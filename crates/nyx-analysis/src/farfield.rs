// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Stratton-Chu near-to-far-field transformation (NFFFT) of 3D frequency-domain electromagnetic results (FR-36, ADR 138).
//!
//! The analysis works only on a [`SolveResult`] and the [`Mesh`] it was solved on: the per-cell complex nodal fields
//! `Ex, Ey, Ez` (V/m) and `Hx, Hy, Hz` (A/m) of one solution (P1 or P2 per-cell layout on `discretise(mesh)`, field names
//! selectable, e.g. the scattered `Esx .. Hsz`) are interpolated onto a closed (or planar) surface, turned into equivalent
//! surface currents `J_s = n x H`, `M_s = -n x E`, and radiated with the homogeneous-medium Green's function.
//!
//! Conventions: time dependence `exp(+j w t)`, outgoing waves `exp(-j k r)`, `H = j/(w mu) curl E`, SI units throughout,
//! background medium of real relative permittivity `eps_r` and permeability `mu_r`. Far fields are reported as `r E`
//! (volts) with the factor `exp(-j k r)` omitted.
//!
//! * surfaces ([`SurfaceSpec`]): labelled boundary facets, the boundary of a set of domains, an axis-aligned box or a
//!   single plane (half-space aperture); mesh faces lying on the surface are integrated with a 7-point degree-5 triangle
//!   rule, otherwise a 2x2 Gauss grid of located points is used.
//! * [`Radiator::far`]: radiation vectors `N`, `L` and `r E = (j k / 4 pi) (r^ x L - eta N_perp)`.
//! * [`Radiator::near`]: the full dyadic (Franz) kernel, valid at any distance from the surface (focal planes, PSF).
//! * totals: radiated power from the far field (Gauss-Legendre x uniform azimuth) and from the surface Poynting flux,
//!   directivity, gain, bistatic RCS, `(kx, ky)` maps, focal-plane intensity with Strehl ratio against the Airy reference,
//!   and a first-version planar-interface (Fresnel) projection.
//!
//! Every parallel loop computes independent outputs with a fixed sequential inner sum, so results are bit-identical for
//! any thread count.

#![allow(clippy::needless_range_loop)]

use nyx_physics::{CellField, Discretisation, Error, FieldFrame, Mesh, SolveResult};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

/// Vacuum permeability (H/m), the value the `ems` kernels use.
pub const MU0: f64 = 4.0 * PI * 1e-7;
/// Speed of light in vacuum (m/s).
pub const C0: f64 = 299_792_458.0;
/// Vacuum permittivity (F/m), `1 / (mu0 c0^2)`.
pub const EPS0: f64 = 1.0 / (MU0 * C0 * C0);

pub(crate) const AXIS_NAMES: [&str; 3] = ["x", "y", "z"];

// ------------------------------------------------------------------------------------------------ complex arithmetic

/// A complex number (the analysis crate has no complex dependency of its own).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }
    /// `exp(j phase)`.
    pub fn cis(phase: f64) -> C64 {
        let (s, c) = phase.sin_cos();
        C64 { re: c, im: s }
    }
    pub fn conj(self) -> C64 {
        C64 { re: self.re, im: -self.im }
    }
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn scale(self, s: f64) -> C64 {
        C64 { re: self.re * s, im: self.im * s }
    }
}

impl Add for C64 {
    type Output = C64;
    fn add(self, o: C64) -> C64 {
        C64 { re: self.re + o.re, im: self.im + o.im }
    }
}
impl Sub for C64 {
    type Output = C64;
    fn sub(self, o: C64) -> C64 {
        C64 { re: self.re - o.re, im: self.im - o.im }
    }
}
impl Mul for C64 {
    type Output = C64;
    fn mul(self, o: C64) -> C64 {
        C64 { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}
impl Mul<f64> for C64 {
    type Output = C64;
    fn mul(self, s: f64) -> C64 {
        self.scale(s)
    }
}
impl Div for C64 {
    type Output = C64;
    fn div(self, o: C64) -> C64 {
        let d = o.norm_sqr();
        C64 { re: (self.re * o.re + self.im * o.im) / d, im: (self.im * o.re - self.re * o.im) / d }
    }
}
impl Neg for C64 {
    type Output = C64;
    fn neg(self) -> C64 {
        C64 { re: -self.re, im: -self.im }
    }
}
impl AddAssign for C64 {
    fn add_assign(&mut self, o: C64) {
        self.re += o.re;
        self.im += o.im;
    }
}

/// A complex Cartesian vector.
pub type CVec3 = [C64; 3];

pub(crate) fn rdot(a: [f64; 3], b: &CVec3) -> C64 {
    b[0] * a[0] + b[1] * a[1] + b[2] * a[2]
}
pub(crate) fn rcross(a: [f64; 3], b: &CVec3) -> CVec3 {
    [b[2] * a[1] - b[1] * a[2], b[0] * a[2] - b[2] * a[0], b[1] * a[0] - b[0] * a[1]]
}
pub(crate) fn ccross(a: &CVec3, b: &CVec3) -> CVec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub(crate) fn cconj(a: &CVec3) -> CVec3 {
    [a[0].conj(), a[1].conj(), a[2].conj()]
}
pub(crate) fn cnorm2(a: &CVec3) -> f64 {
    a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr()
}
pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub(crate) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub(crate) fn unit3(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm3(a);
    (n.is_finite() && n > 0.0).then(|| scale3(a, 1.0 / n))
}
/// Two unit vectors completing `p` to a right-handed orthonormal frame `(e1, e2, p)`.
pub(crate) fn frame_of(p: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let a = if p[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let e1 = unit3(sub3(a, scale3(p, dot3(a, p)))).unwrap_or([1.0, 0.0, 0.0]);
    (e1, cross3(p, e1))
}
/// The in-plane axes `(u, v)` of a plane normal to `axis`, with `e_u x e_v = e_axis`.
pub fn uv_axes(axis: usize) -> (usize, usize) {
    match axis {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    }
}
fn unit_axis(axis: usize) -> [f64; 3] {
    let mut e = [0.0; 3];
    e[axis.min(2)] = 1.0;
    e
}

/// `f(0), .., f(n-1)` in order, on `threads` scoped threads with fixed 64-item chunks dealt round-robin. Every output is
/// computed by the same sequential code, so the result does not depend on `threads`.
pub(crate) fn par_map<R, F>(n: usize, threads: usize, f: F) -> Vec<R>
where
    R: Send + Default + Clone,
    F: Fn(usize) -> R + Sync,
{
    const CHUNK: usize = 64;
    let mut out = vec![R::default(); n];
    let nchunks = n.div_ceil(CHUNK);
    if threads <= 1 || nchunks < 2 {
        for (i, o) in out.iter_mut().enumerate() {
            *o = f(i);
        }
        return out;
    }
    let nt = threads.min(nchunks);
    let mut work: Vec<Vec<(usize, &mut [R])>> = (0..nt).map(|_| Vec::new()).collect();
    for (ci, chunk) in out.chunks_mut(CHUNK).enumerate() {
        work[ci % nt].push((ci, chunk));
    }
    std::thread::scope(|s| {
        let f = &f;
        for w in work {
            s.spawn(move || {
                for (ci, chunk) in w {
                    for (k, o) in chunk.iter_mut().enumerate() {
                        *o = f(ci * CHUNK + k);
                    }
                }
            });
        }
    });
    out
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`, ascending.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n.max(1);
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut pp = 1.0;
        for _ in 0..100 {
            let (mut p1, mut p2) = (1.0, 0.0);
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp;
            if (z - z1).abs() < 1e-15 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        let wi = 2.0 / ((1.0 - z * z) * pp * pp);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}

// ------------------------------------------------------------------------------------------------ medium

/// The homogeneous background medium the surface radiates into, at one frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Medium {
    pub eps_r: f64,
    pub mu_r: f64,
    pub freq_hz: f64,
}

impl Medium {
    pub fn new(eps_r: f64, mu_r: f64, freq_hz: f64) -> Result<Medium, Error> {
        if !(eps_r.is_finite() && eps_r > 0.0 && mu_r.is_finite() && mu_r > 0.0 && freq_hz.is_finite() && freq_hz > 0.0)
        {
            return Err(Error::invalid(format!(
                "the background medium needs eps_r > 0, mu_r > 0 and freq_hz > 0 (got {eps_r}, {mu_r}, {freq_hz})"
            )));
        }
        Ok(Medium { eps_r, mu_r, freq_hz })
    }
    pub fn omega(&self) -> f64 {
        2.0 * PI * self.freq_hz
    }
    /// Refractive index `sqrt(eps_r mu_r)`.
    pub fn index(&self) -> f64 {
        (self.eps_r * self.mu_r).sqrt()
    }
    /// Vacuum wavenumber (1/m).
    pub fn k0(&self) -> f64 {
        self.omega() / C0
    }
    /// Wavenumber in the medium (1/m).
    pub fn k(&self) -> f64 {
        self.k0() * self.index()
    }
    /// Wave impedance of the medium (ohm).
    pub fn eta(&self) -> f64 {
        (MU0 * self.mu_r / (EPS0 * self.eps_r)).sqrt()
    }
    /// Wavelength in the medium (m).
    pub fn wavelength(&self) -> f64 {
        2.0 * PI / self.k()
    }
}

// ------------------------------------------------------------------------------------------------ field sampling

/// Which solution and which per-cell fields carry `E` and `H`.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldSelection {
    pub solution: usize,
    pub e_names: [String; 3],
    pub h_names: [String; 3],
}

impl FieldSelection {
    /// `E{prefix}x .. H{prefix}z` of solution `solution` (`""` gives `Ex .. Hz`, `"s"` the scattered `Esx .. Hsz`).
    pub fn with_prefix(solution: usize, prefix: &str) -> FieldSelection {
        let n = |f: &str, c: &str| format!("{f}{prefix}{c}");
        FieldSelection {
            solution,
            e_names: [n("E", "x"), n("E", "y"), n("E", "z")],
            h_names: [n("H", "x"), n("H", "y"), n("H", "z")],
        }
    }

    /// Reads `solution`, `field_prefix`, `e_fields`, `h_fields` from analysis parameters.
    pub fn from_params(params: &Value) -> Result<FieldSelection, Error> {
        let solution = usize_param(params, "solution", 0)?;
        let prefix = match params.get("field_prefix") {
            None | Some(Value::Null) => "",
            Some(v) => v.as_str().ok_or_else(|| Error::invalid("params.field_prefix must be a string"))?,
        };
        let mut sel = FieldSelection::with_prefix(solution, prefix);
        for (key, slot) in [("e_fields", &mut sel.e_names), ("h_fields", &mut sel.h_names)] {
            if let Some(v) = params.get(key).filter(|v| !v.is_null()) {
                let names: Option<Vec<String>> =
                    v.as_array().map(|a| a.iter().filter_map(|s| s.as_str().map(str::to_string)).collect());
                match names {
                    Some(n) if n.len() == 3 => *slot = [n[0].clone(), n[1].clone(), n[2].clone()],
                    _ => return Err(Error::invalid(format!("params.{key} must be an array of three field names"))),
                }
            }
        }
        Ok(sel)
    }
}

/// Interpolates the selected per-cell complex fields at any point of a cell of `discretise(mesh)`.
pub struct CellSampler<'a> {
    disc: Discretisation,
    nloc: usize,
    e: [&'a CellField; 3],
    h: Option<[&'a CellField; 3]>,
    value: Option<f64>,
}

fn shape_fn(nloc: usize, l: &[f64; 4], out: &mut [f64; 10]) -> usize {
    if nloc == 4 {
        out[..4].copy_from_slice(l);
        return 4;
    }
    for i in 0..4 {
        out[i] = l[i] * (2.0 * l[i] - 1.0);
    }
    let mut k = 4;
    for i in 0..4 {
        for j in i + 1..4 {
            out[k] = 4.0 * l[i] * l[j];
            k += 1;
        }
    }
    10
}

impl<'a> CellSampler<'a> {
    /// The sampler of solution `sel.solution` of `result`; `H` is optional (only the magnetic-current aperture and the
    /// E-only overlap work without it).
    pub fn new(mesh: &Mesh, result: &'a SolveResult, sel: &FieldSelection) -> Result<CellSampler<'a>, Error> {
        if mesh.dim() != 3 {
            return Err(Error::Unsupported("the far-field and beam-overlap analyses need a 3D result".into()));
        }
        let sol = result.solutions.get(sel.solution).ok_or_else(|| {
            Error::invalid(format!(
                "solution {} does not exist (the result has {})",
                sel.solution,
                result.solutions.len()
            ))
        })?;
        let names: Vec<&str> = sol.cell_fields.iter().map(|f| f.name.as_str()).collect();
        let find = |n: &str| sol.cell_fields.iter().find(|f| f.name == n);
        let need = |n: &str| {
            find(n).ok_or_else(|| {
                Error::invalid(format!(
                    "solution {} has no per-cell field '{n}' (it has: {})",
                    sel.solution,
                    names.join(", ")
                ))
            })
        };
        let e = [need(&sel.e_names[0])?, need(&sel.e_names[1])?, need(&sel.e_names[2])?];
        let h = match (find(&sel.h_names[0]), find(&sel.h_names[1]), find(&sel.h_names[2])) {
            (Some(a), Some(b), Some(c)) => Some([a, b, c]),
            _ => None,
        };
        let disc = nyx_physics::discretise(mesh)?;
        let nloc = e[0].nloc;
        if nloc != 4 && nloc != 10 {
            return Err(Error::invalid(format!(
                "per-cell field '{}' has nloc {nloc}; P1 (4) or P2 (10) expected",
                e[0].name
            )));
        }
        let want = disc.n_cells() * nloc;
        for f in e.iter().chain(h.iter().flatten()) {
            if f.nloc != nloc || f.re.len() != want || !(f.im.is_empty() || f.im.len() == want) {
                return Err(Error::invalid(format!(
                    "per-cell field '{}' does not match the mesh ({} cells x nloc {nloc} = {want} values, got {} with nloc {})",
                    f.name,
                    disc.n_cells(),
                    f.re.len(),
                    f.nloc
                )));
            }
        }
        Ok(CellSampler { disc, nloc, e, h, value: sol.value })
    }

    /// The discretised mesh the per-cell fields live on.
    pub fn disc(&self) -> &Discretisation {
        &self.disc
    }
    pub fn has_h(&self) -> bool {
        self.h.is_some()
    }
    /// The real value of the solution (its frequency in Hz for a driven study).
    pub fn solution_value(&self) -> Option<f64> {
        self.value
    }
    /// The names of the E and H fields used.
    pub fn field_names(&self) -> (Vec<String>, Vec<String>) {
        (self.e.iter().map(|f| f.name.clone()).collect(), self.h.iter().flatten().map(|f| f.name.clone()).collect())
    }

    /// `(E, H)` in cell `cell` at barycentric coordinates `lam` (H is zero when the result carries none).
    pub fn eval(&self, cell: usize, lam: &[f64; 4]) -> (CVec3, CVec3) {
        let mut nf = [0.0; 10];
        let n = shape_fn(self.nloc, lam, &mut nf);
        let base = cell * self.nloc;
        let val = |f: &CellField| {
            let mut z = C64::ZERO;
            for a in 0..n {
                let i = base + a;
                z.re += nf[a] * f.re.get(i).copied().unwrap_or(0.0);
                z.im += nf[a] * f.im.get(i).copied().unwrap_or(0.0);
            }
            z
        };
        let e = [val(self.e[0]), val(self.e[1]), val(self.e[2])];
        let h = match &self.h {
            Some(h) => [val(h[0]), val(h[1]), val(h[2])],
            None => [C64::ZERO; 3],
        };
        (e, h)
    }

    fn vtx(&self, i: u32) -> [f64; 3] {
        let p = self.disc.point(i as usize);
        [p[0], p[1], p[2]]
    }
}

/// One quadrature point of the radiating surface: position, weight (m^2), unit normal (pointing to the radiation side)
/// and the fields there.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceSample {
    pub x: [f64; 3],
    pub w: f64,
    pub n: [f64; 3],
    pub e: CVec3,
    pub h: CVec3,
}

/// The surface the far field is computed from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sampling {
    /// mesh faces on the surface when they cover it (within 0.1 %), otherwise the located grid
    Auto,
    /// mesh faces lying on the planes only
    Faces,
    /// a 2x2-Gauss grid of points located in the mesh
    Grid,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceSpec {
    /// boundary facets carrying one of `labels`; the normal points out of the cells of `inside_domains` (or out of the
    /// mesh for exterior facets, or away from the surface centroid)
    Labels { labels: Vec<i32>, inside_domains: Vec<i32> },
    /// the boundary of the union of the cells of the domains `inside` (e.g. every non-absorbing domain)
    Domains { inside: Vec<i32> },
    /// the six faces of the axis-aligned box `[lo, hi]`, normals outward
    Box { lo: [f64; 3], hi: [f64; 3], sampling: Sampling, h: Option<f64> },
    /// one plane `x[axis] = position`, normal `normal_sign * e_axis`, over the rectangle `[lo, hi]` of its `(u, v)` axes
    /// (default: the mesh extent)
    Plane {
        axis: usize,
        position: f64,
        normal_sign: f64,
        lo: Option<[f64; 2]>,
        hi: Option<[f64; 2]>,
        sampling: Sampling,
        h: Option<f64>,
    },
}

/// The sampled surface.
#[derive(Clone, Debug, Default)]
pub struct SampledSurface {
    pub samples: Vec<SurfaceSample>,
    /// integrated area of the samples (m^2)
    pub area: f64,
    /// area the surface should have (box, plane), when known
    pub target_area: Option<f64>,
    /// `faces` or `grid`
    pub sampling: String,
    /// number of triangles or grid cells
    pub n_elements: usize,
    /// grid points that fell outside the mesh (left out)
    pub n_outside: usize,
    /// the plane normal of a planar aperture (the radiating half-space)
    pub half_space_normal: Option<[f64; 3]>,
}

/// Dunavant 7-point degree-5 triangle rule: barycentric point, weight (sum 1).
const TRI_Q: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], 0.225),
    ([0.059715871789770, 0.470142064105115, 0.470142064105115], 0.132394152788506),
    ([0.470142064105115, 0.059715871789770, 0.470142064105115], 0.132394152788506),
    ([0.470142064105115, 0.470142064105115, 0.059715871789770], 0.132394152788506),
    ([0.797426985353087, 0.101286507323456, 0.101286507323456], 0.125939180544827),
    ([0.101286507323456, 0.797426985353087, 0.101286507323456], 0.125939180544827),
    ([0.101286507323456, 0.101286507323456, 0.797426985353087], 0.125939180544827),
];

const NONE: u32 = u32::MAX;

fn face_key(v: [u32; 3]) -> [u32; 3] {
    let mut k = v;
    k.sort_unstable();
    k
}

/// Local vertex triples of the four faces of a tetrahedron (face `k` is opposite vertex `k`).
const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Cells sharing every face (`NONE` for a missing second cell).
fn face_map(m: &Discretisation) -> HashMap<[u32; 3], [u32; 2]> {
    let mut map: HashMap<[u32; 3], [u32; 2]> = HashMap::with_capacity(m.n_cells() * 2);
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for f in TET_FACES {
            let e = map.entry(face_key([v[f[0]], v[f[1]], v[f[2]]])).or_insert([NONE, NONE]);
            if e[0] == NONE {
                e[0] = c as u32;
            } else if e[1] == NONE {
                e[1] = c as u32;
            }
        }
    }
    map
}

fn local_positions(cell: &[u32], verts: &[u32; 3]) -> Option<[usize; 3]> {
    let mut out = [0; 3];
    for (k, v) in verts.iter().enumerate() {
        out[k] = cell.iter().position(|x| x == v)?;
    }
    Some(out)
}

impl CellSampler<'_> {
    /// Quadrature samples of triangle `verts` with normal `n`, the fields averaged over `cells`; returns the area.
    fn push_face(&self, out: &mut Vec<SurfaceSample>, verts: [u32; 3], cells: &[usize], n: [f64; 3]) -> f64 {
        let p = [self.vtx(verts[0]), self.vtx(verts[1]), self.vtx(verts[2])];
        let flat_area = 0.5 * norm3(cross3(sub3(p[1], p[0]), sub3(p[2], p[0])));
        let locs: Vec<(usize, [usize; 3])> =
            cells.iter().filter_map(|&c| local_positions(self.disc.cell(c), &verts).map(|l| (c, l))).collect();
        if locs.is_empty() || flat_area <= 0.0 {
            return 0.0;
        }
        let curved = self
            .disc
            .cell_is_curved(locs[0].0)
            .then(|| (femcore::fem::CurvedTet::new(self.disc.cell_geometry_nodes(locs[0].0)), locs[0].1));
        let inv = 1.0 / locs.len() as f64;
        let mut area_acc = 0.0;
        for (xi, wq) in TRI_Q {
            let (x, area_q, n_q) = match &curved {
                Some((ct, lf)) => {
                    let mut lam4 = [0.0; 4];
                    for i in 0..3 {
                        lam4[lf[i]] = xi[i];
                    }
                    let d = ct.dx_dlam(&lam4);
                    let nr = cross3(sub3(d[lf[1]], d[lf[0]]), sub3(d[lf[2]], d[lf[0]]));
                    let nl = norm3(nr);
                    if nl > 1e-30 {
                        let mut nq = [nr[0] / nl, nr[1] / nl, nr[2] / nl];
                        if dot3(nq, n) < 0.0 {
                            nq = scale3(nq, -1.0);
                        }
                        (ct.position(&lam4), 0.5 * nl, nq)
                    } else {
                        (ct.position(&lam4), flat_area, n)
                    }
                }
                None => {
                    let mut x = [0.0; 3];
                    for i in 0..3 {
                        for k in 0..3 {
                            x[k] += xi[i] * p[i][k];
                        }
                    }
                    (x, flat_area, n)
                }
            };
            let mut e = [C64::ZERO; 3];
            let mut h = [C64::ZERO; 3];
            for &(c, loc) in &locs {
                let mut lam = [0.0; 4];
                for i in 0..3 {
                    lam[loc[i]] = xi[i];
                }
                let (ec, hc) = self.eval(c, &lam);
                for k in 0..3 {
                    e[k] += ec[k] * inv;
                    h[k] += hc[k] * inv;
                }
            }
            let w = wq * area_q;
            area_acc += w;
            out.push(SurfaceSample { x, w, n: n_q, e, h });
        }
        area_acc
    }

    /// The unit normal of triangle `verts` oriented so that it points away from `from`.
    fn oriented_normal(&self, verts: [u32; 3], from: [f64; 3]) -> [f64; 3] {
        let p = [self.vtx(verts[0]), self.vtx(verts[1]), self.vtx(verts[2])];
        let n = unit3(cross3(sub3(p[1], p[0]), sub3(p[2], p[0]))).unwrap_or([0.0, 0.0, 1.0]);
        let c =
            scale3([p[0][0] + p[1][0] + p[2][0], p[0][1] + p[1][1] + p[2][1], p[0][2] + p[1][2] + p[2][2]], 1.0 / 3.0);
        if dot3(n, sub3(c, from)) < 0.0 {
            scale3(n, -1.0)
        } else {
            n
        }
    }

    fn opposite_vertex(&self, cell: usize, verts: &[u32; 3]) -> [f64; 3] {
        let v = self.disc.cell(cell);
        let o = v.iter().copied().find(|x| !verts.contains(x)).unwrap_or(v[0]);
        self.vtx(o)
    }

    /// Mesh faces on the plane `x[axis] = pos` inside the rectangle `[lo, hi]` of `(u, v)`; returns (area, faces).
    #[allow(clippy::too_many_arguments)]
    fn plane_faces(
        &self,
        fm: &HashMap<[u32; 3], [u32; 2]>,
        axis: usize,
        pos: f64,
        lo: [f64; 2],
        hi: [f64; 2],
        n_sign: f64,
        tol: f64,
        out: &mut Vec<SurfaceSample>,
    ) -> (f64, usize) {
        let (iu, iv) = uv_axes(axis);
        let n = scale3(unit_axis(axis), n_sign);
        let mut seen: HashSet<[u32; 3]> = HashSet::new();
        let (mut area, mut count) = (0.0, 0);
        for c in 0..self.disc.n_cells() {
            let v = self.disc.cell(c);
            for f in TET_FACES {
                let verts = [v[f[0]], v[f[1]], v[f[2]]];
                let p = [self.vtx(verts[0]), self.vtx(verts[1]), self.vtx(verts[2])];
                if p.iter().any(|q| (q[axis] - pos).abs() > tol) {
                    continue;
                }
                let cu = (p[0][iu] + p[1][iu] + p[2][iu]) / 3.0;
                let cv = (p[0][iv] + p[1][iv] + p[2][iv]) / 3.0;
                if cu < lo[0] - tol || cu > hi[0] + tol || cv < lo[1] - tol || cv > hi[1] + tol {
                    continue;
                }
                let key = face_key(verts);
                if !seen.insert(key) {
                    continue;
                }
                let cells: Vec<usize> = fm
                    .get(&key)
                    .map(|e| e.iter().filter(|&&x| x != NONE).map(|&x| x as usize).collect())
                    .unwrap_or_default();
                area += self.push_face(out, verts, &cells, n);
                count += 1;
            }
        }
        (area, count)
    }

    /// A 2x2-Gauss grid over the rectangle of a plane, points located in the mesh; returns (area, cells, outside).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn plane_grid(
        &self,
        loc: &Locator,
        axis: usize,
        pos: f64,
        lo: [f64; 2],
        hi: [f64; 2],
        nuv: [usize; 2],
        n_sign: f64,
        out: &mut Vec<SurfaceSample>,
    ) -> (f64, usize, usize) {
        let n = scale3(unit_axis(axis), n_sign);
        let pts = plane_grid_points(axis, pos, lo, hi, nuv);
        let mut outside = 0;
        let mut area = 0.0;
        for (x, w) in pts {
            match loc.locate(x) {
                Some((c, lam)) => {
                    let (e, h) = self.eval(c, &lam);
                    out.push(SurfaceSample { x, w, n, e, h });
                    area += w;
                }
                None => outside += 1,
            }
        }
        (area, nuv[0] * nuv[1], outside)
    }

    /// A default grid spacing: half the mean edge length of the mesh (the P2 node spacing).
    pub(crate) fn default_h(&self) -> f64 {
        let edges = self.disc.edges();
        if edges.is_empty() {
            return 1.0;
        }
        let mut s = 0.0;
        for &(a, b) in edges {
            s += norm3(sub3(self.vtx(a), self.vtx(b)));
        }
        0.5 * s / edges.len() as f64
    }
}

/// The 2x2-Gauss points `(x, weight)` of an `nu x nv` grid over the rectangle `[lo, hi]` of the plane `x[axis] = pos`.
pub fn plane_grid_points(axis: usize, pos: f64, lo: [f64; 2], hi: [f64; 2], nuv: [usize; 2]) -> Vec<([f64; 3], f64)> {
    let (iu, iv) = uv_axes(axis);
    let (nu, nv) = (nuv[0].max(1), nuv[1].max(1));
    let du = (hi[0] - lo[0]) / nu as f64;
    let dv = (hi[1] - lo[1]) / nv as f64;
    let g = 0.5 / 3f64.sqrt();
    let w = (du * dv).abs() / 4.0;
    let mut out = Vec::with_capacity(nu * nv * 4);
    for j in 0..nv {
        for gj in [-g, g] {
            let v = lo[1] + (j as f64 + 0.5 + gj) * dv;
            for i in 0..nu {
                for gi in [-g, g] {
                    let u = lo[0] + (i as f64 + 0.5 + gi) * du;
                    let mut x = [0.0; 3];
                    x[axis] = pos;
                    x[iu] = u;
                    x[iv] = v;
                    out.push((x, w));
                }
            }
        }
    }
    out
}

/// Point location in a tetrahedral mesh through a uniform grid of bins.
pub(crate) struct Locator<'m> {
    m: &'m Discretisation,
    lo: [f64; 3],
    size: [f64; 3],
    dims: [usize; 3],
    start: Vec<usize>,
    items: Vec<u32>,
}

impl<'m> Locator<'m> {
    pub(crate) fn new(m: &'m Discretisation) -> Locator<'m> {
        let (lo, hi) = m.bbox();
        let ext = [(hi[0] - lo[0]).max(1e-300), (hi[1] - lo[1]).max(1e-300), (hi[2] - lo[2]).max(1e-300)];
        let nb = m.n_cells().max(1) as f64;
        let s = (ext[0] * ext[1] * ext[2] / nb).cbrt().max(1e-300);
        let dims = [0, 1, 2].map(|k| ((ext[k] / s).ceil() as usize).clamp(1, 512));
        let size = [0, 1, 2].map(|k| ext[k] / dims[k] as f64);
        let mut loc = Locator { m, lo, size, dims, start: Vec::new(), items: Vec::new() };
        let nbins = dims[0] * dims[1] * dims[2];
        let ranges: Vec<([usize; 3], [usize; 3])> = (0..m.n_cells()).map(|c| loc.cell_range(c)).collect();
        let mut count = vec![0usize; nbins + 1];
        for (a, b) in &ranges {
            for k in a[2]..=b[2] {
                for j in a[1]..=b[1] {
                    for i in a[0]..=b[0] {
                        count[loc.bin(i, j, k) + 1] += 1;
                    }
                }
            }
        }
        for i in 0..nbins {
            count[i + 1] += count[i];
        }
        let mut cursor = count.clone();
        let mut items = vec![0u32; count[nbins]];
        for (c, (a, b)) in ranges.iter().enumerate() {
            for k in a[2]..=b[2] {
                for j in a[1]..=b[1] {
                    for i in a[0]..=b[0] {
                        let bi = loc.bin(i, j, k);
                        items[cursor[bi]] = c as u32;
                        cursor[bi] += 1;
                    }
                }
            }
        }
        loc.start = count;
        loc.items = items;
        loc
    }

    fn bin(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.dims[1] + j) * self.dims[0] + i
    }

    fn index(&self, x: f64, k: usize) -> usize {
        let t = ((x - self.lo[k]) / self.size[k]).floor();
        if t <= 0.0 {
            0
        } else {
            (t as usize).min(self.dims[k] - 1)
        }
    }

    fn cell_range(&self, c: usize) -> ([usize; 3], [usize; 3]) {
        let mut a = [usize::MAX; 3];
        let mut b = [0; 3];
        for &v in self.m.cell(c) {
            let p = self.m.point(v as usize);
            for k in 0..3 {
                let i = self.index(p[k], k);
                a[k] = a[k].min(i);
                b[k] = b[k].max(i);
            }
        }
        (a, b)
    }

    fn bary(&self, c: usize, x: [f64; 3]) -> [f64; 4] {
        let v = self.m.cell(c);
        let p = |i: usize| {
            let q = self.m.point(v[i] as usize);
            [q[0], q[1], q[2]]
        };
        let p0 = p(0);
        let (a, b, cc) = (sub3(p(1), p0), sub3(p(2), p0), sub3(p(3), p0));
        let det = dot3(a, cross3(b, cc));
        if det.abs() < 1e-300 {
            return [f64::NEG_INFINITY; 4];
        }
        let r = sub3(x, p0);
        let l1 = dot3(r, cross3(b, cc)) / det;
        let l2 = dot3(a, cross3(r, cc)) / det;
        let l3 = dot3(a, cross3(b, r)) / det;
        [1.0 - l1 - l2 - l3, l1, l2, l3]
    }

    /// The cell containing `x` (the one with the largest smallest barycentric coordinate) and its barycentric
    /// coordinates, clamped to the cell; `None` outside the mesh.
    pub(crate) fn locate(&self, x: [f64; 3]) -> Option<(usize, [f64; 4])> {
        let b = self.bin(self.index(x[0], 0), self.index(x[1], 1), self.index(x[2], 2));
        let mut best: Option<(usize, [f64; 4], f64)> = None;
        for &c in &self.items[self.start[b]..self.start[b + 1]] {
            let lam = self.bary(c as usize, x);
            let m = lam.iter().copied().fold(f64::INFINITY, f64::min);
            if best.as_ref().is_none_or(|bb| m > bb.2) {
                best = Some((c as usize, lam, m));
            }
        }
        let (c, lam, m) = best?;
        if m < -1e-8 {
            return None;
        }
        let mut l = lam.map(|x| x.max(0.0));
        let s: f64 = l.iter().sum();
        if s > 0.0 {
            l = l.map(|x| x / s);
        }
        Some((c, l))
    }
}

/// Samples the surface `spec` with the fields of `sampler`.
pub fn sample_surface(sampler: &CellSampler, spec: &SurfaceSpec) -> Result<SampledSurface, Error> {
    let disc = sampler.disc();
    let mut out = SampledSurface { sampling: "faces".into(), ..Default::default() };
    match spec {
        SurfaceSpec::Labels { labels, inside_domains } => {
            let fm = face_map(disc);
            let sel: HashSet<i32> = labels.iter().copied().collect();
            let inside: HashSet<i32> = inside_domains.iter().copied().collect();
            let facets: Vec<usize> = (0..disc.n_facets()).filter(|&f| sel.contains(&disc.facet_labels[f])).collect();
            if facets.is_empty() {
                return Err(Error::invalid(format!(
                    "no boundary facet carries the labels {labels:?} (the mesh has {:?})",
                    disc.boundary_labels()
                )));
            }
            let mut cen = [0.0; 3];
            let mut atot = 0.0;
            for &f in &facets {
                let v = disc.facet(f);
                let p = [sampler.vtx(v[0]), sampler.vtx(v[1]), sampler.vtx(v[2])];
                let a = 0.5 * norm3(cross3(sub3(p[1], p[0]), sub3(p[2], p[0])));
                for k in 0..3 {
                    cen[k] += a * (p[0][k] + p[1][k] + p[2][k]) / 3.0;
                }
                atot += a;
            }
            let cen = scale3(cen, 1.0 / atot.max(1e-300));
            for &f in &facets {
                let v = disc.facet(f);
                let verts = [v[0], v[1], v[2]];
                let cells: Vec<usize> = fm
                    .get(&face_key(verts))
                    .map(|e| e.iter().filter(|&&x| x != NONE).map(|&x| x as usize).collect())
                    .unwrap_or_default();
                if cells.is_empty() {
                    return Err(Error::invalid(format!("boundary facet {f} is not a face of any cell")));
                }
                let in_cell = cells.iter().copied().find(|&c| inside.contains(&disc.cell_labels[c]));
                let reference = in_cell.or((cells.len() == 1).then_some(cells[0]));
                let n = match reference {
                    Some(c) => sampler.oriented_normal(verts, sampler.opposite_vertex(c, &verts)),
                    None => sampler.oriented_normal(verts, cen),
                };
                out.area += sampler.push_face(&mut out.samples, verts, &cells, n);
                out.n_elements += 1;
            }
        }
        SurfaceSpec::Domains { inside } => {
            let fm = face_map(disc);
            let set: HashSet<i32> = inside.iter().copied().collect();
            for c in 0..disc.n_cells() {
                if !set.contains(&disc.cell_labels[c]) {
                    continue;
                }
                let v = disc.cell(c);
                for (k, f) in TET_FACES.iter().enumerate() {
                    let verts = [v[f[0]], v[f[1]], v[f[2]]];
                    let on_surface = fm.get(&face_key(verts)).is_none_or(|e| {
                        e.iter().all(|&o| o == NONE || o as usize == c || !set.contains(&disc.cell_labels[o as usize]))
                    });
                    if on_surface {
                        let n = sampler.oriented_normal(verts, sampler.vtx(v[k]));
                        out.area += sampler.push_face(&mut out.samples, verts, &[c], n);
                        out.n_elements += 1;
                    }
                }
            }
            if out.samples.is_empty() {
                return Err(Error::invalid(format!("the domains {inside:?} have no boundary in the mesh")));
            }
        }
        SurfaceSpec::Box { lo, hi, sampling, h } => {
            if (0..3).any(|k| !lo[k].is_finite() || !hi[k].is_finite() || hi[k] <= lo[k]) {
                return Err(Error::invalid("surface box needs hi > lo on every axis"));
            }
            let ext = sub3(*hi, *lo);
            let target = 2.0 * (ext[0] * ext[1] + ext[1] * ext[2] + ext[2] * ext[0]);
            out.target_area = Some(target);
            let tol = 1e-8 * norm3(ext);
            let fm = if *sampling != Sampling::Grid { face_map(disc) } else { HashMap::new() };
            let planes = |out: &mut SampledSurface, grid: Option<(&Locator, f64)>| {
                for axis in 0..3 {
                    let (iu, iv) = uv_axes(axis);
                    let (plo, phi) = ([lo[iu], lo[iv]], [hi[iu], hi[iv]]);
                    for (pos, sign) in [(lo[axis], -1.0), (hi[axis], 1.0)] {
                        match grid {
                            None => {
                                let (a, n) = sampler.plane_faces(&fm, axis, pos, plo, phi, sign, tol, &mut out.samples);
                                out.area += a;
                                out.n_elements += n;
                            }
                            Some((loc, hh)) => {
                                let nuv = [
                                    ((phi[0] - plo[0]) / hh).ceil() as usize,
                                    ((phi[1] - plo[1]) / hh).ceil() as usize,
                                ];
                                let (a, n, o) =
                                    sampler.plane_grid(loc, axis, pos, plo, phi, nuv, sign, &mut out.samples);
                                out.area += a;
                                out.n_elements += n;
                                out.n_outside += o;
                            }
                        }
                    }
                }
            };
            if *sampling != Sampling::Grid {
                planes(&mut out, None);
            }
            let covered = (out.area - target).abs() <= 1e-3 * target;
            if *sampling == Sampling::Grid || (*sampling == Sampling::Auto && !covered) {
                out = SampledSurface { target_area: Some(target), sampling: "grid".into(), ..Default::default() };
                let loc = Locator::new(disc);
                let hh = h.unwrap_or_else(|| sampler.default_h()).max(1e-300);
                planes(&mut out, Some((&loc, hh)));
            }
        }
        SurfaceSpec::Plane { axis, position, normal_sign, lo, hi, sampling, h } => {
            let axis = (*axis).min(2);
            let (iu, iv) = uv_axes(axis);
            let (blo, bhi) = disc.bbox();
            let plo = lo.unwrap_or([blo[iu], blo[iv]]);
            let phi = hi.unwrap_or([bhi[iu], bhi[iv]]);
            if !(phi[0] > plo[0] && phi[1] > plo[1]) {
                return Err(Error::invalid("surface plane needs hi > lo on both in-plane axes"));
            }
            let target = (phi[0] - plo[0]) * (phi[1] - plo[1]);
            let sign = if *normal_sign < 0.0 { -1.0 } else { 1.0 };
            out.target_area = Some(target);
            out.half_space_normal = Some(scale3(unit_axis(axis), sign));
            let tol = 1e-8 * norm3(sub3(bhi, blo));
            if *sampling != Sampling::Grid {
                let (a, n) =
                    sampler.plane_faces(&face_map(disc), axis, *position, plo, phi, sign, tol, &mut out.samples);
                out.area = a;
                out.n_elements = n;
            }
            let covered = (out.area - target).abs() <= 1e-3 * target;
            if *sampling == Sampling::Grid || (*sampling == Sampling::Auto && !covered) {
                out.samples.clear();
                out.sampling = "grid".into();
                let loc = Locator::new(disc);
                let hh = h.unwrap_or_else(|| sampler.default_h()).max(1e-300);
                let nuv = [((phi[0] - plo[0]) / hh).ceil() as usize, ((phi[1] - plo[1]) / hh).ceil() as usize];
                let (a, n, o) = sampler.plane_grid(&loc, axis, *position, plo, phi, nuv, sign, &mut out.samples);
                out.area = a;
                out.n_elements = n;
                out.n_outside = o;
            }
        }
    }
    if out.samples.is_empty() {
        return Err(Error::invalid("the far-field surface has no quadrature point inside the mesh"));
    }
    Ok(out)
}

/// Time-averaged power through the samples, `1/2 Re sum w (E x H*) . n` (W).
pub fn surface_flux(samples: &[SurfaceSample]) -> f64 {
    let mut p = 0.0;
    for s in samples {
        let sv = ccross(&s.e, &cconj(&s.h));
        p += s.w * (sv[0].re * s.n[0] + sv[1].re * s.n[1] + sv[2].re * s.n[2]);
    }
    0.5 * p
}

// ------------------------------------------------------------------------------------------------ radiation

/// Which equivalent currents radiate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Equivalence {
    /// `J = n x H` and `M = -n x E` (Love; exact for a closed surface)
    Love,
    /// `2 M` only (PEC image of an infinite plane, radiating into the half-space of the normal)
    Magnetic,
    /// `2 J` only (PMC image of an infinite plane)
    Electric,
}

impl Equivalence {
    pub fn parse(s: &str) -> Result<Equivalence, Error> {
        match s {
            "love" => Ok(Equivalence::Love),
            "magnetic" => Ok(Equivalence::Magnetic),
            "electric" => Ok(Equivalence::Electric),
            _ => Err(Error::invalid(format!("equivalence '{s}' is not one of love, magnetic, electric"))),
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Equivalence::Love => "love",
            Equivalence::Magnetic => "magnetic",
            Equivalence::Electric => "electric",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Src {
    x: [f64; 3],
    j: CVec3,
    m: CVec3,
}

/// The equivalent surface currents of a sampled surface in a homogeneous medium, ready to radiate.
pub struct Radiator {
    k: f64,
    eta: f64,
    wmu: f64,
    weps: f64,
    srcs: Vec<Src>,
    half_space: Option<[f64; 3]>,
    threads: usize,
    min_r: f64,
}

impl Radiator {
    /// The radiator of `samples`; `half_space` limits the far field to the directions with `dir . n >= 0`.
    pub fn new(
        samples: &[SurfaceSample],
        medium: &Medium,
        equivalence: Equivalence,
        half_space: Option<[f64; 3]>,
        threads: usize,
    ) -> Radiator {
        let (fj, fm) = match equivalence {
            Equivalence::Love => (1.0, 1.0),
            Equivalence::Magnetic => (0.0, 2.0),
            Equivalence::Electric => (2.0, 0.0),
        };
        let srcs = samples
            .iter()
            .map(|s| {
                let j = rcross(s.n, &s.h);
                let m = rcross(s.n, &s.e);
                Src { x: s.x, j: j.map(|z| z * (fj * s.w)), m: m.map(|z| z * (-fm * s.w)) }
            })
            .collect();
        let k = medium.k();
        Radiator {
            k,
            eta: medium.eta(),
            wmu: medium.omega() * MU0 * medium.mu_r,
            weps: medium.omega() * EPS0 * medium.eps_r,
            srcs,
            half_space,
            threads: threads.max(1),
            min_r: 1e-9 * 2.0 * PI / k,
        }
    }

    /// Wavenumber of the medium.
    pub fn k(&self) -> f64 {
        self.k
    }
    /// Wave impedance of the medium.
    pub fn eta(&self) -> f64 {
        self.eta
    }

    /// `r E(r^)` (V, Cartesian, the factor `exp(-j k r)` left out) in the unit direction `dir`.
    pub fn far(&self, dir: [f64; 3]) -> CVec3 {
        if let Some(n0) = self.half_space {
            if dot3(dir, n0) < 0.0 {
                return [C64::ZERO; 3];
            }
        }
        let mut n = [C64::ZERO; 3];
        let mut l = [C64::ZERO; 3];
        for s in &self.srcs {
            let ph = C64::cis(self.k * dot3(dir, s.x));
            for i in 0..3 {
                n[i] += s.j[i] * ph;
                l[i] += s.m[i] * ph;
            }
        }
        let rn = rdot(dir, &n);
        let rxl = rcross(dir, &l);
        let pref = C64::new(0.0, self.k / (4.0 * PI));
        let mut out = [C64::ZERO; 3];
        for i in 0..3 {
            let nperp = n[i] - rn * dir[i];
            out[i] = pref * (rxl[i] - nperp * self.eta);
        }
        out
    }

    /// `far` over many directions (parallel, order kept).
    pub fn far_many(&self, dirs: &[[f64; 3]]) -> Vec<CVec3> {
        par_map(dirs.len(), self.threads, |i| self.far(dirs[i]))
    }

    /// Radiation intensity `U = |r E|^2 / (2 eta)` (W/sr).
    pub fn intensity(&self, r_e: &CVec3) -> f64 {
        cnorm2(r_e) / (2.0 * self.eta)
    }

    /// `(E, H)` at the point `x` with the full dyadic kernel (any distance; `x` must not lie on the surface).
    pub fn near(&self, x: [f64; 3]) -> (CVec3, CVec3) {
        let k = self.k;
        let mut e = [C64::ZERO; 3];
        let mut h = [C64::ZERO; 3];
        let mjwmu = C64::new(0.0, -self.wmu);
        let mjweps = C64::new(0.0, -self.weps);
        for s in &self.srcs {
            let rv = sub3(x, s.x);
            let r = norm3(rv);
            if r < self.min_r {
                continue;
            }
            let rh = scale3(rv, 1.0 / r);
            let kr = k * r;
            let g = C64::cis(-kr) * (1.0 / (4.0 * PI * r));
            let a = C64::new(1.0 - 1.0 / (kr * kr), -1.0 / kr) * g;
            let b = C64::new(-1.0 + 3.0 / (kr * kr), 3.0 / kr) * g;
            let c = g * C64::new(1.0 / r, k);
            let rj = rdot(rh, &s.j);
            let rm = rdot(rh, &s.m);
            let rxm = rcross(rh, &s.m);
            let rxj = rcross(rh, &s.j);
            for i in 0..3 {
                let gj = a * s.j[i] + b * rj * rh[i];
                let gm = a * s.m[i] + b * rm * rh[i];
                e[i] += mjwmu * gj + c * rxm[i];
                h[i] += mjweps * gm - c * rxj[i];
            }
        }
        (e, h)
    }

    /// `near` at many points (parallel, order kept).
    pub fn near_many(&self, xs: &[[f64; 3]]) -> Vec<(CVec3, CVec3)> {
        par_map(xs.len(), self.threads, |i| self.near(xs[i]))
    }

    /// Radiated power (W): `U` integrated over the sphere (or the hemisphere around `pole`) with `n_theta` Gauss-Legendre
    /// nodes in `cos theta` and `n_phi` uniform azimuths.
    pub fn radiated_power(&self, pole: [f64; 3], hemisphere: bool, n_theta: usize, n_phi: usize) -> f64 {
        let pole = unit3(pole).unwrap_or([0.0, 0.0, 1.0]);
        let (e1, e2) = frame_of(pole);
        let (x, w) = gauss_legendre(n_theta);
        let n_phi = n_phi.max(1);
        let (mu, wmu): (Vec<f64>, Vec<f64>) = if hemisphere {
            (x.iter().map(|t| 0.5 * (t + 1.0)).collect(), w.iter().map(|t| 0.5 * t).collect())
        } else {
            (x, w)
        };
        let mut dirs = Vec::with_capacity(mu.len() * n_phi);
        for &m in &mu {
            let st = (1.0 - m * m).max(0.0).sqrt();
            for j in 0..n_phi {
                let (sp, cp) = (2.0 * PI * j as f64 / n_phi as f64).sin_cos();
                dirs.push([0, 1, 2].map(|k| st * (cp * e1[k] + sp * e2[k]) + m * pole[k]));
            }
        }
        let fields = self.far_many(&dirs);
        let dphi = 2.0 * PI / n_phi as f64;
        let mut p = 0.0;
        for (i, wi) in wmu.iter().enumerate() {
            for j in 0..n_phi {
                p += wi * dphi * self.intensity(&fields[i * n_phi + j]);
            }
        }
        p
    }
}

/// Spherical unit vectors `(r^, theta^, phi^)` of the global frame (pole `+z`).
pub fn spherical_frame(theta: f64, phi: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    ([st * cp, st * sp, ct], [ct * cp, ct * sp, -st], [-sp, cp, 0.0])
}

/// First zero of `J1` (the Airy first dark ring is at `AIRY_ZERO / (k0 NA)`).
pub const AIRY_ZERO: f64 = 3.8317059702075;

/// The Fresnel transmission of a far field into a second half-space (`stratified`, first version): the direction
/// `dir2` in medium 2 (index `n2`) is fed by the medium-1 direction with the same transverse wavevector; the s and p
/// components are scaled by `t_s`, `t_p` and the stationary-phase factor `k2 cos2 / (k1 cos1)`, with the phase of an
/// interface at `x[axis] = z0`. Directions fed by evanescent medium-1 waves give zero.
#[allow(clippy::too_many_arguments)]
pub fn transmit_far_field(
    rad: &Radiator,
    medium: &Medium,
    axis: usize,
    z0: f64,
    side: f64,
    eps2: f64,
    dir2: [f64; 3],
) -> CVec3 {
    let ea = scale3(unit_axis(axis), if side < 0.0 { -1.0 } else { 1.0 });
    let n1 = medium.index();
    let n2 = (eps2 * medium.mu_r).sqrt();
    let k0 = medium.k0();
    let c2 = dot3(dir2, ea);
    if c2 <= 0.0 {
        return [C64::ZERO; 3];
    }
    let t = sub3(dir2, scale3(ea, c2));
    let s2 = norm3(t);
    let s1 = n2 * s2 / n1;
    if s1 >= 1.0 {
        return [C64::ZERO; 3];
    }
    let c1 = (1.0 - s1 * s1).sqrt();
    let th = if s2 > 0.0 { scale3(t, 1.0 / s2) } else { frame_of(ea).0 };
    let dir1 = [0, 1, 2].map(|k| s1 * th[k] + c1 * ea[k]);
    let e1 = rad.far(dir1);
    let sh = unit3(cross3(dir1, ea)).unwrap_or(frame_of(ea).1);
    let p1 = cross3(dir1, sh);
    let p2 = cross3(dir2, sh);
    let es = rdot(sh, &e1);
    let ep = rdot(p1, &e1);
    let ts = 2.0 * n1 * c1 / (n1 * c1 + n2 * c2);
    let tp = 2.0 * n1 * c1 / (n2 * c1 + n1 * c2);
    let factor = (n2 * c2) / (n1 * c1);
    let zs = z0 * if side < 0.0 { -1.0 } else { 1.0 };
    let phase = C64::cis(-(k0 * n1 * c1 - k0 * n2 * c2) * zs) * factor;
    let mut out = [C64::ZERO; 3];
    for k in 0..3 {
        out[k] = (es * (ts * sh[k]) + ep * (tp * p2[k])) * phase;
    }
    out
}

// ------------------------------------------------------------------------------------------------ parameters

pub(crate) fn f64_param(p: &Value, key: &str, default: f64) -> Result<f64, Error> {
    Ok(opt_f64(p, key)?.unwrap_or(default))
}
pub(crate) fn opt_f64(p: &Value, key: &str) -> Result<Option<f64>, Error> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => v
            .as_f64()
            .filter(|x| x.is_finite())
            .map(Some)
            .ok_or_else(|| Error::invalid(format!("params.{key} must be a finite number"))),
    }
}
pub(crate) fn usize_param(p: &Value, key: &str, default: usize) -> Result<usize, Error> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(v) => v
            .as_u64()
            .map(|x| x as usize)
            .ok_or_else(|| Error::invalid(format!("params.{key} must be a non-negative integer"))),
    }
}
fn f64_array(p: &Value, key: &str) -> Result<Option<Vec<f64>>, Error> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => {
            let a = v.as_array().ok_or_else(|| Error::invalid(format!("params.{key} must be an array of numbers")))?;
            let out: Option<Vec<f64>> = a.iter().map(|x| x.as_f64().filter(|y| y.is_finite())).collect();
            out.map(Some).ok_or_else(|| Error::invalid(format!("params.{key} must be an array of finite numbers")))
        }
    }
}
pub(crate) fn vec3_param(p: &Value, key: &str) -> Result<Option<[f64; 3]>, Error> {
    match f64_array(p, key)? {
        None => Ok(None),
        Some(a) if a.len() == 3 => Ok(Some([a[0], a[1], a[2]])),
        Some(_) => Err(Error::invalid(format!("params.{key} must have three numbers"))),
    }
}
pub(crate) fn vec2_param(p: &Value, key: &str) -> Result<Option<[f64; 2]>, Error> {
    match f64_array(p, key)? {
        None => Ok(None),
        Some(a) if a.len() == 2 => Ok(Some([a[0], a[1]])),
        Some(_) => Err(Error::invalid(format!("params.{key} must have two numbers"))),
    }
}
pub(crate) fn axis_param(p: &Value, key: &str, default: usize) -> Result<usize, Error> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::String(s)) => match s.as_str() {
            "x" => Ok(0),
            "y" => Ok(1),
            "z" => Ok(2),
            _ => Err(Error::invalid(format!("params.{key} must be x, y or z"))),
        },
        Some(v) => v
            .as_u64()
            .filter(|&a| a < 3)
            .map(|a| a as usize)
            .ok_or_else(|| Error::invalid(format!("params.{key} must be x, y or z"))),
    }
}
fn i32_list(p: &Value, key: &str) -> Result<Vec<i32>, Error> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(vec![]),
        Some(v) => {
            let a = v.as_array().ok_or_else(|| Error::invalid(format!("params.{key} must be an array of integers")))?;
            a.iter()
                .map(|x| {
                    x.as_i64()
                        .map(|y| y as i32)
                        .ok_or_else(|| Error::invalid(format!("params.{key} must hold integers")))
                })
                .collect()
        }
    }
}
fn sampling_param(p: &Value) -> Result<Sampling, Error> {
    match p.get("sampling").and_then(Value::as_str) {
        None | Some("auto") => Ok(Sampling::Auto),
        Some("faces") => Ok(Sampling::Faces),
        Some("grid") => Ok(Sampling::Grid),
        Some(s) => Err(Error::invalid(format!("sampling '{s}' is not one of auto, faces, grid"))),
    }
}
/// `[min, max, n]` -> `n` points from `min` to `max` (one point `min` when `n == 1`).
fn range_param(p: &Value, key: &str) -> Result<Option<Vec<f64>>, Error> {
    match f64_array(p, key)? {
        None => Ok(None),
        Some(a) if a.len() == 3 && a[2] >= 1.0 => {
            let n = a[2] as usize;
            Ok(Some(linspace(a[0], a[1], n)))
        }
        Some(_) => Err(Error::invalid(format!("params.{key} must be [min, max, n] with n >= 1"))),
    }
}
fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![a];
    }
    (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect()
}

/// The background medium of analysis parameters (`eps_r`, `mu_r`, `freq_hz`; the frequency defaults to the solution
/// value).
pub fn medium_from_params(params: &Value, solution_value: Option<f64>) -> Result<Medium, Error> {
    let f = match opt_f64(params, "freq_hz")? {
        Some(f) => f,
        None => solution_value
            .ok_or_else(|| Error::invalid("params.freq_hz is needed: the solution carries no frequency"))?,
    };
    Medium::new(f64_param(params, "eps_r", 1.0)?, f64_param(params, "mu_r", 1.0)?, f)
}

/// The surface of analysis parameters; without one, the boundary of the non-absorbing domains the result's frame
/// states.
pub fn surface_from_params(
    surface: Option<&Value>,
    result: &SolveResult,
    disc: &Discretisation,
) -> Result<SurfaceSpec, Error> {
    let s = match surface {
        None | Some(Value::Null) => {
            let outside = FieldFrame::of(result).map(|f| f.outside).unwrap_or_default();
            if outside.is_empty() {
                return Err(Error::invalid(
                    "params.surface is needed (labels, domains, box or plane): the result states no absorbing domains",
                ));
            }
            let inside: Vec<i32> = disc.domain_labels().into_iter().filter(|d| !outside.contains(d)).collect();
            return Ok(SurfaceSpec::Domains { inside });
        }
        Some(s) => s,
    };
    let kind = s.get("kind").and_then(Value::as_str).unwrap_or("");
    match kind {
        "labels" => {
            Ok(SurfaceSpec::Labels { labels: i32_list(s, "labels")?, inside_domains: i32_list(s, "inside_domains")? })
        }
        "domains" => Ok(SurfaceSpec::Domains { inside: i32_list(s, "inside")? }),
        "box" => Ok(SurfaceSpec::Box {
            lo: vec3_param(s, "lo")?.ok_or_else(|| Error::invalid("surface box needs lo"))?,
            hi: vec3_param(s, "hi")?.ok_or_else(|| Error::invalid("surface box needs hi"))?,
            sampling: sampling_param(s)?,
            h: opt_f64(s, "h")?,
        }),
        "plane" => Ok(SurfaceSpec::Plane {
            axis: axis_param(s, "axis", 2)?,
            position: f64_param(s, "position", 0.0)?,
            normal_sign: f64_param(s, "normal_sign", 1.0)?,
            lo: vec2_param(s, "lo")?,
            hi: vec2_param(s, "hi")?,
            sampling: sampling_param(s)?,
            h: opt_f64(s, "h")?,
        }),
        _ => Err(Error::invalid(format!("surface kind '{kind}' is not one of labels, domains, box, plane"))),
    }
}

fn surface_json(s: &SampledSurface, kind: &str) -> Value {
    json!({
        "kind": kind, "sampling": s.sampling, "n_samples": s.samples.len(), "n_elements": s.n_elements,
        "area_m2": s.area, "target_area_m2": s.target_area,
        "coverage": s.target_area.map(|t| s.area / t), "n_outside": s.n_outside
    })
}

fn kind_of(spec: &SurfaceSpec) -> &'static str {
    match spec {
        SurfaceSpec::Labels { .. } => "labels",
        SurfaceSpec::Domains { .. } => "domains",
        SurfaceSpec::Box { .. } => "box",
        SurfaceSpec::Plane { .. } => "plane",
    }
}

fn split(z: &[C64]) -> (Vec<f64>, Vec<f64>) {
    (z.iter().map(|c| c.re).collect(), z.iter().map(|c| c.im).collect())
}

fn rows(v: &[f64], ncol: usize) -> Vec<Vec<f64>> {
    v.chunks(ncol.max(1)).map(|r| r.to_vec()).collect()
}

/// Whether `result` carries 3D per-cell `E` fields (and so the far-field analysis applies).
pub fn farfield_applies(result: &SolveResult) -> bool {
    result.dim == 3
        && result.solutions.iter().any(|s| {
            s.cell_fields.iter().any(|f| f.name == "Ex" || f.name == "Esx")
                && s.cell_fields.iter().any(|f| f.name.starts_with('H'))
        })
}

/// The far-field analysis with JSON parameters (see [`FarFieldAnalysis`]).
pub fn farfield_with_params(mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
    let sel = FieldSelection::from_params(params)?;
    let sampler = CellSampler::new(mesh, result, &sel)?;
    let threads = usize_param(params, "threads", femcore::par::threads())?.max(1);
    let medium = medium_from_params(params, sampler.solution_value())?;
    let spec = surface_from_params(params.get("surface"), result, sampler.disc())?;
    let eq = match params.get("equivalence").and_then(Value::as_str) {
        Some(s) => Equivalence::parse(s)?,
        None if matches!(spec, SurfaceSpec::Plane { .. }) => Equivalence::Magnetic,
        None => Equivalence::Love,
    };
    if eq != Equivalence::Magnetic && !sampler.has_h() {
        return Err(Error::invalid(format!(
            "the {} equivalence needs the H fields {:?}; use equivalence \"magnetic\" on a plane without them",
            eq.label(),
            sel.h_names
        )));
    }
    let surf = sample_surface(&sampler, &spec)?;
    let half = if eq == Equivalence::Love { None } else { surf.half_space_normal };
    let rad = Radiator::new(&surf.samples, &medium, eq, half, threads);
    let p_flux = surface_flux(&surf.samples);
    let pq = params.get("power_quadrature").cloned().unwrap_or(Value::Null);
    let (pole, hemi) = match half {
        Some(n) => (n, true),
        None => ([0.0, 0.0, 1.0], false),
    };
    let p_ff =
        rad.radiated_power(pole, hemi, usize_param(&pq, "n_theta", 48)?.max(2), usize_param(&pq, "n_phi", 96)?.max(4));

    // angular grid
    let theta_deg = match f64_array(params, "theta_deg")? {
        Some(t) => t,
        None => linspace(0.0, 180.0, usize_param(params, "n_theta", 37)?.max(1)),
    };
    let phi_deg = match f64_array(params, "phi_deg")? {
        Some(p) => p,
        None => {
            let n = usize_param(params, "n_phi", 72)?.max(1);
            (0..n).map(|j| 360.0 * j as f64 / n as f64).collect()
        }
    };
    let (nt, np) = (theta_deg.len(), phi_deg.len());
    let mut dirs = Vec::with_capacity(nt * np);
    let mut frames = Vec::with_capacity(nt * np);
    for &t in &theta_deg {
        for &p in &phi_deg {
            let f = spherical_frame(t.to_radians(), p.to_radians());
            dirs.push(f.0);
            frames.push(f);
        }
    }
    let re = rad.far_many(&dirs);
    let u: Vec<f64> = re.iter().map(|e| rad.intensity(e)).collect();
    let e_th: Vec<C64> = re.iter().zip(&frames).map(|(e, f)| rdot(f.1, e)).collect();
    let e_ph: Vec<C64> = re.iter().zip(&frames).map(|(e, f)| rdot(f.2, e)).collect();
    let d: Vec<f64> = u.iter().map(|x| if p_ff > 0.0 { 4.0 * PI * x / p_ff } else { 0.0 }).collect();
    let (imax, dmax) =
        d.iter().copied().enumerate().fold((0, f64::NEG_INFINITY), |a, (i, x)| if x > a.1 { (i, x) } else { a });
    let (eth_re, eth_im) = split(&e_th);
    let (eph_re, eph_im) = split(&e_ph);
    let (en, hn) = sampler.field_names();
    let mut out = json!({
        "analysis": "farfield",
        "solution": sel.solution, "e_fields": en, "h_fields": hn,
        "freq_hz": medium.freq_hz, "eps_r": medium.eps_r, "mu_r": medium.mu_r,
        "wavelength_m": medium.wavelength(), "k_per_m": medium.k(), "eta_ohm": medium.eta(),
        "equivalence": eq.label(), "surface": surface_json(&surf, kind_of(&spec)),
        "P_rad_farfield_w": p_ff, "P_rad_flux_w": p_flux,
        "power_balance_rel": if p_flux != 0.0 { (p_ff - p_flux).abs() / p_flux.abs() } else { f64::NAN },
        "D_max": dmax, "D_max_dbi": 10.0 * dmax.log10(),
        "D_max_theta_deg": theta_deg.get(imax / np.max(1)), "D_max_phi_deg": phi_deg.get(imax % np.max(1)),
        "theta_deg": theta_deg, "phi_deg": phi_deg,
        "U_w_per_sr": rows(&u, np), "directivity": rows(&d, np),
        "r_e_theta_re": rows(&eth_re, np), "r_e_theta_im": rows(&eth_im, np),
        "r_e_phi_re": rows(&eph_re, np), "r_e_phi_im": rows(&eph_im, np),
        "note": "far fields are r*E in volts with exp(-j k r) left out; exp(+j w t) convention"
    });
    let o = out.as_object_mut().ok_or_else(|| Error::invalid("internal: report is not an object"))?;
    if let Some(pin) = opt_f64(params, "input_power_w")? {
        if pin > 0.0 {
            let g: Vec<f64> = u.iter().map(|x| 4.0 * PI * x / pin).collect();
            let gmax = g.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            o.insert("gain".into(), json!(rows(&g, np)));
            o.insert("G_max".into(), json!(gmax));
            o.insert("G_max_dbi".into(), json!(10.0 * gmax.log10()));
            o.insert("radiation_efficiency".into(), json!(p_ff / pin));
        }
    }
    if let Some(einc) = opt_f64(params, "incident_amplitude_v_per_m")? {
        if einc > 0.0 {
            let rcs: Vec<f64> = re.iter().map(|e| 4.0 * PI * cnorm2(e) / (einc * einc)).collect();
            o.insert("rcs_m2".into(), json!(rows(&rcs, np)));
            o.insert("scattering_cross_section_m2".into(), json!(2.0 * medium.eta() * p_ff / (einc * einc)));
        }
    }
    if let Some(ks) = params.get("kspace").filter(|v| !v.is_null()) {
        o.insert("kspace".into(), kspace_json(&rad, ks)?);
    }
    if let Some(fp) = params.get("focal").filter(|v| !v.is_null()) {
        o.insert("focal".into(), focal_json(&rad, &medium, fp, p_flux)?);
    }
    if let Some(st) = params.get("stratified").filter(|v| !v.is_null()) {
        o.insert("stratified".into(), stratified_json(&rad, &medium, st, &phi_deg, nt)?);
    }
    Ok(out)
}

fn kspace_json(rad: &Radiator, p: &Value) -> Result<Value, Error> {
    let axis = axis_param(p, "axis", 2)?;
    let sign = if f64_param(p, "sign", 1.0)? < 0.0 { -1.0 } else { 1.0 };
    let n = usize_param(p, "n", 41)?.max(2);
    let (iu, iv) = uv_axes(axis);
    let us = linspace(-1.0, 1.0, n);
    let mut dirs = Vec::new();
    let mut idx = Vec::new();
    for (j, &uy) in us.iter().enumerate() {
        for (i, &ux) in us.iter().enumerate() {
            let s = ux * ux + uy * uy;
            if s <= 1.0 {
                let mut d = [0.0; 3];
                d[iu] = ux;
                d[iv] = uy;
                d[axis] = sign * (1.0 - s).sqrt();
                dirs.push(d);
                idx.push(j * n + i);
            }
        }
    }
    let f = rad.far_many(&dirs);
    let mut grid = vec![Value::Null; n * n];
    for (e, &i) in f.iter().zip(&idx) {
        grid[i] = json!(rad.intensity(e));
    }
    let rows: Vec<Value> = grid.chunks(n).map(|r| Value::Array(r.to_vec())).collect();
    Ok(json!({
        "axis": AXIS_NAMES[axis], "sign": sign, "ux": us, "uy": us,
        "kx_per_m": us.iter().map(|u| u * rad.k()).collect::<Vec<_>>(), "U_w_per_sr": rows
    }))
}

fn focal_json(rad: &Radiator, medium: &Medium, p: &Value, p_flux: f64) -> Result<Value, Error> {
    let axis = axis_param(p, "axis", 2)?;
    let pos = f64_param(p, "position", 0.0)?;
    let (iu, iv) = uv_axes(axis);
    let us = range_param(p, "u")?.ok_or_else(|| Error::invalid("focal needs u: [min, max, n]"))?;
    let vs = range_param(p, "v")?.unwrap_or_else(|| vec![0.0]);
    let center = vec2_param(p, "center")?.unwrap_or([0.0, 0.0]);
    let mut pts = Vec::with_capacity(us.len() * vs.len());
    for &v in &vs {
        for &u in &us {
            let mut x = [0.0; 3];
            x[axis] = pos;
            x[iu] = center[0] + u;
            x[iv] = center[1] + v;
            pts.push(x);
        }
    }
    let f = rad.near_many(&pts);
    let ea = unit_axis(axis);
    let inten: Vec<f64> = f
        .iter()
        .map(|(e, h)| {
            let s = ccross(e, &cconj(h));
            0.5 * (s[0].re * ea[0] + s[1].re * ea[1] + s[2].re * ea[2]).abs()
        })
        .collect();
    let e2: Vec<f64> = f.iter().map(|(e, _)| cnorm2(e)).collect();
    let nu = us.len();
    let (ipk, peak) =
        inten.iter().copied().enumerate().fold((0, f64::NEG_INFINITY), |a, (i, x)| if x > a.1 { (i, x) } else { a });
    let (pr, pc) = (ipk / nu, ipk % nu);
    // first minimum along +u from the peak, parabolic refinement
    let row = &inten[pr * nu..(pr + 1) * nu];
    let mut first_zero = None;
    let mut c = pc;
    while c + 1 < nu && row[c + 1] <= row[c] {
        c += 1;
    }
    if c > pc && c + 1 < nu {
        let (a, b, d) = (row[c - 1], row[c], row[c + 1]);
        let den = a - 2.0 * b + d;
        let off = if den.abs() > 0.0 { 0.5 * (a - d) / den } else { 0.0 };
        let du = if nu > 1 { us[1] - us[0] } else { 0.0 };
        first_zero = Some(us[c] + off * du - us[pc]);
    }
    let mut out = json!({
        "axis": AXIS_NAMES[axis], "position_m": pos, "center_m": center, "u_m": us, "v_m": vs,
        "intensity_w_per_m2": rows(&inten, nu), "e2_v2_per_m2": rows(&e2, nu),
        "peak_w_per_m2": peak, "peak_u_m": us[pc] + center[0], "peak_v_m": vs[pr] + center[1],
        "measured_first_zero_m": first_zero
    });
    if let Some(na) = opt_f64(p, "na")? {
        if na > 0.0 {
            let pref = f64_param(p, "power_ref_w", p_flux.abs())?;
            let k0na = medium.k0() * na;
            let airy_peak = pref * k0na * k0na / (4.0 * PI);
            if let Some(o) = out.as_object_mut() {
                o.insert("na".into(), json!(na));
                o.insert("power_ref_w".into(), json!(pref));
                o.insert("airy_peak_w_per_m2".into(), json!(airy_peak));
                o.insert("strehl".into(), json!(if airy_peak > 0.0 { peak / airy_peak } else { f64::NAN }));
                o.insert("airy_first_zero_m".into(), json!(AIRY_ZERO / k0na));
            }
        }
    }
    Ok(out)
}

fn stratified_json(rad: &Radiator, medium: &Medium, p: &Value, phi_deg: &[f64], nt: usize) -> Result<Value, Error> {
    let axis = axis_param(p, "axis", 2)?;
    let z0 = f64_param(p, "z0", 0.0)?;
    let eps2 = f64_param(p, "eps2", 1.0)?;
    if eps2 <= 0.0 {
        return Err(Error::invalid("stratified.eps2 must be > 0"));
    }
    let side = match p.get("side").and_then(Value::as_str) {
        None | Some("below") => -1.0,
        Some("above") => 1.0,
        Some(s) => return Err(Error::invalid(format!("stratified.side '{s}' is not below or above"))),
    };
    let ea = scale3(unit_axis(axis), side);
    let (e1, e2) = frame_of(ea);
    let n2 = (eps2 * medium.mu_r).sqrt();
    let eta2 = medium.eta() * medium.index() / n2;
    let ntheta = usize_param(p, "n_theta", nt.div_ceil(2).max(2))?.max(2);
    let theta = linspace(0.0, 90.0, ntheta);
    let np = phi_deg.len();
    let mut u = Vec::with_capacity(ntheta * np);
    for &t in &theta {
        let (st, ct) = t.to_radians().sin_cos();
        for &ph in phi_deg {
            let (sp, cp) = ph.to_radians().sin_cos();
            let d = [0, 1, 2].map(|k| st * (cp * e1[k] + sp * e2[k]) + ct * ea[k]);
            u.push(cnorm2(&transmit_far_field(rad, medium, axis, z0, side, eps2, d)) / (2.0 * eta2));
        }
    }
    let (x, w) = gauss_legendre(usize_param(p, "n_quad", 64)?.max(2));
    let nphi_q = 64;
    let mut ptx = 0.0;
    for (xi, wi) in x.iter().zip(&w) {
        let m = 0.5 * (xi + 1.0);
        let st = (1.0 - m * m).max(0.0).sqrt();
        for j in 0..nphi_q {
            let (sp, cp) = (2.0 * PI * j as f64 / nphi_q as f64).sin_cos();
            let d = [0, 1, 2].map(|k| st * (cp * e1[k] + sp * e2[k]) + m * ea[k]);
            let e = transmit_far_field(rad, medium, axis, z0, side, eps2, d);
            ptx += 0.5 * wi * (2.0 * PI / nphi_q as f64) * cnorm2(&e) / (2.0 * eta2);
        }
    }
    Ok(json!({
        "axis": AXIS_NAMES[axis], "z0_m": z0, "eps2": eps2, "side": if side < 0.0 { "below" } else { "above" },
        "theta_deg": theta, "phi_deg": phi_deg, "U_w_per_sr": rows(&u, np), "P_transmitted_w": ptx,
        "note": "first version: plane-wave Fresnel transmission of the medium-1 far field through an interface outside the surface; \
                 supercritical (evanescent-fed) directions and multiple reflections are left out; theta from the side normal"
    }))
}

/// The Stratton-Chu near-to-far-field analysis behind the `Analysis` trait.
///
/// Parameters: `solution` (0), `field_prefix` ("" or "s"), `e_fields`/`h_fields`, `freq_hz` (solution value), `eps_r`,
/// `mu_r` (1), `surface` ({kind: labels|domains|box|plane, ...}; default: the boundary of the non-absorbing domains),
/// `equivalence` (love; magnetic for a plane), `n_theta` (37), `n_phi` (72) or `theta_deg`/`phi_deg`,
/// `power_quadrature` {n_theta 48, n_phi 96}, `input_power_w`, `incident_amplitude_v_per_m`, `kspace` {axis, sign, n},
/// `focal` {axis, position, u: [min,max,n], v, center, na, power_ref_w}, `stratified` {axis, z0, eps2, side, n_theta},
/// `threads`.
pub struct FarFieldAnalysis;

impl crate::Analysis for FarFieldAnalysis {
    fn id(&self) -> &'static str {
        "farfield"
    }
    fn name(&self) -> &'static str {
        "Near-to-far-field transformation (Stratton-Chu)"
    }
    fn applies_to(&self, result: &SolveResult) -> bool {
        farfield_applies(result)
    }
    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        farfield_with_params(mesh, result, &Value::Null)
    }
    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        farfield_with_params(mesh, result, params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_legendre_integrates_polynomials() {
        let (x, w) = gauss_legendre(7);
        let s: f64 = w.iter().sum();
        assert!((s - 2.0).abs() < 1e-14);
        let i6: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(6)).sum();
        assert!((i6 - 2.0 / 7.0).abs() < 1e-14);
    }

    #[test]
    fn triangle_rule_is_degree_five() {
        // integral of l0^2 l1^3 over the reference triangle (area 1/2) = 2! 3! 2 / 7! * area... : 2!3!/(7!) * 2 * (1/2)
        let exact = 2.0 * 6.0 / 5040.0;
        let q: f64 = TRI_Q.iter().map(|(l, w)| w * l[0].powi(2) * l[1].powi(3)).sum::<f64>();
        assert!((q - exact * 2.0).abs() < 1e-14, "{q} vs {}", exact * 2.0);
    }

    #[test]
    fn par_map_is_thread_independent() {
        let f = |i: usize| (i as f64).sqrt().sin();
        let a = par_map(1000, 1, f);
        let b = par_map(1000, 7, f);
        assert_eq!(a, b);
    }
}
