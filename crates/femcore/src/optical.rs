//! Optical modes of a 2D cross-section as a field-evaluable set: the modes of the `nyx-ems` mode solvers together with the material data
//! needed to evaluate their E, D and H on the mesh, the vacuum constants, the 24-point quadrature rule used for field
//! integrals, and the quadrature rule the polarization analysis (`nyx-analysis`) samples the modes with. Everything works on the native
//! mesh with the finite-element basis functions of the mode solvers. (Split out of `sbs`, which held it next to the Brillouin
//! coupling integrals; the Brillouin code is gone, see ADR 24.) The solvers that produce a set live in `nyx-ems` since restructure
//! phase 4 (ADR 45); the complex material tensors they and the field sampler share (`EpsC`, `MuC`, `cell_materials`) stay here.

use crate::linalg::csr::{c, C};
use crate::edge::{cell_data, edge_basis, tri_g, CellData};
use crate::pml::stretch;
use crate::fem::*;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use std::collections::BTreeSet;

pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const MU0: f64 = 1.256_637_062_12e-6;
/// Speed of light in vacuum (m/s).
pub const C0: f64 = 299_792_458.0;

#[derive(Clone, Copy, Debug)]
pub struct EpsC {
    pub xx: C,
    pub xy: C,
    pub yy: C,
    pub zz: C,
    /// couplings of the transverse components x, y to the propagation axis
    pub xz: C,
    pub yz: C,
}

#[derive(Clone, Copy, Debug)]
pub struct MuC {
    pub x: C,
    pub y: C,
    pub z: C,
}


/// Complex material tensors of every cell including the PML stretching, in the local axes of the mode solvers (`section::local_tensor`).
/// The layers of a cross-section absorb along the wafer axes x and z.
pub fn cell_materials(mesh: &SimplexMesh, m: &Materials) -> Result<(Vec<EpsC>, Vec<MuC>), String> {
    cell_materials_signed(mesh, m, crate::section::PROP_SIGN)
}

/// As `cell_materials` for a propagation axis along `prop_sign * y` (only the reciprocity test uses a sign other than `PROP_SIGN`).
pub fn cell_materials_signed(mesh: &SimplexMesh, m: &Materials, prop_sign: f64) -> Result<(Vec<EpsC>, Vec<MuC>), String> {
    let tens: Vec<_> = m.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?.iter().map(|t| crate::section::local_tensor_signed(t, prop_sign)).collect();
    let mu = m.cell_values(mesh, "mu_r", Some(1.0))?;
    // bounding boxes of the PML domains
    let mut boxes: std::collections::BTreeMap<i32, ([f64; 2], [f64; 2])> = Default::default();
    for c in 0..mesh.n_cells() {
        let l = mesh.cell_labels[c];
        if m.pml.get(&l).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off) {
            let e = boxes.entry(l).or_insert(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]));
            for &v in mesh.cell(c) {
                let p = mesh.point(v as usize);
                for k in 0..2 {
                    e.0[k] = e.0[k].min(p[k]);
                    e.1[k] = e.1[k].max(p[k]);
                }
            }
        }
    }
    let mut eps = Vec::with_capacity(mesh.n_cells());
    let mut mus = Vec::with_capacity(mesh.n_cells());
    for ci in 0..mesh.n_cells() {
        let t = tens[ci];
        let l = mesh.cell_labels[ci];
        let (mut sx, mut sy) = (c(1.0, 0.0), c(1.0, 0.0));
        if let (Some(p), Some((lo, hi))) = (m.pml.get(&l), boxes.get(&l)) {
            let v = mesh.cell(ci);
            let mut cen = [0.0; 2];
            for &i in v {
                for k in 0..2 {
                    cen[k] += mesh.point(i as usize)[k] / 3.0;
                }
            }
            let s = stretch(p, [p.x, p.z], cen, *lo, *hi);
            sx = s.0;
            sy = s.1;
        }
        let (lx, ly, lz) = (sy / sx, sx / sy, sx * sy);
        // eps_ij -> det(S) eps_ij / (s_i s_j) with S = diag(sx, sy, 1): the diagonal gets lx, ly, lz, xy is unchanged,
        // xz and yz are multiplied by sy and sx
        eps.push(EpsC { xx: lx * t[0][0], xy: c(t[0][1], 0.0), yy: ly * t[1][1], zz: lz * t[2][2], xz: sy * t[0][2], yz: sx * t[1][2] });
        mus.push(MuC { x: lx * mu[ci], y: ly * mu[ci], z: lz * mu[ci] });
    }
    Ok((eps, mus))
}

// ------------------------------------------------------------------ optical modes

/// One guided mode: effective index and the edge-element / Lagrange coefficients.
#[derive(Clone, Debug)]
pub struct OptMode {
    pub neff: C,
    pub et: Vec<C>,
    pub phi: Vec<C>,
}

/// The optical modes of a cross-section with the material data needed to evaluate their fields.
#[derive(Clone, Debug)]
pub struct OpticalSet {
    pub freq: f64,
    pub order: usize,
    pub eps: Vec<EpsC>,
    pub mu: Vec<MuC>,
    pub modes: Vec<OptMode>,
    pub converged: bool,
    /// cells inside a PML (left out of power integrals)
    pub in_pml: Vec<bool>,
    /// the local axes (lateral, normal, propagation) of the solver form a left-handed triple (propagation along +y), which reverses the
    /// sense of rotation of the polarization
    pub left_handed: bool,
}

#[derive(Clone, Copy, Debug, Default)]
struct Em {
    e: [C; 3],
}

struct FieldCtx<'a> {
    mesh: &'a SimplexMesh,
    geo: Geometry,
    cd: CellData,
    order: usize,
    omega: f64,
}

impl<'a> FieldCtx<'a> {
    fn new(mesh: &'a SimplexMesh, order: usize, freq: f64) -> Self {
        FieldCtx { mesh, geo: geometry(mesh), cd: cell_data(mesh, order), order, omega: 2.0 * std::f64::consts::PI * freq }
    }

    /// E of `mode` at barycentric point `l` of cell `ci` (envelope field without the `exp(-j beta z)` factor).
    fn at(&self, mode: &OptMode, ci: usize, l: &[f64; 3]) -> Em {
        let nloc_e = if self.order == 1 { 3 } else { 8 };
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        let g = tri_g(&self.geo, ci);
        edge_basis(l, &g, &self.cd.sign[ci], self.order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (mut ex, mut ey) = (c(0.0, 0.0), c(0.0, 0.0));
        for s in 0..nloc_e {
            let cf = mode.et[self.cd.gdofs[ci][s]];
            ex += cf * w[s][0];
            ey += cf * w[s][1];
        }
        let (vals, _) = shape(l, self.order);
        let mut pd = Vec::with_capacity(6);
        self.mesh.cell_dofs(ci, self.order, &mut pd);
        let mut phi = c(0.0, 0.0);
        for a in 0..vals.len() {
            phi += mode.phi[pd[a]] * vals[a];
        }
        let k0 = self.omega / C0;
        let beta = mode.neff * k0;
        let ez = c(0.0, 1.0) * beta * phi;
        Em { e: [ex, ey, ez] }
    }
}

/// The integrals of `|Ex|^2`, `|Ey|^2` and `|Ez|^2` of optical mode `j` over the cells outside the PML (restricted to `domains` if given).
/// In the cross-section x is lateral, y the surface normal and z the propagation axis, so `Ex` dominant means a quasi-TE mode and `Ey` dominant a quasi-TM mode.
pub fn mode_energy_parts(mesh: &SimplexMesh, set: &OpticalSet, j: usize, domains: Option<&BTreeSet<i32>>) -> [f64; 3] {
    let fc = FieldCtx::new(mesh, set.order, set.freq);
    let (pts, wts) = fine_rule();
    let mut acc = [0.0; 3];
    for ci in 0..mesh.n_cells() {
        if set.in_pml[ci] || domains.is_some_and(|d| !d.contains(&mesh.cell_labels[ci])) {
            continue;
        }
        for (l, w) in pts.iter().zip(&wts) {
            let f = fc.at(&set.modes[j], ci, l);
            for k in 0..3 {
                acc[k] += f.e[k].norm_sqr() * w * fc.geo.vol[ci];
            }
        }
    }
    acc
}

/// The electric field `(Ex, Ey, Ez)` of optical mode `j` at the points `pts` (x, y) of the cross-section; `None` outside the mesh.
pub fn sample_mode_e(mesh: &SimplexMesh, set: &OpticalSet, j: usize, pts: &[[f64; 2]]) -> Vec<Option<[C; 3]>> {
    let fc = FieldCtx::new(mesh, set.order, set.freq);
    let tris: Vec<[[f64; 2]; 3]> = (0..mesh.n_cells())
        .map(|ci| {
            let v = mesh.cell(ci);
            let p = |i: usize| {
                let q = mesh.point(v[i] as usize);
                [q[0], q[1]]
            };
            [p(0), p(1), p(2)]
        })
        .collect();
    pts.iter()
        .map(|&[x, y]| {
            for (ci, t) in tris.iter().enumerate() {
                let det = (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[2][0] - t[0][0]) * (t[1][1] - t[0][1]);
                let l1 = ((x - t[0][0]) * (t[2][1] - t[0][1]) - (t[2][0] - t[0][0]) * (y - t[0][1])) / det;
                let l2 = ((t[1][0] - t[0][0]) * (y - t[0][1]) - (x - t[0][0]) * (t[1][1] - t[0][1])) / det;
                let l0 = 1.0 - l1 - l2;
                let eps = -1e-9;
                if l0 >= eps && l1 >= eps && l2 >= eps {
                    return Some(fc.at(&set.modes[j], ci, &[l0, l1, l2]).e);
                }
            }
            None
        })
        .collect()
}

/// The 4 x 6-point rule on the four sub-triangles of a triangle (24 points, degree 4 on each part).
pub fn fine_rule() -> (Vec<[f64; 3]>, Vec<f64>) {
    let q = quad_rule(2, 4);
    let corners = [[[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]], [[0.0, 1.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5]], [[0.0, 0.0, 1.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]];
    let (mut pts, mut ws) = (vec![], vec![]);
    for tri in &corners {
        for (lam, w) in q.lam.iter().zip(&q.w) {
            let mut p = [0.0; 3];
            for k in 0..3 {
                p[k] = lam[0] * tri[0][k] + lam[1] * tri[1][k] + lam[2] * tri[2][k];
            }
            pts.push(p);
            ws.push(w / 4.0);
        }
    }
    (pts, ws)
}

