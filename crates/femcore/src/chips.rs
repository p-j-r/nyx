//! The silicon-on-sapphire rib waveguides of the four ring chips, as Brillouin problems in the
//! wafer frame: x lateral, y the direction of propagation, z the surface normal; the cross-section is the x-z plane.
//!
//! Layers (domain numbers): 1 sapphire, 2 silicon slab, 3 rib, 4 air; the absorbing layers are 11-15 (sapphire: left,
//! right, below, below-left, below-right), 21 and 22 (slab, left and right) and 41-45 (air: left, right, above, above-left,
//! above-right). The physical region is 2 um (3 um for the 500TE chip) wide on each
//! side of the rib, 2 um of sapphire below, 1.2 um of air above, with 1 um absorbing layers around it.
//!
//! The crystals are oriented as follows: the silicon axes are `[p, q, z]` and the sapphire axes
//! `[a, c x a, c]` with `p = (sin psi, -cos psi, 0)` (the projection of the sapphire c axis on the wafer; at `psi = 0` it
//! lies along the waveguide), `q = (cos psi, sin psi, 0)` the in-plane a axis and `c = sin(57.6 deg) p + cos(57.6 deg) z`
//! (r-plane sapphire). The wafer frame is (lateral x, propagation y, normal z), right-handed; a wave along `-y` is as good as one along `+y`: the
//! elastic and optical modes have the same frequencies and gains.
//! The sapphire permittivity is the uniaxial tensor `diag(n_o^2, n_o^2, n_e^2)` of the same crystal.

use crate::material::*;
use crate::mesh::SimplexMesh;
use crate::mesher::*;
use crate::problem::*;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub struct ChipSpec {
    pub name: &'static str,
    /// rib height (total thickness of the silicon) in um
    pub htot: f64,
    /// etch depth in um (the slab is `htot - etch` thick)
    pub etch: f64,
    /// rib width in um
    pub width: f64,
    /// half width of the physical region in um
    pub xin: f64,
    /// effective index of the design mode (the search target for the optical modes)
    pub n_eff: f64,
    pub polarisation: &'static str,
}

pub const CHIPS: [ChipSpec; 4] = [
    ChipSpec { name: "300TE", htot: 0.3, etch: 0.17, width: 0.5, xin: 2.0, n_eff: 2.784, polarisation: "TE" },
    ChipSpec { name: "500TE", htot: 0.5, etch: 0.15, width: 1.0, xin: 3.0, n_eff: 3.2223, polarisation: "TE" },
    ChipSpec { name: "300TM", htot: 0.3, etch: 0.235, width: 0.5, xin: 2.0, n_eff: 2.3677, polarisation: "TM" },
    ChipSpec { name: "500TM", htot: 0.5, etch: 0.325, width: 0.5, xin: 2.0, n_eff: 2.9401, polarisation: "TM" },
];

pub fn chip_by_name(name: &str) -> Option<ChipSpec> {
    CHIPS.iter().copied().find(|c| c.name.eq_ignore_ascii_case(name))
}

/// Optical constants at 1550 nm.
pub const N_SILICON: f64 = 3.476;
pub const N_SAPPHIRE_ORDINARY: f64 = 1.7462;
pub const N_SAPPHIRE_EXTRAORDINARY: f64 = 1.7383;

/// The shapes of the cross-section (lengths in m) and the sizes of the mesh (`scale` multiplies the element sizes).
pub fn chip_shapes(c: &ChipSpec, scale: f64) -> Spec2d {
    let um = 1e-6;
    let (w, htot) = (c.width * um, c.htot * um);
    let hslab = (c.htot - c.etch) * um;
    let (xin, hin, lp, ha) = (c.xin * um, 2.0 * um, 1.0 * um, 1.2 * um);
    let x = xin + lp;
    let (h_si, h_sap, h_pml, h_far) = (0.028 * scale * um, 0.05 * scale * um, 0.07 * scale * um, 0.12 * um);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain: d, h };
    let ytop = htot + ha;
    Spec2d {
        parts: vec![
            rect(-x, -hin - lp, x, ytop + lp, 4, None),
            rect(-x, hslab, -xin, ytop, 41, Some(h_pml)),
            rect(xin, hslab, x, ytop, 42, Some(h_pml)),
            rect(-xin, ytop, xin, ytop + lp, 43, Some(h_pml)),
            rect(-x, ytop, -xin, ytop + lp, 44, Some(h_pml)),
            rect(xin, ytop, x, ytop + lp, 45, Some(h_pml)),
            rect(-xin, -hin, xin, 0.0, 1, Some(h_sap)),
            rect(-x, -hin, -xin, 0.0, 11, Some(h_pml)),
            rect(xin, -hin, x, 0.0, 12, Some(h_pml)),
            rect(-xin, -hin - lp, xin, -hin, 13, Some(h_pml)),
            rect(-x, -hin - lp, -xin, -hin, 14, Some(h_pml)),
            rect(xin, -hin - lp, x, -hin, 15, Some(h_pml)),
            rect(-xin, 0.0, xin, hslab, 2, Some(h_si)),
            rect(-x, 0.0, -xin, hslab, 21, Some(h_pml)),
            rect(xin, 0.0, x, hslab, 22, Some(h_pml)),
            rect(-w / 2.0, hslab, w / 2.0, htot, 3, Some(h_si)),
        ],
        h: h_far,
        min_angle_deg: 25.0,
    }
}

pub const SAPPHIRE_DOMAINS: [i32; 6] = [1, 11, 12, 13, 14, 15];
pub const SILICON_DOMAINS: [i32; 4] = [2, 3, 21, 22];
pub const AIR_DOMAINS: [i32; 6] = [4, 41, 42, 43, 44, 45];

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
}

fn columns(c1: [f64; 3], c2: [f64; 3], c3: [f64; 3]) -> M3 {
    [[c1[0], c2[0], c3[0]], [c1[1], c2[1], c3[1]], [c1[2], c2[2], c3[2]]]
}

/// Orientation matrices (columns = crystal axes in the wafer frame) of silicon and sapphire at the azimuth `psi`.
pub fn crystal_rotations(psi_deg: f64) -> (M3, M3) {
    let (ps, b) = (psi_deg.to_radians(), 57.6f64.to_radians());
    let p = [ps.sin(), -ps.cos(), 0.0];
    let q = [ps.cos(), ps.sin(), 0.0];
    let z = [0.0, 0.0, 1.0];
    let c = [b.sin() * p[0], b.sin() * p[1], b.cos()];
    (columns(p, q, z), columns(q, cross(c, q), c))
}

/// The Brillouin problem of a chip (backward scattering of the design mode with itself): `psi_deg` is the crystal azimuth,
/// `scale` the element size factor, `n_modes` and `f_search` (Hz) set the elastic search.
pub fn chip_problem(c: &ChipSpec, psi_deg: f64, n_modes: usize, f_search: f64, outer_labels: Vec<i32>) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Brillouin);
    s.order = 2;
    s.materials.frame = Frame::Wafer;
    s.inactive = AIR_DOMAINS.to_vec();
    let (r_si, r_sap) = crystal_rotations(psi_deg);
    let (eul_si, eul_sap) = (euler_from_rotation(&r_si), euler_from_rotation(&r_sap));
    for d in 1..=45 {
        s.materials.set(d, "E", 100e9);
        s.materials.set(d, "nu", 0.25);
        s.materials.set(d, "rho", 1.0);
    }
    for d in SAPPHIRE_DOMAINS {
        s.materials.aniso.insert(
            d,
            Anisotropy {
                stiffness: Some(Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 }),
                permittivity: Some(Tensor::Diagonal { x: N_SAPPHIRE_ORDINARY.powi(2), y: N_SAPPHIRE_ORDINARY.powi(2), z: N_SAPPHIRE_EXTRAORDINARY.powi(2) }),
                euler_deg: eul_sap,
                ..Default::default()
            },
        );
        s.materials.set(d, "rho", 3980.0);
    }
    for d in SILICON_DOMAINS {
        s.materials.aniso.insert(d, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: eul_si, ..Default::default() });
        s.materials.set(d, "rho", 2329.0);
    }
    let pml = |x: PmlDir, y: PmlDir| Pml { x, y, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 };
    for (d, pm) in [
        (11, pml(PmlDir::Minus, PmlDir::Off)),
        (12, pml(PmlDir::Plus, PmlDir::Off)),
        (21, pml(PmlDir::Minus, PmlDir::Off)),
        (22, pml(PmlDir::Plus, PmlDir::Off)),
        (13, pml(PmlDir::Off, PmlDir::Minus)),
        (14, pml(PmlDir::Minus, PmlDir::Minus)),
        (15, pml(PmlDir::Plus, PmlDir::Minus)),
    ] {
        s.materials.pml.insert(d, pm);
    }
    s.bcs = vec![Bc { labels: outer_labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s.study = Study::Eigenfrequency { n_modes, f_search };
    let b = s.brillouin.as_mut().unwrap();
    b.wavelength = 1.55e-6;
    b.n_target = c.n_eff;
    b.n_optical_modes = 1;
    b.backward = true;
    b.pump_modes = vec![1];
    b.stokes_modes = vec![1];
    for d in SILICON_DOMAINS {
        b.refractive_index.insert(d, N_SILICON);
    }
    for d in SAPPHIRE_DOMAINS {
        b.refractive_index.insert(d, N_SAPPHIRE_ORDINARY);
    }
    for d in AIR_DOMAINS {
        b.refractive_index.insert(d, 1.0);
    }
    for d in [2, 3] {
        b.photoelastic.insert(d, PhotoSpec { kind: PhotoKind::Cubic { p11: -0.094, p12: 0.017, p44: -0.051 }, euler_deg: [0.0; 3], follow_stiffness: true });
    }
    b.photoelastic.insert(1, PhotoSpec { kind: PhotoKind::Trigonal { p11: -0.23, p12: 0.03, p13: 0.02, p14: 0.0, p33: -0.2, p44: -0.1 }, euler_deg: [0.0; 3], follow_stiffness: true });
    b.moving = vec![2, 3];
    s
}

/// Mesh and problem of a chip.
pub fn chip_model(c: &ChipSpec, psi_deg: f64, scale: f64, n_modes: usize, f_search: f64) -> Result<(SimplexMesh, ProblemSpec), String> {
    let spec2 = chip_shapes(c, scale);
    let mesh = mesh_2d(&spec2)?;
    let um = 1e-6;
    let x = (c.xin + 1.0) * um;
    let hslab = (c.htot - c.etch) * um;
    let ybottom = -3.0 * um;
    let mut labels = vec![];
    for f in 0..mesh.n_facets() {
        let (pa, pb) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
        let (mx, my) = (0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]));
        let outer = (mx.abs() - x).abs() < 1e-12 || (my - ybottom).abs() < 1e-12;
        if outer && my < hslab + 1e-12 {
            labels.push(mesh.facet_labels[f]);
        }
    }
    labels.sort();
    labels.dedup();
    let spec = chip_problem(c, psi_deg, n_modes, f_search, labels);
    Ok((mesh, spec))
}

#[allow(dead_code)]
fn _unused(_: &BTreeMap<i32, f64>) {}
