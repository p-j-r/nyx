//! Optical modes of a 2D cross-section as a field-evaluable set: the modes of `ems`/`ems_pml` together with the material data
//! needed to evaluate their E, D and H on the mesh, the vacuum constants, the 24-point quadrature rule used for field
//! integrals, and the samplers the polarization analysis (`analysis`, `polarization`) is built on. Everything works on the native
//! mesh with the finite-element basis functions of the mode solvers. (Split out of `sbs`, which held it next to the Brillouin
//! coupling integrals; the Brillouin code is gone, see ADR 24.)

use crate::solvers::csr::{c, C};
use crate::ems::{cell_data, edge_basis, tri_g, CellData, EpsCell, Modes, C0};
use crate::ems_pml::{cell_materials, cell_materials_signed, solve_modes_complex, EpsC, ModesC, MuC};
use crate::fem::*;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use std::collections::BTreeSet;

pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const MU0: f64 = 1.256_637_062_12e-6;

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

fn from_real(m: &Modes) -> Vec<OptMode> {
    (0..m.neff.len())
        .map(|j| OptMode { neff: c(m.neff[j].0, m.neff[j].1), et: m.et[j].iter().map(|&x| c(x, 0.0)).collect(), phi: m.phi[j].iter().map(|&x| c(x, 0.0)).collect() })
        .collect()
}

fn from_complex(m: &ModesC) -> Vec<OptMode> {
    (0..m.neff.len()).map(|j| OptMode { neff: m.neff[j], et: m.et[j].clone(), phi: m.phi[j].clone() }).collect()
}

/// Solves the optical modes of `spec` (a `Physics::Emw` problem with a `Modes` study in wafer axes), propagating along `+y`; complex
/// materials and PMLs go through the complex solver.
pub fn solve_optical(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<OpticalSet, String> {
    solve_optical_wafer(mesh, spec, crate::section::PROP_SIGN)
}

/// Like `solve_optical` for a propagation axis along `prop_sign * y` (only the reciprocity test uses a sign other than `PROP_SIGN`).
pub fn solve_optical_wafer(mesh: &SimplexMesh, spec: &ProblemSpec, prop_sign: f64) -> Result<OpticalSet, String> {
    let Study::Modes { freq, n_target, evanescent, n_modes } = spec.study.clone() else {
        return Err("the optical problem needs a Modes study".into());
    };
    if freq <= 0.0 || n_modes < 1 {
        return Err("frequency must be positive and at least one mode is needed".into());
    }
    let (eps, mu) = cell_materials_signed(mesh, &spec.materials, prop_sign)?;
    let in_pml: Vec<bool> = (0..mesh.n_cells()).map(|ci| spec.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off)).collect();
    let coupled = eps.iter().any(|e| e.xz.norm() > 0.0 || e.yz.norm() > 0.0);
    let (modes, converged) = if spec.materials.has_pml() || coupled {
        let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
        let m = solve_modes_complex(mesh, freq, &eps, &mu, target, n_modes, spec.order)?;
        (from_complex(&m), m.converged)
    } else {
        let tens: Vec<_> = spec.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?.iter().map(|t| crate::section::local_tensor_signed(t, prop_sign)).collect();
        let e: Vec<EpsCell> = tens.iter().map(|t| EpsCell { xx: t[0][0], xy: t[0][1], yy: t[1][1], zz: t[2][2] }).collect();
        let mur = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
        let m = crate::ems::solve_modes_tensor(mesh, freq, &e, &mur, n_target, evanescent, n_modes, spec.order)?;
        (from_real(&m), m.converged)
    };
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes, converged, in_pml, left_handed: crate::section::is_left_handed(prop_sign) })
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

/// The electric field of optical mode `j` at the quadrature points of every cell outside the PML, as the samples of a
/// polarization analysis (`polarization::Sample`: `(Ex, Ey, Ez)` with `x` lateral, `y` the surface normal, `z` the propagation axis, and
/// the area each point stands for). The 24-point rule of `fine_rule` is exact for `|E|^2` of the quadratic elements.
pub fn mode_samples(mesh: &SimplexMesh, set: &OpticalSet, j: usize) -> Vec<crate::polarization::Sample> {
    let fc = FieldCtx::new(mesh, set.order, set.freq);
    let (pts, wts) = fine_rule();
    let mut out = Vec::with_capacity(mesh.n_cells() * pts.len());
    for ci in 0..mesh.n_cells() {
        if set.in_pml[ci] {
            continue;
        }
        for (l, w) in pts.iter().zip(&wts) {
            let f = fc.at(&set.modes[j], ci, l);
            let cx = |z: C| [z.re, z.im];
            out.push(crate::polarization::Sample { ex: cx(f.e[0]), ey: cx(f.e[1]), ez: cx(f.e[2]), w: w * fc.geo.vol[ci] });
        }
    }
    out
}

/// The optical modes found by `ems::solve_modes_tensor` for `spec` as an `OpticalSet`, to evaluate their fields.
pub(crate) fn optical_set_real(mesh: &SimplexMesh, spec: &ProblemSpec, m: &Modes, freq: f64) -> Result<OpticalSet, String> {
    let (eps, mu) = cell_materials(mesh, &spec.materials)?;
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes: from_real(m), converged: m.converged, in_pml: vec![false; mesh.n_cells()], left_handed: crate::section::is_left_handed(crate::section::PROP_SIGN) })
}

/// The modes of `ems_pml::solve_modes_complex` for `spec` as an `OpticalSet`.
pub(crate) fn optical_set_complex(mesh: &SimplexMesh, spec: &ProblemSpec, m: &ModesC, freq: f64) -> Result<OpticalSet, String> {
    let (eps, mu) = cell_materials(mesh, &spec.materials)?;
    let in_pml = (0..mesh.n_cells()).map(|ci| spec.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off)).collect();
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes: from_complex(m), converged: m.converged, in_pml, left_handed: crate::section::is_left_handed(crate::section::PROP_SIGN) })
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

