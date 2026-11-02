//! Stimulated Brillouin scattering of a 2D waveguide cross-section: the optical modes (pump and
//! Stokes, from `emw`/`emw_pml`) and the elastic modes with the phase-matched out-of-plane wave
//! number (from `solid_kz`) are combined in the two coupling integrals of the standard theory
//! (Wolff et al., Sturmberg et al.):
//!
//! * the photoelastic overlap `f_PE = -eps0 eps_r^2 * Int conj(E_s) p:S E_p dA`, over the domains
//!   that carry a photoelastic tensor,
//! * the moving-boundary overlap along the interfaces between the body and its surroundings,
//!   `f_MB = Int (u.n) [ (eps_a - eps_b) eps0 conj(E_s)|| E_p|| - (1/eps_a - 1/eps_b) conj(D_s)_n D_p,n / eps0 ] dl`,
//!   with `n` pointing from the body `a` into the surroundings `b`,
//!
//! both divided by the square root of the product of the optical powers, and turned into a gain with
//! the usual conventions: `gain = 2 omega_p (Q / k_eff)
//! |f|^2` with `k_eff = m_eff Omega^2`, power `P = 2 Int Re(conj(Ex) Hy - conj(Ey) Hx)`, and
//! `g0 = sqrt(|vg_p vg_s| hbar omega_p Omega / (4 L) gain / Q)`. For backward scattering the Stokes mode
//! is the mirror image of a forward mode (`Ez`, `Dz` and `Ht` change sign) and `q = beta_p + beta_s`;
//! for forward scattering `q = beta_p - beta_s`.
//!
//! Everything works on the native mesh with the finite-element basis functions of the solvers; the
//! coupling is evaluated once for an optical pair and then applied to every elastic mode.

use crate::csparse::{c, C};
use crate::emw::{cell_data, edge_basis, tri_g, CellData, EpsCell, Modes, C0};
use crate::emw_pml::{cell_materials_signed, solve_modes_complex, EpsC, ModesC, MuC};
use crate::fem::*;
use crate::material::M3;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::solid_kz::KzModes;
use std::collections::{BTreeMap, BTreeSet};

pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const MU0: f64 = 1.256_637_062_12e-6;
pub const HBAR: f64 = 1.054_571_817e-34;

// ------------------------------------------------------------------ photoelastic tensors

/// A photoelastic tensor `p_abkl` in the global frame.
#[derive(Clone, Debug)]
pub struct PhotoElastic {
    pub p: [[[[f64; 3]; 3]; 3]; 3],
}

const VP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];

impl PhotoElastic {
    /// From the 6x6 matrix in standard Voigt order (xx, yy, zz, yz, xz, xy) of the crystal frame and the
    /// rotation `r` whose columns are the crystal axes in the global frame.
    pub fn from_voigt(p6: &[[f64; 6]; 6], r: &M3) -> PhotoElastic {
        let mut t = [[[[0.0; 3]; 3]; 3]; 3];
        for (i, &(a, b)) in VP.iter().enumerate() {
            for (j, &(cc, d)) in VP.iter().enumerate() {
                for (pp, q) in [(a, b), (b, a)] {
                    for (s, u) in [(cc, d), (d, cc)] {
                        t[pp][q][s][u] = p6[i][j];
                    }
                }
            }
        }
        let mut rot = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                for cc in 0..3 {
                                    for d in 0..3 {
                                        s += r[i][a] * r[j][b] * r[k][cc] * r[l][d] * t[a][b][cc][d];
                                    }
                                }
                            }
                        }
                        rot[i][j][k][l] = s;
                    }
                }
            }
        }
        PhotoElastic { p: rot }
    }

    /// Cubic crystal (Voigt constants `p11, p12, p44`).
    pub fn cubic(p11: f64, p12: f64, p44: f64, r: &M3) -> PhotoElastic {
        let mut m = [[0.0; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] = if i == j { p11 } else { p12 };
            }
            m[3 + i][3 + i] = p44;
        }
        PhotoElastic::from_voigt(&m, r)
    }

    /// Trigonal crystal (3m, sapphire), constants as in the literature (`p66 = (p11 - p12) / 2`).
    pub fn trigonal(p11: f64, p12: f64, p13: f64, p14: f64, p33: f64, p44: f64, r: &M3) -> PhotoElastic {
        let p66 = (p11 - p12) / 2.0;
        let m = [
            [p11, p12, p13, p14, 0.0, 0.0],
            [p12, p11, p13, -p14, 0.0, 0.0],
            [p13, p13, p33, 0.0, 0.0, 0.0],
            [p14, -p14, 0.0, p44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, p44, p14],
            [0.0, 0.0, 0.0, 0.0, p14, p66],
        ];
        PhotoElastic::from_voigt(&m, r)
    }

    /// Isotropic material with `p11, p12` (`p44 = (p11 - p12) / 2`).
    pub fn isotropic(p11: f64, p12: f64) -> PhotoElastic {
        PhotoElastic::cubic(p11, p12, (p11 - p12) / 2.0, &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
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

fn from_real(m: &Modes) -> Vec<OptMode> {
    (0..m.neff.len())
        .map(|j| OptMode { neff: c(m.neff[j].0, m.neff[j].1), et: m.et[j].iter().map(|&x| c(x, 0.0)).collect(), phi: m.phi[j].iter().map(|&x| c(x, 0.0)).collect() })
        .collect()
}

fn from_complex(m: &ModesC) -> Vec<OptMode> {
    (0..m.neff.len()).map(|j| OptMode { neff: m.neff[j], et: m.et[j].clone(), phi: m.phi[j].clone() }).collect()
}

/// Solves the optical modes of `spec` (a `Physics::Emw` problem with a `Modes` study) given in the model frame of the Brillouin
/// kernel. Transitional: the Brillouin kernels still work in model axes (phase 5 of the kernel port), while the optical mode solvers
/// read wafer axes, so the problem is turned back into wafer axes here (`frame_adapter::to_wafer_frame`) and the modes propagate along
/// the model's `z = -y` (`section::LEGACY_PROP_SIGN`) until the Brillouin kernels move to `+y` with everything else.
pub fn solve_optical(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<OpticalSet, String> {
    solve_optical_wafer(mesh, &crate::frame_adapter::to_wafer_frame(spec), crate::section::LEGACY_PROP_SIGN)
}

/// Solves the optical modes of `spec` (a `Physics::Emw` problem with a `Modes` study in wafer axes) for a propagation axis along
/// `prop_sign * y`; complex materials and PMLs go through the complex solver.
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
        let m = crate::emw::solve_modes_tensor(mesh, freq, &e, &mur, n_target, evanescent, n_modes, spec.order)?;
        (from_real(&m), m.converged)
    };
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes, converged, in_pml, left_handed: crate::section::is_left_handed(prop_sign) })
}

#[derive(Clone, Copy, Debug, Default)]
struct Em {
    e: [C; 3],
    d: [C; 3],
    h: [C; 3],
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

    /// E, D and H of `mode` at barycentric point `l` of cell `ci` (envelope fields without the
    /// `exp(-j beta z)` factor).
    fn at(&self, set: &OpticalSet, mode: &OptMode, ci: usize, l: &[f64; 3]) -> Em {
        let nloc_e = if self.order == 1 { 3 } else { 8 };
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        let g = tri_g(&self.geo, ci);
        edge_basis(l, &g, &self.cd.sign[ci], self.order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (mut ex, mut ey, mut ct) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        for s in 0..nloc_e {
            let cf = mode.et[self.cd.gdofs[ci][s]];
            ex += cf * w[s][0];
            ey += cf * w[s][1];
            ct += cf * curl[s];
        }
        let (vals, dn) = shape(l, self.order);
        let mut pd = Vec::with_capacity(6);
        self.mesh.cell_dofs(ci, self.order, &mut pd);
        let (mut phi, mut dpx, mut dpy) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        let gb = ci * 6;
        for a in 0..vals.len() {
            let pa = mode.phi[pd[a]];
            phi += pa * vals[a];
            let (mut gx, mut gy) = (0.0, 0.0);
            for b in 0..3 {
                gx += dn[a * 3 + b] * self.geo.g[gb + b * 2];
                gy += dn[a * 3 + b] * self.geo.g[gb + b * 2 + 1];
            }
            dpx += pa * gx;
            dpy += pa * gy;
        }
        let k0 = self.omega / C0;
        let beta = mode.neff * k0;
        let jb = c(0.0, 1.0) * beta;
        let ez = jb * phi;
        let curl_x = jb * dpy + jb * ey;
        let curl_y = -(jb * ex) - jb * dpx;
        let curl_z = ct;
        let e = [ex, ey, ez];
        let m = &set.mu[ci];
        let q = c(0.0, 1.0) / (self.omega * MU0);
        let h = [curl_x * q / m.x, curl_y * q / m.y, curl_z * q / m.z];
        let ep = &set.eps[ci];
        let d = [(ep.xx * ex + ep.xy * ey + ep.xz * ez) * EPS0, (ep.xy * ex + ep.yy * ey + ep.yz * ez) * EPS0, (ep.xz * ex + ep.yz * ey + ep.zz * ez) * EPS0];
        Em { e, d, h }
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
            let f = fc.at(set, &set.modes[j], ci, l);
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
                    return Some(fc.at(set, &set.modes[j], ci, &[l0, l1, l2]).e);
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
            let f = fc.at(set, &set.modes[j], ci, l);
            let cx = |z: C| [z.re, z.im];
            out.push(crate::polarization::Sample { ex: cx(f.e[0]), ey: cx(f.e[1]), ez: cx(f.e[2]), w: w * fc.geo.vol[ci] });
        }
    }
    out
}

/// The optical modes found by `emw::solve_modes_tensor` for `spec` as an `OpticalSet`, to evaluate their fields.
pub(crate) fn optical_set_real(mesh: &SimplexMesh, spec: &ProblemSpec, m: &Modes, freq: f64) -> Result<OpticalSet, String> {
    let prop = crate::section::PROP_SIGN;
    let (eps, mu) = cell_materials_signed(mesh, &spec.materials, prop)?;
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes: from_real(m), converged: m.converged, in_pml: vec![false; mesh.n_cells()], left_handed: crate::section::is_left_handed(prop) })
}

/// The modes of `emw_pml::solve_modes_complex` for `spec` as an `OpticalSet`.
pub(crate) fn optical_set_complex(mesh: &SimplexMesh, spec: &ProblemSpec, m: &ModesC, freq: f64) -> Result<OpticalSet, String> {
    let prop = crate::section::PROP_SIGN;
    let (eps, mu) = cell_materials_signed(mesh, &spec.materials, prop)?;
    let in_pml = (0..mesh.n_cells()).map(|ci| spec.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off)).collect();
    Ok(OpticalSet { freq, order: spec.order, eps, mu, modes: from_complex(m), converged: m.converged, in_pml, left_handed: crate::section::is_left_handed(prop) })
}


/// The relative permittivity of a cell as a 3 x 3 matrix in the order x, y, z (the propagation axis).
pub(crate) fn eps_matrix(e: &EpsC) -> [[C; 3]; 3] {
    [[e.xx, e.xy, e.xz], [e.xy, e.yy, e.yz], [e.xz, e.yz, e.zz]]
}

/// The moving-boundary form for an interface with tangent `t` and normal `n` (in the cross-section plane) between a
/// body `a` and its surroundings `b` (`n` points from `a` to `b`), for anisotropic media. With `E_T = (E.t, E_p)` and
/// `D_n = D.n` continuous across the interface, the energy functional is built from
/// `M = [[eps_TT - eps_Tn eps_nn^-1 eps_nT, eps_Tn eps_nn^-1], [eps_nn^-1 eps_nT, -eps_nn^-1]]`,
/// and the moving-boundary integrand is `conj(Es) (M_a - M_b) Ep` in those variables (for isotropic media this is
/// `eps0 (eps_a - eps_b) conj(Es_T) Ep_T - (1/eps0) (1/eps_a - 1/eps_b) conj(Ds_n) Dp_n`). Returns the blocks
/// `(dA 2x2, dB 2x1, dBt 1x2, dC)` of `M_a - M_b`.
pub(crate) fn mb_blocks(ea: &[[C; 3]; 3], eb: &[[C; 3]; 3], t: [f64; 2], n: [f64; 2]) -> ([[C; 2]; 2], [C; 2], [C; 2], C) {
    let (et, ep, en) = ([t[0], t[1], 0.0], [0.0, 0.0, 1.0], [n[0], n[1], 0.0]);
    let proj = |m: &[[C; 3]; 3], u: &[f64; 3], v: &[f64; 3]| -> C {
        let mut acc = c(0.0, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                acc += m[i][j] * (u[i] * v[j]);
            }
        }
        acc
    };
    let blocks = |m: &[[C; 3]; 3]| -> ([[C; 2]; 2], [C; 2], [C; 2], C) {
        let tt = [[proj(m, &et, &et), proj(m, &et, &ep)], [proj(m, &ep, &et), proj(m, &ep, &ep)]];
        let tn = [proj(m, &et, &en), proj(m, &ep, &en)];
        let nt = [proj(m, &en, &et), proj(m, &en, &ep)];
        let nn = proj(m, &en, &en);
        let inv = c(1.0, 0.0) / nn;
        let mut a = tt;
        for i in 0..2 {
            for j in 0..2 {
                a[i][j] -= tn[i] * nt[j] * inv;
            }
        }
        ([[a[0][0], a[0][1]], [a[1][0], a[1][1]]], [tn[0] * inv, tn[1] * inv], [nt[0] * inv, nt[1] * inv], -inv)
    };
    let (a, b) = (blocks(ea), blocks(eb));
    let mut da = [[c(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            da[i][j] = a.0[i][j] - b.0[i][j];
        }
    }
    (da, [a.1[0] - b.1[0], a.1[1] - b.1[1]], [a.2[0] - b.2[0], a.2[1] - b.2[1]], a.3 - b.3)
}

/// The moving-boundary integrand `conj(Es) (M_a - M_b) Ep` from the fields (E and D as 3-vectors, in-plane
/// components first) on the interface with the blocks of `mb_blocks`.
pub(crate) fn mb_integrand(blocks: &([[C; 2]; 2], [C; 2], [C; 2], C), es_e: &[C; 3], es_d: &[C; 3], ep_e: &[C; 3], ep_d: &[C; 3], t: [f64; 2], n: [f64; 2]) -> C {
    let (da, db, dbt, dc) = blocks;
    let tang = |e: &[C; 3]| [e[0] * t[0] + e[1] * t[1], e[2]];
    let (est, ept) = (tang(es_e), tang(ep_e));
    let (dsn, dpn) = (es_d[0] * n[0] + es_d[1] * n[1], ep_d[0] * n[0] + ep_d[1] * n[1]);
    let mut v = c(0.0, 0.0);
    for i in 0..2 {
        for j in 0..2 {
            v += est[i].conj() * da[i][j] * ept[j] * EPS0;
        }
        v += est[i].conj() * db[i] * dpn;
        v += dsn.conj() * dbt[i] * ept[i];
    }
    v + dsn.conj() * *dc * dpn / EPS0
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

// ------------------------------------------------------------------ coupling samples

#[derive(Clone, Debug)]
struct PeSample {
    dofs: [u32; 6],
    nv: [f64; 6],
    grad: [[f64; 2]; 6],
    /// weight (quadrature weight times cell area)
    w: f64,
    /// `-eps0 eps_r^2 conj(E_s)_a E_p,b p_abkl`, in the order xx, yy, zz, yz, xz, xy with the off-diagonal terms doubled
    a: [C; 6],
}

#[derive(Clone, Debug)]
struct MbSample {
    dofs: [u32; 6],
    nv: [f64; 6],
    w: f64,
    n: [f64; 2],
    m: C,
}

/// How the elastic strain enters the photoelastic overlap.
///
/// With pump `exp(-j beta_p z)`, Stokes `exp(-+j beta_s z)` and phonon `u exp(-j q z)` the Stokes wave is
/// driven by `conj(dEps[u]) E_p`, so the phase-matched overlap is `Int conj(E_s) dEps[conj(S)] E_p`
/// (`Conjugate`, the default; the moving-boundary term uses `conj(u.n)` in the same way). `AsIs` evaluates
/// `conj(E_s) dEps[S] E_p` with the strain as it is (and `conj(u.n)` in the boundary term), the same expression
/// for a phonon of the opposite wave-number sign. The two differ by a large factor whenever the strain
/// components have different phases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StrainConvention {
    /// `conj(E_s) dEps[S] E_p` with the strain as it is
    AsIs,
    /// with the complex conjugate of the strain (the phase-matched overlap)
    Conjugate,
}

#[derive(Clone, Debug)]
pub struct GainOptions {
    /// the mechanical Q used in the gain is `min(Q, q_cap)`
    pub q_cap: f64,
    /// factor in the effective mass
    pub mass_factor: f64,
    /// waveguide length for the single-photon coupling rate (m)
    pub length: f64,
    pub strain: StrainConvention,
    /// use the complex conjugate of the normal displacement in the moving-boundary overlap
    pub mb_conjugate: bool,
    /// additional loss `1/Q` of every elastic mode (for example from `viscous_inverse_q`), added to the
    /// radiation loss that comes from the imaginary part of the frequency
    pub inv_q_extra: Vec<f64>,
    /// classify the modes by their parity under the mirror x -> -x (axis 0) or y -> -y (axis 1) about the middle of
    /// the mesh (for symmetric cross-sections; stays `None` if the mesh is not symmetric)
    pub mirror_axis: Option<usize>,
    /// domains whose elastic energy counts as "body" and as "PML" in the energy fractions
    pub body_domains: BTreeSet<i32>,
    pub pml_domains: BTreeSet<i32>,
}

impl Default for GainOptions {
    fn default() -> Self {
        GainOptions { q_cap: 1000.0, mass_factor: 1.0, length: 0.02, strain: StrainConvention::Conjugate, mb_conjugate: true, inv_q_extra: vec![], mirror_axis: None, body_domains: BTreeSet::new(), pml_domains: BTreeSet::new() }
    }
}

/// What the optical pair contributes to every gain: power-normalised overlap data.
pub struct Coupling {
    pe: Vec<PeSample>,
    mb: Vec<MbSample>,
    pub omega_pump: f64,
    pub omega_stokes: f64,
    pub beta_pump: f64,
    pub beta_stokes: f64,
    pub p_pump: f64,
    pub p_stokes: f64,
    pub vg_pump: f64,
    pub vg_stokes: f64,
    pub backward: bool,
    pub mirror_pump: bool,
    norm: f64,
    order: usize,
}

impl Coupling {
    /// Out-of-plane wave number of the phase-matched elastic wave.
    pub fn q(&self) -> f64 {
        let q = if self.backward { self.beta_pump + self.beta_stokes } else { self.beta_pump - self.beta_stokes };
        if self.mirror_pump { -q } else { q }
    }
    pub fn group_index(&self) -> (f64, f64) {
        (C0 / self.vg_pump.abs(), C0 / self.vg_stokes.abs())
    }
}

/// Inputs of the optical part of the coupling.
pub struct PairInput<'a> {
    pub mesh: &'a SimplexMesh,
    pub pump: (&'a OpticalSet, usize),
    pub stokes: (&'a OpticalSet, usize),
    pub backward: bool,
    /// mirror image of the backward process: the pump is the mirrored (counter-propagating) mode and the
    /// Stokes wave runs forward, so the phonon has the wave number `-(beta_p + beta_s)`; the gains must equal
    /// those of the original process in a waveguide with mirror symmetry along z
    pub mirror_pump: bool,
    /// photoelastic tensors by domain
    pub photoelastic: &'a BTreeMap<i32, PhotoElastic>,
    /// the domains whose boundary moves (the body)
    pub moving: &'a BTreeSet<i32>,
}

/// Evaluates the optical powers, group velocities and the sample points of both overlaps.
pub fn coupling(inp: &PairInput) -> Result<Coupling, String> {
    let mesh = inp.mesh;
    if mesh.dim != 2 {
        return Err("the Brillouin coupling works on a 2D cross-section".into());
    }
    let (set_p, jp) = inp.pump;
    let (set_s, js) = inp.stokes;
    if set_p.order != set_s.order {
        return Err("pump and Stokes modes must use the same element order".into());
    }
    let order = set_p.order;
    let (mode_p, mode_s) = (&set_p.modes[jp], &set_s.modes[js]);
    let fp = FieldCtx::new(mesh, order, set_p.freq);
    let fs = FieldCtx::new(mesh, order, set_s.freq);
    let sgn = if inp.backward { -1.0 } else { 1.0 };
    // the Stokes envelope as seen by the pump: mirrored for backward scattering (or the pump, in the mirror image)
    let flip_with = |mut f: Em, sgn: f64| -> Em {
        f.e[2] *= sgn;
        f.d[2] *= sgn;
        f.h[0] *= sgn;
        f.h[1] *= sgn;
        f
    };
    let flip = |f: Em| -> Em { if inp.mirror_pump { f } else { flip_with(f, sgn) } };
    let flip_pump = |f: Em| -> Em { if inp.mirror_pump { flip_with(f, sgn) } else { f } };
    let (pts, wts) = fine_rule();
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();
    // powers and stored energies over the cells outside the PML
    let (mut pw_p, mut pw_s, mut en_p, mut en_s) = (0.0, 0.0, 0.0, 0.0);
    for ci in 0..nc {
        if set_p.in_pml[ci] {
            continue;
        }
        let area = fp.geo.vol[ci];
        for (l, w) in pts.iter().zip(&wts) {
            let ep = flip_pump(fp.at(set_p, mode_p, ci, l));
            let es = flip(fs.at(set_s, mode_s, ci, l));
            let sz = |f: &Em| (f.e[0].conj() * f.h[1] - f.e[1].conj() * f.h[0]).re;
            pw_p += 2.0 * sz(&ep) * w * area;
            pw_s += 2.0 * sz(&es) * w * area;
            let eps = set_p.eps[ci];
            let en = |f: &Em| (f.e[0].conj() * (eps.xx * f.e[0] + eps.xy * f.e[1] + eps.xz * f.e[2]) + f.e[1].conj() * (eps.xy * f.e[0] + eps.yy * f.e[1] + eps.yz * f.e[2]) + f.e[2].conj() * (eps.xz * f.e[0] + eps.yz * f.e[1] + eps.zz * f.e[2])).re * EPS0;
            en_p += 2.0 * en(&ep) * w * area;
            en_s += 2.0 * en(&es) * w * area;
        }
    }
    if pw_p.abs() < 1e-300 || pw_s.abs() < 1e-300 {
        return Err("an optical mode carries no power (is it evanescent or inside the PML?)".into());
    }
    let norm = 1.0 / (pw_p * pw_s).abs().sqrt();
    // photoelastic samples
    let mut pe = vec![];
    let mut pd = Vec::with_capacity(6);
    for ci in 0..nc {
        let Some(pt) = inp.photoelastic.get(&mesh.cell_labels[ci]) else { continue };
        mesh.cell_dofs(ci, order, &mut pd);
        let area = fp.geo.vol[ci];
        // delta eps = -eps (p:S) eps, so conj(Es) delta eps Ep = -(1/eps0) conj(Ds) (p:S) Dp for any (symmetric) permittivity
        let pref = -1.0 / EPS0;
        for (l, w) in pts.iter().zip(&wts) {
            let ep = flip_pump(fp.at(set_p, mode_p, ci, l));
            let es = flip(fs.at(set_s, mode_s, ci, l));
            // T_kl = sum_ab conj(Es_a) Ep_b p_abkl
            let mut a = [c(0.0, 0.0); 6];
            for (m, &(k, ll)) in VP_ORDER.iter().enumerate() {
                let mut acc = c(0.0, 0.0);
                for aa in 0..3 {
                    for bb in 0..3 {
                        let coef = pt.p[aa][bb][k][ll];
                        if coef != 0.0 {
                            acc += es.d[aa].conj() * ep.d[bb] * coef;
                        }
                    }
                }
                // the shear terms stand for both (k,l) and (l,k)
                a[m] = acc * if k == ll { pref } else { 2.0 * pref };
            }
            let (vals, dn) = shape(l, order);
            let mut sample = PeSample { dofs: [0; 6], nv: [0.0; 6], grad: [[0.0; 2]; 6], w: w * area, a };
            for aa in 0..nloc {
                sample.dofs[aa] = pd[aa] as u32;
                sample.nv[aa] = vals[aa];
                for k in 0..2 {
                    let mut gk = 0.0;
                    for b in 0..3 {
                        gk += dn[aa * 3 + b] * fp.geo.g[ci * 6 + b * 2 + k];
                    }
                    sample.grad[aa][k] = gk;
                }
            }
            pe.push(sample);
        }
    }
    // moving-boundary samples on edges between a moving domain and its surroundings
    let mut edge_cells: Vec<Vec<(u32, u8)>> = vec![vec![]; mesh.n_edges()];
    for ci in 0..nc {
        for (loc, &e) in mesh.cell_edge_ids(ci).iter().enumerate() {
            edge_cells[e as usize].push((ci as u32, loc as u8));
        }
    }
    const LOC: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];
    let gx = [0.5 - 0.5 * 0.861_136_311_594_052_6, 0.5 - 0.5 * 0.339_981_043_584_856_3, 0.5 + 0.5 * 0.339_981_043_584_856_3, 0.5 + 0.5 * 0.861_136_311_594_052_6];
    let gw = [0.5 * 0.347_854_845_137_453_9, 0.5 * 0.652_145_154_862_546_1, 0.5 * 0.652_145_154_862_546_1, 0.5 * 0.347_854_845_137_453_9];
    let mut mb = vec![];
    for cells in &edge_cells {
        if cells.len() != 2 {
            continue;
        }
        let (a, b) = {
            let (c0, c1) = (cells[0], cells[1]);
            let (m0, m1) = (inp.moving.contains(&mesh.cell_labels[c0.0 as usize]), inp.moving.contains(&mesh.cell_labels[c1.0 as usize]));
            match (m0, m1) {
                (true, false) => (c0, c1),
                (false, true) => (c1, c0),
                _ => continue,
            }
        };
        let (ca, cb) = (a.0 as usize, b.0 as usize);
        if set_p.in_pml[ca] || set_p.in_pml[cb] {
            continue;
        }
        let la = LOC[a.1 as usize];
        let va = mesh.cell(ca);
        let vb = mesh.cell(cb);
        let (v0, v1) = (va[la.0], va[la.1]);
        let (p0, p1) = (mesh.point(v0 as usize), mesh.point(v1 as usize));
        let len = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2)).sqrt();
        let t = [(p1[0] - p0[0]) / len, (p1[1] - p0[1]) / len];
        // the local position of the edge vertices in cell b
        let pos_b = |v: u32| vb.iter().position(|&x| x == v).unwrap();
        let (ib0, ib1) = (pos_b(v0), pos_b(v1));
        // normal from a to b
        let cen = |cc: usize| {
            let mut s = [0.0; 2];
            for &v in mesh.cell(cc) {
                for k in 0..2 {
                    s[k] += mesh.point(v as usize)[k] / 3.0;
                }
            }
            s
        };
        let (cen_a, cen_b) = (cen(ca), cen(cb));
        let mut n = [t[1], -t[0]];
        if n[0] * (cen_b[0] - cen_a[0]) + n[1] * (cen_b[1] - cen_a[1]) < 0.0 {
            n = [-n[0], -n[1]];
        }
        let blocks = mb_blocks(&eps_matrix(&set_p.eps[ca]), &eps_matrix(&set_p.eps[cb]), t, n);
        for q in 0..4 {
            let s = gx[q];
            let mut l_a = [0.0; 3];
            l_a[la.0] = 1.0 - s;
            l_a[la.1] = s;
            let mut l_b = [0.0; 3];
            l_b[ib0] = 1.0 - s;
            l_b[ib1] = s;
            let avg = |f: &dyn Fn(usize, &[f64; 3]) -> Em| -> Em {
                let (x, y) = (f(ca, &l_a), f(cb, &l_b));
                let mut r = Em::default();
                for k in 0..3 {
                    r.e[k] = (x.e[k] + y.e[k]) * 0.5;
                    r.d[k] = (x.d[k] + y.d[k]) * 0.5;
                    r.h[k] = (x.h[k] + y.h[k]) * 0.5;
                }
                r
            };
            let ep = flip_pump(avg(&|cc, l| fp.at(set_p, mode_p, cc, l)));
            let es = flip(avg(&|cc, l| fs.at(set_s, mode_s, cc, l)));
            let m = mb_integrand(&blocks, &es.e, &es.d, &ep.e, &ep.d, t, n);
            let (vals, _) = shape(&l_a, order);
            let mut pd_a = Vec::with_capacity(6);
            mesh.cell_dofs(ca, order, &mut pd_a);
            let mut sample = MbSample { dofs: [0; 6], nv: [0.0; 6], w: gw[q] * len, n, m };
            for aa in 0..nloc {
                sample.dofs[aa] = pd_a[aa] as u32;
                sample.nv[aa] = vals[aa];
            }
            mb.push(sample);
        }
    }
    let omega_p = 2.0 * std::f64::consts::PI * set_p.freq;
    let omega_s = 2.0 * std::f64::consts::PI * set_s.freq;
    Ok(Coupling {
        pe,
        mb,
        omega_pump: omega_p,
        omega_stokes: omega_s,
        beta_pump: (mode_p.neff * (omega_p / C0)).re,
        beta_stokes: (mode_s.neff * (omega_s / C0)).re,
        p_pump: pw_p,
        p_stokes: pw_s,
        vg_pump: pw_p / en_p,
        vg_stokes: pw_s / en_s,
        backward: inp.backward,
        mirror_pump: inp.mirror_pump,
        norm,
        order,
    })
}

const VP_ORDER: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];

// ------------------------------------------------------------------ gains

#[derive(Clone, Debug)]
pub struct GainRow {
    pub mode: usize,
    /// complex frequency (Hz)
    pub freq: C,
    pub q_eig: f64,
    pub q_used: f64,
    pub f_mb: C,
    pub f_pe: C,
    pub gain_mb: f64,
    pub gain_pe: f64,
    pub gain_total: f64,
    pub g0_mb: f64,
    pub g0_pe: f64,
    pub g0_total: f64,
    pub frac_body: f64,
    pub frac_pml: f64,
    /// effective mass per length (kg/m) for displacement amplitude 1 (raw eigenvector)
    pub m_eff: f64,
    /// `+1` even, `-1` odd under the mirror of `GainOptions::mirror_axis` (in between for mixed modes)
    pub parity: Option<f64>,
}

/// Gains of every elastic mode in `mech` for the optical pair `cp`.
pub fn gains(mesh: &SimplexMesh, cp: &Coupling, mech: &KzModes, opts: &GainOptions) -> Vec<GainRow> {
    let order = cp.order;
    let n1 = mesh.n_dofs(order);
    let nloc = mesh.nloc(order);
    let kz = mech.kz;
    let (qr, wq) = {
        let q = quad_rule(2, 4);
        (q.lam, q.w)
    };
    let geo = geometry(mesh);
    let vals: Vec<Vec<f64>> = qr.iter().map(|l| shape(l, order).0).collect();
    let nc = mesh.n_cells();
    let mut pd = Vec::with_capacity(6);
    let mut rows = vec![];
    let mirror = opts.mirror_axis.and_then(|ax| mirror_map(mesh, order, ax));
    for (j, fz) in mech.freqs.iter().enumerate() {
        let u = &mech.vectors[j];
        // elastic mass and energy fractions
        let (mut mass, mut body, mut pmlm) = (0.0, 0.0, 0.0);
        for ci in 0..nc {
            let rho = mech.rho[ci];
            if rho <= 0.0 {
                continue;
            }
            mesh.cell_dofs(ci, order, &mut pd);
            let mut cell_mass = 0.0;
            for (q, w) in wq.iter().enumerate() {
                let mut s = 0.0;
                for k in 0..3 {
                    let mut v = c(0.0, 0.0);
                    for a in 0..nloc {
                        v += u[k * n1 + pd[a]] * vals[q][a];
                    }
                    s += v.norm_sqr();
                }
                cell_mass += s * w;
            }
            let m = cell_mass * rho * geo.vol[ci];
            mass += m;
            let lab = mesh.cell_labels[ci];
            if opts.body_domains.contains(&lab) {
                body += m;
            }
            if opts.pml_domains.contains(&lab) {
                pmlm += m;
            }
        }
        let minus_ikz = c(0.0, -kz);
        let mut f_pe = c(0.0, 0.0);
        for s in &cp.pe {
            let (mut uxx, mut uxy, mut uyx, mut uyy, mut uzx, mut uzy) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            let (mut ux, mut uy, mut uz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nloc {
                let d = s.dofs[a] as usize;
                let (ax, ay, az) = (u[d], u[n1 + d], u[2 * n1 + d]);
                uxx += ax * s.grad[a][0];
                uxy += ax * s.grad[a][1];
                uyx += ay * s.grad[a][0];
                uyy += ay * s.grad[a][1];
                uzx += az * s.grad[a][0];
                uzy += az * s.grad[a][1];
                ux += ax * s.nv[a];
                uy += ay * s.nv[a];
                uz += az * s.nv[a];
            }
            // tensor strains in the order xx, yy, zz, yz, xz, xy
            let strain = [uxx, uyy, minus_ikz * uz, (uzy + minus_ikz * uy) * 0.5, (uzx + minus_ikz * ux) * 0.5, (uxy + uyx) * 0.5];
            let mut acc = c(0.0, 0.0);
            for m in 0..6 {
                let e = if opts.strain == StrainConvention::Conjugate { strain[m].conj() } else { strain[m] };
                acc += s.a[m] * e;
            }
            f_pe += acc * s.w;
        }
        let mut f_mb = c(0.0, 0.0);
        for s in &cp.mb {
            let (mut ux, mut uy) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nloc {
                let d = s.dofs[a] as usize;
                ux += u[d] * s.nv[a];
                uy += u[n1 + d] * s.nv[a];
            }
            let un = ux * s.n[0] + uy * s.n[1];
            f_mb += s.m * if opts.mb_conjugate { un.conj() } else { un } * s.w;
        }
        f_pe *= cp.norm;
        f_mb *= cp.norm;
        let omega_m = 2.0 * std::f64::consts::PI * fz.re;
        let inv_rad = if fz.im.abs() > 0.0 { 2.0 * fz.im.abs() / fz.re } else { 0.0 };
        let inv_q = inv_rad + opts.inv_q_extra.get(j).copied().unwrap_or(0.0);
        let q_eig = if inv_q > 0.0 { 1.0 / inv_q } else { f64::INFINITY };
        let q_used = q_eig.min(opts.q_cap);
        let k_eff = opts.mass_factor * mass * omega_m * omega_m;
        let gain = |f: C| 2.0 * cp.omega_pump * (q_used / k_eff) * f.norm_sqr();
        let g0 = |f: C, g: f64| -> f64 {
            let v = (cp.vg_pump * cp.vg_stokes).abs() * (HBAR * cp.omega_pump) * omega_m / (4.0 * opts.length) * (g / q_used);
            v.sqrt() * if f.re > 0.0 { -1.0 } else { 1.0 }
        };
        let (gm, gp, gt) = (gain(f_mb), gain(f_pe), gain(f_pe + f_mb));
        let (g0m, g0p) = (g0(f_mb, gm), g0(f_pe, gp));
        rows.push(GainRow {
            mode: j,
            freq: *fz,
            q_eig,
            q_used,
            f_mb,
            f_pe,
            gain_mb: gm,
            gain_pe: gp,
            gain_total: gt,
            g0_mb: g0m,
            g0_pe: g0p,
            g0_total: (g0m + g0p).abs(),
            frac_body: body / mass.max(1e-300),
            frac_pml: pmlm / mass.max(1e-300),
            m_eff: mass,
            parity: mirror.as_ref().map(|mp| mirror_parity(u, n1, mp, opts.mirror_axis.unwrap_or(0))),
        });
    }
    rows
}

// ------------------------------------------------------------------ symmetry

/// For every vertex the vertex at the mirrored position (about the middle of the mesh along `axis`), or
/// `None` if the vertices are not mirror symmetric (the triangulation may differ: only vertex values are used).
pub fn mirror_map(mesh: &SimplexMesh, _order: usize, axis: usize) -> Option<Vec<usize>> {
    if axis >= mesh.dim {
        return None;
    }
    let x = mesh.dof_coordinates(1);
    let d = mesh.dim;
    let n = mesh.n_points();
    let (lo, hi) = mesh.bbox();
    let mid = 0.5 * (lo[axis] + hi[axis]);
    let extent = (0..d).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let cell = 1e-7 * extent;
    let key = |p: &[f64]| -> Vec<i64> { p.iter().map(|v| (v / cell).round() as i64).collect() };
    let mut table: std::collections::HashMap<Vec<i64>, usize> = std::collections::HashMap::with_capacity(n);
    for i in 0..n {
        table.insert(key(&x[i * d..(i + 1) * d]), i);
    }
    let mut map = Vec::with_capacity(n);
    for i in 0..n {
        let mut p = x[i * d..(i + 1) * d].to_vec();
        p[axis] = 2.0 * mid - p[axis];
        let k = key(&p);
        // allow one cell of rounding
        let found = table.get(&k).copied().or_else(|| {
            for dk in [-1i64, 1] {
                for ax in 0..d {
                    let mut kk = k.clone();
                    kk[ax] += dk;
                    if let Some(&j) = table.get(&kk) {
                        return Some(j);
                    }
                }
            }
            None
        });
        map.push(found?);
    }
    Some(map)
}

/// `Re <u, M u> / <u, u>` over the vertices for the mirror `M` that flips the `axis` component of the
/// displacement and moves the vertices by `map`; the components are stored in blocks of `n1` dofs.
pub fn mirror_parity(u: &[C], n1: usize, map: &[usize], axis: usize) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for k in 0..3 {
        let sign = if k == axis { -1.0 } else { 1.0 };
        for i in 0..map.len() {
            num += sign * (u[k * n1 + i].conj() * u[k * n1 + map[i]]).re;
            den += u[k * n1 + i].norm_sqr();
        }
    }
    num / den.max(1e-300)
}

// ------------------------------------------------------------------ viscous loss

/// Viscosity tensors (Pa s) per cell in the solver's Voigt order from per-domain anisotropic data
/// (`stiffness` holds the viscosity constants in the same classes as the stiffness: cubic, hexagonal,
/// trigonal, orthotropic).
pub fn cell_viscosity(mesh: &SimplexMesh, viscosity: &BTreeMap<i32, crate::material::Anisotropy>) -> Vec<crate::material::M6> {
    (0..mesh.n_cells())
        .map(|ci| match viscosity.get(&mesh.cell_labels[ci]) {
            Some(a) => match &a.stiffness {
                Some(s) => crate::material::stiffness_global(s, a.euler_deg),
                None => [[0.0; 6]; 6],
            },
            None => [[0.0; 6]; 6],
        })
        .collect()
}

/// `1/Q` of every elastic mode from a viscosity tensor, to first order. The stress is `C S + eta dS/dt`,
/// so `(K + j Omega K_eta - Omega^2 M) u = 0` and `Im(Omega) = u^H K_eta u / (2 u^H M u)`, that is
/// `1/Q = 2 Im(Omega) / Omega = u^H K_eta u / (Omega u^H M u)`. For `eta = tau C` this is exactly `Omega tau`.
pub fn viscous_inverse_q(mesh: &SimplexMesh, order: usize, mech: &KzModes, eta: &[crate::material::M6]) -> Vec<f64> {
    let k = crate::solid_kz::assemble_stiffness_kz(mesh, eta, order, mech.kz, None);
    let n1 = mesh.n_dofs(order);
    let nloc = mesh.nloc(order);
    let q = quad_rule(2, 4);
    let vals: Vec<Vec<f64>> = q.lam.iter().map(|l| shape(l, order).0).collect();
    let geo = geometry(mesh);
    let mut pd = Vec::with_capacity(6);
    mech.vectors
        .iter()
        .zip(&mech.freqs)
        .map(|(u, f)| {
            // u^H K_eta u
            let mut ku = vec![c(0.0, 0.0); u.len()];
            for r in 0..k.nrows {
                let mut acc = c(0.0, 0.0);
                for p in k.indptr[r]..k.indptr[r + 1] {
                    acc += k.data[p] * u[k.indices[p]];
                }
                ku[r] = acc;
            }
            let ke: f64 = u.iter().zip(&ku).map(|(a, b)| (a.conj() * *b).re).sum();
            let mut mass = 0.0;
            for ci in 0..mesh.n_cells() {
                let rho = mech.rho[ci];
                if rho <= 0.0 {
                    continue;
                }
                mesh.cell_dofs(ci, order, &mut pd);
                let mut cell = 0.0;
                for (qi, w) in q.w.iter().enumerate() {
                    let mut s = 0.0;
                    for kk in 0..3 {
                        let mut v = c(0.0, 0.0);
                        for a in 0..nloc {
                            v += u[kk * n1 + pd[a]] * vals[qi][a];
                        }
                        s += v.norm_sqr();
                    }
                    cell += s * w;
                }
                mass += cell * rho * geo.vol[ci];
            }
            let omega = 2.0 * std::f64::consts::PI * f.re;
            (ke / (omega * mass).max(1e-300)).max(0.0)
        })
        .collect()
}

// ------------------------------------------------------------------ mode tracking

/// Matches the modes of `next` to those of `prev` by the overlap of their mode shapes
/// (`|<a, b>| / (|a| |b|)` over the displacement coefficients, insensitive to the arbitrary phase).
/// Returns for every mode of `prev` the index of the best partner in `next` (and the overlap), or
/// `None` if the best overlap is below `threshold` or the partner is already taken by a better match.
pub fn track_modes(prev: &KzModes, next: &KzModes, threshold: f64) -> Vec<Option<(usize, f64)>> {
    let norm = |v: &Vec<C>| v.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt().max(1e-300);
    let nprev: Vec<f64> = prev.vectors.iter().map(norm).collect();
    let nnext: Vec<f64> = next.vectors.iter().map(norm).collect();
    let mut cand: Vec<(f64, usize, usize)> = vec![];
    for (i, a) in prev.vectors.iter().enumerate() {
        for (j, b) in next.vectors.iter().enumerate() {
            if a.len() != b.len() {
                continue;
            }
            let mut dot = c(0.0, 0.0);
            for k in 0..a.len() {
                dot += a[k].conj() * b[k];
            }
            cand.push((dot.norm() / (nprev[i] * nnext[j]), i, j));
        }
    }
    cand.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    let (mut used_i, mut used_j) = (vec![false; prev.vectors.len()], vec![false; next.vectors.len()]);
    let mut out = vec![None; prev.vectors.len()];
    for (ov, i, j) in cand {
        if ov < threshold {
            break;
        }
        if !used_i[i] && !used_j[j] {
            used_i[i] = true;
            used_j[j] = true;
            out[i] = Some((j, ov));
        }
    }
    out
}

/// Elastic modes at several wave numbers with the modes followed from one wave number to the next: returns the
/// modes at every `q` and, for every mode at `qs[0]`, the index of the mode that continues it at each `q`
/// (`None` where the overlap with every candidate is below `threshold`).
/// For every mode of one wave number, the index of the matching mode at the next one (if any).
pub type ModeTracks = Vec<Vec<Option<usize>>>;

pub fn elastic_dispersion(mesh: &SimplexMesh, mech: &ProblemSpec, qs: &[f64], threshold: f64) -> Result<(Vec<KzModes>, ModeTracks), String> {
    let mut all: Vec<KzModes> = vec![];
    for &q in qs {
        all.push(crate::solid_kz::solve_solid_kz_full(mesh, mech, q)?.1);
    }
    if all.is_empty() {
        return Ok((all, vec![]));
    }
    let n0 = all[0].vectors.len();
    let mut tracks: Vec<Vec<Option<usize>>> = (0..n0).map(|i| vec![Some(i)]).collect();
    for k in 1..all.len() {
        // current mode of every track, matched against the modes at the next q
        let prev: Vec<Option<usize>> = tracks.iter().map(|t| *t.last().unwrap()).collect();
        let live: Vec<usize> = (0..n0).filter(|&t| prev[t].is_some()).collect();
        let sub = KzModes {
            kz: all[k - 1].kz,
            freqs: live.iter().map(|&t| all[k - 1].freqs[prev[t].unwrap()]).collect(),
            vectors: live.iter().map(|&t| all[k - 1].vectors[prev[t].unwrap()].clone()).collect(),
            c6: vec![],
            rho: vec![],
            converged: true,
        };
        let m = track_modes(&sub, &all[k], threshold);
        for (pos, &t) in live.iter().enumerate() {
            tracks[t].push(m[pos].map(|x| x.0));
        }
        for t in 0..n0 {
            if prev[t].is_none() {
                tracks[t].push(None);
            }
        }
    }
    Ok((all, tracks))
}

// ------------------------------------------------------------------ study

/// A complete Brillouin calculation: optical modes, pairs of pump and Stokes modes, the elastic
/// modes at each phase-matched wave number and all gains.
pub struct BrillouinStudy<'a> {
    pub mesh: &'a SimplexMesh,
    /// optical problem (`Physics::Emw`, a `Modes` study at the pump wavelength; PML allowed)
    pub optics: ProblemSpec,
    /// elastic problem (`Physics::Solid`, an `Eigenfrequency` study whose `f_search` and `n_modes` are used)
    pub mech: ProblemSpec,
    pub photoelastic: BTreeMap<i32, PhotoElastic>,
    pub moving: BTreeSet<i32>,
    pub backward: bool,
    /// indices into the optical modes; empty means every mode
    pub pump_modes: Vec<usize>,
    pub stokes_modes: Vec<usize>,
    pub options: GainOptions,
}

pub struct PairResult {
    pub pump: usize,
    pub stokes: usize,
    pub q: f64,
    pub neff_pump: f64,
    pub neff_stokes: f64,
    pub rows: Vec<GainRow>,
    pub log: String,
}

pub struct StudyResult {
    pub optical: OpticalSet,
    pub pairs: Vec<PairResult>,
}

impl BrillouinStudy<'_> {
    pub fn run(&self) -> Result<StudyResult, String> {
        let set = solve_optical(self.mesh, &self.optics)?;
        let all: Vec<usize> = (0..set.modes.len()).collect();
        let pumps = if self.pump_modes.is_empty() { all.clone() } else { self.pump_modes.clone() };
        let stokes = if self.stokes_modes.is_empty() { all } else { self.stokes_modes.clone() };
        for &j in pumps.iter().chain(&stokes) {
            if j >= set.modes.len() {
                return Err(format!("optical mode {j} does not exist: {} modes were found", set.modes.len()));
            }
        }
        let mut pairs = vec![];
        let mut solved: Vec<(f64, SolveResult, KzModes)> = vec![];
        for &p in &pumps {
            for &s in &stokes {
                let cp = coupling(&PairInput { mesh: self.mesh, pump: (&set, p), stokes: (&set, s), backward: self.backward, mirror_pump: false, photoelastic: &self.photoelastic, moving: &self.moving });
                let cp = match cp {
                    Ok(cp) => cp,
                    Err(e) => {
                        pairs.push(PairResult { pump: p, stokes: s, q: 0.0, neff_pump: set.modes[p].neff.re, neff_stokes: set.modes[s].neff.re, rows: vec![], log: e });
                        continue;
                    }
                };
                let q = cp.q();
                let at = match solved.iter().position(|(k, _, _)| (k - q).abs() <= 1e-9 * q.abs()) {
                    Some(i) => i,
                    None => match crate::solid_kz::solve_solid_kz_full(self.mesh, &self.mech, q) {
                        Ok((r, m)) => {
                            solved.push((q, r, m));
                            solved.len() - 1
                        }
                        Err(e) => {
                            pairs.push(PairResult { pump: p, stokes: s, q, neff_pump: set.modes[p].neff.re, neff_stokes: set.modes[s].neff.re, rows: vec![], log: e });
                            continue;
                        }
                    },
                };
                let (res, modes) = (&solved[at].1, &solved[at].2);
                let rows = gains(self.mesh, &cp, modes, &self.options);
                pairs.push(PairResult { pump: p, stokes: s, q, neff_pump: set.modes[p].neff.re, neff_stokes: set.modes[s].neff.re, rows, log: res.log.join(" | ") });
            }
        }
        Ok(StudyResult { optical: set, pairs })
    }
}

// ------------------------------------------------------------------ the study as a physics

pub fn photo_from_spec(p: &PhotoSpec, stiffness_euler: Option<[f64; 3]>) -> PhotoElastic {
    let r = crate::material::rotation(if p.follow_stiffness { stiffness_euler.unwrap_or(p.euler_deg) } else { p.euler_deg });
    match &p.kind {
        PhotoKind::Isotropic { p11, p12 } => PhotoElastic::cubic(*p11, *p12, (p11 - p12) / 2.0, &r),
        PhotoKind::Cubic { p11, p12, p44 } => PhotoElastic::cubic(*p11, *p12, *p44, &r),
        PhotoKind::Trigonal { p11, p12, p13, p14, p33, p44 } => PhotoElastic::trigonal(*p11, *p12, *p13, *p14, *p33, *p44, &r),
    }
}

/// Mean temperature of every domain from a nodal temperature field (vertex values first).
fn domain_temperatures(mesh: &SimplexMesh, t: &[f64]) -> BTreeMap<i32, f64> {
    let mut sum: BTreeMap<i32, (f64, usize)> = BTreeMap::new();
    for ci in 0..mesh.n_cells() {
        let e = sum.entry(mesh.cell_labels[ci]).or_insert((0.0, 0));
        for &v in mesh.cell(ci) {
            if let Some(x) = t.get(v as usize) {
                e.0 += x;
                e.1 += 1;
            }
        }
    }
    sum.into_iter().filter(|(_, (_, n))| *n > 0).map(|(d, (s, n))| (d, s / n as f64)).collect()
}

/// Adds the table rows and the mode shapes of one pump/Stokes pair to the report.
fn report_rows(tag: &str, rows: &[GainRow], res: &SolveResult, scalars: &mut Vec<(String, Scalar)>, solutions: &mut Vec<(f64, Solution)>, best: &mut Option<(f64, String)>) {
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &c| rows[c].gain_total.partial_cmp(&rows[a].gain_total).unwrap_or(std::cmp::Ordering::Equal));
    for (rank, &i) in order.iter().enumerate() {
        let r = &rows[i];
        let qtxt = if r.q_eig.is_finite() { format!("{:.3e}", r.q_eig) } else { "inf".into() };
        let par = r.parity.map_or(String::new(), |p| format!(" | parity {:+.2}", p));
        let line = format!("{:.5} GHz | Q {} | body {:.2} pml {:.2}{par} | gain MB {:.4e} PE {:.4e} total {:.4e} /W/m | g0 {:.3e} Hz", r.freq.re * 1e-9, qtxt, r.frac_body, r.frac_pml, r.gain_mb, r.gain_pe, r.gain_total, r.g0_total / (2.0 * std::f64::consts::PI));
        if rank < 12 {
            scalars.push((format!("{tag} #{}", rank + 1), Scalar::Text(line)));
        }
        if best.as_ref().is_none_or(|(g, _)| r.gain_total > *g) {
            *best = Some((r.gain_total, format!("{tag}: {:.4} GHz", r.freq.re * 1e-9)));
        }
        let mut sol = res.solutions[i].clone();
        sol.label = format!("{tag}, {:.4} GHz, gain {:.3e} /W/m", r.freq.re * 1e-9, r.gain_total);
        solutions.push((r.gain_total, sol));
    }
}

/// Brillouin study of a ring resonator (`brillouin.ring_radius`): the mesh is the (r, z) half plane.
fn solve_ring(mesh: &SimplexMesh, spec: &ProblemSpec, state: &SolveState, r_ref: f64) -> Result<SolveResult, String> {
    use crate::emw_axisym::*;
    use crate::ring_brillouin::*;
    let t0 = std::time::Instant::now();
    let b = spec.brillouin.as_ref().ok_or("this problem has no Brillouin settings")?;
    let mut log = vec![format!("ring resonator, reference radius {:.4e} m; the crystal data are averaged over rotations about the ring axis", r_ref)];
    let temps = state.heat_t.as_ref().map(|t| domain_temperatures(mesh, t));
    let mut index = BTreeMap::new();
    for d in mesh.domain_labels() {
        let mut n = b.refractive_index.get(&d).copied().unwrap_or(1.0);
        if let (Some(dn), Some(tm)) = (b.dn_dt.get(&d), temps.as_ref().and_then(|t| t.get(&d))) {
            n += dn * (tm - b.t_ref);
        }
        index.insert(d, n);
    }
    let mat = AxiMaterial::from_indices(mesh, &index).with_pml(mesh, &b.optical_pml);
    let f0 = C0 / b.wavelength;
    let k0 = 2.0 * std::f64::consts::PI * f0 / C0;
    let m_pump = (b.n_target * k0 * r_ref).round();
    let modes = solve_axisym_modes(mesh, &mat, m_pump, f0, b.n_optical_modes.max(1), spec.order)?;
    let mut scalars = vec![("azimuthal order of the optical modes".to_string(), Scalar::Num(m_pump))];
    for (j, f) in modes.freqs.iter().enumerate() {
        scalars.push((format!("optical mode {}: frequency (THz) | neff at r_ref", j + 1), Scalar::Text(format!("{:.5} | {:.5}", f * 1e-12, m_pump * C0 / (2.0 * std::f64::consts::PI * f * r_ref)))));
    }
    let pick = |v: &Vec<usize>| -> Result<Vec<usize>, String> {
        if v.is_empty() {
            return Ok((0..modes.freqs.len()).collect());
        }
        v.iter().map(|&m| if m >= 1 && m <= modes.freqs.len() { Ok(m - 1) } else { Err(format!("optical mode {m} does not exist: {} were found", modes.freqs.len())) }).collect()
    };
    let (pumps, stokes) = (pick(&b.pump_modes)?, pick(&b.stokes_modes)?);
    let mut mech = spec.clone();
    mech.physics = Physics::Solid;
    mech.brillouin = None;
    mech.out_of_plane = None;
    let photoelastic: BTreeMap<i32, PhotoElastic> = b.photoelastic.iter().map(|(d, p)| (*d, photo_from_spec(p, spec.materials.aniso.get(d).map(|a| a.euler_deg)).average_about_axis(1, 16))).collect();
    let moving: BTreeSet<i32> = b.moving.iter().copied().collect();
    let pml_domains: BTreeSet<i32> = spec.materials.pml.keys().copied().collect();
    let opts = GainOptions { q_cap: b.q_cap, length: 2.0 * std::f64::consts::PI * r_ref, strain: if b.conjugate_strain { StrainConvention::Conjugate } else { StrainConvention::AsIs }, body_domains: moving.clone(), pml_domains, ..Default::default() };
    let mut solutions = vec![];
    let mut best: Option<(f64, String)> = None;
    let mut solved: Vec<(f64, SolveResult, KzModes)> = vec![];
    for &p in &pumps {
        for &s in &stokes {
            let cp = ring_coupling(&RingInput { mesh, material: &mat, pump: (&modes, p), stokes: (&modes, s), backward: b.backward, photoelastic: &photoelastic, moving: &moving, r_ref })?;
            let m_el = cp.elastic_order();
            let tag = format!("pump {} -> Stokes {}", p + 1, s + 1);
            if m_el.abs() < 0.5 {
                log.push(format!("{tag}: azimuthal order of the phonon is zero; nothing to compute"));
                continue;
            }
            let at = match solved.iter().position(|(k, _, _)| (k - m_el).abs() < 1e-9) {
                Some(i) => i,
                None => {
                    let (r, m) = crate::solid_kz::solve_solid_axisym(mesh, &mech, m_el)?;
                    solved.push((m_el, r, m));
                    solved.len() - 1
                }
            };
            let (res, emodes) = (&solved[at].1, &solved[at].2);
            let rows = ring_gains(mesh, &cp, emodes, &opts);
            scalars.push((format!("{tag}: azimuthal order of the phonon"), Scalar::Num(m_el)));
            let (ngp, ngs) = cp.group_index();
            scalars.push((format!("{tag}: group indices"), Scalar::Text(format!("{ngp:.4} / {ngs:.4}"))));
            report_rows(&tag, &rows, res, &mut scalars, &mut solutions, &mut best);
            log.push(format!("{tag}: M = {m_el}, {} elastic modes; {}", rows.len(), res.log.join(" | ")));
        }
    }
    solutions.sort_by(|a, c| c.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    if let Some((g, what)) = best {
        scalars.insert(0, ("largest total gain (1/W/m)".into(), Scalar::Num(g)));
        scalars.insert(1, ("... for".into(), Scalar::Text(what)));
    }
    let sols: Vec<Solution> = solutions.into_iter().map(|x| x.1).collect();
    if sols.is_empty() {
        return Err("no pair of optical modes gave an elastic problem to solve".into());
    }
    let mut result = SolveResult { order: spec.order, dim: 2, solutions: sols, scalars, log, seconds: t0.elapsed().as_secs_f64(), dofs: 3 * mesh.n_dofs(spec.order), ..Default::default() };
    // polarization of the whispering-gallery modes: (E_r, E_z, E_theta) as (lateral, normal, propagation); neff at the reference radius
    let in_pml: Vec<bool> = (0..mesh.n_cells()).map(|ci| b.optical_pml.contains_key(&mesh.cell_labels[ci])).collect();
    let pol: Vec<crate::analysis::ModeInput> = (0..modes.freqs.len())
        .map(|j| crate::analysis::ModeInput {
            label: format!("mode {}", j + 1),
            neff: (m_pump * C0 / (2.0 * std::f64::consts::PI * modes.freqs[j] * r_ref), 0.0),
            samples: crate::emw_axisym::mode_samples(mesh, &mat, &modes, j, &in_pml),
        })
        .collect();
    let mut opts = crate::analysis::AnalysisOptions::default();
    opts.class.time = crate::polarization::TimeConvention::MinusI;
    crate::analysis::analyze_modes(&pol, &opts).attach(&mut result);
    Ok(result)
}

/// Brillouin study from a `Physics::Brillouin` problem: the optical modes from `brillouin.refractive_index`
/// (with `dn_dt` and a temperature field from a heat study in `state`), the elastic modes from the
/// solid-mechanics data of the problem at the phase-matched wave number of every pump/Stokes pair, and all gains.
pub fn solve_brillouin(mesh: &SimplexMesh, spec: &ProblemSpec, state: &SolveState) -> Result<SolveResult, String> {
    let t0 = std::time::Instant::now();
    let b = spec.brillouin.as_ref().ok_or("this problem has no Brillouin settings")?;
    if mesh.dim != 2 {
        return Err("the Brillouin study works on a 2D cross-section".into());
    }
    if !matches!(spec.study, Study::Eigenfrequency { .. }) {
        return Err("the Brillouin study needs an Eigenfrequency study (number of elastic modes and search frequency)".into());
    }
    if b.wavelength.is_nan() || b.wavelength <= 0.0 {
        return Err("the wavelength must be positive".into());
    }
    if let Some(r_ref) = b.ring_radius {
        if r_ref.is_nan() || r_ref <= 0.0 {
            return Err("the ring radius must be positive".into());
        }
        return solve_ring(mesh, spec, state, r_ref);
    }
    let study = study_straight(mesh, spec, state)?;
    let mut scalars = vec![];
    for (j, z) in study.neff.iter().enumerate() {
        let txt = if z.im.abs() > 1e-9 { format!("{:.6} {:+.3e}i", z.re, z.im) } else { format!("{:.6}", z.re) };
        scalars.push((format!("optical neff {}", j + 1), Scalar::Text(txt)));
    }
    let mut solutions = vec![];
    let mut best: Option<(f64, String)> = None;
    for pr in &study.pairs {
        scalars.push((format!("{}: q (1/m)", pr.tag), Scalar::Num(pr.q)));
        scalars.push((format!("{}: group indices", pr.tag), Scalar::Text(format!("{:.4} / {:.4}", pr.group_index.0, pr.group_index.1))));
        report_rows(&pr.tag, &pr.rows, &pr.res, &mut scalars, &mut solutions, &mut best);
    }
    solutions.sort_by(|a, c| c.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    if let Some((g, what)) = best {
        scalars.insert(0, ("largest total gain (1/W/m)".into(), Scalar::Num(g)));
        scalars.insert(1, ("... for".into(), Scalar::Text(what)));
    }
    let sols: Vec<Solution> = solutions.into_iter().map(|x| x.1).collect();
    if sols.is_empty() {
        return Err("no pair of optical modes gave an elastic problem to solve".into());
    }
    let mut result = SolveResult { order: spec.order, dim: 2, solutions: sols, scalars, log: study.log, seconds: t0.elapsed().as_secs_f64(), dofs: 3 * mesh.n_dofs(spec.order), ..Default::default() };
    if let Some(a) = study.analysis {
        // the analysis describes the optical modes the pump and Stokes waves are taken from
        a.attach(&mut result);
    }
    Ok(result)
}

/// One pump/Stokes pair of a straight-waveguide study: the elastic modes at the phase-matched wave number and their gains.
pub struct PairOutput {
    pub tag: String,
    pub q: f64,
    pub group_index: (f64, f64),
    pub rows: Vec<GainRow>,
    /// the same modes with the opposite strain convention (the strain as it is if `rows` is phase matched, and vice versa)
    pub rows_alt: Vec<GainRow>,
    /// the elastic solution (mode shapes in the order of `rows`)
    pub res: SolveResult,
}

pub struct StraightStudy {
    pub neff: Vec<C>,
    pub pairs: Vec<PairOutput>,
    pub log: Vec<String>,
    /// polarization of the optical modes (pump and Stokes)
    pub analysis: Option<crate::analysis::Analysis>,
}

/// The straight-waveguide Brillouin study of a `Physics::Brillouin` problem (a 2D mesh, `ring_radius` ignored), with the
/// numbers of every elastic mode instead of a report: optical modes, then for every pump/Stokes pair the elastic modes at
/// `q` and their gains.
pub fn study_straight(mesh: &SimplexMesh, spec: &ProblemSpec, state: &SolveState) -> Result<StraightStudy, String> {
    let b = spec.brillouin.as_ref().ok_or("this problem has no Brillouin settings")?;
    // optics
    let mut opt = ProblemSpec::new(Physics::Emw);
    opt.order = spec.order;
    let temps = state.heat_t.as_ref().map(|t| domain_temperatures(mesh, t));
    let mut log = vec![];
    for d in mesh.domain_labels() {
        let mut n = b.refractive_index.get(&d).copied().unwrap_or(1.0);
        if let (Some(dn), Some(tm)) = (b.dn_dt.get(&d), temps.as_ref().and_then(|t| t.get(&d))) {
            n += dn * (tm - b.t_ref);
        }
        opt.materials.set(d, "eps_r", n * n);
        opt.materials.set(d, "mu_r", 1.0);
        // a permittivity tensor of the crystal data of the domain (with its orientation) replaces the refractive index
        if let Some(a) = spec.materials.aniso.get(&d) {
            if a.permittivity.is_some() {
                opt.materials.aniso.insert(d, crate::material::Anisotropy { permittivity: a.permittivity.clone(), euler_deg: a.euler_deg, ..Default::default() });
            }
        }
    }
    if temps.is_some() && !b.dn_dt.is_empty() {
        log.push("refractive indices shifted by dn/dT with the temperature of the heat study".into());
    }
    opt.materials.pml = b.optical_pml.clone();
    let freq = C0 / b.wavelength;
    opt.study = Study::Modes { freq, n_target: b.n_target, evanescent: false, n_modes: b.n_optical_modes.max(1) };
    let set = {
        let _s = crate::progress::stage("optical modes of the cross-section");
        solve_optical(mesh, &opt)?
    };
    if set.modes.is_empty() {
        return Err("no optical mode was found near the target effective index".into());
    }
    crate::progress::note(format!("{} optical mode(s) found, effective indices {}", set.modes.len(), set.modes.iter().map(|m| format!("{:.5}", m.neff.re)).collect::<Vec<_>>().join(", ")));
    let pick = |v: &Vec<usize>| -> Result<Vec<usize>, String> {
        if v.is_empty() {
            return Ok((0..set.modes.len()).collect());
        }
        v.iter().map(|&m| if m >= 1 && m <= set.modes.len() { Ok(m - 1) } else { Err(format!("optical mode {m} does not exist: {} were found", set.modes.len())) }).collect()
    };
    let (pumps, stokes) = (pick(&b.pump_modes)?, pick(&b.stokes_modes)?);
    // elastic problem
    let mut mech = spec.clone();
    mech.physics = Physics::Solid;
    mech.brillouin = None;
    mech.out_of_plane = None;
    let photoelastic: BTreeMap<i32, PhotoElastic> = b.photoelastic.iter().map(|(d, p)| (*d, photo_from_spec(p, spec.materials.aniso.get(d).map(|a| a.euler_deg)))).collect();
    let moving: BTreeSet<i32> = b.moving.iter().copied().collect();
    let pml_domains: BTreeSet<i32> = spec.materials.pml.keys().copied().collect();
    let mut opts = GainOptions { q_cap: b.q_cap, length: b.length, strain: if b.conjugate_strain { StrainConvention::Conjugate } else { StrainConvention::AsIs }, body_domains: moving.clone(), pml_domains, mirror_axis: b.mirror_axis, ..Default::default() };
    let eta = if b.viscosity.is_empty() { None } else { Some(cell_viscosity(mesh, &b.viscosity)) };
    if !set.converged {
        log.push("warning: the optical eigensolver did not converge every mode".into());
    }
    let mut pairs = vec![];
    // elastic problems already solved, by wave number (pairs with the same q share them)
    let mut solved: Vec<(f64, SolveResult, KzModes)> = vec![];
    for &p in &pumps {
        for &s in &stokes {
            crate::progress::check()?;
            let tag = format!("pump {} -> Stokes {}", p + 1, s + 1);
            let _pair_stage = crate::progress::stage(format!("{tag}: phase matching, elastic modes and gains"));
            let cp = coupling(&PairInput { mesh, pump: (&set, p), stokes: (&set, s), backward: b.backward, mirror_pump: false, photoelastic: &photoelastic, moving: &moving })?;
            let q = cp.q();
            if q.abs() < 1.0 {
                log.push(format!("{tag}: the phase-matched wave number is zero (forward scattering of a mode with itself); nothing to compute"));
                continue;
            }
            let at = match solved.iter().position(|(k, _, _)| (k - q).abs() <= 1e-9 * q.abs()) {
                Some(i) => i,
                None => {
                    let (r, m) = crate::solid_kz::solve_solid_kz_full(mesh, &mech, q)?;
                    solved.push((q, r, m));
                    solved.len() - 1
                }
            };
            let (res, modes) = (&solved[at].1, &solved[at].2);
            if let Some(eta) = &eta {
                opts.inv_q_extra = viscous_inverse_q(mesh, spec.order, modes, eta);
            }
            let rows = gains(mesh, &cp, modes, &opts);
            let alt = GainOptions { strain: if opts.strain == StrainConvention::Conjugate { StrainConvention::AsIs } else { StrainConvention::Conjugate }, ..opts.clone() };
            let rows_alt = gains(mesh, &cp, modes, &alt);
            log.push(format!("{tag}: q = {q:.4e} /m, {} elastic modes; {}", rows.len(), res.log.join(" | ")));
            pairs.push(PairOutput { tag, q, group_index: cp.group_index(), rows, rows_alt, res: res.clone() });
        }
    }
    let analysis = Some(crate::analysis::analyze_optical(mesh, &set, &Default::default()));
    Ok(StraightStudy { neff: set.modes.iter().map(|m| m.neff).collect(), pairs, log, analysis })
}
