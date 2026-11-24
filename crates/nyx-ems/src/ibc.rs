// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Kinetic-inductance surface-impedance boundary conditions (`surface_impedance`, KICQ IBC), explicit interior/exterior
//! `pec` sheets, rank-1 and distributed-sheet lumped Josephson junction (`lumped_jj`) operators, single-facet deduplication,
//! and symmetric Ruiz diagonal equilibration for 3D `cavity` and 2D `modes` / `axisymmetric_modes`.

use crate::cavity::{basis3, safe_curved_jacobian, Ned3};
use crate::modes::C0;
use crate::problem::{LumpedJjBc, SurfaceImpedanceBc};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{facet_rule, geometry, quad_rule, shape, CurvedTet};
use femcore::linalg::csr::{c, CCoo, CCsr, Coo, RCsr, C};
use nyx_physics::Discretisation;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;

/// Vacuum magnetic permeability `mu_0` in H/m.
pub const MU0: f64 = 4.0 * PI * 1e-7;
/// Vacuum electric permittivity `eps_0 = 1 / (mu_0 c_0^2)` in F/m, consistent with `MU0` and `C0` to machine precision.
pub const EPS0: f64 = 1.0 / (MU0 * C0 * C0);
/// Vacuum wave impedance `eta_0 = mu_0 c_0` in Ohm.
pub const ETA0: f64 = MU0 * C0;
/// Superconducting magnetic flux quantum `Phi_0 = h / (2e)` in Wb.
pub const PHI0: f64 = 2.067_833_848_461_929_5e-15;
/// Planck constant `h` in J*s.
pub const PLANCK_H: f64 = 6.626_070_15e-34;

const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Conversion scale `(Phi_0 / (2 pi))^2 / h` relating `L_J` (H) and `E_J / h` (Hz): `L_J * E_J_Hz = jj_flux_scale()`.
#[inline]
pub fn jj_flux_scale() -> f64 {
    let phi0_2pi = PHI0 / (2.0 * PI);
    (phi0_2pi * phi0_2pi) / PLANCK_H
}

/// Resolves `(L_j, E_j_Hz)` when one is given and the other is `<= 0.0`.
#[inline]
pub fn resolve_jj_l_and_ej(l_j: f64, e_j_hz: f64) -> (f64, f64) {
    let scale = jj_flux_scale();
    if l_j <= 0.0 && e_j_hz > 0.0 {
        (scale / e_j_hz, e_j_hz)
    } else if l_j > 0.0 && e_j_hz <= 0.0 {
        (l_j, scale / l_j)
    } else {
        (l_j, e_j_hz)
    }
}

/// Total inverse linear inductance `1 / L_tot = 1 / L_J + 1 / L_p` (1/H) of a lumped JJ port.
#[inline]
pub fn jj_inv_l_tot(l_j: f64, e_l_hz: f64) -> f64 {
    let inv_lj = if l_j > 0.0 { 1.0 / l_j } else { 0.0 };
    let inv_lp = if e_l_hz > 0.0 { e_l_hz / jj_flux_scale() } else { 0.0 };
    inv_lj + inv_lp
}

/// Series-to-parallel conversion of `Z_s(omega_ref) = R_s + i omega_ref L_k_sq`:
/// returns `(1 / L_eff, G_eff)` where `L_eff = L_k_sq + R_s^2 / (omega_ref^2 L_k_sq)` and
/// `G_eff = R_s / (R_s^2 + omega_ref^2 L_k_sq^2)`.
#[inline]
pub fn ibc_parallel_params(r_s: f64, l_k_sq: f64, omega_ref: f64) -> (f64, f64) {
    if l_k_sq > 0.0 {
        let wl = omega_ref * l_k_sq;
        let denom = r_s * r_s + wl * wl;
        let inv_l_eff = if denom > 0.0 { (omega_ref * wl) / denom } else { 1.0 / l_k_sq };
        let g_eff = if denom > 0.0 { r_s / denom } else { 0.0 };
        (inv_l_eff, g_eff)
    } else if r_s > 0.0 {
        (0.0, 1.0 / r_s)
    } else {
        (0.0, 0.0)
    }
}

#[inline]
pub fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn pt3(mesh: &Discretisation, v: u32) -> [f64; 3] {
    let p = mesh.point(v as usize);
    [p[0], p[1], p[2]]
}

/// Assembled operators and geometric metadata for one `lumped_jj` terminal port.
#[derive(Clone, Debug)]
pub struct AssembledJjPort {
    pub bc: LumpedJjBc,
    pub d_hat: [f64; 3],
    pub w_hat: [f64; 3],
    pub n_hat: [f64; 3],
    pub length_m: f64,
    pub width_m: f64,
    pub u_min: f64,
    pub u_max: f64,
    pub inv_l_tot: f64,
    pub triangles: Vec<[u32; 3]>,
    /// Terminal voltage functional `b_j` on all full Nedelec DOFs (`n`).
    pub b_full: Vec<f64>,
    /// Terminal voltage functional `b_j` on the reduced free Nedelec DOFs (`nred`).
    pub b_red: Vec<f64>,
    /// Reduced stiffness contribution `K_JJ,j` (`nred x nred`).
    pub k_jj_red: RCsr,
    /// Reduced mass contribution `M_JJ,j` (`nred x nred`).
    pub m_jj_red: RCsr,
}

/// Sheet surface-current evaluator rule attached to a 3D cell face.
#[derive(Clone, Copy, Debug)]
pub enum SheetCellKind {
    Ibc { inv_l_eff: f64, g_eff: f64, n_hat: [f64; 3] },
    Jj { port_index: usize, d_hat: [f64; 3], width_m: f64 },
    Pec { n_hat: [f64; 3] },
}

/// One boundary triangle handled by perturbative surface-current IBC (`ibc_mode = "perturbative"` or `"auto"`, FR-30).
#[derive(Clone, Copy, Debug)]
pub struct PerturbativeIbcTri {
    pub tri: [u32; 3],
    pub l_k_sq: f64,
    pub r_s: f64,
    pub f_ref: f64,
}

/// Combined 3D assembly of `pec`, `surface_impedance`, and `lumped_jj` boundary conditions.
pub struct IbcJjAssembly3D {
    /// Unique sorted vertex triples of explicit `pec` facets.
    pub pec_tris: Vec<[u32; 3]>,
    /// Unique sorted vertex triples of direct `surface_impedance` facets (Nedelec DOFs free).
    pub ibc_tris: Vec<[u32; 3]>,
    /// Unique sorted vertex triples of perturbative `surface_impedance` facets (strong PEC in Stage 1/2).
    pub pert_ibc_tris: Vec<PerturbativeIbcTri>,
    /// Boundary labels belonging to direct (non-perturbative) `surface_impedance` conditions.
    pub direct_ibc_labels: Vec<i32>,
    /// Unique sorted vertex triples of all `lumped_jj` facets.
    pub jj_tris: Vec<[u32; 3]>,
    /// Full-DOF real kinetic inductance matrix `K_IBC` (`n x n`).
    pub k_ibc_full: RCsr,
    /// Full-DOF real surface loss matrix `C_IBC` (`n x n`).
    pub c_ibc_full: RCsr,
    /// Full-DOF complex IBC contribution `K_IBC + i k_{0,ref} C_IBC` (`n x n`).
    pub kc_ibc_full: CCsr,
    /// Full-DOF real JJ inductive matrix `sum_j K_JJ,j` (`n x n`).
    pub k_jj_full: RCsr,
    /// Full-DOF real JJ capacitive mass matrix `sum_j M_JJ,j` (`n x n`).
    pub m_jj_full: RCsr,
    /// Per-port assembled JJ data.
    pub jj_ports: Vec<AssembledJjPort>,
    /// Per cell `ci`, optional sheet face rule for exporting `Kx, Ky, Kz` in `"A/m"`.
    pub cell_sheet: Vec<Option<SheetCellKind>>,
    /// Diagnostic warnings generated during IBC assembly (e.g. coarse P1 direct IBC resolution).
    pub warnings: Vec<String>,
}

pub use femcore::fem::build_face_to_cell_3d;

/// Maps each sorted 3D triangle `[v0, v1, v2]` to all adjacent tetrahedra `(ci, lf)` (1 for exterior faces,
/// 2 for interior sheet faces), with `lf` ordered so `mesh.cell(ci)[lf[k]] == tri[k]` for `k = 0..3`.
pub fn build_face_to_cells_3d(mesh: &Discretisation) -> HashMap<[u32; 3], Vec<(usize, [usize; 3])>> {
    let mut map: HashMap<[u32; 3], Vec<(usize, [usize; 3])>> = HashMap::with_capacity(mesh.n_cells() * 2);
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for f in FACES {
            let key = sorted3(v[f[0]], v[f[1]], v[f[2]]);
            let mut lf = [0usize; 3];
            for (k, &tv) in key.iter().enumerate() {
                lf[k] = v.iter().position(|&cv| cv == tv).expect("face vertex belongs to cell");
            }
            map.entry(key).or_default().push((ci, lf));
        }
    }
    map
}

/// Precomputed quadrature-point surface-current basis on perturbative `surface_impedance` triangles (FR-30).
/// Each quadrature point evaluates `K_tilde = n_hat x [mu_r^-1 curl E]` (jump across both sides on interior
/// 2-sided sheets, single-sided `n_in x (mu_r^-1 curl E)` on exterior boundary faces).
pub struct PerturbativeIbcEvaluator {
    /// Per quadrature point across all perturbative triangles: `(weight * l_k_sq / mu_0, weight * r_s, dof_vecs)`.
    points: Vec<(f64, f64, Vec<(usize, [f64; 3])>)>,
}

impl PerturbativeIbcEvaluator {
    pub fn new(
        mesh: &Discretisation,
        ned: &Ned3,
        geo: &femcore::fem::Geometry,
        mur_phys: &[f64],
        pert_tris: &[PerturbativeIbcTri],
    ) -> Result<Self, String> {
        let face_to_cells = build_face_to_cells_3d(mesh);
        let quad2 = quad_rule(2, 4);
        let nq = quad2.w.len();
        let nb = ned.nb();
        let mut points = Vec::with_capacity(pert_tris.len() * nq);
        let (mut w, mut cu) = (Vec::new(), Vec::new());

        for pt in pert_tris {
            let adj = face_to_cells
                .get(&pt.tri)
                .ok_or_else(|| format!("perturbative IBC facet {:?} is not a face of any tetrahedron", pt.tri))?;
            let (p0, p1, p2) = (pt3(mesh, pt.tri[0]), pt3(mesh, pt.tri[1]), pt3(mesh, pt.tri[2]));
            let n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
            let n_len = norm3(n_raw);
            let area = 0.5 * n_len;
            if area <= 0.0 {
                continue;
            }
            let n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];

            // Precompute per-adjacent-cell geometry & orientation sign (`sign * n_hat` points INTO cell `ci`)
            struct AdjCell {
                ci: usize,
                lf: [usize; 3],
                sign_inv_mur: f64,
                gid: [u32; 4],
                gs: [[f64; 3]; 4],
                curved: Option<CurvedTet>,
            }
            let mut adj_cells = Vec::with_capacity(adj.len());
            for &(ci, lf) in adj {
                let v = mesh.cell(ci);
                let opp = 6 - lf[0] - lf[1] - lf[2];
                let p_opp = pt3(mesh, v[opp]);
                let sign = if dot3(sub3(p_opp, p0), n_hat) >= 0.0 { 1.0 } else { -1.0 };
                let gid = ned.keys(mesh, ci);
                let mut gs = [[0.0; 3]; 4];
                for b in 0..4 {
                    for k in 0..3 {
                        gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                    }
                }
                let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
                adj_cells.push(AdjCell { ci, lf, sign_inv_mur: sign / mur_phys[ci], gid, gs, curved });
            }

            let lk_over_mu0 = pt.l_k_sq / MU0;
            for q in 0..nq {
                let l2 = &quad2.lam[q];
                let mut area_q = area;
                let mut n_q = n_hat;
                let mut dof_map: BTreeMap<usize, [f64; 3]> = BTreeMap::new();
                for (ac_idx, ac) in adj_cells.iter().enumerate() {
                    let mut lam4 = [0.0; 4];
                    lam4[ac.lf[0]] = l2[0];
                    lam4[ac.lf[1]] = l2[1];
                    lam4[ac.lf[2]] = l2[2];
                    match &ac.curved {
                        Some(ct) => {
                            let gq = safe_curved_jacobian(ct, &lam4, &ac.gs);
                            basis3(ned.order, &lam4, &gq, &ac.gid, &mut w, &mut cu);
                            if ac_idx == 0 {
                                let d = ct.dx_dlam(&lam4);
                                let nr = cross3(sub3(d[ac.lf[1]], d[ac.lf[0]]), sub3(d[ac.lf[2]], d[ac.lf[0]]));
                                let nl = norm3(nr);
                                if nl > 1e-30 {
                                    area_q = 0.5 * nl;
                                    n_q = [nr[0] / nl, nr[1] / nl, nr[2] / nl];
                                    if dot3(n_q, n_hat) < 0.0 {
                                        n_q = [-n_q[0], -n_q[1], -n_q[2]];
                                    }
                                }
                            }
                        }
                        None => basis3(ned.order, &lam4, &ac.gs, &ac.gid, &mut w, &mut cu),
                    }
                    let cdofs = &ned.cell_dofs[ac.ci];
                    for s in 0..nb {
                        let kx = cross3(n_q, cu[s]);
                        let entry = dof_map.entry(cdofs[s]).or_insert([0.0; 3]);
                        entry[0] += ac.sign_inv_mur * kx[0];
                        entry[1] += ac.sign_inv_mur * kx[1];
                        entry[2] += ac.sign_inv_mur * kx[2];
                    }
                }
                let wq = quad2.w[q] * area_q;
                let dof_vecs: Vec<(usize, [f64; 3])> = dof_map.into_iter().collect();
                points.push((wq * lk_over_mu0, wq * pt.r_s, dof_vecs));
            }
        }

        Ok(Self { points })
    }

    /// Evaluates `(sum_T (L_k_sq / mu_0) int_T |K_tilde|^2 dS, sum_T R_s int_T |K_tilde|^2 dS)`
    /// for a full-DOF complex electric field vector `x_full`.
    pub fn eval_mode(&self, x_full: &[C]) -> (f64, f64) {
        let mut sum_lk = 0.0_f64;
        let mut sum_rs = 0.0_f64;
        for (w_lk, w_rs, dof_vecs) in &self.points {
            let mut kt = [c(0.0, 0.0); 3];
            for &(d, kv) in dof_vecs {
                let xd = x_full[d];
                kt[0] += xd * kv[0];
                kt[1] += xd * kv[1];
                kt[2] += xd * kv[2];
            }
            let norm2 = kt[0].norm_sqr() + kt[1].norm_sqr() + kt[2].norm_sqr();
            sum_lk += w_lk * norm2;
            sum_rs += w_rs * norm2;
        }
        (sum_lk, sum_rs)
    }

    /// Evaluates the symmetric `r x r` Gram matrices `(dk_kin, g_r_mat)` on full-DOF real basis columns `q_full_cols`:
    /// - `dk_kin[a][b] = sum_T (L_k_sq / mu_0) int_T K_tilde_a . K_tilde_b dS`
    /// - `g_r_mat[a][b] = sum_T R_s int_T K_tilde_a . K_tilde_b dS`
    pub fn eval_rom_grams(&self, q_full_cols: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let r = q_full_cols.len();
        let mut dk_kin = vec![vec![0.0_f64; r]; r];
        let mut g_r = vec![vec![0.0_f64; r]; r];
        if r == 0 || self.points.is_empty() {
            return (dk_kin, g_r);
        }
        let mut kts = vec![[0.0_f64; 3]; r];
        for (w_lk, w_rs, dof_vecs) in &self.points {
            if *w_lk == 0.0 && *w_rs == 0.0 {
                continue;
            }
            for (a, q_full) in q_full_cols.iter().enumerate() {
                let mut kt = [0.0_f64; 3];
                for &(d, kv) in dof_vecs {
                    let qd = q_full[d];
                    kt[0] += qd * kv[0];
                    kt[1] += qd * kv[1];
                    kt[2] += qd * kv[2];
                }
                kts[a] = kt;
            }
            for a in 0..r {
                for b in a..r {
                    let dp = dot3(kts[a], kts[b]);
                    let dlk = w_lk * dp;
                    let drs = w_rs * dp;
                    dk_kin[a][b] += dlk;
                    g_r[a][b] += drs;
                    if a != b {
                        dk_kin[b][a] += dlk;
                        g_r[b][a] += drs;
                    }
                }
            }
        }
        (dk_kin, g_r)
    }
}

/// Evaluates the tangential Nedelec basis vectors `psi_{s,||} = w_s - (w_s . n) n` at the 6-point Dunavant
/// quadrature points on a single triangle `tri_sorted = [v0, v1, v2]`.
/// Returns `(area, n_hat, bnd_dofs, quad_weights, psi_parallel_at_q)` where `psi_parallel_at_q[q][k]` is the
/// tangential basis vector of `bnd_dofs[k]` at quadrature point `q`.
fn eval_triangle_tangential_basis(
    mesh: &Discretisation,
    ned: &Ned3,
    geo: &femcore::fem::Geometry,
    face_to_cell: &HashMap<[u32; 3], (usize, [usize; 3])>,
    tri_sorted: [u32; 3],
) -> Result<(f64, [f64; 3], Vec<usize>, Vec<f64>, Vec<Vec<[f64; 3]>>), String> {
    let &(ci, lf) = face_to_cell
        .get(&tri_sorted)
        .ok_or_else(|| format!("facet {:?} is not a face of any tetrahedron in the mesh", tri_sorted))?;
    let v = mesh.cell(ci);
    let (p0, p1, p2) = (pt3(mesh, v[lf[0]]), pt3(mesh, v[lf[1]]), pt3(mesh, v[lf[2]]));
    let n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
    let n_len = norm3(n_raw);
    let area = 0.5 * n_len;
    if area <= 0.0 {
        return Err("degenerate zero-area boundary triangle".into());
    }
    let n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];

    let bnd_dofs = ned.boundary_dofs(mesh, tri_sorted);
    let cell_dofs = &ned.cell_dofs[ci];
    // Map each boundary dof to its local basis index `s` in cell `ci`
    let local_slots: Vec<usize> = bnd_dofs
        .iter()
        .map(|&gd| {
            cell_dofs
                .iter()
                .position(|&cd| cd == gd)
                .expect("every boundary dof of a face belongs to its adjacent tetrahedron")
        })
        .collect();

    let gid = ned.keys(mesh, ci);
    let mut gs = [[0.0; 3]; 4];
    for b in 0..4 {
        for k in 0..3 {
            gs[b][k] = geo.g[ci * 12 + b * 3 + k];
        }
    }
    let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
    let quad2 = quad_rule(2, 4);
    let nq = quad2.w.len();
    let nb_face = bnd_dofs.len();
    let mut weights = Vec::with_capacity(nq);
    let mut psi_q = Vec::with_capacity(nq);
    let (mut w, mut cu) = (Vec::new(), Vec::new());

    for q in 0..nq {
        let l2 = &quad2.lam[q];
        let mut lam4 = [0.0; 4];
        lam4[lf[0]] = l2[0];
        lam4[lf[1]] = l2[1];
        lam4[lf[2]] = l2[2];
        let (area_q, n_q) = match &curved {
            Some(ct) => {
                let gq = safe_curved_jacobian(ct, &lam4, &gs);
                basis3(ned.order, &lam4, &gq, &gid, &mut w, &mut cu);
                let d = ct.dx_dlam(&lam4);
                let nr = cross3(sub3(d[lf[1]], d[lf[0]]), sub3(d[lf[2]], d[lf[0]]));
                let nl = norm3(nr);
                if nl > 1e-30 {
                    let mut nq = [nr[0] / nl, nr[1] / nl, nr[2] / nl];
                    if dot3(nq, n_hat) < 0.0 {
                        nq = [-nq[0], -nq[1], -nq[2]];
                    }
                    (0.5 * nl, nq)
                } else {
                    (area, n_hat)
                }
            }
            None => {
                basis3(ned.order, &lam4, &gs, &gid, &mut w, &mut cu);
                (area, n_hat)
            }
        };
        let mut row = Vec::with_capacity(nb_face);
        for &s in &local_slots {
            let ws = w[s];
            let wn = dot3(ws, n_q);
            row.push([ws[0] - wn * n_q[0], ws[1] - wn * n_q[1], ws[2] - wn * n_q[2]]);
        }
        weights.push(quad2.w[q] * area_q);
        psi_q.push(row);
    }

    Ok((area, n_hat, bnd_dofs, weights, psi_q))
}

/// Deduplicates and assembles all 3D `pec`, `surface_impedance`, and `lumped_jj` boundary conditions
/// following the Mandatory Single-Facet Assembly Rule (direct IBC assembly, for `driven` and `bloch`).
pub fn assemble_ibc_and_jj_3d(
    mesh: &Discretisation,
    ned: &Ned3,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f0: f64,
) -> Result<IbcJjAssembly3D, String> {
    assemble_ibc_and_jj_3d_with_mode(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f0, false)
}

/// Deduplicates and assembles all 3D `pec`, `surface_impedance`, and `lumped_jj` boundary conditions
/// for `cavity` studies, supporting `ibc_mode` (`"direct" | "perturbative" | "auto"`, FR-30).
pub fn assemble_ibc_and_jj_3d_cavity(
    mesh: &Discretisation,
    ned: &Ned3,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f0: f64,
) -> Result<IbcJjAssembly3D, String> {
    assemble_ibc_and_jj_3d_with_mode(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f0, true)
}

fn assemble_ibc_and_jj_3d_with_mode(
    mesh: &Discretisation,
    ned: &Ned3,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f0: f64,
    allow_perturbative: bool,
) -> Result<IbcJjAssembly3D, String> {
    let n = ned.n_dofs();
    let geo = geometry(mesh);
    let face_to_cell = build_face_to_cell_3d(mesh);

    let mut seen_pec = BTreeSet::new();
    let mut pec_tris = Vec::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if seen_pec.insert(key) {
                pec_tris.push(key);
            }
        }
    }

    let mut seen_ibc = BTreeSet::new();
    let mut ibc_tris = Vec::new();
    let mut pert_ibc_tris = Vec::new();
    let mut direct_ibc_labels = Vec::new();
    let mut warnings = Vec::new();
    let mut k_ibc_coo = Coo::new(n, n);
    let mut c_ibc_coo = Coo::new(n, n);
    let mut kc_ibc_coo = CCoo::new(n, n);
    let mut face_sheet_map: HashMap<[u32; 3], SheetCellKind> = HashMap::new();

    for &tri in &pec_tris {
        if let Some(&(ci, lf)) = face_to_cell.get(&tri) {
            let v = mesh.cell(ci);
            let (p0, p1, p2) = (pt3(mesh, v[lf[0]]), pt3(mesh, v[lf[1]]), pt3(mesh, v[lf[2]]));
            let opp = 6 - lf[0] - lf[1] - lf[2];
            let p_opp = pt3(mesh, v[opp]);
            let mut n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
            if dot3(n_raw, sub3(p_opp, p0)) < 0.0 {
                n_raw = [-n_raw[0], -n_raw[1], -n_raw[2]];
            }
            let nl = norm3(n_raw).max(1e-300);
            face_sheet_map.insert(tri, SheetCellKind::Pec { n_hat: [n_raw[0] / nl, n_raw[1] / nl, n_raw[2] / nl] });
        }
    }

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f0 };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha_k = MU0 * inv_l_eff;
        let alpha_c = MU0 * C0 * g_eff;

        // Collect unique triangles belonging to this surface_impedance BC first so we can measure
        // sheet mesh resolution for `ibc_mode = "auto"` and direct-mode diagnostics (FR-30).
        let mut bc_tris = Vec::new();
        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            // Mandatory Single-Facet Rule: integrate each unique triangle at most once
            if seen_pec.contains(&key) || !seen_ibc.insert(key) {
                continue;
            }
            bc_tris.push(key);
        }
        if bc_tris.is_empty() {
            continue;
        }

        let mut edge_len_sum = 0.0_f64;
        let mut area_tot = 0.0_f64;
        let mut edge_counts: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        for &tri in &bc_tris {
            let (p0, p1, p2) = (pt3(mesh, tri[0]), pt3(mesh, tri[1]), pt3(mesh, tri[2]));
            let l01 = norm3(sub3(p1, p0));
            let l12 = norm3(sub3(p2, p1));
            let l20 = norm3(sub3(p0, p2));
            edge_len_sum += (l01 + l12 + l20) / 3.0;
            area_tot += 0.5 * norm3(cross3(sub3(p1, p0), sub3(p2, p0)));
            for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
                *edge_counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let h_sheet_mean = edge_len_sum / (bc_tris.len() as f64);
        let perim_tot: f64 = edge_counts
            .iter()
            .filter(|(_, &cnt)| cnt == 1)
            .map(|(&(a, b), _)| norm3(sub3(pt3(mesh, b), pt3(mesh, a))))
            .sum();
        let est_tris_across_width = if perim_tot > 0.0 && h_sheet_mean > 0.0 {
            (2.0 * area_tot / perim_tot) / h_sheet_mean
        } else {
            f64::INFINITY
        };
        let lambda_k = bc.l_k_sq / MU0;
        let z_s_mag = bc.r_s.hypot(omega_ref * bc.l_k_sq);
        let use_perturbative = allow_perturbative
            && match bc.ibc_mode.as_str() {
                "perturbative" => true,
                "auto" => ned.order == 1 && h_sheet_mean > 10.0 * lambda_k && z_s_mag < 0.05 * ETA0,
                _ => false,
            };

        if use_perturbative {
            for &key in &bc_tris {
                pert_ibc_tris.push(PerturbativeIbcTri { tri: key, l_k_sq: bc.l_k_sq, r_s: bc.r_s, f_ref });
                if let Some(&(ci, lf)) = face_to_cell.get(&key) {
                    let v = mesh.cell(ci);
                    let (p0, p1, p2) = (pt3(mesh, v[lf[0]]), pt3(mesh, v[lf[1]]), pt3(mesh, v[lf[2]]));
                    let opp = 6 - lf[0] - lf[1] - lf[2];
                    let p_opp = pt3(mesh, v[opp]);
                    let mut n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
                    if dot3(n_raw, sub3(p_opp, p0)) < 0.0 {
                        n_raw = [-n_raw[0], -n_raw[1], -n_raw[2]];
                    }
                    let nl = norm3(n_raw).max(1e-300);
                    face_sheet_map
                        .insert(key, SheetCellKind::Pec { n_hat: [n_raw[0] / nl, n_raw[1] / nl, n_raw[2] / nl] });
                }
            }
            continue;
        }

        direct_ibc_labels.extend(bc.labels.iter().copied());
        if allow_perturbative
            && ned.order == 1
            && ((bc.l_k_sq > 0.0 && h_sheet_mean > 20.0 * lambda_k) || est_tris_across_width < 3.0)
        {
            warnings.push(format!(
                "warning: surface_impedance in 'direct' mode at order 1 has coarse strip resolution (h_tri = {h_sheet_mean:.3e} m, L_k_sq/mu_0 = {lambda_k:.3e} m, est_tris_across_width = {est_tris_across_width:.2}); consider ibc_mode = 'perturbative' or 'auto' or order = 2"
            ));
        }

        for key in bc_tris {
            ibc_tris.push(key);
            let (_, n_hat, bnd_dofs, weights, psi_q) =
                eval_triangle_tangential_basis(mesh, ned, &geo, &face_to_cell, key)?;
            face_sheet_map.insert(key, SheetCellKind::Ibc { inv_l_eff, g_eff, n_hat });

            let nb = bnd_dofs.len();
            let mut loc_g = vec![0.0; nb * nb];
            for (q, &wq) in weights.iter().enumerate() {
                let pq = &psi_q[q];
                for a in 0..nb {
                    for b in 0..nb {
                        loc_g[a * nb + b] += wq * dot3(pq[a], pq[b]);
                    }
                }
            }
            if alpha_k > 0.0 {
                let ke: Vec<f64> = loc_g.iter().map(|&g| alpha_k * g).collect();
                k_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &ke);
            }
            if alpha_c > 0.0 {
                let ce: Vec<f64> = loc_g.iter().map(|&g| alpha_c * g).collect();
                c_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &ce);
            }
            if alpha_k > 0.0 || alpha_c > 0.0 {
                let kce: Vec<C> = loc_g.iter().map(|&g| c(alpha_k * g, k0_ref * alpha_c * g)).collect();
                kc_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &kce);
            }
        }
    }

    let mut seen_jj_all = BTreeSet::new();
    let mut jj_tris = Vec::new();
    let mut k_jj_coo = Coo::new(n, n);
    let mut m_jj_coo = Coo::new(n, n);
    let mut jj_ports = Vec::with_capacity(jj_bcs.len());

    for (port_index, bc) in jj_bcs.iter().enumerate() {
        let mut port_tris = Vec::new();
        let mut seen_port = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if seen_pec.contains(&key) || !seen_port.insert(key) {
                continue;
            }
            port_tris.push(key);
            if seen_jj_all.insert(key) {
                jj_tris.push(key);
            }
        }
        if port_tris.is_empty() {
            return Err(format!("lumped_jj '{}' has no active boundary triangles", bc.name));
        }

        // Compute consistent sheet normal n_hat across the port triangles
        let mut n_acc = [0.0; 3];
        let mut n_ref: Option<[f64; 3]> = None;
        let mut port_verts = BTreeSet::new();
        for &tri in &port_tris {
            for v in tri {
                port_verts.insert(v);
            }
            let (p0, p1, p2) = (pt3(mesh, tri[0]), pt3(mesh, tri[1]), pt3(mesh, tri[2]));
            let mut nr = cross3(sub3(p1, p0), sub3(p2, p0));
            if let Some(r) = n_ref {
                if dot3(nr, r) < 0.0 {
                    nr = [-nr[0], -nr[1], -nr[2]];
                }
            } else if norm3(nr) > 0.0 {
                n_ref = Some(nr);
            }
            n_acc[0] += nr[0];
            n_acc[1] += nr[1];
            n_acc[2] += nr[2];
        }
        let nl = norm3(n_acc);
        if nl <= 0.0 {
            return Err(format!("lumped_jj '{}' triangles have zero normal", bc.name));
        }
        let n_hat = [n_acc[0] / nl, n_acc[1] / nl, n_acc[2] / nl];

        // Orthonormal in-plane frame (d_hat, w_hat = n_hat x d_hat) preserving the sign of bc.direction
        let dn = dot3(bc.direction, n_hat);
        let d_proj =
            [bc.direction[0] - dn * n_hat[0], bc.direction[1] - dn * n_hat[1], bc.direction[2] - dn * n_hat[2]];
        let dl = norm3(d_proj);
        if dl <= 1e-12 {
            return Err(format!(
                "lumped_jj '{}': direction {:?} is normal to the sheet plane (normal {:?})",
                bc.name, bc.direction, n_hat
            ));
        }
        let d_hat = [d_proj[0] / dl, d_proj[1] / dl, d_proj[2] / dl];
        let w_raw = cross3(n_hat, d_hat);
        let wl = norm3(w_raw).max(1e-300);
        let w_hat = [w_raw[0] / wl, w_raw[1] / wl, w_raw[2] / wl];

        let (mut u_min, mut u_max) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut v_min, mut v_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &port_verts {
            let p = pt3(mesh, v);
            let u = dot3(p, d_hat);
            let w = dot3(p, w_hat);
            u_min = u_min.min(u);
            u_max = u_max.max(u);
            v_min = v_min.min(w);
            v_max = v_max.max(w);
        }
        let length_m = (u_max - u_min).max(1e-15);
        let width_m = (v_max - v_min).max(1e-15);

        let inv_l_tot = jj_inv_l_tot(bc.l_j, bc.e_l_hz);
        let alpha_k = MU0 * inv_l_tot;
        let alpha_m = bc.c_j / EPS0;

        let mut b_full = vec![0.0; n];
        let mut dist_k_coo = Coo::new(n, n);
        let mut dist_m_coo = Coo::new(n, n);

        for &tri in &port_tris {
            face_sheet_map.insert(tri, SheetCellKind::Jj { port_index, d_hat, width_m });
            let (_, _, bnd_dofs, weights, psi_q) = eval_triangle_tangential_basis(mesh, ned, &geo, &face_to_cell, tri)?;
            let nb = bnd_dofs.len();
            let mut loc_dd = vec![0.0; nb * nb];
            for (q, &wq) in weights.iter().enumerate() {
                let pq = &psi_q[q];
                for a in 0..nb {
                    let pa_d = dot3(pq[a], d_hat);
                    b_full[bnd_dofs[a]] += (wq / width_m) * pa_d;
                    if bc.formulation == "distributed_sheet" {
                        for b in 0..nb {
                            let pb_d = dot3(pq[b], d_hat);
                            loc_dd[a * nb + b] += wq * (length_m / width_m) * pa_d * pb_d;
                        }
                    }
                }
            }
            if bc.formulation == "distributed_sheet" {
                if alpha_k > 0.0 {
                    let ke: Vec<f64> = loc_dd.iter().map(|&g| alpha_k * g).collect();
                    dist_k_coo.push_block(&bnd_dofs, &bnd_dofs, &ke);
                    k_jj_coo.push_block(&bnd_dofs, &bnd_dofs, &ke);
                }
                if alpha_m > 0.0 {
                    let me: Vec<f64> = loc_dd.iter().map(|&g| alpha_m * g).collect();
                    dist_m_coo.push_block(&bnd_dofs, &bnd_dofs, &me);
                    m_jj_coo.push_block(&bnd_dofs, &bnd_dofs, &me);
                }
            }
        }

        let (k_port_full, m_port_full) = if bc.formulation == "rank1" {
            let nz: Vec<(usize, f64)> = b_full.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 0.0).collect();
            let mut pk = Coo::new(n, n);
            let mut pm = Coo::new(n, n);
            for &(r, br) in &nz {
                for &(cc, bc_val) in &nz {
                    let bb = br * bc_val;
                    if alpha_k > 0.0 {
                        pk.push(r, cc, alpha_k * bb);
                        k_jj_coo.push(r, cc, alpha_k * bb);
                    }
                    if alpha_m > 0.0 {
                        pm.push(r, cc, alpha_m * bb);
                        m_jj_coo.push(r, cc, alpha_m * bb);
                    }
                }
            }
            (pk.to_csr(), pm.to_csr())
        } else {
            (dist_k_coo.to_csr(), dist_m_coo.to_csr())
        };

        jj_ports.push(AssembledJjPort {
            bc: bc.clone(),
            d_hat,
            w_hat,
            n_hat,
            length_m,
            width_m,
            u_min,
            u_max,
            inv_l_tot,
            triangles: port_tris,
            b_full,
            b_red: Vec::new(),
            k_jj_red: k_port_full,
            m_jj_red: m_port_full,
        });
    }

    let mut cell_sheet = vec![None; mesh.n_cells()];
    if !face_sheet_map.is_empty() {
        for ci in 0..mesh.n_cells() {
            let v = mesh.cell(ci);
            for f in FACES {
                let key = sorted3(v[f[0]], v[f[1]], v[f[2]]);
                if let Some(&kind) = face_sheet_map.get(&key) {
                    cell_sheet[ci] = Some(kind);
                    break;
                }
            }
        }
    }

    Ok(IbcJjAssembly3D {
        pec_tris,
        ibc_tris,
        pert_ibc_tris,
        direct_ibc_labels,
        jj_tris,
        k_ibc_full: k_ibc_coo.to_csr(),
        c_ibc_full: c_ibc_coo.to_csr(),
        kc_ibc_full: kc_ibc_coo.to_csr(),
        k_jj_full: k_jj_coo.to_csr(),
        m_jj_full: m_jj_coo.to_csr(),
        jj_ports,
        cell_sheet,
        warnings,
    })
}

/// Real congruence transform `P^H A P` for a real matrix `A` under a periodic/Dirichlet reduction `red`.
pub fn reduce_real_csr(red: &femcore::periodic::Reduction, pmat: &CCsr, a: &RCsr) -> RCsr {
    if a.nnz() == 0 {
        return RCsr::zeros(red.n_reduced(), red.n_reduced());
    }
    femcore::periodic::congruence(pmat, &a.to_complex()).real_part()
}

/// Restricts a full real load vector `b_full` (`n`) to the reduced space `b_red = P^T b_full` (`nred`).
pub fn reduce_real_vec(red: &femcore::periodic::Reduction, b_full: &[f64]) -> Vec<f64> {
    let mut b_red = vec![0.0; red.n_reduced()];
    for (i, &val) in b_full.iter().enumerate() {
        let r = red.index_of[i];
        if r != usize::MAX && val != 0.0 {
            b_red[r] += red.ratio[i].re * val;
        }
    }
    b_red
}

/// Computes iterative symmetric Ruiz diagonal scaling factors `d` on `|K| + sigma |M|`,
/// returning `(d, K' = D K D, M' = D M D, G' = D^-1 G, G'^H)`.
///
/// Note on nullspace consistency: since `K G = 0` and `K' = D K D`, the nullspace of `K'` is
/// `G' = D^-1 G` (`K' G' = D K D D^-1 G = D K G = 0`), and `G'^H M' G' = G^H M G`.
pub fn ruiz_equilibrate(kc: &CCsr, mc: &CCsr, g_red: &CCsr, sigma: f64) -> (Vec<f64>, CCsr, CCsr, CCsr, CCsr) {
    let n = kc.nrows();
    let mut d = vec![1.0_f64; n];
    let sig = sigma.abs().max(1.0);
    for _ in 0..8 {
        let mut rmax = vec![0.0_f64; n];
        for i in 0..n {
            let di = d[i];
            let mut m = 0.0_f64;
            let (ck, vk) = kc.row(i);
            for (&j, &v) in ck.iter().zip(vk) {
                m = m.max(di * v.norm() * d[j as usize]);
            }
            let (cm, vm) = mc.row(i);
            for (&j, &v) in cm.iter().zip(vm) {
                m = m.max(di * sig * v.norm() * d[j as usize]);
            }
            rmax[i] = m;
        }
        for i in 0..n {
            if rmax[i] > 1e-300 {
                d[i] = (d[i] / rmax[i].sqrt()).clamp(1e-14, 1e14);
            }
        }
    }

    let scale_sym = |mat: &CCsr| -> CCsr {
        let mut out = mat.clone();
        for i in 0..n {
            let di = d[i];
            let range = out.row_range(i);
            let cols: Vec<u32> = out.indices()[range.clone()].to_vec();
            let vals = &mut out.values_mut()[range];
            for (v, &col) in vals.iter_mut().zip(&cols) {
                *v *= di * d[col as usize];
            }
        }
        out
    };

    let kc_scaled = scale_sym(kc);
    let mc_scaled = scale_sym(mc);

    let mut g_scaled = g_red.clone();
    for i in 0..g_scaled.nrows() {
        let inv_di = 1.0 / d[i];
        let range = g_scaled.row_range(i);
        for v in &mut g_scaled.values_mut()[range] {
            *v *= inv_di;
        }
    }
    let gh_scaled = femcore::periodic::ccsr_adjoint(&g_scaled);
    (d, kc_scaled, mc_scaled, g_scaled, gh_scaled)
}

/// Unscales a Ruiz-equilibrated eigenvector `x = D y` and computes its exact Rayleigh-quotient
/// complex frequency `f = (c_0 / (2 pi)) sqrt((x^H K x) / (x^H M x))` on the unscaled matrices `(kc, mc)`.
pub fn ruiz_unscale_and_rayleigh(d: &[f64], y: &[C], kc: &CCsr, mc: &CCsr) -> (Vec<C>, C) {
    let x: Vec<C> = y.iter().zip(d).map(|(&yi, &di)| yi * di).collect();
    let kx = kc.matvec(&x);
    let mx = mc.matvec(&x);
    let num: C = x.iter().zip(&kx).fold(c(0.0, 0.0), |acc, (&xi, &kxi)| acc + xi.conj() * kxi);
    let den: C = x.iter().zip(&mx).fold(c(0.0, 0.0), |acc, (&xi, &mxi)| acc + xi.conj() * mxi);
    let lam = if den.norm() > 1e-300 { num / den } else { c(0.0, 0.0) };
    let mut f = lam.sqrt() * (C0 / (2.0 * PI));
    if f.re < 0.0 {
        f = -f;
    }
    (x, f)
}

/// Unscales a Ruiz-equilibrated eigenvector `x = D y` and evaluates the mode-exact damped Rayleigh quotient
/// `k_R^2 = (x^H K_no_damp x) / (x^H M x)` and `gamma_C = (x^H C_damp x) / (x^H M x)`, returning
/// `f = (c_0 / (2 pi)) sqrt(k_R^2 + i k_R gamma_C)` so `Q_m = k_R / gamma_C` and `kappa_m / (2 pi) = (c_0 / (2 pi)) gamma_C`
/// are independent of `f_ref`/`f_search` and preserve `Re(k_0^2) = k_R^2` energy equipartition.
pub fn ruiz_unscale_and_rayleigh_damped(
    d: &[f64],
    y: &[C],
    kc_no_damp: &CCsr,
    c_damp: &RCsr,
    mc: &CCsr,
) -> (Vec<C>, C) {
    let x: Vec<C> = y.iter().zip(d).map(|(&yi, &di)| yi * di).collect();
    let kx = kc_no_damp.matvec(&x);
    let mx = mc.matvec(&x);
    let num_k: C = x.iter().zip(&kx).fold(c(0.0, 0.0), |acc, (&xi, &kxi)| acc + xi.conj() * kxi);
    let den: C = x.iter().zip(&mx).fold(c(0.0, 0.0), |acc, (&xi, &mxi)| acc + xi.conj() * mxi);
    if den.norm() <= 1e-300 {
        return (x, c(0.0, 0.0));
    }
    let q_k = num_k / den;
    let lam = if c_damp.nnz() > 0 {
        let mut num_c = c(0.0, 0.0);
        for i in 0..c_damp.nrows() {
            let xi_conj = x[i].conj();
            let (cols, vals) = c_damp.row(i);
            let mut rs = c(0.0, 0.0);
            for (&col, &v) in cols.iter().zip(vals) {
                rs += x[col as usize] * v;
            }
            num_c += xi_conj * rs;
        }
        let gamma_c = num_c / den;
        let k_r = q_k.re.max(0.0).sqrt();
        // Overdamped non-oscillatory RC relaxation pole (gamma_c > 2 k_r, e.g. floating feedline on lumped_port):
        // its exact QEP root k_0^2 - i gamma_c k_0 - k_r^2 = 0 is purely imaginary (Re(f) = 0).
        if gamma_c.re > 2.0 * k_r && q_k.im.abs() <= 1e-6 * q_k.re.abs().max(1.0) {
            let disc = (gamma_c.re * gamma_c.re - 4.0 * q_k.re.max(0.0)).max(0.0).sqrt();
            let k_im = 0.5 * (gamma_c.re + disc);
            return (x, c(0.0, k_im * (C0 / (2.0 * PI))));
        }
        q_k + c(0.0, 1.0) * (gamma_c * k_r)
    } else {
        q_k
    };
    let mut f = lam.sqrt() * (C0 / (2.0 * PI));
    if f.re < 0.0 {
        f = -f;
    }
    (x, f)
}

// ------------------------------------------------------------------ 2D PEC and surface_impedance assembly

/// Builds a map from each sorted 1D edge `(v0, v1)` in a 2D triangular mesh to one adjacent triangle `ci`
/// and its 2 local vertex indices `[la, lb]` in `0..3`.
fn build_edge_to_cell_2d(mesh: &Discretisation) -> HashMap<(u32, u32), (usize, [usize; 2])> {
    let mut map = HashMap::with_capacity(mesh.n_edges());
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for [la, lb] in [[0, 1], [1, 2], [0, 2]] {
            let key = (v[la].min(v[lb]), v[la].max(v[lb]));
            map.entry(key).or_insert((ci, [la, lb]));
        }
    }
    map
}

/// Collects PEC boundary edge indices `be` and vertex indices `bn` for a 2D mesh, honoring explicit `pec`
/// facets and excluding `surface_impedance` and `pmc` facets from exterior default PEC walls.
pub fn pec_edges_and_nodes_2d(
    mesh: &Discretisation,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
) -> (Vec<usize>, Vec<usize>) {
    let exterior = mesh.exterior_facets();
    let ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let is_pec_facet = |f: usize| -> bool {
        let l = mesh.facet_labels[f];
        pec_labels.contains(&l) || (exterior[f] && !ibc_labels.contains(&l) && !pmc_labels.contains(&l))
    };
    let mut be: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
        .collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
        .collect();
    bn.sort_unstable();
    bn.dedup();
    (be, bn)
}

/// Assembles 2D `surface_impedance` contributions into `(a_coo, b_coo)` for `modes` / `modes_complex`.
/// In the 2D mixed waveguide formulation:
/// - `A_{ee} = S - k0^2 T_e` gets `+(mu_0 / L_eff + i k0 mu_0 c0 G_eff) int_gamma (w_a . tau)(w_b . tau) dl`,
/// - `B_{phi,phi} = K_p - k0^2 T_p` gets `+(mu_0 / L_eff + i k0 mu_0 c0 G_eff) int_gamma phi_a phi_b dl`.
pub fn add_ibc_2d_waveguide(
    mesh: &Discretisation,
    order: usize,
    freq: f64,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    a_coo: &mut CCoo,
    b_coo: &mut CCoo,
) -> Result<(), String> {
    if ibc_bcs.is_empty() {
        return Ok(());
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let rule = facet_rule(1);
    let edge_to_cell = build_edge_to_cell_2d(mesh);

    let mut seen_pec = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            seen_pec.insert((fv[0].min(fv[1]), fv[0].max(fv[1])));
        }
    }

    let mut seen_ibc = BTreeSet::new();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { freq };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha = c(MU0 * inv_l_eff, k0_ref * MU0 * C0 * g_eff);
        if alpha.norm() == 0.0 {
            continue;
        }

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
            if seen_pec.contains(&key) || !seen_ibc.insert(key) {
                continue;
            }
            let &(ci, [la, lb]) = edge_to_cell
                .get(&key)
                .ok_or_else(|| format!("1D facet {:?} is not an edge of any 2D triangle", key))?;
            let (p0, p1) = (mesh.point(fv[0] as usize), mesh.point(fv[1] as usize));
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 0.0 {
                continue;
            }
            let tau = [dx / len, dy / len];
            let g = tri_g(&geo, ci);
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();

            let mut ke_e = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut ke_p = vec![c(0.0, 0.0); nl * nl];
            for (q, l1) in rule.lam.iter().enumerate() {
                let mut l2 = [0.0; 3];
                l2[la] = l1[0];
                l2[lb] = l1[1];
                edge_basis(&l2, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv, _) = shape(&l2, order);
                let wq = rule.w[q] * len;
                for s in 0..nloc_e {
                    let ws_t = w[s][0] * tau[0] + w[s][1] * tau[1];
                    for t in 0..nloc_e {
                        let wt_t = w[t][0] * tau[0] + w[t][1] * tau[1];
                        ke_e[s * nloc_e + t] += alpha * (wq * ws_t * wt_t);
                    }
                }
                for a in 0..nl {
                    for b in 0..nl {
                        ke_p[a * nl + b] += alpha * (wq * nv[a] * nv[b]);
                    }
                }
            }
            a_coo.push_block(rows_e, rows_e, &ke_e);
            b_coo.push_block(&cols_p, &cols_p, &ke_p);
        }
    }
    Ok(())
}

/// Assembles 2D `surface_impedance` contributions into `a_coo` for `axisymmetric_modes` (`(r, z)` half-plane).
pub fn add_ibc_2d_axisym(
    mesh: &Discretisation,
    order: usize,
    f_search: f64,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    a_coo: &mut CCoo,
) -> Result<(), String> {
    if ibc_bcs.is_empty() {
        return Ok(());
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let rule = facet_rule(1);
    let edge_to_cell = build_edge_to_cell_2d(mesh);

    let mut seen_pec = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            seen_pec.insert((fv[0].min(fv[1]), fv[0].max(fv[1])));
        }
    }

    let mut seen_ibc = BTreeSet::new();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f_search };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha = c(MU0 * inv_l_eff, k0_ref * MU0 * C0 * g_eff);
        if alpha.norm() == 0.0 {
            continue;
        }

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
            if seen_pec.contains(&key) || !seen_ibc.insert(key) {
                continue;
            }
            let &(ci, [la, lb]) = edge_to_cell
                .get(&key)
                .ok_or_else(|| format!("1D facet {:?} is not an edge of any 2D triangle", key))?;
            let (p0, p1) = (mesh.point(fv[0] as usize), mesh.point(fv[1] as usize));
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 0.0 {
                continue;
            }
            let tau = [dx / len, dy / len];
            let g = tri_g(&geo, ci);
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
            let v_cell = mesh.cell(ci);

            let mut ke_e = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut ke_p = vec![c(0.0, 0.0); nl * nl];
            for (q, l1) in rule.lam.iter().enumerate() {
                let mut l2 = [0.0; 3];
                l2[la] = l1[0];
                l2[lb] = l1[1];
                let r: f64 = (0..3).map(|k| l2[k] * mesh.point(v_cell[k] as usize)[0]).sum();
                edge_basis(&l2, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv, _) = shape(&l2, order);
                let wr = rule.w[q] * len * r;
                for s in 0..nloc_e {
                    let ws_t = w[s][0] * tau[0] + w[s][1] * tau[1];
                    for t in 0..nloc_e {
                        let wt_t = w[t][0] * tau[0] + w[t][1] * tau[1];
                        ke_e[s * nloc_e + t] += alpha * (wr * ws_t * wt_t);
                    }
                }
                for a in 0..nl {
                    for b in 0..nl {
                        ke_p[a * nl + b] += alpha * (wr * nv[a] * nv[b]);
                    }
                }
            }
            a_coo.push_block(rows_e, rows_e, &ke_e);
            a_coo.push_block(&cols_p, &cols_p, &ke_p);
        }
    }
    Ok(())
}
