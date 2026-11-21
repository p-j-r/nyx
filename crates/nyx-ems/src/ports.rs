// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 3D lumped sheet ports (`lumped_port`) and 2D boundary-mode waveguide ports (`waveport`) for the
//! driven frequency-domain solver (`Study::Driven`).
//!
//! * `lumped_port`: rank-1 sheet voltage/current terminal operator with characteristic impedance `Z_0`
//!   and current-flow direction `direction`.
//! * `waveport`: extracts a 2D `SimplexMesh` on the port boundary plane `(u_hat, v_hat, n_out)`,
//!   inherits adjacent 3D cell `eps_r` and `mu_r` and 1D PEC/PMC/IBC conductor edges, solves the 2D
//!   waveguide eigenvalue problem (`nyx_ems::modes`) to obtain the `mode`-th guided/TEM mode profile
//!   `e_{t,p}`, dispersive propagation constant `beta_p(omega)` (with smooth phase alignment across
//!   anchors), and modal wave impedance `Z_{0,p}(omega) = omega mu_0 / beta_p(omega)`, and assembles
//!   the sparse 3D waveport boundary matrix `W_p` and modal vector `b_p`.

use crate::cavity::{basis3, Ned3};
use crate::ibc::{build_face_to_cell_3d, sorted3, MU0};
use crate::modes::{assemble_modes_with_bcs, solve_mode_pairs, EpsCell, C0};
use crate::problem::{LumpedPortBc, SurfaceImpedanceBc, WavePortBc};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, quad_rule, CurvedTet};
use femcore::linalg::csr::{Coo, RCsr};
use femcore::mesh::SimplexMesh;
use nyx_physics::Discretisation;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;

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

/// Assembled operators and geometric metadata for one `lumped_port` terminal port.
#[derive(Clone, Debug)]
pub struct AssembledLumpedPort {
    pub bc: LumpedPortBc,
    pub d_hat: [f64; 3],
    pub w_hat: [f64; 3],
    pub n_hat: [f64; 3],
    pub length_m: f64,
    pub width_m: f64,
    pub triangles: Vec<[u32; 3]>,
    /// Terminal voltage functional `b_p` on all full Nedelec DOFs (`n`).
    pub b_full: Vec<f64>,
    /// Terminal voltage functional `b_p` on reduced free Nedelec DOFs (`nred`).
    pub b_red: Vec<f64>,
    /// Full-DOF rank-1 port damping contribution `C_port,p = (mu_0 c_0 / Z_{0,p}) b_full b_full^T` (`n x n`).
    pub c_port_full: RCsr,
}

/// Assembled operators, 2D modal dispersion, and geometric metadata for one `waveport`.
#[derive(Clone, Debug)]
pub struct AssembledWavePort {
    pub bc: WavePortBc,
    pub n_out: [f64; 3],
    pub u_hat: [f64; 3],
    pub v_hat: [f64; 3],
    pub triangles: Vec<[u32; 3]>,
    /// Modal excitation/projection functional `(b_p)_a = int_{Gamma_p} mu_r^-1 (psi_a)_|| . e_hat_{t,p} dS` on full DOFs (`n`).
    pub b_full: Vec<f64>,
    /// Modal excitation/projection functional `b_p` on reduced DOFs (`nred`).
    pub b_red: Vec<f64>,
    /// Full-DOF waveport tangential boundary matrix `(W_p)_{ab} = int_{Gamma_p} mu_r^-1 (psi_a)_|| . (psi_b)_|| dS` (`n x n`).
    pub w_full: RCsr,
    /// Reduced-DOF waveport tangential boundary matrix `W_{p,red}` (`nred x nred`).
    pub w_red: RCsr,
    /// Quadratic dispersion coefficients `[c0, c1, c2]` for `beta_p^2(k0^2) = c0 + c1 * k0^2 + c2 * k0^4`.
    pub beta2_poly: [f64; 3],
}

impl AssembledWavePort {
    /// Evaluates `(beta_p(omega), d(beta_p)/d(omega), d^2(beta_p)/d(omega)^2, z0_p(omega))` at angular frequency `omega` (rad/s).
    #[inline]
    pub fn beta_and_derivs(&self, omega: f64) -> (f64, f64, f64, f64) {
        let k0 = omega / C0;
        let s = k0 * k0;
        let ds_dw = 2.0 * omega / (C0 * C0);
        let d2s_dw2 = 2.0 / (C0 * C0);
        let [c0, c1, c2] = self.beta2_poly;
        let beta2 = (c0 + c1 * s + c2 * s * s).max(1e-12);
        let dbeta2_ds = c1 + 2.0 * c2 * s;
        let d2beta2_ds2 = 2.0 * c2;
        let dbeta2_dw = dbeta2_ds * ds_dw;
        let d2beta2_dw2 = d2beta2_ds2 * ds_dw * ds_dw + dbeta2_ds * d2s_dw2;
        let beta = beta2.sqrt();
        let dbeta_dw = 0.5 * dbeta2_dw / beta;
        let d2beta_dw2 = 0.5 * d2beta2_dw2 / beta - (dbeta_dw * dbeta_dw) / beta;
        let z0 = (omega * MU0 / beta).max(1e-12);
        (beta, dbeta_dw, d2beta_dw2, z0)
    }
}

/// Evaluates the 3D tangential Nedelec basis vectors `psi_{s,||} = w_s - (w_s . n) n` at the 6-point Dunavant
/// quadrature points on a single triangle `tri_sorted = [v0, v1, v2]` ($v_0 < v_1 < v_2$), with local tetrahedron
/// vertex indices ordered to match `[tri_sorted[0], tri_sorted[1], tri_sorted[2]]` so that 2D and 3D quadrature points
/// coincide point-for-point.
/// Returns `(adjacent_cell_ci, area, n_out, bnd_dofs, quad_weights, psi_parallel_at_q)`.
pub fn eval_triangle_tangential_basis_ordered(
    mesh: &Discretisation,
    ned: &Ned3,
    geo: &femcore::fem::Geometry,
    face_to_cell: &HashMap<[u32; 3], (usize, [usize; 3])>,
    tri_sorted: [u32; 3],
) -> Result<(usize, f64, [f64; 3], Vec<usize>, Vec<f64>, Vec<Vec<[f64; 3]>>), String> {
    let &(ci, _) = face_to_cell
        .get(&tri_sorted)
        .ok_or_else(|| format!("facet {:?} is not a face of any tetrahedron in the mesh", tri_sorted))?;
    let v = mesh.cell(ci);
    let mut lf = [0usize; 3];
    for (k, &tv) in tri_sorted.iter().enumerate() {
        lf[k] = v.iter().position(|&cv| cv == tv).expect("triangle vertex belongs to adjacent tetrahedron");
    }
    let opp = 6 - lf[0] - lf[1] - lf[2];
    let (p0, p1, p2, p_opp) = (pt3(mesh, v[lf[0]]), pt3(mesh, v[lf[1]]), pt3(mesh, v[lf[2]]), pt3(mesh, v[opp]));
    let mut n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
    let n_len = norm3(n_raw);
    let area = 0.5 * n_len;
    if area <= 0.0 {
        return Err("degenerate zero-area boundary triangle".into());
    }
    if dot3(n_raw, sub3(p0, p_opp)) < 0.0 {
        n_raw = [-n_raw[0], -n_raw[1], -n_raw[2]];
    }
    let n_out = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];

    let bnd_dofs = ned.boundary_dofs(mesh, tri_sorted);
    let cell_dofs = &ned.cell_dofs[ci];
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
        match &curved {
            Some(ct) => basis3(ned.order, &lam4, &ct.jacobian(&lam4).0, &gid, &mut w, &mut cu),
            None => basis3(ned.order, &lam4, &gs, &gid, &mut w, &mut cu),
        }
        let mut row = Vec::with_capacity(nb_face);
        for &s in &local_slots {
            let ws = w[s];
            let wn = dot3(ws, n_out);
            row.push([ws[0] - wn * n_out[0], ws[1] - wn * n_out[1], ws[2] - wn * n_out[2]]);
        }
        weights.push(quad2.w[q] * area);
        psi_q.push(row);
    }

    Ok((ci, area, n_out, bnd_dofs, weights, psi_q))
}

/// Assembles all 3D `lumped_port` boundary conditions following the single-facet rule.
pub fn assemble_lumped_ports_3d(
    mesh: &Discretisation,
    ned: &Ned3,
    pec_labels: &[i32],
    lumped_bcs: &[LumpedPortBc],
) -> Result<Vec<AssembledLumpedPort>, String> {
    if lumped_bcs.is_empty() {
        return Ok(Vec::new());
    }
    let n = ned.n_dofs();
    let geo = geometry(mesh);
    let face_to_cell = build_face_to_cell_3d(mesh);

    let mut seen_pec = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            seen_pec.insert(sorted3(fv[0], fv[1], fv[2]));
        }
    }

    let mut ports = Vec::with_capacity(lumped_bcs.len());
    for bc in lumped_bcs {
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
        }
        if port_tris.is_empty() {
            return Err(format!("lumped_port '{}' has no active boundary triangles", bc.name));
        }

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
            return Err(format!("lumped_port '{}' triangles have zero normal", bc.name));
        }
        let n_hat = [n_acc[0] / nl, n_acc[1] / nl, n_acc[2] / nl];

        let dn = dot3(bc.direction, n_hat);
        let d_proj =
            [bc.direction[0] - dn * n_hat[0], bc.direction[1] - dn * n_hat[1], bc.direction[2] - dn * n_hat[2]];
        let dl = norm3(d_proj);
        if dl <= 1e-12 {
            return Err(format!(
                "lumped_port '{}': direction {:?} is normal to the sheet plane (normal {:?})",
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

        let mut b_full = vec![0.0; n];
        for &tri in &port_tris {
            let (_, _, _, bnd_dofs, weights, psi_q) =
                eval_triangle_tangential_basis_ordered(mesh, ned, &geo, &face_to_cell, tri)?;
            for (q, &wq) in weights.iter().enumerate() {
                let pq = &psi_q[q];
                for (a, &dof) in bnd_dofs.iter().enumerate() {
                    b_full[dof] += (wq / width_m) * dot3(pq[a], d_hat);
                }
            }
        }

        let alpha_c = MU0 * C0 / bc.z0;
        let nz: Vec<(usize, f64)> = b_full.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 0.0).collect();
        let mut c_coo = Coo::new(n, n);
        for &(r, br) in &nz {
            for &(cc, bc_val) in &nz {
                c_coo.push(r, cc, alpha_c * br * bc_val);
            }
        }

        ports.push(AssembledLumpedPort {
            bc: bc.clone(),
            d_hat,
            w_hat,
            n_hat,
            length_m,
            width_m,
            triangles: port_tris,
            b_full,
            b_red: Vec::new(),
            c_port_full: c_coo.to_csr(),
        });
    }
    Ok(ports)
}

/// Builds a canonical in-plane orthonormal frame `(u_hat, v_hat)` for a plane with unit normal `n_out`.
/// The frame depends only on the plane orientation (via `n_canon` with first dominant coordinate positive),
/// so parallel opposite port planes (e.g. `y = 0` and `y = L`) share identical `(u_hat, v_hat)` axes.
fn canonical_plane_axes(n_out: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut sign = 1.0;
    for k in 0..3 {
        if n_out[k].abs() > 1e-6 {
            if n_out[k] < 0.0 {
                sign = -1.0;
            }
            break;
        }
    }
    let n_canon = [sign * n_out[0], sign * n_out[1], sign * n_out[2]];
    let ref_axis = if n_canon[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let rn = dot3(ref_axis, n_canon);
    let u_raw = [ref_axis[0] - rn * n_canon[0], ref_axis[1] - rn * n_canon[1], ref_axis[2] - rn * n_canon[2]];
    let ul = norm3(u_raw).max(1e-300);
    let u_hat = [u_raw[0] / ul, u_raw[1] / ul, u_raw[2] / ul];
    let v_raw = cross3(n_canon, u_hat);
    let vl = norm3(v_raw).max(1e-300);
    let v_hat = [v_raw[0] / vl, v_raw[1] / vl, v_raw[2] / vl];
    (u_hat, v_hat)
}

/// Evaluates the 3D transverse electric field `e_t(x_{c,q}) = E_u u_hat + E_v v_hat` at the 6 Dunavant quadrature
/// points of every triangle `c` in `mesh_2d` for a 2D edge-element coefficient vector `et`, normalizes it so
/// `int_{Gamma_p} mu_r^-1 |e_t|^2 dS = 1`, and returns `(e_t_at_cq, canonical_orientation_score)`.
fn evaluate_and_normalize_2d_mode(
    mesh_2d: &SimplexMesh,
    order: usize,
    et: &[f64],
    u_hat: [f64; 3],
    v_hat: [f64; 3],
    mu_r_2d: &[f64],
    tri_weights: &[Vec<f64>],
) -> (Vec<Vec<[f64; 3]>>, f64) {
    let nc = mesh_2d.n_cells();
    let geo_2d = geometry(mesh_2d);
    let cd_2d = cell_data(mesh_2d, order);
    let quad2 = quad_rule(2, 4);
    let nq = quad2.w.len();
    let nloc_e = if order == 1 { 3 } else { 8 };
    let mut w2 = vec![[0.0; 2]; nloc_e];
    let mut curl2 = vec![0.0; nloc_e];

    let mut fields = Vec::with_capacity(nc);
    let mut norm2_acc = 0.0;
    let mut u_c = 0.0;
    let mut v_c = 0.0;
    let mut area_tot = 0.0;

    for c in 0..nc {
        let v = mesh_2d.cell(c);
        let (p0, p1, p2) = (mesh_2d.point(v[0] as usize), mesh_2d.point(v[1] as usize), mesh_2d.point(v[2] as usize));
        let g2 = tri_g(&geo_2d, c);
        let imu = 1.0 / mu_r_2d[c];
        let mut cell_q = Vec::with_capacity(nq);
        for q in 0..nq {
            let l = &quad2.lam[q];
            let uq = l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0];
            let vq = l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1];
            let wq = tri_weights[c][q];
            u_c += wq * uq;
            v_c += wq * vq;
            area_tot += wq;

            edge_basis(l, &g2, &cd_2d.sign[c], order, &mut w2, &mut curl2);
            let (mut eu, mut ev) = (0.0, 0.0);
            for s in 0..nloc_e {
                let coef = et[cd_2d.gdofs[c][s]];
                eu += coef * w2[s][0];
                ev += coef * w2[s][1];
            }
            let e3 = [eu * u_hat[0] + ev * v_hat[0], eu * u_hat[1] + ev * v_hat[1], eu * u_hat[2] + ev * v_hat[2]];
            norm2_acc += wq * imu * dot3(e3, e3);
            cell_q.push(e3);
        }
        fields.push(cell_q);
    }

    let inv_norm = 1.0 / norm2_acc.max(1e-300).sqrt();
    for cell_q in &mut fields {
        for e3 in cell_q {
            e3[0] *= inv_norm;
            e3[1] *= inv_norm;
            e3[2] *= inv_norm;
        }
    }

    u_c /= area_tot.max(1e-300);
    v_c /= area_tot.max(1e-300);
    let span = area_tot.sqrt().max(1e-12);

    // Deterministic orientation moment to fix the global +/- sign on parallel planes
    let pref_dir = [
        u_hat[0] + std::f64::consts::FRAC_1_PI * v_hat[0],
        u_hat[1] + std::f64::consts::FRAC_1_PI * v_hat[1],
        u_hat[2] + std::f64::consts::FRAC_1_PI * v_hat[2],
    ];
    let mut m0 = 0.0;
    let mut m1 = 0.0;
    let mut peak_val = 0.0;
    let mut peak_mag2 = -1.0;

    for c in 0..nc {
        let v = mesh_2d.cell(c);
        let (p0, p1, p2) = (mesh_2d.point(v[0] as usize), mesh_2d.point(v[1] as usize), mesh_2d.point(v[2] as usize));
        for q in 0..nq {
            let l = &quad2.lam[q];
            let uq = (l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0] - u_c) / span;
            let vq = (l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1] - v_c) / span;
            let wq = tri_weights[c][q];
            let e3 = fields[c][q];
            let proj = dot3(e3, pref_dir);
            m0 += wq * proj;
            m1 += wq * (uq + 0.2718281828 * vq) * proj;
            let mag2 = dot3(e3, e3);
            if mag2 > peak_mag2 + 1e-12 && proj.abs() > 1e-6 * mag2.sqrt() {
                peak_mag2 = mag2;
                peak_val = proj;
            }
        }
    }

    let score = if m0.abs() > 1e-4 * area_tot.sqrt() {
        m0
    } else if m1.abs() > 1e-5 * area_tot.sqrt() {
        m1
    } else {
        peak_val
    };

    if score < 0.0 {
        for cell_q in &mut fields {
            for e3 in cell_q {
                e3[0] = -e3[0];
                e3[1] = -e3[1];
                e3[2] = -e3[2];
            }
        }
    }

    (fields, score.abs())
}

/// Computes the $\mu_r^{-1}$-weighted inner product $\langle \mathbf{e}_1, \mathbf{e}_2 \rangle_{\Gamma_p}$.
fn mode_overlap(e1: &[Vec<[f64; 3]>], e2: &[Vec<[f64; 3]>], mu_r_2d: &[f64], tri_weights: &[Vec<f64>]) -> f64 {
    let mut acc = 0.0;
    for (c, (c1, c2)) in e1.iter().zip(e2).enumerate() {
        let imu = 1.0 / mu_r_2d[c];
        for (q, (&v1, &v2)) in c1.iter().zip(c2).enumerate() {
            acc += tri_weights[c][q] * imu * dot3(v1, v2);
        }
    }
    acc
}

/// Assembles all 3D `waveport` boundary conditions by extracting the 2D boundary mesh on each waveport plane
/// and solving the 2D waveguide mode problem (`nyx_ems::modes`) with smooth phase alignment across anchor frequencies.
#[allow(clippy::too_many_arguments)]
pub fn assemble_waveports_3d(
    mesh: &Discretisation,
    ned: &Ned3,
    pec_labels: &[i32],
    pmc_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    waveport_bcs: &[WavePortBc],
    cell_eps_3d: &[[[f64; 3]; 3]],
    cell_mur_3d: &[f64],
    f_min: f64,
    f_max: f64,
) -> Result<Vec<AssembledWavePort>, String> {
    if waveport_bcs.is_empty() {
        return Ok(Vec::new());
    }
    let n = ned.n_dofs();
    let order = ned.order;
    let geo = geometry(mesh);
    let face_to_cell = build_face_to_cell_3d(mesh);
    let exterior = mesh.exterior_facets();

    let all_wp_labels: Vec<i32> = waveport_bcs.iter().flat_map(|w| w.labels.iter().copied()).collect();
    let all_ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();

    // Collect 3D edges belonging to explicit PEC facets, exterior default PEC walls, PMC walls, and IBC sheets
    let mut pec_edges_3d = BTreeSet::new();
    let mut pmc_edges_3d = BTreeSet::new();
    let mut ibc_edges_3d: Vec<BTreeSet<(u32, u32)>> = vec![BTreeSet::new(); ibc_bcs.len()];
    let mut seen_pec_tris = BTreeSet::new();

    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        let tri = sorted3(fv[0], fv[1], fv[2]);
        let is_explicit_pec = pec_labels.contains(&l);
        if is_explicit_pec {
            seen_pec_tris.insert(tri);
        }
        let is_ext_pec =
            exterior[f] && !pmc_labels.contains(&l) && !all_ibc_labels.contains(&l) && !all_wp_labels.contains(&l);
        let edges = [
            (fv[0].min(fv[1]), fv[0].max(fv[1])),
            (fv[1].min(fv[2]), fv[1].max(fv[2])),
            (fv[0].min(fv[2]), fv[0].max(fv[2])),
        ];
        if is_explicit_pec || is_ext_pec {
            for e in edges {
                pec_edges_3d.insert(e);
            }
        } else if pmc_labels.contains(&l) {
            for e in edges {
                pmc_edges_3d.insert(e);
            }
        }
        for (ibc_idx, ibc) in ibc_bcs.iter().enumerate() {
            if ibc.labels.contains(&l) {
                for e in edges {
                    ibc_edges_3d[ibc_idx].insert(e);
                }
            }
        }
    }

    let mut ports = Vec::with_capacity(waveport_bcs.len());
    let f_mid = 0.5 * (f_min + f_max);

    for bc in waveport_bcs {
        let mut port_tris = Vec::new();
        let mut seen_port = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if seen_pec_tris.contains(&key) || !seen_port.insert(key) {
                continue;
            }
            port_tris.push(key);
        }
        if port_tris.is_empty() {
            return Err(format!("waveport '{}' has no active boundary triangles", bc.name));
        }

        // Evaluate 3D tangential basis and outward normals on every triangle of the waveport
        let mut tri_evals = Vec::with_capacity(port_tris.len());
        let mut n_acc = [0.0; 3];
        let mut verts_set = BTreeSet::new();
        for &tri in &port_tris {
            for v in tri {
                verts_set.insert(v);
            }
            let eval = eval_triangle_tangential_basis_ordered(mesh, ned, &geo, &face_to_cell, tri)?;
            let area = eval.1;
            let n_tri = eval.2;
            n_acc[0] += area * n_tri[0];
            n_acc[1] += area * n_tri[1];
            n_acc[2] += area * n_tri[2];
            tri_evals.push(eval);
        }
        let nl = norm3(n_acc);
        if nl <= 0.0 {
            return Err(format!("waveport '{}' triangles have zero normal", bc.name));
        }
        let n_out = [n_acc[0] / nl, n_acc[1] / nl, n_acc[2] / nl];
        for (_, _, n_tri, _, _, _) in &tri_evals {
            if dot3(n_out, *n_tri) < 0.95 {
                return Err(format!("waveport '{}' facets are not coplanar", bc.name));
            }
        }

        let (u_hat, v_hat) = canonical_plane_axes(n_out);

        // Build 2D SimplexMesh of the waveport cross-section
        let verts_3d: Vec<u32> = verts_set.into_iter().collect();
        let mut points_2d = Vec::with_capacity(verts_3d.len() * 2);
        for &v3 in &verts_3d {
            let p = pt3(mesh, v3);
            points_2d.push(dot3(p, u_hat));
            points_2d.push(dot3(p, v_hat));
        }
        let v2d_of = |v3: u32| -> u32 { verts_3d.binary_search(&v3).expect("vertex in verts_3d") as u32 };

        let nc_2d = port_tris.len();
        let mut cells_2d = Vec::with_capacity(nc_2d * 3);
        let mut cell_labels_2d = Vec::with_capacity(nc_2d);
        let mut eps_2d = Vec::with_capacity(nc_2d);
        let mut mur_2d = Vec::with_capacity(nc_2d);
        let mut tri_weights = Vec::with_capacity(nc_2d);
        let mut edge_use_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();

        for (c_idx, (&tri, eval)) in port_tris.iter().zip(&tri_evals).enumerate() {
            let (i0, i1, i2) = (v2d_of(tri[0]), v2d_of(tri[1]), v2d_of(tri[2]));
            cells_2d.extend_from_slice(&[i0, i1, i2]);
            cell_labels_2d.push(c_idx as i32 + 1);
            for (a, b) in [(i0, i1), (i1, i2), (i0, i2)] {
                *edge_use_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }

            let ci = eval.0;
            let t3 = &cell_eps_3d[ci];
            // Project 3D permittivity tensor onto (u_hat, v_hat, n_out)
            let project_eps = |a: [f64; 3], b: [f64; 3]| -> f64 {
                let mut s = 0.0;
                for r in 0..3 {
                    for cc in 0..3 {
                        s += a[r] * t3[r][cc] * b[cc];
                    }
                }
                s
            };
            eps_2d.push(EpsCell {
                xx: project_eps(u_hat, u_hat),
                xy: project_eps(u_hat, v_hat),
                yy: project_eps(v_hat, v_hat),
                zz: project_eps(n_out, n_out),
            });
            mur_2d.push(cell_mur_3d[ci]);
            tri_weights.push(eval.4.clone());
        }

        // Build 1D facets of the 2D waveport mesh:
        // label 1 = PEC, label 2 = PMC, label (3 + ibc_idx) = surface_impedance
        let mut facets_2d = Vec::new();
        let mut facet_labels_2d = Vec::new();
        for (&(i_a, i_b), &count) in &edge_use_count {
            let (v_a, v_b) = (verts_3d[i_a as usize], verts_3d[i_b as usize]);
            let e3 = (v_a.min(v_b), v_a.max(v_b));
            let is_perim = count == 1;
            if pec_edges_3d.contains(&e3) || (is_perim && !pmc_edges_3d.contains(&e3)) {
                facets_2d.extend_from_slice(&[i_a, i_b]);
                facet_labels_2d.push(1);
            } else if is_perim && pmc_edges_3d.contains(&e3) {
                facets_2d.extend_from_slice(&[i_a, i_b]);
                facet_labels_2d.push(2);
            } else {
                for (ibc_idx, set) in ibc_edges_3d.iter().enumerate() {
                    if set.contains(&e3) {
                        facets_2d.extend_from_slice(&[i_a, i_b]);
                        facet_labels_2d.push(3 + ibc_idx as i32);
                        break;
                    }
                }
            }
        }

        let mesh_2d = SimplexMesh::new(2, points_2d, cells_2d, cell_labels_2d, facets_2d, facet_labels_2d);
        let mut ibc_2d = Vec::new();
        for (ibc_idx, ibc) in ibc_bcs.iter().enumerate() {
            let lbl = 3 + ibc_idx as i32;
            if mesh_2d.facet_labels.contains(&lbl) {
                let mut b2 = ibc.clone();
                b2.labels = vec![lbl];
                ibc_2d.push(b2);
            }
        }

        let n_max =
            eps_2d.iter().zip(&mur_2d).map(|(e, &m)| ((e.xx.max(e.yy)) * m).max(1.0).sqrt()).fold(1.0_f64, f64::max);
        let n_search = 1.02 * n_max;
        let k_ask = (bc.mode + 3).max(4);

        let solve_2d_at = |freq: f64| -> Result<Vec<(f64, Vec<Vec<[f64; 3]>>)>, String> {
            let asm2 = assemble_modes_with_bcs(&mesh_2d, freq, &eps_2d, &mur_2d, order, &[1], &ibc_2d, &[2])?;
            let k_actual = k_ask.min(asm2.free.len().saturating_sub(2)).max(1);
            let modes = solve_mode_pairs(&asm2, n_search, false, k_actual)?;
            if modes.neff.is_empty() {
                return Err(format!("waveport '{}': no 2D guided modes converged at {:.4} GHz", bc.name, freq / 1e9));
            }
            let k0 = 2.0 * PI * freq / C0;
            let mut candidates: Vec<(f64, f64, Vec<Vec<[f64; 3]>>)> = Vec::with_capacity(modes.neff.len());
            for (m_i, &(n_re, n_im)) in modes.neff.iter().enumerate() {
                let beta2 = if n_re > 0.0 { (n_re * k0).powi(2) } else { -(n_im * k0).powi(2) };
                let (e_field, score) = evaluate_and_normalize_2d_mode(
                    &mesh_2d,
                    order,
                    &modes.et[m_i],
                    u_hat,
                    v_hat,
                    &mur_2d,
                    &tri_weights,
                );
                candidates.push((beta2, score, e_field));
            }
            // Sort by descending beta^2 (fundamental guided mode first), breaking ties deterministically by orientation score
            candidates.sort_by(|a, b| {
                let scale = a.0.abs().max(b.0.abs()).max(1.0);
                if (a.0 - b.0).abs() > 1e-5 * scale {
                    b.0.partial_cmp(&a.0).unwrap()
                } else {
                    b.1.partial_cmp(&a.1).unwrap()
                }
            });
            Ok(candidates.into_iter().map(|(b2, _, ef)| (b2, ef)).collect())
        };

        let mid_modes = solve_2d_at(f_mid)?;
        if bc.mode > mid_modes.len() {
            return Err(format!(
                "waveport '{}' requested mode {}, but only {} modes were found",
                bc.name,
                bc.mode,
                mid_modes.len()
            ));
        }
        let (beta2_mid, ref_e_field) = mid_modes[bc.mode - 1].clone();
        if beta2_mid <= 0.0 {
            return Err(format!(
                "waveport '{}' mode {} is evanescent (below cutoff) at midband {:.4} GHz",
                bc.name,
                bc.mode,
                f_mid / 1e9
            ));
        }

        // Evaluate dispersion at f_min and f_max with phase/mode tracking against ref_e_field
        let track_beta2_at = |freq: f64| -> Result<f64, String> {
            let cands = solve_2d_at(freq)?;
            let mut best_ov = -1.0_f64;
            let mut best_b2 = beta2_mid;
            for (b2, mut ef) in cands {
                let ov = mode_overlap(&ef, &ref_e_field, &mur_2d, &tri_weights);
                // Enforce smooth phase alignment across anchors: flip sign if inner product < 0
                if ov < 0.0 {
                    for cq in &mut ef {
                        for v in cq {
                            v[0] = -v[0];
                            v[1] = -v[1];
                            v[2] = -v[2];
                        }
                    }
                }
                if ov.abs() > best_ov {
                    best_ov = ov.abs();
                    best_b2 = b2;
                }
            }
            Ok(best_b2)
        };

        let beta2_poly = if (f_max - f_min).abs() <= 1e-9 * f_mid {
            let s_mid = (2.0 * PI * f_mid / C0).powi(2);
            [0.0, beta2_mid / s_mid.max(1e-300), 0.0]
        } else {
            let b2_min = track_beta2_at(f_min)?;
            let b2_max = track_beta2_at(f_max)?;
            let s0 = (2.0 * PI * f_min / C0).powi(2);
            let s1 = (2.0 * PI * f_mid / C0).powi(2);
            let s2 = (2.0 * PI * f_max / C0).powi(2);
            // 3-point Lagrange quadratic fit of beta^2(s) through (s0, b2_min), (s1, beta2_mid), (s2, b2_max)
            let l0_denom = (s0 - s1) * (s0 - s2);
            let l1_denom = (s1 - s0) * (s1 - s2);
            let l2_denom = (s2 - s0) * (s2 - s1);
            let w0 = b2_min / l0_denom;
            let w1 = beta2_mid / l1_denom;
            let w2 = b2_max / l2_denom;
            let c2 = w0 + w1 + w2;
            let c1 = -(w0 * (s1 + s2) + w1 * (s0 + s2) + w2 * (s0 + s1));
            let c0 = w0 * s1 * s2 + w1 * s0 * s2 + w2 * s0 * s1;
            [c0, c1, c2]
        };

        // Assemble 3D waveport boundary matrix W_p and modal vector b_p
        let mut w_coo = Coo::new(n, n);
        let mut b_full = vec![0.0; n];
        for (c_idx, eval) in tri_evals.iter().enumerate() {
            let bnd_dofs = &eval.3;
            let weights = &eval.4;
            let psi_q = &eval.5;
            let imu = 1.0 / mur_2d[c_idx];
            let nb = bnd_dofs.len();
            let mut loc_w = vec![0.0; nb * nb];
            for (q, &wq) in weights.iter().enumerate() {
                let pq = &psi_q[q];
                let eq = ref_e_field[c_idx][q];
                let w_imu = wq * imu;
                for a in 0..nb {
                    b_full[bnd_dofs[a]] += w_imu * dot3(pq[a], eq);
                    for b in 0..nb {
                        loc_w[a * nb + b] += w_imu * dot3(pq[a], pq[b]);
                    }
                }
            }
            w_coo.push_block(bnd_dofs, bnd_dofs, &loc_w);
        }

        ports.push(AssembledWavePort {
            bc: bc.clone(),
            n_out,
            u_hat,
            v_hat,
            triangles: port_tris,
            b_full,
            b_red: Vec::new(),
            w_full: w_coo.to_csr(),
            w_red: RCsr::zeros(0, 0),
            beta2_poly,
        });
    }

    Ok(ports)
}
