// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Modal energy partitioning (`U_elec = 1 J` normalization, bulk magnetic energy, kinetic IBC energy,
//! lumped JJ inductive energy, port voltage/current phasors, equipartition verification), 3D island-tied
//! `H^1` port capacitance extraction (`C_port`), passive Craig-Bampton reduced-order model (`CircuitRom`)
//! assembly, and dielectric `TlsMeta` generation.

use crate::cavity::{gradient_full, EmwAssembly};
use crate::ibc::{jj_inv_l_tot, sorted3, EPS0, MU0};
use crate::modes::C0;
use crate::problem::Problem;
use femcore::fem::solve_dirichlet;
use femcore::linalg::amg_elastic::spmm;
use femcore::linalg::csr::{c, RCsr, C};
use femcore::linalg::dense::Dense;
use femcore::linalg::driver::{dense_gen_eigh, DirectSolver};
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::periodic::Reduction;
use nyx_physics::{CircuitRom, Discretisation, ModeEnergyRow, RomPort, TlsMeta};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

#[inline]
fn dot_real(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

#[inline]
fn quad_form_real_mat_complex_vec(mat: &RCsr, x: &[C]) -> f64 {
    let mut acc = 0.0;
    for i in 0..mat.nrows() {
        let xi = x[i];
        let (cols, vals) = mat.row(i);
        let mut row_sum = c(0.0, 0.0);
        for (&col, &v) in cols.iter().zip(vals) {
            row_sum += x[col as usize] * v;
        }
        acc += (xi.conj() * row_sum).re;
    }
    acc
}

/// Computes the physical scaling factor `s_E > 0` such that `U_elec,m = (1/4) eps_0 (s_E x_m)^H M_tot (s_E x_m) = 1 J`,
/// the exact modal energy row (`ModeEnergyRow`), and the relative equipartition error
/// `|U_elec - (U_mag^vol + U_kin^IBC + sum_j U_ind,j^JJ)| / U_elec`.
pub fn compute_mode_energy_and_scale(f: C, x_red: &[C], rot: C, asm: &EmwAssembly) -> (f64, ModeEnergyRow, f64) {
    let mut m_tot_phys = asm.m_phys_red.clone();
    for port in &asm.jj_ports {
        if port.m_jj_red.nnz() > 0 {
            m_tot_phys = m_tot_phys.add_scaled(&port.m_jj_red, 1.0);
        }
    }
    let q_m = quad_form_real_mat_complex_vec(&m_tot_phys, x_red).max(1e-300);
    let s_e = 2.0 / (EPS0 * q_m).sqrt();
    let scale = rot * s_e;
    let x_phys: Vec<C> = x_red.iter().map(|&z| z * scale).collect();

    let u_elec_j = 0.25 * EPS0 * quad_form_real_mat_complex_vec(&m_tot_phys, &x_phys);
    let omega_m = (2.0 * PI * f.re).max(1e-300);
    let re_k0_sq = ((2.0 * PI / C0).powi(2) * (f.re * f.re - f.im * f.im)).max(1e-300);
    let mag_factor = EPS0 / (4.0 * re_k0_sq);

    let u_mag_vol_j = mag_factor * quad_form_real_mat_complex_vec(&asm.k_phys_red, &x_phys);
    let u_kin_ibc_j = if asm.k_ibc_red.nnz() > 0 {
        mag_factor * quad_form_real_mat_complex_vec(&asm.k_ibc_red, &x_phys)
    } else {
        0.0
    };

    let n_ports = asm.jj_ports.len() + asm.lumped_ports.len();
    let mut u_ind_jj_j = Vec::with_capacity(n_ports);
    let mut v_mj = Vec::with_capacity(n_ports);
    let mut i_mj = Vec::with_capacity(n_ports);

    for port in &asm.jj_ports {
        let mut v_port = c(0.0, 0.0);
        for (&b_a, &x_a) in port.b_red.iter().zip(&x_phys) {
            if b_a != 0.0 {
                v_port += x_a * b_a;
            }
        }
        // i_mj = v_mj / (i omega_m L_tot,j) = -i (inv_l_tot / omega_m) v_mj
        let i_port = c(0.0, -1.0) * (port.inv_l_tot / omega_m) * v_port;
        let u_ind = if port.k_jj_red.nnz() > 0 {
            mag_factor * quad_form_real_mat_complex_vec(&port.k_jj_red, &x_phys)
        } else {
            0.0
        };
        v_mj.push([v_port.re, v_port.im]);
        i_mj.push([i_port.re, i_port.im]);
        u_ind_jj_j.push(u_ind);
    }
    for lp in &asm.lumped_ports {
        let mut v_port = c(0.0, 0.0);
        for (&b_a, &x_a) in lp.b_red.iter().zip(&x_phys) {
            if b_a != 0.0 {
                v_port += x_a * b_a;
            }
        }
        let i_port = v_port / lp.bc.z0.max(1e-30);
        v_mj.push([v_port.re, v_port.im]);
        i_mj.push([i_port.re, i_port.im]);
        u_ind_jj_j.push(0.0);
    }

    let u_mag_tot = u_mag_vol_j + u_kin_ibc_j + u_ind_jj_j.iter().sum::<f64>();
    let equipartition_error = (u_elec_j - u_mag_tot).abs() / u_elec_j.max(1e-300);
    let q_factor = (f.im.abs() > 1e-12 * f.re.abs().max(1.0)).then(|| f.re / (2.0 * f.im.abs()));

    let row = ModeEnergyRow { freq_hz: f.re, q_factor, u_elec_j, u_mag_vol_j, u_kin_ibc_j, u_ind_jj_j, v_mj, i_mj };
    (s_e, row, equipartition_error)
}

/// Simple Union-Find on vertex indices `0..nv`.
struct VertexDsu {
    parent: Vec<usize>,
}

impl VertexDsu {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect() }
    }
    fn find(&mut self, mut i: usize) -> usize {
        while self.parent[i] != i {
            self.parent[i] = self.parent[self.parent[i]];
            i = self.parent[i];
        }
        i
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            let (lo, hi) = (ra.min(rb), ra.max(rb));
            self.parent[hi] = lo;
        }
    }
}

/// Computes the 3D passive port capacitance matrix `C_port` (`J x J`, in F) on the island-tied `H^1` nodal
/// subspace and returns `(c_port_f, q_cap_red)` where `q_cap_red[j]` is the reduced Nedelec electrostatic
/// potential field `(G v_j)_red` across JJ / lumped port `j`.
fn compute_port_capacitance_and_static_cap(
    mesh: &Discretisation,
    problem: &Problem,
    asm: &EmwAssembly,
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    let n_ports = asm.jj_ports.len() + asm.lumped_ports.len();
    if n_ports == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    let order = asm.order;
    let n_nodes = if order == 2 { nv + ne } else { nv };

    // Identify Terminal - and Terminal + vertices for each JJ and lumped_port
    let mut term_minus: Vec<Vec<usize>> = Vec::with_capacity(n_ports);
    let mut term_plus: Vec<Vec<usize>> = Vec::with_capacity(n_ports);
    let mut port_b_reds: Vec<&[f64]> = Vec::with_capacity(n_ports);
    let mut jj_cross_edges: BTreeSet<(u32, u32)> = BTreeSet::new();

    let mut collect_single_terminals =
        |triangles: &[[u32; 3]], d_hat: [f64; 3], u_min: f64, u_max: f64, length_m: f64| -> (Vec<usize>, Vec<usize>) {
            let mut verts = BTreeSet::new();
            for &tri in triangles {
                for v in tri {
                    verts.insert(v);
                }
                for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
                    let pa = mesh.point(a as usize);
                    let pb = mesh.point(b as usize);
                    let ua = pa[0] * d_hat[0] + pa[1] * d_hat[1] + pa[2] * d_hat[2];
                    let ub = pb[0] * d_hat[0] + pb[1] * d_hat[1] + pb[2] * d_hat[2];
                    if (ua - ub).abs() > 0.2 * length_m {
                        jj_cross_edges.insert((a.min(b), a.max(b)));
                    }
                }
            }
            let mut tm = Vec::new();
            let mut tp = Vec::new();
            for &v in &verts {
                let p = mesh.point(v as usize);
                let u = p[0] * d_hat[0] + p[1] * d_hat[1] + p[2] * d_hat[2];
                if u <= u_min + 0.15 * length_m {
                    tm.push(v as usize);
                } else if u >= u_max - 0.15 * length_m {
                    tp.push(v as usize);
                }
            }
            if tm.is_empty() {
                if let Some(&v) = verts.first() {
                    tm.push(v as usize);
                }
            }
            if tp.is_empty() {
                if let Some(&v) = verts.last() {
                    tp.push(v as usize);
                }
            }
            (tm, tp)
        };

    for port in &asm.jj_ports {
        let (tm, tp) = collect_single_terminals(&port.triangles, port.d_hat, port.u_min, port.u_max, port.length_m);
        term_minus.push(tm);
        term_plus.push(tp);
        port_b_reds.push(&port.b_red);
    }
    for lp in &asm.lumped_ports {
        if let Some((ref comp_l, ref comp_r)) = lp.cpw_slots {
            // Balanced CPW even mode: center trace (inner edges u_L_max, u_R_min) is +, ground planes (outer edges) are -
            let span_of = |tris: &[[u32; 3]]| -> (f64, f64) {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for &tri in tris {
                    for v in tri {
                        let p = mesh.point(v as usize);
                        let u = p[0] * lp.d_hat[0] + p[1] * lp.d_hat[1] + p[2] * lp.d_hat[2];
                        lo = lo.min(u);
                        hi = hi.max(u);
                    }
                }
                (lo, hi)
            };
            let (ul0, ul1) = span_of(comp_l);
            let (ur0, ur1) = span_of(comp_r);
            let (tm_l, tp_l) = collect_single_terminals(comp_l, lp.d_hat, ul0, ul1, (ul1 - ul0).max(1e-15));
            let (tp_r, tm_r) = collect_single_terminals(comp_r, lp.d_hat, ur0, ur1, (ur1 - ur0).max(1e-15));
            let mut tm = tm_l;
            tm.extend(tm_r);
            let mut tp = tp_l;
            tp.extend(tp_r);
            term_minus.push(tm);
            term_plus.push(tp);
        } else {
            let (tm, tp) = collect_single_terminals(&lp.triangles, lp.d_hat, lp.u_min, lp.u_max, lp.length_m);
            term_minus.push(tm);
            term_plus.push(tp);
        }
        port_b_reds.push(&lp.b_red);
    }

    let mut dsu = VertexDsu::new(nv);
    for j in 0..n_ports {
        let m0 = term_minus[j][0];
        for &v in &term_minus[j][1..] {
            dsu.union(m0, v);
        }
        let p0 = term_plus[j][0];
        for &v in &term_plus[j][1..] {
            dsu.union(p0, v);
        }
    }

    let would_short_port = |dsu: &mut VertexDsu, a: usize, b: usize| -> bool {
        let (ra, rb) = (dsu.find(a), dsu.find(b));
        if ra == rb {
            return false;
        }
        for j in 0..n_ports {
            let rm = dsu.find(term_minus[j][0]);
            let rp = dsu.find(term_plus[j][0]);
            if (ra == rm && rb == rp) || (ra == rp && rb == rm) {
                return true;
            }
        }
        false
    };

    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();

    let mut sheet_cond_edges = Vec::new();
    let mut ext_pec_edges = Vec::new();
    let mut ext_pec_verts = BTreeSet::new();
    let mut cond_verts = vec![false; nv];

    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let is_sheet_cond = problem.pec.contains(&l) || ibc_labels.contains(&l);
        let is_ext_pec = exterior[f]
            && !periodic_labels.contains(&l)
            && !problem.pmc.contains(&l)
            && !ibc_labels.contains(&l)
            && !jj_labels.contains(&l)
            && !lp_labels.contains(&l);
        if !is_sheet_cond && !is_ext_pec {
            continue;
        }
        let fv = mesh.facet(f);
        for &v in fv {
            cond_verts[v as usize] = true;
            if is_ext_pec {
                ext_pec_verts.insert(v as usize);
            }
        }
        for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
            let key = (a.min(b), a.max(b));
            if jj_cross_edges.contains(&key) {
                continue;
            }
            if is_sheet_cond {
                sheet_cond_edges.push((key.0 as usize, key.1 as usize));
            } else {
                ext_pec_edges.push((key.0 as usize, key.1 as usize));
            }
        }
    }
    for &v in term_minus.iter().flatten().chain(term_plus.iter().flatten()) {
        cond_verts[v] = true;
    }

    for &(a, b) in &sheet_cond_edges {
        if !would_short_port(&mut dsu, a, b) {
            dsu.union(a, b);
        }
    }
    for &(a, b) in &ext_pec_edges {
        if !would_short_port(&mut dsu, a, b) {
            dsu.union(a, b);
        }
    }

    // Build island-tied nodal Reduction
    let mut rels = asm.periodic_node_relations.clone();
    for v in 0..nv {
        let r = dsu.find(v);
        if r != v {
            rels.push((v, r, c(1.0, 0.0)));
        }
    }
    let mut fixed_bubbles = Vec::new();
    if order == 2 {
        for (e, &(a, b)) in mesh.edges().iter().enumerate() {
            let (ua, ub) = (a as usize, b as usize);
            if cond_verts[ua] && cond_verts[ub] && dsu.find(ua) == dsu.find(ub) {
                fixed_bubbles.push(nv + e);
            }
        }
    }

    let red_island = Reduction::build(n_nodes, &rels, &fixed_bubbles);
    let p_island = red_island.prolong().real_part();
    let g_full = gradient_full(mesh, &asm.ned);
    let g_island = spmm(&g_full, &p_island);
    let g_island_t = g_island.transpose();
    let m0_g = spmm(&asm.m0_full, &g_island);
    let mut l_isl = spmm(&g_island_t, &m0_g);

    // Stabilize any zero-diagonal entries of l_isl
    let n_isl = l_isl.nrows();
    let diag = l_isl.diagonal();
    let max_diag = diag.iter().copied().fold(0.0_f64, f64::max).max(1.0);
    if diag.iter().any(|&d| d <= 1e-18 * max_diag) {
        let reg_diag: Vec<f64> =
            diag.iter().map(|&d| if d <= 1e-18 * max_diag { 1e-15 * max_diag } else { 0.0 }).collect();
        l_isl = l_isl.add_scaled(&RCsr::from_diag(&reg_diag), 1.0);
    }

    let mut ground_roots = BTreeSet::new();
    for &v in &ext_pec_verts {
        let idx = red_island.index_of[dsu.find(v)];
        if idx != usize::MAX {
            ground_roots.insert(idx);
        }
    }

    let mut port_minus_idx = Vec::with_capacity(n_ports);
    let mut port_plus_idx = Vec::with_capacity(n_ports);
    for j in 0..n_ports {
        let im = red_island.index_of[dsu.find(term_minus[j][0])];
        let ip = red_island.index_of[dsu.find(term_plus[j][0])];
        port_minus_idx.push(im);
        port_plus_idx.push(ip);
        // Ensure at least one terminal of port j is not locked to ground if both touched exterior facets
        if ground_roots.contains(&im) && ground_roots.contains(&ip) {
            ground_roots.remove(&ip);
        }
    }

    let mut phi_ports: Vec<Vec<f64>> = Vec::with_capacity(n_ports);
    let mut q_cap_red: Vec<Vec<f64>> = Vec::with_capacity(n_ports);

    for j in 0..n_ports {
        let im = port_minus_idx[j];
        let ip = port_plus_idx[j];
        let mut base_dirichlet: BTreeMap<usize, f64> = BTreeMap::new();
        for &g in &ground_roots {
            base_dirichlet.insert(g, 0.0);
        }
        // Short other ports k != j (V_{k,+} = V_{k,-} = 0 if grounded, or 0)
        for k in 0..n_ports {
            if k == j {
                continue;
            }
            let (km, kp) = (port_minus_idx[k], port_plus_idx[k]);
            if km != im && km != ip {
                base_dirichlet.entry(km).or_insert(0.0);
            }
            if kp != im && kp != ip {
                base_dirichlet.entry(kp).or_insert(0.0);
            }
        }

        let solve_with = |dir: &BTreeMap<usize, f64>| -> Result<Vec<f64>, String> {
            let fixed: Vec<usize> = dir.keys().copied().collect();
            let vals: Vec<f64> = dir.values().copied().collect();
            let rhs = vec![0.0; n_isl];
            solve_dirichlet(&l_isl, &rhs, &fixed, &vals)
        };

        let m_grounded = base_dirichlet.contains_key(&im) || base_dirichlet.is_empty();
        let p_grounded = base_dirichlet.contains_key(&ip);

        let phi_j = if m_grounded && !p_grounded {
            let mut dir = base_dirichlet.clone();
            dir.insert(im, 0.0);
            dir.insert(ip, 1.0);
            solve_with(&dir)?
        } else if p_grounded && !m_grounded {
            let mut dir = base_dirichlet.clone();
            dir.insert(ip, 0.0);
            dir.insert(im, -1.0);
            solve_with(&dir)?
        } else {
            // Both terminals float relative to the grounded exterior shield:
            // minimize electrostatic energy subject to V_+ - V_- = 1 (zero net charge on the floating dipole)
            let mut dir_p = base_dirichlet.clone();
            dir_p.insert(ip, 1.0);
            dir_p.insert(im, 0.0);
            let phi_p = solve_with(&dir_p)?;

            let mut dir_m = base_dirichlet.clone();
            dir_m.insert(ip, 0.0);
            dir_m.insert(im, 1.0);
            let phi_m = solve_with(&dir_m)?;

            let l_phi_p = l_isl.matvec(&phi_p);
            let l_phi_m = l_isl.matvec(&phi_m);
            let l_pp = dot_real(&phi_p, &l_phi_p);
            let l_mm = dot_real(&phi_m, &l_phi_m);
            let l_pm = dot_real(&phi_p, &l_phi_m);
            let denom = l_pp + l_mm + 2.0 * l_pm;
            let gamma = if denom > 1e-18 * (l_pp + l_mm) { ((l_mm + l_pm) / denom).clamp(0.0, 1.0) } else { 0.5 };
            phi_p.iter().zip(&phi_m).map(|(&vp, &vm)| gamma * vp - (1.0 - gamma) * vm).collect()
        };

        let e_full = g_island.matvec(&phi_j);
        let mut e_red: Vec<f64> = asm.red.roots.iter().map(|&r| e_full[r]).collect();
        let vj = dot_real(port_b_reds[j], &e_red);
        if vj.abs() > 1e-12 {
            for x in &mut e_red {
                *x /= vj;
            }
        }
        phi_ports.push(phi_j);
        q_cap_red.push(e_red);
    }

    let mut c_port_f = vec![vec![0.0; n_ports]; n_ports];
    let l_phis: Vec<Vec<f64>> = phi_ports.iter().map(|p| l_isl.matvec(p)).collect();
    for j in 0..n_ports {
        for k in j..n_ports {
            let c_jk = 0.5 * EPS0 * (dot_real(&phi_ports[j], &l_phis[k]) + dot_real(&phi_ports[k], &l_phis[j]));
            c_port_f[j][k] = c_jk;
            c_port_f[k][j] = c_jk;
        }
    }

    Ok((c_port_f, q_cap_red))
}

/// Builds a real `M_0`-orthogonal projector off the gradient matrix `g_proj_red` (`nred x n_grad`).
fn project_off_gradients(m0_red: &RCsr, g_red: &RCsr, vecs: &mut [Vec<f64>]) -> Result<(), String> {
    if g_red.ncols() == 0 || vecs.is_empty() {
        return Ok(());
    }
    let gt = g_red.transpose();
    let s_mat = spmm(&gt, &spmm(m0_red, g_red));
    let solver = DirectSolver::new(&s_mat, true)?;
    for x in vecs.iter_mut() {
        let mx = m0_red.matvec(x);
        let g = gt.matvec(&mx);
        let y = solver.solve(&g);
        let gy = g_red.matvec(&y);
        for (xi, &gyi) in x.iter_mut().zip(&gy) {
            *xi -= gyi;
        }
    }
    Ok(())
}

/// Two-pass Modified Gram-Schmidt orthonormalization in the `M_0` inner product.
fn append_m0_orthonormal(m0: &RCsr, basis: &mut Vec<Vec<f64>>, m0_basis: &mut Vec<Vec<f64>>, cand: &[f64]) -> bool {
    let mut v = cand.to_vec();
    let mv0 = m0.matvec(&v);
    let norm0 = dot_real(&v, &mv0).max(0.0).sqrt();
    if norm0 <= 1e-30 {
        return false;
    }
    for _pass in 0..2 {
        for (q, mq) in basis.iter().zip(m0_basis.iter()) {
            let proj = dot_real(&v, mq);
            for (vi, &qi) in v.iter_mut().zip(q) {
                *vi -= proj * qi;
            }
        }
    }
    let mv = m0.matvec(&v);
    let norm = dot_real(&v, &mv).max(0.0).sqrt();
    if norm <= 1e-11 * norm0 {
        return false;
    }
    let inv = 1.0 / norm;
    for vi in &mut v {
        *vi *= inv;
    }
    let mv_unit: Vec<f64> = mv.into_iter().map(|x| x * inv).collect();
    basis.push(v);
    m0_basis.push(mv_unit);
    true
}

/// Assembles the passive Craig-Bampton reduced-order model `CircuitRom` (`Q = orth_{M_0}([Q_static, Q_dyn])`).
pub fn build_circuit_rom(
    mesh: &Discretisation,
    problem: &Problem,
    asm: &EmwAssembly,
    solved_vectors: &[Vec<C>],
    freq_indices: &[(C, usize)],
    modes: Vec<ModeEnergyRow>,
) -> Result<CircuitRom, String> {
    let mut ports: Vec<RomPort> = asm
        .jj_ports
        .iter()
        .map(|p| RomPort {
            name: p.bc.name.clone(),
            kind: "lumped_jj".into(),
            l_j_h: (p.bc.l_j > 0.0).then_some(p.bc.l_j),
            e_j_hz: (p.bc.e_j_hz > 0.0).then_some(p.bc.e_j_hz),
            c_j_f: p.bc.c_j,
            e_l_hz: (p.bc.e_l_hz > 0.0).then_some(p.bc.e_l_hz),
            phi_ext_rad: p.bc.phi_ext_rad,
            n_g: p.bc.n_g,
            squid_asymmetry: p.bc.squid_asymmetry,
            z0_ohm: None,
            direction: p.bc.direction,
            formulation: p.bc.formulation.clone(),
        })
        .collect();
    for lp in &asm.lumped_ports {
        ports.push(RomPort {
            name: lp.bc.name.clone(),
            kind: "lumped_port".into(),
            l_j_h: None,
            e_j_hz: None,
            c_j_f: 0.0,
            e_l_hz: None,
            phi_ext_rad: 0.0,
            n_g: 0.0,
            squid_asymmetry: None,
            z0_ohm: Some(lp.bc.z0),
            direction: lp.bc.direction,
            formulation: if lp.bc.mode == "cpw_even" { "cpw_even".into() } else { "rank1".into() },
        });
    }

    let (c_port_f, mut q_cap_list) = compute_port_capacitance_and_static_cap(mesh, problem, asm)?;

    // Regularized full-domain (including PML) matrices for non-singular gradient projection & static lifts
    let k_passive_reg =
        if asm.k_ibc_red.nnz() > 0 { asm.k0_red.add_scaled(&asm.k_ibc_red, 1.0) } else { asm.k0_red.clone() };
    let m0_red = &asm.m0_red;
    let g_red_r = asm.g_red.real_part();

    // Project electrostatic fields off interior gradients g_red
    let _ = project_off_gradients(m0_red, &g_red_r, &mut q_cap_list);

    // Compute static inductive lifts u_j = (K_passive + eps_reg M_0)^-1 b_j projected off gradients
    let mut q_ind_list: Vec<Vec<f64>> = Vec::new();
    if !asm.jj_ports.is_empty() || !asm.lumped_ports.is_empty() {
        let eps_reg = (1e-4 * asm.sigma).max(1e-6);
        let a_reg = k_passive_reg.add_scaled(m0_red, eps_reg);
        if let Ok(reg_solver) = DirectSolver::new(&a_reg, true) {
            for port in &asm.jj_ports {
                if port.b_red.iter().any(|&v| v.abs() > 0.0) {
                    q_ind_list.push(reg_solver.solve(&port.b_red));
                }
            }
            for lp in &asm.lumped_ports {
                if lp.b_red.iter().any(|&v| v.abs() > 0.0) {
                    q_ind_list.push(reg_solver.solve(&lp.b_red));
                }
            }
            let _ = project_off_gradients(m0_red, &g_red_r, &mut q_ind_list);
        }
    }

    let mut q_cols: Vec<Vec<f64>> = Vec::new();
    let mut m0_q_cols: Vec<Vec<f64>> = Vec::new();

    for cand in q_cap_list.iter().chain(q_ind_list.iter()) {
        append_m0_orthonormal(m0_red, &mut q_cols, &mut m0_q_cols, cand);
    }
    let r_static = q_cols.len();

    for &(_, idx) in freq_indices {
        let vec_c = &solved_vectors[idx];
        let big = vec_c.iter().copied().fold(c(0.0, 0.0), |m, z| if z.norm_sqr() > m.norm_sqr() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let re_vec: Vec<f64> = vec_c.iter().map(|&z| (z * rot).re).collect();
        let im_vec: Vec<f64> = vec_c.iter().map(|&z| (z * rot).im).collect();
        append_m0_orthonormal(m0_red, &mut q_cols, &mut m0_q_cols, &re_vec);
        let im_norm2: f64 = im_vec.iter().map(|x| x * x).sum();
        let re_norm2: f64 = re_vec.iter().map(|x| x * x).sum();
        if im_norm2 > 1e-12 * re_norm2 {
            append_m0_orthonormal(m0_red, &mut q_cols, &mut m0_q_cols, &im_vec);
        }
    }
    let r_tot = q_cols.len();
    let r_dyn = r_tot - r_static;

    // Project K_r and M_r strictly from physical (non-PML) cells (`k_phys_red` + `k_ibc_red`, `m_phys_red`)
    let k_passive_phys =
        if asm.k_ibc_red.nnz() > 0 { asm.k_phys_red.add_scaled(&asm.k_ibc_red, 1.0) } else { asm.k_phys_red.clone() };
    let kq_cols: Vec<Vec<f64>> = q_cols.iter().map(|q| k_passive_phys.matvec(q)).collect();
    let mq_cols: Vec<Vec<f64>> = q_cols.iter().map(|q| asm.m_phys_red.matvec(q)).collect();
    let mut k_r = vec![vec![0.0; r_tot]; r_tot];
    let mut m_r = vec![vec![0.0; r_tot]; r_tot];
    for a in 0..r_tot {
        for b in a..r_tot {
            let kab = 0.5 * (dot_real(&q_cols[a], &kq_cols[b]) + dot_real(&q_cols[b], &kq_cols[a]));
            let mab = 0.5 * (dot_real(&q_cols[a], &mq_cols[b]) + dot_real(&q_cols[b], &mq_cols[a]));
            k_r[a][b] = kab;
            k_r[b][a] = kab;
            m_r[a][b] = mab;
            m_r[b][a] = mab;
        }
    }

    let c_tot_rom = if asm.c_ibc_red.nnz() > 0 && asm.c_port_red.nnz() > 0 {
        asm.c_ibc_red.add_scaled(&asm.c_port_red, 1.0)
    } else if asm.c_ibc_red.nnz() > 0 {
        asm.c_ibc_red.clone()
    } else {
        asm.c_port_red.clone()
    };
    let c_r = if c_tot_rom.nnz() > 0 {
        let cq_cols: Vec<Vec<f64>> = q_cols.iter().map(|q| c_tot_rom.matvec(q)).collect();
        let mut cr = vec![vec![0.0; r_tot]; r_tot];
        for a in 0..r_tot {
            for b in a..r_tot {
                let cab = 0.5 * (dot_real(&q_cols[a], &cq_cols[b]) + dot_real(&q_cols[b], &cq_cols[a]));
                cr[a][b] = cab;
                cr[b][a] = cab;
            }
        }
        cr
    } else {
        Vec::new()
    };

    let n_jj = asm.jj_ports.len();
    let n_ports = n_jj + asm.lumped_ports.len();
    let mut b_r = vec![vec![0.0; n_ports]; r_tot];
    for (a, qa) in q_cols.iter().enumerate() {
        for (j, port) in asm.jj_ports.iter().enumerate() {
            b_r[a][j] = dot_real(qa, &port.b_red);
        }
        for (j, lp) in asm.lumped_ports.iter().enumerate() {
            b_r[a][n_jj + j] = dot_real(qa, &lp.b_red);
        }
    }

    Ok(CircuitRom { ports, r_static, r_dyn, k_r, c_r, m_r, b_r, c_port_f, modes })
}

/// Evaluates the linearized resonant eigenfrequencies (Hz, ascending) of a `CircuitRom` by solving the
/// `r x r` generalized eigenvalue problem
/// `(K_r + sum_j (mu_0 / L_tot,j) b_{r,j} b_{r,j}^T) z = k_0^2 (M_r + sum_j (C_{J,j} / eps_0) b_{r,j} b_{r,j}^T) z`.
pub fn rom_eigenfrequencies(rom: &CircuitRom) -> Result<Vec<f64>, String> {
    let r = rom.k_r.len();
    if r == 0 {
        return Ok(Vec::new());
    }
    let mut k_tot = vec![0.0; r * r];
    let mut m_tot = vec![0.0; r * r];
    for i in 0..r {
        for j in 0..r {
            k_tot[i * r + j] = rom.k_r[i][j];
            m_tot[i * r + j] = rom.m_r[i][j];
        }
    }
    for (p, port) in rom.ports.iter().enumerate() {
        let l_j = port.l_j_h.unwrap_or(0.0);
        let e_l = port.e_l_hz.unwrap_or(0.0);
        let inv_l_tot = jj_inv_l_tot(l_j, e_l);
        let alpha_k = MU0 * inv_l_tot;
        let alpha_m = port.c_j_f / EPS0;
        for i in 0..r {
            let bi = rom.b_r[i][p];
            for j in 0..r {
                let bj = rom.b_r[j][p];
                k_tot[i * r + j] += alpha_k * bi * bj;
                m_tot[i * r + j] += alpha_m * bi * bj;
            }
        }
    }
    let (evals, _) = dense_gen_eigh(&k_tot, &m_tot, r)?;
    let mut freqs: Vec<f64> =
        evals.into_iter().filter(|&lam| lam > 0.0).map(|lam| lam.sqrt() * (C0 / (2.0 * PI))).collect();
    freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Ok(freqs)
}

/// Evaluates the complex damped resonant eigenfrequencies (Hz, positive real part, ascending by `Re(f)`)
/// of a `CircuitRom` by solving the exact `2r x 2r` quadratic eigenvalue problem
/// `(K_tot + i k_0 C_tot - k_0^2 M_tot) z = 0`.
pub fn rom_complex_eigenfrequencies(rom: &CircuitRom) -> Result<Vec<C>, String> {
    let r = rom.k_r.len();
    if r == 0 {
        return Ok(Vec::new());
    }
    let mut k_tot = vec![0.0; r * r];
    let mut m_tot = vec![0.0; r * r];
    let mut c_tot = vec![0.0; r * r];
    let has_c_r = rom.c_r.len() == r && rom.c_r.iter().all(|row| row.len() == r);
    for i in 0..r {
        for j in 0..r {
            k_tot[i * r + j] = rom.k_r[i][j];
            m_tot[i * r + j] = rom.m_r[i][j];
            if has_c_r {
                c_tot[i * r + j] = rom.c_r[i][j];
            }
        }
    }
    for (p, port) in rom.ports.iter().enumerate() {
        let l_j = port.l_j_h.unwrap_or(0.0);
        let e_l = port.e_l_hz.unwrap_or(0.0);
        let inv_l_tot = jj_inv_l_tot(l_j, e_l);
        let alpha_k = MU0 * inv_l_tot;
        let alpha_m = port.c_j_f / EPS0;
        let alpha_c = if !has_c_r { port.z0_ohm.filter(|&z| z > 0.0).map_or(0.0, |z0| MU0 * C0 / z0) } else { 0.0 };
        if alpha_k == 0.0 && alpha_m == 0.0 && alpha_c == 0.0 {
            continue;
        }
        for i in 0..r {
            let bi = rom.b_r[i][p];
            for j in 0..r {
                let bj = rom.b_r[j][p];
                k_tot[i * r + j] += alpha_k * bi * bj;
                m_tot[i * r + j] += alpha_m * bi * bj;
                c_tot[i * r + j] += alpha_c * bi * bj;
            }
        }
    }

    let fro_k = k_tot.iter().map(|v| v * v).sum::<f64>().sqrt();
    let fro_m = m_tot.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
    let s_k = (fro_k / fro_m).sqrt().max(1.0);
    let inv_sk = 1.0 / s_k;
    let inv_sk2 = inv_sk * inv_sk;

    let n2 = 2 * r;
    let mut a = Dense::<C>::zeros(n2, n2);
    let mut b = Dense::<C>::zeros(n2, n2);
    for i in 0..r {
        a[(i, r + i)] = c(1.0, 0.0);
        b[(i, i)] = c(1.0, 0.0);
        for j in 0..r {
            a[(r + i, j)] = c(k_tot[i * r + j] * inv_sk2, 0.0);
            a[(r + i, r + j)] = c(0.0, c_tot[i * r + j] * inv_sk);
            b[(r + i, r + j)] = c(m_tot[i * r + j], 0.0);
        }
    }

    let (eigs, _) =
        generalized_eig_c64(&a, &b).map_err(|e| format!("CircuitRom complex QEP eigensolve failed: {e}"))?;
    let scale_hz = s_k * (C0 / (2.0 * PI));
    let mut freqs: Vec<C> = eigs
        .into_iter()
        .filter(|z| z.re.is_finite() && z.im.is_finite() && z.re > 0.0)
        .map(|mu| mu * scale_hz)
        .collect();
    freqs.sort_by(|a, b| a.re.total_cmp(&b.re));
    Ok(freqs)
}

/// Builds `TlsMeta` for dielectric TLS surface participation analysis on a 3D `ems` result.
pub fn build_tls_meta(mesh: &Discretisation, problem: &Problem) -> TlsMeta {
    let mut domain_eps_r = BTreeMap::new();
    for d in mesh.domain_labels() {
        let eps = problem
            .materials
            .get(d, "eps_r")
            .or_else(|| {
                problem.materials.aniso.get(&d).and_then(|a| {
                    a.permittivity.as_ref().map(|t| {
                        let g = femcore::material::tensor_global(t, a.euler_deg);
                        (g[0][0] + g[1][1] + g[2][2]) / 3.0
                    })
                })
            })
            .unwrap_or(1.0);
        domain_eps_r.insert(d, eps);
    }

    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let wp_labels: Vec<i32> = problem.waveports.iter().flat_map(|b| b.labels.iter().copied()).collect();

    let mut pec_set = BTreeSet::new();
    let mut cond_set = BTreeSet::new();
    let mut sheet_set = BTreeSet::new();

    for &l in &problem.pec {
        pec_set.insert(l);
        cond_set.insert(l);
    }
    for &l in &ibc_labels {
        cond_set.insert(l);
    }

    let mut seen_interior_tris = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        if exterior[f] {
            if !periodic_labels.contains(&l)
                && !problem.pmc.contains(&l)
                && !ibc_labels.contains(&l)
                && !jj_labels.contains(&l)
                && !lp_labels.contains(&l)
                && !wp_labels.contains(&l)
            {
                pec_set.insert(l);
                cond_set.insert(l);
            }
        } else {
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if seen_interior_tris.insert(key) && (problem.pec.contains(&l) || ibc_labels.contains(&l)) {
                sheet_set.insert(l);
            }
        }
    }

    let mut sheet_thickness_m = BTreeMap::new();
    for bc in &problem.surface_impedance {
        if bc.thickness > 0.0 {
            for &l in &bc.labels {
                sheet_thickness_m.insert(l, bc.thickness);
            }
        }
    }

    TlsMeta {
        domain_eps_r,
        domain_names: problem.domain_names.clone(),
        conductor_labels: cond_set.into_iter().collect(),
        sheet_labels: sheet_set.into_iter().collect(),
        sheet_thickness_m,
        pec_labels: pec_set.into_iter().collect(),
        interfaces: problem.tls_interfaces.clone(),
    }
}
