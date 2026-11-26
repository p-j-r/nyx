// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Driven frequency-domain multi-port piezoelectric solver (`Study::Driven`), computing
//! exact Schur-complement admittance $\mathbf{Y}(\omega)$, impedance $\mathbf{Z}(\omega)$,
//! scattering $\mathbf{S}(\omega)$, and clamped capacitance $\mathbf{C}_0$ across both
//! [`SimplexMesh`] and native [`MixedDiscretisation`] (`Tri`/`Quad`, `Tet`/`Prism`/`Hex`),
//! accelerated by multi-port Hermite Tangential Galerkin PROM (`femcore::linalg::prom`).
//!
//! Eliminating interior mechanical + potential DOFs $I = (u_I, \phi_I)$ for unit terminal
//! voltage excitations $\mathbf{\Phi}_P = \mathbf{e}_q$ yields
//! $$\mathbf{Y}(\omega) = -i\omega\,\mathbf{S}_P(\omega)
//! = -i\omega\big(\mathbf{K}_{PP} - \mathbf{T}_{PI}\,\mathbf{A}_{II}(\omega)^{-1}\,\mathbf{T}_{IP}\big),$$
//! using non-conjugate $\mathbf{T}_{PI} = \mathbf{T}_{IP}^T$ for symmetric problems so
//! dielectric $\tan\delta$, viscoelasticity, and PML preserve exact reciprocity and passivity.
#![allow(clippy::needless_range_loop)]

use crate::piezo::{push_eigen_mode_fields_mixed, push_eigen_mode_fields_simplex};
use crate::piezo_asm::{
    assemble_piezo_driven_mixed, assemble_piezo_driven_simplex, sub_csr, PiezoConstraints, PiezoDrivenFull,
};
use crate::problem::{Problem, Study};
use femcore::fem_mixed::{LagrangeSpace, MixedLagrangeDofMap};
use femcore::linalg::csr::{c, CCoo, CCsr, RCsr, C};
use femcore::linalg::driver::DirectSolver;
use femcore::linalg::equilibrate::{ruiz_symmetric_terms, scale_symmetric, RuizTerm, DEFAULT_RUIZ_PASSES};
use femcore::linalg::prom::{
    run_greedy_prom_sweep, y_to_z_and_s, ImageProm, PromLinearSolver, PromOperator, SparseDrivenSolver,
};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::mesh::SimplexMesh;
use femcore::periodic::congruence_fast;
use nyx_physics::{
    cell_fields_of, DrivenPortData, FieldFrame, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use std::f64::consts::PI;
use std::ops::Range;
use std::time::Instant;

fn pad_csr_to(a: &CCsr, n_tot: usize) -> CCsr {
    if a.nnz() == 0 {
        return CCsr::zeros(n_tot, n_tot);
    }
    let nr_old = a.nrows();
    let mut indptr = Vec::with_capacity(n_tot + 1);
    indptr.extend_from_slice(a.indptr());
    let nnz = a.nnz();
    for _ in nr_old..n_tot {
        indptr.push(nnz);
    }
    let indices: Vec<usize> = a.indices().iter().map(|&c| c as usize).collect();
    CCsr::from_raw(n_tot, n_tot, indptr, indices, a.values().to_vec())
}

fn transpose_ccsr_full(a: &CCsr) -> CCsr {
    let mut coo = CCoo::new(a.ncols(), a.nrows());
    for i in 0..a.nrows() {
        let (cols, vals) = a.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            coo.push(j as usize, i, v);
        }
    }
    coo.to_csr()
}

fn restrict_transpose_rhs(red: &femcore::periodic::Reduction, f: &[C]) -> Vec<C> {
    let nr = red.n_reduced();
    let mut out = vec![c(0.0, 0.0); nr];
    for i in 0..red.n {
        let p = red.index_of[i];
        if p != usize::MAX {
            out[p] += red.ratio[i] * f[i];
        }
    }
    out
}

/// Reduced equilibrated affine piezoelectric frequency-domain problem ready for PROM or HDM sweeps.
pub struct PiezoDrivenAssembly {
    pub constraints: PiezoConstraints,
    pub nred: usize,
    pub f_min: f64,
    pub f_max: f64,
    pub n_points: usize,
    pub freqs_hz: Option<Vec<f64>>,
    pub prom_tol: Option<f64>,
    pub max_hdm_solves: usize,
    pub strict_prom_tol: bool,
    pub seed_freqs_hz: Vec<f64>,
    pub solver: String,
    pub blr_tol: f64,
    pub port_names: Vec<String>,
    pub z0_ports: Vec<f64>,
    pub d_ruiz: Vec<f64>,
    pub terms_hat: Vec<CCsr>,
    pub k_pp: Vec<Vec<C>>,
    pub k_pp_neg_sub: CCsr,
    pub b_drive_unscaled: Vec<Vec<C>>,
    pub b_recv_unscaled: Option<Vec<Vec<C>>>,
    pub b_drive_hat: Vec<Vec<C>>,
    pub b_recv_hat: Option<Vec<Vec<C>>>,
    pub is_sym: bool,
    pub is_complex: bool,
    pub has_pml: bool,
    pub has_damping: bool,
    pub has_interior: bool,
    pub pml_domains: Vec<i32>,
    pub m1_phys_all: Option<RCsr>,
    pub m1_phys_pml: Option<RCsr>,
    pub block_splits: Vec<Range<usize>>,
    pub diagnostics: Vec<String>,
}

impl PiezoDrivenAssembly {
    fn from_full(full: PiezoDrivenFull, spec: &Problem) -> Result<Self, String> {
        let Study::Driven {
            f_min: study_f_min,
            f_max: study_f_max,
            n_points: study_n_points,
            freqs_hz,
            prom_tol,
            max_hdm_solves,
            strict_prom_tol,
            seed_freqs_hz,
            solver,
            blr_tol,
        } = spec.study.clone()
        else {
            return Err("assemble_piezo_driven requires Study::Driven".into());
        };
        let (f_min, f_max, n_points) = if let Some(fs) = freqs_hz.as_ref().filter(|v| !v.is_empty()) {
            (fs[0], fs[fs.len() - 1], fs.len())
        } else {
            (study_f_min, study_f_max, study_n_points)
        };
        let f0 = 0.5 * (f_min + f_max);

        let constraints = full.constraints;
        let n_u_red = constraints.n_u_red;
        let nred = constraints.red_all.n_reduced();
        let p_ports = constraints.terminals.len();
        if p_ports == 0 {
            return Err("driven study requires at least one terminal electrode".into());
        }

        let k_ii_red = congruence_fast(&constraints.red_all, &full.k_full);
        let m_uu_red = congruence_fast(&constraints.red_u, &full.m_uu_full);
        let m_ii_red = pad_csr_to(&m_uu_red, nred);
        let c_ii_red = match &full.c_uu_full {
            Some(cuu) => pad_csr_to(&congruence_fast(&constraints.red_u, cuu), nred),
            None => CCsr::zeros(0, 0),
        };

        let is_sym = constraints.red_all.ratio.iter().all(|r| r.im.abs() <= 1e-12);
        let k_full_t = (!is_sym).then(|| transpose_ccsr_full(&full.k_full));

        let mut k_pp = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
        let mut b_drive_unscaled = Vec::with_capacity(p_ports);
        let mut b_recv_unscaled_vec = (!is_sym).then(|| Vec::with_capacity(p_ports));

        for q in 0..p_ports {
            let g_q = &constraints.terminal_lifts[q];
            let col_q = full.k_full.matvec(g_q);
            let mut bq = constraints.red_all.restrict_rhs(&col_q);
            for v in &mut bq {
                *v = -*v;
            }
            b_drive_unscaled.push(bq);

            if let Some(ref mut b_recv_out) = b_recv_unscaled_vec {
                let g_q_conj: Vec<C> = g_q.iter().map(|z| z.conj()).collect();
                let col_t_q = k_full_t.as_ref().unwrap().matvec(&g_q_conj);
                let mut br = restrict_transpose_rhs(&constraints.red_all, &col_t_q);
                for v in &mut br {
                    *v = -*v;
                }
                b_recv_out.push(br);
            }

            for p in 0..p_ports {
                let g_p = &constraints.terminal_lifts[p];
                let mut dot = c(0.0, 0.0);
                for (&gp_i, &cq_i) in g_p.iter().zip(&col_q) {
                    if gp_i.re != 0.0 || gp_i.im != 0.0 {
                        dot += gp_i.conj() * cq_i;
                    }
                }
                k_pp[p][q] = dot;
            }
        }

        let k_pp_neg_sub = if constraints.n_phi_red > 0 {
            sub_csr(&k_ii_red, n_u_red, nred, n_u_red, nred)
        } else {
            CCsr::zeros(0, 0)
        };

        let omega_mid = (2.0 * PI * f0).max(1.0);
        let mut ruiz_terms =
            vec![RuizTerm::Complex(&k_ii_red, 1.0), RuizTerm::Complex(&m_ii_red, omega_mid * omega_mid)];
        if c_ii_red.nnz() > 0 {
            ruiz_terms.push(RuizTerm::Complex(&c_ii_red, omega_mid));
        }
        let d_ruiz = ruiz_symmetric_terms(nred, &ruiz_terms, DEFAULT_RUIZ_PASSES);

        let k_hat = scale_symmetric(&k_ii_red, &d_ruiz);
        let c_hat = if c_ii_red.nnz() > 0 { scale_symmetric(&c_ii_red, &d_ruiz) } else { CCsr::zeros(0, 0) };
        let m_hat = scale_symmetric(&m_ii_red, &d_ruiz);
        let terms_hat = vec![k_hat, c_hat, m_hat];

        let scale_cvec = |v: &[C]| -> Vec<C> { v.iter().zip(&d_ruiz).map(|(&vi, &di)| vi * di).collect() };
        let b_drive_hat: Vec<Vec<C>> = b_drive_unscaled.iter().map(|v| scale_cvec(v)).collect();
        let b_recv_hat: Option<Vec<Vec<C>>> =
            b_recv_unscaled_vec.as_ref().map(|vs| vs.iter().map(|v| scale_cvec(v)).collect());

        let port_names: Vec<String> = constraints.terminals.iter().map(|t| t.name.clone()).collect();
        let z0_ports: Vec<f64> = constraints.terminals.iter().map(|t| t.z0_ohm).collect();

        let mut block_splits = Vec::new();
        if n_u_red > 0 {
            block_splits.push(0..n_u_red);
        }
        if constraints.n_phi_red > 0 {
            block_splits.push(n_u_red..nred);
        }

        let pml_domains: Vec<i32> = full.resolved_pml.domains.keys().copied().collect();

        Ok(Self {
            constraints,
            nred,
            f_min,
            f_max,
            n_points,
            freqs_hz,
            prom_tol,
            max_hdm_solves,
            strict_prom_tol,
            seed_freqs_hz,
            solver,
            blr_tol,
            port_names,
            z0_ports,
            d_ruiz,
            terms_hat,
            k_pp,
            k_pp_neg_sub,
            b_drive_unscaled,
            b_recv_unscaled: b_recv_unscaled_vec,
            b_drive_hat,
            b_recv_hat,
            is_sym,
            is_complex: full.is_complex,
            has_pml: full.has_pml,
            has_damping: full.has_damping,
            has_interior: full.has_interior,
            pml_domains,
            m1_phys_all: full.m1_phys_all,
            m1_phys_pml: full.m1_phys_pml,
            block_splits,
            diagnostics: full.diagnostics,
        })
    }

    fn compute_clamped_c0(&self) -> Result<Vec<Vec<f64>>, String> {
        let p_ports = self.port_names.len();
        let n_u_red = self.constraints.n_u_red;
        let n_phi_red = self.constraints.n_phi_red;
        let mut c0 = vec![vec![0.0; p_ports]; p_ports];

        if n_phi_red == 0 {
            for p in 0..p_ports {
                for q in 0..p_ports {
                    c0[p][q] = -self.k_pp[p][q].re;
                }
            }
            return Ok(c0);
        }

        let phi_solver = DirectSolver::new(&self.k_pp_neg_sub, false)?;
        let rhs_phi: Vec<Vec<C>> = self.b_drive_unscaled.iter().map(|b| b[n_u_red..].to_vec()).collect();
        let phi_sols = phi_solver.solve_many(&rhs_phi);
        for p in 0..p_ports {
            let b_recv_p = self.b_recv_unscaled.as_ref().map_or(&self.b_drive_unscaled[p], |v| &v[p]);
            let b_phi = &b_recv_p[n_u_red..];
            for q in 0..p_ports {
                let mut inner = c(0.0, 0.0);
                for (&bp, &pq) in b_phi.iter().zip(&phi_sols[q]) {
                    inner += bp * pq;
                }
                c0[p][q] = -(self.k_pp[p][q] - inner).re;
            }
        }
        Ok(c0)
    }
}

impl PromOperator for PiezoDrivenAssembly {
    fn n_dofs(&self) -> usize {
        self.nred
    }

    fn n_ports(&self) -> usize {
        self.port_names.len()
    }

    fn terms_hat(&self) -> &[CCsr] {
        &self.terms_hat
    }

    fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)> {
        vec![
            (c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
            (c(0.0, omega), c(0.0, 1.0), c(0.0, 0.0)),
            (c(-omega * omega, 0.0), c(-2.0 * omega, 0.0), c(-2.0, 0.0)),
        ]
    }

    fn d_ruiz(&self) -> &[f64] {
        &self.d_ruiz
    }

    fn port_b_hat(&self, port: usize) -> Vec<C> {
        self.b_drive_hat[port].clone()
    }

    fn port_b_recv_hat(&self, port: usize) -> Option<Vec<C>> {
        self.b_recv_hat.as_ref().map(|v| v[port].clone())
    }

    fn port_gamma(&self, _port: usize, _omega: f64) -> (C, C, C) {
        (c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0))
    }

    fn is_symmetric(&self) -> bool {
        self.is_sym
    }

    fn block_splits(&self) -> Option<&[Range<usize>]> {
        (!self.block_splits.is_empty()).then_some(&self.block_splits)
    }
}

struct SweepCoreResult {
    driven_data: DrivenPortData,
    c0_matrix: Vec<Vec<f64>>,
    symbolic_count: usize,
    backend_name: String,
    rom_sweep_seconds: f64,
    f_rep_hz: f64,
    rep_x_red: Vec<C>,
}

fn run_piezo_driven_core(asm: &PiezoDrivenAssembly) -> Result<SweepCoreResult, String> {
    let freqs_hz: Vec<f64> = if let Some(fs) = asm.freqs_hz.as_ref().filter(|v| !v.is_empty()) {
        fs.to_vec()
    } else {
        let n_pts = asm.n_points.max(1);
        if n_pts == 1 {
            vec![asm.f_min]
        } else {
            (0..n_pts).map(|k| asm.f_min + (asm.f_max - asm.f_min) * (k as f64) / ((n_pts - 1) as f64)).collect()
        }
    };
    let n_points = freqs_hz.len();
    let omegas: Vec<f64> = freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let p_ports = asm.port_names.len();

    let c0_matrix = asm.compute_clamped_c0()?;

    let mut z0_ohm = Vec::with_capacity(n_points);
    let mut s_matrix = Vec::with_capacity(n_points);
    let mut z_matrix = Vec::with_capacity(n_points);
    let mut y_matrix = Vec::with_capacity(n_points);

    let z0_pairs_row: Vec<[f64; 2]> = asm.z0_ports.iter().map(|&z| [z, 0.0]).collect();

    let hdm_solves;
    let symbolic_count;
    let backend_name;
    let rom_dim;
    let anchor_freqs_hz;
    let max_residual;
    let mut rom_sweep_seconds = 0.0_f64;
    let mut rep_freq_idx = n_points / 2;
    let mut rep_x_red = vec![c(0.0, 0.0); asm.nred];
    let mut max_ceff_metric = -1.0_f64;

    if let Some(tol) = asm.prom_tol {
        let mut solver = SparseDrivenSolver::new(&asm.solver, asm.blr_tol);
        const SWEEP_CHUNK: usize = 32;
        let summary = run_greedy_prom_sweep(
            asm,
            &mut solver,
            &freqs_hz,
            tol,
            asm.max_hdm_solves,
            asm.strict_prom_tol,
            &asm.seed_freqs_hz,
            |band, b_omegas, prom: &ImageProm<'_, PiezoDrivenAssembly>, xr_band| {
                let nb_pts = b_omegas.len();
                let t_post = Instant::now();
                let b_dim = prom.dim();
                let b_recv = prom.b_r_recv();
                let indices: Vec<usize> = (0..nb_pts).collect();
                let post_chunks = femcore::par::map_chunks(&indices, SWEEP_CHUNK, |chunk| {
                    chunk
                        .iter()
                        .map(|&k| {
                            let w = b_omegas[k];
                            let xr_k = &xr_band[k];
                            let mut y_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
                            for p in 0..p_ports {
                                for q in 0..p_ports {
                                    let mut inner = c(0.0, 0.0);
                                    for i in 0..b_dim {
                                        inner += b_recv[i][p] * xr_k[q][i];
                                    }
                                    y_k[p][q] = c(0.0, -w) * (asm.k_pp[p][q] - inner);
                                }
                            }
                            let metric = y_k[0][0].norm() / w.max(1.0);
                            let (s_out, z_out, y_out) = y_to_z_and_s(&y_k, &asm.z0_ports);
                            (s_out, z_out, y_out, metric)
                        })
                        .collect::<Vec<_>>()
                });
                let mut band_best_local = None;
                for (k_loc, (s_out, z_out, y_out, metric)) in post_chunks.into_iter().flatten().enumerate() {
                    if metric > max_ceff_metric {
                        max_ceff_metric = metric;
                        rep_freq_idx = band.start + k_loc;
                        band_best_local = Some(k_loc);
                    }
                    z0_ohm.push(z0_pairs_row.clone());
                    s_matrix.push(s_out);
                    z_matrix.push(z_out);
                    y_matrix.push(y_out);
                }
                if let Some(kl) = band_best_local {
                    rep_x_red = prom.reconstruct_unscaled(&xr_band[kl][0]);
                }
                rom_sweep_seconds += t_post.elapsed().as_secs_f64();
                Ok(())
            },
        )?;
        symbolic_count = solver.symbolic_count();
        backend_name = solver.backend_name();
        drop(solver);
        hdm_solves = summary.hdm_solves;
        rom_dim = summary.rom_dim;
        anchor_freqs_hz = summary.anchor_freqs_hz;
        max_residual = summary.max_residual;
        rom_sweep_seconds += summary.rom_sweep_seconds;
    } else {
        let mut solver = SparseDrivenSolver::new(&asm.solver, asm.blr_tol);
        let rhs_hat = asm.b_drive_hat.clone();
        let mut anchors = Vec::with_capacity(n_points);
        for (k, &w) in omegas.iter().enumerate() {
            femcore::progress::check()?;
            let a_mat = asm.assemble_matrix_at(w);
            let x_hat_cols = solver.update_and_solve(w, &a_mat, &rhs_hat)?;
            anchors.push(freqs_hz[k]);

            let mut y_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
            for p in 0..p_ports {
                let b_recv_p = asm.b_recv_hat.as_ref().map_or(&asm.b_drive_hat[p], |v| &v[p]);
                for q in 0..p_ports {
                    let mut inner = c(0.0, 0.0);
                    for (&bp, &xq) in b_recv_p.iter().zip(&x_hat_cols[q]) {
                        inner += bp * xq;
                    }
                    y_k[p][q] = c(0.0, -w) * (asm.k_pp[p][q] - inner);
                }
            }
            let metric = y_k[0][0].norm() / w.max(1.0);
            if metric > max_ceff_metric {
                max_ceff_metric = metric;
                rep_freq_idx = k;
                rep_x_red = x_hat_cols[0].iter().zip(&asm.d_ruiz).map(|(&xi, &di)| xi * di).collect();
            }
            let (s_out, z_out, y_out) = y_to_z_and_s(&y_k, &asm.z0_ports);
            z0_ohm.push(z0_pairs_row.clone());
            s_matrix.push(s_out);
            z_matrix.push(z_out);
            y_matrix.push(y_out);
        }
        hdm_solves = n_points;
        rom_dim = 0;
        anchor_freqs_hz = anchors;
        max_residual = 0.0;
        symbolic_count = solver.symbolic_count();
        backend_name = solver.backend_name();
        drop(solver);
    }

    let f_rep_hz = freqs_hz[rep_freq_idx];
    let driven_data = DrivenPortData {
        freqs_hz,
        ports: asm.port_names.clone(),
        z0_ohm,
        s_matrix,
        z_matrix,
        y_matrix,
        hdm_solves,
        rom_dim,
        anchor_freqs_hz,
        max_residual,
    };

    Ok(SweepCoreResult { driven_data, c0_matrix, symbolic_count, backend_name, rom_sweep_seconds, f_rep_hz, rep_x_red })
}

fn normalize_and_split_mode(x_full: &[C], n1: usize, d: usize, n_mech: usize) -> (Vec<C>, Vec<C>, bool) {
    let mut big_u = c(0.0, 0.0);
    for &v in &x_full[..n_mech] {
        if v.norm_sqr() > big_u.norm_sqr() {
            big_u = v;
        }
    }
    let big = if big_u.norm() > 0.0 {
        big_u
    } else {
        x_full.iter().copied().fold(c(0.0, 0.0), |m, v| if v.norm_sqr() > m.norm_sqr() { v } else { m })
    };
    let scale = if big.norm() > 0.0 { big.conj() / big.norm_sqr() } else { c(1.0, 0.0) };
    let u_full: Vec<C> = x_full[..n_mech].iter().map(|&v| v * scale).collect();
    let phi_full: Vec<C> = x_full[n_mech..n_mech + n1].iter().map(|&v| v * scale).collect();
    let max_norm = u_full.iter().chain(&phi_full).map(|z| z.norm()).fold(0.0, f64::max).max(1e-300);
    let is_complex = u_full.iter().chain(&phi_full).any(|z| z.im.abs() > 1e-12 * max_norm);
    let _ = d;
    (u_full, phi_full, is_complex)
}

fn quad_real_blocks(m1: &RCsr, u_full: &[C], d: usize, n1: usize) -> f64 {
    let mut tot = 0.0;
    for comp in 0..d {
        let u = &u_full[comp * n1..(comp + 1) * n1];
        for r in 0..m1.nrows() {
            let (cols, vals) = m1.row(r);
            let mut acc = c(0.0, 0.0);
            for (&col, &val) in cols.iter().zip(vals) {
                acc += u[col as usize] * val;
            }
            tot += (u[r].conj() * acc).re;
        }
    }
    tot.max(0.0)
}

fn finish_driven_result(
    asm: &PiezoDrivenAssembly,
    out: SweepCoreResult,
    rep_fields: Vec<nyx_physics::Field>,
    rep_cells: Vec<nyx_physics::CellField>,
    u_full_unscaled: &[C],
    t0: Instant,
    order: usize,
) -> SolveResult {
    let d = asm.constraints.d;
    let n1 = asm.constraints.n1;
    let mut scalars = vec![
        out.driven_data.scalar(),
        ("C_0 (F)".into(), Scalar::Num(out.c0_matrix[0][0])),
        ("hdm_solves".into(), Scalar::Num(out.driven_data.hdm_solves as f64)),
        ("symbolic_count".into(), Scalar::Num(out.symbolic_count as f64)),
        ("rom_dim".into(), Scalar::Num(out.driven_data.rom_dim as f64)),
        ("max_prom_residual".into(), Scalar::Num(out.driven_data.max_residual)),
        ("rom_sweep_seconds".into(), Scalar::Num(out.rom_sweep_seconds)),
        ("f_res_Hz".into(), Scalar::Num(out.f_rep_hz)),
    ];
    if let (Some(ma), Some(mp)) = (&asm.m1_phys_all, &asm.m1_phys_pml) {
        let tot = quad_real_blocks(ma, u_full_unscaled, d, n1);
        let in_layer = quad_real_blocks(mp, u_full_unscaled, d, n1);
        let lf = if tot > 0.0 { (in_layer / tot).clamp(0.0, 1.0) } else { 0.0 };
        scalars.push(("layer fraction 1".into(), Scalar::Num(lf)));
    }
    scalars.push(
        FieldFrame { lateral: 0, normal: 2, propagation: 1, left_handed: false, outside: asm.pml_domains.clone() }
            .scalar(),
    );

    let rep_sol = Solution {
        label: format!("driven ({}): {:.6} MHz", asm.port_names[0], out.f_rep_hz / 1e6),
        fields: rep_fields,
        cell_fields: rep_cells,
        value: Some(out.f_rep_hz),
        eigenvalue: None,
    };

    let mut log = asm.diagnostics.clone();
    log.push(format!(
        "piezo driven sweep [{:.4}..{:.4} MHz, {} points, {} port(s), solver={}, backend={}]: {} HDM solve(s), ROM dim {}, max residual {:.2e}, C_0 = {:.4e} F",
        asm.f_min / 1e6,
        asm.f_max / 1e6,
        out.driven_data.freqs_hz.len(),
        asm.port_names.len(),
        asm.solver,
        out.backend_name,
        out.driven_data.hdm_solves,
        out.driven_data.rom_dim,
        out.driven_data.max_residual,
        out.c0_matrix[0][0],
    ));

    SolveResult {
        physics: "sms".into(),
        study: "driven".into(),
        order: order as u8,
        dim: d as u8,
        solutions: vec![rep_sol],
        scalars,
        log,
        dofs: asm.constraints.n_tot,
        seconds: t0.elapsed().as_secs_f64(),
        diagnostics: None,
        error_budget: None,
    }
}

/// Solves a coupled piezoelectric `Study::Driven` frequency sweep on [`SimplexMesh`].
pub fn solve_piezo_driven(mesh: &SimplexMesh, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let full = assemble_piezo_driven_simplex(mesh, spec)?;
    let asm = PiezoDrivenAssembly::from_full(full, spec)?;
    let out = run_piezo_driven_core(&asm)?;

    let x_full = asm.constraints.red_all.expand_with(&out.rep_x_red, &asm.constraints.terminal_lifts[0]);
    let (u_norm, p_norm, complex_mode) =
        normalize_and_split_mode(&x_full, asm.constraints.n1, asm.constraints.d, asm.constraints.n_mech);
    let (fields, cells) = push_eigen_mode_fields_simplex(
        mesh,
        spec.order,
        asm.constraints.d,
        asm.constraints.n1,
        &u_norm,
        &p_norm,
        complex_mode,
    );
    Ok(finish_driven_result(
        &asm,
        out,
        fields,
        cell_fields_of(&cells),
        &x_full[..asm.constraints.n_mech],
        t0,
        spec.order,
    ))
}

/// Solves a coupled piezoelectric `Study::Driven` frequency sweep natively on [`MixedDiscretisation`].
pub fn solve_piezo_driven_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let space = if spec.order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let full = assemble_piezo_driven_mixed(mixed, &dof_map, spec)?;
    let asm = PiezoDrivenAssembly::from_full(full, spec)?;
    let out = run_piezo_driven_core(&asm)?;

    let x_full = asm.constraints.red_all.expand_with(&out.rep_x_red, &asm.constraints.terminal_lifts[0]);
    let (u_norm, p_norm, complex_mode) =
        normalize_and_split_mode(&x_full, asm.constraints.n1, asm.constraints.d, asm.constraints.n_mech);
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let (fields, cell_fields) = push_eigen_mode_fields_mixed(
        mixed,
        &dof_map,
        &proj,
        spec.order,
        asm.constraints.d,
        asm.constraints.n1,
        &u_norm,
        &p_norm,
        complex_mode,
    );
    Ok(finish_driven_result(&asm, out, fields, cell_fields, &x_full[..asm.constraints.n_mech], t0, spec.order))
}
