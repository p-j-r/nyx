// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Stage 2 & Stage 3 coupled piezoelectric `(u, phi)` solvers for `Study::Static` and
//! `Study::Eigenfrequency`, including block power-of-2 + symmetric Ruiz equilibration,
//! Schur shift-and-invert operators (`SchurShiftInvertOp<T>` and dimensionless
//! `SchurQepShiftInvertOp`), static/modal field recovery, and energy/charge outputs.

use crate::lrb::AMPLIFYING_IM_TOL;
use crate::piezo_asm::{
    assemble_piezo_eigen_mixed, assemble_piezo_eigen_simplex, assemble_piezo_static_mixed,
    assemble_piezo_static_simplex, PiezoEigenAssembly, PiezoStaticAssembly,
};
use crate::piezo_mat::{
    cell_constitutive_blocks, compute_piezo_complex_scaling, compute_piezo_real_scaling, equilibrate_real,
    eval_piezo_point, unscale_solution, PiezoCellMat,
};
use crate::pml::{real_quad_norm, ResolvedPmlDomains};
use crate::problem::{Plane, Problem};
use crate::solid::{lap, von_mises};
use femcore::fem::{geometry, shape, shape_gradients, AmgSystem, RefElement};
use femcore::fem_mixed::{
    lagrange_ref_nodes, lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::{node_average, CellFields};
use femcore::linalg::amg_elastic::AmgOptions;
use femcore::linalg::csr::{c, CCsr, Csr, RCsr, C};
use femcore::linalg::direct::{Definiteness, FactorBackend, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::driver::{cancel_point, catching_cancel, DirectSolver};
use femcore::linalg::eps::{krylov_schur, EpsOptions, EpsReason, Which};
use femcore::linalg::equilibrate::{ruiz_symmetric, scale_symmetric, scale_vec};
use femcore::linalg::faer_glue::FaerScalar;
use femcore::linalg::ksp::{Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
use femcore::linalg::linop::{LinearOperator, Preconditioner};
use femcore::linalg::scalar::{Scalar as LinalgScalar, C64};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::mesh::SimplexMesh;
use femcore::progress;
use femcore::props::LinearSolver;
use femcore::section::{component_names, stress_name};
use nyx_physics::{
    cell_fields_of, CellField, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

// ------------------------------------------------------------------ Stage 2 Schur operators

/// Shift-and-invert Schur operator acting strictly on reduced mechanical displacement DOFs `u` (`dim = n_u`),
/// internally solving the equilibrated monolithic `(u, phi)` saddle-point factorization `(K - sigma M) y = [M_u x; 0]`.
pub struct SchurShiftInvertOp<'a, T: FaerScalar> {
    pub n_u: usize,
    pub n_tot: usize,
    pub m_uu_scaled: &'a Csr<T>,
    pub fac: &'a ShiftedFactorization<T>,
    pub d_scale: &'a [f64],
}

impl<T: FaerScalar> SchurShiftInvertOp<'_, T> {
    fn solve_full_rhs(&self, rhs: &mut [T], nrhs: usize) {
        if nrhs <= 1 || self.n_tot < 10_000 || femcore::par::threads() <= 1 {
            self.fac.solve_block_in_place(rhs, nrhs).expect("SchurShiftInvertOp solver must be factored before apply");
        } else {
            rhs.par_chunks_mut(self.n_tot).for_each(|col| {
                self.fac.solve_block_in_place(col, 1).expect("SchurShiftInvertOp solver must be factored before apply");
            });
        }
    }
}

impl<T: FaerScalar> LinearOperator<T> for SchurShiftInvertOp<'_, T> {
    fn dim(&self) -> usize {
        self.n_u
    }

    fn apply(&self, x: &[T], y: &mut [T]) {
        self.apply_block(x, y, 1);
    }

    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        cancel_point();
        let mut bx = vec![T::ZERO; self.n_u * nrhs];
        for c_idx in 0..nrhs {
            self.m_uu_scaled
                .spmv(&x[c_idx * self.n_u..(c_idx + 1) * self.n_u], &mut bx[c_idx * self.n_u..(c_idx + 1) * self.n_u]);
        }
        self.apply_block_with_bx(x, &bx, y, nrhs);
    }

    fn apply_block_with_bx(&self, _x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        cancel_point();
        let mut rhs = vec![T::ZERO; self.n_tot * nrhs];
        for c_idx in 0..nrhs {
            let src = &bx[c_idx * self.n_u..(c_idx + 1) * self.n_u];
            let dst = &mut rhs[c_idx * self.n_tot..c_idx * self.n_tot + self.n_u];
            for i in 0..self.n_u {
                dst[i] = src[i].scale(self.d_scale[i]);
            }
        }
        self.solve_full_rhs(&mut rhs, nrhs);
        for c_idx in 0..nrhs {
            let src = &rhs[c_idx * self.n_tot..c_idx * self.n_tot + self.n_u];
            let dst = &mut y[c_idx * self.n_u..(c_idx + 1) * self.n_u];
            for i in 0..self.n_u {
                dst[i] = src[i].scale(self.d_scale[i]);
            }
        }
    }
}

/// Companion shift-and-invert operator in **dimensionless** `s = omega / omega_0` on `C^{2 n_u}` for
/// damped/LRB piezoelectric quadratic eigenproblems `(A_0 + s A_1 + s^2 A_2) [u; phi] = 0`:
///
/// `T [u; w] = [w; -P(s_0)^{-1} (A_2 u + P'(s_0) w)]`, eliminating `phi` internally via one equilibrated
/// monolithic saddle-point solve against `[-r_u; 0]`.
pub struct SchurQepShiftInvertOp<'a> {
    pub n_u: usize,
    pub n_tot: usize,
    pub a2_uu: &'a CCsr,
    pub p_prime_uu: &'a CCsr,
    pub fac: &'a ShiftedFactorization<C64>,
    pub d_scale: &'a [f64],
}

impl LinearOperator<C64> for SchurQepShiftInvertOp<'_> {
    fn dim(&self) -> usize {
        2 * self.n_u
    }

    fn apply(&self, x: &[C64], y: &mut [C64]) {
        self.apply_block(x, y, 1);
    }

    fn apply_block(&self, x: &[C64], y: &mut [C64], nrhs: usize) {
        cancel_point();
        let n2 = 2 * self.n_u;
        let mut rhs = vec![C64::ZERO; self.n_tot * nrhs];
        let mut a2u = vec![C64::ZERO; self.n_u];
        let mut ppw = vec![C64::ZERO; self.n_u];
        for c_idx in 0..nrhs {
            let xc = &x[c_idx * n2..(c_idx + 1) * n2];
            let (u, w) = (&xc[..self.n_u], &xc[self.n_u..n2]);
            y[c_idx * n2..c_idx * n2 + self.n_u].copy_from_slice(w);
            self.a2_uu.spmv(u, &mut a2u);
            self.p_prime_uu.spmv(w, &mut ppw);
            let dst = &mut rhs[c_idx * self.n_tot..c_idx * self.n_tot + self.n_u];
            for i in 0..self.n_u {
                dst[i] = -(a2u[i] + ppw[i]) * self.d_scale[i];
            }
        }
        if nrhs <= 1 || self.n_tot < 10_000 || femcore::par::threads() <= 1 {
            self.fac
                .solve_block_in_place(&mut rhs, nrhs)
                .expect("SchurQepShiftInvertOp solver must be factored before apply");
        } else {
            rhs.par_chunks_mut(self.n_tot).for_each(|col| {
                self.fac
                    .solve_block_in_place(col, 1)
                    .expect("SchurQepShiftInvertOp solver must be factored before apply");
            });
        }
        for c_idx in 0..nrhs {
            let src = &rhs[c_idx * self.n_tot..c_idx * self.n_tot + self.n_u];
            let dst = &mut y[c_idx * n2 + self.n_u..(c_idx + 1) * n2];
            for i in 0..self.n_u {
                dst[i] = src[i] * self.d_scale[i];
            }
        }
    }
}

fn estimate_omega2_scale<T: LinalgScalar>(k_red: &Csr<T>, m_uu_red: &Csr<T>, n_u: usize) -> f64 {
    let mut ratios = Vec::with_capacity(n_u);
    let kd = k_red.diagonal();
    let md = m_uu_red.diagonal();
    for i in 0..n_u {
        let m_ii = md[i].modulus();
        let k_ii = kd[i].modulus();
        if m_ii > 1e-300 && k_ii > 0.0 && k_ii.is_finite() && m_ii.is_finite() {
            ratios.push(k_ii / m_ii);
        }
    }
    if ratios.is_empty() {
        return 1.0;
    }
    ratios.sort_by(|a, b| a.total_cmp(b));
    ratios[ratios.len() / 2].max(1.0)
}

fn norm2_vec<T: LinalgScalar>(v: &[T]) -> f64 {
    v.iter().map(|x| x.modulus().powi(2)).sum::<f64>().sqrt()
}

fn recover_phi_real(k_pp: &RCsr, k_pu: &RCsr, u_vecs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n_phi = k_pp.nrows();
    if n_phi == 0 {
        return Ok(vec![Vec::new(); u_vecs.len()]);
    }
    let d_p = ruiz_symmetric(k_pp);
    let k_pp_hat = scale_symmetric(k_pp, &d_p);
    let solver = DirectSolver::new(&k_pp_hat, false).map_err(|e| format!("electrical Schur solve failed: {e}"))?;
    let mut out = Vec::with_capacity(u_vecs.len());
    for u in u_vecs {
        let pu = k_pu.matvec(u);
        let rhs: Vec<f64> = pu.iter().map(|&v| -v).collect();
        let rhs_hat = scale_vec(&rhs, &d_p);
        let p_hat = solver.solve(&rhs_hat);
        out.push(scale_vec(&p_hat, &d_p));
    }
    Ok(out)
}

pub(crate) fn recover_phi_complex(k_pp: &CCsr, k_pu: &CCsr, u_vecs: &[Vec<C64>]) -> Result<Vec<Vec<C64>>, String> {
    let n_phi = k_pp.nrows();
    if n_phi == 0 {
        return Ok(vec![Vec::new(); u_vecs.len()]);
    }
    let d_p = ruiz_symmetric(k_pp);
    let k_pp_hat = scale_symmetric(k_pp, &d_p);
    let solver = DirectSolver::new(&k_pp_hat, false).map_err(|e| format!("electrical Schur solve failed: {e}"))?;
    let mut out = Vec::with_capacity(u_vecs.len());
    for u in u_vecs {
        let pu = k_pu.matvec(u);
        let rhs: Vec<C64> = pu.iter().map(|&v| -v).collect();
        let rhs_hat = scale_vec(&rhs, &d_p);
        let p_hat = solver.solve(&rhs_hat);
        out.push(scale_vec(&p_hat, &d_p));
    }
    Ok(out)
}

pub(crate) struct PiezoComplexEigenInput<'a> {
    pub n_u: usize,
    pub n_tot: usize,
    pub n_modes: usize,
    pub f_search: f64,
    pub is_hermitian: bool,
    pub has_pml: bool,
    pub has_dielectric_loss: bool,
    pub n_pairs: usize,
    pub k_red: &'a CCsr,
    pub m_uu_red: &'a CCsr,
    pub c_uu_red: Option<&'a CCsr>,
    pub k_pp_red: &'a CCsr,
    pub k_pu_red: &'a CCsr,
}

pub(crate) struct PiezoComplexEigenSolved {
    pub freqs: Vec<C>,
    pub reduced_vecs: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
    pub backend: FactorBackend,
}

pub(crate) fn solve_piezo_complex_schur(input: &PiezoComplexEigenInput<'_>) -> Result<PiezoComplexEigenSolved, String> {
    let (n_u, n_tot, nev, f_search) = (input.n_u, input.n_tot, input.n_modes, input.f_search);
    if n_u == 0 {
        return Err("no free mechanical degrees of freedom remain".into());
    }
    if nev < 1 || nev + 1 >= n_u {
        return Err(format!("number of modes ({nev}) must be between 1 and free mechanical dofs - 2 ({n_u})"));
    }
    let k_red = input.k_red;
    let m_uu = input.m_uu_red;
    let k_pp = input.k_pp_red;
    let k_pu = input.k_pu_red;

    let Some(c_uu) = input.c_uu_red else {
        let sigma_phys = if f_search > 0.0 { (2.0 * PI * f_search).powi(2) } else { -1.0 };
        let lam_scale = if f_search > 0.0 {
            (2.0 * PI * f_search).powi(2).max(1.0)
        } else {
            estimate_omega2_scale(k_red, m_uu, n_u)
        };
        let sigma_dimless = if f_search > 0.0 { c(sigma_phys / lam_scale, 0.0) } else { c(-1e-4, 0.0) };
        let m_uu_scaled = m_uu.map(|v| v * c(lam_scale, 0.0));
        let mut m_full_coo = femcore::linalg::csr::CCoo::new(n_tot, n_tot);
        for r in 0..n_u {
            let (cols, vals) = m_uu_scaled.row(r);
            for (&col, &val) in cols.iter().zip(vals) {
                m_full_coo.push(r, col as usize, val);
            }
        }
        let m_full = m_full_coo.to_csr();
        let d_scale =
            compute_piezo_complex_scaling(n_tot, n_u, k_red, &[(k_red, 1.0), (&m_full, sigma_dimless.norm().max(1.0))]);
        let k_hat = scale_symmetric(k_red, &d_scale);
        let m_hat = scale_symmetric(&m_full, &d_scale);

        let mut fac = ShiftedFactorization::new(&k_hat, Some(&m_hat), Definiteness::General)
            .map_err(|e| format!("complex coupled piezo factorisation failed: {e}"))?;
        fac.set_shift(sigma_dimless).map_err(|e| format!("complex coupled piezo shifted factorisation failed: {e}"))?;

        let block = nev.min(8).min(n_u).max(1);
        let eps = EpsOptions {
            which: Which::LargestMagnitude,
            tol: 1e-9,
            block_size: block,
            ..if input.is_hermitian { EpsOptions::hermitian(nev) } else { EpsOptions::general(nev) }
        };
        let op = SchurShiftInvertOp { n_u, n_tot, m_uu_scaled: &m_uu_scaled, fac: &fac, d_scale: &d_scale };
        let weight: Option<&dyn LinearOperator<C64>> = if input.is_hermitian { Some(&m_uu_scaled) } else { None };
        let res = catching_cancel(|| krylov_schur(&op, weight, None, &eps))?
            .map_err(|e| format!("complex coupled piezo eigensolver failed: {e}"))?;

        let th_max = res.values.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
        let th_cut = 1e-15 * th_max;
        let mut ku_buf = vec![C64::ZERO; n_tot];
        let mut mu_buf = vec![C64::ZERO; n_u];
        let mut pairs: Vec<(C64, Vec<C64>)> = Vec::with_capacity(res.values.len());

        let u_raws: Vec<Vec<C64>> = (0..res.values.len()).map(|i| res.vectors.col(i).to_vec()).collect();
        let p_vecs = recover_phi_complex(k_pp, k_pu, &u_raws)?;
        let k_diag_scale = k_red.diagonal()[..n_u].iter().map(|z| z.norm()).fold(0.0_f64, f64::max).max(1.0);

        for (i, &th) in res.values.iter().enumerate() {
            if th.norm() == 0.0 || th.norm() < th_cut {
                continue;
            }
            let mu = sigma_dimless + C64::ONE / th;
            if !mu.is_finite() {
                continue;
            }
            let lam = mu * lam_scale;
            let mut f = lam.sqrt() / (2.0 * PI);
            if f.re < 0.0 {
                f = -f;
            }
            if input.n_pairs > 0 || input.has_pml {
                let mut x_red = Vec::with_capacity(n_tot);
                x_red.extend_from_slice(&u_raws[i]);
                x_red.extend_from_slice(&p_vecs[i]);
                k_red.spmv(&x_red, &mut ku_buf);
                m_uu.spmv(&u_raws[i], &mut mu_buf);
                let r_norm: f64 = (0..n_u).map(|r| (ku_buf[r] - lam * mu_buf[r]).norm_sqr()).sum::<f64>().sqrt();
                let u_norm = norm2_vec(&u_raws[i]);
                let denom =
                    norm2_vec(&ku_buf[..n_u]) + lam.norm() * norm2_vec(&mu_buf) + k_diag_scale * u_norm * 1e-4 + 1e-300;
                if input.n_pairs > 0 && r_norm / denom >= 1e-5 {
                    continue;
                }
            }
            let mut x_red = Vec::with_capacity(n_tot);
            x_red.extend_from_slice(&u_raws[i]);
            x_red.extend_from_slice(&p_vecs[i]);
            pairs.push((f, x_red));
        }
        if f_search > 0.0 {
            pairs.sort_by(|a, b| (a.0.re - f_search).abs().total_cmp(&(b.0.re - f_search).abs()));
        } else {
            pairs.sort_by(|a, b| a.0.re.total_cmp(&b.0.re));
        }
        pairs.truncate(nev);
        pairs.sort_by(|a, b| a.0.re.total_cmp(&b.0.re));
        if pairs.is_empty() {
            return Err("no eigenfrequencies converged near the shift".into());
        }
        return Ok(PiezoComplexEigenSolved {
            freqs: pairs.iter().map(|(f, _)| *f).collect(),
            reduced_vecs: pairs.into_iter().map(|(_, xr)| xr).collect(),
            dim: res.ncv,
            converged: res.reason == EpsReason::Converged,
            backend: fac.stats().backend,
        });
    };

    // Damped / LRB piezoelectric quadratic eigenproblem in dimensionless s = omega / omega_0
    let omega_0 = if f_search > 0.0 { 2.0 * PI * f_search } else { estimate_omega2_scale(k_red, m_uu, n_u).sqrt() };
    let s0 = if f_search > 0.0 { c(1.0, 0.0) } else { c(1e-3, 0.0) };

    let a1_uu = c_uu.map(|v| v * c(0.0, omega_0));
    let a2_uu = m_uu.map(|v| v * c(-omega_0 * omega_0, 0.0));
    let p_prime_uu = Csr::axpby(c(1.0, 0.0), &a1_uu, s0 * 2.0, &a2_uu);

    let pad_uu = |blk: &CCsr| -> CCsr {
        let mut coo = femcore::linalg::csr::CCoo::new(n_tot, n_tot);
        for r in 0..n_u {
            let (cols, vals) = blk.row(r);
            for (&col, &val) in cols.iter().zip(vals) {
                coo.push(r, col as usize, val);
            }
        }
        coo.to_csr()
    };
    let a1_full = pad_uu(&a1_uu);
    let a2_full = pad_uu(&a2_uu);
    let p_s0 = Csr::axpby(c(1.0, 0.0), &Csr::axpby(c(1.0, 0.0), k_red, s0, &a1_full), s0 * s0, &a2_full);

    let d_scale = compute_piezo_complex_scaling(
        n_tot,
        n_u,
        &p_s0,
        &[(k_red, 1.0), (&a1_full, s0.norm().max(1.0)), (&a2_full, (s0 * s0).norm().max(1.0))],
    );
    let p_s0_hat = scale_symmetric(&p_s0, &d_scale);

    let mut fac = ShiftedFactorization::new(&p_s0_hat, None, Definiteness::General)
        .map_err(|e| format!("damped piezo QEP factorisation failed: {e}"))?;
    fac.set_shift(C64::ZERO).map_err(|e| format!("damped piezo QEP shift failed: {e}"))?;

    let nev_q = if f_search > 0.0 { nev } else { (2 * nev).min(2 * n_u - 1).max(nev) };
    let block = nev_q.min(8).min(2 * n_u).max(1);
    let eps = EpsOptions { which: Which::LargestMagnitude, tol: 1e-9, block_size: block, ..EpsOptions::general(nev_q) };
    let op = SchurQepShiftInvertOp { n_u, n_tot, a2_uu: &a2_uu, p_prime_uu: &p_prime_uu, fac: &fac, d_scale: &d_scale };
    let res = catching_cancel(|| krylov_schur(&op, None, None, &eps))?
        .map_err(|e| format!("damped piezo QEP eigensolver failed: {e}"))?;

    let th_max = res.values.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
    let th_cut = 1e-15 * th_max;
    let mut r0 = vec![C64::ZERO; n_u];
    let mut r1 = vec![C64::ZERO; n_u];
    let mut r2 = vec![C64::ZERO; n_u];
    let k_uu = crate::piezo_asm::sub_csr(k_red, 0, n_u, 0, n_u);
    let k_up = crate::piezo_asm::sub_csr(k_red, 0, n_u, n_u, n_tot);

    let mut candidates: Vec<(C64, Vec<C64>)> = Vec::with_capacity(res.values.len());
    for (i, &th) in res.values.iter().enumerate() {
        if th.norm() == 0.0 || th.norm() < th_cut {
            continue;
        }
        let s = s0 + C64::ONE / th;
        if !s.is_finite() {
            continue;
        }
        let u_raw = &res.vectors.col(i)[..n_u];
        let un = norm2_vec(u_raw);
        if un <= 1e-300 {
            continue;
        }
        let mut u: Vec<C64> = u_raw.iter().map(|&z| z / un).collect();
        let p_vec = recover_phi_complex(k_pp, k_pu, &[u.clone()])?.pop().unwrap_or_default();

        k_uu.spmv(&u, &mut r0);
        if k_up.ncols() > 0 {
            let mut ku_p = vec![C64::ZERO; n_u];
            k_up.spmv(&p_vec, &mut ku_p);
            for k in 0..n_u {
                r0[k] += ku_p[k];
            }
        }
        a1_uu.spmv(&u, &mut r1);
        a2_uu.spmv(&u, &mut r2);
        let s2 = s * s;
        let res_num: f64 = (0..n_u).map(|k| (r0[k] + s * r1[k] + s2 * r2[k]).norm_sqr()).sum::<f64>().sqrt();
        let res_den = norm2_vec(&r0) + s.norm() * norm2_vec(&r1) + s2.norm() * norm2_vec(&r2) + 1e-300;
        if res_num / res_den > 1e-4 {
            continue;
        }

        let mut f = s * omega_0 / (2.0 * PI);
        let mut p_fin = p_vec;
        if f.re < -1e-6 * f.im.abs().max(1.0) {
            if !input.has_pml && !input.has_dielectric_loss {
                f = -f.conj();
                u.iter_mut().for_each(|z| *z = z.conj());
                p_fin.iter_mut().for_each(|z| *z = z.conj());
            } else {
                continue;
            }
        }
        if f.re < 0.0 && f.re.abs() < 1e-9 * f.norm().max(1.0) {
            f.re = 0.0;
        }

        let dup = candidates.iter().any(|(fc, _)| (*fc - f).norm() <= 1e-6 * f.norm().max(1.0));
        if !dup {
            let mut x_red = Vec::with_capacity(n_tot);
            x_red.extend(u);
            x_red.extend(p_fin);
            candidates.push((f, x_red));
        }
    }

    if candidates.is_empty() {
        return Err("no coupled piezoelectric eigenfrequencies converged near the shift".into());
    }
    if f_search > 0.0 {
        candidates.sort_by(|a, b| (a.0.re - f_search).abs().total_cmp(&(b.0.re - f_search).abs()));
    } else {
        candidates.sort_by(|a, b| a.0.re.total_cmp(&b.0.re));
    }
    candidates.truncate(nev);
    candidates.sort_by(|a, b| a.0.re.total_cmp(&b.0.re));

    Ok(PiezoComplexEigenSolved {
        freqs: candidates.iter().map(|(f, _)| *f).collect(),
        reduced_vecs: candidates.into_iter().map(|(_, xr)| xr).collect(),
        dim: res.ncv,
        converged: res.reason == EpsReason::Converged,
        backend: fac.stats().backend,
    })
}

pub(crate) struct PiezoModeSolved {
    pub freqs: Vec<C>,
    pub full_vecs: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
    pub backend: FactorBackend,
}

fn solve_piezo_eigen_modes(asm: &PiezoEigenAssembly) -> Result<PiezoModeSolved, String> {
    let n_u = asm.constraints.n_u_red;
    let n_tot = asm.constraints.red_all.n_reduced();
    let nev = asm.n_modes;
    let f_search = asm.f_search;

    if !asm.is_complex {
        let k_red = asm.k_red_real.as_ref().expect("real K_red");
        let m_uu = asm.m_uu_red_real.as_ref().expect("real M_uu_red");
        let sigma_phys = if f_search > 0.0 { (2.0 * PI * f_search).powi(2) } else { -1.0 };
        let lam_scale = if f_search > 0.0 {
            (2.0 * PI * f_search).powi(2).max(1.0)
        } else {
            estimate_omega2_scale(k_red, m_uu, n_u)
        };
        let sigma_dimless = if f_search > 0.0 { sigma_phys / lam_scale } else { -1e-4 };
        let m_uu_scaled = m_uu.map(|v| v * lam_scale);
        let mut m_full_coo = femcore::linalg::csr::Coo::new(n_tot, n_tot);
        for r in 0..n_u {
            let (cols, vals) = m_uu_scaled.row(r);
            for (&col, &val) in cols.iter().zip(vals) {
                m_full_coo.push(r, col as usize, val);
            }
        }
        let m_full = m_full_coo.to_csr();
        let d_scale =
            compute_piezo_real_scaling(n_tot, n_u, k_red, &[(k_red, 1.0), (&m_full, sigma_dimless.abs().max(1.0))]);
        let k_hat = scale_symmetric(k_red, &d_scale);
        let m_hat = scale_symmetric(&m_full, &d_scale);

        let mut fac = ShiftedFactorization::new(&k_hat, Some(&m_hat), Definiteness::General)
            .map_err(|e| format!("coupled piezo factorisation failed: {e}"))?;
        fac.set_shift(sigma_dimless).map_err(|e| format!("coupled piezo shifted factorisation failed: {e}"))?;

        let block = nev.min(8).min(n_u);
        let eps =
            EpsOptions { which: Which::LargestMagnitude, tol: 1e-10, block_size: block, ..EpsOptions::hermitian(nev) };
        let op = SchurShiftInvertOp { n_u, n_tot, m_uu_scaled: &m_uu_scaled, fac: &fac, d_scale: &d_scale };
        let weight: Option<&dyn LinearOperator<f64>> = Some(&m_uu_scaled);
        let res = catching_cancel(|| krylov_schur(&op, weight, None, &eps))?
            .map_err(|e| format!("coupled piezo eigensolver failed: {e}"))?;

        let th_max = res.values.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
        let th_cut = 1e-15 * th_max;
        let mut modes: Vec<(f64, Vec<f64>)> = Vec::with_capacity(res.values.len());
        for (i, &th) in res.values.iter().enumerate() {
            if th.norm() == 0.0 || th.norm() < th_cut {
                continue;
            }
            let mu = sigma_dimless + 1.0 / th.re;
            if !mu.is_finite() {
                continue;
            }
            let lam = mu * lam_scale;
            let f = lam.max(0.0).sqrt() / (2.0 * PI);
            modes.push((f, res.vectors.col(i).to_vec()));
        }
        if f_search > 0.0 {
            modes.sort_by(|a, b| (a.0 - f_search).abs().total_cmp(&(b.0 - f_search).abs()));
        } else {
            modes.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
        modes.truncate(nev);
        modes.sort_by(|a, b| a.0.total_cmp(&b.0));
        if modes.is_empty() {
            return Err("no coupled piezoelectric eigenfrequencies converged near the shift".into());
        }
        let u_vecs: Vec<Vec<f64>> = modes.iter().map(|(_, u)| u.clone()).collect();
        let p_vecs =
            recover_phi_real(asm.k_pp_red_real.as_ref().unwrap(), asm.k_pu_red_real.as_ref().unwrap(), &u_vecs)?;
        let mut full_vecs = Vec::with_capacity(modes.len());
        for (u, p) in u_vecs.into_iter().zip(p_vecs) {
            let mut xr = Vec::with_capacity(n_tot);
            xr.extend(u);
            xr.extend(p);
            let full = asm
                .constraints
                .red_all
                .expand_with(&xr.into_iter().map(|v| c(v, 0.0)).collect::<Vec<_>>(), &asm.constraints.g_all);
            full_vecs.push(full);
        }
        return Ok(PiezoModeSolved {
            freqs: modes.into_iter().map(|(f, _)| c(f, 0.0)).collect(),
            full_vecs,
            dim: res.ncv,
            converged: res.reason == EpsReason::Converged,
            backend: fac.stats().backend,
        });
    }

    let solved = solve_piezo_complex_schur(&PiezoComplexEigenInput {
        n_u,
        n_tot,
        n_modes: nev,
        f_search,
        is_hermitian: asm.is_hermitian,
        has_pml: asm.has_pml,
        has_dielectric_loss: asm.mats.iter().any(|m| m.eps_im.iter().flatten().any(|&v| v != 0.0)),
        n_pairs: asm.n_pairs,
        k_red: asm.k_red_c.as_ref().expect("complex K_red"),
        m_uu_red: asm.m_uu_red_c.as_ref().expect("complex M_uu_red"),
        c_uu_red: if asm.has_damping { asm.c_uu_red_c.as_ref() } else { None },
        k_pp_red: asm.k_pp_red_c.as_ref().unwrap(),
        k_pu_red: asm.k_pu_red_c.as_ref().unwrap(),
    })?;
    let full_vecs =
        solved.reduced_vecs.iter().map(|xr| asm.constraints.red_all.expand_with(xr, &asm.constraints.g_all)).collect();
    Ok(PiezoModeSolved {
        freqs: solved.freqs,
        full_vecs,
        dim: solved.dim,
        converged: solved.converged,
        backend: solved.backend,
    })
}

// ------------------------------------------------------------------ Stage 3 physical energy & field evaluation

fn cell_energies_simplex(
    mesh: &SimplexMesh,
    mats: &[PiezoCellMat],
    order: usize,
    plane: Plane,
    thickness: f64,
    u_full: &[C],
    phi_full: &[C],
    resolved_pml: &ResolvedPmlDomains,
    n1: usize,
) -> (f64, f64) {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let thk = if d == 2 { thickness } else { 1.0 };
    let mut grad = vec![0.0; nloc * d];
    let mut scratch = Vec::new();
    let (mut ue, mut ud) = (0.0, 0.0);

    for ci in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[ci];
        if resolved_pml.contains(dom) {
            continue;
        }
        let mat = &mats[ci];
        if !mat.mech_active && !mat.elec_active {
            continue;
        }
        mesh.cell_dofs(ci, order, &mut scratch);
        let vol_thk = geo.vol[ci] * thk;
        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            let w_vol = re.w[q] * vol_thk;
            let mut g_u = [[c(0.0, 0.0); 3]; 3];
            let mut g_p = [c(0.0, 0.0); 3];
            for (a, &gd) in scratch[..nloc].iter().enumerate() {
                let pa = phi_full[gd];
                for k in 0..d {
                    let gka = grad[a * d + k];
                    g_p[k] += pa * gka;
                    for comp in 0..d {
                        g_u[comp][k] += u_full[comp * n1 + gd] * gka;
                    }
                }
            }
            let (ce, de) = eval_complex_point_energies(mat, d, plane, &g_u, &g_p);
            ue += w_vol * ce;
            ud += w_vol * de;
        }
    }
    (ue, ud)
}

fn cell_energies_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    mats: &[PiezoCellMat],
    plane: Plane,
    thickness: f64,
    u_full: &[C],
    phi_full: &[C],
    resolved_pml: &ResolvedPmlDomains,
) -> (f64, f64) {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let thk = if d == 2 { thickness } else { 1.0 };
    let (mut vals, mut dxi, mut sn, mut sdn) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut ue, mut ud) = (0.0, 0.0);

    for ci in 0..mixed.n_cells() {
        let dom = mixed.cell_labels[ci];
        if resolved_pml.contains(dom) {
            continue;
        }
        let mat = &mats[ci];
        if !mat.mech_active && !mat.elec_active {
            continue;
        }
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let deg = if iso.curved { 2 * order + 4 } else { 2 * order + 2 };
        let rule = quad_rule_mixed(ct, deg);
        let cdofs = dof_map.cell_dofs(ci);

        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let w_vol = wq * jg.det.abs() * thk;
            let mut g_u = [[c(0.0, 0.0); 3]; 3];
            let mut g_p = [c(0.0, 0.0); 3];
            for (a, &gd) in cdofs.iter().enumerate() {
                let gp = jg.map_grad(dxi[a]);
                let pa = phi_full[gd];
                for k in 0..d {
                    let gka = gp[k];
                    g_p[k] += pa * gka;
                    for comp in 0..d {
                        g_u[comp][k] += u_full[comp * n1 + gd] * gka;
                    }
                }
            }
            let (ce, de) = eval_complex_point_energies(mat, d, plane, &g_u, &g_p);
            ue += w_vol * ce;
            ud += w_vol * de;
        }
    }
    (ue, ud)
}

#[inline]
fn eval_complex_point_energies(
    mat: &PiezoCellMat,
    d: usize,
    plane: Plane,
    g_u: &[[C; 3]; 3],
    g_p: &[C; 3],
) -> (f64, f64) {
    let mut ue = 0.0;
    let mut ud = 0.0;
    if d == 2 {
        let s = [g_u[0][0], g_u[1][1], g_u[0][1] + g_u[1][0]];
        let blk = cell_constitutive_blocks(mat, 2, plane);
        if mat.mech_active {
            for a in 0..3 {
                let mut cs = c(0.0, 0.0);
                for b in 0..3 {
                    cs += s[b] * blk.d_uu[a * 3 + b];
                }
                ue += 0.5 * (s[a].conj() * cs).re;
            }
        }
        if mat.elec_active {
            for p in 0..2 {
                let mut ee = c(0.0, 0.0);
                for q in 0..2 {
                    ee += g_p[q] * blk.eps_ps_re[p * 2 + q];
                }
                ud += 0.5 * (g_p[p].conj() * ee).re;
            }
        }
    } else {
        let s = [g_u[0][0], g_u[1][1], g_u[2][2], g_u[0][1] + g_u[1][0], g_u[1][2] + g_u[2][1], g_u[0][2] + g_u[2][0]];
        if mat.mech_active {
            for a in 0..6 {
                let mut cs = c(0.0, 0.0);
                for b in 0..6 {
                    cs += s[b] * mat.c6[a][b];
                }
                ue += 0.5 * (s[a].conj() * cs).re;
            }
        }
        if mat.elec_active {
            for p in 0..3 {
                let mut ee = c(0.0, 0.0);
                for q in 0..3 {
                    ee += g_p[q] * mat.eps_re[p][q];
                }
                ud += 0.5 * (g_p[p].conj() * ee).re;
            }
        }
    }
    (ue, ud)
}

fn push_static_charge_and_energy_scalars(
    asm: &PiezoStaticAssembly,
    x_full: &[f64],
    scalars: &mut Vec<(String, Scalar)>,
    ue: f64,
    ud: f64,
) {
    let utot = ue + ud;
    scalars.push(("elastic energy (J)".into(), Scalar::Num(ue)));
    scalars.push(("dielectric energy (J)".into(), Scalar::Num(ud)));
    if utot > 0.0 {
        scalars.push(("mechanical energy share".into(), Scalar::Num(ue / utot)));
        scalars.push(("electrical energy share".into(), Scalar::Num(ud / utot)));
    }

    // Reaction charge per electrode: `-(K_full * x_full - f_th)[n_mech + dof]`
    let kx = asm.k_full.matvec(x_full);
    let n_mech = asm.constraints.n_mech;
    let mut driven_pairs: Vec<(f64, f64)> = Vec::new();
    let mut has_ground = false;

    for (k, el) in asm.constraints.electrodes.iter().enumerate() {
        let q: f64 = el.dofs.iter().map(|&dof| -(kx[n_mech + dof] - asm.f_th_only[n_mech + dof])).sum();
        scalars.push((format!("charge {} (C)", k + 1), Scalar::Num(q)));
        if let Some(root) = el.root_dof {
            scalars.push((format!("V {} (V)", k + 1), Scalar::Num(x_full[n_mech + root])));
        }
        if let Some(v) = el.fixed_voltage {
            if v == 0.0 {
                has_ground = true;
            } else {
                driven_pairs.push((v, q));
            }
        }
    }
    if driven_pairs.len() == 1 && has_ground && driven_pairs[0].0.abs() > 0.0 {
        scalars.push(("capacitance (F)".into(), Scalar::Num(driven_pairs[0].1 / driven_pairs[0].0)));
    }
}

fn static_fields_simplex(
    mesh: &SimplexMesh,
    order: usize,
    asm: &PiezoStaticAssembly,
    x_full: &[f64],
    plane: Plane,
) -> (Vec<Field>, Vec<CellField>, f64, f64, f64) {
    let d = asm.constraints.d;
    let n1 = asm.constraints.n1;
    let n_mech = asm.constraints.n_mech;
    let u = &x_full[..n_mech];
    let phi = &x_full[n_mech..];
    let names = component_names(d);

    let u_comps: Vec<Vec<f64>> = (0..d).map(|k| u[k * n1..(k + 1) * n1].to_vec()).collect();
    let mut fields = vec![Field::vector("u", "m", names, &u_comps), Field::scalar("phi", "V", phi.to_vec())];
    let mut cells = CellFields::new(mesh, order);
    // `ux`, `uz`, [`uy`], `phi` first in cell_fields
    for (k, ax) in names.iter().enumerate() {
        cells.push_real(
            format!("u{ax}"),
            "m",
            femcore::fem::cell_values_at_nodes(mesh, &u[k * n1..(k + 1) * n1], order),
        );
    }
    cells.push_real("phi", "V", femcore::fem::cell_values_at_nodes(mesh, phi, order));

    let nc = mesh.n_cells();
    let nloc = mesh.nloc(order);
    let geo = geometry(mesh);
    let (mut lam_nodes, _) = match d {
        2 => (
            if order == 1 {
                vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]]
            } else {
                vec![
                    vec![1.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 1.0],
                    vec![0.5, 0.5, 0.0],
                    vec![0.5, 0.0, 0.5],
                    vec![0.0, 0.5, 0.5],
                ]
            },
            0,
        ),
        _ => (
            if order == 1 {
                vec![
                    vec![1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0],
                ]
            } else {
                vec![
                    vec![1.0, 0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0, 1.0],
                    vec![0.5, 0.5, 0.0, 0.0],
                    vec![0.5, 0.0, 0.5, 0.0],
                    vec![0.5, 0.0, 0.0, 0.5],
                    vec![0.0, 0.5, 0.5, 0.0],
                    vec![0.0, 0.5, 0.0, 0.5],
                    vec![0.0, 0.0, 0.5, 0.5],
                ]
            },
            0,
        ),
    };
    let _ = &mut lam_nodes;
    let shape_grads: Vec<Vec<f64>> = lam_nodes.iter().map(|l| shape(l, order).1).collect();

    let mut s_all = vec![[0.0; 6]; nc * nloc];
    let mut e_cell = vec![Vec::with_capacity(nc * nloc); d];
    let mut d_cell = vec![Vec::with_capacity(nc * nloc); d];
    let mut dofs = Vec::new();
    let mut max_e = 0.0_f64;

    for ci in 0..nc {
        mesh.cell_dofs(ci, order, &mut dofs);
        let mat = &asm.mats[ci];
        let base = ci * (d + 1) * d;
        for a_idx in 0..nloc {
            let dn = &shape_grads[a_idx];
            let mut g_u = [[0.0; 3]; 3];
            let mut g_p = [0.0; 3];
            for node_i in 0..nloc {
                let mut gvec = [0.0; 3];
                for k in 0..d {
                    let mut ga = 0.0;
                    for v_idx in 0..=d {
                        ga += dn[node_i * (d + 1) + v_idx] * geo.g[base + v_idx * d + k];
                    }
                    gvec[k] = ga;
                }
                let gd = dofs[node_i];
                let pv = phi[gd];
                for k in 0..d {
                    g_p[k] += pv * gvec[k];
                    for comp in 0..d {
                        g_u[comp][k] += u[comp * n1 + gd] * gvec[k];
                    }
                }
            }
            let th_pt = asm.thermal.as_ref().map(|th| (&th.alpha[ci], th.t[dofs[a_idx]] - th.tref));
            let ev = eval_piezo_point(mat, d, plane, &g_u, &g_p, th_pt);
            s_all[ci * nloc + a_idx] = ev.stress;
            let emag = (0..d).map(|k| ev.e_field[k] * ev.e_field[k]).sum::<f64>().sqrt();
            max_e = max_e.max(emag);
            for k in 0..d {
                e_cell[k].push(ev.e_field[k]);
                d_cell[k].push(ev.d_field[k]);
            }
        }
    }

    let mut e_nodal = Vec::with_capacity(d);
    let mut d_nodal = Vec::with_capacity(d);
    for (k, ax) in names.iter().enumerate() {
        e_nodal.push(node_average(mesh, &e_cell[k], 1, order));
        cells.push_real(format!("E{ax}"), "V/m", e_cell[k].clone());
    }
    for (k, ax) in names.iter().enumerate() {
        d_nodal.push(node_average(mesh, &d_cell[k], 1, order));
        cells.push_real(format!("D{ax}"), "C/m^2", d_cell[k].clone());
    }

    let vm: Vec<f64> = s_all.iter().map(von_mises).collect();
    let vm_nodal = node_average(mesh, &vm, 1, order);
    let vmax = vm_nodal.iter().cloned().fold(0.0, f64::max);
    fields.push(Field::scalar("von Mises", "Pa", vm_nodal));
    cells.push_real("von Mises", "Pa", vm);

    let idxs: &[usize] = if d == 2 { &[0, 1, 3, 2] } else { &[0, 1, 2, 3, 4, 5] };
    for &idx in idxs {
        let v: Vec<f64> = s_all.iter().map(|x| x[idx]).collect();
        let sname = stress_name(d, idx);
        fields.push(Field::scalar(sname, "Pa", node_average(mesh, &v, 1, order)));
        cells.push_real(sname, "Pa", v);
    }

    fields.push(Field::vector("E", "V/m", names, &e_nodal));
    fields.push(Field::vector("D", "C/m^2", names, &d_nodal));

    let umax = (0..n1).map(|i| (0..d).map(|k| u[k * n1 + i].powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max);
    (fields, cell_fields_of(&cells), umax, vmax, max_e)
}

#[allow(clippy::too_many_arguments)]
fn static_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    order: usize,
    asm: &PiezoStaticAssembly,
    x_full: &[f64],
    plane: Plane,
) -> (Vec<Field>, Vec<CellField>, f64, f64, f64) {
    let d = asm.constraints.d;
    let n1 = asm.constraints.n1;
    let n_mech = asm.constraints.n_mech;
    let u = &x_full[..n_mech];
    let phi = &x_full[n_mech..];
    let space = dof_map.space;
    let nc = mixed.n_cells();
    let names = component_names(d);

    let mut fields = vec![
        proj.project_vector_field(mixed, dof_map, "u", "m", names, u, order),
        proj.project_scalar_field(mixed, dof_map, "phi", "V", phi, order),
    ];

    let mut m_offsets = Vec::with_capacity(nc + 1);
    m_offsets.push(0);
    for mc in 0..nc {
        m_offsets.push(m_offsets[mc] + dof_map.cell_dofs(mc).len());
    }
    let total_cell_nodes = *m_offsets.last().unwrap();
    let mut cell_fields = Vec::new();

    for (k, ax) in names.iter().enumerate() {
        let mut uk_mixed = Vec::with_capacity(total_cell_nodes);
        for mc in 0..nc {
            for &gd in dof_map.cell_dofs(mc) {
                uk_mixed.push(u[k * n1 + gd]);
            }
        }
        cell_fields.push(proj.project_cell_field(
            mixed,
            space,
            format!("u{ax}"),
            "m",
            &uk_mixed,
            &[],
            &m_offsets,
            order,
        ));
    }
    let mut phi_mixed = Vec::with_capacity(total_cell_nodes);
    for mc in 0..nc {
        for &gd in dof_map.cell_dofs(mc) {
            phi_mixed.push(phi[gd]);
        }
    }
    cell_fields.push(proj.project_cell_field(mixed, space, "phi", "V", &phi_mixed, &[], &m_offsets, order));

    let mut s_all = Vec::with_capacity(total_cell_nodes);
    let mut e_cell = vec![Vec::with_capacity(total_cell_nodes); d];
    let mut d_cell = vec![Vec::with_capacity(total_cell_nodes); d];
    let (mut vals, mut dxi, mut sn, mut sdn) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut max_e = 0.0_f64;

    for mc in 0..nc {
        let dom = mixed.cell_labels[mc];
        let mat = &asm.mats[mc];
        let cdofs = dof_map.cell_dofs(mc);
        let ct = mixed.cell_type(mc);
        let iso = IsoparametricCell::from_mesh(mixed, mc);
        let in_ie = asm.ie_domains.contains(dom);
        let ref_nodes = lagrange_ref_nodes(ct, space);

        for (a, xi) in ref_nodes.iter().enumerate() {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let mut g_u = [[0.0; 3]; 3];
            let mut g_p = [0.0; 3];
            for (b_idx, &gd) in cdofs.iter().enumerate() {
                let gp = jg.map_grad(dxi[b_idx]);
                for k in 0..d {
                    g_p[k] += phi[gd] * gp[k];
                    for comp in 0..d {
                        g_u[comp][k] += u[comp * n1 + gd] * gp[k];
                    }
                }
            }
            if in_ie {
                let ie_j = asm.ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                let mut g_u_ie = [[0.0; 3]; 3];
                let mut g_p_ie = [0.0; 3];
                for m in 0..d {
                    for k in 0..d {
                        g_p_ie[m] += g_p[k] * ie_j.j_inv[k][m];
                        for comp in 0..d {
                            g_u_ie[comp][m] += g_u[comp][k] * ie_j.j_inv[k][m];
                        }
                    }
                }
                g_u = g_u_ie;
                g_p = g_p_ie;
            }
            let th_pt = asm.thermal.as_ref().map(|th| (&th.alpha[mc], th.t[cdofs[a]] - th.tref));
            let ev = eval_piezo_point(mat, d, plane, &g_u, &g_p, th_pt);
            s_all.push(ev.stress);
            let emag = (0..d).map(|k| ev.e_field[k] * ev.e_field[k]).sum::<f64>().sqrt();
            max_e = max_e.max(emag);
            for k in 0..d {
                e_cell[k].push(ev.e_field[k]);
                d_cell[k].push(ev.d_field[k]);
            }
        }
    }

    let mut e_nodal = Vec::with_capacity(d);
    let mut d_nodal = Vec::with_capacity(d);
    for (k, ax) in names.iter().enumerate() {
        let cf = proj.project_cell_field(mixed, space, format!("E{ax}"), "V/m", &e_cell[k], &[], &m_offsets, order);
        e_nodal.push(node_average(simplex, &cf.re, 1, order));
        cell_fields.push(cf);
    }
    for (k, ax) in names.iter().enumerate() {
        let cf = proj.project_cell_field(mixed, space, format!("D{ax}"), "C/m^2", &d_cell[k], &[], &m_offsets, order);
        d_nodal.push(node_average(simplex, &cf.re, 1, order));
        cell_fields.push(cf);
    }

    let vm_mixed: Vec<f64> = s_all.iter().map(von_mises).collect();
    let vm_cf = proj.project_cell_field(mixed, space, "von Mises", "Pa", &vm_mixed, &[], &m_offsets, order);
    let vm_nodal = node_average(simplex, &vm_cf.re, 1, order);
    let vmax = vm_nodal.iter().cloned().fold(0.0, f64::max);
    fields.push(Field::scalar("von Mises", "Pa", vm_nodal));
    cell_fields.push(vm_cf);

    let idxs: &[usize] = if d == 2 { &[0, 1, 3, 2] } else { &[0, 1, 2, 3, 4, 5] };
    for &idx in idxs {
        let sname = stress_name(d, idx);
        let comp_mixed: Vec<f64> = s_all.iter().map(|x| x[idx]).collect();
        let s_cf = proj.project_cell_field(mixed, space, sname, "Pa", &comp_mixed, &[], &m_offsets, order);
        fields.push(Field::scalar(sname, "Pa", node_average(simplex, &s_cf.re, 1, order)));
        cell_fields.push(s_cf);
    }

    fields.push(Field::vector("E", "V/m", names, &e_nodal));
    fields.push(Field::vector("D", "C/m^2", names, &d_nodal));

    let umax = (0..n1).map(|i| (0..d).map(|k| u[k * n1 + i].powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max);
    (fields, cell_fields, umax, vmax, max_e)
}

pub(crate) fn push_eigen_mode_fields_simplex(
    mesh: &SimplexMesh,
    order: usize,
    d: usize,
    n1: usize,
    u_full: &[C],
    phi_full: &[C],
    complex_mode: bool,
) -> (Vec<Field>, CellFields) {
    let names = component_names(d);
    let mut cells = CellFields::new(mesh, order);
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();

    let mut fields = Vec::new();
    if !complex_mode {
        let u_r: Vec<f64> = u_full.iter().map(|z| z.re).collect();
        let p_r: Vec<f64> = phi_full.iter().map(|z| z.re).collect();
        let u_comps: Vec<Vec<f64>> = (0..d).map(|k| u_r[k * n1..(k + 1) * n1].to_vec()).collect();
        fields.push(Field::vector("u", "a.u.", names, &u_comps));
        fields.push(Field::scalar("phi", "a.u.", p_r.clone()));
        for (k, ax) in names.iter().enumerate() {
            cells.push_real(
                format!("u{ax}"),
                "a.u.",
                femcore::fem::cell_values_at_nodes(mesh, &u_r[k * n1..(k + 1) * n1], order),
            );
        }
        cells.push_real("phi", "a.u.", femcore::fem::cell_values_at_nodes(mesh, &p_r, order));
    } else {
        let mag: Vec<f64> =
            (0..n1).map(|a| (0..d).map(|k| u_full[k * n1 + a].norm_sqr()).sum::<f64>().sqrt()).collect();
        fields.push(Field::scalar("|u|", "a.u.", mag));
        for (comp, ax) in names.iter().enumerate() {
            let v: Vec<f64> = (0..n1).map(|a| u_full[comp * n1 + a].re).collect();
            fields.push(Field::scalar(format!("Re u{ax}"), "a.u.", v));
            let (mut re_c, mut im_c) = (Vec::with_capacity(nc * nloc), Vec::with_capacity(nc * nloc));
            let mut dofs = Vec::new();
            for cell in 0..nc {
                mesh.cell_dofs(cell, order, &mut dofs);
                for &dof in dofs.iter().take(nloc) {
                    let z = u_full[comp * n1 + dof];
                    re_c.push(z.re);
                    im_c.push(z.im);
                }
            }
            cells.push_complex(format!("u{ax}"), "a.u.", re_c, im_c);
        }
        let phi_mag: Vec<f64> = phi_full.iter().map(|z| z.norm()).collect();
        let phi_re: Vec<f64> = phi_full.iter().map(|z| z.re).collect();
        fields.push(Field::scalar("|phi|", "a.u.", phi_mag));
        fields.push(Field::scalar("Re phi", "a.u.", phi_re));
        let (mut re_p, mut im_p) = (Vec::with_capacity(nc * nloc), Vec::with_capacity(nc * nloc));
        let mut dofs = Vec::new();
        for cell in 0..nc {
            mesh.cell_dofs(cell, order, &mut dofs);
            for &dof in dofs.iter().take(nloc) {
                let z = phi_full[dof];
                re_p.push(z.re);
                im_p.push(z.im);
            }
        }
        cells.push_complex("phi", "a.u.", re_p, im_p);
    }
    (fields, cells)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn push_eigen_mode_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    proj: &MixedToSimplexProjector,
    order: usize,
    d: usize,
    n1: usize,
    u_full: &[C],
    phi_full: &[C],
    complex_mode: bool,
) -> (Vec<Field>, Vec<CellField>) {
    let names = component_names(d);
    let space = dof_map.space;
    let nc = mixed.n_cells();
    let mut m_offsets = Vec::with_capacity(nc + 1);
    m_offsets.push(0);
    for mc in 0..nc {
        m_offsets.push(m_offsets[mc] + dof_map.cell_dofs(mc).len());
    }
    let total_nodes = *m_offsets.last().unwrap();
    let mut fields = Vec::new();
    let mut cell_fields = Vec::with_capacity(d + 1);

    if !complex_mode {
        let u_r: Vec<f64> = u_full.iter().map(|z| z.re).collect();
        let p_r: Vec<f64> = phi_full.iter().map(|z| z.re).collect();
        fields.push(proj.project_vector_field(mixed, dof_map, "u", "a.u.", names, &u_r, order));
        fields.push(proj.project_scalar_field(mixed, dof_map, "phi", "a.u.", &p_r, order));
        for (comp, ax) in names.iter().enumerate() {
            let mut c_u = Vec::with_capacity(total_nodes);
            for mc in 0..nc {
                for &gd in dof_map.cell_dofs(mc) {
                    c_u.push(u_r[comp * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                space,
                format!("u{ax}"),
                "a.u.",
                &c_u,
                &[],
                &m_offsets,
                order,
            ));
        }
        let mut c_p = Vec::with_capacity(total_nodes);
        for mc in 0..nc {
            for &gd in dof_map.cell_dofs(mc) {
                c_p.push(p_r[gd]);
            }
        }
        cell_fields.push(proj.project_cell_field(mixed, space, "phi", "a.u.", &c_p, &[], &m_offsets, order));
    } else {
        let re_u_blocks: Vec<f64> = u_full.iter().map(|z| z.re).collect();
        let im_u_blocks: Vec<f64> = u_full.iter().map(|z| z.im).collect();
        let re_u_proj = proj.project_lagrange_blocks(mixed, dof_map, &re_u_blocks, d, order);
        let im_u_proj = proj.project_lagrange_blocks(mixed, dof_map, &im_u_blocks, d, order);
        let n_out = re_u_proj.len() / d;
        let mag: Vec<f64> = (0..n_out)
            .map(|a| {
                (0..d).map(|k| re_u_proj[k * n_out + a].powi(2) + im_u_proj[k * n_out + a].powi(2)).sum::<f64>().sqrt()
            })
            .collect();
        fields.push(Field::scalar("|u|", "a.u.", mag));
        for (comp, ax) in names.iter().enumerate() {
            fields.push(Field::scalar(
                format!("Re u{ax}"),
                "a.u.",
                re_u_proj[comp * n_out..(comp + 1) * n_out].to_vec(),
            ));
            let (mut re_c, mut im_c) = (Vec::with_capacity(total_nodes), Vec::with_capacity(total_nodes));
            for mc in 0..nc {
                for &gd in dof_map.cell_dofs(mc) {
                    re_c.push(re_u_blocks[comp * n1 + gd]);
                    im_c.push(im_u_blocks[comp * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                space,
                format!("u{ax}"),
                "a.u.",
                &re_c,
                &im_c,
                &m_offsets,
                order,
            ));
        }

        let re_p_blocks: Vec<f64> = phi_full.iter().map(|z| z.re).collect();
        let im_p_blocks: Vec<f64> = phi_full.iter().map(|z| z.im).collect();
        let re_p_proj = proj.project_lagrange_blocks(mixed, dof_map, &re_p_blocks, 1, order);
        let im_p_proj = proj.project_lagrange_blocks(mixed, dof_map, &im_p_blocks, 1, order);
        let phi_mag: Vec<f64> = re_p_proj.iter().zip(&im_p_proj).map(|(&r, &i)| r.hypot(i)).collect();
        fields.push(Field::scalar("|phi|", "a.u.", phi_mag));
        fields.push(Field::scalar("Re phi", "a.u.", re_p_proj));
        let (mut re_pc, mut im_pc) = (Vec::with_capacity(total_nodes), Vec::with_capacity(total_nodes));
        for mc in 0..nc {
            for &gd in dof_map.cell_dofs(mc) {
                re_pc.push(re_p_blocks[gd]);
                im_pc.push(im_p_blocks[gd]);
            }
        }
        cell_fields.push(proj.project_cell_field(mixed, space, "phi", "a.u.", &re_pc, &im_pc, &m_offsets, order));
    }
    (fields, cell_fields)
}

/// Builds SPD `M_phi = K_pp + diag_add` from negated `-K_pp` (`k_pp_neg`).
fn negate_add_diag_csr(k_pp_neg: &RCsr, diag_add: &[f64]) -> RCsr {
    let n = k_pp_neg.nrows();
    let mut indptr = Vec::with_capacity(n + 1);
    indptr.push(0);
    let mut indices = Vec::with_capacity(k_pp_neg.nnz() + n);
    let mut values = Vec::with_capacity(k_pp_neg.nnz() + n);
    for (r, &add) in diag_add.iter().enumerate().take(n) {
        let (cols, vals) = k_pp_neg.row(r);
        let mut inserted = false;
        for (&col_u32, &val) in cols.iter().zip(vals) {
            let col = col_u32 as usize;
            if !inserted && col == r {
                indices.push(col);
                values.push(-val + add);
                inserted = true;
            } else {
                if !inserted && col > r && add != 0.0 {
                    indices.push(r);
                    values.push(add);
                    inserted = true;
                }
                indices.push(col);
                values.push(-val);
            }
        }
        if !inserted && add != 0.0 {
            indices.push(r);
            values.push(add);
        }
        indptr.push(indices.len());
    }
    RCsr::from_raw(n, k_pp_neg.ncols(), indptr, indices, values)
}

/// Builds 1-hop SPD Schur approximation `M_phi = K_pp + diag(sum_k K_{pu,ik}^2 / K_{uu,kk})`
/// along with reduced `(K_{u phi}, M_phi, diag_add)`.
fn build_m_phi_with_coupling(k_red: &RCsr, n_u_red: usize, n_phi_red: usize) -> (RCsr, RCsr, Vec<f64>) {
    let n_tot_red = n_u_red + n_phi_red;
    let k_up_red = crate::piezo_asm::sub_csr(k_red, 0, n_u_red, n_u_red, n_tot_red);
    let k_pu_red = crate::piezo_asm::sub_csr(k_red, n_u_red, n_tot_red, 0, n_u_red);
    let k_pp_neg = crate::piezo_asm::sub_csr(k_red, n_u_red, n_tot_red, n_u_red, n_tot_red);

    let mut diag_uu = vec![0.0; n_u_red];
    for (r, d_slot) in diag_uu.iter_mut().enumerate() {
        let (cols, vals) = k_red.row(r);
        for (&col_u32, &val) in cols.iter().zip(vals) {
            let col = col_u32 as usize;
            if col == r {
                *d_slot = val;
                break;
            }
            if col > r {
                break;
            }
        }
    }

    let mut diag_add = vec![0.0; n_phi_red];
    for (p, add_slot) in diag_add.iter_mut().enumerate() {
        let (cols, vals) = k_pu_red.row(p);
        let mut sum = 0.0;
        for (&col_u32, &val) in cols.iter().zip(vals) {
            let d_kk = diag_uu[col_u32 as usize];
            if d_kk > 0.0 {
                sum += (val * val) / d_kk;
            }
        }
        *add_slot = sum;
    }

    let m_phi_red = negate_add_diag_csr(&k_pp_neg, &diag_add);
    (k_up_red, m_phi_red, diag_add)
}

enum PhiPc {
    Direct(Box<DirectSolver<f64>>),
    Amg(Box<AmgSystem>),
}

impl PhiPc {
    fn solve_neg(&self, r_p: &[f64]) -> Vec<f64> {
        let neg_rp: Vec<f64> = r_p.iter().map(|&v| -v).collect();
        match self {
            PhiPc::Direct(ds) => ds.solve(&neg_rp),
            PhiPc::Amg(sys) => sys.precondition_free(&neg_rp),
        }
    }
}

/// Block upper-triangular Schur preconditioner for equilibrated static piezoelectricity
/// (`K = [[K_uu, K_{u phi}], [K_{phi u}, -K_{phi phi}]]`, `P = [[K_uu, K_{u phi}], [0, -M_phi]]`).
struct PiezoSchurPc<'a> {
    n_u_red: usize,
    n_phi_red: usize,
    k_up_red: RCsr,
    sys_u: AmgSystem,
    phi_pc: PhiPc,
    d_scale: &'a [f64],
}

impl Preconditioner<f64> for PiezoSchurPc<'_> {
    fn apply(&self, r_hat: &[f64], z_hat: &mut [f64]) {
        cancel_point();
        let n_u = self.n_u_red;
        let n_p = self.n_phi_red;
        let mut r_u = vec![0.0; n_u];
        let mut r_p = vec![0.0; n_p];
        for i in 0..n_u {
            r_u[i] = r_hat[i] / self.d_scale[i];
        }
        for p in 0..n_p {
            r_p[p] = r_hat[n_u + p] / self.d_scale[n_u + p];
        }
        let z_p = self.phi_pc.solve_neg(&r_p);
        let mut ku_zp = vec![0.0; n_u];
        self.k_up_red.spmv(&z_p, &mut ku_zp);
        let mut rhs_u = vec![0.0; n_u];
        for i in 0..n_u {
            rhs_u[i] = r_u[i] - ku_zp[i];
        }
        let z_u = self.sys_u.precondition_free(&rhs_u);
        for i in 0..n_u {
            z_hat[i] = z_u[i] / self.d_scale[i];
        }
        for p in 0..n_p {
            z_hat[n_u + p] = z_p[p] / self.d_scale[n_u + p];
        }
    }
}

fn solve_static_equilibrated(asm: &PiezoStaticAssembly) -> Result<(Vec<f64>, FactorBackend), String> {
    let n_u_red = asm.constraints.n_u_red;
    if asm.k_red.nrows() == 0 {
        return Ok((asm.constraints.g_all.iter().map(|z| z.re).collect(), FactorBackend::None));
    }
    let (k_hat, rhs_hat, d_scale) = equilibrate_real(&asm.k_red, &asm.rhs_red, n_u_red);
    let solver = DirectSolver::new(&k_hat, false)?;
    let x_hat = solver.solve(&rhs_hat);
    let x_red = unscale_solution(&x_hat, &d_scale);
    let g_real: Vec<f64> = asm.constraints.g_all.iter().map(|z| z.re).collect();
    let x_full = asm.constraints.red_all.expand_with(&x_red, &g_real);
    Ok((x_full, solver.stats().backend))
}

/// Solves a 3D simplex static coupled system with Right-FGMRES + [`PiezoSchurPc`].
/// Returns `Ok(None)` when `spec.solver == LinearSolver::Auto` should fall back to direct
/// (non-convergence or AMG setup issue), always propagating cancellation.
fn solve_static_iterative(
    mesh: &SimplexMesh,
    spec: &Problem,
    asm: &PiezoStaticAssembly,
) -> Result<Option<(Vec<f64>, String)>, String> {
    solve_static_iterative_with_phi_threshold(mesh, spec, asm, 50_000)
}

fn solve_static_iterative_with_phi_threshold(
    mesh: &SimplexMesh,
    spec: &Problem,
    asm: &PiezoStaticAssembly,
    phi_amg_threshold: usize,
) -> Result<Option<(Vec<f64>, String)>, String> {
    let n_u_red = asm.constraints.n_u_red;
    let n_phi_red = asm.constraints.n_phi_red;
    let n_tot_red = asm.k_red.nrows();
    if n_tot_red == 0 || n_u_red == 0 || n_phi_red == 0 {
        return Ok(None);
    }

    let n1 = asm.constraints.n1;
    let n_mech = asm.constraints.n_mech;
    let n_tot = asm.constraints.n_tot;
    let fixed_u_full: Vec<usize> = (0..n_mech).filter(|&i| asm.constraints.red_u.index_of[i] == usize::MAX).collect();
    let free_u_expected: Vec<usize> =
        (0..n_mech).filter(|&i| asm.constraints.red_u.index_of[i] != usize::MAX).collect();
    if asm.constraints.red_u.roots != free_u_expected {
        if spec.solver == LinearSolver::Auto {
            return Ok(None);
        }
        return Err("iterative coupled piezoelectric solver requires uncoupled mechanical DOF ordering".into());
    }
    let k_uu_full = crate::piezo_asm::sub_csr(&asm.k_full, 0, n_mech, 0, n_mech);
    let sys_u = match AmgSystem::new(mesh, spec.order, 3, &k_uu_full, &fixed_u_full) {
        Ok(sys) => sys,
        Err(e) if e == femcore::progress::CANCELLED => return Err(e),
        Err(e) => {
            if spec.solver == LinearSolver::Auto {
                return Ok(None);
            }
            return Err(format!("iterative displacement AMG setup failed: {e}"));
        }
    };

    let (k_up_red, m_phi_red, diag_add_red) = build_m_phi_with_coupling(&asm.k_red, n_u_red, n_phi_red);
    let has_hubs = (n_mech..n_tot).any(|g| {
        let r = asm.constraints.red_all.index_of[g];
        r != usize::MAX && asm.constraints.red_all.roots[r] != g
    });
    let (phi_pc, phi_desc) = if !has_hubs && n_phi_red > phi_amg_threshold {
        let fixed_p_full: Vec<usize> = (n_mech..n_tot)
            .filter(|&g| asm.constraints.red_all.index_of[g] == usize::MAX)
            .map(|g| g - n_mech)
            .collect();
        let free_p_expected: Vec<usize> = (n_mech..n_tot)
            .filter(|&g| asm.constraints.red_all.index_of[g] != usize::MAX)
            .map(|g| g - n_mech)
            .collect();
        let can_amg = free_p_expected.len() == n_phi_red
            && (0..n_phi_red).all(|p| free_p_expected[p] == asm.constraints.red_all.roots[n_u_red + p] - n_mech);
        if can_amg {
            let mut diag_add_full = vec![0.0; n1];
            for (p_red, &add) in diag_add_red.iter().enumerate() {
                let g = asm.constraints.red_all.roots[n_u_red + p_red] - n_mech;
                diag_add_full[g] = add;
            }
            let k_pp_full_neg = crate::piezo_asm::sub_csr(&asm.k_full, n_mech, n_tot, n_mech, n_tot);
            let m_phi_full = negate_add_diag_csr(&k_pp_full_neg, &diag_add_full);
            match AmgSystem::new_with(
                mesh,
                spec.order,
                1,
                &m_phi_full,
                &fixed_p_full,
                AmgOptions { line_smooth: true, ..Default::default() },
            ) {
                Ok(sys_p) => (PhiPc::Amg(Box::new(sys_p)), "phi=AMG"),
                Err(e) if e == femcore::progress::CANCELLED => return Err(e),
                _ => (PhiPc::Direct(Box::new(DirectSolver::new(&m_phi_red, true)?)), "phi=direct-SPD"),
            }
        } else {
            (PhiPc::Direct(Box::new(DirectSolver::new(&m_phi_red, true)?)), "phi=direct-SPD")
        }
    } else {
        (PhiPc::Direct(Box::new(DirectSolver::new(&m_phi_red, true)?)), "phi=direct-SPD")
    };

    let (k_hat, rhs_hat, d_scale) = equilibrate_real(&asm.k_red, &asm.rhs_red, n_u_red);
    let rhs_norm = rhs_hat.iter().map(|&v| v * v).sum::<f64>().sqrt();
    if rhs_norm == 0.0 {
        let g_real: Vec<f64> = asm.constraints.g_all.iter().map(|z| z.re).collect();
        let x_full = asm.constraints.red_all.expand_with(&vec![0.0; n_tot_red], &g_real);
        return Ok(Some((x_full, format!("Schur-AMG Right-FGMRES (0 iterations, {phi_desc})"))));
    }

    let pc = PiezoSchurPc { n_u_red, n_phi_red, k_up_red, sys_u, phi_pc, d_scale: &d_scale };
    let opts = GmresOptions {
        tol: Tolerances { rtol: 1e-10, atol: 1e-30, dtol: 1e4, max_it: 400 },
        restart: 50,
        side: PcSide::Right,
        orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
        flexible: true,
        initial_guess_nonzero: false,
        record_history: false,
    };
    let mut gmres = Gmres::<f64>::new(opts);
    let mut x_hat = vec![0.0; n_tot_red];
    let info = catching_cancel(|| gmres.solve(&k_hat, &pc, &rhs_hat, &mut x_hat))?;

    let mut kx = vec![0.0; n_tot_red];
    k_hat.spmv(&x_hat, &mut kx);
    let r_norm = kx.iter().zip(&rhs_hat).map(|(&a, &b)| (a - b).powi(2)).sum::<f64>().sqrt();
    let rel_res = r_norm / rhs_norm.max(1e-300);
    if !info.reason.is_converged() || !rel_res.is_finite() || rel_res > 1e-8 {
        if spec.solver == LinearSolver::Auto {
            return Ok(None);
        }
        return Err(format!(
            "iterative coupled piezoelectric solver did not converge ({:?}, {} iterations, rel_res {:.2e}); choose solver = \"direct\"",
            info.reason, info.iterations, rel_res
        ));
    }

    let x_red = unscale_solution(&x_hat, &d_scale);
    let g_real: Vec<f64> = asm.constraints.g_all.iter().map(|z| z.re).collect();
    let x_full = asm.constraints.red_all.expand_with(&x_red, &g_real);
    Ok(Some((
        x_full,
        format!("Schur-AMG Right-FGMRES ({} iterations, rel_res {:.2e}, {phi_desc})", info.iterations, rel_res),
    )))
}

fn solve_static_simplex(
    mesh: &SimplexMesh,
    spec: &Problem,
    asm: &PiezoStaticAssembly,
) -> Result<(Vec<f64>, String), String> {
    if mesh.dim != 3 || !spec.periodic.is_empty() {
        if spec.solver == LinearSolver::Iterative {
            return Err(
                "iterative Schur-AMG solver is only available for 3D simplex meshes without periodic boundaries; choose solver = \"direct\" or \"auto\""
                    .into(),
            );
        }
        let (x_full, backend) = solve_static_equilibrated(asm)?;
        return Ok((x_full, format!("backend {backend:?}")));
    }
    let want_iter = spec.solver == LinearSolver::Iterative
        || (spec.solver == LinearSolver::Auto && asm.constraints.n_tot > 120_000);
    if want_iter {
        if let Some((x_full, desc)) = solve_static_iterative(mesh, spec, asm)? {
            return Ok((x_full, desc));
        }
    }
    let (x_full, backend) = solve_static_equilibrated(asm)?;
    if want_iter {
        Ok((x_full, format!("backend {backend:?} (iterative fallback)")))
    } else {
        Ok((x_full, format!("backend {backend:?}")))
    }
}

#[allow(clippy::too_many_arguments)]
fn build_eigen_result(
    order: usize,
    d: usize,
    asm: &PiezoEigenAssembly,
    solved: &PiezoModeSolved,
    t0: Instant,
    energy_fn: &dyn Fn(&[C], &[C]) -> (f64, f64),
    fields_fn: &dyn Fn(&[C], &[C], bool) -> (Vec<Field>, Vec<CellField>),
) -> SolveResult {
    let n1 = asm.constraints.n1;
    let n_mech = asm.constraints.n_mech;
    let mut solutions = Vec::with_capacity(solved.freqs.len());
    let mut scalars = Vec::new();
    let mut unstable_count = 0usize;
    let report_complex = asm.has_pml || asm.has_damping || (asm.is_complex && !asm.is_hermitian);

    for (j, (&raw_f, vec_full)) in solved.freqs.iter().zip(&solved.full_vecs).enumerate() {
        let u_raw = &vec_full[..n_mech];
        let p_raw = &vec_full[n_mech..];
        let big = u_raw.iter().copied().chain(p_raw.iter().copied()).fold(c(0.0, 0.0), |m, z| {
            if z.norm() > m.norm() {
                z
            } else {
                m
            }
        });
        let big_u = u_raw.iter().copied().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let anchor = if big_u.norm() > 0.0 { big_u } else { big };
        let rot = if anchor.norm() > 0.0 { anchor.conj() / anchor.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / anchor.norm().max(1e-300);
        let u_norm: Vec<C> = u_raw.iter().map(|&z| z * rot * scale).collect();
        let p_norm: Vec<C> = p_raw.iter().map(|&z| z * rot * scale).collect();

        let (ue, ud) = energy_fn(&u_norm, &p_norm);
        let utot = ue + ud;
        let (mech_share, elec_share) = if utot > 0.0 { (ue / utot, ud / utot) } else { (1.0, 0.0) };

        let (fields, cell_fields) = fields_fn(&u_norm, &p_norm, asm.is_complex);

        let (label, val_f) = if !report_complex {
            let f_val = raw_f.re;
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(f_val)));
            (format!("mode {}: {:.6} Hz", j + 1, f_val), f_val)
        } else if asm.has_pml {
            let scale2 = scale * scale;
            let (layer_fraction, phys_mass_norm) = match (&asm.m1_phys_all, &asm.m1_phys_pml) {
                (Some(ma), Some(mp)) => {
                    let tot = real_quad_norm(ma, u_raw, n1, d) * scale2;
                    let in_pml = real_quad_norm(mp, u_raw, n1, d) * scale2;
                    let lf = if tot > 0.0 { (in_pml / tot).clamp(0.0, 1.0) } else { 0.0 };
                    (lf, (tot - in_pml).max(0.0))
                }
                _ => (0.0, 0.0),
            };
            crate::pml::push_pml_mode_scalars_and_label(
                j,
                raw_f,
                layer_fraction,
                phys_mass_norm,
                asm.has_pml,
                asm.has_interior,
                &mut scalars,
                &mut unstable_count,
            )
        } else {
            let is_amplifying = raw_f.im < -AMPLIFYING_IM_TOL * raw_f.re.abs().max(1e-30);
            if is_amplifying {
                unstable_count += 1;
            }
            let q = if raw_f.im.abs() > 0.0 { raw_f.re / (2.0 * raw_f.im) } else { f64::INFINITY };
            let lbl = if is_amplifying {
                format!(
                    "mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]",
                    j + 1,
                    raw_f.re,
                    raw_f.im,
                    q
                )
            } else {
                format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, raw_f.re, raw_f.im, q)
            };
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", raw_f.re, raw_f.im))));
            if q.is_finite() {
                scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
            }
            scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(raw_f.re)));
            scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(raw_f.im)));
            (lbl, raw_f.re)
        };

        scalars.push((format!("mechanical energy share {}", j + 1), Scalar::Num(mech_share)));
        scalars.push((format!("electrical energy share {}", j + 1), Scalar::Num(elec_share)));

        solutions.push(Solution { label, fields, cell_fields, value: Some(val_f), ..Default::default() });
    }

    let mut log = asm.diagnostics.clone();
    if asm.n_pairs > 0 {
        log.push(format!(
            "coupled piezoelectric Floquet: {} pair(s), {} total dofs -> {} reduced (u={}, phi={}), Krylov dim {}{}",
            asm.n_pairs,
            asm.constraints.n_tot,
            asm.constraints.red_all.n_reduced(),
            asm.constraints.n_u_red,
            asm.constraints.n_phi_red,
            solved.dim,
            if solved.converged { "" } else { " (not every pair converged)" }
        ));
    } else {
        log.push(format!(
            "coupled piezoelectric eigenfrequency ({}, Krylov dimension {}{}): {} total dofs -> {} reduced (u={}, phi={}), {} gauge pin(s)",
            if asm.has_damping {
                "damped QEP Schur"
            } else if asm.has_pml {
                "PML Schur"
            } else if asm.is_complex {
                "lossy complex Schur"
            } else {
                "Hermitian Schur"
            },
            solved.dim,
            if solved.converged { "" } else { ", incomplete" },
            asm.constraints.n_tot,
            asm.constraints.red_all.n_reduced(),
            asm.constraints.n_u_red,
            asm.constraints.n_phi_red,
            asm.constraints.gauge_pins.len()
        ));
    }
    log.push(format!("coupled piezo factorization backend: {:?}", solved.backend));
    if report_complex && (unstable_count > 0 || (asm.has_interior && asm.has_pml)) {
        scalars.push(("unstable_modes".into(), Scalar::Num(unstable_count as f64)));
    }
    if unstable_count > 0 {
        log.push(format!("WARNING: {unstable_count} amplifying mode(s) detected (Im f < -1e-9 Re f)"));
    }

    SolveResult {
        order: order as u8,
        dim: d as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.constraints.n_tot,
        ..Default::default()
    }
}

/// Solves a coupled piezoelectric static or eigenfrequency problem on [`SimplexMesh`] (`Discretisation`).
pub fn solve_piezo(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let mut tl = t0;
    let order = spec.order;
    let d = mesh.dim;
    match &spec.study {
        crate::problem::Study::Static => {
            let asm = {
                let _s = progress::stage("coupled piezoelectric static assembly (simplex)");
                assemble_piezo_static_simplex(mesh, spec)?
            };
            lap(&mut tl, "assembly (piezo static)");
            let (x_full, solver_desc) = {
                let _s = progress::stage("coupled piezoelectric linear solve");
                solve_static_simplex(mesh, spec, &asm)?
            };
            lap(&mut tl, "linear solve");
            let u_c: Vec<C> = x_full[..asm.constraints.n_mech].iter().map(|&v| c(v, 0.0)).collect();
            let p_c: Vec<C> = x_full[asm.constraints.n_mech..].iter().map(|&v| c(v, 0.0)).collect();
            let (ue, ud) = cell_energies_simplex(
                mesh,
                &asm.mats,
                order,
                spec.plane,
                spec.thickness,
                &u_c,
                &p_c,
                &ResolvedPmlDomains::default(),
                asm.constraints.n1,
            );
            let (fields, cell_fields, umax, vmax, max_e) =
                static_fields_simplex(mesh, order, &asm, &x_full, spec.plane);
            let phi_slice = &x_full[asm.constraints.n_mech..];
            let phi_min = phi_slice.iter().cloned().fold(f64::INFINITY, f64::min);
            let phi_max = phi_slice.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let mut scalars = vec![
                ("max |u| (m)".into(), Scalar::Num(umax)),
                ("max von Mises (Pa)".into(), Scalar::Num(vmax)),
                ("min phi (V)".into(), Scalar::Num(phi_min)),
                ("max phi (V)".into(), Scalar::Num(phi_max)),
                ("max |E| (V/m)".into(), Scalar::Num(max_e)),
            ];
            push_static_charge_and_energy_scalars(&asm, &x_full, &mut scalars, ue, ud);
            let mut log = asm.diagnostics.clone();
            log.push(format!(
                "coupled piezoelectric static: {} dofs -> {} reduced (u={}, phi={}), {} gauge pin(s), {}",
                asm.constraints.n_tot,
                asm.constraints.red_all.n_reduced(),
                asm.constraints.n_u_red,
                asm.constraints.n_phi_red,
                asm.constraints.gauge_pins.len(),
                solver_desc
            ));
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions: vec![Solution { label: "static".into(), fields, cell_fields, ..Default::default() }],
                scalars,
                log,
                seconds: t0.elapsed().as_secs_f64(),
                dofs: asm.constraints.n_tot,
                ..Default::default()
            })
        }
        crate::problem::Study::Eigenfrequency { .. } => {
            let asm = assemble_piezo_eigen_simplex(mesh, spec)?;
            lap(&mut tl, "assembly (piezo eigenfrequency)");
            let solved = solve_piezo_eigen_modes(&asm)?;
            lap(&mut tl, "Schur eigensolve");
            let n1 = asm.constraints.n1;
            let energy_fn = |u: &[C], p: &[C]| {
                cell_energies_simplex(mesh, &asm.mats, order, spec.plane, spec.thickness, u, p, &asm.resolved_pml, n1)
            };
            let fields_fn = |u: &[C], p: &[C], complex_mode: bool| {
                let (f, c) = push_eigen_mode_fields_simplex(mesh, order, d, n1, u, p, complex_mode);
                (f, cell_fields_of(&c))
            };
            Ok(build_eigen_result(order, d, &asm, &solved, t0, &energy_fn, &fields_fn))
        }
        _ => Err(format!("piezoelectric solver does not handle study '{}'", spec.study.kind())),
    }
}

/// Solves a coupled piezoelectric static or eigenfrequency problem natively on [`MixedDiscretisation`].
pub fn solve_piezo_mixed(raw_mesh: &Mesh, mixed: &MixedDiscretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let mut tl = t0;
    let order = spec.order;
    let d = mixed.dim;
    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let (simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;

    match &spec.study {
        crate::problem::Study::Static => {
            if spec.solver == LinearSolver::Iterative {
                return Err(
                    "iterative Schur-AMG solver does not support mixed-cell meshes or infinite elements; choose solver = \"direct\" or \"auto\""
                        .into(),
                );
            }
            let asm = {
                let _s = progress::stage("coupled piezoelectric static assembly (mixed)");
                assemble_piezo_static_mixed(mixed, &dof_map, spec)?
            };
            lap(&mut tl, "assembly (piezo mixed static)");
            let (x_full, backend) = {
                let _s = progress::stage("coupled piezoelectric linear solve");
                solve_static_equilibrated(&asm)?
            };
            lap(&mut tl, "linear solve");
            let u_c: Vec<C> = x_full[..asm.constraints.n_mech].iter().map(|&v| c(v, 0.0)).collect();
            let p_c: Vec<C> = x_full[asm.constraints.n_mech..].iter().map(|&v| c(v, 0.0)).collect();
            let (ue, ud) = cell_energies_mixed(
                mixed,
                &dof_map,
                &asm.mats,
                spec.plane,
                spec.thickness,
                &u_c,
                &p_c,
                &ResolvedPmlDomains::default(),
            );
            let (fields, cell_fields, umax, vmax, max_e) =
                static_fields_mixed(mixed, &dof_map, &simplex, &proj, order, &asm, &x_full, spec.plane);
            let phi_slice = &x_full[asm.constraints.n_mech..];
            let phi_min = phi_slice.iter().cloned().fold(f64::INFINITY, f64::min);
            let phi_max = phi_slice.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let mut scalars = vec![
                ("max |u| (m)".into(), Scalar::Num(umax)),
                ("max von Mises (Pa)".into(), Scalar::Num(vmax)),
                ("min phi (V)".into(), Scalar::Num(phi_min)),
                ("max phi (V)".into(), Scalar::Num(phi_max)),
                ("max |E| (V/m)".into(), Scalar::Num(max_e)),
            ];
            push_static_charge_and_energy_scalars(&asm, &x_full, &mut scalars, ue, ud);
            let mut log = asm.diagnostics.clone();
            log.push(format!(
                "coupled piezoelectric mixed static: {} dofs -> {} reduced (u={}, phi={}), {} gauge pin(s), backend {:?}",
                asm.constraints.n_tot,
                asm.constraints.red_all.n_reduced(),
                asm.constraints.n_u_red,
                asm.constraints.n_phi_red,
                asm.constraints.gauge_pins.len(),
                backend
            ));
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions: vec![Solution { label: "static".into(), fields, cell_fields, ..Default::default() }],
                scalars,
                log,
                seconds: t0.elapsed().as_secs_f64(),
                dofs: asm.constraints.n_tot,
                ..Default::default()
            })
        }
        crate::problem::Study::Eigenfrequency { .. } => {
            let asm = assemble_piezo_eigen_mixed(mixed, &dof_map, spec)?;
            lap(&mut tl, "assembly (piezo mixed eigenfrequency)");
            let solved = solve_piezo_eigen_modes(&asm)?;
            lap(&mut tl, "Schur eigensolve");
            let n1 = asm.constraints.n1;
            let energy_fn = |u: &[C], p: &[C]| {
                cell_energies_mixed(mixed, &dof_map, &asm.mats, spec.plane, spec.thickness, u, p, &asm.resolved_pml)
            };
            let fields_fn = |u: &[C], p: &[C], complex_mode: bool| {
                push_eigen_mode_fields_mixed(mixed, &dof_map, &proj, order, d, n1, u, p, complex_mode)
            };
            Ok(build_eigen_result(order, d, &asm, &solved, t0, &energy_fn, &fields_fn))
        }
        _ => Err(format!("piezoelectric solver does not handle study '{}'", spec.study.kind())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::problem::{Bc, BcKind};
    use femcore::mesh::cells::{CellMesh, CellType};

    #[test]
    fn piezo_schur_pc_phi_amg_matches_direct() {
        let mut cm = CellMesh::<3>::new();
        let (nx, ny, nz) = (4, 4, 4);
        let vid = |i: usize, j: usize, k: usize| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    cm.points.push([
                        1e-3 * i as f64 / nx as f64,
                        1e-3 * j as f64 / ny as f64,
                        0.5e-3 * k as f64 / nz as f64,
                    ]);
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    cm.push_cell(
                        CellType::Hex,
                        &[
                            vid(i, j, k),
                            vid(i + 1, j, k),
                            vid(i + 1, j + 1, k),
                            vid(i, j + 1, k),
                            vid(i, j, k + 1),
                            vid(i + 1, j, k + 1),
                            vid(i + 1, j + 1, k + 1),
                            vid(i, j + 1, k + 1),
                        ],
                        1,
                    );
                }
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                cm.push_facet(&[vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0)], 1);
                cm.push_facet(&[vid(i, j, nz), vid(i + 1, j, nz), vid(i + 1, j + 1, nz), vid(i, j + 1, nz)], 2);
            }
        }
        let tet = nyx_physics::discretise(&femcore::mesh::Mesh::D3(cm)).unwrap();
        let mut spec = Problem::new(crate::problem::Study::Static, 1);
        spec.solver = LinearSolver::Iterative;
        spec.materials.set(1, "E", 300e9);
        spec.materials.set(1, "nu", 0.25);
        spec.materials.set(1, "eps_r", 9.0);
        spec.materials.aniso.insert(
            1,
            femcore::material::Anisotropy {
                piezo_e: Some([
                    [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
                    [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
                    [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0],
                ]),
                ..Default::default()
            },
        );
        spec.bcs = vec![
            Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
            Bc { labels: vec![1], kind: BcKind::Ground },
            Bc { labels: vec![2], kind: BcKind::Potential(3.0) },
        ];
        let asm = assemble_piezo_static_simplex(&tet, &spec).unwrap();
        let (x_dir, _) = solve_static_equilibrated(&asm).unwrap();
        let (x_amg, desc) = solve_static_iterative_with_phi_threshold(&tet, &spec, &asm, 10).unwrap().unwrap();
        assert!(desc.contains("phi=AMG"), "expected phi=AMG path, got {desc}");
        let max_u_diff = x_dir[..asm.constraints.n_mech]
            .iter()
            .zip(&x_amg[..asm.constraints.n_mech])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        let u_scale = x_dir[..asm.constraints.n_mech].iter().map(|v| v.abs()).fold(0.0_f64, f64::max).max(1e-30);
        assert!(max_u_diff / u_scale < 1e-7, "max_u_diff={max_u_diff} scale={u_scale}");
        let max_p_diff = x_dir[asm.constraints.n_mech..]
            .iter()
            .zip(&x_amg[asm.constraints.n_mech..])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_p_diff < 1e-8, "max_p_diff={max_p_diff}");
    }
}
