// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `SolverDiagnostics`, `NumericalErrorBudget`, and checked linear solver helpers for `nyx-dfs`.

use femcore::fem::{amg_wanted, free_dofs, AmgSystem};
use femcore::linalg::amg_elastic::AmgOptions;
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{backward_residual_many, DirectSolveCheck, DirectSolver};
use femcore::props::LinearSolver;
use nyx_physics::{
    DiagnosticSeverity, DiagnosticWarning, Discretisation, NumericalErrorBudget, SolverDiagnostics,
    SpatialErrorEstimate,
};

/// Accumulator for solver telemetry and error-budget construction across `nyx-dfs` studies.
#[derive(Clone, Debug)]
pub(crate) struct DfsDiag {
    pub solver_kind: String,
    pub converged: bool,
    pub iterations: Option<usize>,
    pub restarts: Option<usize>,
    pub max_algebraic_residual: f64,
    pub per_mode_residuals: Vec<f64>,
    pub rcond_estimate: Option<f64>,
    pub cholesky_fallbacks: usize,
    pub warnings: Vec<DiagnosticWarning>,
    pub spatial_discretization: Option<SpatialErrorEstimate>,
    pub temporal_rel: Option<f64>,
    pub model_truncation_rel: Option<f64>,
}

impl DfsDiag {
    pub(crate) fn new(solver_kind: impl Into<String>) -> Self {
        Self {
            solver_kind: solver_kind.into(),
            converged: true,
            iterations: None,
            restarts: None,
            max_algebraic_residual: 0.0,
            per_mode_residuals: Vec::new(),
            rcond_estimate: None,
            cholesky_fallbacks: 0,
            warnings: Vec::new(),
            spatial_discretization: None,
            temporal_rel: None,
            model_truncation_rel: None,
        }
    }

    pub(crate) fn push_warning(&mut self, code: &str, severity: DiagnosticSeverity, message: impl Into<String>) {
        self.warnings.push(DiagnosticWarning { code: code.to_string(), severity, message: message.into() });
    }

    pub(crate) fn absorb_direct(&mut self, chk: &DirectSolveCheck) {
        self.solver_kind = chk.solver_kind.to_string();
        if !chk.algebraic_residual.is_finite() {
            self.max_algebraic_residual = f64::INFINITY;
            self.converged = false;
        } else {
            self.max_algebraic_residual = self.max_algebraic_residual.max(chk.algebraic_residual);
        }
        self.cholesky_fallbacks += chk.cholesky_fallbacks;
        if let Some(rc) = chk.rcond_estimate {
            self.rcond_estimate = Some(self.rcond_estimate.map_or(rc, |prev| prev.min(rc)));
            if rc < 1e-12 {
                self.push_warning(
                    "ill_conditioned_system",
                    DiagnosticSeverity::Warning,
                    format!("estimated 1-norm rcond = {rc:.2e} (< 1e-12), pivot growth = {:.2e}", chk.pivot_growth),
                );
            }
        }
    }

    pub(crate) fn record_amg(&mut self, its: usize, rel: f64, converged: bool) {
        self.solver_kind = "amg_pcg".to_string();
        self.iterations = Some(self.iterations.map_or(its, |prev| prev + its));
        if !rel.is_finite() {
            self.max_algebraic_residual = f64::INFINITY;
        } else {
            self.max_algebraic_residual = self.max_algebraic_residual.max(rel);
        }
        self.converged &= converged;
    }

    pub(crate) fn set_spatial_estimate(&mut self, est: SpatialErrorEstimate) {
        self.spatial_discretization = Some(est);
    }

    pub(crate) fn set_temporal_rel(&mut self, rel: f64) {
        self.temporal_rel = Some(rel);
    }

    pub(crate) fn finish(mut self) -> (SolverDiagnostics, NumericalErrorBudget) {
        let algebraic_rel =
            if self.max_algebraic_residual.is_finite() { self.max_algebraic_residual.max(0.0) } else { f64::INFINITY };

        let mut disc_sq = 0.0_f64;
        let mut has_disc = false;
        if let Some(sp) = &self.spatial_discretization {
            if sp.rel_bound.is_finite() {
                disc_sq += sp.rel_bound * sp.rel_bound;
                has_disc = true;
            }
        }
        if let Some(t) = self.temporal_rel {
            if t.is_finite() {
                disc_sq += t * t;
                has_disc = true;
            }
        }
        if let Some(m) = self.model_truncation_rel {
            if m.is_finite() {
                disc_sq += m * m;
                has_disc = true;
            }
        }
        let disc_rel = has_disc.then(|| disc_sq.sqrt());
        let combined_num_uncertainty_rel = disc_rel.map(|d| d.hypot(algebraic_rel));
        let alg_below_disc_ratio = disc_rel.filter(|&d| d > 0.0).map(|d| algebraic_rel / d);

        if let Some(ratio) = alg_below_disc_ratio {
            if ratio > 0.01 || !ratio.is_finite() {
                let d_val = disc_rel.unwrap_or(0.0);
                self.push_warning(
                    "algebraic_floor_exceeded",
                    DiagnosticSeverity::Warning,
                    format!(
                        "algebraic residual ({algebraic_rel:.3e}) exceeds 1% of discretization/temporal error ({d_val:.3e}, ratio {ratio:.3e})"
                    ),
                );
            }
        }
        if !self.converged && !self.warnings.iter().any(|w| w.severity == DiagnosticSeverity::Error) {
            self.push_warning(
                "solver_not_converged",
                DiagnosticSeverity::Error,
                format!("solver '{}' did not converge (max residual = {algebraic_rel:.3e})", self.solver_kind),
            );
        } else if !algebraic_rel.is_finite() || algebraic_rel > 1e-6 {
            let sev = if !algebraic_rel.is_finite() || algebraic_rel > 1e-3 {
                DiagnosticSeverity::Error
            } else {
                DiagnosticSeverity::Warning
            };
            self.push_warning(
                "high_algebraic_residual",
                sev,
                format!("algebraic relative residual is {algebraic_rel:.3e} (> 1e-6)"),
            );
        }

        let diag = SolverDiagnostics {
            solver_kind: self.solver_kind,
            converged: self.converged,
            iterations: self.iterations,
            restarts: self.restarts,
            max_algebraic_residual: algebraic_rel,
            per_mode_residuals: self.per_mode_residuals,
            rcond_estimate: self.rcond_estimate,
            cholesky_fallbacks: self.cholesky_fallbacks,
            warnings: self.warnings,
        };
        let budget = NumericalErrorBudget {
            algebraic_rel,
            spatial_discretization: self.spatial_discretization,
            temporal_rel: self.temporal_rel,
            model_truncation_rel: self.model_truncation_rel,
            combined_num_uncertainty_rel,
            alg_below_disc_ratio,
        };
        (diag, budget)
    }
}

/// Solves a Dirichlet-constrained SPD linear system via `DirectSolver::solve_many_checked`.
pub(crate) fn solve_spd_direct_checked(
    k: &Csr,
    f: &[f64],
    fixed: &[usize],
    vals: &[f64],
) -> Result<(Vec<f64>, DirectSolveCheck, usize), String> {
    let n = k.nrows();
    let mut u = vec![0.0; n];
    for (&d, &v) in fixed.iter().zip(vals) {
        u[d] = v;
    }
    let free = free_dofs(n, fixed);
    let n_free = free.len();
    if n_free == 0 {
        let zero_a = Csr::zeros(0, 0);
        let dummy_chk = DirectSolveCheck {
            algebraic_residual: 0.0,
            rcond_estimate: None,
            stats: DirectSolver::new(&Csr::from_triplets(1, 1, &[(0, 0, 1.0)]).unwrap(), true)?.stats(),
            diagnostics: Default::default(),
            cholesky_fallbacks: 0,
            pivot_growth: 1.0,
            solver_kind: "direct_chol",
        };
        let _ = zero_a;
        return Ok((u, dummy_chk, 0));
    }
    let kff = k.submatrix(&free, &free);
    let rhs: Vec<f64> = if fixed.is_empty() || vals.iter().all(|&v| v == 0.0) {
        free.iter().map(|&i| f[i]).collect()
    } else {
        let kfx = k.submatrix(&free, fixed);
        let kx = kfx.matvec(vals);
        free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect()
    };
    let solver = DirectSolver::new(&kff, true)?;
    let (u_free, chk) = solver.solve_checked(&kff, &rhs);
    for (a, &i) in free.iter().enumerate() {
        u[i] = u_free[a];
    }
    Ok((u, chk, n_free))
}

/// Solves `K u = f` with `u[fixed] = values` via AMG-CG or direct sparse Cholesky, enforcing
/// `worst_rel <= tol * 2.0` on iterative solves (falling back to direct if `Auto` or affordable,
/// else marking `converged = false` + Error warning) and populating [`DfsDiag`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_spd_checked(
    mesh: &Discretisation,
    order: usize,
    ncomp: usize,
    solver_pref: LinearSolver,
    k: &Csr,
    f: &[f64],
    fixed: &[usize],
    values: &[f64],
    tol: f64,
    line_smooth: bool,
) -> Result<(Vec<f64>, Option<String>, DfsDiag), String> {
    let n = k.nrows();
    let free = free_dofs(n, fixed);
    let n_free = free.len();
    if n_free == 0 {
        let mut u = vec![0.0; n];
        for (&d, &v) in fixed.iter().zip(values) {
            u[d] = v;
        }
        return Ok((u, Some("0 free dofs: determined by boundary constraints".into()), DfsDiag::new("direct_chol")));
    }
    let use_amg = match solver_pref {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => true,
        LinearSolver::Auto => amg_wanted(mesh.dim, n_free),
    } && !fixed.is_empty();

    if use_amg {
        let opts = AmgOptions { line_smooth, ..Default::default() };
        let sys = AmgSystem::new_with(mesh, order, ncomp, k, fixed, opts)?;
        let mut uvec = vec![0.0; n];
        for (&d, &v) in fixed.iter().zip(values) {
            uvec[d] = v;
        }
        let ku = k.matvec(&uvec);
        let rhs_full: Vec<f64> = (0..n).map(|i| f[i] - ku[i]).collect();
        let fv: Vec<(usize, f64)> = fixed.iter().copied().zip(values.iter().copied()).collect();
        let (mut x, its, rel) = sys.solve_with(None, &rhs_full, &fv, tol, 800);
        femcore::progress::check()?;
        for (&d, &v) in fixed.iter().zip(values) {
            x[d] = v;
        }
        let kff = k.submatrix(&free, &free);
        let rhs_free: Vec<f64> = if fixed.is_empty() || values.iter().all(|&v| v == 0.0) {
            free.iter().map(|&i| f[i]).collect()
        } else {
            let kfx = k.submatrix(&free, fixed);
            let kx = kfx.matvec(values);
            free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect()
        };
        let x_free: Vec<f64> = free.iter().map(|&i| x[i]).collect();
        let bw_rel = backward_residual_many(&kff, std::slice::from_ref(&rhs_free), std::slice::from_ref(&x_free));
        let worst_rel = rel.max(bw_rel);
        let mut diag = DfsDiag::new("amg_pcg");
        if worst_rel.is_finite() && worst_rel <= tol * 2.0 {
            diag.record_amg(its, worst_rel, true);
            let note = format!(
                "AMG-CG ({} levels {:?}, {} iterations, residual {:.1e})",
                sys.n_levels(),
                sys.sizes(),
                its,
                worst_rel
            );
            return Ok((x, Some(note), diag));
        }
        if solver_pref == LinearSolver::Auto || n_free <= 50_000 {
            let direct_solver = DirectSolver::new(&kff, true)?;
            let (u_free, chk) = direct_solver.solve_checked(&kff, &rhs_free);
            let mut u = vec![0.0; n];
            for (&d, &v) in fixed.iter().zip(values) {
                u[d] = v;
            }
            for (a, &i) in free.iter().enumerate() {
                u[i] = u_free[a];
            }
            diag.absorb_direct(&chk);
            diag.push_warning(
                "amg_fallback_to_direct",
                DiagnosticSeverity::Warning,
                format!(
                    "AMG-CG residual {rel:.2e} exceeded tight bound {:.2e}; fell back to direct Cholesky",
                    tol * 2.0
                ),
            );
            let note =
                format!("AMG-CG residual {rel:.2e} exceeded 2*tol; direct Cholesky fallback ({n_free} free dofs)");
            return Ok((u, Some(note), diag));
        }
        diag.record_amg(its, worst_rel, false);
        diag.push_warning(
            "amg_not_converged",
            DiagnosticSeverity::Error,
            format!("AMG-CG failed to reach 2*tol ({:.2e}): residual = {rel:.2e} after {its} iterations", tol * 2.0),
        );
        let note = format!("AMG-CG NOT converged ({its} iterations, residual {rel:.2e} > {:.2e})", tol * 2.0);
        return Ok((x, Some(note), diag));
    }

    let (u, chk, n_free) = solve_spd_direct_checked(k, f, fixed, values)?;
    let mut diag = DfsDiag::new(chk.solver_kind);
    diag.absorb_direct(&chk);
    let note = format!("direct Cholesky ({n_free} free dofs, eps_alg = {:.2e})", chk.algebraic_residual);
    Ok((u, Some(note), diag))
}

/// Checked variant of `common::solve_spd_mixed`.
pub(crate) fn solve_spd_mixed_checked(
    _solver_pref: LinearSolver,
    k: &Csr,
    f: &[f64],
    fixed: &[usize],
    vals: &[f64],
) -> Result<(Vec<f64>, Option<String>, DfsDiag), String> {
    let (u, chk, n_free) = solve_spd_direct_checked(k, f, fixed, vals)?;
    let mut diag = DfsDiag::new(chk.solver_kind);
    diag.absorb_direct(&chk);
    let note = if n_free == 0 {
        Some("0 free dofs: determined by boundary constraints".into())
    } else {
        Some(format!("direct Cholesky ({n_free} free mixed-cell dofs, eps_alg = {:.2e})", chk.algebraic_residual))
    };
    Ok((u, note, diag))
}

/// Checked replacement for `AmgSystem::solve_free` accepting `matrix: Option<&Csr>`, `n`, and `free`,
/// returning `(u_free, iterations, final_rel_residual)`.
pub(crate) fn amg_solve_free_checked(
    sys: &AmgSystem,
    matrix: Option<&Csr>,
    n: usize,
    free: &[usize],
    b_free: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, usize, f64), String> {
    let mut b = vec![0.0; n];
    for (a, &i) in free.iter().enumerate() {
        b[i] = b_free[a];
    }
    let (x, its, rel) = sys.solve_with(matrix, &b, &[], tol, max_iter);
    femcore::progress::check()?;
    Ok((free.iter().map(|&i| x[i]).collect(), its, rel))
}
