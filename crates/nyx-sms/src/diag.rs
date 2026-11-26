// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `SolverDiagnostics` and `NumericalErrorBudget` builders for `nyx-sms`.

use femcore::linalg::driver::{CEigPairs, DirectSolveCheck, EigPairs};
use nyx_physics::{
    DiagnosticSeverity, DiagnosticWarning, NumericalErrorBudget, SolverDiagnostics, SpatialErrorEstimate,
};

/// Accumulator for solver telemetry and error-budget construction across `nyx-sms` studies.
#[derive(Clone, Debug)]
pub(crate) struct SmsDiag {
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

impl SmsDiag {
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

    pub(crate) fn from_eig_pairs(solver_kind: impl Into<String>, pairs: &EigPairs) -> Self {
        let max_res = pairs.residuals.iter().copied().fold(0.0_f64, |a, b| {
            if a.is_finite() && b.is_finite() {
                a.max(b)
            } else {
                f64::INFINITY
            }
        });
        Self {
            solver_kind: solver_kind.into(),
            converged: pairs.converged,
            iterations: None,
            restarts: Some(pairs.restarts),
            max_algebraic_residual: max_res,
            per_mode_residuals: pairs.residuals.clone(),
            rcond_estimate: None,
            cholesky_fallbacks: 0,
            warnings: Vec::new(),
            spatial_discretization: None,
            temporal_rel: None,
            model_truncation_rel: None,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn from_complex_eig_pairs(solver_kind: impl Into<String>, pairs: &CEigPairs) -> Self {
        let max_res = pairs.residuals.iter().copied().fold(0.0_f64, |a, b| {
            if a.is_finite() && b.is_finite() {
                a.max(b)
            } else {
                f64::INFINITY
            }
        });
        Self {
            solver_kind: solver_kind.into(),
            converged: pairs.converged,
            iterations: None,
            restarts: Some(pairs.restarts),
            max_algebraic_residual: max_res,
            per_mode_residuals: pairs.residuals.clone(),
            rcond_estimate: None,
            cholesky_fallbacks: 0,
            warnings: Vec::new(),
            spatial_discretization: None,
            temporal_rel: None,
            model_truncation_rel: None,
        }
    }

    pub(crate) fn record_mode_residuals(
        &mut self,
        solver_kind: impl Into<String>,
        residuals: Vec<f64>,
        converged: bool,
        restarts: Option<usize>,
    ) {
        self.solver_kind = solver_kind.into();
        let max_res = residuals.iter().copied().fold(0.0_f64, |a, b| {
            if a.is_finite() && b.is_finite() {
                a.max(b)
            } else {
                f64::INFINITY
            }
        });
        self.max_algebraic_residual = self.max_algebraic_residual.max(max_res);
        self.per_mode_residuals = residuals;
        self.converged &= converged;
        if let Some(r) = restarts {
            self.restarts = Some(self.restarts.map_or(r, |prev| prev + r));
        }
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
                        "algebraic residual ({algebraic_rel:.3e}) exceeds 1% of discretization/truncation error ({d_val:.3e}, ratio {ratio:.3e})"
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
            let sev = if !algebraic_rel.is_finite() || algebraic_rel > 1e-4 {
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

/// Rigal-Gaches normwise generalized eigenproblem backward error
/// `||A x - lambda B x||_2 / ((||A||_F + |lambda| ||B||_F) ||x||_2 + 1e-300)`.
/// Remains small (`O(eps_mach)`) for exact zero / rigid-body modes where vectorwise `||A x||` vanishes.
pub(crate) fn rigal_gaches_complex_residuals(
    a: &femcore::linalg::csr::CCsr,
    b: &femcore::linalg::csr::CCsr,
    values: &[femcore::linalg::csr::C],
    vectors: &[Vec<femcore::linalg::csr::C>],
) -> Vec<f64> {
    let n = a.nrows();
    if values.is_empty() || n == 0 {
        return Vec::new();
    }
    let a_fro = a.values().iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    let b_fro = b.values().iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    let mut ax = vec![femcore::linalg::csr::C::new(0.0, 0.0); n];
    let mut bx = vec![femcore::linalg::csr::C::new(0.0, 0.0); n];
    values
        .iter()
        .zip(vectors)
        .map(|(&lam, x)| {
            let x_norm = x.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            if x_norm <= 0.0 || !x_norm.is_finite() || !lam.is_finite() {
                return f64::INFINITY;
            }
            a.spmv(x, &mut ax);
            b.spmv(x, &mut bx);
            let rr2: f64 = ax.iter().zip(&bx).map(|(&ai, &bi)| (ai - lam * bi).norm_sqr()).sum();
            let denom = (a_fro + lam.norm() * b_fro) * x_norm + 1e-300;
            rr2.sqrt() / denom
        })
        .collect()
}
