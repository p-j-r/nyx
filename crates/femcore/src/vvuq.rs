// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared VVUQ Layer 2 & Layer 3 numerical verification primitives:
//! zero-sized monomorphized [`SourceField`] / [`UniformSource`] for true PDE-level
//! Method of Manufactured Solutions (MMS) without Cargo feature thrashing, and
//! [`assert_convergence_rate`] for formal h-, p-, and dt-convergence-order checks.

/// Zero-sized default compiles away to zero overhead in production builds;
/// integration tests pass a concrete source (or [`FnSource`]) sharing 100% of
/// production quadrature, basis functions, and Piola transforms.
pub trait SourceField<const D: usize>: Sync {
    /// Real volumetric forcing vector at point `x`, time `t`, and domain label `domain`.
    /// Scalar physics (`dfs::electro`, `dfs::heat`) read component `[0]`.
    fn volume_rhs_re(&self, x: [f64; D], t: f64, domain: u32) -> [f64; 3];

    /// Imaginary volumetric forcing vector (complex harmonic / Maxwell MMS).
    fn volume_rhs_im(&self, _x: [f64; D], _t: f64, _domain: u32) -> [f64; 3] {
        [0.0; 3]
    }

    /// Real Dirichlet trace at boundary point `x`, time `t`, and boundary label `bnd`.
    fn dirichlet_re(&self, x: [f64; D], t: f64, bnd: u32) -> [f64; 3];

    /// Imaginary Dirichlet trace.
    fn dirichlet_im(&self, _x: [f64; D], _t: f64, _bnd: u32) -> [f64; 3] {
        [0.0; 3]
    }

    /// Real Neumann / Robin boundary datum at boundary point `x`, outward unit normal `n`,
    /// time `t`, and boundary label `bnd`.
    fn neumann_robin_re(&self, x: [f64; D], n: [f64; D], t: f64, bnd: u32) -> [f64; 3];

    /// Imaginary Neumann / Robin boundary datum.
    fn neumann_robin_im(&self, _x: [f64; D], _n: [f64; D], _t: f64, _bnd: u32) -> [f64; 3] {
        [0.0; 3]
    }

    /// Returns `true` for the zero-sized [`UniformSource`] production fast path so callers
    /// can preserve bit-identical production behavior when no manufactured source is supplied.
    fn is_uniform(&self) -> bool {
        false
    }
}

/// Production zero-sized source (all MMS callbacks compile away).
#[derive(Clone, Copy, Debug, Default)]
pub struct UniformSource;

impl<const D: usize> SourceField<D> for UniformSource {
    #[inline]
    fn volume_rhs_re(&self, _x: [f64; D], _t: f64, _domain: u32) -> [f64; 3] {
        [0.0; 3]
    }

    #[inline]
    fn dirichlet_re(&self, _x: [f64; D], _t: f64, _bnd: u32) -> [f64; 3] {
        [0.0; 3]
    }

    #[inline]
    fn neumann_robin_re(&self, _x: [f64; D], _n: [f64; D], _t: f64, _bnd: u32) -> [f64; 3] {
        [0.0; 3]
    }

    #[inline]
    fn is_uniform(&self) -> bool {
        true
    }
}

/// Closure-backed [`SourceField`] adapter for MMS integration tests.
pub struct FnSource<FV, FD, FN> {
    pub volume_rhs: FV,
    pub dirichlet: FD,
    pub neumann_robin: FN,
}

impl<const D: usize, FV, FD, FN> SourceField<D> for FnSource<FV, FD, FN>
where
    FV: Fn([f64; D], f64, u32) -> [f64; 3] + Sync,
    FD: Fn([f64; D], f64, u32) -> [f64; 3] + Sync,
    FN: Fn([f64; D], [f64; D], f64, u32) -> [f64; 3] + Sync,
{
    fn volume_rhs_re(&self, x: [f64; D], t: f64, domain: u32) -> [f64; 3] {
        (self.volume_rhs)(x, t, domain)
    }

    fn dirichlet_re(&self, x: [f64; D], t: f64, bnd: u32) -> [f64; 3] {
        (self.dirichlet)(x, t, bnd)
    }

    fn neumann_robin_re(&self, x: [f64; D], n: [f64; D], t: f64, bnd: u32) -> [f64; 3] {
        (self.neumann_robin)(x, n, t, bnd)
    }
}

/// Observed convergence summary across `>= 3` grids or step sizes.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvergenceReport {
    /// Pairwise slopes `ln(e_{k-1} / e_k) / ln(h_{k-1} / h_k)`.
    pub pairwise_slopes: Vec<f64>,
    /// Ordinary least-squares slope of `ln(err)` vs `ln(h)`.
    pub ols_slope: f64,
    /// Finest-pair slope `pairwise_slopes.last()`.
    pub asymptotic_slope: f64,
}

/// Computes pairwise and OLS log-log convergence rates over strictly decreasing positive
/// mesh/step sizes `hs` (`hs.len() >= 3`) and positive finite errors `errs`.
pub fn convergence_rate(hs: &[f64], errs: &[f64]) -> Result<ConvergenceReport, String> {
    if hs.len() != errs.len() {
        return Err(format!(
            "hs.len() ({}) != errs.len() ({})",
            hs.len(),
            errs.len()
        ));
    }
    if hs.len() < 3 {
        return Err(format!(
            "convergence_rate requires at least 3 refinement levels, got {}",
            hs.len()
        ));
    }
    for (i, (&h, &e)) in hs.iter().zip(errs.iter()).enumerate() {
        if !h.is_finite() || h <= 0.0 {
            return Err(format!("hs[{i}] = {h} must be finite and > 0"));
        }
        if !e.is_finite() || e <= 0.0 {
            return Err(format!("errs[{i}] = {e} must be finite and > 0"));
        }
        if i > 0 && hs[i] >= hs[i - 1] {
            return Err(format!(
                "hs must be strictly decreasing: hs[{}]={} >= hs[{}]={}",
                i,
                hs[i],
                i - 1,
                hs[i - 1]
            ));
        }
    }

    let mut pairwise = Vec::with_capacity(hs.len() - 1);
    for i in 1..hs.len() {
        let p = (errs[i - 1] / errs[i]).ln() / (hs[i - 1] / hs[i]).ln();
        pairwise.push(p);
    }

    let n = hs.len() as f64;
    let xs: Vec<f64> = hs.iter().map(|&h| h.ln()).collect();
    let ys: Vec<f64> = errs.iter().map(|&e| e.ln()).collect();
    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        num += dx * (y - mean_y);
        den += dx * dx;
    }
    if den <= 0.0 || !den.is_finite() {
        return Err("degenerate step-size variance in OLS slope".into());
    }
    let ols_slope = num / den;
    let asymptotic_slope = *pairwise.last().unwrap_or(&ols_slope);
    Ok(ConvergenceReport {
        pairwise_slopes: pairwise,
        ols_slope,
        asymptotic_slope,
    })
}

/// Asserts strictly decreasing positive `hs`, finite positive `errs`, and that both
/// the OLS log-log slope and the finest pairwise slope fall within `band` around `expected_p`.
pub fn assert_convergence_rate(
    hs: &[f64],
    errs: &[f64],
    expected_p: f64,
    band: (f64, f64),
) -> ConvergenceReport {
    let report = convergence_rate(hs, errs).unwrap_or_else(|e| {
        panic!("invalid convergence inputs for expected p={expected_p}: {e} (hs={hs:?}, errs={errs:?})")
    });
    assert!(
        report.ols_slope >= band.0 && report.ols_slope <= band.1,
        "OLS convergence order {:.4} outside [{:.3}, {:.3}] (expected {expected_p}, pairwise={:?}, hs={hs:?}, errs={errs:?})",
        report.ols_slope,
        band.0,
        band.1,
        report.pairwise_slopes,
    );
    assert!(
        report.asymptotic_slope >= band.0 && report.asymptotic_slope <= band.1,
        "asymptotic convergence order {:.4} outside [{:.3}, {:.3}] (expected {expected_p}, pairwise={:?}, hs={hs:?}, errs={errs:?})",
        report.asymptotic_slope,
        band.0,
        band.1,
        report.pairwise_slopes,
    );
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_source_is_zero_sized_and_reports_uniform() {
        assert_eq!(std::mem::size_of::<UniformSource>(), 0);
        let u = UniformSource;
        assert!(<UniformSource as SourceField<3>>::is_uniform(&u));
        assert_eq!(<UniformSource as SourceField<3>>::volume_rhs_re(&u, [0.0; 3], 0.0, 1), [0.0; 3]);
    }

    #[test]
    fn assert_convergence_rate_recovers_exact_order_and_rejects_bad_inputs() {
        let hs = [0.25, 0.125, 0.0625, 0.03125];
        let errs: Vec<f64> = hs.iter().map(|&h| 0.4 * h * h).collect();
        let rep = assert_convergence_rate(&hs, &errs, 2.0, (1.85, 2.25));
        assert!((rep.ols_slope - 2.0).abs() < 1e-12);
        assert!((rep.asymptotic_slope - 2.0).abs() < 1e-12);
        assert!(convergence_rate(&[0.2, 0.1], &[0.1, 0.05]).is_err());
        assert!(convergence_rate(&[0.1, 0.2, 0.05], &[0.1, 0.05, 0.01]).is_err());
    }
}
