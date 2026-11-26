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

/// Convergence regime classified by Roache's 3-grid Grid Convergence Index (VVUQ Layer 3, §5.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RoacheGciRegime {
    /// `0 < R < 1`: monotonic convergence toward the continuum limit.
    Monotonic,
    /// `-1 < R < 0`: oscillatory damped convergence.
    Oscillatory,
    /// `|R| >= 1` (or degenerate zero changes): non-convergent / unresolved.
    NonConvergent,
}

/// Result of a 3-grid Roache Grid Convergence Index calculation on grids `h1 < h2 < h3`.
#[derive(Clone, Debug, PartialEq)]
pub struct RoacheGciReport {
    /// Classification of the convergence ratio `R = (phi2 - phi1) / (phi3 - phi2)`.
    pub regime: RoacheGciRegime,
    /// Convergence ratio `R = (phi2 - phi1) / (phi3 - phi2)`.
    pub convergence_ratio: f64,
    /// Fine-to-medium effective grid refinement ratio `r21 = (n1 / n2)^(1/d)`.
    pub r21: f64,
    /// Medium-to-coarse effective grid refinement ratio `r32 = (n2 / n3)^(1/d)`.
    pub r32: f64,
    /// Observed order of accuracy `p_obs` (`None` when non-convergent or oscillatory).
    pub observed_order: Option<f64>,
    /// Fine-grid GCI relative uncertainty (`GCI_21` or oscillatory envelope).
    pub gci_21: f64,
    /// Medium-grid GCI (`GCI_32`), if monotonic.
    pub gci_32: Option<f64>,
    /// Applied factor of safety (`1.25` when in asymptotic range, else `3.0`).
    pub safety_factor: f64,
    /// `true` when `|GCI_32 / (r21^p_obs * GCI_21) - 1| <= 0.15` and `p_obs` lies in `[0.5*p_th, 1.5*p_th]`.
    pub in_asymptotic_range: bool,
    /// Generalized Richardson extrapolated continuum value (when monotonic).
    pub extrapolated_value: Option<f64>,
}

/// Computes Roache's 3-grid GCI for scalar QoIs `[phi1, phi2, phi3]` ordered fine-to-coarse (`h1 < h2 < h3`)
/// with effective refinement ratios `r21 > 1.0`, `r32 > 1.0`, and theoretical order `theoretical_order > 0.0`.
pub fn roache_gci(
    phi: [f64; 3],
    r21: f64,
    r32: f64,
    theoretical_order: f64,
) -> Result<RoacheGciReport, String> {
    let [phi1, phi2, phi3] = phi;
    if !phi1.is_finite() || !phi2.is_finite() || !phi3.is_finite() {
        return Err("roache_gci inputs must be finite".into());
    }
    if !r21.is_finite() || r21 <= 1.0 || !r32.is_finite() || r32 <= 1.0 {
        return Err(format!("refinement ratios r21={r21}, r32={r32} must be finite and > 1.0"));
    }
    if !theoretical_order.is_finite() || theoretical_order <= 0.0 {
        return Err(format!("theoretical_order={theoretical_order} must be > 0"));
    }

    let eps21 = phi2 - phi1;
    let eps32 = phi3 - phi2;
    let denom = phi1.abs().max(1e-30);
    if eps21.abs() == 0.0 && eps32.abs() == 0.0 {
        return Ok(RoacheGciReport {
            regime: RoacheGciRegime::Monotonic,
            convergence_ratio: 0.0,
            r21,
            r32,
            observed_order: Some(theoretical_order),
            gci_21: 0.0,
            gci_32: Some(0.0),
            safety_factor: 1.25,
            in_asymptotic_range: true,
            extrapolated_value: Some(phi1),
        });
    }
    if eps32.abs() <= 1e-300 {
        return Ok(RoacheGciReport {
            regime: RoacheGciRegime::NonConvergent,
            convergence_ratio: f64::INFINITY,
            r21,
            r32,
            observed_order: None,
            gci_21: eps21.abs() / denom,
            gci_32: None,
            safety_factor: 3.0,
            in_asymptotic_range: false,
            extrapolated_value: None,
        });
    }

    let r_conv = eps21 / eps32;
    if r_conv > 0.0 && r_conv < 1.0 {
        // Monotonic convergence: fixed-point solve for observed order p_obs
        let ln_r21 = r21.ln();
        let ln_ratio = (eps32 / eps21).abs().ln();
        let s = (eps32 / eps21).signum();
        let mut p = theoretical_order.clamp(0.25, 8.0);
        for _ in 0..80 {
            let num = r21.powf(p) - s;
            let den = r32.powf(p) - s;
            if num <= 0.0 || den <= 0.0 {
                break;
            }
            let q = (num / den).ln();
            let p_next = ((ln_ratio + q).abs() / ln_r21).clamp(0.1, 12.0);
            let done = (p_next - p).abs() <= 1e-10;
            p = p_next;
            if done {
                break;
            }
        }
        let g21_raw = eps21.abs() / (denom * (r21.powf(p) - 1.0));
        let g32_raw = eps32.abs() / (phi2.abs().max(1e-30) * (r32.powf(p) - 1.0));
        let asym_err = if g21_raw > 1e-300 {
            (g32_raw / (r21.powf(p) * g21_raw) - 1.0).abs()
        } else {
            0.0
        };
        let p_in_band = p >= 0.5 * theoretical_order && p <= 1.5 * theoretical_order;
        let in_asymptotic_range = asym_err <= 0.15 && p_in_band;
        let safety_factor = if in_asymptotic_range { 1.25 } else { 3.0 };
        let gci_21 = safety_factor * g21_raw;
        let gci_32 = safety_factor * g32_raw;
        let extrapolated_value = Some(phi1 - eps21 / (r21.powf(p) - 1.0));
        Ok(RoacheGciReport {
            regime: RoacheGciRegime::Monotonic,
            convergence_ratio: r_conv,
            r21,
            r32,
            observed_order: Some(p),
            gci_21,
            gci_32: Some(gci_32),
            safety_factor,
            in_asymptotic_range,
            extrapolated_value,
        })
    } else if r_conv > -1.0 && r_conv <= 0.0 {
        // Oscillatory convergence
        let p_max = phi1.max(phi2).max(phi3);
        let p_min = phi1.min(phi2).min(phi3);
        let gci_21 = 3.0 * (0.5 * (p_max - p_min)) / denom;
        Ok(RoacheGciReport {
            regime: RoacheGciRegime::Oscillatory,
            convergence_ratio: r_conv,
            r21,
            r32,
            observed_order: None,
            gci_21,
            gci_32: None,
            safety_factor: 3.0,
            in_asymptotic_range: false,
            extrapolated_value: None,
        })
    } else {
        // |R| >= 1: non-convergent
        let p_max = phi1.max(phi2).max(phi3);
        let p_min = phi1.min(phi2).min(phi3);
        let rel = (p_max - p_min).max(eps21.abs()) / denom;
        Ok(RoacheGciReport {
            regime: RoacheGciRegime::NonConvergent,
            convergence_ratio: r_conv,
            r21,
            r32,
            observed_order: None,
            gci_21: rel,
            gci_32: None,
            safety_factor: 3.0,
            in_asymptotic_range: false,
            extrapolated_value: None,
        })
    }
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

    #[test]
    fn roache_gci_classifies_monotonic_oscillatory_and_non_convergent() {
        // Pure O(h^2) sequence with exact continuum value 1.0, r=2.0:
        let phi = [1.0 + 0.04, 1.0 + 0.16, 1.0 + 0.64];
        let mono = roache_gci(phi, 2.0, 2.0, 2.0).unwrap();
        assert_eq!(mono.regime, RoacheGciRegime::Monotonic);
        assert!((mono.observed_order.unwrap() - 2.0).abs() < 1e-9);
        assert!(mono.in_asymptotic_range);
        assert_eq!(mono.safety_factor, 1.25);
        assert!((mono.extrapolated_value.unwrap() - 1.0).abs() < 1e-12);

        // Oscillatory sequence: R = (0.1 - (-0.05)) / (-0.3 - 0.1) = -0.375
        let osc = roache_gci([0.95, 1.1, 0.7], 2.0, 2.0, 2.0).unwrap();
        assert_eq!(osc.regime, RoacheGciRegime::Oscillatory);
        assert_eq!(osc.safety_factor, 3.0);
        assert!(osc.observed_order.is_none());

        // Non-convergent sequence: |R| >= 1
        let nonconv = roache_gci([1.0, 1.5, 1.6], 2.0, 2.0, 2.0).unwrap();
        assert_eq!(nonconv.regime, RoacheGciRegime::NonConvergent);
    }
}
