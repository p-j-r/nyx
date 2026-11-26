// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Deterministic Latin Hypercube Sampling (LHS) and low-discrepancy Sobol Quasi-Monte Carlo (QMC) point generators
//! seeded with `femcore::linalg::rng::SplitMix64` (§6.4).

use femcore::linalg::rng::SplitMix64;
use nyx_physics::Error;
#[cfg(test)]
use std::f64::consts::FRAC_1_SQRT_2;

/// Maximum samples permitted in a single streaming QMC `sweep_run` request.
pub const MAX_QMC_SAMPLES: usize = 4096;

/// Sampling strategy for parametric design-of-experiments in `sweep_run`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QmcMethod {
    /// Stratified Latin Hypercube Sampling driven by [`SplitMix64`].
    Lhs,
    /// Low-discrepancy base-2 Sobol sequence (with optional `SplitMix64` Owen-style XOR digital shift).
    Sobol,
}

impl QmcMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Lhs => "lhs",
            Self::Sobol => "sobol",
        }
    }
}

/// Marginal distribution mapping one unit-hypercube coordinate `u in (0, 1)` onto a scalar parameter value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MarginDistribution {
    Uniform { min: f64, max: f64 },
    Normal { mean: f64, sigma: f64, min: Option<f64>, max: Option<f64> },
}

impl MarginDistribution {
    /// Maps a unit-hypercube coordinate `u in (0, 1)` deterministically to the physical parameter value.
    pub fn map_unit(self, u: f64) -> f64 {
        let u_clamped = u.clamp(1.0e-12, 1.0 - 1.0e-12);
        match self {
            Self::Uniform { min, max } => min + u_clamped * (max - min),
            Self::Normal { mean, sigma, min, max } => {
                let mut val = mean + sigma * inv_normal_cdf(u_clamped);
                if let Some(lo) = min {
                    val = val.max(lo);
                }
                if let Some(hi) = max {
                    val = val.min(hi);
                }
                val
            }
        }
    }
}

/// Generates `n_samples` unit-hypercube vectors in `(0, 1)^dim` via stratified Latin Hypercube Sampling (`SplitMix64`).
#[allow(clippy::needless_range_loop)]
pub fn generate_lhs_unit(n_samples: usize, dim: usize, seed: u64) -> Result<Vec<Vec<f64>>, Error> {
    if n_samples == 0 || n_samples > MAX_QMC_SAMPLES {
        return Err(Error::invalid(format!(
            "LHS sampling requires 1 <= n_samples <= {MAX_QMC_SAMPLES}, got {n_samples}"
        )));
    }
    if dim == 0 {
        return Err(Error::invalid("LHS sampling requires at least one parameter dimension"));
    }

    let mut rng = SplitMix64::new(seed);
    let mut table = vec![vec![0.0_f64; dim]; n_samples];
    let n_f = n_samples as f64;
    let mut perm: Vec<usize> = (0..n_samples).collect();

    for d in 0..dim {
        // Fisher-Yates shuffle on perm[0..n_samples] using SplitMix64
        for i in (1..n_samples).rev() {
            let j = (rng.next_u64() % ((i + 1) as u64)) as usize;
            perm.swap(i, j);
        }
        for (dst_row, &stratum) in table.iter_mut().zip(&perm) {
            let jitter = uniform_open01(&mut rng);
            dst_row[d] = ((stratum as f64) + jitter) / n_f;
        }
    }
    Ok(table)
}

/// Maximum number of dimensions supported by the embedded primitive polynomial / direction-number table.
pub const MAX_SOBOL_DIM: usize = 21;

struct SobolDimSpec {
    degree: usize,
    poly_a: u32,
    init_m: &'static [u32],
}

// Primitive polynomials over GF(2) paired with verified initial direction odd integers (`m_k < 2^k`)
// covering dimensions 1..20 alongside Van der Corput dimension 0 (total 21 dimensions).
const SOBOL_TABLE: &[SobolDimSpec] = &[
    SobolDimSpec { degree: 1, poly_a: 0, init_m: &[1] },
    SobolDimSpec { degree: 2, poly_a: 1, init_m: &[1, 3] },
    SobolDimSpec { degree: 3, poly_a: 1, init_m: &[1, 3, 1] },
    SobolDimSpec { degree: 3, poly_a: 2, init_m: &[1, 1, 1] },
    SobolDimSpec { degree: 4, poly_a: 1, init_m: &[1, 1, 3, 3] },
    SobolDimSpec { degree: 4, poly_a: 4, init_m: &[1, 3, 5, 13] },
    SobolDimSpec { degree: 5, poly_a: 2, init_m: &[1, 1, 3, 7, 9] },
    SobolDimSpec { degree: 5, poly_a: 4, init_m: &[1, 3, 5, 3, 29] },
    SobolDimSpec { degree: 5, poly_a: 7, init_m: &[1, 1, 3, 9, 1] },
    SobolDimSpec { degree: 5, poly_a: 11, init_m: &[1, 3, 1, 9, 19] },
    SobolDimSpec { degree: 5, poly_a: 13, init_m: &[1, 1, 7, 11, 1] },
    SobolDimSpec { degree: 5, poly_a: 14, init_m: &[1, 3, 5, 5, 31] },
    SobolDimSpec { degree: 6, poly_a: 1, init_m: &[1, 3, 3, 9, 7, 49] },
    SobolDimSpec { degree: 6, poly_a: 13, init_m: &[1, 1, 1, 15, 21, 21] },
    SobolDimSpec { degree: 6, poly_a: 16, init_m: &[1, 3, 1, 13, 27, 49] },
    SobolDimSpec { degree: 6, poly_a: 19, init_m: &[1, 1, 1, 15, 7, 5] },
    SobolDimSpec { degree: 6, poly_a: 22, init_m: &[1, 3, 1, 15, 13, 25] },
    SobolDimSpec { degree: 6, poly_a: 25, init_m: &[1, 1, 5, 5, 19, 61] },
    SobolDimSpec { degree: 7, poly_a: 1, init_m: &[1, 3, 7, 11, 23, 15, 103] },
    SobolDimSpec { degree: 7, poly_a: 4, init_m: &[1, 3, 7, 13, 13, 15, 69] },
];

fn direction_numbers_32(dim: usize) -> [u32; 32] {
    let mut v = [0u32; 32];
    if dim == 0 {
        for (i, slot) in v.iter_mut().enumerate() {
            *slot = 1u32 << (31 - i);
        }
        return v;
    }
    let spec = &SOBOL_TABLE[dim - 1];
    let s = spec.degree;
    for (i, &m_i) in spec.init_m.iter().enumerate().take(s) {
        v[i] = m_i << (31 - i);
    }
    for i in s..32 {
        let mut vi = v[i - s] ^ (v[i - s] >> s);
        for k in 1..s {
            if ((spec.poly_a >> (s - 1 - k)) & 1) != 0 {
                vi ^= v[i - k];
            }
        }
        v[i] = vi;
    }
    v
}

/// Generates `n_samples` unit-hypercube vectors in `(0, 1)^dim` via a low-discrepancy base-2 Sobol sequence
/// (`offset >= 1` skips the origin corner; when `scramble` is true, applies a deterministic `SplitMix64` XOR
/// digital shift per dimension that preserves binary `(t, m, s)` dyadic stratification).
pub fn generate_sobol_unit(
    n_samples: usize,
    dim: usize,
    seed: u64,
    scramble: bool,
    offset: usize,
) -> Result<Vec<Vec<f64>>, Error> {
    if n_samples == 0 || n_samples > MAX_QMC_SAMPLES {
        return Err(Error::invalid(format!(
            "Sobol QMC sampling requires 1 <= n_samples <= {MAX_QMC_SAMPLES}, got {n_samples}"
        )));
    }
    if dim == 0 || dim > MAX_SOBOL_DIM {
        return Err(Error::invalid(format!(
            "Sobol QMC sampling supports 1 <= dim <= {MAX_SOBOL_DIM} parameters, got {dim}"
        )));
    }

    let dirs: Vec<[u32; 32]> = (0..dim).map(direction_numbers_32).collect();
    let masks: Vec<u32> = if scramble {
        let mut rng = SplitMix64::new(seed);
        (0..dim).map(|_| (rng.next_u64() >> 32) as u32).collect()
    } else {
        vec![0u32; dim]
    };

    let scale = 1.0 / ((1u64 << 32) as f64);
    let mut table = Vec::with_capacity(n_samples);
    let mut x = vec![0u32; dim];

    // Advance Gray-code state cleanly up to `offset` so sequences always start reproducibly.
    for i in 1..offset {
        let c = (i - 1).trailing_ones() as usize;
        for d in 0..dim {
            x[d] ^= dirs[d][c.min(31)];
        }
    }

    for step in 0..n_samples {
        let seq_idx = offset + step;
        if seq_idx > 0 {
            let c = (seq_idx - 1).trailing_ones() as usize;
            for d in 0..dim {
                x[d] ^= dirs[d][c.min(31)];
            }
        }
        let row: Vec<f64> = (0..dim)
            .map(|d| {
                let bits = x[d] ^ masks[d];
                ((bits as f64) + 0.5) * scale
            })
            .collect();
        table.push(row);
    }

    Ok(table)
}

#[inline]
fn uniform_open01(rng: &mut SplitMix64) -> f64 {
    let bits = rng.next_u64() >> 11;
    ((bits as f64) + 0.5) * (1.0 / ((1u64 << 53) as f64))
}

/// Acklam minimax rational approximation to the inverse standard normal CDF `Phi^(-1)(p)` (`< 1.15e-9` max rel err).
pub fn inv_normal_cdf(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e+01,
        2.209_460_984_245_205e+02,
        -2.759_285_104_469_687e+02,
        1.383_577_518_672_69e+02,
        -3.066_479_806_614_716e+01,
        2.506_628_277_459_239e+00,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e+01,
        1.615_858_368_580_409e+02,
        -1.556_989_798_598_866e+02,
        6.680_131_188_771_972e+01,
        -1.328_068_155_288_572e+01,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-03,
        -3.223_964_580_411_365e-01,
        -2.400_758_277_161_838e+00,
        -2.549_732_539_343_734e+00,
        4.374_664_141_464_968e+00,
        2.938_163_982_698_783e+00,
    ];
    const D: [f64; 4] =
        [7.784_695_709_041_462e-03, 3.224_671_290_700_398e-01, 2.445_134_137_142_996e+00, 3.754_408_661_907_416e+00];
    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;

    let p = p.clamp(1.0e-15, 1.0 - 1.0e-15);
    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Standard normal CDF approximation (`erfc(-z / sqrt(2)) / 2`) helper for unit testing `inv_normal_cdf`.
#[cfg(test)]
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + libm_erf(z * FRAC_1_SQRT_2))
}

#[cfg(test)]
fn libm_erf(x: f64) -> f64 {
    let s = x.signum();
    let ax = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * ax);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t + 0.254_829_592)
            * t
            * (-ax * ax).exp();
    s * y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sobol_polynomials_are_primitive_and_equidistributed() {
        // 1. Verify GF(2) companion LFSR period == 2^s - 1 for every polynomial in SOBOL_TABLE
        for (idx, spec) in SOBOL_TABLE.iter().enumerate() {
            let s = spec.degree;
            let mask = (1u32 << s) - 1;
            // Feedback taps: x^s + sum_{k=1}^{s-1} a_k x^{s-k} + 1
            let mut full_taps = 1u32 << (s - 1); // constant term 1 shifted for Fibonacci LFSR
            for k in 1..s {
                if ((spec.poly_a >> (s - 1 - k)) & 1) != 0 {
                    full_taps |= 1u32 << (s - 1 - k);
                }
            }
            let _ = full_taps;
            // Galois LFSR representation: tap on polynomial x^s + c_{s-1} x^{s-1} + ... + c_0
            let poly = (1u32 << s) | (spec.poly_a << 1) | 1;
            let mut state = 1u32;
            let mut period = 0usize;
            for step in 1..=(1usize << s) {
                let lsb = state & 1;
                state >>= 1;
                if lsb != 0 {
                    state ^= (poly >> 1) & mask;
                }
                if state == 1 {
                    period = step;
                    break;
                }
            }
            assert_eq!(period, (1usize << s) - 1, "dim {} poly not primitive (s={s}, a={})", idx + 1, spec.poly_a);
        }

        // 2. N = 64 points (with offset = 0 so dyadic set 0..64 is complete): every 1D bin [k/64, (k+1)/64) has 1 point
        let pts = generate_sobol_unit(64, MAX_SOBOL_DIM, 0, false, 0).unwrap();
        for d in 0..MAX_SOBOL_DIM {
            let mut seen = [false; 64];
            for row in &pts {
                let bin = ((row[d] * 64.0).floor() as usize).min(63);
                assert!(!seen[bin], "1D Sobol collision in dim {d} bin {bin}");
                seen[bin] = true;
            }
        }

        // 3. Pairwise 2D stratification on 4x4 grid (N = 64 -> 4 points per cell!)
        for d1 in 0..8 {
            for d2 in (d1 + 1)..8 {
                let mut counts = [[0usize; 4]; 4];
                for row in &pts {
                    let b1 = ((row[d1] * 4.0).floor() as usize).min(3);
                    let b2 = ((row[d2] * 4.0).floor() as usize).min(3);
                    counts[b1][b2] += 1;
                }
                for c in counts.iter().flatten() {
                    assert_eq!(*c, 4, "2D 4x4 net violation between dims ({d1}, {d2})");
                }
            }
        }

        // 4. Scrambled Sobol preserves exact 1D stratification while varying with seed
        let scr_a = generate_sobol_unit(64, 5, 42, true, 0).unwrap();
        let scr_b = generate_sobol_unit(64, 5, 42, true, 0).unwrap();
        let scr_c = generate_sobol_unit(64, 5, 99, true, 0).unwrap();
        assert_eq!(scr_a, scr_b);
        assert_ne!(scr_a, scr_c);
        for d in 0..5 {
            let mut seen = [false; 64];
            for row in &scr_a {
                let bin = ((row[d] * 64.0).floor() as usize).min(63);
                assert!(!seen[bin]);
                seen[bin] = true;
            }
        }
    }

    #[test]
    fn lhs_is_deterministic_and_strictly_stratified() {
        let a = generate_lhs_unit(64, 6, 12345).unwrap();
        let b = generate_lhs_unit(64, 6, 12345).unwrap();
        let c = generate_lhs_unit(64, 6, 54321).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, c);
        for d in 0..6 {
            let mut seen = [false; 64];
            for row in &a {
                assert!(row[d] > 0.0 && row[d] < 1.0);
                let bin = ((row[d] * 64.0).floor() as usize).min(63);
                assert!(!seen[bin], "LHS stratum collision at dim {d} bin {bin}");
                seen[bin] = true;
            }
        }
    }

    #[test]
    fn inv_normal_cdf_and_margins_match_target_quantiles() {
        for &p in &[1.0e-5, 0.001, 0.025, 0.15865525, 0.5, 0.84134475, 0.975, 0.99999] {
            let z = inv_normal_cdf(p);
            let p_back = normal_cdf(z);
            assert!((p_back - p).abs() < 1.0e-6, "p={p}, z={z}, back={p_back}");
        }
        let u = MarginDistribution::Uniform { min: 2.0, max: 6.0 };
        assert!((u.map_unit(0.25) - 3.0).abs() < 1e-12);
        let n = MarginDistribution::Normal { mean: 10.0, sigma: 2.0, min: Some(7.0), max: None };
        assert!((n.map_unit(0.5) - 10.0).abs() < 1e-9);
        assert_eq!(n.map_unit(1.0e-6), 7.0);
    }
}
