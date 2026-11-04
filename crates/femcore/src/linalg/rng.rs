//! Small deterministic generator for start vectors and tests (SplitMix64).

use crate::linalg::scalar::Scalar;

#[derive(Debug, Clone)]
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in `[-1, 1)`.
    pub fn next_signed(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
}

/// Deterministic pseudo-random vector with entries in `[-1,1)` (real and imaginary parts for
/// complex scalars).
pub fn random_vector<T: Scalar>(n: usize, seed: u64) -> Vec<T> {
    let mut r = SplitMix64::new(seed);
    (0..n)
        .map(|_| {
            let re = r.next_signed();
            let im = if T::IS_COMPLEX { r.next_signed() } else { 0.0 };
            T::from_parts(re, im)
        })
        .collect()
}
