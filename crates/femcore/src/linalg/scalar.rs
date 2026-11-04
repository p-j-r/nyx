//! The scalar abstraction: every solver in this crate is written once for `f64` and
//! `Complex<f64>`.
//!
//! Only these two types implement [`Scalar`] (the trait is sealed). The hot loops live in
//! [`Kernels`] so that a type can later get SIMD-specialised versions
//! without touching callers.

use crate::linalg::kernels::Kernels;
use num_complex::Complex;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// Complex double, the same type as `faer::c64`.
pub type C64 = Complex<f64>;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

/// A real or complex double-precision number.
pub trait Scalar:
    sealed::Sealed
    + Kernels
    + Copy
    + Debug
    + Default
    + PartialEq
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
{
    const IS_COMPLEX: bool;
    const ZERO: Self;
    const ONE: Self;

    fn from_f64(x: f64) -> Self;
    /// `re + i·im`; the imaginary part is ignored for real scalars.
    fn from_parts(re: f64, im: f64) -> Self;
    fn re(self) -> f64;
    fn im(self) -> f64;
    fn conj(self) -> Self;
    /// Modulus `|x|`.
    fn modulus(self) -> f64;
    /// `|x|²`.
    fn modulus_sq(self) -> f64;
    /// Multiplication by a real number.
    fn scale(self, r: f64) -> Self;
    fn is_finite(self) -> bool;
    /// Principal square root (`NaN` for a negative real number).
    fn sqrt_principal(self) -> Self;
    /// Widening to complex.
    fn to_c64(self) -> C64 {
        C64::new(self.re(), self.im())
    }
}

impl Scalar for f64 {
    const IS_COMPLEX: bool = false;
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    #[inline]
    fn from_f64(x: f64) -> Self {
        x
    }
    #[inline]
    fn from_parts(re: f64, _im: f64) -> Self {
        re
    }
    #[inline]
    fn re(self) -> f64 {
        self
    }
    #[inline]
    fn im(self) -> f64 {
        0.0
    }
    #[inline]
    fn conj(self) -> Self {
        self
    }
    #[inline]
    fn modulus(self) -> f64 {
        self.abs()
    }
    #[inline]
    fn modulus_sq(self) -> f64 {
        self * self
    }
    #[inline]
    fn scale(self, r: f64) -> Self {
        self * r
    }
    #[inline]
    fn sqrt_principal(self) -> Self {
        self.sqrt()
    }
    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

impl Scalar for C64 {
    const IS_COMPLEX: bool = true;
    const ZERO: Self = Complex::new(0.0, 0.0);
    const ONE: Self = Complex::new(1.0, 0.0);
    #[inline]
    fn from_f64(x: f64) -> Self {
        Complex::new(x, 0.0)
    }
    #[inline]
    fn from_parts(re: f64, im: f64) -> Self {
        Complex::new(re, im)
    }
    #[inline]
    fn re(self) -> f64 {
        self.re
    }
    #[inline]
    fn im(self) -> f64 {
        self.im
    }
    #[inline]
    fn conj(self) -> Self {
        Complex::conj(&self)
    }
    #[inline]
    fn modulus(self) -> f64 {
        self.norm()
    }
    #[inline]
    fn modulus_sq(self) -> f64 {
        self.norm_sqr()
    }
    #[inline]
    fn scale(self, r: f64) -> Self {
        Complex::new(self.re * r, self.im * r)
    }
    #[inline]
    fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    #[inline]
    fn sqrt_principal(self) -> Self {
        self.sqrt()
    }
}
