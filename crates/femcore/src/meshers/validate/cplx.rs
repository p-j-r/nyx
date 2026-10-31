//! Minimal complex number type (the harness needs nothing beyond the field operations).

use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };
pub const I: C64 = C64 { re: 0.0, im: 1.0 };

impl C64 {
    pub fn new(re: f64, im: f64) -> Self {
        C64 { re, im }
    }
    pub fn real(re: f64) -> Self {
        C64 { re, im: 0.0 }
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl Add for C64 {
    type Output = C64;
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C64 {
    type Output = C64;
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C64 {
    type Output = C64;
    fn mul(self, o: C64) -> C64 {
        C64::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
impl Mul<f64> for C64 {
    type Output = C64;
    fn mul(self, o: f64) -> C64 {
        C64::new(self.re * o, self.im * o)
    }
}
impl Div for C64 {
    type Output = C64;
    fn div(self, o: C64) -> C64 {
        let d = o.norm_sqr();
        C64::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}
impl Neg for C64 {
    type Output = C64;
    fn neg(self) -> C64 {
        C64::new(-self.re, -self.im)
    }
}
impl AddAssign for C64 {
    fn add_assign(&mut self, o: C64) {
        *self = *self + o;
    }
}
impl SubAssign for C64 {
    fn sub_assign(&mut self, o: C64) {
        *self = *self - o;
    }
}
impl MulAssign for C64 {
    fn mul_assign(&mut self, o: C64) {
        *self = *self * o;
    }
}
