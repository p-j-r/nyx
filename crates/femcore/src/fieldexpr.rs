//! A typed expression tree over complex field values, evaluated at one point at a time (productivity plan N1, ADR 90).
//!
//! The tree is what "handle arithmetic" builds: leaves are numbered slots (`Leaf(i)`, the caller says what slot `i` holds: a field of a
//! solution), nodes are the four operations, negation, integer powers and the functions `abs abs2 re im conj sqrt`. Nothing here knows
//! meshes, results or JSON; `nyx-engine::reduce` builds trees from requests and feeds the leaf values. Every operation is written out
//! (no `Complex::norm`, `Complex::powi`) so the order of the floating-point operations is fixed by this file and a hand-written loop can
//! reproduce the bits: `abs2(z) = re*re + im*im`, `z^n = z*z*...*z` left to right, `a*b` and `a/b` are the formulas below.
//!
//! The text formulas of N2 (`femcore::expr` grows to parse them) produce this same tree, so a formula and a handle expression are one
//! request and share one evaluator.

use crate::linalg::csr::C;

/// Binary operations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Add,
    Sub,
    Mul,
    Div,
}

/// Functions of one value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fun {
    /// the modulus, `sqrt(re*re + im*im)` (real)
    Abs,
    /// the squared modulus, `re*re + im*im` (real, no square root: the exact form of `z * conj(z)`)
    Abs2,
    Re,
    Im,
    Conj,
    /// the principal square root
    Sqrt,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Tree {
    Num(C),
    /// the value of leaf slot `i`
    Leaf(usize),
    Neg(Box<Tree>),
    Bin(Op, Box<Tree>, Box<Tree>),
    Call(Fun, Box<Tree>),
    /// an integer power, `n >= 1`, as repeated multiplication
    Pow(Box<Tree>, u32),
}

/// Largest integer power a tree may carry.
pub const MAX_POWER: u32 = 16;

fn mul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

fn div(a: C, b: C) -> C {
    let d = b.re * b.re + b.im * b.im;
    C::new((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d)
}

/// The principal square root: `sqrt((|z| + re) / 2)` for the real part with the sign of the imaginary part on the imaginary one.
fn sqrt(z: C) -> C {
    if z.re == 0.0 && z.im == 0.0 {
        return C::new(0.0, 0.0);
    }
    let m = (z.re * z.re + z.im * z.im).sqrt();
    let re = ((m + z.re) / 2.0).sqrt();
    let im = ((m - z.re) / 2.0).sqrt();
    C::new(re, if z.im < 0.0 { -im } else { im })
}

impl Tree {
    pub fn num(x: f64) -> Tree {
        Tree::Num(C::new(x, 0.0))
    }

    /// The number of leaf slots the tree reads (one more than the largest slot).
    pub fn leaves(&self) -> usize {
        match self {
            Tree::Num(_) => 0,
            Tree::Leaf(i) => i + 1,
            Tree::Neg(a) | Tree::Call(_, a) | Tree::Pow(a, _) => a.leaves(),
            Tree::Bin(_, a, b) => a.leaves().max(b.leaves()),
        }
    }

    /// The value at one point; `leaf[i]` holds slot `i`.
    pub fn eval(&self, leaf: &[C]) -> C {
        match self {
            Tree::Num(z) => *z,
            Tree::Leaf(i) => leaf[*i],
            Tree::Neg(a) => {
                let z = a.eval(leaf);
                C::new(-z.re, -z.im)
            }
            Tree::Bin(op, a, b) => {
                let (x, y) = (a.eval(leaf), b.eval(leaf));
                match op {
                    Op::Add => C::new(x.re + y.re, x.im + y.im),
                    Op::Sub => C::new(x.re - y.re, x.im - y.im),
                    Op::Mul => mul(x, y),
                    Op::Div => div(x, y),
                }
            }
            Tree::Call(f, a) => {
                let z = a.eval(leaf);
                match f {
                    Fun::Abs => C::new((z.re * z.re + z.im * z.im).sqrt(), 0.0),
                    Fun::Abs2 => C::new(z.re * z.re + z.im * z.im, 0.0),
                    Fun::Re => C::new(z.re, 0.0),
                    Fun::Im => C::new(z.im, 0.0),
                    Fun::Conj => C::new(z.re, -z.im),
                    Fun::Sqrt => sqrt(z),
                }
            }
            Tree::Pow(a, n) => {
                let z = a.eval(leaf);
                let mut r = z;
                for _ in 1..*n {
                    r = mul(r, z);
                }
                r
            }
        }
    }

    /// Whether the value is real for real leaves (no leaf is complex): used to report a real answer as a number, not a pair.
    /// `true` when the tree ends in `abs`, `abs2`, `re` or `im`, or reads no imaginary part.
    pub fn is_real_valued(&self, leaf_is_complex: &[bool]) -> bool {
        match self {
            Tree::Num(z) => z.im == 0.0,
            Tree::Leaf(i) => !leaf_is_complex[*i],
            Tree::Neg(a) | Tree::Pow(a, _) => a.is_real_valued(leaf_is_complex),
            Tree::Bin(_, a, b) => a.is_real_valued(leaf_is_complex) && b.is_real_valued(leaf_is_complex),
            Tree::Call(Fun::Abs | Fun::Abs2 | Fun::Re | Fun::Im, _) => true,
            Tree::Call(Fun::Conj, a) => a.is_real_valued(leaf_is_complex),
            // the root of a negative real is imaginary
            Tree::Call(Fun::Sqrt, _) => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn l(i: usize) -> Box<Tree> {
        Box::new(Tree::Leaf(i))
    }

    #[test]
    fn operations_follow_their_written_formulas() {
        let (a, b) = (C::new(1.5, -2.0), C::new(0.25, 3.0));
        let leaf = [a, b];
        assert_eq!(Tree::Bin(Op::Mul, l(0), l(1)).eval(&leaf), C::new(1.5 * 0.25 - (-2.0) * 3.0, 1.5 * 3.0 + (-2.0) * 0.25));
        assert_eq!(Tree::Call(Fun::Abs2, l(0)).eval(&leaf), C::new(1.5 * 1.5 + 4.0, 0.0));
        // |z|^2 as z * conj(z) is the same number, bit for bit, with a zero imaginary part
        let zz = Tree::Bin(Op::Mul, l(0), Box::new(Tree::Call(Fun::Conj, l(0)))).eval(&leaf);
        assert_eq!(zz.re, Tree::Call(Fun::Abs2, l(0)).eval(&leaf).re);
        assert_eq!(zz.im, 0.0);
        // a power is repeated multiplication
        assert_eq!(Tree::Pow(l(0), 3).eval(&leaf), mul(mul(a, a), a));
        // the quotient undoes the product
        let q = Tree::Bin(Op::Div, Box::new(Tree::Bin(Op::Mul, l(0), l(1))), l(1)).eval(&leaf);
        assert!((q - a).norm() < 1e-15);
        let s = Tree::Call(Fun::Sqrt, Box::new(Tree::Num(C::new(-4.0, 0.0)))).eval(&leaf);
        assert_eq!(s, C::new(0.0, 2.0));
        assert_eq!(Tree::Bin(Op::Add, l(0), Box::new(Tree::num(1.0))).leaves(), 1);
    }

    #[test]
    fn realness_is_known_from_the_tree() {
        let t = Tree::Call(Fun::Abs2, l(0));
        assert!(t.is_real_valued(&[true]));
        assert!(!Tree::Leaf(0).is_real_valued(&[true]));
        assert!(Tree::Bin(Op::Add, l(0), Box::new(Tree::num(2.0))).is_real_valued(&[false]));
    }
}
