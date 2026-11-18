//! A typed expression tree over complex field values, evaluated at one point at a time (productivity plan N1 and N2, ADR 90).
//!
//! The tree is what "handle arithmetic" and the text formulas of `femcore::formula` build: leaves are numbered slots (`Leaf(i)`, the caller
//! says what slot `i` holds: a field of a solution, a coordinate, a material value), nodes are the four operations, negation, integer
//! powers, the comparisons, `Where`, and the functions below. Nothing here knows meshes, results or JSON; `nyx-engine::reduce` builds
//! trees from requests and feeds the leaf values. Every operation is written out (no `Complex::norm`, `Complex::powi`) so the order of the
//! floating-point operations is fixed by this file and a hand-written loop can reproduce the bits: `abs2(z) = re*re + im*im`,
//! `z^n = z*z*...*z` left to right, `a*b` and `a/b` are the formulas below.
//!
//! Complex values flow through `+ - * /`, integer powers, `abs abs2 re im conj sqrt exp sin cos`. What has no meaning for a complex value
//! (an order, `log`, `tan`, a non-integer power, a comparison) takes real operands: [`Tree::check`] refuses a complex operand before any
//! point is evaluated, with a message that says what to wrap it in. A comparison and a logical operation answer 1 or 0 (real).

use crate::linalg::csr::C;

/// Binary operations. `And` and `Or` read the truth of the real parts (nonzero is true) and answer 1 or 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Add,
    Sub,
    Mul,
    Div,
    And,
    Or,
}

/// Comparisons of two real values; the answer is 1 or 0 (false when either side is NaN).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
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
    /// `exp(re) (cos im + i sin im)`
    Exp,
    /// `sin(re) cosh(im) + i cos(re) sinh(im)`
    Sin,
    /// `cos(re) cosh(im) - i sin(re) sinh(im)`
    Cos,
    /// the argument (phase) of a complex value, `atan2(im, re)` (real, -pi to pi)
    Arg,
    /// the square root of a real value that is not negative (real; NaN below zero): the root of a norm
    RealSqrt,
    /// natural logarithm of a real value (NaN below zero)
    Log,
    Log10,
    Tan,
    Sinh,
    Cosh,
    Tanh,
    Atan,
    /// 1 for a value of at least 0, else 0
    Step,
    /// 1 for 0, else 0
    Not,
}

impl Fun {
    /// The name a formula uses.
    pub fn name(self) -> &'static str {
        match self {
            Fun::Abs => "abs",
            Fun::Abs2 => "abs2",
            Fun::Re => "re",
            Fun::Im => "im",
            Fun::Conj => "conj",
            Fun::Sqrt => "sqrt",
            Fun::Exp => "exp",
            Fun::Sin => "sin",
            Fun::Cos => "cos",
            Fun::Arg => "arg",
            Fun::RealSqrt => "norm",
            Fun::Log => "log",
            Fun::Log10 => "log10",
            Fun::Tan => "tan",
            Fun::Sinh => "sinh",
            Fun::Cosh => "cosh",
            Fun::Tanh => "tanh",
            Fun::Atan => "atan",
            Fun::Step => "step",
            Fun::Not => "not",
        }
    }

    /// Whether the function is defined on a complex value (the others read a real one).
    fn takes_complex(self) -> bool {
        matches!(
            self,
            Fun::Abs | Fun::Abs2 | Fun::Re | Fun::Im | Fun::Conj | Fun::Sqrt | Fun::Exp | Fun::Sin | Fun::Cos | Fun::Arg
        )
    }
}

/// Functions of two real values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fun2 {
    Min,
    Max,
    /// `a.powf(b)`, for a power that is not a small whole number
    Pow,
    /// `atan2(y, x)`
    Atan2,
}

impl Fun2 {
    pub fn name(self) -> &'static str {
        match self {
            Fun2::Min => "min",
            Fun2::Max => "max",
            Fun2::Pow => "^",
            Fun2::Atan2 => "atan2",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Tree {
    Num(C),
    /// the value of leaf slot `i`
    Leaf(usize),
    Neg(Box<Tree>),
    Bin(Op, Box<Tree>, Box<Tree>),
    Call(Fun, Box<Tree>),
    Call2(Fun2, Box<Tree>, Box<Tree>),
    /// an integer power, `n >= 1`, as repeated multiplication
    Pow(Box<Tree>, u32),
    Cmp(Cmp, Box<Tree>, Box<Tree>),
    /// `Where(c, a, b)`: `a` where the real part of `c` is not 0, else `b`; only the branch taken is evaluated
    Where(Box<Tree>, Box<Tree>, Box<Tree>),
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

fn real(x: f64) -> C {
    C::new(x, 0.0)
}

fn truth(b: bool) -> C {
    real(if b { 1.0 } else { 0.0 })
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
            Tree::Bin(_, a, b) | Tree::Call2(_, a, b) | Tree::Cmp(_, a, b) => a.leaves().max(b.leaves()),
            Tree::Where(c, a, b) => c.leaves().max(a.leaves()).max(b.leaves()),
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
                    Op::And => truth(x.re != 0.0 && y.re != 0.0),
                    Op::Or => truth(x.re != 0.0 || y.re != 0.0),
                }
            }
            Tree::Call(f, a) => {
                let z = a.eval(leaf);
                match f {
                    Fun::Abs => real((z.re * z.re + z.im * z.im).sqrt()),
                    Fun::Abs2 => real(z.re * z.re + z.im * z.im),
                    Fun::Re => real(z.re),
                    Fun::Im => real(z.im),
                    Fun::Conj => C::new(z.re, -z.im),
                    Fun::Sqrt => sqrt(z),
                    Fun::Exp => {
                        let m = z.re.exp();
                        C::new(m * z.im.cos(), m * z.im.sin())
                    }
                    Fun::Sin => C::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
                    Fun::Cos => C::new(z.re.cos() * z.im.cosh(), -(z.re.sin() * z.im.sinh())),
                    Fun::Arg => real(z.im.atan2(z.re)),
                    Fun::RealSqrt => real(z.re.sqrt()),
                    Fun::Log => real(z.re.ln()),
                    Fun::Log10 => real(z.re.log10()),
                    Fun::Tan => real(z.re.tan()),
                    Fun::Sinh => real(z.re.sinh()),
                    Fun::Cosh => real(z.re.cosh()),
                    Fun::Tanh => real(z.re.tanh()),
                    Fun::Atan => real(z.re.atan()),
                    Fun::Step => truth(z.re >= 0.0),
                    Fun::Not => truth(z.re == 0.0),
                }
            }
            Tree::Call2(f, a, b) => {
                let (x, y) = (a.eval(leaf).re, b.eval(leaf).re);
                real(match f {
                    Fun2::Min => x.min(y),
                    Fun2::Max => x.max(y),
                    Fun2::Pow => x.powf(y),
                    Fun2::Atan2 => x.atan2(y),
                })
            }
            Tree::Pow(a, n) => {
                let z = a.eval(leaf);
                let mut r = z;
                for _ in 1..*n {
                    r = mul(r, z);
                }
                r
            }
            Tree::Cmp(c, a, b) => {
                let (x, y) = (a.eval(leaf).re, b.eval(leaf).re);
                truth(match c {
                    Cmp::Lt => x < y,
                    Cmp::Le => x <= y,
                    Cmp::Gt => x > y,
                    Cmp::Ge => x >= y,
                    Cmp::Eq => x == y,
                    Cmp::Ne => x != y,
                })
            }
            Tree::Where(c, a, b) => {
                if c.eval(leaf).re != 0.0 {
                    a.eval(leaf)
                } else {
                    b.eval(leaf)
                }
            }
        }
    }

    /// Whether the value is real for these leaves: `true` when the tree ends in a function that answers a real value, or reads no imaginary
    /// part. Used to report a real answer as a number, not a pair.
    pub fn is_real_valued(&self, leaf_is_complex: &[bool]) -> bool {
        match self {
            Tree::Num(z) => z.im == 0.0,
            Tree::Leaf(i) => !leaf_is_complex[*i],
            Tree::Neg(a) | Tree::Pow(a, _) => a.is_real_valued(leaf_is_complex),
            Tree::Bin(Op::And | Op::Or, _, _) | Tree::Cmp(..) | Tree::Call2(..) => true,
            Tree::Bin(_, a, b) => a.is_real_valued(leaf_is_complex) && b.is_real_valued(leaf_is_complex),
            Tree::Call(Fun::Conj | Fun::Exp | Fun::Sin | Fun::Cos, a) => a.is_real_valued(leaf_is_complex),
            // the root of a negative real is imaginary
            Tree::Call(Fun::Sqrt, _) => false,
            Tree::Call(..) => true,
            Tree::Where(_, a, b) => a.is_real_valued(leaf_is_complex) && b.is_real_valued(leaf_is_complex),
        }
    }

    /// Refuses what has no meaning for the leaves' values: an order, a comparison, a logical operation, `log`, `tan`, a non-integer power
    /// of a complex value. The message names the operation and says what to wrap the value in. Run it once per request, before any point.
    pub fn check(&self, leaf_is_complex: &[bool]) -> Result<(), String> {
        let need_real = |what: &str, t: &Tree| -> Result<(), String> {
            if t.is_real_valued(leaf_is_complex) {
                Ok(())
            } else {
                Err(format!("`{what}` needs a real value, and this one is complex: wrap it in abs, abs2, re or im"))
            }
        };
        match self {
            Tree::Num(_) | Tree::Leaf(_) => Ok(()),
            Tree::Neg(a) | Tree::Pow(a, _) => a.check(leaf_is_complex),
            Tree::Bin(op, a, b) => {
                if matches!(op, Op::And | Op::Or) {
                    let w = if *op == Op::And { "and" } else { "or" };
                    need_real(w, a)?;
                    need_real(w, b)?;
                }
                a.check(leaf_is_complex)?;
                b.check(leaf_is_complex)
            }
            Tree::Call(f, a) => {
                if !f.takes_complex() {
                    need_real(f.name(), a)?;
                }
                a.check(leaf_is_complex)
            }
            Tree::Call2(f, a, b) => {
                need_real(f.name(), a)?;
                need_real(f.name(), b)?;
                a.check(leaf_is_complex)?;
                b.check(leaf_is_complex)
            }
            Tree::Cmp(_, a, b) => {
                need_real("a comparison", a)?;
                need_real("a comparison", b)?;
                a.check(leaf_is_complex)?;
                b.check(leaf_is_complex)
            }
            Tree::Where(c, a, b) => {
                need_real("where", c)?;
                c.check(leaf_is_complex)?;
                a.check(leaf_is_complex)?;
                b.check(leaf_is_complex)
            }
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
        assert!(Tree::Cmp(Cmp::Lt, l(0), l(0)).is_real_valued(&[true]));
        assert!(!Tree::Call(Fun::Exp, l(0)).is_real_valued(&[true]));
        assert!(Tree::Call(Fun::Exp, l(0)).is_real_valued(&[false]));
    }

    #[test]
    fn complex_functions_reduce_to_the_real_ones_on_real_values() {
        let x = C::new(0.7, 0.0);
        let leaf = [x];
        assert_eq!(Tree::Call(Fun::Exp, l(0)).eval(&leaf), C::new(0.7_f64.exp(), 0.0));
        assert_eq!(Tree::Call(Fun::Sin, l(0)).eval(&leaf), C::new(0.7_f64.sin(), 0.0));
        assert_eq!(Tree::Call(Fun::Cos, l(0)).eval(&leaf), C::new(0.7_f64.cos(), 0.0));
        // e^(i pi) = -1
        let z = Tree::Call(Fun::Exp, Box::new(Tree::Num(C::new(0.0, std::f64::consts::PI)))).eval(&leaf);
        assert!((z.re + 1.0).abs() < 1e-15 && z.im.abs() < 1e-15);
        // sin(i) = i sinh(1)
        let s = Tree::Call(Fun::Sin, Box::new(Tree::Num(C::new(0.0, 1.0)))).eval(&leaf);
        assert!(s.re.abs() < 1e-15 && (s.im - 1.0_f64.sinh()).abs() < 1e-15);
        assert_eq!(Tree::Call(Fun::Arg, Box::new(Tree::Num(C::new(0.0, 2.0)))).eval(&leaf).re, std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn where_takes_one_branch_and_comparisons_answer_one_or_zero() {
        let leaf = [C::new(3.0, 0.0)];
        let cond = Tree::Cmp(Cmp::Gt, l(0), Box::new(Tree::num(2.0)));
        assert_eq!(cond.eval(&leaf), C::new(1.0, 0.0));
        let w = Tree::Where(Box::new(cond), Box::new(Tree::num(10.0)), Box::new(Tree::num(20.0)));
        assert_eq!(w.eval(&leaf), C::new(10.0, 0.0));
        assert_eq!(w.eval(&[C::new(1.0, 0.0)]), C::new(20.0, 0.0));
        // the branch not taken is not evaluated: a division by zero there leaves no NaN
        let guarded = Tree::Where(
            Box::new(Tree::Cmp(Cmp::Ne, l(0), Box::new(Tree::num(0.0)))),
            Box::new(Tree::Bin(Op::Div, Box::new(Tree::num(1.0)), l(0))),
            Box::new(Tree::num(0.0)),
        );
        assert_eq!(guarded.eval(&[C::new(0.0, 0.0)]), C::new(0.0, 0.0));
        assert_eq!(guarded.eval(&[C::new(4.0, 0.0)]), C::new(0.25, 0.0));
        // NaN compares false
        assert_eq!(Tree::Cmp(Cmp::Eq, l(0), l(0)).eval(&[C::new(f64::NAN, 0.0)]), C::new(0.0, 0.0));
    }

    #[test]
    fn check_refuses_a_complex_operand_where_an_order_is_needed() {
        let t = Tree::Cmp(Cmp::Gt, l(0), Box::new(Tree::num(0.0)));
        assert!(t.check(&[false]).is_ok());
        let e = t.check(&[true]).unwrap_err();
        assert!(e.contains("complex") && e.contains("abs"), "{e}");
        let wrapped = Tree::Cmp(Cmp::Gt, Box::new(Tree::Call(Fun::Abs, l(0))), Box::new(Tree::num(0.0)));
        assert!(wrapped.check(&[true]).is_ok());
        assert!(Tree::Call(Fun::Log, l(0)).check(&[true]).is_err());
        assert!(Tree::Call(Fun::Exp, l(0)).check(&[true]).is_ok());
        assert!(Tree::Call2(Fun2::Max, l(0), l(0)).check(&[true]).is_err());
    }
}
