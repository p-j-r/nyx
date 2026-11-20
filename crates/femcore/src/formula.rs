// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Text formulas over fields, coordinates and material values (productivity plan N2, ADR 90): `abs(Ex)^2 * eps_r`, `where(x > 0, T, 0)`,
//! `dot(E, conj(E))`. The text is parsed once into a `fieldexpr::Tree`, the same tree that handle arithmetic builds, so a formula and a
//! handle expression are one request with one evaluator. The grammar allows nothing but arithmetic: no input can run code, name a file or
//! call anything that is not listed here.
//!
//! * **Numbers** `1`, `.5`, `1.5e-7`; a trailing `j` makes an imaginary number (`2j`). Constants `pi` and `e`.
//! * **Operators** `+ - * /`, `^` or `**` (right to left; `-2^2` is -4), the comparisons `< <= > >= == !=`, and `and`, `or`, `not`
//!   (comparisons answer 1 or 0, any nonzero real is true). A power with a whole constant exponent from -16 to 16 is repeated multiplication
//!   (complex values allowed); any other power reads real values.
//! * **Functions of one value** `abs abs2 re im conj sqrt exp sin cos arg` (complex values allowed) and `log ln log10 tan sinh cosh tanh atan
//!   step` (real values); of two `min max atan2` (real); of three `where(cond, a, b)`: only the branch that applies is evaluated, so
//!   `where(x != 0, 1/x, 0)` never divides by zero.
//! * **Names** come from the caller ([`Names`]): coordinates `x y z`, fields, material values. A name that is not a plain identifier is written
//!   in backticks (`` `|E|` ``, `` `Re Ex` ``). `material("rho")` and `region("Si")` say what they are, for a name a field also has.
//! * **Vectors** are the components a name stands for (`E` is `Ex Ey Ez` when a solution has them) or a list `[a, b, c]`. They add and
//!   subtract, multiply and divide by a scalar, and give `dot(a, b)` (no conjugate), `inner(a, b)` (conjugates `a`), `cross(a, b)`,
//!   `norm2(v)` (the sum of the squared moduli) and `norm(v)` (its square root); `v.x`, `v.y`, `v.z` take a component. A formula ends in a
//!   scalar.
//!
//! `femcore::expr` stays what it is: the real formula grammar of geometry and initial fields that `geom::expr` copies (`tests/expr_grammar.rs`).

use crate::fieldexpr::{Cmp, Fun, Fun2, Op, Tree, MAX_POWER};
use crate::linalg::csr::C;

/// What a name stands for.
pub enum Named {
    Scalar(Tree),
    /// the components, each with its axis (0 x, 1 y, 2 z)
    Vector(Vec<(usize, Tree)>),
}

/// The caller's side of a formula: what its names mean.
pub trait Names {
    /// `None` when the name is not known (the error then lists [`Names::known`]).
    fn name(&mut self, name: &str) -> Result<Option<Named>, String>;
    /// `material("name")`: a value of the material of each cell.
    fn material(&mut self, name: &str) -> Result<Tree, String>;
    /// `region("name")`: 1 in the cells of the region, else 0.
    fn region(&mut self, name: &str) -> Result<Tree, String>;
    /// A short list of the names the caller has, for an error message.
    fn known(&self) -> String {
        String::new()
    }
}

enum V {
    S(Tree),
    Vec(Vec<(usize, Tree)>),
}

const MAX_DEPTH: usize = 64;
const AXES: [&str; 3] = ["x", "y", "z"];

struct P<'a, 'n> {
    s: &'a [u8],
    i: usize,
    depth: usize,
    names: &'n mut dyn Names,
}

fn b(t: Tree) -> Box<Tree> {
    Box::new(t)
}

fn constant(t: &Tree) -> Option<f64> {
    match t {
        Tree::Num(z) if z.im == 0.0 => Some(z.re),
        Tree::Neg(a) => constant(a).map(|x| -x),
        _ => None,
    }
}

fn axes_of(v: &[(usize, Tree)]) -> String {
    v.iter().map(|(a, _)| AXES[*a]).collect::<Vec<_>>().join(", ")
}

/// Parses `text` into a scalar tree; `names` says what its names are.
pub fn parse(text: &str, names: &mut dyn Names) -> Result<Tree, String> {
    if text.trim().is_empty() {
        return Err("empty formula".into());
    }
    let mut p = P { s: text.as_bytes(), i: 0, depth: 0, names };
    let v = p.or()?;
    p.ws();
    if p.i != p.s.len() {
        return p.fail(format!("unexpected '{}'", p.rest_char()));
    }
    match v {
        V::S(t) => Ok(t),
        V::Vec(c) => Err(format!(
            "the formula is a vector (components {}); end it in a scalar: a component (`.x`), `norm(v)`, `norm2(v)` or `dot(a, b)`",
            axes_of(&c)
        )),
    }
}

impl<'a, 'n> P<'a, 'n> {
    fn fail<T>(&self, msg: impl Into<String>) -> Result<T, String> {
        Err(format!("{} (at character {})", msg.into(), self.i + 1))
    }

    fn rest_char(&self) -> char {
        String::from_utf8_lossy(&self.s[self.i..]).chars().next().unwrap_or('?')
    }

    fn ws(&mut self) {
        while self.i < self.s.len() && (self.s[self.i] as char).is_whitespace() {
            self.i += 1;
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.s.get(self.i).copied()
    }

    fn at(&self, off: usize) -> Option<u8> {
        self.s.get(self.i + off).copied()
    }

    /// The identifier at the cursor, without consuming it.
    fn peek_word(&mut self) -> Option<&'a str> {
        self.ws();
        let s: &'a [u8] = self.s;
        let c = *s.get(self.i)?;
        if !(c.is_ascii_alphabetic() || c == b'_') {
            return None;
        }
        let mut j = self.i;
        while j < s.len() && (s[j].is_ascii_alphanumeric() || s[j] == b'_') {
            j += 1;
        }
        std::str::from_utf8(&s[self.i..j]).ok()
    }

    fn word(&mut self, w: &str) -> bool {
        if self.peek_word() == Some(w) {
            self.i += w.len();
            true
        } else {
            false
        }
    }

    fn scalar(&self, v: V, what: &str) -> Result<Tree, String> {
        match v {
            V::S(t) => Ok(t),
            V::Vec(c) => self.fail(format!("{what} needs a scalar, and this is a vector (components {})", axes_of(&c))),
        }
    }

    fn enter(&mut self) -> Result<(), String> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return self.fail("the formula is nested too deeply");
        }
        Ok(())
    }

    // ---------------------------------------------------------------------------------------------------------- the grammar

    fn or(&mut self) -> Result<V, String> {
        let mut l = self.and()?;
        while self.word("or") {
            let r = self.and()?;
            let (l_, r_) = (self.scalar(l, "`or`")?, self.scalar(r, "`or`")?);
            l = V::S(Tree::Bin(Op::Or, b(l_), b(r_)));
        }
        Ok(l)
    }

    fn and(&mut self) -> Result<V, String> {
        let mut l = self.not()?;
        while self.word("and") {
            let r = self.not()?;
            let (l_, r_) = (self.scalar(l, "`and`")?, self.scalar(r, "`and`")?);
            l = V::S(Tree::Bin(Op::And, b(l_), b(r_)));
        }
        Ok(l)
    }

    fn not(&mut self) -> Result<V, String> {
        if self.word("not") {
            self.enter()?;
            let v = self.not()?;
            self.depth -= 1;
            let t = self.scalar(v, "`not`")?;
            return Ok(V::S(Tree::Call(Fun::Not, b(t))));
        }
        self.cmp()
    }

    fn cmp_op(&mut self) -> Option<Cmp> {
        self.ws();
        let (c, d) = (self.at(0)?, self.at(1));
        let (op, len) = match (c, d) {
            (b'<', Some(b'=')) => (Cmp::Le, 2),
            (b'>', Some(b'=')) => (Cmp::Ge, 2),
            (b'=', Some(b'=')) => (Cmp::Eq, 2),
            (b'!', Some(b'=')) => (Cmp::Ne, 2),
            (b'<', _) => (Cmp::Lt, 1),
            (b'>', _) => (Cmp::Gt, 1),
            _ => return None,
        };
        self.i += len;
        Some(op)
    }

    fn cmp(&mut self) -> Result<V, String> {
        let l = self.sum()?;
        let Some(op) = self.cmp_op() else {
            if self.peek() == Some(b'=') {
                return self.fail("a single '=' assigns; compare with '=='");
            }
            return Ok(l);
        };
        let r = self.sum()?;
        let (l, r) = (self.scalar(l, "a comparison")?, self.scalar(r, "a comparison")?);
        if self.cmp_op().is_some() {
            return self.fail("a comparison does not chain: join two with `and`");
        }
        Ok(V::S(Tree::Cmp(op, b(l), b(r))))
    }

    fn sum(&mut self) -> Result<V, String> {
        let mut l = self.term()?;
        while let Some(c) = self.peek() {
            if c != b'+' && c != b'-' {
                break;
            }
            self.i += 1;
            let r = self.term()?;
            l = self.addsub(if c == b'+' { Op::Add } else { Op::Sub }, l, r)?;
        }
        Ok(l)
    }

    fn addsub(&self, op: Op, l: V, r: V) -> Result<V, String> {
        match (l, r) {
            (V::S(x), V::S(y)) => Ok(V::S(Tree::Bin(op, b(x), b(y)))),
            (V::Vec(x), V::Vec(y)) => {
                if x.len() != y.len() || x.iter().zip(&y).any(|(p, q)| p.0 != q.0) {
                    return self.fail(format!(
                        "vectors with components {} and {} do not add; they need the same components",
                        axes_of(&x),
                        axes_of(&y)
                    ));
                }
                Ok(V::Vec(x.into_iter().zip(y).map(|((a, p), (_, q))| (a, Tree::Bin(op, b(p), b(q)))).collect()))
            }
            _ => self.fail("a vector and a scalar do not add"),
        }
    }

    fn term(&mut self) -> Result<V, String> {
        let mut l = self.unary()?;
        while let Some(c) = self.peek() {
            if !(c == b'/' || c == b'*' && self.at(1) != Some(b'*')) {
                break;
            }
            self.i += 1;
            let r = self.unary()?;
            l = match (c, l, r) {
                (b'*', V::S(x), V::S(y)) => V::S(Tree::Bin(Op::Mul, b(x), b(y))),
                (b'/', V::S(x), V::S(y)) => V::S(Tree::Bin(Op::Div, b(x), b(y))),
                (b'*', V::S(s), V::Vec(v)) => {
                    V::Vec(v.into_iter().map(|(a, t)| (a, Tree::Bin(Op::Mul, b(s.clone()), b(t)))).collect())
                }
                (b'*', V::Vec(v), V::S(s)) => {
                    V::Vec(v.into_iter().map(|(a, t)| (a, Tree::Bin(Op::Mul, b(t), b(s.clone())))).collect())
                }
                (b'/', V::Vec(v), V::S(s)) => {
                    V::Vec(v.into_iter().map(|(a, t)| (a, Tree::Bin(Op::Div, b(t), b(s.clone())))).collect())
                }
                (b'*', V::Vec(_), V::Vec(_)) => {
                    return self.fail("two vectors do not multiply with '*': use dot(a, b), inner(a, b) or cross(a, b)")
                }
                _ => return self.fail("a scalar does not divide by a vector"),
            };
        }
        Ok(l)
    }

    fn unary(&mut self) -> Result<V, String> {
        match self.peek() {
            Some(b'-') => {
                self.i += 1;
                self.enter()?;
                let v = self.unary()?;
                self.depth -= 1;
                Ok(match v {
                    V::S(t) => V::S(Tree::Neg(b(t))),
                    V::Vec(c) => V::Vec(c.into_iter().map(|(a, t)| (a, Tree::Neg(b(t)))).collect()),
                })
            }
            Some(b'+') => {
                self.i += 1;
                self.enter()?;
                let v = self.unary()?;
                self.depth -= 1;
                Ok(v)
            }
            _ => self.power(),
        }
    }

    fn power(&mut self) -> Result<V, String> {
        let base = self.postfix()?;
        let star2 = self.peek() == Some(b'*') && self.at(1) == Some(b'*');
        if self.peek() == Some(b'^') || star2 {
            self.i += if star2 { 2 } else { 1 };
            self.enter()?;
            let e = self.unary()?;
            self.depth -= 1;
            let (base, e) = (self.scalar(base, "a power")?, self.scalar(e, "an exponent")?);
            return Ok(V::S(self.pow(base, e)));
        }
        Ok(base)
    }

    /// `base ^ e`: a whole constant exponent is repeated multiplication (and 1 / that for a negative one), any other power reads real values.
    fn pow(&self, base: Tree, e: Tree) -> Tree {
        if let Some(n) = constant(&e) {
            if n.fract() == 0.0 && n.abs() <= MAX_POWER as f64 {
                let k = n as i64;
                return match k {
                    0 => Tree::num(1.0),
                    1 => base,
                    k if k > 1 => Tree::Pow(b(base), k as u32),
                    k => Tree::Bin(Op::Div, b(Tree::num(1.0)), b(Tree::Pow(b(base), (-k) as u32))),
                };
            }
        }
        Tree::Call2(Fun2::Pow, b(base), b(e))
    }

    fn postfix(&mut self) -> Result<V, String> {
        let mut v = self.atom()?;
        while self.peek() == Some(b'.') && self.at(1).is_some_and(|c| c.is_ascii_alphabetic()) {
            self.i += 1;
            let name = self.peek_word().unwrap_or("");
            let axis = AXES.iter().position(|a| *a == name);
            let (Some(axis), V::Vec(c)) = (axis, &v) else {
                return self.fail(format!("'.{name}' takes a component x, y or z of a vector"));
            };
            self.i += name.len();
            let Some((_, t)) = c.iter().find(|(a, _)| *a == axis) else {
                return self.fail(format!("this vector has the components {}, not {}", axes_of(c), AXES[axis]));
            };
            v = V::S(t.clone());
        }
        Ok(v)
    }

    fn number(&mut self) -> Result<V, String> {
        let st = self.i;
        while self.i < self.s.len() && (self.s[self.i].is_ascii_digit() || self.s[self.i] == b'.') {
            self.i += 1;
        }
        if self.i < self.s.len() && (self.s[self.i] == b'e' || self.s[self.i] == b'E') {
            let save = self.i;
            self.i += 1;
            if self.i < self.s.len() && (self.s[self.i] == b'+' || self.s[self.i] == b'-') {
                self.i += 1;
            }
            if self.i < self.s.len() && self.s[self.i].is_ascii_digit() {
                while self.i < self.s.len() && self.s[self.i].is_ascii_digit() {
                    self.i += 1;
                }
            } else {
                self.i = save;
            }
        }
        let t = std::str::from_utf8(&self.s[st..self.i]).unwrap_or("");
        let x: f64 = match t.parse() {
            Ok(x) => x,
            Err(_) => return self.fail(format!("bad number '{t}'")),
        };
        if self.at(0) == Some(b'j') && !self.at(1).is_some_and(|c| c.is_ascii_alphanumeric() || c == b'_') {
            self.i += 1;
            return Ok(V::S(Tree::Num(C::new(0.0, x))));
        }
        Ok(V::S(Tree::num(x)))
    }

    fn string(&mut self) -> Result<String, String> {
        let q = match self.peek() {
            Some(q @ (b'"' | b'\'')) => q,
            _ => return self.fail("expected a quoted name"),
        };
        self.i += 1;
        let st = self.i;
        while self.i < self.s.len() && self.s[self.i] != q {
            self.i += 1;
        }
        if self.i >= self.s.len() {
            return self.fail("a quoted name is not closed");
        }
        let t = String::from_utf8_lossy(&self.s[st..self.i]).into_owned();
        self.i += 1;
        Ok(t)
    }

    fn named(&mut self, name: &str) -> Result<V, String> {
        match name {
            "pi" => return Ok(V::S(Tree::num(std::f64::consts::PI))),
            "e" => return Ok(V::S(Tree::num(std::f64::consts::E))),
            _ => {}
        }
        match self.names.name(name) {
            Ok(Some(Named::Scalar(t))) => Ok(V::S(t)),
            Ok(Some(Named::Vector(c))) => Ok(V::Vec(c)),
            Ok(None) => {
                let known = self.names.known();
                self.fail(format!(
                    "unknown name '{name}'{}",
                    if known.is_empty() { String::new() } else { format!(" (the names are: {known}, pi, e)") }
                ))
            }
            Err(m) => self.fail(m),
        }
    }

    fn atom(&mut self) -> Result<V, String> {
        let Some(c) = self.peek() else {
            return self.fail("the formula ends where a value was expected");
        };
        if c == b'(' {
            self.i += 1;
            self.enter()?;
            let v = self.or()?;
            self.depth -= 1;
            if self.peek() != Some(b')') {
                return self.fail("missing ')'");
            }
            self.i += 1;
            return Ok(v);
        }
        if c == b'[' {
            self.i += 1;
            self.enter()?;
            let mut items = vec![];
            loop {
                let v = self.or()?;
                items.push(self.scalar(v, "a vector component")?);
                match self.peek() {
                    Some(b',') => self.i += 1,
                    Some(b']') => break,
                    _ => return self.fail("expected ',' or ']'"),
                }
            }
            self.i += 1;
            self.depth -= 1;
            if items.len() > 3 {
                return self.fail("a vector has at most three components (x, y, z)");
            }
            return Ok(V::Vec(items.into_iter().enumerate().collect()));
        }
        if c.is_ascii_digit() || c == b'.' {
            return self.number();
        }
        if c == b'`' {
            self.i += 1;
            let st = self.i;
            while self.i < self.s.len() && self.s[self.i] != b'`' {
                self.i += 1;
            }
            if self.i >= self.s.len() {
                return self.fail("a name in backticks is not closed");
            }
            let name = String::from_utf8_lossy(&self.s[st..self.i]).into_owned();
            self.i += 1;
            return self.named(&name);
        }
        if let Some(name) = self.peek_word() {
            self.i += name.len();
            if self.peek() == Some(b'(') {
                self.i += 1;
                return self.call(name);
            }
            return self.named(name);
        }
        self.fail(format!("unexpected '{}'", self.rest_char()))
    }

    fn args(&mut self) -> Result<Vec<V>, String> {
        let mut args = vec![];
        if self.peek() == Some(b')') {
            self.i += 1;
            return Ok(args);
        }
        self.enter()?;
        loop {
            args.push(self.or()?);
            match self.peek() {
                Some(b',') => self.i += 1,
                Some(b')') => break,
                _ => return self.fail("expected ',' or ')'"),
            }
        }
        self.i += 1;
        self.depth -= 1;
        Ok(args)
    }

    fn call(&mut self, name: &str) -> Result<V, String> {
        if name == "material" || name == "region" {
            let q = self.string()?;
            if self.peek() != Some(b')') {
                return self.fail(format!("{name}(\"name\") takes one quoted name"));
            }
            self.i += 1;
            let r = if name == "material" { self.names.material(&q) } else { self.names.region(&q) };
            return match r {
                Ok(t) => Ok(V::S(t)),
                Err(m) => self.fail(m),
            };
        }
        let args = self.args()?;
        let n = args.len();
        let arity = |want: usize, this: &Self| -> Result<(), String> {
            if n == want {
                Ok(())
            } else {
                this.fail(format!("{name} takes {want} argument{}, not {n}", if want == 1 { "" } else { "s" }))
            }
        };
        let fun = match name {
            "abs" => Some(Fun::Abs),
            "abs2" => Some(Fun::Abs2),
            "re" | "real" => Some(Fun::Re),
            "im" | "imag" => Some(Fun::Im),
            "conj" => Some(Fun::Conj),
            "sqrt" => Some(Fun::Sqrt),
            "exp" => Some(Fun::Exp),
            "sin" => Some(Fun::Sin),
            "cos" => Some(Fun::Cos),
            "arg" => Some(Fun::Arg),
            "log" | "ln" => Some(Fun::Log),
            "log10" => Some(Fun::Log10),
            "tan" => Some(Fun::Tan),
            "sinh" => Some(Fun::Sinh),
            "cosh" => Some(Fun::Cosh),
            "tanh" => Some(Fun::Tanh),
            "atan" | "arctan" => Some(Fun::Atan),
            "step" => Some(Fun::Step),
            _ => None,
        };
        if let Some(f) = fun {
            arity(1, self)?;
            let a = args.into_iter().next().expect("one argument");
            if matches!(a, V::Vec(_)) && matches!(f, Fun::Abs | Fun::Abs2) {
                return self.fail(format!("{name} of a vector: use norm(v) or norm2(v)"));
            }
            let t = self.scalar(a, name)?;
            return Ok(V::S(Tree::Call(f, b(t))));
        }
        let f2 = match name {
            "min" => Some(Fun2::Min),
            "max" => Some(Fun2::Max),
            "atan2" => Some(Fun2::Atan2),
            _ => None,
        };
        if let Some(f) = f2 {
            arity(2, self)?;
            let mut it = args.into_iter();
            let (x, y) = (it.next().expect("two"), it.next().expect("two"));
            let (x, y) = (self.scalar(x, name)?, self.scalar(y, name)?);
            return Ok(V::S(Tree::Call2(f, b(x), b(y))));
        }
        match name {
            "where" => {
                arity(3, self)?;
                let mut it = args.into_iter();
                let (c, x, y) = (it.next().expect("three"), it.next().expect("three"), it.next().expect("three"));
                let c = self.scalar(c, "the condition of where")?;
                match (x, y) {
                    (V::S(x), V::S(y)) => Ok(V::S(Tree::Where(b(c), b(x), b(y)))),
                    (V::Vec(x), V::Vec(y)) => {
                        if x.len() != y.len() || x.iter().zip(&y).any(|(p, q)| p.0 != q.0) {
                            return self.fail("the two branches of where are vectors with different components");
                        }
                        Ok(V::Vec(
                            x.into_iter()
                                .zip(y)
                                .map(|((a, p), (_, q))| (a, Tree::Where(b(c.clone()), b(p), b(q))))
                                .collect(),
                        ))
                    }
                    _ => self.fail("the two branches of where are a scalar and a vector"),
                }
            }
            "dot" | "inner" => {
                arity(2, self)?;
                let mut it = args.into_iter();
                let (x, y) = (it.next().expect("two"), it.next().expect("two"));
                let (V::Vec(x), V::Vec(y)) = (x, y) else {
                    return self.fail(format!("{name} takes two vectors"));
                };
                if x.len() != y.len() || x.iter().zip(&y).any(|(p, q)| p.0 != q.0) {
                    return self.fail(format!(
                        "{name} of vectors with components {} and {}: they need the same components",
                        axes_of(&x),
                        axes_of(&y)
                    ));
                }
                let terms = x.into_iter().zip(y).map(|((_, p), (_, q))| {
                    let p = if name == "inner" { Tree::Call(Fun::Conj, b(p)) } else { p };
                    Tree::Bin(Op::Mul, b(p), b(q))
                });
                Ok(V::S(sum_terms(terms)))
            }
            "cross" => {
                arity(2, self)?;
                let mut it = args.into_iter();
                let (x, y) = (it.next().expect("two"), it.next().expect("two"));
                let (V::Vec(x), V::Vec(y)) = (x, y) else {
                    return self.fail("cross takes two vectors");
                };
                if x.len() != 3 || y.len() != 3 || x.iter().chain(&y).enumerate().any(|(k, (a, _))| *a != k % 3) {
                    return self.fail("cross takes two vectors with the components x, y and z");
                }
                let m = |p: usize, q: usize| Tree::Bin(Op::Mul, b(x[p].1.clone()), b(y[q].1.clone()));
                let comp = |p: usize, q: usize| Tree::Bin(Op::Sub, b(m(p, q)), b(m(q, p)));
                Ok(V::Vec(vec![(0, comp(1, 2)), (1, comp(2, 0)), (2, comp(0, 1))]))
            }
            "norm2" | "norm" => {
                arity(1, self)?;
                let V::Vec(x) = args.into_iter().next().expect("one argument") else {
                    return self.fail(format!("{name} takes a vector (for a scalar use abs)"));
                };
                let s = sum_terms(x.into_iter().map(|(_, t)| Tree::Call(Fun::Abs2, b(t))));
                Ok(V::S(if name == "norm" { Tree::Call(Fun::RealSqrt, b(s)) } else { s }))
            }
            _ => self.fail(format!(
                "unknown function '{name}' (they are: abs abs2 re im conj sqrt exp sin cos arg log log10 tan sinh cosh tanh atan step min max atan2 where dot inner cross norm norm2 material region)"
            )),
        }
    }
}

/// `t0 + t1 + t2`, left to right.
fn sum_terms(mut it: impl Iterator<Item = Tree>) -> Tree {
    let first = it.next().expect("a term");
    it.fold(first, |acc, t| Tree::Bin(Op::Add, b(acc), b(t)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Leaves: 0 x, 1 y, 2 Ex, 3 Ez, 4 rho (a material value), 5 the region "core".
    struct Map;

    impl Names for Map {
        fn name(&mut self, name: &str) -> Result<Option<Named>, String> {
            Ok(match name {
                "x" => Some(Named::Scalar(Tree::Leaf(0))),
                "y" => Some(Named::Scalar(Tree::Leaf(1))),
                "Ex" => Some(Named::Scalar(Tree::Leaf(2))),
                "Ez" => Some(Named::Scalar(Tree::Leaf(3))),
                "E" => Some(Named::Vector(vec![(0, Tree::Leaf(2)), (2, Tree::Leaf(3))])),
                "rho" => Some(Named::Scalar(Tree::Leaf(4))),
                "|E|" => Some(Named::Scalar(Tree::Leaf(6))),
                _ => None,
            })
        }
        fn material(&mut self, name: &str) -> Result<Tree, String> {
            if name == "rho" {
                Ok(Tree::Leaf(4))
            } else {
                Err(format!("no material value '{name}'"))
            }
        }
        fn region(&mut self, name: &str) -> Result<Tree, String> {
            if name == "core" {
                Ok(Tree::Leaf(5))
            } else {
                Err(format!("no region '{name}'"))
            }
        }
        fn known(&self) -> String {
            "x, y, Ex, Ez, E, rho".into()
        }
    }

    fn eval(text: &str, leaf: &[C]) -> C {
        parse(text, &mut Map).unwrap_or_else(|e| panic!("{text}: {e}")).eval(leaf)
    }

    fn real_leaves() -> [C; 7] {
        [C::new(0.5, 0.0), C::new(-2.0, 0.0), C::new(3.0, 4.0), C::new(1.0, -1.0), C::new(2330.0, 0.0), C::new(1.0, 0.0), C::new(7.0, 0.0)]
    }

    #[test]
    fn arithmetic_keeps_the_order_and_precedence_of_the_real_grammar() {
        let l = real_leaves();
        assert_eq!(eval("2 + 3 * 4", &l).re, 14.0);
        assert_eq!(eval("(2 + 3) * 4", &l).re, 20.0);
        assert_eq!(eval("2 ^ 3 ^ 2", &l).re, 512.0);
        assert_eq!(eval("-2^2", &l).re, -4.0);
        assert_eq!(eval("2 ** 3", &l).re, 8.0);
        assert_eq!(eval("x - y - 1", &l).re, 0.5 + 2.0 - 1.0);
        assert_eq!(eval("8 / 4 / 2", &l).re, 1.0);
        assert_eq!(eval("2^-1", &l).re, 0.5);
        assert_eq!(eval("1.5e3 * x", &l).re, 750.0);
        assert_eq!(eval("pi", &l).re, std::f64::consts::PI);
        assert_eq!(eval("e", &l).re, std::f64::consts::E);
        assert_eq!(eval("  max( x , y )  ", &l).re, 0.5);
        assert_eq!(eval("2^0.5", &l).re, 2.0_f64.powf(0.5));
    }

    #[test]
    fn complex_values_and_the_functions_of_a_complex_value() {
        let l = real_leaves();
        // |Ex|^2 = 25, the same bits as abs2 and as the repeated product
        assert_eq!(eval("abs(Ex)^2", &l).re, 25.0);
        assert_eq!(eval("abs2(Ex)", &l).re, 25.0);
        assert_eq!(eval("Ex * conj(Ex)", &l), C::new(25.0, 0.0));
        assert_eq!(eval("re(Ex) + im(Ex)", &l).re, 7.0);
        assert_eq!(eval("2j", &l), C::new(0.0, 2.0));
        assert_eq!(eval("2j * 2j", &l), C::new(-4.0, 0.0));
        assert_eq!(eval("Ex^2", &l), C::new(3.0 * 3.0 - 4.0 * 4.0, 3.0 * 4.0 + 4.0 * 3.0));
        assert!(eval("Ex^-1 * Ex", &l).im.abs() < 1e-15);
        assert_eq!(eval("arg(1j)", &l).re, std::f64::consts::FRAC_PI_2);
        assert_eq!(eval("`|E|`", &l).re, 7.0);
    }

    #[test]
    fn where_and_comparisons() {
        let l = real_leaves();
        assert_eq!(eval("where(x > 0, 10, 20)", &l).re, 10.0);
        assert_eq!(eval("where(x < 0, 10, 20)", &l).re, 20.0);
        assert_eq!(eval("where(x > 0 and y < 0, 1, 2)", &l).re, 1.0);
        assert_eq!(eval("where(x < 0 or y < 0, 1, 2)", &l).re, 1.0);
        assert_eq!(eval("where(not (x < 0), 1, 2)", &l).re, 1.0);
        assert_eq!(eval("(x >= 0.5) + (x <= 0.5) + (x == 0.5) + (x != 0.5)", &l).re, 3.0);
        // the branch not taken is not evaluated
        assert_eq!(eval("where(x != 0.5, 1 / (x - 0.5), -1)", &l).re, -1.0);
        assert_eq!(eval("where(region(\"core\"), rho, 0)", &l).re, 2330.0);
        assert_eq!(eval("material(\"rho\") / 2", &l).re, 1165.0);
    }

    #[test]
    fn vectors_dot_cross_norm_and_components() {
        let l = real_leaves();
        // E = (3+4i, 0 in y, 1-i)
        assert_eq!(eval("E.x", &l), C::new(3.0, 4.0));
        assert_eq!(eval("E.z", &l), C::new(1.0, -1.0));
        assert_eq!(eval("norm2(E)", &l).re, 25.0 + 2.0);
        assert_eq!(eval("norm(E)", &l).re, 27.0_f64.sqrt());
        // inner conjugates the first vector: sum |E_k|^2
        assert_eq!(eval("inner(E, E)", &l), C::new(27.0, 0.0));
        // dot does not: E.E = (3+4i)^2 + (1-i)^2
        assert_eq!(eval("dot(E, E)", &l), C::new(-7.0 + 0.0, 24.0 - 2.0));
        assert_eq!(eval("(2 * E).x", &l), C::new(6.0, 8.0));
        assert_eq!(eval("(E / 2).z", &l), C::new(0.5, -0.5));
        assert_eq!(eval("(E + E).x", &l), C::new(6.0, 8.0));
        assert_eq!(eval("(-E).x", &l), C::new(-3.0, -4.0));
        assert_eq!(eval("dot([1, 2, 3], [4, 5, 6])", &l).re, 32.0);
        assert_eq!(eval("cross([1, 0, 0], [0, 1, 0]).z", &l).re, 1.0);
        assert_eq!(eval("cross([1, 2, 3], [4, 5, 6]).x", &l).re, 2.0 * 6.0 - 3.0 * 5.0);
        assert_eq!(eval("cross([1, 2, 3], [4, 5, 6]).y", &l).re, 3.0 * 4.0 - 1.0 * 6.0);
        assert_eq!(eval("where(x > 0, E, 0 * E).x", &l), C::new(3.0, 4.0));
    }

    #[test]
    fn refused_inputs_say_what_is_wrong() {
        for (bad, why) in [
            ("", "empty"),
            ("   ", "empty"),
            ("1 +", "ends"),
            ("(1 + 2", "missing ')'"),
            ("1 + 2)", "unexpected ')'"),
            ("foo", "unknown name 'foo'"),
            ("foo(1)", "unknown function 'foo'"),
            ("sin()", "takes 1 argument"),
            ("min(1)", "takes 2 arguments"),
            ("where(1, 2)", "takes 3 arguments"),
            ("1 2", "unexpected '2'"),
            ("x y", "unexpected 'y'"),
            ("x = 1", "single '='"),
            ("0 < x < 1", "does not chain"),
            ("__import__('os')", "unexpected"),
            ("x.real", "component"),
            ("E", "is a vector"),
            ("E.y", "has the components x, z"),
            ("E * E", "dot(a, b)"),
            ("E + 1", "do not add"),
            ("1 / E", "does not divide"),
            ("E ^ 2", "needs a scalar"),
            ("abs(E)", "norm(v)"),
            ("dot(1, 2)", "two vectors"),
            ("cross(E, E)", "x, y and z"),
            ("norm(3)", "takes a vector"),
            ("material(\"nope\")", "no material value"),
            ("region(\"nope\")", "no region"),
            ("material(rho)", "quoted name"),
            ("[1, 2, 3, 4]", "at most three"),
            ("1 $ 2", "unexpected '$'"),
            ("`Ex", "not closed"),
            ("1..2", "bad number"),
        ] {
            let e = parse(bad, &mut Map).expect_err(bad);
            assert!(e.contains(why), "{bad:?}: {e} (wanted {why:?})");
        }
        // an unknown name lists the names the caller has
        assert!(parse("foo", &mut Map).unwrap_err().contains("x, y, Ex, Ez, E, rho"));
        // deep nesting is refused, not a stack overflow
        let deep = format!("{}1{}", "(".repeat(200), ")".repeat(200));
        assert!(parse(&deep, &mut Map).unwrap_err().contains("nested too deeply"));
        let neg = format!("{}1", "-".repeat(200));
        assert!(parse(&neg, &mut Map).unwrap_err().contains("nested too deeply"));
    }

    #[test]
    fn a_power_with_a_whole_exponent_is_the_repeated_product() {
        let t = parse("Ex^3", &mut Map).unwrap();
        assert_eq!(t, Tree::Pow(Box::new(Tree::Leaf(2)), 3));
        assert_eq!(parse("Ex^1", &mut Map).unwrap(), Tree::Leaf(2));
        assert_eq!(parse("Ex^0", &mut Map).unwrap(), Tree::num(1.0));
        // abs(Ex)^2 is the tree the handle `pow` builds
        assert_eq!(
            parse("abs(Ex)^2", &mut Map).unwrap(),
            Tree::Pow(Box::new(Tree::Call(Fun::Abs, Box::new(Tree::Leaf(2)))), 2)
        );
        assert!(matches!(parse("Ex^17", &mut Map).unwrap(), Tree::Call2(Fun2::Pow, _, _)));
    }

    #[test]
    fn real_only_operations_are_refused_on_complex_values_before_any_point() {
        let t = parse("where(Ex > 0, 1, 2)", &mut Map).unwrap();
        let e = t.check(&[false, false, true, false, false, false, false]).unwrap_err();
        assert!(e.contains("complex") && e.contains("re"), "{e}");
        let ok = parse("where(re(Ex) > 0, 1, 2)", &mut Map).unwrap();
        assert!(ok.check(&[false, false, true, false, false, false, false]).is_ok());
        assert!(parse("log(Ex)", &mut Map).unwrap().check(&[false, false, true, false, false, false, false]).is_err());
    }
}
