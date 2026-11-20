// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A small, safe expression evaluator for user input such as `293.15 + 50*exp(-((x-0.5)^2)/0.01)` (initial temperature
//! fields, for example). Numbers, `+ - * / ^`, parentheses, the functions `sin cos tan exp log ln log10 sqrt abs sinh cosh tanh atan min max step`
//! and the names `x y z pi e`; nothing else, so no input can run code.

#[derive(Clone, Debug)]
enum Node {
    Num(f64),
    Var(u8), // 0 x, 1 y, 2 z
    Neg(Box<Node>),
    Bin(char, Box<Node>, Box<Node>),
    Call(Func, Vec<Node>),
}

#[derive(Clone, Copy, Debug)]
enum Func {
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
    Log10,
    Sqrt,
    Abs,
    Tanh,
    Sinh,
    Cosh,
    Atan,
    Min,
    Max,
    Step,
}

#[derive(Clone, Debug)]
pub struct Expr {
    root: Node,
    pub source: String,
}

struct P<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> P<'a> {
    fn ws(&mut self) {
        while self.i < self.s.len() && (self.s[self.i] as char).is_whitespace() {
            self.i += 1;
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.s.get(self.i).copied()
    }
    fn expr(&mut self) -> Result<Node, String> {
        let mut l = self.term()?;
        while let Some(c) = self.peek() {
            if c == b'+' || c == b'-' {
                self.i += 1;
                let r = self.term()?;
                l = Node::Bin(c as char, Box::new(l), Box::new(r));
            } else {
                break;
            }
        }
        Ok(l)
    }
    fn term(&mut self) -> Result<Node, String> {
        let mut l = self.unary()?;
        while let Some(c) = self.peek() {
            if c == b'*' && self.s.get(self.i + 1) != Some(&b'*') || c == b'/' {
                self.i += 1;
                let r = self.unary()?;
                l = Node::Bin(c as char, Box::new(l), Box::new(r));
            } else {
                break;
            }
        }
        Ok(l)
    }
    fn unary(&mut self) -> Result<Node, String> {
        match self.peek() {
            Some(b'-') => {
                self.i += 1;
                Ok(Node::Neg(Box::new(self.unary()?)))
            }
            Some(b'+') => {
                self.i += 1;
                self.unary()
            }
            _ => self.power(),
        }
    }
    fn power(&mut self) -> Result<Node, String> {
        let base = self.atom()?;
        match self.peek() {
            Some(b'^') => {
                self.i += 1;
                let e = self.unary()?;
                Ok(Node::Bin('^', Box::new(base), Box::new(e)))
            }
            Some(b'*') if self.s.get(self.i + 1) == Some(&b'*') => {
                self.i += 2;
                let e = self.unary()?;
                Ok(Node::Bin('^', Box::new(base), Box::new(e)))
            }
            _ => Ok(base),
        }
    }
    fn atom(&mut self) -> Result<Node, String> {
        let c = self.peek().ok_or("unexpected end of expression")?;
        if c == b'(' {
            self.i += 1;
            let e = self.expr()?;
            if self.peek() != Some(b')') {
                return Err("missing ')'".into());
            }
            self.i += 1;
            return Ok(e);
        }
        if c.is_ascii_digit() || c == b'.' {
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
            let t = std::str::from_utf8(&self.s[st..self.i]).unwrap();
            return t.parse::<f64>().map(Node::Num).map_err(|_| format!("bad number '{t}'"));
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let st = self.i;
            while self.i < self.s.len() && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'_') {
                self.i += 1;
            }
            let name = std::str::from_utf8(&self.s[st..self.i]).unwrap().to_string();
            if self.peek() == Some(b'(') {
                self.i += 1;
                let mut args = vec![];
                if self.peek() != Some(b')') {
                    loop {
                        args.push(self.expr()?);
                        match self.peek() {
                            Some(b',') => self.i += 1,
                            Some(b')') => break,
                            _ => return Err("expected ',' or ')'".into()),
                        }
                    }
                }
                self.i += 1;
                let (f, arity) = match name.as_str() {
                    "sin" => (Func::Sin, 1),
                    "cos" => (Func::Cos, 1),
                    "tan" => (Func::Tan, 1),
                    "exp" => (Func::Exp, 1),
                    "log" | "ln" => (Func::Log, 1),
                    "log10" => (Func::Log10, 1),
                    "sqrt" => (Func::Sqrt, 1),
                    "abs" => (Func::Abs, 1),
                    "tanh" => (Func::Tanh, 1),
                    "sinh" => (Func::Sinh, 1),
                    "cosh" => (Func::Cosh, 1),
                    "atan" | "arctan" => (Func::Atan, 1),
                    "min" => (Func::Min, 2),
                    "max" => (Func::Max, 2),
                    "step" => (Func::Step, 1),
                    _ => return Err(format!("unknown function '{name}'")),
                };
                if args.len() != arity {
                    return Err(format!("{name} takes {arity} argument(s)"));
                }
                return Ok(Node::Call(f, args));
            }
            return match name.as_str() {
                "x" => Ok(Node::Var(0)),
                "y" => Ok(Node::Var(1)),
                "z" => Ok(Node::Var(2)),
                "pi" => Ok(Node::Num(std::f64::consts::PI)),
                "e" => Ok(Node::Num(std::f64::consts::E)),
                _ => Err(format!("unknown name '{name}' (use x, y, z, pi, e)")),
            };
        }
        Err(format!("unexpected '{}'", c as char))
    }
}

impl Expr {
    pub fn parse(text: &str) -> Result<Expr, String> {
        let mut p = P { s: text.as_bytes(), i: 0 };
        if text.trim().is_empty() {
            return Err("empty expression".into());
        }
        let root = p.expr()?;
        p.ws();
        if p.i != p.s.len() {
            return Err(format!("unexpected '{}' at position {}", p.s[p.i] as char, p.i + 1));
        }
        Ok(Expr { root, source: text.to_string() })
    }

    pub fn eval(&self, x: f64, y: f64, z: f64) -> f64 {
        fn ev(n: &Node, v: [f64; 3]) -> f64 {
            match n {
                Node::Num(c) => *c,
                Node::Var(i) => v[*i as usize],
                Node::Neg(a) => -ev(a, v),
                Node::Bin(op, a, b) => {
                    let (l, r) = (ev(a, v), ev(b, v));
                    match op {
                        '+' => l + r,
                        '-' => l - r,
                        '*' => l * r,
                        '/' => l / r,
                        _ => l.powf(r),
                    }
                }
                Node::Call(f, a) => {
                    let x = ev(&a[0], v);
                    match f {
                        Func::Sin => x.sin(),
                        Func::Cos => x.cos(),
                        Func::Tan => x.tan(),
                        Func::Exp => x.exp(),
                        Func::Log => x.ln(),
                        Func::Log10 => x.log10(),
                        Func::Sqrt => x.sqrt(),
                        Func::Abs => x.abs(),
                        Func::Tanh => x.tanh(),
                        Func::Sinh => x.sinh(),
                        Func::Cosh => x.cosh(),
                        Func::Atan => x.atan(),
                        Func::Step => if x >= 0.0 { 1.0 } else { 0.0 },
                        Func::Min => x.min(ev(&a[1], v)),
                        Func::Max => x.max(ev(&a[1], v)),
                    }
                }
            }
        }
        ev(&self.root, [x, y, z])
    }
}

/// Evaluates `text` at dof coordinates (flat, `dim` per point) in the wafer axes: `x`, `y`, `z`. On a 2D cross-section the mesh is the x-z
/// plane, so the second mesh axis is `z` and `y` is 0.
pub fn eval_at(text: &str, coords: &[f64], dim: usize) -> Result<Vec<f64>, String> {
    let e = Expr::parse(text)?;
    Ok(coords
        .chunks(dim)
        .map(|p| {
            let w = crate::section::wafer_point(p, dim);
            e.eval(w[0], w[1], w[2])
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates() {
        let e = Expr::parse("sin(pi*x)+2*x^2 - -y").unwrap();
        let v = e.eval(0.5, 3.0, 0.0);
        assert!((v - (1.0 + 0.5 + 3.0)).abs() < 1e-14);
        assert_eq!(Expr::parse("2^3^2").unwrap().eval(0.0, 0.0, 0.0), 512.0);
        assert_eq!(Expr::parse("-2^2").unwrap().eval(0.0, 0.0, 0.0), -4.0);
        assert!((Expr::parse("1.5e3*x").unwrap().eval(2.0, 0.0, 0.0) - 3000.0).abs() < 1e-9);
        assert_eq!(Expr::parse("2**3").unwrap().eval(0.0, 0.0, 0.0), 8.0);
        assert_eq!(Expr::parse("max(x, 3)").unwrap().eval(1.0, 0.0, 0.0), 3.0);
        for bad in ["__import__('os')", "x.real", "foo(1)", "1 +", "(1", "x y", ""] {
            assert!(Expr::parse(bad).is_err(), "{bad}");
        }
    }
}
