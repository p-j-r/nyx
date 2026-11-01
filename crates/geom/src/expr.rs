//! Parameter expressions: a number or a small safe formula over named parameters. The grammar mirrors nyx `femcore::expr`
//! (numbers, `+ - * / ^ **`, parentheses, the functions `sin cos tan exp log ln log10 sqrt abs sinh cosh tanh atan min max step`,
//! and the constants `pi e`) with parameter names in place of `x y z`. Nothing in an expression can run code.

use crate::error::{GeomError, GeomResult};
use serde::{Deserialize, Serialize};

/// A number or an expression string; serialises as a bare JSON number or string (`untagged`), which `xmlio` writes as `t="n"` / `t="s"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Expr {
    Num(f64),
    Text(String),
}

impl From<f64> for Expr {
    fn from(v: f64) -> Expr {
        Expr::Num(v)
    }
}
impl From<i32> for Expr {
    fn from(v: i32) -> Expr {
        Expr::Num(v as f64)
    }
}
impl From<&str> for Expr {
    fn from(s: &str) -> Expr {
        Expr::Text(s.to_string())
    }
}
impl From<String> for Expr {
    fn from(s: String) -> Expr {
        Expr::Text(s)
    }
}
impl From<&Expr> for Expr {
    fn from(e: &Expr) -> Expr {
        e.clone()
    }
}

impl Expr {
    /// Source text for printing (numbers use Rust's shortest round-trip form).
    pub fn source(&self) -> String {
        match self {
            Expr::Num(v) => format!("{v:?}"),
            Expr::Text(s) => s.clone(),
        }
    }
}

#[derive(Clone, Debug)]
enum Node {
    Num(f64),
    Var(String),
    Neg(Box<Node>),
    Bin(u8, Box<Node>, Box<Node>),
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

struct P<'a> {
    s: &'a [u8],
    i: usize,
}

impl P<'_> {
    fn ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
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
                l = Node::Bin(c, Box::new(l), Box::new(r));
            } else {
                break;
            }
        }
        Ok(l)
    }
    fn term(&mut self) -> Result<Node, String> {
        let mut l = self.unary()?;
        while let Some(c) = self.peek() {
            if (c == b'*' && self.s.get(self.i + 1) != Some(&b'*')) || c == b'/' {
                self.i += 1;
                let r = self.unary()?;
                l = Node::Bin(c, Box::new(l), Box::new(r));
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
                Ok(Node::Bin(b'^', Box::new(base), Box::new(e)))
            }
            Some(b'*') if self.s.get(self.i + 1) == Some(&b'*') => {
                self.i += 2;
                let e = self.unary()?;
                Ok(Node::Bin(b'^', Box::new(base), Box::new(e)))
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
            let t = std::str::from_utf8(&self.s[st..self.i]).unwrap_or("");
            return t.parse::<f64>().map(Node::Num).map_err(|_| format!("bad number '{t}'"));
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let st = self.i;
            while self.i < self.s.len() && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'_') {
                self.i += 1;
            }
            let name = std::str::from_utf8(&self.s[st..self.i]).unwrap_or("").to_string();
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
            return Ok(match name.as_str() {
                "pi" => Node::Num(std::f64::consts::PI),
                "e" => Node::Num(std::f64::consts::E),
                _ => Node::Var(name),
            });
        }
        Err(format!("unexpected '{}'", c as char))
    }
}

fn parse(text: &str) -> Result<Node, String> {
    if text.trim().is_empty() {
        return Err("empty expression".into());
    }
    let mut p = P { s: text.as_bytes(), i: 0 };
    let root = p.expr()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(format!("unexpected '{}' at position {}", p.s[p.i] as char, p.i + 1));
    }
    Ok(root)
}

fn ev(n: &Node, env: &dyn Fn(&str) -> Option<f64>) -> Result<f64, String> {
    Ok(match n {
        Node::Num(c) => *c,
        Node::Var(v) => env(v).ok_or_else(|| format!("unknown name '{v}'"))?,
        Node::Neg(a) => -ev(a, env)?,
        Node::Bin(op, a, b) => {
            let (l, r) = (ev(a, env)?, ev(b, env)?);
            match op {
                b'+' => l + r,
                b'-' => l - r,
                b'*' => l * r,
                b'/' => l / r,
                _ => l.powf(r),
            }
        }
        Node::Call(f, a) => {
            let x = ev(&a[0], env)?;
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
                Func::Step => (x >= 0.0) as u8 as f64,
                Func::Min => x.min(ev(&a[1], env)?),
                Func::Max => x.max(ev(&a[1], env)?),
            }
        }
    })
}

/// Evaluates `text` with names resolved by `env`; the result must be finite.
pub fn eval_text(text: &str, env: &dyn Fn(&str) -> Option<f64>) -> Result<f64, String> {
    let v = ev(&parse(text).map_err(|e| format!("in expression '{text}': {e}"))?, env)
        .map_err(|e| format!("in expression '{text}': {e}"))?;
    if v.is_finite() {
        Ok(v)
    } else {
        Err(format!("expression '{text}' is not a finite number ({v})"))
    }
}

/// Named parameters evaluated in order; later parameters may use earlier ones.
#[derive(Clone, Debug, Default)]
pub struct Env {
    pub vars: Vec<(String, f64)>,
}

impl Env {
    pub fn get(&self, name: &str) -> Option<f64> {
        self.vars.iter().rev().find(|(n, _)| n == name).map(|(_, v)| *v)
    }

    pub fn eval(&self, e: &Expr, what: &str) -> GeomResult<f64> {
        match e {
            Expr::Num(v) if v.is_finite() => Ok(*v),
            Expr::Num(v) => Err(GeomError(format!("{what} is not a finite number ({v})"))),
            Expr::Text(s) => eval_text(s, &|n| self.get(n)).map_err(|m| GeomError(format!("{what}: {m}"))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates() {
        let env = Env { vars: vec![("w".into(), 0.5), ("t".into(), 3.0)] };
        let f = |s: &str| env.eval(&Expr::Text(s.into()), "x");
        assert!((f("-w/2-1e-6").unwrap() + 0.250001).abs() < 1e-15);
        assert_eq!(f("2^3^2").unwrap(), 512.0);
        assert_eq!(f("-2^2").unwrap(), -4.0);
        assert_eq!(f("max(w, t)").unwrap(), 3.0);
        assert_eq!(f("2**3").unwrap(), 8.0);
        for bad in ["__import__('os')", "w.real", "foo(1)", "1 +", "(1", "w t", "", "zz", "1/0", "sqrt(-1)"] {
            assert!(f(bad).is_err(), "{bad}");
        }
        assert!(env.eval(&Expr::Num(f64::NAN), "radius of core").unwrap_err().0.contains("radius of core"));
    }
}
