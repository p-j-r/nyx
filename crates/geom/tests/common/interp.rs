//! A tiny interpreter for the code `FeatureList::to_rust` prints (builder shorthand only), so the round trip
//! `build(print(list)) == list` can be checked as pure data without compiling the printed code.
#![allow(dead_code)]

use geom::expr::Expr;
use geom::model::*;
use geom::select::Selector;

#[derive(Debug, Clone)]
enum A {
    Str(String),
    Num(f64),
    Id(String),
    Path(String),
    List(Vec<A>),
    Call(String, Vec<A>),
}

struct P<'a> {
    s: &'a [u8],
    i: usize,
}

impl P<'_> {
    fn ws(&mut self) {
        while self.s.get(self.i).is_some_and(|c| c.is_ascii_whitespace()) {
            self.i += 1;
        }
    }
    fn arg(&mut self) -> A {
        self.ws();
        let c = self.s[self.i];
        match c {
            b'"' => {
                self.i += 1;
                let mut out = String::new();
                while self.s[self.i] != b'"' {
                    if self.s[self.i] == b'\\' {
                        self.i += 1;
                        out.push(match self.s[self.i] {
                            b'n' => '\n',
                            b't' => '\t',
                            o => o as char,
                        });
                    } else {
                        let ch = std::str::from_utf8(&self.s[self.i..]).unwrap().chars().next().unwrap();
                        out.push(ch);
                        self.i += ch.len_utf8() - 1;
                    }
                    self.i += 1;
                }
                self.i += 1;
                A::Str(out)
            }
            b'(' | b'[' => {
                self.i += 1;
                A::List(self.items(if c == b'(' { b')' } else { b']' }))
            }
            b'&' => {
                self.i += 1;
                self.arg()
            }
            c if c.is_ascii_digit() || c == b'-' => {
                let st = self.i;
                self.i += 1;
                while self.s.get(self.i).is_some_and(|c| c.is_ascii_digit() || b".eE+-".contains(c)) {
                    self.i += 1;
                }
                A::Num(std::str::from_utf8(&self.s[st..self.i]).unwrap().parse().unwrap())
            }
            _ => {
                let st = self.i;
                while self.s.get(self.i).is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b':') {
                    self.i += 1;
                }
                let name = std::str::from_utf8(&self.s[st..self.i]).unwrap().to_string();
                assert!(!name.is_empty(), "unexpected character {:?} in printed code", self.s[self.i] as char);
                if self.s.get(self.i) == Some(&b'(') {
                    self.i += 1;
                    let c = A::Call(name, self.items(b')'));
                    if self.s.get(self.i) == Some(&b'?') {
                        self.i += 1;
                    }
                    c
                } else if self.s.get(self.i) == Some(&b'!') {
                    self.i += 2; // vec![
                    A::List(self.items(b']'))
                } else if name.contains("::") {
                    A::Path(name)
                } else {
                    A::Id(name)
                }
            }
        }
    }
    fn items(&mut self, close: u8) -> Vec<A> {
        let mut v = vec![];
        loop {
            self.ws();
            if self.s[self.i] == close {
                self.i += 1;
                return v;
            }
            v.push(self.arg());
            self.ws();
            if self.s[self.i] == b',' {
                self.i += 1;
            }
        }
    }
}

fn expr(a: &A) -> Expr {
    match a {
        A::Num(v) => Expr::Num(*v),
        A::Str(s) => Expr::Text(s.clone()),
        A::Call(n, x) if n == "Expr::from" => expr(&x[0]),
        o => panic!("not an expression: {o:?}"),
    }
}
fn e2(a: &A) -> (Expr, Expr) {
    let A::List(v) = a else { panic!("tuple expected: {a:?}") };
    (expr(&v[0]), expr(&v[1]))
}
fn e3(a: &A) -> (Expr, Expr, Expr) {
    let A::List(v) = a else { panic!("tuple expected: {a:?}") };
    (expr(&v[0]), expr(&v[1]), expr(&v[2]))
}
fn st(a: &A) -> String {
    match a {
        A::Str(s) => s.clone(),
        o => panic!("string expected: {o:?}"),
    }
}
fn num(a: &A) -> f64 {
    match a {
        A::Num(v) => *v,
        o => panic!("number expected: {o:?}"),
    }
}
fn fid(a: &A, ids: &std::collections::HashMap<String, FeatureId>) -> FeatureId {
    match a {
        A::Id(n) => ids[n],
        o => panic!("feature expected: {o:?}"),
    }
}
fn fids(a: &A, ids: &std::collections::HashMap<String, FeatureId>) -> Vec<FeatureId> {
    let A::List(v) = a else { panic!() };
    v.iter().map(|x| fid(x, ids)).collect()
}
fn sel(a: &A) -> Selector {
    let A::Call(_, v) = a else { panic!("selector expected") };
    Selector::parse(&st(&v[0])).unwrap()
}
fn opt_expr(a: &A) -> Option<Expr> {
    match a {
        A::Id(n) if n == "None" => None,
        A::Call(n, v) if n == "Some" => Some(expr(&v[0])),
        o => panic!("option expected: {o:?}"),
    }
}

pub fn build(code: &str) -> Model {
    let mut m: Option<Model> = None;
    let mut ids = std::collections::HashMap::new();
    for line in code.lines() {
        let line = line.trim().trim_end_matches(';');
        if let Some(r) = line.strip_prefix("let mut m = Model::new_") {
            m = Some(if r.starts_with('2') { Model::new_2d() } else { Model::new_3d() });
            continue;
        }
        let m = m.as_mut().unwrap();
        let (var, call) = match line.strip_prefix("let ") {
            Some(r) => {
                let (v, c) = r.split_once(" = ").unwrap();
                (Some(v.to_string()), c)
            }
            None => (None, line),
        };
        let call = call.strip_prefix("m.").unwrap();
        let mut p = P { s: call.as_bytes(), i: 0 };
        let A::Call(name, a) = p.arg() else { panic!("call expected: {call}") };
        let id = match name.as_str() {
            "param" => {
                m.param(&st(&a[0]), expr(&a[1]));
                None
            }
            "rect" => Some(m.rect(&st(&a[0]), e2(&a[1]), e2(&a[2]))),
            "rounded_rect" => Some(m.rounded_rect(&st(&a[0]), e2(&a[1]), e2(&a[2]), expr(&a[3]))),
            "circle" => Some(m.circle(&st(&a[0]), e2(&a[1]), expr(&a[2]))),
            "revolve_along" => {
                let A::Path(k) = &a[1] else { panic!() };
                Some(m.revolve_along(&st(&a[0]), axis_of(k), fid(&a[2], &ids), e3(&a[3])))
            }
            "revolve" => Some(m.revolve(&st(&a[0]), fid(&a[1], &ids), e3(&a[2]))),
            "cone_along" => {
                let A::Path(k) = &a[1] else { panic!() };
                Some(m.cone_along(&st(&a[0]), axis_of(k), e3(&a[2]), expr(&a[3]), expr(&a[4]), expr(&a[5])))
            }
            "cone" => Some(m.cone(&st(&a[0]), e3(&a[1]), expr(&a[2]), expr(&a[3]), expr(&a[4]))),
            "sphere" => Some(m.sphere(&st(&a[0]), e3(&a[1]), expr(&a[2]))),
            "torus" => Some(m.torus(&st(&a[0]), e3(&a[1]), expr(&a[2]), expr(&a[3]))),
            "box_" => Some(m.box_(&st(&a[0]), e3(&a[1]), e3(&a[2]))),
            "cylinder" => Some(m.cylinder(&st(&a[0]), e3(&a[1]), expr(&a[2]), expr(&a[3]))),
            "extrude_tapered" => {
                Some(m.extrude_tapered(&st(&a[0]), fid(&a[1], &ids), expr(&a[2]), expr(&a[3]), expr(&a[4])))
            }
            "cylinder_along" => {
                let A::Path(k) = &a[1] else { panic!() };
                Some(m.cylinder_along(&st(&a[0]), axis_of(k), e3(&a[2]), expr(&a[3]), expr(&a[4])))
            }
            "extrude_along" => {
                let A::Path(k) = &a[1] else { panic!() };
                Some(m.extrude_along(&st(&a[0]), axis_of(k), fid(&a[2], &ids), expr(&a[3]), expr(&a[4])))
            }
            "extrude" => Some(m.extrude(&st(&a[0]), fid(&a[1], &ids), expr(&a[2]), expr(&a[3]))),
            "ellipse" => Some(m.ellipse(&st(&a[0]), e2(&a[1]), e2(&a[2]), expr(&a[3]))),
            "polygon" | "polygon_named" => {
                let A::List(pts) = &a[1] else { panic!() };
                let pts: Vec<(Expr, Expr)> = pts.iter().map(e2).collect();
                if name == "polygon" {
                    Some(m.polygon(&st(&a[0]), &pts))
                } else {
                    let A::List(nm) = &a[2] else { panic!() };
                    let nm: Vec<String> = nm.iter().map(st).collect();
                    let nm: Vec<&str> = nm.iter().map(String::as_str).collect();
                    Some(m.polygon_named(&st(&a[0]), &pts, &nm))
                }
            }
            "translate" => Some(m.translate(&st(&a[0]), fid(&a[1], &ids), e3(&a[2]))),
            "linear_pattern" => Some(m.linear_pattern(&st(&a[0]), fid(&a[1], &ids), e3(&a[2]), expr(&a[3]))),
            "union" => Some(m.union(&st(&a[0]), &fids(&a[1], &ids))),
            "intersect" => Some(m.intersect(&st(&a[0]), &fids(&a[1], &ids))),
            "difference" => Some(m.difference(&st(&a[0]), fid(&a[1], &ids), &fids(&a[2], &ids))),
            "overlay" => {
                let A::List(r) = &a[0] else { panic!() };
                let r: Vec<(String, FeatureId)> = r
                    .iter()
                    .map(|x| {
                        let A::List(t) = x else { panic!() };
                        (st(&t[0]), fid(&t[1], &ids))
                    })
                    .collect();
                let r: Vec<(&str, FeatureId)> = r.iter().map(|(n, i)| (n.as_str(), *i)).collect();
                Some(m.overlay(&r))
            }
            "overlay_full" => {
                let A::List(r) = &a[0] else { panic!() };
                let r: Vec<(String, FeatureId, RegionKind, Option<i32>)> = r
                    .iter()
                    .map(|x| {
                        let A::List(t) = x else { panic!() };
                        let A::Path(k) = &t[2] else { panic!() };
                        let kind = match k.as_str() {
                            "RegionKind::Material" => RegionKind::Material,
                            "RegionKind::Hole" => RegionKind::Hole,
                            "RegionKind::Pml" => RegionKind::Pml,
                            _ => RegionKind::Void,
                        };
                        let pin = match &t[3] {
                            A::Call(_, v) => Some(num(&v[0]) as i32),
                            _ => None,
                        };
                        (st(&t[0]), fid(&t[1], &ids), kind, pin)
                    })
                    .collect();
                let r: Vec<(&str, FeatureId, RegionKind, Option<i32>)> =
                    r.iter().map(|(n, i, k, p)| (n.as_str(), *i, *k, *p)).collect();
                let A::List(pn) = &a[1] else { panic!() };
                let pn: Vec<(String, i32)> = pn
                    .iter()
                    .map(|x| {
                        let A::List(t) = x else { panic!() };
                        (st(&t[0]), num(&t[1]) as i32)
                    })
                    .collect();
                let pn: Vec<(&str, i32)> = pn.iter().map(|(s, i)| (s.as_str(), *i)).collect();
                Some(m.overlay_full(&r, &pn))
            }
            "pml" => {
                let A::List(sd) = &a[2] else { panic!() };
                let sides: Vec<Side> = sd
                    .iter()
                    .map(|x| {
                        let A::Path(p) = x else { panic!() };
                        match p.as_str() {
                            "Side::XMin" => Side::XMin,
                            "Side::XMax" => Side::XMax,
                            "Side::YMin" => Side::YMin,
                            "Side::YMax" => Side::YMax,
                            "Side::ZMin" => Side::ZMin,
                            _ => Side::ZMax,
                        }
                    })
                    .collect();
                Some(m.pml(&st(&a[0]), fid(&a[1], &ids), &sides, expr(&a[3]), expr(&a[4]), expr(&a[5])))
            }
            "cut" => {
                let A::List(nv) = &a[1] else { panic!() };
                let A::Path(k) = &a[3] else { panic!() };
                Some(m.cut(
                    &st(&a[0]),
                    (num(&nv[0]), num(&nv[1]), num(&nv[2])),
                    expr(&a[2]),
                    if k == "Keep::Above" { Keep::Above } else { Keep::Below },
                ))
            }
            "periodic" => Some(m.periodic(&st(&a[0]), sel(&a[1]), sel(&a[2]), e3(&a[3]))),
            "size" => Some(m.size(&st(&a[0]), sel(&a[1]), opt_expr(&a[2]), opt_expr(&a[3]), opt_expr(&a[4]))),
            "tag" => Some(m.tag(&st(&a[0]), sel(&a[1]), &st(&a[2]))),
            "set_naming" => {
                let A::Path(k) = &a[0] else { panic!() };
                m.set_naming(match k.as_str() {
                    "NameScheme::Region" => NameScheme::Region,
                    "NameScheme::Numbered" => NameScheme::Numbered,
                    _ => NameScheme::Feature,
                });
                None
            }
            "suppress" => {
                m.suppress(fid(&a[0], &ids), true);
                None
            }
            "tag_feature" => {
                m.tag_feature(fid(&a[0], &ids), &st(&a[1]));
                None
            }
            o => panic!("unsupported call {o}"),
        };
        if let (Some(v), Some(i)) = (var, id) {
            ids.insert(v, i);
        }
    }
    m.unwrap()
}

fn axis_of(p: &str) -> Axis {
    match p {
        "Axis::X" => Axis::X,
        "Axis::Y" => Axis::Y,
        _ => Axis::Z,
    }
}
