//! Selectors: small serialisable rules that resolve against an evaluated complex to typed entity sets (regions, facets,
//! vertices). The text form (`name(core.top) | interface(Si, BOX)`) parses to the serde form and prints back deterministically.
//! `suggest` goes the other way: for a picked set it returns the simplest selectors that resolve to exactly that set.
//!
//! Facets are the codimension-1 entities: edges in a 2D model, faces in a 3D model (reported with kind `Edges` / `Faces`).
//! `A & extreme(+z)` evaluates `extreme` inside the set `A` (so "the top of the substrate" is `boundary(region(S)) & extreme(+z)`).
//! Coordinates are in the wafer frame: x lateral, y along the propagation direction, z the surface normal. A 2D model is the x-z
//! cross-section (its entities have y = 0).
//! A selector that matches nothing is not an error: the result carries a warning naming the empty atom.

use crate::error::{GeomError, GeomResult};
use crate::eval::Evaluated;
use crate::expr::Expr;
use crate::math::*;
use crate::model::Axis;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoxMode {
    Within,
    Touching,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurfKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Line,
    Circle,
    Ellipse,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum NormalSpec {
    /// `+x`, `-z`, ...
    Axis {
        axis: Axis,
        positive: bool,
    },
    Vector([Expr; 3]),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Selector {
    Name(String),
    Tag(String),
    Region(String),
    Boundary(Box<Selector>),
    Outer,
    /// `None`: any neighbour (including the void)
    Interface(String, Option<String>),
    From(String, String),
    Kind(SurfKind),
    Radius(Expr, Option<Expr>),
    Area(Cmp, Expr),
    Plane(Axis, Expr, Option<Expr>),
    Normal(NormalSpec, Option<Expr>),
    Extreme(Axis, bool),
    InBox([Expr; 3], [Expr; 3], BoxMode),
    Adjacent(Box<Selector>),
    Id(i32),
    Union(Box<Selector>, Box<Selector>),
    Inter(Box<Selector>, Box<Selector>),
    Diff(Box<Selector>, Box<Selector>),
    Not(Box<Selector>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EntityKind {
    Regions,
    Faces,
    Edges,
    Vertices,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EntitySet {
    pub kind: EntityKind,
    /// sorted, unique
    pub ids: Vec<u32>,
}

impl EntitySet {
    pub fn new(kind: EntityKind, mut ids: Vec<u32>) -> EntitySet {
        ids.sort_unstable();
        ids.dedup();
        EntitySet { kind, ids }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Resolved {
    pub set: EntitySet,
    pub warnings: Vec<String>,
}

// ------------------------------------------------------------------------------------------------ text form

struct Parser<'a> {
    s: &'a str,
    i: usize,
}

fn bad<T>(msg: impl Into<String>) -> Result<T, GeomError> {
    Err(GeomError(msg.into()))
}

fn parse_expr(s: &str) -> Expr {
    let t = s.trim();
    match t.parse::<f64>() {
        Ok(v)
            if v.is_finite()
                && !t.is_empty()
                && !t.eq_ignore_ascii_case("nan")
                && !t.to_ascii_lowercase().contains("inf") =>
        {
            Expr::Num(v)
        }
        _ => Expr::Text(t.to_string()),
    }
}

fn parse_vec3(s: &str) -> GeomResult<[Expr; 3]> {
    let t = s.trim();
    let inner = t
        .strip_prefix('[')
        .and_then(|x| x.strip_suffix(']'))
        .ok_or_else(|| GeomError(format!("expected a vector like [x, y, z], found '{t}'")))?;
    let parts = split_top(inner);
    if parts.len() != 3 {
        return bad(format!("a vector needs three components, found '{t}'"));
    }
    Ok([parse_expr(&parts[0]), parse_expr(&parts[1]), parse_expr(&parts[2])])
}

/// Splits at top-level commas (not inside () or []).
fn split_top(s: &str) -> Vec<String> {
    let (mut depth, mut cur, mut out) = (0i32, String::new(), vec![]);
    for c in s.chars() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            _ => {}
        }
        if c == ',' && depth == 0 {
            out.push(cur.trim().to_string());
            cur.clear();
        } else {
            cur.push(c);
        }
    }
    if !cur.trim().is_empty() || !out.is_empty() {
        out.push(cur.trim().to_string());
    }
    out
}

fn axis_of(s: &str) -> GeomResult<Axis> {
    match s.trim() {
        "x" => Ok(Axis::X),
        "y" => Ok(Axis::Y),
        "z" => Ok(Axis::Z),
        o => bad(format!("expected an axis x, y or z, found '{o}'")),
    }
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.s[self.i..].starts_with(char::is_whitespace) {
            self.i += self.s[self.i..].chars().next().map_or(1, char::len_utf8);
        }
    }
    fn peek(&mut self) -> Option<char> {
        self.ws();
        self.s[self.i..].chars().next()
    }
    fn sel(&mut self) -> GeomResult<Selector> {
        let mut l = self.unary()?;
        while let Some(c) = self.peek() {
            let ctor: fn(Box<Selector>, Box<Selector>) -> Selector = match c {
                '|' => Selector::Union,
                '&' => Selector::Inter,
                '-' => Selector::Diff,
                _ => break,
            };
            self.i += 1;
            let r = self.unary()?;
            l = ctor(Box::new(l), Box::new(r));
        }
        Ok(l)
    }
    fn unary(&mut self) -> GeomResult<Selector> {
        match self.peek() {
            Some('!') => {
                self.i += 1;
                Ok(Selector::Not(Box::new(self.unary()?)))
            }
            Some('(') => {
                self.i += 1;
                let s = self.sel()?;
                if self.peek() != Some(')') {
                    return bad("missing ')' in selector");
                }
                self.i += 1;
                Ok(s)
            }
            Some(_) => self.atom(),
            None => bad("a selector is missing after an operator"),
        }
    }
    fn atom(&mut self) -> GeomResult<Selector> {
        let st = self.i;
        while self.s[self.i..].starts_with(|c: char| c.is_ascii_alphanumeric() || c == '_') {
            self.i += 1;
        }
        let word = &self.s[st..self.i];
        if word.is_empty() {
            return bad(format!(
                "unexpected '{}' in selector at position {}",
                self.s[self.i..].chars().next().unwrap_or('?'),
                self.i + 1
            ));
        }
        if word == "outer" && self.peek() != Some('(') {
            return Ok(Selector::Outer);
        }
        if self.s[self.i..].starts_with(":edge")
            && !self.s[self.i + 5..].starts_with(|c: char| c.is_ascii_alphanumeric() || c == '_')
        {
            self.i += 5;
            return Ok(Selector::Name(format!("{word}:edge")));
        }
        if self.peek() != Some('(') {
            return bad(format!("'{word}' needs arguments in parentheses, like {word}(...)"));
        }
        self.i += 1;
        let (st, mut depth) = (self.i, 1);
        while self.i < self.s.len() {
            match self.s.as_bytes()[self.i] {
                b'(' | b'[' => depth += 1,
                b')' | b']' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                _ => {}
            }
            self.i += 1;
        }
        if depth != 0 {
            return bad(format!("missing ')' after '{word}('"));
        }
        let raw = &self.s[st..self.i];
        self.i += 1;
        let args = split_top(raw);
        let n = args.len();
        let need = |lo: usize, hi: usize| -> GeomResult<()> {
            if n < lo || n > hi {
                bad(format!(
                    "{word}(...) takes {} argument(s), found {n}",
                    if lo == hi { lo.to_string() } else { format!("{lo} to {hi}") }
                ))
            } else {
                Ok(())
            }
        };
        Ok(match word {
            "name" => {
                need(1, 1)?;
                Selector::Name(args[0].clone())
            }
            "tag" => {
                need(1, 1)?;
                Selector::Tag(args[0].clone())
            }
            "region" => {
                need(1, 1)?;
                Selector::Region(args[0].clone())
            }
            "boundary" | "adjacent" => {
                need(1, 1)?;
                let inner = parse(&args[0])?;
                if word == "boundary" {
                    Selector::Boundary(Box::new(inner))
                } else {
                    Selector::Adjacent(Box::new(inner))
                }
            }
            "interface" => {
                need(2, 2)?;
                Selector::Interface(args[0].clone(), if args[1] == "*" { None } else { Some(args[1].clone()) })
            }
            "from" => {
                need(2, 2)?;
                Selector::From(args[0].clone(), args[1].clone())
            }
            "kind" => {
                need(1, 1)?;
                Selector::Kind(match args[0].as_str() {
                    "plane" => SurfKind::Plane,
                    "cylinder" => SurfKind::Cylinder,
                    "cone" => SurfKind::Cone,
                    "sphere" => SurfKind::Sphere,
                    "torus" => SurfKind::Torus,
                    "line" => SurfKind::Line,
                    "circle" => SurfKind::Circle,
                    "ellipse" => SurfKind::Ellipse,
                    o => {
                        return bad(format!(
                            "unknown kind '{o}' (plane, cylinder, cone, sphere, torus, line, circle, ellipse)"
                        ))
                    }
                })
            }
            "radius" => {
                need(1, 2)?;
                Selector::Radius(parse_expr(&args[0]), args.get(1).map(|a| parse_expr(a)))
            }
            "area" => {
                need(1, 1)?;
                let a = &args[0];
                let (cmp, rest) = if let Some(r) = a.strip_prefix("<=") {
                    (Cmp::Le, r)
                } else if let Some(r) = a.strip_prefix(">=") {
                    (Cmp::Ge, r)
                } else if let Some(r) = a.strip_prefix('<') {
                    (Cmp::Lt, r)
                } else if let Some(r) = a.strip_prefix('>') {
                    (Cmp::Gt, r)
                } else {
                    return bad(format!("area(...) needs a comparison like area(>1e-12), found '{a}'"));
                };
                Selector::Area(cmp, parse_expr(rest))
            }
            "plane" => {
                need(1, 2)?;
                let (ax, val) = args[0]
                    .split_once('=')
                    .ok_or_else(|| GeomError(format!("plane(...) needs 'axis = value', found '{}'", args[0])))?;
                Selector::Plane(axis_of(ax)?, parse_expr(val), args.get(1).map(|a| parse_expr(a)))
            }
            "normal" => {
                need(1, 2)?;
                let a = args[0].as_str();
                let spec = if a.starts_with('[') {
                    NormalSpec::Vector(parse_vec3(a)?)
                } else {
                    let (positive, ax) = match a.as_bytes().first() {
                        Some(b'+') => (true, &a[1..]),
                        Some(b'-') => (false, &a[1..]),
                        _ => (true, a),
                    };
                    NormalSpec::Axis { axis: axis_of(ax)?, positive }
                };
                Selector::Normal(spec, args.get(1).map(|a| parse_expr(a)))
            }
            "extreme" => {
                need(1, 1)?;
                let a = args[0].as_str();
                let (positive, ax) = match a.as_bytes().first() {
                    Some(b'+') => (true, &a[1..]),
                    Some(b'-') => (false, &a[1..]),
                    _ => (true, a),
                };
                Selector::Extreme(axis_of(ax)?, positive)
            }
            "in_box" => {
                need(2, 3)?;
                let mode = match args.get(2).map(String::as_str) {
                    None | Some("within") => BoxMode::Within,
                    Some("touching") => BoxMode::Touching,
                    Some(o) => return bad(format!("in_box mode must be within or touching, found '{o}'")),
                };
                Selector::InBox(parse_vec3(&args[0])?, parse_vec3(&args[1])?, mode)
            }
            "id" => {
                need(1, 1)?;
                Selector::Id(
                    args[0]
                        .parse::<i32>()
                        .map_err(|_| GeomError(format!("id(...) needs an integer, found '{}'", args[0])))?,
                )
            }
            o => return bad(format!("unknown selector term '{o}'")),
        })
    }
}

fn parse(s: &str) -> GeomResult<Selector> {
    let mut p = Parser { s, i: 0 };
    let r = p.sel()?;
    if p.peek().is_some() {
        return bad(format!("unexpected '{}' in selector at position {}", &s[p.i..p.i + 1], p.i + 1));
    }
    Ok(r)
}

fn pe(e: &Expr) -> String {
    e.source()
}

fn pv(v: &[Expr; 3]) -> String {
    format!("[{}, {}, {}]", pe(&v[0]), pe(&v[1]), pe(&v[2]))
}

fn sign_axis(a: Axis, positive: bool) -> String {
    format!("{}{}", if positive { "+" } else { "-" }, a.name())
}

impl std::fmt::Display for Selector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bin = |f: &mut std::fmt::Formatter<'_>, a: &Selector, op: &str, b: &Selector| {
            write!(f, "{a} {op} ")?;
            if matches!(b, Selector::Union(..) | Selector::Inter(..) | Selector::Diff(..)) {
                write!(f, "({b})")
            } else {
                write!(f, "{b}")
            }
        };
        match self {
            Selector::Name(n) => write!(f, "name({n})"),
            Selector::Tag(n) => write!(f, "tag({n})"),
            Selector::Region(n) => write!(f, "region({n})"),
            Selector::Boundary(s) => write!(f, "boundary({s})"),
            Selector::Outer => write!(f, "outer"),
            Selector::Interface(a, b) => write!(f, "interface({a}, {})", b.as_deref().unwrap_or("*")),
            Selector::From(a, b) => write!(f, "from({a}, {b})"),
            Selector::Kind(k) => write!(f, "kind({})", format!("{k:?}").to_lowercase()),
            Selector::Radius(r, t) => match t {
                Some(t) => write!(f, "radius({}, {})", pe(r), pe(t)),
                None => write!(f, "radius({})", pe(r)),
            },
            Selector::Area(c, v) => {
                let op = match c {
                    Cmp::Lt => "<",
                    Cmp::Le => "<=",
                    Cmp::Gt => ">",
                    Cmp::Ge => ">=",
                };
                write!(f, "area({op}{})", pe(v))
            }
            Selector::Plane(a, v, t) => match t {
                Some(t) => write!(f, "plane({} = {}, {})", a.name(), pe(v), pe(t)),
                None => write!(f, "plane({} = {})", a.name(), pe(v)),
            },
            Selector::Normal(n, t) => {
                let ns = match n {
                    NormalSpec::Axis { axis, positive } => sign_axis(*axis, *positive),
                    NormalSpec::Vector(v) => pv(v),
                };
                match t {
                    Some(t) => write!(f, "normal({ns}, {})", pe(t)),
                    None => write!(f, "normal({ns})"),
                }
            }
            Selector::Extreme(a, p) => write!(f, "extreme({})", sign_axis(*a, *p)),
            Selector::InBox(lo, hi, m) => {
                write!(
                    f,
                    "in_box({}, {}, {})",
                    pv(lo),
                    pv(hi),
                    if *m == BoxMode::Within { "within" } else { "touching" }
                )
            }
            Selector::Adjacent(s) => write!(f, "adjacent({s})"),
            Selector::Id(i) => write!(f, "id({i})"),
            Selector::Union(a, b) => bin(f, a, "|", b),
            Selector::Inter(a, b) => bin(f, a, "&", b),
            Selector::Diff(a, b) => bin(f, a, "-", b),
            Selector::Not(s) => {
                if matches!(**s, Selector::Union(..) | Selector::Inter(..) | Selector::Diff(..)) {
                    write!(f, "!({s})")
                } else {
                    write!(f, "!{s}")
                }
            }
        }
    }
}

impl Selector {
    pub fn parse(text: &str) -> GeomResult<Selector> {
        parse(text)
    }
    pub fn atoms(&self) -> usize {
        match self {
            Selector::Union(a, b) | Selector::Inter(a, b) | Selector::Diff(a, b) => a.atoms() + b.atoms(),
            Selector::Not(a) | Selector::Boundary(a) | Selector::Adjacent(a) => 1 + a.atoms(),
            _ => 1,
        }
    }
    pub fn or(self, o: Selector) -> Selector {
        Selector::Union(Box::new(self), Box::new(o))
    }
    pub fn and(self, o: Selector) -> Selector {
        Selector::Inter(Box::new(self), Box::new(o))
    }
    pub fn minus(self, o: Selector) -> Selector {
        Selector::Diff(Box::new(self), Box::new(o))
    }
}

// ------------------------------------------------------------------------------------------------ resolution

struct Ctx<'a> {
    ev: &'a Evaluated,
    warnings: Vec<String>,
    depth: usize,
    no_id: bool,
}

fn num(ev: &Evaluated, e: &Expr, what: &str) -> GeomResult<f64> {
    ev.scene().env.eval(e, what)
}

fn norm3v(v: V3) -> V3 {
    unit3(v)
}

impl Ctx<'_> {
    fn facets(&self) -> EntityKind {
        if self.ev.list.dim == 2 {
            EntityKind::Edges
        } else {
            EntityKind::Faces
        }
    }

    fn empty_warn(&mut self, atom: &Selector, set: EntitySet) -> EntitySet {
        if set.ids.is_empty() {
            self.warnings.push(format!("'{atom}' matches nothing"));
        }
        set
    }

    fn all_facets(&self) -> Vec<u32> {
        (0..self.ev.complex().facet_count() as u32).collect()
    }

    fn all_of_kind(&self, kind: EntityKind) -> EntitySet {
        let cx = self.ev.complex();
        let n = match kind {
            EntityKind::Regions => cx.regions.len(),
            EntityKind::Faces => cx.faces.len(),
            EntityKind::Edges => {
                if self.ev.list.dim == 2 {
                    cx.edges.len()
                } else {
                    cx.edges3.len()
                }
            }
            EntityKind::Vertices => cx.vertices.len(),
        } as u32;
        EntitySet::new(kind, (0..n).collect())
    }

    fn is_polymorphic_box(s: &Selector) -> bool {
        match s {
            Selector::InBox(..) => true,
            Selector::Not(a) => Self::is_polymorphic_box(a),
            Selector::Union(a, b) | Selector::Inter(a, b) | Selector::Diff(a, b) => {
                Self::is_polymorphic_box(a) && Self::is_polymorphic_box(b)
            }
            _ => false,
        }
    }

    fn region_id(&self, name: &str) -> Option<u32> {
        self.ev.complex().regions.iter().position(|r| r.name == name).map(|i| i as u32)
    }

    fn filter(&self, universe: Option<&EntitySet>, f: impl Fn(u32) -> bool) -> Vec<u32> {
        match universe {
            Some(u) if u.kind == self.facets() => u.ids.iter().copied().filter(|&i| f(i)).collect(),
            _ => self.all_facets().into_iter().filter(|&i| f(i)).collect(),
        }
    }

    fn gtol(&self) -> f64 {
        1e-7 * self.ev.bounds().diagonal()
    }

    fn feature_tags(&mut self, tag: &str, out: &mut Vec<EntitySet>) -> GeomResult<()> {
        let list = &self.ev.list;
        let cx = self.ev.complex();
        // tags on features: all facets that descend from the feature (through transforms, patterns and Booleans)
        let mut feats: Vec<u32> = vec![];
        for f in list.features.iter().filter(|f| f.tags.iter().any(|t| t == tag)) {
            let mut stack = vec![f.id];
            while let Some(id) = stack.pop() {
                feats.push(id);
                if let Some(ff) = list.feature(id) {
                    stack.extend(ff.op.inputs());
                }
            }
        }
        for (fid, role, t) in &self.ev.scene().auto_tags {
            if t == tag {
                let ids: Vec<u32> = (0..cx.facet_count() as u32)
                    .filter(|&i| cx.facet_prov(i).iter().any(|p| p.feature == *fid && &p.role == role))
                    .collect();
                out.push(EntitySet::new(self.facets(), ids));
            }
        }
        if !feats.is_empty() {
            let ids: Vec<u32> = (0..cx.facet_count() as u32)
                .filter(|&i| cx.facet_prov(i).iter().any(|p| feats.contains(&p.feature)))
                .collect();
            out.push(EntitySet::new(self.facets(), ids));
        }
        Ok(())
    }

    fn go(&mut self, s: &Selector, universe: Option<&EntitySet>) -> GeomResult<EntitySet> {
        let cx = self.ev.complex();
        let facets = self.facets();
        Ok(match s {
            Selector::Name(n) => {
                if let Some(sheet) = n.strip_suffix(":edge") {
                    let es = cx.sheet_perimeter_edges(sheet);
                    return Ok(self.empty_warn(s, EntitySet::new(EntityKind::Edges, es)));
                }
                let ids: Vec<u32> = self.all_facets().into_iter().filter(|&i| cx.facet_name(i) == n).collect();
                if ids.is_empty() {
                    let es: Vec<u32> =
                        (0..cx.edges3.len() as u32).filter(|&e| cx.edges3[e as usize].name == *n).collect();
                    if !es.is_empty() {
                        return Ok(EntitySet::new(EntityKind::Edges, es));
                    }
                    let vs: Vec<u32> =
                        (0..cx.vertices.len() as u32).filter(|&v| cx.vertices[v as usize].name == *n).collect();
                    self.empty_warn(s, EntitySet::new(if vs.is_empty() { facets } else { EntityKind::Vertices }, vs))
                } else {
                    EntitySet::new(facets, ids)
                }
            }
            Selector::Tag(t) => {
                let mut sets = vec![];
                self.feature_tags(t, &mut sets)?;
                let list = &self.ev.list;
                if self.depth > 8 {
                    return bad(format!("tag '{t}' is defined in terms of itself"));
                }
                for f in &list.features {
                    if let crate::model::Op::Tag { sel, name } = &f.op {
                        if name == t && !f.suppressed {
                            self.depth += 1;
                            let r = self.go(sel, None);
                            self.depth -= 1;
                            sets.push(r?);
                        }
                    }
                }
                let kind = sets.first().map_or(facets, |x| x.kind);
                if sets.iter().any(|x| x.kind != kind) {
                    return bad(format!("tag '{t}' is attached to entities of different kinds"));
                }
                let set = EntitySet::new(kind, sets.into_iter().flat_map(|x| x.ids).collect());
                self.empty_warn(s, set)
            }
            Selector::Region(n) => {
                let ids = self.region_id(n).into_iter().collect();
                self.empty_warn(s, EntitySet::new(EntityKind::Regions, ids))
            }
            Selector::Boundary(inner) => {
                let r = self.go(inner, None)?;
                if r.kind != EntityKind::Regions {
                    return bad(format!("boundary(...) needs regions, but '{inner}' selects {:?}", r.kind));
                }
                let ids = self.all_facets().into_iter().filter(|&i| {
                    let (a, b) = cx.facet_sides(i);
                    a.is_some_and(|x| r.ids.contains(&x)) || b.is_some_and(|x| r.ids.contains(&x))
                });
                let set = EntitySet::new(facets, ids.collect());
                self.empty_warn(s, set)
            }
            Selector::Outer => {
                let ids = self.filter(universe, |i| {
                    let (a, b) = cx.facet_sides(i);
                    a.is_none() != b.is_none()
                });
                EntitySet::new(facets, ids)
            }
            Selector::Interface(a, b) => {
                let (ra, rb) = (self.region_id(a), b.as_ref().map(|b| self.region_id(b)));
                if ra.is_none() || matches!(rb, Some(None)) {
                    self.warnings.push(format!("'{s}': no region with that name"));
                    return Ok(EntitySet::new(facets, vec![]));
                }
                let ra = ra.unwrap();
                let ids = self.filter(None, |i| {
                    let (l, r) = cx.facet_sides(i);
                    match rb {
                        None => l == Some(ra) || r == Some(ra),
                        Some(rb) => (l == Some(ra) && r == rb) || (r == Some(ra) && l == rb),
                    }
                });
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::From(f, role) => {
                let list = &self.ev.list;
                let ids = self.filter(None, |i| {
                    cx.facet_prov(i).iter().any(|p| {
                        &p.role == role && (&p.instance == f || list.feature(p.feature).is_some_and(|ff| &ff.name == f))
                    })
                });
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Kind(k) => {
                let want = format!("{k:?}").to_lowercase();
                let ids = self.filter(universe, |i| cx.facet_kind(i) == want);
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Radius(r, t) => {
                let rv = num(self.ev, r, "radius()")?;
                let tv = match t {
                    Some(t) => num(self.ev, t, "radius() tolerance")?,
                    None => 1e-6 * rv.abs(),
                };
                let ids = self.filter(universe, |i| cx.facet_radius(i).is_some_and(|x| (x - rv).abs() <= tv));
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Area(c, v) => {
                let v = num(self.ev, v, "area()")?;
                let ids = self.filter(universe, |i| {
                    let m = cx.facet_measure(i);
                    match c {
                        Cmp::Lt => m < v,
                        Cmp::Le => m <= v,
                        Cmp::Gt => m > v,
                        Cmp::Ge => m >= v,
                    }
                });
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Plane(ax, v, t) => {
                let (v, t) = (
                    num(self.ev, v, "plane()")?,
                    match t {
                        Some(t) => num(self.ev, t, "plane() tolerance")?,
                        None => self.gtol(),
                    },
                );
                let k = ax.index();
                let ids = self.filter(universe, |i| {
                    let b = cx.facet_bbox(i);
                    (b.lo[k] - v).abs() <= t && (b.hi[k] - v).abs() <= t
                });
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Normal(n, t) => {
                let dir = match n {
                    NormalSpec::Axis { axis, positive } => {
                        let mut d = [0.0; 3];
                        d[axis.index()] = if *positive { 1.0 } else { -1.0 };
                        d
                    }
                    NormalSpec::Vector(v) => norm3v([
                        num(self.ev, &v[0], "normal()")?,
                        num(self.ev, &v[1], "normal()")?,
                        num(self.ev, &v[2], "normal()")?,
                    ]),
                };
                let ang = match t {
                    Some(t) => num(self.ev, t, "normal() tolerance")?,
                    None => 0.5,
                };
                let cos_tol = ang.to_radians().cos();
                let ids = self.filter(universe, |i| {
                    let interface = {
                        let (a, b) = cx.facet_sides(i);
                        a.is_some() && b.is_some()
                    };
                    cx.facet_normals(i).iter().all(|n| {
                        let d = dot3(*n, dir);
                        d >= cos_tol || (interface && -d >= cos_tol)
                    })
                });
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Extreme(ax, pos) => {
                let k = ax.index();
                let t = self.gtol();
                let pool = self.filter(universe, |_| true);
                let val = |i: u32| if *pos { cx.facet_bbox(i).hi[k] } else { cx.facet_bbox(i).lo[k] };
                let best =
                    pool.iter().map(|&i| val(i)).fold(if *pos { f64::NEG_INFINITY } else { f64::INFINITY }, |a, b| {
                        if *pos {
                            a.max(b)
                        } else {
                            a.min(b)
                        }
                    });
                let ids = pool.into_iter().filter(|&i| {
                    let b = cx.facet_bbox(i);
                    (val(i) - best).abs() <= t && b.hi[k] - b.lo[k] <= t
                });
                self.empty_warn(s, EntitySet::new(facets, ids.collect()))
            }
            Selector::InBox(lo, hi, mode) => {
                let ev3 = |v: &[Expr; 3]| -> GeomResult<V3> {
                    Ok([
                        num(self.ev, &v[0], "in_box()")?,
                        num(self.ev, &v[1], "in_box()")?,
                        num(self.ev, &v[2], "in_box()")?,
                    ])
                };
                let (a, b) = (ev3(lo)?, ev3(hi)?);
                let (lo, hi) = (
                    [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])],
                    [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])],
                );
                let t = self.gtol();
                let bx = Aabb::new(lo, hi).pad(t);
                let two_d = self.ev.list.dim == 2;
                let kind = universe.map_or(facets, |u| u.kind);
                let pool: Vec<u32> = match universe {
                    Some(u) if u.kind == kind => u.ids.clone(),
                    _ => self.all_of_kind(kind).ids,
                };
                let matches_box = |mut b: Aabb| -> bool {
                    if two_d {
                        // y is out of the x-z cross-section: ignore it
                        b.lo[1] = lo[1];
                        b.hi[1] = hi[1];
                    }
                    match mode {
                        BoxMode::Within => (0..3).all(|k| b.lo[k] >= bx.lo[k] && b.hi[k] <= bx.hi[k]),
                        BoxMode::Touching => b.overlaps(&bx),
                    }
                };
                let ids: Vec<u32> = pool
                    .into_iter()
                    .filter(|&i| match kind {
                        k if k == facets => matches_box(cx.facet_bbox(i)),
                        EntityKind::Edges => {
                            let e = &cx.edges3[i as usize];
                            let p0 = cx.vertices[e.v0 as usize].p;
                            let p1 = cx.vertices[e.v1 as usize].p;
                            let pm = cx.edge3_point(i, 0.5);
                            let mut b = e.bbox;
                            b.grow(p0);
                            b.grow(p1);
                            b.grow(pm);
                            match mode {
                                BoxMode::Within => {
                                    (0..3).all(|k| b.lo[k] >= bx.lo[k] && b.hi[k] <= bx.hi[k])
                                        && bx.contains(p0)
                                        && bx.contains(p1)
                                        && bx.contains(pm)
                                }
                                BoxMode::Touching => {
                                    b.overlaps(&bx)
                                        && (bx.contains(p0)
                                            || bx.contains(p1)
                                            || bx.contains(pm)
                                            || bx.contains(cx.edge3_closest(i, bx.center()).0))
                                }
                            }
                        }
                        EntityKind::Regions => matches_box(cx.regions[i as usize].bbox),
                        EntityKind::Vertices => {
                            let mut p = cx.vertices[i as usize].p;
                            if two_d {
                                p[1] = lo[1];
                            }
                            bx.contains(p)
                        }
                        _ => false,
                    })
                    .collect();
                self.empty_warn(s, EntitySet::new(kind, ids))
            }
            Selector::Adjacent(inner) => {
                let r = self.go(inner, None)?;
                match r.kind {
                    EntityKind::Regions => {
                        let ids = (0..cx.regions.len() as u32).filter(|x| !r.ids.contains(x)).filter(|&x| {
                            self.all_facets().into_iter().any(|i| {
                                let (a, b) = cx.facet_sides(i);
                                (a == Some(x) && b.is_some_and(|b| r.ids.contains(&b)))
                                    || (b == Some(x) && a.is_some_and(|a| r.ids.contains(&a)))
                            })
                        });
                        EntitySet::new(EntityKind::Regions, ids.collect())
                    }
                    k if k == facets => {
                        let vs: Vec<u32> = r.ids.iter().flat_map(|&i| cx.facet_vertices(i)).collect();
                        let ids = self
                            .all_facets()
                            .into_iter()
                            .filter(|i| !r.ids.contains(i) && cx.facet_vertices(*i).iter().any(|v| vs.contains(v)));
                        EntitySet::new(facets, ids.collect())
                    }
                    k => {
                        return bad(format!(
                            "adjacent(...) is defined for regions and {facets:?}, but '{inner}' selects {k:?}"
                        ))
                    }
                }
            }
            Selector::Id(n) => {
                if self.no_id {
                    return bad("id(...) cannot be used inside a label pin");
                }
                let lm = crate::labels::compile(self.ev)?;
                let ids = lm.edges_of_label(*n);
                self.empty_warn(s, EntitySet::new(facets, ids))
            }
            Selector::Union(a, b) | Selector::Inter(a, b) | Selector::Diff(a, b) => {
                let (poly_a, poly_b) = (Self::is_polymorphic_box(a), Self::is_polymorphic_box(b));
                let (x, y) = match s {
                    Selector::Inter(..) => {
                        if poly_a && !poly_b {
                            let y = self.go(b, universe)?;
                            let x = self.go(a, Some(&y))?;
                            (x, y)
                        } else {
                            let x = self.go(a, universe)?;
                            let y = self.go(b, Some(&x))?;
                            (x, y)
                        }
                    }
                    Selector::Diff(..) => {
                        if poly_a && !poly_b {
                            let y = self.go(b, universe)?;
                            let full = self.all_of_kind(y.kind);
                            let x = self.go(a, Some(universe.unwrap_or(&full)))?;
                            (x, y)
                        } else {
                            let x = self.go(a, universe)?;
                            let y = if poly_b { self.go(b, Some(&x))? } else { self.go(b, universe)? };
                            (x, y)
                        }
                    }
                    Selector::Union(..) => {
                        if poly_a && !poly_b {
                            let y = self.go(b, universe)?;
                            let full = self.all_of_kind(y.kind);
                            let x = self.go(a, Some(universe.unwrap_or(&full)))?;
                            (x, y)
                        } else if !poly_a && poly_b {
                            let x = self.go(a, universe)?;
                            let full = self.all_of_kind(x.kind);
                            let y = self.go(b, Some(universe.unwrap_or(&full)))?;
                            (x, y)
                        } else {
                            (self.go(a, universe)?, self.go(b, universe)?)
                        }
                    }
                    _ => unreachable!(),
                };
                if x.kind != y.kind && !x.ids.is_empty() && !y.ids.is_empty() {
                    return bad(format!(
                        "cannot combine {:?} ('{a}') with {:?} ('{b}') in one selector",
                        x.kind, y.kind
                    ));
                }
                let kind = if x.ids.is_empty() { y.kind } else { x.kind };
                let ids = match s {
                    Selector::Union(..) => x.ids.iter().chain(&y.ids).copied().collect(),
                    Selector::Inter(..) => x.ids.iter().copied().filter(|i| y.ids.contains(i)).collect(),
                    _ => x.ids.iter().copied().filter(|i| !y.ids.contains(i)).collect(),
                };
                EntitySet::new(kind, ids)
            }
            Selector::Not(a) => {
                let target_kind = if Self::is_polymorphic_box(a) { universe.map(|u| u.kind) } else { None };
                let full_opt = target_kind.map(|k| self.all_of_kind(k));
                let x = self.go(a, full_opt.as_ref())?;
                let full = full_opt.unwrap_or_else(|| self.all_of_kind(x.kind));
                EntitySet::new(x.kind, full.ids.into_iter().filter(|i| !x.ids.contains(i)).collect())
            }
        })
    }
}

impl Evaluated {
    /// Resolves a selector to an entity set plus warnings (an empty result is a warning, not an error).
    pub fn resolve(&self, sel: &Selector) -> GeomResult<Resolved> {
        self.resolve_inner(sel, false)
    }

    pub(crate) fn resolve_inner(&self, sel: &Selector, no_id: bool) -> GeomResult<Resolved> {
        let mut c = Ctx { ev: self, warnings: vec![], depth: 0, no_id };
        let set = c.go(sel, None)?;
        let mut warnings = c.warnings;
        if set.ids.is_empty() && warnings.is_empty() {
            warnings.push(format!("'{sel}' matches nothing"));
        }
        Ok(Resolved { set, warnings })
    }

    /// Selectors, simplest first, that each resolve to exactly `picked`; empty if none does.
    pub fn suggest(&self, picked: &EntitySet) -> Vec<Selector> {
        let cx = self.complex();
        let mut cands: Vec<Selector> = vec![];
        let facets = if self.list.dim == 2 { EntityKind::Edges } else { EntityKind::Faces };
        if picked.ids.is_empty() {
            return vec![];
        }
        match picked.kind {
            EntityKind::Regions => {
                let names: Vec<&str> = picked.ids.iter().map(|&i| cx.regions[i as usize].name.as_str()).collect();
                let mut it = names.iter().map(|n| Selector::Region(n.to_string()));
                let first = it.next().unwrap();
                cands.push(it.fold(first, |a, b| a.or(b)));
            }
            EntityKind::Vertices => {
                let mut names: Vec<&str> = picked.ids.iter().map(|&i| cx.vertices[i as usize].name.as_str()).collect();
                names.sort_unstable();
                names.dedup();
                cands.extend(names.iter().map(|n| Selector::Name(n.to_string())));
            }
            k if k == facets => {
                let mut names: Vec<&str> = picked.ids.iter().map(|&i| cx.facet_name(i)).collect();
                names.sort_unstable();
                names.dedup();
                let tags: Vec<String> = {
                    let mut t: Vec<String> = self.list.features.iter().flat_map(|f| f.tags.iter().cloned()).collect();
                    t.extend(self.list.features.iter().filter_map(|f| {
                        if let crate::model::Op::Tag { name, .. } = &f.op {
                            Some(name.clone())
                        } else {
                            None
                        }
                    }));
                    t.extend(self.scene().auto_tags.iter().map(|x| x.2.clone()));
                    t.sort();
                    t.dedup();
                    t
                };
                cands.extend(tags.into_iter().map(Selector::Tag));
                cands.extend(names.iter().map(|n| Selector::Name(n.to_string())));
                if names.len() > 1 && names.len() <= 4 {
                    let mut it = names.iter().map(|n| Selector::Name(n.to_string()));
                    let first = it.next().unwrap();
                    cands.push(it.fold(first, |a, b| a.or(b)));
                }
                // region based
                let sides: Vec<(Option<u32>, Option<u32>)> = picked.ids.iter().map(|&i| cx.facet_sides(i)).collect();
                for (ri, r) in cx.regions.iter().enumerate() {
                    cands.push(Selector::Boundary(Box::new(Selector::Region(r.name.clone()))));
                    cands.push(Selector::Interface(r.name.clone(), None));
                    for (rj, o) in cx.regions.iter().enumerate().skip(ri + 1) {
                        cands.push(Selector::Interface(r.name.clone(), Some(o.name.clone())));
                        let _ = rj;
                    }
                    let _ = &sides;
                }
                cands.push(Selector::Outer);
                // from(feature, role)
                let mut froms: Vec<(String, String)> = picked
                    .ids
                    .iter()
                    .flat_map(|&i| {
                        cx.facet_prov(i).iter().map(|p| {
                            (
                                self.list.feature(p.feature).map_or(p.instance.clone(), |f| f.name.clone()),
                                p.role.clone(),
                            )
                        })
                    })
                    .collect();
                froms.sort();
                froms.dedup();
                cands.extend(froms.into_iter().map(|(f, r)| Selector::From(f, r)));
                // geometric predicates on the picked set's extent
                let mut bb = Aabb::EMPTY;
                for &i in &picked.ids {
                    bb = bb.union(&cx.facet_bbox(i));
                }
                let t = 1e-7 * self.bounds().diagonal();
                let bx = |m: BoxMode| {
                    Selector::InBox(
                        [0, 1, 2].map(|k| Expr::Num(bb.lo[k] - t)),
                        [0, 1, 2].map(|k| Expr::Num(bb.hi[k] + t)),
                        m,
                    )
                };
                let mut geom: Vec<Selector> = vec![bx(BoxMode::Within)];
                for ax in [Axis::X, Axis::Y, Axis::Z] {
                    let k = ax.index();
                    if bb.hi[k] - bb.lo[k] <= t && (self.list.dim == 3 || k != 1) {
                        geom.push(Selector::Plane(ax, Expr::Num(bb.lo[k]), None));
                    }
                    for pos in [true, false] {
                        geom.push(Selector::Normal(NormalSpec::Axis { axis: ax, positive: pos }, None));
                        geom.push(Selector::Extreme(ax, pos));
                    }
                }
                let base: Vec<Selector> = names.iter().map(|n| Selector::Name(n.to_string())).collect();
                for g in &geom {
                    cands.push(g.clone());
                    if base.len() == 1 {
                        cands.push(base[0].clone().and(g.clone()));
                    }
                    for r in &cx.regions {
                        cands.push(Selector::Boundary(Box::new(Selector::Region(r.name.clone()))).and(g.clone()));
                    }
                }
                for g in &geom {
                    for n in &names {
                        cands.push(Selector::Name(n.to_string()).and(g.clone()));
                    }
                }
                if picked.ids.len() == 1 {
                    if let Ok(lm) = crate::labels::compile(self) {
                        if let Some(l) = lm.label_of_edge(picked.ids[0]) {
                            cands.push(Selector::Id(l));
                        }
                    }
                }
            }
            _ => {}
        }
        let mut out: Vec<(usize, String, Selector)> = vec![];
        for c in cands {
            if let Ok(r) = self.resolve(&c) {
                if r.set == *picked {
                    let t = c.to_string();
                    if !out.iter().any(|o| o.1 == t) {
                        out.push((c.atoms(), t, c));
                    }
                }
            }
        }
        out.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        out.into_iter().map(|o| o.2).collect()
    }
}
