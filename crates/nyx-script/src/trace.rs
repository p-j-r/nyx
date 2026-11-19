//! Tracing of script functions (productivity plan S6, N3, ADR 90): a Python function given where a formula is expected is called once with a
//! symbolic point and becomes formula text for the engine.
//!
//! * `prelude`: the Python text fed to the worker before the script: `prelude.nys` (the symbolic value, the point, the tracer) plus the
//!   functions of the formula grammar (`where`, `sqrt`, ...) and the `nyx` object, which wraps every engine command so that a function (or a
//!   symbolic value) in a formula argument is replaced by its text. Written from `Engine::catalog()` and the lists below, so the stub and the
//!   prelude cannot disagree.
//! * `stub`: the same names as Python declarations, for the static check and for an editor.
//! * `guard`: the static check that keeps the plan's rule. Monty 1.6.0 does not call `__bool__`, `__lt__` and its siblings, `__ne__` or
//!   `__getattr__` on an instance, so `if sym:` would silently take the true branch. The guard parses the script, finds the functions passed to
//!   a `nyx.` call, follows which names depend on the point and refuses a branch or a comparison on one, naming `where(...)`.

use std::collections::{HashMap, HashSet};

use nyx_engine::Engine;
use ruff_python_ast::visitor::{self, Visitor};
use ruff_python_ast::{Comprehension, Expr, Parameters, Stmt, StmtFunctionDef, UnaryOp};
use ruff_text_size::Ranged;
use serde_json::Value;

/// The arguments of the engine commands that take a formula: a function (or a symbolic value) given for one of them is traced.
pub const FORMULA_ARGS: &[&str] = &["expr", "a", "b", "fields", "weight"];
/// One value in, the formula function of the same name out (`abs` is Python's own, through `__abs__`).
pub const FUNCS_1: &[&str] = &["sqrt", "exp", "sin", "cos", "arg", "log", "log10", "tan", "sinh", "cosh", "tanh", "atan", "step", "re", "im", "conj", "abs2", "norm", "norm2"];
/// Two values in: the name in a script, the name in the formula.
pub const FUNCS_2: &[(&str, &str)] = &[("atan2", "atan2"), ("minimum", "min"), ("maximum", "max"), ("dot", "dot"), ("inner", "inner"), ("cross", "cross")];
/// Conditions for `where`: the name in a script, the operator in the formula.
pub const CONDS: &[(&str, &str)] = &[("lt", "<"), ("gt", ">"), ("le", "<="), ("ge", ">="), ("eq", "=="), ("ne", "!="), ("land", "and"), ("lor", "or")];

const NO_BRANCH: &str = "write where(cond, a, b) with the conditions lt, gt, le, ge, eq, ne, land, lor, lnot (a function given as a formula runs once, not at every point)";

/// Names a script gets besides `nyx`: for the stub.
fn all_names() -> Vec<(String, usize)> {
    let mut v: Vec<(String, usize)> = vec![("where".into(), 3), ("material".into(), 1), ("region".into(), 1), ("trace".into(), 1), ("lnot".into(), 1)];
    v.extend(FUNCS_1.iter().map(|n| (n.to_string(), 1)));
    v.extend(FUNCS_2.iter().map(|(n, _)| (n.to_string(), 2)));
    v.extend(CONDS.iter().map(|(n, _)| (n.to_string(), 2)));
    v
}

/// The Python text fed to the worker before the script, with `nyx` (the engine's module) as its input: it leaves `nyx` bound to a wrapper.
pub fn prelude(engine: &Engine) -> String {
    let mut s = String::from(include_str!("prelude.nys"));
    s.push_str("\n\n");
    for n in FUNCS_1 {
        s.push_str(&format!("def {n}(x):\n    return _nyx_fn(\"{n}\", x)\n\n\n"));
    }
    for (n, f) in FUNCS_2 {
        s.push_str(&format!("def {n}(a, b):\n    return _nyx_fn(\"{f}\", a, b)\n\n\n"));
    }
    for (n, op) in CONDS {
        s.push_str(&format!("def {n}(a, b):\n    return Sym(\"(\" + _nyx_t(a) + \" {op} \" + _nyx_t(b) + \")\")\n\n\n"));
    }
    s.push_str("def lnot(a):\n    return Sym(\"(not \" + _nyx_t(a) + \")\")\n\n\n");
    s.push_str("class _Nyx:\n");
    s.push_str("    def call(self, name, *a, **kw):\n        return _nyx_host.call(name, *_nyx_pos(a), **_nyx_kw(kw))\n\n");
    for c in engine.catalog() {
        let takes_formula = c
            .schema
            .get("properties")
            .and_then(Value::as_object)
            .is_some_and(|p| FORMULA_ARGS.iter().any(|a| p.contains_key(*a)));
        if takes_formula {
            s.push_str(&format!("    def {n}(self, *a, **kw):\n        return _nyx_host.{n}(*_nyx_pos(a), **_nyx_kw(kw))\n\n", n = c.name));
        } else {
            s.push_str(&format!("    def {n}(self, *a, **kw):\n        return _nyx_host.{n}(*a, **kw)\n\n", n = c.name));
        }
    }
    s.push_str("\nnyx = _Nyx()\n");
    s
}

/// The names of the prelude as Python declarations, appended to the stub of `nyx`. They are declared as `Any` variables, not as functions: a
/// script is free to use `region`, `material`, `step` or `norm` as a variable name, and the static check would refuse to assign a text to a
/// declared function.
pub fn stub() -> String {
    let mut s = String::from("\n# formula functions of a traced script function (nyx_script::trace); each is a function of one to three values\n");
    for (n, arity) in all_names() {
        let doc = match n.as_str() {
            "where" => "where(cond, a, b): a if cond else b, at every point (conditions: lt, gt, le, ge, eq, ne, land, lor, lnot)".to_string(),
            "trace" => "trace(f): the formula text of a function of the point (called once)".to_string(),
            "material" => "material(\"rho\"): a material value, by its name".to_string(),
            "region" => "region(\"Si\"): 1 in the cells of the region, else 0".to_string(),
            _ => format!("{n}: a function of the formula language with {arity} value(s), numbers or values of the point"),
        };
        s.push_str(&format!("{n}: Any  # {doc}\n"));
    }
    s
}

// ---- the guard ----

/// What a branch on the point would look like to a reader of the error.
fn problem(name: &str, src: &str, offset: usize, what: &str) -> String {
    let line = src[..offset.min(src.len())].bytes().filter(|&b| b == b'\n').count() + 1;
    format!("{name}:{line}: {what} on a value of the point cannot be traced: {NO_BRANCH}")
}

/// One function that may be traced: a `def` or a `lambda`, with the names of its parameters.
struct Func<'a> {
    params: Vec<&'a str>,
    body: Body<'a>,
}

enum Body<'a> {
    Stmts(&'a [Stmt]),
    Expr(&'a Expr),
}

fn param_names(p: &Parameters) -> Vec<&str> {
    let mut v: Vec<&str> = p.posonlyargs.iter().chain(p.args.iter()).map(|a| a.parameter.name.as_str()).collect();
    if let Some(a) = &p.vararg {
        v.push(a.name.as_str());
    }
    v.extend(p.kwonlyargs.iter().map(|a| a.parameter.name.as_str()));
    if let Some(a) = &p.kwarg {
        v.push(a.name.as_str());
    }
    v
}

fn defs_of(body: &[Stmt]) -> HashMap<String, Vec<&StmtFunctionDef>> {
    struct Collect<'a>(HashMap<String, Vec<&'a StmtFunctionDef>>);
    impl<'a> Visitor<'a> for Collect<'a> {
        fn visit_stmt(&mut self, stmt: &'a Stmt) {
            if let Stmt::FunctionDef(f) = stmt {
                self.0.entry(f.name.to_string()).or_default().push(f);
            }
            visitor::walk_stmt(self, stmt);
        }
    }
    let mut c = Collect(HashMap::new());
    for s in body {
        c.visit_stmt(s);
    }
    c.0
}

/// The names an expression reads that are in `set`.
fn mentions(e: &Expr, set: &HashSet<String>) -> bool {
    struct M<'s> {
        set: &'s HashSet<String>,
        hit: bool,
    }
    impl<'a> Visitor<'a> for M<'_> {
        fn visit_expr(&mut self, e: &'a Expr) {
            if let Expr::Name(n) = e {
                if self.set.contains(n.id.as_str()) {
                    self.hit = true;
                }
            }
            if !self.hit {
                visitor::walk_expr(self, e);
            }
        }
    }
    let mut m = M { set, hit: false };
    m.visit_expr(e);
    m.hit
}

/// The names a target binds: `a`, `(a, b)`, `a.b` and `a[i]` (which change what `a` holds).
fn target_names(t: &Expr, out: &mut Vec<String>) {
    match t {
        Expr::Name(n) => out.push(n.id.to_string()),
        Expr::Tuple(x) => x.elts.iter().for_each(|e| target_names(e, out)),
        Expr::List(x) => x.elts.iter().for_each(|e| target_names(e, out)),
        Expr::Starred(x) => target_names(&x.value, out),
        Expr::Attribute(x) => target_names(&x.value, out),
        Expr::Subscript(x) => target_names(&x.value, out),
        _ => {}
    }
}

/// Analysis of one function (a body with a set of names that depend on the point): grows the set to a fixed point, collects the problems found
/// and the calls of other functions of the script with the arguments that carry the point.
struct Scan<'a, 's> {
    src: &'s str,
    script: &'s str,
    tainted: HashSet<String>,
    problems: Vec<String>,
    /// (function name, positions of arguments that depend on the point, keywords that do, a starred argument depends)
    calls: Vec<(String, Vec<usize>, Vec<String>, bool)>,
    changed: bool,
    defs: &'s HashMap<String, Vec<&'a StmtFunctionDef>>,
}

impl Scan<'_, '_> {
    fn taint(&mut self, names: Vec<String>) {
        for n in names {
            if self.tainted.insert(n) {
                self.changed = true;
            }
        }
    }

    fn flow(&mut self, targets: &[&Expr], value: &Expr) {
        if mentions(value, &self.tainted) {
            let mut v = vec![];
            targets.iter().for_each(|t| target_names(t, &mut v));
            self.taint(v);
        }
    }

    fn report(&mut self, offset: usize, what: &str) {
        let m = problem(self.script, self.src, offset, what);
        if !self.problems.contains(&m) {
            self.problems.push(m);
        }
    }

    fn test(&mut self, test: &Expr, offset: usize, what: &str) {
        if mentions(test, &self.tainted) {
            self.report(offset, what);
        }
    }

    fn comp(&mut self, c: &Comprehension) {
        self.flow(&[&c.target], &c.iter);
        for cond in &c.ifs {
            self.test(cond, cond.range().start().to_usize(), "a condition in a comprehension");
        }
    }
}

impl<'a> Visitor<'a> for Scan<'a, '_> {
    fn visit_stmt(&mut self, stmt: &'a Stmt) {
        match stmt {
            // a function nested in this one is read with this one's names (a closure over the point); one defined elsewhere is read when it is
            // called with the point (see `guard`)
            Stmt::Assign(a) => self.flow(&a.targets.iter().collect::<Vec<_>>(), &a.value),
            Stmt::AugAssign(a) => self.flow(&[&a.target], &a.value),
            Stmt::AnnAssign(a) => {
                if let Some(v) = &a.value {
                    self.flow(&[&a.target], v);
                }
            }
            Stmt::For(f) => self.flow(&[&f.target], &f.iter),
            Stmt::With(w) => {
                for item in &w.items {
                    if let Some(v) = &item.optional_vars {
                        self.flow(&[v], &item.context_expr);
                    }
                }
            }
            Stmt::If(i) => {
                self.test(&i.test, i.range.start().to_usize(), "an `if`");
                for c in &i.elif_else_clauses {
                    if let Some(t) = &c.test {
                        self.test(t, c.range.start().to_usize(), "an `elif`");
                    }
                }
            }
            Stmt::While(w) => self.test(&w.test, w.range.start().to_usize(), "a `while`"),
            Stmt::Assert(a) => self.test(&a.test, a.range.start().to_usize(), "an `assert`"),
            Stmt::Match(m) => self.test(&m.subject, m.range.start().to_usize(), "a `match`"),
            _ => {}
        }
        visitor::walk_stmt(self, stmt);
    }

    fn visit_expr(&mut self, e: &'a Expr) {
        let at = e.range().start().to_usize();
        match e {
            Expr::If(x) => self.test(&x.test, at, "a conditional expression"),
            Expr::BoolOp(_) => self.test(e, at, "`and` / `or`"),
            Expr::UnaryOp(u) if u.op == UnaryOp::Not => self.test(&u.operand, at, "`not`"),
            Expr::Compare(_) => self.test(e, at, "a comparison"),
            Expr::Named(n) => self.flow(&[&n.target], &n.value),
            Expr::ListComp(c) => c.generators.iter().for_each(|g| self.comp(g)),
            Expr::SetComp(c) => c.generators.iter().for_each(|g| self.comp(g)),
            Expr::DictComp(c) => c.generators.iter().for_each(|g| self.comp(g)),
            Expr::Generator(c) => c.generators.iter().for_each(|g| self.comp(g)),
            Expr::Call(c) => {
                if let Expr::Attribute(a) = &*c.func {
                    if matches!(&*a.value, Expr::Name(n) if n.id.as_str() == "math") && c.arguments.args.iter().any(|x| mentions(x, &self.tainted)) {
                        self.report(at, &format!("`math.{}`", a.attr));
                    }
                }
                if let Expr::Name(n) = &*c.func {
                    let f = n.id.as_str();
                    if (f == "min" || f == "max") && c.arguments.args.iter().any(|x| mentions(x, &self.tainted)) {
                        self.report(at, &format!("`{f}()`"));
                    }
                    if self.defs.contains_key(f) {
                        let pos: Vec<usize> = c.arguments.args.iter().enumerate().filter(|(_, x)| mentions(x, &self.tainted)).map(|(i, _)| i).collect();
                        let kws: Vec<String> =
                            c.arguments.keywords.iter().filter(|k| mentions(&k.value, &self.tainted)).filter_map(|k| k.arg.as_ref().map(|a| a.to_string())).collect();
                        let star = c.arguments.args.iter().any(|x| matches!(x, Expr::Starred(_)) && mentions(x, &self.tainted))
                            || c.arguments.keywords.iter().any(|k| k.arg.is_none() && mentions(&k.value, &self.tainted));
                        if !pos.is_empty() || !kws.is_empty() || star {
                            self.calls.push((f.to_string(), pos, kws, star));
                        }
                    }
                }
            }
            _ => {}
        }
        visitor::walk_expr(self, e);
    }
}

/// Finds the functions a script hands over as formulas: the arguments of `nyx.<anything>(...)` and of `trace(...)` that are the name of a `def`,
/// a `lambda`, or a list of them.
fn roots<'a>(body: &'a [Stmt], defs: &HashMap<String, Vec<&'a StmtFunctionDef>>) -> Vec<Func<'a>> {
    struct Find<'a, 'd> {
        defs: &'d HashMap<String, Vec<&'a StmtFunctionDef>>,
        out: Vec<Func<'a>>,
    }
    impl<'a> Find<'a, '_> {
        fn arg(&mut self, e: &'a Expr) {
            match e {
                Expr::Name(n) => {
                    for f in self.defs.get(n.id.as_str()).into_iter().flatten() {
                        self.out.push(Func { params: param_names(&f.parameters), body: Body::Stmts(&f.body) });
                    }
                }
                Expr::Lambda(l) => {
                    let params = l.parameters.as_deref().map(param_names).unwrap_or_default();
                    self.out.push(Func { params, body: Body::Expr(&l.body) });
                }
                Expr::List(l) => l.elts.iter().for_each(|x| self.arg(x)),
                Expr::Tuple(l) => l.elts.iter().for_each(|x| self.arg(x)),
                _ => {}
            }
        }
    }
    impl<'a> Visitor<'a> for Find<'a, '_> {
        fn visit_expr(&mut self, e: &'a Expr) {
            if let Expr::Call(c) = e {
                let to_engine = matches!(&*c.func, Expr::Attribute(a) if matches!(&*a.value, Expr::Name(n) if n.id.as_str() == "nyx"));
                let traces = matches!(&*c.func, Expr::Name(n) if n.id.as_str() == "trace");
                if to_engine || traces {
                    c.arguments.args.iter().for_each(|x| self.arg(x));
                    c.arguments.keywords.iter().for_each(|k| self.arg(&k.value));
                }
            }
            visitor::walk_expr(self, e);
        }
    }
    let mut f = Find { defs, out: vec![] };
    for s in body {
        f.visit_stmt(s);
    }
    f.out
}

/// Refuses a script in which a function given as a formula branches or compares on a value of the point. `Ok` for a script the parser does not
/// understand (Monty reports a syntax error better) and for one with no function given as a formula.
pub fn guard(name: &str, src: &str) -> Result<(), String> {
    let Ok(parsed) = ruff_python_parser::parse_module(src) else {
        return Ok(());
    };
    let body: &[Stmt] = &parsed.syntax().body;
    let defs = defs_of(body);
    let mut problems: Vec<String> = vec![];
    // functions of the script that receive the point, with the parameters that carry it
    let mut carried: HashMap<String, HashSet<String>> = HashMap::new();
    let mut work: Vec<(Func, HashSet<String>)> = roots(body, &defs).into_iter().map(|f| {
        let t = f.params.iter().map(|p| p.to_string()).collect();
        (f, t)
    }).collect();
    let mut rounds = 0;
    while let Some((f, start)) = work.pop() {
        rounds += 1;
        if rounds > 500 {
            break;
        }
        let mut scan = Scan { src, script: name, tainted: start, problems: vec![], calls: vec![], changed: true, defs: &defs };
        while scan.changed {
            scan.changed = false;
            scan.calls.clear();
            match f.body {
                Body::Stmts(b) => b.iter().for_each(|s| scan.visit_stmt(s)),
                Body::Expr(e) => scan.visit_expr(e),
            }
        }
        for p in scan.problems {
            if !problems.contains(&p) {
                problems.push(p);
            }
        }
        // a helper called with a value of the point: its parameters carry it, and it is analysed with them
        for (callee, pos, kws, star) in scan.calls {
            for d in defs.get(&callee).into_iter().flatten() {
                let names = param_names(&d.parameters);
                let mut grow: Vec<String> = vec![];
                if star {
                    grow.extend(names.iter().map(|n| n.to_string()));
                }
                grow.extend(pos.iter().filter_map(|&i| names.get(i)).map(|n| n.to_string()));
                grow.extend(kws.iter().filter(|k| names.contains(&k.as_str())).cloned());
                let have = carried.entry(callee.clone()).or_default();
                let before = have.len();
                have.extend(grow);
                if have.len() > before {
                    work.push((Func { params: names, body: Body::Stmts(&d.body) }, have.clone()));
                }
            }
        }
    }
    if problems.is_empty() {
        Ok(())
    } else {
        problems.truncate(5);
        Err(problems.join("\n"))
    }
}
