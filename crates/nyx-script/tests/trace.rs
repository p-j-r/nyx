//! Tracing of script functions (productivity plan S6, N3, ADR 90): a function given where a formula is expected is called once, becomes formula
//! text, and gives the number of the same formula written by hand, to the bit. A branch on a value of the point is refused with the line and the
//! `where(...)` spelling.

use std::path::PathBuf;

use nyx_engine::Engine;
use nyx_script::{host, trace, Limits, Options};


fn root() -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!("nyx-trace-test-{}-{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(src: &str) -> Result<nyx_script::Outcome, String> {
    let r = root();
    let mut e = Engine::new();
    let o = Options { root: r.clone(), worker: Some(PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))) };
    let out = host::run(&mut e, "test.nys", src, &Limits { seconds: 60.0, ..Limits::default() }, &o);
    let _ = std::fs::remove_dir_all(&r);
    out
}

/// The strip of the S4 parity test: an SOI wire with TE0 and TM0, solved once.
const MODEL: &str = r#"
nyx.project_new(dim=2, title="trace")
nyx.set_param(name="w", value=5e-7)
nyx.set_param(name="t", value=2.2e-7)
nyx.add_rect(name="box", min=[-2.0e-6, -2.0e-6], max=[2.0e-6, 0.0])
nyx.add_rect(name="clad", min=[-2.0e-6, 0.0], max=[2.0e-6, 2.0e-6])
nyx.add_rect(name="core", min=["-w/2", 0.0], max=["w/2", "t"])
nyx.add_overlay(regions=[{"name": "BOX", "input": "box"}, {"name": "Clad", "input": "clad"}, {"name": "Si", "input": "core"}])
nyx.add_physics(physics="ems", study="modes", study_params={"freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 4})
for reg, mat in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")]:
    nyx.set_material(physics="ems", region=reg, material=mat)
nyx.set_mesh(h=4e-7, adapt_scale=0.4)
nyx.solve_run(physics="ems", wait=True)
"#;

// ------------------------------------------------------------------------------------------------------------------------ the gate

/// A traced function equals the formula written by hand to the bit, and is called once, not at every element.
#[test]
fn a_traced_function_equals_the_written_formula_to_the_bit_and_runs_once() {
    let body = r#"
calls = 0

def integrand(p):
    global calls
    calls += 1
    e = p["E"]
    return p["eps_r"] * (abs(e.x) ** 2 + abs(e.z) ** 2)

def helper(e):
    return abs(e.x) ** 2

def with_helper(p):
    total = 0
    for k in range(3):
        total = total + (k + 1) * helper(p["E"])
    return total

cases = [
    ["energy", integrand, "eps_r * (abs(E.x)^2 + abs(E.z)^2)"],
    ["lambda", lambda p: abs(p["Ex"]) ** 2 + abs(p["Ez"]) ** 2, "abs(Ex)^2 + abs(Ez)^2"],
    ["reflected", lambda p: abs(2 * p["Ex"] - 1 / p["eps_r"]), "abs(2 * Ex - 1 / eps_r)"],
    ["where", lambda p: where(gt(p.z, 0), p["eps_r"] * norm2(p["E"]), 0), "where(z > 0, eps_r * norm2(E), 0)"],
    ["region", lambda p: region("Si") * p["eps_r"] * abs2(p["Ey"]), 'region("Si") * eps_r * abs2(Ey)'],
    ["numbers", lambda p: 0.1 * p["eps_r"] - 1e-6 * p.x + 3, "0.1 * eps_r - 1e-6 * x + 3"],
    ["sqrt", lambda p: sqrt(abs2(p["Ex"]) + 1e-12) / (1 + p["eps_r"]), "sqrt(abs2(Ex) + 1e-12) / (1 + eps_r)"],
    ["negative", lambda p: -abs(p["Ex"]) + 2 ** 0.5 * p["eps_r"], "-abs(Ex) + 1.4142135623730951 * eps_r"],
    ["helper", with_helper, "0 + 1 * abs(E.x)^2 + 2 * abs(E.x)^2 + 3 * abs(E.x)^2"],
    ["land", lambda p: where(land(gt(p.z, 0), lt(p.x, 1e-7)), p["eps_r"], lnot(gt(p.x, 0))), "where(z > 0 and x < 1e-7, eps_r, not x > 0)"],
]
pairs = []
for label, f, text in cases:
    a = nyx.integrate(expr=f, physics="ems")["value"]
    b = nyx.integrate(expr=text, physics="ems")["value"]
    pairs.append([label, a, b])

calls_after_energy = calls

w1 = nyx.confinement(fields=["Ex", "Ez"], weight=lambda p: p["eps_r"], region="Si", physics="ems")["fraction"]
w2 = nyx.confinement(fields=["Ex", "Ez"], weight="eps_r", region="Si", physics="ems")["fraction"]
pts = [[0.0, 1.1e-7], [1e-7, 1e-7]]
pr1 = nyx.probe(expr=lambda p: abs2(p["Ex"]) * p["eps_r"], points=pts, physics="ems")["values"]
pr2 = nyx.probe(expr="abs2(Ex) * eps_r", points=pts, physics="ems")["values"]
lc1 = nyx.line_cut(expr=lambda p: norm(p["E"]), n=5, physics="ems", **{"from": [-1e-6, 1.1e-7], "to": [1e-6, 1.1e-7]})["values"]
lc2 = nyx.line_cut(expr="norm(E)", n=5, physics="ems", **{"from": [-1e-6, 1.1e-7], "to": [1e-6, 1.1e-7]})["values"]
# a symbolic value, not a function, is a formula too; a call through nyx.call and through a positional dict as well
sv = nyx.call("integrate", expr=lambda p: abs2(p["Ex"]), physics="ems")["value"]
sd = nyx.call("integrate", {"expr": lambda p: abs2(p["Ex"]), "physics": "ems"})["value"]
sx = nyx.integrate(expr="abs2(Ex)", physics="ems")["value"]

def message(thunk):
    try:
        thunk()
    except Exception as ex:
        return str(ex)
    return "no error"

errors = [
    message(lambda: nyx.integrate(expr=lambda e, eps: e, physics="ems")),
    message(lambda: nyx.integrate(expr=lambda p: p["no_such_name"] + 1, physics="ems")),
    message(lambda: trace(lambda p: float(p["Ex"]))),
    message(lambda: trace(lambda p: p["Ex"].__eq__(1))),
    message(lambda: trace(lambda p: int(p["Ex"]))),
    message(lambda: trace(lambda p: "text")),
    message(lambda: trace(lambda p: 1e308 * 10 * p["Ex"])),
]
{"pairs": pairs, "calls": calls_after_energy, "weight": [w1, w2], "probe": [pr1, pr2], "line": [lc1, lc2], "forms": [sv, sd, sx], "errors": errors,
 "text": trace(lambda p: 2 * p["Ex"] - 1)}
"#;
    let out = run(&format!("{MODEL}{body}")).unwrap();
    let v = &out.value;
    for pair in v["pairs"].as_array().unwrap() {
        // the text of a JSON number is its shortest round-trip form: equal text is equal bits (a complex value is a pair)
        let label = pair[0].as_str().unwrap();
        let (a, b) = (pair[1].to_string(), pair[2].to_string());
        assert!(a != "0.0" && a != "[0.0,0.0]", "{label}: the integral is {a}");
        assert_eq!(a, b, "{label}: traced against written");
    }
    assert_eq!(v["pairs"].as_array().unwrap().len(), 10);
    assert_eq!(v["calls"], 1, "the function ran once for the command, not at every element");
    assert_eq!(v["weight"][0], v["weight"][1], "confinement weight");
    assert_eq!(v["probe"][0], v["probe"][1], "probe");
    assert_eq!(v["line"][0], v["line"][1], "line_cut");
    assert_eq!(v["forms"][0], v["forms"][2], "nyx.call");
    assert_eq!(v["forms"][1], v["forms"][2], "a positional dict");
    assert_eq!(v["text"], "((2 * Ex) - 1)");
    let errors: Vec<&str> = v["errors"].as_array().unwrap().iter().map(|e| e.as_str().unwrap()).collect();
    assert!(errors[0].contains("takes one argument, the point"), "{}", errors[0]);
    assert!(errors[1].contains("no_such_name") && errors[1].contains("integrate"), "{}", errors[1]);
    assert!(errors[2].contains("not a number here"), "{}", errors[2]);
    assert!(errors[3].contains("where("), "{}", errors[3]);
    assert!(errors[4].contains("not a number here"), "{}", errors[4]);
    assert!(errors[5].contains("str"), "{}", errors[5]);
    assert!(errors[6].contains("finite"), "{}", errors[6]);
}

/// The text a traced function gives: every operation parenthesised, operand order kept, numbers by repr.
#[test]
fn the_text_keeps_operand_order_and_parenthesises_every_operation() {
    let src = r#"
[
    trace(lambda p: 1 - p["Ex"] ** 2 / -p.x),
    trace(lambda p: where(gt(p.z, 0), sqrt(p["Ex"]), 0)),
    trace(lambda p: 2 * p["Ex"] + p["Ex"] * 2),
    trace(lambda p: 1e-7 * p["Ex"] - 2.5e20 / p.x + 5),
    trace(lambda p: p["Re Ex"] + p["E"].z + abs(p["u"]) ** 3 + 2 ** p.x),
    trace(lambda p: material("rho") * region("Si") + minimum(p.x, 1.5) - 0.1),
    trace(lambda p: 7),
]
"#;
    let out = run(src).unwrap();
    let t: Vec<&str> = out.value.as_array().unwrap().iter().map(|x| x.as_str().unwrap()).collect();
    assert_eq!(t[0], "(1 - ((Ex ^ 2) / (-x)))");
    assert_eq!(t[1], "where((z > 0), sqrt(Ex), 0)");
    assert_eq!(t[2], "((2 * Ex) + (Ex * 2))");
    assert_eq!(t[3], "(((1e-07 * Ex) - (2.5e+20 / x)) + 5)");
    assert_eq!(t[4], "(((`Re Ex` + E.z) + (abs(u) ^ 3)) + (2 ^ x))");
    assert_eq!(t[5], "(((material(\"rho\") * region(\"Si\")) + min(x, 1.5)) - 0.1)");
    assert_eq!(t[6], "7");
}

// ------------------------------------------------------------------------------------------------------------------------ the rule

fn refused(src: &str, line: usize) {
    let e = trace::guard("test.nys", src).expect_err(&format!("this script should be refused:\n{src}"));
    assert!(e.contains(&format!("test.nys:{line}:")), "line {line} expected: {e}");
    assert!(e.contains("where(cond, a, b)"), "{e}");
}

#[test]
fn a_branch_on_a_value_of_the_point_is_refused_with_its_line() {
    refused("def f(p):\n    if p[\"Ex\"]:\n        return 1\n    return 2\nnyx.integrate(expr=f)", 2);
    refused("nyx.integrate(expr=lambda p: 1 if p.x > 0 else 2)", 1);
    refused("def h(q):\n    return q.x if q.x > 0 else 0\ndef f(p):\n    return h(p[\"E\"])\nnyx.integrate(expr=f)", 2);
    refused("def f(p):\n    e = p[\"E\"]\n    a = e.x\n    if a:\n        return 1\n    return 2\nnyx.integrate(expr=f)", 4);
    refused("def f(p):\n    d = {}\n    d[\"a\"] = p.x\n    while d[\"a\"]:\n        break\n    return 1\nnyx.integrate(expr=f)", 4);
    refused("nyx.integrate(expr=lambda p: p.x and 1)", 1);
    refused("nyx.integrate(expr=lambda p: not p.x)", 1);
    refused("nyx.integrate(expr=lambda p: p.x != 0)", 1);
    refused("nyx.probe(expr=[lambda p: p.x < 1], points=[[0, 0]])", 1);
    refused("import math\nnyx.integrate(expr=lambda p: math.sqrt(p.x))", 2);
    refused("nyx.integrate(expr=lambda p: min(p.x, 1))", 1);
    refused("nyx.integrate(expr=lambda p: sum([v for v in [p.x, p.z] if v]))", 1);
    refused("def f(p):\n    a, b = p.x, 1\n    assert b\n    for v in [a]:\n        if v:\n            pass\n    return a\nnyx.integrate(expr=f)", 5);
    refused("x = trace(lambda p: 1 if p.x else 2)", 1);
    refused("def f(p, k):\n    g = lambda v: v\n    def inner():\n        return 1 if p.x else 2\n    return inner()\nnyx.integrate(expr=f)", 4);
}

#[test]
fn a_branch_on_something_else_is_allowed_as_in_jax() {
    let ok = [
        // a plain constant, a loop bound, a helper called with numbers only
        "USE = True\ndef f(p):\n    t = 0\n    for k in range(3):\n        if k != 1 and USE:\n            t = t + k * p[\"Ex\"]\n    return t\nnyx.integrate(expr=f)",
        "def g(a):\n    return 1 if a > 0 else 2\ndef f(p):\n    return g(3) * p.x\nnyx.integrate(expr=f)",
        // branches in functions that are not given as formulas
        "def g(a):\n    if a:\n        return 1\n    return 2\nx = g(3)\nnyx.project_info()",
        "n = 3\nif n > 2:\n    nyx.project_info()",
        // a script the parser does not understand is Monty's to report
        "def f(p:\n",
    ];
    for s in ok {
        trace::guard("test.nys", s).unwrap_or_else(|e| panic!("this script is fine:\n{s}\n{e}"));
    }
}

#[test]
fn check_refuses_a_branch_and_runs_nothing_and_knows_the_formula_functions() {
    let r = root();
    let mut e = Engine::new();
    let o = Options { root: r.clone(), worker: Some(PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))) };
    let bad = "nyx.project_new(dim=2)\nnyx.integrate(expr=lambda p: p[\"Ex\"] if p.x else 0)\n";
    let err = host::check(&mut e, "test.nys", bad, &Limits::default(), &o).unwrap_err();
    assert!(err.contains("test.nys:2:") && err.contains("where("), "{err}");
    let good = "nyx.integrate(expr=lambda p: where(gt(p.z, 0), sqrt(abs2(p[\"Ex\"])) + exp(p.x) * material(\"rho\"), region(\"Si\")), physics=\"ems\")\nx = trace(lambda p: dot(p[\"E\"], conj(p[\"E\"])))\n";
    host::check(&mut e, "test.nys", good, &Limits::default(), &o).unwrap_or_else(|m| panic!("the stub should know the formula functions: {m}"));
    let _ = std::fs::remove_dir_all(&r);
}
