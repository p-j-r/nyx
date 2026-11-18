//! S4 spike runner: `spike sample star|rhai [file]` runs the strip sweep sample; `spike bench` measures the cost of one call and of the loop around it;
//! `spike limits` runs the runaway scripts.

use std::time::Instant;

use nyx_engine::Engine;
use nyx_script::{rhai_host, star, Limits, Outcome};
use serde_json::{json, Value};

fn show(label: &str, r: Result<Outcome, String>, t: Instant) {
    match r {
        Ok(o) => {
            for l in &o.printed {
                println!("  | {l}");
            }
            println!("{label}: ok in {:.2} s, {} engine calls, value {}", t.elapsed().as_secs_f64(), o.calls, o.value);
        }
        Err(e) => println!("{label}: FAILED after {:.2} s\n{e}", t.elapsed().as_secs_f64()),
    }
}

fn fresh() -> Engine {
    let mut e = Engine::new();
    e.call("project_new", &json!({ "dim": 2, "title": "spike" })).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 5e-7 })).unwrap();
    e
}

/// seconds per iteration of `f` run `n` times
fn per_call(n: u64, mut f: impl FnMut()) -> f64 {
    let t = Instant::now();
    for _ in 0..n {
        f();
    }
    t.elapsed().as_secs_f64() / n as f64
}

fn bench() {
    let n: u64 = 20_000;
    let limits = Limits { steps: u64::MAX / 2, ..Limits::default() };
    println!("{:<34} {:>12} {:>12} {:>12}", "per call, microseconds", "rust", "starlark", "rhai");
    for (name, args_rust, star_call, rhai_call) in [
        ("set_param (edit)", json!({ "name": "w", "value": 5e-7 }), "nyx.set_param(name = \"w\", value = 5e-7)", "nyx::set_param(#{ name: \"w\", value: 5e-7 })"),
        ("project_info (read, reply of a map)", json!({}), "nyx.project_info()", "nyx::project_info(#{})"),
        ("list_materials (read, reply of a list)", json!({}), "nyx.list_materials()", "nyx::list_materials(#{})"),
    ] {
        let cmd = name.split_whitespace().next().unwrap().to_string();
        let mut e = fresh();
        let rust = per_call(n, || {
            e.call(&cmd, &args_rust).unwrap();
        });
        let mut e = fresh();
        let src = format!("def go():\n    for i in range({n}):\n        {star_call}\ngo()\n");
        let t = Instant::now();
        star::run(&mut e, "bench.star", &src, &limits).unwrap();
        let st = t.elapsed().as_secs_f64() / n as f64;
        let mut e = fresh();
        let src = format!("for i in 0..{n} {{ {rhai_call}; }}");
        let t = Instant::now();
        rhai_host::run(&mut e, "bench.rhai", &src, &limits).unwrap();
        let rh = t.elapsed().as_secs_f64() / n as f64;
        println!("{name:<34} {:>12.2} {:>12.2} {:>12.2}", rust * 1e6, st * 1e6, rh * 1e6);
    }
    // the loop alone: how much does the interpreter cost per iteration with no engine call
    let mut e = fresh();
    let m: u64 = 2_000_000;
    let t = Instant::now();
    star::run(&mut e, "loop.star", &format!("def go():\n    s = 0.0\n    for i in range({m}):\n        s += i * 0.5\n    return s\nresult = go()\n"), &limits).unwrap();
    let st = t.elapsed().as_secs_f64() / m as f64;
    let t = Instant::now();
    rhai_host::run(&mut e, "loop.rhai", &format!("let s = 0.0; for i in 0..{m} {{ s += i.to_float() * 0.5; }} s"), &limits).unwrap();
    let rh = t.elapsed().as_secs_f64() / m as f64;
    println!("{:<34} {:>12} {:>12.3} {:>12.3}   (loop of arithmetic, no engine call)", "per iteration, microseconds", "-", st * 1e6, rh * 1e6);
    // start-up: a script that does nothing
    let mut e = fresh();
    let st = per_call(50, || {
        star::run(&mut e, "empty.star", "result = 1", &limits).unwrap();
    });
    let rh = per_call(50, || {
        rhai_host::run(&mut e, "empty.rhai", "1", &limits).unwrap();
    });
    println!("{:<34} {:>12} {:>12.1} {:>12.1}   (an empty script, bindings of {} commands)", "start-up, microseconds", "-", st * 1e6, rh * 1e6, e.catalog().len());
}

fn limits_demo() {
    let tight = Limits { steps: 2_000_000, heap_bytes: 64 << 20, depth: 50, ..Limits::default() };
    let mut e = fresh();
    for (label, lang, src) in [
        ("starlark: huge loop", "star", "def go():\n    s = 0\n    for i in range(100000000):\n        s += i\n    return s\nresult = go()\n"),
        ("starlark: unbounded list", "star", "def go():\n    x = [0] * 1000\n    for i in range(100000):\n        x = x + x\n    return len(x)\nresult = go()\n"),
        ("starlark: one huge list", "star", "result = len([0] * 10000000000000)
"),
        ("starlark: unbounded string", "star", "def go():\n    s = \"x\"\n    for i in range(100):\n        s = s + s\n    return len(s)\nresult = go()\n"),
        ("rhai: while true", "rhai", "let i = 0; while true { i += 1; } i"),
        ("rhai: unbounded array", "rhai", "let x = [0]; loop { x += x; } x.len()"),
        ("rhai: unbounded string", "rhai", "let s = \"x\"; loop { s += s; } s.len()"),
        ("rhai: deep recursion", "rhai", "fn f(n) { f(n + 1) } f(0)"),
        ("starlark: recursion", "star", "def f(n):\n    return f(n + 1)\nresult = f(0)\n"),
    ] {
        let t = Instant::now();
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| if lang == "star" { star::run(&mut e, "x.star", src, &tight) } else { rhai_host::run(&mut e, "x.rhai", src, &tight) })).unwrap_or_else(|_| Err("PANIC inside the interpreter".into()));
        let msg = match r {
            Ok(o) => format!("NOT STOPPED, value {}", o.value),
            Err(m) => m.lines().rfind(|l| !l.trim().is_empty() && !l.contains("Traceback") && !l.starts_with(" ")).unwrap_or("").chars().take(110).collect(),
        };
        println!("{label:<30} stopped in {:>6.2} s: {msg}", t.elapsed().as_secs_f64());
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("bench") => bench(),
        Some("limits") => limits_demo(),
        Some("direct") => println!("{}", direct()),
        Some("sample") => {
            let lang = args.get(1).map(String::as_str).unwrap_or("star");
            let file = args.get(2).cloned().unwrap_or_else(|| format!("docs/scripting-samples/strip_sweep.{}", if lang == "star" { "star" } else { "rhai" }));
            let src = std::fs::read_to_string(&file).unwrap_or_else(|e| panic!("{file}: {e}"));
            let mut e = Engine::new();
            let t = Instant::now();
            let limits = Limits { steps: 2_000_000_000, ..Limits::default() };
            let r = if lang == "star" { star::run(&mut e, &file, &src, &limits) } else { rhai_host::run(&mut e, &file, &src, &limits) };
            show(lang, r, t);
        }
        _ => eprintln!("usage: spike bench | limits | sample star|rhai [file]"),
    }
    let _ = Value::Null;
}

/// The same workflow as the samples, written as plain engine calls (what an agent does over MCP, call by call). The numbers must equal the scripts'.
fn direct() -> Value {
    let mut e = Engine::new();

    let call = |e: &mut Engine, name: &str, args: Value| -> Value {
        CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        e.call(name, &args).unwrap_or_else(|err| panic!("{name}: {err}"))
    };
    let (c0, lambda) = (299792458.0_f64, 1.55e-6_f64);
    call(&mut e, "project_new", json!({ "dim": 2, "title": "SOI strip width sweep" }));
    call(&mut e, "set_param", json!({ "name": "w", "value": 500e-9 }));
    call(&mut e, "set_param", json!({ "name": "t", "value": 220e-9 }));
    call(&mut e, "add_rect", json!({ "name": "box", "min": [-2.0e-6, -2.0e-6], "max": [2.0e-6, 0.0] }));
    call(&mut e, "add_rect", json!({ "name": "clad", "min": [-2.0e-6, 0.0], "max": [2.0e-6, 2.0e-6] }));
    call(&mut e, "add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "BOX", "input": "box" }, { "name": "Clad", "input": "clad" }, { "name": "Si", "input": "core" }] }));
    call(&mut e, "add_physics", json!({ "physics": "ems", "study": "modes", "study_params": { "freq": c0 / lambda, "n_target": 2.9, "n_modes": 4 } }));
    for (r, m) in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")] {
        call(&mut e, "set_material", json!({ "physics": "ems", "region": r, "material": m }));
    }
    let solve = |e: &mut Engine, do_solve: bool| -> (Value, Value) {
        if do_solve {
            call(e, "solve_run", json!({ "physics": "ems", "wait": true }));
        }
        let pol = call(e, "analyse", json!({ "analysis": "polarization", "physics": "ems" }));
        let (mut te, mut tm) = (Value::Null, Value::Null);
        for row in pol["rows"].as_array().unwrap() {
            let cls = row["report"]["class"].as_str().unwrap();
            if cls.starts_with("quasi-TE") && te.is_null() {
                te = row["value"][0].clone();
            }
            if cls.starts_with("quasi-TM") && tm.is_null() {
                tm = row["value"][0].clone();
            }
        }
        (te, tm)
    };
    let mut ladder = vec![];
    let mut h = 200e-9_f64;
    for _ in 0..8 {
        if h < 25e-9 {
            break;
        }
        ladder.push(h);
        h /= 2.0;
    }
    let mut prev: Option<f64> = None;
    let mut conv = None;
    for (i, h) in ladder.iter().enumerate() {
        call(&mut e, "set_mesh", json!({ "h": h, "adapt_scale": 0.4 }));
        let cur = solve(&mut e, true).0.as_f64().unwrap();
        if let Some(p) = prev {
            if (cur - p).abs() < 2e-3 {
                conv = Some(ladder[i - 1]);
                break;
            }
        }
        prev = Some(cur);
    }
    let h = conv.expect("not converged");
    call(&mut e, "set_mesh", json!({ "h": h, "adapt_scale": 0.4 }));
    let mut rows = vec![];
    for w_nm in (400..701).step_by(100) {
        call(&mut e, "set_param", json!({ "name": "w", "value": w_nm as f64 * 1e-9 }));
        call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true }));
        let (te, tm) = solve(&mut e, false);
        rows.push(json!({ "w_nm": w_nm, "te0": te, "tm0": tm }));
    }
    println!("direct: {} engine calls, converged h = {} nm", CALLS.load(std::sync::atomic::Ordering::Relaxed), h * 1e9);
    Value::Array(rows)
}
static CALLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
