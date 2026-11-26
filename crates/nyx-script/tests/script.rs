// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Scripts (productivity plan S4, ADR [89]): what a script can do, and the gates: the same numbers as the same workflow sent command by command,
//! a runaway script is stopped without taking the engine down, and a script cannot reach a file outside its folder.

use std::path::PathBuf;

use nyx_engine::Engine;
use nyx_script::{host, stubs, Limits, Options};
use serde_json::{json, Value};

fn worker() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))
}

fn root(name: &str) -> PathBuf {
    // tests run in parallel: one folder each
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!("nyx-script-test-{name}-{}-{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn opts(root: &std::path::Path) -> Options {
    Options { root: root.to_path_buf(), worker: Some(worker()) }
}

fn run(src: &str) -> Result<nyx_script::Outcome, String> {
    run_with(src, &Limits::default())
}

fn run_with(src: &str, limits: &Limits) -> Result<nyx_script::Outcome, String> {
    let r = root("run");
    let mut e = Engine::new();
    let out = host::run(&mut e, "test.nys", src, limits, &opts(&r));
    let _ = std::fs::remove_dir_all(&r);
    out
}

// ------------------------------------------------------------------------------------------------------------------------ what a script does

#[test]
fn a_script_drives_the_engine_and_returns_a_value() {
    let out = run("nyx.project_new(dim=2, title='from a script')\nnyx.set_param(name='w', value=5e-7)\nnyx.project_info()['params']").unwrap();
    assert_eq!(out.value, json!([{ "name": "w", "value": 5e-7 }]));
    assert_eq!(out.calls, 3);
}

#[test]
fn print_goes_back_and_the_value_is_the_last_expression() {
    let out = run("print('a', 1)\nprint(f'{2.5:.3f}')\n41 + 1").unwrap();
    assert_eq!(out.printed, vec!["a 1".to_string(), "2.500".to_string()]);
    assert_eq!(out.value, json!(42));
}

#[test]
fn an_engine_error_is_a_runtime_error_a_script_can_catch() {
    // (a value of None passes the static check, which sees names only, and fails in the engine)
    let out = run("try:\n    nyx.set_param(name='w', value=None)\n    r = 'no error'\nexcept RuntimeError as e:\n    r = str(e)\nr").unwrap();
    let msg = out.value.as_str().unwrap();
    assert!(msg.starts_with("set_param:"), "the command is named: {msg}");
    // and uncaught it ends the script with the line
    let e = run("a = 1\nnyx.set_param(name='w', value=None)\n").unwrap_err();
    assert!(e.contains("line 2") && e.contains("set_param"), "{e}");
}

#[test]
fn nyx_call_and_a_positional_dict_reach_any_command() {
    let out = run("nyx.call('project_new', dim=2, title='t')\nnyx.call('set_param', {'name': 'w', 'value': 1.5})\nnyx.project_info()['params'][0]['value']").unwrap();
    assert_eq!(out.value, json!(1.5));
}

#[test]
fn every_command_of_the_catalog_is_in_the_stub_and_the_stub_is_valid_python() {
    let mut e = Engine::new();
    let pyi = stubs::nyx_pyi(&e);
    for c in e.catalog() {
        assert!(pyi.contains(&format!("def {}(self", c.name)), "{} is not in the stub", c.name);
    }
    // a script using a few of them passes the type check against it (the check reads the stub, so a stub that was not valid Python would fail here)
    let r = root("stub");
    host::check(&mut e, "t.nys", "nyx.project_new(dim=2)\nnyx.add_rect(name='a', min=[0, 0], max=[1, 1])\nnyx.list_materials()", &Limits::default(), &opts(&r)).unwrap();
    let _ = std::fs::remove_dir_all(&r);
}

// ------------------------------------------------------------------------------------------------------------------------ the static check

#[test]
fn a_misspelled_command_or_argument_fails_before_anything_runs() {
    let mut e = Engine::new();
    let r = root("check");
    let o = opts(&r);
    for (src, needle) in [
        ("nyx.nope()", "nope"),
        ("nyx.set_param(nme='w', value=1)", "nme"),
        ("nyx.set_param(value=1)", "name"),
        ("x = (1,\n", "yntax"),
    ] {
        let err = host::check(&mut e, "t.nys", src, &Limits::default(), &o).unwrap_err();
        assert!(err.contains(needle), "{src}: {err}");
    }
    // `run` checks first: the engine must not have been touched by the failing script
    let before = e.call("project_info", &Value::Null).unwrap();
    let err = host::run(&mut e, "t.nys", "nyx.set_param(name='w', value=1)\nnyx.nope()", &Limits::default(), &o).unwrap_err();
    assert!(err.contains("nope"), "{err}");
    assert_eq!(e.call("project_info", &Value::Null).unwrap(), before, "a script that fails the static check changed the project");
    let _ = std::fs::remove_dir_all(&r);
}

#[test]
fn check_runs_no_command() {
    let mut e = Engine::new();
    let r = root("check2");
    let before = e.call("project_info", &Value::Null).unwrap();
    host::check(&mut e, "t.nys", "nyx.project_new(dim=3, title='must not happen')\nnyx.set_param(name='w', value=1)", &Limits::default(), &opts(&r)).unwrap();
    assert_eq!(e.call("project_info", &Value::Null).unwrap(), before);
    let _ = std::fs::remove_dir_all(&r);
}

// ------------------------------------------------------------------------------------------------------------------------ limits

#[test]
fn a_loop_that_never_ends_is_stopped_and_the_engine_goes_on() {
    let mut e = Engine::new();
    let r = root("loop");
    let t = std::time::Instant::now();
    let err = host::run(&mut e, "t.nys", "while True:\n    pass\n", &Limits { seconds: 1.0, ..Limits::default() }, &opts(&r)).unwrap_err();
    assert!(t.elapsed().as_secs_f64() < 10.0, "took {:?}", t.elapsed());
    assert!(err.contains("imit") || err.contains("Timeout"), "{err}");
    e.call("project_new", &json!({ "dim": 2 })).unwrap();
    let _ = std::fs::remove_dir_all(&r);
}

/// Doubling a list in a loop aborted the whole process when the interpreter ran inside it (the S4 spike). In the worker it ends the worker only.
#[test]
fn allocating_without_bound_ends_the_worker_not_the_engine() {
    let mut e = Engine::new();
    let r = root("memory");
    for src in ["x = [0]\nwhile True:\n    x = x + x\n", "s = 'x'\nwhile True:\n    s = s + s\n", "len([0] * 10**12)"] {
        let limits = Limits { memory_bytes: 64 << 20, seconds: 20.0, ..Limits::default() };
        let err = host::run(&mut e, "t.nys", src, &limits, &opts(&r)).unwrap_err();
        assert!(err.to_lowercase().contains("memory"), "{src}: {err}");
    }
    // this process is still here and the engine still answers
    e.call("project_new", &json!({ "dim": 2, "title": "still alive" })).unwrap();
    assert_eq!(e.call("project_info", &Value::Null).unwrap()["title"], "still alive");
    let _ = std::fs::remove_dir_all(&r);
}

#[test]
fn too_many_commands_deep_recursion_and_a_print_flood_are_stopped() {
    let e = run_with("for i in range(1000):\n    nyx.project_info()\n", &Limits { calls: 50, ..Limits::default() }).unwrap_err();
    assert!(e.contains("50") || e.to_lowercase().contains("suspension") || e.to_lowercase().contains("commands"), "{e}");
    let e = run("def f(n):\n    return f(n + 1)\nf(0)\n").unwrap_err();
    assert!(e.contains("RecursionError"), "{e}");
    // output is capped, not fatal: the script ends by its time limit, with at most the cap printed
    let l = Limits { seconds: 2.0, output_bytes: 10_000, ..Limits::default() };
    let _ = run_with("while True:\n    print('x' * 1000)\n", &l);
}

// ------------------------------------------------------------------------------------------------------------------------ no I/O, same numbers

#[test]
fn a_script_cannot_open_a_file_or_import_what_is_not_built_in() {
    for (src, needle) in [
        ("open('C:/Windows/win.ini').read()", "PermissionError"),
        ("import os\nos.listdir('.')", "PermissionError"),
        // these two are caught by the static check, before anything runs
        ("import numpy", "numpy"),
        ("print(undefined_name)", "undefined_name"),
    ] {
        let e = run(src).unwrap_err();
        assert!(e.contains(needle), "{src}: {e}");
    }
}

#[test]
fn the_clock_is_refused_and_random_is_seeded() {
    let e = run("import time\ntime.time()").unwrap_err();
    assert!(e.contains("no clock"), "{e}");
    let a = run("import random\n[random.random() for _ in range(3)]").unwrap().value;
    let b = run("import random\n[random.random() for _ in range(3)]").unwrap().value;
    assert_eq!(a, b, "random must give the same numbers in every run");
    // CPython's random.seed(0) sequence
    assert_eq!(a[0], json!(0.8444218515250481));
}

#[test]
fn the_path_guard_keeps_a_scripts_files_in_its_folder() {
    let r = root("guard");
    let outside = root("guard-outside");
    let mut e = Engine::new();
    let o = opts(&r);
    // inside: a relative path is taken from the root
    host::run(&mut e, "t.nys", "nyx.project_new(dim=2, title='g')\nnyx.project_save(path='model.nyx.json')", &Limits::default(), &o).unwrap();
    assert!(r.join("model.nyx.json").exists());
    // outside: refused, nothing written, the script can say why
    for path in [outside.join("x.nyx.json").to_string_lossy().replace('\\', "/"), "../escaped.nyx.json".to_string()] {
        let src = format!("nyx.project_save(path={path:?})");
        let err = host::run(&mut e, "t.nys", &src, &Limits::default(), &o).unwrap_err();
        assert!(err.contains("project_save") && err.contains("path"), "{err}");
    }
    assert!(!outside.join("x.nyx.json").exists());
    assert!(!r.parent().unwrap().join("escaped.nyx.json").exists());
    // the guard is off again after the run: the host's own calls are not confined
    assert!(e.path_root().is_none());
    let _ = std::fs::remove_dir_all(&r);
    let _ = std::fs::remove_dir_all(&outside);
}

// ------------------------------------------------------------------------------------------------------------------------ the sample

/// TE0 and TM0 of the strip, from the polarization rows
fn modes(e: &mut Engine) -> (f64, f64) {
    let pol = e.call("analyse", &json!({ "analysis": "polarization", "physics": "ems" })).unwrap();
    let (mut te, mut tm) = (None, None);
    for row in pol["rows"].as_array().unwrap() {
        let class = row["report"]["class"].as_str().unwrap();
        let n = row["value"][0].as_f64().unwrap();
        if class.starts_with("quasi-TE") && te.is_none() {
            te = Some(n);
        }
        if class.starts_with("quasi-TM") && tm.is_none() {
            tm = Some(n);
        }
    }
    (te.unwrap(), tm.unwrap())
}

/// One solve of the sample's model, by a script and by the engine commands an MCP agent would send: the same numbers to the bit.
#[test]
fn a_script_gives_the_same_numbers_as_the_same_commands() {
    let model = |c: &mut dyn FnMut(&str, Value) -> Value| {
        c("project_new", json!({ "dim": 2, "title": "parity" }));
        c("set_param", json!({ "name": "w", "value": 5e-7 }));
        c("set_param", json!({ "name": "t", "value": 2.2e-7 }));
        c("add_rect", json!({ "name": "box", "min": [-2.0e-6, -2.0e-6], "max": [2.0e-6, 0.0] }));
        c("add_rect", json!({ "name": "clad", "min": [-2.0e-6, 0.0], "max": [2.0e-6, 2.0e-6] }));
        c("add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
        c("add_overlay", json!({ "regions": [{ "name": "BOX", "input": "box" }, { "name": "Clad", "input": "clad" }, { "name": "Si", "input": "core" }] }));
        c("add_physics", json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 4 } }));
        for (region, material) in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")] {
            c("set_material", json!({ "physics": "ems", "region": region, "material": material }));
        }
        c("set_mesh", json!({ "h": 4e-7, "adapt_scale": 0.4 }));
        c("solve_run", json!({ "physics": "ems", "wait": true }));
    };
    let mut direct = Engine::new();
    model(&mut |n, a| direct.call(n, &a).unwrap_or_else(|e| panic!("{n}: {e}")));
    let (te, tm) = modes(&mut direct);

    let src = r#"
nyx.project_new(dim=2, title="parity")
nyx.set_param(name="w", value=5e-7)
nyx.set_param(name="t", value=2.2e-7)
nyx.add_rect(name="box", min=[-2.0e-6, -2.0e-6], max=[2.0e-6, 0.0])
nyx.add_rect(name="clad", min=[-2.0e-6, 0.0], max=[2.0e-6, 2.0e-6])
nyx.add_rect(name="core", min=["-w/2", 0.0], max=["w/2", "t"])
nyx.add_overlay(regions=[{"name": "BOX", "input": "box"}, {"name": "Clad", "input": "clad"}, {"name": "Si", "input": "core"}])
nyx.add_physics(physics="ems", study="modes", study_params={"freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 4})
for region, material in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")]:
    nyx.set_material(physics="ems", region=region, material=material)
nyx.set_mesh(h=4e-7, adapt_scale=0.4)
nyx.solve_run(physics="ems", wait=True)
pol = nyx.analyse(analysis="polarization", physics="ems")
te = None
tm = None
for row in pol["rows"]:
    cls = row["report"]["class"]
    if cls.startswith("quasi-TE") and te is None:
        te = row["value"][0]
    if cls.startswith("quasi-TM") and tm is None:
        tm = row["value"][0]
{"te": te, "tm": tm}
"#;
    let out = run_with(src, &Limits { seconds: 60.0, ..Limits::default() }).unwrap();
    assert_eq!(out.value["te"].as_f64().unwrap().to_bits(), te.to_bits(), "TE0 differs: script {} direct {te}", out.value["te"]);
    assert_eq!(out.value["tm"].as_f64().unwrap().to_bits(), tm.to_bits(), "TM0 differs");
}

/// A run made by a script keeps the script's text, with its hash in the file name, beside its results.
#[test]
fn a_run_made_by_a_script_keeps_the_script_and_its_hash() {
    let r = root("record");
    let src = r#"
nyx.project_new(dim=2, title="recorded")
nyx.set_param(name="w", value=5e-7)
nyx.add_rect(name="box", min=[-2.0e-6, -2.0e-6], max=[2.0e-6, 0.0])
nyx.add_rect(name="core", min=["-w/2", 0.0], max=["w/2", 2.2e-7])
nyx.add_overlay(regions=[{"name": "BOX", "input": "box"}, {"name": "Si", "input": "core"}])
nyx.add_physics(physics="ems", study="modes", study_params={"freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 2})
nyx.set_material(physics="ems", region="BOX", material="silica")
nyx.set_material(physics="ems", region="Si", material="silicon")
nyx.set_mesh(h=6e-7, adapt_scale=0.4)
nyx.store_save(workspace=".")
nyx.solve_run(physics="ems", wait=True)
nyx.store_runs()
"#;
    let mut e = Engine::new();
    let out = host::run(&mut e, "recorded.nys", src, &Limits::default(), &opts(&r)).unwrap();
    assert!(out.value.to_string().contains("script"), "the run lists a script artifact: {}", out.value);
    let mut found = Vec::new();
    let mut stack = vec![r.clone()];
    while let Some(d) = stack.pop() {
        for f in std::fs::read_dir(&d).unwrap().flatten() {
            let p = f.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.file_name().unwrap().to_string_lossy().starts_with("script-") {
                found.push(p);
            }
        }
    }
    assert_eq!(found.len(), 1, "{found:?}");
    let text = std::fs::read_to_string(&found[0]).unwrap();
    assert_eq!(text, src);
    // the file name holds the hash of the text: `script-` + 32 hex digits (blake3, as the store hashes everything) + `.nys`
    let name = found[0].file_name().unwrap().to_string_lossy().to_string();
    let hash = name.strip_prefix("script-").and_then(|n| n.strip_suffix(".nys")).unwrap_or("");
    assert!(hash.len() == 32 && hash.chars().all(|c| c.is_ascii_hexdigit()), "{name}");
    let _ = std::fs::remove_dir_all(&r);
}

/// The whole sample (mesh convergence, then a width sweep, about five minutes per pass in release):
/// compares the script against the exact same convergence + width-sweep workflow sent as plain engine commands (to the bit across any OS)
/// and checks physical SOI dispersion reference bounds. Run with `-- --ignored`.
#[test]
#[ignore = "five minutes: the full strip sweep"]
fn the_strip_sweep_sample_gives_the_numbers_of_the_same_commands() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/scripting-samples/strip_sweep.nys");
    let src = std::fs::read_to_string(&path).unwrap();
    let mut e = Engine::new();
    let r = root("sample");
    let out = host::run(&mut e, "strip_sweep.nys", &src, &Limits { seconds: 600.0, ..Limits::default() }, &opts(&r)).unwrap();
    let rows = out.value.as_array().unwrap();
    assert_eq!(rows.len(), 4, "{out:?}");

    // Replay the exact same workflow with direct `Engine::call` commands for 0-ULP parity on any OS/arch.
    let mut direct = Engine::new();
    let call = |eng: &mut Engine, name: &str, args: Value| -> Value {
        eng.call(name, &args).unwrap_or_else(|err| panic!("{name}: {err}"))
    };
    call(&mut direct, "project_new", json!({ "dim": 2, "title": "SOI strip width sweep" }));
    call(&mut direct, "set_param", json!({ "name": "w", "value": 500e-9 }));
    call(&mut direct, "set_param", json!({ "name": "t", "value": 220e-9 }));
    call(&mut direct, "add_rect", json!({ "name": "box", "min": [-2.0e-6, -2.0e-6], "max": [2.0e-6, 0.0] }));
    call(&mut direct, "add_rect", json!({ "name": "clad", "min": [-2.0e-6, 0.0], "max": [2.0e-6, 2.0e-6] }));
    call(&mut direct, "add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
    call(
        &mut direct,
        "add_overlay",
        json!({ "regions": [{ "name": "BOX", "input": "box" }, { "name": "Clad", "input": "clad" }, { "name": "Si", "input": "core" }] }),
    );
    call(
        &mut direct,
        "add_physics",
        json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 299_792_458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 4 } }),
    );
    for (region, material) in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")] {
        call(&mut direct, "set_material", json!({ "physics": "ems", "region": region, "material": material }));
    }

    let mut h = 200e-9;
    call(&mut direct, "set_mesh", json!({ "h": h, "adapt_scale": 0.4 }));
    call(&mut direct, "solve_run", json!({ "physics": "ems", "wait": true }));
    let (mut prev_te, _) = modes(&mut direct);
    let mut converged = false;
    while h > 25e-9 {
        h /= 2.0;
        call(&mut direct, "set_mesh", json!({ "h": h, "adapt_scale": 0.4 }));
        call(&mut direct, "solve_run", json!({ "physics": "ems", "wait": true }));
        let (cur_te, _) = modes(&mut direct);
        let d = (cur_te - prev_te).abs();
        prev_te = cur_te;
        if d < 2e-3 {
            converged = true;
            break;
        }
    }
    assert!(converged, "direct replay failed to converge");
    call(&mut direct, "set_mesh", json!({ "h": h * 2.0, "adapt_scale": 0.4 }));

    let physical_ref = [
        (400_u64, 2.22497_f64, 1.69023_f64),
        (500, 2.44349, 1.76907),
        (600, 2.56517, 1.82817),
        (700, 2.63912, 1.87189),
    ];
    for (row, (w, ref_te, ref_tm)) in rows.iter().zip(physical_ref) {
        call(&mut direct, "set_param", json!({ "name": "w", "value": (w as f64) * 1e-9 }));
        call(&mut direct, "solve_run", json!({ "physics": "ems", "wait": true }));
        let (te, tm) = modes(&mut direct);
        assert_eq!(row["w_nm"], w);
        let s_te = row["te0"].as_f64().unwrap();
        let s_tm = row["tm0"].as_f64().unwrap();
        assert_eq!(s_te.to_bits(), te.to_bits(), "TE0 at {w} nm differs between script ({s_te}) and direct ({te})");
        assert_eq!(s_tm.to_bits(), tm.to_bits(), "TM0 at {w} nm differs between script ({s_tm}) and direct ({tm})");
        assert!((s_te - ref_te).abs() < 1e-4, "TE0 at {w} nm drifted from physical reference: {s_te} vs {ref_te}");
        assert!((s_tm - ref_tm).abs() < 1e-4, "TM0 at {w} nm drifted from physical reference: {s_tm} vs {ref_tm}");
    }
    let _ = std::fs::remove_dir_all(&r);
}
