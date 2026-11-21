// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Resilience tests for `nyx-script` (FR-7): per-turn Monty suspension budget across >25k cumulative calls,
//! clean snapshot recovery after per-turn call limit or MemoryError aborts, `Kernel::reset`, UTF-8 BOM stripping,
//! and dict recursion in formula tracing and static branch guards.

use std::path::PathBuf;

use nyx_engine::commandlist;
use nyx_engine::Engine;
use nyx_script::kernel::{Kernel, Step};
use nyx_script::{host, trace, Limits, Options};
use serde_json::{json, Value};

fn root(name: &str) -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!(
        "nyx-resilience-test-{name}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn opts(root: &std::path::Path) -> Options {
    Options { root: root.to_path_buf(), worker: Some(PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))) }
}

fn kernel(e: &Engine, limits: Limits) -> (Kernel, PathBuf) {
    let r = root("k");
    let k = Kernel::new(e, limits, opts(&r)).unwrap();
    (k, r)
}

fn go(k: &mut Kernel, e: &mut Engine, src: &str) -> Result<Value, String> {
    match k.run(e, "t.nys", src, None, None, true)? {
        Step::Done(o) => Ok(o.value),
        Step::Running(p) => panic!("still running, waiting for {}", p.waiting),
    }
}

#[test]
fn multi_turn_over_25k_cumulative_calls_with_small_per_turn_limit_succeeds_and_preserves_state() {
    let mut e = Engine::new();
    let limits = Limits { calls: 1_000, seconds: 60.0, ..Limits::default() };
    let (mut k, r) = kernel(&e, limits);

    go(&mut k, &mut e, "total = 0\ndef bump(n):\n    global total\n    for _ in range(n):\n        nyx.project_info()\n    total += n\n    return total\n").unwrap();

    // 30 turns * 900 calls/turn = 27,000 cumulative calls (>25k cumulative, each turn <= 1,000 per-turn limit)
    for turn in 1..=30 {
        let v = go(&mut k, &mut e, "bump(900)").unwrap();
        assert_eq!(v, json!(turn * 900), "turn {turn}");
    }

    assert_eq!(go(&mut k, &mut e, "total").unwrap(), json!(27_000));
    let _ = std::fs::remove_dir_all(&r);
}

#[test]
fn single_turn_over_call_limit_fails_cleanly_and_next_turn_restores_clean_state() {
    let mut e = Engine::new();
    let limits = Limits { calls: 80, seconds: 30.0, ..Limits::default() };
    let (mut k, r) = kernel(&e, limits);

    go(&mut k, &mut e, "saved = [42]\ncounter = 10\n").unwrap();

    // Turn exceeds per-turn limit of 80 calls after mutating REPL variables
    let err = go(
        &mut k,
        &mut e,
        "saved.append(999)\ndirty = 'leaked'\nfor _ in range(120):\n    nyx.project_info()\n",
    )
    .unwrap_err();
    assert!(err.contains("call limit 80 exceeded for this turn"), "unexpected error: {err}");

    // Next turn succeeds cleanly and sees the snapshot from the last completed turn
    assert_eq!(go(&mut k, &mut e, "saved[0] + counter + len(saved)").unwrap(), json!(53));
    assert!(go(&mut k, &mut e, "dirty").is_err(), "dirty variable from aborted turn must not survive");

    // Kernel reset clears REPL state without touching the engine
    e.call("project_new", &json!({ "dim": 2, "title": "kept_in_engine" })).unwrap();
    k.reset().unwrap();
    assert!(go(&mut k, &mut e, "saved").is_err(), "reset must clear REPL globals");
    assert_eq!(go(&mut k, &mut e, "nyx.project_info()['title']").unwrap(), json!("kept_in_engine"));

    let _ = std::fs::remove_dir_all(&r);
}

#[test]
fn utf8_bom_scripts_and_commandlists_pass_check_and_run() {
    let mut e = Engine::new();
    let (mut k, r) = kernel(&e, Limits::default());

    let bom_script = "\u{feff}nyx.project_new(dim=2, title='bom_ok')\nx = 19\nx + 23\n";
    k.check(&mut e, "bom.nys", bom_script).unwrap();
    assert_eq!(go(&mut k, &mut e, bom_script).unwrap(), json!(42));

    // Also via host::check and host::run with traced formula
    let o = opts(&r);
    let bom_host = "\u{feff}nyx.project_new(dim=2, title='bom_host')\ntrace(lambda p: abs2(p['Ex']) + p.x)\n";
    host::check(&mut e, "bom_host.nys", bom_host, &Limits::default(), &o).unwrap();
    let out = host::run(&mut e, "bom_host.nys", bom_host, &Limits::default(), &o).unwrap();
    assert_eq!(out.value, json!("(abs2(Ex) + x)"));

    // And commandlist::parse strips leading UTF-8 BOM
    let bom_jsonl = "\u{feff}{\"cmd\": \"project_new\", \"args\": {\"dim\": 2, \"title\": \"bom_list\"}}\n";
    let (lines, problems) = commandlist::parse(bom_jsonl);
    assert!(problems.is_empty(), "{problems:?}");
    assert_eq!(lines.len(), 1);
    assert_eq!(lines[0].cmd, "project_new");

    let _ = std::fs::remove_dir_all(&r);
}

#[test]
fn dict_traced_lambdas_convert_and_static_branch_guard_checks_dict_values() {
    // 1. Static guard rejects a branch on a point value inside a dict formula argument
    let bad_dict = "nyx.integrate(expr={'ex': lambda p: 1 if p.x > 0 else 2, 'ez': lambda p: abs2(p['Ez'])})";
    let err = trace::guard("bad_dict.nys", bad_dict).unwrap_err();
    assert!(err.contains("bad_dict.nys:1:") && err.contains("where(cond, a, b)"), "{err}");

    let bad_def_in_dict = "def f(p):\n    if p['Ex']:\n        return 1\n    return 0\nnyx.integrate(expr={'a': f})";
    let err = trace::guard("bad_def.nys", bad_def_in_dict).unwrap_err();
    assert!(err.contains("bad_def.nys:2:") && err.contains("where(cond, a, b)"), "{err}");

    // 2. Runtime _nyx_conv recursively converts functions and Sym values inside dicts
    let mut e = Engine::new();
    let (mut k, r) = kernel(&e, Limits::default());
    let src = r#"
d = trace({
    "ex": lambda p: abs2(p["Ex"]),
    "ez": lambda p: where(gt(p.z, 0), p["eps_r"] * abs2(p["Ez"]), 0),
    "raw": "norm2(E)",
})
d
"#;
    let v = go(&mut k, &mut e, src).unwrap();
    assert_eq!(v["ex"], json!("abs2(Ex)"));
    assert_eq!(v["ez"], json!("where((z > 0), (eps_r * abs2(Ez)), 0)"));
    assert_eq!(v["raw"], json!("norm2(E)"));

    let _ = std::fs::remove_dir_all(&r);
}
