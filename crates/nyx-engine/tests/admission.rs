// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Estimates and admission (productivity plan Q2, ADR 86): what solves cost is learned from the runs of workers; a solve whose estimate does not
//! fit the memory free now waits and says why, one the machine could never hold is refused unless forced, a cancel ends the wait, and a
//! solve with no estimate (fewer than two sizes seen) starts as before. The machine's memory is played by a probe, so nothing here depends on
//! what the test machine has free.

use nyx_engine::estimate::Row;
use nyx_engine::{Engine, JobState};
use serde_json::json;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn open(name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    e
}

/// An engine that has solved `name` once (so it holds the mesh), with two history rows around the size of that mesh whose memory is `gb_small`
/// and `gb_large` times 1e9 bytes: the estimate for this project is then in between, whatever the machine.
fn engine_with_history(name: &str, gb_small: f64, gb_large: f64) -> (Engine, String, u8, u64) {
    let mut e = open(name);
    e.call("solve_run", &json!({ "isolate": false })).unwrap();
    let cells = e.model().unwrap().mesh().n_cells() as u64;
    let dim = e.model().unwrap().mesh().dim() as u8;
    let chain = e.project().default_targets().join(",");
    for (c, gb) in [(cells / 2, gb_small), (cells * 2, gb_large)] {
        e.history().add(Row {
            chain: chain.clone(),
            dim,
            cells: c,
            dofs: c * 3,
            peak_working_set: (gb * 1e9) as u64,
            peak_commit: (gb * 1.2e9) as u64,
            seconds: 10.0,
        });
    }
    (e, chain, dim, cells)
}

fn machine(total_gb: f64, free_gb: f64) -> (Arc<AtomicU64>, nyx_engine::estimate::MemoryProbe) {
    let free = Arc::new(AtomicU64::new((free_gb * 1e9) as u64));
    let f = free.clone();
    let total = (total_gb * 1e9) as u64;
    (free, Arc::new(move || Some((total, f.load(Ordering::SeqCst)))))
}

#[test]
fn a_solve_with_no_estimate_starts_at_once() {
    // no history, or one size only: nothing is invented
    let mut e = open("02_cantilever");
    let (_, probe) = machine(16.0, 0.001);
    e.set_memory_probe(probe);
    let r = e.call("solve_run", &json!({ "isolate": false })).unwrap();
    assert_eq!(r["state"], "done");
    let cells = e.model().unwrap().mesh().n_cells() as u64;
    let dim = e.model().unwrap().mesh().dim() as u8;
    let chain = e.project().default_targets().join(",");
    e.history().add(Row {
        chain,
        dim,
        cells,
        dofs: 1,
        peak_working_set: 900_000_000_000,
        peak_commit: 1,
        seconds: 1.0,
    });
    // still one size
    assert_eq!(e.call("solve_run", &json!({ "isolate": false })).unwrap()["state"], "done");
}

#[test]
fn a_solve_waits_for_memory_says_why_and_starts_when_it_is_free() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    let (free, probe) = machine(2000.0, 4.0);
    e.set_memory_probe(probe);
    let id = e.call("solve_run", &json!({ "isolate": false, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let start = Instant::now();
    let mut stage = String::new();
    while start.elapsed() < Duration::from_secs(20) {
        stage = e.jobs().info(id).unwrap().stage;
        if stage.starts_with("waiting for memory") {
            break;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    assert!(stage.starts_with("waiting for memory: needs about") && stage.contains("4.0 GB free"), "{stage}");
    // the entry shows the estimate and its basis
    let entry = e.call("queue_get", &json!({ "id": id })).unwrap();
    assert_eq!(entry["state"], "running");
    assert_eq!(entry["estimate"]["runs"], 2, "{entry}");
    assert!(entry["estimate"]["bytes"].as_u64().unwrap() > 50_000_000_000, "{entry}");
    // memory frees up: it starts and finishes
    free.store(1_900_000_000_000, Ordering::SeqCst);
    let r = e.wait(id).unwrap();
    assert_eq!(r["state"], "done");
    assert_eq!(e.jobs().info(id).unwrap().state, JobState::Done);
}

#[test]
fn a_cancel_ends_the_wait() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    let (_, probe) = machine(2000.0, 4.0);
    e.set_memory_probe(probe);
    let id = e.call("solve_run", &json!({ "isolate": false, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let start = Instant::now();
    while !e.jobs().info(id).unwrap().stage.starts_with("waiting for memory")
        && start.elapsed() < Duration::from_secs(20)
    {
        std::thread::sleep(Duration::from_millis(20));
    }
    let asked = Instant::now();
    assert!(e.jobs().cancel(id));
    assert!(matches!(e.wait(id).unwrap_err(), nyx_physics::Error::Cancelled));
    assert!(asked.elapsed() < Duration::from_secs(2));
}

#[test]
fn what_the_machine_cannot_hold_is_refused_unless_forced() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    // the estimate is about 100 to 400 GB; this machine has 64 GB
    let (_, probe) = machine(64.0, 64.0);
    e.set_memory_probe(probe);
    let err = e.call("solve_run", &json!({ "isolate": false })).unwrap_err().to_string();
    assert!(err.contains("more than this machine can hold") && err.contains("force"), "{err}");
    assert!(err.contains("from 2 earlier runs"), "{err}");
    // forced, it runs (here the estimate is wrong by design; the test only shows the way through)
    let r = e.call("solve_run", &json!({ "isolate": false, "force": true })).unwrap();
    assert_eq!(r["state"], "done");
}

#[test]
fn a_worker_run_teaches_the_next_estimate_and_the_history_survives_in_a_file() {
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"));
    let file = std::env::temp_dir().join(format!("nyx-admission-{}", std::process::id())).join("solve_history.jsonl");
    let mut e = open("14_box_cavity");
    e.set_worker(exe.clone(), false);
    e.set_history_file(&file);
    // two sizes
    for h in [0.006, 0.0045] {
        e.call("set_mesh", &json!({ "h": h })).unwrap();
        e.call("solve_run", &json!({ "isolate": true })).unwrap();
    }
    assert_eq!(e.history().rows().len(), 2);
    // a third size gets an estimate from them, in the entry of its job
    e.call("set_mesh", &json!({ "h": 0.005 })).unwrap();
    let id = e.call("solve_run", &json!({ "isolate": true, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let r = e.wait(id).unwrap();
    assert_eq!(r["state"], "done");
    let entry = e.queue_entry(id).unwrap();
    let est = &entry["estimate"];
    assert_eq!(est["runs"], 2, "{entry}");
    assert_eq!(est["extrapolated"], false, "between the two sizes");
    let actual = r["worker"]["peak_working_set"].as_f64().unwrap();
    let predicted = est["bytes"].as_f64().unwrap();
    assert!(predicted / actual > 0.5 && predicted / actual < 2.0, "predicted {predicted}, measured {actual}");
    // a new engine reads what was learned
    e.history().flush();
    let again = Engine::new();
    again.history().set_file(&file);
    assert_eq!(again.history().rows().len(), 3);
    let _ = std::fs::remove_dir_all(file.parent().unwrap());
}

/// A line of the history is written by a thread of the history (the file open and write took 10 ms on the thread that draws, plan P3): the rows in
/// memory are there at once, the file once `flush` returns or the history is dropped, and in the order the rows were added.
#[test]
fn the_history_file_is_written_beside_the_caller_in_order() {
    use nyx_engine::estimate::{History, Row};
    let dir = std::env::temp_dir().join(format!("nyx-history-writer-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let file = dir.join("solve_history.jsonl");
    let row = |cells: u64| Row {
        chain: "sms".into(),
        dim: 2,
        cells,
        dofs: cells * 2,
        peak_working_set: 1 << 20,
        peak_commit: 1 << 20,
        seconds: 0.1,
    };
    let h = History::default();
    h.set_file(&file);
    for n in 1..=50 {
        h.add(row(n * 100));
    }
    assert_eq!(h.rows().len(), 50, "the rows an estimate reads are in memory at once");
    h.flush();
    let text = std::fs::read_to_string(&file).unwrap();
    let cells: Vec<u64> = text.lines().map(|l| serde_json::from_str::<Row>(l).unwrap().cells).collect();
    assert_eq!(cells, (1..=50).map(|n| n * 100).collect::<Vec<_>>(), "every line, in order");
    // dropping flushes: a row added just before the end of the program is not lost
    h.add(row(9_999));
    drop(h);
    let again = History::default();
    again.set_file(&file);
    assert_eq!(again.rows().len(), 51);
    assert_eq!(again.rows().last().unwrap().cells, 9_999);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn nyx_memory_wait_timeout_env_child_probe() {
    let Ok(mode) = std::env::var("NYX_TEST_MEM_TIMEOUT_MODE") else {
        return;
    };
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    let (_, probe) = machine(2000.0, 4.0);
    e.set_memory_probe(probe);
    match mode.as_str() {
        "env_beats_caller_default" => {
            // Caller default is 10.0 s, NYX_MEMORY_WAIT_TIMEOUT_S is 0.35 s -> aborts in < 1.2 s naming 0.35
            e.set_memory_wait_timeout(Some(10.0));
            let t0 = Instant::now();
            let err = e.call("solve_run", &json!({ "isolate": false })).unwrap_err().to_string();
            let dt = t0.elapsed().as_secs_f64();
            assert!((0.30..=0.85).contains(&dt), "env override timing {dt} s: {err}");
            assert!(
                err.contains("0.35") && err.contains("need") && err.contains("free") && err.contains("total"),
                "{err}"
            );
        }
        "explicit_beats_env" => {
            // NYX_MEMORY_WAIT_TIMEOUT_S is 10.0 s, explicit memory_wait_timeout_s is 0.2 s -> aborts in < 0.8 s naming 0.2
            e.set_memory_wait_timeout(Some(20.0));
            let t0 = Instant::now();
            let err = e
                .call("solve_run", &json!({ "isolate": false, "memory_wait_timeout_s": 0.2 }))
                .unwrap_err()
                .to_string();
            let dt = t0.elapsed().as_secs_f64();
            assert!((0.15..=0.70).contains(&dt), "explicit override timing {dt} s: {err}");
            assert!(
                err.contains("0.2") && err.contains("need") && err.contains("free") && err.contains("total"),
                "{err}"
            );
        }
        other => panic!("unknown probe mode {other}"),
    }
}

#[test]
fn admit_timeout_aborts_promptly_names_need_free_total_timeout_and_obeys_priority() {
    // 1. Direct `estimate::admit` with `timeout_s = Some(0.5)` completes within `<= 0.5 + 0.5 s` (< 1.5 s)
    let (mut e, chain, dim, cells) = engine_with_history("02_cantilever", 100.0, 400.0);
    let est = e.history().estimate(&chain, dim, cells).expect("estimate from 2 history rows");
    let (_, probe) = machine(2000.0, 4.0);
    let progress = femcore::progress::Progress::new();

    let t0 = Instant::now();
    let err = nyx_engine::estimate::admit(&est, false, &probe, &progress, Some(0.5)).unwrap_err();
    let elapsed = t0.elapsed().as_secs_f64();
    assert!(
        (0.45..=1.0).contains(&elapsed),
        "admit(Some(0.5)) must abort within [0.45, 1.0] s (< 1.5 s), took {elapsed:.3} s"
    );
    assert!(matches!(err, nyx_physics::Error::Invalid(_)), "expected Error::Invalid, got {err:?}");
    let msg = err.to_string();
    assert!(
        msg.contains("need")
            && msg.contains("free")
            && msg.contains("total")
            && msg.contains("timeout_s")
            && msg.contains("0.5"),
        "timeout message must name need, free, total, and timeout_s: {msg}"
    );

    // 2. Zero timeout aborts immediately (< 0.1 s), and invalid negative/NaN timeouts are rejected upfront
    let t_zero = Instant::now();
    assert!(nyx_engine::estimate::admit(&est, false, &probe, &progress, Some(0.0)).is_err());
    assert!(t_zero.elapsed() < Duration::from_millis(100));
    assert!(nyx_engine::estimate::admit(&est, false, &probe, &progress, Some(-1.0)).is_err());
    assert!(nyx_engine::estimate::admit(&est, false, &probe, &progress, Some(f64::NAN)).is_err());
    assert!(nyx_engine::estimate::admit(&est, false, &probe, &progress, Some(f64::INFINITY)).is_err());

    // 3. Engine `solve_run` and `set_memory_wait_timeout` plumbing
    e.set_memory_probe(probe);
    assert!(e.call("solve_run", &json!({ "isolate": false, "memory_wait_timeout_s": -0.5 })).is_err());
    let t_eng = Instant::now();
    let eng_err =
        e.call("solve_run", &json!({ "isolate": false, "memory_wait_timeout_s": 0.5 })).unwrap_err().to_string();
    assert!(t_eng.elapsed().as_secs_f64() < 1.5, "solve_run timeout took {:?}", t_eng.elapsed());
    assert!(eng_err.contains("need") && eng_err.contains("free") && eng_err.contains("total"), "{eng_err}");

    // 4. Priority order: explicit `timeout_s` > `NYX_MEMORY_WAIT_TIMEOUT_S` > caller default (`set_memory_wait_timeout`)
    let cur_exe = std::env::current_exe().unwrap();
    for (mode, env_val) in [("env_beats_caller_default", "0.35"), ("explicit_beats_env", "10.0")] {
        let out = std::process::Command::new(&cur_exe)
            .args(["--exact", "nyx_memory_wait_timeout_env_child_probe", "--nocapture"])
            .env("NYX_TEST_MEM_TIMEOUT_MODE", mode)
            .env("NYX_MEMORY_WAIT_TIMEOUT_S", env_val)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "priority probe {mode} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
}
