//! What a front-end pays per frame, idle and while a solve runs (productivity plan P0). A project is opened and solved once (so that there is a
//! result to look at), then a loop at the cadence of a window (one tick per `--tick-ms`, default 16) makes the engine calls the window makes
//! on every repaint and times each one: `result()` (which hashes the whole project through `stamp()`), `stamp()`, `mesh_key()` and `sizing()`
//! alone, `model()`, `jobs().list()`, the progress status of the job, `poll()` and a geometry evaluation (`geometry_names`). The same loop
//! runs first with no job (idle) and then while the project is solved again on a job. The loop also records how late each tick woke up: that is
//! the delay a drawing thread would see from the solve's threads.
//!
//!   cargo run --release -p nyx-engine --example latency -- FILE.nyx.json [--physics ems] [--mesh '{"h":0.0025}'] [--pad N] [--idle-ticks 300] [--tick-ms 16]
//!
//! This measures the engine calls and the scheduling delay of a thread beside the solver's pool. It is not a measurement of window frames:
//! egui drawing, the viewport and the GPU are not in it (`NYX_FRAME_LOG=frames.csv` on the window records those). Run from the repository root.

use nyx_engine::Engine;
use serde_json::json;
use std::time::{Duration, Instant};

const NAMES: [&str; 10] = [
    "result()",
    "stamp()",
    "mesh_key()",
    "sizing()",
    "model()",
    "jobs().list()",
    "progress.status()",
    "poll()",
    "geometry_names",
    "tick wake delay",
];

fn pct(v: &[f64], p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    s[(((s.len() - 1) as f64) * p).round() as usize]
}

/// One tick of what a window does on a repaint; pushes the microseconds of each call into `series`.
fn tick(e: &mut Engine, physics: &str, job: Option<u64>, series: &mut [Vec<f64>]) {
    let mut t = |i: usize, f: &mut dyn FnMut()| {
        let t0 = Instant::now();
        f();
        series[i].push(t0.elapsed().as_secs_f64() * 1e6);
    };
    t(0, &mut || {
        std::hint::black_box(e.result(physics).map(|(_, stale)| stale));
    });
    t(1, &mut || {
        std::hint::black_box(e.stamp());
    });
    t(2, &mut || {
        std::hint::black_box(e.mesh_key());
    });
    t(3, &mut || {
        std::hint::black_box(e.sizing().is_some());
    });
    t(4, &mut || {
        std::hint::black_box(e.model().is_some());
    });
    t(5, &mut || {
        std::hint::black_box(e.jobs().list().len());
    });
    if let Some(id) = job {
        t(6, &mut || {
            std::hint::black_box(e.jobs().progress(id).map(|p| p.status().stage));
        });
    }
    t(7, &mut || {
        std::hint::black_box(e.poll().len());
    });
    t(8, &mut || {
        std::hint::black_box(e.call("geometry_names", &json!({})).is_ok());
    });
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let Some(file) = a.first() else {
        eprintln!(
            "usage: latency FILE.nyx.json [--physics CODE] [--mesh JSON] [--pad N] [--idle-ticks N] [--tick-ms N]"
        );
        std::process::exit(2);
    };
    let opt = |name: &str| a.iter().position(|x| x == name).and_then(|i| a.get(i + 1)).cloned();
    let idle_ticks: usize = opt("--idle-ticks").and_then(|s| s.parse().ok()).unwrap_or(300);
    let tick_ms: u64 = opt("--tick-ms").and_then(|s| s.parse().ok()).unwrap_or(16);
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": file })).unwrap_or_else(|er| panic!("project_open: {er}"));
    if let Some(m) = opt("--mesh") {
        e.call("set_mesh", &serde_json::from_str(&m).expect("--mesh is JSON"))
            .unwrap_or_else(|er| panic!("set_mesh: {er}"));
    }
    // a stand-in for a feature-heavy model: `--pad N` adds N rectangles (names `pad_i`, far from the model, not part of any region) before the first solve;
    // the cost of `stamp()` grows with the size of the project JSON, not with the mesh
    for i in 0..opt("--pad").and_then(|s| s.parse::<usize>().ok()).unwrap_or(0) {
        let x = 100.0 + i as f64;
        e.call("add_rect", &json!({ "name": format!("pad_{i}"), "min": [x, 100.0], "max": [x + 0.5, 100.5] }))
            .unwrap_or_else(|er| panic!("add_rect: {er}"));
    }
    let physics = opt("--physics")
        .or_else(|| e.project().active.clone())
        .unwrap_or_else(|| e.project().physics[0].physics.clone());
    let t0 = Instant::now();
    e.solve_now(None).unwrap_or_else(|er| panic!("first solve: {er}"));
    let first = t0.elapsed().as_secs_f64();
    let cells = e.model().map(|m| m.mesh().n_cells()).unwrap_or(0);
    println!(
        "{file}: physics {physics}, {cells} cells, first solve (mesh included) {first:.2} s, project JSON {} bytes",
        e.project().to_json().len()
    );

    let cadence = Duration::from_millis(tick_ms);
    let run = |e: &mut Engine, job: Option<u64>, ticks: Option<usize>| -> Vec<Vec<f64>> {
        let mut series: Vec<Vec<f64>> = vec![vec![]; NAMES.len()];
        let start = Instant::now();
        let mut k = 0u32;
        loop {
            let want = start + cadence * k;
            let now = Instant::now();
            if now < want {
                std::thread::sleep(want - now);
            }
            series[9].push(Instant::now().saturating_duration_since(want).as_secs_f64() * 1e6);
            tick(e, &physics, job, &mut series);
            k += 1;
            match (ticks, job) {
                (Some(n), _) if k as usize >= n => break,
                (None, Some(id)) if e.jobs().info(id).is_none_or(|i| i.state != nyx_engine::JobState::Running) => break,
                _ => {}
            }
        }
        // the finishing tick left the outcome in the job table: apply it (the cost of `poll()` that takes a result is in the series above)
        series
    };
    let idle = run(&mut e, None, Some(idle_ticks));
    let id = e.start_solve(None).expect("start the second solve");
    let t1 = Instant::now();
    let load = run(&mut e, Some(id), None);
    let wall = t1.elapsed().as_secs_f64();
    let _ = e.poll();
    println!(
        "second solve on a job: {wall:.2} s wall, {} ticks of {tick_ms} ms while it ran ({} idle ticks before)\n",
        load[9].len(),
        idle[9].len()
    );
    println!("{:<20} {:>26}   {:>26}", "", "idle (microseconds)", "while the solve runs (microseconds)");
    println!("{:<20} {:>8} {:>8} {:>8}   {:>8} {:>8} {:>8}", "call", "p50", "p95", "max", "p50", "p95", "max");
    for (i, n) in NAMES.iter().enumerate() {
        let f = |v: &[f64]| (pct(v, 0.5), pct(v, 0.95), pct(v, 1.0));
        let ((a, b, c), (d, ee, ff)) = (f(&idle[i]), f(&load[i]));
        println!("{n:<20} {a:>8.1} {b:>8.1} {c:>8.1}   {d:>8.1} {ee:>8.1} {ff:>8.1}");
    }
}
