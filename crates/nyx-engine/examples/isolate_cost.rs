// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What running a solve in a worker process costs (productivity plan P4): the same project is solved in this process and in a worker, each in a
//! fresh engine (so both mesh), `--reps` times, and the wall time of each is printed with the stages that belong to the worker (the results
//! read back) and the memory the worker reports. The difference of the best times is the price of isolation: starting the program, meshing
//! twice (once here for the mesh the window draws, once in the worker), writing and reading the results file.
//!
//!   cargo run --release -p nyx-engine --example isolate_cost -- FILE.nyx.json [--physics ems] [--mesh '{"h":0.0025}'] [--reps 3] [--worker PROGRAM]
//!
//! `--worker` names the program to run as the worker (`nyx.exe` or `nyx-mcp.exe`, which answer `--job`); the default is `nyx-worker` next to the examples.
//!
//! Run from the repository root, with `nyx-worker` built (`cargo build --release -p nyx-engine --bin nyx-worker`).

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::time::Instant;

fn main() {
    nyx_engine::process::opt_out_of_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut file = None;
    let mut worker_arg = None::<String>;
    let (mut physics, mut mesh, mut reps) = (None::<String>, None::<Value>, 3usize);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--physics" => {
                i += 1;
                physics = args.get(i).cloned();
            }
            "--mesh" => {
                i += 1;
                mesh = args.get(i).and_then(|t| serde_json::from_str(t).ok());
            }
            "--worker" => {
                i += 1;
                worker_arg = args.get(i).cloned();
            }
            "--reps" => {
                i += 1;
                reps = args.get(i).and_then(|t| t.parse().ok()).unwrap_or(3);
            }
            a => file = Some(a.to_string()),
        }
        i += 1;
    }
    let file = file.expect("usage: isolate_cost FILE.nyx.json [--physics ems] [--mesh JSON] [--reps N]");
    let exe = std::env::current_exe().expect("exe");
    let worker = worker_arg.map(std::path::PathBuf::from).unwrap_or_else(|| {
        exe.parent()
            .and_then(|d| d.parent())
            .map(|d| d.join(if cfg!(windows) { "nyx-worker.exe" } else { "nyx-worker" }))
            .expect("path")
    });
    assert!(worker.is_file(), "build the worker first: {}", worker.display());
    let mut rows: Vec<(bool, f64, String)> = vec![];
    for rep in 0..reps {
        for isolate in [false, true] {
            let mut e = Engine::new();
            e.call("project_open", &json!({ "path": file })).expect("open");
            if let Some(m) = &mesh {
                e.call("set_mesh", m).expect("mesh settings");
            }
            e.set_worker(worker.clone(), false);
            let mut a = json!({ "isolate": isolate });
            if let Some(p) = &physics {
                a["physics"] = json!(p);
            }
            let t = Instant::now();
            let r = e.call("solve_run", &a).unwrap_or_else(|e| panic!("solve: {e}"));
            let wall = t.elapsed().as_secs_f64();
            let solve_s = r["seconds"].as_f64().unwrap_or(f64::NAN);
            let mut note = format!("rep {rep}: job {solve_s:.2} s");
            if isolate {
                let w = &r["worker"];
                let gb = |k: &str| w[k].as_u64().map(|b| format!("{:.2}", b as f64 / 1e9)).unwrap_or("?".into());
                note +=
                    &format!(", worker peak {} GB in RAM / {} GB committed", gb("peak_working_set"), gb("peak_commit"));
            }
            let dofs: Vec<String> = r["results"]
                .as_array()
                .map(|a| a.iter().map(|x| format!("{} {} unknowns", x["physics"], x["dofs"])).collect())
                .unwrap_or_default();
            println!(
                "{} {:>7.3} s  {note}  [{}]",
                if isolate { "worker    " } else { "in process" },
                wall,
                dofs.join(", ")
            );
            rows.push((isolate, wall, note));
        }
    }
    let best = |iso: bool| rows.iter().filter(|r| r.0 == iso).map(|r| r.1).fold(f64::INFINITY, f64::min);
    let (a, b) = (best(false), best(true));
    println!(
        "best of {reps}: in process {a:.3} s, in a worker {b:.3} s, difference {:+.3} s ({:+.1} %)",
        b - a,
        (b / a - 1.0) * 100.0
    );
}
