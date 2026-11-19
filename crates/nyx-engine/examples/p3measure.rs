//! What the interactions of productivity plan P3 cost, per example project: the geometry preview (`geometry_names`, what the window evaluates on
//! every edit), the input check (`check_project`), the mesh build (`mesh_build`, to decide whether a small 2D mesh preview pays), and, with
//! `--polarization`, the polarization analysis of a solved mode study (what the window computes when a solve is taken).
//!
//!   cargo run --release -p nyx-engine --example p3measure -- DIR_OR_FILE... [--reps 20] [--polarization]
//!
//! Times are milliseconds, p50 and max over `--reps` repetitions. Run from the repository root; no file is written.

use nyx_engine::Engine;
use serde_json::json;
use std::time::Instant;

fn stats(mut v: Vec<f64>) -> (f64, f64) {
    v.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    (v[v.len() / 2], *v.last().unwrap_or(&f64::NAN))
}

fn time(reps: usize, mut f: impl FnMut()) -> (f64, f64) {
    stats(
        (0..reps)
            .map(|_| {
                let t = Instant::now();
                f();
                t.elapsed().as_secs_f64() * 1e3
            })
            .collect(),
    )
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let reps: usize =
        a.iter().position(|x| x == "--reps").and_then(|i| a.get(i + 1)).and_then(|s| s.parse().ok()).unwrap_or(20);
    let pol = a.iter().any(|x| x == "--polarization");
    let mut files: Vec<std::path::PathBuf> = vec![];
    let mut skip = false;
    for x in &a {
        if skip {
            skip = false;
            continue;
        }
        if x == "--reps" || x == "--mesh" {
            skip = true;
            continue;
        }
        if x.starts_with("--") {
            continue;
        }
        let p = std::path::PathBuf::from(x);
        if p.is_dir() {
            let mut v: Vec<_> = std::fs::read_dir(&p)
                .expect("read dir")
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.to_string_lossy().ends_with(".nyx.json"))
                .collect();
            v.sort();
            files.extend(v);
        } else {
            files.push(p);
        }
    }
    println!(
        "{:<28} {:>3} {:>8} | {:>16} | {:>16} | {:>16} | {:>16}",
        "example", "dim", "cells", "geometry cold ms", "check cold ms", "mesh cold ms", "polarization ms"
    );
    for f in files {
        let mut e = Engine::new();
        if let Err(er) = e.call("project_open", &json!({ "path": f })) {
            println!("{:<28} project_open: {er}", f.file_name().and_then(|s| s.to_str()).unwrap_or("?"));
            continue;
        }
        // `--mesh '{"h":1e-7}'` refines every example, to see how a cost grows with the mesh
        if let Some(m) = a.iter().position(|x| x == "--mesh").and_then(|i| a.get(i + 1)) {
            e.call("set_mesh", &serde_json::from_str(m).expect("--mesh is JSON"))
                .unwrap_or_else(|er| panic!("set_mesh: {er}"));
        }
        let name = f.file_name().and_then(|s| s.to_str()).unwrap_or("?").trim_end_matches(".nyx.json").to_string();
        // COLD times: every repetition is a fresh engine (the open and the mesh settings are not timed), so a cache of the engine cannot make
        // a repeated call free. `warm` is the same call again on the same engine, which is what a window pays per frame.
        let fresh = || -> Engine {
            let mut e = Engine::new();
            e.call("project_open", &json!({ "path": f })).expect("project_open");
            if let Some(m) = a.iter().position(|x| x == "--mesh").and_then(|i| a.get(i + 1)) {
                e.call("set_mesh", &serde_json::from_str(m).expect("--mesh is JSON")).expect("set_mesh");
            }
            e
        };
        let cold = |n: usize, call: &str| -> (f64, f64) {
            stats(
                (0..n)
                    .map(|_| {
                        let mut e = fresh();
                        let t = Instant::now();
                        std::hint::black_box(e.call(call, &json!({})).is_ok());
                        t.elapsed().as_secs_f64() * 1e3
                    })
                    .collect(),
            )
        };
        let geo = cold(reps, "geometry_names");
        let chk = cold(reps, "check_project");
        let mesh = cold(reps.min(10), "mesh_build");
        let warm = time(reps, || {
            std::hint::black_box(e.call("geometry_names", &json!({})).is_ok());
        });
        let _ = e.call("mesh_build", &json!({}));
        let (cells, dim) = e.model().map(|m| (m.mesh().n_cells(), m.mesh().dim() as u8)).unwrap_or((0, 0));
        let polar = if pol {
            let active = e.project().active.clone().or_else(|| e.project().physics.first().map(|p| p.physics.clone()));
            match active {
                // only the 2D electromagnetic examples are solved here (a 3D solve is seconds to minutes and has no polarization analysis)
                Some(p) if p == "ems" && dim == 2 => match e.solve_now(None) {
                    Ok(_) if e.polarization_applies(&p) => {
                        let t = time(reps, || {
                            std::hint::black_box(e.polarization(&p).is_ok());
                        });
                        format!("{:.1}/{:.1}", t.0, t.1)
                    }
                    Ok(_) => "n/a".into(),
                    Err(er) => format!("solve: {er}"),
                },
                _ => "-".into(),
            }
        } else {
            "-".into()
        };
        println!(
            "{:<28} {:>3} {:>8} | {:>7.2}/{:>8.2} | {:>7.2}/{:>8.2} | {:>7.1}/{:>8.1} | {:>16} (geometry warm p50 {:.3})",
            name, dim, cells, geo.0, geo.1, chk.0, chk.1, mesh.0, mesh.1, polar, warm.0
        );
    }
}
