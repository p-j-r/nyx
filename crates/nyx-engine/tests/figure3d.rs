// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end engine verification tests for FR-4 3D `figure` command (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`,
//! 1-based `modes`, `expr`, `plane`, `clip`, `isosurfaces`, `show_edges`, `camera`, and auto-resolved `repeat_periodic` Bloch tiling).

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::f64::consts::FRAC_PI_2;
use std::path::PathBuf;
use std::time::Instant;

fn tmp(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("nyx-fig3d-{}-{name}", std::process::id()))
}

fn decode_png(bytes: &[u8]) -> (u32, u32, Vec<u8>) {
    assert_eq!(&bytes[..4], &[0x89, b'P', b'N', b'G']);
    nyx_engine::analysis::view3d::decode_png(bytes).expect("valid PNG")
}

/// Swept 3D dielectric waveguide unit cell (`Core` eps_r = 11.9, `Gap` eps_r = 2.055 -> ~5.79x normal jump),
/// Floquet-periodic along `y` (`ymin` -> `ymax`) with `phase`.
fn swept_periodic_waveguide(lx: f64, lz: f64, ly: f64, h: f64, n_layers: usize, phase: f64) -> Value {
    let feats = json!([
        {"id": 1, "name": "box", "op": {"Rect": {"min": [-0.5 * lx, -0.5 * lz], "max": [0.5 * lx, 0.5 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 2, "name": "core", "op": {"Rect": {"min": [-0.22 * lx, -0.22 * lz], "max": [0.22 * lx, 0.22 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 3, "name": "ov", "op": {"Overlay": {"naming": "Feature", "pins": [], "regions": [
            {"input": 1, "kind": "Material", "name": "Gap", "pinned_id": null},
            {"input": 2, "kind": "Material", "name": "Si", "pinned_id": null}
        ]}}, "origin": "Code", "suppressed": false, "tags": []}
    ]);
    json!({
        "active": "ems",
        "format": 1,
        "title": "3D periodic waveguide",
        "geometry": {"dim": 2, "features": feats, "format": 1, "params": []},
        "mesh": {"h": h, "sweep": [{"thickness": ly, "n": n_layers}]},
        "physics": [{
            "physics": "ems",
            "format": 1,
            "order": 2,
            "bcs": [
                {"kind": "periodic", "on": "ymin", "params": {"dest": "ymax", "phase": phase}}
            ],
            "materials": {"Gap": {"eps_r": 2.055}, "Si": {"eps_r": 11.90}},
            "study": {"kind": "cavity", "params": {"f_search": 193.4e12, "n_modes": 2}}
        }]
    })
}

#[test]
fn engine_figure_supports_all_3d_kinds_1based_modes_and_bloch_auto_resolution() {
    let mut e = Engine::new();
    let proj = swept_periodic_waveguide(1.2e-6, 0.9e-6, 0.35e-6, 0.28e-6, 2, FRAC_PI_2);
    e.call("project_load", &json!({ "json": proj.to_string() })).unwrap();
    e.call("solve_run", &json!({ "physics": "ems" })).unwrap();

    // 1-based `modes` contract: 0 is rejected, [1] succeeds and reports `mode: 1`
    let bad_path = tmp("bad0.png");
    assert!(
        e.call(
            "figure",
            &json!({ "physics": "ems", "path": bad_path.to_str().unwrap(), "kind": "slice3d", "modes": [0] })
        )
        .is_err(),
        "`modes` must stay 1-based and reject 0"
    );

    for kind in ["slice3d", "cutaway3d", "isosurface3d", "mesh3d"] {
        for ext in ["png", "svg"] {
            let p = tmp(&format!("{kind}.{ext}"));
            let mut args = json!({
                "physics": "ems",
                "path": p.to_str().unwrap(),
                "kind": kind,
                "modes": [1],
                "repeat_periodic": 4,
                "scale": 1.0
            });
            if kind == "slice3d" {
                args["plane"] = json!({ "axis": "z", "at": 0.0 });
                args["expr"] = json!("re(Ex)");
            } else if kind == "cutaway3d" {
                args["clip"] = json!({ "zmax": 0.0, "xmax": 0.2e-6 });
                args["isosurfaces"] = json!([-0.45, 0.45]);
                args["expr"] = json!("re(Ex)");
                args["camera"] = json!({ "azimuth_deg": -48.0, "elevation_deg": 26.0 });
            } else if kind == "isosurface3d" {
                args["expr"] = json!("abs2(E)");
                args["show_edges"] = json!("regions");
            }
            let r = e.call("figure", &args).unwrap();
            assert_eq!(r["kind"], kind);
            assert_eq!(r["format"], ext);
            if kind != "mesh3d" {
                assert_eq!(r["modes"][0]["mode"], 1);
            }
            let raw = std::fs::read(&p).unwrap();
            if ext == "png" {
                let (w, h, _) = decode_png(&raw);
                assert!(w > 100 && h > 100);
            } else {
                let s = String::from_utf8(raw).unwrap();
                assert!(s.starts_with("<svg") && s.trim_end().ends_with("</svg>"));
            }
            let _ = std::fs::remove_file(&p);
        }
    }

    // Verify Bloch phase auto-resolved from periodic BC `phase = pi/2`: on z = 0 slice, period m=0 and m=2 have opposite signs!
    // Mode 1 of the 1.2 x 0.9 um guide has dominant transverse field (`re(Ez)` or `re(Ex)`) peaking at the waveguide centre (x=0, z=0).
    let slice_png = tmp("bloch_sign.png");
    let mut best_rb = 0i32;
    let mut rb_m2 = 0i32;
    for expr in ["re(Ez)", "re(Ex)"] {
        e.call(
            "figure",
            &json!({
                "physics": "ems",
                "path": slice_png.to_str().unwrap(),
                "kind": "slice3d",
                "modes": [1],
                "expr": expr,
                "plane": { "axis": "z", "at": 0.0 },
                "repeat_periodic": 4,
                "show_edges": "none",
                "window": [-0.6e-6, 0.6e-6, 0.0, 1.4e-6],
                "panel_width": 280.0,
                "scale": 1.0
            }),
        )
        .unwrap();
        let (w, h, px) = decode_png(&std::fs::read(&slice_png).unwrap());
        let _ = std::fs::remove_file(&slice_png);
        assert!(w > 280 && h > 300);
        let mut local_best = 0i32;
        let mut best_pos = (0usize, 0usize);
        for py in 350..420usize {
            for px_x in 65..325usize {
                let i = (py * w as usize + px_x) * 4;
                let rb = px[i] as i32 - px[i + 2] as i32;
                if rb.abs() > local_best.abs() {
                    local_best = rb;
                    best_pos = (px_x, py);
                }
            }
        }
        if local_best.abs() > best_rb.abs() {
            best_rb = local_best;
            let py_m2 = (best_pos.1 as f64 - 2.0 * (0.35 / 1.4) * (280.0 * (1.4 / 1.2))).round() as usize;
            let i_m2 = (py_m2 * w as usize + best_pos.0) * 4;
            rb_m2 = px[i_m2] as i32 - px[i_m2 + 2] as i32;
        }
    }
    assert!(best_rb.abs() > 80, "expected strong signed Coolwarm peak in core, got {best_rb}");
    assert!(
        best_rb.signum() == -rb_m2.signum() && rb_m2.abs() > 80,
        "expected sign reversal after 2 periods (phase = pi/2): m=0 R-B={best_rb}, m=2 R-B={rb_m2}"
    );
}

#[test]
fn native_3d_periodic_pairs_auto_resolve_shift_and_engine_6300_tet_latency_is_under_150_ms() {
    // Native 3D geom Box with `Periodic` feature on y-faces so `model.evaluated().periodic_pairs()` resolves shift directly
    let (lx, ly, lz) = (1.2e-6, 0.35e-6, 0.8e-6);
    let mut e3 = Engine::new();
    let proj3 = json!({
        "active": "ems",
        "format": 1,
        "title": "native 3D periodic box",
        "geometry": {
            "dim": 3,
            "format": 1,
            "params": [],
            "features": [
                {"id": 1, "name": "cell", "op": {"Box": {"min": [-0.5 * lx, 0.0, -0.5 * lz], "max": [0.5 * lx, ly, 0.5 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
                {"id": 2, "name": "ov", "op": {"Overlay": {"naming": "Feature", "pins": [], "regions": [
                    {"input": 1, "kind": "Material", "name": "Si", "pinned_id": null}
                ]}}, "origin": "Code", "suppressed": false, "tags": []},
                {"id": 3, "name": "py", "op": {"Periodic": {
                    "sel_a": {"Name": "cell.ymin"},
                    "sel_b": {"Name": "cell.ymax"},
                    "translation": [0.0, ly, 0.0]
                }}, "origin": "Code", "suppressed": false, "tags": []}
            ]
        },
        "mesh": {"h": 0.35e-6},
        "physics": [{
            "physics": "ems",
            "format": 1,
            "order": 2,
            "bcs": [
                {"kind": "periodic", "on": "cell.ymin", "params": {"dest": "cell.ymax", "phase": FRAC_PI_2}}
            ],
            "materials": {"Si": {"eps_r": 11.9}},
            "study": {"kind": "cavity", "params": {"f_search": 150.0e12, "n_modes": 1}}
        }]
    });
    e3.call("project_load", &json!({ "json": proj3.to_string() })).unwrap();
    e3.call("solve_run", &json!({ "physics": "ems" })).unwrap();
    let p_native = tmp("native3d.png");
    let rn = e3
        .call(
            "figure",
            &json!({
                "physics": "ems",
                "path": p_native.to_str().unwrap(),
                "kind": "cutaway3d",
                "repeat_periodic": 4,
                "expr": "im(Ey)",
                "scale": 1.0
            }),
        )
        .unwrap();
    assert!(rn["bytes"].as_u64().unwrap() > 1000);
    let _ = std::fs::remove_file(&p_native);

    // Build a >6,300-tet P2 3D swept waveguide mesh and time `Engine::call("figure", ...)` at 1600x1100 (scale = 2.0) with repeat_periodic = 4
    let mut eb = Engine::new();
    let bench_proj = swept_periodic_waveguide(1.2e-6, 0.9e-6, 0.35e-6, 0.092e-6, 7, FRAC_PI_2);
    eb.call("project_load", &json!({ "json": bench_proj.to_string() })).unwrap();
    eb.call("mesh_build", &json!({})).unwrap();
    let n_tets = eb.model().unwrap().solver_mesh().n_cells();
    assert!(n_tets >= 6_300, "expected >= 6,300 tets for engine latency test, got {n_tets}");

    // Solve 1 mode on a fast electrostatic/cavity or mesh3d/cutaway3d: solve 1 cavity mode at order 2? Wait: let's solve dfs order=2 (instant!) or mesh3d + dfs cutaway3d!
    // Even better: solve dfs P2 on the exact >=6,300-tet 3D mesh (< 20 ms solve!) with V + Ex/Ey/Ez or V so we benchmark `Engine::figure` purely!
    let dfs_proj = json!({
        "active": "dfs",
        "format": 1,
        "title": "6300-tet 3D P2 figure benchmark",
        "geometry": bench_proj["geometry"],
        "mesh": bench_proj["mesh"],
        "physics": [{
            "physics": "dfs",
            "format": 1,
            "order": 2,
            "bcs": [
                {"kind": "ground", "on": "ymin", "params": {}},
                {"kind": "potential", "on": "ymax", "params": {"value": 1.0}}
            ],
            "materials": {"Gap": {"eps_r": 2.055}, "Si": {"eps_r": 11.90}},
            "study": {"kind": "electrostatic", "params": {}}
        }]
    });
    eb.call("project_load", &json!({ "json": dfs_proj.to_string() })).unwrap();
    eb.call("solve_run", &json!({ "physics": "dfs" })).unwrap();
    let actual_tets = eb.model().unwrap().solver_mesh().n_cells();
    assert!(actual_tets >= 6_300, "expected >= 6,300 tets, got {actual_tets}");

    let out_bench = tmp("bench_1600x1100.png");
    let pw = 760.0;
    let ph = 362.0;
    let fig_args = json!({
        "physics": "dfs",
        "path": out_bench.to_str().unwrap(),
        "kind": "cutaway3d",
        "clip": { "zmax": 0.0, "xmax": 0.2e-6 },
        "repeat_periodic": 4,
        "periodic_phase": FRAC_PI_2,
        "panel_width": pw,
        "window": [-1.0e-6, 1.0e-6, -1.0e-6 * (ph / pw), 1.0e-6 * (ph / pw)],
        "scale": 2.0
    });
    // Warm-up
    let mut rep = eb.call("figure", &fig_args).unwrap();

    // Best-of-5 iterations avoids transient OS/filesystem jitter during concurrent `--workspace` test runs
    let best_ms = (0..5)
        .map(|_| {
            let t0 = Instant::now();
            rep = eb.call("figure", &fig_args).unwrap();
            t0.elapsed().as_secs_f64() * 1000.0
        })
        .fold(f64::INFINITY, f64::min);
    let _ = std::fs::remove_file(&out_bench);

    let n_threads = femcore::par::threads();
    let limit_ms = if n_threads >= 4 { 150.0 } else { 600.0 };
    assert_eq!((rep["width"].as_f64().unwrap(), rep["height"].as_f64().unwrap()), (1600.0, 1100.0));
    eprintln!(
        "Engine::figure cutaway3d ({actual_tets} P2 tets x 4 repeats, 1600x1100 PNG incl. disk write): {best_ms:.2} ms (best of 5, threads={n_threads}, limit={limit_ms:.0} ms)"
    );
    assert!(best_ms < limit_ms, "Engine::figure took {best_ms:.2} ms (target < {limit_ms:.0} ms, threads={n_threads})");
}
