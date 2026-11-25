// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Fast, bounded core-engine bug reproducers on unmodified `HEAD` (`9480093a`) for `squadds_bench bug-checks`
//! (`#41`): records explicit `status: "failed"` JSONL records with exact engine errors and distortion metrics
//! for Bug 1 (`tls.rs:789-797`), Bug 2 (`cavity.rs:646-694`), and Bug 3 (`sheet3d.rs` `z_coarsen: true`).

use crate::squadds_bench::schema::{
    parse_json_lenient_nan, InvariantsRecord, LayoutSpec, MemoryBytes, RunRecord, StageTimings,
};
use nyx_engine::facade::physics::CapacitanceMatrixData;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::mpsc;
use std::time::{Duration, Instant};

fn rect_poly(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<[f64; 2]> {
    vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
}

fn base_bug_record(case_id: &str, mode: &str) -> RunRecord {
    RunRecord {
        case_id: case_id.to_string(),
        suite_id: "BUG_CHECKS".to_string(),
        tier: "bug_checks".to_string(),
        mode: mode.to_string(),
        cluster_id: case_id.to_string(),
        stratum_weight: 1.0,
        status: "ok".to_string(),
        failure_detail: None,
        cells: 0,
        dofs: 0,
        order: 1,
        threads: femcore::par::threads(),
        timings_s: StageTimings::default(),
        memory_bytes: MemoryBytes::default(),
        invariants: InvariantsRecord::default(),
        predictions_si: BTreeMap::new(),
        references_si: BTreeMap::new(),
        rel_errors_pct: BTreeMap::new(),
        abs_errors_unit: BTreeMap::new(),
        extra: Value::Null,
    }
}

/// Bug 1 (`crates/nyx-analysis/src/tls.rs:789-797`):
/// On any finite-die 3D stack (`sub_box` confined inside outer vacuum), the 4 outer vertical `Si|Vac` die
/// sidewalls trigger `has_vertical_sa_trench = true` -> `has_resolved_3d_step = true`, clearing `sheet_edge_count`
/// and causing `FR-11` `analyse("capacitance_matrix", {"metal_thickness_m": 200e-9})` to fail.
fn check_bug1_tls_fr11_finite_die(inputs_dir: Option<&Path>) -> RunRecord {
    let t0 = Instant::now();
    let mut rec = base_bug_record("bug1_tls_fr11_finite_die", "fr11_finite_die");
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };

    // Load real `s1_1_smoke_00_r1063.layout.json` from `--inputs` if present, else fall back to synthetic finite-die CPW.
    let loaded_layout: Option<LayoutSpec> = inputs_dir.and_then(|dir| {
        let candidates = [
            dir.join("layouts").join("smoke").join("s1_1_smoke_00_r1063.layout.json"),
            dir.join("smoke").join("s1_1_smoke_00_r1063.layout.json"),
            dir.join("s1_1_smoke_00_r1063.layout.json"),
        ];
        for c in candidates {
            if let Ok(txt) = std::fs::read_to_string(&c) {
                if let Ok(lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                    return Some(lay);
                }
            }
        }
        None
    });

    let setup_res: Result<(u64, u64, String), String> = (|| {
        let (x0, y0, x1, y1, h_sub, eps_sub, pad_m, used_real) = if let Some(ref lay) = loaded_layout {
            let [bx0, by0, bx1, by1] = lay.die_bbox_m;
            let l_chip = lay.chip_length_m();
            (bx0, by0, bx1, by1, lay.stack.h_sub_m, lay.stack.eps_sub, 0.5 * l_chip, true)
        } else {
            (-150.0e-6, -150.0e-6, 150.0e-6, 150.0e-6, 300.0e-6, 11.45, 150.0e-6, false)
        };
        let h_base = ((x1 - x0).max(y1 - y0) / 8.5).clamp(60.0e-6, 145.0e-6);
        call(&mut e, "project_new", json!({ "title": "bug1_fr11_finite_die", "dim": 3, "h": h_base }))?;
        call(
            &mut e,
            "add_box",
            json!({ "name": "vac_bot", "min": [x0 - pad_m, y0 - pad_m, -h_sub - pad_m], "max": [x1 + pad_m, y1 + pad_m, -h_sub] }),
        )?;
        call(
            &mut e,
            "add_box",
            json!({ "name": "vac_mid", "min": [x0 - pad_m, y0 - pad_m, -h_sub], "max": [x1 + pad_m, y1 + pad_m, 0.0] }),
        )?;
        call(
            &mut e,
            "add_box",
            json!({ "name": "sub_box", "min": [x0 - 20.0e-6, y0 - 20.0e-6, -h_sub], "max": [x1 + 20.0e-6, y1 + 20.0e-6, 0.0] }),
        )?;
        call(
            &mut e,
            "add_box",
            json!({ "name": "vac_top", "min": [x0 - pad_m, y0 - pad_m, 0.0], "max": [x1 + pad_m, y1 + pad_m, pad_m] }),
        )?;

        let mut term_sels: Vec<(String, String)> = Vec::new();
        let mut sig_sels: Vec<String> = Vec::new();
        if let Some(ref lay) = loaded_layout {
            for sh in &lay.sheets {
                let p_name = format!("poly_{}", sh.name);
                call(
                    &mut e,
                    "add_polygons",
                    json!({ "name": p_name, "polygons": sh.polygons, "holes": sh.holes, "simplify_tol": sh.simplify_tol_m }),
                )?;
                call(&mut e, "add_sheet", json!({ "name": sh.name, "input": p_name, "z": sh.z_m }))?;
                let tn = sh.terminal.clone().unwrap_or_else(|| sh.name.clone());
                if !tn.to_ascii_lowercase().contains("ground") {
                    sig_sels.push(format!("{}:edge", sh.name));
                }
                term_sels.push((tn, sh.name.clone()));
            }
        } else {
            call(
                &mut e,
                "add_polygons",
                json!({ "name": "p_a", "polygons": [rect_poly(-60e-6, -60e-6, -10e-6, 60e-6)] }),
            )?;
            call(
                &mut e,
                "add_polygons",
                json!({ "name": "p_b", "polygons": [rect_poly(10e-6, -60e-6, 60e-6, 60e-6)] }),
            )?;
            call(&mut e, "add_sheet", json!({ "name": "pad_a", "input": "p_a", "z": 0.0 }))?;
            call(&mut e, "add_sheet", json!({ "name": "pad_b", "input": "p_b", "z": 0.0 }))?;
            sig_sels.push("pad_a:edge".into());
            sig_sels.push("pad_b:edge".into());
            term_sels.push(("pad_a".into(), "pad_a".into()));
            term_sels.push(("pad_b".into(), "pad_b".into()));
        }

        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                    { "name": "vacuum_mid", "input": "vac_mid", "kind": "Material" },
                    { "name": "substrate", "input": "sub_box", "kind": "Material" },
                    { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                ],
                "pins": [{ "sel": "outer", "id": 100 }]
            }),
        )?;
        for (idx, sel) in sig_sels.iter().enumerate() {
            call(&mut e, "add_size", json!({ "name": format!("sz_{idx}"), "sel": sel, "h": 14.0e-6 }))?;
        }
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": h_base,
                "structured_pml": false,
                "z_coarsen": false,
                "dz_min": 20.0e-6,
                "z_growth_ratio": 2.6,
                "max_z_slabs": 8,
                "max_cells": 250_000,
                "max_steiner_multiplier": 20,
                "strict_mesh": false,
            }),
        )?;
        let mrep = call(&mut e, "mesh_build", json!({ "wait": true }))?;
        let cells = mrep["mesh"]["cells"].as_u64().unwrap_or(0);

        call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }))?;
        call(
            &mut e,
            "set_physics_params",
            json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1e-9 } }),
        )?;
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": eps_sub } }),
        )?;
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }),
        )?;
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "vacuum_mid", "material": { "eps_r": 1.0 } }),
        )?;
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "vacuum_top", "material": { "eps_r": 1.0 } }),
        )?;
        for (tn, sel) in &term_sels {
            call(
                &mut e,
                "add_condition",
                json!({ "physics": "dfs", "on": sel, "kind": "terminal", "params": { "name": tn } }),
            )?;
        }
        let cx = 0.5 * (x0 + x1);
        let cy = 0.5 * (y0 + y1);
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2, "center": [cx, cy, 0.0] } }),
        )?;
        call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": true }))?;
        let (res, _) = e.result("dfs").ok_or("missing dfs result")?;
        let dofs = res.dofs as u64;
        let source = if used_real { "s1_1_smoke_00_r1063.layout.json" } else { "synthetic_finite_die" };
        Ok((cells, dofs, source.to_string()))
    })();

    match setup_res {
        Ok((cells, dofs, source)) => {
            rec.cells = cells;
            rec.dofs = dofs;
            match e.call(
                "analyse",
                &json!({
                    "analysis": "capacitance_matrix",
                    "physics": "dfs",
                    "params": { "metal_thickness_m": 200.0e-9 }
                }),
            ) {
                Ok(out) => {
                    let has_delta = out["thickness_correction"]["delta_maxwell_F"].is_array();
                    if has_delta {
                        rec.status = "ok".into();
                        rec.extra = json!({
                            "bug_id": "Bug 1 (tls.rs:789-797)",
                            "source_layout": source,
                            "reproduced": false,
                            "thickness_correction": out["thickness_correction"],
                        });
                    } else {
                        rec.status = "failed".into();
                        rec.failure_detail =
                            Some("FR-11 metal_thickness_m returned no delta_maxwell_F on finite die".into());
                        rec.extra =
                            json!({ "bug_id": "Bug 1 (tls.rs:789-797)", "source_layout": source, "reproduced": true });
                    }
                }
                Err(err) => {
                    let msg = err.to_string();
                    rec.status = "failed".into();
                    rec.failure_detail =
                        Some(format!("Bug 1 reproduced (crates/nyx-analysis/src/tls.rs:789-797): {msg}"));
                    rec.predictions_si.insert("bug_reproduced".into(), 1.0);
                    rec.extra = json!({
                        "bug_id": "Bug 1 (tls.rs:789-797)",
                        "source_layout": source,
                        "reproduced": true,
                        "engine_error": msg,
                    });
                }
            }
        }
        Err(err) => {
            rec.status = "failed".into();
            rec.failure_detail = Some(format!("bug1 setup error: {err}"));
        }
    }

    rec.timings_s.total = t0.elapsed().as_secs_f64();
    if let Some((ws, commit)) = nyx_engine::process::peak_memory() {
        rec.memory_bytes = MemoryBytes { peak_ws: ws, peak_commit: commit };
    }
    rec
}

/// Runs a small 3D planar CPW cavity with a `lumped_port` (`20 x 20 um`, `Z_0 = 50 Ohm`), optionally refining the
/// port sheet edge to `h = 6 um` (creating interior non-PEC vertices on its lateral vacuum edges).
fn run_bug2_cavity_case(refine_port: bool) -> Result<(u64, u64, f64, Vec<f64>), String> {
    let t0 = Instant::now();
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };

    // Compact 800 x 500 um cavity box with a 600 um long CPW strip and 20 x 20 um axial lumped_port
    call(&mut e, "project_new", json!({ "title": "bug2_cavity_port", "dim": 3, "h": 95.0e-6 }))?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "sub_box", "min": [-400e-6, -250e-6, -300e-6], "max": [400e-6, 250e-6, 0.0] }),
    )?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "vac_top", "min": [-400e-6, -250e-6, 0.0], "max": [400e-6, 250e-6, 300e-6] }),
    )?;

    // Signal trace `x in [-300, 300] um, y in [-10, 10] um`; disjoint ground return bar `x in [320, 400] um, y in [-250, 250] um`
    // connected to `sig` across a `20 x 20 um` axial `lumped_port` (`x in [300, 320], y in [-10, 10]`).
    call(&mut e, "add_polygons", json!({ "name": "p_sig", "polygons": [rect_poly(-300e-6, -10e-6, 300e-6, 10e-6)] }))?;
    call(&mut e, "add_polygons", json!({ "name": "p_gnd", "polygons": [rect_poly(320e-6, -250e-6, 400e-6, 250e-6)] }))?;
    call(&mut e, "add_polygons", json!({ "name": "p_port", "polygons": [rect_poly(300e-6, -10e-6, 320e-6, 10e-6)] }))?;
    call(&mut e, "add_sheet", json!({ "name": "sig", "input": "p_sig", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "gnd", "input": "p_gnd", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "sh_port", "input": "p_port", "z": 0.0 }))?;

    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "substrate", "input": "sub_box", "kind": "Material" },
                { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
            ]
        }),
    )?;
    if refine_port {
        // Subdivide the 20 x 20 um port rectangle edges to 5 um -> creates interior non-PEC vertices on lateral vacuum edges
        call(&mut e, "add_size", json!({ "name": "sz_port", "sel": "sh_port:edge", "h": 5.0e-6 }))?;
    }
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 95.0e-6,
            "structured_pml": false,
            "z_coarsen": false,
            "dz_min": 30.0e-6,
            "z_growth_ratio": 2.4,
            "max_z_slabs": 8,
            "max_cells": 60_000,
            "max_steiner_multiplier": 6,
            "strict_mesh": false,
        }),
    )?;
    let mrep = call(&mut e, "mesh_build", json!({ "wait": true }))?;
    let cells = mrep["mesh"]["cells"].as_u64().unwrap_or(0);

    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 1,
            "study_params": { "n_modes": 4, "f_search": 5.0e9 }
        }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "substrate", "material": { "eps_r": 11.45 } }))?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vacuum_top", "material": { "eps_r": 1.0 } }))?;
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "sig", "kind": "pec", "params": {} }))?;
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "gnd", "kind": "pec", "params": {} }))?;
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "outer", "kind": "pec", "params": {} }))?;
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "sh_port",
            "kind": "lumped_port",
            "params": { "name": "port_1", "Z_0": 50.0, "direction": [1.0, 0.0, 0.0] }
        }),
    )?;
    call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true, "isolate": false, "save_fields": false }))?;
    let (res, _) = e.result("ems").ok_or("missing ems result")?;
    let freqs: Vec<f64> = res.solutions.iter().filter_map(|s| s.value).collect();
    Ok((cells, res.dofs as u64, t0.elapsed().as_secs_f64(), freqs))
}

/// Bug 2 (`crates/nyx-ems/src/cavity.rs:646-694`):
/// Unconditionally adding all vertices of a `lumped_port`/`lumped_jj` sheet to `cond_vertices`/`pec_vertex`
/// drops interior lateral-vacuum-edge vertex gradient columns from `g_red`, spawning spurious ~0 Hz modes that
/// stall or non-converge `krylov_schur`. Bounded via 20 s worker timeout so the check never hangs.
fn check_bug2_cavity_refined_lumped_port() -> RunRecord {
    let t0 = Instant::now();
    let mut rec = base_bug_record("bug2_cavity_refined_lumped_port", "ems_cavity_port_refine");

    let base_res = run_bug2_cavity_case(false);
    let (base_cells, base_dofs, base_s, base_freqs) = match base_res {
        Ok(v) => v,
        Err(err) => {
            rec.status = "failed".into();
            rec.failure_detail = Some(format!("bug2 baseline unrefined port failed unexpectedly: {err}"));
            rec.timings_s.total = t0.elapsed().as_secs_f64();
            return rec;
        }
    };

    // Run the refined-port case on a worker thread with a hard 20 s bound so a Krylov-Schur stall never hangs the harness.
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let r = run_bug2_cavity_case(true);
        let _ = tx.send(r);
    });

    let timeout_s = 20.0_f64;
    match rx.recv_timeout(Duration::from_secs(20)) {
        Err(_) => {
            rec.status = "failed".into();
            rec.cells = base_cells;
            rec.dofs = base_dofs;
            rec.failure_detail = Some(format!(
                "Bug 2 reproduced (crates/nyx-ems/src/cavity.rs:646-694): refined lumped_port stalled krylov_schur > {timeout_s:.0}s (vs unrefined baseline {base_s:.2}s)"
            ));
            rec.predictions_si.insert("baseline_solve_s".into(), base_s);
            rec.predictions_si.insert("refined_timeout_s".into(), timeout_s);
            rec.predictions_si.insert("bug_reproduced".into(), 1.0);
            rec.extra = json!({
                "bug_id": "Bug 2 (cavity.rs:646-694)",
                "reproduced": true,
                "stalled_timeout": true,
                "baseline_cells": base_cells,
                "baseline_solve_s": base_s,
                "baseline_freqs_hz": base_freqs,
            });
        }
        Ok(Err(err)) => {
            rec.status = "failed".into();
            rec.cells = base_cells;
            rec.dofs = base_dofs;
            rec.failure_detail = Some(format!(
                "Bug 2 reproduced (crates/nyx-ems/src/cavity.rs:646-694): refined lumped_port solver error: {err} (vs unrefined baseline {base_s:.2}s)"
            ));
            rec.predictions_si.insert("baseline_solve_s".into(), base_s);
            rec.predictions_si.insert("bug_reproduced".into(), 1.0);
            rec.extra = json!({
                "bug_id": "Bug 2 (cavity.rs:646-694)",
                "reproduced": true,
                "engine_error": err,
                "baseline_cells": base_cells,
                "baseline_solve_s": base_s,
                "baseline_freqs_hz": base_freqs,
            });
        }
        Ok(Ok((ref_cells, ref_dofs, ref_s, ref_freqs))) => {
            rec.cells = ref_cells;
            rec.dofs = ref_dofs;
            let spurious_count = ref_freqs.iter().filter(|&&f| f < 1.0e9).count();
            let slowdown = ref_s / base_s.max(0.05);
            rec.predictions_si.insert("baseline_solve_s".into(), base_s);
            rec.predictions_si.insert("refined_solve_s".into(), ref_s);
            rec.predictions_si.insert("slowdown_ratio".into(), slowdown);
            rec.predictions_si.insert("spurious_modes_lt_1ghz".into(), spurious_count as f64);
            if spurious_count > 0 || slowdown > 8.0 {
                rec.status = "failed".into();
                rec.failure_detail = Some(format!(
                    "Bug 2 reproduced (crates/nyx-ems/src/cavity.rs:646-694): refined lumped_port produced {spurious_count} spurious <1 GHz mode(s) and {slowdown:.1}x slowdown ({ref_s:.2}s vs {base_s:.2}s)"
                ));
                rec.predictions_si.insert("bug_reproduced".into(), 1.0);
            } else {
                rec.status = "ok".into();
                rec.predictions_si.insert("bug_reproduced".into(), 0.0);
            }
            rec.extra = json!({
                "bug_id": "Bug 2 (cavity.rs:646-694)",
                "reproduced": rec.status == "failed",
                "baseline_cells": base_cells,
                "baseline_solve_s": base_s,
                "baseline_freqs_hz": base_freqs,
                "refined_cells": ref_cells,
                "refined_solve_s": ref_s,
                "refined_freqs_hz": ref_freqs,
                "spurious_modes_lt_1ghz": spurious_count,
                "slowdown_ratio": slowdown,
            });
        }
    }

    rec.timings_s.total = t0.elapsed().as_secs_f64();
    if let Some((ws, commit)) = nyx_engine::process::peak_memory() {
        rec.memory_bytes = MemoryBytes { peak_ws: ws, peak_commit: commit };
    }
    rec
}

fn run_bug3_z_coarsen_case(_lay: Option<&LayoutSpec>, z_coarsen: bool) -> Result<(u64, u64, f64), String> {
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Result<Value, String> {
        eng.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
    };
    // Multi-layer finite-die TransmonCross+Claw stack (`P = 1.2 L`, no `inner_vac_box` workaround):
    // `vac_box` + `sub_box` with `z_coarsen: true` triggers `z_max` cap collapse
    // (`boundary number(s) [...] do not exist`) and/or >20% capacitance distortion on unmodified HEAD.
    let half = 500.0e-6;
    let pad = 1200.0e-6;
    let h_sub = 750.0e-6;
    let h_air = 890.0e-6;
    let h_bot = 900.0e-6;
    let h_global = 340.0e-6;
    call(&mut e, "project_new", json!({ "title": "bug3_z_coarsen_c02", "dim": 3, "h": h_global }))?;
    call(
        &mut e,
        "add_box",
        json!({
            "name": "vac_box",
            "min": [-half - pad, -half - pad, -h_sub - h_bot - pad],
            "max": [half + pad, half + pad, h_air + pad]
        }),
    )?;
    call(&mut e, "add_box", json!({ "name": "sub_box", "min": [-half, -half, -h_sub], "max": [half, half, 0.0] }))?;
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "cross_poly",
            "polygons": [[[-15e-6, -120e-6], [15e-6, -120e-6], [15e-6, -15e-6], [120e-6, -15e-6], [120e-6, 15e-6], [15e-6, 15e-6], [15e-6, 120e-6], [-15e-6, 120e-6], [-15e-6, 15e-6], [-120e-6, 15e-6], [-120e-6, -15e-6], [-15e-6, -15e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "claw_poly",
            "polygons": [[[-45e-6, 132e-6], [45e-6, 132e-6], [45e-6, 190e-6], [30e-6, 190e-6], [30e-6, 147e-6], [-30e-6, 147e-6], [-30e-6, 190e-6], [-45e-6, 190e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "gnd_poly",
            "polygons": [[[-half, -half], [half, -half], [half, half], [-half, half]]],
            "holes": [[[-140e-6, -140e-6], [140e-6, -140e-6], [140e-6, 210e-6], [-140e-6, 210e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "cross", "input": "cross_poly", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "claw", "input": "claw_poly", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "gnd", "input": "gnd_poly", "z": 0.0 }))?;
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "Vac", "input": "vac_box", "kind": "Material" },
                { "name": "Si", "input": "sub_box", "kind": "Material" }
            ]
        }),
    )?;
    call(&mut e, "add_size", json!({ "name": "sz_c", "sel": "cross:edge", "h": 10.0e-6 }))?;
    call(&mut e, "add_size", json!({ "name": "sz_k", "sel": "claw:edge", "h": 8.0e-6 }))?;
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": h_global,
            "structured_pml": false,
            "z_coarsen": z_coarsen,
            "dz_min": 5.0e-6,
            "z_growth_ratio": 1.6,
            "max_z_slabs": 22,
            "max_cells": 350_000,
            "max_steiner_multiplier": 8,
            "strict_mesh": false,
        }),
    )?;
    let mrep = call(&mut e, "mesh_build", json!({ "wait": true }))?;
    let cells = mrep["mesh"]["cells"].as_u64().unwrap_or(0);

    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }))?;
    call(
        &mut e,
        "set_physics_params",
        json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1e-9 } }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }))?;
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "cross", "kind": "terminal", "params": { "name": "cross" } }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "claw", "kind": "terminal", "params": { "name": "claw" } }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "gnd", "kind": "terminal", "params": { "name": "gnd" } }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2, "center": [0.0, 0.0, 0.0] } }),
    )?;
    call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": false }))?;
    let (res, _) = e.result("dfs").ok_or("missing dfs result")?;
    let cap = CapacitanceMatrixData::of(res).ok_or("missing CapacitanceMatrixData")?;
    let c_diag = cap.maxwell_matrix_f[1][1];
    Ok((cells, res.dofs as u64, c_diag))
}

/// Bug 3 (`crates/femcore/src/meshing/sheet3d.rs` `z_coarsen: true`):
/// Multi-layer chip stack with `z_coarsen: true` triggers `z_max` cap collapse
/// (`boundary number(s) [...] do not exist`) and/or `> 20%` coplanar capacitance distortion.
fn check_bug3_sheet3d_z_coarsen(inputs_dir: Option<&Path>) -> RunRecord {
    let t0 = Instant::now();
    let mut rec = base_bug_record("bug3_sheet3d_z_coarsen", "sheet3d_z_coarsen");

    let loaded_layout: Option<LayoutSpec> = inputs_dir.and_then(|dir| {
        let candidates = [
            dir.join("layouts").join("smoke").join("s1_1_smoke_00_r1063.layout.json"),
            dir.join("smoke").join("s1_1_smoke_00_r1063.layout.json"),
            dir.join("s1_1_smoke_00_r1063.layout.json"),
        ];
        for c in candidates {
            if let Ok(txt) = std::fs::read_to_string(&c) {
                if let Ok(lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                    return Some(lay);
                }
            }
        }
        None
    });

    let ref_res = run_bug3_z_coarsen_case(loaded_layout.as_ref(), false);
    let on_res = run_bug3_z_coarsen_case(loaded_layout.as_ref(), true);

    let mut issues: Vec<String> = Vec::new();
    let mut c_ref_f: Option<f64> = None;
    if let Ok((cells, dofs, c_ref)) = ref_res {
        rec.cells = cells;
        rec.dofs = dofs;
        c_ref_f = Some(c_ref);
        rec.predictions_si.insert("c_diag_z_coarsen_false_f".into(), c_ref);
    } else if let Err(err) = &ref_res {
        issues.push(format!("baseline z_coarsen=false error: {err}"));
    }

    let mut on_err: Option<String> = None;
    let mut delta_pct: Option<f64> = None;
    match on_res {
        Ok((cells_on, _dofs_on, c_on)) => {
            rec.predictions_si.insert("c_diag_z_coarsen_true_f".into(), c_on);
            rec.predictions_si.insert("cells_z_coarsen_true".into(), cells_on as f64);
            if let Some(cr) = c_ref_f {
                let dpct = 100.0 * (c_on - cr).abs() / cr.max(1.0e-30);
                delta_pct = Some(dpct);
                rec.predictions_si.insert("z_coarsen_rel_distortion_pct".into(), dpct);
                if dpct > 20.0 {
                    issues.push(format!(
                        "z_coarsen=true distorted coplanar capacitance by {dpct:.1}% ({:.2} fF vs {:.2} fF baseline)",
                        c_on * 1e15,
                        cr * 1e15
                    ));
                }
            }
        }
        Err(err) => {
            on_err = Some(err.clone());
            issues.push(format!("z_coarsen=true failed (cap collapse / missing boundary): {err}"));
        }
    }

    if !issues.is_empty() {
        rec.status = "failed".into();
        rec.failure_detail = Some(format!(
            "Bug 3 reproduced (crates/femcore/src/meshing/sheet3d.rs z_coarsen=true): {}",
            issues.join("; ")
        ));
        rec.predictions_si.insert("bug_reproduced".into(), 1.0);
    } else {
        rec.status = "ok".into();
        rec.predictions_si.insert("bug_reproduced".into(), 0.0);
    }
    rec.extra = json!({
        "bug_id": "Bug 3 (sheet3d.rs z_coarsen=true)",
        "reproduced": rec.status == "failed",
        "delta_pct": delta_pct,
        "engine_error": on_err,
        "issues": issues,
    });
    rec.timings_s.total = t0.elapsed().as_secs_f64();
    if let Some((ws, commit)) = nyx_engine::process::peak_memory() {
        rec.memory_bytes = MemoryBytes { peak_ws: ws, peak_commit: commit };
    }
    rec
}

/// Executes all 3 bounded unmodified-HEAD core-engine bug reproducers and returns their `RunRecord`s.
pub fn run_all_bug_checks(inputs_dir: Option<&Path>) -> Vec<RunRecord> {
    let r1 = check_bug1_tls_fr11_finite_die(inputs_dir);
    println!(
        "[bug-check] {} status={} t={:.2}s detail={}",
        r1.case_id,
        r1.status,
        r1.timings_s.total,
        r1.failure_detail.as_deref().unwrap_or("none")
    );
    let r2 = check_bug2_cavity_refined_lumped_port();
    println!(
        "[bug-check] {} status={} t={:.2}s detail={}",
        r2.case_id,
        r2.status,
        r2.timings_s.total,
        r2.failure_detail.as_deref().unwrap_or("none")
    );
    let r3 = check_bug3_sheet3d_z_coarsen(inputs_dir);
    println!(
        "[bug-check] {} status={} t={:.2}s detail={}",
        r3.case_id,
        r3.status,
        r3.timings_s.total,
        r3.failure_detail.as_deref().unwrap_or("none")
    );
    vec![r1, r2, r3]
}
