// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 5 VVUQ verification tests for deterministic HDF5/archive provenance,
//! `/provenance/materials_json` drift warnings on `open_results`, and `RunRecord` v1.1 persistence.
#![cfg(feature = "hdf5")]

use hdf5_pure::{AttrValue, File};
use nyx_engine::provenance::{build_provenance, resolved_nyx_factor_threads, resolved_nyx_threads};
use nyx_engine::store_link::{RunOutcome, StoreLink};
use nyx_engine::Engine;
use nyx_physics::DiagnosticSeverity;
use nyx_store::archive::{export_bytes, import, inspect, ExportOptions, Limits};
use nyx_store::runs::RunStatus;
use nyx_store::{ShareLevel, WorkflowDir};
use serde_json::{json, Value};
use std::io::Read;
use std::path::{Path, PathBuf};
use zip::ZipArchive;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn tmp(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-vvuq-p5-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn attr_string(attrs: &std::collections::HashMap<String, AttrValue>, key: &str) -> String {
    let v = attrs.get(key).unwrap_or_else(|| panic!("missing HDF5 attr '{key}' in {attrs:?}"));
    v.as_str().map(str::to_string).or_else(|| v.as_strings().and_then(|a| a.first().cloned())).expect("string attr")
}

fn assert_no_created_utc_anywhere(path: &Path) {
    let file = File::open_streaming(path).unwrap();
    let root_attrs = file.root().attrs().unwrap();
    assert!(!root_attrs.contains_key("created_utc"), "root must not contain created_utc");
    if let Ok(pg) = file.group("provenance") {
        let p_attrs = pg.attrs().unwrap();
        assert!(!p_attrs.contains_key("created_utc"), "/provenance must not contain created_utc");
    }
    if let Ok(rg) = file.group("results") {
        for code in rg.groups().unwrap() {
            let g = rg.group(&code).unwrap();
            let g_attrs = g.attrs().unwrap();
            assert!(!g_attrs.contains_key("created_utc"), "/results/{code} must not contain created_utc");
        }
    }
}

#[test]
fn repeated_save_results_and_nyx_archives_are_byte_identical_and_omit_created_utc() {
    let root = tmp("deterministic");
    let ws = root.join("ws");
    std::fs::create_dir_all(&ws).unwrap();

    let mut e = Engine::new();
    e.set_workspace(ws.clone(), true);
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));

    let h5_a = root.join("export_a.h5");
    let h5_b = root.join("export_b.h5");
    call(&mut e, "export_results", json!({ "path": h5_a.to_str().unwrap() }));
    std::thread::sleep(std::time::Duration::from_millis(25));
    call(&mut e, "export_results", json!({ "path": h5_b.to_str().unwrap() }));

    let bytes_a = std::fs::read(&h5_a).unwrap();
    let bytes_b = std::fs::read(&h5_b).unwrap();
    assert_eq!(bytes_a, bytes_b, "repeated export_results must produce 0-byte-diff HDF5");
    assert_no_created_utc_anywhere(&h5_a);

    // Check deterministic root provenance attributes.
    let f = File::open_streaming(&h5_a).unwrap();
    let attrs = f.root().attrs().unwrap();
    assert_eq!(attr_string(&attrs, "nyx_version"), env!("CARGO_PKG_VERSION"));
    assert_eq!(attr_string(&attrs, "git_commit"), env!("NYX_BUILD_GIT_COMMIT"));
    assert_eq!(attr_string(&attrs, "rustc_version"), env!("NYX_BUILD_RUSTC_VERSION"));
    assert_eq!(attrs["nyx_threads"].as_i64(), Some(resolved_nyx_threads() as i64));
    assert_eq!(attrs["nyx_factor_threads"].as_i64(), Some(resolved_nyx_factor_threads() as i64));
    assert!(!attrs.contains_key("git_dirty"), "git_dirty should not sit on deterministic root attrs");

    // Repeated `.nyx` archive creation of the recorded workflow must also be 0-byte-diff identical.
    let store_info = call(&mut e, "project_info", json!({}))["store"].clone();
    let wf_root = PathBuf::from(store_info["folder"].as_str().unwrap());
    let wf = WorkflowDir::open(&wf_root).unwrap();
    let opt = ExportOptions { level: ShareLevel::Results, only: None, strip_paths: true, exclude: &[] };
    let nyx_a = export_bytes(&wf, &opt).unwrap();
    std::thread::sleep(std::time::Duration::from_millis(25));
    let nyx_b = export_bytes(&wf, &opt).unwrap();
    assert_eq!(nyx_a, nyx_b, "repeated .nyx archive export must be 0-byte-diff identical");

    // Verify archived results.h5 inside .nyx also has no created_utc and matches the run's results.h5.
    let nyx_file = root.join("workflow.nyx");
    std::fs::write(&nyx_file, &nyx_a).unwrap();
    let ins = inspect(&nyx_file, &Limits::default());
    assert!(ins.intact && ins.problems.is_empty(), "{:?}", ins.problems);

    let mut za = ZipArchive::new(std::io::Cursor::new(&nyx_a)).unwrap();
    let h5_entry_name = (0..za.len())
        .map(|i| za.by_index(i).unwrap().name().unwrap().into_owned())
        .find(|n| n.ends_with("/results.h5"))
        .expect("results.h5 inside .nyx archive");
    let mut archived_h5 = Vec::new();
    za.by_name(&h5_entry_name).unwrap().read_to_end(&mut archived_h5).unwrap();
    let extracted_h5 = root.join("from_archive.h5");
    std::fs::write(&extracted_h5, &archived_h5).unwrap();
    assert_no_created_utc_anywhere(&extracted_h5);

    let inbox = root.join("inbox");
    std::fs::create_dir_all(&inbox).unwrap();
    let imported = import(&nyx_file, &inbox, &Limits::default()).unwrap();
    assert_eq!(imported.runs("e1").0.len(), 1);

    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn provenance_materials_json_round_trips_and_open_results_warns_on_local_library_drift() {
    let root = tmp("mat-drift");
    let h5_path = root.join("soi_wire.h5");

    let mut e = Engine::new();
    call(&mut e, "project_open", json!({ "path": example("15_soi_wire") }));
    call(&mut e, "solve_run", json!({}));
    call(&mut e, "export_results", json!({ "path": h5_path.to_str().unwrap() }));

    // 1. `/provenance/materials_json` round-trips cleanly via h5io and results_file_info.
    let stored = nyx_engine::h5io::read_results(&h5_path).unwrap();
    assert_eq!(stored.materials.len(), 2, "15_soi_wire uses silicon + silica: {:?}", stored.materials);
    for mat in ["silicon", "silica"] {
        let entry = stored.materials.get(mat).unwrap_or_else(|| panic!("missing {mat} in {:?}", stored.materials));
        assert!(entry.blake3.starts_with("b3:"), "{entry:?}");
        assert_eq!(entry.source, "library");
        assert!(!entry.citations.is_empty(), "expected citations for {mat}");
        assert!(entry.params_snapshot.is_object());
    }
    let summary = nyx_engine::h5io::read_results_summary(&h5_path).unwrap();
    assert_eq!(summary.materials, stored.materials);

    let info = call(&mut e, "results_file_info", json!({ "path": h5_path.to_str().unwrap() }));
    assert!(info["materials"]["silicon"]["blake3"].as_str().unwrap().starts_with("b3:"));
    assert_eq!(info["material_warnings"], json!([]), "no drift against matching library");

    // 2. Opening against matching library succeeds with zero material_drift warnings.
    let mut clean_open = Engine::new();
    let clean_reply = call(&mut clean_open, "open_results", json!({ "path": h5_path.to_str().unwrap() }));
    let clean_warns = clean_reply["results"][0]["diagnostics"]["warnings"].as_array().cloned().unwrap_or_default();
    assert!(
        !clean_warns.iter().any(|w| w["code"] == "material_drift"),
        "matching library must emit no material_drift warning: {clean_warns:?}"
    );

    // 3. Modify local `silicon` in a fresh Engine's base library -> `open_results` still opens cleanly and warns.
    let mut drifted_texts = nyx_engine::worker::material_texts(&e.library().unwrap());
    drifted_texts.insert(
        "silicon".into(),
        r#"format = 2
name = "silicon"
source = "Modified synthetic silicon for drift detection"
[optical]
model = "constant"
n = 3.55
"#
        .to_string(),
    );
    let drifted_lib = femcore::material::Library::default().with_texts(&drifted_texts).unwrap();

    let mut drifted_engine = Engine::new();
    drifted_engine.set_base_library(drifted_lib);

    let info_drifted = call(&mut drifted_engine, "results_file_info", json!({ "path": h5_path.to_str().unwrap() }));
    let info_warns = info_drifted["material_warnings"].as_array().unwrap();
    assert!(
        info_warns.iter().any(|w| w.as_str().is_some_and(|s| s.contains("silicon"))),
        "results_file_info must warn about drifted silicon: {info_drifted}"
    );

    let opened_drifted = call(&mut drifted_engine, "open_results", json!({ "path": h5_path.to_str().unwrap() }));
    let open_warns = opened_drifted["results"][0]["diagnostics"]["warnings"]
        .as_array()
        .expect("diagnostics.warnings surfaced on open_results");
    assert!(
        open_warns.iter().any(|w| {
            w["code"] == "material_drift"
                && w["severity"] == "warning"
                && w["message"].as_str().is_some_and(|m| m.contains("silicon"))
        }),
        "open_results must surface material_drift warning for silicon: {opened_drifted}"
    );

    let (res_after, stale_after) = drifted_engine.result("ems").unwrap();
    assert!(!stale_after);
    assert!(res_after.diagnostics.as_ref().unwrap().warnings.iter().any(|w| w.code == "material_drift"
        && w.severity == DiagnosticSeverity::Warning
        && w.message.contains("silicon")));
    assert!(res_after.log.iter().any(|l| l.contains("material drift") && l.contains("silicon")));

    let rinfo = call(&mut drifted_engine, "result_info", json!({ "physics": "ems" }));
    assert!(rinfo[0]["diagnostics"]["warnings"].as_array().unwrap().iter().any(|w| w["code"] == "material_drift"));

    // 4. Empty local library (missing materials) also opens cleanly with explicit drift warnings.
    let mut empty_engine = Engine::new();
    empty_engine.set_base_library(femcore::material::Library::default());
    let opened_empty = call(&mut empty_engine, "open_results", json!({ "path": h5_path.to_str().unwrap() }));
    let empty_warns = opened_empty["results"][0]["diagnostics"]["warnings"].as_array().unwrap();
    assert!(
        empty_warns
            .iter()
            .any(|w| w["code"] == "material_drift" && w["message"].as_str().is_some_and(|m| m.contains("missing"))),
        "missing local library materials must emit material_drift warning: {opened_empty}"
    );

    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn store_link_persists_run_record_v1_1_provenance_and_gci_convergence() {
    let root = tmp("store-v11");
    let ws = root.join("ws");
    std::fs::create_dir_all(&ws).unwrap();

    // Use `silicon` in `dfs` electrostatic so SolverDiagnostics contains `material_approximate` warning
    // AND supports 3-grid GCI `convergence_check`.
    let mut e = Engine::new();
    e.set_workspace(ws.clone(), true);
    call(&mut e, "project_new", json!({ "title": "p5-prov-store", "dim": 2, "h": 0.25 }));
    call(&mut e, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "body", "input": "box" }] }));
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "body", "material": "silicon" }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "ground", "params": {} }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "potential", "params": { "value": 1.0 } }),
    );
    call(&mut e, "solve_run", json!({}));

    let runs_val = call(&mut e, "store_runs", json!({}));
    assert_eq!(runs_val.as_array().unwrap().len(), 1);
    let rec: nyx_store::runs::RunRecord = serde_json::from_value(runs_val[0].clone()).unwrap();
    assert_eq!(rec.status, RunStatus::Ok);

    // BuildProvenance matches compile-time build.rs contract.
    let bp = rec.build.expect("RunRecord.build populated");
    assert_eq!(bp, build_provenance());
    assert_eq!(bp.nyx_version, env!("CARGO_PKG_VERSION"));
    assert_eq!(bp.git_commit, env!("NYX_BUILD_GIT_COMMIT"));
    assert_eq!(bp.git_dirty, env!("NYX_BUILD_GIT_DIRTY") == "true");
    assert_eq!(bp.rustc_version, env!("NYX_BUILD_RUSTC_VERSION"));
    assert_eq!(bp.target_triple, env!("NYX_BUILD_TARGET"));
    assert_eq!(bp.profile, env!("NYX_BUILD_PROFILE"));
    assert_eq!(bp.cargo_lock_blake3, env!("NYX_BUILD_CARGO_LOCK_BLAKE3"));
    assert!(bp.cargo_lock_blake3.starts_with("b3:"));

    // RuntimeProvenance keeps created_utc in RunRecord only (never in HDF5).
    let rt = rec.runtime.expect("RunRecord.runtime populated");
    assert_eq!(rt.os, std::env::consts::OS);
    assert_eq!(rt.cpu_arch, std::env::consts::ARCH);
    assert_eq!(rt.nyx_threads, resolved_nyx_threads());
    assert_eq!(rt.nyx_factor_threads, resolved_nyx_factor_threads());
    assert_eq!(rt.created_utc, rec.t);
    assert!(nyx_store::util::parse_iso(&rt.created_utc).is_some());

    // RunDiagnosticsSummary captures counts, worst severity, algebraic residual, and warning summary.
    let diag = rec.diagnostics.expect("RunRecord.diagnostics populated");
    assert_eq!(diag.results_with_diagnostics, 1);
    assert!(diag.all_converged);
    assert!(diag.warning_count >= 1);
    assert_eq!(diag.worst_severity, "warning");
    assert!(diag.max_algebraic_residual.is_finite());
    assert!(diag.warning_codes.iter().any(|c| c == "material_approximate"), "{diag:?}");
    assert!(!diag.warnings.is_empty());

    // material_hashes matches `/provenance/materials_json` in the run's `results.h5`.
    let mh = rec.material_hashes.expect("RunRecord.material_hashes populated");
    let si_hash = mh.get("silicon").expect("silicon hash in RunRecord");
    assert!(si_hash.starts_with("b3:"));

    let store_folder = PathBuf::from(call(&mut e, "project_info", json!({}))["store"]["folder"].as_str().unwrap());
    let run_h5_rel =
        rec.artifacts.iter().find(|a| a.path.ends_with("results.h5")).expect("results.h5 artifact").path.clone();
    let run_h5_path = store_folder.join("e1").join(&run_h5_rel);
    assert_no_created_utc_anywhere(&run_h5_path);
    let stored_run_h5 = nyx_engine::h5io::read_results(&run_h5_path).unwrap();
    assert_eq!(stored_run_h5.materials.get("silicon").map(|m| &m.blake3), Some(si_hash));

    // 3-grid convergence check patches optional GCI fields on RunRecord.convergence.
    call(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.6, "grids": 3, "tolerance": 0.25 }));
    let runs_after = call(&mut e, "store_runs", json!({}));
    let rec_after: nyx_store::runs::RunRecord = serde_json::from_value(runs_after[0].clone()).unwrap();
    let conv = rec_after.convergence.expect("convergence patch persisted");
    assert_eq!(conv.grids, Some(3));
    assert!(conv.r_eff.is_some());
    assert!(conv.gci_fine.is_some());
    assert!(conv.gci_p_obs.is_some());
    assert!(conv.gci_regime.is_some());
    assert_eq!(conv.tier.as_deref(), Some("gci_verified"));

    // Direct StoreLink::record also populates failed runs with build/runtime/material_hashes.
    let direct_ws = root.join("direct_ws");
    std::fs::create_dir_all(&direct_ws).unwrap();
    let mut link = StoreLink::create(&direct_ws, e.project(), "cli").unwrap();
    let err = nyx_physics::Error::Solver("synthetic failure".into());
    let failed_rec = link.record(e.project(), RunOutcome::Failed(&err), 0.05, &["stopped".into()], None).unwrap();
    assert_eq!(failed_rec.status, RunStatus::Failed);
    assert!(failed_rec.build.is_some());
    assert!(failed_rec.runtime.is_some());
    assert!(failed_rec.diagnostics.is_none());
    assert_eq!(failed_rec.material_hashes.as_ref().and_then(|m| m.get("silicon")), Some(si_hash));

    let _ = std::fs::remove_dir_all(&root);
}
