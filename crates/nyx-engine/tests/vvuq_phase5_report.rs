// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 5 VVUQ verification tests for `vvuq_export` and `vvuq_report.v1.json`:
//! - Strict JSON Schema (`vvuq_report_schema()`) validation & non-vacuity checks;
//! - End-to-end `vvuq_export` on a solved & convergence-checked project with active material
//!   warnings (`material_approximate`, `temperature_clamped`) and HDF5 `results.h5` export/linking;
//! - Independent BLAKE3 digest verification of canonical project JSON, solver mesh, resolved
//!   materials, and `results.h5`;
//! - `path_guard` confinement refusal outside `path_root` (`..`, absolute, and symlink/junction escapes);
//! - Contract snapshot check (`NYX_WRITE_CONTRACT=1` writes `docs/architecture/contracts/vvuq_report.v1.{json,schema.json}`).

use nyx_engine::vvuq_report::{
    canonical_mesh_blake3, canonical_project_blake3, vvuq_report_exemplar, vvuq_report_schema, VvuqReportV1,
    VVUQ_REPORT_FORMAT,
};
use nyx_engine::Engine;
use nyx_physics::{canonical_json, DiagnosticSeverity};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

fn validate(root: &Value, schema: &Value, inst: &Value, at: &str, errs: &mut Vec<String>) {
    const KNOWN: [&str; 15] = [
        "$schema",
        "$id",
        "title",
        "description",
        "$defs",
        "$ref",
        "oneOf",
        "const",
        "enum",
        "type",
        "required",
        "properties",
        "additionalProperties",
        "items",
        "minimum",
    ];
    let obj = schema.as_object().expect("schema node must be an object");
    for k in obj.keys() {
        assert!(KNOWN.contains(&k.as_str()), "unsupported JSON Schema keyword `{k}`");
    }
    if let Some(r) = obj.get("$ref") {
        let name = r.as_str().unwrap().strip_prefix("#/$defs/").expect("only #/$defs refs");
        return validate(root, &root["$defs"][name], inst, at, errs);
    }
    if let Some(alts) = obj.get("oneOf") {
        let n = alts
            .as_array()
            .unwrap()
            .iter()
            .filter(|a| {
                let mut e = Vec::new();
                validate(root, a, inst, at, &mut e);
                e.is_empty()
            })
            .count();
        if n != 1 {
            errs.push(format!("{at}: matched {n} of oneOf alternatives (expected 1)"));
        }
    }
    if let Some(c) = obj.get("const") {
        if c != inst {
            errs.push(format!("{at}: expected const {c}, got {inst}"));
        }
    }
    if let Some(e) = obj.get("enum") {
        if !e.as_array().unwrap().contains(inst) {
            errs.push(format!("{at}: {inst} not in enum {e}"));
        }
    }
    if let Some(t) = obj.get("type") {
        let types: Vec<&str> = match t {
            Value::String(s) => vec![s.as_str()],
            Value::Array(a) => a.iter().map(|v| v.as_str().unwrap()).collect(),
            _ => panic!("bad schema type keyword"),
        };
        let ok = types.iter().any(|t| match *t {
            "object" => inst.is_object(),
            "array" => inst.is_array(),
            "string" => inst.is_string(),
            "boolean" => inst.is_boolean(),
            "null" => inst.is_null(),
            "integer" => inst.is_i64() || inst.is_u64(),
            "number" => inst.is_number(),
            other => panic!("unknown schema type `{other}`"),
        });
        if !ok {
            errs.push(format!("{at}: expected {types:?}, got {inst}"));
        }
    }
    if let Some(min) = obj.get("minimum").and_then(Value::as_f64) {
        if let Some(v) = inst.as_f64() {
            if v < min {
                errs.push(format!("{at}: {v} < minimum {min}"));
            }
        }
    }
    if let Some(m) = inst.as_object() {
        for k in obj.get("required").and_then(Value::as_array).into_iter().flatten() {
            if !m.contains_key(k.as_str().unwrap()) {
                errs.push(format!("{at}: missing required key {k}"));
            }
        }
        let props = obj.get("properties").and_then(Value::as_object);
        let addl = obj.get("additionalProperties");
        for (k, v) in m {
            match props.and_then(|p| p.get(k)) {
                Some(s) => validate(root, s, v, &format!("{at}/{k}"), errs),
                None => match addl {
                    Some(Value::Bool(false)) => errs.push(format!("{at}: unexpected key `{k}`")),
                    Some(sub) if sub.is_object() => validate(root, sub, v, &format!("{at}/{k}"), errs),
                    _ => {}
                },
            }
        }
    }
    if let (Some(items), Some(arr)) = (obj.get("items"), inst.as_array()) {
        for (i, v) in arr.iter().enumerate() {
            validate(root, items, v, &format!("{at}/{i}"), errs);
        }
    }
}

fn assert_valid_report(val: &Value, label: &str) {
    let schema = vvuq_report_schema();
    let mut errs = Vec::new();
    validate(&schema, &schema, val, "", &mut errs);
    assert!(errs.is_empty(), "{label} failed vvuq_report_schema:\n{}", errs.join("\n"));
}

fn temp_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-vvuq-p5-report-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn setup_thermal_copper_project() -> Engine {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "P5 Thermal Report", "dim": 2, "h": 0.25 }));
    call(&mut e, "add_rect", json!({ "name": "bar", "min": [0.0, 0.0], "max": [1.0, 0.5] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "bar", "input": "bar" }] }));
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "bar", "material": "copper_4k" }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": 10.0 } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "temperature", "params": { "value": 420.0 } }),
    );
    e
}

#[test]
fn vvuq_report_v1_schema_exemplar_and_non_vacuity() {
    let schema = vvuq_report_schema();
    let exemplar = serde_json::to_value(vvuq_report_exemplar()).unwrap();
    assert_valid_report(&exemplar, "exemplar");

    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        let s_path = contract_file("vvuq_report.v1.schema.json");
        let e_path = contract_file("vvuq_report.v1.json");
        std::fs::create_dir_all(s_path.parent().unwrap()).unwrap();
        std::fs::write(&s_path, serde_json::to_string_pretty(&canonical_json(&schema)).unwrap() + "\n").unwrap();
        std::fs::write(&e_path, serde_json::to_string_pretty(&canonical_json(&exemplar)).unwrap() + "\n").unwrap();
    } else {
        let s_path = contract_file("vvuq_report.v1.schema.json");
        let e_path = contract_file("vvuq_report.v1.json");
        if s_path.is_file() && e_path.is_file() {
            let want_schema: Value = serde_json::from_str(&std::fs::read_to_string(&s_path).unwrap()).unwrap();
            let want_ex: Value = serde_json::from_str(&std::fs::read_to_string(&e_path).unwrap()).unwrap();
            assert_eq!(canonical_json(&want_schema), canonical_json(&schema), "vvuq_report.v1.schema.json drifted");
            assert_eq!(canonical_json(&want_ex), canonical_json(&exemplar), "vvuq_report.v1.json drifted");
        }
    }

    // Non-vacuity: corrupted reports must be rejected.
    let mut bad_format = exemplar.clone();
    bad_format["format"] = json!("vvuq_report.v99");
    let mut bad_missing = exemplar.clone();
    bad_missing.as_object_mut().unwrap().remove("project_blake3");
    let mut bad_signed = exemplar.clone();
    bad_signed["signature"] = json!("must-stay-unsigned");
    let mut bad_tier = exemplar.clone();
    bad_tier["validated"]["tier"] = json!("imaginary_tier");
    let mut bad_sev = exemplar.clone();
    bad_sev["material_warnings"][0]["severity"] = json!("critical");
    let mut bad_dim = exemplar.clone();
    bad_dim["mesh_check"]["dim"] = json!(4);

    for (why, inst) in [
        ("wrong format", bad_format),
        ("missing project_blake3", bad_missing),
        ("unexpected signature key", bad_signed),
        ("bad validation tier", bad_tier),
        ("bad material warning severity", bad_sev),
        ("invalid mesh dim", bad_dim),
    ] {
        let mut errs = Vec::new();
        validate(&schema, &schema, &inst, "", &mut errs);
        assert!(!errs.is_empty(), "validator erroneously accepted {why}");
    }
}

#[test]
fn vvuq_export_end_to_end_digests_diagnostics_warnings_and_h5() {
    let dir = temp_dir("e2e");

    // 1. Unmeshed, unsolved project with inline material still exports valid vvuq_report.v1.
    let mut e_unmeshed = Engine::new();
    call(&mut e_unmeshed, "project_new", json!({ "title": "Unmeshed", "dim": 2, "h": 0.25 }));
    call(&mut e_unmeshed, "add_rect", json!({ "name": "dom", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    call(&mut e_unmeshed, "add_overlay", json!({ "regions": [{ "name": "dom", "input": "dom" }] }));
    call(&mut e_unmeshed, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary" }));
    call(&mut e_unmeshed, "set_material", json!({ "physics": "dfs", "region": "dom", "material": { "k": 15.0 } }));

    let unmeshed_path = dir.join("unmeshed.v1.json");
    let unmeshed_reply = call(&mut e_unmeshed, "vvuq_export", json!({ "path": unmeshed_path.to_str().unwrap() }));
    assert_eq!(unmeshed_reply["format"], VVUQ_REPORT_FORMAT);
    let unmeshed_disk: Value = serde_json::from_str(&std::fs::read_to_string(&unmeshed_path).unwrap()).unwrap();
    assert_valid_report(&unmeshed_disk, "unmeshed");
    assert!(unmeshed_disk["mesh_blake3"].is_null());
    assert!(unmeshed_disk["mesh_check"].is_null());
    assert!(unmeshed_disk["results_h5_blake3"].is_null());
    assert!(unmeshed_disk["materials_blake3"]["inline:dfs:dom"].is_string());
    assert_eq!(unmeshed_disk["physics"]["dfs"]["solved"], json!(false));
    assert_eq!(unmeshed_disk["validated"]["state"], json!("exploratory"));

    // 2. Full solved & 3-grid GCI convergence-checked thermal copper_4k project with warnings + HDF5.
    let mut e = setup_thermal_copper_project();
    call(&mut e, "solve_run", json!({ "physics": "dfs" }));
    let cc =
        call(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.6, "grids": 3, "tolerance": 0.1 }));
    assert_eq!(cc["converged"], json!(true), "{cc}");

    // Export report without HDF5 first -> results_h5_blake3 is null.
    let rep_no_h5 = dir.join("nested/report_no_h5.v1.json");
    call(&mut e, "vvuq_export", json!({ "path": rep_no_h5.to_str().unwrap() }));
    let parsed_no_h5: VvuqReportV1 = serde_json::from_str(&std::fs::read_to_string(&rep_no_h5).unwrap()).unwrap();
    assert!(parsed_no_h5.results_h5_blake3.is_none());

    // Now export with `export_h5: true` -> companion results.h5 created and BLAKE3-linked!
    let rep_with_h5 = dir.join("companion/vvuq_report.v1.json");
    let reply_h5 = call(&mut e, "vvuq_export", json!({ "path": rep_with_h5.to_str().unwrap(), "export_h5": true }));
    assert_valid_report(&reply_h5["report"], "reply report");
    let disk_raw = std::fs::read_to_string(&rep_with_h5).unwrap();
    let disk_val: Value = serde_json::from_str(&disk_raw).unwrap();
    assert_valid_report(&disk_val, "disk report");
    let rep: VvuqReportV1 = serde_json::from_value(disk_val.clone()).unwrap();

    // Check format, timestamp, and independent BLAKE3 digests.
    assert_eq!(rep.format, VVUQ_REPORT_FORMAT);
    assert!(nyx_store::util::parse_iso(&rep.generated_utc).is_some(), "{}", rep.generated_utc);
    assert_eq!(rep.project_blake3, canonical_project_blake3(e.project()));

    let solver_mesh = e.model().unwrap().solver_mesh();
    assert_eq!(rep.mesh_blake3.as_deref(), Some(canonical_mesh_blake3(&solver_mesh).as_str()));

    let lib = e.library().unwrap();
    let cu_toml = lib.get("copper_4k").expect("copper_4k in library").to_toml();
    assert_eq!(
        rep.materials_blake3.get("copper_4k").map(String::as_str),
        Some(nyx_store::util::hash_bytes(cu_toml.as_bytes()).as_str())
    );

    if cfg!(feature = "hdf5") {
        let companion_h5 = PathBuf::from(reply_h5["results_h5_path"].as_str().expect("results_h5_path"));
        assert!(companion_h5.is_file(), "companion results.h5 missing at {}", companion_h5.display());
        let expected_h5_b3 = nyx_store::util::hash_file(&companion_h5).unwrap();
        assert_eq!(rep.results_h5_blake3.as_deref(), Some(expected_h5_b3.as_str()));
    }

    // Explicit `export_results` path updates linked `results_h5_blake3`.
    if cfg!(feature = "hdf5") {
        let explicit_h5 = dir.join("explicit_out.h5");
        call(&mut e, "export_results", json!({ "path": explicit_h5.to_str().unwrap() }));
        let rep_exp = dir.join("explicit_report.v1.json");
        call(&mut e, "vvuq_export", json!({ "path": rep_exp.to_str().unwrap() }));
        let parsed_exp: VvuqReportV1 = serde_json::from_str(&std::fs::read_to_string(&rep_exp).unwrap()).unwrap();
        assert_eq!(
            parsed_exp.results_h5_blake3.as_deref(),
            Some(nyx_store::util::hash_file(&explicit_h5).unwrap().as_str())
        );
    }

    // Build metadata checks (fallbacks when NYX_BUILD_* not set, or exact values when set).
    assert_eq!(rep.build.nyx_version, env!("CARGO_PKG_VERSION"));
    assert_eq!(rep.build.os, std::env::consts::OS);
    assert_eq!(rep.build.arch, std::env::consts::ARCH);
    assert!(!rep.build.target_triple.is_empty());
    assert!(!rep.build.profile.is_empty());
    assert!(rep.build.nyx_threads >= 1 && rep.build.nyx_factor_threads >= 1);
    if let Some(gc) = option_env!("NYX_BUILD_GIT_COMMIT").filter(|s| !s.is_empty()) {
        assert_eq!(rep.build.git_commit.as_deref(), Some(gc));
    }

    // SimplexMesh::check() summary matches solved mesh.
    let mc = rep.mesh_check.as_ref().expect("mesh_check");
    assert!(mc.ok, "{:?}", mc.error);
    assert!(mc.error.is_none());
    assert_eq!(mc.dim, 2);
    assert_eq!(mc.n_points, solver_mesh.n_points());
    assert_eq!(mc.n_cells, solver_mesh.n_cells());
    assert_eq!(mc.n_facets, solver_mesh.n_facets());
    assert_eq!(mc.n_edges, solver_mesh.n_edges());

    // Material warnings include both approximate provenance and temperature clamping.
    assert!(
        rep.material_warnings.iter().any(|w| w.physics == "dfs"
            && w.code == "material_approximate"
            && w.severity == DiagnosticSeverity::Warning),
        "{:?}",
        rep.material_warnings
    );
    assert!(
        rep.material_warnings
            .iter()
            .any(|w| w.physics == "dfs" && w.code == "temperature_clamped" && w.message.contains("copper_4k")),
        "{:?}",
        rep.material_warnings
    );

    // Per-physics diagnostics, error budget, and GCI validation tier.
    let dfs = rep.physics.get("dfs").expect("dfs entry");
    assert!(dfs.solved && !dfs.stale);
    assert!(dfs.dofs.is_some_and(|d| d > 0));
    assert!(dfs.diagnostics.as_ref().is_some_and(|d| d.converged && d.max_algebraic_residual < 1e-10));
    let eb = dfs.error_budget.as_ref().expect("error_budget");
    let sp = eb.spatial_discretization.as_ref().expect("spatial_discretization");
    assert_eq!(sp.method, "roache_gci_3grid");
    assert_eq!(dfs.validated.state, "checked");
    assert_eq!(dfs.validated.tier, "gci_verified");
    assert_eq!(dfs.validated.check.as_ref().map(|c| c.grids), Some(3));
    assert_eq!(rep.validated.state, "checked");
    assert_eq!(rep.validated.tier, "gci_verified");

    // 3. Store-linked run (`store_save` + `solve_run`) links runs/<id>/results.h5 automatically.
    if cfg!(feature = "hdf5") {
        let ws = dir.join("workspace");
        let mut es = setup_thermal_copper_project();
        call(&mut es, "store_save", json!({ "workspace": ws.to_str().unwrap() }));
        let solve_ans = call(&mut es, "solve_run", json!({ "physics": "dfs" }));
        let run_id = solve_ans["run"]["id"].as_str().unwrap();
        let store_rep = dir.join("store_report.v1.json");
        let store_reply = call(&mut es, "vvuq_export", json!({ "path": store_rep.to_str().unwrap() }));
        let linked_h5 = Path::new(store_reply["results_h5_path"].as_str().unwrap());
        assert!(linked_h5.to_string_lossy().contains(run_id), "{}", linked_h5.display());
        let expected_b3 = nyx_store::util::hash_file(linked_h5).unwrap();
        assert_eq!(store_reply["report"]["results_h5_blake3"], json!(expected_b3));
    }

    // 4. Editing the project invalidates the linked HDF5 hash and marks stored results stale.
    call(&mut e, "set_mesh", json!({ "h": 0.20 }));
    let stale_rep_path = dir.join("stale_report.v1.json");
    call(&mut e, "vvuq_export", json!({ "path": stale_rep_path.to_str().unwrap() }));
    let stale_rep: VvuqReportV1 = serde_json::from_str(&std::fs::read_to_string(&stale_rep_path).unwrap()).unwrap();
    assert!(stale_rep.results_h5_blake3.is_none());
    assert!(stale_rep.mesh_blake3.is_none());
    assert!(stale_rep.physics["dfs"].stale);
    assert_eq!(stale_rep.validated.state, "exploratory");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn vvuq_export_path_guard_enforces_root_confinement() {
    let base = temp_dir("guard");
    let root = base.join("sandbox");
    let outside = base.join("outside");
    std::fs::create_dir_all(root.join("sub")).unwrap();
    std::fs::create_dir_all(&outside).unwrap();

    let mut e = setup_thermal_copper_project();
    call(&mut e, "solve_run", json!({ "physics": "dfs" }));
    e.set_path_root(Some(&root)).unwrap();

    // Parent traversal and absolute paths outside root are refused without writing files.
    let err_up = e.call("vvuq_export", &json!({ "path": "../escape.v1.json" })).unwrap_err().to_string();
    assert!(err_up.contains("vvuq_export") && err_up.contains("path"), "{err_up}");
    assert!(!base.join("escape.v1.json").exists());

    let abs_out = outside.join("abs.v1.json");
    assert!(e.call("vvuq_export", &json!({ "path": abs_out.to_str().unwrap() })).is_err());
    assert!(!abs_out.exists());

    // Direct Engine::vvuq_export call is also guarded.
    assert!(e.vvuq_export(abs_out.to_str().unwrap(), false).is_err());
    assert!(!abs_out.exists());

    // Symlink / junction escape is refused.
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(&outside, root.join("link_out")).unwrap();
        assert!(e.call("vvuq_export", &json!({ "path": "link_out/evil.v1.json" })).is_err());
        assert!(!outside.join("evil.v1.json").exists());
    }
    #[cfg(windows)]
    {
        let status = std::process::Command::new("cmd")
            .args(["/C", "mklink", "/J"])
            .arg(root.join("link_out"))
            .arg(&outside)
            .output()
            .unwrap();
        assert!(status.status.success());
        assert!(e.call("vvuq_export", &json!({ "path": "link_out/evil.v1.json" })).is_err());
        assert!(!outside.join("evil.v1.json").exists());
    }

    // Relative path inside root succeeds and writes both report and companion results.h5 inside root.
    let ok_reply = call(
        &mut e,
        "vvuq_export",
        json!({ "path": "reports/vvuq_report.v1.json", "export_h5": cfg!(feature = "hdf5") }),
    );
    let inside_file = root.join("reports/vvuq_report.v1.json");
    assert!(inside_file.is_file());
    assert_valid_report(&ok_reply["report"], "confined inside root");
    if cfg!(feature = "hdf5") {
        assert!(root.join("reports/results.h5").is_file());
    }

    let _ = std::fs::remove_dir_all(&base);
}
