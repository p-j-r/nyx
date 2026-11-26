// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Interface contracts of nyx-store (docs/architecture/adr/0014 and 0019): the files other programs and later versions of nyx read.
//!
//! * `store_workflow` / `store_run`: exemplar `workflow.json` and `runs.jsonl` line (docs/architecture/contracts/store_*.v1.json).
//! * `store_format_gate`: a manifest from a newer nyx is refused, not misread.
//!
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-store --test contracts` rewrites the exemplars (only with a version bump, ADR 19).
use nyx_store::model::{Card, Digest, Experiment, Status, Workflow, FORMAT};
use nyx_store::runs::{Artifact, BuildProvenance, Convergence, RunDiagnosticsSummary, RunRecord, RunStatus, RuntimeProvenance};
use nyx_store::WorkflowDir;
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

fn check_exemplar(name: &str, got: &Value, pretty: bool) {
    let path = contract_file(name);
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let text = if pretty { serde_json::to_string_pretty(got).unwrap() } else { serde_json::to_string(got).unwrap() };
        std::fs::write(&path, format!("{text}\n")).unwrap();
        return;
    }
    let want_text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let want: Value = serde_json::from_str(&want_text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(&want, got, "{name} no longer matches the contract exemplar; this is a change of a saved format (ADR 19, ADR 26)");
}

/// A project as the app saves it today: `Project { version, geometry: GeometrySpec { list, h, extrude, .. }, source, active, specs, colormap }`.
fn current_project() -> Value {
    json!({
        "version": 1,
        "geometry": { "list": { "dim": 2, "params": [], "features": [] }, "h": 5e-8, "extrude": null },
        "source": "Generated",
        "active": "Emw",
        "specs": [{
            "physics": "Emw", "order": 2,
            "materials": { "default": { "eps_r": 1.0 }, "domains": { "1": { "eps_r": 12.08 } } },
            "study": { "Modes": { "freq": 1.9341e14, "n_target": 2.4, "evanescent": false, "n_modes": 4 } }
        }],
        "colormap": "Viridis"
    })
}

fn tmp(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-store-contract-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn store_workflow_manifest_exemplar_is_stable() {
    let mut wf = Workflow::new("Exemplar", "gui");
    wf.id = "wf-example".into();
    wf.created = "2026-11-03T00:00:00Z".into();
    wf.touched = "2026-11-03T00:00:00Z".into();
    wf.tags = vec!["soi".into()];
    wf.version = 1;
    wf.experiments.push(Experiment {
        id: "e1".into(),
        title: "TE modes".into(),
        parent: None,
        links: vec![],
        spec_hash: "b3:0000".into(),
        card: Card::default(),
        runs: Digest::default(),
        model: None,
    });
    assert_eq!(wf.nyx, FORMAT);
    assert_eq!(wf.status, Status::Active);
    check_exemplar("store_workflow.v1.json", &serde_json::to_value(&wf).unwrap(), true);
    let back: Workflow = serde_json::from_str(&std::fs::read_to_string(contract_file("store_workflow.v1.json")).unwrap()).unwrap();
    assert_eq!(back, wf);
}

#[test]
fn store_run_record_exemplar_is_stable_and_old_lines_still_read() {
    let r = RunRecord {
        id: "r0001_abcd".into(),
        t: "2026-11-03T00:00:00Z".into(),
        sec: 0.25,
        status: RunStatus::Ok,
        spec: "b3:0000".into(),
        femcore: "0.1.0".into(),
        version: 2,
        params: [("w".to_string(), 5e-7)].into(),
        metrics: [("neff_re".to_string(), 2.44)].into(),
        tags: vec!["baseline".into()],
        star: true,
        note: "first".into(),
        artifacts: vec![Artifact { kind: "fields".into(), path: "r0001_abcd/fields.h5".into(), bytes: 1024 }],
        error: None,
        convergence: Some(Convergence {
            factor: 0.6,
            tolerance: 1e-3,
            observable: "2 quantities (eigenvalues, values, numeric scalars), solutions matched by index, 3-grid GCI".into(),
            passed: true,
            worst_relative: 2.5e-4,
            grids: Some(3),
            r_eff: Some(1.67),
            gci_fine: Some(1.2e-4),
            gci_p_obs: Some(2.0),
            gci_regime: Some("monotonic".into()),
            tier: Some("gci_verified".into()),
        }),
        reuse: None,
        build: Some(BuildProvenance {
            nyx_version: "0.1.0".into(),
            git_commit: "000000000000".into(),
            git_dirty: false,
            rustc_version: "rustc 1.99.1".into(),
            target_triple: "x86_64-unknown-linux-gnu".into(),
            profile: "release".into(),
            cargo_lock_blake3: "b3:0000".into(),
        }),
        runtime: Some(RuntimeProvenance {
            os: "linux".into(),
            cpu_arch: "x86_64".into(),
            nyx_threads: 4,
            nyx_factor_threads: 4,
            created_utc: "2026-11-03T00:00:00Z".into(),
        }),
        diagnostics: Some(RunDiagnosticsSummary {
            results_with_diagnostics: 1,
            all_converged: true,
            info_count: 0,
            warning_count: 0,
            error_count: 0,
            worst_severity: "none".into(),
            max_algebraic_residual: 1.0e-14,
            min_rcond_estimate: Some(1.0e-4),
            cholesky_fallbacks: 0,
            warning_codes: vec![],
            warnings: vec![],
        }),
        material_hashes: Some([("silicon".to_string(), "b3:1111".to_string())].into()),
    };
    check_exemplar("store_run.v1.json", &serde_json::to_value(&r).unwrap(), false);
    // a line written before `version`, `tags`, `star`, `note` and `artifacts` existed still reads (they are `#[serde(default)]`)
    let old = r#"{"id":"r1","t":"2026-10-31T00:00:00Z","sec":1.0,"status":"ok","spec":"b3:1","femcore":"0.1.0","params":{},"metrics":{}}"#;
    let rec: RunRecord = serde_json::from_str(old).unwrap();
    assert_eq!((rec.version, rec.star, rec.tags.len(), rec.note.as_str(), rec.artifacts.len()), (0, false, 0, "", 0));
    assert!(rec.build.is_none() && rec.runtime.is_none() && rec.diagnostics.is_none() && rec.material_hashes.is_none());
    // a store_run 1.0 record still deserializes cleanly under schema 1.1
    let v1_0 = r#"{"id":"r0001_abcd","t":"2026-11-03T00:00:00Z","sec":0.25,"status":"ok","spec":"b3:0000","femcore":"0.1.0","version":2,"params":{"w":5e-7},"metrics":{"neff_re":2.44},"tags":["baseline"],"star":true,"note":"first","artifacts":[{"kind":"fields","path":"r0001_abcd/fields.h5","bytes":1024}]}"#;
    let rec_v1_0: RunRecord = serde_json::from_str(v1_0).unwrap();
    assert_eq!(rec_v1_0.id, "r0001_abcd");
    assert!(rec_v1_0.build.is_none() && rec_v1_0.runtime.is_none() && rec_v1_0.diagnostics.is_none() && rec_v1_0.material_hashes.is_none());
    assert_eq!(serde_json::to_value(RunStatus::Cancelled).unwrap(), json!("cancelled"));
}

#[test]
fn store_format_gate_refuses_a_newer_manifest() {
    let ws = tmp("gate");
    let d = WorkflowDir::create(&ws, "gate", "gui", &current_project().to_string()).unwrap();
    let text = std::fs::read_to_string(d.manifest_path()).unwrap();
    assert!(text.contains(&format!("\"nyx\": {FORMAT}")), "the manifest carries its format under the key `nyx`");
    std::fs::write(d.manifest_path(), text.replace(&format!("\"nyx\": {FORMAT}"), "\"nyx\": 99")).unwrap();
    assert!(WorkflowDir::open(&d.root).unwrap_err().contains("newer nyx"));
    // the files a workflow is made of (ADR 14)
    assert!(d.spec_path("e1").ends_with("e1/experiment.nyx.json") || d.spec_path("e1").to_string_lossy().replace('\\', "/").ends_with("e1/experiment.nyx.json"));
    let _ = std::fs::remove_dir_all(&ws);
}

