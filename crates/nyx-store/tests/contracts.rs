//! Interface contracts of nyx-store (docs/architecture/adr/0014 and 0019): the files other programs and later versions of nyx read.
//!
//! * `store_workflow` / `store_run`: exemplar `workflow.json` and `runs.jsonl` line (docs/architecture/contracts/store_*.v1.json).
//! * `store_format_gate`: a manifest from a newer nyx is refused, not misread.
//! * `store_metrics`: the names under which a `SolveResult` appears in the runs table.
//! * `store_card`: the concept card is read from the *current* saved project. One test is `#[ignore]`d until ROUGH_SPOTS R25 is fixed.
//!
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-store --test contracts` rewrites the exemplars (only with a version bump, ADR 19).
use femcore::problem::{Field, Scalar, Solution, SolveResult};
use nyx_store::model::{Card, Digest, Experiment, Status, Workflow, FORMAT};
use nyx_store::runs::{Artifact, RunRecord, RunStatus};
use nyx_store::WorkflowDir;
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

fn check_exemplar(name: &str, got: &str) {
    let path = contract_file(name);
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, format!("{got}\n")).unwrap();
        return;
    }
    let want = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(want.replace("\r\n", "\n").trim_end(), got.trim_end(), "{name} no longer matches the contract exemplar; this is a change of a saved format (ADR 19)");
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
    check_exemplar("store_workflow.v1.json", &serde_json::to_string_pretty(&wf).unwrap());
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
    };
    check_exemplar("store_run.v1.json", &serde_json::to_string(&r).unwrap());
    // a line written before `version`, `tags`, `star`, `note` and `artifacts` existed still reads (they are `#[serde(default)]`)
    let old = r#"{"id":"r1","t":"2026-10-31T00:00:00Z","sec":1.0,"status":"ok","spec":"b3:1","femcore":"0.1.0","params":{},"metrics":{}}"#;
    let rec: RunRecord = serde_json::from_str(old).unwrap();
    assert_eq!((rec.version, rec.star, rec.tags.len(), rec.note.as_str(), rec.artifacts.len()), (0, false, 0, "", 0));
    assert_eq!(serde_json::to_value(RunStatus::Cancelled).unwrap(), json!("cancelled"));
}

#[test]
fn store_format_gate_refuses_a_newer_manifest() {
    let ws = tmp("gate");
    let xml = femcore::xmlio::value_to_xml("nyx_project", &current_project()).unwrap();
    let d = WorkflowDir::create(&ws, "gate", "gui", &xml).unwrap();
    let text = std::fs::read_to_string(d.manifest_path()).unwrap();
    assert!(text.contains(&format!("\"nyx\": {FORMAT}")), "the manifest carries its format under the key `nyx`");
    std::fs::write(d.manifest_path(), text.replace(&format!("\"nyx\": {FORMAT}"), "\"nyx\": 99")).unwrap();
    assert!(WorkflowDir::open(&d.root).unwrap_err().contains("newer nyx"));
    // the files a workflow is made of (ADR 14)
    assert!(d.spec_path("e1").ends_with("e1/experiment.nyx.xml") || d.spec_path("e1").to_string_lossy().replace('\\', "/").ends_with("e1/experiment.nyx.xml"));
    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn store_metrics_are_named_from_the_solve_result() {
    let mut r = SolveResult { dofs: 1234, ..Default::default() };
    r.solutions.push(Solution { neff: Some((2.4, 1e-9)), value: Some(5.0), fields: vec![Field { name: "Ex".into(), unit: String::new(), values: vec![0.0] }], ..Default::default() });
    r.scalars.push(("Gain G (1/W/m)".into(), Scalar::Num(1.5)));
    let m = nyx_store::runs::metrics_of(&r);
    let keys: Vec<&str> = m.keys().map(|s| s.as_str()).collect();
    assert_eq!(keys, vec!["dofs", "gain_g_1_w_m", "n_solutions", "neff_im", "neff_re", "value_1"]);
    assert_eq!(m["dofs"], 1234.0);
}

#[test]
fn store_card_still_gets_physics_study_and_wavelength_from_the_current_project() {
    let card = nyx_store::card::from_value(&current_project());
    assert_eq!(card.dim, 2);
    let w = card.wavelength_m.expect("wavelength is C0 / freq of the Modes study");
    assert!((w - 299_792_458.0 / 1.9341e14).abs() < 1e-12, "{w}");
    assert!(card.summary.contains("1550 nm") || card.summary.contains("1551 nm"), "{}", card.summary);
}

#[test]
#[ignore = "ROUGH_SPOTS R25: card::geometry reads the removed layout geometry.geometry.{Rectangle,Block,Shapes,Solids}; a current project has geometry.list"]
fn store_card_reads_current_project() {
    let card = nyx_store::card::from_value(&current_project());
    assert_ne!(card.geometry.family, "unknown", "the geometry of a current project (a geom FeatureList) must be recognised");
    assert_ne!(card.geometry.summary, "geometry not recognised");
}
