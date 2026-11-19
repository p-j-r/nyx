//! Export as script (productivity plan S7b, ADR [94]): `project_script` writes a script that rebuilds the project. The gate: for every example
//! project and for a project shaped like one built in the window, the script type-checks, runs on a fresh engine, and gives the same project;
//! the readable form (`commands`) is required for all of them, so the fallbacks cannot hide a bug; `params` gives the variant.

use std::path::PathBuf;

use nyx_engine::Engine;
use nyx_script::kernel::{Kernel, Step};
use nyx_script::{host, Limits, Options};
use serde_json::{json, Value};

fn opts() -> Options {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!("nyx-export-test-{}-{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    std::fs::create_dir_all(&d).unwrap();
    Options { root: d, worker: Some(PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))) }
}

fn project_value(e: &mut Engine) -> Value {
    serde_json::from_str(e.call("project_json", &json!({})).unwrap()["json"].as_str().unwrap()).unwrap()
}

/// The same project: the same description of the geometry (names, not ids; no origin) and the same everything else.
fn assert_same(a: &mut Engine, b: &mut Engine, what: &str) {
    assert_eq!(a.call("list_features", &json!({})).unwrap(), b.call("list_features", &json!({})).unwrap(), "{what}: the geometry");
    let rest = |e: &mut Engine| {
        let mut v = project_value(e);
        v.as_object_mut().unwrap().remove("geometry");
        v
    };
    assert_eq!(rest(a), rest(b), "{what}: the rest of the project");
}

/// Exports `a`, checks the script, runs it on a fresh engine and compares. Returns the form and the script.
fn roundtrip(a: &mut Engine, what: &str) -> (String, String) {
    let r = a.call("project_script", &json!({})).unwrap_or_else(|e| panic!("{what}: {e}"));
    let (script, form) = (r["script"].as_str().unwrap().to_string(), r["form"].as_str().unwrap().to_string());
    let mut b = Engine::new();
    let o = opts();
    host::check(&mut b, "exported.nys", &script, &Limits::default(), &o).unwrap_or_else(|e| panic!("{what}: the script fails its check: {e}\n{script}"));
    host::run(&mut b, "exported.nys", &script, &Limits::default(), &o).unwrap_or_else(|e| panic!("{what}: the script fails: {e}\n{script}"));
    assert_same(a, &mut b, what);
    (form, script)
}

fn examples() -> Vec<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples");
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir).unwrap().filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.to_string_lossy().ends_with(".nyx.json")).collect();
    v.sort();
    v
}

#[test]
fn every_example_exports_as_commands_and_rebuilds_the_same_project() {
    let list = examples();
    assert!(list.len() >= 5, "the example projects were not found");
    let mut not_readable = vec![];
    for p in &list {
        let mut a = Engine::new();
        a.call("project_open", &json!({ "path": p.to_string_lossy() })).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
        let (form, script) = roundtrip(&mut a, &p.display().to_string());
        if form != "commands" {
            not_readable.push(format!("{} ({form})", p.file_name().unwrap().to_string_lossy()));
        }
        assert!(!script.contains("nyx.json") && !script.contains(":\\\\"), "the script names no file: {}", p.display());
    }
    assert!(not_readable.is_empty(), "these examples were exported in a fallback form, not as commands: {not_readable:?}");
}

/// A model as the window leaves it: features of origin `Gui`, an id gap from a removed feature, a suppressed and a tagged feature, a material of
/// the project, two physics with the active one not first, a mesh with every option.
fn window_shaped(sweep: bool) -> Engine {
    let mut e = Engine::new();
    let c = |e: &mut Engine, cmd: &str, args: Value| {
        e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}: {er}"));
    };
    c(&mut e, "project_new", json!({ "title": "built in the window", "dim": 2 }));
    c(&mut e, "set_param", json!({ "name": "w", "value": 5e-7 }));
    c(&mut e, "set_param", json!({ "name": "t", "value": "2.2e-7 + 0" }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [-2e-6, -2e-6], "max": [2e-6, 0.0] }));
    c(&mut e, "add_rect", json!({ "name": "junk", "min": [0.0, 0.0], "max": [1e-7, 1e-7] }));
    c(&mut e, "add_rect", json!({ "name": "clad", "min": [-2e-6, 0.0], "max": [2e-6, 2e-6] }));
    c(&mut e, "add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
    c(&mut e, "add_rect", json!({ "name": "spare", "min": [0.0, 0.0], "max": [1e-7, 1e-7] }));
    c(&mut e, "remove_feature", json!({ "name": "junk" }));
    c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "BOX", "input": "box" }, { "name": "Clad", "input": "clad" }, { "name": "Si", "input": "core" }] }));
    c(&mut e, "suppress_feature", json!({ "name": "spare" }));
    c(&mut e, "tag_feature", json!({ "name": "core", "tag": "wire" }));
    c(&mut e, "add_material", json!({ "name": "doped", "toml": "extends = \"silica\"\nname = \"doped\"\n[optical]\nmodel = \"constant\"\nn = 1.50\n" }));
    c(&mut e, "add_physics", json!({ "physics": "dfs" }));
    c(&mut e, "add_physics", json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 193.4e12, "n_target": 2.9, "n_modes": 2 } }));
    for (reg, mat) in [("BOX", "silica"), ("Clad", "doped"), ("Si", "silicon")] {
        c(&mut e, "set_material", json!({ "physics": "ems", "region": reg, "material": mat }));
    }
    c(&mut e, "set_active", json!({ "physics": "ems" }));
    c(&mut e, "set_mesh", json!({ "h": 2e-7, "adapt_scale": 0.4, "physics_sizing": true, "points_per_wavelength": 10.0 }));
    if sweep {
        c(&mut e, "set_mesh", json!({ "sweep": [{ "thickness": 1e-6, "n": 4 }] }));
    }
    // what the window writes: every feature is of origin Gui
    let mut v = project_value(&mut e);
    for f in v["geometry"]["features"].as_array_mut().unwrap() {
        f["origin"] = json!("Gui");
    }
    e.call("project_load", &json!({ "json": v.to_string() })).unwrap();
    e
}

#[test]
fn a_model_built_in_the_window_exports_as_commands_and_rebuilds_the_same_project() {
    let mut a = window_shaped(true);
    let ids: Vec<u64> = project_value(&mut a)["geometry"]["features"].as_array().unwrap().iter().map(|f| f["id"].as_u64().unwrap()).collect();
    assert!(ids.windows(2).any(|w| w[1] > w[0] + 1), "the model has an id gap: {ids:?}");
    let (form, script) = roundtrip(&mut a, "window-shaped");
    assert_eq!(form, "commands", "{script}");
    for want in ["nyx.project_new(", "nyx.set_param(name=\"w\", value=params.get(\"w\", 5e-7))", "nyx.add_rect(name=\"core\"", "nyx.suppress_feature(", "nyx.tag_feature(", "nyx.add_material(", "nyx.set_problem(", "nyx.set_active(physics=\"ems\")", "nyx.set_mesh("] {
        assert!(script.contains(want), "{want} is not in the script:\n{script}");
    }
    // the same model exports to the same text every time
    assert_eq!(a.call("project_script", &json!({})).unwrap()["script"], json!(script));
}

#[test]
fn params_give_the_variant_with_nothing_edited() {
    let mut a = window_shaped(true);
    let script = a.call("project_script", &json!({})).unwrap()["script"].as_str().unwrap().to_string();
    let mut b = Engine::new();
    let mut k = Kernel::new(&b, Limits::default(), opts()).unwrap();
    let Step::Done(_) = k.run(&mut b, "exported.nys", &script, Some(&json!({ "w": 3e-7 })), None, true).unwrap() else { panic!() };
    let params = b.call("project_info", &json!({})).unwrap()["params"].clone();
    let w = params.as_array().unwrap().iter().find(|p| p["name"] == "w").unwrap()["value"].as_f64().unwrap();
    assert_eq!(w, 3e-7, "{params}");
    // and the model still builds: the geometry names resolve with the new width
    b.call("check_project", &json!({})).unwrap();
}

#[test]
fn the_exported_script_can_be_run_by_an_agent_and_solved() {
    // the point of the export: a model built by hand, handed on, then solved by the script that follows it
    let mut a = window_shaped(false);
    let script = a.call("project_script", &json!({})).unwrap()["script"].as_str().unwrap().to_string();
    let mut b = Engine::new();
    let mut k = Kernel::new(&b, Limits::default(), opts()).unwrap();
    let src = format!("{script}\nnyx.solve_run(physics='ems')\nnyx.result_info(physics='ems')[0]['solutions'][0]['eigenvalue'][0]");
    let Step::Done(o) = k.run(&mut b, "exported.nys", &src, None, None, true).unwrap() else { panic!() };
    let neff = o.value.as_f64().unwrap();
    assert!(neff > 1.4 && neff < 3.5, "an effective index: {neff}");
}
