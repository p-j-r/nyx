// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Materials, analysis and a sized mesh through engine commands alone: a silicon wire in oxide is built from geometry commands, its
//! materials are named library files (silicon and silica, Sellmeier or tabulated models evaluated at the study wavelength), it is solved
//! with the real `ems` physics and its modes are analysed. Nothing here names a physics type or a solver.

use nyx_engine::Engine;
use serde_json::{json, Value};

fn wire(e: &mut Engine) {
    for (cmd, args) in [
        ("project_new", json!({ "title": "wire", "h": 1.0e-7 })),
        ("add_rect", json!({ "name": "box", "min": [-2.5e-6, -2e-6], "max": [2.5e-6, 0.0] })),
        ("add_rect", json!({ "name": "clad", "min": [-2.5e-6, 0.0], "max": [2.5e-6, 1.5e-6] })),
        ("add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 3e-7] })),
        (
            "add_overlay",
            json!({ "name": "regions", "regions": [
            { "name": "BOX", "input": "box" }, { "name": "clad", "input": "clad" }, { "name": "Si", "input": "core" } ] }),
        ),
        (
            "add_physics",
            json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 3 }, "order": 2 }),
        ),
        ("set_material", json!({ "physics": "ems", "region": "BOX", "material": "silica" })),
        ("set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" })),
        ("set_material", json!({ "physics": "ems", "region": "Si", "material": "silicon" })),
    ] {
        e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"));
    }
}

#[test]
fn a_wire_with_library_materials_solves_and_its_modes_are_analysed() {
    let mut e = Engine::new();
    assert!(
        e.library().unwrap().names().contains(&"silicon".to_string()),
        "the library of the repository is found: {:?}",
        e.library_notes
    );
    wire(&mut e);
    e.call("set_mesh", &json!({ "h": 1.2e-7, "adapt_scale": 0.4 })).unwrap();
    let r = e.call("solve_run", &json!({})).unwrap_or_else(|err| panic!("solve: {err}"));
    assert_eq!(r["results"][0]["physics"], "ems", "{r}");
    let info = e.call("result_info", &json!({ "physics": "ems" })).unwrap();
    let neff = info[0]["solutions"][0]["eigenvalue"][0].as_f64().expect("a mode with an effective index");
    // 500 x 300 nm silicon wire in oxide at 1550 nm: the documented TE0 index is 2.67 (silicon 3.4757, silica 1.444)
    assert!((neff - 2.67).abs() < 0.03, "TE0 effective index {neff}");

    let a = e.call("analyse", &json!({ "analysis": "polarization", "physics": "ems" })).unwrap();
    assert_eq!(a["analysis"], "polarization");
    let first = &a["rows"][0]["report"];
    assert!(first["class"].as_str().unwrap().contains("TE"), "the fundamental mode is quasi-TE: {first}");
    let typed = e.polarization("ems").unwrap();
    assert_eq!(typed.rows.len(), a["rows"].as_array().unwrap().len());
}

#[test]
fn a_material_the_project_carries_can_extend_a_library_material() {
    let mut e = Engine::new();
    let text = "extends = \"silica\"\nname = \"doped\"\n[optical]\nmodel = \"constant\"\nn = 1.50\n";
    let m = e.add_material("doped", text).unwrap();
    assert_eq!(m["name"], "doped");
    assert!(e.library().unwrap().names().contains(&"doped".to_string()));
    assert!(e.add_material("bad", "[optical]\nmodel = \"nonsense\"\n").is_err());
    e.remove_material("doped").unwrap();
    assert!(!e.library().unwrap().names().contains(&"doped".to_string()));
}

#[test]
fn a_physics_sized_mesh_is_finer_where_the_waves_are_short() {
    let mut e = Engine::new();
    wire(&mut e);
    e.call("set_mesh", &json!({ "h": 4e-7 })).unwrap();
    let coarse = e.call("mesh_build", &json!({})).unwrap()["mesh"]["cells"].as_u64().unwrap();
    e.call("set_mesh", &json!({ "physics_sizing": true })).unwrap();
    let sized: Value = e.call("mesh_build", &json!({})).unwrap();
    assert!(sized["mesh"]["cells"].as_u64().unwrap() > coarse, "sizing from the waves refines the mesh: {sized}");
}

#[test]
fn library_list_and_describe_expose_comparative_descriptions() {
    let mut e = Engine::new();
    let list = e.call("list_materials", &json!({})).unwrap();
    let rows = list.as_array().unwrap();
    let si = rows.iter().find(|r| r["name"] == "silicon").expect("silicon in list");
    assert!(si["description"].as_str().is_some_and(|s| !s.is_empty()), "{si}");
    let salz = e.call("describe_material", &json!({ "name": "silicon_salzberg" })).unwrap();
    assert!(salz["description"].as_str().is_some_and(|s| s.contains("Salzberg")), "{salz}");
}

#[test]
fn remote_catalog_and_offline_import_material_work_and_respect_path_guard() {
    let mut e = Engine::new();
    let all = e.call("list_remote_materials", &json!({})).unwrap();
    let items = all.as_array().unwrap();
    assert!(items.len() >= 20, "expected curated catalog entries: {items:?}");
    let photonics = e.call("list_remote_materials", &json!({ "domain": "photonics" })).unwrap();
    assert!(photonics.as_array().unwrap().iter().all(|r| r["domain"] == "photonics"));
    assert!(e.call("list_remote_materials", &json!({ "domain": "nonsense" })).is_err());

    let root = std::env::temp_dir().join(format!("nyx-engine-import-test-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let save_dir = root.join("saved");
    std::fs::create_dir_all(&save_dir).unwrap();
    e.set_path_root(Some(&root)).unwrap();

    // Escaping folder must be refused by path_guard before any file is written.
    let rii = "REFERENCES: \"Synthetic test\"\nDATA:\n  - type: formula 1\n    wavelength_range: 1.0 2.0\n    coefficients: 0 10.6684293 0.301516485 0.0030434748 1.13475115\n";
    assert!(
        e.call(
            "import_material",
            &json!({ "name": "si_synth", "format": "rii", "content": rii, "extends": "silicon", "folder": "../escaped" }),
        )
        .is_err()
    );

    // Valid import extending `silicon`, added to project and written inside `root`.
    let out = e
        .call(
            "import_material",
            &json!({
                "name": "si_synth",
                "format": "rii",
                "content": rii,
                "description": "Synthetic Sellmeier Si variant extending silicon",
                "extends": "silicon",
                "add_to_project": true,
                "folder": "saved"
            }),
        )
        .unwrap();
    assert_eq!(out["name"], "si_synth");
    assert_eq!(out["added_to_project"], true);
    assert!(save_dir.join("si_synth.toml").is_file(), "TOML should be written into root/saved");
    assert!(e.library().unwrap().names().contains(&"si_synth".to_string()));
    let _ = std::fs::remove_dir_all(&root);
}
