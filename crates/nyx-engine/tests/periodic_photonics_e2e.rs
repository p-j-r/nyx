// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end integration tests for periodic photonics (`PFR-1..PFR-7`) through `Engine::call`:
//! * PFR-1: Fixed-frequency 3D `bloch` study emitting `n_B`, `alpha`, `k_y`, and `phase`.
//! * PFR-2: Per-layer `sweep[*].regions` remap and parent->child 3D periodic lattice conformance.
//! * PFR-3: Single-solve dispersive Poynting `v_g` / `n_g`, component shares, and region shares.
//! * PFR-4: 3D `figure` (`slice3d`, `isosurface3d`) with `repeat_periodic` Bloch phase tiling.
//! * PFR-5: Single-pass bit-exact batch `integrate` and `confinement` (`solutions`, `by_region`, `regions`).
//! * PFR-6: `convergence_check` with `scale_sizes` and physical mode-signature `track`.
//! * PFR-7: Inline diagonal & symmetric tensor `set_material` and UTF-8 BOM-prefixed project JSON & material TOML.

use femcore::material::{Library, MaterialSpec};
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn temp_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-periodic-photonics-e2e-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn scalar_num(scalars: &Value, name: &str) -> Option<f64> {
    scalars.as_array()?.iter().find_map(|pair| {
        let arr = pair.as_array()?;
        (arr.first()?.as_str()? == name).then(|| arr.get(1)?.as_f64()).flatten()
    })
}

#[test]
fn periodic_swg_bloch_figure_batch_reduce_convergence_and_bom_e2e() {
    let dir = temp_dir("swg");
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Value {
        eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"))
    };

    let c0 = 299_792_458.0_f64;
    let lambda0 = 1.55e-6_f64;
    let freq = c0 / lambda0;
    let pitch = 0.24e-6_f64;

    call(&mut e, "project_new", json!({ "title": "3D SWG Bloch E2E", "dim": 2, "h": 0.20e-6 }));
    call(&mut e, "add_rect", json!({ "name": "clad", "min": [-0.85e-6, -0.65e-6], "max": [0.85e-6, 0.65e-6] }));
    call(&mut e, "add_rect", json!({ "name": "core", "min": [-0.225e-6, -0.11e-6], "max": [0.225e-6, 0.11e-6] }));
    call(&mut e, "add_rect", json!({ "name": "tooth", "min": [0.225e-6, -0.11e-6], "max": [0.305e-6, 0.11e-6] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "SiO2", "input": "clad" },
                { "name": "Si_Core", "input": "core" },
                { "name": "Si_Tooth", "input": "tooth" }
            ]
        }),
    );
    call(&mut e, "add_size", json!({ "sel": "region(Si_Core) | region(Si_Tooth)", "h": 0.085e-6 }));

    // PFR-2: Sweep along +y with per-layer `regions` remap (`Si_Tooth -> SiO2` in second half-period)
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 0.20e-6,
            "sweep": [
                { "thickness": 0.5 * pitch, "n": 1 },
                { "thickness": 0.5 * pitch, "n": 1, "regions": { "Si_Tooth": "SiO2" } }
            ]
        }),
    );

    // PFR-7: UTF-8 BOM-prefixed material TOML via `add_material`, `MaterialSpec::from_toml_named`, and `Library::load_dir`
    let bom_toml = "\u{feff}format = 2\nname = \"sio2_bom\"\ndescription = \"BOM-prefixed silica cladding\"\nkind = \"dielectric\"\n\n[optical]\nmodel = \"constant\"\nn = 1.444\n";
    let parsed_spec = MaterialSpec::from_toml_named(bom_toml, "sio2_bom").expect("BOM material TOML must parse");
    assert_eq!(parsed_spec.name, "sio2_bom");
    let mat_dir = dir.join("mats");
    std::fs::create_dir_all(&mat_dir).unwrap();
    std::fs::write(mat_dir.join("sio2_bom.toml"), bom_toml.as_bytes()).unwrap();
    let loaded_lib = Library::load_dir(&mat_dir).expect("BOM material folder must load");
    assert!(loaded_lib.get("sio2_bom").is_some());

    call(&mut e, "add_material", json!({ "name": "sio2_bom", "toml": bom_toml }));

    // PFR-1 & PFR-7: 3D `bloch` study with inline anisotropic diagonal & 6-component symmetric tensor `eps_r`
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "bloch",
            "order": 2,
            "study_params": {
                "freq": freq,
                "n_target": 2.40,
                "n_modes": 3,
                "brillouin_zone": 1
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "SiO2", "material": "sio2_bom" }));
    call(
        &mut e,
        "set_material",
        json!({ "physics": "ems", "region": "Si_Core", "material": { "eps_r": [12.08, 12.08, 12.02] } }),
    );
    call(
        &mut e,
        "set_material",
        json!({ "physics": "ems", "region": "Si_Tooth", "material": { "eps_r": [12.08, 12.08, 12.02, 0.0, 0.0, 0.0] } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "ems", "on": "ymin", "kind": "periodic", "params": { "dest": "ymax", "phase": 0.0 } }),
    );
    call(&mut e, "check_project", json!({}));

    // PFR-7: UTF-8 BOM-prefixed project JSON round-trip via `project_load` and `project_open`
    let proj_json = call(&mut e, "project_json", json!({}))["json"].as_str().unwrap().to_string();
    let bom_proj = format!("\u{feff}{proj_json}");
    let bom_proj_path = dir.join("swg_bom.nyx.json");
    std::fs::write(&bom_proj_path, bom_proj.as_bytes()).unwrap();
    call(&mut e, "project_load", json!({ "json": bom_proj }));
    call(&mut e, "project_open", json!({ "path": bom_proj_path.to_str().unwrap() }));

    // Solve 3D Bloch study
    call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true }));
    let info = call(&mut e, "result_info", json!({ "physics": "ems" }));
    let row = &info.as_array().unwrap()[0];
    assert_eq!(row["study"], json!("bloch"));
    assert_eq!(row["dim"], json!(3));
    assert_eq!(row["order"], json!(2));
    assert_eq!(row["validated"]["state"], json!("exploratory"));

    let sols = row["solutions"].as_array().unwrap();
    assert!(sols.len() >= 2, "expected at least 2 Bloch modes, got {}", sols.len());
    let nb_1 = sols[0]["eigenvalue"].as_array().unwrap()[0].as_f64().unwrap();
    assert!(nb_1 > 1.6 && nb_1 < 3.48, "unexpected fundamental SWG Bloch index n_B = {nb_1}");

    // PFR-3: Single-solve Poynting v_g, n_g, component shares, and domain-named region shares
    let scalars = &row["scalars"];
    let vg_1 = scalar_num(scalars, "v_g 1 (m/s)").expect("v_g 1 (m/s) scalar");
    let ng_1 = scalar_num(scalars, "n_g 1").expect("n_g 1 scalar");
    assert!(vg_1 > 0.0 && vg_1 < c0, "expected forward positive group velocity, got {vg_1:e}");
    assert!(ng_1 > nb_1 && ng_1 < 6.5, "expected waveguide group index n_g > n_B, got n_g={ng_1}, n_B={nb_1}");

    let ex_sh = scalar_num(scalars, "Ex share 1").expect("Ex share 1");
    let ey_sh = scalar_num(scalars, "Ey share 1").expect("Ey share 1");
    let ez_sh = scalar_num(scalars, "Ez share 1").expect("Ez share 1");
    assert!(((ex_sh + ey_sh + ez_sh) - 1.0).abs() < 1e-10, "component shares must sum to 1");
    assert!(ex_sh > 0.40, "fundamental wide-core mode should be quasi-TE (Ex share = {ex_sh})");

    let sh_core = scalar_num(scalars, "share in Si_Core 1").expect("share in Si_Core 1");
    let sh_tooth = scalar_num(scalars, "share in Si_Tooth 1").expect("share in Si_Tooth 1");
    let sh_clad = scalar_num(scalars, "share in SiO2 1").expect("share in SiO2 1");
    assert!(sh_core > 0.35, "expected guided mode confined in Si_Core, got {sh_core}");
    assert!(((sh_core + sh_tooth + sh_clad) - 1.0).abs() < 1e-10, "region shares must sum to 1");
    assert_eq!(scalar_num(scalars, "layer fraction 1"), Some(0.0));

    // PFR-4: 3D `figure` (`slice3d` PNG + `isosurface3d` SVG) with `repeat_periodic` Bloch phase tiling
    let png_path = dir.join("swg_slice3d.png");
    let fig_slice = call(
        &mut e,
        "figure",
        json!({
            "physics": "ems",
            "path": png_path.to_str().unwrap(),
            "kind": "slice3d",
            "modes": [1],
            "expr": "re(Ex)",
            "plane": { "axis": "z", "at": 0.0 },
            "repeat_periodic": 3,
            "cmap": "coolwarm",
            "show_edges": "regions",
            "camera": { "azimuth_deg": 35.0, "elevation_deg": 28.0 }
        }),
    );
    assert_eq!(fig_slice["kind"], json!("slice3d"));
    let png_bytes = std::fs::read(&png_path).unwrap();
    assert!(png_bytes.starts_with(b"\x89PNG\r\n\x1a\n") && png_bytes.len() > 500);

    let svg_path = dir.join("swg_iso3d.svg");
    let fig_iso = call(
        &mut e,
        "figure",
        json!({
            "physics": "ems",
            "path": svg_path.to_str().unwrap(),
            "kind": "isosurface3d",
            "modes": [1],
            "expr": "abs2(E)",
            "isosurfaces": [0.3, 0.7],
            "repeat_periodic": 2
        }),
    );
    assert_eq!(fig_iso["kind"], json!("isosurface3d"));
    let svg_text = std::fs::read_to_string(&svg_path).unwrap();
    assert!(svg_text.contains("<svg") && svg_text.len() > 200);

    // PFR-5: Single-pass bit-exact batch `integrate` and `confinement`
    let batch_int = call(
        &mut e,
        "integrate",
        json!({
            "physics": "ems",
            "solutions": "all",
            "by_region": true,
            "expr": {
                "ex2": "abs2(Ex)",
                "ez2": "abs2(Ez)",
                "e_norm2": "abs2(Ex) + abs2(Ey) + abs2(Ez)"
            }
        }),
    );
    assert_eq!(batch_int["batch"], json!(true));
    assert_eq!(batch_int["by_region"], json!(true));
    let single_int =
        call(&mut e, "integrate", json!({ "physics": "ems", "solution": 0, "region": "Si_Core", "expr": "abs2(Ex)" }));
    assert_eq!(
        batch_int["solutions"][0]["regions"]["Si_Core"]["ex2"].as_f64().unwrap().to_bits(),
        single_int["value"].as_f64().unwrap().to_bits(),
        "batch integrate must match standalone integrate bit-for-bit"
    );

    let batch_conf = call(
        &mut e,
        "confinement",
        json!({
            "physics": "ems",
            "solutions": "all",
            "regions": ["Si_Core", "Si_Tooth", "SiO2"],
            "fields": ["Ex", "Ey", "Ez"],
            "weight": "eps_r_xx"
        }),
    );
    assert_eq!(batch_conf["batch"], json!(true));
    let single_conf = call(
        &mut e,
        "confinement",
        json!({
            "physics": "ems",
            "solution": 0,
            "region": "Si_Core",
            "fields": ["Ex", "Ey", "Ez"],
            "weight": "eps_r_xx"
        }),
    );
    assert_eq!(
        batch_conf["solutions"][0]["regions"]["Si_Core"].as_f64().unwrap().to_bits(),
        single_conf["fraction"].as_f64().unwrap().to_bits(),
        "batch confinement must match standalone confinement bit-for-bit"
    );

    // PFR-6: `convergence_check` dry_run + full run with `scale_sizes = true` and physical mode-signature `track`
    let dry = call(
        &mut e,
        "convergence_check",
        json!({ "physics": "ems", "factor": 0.78, "scale_sizes": true, "dry_run": true }),
    );
    assert_eq!(dry["dry_run"], json!(true));
    assert_eq!(dry["scale_sizes"], json!(true));
    assert!(dry["cells"]["refined"].as_u64().unwrap() > dry["cells"]["coarse"].as_u64().unwrap());

    let conv = call(
        &mut e,
        "convergence_check",
        json!({
            "physics": "ems",
            "factor": 0.80,
            "scale_sizes": true,
            "tolerance": 0.10,
            "track": {
                "n_b": nb_1,
                "region": "Si_Core",
                "min_share": 0.35,
                "component": "Ex",
                "min_component": 0.40,
                "max_layer_fraction": 0.05
            }
        }),
    );
    assert_eq!(conv["matching"], json!("track"));
    assert_eq!(conv["converged"], json!(true), "SWG TE0 Bloch mode convergence failed: {conv}");
    assert_eq!(conv["validated"]["state"], json!("checked"));

    let info_after = call(&mut e, "result_info", json!({ "physics": "ems" }));
    assert_eq!(info_after.as_array().unwrap()[0]["validated"]["state"], json!("checked"));

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn periodic_3d_lattice_mesher_conformance_with_bloch_e2e() {
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Value {
        eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"))
    };

    let c0 = 299_792_458.0_f64;
    let pitch = 0.22e-6_f64;
    call(&mut e, "project_new", json!({ "title": "3D Lattice Periodic Bloch", "dim": 3, "h": 0.22e-6 }));
    call(
        &mut e,
        "add_box",
        json!({ "name": "clad", "min": [-0.65e-6, 0.0, -0.50e-6], "max": [0.65e-6, pitch, 0.50e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "core", "min": [-0.22e-6, 0.0, -0.11e-6], "max": [0.22e-6, pitch, 0.11e-6] }),
    );
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "SiO2", "input": "clad" },
                { "name": "Si", "input": "core" }
            ]
        }),
    );
    call(&mut e, "add_size", json!({ "sel": "region(Si)", "h": 0.10e-6 }));
    call(
        &mut e,
        "add_periodic",
        json!({
            "sel_a": "outer & plane(y = 0)",
            "sel_b": format!("outer & plane(y = {pitch:e})"),
            "translation": [0.0, pitch, 0.0]
        }),
    );
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "bloch",
            "order": 2,
            "study_params": { "freq": c0 / 1.55e-6, "n_target": 2.35, "n_modes": 2 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "SiO2", "material": { "eps_r": 2.085 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "Si", "material": { "eps_r": 12.08 } }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "outer & plane(y = 0)",
            "kind": "periodic",
            "params": { "dest": format!("outer & plane(y = {pitch:e})"), "phase": 0.0 }
        }),
    );

    call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true }));
    let info = call(&mut e, "result_info", json!({ "physics": "ems" }));
    let nb = info[0]["solutions"][0]["eigenvalue"][0].as_f64().unwrap();
    assert!(nb > 1.8 && nb < 3.45, "3D parent->child lattice Bloch solve produced n_B = {nb}");
}
