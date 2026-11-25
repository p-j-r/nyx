// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Monty (`nyx-script`) end-to-end verification of Phase A accuracy features (`AF-1..AF-6`) through `.nys` scripts.

use nyx_engine::Engine;
use nyx_script::{host, Limits, Options};
use std::path::PathBuf;

fn worker() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-script-worker"))
}

fn temp_root(name: &str) -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!(
        "nyx-af-script-{name}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn monty_script_drives_af1_through_af6_commands_and_static_check() {
    let script = r#"
nyx.project_new(title="Monty AF1..AF6", dim=3)
nyx.add_polygons(
    name="p1",
    polygons=[[[-40e-6, -25e-6], [-8e-6, -25e-6], [-8e-6, 25e-6], [-40e-6, 25e-6]]],
    arc_tol=150e-9,
)
nyx.add_polygons(
    name="p2",
    polygons=[[[8e-6, -25e-6], [40e-6, -25e-6], [40e-6, 25e-6], [8e-6, 25e-6]]],
    arc_tol=150e-9,
)
nyx.add_sheet(name="s1", input="p1", z=0.0, thickness=0.0, etch_depth=0.0)
nyx.add_sheet(name="s2", input="p2", z=0.0, thickness=0.0, etch_depth=0.0)
stack = nyx.add_chip_stack(lateral_pad=55e-6, t_sub=70e-6, h_vac_top=70e-6, h_vac_bot=0.0)
nyx.set_mesh(
    h=42e-6,
    curved=True,
    gap_cells=2.0,
    adapt_grade=0.35,
    adapt_max_cells=50000,
    corner_grading=False,
)
nyx.add_physics(
    physics="dfs",
    study="electrostatic",
    order=2,
    study_params={
        "adapt_passes": 2,
        "adapt_rel_tol": 1e-3,
        "adapt_mode": "sheet_2d_z",
        "adapt_estimator": "flux_jump_goal",
        "flux_bounds": "midpoint",
        "freq": 13e9,
    },
)
nyx.set_material(physics="dfs", region="Si", material="silicon")
nyx.set_material(physics="dfs", region="VacTop", material={"eps_r": 1.0})
nyx.add_condition(physics="dfs", on="s1", kind="terminal", params={"name": "c1"})
nyx.add_condition(physics="dfs", on="s2", kind="terminal", params={"name": "c2"})
nyx.add_condition(physics="dfs", on="outer", kind="far_field", params={"order": 1})
nyx.solve_run()
nyx.analyse(analysis="capacitance_matrix", physics="dfs")
"#;

    let root = temp_root("e2e");
    let mut e = Engine::new();
    let opts = Options { root: root.clone(), worker: Some(worker()) };
    let limits = Limits::default();
    host::check(&mut e, "af_e2e.nys", script, &limits, &opts)
        .expect("Monty static check must accept all AF-1..AF-6 commands and kwargs");

    let out = host::run(&mut e, "af_e2e.nys", script, &limits, &opts).unwrap();
    let v = &out.value;
    assert_eq!(v["terminals"], serde_json::json!(["c1", "c2"]));
    assert!(v["maxwell_lower_matrix_F"].is_array(), "missing Prager-Synge lower matrix: {v}");
    assert!(v["maxwell_upper_matrix_F"].is_array(), "missing Prager-Synge upper matrix: {v}");
    assert!(v["conductance_matrix_S"].is_array(), "missing 13 GHz conductance matrix: {v}");
    assert!(v["q_dielectric_matrix"].is_array(), "missing 13 GHz Q matrix: {v}");
    assert_eq!(v["adapt_history"].as_array().unwrap().len(), 2, "expected 2 AMR passes: {v}");

    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn monty_script_drives_outward_absorbing_layers_and_pml_diagnostics() {
    let script = r#"
nyx.project_new(title="Monty Outward PML", dim=2)
nyx.add_circle(name="clad_c", center=[0.0, 0.0], r=2.0e-6)
nyx.add_circle(name="core_c", center=[0.0, 0.0], r=1.0e-6)
nyx.add_overlay(
    name="regions",
    regions=[
        {"name": "cladding", "input": "clad_c"},
        {"name": "core", "input": "core_c"},
    ],
)
nyx.add_absorbing_layer(
    name="pml",
    coord="Cylindrical",
    thickness=0.8e-6,
    layers=5,
    grading=2.0,
    layer_grading=2.0,
)
nyx.set_mesh(h=0.4e-6, pml_mode="auto")
nyx.add_physics(
    physics="ems",
    study="modes",
    order=2,
    study_params={"freq": 299792458.0 / 1.55e-6, "n_target": 2.4, "n_modes": 2},
)
nyx.set_material(physics="ems", region="core", material={"eps_r": 12.0, "mu_r": 1.0})
nyx.set_material(physics="ems", region="cladding*", material={"eps_r": 1.0, "mu_r": 1.0})
nyx.add_pml_conditions(physics="ems", auto_scale=True)
nyx.mesh_build()
nyx.pml_diagnostics()
"#;

    let root = temp_root("pml_e2e");
    let mut e = Engine::new();
    let opts = Options { root: root.clone(), worker: Some(worker()) };
    let limits = Limits::default();
    host::check(&mut e, "pml_e2e.nys", script, &limits, &opts)
        .expect("Monty static check must accept add_absorbing_layer, set_mesh(pml_mode=...), add_pml_conditions, and pml_diagnostics");

    let out = host::run(&mut e, "pml_e2e.nys", script, &limits, &opts).unwrap();
    let v = &out.value;
    assert!(v["element_counts"]["quad"].as_u64().unwrap_or(0) > 0, "expected structured Quad elements: {v}");
    assert!(v["min_pml_scaled_jacobian"].as_f64().unwrap_or(0.0) > 0.5, "expected positive PML Jacobian: {v}");
    assert_eq!(v["s_dist_bounds"], serde_json::json!([0.0, 1.0]), "{v}");

    let _ = std::fs::remove_dir_all(&root);
}
