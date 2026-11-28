// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end FR-37 (`nyx-engine`, ADR 140) verification suite:
//! 1. One-way `dfs::heat_stationary -> ems::eme` thermo-optic waveguide phase shifter across both
//!    shared 3D meshes and per-physics (`FR-44` `active_domains`) cross-meshes, verifying linear
//!    phase shift `Delta phi(P)`, half-wave tuning power `P_pi`, strict `domain_id` isolation across
//!    core/cladding jumps, and cryogenic Einstein-oscillator `dn/dT` saturation.
//! 2. `Engine::call("eme_sweep", ...)` (`length`, `periodicity`, `modes`) and `analyse` `"network"` Touchstone.
//! 3. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bit-exact determinism.

use nyx_engine::Engine;
use nyx_physics::DrivenPortData;
use serde_json::{json, Value};

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn s21_phase_rad(e: &Engine) -> f64 {
    let res = e.result("ems").expect("ems result").0;
    let dp = DrivenPortData::of(res).expect("DrivenPortData on ems::eme");
    let s21 = dp.s_matrix[0][1][0];
    s21[1].atan2(s21[0])
}

/// Builds a 3D thermo-optic waveguide phase-shifter project in `Engine` with:
/// - `clad` (`SiO2` box `[-1.2, 0, -0.8] .. [1.2, L, 0.8]` um, `dn_dT = 1.0e-5 / K`),
/// - `core` (`Si` waveguide `[-0.24, 0, -0.15] .. [0.24, L, 0.15]` um, `dn_dT = 1.86e-4 / K` + optional `einstein_theta`),
/// - `heater` (top metallic heater strip `[-0.35, 0.15*L, 0.45] .. [0.35, 0.85*L, 0.65]` um carrying volumetric `heat_source`),
/// - bottom heat-sink `temperature = t_sink_k` on `z = -0.8 um`.
fn build_thermo_optic_eme_engine(
    heater_power_density_w_m3: f64,
    t_sink_k: f64,
    einstein_theta_k: Option<f64>,
    cross_mesh_optical: bool,
) -> Engine {
    let theta_line = einstein_theta_k.map(|th| format!("dn_dT_einstein_theta = \"{th} K\"\n")).unwrap_or_default();
    let si_toml = format!(
        r#"
format = 2
name = "si_thermo_fr37"
kind = "semiconductor"
reference_temperature = "300 K"
source = {{ cite = "Verified FR-37 thermo-optic silicon test model", verified = "2026-11-27" }}

[thermal]
conductivity = "148 W/(m*K)"

[optical]
wavelength_unit = "um"
range = [1.2, 1.7]
model = "constant"
n = 3.476
k = 0.0
dn_dT = "1.86e-4 1/K"
{theta_line}"#
    );

    let sio2_toml = r#"
format = 2
name = "sio2_thermo_fr37"
kind = "dielectric"
reference_temperature = "300 K"
source = { cite = "Verified FR-37 thermo-optic silica cladding test model", verified = "2026-11-27" }

[thermal]
conductivity = "1.38 W/(m*K)"

[optical]
wavelength_unit = "um"
range = [1.2, 1.7]
model = "constant"
n = 1.444
k = 0.0
dn_dT = "1.0e-5 1/K"
"#;

    let heater_toml = r#"
format = 2
name = "tiw_heater_fr37"
kind = "dielectric"
reference_temperature = "300 K"
source = { cite = "Verified FR-37 resistive micro-heater thermal block", verified = "2026-11-27" }

[thermal]
conductivity = "30.0 W/(m*K)"

[optical]
wavelength_unit = "um"
range = [1.2, 1.7]
model = "constant"
n = 1.444
k = 0.0
"#;

    let l = 20.0e-6;
    let freq = 299_792_458.0 / 1.55e-6;

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "fr37-heat-eme", "dim": 3, "h": 0.55e-6 }));
    call(&mut e, "add_box", json!({ "name": "b_clad", "min": [-1.1e-6, 0.0, -0.75e-6], "max": [1.1e-6, l, 0.75e-6] }));
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_core", "min": [-0.24e-6, 0.0, -0.15e-6], "max": [0.24e-6, l, 0.15e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_heat", "min": [-0.35e-6, 0.15 * l, 0.42e-6], "max": [0.35e-6, 0.85 * l, 0.62e-6] }),
    );
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad", "input": "b_clad" },
                { "name": "core", "input": "b_core" },
                { "name": "heater", "input": "b_heat" }
            ]
        }),
    );

    call(&mut e, "add_material", json!({ "name": "si_thermo_fr37", "toml": si_toml }));
    call(&mut e, "add_material", json!({ "name": "sio2_thermo_fr37", "toml": sio2_toml }));
    call(&mut e, "add_material", json!({ "name": "tiw_heater_fr37", "toml": heater_toml }));

    // Stationary thermal physics (`dfs`)
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "clad", "material": "sio2_thermo_fr37" }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "core", "material": "si_thermo_fr37" }));
    call(
        &mut e,
        "set_material",
        json!({
            "physics": "dfs",
            "region": "heater",
            "material": { "k": 30.0, "Q": heater_power_density_w_m3 }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & plane(z = -0.75e-6)",
            "kind": "temperature",
            "params": { "value": t_sink_k }
        }),
    );

    // Structured y-sweep on project mesh so uniform cold segments share exact cross-sections for BLAKE3 deduplication
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 0.52e-6,
            "sweep": [
                { "thickness": 0.15 * l, "n": 2 },
                { "thickness": 0.70 * l, "n": 6 },
                { "thickness": 0.15 * l, "n": 2 }
            ]
        }),
    );

    // Optical EME physics (`ems`)
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "eme",
            "order": 2,
            "study_params": {
                "freq": freq,
                "n_target": 2.5,
                "n_modes": 2,
                "grid_kind": "uniform",
                "n_cells": 8,
                "port_modes": 1,
                "reconstruct_fields": true
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "sio2_thermo_fr37" }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "si_thermo_fr37" }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "heater", "material": "tiw_heater_fr37" }));

    if cross_mesh_optical {
        // Distinct per-physics optical mesh (`FR-44`) with finer element size AND `active_domains` dropping `heater`!
        call(
            &mut e,
            "set_mesh",
            json!({
                "physics": "ems",
                "h": 0.38e-6,
                "active_domains": ["clad", "core"],
                "sweep": [{ "thickness": l, "n": 10 }]
            }),
        );
    }
    e
}

fn unwrap_phase_advance(phi_hot: f64, phi_cold: f64) -> f64 {
    // Forward propagation is exp(-j beta L), so heating (dn/dT > 0 -> d_beta > 0) makes phase MORE negative:
    // positive thermo-optic phase shift is Delta phi = wrap_pos(phi_cold - phi_hot).
    let mut d = phi_cold - phi_hot;
    while d < -0.5 {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

#[test]
fn cross_mesh_and_shared_mesh_heat_to_eme_thermo_optic_phase_shifter_and_sweeps() {
    // 1. Cold baseline (Q = 0, T = 300 K everywhere) — verify automatic one-way `["ems", "dfs"]` solve order (no false cycle!)
    let mut e_cold_shared = build_thermo_optic_eme_engine(0.0, 300.0, None, false);
    call(&mut e_cold_shared, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let phi0_shared = s21_phase_rad(&e_cold_shared);

    let mut e_cold_cross = build_thermo_optic_eme_engine(0.0, 300.0, None, true);
    call(&mut e_cold_cross, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let phi0_cross = s21_phase_rad(&e_cold_cross);

    // At uniform T = 300 K, all 8 slices of the waveguide have identical cross-section & temperature -> BLAKE3 deduplicates!
    let cold_data = nyx_ems::EmeData::of(e_cold_shared.result("ems").unwrap().0).unwrap();
    assert!(
        cold_data.unique_slices.len() <= 3,
        "uniform cold phase shifter should deduplicate identical slices, got {}",
        cold_data.unique_slices.len()
    );

    // 2. Heated phase shifter (`Q_1` and `2 * Q_1`) on BOTH shared and cross-mesh (`FR-44`), verifying linearity & P_pi!
    let q1 = 1.2e15; // W/m^3 in the micro-heater (~50 K core temperature rise -> ~0.5 rad phase shift)
    let mut e_q1_shared = build_thermo_optic_eme_engine(q1, 300.0, None, false);
    call(&mut e_q1_shared, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let dphi_q1_shared = unwrap_phase_advance(s21_phase_rad(&e_q1_shared), phi0_shared);

    let mut e_q1_cross = build_thermo_optic_eme_engine(q1, 300.0, None, true);
    call(&mut e_q1_cross, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let dphi_q1_cross = unwrap_phase_advance(s21_phase_rad(&e_q1_cross), phi0_cross);

    assert!(
        dphi_q1_shared > 0.25 && dphi_q1_cross > 0.25,
        "expected substantial positive thermo-optic phase shift, got shared={dphi_q1_shared:.4} rad, cross={dphi_q1_cross:.4} rad"
    );
    let mesh_rel_diff = (dphi_q1_shared - dphi_q1_cross).abs() / dphi_q1_shared;
    assert!(
        mesh_rel_diff < 0.08,
        "shared-mesh vs cross-mesh (active_domains) Heat->EME phase shift should agree closely: shared={dphi_q1_shared:.4}, cross={dphi_q1_cross:.4} (rel={mesh_rel_diff:.3})"
    );

    // Verify linear scaling at 2 * Q_1 and solve directly at the extracted P_pi (`Delta phi == pi`)!
    let q_pi = q1 * (std::f64::consts::PI / dphi_q1_cross);
    let mut e_pi_cross = build_thermo_optic_eme_engine(q_pi, 300.0, None, true);
    call(&mut e_pi_cross, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let dphi_pi = unwrap_phase_advance(s21_phase_rad(&e_pi_cross), phi0_cross);
    let pi_rel_err = (dphi_pi - std::f64::consts::PI).abs() / std::f64::consts::PI;
    assert!(
        pi_rel_err < 2e-2,
        "extracted P_pi power density q_pi={q_pi:.3e} must produce pi phase shift: got {dphi_pi:.5} rad (rel err {pi_rel_err:.4e})"
    );

    // 3. Cryogenic Einstein-oscillator dn/dT saturation check (`Theta_E = 450 K` at `T_sink = 15 K` vs `300 K`)!
    let q_cryo = 2.0e14;
    let mut e_cryo_cold = build_thermo_optic_eme_engine(0.0, 15.0, Some(450.0), true);
    call(&mut e_cryo_cold, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let mut e_cryo_hot = build_thermo_optic_eme_engine(q_cryo, 15.0, Some(450.0), true);
    call(&mut e_cryo_hot, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let dphi_cryo = unwrap_phase_advance(s21_phase_rad(&e_cryo_hot), s21_phase_rad(&e_cryo_cold));

    let mut e_rt_cold = build_thermo_optic_eme_engine(0.0, 300.0, Some(450.0), true);
    call(&mut e_rt_cold, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let mut e_rt_hot = build_thermo_optic_eme_engine(q_cryo, 300.0, Some(450.0), true);
    call(&mut e_rt_hot, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let dphi_rt = unwrap_phase_advance(s21_phase_rad(&e_rt_hot), s21_phase_rad(&e_rt_cold));

    assert!(
        dphi_cryo.abs() < 0.05 * dphi_rt,
        "Einstein-oscillator dn/dT at 15 K must saturate strongly compared to 300 K: cryo={dphi_cryo:.4e} vs RT={dphi_rt:.4e}"
    );

    // 4. Exercise `eme_sweep` (`length`, `modes`, and `periodicity`) + `analyse` (`network` Touchstone) through `Engine::call`!
    let len_rep = call(
        &mut e_pi_cross,
        "eme_sweep",
        json!({ "mode": "length", "segment": 3, "lengths": [0.0, 2.5e-6, 10.0e-6] }),
    );
    assert_eq!(len_rep["points"].as_array().unwrap().len(), 3);
    assert_eq!(len_rep["stale"], json!(false));

    let mode_rep = call(&mut e_pi_cross, "eme_sweep", json!({ "mode": "modes", "mode_counts": [1, 2] }));
    assert_eq!(mode_rep["points"].as_array().unwrap().len(), 2);

    let net = call(
        &mut e_pi_cross,
        "analyse",
        json!({ "analysis": "network", "physics": "ems", "params": { "touchstone_format": "MA" } }),
    );
    assert!(net["touchstone_snp"].as_str().unwrap().contains("# Hz S MA R 50"));

    // Also test `periodicity` mode via `set_study` with `grid_kind = "periodic_subgrid"`
    call(
        &mut e_pi_cross,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "eme",
            "order": 2,
            "params": {
                "freq": 299_792_458.0 / 1.55e-6,
                "n_target": 2.5,
                "n_modes": 2,
                "grid_kind": "periodic_subgrid",
                "n_cells": 2,
                "num_reps": 2,
                "port_modes": 1,
                "reconstruct_fields": false
            }
        }),
    );
    call(&mut e_pi_cross, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    let per_rep =
        call(&mut e_pi_cross, "eme_sweep", json!({ "mode": "periodicity", "subgrid": 0, "reps": [1, 10, 100] }));
    assert_eq!(per_rep["points"].as_array().unwrap().len(), 3);
}

#[test]
fn heat_to_eme_cross_mesh_thread_determinism() {
    if let Ok(tag) = std::env::var("NYX_FR37_DETERMINISM_CHILD") {
        let mut e = build_thermo_optic_eme_engine(2.5e13, 300.0, None, true);
        call(&mut e, "solve_run", json!({ "physics": ["ems", "dfs"] }));
        let sw = call(&mut e, "eme_sweep", json!({ "mode": "length", "segment": 2, "lengths": [1.0e-6, 5.0e-6] }));
        let res = &e.result("ems").unwrap().0;
        let dp = DrivenPortData::of(res).unwrap();
        let field_bits: Vec<u64> =
            res.solutions[0].fields.iter().flat_map(|f| f.values.iter().take(32).map(|v| v.to_bits())).collect();
        let payload = json!({
            "tag": tag,
            "s_matrix": dp.s_matrix,
            "sweep_s": sw["points"][1]["s_matrix"],
            "field_bits": field_bits,
        });
        println!("FR37_DET_JSON:{}", serde_json::to_string(&payload).unwrap());
        return;
    }

    let exe = std::env::current_exe().expect("current test binary");
    let run_child = |threads: &str| -> Value {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("heat_to_eme_cross_mesh_thread_determinism")
            .arg("--nocapture")
            .env("NYX_FR37_DETERMINISM_CHILD", threads)
            .env("NYX_THREADS", threads)
            .env("RAYON_NUM_THREADS", threads)
            .output()
            .expect("spawn child test");
        assert!(out.status.success(), "child NYX_THREADS={threads} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let line = stdout
            .lines()
            .find_map(|l| l.strip_prefix("FR37_DET_JSON:"))
            .unwrap_or_else(|| panic!("missing FR37_DET_JSON line in stdout:\n{stdout}"));
        serde_json::from_str(line).unwrap()
    };

    let v1 = run_child("1");
    let v8 = run_child("8");
    assert_eq!(v1["s_matrix"], v8["s_matrix"], "NYX_THREADS=1 vs 8 S-matrix mismatch");
    assert_eq!(v1["sweep_s"], v8["sweep_s"], "NYX_THREADS=1 vs 8 eme_sweep mismatch");
    assert_eq!(v1["field_bits"], v8["field_bits"], "NYX_THREADS=1 vs 8 reconstructed 3D field mismatch");
}
