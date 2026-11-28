// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 1 cross-FR integration verification suite (`nyx-engine`, ADRs 135–139):
//! 1. Cross-mesh two-way photothermal Anderson loop (`FR-44` + `FR-38`): `ems::driven` on a fine
//!    `active_domains = ["core"]` optical mesh coupled with `dfs::heat_stationary` (`RadiationBC` +
//!    `Convection`) on a broader/coarser `["core", "clad"]` thermal mesh, verifying `project_scalar_conservative`
//!    global power balance (`< 5e-4`) and physical agreement with a single-mesh solve.
//! 2. Cross-mesh electrostatically gated graphene waveguide (`FR-44` + `FR-35` + `FR-38`): `dfs::electrostatic`
//!    (`P2` potential `dfs.V`) on its own mesh transferred onto `ems::modes` with Kubo graphene `sheet_2d`,
//!    verifying Pauli-blocking modulation `Delta alpha(V_g)`, conservative `Q_opt` sheet dissipation, and
//!    single-mesh eigenvalue agreement.
//! 3. Combined bulk gyrotropy (`FR-34`) + Hall-biased `sheet_2d` (`FR-35`) + scattered-field excitation (`FR-36`)
//!    and `Q_opt` export (`FR-38`) in 3D `ems::driven`, verifying `lblt` rejection, general-LU PROM accuracy,
//!    and exact total-field absorption conservation.
//! 4. Subprocess bit-exact determinism across `NYX_THREADS=1` vs `NYX_THREADS=8`.

use nyx_engine::Engine;
use nyx_physics::DrivenPortData;
use serde_json::{json, Value};

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn scalar_num(e: &Engine, physics: &str, name: &str) -> f64 {
    let res = &e.result(physics).expect("result exists").0;
    for (k, v) in &res.scalars {
        if k == name {
            if let nyx_physics::Scalar::Num(x) = v {
                return *x;
            }
        }
    }
    panic!("scalar '{name}' not found on '{physics}': {:?}", res.scalars);
}

fn build_photothermal_engine(cross_mesh: bool) -> Engine {
    let core_toml = r#"
format = 2
name = "thermo_core_p1"
kind = "dielectric"
reference_temperature = "300 K"
density = "3200 kg/m^3"
source = { cite = "Synthetic Phase 1 cross-mesh photothermal resonator core", verified = "2026-11-27" }

[thermal]
conductivity = "5.0 W/(m*K)"
specific_heat = "700 J/(kg*K)"

[dielectric]
eps_free = 10.0
eps_T = { form = "tabulated_T", T = [250.0, 300.0, 500.0], values = [9.0, 10.0, 14.5] }
loss_tangent = 0.05
"#;

    let clad_toml = r#"
format = 2
name = "thermal_cap_p1"
kind = "dielectric"
reference_temperature = "300 K"
density = "2200 kg/m^3"
source = { cite = "Synthetic Phase 1 thermal heat-spreader cladding", verified = "2026-11-27" }

[thermal]
conductivity = "12.0 W/(m*K)"
specific_heat = "750 J/(kg*K)"

[dielectric]
eps_free = 2.25
"#;

    let (a, l, b) = (6.0e-3, 8.0e-3, 3.0e-3);
    let freq = 6.5e9;
    let source_scale = 140.0;

    let mut e = Engine::new();
    let h_proj = if cross_mesh { 2.6e-3 } else { 2.0e-3 };
    call(&mut e, "project_new", json!({ "title": "p1-photothermal", "dim": 3, "h": h_proj }));
    call(&mut e, "add_box", json!({ "name": "core_box", "min": [0.0, 0.0, 0.0], "max": [a, l, b] }));
    call(&mut e, "add_box", json!({ "name": "clad_box", "min": [0.0, l, 0.0], "max": [a, 2.0 * l, b] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "core", "input": "core_box" },
                { "name": "clad", "input": "clad_box" }
            ]
        }),
    );
    call(&mut e, "add_material", json!({ "name": "thermo_core_p1", "toml": core_toml }));
    call(&mut e, "add_material", json!({ "name": "thermal_cap_p1", "toml": clad_toml }));

    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 1,
            "study_params": { "freqs_hz": [freq], "prom_tol": 0.0 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "thermo_core_p1" }));
    if !cross_mesh {
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "thermal_cap_p1" }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "interface(core, clad)", "kind": "pec" }));
    }
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(-y)",
            "kind": "lumped_port",
            "params": { "name": "p1", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }
        }),
    );

    if cross_mesh {
        call(
            &mut e,
            "set_mesh",
            json!({
                "physics": "ems",
                "h": 2.0e-3,
                "active_domains": ["core"]
            }),
        );
    }

    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "core", "material": "thermo_core_p1" }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "clad", "material": "thermal_cap_p1" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer",
            "kind": "convection",
            "params": { "h": 2200.0, "t_ext": 300.0 }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer",
            "kind": "radiation",
            "params": { "emissivity": 0.85, "t_ambient": 300.0, "view_factor": 1.0 }
        }),
    );

    call(
        &mut e,
        "set_coupled_loop",
        json!({ "max_iters": 12, "rtol": 1e-7, "anderson_m": 5, "source_scale": source_scale }),
    );
    call(&mut e, "solve_run", json!({ "physics": ["ems", "dfs"] }));
    e
}

#[test]
fn cross_mesh_photothermal_anderson_loop_conserves_power_and_matches_single_mesh() {
    let source_scale = 140.0;
    let mut e_cross = build_photothermal_engine(true);
    let mut e_single = build_photothermal_engine(false);

    assert!(e_cross.project().has_own_mesh("ems") && !e_cross.project().has_own_mesh("dfs"));
    let em_q_len = e_cross.result("ems").unwrap().0.solutions[0].field("Q_opt").unwrap().values.len();
    let th_t_len = e_cross.result("dfs").unwrap().0.solutions[0].field("T").unwrap().values.len();
    assert_ne!(
        em_q_len, th_t_len,
        "cross-mesh run must use distinct optical (core-only) and thermal (core+clad) meshes"
    );
    let cm_check = call(
        &mut e_cross,
        "integrate",
        json!({
            "expr": {
                "op": "mul",
                "args": [
                    { "field": "Q_opt", "physics": "ems", "solution": 0 },
                    { "field": "T", "physics": "dfs", "solution": 0 }
                ]
            }
        }),
    );
    assert_eq!(cm_check["cross_mesh_transfer"], json!(true));

    let iters = scalar_num(&e_cross, "dfs", "coupled_iterations");
    let rel_res = scalar_num(&e_cross, "dfs", "coupled_rel_residual");
    let p_loss = scalar_num(&e_cross, "ems", "P_loss (W)");
    let t_max_cross = scalar_num(&e_cross, "dfs", "T max (K)");
    assert!(iters <= 8.0, "cross-mesh Anderson loop must converge in <= 8 iterations, got {iters}");
    assert!(rel_res <= 1e-7, "coupled_rel_residual {rel_res:.3e} > 1e-7");
    assert!(t_max_cross > 304.0, "expected measurable photothermal self-heating above 300 K, got {t_max_cross}");

    // 1. Optical Q_opt integral on fine core-only EM mesh equals P_loss (W)
    let q_opt_int =
        call(&mut e_cross, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }))["value"].as_f64().unwrap();
    assert!(((q_opt_int - p_loss) / p_loss).abs() < 1e-10);

    // 2. Conservative L2 transfer of Q_opt onto coarser core+clad thermal mesh preserves total power:
    //    convective + Stefan-Boltzmann radiative rejection on outer thermal boundary == source_scale * P_loss!
    let p_expected = source_scale * p_loss;
    let p_out = call(
        &mut e_cross,
        "integrate",
        json!({
            "physics": "dfs",
            "boundary": "outer",
            "expr": "2200.0 * (T - 300.0) + 0.85 * 5.670374419e-8 * (T^4 - 300.0^4)"
        }),
    )["value"]
        .as_f64()
        .unwrap();
    let bal_err = (p_out - p_expected).abs() / p_expected;
    assert!(
        bal_err < 5e-4,
        "cross-mesh photothermal boundary power balance error {bal_err:.3e} >= 5e-4: P_out = {p_out:.6e} vs {p_expected:.6e}"
    );

    // 3. Compare cross-mesh converged temperature rise and S11 against single-mesh reference
    let t_max_single = scalar_num(&e_single, "dfs", "T max (K)");
    let dt_cross = t_max_cross - 300.0;
    let dt_single = t_max_single - 300.0;
    let dt_rel = (dt_cross - dt_single).abs() / dt_single.max(1.0);
    assert!(
        dt_rel < 0.05,
        "cross-mesh vs single-mesh temperature rise mismatch {:.3}%: dT_cross = {dt_cross:.3} K vs dT_single = {dt_single:.3} K",
        dt_rel * 100.0
    );

    let s11_cross = DrivenPortData::of(e_cross.result("ems").unwrap().0).unwrap().s_matrix[0][0][0];
    let s11_single = DrivenPortData::of(e_single.result("ems").unwrap().0).unwrap().s_matrix[0][0][0];
    let s11_diff = (s11_cross[0] - s11_single[0]).hypot(s11_cross[1] - s11_single[1]);
    assert!(s11_diff < 0.04, "cross-mesh vs single-mesh |S11_cross - S11_single| = {s11_diff:.4} >= 0.04");
    let _ = call(&mut e_single, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }));
}

fn solve_gated_graphene_waveguide(v_top: f64, cross_mesh: bool) -> (Engine, [f64; 2], f64) {
    let si_toml = r#"
format = 2
name = "si_wg_p1"
kind = "dielectric"
source = { cite = "Synthetic Phase 1 silicon waveguide core", verified = "2026-11-27" }

[dielectric]
eps_free = 12.08
"#;
    let sio2_toml = r#"
format = 2
name = "sio2_clad_p1"
kind = "dielectric"
source = { cite = "Synthetic Phase 1 oxide cladding", verified = "2026-11-27" }

[dielectric]
eps_free = 2.085
"#;

    let lambda0 = 1550.0e-9;
    let freq_hz = 299_792_458.0 / lambda0;
    let k0 = 2.0 * std::f64::consts::PI / lambda0;

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p1-gated-graphene", "dim": 2, "h": 0.11e-6 }));
    call(&mut e, "add_rect", json!({ "name": "clad_r", "min": [-1.0e-6, -0.6e-6], "max": [1.0e-6, 0.6e-6] }));
    call(&mut e, "add_rect", json!({ "name": "bot_r", "min": [-0.25e-6, -0.11e-6], "max": [0.25e-6, 0.0] }));
    call(&mut e, "add_rect", json!({ "name": "top_r", "min": [-0.25e-6, 0.0], "max": [0.25e-6, 0.11e-6] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad", "input": "clad_r" },
                { "name": "core_bot", "input": "bot_r" },
                { "name": "core_top", "input": "top_r" }
            ]
        }),
    );
    call(&mut e, "add_material", json!({ "name": "si_wg_p1", "toml": si_toml }));
    call(&mut e, "add_material", json!({ "name": "sio2_clad_p1", "toml": sio2_toml }));

    // Electrostatic gating study (`dfs.V` -> `ems.potential`)
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "clad", "material": "sio2_clad_p1" }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "core_bot", "material": "si_wg_p1" }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "core_top", "material": "si_wg_p1" }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "extreme(-z)", "kind": "ground" }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "extreme(+z)", "kind": "potential", "params": { "value": v_top } }),
    );
    if cross_mesh {
        call(&mut e, "set_mesh", json!({ "physics": "dfs", "h": 0.17e-6 }));
    }

    // Optical waveguide mode study (`ems::modes`) with monolayer Kubo graphene on `interface(core_bot, core_top)`
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "modes",
            "order": 2,
            "study_params": { "freq": freq_hz, "n_target": 2.45, "evanescent": false, "n_modes": 1 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "sio2_clad_p1" }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core_bot", "material": "si_wg_p1" }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core_top", "material": "si_wg_p1" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "interface(core_bot, core_top)",
            "kind": "sheet_2d",
            "params": {
                "material": "graphene",
                "tau_s": 0.2e-12,
                "temperature": 300.0,
                "eps_ox": 9.0,
                "d_ox": 7.0e-9,
                "v_cnp": -0.15
            }
        }),
    );

    call(&mut e, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let eig = e.result("ems").unwrap().0.solutions[0].eigenvalue.expect("complex neff");
    let alpha_np_m = (-2.0 * k0 * eig[1]).max(0.0);
    let alpha_db_um = 10.0 / std::f64::consts::LN_10 * alpha_np_m * 1e-6;
    (e, eig, alpha_db_um)
}

#[test]
fn cross_mesh_electrostatic_gated_graphene_waveguide_and_hall_gyrotropy_composition() {
    // Part A: Cross-mesh `dfs.V` -> `ems.potential` gated graphene waveguide
    let (mut e_low, eig_low, alpha_low) = solve_gated_graphene_waveguide(0.0, true);
    let (e_high, eig_high, alpha_high) = solve_gated_graphene_waveguide(10.0, true);
    let (_e_ref, eig_ref, _alpha_ref) = solve_gated_graphene_waveguide(10.0, false);

    assert!(e_low.project().has_own_mesh("dfs") && !e_low.project().has_own_mesh("ems"));
    assert_ne!(
        e_low.result("dfs").unwrap().0.dofs,
        e_low.result("ems").unwrap().0.solutions[0].fields[0].values.len(),
        "dfs electrostatic and ems optical mode meshes must be independent"
    );
    assert!(alpha_low > 0.02, "expected strong interband graphene absorption at V_top = 0 V, got {alpha_low} dB/um");
    assert!(
        alpha_high < 0.15 * alpha_low,
        "expected Pauli-blocked suppression at V_top = 10 V ({alpha_high} dB/um vs {alpha_low} dB/um)"
    );
    let neff_rel = ((eig_high[0] - eig_ref[0]).hypot(eig_high[1] - eig_ref[1])) / eig_ref[0].hypot(eig_ref[1]);
    assert!(neff_rel < 5e-3, "cross-mesh vs single-mesh gated graphene neff mismatch {neff_rel:.3e} >= 5e-3");
    let _ = eig_low;

    // Verify 2D `sheet_2d` graphene dissipation is exported in `Q_opt` and satisfies waveguide Poynting balance:
    // For a mode normalized to 1 W guided Poynting flux, `int_Omega Q_opt dA == alpha_np_m = -2 * k0 * Im(neff)` (W/m)!
    let k0 = 2.0 * std::f64::consts::PI / 1550.0e-9;
    let alpha_np_low = (-2.0 * k0 * eig_low[1]).max(0.0);
    let q_wg_int =
        call(&mut e_low, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }))["value"].as_f64().unwrap();
    assert!(q_wg_int > 0.0, "expected positive Q_opt integral from absorbing graphene sheet_2d, got {q_wg_int}");
    let wg_poynting_err = (q_wg_int - alpha_np_low).abs() / alpha_np_low;
    assert!(
        wg_poynting_err < 2e-2,
        "2D waveguide graphene Poynting balance mismatch {wg_poynting_err:.3e}: int Q_opt dA = {q_wg_int:.6e} W/m vs 2*k0*|Im(neff)| = {alpha_np_low:.6e} W/m"
    );
    let _ = e_high;

    // Part B: 3D `ems::driven` composing bulk gyrotropy (`FR-34`) + Hall-biased graphene `sheet_2d` (`FR-35`) + `Q_opt` (`FR-38`)
    let gyro_toml = r#"
format = 2
name = "gyro_substrate_p1"
source = { cite = "Synthetic Phase 1 lossy gyrotropic THz medium", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [50.0, 200.0]
model = "constant"
n = 2.0
k = 0.02

[optical.gyrotropy]
model = "constant"
g = [0.0, 0.0, 0.18]
"#;
    let mut eg = Engine::new();
    let (w, l) = (20.0e-6, 15.0e-6);
    let freq_thz = 3.0e12;
    call(&mut eg, "project_new", json!({ "title": "p1-gyro-hall", "dim": 3, "h": 14.0e-6 }));
    call(&mut eg, "add_box", json!({ "name": "b_in", "min": [0.0, 0.0, -l], "max": [w, w, 0.0] }));
    call(&mut eg, "add_box", json!({ "name": "b_out", "min": [0.0, 0.0, 0.0], "max": [w, w, l] }));
    call(
        &mut eg,
        "add_overlay",
        json!({ "regions": [{ "name": "sub_in", "input": "b_in" }, { "name": "sub_out", "input": "b_out" }] }),
    );
    call(&mut eg, "add_material", json!({ "name": "gyro_substrate_p1", "toml": gyro_toml }));

    // First verify `solver = "lblt"` is cleanly rejected when gyrotropy + Hall sheet are active
    call(
        &mut eg,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 1,
            "study_params": { "f_min": 2.8e12, "f_max": 3.2e12, "n_points": 3, "prom_tol": 1e-4, "solver": "lblt" }
        }),
    );
    call(&mut eg, "set_material", json!({ "physics": "ems", "region": "sub_in", "material": "gyro_substrate_p1" }));
    call(&mut eg, "set_material", json!({ "physics": "ems", "region": "sub_out", "material": "gyro_substrate_p1" }));
    call(&mut eg, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pec" }));
    call(&mut eg, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pec" }));
    call(&mut eg, "add_condition", json!({ "physics": "ems", "on": "extreme(-y)", "kind": "pmc" }));
    call(&mut eg, "add_condition", json!({ "physics": "ems", "on": "extreme(+y)", "kind": "pmc" }));
    call(
        &mut eg,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(-z)",
            "kind": "lumped_port",
            "params": { "name": "p1", "Z_0": 188.0, "direction": [1.0, 0.0, 0.0] }
        }),
    );
    call(
        &mut eg,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(+z)",
            "kind": "lumped_port",
            "params": { "name": "p2", "Z_0": 188.0, "direction": [1.0, 0.0, 0.0] }
        }),
    );
    call(
        &mut eg,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "interface(sub_in, sub_out)",
            "kind": "sheet_2d",
            "params": {
                "material": "graphene",
                "mu_c_eV": 0.35,
                "tau_s": 0.5e-12,
                "temperature": 300.0,
                "b_0": 6.0
            }
        }),
    );

    let lblt_err = eg.call("solve_run", &json!({ "physics": "ems" })).unwrap_err().to_string();
    assert!(lblt_err.contains("lblt"), "expected lblt rejection for gyrotropic + Hall system, got: {lblt_err}");

    // Switch to `solver = "auto"` (general LU) and compare PROM (`prom_tol = 1e-4`) vs direct HDM (`prom_tol = 0.0`)
    call(
        &mut eg,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "driven",
            "params": { "f_min": 2.8e12, "f_max": 3.2e12, "n_points": 3, "prom_tol": 1e-4, "solver": "auto" }
        }),
    );
    call(&mut eg, "solve_run", json!({ "physics": "ems" }));
    let s_prom = DrivenPortData::of(eg.result("ems").unwrap().0).unwrap().s_matrix.clone();

    call(
        &mut eg,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "driven",
            "params": { "f_min": 2.8e12, "f_max": 3.2e12, "n_points": 3, "prom_tol": 0.0, "solver": "auto" }
        }),
    );
    call(&mut eg, "solve_run", json!({ "physics": "ems" }));
    let s_hdm = DrivenPortData::of(eg.result("ems").unwrap().0).unwrap().s_matrix.clone();
    for k in 0..3 {
        let diff = (s_prom[k][1][0][0] - s_hdm[k][1][0][0]).hypot(s_prom[k][1][0][1] - s_hdm[k][1][0][1]);
        assert!(diff < 1e-4, "non-symmetric gyrotropy + Hall PROM vs HDM S21 mismatch {diff:.3e} at point {k}");
    }
    let p_loss_tot = scalar_num(&eg, "ems", "P_loss (W)");
    let q_tot_int = call(&mut eg, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }))["value"].as_f64().unwrap();
    assert!(p_loss_tot > 0.0 && ((q_tot_int - p_loss_tot) / p_loss_tot).abs() < 1e-8);

    // Part C: Scattered-field `excitation` (`FR-36`) + `Q_opt` (`FR-38`) total-field heating consistency
    let lossy_toml = r#"
format = 2
name = "lossy_scat_p1"
kind = "dielectric"
source = { cite = "Synthetic Phase 1 scattered-field absorber", verified = "2026-11-27" }

[dielectric]
eps_free = 4.0
loss_tangent = 0.08
"#;
    let vac_toml = r#"
format = 2
name = "vac_bg_p1"
kind = "dielectric"
source = { cite = "Vacuum background", verified = "2026-11-27" }

[dielectric]
eps_free = 1.0
"#;
    let mut es = Engine::new();
    call(&mut es, "project_new", json!({ "title": "p1-scat-heat", "dim": 3, "h": 18.0e-6 }));
    call(
        &mut es,
        "add_box",
        json!({ "name": "bg_b", "min": [-20.0e-6, -20.0e-6, -20.0e-6], "max": [20.0e-6, 20.0e-6, 20.0e-6] }),
    );
    call(
        &mut es,
        "add_box",
        json!({ "name": "obj_b", "min": [-8.0e-6, -8.0e-6, -8.0e-6], "max": [8.0e-6, 8.0e-6, 8.0e-6] }),
    );
    call(
        &mut es,
        "add_overlay",
        json!({ "regions": [{ "name": "bg", "input": "bg_b" }, { "name": "scatterer", "input": "obj_b" }] }),
    );
    call(&mut es, "add_material", json!({ "name": "lossy_scat_p1", "toml": lossy_toml }));
    call(&mut es, "add_material", json!({ "name": "vac_bg_p1", "toml": vac_toml }));
    call(
        &mut es,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 1,
            "study_params": { "freqs_hz": [freq_thz], "prom_tol": 0.0 }
        }),
    );
    call(&mut es, "set_material", json!({ "physics": "ems", "region": "bg", "material": "vac_bg_p1" }));
    call(&mut es, "set_material", json!({ "physics": "ems", "region": "scatterer", "material": "lossy_scat_p1" }));
    call(
        &mut es,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "scatterer",
            "kind": "plane_wave",
            "params": {
                "name": "pw1",
                "direction": [0.0, 0.0, 1.0],
                "polarization": [1.0, 0.0, 0.0],
                "amplitude": 1000.0
            }
        }),
    );
    call(&mut es, "solve_run", json!({ "physics": "ems" }));
    let p_abs_pw = scalar_num(&es, "ems", "pw1: P_abs (W)");
    let p_loss_pw = scalar_num(&es, "ems", "P_loss (W)");
    let q_pw_int = call(&mut es, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }))["value"].as_f64().unwrap();
    assert!(p_abs_pw > 0.0);
    assert!(
        ((p_loss_pw - p_abs_pw) / p_abs_pw).abs() < 1e-8,
        "scattered-field total-field P_loss {p_loss_pw:.6e} must match pw1: P_abs {p_abs_pw:.6e}"
    );
    assert!(((q_pw_int - p_loss_pw) / p_loss_pw).abs() < 1e-8);
}

#[test]
fn phase1_cross_mesh_coupled_workflows_are_deterministic_across_nyx_threads_1_and_8() {
    const ENV_CHILD: &str = "NYX_PHASE1_INT_DET_CHILD";
    if std::env::var_os(ENV_CHILD).is_some() {
        let e1 = build_photothermal_engine(true);
        let t_bits: Vec<u64> =
            e1.result("dfs").unwrap().0.solutions[0].field("T").unwrap().values.iter().map(|v| v.to_bits()).collect();
        let q_bits: Vec<u64> = e1.result("ems").unwrap().0.solutions[0]
            .field("Q_opt")
            .unwrap()
            .values
            .iter()
            .map(|v| v.to_bits())
            .collect();
        let (_e2, eig, _alpha) = solve_gated_graphene_waveguide(6.0, true);
        let payload = json!({
            "t": t_bits,
            "q": q_bits,
            "eig": [eig[0].to_bits(), eig[1].to_bits()],
        });
        println!("PHASE1_DET={}", serde_json::to_string(&payload).unwrap());
        return;
    }

    let exe = std::env::current_exe().expect("test binary path");
    let run_with_threads = |nt: &str| -> String {
        let out = std::process::Command::new(&exe)
            .args([
                "--exact",
                "phase1_cross_mesh_coupled_workflows_are_deterministic_across_nyx_threads_1_and_8",
                "--nocapture",
            ])
            .env(ENV_CHILD, "1")
            .env("NYX_THREADS", nt)
            .output()
            .expect("spawn child test");
        assert!(
            out.status.success(),
            "child NYX_THREADS={nt} failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        stdout
            .lines()
            .find_map(|l| l.strip_prefix("PHASE1_DET="))
            .unwrap_or_else(|| panic!("missing PHASE1_DET in child output:\n{stdout}"))
            .to_string()
    };

    let b1 = run_with_threads("1");
    let b8 = run_with_threads("8");
    assert_eq!(b1, b8, "Phase 1 cross-mesh coupled solves must be 0-ULP bit-identical between NYX_THREADS=1 and 8");
}
