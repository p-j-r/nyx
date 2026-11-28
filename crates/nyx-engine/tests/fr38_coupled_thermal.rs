// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end FR-38 verification suite (`nyx-engine`, ADR 139):
//! 1. Default cycle rejection in `Registry::solve_order` / `solve_run` when `coupled_loop` is absent.
//! 2. Two-way Anderson electrothermal (`dfs::conduction <-> dfs::heat_stationary`) self-consistent loop
//!    with temperature-dependent $\sigma(T)$, $k(T)$, and `Radiation` + `Convection` BCs converging in
//!    `<= 6` outer iterations with exact boundary flux vs Joule power conservation (`< 1e-4`).
//! 3. Two-way Anderson photothermal (`ems::driven <-> dfs::heat_stationary`) self-consistent loop on a
//!    3D thermo-dielectric lossy resonator converging in `<= 6` outer iterations with exact electromagnetic
//!    and thermal energy conservation (`< 1e-4`).
//! 4. Determinism across `NYX_THREADS=1` vs `NYX_THREADS=8`.

use nyx_engine::{Engine, Project};
use nyx_physics::{DrivenPortData, ProblemSpec};
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

#[test]
fn anderson_two_way_electrothermal_loop_converges_and_conserves_energy() {
    // 2D metallic heater strip [0, L] x [0, W] with strong negative d(sigma)/dT (TCR):
    // sigma(T) drops from ~2.05e5 S/m at 300 K to 1.0e5 S/m at 600 K, while thermal conductivity k(T)
    // varies mildly and heat leaves all outer faces via Convection + Stefan-Boltzmann T^4 Radiation!
    let mat_toml = r#"
format = 2
name = "tcr_heater"
kind = "metal"
reference_temperature = "300 K"
density = "8000 kg/m^3"
source = { cite = "Synthetic verified FR-38 electrothermal test material", verified = "2026-11-27" }

[thermal]
conductivity = { form = "tabulated_T", T = [250.0, 900.0], values = [40.0, 53.0], unit = "W/(m*K)" }
specific_heat = "450 J/(kg*K)"

[electrical]
conductivity = { form = "tabulated_T", T = [250.0, 600.0, 900.0], values = [2.2e5, 1.0e5, 0.7e5], unit = "S/m" }
"#;

    let (l, w) = (2.0e-3, 0.5e-3);
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "electrothermal-heater", "dim": 2, "h": 2.5e-4 }));
    call(&mut e, "add_rect", json!({ "name": "bar", "min": [0.0, 0.0], "max": [l, w] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "heater", "input": "bar" }] }));
    call(&mut e, "add_material", json!({ "name": "tcr_heater", "toml": mat_toml }));

    // Primary `dfs` problem: nonlinear `heat_stationary` with Convection + RadiationBC on outer boundary
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "heater", "material": "tcr_heater" }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer", "kind": "convection", "params": { "h": 800.0, "t_ext": 300.0 } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer",
            "kind": "radiation",
            "params": { "emissivity": 0.8, "t_ambient": 300.0, "view_factor": 1.0 }
        }),
    );

    // Auxiliary `dfs` `conduction` problem driving Joule self-heating V <-> T
    let mut cond_spec = ProblemSpec::new("dfs", 2, "conduction");
    cond_spec.materials.insert("heater".into(), nyx_physics::MaterialAssign::Named("tcr_heater".into()));
    cond_spec.bcs.push(nyx_physics::BcSpec {
        on: "outer & plane(x = 0)".into(),
        kind: "ground".into(),
        params: json!({}),
    });
    cond_spec.bcs.push(nyx_physics::BcSpec {
        on: format!("outer & plane(x = {l})"),
        kind: "voltage".into(),
        params: json!({ "value": 0.12 }),
    });

    call(
        &mut e,
        "set_coupled_loop",
        json!({
            "max_iters": 12,
            "rtol": 1e-8,
            "anderson_m": 5,
            "auxiliary": serde_json::to_value(&cond_spec).unwrap()
        }),
    );

    // Project JSON round-trip check with CoupledLoopConfig
    let json_text = e.project().to_json();
    let rt_proj = Project::from_json(&json_text).expect("Project round-trip with coupled_loop");
    assert_eq!(&rt_proj, e.project());

    call(&mut e, "solve_run", json!({ "physics": "dfs" }));

    let iters = scalar_num(&e, "dfs", "coupled_iterations");
    let rel_res = scalar_num(&e, "dfs", "coupled_rel_residual");
    let p_joule = scalar_num(&e, "dfs", "P_joule (W)");
    let t_max = scalar_num(&e, "dfs", "T max (K)");
    assert!(
        iters <= 6.0,
        "electrothermal Anderson loop must converge in <= 6 iterations, got {iters} (rel_res = {rel_res:.3e})"
    );
    assert!(rel_res <= 1e-8, "coupled_rel_residual {rel_res:.3e} > 1e-8");
    assert!(t_max > 350.0 && t_max < 700.0, "unexpected self-heated T_max = {t_max} K");

    // Verify exact global thermal-electrical energy conservation:
    // Surface convective + Stefan-Boltzmann radiative heat rejection across `"outer"` equals `P_joule (W)`!
    let q_int = call(
        &mut e,
        "integrate",
        json!({
            "physics": "dfs",
            "boundary": "outer",
            "expr": "800.0 * (T - 300.0) + 0.8 * 5.670374419e-8 * (T^4 - 300.0^4)"
        }),
    );
    let p_out = q_int["value"].as_f64().expect("boundary convective + radiative heat rejection");
    let energy_err = (p_out - p_joule).abs() / p_joule;
    assert!(
        energy_err < 1e-4,
        "electrothermal power conservation error {energy_err:.3e} >= 1e-4: P_out = {p_out:.6e} W vs P_joule = {p_joule:.6e} W"
    );

    // Also verify Q_joule volume integral matches P_joule
    let q_vol_int = call(&mut e, "integrate", json!({ "physics": "dfs", "expr": "Q_joule" }));
    let p_vol = q_vol_int["value"].as_f64().unwrap();
    assert!(((p_vol - p_joule) / p_joule).abs() < 1e-8, "volume Q_joule integral {p_vol} vs {p_joule}");
}

#[test]
fn anderson_two_way_photothermal_driven_em_heat_loop_converges_and_conserves_energy() {
    // 3D lossy thermo-dielectric stub resonator [0, a] x [0, L] x [0, b] driven at f = 6.5 GHz via a lumped port:
    // - Relative permittivity eps_T(T) increases strongly with T (10.0 at 300 K -> 14.0 at 500 K), shifting
    //   the electrical length and S11 reflection (and thus absorbed RF/optical power P_loss(T))!
    // - Loss tangent tan_delta = 0.04 produces volumetric Q_opt ("W/m^3").
    // - Scaled by `source_scale` (incident drive power multiplier), `dfs::heat_stationary` heats the resonator
    //   against bottom heat-sink + top radiative cooling until (E, T) reaches self-consistency!
    let thermo_diel_toml = r#"
format = 2
name = "thermo_dielectric_res"
kind = "dielectric"
reference_temperature = "300 K"
density = "3200 kg/m^3"
source = { cite = "Synthetic verified FR-38 photothermal resonator test material", verified = "2026-11-27" }

[thermal]
conductivity = "5.0 W/(m*K)"
specific_heat = "700 J/(kg*K)"

[dielectric]
eps_free = 10.0
eps_T = { form = "tabulated_T", T = [250.0, 300.0, 500.0], values = [9.0, 10.0, 14.0] }
loss_tangent = 0.04
"#;

    let (a, l, b) = (6.0e-3, 8.0e-3, 3.0e-3);
    let freq = 6.5e9;
    let source_scale = 120.0;

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "photothermal-stub", "dim": 3, "h": 3.2e-3 }));
    call(&mut e, "add_box", json!({ "name": "cav", "min": [0.0, 0.0, 0.0], "max": [a, l, b] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "cav" }] }));
    call(&mut e, "add_material", json!({ "name": "thermo_dielectric_res", "toml": thermo_diel_toml }));

    // EM driven study at single frequency `freq`
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
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "thermo_dielectric_res" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "outer & plane(x = 0)", "kind": "pmc" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": format!("outer & plane(x = {a})"), "kind": "pmc" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "outer & plane(y = 0)",
            "kind": "lumped_port",
            "params": { "name": "p1", "Z_0": 50.0, "direction": [0.0, 0.0, 1.0] }
        }),
    );

    // Thermal stationary study on the same 3D mesh (P2 temperature with bottom convective heat-sink + top RadiationBC)
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "core", "material": "thermo_dielectric_res" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer",
            "kind": "convection",
            "params": { "h": 2500.0, "t_ext": 300.0 }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer",
            "kind": "radiation",
            "params": { "emissivity": 0.9, "t_ambient": 300.0, "view_factor": 1.0 }
        }),
    );

    // 1. Without `set_coupled_loop`, solving `["ems", "dfs"]` must reject the accidental cycle!
    let cycle_err = e.call("solve_run", &json!({ "physics": ["ems", "dfs"] })).unwrap_err().to_string();
    assert!(cycle_err.contains("cycle"), "expected cycle rejection without coupled_loop, got: {cycle_err}");

    // Solve cold EM alone first so we can verify hot thermo-optic detuning
    call(&mut e, "solve_run", json!({ "physics": "ems" }));
    let cold_em = e.result("ems").unwrap().0.clone();
    let cold_s11 = DrivenPortData::of(&cold_em).unwrap().s_matrix[0][0][0];

    // 2. Enable two-way Anderson coupled loop and solve `["ems", "dfs"]`!
    call(
        &mut e,
        "set_coupled_loop",
        json!({ "max_iters": 12, "rtol": 1e-7, "anderson_m": 5, "source_scale": source_scale }),
    );
    call(&mut e, "solve_run", json!({ "physics": ["ems", "dfs"] }));

    let iters = scalar_num(&e, "dfs", "coupled_iterations");
    let rel_res = scalar_num(&e, "dfs", "coupled_rel_residual");
    let p_loss_unit = scalar_num(&e, "ems", "P_loss (W)");
    let t_max = scalar_num(&e, "dfs", "T max (K)");
    assert!(
        iters <= 6.0,
        "photothermal Anderson loop must converge in <= 6 outer iterations, got {iters} (rel_res = {rel_res:.3e})"
    );
    assert!(rel_res <= 1e-7, "photothermal coupled_rel_residual {rel_res:.3e} > 1e-7");
    assert!(t_max > 305.0, "expected measurable photothermal heating above 300 K, got T_max = {t_max}");

    // Check EM Poynting power balance: P_loss_unit == 0.5 * (1 - |S11|^2)
    let hot_em = e.result("ems").unwrap().0.clone();
    let hot_s11 = DrivenPortData::of(&hot_em).unwrap().s_matrix[0][0][0];
    let s11_sq = hot_s11[0] * hot_s11[0] + hot_s11[1] * hot_s11[1];
    let p_in_em = 0.5 * (1.0 - s11_sq);
    let em_bal_err = (p_loss_unit - p_in_em).abs() / p_in_em.max(1e-30);
    assert!(
        em_bal_err < 1e-4,
        "EM Poynting power balance error {em_bal_err:.3e}: P_loss = {p_loss_unit:.6e} vs 0.5*(1-|S11|^2) = {p_in_em:.6e}"
    );

    // Check thermo-optic shift moved S11 measurably from cold S11
    let ds11 = (hot_s11[0] - cold_s11[0]).hypot(hot_s11[1] - cold_s11[1]);
    assert!(ds11 > 1e-3, "expected measurable thermo-optic S11 shift from self-heating, got |dS11| = {ds11:.3e}");

    // Check Q_opt integral on EM == P_loss_unit
    let q_opt_int = call(&mut e, "integrate", json!({ "physics": "ems", "expr": "Q_opt" }))["value"].as_f64().unwrap();
    assert!(((q_opt_int - p_loss_unit) / p_loss_unit).abs() < 1e-10);

    // Check global thermal energy conservation: total convective + radiative heat rejection across outer boundary == source_scale * P_loss_unit!
    let p_expected_thermal = source_scale * p_loss_unit;
    let p_thermal_out = call(
        &mut e,
        "integrate",
        json!({
            "physics": "dfs",
            "boundary": "outer",
            "expr": "2500.0 * (T - 300.0) + 0.9 * 5.670374419e-8 * (T^4 - 300.0^4)"
        }),
    )["value"]
        .as_f64()
        .unwrap();
    let thermal_bal_err = (p_thermal_out - p_expected_thermal).abs() / p_expected_thermal;
    assert!(
        thermal_bal_err < 1e-4,
        "photothermal boundary heat balance error {thermal_bal_err:.3e} >= 1e-4: P_out = {p_thermal_out:.6e} W vs scaled P_loss = {p_expected_thermal:.6e} W"
    );
}

#[test]
fn fr38_coupled_thermal_and_radiation_are_deterministic_across_nyx_threads() {
    const ENV_CHILD: &str = "NYX_FR38_DET_CHILD";
    if std::env::var_os(ENV_CHILD).is_some() {
        let mat_toml = r#"
format = 2
name = "det_mat"
kind = "dielectric"
reference_temperature = "300 K"
density = "2500 kg/m^3"
source = { cite = "Synthetic verified FR-38 determinism test material", verified = "2026-11-27" }

[thermal]
conductivity = { form = "tabulated_T", T = [250.0, 800.0], values = [15.0, 25.0], unit = "W/(m*K)" }

[electrical]
conductivity = { form = "tabulated_T", T = [250.0, 800.0], values = [1.5e5, 0.8e5], unit = "S/m" }
"#;
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "det", "dim": 2, "h": 2.5e-4 }));
        call(&mut e, "add_rect", json!({ "name": "bar", "min": [0.0, 0.0], "max": [1.5e-3, 0.5e-3] }));
        call(&mut e, "add_overlay", json!({ "regions": [{ "name": "r", "input": "bar" }] }));
        call(&mut e, "add_material", json!({ "name": "det_mat", "toml": mat_toml }));
        call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": "r", "material": "det_mat" }));
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
            "add_condition",
            json!({ "physics": "dfs", "on": "outer", "kind": "convection", "params": { "h": 600.0, "t_ext": 300.0 } }),
        );
        let mut cond = ProblemSpec::new("dfs", 2, "conduction");
        cond.materials.insert("r".into(), nyx_physics::MaterialAssign::Named("det_mat".into()));
        cond.bcs.push(nyx_physics::BcSpec {
            on: "outer & plane(x = 0)".into(),
            kind: "ground".into(),
            params: json!({}),
        });
        cond.bcs.push(nyx_physics::BcSpec {
            on: "outer & plane(x = 0.0015)".into(),
            kind: "voltage".into(),
            params: json!({ "value": 0.3 }),
        });
        call(
            &mut e,
            "set_coupled_loop",
            json!({ "max_iters": 10, "rtol": 1e-9, "anderson_m": 5, "auxiliary": serde_json::to_value(&cond).unwrap() }),
        );
        call(&mut e, "solve_run", json!({ "physics": "dfs" }));
        let t_vals = &e.result("dfs").unwrap().0.solutions[0].field("T").unwrap().values;
        let bits: Vec<u64> = t_vals.iter().map(|v| v.to_bits()).collect();
        println!("FR38_BITS={}", serde_json::to_string(&bits).unwrap());
        return;
    }

    let exe = std::env::current_exe().expect("test binary path");
    let run_threads = |nt: &str| -> Vec<u64> {
        let out = std::process::Command::new(&exe)
            .env(ENV_CHILD, "1")
            .env("NYX_THREADS", nt)
            .args(["--exact", "fr38_coupled_thermal_and_radiation_are_deterministic_across_nyx_threads", "--nocapture"])
            .output()
            .expect("spawn child");
        assert!(out.status.success(), "child NYX_THREADS={nt} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let line = stdout.lines().find_map(|l| l.strip_prefix("FR38_BITS=")).expect("FR38_BITS line");
        serde_json::from_str(line).unwrap()
    };

    let b1 = run_threads("1");
    let b8 = run_threads("8");
    assert_eq!(b1, b8, "0-ULP thread determinism violated between NYX_THREADS=1 and NYX_THREADS=8");
}
