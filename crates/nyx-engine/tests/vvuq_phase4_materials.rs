// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 VVUQ verification tests for active-study material provenance & temperature-range guards:
//! - Active-study `Source::flags_uncertainty()` (`material_approximate`, `material_unverified`) emitted ONLY for queried tables;
//! - Active-study `TempCurve::eval_checked` and `AxisModel::dn_dt_range` emitting `"temperature_clamped"` in `SolveResult::diagnostics.warnings`;
//! - Opt-in `strict_material_range` refusing out-of-range temperature evaluations via both `Project` and `solve_run`;
//! - Worker process (`isolate: true`) round-trip persistence of warnings and strict refusal.

use nyx_engine::{Engine, Project};
use nyx_physics::DiagnosticSeverity;
use serde_json::{json, Value};
use std::path::PathBuf;

fn worker_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"))
}

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn warnings_of(e: &Engine, physics: &str) -> Vec<(String, DiagnosticSeverity, String)> {
    e.result(physics)
        .expect("result exists")
        .0
        .diagnostics
        .as_ref()
        .map(|d| d.warnings.iter().map(|w| (w.code.clone(), w.severity, w.message.clone())).collect())
        .unwrap_or_default()
}

fn warning_codes(e: &Engine, physics: &str) -> Vec<String> {
    warnings_of(e, physics).into_iter().map(|(c, _, _)| c).collect()
}

fn setup_unit_square(physics: &str, study: &str, study_params: Value, order: u8, region_mat: &str) -> Engine {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p4-mat-test", "dim": 2, "h": 0.25 }));
    call(&mut e, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "body", "input": "box" }] }));
    call(
        &mut e,
        "add_physics",
        json!({ "physics": physics, "study": study, "study_params": study_params, "order": order }),
    );
    call(&mut e, "set_material", json!({ "physics": physics, "region": "body", "material": region_mat }));
    e
}

#[test]
fn source_uncertainty_warnings_are_scoped_to_actively_queried_tables() {
    // 1. Optical `ems` waveguide using `silicon` (verified [optical], APPROXIMATE [thermal] & [dielectric]):
    // must NOT emit `material_approximate` for unused thermal/dielectric tables.
    let mut e_opt = Engine::new();
    call(&mut e_opt, "project_new", json!({ "title": "opt-si", "dim": 2, "h": 1.2e-7 }));
    call(&mut e_opt, "add_rect", json!({ "name": "clad", "min": [-1.2e-6, -0.8e-6], "max": [1.2e-6, 0.8e-6] }));
    call(&mut e_opt, "add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 2.2e-7] }));
    call(
        &mut e_opt,
        "add_overlay",
        json!({ "regions": [{ "name": "air", "input": "clad" }, { "name": "si", "input": "core" }] }),
    );
    call(
        &mut e_opt,
        "add_physics",
        json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 1.934e14, "n_modes": 1 }, "order": 1 }),
    );
    call(&mut e_opt, "set_material", json!({ "physics": "ems", "region": "air", "material": "air" }));
    call(&mut e_opt, "set_material", json!({ "physics": "ems", "region": "si", "material": "silicon" }));
    call(&mut e_opt, "solve_run", json!({}));

    let opt_codes = warning_codes(&e_opt, "ems");
    assert!(
        !opt_codes.iter().any(|c| c == "material_approximate" || c == "material_unverified"),
        "optical ems on silicon/air must not warn about unused approximate thermal/dielectric tables: {opt_codes:?}"
    );

    // 2. Same `silicon` in `dfs` electrostatic DOES query [dielectric] (which is APPROXIMATE):
    let mut e_es = setup_unit_square("dfs", "electrostatic", json!({}), 1, "silicon");
    call(
        &mut e_es,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "ground", "params": {} }),
    );
    call(
        &mut e_es,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "potential", "params": { "value": 1.0 } }),
    );
    call(&mut e_es, "solve_run", json!({}));

    let es_warns = warnings_of(&e_es, "dfs");
    assert!(
        es_warns.iter().any(|(c, sev, msg)| c == "material_approximate"
            && *sev == DiagnosticSeverity::Warning
            && msg.contains("[dielectric]")),
        "expected material_approximate warning for silicon [dielectric], got: {es_warns:?}"
    );
    assert!(
        !es_warns.iter().any(|(c, _, _)| c == "material_unverified"),
        "unexpected material_unverified warning: {es_warns:?}"
    );

    // 3. `sms` eigenfrequency with verified `silicon` emits 0 uncertainty warnings, while `silicon_nitride` (`density` APPROXIMATE) emits `material_approximate`.
    // 3. Mechanical studies: `silicon` has verified [elastic] and density, plus APPROXIMATE [viscoelastic].
    // `sms::static` does NOT use [viscoelastic], so `silicon` produces 0 material_* warnings;
    // `sms::eigenfrequency` DOES load `eta` ([viscoelastic]) for Q damping, so it flags ONLY [viscoelastic].
    let mut e_sms_static = setup_unit_square("sms", "static", json!({}), 1, "silicon");
    call(
        &mut e_sms_static,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    call(&mut e_sms_static, "solve_run", json!({}));
    assert!(
        !warning_codes(&e_sms_static, "sms").iter().any(|c| c.starts_with("material_")),
        "verified silicon static elastic must not warn about unused viscoelastic/dielectric/thermal: {:?}",
        warnings_of(&e_sms_static, "sms")
    );

    let mut e_sms_si =
        setup_unit_square("sms", "eigenfrequency", json!({ "n_modes": 1, "f_search": 1000.0 }), 1, "silicon");
    call(
        &mut e_sms_si,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    call(&mut e_sms_si, "solve_run", json!({}));
    let si_eig_warns = warnings_of(&e_sms_si, "sms");
    assert!(
        si_eig_warns.iter().any(|(c, _, msg)| c == "material_approximate" && msg.contains("[viscoelastic]")),
        "sms eigenfrequency uses silicon [viscoelastic] eta for Q: {si_eig_warns:?}"
    );
    assert!(
        !si_eig_warns
            .iter()
            .any(|(c, _, msg)| c == "material_approximate" && (msg.contains("[elastic]") || msg.contains("[density]"))),
        "verified silicon [elastic] and [density] must not warn: {si_eig_warns:?}"
    );

    let mut e_sms_sin =
        setup_unit_square("sms", "eigenfrequency", json!({ "n_modes": 1, "f_search": 1000.0 }), 1, "silicon_nitride");
    call(
        &mut e_sms_sin,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    call(&mut e_sms_sin, "solve_run", json!({}));
    assert!(
        warning_codes(&e_sms_sin, "sms").iter().any(|c| c == "material_approximate"),
        "expected material_approximate for silicon_nitride density/mechanics, got: {:?}",
        warnings_of(&e_sms_sin, "sms")
    );

    // 4. Superconducting BC (`surface_impedance`) on `aluminium` (APPROXIMATE sigma_n) flags `material_approximate`, while `niobium` does not.
    let run_sc_cavity = |sc_mat: &str| {
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p4-sc-cavity", "dim": 3, "h": 0.5 }));
        call(&mut e, "add_box", json!({ "name": "box", "min": [0.0, 0.0, 0.0], "max": [1.0, 1.0, 1.0] }));
        call(&mut e, "add_overlay", json!({ "regions": [{ "name": "body", "input": "box" }] }));
        call(
            &mut e,
            "add_physics",
            json!({ "physics": "ems", "study": "cavity", "study_params": { "n_modes": 1, "f_search": 2.1e8 }, "order": 1 }),
        );
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "body", "material": "air" }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "outer & !plane(z = 0)", "kind": "pec", "params": {} }),
        );
        call(
            &mut e,
            "add_condition",
            json!({
                "physics": "ems",
                "on": "outer & plane(z = 0)",
                "kind": "surface_impedance",
                "params": { "material": sc_mat, "thickness": 100.0e-9, "temperature": 0.2, "lossless": true }
            }),
        );
        call(&mut e, "solve_run", json!({}));
        e
    };
    let e_sc_al = run_sc_cavity("aluminium");
    assert!(
        warnings_of(&e_sc_al, "ems").iter().any(|(c, sev, msg)| c == "material_approximate"
            && *sev == DiagnosticSeverity::Warning
            && msg.contains("aluminium")
            && msg.contains("[superconducting]")),
        "expected surface_impedance aluminium warning: {:?}",
        warnings_of(&e_sc_al, "ems")
    );
    let e_sc_nb = run_sc_cavity("niobium");
    assert!(
        !warning_codes(&e_sc_nb, "ems").iter().any(|c| c.starts_with("material_")),
        "verified niobium [superconducting] must not warn: {:?}",
        warnings_of(&e_sc_nb, "ems")
    );

    // 5. Unverified custom material in active `[thermal]` table emits `material_unverified`.
    let mut e_unver = setup_unit_square("dfs", "heat_stationary", json!({}), 1, "unver_mat");
    call(
        &mut e_unver,
        "add_material",
        json!({
            "name": "unver_mat",
            "toml": r#"
format = 2
name = "unver_mat"
source = { cite = "Verified wrapper", verified = "2026-11-26" }
[thermal]
conductivity = "12.0 W/(m*K)"
source = { cite = "UNVERIFIED synthetic test thermal value" }
"#
        }),
    );
    call(
        &mut e_unver,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": 290.0 } }),
    );
    call(
        &mut e_unver,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "temperature", "params": { "value": 300.0 } }),
    );
    call(&mut e_unver, "solve_run", json!({}));

    let unver_warns = warnings_of(&e_unver, "dfs");
    assert!(
        unver_warns.iter().any(|(c, sev, _)| c == "material_unverified" && *sev == DiagnosticSeverity::Warning),
        "expected material_unverified warning, got: {unver_warns:?}"
    );
}

#[test]
fn thermal_solve_temperature_clamping_warning_and_strict_material_range_refusal() {
    // 1. In-range copper_4k solve (T in [10, 200] K within [4, 293.15] K) passes even with strict_material_range = true.
    let mut e_ok = setup_unit_square("dfs", "heat_stationary", json!({}), 1, "copper_4k");
    call(
        &mut e_ok,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": 10.0 } }),
    );
    call(
        &mut e_ok,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "temperature", "params": { "value": 200.0 } }),
    );
    call(&mut e_ok, "solve_run", json!({ "strict_material_range": true }));
    let ok_codes = warning_codes(&e_ok, "dfs");
    assert!(!ok_codes.iter().any(|c| c == "temperature_clamped"), "in-range solve must not warn: {ok_codes:?}");

    // 2. Out-of-range solve (420 K > 293.15 K) warns by default (`temperature_clamped`) without failing.
    let mut e_clamp = setup_unit_square("dfs", "heat_stationary", json!({}), 1, "copper_4k");
    call(
        &mut e_clamp,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": 10.0 } }),
    );
    call(
        &mut e_clamp,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "temperature", "params": { "value": 420.0 } }),
    );
    call(&mut e_clamp, "solve_run", json!({}));

    let clamp_warns = warnings_of(&e_clamp, "dfs");
    assert!(
        clamp_warns.iter().any(|(c, sev, msg)| c == "temperature_clamped"
            && *sev == DiagnosticSeverity::Warning
            && msg.contains("copper_4k")),
        "expected temperature_clamped warning, got: {clamp_warns:?}"
    );
    assert!(
        clamp_warns.iter().any(|(c, _, _)| c == "material_approximate"),
        "expected material_approximate warning alongside clamp, got: {clamp_warns:?}"
    );

    // `result_info` exposes these warnings cleanly through the engine JSON surface.
    let info = call(&mut e_clamp, "result_info", json!({ "physics": "dfs" }));
    let info_codes: Vec<&str> =
        info[0]["diagnostics"]["warnings"].as_array().unwrap().iter().filter_map(|w| w["code"].as_str()).collect();
    assert!(info_codes.contains(&"temperature_clamped"), "{info}");

    // 3. Opt-in strict_material_range via `solve_run` option refuses with Err.
    let err_opt = e_clamp.call("solve_run", &json!({ "strict_material_range": true })).unwrap_err().to_string();
    assert!(err_opt.contains("strict_material_range") && err_opt.contains("copper_4k"), "{err_opt}");

    // 4. Persistent `Project::strict_material_range` round-trips through JSON and refuses solves.
    let mut proj = e_clamp.project().clone();
    proj.strict_material_range = true;
    let roundtrip = Project::from_json(&proj.to_json()).unwrap();
    assert!(roundtrip.strict_material_range);
    let mut e_proj = Engine::new();
    e_proj.set_project(roundtrip);
    let err_proj = e_proj.call("solve_run", &json!({})).unwrap_err().to_string();
    assert!(err_proj.contains("strict_material_range"), "{err_proj}");

    // 5. Worker-isolated solve (`isolate: true`) persists warnings AND respects strict_material_range refusal.
    let mut e_worker = setup_unit_square("dfs", "heat_stationary", json!({}), 1, "copper_4k");
    e_worker.set_worker(worker_bin(), false);
    call(
        &mut e_worker,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": 10.0 } }),
    );
    call(
        &mut e_worker,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "temperature", "params": { "value": 420.0 } }),
    );
    call(&mut e_worker, "solve_run", json!({ "isolate": true }));
    let w_codes = warning_codes(&e_worker, "dfs");
    assert!(w_codes.iter().any(|c| c == "temperature_clamped"), "worker solve lost temperature_clamped: {w_codes:?}");

    let err_worker =
        e_worker.call("solve_run", &json!({ "isolate": true, "strict_material_range": true })).unwrap_err().to_string();
    assert!(err_worker.contains("strict_material_range"), "{err_worker}");
}

#[test]
fn sms_thermal_expansion_and_coupled_optical_dn_dt_range_are_checked_only_when_active() {
    let custom_mat = r#"
format = 2
name = "bounded_coupled"
source = { cite = "Verified unit test material", verified = "2026-11-26" }
density = "2500 kg/m^3"

[elastic]
unit = "GPa"
E = 100.0
nu = 0.25
source = { cite = "Verified unit test elastic", verified = "2026-11-26" }

[thermal]
conductivity = "50.0 W/(m*K)"
specific_heat = "500.0 J/(kg*K)"
expansion_T = { form = "log10_poly", T_range = [100.0, 400.0], coeffs = [-5.5, 0.2], unit = "1/K" }
source = { cite = "Verified unit test thermal", verified = "2026-11-26" }

[optical]
wavelength_unit = "um"
reference_temperature = "293.15 K"
model = "constant"
n = 2.0
dn_dT = "1.0e-4 1/K"
dn_dT_range = [200.0, 400.0]
source = { cite = "Verified unit test optical", verified = "2026-11-26" }
"#;

    // 1. Sms static with thermal="none": expansion_T is unused -> 0 warnings even with strict_material_range.
    let mut e_sms = setup_unit_square("sms", "static", json!({}), 1, "bounded_coupled");
    call(&mut e_sms, "add_material", json!({ "name": "bounded_coupled", "toml": custom_mat }));
    call(
        &mut e_sms,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    call(&mut e_sms, "solve_run", json!({ "strict_material_range": true }));
    assert!(warning_codes(&e_sms, "sms").is_empty());

    // 2. Sms static with thermal="uniform" and temperature = 50 K (below [100, 400] K) warns & strict-fails.
    call(
        &mut e_sms,
        "set_physics_params",
        json!({ "physics": "sms", "params": { "plane": "stress", "thermal": "uniform", "temperature": 50.0, "t_ref": 293.15 } }),
    );
    call(&mut e_sms, "solve_run", json!({}));
    let sms_codes = warning_codes(&e_sms, "sms");
    assert!(sms_codes.iter().any(|c| c == "temperature_clamped"), "{sms_codes:?}");

    let err_sms = e_sms.call("solve_run", &json!({ "strict_material_range": true })).unwrap_err().to_string();
    assert!(err_sms.contains("strict_material_range") && err_sms.contains("expansion_T"), "{err_sms}");

    // 3. Coupled dfs heat (500 K) -> ems optical (1.55 um, dn_dT_range = [200, 400] K): ems warns and strict-refuses.
    let mut e_chain = Engine::new();
    call(&mut e_chain, "project_new", json!({ "title": "thermo-optic-guard", "dim": 2, "h": 2.5e-7 }));
    call(&mut e_chain, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0e-6, 1.0e-6] }));
    call(&mut e_chain, "add_overlay", json!({ "regions": [{ "name": "body", "input": "box" }] }));
    call(&mut e_chain, "add_material", json!({ "name": "bounded_coupled", "toml": custom_mat }));
    call(&mut e_chain, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 1 }));
    call(&mut e_chain, "set_material", json!({ "physics": "dfs", "region": "body", "material": "bounded_coupled" }));
    call(
        &mut e_chain,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer", "kind": "temperature", "params": { "value": 500.0 } }),
    );
    call(
        &mut e_chain,
        "add_physics",
        json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 1.934e14, "n_modes": 1 }, "order": 1 }),
    );
    call(&mut e_chain, "set_material", json!({ "physics": "ems", "region": "body", "material": "bounded_coupled" }));
    // Without explicit opt-in `set_coupled_loop`, bidirectional `dfs <-> ems` (`temperature` + `heat_source`) is rejected as a cycle.
    let err_cycle = e_chain.call("solve_run", &json!({ "physics": ["dfs", "ems"] })).unwrap_err().to_string();
    assert!(err_cycle.contains("cycle"), "expected unconfigured cycle rejection, got: {err_cycle}");
    call(&mut e_chain, "set_coupled_loop", json!({ "max_iters": 10, "rtol": 1e-6, "anderson_m": 5 }));
    call(&mut e_chain, "solve_run", json!({ "physics": ["dfs", "ems"] }));

    let ems_warns = warnings_of(&e_chain, "ems");
    assert!(
        ems_warns.iter().any(|(c, _, msg)| c == "temperature_clamped" && msg.contains("optical")),
        "expected optical dn_dT_range clamp warning from coupled heat temperature: {ems_warns:?}"
    );

    let err_ems = e_chain
        .call("solve_run", &json!({ "physics": ["dfs", "ems"], "strict_material_range": true }))
        .unwrap_err()
        .to_string();
    assert!(err_ems.contains("strict_material_range") && err_ems.contains("optical"), "{err_ems}");
}
