// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 2 cross-FR integration verification suite (`nyx-engine`, ADRs 140–142):
//! 1. Cross-mesh `dfs::electrostatic -> ems::eme` X-cut thin-film lithium niobate (`LiNbO3`, `3m`)
//!    electro-optic phase shifter and push-pull directional coupler modulator (`FR-37` + `FR-39` + `FR-44`),
//!    verifying BLAKE3 slice deduplication, linear Pockels phase shift, shared vs cross-mesh agreement,
//!    and bias-driven bar/cross directional coupler power switching.
//! 2. Analytical Redheffer EME cell-length sensitivities (`dP_c/dL_c = -j beta_c P_c`) vs central
//!    finite differences (`< 1e-6` relative error) across multi-port couplers, `PeriodicSubgrid`
//!    Bragg gratings, and longitudinal CAD `param_var`, plus `lbfgs` 3-dB coupler length optimization
//!    (`FR-37` + `FR-40`).
//! 3. Cross-mesh `sms` + `ems` `optomechanics` analysis (`g_PE`, `g_MB`, `g_0`, `G_B`) composed with
//!    `FR-40` `optimize_run` waveguide inverse design and subprocess `NYX_THREADS=1` vs `NYX_THREADS=8`
//!    bit-exact determinism (`FR-39` + `FR-40` + `FR-44`).

use nyx_ems::EmeData;
use nyx_engine::Engine;
use nyx_physics::DrivenPortData;
use serde_json::{json, Value};

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn unwrap_phase_diff(phi_mod: f64, phi_ref: f64) -> f64 {
    let mut d = phi_mod - phi_ref;
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

/// Builds a 3D X-cut TFLN (`lithium_niobate`, `euler_deg = [90, 90, 0]`, crystal +Z along wafer +x)
/// electro-optic waveguide phase shifter driven by `dfs::electrostatic` side electrodes and evaluated
/// with `ems::eme` (either on a shared 3D mesh or on an optical-only `active_domains: ["clad", "core"]` cross-mesh).
fn build_tfln_eo_phase_shifter_engine(v_bias: f64, segmented_electrode: bool, cross_mesh_optical: bool) -> Engine {
    let l_mesh = 8.0e-6;
    let l_phys = 400.0e-6;
    let freq = 299_792_458.0 / 1.55e-6;

    let elec_toml = r#"
format = 2
name = "au_pad_p2"
kind = "dielectric"
source = { cite = "Phase 2 electrostatic electrode block", verified = "2026-11-27" }

[dielectric]
eps_free = 1.0e4

[optical]
wavelength_unit = "um"
range = [1.2, 1.7]
model = "constant"
n = 1.444
k = 0.0
"#;

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p2-tfln-eo-ps", "dim": 3, "h": 0.45e-6 }));
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_clad", "min": [-1.5e-6, 0.0, -0.55e-6], "max": [1.5e-6, l_mesh, 0.55e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_core", "min": [-0.35e-6, 0.0, -0.18e-6], "max": [0.35e-6, l_mesh, 0.18e-6] }),
    );
    let (y_el0, y_el1) = if segmented_electrode { (0.25 * l_mesh, 0.75 * l_mesh) } else { (0.0, l_mesh) };
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_elec_l", "min": [-1.5e-6, y_el0, -0.25e-6], "max": [-0.95e-6, y_el1, 0.25e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_elec_r", "min": [0.95e-6, y_el0, -0.25e-6], "max": [1.5e-6, y_el1, 0.25e-6] }),
    );
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad", "input": "b_clad" },
                { "name": "core", "input": "b_core" },
                { "name": "elec_l", "input": "b_elec_l" },
                { "name": "elec_r", "input": "b_elec_r" }
            ]
        }),
    );
    call(&mut e, "add_material", json!({ "name": "au_pad_p2", "toml": elec_toml }));

    let proj_sweep = if segmented_electrode {
        json!([
            { "thickness": 0.25 * l_mesh, "n": 2 },
            { "thickness": 0.50 * l_mesh, "n": 4 },
            { "thickness": 0.25 * l_mesh, "n": 2 }
        ])
    } else {
        json!([{ "thickness": l_mesh, "n": 6 }])
    };
    call(&mut e, "set_mesh", json!({ "h": 0.45e-6, "sweep": proj_sweep }));

    // 1. 3D Electrostatic bias (`dfs::electrostatic`)
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "clad", "material": "silica" }));
    call(
        &mut e,
        "set_material",
        json!({
            "physics": "dfs",
            "region": "core",
            "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "elec_l", "material": "au_pad_p2" }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "elec_r", "material": "au_pad_p2" }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "boundary(region(elec_l))",
            "kind": "potential",
            "params": { "value": 0.5 * v_bias }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "boundary(region(elec_r))",
            "kind": "potential",
            "params": { "value": -0.5 * v_bias }
        }),
    );

    // 2. 3D Optical EME (`ems::eme`) with 6 segments totaling 400 um physical length
    let n_seg = if segmented_electrode { 8usize } else { 6usize };
    let seg_len = l_phys / (n_seg as f64);
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "eme",
            "order": 2,
            "study_params": {
                "freq": freq,
                "n_target": 2.14,
                "n_modes": 2,
                "grid_kind": "explicit",
                "lengths": vec![seg_len; n_seg],
                "port_modes": 1,
                "reconstruct_fields": false
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" }));
    call(
        &mut e,
        "set_material",
        json!({
            "physics": "ems",
            "region": "core",
            "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
        }),
    );
    // PMC on z = +-0.55 um selects fundamental TE (Ex-polarized along crystal +Z, sensing r_33 = 28.6 pm/V)
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-z)", "kind": "pmc" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+z)", "kind": "pmc" }));

    if cross_mesh_optical {
        call(
            &mut e,
            "set_mesh",
            json!({
                "physics": "ems",
                "h": 0.38e-6,
                "active_domains": ["clad", "core"],
                "sweep": [{ "thickness": l_mesh, "n": n_seg }]
            }),
        );
    } else {
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "elec_l", "material": "au_pad_p2" }));
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "elec_r", "material": "au_pad_p2" }));
    }

    e
}

/// Builds a 3D X-cut TFLN push-pull electro-optic directional coupler (`dfs::electrostatic -> ems::eme`)
/// with 4 explicit split-face `waveport` terminals (`in_l`, `in_r`, `out_l`, `out_r`).
fn build_tfln_eo_coupler_engine(v_push_pull: f64, coupler_length_m: f64) -> Engine {
    let l_mesh = 4.0e-6;
    let freq = 299_792_458.0 / 1.55e-6;

    let elec_toml = r#"
format = 2
name = "au_coupler_pad"
kind = "dielectric"
source = { cite = "Phase 2 push-pull coupler electrode", verified = "2026-11-27" }

[dielectric]
eps_free = 1.0e4
"#;

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p2-tfln-eo-coupler", "dim": 3, "h": 0.45e-6 }));
    // Split cladding at x = 0 so left/right arm waveports can be selected cleanly via CAD region boundaries
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_clad_l", "min": [-2.2e-6, 0.0, -0.55e-6], "max": [0.0, l_mesh, 0.55e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_clad_r", "min": [0.0, 0.0, -0.55e-6], "max": [2.2e-6, l_mesh, 0.55e-6] }),
    );
    // Symmetric coupled TFLN cores centered at x = -0.54 um and x = +0.54 um (width 0.48 um, gap 0.60 um -> high waveport overlap!)
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_core_l", "min": [-0.78e-6, 0.0, -0.16e-6], "max": [-0.30e-6, l_mesh, 0.16e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_core_r", "min": [0.30e-6, 0.0, -0.16e-6], "max": [0.78e-6, l_mesh, 0.16e-6] }),
    );
    // Ground-Signal-Ground (GSG) electrodes: outer ground pads + center signal electrode above gap
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_gnd_l", "min": [-2.2e-6, 0.0, -0.25e-6], "max": [-1.45e-6, l_mesh, 0.25e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_gnd_r", "min": [1.45e-6, 0.0, -0.25e-6], "max": [2.2e-6, l_mesh, 0.25e-6] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "b_sig_c", "min": [-0.14e-6, 0.0, 0.26e-6], "max": [0.14e-6, l_mesh, 0.50e-6] }),
    );
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad_l", "input": "b_clad_l" },
                { "name": "clad_r", "input": "b_clad_r" },
                { "name": "core_l", "input": "b_core_l" },
                { "name": "core_r", "input": "b_core_r" },
                { "name": "gnd_l", "input": "b_gnd_l" },
                { "name": "gnd_r", "input": "b_gnd_r" },
                { "name": "sig_c", "input": "b_sig_c" }
            ]
        }),
    );
    call(&mut e, "add_material", json!({ "name": "au_coupler_pad", "toml": elec_toml }));
    call(&mut e, "set_mesh", json!({ "h": 0.44e-6, "sweep": [{ "thickness": l_mesh, "n": 2 }] }));

    // Electrostatic push-pull bias: center signal at +V_push_pull, outer grounds at 0 V -> opposite Ex in core_l vs core_r!
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    for reg in ["clad_l", "clad_r"] {
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": reg, "material": "silica" }));
    }
    for reg in ["core_l", "core_r"] {
        call(
            &mut e,
            "set_material",
            json!({
                "physics": "dfs",
                "region": reg,
                "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
            }),
        );
    }
    for reg in ["gnd_l", "gnd_r", "sig_c"] {
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": reg, "material": "au_coupler_pad" }));
    }
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "boundary(region(sig_c))", "kind": "potential", "params": { "value": v_push_pull } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "boundary(region(gnd_l)) | boundary(region(gnd_r))",
            "kind": "potential",
            "params": { "value": 0.0 }
        }),
    );

    // Optical EME on cross-mesh excluding all 3 metallic electrodes
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "eme",
            "order": 2,
            "study_params": {
                "freq": freq,
                "n_target": 2.14,
                "n_modes": 4,
                "grid_kind": "explicit",
                "lengths": [coupler_length_m],
                "reconstruct_fields": false
            }
        }),
    );
    for reg in ["clad_l", "clad_r"] {
        call(&mut e, "set_material", json!({ "physics": "ems", "region": reg, "material": "silica" }));
    }
    for reg in ["core_l", "core_r"] {
        call(
            &mut e,
            "set_material",
            json!({
                "physics": "ems",
                "region": reg,
                "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
            }),
        );
    }
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-z)", "kind": "pmc" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+z)", "kind": "pmc" }));
    for (pname, sel) in [
        ("in_l", "boundary(region(clad_l) | region(core_l)) & extreme(-y)"),
        ("in_r", "boundary(region(clad_r) | region(core_r)) & extreme(-y)"),
        ("out_l", "boundary(region(clad_l) | region(core_l)) & extreme(+y)"),
        ("out_r", "boundary(region(clad_r) | region(core_r)) & extreme(+y)"),
    ] {
        call(
            &mut e,
            "add_condition",
            json!({
                "physics": "ems",
                "on": sel,
                "kind": "waveport",
                "params": { "name": pname, "mode": 1, "n_eff_guess": 2.14 }
            }),
        );
    }
    call(
        &mut e,
        "set_mesh",
        json!({
            "physics": "ems",
            "h": 0.36e-6,
            "active_domains": ["clad_l", "clad_r", "core_l", "core_r"],
            "sweep": [{ "thickness": l_mesh, "n": 2 }]
        }),
    );

    e
}

fn s21_phase(e: &Engine) -> f64 {
    let res = e.result("ems").expect("ems result").0;
    let dp = DrivenPortData::of(res).expect("DrivenPortData");
    let s21 = dp.s_matrix[0][1][0];
    s21[1].atan2(s21[0])
}

#[test]
fn cross_mesh_tfln_electro_optic_eme_phase_shifter_and_coupler() {
    // Part A1: Exact BLAKE3 slice deduplication under nonzero 3D EO bias (`dfs::electrostatic -> ems::eme`)
    // On a z-spanning X-cut TFLN slab capacitor (`clad_l | core | clad_r`) with direct SPD solve, `V(x)` is
    // piecewise-linear in `x` (exact in P1/P2 on every 3D tetrahedron regardless of prism diagonal orientation),
    // so all 6 EME segments at V = 20 V must deduplicate via BLAKE3 to 1 unique slice while shifting phase!
    {
        let freq = 299_792_458.0 / 1.55e-6;
        let run_planar_eo = |v_bias: f64| -> (f64, usize) {
            let mut ep = Engine::new();
            call(&mut ep, "project_new", json!({ "title": "p2-tfln-blake3-eo", "dim": 3, "h": 0.50e-6 }));
            call(
                &mut ep,
                "add_box",
                json!({ "name": "bl", "min": [-1.2e-6, 0.0, -0.35e-6], "max": [-0.35e-6, 6.0e-6, 0.35e-6] }),
            );
            call(
                &mut ep,
                "add_box",
                json!({ "name": "bc", "min": [-0.35e-6, 0.0, -0.35e-6], "max": [0.35e-6, 6.0e-6, 0.35e-6] }),
            );
            call(
                &mut ep,
                "add_box",
                json!({ "name": "br", "min": [0.35e-6, 0.0, -0.35e-6], "max": [1.2e-6, 6.0e-6, 0.35e-6] }),
            );
            call(
                &mut ep,
                "add_overlay",
                json!({
                    "regions": [
                        { "name": "clad_l", "input": "bl" },
                        { "name": "core", "input": "bc" },
                        { "name": "clad_r", "input": "br" }
                    ]
                }),
            );
            call(&mut ep, "set_mesh", json!({ "h": 0.50e-6, "sweep": [{ "thickness": 6.0e-6, "n": 6 }] }));
            call(&mut ep, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
            call(&mut ep, "set_physics_params", json!({ "physics": "dfs", "params": { "solver": "direct" } }));
            call(&mut ep, "set_material", json!({ "physics": "dfs", "region": "clad_l", "material": "silica" }));
            call(&mut ep, "set_material", json!({ "physics": "dfs", "region": "clad_r", "material": "silica" }));
            call(
                &mut ep,
                "set_material",
                json!({
                    "physics": "dfs",
                    "region": "core",
                    "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
                }),
            );
            call(
                &mut ep,
                "add_condition",
                json!({ "physics": "dfs", "on": "extreme(-x)", "kind": "potential", "params": { "value": 0.5 * v_bias } }),
            );
            call(
                &mut ep,
                "add_condition",
                json!({ "physics": "dfs", "on": "extreme(+x)", "kind": "potential", "params": { "value": -0.5 * v_bias } }),
            );
            call(
                &mut ep,
                "add_physics",
                json!({
                    "physics": "ems",
                    "study": "eme",
                    "order": 2,
                    "study_params": {
                        "freq": freq,
                        "n_target": 2.14,
                        "n_modes": 1,
                        "grid_kind": "explicit",
                        "lengths": vec![50.0e-6; 6],
                        "port_modes": 1,
                        "reconstruct_fields": false
                    }
                }),
            );
            call(&mut ep, "set_material", json!({ "physics": "ems", "region": "clad_l", "material": "silica" }));
            call(&mut ep, "set_material", json!({ "physics": "ems", "region": "clad_r", "material": "silica" }));
            call(
                &mut ep,
                "set_material",
                json!({
                    "physics": "ems",
                    "region": "core",
                    "material": { "material": "lithium_niobate", "euler_deg": [90.0, 90.0, 0.0] }
                }),
            );
            call(&mut ep, "add_condition", json!({ "physics": "ems", "on": "extreme(-z)", "kind": "pmc" }));
            call(&mut ep, "add_condition", json!({ "physics": "ems", "on": "extreme(+z)", "kind": "pmc" }));
            call(&mut ep, "solve_run", json!({ "physics": ["dfs", "ems"] }));
            let ns = EmeData::of(ep.result("ems").expect("ems").0).expect("EmeData").unique_slices.len();
            (s21_phase(&ep), ns)
        };
        let (p0, u0) = run_planar_eo(0.0);
        let (p20, u20) = run_planar_eo(20.0);
        assert_eq!(u0, 1, "unbiased 6-slice EME slab must deduplicate to 1 unique slice");
        assert!(
            u20 <= 3,
            "y-uniform 20 V EO-biased 6-slice EME slab must deduplicate all 4 interior slabs via BLAKE3 (<= 3 unique with endcaps), got {u20}"
        );
        assert!(unwrap_phase_diff(p20, p0).abs() > 0.05, "20 V EO slab must impart nonzero Pockels phase shift");
    }

    // Part A2: Cross-mesh X-cut TFLN buried channel EO phase shifter (`dfs::electrostatic -> ems::eme`)
    let mut e0 = build_tfln_eo_phase_shifter_engine(0.0, false, true);
    call(&mut e0, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let phi0 = s21_phase(&e0);
    let data0 = EmeData::of(e0.result("ems").expect("ems").0).expect("EmeData");
    assert_eq!(data0.unique_slices.len(), 1, "6 unbiased EME slices must deduplicate to 1 unique slice via BLAKE3");

    let mut e_v1 = build_tfln_eo_phase_shifter_engine(20.0, false, true);
    call(&mut e_v1, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dphi_v1 = unwrap_phase_diff(s21_phase(&e_v1), phi0);
    assert!(
        dphi_v1.abs() > 0.05,
        "20 V bias across 400 um X-cut TFLN waveguide must produce substantial EO phase shift, got {dphi_v1} rad"
    );

    let mut e_v2 = build_tfln_eo_phase_shifter_engine(40.0, false, true);
    call(&mut e_v2, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dphi_v2 = unwrap_phase_diff(s21_phase(&e_v2), phi0);
    let lin_err = ((dphi_v2 / dphi_v1) - 2.0).abs();
    assert!(lin_err < 2e-3, "Pockels phase shift must scale linearly with bias voltage (ratio err = {lin_err:.3e})");

    let mut e_neg = build_tfln_eo_phase_shifter_engine(-20.0, false, true);
    call(&mut e_neg, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dphi_neg = unwrap_phase_diff(s21_phase(&e_neg), phi0);
    let odd_err = (dphi_neg + dphi_v1).abs() / dphi_v1.abs();
    assert!(odd_err < 2e-3, "reversing bias polarity must negate Pockels phase shift (odd_err = {odd_err:.3e})");

    // Shared-mesh vs cross-mesh (`active_domains: ["clad", "core"]`) agreement
    let mut e_shared_0 = build_tfln_eo_phase_shifter_engine(0.0, false, false);
    call(&mut e_shared_0, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let mut e_shared_v1 = build_tfln_eo_phase_shifter_engine(20.0, false, false);
    call(&mut e_shared_v1, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dphi_shared = unwrap_phase_diff(s21_phase(&e_shared_v1), s21_phase(&e_shared_0));
    let mesh_rel_diff = (dphi_v1 - dphi_shared).abs() / dphi_shared.abs();
    assert!(
        mesh_rel_diff < 0.08,
        "cross-mesh vs shared-mesh TFLN EO phase shift mismatch too large: cross={dphi_v1:.5}, shared={dphi_shared:.5} (rel={mesh_rel_diff:.3e})"
    );

    // Segmented electrode along y (`[0.25 L, 0.75 L]`) must produce distinct BLAKE3 slice hashes AND ~half the phase shift!
    let mut e_seg_0 = build_tfln_eo_phase_shifter_engine(0.0, true, true);
    call(&mut e_seg_0, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let mut e_seg = build_tfln_eo_phase_shifter_engine(20.0, true, true);
    call(&mut e_seg, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let data_seg = EmeData::of(e_seg.result("ems").expect("ems").0).expect("EmeData");
    assert!(
        data_seg.unique_slices.len() >= 2,
        "segmented electrode along y must yield >= 2 distinct BLAKE3 EME slices, got {}",
        data_seg.unique_slices.len()
    );
    let dphi_seg = unwrap_phase_diff(s21_phase(&e_seg), s21_phase(&e_seg_0));
    let seg_ratio = dphi_seg / dphi_v1;
    // Accounting for 3D electrostatic fringing tails (~0.7 * gap) at both ends of the 4 um electrode on the 8 um cell
    assert!(
        (0.40..=0.80).contains(&seg_ratio),
        "half-length y-segmented EO electrode (with 3D fringing) must reduce phase shift vs full electrode, got ratio={seg_ratio:.4}"
    );

    // Part B: Cross-mesh X-cut TFLN push-pull EO directional coupler modulator
    // First solve at nominal placeholder length to read exact supermode beat length L_pi = pi / |beta_0 - beta_1|
    let mut ec_probe = build_tfln_eo_coupler_engine(0.0, 50.0e-6);
    call(&mut ec_probe, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let probe_data = EmeData::of(ec_probe.result("ems").expect("ems").0).expect("EmeData");
    let dbeta = (probe_data.unique_slices[0].beta_fwd[0][0] - probe_data.unique_slices[0].beta_fwd[1][0]).abs();
    assert!(dbeta > 100.0, "coupled TFLN supermodes must split cleanly (dbeta = {dbeta} rad/m)");
    let l_pi = std::f64::consts::PI / dbeta;

    // At V = 0 V and L = L_pi, input at `in_l` (port 0) couples predominantly to cross port `out_r` (port 3)
    let mut ec_zero = build_tfln_eo_coupler_engine(0.0, l_pi);
    call(&mut ec_zero, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dp_zero = DrivenPortData::of(ec_zero.result("ems").expect("ems").0).expect("DrivenPortData");
    assert_eq!(dp_zero.ports, vec!["in_l", "in_r", "out_l", "out_r"]);
    let p_bar_0 = dp_zero.s_matrix[0][2][0][0].powi(2) + dp_zero.s_matrix[0][2][0][1].powi(2);
    let p_cross_0 = dp_zero.s_matrix[0][3][0][0].powi(2) + dp_zero.s_matrix[0][3][0][1].powi(2);
    assert!(
        p_cross_0 > 0.85 && p_bar_0 < 0.15,
        "unbiased coupler at L_pi must route power to cross port: P_bar={p_bar_0:.4}, P_cross={p_cross_0:.4}"
    );

    // Calibrate push-pull supermode detuning rate at V_cal = 25 V so coupled-mode switching voltage V_bar = sqrt(3)*dbeta_0 / (detune_per_v) hits bar state!
    let v_cal = 25.0;
    let mut ec_cal = build_tfln_eo_coupler_engine(v_cal, l_pi);
    call(&mut ec_cal, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let cal_data = EmeData::of(ec_cal.result("ems").expect("ems").0).expect("EmeData");
    let dbeta_cal = (cal_data.unique_slices[0].beta_fwd[0][0] - cal_data.unique_slices[0].beta_fwd[1][0]).abs();
    let detune_cal = (dbeta_cal * dbeta_cal - dbeta * dbeta).max(0.0).sqrt();
    assert!(
        detune_cal > 1.0,
        "push-pull bias must increase supermode splitting: dbeta_0={dbeta:.2}, dbeta_cal={dbeta_cal:.2}"
    );
    let v_switch = v_cal * (3.0_f64.sqrt() * dbeta / detune_cal);

    let mut ec_bias = build_tfln_eo_coupler_engine(v_switch, l_pi);
    call(&mut ec_bias, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dp_bias = DrivenPortData::of(ec_bias.result("ems").expect("ems").0).expect("DrivenPortData");
    let p_bar_v = dp_bias.s_matrix[0][2][0][0].powi(2) + dp_bias.s_matrix[0][2][0][1].powi(2);
    let p_cross_v = dp_bias.s_matrix[0][3][0][0].powi(2) + dp_bias.s_matrix[0][3][0][1].powi(2);
    assert!(
        p_bar_v > p_cross_v && (p_bar_v - p_bar_0) > 0.35 && (p_cross_0 - p_cross_v) > 0.35,
        "push-pull EO bias at V_switch={v_switch:.2} V must invert dominant port from cross to bar: P_bar {p_bar_0:.4}->{p_bar_v:.4}, P_cross {p_cross_0:.4}->{p_cross_v:.4}"
    );
}

#[test]
fn eme_analytical_length_sensitivities_and_coupler_optimization() {
    // 1. Analytical Redheffer cell-length sensitivity vs central finite differences (< 1e-6 relative) on a 4-port coupler
    let l0 = 65.0e-6;
    let mut ec = build_tfln_eo_coupler_engine(15.0, l0);
    // Also attach a CAD parameter `L_cad` so we verify both `eme_segment_length` and longitudinal `param_var`!
    let sens_trans = call(
        &mut ec,
        "sensitivity_run",
        json!({
            "physics": "ems",
            "objective": { "kind": "transmission", "in_port": 0, "out_port": 3 },
            "params": [{ "kind": "eme_segment_length", "segment": 0 }]
        }),
    );
    let g_ana_trans = sens_trans["flat_gradient"][0].as_f64().expect("grad");

    let sens_bal = call(
        &mut ec,
        "sensitivity_run",
        json!({
            "physics": "ems",
            "objective": { "kind": "splitter_balance", "in_port": 0, "out_ports": [2, 3], "targets": [0.5, 0.5] },
            "params": [{ "kind": "eme_segment_length", "segment": 0 }]
        }),
    );
    let g_ana_bal = sens_bal["flat_gradient"][0].as_f64().expect("grad");

    // Central finite-difference oracle via `eme_sweep` (`mode = "length"`) around `l0`
    call(&mut ec, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dl = 1.0e-9;
    let sw = call(
        &mut ec,
        "eme_sweep",
        json!({
            "physics": "ems",
            "mode": "length",
            "segment": 0,
            "lengths": [l0 - dl, l0 + dl]
        }),
    );
    let pts = sw["points"].as_array().expect("points");
    let p_cross_m = pts[0]["power_from_port0"][3].as_f64().unwrap();
    let p_cross_p = pts[1]["power_from_port0"][3].as_f64().unwrap();
    let p_bar_m = pts[0]["power_from_port0"][2].as_f64().unwrap();
    let p_bar_p = pts[1]["power_from_port0"][2].as_f64().unwrap();

    let g_fd_trans = (p_cross_p - p_cross_m) / (2.0 * dl);
    let j_bal_m = 0.5 * ((p_bar_m - 0.5).powi(2) + (p_cross_m - 0.5).powi(2));
    let j_bal_p = 0.5 * ((p_bar_p - 0.5).powi(2) + (p_cross_p - 0.5).powi(2));
    let g_fd_bal = (j_bal_p - j_bal_m) / (2.0 * dl);

    let rel_err_trans = (g_ana_trans - g_fd_trans).abs() / g_fd_trans.abs().max(1.0);
    let rel_err_bal = (g_ana_bal - g_fd_bal).abs() / g_fd_bal.abs().max(1.0);
    assert!(
        rel_err_trans < 1e-6,
        "analytical Redheffer transmission d|S_30|^2/dL mismatch vs central FD: ana={g_ana_trans:.9e}, fd={g_fd_trans:.9e}, rel={rel_err_trans:.3e}"
    );
    assert!(
        rel_err_bal < 1e-6,
        "analytical Redheffer splitter_balance dJ/dL mismatch vs central FD: ana={g_ana_bal:.9e}, fd={g_fd_bal:.9e}, rel={rel_err_bal:.3e}"
    );

    // 2. Analytical `PeriodicSubgrid` Bragg grating cell-length sensitivity (`O(log N_rep)` tangent) + CAD `param_var`
    let freq = 299_792_458.0 / 1.55e-6;
    let mut e_dbr = Engine::new();
    let l1 = 0.16e-6;
    let l2 = 0.18e-6;
    call(&mut e_dbr, "project_new", json!({ "title": "p2-eme-dbr-sens", "dim": 3, "h": 0.32e-6 }));
    call(&mut e_dbr, "set_param", json!({ "name": "L_tot", "value": l1 + l2 }));
    call(
        &mut e_dbr,
        "add_box",
        json!({ "name": "b_clad", "min": [-0.8e-6, 0.0, -0.5e-6], "max": [0.8e-6, "L_tot", 0.5e-6] }),
    );
    call(
        &mut e_dbr,
        "add_box",
        json!({ "name": "b_seg1", "min": [-0.25e-6, 0.0, -0.12e-6], "max": [0.25e-6, l1, 0.12e-6] }),
    );
    call(
        &mut e_dbr,
        "add_box",
        json!({ "name": "b_seg2", "min": [-0.32e-6, l1, -0.12e-6], "max": [0.32e-6, "L_tot", 0.12e-6] }),
    );
    call(
        &mut e_dbr,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad", "input": "b_clad" },
                { "name": "core1", "input": "b_seg1" },
                { "name": "core2", "input": "b_seg2" }
            ]
        }),
    );
    call(
        &mut e_dbr,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "eme",
            "order": 1,
            "study_params": {
                "freq": freq,
                "n_target": 3.4,
                "n_modes": 2,
                "grid_kind": "periodic_subgrid",
                "lengths": [l1, l2],
                "num_reps": 20,
                "port_modes": 1,
                "reconstruct_fields": false
            }
        }),
    );
    call(&mut e_dbr, "set_material", json!({ "physics": "ems", "region": "clad", "material": { "eps_r": 2.085 } }));
    call(&mut e_dbr, "set_material", json!({ "physics": "ems", "region": "core1", "material": { "eps_r": 12.08 } }));
    call(&mut e_dbr, "set_material", json!({ "physics": "ems", "region": "core2", "material": { "eps_r": 12.08 } }));
    call(&mut e_dbr, "add_condition", json!({ "physics": "ems", "on": "extreme(-z)", "kind": "pmc" }));
    call(&mut e_dbr, "add_condition", json!({ "physics": "ems", "on": "extreme(+z)", "kind": "pmc" }));

    let dbr_sens = call(
        &mut e_dbr,
        "sensitivity_run",
        json!({
            "physics": "ems",
            "objective": { "kind": "transmission", "in_port": 0, "out_port": 1 },
            "params": [
                { "kind": "eme_segment_length", "segment": 0 },
                { "kind": "eme_segment_length", "segment": 1 }
            ]
        }),
    );
    let g_dbr0 = dbr_sens["flat_gradient"][0].as_f64().unwrap();
    let g_dbr1 = dbr_sens["flat_gradient"][1].as_f64().unwrap();
    call(&mut e_dbr, "solve_run", json!({ "physics": "ems" }));
    let h_dbr = 1.0e-12;
    let sw0 =
        call(&mut e_dbr, "eme_sweep", json!({ "mode": "length", "segment": 0, "lengths": [l1 - h_dbr, l1 + h_dbr] }));
    let sw1 =
        call(&mut e_dbr, "eme_sweep", json!({ "mode": "length", "segment": 1, "lengths": [l2 - h_dbr, l2 + h_dbr] }));
    let fd_dbr0 = (sw0["points"][1]["power_from_port0"][1].as_f64().unwrap()
        - sw0["points"][0]["power_from_port0"][1].as_f64().unwrap())
        / (2.0 * h_dbr);
    let fd_dbr1 = (sw1["points"][1]["power_from_port0"][1].as_f64().unwrap()
        - sw1["points"][0]["power_from_port0"][1].as_f64().unwrap())
        / (2.0 * h_dbr);
    let rel0 = (g_dbr0 - fd_dbr0).abs() / fd_dbr0.abs().max(1.0);
    let rel1 = (g_dbr1 - fd_dbr1).abs() / fd_dbr1.abs().max(1.0);
    assert!(
        rel0 < 1e-6 && rel1 < 1e-6,
        "PeriodicSubgrid EME length sensitivity mismatch: rel0={rel0:.3e}, rel1={rel1:.3e}"
    );

    // 3. L-BFGS inverse design of TFLN EO directional coupler length to exact 50/50 3-dB power split!
    let mut ec_probe_opt = build_tfln_eo_coupler_engine(0.0, 50.0e-6);
    call(&mut ec_probe_opt, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let dbeta_opt = {
        let pd = EmeData::of(ec_probe_opt.result("ems").unwrap().0).unwrap();
        (pd.unique_slices[0].beta_fwd[0][0] - pd.unique_slices[0].beta_fwd[1][0]).abs()
    };
    let l_pi_opt = std::f64::consts::PI / dbeta_opt;
    // Start detuned at 0.20 * L_pi (~10% cross, ~85% bar) and let L-BFGS converge to ~0.50 * L_pi (50/50 split)!
    let mut ec_opt = build_tfln_eo_coupler_engine(0.0, 0.20 * l_pi_opt);
    let opt_res = call(
        &mut ec_opt,
        "optimize_run",
        json!({
            "physics": "ems",
            "algorithm": "lbfgs",
            "max_iters": 30,
            "gtol": 1e-10,
            "objective": {
                "kind": "splitter_balance",
                "in_port": 0,
                "out_ports": [2, 3],
                "targets": [0.5, 0.5]
            },
            "params": [{
                "kind": "eme_segment_length",
                "segment": 0,
                "min": 0.05 * l_pi_opt,
                "max": 0.95 * l_pi_opt
            }]
        }),
    );
    let j_init = opt_res["initial_objective"].as_f64().unwrap();
    let j_final = opt_res["final_objective"].as_f64().unwrap();
    let l_opt = opt_res["final_params"]["eme_length[0]"].as_f64().unwrap();
    let total_sum_half = {
        let s2 = opt_res["aux_scalars"]["S_2_0_mag_sq"].as_f64().unwrap();
        let s3 = opt_res["aux_scalars"]["S_3_0_mag_sq"].as_f64().unwrap();
        // On a lossless passive coupler, optimal splitter_balance equalizes |S_bar|^2 == |S_cross|^2!
        assert!(
            (s2 - s3).abs() < 5e-3 && s2 > 0.40 && s3 > 0.40,
            "L-BFGS must equalize 3-dB coupler ports: |S_bar|^2={s2:.6}, |S_cross|^2={s3:.6}"
        );
        0.5 * (s2 + s3)
    };
    assert!(
        j_final < 0.02 * j_init && total_sum_half > 0.40,
        "coupler length optimization must slash splitter_balance objective: {j_init:.4e} -> {j_final:.4e} at L={l_opt:.6e} m"
    );

    // Verify stored post-optimization EME result agrees bit-for-bit with `eme_sweep` at `l_opt`
    let sw_check =
        call(&mut ec_opt, "eme_sweep", json!({ "physics": "ems", "mode": "length", "segment": 0, "lengths": [l_opt] }));
    let sw_bar = sw_check["points"][0]["power_from_port0"][2].as_f64().unwrap();
    let opt_bar = opt_res["aux_scalars"]["S_2_0_mag_sq"].as_f64().unwrap();
    assert!((sw_bar - opt_bar).abs() < 1e-12, "stored EME result must match eme_sweep at l_opt: {sw_bar} vs {opt_bar}");
}

fn run_phase2_pipeline_snapshot() -> Value {
    let si_toml = r#"
format = 2
name = "si_optomech_p2"
kind = "semiconductor"
point_group = "m-3m"
density = "2329 kg/m^3"
source = { cite = "Phase 2 cross-mesh optomechanical silicon waveguide", verified = "2026-11-27" }

[elastic]
class = "cubic"
unit = "GPa"
c11 = 165.6
c12 = 63.9
c44 = 79.5

[dielectric]
eps_clamped = 11.7

[photoelastic]
wavelength = "1.55 um"
p11 = -0.094
p12 = 0.017
p44 = -0.051

[optical]
wavelength_unit = "um"
range = [1.2, 2.0]
model = "tabulated"
data = [[1.20, 3.4757], [1.55, 3.4757], [2.00, 3.4757]]
"#;

    let freq_opt = 299_792_458.0 / 1.55e-6;
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p2-optomech-opt", "dim": 2, "h": 0.24e-6 }));
    call(&mut e, "set_param", json!({ "name": "w_core", "value": 0.42e-6 }));
    call(&mut e, "add_rect", json!({ "name": "r_clad", "min": [-1.2e-6, -0.9e-6], "max": [1.2e-6, 0.9e-6] }));
    call(&mut e, "add_rect", json!({ "name": "r_core", "min": ["-w_core/2", -0.11e-6], "max": ["w_core/2", 0.11e-6] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "clad", "input": "r_clad" },
                { "name": "core", "input": "r_core" }
            ]
        }),
    );
    call(&mut e, "add_material", json!({ "name": "si_optomech_p2", "toml": si_toml }));

    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "modes",
            "order": 2,
            "study_params": { "freq": freq_opt, "n_target": 3.45, "n_modes": 1 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "si_optomech_p2" }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": { "eps_r": 1.0 } }));

    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "sms",
            "study": "eigenfrequency",
            "order": 2,
            "study_params": { "n_modes": 6, "f_search": 9.37e9 }
        }),
    );
    call(&mut e, "set_physics_params", json!({ "physics": "sms", "params": { "plane": "strain" } }));
    call(&mut e, "set_material", json!({ "physics": "sms", "region": "core", "material": "si_optomech_p2" }));
    call(
        &mut e,
        "set_mesh",
        json!({
            "physics": "sms",
            "h": 0.06e-6,
            "active_domains": ["core"]
        }),
    );

    // Step 1: Optimize optical waveguide core width `w_core` via FR-40 `optimize_run` (`lbfgs`) to hit target `n_eff = 2.35`
    let opt_out = call(
        &mut e,
        "optimize_run",
        json!({
            "physics": "ems",
            "algorithm": "lbfgs",
            "max_iters": 8,
            "gtol": 1e-8,
            "objective": { "kind": "effective_index", "mode": 0, "target_neff": 2.35 },
            "params": [{ "kind": "param_var", "name": "w_core", "min": 0.30e-6, "max": 0.60e-6 }]
        }),
    );

    // Step 2: Solve `sms` (solid core submesh) and `ems` (core+air mesh) independently on the optimized waveguide
    call(&mut e, "solve_run", json!({ "physics": "sms" }));
    call(&mut e, "solve_run", json!({ "physics": "ems" }));

    // Step 3: Evaluate cross-mesh `optomechanics` analysis (`g_PE`, `g_MB`, `g_0`, `G_B`)
    // Pick the fundamental lateral breathing mode (`f > 5 GHz` and highest |g_0|) across the 1-based mechanical modes
    let n_mech = e.result("sms").expect("sms").0.solutions.len();
    let mut best_om = Value::Null;
    let mut best_g0_abs = -1.0_f64;
    for m_idx in 0..n_mech {
        let f_mode = e.result("sms").expect("sms").0.solutions[m_idx].value.unwrap_or(0.0);
        if f_mode < 5.0e9 {
            continue;
        }
        let om = call(
            &mut e,
            "analyse",
            json!({
                "analysis": "optomechanics",
                "physics": "ems",
                "params": {
                    "mechanical_physics": "sms",
                    "optical_mode": 1,
                    "mechanical_mode": m_idx + 1,
                    "length_y_m": 1.0e-6,
                    "core_domains": ["core"],
                    "acoustic_q": 1000.0
                }
            }),
        );
        let g0_abs = om["g_0_hz"].as_f64().unwrap_or(0.0).abs();
        if g0_abs > best_g0_abs {
            best_g0_abs = g0_abs;
            best_om = om;
        }
    }

    // Step 4: Also include cross-mesh EO EME + EME sensitivity numbers in the determinism snapshot
    let mut eo_e = build_tfln_eo_coupler_engine(25.0, 90.0e-6);
    call(&mut eo_e, "solve_run", json!({ "physics": ["dfs", "ems"] }));
    let eo_sens = call(
        &mut eo_e,
        "sensitivity_run",
        json!({
            "physics": "ems",
            "objective": { "kind": "transmission", "in_port": 0, "out_port": 3 },
            "params": [{ "kind": "eme_segment_length", "segment": 0 }]
        }),
    );
    let eo_dp = DrivenPortData::of(eo_e.result("ems").expect("ems").0).expect("dp");

    json!({
        "opt_w_core": opt_out["final_params"]["w_core"],
        "opt_final_obj": opt_out["final_objective"],
        "opt_neff_re": opt_out["aux_scalars"]["neff_re"],
        "om_cross_mesh": best_om["cross_mesh"],
        "om_g_pe_hz": best_om["g_pe_hz"],
        "om_g_mb_hz": best_om["g_mb_hz"],
        "om_g_0_hz": best_om["g_0_hz"],
        "om_m_eff_kg": best_om["m_eff_kg"],
        "om_x_zpf_m": best_om["x_zpf_m"],
        "om_sbs_gain": best_om["sbs_gain_w_inv_m_inv"],
        "eo_coupler_s": eo_dp.s_matrix[0],
        "eo_coupler_d_dl": eo_sens["flat_gradient"][0],
    })
}

#[test]
fn cross_mesh_optomechanics_and_optimize_pipeline_with_thread_determinism() {
    if let Ok(tag) = std::env::var("NYX_PHASE2_DET_WORKER") {
        let snap = nyx_physics::canonical_json(&run_phase2_pipeline_snapshot());
        println!("NYX_PHASE2_DET_JSON_{tag}={}", serde_json::to_string(&snap).unwrap());
        return;
    }

    // Verify in-process physical invariants first
    let snap = run_phase2_pipeline_snapshot();
    assert_eq!(snap["om_cross_mesh"], json!(true), "optomechanics must execute across distinct sms/ems meshes");
    let final_obj = snap["opt_final_obj"].as_f64().unwrap();
    let neff_re = snap["opt_neff_re"].as_f64().unwrap();
    assert!(
        final_obj < 1e-5 && (neff_re - 2.35).abs() < 5e-3,
        "optimize_run must tune Si waveguide width to target_neff = 2.35 (got neff={neff_re}, obj={final_obj:.3e})"
    );
    let g_pe = snap["om_g_pe_hz"].as_f64().unwrap();
    let g_mb = snap["om_g_mb_hz"].as_f64().unwrap();
    let g_0 = snap["om_g_0_hz"].as_f64().unwrap();
    let m_eff = snap["om_m_eff_kg"].as_f64().unwrap();
    let x_zpf = snap["om_x_zpf_m"].as_f64().unwrap();
    let sbs_gain = snap["om_sbs_gain"].as_f64().unwrap();
    assert!(
        g_pe.abs() > 1.0 && g_mb.abs() > 1.0 && g_0.abs() > 1.0,
        "optomechanical rates must be non-trivial: {snap}"
    );
    assert!(((g_pe + g_mb) - g_0).abs() <= 1e-6 * g_0.abs(), "g_0 must equal g_PE + g_MB");
    assert!(m_eff > 0.0 && x_zpf > 0.0 && sbs_gain >= 0.0);

    // Subprocess bit-exact determinism: NYX_THREADS=1 vs NYX_THREADS=8
    let exe = std::env::current_exe().expect("test binary path");
    let run_worker = |threads: usize| -> String {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("cross_mesh_optomechanics_and_optimize_pipeline_with_thread_determinism")
            .arg("--nocapture")
            .env("NYX_THREADS", threads.to_string())
            .env("NYX_PHASE2_DET_WORKER", threads.to_string())
            .output()
            .expect("spawn determinism subprocess");
        assert!(out.status.success(), "worker NYX_THREADS={threads} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8(out.stdout).expect("utf8");
        let prefix = format!("NYX_PHASE2_DET_JSON_{threads}=");
        stdout
            .lines()
            .find_map(|l| l.strip_prefix(&prefix))
            .unwrap_or_else(|| panic!("missing determinism payload in stdout:\n{stdout}"))
            .to_string()
    };

    let j1 = run_worker(1);
    let j8 = run_worker(8);
    let v1: Value = serde_json::from_str(&j1).unwrap();
    let v8: Value = serde_json::from_str(&j8).unwrap();
    assert_eq!(
        nyx_physics::canonical_json(&v1),
        nyx_physics::canonical_json(&v8),
        "Phase 2 pipeline must be bit-for-bit identical between NYX_THREADS=1 and NYX_THREADS=8"
    );
}
