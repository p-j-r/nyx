// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 3 cross-FR integration verification suite (`nyx-engine`, ADRs 143–145):
//! 1. 3-way broadband S-parameter parity (`< 1e-3` in `|S_11|^2` and `|S_21|^2`) across frequency-domain
//!    `Study::Driven`, high-order `Study::Dgtd`, and FEM-seeded `Study::Fdtd` (`ModeSource` + `ModeMonitor`)
//!    on a 3D dielectric Fabry–Perot waveguide filter through `Engine::call`.
//! 2. Shared `[nonlinear]` format-2 TOML material parity across `Study::NonlinearDriven` (Harmonic Balance FEM),
//!    `Study::Dgtd`, and `Study::Fdtd`:
//!    - $\chi^{(2)}$ second-harmonic generation (`d33`, quadratic $P_{2\omega} \propto P_{\text{in}}^2$ scaling,
//!      quantitative conversion efficiency parity with CW Harmonic Balance and Boyd $\Gamma^2 L^2$, plus
//!      `[nonlinear.qpm]` spatial poling suppression when $\Delta k = 0$),
//!    - $\chi^{(3)}$ Kerr self-phase modulation (`silicon_nitride`) phase-shift parity ($\Delta\phi_{\text{pulse}}(f_0)
//!      = \Delta\phi_{\text{CW}} / \sqrt{3}$) and crystalline semiconductor TPA + free-carrier absorption (`drude_carrier_step`).
//! 3. Time-domain DFT `FieldMonitor` (`Study::Fdtd` and `Study::Dgtd`) feeding `nyx-analysis` Stratton–Chu
//!    `farfield` NFFFT (`D_max = 1.50`) and `overlap_beam` vectorial Gaussian coupling through `Engine::call("analyse", ...)`.
//! 4. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bitwise `0-ULP` determinism across all integrated pipelines.

use nyx_ems::{DgtdData, FdtdMonitorData, C0};
use nyx_engine::Engine;
use nyx_physics::{canonical_json, DrivenPortData, Scalar, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;
use std::process::Command;

const MU0: f64 = 4.0 * PI * 1e-7;
const EPS0: f64 = 1.0 / (MU0 * C0 * C0);

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn num_scalar(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(v))) => *v,
        other => panic!("missing numeric scalar '{name}', got {other:?} in {:?}", r.scalars),
    }
}

fn unwrap_phase_diff(phi_mod: f64, phi_ref: f64) -> f64 {
    let mut d = phi_mod - phi_ref;
    while d > PI {
        d -= 2.0 * PI;
    }
    while d < -PI {
        d += 2.0 * PI;
    }
    d
}

/// Exact normal-incidence Airy `(R, T)` for a slab of refractive index `n_slab` and thickness `d_m`
/// embedded between semi-infinite leads of refractive index `n_lead`.
fn airy_filter_rt(freq_hz: f64, d_m: f64, n_lead: f64, n_slab: f64) -> (f64, f64) {
    let r12 = (n_lead - n_slab) / (n_lead + n_slab);
    let r21 = -r12;
    let t12_t21 = 1.0 - r12 * r12;
    let phi = 2.0 * PI * freq_hz * n_slab * d_m / C0;
    let (s1, c1) = (-phi).sin_cos();
    let (s2, c2) = (-2.0 * phi).sin_cos();
    let den_re = 1.0 + r12 * r21 * c2;
    let den_im = r12 * r21 * s2;
    let den2 = den_re * den_re + den_im * den_im;
    let num_r_re = r12 + r21 * c2;
    let num_r_im = r21 * s2;
    let r_sq = (num_r_re * num_r_re + num_r_im * num_r_im) / den2;
    let t_sq = (t12_t21 * t12_t21 * (c1 * c1 + s1 * s1)) / den2;
    (r_sq, t_sq)
}

#[test]
fn three_way_broadband_sparameter_parity_driven_dgtd_fdtd() {
    // 3D sub-cutoff dielectric Fabry-Perot waveguide filter:
    //   Cross-section: x in [0, 0.16 um], z in [0, 0.16 um], PMC on x = 0, w_x, PEC on z = 0, h_z.
    //   Leads: eps_r = 2.25 (n_lead = 1.50); Resonator slab: eps_r = 6.25 (n_slab = 2.50), thickness d_slab = 0.40 um.
    //   Reference port planes: y0 = 0.0, y1 = 1.10 um (lead_in = 0.35 um, slab = 0.40 um, lead_out = 0.35 um).
    let w_x = 0.16e-6;
    let h_z = 0.16e-6;
    let y_slab0 = 0.35e-6;
    let y_slab1 = 0.75e-6;
    let y_end = 1.10e-6;
    let eps_lead = 2.25_f64;
    let eps_slab = 6.25_f64;

    let lambdas_m = [1.65e-6, 1.55e-6, 1.45e-6, 1.35e-6];
    let freqs_hz: Vec<f64> = lambdas_m.iter().map(|&lam| C0 / lam).collect();
    let f0 = 0.5 * (freqs_hz[0] + freqs_hz[freqs_hz.len() - 1]);
    let fwidth = 0.32 * f0;

    let build_filter_geom = |e: &mut Engine, y_lo: f64, y_hi: f64, sweep_layers: Value| {
        call(e, "project_new", json!({ "title": "p3-3way-filter", "dim": 3, "h": 0.16e-6 }));
        call(e, "add_box", json!({ "name": "b_in", "min": [0.0, y_lo, 0.0], "max": [w_x, y_slab0, h_z] }));
        call(e, "add_box", json!({ "name": "b_slab", "min": [0.0, y_slab0, 0.0], "max": [w_x, y_slab1, h_z] }));
        call(e, "add_box", json!({ "name": "b_out", "min": [0.0, y_slab1, 0.0], "max": [w_x, y_hi, h_z] }));
        call(
            e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "lead_in", "input": "b_in" },
                    { "name": "slab", "input": "b_slab" },
                    { "name": "lead_out", "input": "b_out" }
                ]
            }),
        );
        call(e, "set_mesh", json!({ "h": 0.16e-6, "sweep": sweep_layers }));
    };

    let sweep_fem = json!([
        { "thickness": y_slab0, "n": 7 },
        { "thickness": y_slab1 - y_slab0, "n": 12 },
        { "thickness": y_end - y_slab1, "n": 7 }
    ]);

    // 1. Frequency-domain `Study::Driven` via Engine
    let mut e_drv = Engine::new();
    build_filter_geom(&mut e_drv, 0.0, y_end, sweep_fem.clone());
    call(
        &mut e_drv,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 2,
            "study_params": {
                "freqs_hz": freqs_hz,
                "prom_tol": 0.0,
                "solver": "direct"
            }
        }),
    );
    for reg in ["lead_in", "lead_out"] {
        call(&mut e_drv, "set_material", json!({ "physics": "ems", "region": reg, "material": { "eps_r": eps_lead } }));
    }
    call(&mut e_drv, "set_material", json!({ "physics": "ems", "region": "slab", "material": { "eps_r": eps_slab } }));
    call(&mut e_drv, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e_drv, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e_drv,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
    );
    call(
        &mut e_drv,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
    );
    call(&mut e_drv, "solve_run", json!({}));
    let dp_drv = DrivenPortData::of(e_drv.result("ems").expect("driven result").0).expect("driven DrivenPortData");

    // 2. High-order `Study::Dgtd` via Engine (`order = 2`, LTS)
    let mut e_dgtd = Engine::new();
    build_filter_geom(&mut e_dgtd, 0.0, y_end, sweep_fem);
    call(
        &mut e_dgtd,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "dgtd",
            "order": 2,
            "study_params": {
                "t_final": 95.0e-15,
                "cfl": 0.75,
                "flux": "upwind",
                "time_stepper": "lts",
                "max_lts_tiers": 4,
                "f0": f0,
                "fwidth": fwidth,
                "freqs_hz": freqs_hz,
                "excited_port": 1
            }
        }),
    );
    for reg in ["lead_in", "lead_out"] {
        call(
            &mut e_dgtd,
            "set_material",
            json!({ "physics": "ems", "region": reg, "material": { "eps_r": eps_lead } }),
        );
    }
    call(&mut e_dgtd, "set_material", json!({ "physics": "ems", "region": "slab", "material": { "eps_r": eps_slab } }));
    call(&mut e_dgtd, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e_dgtd, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e_dgtd,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
    );
    call(
        &mut e_dgtd,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
    );
    call(&mut e_dgtd, "solve_run", json!({}));
    let dp_dgtd = DrivenPortData::of(e_dgtd.result("ems").expect("dgtd result").0).expect("dgtd DrivenPortData");

    // 3. Subpixel `Study::Fdtd` via Engine (CPML y-leads + FEM-seeded `ModeSource` & `ModeMonitor`s at y = 0 and y = y_end)
    let y_fdtd_lo = -0.90e-6;
    let y_fdtd_hi = y_end + 0.90e-6;
    let sweep_fdtd = json!([
        { "thickness": y_slab0 - y_fdtd_lo, "n": 4 },
        { "thickness": y_slab1 - y_slab0, "n": 4 },
        { "thickness": y_fdtd_hi - y_slab1, "n": 4 }
    ]);
    let fdtd_cfg = json!({
        "courant": 0.95,
        "shutoff_energy_ratio": 1e-10,
        "subpixel": "polarized",
        "grid": {
            "auto": { "wavelength_min_m": 1.30e-6, "min_steps_per_wave": 52.0 },
            "uniform_dl": [0.08e-6, 0.0, 0.08e-6],
            "snap_planes": [[1, 0.0], [1, y_slab0], [1, y_slab1], [1, y_end]]
        },
        "boundaries": {
            "faces": [
                { "kind": "pmc" },
                { "kind": "pmc" },
                { "kind": "cpml", "layers": 20, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                { "kind": "cpml", "layers": 20, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                { "kind": "pec" },
                { "kind": "pec" }
            ]
        },
        "sources": [{
            "kind": "mode_source",
            "name": "src_p1",
            "axis": 1,
            "position_m": -0.35e-6,
            "direction": 1.0,
            "mode_index": 1,
            "waveform": {
                "kind": "gaussian_pulse",
                "f0_hz": f0,
                "fwidth_hz": fwidth,
                "amplitude": 1.0,
                "phase_rad": 0.0,
                "peak_offset_sigma": 5.0
            }
        }],
        "monitors": [
            {
                "kind": "mode_monitor",
                "name": "p1",
                "axis": 1,
                "position_m": 0.0,
                "forward_sign": 1.0,
                "mode_index": 1,
                "freqs_hz": freqs_hz
            },
            {
                "kind": "mode_monitor",
                "name": "p2",
                "axis": 1,
                "position_m": y_end,
                "forward_sign": 1.0,
                "mode_index": 1,
                "freqs_hz": freqs_hz
            }
        ]
    });

    let mut e_fdtd = Engine::new();
    build_filter_geom(&mut e_fdtd, y_fdtd_lo, y_fdtd_hi, sweep_fdtd);
    call(
        &mut e_fdtd,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "fdtd",
            "order": 2,
            "study_params": {
                "config_json": fdtd_cfg.to_string()
            }
        }),
    );
    for reg in ["lead_in", "lead_out"] {
        call(
            &mut e_fdtd,
            "set_material",
            json!({ "physics": "ems", "region": reg, "material": { "eps_r": eps_lead } }),
        );
    }
    call(&mut e_fdtd, "set_material", json!({ "physics": "ems", "region": "slab", "material": { "eps_r": eps_slab } }));
    call(&mut e_fdtd, "solve_run", json!({}));
    let dp_fdtd = DrivenPortData::of(e_fdtd.result("ems").expect("fdtd result").0).expect("fdtd DrivenPortData");

    for (fi, &f) in freqs_hz.iter().enumerate() {
        let s11_dr = dp_drv.s_matrix[fi][0][0][0].hypot(dp_drv.s_matrix[fi][0][0][1]).powi(2);
        let s21_dr = dp_drv.s_matrix[fi][1][0][0].hypot(dp_drv.s_matrix[fi][1][0][1]).powi(2);
        let s11_dg = dp_dgtd.s_matrix[fi][0][0][0].hypot(dp_dgtd.s_matrix[fi][0][0][1]).powi(2);
        let s21_dg = dp_dgtd.s_matrix[fi][1][0][0].hypot(dp_dgtd.s_matrix[fi][1][0][1]).powi(2);
        let s11_fd = dp_fdtd.s_matrix[fi][0][0][0].hypot(dp_fdtd.s_matrix[fi][0][0][1]).powi(2);
        let s21_fd = dp_fdtd.s_matrix[fi][1][0][0].hypot(dp_fdtd.s_matrix[fi][1][0][1]).powi(2);
        let (r_airy, t_airy) = airy_filter_rt(f, y_slab1 - y_slab0, eps_lead.sqrt(), eps_slab.sqrt());

        assert!(
            (s11_dr - r_airy).abs() < 5e-4 && (s21_dr - t_airy).abs() < 5e-4,
            "Driven vs Airy mismatch at f={:.2} THz: Driven ({s11_dr:.5}, {s21_dr:.5}) vs Airy ({r_airy:.5}, {t_airy:.5})",
            f * 1e-12
        );
        assert!(
            (s11_dg - s11_dr).abs() < 1e-3 && (s21_dg - s21_dr).abs() < 1e-3,
            "DGTD vs Driven mismatch >= 1e-3 at f={:.2} THz: DGTD ({s11_dg:.5}, {s21_dg:.5}) vs Driven ({s11_dr:.5}, {s21_dr:.5})",
            f * 1e-12
        );
        assert!(
            (s11_fd - s11_dr).abs() < 1e-3 && (s21_fd - s21_dr).abs() < 1e-3,
            "FDTD vs Driven mismatch >= 1e-3 at f={:.2} THz: FDTD ({s11_fd:.5}, {s21_fd:.5}) vs Driven ({s11_dr:.5}, {s21_dr:.5})",
            f * 1e-12
        );
        assert!(
            (s11_dg - s11_fd).abs() < 1e-3 && (s21_dg - s21_fd).abs() < 1e-3,
            "DGTD vs FDTD mismatch >= 1e-3 at f={:.2} THz: DGTD ({s11_dg:.5}, {s21_dg:.5}) vs FDTD ({s11_fd:.5}, {s21_fd:.5})",
            f * 1e-12
        );
        assert!(
            s11_dg + s21_dg <= 1.0 + 1e-3 && s11_fd + s21_fd <= 1.0 + 1e-3,
            "passivity violated at f={:.2} THz: DGTD sum={:.6}, FDTD sum={:.6}",
            f * 1e-12,
            s11_dg + s21_dg,
            s11_fd + s21_fd
        );
    }
}

#[test]
fn shared_nonlinear_toml_shg_qpm_and_kerr_tpa_parity_across_hb_dgtd_fdtd() {
    let lambda1 = 1.55e-6;
    let freq1 = C0 / lambda1;
    let omega1 = 2.0 * PI * freq1;
    let eta0 = (MU0 / EPS0).sqrt();

    // -------------------------------------------------------------------------
    // Part A: Shared `[nonlinear.chi2]` + `[nonlinear.qpm]` TOML across
    //         `Study::NonlinearDriven`, `Study::Dgtd`, and `Study::Fdtd`
    // -------------------------------------------------------------------------
    let n0_ln = 2.14_f64;
    let d33_ln = -25.2e-12_f64; // congruent LiNbO3 d33 = -25.2 pm/V
    let w_x = 0.20e-6;
    let h_z = 0.20e-6;
    let area = w_x * h_z;
    let length_shg = 1.55e-6; // 2 * 0.775 um
    let z0_ln = eta0 / n0_ln;

    let ln_shg_toml = r#"
format = 2
name = "ln_p3_shg"
kind = "crystalline"
point_group = "3m"
source = { cite = "Phase 3 nondispersive LiNbO3-scale chi2 benchmark", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [0.5, 2.5]
model = "constant"
n = 2.14
k = 0.0

[nonlinear.chi2]
wavelength = "1.55 um"
unit = "pm/V"
d33 = -25.2
d31 = -4.6
d22 = 2.1
d15 = -4.6
"#;

    let ln_qpm_cancel_toml = r#"
format = 2
name = "ln_p3_qpm_cancel"
kind = "crystalline"
point_group = "3m"
source = { cite = "Phase 3 QPM periodic domain-cancellation benchmark", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [0.5, 2.5]
model = "constant"
n = 2.14
k = 0.0

[nonlinear.chi2]
wavelength = "1.55 um"
unit = "pm/V"
d33 = -25.2
d31 = -4.6
d22 = 2.1
d15 = -4.6

[nonlinear.qpm]
sign = 1.0
Lambda_qpm = "0.775 um"
axis = [0.0, 1.0, 0.0]
duty_cycle = 0.5
phase_offset = "0.0 um"
"#;

    // Choose weak-depletion peak power `P_lo` (`Gamma * L ~ 0.045 -> eta_CW ~ 2.0e-3`) and `P_hi = 4 * P_lo`:
    let p_shg_lo = 0.25_f64;
    let p_shg_hi = 4.0 * p_shg_lo;
    let gamma_lo =
        ((2.0 * omega1 * omega1 * d33_ln * d33_ln * p_shg_lo) / (n0_ln.powi(3) * EPS0 * C0.powi(3) * area)).sqrt();
    let eta_boyd_lo = (gamma_lo * length_shg).tanh().powi(2);

    // 1. Harmonic Balance `Study::NonlinearDriven` (`process = "shg"`) via Engine
    let mut e_hb_shg = Engine::new();
    call(&mut e_hb_shg, "project_new", json!({ "title": "p3-hb-shg", "dim": 3, "h": 0.20e-6 }));
    call(&mut e_hb_shg, "add_box", json!({ "name": "b_core", "min": [0.0, 0.0, 0.0], "max": [w_x, length_shg, h_z] }));
    call(&mut e_hb_shg, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_core" }] }));
    call(&mut e_hb_shg, "add_material", json!({ "name": "ln_p3_shg", "toml": ln_shg_toml }));
    call(&mut e_hb_shg, "set_mesh", json!({ "h": 0.20e-6, "sweep": [{ "thickness": length_shg, "n": 36 }] }));
    call(
        &mut e_hb_shg,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "nonlinear_driven",
            "order": 2,
            "study_params": {
                "process": "shg",
                "freq": freq1,
                "powers_w": [p_shg_lo, p_shg_hi],
                "drive_port": 0,
                "coupling": "depleted",
                "rtol": 1e-10
            }
        }),
    );
    call(&mut e_hb_shg, "set_material", json!({ "physics": "ems", "region": "core", "material": "ln_p3_shg" }));
    call(&mut e_hb_shg, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e_hb_shg, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e_hb_shg,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(-y)",
            "kind": "lumped_port",
            "params": { "name": "in", "Z_0": z0_ln, "direction": [0.0, 0.0, 1.0] }
        }),
    );
    call(
        &mut e_hb_shg,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(+y)",
            "kind": "lumped_port",
            "params": { "name": "out", "Z_0": z0_ln, "direction": [0.0, 0.0, 1.0] }
        }),
    );
    call(&mut e_hb_shg, "solve_run", json!({}));
    let res_hb_shg = e_hb_shg.result("ems").expect("HB SHG result").0;
    let eta_hb_lo = num_scalar(res_hb_shg, "shg_efficiency step 1");
    let eta_hb_hi = num_scalar(res_hb_shg, "shg_efficiency step 2");
    assert!(
        (eta_hb_lo - eta_boyd_lo).abs() / eta_boyd_lo < 0.02,
        "HB SHG efficiency {eta_hb_lo:.6e} must match Boyd tanh^2 {eta_boyd_lo:.6e}"
    );

    // 2. Time-domain `Study::Dgtd` driven by the EXACT same `ln_p3_shg` & `ln_p3_qpm_cancel` TOMLs:
    //    For TEM port with int |e_t|^2 dS = 1 and wave admittance Y_0 = 1 / z0_ln,
    //    peak cycle-averaged power P_peak = 0.5 * Y_0 * A^2 => A = sqrt(2 * z0_ln * P_peak).
    let run_dgtd_shg = |mat_name: &str, p_peak_w: f64| -> (f64, f64) {
        let amp = (2.0 * z0_ln * p_peak_w).sqrt();
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-dgtd-shg", "dim": 3, "h": 0.20e-6 }));
        call(&mut e, "add_box", json!({ "name": "b_core", "min": [0.0, 0.0, 0.0], "max": [w_x, length_shg, h_z] }));
        call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_core" }] }));
        call(&mut e, "add_material", json!({ "name": "ln_p3_shg", "toml": ln_shg_toml }));
        call(&mut e, "add_material", json!({ "name": "ln_p3_qpm_cancel", "toml": ln_qpm_cancel_toml }));
        call(&mut e, "set_mesh", json!({ "h": 0.20e-6, "sweep": [{ "thickness": length_shg, "n": 24 }] }));
        call(
            &mut e,
            "add_physics",
            json!({
                "physics": "ems",
                "study": "dgtd",
                "order": 2,
                "study_params": {
                    "t_final": 72.0e-15,
                    "cfl": 0.75,
                    "flux": "upwind",
                    "time_stepper": "lts",
                    "max_lts_tiers": 4,
                    "f0": freq1,
                    "fwidth": 0.13 * freq1,
                    "pulse_amplitude": amp,
                    "freqs_hz": [freq1, 2.0 * freq1],
                    "excited_port": 1,
                    "flux_monitors_on": "extreme(+y)",
                    "nonlinear": true
                }
            }),
        );
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": mat_name }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
        );
        call(&mut e, "solve_run", json!({}));
        let td = DgtdData::of(e.result("ems").expect("dgtd shg").0).expect("DgtdData");
        let fm = &td.flux_monitors[0];
        (fm.raw_flux_w_s2[0], fm.raw_flux_w_s2[1])
    };

    let (f1_dg_lo, f2_dg_lo) = run_dgtd_shg("ln_p3_shg", p_shg_lo);
    let (_f1_dg_hi, f2_dg_hi) = run_dgtd_shg("ln_p3_shg", p_shg_hi);
    let (_f1_dg_qpm, f2_dg_qpm) = run_dgtd_shg("ln_p3_qpm_cancel", p_shg_lo);

    let quad_ratio_dg = f2_dg_hi / f2_dg_lo.max(1e-36);
    // Gaussian pulse envelope g(t)^2 has spectral energy |s_2(2 f_0)|^2 = 0.5 * |s_1(f_0)|^2,
    // so CW peak conversion efficiency eta_CW(P_peak) = 2 * raw_flux(2 f_0) / raw_flux(f_0):
    let eta_dg_lo_cw_equiv = 2.0 * (f2_dg_lo / f1_dg_lo.max(1e-36));
    assert!(
        (quad_ratio_dg - 16.0).abs() / 16.0 < 0.06,
        "DGTD SHG power must scale quadratically (4x pump -> 16x SHG), got ratio {quad_ratio_dg:.3}"
    );
    assert!(
        (eta_dg_lo_cw_equiv - eta_hb_lo).abs() / eta_hb_lo < 0.12,
        "DGTD SHG CW-equivalent efficiency {eta_dg_lo_cw_equiv:.5e} must agree with HB {eta_hb_lo:.5e}"
    );
    assert!(
        f2_dg_lo > 25.0 * f2_dg_qpm,
        "DGTD QPM alternating domains at Delta_k = 0 must cancel SHG by >25x: unpoled={f2_dg_lo:.4e} vs poled={f2_dg_qpm:.4e}"
    );

    // 3. Time-domain `Study::Fdtd` driven by the EXACT same `ln_p3_shg` & `ln_p3_qpm_cancel` TOMLs:
    //    `ModeSource` normalizes modal cross-power to 1 W at unit amplitude, so `amplitude = sqrt(P_peak)`.
    let run_fdtd_shg = |mat_name: &str, p_peak_w: f64| -> (f64, f64) {
        let amp = p_peak_w.sqrt();
        let y_lo = -0.65e-6;
        let y_hi = length_shg + 0.65e-6;
        let cfg = json!({
            "courant": 0.95,
            "shutoff_energy_ratio": 1e-10,
            "subpixel": "polarized",
            "nonlinear": true,
            "grid": {
                "auto": { "wavelength_min_m": 0.775e-6, "min_steps_per_wave": 28.0 },
                "uniform_dl": [0.10e-6, 0.0, 0.10e-6],
                "snap_planes": [[1, 0.0], [1, 0.5 * length_shg], [1, length_shg]]
            },
            "boundaries": {
                "faces": [
                    { "kind": "pmc" },
                    { "kind": "pmc" },
                    { "kind": "cpml", "layers": 14, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                    { "kind": "cpml", "layers": 14, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                    { "kind": "pec" },
                    { "kind": "pec" }
                ]
            },
            "sources": [{
                "kind": "mode_source",
                "name": "src_pump",
                "axis": 1,
                "position_m": -0.20e-6,
                "direction": 1.0,
                "mode_index": 1,
                "waveform": {
                    "kind": "gaussian_pulse",
                    "f0_hz": freq1,
                    "fwidth_hz": 0.13 * freq1,
                    "amplitude": amp,
                    "phase_rad": 0.0,
                    "peak_offset_sigma": 4.2
                }
            }],
            "monitors": [{
                "kind": "flux_monitor",
                "name": "out_flux",
                "axis": 1,
                "position_m": length_shg + 0.10e-6,
                "normal_sign": 1.0,
                "freqs_hz": [freq1, 2.0 * freq1]
            }]
        });

        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-fdtd-shg", "dim": 3, "h": 0.20e-6 }));
        call(&mut e, "add_box", json!({ "name": "b_in", "min": [0.0, y_lo, 0.0], "max": [w_x, 0.0, h_z] }));
        call(&mut e, "add_box", json!({ "name": "b_nl", "min": [0.0, 0.0, 0.0], "max": [w_x, length_shg, h_z] }));
        call(&mut e, "add_box", json!({ "name": "b_out", "min": [0.0, length_shg, 0.0], "max": [w_x, y_hi, h_z] }));
        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "lead_in", "input": "b_in" },
                    { "name": "nl_core", "input": "b_nl" },
                    { "name": "lead_out", "input": "b_out" }
                ]
            }),
        );
        call(&mut e, "add_material", json!({ "name": "ln_p3_shg", "toml": ln_shg_toml }));
        call(&mut e, "add_material", json!({ "name": "ln_p3_qpm_cancel", "toml": ln_qpm_cancel_toml }));
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": 0.20e-6,
                "sweep": [
                    { "thickness": -y_lo, "n": 3 },
                    { "thickness": length_shg, "n": 4 },
                    { "thickness": y_hi - length_shg, "n": 3 }
                ]
            }),
        );
        call(
            &mut e,
            "add_physics",
            json!({ "physics": "ems", "study": "fdtd", "order": 2, "study_params": { "config_json": cfg.to_string() } }),
        );
        for reg in ["lead_in", "lead_out"] {
            call(
                &mut e,
                "set_material",
                json!({ "physics": "ems", "region": reg, "material": { "eps_r": n0_ln * n0_ln } }),
            );
        }
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "nl_core", "material": mat_name }));
        call(&mut e, "solve_run", json!({}));
        let md = FdtdMonitorData::of(e.result("ems").expect("fdtd shg").0).expect("FdtdMonitorData");
        let fm = &md.flux_monitors[0];
        (fm.raw_flux_w_s2[0], fm.raw_flux_w_s2[1])
    };

    let (f1_fd_lo, f2_fd_lo) = run_fdtd_shg("ln_p3_shg", p_shg_lo);
    let (_f1_fd_hi, f2_fd_hi) = run_fdtd_shg("ln_p3_shg", p_shg_hi);
    let (_f1_fd_qpm, f2_fd_qpm) = run_fdtd_shg("ln_p3_qpm_cancel", p_shg_lo);

    let quad_ratio_fd = f2_fd_hi / f2_fd_lo.max(1e-36);
    let eta_fd_lo_cw_equiv = 2.0 * (f2_fd_lo / f1_fd_lo.max(1e-36));
    assert!(
        (quad_ratio_fd - 16.0).abs() / 16.0 < 0.06 && eta_hb_hi > 3.8 * eta_hb_lo,
        "FDTD SHG power must scale quadratically (4x pump -> 16x SHG), got {quad_ratio_fd:.3}"
    );
    assert!(
        (eta_fd_lo_cw_equiv - eta_hb_lo).abs() / eta_hb_lo < 0.12
            && (eta_fd_lo_cw_equiv - eta_dg_lo_cw_equiv).abs() / eta_hb_lo < 0.12,
        "FDTD SHG efficiency {eta_fd_lo_cw_equiv:.5e} must agree with DGTD {eta_dg_lo_cw_equiv:.5e} and HB {eta_hb_lo:.5e}"
    );
    assert!(
        f2_fd_lo > 20.0 * f2_fd_qpm,
        "FDTD QPM alternating domains at Delta_k = 0 must cancel SHG by >20x: unpoled={f2_fd_lo:.4e} vs poled={f2_fd_qpm:.4e}"
    );

    // -------------------------------------------------------------------------
    // Part B: Shared `[nonlinear.chi3]` Kerr SPM (`silicon_nitride`) + TPA/FCA
    //         across `Study::NonlinearDriven`, `Study::Dgtd`, and `Study::Fdtd`
    // -------------------------------------------------------------------------
    let p_kerr_lo = 1.0_f64;
    let p_kerr_mid = 200.0_f64;
    let p_kerr_hi = 400.0_f64;
    let length_kerr = 1.55e-6;
    // Luke 2015 Si3N4 refractive index at 1.55 um is n0 = 1.9963, n_2 = 2.4e-19 m^2/W:
    let n0_sin = 1.9963_f64;
    let z0_sin = eta0 / n0_sin;

    // 1. Harmonic Balance `Study::NonlinearDriven` (`process = "chi3"`, shipped `silicon_nitride`)
    let mut e_hb_kerr = Engine::new();
    call(&mut e_hb_kerr, "project_new", json!({ "title": "p3-hb-kerr", "dim": 3, "h": 0.20e-6 }));
    call(
        &mut e_hb_kerr,
        "add_box",
        json!({ "name": "b_core", "min": [0.0, 0.0, 0.0], "max": [w_x, length_kerr, h_z] }),
    );
    call(&mut e_hb_kerr, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_core" }] }));
    call(&mut e_hb_kerr, "set_mesh", json!({ "h": 0.20e-6, "sweep": [{ "thickness": length_kerr, "n": 24 }] }));
    call(
        &mut e_hb_kerr,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "nonlinear_driven",
            "order": 2,
            "study_params": {
                "process": "chi3",
                "freq": freq1,
                "powers_w": [p_kerr_lo, p_kerr_mid, p_kerr_hi],
                "drive_port": 0,
                "rtol": 1e-11
            }
        }),
    );
    call(&mut e_hb_kerr, "set_material", json!({ "physics": "ems", "region": "core", "material": "silicon_nitride" }));
    call(&mut e_hb_kerr, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e_hb_kerr, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e_hb_kerr,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(-y)",
            "kind": "lumped_port",
            "params": { "name": "in", "Z_0": z0_sin, "direction": [0.0, 0.0, 1.0] }
        }),
    );
    call(
        &mut e_hb_kerr,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "extreme(+y)",
            "kind": "lumped_port",
            "params": { "name": "out", "Z_0": z0_sin, "direction": [0.0, 0.0, 1.0] }
        }),
    );
    call(&mut e_hb_kerr, "solve_run", json!({}));
    let res_hb_kerr = e_hb_kerr.result("ems").expect("HB Kerr").0;
    let phi_hb = |step: usize| {
        let bre = num_scalar(res_hb_kerr, &format!("harmonic 1 port 'out' b_out_re step {step} (W^0.5)"));
        let bim = num_scalar(res_hb_kerr, &format!("harmonic 1 port 'out' b_out_im step {step} (W^0.5)"));
        bim.atan2(bre)
    };
    let dphi_hb_hi = unwrap_phase_diff(phi_hb(3), phi_hb(1));
    let dphi_hb_mid = unwrap_phase_diff(phi_hb(2), phi_hb(1));
    let dphi_exact_hi = -(2.0 * PI / lambda1) * 2.4e-19 * ((p_kerr_hi - p_kerr_lo) / area) * length_kerr;
    assert!(
        dphi_hb_hi < 0.0 && (dphi_hb_hi - dphi_exact_hi).abs() / dphi_exact_hi.abs() < 0.02,
        "HB Kerr phase shift {dphi_hb_hi:.6e} must match -k0*n2*DeltaI*L = {dphi_exact_hi:.6e}"
    );

    // 2. Time-domain `Study::Dgtd` Kerr SPM on shipped `silicon_nitride`
    let run_dgtd_s21 = |mat_name: &str, custom_toml: Option<(&str, &str)>, z0_port: f64, p_peak_w: f64| -> [f64; 2] {
        let amp = (2.0 * z0_port * p_peak_w).sqrt();
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-dgtd-chi3", "dim": 3, "h": 0.20e-6 }));
        call(&mut e, "add_box", json!({ "name": "b_core", "min": [0.0, 0.0, 0.0], "max": [w_x, length_kerr, h_z] }));
        call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_core" }] }));
        if let Some((nm, tml)) = custom_toml {
            call(&mut e, "add_material", json!({ "name": nm, "toml": tml }));
        }
        call(&mut e, "set_mesh", json!({ "h": 0.20e-6, "sweep": [{ "thickness": length_kerr, "n": 16 }] }));
        call(
            &mut e,
            "add_physics",
            json!({
                "physics": "ems",
                "study": "dgtd",
                "order": 2,
                "study_params": {
                    "t_final": 78.0e-15,
                    "cfl": 0.75,
                    "flux": "upwind",
                    "time_stepper": "lts",
                    "max_lts_tiers": 4,
                    "f0": freq1,
                    "fwidth": 0.13 * freq1,
                    "pulse_amplitude": amp,
                    "freqs_hz": [freq1],
                    "excited_port": 1,
                    "nonlinear": true
                }
            }),
        );
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": mat_name }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
        );
        call(&mut e, "solve_run", json!({}));
        DrivenPortData::of(e.result("ems").expect("dgtd chi3").0).expect("DrivenPortData").s_matrix[0][1][0]
    };

    let s21_dg_lo = run_dgtd_s21("silicon_nitride", None, z0_sin, p_kerr_lo);
    let s21_dg_mid = run_dgtd_s21("silicon_nitride", None, z0_sin, p_kerr_mid);
    let s21_dg_hi = run_dgtd_s21("silicon_nitride", None, z0_sin, p_kerr_hi);
    let dphi_dg_mid = unwrap_phase_diff(s21_dg_mid[1].atan2(s21_dg_mid[0]), s21_dg_lo[1].atan2(s21_dg_lo[0]));
    let dphi_dg_hi = unwrap_phase_diff(s21_dg_hi[1].atan2(s21_dg_hi[0]), s21_dg_lo[1].atan2(s21_dg_lo[0]));
    // Gaussian pulse carrier-center spectral Kerr phase shift satisfies dphi_CW = sqrt(3) * dphi_pulse(f0):
    let dphi_dg_hi_cw_equiv = 3.0_f64.sqrt() * dphi_dg_hi;
    assert!(
        dphi_dg_hi < 0.0 && ((dphi_dg_hi / dphi_dg_mid) - (dphi_hb_hi / dphi_hb_mid)).abs() < 0.08,
        "DGTD Kerr phase shift must be negative and linear in power: mid={dphi_dg_mid:.5e}, hi={dphi_dg_hi:.5e}"
    );
    assert!(
        (dphi_dg_hi_cw_equiv - dphi_hb_hi).abs() / dphi_hb_hi.abs() < 0.10,
        "DGTD Kerr CW-equivalent phase shift {dphi_dg_hi_cw_equiv:.5e} must agree with HB {dphi_hb_hi:.5e}"
    );

    // 3. Time-domain `Study::Fdtd` Kerr SPM on shipped `silicon_nitride`
    let run_fdtd_s21 = |mat_name: &str,
                        custom_toml: Option<(&str, &str)>,
                        eps_lead_val: f64,
                        p_peak_w: f64|
     -> [f64; 2] {
        let amp = p_peak_w.sqrt();
        let y_lo = -0.65e-6;
        let y_hi = length_kerr + 0.65e-6;
        let cfg = json!({
            "courant": 0.95,
            "shutoff_energy_ratio": 1e-10,
            "subpixel": "polarized",
            "nonlinear": true,
            "grid": {
                "auto": { "wavelength_min_m": 1.50e-6, "min_steps_per_wave": 30.0 },
                "uniform_dl": [0.10e-6, 0.0, 0.10e-6],
                "snap_planes": [[1, 0.0], [1, length_kerr]]
            },
            "boundaries": {
                "faces": [
                    { "kind": "pmc" },
                    { "kind": "pmc" },
                    { "kind": "cpml", "layers": 14, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                    { "kind": "cpml", "layers": 14, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-8 },
                    { "kind": "pec" },
                    { "kind": "pec" }
                ]
            },
            "sources": [{
                "kind": "mode_source",
                "name": "src_in",
                "axis": 1,
                "position_m": -0.20e-6,
                "direction": 1.0,
                "mode_index": 1,
                "waveform": {
                    "kind": "gaussian_pulse",
                    "f0_hz": freq1,
                    "fwidth_hz": 0.13 * freq1,
                    "amplitude": amp,
                    "phase_rad": 0.0,
                    "peak_offset_sigma": 4.2
                }
            }],
            "monitors": [
                {
                    "kind": "mode_monitor",
                    "name": "in",
                    "axis": 1,
                    "position_m": 0.0,
                    "forward_sign": 1.0,
                    "mode_index": 1,
                    "freqs_hz": [freq1]
                },
                {
                    "kind": "mode_monitor",
                    "name": "out",
                    "axis": 1,
                    "position_m": length_kerr,
                    "forward_sign": 1.0,
                    "mode_index": 1,
                    "freqs_hz": [freq1]
                }
            ]
        });

        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-fdtd-chi3", "dim": 3, "h": 0.20e-6 }));
        call(&mut e, "add_box", json!({ "name": "b_in", "min": [0.0, y_lo, 0.0], "max": [w_x, 0.0, h_z] }));
        call(&mut e, "add_box", json!({ "name": "b_nl", "min": [0.0, 0.0, 0.0], "max": [w_x, length_kerr, h_z] }));
        call(&mut e, "add_box", json!({ "name": "b_out", "min": [0.0, length_kerr, 0.0], "max": [w_x, y_hi, h_z] }));
        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "lead_in", "input": "b_in" },
                    { "name": "nl_core", "input": "b_nl" },
                    { "name": "lead_out", "input": "b_out" }
                ]
            }),
        );
        if let Some((nm, tml)) = custom_toml {
            call(&mut e, "add_material", json!({ "name": nm, "toml": tml }));
        }
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": 0.20e-6,
                "sweep": [
                    { "thickness": -y_lo, "n": 3 },
                    { "thickness": length_kerr, "n": 4 },
                    { "thickness": y_hi - length_kerr, "n": 3 }
                ]
            }),
        );
        call(
            &mut e,
            "add_physics",
            json!({ "physics": "ems", "study": "fdtd", "order": 2, "study_params": { "config_json": cfg.to_string() } }),
        );
        for reg in ["lead_in", "lead_out"] {
            call(
                &mut e,
                "set_material",
                json!({ "physics": "ems", "region": reg, "material": { "eps_r": eps_lead_val } }),
            );
        }
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "nl_core", "material": mat_name }));
        call(&mut e, "solve_run", json!({}));
        DrivenPortData::of(e.result("ems").expect("fdtd chi3").0).expect("DrivenPortData").s_matrix[0][1][0]
    };

    let s21_fd_lo = run_fdtd_s21("silicon_nitride", None, n0_sin * n0_sin, p_kerr_lo);
    let s21_fd_mid = run_fdtd_s21("silicon_nitride", None, n0_sin * n0_sin, p_kerr_mid);
    let s21_fd_hi = run_fdtd_s21("silicon_nitride", None, n0_sin * n0_sin, p_kerr_hi);
    let dphi_fd_mid = unwrap_phase_diff(s21_fd_mid[1].atan2(s21_fd_mid[0]), s21_fd_lo[1].atan2(s21_fd_lo[0]));
    let dphi_fd_hi = unwrap_phase_diff(s21_fd_hi[1].atan2(s21_fd_hi[0]), s21_fd_lo[1].atan2(s21_fd_lo[0]));
    let dphi_fd_hi_cw_equiv = 3.0_f64.sqrt() * dphi_fd_hi;
    assert!(
        dphi_fd_hi < 0.0 && ((dphi_fd_hi / dphi_fd_mid) - (dphi_hb_hi / dphi_hb_mid)).abs() < 0.08,
        "FDTD Kerr phase shift must be negative and linear in power: mid={dphi_fd_mid:.5e}, hi={dphi_fd_hi:.5e}"
    );
    assert!(
        (dphi_fd_hi_cw_equiv - dphi_hb_hi).abs() / dphi_hb_hi.abs() < 0.10
            && (dphi_fd_hi - dphi_dg_hi).abs() / dphi_dg_hi.abs() < 0.10,
        "FDTD Kerr CW-equivalent phase shift {dphi_fd_hi_cw_equiv:.5e} must agree with DGTD {dphi_dg_hi_cw_equiv:.5e} and HB {dphi_hb_hi:.5e}"
    );

    // 4. Semiconductor TPA + Free-Carrier Absorption (`drude_carrier_step`) in BOTH `Study::Dgtd` and `Study::Fdtd`
    let si_tpa_only_toml = r#"
format = 2
name = "si_p3_tpa_only"
kind = "semiconductor"
point_group = "m-3m"
source = { cite = "Phase 3 TPA-only silicon benchmark", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [1.2, 2.0]
model = "constant"
n = 3.4757
k = 0.0

[nonlinear.chi3]
wavelength = "1.55 um"
n_2 = "4.5e-18 m^2/W"
beta_tpa = "5.0e-11 m/W"
sigma_fca = "0.0 m^2"
tau_c = "40 fs"
"#;

    let si_tpa_fca_toml = r#"
format = 2
name = "si_p3_tpa_fca"
kind = "semiconductor"
point_group = "m-3m"
source = { cite = "Phase 3 TPA + FCA/FCD carrier-dynamics silicon benchmark", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [1.2, 2.0]
model = "constant"
n = 3.4757
k = 0.0

[nonlinear.chi3]
wavelength = "1.55 um"
n_2 = "4.5e-18 m^2/W"
beta_tpa = "5.0e-11 m/W"
sigma_fca = "5.0e-20 m^2"
k_fcd = "-1.35e-27 m^3"
tau_c = "40 fs"
"#;

    let n0_si = 3.4757_f64;
    let z0_si = eta0 / n0_si;
    let p_si_lo = 1.0e-3_f64;
    let p_si_hi = 36.0_f64;

    let s21_dg_si_lin = run_dgtd_s21("si_p3_tpa_only", Some(("si_p3_tpa_only", si_tpa_only_toml)), z0_si, p_si_lo);
    let s21_dg_si_tpa = run_dgtd_s21("si_p3_tpa_only", Some(("si_p3_tpa_only", si_tpa_only_toml)), z0_si, p_si_hi);
    let s21_dg_si_fca = run_dgtd_s21("si_p3_tpa_fca", Some(("si_p3_tpa_fca", si_tpa_fca_toml)), z0_si, p_si_hi);
    let t_dg_lin = s21_dg_si_lin[0].hypot(s21_dg_si_lin[1]).powi(2);
    let t_dg_tpa = s21_dg_si_tpa[0].hypot(s21_dg_si_tpa[1]).powi(2);
    let t_dg_fca = s21_dg_si_fca[0].hypot(s21_dg_si_fca[1]).powi(2);
    assert!(
        t_dg_tpa < 0.975 * t_dg_lin && t_dg_fca < 0.975 * t_dg_tpa,
        "DGTD TPA + drude_carrier_step FCA must attenuate transmission progressively: T_lin={t_dg_lin:.4}, T_tpa={t_dg_tpa:.4}, T_tpa_fca={t_dg_fca:.4}"
    );

    let s21_fd_si_lin =
        run_fdtd_s21("si_p3_tpa_only", Some(("si_p3_tpa_only", si_tpa_only_toml)), n0_si * n0_si, p_si_lo);
    let s21_fd_si_tpa =
        run_fdtd_s21("si_p3_tpa_only", Some(("si_p3_tpa_only", si_tpa_only_toml)), n0_si * n0_si, p_si_hi);
    let s21_fd_si_fca = run_fdtd_s21("si_p3_tpa_fca", Some(("si_p3_tpa_fca", si_tpa_fca_toml)), n0_si * n0_si, p_si_hi);
    let t_fd_lin = s21_fd_si_lin[0].hypot(s21_fd_si_lin[1]).powi(2);
    let t_fd_tpa = s21_fd_si_tpa[0].hypot(s21_fd_si_tpa[1]).powi(2);
    let t_fd_fca = s21_fd_si_fca[0].hypot(s21_fd_si_fca[1]).powi(2);
    assert!(
        t_fd_tpa < 0.975 * t_fd_lin && t_fd_fca < 0.975 * t_fd_tpa,
        "FDTD TPA + drude_carrier_step FCA must attenuate transmission progressively: T_lin={t_fd_lin:.4}, T_tpa={t_fd_tpa:.4}, T_tpa_fca={t_fd_fca:.4}"
    );
}

#[test]
fn time_domain_dft_to_nffft_farfield_and_beam_overlap_through_engine() {
    let f_dip = C0; // lambda = 1.0 m in vacuum

    // 1. 3D `Study::Fdtd` Hertzian z-dipole -> DFT `FieldMonitor` -> `Engine::call("analyse", "farfield" & "overlap_beam")`
    let fdtd_dip_cfg = json!({
        "courant": 0.95,
        "shutoff_energy_ratio": 1e-8,
        "grid": {
            "auto": { "wavelength_min_m": 1.0, "min_steps_per_wave": 20.0 }
        },
        "sources": [{
            "kind": "point_dipole",
            "name": "ez_dipole",
            "position": [0.0, 0.0, 0.0],
            "component": 2,
            "magnetic": false,
            "waveform": {
                "kind": "gaussian_pulse",
                "f0_hz": f_dip,
                "fwidth_hz": 0.25 * f_dip,
                "amplitude": 1.0,
                "phase_rad": 0.0,
                "peak_offset_sigma": 4.5
            }
        }],
        "monitors": [{
            "kind": "field_monitor",
            "name": "vol_dft",
            "freqs_hz": [f_dip]
        }]
    });

    let mut e_fdtd = Engine::new();
    call(&mut e_fdtd, "project_new", json!({ "title": "p3-fdtd-nffft", "dim": 3, "h": 0.16 }));
    call(&mut e_fdtd, "add_box", json!({ "name": "b_vac", "min": [-0.90, -0.90, -0.90], "max": [0.90, 0.90, 0.90] }));
    call(&mut e_fdtd, "add_overlay", json!({ "regions": [{ "name": "vac", "input": "b_vac" }] }));
    call(&mut e_fdtd, "set_mesh", json!({ "h": 0.16, "sweep": [{ "thickness": 1.80, "n": 12 }] }));
    call(
        &mut e_fdtd,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "fdtd",
            "order": 2,
            "study_params": { "config_json": fdtd_dip_cfg.to_string() }
        }),
    );
    call(&mut e_fdtd, "set_material", json!({ "physics": "ems", "region": "vac", "material": { "eps_r": 1.0 } }));
    call(&mut e_fdtd, "solve_run", json!({}));

    let ff_fdtd = call(
        &mut e_fdtd,
        "analyse",
        json!({
            "analysis": "farfield",
            "physics": "ems",
            "params": {
                "surface": { "kind": "box", "lo": [-0.45, -0.45, -0.45], "hi": [0.45, 0.45, 0.45], "sampling": "grid" },
                "n_theta": 37,
                "n_phi": 24
            }
        }),
    );
    let d_max_fdtd = ff_fdtd["D_max"].as_f64().expect("FDTD NFFFT D_max");
    let th_fdtd = ff_fdtd["D_max_theta_deg"].as_f64().expect("FDTD NFFFT D_max_theta_deg");
    let p_rad_ff = ff_fdtd["P_rad_farfield_w"].as_f64().expect("FDTD P_rad_farfield_w");
    let p_rad_flux = ff_fdtd["P_rad_flux_w"].as_f64().expect("FDTD P_rad_flux_w");
    assert!(
        (d_max_fdtd - 1.50).abs() < 0.06,
        "FDTD -> Engine farfield expected Hertzian dipole D_max ~ 1.50, got {d_max_fdtd:.4}"
    );
    assert!(
        (th_fdtd - 90.0).abs() < 1e-6 && p_rad_ff > 0.0 && p_rad_flux > 0.0,
        "FDTD -> Engine farfield expected equatorial peak (90 deg) and positive radiated power, got theta={th_fdtd}, P_ff={p_rad_ff:e}, P_flux={p_rad_flux:e}"
    );
    assert!(
        (p_rad_ff - p_rad_flux).abs() / p_rad_flux < 0.10,
        "FDTD far-field power {p_rad_ff:e} and surface Poynting flux {p_rad_flux:e} must agree within 10%"
    );

    let ob_fdtd = call(
        &mut e_fdtd,
        "analyse",
        json!({
            "analysis": "overlap_beam",
            "physics": "ems",
            "params": {
                "plane": {
                    "axis": "x",
                    "position": 0.45,
                    "lo": [-0.45, -0.45],
                    "hi": [0.45, 0.45],
                    "n": [16, 16],
                    "normal_sign": 1.0
                },
                "beam": {
                    "kind": "gaussian",
                    "waist_m": 0.40,
                    "focus": [0.45, 0.0, 0.0],
                    "direction": [1.0, 0.0, 0.0],
                    "polarization": [0.0, 0.0, 1.0]
                }
            }
        }),
    );
    let eta_beam = ob_fdtd["eta"].as_f64().expect("FDTD overlap_beam eta");
    let p_field_w = ob_fdtd["P_field_w"].as_f64().expect("FDTD overlap_beam P_field_w");
    assert!(
        eta_beam > 0.20 && eta_beam <= 1.0 && p_field_w > 0.0,
        "FDTD -> Engine overlap_beam expected valid coupling efficiency in (0.20, 1.0], got eta={eta_beam:.4}, P_field_w={p_field_w:e}"
    );

    // 2. 3D `Study::Dgtd` Hertzian z-dipole -> DFT `CellField`s -> `Engine::call("analyse", "farfield")`
    //    (odd sweep count n = 9 places [0, 0, 0] strictly inside the central y-layer away from element interfaces)
    let mut e_dgtd = Engine::new();
    call(&mut e_dgtd, "project_new", json!({ "title": "p3-dgtd-nffft", "dim": 3, "h": 0.34 }));
    call(&mut e_dgtd, "add_box", json!({ "name": "b_vac", "min": [-0.85, -0.85, -0.85], "max": [0.85, 0.85, 0.85] }));
    call(&mut e_dgtd, "add_overlay", json!({ "regions": [{ "name": "vac", "input": "b_vac" }] }));
    call(&mut e_dgtd, "set_mesh", json!({ "h": 0.34, "sweep": [{ "thickness": 1.70, "n": 9 }] }));
    call(
        &mut e_dgtd,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "dgtd",
            "order": 2,
            "study_params": {
                "t_final": 15.0e-9,
                "cfl": 0.75,
                "flux": "upwind",
                "time_stepper": "lts",
                "max_lts_tiers": 4,
                "f0": f_dip,
                "fwidth": 0.25 * f_dip,
                "freqs_hz": [f_dip],
                "absorbing_on": "extreme(-x), extreme(+x), extreme(-y), extreme(+y), extreme(-z), extreme(+z)"
            }
        }),
    );
    call(&mut e_dgtd, "set_material", json!({ "physics": "ems", "region": "vac", "material": { "eps_r": 1.0 } }));
    call(
        &mut e_dgtd,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "vac",
            "kind": "dipole",
            "params": {
                "name": "ez_dip",
                "position": [0.0, 0.0, 0.0],
                "moment": [0.0, 0.0, 1.0e-9],
                "magnetic": false
            }
        }),
    );
    call(&mut e_dgtd, "solve_run", json!({}));

    let ff_dgtd = call(
        &mut e_dgtd,
        "analyse",
        json!({
            "analysis": "farfield",
            "physics": "ems",
            "params": {
                "surface": { "kind": "box", "lo": [-0.45, -0.45, -0.45], "hi": [0.45, 0.45, 0.45], "sampling": "grid" },
                "n_theta": 37,
                "n_phi": 24
            }
        }),
    );
    let d_max_dgtd = ff_dgtd["D_max"].as_f64().expect("DGTD NFFFT D_max");
    let th_dgtd = ff_dgtd["D_max_theta_deg"].as_f64().expect("DGTD NFFFT D_max_theta_deg");
    let p_ff_dgtd = ff_dgtd["P_rad_farfield_w"].as_f64().expect("DGTD P_rad_farfield_w");
    let p_flux_dgtd = ff_dgtd["P_rad_flux_w"].as_f64().expect("DGTD P_rad_flux_w");
    assert!(
        (d_max_dgtd - 1.50).abs() < 0.22
            && (th_dgtd - 90.0).abs() < 1e-6
            && p_ff_dgtd > 0.0
            && (p_ff_dgtd - p_flux_dgtd).abs() / p_flux_dgtd < 0.02,
        "DGTD -> Engine farfield expected equatorial Hertzian dipole & <2% Stratton-Chu power balance, got D_max={d_max_dgtd:.4}, theta={th_dgtd}, P_ff={p_ff_dgtd:e}, P_flux={p_flux_dgtd:e}"
    );
}

fn run_determinism_payload() -> String {
    let nl_toml = r#"
format = 2
name = "nl_det_p3"
kind = "crystalline"
point_group = "3m"
source = { cite = "Phase 3 multi-thread determinism test", verified = "2026-11-27" }

[optical]
wavelength_unit = "um"
range = [0.5, 2.5]
model = "constant"
n = 2.2
k = 0.0

[nonlinear.chi2]
wavelength = "1.55 um"
unit = "pm/V"
d33 = -25.0
d31 = -4.5
d22 = 2.0
d15 = -4.5

[nonlinear.qpm]
sign = 1.0
Lambda_qpm = "0.60 um"
axis = [0.0, 1.0, 0.0]
duty_cycle = 0.5
phase_offset = "0.0 um"

[nonlinear.chi3]
wavelength = "1.55 um"
n_2 = "2.5e-18 m^2/W"
beta_tpa = "2.0e-11 m/W"
sigma_fca = "2.0e-20 m^2"
k_fcd = "-1.35e-27 m^3"
tau_c = "30 fs"
"#;

    let f0 = C0 / 1.55e-6;
    let mut bits: Vec<u64> = Vec::new();

    // 1. Nonlinear DGTD (`chi2` + `qpm` + `chi3` + `TPA/FCA`) + `farfield` analysis
    {
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-det-dgtd", "dim": 3, "h": 0.25e-6 }));
        call(&mut e, "add_box", json!({ "name": "b", "min": [0.0, 0.0, 0.0], "max": [0.25e-6, 0.90e-6, 0.25e-6] }));
        call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b" }] }));
        call(&mut e, "add_material", json!({ "name": "nl_det_p3", "toml": nl_toml }));
        call(&mut e, "set_mesh", json!({ "h": 0.25e-6, "sweep": [{ "thickness": 0.90e-6, "n": 12 }] }));
        call(
            &mut e,
            "add_physics",
            json!({
                "physics": "ems",
                "study": "dgtd",
                "order": 2,
                "study_params": {
                    "t_final": 45.0e-15,
                    "cfl": 0.70,
                    "flux": "upwind",
                    "time_stepper": "lts",
                    "max_lts_tiers": 4,
                    "f0": f0,
                    "fwidth": 0.15 * f0,
                    "pulse_amplitude": 35.0,
                    "freqs_hz": [f0, 2.0 * f0],
                    "excited_port": 1,
                    "flux_monitors_on": "extreme(+y)",
                    "nonlinear": true
                }
            }),
        );
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "nl_det_p3" }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
        );
        call(&mut e, "solve_run", json!({}));
        let res = e.result("ems").expect("det dgtd").0;
        let dp = DrivenPortData::of(res).expect("det dgtd dp");
        for f_mat in &dp.s_matrix {
            for row in f_mat {
                for s in row {
                    bits.push(s[0].to_bits());
                    bits.push(s[1].to_bits());
                }
            }
        }
        let td = DgtdData::of(res).expect("det dgtd data");
        for &v in &td.flux_monitors[0].raw_flux_w_s2 {
            bits.push(v.to_bits());
        }
    }

    // 2. Nonlinear FDTD (`chi2` + `qpm` + `chi3` + `TPA/FCA`) + `farfield` & `overlap_beam`
    let ff_canon = {
        let cfg = json!({
            "courant": 0.95,
            "max_steps": 160,
            "shutoff_energy_ratio": 0.0,
            "subpixel": "polarized",
            "nonlinear": true,
            "grid": {
                "auto": { "wavelength_min_m": 0.775e-6, "min_steps_per_wave": 24.0 },
                "snap_planes": [[1, 0.0], [1, 0.90e-6]]
            },
            "boundaries": {
                "faces": [
                    { "kind": "pmc" },
                    { "kind": "pmc" },
                    { "kind": "cpml", "layers": 10, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-7 },
                    { "kind": "cpml", "layers": 10, "order": 3.5, "kappa_max": 4.0, "alpha_max": 0.05, "target_reflection": 1e-7 },
                    { "kind": "pec" },
                    { "kind": "pec" }
                ]
            },
            "sources": [{
                "kind": "mode_source",
                "name": "src",
                "axis": 1,
                "position_m": -0.20e-6,
                "direction": 1.0,
                "mode_index": 1,
                "waveform": {
                    "kind": "gaussian_pulse",
                    "f0_hz": f0,
                    "fwidth_hz": 0.15 * f0,
                    "amplitude": 2.0,
                    "phase_rad": 0.0,
                    "peak_offset_sigma": 4.0
                }
            }],
            "monitors": [
                {
                    "kind": "flux_monitor",
                    "name": "out",
                    "axis": 1,
                    "position_m": 0.90e-6,
                    "normal_sign": 1.0,
                    "freqs_hz": [f0, 2.0 * f0]
                },
                {
                    "kind": "field_monitor",
                    "name": "vol",
                    "freqs_hz": [f0]
                }
            ]
        });

        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "p3-det-fdtd", "dim": 3, "h": 0.25e-6 }));
        call(&mut e, "add_box", json!({ "name": "bi", "min": [0.0, -0.60e-6, 0.0], "max": [0.25e-6, 0.0, 0.25e-6] }));
        call(&mut e, "add_box", json!({ "name": "bc", "min": [0.0, 0.0, 0.0], "max": [0.25e-6, 0.90e-6, 0.25e-6] }));
        call(
            &mut e,
            "add_box",
            json!({ "name": "bo", "min": [0.0, 0.90e-6, 0.0], "max": [0.25e-6, 1.50e-6, 0.25e-6] }),
        );
        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "in_reg", "input": "bi" },
                    { "name": "core", "input": "bc" },
                    { "name": "out_reg", "input": "bo" }
                ]
            }),
        );
        call(&mut e, "add_material", json!({ "name": "nl_det_p3", "toml": nl_toml }));
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": 0.25e-6,
                "sweep": [
                    { "thickness": 0.60e-6, "n": 2 },
                    { "thickness": 0.90e-6, "n": 4 },
                    { "thickness": 0.60e-6, "n": 2 }
                ]
            }),
        );
        call(
            &mut e,
            "add_physics",
            json!({ "physics": "ems", "study": "fdtd", "order": 2, "study_params": { "config_json": cfg.to_string() } }),
        );
        for reg in ["in_reg", "out_reg"] {
            call(&mut e, "set_material", json!({ "physics": "ems", "region": reg, "material": { "eps_r": 4.84 } }));
        }
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": "nl_det_p3" }));
        call(&mut e, "solve_run", json!({}));
        let res = e.result("ems").expect("det fdtd").0;
        let md = FdtdMonitorData::of(res).expect("det fdtd md");
        for &v in &md.flux_monitors[0].raw_flux_w_s2 {
            bits.push(v.to_bits());
        }
        for cf in &res.solutions[0].cell_fields {
            for &v in &cf.re {
                bits.push(v.to_bits());
            }
            for &v in &cf.im {
                bits.push(v.to_bits());
            }
        }
        let mut ff = call(
            &mut e,
            "analyse",
            json!({
                "analysis": "farfield",
                "physics": "ems",
                "params": {
                    "surface": { "kind": "box", "lo": [0.04e-6, 0.15e-6, 0.04e-6], "hi": [0.21e-6, 0.75e-6, 0.21e-6], "sampling": "grid" },
                    "n_theta": 13,
                    "n_phi": 12
                }
            }),
        );
        if let Some(o) = ff.as_object_mut() {
            o.remove("age_s");
        }
        canonical_json(&ff).to_string()
    };

    canonical_json(&json!({ "bits": bits, "farfield": ff_canon })).to_string()
}

#[test]
fn phase3_subprocess_thread_determinism_1_vs_8() {
    if std::env::var("NYX_P3_DET_CHILD").as_deref() == Ok("1") {
        println!("NYX_P3_DET_PAYLOAD={}", run_determinism_payload());
        return;
    }

    let exe = std::env::current_exe().expect("current test executable");
    let run_with_threads = |nt: &str| -> String {
        let out = Command::new(&exe)
            .arg("--exact")
            .arg("phase3_subprocess_thread_determinism_1_vs_8")
            .arg("--nocapture")
            .env("NYX_P3_DET_CHILD", "1")
            .env("NYX_THREADS", nt)
            .output()
            .expect("spawn determinism child");
        assert!(
            out.status.success(),
            "determinism child NYX_THREADS={nt} failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_P3_DET_PAYLOAD="))
            .expect("missing NYX_P3_DET_PAYLOAD line")
            .to_string()
    };

    let p1 = run_with_threads("1");
    let p8 = run_with_threads("8");
    assert_eq!(
        p1, p8,
        "Phase 3 integrated nonlinear DGTD + FDTD + NFFFT outputs differed between NYX_THREADS=1 and NYX_THREADS=8"
    );
}
