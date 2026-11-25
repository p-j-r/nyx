// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! SQuADDS PR CI end-to-end `Engine::call` smoke suite for `nyx-engine` (Plan §8.3, `< 90 s` release budget):
//! - `S2.1 / C0.3`: Quarter-wave (`CLT`) hanger resonator real `CircuitRom` QEP `(f_r, kappa)` vs.
//!   2-port `ems::driven` Taubin circle fit (`s21_fit`) on the identical 3D `SimplexMesh`, plus S-matrix
//!   reciprocity and lossless unitarity invariants.
//! - `S2.3 + S3.2`: Unified `QubitCavity` equivalence across Method M1 (closed-port `craig_bampton`),
//!   Method M2 (passive-port `CircuitRom` QEP + BBQ), and Method M3 (3-port `ems::driven` Schur BBQ +
//!   exact Two-Pole Kubo `Y_Kubo^(g/e)(omega)` dispersive pull `chi = f_r^(e) - f_r^(g)` with TRK sum rule).
//! - `S3.4`: 3D CPW airbridge crossovers via `Op::Box` (`kind: "Hole"`) landing posts + elevated `Op::Sheet`
//!   bridge decks (`structured_pml: false`), verifying facet preservation and downward capacitive pull.
//! - `ADR 11`: Fresh-process `NYX_THREADS=1` vs. `NYX_THREADS=8` determinism (`|y_8 - y_1| / |y_1| < 1e-13`)
//!   across both 3D `dfs::electrostatic` (`far_field` + `I_edge`) and 3D `ems::cavity` (`CircuitRom`).

use femcore::linalg::driver::dense_eigh;
use femcore::linalg::scalar::{Scalar, C64};
use nyx_analysis::epr::{E_CHARGE, PLANCK_H};
use nyx_analysis::network_with_params;
use nyx_engine::Engine;
use nyx_physics::{CircuitRom, DrivenPortData, LabelResolver, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;
use std::sync::Arc;

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}): {err}"))
}

fn c_inv(z: C64) -> C64 {
    let d = (z.re * z.re + z.im * z.im).max(1.0e-300);
    C64::new(z.re / d, -z.im / d)
}

/// Checks `max |S_ij - S_ji|` and `max ||S^H S - I||_F` across all driven frequency points.
fn s_matrix_reciprocity_and_unitarity(dp: &DrivenPortData) -> (f64, f64) {
    let mut max_recip = 0.0_f64;
    let mut max_unit = 0.0_f64;
    for s_f in &dp.s_matrix {
        let p = s_f.len();
        for (i, row_i) in s_f.iter().enumerate().take(p) {
            for (j, row_j) in s_f.iter().enumerate().take(p).skip(i + 1) {
                let dr = row_i[j][0] - row_j[i][0];
                let di = row_i[j][1] - row_j[i][1];
                max_recip = max_recip.max(dr.hypot(di));
            }
        }
        let mut frob_sq = 0.0_f64;
        for i in 0..p {
            for j in 0..p {
                let mut acc = C64::ZERO;
                for row in s_f.iter().take(p) {
                    let ski = C64::new(row[i][0], -row[i][1]);
                    let skj = C64::new(row[j][0], row[j][1]);
                    acc += ski * skj;
                }
                if i == j {
                    acc = C64::new(acc.re - 1.0, acc.im);
                }
                frob_sq += acc.re * acc.re + acc.im * acc.im;
            }
        }
        max_unit = max_unit.max(frob_sq.sqrt());
    }
    (max_recip, max_unit)
}

/// Builds a high-Q coplanar-waveguide 3D quarter-wave (`CLT`) hanger resonator on a `2.4 x 3.8 mm` chip
/// (matching `C0.3` in `calibrate.rs`: `ground_plane` with feedline slot + `lambda/4` hanger trench, `Q_c >> 1000`),
/// optionally adding a floating 2-pad transmon (`pad_a`, `pad_b`, `jj_sheet`) in an open-end coupling pocket.
fn build_clt_chip(e: &mut Engine, with_qubit: bool) {
    let (wx, wy) = (2.4e-3, 3.8e-3);
    let (h_sub, h_top, h_bot) = (750.0e-6, 890.0e-6, 900.0e-6);
    call(e, "project_new", json!({ "title": "CLT C0.3 Smoke", "dim": 3, "h": 420.0e-6 }));
    call(
        e,
        "add_box",
        json!({ "name": "vac_bot", "min": [-wx / 2.0, -wy / 2.0, -h_sub - h_bot], "max": [wx / 2.0, wy / 2.0, -h_sub] }),
    );
    call(
        e,
        "add_box",
        json!({ "name": "si_box", "min": [-wx / 2.0, -wy / 2.0, -h_sub], "max": [wx / 2.0, wy / 2.0, 0.0] }),
    );
    call(
        e,
        "add_box",
        json!({ "name": "air_box", "min": [-wx / 2.0, -wy / 2.0, 0.0], "max": [wx / 2.0, wy / 2.0, h_top] }),
    );

    let mut gnd_holes = vec![
        vec![[-655.0e-6, -1850.0e-6], [-525.0e-6, -1850.0e-6], [-525.0e-6, 1850.0e-6], [-655.0e-6, 1850.0e-6]],
        vec![[-505.0e-6, -1450.0e-6], [-385.0e-6, -1450.0e-6], [-385.0e-6, 1630.0e-6], [-505.0e-6, 1630.0e-6]],
    ];
    if with_qubit {
        gnd_holes.push(vec![
            [-345.0e-6, 1260.0e-6],
            [55.0e-6, 1260.0e-6],
            [55.0e-6, 1780.0e-6],
            [-345.0e-6, 1780.0e-6],
        ]);
    }
    call(
        e,
        "add_polygons",
        json!({
            "name": "gnd_poly",
            "polygons": [[[-wx / 2.0, -wy / 2.0], [wx / 2.0, -wy / 2.0], [wx / 2.0, wy / 2.0], [-wx / 2.0, wy / 2.0]]],
            "holes": gnd_holes
        }),
    );
    call(e, "add_sheet", json!({ "name": "ground_plane", "input": "gnd_poly", "z": 0.0 }));

    // Through-feedline + two axial end-cap 50-Ohm port sheets touching ground_plane at y = +/-1.85 mm
    call(e, "add_rect", json!({ "name": "pin_2d", "min": [-620.0e-6, -1850.0e-6], "max": [-560.0e-6, -1750.0e-6] }));
    call(e, "add_sheet", json!({ "name": "port_in", "input": "pin_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "feed_2d", "min": [-620.0e-6, -1750.0e-6], "max": [-560.0e-6, 1750.0e-6] }));
    call(e, "add_sheet", json!({ "name": "feed_trace", "input": "feed_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "pout_2d", "min": [-620.0e-6, 1750.0e-6], "max": [-560.0e-6, 1850.0e-6] }));
    call(e, "add_sheet", json!({ "name": "port_out", "input": "pout_2d", "z": 0.0 }));

    // Quarter-wave CLT resonator shorted to ground_plane at y = -1.45 mm and open at y = +1.55 mm (L_res = 3.0 mm)
    call(e, "add_rect", json!({ "name": "res_2d", "min": [-475.0e-6, -1450.0e-6], "max": [-415.0e-6, 1550.0e-6] }));
    call(e, "add_sheet", json!({ "name": "clt_res", "input": "res_2d", "z": 0.0 }));

    if with_qubit {
        call(e, "add_rect", json!({ "name": "pa_2d", "min": [-305.0e-6, 1300.0e-6], "max": [15.0e-6, 1490.0e-6] }));
        call(e, "add_sheet", json!({ "name": "pad_a", "input": "pa_2d", "z": 0.0 }));
        call(e, "add_rect", json!({ "name": "pb_2d", "min": [-305.0e-6, 1550.0e-6], "max": [15.0e-6, 1740.0e-6] }));
        call(e, "add_sheet", json!({ "name": "pad_b", "input": "pb_2d", "z": 0.0 }));
        call(e, "add_rect", json!({ "name": "jj_2d", "min": [-170.0e-6, 1490.0e-6], "max": [-120.0e-6, 1550.0e-6] }));
        call(e, "add_sheet", json!({ "name": "jj_sheet", "input": "jj_2d", "z": 0.0 }));
    }

    call(
        e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "vac_bot", "input": "vac_bot" },
                { "name": "si", "input": "si_box" },
                { "name": "air", "input": "air_box" }
            ]
        }),
    );
    let excl = "!in_box([-650e-6, -1880e-6, -5e-6], [-530e-6, -1720e-6, 5e-6], touching) & !in_box([-650e-6, 1720e-6, -5e-6], [-530e-6, 1880e-6, 5e-6], touching) & !in_box([-200e-6, 1460e-6, -5e-6], [-90e-6, 1580e-6, 5e-6], touching)";
    call(e, "add_size", json!({ "sel": format!("clt_res:edge & {excl}"), "h": 140.0e-6 }));
    call(e, "add_size", json!({ "sel": format!("feed_trace:edge & {excl}"), "h": 160.0e-6 }));
    if with_qubit {
        call(e, "add_size", json!({ "sel": format!("(pad_a:edge | pad_b:edge) & {excl}"), "h": 85.0e-6 }));
    }
    call(
        e,
        "set_mesh",
        json!({
            "h": 450.0e-6,
            "dz_min": 140.0e-6,
            "max_z_slabs": 5,
            "z_growth_ratio": 2.2,
            "z_coarsen": true,
            "max_cells": 65_000,
            "structured_pml": false
        }),
    );
    call(e, "mesh_build", json!({ "wait": true }));
}

#[test]
fn clt_quarter_wave_qep_vs_driven_taubin_fit_smoke() {
    let mut e = Engine::new();
    build_clt_chip(&mut e, false);
    let mesh_arc = e.model().unwrap().mesh().clone();

    // 1. Real ems::cavity + passive lumped_port CircuitRom QEP
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 1,
            "study_params": { "n_modes": 2, "f_search": 8.8e9 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vac_bot", "material": { "eps_r": 1.0 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
    for sh in ["ground_plane", "feed_trace", "clt_res"] {
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": sh, "kind": "pec", "params": {} }));
    }
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port_in",
            "kind": "lumped_port",
            "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "complex_damping": false }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port_out",
            "kind": "lumped_port",
            "params": { "name": "port_out", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0], "complex_damping": false }
        }),
    );
    call(&mut e, "solve_run", json!({ "isolate": false }));
    assert!(Arc::ptr_eq(e.model().unwrap().mesh(), &mesh_arc));

    let net_qep = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": { "s21_ports": ["port_in", "port_out"], "resonator_branch": "hanger" }
        }),
    );
    let qep_modes = net_qep["qep_modes"].as_array().expect("qep_modes");
    let res_mode = qep_modes
        .iter()
        .filter(|m| {
            let f = m["f_re_Hz"].as_f64().unwrap_or(0.0);
            let k = m["kappa_Hz"].as_f64().unwrap_or(0.0);
            (2.0e9..14.0e9).contains(&f) && k > 5.0e3 && k < 0.2 * f
        })
        .min_by(|a, b| a["f_re_Hz"].as_f64().unwrap().total_cmp(&b["f_re_Hz"].as_f64().unwrap()))
        .unwrap_or_else(|| panic!("no fundamental CLT QEP mode in {net_qep}"));
    let f_qep = res_mode["f_re_Hz"].as_f64().unwrap();
    let kappa_qep = res_mode["kappa_Hz"].as_f64().unwrap();
    assert!(kappa_qep > 5.0e3 && kappa_qep < 50.0e6, "expected high-Q CLT kappa, got {kappa_qep:e}");

    // 2. 2-port ems::driven + Probst-Taubin circle fit on the exact same SimplexMesh
    let span = (5.0 * kappa_qep).clamp(0.5e6, 150.0e6);
    call(
        &mut e,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "driven",
            "order": 1,
            "params": {
                "f_min": f_qep - span,
                "f_max": f_qep + span,
                "n_points": 51,
                "prom_tol": 2.0e-4,
                "max_hdm_solves": 6
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "isolate": false }));
    assert!(Arc::ptr_eq(e.model().unwrap().mesh(), &mesh_arc), "driven solve must reuse identical mesh");

    let dp = DrivenPortData::of(e.result("ems").unwrap().0).expect("DrivenPortData present");
    assert!(dp.max_residual < 2.0e-4, "PROM residual {:.3e} >= 2e-4", dp.max_residual);
    let (recip_err, unit_err) = s_matrix_reciprocity_and_unitarity(&dp);
    assert!(recip_err < 1.0e-10, "driven S reciprocity error {recip_err:e} >= 1e-10");
    assert!(unit_err < 1.0e-8, "lossless driven S unitarity error {unit_err:e} >= 1e-8");

    let net_drv = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": { "s21_ports": ["port_in", "port_out"], "resonator_branch": "hanger" }
        }),
    );
    let fit = &net_drv["s21_fit"];
    let f_drv = fit["f_r_Hz"].as_f64().unwrap();
    let kappa_drv = fit["kappa_Hz"].as_f64().unwrap();

    let rel_f = (f_drv - f_qep).abs() / f_qep;
    let rel_kappa = (kappa_drv - kappa_qep).abs() / kappa_qep;
    assert!(rel_f <= 1.0e-3, "CLT QEP vs Driven f_r mismatch {:.4}% > 0.1% ({f_qep:e} vs {f_drv:e})", rel_f * 100.0);
    assert!(
        rel_kappa <= 0.02,
        "CLT QEP vs Driven kappa mismatch {:.2}% > 2.0% ({kappa_qep:e} vs {kappa_drv:e})",
        rel_kappa * 100.0
    );
}

/// Exact 41-state charge-basis transmon Kubo admittance `Y_Kubo^(g/e)(omega)` (`Plan Appendix A.2`,
/// with exact `16 E_C` oscillator-strength prefactor satisfying the Thomas–Reiche–Kuhn sum rule).
struct KuboTransmon {
    c_env_f: f64,
    c_sigma_f: f64,
    f01_bare_hz: f64,
    alpha_bare_hz: f64,
    poles_g: Vec<(f64, f64)>,
    poles_e_up: Vec<(f64, f64)>,
    omega_01: f64,
    w_01: f64,
    trk_sum_g: f64,
    trk_sum_e: f64,
}

impl KuboTransmon {
    fn new(c_env_f: f64, c_j_f: f64, e_j_hz: f64) -> Self {
        let c_sigma_f = c_env_f + c_j_f;
        let e_c_hz = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * c_sigma_f);
        let n_max: isize = 20;
        let dim = (2 * n_max + 1) as usize;
        let charges: Vec<f64> = (-n_max..=n_max).map(|n| n as f64).collect();
        let mut h_flat = vec![0.0_f64; dim * dim];
        for i in 0..dim {
            h_flat[i * dim + i] = 4.0 * e_c_hz * charges[i] * charges[i];
            if i + 1 < dim {
                h_flat[i * dim + (i + 1)] = -0.5 * e_j_hz;
                h_flat[(i + 1) * dim + i] = -0.5 * e_j_hz;
            }
        }
        let (evals, evecs_col) = dense_eigh(&h_flat, dim).expect("dense_eigh transmon charge basis");
        // dense_eigh returns eigenvectors as Vec<Vec<f64>> indexed [col][row]
        let n_elem =
            |k: usize, m: usize| -> f64 { (0..dim).map(|r| charges[r] * evecs_col[k][r] * evecs_col[m][r]).sum() };

        let f01 = evals[1] - evals[0];
        let f12 = evals[2] - evals[1];
        let n01 = n_elem(0, 1);
        let w_01 = (16.0 * e_c_hz * n01 * n01) / f01;

        let mut poles_g = Vec::new();
        let mut trk_sum_g = 0.0;
        for m in 1..16 {
            let f0m = evals[m] - evals[0];
            let n0m = n_elem(0, m);
            let w0m = (16.0 * e_c_hz * n0m * n0m) / f0m;
            trk_sum_g += w0m;
            if w0m > 1.0e-15 {
                poles_g.push((2.0 * PI * f0m, w0m));
            }
        }

        let mut poles_e_up = Vec::new();
        let mut trk_up = 0.0;
        for m in 2..17 {
            let f1m = evals[m] - evals[1];
            let n1m = n_elem(1, m);
            let w1m = (16.0 * e_c_hz * n1m * n1m) / f1m;
            trk_up += w1m;
            if w1m > 1.0e-15 {
                poles_e_up.push((2.0 * PI * f1m, w1m));
            }
        }

        Self {
            c_env_f,
            c_sigma_f,
            f01_bare_hz: f01,
            alpha_bare_hz: f12 - f01,
            poles_g,
            poles_e_up,
            omega_01: 2.0 * PI * f01,
            w_01,
            trk_sum_g,
            trk_sum_e: trk_up - w_01,
        }
    }

    fn y_im(&self, freq_hz: f64, excited: bool) -> f64 {
        let w = 2.0 * PI * freq_hz;
        let w2 = w * w;
        let sum_term = if !excited {
            self.poles_g.iter().map(|&(wm, wt)| wt / (wm * wm - w2)).sum::<f64>()
        } else {
            self.poles_e_up.iter().map(|&(wm, wt)| wt / (wm * wm - w2)).sum::<f64>()
                - self.w_01 / (self.omega_01 * self.omega_01 - w2)
        };
        -self.c_sigma_f / (w * sum_term) - w * self.c_env_f
    }
}

#[test]
fn unified_qubit_cavity_m1_m2_m3_equivalence_smoke() {
    let mut e = Engine::new();
    build_clt_chip(&mut e, true);
    let mesh_arc = e.model().unwrap().mesh().clone();

    let e_j_hz = 13.5e9;
    let c_j_f = 45.0e-15;

    // --- Method M1: Closed-port ems::cavity + craig_bampton ---
    let pec_sheets = ["ground_plane", "feed_trace", "clt_res", "pad_a", "pad_b"];
    let mut bcs_m1: Vec<Value> = pec_sheets.iter().map(|s| json!({ "on": s, "kind": "pec", "params": {} })).collect();
    bcs_m1.push(json!({
        "on": "jj_sheet",
        "kind": "lumped_jj",
        "params": {
            "name": "Q1_JJ",
            "E_j_Hz": e_j_hz,
            "C_j": c_j_f,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }
    }));
    call(
        &mut e,
        "set_problem",
        json!({
            "physics": "ems",
            "problem": {
                "format": 1,
                "physics": "ems",
                "order": 1,
                "study": { "kind": "cavity", "params": { "n_modes": 3, "f_search": 6.5e9 } },
                "materials": { "vac_bot": { "eps_r": 1.0 }, "si": { "eps_r": 11.45 }, "air": { "eps_r": 1.0 } },
                "bcs": bcs_m1,
                "params": {}
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "isolate": false }));
    assert!(Arc::ptr_eq(e.model().unwrap().mesh(), &mesh_arc));
    let ham_m1 = call(
        &mut e,
        "analyse",
        json!({ "analysis": "hamiltonian", "physics": "ems", "params": { "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [5, 4] } }),
    );
    let f01_m1 = ham_m1["f01_Hz"][0].as_f64().unwrap();
    let alpha_m1 = ham_m1["alpha_Hz"][0].as_f64().unwrap();
    let chi_m1 = ham_m1["chi_matrix_Hz"][0][1].as_f64().unwrap();

    // --- Method M2: Passive-port CircuitRom QEP + craig_bampton + BBQ on identical mesh ---
    let mut bcs_m2 = bcs_m1.clone();
    bcs_m2.push(json!({
        "on": "port_in",
        "kind": "lumped_port",
        "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "complex_damping": false }
    }));
    bcs_m2.push(json!({
        "on": "port_out",
        "kind": "lumped_port",
        "params": { "name": "port_out", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0], "complex_damping": false }
    }));
    call(
        &mut e,
        "set_problem",
        json!({
            "physics": "ems",
            "problem": {
                "format": 1,
                "physics": "ems",
                "order": 1,
                "study": { "kind": "cavity", "params": { "n_modes": 3, "f_search": 6.5e9 } },
                "materials": { "vac_bot": { "eps_r": 1.0 }, "si": { "eps_r": 11.45 }, "air": { "eps_r": 1.0 } },
                "bcs": bcs_m2,
                "params": {}
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "isolate": false }));
    assert!(Arc::ptr_eq(e.model().unwrap().mesh(), &mesh_arc));
    let c_env_rom = CircuitRom::of(e.result("ems").unwrap().0).map(|rom| rom.c_port_f[0][0]);
    let ham_m2 = call(
        &mut e,
        "analyse",
        json!({ "analysis": "hamiltonian", "physics": "ems", "params": { "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [5, 4] } }),
    );
    let net_m2 = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": {
                "active_ports": ["Q1_JJ"],
                "terminations_ohm": { "port_in": 50.0, "port_out": 50.0 },
                "s21_ports": ["port_in", "port_out"],
                "resonator_branch": "hanger"
            }
        }),
    );
    let f01_m2 = ham_m2["f01_Hz"][0].as_f64().unwrap();
    let alpha_m2 = ham_m2["alpha_Hz"][0].as_f64().unwrap();
    let chi_m2 = ham_m2["chi_matrix_Hz"][0][1].as_f64().unwrap();

    // Find resonator QEP mode near ham_m2["f_cav_Hz"][0]
    let f_cav_bare = ham_m2["f_cav_Hz"][0].as_f64().unwrap();
    let qep_res = net_m2["qep_modes"]
        .as_array()
        .unwrap()
        .iter()
        .min_by(|a, b| {
            let da = (a["f_re_Hz"].as_f64().unwrap() - f_cav_bare).abs();
            let db = (b["f_re_Hz"].as_f64().unwrap() - f_cav_bare).abs();
            da.partial_cmp(&db).unwrap()
        })
        .unwrap();
    let f_r_m2 = qep_res["f_re_Hz"].as_f64().unwrap();
    let kappa_m2 = qep_res["kappa_Hz"].as_f64().unwrap().max(10.0e3);

    assert!((f01_m2 - f01_m1).abs() / f01_m1 <= 0.005, "M1 vs M2 f01 drift: {f01_m1:e} vs {f01_m2:e}");
    assert!((alpha_m2 - alpha_m1).abs() / alpha_m1.abs() <= 0.02, "M1 vs M2 alpha drift: {alpha_m1:e} vs {alpha_m2:e}");
    assert!((chi_m2 - chi_m1).abs() / chi_m1.abs().max(25.0e3) <= 0.05, "M1 vs M2 chi drift: {chi_m1:e} vs {chi_m2:e}");

    // --- Method M3: 3-Port ems::driven Schur + Exact Two-Pole Kubo Y_Kubo^(g/e) on identical mesh ---
    // Include tight local windows around `f_r^(g)` and `f_r^(e) = f_r + chi` (matching `em_hw.rs:815-827`)
    let mut freqs_m3 = vec![1.0e9, 1.8e9];
    let kp = kappa_m2.max(10.0e3);
    for center in [f_r_m2, f_r_m2 + chi_m2] {
        for k in 0..=24 {
            freqs_m3.push((center - 6.0 * kp) + 12.0 * kp * (k as f64) / 24.0);
        }
    }
    freqs_m3.sort_by(|a, b| a.total_cmp(b));
    let mut dedup_m3: Vec<f64> = Vec::with_capacity(freqs_m3.len());
    for f in freqs_m3 {
        if dedup_m3.last().is_none_or(|&prev| f - prev > (0.02 * kp).clamp(50.0, 500.0)) {
            dedup_m3.push(f);
        }
    }
    let mut bcs_m3: Vec<Value> = pec_sheets.iter().map(|s| json!({ "on": s, "kind": "pec", "params": {} })).collect();
    bcs_m3.push(json!({
        "on": "port_in",
        "kind": "lumped_port",
        "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
    }));
    bcs_m3.push(json!({
        "on": "port_out",
        "kind": "lumped_port",
        "params": { "name": "port_out", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }
    }));
    bcs_m3.push(json!({
        "on": "jj_sheet",
        "kind": "lumped_port",
        "params": { "name": "Q1_JJ", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
    }));
    call(
        &mut e,
        "set_problem",
        json!({
            "physics": "ems",
            "problem": {
                "format": 1,
                "physics": "ems",
                "order": 1,
                "study": {
                    "kind": "driven",
                    "params": {
                        "freqs_hz": dedup_m3,
                        "prom_tol": 2.0e-4,
                        "max_hdm_solves": 8,
                        "seed_freqs_hz": [f01_m2, f_r_m2]
                    }
                },
                "materials": { "vac_bot": { "eps_r": 1.0 }, "si": { "eps_r": 11.45 }, "air": { "eps_r": 1.0 } },
                "bcs": bcs_m3,
                "params": {}
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "isolate": false }));
    assert!(Arc::ptr_eq(e.model().unwrap().mesh(), &mesh_arc));

    let dp3 = DrivenPortData::of(e.result("ems").unwrap().0).unwrap();
    let idx_in = dp3.ports.iter().position(|p| p == "port_in").unwrap();
    let idx_out = dp3.ports.iter().position(|p| p == "port_out").unwrap();
    let idx_jj = dp3.ports.iter().position(|p| p == "Q1_JJ").unwrap();

    // Use exact `CircuitRom::c_port_f[0][0]` (matching `em_hw.rs:1535`), falling back to low-f Y_env extrapolation
    let w1 = 2.0 * PI * dp3.freqs_hz[0];
    let w2 = 2.0 * PI * dp3.freqs_hz[1];
    let c1 = dp3.y_matrix[0][idx_jj][idx_jj][1] / w1;
    let c2 = dp3.y_matrix[1][idx_jj][idx_jj][1] / w2;
    let c_env_lowf = (w2 * w2 * c1 - w1 * w1 * c2) / (w2 * w2 - w1 * w1);
    let c_env = c_env_rom.unwrap_or(c_env_lowf);
    assert!(c_env > 1.0e-15, "expected positive C_env, got {c_env:e}");

    let kubo = KuboTransmon::new(c_env, c_j_f, e_j_hz);
    assert!((kubo.trk_sum_g - 1.0).abs() < 1.0e-10, "TRK sum (g) violated: {}", kubo.trk_sum_g);
    assert!((kubo.trk_sum_e - 1.0).abs() < 1.0e-10, "TRK sum (e) violated: {}", kubo.trk_sum_e);

    // Synthesize 2-port feedline S-parameters over separate local ground/excited dip windows
    let synth_band = |center_hz: f64, excited: bool| -> (f64, f64) {
        let y0 = 1.0 / 50.0;
        let half_w = (8.0 * kp).max(40.0e3);
        let idxs: Vec<usize> = dp3
            .freqs_hz
            .iter()
            .enumerate()
            .filter_map(|(k, &f)| (k >= 2 && (f - center_hz).abs() <= half_w).then_some(k))
            .collect();
        assert!(idxs.len() >= 9, "insufficient points near {center_hz:e} (len={})", idxs.len());
        let band_freqs: Vec<f64> = idxs.iter().map(|&k| dp3.freqs_hz[k]).collect();
        let nf = band_freqs.len();
        let mut s_mat = Vec::with_capacity(nf);
        let mut y_mat = Vec::with_capacity(nf);
        let mut z_mat = Vec::with_capacity(nf);
        for &k in &idxs {
            let f_hz = dp3.freqs_hz[k];
            let y3 = &dp3.y_matrix[k];
            let y_jj = C64::new(y3[idx_jj][idx_jj][0], y3[idx_jj][idx_jj][1] + kubo.y_im(f_hz, excited));
            let inv_yjj = c_inv(y_jj);
            let feed = [idx_in, idx_out];
            let mut y_ff = [[C64::ZERO; 2]; 2];
            for a in 0..2 {
                let y_aj = C64::new(y3[feed[a]][idx_jj][0], y3[feed[a]][idx_jj][1]);
                for b in 0..2 {
                    let y_jb = C64::new(y3[idx_jj][feed[b]][0], y3[idx_jj][feed[b]][1]);
                    let y_ab = C64::new(y3[feed[a]][feed[b]][0], y3[feed[a]][feed[b]][1]);
                    y_ff[a][b] = y_ab - y_aj * inv_yjj * y_jb;
                }
            }
            let a00 = C64::new(y0 + y_ff[0][0].re, y_ff[0][0].im);
            let a01 = y_ff[0][1];
            let a10 = y_ff[1][0];
            let a11 = C64::new(y0 + y_ff[1][1].re, y_ff[1][1].im);
            let inv_det = c_inv(a00 * a11 - a01 * a10);
            let inv_a = [[a11 * inv_det, a01.scale(-1.0) * inv_det], [a10.scale(-1.0) * inv_det, a00 * inv_det]];
            let n = [
                [C64::new(y0 - y_ff[0][0].re, -y_ff[0][0].im), y_ff[0][1].scale(-1.0)],
                [y_ff[1][0].scale(-1.0), C64::new(y0 - y_ff[1][1].re, -y_ff[1][1].im)],
            ];
            let mut s2 = [[0.0_f64; 2]; 4];
            for r in 0..2 {
                for c in 0..2 {
                    let v = n[r][0] * inv_a[0][c] + n[r][1] * inv_a[1][c];
                    s2[r * 2 + c] = [v.re, v.im];
                }
            }
            let inv_dy = c_inv(y_ff[0][0] * y_ff[1][1] - y_ff[0][1] * y_ff[1][0]);
            let z00 = y_ff[1][1] * inv_dy;
            let z01 = y_ff[0][1].scale(-1.0) * inv_dy;
            let z10 = y_ff[1][0].scale(-1.0) * inv_dy;
            let z11 = y_ff[0][0] * inv_dy;
            s_mat.push(vec![vec![s2[0], s2[1]], vec![s2[2], s2[3]]]);
            y_mat.push(vec![
                vec![[y_ff[0][0].re, y_ff[0][0].im], [y_ff[0][1].re, y_ff[0][1].im]],
                vec![[y_ff[1][0].re, y_ff[1][0].im], [y_ff[1][1].re, y_ff[1][1].im]],
            ]);
            z_mat.push(vec![vec![[z00.re, z00.im], [z01.re, z01.im]], vec![[z10.re, z10.im], [z11.re, z11.im]]]);
        }
        let dp2 = DrivenPortData {
            freqs_hz: band_freqs,
            ports: vec!["port_in".into(), "port_out".into()],
            z0_ohm: vec![vec![[50.0, 0.0], [50.0, 0.0]]; nf],
            s_matrix: s_mat,
            z_matrix: z_mat,
            y_matrix: y_mat,
            hdm_solves: dp3.hdm_solves,
            rom_dim: dp3.rom_dim,
            anchor_freqs_hz: dp3.anchor_freqs_hz.clone(),
            max_residual: dp3.max_residual,
        };
        let sr = SolveResult { scalars: vec![dp2.scalar()], ..Default::default() };
        let fit =
            network_with_params(&sr, &json!({ "s21_ports": ["port_in", "port_out"], "resonator_branch": "hanger" }))
                .unwrap()
                .s21_fit
                .unwrap();
        (fit.f_r_hz, fit.kappa_hz)
    };

    let (fr_g, _kappa_g) = synth_band(f_r_m2, false);
    let (fr_e, _kappa_e) = synth_band(f_r_m2 + chi_m2, true);
    let chi_m3 = fr_e - fr_g;

    assert!(
        (kubo.f01_bare_hz - f01_m2).abs() / f01_m2 <= 0.02,
        "M3 vs M2 f01 mismatch: {:e} vs {f01_m2:e}",
        kubo.f01_bare_hz
    );
    assert!(
        (kubo.alpha_bare_hz - alpha_m2).abs() / alpha_m2.abs() <= 0.05,
        "M3 vs M2 alpha mismatch: {:e} vs {alpha_m2:e}",
        kubo.alpha_bare_hz
    );
    assert!((fr_g - f_r_m2).abs() / f_r_m2 <= 0.002, "M3 vs M2 f_r mismatch: {fr_g:e} vs {f_r_m2:e}");
    assert!(chi_m3 * chi_m2 > 0.0, "M3 Kubo chi ({chi_m3:e}) and M2 CB chi ({chi_m2:e}) must have identical sign");
    let rel_chi = (chi_m3 - chi_m2).abs() / chi_m2.abs().max(25.0e3);
    assert!(
        rel_chi <= 0.15,
        "M3 Two-Pole Kubo chi ({chi_m3:e}) vs M2 CB chi ({chi_m2:e}): rel={:.2}%",
        rel_chi * 100.0
    );
}

#[test]
fn airbridge_hole_post_and_deck_sheet_facet_preservation_smoke() {
    let solve_cpw_with_bridges = |with_bridges: bool| -> (f64, Engine) {
        let mut e = Engine::new();
        call(&mut e, "project_new", json!({ "title": "S3.4 Airbridge Smoke", "dim": 3, "h": 180.0e-6 }));
        let (w_box, l_box, h_sub, h_air) = (500.0e-6, 2.4e-3, 200.0e-6, 200.0e-6);
        call(
            &mut e,
            "add_box",
            json!({ "name": "si_box", "min": [-w_box / 2.0, -l_box / 2.0, -h_sub], "max": [w_box / 2.0, l_box / 2.0, 0.0] }),
        );
        call(
            &mut e,
            "add_box",
            json!({ "name": "air_box", "min": [-w_box / 2.0, -l_box / 2.0, 0.0], "max": [w_box / 2.0, l_box / 2.0, h_air] }),
        );

        // Open-ended lambda/2 CPW center trace ([-20, 20] um x [-900, 900] um) + flanking ground planes at |x| >= 45 um
        call(&mut e, "add_rect", json!({ "name": "sig_2d", "min": [-20.0e-6, -900.0e-6], "max": [20.0e-6, 900.0e-6] }));
        call(&mut e, "add_sheet", json!({ "name": "cpw_sig", "input": "sig_2d", "z": 0.0 }));
        call(&mut e, "add_rect", json!({ "name": "gl_2d", "min": [-250.0e-6, -1.2e-3], "max": [-45.0e-6, 1.2e-3] }));
        call(&mut e, "add_sheet", json!({ "name": "gnd_l", "input": "gl_2d", "z": 0.0 }));
        call(&mut e, "add_rect", json!({ "name": "gr_2d", "min": [45.0e-6, -1.2e-3], "max": [250.0e-6, 1.2e-3] }));
        call(&mut e, "add_sheet", json!({ "name": "gnd_r", "input": "gr_2d", "z": 0.0 }));

        let h_bridge = 4.0e-6;
        let bridge_ys = [-350.0e-6, 350.0e-6];
        let mut regions = vec![
            json!({ "name": "si", "input": "si_box", "kind": "Material" }),
            json!({ "name": "air", "input": "air_box", "kind": "Material" }),
        ];

        if with_bridges {
            for (idx, &yc) in bridge_ys.iter().enumerate() {
                let pl = format!("post_l_{idx}");
                let pr = format!("post_r_{idx}");
                let bd = format!("bridge_{idx}");
                let bd_2d = format!("bridge_2d_{idx}");
                call(
                    &mut e,
                    "add_box",
                    json!({ "name": pl, "min": [-95.0e-6, yc - 20.0e-6, 0.0], "max": [-55.0e-6, yc + 20.0e-6, h_bridge] }),
                );
                call(
                    &mut e,
                    "add_box",
                    json!({ "name": pr, "min": [55.0e-6, yc - 20.0e-6, 0.0], "max": [95.0e-6, yc + 20.0e-6, h_bridge] }),
                );
                call(
                    &mut e,
                    "add_rect",
                    json!({ "name": bd_2d, "min": [-95.0e-6, yc - 20.0e-6], "max": [95.0e-6, yc + 20.0e-6] }),
                );
                call(&mut e, "add_sheet", json!({ "name": bd, "input": bd_2d, "z": h_bridge }));
                regions.push(json!({ "name": pl, "input": pl, "kind": "Hole" }));
                regions.push(json!({ "name": pr, "input": pr, "kind": "Hole" }));
            }
        }

        call(&mut e, "add_overlay", json!({ "regions": regions }));
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": 180.0e-6,
                "dz_min": 4.0e-6,
                "max_z_slabs": 12,
                "z_growth_ratio": 1.75,
                "structured_pml": false
            }),
        );

        call(
            &mut e,
            "add_physics",
            json!({
                "physics": "ems",
                "study": "cavity",
                "order": 1,
                "study_params": { "n_modes": 2, "f_search": 30.0e9 }
            }),
        );
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
        call(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
        for s in ["cpw_sig", "gnd_l", "gnd_r"] {
            call(&mut e, "add_condition", json!({ "physics": "ems", "on": s, "kind": "pec", "params": {} }));
        }
        if with_bridges {
            for idx in 0..bridge_ys.len() {
                let bd = format!("bridge_{idx}");
                call(&mut e, "add_condition", json!({ "physics": "ems", "on": bd, "kind": "pec", "params": {} }));
                for side in ["l", "r"] {
                    let p = format!("post_{side}_{idx}");
                    let wall_sel = format!("(name({p}.xmin) | name({p}.xmax) | name({p}.ymin) | name({p}.ymax))");
                    call(
                        &mut e,
                        "add_condition",
                        json!({ "physics": "ems", "on": wall_sel, "kind": "pec", "params": {} }),
                    );
                }
            }
        }

        call(&mut e, "solve_run", json!({ "isolate": false }));
        let (res, _) = e.result("ems").unwrap();
        let ev0 = res.solutions[0].eigenvalue.expect("cavity eigenvalue");
        assert_eq!(ev0[1], 0.0, "PEC airbridge cavity must have purely real eigenfrequency");
        (ev0[0], e)
    };

    let (f_nb0, _) = solve_cpw_with_bridges(false);
    let (f_nb2, e_br) = solve_cpw_with_bridges(true);

    // Verify every bridge deck and Hole post lateral boundary retained physical facets in the 3D SimplexMesh
    let model = e_br.model().unwrap();
    let sm = model.solver_mesh();
    for idx in 0..2 {
        let deck_labs = model.boundaries(&format!("bridge_{idx}")).expect("bridge deck boundary labels");
        let deck_facets = sm.facet_labels.iter().filter(|l| deck_labs.contains(l)).count();
        assert!(deck_facets >= 2, "bridge_{idx} deck lost mesh facets (count={deck_facets})");

        for side in ["l", "r"] {
            let p = format!("post_{side}_{idx}");
            let wall_sel = format!("(name({p}.xmin) | name({p}.xmax) | name({p}.ymin) | name({p}.ymax))");
            let post_labs = model.boundaries(&wall_sel).expect("Hole post lateral wall labels");
            let post_facets = sm.facet_labels.iter().filter(|l| post_labs.contains(l)).count();
            assert!(post_facets >= 8, "{p} lateral walls lost mesh facets (count={post_facets})");
        }
    }

    assert!(
        f_nb2 < f_nb0,
        "3D capacitive airbridge crossovers must pull CPW resonance downward: f(Nb=2)={f_nb2:e} >= f(Nb=0)={f_nb0:e}"
    );
}

#[test]
fn nyx_threads_1_vs_8_determinism_smoke() {
    if std::env::var("NYX_SMOKE_SUBPROC_DFS_DET").as_deref() == Ok("1") {
        // 3D P2 dfs::electrostatic with AF-3 far_field + FR-11 200 nm thickness correction
        let mut e_dfs = Engine::new();
        call(&mut e_dfs, "project_new", json!({ "title": "Det DFS", "dim": 3, "h": 120.0e-6 }));
        call(
            &mut e_dfs,
            "add_rect",
            json!({ "name": "p1_2d", "min": [-90.0e-6, -50.0e-6], "max": [-12.0e-6, 50.0e-6] }),
        );
        call(&mut e_dfs, "add_rect", json!({ "name": "p2_2d", "min": [12.0e-6, -50.0e-6], "max": [90.0e-6, 50.0e-6] }));
        call(&mut e_dfs, "add_sheet", json!({ "name": "pad1", "input": "p1_2d", "z": 0.0 }));
        call(&mut e_dfs, "add_sheet", json!({ "name": "pad2", "input": "p2_2d", "z": 0.0 }));
        call(
            &mut e_dfs,
            "add_chip_stack",
            json!({ "lateral_pad": 120.0e-6, "t_sub": 120.0e-6, "h_vac_top": 120.0e-6, "h_vac_bot": 0.0 }),
        );
        call(&mut e_dfs, "add_size", json!({ "sel": "pad1:edge", "h": 12.0e-6 }));
        call(&mut e_dfs, "add_size", json!({ "sel": "pad2:edge", "h": 12.0e-6 }));
        call(&mut e_dfs, "set_mesh", json!({ "h": 85.0e-6, "dz_min": 12.0e-6, "max_z_slabs": 10 }));
        call(&mut e_dfs, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
        call(
            &mut e_dfs,
            "set_physics_params",
            json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1.0e-11 } }),
        );
        call(&mut e_dfs, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } }));
        call(&mut e_dfs, "set_material", json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } }));
        call(
            &mut e_dfs,
            "add_condition",
            json!({ "physics": "dfs", "on": "pad1", "kind": "terminal", "params": { "name": "c1" } }),
        );
        call(
            &mut e_dfs,
            "add_condition",
            json!({ "physics": "dfs", "on": "pad2", "kind": "terminal", "params": { "name": "c2" } }),
        );
        call(
            &mut e_dfs,
            "add_condition",
            json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2 } }),
        );
        call(&mut e_dfs, "solve_run", json!({ "isolate": false }));
        let cap = call(
            &mut e_dfs,
            "analyse",
            json!({ "analysis": "capacitance_matrix", "physics": "dfs", "params": { "metal_thickness_m": 200.0e-9 } }),
        );
        let c00 = cap["maxwell_matrix_F"][0][0].as_f64().unwrap();
        let c01 = cap["maxwell_matrix_F"][0][1].as_f64().unwrap();
        println!("NYX_SQUADDS_DFS_DET:{}", json!([c00, c01]));
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |threads: &str| -> Vec<f64> {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("nyx_threads_1_vs_8_determinism_smoke")
            .arg("--nocapture")
            .env("NYX_THREADS", threads)
            .env("NYX_SMOKE_SUBPROC_DFS_DET", "1")
            .output()
            .expect("spawn dfs determinism subprocess");
        assert!(out.status.success(), "NYX_THREADS={threads} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let payload = stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_SQUADDS_DFS_DET:"))
            .unwrap_or_else(|| panic!("missing NYX_SQUADDS_DFS_DET line in stdout:\n{stdout}"));
        serde_json::from_str(payload).unwrap()
    };

    let y1 = run_with_threads("1");
    let y8 = run_with_threads("8");
    let names = ["C_00", "C_01"];
    for ((&v1, &v8), name) in y1.iter().zip(&y8).zip(&names) {
        let rel = (v8 - v1).abs() / v1.abs().max(1.0e-30);
        assert!(
            rel < 1.0e-13,
            "ADR 11 DFS thread-determinism drift on {name}: NYX_THREADS=1 ({v1:.16e}) vs NYX_THREADS=8 ({v8:.16e}), rel={rel:e}"
        );
    }
}

/// Core Bug 2 Regression (`crates/nyx-ems/src/cavity.rs:646-694`):
/// When a `rank1` `lumped_port` or `lumped_jj` sheet has subdivided interior or lateral-vacuum vertices,
/// `cavity.rs:648-652, 676-681` unconditionally inserts ALL vertices of `port.triangles` into `cond_vertices`,
/// dropping their zero-frequency gradient columns `(G^T b_j)_v == 0` from `g_red` and creating spurious
/// near-zero modes that stall `krylov_schur`. Expected to FAIL on unmodified `HEAD` (`9480093a`).
#[test]
fn bug2_cavity_refined_lumped_port_stalls() {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "Bug2 Refined Port Cavity", "dim": 3, "h": 280.0e-6 }));
    call(
        &mut e,
        "add_box",
        json!({ "name": "si", "min": [-500.0e-6, -900.0e-6, -240.0e-6], "max": [500.0e-6, 900.0e-6, 0.0] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "air", "min": [-500.0e-6, -900.0e-6, 0.0], "max": [500.0e-6, 900.0e-6, 240.0e-6] }),
    );
    call(&mut e, "add_rect", json!({ "name": "feed_r", "min": [-80.0e-6, -650.0e-6], "max": [80.0e-6, 650.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "feed", "input": "feed_r", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "port_r", "min": [-80.0e-6, -900.0e-6], "max": [80.0e-6, -650.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "port_in", "input": "port_r", "z": 0.0 }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "si", "input": "si" }, { "name": "air", "input": "air" }] }),
    );
    // Explicitly refine the lumped_port rectangle edges so it acquires interior + lateral-vacuum vertices
    call(&mut e, "add_size", json!({ "name": "sz_port", "sel": "port_in:edge", "h": 40.0e-6 }));
    call(
        &mut e,
        "set_mesh",
        json!({ "h": 280.0e-6, "dz_min": 45.0e-6, "max_z_slabs": 8, "z_growth_ratio": 1.6, "max_cells": 60_000 }),
    );
    call(
        &mut e,
        "add_physics",
        json!({ "physics": "ems", "study": "cavity", "order": 1, "study_params": { "n_modes": 2, "f_search": 8.0e9 } }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "feed", "kind": "pec", "params": {} }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port_in",
            "kind": "lumped_port",
            "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
        }),
    );

    let reply = e
        .call("solve_run", &json!({ "physics": "ems", "isolate": false }))
        .expect("Core Bug 2 (cavity.rs:646-694): refined lumped_port must not stall krylov_schur");
    let (res, _) = e.result("ems").expect("ems result missing");
    assert!(
        res.solutions.iter().all(|s| s.eigenvalue.is_some_and(|ev| ev[0] > 1.0e8)),
        "Core Bug 2 (cavity.rs:646-694): refined lumped_port emitted spurious near-zero mode(s): {:?} ({reply})",
        res.solutions.iter().map(|s| s.eigenvalue).collect::<Vec<_>>()
    );
}

/// Core Bug 3 Regression (`crates/femcore/src/meshing/sheet3d.rs` `z_coarsen: true`):
/// On a multi-layer finite-die electrostatic stack (`vac_box` enclosing confined `sub_box`, `P = 1.2 L`),
/// enabling `z_coarsen: true` either collapses exterior cap boundary labels (`boundary number(s) [...] do not exist`)
/// or distorts self/mutual capacitances by `+20%..+60%` vs `z_coarsen: false`.
/// Expected to FAIL on unmodified `HEAD` (`9480093a`).
#[test]
fn bug3_sheet3d_z_coarsen_cap_collapse_and_distortion() {
    let run_c02_box = |z_coarsen: bool| -> Result<(f64, f64), String> {
        let mut e = Engine::new();
        let (half, pad, h_sub, h_air, h_bot) = (500.0e-6, 1200.0e-6, 750.0e-6, 890.0e-6, 900.0e-6);
        call(&mut e, "project_new", json!({ "title": "Bug3 z_coarsen", "dim": 3, "h": 340.0e-6 }));
        call(
            &mut e,
            "add_box",
            json!({
                "name": "vac_box",
                "min": [-half - pad, -half - pad, -h_sub - h_bot - pad],
                "max": [half + pad, half + pad, h_air + pad]
            }),
        );
        call(&mut e, "add_box", json!({ "name": "sub_box", "min": [-half, -half, -h_sub], "max": [half, half, 0.0] }));
        call(
            &mut e,
            "add_polygons",
            json!({
                "name": "cross_poly",
                "polygons": [[[-15e-6, -120e-6], [15e-6, -120e-6], [15e-6, -15e-6], [120e-6, -15e-6], [120e-6, 15e-6], [15e-6, 15e-6], [15e-6, 120e-6], [-15e-6, 120e-6], [-15e-6, 15e-6], [-120e-6, 15e-6], [-120e-6, -15e-6], [-15e-6, -15e-6]]],
                "simplify_tol": 50e-9
            }),
        );
        call(
            &mut e,
            "add_polygons",
            json!({
                "name": "claw_poly",
                "polygons": [[[-45e-6, 132e-6], [45e-6, 132e-6], [45e-6, 190e-6], [30e-6, 190e-6], [30e-6, 147e-6], [-30e-6, 147e-6], [-30e-6, 190e-6], [-45e-6, 190e-6]]],
                "simplify_tol": 50e-9
            }),
        );
        call(
            &mut e,
            "add_polygons",
            json!({
                "name": "gnd_poly",
                "polygons": [[[-half, -half], [half, -half], [half, half], [-half, half]]],
                "holes": [[[-140e-6, -140e-6], [140e-6, -140e-6], [140e-6, 210e-6], [-140e-6, 210e-6]]],
                "simplify_tol": 50e-9
            }),
        );
        call(&mut e, "add_sheet", json!({ "name": "cross", "input": "cross_poly", "z": 0.0 }));
        call(&mut e, "add_sheet", json!({ "name": "claw", "input": "claw_poly", "z": 0.0 }));
        call(&mut e, "add_sheet", json!({ "name": "ground_plane", "input": "gnd_poly", "z": 0.0 }));
        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "vacuum", "input": "vac_box", "kind": "Material" },
                    { "name": "substrate", "input": "sub_box", "kind": "Material" }
                ]
            }),
        );
        call(&mut e, "add_size", json!({ "sel": "cross:edge", "h": 10.0e-6 }));
        call(&mut e, "add_size", json!({ "sel": "claw:edge", "h": 8.0e-6 }));
        call(
            &mut e,
            "set_mesh",
            json!({
                "h": 340.0e-6,
                "z_coarsen": z_coarsen,
                "dz_min": 5.0e-6,
                "z_growth_ratio": 1.6,
                "max_z_slabs": 22,
                "max_cells": 350_000,
                "max_steiner_multiplier": 8,
                "strict_mesh": false,
                "structured_pml": false
            }),
        );
        call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
        call(
            &mut e,
            "set_physics_params",
            json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1.0e-9 } }),
        );
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }),
        );
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum", "material": { "eps_r": 1.0 } }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "cross", "kind": "terminal", "params": { "name": "cross" } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "claw", "kind": "terminal", "params": { "name": "claw" } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "ground_plane", "kind": "terminal", "params": { "name": "ground" } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2, "center": [0.0, 0.0, 0.0] } }),
        );

        e.call("solve_run", &json!({ "physics": "dfs", "isolate": false, "save_fields": false }))
            .map_err(|err| err.to_string())?;
        let cap = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
        Ok((cap["maxwell_matrix_F"][0][0].as_f64().unwrap(), -cap["maxwell_matrix_F"][0][1].as_f64().unwrap()))
    };

    let (c_xx_ref, c_xc_ref) = run_c02_box(false).expect("baseline z_coarsen=false solve must succeed");
    let (c_xx_zc, c_xc_zc) = run_c02_box(true)
        .expect("Core Bug 3 (sheet3d.rs z_coarsen=true): cap boundary collapse during multi-layer solve");

    let rel_xx = (c_xx_zc - c_xx_ref).abs() / c_xx_ref;
    let rel_xc = (c_xc_zc - c_xc_ref).abs() / c_xc_ref;
    assert!(
        rel_xx < 0.02 && rel_xc < 0.02,
        "Core Bug 3 (sheet3d.rs z_coarsen=true): capacitance distorted vs z_coarsen=false: C_xx {:.2}% ({c_xx_zc:e} vs {c_xx_ref:e}), C_xc {:.2}% ({c_xc_zc:e} vs {c_xc_ref:e})",
        rel_xx * 100.0,
        rel_xc * 100.0
    );
}

/// Core Bug 4 Regression (`crates/femcore/src/meshing/mesher.rs` / `sheet3d.rs` narrow-gap PSLG Steiner locking):
/// On `S1.4` (`GeneralizedCapNInterdigital`) and `Ref.0` (`QCap` GDS) geometries with dense `~200-300 nm`
/// rounded-tip / corner-fillet chords adjacent to `2 um` interdigital gaps inside a `500 um` ground pocket,
/// Ruppert 2D refinement exhausts `max_steiner_multiplier` near the short fillet chords and leaves unrefined
/// single-prism cross-gap triangles spanning the `2 um` interdigital gap (`+100%..+400%` capacitance locking),
/// whereas `planar_bem` (`QFR-4`) solves the identical polygons without volumetric gap locking.
/// Expected to FAIL on unmodified `HEAD` (`9480093a`).
#[test]
fn bug4_sheet3d_narrow_gap_steiner_locking() {
    type PolyRings = Vec<Vec<[f64; 2]>>;
    let build_rounded_interdigital_polys = || -> (PolyRings, PolyRings, Vec<[f64; 2]>) {
        // Interdigital fingers along x in [-60, 60] um with 2.0 um gap, 4.0 um width, and 20-chord semicircular
        // tips (~314 nm chords) that exhaust Ruppert's default Steiner budget before the 2 um parallel gaps refine.
        let mut bus_a = vec![[-75.0e-6, -22.0e-6], [-60.0e-6, -22.0e-6]];
        let mut bus_b = vec![[75.0e-6, 22.0e-6], [60.0e-6, 22.0e-6]];
        let w = 4.0e-6;
        let g = 2.0e-6;
        let pitch = 2.0 * (w + g); // 12 um per pair
        let r = w / 2.0; // 2.0 um tip radius -> pi*r / 20 = 314 nm chord!
        let n_arc = 20usize;

        // 3 fingers extending from left bus (-60 um -> +48 um)
        for k in 0..3 {
            let yc = -15.0e-6 + (k as f64) * pitch;
            bus_a.push([-60.0e-6, yc - r]);
            for i in 0..=n_arc {
                let th = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc as f64);
                bus_a.push([48.0e-6 + r * th.cos(), yc + r * th.sin()]);
            }
            bus_a.push([-60.0e-6, yc + r]);
        }
        bus_a.push([-60.0e-6, 22.0e-6]);
        bus_a.push([-75.0e-6, 22.0e-6]);

        // 3 interleaved fingers extending from right bus (+60 um -> -48 um), separated by 2.0 um gaps
        for k in (0..3).rev() {
            let yc = -15.0e-6 + (w + g) + (k as f64) * pitch;
            bus_b.push([60.0e-6, yc + r]);
            for i in 0..=n_arc {
                let th = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc as f64);
                bus_b.push([-48.0e-6 + r * th.cos(), yc + r * th.sin()]);
            }
            bus_b.push([60.0e-6, yc - r]);
        }
        bus_b.push([60.0e-6, -22.0e-6]);
        bus_b.push([75.0e-6, -22.0e-6]);

        // Filleted pocket hole [-85, 85] x [-30, 30] um with 20-chord R=4 um corners (~314 nm chords)
        let mut hole = Vec::new();
        let (hx, hy, hr) = (85.0e-6, 30.0e-6, 4.0e-6);
        let corners = [
            (hx - hr, -hy + hr, -std::f64::consts::FRAC_PI_2),
            (hx - hr, hy - hr, 0.0),
            (-hx + hr, hy - hr, std::f64::consts::FRAC_PI_2),
            (-hx + hr, -hy + hr, std::f64::consts::PI),
        ];
        for &(cx, cy, th0) in &corners {
            for i in 0..=n_arc {
                let th = th0 + std::f64::consts::FRAC_PI_2 * (i as f64) / (n_arc as f64);
                hole.push([cx + hr * th.cos(), cy + hr * th.sin()]);
            }
        }
        (vec![bus_a], vec![bus_b], hole)
    };

    let solve_idc = |planar_bem: bool| -> (f64, f64) {
        let (pa, pb, pocket_hole) = build_rounded_interdigital_polys();
        let mut e = Engine::new();
        let half = 260.0e-6;
        call(&mut e, "project_new", json!({ "title": "Bug4 S1.4 Narrow Gap", "dim": 3, "h": 95.0e-6 }));
        call(
            &mut e,
            "add_box",
            json!({ "name": "vac_bot", "min": [-half, -half, -500.0e-6], "max": [half, half, -250.0e-6] }),
        );
        call(
            &mut e,
            "add_box",
            json!({ "name": "sub_box", "min": [-half, -half, -250.0e-6], "max": [half, half, 0.0] }),
        );
        call(
            &mut e,
            "add_box",
            json!({ "name": "vac_top", "min": [-half, -half, 0.0], "max": [half, half, 250.0e-6] }),
        );
        call(&mut e, "add_polygons", json!({ "name": "p_a", "polygons": pa, "simplify_tol": 1.0e-9 }));
        call(&mut e, "add_polygons", json!({ "name": "p_b", "polygons": pb, "simplify_tol": 1.0e-9 }));
        call(
            &mut e,
            "add_polygons",
            json!({
                "name": "p_gnd",
                "polygons": [[[-half, -half], [half, -half], [half, half], [-half, half]]],
                "holes": [pocket_hole],
                "simplify_tol": 1.0e-9
            }),
        );
        call(&mut e, "add_sheet", json!({ "name": "cap_a", "input": "p_a", "z": 0.0 }));
        call(&mut e, "add_sheet", json!({ "name": "cap_b", "input": "p_b", "z": 0.0 }));
        call(&mut e, "add_sheet", json!({ "name": "gnd", "input": "p_gnd", "z": 0.0 }));
        call(
            &mut e,
            "add_overlay",
            json!({
                "regions": [
                    { "name": "vacuum_bot", "input": "vac_bot", "kind": "Material" },
                    { "name": "substrate", "input": "sub_box", "kind": "Material" },
                    { "name": "vacuum_top", "input": "vac_top", "kind": "Material" }
                ]
            }),
        );

        if planar_bem {
            call(&mut e, "add_size", json!({ "sel": "cap_a:edge", "h": 4.0e-6 }));
            call(&mut e, "add_size", json!({ "sel": "cap_b:edge", "h": 4.0e-6 }));
            call(
                &mut e,
                "set_mesh",
                json!({
                    "h": 95.0e-6,
                    "z_coarsen": false,
                    "dz_min": 50.0e-6,
                    "max_z_slabs": 4,
                    "max_steiner_multiplier": 4,
                    "strict_mesh": false,
                    "structured_pml": false
                }),
            );
            call(
                &mut e,
                "add_physics",
                json!({
                    "physics": "dfs",
                    "study": "electrostatic",
                    "order": 1,
                    "study_params": { "formulation": "planar_bem", "richardson": false }
                }),
            );
        } else {
            call(
                &mut e,
                "set_mesh",
                json!({
                    "h": 95.0e-6,
                    "z_coarsen": false,
                    "dz_min": 2.0e-6,
                    "z_growth_ratio": 2.4,
                    "max_z_slabs": 12,
                    "max_steiner_multiplier": 8,
                    "singular_edges": "quarter_point",
                    "quarter_point_xi": 0.25,
                    "corner_grading": true,
                    "strict_mesh": false,
                    "structured_pml": false
                }),
            );
            call(
                &mut e,
                "add_physics",
                json!({
                    "physics": "dfs",
                    "study": "electrostatic",
                    "order": 2,
                    "study_params": { "formulation": "volumetric_fem", "singular_edges": "quarter_point" }
                }),
            );
        }

        call(&mut e, "set_physics_params", json!({ "physics": "dfs", "params": { "solver": "direct" } }));
        call(
            &mut e,
            "set_material",
            json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }),
        );
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum_top", "material": { "eps_r": 1.0 } }));
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }));
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "cap_a", "kind": "terminal", "params": { "name": "a" } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "cap_b", "kind": "terminal", "params": { "name": "b" } }),
        );
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "gnd", "kind": "terminal", "params": { "name": "gnd" } }),
        );
        if !planar_bem {
            call(
                &mut e,
                "add_condition",
                json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2 } }),
            );
        }
        call(&mut e, "solve_run", json!({ "isolate": false, "save_fields": false }));
        let cap = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
        let c_aa = cap["maxwell_matrix_F"][0][0].as_f64().unwrap();
        let c_ab = -cap["maxwell_matrix_F"][0][1].as_f64().unwrap();
        (c_aa, c_ab)
    };

    let (c_aa_bem, c_ab_bem) = solve_idc(true);
    let (c_aa_fem, c_ab_fem) = solve_idc(false);
    let rel_aa = (c_aa_fem - c_aa_bem).abs() / c_aa_bem;
    let rel_ab = (c_ab_fem - c_ab_bem).abs() / c_ab_bem;
    assert!(
        rel_aa <= 0.15 && rel_ab <= 0.15,
        "Core Bug 4 (mesher.rs / sheet3d.rs narrow-gap PSLG Steiner locking): volumetric FEM locked vs planar_bem on 2 um filleted IDC: C_aa {:.1}% ({c_aa_fem:e} vs {c_aa_bem:e}), C_ab {:.1}% ({c_ab_fem:e} vs {c_ab_bem:e})",
        rel_aa * 100.0,
        rel_ab * 100.0
    );
}

/// Core Bug 5 Regression (`crates/femcore` / `nyx-ems` parallel `faer` sparse factorization non-determinism):
/// ADR 11 requires identical inputs to produce identical outputs (`< 1e-13` relative drift) regardless of
/// `NYX_THREADS`. While `dfs` AMG is bitwise deterministic (`0.000e0`), 3D `ems::cavity` with parallel `faer`
/// sparse direct factorization exhibits `~1e-10 .. 1e-7` thread-count drift between `NYX_THREADS=1` and `8`.
/// Expected to FAIL on unmodified `HEAD` (`9480093a`).
#[test]
fn bug5_ems_thread_determinism_drift() {
    if std::env::var("NYX_SMOKE_SUBPROC_EMS_DET").as_deref() == Ok("1") {
        let mut e_ems = Engine::new();
        build_clt_chip(&mut e_ems, true);
        call(
            &mut e_ems,
            "add_physics",
            json!({ "physics": "ems", "study": "cavity", "order": 1, "study_params": { "n_modes": 3, "f_search": 6.5e9 } }),
        );
        call(
            &mut e_ems,
            "set_material",
            json!({ "physics": "ems", "region": "vac_bot", "material": { "eps_r": 1.0 } }),
        );
        call(&mut e_ems, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
        call(&mut e_ems, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
        for s in ["ground_plane", "feed_trace", "clt_res", "pad_a", "pad_b"] {
            call(&mut e_ems, "add_condition", json!({ "physics": "ems", "on": s, "kind": "pec", "params": {} }));
        }
        call(
            &mut e_ems,
            "add_condition",
            json!({
                "physics": "ems",
                "on": "jj_sheet",
                "kind": "lumped_jj",
                "params": { "name": "Q1_JJ", "E_j_Hz": 14.0e9, "C_j": 40.0e-15, "direction": [0.0, 1.0, 0.0], "formulation": "rank1" }
            }),
        );
        call(
            &mut e_ems,
            "add_condition",
            json!({
                "physics": "ems",
                "on": "port_in",
                "kind": "lumped_port",
                "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
            }),
        );
        call(&mut e_ems, "solve_run", json!({ "isolate": false }));
        let ham = call(
            &mut e_ems,
            "analyse",
            json!({ "analysis": "hamiltonian", "physics": "ems", "params": { "basis": "craig_bampton" } }),
        );
        let net = call(&mut e_ems, "analyse", json!({ "analysis": "network", "physics": "ems" }));
        let f01 = ham["f01_Hz"][0].as_f64().unwrap();
        let alpha = ham["alpha_Hz"][0].as_f64().unwrap();
        let f_qep = net["qep_modes"][0]["f_re_Hz"].as_f64().unwrap();
        let kappa = net["qep_modes"][0]["kappa_Hz"].as_f64().unwrap();
        println!("NYX_SQUADDS_EMS_DET:{}", json!([f01, alpha, f_qep, kappa]));
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |threads: &str| -> Vec<f64> {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("bug5_ems_thread_determinism_drift")
            .arg("--nocapture")
            .env("NYX_THREADS", threads)
            .env("NYX_SMOKE_SUBPROC_EMS_DET", "1")
            .output()
            .expect("spawn ems determinism subprocess");
        assert!(out.status.success(), "NYX_THREADS={threads} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let payload = stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_SQUADDS_EMS_DET:"))
            .unwrap_or_else(|| panic!("missing NYX_SQUADDS_EMS_DET line in stdout:\n{stdout}"));
        serde_json::from_str(payload).unwrap()
    };

    let y1 = run_with_threads("1");
    let y8 = run_with_threads("8");
    let names = ["f01_Hz", "alpha_Hz", "f_qep_Hz", "kappa_Hz"];
    for ((&v1, &v8), name) in y1.iter().zip(&y8).zip(&names) {
        let rel = (v8 - v1).abs() / v1.abs().max(1.0e-30);
        assert!(
            rel < 1.0e-13,
            "Core Bug 5 (ADR 11 EMS parallel faer factorization drift) on {name}: NYX_THREADS=1 ({v1:.16e}) vs NYX_THREADS=8 ({v8:.16e}), rel={rel:e} >= 1e-13"
        );
    }
}
