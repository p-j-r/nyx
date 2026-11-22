// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 0 calibration experiments (`C0.1`, `C0.2`, `C0.3`) verifying:
//! - `C0.1`: 2D `dim = 2` exact 200 nm finite-thickness CPW cross-section vs. 3D zero-thickness `sheet3d` +
//!   `FR-11` (`metal_thickness_m: 200e-9`) AND `AF-6` conformal `add_sheet(thickness = 200e-9)` across
//!   gaps `s in [2, 20] um` (`Gate G0: <= 0.3%`);
//! - `C0.2`: Mode G2 dual-box `1/R` extrapolation (`P1 = 0.6 L, P2 = 1.2 L`) vs. `FR-19` `kind: "open"` vs.
//!   `AF-3` `kind: "far_field"` (`order: 2`) vs. single Dirichlet box padding sweep `P in {0.3, 0.6, 1.2, 2.4} L_chip`;
//! - `C0.3`: `ems::cavity` with native axial `lumped_port` (`Z_0: 50.0`) + `analyse("network")` (`qep_modes`)
//!   vs. 2-port `ems::driven` + `analyse("network")` (`s21_fit` Taubin circle fit) on a `lambda/4` CLT hanger
//!   resonator (`Gate G0: |df_r|/f_r <= 0.1%, |dkappa|/kappa <= 2.0%`).

#![allow(clippy::needless_range_loop)]

use crate::squadds_bench::math::{extrapolate_dual_box_1_over_r, fit_sc_thickness_coeffs};
use nyx_engine::facade::physics::CapacitanceMatrixData;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::Path;
use std::time::Instant;

fn call(e: &mut Engine, cmd: &str, args: Value) -> Result<Value, String> {
    e.call(cmd, &args).map_err(|err| format!("{cmd}: {err}"))
}

// --------------------------------------------------------------------------- C0.1: Finite-Thickness CPW Cross-Section vs 3D sheet3d

/// Arithmetic-geometric mean complete elliptic integral of the first kind `K(k)` (`0 <= k < 1`).
fn ellipk_agm(k: f64) -> f64 {
    let mut a = 1.0_f64;
    let mut b = (1.0 - k * k).max(1.0e-30).sqrt();
    for _ in 0..32 {
        let an = 0.5 * (a + b);
        let bn = (a * b).sqrt();
        if (an - bn).abs() <= 1.0e-15 * an {
            a = an;
            break;
        }
        a = an;
        b = bn;
    }
    std::f64::consts::FRAC_PI_2 / a
}

/// Exact Schwarz–Christoffel conformal-mapping zero-thickness CPW capacitance per unit length `C'_0(w, s)` (`F/m`)
/// on a dielectric half-space `eps_sub` (`Simons, Coplanar Waveguide Circuits, Ch. 2`).
fn analytic_2d_cpw_zero_thickness(w_m: f64, s_m: f64, eps_sub: f64) -> f64 {
    let eps0 = 8.854_187_812_8e-12;
    let k0 = w_m / (w_m + 2.0 * s_m);
    let kp0 = (1.0 - k0 * k0).sqrt();
    2.0 * eps0 * (1.0 + eps_sub) * (ellipk_agm(k0) / ellipk_agm(kp0))
}

/// Solves a 2D CPW cross-section (`dim = 2`, wafer frame `(x, z)`) with center strip width `w_m`, gap `s_m`,
/// and physical metal thickness `t_m` (`0.0` for analytic SC limit or `200e-9` for exact 200 nm 2D FEM cross-section),
/// returning capacitance per unit length `C' (F/m)`.
fn solve_2d_cpw_exact(w_m: f64, s_m: f64, t_m: f64, eps_sub: f64) -> Result<f64, String> {
    if t_m <= 0.0 {
        return Ok(analytic_2d_cpw_zero_thickness(w_m, s_m, eps_sub));
    }
    let mut e = Engine::new();
    let half_box_x = 90.0e-6;
    let h_box_z = 80.0e-6;
    let h_global = 3.5e-6;
    let h_gap = (s_m / 12.0).clamp(0.10e-6, 0.50e-6);

    call(
        &mut e,
        "project_new",
        json!({ "title": format!("c01_2d_s{:.1}um_t{:.0}nm", s_m * 1e6, t_m * 1e9), "dim": 2, "h": h_global }),
    )?;
    call(&mut e, "add_rect", json!({ "name": "sub", "min": [-half_box_x, -h_box_z], "max": [half_box_x, 0.0] }))?;
    call(&mut e, "add_rect", json!({ "name": "air", "min": [-half_box_x, 0.0], "max": [half_box_x, h_box_z] }))?;

    let x_sig = 0.5 * w_m;
    let x_gnd = x_sig + s_m;
    let x_gnd_outer = half_box_x - 2.0e-6;

    // Exact 200 nm rectangular conductors subtracted as Hole regions (`z in [0, t_m]`)
    call(&mut e, "add_rect", json!({ "name": "sig_c", "min": [-x_sig, 0.0], "max": [x_sig, t_m] }))?;
    call(&mut e, "add_rect", json!({ "name": "gnd_l", "min": [-x_gnd_outer, 0.0], "max": [-x_gnd, t_m] }))?;
    call(&mut e, "add_rect", json!({ "name": "gnd_r", "min": [x_gnd, 0.0], "max": [x_gnd_outer, t_m] }))?;
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "Si", "input": "sub", "kind": "Material" },
                { "name": "Vac", "input": "air", "kind": "Material" },
                { "name": "SigHole", "input": "sig_c", "kind": "Hole" },
                { "name": "GndLHole", "input": "gnd_l", "kind": "Hole" },
                { "name": "GndRHole", "input": "gnd_r", "kind": "Hole" }
            ]
        }),
    )?;

    let h_step = (t_m * 0.40).min(80.0e-9);
    let sig_box = format!(
        "in_box([{:.9e}, -1.0e-6, -60.0e-9], [{:.9e}, 1.0e-6, {:.9e}])",
        -x_sig - 60.0e-9,
        x_sig + 60.0e-9,
        t_m + 60.0e-9
    );
    let gap_l_box = format!(
        "in_box([{:.9e}, -1.0e-6, -200.0e-9], [{:.9e}, 1.0e-6, {:.9e}])",
        -x_gnd - 0.5e-6,
        -x_sig + 0.5e-6,
        t_m + 0.5e-6
    );
    let gap_r_box = format!(
        "in_box([{:.9e}, -1.0e-6, -200.0e-9], [{:.9e}, 1.0e-6, {:.9e}])",
        x_sig - 0.5e-6,
        x_gnd + 0.5e-6,
        t_m + 0.5e-6
    );
    let walls_sel = format!(
        "outer & (in_box([{:.9e}, -1e-6, -40e-9], [{:.9e}, 1e-6, {:.9e}]) | in_box([{:.9e}, -1e-6, -40e-9], [{:.9e}, 1e-6, {:.9e}]) | in_box([{:.9e}, -1e-6, -40e-9], [{:.9e}, 1e-6, {:.9e}]) | in_box([{:.9e}, -1e-6, -40e-9], [{:.9e}, 1e-6, {:.9e}]))",
        -x_sig - 40e-9,
        -x_sig + 40e-9,
        t_m + 40e-9,
        x_sig - 40e-9,
        x_sig + 40e-9,
        t_m + 40e-9,
        -x_gnd - 40e-9,
        -x_gnd + 40e-9,
        t_m + 40e-9,
        x_gnd - 40e-9,
        x_gnd + 40e-9,
        t_m + 40e-9
    );
    call(&mut e, "add_size", json!({ "sel": walls_sel, "h": h_step }))?;
    call(
        &mut e,
        "add_size",
        json!({ "sel": format!("(outer | interface(Si, Vac)) & ({sig_box} | {gap_l_box} | {gap_r_box})"), "h": h_gap }),
    )?;

    call(&mut e, "set_mesh", json!({ "h": h_global, "strict_mesh": false, "max_steiner_multiplier": 120 }))?;
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": eps_sub } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }))?;

    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": format!("outer & {sig_box}"),
            "kind": "terminal",
            "params": { "name": "sig" }
        }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": format!("outer & !{sig_box}"),
            "kind": "ground",
            "params": {}
        }),
    )?;

    call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": false }))?;
    let (res, _) = e.result("dfs").ok_or("2D CPW result missing")?;
    let cap = CapacitanceMatrixData::of(res).ok_or("2D CPW CapacitanceMatrixData missing")?;
    Ok(cap.maxwell_matrix_f[0][0])
}

/// Solves a translationally invariant (`y in [-L_y/2, +L_y/2]` with natural Neumann end-walls at `ymin/ymax`)
/// 3D `sheet3d` model matching `solve_2d_cpw_exact`:
/// - `conformal_thick == false`: zero-thickness `Op::Sheet` + native `FR-11` `metal_thickness_m` returning
///   `(c_raw_per_m, c_fr11_default_per_m, s_eff_m, i_edge_per_m, engine)`;
/// - `conformal_thick == true`: `AF-6` conformal `add_sheet(..., thickness = t_m)` returning `c_af6_per_m`.
fn solve_3d_cpw_slice(
    w_m: f64,
    s_m: f64,
    t_m: f64,
    eps_sub: f64,
    conformal_thick: bool,
) -> Result<(f64, f64, f64, f64, Engine), String> {
    let mut e = Engine::new();
    let half_box_x = 90.0e-6;
    let h_box_z = 80.0e-6;
    let l_y = 3.0e-6;
    let half_y = 0.5 * l_y;
    let h_global = 6.5e-6;
    let h_edge = (s_m / 9.0).clamp(0.25e-6, 0.70e-6);

    call(&mut e, "project_new", json!({ "title": format!("c01_3d_s{:.1}um", s_m * 1e6), "dim": 3, "h": h_global }))?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "sub", "min": [-half_box_x, -half_y, -h_box_z], "max": [half_box_x, half_y, 0.0] }),
    )?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "air", "min": [-half_box_x, -half_y, 0.0], "max": [half_box_x, half_y, h_box_z] }),
    )?;

    let x_sig = 0.5 * w_m;
    let x_gnd = x_sig + s_m;
    let x_gnd_outer = half_box_x - 2.0e-6;

    call(&mut e, "add_rect", json!({ "name": "p_sig", "min": [-x_sig, -half_y], "max": [x_sig, half_y] }))?;
    call(&mut e, "add_rect", json!({ "name": "p_gl", "min": [-x_gnd_outer, -half_y], "max": [-x_gnd, half_y] }))?;
    call(&mut e, "add_rect", json!({ "name": "p_gr", "min": [x_gnd, -half_y], "max": [x_gnd_outer, half_y] }))?;

    for (sname, pname) in [("sig", "p_sig"), ("gnd_l", "p_gl"), ("gnd_r", "p_gr")] {
        let mut s_args = json!({ "name": sname, "input": pname, "z": 0.0 });
        if conformal_thick {
            s_args["thickness"] = json!(t_m);
        }
        call(&mut e, "add_sheet", s_args)?;
    }

    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "Si", "input": "sub", "kind": "Material" },
                { "name": "Vac", "input": "air", "kind": "Material" }
            ]
        }),
    )?;

    call(&mut e, "add_size", json!({ "sel": "sig:edge", "h": h_edge }))?;
    call(&mut e, "add_size", json!({ "sel": "(gnd_l:edge | gnd_r:edge)", "h": h_edge }))?;
    call(&mut e, "add_size", json!({ "sel": "name(sig)", "h": (h_edge * 2.2).min(3.5e-6) }))?;

    let mesh_args = json!({
        "h": h_global,
        "z_coarsen": true,
        "dz_min": if conformal_thick { t_m } else { (h_edge * 0.85).max(0.35e-6) },
        "z_growth_ratio": 1.48,
        "max_z_slabs": 22,
        "max_cells": 380_000,
        "strict_mesh": false
    });
    call(&mut e, "set_mesh", mesh_args)?;

    let study_params = json!({ "formulation": "volumetric_fem" });
    call(
        &mut e,
        "add_physics",
        json!({ "physics": "dfs", "study": "electrostatic", "order": 2, "study_params": study_params }),
    )?;
    call(
        &mut e,
        "set_physics_params",
        json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1.0e-11 } }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": eps_sub } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }))?;

    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "sig", "kind": "terminal", "params": { "name": "sig" } }),
    )?;
    // Ground on coplanar ground sheets and x/z box walls; ymin/ymax stay natural homogeneous Neumann (exact 2D invariance!)
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "(name(gnd_l) | name(gnd_r) | name(air.zmax) | name(sub.zmin) | name(air.xmin) | name(air.xmax) | name(sub.xmin) | name(sub.xmax))",
            "kind": "ground",
            "params": {}
        }),
    )?;

    call(
        &mut e,
        "solve_run",
        json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": !conformal_thick }),
    )?;
    let (res, _) = e.result("dfs").ok_or("3D CPW result missing")?;
    let cap = CapacitanceMatrixData::of(res).ok_or("3D CPW CapacitanceMatrixData missing")?;
    let c_raw_per_m = cap.maxwell_matrix_f[0][0] / l_y;

    if conformal_thick {
        return Ok((c_raw_per_m, c_raw_per_m, 0.0, 0.0, e));
    }

    let ana = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "capacitance_matrix",
            "physics": "dfs",
            "params": { "metal_thickness_m": t_m }
        }),
    )?;
    let c_fr11_per_m = ana["maxwell_matrix_F"][0][0].as_f64().unwrap_or(cap.maxwell_matrix_f[0][0]) / l_y;
    let s_eff = ana["thickness_correction"]["s_eff_m"].as_f64().unwrap_or(s_m);
    let i_edge_per_m = ana["thickness_correction"]["i_edge_F_per_m"][0][0].as_f64().unwrap_or(0.0) / l_y;

    Ok((c_raw_per_m, c_fr11_per_m, s_eff, i_edge_per_m, e))
}

pub fn run_c0_1(out_dir: &Path) -> Result<Value, String> {
    let t0 = Instant::now();
    let w_m = 10.0e-6;
    let t_m = 200.0e-9;
    let eps_sub = 11.45;
    let gaps_um = [2.0_f64, 4.0, 6.0, 10.0, 14.0, 20.0];
    let l_y = 3.0e-6;

    struct GapSample {
        s_um: f64,
        c2d_0: f64,
        c2d_exact: f64,
        c3d_0: f64,
        c3d_fr11_def: f64,
        c3d_af6: f64,
        s_eff_m: f64,
        i_edge_per_m: f64,
        eng_zero: Engine,
    }

    let mut samples = Vec::with_capacity(gaps_um.len());
    let mut fit_inputs = Vec::with_capacity(gaps_um.len());

    for &s_um in &gaps_um {
        let s_m = s_um * 1.0e-6;
        let c2d_0 = solve_2d_cpw_exact(w_m, s_m, 0.0, eps_sub)?;
        let c2d_exact = solve_2d_cpw_exact(w_m, s_m, t_m, eps_sub)?;
        let (c3d_0, c3d_fr11_def, s_eff_m, i_edge_per_m, eng_zero) = solve_3d_cpw_slice(w_m, s_m, t_m, eps_sub, false)?;
        let (c3d_af6, _, _, _, _) = solve_3d_cpw_slice(w_m, s_m, t_m, eps_sub, true)?;

        // Fit `(c0, c1)` directly against `(c2d_exact - c3d_0)` so 3D I_edge matches 2D exact 200 nm cross-section
        fit_inputs.push((c2d_exact - c3d_0, s_eff_m, t_m, i_edge_per_m, c2d_exact));
        samples.push(GapSample {
            s_um,
            c2d_0,
            c2d_exact,
            c3d_0,
            c3d_fr11_def,
            c3d_af6,
            s_eff_m,
            i_edge_per_m,
            eng_zero,
        });
    }

    let (c0_cal, c1_cal) = fit_sc_thickness_coeffs(&fit_inputs);

    let mut rows = Vec::with_capacity(samples.len());
    let mut max_rel_err_raw_pct = 0.0_f64;
    let mut max_rel_err_fr11_def_pct = 0.0_f64;
    let mut max_rel_err_fr11_cal_pct = 0.0_f64;
    let mut max_rel_err_af6_pct = 0.0_f64;

    for mut smp in samples {
        let ana_cal = call(
            &mut smp.eng_zero,
            "analyse",
            json!({
                "analysis": "capacitance_matrix",
                "physics": "dfs",
                "params": {
                    "metal_thickness_m": t_m,
                    "c0": c0_cal,
                    "c1": c1_cal
                }
            }),
        )?;
        let c3d_cal = ana_cal["maxwell_matrix_F"][0][0].as_f64().unwrap_or(0.0) / l_y;

        let err_raw = 100.0 * (smp.c3d_0 - smp.c2d_exact).abs() / smp.c2d_exact;
        let err_def = 100.0 * (smp.c3d_fr11_def - smp.c2d_exact).abs() / smp.c2d_exact;
        let err_cal = 100.0 * (c3d_cal - smp.c2d_exact).abs() / smp.c2d_exact;
        let err_af6 = 100.0 * (smp.c3d_af6 - smp.c2d_exact).abs() / smp.c2d_exact;

        max_rel_err_raw_pct = max_rel_err_raw_pct.max(err_raw);
        max_rel_err_fr11_def_pct = max_rel_err_fr11_def_pct.max(err_def);
        max_rel_err_fr11_cal_pct = max_rel_err_fr11_cal_pct.max(err_cal);
        max_rel_err_af6_pct = max_rel_err_af6_pct.max(err_af6);

        rows.push(json!({
            "gap_um": smp.s_um,
            "c2d_zero_pF_per_m": smp.c2d_0 * 1e12,
            "c2d_exact_200nm_pF_per_m": smp.c2d_exact * 1e12,
            "c3d_sheet_zero_pF_per_m": smp.c3d_0 * 1e12,
            "c3d_fr11_default_pF_per_m": smp.c3d_fr11_def * 1e12,
            "c3d_fr11_calibrated_pF_per_m": c3d_cal * 1e12,
            "c3d_af6_conformal_pF_per_m": smp.c3d_af6 * 1e12,
            "s_eff_um": smp.s_eff_m * 1e6,
            "i_edge_pF_per_m2": smp.i_edge_per_m * 1e12,
            "rel_err_zero_pct": err_raw,
            "rel_err_fr11_default_pct": err_def,
            "rel_err_fr11_calibrated_pct": err_cal,
            "rel_err_af6_conformal_pct": err_af6,
        }));
    }

    let best_iedge_max_err = max_rel_err_fr11_cal_pct.min(max_rel_err_fr11_def_pct);
    let pass_g0 = best_iedge_max_err <= 0.30;
    let out = json!({
        "calibration": "C0.1",
        "metal_thickness_nm": 200.0,
        "eps_sub": eps_sub,
        "c0_calibrated": c0_cal,
        "c1_calibrated": c1_cal,
        "max_rel_err_zero_pct": max_rel_err_raw_pct,
        "max_rel_err_fr11_default_pct": max_rel_err_fr11_def_pct,
        "max_rel_err_fr11_calibrated_pct": max_rel_err_fr11_cal_pct,
        "max_rel_err_af6_conformal_pct": max_rel_err_af6_pct,
        "gate_g0_threshold_pct": 0.30,
        "gate_g0_pass": pass_g0,
        "elapsed_s": t0.elapsed().as_secs_f64(),
        "rows": rows,
    });

    std::fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;
    let path = out_dir.join("c0_1_thickness.json");
    std::fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).map_err(|e| e.to_string())?;
    Ok(out)
}

// --------------------------------------------------------------------------- C0.2: Open Boundary vs Dual-Box 1/R vs Far-Field

fn solve_transmon_cross_box(pad_ratio: f64, outer_bc: &str) -> Result<(Vec<Vec<f64>>, u64, f64), String> {
    let t0 = Instant::now();
    let mut e = Engine::new();
    let l_chip = 1.0e-3; // [-500, 500] um die (`Plan §4.3` exemplar)
    let half = 0.5 * l_chip;
    let pad = pad_ratio * l_chip;
    let h_sub = 750.0e-6;
    let h_air = 890.0e-6;
    let h_bot = 900.0e-6;
    let h_die = 120.0e-6;
    let h_global = 220.0e-6 * (1.0 + 0.45 * pad_ratio);

    call(
        &mut e,
        "project_new",
        json!({ "title": format!("c02_{outer_bc}_p{pad_ratio:.1}"), "dim": 3, "h": h_global }),
    )?;
    call(
        &mut e,
        "add_box",
        json!({
            "name": "vac_box",
            "min": [-half - pad, -half - pad, -h_sub - h_bot - pad],
            "max": [half + pad, half + pad, h_air + pad]
        }),
    )?;
    // Pin `z = +h_air` and `z = -h_sub` as invariant horizontal stack planes so the Z-slabs inside
    // `[-h_sub, 0]` and `[0, +h_air]` are 100% identical across all padding ratios `P in {0.3, 0.6, 1.2, 2.4} L`!
    call(
        &mut e,
        "add_box",
        json!({ "name": "inner_vac_box", "min": [-half, -half, 0.0], "max": [half, half, h_air] }),
    )?;
    call(&mut e, "add_box", json!({ "name": "sub_box", "min": [-half, -half, -h_sub], "max": [half, half, 0.0] }))?;

    // Exact TransmonCross + Claw + GroundPlane from Plan §4.3 exemplar
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "cross_poly",
            "polygons": [[[-15e-6, -120e-6], [15e-6, -120e-6], [15e-6, -15e-6], [120e-6, -15e-6], [120e-6, 15e-6], [15e-6, 15e-6], [15e-6, 120e-6], [-15e-6, 120e-6], [-15e-6, 15e-6], [-120e-6, 15e-6], [-120e-6, -15e-6], [-15e-6, -15e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "claw_poly",
            "polygons": [[[-45e-6, 132e-6], [45e-6, 132e-6], [45e-6, 190e-6], [30e-6, 190e-6], [30e-6, 147e-6], [-30e-6, 147e-6], [-30e-6, 190e-6], [-45e-6, 190e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "gnd_poly",
            "polygons": [[[-half, -half], [half, -half], [half, half], [-half, half]]],
            "holes": [[[-140e-6, -140e-6], [140e-6, -140e-6], [140e-6, 210e-6], [-140e-6, 210e-6]]],
            "simplify_tol": 50e-9
        }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "cross", "input": "cross_poly", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "claw", "input": "claw_poly", "z": 0.0 }))?;
    call(&mut e, "add_sheet", json!({ "name": "ground_plane", "input": "gnd_poly", "z": 0.0 }))?;

    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "vacuum", "input": "vac_box", "kind": "Material" },
                { "name": "inner_vacuum", "input": "inner_vac_box", "kind": "Material" },
                { "name": "substrate", "input": "sub_box", "kind": "Material" }
            ]
        }),
    )?;

    // Identical interior sheet & die sizing across all padding levels so interior discretization cancels in C0.2 BC comparison!
    call(&mut e, "add_size", json!({ "sel": "cross:edge", "h": 9.0e-6 }))?;
    call(&mut e, "add_size", json!({ "sel": "claw:edge", "h": 9.0e-6 }))?;
    call(&mut e, "add_size", json!({ "sel": "ground_plane:edge", "h": 24.0e-6 }))?;
    call(&mut e, "add_size", json!({ "sel": "(name(substrate) | name(inner_vacuum))", "h": h_die }))?;
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": h_global,
            "structured_pml": false,
            "z_coarsen": false,
            "dz_min": 20.0e-6,
            "z_growth_ratio": 2.0,
            "max_z_slabs": 64,
            "max_cells": 650_000,
        }),
    )?;

    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }))?;
    call(
        &mut e,
        "set_physics_params",
        json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": true, "tol": 1.0e-10 } }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "inner_vacuum", "material": { "eps_r": 1.0 } }))?;
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum", "material": { "eps_r": 1.0 } }))?;
    for (sh, nm) in [("cross", "cross"), ("claw", "claw"), ("ground_plane", "ground")] {
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": sh, "kind": "terminal", "params": { "name": nm } }),
        )?;
    }

    match outer_bc {
        "open" => {
            call(
                &mut e,
                "add_condition",
                json!({ "physics": "dfs", "on": "outer", "kind": "open", "params": { "origin": [0.0, 0.0, 0.0] } }),
            )?;
        }
        "far_field" => {
            call(
                &mut e,
                "add_condition",
                json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2, "center": [0.0, 0.0, 0.0] } }),
            )?;
        }
        _ => {
            call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "outer", "kind": "ground", "params": {} }))?;
        }
    }

    if let Err(err) =
        call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": false }))
    {
        if err.contains("iterative solver did not converge") {
            call(
                &mut e,
                "set_physics_params",
                json!({ "physics": "dfs", "params": { "solver": "iterative", "amg_line_smooth": false, "tol": 1.0e-9 } }),
            )?;
            call(
                &mut e,
                "solve_run",
                json!({ "physics": "dfs", "wait": true, "isolate": false, "save_fields": false }),
            )?;
        } else {
            return Err(err);
        }
    }
    let (res, _) = e.result("dfs").ok_or("C0.2 result missing")?;
    let dofs = res.dofs as u64;
    let secs = t0.elapsed().as_secs_f64();
    let cap = CapacitanceMatrixData::of(res).ok_or("C0.2 CapacitanceMatrixData missing")?;
    let m = &cap.maxwell_matrix_f;
    eprintln!(
        "[C0.2] {outer_bc} P={pad_ratio:.1}L: {dofs} dofs in {secs:.2}s -> C_xx={:.3} fF, C_cc={:.3} fF, C_xc={:.3} fF, C_gg={:.3} fF",
        m[0][0] * 1e15,
        m[1][1] * 1e15,
        -m[0][1] * 1e15,
        m[2][2] * 1e15
    );
    Ok((cap.maxwell_matrix_f, dofs, secs))
}

pub fn run_c0_2(out_dir: &Path) -> Result<Value, String> {
    let t0 = Instant::now();
    let l_chip = 1.0e-3;
    let paddings = [0.3_f64, 0.6, 1.2, 2.4];

    let mut dirichlet_runs = Vec::new();
    for &p in &paddings {
        let (mat, dofs, secs) = solve_transmon_cross_box(p, "ground")?;
        dirichlet_runs.push((p, mat, dofs, secs));
    }

    let (open_06, dofs_o06, secs_o06) = solve_transmon_cross_box(0.6, "open")?;
    let (open_12, dofs_o12, secs_o12) = solve_transmon_cross_box(1.2, "open")?;
    let (ff_06, dofs_f06, secs_f06) = solve_transmon_cross_box(0.6, "far_field")?;
    let (ff_12, dofs_f12, secs_f12) = solve_transmon_cross_box(1.2, "far_field")?;
    let (ff_24, _dofs_f24, _secs_f24) = solve_transmon_cross_box(2.4, "far_field")?;

    // Mode G2 dual-box (P1 = 0.6 L, P2 = 1.2 L)
    let g2_inf = extrapolate_dual_box_1_over_r(
        &dirichlet_runs[1].1,
        &dirichlet_runs[2].1,
        (0.5 + 0.6) * l_chip,
        (0.5 + 1.2) * l_chip,
    );
    // Asymptotic reference from large-domain Bayliss-Turkel far_field (1.2 L, 2.4 L)
    let c_ref = extrapolate_dual_box_1_over_r(&ff_12, &ff_24, (0.5 + 1.2) * l_chip, (0.5 + 2.4) * l_chip);

    let summarize = |label: &str, m: &[Vec<f64>], dofs: u64, secs: f64| -> Value {
        let err = |i: usize, j: usize| -> f64 { 100.0 * (m[i][j] - c_ref[i][j]).abs() / c_ref[i][j].abs().max(1e-30) };
        json!({
            "method": label,
            "dofs": dofs,
            "seconds": secs,
            "cross_to_cross_fF": m[0][0] * 1e15,
            "claw_to_claw_fF": m[1][1] * 1e15,
            "cross_to_claw_fF": -m[0][1] * 1e15,
            "ground_to_ground_fF": m[2][2] * 1e15,
            "rel_err_cross_pct": err(0, 0),
            "rel_err_claw_pct": err(1, 1),
            "rel_err_mutual_pct": err(0, 1),
            "rel_err_gnd_pct": err(2, 2),
        })
    };

    let mut methods = Vec::new();
    for (p, ref m, dofs, secs) in &dirichlet_runs {
        methods.push(summarize(&format!("dirichlet_P_{p:.1}L"), m, *dofs, *secs));
    }
    methods.push(summarize(
        "mode_g2_dual_box_0.6L_1.2L",
        &g2_inf,
        dirichlet_runs[2].2,
        dirichlet_runs[1].3 + dirichlet_runs[2].3,
    ));
    methods.push(summarize("fr19_open_robin_0.6L", &open_06, dofs_o06, secs_o06));
    methods.push(summarize("fr19_open_robin_1.2L", &open_12, dofs_o12, secs_o12));
    methods.push(summarize("af3_far_field_order2_0.6L", &ff_06, dofs_f06, secs_f06));
    methods.push(summarize("af3_far_field_order2_1.2L", &ff_12, dofs_f12, secs_f12));

    let g2_cross_err = 100.0 * (g2_inf[0][0] - c_ref[0][0]).abs() / c_ref[0][0];
    let ff_cross_err = 100.0 * (ff_06[0][0] - c_ref[0][0]).abs() / c_ref[0][0];

    let out = json!({
        "calibration": "C0.2",
        "reference": "extrapolated_far_field_order2_1.2L_2.4L",
        "reference_fF": {
            "cross_to_cross_fF": c_ref[0][0] * 1e15,
            "claw_to_claw_fF": c_ref[1][1] * 1e15,
            "cross_to_claw_fF": -c_ref[0][1] * 1e15,
            "ground_to_ground_fF": c_ref[2][2] * 1e15,
        },
        "g2_cross_rel_err_pct": g2_cross_err,
        "af3_far_field_0.6L_cross_rel_err_pct": ff_cross_err,
        "elapsed_s": t0.elapsed().as_secs_f64(),
        "methods": methods,
    });

    std::fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;
    let path = out_dir.join("c0_2_open_boundary.json");
    std::fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).map_err(|e| e.to_string())?;
    Ok(out)
}

// --------------------------------------------------------------------------- C0.3: Real Cavity QEP vs 2-Port Driven Taubin Fit

pub fn run_c0_3(out_dir: &Path) -> Result<Value, String> {
    let t0 = Instant::now();
    let mut e = Engine::new();

    // Compact 3D lambda/4 CLT hanger resonator coupled to a 2-port 50-Ohm feedline with axial end-cap lumped_ports
    // Shorted at one end directly to the surrounding ground plane (`lambda/4`), open claw/tip at the other end!
    let (wx, wy) = (2.4e-3, 3.8e-3);
    let (h_sub, h_top, h_bot) = (750.0e-6, 890.0e-6, 900.0e-6);

    call(&mut e, "project_new", json!({ "title": "c03_clt_qep_vs_driven", "dim": 3, "h": 450.0e-6 }))?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "vac_bot", "min": [-wx / 2.0, -wy / 2.0, -h_sub - h_bot], "max": [wx / 2.0, wy / 2.0, -h_sub] }),
    )?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "sub_box", "min": [-wx / 2.0, -wy / 2.0, -h_sub], "max": [wx / 2.0, wy / 2.0, 0.0] }),
    )?;
    call(
        &mut e,
        "add_box",
        json!({ "name": "vac_top", "min": [-wx / 2.0, -wy / 2.0, 0.0], "max": [wx / 2.0, wy / 2.0, h_top] }),
    )?;

    // Ground plane covering [-wx/2, wx/2] x [-wy/2, wy/2] with feedline channel + hanger resonator trench
    // Feedline trace: x in [-620, -560] um (`w = 60 um`), y in [-1750, +1750] um; channel: x in [-655, -525] um, y in [-1850, +1850] um
    // Axial port end-caps: port1 y in [-1850, -1750] um, port2 y in [+1750, +1850] um (touching ground at y = +/-1850 um!)
    // Hanger lambda/4 resonator: coupled parallel section x in [-475, -415] um, y in [-1450, +1550] um (`L ~ 3.0 mm`),
    // shorted into ground plane at y = -1450 um, open at y = +1550 um (trench y in [-1450, +1630] um, x in [-505, -385] um)!
    call(
        &mut e,
        "add_polygons",
        json!({
            "name": "gnd_poly",
            "polygons": [[[-wx / 2.0, -wy / 2.0], [wx / 2.0, -wy / 2.0], [wx / 2.0, wy / 2.0], [-wx / 2.0, wy / 2.0]]],
            "holes": [
                [[-655.0e-6, -1850.0e-6], [-525.0e-6, -1850.0e-6], [-525.0e-6, 1850.0e-6], [-655.0e-6, 1850.0e-6]],
                [[-505.0e-6, -1450.0e-6], [-385.0e-6, -1450.0e-6], [-385.0e-6, 1630.0e-6], [-505.0e-6, 1630.0e-6]]
            ]
        }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "ground_plane", "input": "gnd_poly", "z": 0.0 }))?;

    call(
        &mut e,
        "add_rect",
        json!({ "name": "feed_r", "min": [-620.0e-6, -1750.0e-6], "max": [-560.0e-6, 1750.0e-6] }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "feed_trace", "input": "feed_r", "z": 0.0 }))?;

    call(
        &mut e,
        "add_rect",
        json!({ "name": "p1_r", "min": [-620.0e-6, -1850.0e-6], "max": [-560.0e-6, -1750.0e-6] }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "port1_sheet", "input": "p1_r", "z": 0.0 }))?;

    call(&mut e, "add_rect", json!({ "name": "p2_r", "min": [-620.0e-6, 1750.0e-6], "max": [-560.0e-6, 1850.0e-6] }))?;
    call(&mut e, "add_sheet", json!({ "name": "port2_sheet", "input": "p2_r", "z": 0.0 }))?;

    // lambda/4 hanger trace: touches ground_plane along y = -1450 um (`short`), open at y = +1550 um!
    call(
        &mut e,
        "add_rect",
        json!({ "name": "res_r", "min": [-475.0e-6, -1450.0e-6], "max": [-415.0e-6, 1550.0e-6] }),
    )?;
    call(&mut e, "add_sheet", json!({ "name": "clt_resonator", "input": "res_r", "z": 0.0 }))?;

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
    )?;

    call(&mut e, "add_size", json!({ "sel": "clt_resonator:edge", "h": 140.0e-6 }))?;
    call(&mut e, "add_size", json!({ "sel": "feed_trace:edge", "h": 160.0e-6 }))?;
    call(&mut e, "add_size", json!({ "sel": "(port1_sheet:edge | port2_sheet:edge)", "h": 85.0e-6 }))?;
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 450.0e-6,
            "structured_pml": false,
            "z_coarsen": true,
            "dz_min": 140.0e-6,
            "z_growth_ratio": 2.2,
            "max_z_slabs": 5,
            "max_cells": 65_000,
        }),
    )?;

    let mesh_info = call(&mut e, "mesh_build", json!({ "wait": true }))?;
    let cells = mesh_info["mesh"]["cells"].as_u64().unwrap_or(0);
    eprintln!("[C0.3] mesh built: {cells} tets");

    // 1. Real `ems::cavity` with passive `lumped_port` (`Z_0: 50.0`, `complex_damping: false`) -> `analyse("network")` QEP
    let t_cav0 = Instant::now();
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 1,
            "study_params": { "n_modes": 2, "f_search": 9.0e9 }
        }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "substrate", "material": { "eps_r": 11.45 } }))?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vacuum_top", "material": { "eps_r": 1.0 } }))?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }))?;
    for sh in ["ground_plane", "feed_trace", "clt_resonator"] {
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": sh, "kind": "pec", "params": {} }))?;
    }
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port1_sheet",
            "kind": "lumped_port",
            "params": { "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "complex_damping": false }
        }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port2_sheet",
            "kind": "lumped_port",
            "params": { "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0], "complex_damping": false }
        }),
    )?;
    call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true, "isolate": false, "save_fields": false }))?;
    let cav_s = t_cav0.elapsed().as_secs_f64();

    let net_cav = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": {
                "s21_ports": ["port1", "port2"],
                "resonator_branch": "hanger",
                "touchstone_format": "RI"
            }
        }),
    )?;
    let qep_modes = net_cav["qep_modes"].as_array().ok_or("C0.3: empty qep_modes")?;
    // Pick the hanger resonator mode (moderate kappa: 50 kHz < kappa < 500 MHz, lowest frequency)
    let res_mode = qep_modes
        .iter()
        .filter(|m| {
            let k = m["kappa_Hz"].as_f64().unwrap_or(0.0);
            let f = m["f_re_Hz"].as_f64().unwrap_or(0.0);
            f > 2.0e9 && k > 10.0e3 && k < 0.2 * f
        })
        .min_by(|a, b| a["f_re_Hz"].as_f64().unwrap().total_cmp(&b["f_re_Hz"].as_f64().unwrap()))
        .or_else(|| qep_modes.first())
        .ok_or("C0.3: no physical resonator mode in qep_modes")?;

    let fr_qep = res_mode["f_re_Hz"].as_f64().unwrap();
    let kappa_qep = res_mode["kappa_Hz"].as_f64().unwrap();
    let qc_qep = res_mode["q_c"].as_f64().unwrap_or(0.0);
    eprintln!("[C0.3] cavity QEP done in {cav_s:.2}s: f_r={:.6} GHz, kappa={:.3} MHz", fr_qep * 1e-9, kappa_qep * 1e-6);

    // 2. 2-port `ems::driven` on the EXACT SAME 3D mesh (`cells` unchanged!) -> Probst-Taubin `s21_fit`
    let half_span = (5.0 * kappa_qep).clamp(0.5e6, 150.0e6);
    let t_drv0 = Instant::now();
    call(&mut e, "remove_physics", json!({ "physics": "ems" }))?;
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 1,
            "study_params": {
                "f_min": fr_qep - half_span,
                "f_max": fr_qep + half_span,
                "n_points": 51,
                "prom_tol": 2.0e-4,
                "max_hdm_solves": 6,
            }
        }),
    )?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "substrate", "material": { "eps_r": 11.45 } }))?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vacuum_top", "material": { "eps_r": 1.0 } }))?;
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "vacuum_bot", "material": { "eps_r": 1.0 } }))?;
    for sh in ["ground_plane", "feed_trace", "clt_resonator"] {
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": sh, "kind": "pec", "params": {} }))?;
    }
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port1_sheet",
            "kind": "lumped_port",
            "params": { "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
        }),
    )?;
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port2_sheet",
            "kind": "lumped_port",
            "params": { "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }
        }),
    )?;
    call(&mut e, "solve_run", json!({ "physics": "ems", "wait": true, "isolate": false, "save_fields": false }))?;
    let drv_s = t_drv0.elapsed().as_secs_f64();

    let net_drv = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": {
                "s21_ports": ["port1", "port2"],
                "resonator_branch": "hanger",
                "touchstone_format": "RI"
            }
        }),
    )?;
    let fit = net_drv.get("s21_fit").ok_or("C0.3: driven network analysis returned no s21_fit")?;
    let fr_drv = fit["f_r_Hz"].as_f64().ok_or("C0.3: missing driven f_r_Hz")?;
    let kappa_drv = fit["kappa_Hz"].as_f64().ok_or("C0.3: missing driven kappa_Hz")?;
    let qc_drv = fit["q_c"].as_f64().unwrap_or(0.0);

    let dfr_pct = 100.0 * (fr_drv - fr_qep).abs() / fr_qep;
    let dkappa_pct = 100.0 * (kappa_drv - kappa_qep).abs() / kappa_qep.max(1.0);
    let pass_g0 = dfr_pct <= 0.10 && dkappa_pct <= 2.00;

    let out = json!({
        "calibration": "C0.3",
        "cells": cells,
        "cavity_qep": {
            "f_r_Hz": fr_qep,
            "kappa_Hz": kappa_qep,
            "q_c": qc_qep,
            "seconds": cav_s,
            "all_qep_modes": qep_modes,
        },
        "driven_taubin": {
            "f_r_Hz": fr_drv,
            "kappa_Hz": kappa_drv,
            "q_c": qc_drv,
            "seconds": drv_s,
            "s21_fit": fit,
        },
        "delta_f_r_pct": dfr_pct,
        "delta_kappa_pct": dkappa_pct,
        "gate_g0_thresholds": { "delta_f_r_pct_max": 0.10, "delta_kappa_pct_max": 2.00 },
        "gate_g0_pass": pass_g0,
        "elapsed_s": t0.elapsed().as_secs_f64(),
    });

    std::fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;
    let path = out_dir.join("c0_3_qep_vs_driven.json");
    std::fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).map_err(|e| e.to_string())?;
    Ok(out)
}

pub fn run_calibrations(which: &str, out_dir: &Path) -> Result<Value, String> {
    let w = which.trim().to_ascii_uppercase();
    let mut summary = serde_json::Map::new();
    for (key, fname) in
        [("C0.1", "c0_1_thickness.json"), ("C0.2", "c0_2_open_boundary.json"), ("C0.3", "c0_3_qep_vs_driven.json")]
    {
        if let Ok(txt) = std::fs::read_to_string(out_dir.join(fname)) {
            if let Ok(v) = serde_json::from_str::<Value>(&txt) {
                summary.insert(key.into(), v);
            }
        }
    }
    let mut ran_any = false;
    if w == "C0.1" || w == "ALL" {
        summary.insert("C0.1".into(), run_c0_1(out_dir)?);
        ran_any = true;
    }
    if w == "C0.2" || w == "ALL" {
        summary.insert("C0.2".into(), run_c0_2(out_dir)?);
        ran_any = true;
    }
    if w == "C0.3" || w == "ALL" {
        summary.insert("C0.3".into(), run_c0_3(out_dir)?);
        ran_any = true;
    }
    if !ran_any {
        return Err(format!("unknown calibration '{which}' (expected C0.1, C0.2, C0.3, or all)"));
    }
    let combined = Value::Object(summary);
    let comb_path = out_dir.join("calibrations_summary.json");
    std::fs::write(&comb_path, serde_json::to_string_pretty(&combined).unwrap()).map_err(|e| e.to_string())?;
    Ok(combined)
}
