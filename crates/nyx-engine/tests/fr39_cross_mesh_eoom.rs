// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-39 (ADR 141) end-to-end integration tests for Cross-Mesh Electro-Optic, Photoelastic,
//! and Optomechanical (`g_0 = g_PE + g_MB`) coupling in `nyx-engine`:
//!
//! 1. Cross-mesh Silicon waveguide optomechanical coupling (`sms` on solid-only P2 mesh with
//!    `active_domains: ["core"]` vs `ems` on core+air P2 mesh):
//!    - Agreement within `< 1%` against a direct centered finite-difference mesh-deformation +
//!      photoelastic-permittivity perturbation oracle $\frac{\omega_{\text{opt}}(+\delta x) - \omega_{\text{opt}}(-\delta x)}{2\,\delta x}$,
//!    - Agreement within `< 0.5%` against a fine unified shared mesh while using `>= 2x` fewer mechanical DOFs.
//! 2. End-to-end `Engine` cross-mesh electro-optic (`dfs -> ems`) and photoelastic (`sms -> ems`) chains.
//! 3. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bitwise determinism.

use femcore::optical::C0;
use geom::model::Model;
use nyx_engine::project::{MeshSettings, Project};
use nyx_engine::Engine;
use nyx_physics::{BcSpec, MaterialAssign, Physics, ProblemSpec, PropValue, Scalar};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Builds a 2D Silicon-in-air waveguide project with:
/// - `"clad"`: `[-w_clad/2, +w_clad/2] x [-h_clad/2, +h_clad/2]` (air, `eps_r = 1.0`)
/// - `"core"`: `[-w_core/2, +w_core/2] x [-h_core/2, +h_core/2]` (`silicon`, `n = 3.476`, `p11 = -0.094, p12 = 0.017, p44 = -0.051`)
#[allow(clippy::too_many_arguments)]
fn build_si_waveguide_project(
    w_core: f64,
    h_core: f64,
    w_clad: f64,
    h_clad: f64,
    h_opt: f64,
    h_core_opt: f64,
    adapt_scale_opt: Option<f64>,
    h_mech: Option<f64>,
    freq_opt: f64,
    f_search_mech: f64,
) -> Project {
    let mut m = Model::new_2d();
    let clad = m.rect("clad_r", (-0.5 * w_clad, -0.5 * h_clad), (0.5 * w_clad, 0.5 * h_clad));
    let core = m.rect("core_r", (-0.5 * w_core, -0.5 * h_core), (0.5 * w_core, 0.5 * h_core));
    m.overlay(&[("clad", clad), ("core", core)]);
    m.size("core_fine", geom::select::Selector::parse("region(core)").unwrap(), Some(h_core_opt.into()), None, None);

    let mut proj = Project::new("Si Optomechanical Waveguide", 2);
    proj.geometry = m.list;
    proj.materials.insert(
        "si_optomech".into(),
        r#"
format = 2
name = "si_optomech"
kind = "semiconductor"
point_group = "m-3m"
density = "2329 kg/m^3"

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
"#
        .trim()
        .into(),
    );
    proj.mesh = MeshSettings { h: h_opt, adapt_scale: adapt_scale_opt, ..Default::default() };

    if let Some(hm) = h_mech {
        proj.mesh_by_physics.insert(
            "sms".into(),
            MeshSettings { h: hm, adapt_scale: None, active_domains: Some(vec!["core".into()]), ..Default::default() },
        );
    }

    // Optical `ems::modes` problem
    let mut ems = ProblemSpec::new("ems", 2, "modes");
    ems.study.params = json!({
        "freq": freq_opt,
        "n_modes": 1,
        "n_target": 3.48
    });
    ems.materials.insert("core".into(), MaterialAssign::Named("si_optomech".into()));
    ems.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));

    // Mechanical `sms::eigenfrequency` problem on the Silicon core
    let mut sms = ProblemSpec::new("sms", 2, "eigenfrequency");
    sms.params = if h_mech.is_some() {
        json!({
            "plane": "strain"
        })
    } else {
        json!({
            "plane": "strain",
            "inactive": ["clad"]
        })
    };
    sms.study.params = json!({
        "n_modes": 6,
        "f_search": f_search_mech
    });
    sms.materials.insert("core".into(), MaterialAssign::Named("si_optomech".into()));
    if h_mech.is_none() {
        // On the unified shared mesh, `"clad"` is inactive for mechanics (`inactive: ["clad"]`)
        sms.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    }

    proj.physics = vec![sms, ems];
    proj
}

fn solve_project_sync(proj: Project, targets: &[&str]) -> Engine {
    let mut eng = Engine::new();
    eng.set_project(proj);
    let t: Vec<String> = targets.iter().map(|s| (*s).to_string()).collect();
    let jid = eng.start_solve(Some(t)).expect("start solve");
    eng.wait(jid).expect("wait solve");
    eng
}

#[test]
fn cross_mesh_si_waveguide_optomechanics_vs_fd_oracle_and_unified_mesh() {
    let (w_core, h_core) = (0.45e-6, 0.22e-6);
    let (w_clad, h_clad) = (1.60e-6, 1.10e-6);
    let lambda0 = 1.55e-6;
    let freq_opt = C0 / lambda0;
    // Fundamental lateral breathing frequency of a Si core of width 0.45 um with uz=0:
    // v_L = sqrt(C_11 / rho) = sqrt(165.7e9 / 2329) = 8434.8 m/s => f_1 = v_L / (2 * w_core) = 9.372 GHz
    let f_search_mech = 9.37e9;

    // 1. Cross-mesh solve: `sms` on solid-only P2 mesh (`active_domains: ["core"]`, h = 0.044 um)
    //    and `ems` on core+air adaptive P2 mesh (`h = 0.088 um` in clad, `0.056 um` in core, `adapt_scale = 0.62`).
    let proj_cross = build_si_waveguide_project(
        w_core,
        h_core,
        w_clad,
        h_clad,
        0.088e-6,
        0.056e-6,
        Some(0.62),
        Some(0.044e-6),
        freq_opt,
        f_search_mech,
    );
    // Solve `sms` and `ems` independently (without coupling `sms.S -> ems.strain` during the baseline solve)
    let mut eng_cross = Engine::new();
    eng_cross.set_project(proj_cross);
    let j_sms = eng_cross.start_solve(Some(vec!["sms".into()])).unwrap();
    eng_cross.wait(j_sms).unwrap();
    let j_ems = eng_cross.start_solve(Some(vec!["ems".into()])).unwrap();
    eng_cross.wait(j_ems).unwrap();

    let sms_dofs_cross = eng_cross.result("sms").unwrap().0.dofs;
    let ems_cells_cross = eng_cross.result_model("ems").unwrap().mesh().n_cells();
    let sms_cells_cross = eng_cross.result_model("sms").unwrap().mesh().n_cells();

    // Find the fundamental lateral breathing mode index (elastic mode with f > 5 GHz and largest |g_0_hz|)
    let sms_res = eng_cross.result("sms").unwrap().0.clone();
    let mut breath_idx = 0usize;
    let mut best_g0 = -1.0_f64;
    for idx in 0..sms_res.solutions.len() {
        let f_mode = sms_res.solutions[idx].value.unwrap_or(0.0);
        if f_mode < 5.0e9 {
            continue;
        }
        let r = eng_cross
            .analyse_with(
                "optomechanics",
                "ems",
                &json!({
                    "optical_mode": 1,
                    "mechanical_mode": idx + 1,
                    "length_y_m": 1.0e-6
                }),
            )
            .expect("optomechanics mode scan");
        let g0_abs = r["g_0_hz"].as_f64().unwrap().abs();
        if g0_abs > best_g0 {
            best_g0 = g0_abs;
            breath_idx = idx;
        }
    }
    let f_breath = sms_res.solutions[breath_idx].value.unwrap();
    assert!(
        f_breath > 7.0e9 && f_breath < 11.0e9,
        "expected Si lateral breathing frequency in 7..11 GHz, got {f_breath:.6e} Hz"
    );

    let rep_cross = eng_cross
        .analyse_with(
            "optomechanics",
            "ems",
            &json!({
                "optical_mode": 1,
                "mechanical_mode": breath_idx + 1,
                "length_y_m": 1.0e-6
            }),
        )
        .expect("cross-mesh optomechanics analysis");

    assert_eq!(rep_cross["cross_mesh"], json!(true));
    let g_om_pe_cross = rep_cross["g_om_pe_rad_s_per_m"].as_f64().unwrap();
    let g_om_mb_cross = rep_cross["g_om_mb_rad_s_per_m"].as_f64().unwrap();
    let g_om_tot_cross = rep_cross["g_om_rad_s_per_m"].as_f64().unwrap();
    let g_pe_hz_cross = rep_cross["g_pe_hz"].as_f64().unwrap();
    let g_mb_hz_cross = rep_cross["g_mb_hz"].as_f64().unwrap();
    let g_0_hz_cross = rep_cross["g_0_hz"].as_f64().unwrap();

    // 2. Direct Centered Finite-Difference Oracle for G_MB, G_PE, and Total G_OM = G_PE + G_MB!
    let ems_res_base = eng_cross.result("ems").unwrap().0;
    let n_g = ems_res_base
        .scalars
        .iter()
        .find_map(|(k, v)| match (k.as_str(), v) {
            ("n_g 1", Scalar::Num(x)) => Some(*x),
            _ => None,
        })
        .unwrap();
    let omega_0 = 2.0 * PI * freq_opt;
    let dx = 0.25e-9_f64;

    {
        use femcore::transfer::{SourceMesh, TransferIndex, TransferOptions};

        let s_sol = &sms_res.solutions[breath_idx];
        let u_vals = &s_sol.field("u").unwrap().values;
        let sms_model = eng_cross.result_model("sms").unwrap().clone();
        let ems_model = eng_cross.result_model("ems").unwrap().clone();
        let sms_disc = sms_model.solver_mesh();

        // Compute the exact same mode normalization factor `scale_unit` as `optomechanics.rs`:
        // peak |u| = 1 with positive radial divergence sign `int rho (x - x_cm) . u dV > 0`.
        let mut max_u = 0.0_f64;
        for c in u_vals.chunks(2) {
            max_u = max_u.max(c[0].hypot(c[1]));
        }
        let mut i_div = 0.0_f64;
        for ci in 0..sms_disc.n_cells() {
            let cv = sms_disc.cell(ci);
            let mut xc = [0.0_f64; 2];
            let mut uc = [0.0_f64; 2];
            for &v in cv {
                let p = sms_disc.point(v as usize);
                xc[0] += p[0] / 3.0;
                xc[1] += p[1] / 3.0;
                uc[0] += u_vals[2 * (v as usize)] / 3.0;
                uc[1] += u_vals[2 * (v as usize) + 1] / 3.0;
            }
            i_div += xc[0] * uc[0] + xc[1] * uc[1];
        }
        let sign_flip = if i_div < 0.0 { -1.0 } else { 1.0 };
        let scale_unit = sign_flip / max_u;

        let u_norm: Vec<f64> = u_vals.iter().map(|&v| v * scale_unit).collect();
        let s_raw = &s_sol.field("S").unwrap().values;

        // Build TransferIndex over `sms_disc` to evaluate normalized displacement `u(x, z)` at any point in/on `"core"`
        let sms_dom_names = sms_model.domain_name_map();
        let tidx =
            TransferIndex::new(SourceMesh::Simplex(&sms_disc), &sms_dom_names, TransferOptions::default()).unwrap();

        let ems_spec = eng_cross.project().spec("ems").unwrap().clone();
        let lib = eng_cross.library().unwrap();

        // Helper to deform the optical mesh coordinates smoothly by `amp * w(x, z)` (preserving exact mesh topology!):
        // Inside/on `"core"`, `w(x, z) = u_norm(x, z)`; outside `"core"` in `"clad"`, `w(x, z)` decays smoothly from
        // the closest point `(qx, qz)` on the `"core"` boundary to 0 at the `"clad"` outer wall.
        let deform_optical_mesh = |amp: f64| -> nyx_physics::Mesh {
            let mut m = (**ems_model.mesh()).clone();
            if amp == 0.0 {
                return m;
            }
            if let nyx_physics::Mesh::D2(cm) = &mut m {
                let hx = 0.5 * w_core;
                let hz = 0.5 * h_core;
                let d_buf = 0.22e-6;
                let cx_max = (0.5 * w_clad - hx - d_buf).max(1e-30);
                let cz_max = (0.5 * h_clad - hz - d_buf).max(1e-30);
                for p in &mut cm.points {
                    let qx = p[0].clamp(-hx, hx);
                    let qz = p[1].clamp(-hz, hz);
                    let dx_out = ((p[0] - qx).abs() - d_buf).max(0.0);
                    let dz_out = ((p[1] - qz).abs() - d_buf).max(0.0);
                    let ax = (1.0 - dx_out / cx_max).clamp(0.0, 1.0);
                    let az = (1.0 - dz_out / cz_max).clamp(0.0, 1.0);
                    let decay = ax * az;
                    if let Some(loc) = tidx.locate(&[qx, qz], Some("core")) {
                        let (n_sh, _) = femcore::fem::shape(&loc.coords[..=2], 2);
                        let mut dofs = Vec::new();
                        sms_disc.cell_dofs(loc.cell, 2, &mut dofs);
                        let mut u_pt = [0.0_f64; 2];
                        for (a, &gd) in dofs.iter().enumerate() {
                            u_pt[0] += n_sh[a] * u_norm[2 * gd];
                            u_pt[1] += n_sh[a] * u_norm[2 * gd + 1];
                        }
                        p[0] += amp * decay * u_pt[0];
                        p[1] += amp * decay * u_pt[1];
                    }
                }
            }
            m
        };

        let solve_optical_perturbed = |amp_mb: f64, amp_pe: f64| -> f64 {
            let m_pert = deform_optical_mesh(amp_mb);
            let mut ctx = nyx_physics::Ctx::new(ems_model.as_ref());
            ctx.library = Some(lib.clone());
            if amp_pe != 0.0 {
                let s_scaled: Vec<f64> = s_raw.iter().map(|&v| v * scale_unit * amp_pe).collect();
                let tf = nyx_physics::TransferredField {
                    physics: "sms".into(),
                    field: "S".into(),
                    mesh: sms_model.mesh().clone(),
                    discretisation: Some(sms_model.solver_mesh()),
                    domain_names: sms_model.domain_name_map(),
                    boundary_names: sms_model.boundary_name_map(),
                    order: 2,
                    components: 6,
                    values: s_scaled,
                    cell_fields: Vec::new(),
                };
                ctx.transferred_inputs.insert("strain".into(), tf);
            }
            let r = nyx_ems::Ems.solve(&m_pert, &ems_spec, &mut ctx).unwrap();
            r.solutions[0].eigenvalue.unwrap()[0]
        };

        // 2a. Moving-Boundary Finite-Difference Oracle (`amp_mb = +/- dx, amp_pe = 0`)
        let neff_mb_plus = solve_optical_perturbed(dx, 0.0);
        let neff_mb_minus = solve_optical_perturbed(-dx, 0.0);
        let d_neff_mb_dx = (neff_mb_plus - neff_mb_minus) / (2.0 * dx);
        let g_om_mb_fd = -(omega_0 / n_g) * d_neff_mb_dx;

        let rel_err_mb_fd = (g_om_mb_cross - g_om_mb_fd).abs() / g_om_mb_fd.abs();
        assert!(
            rel_err_mb_fd < 0.01,
            "G_MB cross-mesh ({g_om_mb_cross:.6e}) vs finite-difference mesh-perturbation oracle ({g_om_mb_fd:.6e}): rel_err = {rel_err_mb_fd:.4} >= 1%"
        );

        // 2b. Photoelastic Finite-Difference Oracle (`amp_mb = 0, amp_pe = +/- dx`)
        let neff_pe_plus = solve_optical_perturbed(0.0, dx);
        let neff_pe_minus = solve_optical_perturbed(0.0, -dx);
        let d_neff_pe_dx = (neff_pe_plus - neff_pe_minus) / (2.0 * dx);
        let g_om_pe_fd = -(omega_0 / n_g) * d_neff_pe_dx;

        let rel_err_pe_fd = (g_om_pe_cross - g_om_pe_fd).abs() / g_om_pe_fd.abs();
        assert!(
            rel_err_pe_fd < 0.005,
            "G_PE cross-mesh ({g_om_pe_cross:.6e}) vs finite-difference photoelastic solve oracle ({g_om_pe_fd:.6e}): rel_err = {rel_err_pe_fd:.4} >= 0.5%"
        );

        // 2c. Simultaneous Moving-Boundary + Photoelastic Finite-Difference Oracle (`amp_mb = +/- dx, amp_pe = +/- dx`)
        let neff_tot_plus = solve_optical_perturbed(dx, dx);
        let neff_tot_minus = solve_optical_perturbed(-dx, -dx);
        let d_neff_tot_dx = (neff_tot_plus - neff_tot_minus) / (2.0 * dx);
        let g_om_tot_fd = -(omega_0 / n_g) * d_neff_tot_dx;

        let rel_err_tot_fd = (g_om_tot_cross - g_om_tot_fd).abs() / g_om_tot_fd.abs();
        assert!(
            rel_err_tot_fd < 0.01,
            "Total G_OM = G_PE + G_MB ({g_om_tot_cross:.6e}) vs combined FD oracle ({g_om_tot_fd:.6e}): rel_err = {rel_err_tot_fd:.4} >= 1%"
        );
    }

    // 3. Compare cross-mesh `g_0 = g_PE + g_MB` against a fine unified shared mesh (`h_mech = None`, `h = 0.085 um`, `h_core = 0.054 um`, `adapt_scale = 0.62`)
    //    and verify `< 0.5%` agreement while `sms` on the solid-only mesh uses `>= 2x` fewer cells/DOFs!
    let proj_unified = build_si_waveguide_project(
        w_core,
        h_core,
        w_clad,
        h_clad,
        0.085e-6,
        0.054e-6,
        Some(0.62),
        None,
        freq_opt,
        f_search_mech,
    );
    let mut eng_uni = Engine::new();
    eng_uni.set_project(proj_unified);
    let j1 = eng_uni.start_solve(Some(vec!["sms".into()])).unwrap();
    eng_uni.wait(j1).unwrap();
    let j2 = eng_uni.start_solve(Some(vec!["ems".into()])).unwrap();
    eng_uni.wait(j2).unwrap();

    let uni_mesh_cells = eng_uni.result_model("sms").unwrap().mesh().n_cells();
    let uni_full_mech_dofs = 2 * eng_uni.result_model("sms").unwrap().solver_mesh().n_dofs(2);
    assert!(
        uni_mesh_cells >= 2 * sms_cells_cross && uni_full_mech_dofs >= 2 * sms_dofs_cross,
        "expected unified mesh ({uni_mesh_cells} cells, {uni_full_mech_dofs} full mech DOFs) >= 2x cross-mesh solid-only sms ({sms_cells_cross} cells, {sms_dofs_cross} DOFs); opt cells = {ems_cells_cross}"
    );

    let sms_uni = eng_uni.result("sms").unwrap().0.clone();
    let breath_idx_uni = (0..sms_uni.solutions.len())
        .min_by(|&a, &b| {
            let fa = sms_uni.solutions[a].value.unwrap_or(0.0);
            let fb = sms_uni.solutions[b].value.unwrap_or(0.0);
            (fa - f_breath).abs().partial_cmp(&(fb - f_breath).abs()).unwrap()
        })
        .unwrap();

    let rep_uni = eng_uni
        .analyse_with(
            "optomechanics",
            "ems",
            &json!({
                "optical_mode": 1,
                "mechanical_mode": breath_idx_uni + 1,
                "length_y_m": 1.0e-6
            }),
        )
        .expect("unified shared-mesh optomechanics analysis");
    assert_eq!(rep_uni["cross_mesh"], json!(false));

    let g_pe_hz_uni = rep_uni["g_pe_hz"].as_f64().unwrap();
    let g_mb_hz_uni = rep_uni["g_mb_hz"].as_f64().unwrap();
    let g_0_hz_uni = rep_uni["g_0_hz"].as_f64().unwrap();

    let rel_pe_uni = (g_pe_hz_cross - g_pe_hz_uni).abs() / g_pe_hz_uni.abs();
    let rel_mb_uni = (g_mb_hz_cross - g_mb_hz_uni).abs() / g_mb_hz_uni.abs();
    let rel_g0_uni = (g_0_hz_cross - g_0_hz_uni).abs() / g_0_hz_uni.abs();
    assert!(
        rel_g0_uni < 0.005 && rel_pe_uni < 0.005 && rel_mb_uni < 0.01,
        "cross-mesh vs fine unified mesh disagreement: rel_pe={rel_pe_uni:.4}, rel_mb={rel_mb_uni:.4}, rel_g0={rel_g0_uni:.4}"
    );
}

fn run_eoom_chain_snapshot() -> Value {
    let (w_core, h_core) = (0.45e-6, 0.22e-6);
    let (w_clad, h_clad) = (1.60e-6, 1.10e-6);
    let freq_opt = C0 / 1.55e-6;

    let mut proj = build_si_waveguide_project(
        w_core,
        h_core,
        w_clad,
        h_clad,
        0.16e-6,
        0.065e-6,
        None,
        Some(0.065e-6),
        freq_opt,
        9.37e9,
    );

    // Add `dfs::electrostatic` and give `"core"` both Pockels `r33` and photoelastic `p11, p12, p44`
    let mut dfs = ProblemSpec::new("dfs", 2, "electrostatic");
    dfs.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".into(), 3.9)].into()));
    dfs.materials.insert("core".into(), MaterialAssign::Props([("eps_r".into(), 12.082576)].into()));
    dfs.bcs = vec![
        BcSpec { on: "extreme(-x)".into(), kind: "potential".into(), params: json!({ "value": 15.0 }) },
        BcSpec { on: "extreme(+x)".into(), kind: "ground".into(), params: json!({}) },
    ];

    let mut core_props = BTreeMap::new();
    core_props.insert("eps_r".into(), PropValue::Scalar(12.082576));
    core_props.insert("r33".into(), PropValue::Scalar(30.8e-12));
    core_props.insert("p11".into(), PropValue::Scalar(-0.094));
    core_props.insert("p12".into(), PropValue::Scalar(0.017));
    core_props.insert("p44".into(), PropValue::Scalar(-0.051));
    proj.spec_mut("ems").unwrap().materials.insert("core".into(), MaterialAssign::PropsTensor(core_props));
    proj.physics.insert(0, dfs);

    let eng = solve_project_sync(proj, &["dfs", "sms", "ems"]);
    let neff = eng.result("ems").unwrap().0.solutions[0].eigenvalue.unwrap()[0];
    let sms_sols = &eng.result("sms").unwrap().0.solutions;
    let mech_idx = sms_sols.iter().position(|s| s.value.unwrap_or(0.0) > 5.0e9).unwrap_or(0);
    let rep = eng
        .analyse_with(
            "optomechanics",
            "ems",
            &json!({
                "optical_mode": 1,
                "mechanical_mode": mech_idx + 1,
                "length_y_m": 1.0e-6
            }),
        )
        .unwrap();

    json!({
        "neff_bits": neff.to_bits(),
        "g_pe_bits": rep["g_pe_hz"].as_f64().unwrap().to_bits(),
        "g_mb_bits": rep["g_mb_hz"].as_f64().unwrap().to_bits(),
        "g_0_bits": rep["g_0_hz"].as_f64().unwrap().to_bits(),
        "m_eff_bits": rep["m_eff_kg"].as_f64().unwrap().to_bits(),
    })
}

#[test]
fn cross_mesh_eoom_thread_determinism_1_vs_8() {
    if let Ok(tag) = std::env::var("NYX_FR39_CHILD") {
        let snap = run_eoom_chain_snapshot();
        println!("NYX_FR39_SNAP_{tag}={}", serde_json::to_string(&snap).unwrap());
        return;
    }
    let exe = std::env::current_exe().expect("current test binary");
    let run_with_threads = |nt: usize| -> Value {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("cross_mesh_eoom_thread_determinism_1_vs_8")
            .arg("--nocapture")
            .env("NYX_THREADS", nt.to_string())
            .env("RAYON_NUM_THREADS", nt.to_string())
            .env("NYX_FR39_CHILD", nt.to_string())
            .output()
            .expect("spawn child test");
        assert!(out.status.success(), "child failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let prefix = format!("NYX_FR39_SNAP_{nt}=");
        let line = stdout.lines().find(|l| l.starts_with(&prefix)).expect("find snapshot line");
        serde_json::from_str(&line[prefix.len()..]).unwrap()
    };

    let s1 = run_with_threads(1);
    let s8 = run_with_threads(8);
    assert_eq!(s1, s8, "NYX_THREADS=1 vs 8 produced non-identical bits: {s1} vs {s8}");
}
