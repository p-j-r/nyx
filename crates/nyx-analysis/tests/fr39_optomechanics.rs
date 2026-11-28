// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-39 (ADR 141) verification suite for pure `nyx-analysis` Cross-Mesh `optomechanics`
//! (`optomechanics_with_meshes`, `OptomechanicsAnalysis`) without depending on `nyx-ems` or `nyx-sms`.

use femcore::material::iso_tensor;
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::optical::{C0, EPS0};
use nyx_analysis::{
    optomechanics_with_meshes, Analysis, OptomechanicsAnalysis, OptomechanicsDomainMat, OptomechanicsOptions, HBAR,
};
use nyx_physics::{discretise, CellField, Field, Mesh, Scalar, Solution, SolveResult};
use serde_json::json;
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Builds a 2D two-domain mesh `[-w_clad/2, +w_clad/2] x [-h_clad/2, +h_clad/2]` with core
/// `[-w_core/2, +w_core/2] x [-h_core/2, +h_core/2]` (`domain 1 = "core"`, `domain 2 = "clad"`).
/// When `core_only == true`, only `"core"` cells are emitted (for solid-only mechanical mesh testing).
/// Facets on `x = +/- w_core/2` carry boundary label `10` (`"core_wall"`).
#[allow(clippy::too_many_arguments)]
fn build_rect_waveguide_mesh(
    w_core: f64,
    h_core: f64,
    w_clad: f64,
    h_clad: f64,
    nx_core: usize,
    nz_core: usize,
    nx_pad: usize,
    nz_pad: usize,
    core_only: bool,
) -> Mesh {
    let mut xs = Vec::new();
    if !core_only {
        for i in 0..nx_pad {
            xs.push(-0.5 * w_clad + (0.5 * (w_clad - w_core)) * (i as f64) / (nx_pad as f64));
        }
    }
    let ix_core_lo = xs.len();
    for i in 0..=nx_core {
        xs.push(-0.5 * w_core + w_core * (i as f64) / (nx_core as f64));
    }
    let ix_core_hi = xs.len() - 1;
    if !core_only {
        for i in 1..=nx_pad {
            xs.push(0.5 * w_core + (0.5 * (w_clad - w_core)) * (i as f64) / (nx_pad as f64));
        }
    }

    let mut zs = Vec::new();
    if !core_only {
        for j in 0..nz_pad {
            zs.push(-0.5 * h_clad + (0.5 * (h_clad - h_core)) * (j as f64) / (nz_pad as f64));
        }
    }
    let iz_core_lo = zs.len();
    for j in 0..=nz_core {
        zs.push(-0.5 * h_core + h_core * (j as f64) / (nz_core as f64));
    }
    let iz_core_hi = zs.len() - 1;
    if !core_only {
        for j in 1..=nz_pad {
            zs.push(0.5 * h_core + (0.5 * (h_clad - h_core)) * (j as f64) / (nz_pad as f64));
        }
    }

    let nx = xs.len() - 1;
    let nz = zs.len() - 1;
    let mut cm = CellMesh::<2>::new();
    let vid = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    for &z in &zs {
        for &x in &xs {
            cm.points.push([x, z]);
        }
    }
    for j in 0..nz {
        for i in 0..nx {
            let in_core = i >= ix_core_lo && i < ix_core_hi && j >= iz_core_lo && j < iz_core_hi;
            let dom = if in_core { 1 } else { 2 };
            let (v00, v10, v11, v01) = (vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1));
            cm.push_cell(CellType::Tri, &[v00, v10, v11], dom);
            cm.push_cell(CellType::Tri, &[v00, v11, v01], dom);
        }
    }
    for j in iz_core_lo..iz_core_hi {
        cm.push_facet(&[vid(ix_core_lo, j), vid(ix_core_lo, j + 1)], 10);
        cm.push_facet(&[vid(ix_core_hi, j), vid(ix_core_hi, j + 1)], 10);
    }
    Mesh::D2(cm)
}

/// Builds synthetic optical (`ems`) and mechanical (`sms`) `SolveResult`s with exact polynomial fields
/// so that $m_{\text{eff}}, x_{\text{zpf}}, G_{\text{PE}}, G_{\text{MB}}, g_{\text{PE}}, g_{\text{MB}}, g_0, G_B$
/// have exact closed-form integrals.
fn build_synthetic_results(
    opt_mesh: &Mesh,
    mech_mesh: &Mesh,
    w_core: f64,
    _f_opt: f64,
    f_mech: f64,
    v_g: f64,
    q_mech: f64,
) -> (SolveResult, SolveResult) {
    let opt_disc = discretise(opt_mesh).unwrap();
    let mech_disc = discretise(mech_mesh).unwrap();
    let order = 2usize;
    let nloc_opt = opt_disc.nloc(order);
    let nloc_mech = mech_disc.nloc(order);

    // Optical field: constant E_z = E_0 = 1.0 V/m (tangential to the vertical walls x = +/- w_core/2)
    // and E_x = 0.5 V/m inside core, continuous D_x = eps_1 * 0.5 across x = +/- w_core/2 (`E_x_clad = 0.5 * eps_1 / eps_2`)
    let eps1 = 12.0_f64;
    let eps2 = 1.0_f64;
    let ex_core = 0.5_f64;
    let ex_clad = ex_core * (eps1 / eps2);
    let ez_all = 1.0_f64;

    let nc_opt = opt_disc.n_cells();
    let mut cf_ex = Vec::with_capacity(nc_opt * nloc_opt);
    let mut cf_ey = Vec::with_capacity(nc_opt * nloc_opt);
    let mut cf_ez = Vec::with_capacity(nc_opt * nloc_opt);
    for ci in 0..nc_opt {
        let ex = if opt_disc.cell_labels[ci] == 1 { ex_core } else { ex_clad };
        for _ in 0..nloc_opt {
            cf_ex.push(ex);
            cf_ey.push(0.0);
            cf_ez.push(ez_all);
        }
    }
    let opt_res = SolveResult {
        physics: "ems".into(),
        study: "modes".into(),
        order: order as u8,
        dim: 2,
        dofs: opt_disc.n_dofs(order),
        seconds: 0.0,
        solutions: vec![Solution {
            label: "mode 1".into(),
            value: Some(2.65),
            fields: vec![],
            cell_fields: vec![
                CellField { name: "Ex".into(), unit: "V/m".into(), nloc: nloc_opt, re: cf_ex, im: vec![] },
                CellField { name: "Ey".into(), unit: "V/m".into(), nloc: nloc_opt, re: cf_ey, im: vec![] },
                CellField { name: "Ez".into(), unit: "V/m".into(), nloc: nloc_opt, re: cf_ez, im: vec![] },
            ],
            ..Default::default()
        }],
        scalars: vec![("v_g 1 (m/s)".into(), Scalar::Num(v_g)), ("n_g 1".into(), Scalar::Num(C0 / v_g))],
        log: vec![],
        ..Default::default()
    };

    // Mechanical mode on `mech_mesh`: lateral breathing mode `u_x(x, z) = x / (0.5 * w_core)`, `u_y = u_z = 0`.
    // Peak displacement inside core `[-w_core/2, +w_core/2]` is `max ||u|| = 1.0` at `x = +/- w_core/2`,
    // and strain `S_xx = du_x / dx = 1 / (0.5 * w_core) = 2 / w_core` everywhere in core!
    let coords = mech_disc.dof_coordinates(order);
    let n1_mech = mech_disc.n_dofs(order);
    let mut ux_nodal = Vec::with_capacity(n1_mech);
    let mut uz_nodal = Vec::with_capacity(n1_mech);
    let s_xx = 2.0 / w_core;
    let sxx_nodal = vec![s_xx; n1_mech];
    let zero_nodal = vec![0.0_f64; n1_mech];
    for i in 0..n1_mech {
        let x = coords[2 * i];
        let ux = x / (0.5 * w_core);
        ux_nodal.push(ux);
        uz_nodal.push(0.0);
    }
    let nc_mech = mech_disc.n_cells();
    let mut cf_ux = Vec::with_capacity(nc_mech * nloc_mech);
    let mut cf_uz = Vec::with_capacity(nc_mech * nloc_mech);
    let mut dofs = Vec::new();
    for ci in 0..nc_mech {
        mech_disc.cell_dofs(ci, order, &mut dofs);
        for &gd in dofs.iter().take(nloc_mech) {
            cf_ux.push(ux_nodal[gd]);
            cf_uz.push(0.0);
        }
    }

    let mech_res = SolveResult {
        physics: "sms".into(),
        study: "eigenfrequency".into(),
        order: order as u8,
        dim: 2,
        dofs: 2 * n1_mech,
        seconds: 0.0,
        solutions: vec![Solution {
            label: format!("mode 1: {f_mech:.1} Hz"),
            value: Some(f_mech),
            fields: vec![
                Field::vector("u", "m", &["x", "z"], &[ux_nodal, uz_nodal]),
                Field::vector(
                    "S",
                    "m/m",
                    &["xx", "yy", "zz", "xy", "yz", "xz"],
                    &[
                        sxx_nodal,
                        zero_nodal.clone(),
                        zero_nodal.clone(),
                        zero_nodal.clone(),
                        zero_nodal.clone(),
                        zero_nodal,
                    ],
                ),
            ],
            cell_fields: vec![
                CellField { name: "ux".into(), unit: "m".into(), nloc: nloc_mech, re: cf_ux, im: vec![] },
                CellField { name: "uz".into(), unit: "m".into(), nloc: nloc_mech, re: cf_uz, im: vec![] },
            ],
            ..Default::default()
        }],
        scalars: vec![("f1 (Hz)".into(), Scalar::Num(f_mech)), ("Q1".into(), Scalar::Num(q_mech))],
        log: vec![],
        ..Default::default()
    };

    (opt_res, mech_res)
}

#[test]
fn optomechanics_cross_mesh_vs_exact_closed_form_integrals() {
    let (w_core, h_core) = (0.50e-6, 0.25e-6);
    let (w_clad, h_clad) = (1.50e-6, 0.85e-6);
    let ly = 1.0_f64;
    let f_opt = C0 / 1.55e-6;
    let omega_opt = 2.0 * PI * f_opt;
    let f_mech = 8.0e9;
    let omega_mech = 2.0 * PI * f_mech;
    let v_g = C0 / 4.0;
    let q_mech = 1000.0;
    let rho_si = 2329.0;
    let eps1 = 12.0_f64;
    let eps2 = 1.0_f64;
    let p11 = -0.094_f64;
    let p12 = 0.017_f64;
    let p44 = -0.051_f64;

    let mut pe_mat = [[0.0; 6]; 6];
    pe_mat[0][0] = p11;
    pe_mat[1][1] = p11;
    pe_mat[2][2] = p11;
    pe_mat[0][1] = p12;
    pe_mat[1][0] = p12;
    pe_mat[0][2] = p12;
    pe_mat[2][0] = p12;
    pe_mat[1][2] = p12;
    pe_mat[2][1] = p12;
    pe_mat[3][3] = p44;
    pe_mat[4][4] = p44;
    pe_mat[5][5] = p44;

    // Build optical mesh (core + air cladding, 6x4 core grid + 3x2 pad grid)
    let opt_mesh = build_rect_waveguide_mesh(w_core, h_core, w_clad, h_clad, 6, 4, 3, 2, false);
    // Build distinct solid-only mechanical mesh (`core_only = true`, different grid density 5x3)
    let mech_mesh = build_rect_waveguide_mesh(w_core, h_core, w_clad, h_clad, 5, 3, 0, 0, true);

    assert!(
        mech_mesh.n_cells() * 3 <= opt_mesh.n_cells(),
        "solid-only mechanical mesh ({} cells) should have >= 3x fewer cells than optical core+air mesh ({} cells)",
        mech_mesh.n_cells(),
        opt_mesh.n_cells()
    );

    let (opt_res, mech_res) = build_synthetic_results(&opt_mesh, &mech_mesh, w_core, f_opt, f_mech, v_g, q_mech);

    let mut opts = OptomechanicsOptions {
        optical_freq_hz: Some(f_opt),
        length_m: ly,
        opt_domain_names: BTreeMap::from([(1, "core".into()), (2, "clad".into())]),
        mech_domain_names: BTreeMap::from([(1, "core".into())]),
        opt_boundary_names: BTreeMap::from([(10, "core_wall".into())]),
        mech_boundary_names: BTreeMap::from([(10, "core_wall".into())]),
        ..Default::default()
    };
    opts.domain_materials.insert(
        "core".into(),
        OptomechanicsDomainMat { eps_r: iso_tensor(eps1), photoelastic: Some(pe_mat), rho: rho_si },
    );
    opts.domain_materials
        .insert("clad".into(), OptomechanicsDomainMat { eps_r: iso_tensor(eps2), photoelastic: None, rho: 0.0 });

    let rep = optomechanics_with_meshes(&opt_mesh, &opt_res, &mech_mesh, &mech_res, &opts)
        .expect("cross-mesh optomechanics analysis");
    assert!(rep.cross_mesh, "expected cross_mesh == true for distinct optical and mechanical meshes");

    // Exact closed-form analytical values:
    // 1. m_eff = int_{-h_core/2}^{h_core/2} int_{-w_core/2}^{w_core/2} rho * (x / (w_core / 2))^2 dx dz * ly
    //          = (1 / 3) * rho * w_core * h_core * ly
    let a_core = w_core * h_core;
    let a_clad = w_clad * h_clad - a_core;
    let m_eff_exact = (1.0 / 3.0) * rho_si * a_core * ly;
    let x_zpf_exact = (HBAR / (2.0 * m_eff_exact * omega_mech)).sqrt();
    assert!(
        (rep.m_eff_kg - m_eff_exact).abs() / m_eff_exact < 1e-12,
        "m_eff ({:.6e}) vs exact ({m_eff_exact:.6e})",
        rep.m_eff_kg
    );
    assert!(
        (rep.x_zpf_m - x_zpf_exact).abs() / x_zpf_exact < 1e-12,
        "x_zpf ({:.6e}) vs exact ({x_zpf_exact:.6e})",
        rep.x_zpf_m
    );

    // 2. Optical stored energy U_opt:
    //    In core: eps_1 * (ex_core^2 + ez_all^2) * a_core
    //    In clad: eps_2 * (ex_clad^2 + ez_all^2) * a_clad
    let ex_core = 0.5_f64;
    let ex_clad = ex_core * (eps1 / eps2);
    let ez_all = 1.0_f64;
    let u_opt_exact = EPS0
        * ly
        * (eps1 * (ex_core * ex_core + ez_all * ez_all) * a_core
            + eps2 * (ex_clad * ex_clad + ez_all * ez_all) * a_clad);
    assert!(
        (rep.u_opt_j - u_opt_exact).abs() / u_opt_exact < 1e-12,
        "U_opt ({:.6e}) vs exact ({u_opt_exact:.6e})",
        rep.u_opt_j
    );

    // 3. Cross-mesh photoelastic coupling G_PE and g_PE:
    //    S_xx = 2 / w_core => Delta eps_xx = -eps1^2 * p11 * S_xx, Delta eps_zz = -eps1^2 * p12 * S_xx
    //    pe_integral = EPS0 * ly * a_core * (Delta eps_xx * ex_core^2 + Delta eps_zz * ez_all^2)
    let s_xx = 2.0 / w_core;
    let d_eps_xx = -eps1 * eps1 * p11 * s_xx;
    let d_eps_zz = -eps1 * eps1 * p12 * s_xx;
    let pe_int_exact = EPS0 * ly * a_core * (d_eps_xx * ex_core * ex_core + d_eps_zz * ez_all * ez_all);
    let g_om_pe_exact = -omega_opt * pe_int_exact / (2.0 * u_opt_exact);
    let g_pe_hz_exact = (g_om_pe_exact * x_zpf_exact) / (2.0 * PI);
    assert!(
        (rep.g_pe_hz - g_pe_hz_exact).abs() / g_pe_hz_exact.abs() < 1e-12,
        "cross-mesh g_PE ({:.6e} Hz) vs exact ({g_pe_hz_exact:.6e} Hz)",
        rep.g_pe_hz
    );

    // 4. Cross-mesh moving-boundary coupling G_MB and g_MB:
    //    On the two vertical walls x = +/- w_core/2 (total area = 2 * h_core * ly),
    //    u . n_hat = +1 on both left (u_x = -1, n_x = -1) and right (u_x = +1, n_x = +1) walls!
    //    On the horizontal walls z = +/- h_core/2, u . n_hat = u_z = 0.
    //    On vertical walls: E_par = (0, 0, ez_all), D_perp = EPS0 * eps1 * ex_core = EPS0 * eps2 * ex_clad.
    //    surface_density = EPS0 * (eps1 - eps2) * ez_all^2 - (1/(EPS0*eps1) - 1/(EPS0*eps2)) * D_perp^2
    let d_perp = EPS0 * eps1 * ex_core;
    let surf_density_exact =
        EPS0 * (eps1 - eps2) * ez_all * ez_all - (1.0 / (EPS0 * eps1) - 1.0 / (EPS0 * eps2)) * d_perp * d_perp;
    let mb_int_exact = (2.0 * h_core * ly) * 1.0 * surf_density_exact;
    let g_om_mb_exact = -omega_opt * mb_int_exact / (2.0 * u_opt_exact);
    let g_mb_hz_exact = (g_om_mb_exact * x_zpf_exact) / (2.0 * PI);
    assert!(
        (rep.g_mb_hz - g_mb_hz_exact).abs() / g_mb_hz_exact.abs() < 1e-12,
        "cross-mesh g_MB ({:.6e} Hz) vs exact ({g_mb_hz_exact:.6e} Hz)",
        rep.g_mb_hz
    );

    // 5. Total g_0 and SBS gain G_B:
    let g_0_rad_s_exact = (g_om_pe_exact + g_om_mb_exact) * x_zpf_exact;
    let gb_exact =
        (16.0 * q_mech * g_0_rad_s_exact * g_0_rad_s_exact * ly) / (HBAR * omega_opt * omega_mech * v_g * v_g);
    assert!(
        (rep.g_0_rad_s - g_0_rad_s_exact).abs() / g_0_rad_s_exact.abs() < 1e-12,
        "g_0 ({:.6e}) vs exact ({g_0_rad_s_exact:.6e})",
        rep.g_0_rad_s
    );
    let gb_sim = rep.sbs_gain_w_inv_m_inv.expect("SBS gain G_B");
    assert!((gb_sim - gb_exact).abs() / gb_exact < 1e-12, "G_B ({gb_sim:.6e}) vs exact ({gb_exact:.6e})");

    // 6. Also verify `OptomechanicsAnalysis::run_with_params` trait dispatch on shared mesh
    let (_, mech_res_same) = build_synthetic_results(&opt_mesh, &opt_mesh, w_core, f_opt, f_mech, v_g, q_mech);
    let params = json!({
        "optical_freq_hz": f_opt,
        "length_y_m": ly,
        "mechanical_result": mech_res_same,
        "domain_materials": {
            "1": {
                "eps_r": eps1,
                "p11": p11,
                "p12": p12,
                "p44": p44,
                "rho": rho_si
            },
            "2": {
                "eps_r": eps2,
                "rho": 0.0
            }
        }
    });
    let val = OptomechanicsAnalysis
        .run_with_params(&opt_mesh, &opt_res, &params)
        .expect("OptomechanicsAnalysis::run_with_params");
    assert_eq!(val["cross_mesh"], json!(false));
    let g0_same = val["g_0_hz"].as_f64().unwrap();
    assert!((g0_same - rep.g_0_hz).abs() / rep.g_0_hz.abs() < 1e-12);
}
