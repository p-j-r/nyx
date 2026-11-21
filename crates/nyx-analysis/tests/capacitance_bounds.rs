// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of Prager-Synge complementary capacitance bounds (`AF-4`) in `nyx-analysis`:
//! - Exact legacy JSON preservation when `maxwell_{lower,upper,mid}_matrix_f` are `None`;
//! - Rigorous Loewner propagation `[C_Sigma_low_F, C_Sigma_up_F]` and `[E_C_low_Hz, E_C_up_Hz]` for both 1-pad grounded
//!   (`Xmon`) and 2-pad differential (`TransmonPocket`, Eq. 18) qubits under `"eliminate": "drop"` and
//!   `"eliminate": "floating_schur"` (including off-diagonal sign-shift stress cases);
//! - End-to-end 3D `nyx-dfs` + `femcore::rt::equilibrate_electro_flux_bounds` + `nyx-analysis` transmon qubit bounds for
//!   `order = 1` (`RT_0`) and `order = 2` (`RT_1`), verifying certified `E_C` bracketing and `> 2.5x` midpoint accuracy gain.
#![allow(clippy::needless_range_loop)]

use femcore::meshing::{mesh_model, MeshPlan};
use femcore::rt::{equilibrate_electro_flux_bounds, CellPermittivity, RobinBoundary};
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_analysis::capacitance::capacitance_matrix_with_params;
use nyx_analysis::epr::{E_CHARGE, PLANCK_H};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study, EPS0};
use nyx_physics::{CapacitanceMatrixData, SolveResult};
use serde_json::{json, Value};

fn mutual_from_maxwell(c: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = c.len();
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut s = 0.0;
        for (k, &v) in c[i].iter().enumerate() {
            s += v;
            if i != k {
                m[i][k] = -v;
            }
        }
        m[i][i] = s;
    }
    m
}

fn mat_add_rank1(base: &[Vec<f64>], v: &[f64], scale: f64) -> Vec<Vec<f64>> {
    let n = base.len();
    let mut out = base.to_vec();
    for i in 0..n {
        for j in 0..n {
            out[i][j] += scale * v[i] * v[j];
        }
    }
    out
}

#[test]
fn bounds_absent_preserves_legacy_json_and_present_propagates_1pad_and_2pad_brackets() {
    let f = 1.0e-15;
    // Exact 4-terminal network: c0 = pocket shield, c1 = readout claw, c2 = pad S, c3 = pad B
    let c_exact = vec![
        vec![190.0 * f, -35.0 * f, -55.0 * f, -60.0 * f],
        vec![-35.0 * f, 48.0 * f, -8.5 * f, -0.8 * f],
        vec![-55.0 * f, -8.5 * f, 85.0 * f, -20.0 * f],
        vec![-60.0 * f, -0.8 * f, -20.0 * f, 83.0 * f],
    ];
    // Construct guaranteed PSD perturbations: C^low = C^exact - Delta_low, C^up = C^exact + Delta_up
    let mut c_low = mat_add_rank1(&c_exact, &[1.0, -0.4, 0.8, -0.7], -6.0 * f);
    c_low = mat_add_rank1(&c_low, &[0.3, 0.5, 0.6, 0.5], -4.0 * f);
    let mut c_up = mat_add_rank1(&c_exact, &[0.9, 0.3, -0.7, 0.8], 5.5 * f);
    c_up = mat_add_rank1(&c_up, &[0.4, -0.6, 0.5, 0.5], 3.5 * f);
    let c_mid: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| 0.5 * (c_low[i][j] + c_up[i][j])).collect()).collect();
    let c_err: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| 0.5 * ((c_up[i][i] - c_low[i][i]) * (c_up[j][j] - c_low[j][j])).sqrt()).collect())
        .collect();

    // 1. Legacy payload (bounds = None) omits all bound keys
    let cap_legacy = CapacitanceMatrixData {
        terminals: vec!["c0".into(), "c1".into(), "c2".into(), "c3".into()],
        maxwell_matrix_f: c_exact.clone(),
        mutual_matrix_f: mutual_from_maxwell(&c_exact),
        ec_matrix_hz: vec![vec![0.0; 4]; 4],
        ..Default::default()
    };
    let res_legacy = SolveResult { scalars: vec![cap_legacy.scalar()], ..Default::default() };
    let out_legacy_base = capacitance_matrix_with_params(&res_legacy, &Value::Null).unwrap();
    assert_eq!(out_legacy_base.as_object().unwrap().len(), 4);
    let out_legacy_q = capacitance_matrix_with_params(
        &res_legacy,
        &json!({
            "ground": "c0",
            "qubits": [
                { "name": "Q_diff", "pads": ["c2", "c3"] },
                { "name": "Q_gnd", "pads": ["c1"] }
            ]
        }),
    )
    .unwrap();
    assert!(out_legacy_q.get("maxwell_lower_matrix_F").is_none());
    assert!(out_legacy_q["qubits"][0].get("C_Sigma_low_F").is_none());
    assert!(out_legacy_q["qubits"][0].get("E_C_low_Hz").is_none());
    assert!(out_legacy_q["qubits"][1].get("C_Sigma_low_F").is_none());

    // 2. Payload with Prager-Synge bounds propagates exact brackets under BOTH `drop` and `floating_schur`
    let cap_bounded = CapacitanceMatrixData {
        terminals: vec!["c0".into(), "c1".into(), "c2".into(), "c3".into()],
        maxwell_matrix_f: c_mid.clone(),
        mutual_matrix_f: mutual_from_maxwell(&c_mid),
        ec_matrix_hz: vec![vec![0.0; 4]; 4],
        maxwell_lower_matrix_f: Some(c_low),
        maxwell_upper_matrix_f: Some(c_up),
        maxwell_mid_matrix_f: Some(c_mid),
        maxwell_error_bound_f: Some(c_err),
    };
    let res_bounded = SolveResult { scalars: vec![cap_bounded.scalar()], ..Default::default() };
    let out_base = capacitance_matrix_with_params(&res_bounded, &Value::Null).unwrap();
    assert!(out_base.get("maxwell_lower_matrix_F").is_some());
    assert!(out_base.get("maxwell_upper_matrix_F").is_some());
    assert!(out_base.get("maxwell_mid_matrix_F").is_some());
    assert!(out_base.get("maxwell_error_bound_F").is_some());

    for elim in ["drop", "floating_schur"] {
        let query = json!({
            "ground": "c0",
            "eliminate": elim,
            "qubits": [
                { "name": "Q_diff", "pads": ["c2", "c3"], "E_J_Hz": 16.0e9 },
                { "name": "Q_gnd", "pads": ["c1"], "E_J_Hz": 14.0e9 }
            ]
        });
        let exact_out = capacitance_matrix_with_params(&res_legacy, &query).unwrap();
        let bound_out = capacitance_matrix_with_params(&res_bounded, &query).unwrap();
        assert!(bound_out.get("maxwell_reduced_lower_F").is_some());
        assert!(bound_out.get("maxwell_reduced_upper_F").is_some());

        for qi in 0..2 {
            let cs_exact = exact_out["qubits"][qi]["C_Sigma_F"].as_f64().unwrap();
            let ec_exact = exact_out["qubits"][qi]["E_C_Hz"].as_f64().unwrap();
            let q = &bound_out["qubits"][qi];
            let cs_low = q["C_Sigma_low_F"].as_f64().unwrap();
            let cs_up = q["C_Sigma_up_F"].as_f64().unwrap();
            let ec_low = q["E_C_low_Hz"].as_f64().unwrap();
            let ec_up = q["E_C_up_Hz"].as_f64().unwrap();

            assert!(
                cs_low < cs_exact && cs_exact < cs_up,
                "elim={elim} q={qi}: expected cs_low ({cs_low:e}) < cs_exact ({cs_exact:e}) < cs_up ({cs_up:e})"
            );
            assert!(
                ec_low < ec_exact && ec_exact < ec_up,
                "elim={elim} q={qi}: expected ec_low ({ec_low:e}) < ec_exact ({ec_exact:e}) < ec_up ({ec_up:e})"
            );
        }
    }

    // 3. Weak-coupling stress test where off-diagonal `C_sb^low > 0` (turns slightly positive while `C^low` stays SPD!)
    let c_weak_exact = vec![vec![60.0 * f, -f], vec![-f, 62.0 * f]];
    let c_weak_low = vec![vec![55.0 * f, 1.5 * f], vec![1.5 * f, 57.0 * f]]; // C_01^low = +1.5 fF > 0, C^exact - C^low = [[5, -2.5], [-2.5, 5]] > 0!
    let c_weak_up = vec![vec![65.0 * f, -3.5 * f], vec![-3.5 * f, 67.0 * f]]; // C^up - C^exact = [[5, -2.5], [-2.5, 5]] > 0!
    let res_weak_ex = SolveResult {
        scalars: vec![CapacitanceMatrixData {
            terminals: vec!["s".into(), "b".into()],
            maxwell_matrix_f: c_weak_exact.clone(),
            mutual_matrix_f: mutual_from_maxwell(&c_weak_exact),
            ec_matrix_hz: vec![vec![0.0; 2]; 2],
            ..Default::default()
        }
        .scalar()],
        ..Default::default()
    };
    let res_weak_bd = SolveResult {
        scalars: vec![CapacitanceMatrixData {
            terminals: vec!["s".into(), "b".into()],
            maxwell_matrix_f: c_weak_exact,
            mutual_matrix_f: vec![vec![0.0; 2]; 2],
            ec_matrix_hz: vec![vec![0.0; 2]; 2],
            maxwell_lower_matrix_f: Some(c_weak_low),
            maxwell_upper_matrix_f: Some(c_weak_up),
            ..Default::default()
        }
        .scalar()],
        ..Default::default()
    };
    let q_weak = json!({ "qubits": [{ "name": "Qw", "pads": ["s", "b"] }] });
    let ex_w = capacitance_matrix_with_params(&res_weak_ex, &q_weak).unwrap();
    let bd_w = capacitance_matrix_with_params(&res_weak_bd, &q_weak).unwrap();
    let ec_ex_w = ex_w["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let ec_lo_w = bd_w["qubits"][0]["E_C_low_Hz"].as_f64().unwrap();
    let ec_hi_w = bd_w["qubits"][0]["E_C_up_Hz"].as_f64().unwrap();
    assert!(ec_lo_w < ec_ex_w && ec_ex_w < ec_hi_w);
}

fn solve_3d_transmon_with_flux_bounds(global_h: f64, edge_h: f64, order: usize, use_midpoint: bool) -> SolveResult {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-180e-6, -180e-6, -100e-6), (180e-6, 180e-6, 0.0));
    let vac = m.box_("vac", (-180e-6, -180e-6, 0.0), (180e-6, 180e-6, 100e-6));

    // Coplanar ground frame `gnd` + two differential pads `pad_s`, `pad_b` at z = 0
    let gnd_outer = vec![pt(-160e-6, -160e-6), pt(160e-6, -160e-6), pt(160e-6, 160e-6), pt(-160e-6, 160e-6)];
    let gnd_hole = vec![pt(-115e-6, -85e-6), pt(115e-6, -85e-6), pt(115e-6, 85e-6), pt(-115e-6, 85e-6)];
    let pads_outer = vec![pt(-95e-6, -60e-6), pt(-15e-6, -60e-6), pt(-15e-6, 60e-6), pt(-95e-6, 60e-6)];
    let padb_outer = vec![pt(15e-6, -60e-6), pt(95e-6, -60e-6), pt(95e-6, 60e-6), pt(15e-6, 60e-6)];

    let g_poly = m.polygons("g_poly", vec![gnd_outer], vec![gnd_hole]);
    let s_poly = m.polygons("s_poly", vec![pads_outer], vec![]);
    let b_poly = m.polygons("b_poly", vec![padb_outer], vec![]);

    m.sheet("c0", g_poly, 0.0);
    m.sheet("c1", s_poly, 0.0);
    m.sheet("c2", b_poly, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for name in ["c0", "c1", "c2"] {
        m.size(
            &format!("sz_{name}"),
            Selector::parse(&format!("{name}:edge")).unwrap(),
            Some(Expr::Num(edge_h)),
            None,
            None,
        );
    }

    let gm = mesh_model(&m.list, &MeshPlan::simple(global_h)).unwrap();
    let sub_dom = gm.domain("substrate").unwrap();
    let vac_dom = gm.domain("vacuum").unwrap();
    let c0_labs = gm.boundary("c0");
    let c1_labs = gm.boundary("c1");
    let c2_labs = gm.boundary("c2");
    let outer_labs = gm.select("outer").unwrap();

    let mut prob = Problem::new(Study::Electrostatic, order);
    prob.materials.set(sub_dom, "eps_r", 11.45);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.bcs = vec![
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
        Bc { labels: c0_labs.clone(), kind: BcKind::Terminal("c0".into()) },
        Bc { labels: c1_labs.clone(), kind: BcKind::Terminal("c1".into()) },
        Bc { labels: c2_labs.clone(), kind: BcKind::Terminal("c2".into()) },
    ];
    let mut res = solve_electrostatics(&gm.mesh, &prob).unwrap();

    let eps_cells: Vec<f64> =
        gm.mesh.cell_labels.iter().map(|&d| if d == sub_dom { 11.45 * EPS0 } else { EPS0 }).collect();
    let v_all: Vec<Vec<f64>> = res.solutions.iter().map(|sol| sol.field("V").unwrap().values.clone()).collect();
    let mut dir_labs = outer_labs;
    dir_labs.extend(c0_labs);
    dir_labs.extend(c1_labs);
    dir_labs.extend(c2_labs);
    let mut cap = CapacitanceMatrixData::of(&res).unwrap();
    let eq = equilibrate_electro_flux_bounds(
        &gm.mesh,
        order,
        &v_all,
        CellPermittivity::Scalar(&eps_cells),
        &dir_labs,
        RobinBoundary::None,
        Some(&cap.maxwell_matrix_f),
    )
    .unwrap();

    cap.maxwell_lower_matrix_f = Some(eq.c_lower_matrix);
    cap.maxwell_upper_matrix_f = Some(eq.c_upper_matrix);
    cap.maxwell_mid_matrix_f = Some(eq.c_mid_matrix.clone());
    cap.maxwell_error_bound_f = Some(eq.c_error_bound_matrix);
    if use_midpoint {
        cap.maxwell_matrix_f = eq.c_mid_matrix.clone();
        cap.mutual_matrix_f = mutual_from_maxwell(&eq.c_mid_matrix);
    }
    let new_scalar = cap.scalar();
    if let Some(slot) = res.scalars.iter_mut().find(|(k, _)| k == &new_scalar.0) {
        *slot = new_scalar;
    }
    res
}

#[test]
fn end_to_end_3d_transmon_pocket_qubit_ec_is_bracketed_and_midpoint_cuts_error_over_2_5x() {
    let query = json!({
        "ground": "c0",
        "qubits": [
            { "name": "Q_pocket", "pads": ["c1", "c2"], "E_J_Hz": 18.0e9 },
            { "name": "Q_single", "pads": ["c1"] }
        ]
    });

    // Converged P2 midpoint reference on fine sheet3d mesh
    let res_ref = solve_3d_transmon_with_flux_bounds(24e-6, 4.5e-6, 2, true);
    let ref_json = capacitance_matrix_with_params(&res_ref, &query).unwrap();
    let ec_ref_diff = ref_json["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let cs_ref_diff = ref_json["qubits"][0]["C_Sigma_F"].as_f64().unwrap();
    let ec_ref_1pad = ref_json["qubits"][1]["E_C_Hz"].as_f64().unwrap();

    for order in [1usize, 2usize] {
        let res_primal = solve_3d_transmon_with_flux_bounds(60e-6, 18e-6, order, false);
        let res_mid = solve_3d_transmon_with_flux_bounds(60e-6, 18e-6, order, true);
        let out_primal = capacitance_matrix_with_params(&res_primal, &query).unwrap();
        let out_mid = capacitance_matrix_with_params(&res_mid, &query).unwrap();

        // Differential 2-pad qubit bounds
        let q0 = &out_mid["qubits"][0];
        let cs_low = q0["C_Sigma_low_F"].as_f64().unwrap();
        let cs_up = q0["C_Sigma_up_F"].as_f64().unwrap();
        let ec_low = q0["E_C_low_Hz"].as_f64().unwrap();
        let ec_up = q0["E_C_up_Hz"].as_f64().unwrap();
        assert!(
            cs_low < cs_ref_diff && cs_ref_diff < cs_up,
            "order={order} 2-pad C_Sigma not bracketed: [{cs_low:e}, {cs_up:e}] vs ref {cs_ref_diff:e}"
        );
        assert!(
            ec_low < ec_ref_diff && ec_ref_diff < ec_up,
            "order={order} 2-pad E_C not bracketed: [{ec_low:e}, {ec_up:e}] vs ref {ec_ref_diff:e}"
        );

        // Check E_C = e^2 / (2 h C_Sigma) exact reciprocal relationship on bounds
        let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
        assert!(((scale / cs_up) - ec_low).abs() < 1.0e-12 * ec_low);
        assert!(((scale / cs_low) - ec_up).abs() < 1.0e-12 * ec_up);

        // Grounded 1-pad qubit bounds
        let q1 = &out_mid["qubits"][1];
        let ec1_low = q1["E_C_low_Hz"].as_f64().unwrap();
        let ec1_up = q1["E_C_up_Hz"].as_f64().unwrap();
        assert!(
            ec1_low < ec_ref_1pad && ec_ref_1pad < ec1_up,
            "order={order} 1-pad E_C not bracketed: [{ec1_low:e}, {ec1_up:e}] vs ref {ec_ref_1pad:e}"
        );

        if order == 2 {
            // Midpoint `E_C` accuracy gain `> 2.5x` over one-sided primal `E_C` at P2 (`RT_1`)
            let ec_primal = out_primal["qubits"][0]["E_C_Hz"].as_f64().unwrap();
            let ec_mid = q0["E_C_Hz"].as_f64().unwrap();
            let err_primal = (ec_primal - ec_ref_diff).abs();
            let err_mid = (ec_mid - ec_ref_diff).abs();
            let gain = err_primal / err_mid.max(1.0e-30);
            eprintln!("3D TransmonPocket P2/RT_1 qubit E_C midpoint gain = {gain:.2}x (primal={:.2} MHz, mid={:.2} MHz, ref={:.2} MHz, bracket=[{:.2}, {:.2}] MHz)",
                ec_primal * 1e-6, ec_mid * 1e-6, ec_ref_diff * 1e-6, ec_low * 1e-6, ec_up * 1e-6);
            assert!(
                gain > 2.5,
                "order={order}: expected >2.5x qubit E_C accuracy gain from C^mid, got {gain:.2}x (primal err={err_primal:e}, mid err={err_mid:e})"
            );
        }
    }
}
