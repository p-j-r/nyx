// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analytical verification of ground/Schur elimination, single-pad (`Xmon`) and differential two-pad (`TransmonPocket`,
//! Eqs. 17–18 of `arXiv:2611.22057v1`) lumped cQED reduction, resonator ZPF coupling (`g_Hz`, `chi_Hz`), and gauge
//! null-mode handling in `nyx_analysis::capacitance::capacitance_matrix_with_params` (QFR-6).

use nyx_analysis::capacitance::capacitance_matrix_with_params;
use nyx_analysis::epr::{E_CHARGE, PLANCK_H};
use nyx_analysis::{Analysis, CapacitanceMatrixAnalysis};
use nyx_physics::{CapacitanceMatrixData, Mesh, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;

fn dummy_mesh() -> Mesh {
    let sm = femcore::testkit::rect(1.0, 0.5, 2, 2);
    Mesh::from_solver_mesh(&sm).unwrap()
}

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

fn make_result(terminals: &[&str], maxwell_f: Vec<Vec<f64>>) -> SolveResult {
    let n = terminals.len();
    let mutual_matrix_f = mutual_from_maxwell(&maxwell_f);
    let ec_matrix_hz = vec![vec![0.0; n]; n];
    let cap = CapacitanceMatrixData {
        terminals: terminals.iter().map(|&s| s.to_string()).collect(),
        maxwell_matrix_f: maxwell_f,
        mutual_matrix_f,
        ec_matrix_hz,
        ..Default::default()
    };
    SolveResult {
        physics: "dfs".into(),
        study: "electrostatic".into(),
        scalars: vec![cap.scalar()],
        ..Default::default()
    }
}

#[inline]
fn rel_err(got: f64, want: f64) -> f64 {
    (got - want).abs() / want.abs().max(1.0e-30)
}

#[test]
fn empty_and_null_params_preserve_exact_legacy_json() {
    let mesh = dummy_mesh();
    let c = vec![vec![120.0e-15, -35.0e-15], vec![-35.0e-15, 95.0e-15]];
    let res = make_result(&["p1", "p2"], c);

    let base = CapacitanceMatrixAnalysis.run(&mesh, &res).unwrap();
    let via_null = capacitance_matrix_with_params(&res, &Value::Null).unwrap();
    let via_empty = capacitance_matrix_with_params(&res, &json!({})).unwrap();
    let via_trait = CapacitanceMatrixAnalysis.run_with_params(&mesh, &res, &json!({})).unwrap();

    assert_eq!(base, via_null);
    assert_eq!(base, via_empty);
    assert_eq!(base, via_trait);
    let keys: Vec<&str> = via_empty.as_object().unwrap().keys().map(String::as_str).collect();
    assert_eq!(keys.len(), 4);
    assert!(via_empty.get("reduced_terminals").is_none());
}

#[test]
fn transmon_pocket_2pad_and_coupler_resonators_match_eq18_to_machine_precision() {
    // Realistic 5-terminal TransmonPocket (c0 = ground, c1 = readout claw, c2 = pad S, c3 = pad B, c4 = bus coupler)
    // grounded-box Maxwell matrix (in Farads, modeled after SQuADDS q3d_cap_0470 / q3d_cap_0361).
    let f = 1.0e-15;
    let c_m = vec![
        vec![195.40 * f, -40.18 * f, -60.37 * f, -65.90 * f, -16.20 * f],
        vec![-40.18 * f, 50.56 * f, -9.53 * f, -0.64 * f, -0.05 * f],
        vec![-60.37 * f, -9.53 * f, 86.8832807652342 * f, -24.3704459320443 * f, -1.12 * f], // note row sum > 0 adjusted below
        vec![-65.90 * f, -0.64 * f, -24.3704459320443 * f, 84.0494025767063 * f, -4.85 * f],
        vec![-16.20 * f, -0.05 * f, -1.12 * f, -4.85 * f, 24.50 * f],
    ];
    // Adjust c0, c2, c3 diagonals so every conductor has strictly positive shunt capacitance to box ground (strictly diagonally dominant SPD)
    let mut c = c_m;
    c[2][0] = -51.40 * f;
    c[0][2] = -51.40 * f;
    c[3][0] = -52.80 * f;
    c[0][3] = -52.80 * f;

    let res = make_result(&["c0", "c1", "c2", "c3", "c4"], c.clone());
    let e_j_hz = 18.0e9;
    let f_r1_hz = 6.8e9;
    let z0_r1 = 50.0;
    let f_r2_hz = 7.4e9;
    let z0_r2 = 55.0;
    let f_r3_hz = 5.2e9;
    let c_r3_lumped = 420.0e-15;

    let out = capacitance_matrix_with_params(
        &res,
        &json!({
            "ground": "c0",
            "qubits": [{ "name": "Q1", "pads": ["c2", "c3"], "E_J_Hz": e_j_hz }],
            "resonators": [
                { "name": "R_ro", "terminal": "c1", "f_r_Hz": f_r1_hz, "Z0_ohm": z0_r1, "resonator_type": "quarter_wave" },
                { "name": "R_bus", "terminal": "c4", "f_r_Hz": f_r2_hz, "Z0_ohm": z0_r2, "resonator_type": "half_wave" },
                { "name": "R_lump", "terminal": "c1", "f_r_Hz": f_r3_hz, "resonator_type": "lumped", "C_r_F": c_r3_lumped }
            ]
        }),
    )
    .unwrap();

    assert_eq!(out["eliminate"], json!("drop"));
    assert_eq!(out["reduced_terminals"], json!(["c1", "c2", "c3", "c4"]));
    assert_eq!(out["has_null_mode"], json!(false));
    assert_eq!(out["null_modes"], json!(0));
    let cond = out["condition_number"].as_f64().unwrap();
    assert!(cond > 1.0 && cond < 100.0, "unexpected condition_number = {cond}");

    // Analytical Eq. (18) check
    let c_ss = c[2][2];
    let c_bb = c[3][3];
    let c_sb = c[2][3]; // negative Maxwell off-diagonal
    assert!(c_sb < 0.0);
    let denom = c_ss + c_bb + 2.0 * c_sb;
    let want_c_sigma = (c_ss * c_bb - c_sb * c_sb) / denom;
    let want_ec_hz = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * want_c_sigma);
    let want_ceff_c1 = ((c_bb + c_sb) * c[2][1].abs() - (c_ss + c_sb) * c[3][1].abs()).abs() / denom;
    let want_ceff_c4 = ((c_bb + c_sb) * c[2][4].abs() - (c_ss + c_sb) * c[3][4].abs()).abs() / denom;
    let want_beta_c1 = want_ceff_c1 / want_c_sigma;
    let want_beta_c4 = want_ceff_c4 / want_c_sigma;

    // Cross-check Eq. (18) against 2x2 pad-submatrix inverse differential charging energy
    let det2 = c_ss * c_bb - c_sb * c_sb;
    let inv_diff = (c_bb + c_ss + 2.0 * c_sb) / det2; // (1, -1) [[c_bb, -c_sb], [-c_sb, c_ss]] (1, -1)^T / det2
    let ec_2pad_inv = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H) * inv_diff;
    assert!(rel_err(want_ec_hz, ec_2pad_inv) < 1.0e-15);

    let q1 = &out["qubits"][0];
    assert_eq!(q1["name"], json!("Q1"));
    assert_eq!(q1["pads"], json!(["c2", "c3"]));
    assert!(rel_err(q1["C_Sigma_F"].as_f64().unwrap(), want_c_sigma) < 1.0e-14);
    assert!(rel_err(q1["E_C_Hz"].as_f64().unwrap(), want_ec_hz) < 1.0e-14);
    assert!(rel_err(q1["C_eff_F"]["c1"].as_f64().unwrap(), want_ceff_c1) < 1.0e-14);
    assert!(rel_err(q1["C_eff_F"]["c4"].as_f64().unwrap(), want_ceff_c4) < 1.0e-14);
    assert!(rel_err(q1["beta"]["c1"].as_f64().unwrap(), want_beta_c1) < 1.0e-14);
    assert!(rel_err(q1["beta"]["c4"].as_f64().unwrap(), want_beta_c4) < 1.0e-14);
    assert!(q1["C_eff_F"].get("c0").is_none(), "ground terminal c0 must not appear as a coupler");

    let want_ej_ec = e_j_hz / want_ec_hz;
    let want_f01 = (8.0 * e_j_hz * want_ec_hz).sqrt() - want_ec_hz;
    let want_alpha = -want_ec_hz;
    let want_nzpf = (e_j_hz / (32.0 * want_ec_hz)).powf(0.25);
    assert!(rel_err(q1["ej_over_ec"].as_f64().unwrap(), want_ej_ec) < 1.0e-14);
    assert!(rel_err(q1["f01_Hz"].as_f64().unwrap(), want_f01) < 1.0e-14);
    assert!(rel_err(q1["anharmonicity_Hz"].as_f64().unwrap(), want_alpha) < 1.0e-14);
    assert!(rel_err(q1["n_zpf"].as_f64().unwrap(), want_nzpf) < 1.0e-14);

    // Resonator & dispersive coupling checks (quarter_wave, half_wave, lumped)
    let hbar = PLANCK_H / (2.0 * PI);
    let check_pair = |idx: usize, f_r: f64, want_cr: f64, want_beta: f64| {
        let r_obj = &out["resonators"][idx];
        let c_obj = &out["couplings"][idx];
        let omega_r = 2.0 * PI * f_r;
        let want_vrms = (hbar * omega_r / (2.0 * want_cr)).sqrt();
        let want_g = 2.0 * E_CHARGE * want_beta * want_nzpf * want_vrms / PLANCK_H;
        let delta = want_f01 - f_r;
        let want_chi = -(want_g * want_g) / delta * (want_ec_hz / (delta - want_ec_hz));

        assert!(rel_err(r_obj["C_r_F"].as_f64().unwrap(), want_cr) < 1.0e-14);
        assert!(rel_err(r_obj["V_rms_V"].as_f64().unwrap(), want_vrms) < 1.0e-14);
        assert!(rel_err(c_obj["beta"].as_f64().unwrap(), want_beta) < 1.0e-14);
        assert!(rel_err(c_obj["n_zpf"].as_f64().unwrap(), want_nzpf) < 1.0e-14);
        assert!(rel_err(c_obj["g_Hz"].as_f64().unwrap(), want_g) < 1.0e-14);
        assert!(rel_err(c_obj["chi_Hz"].as_f64().unwrap(), want_chi) < 1.0e-14);
    };

    check_pair(0, f_r1_hz, PI / (4.0 * 2.0 * PI * f_r1_hz * z0_r1), want_beta_c1);
    check_pair(1, f_r2_hz, PI / (2.0 * 2.0 * PI * f_r2_hz * z0_r2), want_beta_c4);
    check_pair(2, f_r3_hz, c_r3_lumped, want_beta_c1);
}

#[test]
fn xmon_1pad_grounded_matches_eq17_to_machine_precision() {
    // 3-terminal Xmon from SQuADDS q3d_cap_0554: c0 = G, c1 = N (claw), c2 = S (cross)
    let f = 1.0e-15;
    let c = vec![
        vec![266.2751602122418 * f, -67.06011354094588 * f, -172.31834463911892 * f],
        vec![-67.06011354094588 * f, 69.17479232898475 * f, -1.983533590526942 * f],
        vec![-172.31834463911892 * f, -1.983533590526942 * f, 175.60664203508705 * f],
    ];
    let res = make_result(&["c0", "c1", "c2"], c.clone());
    let e_j_hz = 15.0e9;
    let f_r_hz = 6.5e9;

    let out = capacitance_matrix_with_params(
        &res,
        &json!({
            "ground": ["c0"],
            "qubits": [{ "name": "X1", "pads": ["c2"], "E_J_Hz": e_j_hz }],
            "resonators": [{ "name": "R1", "terminal": "c1", "f_r_Hz": f_r_hz }]
        }),
    )
    .unwrap();

    let want_c_sigma = c[2][2];
    let want_ec = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * want_c_sigma);
    let want_ceff = c[1][2].abs();
    let want_beta = want_ceff / want_c_sigma;

    let q = &out["qubits"][0];
    assert_eq!(q["topology"], json!("grounded_1pad"));
    assert!(rel_err(q["C_Sigma_F"].as_f64().unwrap(), want_c_sigma) < 1.0e-14);
    assert!(rel_err(q["E_C_Hz"].as_f64().unwrap(), want_ec) < 1.0e-14);
    assert!(rel_err(q["C_eff_F"]["c1"].as_f64().unwrap(), want_ceff) < 1.0e-14);
    assert!(rel_err(q["beta"]["c1"].as_f64().unwrap(), want_beta) < 1.0e-14);

    let coup = &out["couplings"][0];
    assert!(coup["g_Hz"].as_f64().unwrap() > 5.0e6);
    assert!(coup["chi_Hz"].as_f64().unwrap() < 0.0);
}

#[test]
fn floating_schur_and_null_mode_detection_work_on_grounded_and_ungrounded_networks() {
    let f = 1.0e-15;
    // 1) Closed ungrounded 4-conductor network: exact Kirchhoff row sums == 0 (constant-V null mode)
    // Mutuals: c01 = 40, c02 = 60, c03 = 65, c12 = 10, c13 = 1, c23 = 25 (fF)
    let c_closed = vec![
        vec![165.0 * f, -40.0 * f, -60.0 * f, -65.0 * f],
        vec![-40.0 * f, 51.0 * f, -10.0 * f, -f],
        vec![-60.0 * f, -10.0 * f, 95.0 * f, -25.0 * f],
        vec![-65.0 * f, -f, -25.0 * f, 91.0 * f],
    ];
    let res_closed = make_result(&["c0", "c1", "c2", "c3"], c_closed.clone());

    // Without dropping ground: detects exact 1 gauge null mode, finite capped condition_number, finite Moore-Penrose ec_reduced_Hz
    let raw_check = capacitance_matrix_with_params(&res_closed, &json!({ "eliminate": "drop" })).unwrap();
    assert_eq!(raw_check["has_null_mode"], json!(true));
    assert_eq!(raw_check["null_modes"], json!(1));
    assert!(raw_check["condition_number"].as_f64().unwrap() >= 1.0e10);

    // Dropping c0 fixes V_c0 = 0 and eliminates the null mode completely
    let dropped = capacitance_matrix_with_params(&res_closed, &json!({ "ground": "c0", "eliminate": "drop" })).unwrap();
    assert_eq!(dropped["has_null_mode"], json!(false));
    assert_eq!(dropped["null_modes"], json!(0));
    assert!(dropped["condition_number"].as_f64().unwrap() < 20.0);

    // Verify exact charge-neutral subspace invariant:
    // For any signal node a in {1, 2, 3}, (e_a - e_0)^T E_{C,closed}^+ (e_a - e_0) == (E_{C,dropped})[a-1][a-1]
    let ec_pinv = &raw_check["ec_reduced_Hz"];
    let ec_drop = &dropped["ec_reduced_Hz"];
    for a in 1..4 {
        let dipole_ec =
            ec_pinv[a][a].as_f64().unwrap() + ec_pinv[0][0].as_f64().unwrap() - 2.0 * ec_pinv[a][0].as_f64().unwrap();
        let drop_ec = ec_drop[a - 1][a - 1].as_f64().unwrap();
        assert!(rel_err(dipole_ec, drop_ec) < 1.0e-12, "a={a}: dipole={dipole_ec}, drop={drop_ec}");
    }

    // Floating Schur elimination on a closed ungrounded network preserves the constant-potential null mode on S
    let closed_schur =
        capacitance_matrix_with_params(&res_closed, &json!({ "ground": "c0", "eliminate": "floating_schur" })).unwrap();
    assert_eq!(closed_schur["has_null_mode"], json!(true));
    assert_eq!(closed_schur["null_modes"], json!(1));
    for i in 0..3 {
        let row_sum: f64 = (0..3).map(|j| closed_schur["maxwell_reduced_F"][i][j].as_f64().unwrap()).sum();
        assert!(row_sum.abs() < 1.0e-27, "closed Schur row {i} sum = {row_sum:e}");
    }

    // 2) Grounded-box network where `c0` is an ungrounded floating pocket/shield (Q_c0 = 0, box has stray cap)
    let mut c_box = c_closed.clone();
    c_box[0][0] += 30.0 * f;
    c_box[1][1] += 5.0 * f;
    c_box[2][2] += 8.0 * f;
    c_box[3][3] += 8.0 * f;
    let res_box = make_result(&["c0", "c1", "c2", "c3"], c_box.clone());

    let box_schur =
        capacitance_matrix_with_params(&res_box, &json!({ "ground": "c0", "eliminate": "floating_schur" })).unwrap();
    assert_eq!(box_schur["has_null_mode"], json!(false));
    assert_eq!(box_schur["null_modes"], json!(0));

    // Analytical C_ss - C_sg * C_gg^-1 * C_gs check (1x1 C_gg = c_box[0][0])
    let c00 = c_box[0][0];
    for a in 0..3 {
        for b in 0..3 {
            let want = c_box[a + 1][b + 1] - (c_box[a + 1][0] * c_box[0][b + 1]) / c00;
            let got = box_schur["maxwell_reduced_F"][a][b].as_f64().unwrap();
            assert!(rel_err(got, want) < 1.0e-14, "({a},{b}): got={got:e}, want={want:e}");
        }
    }

    // Verify C_red * ec_reduced_Hz == (e^2 / (2h)) * I_3
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    for i in 0..3 {
        for j in 0..3 {
            let mut prod = 0.0;
            for k in 0..3 {
                prod += box_schur["maxwell_reduced_F"][i][k].as_f64().unwrap()
                    * box_schur["ec_reduced_Hz"][k][j].as_f64().unwrap();
            }
            let target = if i == j { scale } else { 0.0 };
            assert!((prod - target).abs() < 1.0e-13 * scale, "({i},{j}): prod={prod:e}, target={target:e}");
        }
    }
}

#[test]
fn invalid_parameters_return_clean_diagnostic_errors() {
    let f = 1.0e-15;
    let c = vec![
        vec![100.0 * f, -30.0 * f, -20.0 * f],
        vec![-30.0 * f, 80.0 * f, -10.0 * f],
        vec![-20.0 * f, -10.0 * f, 70.0 * f],
    ];
    let res = make_result(&["c0", "c1", "c2"], c);

    // Unknown top-level parameter
    assert!(capacitance_matrix_with_params(&res, &json!({ "bogus": 1 })).is_err());
    // Non-existent ground terminal
    assert!(capacitance_matrix_with_params(&res, &json!({ "ground": "nope" })).is_err());
    // Eliminating all terminals
    assert!(capacitance_matrix_with_params(&res, &json!({ "ground": ["c0", "c1", "c2"] })).is_err());
    // Qubit pad also marked as ground
    assert!(capacitance_matrix_with_params(
        &res,
        &json!({ "ground": "c0", "qubits": [{ "name": "Q1", "pads": ["c0", "c1"] }] })
    )
    .is_err());
    // Lumped resonator missing C_r_F
    assert!(capacitance_matrix_with_params(
        &res,
        &json!({ "resonators": [{ "name": "R1", "terminal": "c1", "f_r_Hz": 6.5e9, "resonator_type": "lumped" }] })
    )
    .is_err());
    // Resonator sharing a qubit pad
    assert!(capacitance_matrix_with_params(
        &res,
        &json!({
            "ground": "c0",
            "qubits": [{ "name": "Q1", "pads": ["c2"] }],
            "resonators": [{ "name": "R1", "terminal": "c2", "f_r_Hz": 6.5e9 }]
        })
    )
    .is_err());
}
