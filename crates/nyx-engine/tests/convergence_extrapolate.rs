// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Two-grid Richardson extrapolation (`richardson_log_spd` and `richardson_linear`) and bidirectional `factor` (`< 1` refining, `> 1`
//! coarsening) in `convergence_check` (QFR-3 part 1).

use femcore::linalg::driver::dense_eigh;
use nyx_engine::convergence::{
    apply_extrapolated_capacitance, extrapolate_capacitance, extrapolated_capacitance, should_apply_extrapolated,
    Limits,
};
use nyx_engine::{Engine, FemcoreMesher};
use nyx_physics::{
    BcDesc, CapacitanceMatrixData, Ctx, Dims, Error, FieldDesc, FieldLink, Mesh, ParamDesc, Physics, ProblemSpec,
    PropDesc, Registry, Scalar, Sizing, Solution, SolveResult, StudyDesc,
};
use serde_json::{json, Value};
use std::sync::Arc;

const E_CHARGE: f64 = 1.602_176_634e-19;
const PLANCK_H: f64 = 6.626_070_15e-34;

#[allow(clippy::needless_range_loop)]
fn build_cap_data(terminals: &[&str], maxwell: Vec<Vec<f64>>) -> CapacitanceMatrixData {
    let n = maxwell.len();
    let mut mutual = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut s = 0.0;
        for k in 0..n {
            s += maxwell[i][k];
            if i != k {
                mutual[i][k] = -maxwell[i][k];
            }
        }
        mutual[i][i] = s;
    }
    let mut flat = Vec::with_capacity(n * n);
    for r in &maxwell {
        flat.extend_from_slice(r);
    }
    let (eigs, vecs) = dense_eigh(&flat, n).unwrap();
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; n]; n];
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                ec[i][j] += (scale / eigs[k]) * vecs[k][i] * vecs[k][j];
            }
        }
    }
    CapacitanceMatrixData {
        terminals: terminals.iter().map(|s| (*s).to_string()).collect(),
        maxwell_matrix_f: maxwell,
        mutual_matrix_f: mutual,
        ec_matrix_hz: ec,
        ..Default::default()
    }
}

/// Synthetic 3D physics whose effective mesh spacing `h_eff = cbrt(6.0 / n_cells)` scales linearly with `mesh.h` on a uniform
/// 3D lattice box (`cells = round(1/h)^3`, `n_cells = 6 * round(1/h)^3` -> `h_eff == h` when `1/h` is an exact integer!).
struct SyntheticCapPhysics {
    exponential_decay: bool,
}

impl SyntheticCapPhysics {
    fn exact_asymptote() -> CapacitanceMatrixData {
        // 3-terminal Maxwell matrix: self_inf = [45 fF, 50 fF, 30 fF], mutuals m01 = 22 fF, m02 = 3 fF, m12 = 14 fF
        let s = [45.0e-15, 50.0e-15, 30.0e-15];
        let (m01, m02, m12) = (22.0e-15, 3.0e-15, 14.0e-15);
        let c = vec![
            vec![s[0] + m01 + m02, -m01, -m02],
            vec![-m01, s[1] + m01 + m12, -m12],
            vec![-m02, -m12, s[2] + m02 + m12],
        ];
        build_cap_data(&["pad1", "pad2", "coupler"], c)
    }

    fn at_h(&self, h: f64) -> CapacitanceMatrixData {
        let s_inf = [45.0e-15, 50.0e-15, 30.0e-15];
        let (m01_inf, m02_inf, m12_inf) = (22.0e-15, 3.0e-15, 14.0e-15);
        let c = if self.exponential_decay {
            // Exact O(h) in log-space: m(h) = m_inf * exp(alpha * h), self_inf(h) = s_inf * exp(gamma * h)
            let s0 = s_inf[0] * (-0.35 * h).exp();
            let s1 = s_inf[1] * (-0.28 * h).exp();
            let s2 = s_inf[2] * (-0.42 * h).exp();
            let m01 = m01_inf * (-0.50 * h).exp();
            let m02 = m02_inf * (-0.65 * h).exp();
            let m12 = m12_inf * (-0.30 * h).exp();
            vec![vec![s0 + m01 + m02, -m01, -m02], vec![-m01, s1 + m01 + m12, -m12], vec![-m02, -m12, s2 + m02 + m12]]
        } else {
            // Exact additive O(h): C(h) = C_inf + h * DeltaC
            let exact = Self::exact_asymptote().maxwell_matrix_f;
            let dc = [[-12.0e-15, 4.0e-15, 1.0e-15], [4.0e-15, -10.0e-15, 3.0e-15], [1.0e-15, 3.0e-15, -8.0e-15]];
            (0..3).map(|i| (0..3).map(|j| exact[i][j] + h * dc[i][j]).collect()).collect()
        };
        build_cap_data(&["pad1", "pad2", "coupler"], c)
    }
}

impl Physics for SyntheticCapPhysics {
    fn id(&self) -> &'static str {
        "dfs"
    }
    fn name(&self) -> &'static str {
        "Synthetic O(h) capacitance"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc {
            kind: "electrostatic".into(),
            name: "Electrostatic".into(),
            doc: "".into(),
            dims: Dims::Both,
            params: vec![],
        }]
    }
    fn bc_kinds(&self) -> Vec<BcDesc> {
        vec![]
    }
    fn properties(&self) -> Vec<PropDesc> {
        vec![]
    }
    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }
    fn fields(&self, _: &str) -> Vec<FieldDesc> {
        vec![]
    }
    fn exports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn validate(&self, _: &ProblemSpec, _: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, spec: &ProblemSpec, _: Option<&Arc<femcore::material::Library>>) -> Sizing {
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        // Estimate exact characteristic spacing h from bounding box / cell count: for our 2D unit square
        // `[0,1]^2` meshed with `FemcoreMesher`, instead of relying on cell count jitter, read `h` from `1 / sqrt(cells)`?
        // Better: on a 1D-extruded 3D mesh (`sweep: [{ thickness: 1.0, n: N }]`) `l.n` scales as `(N / factor)`!
        // With `n = 12` at `h_coarse` and `factor = 1/1.5` (`n = 18`) or `n = 18` at `h_fine` and `factor = 1.5` (`n = 12`),
        // `1.0 / n_z` has exact ratio `18 / 12 == 1.5`!
        let d = nyx_physics::discretise(mesh)?;
        let mut y_set: Vec<i64> = (0..d.n_points()).map(|i| (d.point(i)[1] * 1e9).round() as i64).collect();
        y_set.sort_unstable();
        y_set.dedup();
        let n_layers = (y_set.len().saturating_sub(1)).max(1) as f64;
        let h_eff = 1.0 / n_layers;
        let cap = self.at_h(h_eff);
        let mut scalars = Vec::new();
        let mut solutions = Vec::new();
        for i in 0..3 {
            solutions.push(Solution {
                label: format!("terminal {}", cap.terminals[i]),
                value: Some(cap.maxwell_matrix_f[i][i]),
                ..Default::default()
            });
            for j in 0..3 {
                scalars.push((format!("C_{i}_{j} (F)"), Scalar::Num(cap.maxwell_matrix_f[i][j])));
                scalars.push((format!("E_C_{i}_{j} (Hz)"), Scalar::Num(cap.ec_matrix_hz[i][j])));
            }
        }
        scalars.push(cap.scalar());
        Ok(SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            order: spec.order,
            dim: mesh.dim() as u8,
            solutions,
            scalars,
            log: vec![],
            seconds: 0.001,
            dofs: mesh.n_points(),
        })
    }
}

fn assert_maxwell_spd(cap: &CapacitanceMatrixData) {
    let n = cap.terminals.len();
    assert_eq!(cap.maxwell_matrix_f.len(), n);
    let mut flat = Vec::with_capacity(n * n);
    for i in 0..n {
        assert!(cap.maxwell_matrix_f[i][i] > 0.0, "diagonal C[{i},{i}] must be positive");
        let mut sum_off = 0.0;
        for j in 0..n {
            let cij = cap.maxwell_matrix_f[i][j];
            let cji = cap.maxwell_matrix_f[j][i];
            assert!((cij - cji).abs() <= 1e-28, "symmetry violated at ({i},{j}): {cij} vs {cji}");
            if i != j {
                assert!(cij <= 0.0, "off-diagonal Maxwell C[{i},{j}] must be <= 0, got {cij}");
                sum_off += cij.abs();
            }
            flat.push(cij);
        }
        assert!(
            cap.maxwell_matrix_f[i][i] >= sum_off - 1e-28,
            "diagonal dominance violated on row {i}: {} < {sum_off}",
            cap.maxwell_matrix_f[i][i]
        );
    }
    let (eigs, _) = dense_eigh(&flat, n).expect("dense_eigh succeeds");
    for (k, &lam) in eigs.iter().enumerate() {
        assert!(lam > 0.0, "Maxwell capacitance eigenvalue {k} must be strictly positive, got {lam}");
    }
}

fn max_rel_matrix_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let mut m = 0.0_f64;
    for i in 0..a.len() {
        for j in 0..a[i].len() {
            let scale = a[i][j].abs().max(b[i][j].abs()).max(1e-30);
            m = m.max((a[i][j] - b[i][j]).abs() / scale);
        }
    }
    m
}

fn setup_swept_synthetic_engine(exponential_decay: bool, n_sweep: usize) -> Engine {
    let mut reg = Registry::new();
    reg.register(Arc::new(SyntheticCapPhysics { exponential_decay })).unwrap();
    let mut e = Engine::with(reg, Arc::new(FemcoreMesher));
    e.call("project_new", &json!({ "title": "synth_3d_cap", "dim": 2, "h": 0.5 })).unwrap();
    e.call("add_rect", &json!({ "name": "r", "min": [0.0, 0.0], "max": [1.0, 1.0] })).unwrap();
    e.call("add_overlay", &json!({ "name": "ov", "regions": [{ "name": "sub", "input": "r" }] })).unwrap();
    e.call("set_mesh", &json!({ "sweep": [{ "thickness": 1.0, "n": n_sweep }] })).unwrap();
    e.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
    e
}

#[test]
fn synthetic_oh_capacitance_recovers_asymptote_for_refining_and_coarsening_factors_and_preserves_spd() {
    let exact = SyntheticCapPhysics::exact_asymptote();

    // 1. `richardson_log_spd` with `factor = 1.0 / 1.5` (refining: 12 -> 18 sweep layers, r = 1.5, default order = 1.0, w1 = 3.0)
    let mut e_ref = setup_swept_synthetic_engine(true, 12);
    let ans_ref = e_ref
        .call(
            "convergence_check",
            &json!({
                "physics": "dfs",
                "factor": 1.0 / 1.5,
                "tolerance": 0.25,
                "extrapolate": "richardson_log_spd",
                "apply_extrapolated": true
            }),
        )
        .unwrap();
    let ext_ref = &ans_ref["extrapolated"];
    assert_eq!(ext_ref["method"], json!("richardson_log_spd"));
    assert_eq!(ext_ref["order"], json!(1.0), "3D capacitance mesh defaults to p = 1.0");
    assert!((ext_ref["r"].as_f64().unwrap() - 1.5).abs() < 1e-14);
    assert!((ext_ref["weight_fine"].as_f64().unwrap() - 3.0).abs() < 1e-14);
    assert!(should_apply_extrapolated(&ans_ref));
    let cap_ref = extrapolated_capacitance(&ans_ref).expect("extrapolated capacitance decoded");
    assert_maxwell_spd(&cap_ref);
    assert!(
        max_rel_matrix_diff(&cap_ref.maxwell_matrix_f, &exact.maxwell_matrix_f) < 1e-12,
        "refining log-SPD max_rel = {:e}",
        max_rel_matrix_diff(&cap_ref.maxwell_matrix_f, &exact.maxwell_matrix_f)
    );
    assert!(max_rel_matrix_diff(&cap_ref.mutual_matrix_f, &exact.mutual_matrix_f) < 1e-12);
    assert!(max_rel_matrix_diff(&cap_ref.ec_matrix_hz, &exact.ec_matrix_hz) < 1e-12);
    assert!(ext_ref["max_entry_rel_error"].as_f64().unwrap() > 1e-3);

    // 2. `richardson_log_spd` with `factor = 1.5` (coarsening: present = 18 layers fine, comparison = 12 layers coarse, r = 1.5, w1 = 3.0)
    let e_coarse = setup_swept_synthetic_engine(true, 18);
    let plan_coarse = e_coarse.convergence_job().plan("dfs", 1.5).unwrap();
    assert!(
        plan_coarse["cells"]["coarse"].as_u64().unwrap() < plan_coarse["cells"]["refined"].as_u64().unwrap(),
        "plan with factor > 1 orders coarse < refined: {plan_coarse}"
    );
    assert!(plan_coarse["refined_over_present"].as_f64().unwrap() < 1.0);

    let ans_coarse = e_coarse
        .convergence_job()
        .with_limits(Limits {
            extrapolate: Some("richardson_log_spd".into()),
            apply_extrapolated: false,
            ..Limits::default()
        })
        .run("dfs", 1.5, Some(0.25))
        .unwrap();
    assert!(
        ans_coarse["cells"]["coarse"].as_u64().unwrap() < ans_coarse["cells"]["refined"].as_u64().unwrap(),
        "{ans_coarse}"
    );
    assert!((ans_coarse["extrapolated"]["r"].as_f64().unwrap() - 1.5).abs() < 1e-14);
    assert!((ans_coarse["extrapolated"]["weight_fine"].as_f64().unwrap() - 3.0).abs() < 1e-14);
    assert!(!should_apply_extrapolated(&ans_coarse));
    let cap_coarse = extrapolated_capacitance(&ans_coarse).expect("extrapolated capacitance decoded");
    assert_maxwell_spd(&cap_coarse);
    assert!(
        max_rel_matrix_diff(&cap_coarse.maxwell_matrix_f, &exact.maxwell_matrix_f) < 1e-12,
        "coarsening log-SPD max_rel = {:e}",
        max_rel_matrix_diff(&cap_coarse.maxwell_matrix_f, &exact.maxwell_matrix_f)
    );
    assert!(max_rel_matrix_diff(&cap_coarse.maxwell_matrix_f, &cap_ref.maxwell_matrix_f) < 1e-14);

    // 3. `richardson_linear` on additive O(h) matrices for BOTH `factor = 1.0 / 1.5` and `factor = 1.5`
    let mut e_lin_ref = setup_swept_synthetic_engine(false, 12);
    let ans_lin_ref = e_lin_ref
        .call(
            "convergence_check",
            &json!({ "physics": "dfs", "factor": 1.0 / 1.5, "extrapolate": "richardson_linear" }),
        )
        .unwrap();
    let cap_lin_ref = extrapolated_capacitance(&ans_lin_ref).unwrap();
    assert_maxwell_spd(&cap_lin_ref);
    assert!(max_rel_matrix_diff(&cap_lin_ref.maxwell_matrix_f, &exact.maxwell_matrix_f) < 1e-12);

    let e_lin_coarse = setup_swept_synthetic_engine(false, 18);
    let ans_lin_coarse = e_lin_coarse
        .convergence_job()
        .with_limits(Limits { extrapolate: Some("richardson_linear".into()), ..Limits::default() })
        .run("dfs", 1.5, None)
        .unwrap();
    let cap_lin_coarse = extrapolated_capacitance(&ans_lin_coarse).unwrap();
    assert_maxwell_spd(&cap_lin_coarse);
    assert!(max_rel_matrix_diff(&cap_lin_coarse.maxwell_matrix_f, &exact.maxwell_matrix_f) < 1e-12);

    // 4. Severe coarse coupling overshoot where linear Richardson flips sign (`C_01 > 0` non-physical!) while `richardson_log_spd` preserves SPD & signs
    let c_overshoot_coarse = build_cap_data(&["a", "b"], vec![vec![50.0e-15, -4.0e-15], vec![-4.0e-15, 50.0e-15]]);
    let c_overshoot_fine = build_cap_data(&["a", "b"], vec![vec![48.0e-15, -1.0e-15], vec![-1.0e-15, 48.0e-15]]);
    let (lin_bad, _, _) =
        extrapolate_capacitance(&c_overshoot_coarse, &c_overshoot_fine, "richardson_linear", 1.5, 1.0).unwrap();
    assert!(lin_bad.maxwell_matrix_f[0][1] > 0.0, "linear overshoot flipped off-diagonal sign");
    let (log_good, _, _) =
        extrapolate_capacitance(&c_overshoot_coarse, &c_overshoot_fine, "richardson_log_spd", 1.5, 1.0).unwrap();
    assert_maxwell_spd(&log_good);
    assert!(log_good.maxwell_matrix_f[0][1] < 0.0, "log-SPD strictly preserves negative off-diagonal Maxwell entry");

    // 5. `apply_extrapolated_capacitance` updates `CapacitanceMatrixData`, `C_*` / `E_C_*` scalars, and solution values
    e_ref.call("solve_run", &json!({})).unwrap();
    let mut promoted: SolveResult = e_ref.result("dfs").unwrap().0.as_ref().clone();
    apply_extrapolated_capacitance(&mut promoted, &cap_ref);
    assert_eq!(CapacitanceMatrixData::of(&promoted).unwrap(), cap_ref);
    assert_eq!(promoted.solutions[0].value, Some(cap_ref.maxwell_matrix_f[0][0]));
    let c01_scalar = promoted
        .scalars
        .iter()
        .find_map(|(k, v)| match v {
            Scalar::Num(x) if k == "C_0_1 (F)" => Some(*x),
            _ => None,
        })
        .unwrap();
    assert_eq!(c01_scalar, cap_ref.maxwell_matrix_f[0][1]);
}

fn build_two_pad_sheet_engine(h_global: f64, h_edge: f64) -> Engine {
    let mut e = Engine::new();
    let mut call = |cmd: &str, args: Value| e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}: {er}"));
    call("project_new", json!({ "title": "two_pad_sheet_qcap", "dim": 3, "h": h_global }));
    call("add_box", json!({ "name": "sub", "min": [-120e-6, -120e-6, -60e-6], "max": [120e-6, 120e-6, 0.0] }));
    call("add_box", json!({ "name": "vac", "min": [-120e-6, -120e-6, 0.0], "max": [120e-6, 120e-6, 60e-6] }));
    call(
        "add_polygons",
        json!({
            "name": "p1_poly",
            "polygons": [[[-55e-6, -40e-6], [-10e-6, -40e-6], [-10e-6, 40e-6], [-55e-6, 40e-6]]]
        }),
    );
    call(
        "add_polygons",
        json!({
            "name": "p2_poly",
            "polygons": [[[10e-6, -40e-6], [55e-6, -40e-6], [55e-6, 40e-6], [10e-6, 40e-6]]]
        }),
    );
    call("add_sheet", json!({ "name": "pad1", "input": "p1_poly", "z": 0.0 }));
    call("add_sheet", json!({ "name": "pad2", "input": "p2_poly", "z": 0.0 }));
    call(
        "add_overlay",
        json!({
            "name": "ov",
            "regions": [
                { "name": "substrate", "input": "sub" },
                { "name": "vacuum", "input": "vac" }
            ]
        }),
    );
    call("add_size", json!({ "name": "sheet_rim", "sel": "pad1:edge | pad2:edge", "h": h_edge }));
    call("add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    call("set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }));
    call("set_material", json!({ "physics": "dfs", "region": "vacuum", "material": { "eps_r": 1.0 } }));
    call("add_condition", json!({ "physics": "dfs", "on": "pad1", "kind": "terminal", "params": { "name": "pad1" } }));
    call("add_condition", json!({ "physics": "dfs", "on": "pad2", "kind": "terminal", "params": { "name": "pad2" } }));
    call("add_condition", json!({ "physics": "dfs", "on": "outer", "kind": "ground", "params": {} }));
    e
}

#[test]
fn real_dfs_sheet_electrostatic_extrapolation_reduces_error_against_finer_reference() {
    // Reference fine solve (`h = 12 um`, `h_edge = 3.5 um`)
    let mut e_ref = build_two_pad_sheet_engine(12.0e-6, 3.5e-6);
    e_ref.call("solve_run", &json!({})).unwrap();
    let ref_cap = CapacitanceMatrixData::of(e_ref.result("dfs").unwrap().0).expect("reference capacitance");

    // Two-grid check starting from coarse (`h = 36 um`, `h_edge = 10.5 um`) with `factor = 2/3` (`h_fine = 24 um`, `h_edge = 7 um`),
    // `scale_sizes = true` and `extrapolate = "richardson_log_spd"`
    let mut e = build_two_pad_sheet_engine(36.0e-6, 10.5e-6);
    e.call("solve_run", &json!({})).unwrap();
    let coarse_cap = CapacitanceMatrixData::of(e.result("dfs").unwrap().0).unwrap();

    // Also solve the un-extrapolated fine grid (`24 um`, `7 um`) directly to measure error reduction
    let mut e_mid = build_two_pad_sheet_engine(24.0e-6, 7.0e-6);
    e_mid.call("solve_run", &json!({})).unwrap();
    let mid_cap = CapacitanceMatrixData::of(e_mid.result("dfs").unwrap().0).unwrap();

    let conv = e
        .call(
            "convergence_check",
            &json!({
                "physics": "dfs",
                "factor": 2.0 / 3.0,
                "scale_sizes": true,
                "tolerance": 0.15,
                "extrapolate": "richardson_log_spd"
            }),
        )
        .unwrap();
    let ext_cap = extrapolated_capacitance(&conv).expect("extrapolated capacitance");
    assert_maxwell_spd(&ext_cap);
    assert_eq!(conv["extrapolated"]["order"], json!(1.0));
    assert!(ext_cap.ec_matrix_hz[0][0].is_finite() && ext_cap.ec_matrix_hz[0][0] > 1.0e6);

    let err_coarse = max_rel_matrix_diff(&coarse_cap.maxwell_matrix_f, &ref_cap.maxwell_matrix_f);
    let err_fine = max_rel_matrix_diff(&mid_cap.maxwell_matrix_f, &ref_cap.maxwell_matrix_f);
    let err_ext = max_rel_matrix_diff(&ext_cap.maxwell_matrix_f, &ref_cap.maxwell_matrix_f);
    assert!(
        err_ext < err_fine && err_fine < err_coarse,
        "Richardson log-SPD must reduce capacitance error vs reference: coarse={err_coarse:.4e}, fine={err_fine:.4e}, extrap={err_ext:.4e}"
    );
}

#[test]
fn validation_rejects_invalid_extrapolate_options_and_bad_factors() {
    let mut e = setup_swept_synthetic_engine(true, 8);
    for bad in [
        json!({ "extrapolate": "unknown_method" }),
        json!({ "extrapolate": 42 }),
        json!({ "extrapolate": "richardson_log_spd", "extrapolate_order": 0.0 }),
        json!({ "extrapolate": "richardson_log_spd", "extrapolate_order": -1.5 }),
        json!({ "extrapolate": "richardson_log_spd", "extrapolate_order": "two" }),
        json!({ "extrapolate": "richardson_log_spd", "apply_extrapolated": "yes" }),
    ] {
        let mut a = json!({ "physics": "dfs", "factor": 0.6 });
        for (k, v) in bad.as_object().unwrap() {
            a[k] = v.clone();
        }
        assert!(e.call("convergence_check", &a).is_err(), "expected validation failure for {a}");
    }

    // Direct `ConvergenceJob` factor validation: 0.0, negative, 1.0, NaN, INFINITY all fail
    for bad_factor in [0.0, -0.6, 1.0, f64::NAN, f64::INFINITY] {
        assert!(e.convergence_job().run("dfs", bad_factor, None).is_err(), "bad factor {bad_factor}");
        assert!(e.convergence_job().plan("dfs", bad_factor).is_err(), "bad plan factor {bad_factor}");
    }
    // And non-finite `extrapolate_order` constructed in Rust is rejected too
    assert!(e
        .convergence_job()
        .with_limits(Limits {
            extrapolate: Some("richardson_linear".into()),
            extrapolate_order: Some(f64::NAN),
            ..Limits::default()
        })
        .run("dfs", 0.6, None)
        .is_err());
}
