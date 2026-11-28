// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification suite for FR-38 (`nyx-dfs`, ADR 139):
//! 1. 1D radiative-conductive slab (`BcKind::Radiation` $T^4$) vs semi-analytical ODE (`< 1e-6` in `<= 6` Newton steps, `symbolic_analyses == 1`),
//!    preserving the linear SPD fast path bit-identically when `Radiation` and $k(T)$ are absent.
//! 2. Kirchhoff-transform analytical verification for nonlinear $k(T)$ (`Tabulated`) and transient BDF2 radiative cooling (`heat_transient` + `Cp(T)`).
//! 3. `Study::Conduction` 2D & 3D analytical Joule heating, `CurrentDensity` BC, $\sigma(T)$ harmonic resistance, and exact power conservation (`< 1e-10`).

use femcore::fem::assemble_mass;
use femcore::material::state_t::{TempCurve, TempTensor};
use femcore::testkit::{block, rect};
use nyx_dfs::{run, solve_conduction, solve_heat, Bc, BcKind, Problem, Study, SIGMA_SB};
use nyx_physics::{Mesh, Scalar};

fn num_scalar(res: &nyx_physics::SolveResult, key: &str) -> f64 {
    res.scalars
        .iter()
        .find_map(|(k, v)| match (k == key, v) {
            (true, Scalar::Num(x)) => Some(*x),
            _ => None,
        })
        .unwrap_or_else(|| panic!("missing numeric scalar '{key}' in {:?}", res.scalars))
}

/// Solves `eps * F * sigma_SB * (T_L^4 - T_amb^4) + (k / L) * (T_L - T_0) - 0.5 * Q * L = 0` to machine precision.
fn solve_1d_radiation_slab_tip_temp(
    k: f64,
    length: f64,
    q_vol: f64,
    t0: f64,
    emissivity: f64,
    view_factor: f64,
    t_amb: f64,
) -> f64 {
    let rad = emissivity * view_factor * SIGMA_SB;
    let mut tl = t0;
    for _ in 0..30 {
        let g = rad * (tl.powi(4) - t_amb.powi(4)) + (k / length) * (tl - t0) - 0.5 * q_vol * length;
        let dg = 4.0 * rad * tl.powi(3) + k / length;
        let step = g / dg;
        tl -= step;
        if step.abs() <= 1e-15 * tl.abs() {
            break;
        }
    }
    tl
}

#[test]
fn radiation_bc_1d_slab_matches_ode_within_1e_6_in_at_most_6_newton_steps() {
    let length = 0.01;
    let height = 0.002;
    let mesh = rect(length, height, 20, 2);
    let mut lbl_x0 = None;
    let mut lbl_xl = None;
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        if p0[0].abs() < 1e-12 && p1[0].abs() < 1e-12 {
            lbl_x0 = Some(mesh.facet_labels[f]);
        } else if (p0[0] - length).abs() < 1e-12 && (p1[0] - length).abs() < 1e-12 {
            lbl_xl = Some(mesh.facet_labels[f]);
        }
    }
    let (left_lbl, right_lbl) = (lbl_x0.unwrap(), lbl_xl.unwrap());

    let k_cond = 2.0;
    let q_vol = 5.0e5;
    let t0 = 800.0;
    let t_amb = 300.0;
    let emissivity = 0.85;
    let view_factor = 0.9;

    let mut spec = Problem::new(Study::HeatStationary, 2).with_tol(1e-12);
    spec.materials.default.insert("k".into(), k_cond);
    spec.materials.default.insert("Q".into(), q_vol);
    spec.bcs = vec![
        Bc { labels: vec![left_lbl], kind: BcKind::Temperature(t0) },
        Bc { labels: vec![right_lbl], kind: BcKind::Radiation { emissivity, t_ambient: t_amb, view_factor } },
    ];

    let res = solve_heat(&mesh, &spec).expect("nonlinear radiative slab solve");
    let diag = res.diagnostics.as_ref().expect("diagnostics present");
    let newton_iters = diag.iterations.expect("Newton iterations recorded");
    assert!(
        newton_iters <= 6,
        "expected <= 6 Newton iterations with exact +4*eps*F*sigma*T^3 tangent, got {newton_iters}"
    );
    assert!(
        res.log.iter().any(|l| l.contains("1 symbolic")),
        "expected single symbolic factorization reuse in log: {:?}",
        res.log
    );

    let tl_exact = solve_1d_radiation_slab_tip_temp(k_cond, length, q_vol, t0, emissivity, view_factor, t_amb);
    let slope0 = (tl_exact - t0) / length + 0.5 * q_vol * length / k_cond;
    let t_nodal = &res.solutions[0].field("T").expect("T field").values;

    let mut max_rel_err = 0.0_f64;
    for (v, &tv) in t_nodal.iter().enumerate().take(mesh.n_points()) {
        let x = mesh.point(v)[0];
        let t_exact = t0 + slope0 * x - 0.5 * (q_vol / k_cond) * x * x;
        let rel = (tv - t_exact).abs() / t0;
        max_rel_err = max_rel_err.max(rel);
    }
    assert!(
        max_rel_err < 1e-6,
        "1D radiative-conductive slab relative error {max_rel_err:.3e} >= 1e-6 (T_L exact = {tl_exact:.6} K)"
    );

    // Verify radiative heat flux balance at x = L: k * (-dT/dx)|_{x=L} == eps * F * sigma_SB * (T_L^4 - T_amb^4)
    let qx = res.solutions[0].flat("qx").expect("qx field");
    let mut qx_right_sum = 0.0;
    let mut qx_right_cnt = 0usize;
    for (v, &qxv) in qx.iter().enumerate().take(mesh.n_points()) {
        if (mesh.point(v)[0] - length).abs() < 1e-12 {
            qx_right_sum += qxv;
            qx_right_cnt += 1;
        }
    }
    let qx_right = qx_right_sum / (qx_right_cnt as f64);
    let q_rad_exact = emissivity * view_factor * SIGMA_SB * (tl_exact.powi(4) - t_amb.powi(4));
    let flux_rel_err = (qx_right - q_rad_exact).abs() / q_rad_exact;
    assert!(flux_rel_err < 1e-6, "radiative boundary flux error {flux_rel_err:.3e}: {qx_right} vs {q_rad_exact}");

    // Verify linear fast path stays bit-identical when no Radiation or k(T) curve is present.
    let mut lin_spec = Problem::new(Study::HeatStationary, 2);
    lin_spec.materials.default.insert("k".into(), k_cond);
    lin_spec.materials.default.insert("Q".into(), q_vol);
    lin_spec.bcs = vec![
        Bc { labels: vec![left_lbl], kind: BcKind::Temperature(t0) },
        Bc { labels: vec![right_lbl], kind: BcKind::Convection { h: 25.0, text: t_amb } },
    ];
    assert!(!lin_spec.is_nonlinear_heat());
    let r_lin1 = solve_heat(&mesh, &lin_spec).unwrap();
    let r_lin2 = solve_heat(&mesh, &lin_spec).unwrap();
    assert_eq!(r_lin1.solutions[0].fields[0].values, r_lin2.solutions[0].fields[0].values);
    assert!(!r_lin1.log.iter().any(|l| l.contains("nonlinear stationary heat")));
}

#[test]
fn nonlinear_kirchhoff_conductivity_and_transient_radiation_bdf2() {
    // Part A: Steady nonlinear k(T) = k0 * (1 + beta * (T - T1)) via Tabulated TempCurve vs exact Kirchhoff transform.
    let length = 0.1;
    let mesh = rect(length, 0.01, 24, 2);
    let mut left_lbl = 1;
    let mut right_lbl = 2;
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        if p0[0].abs() < 1e-12 && p1[0].abs() < 1e-12 {
            left_lbl = mesh.facet_labels[f];
        } else if (p0[0] - length).abs() < 1e-12 && (p1[0] - length).abs() < 1e-12 {
            right_lbl = mesh.facet_labels[f];
        }
    }

    let (t1, t2) = (300.0, 600.0);
    let k0 = 10.0;
    let beta = 0.004; // k(300) = 10, k(600) = 22 W/(m*K)
    let curve_k = TempCurve::Tabulated {
        t: vec![200.0, 800.0],
        values: vec![k0 * (1.0 + beta * (200.0 - t1)), k0 * (1.0 + beta * (800.0 - t1))],
        t_ref: Some(300.0),
        source: Vec::new(),
    };
    let dom = mesh.cell_labels[0];
    let mut spec_nl =
        Problem::new(Study::HeatStationary, 2).with_k_curve(dom, TempTensor::Scalar(curve_k), [0.0; 3]).with_tol(1e-12);
    spec_nl.materials.default.insert("k".into(), k0);
    spec_nl.bcs = vec![
        Bc { labels: vec![left_lbl], kind: BcKind::Temperature(t1) },
        Bc { labels: vec![right_lbl], kind: BcKind::Temperature(t2) },
    ];

    let res_nl = solve_heat(&mesh, &spec_nl).expect("nonlinear k(T) Kirchhoff solve");
    let iters_nl = res_nl.diagnostics.as_ref().unwrap().iterations.unwrap();
    assert!(iters_nl <= 6, "Kirchhoff nonlinear k(T) should converge in <= 6 Newton steps, got {iters_nl}");

    let theta_2 = k0 * ((t2 - t1) + 0.5 * beta * (t2 - t1).powi(2));
    let t_vals = &res_nl.solutions[0].field("T").unwrap().values;
    let mut max_kirchhoff_err = 0.0_f64;
    for (v, &tv) in t_vals.iter().enumerate().take(mesh.n_points()) {
        let x = mesh.point(v)[0];
        let theta_x = theta_2 * (x / length);
        let t_exact = t1 + (-1.0 + (1.0 + 2.0 * beta * theta_x / k0).sqrt()) / beta;
        max_kirchhoff_err = max_kirchhoff_err.max((tv - t_exact).abs() / t2);
    }
    assert!(max_kirchhoff_err < 1e-5, "Kirchhoff nonlinear k(T) error {max_kirchhoff_err:.3e} >= 1e-5");

    // Part B: Transient BDF2 radiative cooling of a high-conductivity 3D body in vacuum (pure RadiationBC, 0 Dirichlet BCs!).
    let side = 1.0e-3; // 1 mm cube: V = 1e-9 m^3, A = 6e-6 m^2
    let mesh3 = block(side, side, side, 2, 2, 2);
    let all_outer: Vec<i32> = {
        let mut s: Vec<i32> = mesh3.facet_labels.clone();
        s.sort_unstable();
        s.dedup();
        s
    };
    let rho = 2000.0;
    let cp0 = 500.0; // rho * Cp * V = 1e-3 J/K
    let k_high = 500.0; // Bi_rad = eps * sigma * T0^3 * L / k ~ 2.4e-5 << 1 (lumped limit)
    let (t_init, t_amb) = (600.0, 300.0);
    let eps = 1.0;
    let t_end = 40.0;
    let n_out = 40;

    let cp_curve =
        TempCurve::Tabulated { t: vec![100.0, 1000.0], values: vec![cp0, cp0], t_ref: Some(300.0), source: Vec::new() };
    let dom3 = mesh3.cell_labels[0];
    let mut tr_spec =
        Problem::new(Study::HeatTransient { t_end, n_out, t0: format!("{t_init}") }, 1).with_cp_curve(dom3, cp_curve);
    tr_spec.materials.default.insert("k".into(), k_high);
    tr_spec.materials.default.insert("rho".into(), rho);
    tr_spec.materials.default.insert("Cp".into(), cp0);
    tr_spec.bcs.push(Bc {
        labels: all_outer,
        kind: BcKind::Radiation { emissivity: eps, t_ambient: t_amb, view_factor: 1.0 },
    });

    let tr_res = solve_heat(&mesh3, &tr_spec).expect("transient BDF2 RadiationBC solve");
    assert_eq!(tr_res.solutions.len(), n_out + 1);
    let t_final_vals = &tr_res.solutions.last().unwrap().field("T").unwrap().values;
    let t_final = t_final_vals.iter().sum::<f64>() / (t_final_vals.len() as f64);
    assert!(t_final > t_amb && t_final < t_init, "T_final = {t_final} must lie strictly in ({t_amb}, {t_init})");

    // Exact closed-form radiative cooling time integral t(T) for dT/dt = -(A * eps * sigma / (rho * Cp * V)) * (T^4 - T_amb^4):
    let vol = side * side * side;
    let area = 6.0 * side * side;
    let tau_coeff = (rho * cp0 * vol) / (4.0 * area * eps * SIGMA_SB * t_amb.powi(3));
    let exact_time_for_temp = |tk: f64| -> f64 {
        let log_term = ((tk + t_amb) * (t_init - t_amb) / ((tk - t_amb) * (t_init + t_amb))).ln();
        let atan_term = 2.0 * ((tk / t_amb).atan() - (t_init / t_amb).atan());
        tau_coeff * (log_term + atan_term)
    };
    let t_implied = exact_time_for_temp(t_final);
    let rel_time_err = (t_implied - t_end).abs() / t_end;
    assert!(
        rel_time_err < 2e-3,
        "BDF2 radiative cooling trajectory mismatch: implied t(T_final={t_final:.4} K) = {t_implied:.4} s vs {t_end} s (rel err {rel_time_err:.3e})"
    );
}

#[test]
fn conduction_analytical_joule_heating_and_power_conservation_2d_and_3d() {
    // 1. 2D rectangular bar [0, L] x [0, W] at order 1 and order 2 with Voltage & Ground BCs.
    let (l, w) = (4.0e-3, 1.5e-3);
    let mesh2 = rect(l, w, 12, 5);
    let mut left2 = 1;
    let mut right2 = 2;
    for f in 0..mesh2.n_facets() {
        let fv = mesh2.facet(f);
        let p0 = mesh2.point(fv[0] as usize);
        let p1 = mesh2.point(fv[1] as usize);
        if p0[0].abs() < 1e-12 && p1[0].abs() < 1e-12 {
            left2 = mesh2.facet_labels[f];
        } else if (p0[0] - l).abs() < 1e-12 && (p1[0] - l).abs() < 1e-12 {
            right2 = mesh2.facet_labels[f];
        }
    }
    let sigma = 1.25e4;
    let v_in = 2.4;
    let p_exact_2d = sigma * (w / l) * v_in * v_in;
    let q_exact_2d = sigma * (v_in / l).powi(2);

    for order in [1usize, 2usize] {
        let mut spec2 = Problem::new(Study::Conduction, order);
        spec2.materials.default.insert("sigma".into(), sigma);
        spec2.bcs = vec![
            Bc { labels: vec![left2], kind: BcKind::Ground },
            Bc { labels: vec![right2], kind: BcKind::Voltage(v_in) },
        ];
        let res2 = solve_conduction(&mesh2, &spec2).expect("2D conduction solve");
        let p_scalar = num_scalar(&res2, "P_joule (W)");
        assert!(
            ((p_scalar - p_exact_2d) / p_exact_2d).abs() < 1e-11,
            "2D order {order} P_joule {p_scalar} vs exact {p_exact_2d}"
        );
        let q_nodal = &res2.solutions[0].field("Q_joule").expect("Q_joule field").values;
        assert_eq!(q_nodal.len(), mesh2.n_dofs(order));
        let m_mat = assemble_mass(&mesh2, order, &vec![1.0; mesh2.n_cells()]);
        let p_int: f64 = m_mat.matvec(q_nodal).iter().sum();
        assert!(
            ((p_int - p_exact_2d) / p_exact_2d).abs() < 1e-10,
            "2D order {order} integrated Q_joule {p_int} vs exact {p_exact_2d}"
        );
        for &qv in q_nodal {
            assert!(((qv - q_exact_2d) / q_exact_2d).abs() < 1e-9, "pointwise Q_joule {qv} vs {q_exact_2d}");
        }
    }

    // 2. 3D block [0, Lx] x [0, Ly] x [0, Lz] with CurrentDensity BC on x = Lx and Ground on x = 0,
    //    plus spatially varying imported temperature T(x) driving sigma(T).
    let (lx, ly, lz) = (3.0e-3, 1.0e-3, 0.8e-3);
    let mesh3 = block(lx, ly, lz, 9, 3, 3);
    let mut left3 = 1;
    let mut right3 = 4;
    for f in 0..mesh3.n_facets() {
        let fv = mesh3.facet(f);
        if fv.iter().all(|&v| mesh3.point(v as usize)[0].abs() < 1e-12) {
            left3 = mesh3.facet_labels[f];
        } else if fv.iter().all(|&v| (mesh3.point(v as usize)[0] - lx).abs() < 1e-12) {
            right3 = mesh3.facet_labels[f];
        }
    }
    let j_in = 5.0e5; // A/m^2 inward at x = Lx
    let area = ly * lz;
    let current_i = j_in * area;
    let r_exact_3d = lx / (sigma * area);
    let p_exact_3d = current_i * current_i * r_exact_3d;

    let mut spec3 = Problem::new(Study::Conduction, 2);
    spec3.materials.default.insert("sigma".into(), sigma);
    spec3.bcs = vec![
        Bc { labels: vec![left3], kind: BcKind::Ground },
        Bc { labels: vec![right3], kind: BcKind::CurrentDensity(j_in) },
    ];
    let native_mesh3 = Mesh::from_solver_mesh(&mesh3).expect("valid 3D mesh");
    let prog = femcore::progress::Progress::new();
    let res3 = run(&native_mesh3, &spec3, &prog).expect("3D conduction run");
    let p3_scalar = num_scalar(&res3, "P_joule (W)");
    let v3_max = num_scalar(&res3, "V max (V)");
    assert!(
        ((v3_max - current_i * r_exact_3d) / (current_i * r_exact_3d)).abs() < 1e-10,
        "3D terminal voltage {v3_max} vs I*R = {}",
        current_i * r_exact_3d
    );
    assert!(
        ((p3_scalar - p_exact_3d) / p_exact_3d).abs() < 1e-10,
        "3D CurrentDensity P_joule {p3_scalar} vs I^2*R = {p_exact_3d}"
    );
    let q3_nodal = &res3.solutions[0].field("Q_joule").unwrap().values;
    let m3 = assemble_mass(&mesh3, 2, &vec![1.0; mesh3.n_cells()]);
    let p3_int: f64 = m3.matvec(q3_nodal).iter().sum();
    assert!(((p3_int - p_exact_3d) / p_exact_3d).abs() < 1e-10);

    // 3. Temperature-dependent sigma(T(x)) = sigma_0 * (1 + alpha * x / Lx):
    //    Exact series resistance R = int_0^Lx dx / (A * sigma(x)) = Lx * ln(1 + alpha) / (A * sigma_0 * alpha).
    let alpha = 1.0; // sigma doubles from x=0 (300 K) to x=Lx (400 K)
    let t_field: Vec<f64> = (0..mesh3.n_points()).map(|v| 300.0 + 100.0 * (mesh3.point(v)[0] / lx)).collect();
    let sig_curve = TempCurve::Tabulated {
        t: vec![300.0, 400.0],
        values: vec![sigma, sigma * (1.0 + alpha)],
        t_ref: Some(300.0),
        source: Vec::new(),
    };
    let dom3 = mesh3.cell_labels[0];
    let mut spec3_t = Problem::new(Study::Conduction, 2)
        .with_sigma_curve(dom3, TempTensor::Scalar(sig_curve), [0.0; 3])
        .with_temperature_field(t_field);
    spec3_t.materials.default.insert("sigma".into(), sigma);
    spec3_t.bcs = vec![
        Bc { labels: vec![left3], kind: BcKind::Ground },
        Bc { labels: vec![right3], kind: BcKind::Voltage(v_in) },
    ];
    let res3_t = solve_conduction(&mesh3, &spec3_t).expect("3D sigma(T) conduction");
    let r_exact_t = lx * (1.0 + alpha).ln() / (area * sigma * alpha);
    let p_exact_t = (v_in * v_in) / r_exact_t;
    let p_num_t = num_scalar(&res3_t, "P_joule (W)");
    assert!(((p_num_t - p_exact_t) / p_exact_t).abs() < 1e-5, "3D sigma(T) P_joule {p_num_t} vs exact {p_exact_t}");
}
