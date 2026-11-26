// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 validation suite for NAFEMS T1-T4 thermal benchmarks & coupled thermomechanical:
//! - NAFEMS T1: Membrane with hot spot (quarter symmetry, plane stress, sigma_zz at Point D)
//! - NAFEMS T2: 1D bar with nonlinear Stefan-Boltzmann radiation (Newton Robin iterations on `dfs::heat`)
//! - NAFEMS T3: 1D transient heat conduction with sinusoidal boundary drive (BDF2 over `nyx-dfs` matrices)
//! - NAFEMS T4: 2D steady heat conduction with mixed Dirichlet / insulated / convective cooling
//! - Coupled Timoshenko (1925) Si/Al bimetallic cantilever thermal deflection via `dfs::heat -> sms::static`

mod common;

use common::{cmd, load_validation_json, u_num_three_grid, u_num_two_grid, verify_asme_vv20, ValidationTarget};
use femcore::fem::assemble_mass;
use femcore::linalg::driver::DirectSolver;
use nyx_dfs::heat::{assemble_heat, heat_fields};
use nyx_dfs::problem::{Bc, BcKind, Problem, Study};
use nyx_engine::Engine;
use nyx_physics::discretise;
use serde::Deserialize;
use serde_json::json;
use std::f64::consts::PI;

#[derive(Deserialize)]
struct ThermalDataset {
    nafems_t1_hotspot: NafemsT1,
    nafems_t2_radiation: NafemsT2,
    nafems_t3_transient: NafemsT3,
    nafems_t4_convection: NafemsT4,
    coupled_bimorph_si_al: BimorphSection,
}

#[derive(Deserialize)]
struct NafemsT1 {
    plate_quarter_m: f64,
    hotspot_radius_m: f64,
    e_pa: f64,
    nu: f64,
    alpha_per_k: f64,
    delta_t_k: f64,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct NafemsT2 {
    length_m: f64,
    width_m: f64,
    k_w_mk: f64,
    t_a_k: f64,
    t_inf_k: f64,
    emissivity: f64,
    sigma_sb: f64,
    target: ValidationTarget,
}

#[derive(Deserialize)]
struct NafemsT3 {
    length_m: f64,
    width_m: f64,
    k_w_mk: f64,
    rho_kg_m3: f64,
    cp_j_kgk: f64,
    t_end_s: f64,
    period_half_s: f64,
    amp_c: f64,
    probe_x_m: f64,
    target: ValidationTarget,
}

#[derive(Deserialize)]
struct NafemsT4 {
    width_m: f64,
    height_m: f64,
    k_w_mk: f64,
    t_base_c: f64,
    h_w_m2k: f64,
    t_ext_c: f64,
    probe_e_m: [f64; 2],
    target: ValidationTarget,
}

#[derive(Deserialize)]
struct MatThermMech {
    e_pa: f64,
    nu: f64,
    alpha: f64,
    k: f64,
}

#[derive(Deserialize)]
struct BimorphSection {
    length_m: f64,
    h_si_m: f64,
    h_al_m: f64,
    si: MatThermMech,
    al: MatThermMech,
    t_ref_k: f64,
    t_hot_k: f64,
    target: ValidationTarget,
}

fn solve_t1_membrane(t1: &NafemsT1, h_global: f64, h_hot: f64) -> (f64, f64) {
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "NAFEMS T1 Hot-Spot", "dim": 2, "h": h_global }));
    cmd(
        &mut e,
        "add_rect",
        json!({ "name": "plate_box", "min": [0.0, 0.0], "max": [t1.plate_quarter_m, t1.plate_quarter_m] }),
    );
    cmd(&mut e, "add_circle", json!({ "name": "ring_circle", "center": [0.0, 0.0], "r": 1.35 * t1.hotspot_radius_m }));
    cmd(&mut e, "add_intersect", json!({ "name": "ring_q", "ids": ["plate_box", "ring_circle"] }));
    cmd(&mut e, "add_circle", json!({ "name": "hot_circle", "center": [0.0, 0.0], "r": t1.hotspot_radius_m }));
    cmd(&mut e, "add_intersect", json!({ "name": "hot_q", "ids": ["plate_box", "hot_circle"] }));
    cmd(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "outer_ring", "input": "plate_box" }, { "name": "cold_ring", "input": "ring_q" }, { "name": "hot_spot", "input": "hot_q" }] }),
    );
    cmd(&mut e, "add_size", json!({ "name": "sz_hot", "sel": "region(hot_spot) | region(cold_ring)", "h": h_hot }));
    cmd(&mut e, "set_mesh", json!({ "h": h_global, "max_dev": 2.0e-6 }));
    cmd(&mut e, "add_physics", json!({ "physics": "sms", "study": "static", "order": 2 }));
    cmd(
        &mut e,
        "set_physics_params",
        json!({
            "physics": "sms",
            "params": { "plane": "stress", "thermal": "uniform", "temperature": 293.15 + t1.delta_t_k, "t_ref": 293.15 }
        }),
    );
    cmd(
        &mut e,
        "set_material",
        json!({ "physics": "sms", "region": "hot_spot", "material": { "E": t1.e_pa, "nu": t1.nu, "alpha": t1.alpha_per_k, "rho": 7800.0 } }),
    );
    for reg in ["cold_ring", "outer_ring"] {
        cmd(
            &mut e,
            "set_material",
            json!({ "physics": "sms", "region": reg, "material": { "E": t1.e_pa, "nu": t1.nu, "alpha": 0.0, "rho": 7800.0 } }),
        );
    }
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0 } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(z = 0)", "kind": "fixed", "params": { "uz": 0.0 } }),
    );
    cmd(&mut e, "solve_run", json!({}));

    let r_out = t1.hotspot_radius_m + 3.0e-6;
    let pr = cmd(
        &mut e,
        "probe",
        json!({
            "physics": "sms",
            "expr": "szz",
            "points": [[r_out, 3.0e-6], [0.5 * t1.hotspot_radius_m, 0.5 * t1.hotspot_radius_m]]
        }),
    );
    let szz_out = pr["values"][0].as_f64().unwrap();
    let pr_in = cmd(
        &mut e,
        "probe",
        json!({
            "physics": "sms",
            "expr": "sxx",
            "points": [[0.4 * t1.hotspot_radius_m, 0.4 * t1.hotspot_radius_m]]
        }),
    );
    let sxx_in = pr_in["values"][0].as_f64().unwrap();
    (szz_out, sxx_in)
}

#[test]
fn nafems_t1_membrane_hotspot_vv20() {
    let ds: ThermalDataset = load_validation_json("thermal_nafems_benchmarks.json");
    let t1 = &ds.nafems_t1_hotspot;

    let (d3, in3) = solve_t1_membrane(t1, 1.2e-3, 56.0e-6);
    let (d2, in2) = solve_t1_membrane(t1, 0.9e-3, 36.0e-6);
    let (d1, in1) = solve_t1_membrane(t1, 0.65e-3, 24.0e-6);
    let (u_d, _) = u_num_three_grid(d1, d2, d3, 36.0 / 24.0, 56.0 / 36.0, 1.0);
    let (u_in, _) = u_num_three_grid(in1, in2, in3, 36.0 / 24.0, 56.0 / 36.0, 2.0);

    verify_asme_vv20(&t1.targets[0], d1, u_d);
    verify_asme_vv20(&t1.targets[1], in1, u_in);
}

fn solve_t2_radiation_bar(t2: &NafemsT2, h_mesh: f64) -> f64 {
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "NAFEMS T2 Radiation", "dim": 2, "h": h_mesh }));
    cmd(&mut e, "add_rect", json!({ "name": "bar", "min": [0.0, 0.0], "max": [t2.length_m, t2.width_m] }));
    cmd(&mut e, "add_overlay", json!({ "regions": [{ "name": "bar", "input": "bar" }] }));
    cmd(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "bar", "material": { "k": t2.k_w_mk } }));
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": t2.t_a_k } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": format!("outer & plane(x = {:e})", t2.length_m),
            "kind": "convection",
            "params": { "h": 1.0, "t_ext": t2.t_inf_k }
        }),
    );

    // Newton-Raphson linearization of grey-body Stefan-Boltzmann radiation BC using linear dfs::heat Robin solves
    let mut tb = 900.0_f64;
    for _ in 0..10 {
        let h_rad = 4.0 * t2.emissivity * t2.sigma_sb * tb.powi(3);
        let t_eff = (3.0 * tb.powi(4) + t2.t_inf_k.powi(4)) / (4.0 * tb.powi(3));
        cmd(&mut e, "remove_condition", json!({ "physics": "dfs", "index": 1 }));
        cmd(
            &mut e,
            "add_condition",
            json!({
                "physics": "dfs",
                "on": format!("outer & plane(x = {:e})", t2.length_m),
                "kind": "convection",
                "params": { "h": h_rad, "t_ext": t_eff }
            }),
        );
        cmd(&mut e, "solve_run", json!({}));
        let pr =
            cmd(&mut e, "probe", json!({ "physics": "dfs", "expr": "T", "points": [[t2.length_m, 0.5 * t2.width_m]] }));
        let tb_new = pr["values"][0].as_f64().unwrap();
        let diff = (tb_new - tb).abs();
        tb = tb_new;
        if diff < 1.0e-8 {
            break;
        }
    }
    tb
}

#[test]
fn nafems_t2_radiation_bar_vv20() {
    let ds: ThermalDataset = load_validation_json("thermal_nafems_benchmarks.json");
    let t2 = &ds.nafems_t2_radiation;

    let tb3 = solve_t2_radiation_bar(t2, 0.01);
    let tb2 = solve_t2_radiation_bar(t2, 0.005);
    let tb1 = solve_t2_radiation_bar(t2, 0.0025);
    let (u_num, _) = u_num_three_grid(tb1, tb2, tb3, 2.0, 2.0, 2.0);
    verify_asme_vv20(&t2.target, tb1, u_num);
}

fn solve_t3_sinusoidal_bdf2(t3: &NafemsT3, nx: usize, n_steps: usize) -> f64 {
    let sm = femcore::testkit::rect_split_x(
        t3.length_m,
        t3.width_m,
        t3.probe_x_m,
        (t3.length_m / (nx as f64)).min(0.5 * t3.width_m),
    );
    let mesh_ex = nyx_physics::Mesh::from_solver_mesh(&sm).unwrap();
    let disc = discretise(&mesh_ex).unwrap();

    let mut prob = Problem::new(Study::HeatTransient { t_end: t3.t_end_s, n_out: 1, t0: "0.0".into() }, 2);
    prob.materials.default.insert("k".into(), t3.k_w_mk);
    prob.materials.default.insert("rho".into(), t3.rho_kg_m3);
    prob.materials.default.insert("Cp".into(), t3.cp_j_kgk);
    // Labels in femcore::testkit::rect: 1=left(x=0), 4=right(x=L), 2=bottom, 3=top
    prob.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(0.0) },
        Bc { labels: vec![4], kind: BcKind::Temperature(1.0) },
    ];

    let sys = assemble_heat(&disc, &prob).unwrap();
    let n = sys.k.nrows();
    let rcp = vec![t3.rho_kg_m3 * t3.cp_j_kgk; sm.n_cells()];
    let m = femcore::fem::assemble_mass(&sm, 2, &rcp);

    // Partition DOFs into x=0 fixed (label 1), x=L driven (label 4), and interior free DOFs
    let coords = sm.dof_coordinates(2);
    let mut left_dofs = Vec::new();
    let mut right_dofs = Vec::new();
    let mut free = Vec::new();
    for (i, c) in coords.as_chunks::<2>().0.iter().enumerate() {
        if c[0] < 1.0e-9 {
            left_dofs.push(i);
        } else if (c[0] - t3.length_m).abs() < 1.0e-9 {
            right_dofs.push(i);
        } else {
            free.push(i);
        }
    }
    let nf = free.len();
    let kff = sys.k.submatrix(&free, &free);
    let mff = m.submatrix(&free, &free);
    let kfd = sys.k.submatrix(&free, &right_dofs);
    let mfd = m.submatrix(&free, &right_dofs);

    let dt = t3.t_end_s / (n_steps as f64);
    let a1 = mff.scaled(1.0 / dt).add_scaled(&kff, 1.0);
    let a2 = mff.scaled(1.5 / dt).add_scaled(&kff, 1.0);
    let s1 = DirectSolver::new(&a1, true).unwrap();
    let s2 = DirectSolver::new(&a2, true).unwrap();
    let ones_r = vec![1.0; right_dofs.len()];
    let k_col = kfd.matvec(&ones_r);
    let m_col = mfd.matvec(&ones_r);

    let bc_val = |t: f64| t3.amp_c * (PI * t / t3.period_half_s).sin();
    let mut uf_prev = vec![0.0; nf];
    let mut uf_curr = vec![0.0; nf];
    for step in 1..=n_steps {
        let t_now = (step as f64) * dt;
        let g_now = bc_val(t_now);
        if step == 1 {
            let dg_dt = g_now / dt;
            let mu = mff.matvec(&uf_curr);
            let rhs: Vec<f64> = (0..nf).map(|a| mu[a] / dt - g_now * k_col[a] - dg_dt * m_col[a]).collect();
            uf_prev = uf_curr;
            uf_curr = s1.solve(&rhs);
        } else {
            let g_p1 = bc_val((step as f64 - 1.0) * dt);
            let g_p2 = bc_val((step as f64 - 2.0) * dt);
            let dg_dt = (1.5 * g_now - 2.0 * g_p1 + 0.5 * g_p2) / dt;
            let mu1 = mff.matvec(&uf_curr);
            let mu2 = mff.matvec(&uf_prev);
            let rhs: Vec<f64> =
                (0..nf).map(|a| (2.0 * mu1[a] - 0.5 * mu2[a]) / dt - g_now * k_col[a] - dg_dt * m_col[a]).collect();
            let u_next = s2.solve(&rhs);
            uf_prev = uf_curr;
            uf_curr = u_next;
        }
    }

    let mut full = vec![0.0; n];
    for (a, &i) in free.iter().enumerate() {
        full[i] = uf_curr[a];
    }
    let g_end = bc_val(t3.t_end_s);
    for &i in &right_dofs {
        full[i] = g_end;
    }

    // Build T field through nyx_dfs::heat_fields and read DOF at (x = 0.08, z = W/2)
    let k_mat = vec![[[t3.k_w_mk, 0.0, 0.0], [0.0, t3.k_w_mk, 0.0], [0.0, 0.0, t3.k_w_mk]]; sm.n_cells()];
    let (fields, _) = heat_fields(&disc, &full, &k_mat, 2);
    let _ = assemble_mass;
    let idx = coords
        .as_chunks::<2>()
        .0
        .iter()
        .enumerate()
        .min_by(|a, b| {
            let da = (a.1[0] - t3.probe_x_m).hypot(a.1[1] - 0.5 * t3.width_m);
            let db = (b.1[0] - t3.probe_x_m).hypot(b.1[1] - 0.5 * t3.width_m);
            da.total_cmp(&db)
        })
        .unwrap()
        .0;
    assert!((coords[2 * idx] - t3.probe_x_m).abs() < 1e-12);
    let field_t = fields.iter().find(|f| f.name == "T").unwrap();
    field_t.values[idx]
}

#[test]
fn nafems_t3_transient_conduction_vv20() {
    let ds: ThermalDataset = load_validation_json("thermal_nafems_benchmarks.json");
    let t3 = &ds.nafems_t3_transient;

    // Spatial + temporal refinement ladder (10x20, 20x40, 40x80 steps -> 2x refinement ratio)
    let v3 = solve_t3_sinusoidal_bdf2(t3, 10, 32);
    let v2 = solve_t3_sinusoidal_bdf2(t3, 20, 64);
    let v1 = solve_t3_sinusoidal_bdf2(t3, 40, 128);
    let (u_num, _) = u_num_three_grid(v1, v2, v3, 2.0, 2.0, 2.0);
    verify_asme_vv20(&t3.target, v1, u_num);
}

fn solve_t4_convection(t4: &NafemsT4, h_mesh: f64) -> f64 {
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "NAFEMS T4 Convection", "dim": 2, "h": h_mesh }));
    cmd(&mut e, "add_rect", json!({ "name": "plate", "min": [0.0, 0.0], "max": [t4.width_m, t4.height_m] }));
    cmd(&mut e, "add_overlay", json!({ "regions": [{ "name": "plate", "input": "plate" }] }));
    cmd(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "plate", "material": { "k": t4.k_w_mk } }));
    // y=0 (plane z=0) fixed 100 C; x=0 insulated; x=W and y=H convect with h=750 to 0 C
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(z = 0)", "kind": "temperature", "params": { "value": t4.t_base_c } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": format!("outer & plane(x = {:e})", t4.width_m), "kind": "convection", "params": { "h": t4.h_w_m2k, "t_ext": t4.t_ext_c } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": format!("outer & plane(z = {:e})", t4.height_m), "kind": "convection", "params": { "h": t4.h_w_m2k, "t_ext": t4.t_ext_c } }),
    );
    cmd(&mut e, "solve_run", json!({}));

    let pr = cmd(&mut e, "probe", json!({ "physics": "dfs", "expr": "T", "points": [t4.probe_e_m] }));
    pr["values"][0].as_f64().unwrap()
}

#[test]
fn nafems_t4_2d_convection_vv20() {
    let ds: ThermalDataset = load_validation_json("thermal_nafems_benchmarks.json");
    let t4 = &ds.nafems_t4_convection;

    let t_c3 = solve_t4_convection(t4, 0.06);
    let t_c2 = solve_t4_convection(t4, 0.03);
    let t_c1 = solve_t4_convection(t4, 0.015);
    let (u_num, _) = u_num_three_grid(t_c1, t_c2, t_c3, 2.0, 2.0, 2.0);
    verify_asme_vv20(&t4.target, t_c1, u_num);
}

fn solve_coupled_bimorph(b: &BimorphSection, h_mesh: f64) -> f64 {
    let mut e = Engine::new();
    let h_tot = b.h_si_m + b.h_al_m;
    cmd(&mut e, "project_new", json!({ "title": "Coupled Si/Al Bimorph", "dim": 2, "h": h_mesh }));
    cmd(&mut e, "add_rect", json!({ "name": "bot_si", "min": [0.0, 0.0], "max": [b.length_m, b.h_si_m] }));
    cmd(&mut e, "add_rect", json!({ "name": "top_al", "min": [0.0, b.h_si_m], "max": [b.length_m, h_tot] }));
    cmd(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "si", "input": "bot_si" }, { "name": "al", "input": "top_al" }] }),
    );

    // Coupled thermal stage: uniform steady temperature t_hot_k driving mechanical expansion
    cmd(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 2 }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "si", "material": { "k": b.si.k } }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "al", "material": { "k": b.al.k } }));
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "temperature", "params": { "value": b.t_hot_k } }),
    );

    cmd(&mut e, "add_physics", json!({ "physics": "sms", "study": "static", "order": 2 }));
    cmd(
        &mut e,
        "set_physics_params",
        json!({ "physics": "sms", "params": { "plane": "stress", "thermal": "imported", "t_ref": b.t_ref_k } }),
    );
    cmd(
        &mut e,
        "set_material",
        json!({ "physics": "sms", "region": "si", "material": { "E": b.si.e_pa, "nu": b.si.nu, "alpha": b.si.alpha, "rho": 2330.0 } }),
    );
    cmd(
        &mut e,
        "set_material",
        json!({ "physics": "sms", "region": "al", "material": { "E": b.al.e_pa, "nu": b.al.nu, "alpha": b.al.alpha, "rho": 2700.0 } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    cmd(&mut e, "set_active", json!({ "physics": "sms" }));
    cmd(&mut e, "solve_run", json!({}));

    let pr = cmd(&mut e, "probe", json!({ "physics": "sms", "expr": "uz", "points": [[b.length_m, 0.5 * h_tot]] }));
    pr["values"][0].as_f64().unwrap().abs()
}

#[test]
fn coupled_thermomechanical_timoshenko_bimorph_vv20() {
    let ds: ThermalDataset = load_validation_json("thermal_nafems_benchmarks.json");
    let b = &ds.coupled_bimorph_si_al;

    let w_c = solve_coupled_bimorph(b, 0.35e-3);
    let w_f = solve_coupled_bimorph(b, 0.18e-3);
    let u_num = u_num_two_grid(w_c, w_f, 0.35 / 0.18, 2.0);
    verify_asme_vv20(&b.target, w_f, u_num);
}
