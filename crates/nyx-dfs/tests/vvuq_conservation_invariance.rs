// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME Layer 2 Pillars (d), (e), and R32 regression verification:
//! * Global charge conservation (`< 1e-10`) and 3+ terminal Maxwell C strict SPD, symmetry, row-sum identity;
//! * Stationary heat flux conservation (`< 1e-10`) with zero-net-source denominator safety;
//! * Magnetostatics 3+ terminal L symmetry & SPD, weak div(B)=0 (`< 1e-10`), and discrete surface-current
//!   conservation `||D_s K_s||_2 / ||K_s||_2 < 1e-10` on `planar_stream`;
//! * Rigid Euclidean mesh invariance (`< 1e-10`), geometry rotated remeshing (`< 1e-3`), vertex/cell
//!   permutation invariance (`< 1e-10`), terminal permutation `P M P^T` (`< 1e-12`), and 3D scaling (`< 1e-11`);
//! * R32 explicit `Error::Unsupported` rejection when curved meshes enter `nyx-dfs::{heat, magneto}`.

use femcore::fem::geometry;
use femcore::linalg::driver::dense_eigh;
use femcore::material::{rotation, Anisotropy, Tensor};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::progress::Progress;
use femcore::vvuq::FnSource;
use geom::expr::Expr;
use geom::model::Model;
use nyx_dfs::electro::{assemble_electrostatics_with_source, solve_potential};
use nyx_dfs::heat::{assemble_heat, assemble_heat_with_source, solve_heat_mixed, solve_heat_with_source};
use nyx_dfs::magneto::{solve_magnetostatics_mixed, MU0};
use nyx_dfs::{run, solve_electrostatics, solve_heat, solve_magnetostatics, Bc, BcKind, Problem, Study, EPS0};
use nyx_physics::{CapacitanceMatrixData, Error, InductanceMatrixData, Mesh};
use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::sync::Arc;

fn assert_spd_sym_matrix(label: &str, m: &[Vec<f64>], check_row_sum: bool) {
    let n = m.len();
    assert!(n >= 3, "{label}: need >= 3 terminals, got {n}");
    let mut frob2 = 0.0_f64;
    let mut max_asym = 0.0_f64;
    let mut flat = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            frob2 += m[i][j] * m[i][j];
            max_asym = max_asym.max((m[i][j] - m[j][i]).abs());
            flat[i * n + j] = 0.5 * (m[i][j] + m[j][i]);
        }
    }
    let frob = frob2.sqrt().max(1e-30);
    assert!(max_asym / frob < 1e-12, "{label}: asymmetry {:e} >= 1e-12", max_asym / frob);
    let (eigs, _) = dense_eigh(&flat, n).unwrap();
    assert!(eigs[0] > 0.0, "{label}: lambda_min = {:e} must be > 0", eigs[0]);
    if check_row_sum {
        for (i, row) in m.iter().enumerate() {
            let s: f64 = row.iter().sum();
            assert!(s >= -1e-12 * frob, "{label}: row {i} sum {s:e} < -1e-12 ||C||");
        }
    }
}

#[test]
fn electro_global_charge_conservation_and_3terminal_maxwell_spd() {
    // 1. 3D global charge conservation with non-uniform rho_v, Dirichlet + Neumann + Open BCs (< 1e-10)
    let mesh = femcore::testkit::block(1.0, 0.8, 0.6, 5, 4, 4);
    let source = FnSource {
        volume_rhs: |x: [f64; 3], _t: f64, _d: u32| [EPS0 * (10.0 + 3.0 * (PI * x[0]).sin() * x[1]), 0.0, 0.0],
        dirichlet: |x: [f64; 3], _t: f64, bnd: u32| if bnd == 1 { [1.5 + 0.2 * x[1], 0.0, 0.0] } else { [0.0; 3] },
        neumann_robin: |x: [f64; 3], _n: [f64; 3], _t: f64, _bnd: u32| [EPS0 * (0.5 + x[0] * x[2]), 0.0, 0.0],
    };
    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.bcs = vec![
        Bc { labels: vec![1, 2], kind: BcKind::Potential(0.0) },
        Bc { labels: vec![3], kind: BcKind::SurfaceCharge(0.0) },
        Bc { labels: vec![6], kind: BcKind::Open { origin: Some([0.5, 0.4, -0.5]) } },
    ];
    let asm = assemble_electrostatics_with_source::<3, _>(&mesh, &prob, &source).unwrap();
    let (v, _) = solve_potential(&mesh, &prob, &asm).unwrap();

    // Also assemble the open-free operator (without Open Robin) to isolate Q_open = 1^T (K_open * V)
    let mut prob_no_open = prob.clone();
    prob_no_open.bcs.pop();
    let asm_no_open = assemble_electrostatics_with_source::<3, _>(&mesh, &prob_no_open, &source).unwrap();

    let kv_asm = asm.k.matvec(&v);
    let kv_vol = asm_no_open.k.matvec(&v);
    let fixed_set: BTreeSet<usize> = asm.fixed.iter().copied().collect();
    let mut q_dir = 0.0;
    let mut q_dir_abs = 0.0;
    for &d in &asm.fixed {
        let r_d = kv_asm[d] - asm.f[d];
        q_dir += r_d;
        q_dir_abs += r_d.abs();
    }
    // Check free rows residual ~ 0
    for (i, (&kv_i, &f_i)) in kv_asm.iter().zip(&asm.f).enumerate() {
        if !fixed_set.contains(&i) {
            assert!((kv_i - f_i).abs() < 1e-10 * (q_dir_abs.max(1e-20)));
        }
    }
    let q_open: f64 = kv_asm.iter().zip(&kv_vol).map(|(a, b)| a - b).sum();
    let f_tot: f64 = asm.f.iter().sum();
    let f_abs: f64 = asm.f.iter().map(|x| x.abs()).sum();
    let denom = q_dir_abs.max(f_abs).max(1e-30);
    // Outward flux residual: (-Q_dir) + Q_open = int rho_v dV + int g_N dS (= f_tot)
    let cons_err = (q_open - q_dir - f_tot).abs() / denom;
    eprintln!("Electro 3D global charge conservation residual = {cons_err:e}");
    assert!(cons_err < 1e-10, "Electro 3D charge conservation {cons_err:e} >= 1e-10");

    // 2. 3-terminal Maxwell C matrix: strict SPD, symmetry < 1e-12, row-sum identity
    let mut p3 = Problem::new(Study::Electrostatic, 2);
    p3.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("t1".into()) },
        Bc { labels: vec![2], kind: BcKind::Terminal("t2".into()) },
        Bc { labels: vec![3], kind: BcKind::Terminal("t3".into()) },
        Bc { labels: vec![6], kind: BcKind::Ground },
    ];
    let c3 = CapacitanceMatrixData::of(&solve_electrostatics(&mesh, &p3).unwrap()).unwrap().maxwell_matrix_f;
    assert_spd_sym_matrix("Electro 3-terminal Maxwell C", &c3, true);
}

#[test]
fn heat_stationary_flux_conservation_and_zero_net_source_safety() {
    let mesh = femcore::testkit::block(1.0, 0.8, 0.5, 4, 4, 4);
    let q_src = 24.0;
    let k_cond = 3.0;
    // Exact quadratic profile along z with mixed Dirichlet (z=0), Neumann (z=H), Convection (x,y sides):
    // T(z) = 300 + 20 z - 0.5 (Q / k) z^2  =>  q_z = -k (20 - Q/k * z) = -k*20 + Q z
    let h_val = 0.5;
    let mut prob = Problem::new(Study::HeatStationary, 2);
    prob.materials.default.insert("k".into(), k_cond);
    prob.materials.default.insert("Q".into(), q_src);
    // At z=H: +k dT/dz = k*(20 - 8*0.5) = 48 W/m^2 incoming; sides insulated except bottom T=300
    prob.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Flux(k_cond * 20.0 - q_src * h_val) },
        Bc { labels: vec![1, 2, 5, 6], kind: BcKind::Flux(0.0) },
    ];
    let res = solve_heat(&mesh, &prob).unwrap();
    let sol = &res.solutions[0];
    let q_cells = [
        sol.cell_fields.iter().find(|f| f.name == "qx").unwrap(),
        sol.cell_fields.iter().find(|f| f.name == "qy").unwrap(),
        sol.cell_fields.iter().find(|f| f.name == "qz").unwrap(),
    ];

    // Integrate exact P2 surface flux int_{dOmega} (-q . n) dS? Wait: div q = Q => int_{dOmega} q . n dS = int_Omega Q dV
    // Recall heat_fields defines q = -k grad T, so -div(k grad T) = Q => div(q) = Q => int_{dOmega} q . n dS = int_Omega Q dV
    let opp_map = {
        let mut map = std::collections::HashMap::new();
        for c in 0..mesh.n_cells() {
            let cv = mesh.cell(c);
            for skip in 0..4 {
                let tri = match skip {
                    0 => [cv[1], cv[2], cv[3]],
                    1 => [cv[0], cv[2], cv[3]],
                    2 => [cv[0], cv[1], cv[3]],
                    _ => [cv[0], cv[1], cv[2]],
                };
                let mut key = tri;
                key.sort_unstable();
                map.insert(key, (c, cv[skip]));
            }
        }
        map
    };
    let mut bnd_flux = 0.0;
    let mut bnd_flux_abs = 0.0;
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let mut key = [fv[0], fv[1], fv[2]];
        key.sort_unstable();
        let &(c, opp) = opp_map.get(&key).unwrap();
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        let p2 = mesh.point(fv[2] as usize);
        let po = mesh.point(opp as usize);
        let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let cx = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let nlen = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
        let area = 0.5 * nlen;
        let mut n_hat = [cx[0] / nlen, cx[1] / nlen, cx[2] / nlen];
        if (p0[0] - po[0]) * n_hat[0] + (p0[1] - po[1]) * n_hat[1] + (p0[2] - po[2]) * n_hat[2] < 0.0 {
            n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
        }
        // Evaluate linear-in-tet q at facet barycenter (q is exact linear for quadratic T in P2)
        let cv = mesh.cell(c);
        let mut q_avg = [0.0; 3];
        let nloc = q_cells[0].nloc;
        for &fv_i in fv {
            let loc = cv.iter().position(|&u| u == fv_i).unwrap();
            for k in 0..3 {
                q_avg[k] += q_cells[k].re[c * nloc + loc] / 3.0;
            }
        }
        let qn = q_avg[0] * n_hat[0] + q_avg[1] * n_hat[1] + q_avg[2] * n_hat[2];
        bnd_flux += area * qn;
        bnd_flux_abs += area * qn.abs();
    }
    let q_vol_tot = q_src * (1.0 * 0.8 * 0.5);
    let cons = (bnd_flux - q_vol_tot).abs() / bnd_flux_abs.max(q_vol_tot.abs()).max(1e-30);
    eprintln!("Heat 3D stationary surface flux conservation = {cons:e}");
    assert!(cons < 1e-10, "Heat 3D surface flux conservation {cons:e} >= 1e-10");

    // 2. Variational reaction conservation with arbitrary smooth FnSource + Convection (< 1e-12)
    let source = FnSource {
        volume_rhs: |x: [f64; 3], _t: f64, _d: u32| [5.0 + 3.0 * (PI * x[0]).cos() * x[2], 0.0, 0.0],
        dirichlet: |_x: [f64; 3], _t: f64, _b: u32| [310.0, 0.0, 0.0],
        neumann_robin: |x: [f64; 3], _n: [f64; 3], _t: f64, b: u32| {
            if b == 6 {
                [2.0 * (300.0 + x[0]), 0.0, 0.0]
            } else {
                [0.0; 3]
            }
        },
    };
    let mut p_mms = Problem::new(Study::HeatStationary, 2);
    p_mms.materials.default.insert("k".into(), k_cond);
    p_mms.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Temperature(310.0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h: 2.0, text: 300.0 } },
    ];
    let sys = assemble_heat_with_source::<3, _>(&mesh, &p_mms, &source).unwrap();
    let mut p_vol = p_mms.clone();
    p_vol.bcs.pop(); // remove Convection to get K_vol alone
    let sys_vol = assemble_heat_with_source::<3, _>(&mesh, &p_vol, &source).unwrap();
    let rmms = solve_heat_with_source::<3, _>(&mesh, &p_mms, &source).unwrap();
    let tv = &rmms.solutions[0].field("T").unwrap().values;
    let kt = sys.k.matvec(tv);
    let kt_vol = sys_vol.k.matvec(tv);
    let mut q_dir = 0.0;
    for &d in &sys.fixed {
        q_dir += kt[d] - sys.f[d];
    }
    let q_conv: f64 = kt.iter().zip(&kt_vol).map(|(a, b)| a - b).sum();
    let f_sum: f64 = sys.f.iter().sum();
    // 1^T K_vol T = 0 => (sum_fixed (K T - f)_d) - (1^T K_conv T) + f_sum = 0
    let res_var = (q_dir - q_conv + f_sum).abs() / f_sum.abs().max(1e-30);
    assert!(res_var < 1e-10, "Heat variational reaction residual {res_var:e}");

    // 3. Zero-net-source denominator safety check (Q = 0, uniform T = 300)
    let mut p_zero = Problem::new(Study::HeatStationary, 2);
    p_zero.materials.default.insert("k".into(), k_cond);
    p_zero.bcs = vec![Bc { labels: vec![5, 6], kind: BcKind::Temperature(300.0) }];
    let r_z = solve_heat(&mesh, &p_zero).unwrap();
    let tz = &r_z.solutions[0].field("T").unwrap().values;
    let sys_z = assemble_heat(&mesh, &p_zero).unwrap();
    let ktz = sys_z.k.matvec(tz);
    let num_z: f64 = ktz.iter().sum::<f64>().abs();
    let denom_z = sys_z.f.iter().map(|x| x.abs()).sum::<f64>().max(1e-30);
    assert!(num_z / denom_z < 1e-10 || num_z < 1e-10 * 300.0 * k_cond);
}

#[test]
fn magneto_3terminal_spd_weak_div_b_and_planar_stream_current_conservation() {
    // 1. 3-terminal cross_section SPD & symmetry < 1e-12
    let mut m2 = Model::new_2d();
    let box_d = m2.rect("air", (-40e-6, -40e-6), (40e-6, 40e-6));
    let w1 = m2.circle("w1", (-15e-6, 0.0), 4e-6);
    let w2 = m2.circle("w2", (0.0, 0.0), 4e-6);
    let w3 = m2.circle("w3", (15e-6, 0.0), 4e-6);
    use geom::model::RegionKind;
    m2.overlay_full(
        &[
            ("air", box_d, RegionKind::Material, None),
            ("w1", w1, RegionKind::Hole, None),
            ("w2", w2, RegionKind::Hole, None),
            ("w3", w3, RegionKind::Hole, None),
        ],
        &[],
    );
    let g2 = mesh_geom_2d(evaluate(&m2.list).unwrap(), &GeomMeshOpts::new(3.0e-6)).unwrap();
    let (mut labs1, mut labs2, mut labs3, mut gnd) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let ym = 0.5 * (p0[1] + p1[1]);
        let lab = g2.mesh.facet_labels[f];
        if (xm + 15e-6).hypot(ym) < 6e-6 {
            if !labs1.contains(&lab) {
                labs1.push(lab);
            }
        } else if xm.hypot(ym) < 6e-6 {
            if !labs2.contains(&lab) {
                labs2.push(lab);
            }
        } else if (xm - 15e-6).hypot(ym) < 6e-6 {
            if !labs3.contains(&lab) {
                labs3.push(lab);
            }
        } else if !gnd.contains(&lab) {
            gnd.push(lab);
        }
    }
    let mut p_3w = Problem::new(Study::Magnetostatic { formulation: "cross_section".into() }, 2);
    p_3w.bcs = vec![
        Bc { labels: gnd, kind: BcKind::Ground },
        Bc {
            labels: labs1,
            kind: BcKind::CurrentTerminal { name: "w1".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc {
            labels: labs2,
            kind: BcKind::CurrentTerminal { name: "w2".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc {
            labels: labs3,
            kind: BcKind::CurrentTerminal { name: "w3".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
    ];
    let l_3w = InductanceMatrixData::of(&solve_magnetostatics(&g2.mesh, &p_3w).unwrap()).unwrap().inductance_matrix_h;
    assert_spd_sym_matrix("Magneto 3-terminal cross_section L", &l_3w, false);

    // 2. Weak div(B) = 0 (< 1e-10) AND 3-terminal SPD on mixed_3d
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m3 = Model::new_3d();
    let sub = m3.box_("sub", (-60e-6, -50e-6, -30e-6), (60e-6, 50e-6, 0.0));
    let vac = m3.box_("vac", (-60e-6, -50e-6, 0.0), (60e-6, 50e-6, 30e-6));
    let mk_trace = |m: &mut Model, nm: &str, y0: f64, y1: f64| {
        let pname = format!("{nm}_p");
        let p = m.polygons(&pname, vec![vec![pt(-60e-6, y0), pt(60e-6, y0), pt(60e-6, y1), pt(-60e-6, y1)]], vec![]);
        m.sheet(nm, p, 0.0);
    };
    mk_trace(&mut m3, "t1", -30e-6, -20e-6);
    mk_trace(&mut m3, "t2", -5e-6, 5e-6);
    mk_trace(&mut m3, "t3", 20e-6, 30e-6);
    m3.overlay(&[("sub", sub), ("vac", vac)]);
    let g3 = mesh_model(&m3.list, &MeshPlan::simple(16e-6)).unwrap();
    let t1_l = g3.boundary("t1");
    let t2_l = g3.boundary("t2");
    let t3_l = g3.boundary("t3");
    let mut sh_l = t1_l.clone();
    sh_l.extend_from_slice(&t2_l);
    sh_l.extend_from_slice(&t3_l);
    let sh_set: std::collections::BTreeSet<i32> = sh_l.iter().copied().collect();
    let mut gnd3 = Vec::new();
    for (f, &is_e) in g3.mesh.exterior_facets().iter().enumerate() {
        let lab = g3.mesh.facet_labels[f];
        if is_e && !sh_set.contains(&lab) && !gnd3.contains(&lab) {
            gnd3.push(lab);
        }
    }
    let mut p_m3 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
    p_m3.bcs = vec![
        Bc { labels: sh_l, kind: BcKind::LondonSheet { l_k_sq: 1.0e-12 } },
        Bc { labels: gnd3, kind: BcKind::Ground },
        Bc {
            labels: t1_l,
            kind: BcKind::CurrentTerminal { name: "t1".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
        },
        Bc {
            labels: t2_l,
            kind: BcKind::CurrentTerminal { name: "t2".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
        },
        Bc {
            labels: t3_l,
            kind: BcKind::CurrentTerminal { name: "t3".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
        },
    ];
    let r_m3 = solve_magnetostatics(&g3.mesh, &p_m3).unwrap();
    let l_3m = InductanceMatrixData::of(&r_m3).unwrap().inductance_matrix_h;
    assert_spd_sym_matrix("Magneto 3-terminal mixed_3d L", &l_3m, false);

    // Check weak div(B) = 0 on interior vertices
    let sol0 = &r_m3.solutions[0];
    let bx = sol0.cell_fields.iter().find(|f| f.name == "Bx").unwrap();
    let by = sol0.cell_fields.iter().find(|f| f.name == "By").unwrap();
    let bz = sol0.cell_fields.iter().find(|f| f.name == "Bz").unwrap();
    let nloc = bx.nloc;
    let geo3 = geometry(&g3.mesh);
    let nv3 = g3.mesh.n_points();
    let mut on_bnd3 = vec![false; nv3];
    let ext3 = g3.mesh.exterior_facets();
    for (f, &is_e) in ext3.iter().enumerate() {
        if is_e {
            for &v in g3.mesh.facet(f) {
                on_bnd3[v as usize] = true;
            }
        }
    }
    let mut div_v = vec![0.0; nv3];
    let mut grad_norm2_v = vec![0.0; nv3];
    let mut b_l2_sq = 0.0_f64;
    for c in 0..g3.mesh.n_cells() {
        let cv = g3.mesh.cell(c);
        let vol = geo3.vol[c];
        let bvec = [bx.re[c * nloc], by.re[c * nloc], bz.re[c * nloc]];
        b_l2_sq += vol * (bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
        for (loc, &v) in cv.iter().enumerate() {
            let g_i = [geo3.g[c * 12 + loc * 3], geo3.g[c * 12 + loc * 3 + 1], geo3.g[c * 12 + loc * 3 + 2]];
            div_v[v as usize] += vol * (bvec[0] * g_i[0] + bvec[1] * g_i[1] + bvec[2] * g_i[2]);
            grad_norm2_v[v as usize] += vol * (g_i[0] * g_i[0] + g_i[1] * g_i[1] + g_i[2] * g_i[2]);
        }
    }
    let b_l2 = b_l2_sq.sqrt().max(1e-30);
    let mut max_weak_div = 0.0_f64;
    let mut worst_v = 0usize;
    for v in 0..nv3 {
        if !on_bnd3[v] && grad_norm2_v[v] > 0.0 {
            let norm_rel = div_v[v].abs() / (b_l2 * grad_norm2_v[v].sqrt());
            if norm_rel > max_weak_div {
                max_weak_div = norm_rel;
                worst_v = v;
            }
        }
    }
    eprintln!(
        "Magneto mixed_3d max interior weak div(B) = {max_weak_div:e} at v={worst_v} {:?} (on_bnd3={})",
        g3.mesh.point(worst_v),
        on_bnd3[worst_v]
    );
    assert!(max_weak_div < 1e-10, "mixed_3d weak div(B) {max_weak_div:e} >= 1e-10");

    // 3. planar_stream discrete surface-current conservation ||D_s K_s||_2 / ||K_s||_2 < 1e-10
    let mesh2 = femcore::testkit::rect(40e-6, 15e-6, 12, 6);
    let mut p_ps = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_ps.bcs = vec![
        Bc { labels: vec![2], kind: BcKind::Ground },
        Bc {
            labels: vec![3],
            kind: BcKind::CurrentTerminal { name: "strip".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
        },
    ];
    let r_ps = solve_magnetostatics(&mesh2, &p_ps).unwrap();
    let s0 = &r_ps.solutions[0];
    let kx = s0.cell_fields.iter().find(|f| f.name == "Bx").unwrap();
    let kz = s0.cell_fields.iter().find(|f| f.name == "Bz").unwrap();
    let geo2 = geometry(&mesh2);
    let nv2 = mesh2.n_points();
    let mut on_bnd2 = vec![false; nv2];
    for f in 0..mesh2.n_facets() {
        for &v in mesh2.facet(f) {
            on_bnd2[v as usize] = true;
        }
    }
    let mut ds_ks = vec![0.0; nv2];
    let mut ks_l2_sq = 0.0;
    for c in 0..mesh2.n_cells() {
        let cv = mesh2.cell(c);
        let area = geo2.vol[c];
        let kvec = [kx.re[c * kx.nloc] / MU0, kz.re[c * kz.nloc] / MU0];
        ks_l2_sq += area * (kvec[0] * kvec[0] + kvec[1] * kvec[1]);
        for (loc, &v) in cv.iter().enumerate() {
            let g_i = [geo2.g[c * 6 + loc * 2], geo2.g[c * 6 + loc * 2 + 1]];
            ds_ks[v as usize] += area * (kvec[0] * g_i[0] + kvec[1] * g_i[1]);
        }
    }
    let ds_norm = (0..nv2).filter(|&v| !on_bnd2[v]).map(|v| ds_ks[v] * ds_ks[v]).sum::<f64>().sqrt();
    let ks_rel = ds_norm / ks_l2_sq.sqrt().max(1e-30);
    eprintln!("Magneto planar_stream ||D_s K_s||_2 / ||K_s||_2 = {ks_rel:e}");
    assert!(ks_rel < 1e-10, "planar_stream current conservation {ks_rel:e} >= 1e-10");
}

fn rigid_transform_mesh3d(mesh: &SimplexMesh, euler_deg: [f64; 3], shift: [f64; 3]) -> SimplexMesh {
    let r = rotation(euler_deg);
    let pts: Vec<f64> = mesh
        .points
        .chunks(3)
        .flat_map(|p| {
            [
                r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2] + shift[0],
                r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2] + shift[1],
                r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2] + shift[2],
            ]
        })
        .collect();
    SimplexMesh::new(
        3,
        pts,
        mesh.cells.clone(),
        mesh.cell_labels.clone(),
        mesh.facets.clone(),
        mesh.facet_labels.clone(),
    )
}

fn scale_mesh(mesh: &SimplexMesh, s: f64) -> SimplexMesh {
    let pts: Vec<f64> = mesh.points.iter().map(|&x| x * s).collect();
    SimplexMesh::new(
        mesh.dim,
        pts,
        mesh.cells.clone(),
        mesh.cell_labels.clone(),
        mesh.facets.clone(),
        mesh.facet_labels.clone(),
    )
}

fn permute_mesh(mesh: &SimplexMesh) -> SimplexMesh {
    let d = mesh.dim;
    let nv = mesh.n_points();
    let nc = mesh.n_cells();
    let nloc = d + 1;
    // Deterministic permutation of vertices and cells
    let mut v_perm: Vec<usize> = (0..nv).collect();
    let mut state = 0x517cc1b727220a95u64;
    for i in (1..nv).rev() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let j = (state as usize) % (i + 1);
        v_perm.swap(i, j);
    }
    let mut v_inv = vec![0u32; nv];
    for (new_i, &old_i) in v_perm.iter().enumerate() {
        v_inv[old_i] = new_i as u32;
    }
    let mut new_pts = Vec::with_capacity(nv * d);
    for &old_i in &v_perm {
        new_pts.extend_from_slice(mesh.point(old_i));
    }
    let mut c_perm: Vec<usize> = (0..nc).collect();
    for i in (1..nc).rev() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let j = (state as usize) % (i + 1);
        c_perm.swap(i, j);
    }
    let mut new_cells = Vec::with_capacity(nc * nloc);
    let mut new_clabs = Vec::with_capacity(nc);
    for (k, &old_c) in c_perm.iter().enumerate() {
        let cv = mesh.cell(old_c);
        // Apply even cyclic permutation on local cell vertices preserving positive Jacobian orientation
        if d == 3 {
            let even_rot = [[0, 1, 2, 3], [1, 2, 0, 3], [2, 3, 1, 0], [3, 0, 2, 1]][k % 4];
            for idx in even_rot {
                new_cells.push(v_inv[cv[idx] as usize]);
            }
        } else {
            let sh = k % 3;
            for l in 0..3 {
                new_cells.push(v_inv[cv[(l + sh) % 3] as usize]);
            }
        }
        new_clabs.push(mesh.cell_labels[old_c]);
    }
    let mut new_facets = Vec::with_capacity(mesh.facets.len());
    for f in 0..mesh.n_facets() {
        for &v in mesh.facet(f) {
            new_facets.push(v_inv[v as usize]);
        }
    }
    SimplexMesh::new(d, new_pts, new_cells, new_clabs, new_facets, mesh.facet_labels.clone())
}

#[test]
fn invariances_rigid_euclidean_remesh_permutation_terminal_and_scaling() {
    // 1. Mesh-level rigid Euclidean (< 1e-10): rotation + translation on 3D anisotropic electro & heat + mixed_3d
    let base3 = femcore::testkit::block(1.0, 0.8, 0.5, 3, 3, 3);
    let euler = [37.0, 53.0, 19.0];
    let shift = [1.7, -0.8, 2.3];
    let rot3 = rigid_transform_mesh3d(&base3, euler, shift);
    let rmat = rotation(euler);

    // Electro anisotropic
    let eps_t = Tensor::Diagonal { x: 4.5, y: 7.2, z: 11.3 };
    let mut pe0 = Problem::new(Study::Electrostatic, 2);
    pe0.materials
        .aniso
        .insert(base3.cell_labels[0], Anisotropy { permittivity: Some(eps_t.clone()), ..Default::default() });
    pe0.bcs = vec![Bc { labels: vec![5], kind: BcKind::Potential(3.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    let mut per = pe0.clone();
    per.materials.aniso.get_mut(&base3.cell_labels[0]).unwrap().euler_deg = euler;
    let re0 = solve_electrostatics(&base3, &pe0).unwrap();
    let rer = solve_electrostatics(&rot3, &per).unwrap();
    let c_e0 = re0.scalar("capacitance").unwrap();
    let c_er = rer.scalar("capacitance").unwrap();
    assert!((c_er - c_e0).abs() / c_e0 < 1e-10, "Rigid Euclidean Electro C err");
    // Vector E field transformed by R: E_rot = R * E_0 (Field::vector stores interleaved XYZ per DOF)
    let e0 = re0.solutions[0].field("E").unwrap();
    let er = rer.solutions[0].field("E").unwrap();
    let n_pts = base3.n_dofs(2);
    for k in 0..n_pts {
        let v0 = [e0.values[3 * k], e0.values[3 * k + 1], e0.values[3 * k + 2]];
        let vr = [er.values[3 * k], er.values[3 * k + 1], er.values[3 * k + 2]];
        for i in 0..3 {
            let r_v0_i = rmat[i][0] * v0[0] + rmat[i][1] * v0[1] + rmat[i][2] * v0[2];
            assert!((vr[i] - r_v0_i).abs() < 1e-10 * 10.0, "E vector rigid rotation mismatch at dof {k}");
        }
    }

    // Heat anisotropic rigid Euclidean
    let k_t = Tensor::Diagonal { x: 2.0, y: 5.0, z: 9.0 };
    let mut ph0 = Problem::new(Study::HeatStationary, 2);
    ph0.materials.aniso.insert(base3.cell_labels[0], Anisotropy { conductivity: Some(k_t), ..Default::default() });
    ph0.materials.default.insert("Q".into(), 12.0);
    ph0.bcs = vec![
        Bc { labels: vec![5], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![6], kind: BcKind::Convection { h: 4.0, text: 290.0 } },
    ];
    let mut phr = ph0.clone();
    phr.materials.aniso.get_mut(&base3.cell_labels[0]).unwrap().euler_deg = euler;
    let rh0 = solve_heat(&base3, &ph0).unwrap();
    let rhr = solve_heat(&rot3, &phr).unwrap();
    assert!((rhr.scalar("T max").unwrap() - rh0.scalar("T max").unwrap()).abs() < 1e-10);
    assert!((rhr.scalar("T min").unwrap() - rh0.scalar("T min").unwrap()).abs() < 1e-10);

    // mixed_3d rigid Euclidean rotation + translation (rotating CurrentTerminal.direction!)
    let mut pm0 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
    pm0.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::LondonSheet { l_k_sq: 1.0e-12 } },
        Bc {
            labels: vec![1],
            kind: BcKind::CurrentTerminal { name: "t".into(), current: 1.0, direction: [0.0, 0.0, 1.0] },
        },
        Bc { labels: vec![2, 3, 4, 5, 6], kind: BcKind::Ground },
    ];
    let mut pmr = pm0.clone();
    let d_rot = [rmat[0][2], rmat[1][2], rmat[2][2]];
    pmr.bcs[1] =
        Bc { labels: vec![1], kind: BcKind::CurrentTerminal { name: "t".into(), current: 1.0, direction: d_rot } };
    let l_m0 =
        InductanceMatrixData::of(&solve_magnetostatics(&base3, &pm0).unwrap()).unwrap().inductance_matrix_h[0][0];
    let l_mr = InductanceMatrixData::of(&solve_magnetostatics(&rot3, &pmr).unwrap()).unwrap().inductance_matrix_h[0][0];
    assert!(
        (l_mr - l_m0).abs() / l_m0 < 1e-10,
        "Rigid Euclidean mixed_3d L err: l_m0={l_m0:e}, l_mr={l_mr:e}, rel={:e}",
        (l_mr - l_m0).abs() / l_m0
    );

    // 2. Geometry-level rotated remeshing (< 1e-3): square inner hole inside disk rotated by 37 degrees
    let solve_rot_geom = |ang_deg: f64| -> f64 {
        let th = ang_deg.to_radians();
        let (c, s) = (th.cos(), th.sin());
        let rot_pt = |x: f64, y: f64| (Expr::Num(c * x - s * y), Expr::Num(s * x + c * y));
        let mut m = Model::new_2d();
        let outer = m.circle("diel", (0.0, 0.0), 1.0);
        let hole = m.polygon("core", &[rot_pt(-0.3, -0.2), rot_pt(0.3, -0.2), rot_pt(0.3, 0.2), rot_pt(-0.3, 0.2)]);
        use geom::model::RegionKind;
        m.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", hole, RegionKind::Hole, None)], &[]);
        let mut opts = GeomMeshOpts::new(0.05);
        opts.max_dev = Some(0.001);
        let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
        let (mut inner_l, mut outer_l) = (Vec::new(), Vec::new());
        for f in 0..g.mesh.n_facets() {
            let fv = g.mesh.facet(f);
            let rm = (0.5 * (g.mesh.point(fv[0] as usize)[0] + g.mesh.point(fv[1] as usize)[0]))
                .hypot(0.5 * (g.mesh.point(fv[0] as usize)[1] + g.mesh.point(fv[1] as usize)[1]));
            let lab = g.mesh.facet_labels[f];
            if rm < 0.7 && !inner_l.contains(&lab) {
                inner_l.push(lab);
            } else if rm >= 0.7 && !outer_l.contains(&lab) {
                outer_l.push(lab);
            }
        }
        let mut p = Problem::new(Study::Electrostatic, 2);
        p.bcs =
            vec![Bc { labels: inner_l, kind: BcKind::Potential(1.0) }, Bc { labels: outer_l, kind: BcKind::Ground }];
        solve_electrostatics(&g.mesh, &p).unwrap().scalar("capacitance").unwrap()
    };
    let cg0 = solve_rot_geom(0.0);
    let cg37 = solve_rot_geom(37.0);
    assert!((cg37 - cg0).abs() / cg0 < 1e-3, "Geometry remesh 37deg C err {}", (cg37 - cg0).abs() / cg0);

    // 3. Vertex & cell permutation invariance (< 1e-10)
    let pmesh3 = permute_mesh(&base3);
    let rpe = solve_electrostatics(&pmesh3, &pe0).unwrap().scalar("capacitance").unwrap();
    let rph = solve_heat(&pmesh3, &ph0).unwrap().scalar("T max").unwrap();
    let rpm =
        InductanceMatrixData::of(&solve_magnetostatics(&pmesh3, &pm0).unwrap()).unwrap().inductance_matrix_h[0][0];
    assert!((rpe - c_e0).abs() / c_e0 < 1e-11);
    assert!((rph - rh0.scalar("T max").unwrap()).abs() < 1e-10);
    assert!((rpm - l_m0).abs() / l_m0 < 1e-11);

    // 4. Terminal permutation P M P^T (< 1e-12)
    let mut pt0 = Problem::new(Study::Electrostatic, 2);
    pt0.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("a".into()) },
        Bc { labels: vec![2], kind: BcKind::Terminal("b".into()) },
        Bc { labels: vec![3], kind: BcKind::Terminal("c".into()) },
        Bc { labels: vec![6], kind: BcKind::Ground },
    ];
    let mut ptp = pt0.clone();
    let perm = [2usize, 0, 1];
    ptp.bcs = vec![pt0.bcs[perm[0]].clone(), pt0.bcs[perm[1]].clone(), pt0.bcs[perm[2]].clone(), pt0.bcs[3].clone()];
    let m0 = CapacitanceMatrixData::of(&solve_electrostatics(&base3, &pt0).unwrap()).unwrap().maxwell_matrix_f;
    let mp = CapacitanceMatrixData::of(&solve_electrostatics(&base3, &ptp).unwrap()).unwrap().maxwell_matrix_f;
    for i in 0..3 {
        for j in 0..3 {
            let diff = (mp[i][j] - m0[perm[i]][perm[j]]).abs() / m0[perm[i]][perm[j]].abs().max(1e-30);
            assert!(diff < 1e-12, "Terminal permutation C mismatch ({i},{j}): {diff:e}");
        }
    }

    // 5. Dimensional scaling invariance (< 1e-11 in 3D: C(s*Omega) = s*C, L(s*Omega) = s*L when sheet_lk = 0)
    let scale = 1000.0;
    let smesh3 = scale_mesh(&base3, scale);
    let c_sc = solve_electrostatics(&smesh3, &pe0).unwrap().scalar("capacitance").unwrap();
    assert!((c_sc / (scale * c_e0) - 1.0).abs() < 1e-11, "3D Electro scaling err");

    let mut pm_geom = pm0.clone();
    pm_geom.bcs[0] = Bc { labels: vec![1], kind: BcKind::LondonSheet { l_k_sq: 0.0 } };
    let l_g0 =
        InductanceMatrixData::of(&solve_magnetostatics(&base3, &pm_geom).unwrap()).unwrap().inductance_matrix_h[0][0];
    let l_gs =
        InductanceMatrixData::of(&solve_magnetostatics(&smesh3, &pm_geom).unwrap()).unwrap().inductance_matrix_h[0][0];
    assert!((l_gs / (scale * l_g0) - 1.0).abs() < 1e-11, "3D mixed_3d scaling err");
}

#[test]
fn r32_curved_meshes_rejected_with_unsupported_in_heat_and_magneto() {
    let mut curved_3d = femcore::testkit::block(1.0, 1.0, 1.0, 2, 2, 2);
    let e0 = curved_3d.edges()[0];
    let p0 = curved_3d.point(e0.0 as usize);
    let p1 = curved_3d.point(e0.1 as usize);
    curved_3d.set_edge_nodes(&[(
        [e0.0, e0.1],
        [0.5 * (p0[0] + p1[0]) + 0.01, 0.5 * (p0[1] + p1[1]), 0.5 * (p0[2] + p1[2])],
    )]);
    assert!(curved_3d.has_curved());
    let ex_3d = Mesh::from_solver_mesh(&curved_3d).unwrap();

    let mut curved_2d = femcore::testkit::rect(1.0, 1.0, 2, 2);
    let e2 = curved_2d.edges()[0];
    let q0 = curved_2d.point(e2.0 as usize);
    let q1 = curved_2d.point(e2.1 as usize);
    let ex_2d = {
        let m = Mesh::from_solver_mesh(&curved_2d).unwrap();
        match m {
            Mesh::D2(mut cm) => {
                cm.set_edge_nodes([(
                    [e2.0.min(e2.1), e2.0.max(e2.1)],
                    [0.5 * (q0[0] + q1[0]) + 0.01, 0.5 * (q0[1] + q1[1])],
                )]);
                Mesh::D2(cm)
            }
            _ => unreachable!(),
        }
    };
    curved_2d.set_edge_nodes(&[([e2.0, e2.1], [0.5 * (q0[0] + q1[0]) + 0.01, 0.5 * (q0[1] + q1[1]), 0.0])]);

    let prog = Arc::new(Progress::new());

    // 1. Electrostatics accepts curved 3D Tet positively
    let mut p_es = Problem::new(Study::Electrostatic, 2);
    p_es.bcs = vec![Bc { labels: vec![5], kind: BcKind::Potential(1.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    assert!(run(&ex_3d, &p_es, &prog).is_ok());

    // 2. Heat (Stationary + Transient) rejects curved simplex via run() with Error::Unsupported & direct kernels
    let mut p_hs = Problem::new(Study::HeatStationary, 2);
    p_hs.bcs = vec![Bc { labels: vec![5], kind: BcKind::Temperature(300.0) }];
    assert!(matches!(run(&ex_3d, &p_hs, &prog), Err(Error::Unsupported(_))));
    assert!(matches!(run(&ex_2d, &p_hs, &prog), Err(Error::Unsupported(_))));
    match assemble_heat(&curved_3d, &p_hs) {
        Err(e) => assert!(e.contains("unsupported"), "expected unsupported error, got: {e}"),
        Ok(_) => panic!("assemble_heat must reject curved simplex mesh (R32)"),
    }
    assert!(solve_heat(&curved_3d, &p_hs).unwrap_err().contains("unsupported"));
    assert!(solve_heat_mixed(&ex_3d, &p_hs).unwrap_err().contains("unsupported"));

    let mut p_ht = Problem::new(Study::HeatTransient { t_end: 0.1, n_out: 2, t0: "300.0".into() }, 2);
    p_ht.bcs = vec![Bc { labels: vec![5], kind: BcKind::Temperature(300.0) }];
    assert!(matches!(run(&ex_3d, &p_ht, &prog), Err(Error::Unsupported(_))));

    // 3. Magnetostatics (cross_section, planar_stream, mixed_3d) rejects curved meshes via run() & direct kernels
    for (form, mesh_ex, mesh_d) in
        [("cross_section", &ex_2d, &curved_2d), ("planar_stream", &ex_2d, &curved_2d), ("mixed_3d", &ex_3d, &curved_3d)]
    {
        let mut pm = Problem::new(Study::Magnetostatic { formulation: form.into() }, 1);
        pm.bcs = vec![
            Bc { labels: vec![1], kind: BcKind::Ground },
            Bc {
                labels: vec![2],
                kind: BcKind::CurrentTerminal { name: "t".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
            },
        ];
        assert!(
            matches!(run(mesh_ex, &pm, &prog), Err(Error::Unsupported(_))),
            "{form} run() must reject with Unsupported"
        );
        assert!(solve_magnetostatics(mesh_d, &pm).unwrap_err().contains("unsupported"));
        assert!(solve_magnetostatics_mixed(mesh_ex, &pm).unwrap_err().contains("unsupported"));
    }
}
