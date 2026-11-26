// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Transient heat BDF2 temporal order and energy-balance verification:
//! * Verifies stage-1 and multi-step forcing inside `AmgSystem` BDF2 with opt-in `force_bdf2` bypass
//!   using `T*(x, z, t) = (1 + 0.5 sin(omega t) + 0.25 t^3) sin(pi x) sin(pi z)` across `sub = 4, 8, 16`.
//! * Verifies transient thermal energy balance converges at `O(dt^2)`.

use femcore::fem::{geometry, shape};
use femcore::props::LinearSolver;
use femcore::vvuq::{assert_convergence_rate, FnSource};
use nyx_dfs::heat::{solve_heat_bdf2_with_source, TransientHeatOptions};
use nyx_dfs::{Bc, BcKind, Problem, Study};
use std::f64::consts::PI;

fn rect_2d(w: f64, h: f64, nx: usize, nz: usize) -> femcore::mesh::SimplexMesh {
    femcore::testkit::rect(w, h, nx, nz)
}

fn quad_tri_7() -> (Vec<[f64; 3]>, Vec<f64>) {
    const T4: [f64; 4] =
        [-0.861_136_311_594_052_6, -0.339_981_043_584_856_3, 0.339_981_043_584_856_3, 0.861_136_311_594_052_6];
    const W4: [f64; 4] =
        [0.347_854_845_137_453_8, 0.652_145_154_862_546_1, 0.652_145_154_862_546_1, 0.347_854_845_137_453_8];
    let mut lam = Vec::with_capacity(16);
    let mut w = Vec::with_capacity(16);
    for (&tu, &wu) in T4.iter().zip(&W4) {
        let u = 0.5 * (1.0 + tu);
        let du = 0.5 * wu;
        for (&tv, &wv) in T4.iter().zip(&W4) {
            let v = 0.5 * (1.0 + tv);
            let dv = 0.5 * wv;
            let l1 = u;
            let l2 = (1.0 - u) * v;
            let l0 = 1.0 - l1 - l2;
            lam.push([l0, l1, l2]);
            w.push(2.0 * (1.0 - u) * du * dv);
        }
    }
    (lam, w)
}

fn l2_error_2d<F: Fn(f64, f64) -> f64>(mesh: &femcore::mesh::SimplexMesh, order: usize, u: &[f64], exact: F) -> f64 {
    let (qlam, qw) = quad_tri_7();
    let geo = geometry(mesh);
    let shapes: Vec<Vec<f64>> = qlam.iter().map(|l| shape(l, order).0).collect();
    let mut dofs = Vec::new();
    let mut err2 = 0.0_f64;
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let pts = [mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize)];
        mesh.cell_dofs(c, order, &mut dofs);
        for (q, lam) in qlam.iter().enumerate() {
            let x = lam[0] * pts[0][0] + lam[1] * pts[1][0] + lam[2] * pts[2][0];
            let z = lam[0] * pts[0][1] + lam[1] * pts[1][1] + lam[2] * pts[2][1];
            let uh: f64 = shapes[q].iter().zip(&dofs).map(|(&n, &d)| n * u[d]).sum();
            let diff = uh - exact(x, z);
            err2 += geo.vol[c] * qw[q] * diff * diff;
        }
    }
    err2.sqrt()
}

#[test]
fn heat_bdf2_amg_time_varying_forcing_achieves_second_order_convergence() {
    // T*(x, z, t) = (1 + 0.5 sin(omega t) + 0.25 t^3) sin(pi x) sin(pi z)
    // Spatial P2 on 32x32 grid gives O(h^3) ~ 7e-6 << BDF2 O(dt^2) errors at sub = 4, 8, 16.
    let mesh = rect_2d(1.0, 1.0, 32, 32);
    let k_cond = 0.2;
    let omega = 1.0;
    let t_end = 1.0;

    let g = |t: f64| 1.0 + 0.5 * (omega * t).sin() + 0.25 * t.powi(3);
    let gp = |t: f64| 0.5 * omega * (omega * t).cos() + 0.75 * t * t;
    let t_star = |x: f64, z: f64, t: f64| g(t) * (PI * x).sin() * (PI * z).sin();

    let source = FnSource {
        volume_rhs: |x: [f64; 2], t: f64, _dom: u32| -> [f64; 3] {
            let rhs = (gp(t) + k_cond * 2.0 * PI * PI * g(t)) * (PI * x[0]).sin() * (PI * x[1]).sin();
            [rhs, 0.0, 0.0]
        },
        dirichlet: |_x: [f64; 2], _t: f64, _bnd: u32| -> [f64; 3] { [0.0; 3] },
        neumann_robin: |_x: [f64; 2], _n: [f64; 2], _t: f64, _bnd: u32| -> [f64; 3] { [0.0; 3] },
    };

    let mut spec = Problem::new(Study::HeatTransient { t_end, n_out: 1, t0: "sin(pi*x)*sin(pi*z)".into() }, 2);
    spec.solver = LinearSolver::Iterative;
    spec.tol = 1e-12;
    spec.materials.default.insert("k".into(), k_cond);
    spec.materials.default.insert("rho".into(), 1.0);
    spec.materials.default.insert("Cp".into(), 1.0);
    spec.bcs = vec![Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Temperature(0.0) }];

    let subs = [4usize, 8, 16];
    let mut dts = Vec::new();
    let mut errs = Vec::new();
    for &sub in &subs {
        let res = solve_heat_bdf2_with_source::<2, _>(
            &mesh,
            &spec,
            &source,
            TransientHeatOptions { bdf2_substeps: sub, force_bdf2: true },
        )
        .unwrap();
        assert!(
            res.log.iter().any(|l| l.contains("AMG-CG")),
            "expected AmgSystem path when Iterative is selected, got {:?}",
            res.log
        );
        let t_final = &res.solutions.last().unwrap().field("T").unwrap().values;
        let err = l2_error_2d(&mesh, 2, t_final, |x, z| t_star(x, z, t_end));
        dts.push(t_end / sub as f64);
        errs.push(err);
    }
    let rep = assert_convergence_rate(&dts, &errs, 2.0, (1.90, 2.15));
    eprintln!(
        "heat_bdf2_order: asymptotic_slope = {:.4}, ols = {:.4}, dts = {dts:?}, errs = {errs:?}",
        rep.asymptotic_slope, rep.ols_slope
    );

    // 2. Smooth transient energy balance converges at O(dt^2):
    //    For homogeneous Dirichlet T = 0 on dOmega, testing rho Cp T_t - div(k grad T) = f by T gives
    //    0.5 * rho * Cp * (||T(t_end)||^2 - ||T(0)||^2) + int_0^{t_end} (int k |grad T|^2 dV) dt
    //      = int_0^{t_end} (int f * T dV) dt
    //    Equivalently on each time step with T* smooth (zero boundary-layer singularity), composite
    //    trapezoidal integration of the exact continuous energy balance converges at O(dt^2).
    let (qlam, qw) = quad_tri_7();
    let geo = geometry(&mesh);
    let shapes_v: Vec<Vec<f64>> = qlam.iter().map(|l| shape(l, 2).0).collect();
    let sym = nyx_dfs::heat::assemble_heat_with_source::<2, _>(&mesh, &spec, &source).unwrap();
    let m_mat = femcore::fem::assemble_mass(&mesh, 2, &vec![1.0; mesh.n_cells()]);
    let energy_residual_l2 = |n_out: usize| -> (f64, f64) {
        let mut s_eb = spec.clone();
        s_eb.study = Study::HeatTransient { t_end, n_out, t0: "sin(pi*x)*sin(pi*z)".into() };
        let res = solve_heat_bdf2_with_source::<2, _>(
            &mesh,
            &s_eb,
            &source,
            TransientHeatOptions { bdf2_substeps: 1, force_bdf2: true },
        )
        .unwrap();
        let dt = t_end / n_out as f64;
        let t0_v = &res.solutions[0].field("T").unwrap().values;
        let t1_v = &res.solutions[n_out].field("T").unwrap().values;
        let dot_mat = |u: &[f64], mat: &femcore::linalg::csr::RCsr| -> f64 {
            let mu = mat.matvec(u);
            u.iter().zip(&mu).map(|(a, b)| a * b).sum()
        };
        let d_store = 0.5 * (dot_mat(t1_v, &m_mat) - dot_mat(t0_v, &m_mat));
        let mut diss = 0.0_f64;
        let mut power = 0.0_f64;
        let mut dofs = Vec::new();
        for j in 0..=n_out {
            let tj = j as f64 * dt;
            let wt = if j == 0 || j == n_out { 0.5 * dt } else { dt };
            let tv = &res.solutions[j].field("T").unwrap().values;
            diss += wt * dot_mat(tv, &sym.k);
            // int_\Omega f(x, tj) * T_h(x, tj) dV
            let mut p_j = 0.0_f64;
            for c in 0..mesh.n_cells() {
                let v = mesh.cell(c);
                let pts = [mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize)];
                mesh.cell_dofs(c, 2, &mut dofs);
                for (q, lam) in qlam.iter().enumerate() {
                    let x = lam[0] * pts[0][0] + lam[1] * pts[1][0] + lam[2] * pts[2][0];
                    let z = lam[0] * pts[0][1] + lam[1] * pts[1][1] + lam[2] * pts[2][1];
                    let uh: f64 = shapes_v[q].iter().zip(&dofs).map(|(&n, &d)| n * tv[d]).sum();
                    let fval = (gp(tj) + k_cond * 2.0 * PI * PI * g(tj)) * (PI * x).sin() * (PI * z).sin();
                    p_j += geo.vol[c] * qw[q] * fval * uh;
                }
            }
            power += wt * p_j;
        }
        let residual = (d_store + diss - power).abs() / power.abs();
        (dt, residual)
    };
    let eb_steps = [8usize, 16, 32];
    let mut eb_dts = Vec::new();
    let mut eb_errs = Vec::new();
    for &n_out in &eb_steps {
        let (dt, r) = energy_residual_l2(n_out);
        eb_dts.push(dt);
        eb_errs.push(r);
    }
    let rep_eb = assert_convergence_rate(&eb_dts, &eb_errs, 2.0, (1.80, 2.25));
    eprintln!(
        "heat_bdf2 energy balance: asymptotic_slope = {:.4}, ols = {:.4}, errs = {eb_errs:?}",
        rep_eb.asymptotic_slope, rep_eb.ols_slope
    );
}
