// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME VVUQ Pillars (b) & (c) + Rough-Spot Closure Tests (`R6a`, `R6b`, `R6d`, `R32`):
//! - Zero-sized monomorphized `SourceField<D>` MMS in 2D (`Plane::Strain`, `Plane::Stress`) and 3D
//!   (isotropic, rotated cubic Si, Bunge Z-X-Z rotated sapphire + anisotropic thermal strain)
//! - Explicit `PropValue::Sym6` `[0,1,2,5,4,3]` ordering verification via `Sms::solve` (`validate.rs:144`)
//! - Native `MixedDiscretisation` MMS (`Tri` + `Quad`, curved annulus)
//! - Formal `femcore::vvuq::assert_convergence_rate` suite (P1/P2 `L^2`, `H^1`, symmetric eigenfrequency,
//!   and `axisym` PML `R6d` convergence rates)
//! - `R6a`: `u_r = i m u_theta` algebraic MPC on `r = 0` for `m = +-1` (elastic & piezo, simplex & mixed)
//! - `R6b`: axisymmetric `lrb::damping_matrix_axisym` radial weight `r` unit & energy-flux balance
//! - `R32`: explicit `Error::Unsupported` / refusal when curved or `qp_edges` meshes enter unsupported paths

mod common;

use femcore::fem::{geometry, shape, shape_gradients, RefElement};
use femcore::fem_mixed::{quad_rule_mixed, LagrangeSpace, MixedLagrangeDofMap};
use femcore::linalg::csr::{c, C};
use femcore::material::{stiffness_global, tensor_global, Anisotropy, Stiffness, Tensor, M3, M6};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use femcore::progress::Progress;
use femcore::props::{LinearSolver, Pml, PmlDir};
use femcore::section::plane_stiffness;
use femcore::vvuq::{assert_convergence_rate, FnSource};
use nyx_physics::{
    discretise_mixed, BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, PropValue, Scalar,
};
use nyx_sms::cross_section::{assemble_cross_section, solve_cross_section_modes, CrossSection};
use nyx_sms::lrb::damping_matrix_axisym;
use nyx_sms::solid::{setup_with_source, solve_solid_mixed_with_source, solve_solid_with_source};
use nyx_sms::*;
use serde_json::json;
use std::f64::consts::PI;

// ---- Analytical MMS helpers -------------------------------------------------------------------------

#[inline]
fn sin_cos_pi(x: f64) -> (f64, f64) {
    (PI * x).sin_cos()
}

/// Evaluates 2D `u*([x, z]) = [sin(pi x) cos(pi z), cos(pi x) sin(pi z)]` and its gradient `du[comp][axis]`.
fn mms_u_2d(p: [f64; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
    let (sx, cx) = sin_cos_pi(p[0]);
    let (sz, cz) = sin_cos_pi(p[1]);
    let u = [sx * cz, cx * sz];
    let du = [[PI * cx * cz, -PI * sx * sz], [-PI * sx * sz, PI * cx * cz]];
    (u, du)
}

/// Evaluates `-div(sigma(u*, dT))` and `sigma(u*, dT) * n` for arbitrary 2D constitutive plane matrix `C_plane` and `alpha_plane`.
fn mms_forcing_2d(
    p: [f64; 2],
    c_plane: &M6,
    plane: Plane,
    alpha_plane: &M3,
    dt_fn: &dyn Fn([f64; 2]) -> (f64, [f64; 2]),
) -> ([f64; 2], [[f64; 2]; 2]) {
    let (sx, cx) = sin_cos_pi(p[0]);
    let (sz, cz) = sin_cos_pi(p[1]);
    let p2 = PI * PI;

    // d^2 u_i / dx_a dx_b: symmetric 2x2 per comp -> (xx, xz, zz)
    let d2_u0 = [-p2 * sx * cz, -p2 * cx * sz, -p2 * sx * cz];
    let d2_u1 = [-p2 * cx * sz, -p2 * sx * cz, -p2 * cx * sz];

    // eps = [du0/dx, du1/dz, du0/dz + du1/dx]
    let (_, du) = mms_u_2d(p);
    let eps = [du[0][0], du[1][1], du[0][1] + du[1][0]];
    // deps / dx = [d2_u0/dxx, d2_u1/dxz, d2_u0/dxz + d2_u1/dxx]
    // deps / dz = [d2_u0/dxz, d2_u1/dzz, d2_u0/dzz + d2_u1/dxz]
    let d_eps = [[d2_u0[0], d2_u1[1], d2_u0[1] + d2_u1[0]], [d2_u0[1], d2_u1[2], d2_u0[2] + d2_u1[1]]];

    const INP: [usize; 3] = [0, 1, 3];
    let mut d_mat = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            d_mat[a][b] = match plane {
                Plane::Strain => c_plane[INP[a]][INP[b]],
                Plane::Stress => c_plane[INP[a]][INP[b]] - c_plane[INP[a]][2] * c_plane[2][INP[b]] / c_plane[2][2],
            };
        }
    }

    let (dt, d_dt) = dt_fn(p);
    let eth6 = |v: f64| -> [f64; 6] {
        [
            alpha_plane[0][0] * v,
            alpha_plane[1][1] * v,
            alpha_plane[2][2] * v,
            2.0 * alpha_plane[0][1] * v,
            2.0 * alpha_plane[1][2] * v,
            2.0 * alpha_plane[0][2] * v,
        ]
    };
    let th_stress = |v: f64| -> [f64; 3] {
        let e = eth6(v);
        let mut s = [0.0; 3];
        match plane {
            Plane::Strain => {
                for a in 0..3 {
                    s[a] = (0..6).map(|j| c_plane[INP[a]][j] * e[j]).sum();
                }
            }
            Plane::Stress => {
                for a in 0..3 {
                    s[a] = (0..3).map(|b| d_mat[a][b] * e[INP[b]]).sum();
                }
            }
        }
        s
    };

    let s_th = th_stress(dt);
    let ds_th = [th_stress(d_dt[0]), th_stress(d_dt[1])];

    let mut sig = [0.0; 3];
    let mut dsig = [[0.0; 3]; 2];
    for a in 0..3 {
        sig[a] = (0..3).map(|b| d_mat[a][b] * eps[b]).sum::<f64>() - s_th[a];
        for k in 0..2 {
            dsig[k][a] = (0..3).map(|b| d_mat[a][b] * d_eps[k][b]).sum::<f64>() - ds_th[k][a];
        }
    }
    let f = [-(dsig[0][0] + dsig[1][2]), -(dsig[0][2] + dsig[1][1])];
    let sig_tensor = [[sig[0], sig[2]], [sig[2], sig[1]]];
    (f, sig_tensor)
}

/// Evaluates 3D `u* = [sin(pi x)cos(pi y)cos(pi z), cos(pi x)sin(pi y)cos(pi z), cos(pi x)cos(pi y)sin(pi z)]`
/// and its gradient `du[comp][axis]`.
fn mms_u_3d(p: [f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let (sx, cx) = sin_cos_pi(p[0]);
    let (sy, cy) = sin_cos_pi(p[1]);
    let (sz, cz) = sin_cos_pi(p[2]);
    let u = [sx * cy * cz, cx * sy * cz, cx * cy * sz];
    let du = [
        [PI * cx * cy * cz, -PI * sx * sy * cz, -PI * sx * cy * sz],
        [-PI * sx * sy * cz, PI * cx * cy * cz, -PI * cx * sy * sz],
        [-PI * sx * cy * sz, -PI * cx * sy * sz, PI * cx * cy * cz],
    ];
    (u, du)
}

/// Evaluates `-div(sigma(u*, dT))` and `sigma(u*, dT)` tensor in 3D solver Voigt order `(xx, yy, zz, xy, yz, xz)`.
fn mms_forcing_3d(
    p: [f64; 3],
    c6: &M6,
    alpha: &M3,
    dt_fn: &dyn Fn([f64; 3]) -> (f64, [f64; 3]),
) -> ([f64; 3], [[f64; 3]; 3]) {
    let (sx, cx) = sin_cos_pi(p[0]);
    let (sy, cy) = sin_cos_pi(p[1]);
    let (sz, cz) = sin_cos_pi(p[2]);
    let p2 = PI * PI;
    let (u, du) = mms_u_3d(p);

    // Partial derivatives of du[comp][j] wrt x_m: d2u[comp][j][m]
    let mut d2 = [[[0.0; 3]; 3]; 3];
    // comp 0: u0 = sx * cy * cz
    d2[0][0] = [-p2 * u[0], -p2 * cx * sy * cz, -p2 * cx * cy * sz];
    d2[0][1] = [-p2 * cx * sy * cz, -p2 * u[0], p2 * sx * sy * sz];
    d2[0][2] = [-p2 * cx * cy * sz, p2 * sx * sy * sz, -p2 * u[0]];
    // comp 1: u1 = cx * sy * cz
    d2[1][0] = [-p2 * u[1], -p2 * sx * cy * cz, p2 * sx * sy * sz];
    d2[1][1] = [-p2 * sx * cy * cz, -p2 * u[1], -p2 * cx * cy * sz];
    d2[1][2] = [p2 * sx * sy * sz, -p2 * cx * cy * sz, -p2 * u[1]];
    // comp 2: u2 = cx * cy * sz
    d2[2][0] = [-p2 * u[2], p2 * sx * sy * sz, -p2 * sx * cy * cz];
    d2[2][1] = [p2 * sx * sy * sz, -p2 * u[2], -p2 * cx * sy * cz];
    d2[2][2] = [-p2 * sx * cy * cz, -p2 * cx * sy * cz, -p2 * u[2]];

    let eps = [du[0][0], du[1][1], du[2][2], du[0][1] + du[1][0], du[1][2] + du[2][1], du[0][2] + du[2][0]];
    let mut d_eps = [[0.0; 6]; 3];
    for m in 0..3 {
        d_eps[m] = [
            d2[0][0][m],
            d2[1][1][m],
            d2[2][2][m],
            d2[0][1][m] + d2[1][0][m],
            d2[1][2][m] + d2[2][1][m],
            d2[0][2][m] + d2[2][0][m],
        ];
    }
    let eth6 = |v: f64| -> [f64; 6] {
        [
            alpha[0][0] * v,
            alpha[1][1] * v,
            alpha[2][2] * v,
            2.0 * alpha[0][1] * v,
            2.0 * alpha[1][2] * v,
            2.0 * alpha[0][2] * v,
        ]
    };
    let (dt, d_dt) = dt_fn(p);
    let e_th = eth6(dt);
    let de_th = [eth6(d_dt[0]), eth6(d_dt[1]), eth6(d_dt[2])];

    let mut sig = [0.0; 6];
    let mut dsig = [[0.0; 6]; 3];
    for i in 0..6 {
        sig[i] = (0..6).map(|j| c6[i][j] * (eps[j] - e_th[j])).sum();
        for m in 0..3 {
            dsig[m][i] = (0..6).map(|j| c6[i][j] * (d_eps[m][j] - de_th[m][j])).sum();
        }
    }
    // Solver Voigt: 0:xx, 1:yy, 2:zz, 3:xy, 4:yz, 5:xz
    let f = [
        -(dsig[0][0] + dsig[1][3] + dsig[2][5]),
        -(dsig[0][3] + dsig[1][1] + dsig[2][4]),
        -(dsig[0][5] + dsig[1][4] + dsig[2][2]),
    ];
    let s3 = [[sig[0], sig[3], sig[5]], [sig[3], sig[1], sig[4]], [sig[5], sig[4], sig[2]]];
    (f, s3)
}

/// Evaluates continuous `(L^2, H^1_semi)` errors on a `SimplexMesh` using degree-6 quadrature.
fn simplex_errors<const D: usize>(
    mesh: &SimplexMesh,
    order: usize,
    res: &nyx_physics::SolveResult,
    exact: &dyn Fn([f64; D]) -> ([f64; D], [[f64; D]; D]),
) -> (f64, f64) {
    let uf = res.solutions[0].field("u").expect("displacement field u");
    let u_comps: Vec<Vec<f64>> = (0..D).map(|k| uf.component(k)).collect();
    let geo = geometry(mesh);
    let (ct, ref_inv) = if D == 2 { (CellType::Tri, 2.0) } else { (CellType::Tet, 6.0) };
    let rule = quad_rule_mixed(ct, 6);
    let nloc = mesh.nloc(order);
    let mut re = RefElement::new(D, order, 1);
    let lams: Vec<Vec<f64>> = rule
        .pts
        .iter()
        .map(|xi| {
            if D == 2 {
                vec![1.0 - xi[0] - xi[1], xi[0], xi[1]]
            } else {
                vec![1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]]
            }
        })
        .collect();
    re.w = rule.w.iter().map(|&w| w * ref_inv).collect();
    re.n.clear();
    re.dn.clear();
    for lam in &lams {
        let (nv, dnv) = shape(lam, order);
        re.n.push(nv);
        re.dn.push(dnv);
    }

    let mut l2_sq = 0.0;
    let mut h1_sq = 0.0;
    let mut dofs = Vec::new();
    let mut grad = vec![0.0; nloc * D];
    for ci in 0..mesh.n_cells() {
        let cv = mesh.cell(ci);
        mesh.cell_dofs(ci, order, &mut dofs);
        for (q, lam) in lams.iter().enumerate() {
            let mut x_q = [0.0; D];
            for (i, &v) in cv.iter().enumerate() {
                let pt = mesh.point(v as usize);
                for k in 0..D {
                    x_q[k] += lam[i] * pt[k];
                }
            }
            shape_gradients(&re, &geo, ci, q, &mut grad);
            let (u_ex, du_ex) = exact(x_q);
            let mut u_h = [0.0; D];
            let mut du_h = [[0.0; D]; D];
            for a in 0..nloc {
                let na = re.n[q][a];
                for k in 0..D {
                    let uak = u_comps[k][dofs[a]];
                    u_h[k] += na * uak;
                    for j in 0..D {
                        du_h[k][j] += grad[a * D + j] * uak;
                    }
                }
            }
            let w_phys = re.w[q] * geo.vol[ci];
            for k in 0..D {
                l2_sq += w_phys * (u_h[k] - u_ex[k]).powi(2);
                for j in 0..D {
                    h1_sq += w_phys * (du_h[k][j] - du_ex[k][j]).powi(2);
                }
            }
        }
    }
    (l2_sq.sqrt(), h1_sq.sqrt())
}

// ---- Pillar (b) + Pillar (c) 2D MMS & convergence ---------------------------------------------------

fn run_2d_mms(
    n: usize,
    order: usize,
    plane: Plane,
    c_3d: M6,
    euler_deg: [f64; 3],
    aniso: Option<Anisotropy>,
) -> (f64, (f64, f64)) {
    let mesh = common::rect(1.0, 1.0, n, n);
    let c_plane = plane_stiffness(&c_3d, 2);
    let alpha_zero = [[0.0; 3]; 3];
    let mut sp = Problem::new(Study::Static, order);
    sp.plane = plane;
    if let Some(an) = aniso {
        sp.materials.aniso.insert(1, an);
    } else {
        sp.materials.default.insert("E".into(), 100.0e9);
        sp.materials.default.insert("nu".into(), 0.28);
    }
    // Mix inhomogeneous Dirichlet (x=0, z=0 => labels 1, 2) and Neumann traction (x=1, z=1 => labels 3, 4).
    let _ = euler_deg;
    sp.bcs = vec![
        Bc { labels: vec![1, 2], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } },
        Bc { labels: vec![3, 4], kind: BcKind::Load([0.0; 3]) },
    ];
    let src = FnSource {
        volume_rhs: |x: [f64; 2], _t, _d| {
            let (f, _) = mms_forcing_2d(x, &c_plane, plane, &alpha_zero, &|_| (0.0, [0.0, 0.0]));
            [f[0], f[1], 0.0]
        },
        dirichlet: |x: [f64; 2], _t, _b| {
            let (u, _) = mms_u_2d(x);
            [u[0], u[1], 0.0]
        },
        neumann_robin: |x: [f64; 2], n_hat: [f64; 2], _t, _b| {
            let (_, sig) = mms_forcing_2d(x, &c_plane, plane, &alpha_zero, &|_| (0.0, [0.0, 0.0]));
            [sig[0][0] * n_hat[0] + sig[0][1] * n_hat[1], sig[1][0] * n_hat[0] + sig[1][1] * n_hat[1], 0.0]
        },
    };
    let res = solve_solid_with_source(&mesh, &sp, &src).unwrap();
    let err = simplex_errors(&mesh, order, &res, &mms_u_2d);
    (1.0 / n as f64, err)
}

#[test]
fn mms_2d_plane_strain_and_stress_p1_p2_convergence_rates_and_native_quad() {
    let sap_st =
        Stiffness::Trigonal { c11: 496.8e9, c12: 163.6e9, c13: 110.9e9, c14: -23.5e9, c33: 498.1e9, c44: 147.4e9 };
    let euler_2d = [33.0, 57.6, -19.0];
    let sap_c6 = stiffness_global(&sap_st, euler_2d);
    let sap_an = Anisotropy { stiffness: Some(sap_st), euler_deg: euler_2d, ..Default::default() };

    let grids = [4usize, 8, 16];
    let hs: Vec<f64> = grids.iter().map(|&n| 1.0 / n as f64).collect();

    // 1. 2D Plane Strain rotated sapphire P1 and P2
    for (order, l2_band, h1_band, p_l2, p_h1) in
        [(1usize, (1.85, 2.25), (0.90, 1.15), 2.0, 1.0), (2usize, (2.80, 3.25), (1.85, 2.20), 3.0, 2.0)]
    {
        let mut l2s = Vec::new();
        let mut h1s = Vec::new();
        for &n in &grids {
            let (_, (l2, h1)) = run_2d_mms(n, order, Plane::Strain, sap_c6, euler_2d, Some(sap_an.clone()));
            l2s.push(l2);
            h1s.push(h1);
        }
        let r_l2 = assert_convergence_rate(&hs, &l2s, p_l2, l2_band);
        let r_h1 = assert_convergence_rate(&hs, &h1s, p_h1, h1_band);
        println!("2D PlaneStrain sapphire P{order}: L2={r_l2:?}, H1={r_h1:?}");
    }

    // 2. 2D Plane Stress (isotropic / orthotropic-aligned) P1 and P2
    let si_st = Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 };
    let si_c6 = stiffness_global(&si_st, [0.0; 3]);
    let si_an = Anisotropy { stiffness: Some(si_st), euler_deg: [0.0; 3], ..Default::default() };
    for (order, l2_band, h1_band, p_l2, p_h1) in
        [(1usize, (1.85, 2.25), (0.90, 1.15), 2.0, 1.0), (2usize, (2.80, 3.25), (1.85, 2.20), 3.0, 2.0)]
    {
        let mut l2s = Vec::new();
        let mut h1s = Vec::new();
        for &n in &grids {
            let (_, (l2, h1)) = run_2d_mms(n, order, Plane::Stress, si_c6, [0.0; 3], Some(si_an.clone()));
            l2s.push(l2);
            h1s.push(h1);
        }
        let r_l2 = assert_convergence_rate(&hs, &l2s, p_l2, l2_band);
        let r_h1 = assert_convergence_rate(&hs, &h1s, p_h1, h1_band);
        println!("2D PlaneStress cubic Si P{order}: L2={r_l2:?}, H1={r_h1:?}");
    }

    // 3. Also verify solve_solid_mixed_with_source on a native 2D Quad mesh [0,1]^2 at P2!
    let mut cm = CellMesh::<2>::new();
    let nq = 16usize;
    for j in 0..=nq {
        for i in 0..=nq {
            cm.points.push([i as f64 / nq as f64, j as f64 / nq as f64]);
        }
    }
    let vid = |i: usize, j: usize| (j * (nq + 1) + i) as u32;
    for j in 0..nq {
        for i in 0..nq {
            cm.push_cell(CellType::Quad, &[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)], 1);
        }
    }
    for j in 0..nq {
        cm.push_facet(&[vid(0, j), vid(0, j + 1)], 1);
        cm.push_facet(&[vid(nq, j), vid(nq, j + 1)], 4);
    }
    for i in 0..nq {
        cm.push_facet(&[vid(i, 0), vid(i + 1, 0)], 2);
        cm.push_facet(&[vid(i, nq), vid(i + 1, nq)], 3);
    }
    let raw_q = Mesh::D2(cm);
    let mixed_q = discretise_mixed(&raw_q).unwrap();
    let mut sp_q = Problem::new(Study::Static, 2);
    sp_q.plane = Plane::Strain;
    sp_q.materials.aniso.insert(1, sap_an);
    sp_q.bcs = vec![
        Bc { labels: vec![1, 2], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } },
        Bc { labels: vec![3, 4], kind: BcKind::Load([0.0; 3]) },
    ];
    let alpha_zero = [[0.0; 3]; 3];
    let src_q = FnSource {
        volume_rhs: |x: [f64; 2], _t, _d| {
            let (f, _) =
                mms_forcing_2d(x, &plane_stiffness(&sap_c6, 2), Plane::Strain, &alpha_zero, &|_| (0.0, [0.0, 0.0]));
            [f[0], f[1], 0.0]
        },
        dirichlet: |x: [f64; 2], _t, _b| {
            let (u, _) = mms_u_2d(x);
            [u[0], u[1], 0.0]
        },
        neumann_robin: |x: [f64; 2], n_hat: [f64; 2], _t, _b| {
            let (_, sig) =
                mms_forcing_2d(x, &plane_stiffness(&sap_c6, 2), Plane::Strain, &alpha_zero, &|_| (0.0, [0.0, 0.0]));
            [sig[0][0] * n_hat[0] + sig[0][1] * n_hat[1], sig[1][0] * n_hat[0] + sig[1][1] * n_hat[1], 0.0]
        },
    };
    let r_q = solve_solid_mixed_with_source(&raw_q, &mixed_q, &sp_q, &src_q).unwrap();
    let dm_q = MixedLagrangeDofMap::new(&mixed_q, LagrangeSpace::Order2Serendipity);
    let qcoords = dm_q.dof_coordinates(&mixed_q);
    let ux_proj = r_q.solutions[0].field("u").unwrap().component(0);
    let uz_proj = r_q.solutions[0].field("u").unwrap().component(1);
    let simplex_q = nyx_physics::discretise(&raw_q).unwrap();
    let sc_q = simplex_q.dof_coordinates(2);
    let max_node_err = (0..simplex_q.n_dofs(2))
        .map(|i| {
            let (ue, _) = mms_u_2d([sc_q[2 * i], sc_q[2 * i + 1]]);
            (ux_proj[i] - ue[0]).hypot(uz_proj[i] - ue[1])
        })
        .fold(0.0_f64, f64::max);
    let _ = qcoords;
    assert!(max_node_err < 1e-4, "native Quad2 P2 MMS max error {max_node_err:.3e}");
}

// ---- Pillar (b) + Pillar (c) 3D MMS (isotropic + Bunge rotated Si/sapphire + thermal + Sym6) ---------

fn kuhn_cube(n: usize) -> SimplexMesh {
    let mut pts = Vec::with_capacity(3 * (n + 1).pow(3));
    for k in 0..=n {
        for j in 0..=n {
            for i in 0..=n {
                pts.extend_from_slice(&[i as f64 / n as f64, j as f64 / n as f64, k as f64 / n as f64]);
            }
        }
    }
    let vid = |i: usize, j: usize, k: usize| ((k * (n + 1) + j) * (n + 1) + i) as u32;
    let mut cells = Vec::with_capacity(6 * 4 * n.pow(3));
    let mut clbls = Vec::with_capacity(6 * n.pow(3));
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let c = [
                    vid(i, j, k),
                    vid(i + 1, j, k),
                    vid(i + 1, j + 1, k),
                    vid(i, j + 1, k),
                    vid(i, j, k + 1),
                    vid(i + 1, j, k + 1),
                    vid(i + 1, j + 1, k + 1),
                    vid(i, j + 1, k + 1),
                ];
                let tets = [
                    [c[0], c[1], c[2], c[6]],
                    [c[0], c[2], c[3], c[6]],
                    [c[0], c[3], c[7], c[6]],
                    [c[0], c[7], c[4], c[6]],
                    [c[0], c[4], c[5], c[6]],
                    [c[0], c[5], c[1], c[6]],
                ];
                for t in tets {
                    cells.extend_from_slice(&t);
                    clbls.push(1);
                }
            }
        }
    }
    let mut facets = Vec::new();
    let mut flbls = Vec::new();
    for a in 0..n {
        for b in 0..n {
            facets.extend_from_slice(&[
                vid(0, a, b),
                vid(0, a + 1, b),
                vid(0, a + 1, b + 1),
                vid(0, a, b),
                vid(0, a + 1, b + 1),
                vid(0, a, b + 1),
            ]);
            flbls.extend_from_slice(&[1, 1]);
            facets.extend_from_slice(&[
                vid(n, a, b),
                vid(n, a + 1, b),
                vid(n, a + 1, b + 1),
                vid(n, a, b),
                vid(n, a + 1, b + 1),
                vid(n, a, b + 1),
            ]);
            flbls.extend_from_slice(&[6, 6]);
            facets.extend_from_slice(&[
                vid(a, 0, b),
                vid(a + 1, 0, b),
                vid(a + 1, 0, b + 1),
                vid(a, 0, b),
                vid(a + 1, 0, b + 1),
                vid(a, 0, b + 1),
            ]);
            flbls.extend_from_slice(&[2, 2]);
            facets.extend_from_slice(&[
                vid(a, n, b),
                vid(a + 1, n, b),
                vid(a + 1, n, b + 1),
                vid(a, n, b),
                vid(a + 1, n, b + 1),
                vid(a, n, b + 1),
            ]);
            flbls.extend_from_slice(&[5, 5]);
            facets.extend_from_slice(&[
                vid(a, b, 0),
                vid(a + 1, b, 0),
                vid(a + 1, b + 1, 0),
                vid(a, b, 0),
                vid(a + 1, b + 1, 0),
                vid(a, b + 1, 0),
            ]);
            flbls.extend_from_slice(&[3, 3]);
            facets.extend_from_slice(&[
                vid(a, b, n),
                vid(a + 1, b, n),
                vid(a + 1, b + 1, n),
                vid(a, b, n),
                vid(a + 1, b + 1, n),
                vid(a, b + 1, n),
            ]);
            flbls.extend_from_slice(&[4, 4]);
        }
    }
    SimplexMesh::new(3, pts, cells, clbls, facets, flbls)
}

#[test]
fn mms_3d_isotropic_bunge_rotated_sapphire_si_thermal_and_sym6_convergence() {
    let sap_st =
        Stiffness::Trigonal { c11: 496.8e9, c12: 163.6e9, c13: 110.9e9, c14: -23.5e9, c33: 498.1e9, c44: 147.4e9 };
    let euler_3d = [37.0, 57.6, -22.0];
    let sap_c6 = stiffness_global(&sap_st, euler_3d);
    let alpha_tensor = Tensor::Uniaxial { perp: 5.0e-6, par: 6.2e-6 };
    let alpha_global = tensor_global(&alpha_tensor, euler_3d);
    let sap_an = Anisotropy {
        stiffness: Some(sap_st),
        expansion: Some(alpha_tensor),
        euler_deg: euler_3d,
        ..Default::default()
    };

    // Spatially varying thermal field dT(x,y,z) = 120.0 + 40.0 * x - 30.0 * y^2 + 25.0 * z
    let dt_fn = |p: [f64; 3]| -> (f64, [f64; 3]) {
        (120.0 + 40.0 * p[0] - 30.0 * p[1] * p[1] + 25.0 * p[2], [40.0, -60.0 * p[1], 25.0])
    };
    let t_ref = 293.15;

    let grids = [4usize, 6, 8];
    let hs: Vec<f64> = grids.iter().map(|&n| 1.0 / n as f64).collect();

    for (order, l2_band, h1_band, p_l2, p_h1) in
        [(1usize, (1.85, 2.25), (0.90, 1.15), 2.0, 1.0), (2usize, (2.80, 3.25), (1.85, 2.20), 3.0, 2.0)]
    {
        let mut l2s = Vec::new();
        let mut h1s = Vec::new();
        for &n in &grids {
            let mesh = kuhn_cube(n);
            let coords = mesh.dof_coordinates(order);
            let t_imported: Vec<f64> = (0..mesh.n_dofs(order))
                .map(|i| t_ref + dt_fn([coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]).0)
                .collect();

            let mut sp = Problem::new(Study::Static, order);
            sp.solver = LinearSolver::Direct;
            sp.materials.aniso.insert(1, sap_an.clone());
            sp.thermal = Some(ThermalExpansion { source: ThermalSource::Imported(t_imported), t_ref });
            sp.bcs = vec![
                Bc { labels: vec![1, 2, 3], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
                Bc { labels: vec![4, 5, 6], kind: BcKind::Load([0.0; 3]) },
            ];
            let src = FnSource {
                volume_rhs: |x: [f64; 3], _t, _d| mms_forcing_3d(x, &sap_c6, &alpha_global, &dt_fn).0,
                dirichlet: |x: [f64; 3], _t, _b| mms_u_3d(x).0,
                neumann_robin: |x: [f64; 3], n_hat: [f64; 3], _t, _b| {
                    let s = mms_forcing_3d(x, &sap_c6, &alpha_global, &dt_fn).1;
                    [
                        s[0][0] * n_hat[0] + s[0][1] * n_hat[1] + s[0][2] * n_hat[2],
                        s[1][0] * n_hat[0] + s[1][1] * n_hat[1] + s[1][2] * n_hat[2],
                        s[2][0] * n_hat[0] + s[2][1] * n_hat[1] + s[2][2] * n_hat[2],
                    ]
                },
            };
            let res = solve_solid_with_source(&mesh, &sp, &src).unwrap();
            let (l2, h1) = simplex_errors(&mesh, order, &res, &mms_u_3d);
            l2s.push(l2);
            h1s.push(h1);
        }
        let r_l2 = assert_convergence_rate(&hs, &l2s, p_l2, l2_band);
        let r_h1 = assert_convergence_rate(&hs, &h1s, p_h1, h1_band);
        println!("3D Bunge-rotated Sapphire + thermal P{order}: L2={r_l2:?}, H1={r_h1:?}");
    }

    // Verify rotated cubic Si & isotropic 3D P2 convergence rates (O(h^3) L2, O(h^2) H1):
    let si_st = Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 };
    let si_c6 = stiffness_global(&si_st, [19.0, 31.0, -47.0]);
    let iso_c6 = femcore::material::isotropic_stiffness(165.0e9, 0.25);
    let zero_a = [[0.0; 3]; 3];
    let zero_dt = |_: [f64; 3]| (0.0, [0.0; 3]);
    for (tag, c6_ref, aniso_opt) in [
        ("3D isotropic", iso_c6, None),
        (
            "3D rotated Si",
            si_c6,
            Some(Anisotropy { stiffness: Some(si_st), euler_deg: [19.0, 31.0, -47.0], ..Default::default() }),
        ),
    ] {
        let mut l2s = Vec::new();
        let mut h1s = Vec::new();
        for &n in &grids {
            let m = kuhn_cube(n);
            let mut sp = Problem::new(Study::Static, 2);
            sp.solver = LinearSolver::Direct;
            if let Some(an) = &aniso_opt {
                sp.materials.aniso.insert(1, an.clone());
            } else {
                sp.materials.set(1, "E", 165.0e9);
                sp.materials.set(1, "nu", 0.25);
            }
            sp.bcs = vec![
                Bc { labels: vec![1, 2, 3], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
                Bc { labels: vec![4, 5, 6], kind: BcKind::Load([0.0; 3]) },
            ];
            let src = FnSource {
                volume_rhs: |x: [f64; 3], _t, _d| mms_forcing_3d(x, &c6_ref, &zero_a, &zero_dt).0,
                dirichlet: |x: [f64; 3], _t, _b| mms_u_3d(x).0,
                neumann_robin: |x: [f64; 3], n_hat: [f64; 3], _t, _b| {
                    let s = mms_forcing_3d(x, &c6_ref, &zero_a, &zero_dt).1;
                    [
                        s[0][0] * n_hat[0] + s[0][1] * n_hat[1] + s[0][2] * n_hat[2],
                        s[1][0] * n_hat[0] + s[1][1] * n_hat[1] + s[1][2] * n_hat[2],
                        s[2][0] * n_hat[0] + s[2][1] * n_hat[1] + s[2][2] * n_hat[2],
                    ]
                },
            };
            let r = solve_solid_with_source(&m, &sp, &src).unwrap();
            let (l2, h1) = simplex_errors(&m, 2, &r, &mms_u_3d);
            l2s.push(l2);
            h1s.push(h1);
        }
        let r_l2 = assert_convergence_rate(&hs, &l2s, 3.0, (2.80, 3.25));
        let r_h1 = assert_convergence_rate(&hs, &h1s, 2.0, (1.85, 2.20));
        assert!(l2s[2] < 1.5e-3, "{tag} n=8 P2 L2={:.3e}", l2s[2]);
        println!("{tag} P2: L2={r_l2:?}, H1={r_h1:?}, n=8 err={:.3e}", l2s[2]);
    }

    // Explicitly verify `PropValue::Sym6` [xx, yy, zz, xy, xz, yz] ordering via `Sms::solve` (`validate.rs:144`)!
    // Rotated thermal expansion tensor alpha_global applied through ProblemSpec matches Tensor::Symmetric identically.
    struct BoxResolver;
    impl LabelResolver for BoxResolver {
        fn domains(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "body" => Ok(vec![1]),
                _ => Err(Error::invalid(format!("unknown domain {name}"))),
            }
        }
        fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "xmin" => Ok(vec![1]),
                "ymin" => Ok(vec![2]),
                "zmin" => Ok(vec![3]),
                _ => Err(Error::invalid(format!("unknown boundary {name}"))),
            }
        }
    }
    let test_mesh = common::block(1.0, 0.5, 0.5, 3, 2, 2);
    let raw_ex = Mesh::from_solver_mesh(&test_mesh).unwrap();
    let sym6_vec = [
        alpha_global[0][0],
        alpha_global[1][1],
        alpha_global[2][2],
        alpha_global[0][1],
        alpha_global[0][2],
        alpha_global[1][2],
    ];
    let mut spec = ProblemSpec::new("sms", 2, "static");
    spec.materials.insert(
        "body".into(),
        MaterialAssign::PropsTensor(
            [
                ("E".to_string(), PropValue::Scalar(200.0e9)),
                ("nu".to_string(), PropValue::Scalar(0.28)),
                ("alpha".to_string(), PropValue::Sym6(sym6_vec)),
                // Exercise `validate.rs:154` Sym6 `[0,1,2,5,4,3]` permutation check on `loss_tangent` too!
                ("loss_tangent".to_string(), PropValue::Sym6([0.02, 0.03, 0.04, 0.001, -0.002, 0.003])),
            ]
            .into(),
        ),
    );
    spec.bcs = vec![
        BcSpec { on: "xmin".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
        BcSpec { on: "ymin".into(), kind: "fixed".into(), params: json!({ "uy": 0.0 }) },
        BcSpec { on: "zmin".into(), kind: "fixed".into(), params: json!({ "uz": 0.0 }) },
    ];
    spec.params = json!({ "thermal": "uniform", "temperature": 393.15, "t_ref": 293.15 });
    let res_sym6 = Sms.solve(&raw_ex, &spec, &mut Ctx::new(&BoxResolver)).unwrap();
    // Verify against typed `Tensor::Symmetric { xx, yy, zz, xy, yz, xz }` on identical BCs (`< 1e-13`),
    // plus confirm non-degeneracy when shear entries (xy, yz) are swapped:
    let mut p_ref_sym = Problem::new(Study::Static, 2);
    p_ref_sym.materials.set(1, "E", 200.0e9);
    p_ref_sym.materials.set(1, "nu", 0.28);
    p_ref_sym.materials.aniso.insert(
        1,
        Anisotropy {
            expansion: Some(femcore::material::Tensor::Symmetric {
                xx: alpha_global[0][0],
                yy: alpha_global[1][1],
                zz: alpha_global[2][2],
                xy: alpha_global[0][1],
                yz: alpha_global[1][2],
                xz: alpha_global[0][2],
            }),
            euler_deg: [0.0; 3],
            ..Default::default()
        },
    );
    p_ref_sym.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
    ];
    p_ref_sym.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
    let res_ref = solve_solid(&test_mesh, &p_ref_sym).unwrap();
    let mut spec_swapped = spec.clone();
    if let MaterialAssign::PropsTensor(map) = spec_swapped.materials.get_mut("body").unwrap() {
        // Swap xy (idx 3) and yz (idx 5) in Sym6 to prove shear components are distinct & ordered!
        map.insert(
            "alpha".to_string(),
            PropValue::Sym6([
                alpha_global[0][0],
                alpha_global[1][1],
                alpha_global[2][2],
                alpha_global[1][2],
                alpha_global[0][2],
                alpha_global[0][1],
            ]),
        );
    }
    let res_swapped = Sms.solve(&raw_ex, &spec_swapped, &mut Ctx::new(&BoxResolver)).unwrap();
    let mut swap_diff = 0.0_f64;
    for comp in 0..3 {
        let u_num = res_sym6.solutions[0].field("u").unwrap().component(comp);
        let u_ref = res_ref.solutions[0].field("u").unwrap().component(comp);
        let u_swp = res_swapped.solutions[0].field("u").unwrap().component(comp);
        let scale = u_ref.iter().map(|v| v.abs()).fold(1e-12_f64, f64::max);
        for i in 0..test_mesh.n_dofs(2) {
            assert!(
                (u_num[i] - u_ref[i]).abs() < 1e-13 * scale,
                "Sym6 [0,1,2,5,4,3] ordering mismatch at dof {i} comp {comp}: {} vs {}",
                u_num[i],
                u_ref[i]
            );
            swap_diff = swap_diff.max((u_num[i] - u_swp[i]).abs() / scale);
        }
    }
    assert!(swap_diff > 1e-3, "Sym6 shear permutation sensitivity too small: {swap_diff:.3e}");
}

// ---- Pillar (c): Formal Eigenfrequency p_obs & R6d Axisymmetric PML p_obs ----------------------------

#[test]
fn symmetric_eigenfrequency_and_axisym_pml_r6d_convergence_rates() {
    // 1. 2D Waveguide SH(1, 0) exact frequency on [0, L] x [0, H] -> P1 O(h^2), P2 O(h^4)
    let (w, h_plate, e, nu, rho) = (1.0_f64, 0.25_f64, 169.0e9_f64, 0.22_f64, 2329.0_f64);
    let mu = e / (2.0 * (1.0 + nu));
    let cs = (mu / rho).sqrt();
    let f_sh_exact = cs / (2.0 * w);

    // Build structured nested 2D tri meshes for SH(1, 0) so unstructured Delaunay jitter doesn't distort p_obs:
    let rect_struct = |lx: f64, ly: f64, nx: usize, ny: usize| -> SimplexMesh {
        let mut pts = Vec::with_capacity(2 * (nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                pts.push(lx * i as f64 / nx as f64);
                pts.push(ly * j as f64 / ny as f64);
            }
        }
        let vid = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
        let mut cells = Vec::new();
        let mut clbls = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let (a, b, c, d) = (vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1));
                cells.extend_from_slice(&[a, b, c, a, c, d]);
                clbls.extend_from_slice(&[1, 1]);
            }
        }
        let mut facets = Vec::new();
        let mut flbls = Vec::new();
        for j in 0..ny {
            facets.extend_from_slice(&[vid(0, j), vid(0, j + 1), vid(nx, j), vid(nx, j + 1)]);
            flbls.extend_from_slice(&[1, 4]);
        }
        for i in 0..nx {
            facets.extend_from_slice(&[vid(i, 0), vid(i + 1, 0), vid(i, ny), vid(i + 1, ny)]);
            flbls.extend_from_slice(&[2, 3]);
        }
        SimplexMesh::new(2, pts, cells, clbls, facets, flbls)
    };
    let ns_2d = [4usize, 8, 16];
    let hs_2d: Vec<f64> = ns_2d.iter().map(|&n| 1.0 / n as f64).collect();
    for (order, expected_p, band) in [(1usize, 2.0, (1.85, 2.25)), (2usize, 4.0, (3.60, 4.30))] {
        let mut errs = Vec::new();
        for &n in &ns_2d {
            let mesh = rect_struct(w, h_plate, n, n / 4);
            let mut p = Problem::new(Study::Modes { n_modes: 4, f_search: 0.9 * f_sh_exact, wavevector: 0.0 }, order);
            p.materials.set(1, "E", e);
            p.materials.set(1, "nu", nu);
            p.materials.set(1, "rho", rho);
            p.bcs = vec![Bc { labels: vec![1, 4], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
            let r = common::solve(&mesh, &p);
            let f_got = (1..=r.solutions.len())
                .filter(|&j| match r.scalars.iter().find(|(k, _)| k == &format!("longitudinal share {j}")) {
                    Some((_, Scalar::Num(x))) => *x > 0.99,
                    _ => false,
                })
                .map(|j| r.solutions[j - 1].value.unwrap())
                .min_by(|a, b| (a - f_sh_exact).abs().total_cmp(&(b - f_sh_exact).abs()))
                .unwrap();
            errs.push((f_got - f_sh_exact).abs() / f_sh_exact);
        }
        let rep = assert_convergence_rate(&hs_2d, &errs, expected_p, band);
        println!("SH eigenfrequency P{order}: {rep:?}");
    }

    // 2. 3D Guided-face cube T(1,1,0) = sqrt(2)*c_s/(2L) on structured Kuhn-tet meshes -> P1 O(h^2), P2 O(h^4)
    let (e3, nu3, rho3) = (200.0e9_f64, 0.25_f64, 7850.0_f64);
    let cp3 = (e3 * (1.0 - nu3) / ((1.0 + nu3) * (1.0 - 2.0 * nu3)) / rho3).sqrt();
    let cs3 = (e3 / (2.0 * (1.0 + nu3)) / rho3).sqrt();
    let f_lame = (2.0_f64).sqrt() * cs3 / 2.0; // on unit cube [0,1]^3 -> T(1,1,0)
    let _ = cp3;
    let ns_3d = [4usize, 6, 8];
    let hs_3d: Vec<f64> = ns_3d.iter().map(|&n| 1.0 / n as f64).collect();
    for (order, expected_p, band) in [(1usize, 2.0, (1.85, 2.25)), (2usize, 4.0, (3.60, 4.30))] {
        let mut errs = Vec::new();
        for &n in &ns_3d {
            let mesh = kuhn_cube(n);
            let mut p = Problem::new(Study::Eigenfrequency { n_modes: 6, f_search: 0.9 * f_lame }, order);
            p.solver = LinearSolver::Direct;
            p.materials.set(1, "E", e3);
            p.materials.set(1, "nu", nu3);
            p.materials.set(1, "rho", rho3);
            p.bcs = vec![
                Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
                Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
                Bc { labels: vec![3, 4], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
            ];
            let r = common::solve(&mesh, &p);
            let err =
                r.solutions.iter().map(|s| (s.value.unwrap() - f_lame).abs() / f_lame).fold(f64::INFINITY, f64::min);
            errs.push(err);
        }
        let rep = assert_convergence_rate(&hs_3d, &errs, expected_p, band);
        println!("3D guided cube eigenfrequency P{order}: {rep:?}");
    }

    // 3. R6d: Axisymmetric cylindrical PML with quadrature-point s_theta = r~/r must achieve optimal O(h^2) (P1) / O(h^4) (P2)
    //    convergence rate (no O(h/r) cell-center floor) against the exact physical longer ring!
    //    Use roller-guided ends (u_z = 0 on z = 0, H; u_r = 0 on r = r0) so the radial eigenmode is C^inf smooth
    //    without corner singularities masking P2 O(h^4) accuracy.
    let (r0, w_pml, h_ring, re_stretch) = (1.0_f64, 0.3_f64, 0.1_f64, 1.5_f64);
    let roller_bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, false], values: [0.0; 3] } },
        Bc { labels: vec![2, 3], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
    ];
    let ref_long = {
        let mut m = rect_struct(w_pml * re_stretch, h_ring, 64, 16);
        m.points.chunks_mut(2).for_each(|p| p[0] += r0);
        let mesh_l = SimplexMesh::new(2, m.points, m.cells, m.cell_labels, m.facets, m.facet_labels);
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
        p.materials.set(1, "E", 100.0e9);
        p.materials.set(1, "nu", 0.3);
        p.materials.set(1, "rho", 2500.0);
        p.bcs = roller_bcs.clone();
        let (_, modes) = solve_solid_axisym(&mesh_l, &p, 0.0).unwrap();
        modes.freqs[0].re
    };
    let grids_ax = [4usize, 8, 12];
    let hs_ax: Vec<f64> = grids_ax.iter().map(|&n| 1.0 / n as f64).collect();
    for (order, expected_p, band) in [(1usize, 2.0, (1.85, 2.25)), (2usize, 4.0, (3.60, 4.30))] {
        let mut errs = Vec::new();
        for &n in &grids_ax {
            let mut m = rect_struct(w_pml, h_ring, n, (n / 4).max(1));
            m.points.chunks_mut(2).for_each(|p| p[0] += r0);
            let mesh = SimplexMesh::new(2, m.points, m.cells, m.cell_labels, m.facets, m.facet_labels);
            let mut p = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, order);
            p.materials.set(1, "E", 100.0e9);
            p.materials.set(1, "nu", 0.3);
            p.materials.set(1, "rho", 2500.0);
            p.bcs = roller_bcs.clone();
            p.materials.pml.insert(
                1,
                Pml {
                    x: PmlDir::Plus,
                    strength: 0.0,
                    order: 0.0,
                    re: re_stretch,
                    clamp_outer: false,
                    ..Default::default()
                },
            );
            let (_, modes) = solve_solid_axisym(&mesh, &p, 0.0).unwrap();
            errs.push((modes.freqs[0].re - ref_long).abs() / ref_long);
        }
        let rep = assert_convergence_rate(&hs_ax, &errs, expected_p, band);
        println!("R6d axisym PML P{order} convergence: {rep:?}");
    }
}

// ---- Rough-Spot Closures R6a, R6b, R32 ---------------------------------------------------------------

#[test]
fn r6a_axisym_m_plus_minus_1_enforces_exact_algebraic_mpc_on_r0() {
    let (e, nu, rho, r_cyl, len) = (169.0e9_f64, 0.28_f64, 2329.0_f64, 0.5_f64, 2.0_f64);
    let mesh = common::parts_mesh(&[common::part(0.0, 0.0, r_cyl, len, 1, Some(0.1))], 0.1);
    let bottom_lbls: Vec<i32> = (0..mesh.n_facets())
        .filter(|&f| mesh.facet(f).iter().all(|&v| mesh.point(v as usize)[1].abs() < 1e-12))
        .map(|f| mesh.facet_labels[f])
        .collect();

    let tol = nyx_sms::axisym::axis_tolerance(&mesh);
    let coords = mesh.dof_coordinates(2);
    let n1 = mesh.n_dofs(2);
    let axis_dofs: Vec<usize> = (0..n1).filter(|&d| coords[2 * d] <= tol).collect();
    assert!(!axis_dofs.is_empty());

    // Test both m = +1 and m = -1 via SimplexMesh kernel and MixedDiscretisation + Sms::run!
    let raw = Mesh::from_solver_mesh(&mesh).unwrap();
    let mut freqs_pos_neg = Vec::new();
    for m_ord in [1_i64, -1_i64] {
        let mut p =
            Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search: 0.0, azimuthal_order: m_ord }, 2);
        p.materials.set(1, "E", e);
        p.materials.set(1, "nu", nu);
        p.materials.set(1, "rho", rho);
        p.bcs = vec![Bc { labels: bottom_lbls.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];

        let asm = assemble_cross_section(&mesh, &p, CrossSection::Axisymmetric(m_ord as f64)).unwrap();
        assert!(asm.red.is_some(), "R6a: m={m_ord} touching r=0 must build Reduction MPC");
        let solved = solve_cross_section_modes(&asm).unwrap();
        freqs_pos_neg.push(solved.modes.freqs[0].re);

        // Verify exact algebraic MPC u_r(0, z) - i m u_theta(0, z) == 0 and u_z(0, z) == 0 on EVERY mode & axis DOF!
        let factor = c(0.0, m_ord as f64);
        for vec in &solved.modes.vectors {
            let max_norm = vec.iter().map(|z| z.norm()).fold(0.0_f64, f64::max);
            for &d in &axis_dofs {
                let ur = vec[d];
                let uz = vec[n1 + d];
                let ut = vec[2 * n1 + d];
                assert!(uz.norm() < 1e-13 * max_norm, "m={m_ord}: uz on r=0 not zero: {uz:?}");
                let tie_res = (ur - factor * ut).norm();
                assert!(
                    tie_res < 1e-13 * max_norm,
                    "R6a MPC failed at axis dof {d} (m={m_ord}): |u_r - i*m*u_theta| = {tie_res:.3e}"
                );
            }
        }

        // Also check log documents u_r = i m u_theta regularity and run(&raw, &p) matches.
        let r_run = run(&raw, &p, &Progress::new()).unwrap();
        assert!(
            r_run.log.iter().any(|l| l.contains("u_r = i m u_theta")),
            "expected MPC regularity log message, got {:?}",
            r_run.log
        );
        assert!((r_run.solutions[0].value.unwrap() - solved.modes.freqs[0].re).abs() < 1e-10);
    }
    assert!(
        (freqs_pos_neg[0] - freqs_pos_neg[1]).abs() < 1e-11 * freqs_pos_neg[0],
        "m=+1 ({}) vs m=-1 ({}) degeneracy mismatch",
        freqs_pos_neg[0],
        freqs_pos_neg[1]
    );
}

#[test]
fn r6b_axisym_lrb_applies_radial_jacobian_weight_and_energy_flux_balance() {
    // 1. Unit check of `lrb::damping_matrix_axisym` radial Jacobian integral:
    //    Horizontal facet y = 0 spanning r in [r0, r1], n = (0, -1, 0):
    //    sum of u_theta block rows = int_{r0}^{r1} (rho c_s) r dr = 0.5 * rho * c_s * (r1^2 - r0^2).
    let (r0, r1, h) = (2.0_f64, 5.0_f64, 0.5_f64);
    let (e, nu, rho) = (100.0e9_f64, 0.25_f64, 2500.0_f64);
    let mut base = common::rect(r1 - r0, h, 8, 2);
    base.points.chunks_mut(2).for_each(|p| p[0] += r0);
    let mesh = SimplexMesh::new(2, base.points, base.cells, base.cell_labels, base.facets, base.facet_labels);

    let c_iso = femcore::material::isotropic_stiffness(e, nu);
    let c6 = vec![plane_stiffness(&c_iso, 2); mesh.n_cells()];
    let rho_cells = vec![rho; mesh.n_cells()];
    let (cp, cs) = nyx_sms::lrb::voigt_speeds(&c6[0], rho);

    let n1 = mesh.n_dofs(2);
    let c_ax = damping_matrix_axisym(&mesh, &c6, &rho_cells, &[2], 2, 3);
    // Unit test vectors in block 2 (u_theta, tangential -> c_s) and block 1 (u_z, normal -> c_p):
    let mut u_theta = vec![0.0; 3 * n1];
    let mut u_z = vec![0.0; 3 * n1];
    for d in 0..n1 {
        u_theta[2 * n1 + d] = 1.0;
        u_z[n1 + d] = 1.0;
    }
    let q_t: f64 = u_theta.iter().zip(c_ax.matvec(&u_theta)).map(|(a, b)| a * b).sum();
    let q_n: f64 = u_z.iter().zip(c_ax.matvec(&u_z)).map(|(a, b)| a * b).sum();
    let r_int_exact = 0.5 * (r1 * r1 - r0 * r0);
    let want_t = rho * cs * r_int_exact;
    let want_n = rho * cp * r_int_exact;
    assert!((q_t - want_t).abs() < 1e-12 * want_t, "axisym LRB shear integral {q_t:.6e} vs {want_t:.6e}");
    assert!((q_n - want_n).abs() < 1e-12 * want_n, "axisym LRB normal integral {q_n:.6e} vs {want_n:.6e}");

    // 2. Energy-flux & quality-factor identity on library entry `assemble_cross_section`:
    //    On a ring with outer radial LRB (right boundary = 4), (K + i omega C - omega^2 M) u = 0
    //    implies exact quadratic eigen-energy flux identity Q_eig = Re(w)/(2 Im(w)) == Re(w) * (u^H M u) / (u^H C u).
    let mut p_ax = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 2.0e4 }, 2);
    p_ax.materials.set(1, "E", e);
    p_ax.materials.set(1, "nu", nu);
    p_ax.materials.set(1, "rho", rho);
    p_ax.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![4], kind: BcKind::LowReflecting },
    ];
    let asm = assemble_cross_section(&mesh, &p_ax, CrossSection::Axisymmetric(0.0)).unwrap();
    let cf = asm.cf.as_ref().expect("R6b library entry must assemble axisymmetric cf");
    let solved = solve_cross_section_modes(&asm).unwrap();
    let quad_re = |m: &femcore::linalg::csr::CCsr, v: &[C]| -> f64 {
        let mv = m.matvec(v);
        v.iter().zip(&mv).map(|(a, b)| (a.conj() * b).re).sum()
    };
    let cf_c = cf.to_complex();
    let free = &asm.free;
    for (idx, fz) in solved.modes.freqs.iter().enumerate() {
        let omega = fz * (2.0 * PI);
        let u_full = &solved.modes.vectors[idx];
        let u_free: Vec<C> = free.iter().map(|&d| u_full[d]).collect();
        let m_q = quad_re(&asm.mf, &u_free);
        let c_q = quad_re(&cf_c, &u_free);
        assert!(omega.im > 0.0 && c_q > 0.0, "expected passive outward radiation damping");
        // Radiated energy flux balance: 2 * Re(omega) * Im(omega) * (u^H M u) == Re(omega) * (u^H C u)
        let lhs = 2.0 * omega.re * omega.im * m_q;
        let rhs = omega.re * c_q;
        assert!(
            (lhs - rhs).abs() < 1e-6 * lhs.abs().max(rhs.abs()),
            "mode {} energy flux mismatch: lhs={lhs:.6e} vs rhs={rhs:.6e}",
            idx + 1
        );
    }

    // 3. Also test R6a + R6b combined (`m = 1` solid cylinder on axis + outer radial LRB via library entry)!
    let cyl = common::parts_mesh(&[common::part(0.0, 0.0, 1.0, 0.5, 1, Some(0.125))], 0.125);
    let mut right_lbls: Vec<i32> = (0..cyl.n_facets())
        .filter(|&f| cyl.facet(f).iter().all(|&v| (cyl.point(v as usize)[0] - 1.0).abs() < 1e-9))
        .map(|f| cyl.facet_labels[f])
        .collect();
    right_lbls.sort_unstable();
    right_lbls.dedup();
    let mut p_m1_lrb = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 2.0e3 }, 2);
    p_m1_lrb.materials.set(1, "E", e);
    p_m1_lrb.materials.set(1, "nu", nu);
    p_m1_lrb.materials.set(1, "rho", rho);
    p_m1_lrb.bcs = vec![Bc { labels: right_lbls, kind: BcKind::LowReflecting }];
    let (_, m1_modes) = solve_solid_axisym(&cyl, &p_m1_lrb, 1.0).unwrap();
    assert!(m1_modes.freqs.iter().any(|z| z.re > 100.0 && z.im > 0.0), "{:?}", m1_modes.freqs);
}

#[test]
fn r32_curved_and_qp_meshes_emit_explicit_unsupported_on_incompatible_paths() {
    // 1. Curved SimplexMesh passed directly to affine simplex entry points MUST emit explicit "unsupported" error
    //    rather than silently ignoring edge_nodes!
    let mut curved_tri = common::rect(1.0, 0.5, 4, 2);
    let (a, b) = curved_tri.edges()[0];
    let pa = curved_tri.point(a as usize);
    let pb = curved_tri.point(b as usize);
    curved_tri.set_edge_nodes(&[([a, b], [0.5 * (pa[0] + pb[0]) + 0.01, 0.5 * (pa[1] + pb[1]), 0.0])]);
    assert!(curved_tri.has_curved());

    let mut sp = Problem::new(Study::Static, 2);
    sp.materials.set(1, "E", 100e9);
    sp.materials.set(1, "nu", 0.3);
    sp.materials.set(1, "rho", 2500.0);
    sp.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];

    let e_solid = solve_solid(&curved_tri, &sp).unwrap_err();
    assert!(e_solid.contains("unsupported"), "expected unsupported error from solve_solid, got: {e_solid}");
    let e_setup = setup_with_source::<2, _>(&curved_tri, &sp, &femcore::vvuq::UniformSource)
        .err()
        .expect("expected unsupported error from setup_with_source");
    assert!(e_setup.contains("unsupported"), "got: {e_setup}");

    let mut sp_ax = sp.clone();
    sp_ax.study = Study::AxisymmetricEigenfrequency { n_modes: 3, f_search: 0.0, azimuthal_order: 0 };
    let e_ax = solve_axisymmetric(&curved_tri, &sp_ax).unwrap_err();
    assert!(e_ax.contains("unsupported"), "got: {e_ax}");

    let mut sp_wg = sp.clone();
    sp_wg.study = Study::Modes { n_modes: 3, f_search: 0.0, wavevector: 1.0 };
    let e_wg = solve_waveguide_modes(&curved_tri, &sp_wg).unwrap_err();
    assert!(e_wg.contains("unsupported"), "got: {e_wg}");

    let mut sp_pml = sp.clone();
    sp_pml.study = Study::Eigenfrequency { n_modes: 3, f_search: 0.0 };
    sp_pml.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 1.0, order: 1.0, ..Default::default() });
    let e_pml = nyx_sms::pml::solve_solid_pml(&curved_tri, &sp_pml).unwrap_err();
    assert!(e_pml.contains("unsupported"), "got: {e_pml}");

    let mut sp_pz = sp.clone();
    sp_pz.piezoelectric = PiezoMode::On;
    sp_pz.materials.set(1, "eps_r", 9.0);
    sp_pz.bcs.push(Bc { labels: vec![1], kind: BcKind::Ground });
    let e_pz = nyx_sms::piezo::solve_piezo(&curved_tri, &sp_pz).unwrap_err();
    assert!(e_pz.contains("unsupported"), "got: {e_pz}");

    // 2. Barsoum quarter-point edge (`qp_edges`) passed to `Sms::validate` / `run` MUST return `Error::Unsupported`.
    let mut cm = CellMesh::<2>::new();
    cm.points = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    cm.push_cell(CellType::Tri, &[0, 1, 2], 1);
    cm.push_cell(CellType::Tri, &[0, 2, 3], 1);
    cm.push_facet(&[0, 3], 1);
    cm.push_facet(&[1, 2], 4);
    cm.set_edge_nodes(vec![([0, 1], [0.25, 0.0])]);
    cm.qp_edges = vec![[0, 1]];
    let qp_mesh = Mesh::D2(cm);
    struct R32Res;
    impl LabelResolver for R32Res {
        fn domains(&self, _: &str) -> Result<Vec<i32>, Error> {
            Ok(vec![1])
        }
        fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "left" => Ok(vec![1]),
                _ => Ok(vec![4]),
            }
        }
    }
    let mut spec = ProblemSpec::new("sms", 2, "static");
    spec.materials.insert("dom".into(), MaterialAssign::Props([("E".into(), 100e9), ("nu".into(), 0.3)].into()));
    spec.bcs.push(BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0 }) });
    let err_qp = Sms.solve(&qp_mesh, &spec, &mut Ctx::new(&R32Res)).unwrap_err();
    assert_eq!(err_qp.code(), "unsupported", "expected Error::Unsupported for qp_edges, got {err_qp:?}");

    let err_run_qp = run(&qp_mesh, &sp, &Progress::new()).unwrap_err();
    assert_eq!(err_run_qp.code(), "unsupported", "expected Error::Unsupported from run(), got {err_run_qp:?}");
}
