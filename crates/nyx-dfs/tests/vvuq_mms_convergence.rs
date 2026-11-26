// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME Layer 2 Pillars (b) & (c): Zero-sized monomorphized `SourceField<D>` Manufactured Solutions (MMS)
//! and formal `assert_convergence_rate` verification suite across P1/P2 Electrostatics (2D/3D, iso & rotated
//! aniso, straight & non-QP `CurvedTet`, mixed Dirichlet/Neumann/Robin), Stationary Heat, and `mixed_3d`
//! Nedelec Magnetostatics using degree >= 6 quadrature.

use femcore::fem::{geometry, quad_rule_collapsed, shape, CurvedTet, TET_EDGES};
use femcore::material::{tensor_global, Anisotropy, Tensor, M3};
use femcore::mesh::SimplexMesh;
use femcore::section::plane_tensor;
use femcore::vvuq::{assert_convergence_rate, FnSource};
use nyx_dfs::electro::solve_electrostatics_with_source_and_robin;
use nyx_dfs::heat::solve_heat_with_source;
use nyx_dfs::magneto::solve_mixed_3d_edges_with_source;
use nyx_dfs::{Bc, BcKind, Problem, Study, EPS0, MU0};
use std::f64::consts::PI;

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
            lam.push([1.0 - u - (1.0 - u) * v, u, (1.0 - u) * v]);
            w.push(2.0 * (1.0 - u) * du * dv);
        }
    }
    (lam, w)
}

/// Exact value, gradient, and Hessian of phi*(x, y, z) = sin(pi x) sin(pi y) cos(0.5 pi z) + 0.3 exp(0.4 x - 0.3 y + 0.2 z).
/// In 2D (`y_in = z_mesh, z_in = 0`), reduces cleanly to sin(pi x) sin(pi z) + 0.3 exp(0.4 x - 0.3 z).
fn phi_jet(x: [f64; 3]) -> (f64, [f64; 3], M3) {
    let [px, py, pz] = x;
    let sx = (PI * px).sin();
    let cx = (PI * px).cos();
    let sy = (PI * py).sin();
    let cy = (PI * py).cos();
    let cz = (0.5 * PI * pz).cos();
    let sz = (0.5 * PI * pz).sin();

    let ex = 0.3 * (0.4 * px - 0.3 * py + 0.2 * pz).exp();
    let val = sx * sy * cz + ex;

    let g0 = PI * cx * sy * cz + 0.4 * ex;
    let g1 = PI * sx * cy * cz - 0.3 * ex;
    let g2 = -0.5 * PI * sx * sy * sz + 0.2 * ex;

    let h00 = -PI * PI * sx * sy * cz + 0.16 * ex;
    let h11 = -PI * PI * sx * sy * cz + 0.09 * ex;
    let h22 = -0.25 * PI * PI * sx * sy * cz + 0.04 * ex;
    let h01 = PI * PI * cx * cy * cz - 0.12 * ex;
    let h02 = -0.5 * PI * PI * cx * sy * sz + 0.08 * ex;
    let h12 = -0.5 * PI * PI * sx * cy * sz - 0.06 * ex;
    (val, [g0, g1, g2], [[h00, h01, h02], [h01, h11, h12], [h02, h12, h22]])
}

fn scalar_mms_norms_2d(mesh: &SimplexMesh, order: usize, u: &[f64]) -> (f64, f64) {
    let (qlam, qw) = quad_tri_7();
    let geo = geometry(mesh);
    let shapes: Vec<(Vec<f64>, Vec<f64>)> = qlam.iter().map(|l| shape(l, order)).collect();
    let nloc = mesh.nloc(order);
    let mut dofs = Vec::new();
    let (mut l2_2, mut h1_2) = (0.0_f64, 0.0_f64);
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let pts = [mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize)];
        mesh.cell_dofs(c, order, &mut dofs);
        let g = [
            [geo.g[c * 6], geo.g[c * 6 + 1]],
            [geo.g[c * 6 + 2], geo.g[c * 6 + 3]],
            [geo.g[c * 6 + 4], geo.g[c * 6 + 5]],
        ];
        for (q, lam) in qlam.iter().enumerate() {
            let x = lam[0] * pts[0][0] + lam[1] * pts[1][0] + lam[2] * pts[2][0];
            let z = lam[0] * pts[0][1] + lam[1] * pts[1][1] + lam[2] * pts[2][1];
            let (nv, dn) = &shapes[q];
            let mut uh = 0.0;
            let mut guh = [0.0; 2];
            for a in 0..nloc {
                let ud = u[dofs[a]];
                uh += nv[a] * ud;
                for b in 0..3 {
                    let db = dn[a * 3 + b] * ud;
                    guh[0] += db * g[b][0];
                    guh[1] += db * g[b][1];
                }
            }
            let (ve, ge, _) = phi_jet([x, z, 0.0]);
            let wq = geo.vol[c] * qw[q];
            l2_2 += wq * (uh - ve).powi(2);
            h1_2 += wq * ((guh[0] - ge[0]).powi(2) + (guh[1] - ge[1]).powi(2));
        }
    }
    (l2_2.sqrt(), h1_2.sqrt())
}

fn scalar_mms_norms_3d(mesh: &SimplexMesh, order: usize, u: &[f64]) -> (f64, f64) {
    let qrule = if mesh.has_curved() { quad_rule_collapsed(6) } else { quad_rule_collapsed(4) };
    let geo = geometry(mesh);
    let shapes: Vec<(Vec<f64>, Vec<f64>)> = qrule.lam.iter().map(|l| shape(l, order)).collect();
    let nloc = mesh.nloc(order);
    let mut dofs = Vec::new();
    let (mut l2_2, mut h1_2) = (0.0_f64, 0.0_f64);
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        mesh.cell_dofs(c, order, &mut dofs);
        if mesh.cell_is_curved(c) {
            let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
            for (q, lam_v) in qrule.lam.iter().enumerate() {
                let lam = [lam_v[0], lam_v[1], lam_v[2], lam_v[3]];
                let x = ct.position(&lam);
                let (gq, det) = ct.jacobian(&lam);
                let (nv, dn) = &shapes[q];
                let mut uh = 0.0;
                let mut guh = [0.0; 3];
                for a in 0..nloc {
                    let ud = u[dofs[a]];
                    uh += nv[a] * ud;
                    for b in 0..4 {
                        let db = dn[a * 4 + b] * ud;
                        for k in 0..3 {
                            guh[k] += db * gq[b][k];
                        }
                    }
                }
                let (ve, ge, _) = phi_jet(x);
                let wq = qrule.w[q] * det.abs() / 6.0;
                l2_2 += wq * (uh - ve).powi(2);
                h1_2 += wq * ((guh[0] - ge[0]).powi(2) + (guh[1] - ge[1]).powi(2) + (guh[2] - ge[2]).powi(2));
            }
        } else {
            let pts = [
                mesh.point(v[0] as usize),
                mesh.point(v[1] as usize),
                mesh.point(v[2] as usize),
                mesh.point(v[3] as usize),
            ];
            let g = [
                [geo.g[c * 12], geo.g[c * 12 + 1], geo.g[c * 12 + 2]],
                [geo.g[c * 12 + 3], geo.g[c * 12 + 4], geo.g[c * 12 + 5]],
                [geo.g[c * 12 + 6], geo.g[c * 12 + 7], geo.g[c * 12 + 8]],
                [geo.g[c * 12 + 9], geo.g[c * 12 + 10], geo.g[c * 12 + 11]],
            ];
            for (q, lam) in qrule.lam.iter().enumerate() {
                let mut x = [0.0; 3];
                for b in 0..4 {
                    for k in 0..3 {
                        x[k] += lam[b] * pts[b][k];
                    }
                }
                let (nv, dn) = &shapes[q];
                let mut uh = 0.0;
                let mut guh = [0.0; 3];
                for a in 0..nloc {
                    let ud = u[dofs[a]];
                    uh += nv[a] * ud;
                    for b in 0..4 {
                        let db = dn[a * 4 + b] * ud;
                        for k in 0..3 {
                            guh[k] += db * g[b][k];
                        }
                    }
                }
                let (ve, ge, _) = phi_jet(x);
                let wq = geo.vol[c] * qrule.w[q];
                l2_2 += wq * (uh - ve).powi(2);
                h1_2 += wq * ((guh[0] - ge[0]).powi(2) + (guh[1] - ge[1]).powi(2) + (guh[2] - ge[2]).powi(2));
            }
        }
    }
    (l2_2.sqrt(), h1_2.sqrt())
}

fn mms_2d_tensor_source(
    a2: [[f64; 2]; 2],
    robin_alpha: f64,
    robin_labs: &[i32],
) -> impl femcore::vvuq::SourceField<2> + '_ {
    FnSource {
        volume_rhs: move |x: [f64; 2], _t: f64, _dom: u32| -> [f64; 3] {
            let (_, _, h) = phi_jet([x[0], x[1], 0.0]);
            let rhs = -(a2[0][0] * h[0][0] + 2.0 * a2[0][1] * h[0][1] + a2[1][1] * h[1][1]);
            [rhs, 0.0, 0.0]
        },
        dirichlet: move |x: [f64; 2], _t: f64, _bnd: u32| -> [f64; 3] { [phi_jet([x[0], x[1], 0.0]).0, 0.0, 0.0] },
        neumann_robin: move |x: [f64; 2], n: [f64; 2], _t: f64, bnd: u32| -> [f64; 3] {
            let (val, g, _) = phi_jet([x[0], x[1], 0.0]);
            let flux = n[0] * (a2[0][0] * g[0] + a2[0][1] * g[1]) + n[1] * (a2[1][0] * g[0] + a2[1][1] * g[1]);
            let rob = if robin_labs.contains(&(bnd as i32)) { robin_alpha * val } else { 0.0 };
            [flux + rob, 0.0, 0.0]
        },
    }
}

fn mms_3d_tensor_source(a3: M3, robin_alpha: f64, robin_labs: &[i32]) -> impl femcore::vvuq::SourceField<3> + '_ {
    FnSource {
        volume_rhs: move |x: [f64; 3], _t: f64, _dom: u32| -> [f64; 3] {
            let (_, _, h) = phi_jet(x);
            let mut s = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    s -= a3[i][j] * h[i][j];
                }
            }
            [s, 0.0, 0.0]
        },
        dirichlet: move |x: [f64; 3], _t: f64, _bnd: u32| -> [f64; 3] { [phi_jet(x).0, 0.0, 0.0] },
        neumann_robin: move |x: [f64; 3], n: [f64; 3], _t: f64, bnd: u32| -> [f64; 3] {
            let (val, g, _) = phi_jet(x);
            let mut flux = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    flux += n[i] * a3[i][j] * g[j];
                }
            }
            let rob = if robin_labs.contains(&(bnd as i32)) { robin_alpha * val } else { 0.0 };
            [flux + rob, 0.0, 0.0]
        },
    }
}

#[test]
fn electro_and_heat_2d_mms_p1_p2_iso_and_rotated_aniso_mixed_bcs() {
    let aniso_t = Tensor::Diagonal { x: 1.5, y: 2.2, z: 3.1 };
    let euler = [90.0, 32.0, -90.0]; // rotation in (x, z) wafer plane
    let a3_rel = tensor_global(&aniso_t, euler);
    let a2_full = plane_tensor(&a3_rel, 2);
    let a2_rel = [[a2_full[0][0], a2_full[0][1]], [a2_full[1][0], a2_full[1][1]]];
    let iso2_rel = [[2.5, 0.0], [0.0, 2.5]];

    let ns = &[6usize, 12, 24][..];
    let hs: Vec<f64> = ns.iter().map(|&n| 1.0 / n as f64).collect();

    for (case_name, is_aniso) in [("iso", false), ("rot_aniso", true)] {
        for &order in &[1usize, 2] {
            // Electro 2D (Dirichlet on 1,2; Neumann on 3; Robin on 4)
            let robin_alpha_e = 2.0 * EPS0;
            let a_eff_e = if is_aniso {
                [[EPS0 * a2_rel[0][0], EPS0 * a2_rel[0][1]], [EPS0 * a2_rel[1][0], EPS0 * a2_rel[1][1]]]
            } else {
                [[EPS0 * iso2_rel[0][0], 0.0], [0.0, EPS0 * iso2_rel[1][1]]]
            };
            let src_e = mms_2d_tensor_source(a_eff_e, robin_alpha_e, &[4]);
            let (mut l2s_e, mut h1s_e) = (Vec::new(), Vec::new());
            for &n in ns {
                let mesh = femcore::testkit::rect(1.0, 1.0, n, n);
                let mut p = Problem::new(Study::Electrostatic, order);
                if is_aniso {
                    p.materials.aniso.insert(
                        mesh.cell_labels[0],
                        Anisotropy { permittivity: Some(aniso_t.clone()), euler_deg: euler, ..Default::default() },
                    );
                } else {
                    p.materials.default.insert("eps_r".into(), 2.5);
                }
                p.bcs = vec![
                    Bc { labels: vec![1, 2], kind: BcKind::Potential(0.0) },
                    Bc { labels: vec![3], kind: BcKind::SurfaceCharge(0.0) },
                ];
                let res =
                    solve_electrostatics_with_source_and_robin::<2, _>(&mesh, &p, &src_e, &[(vec![4], robin_alpha_e)])
                        .unwrap();
                let v = &res.solutions[0].field("V").unwrap().values;
                let (e_l2, e_h1) = scalar_mms_norms_2d(&mesh, order, v);
                l2s_e.push(e_l2);
                h1s_e.push(e_h1);
            }
            let (l2_exp, l2_band) = if order == 1 { (2.0, (1.85, 2.25)) } else { (3.0, (2.80, 3.25)) };
            let (h1_exp, h1_band) = if order == 1 { (1.0, (0.90, 1.15)) } else { (2.0, (1.85, 2.20)) };
            let r_l2 = assert_convergence_rate(&hs, &l2s_e, l2_exp, l2_band);
            let r_h1 = assert_convergence_rate(&hs, &h1s_e, h1_exp, h1_band);
            eprintln!(
                "Electro 2D {case_name} P{order}: L2={:.3}, H1={:.3}",
                r_l2.asymptotic_slope, r_h1.asymptotic_slope
            );

            // Heat 2D (Dirichlet on 1,2; Flux on 3; Convection/Robin on 4)
            let robin_alpha_h = 2.0;
            let a_eff_h = if is_aniso { a2_rel } else { iso2_rel };
            let src_h = mms_2d_tensor_source(a_eff_h, robin_alpha_h, &[4]);
            let (mut l2s_h, mut h1s_h) = (Vec::new(), Vec::new());
            for &n in ns {
                let mesh = femcore::testkit::rect(1.0, 1.0, n, n);
                let mut p = Problem::new(Study::HeatStationary, order);
                if is_aniso {
                    p.materials.aniso.insert(
                        mesh.cell_labels[0],
                        Anisotropy { conductivity: Some(aniso_t.clone()), euler_deg: euler, ..Default::default() },
                    );
                } else {
                    p.materials.default.insert("k".into(), 2.5);
                }
                p.bcs = vec![
                    Bc { labels: vec![1, 2], kind: BcKind::Temperature(0.0) },
                    Bc { labels: vec![3], kind: BcKind::Flux(0.0) },
                    Bc { labels: vec![4], kind: BcKind::Convection { h: robin_alpha_h, text: 0.0 } },
                ];
                let res = solve_heat_with_source::<2, _>(&mesh, &p, &src_h).unwrap();
                let t = &res.solutions[0].field("T").unwrap().values;
                let (e_l2, e_h1) = scalar_mms_norms_2d(&mesh, order, t);
                l2s_h.push(e_l2);
                h1s_h.push(e_h1);
            }
            let r_hl2 = assert_convergence_rate(&hs, &l2s_h, l2_exp, l2_band);
            let r_hh1 = assert_convergence_rate(&hs, &h1s_h, h1_exp, h1_band);
            eprintln!(
                "Heat 2D {case_name} P{order}: L2={:.3}, H1={:.3}",
                r_hl2.asymptotic_slope, r_hh1.asymptotic_slope
            );
        }
    }
}

/// Perturbs interior tet edges smoothly `delta_x = c * |e|^2 * bump` (non-QP `CurvedTet`)
/// keeping boundary facets flat and det(J) strictly positive.
fn add_smooth_interior_curvature(mesh: &SimplexMesh, _h: f64) -> SimplexMesh {
    let mut bnd_edges = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
            bnd_edges.insert([a.min(b), a.max(b)]);
        }
    }
    let mut en = Vec::new();
    for &(va, vb) in mesh.edges() {
        let key = [va.min(vb), va.max(vb)];
        if bnd_edges.contains(&key) {
            continue;
        }
        let pa = mesh.point(va as usize);
        let pb = mesh.point(vb as usize);
        let len2 = (pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2) + (pb[2] - pa[2]).powi(2);
        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
        let bump = 0.25 * len2 * (PI * mid[0]).sin() * (PI * mid[1]).sin() * (PI * mid[2]).sin();
        en.push((key, [mid[0] + bump, mid[1] + 0.5 * bump, mid[2] - 0.5 * bump]));
    }
    let mut out = mesh.clone();
    out.set_edge_nodes(&en);
    assert!(out.has_curved() && !out.has_qp());
    out
}

#[test]
fn electro_and_heat_3d_mms_straight_and_curved_tet_iso_and_rotated_aniso() {
    let aniso_t = Tensor::Diagonal { x: 1.4, y: 2.3, z: 3.5 };
    let euler = [35.0, 48.0, 22.0];
    let a3_rel = tensor_global(&aniso_t, euler);
    let iso3_rel = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];

    let ns_p1 = [4usize, 8, 16];
    let ns_p2 = [4usize, 6, 9];

    // 3D Heat + Electro P1 & P2 on straight tets (iso + rotated aniso)
    for (case_name, is_aniso) in [("iso", false), ("rot_aniso", true)] {
        for &order in &[1usize, 2] {
            let ns = if order == 1 { &ns_p1[..] } else { &ns_p2[..] };
            let hs: Vec<f64> = ns.iter().map(|&n| 1.0 / n as f64).collect();
            let robin_alpha_e = 1.5 * EPS0;
            let mut a_e = if is_aniso { a3_rel } else { iso3_rel };
            a_e.iter_mut().for_each(|r| r.iter_mut().for_each(|v| *v *= EPS0));
            let src_e = mms_3d_tensor_source(a_e, robin_alpha_e, &[6]);
            let (mut l2_e, mut h1_e) = (Vec::new(), Vec::new());

            let robin_alpha_h = 1.5;
            let a_h = if is_aniso { a3_rel } else { iso3_rel };
            let src_h = mms_3d_tensor_source(a_h, robin_alpha_h, &[6]);
            let (mut l2_h, mut h1_h) = (Vec::new(), Vec::new());

            for &n in ns {
                let mesh = femcore::testkit::block(1.0, 1.0, 1.0, n, n, n);
                let mut pe = Problem::new(Study::Electrostatic, order);
                let mut ph = Problem::new(Study::HeatStationary, order);
                if is_aniso {
                    pe.materials.aniso.insert(
                        mesh.cell_labels[0],
                        Anisotropy { permittivity: Some(aniso_t.clone()), euler_deg: euler, ..Default::default() },
                    );
                    ph.materials.aniso.insert(
                        mesh.cell_labels[0],
                        Anisotropy { conductivity: Some(aniso_t.clone()), euler_deg: euler, ..Default::default() },
                    );
                } else {
                    pe.materials.default.insert("eps_r".into(), 2.0);
                    ph.materials.default.insert("k".into(), 2.0);
                }
                pe.bcs = vec![
                    Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Potential(0.0) },
                    Bc { labels: vec![5], kind: BcKind::SurfaceCharge(0.0) },
                ];
                let re =
                    solve_electrostatics_with_source_and_robin::<3, _>(&mesh, &pe, &src_e, &[(vec![6], robin_alpha_e)])
                        .unwrap();
                let (el2, eh1) = scalar_mms_norms_3d(&mesh, order, &re.solutions[0].field("V").unwrap().values);
                l2_e.push(el2);
                h1_e.push(eh1);

                ph.bcs = vec![
                    Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Temperature(0.0) },
                    Bc { labels: vec![5], kind: BcKind::Flux(0.0) },
                    Bc { labels: vec![6], kind: BcKind::Convection { h: robin_alpha_h, text: 0.0 } },
                ];
                let rh = solve_heat_with_source::<3, _>(&mesh, &ph, &src_h).unwrap();
                let (hl2, hh1) = scalar_mms_norms_3d(&mesh, order, &rh.solutions[0].field("T").unwrap().values);
                l2_h.push(hl2);
                h1_h.push(hh1);
            }
            let (l2_exp, l2_band) = if order == 1 { (2.0, (1.85, 2.25)) } else { (3.0, (2.80, 3.25)) };
            let (h1_exp, h1_band) = if order == 1 { (1.0, (0.90, 1.15)) } else { (2.0, (1.85, 2.20)) };
            assert_convergence_rate(&hs, &l2_e, l2_exp, l2_band);
            assert_convergence_rate(&hs, &h1_e, h1_exp, h1_band);
            assert_convergence_rate(&hs, &l2_h, l2_exp, l2_band);
            assert_convergence_rate(&hs, &h1_h, h1_exp, h1_band);
            eprintln!("3D Straight {case_name} P{order} Electro & Heat passed");
        }
    }

    // 3D Electro P2 on non-QP CurvedTet meshes (rotated anisotropic + Dirichlet/Neumann/Robin)
    let robin_alpha_c = 1.5 * EPS0;
    let mut a_c = a3_rel;
    a_c.iter_mut().for_each(|r| r.iter_mut().for_each(|v| *v *= EPS0));
    let src_c = mms_3d_tensor_source(a_c, robin_alpha_c, &[6]);
    let (mut l2_c, mut h1_c) = (Vec::new(), Vec::new());
    let hs_p2: Vec<f64> = ns_p2.iter().map(|&n| 1.0 / n as f64).collect();
    for &n in &ns_p2 {
        let base = femcore::testkit::block(1.0, 1.0, 1.0, n, n, n);
        let mesh = add_smooth_interior_curvature(&base, 1.0 / n as f64);
        let mut pe = Problem::new(Study::Electrostatic, 2);
        pe.materials.aniso.insert(
            mesh.cell_labels[0],
            Anisotropy { permittivity: Some(aniso_t.clone()), euler_deg: euler, ..Default::default() },
        );
        pe.bcs = vec![
            Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Potential(0.0) },
            Bc { labels: vec![5], kind: BcKind::SurfaceCharge(0.0) },
        ];
        let re = solve_electrostatics_with_source_and_robin::<3, _>(&mesh, &pe, &src_c, &[(vec![6], robin_alpha_c)])
            .unwrap();
        let (el2, eh1) = scalar_mms_norms_3d(&mesh, 2, &re.solutions[0].field("V").unwrap().values);
        l2_c.push(el2);
        h1_c.push(eh1);
    }
    let rc_l2 = assert_convergence_rate(&hs_p2, &l2_c, 3.0, (2.80, 3.25));
    let rc_h1 = assert_convergence_rate(&hs_p2, &h1_c, 2.0, (1.85, 2.20));
    eprintln!(
        "3D CurvedTet (non-QP) P2 rot_aniso Electro: L2={:.3}, H1={:.3}",
        rc_l2.asymptotic_slope, rc_h1.asymptotic_slope
    );
}

#[test]
fn magneto_mixed_3d_mms_nedelec_p1_hcurl_convergence() {
    // A*(x,y,z) = curl(sin(pi x) sin(pi y) sin(pi z) e_z)
    // J*(x,y,z) = curl(mu^-1 curl(A*)) = (3 pi^2 / (mu0 mu_r)) A*(x,y,z)
    let mu_r_val = 2.0_f64;
    let a_star = |x: [f64; 3]| -> [f64; 3] {
        [
            PI * (PI * x[0]).sin() * (PI * x[1]).cos() * (PI * x[2]).sin(),
            -PI * (PI * x[0]).cos() * (PI * x[1]).sin() * (PI * x[2]).sin(),
            0.0,
        ]
    };
    let curl_a_star = |x: [f64; 3]| -> [f64; 3] {
        [
            PI * PI * (PI * x[0]).cos() * (PI * x[1]).sin() * (PI * x[2]).cos(),
            PI * PI * (PI * x[0]).sin() * (PI * x[1]).cos() * (PI * x[2]).cos(),
            2.0 * PI * PI * (PI * x[0]).sin() * (PI * x[1]).sin() * (PI * x[2]).sin(),
        ]
    };
    let source = FnSource {
        volume_rhs: move |x: [f64; 3], _t: f64, _dom: u32| -> [f64; 3] {
            let a = a_star(x);
            let c = 3.0 * PI * PI / (MU0 * mu_r_val);
            [c * a[0], c * a[1], c * a[2]]
        },
        dirichlet: move |x: [f64; 3], _t: f64, _bnd: u32| -> [f64; 3] { a_star(x) },
        neumann_robin: |_x: [f64; 3], _n: [f64; 3], _t: f64, _bnd: u32| -> [f64; 3] { [0.0; 3] },
    };

    let ns = [6usize, 9, 14];
    let hs: Vec<f64> = ns.iter().map(|&n| 1.0 / n as f64).collect();
    let (mut err_a_l2, mut err_curl_l2) = (Vec::new(), Vec::new());
    let q64 = quad_rule_collapsed(4);

    for &n in &ns {
        let mesh = femcore::testkit::block(1.0, 1.0, 1.0, n, n, n);
        let mut spec = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
        spec.materials.default.insert("mu_r".into(), mu_r_val);
        spec.bcs = vec![Bc { labels: vec![1, 2, 3, 4, 5, 6], kind: BcKind::Ground }];
        let (_res, a_edges) = solve_mixed_3d_edges_with_source(&mesh, &spec, &source).unwrap();

        let geo = geometry(&mesh);
        let (mut l2a, mut l2c) = (0.0_f64, 0.0_f64);
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let eids = mesh.cell_edge_ids(c);
            let pts = [
                mesh.point(v[0] as usize),
                mesh.point(v[1] as usize),
                mesh.point(v[2] as usize),
                mesh.point(v[3] as usize),
            ];
            let mut g = [[0.0; 3]; 4];
            for (b, row) in g.iter_mut().enumerate() {
                row.copy_from_slice(&geo.g[c * 12 + b * 3..c * 12 + b * 3 + 3]);
            }
            for q in 0..q64.w.len() {
                let lam = [q64.lam[q][0], q64.lam[q][1], q64.lam[q][2], q64.lam[q][3]];
                let mut xq = [0.0; 3];
                for b in 0..4 {
                    for k in 0..3 {
                        xq[k] += lam[b] * pts[b][k];
                    }
                }
                let (mut ah, mut bh) = ([0.0; 3], [0.0; 3]);
                for (s, &(i, j)) in TET_EDGES.iter().enumerate() {
                    let (a, b) = if v[i] < v[j] { (i, j) } else { (j, i) };
                    let w = [
                        lam[a] * g[b][0] - lam[b] * g[a][0],
                        lam[a] * g[b][1] - lam[b] * g[a][1],
                        lam[a] * g[b][2] - lam[b] * g[a][2],
                    ];
                    let cu = [
                        2.0 * (g[a][1] * g[b][2] - g[a][2] * g[b][1]),
                        2.0 * (g[a][2] * g[b][0] - g[a][0] * g[b][2]),
                        2.0 * (g[a][0] * g[b][1] - g[a][1] * g[b][0]),
                    ];
                    let coef = a_edges[eids[s] as usize];
                    for k in 0..3 {
                        ah[k] += coef * w[k];
                        bh[k] += coef * cu[k];
                    }
                }
                let ae = a_star(xq);
                let ce = curl_a_star(xq);
                let wq = q64.w[q] * geo.vol[c];
                l2a += wq * ((ah[0] - ae[0]).powi(2) + (ah[1] - ae[1]).powi(2) + (ah[2] - ae[2]).powi(2));
                l2c += wq * ((bh[0] - ce[0]).powi(2) + (bh[1] - ce[1]).powi(2) + (bh[2] - ce[2]).powi(2));
            }
        }
        err_a_l2.push(l2a.sqrt());
        err_curl_l2.push(l2c.sqrt());
    }
    let ra = assert_convergence_rate(&hs, &err_a_l2, 1.0, (0.90, 1.15));
    let rc = assert_convergence_rate(&hs, &err_curl_l2, 1.0, (0.90, 1.15));
    eprintln!("Magneto mixed_3d P1 H(curl): L2(A)={:.3}, L2(curl A)={:.3}", ra.asymptotic_slope, rc.asymptotic_slope);
}
