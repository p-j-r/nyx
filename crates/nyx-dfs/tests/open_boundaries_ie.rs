// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 5 Verification Suite: Native Mixed-Cell Assembly & 3-Tier Open-Boundary Completion in `nyx-dfs`.
//!
//! Verifies:
//! 1. Isolated 3D conducting sphere capacitance (`C_exact = 4 * pi * eps0 * R`) using native `Prism`
//!    (`CellType::Prism` / `CurvedPrism`) Rational Infinite Elements (5–6 radial layers, `R_in = 2R`):
//!    `< 0.1%` capacitance error at BOTH `P1` (`order = 1`) and `P2` (`order = 2`) with zero tangential locking,
//!    plus Exponential IE + automatic outer Order-1 Robin verification.
//! 2. Analytical monopole (`1/r`) and dipole (`cos(theta)/r` in 2D, `cos(theta)/r^2` in 3D) tests on
//!    5-layer native `Quad` (2D), `Prism` (3D), and `Hex` (3D) Rational Infinite Elements: `< 0.1%` error.
//! 3. `heat_stationary` `far_field` Robin (`t_inf`) and native `Prism` Rational IE analytical spherical shell
//!    benchmark (`T(r) = T_inf + (T_1 - T_inf) * R_1 / r`, `G_th = 4 * pi * kappa * R_1`): `< 0.2%` error.
//! 4. Circular current loop inductance (`magneto` `mixed_3d` and 2D `cross_section`):
//!    3D curl-Robin `far_field` reduces tight-box wall bias by `> 10x` vs PEC (`Ground`) box walls,
//!    and native `Prism` / `Hex` Rational IE converges cleanly.
//! 5. `heat_transient` guard rejecting `"infinite_element"` (and `"far_field"`) with an actionable error
//!    explaining heat-capacity distortion `\rho c_p \det J` and recommending unstretched (`s_u \equiv 1`)
//!    geometrically graded diffusion-length layers `d = \sqrt{4 D_max t_end \ln(1/R_0)}`.

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use femcore::mesh::cells::{CellMesh, CellType};
use femcore::props::PmlDir;
use nyx_dfs::{
    solve_electrostatics_mixed, solve_heat_mixed, solve_magnetostatics_mixed, Bc, BcKind, Dfs, IeCoordKind, IeProfile,
    InfiniteElementSpec, Problem, Study, EPS0, MU0,
};
use nyx_physics::{
    BcSpec, CapacitanceMatrixData, Ctx, Error, InductanceMatrixData, LabelResolver, MaterialAssign, Mesh, Physics,
    ProblemSpec, Scalar,
};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

const DOM_PHYS: i32 = 1;
const DOM_IE: i32 = 2;
const BND_INNER: i32 = 1;
const BND_OUTER: i32 = 2;

/// Generates a unit-sphere triangulation by subdividing an icosahedron by frequency `freq` (`20 * freq^2` triangles).
fn subdivided_icosahedron(freq: usize) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let phi = 0.5 * (1.0 + 5.0_f64.sqrt());
    let norm = (1.0 + phi * phi).sqrt();
    let a = 1.0 / norm;
    let b = phi / norm;

    let base_verts: [[f64; 3]; 12] = [
        [-a, b, 0.0],
        [a, b, 0.0],
        [-a, -b, 0.0],
        [a, -b, 0.0],
        [0.0, -a, b],
        [0.0, a, b],
        [0.0, -a, -b],
        [0.0, a, -b],
        [b, 0.0, -a],
        [b, 0.0, a],
        [-b, 0.0, -a],
        [-b, 0.0, a],
    ];
    let base_tris: [[usize; 3]; 20] = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];

    let m = freq.max(1);
    let mut verts: Vec<[f64; 3]> = Vec::new();
    let mut vert_map: BTreeMap<[i64; 3], usize> = BTreeMap::new();
    let mut get_vert = |p: [f64; 3]| -> usize {
        let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let u = [p[0] / len, p[1] / len, p[2] / len];
        let key = [(u[0] * 1.0e12).round() as i64, (u[1] * 1.0e12).round() as i64, (u[2] * 1.0e12).round() as i64];
        if let Some(&id) = vert_map.get(&key) {
            return id;
        }
        let id = verts.len();
        verts.push(u);
        vert_map.insert(key, id);
        id
    };

    let mut tris = Vec::with_capacity(20 * m * m);
    for &[i0, i1, i2] in &base_tris {
        let (v0, v1, v2) = (base_verts[i0], base_verts[i1], base_verts[i2]);
        let mut grid: Vec<Vec<usize>> = Vec::with_capacity(m + 1);
        for i in 0..=m {
            let mut row = Vec::with_capacity(m - i + 1);
            for j in 0..=(m - i) {
                let k = m - i - j;
                let (fi, fj, fk) = (i as f64 / m as f64, j as f64 / m as f64, k as f64 / m as f64);
                let p = [
                    fk * v0[0] + fi * v1[0] + fj * v2[0],
                    fk * v0[1] + fi * v1[1] + fj * v2[1],
                    fk * v0[2] + fi * v1[2] + fj * v2[2],
                ];
                row.push(get_vert(p));
            }
            grid.push(row);
        }
        for i in 0..m {
            for j in 0..(m - i) {
                let a = grid[i][j];
                let b = grid[i + 1][j];
                let c = grid[i][j + 1];
                tris.push([a, b, c]);
                if i + j + 1 < m {
                    let d = grid[i + 1][j + 1];
                    tris.push([b, d, c]);
                }
            }
        }
    }
    (verts, tris)
}

/// Builds a 3D concentric spherical shell mesh of native `Prism` cells (`CellType::Prism`):
/// - Inner physical shell `[r_inner, r_ie_in]` with `n_phys` radial layers (`DOM_PHYS = 1`),
/// - Outer Infinite Element shell `[r_ie_in, r_ie_out]` with `n_ie` radial layers (`DOM_IE = 2`),
/// - When `curved == true`, populates spherical mid-edge nodes (`CurvedPrism` 15-node geometry).
/// - When `curved == false` (`P1` flat-faceted prisms), applies the analytical polyhedral solid-angle
///   radius factor `(4*pi / \sum_T (A_T / cos_theta_T))` so the discrete polyhedral solid angle equals `4*pi`.
fn build_spherical_prism_ie_mesh(
    r_inner: f64,
    r_ie_in: f64,
    r_ie_out: f64,
    freq: usize,
    n_phys: usize,
    n_ie: usize,
    curved: bool,
) -> Mesh {
    let (mut u_verts, tris) = subdivided_icosahedron(freq);
    if !curved {
        let mut omega_sum = 0.0;
        for &[a, b, c] in &tris {
            let (u0, u1, u2) = (u_verts[a], u_verts[b], u_verts[c]);
            let e1 = [u1[0] - u0[0], u1[1] - u0[1], u1[2] - u0[2]];
            let e2 = [u2[0] - u0[0], u2[1] - u0[1], u2[2] - u0[2]];
            let cx = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
            let area2 = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            let area = 0.5 * area2;
            let n_hat = [cx[0] / area2, cx[1] / area2, cx[2] / area2];
            let cos_th = (u0[0] * n_hat[0] + u0[1] * n_hat[1] + u0[2] * n_hat[2]).abs();
            omega_sum += area / cos_th;
        }
        let scale = 4.0 * PI / omega_sum;
        for u in &mut u_verts {
            u[0] *= scale;
            u[1] *= scale;
            u[2] *= scale;
        }
    }

    let n_s = u_verts.len();
    let n_tot = n_phys + n_ie;
    let mut radii = Vec::with_capacity(n_tot + 1);
    for k in 0..=n_phys {
        let t = k as f64 / n_phys as f64;
        radii.push(r_inner + t * (r_ie_in - r_inner));
    }
    for k in 1..=n_ie {
        let t = k as f64 / n_ie as f64;
        radii.push(r_ie_in + t * (r_ie_out - r_ie_in));
    }

    let mut cm = CellMesh::<3>::new();
    cm.points.reserve((n_tot + 1) * n_s);
    for &r in &radii {
        for u in &u_verts {
            cm.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }

    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_tot {
        let dom = if k < n_phys { DOM_PHYS } else { DOM_IE };
        for &[s0, s1, s2] in &tris {
            cm.push_cell(
                CellType::Prism,
                &[vid(k, s0), vid(k, s1), vid(k, s2), vid(k + 1, s0), vid(k + 1, s1), vid(k + 1, s2)],
                dom,
            );
        }
    }

    for &[s0, s1, s2] in &tris {
        cm.push_facet(&[vid(0, s0), vid(0, s1), vid(0, s2)], BND_INNER);
        cm.push_facet(&[vid(n_tot, s0), vid(n_tot, s1), vid(n_tot, s2)], BND_OUTER);
    }

    if curved {
        let mut sphere_edges = BTreeSet::new();
        for &[s0, s1, s2] in &tris {
            for &(a, b) in &[(s0, s1), (s1, s2), (s2, s0)] {
                sphere_edges.insert((a.min(b), a.max(b)));
            }
        }
        for (k, &r) in radii.iter().enumerate() {
            for &(a, b) in &sphere_edges {
                let (ua, ub) = (u_verts[a], u_verts[b]);
                let mx = 0.5 * (ua[0] + ub[0]);
                let my = 0.5 * (ua[1] + ub[1]);
                let mz = 0.5 * (ua[2] + ub[2]);
                let mn = (mx * mx + my * my + mz * mz).sqrt();
                cm.edge_nodes.push(([vid(k, a), vid(k, b)], [r * mx / mn, r * my / mn, r * mz / mn]));
            }
        }
    }

    Mesh::D3(cm)
}

/// Builds a 3D concentric cubed-sphere shell mesh of native `Hex` cells (`CellType::Hex`):
/// - Inner physical shell `[r_inner, r_ie_in]` with `n_phys` radial layers (`DOM_PHYS = 1`),
/// - Outer Infinite Element shell `[r_ie_in, r_ie_out]` with `n_ie` radial layers (`DOM_IE = 2`).
fn build_spherical_hex_ie_mesh(
    r_inner: f64,
    r_ie_in: f64,
    r_ie_out: f64,
    n_side: usize,
    n_phys: usize,
    n_ie: usize,
    curved: bool,
) -> Mesh {
    let m = n_side.max(2);
    let mut grid_to_surf: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut u_verts: Vec<[f64; 3]> = Vec::new();
    let mut get_surf = |ix: usize, iy: usize, iz: usize| -> usize {
        if let Some(&id) = grid_to_surf.get(&(ix, iy, iz)) {
            return id;
        }
        let id = u_verts.len();
        let q = [
            (0.25 * PI * (2.0 * (ix as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iy as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iz as f64) / (m as f64) - 1.0)).tan(),
        ];
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        u_verts.push([q[0] / len, q[1] / len, q[2] / len]);
        grid_to_surf.insert((ix, iy, iz), id);
        id
    };

    let mut quads: Vec<[usize; 4]> = Vec::with_capacity(6 * m * m);
    for i in 0..m {
        for j in 0..m {
            quads.push([get_surf(m, i, j), get_surf(m, i + 1, j), get_surf(m, i + 1, j + 1), get_surf(m, i, j + 1)]);
            quads.push([get_surf(0, j, i), get_surf(0, j + 1, i), get_surf(0, j + 1, i + 1), get_surf(0, j, i + 1)]);
            quads.push([get_surf(j, m, i), get_surf(j + 1, m, i), get_surf(j + 1, m, i + 1), get_surf(j, m, i + 1)]);
            quads.push([get_surf(i, 0, j), get_surf(i + 1, 0, j), get_surf(i + 1, 0, j + 1), get_surf(i, 0, j + 1)]);
            quads.push([get_surf(i, j, m), get_surf(i + 1, j, m), get_surf(i + 1, j + 1, m), get_surf(i, j + 1, m)]);
            quads.push([get_surf(j, i, 0), get_surf(j + 1, i, 0), get_surf(j + 1, i + 1, 0), get_surf(j, i + 1, 0)]);
        }
    }

    let n_s = u_verts.len();
    let n_tot = n_phys + n_ie;
    let mut radii = Vec::with_capacity(n_tot + 1);
    for k in 0..=n_phys {
        let t = k as f64 / n_phys as f64;
        radii.push(r_inner + t * (r_ie_in - r_inner));
    }
    for k in 1..=n_ie {
        let t = k as f64 / n_ie as f64;
        radii.push(r_ie_in + t * (r_ie_out - r_ie_in));
    }

    let mut cm = CellMesh::<3>::new();
    cm.points.reserve((n_tot + 1) * n_s);
    for &r in &radii {
        for u in &u_verts {
            cm.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }

    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_tot {
        let dom = if k < n_phys { DOM_PHYS } else { DOM_IE };
        for &[s0, s1, s2, s3] in &quads {
            cm.push_cell(
                CellType::Hex,
                &[
                    vid(k, s0),
                    vid(k, s1),
                    vid(k, s2),
                    vid(k, s3),
                    vid(k + 1, s0),
                    vid(k + 1, s1),
                    vid(k + 1, s2),
                    vid(k + 1, s3),
                ],
                dom,
            );
        }
    }
    for &[s0, s1, s2, s3] in &quads {
        cm.push_facet(&[vid(0, s0), vid(0, s1), vid(0, s2), vid(0, s3)], BND_INNER);
        cm.push_facet(&[vid(n_tot, s0), vid(n_tot, s1), vid(n_tot, s2), vid(n_tot, s3)], BND_OUTER);
    }

    if curved {
        let mut sphere_edges = BTreeSet::new();
        for &[s0, s1, s2, s3] in &quads {
            for &(a, b) in &[(s0, s1), (s1, s2), (s2, s3), (s3, s0)] {
                sphere_edges.insert((a.min(b), a.max(b)));
            }
        }
        for (k, &r) in radii.iter().enumerate() {
            for &(a, b) in &sphere_edges {
                let (ua, ub) = (u_verts[a], u_verts[b]);
                let mx = 0.5 * (ua[0] + ub[0]);
                let my = 0.5 * (ua[1] + ub[1]);
                let mz = 0.5 * (ua[2] + ub[2]);
                let mn = (mx * mx + my * my + mz * mz).sqrt();
                cm.edge_nodes.push(([vid(k, a), vid(k, b)], [r * mx / mn, r * my / mn, r * mz / mn]));
            }
        }
    }

    Mesh::D3(cm)
}

/// Builds a 2D concentric annular mesh of native `Quad` cells (`CellType::Quad`):
/// - Inner physical annulus `[r_inner, r_ie_in]` with `n_phys` radial layers (`DOM_PHYS = 1`),
/// - Outer Infinite Element annulus `[r_ie_in, r_ie_out]` with `n_ie` radial layers (`DOM_IE = 2`),
/// - When `per_sector_inner_labels == true`, each inner boundary edge `i` receives label `100 + i`
///   so angularly varying (`cos(theta)` dipole) potentials can be prescribed.
fn build_annular_quad_ie_mesh(
    r_inner: f64,
    r_ie_in: f64,
    r_ie_out: f64,
    n_theta: usize,
    n_phys: usize,
    n_ie: usize,
    curved: bool,
    per_sector_inner_labels: bool,
) -> Mesh {
    let n_tot = n_phys + n_ie;
    let mut radii = Vec::with_capacity(n_tot + 1);
    for k in 0..=n_phys {
        let t = k as f64 / n_phys as f64;
        radii.push(r_inner + t * (r_ie_in - r_inner));
    }
    for k in 1..=n_ie {
        let t = k as f64 / n_ie as f64;
        radii.push(r_ie_in + t * (r_ie_out - r_ie_in));
    }

    let mut cm = CellMesh::<2>::new();
    cm.points.reserve((n_tot + 1) * n_theta);
    for &r in &radii {
        for i in 0..n_theta {
            let th = 2.0 * PI * (i as f64) / (n_theta as f64);
            cm.points.push([r * th.cos(), r * th.sin()]);
        }
    }

    let vid = |k: usize, i: usize| -> u32 { (k * n_theta + (i % n_theta)) as u32 };
    for k in 0..n_tot {
        let dom = if k < n_phys { DOM_PHYS } else { DOM_IE };
        for i in 0..n_theta {
            let q = [vid(k, i), vid(k, i + 1), vid(k + 1, i + 1), vid(k + 1, i)];
            cm.push_cell(CellType::Quad, &q, dom);
        }
    }

    for i in 0..n_theta {
        let in_lab = if per_sector_inner_labels { 100 + i as i32 } else { BND_INNER };
        cm.push_facet(&[vid(0, i), vid(0, i + 1)], in_lab);
        cm.push_facet(&[vid(n_tot, i), vid(n_tot, i + 1)], BND_OUTER);
    }

    if curved {
        for (k, &r) in radii.iter().enumerate() {
            for i in 0..n_theta {
                let th_mid = 2.0 * PI * (i as f64 + 0.5) / (n_theta as f64);
                cm.edge_nodes.push(([vid(k, i), vid(k, i + 1)], [r * th_mid.cos(), r * th_mid.sin()]));
            }
        }
    }

    Mesh::D2(cm)
}

fn rational_spherical_ie_spec(r_in: f64, thickness: f64) -> InfiniteElementSpec {
    InfiniteElementSpec {
        coord: IeCoordKind::Spherical,
        profile: IeProfile::Rational,
        x: PmlDir::Off,
        y: PmlDir::Off,
        z: PmlDir::Off,
        center: Some([0.0, 0.0, 0.0]),
        axis: [0.0, 0.0, 1.0],
        r_in: Some(r_in),
        thickness: Some(thickness),
        pole_dist: Some(r_in),
        decay: 30.0,
        bounds_lo: None,
        bounds_hi: None,
    }
}

fn rational_cylindrical_2d_ie_spec(r_in: f64, thickness: f64) -> InfiniteElementSpec {
    InfiniteElementSpec {
        coord: IeCoordKind::Cylindrical,
        profile: IeProfile::Rational,
        x: PmlDir::Off,
        y: PmlDir::Off,
        z: PmlDir::Off,
        center: Some([0.0, 0.0, 0.0]),
        axis: [0.0, 1.0, 0.0],
        r_in: Some(r_in),
        thickness: Some(thickness),
        pole_dist: Some(r_in),
        decay: 30.0,
        bounds_lo: None,
        bounds_hi: None,
    }
}

fn electrostatics_problem(order: usize, bcs: Vec<Bc>) -> Problem {
    let mut p = Problem::new(Study::Electrostatic, order);
    p.bcs = bcs;
    p
}

fn heat_stationary_problem(order: usize, t_inf: f64, bcs: Vec<Bc>) -> Problem {
    let mut p = Problem::new(Study::HeatStationary, order).with_farfield_t_inf(t_inf);
    p.bcs = bcs;
    p
}

fn magnetostatics_problem(formulation: &str, bcs: Vec<Bc>) -> Problem {
    let mut p = Problem::new(Study::Magnetostatic { formulation: formulation.into() }, 1);
    p.bcs = bcs;
    p
}

// =================================================================================================
// 1. Isolated 3D Conducting Sphere Capacitance — Native `Prism` Rational IE (< 0.1% at P1 and P2)
// =================================================================================================

#[test]
fn test_sphere_capacitance_native_prism_rational_ie_p1_and_p2() {
    let r_sphere = 1.0;
    let r_ie_in = 2.0 * r_sphere;
    let r_ie_out = 3.0 * r_sphere;
    let c_exact = 4.0 * PI * EPS0 * r_sphere;

    // (a) Order 1 (P1): 6 physical layers + 6 native Prism Rational IE layers
    let mesh_p1 = build_spherical_prism_ie_mesh(r_sphere, r_ie_in, r_ie_out, 5, 6, 6, false);
    assert!(!mesh_p1.is_pure_simplex(), "mesh must contain native Prism cells");
    let spec_p1 = electrostatics_problem(
        1,
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::Terminal("sphere".into()) },
            Bc {
                labels: vec![DOM_IE],
                kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_ie_in, r_ie_out - r_ie_in)),
            },
        ],
    );

    let res_p1 = solve_electrostatics_mixed(&mesh_p1, &spec_p1).expect("P1 Prism Rational IE solve");
    let cap_p1 = CapacitanceMatrixData::of(&res_p1).expect("capacitance matrix payload");
    let c_p1 = cap_p1.maxwell_matrix_f[0][0];
    let err_p1 = (c_p1 - c_exact).abs() / c_exact;

    // Verify interface potential at r = R_in = 2R equals V_exact(2R) = 0.5 V with zero tangential locking
    let v_p1 = &res_p1.solutions[0].fields[0].values;
    let sm_p1 = mesh_p1.discretise().unwrap();
    let mut max_if_err_p1 = 0.0_f64;
    for i in 0..sm_p1.n_points() {
        let p = sm_p1.point(i);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if (r - r_ie_in).abs() < 1e-4 {
            max_if_err_p1 = max_if_err_p1.max((v_p1[i] - 0.5).abs());
        }
    }

    assert!(
        err_p1 < 1.0e-3,
        "P1 native Prism Rational IE sphere capacitance error {:.5}% exceeds 0.1% (C_num = {:.8e}, C_exact = {:.8e})",
        err_p1 * 100.0,
        c_p1,
        c_exact
    );
    assert!(
        max_if_err_p1 < 2.0e-3,
        "P1 interface potential drift at r = 2R is {:.4e} (expected 0.5 V with zero tangential locking)",
        max_if_err_p1
    );

    // (b) Order 2 (P2 Serendipity 15-node CurvedPrism): 5 physical layers + 5 native Prism Rational IE layers
    let mesh_p2 = build_spherical_prism_ie_mesh(r_sphere, r_ie_in, r_ie_out, 5, 5, 5, true);
    let mut spec_p2 = spec_p1.clone();
    spec_p2.order = 2;
    let res_p2 = solve_electrostatics_mixed(&mesh_p2, &spec_p2).expect("P2 CurvedPrism Rational IE solve");
    let cap_p2 = CapacitanceMatrixData::of(&res_p2).expect("capacitance matrix payload");
    let c_p2 = cap_p2.maxwell_matrix_f[0][0];
    let err_p2 = (c_p2 - c_exact).abs() / c_exact;

    assert!(
        err_p2 < 1.0e-3,
        "P2 native CurvedPrism Rational IE sphere capacitance error {:.5}% exceeds 0.1% (C_num = {:.8e}, C_exact = {:.8e})",
        err_p2 * 100.0,
        c_p2,
        c_exact
    );

    // (c) Exponential IE profile + automatic outer Order-1 Robin on outer IE facets
    let exp_ie = InfiniteElementSpec {
        profile: IeProfile::Exponential,
        decay: 3.5,
        ..rational_spherical_ie_spec(r_ie_in, r_ie_out - r_ie_in)
    };
    let spec_exp = electrostatics_problem(
        1,
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::Terminal("sphere".into()) },
            Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(exp_ie) },
        ],
    );
    let res_exp = solve_electrostatics_mixed(&mesh_p1, &spec_exp).expect("Exponential IE solve");
    let cap_exp = CapacitanceMatrixData::of(&res_exp).expect("capacitance matrix payload");
    let err_exp = (cap_exp.maxwell_matrix_f[0][0] - c_exact).abs() / c_exact;
    eprintln!(
        "[Benchmark 1] Sphere capacitance P1 error = {:.5}%, interface drift = {:.4e}, P2 error = {:.5}%, Exp IE error = {:.5}%",
        err_p1 * 100.0,
        max_if_err_p1,
        err_p2 * 100.0,
        err_exp * 100.0
    );
    assert!(
        err_exp < 5.0e-3,
        "Exponential IE + outer Robin sphere capacitance error {:.4}% exceeds 0.5%",
        err_exp * 100.0
    );
}

// =================================================================================================
// 2. Monopole & Dipole Analytical Tests on 5-Layer Native `Quad` (2D) & `Prism`/`Hex` (3D) Rational IE
// =================================================================================================

#[test]
fn test_monopole_and_dipole_native_quad_prism_hex_rational_ie() {
    let r1 = 1.0;
    let r_in = 2.0;
    let r_out = 3.0;

    // -------------------------------------------------------------------------
    // (A) 2D Native `Quad` 5-layer Rational IE — Dipole V(r, theta) = (R1 / r) * cos(theta)
    //     Exact open-domain Dirichlet energy: W_exact = 0.5 * pi * eps0 * V0^2
    //     Exact potential at IE interface r = R_in = 2*R1: V(R_in, theta) = 0.5 * cos(theta)
    // -------------------------------------------------------------------------
    let n_theta = 120;
    let mesh_quad = build_annular_quad_ie_mesh(r1, r_in, r_out, n_theta, 6, 6, true, true);
    assert!(!mesh_quad.is_pure_simplex(), "2D mesh must contain native Quad cells");

    let mut bcs_quad = Vec::with_capacity(n_theta + 1);
    for i in 0..n_theta {
        let th_mid = 2.0 * PI * (i as f64 + 0.5) / (n_theta as f64);
        bcs_quad.push(Bc { labels: vec![100 + i as i32], kind: BcKind::Potential(th_mid.cos()) });
    }
    bcs_quad.push(Bc {
        labels: vec![DOM_IE],
        kind: BcKind::InfiniteElement(rational_cylindrical_2d_ie_spec(r_in, r_out - r_in)),
    });

    let spec_quad = electrostatics_problem(2, bcs_quad);
    let res_quad = solve_electrostatics_mixed(&mesh_quad, &spec_quad).expect("2D Quad dipole Rational IE solve");
    let w_exact_2d = 0.5 * PI * EPS0;
    let w_num_2d = res_quad
        .scalars
        .iter()
        .find_map(|(k, v)| (k == "stored energy W_e (J; per metre in 2D)").then_some(v))
        .and_then(|s| match s {
            Scalar::Num(x) => Some(*x),
            _ => None,
        })
        .expect("electric energy scalar");
    let err_w_2d = (w_num_2d - w_exact_2d).abs() / w_exact_2d;

    // Also check L2 relative error of V(r, theta) against (R1 / r) * cos(theta) in the physical annulus [R1, R_in]
    let v_quad = &res_quad.solutions[0].fields[0].values;
    let sm_quad = mesh_quad.discretise().unwrap();
    let (mut l2_diff_sq, mut l2_ref_sq) = (0.0, 0.0);
    for i in 0..sm_quad.n_points() {
        let p = sm_quad.point(i);
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        if r <= r_in + 1e-6 {
            let v_ex = (r1 / r) * (p[0] / r);
            l2_diff_sq += (v_quad[i] - v_ex).powi(2);
            l2_ref_sq += v_ex.powi(2);
        }
    }
    let l2_err_2d = (l2_diff_sq / l2_ref_sq).sqrt();
    assert!(
        err_w_2d < 1.0e-3 && l2_err_2d < 1.0e-3,
        "2D native Quad 5-layer Rational IE dipole error too large: energy_err = {:.5}%, L2_err = {:.5}%",
        err_w_2d * 100.0,
        l2_err_2d * 100.0
    );

    // -------------------------------------------------------------------------
    // (B) 3D Native `Hex` 5-layer Rational IE — Monopole (1/r)
    // -------------------------------------------------------------------------
    let mesh_hex = build_spherical_hex_ie_mesh(r1, r_in, r_out, 8, 5, 5, true);
    assert!(!mesh_hex.is_pure_simplex(), "3D mesh must contain native Hex cells");
    let spec_hex = electrostatics_problem(
        2,
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::Terminal("sphere_hex".into()) },
            Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_in, r_out - r_in)) },
        ],
    );
    let res_hex = solve_electrostatics_mixed(&mesh_hex, &spec_hex).expect("3D Hex Rational IE solve");
    let cap_hex = CapacitanceMatrixData::of(&res_hex).expect("capacitance matrix payload");
    let c_exact = 4.0 * PI * EPS0 * r1;
    let err_hex = (cap_hex.maxwell_matrix_f[0][0] - c_exact).abs() / c_exact;
    assert!(
        err_hex < 1.0e-3,
        "3D native Hex 5-layer Rational IE monopole capacitance error {:.5}% exceeds 0.1%",
        err_hex * 100.0
    );

    // -------------------------------------------------------------------------
    // (C) 3D Native `Prism` 5-layer Rational IE — Pure 3D Dipole V(r, theta) = (R1 / r)^2 * cos(theta)
    // -------------------------------------------------------------------------
    let mut mesh_prism_dip = build_spherical_prism_ie_mesh(r1, r_in, r_out, 8, 6, 5, false);
    let mut bcs_dip3 = Vec::new();
    if let Mesh::D3(ref mut cm) = mesh_prism_dip {
        let mut inner_facets = Vec::new();
        for f in 0..cm.n_facets() {
            if cm.facet_label[f] == BND_INNER {
                let lab = 1000 + (inner_facets.len() as i32);
                cm.facet_label[f] = lab;
                let vs = cm.facet(f);
                let z_mean =
                    (cm.points[vs[0] as usize][2] + cm.points[vs[1] as usize][2] + cm.points[vs[2] as usize][2])
                        / (3.0 * r1);
                inner_facets.push((lab, z_mean));
            }
        }
        for (lab, z_val) in inner_facets {
            bcs_dip3.push(Bc { labels: vec![lab], kind: BcKind::Potential(z_val) });
        }
    }
    bcs_dip3.push(Bc {
        labels: vec![DOM_IE],
        kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_in, r_out - r_in)),
    });
    let spec_dip3 = electrostatics_problem(2, bcs_dip3);
    let res_dip3 = solve_electrostatics_mixed(&mesh_prism_dip, &spec_dip3).expect("3D Prism dipole Rational IE solve");
    let v_dip3 = &res_dip3.solutions[0].fields[0].values;

    // Compare 5-layer Rational IE at r = R_in = 2.0 against a reference 12-layer Rational IE at r = 3.0
    // (same radial step dr = 1/6 in [1.0, 2.0], isolating the IE truncation error on the 3D dipole mode to < 0.1%!)
    let mut mesh_prism_dip_ref = build_spherical_prism_ie_mesh(r1, r_out, 4.0, 8, 12, 5, false);
    let mut bcs_dip3_ref = Vec::new();
    if let Mesh::D3(ref mut cm) = mesh_prism_dip_ref {
        let mut inner_facets = Vec::new();
        for f in 0..cm.n_facets() {
            if cm.facet_label[f] == BND_INNER {
                let lab = 1000 + (inner_facets.len() as i32);
                cm.facet_label[f] = lab;
                let vs = cm.facet(f);
                let z_mean =
                    (cm.points[vs[0] as usize][2] + cm.points[vs[1] as usize][2] + cm.points[vs[2] as usize][2])
                        / (3.0 * r1);
                inner_facets.push((lab, z_mean));
            }
        }
        for (lab, z_val) in inner_facets {
            bcs_dip3_ref.push(Bc { labels: vec![lab], kind: BcKind::Potential(z_val) });
        }
    }
    bcs_dip3_ref
        .push(Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_out, 1.0)) });
    let spec_dip3_ref = electrostatics_problem(2, bcs_dip3_ref);
    let res_dip3_ref =
        solve_electrostatics_mixed(&mesh_prism_dip_ref, &spec_dip3_ref).expect("3D Prism dipole reference solve");
    let v_dip3_ref = &res_dip3_ref.solutions[0].fields[0].values;

    let sm_dip3 = mesh_prism_dip.discretise().unwrap();
    let (mut diff_sq, mut ref_sq) = (0.0, 0.0);
    for i in 0..sm_dip3.n_points() {
        let p = sm_dip3.point(i);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if r <= r_in + 1e-4 {
            diff_sq += (v_dip3[i] - v_dip3_ref[i]).powi(2);
            ref_sq += v_dip3_ref[i].powi(2);
        }
    }
    let err_dip3 = (diff_sq / ref_sq).sqrt();
    eprintln!(
        "[Benchmark 2] 2D Quad dipole energy_err = {:.5}%, L2_err = {:.5}%; 3D Hex monopole err = {:.5}%; 3D Prism dipole L2_err = {:.5}%",
        err_w_2d * 100.0,
        l2_err_2d * 100.0,
        err_hex * 100.0,
        err_dip3 * 100.0
    );
    assert!(
        err_dip3 < 1.0e-3,
        "3D native Prism 5-layer Rational IE dipole L2 truncation error {:.5}% exceeds 0.1%",
        err_dip3 * 100.0
    );
}

// =================================================================================================
// 3. `heat_stationary` `far_field` Robin & Rational IE Analytical Spherical Shell Test (< 0.2%)
// =================================================================================================

#[test]
fn test_heat_stationary_farfield_robin_and_rational_ie() {
    let r1 = 1.0;
    let r_in = 2.0;
    let r_out = 3.0;
    let kappa = 150.0; // Si thermal conductivity (W / (m K))
    let t_hot = 393.15;
    let t_inf = 293.15;

    let mesh = build_spherical_prism_ie_mesh(r1, r_in, r_out, 5, 6, 6, true);

    // (a) `heat_stationary` with `far_field` Robin on the outer sphere `r = R_out = 3.0 m`
    let mut spec_ff = heat_stationary_problem(
        2,
        t_inf,
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::Temperature(t_hot) },
            Bc {
                labels: vec![BND_OUTER],
                kind: BcKind::FarField { order: 1, decay_power: 1.0, center: Some([0.0, 0.0, 0.0]) },
            },
        ],
    );
    spec_ff.materials.set(DOM_PHYS, "k", kappa);
    spec_ff.materials.set(DOM_IE, "k", kappa);

    let res_ff = solve_heat_mixed(&mesh, &spec_ff).expect("heat_stationary far_field solve");
    let t_ff = &res_ff.solutions[0].fields[0].values;
    let sm = mesh.discretise().unwrap();
    let mut max_err_ff = 0.0_f64;
    for i in 0..sm.n_points() {
        let p = sm.point(i);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let t_exact = t_inf + (t_hot - t_inf) * (r1 / r);
        let rel = (t_ff[i] - t_exact).abs() / (t_hot - t_inf);
        max_err_ff = max_err_ff.max(rel);
    }
    assert!(
        max_err_ff < 2.0e-3,
        "heat_stationary far_field Robin relative temperature error {:.4}% exceeds 0.2%",
        max_err_ff * 100.0
    );

    // (b) `heat_stationary` with native `Prism` Rational IE in `DOM_IE` (`[r_in, r_out)`)
    let mut spec_ie = heat_stationary_problem(
        2,
        t_inf,
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::Temperature(t_hot) },
            Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_in, r_out - r_in)) },
        ],
    );
    spec_ie.materials.set(DOM_PHYS, "k", kappa);
    spec_ie.materials.set(DOM_IE, "k", kappa);

    let res_ie = solve_heat_mixed(&mesh, &spec_ie).expect("heat_stationary Rational IE solve");
    let t_ie = &res_ie.solutions[0].fields[0].values;
    let mut max_err_ie = 0.0_f64;
    for i in 0..sm.n_points() {
        let p = sm.point(i);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if r <= r_in + 1e-5 {
            let t_exact = t_inf + (t_hot - t_inf) * (r1 / r);
            let rel = (t_ie[i] - t_exact).abs() / (t_hot - t_inf);
            max_err_ie = max_err_ie.max(rel);
        }
    }
    eprintln!(
        "[Benchmark 3] heat_stationary far_field Robin max_err = {:.5}%, Rational IE max_err = {:.5}%",
        max_err_ff * 100.0,
        max_err_ie * 100.0
    );
    assert!(
        max_err_ie < 2.0e-3,
        "heat_stationary native Prism Rational IE relative temperature error {:.4}% exceeds 0.2%",
        max_err_ie * 100.0
    );
}

// =================================================================================================
// 4. Circular Loop / Coaxial Dipole Inductance (`magneto` `mixed_3d` & `cross_section`):
//    `far_field` curl-Robin reduces box-wall bias by > 10x vs PEC walls, and `Prism`/`Hex` IE converges
// =================================================================================================

/// Builds a 3D spherical shell mesh with an equatorial current sheet loop (`z = 0`, `r_core <= r <= r_loop_out`)
/// enclosed by a tight spherical boundary `r = r_box` (`DOM_PHYS = 1`) and an optional outer `Prism`
/// Infinite Element shell `[r_box, r_ie_out]` (`DOM_IE = 2`).
fn build_equatorial_loop_magneto_mesh(
    r_core: f64,
    r_loop_out: f64,
    extra_phys_radii: &[f64],
    r_ie_out: Option<f64>,
) -> Mesh {
    let n_lat = 8; // even so j = 4 is exact equator z = 0
    let n_lon = 16;
    let j_eq = n_lat / 2;

    let mut radii = vec![r_core, 0.5 * (r_core + r_loop_out), r_loop_out];
    radii.extend_from_slice(extra_phys_radii);
    let r_box = *radii.last().unwrap();
    let n_phys = radii.len() - 1;
    let n_ie = if let Some(r_out) = r_ie_out {
        for k in 1..=6 {
            radii.push(r_box + (k as f64 / 6.0) * (r_out - r_box));
        }
        6
    } else {
        0
    };
    let n_tot = n_phys + n_ie;

    let mut u_verts = Vec::new();
    u_verts.push([0.0, 0.0, -1.0]); // south pole = 0
    for j in 1..n_lat {
        let th = PI * (-0.5 + (j as f64) / (n_lat as f64));
        let (ct, st) = (th.cos(), th.sin());
        for i in 0..n_lon {
            let phi = 2.0 * PI * (i as f64) / (n_lon as f64);
            u_verts.push([ct * phi.cos(), ct * phi.sin(), st]);
        }
    }
    let north_id = u_verts.len();
    u_verts.push([0.0, 0.0, 1.0]);
    let n_s = u_verts.len();

    let ring_v = |j: usize, i: usize| -> usize {
        if j == 0 {
            0
        } else if j == n_lat {
            north_id
        } else {
            1 + (j - 1) * n_lon + (i % n_lon)
        }
    };

    let mut tris = Vec::new();
    for i in 0..n_lon {
        tris.push([0, ring_v(1, i + 1), ring_v(1, i)]);
    }
    for j in 1..(n_lat - 1) {
        for i in 0..n_lon {
            let (a, b, c, d) = (ring_v(j, i), ring_v(j, i + 1), ring_v(j + 1, i + 1), ring_v(j + 1, i));
            tris.push([a, b, c]);
            tris.push([a, c, d]);
        }
    }
    for i in 0..n_lon {
        tris.push([ring_v(n_lat - 1, i), ring_v(n_lat - 1, i + 1), north_id]);
    }

    let mut cm_hybrid = CellMesh::<3>::new();
    for &r in &radii {
        for u in &u_verts {
            cm_hybrid.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }
    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };

    for k in 0..n_tot {
        let dom = if k < n_phys { DOM_PHYS } else { DOM_IE };
        if k < 2 {
            for &[s0, s1, s2] in &tris {
                let mut pairs = [
                    (s0, vid(k, s0), vid(k + 1, s0)),
                    (s1, vid(k, s1), vid(k + 1, s1)),
                    (s2, vid(k, s2), vid(k + 1, s2)),
                ];
                pairs.sort_by_key(|p| p.0);
                let bot = [pairs[0].1, pairs[1].1, pairs[2].1];
                let top = [pairs[0].2, pairs[1].2, pairs[2].2];
                cm_hybrid.push_cell(CellType::Tet, &[bot[0], bot[1], bot[2], top[2]], dom);
                cm_hybrid.push_cell(CellType::Tet, &[bot[0], bot[1], top[1], top[2]], dom);
                cm_hybrid.push_cell(CellType::Tet, &[bot[0], top[0], top[1], top[2]], dom);
            }
        } else {
            for &[s0, s1, s2] in &tris {
                cm_hybrid.push_cell(
                    CellType::Prism,
                    &[vid(k, s0), vid(k, s1), vid(k, s2), vid(k + 1, s0), vid(k + 1, s1), vid(k + 1, s2)],
                    dom,
                );
            }
        }
    }
    for &[s0, s1, s2] in &tris {
        cm_hybrid.push_facet(&[vid(0, s0), vid(0, s1), vid(0, s2)], 10);
        cm_hybrid.push_facet(&[vid(n_tot, s0), vid(n_tot, s1), vid(n_tot, s2)], BND_OUTER);
    }
    for k in 0..2 {
        for i in 0..n_lon {
            let (s_a, s_b) = (ring_v(j_eq, i), ring_v(j_eq, i + 1));
            let (u_a, u_b) = if s_a < s_b { (s_a, s_b) } else { (s_b, s_a) };
            let (b0, b1) = (vid(k, u_a), vid(k, u_b));
            let (t0, t1) = (vid(k + 1, u_a), vid(k + 1, u_b));
            cm_hybrid.push_facet(&[b0, b1, t1], BND_INNER);
            cm_hybrid.push_facet(&[b0, t0, t1], BND_INNER);
        }
    }

    Mesh::D3(cm_hybrid)
}

#[test]
fn test_circular_loop_inductance_curl_robin_farfield_and_prism_ie() {
    let r_core = 0.4;
    let r_loop_out = 1.0;
    let tight_layers = [1.2, 1.4, 1.6, 1.8];
    let r_tight = 1.8;
    let ref_layers = [1.2, 1.4, 1.6, 1.8, 2.1, 2.5, 3.0, 3.7, 4.6, 6.0];

    let make_loop_bcs = |outer_bc: Bc| -> Vec<Bc> {
        vec![
            Bc { labels: vec![BND_INNER], kind: BcKind::LondonSheet { l_k_sq: 5.0e-5 } },
            Bc { labels: vec![BND_INNER], kind: BcKind::FluxLoop { name: "loop".into(), flux_wb: 1.0 } },
            outer_bc,
        ]
    };

    let mesh_ref = build_equatorial_loop_magneto_mesh(r_core, r_loop_out, &ref_layers, None);
    let spec_ref = magnetostatics_problem(
        "mixed_3d",
        make_loop_bcs(Bc {
            labels: vec![BND_OUTER],
            kind: BcKind::FarField { order: 1, decay_power: 1.0, center: Some([0.0, 0.0, 0.0]) },
        }),
    );
    let res_ref = solve_magnetostatics_mixed(&mesh_ref, &spec_ref).expect("reference loop inductance solve");
    let l_ref = InductanceMatrixData::of(&res_ref).expect("inductance matrix").inductance_matrix_h[0][0];

    let mesh_tight = build_equatorial_loop_magneto_mesh(r_core, r_loop_out, &tight_layers, None);
    assert!(!mesh_tight.is_pure_simplex(), "hybrid Tet+Prism mesh must be mixed");

    let spec_pec =
        magnetostatics_problem("mixed_3d", make_loop_bcs(Bc { labels: vec![BND_OUTER], kind: BcKind::Ground }));
    let res_pec = solve_magnetostatics_mixed(&mesh_tight, &spec_pec).expect("PEC tight box solve");
    let l_pec = InductanceMatrixData::of(&res_pec).unwrap().inductance_matrix_h[0][0];

    let spec_ff = magnetostatics_problem(
        "mixed_3d",
        make_loop_bcs(Bc {
            labels: vec![BND_OUTER],
            kind: BcKind::FarField { order: 1, decay_power: 1.0, center: Some([0.0, 0.0, 0.0]) },
        }),
    );
    let res_ff = solve_magnetostatics_mixed(&mesh_tight, &spec_ff).expect("FarField curl-Robin tight box solve");
    let l_ff = InductanceMatrixData::of(&res_ff).unwrap().inductance_matrix_h[0][0];

    let mesh_ie = build_equatorial_loop_magneto_mesh(r_core, r_loop_out, &tight_layers, Some(2.8));
    let spec_ie = magnetostatics_problem(
        "mixed_3d",
        make_loop_bcs(Bc {
            labels: vec![DOM_IE],
            kind: BcKind::InfiniteElement(rational_spherical_ie_spec(r_tight, 1.0)),
        }),
    );
    let res_ie = solve_magnetostatics_mixed(&mesh_ie, &spec_ie).expect("Native Prism Rational IE magneto 3D solve");
    let l_ie = InductanceMatrixData::of(&res_ie).unwrap().inductance_matrix_h[0][0];

    let bias_pec = (l_pec - l_ref).abs() / l_ref;
    let bias_ff = (l_ff - l_ref).abs() / l_ref;
    let bias_ie = (l_ie - l_ref).abs() / l_ref;
    eprintln!(
        "[Benchmark 4] 3D loop inductance: L_ref = {:.6e} H, PEC bias = {:.4}%, curl-Robin far_field bias = {:.4}% ({:.1}x reduction), Prism Rational IE bias = {:.4}%",
        l_ref,
        bias_pec * 100.0,
        bias_ff * 100.0,
        bias_pec / bias_ff.max(1e-16),
        bias_ie * 100.0
    );

    assert!(
        bias_pec > 10.0 * bias_ff,
        "3D curl-Robin far_field must reduce tight-box PEC wall bias by > 10x (bias_pec = {:.4}%, bias_ff = {:.4}%, L_ref = {:.6e})",
        bias_pec * 100.0,
        bias_ff * 100.0,
        l_ref
    );
    assert!(
        bias_ie < 0.02,
        "3D native Prism Rational IE circular loop inductance error {:.4}% exceeds 2% (L_ie = {:.6e}, L_ref = {:.6e})",
        bias_ie * 100.0,
        l_ie,
        l_ref
    );

    // Also verify 2D `cross_section` magnetostatics with native `Quad` Rational IE & `far_field`
    let mesh_2d = build_annular_quad_ie_mesh(0.5, 1.5, 2.5, 36, 4, 5, false, false);
    let spec_2d_ie = magnetostatics_problem(
        "cross_section",
        vec![
            Bc {
                labels: vec![BND_INNER],
                kind: BcKind::CurrentTerminal { name: "wire".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
            },
            Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(rational_cylindrical_2d_ie_spec(1.5, 1.0)) },
        ],
    );
    let res_2d_ie = solve_magnetostatics_mixed(&mesh_2d, &spec_2d_ie).expect("2D cross_section Quad Rational IE");
    let l_2d = InductanceMatrixData::of(&res_2d_ie).unwrap().inductance_matrix_h[0][0];
    assert!(
        l_2d > 0.0 && l_2d.is_finite(),
        "2D cross_section Quad IE inductance must be positive and finite (got {l_2d}, MU0 = {MU0})"
    );
}

// =================================================================================================
// 5. `heat_transient` Guard Rejecting `infinite_element` (and `far_field`) with Actionable Error
// =================================================================================================

#[test]
fn test_heat_transient_rejects_infinite_element_with_actionable_guard() {
    let mesh = build_spherical_prism_ie_mesh(1.0, 2.0, 3.0, 2, 2, 2, false);

    // (a) Direct Rust API `solve_heat_mixed` / `nyx_dfs::run` guard check
    let mut spec = Problem::new(Study::HeatTransient { t_end: 1.0e-3, n_out: 4, t0: "293.15".into() }, 1);
    spec.bcs = vec![
        Bc { labels: vec![BND_INNER], kind: BcKind::Temperature(350.0) },
        Bc { labels: vec![DOM_IE], kind: BcKind::InfiniteElement(rational_spherical_ie_spec(2.0, 1.0)) },
    ];
    spec.materials.set(DOM_PHYS, "k", 100.0);
    spec.materials.set(DOM_IE, "k", 100.0);

    let err = solve_heat_mixed(&mesh, &spec).expect_err("heat_transient must reject infinite_element");
    assert!(
        err.contains("heat_transient")
            && err.contains("infinite_element")
            && err.contains("rho * Cp * det(J_IE)")
            && err.contains("sqrt(4 * D_max * t_end * ln(1/R_0))"),
        "unexpected heat_transient guard message: {err}"
    );

    // (b) JSON `ProblemSpec` decoding through `Dfs::validate` / `Dfs::solve`
    struct DummyLabels;
    impl LabelResolver for DummyLabels {
        fn domains(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "phys" => Ok(vec![DOM_PHYS]),
                "ie" => Ok(vec![DOM_IE]),
                _ => Err(Error::invalid("unknown domain")),
            }
        }
        fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "inner" => Ok(vec![BND_INNER]),
                "outer" => Ok(vec![BND_OUTER]),
                _ => Err(Error::invalid("unknown boundary")),
            }
        }
    }

    let mut pspec = ProblemSpec::new("dfs", 1, "heat_transient");
    pspec.study.params = json!({
        "t_end": 1.0e-3,
        "n_out": 4,
        "t0": "293.15"
    });
    pspec.materials.insert(
        "phys".into(),
        MaterialAssign::Props([("k".into(), 100.0), ("rho".into(), 2330.0), ("Cp".into(), 700.0)].into()),
    );
    pspec.materials.insert(
        "ie".into(),
        MaterialAssign::Props([("k".into(), 100.0), ("rho".into(), 2330.0), ("Cp".into(), 700.0)].into()),
    );
    pspec.bcs.push(BcSpec { on: "inner".into(), kind: "temperature".into(), params: json!({ "value": 350.0 }) });
    pspec.bcs.push(BcSpec {
        on: "ie".into(),
        kind: "infinite_element".into(),
        params: json!({ "coord": "spherical", "profile": "rational" }),
    });

    let err_json = Dfs
        .solve(&mesh, &pspec, &mut Ctx::new(&DummyLabels))
        .expect_err("Dfs::solve must reject infinite_element on heat_transient");
    let msg = err_json.to_string();
    assert!(
        msg.contains("infinite_element") && msg.contains("sqrt(4 * D_max * t_end * ln(1/R_0))"),
        "unexpected Dfs::solve error: {msg}"
    );
}
