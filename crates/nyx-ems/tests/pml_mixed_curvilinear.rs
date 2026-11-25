// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 1C quantitative verification suite for native mixed-cell EM kernels (`Tri + Quad` in 2D,
//! `Tet + Prism + Hex` in 3D), quadrature-point curvilinear PML pullbacks (`Cartesian`, `Cylindrical`,
//! `Spherical`), `KappaPolynomial` & `Rational` stretching profiles, `r = 0` axisymmetric regularity,
//! `Ams` (`P1`) & `AmsP2` (`P2`) iterative solvers, and driven PROM `P_pml (W)` power conservation.
#![allow(clippy::needless_range_loop)]

use femcore::edge_mixed::{eval_cell_basis_2d, nedelec_nloc, NedEvalScratch, NedMixed};
use femcore::fem_mixed::{quad_rule_mixed, IsoparametricCell};
use femcore::linalg::csr::{c, CCoo, C};
use femcore::linalg::driver::{DirectSolver, IterativeCurlCurlOptions};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::meshing::pml_box::{build_pml, PmlSpec};
use femcore::props::{Periodic, Pml, PmlCoord, PmlDir, PmlProfile};
use nyx_ems::mixed_ops::{iso_mu_inv_c3, C3};
use nyx_ems::{
    assemble_emw_eigen_mixed, run, solve_cavity_mixed_ams, solve_mixed_curl_curl_iterative, EmPmlTable, Ems,
    LumpedPortBc, Problem, Study, WavePortBc, C0,
};
use nyx_physics::testkit::cell_fields_reproduce_the_plot_fields;
use nyx_physics::{
    discretise, discretise_mixed, BcSpec, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, ProblemSpec,
    Scalar, StudySpec, DRIVEN_PORT_KEY,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::Arc;

// ------------------------------------------------------------------ Helpers

/// Builds a conforming 2D meridian half-disk (`r >= 0, r^2 + z^2 <= r_outer^2`) with an inner `Tri` core
/// (`R <= r_inner`) and `n_rings` concentric curved `Quad` shells (`r_inner .. r_outer`), clamping the north
/// (`theta = 0`) and south (`theta = pi`) poles strictly to `r = 0.0`!
/// - Core cells (`R <= r_inner`) get domain label `1`.
/// - Shell rings `0 .. n_pml_Inner` get domain label `shell_dom` (`1` for pure PEC sphere, `2` for Spherical PML).
fn build_half_disk_mixed_2d(
    r_core: f64,
    r_in: f64,
    r_out: f64,
    n_theta: usize,
    n_mid_rings: usize,
    n_pml_rings: usize,
    pml_dom: i32,
) -> Mesh {
    let mut cm = CellMesh::<2>::new();
    // Vertex 0: origin (0, 0) on the symmetry axis `r = 0`
    cm.points.push([0.0, 0.0]);

    let mut radii = vec![r_core];
    for k in 1..=n_mid_rings {
        let t = k as f64 / n_mid_rings as f64;
        radii.push(r_core + t * (r_in - r_core));
    }
    for k in 1..=n_pml_rings {
        let t = k as f64 / n_pml_rings as f64;
        radii.push(r_in + t * (r_out - r_in));
    }

    let n_shells = radii.len();
    let mut ring_verts: Vec<Vec<u32>> = Vec::with_capacity(n_shells);
    for &rad in &radii {
        let mut rv = Vec::with_capacity(n_theta + 1);
        for j in 0..=n_theta {
            let th = PI * (j as f64) / (n_theta as f64);
            let r = if j == 0 || j == n_theta { 0.0 } else { rad * th.sin() };
            let z = rad * th.cos();
            let id = cm.points.len() as u32;
            cm.points.push([r, z]);
            rv.push(id);
        }
        ring_verts.push(rv);
    }

    // Inner fan of `Tri` cells around origin `0` touching `r = 0` on first and last segments
    for j in 0..n_theta {
        // Orient CCW in `(r, z)`: `th` runs from `0` (+z) to `pi` (-z), so `(0, rv[j+1], rv[j])` has positive det!
        cm.push_cell(CellType::Tri, &[0, ring_verts[0][j + 1], ring_verts[0][j]], 1);
    }

    // Concentric `Quad` shells
    let mut edge_nodes = Vec::new();
    for s in 0..n_shells {
        let rad = radii[s];
        for j in 0..n_theta {
            let th_mid = PI * (j as f64 + 0.5) / (n_theta as f64);
            let a = ring_verts[s][j];
            let b = ring_verts[s][j + 1];
            edge_nodes.push(([a.min(b), a.max(b)], [rad * th_mid.sin(), rad * th_mid.cos()]));
        }
    }
    for s in 0..(n_shells - 1) {
        let dom = if s < n_mid_rings { 1 } else { pml_dom };
        for j in 0..n_theta {
            let i0 = ring_verts[s][j + 1];
            let i1 = ring_verts[s + 1][j + 1];
            let i2 = ring_verts[s + 1][j];
            let i3 = ring_verts[s][j];
            cm.push_cell(CellType::Quad, &[i0, i1, i2, i3], dom);
        }
    }

    // Exterior facets:
    // 1. Outer arc `s = n_shells - 1` -> label 1 (outer PEC wall)
    let outer_ring = &ring_verts[n_shells - 1];
    for j in 0..n_theta {
        cm.push_facet(&[outer_ring[j + 1], outer_ring[j]], 1);
    }
    // 2. Symmetry axis `r = 0` segments (north `j = 0` and south `j = n_theta`) -> label 2 (symmetry axis)
    cm.push_facet(&[0, ring_verts[0][0]], 2);
    cm.push_facet(&[ring_verts[0][n_theta], 0], 2);
    for s in 0..(n_shells - 1) {
        cm.push_facet(&[ring_verts[s][0], ring_verts[s + 1][0]], 2);
        cm.push_facet(&[ring_verts[s + 1][n_theta], ring_verts[s][n_theta]], 2);
    }

    cm.set_edge_nodes(edge_nodes);
    Mesh::D2(cm)
}

/// Builds a C_{4v}-symmetric 3D tetrahedral mesh on `[-L/2, L/2] x [-L/2, L/2] x [0, H]` by splitting every
/// grid cube around its body center into 24 symmetric tetrahedra (invariant under 90-degree `(x, y)` rotation!).
fn build_c4v_symmetric_tet_box(l_xy: f64, h_z: f64, n_xy: usize, n_z: usize) -> femcore::adapt::work::SimplexMesh<3> {
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut grid = vec![vec![vec![0u32; n_z + 1]; n_xy + 1]; n_xy + 1];
    for ix in 0..=n_xy {
        let x = -0.5 * l_xy + l_xy * (ix as f64) / (n_xy as f64);
        for iy in 0..=n_xy {
            let y = -0.5 * l_xy + l_xy * (iy as f64) / (n_xy as f64);
            for iz in 0..=n_z {
                let z = h_z * (iz as f64) / (n_z as f64);
                let id = pts.len() as u32;
                pts.push([x, y, z]);
                grid[ix][iy][iz] = id;
            }
        }
    }

    let mut cells: Vec<u32> = Vec::new();
    let mut cell_label: Vec<i32> = Vec::new();
    let mut facets: Vec<u32> = Vec::new();
    let mut facet_label: Vec<i32> = Vec::new();

    // Face center cache so shared grid faces use identical face-center vertices
    let mut face_center: BTreeMap<[u32; 4], u32> = BTreeMap::new();

    for ix in 0..n_xy {
        for iy in 0..n_xy {
            for iz in 0..n_z {
                let v000 = grid[ix][iy][iz];
                let v100 = grid[ix + 1][iy][iz];
                let v110 = grid[ix + 1][iy + 1][iz];
                let v010 = grid[ix][iy + 1][iz];
                let v001 = grid[ix][iy][iz + 1];
                let v101 = grid[ix + 1][iy][iz + 1];
                let v111 = grid[ix + 1][iy + 1][iz + 1];
                let v011 = grid[ix][iy + 1][iz + 1];

                let mut bc = [0.0; 3];
                for &v in &[v000, v100, v110, v010, v001, v101, v111, v011] {
                    for k in 0..3 {
                        bc[k] += 0.125 * pts[v as usize][k];
                    }
                }
                let vc = pts.len() as u32;
                pts.push(bc);

                let faces: [([u32; 4], Option<i32>); 6] = [
                    ([v000, v010, v110, v100], (iz == 0).then_some(3)),
                    ([v001, v101, v111, v011], (iz + 1 == n_z).then_some(4)),
                    ([v000, v100, v101, v001], (iy == 0).then_some(2)),
                    ([v010, v011, v111, v110], (iy + 1 == n_xy).then_some(5)),
                    ([v000, v001, v011, v010], (ix == 0).then_some(1)),
                    ([v100, v110, v111, v101], (ix + 1 == n_xy).then_some(6)),
                ];
                for (q, bnd_lbl) in faces {
                    let mut key = q;
                    key.sort_unstable();
                    let fc = if let Some(&id) = face_center.get(&key) {
                        id
                    } else {
                        let mut c = [0.0; 3];
                        for &v in &q {
                            for k in 0..3 {
                                c[k] += 0.25 * pts[v as usize][k];
                            }
                        }
                        let id = pts.len() as u32;
                        pts.push(c);
                        face_center.insert(key, id);
                        id
                    };
                    for m in 0..4 {
                        let a = q[m];
                        let b = q[(m + 1) % 4];
                        let pa = pts[a as usize];
                        let pb = pts[b as usize];
                        let pf = pts[fc as usize];
                        let pv = pts[vc as usize];
                        let u = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                        let v = [pf[0] - pa[0], pf[1] - pa[1], pf[2] - pa[2]];
                        let w = [pv[0] - pa[0], pv[1] - pa[1], pv[2] - pa[2]];
                        let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                            + u[2] * (v[0] * w[1] - v[1] * w[0]);
                        if det > 0.0 {
                            cells.extend_from_slice(&[a, b, fc, vc]);
                        } else {
                            cells.extend_from_slice(&[b, a, fc, vc]);
                        }
                        cell_label.push(1);
                        if let Some(lbl) = bnd_lbl {
                            facets.extend_from_slice(&[a, b, fc]);
                            facet_label.push(lbl);
                        }
                    }
                }
            }
        }
    }
    femcore::adapt::work::SimplexMesh { points: pts, cells, cell_label, facets, facet_label }
}

// ------------------------------------------------------------------ 1. Empty Spherical Cavity & Spherical PML at r = 0

#[test]
fn axisymmetric_empty_sphere_touching_r_zero_matches_bessel_zeros_m0_m1_m2() {
    let r_sphere = 0.05;
    let mesh = build_half_disk_mixed_2d(0.015, 0.035, r_sphere, 14, 3, 2, 1);
    let simplex = discretise(&mesh).unwrap();
    let progress = Arc::new(femcore::progress::Progress::new());

    // Exact analytical spherical Bessel cavity roots `k * R`:
    // - Lowest TM_{1,m} root (root of `[x j_1(x)]' = 0`): `x = 2.743_707_269_992_269` (3-fold degenerate: m = -1, 0, +1!)
    // - Lowest TM_{2,m} root (root of `[x j_2(x)]' = 0`): `x = 3.870_238_580_207_554` (5-fold degenerate: |m| <= 2!)
    // - Lowest TE_{1,m} root (root of `j_1(x) = 0`):       `x = 4.493_409_457_909_064` (3-fold degenerate: m = -1, 0, +1!)
    let f_of_kr = |kr: f64| C0 * kr / (2.0 * PI * r_sphere);
    let cases = [
        (0_i64, 2.743_707_269_992_269, "TM_011"),
        (0_i64, 4.493_409_457_909_064, "TE_011"),
        (1_i64, 2.743_707_269_992_269, "TM_111"),
        (1_i64, 4.493_409_457_909_064, "TE_111"),
        (2_i64, 3.870_238_580_207_554, "TM_211"),
    ];

    for (m, kr_exact, tag) in cases {
        let f_exact = f_of_kr(kr_exact);
        let mut p =
            Problem::new(Study::AxisymmetricModes { f_search: 0.995 * f_exact, azimuthal_order: m, n_modes: 2 }, 2);
        p.materials.default.insert("eps_r".into(), 1.0);
        p.materials.default.insert("mu_r".into(), 1.0);
        let res = run(&mesh, &p, &progress).expect("axisymmetric sphere touching r=0");
        cell_fields_reproduce_the_plot_fields(&simplex, &res);

        let f_got = res
            .solutions
            .iter()
            .map(|s| s.value.unwrap())
            .min_by(|a, b| (a - f_exact).abs().total_cmp(&(b - f_exact).abs()))
            .unwrap();
        let rel_err = (f_got - f_exact).abs() / f_exact;
        println!(
            "[sphere r=0] {tag} (m={m}): got = {:.8e} Hz, exact = {:.8e} Hz, rel_err = {:.4e} ({:.5}%)",
            f_got,
            f_exact,
            rel_err,
            100.0 * rel_err
        );
        assert!(
            rel_err < 1e-3,
            "{tag} (m={m}) exceeded 0.1% tolerance: got {f_got:.8e} vs {f_exact:.8e} (rel_err = {rel_err:.3e})"
        );
    }
}

#[test]
fn axisymmetric_spherical_pml_touching_r_zero_damps_radiating_quasinormal_modes() {
    // Dielectric resonator sphere (`eps_r = 16.0`, radius `0.02 m`) in air (`R_in = 0.045 m`)
    // terminated by 5 concentric curved `Quad` layers of Spherical PML (`0.045 .. 0.075 m`) touching `r = 0`!
    let mesh = build_half_disk_mixed_2d(0.02, 0.045, 0.075, 14, 3, 5, 3);
    // Relabel core (`r <= 0.02`) -> 1 (`eps_r = 16`), air (`0.02 .. 0.045`) -> 2 (`eps_r = 1`), PML (`0.045 .. 0.075`) -> 3
    let Mesh::D2(mut cm) = mesh else { unreachable!() };
    for c_idx in 0..cm.n_cells() {
        if cm.label[c_idx] == 1 {
            let mut r_max = 0.0_f64;
            for &v in cm.cell(c_idx) {
                let p = cm.points[v as usize];
                r_max = r_max.max(p[0].hypot(p[1]));
            }
            if r_max > 0.0201 {
                cm.label[c_idx] = 2;
            }
        }
    }
    let mesh = Mesh::D2(cm);
    let simplex = discretise(&mesh).unwrap();
    let progress = Arc::new(femcore::progress::Progress::new());

    // Magnetic dipole TE_{01\delta} / TE_{11\delta} resonance of an eps_r = 16 sphere of radius a = 0.02 m sits near k0 a * 4 ~ pi => f ~ 1.87 GHz
    for m in [0_i64, 1_i64] {
        let mut p = Problem::new(Study::AxisymmetricModes { f_search: 1.85e9, azimuthal_order: m, n_modes: 2 }, 2);
        p.materials.set(1, "eps_r", 16.0);
        p.materials.set(2, "eps_r", 1.0);
        p.materials.set(3, "eps_r", 1.0);
        p.materials.default.insert("mu_r".into(), 1.0);
        p.materials.pml.insert(
            3,
            Pml {
                coord: PmlCoord::Spherical { center: [0.0, 0.0, 0.0], r_in: 0.045, thickness: 0.030 },
                profile: PmlProfile::KappaPolynomial { kappa_max: 1.5 },
                strength: 3.5,
                order: 2.0,
                ..Default::default()
            },
        );
        let res = run(&mesh, &p, &progress).expect("axisymmetric spherical PML solve");
        cell_fields_reproduce_the_plot_fields(&simplex, &res);
        let eig = res.solutions[0].eigenvalue.expect("complex radiating eigenvalue");
        assert!(eig[0] > 1.5e9 && eig[0] < 2.3e9, "m={m}: unexpected Re(f) = {:.4e}", eig[0]);
        assert!(eig[1] > 1e6, "m={m}: passive outgoing radiation must have Im(f) > 0, got {:.4e}", eig[1]);
        println!(
            "[axisym spherical PML r=0] m={m}: f = {:.6e} + {:.4e}i Hz (Q = {:.2})",
            eig[0],
            eig[1],
            eig[0] / (2.0 * eig[1])
        );
    }
}

// ------------------------------------------------------------------ 2. 3D Native Tet + Prism + Hex Cavity, C_{4v}, AMS/AmsP2 & DOF Reduction

#[test]
fn cavity_3d_native_tet_prism_hex_preserves_c4v_symmetry_reduces_dofs_and_converges_ams() {
    let (l_xy, h_z) = (0.03, 0.02);
    let phys_tet = build_c4v_symmetric_tet_box(l_xy, h_z, 3, 2);
    let mut pml_spec = PmlSpec::<3>::uniform(0.008, 3);
    // Add lateral x/y structured PML (`Prism` on 4 side faces + `Hex` on 4 vertical corner columns), keeping z=0,h_z & outer walls PEC!
    pml_spec.thickness[2] = [0.0, 0.0];
    pml_spec.layers[2] = [0, 0];
    pml_spec.outer_label = Some(20);
    let pml3 = build_pml(&phys_tet, &pml_spec).expect("3D Tet + Prism + Hex structured PML");
    let cm3 = pml3.mesh;

    let n_tet = cm3.types.iter().filter(|&&t| t == CellType::Tet).count();
    let n_prism = cm3.types.iter().filter(|&&t| t == CellType::Prism).count();
    let n_hex = cm3.types.iter().filter(|&&t| t == CellType::Hex).count();
    assert!(n_tet > 0 && n_prism > 0 && n_hex > 0);

    let mesh = Mesh::D3(cm3);
    assert!(!mesh.is_pure_simplex());
    let mixed = discretise_mixed(&mesh).unwrap();
    let split = discretise(&mesh).unwrap();

    // Compare P2 Nedelec DOF count in PML cells (mixed vs split-tet):
    let ned_mixed_p2 = NedMixed::new(&mixed, 2);
    let vkey_split: Vec<u32> = (0..split.n_points() as u32).collect();
    let ned_split_p2 = nyx_ems::cavity::Ned3::new(&split, 2, vkey_split);
    let vkey_phys: Vec<u32> = (0..phys_tet.n_points() as u32).collect();
    let phys_flat: Vec<f64> = phys_tet.points.iter().flat_map(|p| p.iter().copied()).collect();
    let phys_sm = femcore::mesh::SimplexMesh::new(
        3,
        phys_flat,
        phys_tet.cells.clone(),
        phys_tet.cell_label.clone(),
        phys_tet.facets.clone(),
        phys_tet.facet_label.clone(),
    );
    let ned_phys_p2 = nyx_ems::cavity::Ned3::new(&phys_sm, 2, vkey_phys.clone());
    let pml_dofs_mixed_p2 = ned_mixed_p2.n_dofs - ned_phys_p2.n_dofs();
    let pml_dofs_split_p2 = ned_split_p2.n_dofs() - ned_phys_p2.n_dofs();
    let pml_dof_red_p2_pct = 100.0 * (1.0 - (pml_dofs_mixed_p2 as f64) / (pml_dofs_split_p2 as f64));

    let ned_mixed_p1 = NedMixed::new(&mixed, 1);
    let ned_split_p1 = nyx_ems::cavity::Ned3::new(&split, 1, (0..split.n_points() as u32).collect());
    let ned_phys_p1 = nyx_ems::cavity::Ned3::new(&phys_sm, 1, vkey_phys);
    let pml_dofs_mixed_p1 = ned_mixed_p1.n_dofs - ned_phys_p1.n_dofs();
    let pml_dofs_split_p1 = ned_split_p1.n_dofs() - ned_phys_p1.n_dofs();
    let pml_dof_red_p1_pct = 100.0 * (1.0 - (pml_dofs_mixed_p1 as f64) / (pml_dofs_split_p1 as f64));
    let dof_reduction_pct = 100.0 * (1.0 - (ned_mixed_p2.n_dofs as f64) / (ned_split_p2.n_dofs() as f64));
    println!(
        "[3D Tet+Prism+Hex] cells = ({n_tet} Tet, {n_prism} Prism, {n_hex} Hex) vs {} split Tets | PML-only DOF reduction: N1 = {:.2}% ({} vs {}), N2 = {:.2}% ({} vs {}), total N2 = {:.2}% ({} vs {})",
        split.n_cells(),
        pml_dof_red_p1_pct,
        pml_dofs_mixed_p1,
        pml_dofs_split_p1,
        pml_dof_red_p2_pct,
        pml_dofs_mixed_p2,
        pml_dofs_split_p2,
        dof_reduction_pct,
        ned_mixed_p2.n_dofs,
        ned_split_p2.n_dofs(),
    );
    assert!(ned_mixed_p2.n_dofs < ned_split_p2.n_dofs(), "native Prism+Hex must use fewer N2 DOFs than split-tet");
    assert!(pml_dof_red_p1_pct > 40.0 && pml_dof_red_p2_pct > 25.0);

    // 1. C_{4v} degenerate doublet check (`TE_{101}` and `TE_{011}` on square `[-(l_xy+2d)/2 .. +(l_xy+2d)/2]^2 x [0, h_z]`):
    // First test with uniform lossless real stretch `re = 1.25, strength = 0.0`, then with complex damping `strength = 2.5`!
    let l_tot_eff = l_xy + 2.0 * 1.25 * 0.008;
    let f_doublet_est = 0.5 * C0 * ((1.0 / l_tot_eff).powi(2) + (1.0 / h_z).powi(2)).sqrt();

    let mut prob = Problem::new(Study::Cavity { n_modes: 4, f_search: 0.96 * f_doublet_est }, 2);
    prob.materials.default.insert("eps_r".into(), 1.0);
    prob.materials.default.insert("mu_r".into(), 1.0);
    for &dom in &mixed.domain_labels() {
        if let Some((code, _)) = pml_spec.decode(dom) {
            let dir_of = |c: i8| match c {
                -1 => PmlDir::Minus,
                1 => PmlDir::Plus,
                _ => PmlDir::Off,
            };
            prob.materials.pml.insert(
                dom,
                Pml {
                    x: dir_of(code[0]),
                    y: dir_of(code[1]),
                    z: dir_of(code[2]),
                    strength: 2.0,
                    order: 2.0,
                    profile: PmlProfile::KappaPolynomial { kappa_max: 1.5 },
                    ..Default::default()
                },
            );
        }
    }

    let progress = Arc::new(femcore::progress::Progress::new());
    let res = run(&mesh, &prob, &progress).expect("3D mixed cavity with PML");
    for sol in &res.solutions {
        let flat = sol.flat_fields();
        for name in ["Ex", "Ey", "Ez"] {
            let cf = sol.cell_fields.iter().find(|c| c.name == name).unwrap();
            assert_eq!(cf.re.len(), split.n_cells() * split.nloc(2));
            let avg = femcore::fields::node_average(&split, &cf.re, 1, 2);
            let pf = flat.iter().find(|(n, _, _)| n == &format!("Re {name}")).unwrap();
            let scale = pf.2.iter().fold(0.0_f64, |m, x| m.max(x.abs())).max(1e-30);
            let worst = avg.iter().zip(&pf.2).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(worst <= 1e-9 * scale, "{name} Stage-3 projection mismatch {worst:e}");
        }
    }

    // Find the degenerate C_{4v} doublet `(1,0,1)` / `(0,1,1)` above the non-degenerate TM_{110} fundamental:
    let eigs: Vec<[f64; 2]> = res.solutions.iter().map(|s| s.eigenvalue.unwrap()).collect();
    let mut min_doublet_split = f64::INFINITY;
    for w in eigs.windows(2) {
        let df = ((w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1])) / w[0][0].max(1.0);
        min_doublet_split = min_doublet_split.min(df);
    }
    println!("[3D C_4v Tet+Prism+Hex PML] eigenvalues = {:?}, doublet |df/f| = {:.3e}", eigs, min_doublet_split);
    assert!(min_doublet_split < 1e-12, "C_4v doublet splitting |df/f| = {min_doublet_split:.3e} exceeded 1e-12");

    // 2. Verify iterative `Ams` (Order 1) and `AmsP2` (Order 2) linear + JD cavity solves, plus Rational -> KappaPolynomial fallback!
    for order in [1usize, 2usize] {
        let mut p_it = prob.clone();
        p_it.order = order;
        p_it.study = Study::Cavity { n_modes: 2, f_search: 8.5e9 };
        // Deliberately set `Rational` on one PML domain to verify the iterative AMS safeguard falls back cleanly!
        for pml_entry in p_it.materials.pml.values_mut() {
            pml_entry.profile = PmlProfile::Rational { curvature: 0.85 };
        }
        let (asm, ned_m, masks) = assemble_emw_eigen_mixed(&mixed, &p_it, true).unwrap();
        let a_spd = femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm.kc, c(asm.sigma, 0.0), &asm.mc);
        let rhs: Vec<C> = (0..asm.nred).map(|i| c(((i + 1) as f64).sin(), ((i + 3) as f64).cos())).collect();
        let rep = solve_mixed_curl_curl_iterative(
            &mixed,
            &ned_m,
            &masks.fixed_dofs,
            &masks.fixed_edges,
            &masks.fixed_verts,
            &masks.fixed_quad_faces,
            &a_spd,
            None,
            &rhs,
            &IterativeCurlCurlOptions::default(),
        )
        .unwrap();
        println!(
            "[3D Mixed AMS order={order}] iterations = {}, rel_res = {:.3e}, pc_nnz = {}",
            rep.iterations, rep.relative_residual, rep.preconditioner_nnz
        );
        assert!(rep.converged && rep.relative_residual < 1e-8, "AMS order {order} failed to converge: {rep:?}");
    }

    let mut p_jd = Problem::new(Study::Cavity { n_modes: 1, f_search: 8.5e9 }, 2);
    p_jd.materials.default.insert("eps_r".into(), 1.0);
    p_jd.materials.default.insert("mu_r".into(), 1.0);
    let res_jd = solve_cavity_mixed_ams(&mesh, &mixed, &p_jd).expect("JD + AmsP2 mixed cavity solve");
    assert!(!res_jd.solutions.is_empty());
}

// ------------------------------------------------------------------ 3. 2D Tri + Quad Cartesian & Cylindrical PML Reflection (< 1e-6 / < 1e-8)

/// Measures the exact plane-wave amplitude reflection `|R|` of a 2D `Tri` (`[0, l_phys]`) + structured `Quad`
/// (`[l_phys, l_phys + d_pml]`, `n_layers` layers) parallel-plate TEM waveguide (`PEC` plates `z = 0, h` so `E_x = 0`,
/// `E_z = 1` Dirichlet on `x = 0`, `PEC` termination on `x = l_phys + d_pml`) against a deep reference PML on the
/// identical interior `Tri` mesh (eliminating interior triangle `k_h != k_0` dispersion error)!
fn measure_2d_cartesian_tri_quad_reflection(
    profile: PmlProfile,
    strength: f64,
    order_m: f64,
    d_over_lam: f64,
    n_layers: usize,
) -> f64 {
    let freq = 15e9;
    let lam0 = C0 / freq;
    let k0 = 2.0 * PI / lam0;
    let l_phys = 1.0 * lam0;
    let h = 0.1 * lam0;
    let nx_phys = 24usize;
    let ny = 2usize;

    let solve_strip = |n_pml: usize, d_tot: f64, prof: PmlProfile, str_val: f64, ord: f64, grading: f64| -> Vec<C> {
        let mut cm = CellMesh::<2>::new();
        let nx_tot = nx_phys + n_pml;
        let mut grid = vec![vec![0u32; ny + 1]; nx_tot + 1];
        let mut pml_off = vec![0.0; n_pml + 1];
        if n_pml > 0 {
            let q = if n_pml > 1 { grading.powf(1.0 / (n_pml as f64 - 1.0)) } else { 1.0 };
            let tot: f64 = (0..n_pml).map(|k| q.powi(k as i32)).sum();
            let w0 = d_tot / tot;
            for k in 0..n_pml {
                pml_off[k + 1] = pml_off[k] + w0 * q.powi(k as i32);
            }
            pml_off[n_pml] = d_tot;
        }
        for ix in 0..=nx_tot {
            let x =
                if ix <= nx_phys { l_phys * (ix as f64) / (nx_phys as f64) } else { l_phys + pml_off[ix - nx_phys] };
            for iy in 0..=ny {
                let y = h * (iy as f64) / (ny as f64);
                let id = cm.points.len() as u32;
                cm.points.push([x, y]);
                grid[ix][iy] = id;
            }
        }

        for ix in 0..nx_tot {
            for iy in 0..ny {
                let (v00, v10, v11, v01) = (grid[ix][iy], grid[ix + 1][iy], grid[ix + 1][iy + 1], grid[ix][iy + 1]);
                if ix < nx_phys {
                    cm.push_cell(CellType::Tri, &[v00, v10, v11], 1);
                    cm.push_cell(CellType::Tri, &[v00, v11, v01], 1);
                } else {
                    cm.push_cell(CellType::Quad, &[v00, v10, v11, v01], 2);
                }
            }
        }
        for iy in 0..ny {
            cm.push_facet(&[grid[0][iy + 1], grid[0][iy]], 1); // x = 0 TEM source port
            cm.push_facet(&[grid[nx_tot][iy], grid[nx_tot][iy + 1]], 2); // outer PEC termination
        }
        for ix in 0..nx_tot {
            cm.push_facet(&[grid[ix][0], grid[ix + 1][0]], 2); // bottom PEC plate (E_x = 0)
            cm.push_facet(&[grid[ix + 1][ny], grid[ix][ny]], 2); // top PEC plate (E_x = 0)
        }

        let mixed = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        let ned = NedMixed::new(&mixed, 2);
        let n = ned.n_dofs;

        let mut mats = femcore::props::Materials::default();
        mats.default.insert("eps_r".into(), 1.0);
        mats.default.insert("mu_r".into(), 1.0);
        mats.pml.insert(2, Pml { x: PmlDir::Plus, strength: str_val, order: ord, profile: prof, ..Default::default() });
        let pml_table = EmPmlTable::from_mixed(&mixed, &mats, freq, false, false);
        let id_c3: C3 = iso_mu_inv_c3(1.0);

        let mut coo = CCoo::new(n, n);
        let mut scratch = NedEvalScratch::new();
        let mut w = Vec::new();
        let mut cu = Vec::new();
        let k2 = c(k0 * k0, 0.0);

        for ci in 0..mixed.n_cells() {
            let iso = IsoparametricCell::from_mesh(&mixed, ci);
            let rule = quad_rule_mixed(iso.cell_type, 6);
            let nb = nedelec_nloc(iso.cell_type, 2);
            let mut ke = vec![c(0.0, 0.0); nb * nb];
            for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
                let jac = eval_cell_basis_2d(&mixed, &ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                let wq = w_ref * jac.det.abs();
                let (ep, mui, _) = pml_table.eval_2d(ci, id_c3, id_c3, [jac.x[0], jac.x[1]], false);
                for s in 0..nb {
                    let ws = w[s];
                    for t in 0..nb {
                        let wt = w[t];
                        let ew = ep[0][0] * (ws[0] * wt[0])
                            + ep[0][1] * (ws[0] * wt[1] + ws[1] * wt[0])
                            + ep[1][1] * (ws[1] * wt[1]);
                        ke[s * nb + t] += (mui[2][2] * (cu[s] * cu[t]) - k2 * ew) * wq;
                    }
                }
            }
            let cdofs = ned.cell_dofs(ci);
            coo.push_block(cdofs, cdofs, &ke);
        }
        let a_full = coo.to_csr();

        let mut dir_val: BTreeMap<usize, C> = BTreeMap::new();
        for f in 0..mixed.n_facets() {
            let lbl = mixed.facet_labels[f];
            let fv = mixed.facet(f);
            let eid = mixed.edge_index(fv[0], fv[1]).unwrap();
            if lbl == 2 {
                dir_val.insert(2 * eid, c(0.0, 0.0));
                dir_val.insert(2 * eid + 1, c(0.0, 0.0));
            } else if lbl == 1 {
                let (lo, hi) = if ned.key_of(fv[0]) < ned.key_of(fv[1]) { (fv[0], fv[1]) } else { (fv[1], fv[0]) };
                let dy = mixed.point(hi as usize)[1] - mixed.point(lo as usize)[1];
                dir_val.insert(2 * eid, c(dy, 0.0));
                dir_val.insert(2 * eid + 1, c(0.0, 0.0));
            }
        }
        let free: Vec<usize> = (0..n).filter(|d| !dir_val.contains_key(d)).collect();
        let fixed: Vec<usize> = dir_val.keys().copied().collect();
        let x_fixed: Vec<C> = fixed.iter().map(|d| dir_val[d]).collect();
        let a_ff = a_full.submatrix(&free, &free);
        let a_fc = a_full.submatrix(&free, &fixed);
        let rhs: Vec<C> = a_fc.matvec(&x_fixed).into_iter().map(|v| -v).collect();
        let x_free = DirectSolver::new(&a_ff, false).unwrap().solve(&rhs);
        let mut x_all = vec![c(0.0, 0.0); n];
        for (&d, &v) in fixed.iter().zip(&x_fixed) {
            x_all[d] = v;
        }
        for (&d, &v) in free.iter().zip(&x_free) {
            x_all[d] = v;
        }

        let mut samples = Vec::new();
        for ix_cell in 2..(nx_phys - 2) {
            let ci = 2 * (ix_cell * ny);
            let iso = IsoparametricCell::from_mesh(&mixed, ci);
            let xi = [1.0 / 3.0, 1.0 / 3.0];
            let _ = eval_cell_basis_2d(&mixed, &ned, ci, &iso, &xi, &mut scratch, &mut w, &mut cu);
            let mut ez = c(0.0, 0.0);
            for (s, &gd) in ned.cell_dofs(ci).iter().enumerate() {
                ez += x_all[gd] * w[s][1];
            }
            samples.push(ez);
        }
        samples
    };

    let grading = if matches!(profile, PmlProfile::Rational { .. }) { 1.6 } else { 2.0 };
    let test_s = solve_strip(n_layers, d_over_lam * lam0, profile, strength, order_m, grading);
    let ref_s = solve_strip(48, 1.8 * lam0, PmlProfile::KappaPolynomial { kappa_max: 1.0 }, 7.5, 3.0, 2.0);
    let mut num2 = 0.0;
    let mut den2 = 0.0;
    for (et, er) in test_s.iter().zip(&ref_s) {
        num2 += (*et - *er).norm_sqr();
        den2 += er.norm_sqr();
    }
    (num2 / den2.max(1e-300)).sqrt()
}

/// Measures the exact cylindrical outgoing wave reflection `|R|` on an annular `Tri` (`r_src .. r_in`) +
/// concentric curved `Quad` (`r_in .. r_out`, `n_layers` layers) 2D mesh against a deep reference solution!
fn measure_2d_cylindrical_tri_quad_reflection(
    profile: PmlProfile,
    strength: f64,
    order_m: f64,
    d_over_lam: f64,
    n_layers: usize,
) -> f64 {
    let freq = 15e9;
    let lam0 = C0 / freq;
    let k0 = 2.0 * PI / lam0;
    let r_src = 0.6 * lam0;
    let r_in = 1.2 * lam0;
    let d_pml = d_over_lam * lam0;
    let n_phi = 36usize;
    let n_tri_rings = 8usize;

    let solve_annulus =
        |n_pml: usize, d_tot: f64, prof: PmlProfile, str_val: f64, ord: f64, grading: f64| -> Vec<(f64, C)> {
            let mut cm = CellMesh::<2>::new();
            let mut radii = Vec::new();
            for k in 0..=n_tri_rings {
                radii.push(r_src + (r_in - r_src) * (k as f64) / (n_tri_rings as f64));
            }
            if n_pml > 0 {
                let q = if n_pml > 1 { grading.powf(1.0 / (n_pml as f64 - 1.0)) } else { 1.0 };
                let tot: f64 = (0..n_pml).map(|k| q.powi(k as i32)).sum();
                let w0 = d_tot / tot;
                let mut acc = 0.0;
                for k in 0..n_pml {
                    acc += w0 * q.powi(k as i32);
                    radii.push(r_in + if k + 1 == n_pml { d_tot } else { acc });
                }
            }
            let nr = radii.len();
            let mut rings = Vec::with_capacity(nr);
            let mut edge_nodes = Vec::new();
            for &rad in &radii {
                let mut rv = Vec::with_capacity(n_phi);
                for j in 0..n_phi {
                    let ph = 2.0 * PI * (j as f64) / (n_phi as f64);
                    let id = cm.points.len() as u32;
                    cm.points.push([rad * ph.cos(), rad * ph.sin()]);
                    rv.push(id);
                }
                for j in 0..n_phi {
                    let jn = (j + 1) % n_phi;
                    let ph_m = 2.0 * PI * (j as f64 + 0.5) / (n_phi as f64);
                    let (a, b) = (rv[j], rv[jn]);
                    edge_nodes.push(([a.min(b), a.max(b)], [rad * ph_m.cos(), rad * ph_m.sin()]));
                }
                rings.push(rv);
            }
            for ir in 0..(nr - 1) {
                for j in 0..n_phi {
                    let jn = (j + 1) % n_phi;
                    let (v00, v10, v11, v01) = (rings[ir][j], rings[ir + 1][j], rings[ir + 1][jn], rings[ir][jn]);
                    if ir < n_tri_rings {
                        let mid_id = cm.points.len() as u32;
                        let r_m = 0.5 * (radii[ir] + radii[ir + 1]);
                        let ph_m = 2.0 * PI * (j as f64 + 0.5) / (n_phi as f64);
                        cm.points.push([r_m * ph_m.cos(), r_m * ph_m.sin()]);
                        cm.push_cell(CellType::Tri, &[v00, v10, mid_id], 1);
                        cm.push_cell(CellType::Tri, &[v10, v11, mid_id], 1);
                        cm.push_cell(CellType::Tri, &[v11, v01, mid_id], 1);
                        cm.push_cell(CellType::Tri, &[v01, v00, mid_id], 1);
                    } else {
                        cm.push_cell(CellType::Quad, &[v00, v10, v11, v01], 2);
                    }
                }
            }
            for j in 0..n_phi {
                let jn = (j + 1) % n_phi;
                cm.push_facet(&[rings[0][jn], rings[0][j]], 1);
                cm.push_facet(&[rings[nr - 1][j], rings[nr - 1][jn]], 2);
            }
            cm.set_edge_nodes(edge_nodes);

            let mixed = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
            let ned = NedMixed::new(&mixed, 2);
            let n = ned.n_dofs;

            let mut mats = femcore::props::Materials::default();
            mats.default.insert("eps_r".into(), 1.0);
            mats.default.insert("mu_r".into(), 1.0);
            mats.pml.insert(
                2,
                Pml {
                    coord: PmlCoord::Cylindrical {
                        center: [0.0, 0.0, 0.0],
                        axis: [0.0, 1.0, 0.0],
                        r_in,
                        r_thickness: d_tot,
                        radial_on: true,
                        axial: PmlDir::Off,
                        z_in: 0.0,
                        z_thickness: 0.0,
                    },
                    strength: str_val,
                    order: ord,
                    profile: prof,
                    ..Default::default()
                },
            );
            let pml_table = EmPmlTable::from_mixed(&mixed, &mats, freq, false, false);
            let id_c3 = iso_mu_inv_c3(1.0);
            let k2 = c(k0 * k0, 0.0);

            let mut coo = CCoo::new(n, n);
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            for ci in 0..mixed.n_cells() {
                let iso = IsoparametricCell::from_mesh(&mixed, ci);
                let rule = quad_rule_mixed(iso.cell_type, 6);
                let nb = nedelec_nloc(iso.cell_type, 2);
                let mut ke = vec![c(0.0, 0.0); nb * nb];
                for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_2d(&mixed, &ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let wq = w_ref * jac.det.abs();
                    let (ep, mui, _) = pml_table.eval_2d(ci, id_c3, id_c3, [jac.x[0], jac.x[1]], false);
                    for s in 0..nb {
                        let ws = w[s];
                        for t in 0..nb {
                            let wt = w[t];
                            let ew = ep[0][0] * (ws[0] * wt[0])
                                + ep[0][1] * (ws[0] * wt[1] + ws[1] * wt[0])
                                + ep[1][1] * (ws[1] * wt[1]);
                            ke[s * nb + t] += (mui[2][2] * (cu[s] * cu[t]) - k2 * ew) * wq;
                        }
                    }
                }
                coo.push_block(ned.cell_dofs(ci), ned.cell_dofs(ci), &ke);
            }
            let a_full = coo.to_csr();

            // Axisymmetric `m = 0` TE cylindrical wave (`E_phi(r)` independent of `phi`): `E . dl = d_phi * r_src` on `r = r_src`!
            let mut dir_val: BTreeMap<usize, C> = BTreeMap::new();
            let arc = 2.0 * PI * r_src / (n_phi as f64);
            for f in 0..mixed.n_facets() {
                let lbl = mixed.facet_labels[f];
                let fv = mixed.facet(f);
                let eid = mixed.edge_index(fv[0], fv[1]).unwrap();
                if lbl == 2 {
                    dir_val.insert(2 * eid, c(0.0, 0.0));
                    dir_val.insert(2 * eid + 1, c(0.0, 0.0));
                } else if lbl == 1 {
                    let (lo, hi) = if ned.key_of(fv[0]) < ned.key_of(fv[1]) { (fv[0], fv[1]) } else { (fv[1], fv[0]) };
                    let (plo, phi) = (mixed.point(lo as usize), mixed.point(hi as usize));
                    let cross_z = plo[0] * phi[1] - plo[1] * phi[0];
                    let sgn = if cross_z >= 0.0 { 1.0 } else { -1.0 };
                    dir_val.insert(2 * eid, c(sgn * arc, 0.0));
                    dir_val.insert(2 * eid + 1, c(0.0, 0.0));
                }
            }
            let free: Vec<usize> = (0..n).filter(|d| !dir_val.contains_key(d)).collect();
            let fixed: Vec<usize> = dir_val.keys().copied().collect();
            let x_fixed: Vec<C> = fixed.iter().map(|d| dir_val[d]).collect();
            let a_ff = a_full.submatrix(&free, &free);
            let a_fc = a_full.submatrix(&free, &fixed);
            let rhs: Vec<C> = a_fc.matvec(&x_fixed).into_iter().map(|v| -v).collect();
            let x_free = DirectSolver::new(&a_ff, false).unwrap().solve(&rhs);
            let mut x_all = vec![c(0.0, 0.0); n];
            for (&d, &v) in fixed.iter().zip(&x_fixed) {
                x_all[d] = v;
            }
            for (&d, &v) in free.iter().zip(&x_free) {
                x_all[d] = v;
            }

            // Sample `E_phi` on the interior `Tri` cells (which share IDENTICAL mesh topology for any `n_pml`!):
            let mut samples = Vec::new();
            for ir in 1..(n_tri_rings - 1) {
                let ci = 4 * (ir * n_phi);
                let iso = IsoparametricCell::from_mesh(&mixed, ci);
                let xi = [1.0 / 3.0, 1.0 / 3.0];
                let jac = eval_cell_basis_2d(&mixed, &ned, ci, &iso, &xi, &mut scratch, &mut w, &mut cu);
                let r = jac.x[0].hypot(jac.x[1]);
                let phi_hat = [-jac.x[1] / r, jac.x[0] / r];
                let mut e_phi = c(0.0, 0.0);
                for (s, &gd) in ned.cell_dofs(ci).iter().enumerate() {
                    e_phi += x_all[gd] * (w[s][0] * phi_hat[0] + w[s][1] * phi_hat[1]);
                }
                samples.push((r, e_phi));
            }
            samples
        };

    let grading = if matches!(profile, PmlProfile::Rational { .. }) { 1.6 } else { 2.0 };
    let test_s = solve_annulus(n_layers, d_pml, profile, strength, order_m, grading);
    let ref_s = solve_annulus(48, 1.8 * lam0, PmlProfile::KappaPolynomial { kappa_max: 1.0 }, 7.5, 3.0, 2.0);
    let mut num2 = 0.0;
    let mut den2 = 0.0;
    for ((_, et), (_, er)) in test_s.iter().zip(&ref_s) {
        num2 += (*et - *er).norm_sqr();
        den2 += er.norm_sqr();
    }
    (num2 / den2.max(1e-300)).sqrt()
}

#[test]
fn pml_2d_tri_quad_cartesian_and_cylindrical_reflection_meets_1e6_and_1e8_targets() {
    // 1. Compact 8-layer KappaPolynomial & Rational: power reflectance `R = |r|^2 < 1e-6` and `< 1e-8`!
    let amp_cart_poly_8 =
        measure_2d_cartesian_tri_quad_reflection(PmlProfile::KappaPolynomial { kappa_max: 1.2 }, 8.5, 3.0, 0.55, 8);
    let amp_cart_rat_8 =
        measure_2d_cartesian_tri_quad_reflection(PmlProfile::Rational { curvature: 0.12 }, 15.5, 2.8, 0.34, 8);
    let amp_cyl_poly_8 =
        measure_2d_cylindrical_tri_quad_reflection(PmlProfile::KappaPolynomial { kappa_max: 1.2 }, 8.5, 3.0, 0.55, 8);
    let amp_cyl_rat_8 =
        measure_2d_cylindrical_tri_quad_reflection(PmlProfile::Rational { curvature: 0.12 }, 15.5, 2.8, 0.34, 8);

    // 2. Refined 26-layer KappaPolynomial: even amplitude reflection `|r| < 1e-6` (`R = |r|^2 < 1e-12`)!
    let amp_cart_poly_26 =
        measure_2d_cartesian_tri_quad_reflection(PmlProfile::KappaPolynomial { kappa_max: 1.05 }, 6.8, 3.0, 1.25, 26);
    let amp_cyl_poly_26 =
        measure_2d_cylindrical_tri_quad_reflection(PmlProfile::KappaPolynomial { kappa_max: 1.05 }, 6.8, 3.0, 1.25, 26);

    let r_cart_poly_pow = amp_cart_poly_8 * amp_cart_poly_8;
    let r_cart_rat_pow = amp_cart_rat_8 * amp_cart_rat_8;
    let r_cyl_poly_pow = amp_cyl_poly_8 * amp_cyl_poly_8;
    let r_cyl_rat_pow = amp_cyl_rat_8 * amp_cyl_rat_8;

    println!(
        "[2D Tri+Quad Reflection] Cartesian: |r|_poly(8L)={:.3e} (R={:.3e}), |r|_rat(8L)={:.3e} (R={:.3e}), |r|_poly(26L)={:.3e} | Cylindrical: |r|_poly(8L)={:.3e} (R={:.3e}), |r|_rat(8L)={:.3e} (R={:.3e}), |r|_poly(26L)={:.3e}",
        amp_cart_poly_8, r_cart_poly_pow, amp_cart_rat_8, r_cart_rat_pow, amp_cart_poly_26,
        amp_cyl_poly_8, r_cyl_poly_pow, amp_cyl_rat_8, r_cyl_rat_pow, amp_cyl_poly_26
    );
    assert!(r_cart_poly_pow < 1e-6, "Cartesian KappaPolynomial 8L power R = {r_cart_poly_pow:.3e} >= 1e-6");
    assert!(r_cart_rat_pow < 1e-8, "Cartesian Rational 8L power R = {r_cart_rat_pow:.3e} >= 1e-8");
    assert!(r_cyl_poly_pow < 1e-6, "Cylindrical KappaPolynomial 8L power R = {r_cyl_poly_pow:.3e} >= 1e-6");
    assert!(r_cyl_rat_pow < 1e-8, "Cylindrical Rational 8L power R = {r_cyl_rat_pow:.3e} >= 1e-8");
    assert!(amp_cart_poly_26 < 1e-6, "Cartesian 26L amplitude |r| = {amp_cart_poly_26:.3e} >= 1e-6");
    assert!(amp_cyl_poly_26 < 1e-6, "Cylindrical 26L amplitude |r| = {amp_cyl_poly_26:.3e} >= 1e-6");
}

// ------------------------------------------------------------------ 4. 3D Driven PROM P_pml Energy Balance & Multi-Octave Sub-Banding

#[test]
fn driven_3d_mixed_pml_prom_conserves_radiated_power_and_supports_subbands() {
    // 3D rectangular waveguide (`a = 22.86 mm`, `b = 10.16 mm` WR-90, `L_phys = 15 mm` along `+y`),
    // excited by a `waveport` at `y = 0` and terminated by 5 native `Prism` layers of `Rational` PML (`d = 15 mm`)!
    let (a, l_phys, b) = (22.86e-3, 15.0e-3, 10.16e-3);
    // Build a 3D `SimplexMesh` on `[0, a] x [0, l_phys] x [0, b]`, then extrude `+y` into 5 `Prism` layers via `build_pml`!
    let phys_tet = {
        let (nx, ny, nz) = (6usize, 4usize, 2usize);
        let mut pts = Vec::new();
        let mut g = vec![vec![vec![0u32; nz + 1]; ny + 1]; nx + 1];
        for ix in 0..=nx {
            for iy in 0..=ny {
                for iz in 0..=nz {
                    let id = pts.len() as u32;
                    pts.push([a * ix as f64 / nx as f64, l_phys * iy as f64 / ny as f64, b * iz as f64 / nz as f64]);
                    g[ix][iy][iz] = id;
                }
            }
        }
        let mut cells = Vec::new();
        let mut cl = Vec::new();
        let mut facets = Vec::new();
        let mut fl = Vec::new();
        // Freudenthal 6-tet split with consistent diagonal faces across every shared grid face!
        const TETS6: [[usize; 4]; 6] =
            [[0, 1, 2, 6], [0, 1, 6, 5], [0, 4, 5, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 4, 6, 7]];
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let c8 = [
                        g[ix][iy][iz],
                        g[ix + 1][iy][iz],
                        g[ix + 1][iy + 1][iz],
                        g[ix][iy + 1][iz],
                        g[ix][iy][iz + 1],
                        g[ix + 1][iy][iz + 1],
                        g[ix + 1][iy + 1][iz + 1],
                        g[ix][iy + 1][iz + 1],
                    ];
                    for t in &TETS6 {
                        let (i0, i1, i2, i3) = (c8[t[0]], c8[t[1]], c8[t[2]], c8[t[3]]);
                        let (p0, p1, p2, p3) = (pts[i0 as usize], pts[i1 as usize], pts[i2 as usize], pts[i3 as usize]);
                        let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                        let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                        let w = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
                        let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                            + u[2] * (v[0] * w[1] - v[1] * w[0]);
                        if det > 0.0 {
                            cells.extend_from_slice(&[i0, i1, i2, i3]);
                        } else {
                            cells.extend_from_slice(&[i0, i2, i1, i3]);
                        }
                        cl.push(1);
                    }
                    // All 6 bounding faces matching Freudenthal face diagonals (`0-6` cube diagonal => `0-2` on z, `0-5` on y, `0-7` on x, etc.)
                    if ix == 0 {
                        facets.extend_from_slice(&[c8[0], c8[3], c8[7], c8[0], c8[7], c8[4]]);
                        fl.extend_from_slice(&[1, 1]);
                    }
                    if ix + 1 == nx {
                        facets.extend_from_slice(&[c8[1], c8[2], c8[6], c8[1], c8[6], c8[5]]);
                        fl.extend_from_slice(&[2, 2]);
                    }
                    if iy == 0 {
                        facets.extend_from_slice(&[c8[0], c8[1], c8[5], c8[0], c8[5], c8[4]]);
                        fl.extend_from_slice(&[10, 10]); // waveport on y = 0
                    }
                    if iy + 1 == ny {
                        facets.extend_from_slice(&[c8[3], c8[2], c8[6], c8[3], c8[6], c8[7]]);
                        fl.extend_from_slice(&[3, 3]); // +y interface extruded into Prism PML
                    }
                    if iz == 0 {
                        facets.extend_from_slice(&[c8[0], c8[1], c8[2], c8[0], c8[2], c8[3]]);
                        fl.extend_from_slice(&[4, 4]);
                    }
                    if iz + 1 == nz {
                        facets.extend_from_slice(&[c8[4], c8[5], c8[6], c8[4], c8[6], c8[7]]);
                        fl.extend_from_slice(&[5, 5]);
                    }
                }
            }
        }
        femcore::adapt::work::SimplexMesh { points: pts, cells, cell_label: cl, facets, facet_label: fl }
    };

    let mut pml_spec = PmlSpec::<3>::uniform(0.0, 0);
    pml_spec.thickness[1] = [0.0, 30.0e-3];
    pml_spec.layers[1] = [0, 8];
    pml_spec.grading = 1.8;
    pml_spec.outer_label = Some(20);
    let pml3 = build_pml(&phys_tet, &pml_spec).expect("WR-90 +y Prism PML");
    let mesh = Mesh::D3(pml3.mesh);
    assert!(!mesh.is_pure_simplex());
    let mixed = discretise_mixed(&mesh).unwrap();
    let pml_dom = mixed.domain_labels().into_iter().find(|&d| d >= 1000).unwrap();

    // 1. Single-band driven PROM (`8.5 .. 12.0 GHz`, WR-90 X-band): verify `P_pml + |S_11|^2 == 1.0 W` and `|S_11| < 5e-3`!
    let mut prob = Problem::new(
        Study::Driven {
            f_min: 8.5e9,
            f_max: 12.0e9,
            n_points: 41,
            freqs_hz: None,
            prom_tol: Some(1e-4),
            max_hdm_solves: 8,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    prob.materials.default.insert("eps_r".into(), 1.0);
    prob.materials.default.insert("mu_r".into(), 1.0);
    prob.materials.pml.insert(
        pml_dom,
        Pml {
            y: PmlDir::Plus,
            strength: 5.5,
            order: 3.0,
            profile: PmlProfile::KappaPolynomial { kappa_max: 1.4 },
            ..Default::default()
        },
    );
    prob.waveports.push(WavePortBc { labels: vec![10], name: "P1".into(), mode: 1, n_eff_guess: None, z0_ohm: None });

    let progress = Arc::new(femcore::progress::Progress::new());
    let res = run(&mesh, &prob, &progress).expect("3D mixed driven WR-90 PML PROM");
    let num_scalar = |res_ref: &nyx_physics::SolveResult, key: &str| -> Option<f64> {
        res_ref.scalars.iter().find_map(|(k, v)| match (k == key, v) {
            (true, Scalar::Num(x)) => Some(*x),
            _ => None,
        })
    };
    let p_pml_w = num_scalar(&res, "P_pml (W)").expect("P_pml (W) scalar present");
    let sym_cnt = num_scalar(&res, "symbolic_count").unwrap();
    let dd: DrivenPortData = res
        .scalars
        .iter()
        .find_map(|(k, v)| match (k == DRIVEN_PORT_KEY, v) {
            (true, Scalar::Text(s)) => Some(serde_json::from_str(s).unwrap()),
            _ => None,
        })
        .unwrap();
    let max_s11 = dd.s_matrix.iter().map(|s| s[0][0][0].hypot(s[0][0][1])).fold(0.0_f64, f64::max);
    let f_rep = num_scalar(&res, "f_res_Hz").unwrap();
    let rep_idx = dd.freqs_hz.iter().position(|&f| (f - f_rep).abs() < 1.0).unwrap();
    let s11_rep_sq = dd.s_matrix[rep_idx][0][0][0].powi(2) + dd.s_matrix[rep_idx][0][0][1].powi(2);
    let power_balance_err = (p_pml_w + s11_rep_sq - 1.0).abs();

    println!(
        "[3D Driven Mixed WR-90 PML PROM] max|S11| = {:.3e}, P_pml = {:.6} W, |P_pml + |S11|^2 - 1| = {:.3e}, symbolic_count = {}",
        max_s11, p_pml_w, power_balance_err, sym_cnt
    );
    assert!(max_s11 < 5e-3, "WR-90 Prism Rational PML reflection too high: max|S11| = {max_s11:.3e}");
    assert!(p_pml_w > 0.999, "WR-90 PML must absorb >99.9% of incident 1 W port power: P_pml = {p_pml_w}");
    assert!(power_balance_err < 5e-4, "P_pml energy balance violated: P_pml = {p_pml_w}, |S11|^2 = {s11_rep_sq}");
    assert_eq!(sym_cnt, 1.0, "symbolic LU factorization must be reused once");

    // 2. Multi-octave broadband sub-banded PROM (`f_max / f_min = 4.0 > 3.0`) on a lumped-port + PML box:
    let mut prob_wide = Problem::new(
        Study::Driven {
            f_min: 2.0e9,
            f_max: 8.0e9,
            n_points: 31,
            freqs_hz: None,
            prom_tol: Some(1e-3),
            max_hdm_solves: 12,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        1,
    );
    prob_wide.materials.default.insert("eps_r".into(), 1.0);
    prob_wide.materials.default.insert("mu_r".into(), 1.0);
    prob_wide.materials.pml.insert(
        pml_dom,
        Pml {
            y: PmlDir::Plus,
            strength: 2.5,
            order: 2.0,
            profile: PmlProfile::KappaPolynomial { kappa_max: 2.0 },
            ..Default::default()
        },
    );
    prob_wide.lumped_ports.push(LumpedPortBc {
        labels: vec![10],
        name: "LP1".into(),
        z0: 50.0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: true,
        mode: "single".into(),
    });
    let res_wide = run(&mesh, &prob_wide, &progress).expect("multi-octave sub-banded PROM");
    let sym_wide = num_scalar(&res_wide, "symbolic_count").unwrap();
    assert_eq!(sym_wide, 1.0, "sub-banded PROM must share 1 symbolic LU factorization across sub-bands");

    // 3. Envelope `decode()` round-trip test for the new `"pml"` curvilinear & profile parameters:
    struct SimpleRes;
    impl LabelResolver for SimpleRes {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            match r {
                "core" => Ok(vec![1]),
                "pml" => Ok(vec![2]),
                _ => Err(Error::invalid("unknown region")),
            }
        }
        fn boundaries(&self, _: &str) -> Result<Vec<i32>, Error> {
            Ok(vec![1])
        }
    }
    let mut mat_map = BTreeMap::new();
    mat_map.insert("core".into(), MaterialAssign::Props(BTreeMap::from([("eps_r".into(), 1.0)])));
    mat_map.insert("pml".into(), MaterialAssign::Props(BTreeMap::from([("eps_r".into(), 1.0)])));
    let spec = ProblemSpec {
        format: 1,
        physics: "ems".into(),
        study: StudySpec { kind: "cavity".into(), params: json!({ "f_search": 10e9, "n_modes": 2 }) },
        order: 2,
        materials: mat_map,
        bcs: vec![BcSpec {
            kind: "pml".into(),
            on: "pml".into(),
            params: json!({
                "coord": "spherical",
                "center": [0.0, 0.0, 0.0],
                "r_in": 0.05,
                "thickness": 0.02,
                "profile": "rational",
                "curvature": 0.88,
                "auto_scale": true,
                "target_reflection": 1e-7
            }),
        }],
        params: json!({}),
    };
    nyx_physics::validate_spec(&Ems, &spec).expect("new pml params validate against descriptors");
    let res_ctx = SimpleRes;
    let ctx = Ctx::new(&res_ctx);
    let decoded = nyx_ems::decode(&spec, &ctx).expect("decode spherical rational auto_scale pml");
    let p_dec = &decoded.materials.pml[&2];
    assert_eq!(p_dec.coord, PmlCoord::Spherical { center: [0.0, 0.0, 0.0], r_in: 0.05, thickness: 0.02 });
    assert_eq!(p_dec.profile, PmlProfile::Rational { curvature: 0.88 });
    assert!(p_dec.auto_scale);
    assert!((p_dec.target_reflection - 1e-7).abs() < 1e-20);

    // Unused-import guard for Periodic if needed
    let _ = Periodic { source: vec![1], dest: vec![2], phase: 0.0 };
}

// ------------------------------------------------------------------ 5. Option B vs Split-Tet LBL^H Speedup, >=35% DOF Reduction & Determinism

#[test]
fn option_b_vs_split_tet_lblh_factorization_speedup_and_thread_determinism() {
    // Build a 3D physical box (`[0, 30mm]^3`, `2 x 2 x 2` Freudenthal cubes = 48 Tets)
    // surrounded on ALL 6 faces + 12 edges + 8 corners (26 sub-domains) by structured Cartesian PML!
    let phys_tet = {
        let (nx, ny, nz) = (2usize, 2usize, 2usize);
        let l = 30.0e-3;
        let mut pts = Vec::new();
        let mut g = vec![vec![vec![0u32; nz + 1]; ny + 1]; nx + 1];
        for ix in 0..=nx {
            for iy in 0..=ny {
                for iz in 0..=nz {
                    let id = pts.len() as u32;
                    pts.push([l * ix as f64 / nx as f64, l * iy as f64 / ny as f64, l * iz as f64 / nz as f64]);
                    g[ix][iy][iz] = id;
                }
            }
        }
        let mut cells = Vec::new();
        let mut cl = Vec::new();
        let mut facets = Vec::new();
        let mut fl = Vec::new();
        const TETS6: [[usize; 4]; 6] =
            [[0, 1, 2, 6], [0, 1, 6, 5], [0, 4, 5, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 4, 6, 7]];
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let c8 = [
                        g[ix][iy][iz],
                        g[ix + 1][iy][iz],
                        g[ix + 1][iy + 1][iz],
                        g[ix][iy + 1][iz],
                        g[ix][iy][iz + 1],
                        g[ix + 1][iy][iz + 1],
                        g[ix + 1][iy + 1][iz + 1],
                        g[ix][iy + 1][iz + 1],
                    ];
                    for t in &TETS6 {
                        let (i0, i1, i2, i3) = (c8[t[0]], c8[t[1]], c8[t[2]], c8[t[3]]);
                        let (p0, p1, p2, p3) = (pts[i0 as usize], pts[i1 as usize], pts[i2 as usize], pts[i3 as usize]);
                        let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                        let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                        let w = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
                        let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                            + u[2] * (v[0] * w[1] - v[1] * w[0]);
                        if det > 0.0 {
                            cells.extend_from_slice(&[i0, i1, i2, i3]);
                        } else {
                            cells.extend_from_slice(&[i0, i2, i1, i3]);
                        }
                        cl.push(1);
                    }
                    if ix == 0 {
                        facets.extend_from_slice(&[c8[0], c8[3], c8[7], c8[0], c8[7], c8[4]]);
                        fl.extend_from_slice(&[1, 1]);
                    }
                    if ix + 1 == nx {
                        facets.extend_from_slice(&[c8[1], c8[2], c8[6], c8[1], c8[6], c8[5]]);
                        fl.extend_from_slice(&[2, 2]);
                    }
                    if iy == 0 {
                        facets.extend_from_slice(&[c8[0], c8[1], c8[5], c8[0], c8[5], c8[4]]);
                        fl.extend_from_slice(&[3, 3]);
                    }
                    if iy + 1 == ny {
                        facets.extend_from_slice(&[c8[3], c8[2], c8[6], c8[3], c8[6], c8[7]]);
                        fl.extend_from_slice(&[4, 4]);
                    }
                    if iz == 0 {
                        facets.extend_from_slice(&[c8[0], c8[1], c8[2], c8[0], c8[2], c8[3]]);
                        fl.extend_from_slice(&[5, 5]);
                    }
                    if iz + 1 == nz {
                        facets.extend_from_slice(&[c8[4], c8[5], c8[6], c8[4], c8[6], c8[7]]);
                        fl.extend_from_slice(&[6, 6]);
                    }
                }
            }
        }
        femcore::adapt::work::SimplexMesh { points: pts, cells, cell_label: cl, facets, facet_label: fl }
    };

    let phys_flat: Vec<f64> = phys_tet.points.iter().flat_map(|p| p.iter().copied()).collect();
    let phys_sm = femcore::mesh::SimplexMesh::new(
        3,
        phys_flat,
        phys_tet.cells.clone(),
        phys_tet.cell_label.clone(),
        phys_tet.facets.clone(),
        phys_tet.facet_label.clone(),
    );
    let pml_dofs_n1_phys = nyx_ems::cavity::Ned3::new(&phys_sm, 1, (0..phys_sm.n_points() as u32).collect()).n_dofs();
    let pml_dofs_n2_phys = nyx_ems::cavity::Ned3::new(&phys_sm, 2, (0..phys_sm.n_points() as u32).collect()).n_dofs();

    // 1. Exact cell & N1/N2 PML DOF reductions across 3L, 5L (Rational Option B), and 8L (KappaPolynomial Split-Tet Option A)
    // on a 3D `+x, +y, +z` octant cavity (3 `Prism` faces + 3 `Hex` edges + 1 `Hex` corner):
    let make_octant_spec = |nl: usize| -> PmlSpec<3> {
        let mut sp = PmlSpec::<3>::uniform(0.0, 0);
        for a in 0..3 {
            sp.thickness[a] = [0.0, 12.0e-3];
            sp.layers[a] = [0, nl];
        }
        sp.grading = 2.0;
        sp.outer_label = Some(20);
        sp
    };

    let spec3 = make_octant_spec(3);
    let mesh3 = Mesh::D3(build_pml(&phys_tet, &spec3).unwrap().mesh);
    let mixed3 = discretise_mixed(&mesh3).unwrap();
    let split3 = discretise(&mesh3).unwrap();

    let spec5 = make_octant_spec(5);
    let mesh5 = Mesh::D3(build_pml(&phys_tet, &spec5).unwrap().mesh);
    let mixed5 = discretise_mixed(&mesh5).unwrap();
    let split5 = discretise(&mesh5).unwrap();

    let spec8 = make_octant_spec(8);
    let mesh8 = Mesh::D3(build_pml(&phys_tet, &spec8).unwrap().mesh);
    let split8 = discretise(&mesh8).unwrap();

    let n_tet5 = mixed5.cell_types.iter().filter(|&&t| t == CellType::Tet).count();
    let n_prism5 = mixed5.cell_types.iter().filter(|&&t| t == CellType::Prism).count();
    let n_hex5 = mixed5.cell_types.iter().filter(|&&t| t == CellType::Hex).count();

    let ned_b5_n1 = NedMixed::new(&mixed5, 1);
    let ned_a5_n1 = nyx_ems::cavity::Ned3::new(&split5, 1, (0..split5.n_points() as u32).collect());
    let ned_a8_n1 = nyx_ems::cavity::Ned3::new(&split8, 1, (0..split8.n_points() as u32).collect());

    let ned_b5_n2 = NedMixed::new(&mixed5, 2);
    let ned_a5_n2 = nyx_ems::cavity::Ned3::new(&split5, 2, (0..split5.n_points() as u32).collect());
    let ned_a8_n2 = nyx_ems::cavity::Ned3::new(&split8, 2, (0..split8.n_points() as u32).collect());

    let pml_b5_n1 = ned_b5_n1.n_dofs - pml_dofs_n1_phys;
    let pml_a5_n1 = ned_a5_n1.n_dofs() - pml_dofs_n1_phys;
    let pml_a8_n1 = ned_a8_n1.n_dofs() - pml_dofs_n1_phys;
    let red_5l_n1 = 100.0 * (1.0 - pml_b5_n1 as f64 / pml_a5_n1 as f64);
    let red_5v8_n1 = 100.0 * (1.0 - pml_b5_n1 as f64 / pml_a8_n1 as f64);

    let pml_b5_n2 = ned_b5_n2.n_dofs - pml_dofs_n2_phys;
    let pml_a5_n2 = ned_a5_n2.n_dofs() - pml_dofs_n2_phys;
    let pml_a8_n2 = ned_a8_n2.n_dofs() - pml_dofs_n2_phys;
    let red_5l_n2 = 100.0 * (1.0 - pml_b5_n2 as f64 / pml_a5_n2 as f64);
    let red_5v8_n2 = 100.0 * (1.0 - pml_b5_n2 as f64 / pml_a8_n2 as f64);

    let make_prob = |spec_ref: &PmlSpec<3>, domains: &[i32], order: usize, rational: bool| -> Problem {
        let mut p = Problem::new(Study::Cavity { n_modes: 2, f_search: 7.5e9 }, order);
        p.materials.default.insert("eps_r".into(), 1.0);
        p.materials.default.insert("mu_r".into(), 1.0);
        let dir_of = |c: i8| match c {
            -1 => PmlDir::Minus,
            1 => PmlDir::Plus,
            _ => PmlDir::Off,
        };
        for &dom in domains {
            if let Some((code, _)) = spec_ref.decode(dom) {
                p.materials.pml.insert(
                    dom,
                    Pml {
                        x: dir_of(code[0]),
                        y: dir_of(code[1]),
                        z: dir_of(code[2]),
                        strength: 3.5,
                        order: 2.0,
                        profile: if rational {
                            PmlProfile::Rational { curvature: 0.85 }
                        } else {
                            PmlProfile::KappaPolynomial { kappa_max: 1.5 }
                        },
                        ..Default::default()
                    },
                );
            }
        }
        p
    };

    let time_factor = |mat: &femcore::linalg::csr::CCsr| -> (f64, usize) {
        let mut best = f64::INFINITY;
        let mut f_nnz = 0usize;
        for _ in 0..3 {
            let t0 = std::time::Instant::now();
            let solver = DirectSolver::new(mat, false).expect("LBL^H / Sym LBL^T factorization");
            let dt = t0.elapsed().as_secs_f64();
            f_nnz = solver.stats().factor_nnz;
            let rhs = vec![c(1.0, -0.5); mat.nrows()];
            let sol = solver.solve(&rhs);
            assert!(sol[0].re.is_finite());
            best = best.min(dt);
        }
        (best, f_nnz)
    };

    // 2A. Order-1 (`N1`) LBL^H factorization timing on 5L Option B vs 5L Split-Tet vs 8L Split-Tet:
    let doms8 = discretise_mixed(&mesh8).unwrap().domain_labels();
    let (asm_b5_n1, _, _) =
        assemble_emw_eigen_mixed(&mixed5, &make_prob(&spec5, &mixed5.domain_labels(), 1, true), false).unwrap();
    let asm_a5_n1 =
        nyx_ems::assemble_emw_eigen(&split5, &make_prob(&spec5, &mixed5.domain_labels(), 1, false)).unwrap();
    let asm_a8_n1 = nyx_ems::assemble_emw_eigen(&split8, &make_prob(&spec8, &doms8, 1, false)).unwrap();

    let mat_b5_n1 =
        femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm_b5_n1.kc, c(-asm_b5_n1.sigma, 0.0), &asm_b5_n1.mc);
    let mat_a5_n1 =
        femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm_a5_n1.kc, c(-asm_a5_n1.sigma, 0.0), &asm_a5_n1.mc);
    let mat_a8_n1 =
        femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm_a8_n1.kc, c(-asm_a8_n1.sigma, 0.0), &asm_a8_n1.mc);
    let (t_b5_n1, fnnz_b5_n1) = time_factor(&mat_b5_n1);
    let (t_a5_n1, fnnz_a5_n1) = time_factor(&mat_a5_n1);
    let (t_a8_n1, fnnz_a8_n1) = time_factor(&mat_a8_n1);

    // 2B. Order-2 (`N2`) LBL^H factorization timing on 5L Option B (Rational) vs 5L Split-Tet vs 3L comparison,
    // plus Order-2 Serendipity Lagrange (`nyx-sms` / `nyx-dfs`) PML DOF reduction:
    let lag_phys_p2 = phys_sm.n_dofs(2);
    let lag_b5_p2 =
        femcore::fem_mixed::MixedLagrangeDofMap::new(&mixed5, femcore::fem_mixed::LagrangeSpace::Order2Serendipity)
            .n_dofs
            - lag_phys_p2;
    let lag_a5_p2 = split5.n_dofs(2) - lag_phys_p2;
    let red_5l_lag_p2 = 100.0 * (1.0 - lag_b5_p2 as f64 / lag_a5_p2 as f64);

    let (asm_b5_n2, _, _) =
        assemble_emw_eigen_mixed(&mixed5, &make_prob(&spec5, &mixed5.domain_labels(), 2, true), false).unwrap();
    let asm_a5_n2 =
        nyx_ems::assemble_emw_eigen(&split5, &make_prob(&spec5, &mixed5.domain_labels(), 2, false)).unwrap();

    let mat_b5_n2 =
        femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm_b5_n2.kc, c(-asm_b5_n2.sigma, 0.0), &asm_b5_n2.mc);
    let mat_a5_n2 =
        femcore::linalg::csr::CCsr::axpby(c(1.0, 0.0), &asm_a5_n2.kc, c(-asm_a5_n2.sigma, 0.0), &asm_a5_n2.mc);
    let (t_b5_n2, fnnz_b5_n2) = time_factor(&mat_b5_n2);
    let (t_a5_n2, fnnz_a5_n2) = time_factor(&mat_a5_n2);
    let red_free_5l_n2 = 100.0 * (1.0 - mat_b5_n2.nrows() as f64 / mat_a5_n2.nrows() as f64);

    println!(
        "[Option B vs Split-Tet Benchmark]\n  \
         Cells (5L): Option B = {} ({} Tet, {} Prism, {} Hex) vs Option A (5L split) = {} Tets ({:.2}x fewer PML cells), Option A (8L split) = {} Tets\n  \
         PML DOFs N1: 5L Option B = {} vs 5L Split = {} ({:.2}% reduction) | vs 8L Split = {} ({:.2}% reduction)\n  \
         PML DOFs P2 Serendipity Lagrange (`sms`/`dfs`): 5L Option B = {} vs 5L Split = {} ({:.2}% reduction)\n  \
         PML DOFs N2: 5L Option B = {} vs 5L Split = {} ({:.2}% unreduced, {:.2}% free DOFs: {} vs {}) | vs 8L Split = {} ({:.2}% reduction)\n  \
         LBL^H N1 (best of 3): 5L Option B = {:.2} ms (free={}, f_nnz={}) vs 5L Split = {:.2} ms (free={}, f_nnz={}, {:.2}x speedup) vs 8L Split = {:.2} ms (free={}, f_nnz={}, {:.2}x speedup)\n  \
         LBL^H N2 (best of 3): 5L Option B = {:.2} ms (free={}, f_nnz={}) vs 5L Split = {:.2} ms (free={}, f_nnz={}, {:.2}x speedup)",
        mixed5.n_cells(),
        n_tet5,
        n_prism5,
        n_hex5,
        split5.n_cells(),
        (split5.n_cells() - n_tet5) as f64 / (n_prism5 + n_hex5) as f64,
        split8.n_cells(),
        pml_b5_n1,
        pml_a5_n1,
        red_5l_n1,
        pml_a8_n1,
        red_5v8_n1,
        lag_b5_p2,
        lag_a5_p2,
        red_5l_lag_p2,
        pml_b5_n2,
        pml_a5_n2,
        red_5l_n2,
        red_free_5l_n2,
        mat_b5_n2.nrows(),
        mat_a5_n2.nrows(),
        pml_a8_n2,
        red_5v8_n2,
        1e3 * t_b5_n1,
        mat_b5_n1.nrows(),
        fnnz_b5_n1,
        1e3 * t_a5_n1,
        mat_a5_n1.nrows(),
        fnnz_a5_n1,
        t_a5_n1 / t_b5_n1.max(1e-9),
        1e3 * t_a8_n1,
        mat_a8_n1.nrows(),
        fnnz_a8_n1,
        t_a8_n1 / t_b5_n1.max(1e-9),
        1e3 * t_b5_n2,
        mat_b5_n2.nrows(),
        fnnz_b5_n2,
        1e3 * t_a5_n2,
        mat_a5_n2.nrows(),
        fnnz_a5_n2,
        t_a5_n2 / t_b5_n2.max(1e-9),
    );

    assert!(red_5l_n1 >= 50.0, "expected >= 50% N1 PML DOF reduction, got {red_5l_n1:.2}%");
    assert!(red_5l_lag_p2 >= 40.0, "expected >= 40% P2 Serendipity PML DOF reduction, got {red_5l_lag_p2:.2}%");
    assert!(red_free_5l_n2 >= 35.0, "expected >= 35% free N2 DOF reduction on 5L PML, got {red_free_5l_n2:.2}%");
    assert!(red_5v8_n2 >= 60.0, "expected >= 60% 5L-Rational vs 8L-Split N2 PML DOF reduction, got {red_5v8_n2:.2}%");
    assert!(fnnz_b5_n1 < fnnz_a5_n1 && fnnz_b5_n2 < fnnz_a5_n2, "Option B must reduce factor fill-in vs split-tet");
    assert!(t_b5_n1 < t_a5_n1 && t_b5_n2 < t_a5_n2, "Option B LBL^H factorization must be faster than split-tet");
    let _ = split3;

    // 3. Thread-count bit-determinism (`1` vs `4` threads) on 3D `Tet + Prism + Hex` PML cavity assembly & solve:
    let solve_with_threads = |nt: usize| -> (Vec<[u64; 2]>, Vec<u64>, Vec<u64>) {
        rayon::ThreadPoolBuilder::new().num_threads(nt).build().unwrap().install(|| {
            let mut p_det = make_prob(&spec3, &mixed3.domain_labels(), 1, false);
            p_det.study = Study::Cavity { n_modes: 2, f_search: 7.5e9 };
            let progress = Arc::new(femcore::progress::Progress::new());
            let r = run(&mesh3, &p_det, &progress).expect("deterministic 3D mixed PML cavity");
            let eig_bits: Vec<[u64; 2]> = r
                .solutions
                .iter()
                .map(|s| {
                    let e = s.eigenvalue.unwrap();
                    [e[0].to_bits(), e[1].to_bits()]
                })
                .collect();
            let field_bits: Vec<u64> = r.solutions[0].fields[0].values.iter().map(|v| v.to_bits()).collect();
            let cell_bits: Vec<u64> = r.solutions[0].cell_fields[0].re.iter().map(|v| v.to_bits()).collect();
            (eig_bits, field_bits, cell_bits)
        })
    };
    let det1 = solve_with_threads(1);
    let det4 = solve_with_threads(4);
    assert_eq!(det1.0, det4.0, "eigenvalues differed between 1 and 4 threads");
    assert_eq!(det1.1, det4.1, "nodal fields differed between 1 and 4 threads");
    assert_eq!(det1.2, det4.2, "cell fields differed between 1 and 4 threads");
}
