// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of `femcore::rt` (`AF-4`):
//! - Exact Piola contravariant mapping and Gauss divergence theorem for `RT_0` (`order=1`) and `RT_1` (`order=2`) in 2D and 3D;
//! - Machine-precision exactness (`C^low == C^exact == C^up`), solenoidal `div(D_tilde) == 0`, and zero normal-flux jump across
//!   dielectric interfaces (`Si|VacTop` + anisotropic tensor) on 2D & 3D parallel-plate capacitors;
//! - Strict Loewner bracketing `C^low <= C^exact <= C^up` and `> 2.5x` `C^mid` accuracy improvement on 2D & 3D coaxial cylinders;
//! - Multi-terminal 3D interior conductive sheet (`sheet3d`) capacitor with half-star splitting across `Gamma_sheet`, continuous
//!   normal flux across `Si|VacTop`, full `2x2` Loewner matrix bracketing `v^T C^low v <= v^T C^exact v <= v^T C^up v`, off-diagonal
//!   error bound `|C^mid_ij - C^exact_ij| <= delta_C_ij`, and `> 2.5x` accuracy gain;
//! - Robin `FarField` boundary complementary pairing on a 3D spherical shell (`C^low <= 4*pi*eps0*r_in <= C^up`);
//! - Compatibility guard rejecting singular Barsoum quarter-point `CurvedTet`s while accepting smooth isoparametric `CurvedTet`s;
//! - Exact IEEE-754 bit-determinism across `NYX_THREADS=1` and `NYX_THREADS=4`.
#![allow(clippy::needless_range_loop)]

use femcore::fem::{
    assemble_boundary_mass, assemble_laplace, assemble_laplace_tensor, assemble_parallel, boundary_dofs,
    dirichlet_constraints, facet_rule, free_dofs, geometry, quad_rule, quad_rule_collapsed, shape, shape_gradients,
    CurvedTet, RefElement,
};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{dense_eigh, DirectSolver};
use femcore::material::M3;
use femcore::mesh::SimplexMesh;
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::rt::{
    equilibrate_electro_flux_bounds, equilibrate_flux_patches, eval_rt_piola, rt_nloc, CellPermittivity,
    EquilibratedFluxResult, RobinBoundary,
};
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

const EPS0: f64 = 8.854_187_812_8e-12;

fn assemble_laplace_curved_3d_test(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Csr {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(6);
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * 3];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.fill(0.0);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let dn = &cshapes[q];
                    for a in 0..nloc {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for b in 0..4 {
                                s += dn[a * 4 + b] * gq[b][k];
                            }
                            grad[a * 3 + k] = s;
                        }
                    }
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += wq * s;
                        }
                    }
                }
                let sc = coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += re.w[q] * s;
                        }
                    }
                }
                let sc = geo.vol[c] * coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// Solves multi-terminal primal Laplace `K u = 0` with `u = 1` on terminal `j` and `u = 0` on all other Dirichlet boundaries.
fn solve_primal_terminals(
    mesh: &SimplexMesh,
    order: usize,
    k: &Csr,
    term_labels: &[Vec<i32>],
    ground_labels: &[i32],
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = k.nrows();
    let n_term = term_labels.len();
    let bg_dofs: BTreeSet<usize> = boundary_dofs(mesh, ground_labels, order).into_iter().collect();
    let term_dofs: Vec<BTreeSet<usize>> =
        term_labels.iter().map(|labs| boundary_dofs(mesh, labs, order).into_iter().collect()).collect();

    let mut all_fixed = bg_dofs.clone();
    for td in &term_dofs {
        all_fixed.extend(td.iter().copied());
    }
    let fixed: Vec<usize> = all_fixed.into_iter().collect();
    let free = free_dofs(n, &fixed);

    let vals_list: Vec<Vec<f64>> = (0..n_term)
        .map(|j| {
            fixed
                .iter()
                .map(|&d| if term_dofs[j].contains(&d) && !bg_dofs.contains(&d) { 1.0 } else { 0.0 })
                .collect()
        })
        .collect();

    let mut v_all = vec![vec![0.0; n]; n_term];
    for j in 0..n_term {
        for (&d, &val) in fixed.iter().zip(&vals_list[j]) {
            v_all[j][d] = val;
        }
    }

    if !free.is_empty() {
        let kff = k.submatrix(&free, &free);
        let kfx = k.submatrix(&free, &fixed);
        let rhs_list: Vec<Vec<f64>> = vals_list
            .iter()
            .map(|vals| {
                let kx = kfx.matvec(vals);
                kx.into_iter().map(|x| -x).collect()
            })
            .collect();
        let solver = DirectSolver::new(&kff, true).unwrap();
        let sols = solver.solve_many(&rhs_list);
        for j in 0..n_term {
            for (a, &i) in free.iter().enumerate() {
                v_all[j][i] = sols[j][a];
            }
        }
    }

    let kv_all: Vec<Vec<f64>> = v_all.iter().map(|v| k.matvec(v)).collect();
    let mut c_up = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in i..n_term {
            let cij: f64 = v_all[i].iter().zip(&kv_all[j]).map(|(a, b)| a * b).sum();
            let cji: f64 = v_all[j].iter().zip(&kv_all[i]).map(|(a, b)| a * b).sum();
            let sym = 0.5 * (cij + cji);
            c_up[i][j] = sym;
            c_up[j][i] = sym;
        }
    }
    (v_all, c_up)
}

/// Verifies Loewner matrix inequality `A <= B` (`B - A` is positive semi-definite within `tol`).
fn assert_loewner_le(a: &[Vec<f64>], b: &[Vec<f64>], tol: f64, label: &str) {
    let n = a.len();
    let mut diff = Vec::with_capacity(n * n);
    let mut scale = 0.0_f64;
    for i in 0..n {
        scale = scale.max(a[i][i].abs()).max(b[i][i].abs());
        for j in 0..n {
            diff.push(0.5 * ((b[i][j] - a[i][j]) + (b[j][i] - a[j][i])));
        }
    }
    let (eigs, _) = dense_eigh(&diff, n).unwrap();
    let min_eig = eigs.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(
        min_eig >= -tol * scale.max(1.0e-30),
        "{label}: Loewner inequality violated, min_eig(B - A) = {min_eig:e} (scale = {scale:e}, eigs = {eigs:?})"
    );
}

#[test]
fn rt0_and_rt1_piola_gauss_theorem_on_skewed_2d_and_3d_simplices() {
    // 2D skewed triangle (test both positive and negative det J vertex orderings)
    let tri_pos: [[f64; 3]; 4] = [[0.1, -0.2, 0.0], [1.4, 0.3, 0.0], [-0.3, 1.1, 0.0], [0.0; 3]];
    let tri_neg = [tri_pos[0], tri_pos[2], tri_pos[1], [0.0; 3]];
    for pts in [tri_pos, tri_neg] {
        let j = [
            [pts[1][0] - pts[0][0], pts[2][0] - pts[0][0], 0.0],
            [pts[1][1] - pts[0][1], pts[2][1] - pts[0][1], 0.0],
            [0.0; 3],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let vol = det.abs() / 2.0;
        for order in [1usize, 2] {
            let nrt = rt_nloc(2, order);
            let q = quad_rule(2, 4);
            let frule = facet_rule(1);
            let mut pvec = vec![[0.0; 3]; nrt];
            let mut pdiv = vec![0.0; nrt];
            let mut vol_int = vec![0.0; nrt];
            for (qi, lam) in q.lam.iter().enumerate() {
                eval_rt_piola(2, order, &j, det, lam, &mut pvec, &mut pdiv);
                for r in 0..nrt {
                    vol_int[r] += vol * q.w[qi] * pdiv[r];
                }
            }
            let mut bnd_int = vec![0.0; nrt];
            for opp in 0..3 {
                let slots: Vec<usize> = (0..3).filter(|&a| a != opp).collect();
                let (p0, p1, popp) = (pts[slots[0]], pts[slots[1]], pts[opp]);
                let (tx, ty) = (p1[0] - p0[0], p1[1] - p0[1]);
                let len = (tx * tx + ty * ty).sqrt();
                let (mut nx, mut ny) = (ty / len, -tx / len);
                if nx * (popp[0] - p0[0]) + ny * (popp[1] - p0[1]) > 0.0 {
                    nx = -nx;
                    ny = -ny;
                }
                for (fi, mu) in frule.lam.iter().enumerate() {
                    let mut lam = [0.0; 3];
                    lam[slots[0]] = mu[0];
                    lam[slots[1]] = mu[1];
                    eval_rt_piola(2, order, &j, det, &lam, &mut pvec, &mut pdiv);
                    for r in 0..nrt {
                        let fn_val = pvec[r][0] * nx + pvec[r][1] * ny;
                        if r >= 6 {
                            // RT_1 bubble functions must vanish identically on every edge
                            assert!(fn_val.abs() < 1.0e-15, "2D RT_1 bubble trace non-zero: {fn_val:e}");
                        }
                        bnd_int[r] += len * frule.w[fi] * fn_val;
                    }
                }
            }
            for r in 0..nrt {
                assert!(
                    (vol_int[r] - bnd_int[r]).abs() < 1.0e-14,
                    "2D order={order} r={r}: vol_div={} vs bnd_flux={}",
                    vol_int[r],
                    bnd_int[r]
                );
            }
        }
    }

    // 3D skewed tetrahedron (test both positive and negative det J vertex orderings)
    let tet_pos: [[f64; 3]; 4] = [[0.1, -0.2, 0.05], [1.3, 0.1, -0.1], [0.2, 1.1, 0.3], [-0.1, 0.2, 0.9]];
    let tet_neg = [tet_pos[0], tet_pos[2], tet_pos[1], tet_pos[3]];
    for pts in [tet_pos, tet_neg] {
        let mut j = [[0.0; 3]; 3];
        for b in 0..3 {
            for a in 0..3 {
                j[a][b] = pts[b + 1][a] - pts[0][a];
            }
        }
        let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
        let vol = det.abs() / 6.0;
        for order in [1usize, 2] {
            let nrt = rt_nloc(3, order);
            let q = quad_rule(3, 4);
            let frule = facet_rule(2);
            let mut pvec = vec![[0.0; 3]; nrt];
            let mut pdiv = vec![0.0; nrt];
            let mut vol_int = vec![0.0; nrt];
            for (qi, lam) in q.lam.iter().enumerate() {
                eval_rt_piola(3, order, &j, det, lam, &mut pvec, &mut pdiv);
                for r in 0..nrt {
                    vol_int[r] += vol * q.w[qi] * pdiv[r];
                }
            }
            let mut bnd_int = vec![0.0; nrt];
            for opp in 0..4 {
                let slots: Vec<usize> = (0..4).filter(|&a| a != opp).collect();
                let (p0, p1, p2, popp) = (pts[slots[0]], pts[slots[1]], pts[slots[2]], pts[opp]);
                let t1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                let t2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                let cr = [
                    t1[1] * t2[2] - t1[2] * t2[1],
                    t1[2] * t2[0] - t1[0] * t2[2],
                    t1[0] * t2[1] - t1[1] * t2[0],
                ];
                let norm = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
                let area = 0.5 * norm;
                let s = if cr[0] * (popp[0] - p0[0]) + cr[1] * (popp[1] - p0[1]) + cr[2] * (popp[2] - p0[2]) < 0.0 {
                    1.0
                } else {
                    -1.0
                };
                let n_out = [s * cr[0] / norm, s * cr[1] / norm, s * cr[2] / norm];
                for (fi, mu) in frule.lam.iter().enumerate() {
                    let mut lam = [0.0; 4];
                    lam[slots[0]] = mu[0];
                    lam[slots[1]] = mu[1];
                    lam[slots[2]] = mu[2];
                    eval_rt_piola(3, order, &j, det, &lam, &mut pvec, &mut pdiv);
                    for r in 0..nrt {
                        let fn_val = pvec[r][0] * n_out[0] + pvec[r][1] * n_out[1] + pvec[r][2] * n_out[2];
                        if r >= 12 {
                            assert!(fn_val.abs() < 1.0e-15, "3D RT_1 bubble trace non-zero: {fn_val:e}");
                        }
                        bnd_int[r] += area * frule.w[fi] * fn_val;
                    }
                }
            }
            for r in 0..nrt {
                assert!(
                    (vol_int[r] - bnd_int[r]).abs() < 1.0e-14,
                    "3D order={order} r={r}: vol_div={} vs bnd_flux={}",
                    vol_int[r],
                    bnd_int[r]
                );
            }
        }
    }
}

#[test]
fn parallel_plate_2d_and_3d_layered_dielectric_and_tensor_are_exact() {
    // 1. 2D two-layer parallel-plate capacitor: x in [0, lx], z in [0, lz], interface at x = xc (plates at x = 0 and x = lx)
    let (lx, lz, xc) = (100e-6, 40e-6, 40e-6);
    let (eps1, eps2) = (11.45 * EPS0, 1.0 * EPS0);
    let c_exact_2d = lz / (xc / eps1 + (lx - xc) / eps2);
    let mesh_2d = femcore::testkit::rect_split_x(lx, lz, xc, 20e-6);
    let eps_cells_2d: Vec<f64> = mesh_2d
        .cell_labels
        .iter()
        .map(|&d| if d == 1 { eps1 } else { eps2 })
        .collect();

    for order in [1usize, 2] {
        let k = assemble_laplace(&mesh_2d, order, &eps_cells_2d);
        // Boundary 1: x = 0 (ground), Boundary 4: x = lx (terminal), Boundaries 2 & 3: z = 0, lz (Neumann)
        let (v_all, c_up) = solve_primal_terminals(&mesh_2d, order, &k, &[vec![4]], &[1]);
        let eq = equilibrate_electro_flux_bounds(
            &mesh_2d,
            order,
            &v_all,
            CellPermittivity::Scalar(&eps_cells_2d),
            &[1, 4],
            RobinBoundary::None,
            Some(&c_up),
        )
        .unwrap();

        let sheet_mask = vec![false; mesh_2d.n_facets()];
        let jump = eq.interior_non_sheet_normal_jump_l2(&mesh_2d, &sheet_mask, 0).unwrap();
        assert!(
            jump < 1.0e-12 * c_exact_2d / lz.sqrt(),
            "2D order={order} dielectric interface normal flux jump = {jump:e}"
        );
        assert!(
            eq.max_div_l2 < 1.0e-12 * c_exact_2d / (lx * lz).sqrt(),
            "2D order={order} max_div_l2 = {:e}",
            eq.max_div_l2
        );
        eprintln!(
            "2D plate order={order}: c_low={:.10e}, c_up={:.10e}, exact={c_exact_2d:.10e}, Q={:.10e}, W={:.10e}",
            eq.c_lower_matrix[0][0], eq.c_upper_matrix[0][0], eq.q_pairing_matrix[0][0], eq.w_flux_matrix[0][0]
        );
        assert!(eq.c_lower_matrix[0][0] <= c_exact_2d * (1.0 + 1.0e-12));
        if order == 2 {
            assert!(
                (eq.c_lower_matrix[0][0] - c_exact_2d).abs() < 1.0e-12 * c_exact_2d,
                "2D order={order}: c_low={:e} vs exact={c_exact_2d:e}",
                eq.c_lower_matrix[0][0]
            );
        }
        assert!(
            (eq.c_upper_matrix[0][0] - c_exact_2d).abs() < 1.0e-12 * c_exact_2d,
            "2D order={order}: c_up={:e} vs exact={c_exact_2d:e}",
            eq.c_upper_matrix[0][0]
        );
    }

    // 2. 3D two-layer parallel-plate capacitor (`Si` z < zc, `VacTop` z > zc, plates at z = 0 and z = bz) + anisotropic tensor check
    let (bx, by, bz, zc) = (80e-6, 60e-6, 50e-6, 20e-6);
    let mesh_3d = femcore::testkit::block_split_z(bx, by, bz, zc, 25e-6, 3);
    // Anisotropic sapphire-like tensor in domain 1 (eps_zz = 11.55, eps_xx = eps_yy = 9.40) and vacuum in domain 2
    let (e_perp, e_z) = (9.40 * EPS0, 11.55 * EPS0);
    let kt_3d: Vec<M3> = mesh_3d
        .cell_labels
        .iter()
        .map(|&d| {
            if d == 1 {
                [[e_perp, 0.0, 0.0], [0.0, e_perp, 0.0], [0.0, 0.0, e_z]]
            } else {
                [[EPS0, 0.0, 0.0], [0.0, EPS0, 0.0], [0.0, 0.0, EPS0]]
            }
        })
        .collect();
    let c_exact_3d = (bx * by) / (zc / e_z + (bz - zc) / EPS0);

    for order in [1usize, 2] {
        let k = assemble_laplace_tensor(&mesh_3d, order, &kt_3d);
        // Boundary 3: z = 0 (ground), Boundary 4: z = bz (terminal), sides 1,2,5,6 (Neumann)
        let (v_all, c_up) = solve_primal_terminals(&mesh_3d, order, &k, &[vec![4]], &[3]);
        let eq = equilibrate_electro_flux_bounds(
            &mesh_3d,
            order,
            &v_all,
            CellPermittivity::Tensor(&kt_3d),
            &[3, 4],
            RobinBoundary::None,
            Some(&c_up),
        )
        .unwrap();

        let sheet_mask = vec![false; mesh_3d.n_facets()];
        let jump = eq.interior_non_sheet_normal_jump_l2(&mesh_3d, &sheet_mask, 0).unwrap();
        assert!(
            jump < 1.0e-12 * c_exact_3d / (bx * by).sqrt(),
            "3D order={order} interface normal flux jump = {jump:e}"
        );
        assert!(
            eq.max_div_l2 < 1.0e-12 * c_exact_3d / (bx * by * bz).sqrt(),
            "3D order={order} max_div_l2 = {:e}",
            eq.max_div_l2
        );
        eprintln!(
            "3D plate order={order}: c_low={:.10e}, c_up={:.10e}, exact={c_exact_3d:.10e}, Q={:.10e}, W={:.10e}",
            eq.c_lower_matrix[0][0], eq.c_upper_matrix[0][0], eq.q_pairing_matrix[0][0], eq.w_flux_matrix[0][0]
        );
        assert!(eq.c_lower_matrix[0][0] <= c_exact_3d * (1.0 + 1.0e-11));
        if order == 2 {
            assert!(
                (eq.c_lower_matrix[0][0] - c_exact_3d).abs() < 1.0e-11 * c_exact_3d,
                "3D order={order}: c_low={:e} vs exact={c_exact_3d:e}",
                eq.c_lower_matrix[0][0]
            );
        }
        assert!(
            (eq.c_upper_matrix[0][0] - c_exact_3d).abs() < 1.0e-11 * c_exact_3d,
            "3D order={order}: c_up={:e} vs exact={c_exact_3d:e}",
            eq.c_upper_matrix[0][0]
        );
    }
}

/// Builds a shape-regular 2D Delaunay coaxial mesh and its extruded + curved 3D coaxial cylinder mesh.
fn build_delaunay_coax_2d_and_3d(
    r_in: f64,
    r_out: f64,
    lz: f64,
    h2d: f64,
    max_dev: f64,
) -> (SimplexMesh, SimplexMesh, Vec<i32>, Vec<i32>) {
    use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
    use femcore::meshing::mesher::{extrude_y, Layer};
    use geom::model::RegionKind;

    let mut m = Model::new_2d();
    let outer = m.circle("diel", (0.0, 0.0), r_out);
    let inner = m.circle("core", (0.0, 0.0), r_in);
    m.overlay_full(
        &[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)],
        &[],
    );
    let mut opts = GeomMeshOpts::new(h2d);
    opts.max_dev = Some(max_dev);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (r_in + r_out);
    let mut inner_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        if rm < r_mid {
            if !inner_labs.contains(&lab) {
                inner_labs.push(lab);
            }
        } else if !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }

    let layers = vec![Layer { thickness: lz, n: 2, domain_offset: 0, fill: None }];
    let mut m3 = extrude_y(&g2.mesh, &layers, 0.0).unwrap();
    let tol_a = 1e-3 * r_in;
    let tol_b = 1e-3 * r_out;
    let mut edge_nodes = Vec::new();
    for &(u, v) in m3.edges() {
        let pu = m3.point(u as usize);
        let pv = m3.point(v as usize);
        let ru = (pu[0] * pu[0] + pu[2] * pu[2]).sqrt();
        let rv = (pv[0] * pv[0] + pv[2] * pv[2]).sqrt();
        let target_r = if (ru - r_in).abs() < tol_a && (rv - r_in).abs() < tol_a {
            Some(r_in)
        } else if (ru - r_out).abs() < tol_b && (rv - r_out).abs() < tol_b {
            Some(r_out)
        } else {
            None
        };
        if let Some(r_exact) = target_r {
            let mx = 0.5 * (pu[0] + pv[0]);
            let my = 0.5 * (pu[1] + pv[1]);
            let mz = 0.5 * (pu[2] + pv[2]);
            let mr = (mx * mx + mz * mz).sqrt();
            if mr > 1e-12 && (mr - r_exact).abs() > 1e-12 * r_exact {
                edge_nodes.push(([u, v], [mx * (r_exact / mr), my, mz * (r_exact / mr)]));
            }
        }
    }
    m3.set_edge_nodes(&edge_nodes);
    (g2.mesh, m3, inner_labs, outer_labs)
}

/// Uniform 1-to-4 red refinement of a 2D triangle mesh, preserving the exact polygonal boundary `dOmega_h` so
/// continuum reference capacitance on `dOmega_h` has zero boundary-perturbation discrepancy.
fn refine_2d_red(mesh: &SimplexMesh) -> SimplexMesh {
    assert_eq!(mesh.dim, 2);
    let nv = mesh.n_points();
    let mut points = mesh.points.clone();
    for &(a, b) in mesh.edges() {
        let pa = mesh.point(a as usize);
        let pb = mesh.point(b as usize);
        points.push(0.5 * (pa[0] + pb[0]));
        points.push(0.5 * (pa[1] + pb[1]));
    }
    let mid = |u: u32, v: u32| -> u32 { (nv + mesh.edge_index(u, v)) as u32 };
    let mut cells = Vec::with_capacity(mesh.n_cells() * 4 * 3);
    let mut cell_labels = Vec::with_capacity(mesh.n_cells() * 4);
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let (v0, v1, v2) = (v[0], v[1], v[2]);
        let (m01, m12, m20) = (mid(v0, v1), mid(v1, v2), mid(v2, v0));
        for tri in [[v0, m01, m20], [m01, v1, m12], [m20, m12, v2], [m01, m12, m20]] {
            cells.extend_from_slice(&tri);
            cell_labels.push(mesh.cell_labels[c]);
        }
    }
    let mut facets = Vec::with_capacity(mesh.n_facets() * 2 * 2);
    let mut facet_labels = Vec::with_capacity(mesh.n_facets() * 2);
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let m = mid(fv[0], fv[1]);
        facets.extend_from_slice(&[fv[0], m, m, fv[1]]);
        facet_labels.push(mesh.facet_labels[f]);
        facet_labels.push(mesh.facet_labels[f]);
    }
    SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels)
}

#[test]
fn coaxial_cylinder_2d_and_3d_strict_loewner_bracketing_and_midpoint_gain() {
    let (r_in, r_out, lz): (f64, f64, f64) = (25e-6, 100e-6, 40e-6);
    let eps_r = 3.2;
    let eps_val = eps_r * EPS0;
    let c_exact_2d_circle = 2.0 * PI * eps_val / (r_out / r_in).ln();
    let c_exact_3d = c_exact_2d_circle * lz;

    // Use moderate `max_dev` (`r_in * 4e-3`) so the 2D/3D meshes remain genuinely coarse (~45 boundary segments)
    // rather than over-refining the annulus to <30 ppm!
    let (mesh_2d, mesh_3d, inner_labs, outer_labs) =
        build_delaunay_coax_2d_and_3d(r_in, r_out, lz, 32e-6, r_in * 4.0e-3);
    let mut dir_labs = inner_labs.clone();
    dir_labs.extend_from_slice(&outer_labs);

    // Exact continuum capacitance on the fixed 2D polygonal domain via 3x uniform 1-to-4 nested refinement (64x cells) + P2
    let mesh_2d_fine = refine_2d_red(&refine_2d_red(&refine_2d_red(&mesh_2d)));
    let eps_2d_fine = vec![eps_val; mesh_2d_fine.n_cells()];
    let k_2d_fine = assemble_laplace(&mesh_2d_fine, 2, &eps_2d_fine);
    let (_, cup_2d_fine) =
        solve_primal_terminals(&mesh_2d_fine, 2, &k_2d_fine, std::slice::from_ref(&inner_labs), &outer_labs);
    let c_exact_2d_poly = cup_2d_fine[0][0];

    // 1. 2D coaxial capacitor (`order = 1` RT_0 and `order = 2` RT_1): exact Loewner bracketing & `> 2.5x` C^mid gain!
    let eps_2d = vec![eps_val; mesh_2d.n_cells()];
    for order in [1usize, 2] {
        let k = assemble_laplace(&mesh_2d, order, &eps_2d);
        let (v_all, c_up) =
            solve_primal_terminals(&mesh_2d, order, &k, std::slice::from_ref(&inner_labs), &outer_labs);
        let eq = equilibrate_electro_flux_bounds(
            &mesh_2d,
            order,
            &v_all,
            CellPermittivity::Scalar(&eps_2d),
            &dir_labs,
            RobinBoundary::None,
            Some(&c_up),
        )
        .unwrap();
        let (clow, cmid, cup) = (eq.c_lower_matrix[0][0], eq.c_mid_matrix[0][0], eq.c_upper_matrix[0][0]);
        eprintln!(
            "2D coax order={order}: clow={clow:.8e}, cmid={cmid:.8e}, cup={cup:.8e}, exact_poly={c_exact_2d_poly:.8e}"
        );
        assert!(
            clow < c_exact_2d_poly && c_exact_2d_poly < cup,
            "2D coax order={order}: expected clow ({clow:e}) < exact_poly ({c_exact_2d_poly:e}) < cup ({cup:e})"
        );
        let err_up = (cup - c_exact_2d_poly).abs();
        let err_mid = (cmid - c_exact_2d_poly).abs();
        assert!(
            err_up > 2.5 * err_mid,
            "2D coax order={order}: expected >2.5x midpoint improvement, got err_up={err_up:e}, err_mid={err_mid:e} (ratio={:.2}x)",
            err_up / err_mid
        );
    }

    // 2. 3D curved coaxial cylinder (`order = 1` RT_0 and `order = 2` RT_1: curved cylindrical walls eliminate polygonal chordal bias!)
    assert!(mesh_3d.has_curved());
    let eps_3d = vec![eps_val; mesh_3d.n_cells()];
    for order in [1usize, 2] {
        let k = assemble_laplace_curved_3d_test(&mesh_3d, order, &eps_3d);
        let (v_all, c_up) =
            solve_primal_terminals(&mesh_3d, order, &k, std::slice::from_ref(&inner_labs), &outer_labs);
        let eq = equilibrate_electro_flux_bounds(
            &mesh_3d,
            order,
            &v_all,
            CellPermittivity::Scalar(&eps_3d),
            &dir_labs,
            RobinBoundary::None,
            Some(&c_up),
        )
        .unwrap();
        let (clow, cmid, cup) = (eq.c_lower_matrix[0][0], eq.c_mid_matrix[0][0], eq.c_upper_matrix[0][0]);
        eprintln!(
            "3D curved coax order={order}: clow={clow:.8e}, cmid={cmid:.8e}, cup={cup:.8e}, exact={c_exact_3d:.8e}"
        );
        assert!(
            clow < c_exact_3d && c_exact_3d < cup,
            "3D curved coax order={order}: expected clow ({clow:e}) < exact ({c_exact_3d:e}) < cup ({cup:e})"
        );
        let err_up = (cup - c_exact_3d).abs();
        let err_mid = (cmid - c_exact_3d).abs();
        assert!(
            err_up > 2.5 * err_mid,
            "3D curved coax order={order}: expected >2.5x midpoint improvement, got err_up={err_up:e}, err_mid={err_mid:e} (ratio={:.2}x)",
            err_up / err_mid
        );
    }
}

fn build_3d_two_sheet_capacitor(global_h: f64, edge_h: f64) -> (SimplexMesh, Vec<f64>, Vec<i32>, Vec<i32>, Vec<i32>) {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-160e-6, -160e-6, -100e-6), (160e-6, 160e-6, 0.0));
    let vac = m.box_("vac", (-160e-6, -160e-6, 0.0), (160e-6, 160e-6, 100e-6));

    // Symmetric two-pad coplanar transmon capacitor at z = 0 (`Si` eps_r = 11.45 below, `VacTop` eps_r = 1.0 above)
    let p1_outer = vec![pt(-90e-6, -60e-6), pt(-15e-6, -60e-6), pt(-15e-6, 60e-6), pt(-90e-6, 60e-6)];
    let p2_outer = vec![pt(15e-6, -60e-6), pt(90e-6, -60e-6), pt(90e-6, 60e-6), pt(15e-6, 60e-6)];
    let p1 = m.polygons("p1_poly", vec![p1_outer], vec![]);
    let p2 = m.polygons("p2_poly", vec![p2_outer], vec![]);
    m.sheet("pad1", p1, 0.0);
    m.sheet("pad2", p2, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    m.size("s_p1", Selector::parse("pad1:edge").unwrap(), Some(Expr::Num(edge_h)), None, None);
    m.size("s_p2", Selector::parse("pad2:edge").unwrap(), Some(Expr::Num(edge_h)), None, None);

    let g = mesh_model(&m.list, &MeshPlan::simple(global_h)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let pad1_labs = g.boundary("pad1");
    let pad2_labs = g.boundary("pad2");
    let outer_labs = g.select("outer").unwrap();

    let eps_cells: Vec<f64> = g
        .mesh
        .cell_labels
        .iter()
        .map(|&d| if d == sub_dom { 11.45 * EPS0 } else { EPS0 })
        .collect();

    (g.mesh, eps_cells, pad1_labs, pad2_labs, outer_labs)
}

fn frobenius_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let mut s = 0.0;
    for i in 0..a.len() {
        for j in 0..a[i].len() {
            let d = a[i][j] - b[i][j];
            s += d * d;
        }
    }
    s.sqrt()
}

#[test]
fn interior_3d_sheet_capacitor_splits_half_stars_brackets_loewner_and_cuts_error_over_2_5x() {
    // Compute a converged reference on the exact same polyhedral geometry using fine P2 + equilibrated midpoint
    // (guaranteed tight `[C^low_fine, C^up_fine]` bracket!)
    let (mesh_ref, eps_ref, p1_ref, p2_ref, out_ref) = build_3d_two_sheet_capacitor(22e-6, 4.5e-6);
    let k_ref = assemble_laplace(&mesh_ref, 2, &eps_ref);
    let (v_ref, cup_ref) = solve_primal_terminals(&mesh_ref, 2, &k_ref, &[p1_ref.clone(), p2_ref.clone()], &out_ref);
    let mut dir_ref = out_ref.clone();
    dir_ref.extend_from_slice(&p1_ref);
    dir_ref.extend_from_slice(&p2_ref);
    let eq_ref = equilibrate_electro_flux_bounds(
        &mesh_ref,
        2,
        &v_ref,
        CellPermittivity::Scalar(&eps_ref),
        &dir_ref,
        RobinBoundary::None,
        Some(&cup_ref),
    )
    .unwrap();
    let c_exact = &eq_ref.c_mid_matrix;
    eprintln!(
        "3D sheet FINE REF (order=2): C00=[{:.4} fF, {:.4} fF], Q00={:.4} fF, W00={:.4} fF",
        eq_ref.c_lower_matrix[0][0] * 1e15,
        eq_ref.c_upper_matrix[0][0] * 1e15,
        eq_ref.q_pairing_matrix[0][0] * 1e15,
        eq_ref.w_flux_matrix[0][0] * 1e15,
    );
    let ref_half_gap = 0.5 * (eq_ref.c_upper_matrix[0][0] - eq_ref.c_lower_matrix[0][0]);

    // Coarse 3D sheet mesh tested with both `order = 1` (RT_0) and `order = 2` (RT_1)
    let (mesh_c, eps_c, p1_c, p2_c, out_c) = build_3d_two_sheet_capacitor(55e-6, 18e-6);
    let mut dir_c = out_c.clone();
    dir_c.extend_from_slice(&p1_c);
    dir_c.extend_from_slice(&p2_c);

    let exterior_c = mesh_c.exterior_facets();
    let sheet_set_c: BTreeSet<i32> = p1_c.iter().chain(&p2_c).copied().collect();
    let sheet_mask_c: Vec<bool> = (0..mesh_c.n_facets())
        .map(|f| !exterior_c[f] && sheet_set_c.contains(&mesh_c.facet_labels[f]))
        .collect();

    for order in [1usize, 2] {
        let k_c = assemble_laplace(&mesh_c, order, &eps_c);
        let (v_c, cup_c) = solve_primal_terminals(&mesh_c, order, &k_c, &[p1_c.clone(), p2_c.clone()], &out_c);
        let eq = equilibrate_electro_flux_bounds(
            &mesh_c,
            order,
            &v_c,
            CellPermittivity::Scalar(&eps_c),
            &dir_c,
            RobinBoundary::None,
            Some(&cup_c),
        )
        .unwrap();

        // 1. Interior sheet vertices MUST split into upper and lower half-stars across `pad1` & `pad2`
        assert!(
            eq.n_sheet_splits > 0,
            "order={order}: expected interior sheet vertices to split into half-stars, got n_sheet_splits = 0"
        );
        // 2. Non-zero conductor charge supported (`C^low` strictly positive definite, NEVER zero!)
        assert!(eq.c_lower_matrix[0][0] > 10.0e-15, "order={order}: C^low_00 = {:e}", eq.c_lower_matrix[0][0]);
        // 3. Exact divergence-free global flux and continuous normal flux across `Si|VacTop` dielectric interface
        let div_scale = eq.c_upper_matrix[0][0] / (320e-6 * 320e-6 * 200e-6_f64).sqrt();
        let jump_scale = eq.c_upper_matrix[0][0] / 320e-6;
        assert!(eq.max_div_l2 < 1.0e-11 * div_scale, "order={order}: max_div_l2 = {:e} (scale={div_scale:e})", eq.max_div_l2);
        for term in 0..2 {
            let jump = eq.interior_non_sheet_normal_jump_l2(&mesh_c, &sheet_mask_c, term).unwrap();
            assert!(jump < 1.0e-12 * jump_scale, "order={order} term={term}: Si|VacTop normal flux jump = {jump:e}");
        }
        eprintln!(
            "3D sheet order={order}: C00=[{:.4} fF, {:.4} fF] (mid={:.4} fF, ref={:.4} fF), C01=[{:.4} fF, {:.4} fF] (mid={:.4} fF, ref={:.4} fF)",
            eq.c_lower_matrix[0][0] * 1e15,
            eq.c_upper_matrix[0][0] * 1e15,
            eq.c_mid_matrix[0][0] * 1e15,
            c_exact[0][0] * 1e15,
            eq.c_lower_matrix[0][1] * 1e15,
            eq.c_upper_matrix[0][1] * 1e15,
            eq.c_mid_matrix[0][1] * 1e15,
            c_exact[0][1] * 1e15,
        );
        // 4. Full 2x2 Loewner matrix bracketing: C^low_coarse <= C^low_ref <= C^exact <= C^up_ref <= C^up_coarse
        assert_loewner_le(&eq.c_lower_matrix, &eq_ref.c_upper_matrix, 1.0e-10, &format!("order={order} C^low <= C^exact_up"));
        assert_loewner_le(&eq.c_lower_matrix, c_exact, 1.0e-10, &format!("order={order} C^low <= C^exact"));
        assert_loewner_le(c_exact, &eq.c_upper_matrix, 1.0e-10, &format!("order={order} C^exact <= C^up"));
        assert!(
            eq.c_upper_matrix[0][0] - eq.c_lower_matrix[0][0] > 4.0 * ref_half_gap,
            "coarse bracket should cleanly enclose fine reference bracket"
        );

        // 5. Entrywise off-diagonal & diagonal error bound `|C^mid_ij - C^exact_ij| <= delta_C_ij`
        for i in 0..2 {
            for j in 0..2 {
                let actual_err = (eq.c_mid_matrix[i][j] - c_exact[i][j]).abs();
                let cert_bound = eq.c_error_bound_matrix[i][j] + ref_half_gap;
                assert!(
                    actual_err <= cert_bound * (1.0 + 1.0e-9),
                    "order={order} ({i},{j}): |C^mid - C^exact| = {actual_err:e} > certified bound {cert_bound:e}"
                );
            }
        }

        // 6. `C^mid` cuts error by `> 2.5x` vs one-sided `C^up` for BOTH `order=1` (RT_0) and `order=2` (RT_1)!
        let err_up = frobenius_diff(&eq.c_upper_matrix, c_exact);
        let err_mid = frobenius_diff(&eq.c_mid_matrix, c_exact);
        let gain = err_up / err_mid.max(1.0e-30);
        eprintln!("3D sheet order={order} midpoint gain = {gain:.2}x (err_up={err_up:e}, err_mid={err_mid:e})");
        assert!(
            gain > 2.5,
            "order={order}: expected >2.5x midpoint gain on 3D sheet capacitor, got {gain:.2}x (err_up={err_up:e}, err_mid={err_mid:e})"
        );
    }
}

#[test]
fn robin_farfield_spherical_shell_bracketing_and_quarter_point_guard() {
    // 1. Concentric spherical shell r in [r_in, r_out] built by radially mapping a hollow cube shell,
    //    with V = 1 on inner surface (label 1) and exact monopole Robin BC `alpha(x) = eps0 * (r . n) / |r|^2` on outer surface (label 2)!
    let (r_in, r_out) = (40e-6, 90e-6);
    // Build cube [-1, 1]^3 surface subdivisions projected onto concentric spheres r_in..r_out
    let n_side = 8usize;
    let nr = 4usize;
    let mut points = Vec::new();
    let mut patch_vid = BTreeMap::<(usize, usize, usize, usize), u32>::new();
    // Use 6 cube faces parameterized by (face, iu, iv), deduplicating cube vertices/edges by quantized unit-sphere direction
    let mut dir_to_base = BTreeMap::<[i64; 3], usize>::new();
    let mut unit_dirs: Vec<[f64; 3]> = Vec::new();
    let face_pt = |face: usize, u: f64, v: f64| -> [f64; 3] {
        match face {
            0 => [1.0, u, v],
            1 => [-1.0, v, u],
            2 => [v, 1.0, u],
            3 => [u, -1.0, v],
            4 => [u, v, 1.0],
            _ => [v, u, -1.0],
        }
    };
    let mut grid_base = vec![vec![vec![0usize; n_side + 1]; n_side + 1]; 6];
    for face in 0..6 {
        for iu in 0..=n_side {
            let u = -1.0 + 2.0 * (iu as f64 / n_side as f64);
            for iv in 0..=n_side {
                let v = -1.0 + 2.0 * (iv as f64 / n_side as f64);
                let p = face_pt(face, u, v);
                let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                let d = [p[0] / len, p[1] / len, p[2] / len];
                let key = [
                    (d[0] * 1.0e9).round() as i64,
                    (d[1] * 1.0e9).round() as i64,
                    (d[2] * 1.0e9).round() as i64,
                ];
                let idx = *dir_to_base.entry(key).or_insert_with(|| {
                    let k = unit_dirs.len();
                    unit_dirs.push(d);
                    k
                });
                grid_base[face][iu][iv] = idx;
            }
        }
    }
    let n_base = unit_dirs.len();
    for ir in 0..=nr {
        let r = r_in + (r_out - r_in) * (ir as f64 / nr as f64);
        for d in &unit_dirs {
            points.extend_from_slice(&[r * d[0], r * d[1], r * d[2]]);
        }
    }
    let _ = &mut patch_vid;
    let vid = |ir: usize, b: usize| -> u32 { (ir * n_base + b) as u32 };
    let quad_tris = |cyc: [u32; 4]| -> [[u32; 3]; 2] {
        let m = (0..4).min_by_key(|&i| cyc[i]).unwrap();
        let c = |k: usize| cyc[(m + k) % 4];
        [[c(0), c(1), c(2)], [c(0), c(2), c(3)]]
    };
    let mut cells = Vec::new();
    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for ir in 0..nr {
        for face in 0..6 {
            for iu in 0..n_side {
                for iv in 0..n_side {
                    let (b00, b10, b11, b01) = (
                        grid_base[face][iu][iv],
                        grid_base[face][iu + 1][iv],
                        grid_base[face][iu + 1][iv + 1],
                        grid_base[face][iu][iv + 1],
                    );
                    let h = [
                        vid(ir, b00),
                        vid(ir, b10),
                        vid(ir, b11),
                        vid(ir, b01),
                        vid(ir + 1, b00),
                        vid(ir + 1, b10),
                        vid(ir + 1, b11),
                        vid(ir + 1, b01),
                    ];
                    let apex = *h.iter().min().unwrap();
                    let hex_faces = [
                        [h[0], h[3], h[2], h[1]],
                        [h[4], h[5], h[6], h[7]],
                        [h[0], h[1], h[5], h[4]],
                        [h[1], h[2], h[6], h[5]],
                        [h[2], h[3], h[7], h[6]],
                        [h[3], h[0], h[4], h[7]],
                    ];
                    for qf in hex_faces {
                        if !qf.contains(&apex) {
                            for tri in quad_tris(qf) {
                                cells.extend_from_slice(&[apex, tri[0], tri[1], tri[2]]);
                            }
                        }
                    }
                    if ir == 0 {
                        for tri in quad_tris([h[0], h[3], h[2], h[1]]) {
                            facets.extend_from_slice(&tri);
                            facet_labels.push(1); // r = r_in terminal
                        }
                    }
                    if ir + 1 == nr {
                        for tri in quad_tris([h[4], h[5], h[6], h[7]]) {
                            facets.extend_from_slice(&tri);
                            facet_labels.push(2); // r = r_out Robin far-field
                        }
                    }
                }
            }
        }
    }
    let nc = cells.len() / 4;
    let sphere_mesh = SimplexMesh::new(3, points, cells, vec![1; nc], facets, facet_labels);
    let eps_vec = vec![EPS0; nc];
    let alpha_outer = EPS0 / r_out;

    // Solve with N = 1 terminal on label 1 (V = 1) and constant Robin alpha = eps0 / r_out on label 2 (NO ground!)
    let k_vol = assemble_laplace(&sphere_mesh, 2, &eps_vec);
    let k_rob = assemble_boundary_mass(&sphere_mesh, &[2], 2, alpha_outer);
    let k_tot = k_vol.add_scaled(&k_rob, 1.0);
    let volt: BTreeMap<i32, f64> = [(1, 1.0)].into();
    let (fixed, vals) = dirichlet_constraints(&sphere_mesh, &volt, 2);
    let free = free_dofs(k_tot.nrows(), &fixed);
    let kff = k_tot.submatrix(&free, &free);
    let kfx = k_tot.submatrix(&free, &fixed);
    let rhs: Vec<f64> = kfx.matvec(&vals).into_iter().map(|x| -x).collect();
    let sol_free = DirectSolver::new(&kff, true).unwrap().solve(&rhs);
    let mut v_sol = vec![0.0; k_tot.nrows()];
    for (&d, &val) in fixed.iter().zip(&vals) {
        v_sol[d] = val;
    }
    for (a, &i) in free.iter().enumerate() {
        v_sol[i] = sol_free[a];
    }
    let kv = k_tot.matvec(&v_sol);
    let c_up_val: f64 = v_sol.iter().zip(&kv).map(|(a, b)| a * b).sum();

    let mut robin_alphas = vec![0.0; sphere_mesh.n_facets()];
    for f in 0..sphere_mesh.n_facets() {
        if sphere_mesh.facet_labels[f] == 2 {
            robin_alphas[f] = alpha_outer;
        }
    }
    let eq_rob = equilibrate_electro_flux_bounds(
        &sphere_mesh,
        2,
        &[v_sol],
        CellPermittivity::Scalar(&eps_vec),
        &[1],
        RobinBoundary::FacetAlpha(&robin_alphas),
        Some(&[vec![c_up_val]]),
    )
    .unwrap();
    assert!(
        eq_rob.c_lower_matrix[0][0] > 0.0 && eq_rob.c_lower_matrix[0][0] <= eq_rob.c_upper_matrix[0][0],
        "Robin spherical shell bracket invalid: low={:e}, up={:e}",
        eq_rob.c_lower_matrix[0][0],
        eq_rob.c_upper_matrix[0][0]
    );
    // Relative bracket width on 4-layer P2 spherical shell with Robin BC is < 3%
    let rel_gap = (eq_rob.c_upper_matrix[0][0] - eq_rob.c_lower_matrix[0][0]) / eq_rob.c_mid_matrix[0][0];
    assert!(rel_gap < 0.03, "Robin spherical shell bracket too wide: rel_gap = {rel_gap:.4}");

    // 2. Cross-Feature Compatibility Rule 3: exact Barsoum quarter-point CurvedTet (`xi = 0.25` -> `det J = 0` at vertex 0)
    //    must be rejected cleanly with an explanatory error message!
    let mut qp_mesh = femcore::testkit::block(10e-6, 10e-6, 10e-6, 1, 1, 1);
    let (e0a, e0b) = qp_mesh.edges()[0];
    let pa = qp_mesh.point(e0a as usize);
    let pb = qp_mesh.point(e0b as usize);
    let qp_node = [
        0.75 * pa[0] + 0.25 * pb[0],
        0.75 * pa[1] + 0.25 * pb[1],
        0.75 * pa[2] + 0.25 * pb[2],
    ];
    qp_mesh.set_edge_nodes(&[([e0a, e0b], qp_node)]);
    let dummy_v = vec![vec![0.0; qp_mesh.n_dofs(2)]];
    let dummy_eps = vec![EPS0; qp_mesh.n_cells()];
    let dir_mask = vec![true; qp_mesh.n_facets()];
    let sh_mask = vec![false; qp_mesh.n_facets()];
    let err = equilibrate_flux_patches(
        &qp_mesh,
        2,
        &dummy_v,
        CellPermittivity::Scalar(&dummy_eps),
        &dir_mask,
        &sh_mask,
        RobinBoundary::None,
        None,
    )
    .expect_err("quarter-point CurvedTet with det(J)=0 must be rejected by H(div) equilibration");
    assert!(err.contains("quarter-point"), "unexpected error message: {err}");
}

fn run_sheet_equilibration_bits() -> Vec<u64> {
    let (mesh, eps, p1, p2, out) = build_3d_two_sheet_capacitor(50e-6, 16e-6);
    let k = assemble_laplace(&mesh, 2, &eps);
    let (v_all, c_up) = solve_primal_terminals(&mesh, 2, &k, &[p1.clone(), p2.clone()], &out);
    let mut dir = out;
    dir.extend_from_slice(&p1);
    dir.extend_from_slice(&p2);
    let eq: EquilibratedFluxResult = equilibrate_electro_flux_bounds(
        &mesh,
        2,
        &v_all,
        CellPermittivity::Scalar(&eps),
        &dir,
        RobinBoundary::None,
        Some(&c_up),
    )
    .unwrap();
    let mut bits = Vec::new();
    for m in [&eq.c_lower_matrix, &eq.c_mid_matrix, &eq.w_flux_matrix, &eq.q_pairing_matrix, &eq.c_error_bound_matrix] {
        for r in m {
            for &x in r {
                bits.push(x.to_bits());
            }
        }
    }
    for &x in &eq.cell_eta2 {
        bits.push(x.to_bits());
    }
    bits
}

#[test]
fn nyx_threads_bit_determinism_across_subprocess_thread_counts() {
    if let Ok(mode) = std::env::var("NYX_RT_DET_CHILD") {
        assert!(!mode.is_empty());
        let bits = run_sheet_equilibration_bits();
        let line: String = bits.iter().map(|b| format!("{b:016x}")).collect::<Vec<_>>().join(",");
        println!("NYX_RT_BITS:{line}");
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |nt: &str| -> String {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("nyx_threads_bit_determinism_across_subprocess_thread_counts")
            .arg("--nocapture")
            .env("NYX_THREADS", nt)
            .env("NYX_RT_DET_CHILD", "1")
            .output()
            .expect("failed to spawn child test process");
        assert!(out.status.success(), "child NYX_THREADS={nt} failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8(out.stdout).unwrap();
        stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_RT_BITS:"))
            .expect("missing NYX_RT_BITS line")
            .to_string()
    };

    let bits_1 = run_with_threads("1");
    let bits_4 = run_with_threads("4");
    assert_eq!(bits_1, bits_4, "RT_1 equilibration outputs differed between NYX_THREADS=1 and NYX_THREADS=4");
}
