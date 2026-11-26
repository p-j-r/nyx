// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pillar (c) Verification: Formal h- and p-Convergence-Order Suite (§4.4).
//!
//! Uses `femcore::vvuq::assert_convergence_rate` across `modes`, `cavity`, and `driven` to certify:
//! - 2D `modes` `H(curl)` L^2 & curl error (P1 `[0.90, 1.15]`, P2 `[1.85, 2.20]`)
//! - 2D `modes` & 3D `cavity` symmetric eigenfrequency error `|lambda_h - lambda|` (P1 `[1.85, 2.25]`, P2 `[3.60, 4.30]`)
//! - 3D `driven` waveport S-parameter `h`- and `p`-convergence

mod common;

use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, quad_rule};
use femcore::linalg::csr::c;
use femcore::mesh::SimplexMesh;
use femcore::vvuq::assert_convergence_rate;
use nyx_ems::modes::solve_modes;
use nyx_ems::problem::{Problem, Study, WavePortBc};
use nyx_ems::C0;
use nyx_physics::DrivenPortData;
use std::f64::consts::PI;

/// Nested structured 2D triangular mesh on `[0, lx] x [0, ly]` (labels: `1: x=0, 2: y=0, 3: y=ly, 4: x=lx`).
/// Eliminates unstructured Delaunay diagonal flips across grid refinements so `p_obs` reflects true FEM order.
fn structured_rect(lx: f64, ly: f64, nx: usize, ny: usize) -> SimplexMesh {
    let mut points = Vec::with_capacity(2 * (nx + 1) * (ny + 1));
    for iy in 0..=ny {
        let y = ly * (iy as f64) / (ny as f64);
        for ix in 0..=nx {
            let x = lx * (ix as f64) / (nx as f64);
            points.push(x);
            points.push(y);
        }
    }
    let vid = |ix: usize, iy: usize| -> u32 { (iy * (nx + 1) + ix) as u32 };
    let mut cells = Vec::with_capacity(6 * nx * ny);
    let mut cell_labels = Vec::with_capacity(2 * nx * ny);
    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for iy in 0..ny {
        for ix in 0..nx {
            let (v00, v10, v01, v11) = (vid(ix, iy), vid(ix + 1, iy), vid(ix, iy + 1), vid(ix + 1, iy + 1));
            cells.extend_from_slice(&[v00, v10, v11, v00, v11, v01]);
            cell_labels.push(1);
            cell_labels.push(1);
            if ix == 0 {
                facets.extend_from_slice(&[v01, v00]);
                facet_labels.push(1);
            }
            if iy == 0 {
                facets.extend_from_slice(&[v00, v10]);
                facet_labels.push(2);
            }
            if iy + 1 == ny {
                facets.extend_from_slice(&[v11, v01]);
                facet_labels.push(3);
            }
            if ix + 1 == nx {
                facets.extend_from_slice(&[v10, v11]);
                facet_labels.push(4);
            }
        }
    }
    SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels)
}

/// Nested structured 6-tet Kuhn tetrahedral mesh on `[0, lx] x [0, ly] x [0, lz]`
/// (labels matching `testkit::block`: `1: x=0, 2: y=0, 3: z=0, 4: z=lz, 5: y=ly, 6: x=lx`).
fn structured_block(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> SimplexMesh {
    let mut points = Vec::with_capacity(3 * (nx + 1) * (ny + 1) * (nz + 1));
    for iz in 0..=nz {
        let z = lz * (iz as f64) / (nz as f64);
        for iy in 0..=ny {
            let y = ly * (iy as f64) / (ny as f64);
            for ix in 0..=nx {
                let x = lx * (ix as f64) / (nx as f64);
                points.extend_from_slice(&[x, y, z]);
            }
        }
    }
    let vid = |ix: usize, iy: usize, iz: usize| -> u32 { ((iz * (ny + 1) + iy) * (nx + 1) + ix) as u32 };
    let mut cells = Vec::with_capacity(24 * nx * ny * nz);
    let mut cell_labels = Vec::with_capacity(6 * nx * ny * nz);
    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let c000 = vid(ix, iy, iz);
                let c100 = vid(ix + 1, iy, iz);
                let c010 = vid(ix, iy + 1, iz);
                let c110 = vid(ix + 1, iy + 1, iz);
                let c001 = vid(ix, iy, iz + 1);
                let c101 = vid(ix + 1, iy, iz + 1);
                let c011 = vid(ix, iy + 1, iz + 1);
                let c111 = vid(ix + 1, iy + 1, iz + 1);
                for tet in [
                    [c000, c100, c110, c111],
                    [c000, c101, c100, c111],
                    [c000, c010, c111, c110],
                    [c000, c010, c011, c111],
                    [c000, c001, c101, c111],
                    [c000, c011, c001, c111],
                ] {
                    cells.extend_from_slice(&tet);
                    cell_labels.push(1);
                }
                if ix == 0 {
                    facets.extend_from_slice(&[c000, c010, c011, c000, c011, c001]);
                    facet_labels.extend_from_slice(&[1, 1]);
                }
                if iy == 0 {
                    facets.extend_from_slice(&[c000, c100, c101, c000, c101, c001]);
                    facet_labels.extend_from_slice(&[2, 2]);
                }
                if iz == 0 {
                    facets.extend_from_slice(&[c000, c100, c110, c000, c110, c010]);
                    facet_labels.extend_from_slice(&[3, 3]);
                }
                if iz + 1 == nz {
                    facets.extend_from_slice(&[c001, c101, c111, c001, c111, c011]);
                    facet_labels.extend_from_slice(&[4, 4]);
                }
                if iy + 1 == ny {
                    facets.extend_from_slice(&[c010, c110, c111, c010, c111, c011]);
                    facet_labels.extend_from_slice(&[5, 5]);
                }
                if ix + 1 == nx {
                    facets.extend_from_slice(&[c100, c110, c111, c100, c111, c101]);
                    facet_labels.extend_from_slice(&[6, 6]);
                }
            }
        }
    }
    SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels)
}

#[test]
fn modes_hcurl_l2_curl_and_eigenvalue_convergence_p1_p2() {
    let (a, b) = (22.86e-3_f64, 10.16e-3_f64);
    let freq = 10.0e9;
    let k0 = 2.0 * PI * freq / C0;
    let neff_exact = (1.0 - (PI / (a * k0)).powi(2)).sqrt();

    // Analytic WR-90 TE10 transverse field E_t = (0, sin(pi x / a)), curl E_t = (pi / a) cos(pi x / a)
    let ks = [4usize, 8, 16];
    let hs: Vec<f64> = ks.iter().map(|&k| a / (2 * k) as f64).collect();
    let quad = quad_rule(2, 6);

    for order in [1usize, 2] {
        let mut l2_errs = Vec::new();
        let mut curl_errs = Vec::new();
        let mut eig_errs = Vec::new();
        let nloc_e = if order == 1 { 3 } else { 8 };
        let mut w = vec![[0.0; 2]; nloc_e];
        let mut cu = vec![0.0; nloc_e];

        for &k in &ks {
            let mesh = structured_rect(a, b, 2 * k, k);
            let eps = vec![1.0; mesh.n_cells()];
            let mur = vec![1.0; mesh.n_cells()];
            let modes = solve_modes(&mesh, freq, &eps, &mur, neff_exact * 1.002, false, 1, order).unwrap();
            let neff_h = modes.neff[0].0;
            let eig_err = ((neff_h * neff_h - neff_exact * neff_exact) / (neff_exact * neff_exact)).abs();
            eig_errs.push(eig_err);

            // Compute continuous L^2 inner products to align eigenvector phase/scale with exact TE10
            let geo = geometry(&mesh);
            let cd = cell_data(&mesh, order);
            let mut dot_he = 0.0_f64;
            let mut ee_norm2 = 0.0_f64;
            for ci in 0..mesh.n_cells() {
                let g = tri_g(&geo, ci);
                let area = geo.vol[ci];
                let vs = mesh.cell(ci);
                let pts = [mesh.point(vs[0] as usize), mesh.point(vs[1] as usize), mesh.point(vs[2] as usize)];
                for (q, l) in quad.lam.iter().enumerate() {
                    let wq = quad.w[q] * area;
                    let xq = l[0] * pts[0][0] + l[1] * pts[1][0] + l[2] * pts[2][0];
                    let ex_y = (PI * xq / a).sin();
                    edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut cu);
                    let mut hy = 0.0;
                    for (s, ws) in w.iter().enumerate().take(nloc_e) {
                        hy += modes.et[0][cd.gdofs[ci][s]] * ws[1];
                    }
                    dot_he += wq * hy * ex_y;
                    ee_norm2 += wq * ex_y * ex_y;
                }
            }
            let scale = dot_he / ee_norm2;

            let (mut err_e2, mut ref_e2, mut err_c2, mut ref_c2) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
            for ci in 0..mesh.n_cells() {
                let g = tri_g(&geo, ci);
                let area = geo.vol[ci];
                let vs = mesh.cell(ci);
                let pts = [mesh.point(vs[0] as usize), mesh.point(vs[1] as usize), mesh.point(vs[2] as usize)];
                for (q, l) in quad.lam.iter().enumerate() {
                    let wq = quad.w[q] * area;
                    let xq = l[0] * pts[0][0] + l[1] * pts[1][0] + l[2] * pts[2][0];
                    let ex = [0.0, scale * (PI * xq / a).sin()];
                    let cx = scale * (PI / a) * (PI * xq / a).cos();
                    edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut cu);
                    let (mut eh, mut ch) = ([0.0; 2], 0.0);
                    for s in 0..nloc_e {
                        let coef = modes.et[0][cd.gdofs[ci][s]];
                        eh[0] += coef * w[s][0];
                        eh[1] += coef * w[s][1];
                        ch += coef * cu[s];
                    }
                    err_e2 += wq * ((eh[0] - ex[0]).powi(2) + (eh[1] - ex[1]).powi(2));
                    ref_e2 += wq * (ex[0].powi(2) + ex[1].powi(2));
                    err_c2 += wq * (ch - cx).powi(2);
                    ref_c2 += wq * cx * cx;
                }
            }
            l2_errs.push((err_e2 / ref_e2).sqrt());
            curl_errs.push((err_c2 / ref_c2).sqrt());
        }

        let field_band = if order == 1 { (0.90, 1.15) } else { (1.85, 2.20) };
        let eig_band = if order == 1 { (1.85, 2.25) } else { (3.60, 4.30) };
        eprintln!("2D modes P{order}: L2={l2_errs:?}, curl={curl_errs:?}, eig={eig_errs:?}");
        assert_convergence_rate(&hs, &l2_errs, order as f64, field_band);
        assert_convergence_rate(&hs, &curl_errs, order as f64, field_band);
        assert_convergence_rate(&hs, &eig_errs, (2 * order) as f64, eig_band);
    }
}

#[test]
fn cavity_eigenvalue_convergence_p1_p2() {
    let (a, b, d): (f64, f64, f64) = (0.04, 0.03, 0.02);
    let f_exact = 0.5 * C0 * ((1.0 / a).powi(2) + (1.0 / b).powi(2)).sqrt();
    let lam_exact = f_exact * f_exact;

    // Nested grids: P1 (4,3,2), (8,6,4), (16,12,8); P2 (3,2,2), (6,4,4), (9,6,6)
    for order in [1usize, 2] {
        let grids: Vec<(usize, usize, usize)> =
            if order == 1 { vec![(4, 3, 2), (8, 6, 4), (16, 12, 8)] } else { vec![(3, 2, 2), (6, 4, 4), (9, 6, 6)] };
        let hs: Vec<f64> = grids.iter().map(|&(nx, _, _)| a / nx as f64).collect();
        let mut errs = Vec::new();
        for &(nx, ny, nz) in &grids {
            let mesh = structured_block(a, b, d, nx, ny, nz);
            let mut p = Problem::new(Study::Cavity { n_modes: 1, f_search: 0.90 * f_exact }, order);
            p.materials.default.insert("eps_r".into(), 1.0);
            let res = common::solve(&mesh, &p).unwrap();
            let f_h = res.solutions[0].value.unwrap();
            errs.push(((f_h * f_h - lam_exact) / lam_exact).abs());
        }
        let eig_band = if order == 1 { (1.85, 2.25) } else { (3.60, 4.30) };
        eprintln!("3D cavity P{order}: hs={hs:?}, lam_errs={errs:?}");
        assert_convergence_rate(&hs, &errs, (2 * order) as f64, eig_band);
    }
}

#[test]
fn driven_waveport_sparameter_h_and_p_convergence() {
    let (a, l, b) = (22.86e-3_f64, 15.0e-3_f64, 10.16e-3_f64);
    let freq = 10.0e9;
    let k0 = 2.0 * PI * freq / C0;
    let beta = (k0 * k0 - (PI / a).powi(2)).sqrt();
    let s21_exact = c(0.0, -beta * l).exp();

    let solve_s21_err = |nx: usize, ny: usize, nz: usize, order: usize| -> f64 {
        let mesh = structured_block(a, l, b, nx, ny, nz);
        let mut p = Problem::new(
            Study::Driven {
                f_min: freq,
                f_max: freq,
                n_points: 1,
                freqs_hz: Some(vec![freq]),
                prom_tol: None,
                max_hdm_solves: 5,
                strict_prom_tol: false,
                seed_freqs_hz: Vec::new(),
                solver: "auto".into(),
                blr_tol: 1e-6,
            },
            order,
        );
        p.materials.default.insert("eps_r".into(), 1.0);
        p.waveports.push(WavePortBc { labels: vec![2], name: "p1".into(), mode: 1, n_eff_guess: None, z0_ohm: None });
        p.waveports.push(WavePortBc { labels: vec![5], name: "p2".into(), mode: 1, n_eff_guess: None, z0_ohm: None });
        let res = common::solve(&mesh, &p).unwrap();
        let data = DrivenPortData::of(&res).unwrap();
        let s21 = c(data.s_matrix[0][1][0][0], data.s_matrix[0][1][0][1]);
        (s21 - s21_exact).norm()
    };

    // h-convergence of S21 for P1 and P2
    let ks_p1 = [4usize, 8, 12];
    let hs_p1: Vec<f64> = ks_p1.iter().map(|&k| a / k as f64).collect();
    let errs_p1: Vec<f64> = ks_p1.iter().map(|&k| solve_s21_err(k, k, (k / 2).max(2), 1)).collect();
    eprintln!("Driven WR-90 S21 P1 h-conv: errs={errs_p1:?}");
    assert_convergence_rate(&hs_p1, &errs_p1, 2.0, (1.75, 2.35));

    let ks_p2 = [4usize, 6, 8];
    let hs_p2: Vec<f64> = ks_p2.iter().map(|&k| a / k as f64).collect();
    let errs_p2: Vec<f64> = ks_p2.iter().map(|&k| solve_s21_err(k, k, (k / 2).max(2), 2)).collect();
    eprintln!("Driven WR-90 S21 P2 h-conv: errs={errs_p2:?}");
    assert_convergence_rate(&hs_p2, &errs_p2, 4.0, (3.40, 4.50));

    // p-convergence on identical mesh (4, 4, 2): P2 error << P1 error
    assert!(
        errs_p2[0] < 0.05 * errs_p1[0],
        "p-refinement at (4,4,2): P2 err {:.3e} not << 0.05 * P1 err {:.3e}",
        errs_p2[0],
        errs_p1[0]
    );
}
