// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Rough-spot closures: R6d (quadrature-point axisymmetric PML `s_theta = r_tilde / r` via `AxiMaterial::with_pml`
//! and convergence under refinement) and R32 (explicit curved-mesh errors in non-cavity simplex/port kernels).

use femcore::linalg::csr::{c, C};
use femcore::mesh::SimplexMesh;
use femcore::optical::{EpsC, MuC};
use femcore::props::{Pml, PmlCoord, PmlDir, PmlProfile};
use femcore::vvuq::assert_convergence_rate;
use nyx_ems::axisym::{assemble_axisym, solve_axisymmetric, AxiMaterial};
use nyx_ems::cavity::Ned3;
use nyx_ems::modes::{solve_modes, solve_waveguide};
use nyx_ems::modes_complex::{solve_modes_complex, solve_waveguide_complex};
use nyx_ems::optical::solve_optical;
use nyx_ems::ports::assemble_waveports_3d;
use nyx_ems::problem::{Problem, Study, WavePortBc};
use nyx_physics::{Discretisation, Mesh};
use std::collections::BTreeMap;

/// Builds a two-domain `(r, z)` axisymmetric meridian rectangle `r in [0, r_phys + d_pml]`, `z in [0, h]`,
/// splitting `r = r_phys` cleanly along vertical grid lines: domain 1 is `[0, r_phys]`, domain 2 is `[r_phys, r_phys + d_pml]`.
fn axisym_split_r_mesh(r_phys: f64, d_pml: f64, h: f64, nr_phys: usize, nr_pml: usize, nz: usize) -> Discretisation {
    let nr = nr_phys + nr_pml;
    let nx1 = nr + 1;
    let ny1 = nz + 1;
    let mut pts = Vec::with_capacity(nx1 * ny1 * 2);
    for iy in 0..ny1 {
        let z = h * (iy as f64) / (nz as f64);
        for ix in 0..nx1 {
            let r = if ix <= nr_phys {
                r_phys * (ix as f64) / (nr_phys as f64)
            } else {
                r_phys + d_pml * ((ix - nr_phys) as f64) / (nr_pml as f64)
            };
            pts.push(r);
            pts.push(z);
        }
    }
    let idx = |ix: usize, iy: usize| -> u32 { (iy * nx1 + ix) as u32 };
    let mut cells = Vec::new();
    let mut cell_labels = Vec::new();
    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for iy in 0..nz {
        for ix in 0..nr {
            let dom = if ix < nr_phys { 1 } else { 2 };
            let (a, b, c, d) = (idx(ix, iy), idx(ix + 1, iy), idx(ix + 1, iy + 1), idx(ix, iy + 1));
            cells.extend_from_slice(&[a, b, c, a, c, d]);
            cell_labels.extend_from_slice(&[dom, dom]);
            if iy == 0 {
                facets.extend_from_slice(&[a, b]);
                facet_labels.push(2);
            }
            if iy + 1 == nz {
                facets.extend_from_slice(&[c, d]);
                facet_labels.push(3);
            }
            if ix == 0 {
                facets.extend_from_slice(&[d, a]);
                facet_labels.push(1);
            }
            if ix + 1 == nr {
                facets.extend_from_slice(&[b, c]);
                facet_labels.push(4);
            }
        }
    }
    SimplexMesh::new(2, pts, cells, cell_labels, facets, facet_labels)
}

#[test]
fn r6d_axisym_with_pml_real_radial_stretch_matches_wider_cylinder() {
    // Physical cylinder R = 0.08 m, H = 0.10 m; stretched outer ring R_in = 0.08 .. 0.10 m with real factor kappa = 3.0
    // so optical radius becomes 0.08 + 3.0 * 0.02 = 0.14 m (m = 0 TM010: f_0 = c0 * p_01 / (2 pi * 0.14)).
    let r_phys = 0.08;
    let d_pml = 0.02;
    let kappa = 3.0;
    let h = 0.10;
    let r_eff = r_phys + kappa * d_pml;
    let p01 = 2.404_825_557_695_773;
    let f_exact = 299_792_458.0 * p01 / (2.0 * std::f64::consts::PI * r_eff);

    let mesh = axisym_split_r_mesh(r_phys, d_pml, h, 8, 4, 4);
    let mut pml = BTreeMap::new();
    pml.insert(
        2,
        Pml {
            coord: PmlCoord::Cartesian,
            x: PmlDir::Plus,
            strength: 0.0, // purely real stretch -> real coordinate scaling kappa
            order: 0.0,
            profile: PmlProfile::KappaPolynomial { kappa_max: kappa },
            ..Pml::default()
        },
    );
    let mut idx = BTreeMap::new();
    idx.insert(1, 1.0);
    idx.insert(2, 1.0);
    let mat = AxiMaterial::from_indices(&mesh, &idx).with_pml(&mesh, &pml);
    assert!(mat.pml.is_some(), "R6d: with_pml must attach EmPmlTable for quadrature evaluation of s_theta");

    let modes = nyx_ems::axisym::solve_axisym_modes(&mesh, &mat, 0.0, f_exact * 0.98, 1, 2).unwrap();
    let f_num = modes.freqs[0];
    let rel = (f_num / f_exact - 1.0).abs();
    assert!(rel < 5e-4, "R6d real radial stretch TM010: got {f_num:.6e} vs {f_exact:.6e} (rel {rel:.2e})");
}

#[test]
fn r6d_axisym_radial_pml_quadrature_stheta_converges() {
    // Open radiating dielectric/empty disk in radial cylindrical PML: as mesh refines (nr_pml = 2, 4, 8),
    // complex frequency f_re and damping Im(f) converge monotonically at optimal order without centroid O(h/r) stagnation.
    let r_phys = 0.08;
    let d_pml = 0.04;
    let h = 0.06;
    let p01 = 2.404_825_557_695_773;
    let f_shift = 299_792_458.0 * p01 / (2.0 * std::f64::consts::PI * (r_phys + d_pml));

    let solve_at = |n: usize| -> C {
        let mesh = axisym_split_r_mesh(r_phys, d_pml, h, 2 * n, n, n);
        let mut p = Problem::new(Study::AxisymmetricModes { f_search: f_shift, azimuthal_order: 1, n_modes: 1 }, 2);
        p.materials.set(1, "eps_r", 4.0);
        p.materials.set(2, "eps_r", 1.0);
        p.materials.pml.insert(
            2,
            Pml { coord: PmlCoord::Cartesian, x: PmlDir::Plus, strength: 4.0, order: 2.0, ..Pml::default() },
        );
        let res = solve_axisymmetric(&mesh, &p).unwrap();
        let ev = res.solutions[0].eigenvalue.expect("lossy PML mode must carry complex frequency");
        c(ev[0], ev[1])
    };

    let f_coarse = solve_at(3);
    let f_mid = solve_at(6);
    let f_fine = solve_at(12);
    let f_ref = solve_at(20);
    assert!(f_fine.im.abs() > 1e-6 * f_fine.re, "PML must produce radiative damping, got {f_fine:?}");
    let errs: Vec<f64> = [f_coarse, f_mid, f_fine].iter().map(|z| (z - f_ref).norm() / f_ref.norm()).collect();
    assert!(
        errs[1] < errs[0] && errs[2] < errs[1],
        "R6d axisymmetric PML errors must decrease monotonically: {errs:?}"
    );
    let hs = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 12.0];
    let rep = assert_convergence_rate(&hs, &errs, 3.0, (2.0, 5.0));
    assert!(rep.asymptotic_slope >= 2.0, "expected order-2 convergence rate >= 2.0, got {}", rep.asymptotic_slope);
}

#[test]
fn r32_curved_meshes_are_refused_by_unsupported_non_cavity_kernels() {
    // 2D curved simplex mesh: simplex waveguide (`modes`, `modes_complex`, `optical`) and simplex `axisym` refuse cleanly.
    let mut m2 = femcore::testkit::rect(0.03, 0.015, 4, 2);
    let e0 = m2.edges()[0];
    let p0 = m2.point(e0.0 as usize);
    let p1 = m2.point(e0.1 as usize);
    m2.set_edge_nodes(&[([e0.0, e0.1], [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]) + 1e-4, 0.0])]);
    assert!(m2.has_curved());

    let wg_prob = Problem::new(Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 1 }, 2);
    let err = solve_waveguide(&m2, &wg_prob).unwrap_err();
    assert!(err.contains("curved"), "solve_waveguide error: {err}");
    let err = solve_modes(&m2, 12e9, &vec![1.0; m2.n_cells()], &vec![1.0; m2.n_cells()], 1.0, false, 1, 2)
        .err()
        .expect("solve_modes must reject curved simplex mesh");
    assert!(err.contains("curved"), "solve_modes error: {err}");
    let err = solve_waveguide_complex(&m2, &wg_prob).unwrap_err();
    assert!(err.contains("curved"), "solve_waveguide_complex error: {err}");
    let eps_c = vec![
        EpsC {
            xx: c(1.0, 0.0),
            xy: c(0.0, 0.0),
            yy: c(1.0, 0.0),
            zz: c(1.0, 0.0),
            xz: c(0.0, 0.0),
            yz: c(0.0, 0.0)
        };
        m2.n_cells()
    ];
    let mu_c = vec![MuC { x: c(1.0, 0.0), y: c(1.0, 0.0), z: c(1.0, 0.0) }; m2.n_cells()];
    let err = solve_modes_complex(&m2, 12e9, &eps_c, &mu_c, c(1.0, 0.0), 1, 2)
        .err()
        .expect("solve_modes_complex must reject curved simplex mesh");
    assert!(err.contains("curved"), "solve_modes_complex error: {err}");
    let err = solve_optical(&m2, &wg_prob).expect_err("solve_optical must reject curved simplex mesh");
    assert!(err.contains("curved"), "solve_optical error: {err}");

    let axi_prob = Problem::new(Study::AxisymmetricModes { f_search: 2e9, azimuthal_order: 0, n_modes: 1 }, 2);
    let mat = AxiMaterial::from_indices(&m2, &BTreeMap::from([(1, 1.0)]));
    let err = assemble_axisym(&m2, &mat, 0.0, 2).err().expect("assemble_axisym must reject curved simplex mesh");
    assert!(err.contains("curved"), "assemble_axisym error: {err}");
    let err = solve_axisymmetric(&m2, &axi_prob).unwrap_err();
    assert!(err.contains("curved"), "solve_axisymmetric error: {err}");

    // 3D curved waveport: `assemble_waveports_3d` & `assemble_waveports_3d_mixed` must explicitly refuse curved waveports.
    let mut m3 = femcore::testkit::block(0.03, 0.02, 0.015, 2, 2, 2);
    // Curve an edge lying on boundary label 2 (`y = 0` waveport face).
    let f_y0 = (0..m3.n_facets()).find(|&f| m3.facet_labels[f] == 2).unwrap();
    let fv = m3.facet(f_y0);
    let pa = m3.point(fv[0] as usize);
    let pb = m3.point(fv[1] as usize);
    m3.set_edge_nodes(&[(
        [fv[0], fv[1]],
        [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2]) + 1e-4],
    )]);
    let vkey: Vec<u32> = (0..m3.n_points() as u32).collect();
    let ned = Ned3::new(&m3, 1, vkey);
    let wp = vec![WavePortBc { name: "wp1".into(), labels: vec![2], mode: 1, n_eff_guess: None, z0_ohm: None }];
    let eps3 = vec![[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]; m3.n_cells()];
    let mur3 = vec![1.0; m3.n_cells()];
    let err = assemble_waveports_3d(&m3, &ned, &[], &[], &[], &wp, &eps3, &mur3, 8e9, 12e9).unwrap_err();
    assert!(err.contains("curved"), "assemble_waveports_3d error: {err}");

    // Also test mixed-path waveport refusal on curved 3D mesh
    let ex = Mesh::from_solver_mesh(&m3).unwrap();
    let mixed = ex.discretise_mixed().unwrap();
    let mut driven_p = Problem::new(
        Study::Driven {
            f_min: 8e9,
            f_max: 12e9,
            n_points: 5,
            freqs_hz: None,
            prom_tol: Some(1e-4),
            max_hdm_solves: 10,
            strict_prom_tol: false,
            seed_freqs_hz: Vec::new(),
            solver: "direct".into(),
            blr_tol: 0.0,
        },
        1,
    );
    driven_p.waveports = wp;
    let err = nyx_ems::run(&ex, &driven_p, &femcore::progress::Progress::new()).unwrap_err().to_string();
    assert!(err.contains("curved"), "driven mixed waveport error: {err}");
    let _ = mixed;
}
