// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Floquet-periodic boundary conditions combined with perfectly matched layers in `nyx-sms`,
//! and inactive domains without mass density in the PML kernel.

mod common;
use femcore::material::{Anisotropy, Stiffness};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir, PmlStabilization};
use nyx_physics::testkit::cell_fields_reproduce_the_plot_fields;
use nyx_physics::SolveResult;
use nyx_sms::*;

fn freqs(r: &SolveResult) -> Vec<f64> {
    r.solutions.iter().map(|s| s.value.unwrap()).collect()
}

fn r_plane_sapphire() -> Anisotropy {
    Anisotropy {
        stiffness: Some(Stiffness::Trigonal {
            c11: 497.35e9,
            c12: 163.97e9,
            c13: 112.20e9,
            c14: -23.58e9,
            c33: 499.11e9,
            c44: 147.39e9,
        }),
        euler_deg: [0.0, 57.6, 0.0],
        ..Default::default()
    }
}

#[test]
fn zero_strength_pml_with_floquet_phase_reproduces_lossless_floquet_2d_and_3d() {
    // 2D periodic cell along x with phase = 0.73
    let mesh2 = common::rect(1.0e-3, 0.3e-3, 8, 3);
    let mut floquet2 = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 1.0e6 }, 2);
    floquet2.materials.aniso.insert(1, r_plane_sapphire());
    floquet2.materials.default.insert("rho".into(), 3980.0);
    floquet2.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 0.73 }];
    let r_floquet2 = common::solve(&mesh2, &floquet2);

    let mut pml2 = floquet2.clone();
    pml2.materials.pml.insert(
        1,
        Pml {
            x: PmlDir::Plus,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 0.0,
            order: 2.0,
            re: 1.0,
            clamp_outer: false,
            stabilization: PmlStabilization::Off,
            ..Default::default()
        },
    );
    let r_pml2 = common::solve(&mesh2, &pml2);
    cell_fields_reproduce_the_plot_fields(&mesh2, &r_pml2);
    assert_eq!(r_floquet2.solutions.len(), r_pml2.solutions.len());
    for (fa, fb) in freqs(&r_floquet2).iter().zip(freqs(&r_pml2)) {
        let rel = (fa - fb).abs() / fa.abs();
        assert!(rel < 1e-9, "2D zero-strength PML vs Floquet: {fb} vs {fa} (rel {rel:e})");
    }

    // 3D periodic cell along y (ymin = 2 -> ymax = 5) with bottom z = 0 (3) fixed and phase = 0.73
    let mesh3 = common::block(1.0e-6, 0.2e-6, 0.4e-6, 3, 1, 2);
    let mut floquet3 = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 1);
    floquet3.materials.aniso.insert(1, r_plane_sapphire());
    floquet3.materials.default.insert("rho".into(), 3980.0);
    floquet3.bcs.push(Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } });
    floquet3.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase: 0.73 }];
    let r_floquet3 = common::solve(&mesh3, &floquet3);

    let mut pml3 = floquet3.clone();
    pml3.materials.pml.insert(
        1,
        Pml {
            x: PmlDir::Plus,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 0.0,
            order: 2.0,
            re: 1.0,
            clamp_outer: false,
            stabilization: PmlStabilization::Off,
            ..Default::default()
        },
    );
    let r_pml3 = common::solve(&mesh3, &pml3);
    cell_fields_reproduce_the_plot_fields(&mesh3, &r_pml3);
    assert_eq!(r_floquet3.solutions.len(), r_pml3.solutions.len());
    for (fa, fb) in freqs(&r_floquet3).iter().zip(freqs(&r_pml3)) {
        let rel = (fa - fb).abs() / fa.abs();
        assert!(rel < 1e-9, "3D zero-strength PML vs Floquet: {fb} vs {fa} (rel {rel:e})");
    }
}

#[test]
fn swept_3d_si_on_sapphire_with_lateral_pml_respects_floquet_phase() {
    // 2D x-z cross-section: R-plane sapphire substrate (domain 2) with XMin/XMax PML strips (domains 3, 4)
    // and a silicon ridge core on top (domain 1), swept 1 layer along y.
    let parts = [
        common::part(-2.2e-6, -1.0e-6, 2.2e-6, 0.0, 2, Some(0.35e-6)),
        common::part(-2.2e-6, -1.0e-6, -1.6e-6, 0.0, 3, Some(0.35e-6)),
        common::part(1.6e-6, -1.0e-6, 2.2e-6, 0.0, 4, Some(0.35e-6)),
        common::part(-0.35e-6, 0.0, 0.35e-6, 0.25e-6, 1, Some(0.18e-6)),
    ];
    let m2 = femcore::meshing::geomdrive::sort_nodes(&common::parts_mesh(&parts, 0.35e-6));
    let nb = m2.facet_labels.iter().copied().max().unwrap();
    let (ymin, ymax) = (nb + 1, nb + 2);
    let mesh = extrude_y(&m2, &[Layer { thickness: 0.05e-6, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();

    // Fix the bottom and outer lateral walls of the sapphire substrate
    let mut fixed_labels: Vec<i32> = (0..mesh.n_facets())
        .filter(|&f| {
            let lab = mesh.facet_labels[f];
            if lab == ymin || lab == ymax {
                return false;
            }
            let fc = mesh.facet(f);
            let z_avg = fc.iter().map(|&v| mesh.point(v as usize)[2]).sum::<f64>() / fc.len() as f64;
            let x_avg = fc.iter().map(|&v| mesh.point(v as usize)[0]).sum::<f64>() / fc.len() as f64;
            z_avg < -1.0e-6 + 1e-9 || x_avg.abs() > 2.2e-6 - 1e-9
        })
        .map(|f| mesh.facet_labels[f])
        .collect();
    fixed_labels.sort_unstable();
    fixed_labels.dedup();

    let make_problem = |phase: f64, with_pml: bool| {
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 1);
        // Domain 1: Silicon waveguide core
        for (k, v) in [("E", 169e9), ("nu", 0.22), ("rho", 2329.0)] {
            p.materials.set(1, k, v);
        }
        // Domains 2, 3, 4: R-plane Sapphire substrate + XMin / XMax PMLs
        for dom in [2, 3, 4] {
            p.materials.aniso.insert(dom, r_plane_sapphire());
            p.materials.set(dom, "rho", 3980.0);
        }
        if with_pml {
            p.materials.pml.insert(
                3,
                Pml {
                    x: PmlDir::Minus,
                    y: PmlDir::Off,
                    z: PmlDir::Off,
                    strength: 1.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
            p.materials.pml.insert(
                4,
                Pml {
                    x: PmlDir::Plus,
                    y: PmlDir::Off,
                    z: PmlDir::Off,
                    strength: 1.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
        }
        p.bcs.push(Bc { labels: fixed_labels.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } });
        p.periodic = vec![Periodic { source: vec![ymin], dest: vec![ymax], phase }];
        p
    };

    let r_pml_phase0 = common::solve(&mesh, &make_problem(0.0, true));
    let r_pml_bloch = common::solve(&mesh, &make_problem(0.281619, true));
    let r_nopml_bloch = common::solve(&mesh, &make_problem(0.281619, false));

    cell_fields_reproduce_the_plot_fields(&mesh, &r_pml_phase0);
    cell_fields_reproduce_the_plot_fields(&mesh, &r_pml_bloch);
    assert!(
        r_pml_bloch.log.iter().any(|l| l.contains("Floquet-periodic pair(s) with PML")),
        "expected Floquet+PML log line, got {:?}",
        r_pml_bloch.log
    );

    let f_phase0 = freqs(&r_pml_phase0);
    let f_pml = freqs(&r_pml_bloch);
    let f_nopml = freqs(&r_nopml_bloch);

    // Phase = 0.281619 vs phase = 0.0 must produce distinct eigenfrequencies (> 1e-2 relative diff)
    for (k, (&f0, &fp)) in f_phase0.iter().zip(&f_pml).enumerate() {
        let rel_shift = (fp - f0).abs() / fp.abs();
        assert!(
            rel_shift > 1e-2,
            "mode {} frequency did not shift with Floquet phase: phase=0 gives {f0}, phase=0.281619 gives {fp} (rel {rel_shift:e})",
            k + 1
        );
    }

    // Bound Si-core modes match the no-PML Floquet solve at phase = 0.281619 within < 1e-4 relative
    for (k, (&fp, &fnop)) in f_pml.iter().zip(&f_nopml).take(2).enumerate() {
        let rel_diff = (fp - fnop).abs() / fnop.abs();
        assert!(
            rel_diff < 1e-4,
            "bound mode {} with XMin/XMax PML ({fp}) should match lossless Floquet ({fnop}) within 1e-4, got {rel_diff:e}",
            k + 1
        );
    }
}

#[test]
fn inactive_domain_without_rho_succeeds_in_solve_solid_pml() {
    // Domain 1: active bar with PML; Domain 2: inactive air cladding without 'rho' (or 'E', 'nu')
    let mesh = common::rect_split_x(1.0, 0.2, 0.8, 0.1);
    let mut p = Problem::new(Study::Eigenfrequency { n_modes: 3, f_search: 0.0 }, 1);
    for (k, v) in [("E", 169e9), ("nu", 0.22), ("rho", 2329.0)] {
        p.materials.set(1, k, v);
    }
    p.inactive = vec![2];
    p.materials.pml.insert(
        1,
        Pml {
            x: PmlDir::Plus,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 0.2,
            order: 2.0,
            re: 1.0,
            ..Default::default()
        },
    );
    p.bcs.push(Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } });

    let r = common::solve(&mesh, &p);
    cell_fields_reproduce_the_plot_fields(&mesh, &r);
    assert!(r.solutions.len() >= 3);

    // An active domain without 'rho' must still report the exact error
    let mut bad = p.clone();
    bad.inactive.clear();
    for (k, v) in [("E", 70e9), ("nu", 0.2)] {
        bad.materials.set(2, k, v);
    }
    let err = solve_solid(&mesh, &bad).unwrap_err();
    assert_eq!(err, "domain 2 has no value for 'rho'");
}
