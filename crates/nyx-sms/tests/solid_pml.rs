// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Perfectly matched layers in solid mechanics eigenfrequencies.
mod common;
use femcore::props::*;
use nyx_physics::{Scalar, SolveResult};
use nyx_sms::*;

fn bar_spec(plane: Plane, nm: usize) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.plane = plane;
    s.study = Study::Eigenfrequency { n_modes: nm, f_search: 0.0 };
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s
}

fn freqs(r: &SolveResult) -> Vec<f64> {
    r.solutions.iter().map(|s| s.value.unwrap()).collect()
}

#[test]
fn pml_with_unit_stretch_reproduces_the_real_solver() {
    let mesh = common::rect(1.0, 0.1, 30, 3);
    let real = solve_solid(&mesh, &bar_spec(Plane::Stress, 4)).unwrap();
    let mut s = bar_spec(Plane::Stress, 4);
    s.materials
        .pml
        .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.0, order: 0.0, re: 1.0 });
    let pml = solve_solid(&mesh, &s).unwrap();
    for (a, b) in freqs(&real).iter().zip(freqs(&pml)) {
        assert!((a / b - 1.0).abs() < 1e-8, "{a} vs {b}");
    }
    assert!(pml.solutions.iter().all(|s| s.field("Re ux").is_some() && s.field("Re uz").is_some()));
}

#[test]
fn real_stretch_equals_a_longer_bar() {
    let short = common::rect(1.0, 0.1, 30, 3);
    // the same mesh stretched in x, so that the two discretizations are the same one
    let mut long = short.clone();
    long.points.chunks_mut(2).for_each(|p| p[0] *= 1.5);
    let long =
        femcore::mesh::SimplexMesh::new(2, long.points, long.cells, long.cell_labels, long.facets, long.facet_labels);
    let reference = solve_solid(&long, &bar_spec(Plane::Strain, 4)).unwrap();
    let mut s = bar_spec(Plane::Strain, 4);
    s.materials
        .pml
        .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.0, order: 0.0, re: 1.5 });
    let r = solve_solid(&short, &s).unwrap();
    for (a, b) in freqs(&reference).iter().zip(freqs(&r)) {
        assert!((a / b - 1.0).abs() < 1e-8, "{a} vs {b}");
    }
}

#[test]
fn complex_stretch_gives_the_analytic_complex_axial_frequencies() {
    // 1D axial modes of a bar clamped at x = 0: f_n = (2n - 1) c / (4 s L), c = sqrt(E / rho); with
    // s = 1 - 0.1 j the frequency is complex and Q = Re f / (2 |Im f|) = 5
    let (e, rho, l) = (200e9f64, 7850.0f64, 1.0f64);
    let mesh = common::rect(l, 0.1, 80, 2);
    let mut s = bar_spec(Plane::Stress, 3);
    s.materials.default.insert("nu".into(), 0.0);
    // roller: uz = 0 on the long sides, so only axial motion exists
    s.bcs.push(Bc { labels: vec![2, 3], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } });
    s.materials
        .pml
        .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.1, order: 0.0, re: 1.0 });
    let r = solve_solid(&mesh, &s).unwrap();
    let cs = (e / rho).sqrt();
    let q = r.scalar("Q1").unwrap();
    assert!((q / 5.0 - 1.0).abs() < 2e-3, "Q = {q}");
    for (k, sol) in r.solutions.iter().enumerate() {
        let n = (2 * (k + 1) - 1) as f64;
        let want = n * cs / (4.0 * l);
        // s = 1 - 0.1j -> 1/s = (1 + 0.1j)/1.01
        let wre = want / 1.01;
        let got = sol.value.unwrap();
        assert!((got / wre - 1.0).abs() < 2e-3, "mode {}: Re f {got} vs {wre}", k + 1);
    }
}

#[test]
fn a_pml_without_eigenfrequency_study_is_refused() {
    let mesh = common::rect(1.0, 0.1, 10, 2);
    let mut s = bar_spec(Plane::Stress, 2);
    s.study = Study::Static;
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, ..Default::default() });
    assert!(solve_solid(&mesh, &s).unwrap_err().contains("Eigenfrequency"));
}

#[test]
fn beam_on_a_pml_substrate_loses_energy_but_keeps_its_frequency() {
    let rect = |x0, z0, x1, z1, d| common::part(x0, z0, x1, z1, d, None);
    // substrate 1 x 0.5 (domain 2) with PML strips on the left/right/bottom (3, 4, 5 and the corners 6, 7),
    // a vertical beam 0.02 x 0.3 (domain 1) on top
    let parts = [
        rect(0.0, 0.0, 1.0, 0.5, 2),
        rect(0.0, 0.15, 0.2, 0.5, 3),
        rect(0.8, 0.15, 1.0, 0.5, 4),
        rect(0.2, 0.0, 0.8, 0.15, 5),
        rect(0.0, 0.0, 0.2, 0.15, 6),
        rect(0.8, 0.0, 1.0, 0.15, 7),
        rect(0.49, 0.5, 0.51, 0.8, 1),
    ];
    let mesh = common::parts_mesh(&parts, 0.03);
    // fix every boundary facet below the substrate surface (the outer edge of the layers)
    let mut outer: Vec<i32> = (0..mesh.n_facets())
        .filter(|&f| {
            let (p, q) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
            (p[1] + q[1]) / 2.0 < 0.5 - 1e-9
        })
        .map(|f| mesh.facet_labels[f])
        .collect();
    outer.sort();
    outer.dedup();
    let mut s = bar_spec(Plane::Strain, 4);
    s.bcs = vec![Bc { labels: outer, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let st = |x, z| Pml { x, y: PmlDir::Off, z, strength: 3.0, order: 2.0, re: 1.0 };
    s.materials.pml.insert(3, st(PmlDir::Minus, PmlDir::Off));
    s.materials.pml.insert(4, st(PmlDir::Plus, PmlDir::Off));
    s.materials.pml.insert(5, st(PmlDir::Off, PmlDir::Minus));
    s.materials.pml.insert(6, st(PmlDir::Minus, PmlDir::Minus));
    s.materials.pml.insert(7, st(PmlDir::Plus, PmlDir::Minus));
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    let r = solve_solid(&mesh, &s).unwrap();
    // the beam alone, clamped at its foot
    let beam = common::parts_mesh(&[rect(0.49, 0.5, 0.51, 0.8, 1)], 0.01);
    // the foot is the boundary at the bottom of the beam
    let foot: Vec<i32> = {
        let mut l: Vec<i32> = (0..beam.n_facets())
            .filter(|&f| {
                beam.point(beam.facet(f)[0] as usize)[1] < 0.5 + 1e-9
                    && beam.point(beam.facet(f)[1] as usize)[1] < 0.5 + 1e-9
            })
            .map(|f| beam.facet_labels[f])
            .collect();
        l.sort();
        l.dedup();
        l
    };
    let mut b = bar_spec(Plane::Strain, 1);
    b.bcs = vec![Bc { labels: foot, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let fb = solve_solid(&beam, &b).unwrap().solutions[0].value.unwrap();
    let f = freqs(&r);
    let near = f.iter().cloned().min_by(|a, c| (a - fb).abs().partial_cmp(&(c - fb).abs()).unwrap()).unwrap();
    assert!((near / fb - 1.0).abs() < 0.08, "beam clamped at its foot {fb} Hz, on the PML substrate {f:?}");
    let q = r
        .scalars
        .iter()
        .filter(|(n, _)| n.starts_with('Q'))
        .filter_map(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
        .fold(f64::INFINITY, f64::min);
    eprintln!("beam alone {fb} Hz; modes on substrate {f:?}; smallest Q {q}");
    assert!(q.is_finite() && q > 1.0, "expected finite quality factors, got {q}");
}

#[test]
fn low_reflecting_boundary_gives_the_analytic_damped_bar_modes() {
    // bar clamped at x = 0, absorbing end at x = L whose adjacent material has impedance ratio zeta = 2
    // (a thin end layer with 2 rho, 2 E: same wave speed). Standing waves: exp(2 i k L) = (zeta - 1)/(zeta + 1),
    // so Re f_n = n c / (2 L) and Q_n = n pi / ln((zeta + 1)/(zeta - 1))
    let (e, rho, l) = (200e9f64, 7850.0f64, 1.0f64);
    let layer = 0.01;
    let mesh = common::rect_split_x(l, 0.05, l - layer, layer / 4.0);
    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Strain;
    for (k, v) in [("E", e), ("nu", 0.0), ("rho", rho)] {
        s.materials.set(1, k, v);
    }
    for (k, v) in [("E", 2.0 * e), ("nu", 0.0), ("rho", 2.0 * rho)] {
        s.materials.set(2, k, v);
    }
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        // rollers along the long sides: uz = 0, so the motion is uniaxial
        Bc { labels: vec![2, 3], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { labels: vec![4], kind: BcKind::LowReflecting },
    ];
    let c = (e / rho).sqrt();
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 1.05 * c / (2.0 * l) };
    let r = solve_solid(&mesh, &s).unwrap();
    let zeta: f64 = 2.0;
    let qn = |n: f64| n * std::f64::consts::PI / ((zeta + 1.0) / (zeta - 1.0)).ln();
    let first = &r.solutions[0];
    let f1 = first.value.unwrap();
    assert!((f1 / (c / (2.0 * l)) - 1.0).abs() < 0.03, "Re f1 = {f1} vs {}", c / (2.0 * l));
    let q1 = r.scalar("Q1").unwrap();
    assert!((q1 / qn(1.0) - 1.0).abs() < 0.06, "Q1 = {q1} vs {}", qn(1.0));
    let f2 = r.solutions[1].value.unwrap();
    assert!((f2 / (2.0 * c / (2.0 * l)) - 1.0).abs() < 0.03, "Re f2 = {f2}");
    let q2 = r.scalar("Q2").unwrap();
    assert!((q2 / qn(2.0) - 1.0).abs() < 0.06, "Q2 = {q2} vs {}", qn(2.0));
}
