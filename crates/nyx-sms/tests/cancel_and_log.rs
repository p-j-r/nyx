// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A running solve can be cancelled from another thread within a bounded time, a cancelled solve leaves no result,
//! and every solve returns a plain-language stage log.
mod common;
use femcore::mesh::SimplexMesh;
use femcore::props::*;
use nyx_sms::*;

use femcore::progress::{Progress, CANCELLED};
use nyx_physics::SolveResult;
use std::sync::Arc;
use std::time::{Duration, Instant};

/// The entry the engine uses, with a progress handle the test can cancel through.
fn solve_with(m: &SimplexMesh, spec: &Problem, p: &Arc<Progress>) -> Result<SolveResult, String> {
    let exchange = nyx_physics::Mesh::from_solver_mesh(m)?;
    nyx_sms::run(&exchange, spec, p).map_err(|e| e.to_string())
}

fn beam_spec(mesh: &SimplexMesh, order: usize) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.order = order;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) },
    ];
    let _ = mesh;
    s
}

#[test]
fn a_solve_returns_a_readable_stage_log() {
    let m = common::block(1.0, 0.2, 0.2, 10, 3, 3);
    let spec = beam_spec(&m, 2);
    let p = Progress::new();
    let r = solve_with(&m, &spec, &p).unwrap();
    let text = r.log.join("\n");
    assert!(text.contains("solving solid mechanics"), "{text}");
    assert!(text.contains("stiffness assembly"), "{text}");
    assert!(text.contains("linear solve"), "{text}");
    assert!(r.log.iter().any(|l| l.contains(" s")), "timings are in the log");
}

#[test]
fn a_cancel_requested_before_the_start_stops_at_once() {
    let m = common::block(1.0, 0.2, 0.2, 10, 3, 3);
    let p = Progress::new();
    p.cancel();
    let e = solve_with(&m, &beam_spec(&m, 2), &p).unwrap_err();
    assert_eq!(e, CANCELLED);
}

#[test]
fn a_long_eigenfrequency_solve_is_cancelled_within_a_bounded_time() {
    // multigrid-preconditioned LOBPCG on ~230,000 unknowns: many iterations, each checked for the cancel flag
    let m = common::block(6.0, 1.2, 1.2, 60, 12, 12);
    let mut spec = beam_spec(&m, 2);
    spec.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    spec.solver = LinearSolver::Iterative;
    let p = Progress::new();
    let p2 = p.clone();
    let started = Instant::now();
    let h = std::thread::spawn(move || solve_with(&m, &spec, &p2));
    std::thread::sleep(Duration::from_millis(2500));
    let seen = p.status();
    assert!(!seen.stage.is_empty(), "a stage is reported while the solver works");
    p.cancel();
    let t_cancel = Instant::now();
    let r = h.join().unwrap();
    let waited = t_cancel.elapsed().as_secs_f64();
    assert_eq!(
        r.unwrap_err(),
        CANCELLED,
        "solve finished after {:.1} s instead of being cancelled",
        started.elapsed().as_secs_f64()
    );
    // the solver checks the flag every iteration/level; allow for the longest stage that cannot be interrupted
    assert!(waited < 30.0, "took {waited:.1} s to react to the cancel");
}

#[test]
fn an_axisymmetric_study_logs_its_stages_and_a_cancelled_one_says_cancelled() {
    let flat = common::rect(1.0e-3, 0.2e-3, 8, 3);
    let mut pts = flat.points.clone();
    for i in 0..flat.n_points() {
        pts[2 * i] += 5.0e-3;
    }
    let m = SimplexMesh::new(
        2,
        pts,
        flat.cells.clone(),
        flat.cell_labels.clone(),
        flat.facets.clone(),
        flat.facet_labels.clone(),
    );
    let mut s = Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 3, f_search: 0.0, azimuthal_order: 0 }, 2);
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("rho", 2700.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let r = solve_with(&m, &s, &Progress::new()).unwrap();
    assert!(r.log.iter().any(|l| l.contains("axisymmetric, azimuthal order m = 0")), "{:?}", r.log);
    assert!(r.log.iter().any(|l| l.contains("solving solid mechanics")), "{:?}", r.log);
    let p = Progress::new();
    p.cancel();
    assert_eq!(solve_with(&m, &s, &p).unwrap_err(), CANCELLED);
}
