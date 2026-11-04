//! The multigrid-preconditioned CG path must give the same answers as the direct solver.
mod common;
use common::maxdiff;
use femcore::props::*;
use nyx_sms::*;

#[test]
fn iterative_elasticity_matches_the_direct_solution() {
    let mesh = common::block(1.0, 0.2, 0.2, 20, 4, 4);
    let mut s = Problem::new(Study::Static, 2);
    s.materials.default.insert("E".into(), 200e9);
    s.materials.default.insert("nu".into(), 0.3);
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) },
    ];
    s.solver = LinearSolver::Direct;
    let d = common::solve(&mesh, &s);
    s.solver = LinearSolver::Iterative;
    let it = common::solve(&mesh, &s);
    assert!(it.log.iter().any(|l| l.starts_with("AMG-CG")), "{:?}", it.log);
    for name in ["ux", "uy", "uz", "von Mises"] {
        let (a, b) = (&d.solutions[0].flat(name).unwrap(), &it.solutions[0].flat(name).unwrap());
        let sc = a.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(maxdiff(a, b) < 1e-6 * sc, "{name}: {} (scale {sc})", maxdiff(a, b));
    }
    eprintln!("{}", it.log[0]);
}

#[test]
fn iterative_thermal_expansion_and_prescribed_displacement() {
    // 2D plate with a prescribed displacement and a mixed constraint (roller)
    let mesh = common::rect(1.0, 0.2, 60, 12);
    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Stress;
    s.materials.default.insert("E".into(), 1e9);
    s.materials.default.insert("nu".into(), 0.25);
    s.materials.default.insert("alpha".into(), 1e-5);
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { labels: vec![4], kind: BcKind::Fixed { fixed: [true, false, false], values: [1e-3, 0.0, 0.0] } },
    ];
    s.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(350.0), t_ref: 293.15 });
    s.solver = LinearSolver::Direct;
    let d = common::solve(&mesh, &s);
    s.solver = LinearSolver::Iterative;
    let it = common::solve(&mesh, &s);
    let (a, b) = (&d.solutions[0].flat("ux").unwrap(), &it.solutions[0].flat("ux").unwrap());
    assert!(maxdiff(a, b) < 1e-6 * 1e-3, "{}", maxdiff(a, b));
}

#[test]
fn eigenfrequencies_by_lobpcg_match_the_factorisation() {
    let mesh = common::block(1.0, 0.2, 0.2, 14, 3, 3);
    let mut s = Problem::new(Study::Static, 2);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    s.solver = LinearSolver::Direct;
    let d = common::solve(&mesh, &s);
    s.solver = LinearSolver::Iterative;
    let it = common::solve(&mesh, &s);
    assert!(it.log.iter().any(|l| l.starts_with("LOBPCG")), "{:?}", it.log);
    for (a, b) in d.solutions.iter().zip(&it.solutions) {
        let (fa, fb) = (a.value.unwrap(), b.value.unwrap());
        assert!((fa - fb).abs() < 1e-6 * fa, "{fa} vs {fb}");
    }
    // the mode shapes agree too (up to sign and a rotation within a degenerate pair, so compare the norm of |u|)
    let na: f64 = d.solutions[0].flat("|u|").unwrap().iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = it.solutions[0].flat("|u|").unwrap().iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!((na - nb).abs() < 1e-5 * na);
}
