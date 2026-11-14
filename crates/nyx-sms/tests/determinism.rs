//! The shared thread pool must not make results depend on scheduling: the same problem solved twice in one process
//! gives bit-identical numbers (work is split into fixed pieces and combined in a fixed order).
mod common;
use femcore::props::*;
use nyx_sms::*;

use nyx_physics::{Scalar, SolveResult};

fn spec(order: usize, study: Study) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.order = order;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) },
    ];
    s.study = study;
    s.solver = LinearSolver::Iterative;
    s
}

fn bits(r: &SolveResult) -> Vec<u64> {
    let mut v = vec![];
    for sol in &r.solutions {
        for f in &sol.fields {
            v.extend(f.values.iter().map(|x| x.to_bits()));
        }
    }
    v.extend(r.scalars.iter().filter_map(|(_, s)| if let Scalar::Num(x) = s { Some(x.to_bits()) } else { None }));
    v
}

#[test]
fn multigrid_cg_gives_identical_results_every_time() {
    let m = common::block(1.0, 0.5, 0.5, 24, 12, 12);
    let s = spec(1, Study::Static);
    let a = common::solve(&m, &s);
    let b = common::solve(&m, &s);
    assert!(a.log.iter().any(|l| l.contains("conjugate gradients") || l.contains("multigrid")), "{:?}", a.log);
    assert_eq!(bits(&a), bits(&b));
}

#[test]
fn lobpcg_gives_identical_results_every_time() {
    let m = common::block(2.0, 0.5, 0.5, 40, 10, 10);
    let s = spec(1, Study::Eigenfrequency { n_modes: 4, f_search: 0.0 });
    let a = common::solve(&m, &s);
    let b = common::solve(&m, &s);
    assert_eq!(bits(&a), bits(&b));
}

#[test]
fn the_axisymmetric_study_gives_identical_results_every_time() {
    // an absorbing layer on a ring section (complex shift-invert), solved twice in one process
    let flat = common::rect(2.0e-3, 0.4e-3, 20, 4);
    let mut pts = flat.points.clone();
    for i in 0..flat.n_points() {
        pts[2 * i] += 5.0e-3;
    }
    let m = femcore::mesh::SimplexMesh::new(
        2,
        pts,
        flat.cells.clone(),
        flat.cell_labels.clone(),
        flat.facets.clone(),
        flat.facet_labels.clone(),
    );
    let mut s = Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search: 0.0, azimuthal_order: 1 }, 2);
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("rho", 2700.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 3.0, order: 2.0, ..Pml::default() });
    let a = common::solve(&m, &s);
    let b = common::solve(&m, &s);
    assert_eq!(bits(&a), bits(&b));
    let eig = |r: &SolveResult| {
        r.solutions.iter().map(|x| x.eigenvalue.map(|e| (e[0].to_bits(), e[1].to_bits()))).collect::<Vec<_>>()
    };
    assert_eq!(eig(&a), eig(&b));
}
