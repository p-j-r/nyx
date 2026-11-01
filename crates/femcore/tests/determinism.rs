//! The shared thread pool must not make results depend on scheduling: the same problem solved twice in one process
//! gives bit-identical numbers (work is split into fixed pieces and combined in a fixed order).
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::solve::solve;

fn spec(order: usize, study: Study) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = order;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) }];
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
    let m = SimplexMesh::block(1.0, 0.5, 0.5, 24, 12, 12);
    let s = spec(1, Study::Static);
    let a = solve(&m, &s, &mut SolveState::default()).unwrap();
    let b = solve(&m, &s, &mut SolveState::default()).unwrap();
    assert!(a.log.iter().any(|l| l.contains("conjugate gradients") || l.contains("multigrid")), "{:?}", a.log);
    assert_eq!(bits(&a), bits(&b));
}

#[test]
fn lobpcg_gives_identical_results_every_time() {
    let m = SimplexMesh::block(2.0, 0.5, 0.5, 40, 10, 10);
    let s = spec(1, Study::Eigenfrequency { n_modes: 4, f_search: 0.0 });
    let a = solve(&m, &s, &mut SolveState::default()).unwrap();
    let b = solve(&m, &s, &mut SolveState::default()).unwrap();
    assert_eq!(bits(&a), bits(&b));
}
