//! Prints the plain-language stage log of a few small solves, to see what the simulator reports for each kind of study
//! (the text of docs/HOW_IT_WORKS.md quotes these):
//!
//!   cargo run --release -p femcore --example stage_log
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::solve::solve;

fn show(title: &str, mesh: &SimplexMesh, spec: &ProblemSpec) {
    println!("=== {title}");
    match solve(mesh, spec, &mut SolveState::default()) {
        Ok(r) => {
            for (n, v) in &r.scalars {
                if let Scalar::Num(x) = v {
                    println!("  result {n} = {x:.6e}");
                }
            }
            for l in &r.log {
                println!("  {l}");
            }
        }
        Err(e) => println!("  error: {e}"),
    }
}

fn main() {
    let rect = SimplexMesh::rectangle(1.0, 0.5, 20, 10);
    let mut es = ProblemSpec::new(Physics::Electrostatics);
    es.materials.default.insert("eps_r".into(), 4.0);
    es.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![3], kind: BcKind::Ground }];
    show("electrostatics, 2D", &rect, &es);

    let mut so = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        so.materials.default.insert(k.into(), v);
    }
    so.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { labels: vec![4], kind: BcKind::Load([0.0, -1e5, 0.0]) }];
    show("solid mechanics, static, 2D", &rect, &so);
    so.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    show("solid mechanics, 4 eigenfrequencies, 2D", &rect, &so);

    let beam = SimplexMesh::block(1.0, 0.1, 0.1, 40, 8, 8);
    so.bcs.truncate(1);
    show("solid mechanics, 6 eigenfrequencies, 3D slender beam (multigrid LOBPCG is tried above 80,000 unknowns)", &beam, &so);
    let big = SimplexMesh::block(1.0, 1.0, 1.0, 14, 14, 14);
    so.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    so.solver = LinearSolver::Iterative;
    show("solid mechanics, 6 eigenfrequencies, 3D cube, multigrid LOBPCG forced", &big, &so);

    let guide = SimplexMesh::rectangle(22.86e-3, 10.16e-3, 40, 18);
    let mut em = ProblemSpec::new(Physics::Emw);
    em.materials.default.insert("eps_r".into(), 1.0);
    em.materials.default.insert("mu_r".into(), 1.0);
    em.study = Study::Modes { freq: 20e9, n_target: 0.9, evanescent: false, n_modes: 3 };
    show("electromagnetic modes of a rectangular waveguide (WR-90), 2D cross-section", &guide, &em);
}
