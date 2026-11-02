//! Timings on larger generated models: cargo run --release -p femcore --example bench
mod benchmesh;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::solve::solve;
use std::time::Instant;

fn fixed() -> BcKind {
    BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }
}

fn run(name: &str, mesh: &SimplexMesh, spec: &ProblemSpec) {
    let t = Instant::now();
    match solve(mesh, spec, &mut SolveState::default()) {
        Ok(r) => println!("{name:<44} {:>8} dofs  {:>7.2} s {}", r.dofs, t.elapsed().as_secs_f64(), r.log.first().cloned().unwrap_or_default()),
        Err(e) => println!("{name:<44} FAILED: {e}"),
    }
}

fn main() {
    let m2 = benchmesh::rect(1.0, 0.5, 300, 150);
    println!("2D mesh: {} cells", m2.n_cells());
    let mut es = ProblemSpec::new(Physics::Electrostatics);
    es.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    run("electrostatics 2D", &m2, &es);
    let mut so = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        so.materials.default.insert(k.into(), v);
    }
    so.bcs = vec![Bc { on: None, labels: vec![1], kind: fixed() }, Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, -1e5, 0.0]) }];
    run("solid static 2D", &m2, &so);
    so.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    run("solid 6 eigenfrequencies 2D", &m2, &so);

    let m3 = benchmesh::block(1.0, 0.1, 0.1, 40, 8, 8);
    println!("3D mesh: {} cells", m3.n_cells());
    let mut s3 = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s3.materials.default.insert(k.into(), v);
    }
    s3.bcs = vec![Bc { on: None, labels: vec![1], kind: fixed() }, Bc { on: None, labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) }];
    run("solid static 3D", &m3, &s3);
    s3.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    run("solid 6 eigenfrequencies 3D", &m3, &s3);

    for (nx, ny, nz) in [(40usize, 8usize, 8usize)] {
        let big = benchmesh::block(1.0, 0.1, 0.1, nx, ny, nz);
        let mut sb = s3.clone();
        sb.study = Study::Static;
        for (name, solver) in [("direct", LinearSolver::Direct), ("AMG-CG", LinearSolver::Iterative)] {
            sb.solver = solver;
            if nx == 80 && matches!(solver, LinearSolver::Direct) {
                continue;
            }
            run(&format!("solid static 3D {nx}x{ny}x{nz} {name}"), &big, &sb);
        }
    }
    let mw = benchmesh::rect(22.86e-3, 10.16e-3, 80, 36);
    let mut em = ProblemSpec::new(Physics::Emw);
    em.study = Study::Modes { freq: 20e9, n_target: 0.9, evanescent: false, n_modes: 6 };
    run("EM modes 2D (Nedelec 2)", &mw, &em);
}
