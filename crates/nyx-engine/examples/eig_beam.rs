//! Eigenfrequencies of a block beam (fixed at x = 0) with a chosen linear-solver setting, to compare the multigrid/LOBPCG
//! path with the sparse factorisation and to time them:
//!
//!   cargo run --release -p nyx-engine --example eig_beam -- LX LY LZ NXxNYxNZ auto|direct|iter [order]
//!   cargo run --release -p nyx-engine --example eig_beam -- 1 0.1 0.1 40x8x8 auto      (the slender 70,000-unknown beam of bench.rs)
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
use femcore::problem::*;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let (lx, ly, lz): (f64, f64, f64) = (a[0].parse().unwrap(), a[1].parse().unwrap(), a[2].parse().unwrap());
    let n: Vec<usize> = a[3].split('x').map(|x| x.parse().unwrap()).collect();
    let solver = match a[4].as_str() {
        "direct" => LinearSolver::Direct,
        "iter" => LinearSolver::Iterative,
        _ => LinearSolver::Auto,
    };
    let order: usize = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(2);
    let m = benchmesh::block(lx, ly, lz, n[0], n[1], n[2]);
    let mut s = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
    s.solver = solver;
    s.order = order;
    let t = Instant::now();
    let r = nyx_engine::legacy::solve(&m, &s, &mut SolveState::default()).unwrap();
    println!("{:.1} s, {} cells, {} dofs, solver {:?}", t.elapsed().as_secs_f64(), m.n_cells(), r.dofs, solver);
    for (name, v) in &r.scalars {
        if let Scalar::Num(x) = v {
            println!("{name} {x:.12e}");
        }
    }
    for l in &r.log {
        println!("log: {l}");
    }
}
