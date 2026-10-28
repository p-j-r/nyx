//! Stage timings of a 3D solid static solve: NYX_PROFILE=1 cargo run --release -p femcore --example profile3d -- NX NY NZ [p1|p2]
use femcore::mesh::SimplexMesh;
use femcore::problem::*;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: Vec<usize> = a.iter().take(3).map(|x| x.parse().unwrap()).collect();
    let mesh = SimplexMesh::block(n[0] as f64, n[1] as f64, n[2] as f64, n[0], n[1], n[2]);
    let mut s = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.order = if a.get(3).map(|x| x.as_str()) == Some("p1") { 1 } else { 2 };
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) }];
    let t = std::time::Instant::now();
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    println!("{} cells, {} dofs, total {:.2} s | {}", mesh.n_cells(), r.dofs, t.elapsed().as_secs_f64(), r.log.first().cloned().unwrap_or_default());
}
