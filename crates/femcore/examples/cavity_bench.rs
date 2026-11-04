//! Spherical PEC cavity timing: cargo run --release -p femcore --example cavity_bench -- h order
#[path = "../tests/common/solid3d/mod.rs"]
mod solid3d;
use solid3d::mesher3d::*;
use femcore::problem::*;
fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (h, order) = (a.first().copied().unwrap_or(0.012), a.get(1).copied().unwrap_or(2.0) as usize);
    let r = 0.05;
    let s3 = Spec3d { bodies: vec![Body3d { solid: Solid::Sphere { c: [0.0; 3], r }, domain: 1, h: None }], h, alpha: 0.3, snap_features: true, min_layers: 2.0 };
    let (mesh, _) = mesh_3d(&s3).unwrap();
    println!("{} tetrahedra", mesh.n_cells());
    let exact = 299792458.0 * 2.7437 / (2.0 * std::f64::consts::PI * r);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = order;
    s.materials.default.insert("eps_r".into(), 1.0);
    s.materials.default.insert("mu_r".into(), 1.0);
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: exact };
    let t = std::time::Instant::now();
    let res = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    println!("{} dofs in {:.1} s: {:?} (exact {:.5e}) | {}", res.dofs, t.elapsed().as_secs_f64(), res.solutions.iter().map(|x| x.value.unwrap() / 1e9).collect::<Vec<_>>(), exact, res.log.join(" | "));
}
