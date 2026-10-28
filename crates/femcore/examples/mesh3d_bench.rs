//! Speed of the 3D mesher: cargo run --release -p femcore --example mesh3d_bench [h ...]
use femcore::mesher3d::*;
use std::time::Instant;

fn main() {
    let hs: Vec<f64> = std::env::args().skip(1).filter_map(|a| a.parse().ok()).collect();
    for h in if hs.is_empty() { vec![0.1, 0.05, 0.03] } else { hs } {
        let s = Spec3d { bodies: vec![Body3d { solid: Solid::Sphere { c: [0.0; 3], r: 1.0 }, domain: 1, h: None }, Body3d { solid: Solid::Cylinder { base: [-2.0, 0.0, 0.0], axis: 0, r: 0.3, len: 4.0 }, domain: 0, h: None }], h, alpha: 0.3, snap_features: true, min_layers: 2.0 };
        let t = Instant::now();
        let (m, rep) = mesh_3d(&s).unwrap();
        let dt = t.elapsed().as_secs_f64();
        println!("h = {h}: {} tetrahedra, {} points, {} warped, {} cut points in {dt:.2} s ({:.2} million tets/s)", m.n_cells(), m.n_points(), rep.warped, rep.cut_points, m.n_cells() as f64 / dt / 1e6);
    }
}
