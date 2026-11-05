//! Spherical PEC cavity timing: cargo run --release -p nyx-engine --example cavity_bench -- h order
mod common;
use common::*;
use serde_json::json;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (h, order) = (a.first().copied().unwrap_or(0.012), a.get(1).copied().unwrap_or(2.0) as u8);
    let r = 0.05;

    let mesh = femcore::meshing::geomdrive3::sphere_cavity(r, h).unwrap().mesh;
    println!("{} tetrahedra", mesh.n_cells());
    let exact = 299792458.0 * 2.7437 / (2.0 * std::f64::consts::PI * r);
    let mut s = spec("ems", order, "cavity", json!({ "n_modes": 4, "f_search": exact }));
    props(&mut s, &[("eps_r", 1.0), ("mu_r", 1.0)]);
    let t = std::time::Instant::now();
    let res = solve(&mesh, &s, None).unwrap();
    println!(
        "{} dofs in {:.1} s: {:?} (exact {:.5e}) | {}",
        res.dofs,
        t.elapsed().as_secs_f64(),
        res.solutions.iter().map(|x| x.value.unwrap() / 1e9).collect::<Vec<_>>(),
        exact,
        res.log.join(" | ")
    );
}
