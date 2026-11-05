//! Eigenfrequencies of a block beam (fixed at x = 0) with a chosen linear-solver setting, to compare the multigrid/LOBPCG
//! path with the sparse factorisation and to time them:
//!
//!   cargo run --release -p nyx-engine --example eig_beam -- LX LY LZ NXxNYxNZ auto|direct|iter [order]
//!   cargo run --release -p nyx-engine --example eig_beam -- 1 0.1 0.1 40x8x8 auto      (the slender 70,000-unknown beam of bench.rs)
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
mod common;
use common::*;
use serde_json::json;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let (lx, ly, lz): (f64, f64, f64) = (a[0].parse().unwrap(), a[1].parse().unwrap(), a[2].parse().unwrap());
    let n: Vec<usize> = a[3].split('x').map(|x| x.parse().unwrap()).collect();
    let linear = match a[4].as_str() {
        "direct" => "direct",
        "iter" => "iterative",
        _ => "auto",
    };
    let order: u8 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(2);
    let m = benchmesh::block(lx, ly, lz, n[0], n[1], n[2]);
    let mut s = spec("sms", order, "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }));
    steel(&mut s);
    solver(&mut s, linear);
    fixed(&mut s, "xmin");
    let t = Instant::now();
    match solve(&m, &s, None) {
        Ok(r) => {
            println!("{} cells, {} dofs, {:.2} s", m.n_cells(), r.dofs, t.elapsed().as_secs_f64());
            println!("frequencies (Hz): {:?}", r.solutions.iter().filter_map(|s| s.value).collect::<Vec<_>>());
            for l in &r.log {
                println!("  {l}");
            }
        }
        Err(e) => println!("FAILED: {e}"),
    }
}
