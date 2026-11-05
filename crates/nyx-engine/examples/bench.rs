//! Timings on larger generated models: cargo run --release -p nyx-engine --example bench
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
mod common;
use common::*;
use nyx_engine::facade::physics::ProblemSpec;
use nyx_engine::facade::SimplexMesh;
use serde_json::{json, Value};
use std::time::Instant;

fn run(name: &str, mesh: &SimplexMesh, spec: &ProblemSpec) {
    let t = Instant::now();
    match solve(mesh, spec, None) {
        Ok(r) => println!("{name:<44} {:>8} dofs  {:>7.2} s {}", r.dofs, t.elapsed().as_secs_f64(), r.log.first().cloned().unwrap_or_default()),
        Err(e) => println!("{name:<44} FAILED: {e}"),
    }
}

fn main() {
    let m2 = benchmesh::rect(1.0, 0.5, 300, 150);
    println!("2D mesh: {} cells", m2.n_cells());
    let mut es = spec("dfs", 2, "electrostatic", Value::Null);
    condition(&mut es, "left", "potential", json!({ "value": 1.0 }));
    condition(&mut es, "right", "ground", Value::Null);
    run("electrostatics 2D", &m2, &es);
    let mut so = spec("sms", 2, "static", Value::Null);
    steel(&mut so);
    fixed(&mut so, "left");
    load(&mut so, "right", [0.0, -1e5, 0.0]);
    run("solid static 2D", &m2, &so);
    // a boundary load belongs to the static study; the eigenfrequency study has the clamp only
    let mut eig = so.clone();
    eig.bcs.truncate(1);
    study(&mut eig, "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }));
    run("solid 6 eigenfrequencies 2D", &m2, &eig);

    let m3 = benchmesh::block(1.0, 0.1, 0.1, 40, 8, 8);
    println!("3D mesh: {} cells", m3.n_cells());
    let mut s3 = spec("sms", 2, "static", Value::Null);
    steel(&mut s3);
    fixed(&mut s3, "xmin");
    load(&mut s3, "xmax", [0.0, 0.0, -1e5]);
    run("solid static 3D", &m3, &s3);
    let mut eig3 = s3.clone();
    eig3.bcs.truncate(1);
    study(&mut eig3, "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }));
    run("solid 6 eigenfrequencies 3D", &m3, &eig3);

    for (nx, ny, nz) in [(40usize, 8usize, 8usize)] {
        let big = benchmesh::block(1.0, 0.1, 0.1, nx, ny, nz);
        let mut sb = s3.clone();
        for (name, s) in [("direct", "direct"), ("AMG-CG", "iterative")] {
            solver(&mut sb, s);
            if nx == 80 && s == "direct" {
                continue;
            }
            run(&format!("solid static 3D {nx}x{ny}x{nz} {name}"), &big, &sb);
        }
    }
    let mw = benchmesh::rect(22.86e-3, 10.16e-3, 80, 36);
    let mut em = spec("ems", 2, "modes", json!({ "freq": 20e9, "n_target": 0.9, "evanescent": false, "n_modes": 6 }));
    props(&mut em, &[("eps_r", 1.0), ("mu_r", 1.0)]);
    run("EM modes 2D (Nedelec 2)", &mw, &em);
}
