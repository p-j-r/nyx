// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Stage timings of a 3D solid static solve: NYX_PROFILE=1 cargo run --release -p nyx-engine --example profile3d -- NX NY NZ [p1|p2]
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
mod common;
use common::*;
use serde_json::Value;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: Vec<usize> = a.iter().take(3).map(|x| x.parse().unwrap()).collect();
    let mesh = benchmesh::block(n[0] as f64, n[1] as f64, n[2] as f64, n[0], n[1], n[2]);
    let order = if a.get(3).map(|x| x.as_str()) == Some("p1") { 1 } else { 2 };
    let mut s = spec("sms", order, "static", Value::Null);
    steel(&mut s);
    fixed(&mut s, "xmin");
    load(&mut s, "xmax", [0.0, 0.0, -1e5]);
    let t = std::time::Instant::now();
    let r = solve(&mesh, &s, None).unwrap();
    println!(
        "{} cells, {} dofs, total {:.2} s | {}",
        mesh.n_cells(),
        r.dofs,
        t.elapsed().as_secs_f64(),
        r.log.first().cloned().unwrap_or_default()
    );
}
