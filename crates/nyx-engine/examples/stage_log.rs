// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Prints the plain-language stage log of a few small solves, to see what the simulator reports for each kind of study
//! (the text of docs/HOW_IT_WORKS.md quotes these):
//!
//!   cargo run --release -p nyx-engine --example stage_log
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
mod common;
use common::*;
use nyx_engine::facade::physics::{ProblemSpec, Scalar};
use nyx_engine::facade::SimplexMesh;
use serde_json::{json, Value};

fn show(title: &str, mesh: &SimplexMesh, spec: &ProblemSpec) {
    println!("=== {title}");
    match solve(mesh, spec, None) {
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
    let rect = benchmesh::rect(1.0, 0.5, 20, 10);
    let mut es = spec("dfs", 2, "electrostatic", Value::Null);
    props(&mut es, &[("eps_r", 4.0)]);
    condition(&mut es, "left", "potential", json!({ "value": 1.0 }));
    condition(&mut es, "top", "ground", Value::Null);
    show("electrostatics, 2D", &rect, &es);

    let mut so = spec("sms", 2, "static", Value::Null);
    steel(&mut so);
    fixed(&mut so, "left");
    load(&mut so, "right", [0.0, -1e5, 0.0]);
    show("solid mechanics, static, 2D", &rect, &so);
    study(&mut so, "eigenfrequency", json!({ "n_modes": 4, "f_search": 0.0 }));
    show("solid mechanics, 4 eigenfrequencies, 2D", &rect, &so);

    let beam = benchmesh::block(1.0, 0.1, 0.1, 40, 8, 8);
    so.bcs.truncate(1);
    so.bcs[0].on = "xmin".into();
    study(&mut so, "eigenfrequency", json!({ "n_modes": 4, "f_search": 0.0 }));
    show(
        "solid mechanics, 4 eigenfrequencies, 3D slender beam (multigrid LOBPCG is tried above 80,000 unknowns)",
        &beam,
        &so,
    );
    let big = benchmesh::block(1.0, 1.0, 1.0, 14, 14, 14);
    study(&mut so, "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }));
    solver(&mut so, "iterative");
    show("solid mechanics, 6 eigenfrequencies, 3D cube, multigrid LOBPCG forced", &big, &so);

    let guide = benchmesh::rect(22.86e-3, 10.16e-3, 40, 18);
    let mut em = spec("ems", 2, "modes", json!({ "freq": 20e9, "n_target": 0.9, "evanescent": false, "n_modes": 3 }));
    props(&mut em, &[("eps_r", 1.0), ("mu_r", 1.0)]);
    show("electromagnetic modes of a rectangular waveguide (WR-90), 2D cross-section", &guide, &em);
}
