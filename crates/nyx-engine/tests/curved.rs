// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A dielectric sphere in a periodic cube at k = 0 with the `curved` mesh setting on and off (ADR 82): the reduced system keeps its size,
//! the three uniform-field modes stay at zero, and the frequencies move the expected way (a curved sphere holds more dielectric than the
//! polyhedron inscribed in it, so its modes are lower). The settings and the report travel through the engine commands.
use nyx_engine::Engine;
use serde_json::{json, Value};

const C: f64 = 299_792_458.0;

fn project(curved: bool, h: f64) -> String {
    let a = 1e-6;
    let f = |i: u32, n: &str, op: Value| json!({ "id": i, "name": n, "op": op, "origin": "Code", "suppressed": false, "tags": [] });
    let mut feats = vec![
        f(1, "cell", json!({ "Box": { "min": [0, 0, 0], "max": [a, a, a] } })),
        f(2, "ball", json!({ "Sphere": { "center": [a / 2.0, a / 2.0, a / 2.0], "r": 0.3 * a } })),
        f(
            3,
            "overlay1",
            json!({ "Overlay": { "naming": "Feature", "pins": [], "regions": [
            { "input": 1, "kind": "Material", "name": "bg", "pinned_id": null }, { "input": 2, "kind": "Material", "name": "sph", "pinned_id": null } ] } }),
        ),
    ];
    for (j, ax) in ["x", "y", "z"].iter().enumerate() {
        let mut t = [0.0; 3];
        t[j] = a;
        feats.push(f(4 + j as u32, &format!("periodic{}", j + 1), json!({ "Periodic": { "sel_a": { "Name": format!("cell.{ax}min") }, "sel_b": { "Name": format!("cell.{ax}max") }, "translation": t } })));
    }
    let bcs: Vec<Value> = ["x", "y", "z"]
        .iter()
        .map(|ax| json!({ "kind": "periodic", "on": format!("name(cell.{ax}max)"), "params": { "dest": format!("name(cell.{ax}min)"), "phase": 0.0 } }))
        .collect();
    let mut mesh = json!({ "h": h * a });
    if curved {
        mesh["curved"] = json!(true);
    }
    json!({ "active": "ems", "format": 1, "geometry": { "dim": 3, "features": feats, "format": 1, "params": [] }, "mesh": mesh,
        "physics": [{ "bcs": bcs, "format": 1, "materials": { "bg": { "eps_r": 1.0 }, "sph": { "eps_r": 12.96 } }, "order": 2, "params": {},
            "physics": "ems", "study": { "kind": "cavity", "params": { "f_search": 0.05 * C / a, "n_modes": 8 } } }],
        "title": "curved periodic sphere" })
    .to_string()
}

fn run(curved: bool) -> (Vec<f64>, i64, Value) {
    let mut e = Engine::new();
    e.call("project_load", &json!({ "json": project(curved, 0.2) })).unwrap();
    let built = e.call("mesh_build", &json!({})).unwrap();
    let r = e.call("solve_run", &json!({})).unwrap();
    let info = e.call("result_info", &json!({ "physics": "ems" })).unwrap();
    let vals: Vec<f64> =
        info[0]["solutions"].as_array().unwrap().iter().map(|s| s["value"].as_f64().unwrap()).collect();
    let dofs = info[0]["dofs"].as_i64().unwrap();
    let _ = r;
    (vals, dofs, built)
}

#[test]
fn a_curved_dielectric_sphere_in_a_periodic_cube_lowers_the_frequencies_and_keeps_the_null_space() {
    let (fs, ds, mesh_s) = run(false);
    let (fc, dc, mesh_c) = run(true);
    eprintln!("straight {fs:?} ({ds} dofs)\ncurved   {fc:?} ({dc} dofs)\nmesh report {mesh_c}");
    assert!(mesh_s.get("curved").is_none() || mesh_s["curved"].is_null(), "{mesh_s}");
    assert_eq!(ds, dc, "the numbering does not change");
    assert_eq!(fs.len(), fc.len());
    // the uniform-field modes sit at (numerically) zero in both, far below the first physical band
    for f in [&fs, &fc] {
        let physical = f.iter().cloned().fold(f64::INFINITY, |m, v| if v > 0.02 * C / 1e-6 { m.min(v) } else { m });
        assert!(f.iter().filter(|&&v| v < 0.01 * physical).count() >= 3, "{f:?}");
    }
    // curved holds the true sphere (more dielectric): every physical mode is lower
    let phys = |f: &Vec<f64>| -> Vec<f64> { f.iter().cloned().filter(|&v| v > 0.02 * C / 1e-6).collect() };
    let (ps, pc) = (phys(&fs), phys(&fc));
    assert!(!ps.is_empty() && ps.len() == pc.len());
    for (s, c) in ps.iter().zip(&pc) {
        assert!(c < s, "curved {c} not below straight {s}");
        assert!((s / c - 1.0) < 0.05, "the change is a few percent at most: {s} vs {c}");
    }
}
