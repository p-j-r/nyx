//! What an agent lost time on in the Antigravity benchmark session (L1 to L3, ADR 97): a wildcard material name that failed only in `mesh_build`
//! with `physics_sizing`, a mesh far larger than needed that was found out by solving it, and a convergence check that could not be asked
//! what it would cost. The tests use small 2D models; nothing here needs a long solve.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

fn c(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}: {er}"))
}

/// A wire in silica with absorbing layers on both sides, the materials given for `BOX*` and `clad*` (a name with a star stands for the region and
/// its absorbing pieces), the mesh sized from the waves.
fn wire_with_pml_and_stars() -> Engine {
    let mut e = Engine::new();
    c(&mut e, "project_new", json!({ "title": "pml and stars", "dim": 2, "h": 3e-7 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [-3e-6, -2e-6], "max": [3e-6, 0.0] }));
    c(&mut e, "add_rect", json!({ "name": "clad", "min": [-3e-6, 0.0], "max": [3e-6, 1.5e-6] }));
    c(&mut e, "add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 3e-7] }));
    c(&mut e, "add_rect", json!({ "name": "window", "min": [-2.5e-6, -2e-6], "max": [2.5e-6, 1.5e-6] }));
    c(
        &mut e,
        "add_overlay",
        json!({ "name": "regions", "regions": [{ "name": "BOX", "input": "box" }, { "name": "clad", "input": "clad" }, { "name": "Si", "input": "core" }] }),
    );
    c(
        &mut e,
        "add_pml",
        json!({ "name": "pml", "window": "window", "sides": ["XMin", "XMax"], "thickness": 5e-7, "layers": 3.0, "grading": 2.0, "regions": ["BOX", "clad"] }),
    );
    c(
        &mut e,
        "add_physics",
        json!({ "physics": "ems", "study": "modes", "order": 2, "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 1 } }),
    );
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "BOX*", "material": "silica" }));
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "clad*", "material": "silica" }));
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "Si", "material": "silicon" }));
    c(&mut e, "add_pml_conditions", json!({ "physics": "ems", "strength": 3.0 }));
    e
}

#[test]
fn a_material_given_for_a_starred_region_sizes_the_mesh_of_a_model_with_absorbing_layers() {
    let mut e = wire_with_pml_and_stars();
    c(&mut e, "set_mesh", json!({ "physics_sizing": true, "adapt_scale": 1.0 }));
    // before: the mesher looked up the region `BOX*` literally and refused ("the sizing names region 'BOX*', which the model does not have")
    let m = c(&mut e, "mesh_build", json!({}));
    assert_eq!(m["state"], "done", "{m}");
    // the answer says what a solve of this mesh is expected to cost: nothing was solved on this machine yet, so no number is invented
    let cost = &m["cost"];
    assert!(cost["cells"].as_u64().unwrap() > 100, "{cost}");
    assert_eq!(cost["dim"], 2);
    assert!(cost["solve"].is_null(), "no solve has been run: no estimate: {cost}");
    assert!(cost.get("warning").is_none(), "a small mesh is not a warning: {cost}");
}

#[test]
fn a_convergence_dry_run_tells_the_cells_of_both_meshes_and_solves_nothing() {
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("02_cantilever") }));
    let p = c(&mut e, "convergence_check", json!({ "physics": "sms", "factor": 0.6, "dry_run": true }));
    assert_eq!(p["dry_run"], true, "{p}");
    let (coarse, refined) = (p["cells"]["coarse"].as_u64().unwrap(), p["cells"]["refined"].as_u64().unwrap());
    assert!(refined > coarse, "{p}");
    assert!(p["refined_over_present"].as_f64().unwrap() > 1.0);
    assert!(p["seconds_estimate"].is_null(), "no solve has been run on this engine: no estimate: {p}");
    assert!(p["cost"]["refined"]["solve"].is_null());
    // nothing was solved and nothing queued
    assert!(e.call("result_info", &json!({})).is_err(), "no result: nothing was solved");
    assert!(c(&mut e, "job_list", json!({})).as_array().unwrap().iter().all(|j| j["kind"] != "convergence"));
}

#[test]
fn a_convergence_check_larger_than_max_cells_is_refused_before_any_solve() {
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("02_cantilever") }));
    let p = c(&mut e, "convergence_check", json!({ "physics": "sms", "factor": 0.6, "dry_run": true }));
    let refined = p["cells"]["refined"].as_u64().unwrap();
    let err = e
        .call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "max_cells": refined - 1 }))
        .unwrap_err()
        .to_string();
    assert!(err.contains("max_cells") && err.contains(&refined.to_string()), "{err}");
    assert!(err.contains("not started"), "{err}");
    // the same limit above the size lets it run
    let ok = c(
        &mut e,
        "convergence_check",
        json!({ "physics": "sms", "factor": 0.8, "max_cells": 10 * refined, "tolerance": 1.0 }),
    );
    assert!(!ok["rows"].as_array().unwrap().is_empty(), "{ok}");
    // bad limits are refused in plain words
    for bad in [
        json!({ "max_cells": 0 }),
        json!({ "max_cells": -3 }),
        json!({ "max_seconds": 0 }),
        json!({ "max_seconds": "soon" }),
    ] {
        let mut a = json!({ "physics": "sms" });
        for (k, v) in bad.as_object().unwrap() {
            a[k] = v.clone();
        }
        assert!(e.call("convergence_check", &a).is_err(), "{a}");
    }
}

#[test]
fn a_lossless_guide_on_a_rotated_crystal_converges_although_its_imaginary_part_is_rounding() {
    let mut e = wire_with_pml_and_stars();
    // sapphire with a tilted axis has off-diagonal permittivity: the complex solver runs, and a lossless guide gets an imaginary part of 1e-11
    // and a quality factor of 1e10 made of it (the benchmark session of the SOS wire)
    c(
        &mut e,
        "set_material",
        json!({ "physics": "ems", "region": "BOX*", "material": { "material": "sapphire", "euler_deg": [0.0, 57.6, 0.0] } }),
    );
    c(&mut e, "set_mesh", json!({ "physics_sizing": true, "adapt_scale": 1.0 }));
    let s = c(&mut e, "solve_run", json!({}));
    let scalars = s["results"][0]["scalars"].to_string();
    assert!(scalars.contains("Q 1 = Re/(2|Im|) of neff"), "the solve reports the noise quality factor: {scalars}");
    let r = c(&mut e, "convergence_check", json!({ "physics": "ems", "factor": 0.7, "tolerance": 5e-2 }));
    assert_eq!(r["noise_modes"], json!([0]), "{r}");
    assert_eq!(r["converged"], true, "the real part converged; the noise does not count: {r}");
    let rows = r["rows"].as_array().unwrap();
    let noisy: Vec<&str> =
        rows.iter().filter(|x| x["noise"] == true).map(|x| x["quantity"].as_str().unwrap()).collect();
    assert!(noisy.contains(&"eigenvalue[0].im"), "{noisy:?}");
    assert!(noisy.iter().any(|q| q.starts_with("scalar Q")), "{noisy:?}");
    let worst = r["worst_quantity"].as_str().unwrap();
    assert!(!worst.contains(".im") && !worst.contains("scalar Q"), "the worst row is a real one: {worst}");
    // a tensor permittivity is refused as a weight with the names of its entries (the guide says to read them)
    let err = e
        .call(
            "confinement",
            &json!({ "fields": ["Ex", "Ey", "Ez"], "weight": "eps_r", "region": "Si", "physics": "ems" }),
        )
        .unwrap_err()
        .to_string();
    assert!(err.contains("eps_r_xx") && err.contains("tensor"), "{err}");
    // and the contraction the optical-modes card gives for the energy works, over the region and over all of it
    let w = "re(conj(Ex)*(eps_r_xx*Ex + eps_r_xy*Ey + eps_r_xz*Ez) + conj(Ey)*(eps_r_xy*Ex + eps_r_yy*Ey + eps_r_yz*Ez) + conj(Ez)*(eps_r_xz*Ex + eps_r_yz*Ey + eps_r_zz*Ez))";
    let part = e.call("integrate", &json!({ "expr": w, "physics": "ems", "region": "Si" }));
    let all = e.call("integrate", &json!({ "expr": w, "physics": "ems" }));
    assert!(part.is_ok() && all.is_ok(), "the card's formula runs: {part:?} {all:?}");
}

#[test]
fn a_sweep_and_a_mesh_say_what_they_are_expected_to_take_once_the_machine_has_solved_a_few() {
    use nyx_engine::estimate::Row;
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("02_cantilever") }));
    let m = c(&mut e, "mesh_build", json!({}));
    let cells = m["cost"]["cells"].as_u64().unwrap();
    assert!(m["cost"]["solve"].is_null(), "no runs yet: {m}");
    // three earlier solves of other sizes, 2 s per thousand cells
    for k in [1u64, 2, 4] {
        e.history().add(Row {
            chain: "sms".into(),
            dim: 2,
            cells: cells * k / 2,
            dofs: cells * k,
            peak_working_set: 1_000_000 * k,
            peak_commit: 1_000_000 * k,
            seconds: 0.001 * (cells * k / 2) as f64,
        });
    }
    let m = c(&mut e, "mesh_build", json!({ "wait": true }));
    let per = m["cost"]["solve"]["seconds"].as_f64().expect("an estimate from three runs");
    assert!(per > 0.0 && m["cost"]["solve"]["runs"] == 3, "{m}");
    // a sweep of three points on this mesh: the answer carries the expected total
    let r = c(
        &mut e,
        "sweep_run",
        json!({ "series": { "cmd": "set_mesh", "args": {}, "field": "h", "values": [4e-3, 4e-3, 4e-3] } }),
    );
    let cost = &r["cost"];
    assert_eq!(cost["points"], 3, "{r}");
    let total = cost["seconds_estimate"].as_f64().expect("a total");
    assert!((total - 3.0 * per).abs() < 1e-9 * total.max(1.0), "{cost}");
}
