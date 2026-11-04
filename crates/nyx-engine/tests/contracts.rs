//! Contracts of `nyx-engine` (ADR 19, ADR 37): the catalog of commands (names, access, schemas: what the MCP server turns into tools) and an
//! example project written by commands alone. The project layout is provisional until phase 8 derives the final one (gate G6), so its
//! exemplar is version 0.x.
//!
//! Rewrite after a deliberate change (with a `VERSIONS.md` line):
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-engine --test contracts`
use nyx_engine::*;
use nyx_physics::contract::check_exemplar;
use serde_json::{json, Value};

#[test]
fn the_command_catalog_is_pinned() {
    let e = Engine::new();
    let cat = serde_json::to_value(e.catalog()).unwrap();
    check_exemplar(env!("CARGO_MANIFEST_DIR"), "engine_commands.v1.json", &cat);
}

#[test]
fn an_example_project_made_by_commands_is_pinned() {
    let mut e = Engine::new();
    let mut call = |name: &str, a: Value| e.call(name, &a).unwrap_or_else(|err| panic!("{name}: {err}"));
    call("project_new", json!({ "title": "SOI strip", "h": 1e-7 }));
    call("set_param", json!({ "name": "w", "value": 5e-7 }));
    call("set_param", json!({ "name": "t", "value": 2.2e-7 }));
    call("add_rect", json!({ "name": "window", "min": [-2e-6, -1.5e-6], "max": [2e-6, 1.5e-6] }));
    call("add_rect", json!({ "name": "box", "min": [-2e-6, -1.5e-6], "max": [2e-6, 0.0] }));
    call("add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
    call(
        "add_overlay",
        json!({ "regions": [
        { "name": "air", "input": "window" }, { "name": "oxide", "input": "box" }, { "name": "silicon", "input": "core" } ] }),
    );
    call(
        "add_physics",
        json!({ "physics": "ems", "study": "modes", "study_params": { "freq": 1.934e14, "n_modes": 2 }, "order": 2 }),
    );
    call("set_material", json!({ "physics": "ems", "region": "silicon", "material": { "eps_r": 12.08 } }));
    call("check_project", json!({}));
    let v: Value = serde_json::from_str(&e.project().to_json()).unwrap();
    check_exemplar(env!("CARGO_MANIFEST_DIR"), "engine_project.v1.json", &v);
    // and the file reads back to the same project
    assert_eq!(Project::from_json(&e.project().to_json()).unwrap(), *e.project());
}
