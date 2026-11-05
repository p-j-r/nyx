//! The layering of the workspace (ADR 26; ROUGH_SPOTS R18, R15, R16): which crate may depend on which. The table below is the target
//! structure of `docs/architecture/RESTRUCTURE_PLAN.md` section 2, plus a short list of edges that exist today and that a named phase removes.
//!
//! * a dependency between workspace crates that is in neither list fails the test: change the structure on purpose, in the same commit as
//!   this table, and say why in the commit message;
//! * an entry of `TEMPORARY` that no longer exists fails too, so the table cannot keep an edge that is gone: delete the row.
//!
//! Only normal dependencies count (not dev- or build-dependencies).

use std::collections::{BTreeMap, BTreeSet};

/// crate -> the workspace crates it may depend on in the finished structure
const ALLOWED: &[(&str, &[&str])] = &[
    ("geom", &[]),
    ("femcore", &["geom"]),
    ("nyx-physics", &["femcore"]),
    // the physics crates keep femcore: the kernels use its mesh, fem and linalg
    ("nyx-dfs", &["nyx-physics", "femcore"]),
    ("nyx-sms", &["nyx-physics", "femcore"]),
    ("nyx-ems", &["nyx-physics", "femcore"]),
    ("nyx-analysis", &["nyx-physics", "femcore"]), // femcore: the mesh geometry and the quadrature rule the modes are sampled with
    ("nyx-store", &[]),
    ("nyx-engine", &["nyx-physics", "nyx-dfs", "nyx-sms", "nyx-ems", "nyx-analysis", "nyx-store", "femcore", "geom"]),
    ("nyx-connect", &[]),
    ("nyx-mcp", &["nyx-engine"]),
    ("nyx", &["nyx-engine", "nyx-connect"]), // crates/app
];

/// edges that exist today and the phase that removes each
const TEMPORARY: &[(&str, &str, &str)] = &[
    ("nyx-store", "femcore", "phase 7: metrics_of and the concept card move to the engine (R21, R25)"),
    ("nyx", "femcore", "phase 7: the window calls engine commands only"),
    ("nyx", "geom", "phase 7: the window calls engine commands only"),
    ("nyx", "nyx-store", "phase 7: persistence goes through the engine"),
];

/// a crate that is not in `ALLOWED` may not exist (it must be put in the table)
const ONLY_TEMPORARY_CRATES: &[&str] = &[];

fn workspace_graph() -> BTreeMap<String, BTreeSet<String>> {
    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());
    let out = std::process::Command::new(cargo)
        .args(["metadata", "--format-version", "1", "--no-deps", "--offline"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("run cargo metadata");
    assert!(out.status.success(), "cargo metadata failed: {}", String::from_utf8_lossy(&out.stderr));
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    let members: BTreeSet<String> =
        v["packages"].as_array().unwrap().iter().map(|p| p["name"].as_str().unwrap().to_string()).collect();
    let mut g = BTreeMap::new();
    for p in v["packages"].as_array().unwrap() {
        let name = p["name"].as_str().unwrap().to_string();
        let deps: BTreeSet<String> = p["dependencies"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|d| d["kind"].is_null() && d["path"].is_string())
            .map(|d| d["name"].as_str().unwrap().to_string())
            .filter(|d| members.contains(d))
            .collect();
        g.insert(name, deps);
    }
    g
}

#[test]
fn the_workspace_follows_the_layering_table() {
    let g = workspace_graph();
    let allowed: BTreeMap<&str, BTreeSet<&str>> =
        ALLOWED.iter().map(|(c, d)| (*c, d.iter().copied().collect())).collect();
    let mut errors = vec![];
    for (name, deps) in &g {
        let known = allowed.contains_key(name.as_str()) || ONLY_TEMPORARY_CRATES.contains(&name.as_str());
        if !known {
            errors.push(format!(
                "crate {name} is not in the layering table (ALLOWED): add it with the dependencies it may have"
            ));
            continue;
        }
        for d in deps {
            let ok = allowed.get(name.as_str()).is_some_and(|a| a.contains(d.as_str()))
                || TEMPORARY.iter().any(|(c, dd, _)| c == name && dd == d);
            if !ok {
                errors.push(format!("{name} depends on {d}, which the layering table does not allow"));
            }
        }
    }
    for (c, d, why) in TEMPORARY {
        match g.get(*c) {
            None => {
                errors.push(format!("TEMPORARY names the crate {c}, which no longer exists: delete the row ({why})"))
            }
            Some(deps) if !deps.contains(*d) => {
                errors.push(format!("{c} no longer depends on {d}: delete the TEMPORARY row ({why})"))
            }
            _ => {}
        }
    }
    for c in allowed.keys() {
        if !g.contains_key(*c) {
            errors.push(format!("the table lists {c}, which is not a workspace crate"));
        }
    }
    assert!(errors.is_empty(), "layering violations:\n{}", errors.join("\n"));
}

#[test]
fn no_physics_depends_on_another_physics_and_only_the_engine_names_them() {
    let g = workspace_graph();
    let physics = ["nyx-dfs", "nyx-sms", "nyx-ems"];
    for (name, deps) in &g {
        for d in deps {
            if physics.contains(&d.as_str()) {
                assert!(name == "nyx-engine", "{name} depends on the physics crate {d}: only nyx-engine may (nothing outside the physics crates names a physics)");
            }
        }
    }
    assert!(g["femcore"].iter().all(|d| !d.starts_with("nyx-")), "femcore must not depend on a nyx crate");
}
