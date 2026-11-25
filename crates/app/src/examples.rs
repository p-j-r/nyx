// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-
//! Starting points: the example projects (files in the `examples` folder, found next to the program or in a parent folder of it, or at
//! `NYX_EXAMPLES`) and a blank simulation built from engine commands. An example is a project file like any other; nothing about a device is
//! in the program.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

/// One example file with the title written in it.
#[derive(Clone, Debug)]
pub struct Example {
    pub path: PathBuf,
    pub title: String,
}

/// The folders searched for examples, first with files wins: `NYX_EXAMPLES`, split/bundle layouts near `current_exe()`
/// (including symlinked executable targets via `canonicalize()`), and ancestor `examples` folders of the program or CWD.
pub fn dirs() -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = vec![];
    let mut push_unique = |p: PathBuf| {
        if !v.contains(&p) {
            v.push(p);
        }
    };
    if let Some(p) = std::env::var_os("NYX_EXAMPLES") {
        for dir in std::env::split_paths(&p) {
            push_unique(dir);
        }
    }
    let mut exe_dirs: Vec<PathBuf> = vec![];
    if let Ok(exe) = std::env::current_exe() {
        for cand in [Some(exe.clone()), exe.canonicalize().ok()] {
            if let Some(d) = cand.as_deref().and_then(Path::parent) {
                let dp = d.to_path_buf();
                if !exe_dirs.contains(&dp) {
                    exe_dirs.push(dp);
                }
            }
        }
    }
    for d in &exe_dirs {
        for rel in [
            d.join("examples"),
            d.join("..").join("examples"),
            d.join("..").join("Resources").join("examples"),
            d.join("..").join("share").join("nyx").join("examples"),
        ] {
            if has_projects(&rel) {
                push_unique(rel);
            }
        }
    }
    let mut starts = exe_dirs;
    starts.extend(std::env::current_dir().ok());
    for s in starts {
        let mut d: Option<&Path> = Some(&s);
        while let Some(dir) = d {
            let e = dir.join("examples");
            if has_projects(&e) {
                push_unique(e);
                break;
            }
            d = dir.parent();
        }
    }
    v
}

fn has_projects(d: &Path) -> bool {
    std::fs::read_dir(d).map(|r| r.filter_map(|e| e.ok()).any(|e| e.file_name().to_string_lossy().ends_with(".nyx.json"))).unwrap_or(false)
}

/// The examples of the first folder that has any, in file-name order.
pub fn list() -> Vec<Example> {
    for d in dirs() {
        let Ok(rd) = std::fs::read_dir(&d) else { continue };
        let mut files: Vec<PathBuf> = rd.filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.to_string_lossy().ends_with(".nyx.json")).collect();
        if files.is_empty() {
            continue;
        }
        files.sort();
        return files
            .into_iter()
            .map(|path| {
                let title = std::fs::read_to_string(&path)
                    .ok()
                    .and_then(|t| serde_json::from_str::<Value>(&t).ok())
                    .and_then(|v| v.get("title").and_then(Value::as_str).map(str::to_string))
                    .unwrap_or_else(|| path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default());
                Example { path, title }
            })
            .collect();
    }
    vec![]
}

/// A value for a property the descriptor gives no default for: a typical number, so that a new simulation runs before anything is edited.
pub fn typical(name: &str) -> f64 {
    match name {
        "eps_r" | "mu_r" => 1.0,
        "rho_v" | "Q" | "alpha" => 0.0,
        "nu" => 0.3,
        "E" => 200e9,
        "k" => 1.0,
        "rho" => 1000.0,
        "Cp" => 500.0,
        _ => 1.0,
    }
}

/// Starts a blank simulation in `engine`: a unit plate (one region `plate`, x-z) or, in 3D, the unit cube made by sweeping it along y,
/// with the physics `physics`, its first study and typical materials, and no conditions.
pub fn blank(engine: &mut Engine, physics: &str, three_d: bool) -> Result<(), String> {
    let mut call = |cmd: &str, args: Value| engine.call(cmd, &args).map(|_| ()).map_err(|e| format!("{cmd}: {e}"));
    call("project_new", json!({ "title": "New simulation", "dim": 2, "h": if three_d { 0.125 } else { 0.05 } }))?;
    call("add_rect", json!({ "name": "plate_box", "min": [0.0, 0.0], "max": [1.0, 1.0] }))?;
    call("add_overlay", json!({ "name": "regions", "regions": [{ "name": "plate", "input": "plate_box" }] }))?;
    if three_d {
        call("set_mesh", json!({ "sweep": [{ "thickness": 1.0, "n": 8 }] }))?;
    }
    let desc = engine.call("describe_physics", &json!({ "physics": physics })).map_err(|e| e.to_string())?;
    let studies = desc["studies"].as_array().cloned().unwrap_or_default();
    // a 3D electromagnetic problem is a cavity, a 2D one a waveguide cross-section: pick the first study that runs in this space
    let want = if three_d { "three" } else { "two" };
    let study = studies
        .iter()
        .find(|s| matches!(s["dims"].as_str(), Some(d) if d == want || d == "both"))
        .or_else(|| studies.first())
        .and_then(|s| s["kind"].as_str())
        .ok_or("this physics lists no study")?
        .to_string();
    let mut params = serde_json::Map::new();
    for p in desc["studies"].as_array().into_iter().flatten().filter(|s| s["kind"] == study.as_str()).flat_map(|s| s["params"].as_array().cloned().unwrap_or_default()) {
        // a required parameter has no usable default; give it the value the old window started with
        if p["required"] == true {
            let name = p["name"].as_str().unwrap_or_default().to_string();
            let v = match name.as_str() {
                "freq" => json!(10e9),
                "f_search" => json!(10e9),
                _ => p["kind"]["default"].clone(),
            };
            params.insert(name, v);
        }
    }
    engine
        .call("add_physics", &json!({ "physics": physics, "study": study, "study_params": params, "order": 2 }))
        .map_err(|e| format!("add_physics: {e}"))?;
    let mut mat = serde_json::Map::new();
    for p in desc["properties"].as_array().into_iter().flatten() {
        let name = p["name"].as_str().unwrap_or_default();
        if matches!(p["kind"].as_str(), Some("matrix3x6" | "dispersive")) {
            continue;
        }
        let v = p["default"].as_f64().unwrap_or_else(|| typical(name));
        mat.insert(name.to_string(), json!(v));
    }
    engine
        .call("set_material", &json!({ "physics": physics, "region": "plate", "material": Value::Object(mat) }))
        .map_err(|e| format!("set_material: {e}"))?;
    Ok(())
}

/// The JSON text of a blank simulation (`blank`) made on a scratch engine: the spec a new experiment starts from.
pub fn blank_json(physics: &str, three_d: bool) -> Result<String, String> {
    let mut e = Engine::new();
    blank(&mut e, physics, three_d)?;
    Ok(e.project().to_json())
}
