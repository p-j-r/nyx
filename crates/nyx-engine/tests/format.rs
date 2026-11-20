// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The file format key (ADR 63): written in the project, its feature list and every problem; a file without it or of another format is
//! refused with a plain message, and nothing is converted.

use nyx_engine::Project;
use nyx_physics::ProblemSpec;
use serde_json::{json, Value};

fn small() -> Project {
    let mut p = Project::new("format test", 2);
    p.physics.push(ProblemSpec::new("dfs", 1, "electrostatic"));
    p
}

fn text_of(f: impl FnOnce(&mut Value)) -> String {
    let mut v: Value = serde_json::from_str(&small().to_json()).unwrap();
    f(&mut v);
    v.to_string()
}

#[test]
fn every_level_carries_the_format_key_and_the_project_reads_back() {
    let p = small();
    let v: Value = serde_json::from_str(&p.to_json()).unwrap();
    assert_eq!(v["format"], json!(1));
    assert_eq!(v["geometry"]["format"], json!(1));
    assert_eq!(v["physics"][0]["format"], json!(1));
    assert_eq!(Project::from_json(&p.to_json()).unwrap(), p);
}

#[test]
fn a_file_without_the_key_is_refused_at_each_level() {
    for (what, text) in [
        ("the project", text_of(|v| drop(v.as_object_mut().unwrap().remove("format")))),
        ("its geometry", text_of(|v| drop(v["geometry"].as_object_mut().unwrap().remove("format")))),
        ("problem 1", text_of(|v| drop(v["physics"][0].as_object_mut().unwrap().remove("format")))),
    ] {
        let e = Project::from_json(&text).unwrap_err().to_string();
        assert!(e.contains(what) && e.contains("\"format\""), "{what}: {e}");
    }
}

#[test]
fn a_file_of_another_format_is_refused_and_not_converted() {
    for (what, text) in [
        ("the project", text_of(|v| v["format"] = json!(2))),
        ("the geometry", text_of(|v| v["geometry"]["format"] = json!(2))),
        ("the 'dfs' problem", text_of(|v| v["physics"][0]["format"] = json!(0))),
    ] {
        let e = Project::from_json(&text).unwrap_err().to_string();
        assert!(e.contains(what) && e.contains("reads format 1"), "{what}: {e}");
    }
}

#[test]
fn the_example_files_are_format_1() {
    let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples");
    let mut n = 0;
    for e in std::fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        if p.to_string_lossy().ends_with(".nyx.json") {
            let t = std::fs::read_to_string(&p).unwrap();
            Project::from_json(&t).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
            n += 1;
        }
    }
    assert!(n >= 16, "found {n} example files");
}
