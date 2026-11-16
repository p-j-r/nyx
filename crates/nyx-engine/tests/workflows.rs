//! The command lists in `tools/workflows/` (productivity plan P0): JSON lines of `{cmd, args}` that the `workflow` example replays and that
//! `nyx check` / `nyx run` will take over (plan step S1). Every line must name a command of the catalog with an object for its arguments and
//! use only the replay keys the format defines; the start of the first list is replayed to show the format runs against an `Engine`.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

const REPLAY_KEYS: [&str; 7] = ["cmd", "args", "note", "job_arg", "after_stage", "delay_s", "expect_error"];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// (file name, line number, parsed line) of every command in every list.
fn lines() -> Vec<(String, usize, Value)> {
    let mut out = vec![];
    let dir = root().join("tools/workflows");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("tools/workflows")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "jsonl"))
        .collect();
    files.sort();
    assert!(!files.is_empty(), "no command lists in tools/workflows");
    for f in files {
        let name = f.file_name().unwrap().to_string_lossy().into_owned();
        for (i, l) in std::fs::read_to_string(&f).unwrap().lines().enumerate() {
            let l = l.trim();
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let v: Value = serde_json::from_str(l).unwrap_or_else(|e| panic!("{name}:{}: not JSON: {e}", i + 1));
            out.push((name.clone(), i + 1, v));
        }
    }
    out
}

#[test]
fn every_workflow_line_is_a_catalog_command_with_known_keys() {
    let catalog: Vec<String> = Engine::new().catalog().into_iter().map(|c| c.name).collect();
    let mut started_by_file: std::collections::BTreeMap<String, bool> = Default::default();
    for (file, n, v) in lines() {
        let at = format!("{file}:{n}");
        let obj = v.as_object().unwrap_or_else(|| panic!("{at}: a line is a JSON object"));
        for k in obj.keys() {
            assert!(REPLAY_KEYS.contains(&k.as_str()), "{at}: unknown key '{k}'");
        }
        let cmd = obj.get("cmd").and_then(Value::as_str).unwrap_or_else(|| panic!("{at}: no `cmd`"));
        assert!(catalog.iter().any(|c| c == cmd), "{at}: '{cmd}' is not an engine command");
        if let Some(a) = obj.get("args") {
            assert!(a.is_object(), "{at}: `args` is an object");
        }
        let started = started_by_file.entry(file.clone()).or_insert(false);
        if obj.contains_key("job_arg") || obj.contains_key("after_stage") {
            assert!(
                *started,
                "{at}: `job_arg` / `after_stage` need an earlier command that starts a job (`wait: false`)"
            );
        }
        if obj.get("args").and_then(|a| a.get("wait")) == Some(&json!(false)) {
            *started = true;
        }
    }
}

#[test]
fn the_start_of_the_first_workflow_runs_against_an_engine() {
    let mut e = Engine::new();
    let mut ran = 0;
    for (file, n, v) in lines().into_iter().filter(|(f, ..)| f == "w1_soi_width.jsonl").take(3) {
        let cmd = v["cmd"].as_str().unwrap();
        let mut args = v.get("args").cloned().unwrap_or(json!({}));
        // the lists name files relative to the repository root, where they are run from
        if let Some(p) = args.get("path").and_then(Value::as_str) {
            args["path"] = json!(root().join(p).to_string_lossy());
        }
        e.call(cmd, &args).unwrap_or_else(|err| panic!("{file}:{n}: {cmd}: {err}"));
        ran += 1;
    }
    assert_eq!(ran, 3);
    assert!(e.result("ems").is_some(), "the solve of the first list gave a result");
}
