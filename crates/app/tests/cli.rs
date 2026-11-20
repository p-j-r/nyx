// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! `nyx check` and `nyx run` (productivity plan S1, ADR 89): the command line checks a list before it runs, names the bad line, exits 2 and runs
//! nothing; a good list is planned without meshing and runs through the engine's commands.

use std::path::PathBuf;
use std::process::{Command, Output};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..").canonicalize().unwrap()
}

fn nyx(args: &[&str]) -> Output {
    // run from the repository root, as the lists' relative example paths expect
    Command::new(env!("CARGO_BIN_EXE_nyx")).args(args).current_dir(root()).output().expect("the nyx program")
}

fn tmp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("nyx-cli-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn text(o: &[u8]) -> String {
    String::from_utf8_lossy(o).into_owned()
}

#[test]
fn check_prints_the_plan_of_a_workflow_and_exits_zero() {
    let o = nyx(&["check", "tools/workflows/w1_soi_width.jsonl"]);
    assert_eq!(o.status.code(), Some(0), "{}", text(&o.stderr));
    let out = text(&o.stdout);
    assert!(out.contains("6 mesh build(s) and 6 solve job(s)"), "{out}");
    assert!(out.contains("nothing was meshed, solved or written"));
}

#[test]
fn check_names_a_bad_line_and_exits_two() {
    let dir = tmp("bad");
    let f = dir.join("bad.jsonl");
    std::fs::write(
        &f,
        "{\"cmd\":\"project_open\",\"args\":{\"path\":\"examples/15_soi_wire.nyx.json\"}}\n{\"cmd\":\"set_param\",\"args\":{\"name\":\"w\"}}\n{\"cmd\":\"solve_runn\"}\n",
    )
    .unwrap();
    let o = nyx(&["check", f.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2));
    let err = text(&o.stderr);
    assert!(err.contains("line 2, set_param") && err.contains("missing 'value'"), "{err}");
    assert!(err.contains("line 3, solve_runn") && err.contains("did you mean solve_run"), "{err}");
    // `--json` is for an agent: the problems are in the report
    let j = nyx(&["check", f.to_str().unwrap(), "--json"]);
    let v: serde_json::Value = serde_json::from_str(&text(&j.stdout)).expect("a JSON report");
    assert_eq!(v["problems"].as_array().unwrap().len(), 2);
    assert_eq!(v["problems"][0]["line"], 2);
}

#[test]
fn run_refuses_a_bad_list_before_running_anything_and_runs_a_good_one() {
    let dir = tmp("run");
    let ws = dir.join("ws");
    let bad = dir.join("bad.jsonl");
    std::fs::write(
        &bad,
        "{\"cmd\":\"project_open\",\"args\":{\"path\":\"examples/15_soi_wire.nyx.json\"}}\n{\"cmd\":\"solve_run\"}\n{\"cmd\":\"set_param\",\"args\":{\"name\":\"w\"}}\n",
    )
    .unwrap();
    let o = nyx(&["run", bad.to_str().unwrap(), "--workspace", ws.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2), "{}", text(&o.stderr));
    assert!(!text(&o.stdout).contains("solve_run"), "a command ran: {}", text(&o.stdout));
    assert!(!ws.exists(), "a workspace was made for a list that never ran");

    let good = dir.join("good.jsonl");
    std::fs::write(
        &good,
        "{\"cmd\":\"project_open\",\"args\":{\"path\":\"examples/15_soi_wire.nyx.json\"}}\n{\"cmd\":\"solve_run\"}\n{\"cmd\":\"result_info\",\"args\":{\"physics\":\"ems\"}}\n",
    )
    .unwrap();
    let o = nyx(&["run", good.to_str().unwrap(), "--workspace", ws.to_str().unwrap(), "--json"]);
    assert_eq!(o.status.code(), Some(0), "{}", text(&o.stderr));
    let rows: Vec<serde_json::Value> = text(&o.stdout).lines().map(|l| serde_json::from_str(l).unwrap()).collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[1]["cmd"], "solve_run");
    assert!(rows[1]["reply"]["run"]["id"].is_string(), "the solve was not recorded: {}", rows[1]);
    assert!(ws.exists());
}

#[test]
fn run_exits_one_when_a_command_fails_and_names_it() {
    let dir = tmp("fail");
    let f = dir.join("fail.jsonl");
    // passes the check (the feature name is just text) and fails when the engine runs it? the plan catches it too: skip the check
    std::fs::write(
        &f,
        "{\"cmd\":\"project_open\",\"args\":{\"path\":\"examples/15_soi_wire.nyx.json\"}}\n{\"cmd\":\"remove_feature\",\"args\":{\"name\":\"nothing_like_this\"}}\n{\"cmd\":\"project_info\"}\n",
    )
    .unwrap();
    let o = nyx(&["run", f.to_str().unwrap(), "--no-check"]);
    assert_eq!(o.status.code(), Some(1));
    assert!(text(&o.stderr).contains("line 2, remove_feature"), "{}", text(&o.stderr));
    assert!(!text(&o.stdout).contains("project_info"), "a command after the failure ran");
    // with the check on, the plan finds the same line and nothing runs
    let o = nyx(&["run", f.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2));
    assert!(text(&o.stderr).contains("line 2, remove_feature"));
}

#[test]
fn the_usage_errors_exit_two() {
    assert_eq!(nyx(&["check"]).status.code(), Some(2));
    assert_eq!(nyx(&["run", "x.jsonl", "--nope"]).status.code(), Some(2));
    assert_eq!(nyx(&["check", "no_such_file.jsonl"]).status.code(), Some(2));
}

#[test]
fn api_list_and_show_print_the_reference_of_the_catalog() {
    let o = nyx(&["api", "list"]);
    assert_eq!(o.status.code(), Some(0));
    let out = text(&o.stdout);
    assert!(out.lines().any(|l| l.starts_with("set_param") && l.contains("edit")), "{out}");
    assert!(out.lines().count() > 80, "every command is listed");
    let j: serde_json::Value = serde_json::from_str(&text(&nyx(&["api", "list", "--json"]).stdout)).unwrap();
    assert_eq!(j.as_array().unwrap().len(), out.lines().count());

    let o = nyx(&["api", "show", "set_param"]);
    assert_eq!(o.status.code(), Some(0));
    // the page is the file in docs/api, byte for byte
    let file = std::fs::read_to_string(root().join("docs/api/set_param.md")).unwrap().replace("\r\n", "\n");
    assert_eq!(text(&o.stdout), file);
    let j: serde_json::Value = serde_json::from_str(&text(&nyx(&["api", "show", "set_param", "--json"]).stdout)).unwrap();
    assert_eq!(j["example"]["cmd"], "set_param");
    assert!(j["schema"]["properties"]["value"].is_object());

    let o = nyx(&["api", "show", "set_parm"]);
    assert_eq!(o.status.code(), Some(2));
    assert!(text(&o.stderr).contains("did you mean set_param"), "{}", text(&o.stderr));
    assert_eq!(nyx(&["api"]).status.code(), Some(2));
    assert_eq!(nyx(&["api", "frob"]).status.code(), Some(2));
}

#[test]
fn the_program_answers_job_as_a_worker_and_never_opens_a_window() {
    // no folder: the usage, exit 2, at once (a window would not exit)
    let o = nyx(&["--job"]);
    assert_eq!(o.status.code(), Some(2), "{}", text(&o.stderr));
    // a real job, as an engine writes it: the project, the materials of its library and the chain; the answer is on stdout, the results in a file
    let dir = tmp("job");
    std::fs::copy(root().join("examples/01_capacitor.nyx.json"), dir.join("project.nyx.json")).unwrap();
    let engine = nyx_engine::Engine::new();
    let lib = engine.library().unwrap();
    let job = serde_json::json!({ "targets": ["dfs"], "materials": nyx_engine::worker::material_texts(&lib) });
    std::fs::write(dir.join("job.json"), job.to_string()).unwrap();
    // the worker ends when its standard input closes, so the parent keeps it open while it reads
    let mut child = Command::new(env!("CARGO_BIN_EXE_nyx"))
        .args(["--job", dir.to_str().unwrap()])
        .current_dir(root())
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("the nyx program");
    let mut out = String::new();
    std::io::Read::read_to_string(child.stdout.as_mut().unwrap(), &mut out).unwrap();
    let status = child.wait().unwrap();
    assert_eq!(status.code(), Some(0), "{out}");
    let done = out
        .lines()
        .filter_map(|l| l.strip_prefix(nyx_engine::worker::PREFIX))
        .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
        .find_map(|m| m.get("done").cloned())
        .unwrap_or_else(|| panic!("no done message in: {out}"));
    assert_eq!(done["ok"], true, "{done}");
    assert_eq!(done["has_results"], true);
    assert!(dir.join("results.h5").is_file());
}

// ------------------------------------------------------------------------------------------------------------------------ scripts (S4, ADR 93)

/// The script worker is built with the workspace and sits next to `nyx`; a run of `-p nyx` alone may not have built it.
fn need_worker() {
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_nyx"));
    let w = exe.parent().unwrap().join(format!("nyx-script-worker{}", std::env::consts::EXE_SUFFIX));
    assert!(w.is_file(), "{} is missing: build it with `cargo build --release -p nyx-script` (a workspace test run builds it)", w.display());
}

#[test]
fn a_script_is_checked_against_the_commands_before_it_runs() {
    need_worker();
    let dir = tmp("script-check");
    let good = dir.join("good.nys");
    std::fs::write(&good, "nyx.project_new(dim=2, title='t')\nnyx.set_param(name='w', value=5e-7)\n").unwrap();
    let o = nyx(&["check", good.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(0), "{}", text(&o.stderr));
    assert!(text(&o.stdout).contains("nothing was run"));

    let bad = dir.join("bad.nys");
    std::fs::write(&bad, "nyx.project_new(dim=2)\nnyx.solve_runn(physics='ems')\n").unwrap();
    let o = nyx(&["check", bad.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2));
    let err = text(&o.stderr);
    assert!(err.contains("solve_runn") && err.contains("bad.nys:2"), "{err}");
    // `run` checks first: the script fails the static check (exit 2) and nothing ran
    let o = nyx(&["run", bad.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(2), "{}", text(&o.stderr));
}

#[test]
fn a_script_runs_prints_returns_and_is_stopped_by_its_limits() {
    need_worker();
    let dir = tmp("script-run");
    let f = dir.join("ok.nys");
    std::fs::write(&f, "nyx.project_new(dim=3, title='t')\nprint('hello from a script')\nnyx.project_info()['dim']\n").unwrap();
    let o = nyx(&["run", f.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(0), "{}", text(&o.stderr));
    let out = text(&o.stdout);
    assert!(out.contains("hello from a script") && out.contains("3") && out.contains("2 engine commands"), "{out}");
    let j = nyx(&["run", f.to_str().unwrap(), "--json"]);
    let v: serde_json::Value = serde_json::from_str(text(&j.stdout).lines().last().unwrap()).expect("a JSON result");
    assert_eq!(v["value"], 3);
    assert_eq!(v["calls"], 2);

    let spin = dir.join("spin.nys");
    std::fs::write(&spin, "while True:\n    pass\n").unwrap();
    let o = nyx(&["run", spin.to_str().unwrap(), "--seconds", "1"]);
    assert_eq!(o.status.code(), Some(1));
    // a script cannot write outside its folder: the guard answers with the command and the argument
    let esc = dir.join("escape.nys");
    std::fs::write(&esc, "nyx.project_new(dim=2)\nnyx.project_save(path='../escaped.nyx.json')\n").unwrap();
    let o = nyx(&["run", esc.to_str().unwrap()]);
    assert_eq!(o.status.code(), Some(1));
    assert!(text(&o.stderr).contains("project_save"), "{}", text(&o.stderr));
    assert!(!dir.parent().unwrap().join("escaped.nyx.json").exists());
}

#[test]
fn api_stubs_prints_the_type_stub_of_the_nyx_object() {
    let o = nyx(&["api", "stubs"]);
    assert_eq!(o.status.code(), Some(0), "{}", text(&o.stderr));
    let out = text(&o.stdout);
    assert!(out.contains("class _Nyx:") && out.contains("def solve_run(self") && out.contains("nyx: _Nyx"), "{}", &out[..out.len().min(400)]);
}
