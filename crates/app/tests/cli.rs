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
