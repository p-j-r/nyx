// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Command lists (productivity plan S1): a JSON-lines list is checked before it runs (every command in the catalog, every argument against its
//! schema, then a dry run that finds what only the project can tell), the check starts no job and changes nothing, a bad line is named with its
//! command, and the workflows of `tools/workflows/` run from their files.

use nyx_engine::commandlist::{parse, run, Report};
use nyx_engine::{Engine, MeshSettings, MeshedModel, Mesher};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

struct Counting {
    inner: nyx_engine::FemcoreMesher,
    calls: AtomicUsize,
}

impl Mesher for Counting {
    fn mesh(
        &self,
        geometry: &geom::model::FeatureList,
        settings: &MeshSettings,
        sizing: Option<&nyx_physics::Sizing>,
    ) -> Result<MeshedModel, Error> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        self.inner.mesh(geometry, settings, sizing)
    }
}

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The text of a list in `tools/workflows/`, with the example paths made absolute so the test does not depend on the working directory.
fn list(name: &str) -> String {
    let text = std::fs::read_to_string(root().join("tools/workflows").join(name)).expect("a workflow list");
    let abs = root().canonicalize().unwrap().to_string_lossy().replace('\\', "/");
    let abs = abs.strip_prefix("//?/").unwrap_or(&abs).to_string();
    text.replace("\"examples/", &format!("\"{abs}/examples/")).replace("\"runs/", &format!("\"{abs}/runs/"))
}

fn counting_engine() -> (Engine, Arc<Counting>) {
    let mut e = Engine::new();
    let m = Arc::new(Counting { inner: nyx_engine::FemcoreMesher, calls: AtomicUsize::new(0) });
    e.set_mesher(m.clone());
    (e, m)
}

fn example(name: &str) -> String {
    root()
        .canonicalize()
        .unwrap()
        .join("examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .replace('\\', "/")
        .trim_start_matches("//?/")
        .to_string()
}

/// A list of lines as text.
fn text(lines: &[Value]) -> String {
    lines.iter().map(|l| l.to_string()).collect::<Vec<_>>().join("\n")
}

fn problems_of(r: &Report) -> Vec<String> {
    r.problems.iter().map(|p| p.to_string()).collect()
}

#[test]
fn every_workflow_list_checks_clean_and_the_check_starts_nothing() {
    let (e, mesher) = counting_engine();
    let stamp = e.stamp();
    let mut files: Vec<String> = std::fs::read_dir(root().join("tools/workflows"))
        .unwrap()
        .map(|f| f.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|f| f.ends_with(".jsonl"))
        .collect();
    files.sort();
    assert!(files.len() >= 7, "{files:?}");
    for f in &files {
        let r = e.check_list(&list(f));
        assert!(r.ok(), "{f}: {:?}", problems_of(&r));
        assert!(r.commands > 0 && r.hash.len() >= 8);
    }
    assert!(e.jobs().list().is_empty(), "a check started a job");
    assert_eq!(mesher.calls.load(Ordering::SeqCst), 0, "a check meshed");
    assert_eq!(e.stamp(), stamp, "a check changed the engine's project");
}

#[test]
fn the_plan_counts_the_meshes_and_solves_a_list_would_make() {
    let (e, _) = counting_engine();
    // w1: the width set five times with a solve each, and a first solve: six meshes, six solves
    let r = e.check_list(&list("w1_soi_width.jsonl"));
    assert!(r.ok(), "{:?}", problems_of(&r));
    assert_eq!((r.meshes, r.solves), (6, 6), "{:?}", r.steps);
    // two solves of one project build one mesh
    let wire = example("15_soi_wire");
    let two = text(&[
        json!({ "cmd": "project_open", "args": { "path": wire } }),
        json!({ "cmd": "solve_run" }),
        json!({ "cmd": "solve_run" }),
    ]);
    let r = e.check_list(&two);
    assert_eq!((r.meshes, r.solves), (1, 2));
    assert!(r.steps[2].summary.contains("reused"), "{:?}", r.steps);
    // a convergence check right after a solve of that physics alone reuses its solve; the finer mesh is new
    let conv = text(&[
        json!({ "cmd": "project_open", "args": { "path": wire } }),
        json!({ "cmd": "solve_run", "args": { "physics": "ems" } }),
        json!({ "cmd": "convergence_check", "args": { "physics": "ems" } }),
    ]);
    let r = e.check_list(&conv);
    assert!(r.ok(), "{:?}", problems_of(&r));
    assert_eq!((r.meshes, r.solves), (2, 2), "{:?}", r.steps);
    // a sweep counts its points
    let sweep = text(&[
        json!({ "cmd": "project_open", "args": { "path": wire } }),
        json!({ "cmd": "sweep_run", "args": { "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": [4.5e-7, 5.5e-7] } } }),
    ]);
    let r = e.check_list(&sweep);
    assert!(r.ok(), "{:?}", problems_of(&r));
    assert_eq!((r.meshes, r.solves), (2, 2), "{:?}", r.steps);
}

#[test]
fn a_bad_line_is_named_with_its_command_and_nothing_runs() {
    let (mut e, mesher) = counting_engine();
    let wire = example("15_soi_wire");
    let cases: Vec<(Value, &str)> = vec![
        (json!({ "cmd": "set_param", "args": { "name": "w" } }), "missing 'value'"),
        (json!({ "cmd": "set_param", "args": { "name": "w", "value": 1e-7, "unit": "m" } }), "no argument 'unit'"),
        (json!({ "cmd": "set_param", "args": { "name": 5, "value": 1e-7 } }), "args.name: must be string"),
        (json!({ "cmd": "set_param", "args": { "name": "w", "value": true } }), "args.value"),
        (json!({ "cmd": "solve_runn" }), "did you mean solve_run"),
        (json!({ "cmd": "set_mesh", "args": { "h": "big" } }), "args.h"),
        (json!({ "cmd": "describe_physics", "args": { "physics": "xyz" } }), "is not one of"),
        (json!({ "cmd": "remove_condition", "args": { "physics": "ems", "index": 1.5 } }), "args.index"),
    ];
    for (bad, word) in cases {
        let t = text(&[
            json!({ "cmd": "project_open", "args": { "path": wire } }),
            json!({ "cmd": "solve_run", "note": "line 2" }),
            bad.clone(),
        ]);
        let r = e.check_list(&t);
        assert_eq!(r.problems.len(), 1, "{bad}: {:?}", problems_of(&r));
        let p = &r.problems[0];
        assert_eq!(p.line, 3, "{bad}");
        assert_eq!(p.cmd, bad["cmd"].as_str().unwrap());
        assert!(p.message.contains(word), "{bad}: {p}");
        assert!(r.steps.is_empty(), "a list with a bad line was planned");
    }
    // `run` refuses the same list before it runs anything: the parse and the static check come first
    let (lines, parse_problems) = parse(&text(&[json!({ "cmd": "solve_run" }), json!({ "cmd": "nope" })]));
    assert!(parse_problems.is_empty());
    assert_eq!(e.check_static(&lines).len(), 1);
    assert!(e.jobs().list().is_empty());
    assert_eq!(mesher.calls.load(Ordering::SeqCst), 0);
    let _ = &mut e;
}

#[test]
fn lines_that_are_not_commands_are_all_reported_with_their_line_numbers() {
    let t = "# a comment\n\n{\"cmd\":\"project_new\"}\nnot json\n[1,2]\n{\"args\":{}}\n{\"cmd\":\"set_param\",\"extra\":1}\n{\"cmd\":\"set_param\",\"args\":3}\n";
    let (lines, problems) = parse(t);
    assert_eq!(lines.len(), 1);
    let at: Vec<usize> = problems.iter().map(|p| p.line).collect();
    assert_eq!(at, vec![4, 5, 6, 7, 8], "{problems:?}");
    assert!(problems[0].message.contains("not JSON"));
    assert!(problems[2].message.contains("no `cmd`"));
    assert!(problems[3].message.contains("unknown key 'extra'") && problems[3].cmd == "set_param");
    let (e, _) = counting_engine();
    assert_eq!(e.check_list(t).problems.len(), 5);
}

#[test]
fn the_plan_finds_what_only_the_project_can_tell() {
    let (e, mesher) = counting_engine();
    let wire = example("15_soi_wire");
    let open = json!({ "cmd": "project_open", "args": { "path": wire } });
    // a solve with no physics; a feature that does not exist; a result nothing solves; a problem that does not check
    let cases: Vec<(Vec<Value>, usize, &str)> = vec![
        (vec![json!({ "cmd": "project_new" }), json!({ "cmd": "solve_run" })], 2, "no physics"),
        (vec![open.clone(), json!({ "cmd": "remove_feature", "args": { "name": "nothing_like_this" } })], 2, ""),
        (
            vec![open.clone(), json!({ "cmd": "result_info", "args": { "physics": "ems" } })],
            2,
            "nothing before this line solves it",
        ),
        (vec![open.clone(), json!({ "cmd": "solve_run", "args": { "physics": "dfs" } })], 2, "no problem for 'dfs'"),
        (
            vec![
                open.clone(),
                json!({ "cmd": "solve_run" }),
                open.clone(),
                json!({ "cmd": "result_info", "args": { "physics": "ems" } }),
            ],
            4,
            "nothing before this line solves it",
        ),
    ];
    for (lines, line, word) in cases {
        let r = e.check_list(&text(&lines));
        assert_eq!(r.problems.len(), 1, "{lines:?}: {:?}", problems_of(&r));
        assert_eq!(r.problems[0].line, line, "{:?}", problems_of(&r));
        assert!(r.problems[0].message.contains(word), "{}", r.problems[0]);
    }
    // a result read after the solve that makes it is fine
    let ok = e.check_list(&text(&[
        open.clone(),
        json!({ "cmd": "solve_run" }),
        json!({ "cmd": "result_info", "args": { "physics": "ems" } }),
    ]));
    assert!(ok.ok(), "{:?}", problems_of(&ok));
    assert_eq!(mesher.calls.load(Ordering::SeqCst), 0);
}

#[test]
fn a_job_arg_line_needs_no_id_to_pass_the_check() {
    let (e, _) = counting_engine();
    let r = e.check_list(&list("w3a_cancel_early.jsonl"));
    assert!(r.ok(), "{:?}", problems_of(&r));
}

#[test]
fn a_run_replays_the_lists_the_same_way_the_engine_would() {
    // w1 from its file: open, six solves, a read after each; the observer sees every command, and the numbers are those of a direct call
    let (mut e, _) = counting_engine();
    let (lines, bad) = parse(&list("w1_soi_width.jsonl"));
    assert!(bad.is_empty());
    let mut seen = vec![];
    run(&mut e, &lines, &mut |d| seen.push((d.line.number, d.line.cmd.clone(), d.reply.as_ref().ok().cloned())))
        .unwrap();
    assert_eq!(seen.len(), lines.len());
    assert_eq!(seen.iter().filter(|s| s.1 == "solve_run").count(), 6);
    let last = seen.iter().rev().find(|s| s.1 == "result_info").unwrap().2.clone().unwrap();
    let mut direct = Engine::new();
    direct.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
    direct.call("set_param", &json!({ "name": "w", "value": 6.5e-7 })).unwrap();
    direct.call("solve_run", &json!({})).unwrap();
    // the numbers: the eigenvalues and the scalars (the reply also carries timings, which differ)
    let numbers = |v: &Value| {
        let eig: Vec<Value> = v[0]["solutions"].as_array().unwrap().iter().map(|s| s["eigenvalue"].clone()).collect();
        (eig, v[0]["scalars"].clone(), v[0]["dofs"].clone())
    };
    let want = direct.call("result_info", &json!({ "physics": "ems" })).unwrap();
    assert_eq!(numbers(&last), numbers(&want));
}

#[test]
fn a_run_stops_at_the_first_failing_line_and_names_it() {
    let (mut e, _) = counting_engine();
    let wire = example("15_soi_wire");
    let (lines, _) = parse(&text(&[
        json!({ "cmd": "project_open", "args": { "path": wire } }),
        json!({ "cmd": "remove_feature", "args": { "name": "nothing_like_this" } }),
        json!({ "cmd": "solve_run" }),
    ]));
    let mut n = 0;
    let err = run(&mut e, &lines, &mut |_| n += 1).unwrap_err();
    assert_eq!((err.line, err.cmd.as_str()), (2, "remove_feature"), "{err}");
    assert_eq!(n, 2, "the failing line is reported to the observer, the next one is not run");
    assert!(e.jobs().list().is_empty(), "the solve after the failure ran");
    // an error that the line expects is the outcome; a success where an error was expected stops the run
    let (lines, _) = parse(&text(&[
        json!({ "cmd": "project_open", "args": { "path": wire } }),
        json!({ "cmd": "remove_feature", "args": { "name": "nothing_like_this" }, "expect_error": true }),
        json!({ "cmd": "project_info", "expect_error": true }),
    ]));
    let err = run(&mut e, &lines, &mut |_| {}).unwrap_err();
    assert_eq!(err.line, 3, "{err}");
    assert!(err.message.contains("expect_error"));
}

#[test]
fn the_cancel_workflow_runs_from_its_file_and_a_cancel_that_does_not_happen_fails() {
    let (mut e, _) = counting_engine();
    let (lines, bad) = parse(&list("w3a_cancel_early.jsonl"));
    assert!(bad.is_empty());
    let mut since = None;
    run(&mut e, &lines, &mut |d| {
        if d.since_cancel.is_some() {
            since = d.since_cancel;
        }
    })
    .unwrap();
    assert!(since.is_some_and(|s| s >= 0.0), "the run did not measure the cancel");
    // a cancel with no job to cancel is an error of the run, not a pass
    let (mut fresh, _) = counting_engine();
    let (lines, _) = parse(&text(&[json!({ "cmd": "job_cancel", "job_arg": "id" })]));
    assert!(run(&mut fresh, &lines, &mut |_| {}).unwrap_err().message.contains("no job was started"));
}
