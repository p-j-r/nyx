// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Drives the real binary over stdio the way an MCP host does (ADR 94: four script tools). The tests that run a script need the script worker
//! program next to `nyx-mcp`: `cargo build --release -p nyx-script` (a workspace build has it).
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

struct Host {
    child: Child,
    tx: ChildStdin,
    rx: BufReader<ChildStdout>,
    next: u64,
    dir: PathBuf,
}

fn worker_present() -> bool {
    let w = PathBuf::from(env!("CARGO_BIN_EXE_nyx-mcp")).with_file_name(format!("nyx-script-worker{}", std::env::consts::EXE_SUFFIX));
    assert!(w.is_file(), "the script worker is not at {}: build it with `cargo build --release -p nyx-script`", w.display());
    true
}

impl Host {
    fn start() -> Host {
        static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!("nyx_mcp_test_{}_{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
        let mut child = Command::new(env!("CARGO_BIN_EXE_nyx-mcp"))
            .arg("--dir")
            .arg(&dir)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn nyx-mcp");
        let tx = child.stdin.take().unwrap();
        let rx = BufReader::new(child.stdout.take().unwrap());
        Host { child, tx, rx, next: 1, dir }
    }
    fn send(&mut self, v: Value) {
        writeln!(self.tx, "{v}").unwrap();
        self.tx.flush().unwrap();
    }
    fn request(&mut self, method: &str, params: Value) -> Value {
        let id = self.next;
        self.next += 1;
        self.send(json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params }));
        let mut line = String::new();
        self.rx.read_line(&mut line).unwrap();
        let v: Value = serde_json::from_str(&line).unwrap_or_else(|e| panic!("bad reply {line:?}: {e}"));
        assert_eq!(v["id"], json!(id));
        v
    }
    fn tool(&mut self, name: &str, args: Value) -> Value {
        self.request("tools/call", json!({ "name": name, "arguments": args }))
    }
    /// Sends a request without waiting for its answer: its id.
    fn send_request(&mut self, method: &str, params: Value) -> u64 {
        let id = self.next;
        self.next += 1;
        self.send(json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params }));
        id
    }
    /// The next line of the server, whatever it is.
    fn read(&mut self) -> Value {
        let mut line = String::new();
        self.rx.read_line(&mut line).unwrap();
        serde_json::from_str(&line).unwrap_or_else(|e| panic!("bad line {line:?}: {e}"))
    }
    fn init(&mut self, proto: &str) -> Value {
        let r = self.request("initialize", json!({ "protocolVersion": proto, "capabilities": {}, "clientInfo": { "name": "test", "version": "0" } }));
        self.send(json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }));
        r
    }
    /// A script that must succeed: the answer's JSON.
    fn script(&mut self, source: &str) -> Value {
        ok(&self.tool("script_run", json!({ "source": source })))
    }
}

impl Drop for Host {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn payload(r: &Value) -> Value {
    serde_json::from_str(r["result"]["content"][0]["text"].as_str().unwrap()).unwrap()
}

fn ok(r: &Value) -> Value {
    assert_eq!(r["result"]["isError"], false, "{r}");
    payload(r)
}

fn failure(r: &Value) -> String {
    assert_eq!(r["result"]["isError"], true, "{r}");
    payload(r)["message"].as_str().unwrap().to_string()
}

#[test]
fn handshake_negotiates_the_protocol_and_lists_exactly_the_four_script_tools() {
    let mut h = Host::start();
    let r = h.init("2024-11-05");
    assert_eq!(r["result"]["protocolVersion"], "2024-11-05");
    assert_eq!(r["result"]["serverInfo"]["name"], "nyx-mcp");
    let instructions = r["result"]["instructions"].as_str().unwrap();
    assert!(instructions.contains("script_run") && instructions.contains("api"));
    assert_eq!(h.request("ping", json!({}))["result"], json!({}));
    let l = h.request("tools/list", json!({}));
    let mut names: Vec<&str> = l["result"]["tools"].as_array().unwrap().iter().map(|t| t["name"].as_str().unwrap()).collect();
    names.sort();
    assert_eq!(names, vec!["api", "script_check", "script_run", "script_wait"], "the one-tool-per-command catalog is gone (Q9)");
    for t in l["result"]["tools"].as_array().unwrap() {
        assert_eq!(t["inputSchema"]["type"], "object", "{}", t["name"]);
    }
}

#[test]
fn unknown_version_falls_back_and_unknown_method_and_tool_are_errors() {
    let mut h = Host::start();
    assert_eq!(h.init("1999-01-01")["result"]["protocolVersion"], "2025-06-18");
    assert_eq!(h.request("nope/nothing", json!({}))["error"]["code"], -32601);
    // the old per-command tools do not exist any more
    assert_eq!(h.tool("add_rect", json!({}))["error"]["code"], -32602);
    assert_eq!(h.tool("get_guide", json!({}))["error"]["code"], -32602);
}

#[test]
fn api_is_the_guide_the_command_list_and_one_command_in_plain_text() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("api", json!({}));
    assert_eq!(r["result"]["isError"], false);
    assert!(r["result"].get("structuredContent").is_none(), "an answer is sent once, as text");
    let guide = r["result"]["content"][0]["text"].as_str().unwrap();
    assert!(guide.contains("script_run") && guide.contains("nyx.set_param"), "{guide}");
    let list = h.tool("api", json!({ "topic": "commands" }));
    let list = list["result"]["content"][0]["text"].as_str().unwrap();
    assert!(list.contains("add_rect") && list.contains("solve_run") && list.contains("integrate"), "{list}");
    let one = h.tool("api", json!({ "command": "add_rect" }));
    let one = one["result"]["content"][0]["text"].as_str().unwrap();
    assert!(one.contains("def add_rect(") && one.contains("min"), "the signature and the arguments: {one}");
    let geometry = h.tool("api", json!({ "topic": "geometry" }));
    assert!(geometry["result"]["content"][0]["text"].as_str().unwrap().contains("add_overlay"));
    // bad arguments are a tool error that names the valid values
    let m = failure(&h.tool("api", json!({ "topic": "bogus" })));
    assert!(m.contains("optical-modes"), "{m}");
    let m = failure(&h.tool("api", json!({ "command": "no_such" })));
    assert!(m.contains("no_such"), "{m}");
}

#[test]
fn a_script_runs_prints_and_returns_one_answer_with_the_project_header() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("script_run", json!({ "allow_destructive": true, "source": "nyx.project_new(dim=2, title='t')\nprint('hello', 1 + 1)\nnyx.project_info()['title']" }));
    assert!(r["result"].get("structuredContent").is_none());
    let a = ok(&r);
    assert_eq!(a["status"], "done", "{a}");
    assert_eq!(a["value"], "t");
    assert_eq!(a["printed"], json!(["hello 2"]));
    assert_eq!(a["calls"], 2);
    assert!(a["nyx"]["session"].is_string() && a["nyx"]["project_hash"].is_string() && a["nyx"]["revision"].as_u64().unwrap() >= 1, "{a}");
    // a script error is a tool error with the line
    let m = failure(&h.tool("script_run", json!({ "source": "a = 1\nnyx.set_param(name='w', value=None)\n" })));
    assert!(m.contains("line 2") && m.contains("set_param"), "{m}");
    // a misspelled command is refused before anything runs
    let m = failure(&h.tool("script_run", json!({ "source": "nyx.set_prm(name='w', value=1)" })));
    assert!(m.contains("set_prm"), "{m}");
    // and the session goes on
    assert_eq!(h.script("1 + 1")["value"], 2);
}

#[test]
fn variables_stay_and_a_saved_script_runs_again_with_other_params() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    h.script("x = 41\ndef f(a):\n    return a + x");
    assert_eq!(h.script("f(1)")["value"], 42);
    // the check sees the session's variables, and runs nothing
    assert_eq!(ok(&h.tool("script_check", json!({ "source": "f(x)" })))["status"], "ok");
    let m = failure(&h.tool("script_check", json!({ "source": "nyx.set_param(nam='w', value=1)" })));
    assert!(m.contains("nam") || m.contains("name"), "{m}");
    // saved once, run again with other values and nothing resent
    let src = "nyx.project_new(dim=2, title='saved')\nnyx.set_param(name='w', value=params.get('w', 5e-7))\nnyx.project_info()['params'][0]['value']";
    assert_eq!(ok(&h.tool("script_run", json!({ "source": src, "save_as": "wire", "allow_destructive": true })))["value"], 5e-7);
    assert!(h.dir.join("scripts").join("wire.nys").is_file());
    let a = ok(&h.tool("script_run", json!({ "path": "wire", "params": { "w": 3e-7 }, "allow_destructive": true })));
    assert_eq!(a["value"], 3e-7, "{a}");
    assert_eq!(ok(&h.tool("script_check", json!({ "path": "wire" })))["status"], "ok");
    // names that reach outside the scripts folder, or no such script, are refused
    for bad in ["../wire", "a/b", "", "wire.py"] {
        let m = failure(&h.tool("script_run", json!({ "path": bad })));
        assert!(m.contains("script name") || m.contains("no saved script") || m.contains("`source`"), "{bad}: {m}");
    }
    let m = failure(&h.tool("script_run", json!({ "path": "nothere" })));
    assert!(m.contains("no saved script"), "{m}");
    let m = failure(&h.tool("script_run", json!({ "source": "1", "path": "wire" })));
    assert!(m.contains("not both"), "{m}");
}

#[test]
fn the_commands_that_discard_the_project_need_the_flag() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let m = failure(&h.tool("script_run", json!({ "source": "nyx.project_new(dim=2)" })));
    assert!(m.contains("allow_destructive"), "{m}");
    assert_eq!(ok(&h.tool("script_run", json!({ "source": "nyx.project_new(dim=2)\n1", "allow_destructive": true })))["status"], "done");
}

#[test]
fn a_script_cannot_name_a_file_outside_the_root() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let outside = std::env::temp_dir().join("nyx_mcp_outside.nyx.json");
    let m = failure(&h.tool("script_run", json!({ "source": format!("nyx.project_open(path={:?})", outside.to_string_lossy()), "allow_destructive": true })));
    assert!(m.to_lowercase().contains("outside") || m.to_lowercase().contains("folder") || m.to_lowercase().contains("not allowed"), "{m}");
}

#[test]
fn a_long_value_is_cut_and_says_how_to_print_instead() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let a = h.script("list(range(200000))");
    assert_eq!(a["value"]["truncated"], true, "{a}");
    assert!(a["value"]["bytes"].as_u64().unwrap() > 48_000);
}

const MODEL: &str = "
nyx.project_new(title='wire from mcp', h=1.2e-7)
nyx.add_rect(name='box', min=[-2.5e-6, -2e-6], max=[2.5e-6, 0.0])
nyx.add_rect(name='clad', min=[-2.5e-6, 0.0], max=[2.5e-6, 1.5e-6])
nyx.add_rect(name='core', min=[-2.5e-7, 0.0], max=[2.5e-7, 3e-7])
nyx.add_overlay(name='regions', regions=[{'name': 'BOX', 'input': 'box'}, {'name': 'clad', 'input': 'clad'}, {'name': 'Si', 'input': 'core'}])
nyx.add_physics(physics='ems', study='modes', order=2, study_params={'freq': 299792458.0 / 1.55e-6, 'n_target': 2.9, 'n_modes': 2})
nyx.set_material(physics='ems', region='BOX', material='silica')
nyx.set_material(physics='ems', region='clad', material='silica')
nyx.set_material(physics='ems', region='Si', material='silicon')
nyx.set_mesh(adapt_scale=0.4)
";

#[test]
fn a_host_builds_a_wire_solves_it_in_one_script_and_the_work_is_a_stored_workflow() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let src = format!(
        "{MODEL}\nsolved = nyx.solve_run()\ninfo = nyx.result_info(physics='ems')\nneff = info[0]['solutions'][0]['eigenvalue'][0]\nprint('neff', neff)\npol = nyx.analyse(analysis='polarization', physics='ems')\n{{'neff': neff, 'class': pol['rows'][0]['report']['class'], 'run': solved['run']['id'], 'runs': len(nyx.store_runs()), 'folder': nyx.project_info()['store']['folder']}}"
    );
    let a = ok(&h.tool("script_run", json!({ "source": src, "allow_destructive": true })));
    assert_eq!(a["status"], "done", "{a}");
    let neff = a["value"]["neff"].as_f64().unwrap();
    assert!((neff - 2.67).abs() < 0.03, "TE0 effective index {neff}");
    assert!(a["value"]["class"].as_str().unwrap().contains("TE"), "{a}");
    assert!(a["value"]["run"].is_string(), "the run is recorded in the workflow's run table: {a}");
    assert_eq!(a["value"]["runs"], 1);
    assert_eq!(a["printed"].as_array().unwrap().len(), 1);
    // the host did not name a folder: the work sits in the server's workspace as a workflow with a run, and the script that made it is kept
    let folder = a["value"]["folder"].as_str().expect("bound to a workflow folder");
    assert!(std::path::Path::new(folder).join("workflow.json").exists(), "{folder}");
}

#[test]
fn a_host_finds_the_axisymmetric_solid_study_and_gets_a_quality_factor_for_a_disk_on_a_post() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    // a disk (r 0 to 6 um) on a post on a substrate, layers to the right and below, as the engine example 17 (x is r); the study is found in the
    // description of sms with its three parameters; the convergence check pairs the modes of the two meshes, or follows one (ADR 80)
    let src = r#"
um = 1e-6
xw = 16 * um
zw = 8 * um
t = 5 * um
d = nyx.describe_physics(physics='sms')
names = []
for s in d['studies']:
    if s['kind'] == 'axisymmetric_eigenfrequency':
        for p in s['params']:
            names.append(p['name'])
nyx.project_new(title='disk on a post', dim=2, h=um)
nyx.add_rect(name='substrate_box', min=[0.0, -zw - 0.5 * um - t], max=[xw + t, -0.5 * um])
nyx.add_rect(name='post_box', min=[0.0, -0.5 * um], max=[0.8 * um, 0.0])
nyx.add_rect(name='disk_box', min=[0.0, 0.0], max=[6.0 * um, um])
nyx.add_rect(name='window', min=[0.0, -zw - 0.5 * um], max=[xw, 3.0 * um])
nyx.add_overlay(name='regions', regions=[{'name': 'substrate', 'input': 'substrate_box'}, {'name': 'post', 'input': 'post_box'}, {'name': 'disk', 'input': 'disk_box'}])
nyx.add_pml(name='pml', window='window', sides=['XMax', 'ZMin'], thickness=t, layers=3.0, grading=2.0, regions=['substrate'])
nyx.add_size(name='size_disk', sel='region(disk) | region(post)', h=0.4 * um)
nyx.add_physics(physics='sms', study='axisymmetric_eigenfrequency', order=2, study_params={'n_modes': 8, 'f_search': 4.6e8, 'azimuthal_order': 0})
steel = {'E': 160e9, 'nu': 0.22, 'rho': 2330.0}
nyx.set_material(physics='sms', region='substrate*', material=steel)
nyx.set_material(physics='sms', region='post', material=steel)
nyx.set_material(physics='sms', region='disk', material=steel)
nyx.add_pml_conditions(physics='sms', strength=3.0)
fixed = {'ux': 0.0, 'uy': 0.0, 'uz': 0.0}
nyx.add_condition(physics='sms', on='outer & plane(x = ' + str(xw + t) + ')', kind='fixed', params=fixed)
nyx.add_condition(physics='sms', on='outer & plane(z = ' + str(-zw - 0.5 * um - t) + ')', kind='fixed', params=fixed)
solved = nyx.solve_run()
res = nyx.result_info(physics='sms')[0]
near = nyx.convergence_check(physics='sms', factor=0.7, match='nearest')
track = nyx.convergence_check(physics='sms', factor=0.7, track={'f': 4.6e8})
bad = ''
try:
    nyx.convergence_check(physics='sms', match='closest')
except RuntimeError as e:
    bad = str(e)
{'names': names, 'solved': solved['results'][0]['physics'], 'study': res['study'], 'scalars': str(res['scalars']),
 'eigen_parts': len(res['solutions'][0]['eigenvalue']), 'near_match': near['match'], 'pairs': len(near['pairs']),
 'unpaired_listed': near['unpaired']['coarse'] is not None, 'track_pairs': len(track['pairs']), 'bad': bad}
"#;
    let a = ok(&h.tool("script_run", json!({ "source": src, "allow_destructive": true, "wait_s": 300 })));
    assert_eq!(a["status"], "done", "{a}");
    let v = &a["value"];
    assert_eq!(v["names"], json!(["n_modes", "f_search", "azimuthal_order"]), "{v}");
    assert_eq!(v["solved"], "sms", "{v}");
    assert_eq!(v["study"], "axisymmetric_eigenfrequency", "{v}");
    let scalars = v["scalars"].as_str().unwrap();
    for name in ["Q1", "layer fraction 1", "f1 (Hz)", "f re 1 (Hz)", "f im 1 (Hz)", "radial share 1", "share in disk 1"] {
        assert!(scalars.contains(name), "the result carries {name}: {scalars}");
    }
    assert_eq!(v["eigen_parts"], 2, "{v}");
    assert_eq!(v["near_match"], "nearest", "{v}");
    assert!(v["pairs"].as_u64().unwrap() > 0 && v["unpaired_listed"] == true, "{v}");
    assert_eq!(v["track_pairs"], 1, "{v}");
    assert!(!v["bad"].as_str().unwrap().is_empty(), "a bad `match` is an engine error the script can catch: {v}");
}

#[test]
fn a_solve_that_outlasts_wait_s_answers_running_and_script_wait_finishes_it() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let a = ok(&h.tool("script_run", json!({ "source": format!("{MODEL}\nprint('before')\nr = nyx.solve_run()\nprint('after')\n'finished'"), "allow_destructive": true, "wait_s": 0 })));
    assert_eq!(a["status"], "running", "{a}");
    assert!(a["waiting_for"].as_str().unwrap().contains("solve_run"), "{a}");
    assert_eq!(a["printed"], json!(["before"]));
    // another script cannot start while this one waits
    let m = failure(&h.tool("script_run", json!({ "source": "1" })));
    assert!(m.contains("script_wait"), "{m}");
    let b = ok(&h.tool("script_wait", json!({ "wait_s": 120 })));
    assert_eq!(b["status"], "done", "{b}");
    assert_eq!(b["value"], "finished");
    assert_eq!(b["printed"], json!(["after"]));
    // the variables of the finished script are there; nothing waits now
    assert!(h.script("r")["value"].is_object());
    let m = failure(&h.tool("script_wait", json!({})));
    assert!(m.contains("no script is running"), "{m}");
}

#[test]
fn a_wait_longer_than_the_server_allows_is_cut_and_the_answer_says_so() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let a = ok(&h.tool("script_run", json!({ "source": "1 + 1", "wait_s": 600 })));
    assert_eq!(a["status"], "done", "{a}");
    assert_eq!(a["wait_s"]["asked"], 600.0, "{a}");
    assert_eq!(a["wait_s"]["used"], 90.0, "{a}");
    // within the limit nothing is said
    assert!(ok(&h.tool("script_run", json!({ "source": "2", "wait_s": 30 })))["wait_s"].is_null());
    // the limit is in the tool description and the schema, and a negative wait is refused
    let l = h.request("tools/list", json!({}));
    let run = l["result"]["tools"].as_array().unwrap().iter().find(|t| t["name"] == "script_run").unwrap().clone();
    assert_eq!(run["inputSchema"]["properties"]["wait_s"]["maximum"], 90.0, "{run}");
    let m = failure(&h.tool("script_run", json!({ "source": "3", "wait_s": -1 })));
    assert!(m.contains("wait_s"), "{m}");
}

/// A model whose solve takes long enough for a host to send three more messages while it runs.
const HEAVY: &str = "
nyx.set_mesh(h=4e-8, adapt_scale=0.4)
";

#[test]
fn a_waiting_call_ends_when_the_host_sends_a_newer_one_or_cancels_and_ping_is_answered_meanwhile() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let run = h.send_request("tools/call", json!({ "name": "script_run", "arguments": { "source": format!("{MODEL}{HEAVY}print('before')\nnyx.solve_run()\n'finished'"), "allow_destructive": true, "wait_s": 60 } }));
    std::thread::sleep(std::time::Duration::from_millis(1500));
    // a ping is answered at once, the long wait goes on (the old server answered it after the wait)
    let ping = h.send_request("ping", json!({}));
    let r = h.read();
    assert_eq!(r["id"], json!(ping), "the ping was answered during the wait: {r}");
    assert_eq!(r["result"], json!({}));
    // the host gave up on the first call and asks again: the first is answered `running` now, with what it printed
    let wait = h.send_request("tools/call", json!({ "name": "script_wait", "arguments": { "wait_s": 60 } }));
    let r = h.read();
    assert_eq!(r["id"], json!(run), "{r}");
    let a = ok(&r);
    assert_eq!(a["status"], "running", "the solve outlasts this test; make HEAVY heavier if this is `done`: {a}");
    assert_eq!(a["printed"], json!(["before"]), "{a}");
    assert!(a["note"].as_str().unwrap().contains("newer request"), "{a}");
    // the host cancels the second: it gets no answer, and the request after it is served at once
    h.send(json!({ "jsonrpc": "2.0", "method": "notifications/cancelled", "params": { "requestId": wait, "reason": "timeout" } }));
    let cancel = h.send_request("tools/call", json!({ "name": "script_wait", "arguments": { "cancel": true } }));
    let r = h.read();
    assert_eq!(r["id"], json!(cancel), "no answer is sent for a cancelled request: {r}");
    assert_eq!(ok(&r)["status"], "cancelled");
    // the server is still usable
    assert_eq!(h.script("1 + 1")["value"], 2);
}

#[test]
fn a_call_with_a_progress_token_gets_progress_notifications_while_it_waits() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    let id = h.send_request("tools/call", json!({ "name": "script_run", "_meta": { "progressToken": "t1" }, "arguments": { "source": format!("{MODEL}{HEAVY}nyx.solve_run()"), "allow_destructive": true, "wait_s": 60 } }));
    let n = h.read();
    assert_eq!(n["method"], "notifications/progress", "{n}");
    assert_eq!(n["params"]["progressToken"], "t1");
    assert!(n["params"]["progress"].as_f64().unwrap() >= 4.0, "{n}");
    // the call itself is answered when its wait ends: here the host abandons it
    let cancel = h.send_request("tools/call", json!({ "name": "script_wait", "arguments": { "cancel": true } }));
    let mut r = h.read();
    if r["id"] == json!(id) {
        r = h.read();
    }
    assert_eq!(r["id"], json!(cancel), "{r}");
}

#[test]
fn script_wait_can_abandon_a_waiting_script() {
    worker_present();
    let mut h = Host::start();
    h.init("2025-06-18");
    h.script("kept = 7");
    let a = ok(&h.tool("script_run", json!({ "source": format!("{MODEL}\nnyx.solve_run()"), "allow_destructive": true, "wait_s": 0 })));
    assert_eq!(a["status"], "running", "{a}");
    assert_eq!(ok(&h.tool("script_wait", json!({ "cancel": true })))["status"], "cancelled");
    assert_eq!(h.script("kept")["value"], 7);
}

// ------------------------------------------------------------------------------------------------------------------------ the gate of S7

/// A server whose scripts may name files under `root` (the repository, for the example projects).
fn start_with_root(root: &std::path::Path) -> Host {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(100);
    let dir = std::env::temp_dir().join(format!("nyx_mcp_gate_{}_{}", std::process::id(), N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    let mut child = Command::new(env!("CARGO_BIN_EXE_nyx-mcp"))
        .arg("--dir")
        .arg(&dir)
        .arg("--root")
        .arg(root)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn nyx-mcp");
    let tx = child.stdin.take().unwrap();
    let rx = BufReader::new(child.stdout.take().unwrap());
    Host { child, tx, rx, next: 1, dir }
}

/// What the old server (one tool per command, ADR 57) cost for a command list [assumed]: one call per command; the agent sent the tool name and
/// the arguments; the answer came back twice (as text and as `structuredContent`) with a header of about 85 bytes.
fn per_command_baseline(list: &str, repo: &std::path::Path) -> (usize, usize, usize, Vec<f64>) {
    let mut e = nyx_engine::Engine::new();
    let (mut calls, mut sent, mut returned, mut neffs) = (0, 0, 0, vec![]);
    for line in list.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let l: Value = serde_json::from_str(line).unwrap();
        let cmd = l["cmd"].as_str().unwrap();
        let mut args = l.get("args").cloned().unwrap_or_else(|| json!({}));
        if let Some(p) = args.get("path").and_then(Value::as_str).map(str::to_string) {
            args["path"] = json!(repo.join(p).to_string_lossy());
        }
        let r = e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}: {er}"));
        calls += 1;
        sent += cmd.len() + args.to_string().len();
        returned += 2 * r.to_string().len() + 85;
        if cmd == "result_info" {
            neffs.push(r[0]["solutions"][0]["eigenvalue"][0].as_f64().unwrap());
        }
    }
    (calls, sent, returned, neffs)
}

/// The plan's gate for S7: the P0 workflows in script mode, tool calls and bytes against the per-command baseline, with the same numbers.
#[test]
fn script_mode_is_one_call_and_a_fraction_of_the_bytes_for_the_p0_workflows_with_the_same_numbers() {
    worker_present();
    let repo = std::fs::canonicalize(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")).unwrap();
    let repo = PathBuf::from(repo.to_string_lossy().trim_start_matches(r"\?\"));
    let mut h = start_with_root(&repo);
    h.init("2025-06-18");
    let defs = h.request("tools/list", json!({}))["result"]["tools"].to_string().len();
    println!("tool definitions: {defs} bytes in script mode");
    for (list, script) in [("tools/workflows/w1_soi_width.jsonl", "docs/scripting-samples/w1_soi_width.nys"), ("tools/workflows/w2_soi_wavelength.jsonl", "docs/scripting-samples/w2_soi_wavelength.nys")] {
        let (calls, sent, returned, want) = per_command_baseline(&std::fs::read_to_string(repo.join(list)).unwrap(), &repo);
        let src = std::fs::read_to_string(repo.join(script)).unwrap();
        let req = json!({ "source": src, "allow_destructive": true, "wait_s": 300 });
        let r = h.tool("script_run", req.clone());
        let a = ok(&r);
        assert_eq!(a["status"], "done", "{a}");
        let got: Vec<f64> = a["value"].as_array().unwrap().iter().map(|row| row[1].as_f64().unwrap()).collect();
        assert_eq!(got.len(), want.len(), "{list}");
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).abs() < 1e-9 * w.abs(), "{list}: script {g} against the commands {w}");
        }
        let (s_sent, s_ret) = (req.to_string().len(), r["result"]["content"][0]["text"].as_str().unwrap().len());
        println!("{list}: per command {calls} calls, {sent} bytes sent, {returned} bytes returned | script mode 1 call, {s_sent} bytes sent, {s_ret} bytes returned");
        assert!(s_ret * 5 < returned, "script mode returns {s_ret} bytes against {returned}");
        assert!(calls >= 15, "{list} is a long list of single calls");
    }
}
