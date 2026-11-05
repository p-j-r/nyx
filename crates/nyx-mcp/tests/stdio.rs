//! Drives the real binary over stdio the way an MCP host does.
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

struct Host {
    child: Child,
    tx: ChildStdin,
    rx: BufReader<ChildStdout>,
    next: u64,
}

impl Host {
    fn start() -> Host {
        let dir = std::env::temp_dir().join(format!("nyx_mcp_test_{}", std::process::id()));
        let mut child = Command::new(env!("CARGO_BIN_EXE_nyx-mcp"))
            .arg("--dir")
            .arg(dir)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn nyx-mcp");
        let tx = child.stdin.take().unwrap();
        let rx = BufReader::new(child.stdout.take().unwrap());
        Host { child, tx, rx, next: 1 }
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
    fn init(&mut self, proto: &str) -> Value {
        let r = self.request("initialize", json!({ "protocolVersion": proto, "capabilities": {}, "clientInfo": { "name": "test", "version": "0" } }));
        self.send(json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }));
        r
    }
}

impl Drop for Host {
    fn drop(&mut self) {
        let _ = self.child.kill();
    }
}

fn payload(r: &Value) -> Value {
    serde_json::from_str(r["result"]["content"][0]["text"].as_str().unwrap()).unwrap()
}

#[test]
fn handshake_negotiates_the_protocol_and_lists_tools() {
    let mut h = Host::start();
    let r = h.init("2024-11-05");
    assert_eq!(r["result"]["protocolVersion"], "2024-11-05");
    assert_eq!(r["result"]["serverInfo"]["name"], "nyx-mcp");
    assert!(r["result"]["instructions"].as_str().unwrap().contains("get_guide"));
    assert_eq!(h.request("ping", json!({}))["result"], json!({}));
    let l = h.request("tools/list", json!({}));
    let names: Vec<&str> = l["result"]["tools"].as_array().unwrap().iter().map(|t| t["name"].as_str().unwrap()).collect();
    assert!(names.contains(&"get_guide") && names.contains(&"add_rect") && names.contains(&"solve_run") && names.contains(&"analyse"));
    for geometry in ["add_circle", "add_polygon", "add_overlay", "add_pml", "add_difference"] {
        assert!(names.contains(&geometry), "the geometry tools are generated from the descriptors: {geometry} is missing");
    }
    for t in l["result"]["tools"].as_array().unwrap() {
        assert_eq!(t["inputSchema"]["type"], "object", "{}", t["name"]);
    }
}

#[test]
fn unknown_version_falls_back_and_unknown_method_is_an_error() {
    let mut h = Host::start();
    assert_eq!(h.init("1999-01-01")["result"]["protocolVersion"], "2025-06-18");
    assert_eq!(h.request("nope/nothing", json!({}))["error"]["code"], -32601);
    assert_eq!(h.tool("no_such_tool", json!({}))["error"]["code"], -32602);
}

#[test]
fn guide_is_plain_text() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("get_guide", json!({ "topic": "nyx-instructions" }));
    assert_eq!(r["result"]["isError"], false);
    assert!(r["result"]["content"][0]["text"].as_str().unwrap().contains("add_overlay"));
}

#[test]
fn bad_arguments_are_a_tool_error_that_names_the_valid_values() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("get_guide", json!({ "topic": "bogus" }));
    assert_eq!(r["result"]["isError"], true);
    let p = payload(&r);
    assert_eq!(p["error"], "invalid");
    assert!(p["message"].as_str().unwrap().contains("optical-modes"), "the error names the valid topics: {p}");
}


fn ok(r: &Value) -> Value {
    assert_eq!(r["result"]["isError"], false, "{r}");
    payload(r)
}

#[test]
fn a_host_builds_a_wire_solves_it_and_the_work_is_a_stored_workflow() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let steps: Vec<(&str, Value)> = vec![
        ("project_new", json!({ "title": "wire from mcp", "h": 1.2e-7 })),
        ("add_rect", json!({ "name": "box", "min": [-2.5e-6, -2e-6], "max": [2.5e-6, 0.0] })),
        ("add_rect", json!({ "name": "clad", "min": [-2.5e-6, 0.0], "max": [2.5e-6, 1.5e-6] })),
        ("add_rect", json!({ "name": "core", "min": [-2.5e-7, 0.0], "max": [2.5e-7, 3e-7] })),
        ("add_overlay", json!({ "name": "regions", "regions": [
            { "name": "BOX", "input": "box" }, { "name": "clad", "input": "clad" }, { "name": "Si", "input": "core" } ] })),
        ("add_physics", json!({ "physics": "ems", "study": "modes", "order": 2,
            "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.9, "n_modes": 2 } })),
        ("set_material", json!({ "physics": "ems", "region": "BOX", "material": "silica" })),
        ("set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" })),
        ("set_material", json!({ "physics": "ems", "region": "Si", "material": "silicon" })),
        ("set_mesh", json!({ "adapt_scale": 0.4 })),
    ];
    let mut last = Value::Null;
    for (tool, args) in steps {
        last = ok(&h.tool(tool, args));
    }
    assert!(last["nyx"]["revision"].as_u64().unwrap() >= 9, "every edit moves the revision: {last}");
    let hash_before = last["nyx"]["project_hash"].clone();

    let solved = ok(&h.tool("solve_run", json!({})));
    assert_eq!(solved["results"][0]["physics"], "ems", "{solved}");
    assert!(solved["run"]["id"].is_string(), "the run is recorded in the workflow's run table: {solved}");
    let info = ok(&h.tool("result_info", json!({ "physics": "ems" })));
    let neff = info["result"][0]["solutions"][0]["eigenvalue"][0].as_f64().unwrap();
    assert!((neff - 2.67).abs() < 0.03, "TE0 effective index {neff}");
    assert_eq!(info["nyx"]["project_hash"], hash_before, "reading results does not change the project");

    let a = ok(&h.tool("analyse", json!({ "analysis": "polarization", "physics": "ems" })));
    assert!(a["rows"][0]["report"]["class"].as_str().unwrap().contains("TE"), "{a}");

    // the host did not name a folder: the work sits in the server's workspace as a workflow with a run
    let runs = ok(&h.tool("store_runs", json!({})));
    assert_eq!(runs["result"].as_array().map(|r| r.len()), Some(1), "{runs}");
    let status = ok(&h.tool("project_info", json!({})));
    let folder = status["store"]["folder"].as_str().expect("bound to a workflow folder");
    assert!(std::path::Path::new(folder).join("workflow.json").exists(), "{folder}");
}
