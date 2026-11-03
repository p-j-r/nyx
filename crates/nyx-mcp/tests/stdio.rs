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
    assert!(names.contains(&"recipe_run") && names.contains(&"get_guide"));
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
    assert!(r["result"]["content"][0]["text"].as_str().unwrap().contains("recipe_run"));
}

#[test]
fn bad_arguments_are_a_tool_error_that_names_the_valid_values() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("get_guide", json!({ "topic": "bogus" }));
    assert_eq!(r["result"]["isError"], true);
    let p = payload(&r);
    assert_eq!(p["error"], "bad_value");
    assert!(p["valid"].as_array().unwrap().iter().any(|v| v == "optical-modes"));
}

#[test]
fn recipe_asks_for_the_width_then_answers() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("recipe_run", json!({ "recipe": "wg_modes_2d", "thickness": 3e-7 }));
    assert_eq!(r["result"]["isError"], false);
    let p = payload(&r);
    assert_eq!(p["status"], "needs_input", "{p}");
    assert_eq!(p["slot"], "width");

    let r = h.tool("recipe_run", json!({ "recipe": "wg_modes_2d", "thickness": 3e-7, "width": 5e-7, "cladding": "oxide", "box": 2e-6 }));
    let p = payload(&r);
    assert_eq!(p["status"], "done", "{p}");
    let neff = p["result"]["neff"].as_f64().expect("result.neff");
    assert!((neff - 2.6718).abs() < 0.01, "TE0 neff {neff} vs documented 2.6718");
    assert!(p["checks"].as_array().unwrap().iter().all(|c| c["pass"] == json!(true)), "{p}");
    assert!(r["result"]["structuredContent"].is_object());
}

#[test]
fn a_wide_tm_wire_is_reported_unverified_not_done() {
    let mut h = Host::start();
    h.init("2025-06-18");
    let r = h.tool("recipe_run", json!({ "recipe": "wg_modes_2d", "thickness": 3e-7, "width": 4e-6, "cladding": "air", "pol": "TM" }));
    let p = payload(&r);
    assert_eq!(p["status"], "unverified", "{p}");
    assert!(p["answer"].as_str().unwrap().starts_with("UNVERIFIED"));
}
