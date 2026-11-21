// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-session MCP tests over stdio (QFR-7 part 2): verifies isolated per-session `Engine` and `Kernel`
//! state (`"s1"`, `"s2"`, `"default"`), session-scoped `reset_session: true`, invalid session name rejection,
//! and `"cqed"` topic exposure in `api`.

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

fn worker_present() {
    let w = PathBuf::from(env!("CARGO_BIN_EXE_nyx-mcp"))
        .with_file_name(format!("nyx-script-worker{}", std::env::consts::EXE_SUFFIX));
    assert!(
        w.is_file(),
        "the script worker is not at {}: build it with `cargo build --release -p nyx-script`",
        w.display()
    );
}

impl Host {
    fn start() -> Host {
        static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "nyx_mcp_multisess_test_{}_{}",
            std::process::id(),
            N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
        ));
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
        let mut h = Host { child, tx, rx, next: 1, dir };
        let _ = h.request(
            "initialize",
            json!({
                "protocolVersion": "2025-06-18",
                "capabilities": {},
                "clientInfo": { "name": "test", "version": "0" }
            }),
        );
        writeln!(h.tx, "{}", json!({ "jsonrpc": "2.0", "method": "notifications/initialized" })).unwrap();
        h.tx.flush().unwrap();
        h
    }

    fn request(&mut self, method: &str, params: Value) -> Value {
        let id = self.next;
        self.next += 1;
        writeln!(self.tx, "{}", json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params })).unwrap();
        self.tx.flush().unwrap();
        let mut line = String::new();
        self.rx.read_line(&mut line).unwrap();
        let v: Value = serde_json::from_str(&line).unwrap();
        assert_eq!(v["id"], json!(id));
        v
    }

    fn tool(&mut self, name: &str, args: Value) -> Value {
        self.request("tools/call", json!({ "name": name, "arguments": args }))
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
fn tools_schema_exposes_session_and_cqed_topic() {
    let mut h = Host::start();
    let tools = h.request("tools/list", json!({}))["result"]["tools"].as_array().unwrap().clone();
    for tool_name in ["script_run", "script_wait", "script_check"] {
        let t = tools.iter().find(|t| t["name"] == tool_name).unwrap();
        assert_eq!(t["inputSchema"]["properties"]["session"]["type"], "string", "{tool_name}");
        assert_eq!(t["inputSchema"]["properties"]["session"]["default"], "default", "{tool_name}");
    }
    let api_tool = tools.iter().find(|t| t["name"] == "api").unwrap();
    let topics: Vec<&str> = api_tool["inputSchema"]["properties"]["topic"]["enum"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert!(topics.contains(&"cqed"), "api topic enum should include cqed: {topics:?}");
}

#[test]
fn interleaved_sessions_maintain_independent_engines_repl_globals_and_scoped_reset() {
    worker_present();
    let mut h = Host::start();

    // Session s1: create project 'proj_s1', param w = 4e-7, REPL global x = 10
    let r1 = ok(&h.tool(
        "script_run",
        json!({
            "session": "s1",
            "allow_destructive": true,
            "source": "nyx.project_new(dim=2, title='proj_s1')\nnyx.set_param(name='w', value=4e-7)\nx = 10\ndef f(v):\n    return v + x\nf(5)"
        }),
    ));
    assert_eq!(r1["value"], json!(15));
    assert_eq!(r1["nyx"]["session_name"], json!("s1"));

    // Session s2: create project 'proj_s2', param w = 8e-7, REPL global x = 200
    let r2 = ok(&h.tool(
        "script_run",
        json!({
            "session": "s2",
            "allow_destructive": true,
            "source": "nyx.project_new(dim=3, title='proj_s2')\nnyx.set_param(name='w', value=8e-7)\nx = 200\ndef f(v):\n    return v * x\nf(3)"
        }),
    ));
    assert_eq!(r2["value"], json!(600));
    assert_eq!(r2["nyx"]["session_name"], json!("s2"));
    assert_ne!(r1["nyx"]["session"], r2["nyx"]["session"]);
    assert_ne!(r1["nyx"]["project_hash"], r2["nyx"]["project_hash"]);

    // Interleaved reads confirm s1 and s2 do not leak projects, parameters, or REPL globals
    let check_s1 = ok(&h.tool(
        "script_run",
        json!({
            "session": "s1",
            "source": "info = nyx.project_info()\n{'title': info['title'], 'dim': info['dim'], 'w': info['params'][0]['value'], 'x': x, 'f': f(1)}"
        }),
    ));
    assert_eq!(
        check_s1["value"],
        json!({ "title": "proj_s1", "dim": 2, "w": 4e-7, "x": 10, "f": 11 })
    );

    let check_s2 = ok(&h.tool(
        "script_run",
        json!({
            "session": "s2",
            "source": "info = nyx.project_info()\n{'title': info['title'], 'dim': info['dim'], 'w': info['params'][0]['value'], 'x': x, 'f': f(2)}"
        }),
    ));
    assert_eq!(
        check_s2["value"],
        json!({ "title": "proj_s2", "dim": 3, "w": 8e-7, "x": 200, "f": 400 })
    );

    // reset_session on s1 clears s1's REPL globals while preserving s1's engine project AND leaving s2 untouched
    let rst1 = ok(&h.tool("script_run", json!({ "session": "s1", "reset_session": true })));
    assert_eq!(rst1["ok"], json!(true));
    assert_eq!(rst1["reset_session"], json!(true));

    // In s1, x is now undefined, but proj_s1 is still in s1's engine
    let err_s1 = failure(&h.tool("script_check", json!({ "session": "s1", "source": "x + 1" })));
    assert!(err_s1.contains("x"), "{err_s1}");
    let s1_proj = ok(&h.tool("script_run", json!({ "session": "s1", "source": "nyx.project_info()['title']" })));
    assert_eq!(s1_proj["value"], json!("proj_s1"));

    // In s2, x and f(v) and proj_s2 are all still intact!
    assert_eq!(ok(&h.tool("script_check", json!({ "session": "s2", "source": "f(x)" })))["status"], json!("ok"));
    let s2_still = ok(&h.tool(
        "script_run",
        json!({
            "session": "s2",
            "source": "[nyx.project_info()['title'], x, f(4)]"
        }),
    ));
    assert_eq!(s2_still["value"], json!(["proj_s2", 200, 800]));
}

#[test]
fn invalid_session_names_are_rejected() {
    let mut h = Host::start();
    for bad in ["", ".", "..", "../escape", "a/b", "a\\b", "bad space"] {
        for tool in ["script_run", "script_wait", "script_check"] {
            let err = failure(&h.tool(tool, json!({ "session": bad, "source": "1" })));
            assert!(err.contains("invalid session name"), "{tool} with {bad:?}: {err}");
        }
    }
}
