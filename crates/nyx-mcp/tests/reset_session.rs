// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! MCP `reset_session` and UTF-8 BOM handling tests over stdio (FR-7): verifies that `script_run` with
//! `reset_session: true` clears REPL globals (standalone or before running a script) while preserving the
//! `Engine` project, and that leading UTF-8 BOMs are stripped in `script_run` and `script_check`.

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
            "nyx_mcp_reset_test_{}_{}",
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
fn reset_session_clears_repl_globals_without_touching_engine_and_handles_bom() {
    worker_present();
    let mut h = Host::start();

    // Schema lists reset_session on script_run
    let tools = h.request("tools/list", json!({}));
    let run_tool = tools["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .find(|t| t["name"] == "script_run")
        .unwrap();
    assert_eq!(run_tool["inputSchema"]["properties"]["reset_session"]["type"], "boolean");

    // Populate REPL globals and engine project
    let a = ok(&h.tool(
        "script_run",
        json!({
            "allow_destructive": true,
            "source": "nyx.project_new(dim=2, title='persist_project')\nnyx.set_param(name='w', value=5e-7)\nrepl_var = 123\nrepl_var"
        }),
    ));
    assert_eq!(a["value"], json!(123));

    // Standalone reset_session without source/path
    let reset_ans = ok(&h.tool("script_run", json!({ "reset_session": true })));
    assert_eq!(reset_ans["ok"], json!(true));
    assert_eq!(reset_ans["reset_session"], json!(true));

    // repl_var is gone, but engine project is intact
    let err = failure(&h.tool("script_run", json!({ "source": "repl_var" })));
    assert!(err.contains("repl_var"), "{err}");
    let title = ok(&h.tool("script_run", json!({ "source": "nyx.project_info()['title']" })));
    assert_eq!(title["value"], json!("persist_project"));

    // Define another variable, then reset_session + run new source in one call
    ok(&h.tool("script_run", json!({ "source": "temp_var = 77" })));
    let combined = ok(&h.tool(
        "script_run",
        json!({
            "reset_session": true,
            "source": "fresh_var = 99\nfresh_var"
        }),
    ));
    assert_eq!(combined["value"], json!(99));
    let err = failure(&h.tool("script_run", json!({ "source": "temp_var" })));
    assert!(err.contains("temp_var"), "{err}");

    // UTF-8 BOM in inline source and saved script
    let bom_chk = ok(&h.tool("script_check", json!({ "source": "\u{feff}fresh_var + 1" })));
    assert_eq!(bom_chk["status"], json!("ok"));
    let bom_run = ok(&h.tool(
        "script_run",
        json!({
            "source": "\u{feff}fresh_var + 1",
            "save_as": "bom_saved"
        }),
    ));
    assert_eq!(bom_run["value"], json!(100));
    let bom_path = ok(&h.tool("script_run", json!({ "path": "bom_saved" })));
    assert_eq!(bom_path["value"], json!(100));
}
