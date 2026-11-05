//! A stand-in MCP server for the tests of `nyx_connect::test_server`: it answers `initialize` and `tools/list` over stdio, one JSON line each.
//! `FAKE_MCP=silent` makes it read and never answer (a hung server); `FAKE_MCP=error` answers `tools/list` with a JSON-RPC error;
//! `FAKE_MCP=noise` writes an unrelated notification before each answer. Built by `cargo test`; never shipped.
use serde_json::{json, Value};
use std::io::{BufRead, Write};

fn main() {
    let mode = std::env::var("FAKE_MCP").unwrap_or_default();
    let stdin = std::io::stdin();
    let mut out = std::io::stdout();
    for line in stdin.lock().lines().map_while(Result::ok) {
        let Ok(req) = serde_json::from_str::<Value>(&line) else { continue };
        if mode == "silent" {
            continue;
        }
        if mode == "noise" {
            let _ = writeln!(out, "{}", json!({ "jsonrpc": "2.0", "method": "notifications/message", "params": { "level": "info", "data": "hello" } }));
        }
        let id = req["id"].clone();
        let reply = match req["method"].as_str() {
            Some("initialize") => json!({ "jsonrpc": "2.0", "id": id, "result": { "protocolVersion": "2025-06-18", "capabilities": {}, "serverInfo": { "name": "fake-nyx", "version": "9.9.9" } } }),
            Some("tools/list") if mode == "error" => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": -32603, "message": "boom" } }),
            Some("tools/list") => json!({ "jsonrpc": "2.0", "id": id, "result": { "tools": [{ "name": "geometry_build" }, { "name": "solve" }] } }),
            _ => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": -32601, "message": "unknown method" } }),
        };
        let _ = writeln!(out, "{reply}");
        let _ = out.flush();
    }
}
