//! Interface contract of nyx-mcp (docs/architecture/adr/0015, 0019 and 0057): what an AI host that has registered `nyx` can rely on.
//!
//! The snapshot `docs/architecture/contracts/mcp_tools.v3.json` lists, for every tool, its name, its annotations, its required
//! arguments and every argument with its JSON type (and enum values). Descriptions and defaults are free text and are not pinned.
//! Adding a tool or an optional argument changes the snapshot (additive: bump the minor in `VERSIONS.md`); removing or renaming a tool or an
//! argument, or making an argument required, is breaking (bump the major and tell hosts).
//!
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-mcp --test contracts` rewrites the snapshot.
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

/// Sends the lines and returns one reply per request (a line with an `id`).
fn session(lines: &[Value]) -> Vec<Value> {
    let dir = std::env::temp_dir().join(format!("nyx_mcp_contract_{}", std::process::id()));
    let mut child = Command::new(env!("CARGO_BIN_EXE_nyx-mcp")).arg("--dir").arg(dir).stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::null()).spawn().unwrap();
    let mut tx = child.stdin.take().unwrap();
    let mut rx = BufReader::new(child.stdout.take().unwrap());
    let mut out = vec![];
    for l in lines {
        writeln!(tx, "{l}").unwrap();
        tx.flush().unwrap();
        if l.get("id").is_some() {
            let mut s = String::new();
            rx.read_line(&mut s).unwrap();
            out.push(serde_json::from_str(&s).unwrap_or_else(|e| panic!("bad reply {s:?}: {e}")));
        }
    }
    let _ = child.kill();
    let _ = child.wait();
    out
}

fn init(proto: &str) -> Value {
    json!({ "jsonrpc": "2.0", "id": 1, "method": "initialize", "params": { "protocolVersion": proto, "capabilities": {}, "clientInfo": { "name": "contract", "version": "0" } } })
}

fn tools_snapshot() -> Value {
    let r = session(&[init("2025-06-18"), json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }), json!({ "jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {} })]);
    let mut tools: Vec<Value> = r[1]["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| {
            let schema = &t["inputSchema"];
            let mut args = serde_json::Map::new();
            if let Some(p) = schema["properties"].as_object() {
                for (k, v) in p {
                    let mut a = serde_json::Map::new();
                    a.insert("type".into(), v.get("type").cloned().unwrap_or(Value::Null));
                    if let Some(e) = v.get("enum") {
                        a.insert("enum".into(), e.clone());
                    }
                    args.insert(k.clone(), Value::Object(a));
                }
            }
            let mut required: Vec<Value> = schema.get("required").and_then(|r| r.as_array()).cloned().unwrap_or_default();
            required.sort_by_key(|v| v.as_str().unwrap_or("").to_string());
            json!({ "name": t["name"], "annotations": t.get("annotations").cloned().unwrap_or(Value::Null), "required": required, "arguments": Value::Object(args) })
        })
        .collect();
    tools.sort_by_key(|t| t["name"].as_str().unwrap().to_string());
    json!({ "server": "nyx-mcp", "tools": tools })
}

#[test]
fn mcp_tool_list_matches_the_contract_snapshot() {
    let got = serde_json::to_string_pretty(&tools_snapshot()).unwrap();
    let path = contract_file("mcp_tools.v3.json");
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, format!("{got}\n")).unwrap();
        return;
    }
    let want = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    // compared as values: the order of the keys of an object is not part of the contract (`serde_json` keeps insertion order or sorts,
    // depending on whether any crate of the build enables `preserve_order`)
    let want_value: Value = serde_json::from_str(&want).unwrap();
    assert_eq!(want_value, tools_snapshot(), "the tool list changed: this is an interface change for every registered host (ADR 19)");
}

#[test]
fn mcp_handshake_contract() {
    for proto in ["2024-11-05", "2025-06-18"] {
        let r = session(&[init(proto)]);
        assert_eq!(r[0]["result"]["protocolVersion"], json!(proto), "a supported protocol version is echoed");
        assert_eq!(r[0]["result"]["serverInfo"]["name"], json!("nyx-mcp"));
        assert!(r[0]["result"]["capabilities"]["tools"].is_object());
        assert!(r[0]["result"]["instructions"].as_str().unwrap().contains("script_run"), "hosts are told to read the api first and that the work is scripts");
    }
    // an unknown version falls back to a version the server speaks
    let r = session(&[init("1999-01-01")]);
    assert!(r[0]["result"]["protocolVersion"].as_str().unwrap().starts_with("20"));
    // unknown method is a JSON-RPC error, not a crash
    let r = session(&[init("2025-06-18"), json!({ "jsonrpc": "2.0", "id": 2, "method": "nope" })]);
    assert_eq!(r[1]["error"]["code"], json!(-32601));
}

