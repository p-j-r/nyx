//! nyx-mcp: a Model Context Protocol server over stdio for nyx. Newline-delimited JSON-RPC 2.0 on stdin/stdout, logs on stderr.
//! No async runtime and no SDK: the protocol surface nyx needs is five methods. All behaviour lives in `nyx-session`.
use std::io::{BufRead, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;

use nyx_session::guide::SERVER_INSTRUCTIONS;
use nyx_session::{call, render, Outcome, Workspace, TOOLS};
use serde_json::{json, Value};

const NAME: &str = "nyx-mcp";
const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Protocol revisions this server speaks, newest first. It answers `initialize` with the client's if listed, else the newest.
const PROTOCOLS: &[&str] = &["2025-06-18", "2025-03-26", "2024-11-05"];

fn default_dir() -> PathBuf {
    let home = std::env::var_os("USERPROFILE").or_else(|| std::env::var_os("HOME")).map(PathBuf::from);
    match home {
        Some(h) => h.join("Documents").join("nyx_runs").join("mcp"),
        None => std::env::temp_dir().join("nyx_runs").join("mcp"),
    }
}

fn tool_list() -> Value {
    let tools: Vec<Value> = TOOLS
        .iter()
        .map(|t| json!({ "name": t.name, "title": t.title, "description": t.description, "inputSchema": t.input_schema(), "annotations": t.annotations() }))
        .collect();
    json!(tools)
}

/// Protocol revisions compare as dates, so string order is chronological.
fn at_least(proto: &str, min: &str) -> bool {
    proto >= min
}

struct Server {
    ws: Workspace,
    proto: &'static str,
}

impl Server {
    fn initialize(&mut self, params: &Value) -> Value {
        let asked = params.get("protocolVersion").and_then(Value::as_str).unwrap_or("");
        self.proto = PROTOCOLS.iter().copied().find(|p| *p == asked).unwrap_or(PROTOCOLS[0]);
        json!({
            "protocolVersion": self.proto,
            "capabilities": { "tools": { "listChanged": false } },
            "serverInfo": { "name": NAME, "title": "nyx photonics", "version": VERSION },
            "instructions": SERVER_INSTRUCTIONS,
        })
    }

    fn tools_call(&mut self, params: &Value) -> Result<Value, (i64, String)> {
        let name = params.get("name").and_then(Value::as_str).ok_or((-32602, "tools/call needs a string 'name'".to_string()))?.to_string();
        let args = match params.get("arguments") {
            None | Some(Value::Null) => json!({}),
            Some(a) => a.clone(),
        };
        let out = catch_unwind(AssertUnwindSafe(|| call(&mut self.ws, &name, &args)));
        let (value, is_error) = match out {
            Ok(Outcome::Ok(v)) | Ok(Outcome::NeedsInput(v)) => (v, false),
            Ok(Outcome::Err(e)) => (e.to_json(), true),
            Ok(Outcome::UnknownTool) => return Err((-32602, format!("unknown tool '{name}'"))),
            Err(_) => (json!({ "error": "internal", "message": "nyx hit an internal fault on this call; the session is still usable" }), true),
        };
        let mut result = json!({ "content": [{ "type": "text", "text": render(&name, &value) }], "isError": is_error });
        if at_least(self.proto, "2025-06-18") && !is_error && name != "get_guide" {
            result["structuredContent"] = value;
        }
        Ok(result)
    }

    /// One incoming message; `None` means no reply (notifications).
    fn handle(&mut self, msg: &Value) -> Option<Value> {
        let id = msg.get("id").cloned();
        let method = msg.get("method").and_then(Value::as_str).unwrap_or("");
        let params = msg.get("params").cloned().unwrap_or(Value::Null);
        // notifications carry no id and get no reply; so do stray responses
        let id = match id {
            Some(i) if !method.is_empty() => i,
            _ => return None,
        };
        let r = match method {
            "initialize" => Ok(self.initialize(&params)),
            "ping" => Ok(json!({})),
            "tools/list" => Ok(json!({ "tools": tool_list() })),
            "tools/call" => self.tools_call(&params),
            m => Err((-32601, format!("method not found: {m}"))),
        };
        Some(match r {
            Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }),
            Err((code, message)) => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } }),
        })
    }
}

fn serve(dir: PathBuf) {
    let mut s = Server { ws: Workspace::new(dir), proto: PROTOCOLS[0] };
    eprintln!("{NAME} {VERSION} ready");
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    for line in stdin.lock().lines() {
        let Ok(line) = line else { break };
        if line.trim().is_empty() {
            continue;
        }
        let reply = match serde_json::from_str::<Value>(&line) {
            Ok(Value::Array(batch)) => {
                // 2025-06-18 dropped batching; answering one is harmless
                let out: Vec<Value> = batch.iter().filter_map(|m| s.handle(m)).collect();
                if out.is_empty() { None } else { Some(Value::Array(out)) }
            }
            Ok(m) => s.handle(&m),
            Err(e) => Some(json!({ "jsonrpc": "2.0", "id": null, "error": { "code": -32700, "message": format!("parse error: {e}") } })),
        };
        if let Some(r) = reply {
            let mut out = stdout.lock();
            if writeln!(out, "{r}").and_then(|_| out.flush()).is_err() {
                break; // the host went away
            }
        }
    }
    eprintln!("{NAME} stdin closed, exiting");
}

fn usage() {
    println!("nyx-mcp {VERSION}\n\nUSAGE\n  nyx-mcp [--dir PATH]                 serve MCP on stdio\n  nyx-mcp --list                       print the tool names\n  nyx-mcp [--dir PATH] --call TOOL JSON  run one tool and print its result\n\nSession files go to PATH (default {}).", default_dir().display());
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut dir = default_dir();
    if let Some(i) = args.iter().position(|a| a == "--dir") {
        if i + 1 >= args.len() {
            eprintln!("--dir needs a path");
            std::process::exit(2);
        }
        dir = PathBuf::from(args.remove(i + 1));
        args.remove(i);
    }
    match args.first().map(String::as_str) {
        None => serve(dir),
        Some("--version") => println!("{NAME} {VERSION}"),
        Some("--help") | Some("-h") => usage(),
        Some("--list") => {
            for t in TOOLS {
                println!("{:16} {}", t.name, t.title);
            }
        }
        Some("--call") => {
            let (Some(tool), json_arg) = (args.get(1), args.get(2)) else {
                eprintln!("--call needs a tool name and optional JSON arguments");
                std::process::exit(2);
            };
            let a: Value = match json_arg.map(|s| serde_json::from_str(s)).unwrap_or(Ok(json!({}))) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("arguments are not JSON: {e}");
                    std::process::exit(2);
                }
            };
            let mut ws = Workspace::new(dir);
            match call(&mut ws, tool, &a) {
                Outcome::Ok(v) | Outcome::NeedsInput(v) => println!("{}", render(tool, &v)),
                Outcome::Err(e) => {
                    println!("{}", e.to_json());
                    std::process::exit(1);
                }
                Outcome::UnknownTool => {
                    eprintln!("unknown tool '{tool}' (see --list)");
                    std::process::exit(2);
                }
            }
        }
        Some(other) => {
            eprintln!("unknown argument '{other}'");
            usage();
            std::process::exit(2);
        }
    }
}
