//! nyx-mcp: a Model Context Protocol server over stdio for nyx. Newline-delimited JSON-RPC 2.0 on stdin/stdout, logs on stderr.
//! No async runtime and no SDK: the protocol surface nyx needs is five methods. The server holds one `nyx_engine::Engine` and nothing else:
//! the tool list is the engine's command catalog (the geometry tools among them are generated from the descriptors of the geom
//! operations) and a tool call is `Engine::call`. What a host does is stored as a workflow of the workspace folder (origin `ai`), so the
//! person sees it in the window's dashboard.
use std::io::{BufRead, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;

use nyx_engine::guide::SERVER_INSTRUCTIONS;
use nyx_engine::{Access, Engine};
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

/// Commands that discard the project or a part of it: a host may ask the person first.
fn destructive(name: &str) -> bool {
    matches!(name, "project_new" | "project_load" | "project_open" | "store_open") || name.starts_with("remove_")
}

fn tool_list(engine: &Engine) -> Value {
    let tools: Vec<Value> = engine
        .catalog()
        .iter()
        .map(|c| {
            json!({
                "name": c.name, "title": c.title, "description": c.description, "inputSchema": c.schema,
                "annotations": {
                    "title": c.title,
                    "readOnlyHint": c.access == Access::Read,
                    "destructiveHint": destructive(&c.name),
                    "idempotentHint": c.access == Access::Read,
                    "openWorldHint": false
                }
            })
        })
        .collect();
    json!(tools)
}

/// Protocol revisions compare as dates, so string order is chronological.
fn at_least(proto: &str, min: &str) -> bool {
    proto >= min
}

/// A short id for this process: the session a result belongs to.
fn session_id() -> String {
    let t = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    format!("s_{:08x}", (t ^ ((std::process::id() as u128) << 20)) as u32)
}

struct Server {
    engine: Engine,
    proto: &'static str,
    session: String,
    /// counts the calls that changed the project, so a host can tell that something moved
    revision: u64,
}

impl Server {
    fn new(dir: PathBuf) -> Server {
        let mut engine = Engine::new();
        engine.set_workspace(dir, true);
        engine.set_owner("ai"); // its requests queue behind a person's at the same priority (ADR 87)
        Server { engine, proto: PROTOCOLS[0], session: session_id(), revision: 0 }
    }

    fn header(&self) -> Value {
        let hash = self.engine.stamp();
        json!({ "session": self.session, "revision": self.revision, "project_hash": hash.get(..16).unwrap_or(&hash) })
    }

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

    /// Runs one command: its JSON with the session header, or the error as JSON. `None` when there is no such tool.
    fn run(&mut self, name: &str, args: &Value) -> Option<(Value, bool)> {
        let access = self.engine.catalog().into_iter().find(|c| c.name == name)?.access;
        let out = catch_unwind(AssertUnwindSafe(|| self.engine.call(name, args)));
        Some(match out {
            Ok(Ok(v)) => {
                if access != Access::Read {
                    self.revision += 1;
                }
                let mut v = match v {
                    Value::Object(_) => v,
                    other => json!({ "result": other }),
                };
                v["nyx"] = self.header();
                (v, false)
            }
            Ok(Err(e)) => (json!({ "error": e.code(), "message": e.to_string(), "nyx": self.header() }), true),
            Err(_) => (
                json!({ "error": "internal", "message": "nyx hit an internal fault on this call; the session is still usable", "nyx": self.header() }),
                true,
            ),
        })
    }

    fn tools_call(&mut self, params: &Value) -> Result<Value, (i64, String)> {
        let name = params.get("name").and_then(Value::as_str).ok_or((-32602, "tools/call needs a string 'name'".to_string()))?.to_string();
        let args = match params.get("arguments") {
            None | Some(Value::Null) => json!({}),
            Some(a) => a.clone(),
        };
        let Some((value, is_error)) = self.run(&name, &args) else {
            return Err((-32602, format!("unknown tool '{name}'")));
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
            "tools/list" => Ok(json!({ "tools": tool_list(&self.engine) })),
            "tools/call" => self.tools_call(&params),
            m => Err((-32601, format!("method not found: {m}"))),
        };
        Some(match r {
            Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }),
            Err((code, message)) => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } }),
        })
    }
}

/// Text for the host: a guide is plain text, everything else compact JSON.
fn render(tool: &str, v: &Value) -> String {
    if tool == "get_guide" {
        if let Some(t) = v.get("text").and_then(Value::as_str) {
            return t.to_string();
        }
    }
    v.to_string()
}

fn serve(dir: PathBuf) {
    let mut s = Server::new(dir);
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
    println!("nyx-mcp {VERSION}\n\nUSAGE\n  nyx-mcp [--dir PATH]                 serve MCP on stdio\n  nyx-mcp --list                       print the tool names\n  nyx-mcp [--dir PATH] --call TOOL JSON  run one tool and print its result\n\nWorkflows go to PATH (default {}).", default_dir().display());
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
            for c in Engine::new().catalog() {
                println!("{:22} {}", c.name, c.title);
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
            let mut s = Server::new(dir);
            match s.run(tool, &a) {
                Some((v, false)) => println!("{}", render(tool, &v)),
                Some((v, true)) => {
                    println!("{v}");
                    std::process::exit(1);
                }
                None => {
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
