//! nyx-mcp: a Model Context Protocol server over stdio for nyx. Newline-delimited JSON-RPC 2.0 on stdin/stdout, logs on stderr.
//! No async runtime and no SDK: the protocol surface nyx needs is five methods. Since S7 (ADR 94) the server offers four tools, not one per engine
//! command: `script_run`, `script_wait`, `script_check` and `api`. The server holds one `nyx_engine::Engine` and one script kernel
//! (`nyx_script::kernel`, created at the first script, so `initialize` and `tools/list` start no process): the agent writes Python that calls the
//! engine commands and only what the script prints or returns comes back. What a host does is stored as a workflow of the workspace folder
//! (origin `ai`), so the person sees it in the window's dashboard, with the script that made each run.
use std::io::{BufRead, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::time::Duration;

use nyx_engine::guide::{self, SERVER_INSTRUCTIONS};
use nyx_engine::Engine;
use nyx_script::kernel::{Kernel, Step};
use nyx_script::{stubs, Limits, Options};
use serde_json::{json, Value};

const NAME: &str = "nyx-mcp";
const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Protocol revisions this server speaks, newest first. It answers `initialize` with the client's if listed, else the newest.
const PROTOCOLS: &[&str] = &["2025-06-18", "2025-03-26", "2024-11-05"];
/// How long `script_run` and `script_wait` give a script before answering `running`, and the most they may be asked to.
const DEFAULT_WAIT_S: f64 = 25.0;
const MAX_WAIT_S: f64 = 900.0;
/// A script's value is cut past this many bytes of JSON: the point of a script is that fields stay in the engine.
const MAX_VALUE_BYTES: usize = 48 << 10;

fn home() -> Option<PathBuf> {
    std::env::var_os("USERPROFILE").or_else(|| std::env::var_os("HOME")).map(PathBuf::from)
}

fn default_dir() -> PathBuf {
    nyx_script::mcp_dir()
}

/// The folder the file arguments of a script's commands must lie in (the engine's path guard, ADR 93): `--root`, else the workspace when
/// `--dir` was given, else the person's Documents folder (their projects are there, and so is the default workspace).
fn default_root(dir_given: bool, dir: &Path) -> PathBuf {
    if dir_given {
        return dir.to_path_buf();
    }
    home().map(|h| h.join("Documents")).unwrap_or_else(|| dir.to_path_buf())
}

fn tools() -> Value {
    let wait = json!({ "type": "number", "description": format!("seconds to let the script run before answering `running` (default {DEFAULT_WAIT_S}, at most {MAX_WAIT_S})") });
    json!([
        {
            "name": "script_run", "title": "Run a script",
            "description": "Runs a Python script that calls the nyx engine commands as nyx.<command>(name=value, ...) and answers {status, value, printed, calls, nyx}. `value` is the script's last expression, `printed` what it printed: print the numbers you need, fields stay in the engine. Variables and functions stay from one call to the next. The script is checked before it runs. status `running` means a solve is still going: call script_wait. Read api first.",
            "inputSchema": { "type": "object", "properties": {
                "source": { "type": "string", "description": "the script text" },
                "path": { "type": "string", "description": "the name of a script saved in the workspace instead of `source` (letters, digits, - and _)" },
                "params": { "type": "object", "description": "values for the script, which reads them as the dict `params`" },
                "wait_s": wait,
                "allow_destructive": { "type": "boolean", "description": "allow project_new, project_load, project_open, store_open and remove_* (they discard the project or a part of it); ask the person first" },
                "save_as": { "type": "string", "description": "also save `source` in the workspace under this name, to run again with `path`" }
            } },
            "annotations": { "title": "Run a script", "readOnlyHint": false, "destructiveHint": false, "idempotentHint": false, "openWorldHint": false }
        },
        {
            "name": "script_wait", "title": "Continue or cancel a waiting script",
            "description": "Continues the script that script_run left `running` (a solve or a convergence check was still going) and answers like script_run. With cancel true it abandons the script instead: its job is cancelled and the variables of the scripts that finished stay.",
            "inputSchema": { "type": "object", "properties": { "wait_s": wait, "cancel": { "type": "boolean", "description": "abandon the waiting script" } } },
            "annotations": { "title": "Continue or cancel a waiting script", "readOnlyHint": false, "destructiveHint": false, "idempotentHint": false, "openWorldHint": false }
        },
        {
            "name": "script_check", "title": "Check a script",
            "description": "Parses and type-checks a script against the nyx commands and the variables of the session and runs nothing: a misspelled command or argument, a missing argument, a syntax error, a branch on a value of a traced function. Answers {status: ok} or the error with its line.",
            "inputSchema": { "type": "object", "properties": { "source": { "type": "string", "description": "the script text" }, "path": { "type": "string", "description": "the name of a saved script instead of `source`" } } },
            "annotations": { "title": "Check a script", "readOnlyHint": true, "destructiveHint": false, "idempotentHint": true, "openWorldHint": false }
        },
        {
            "name": "api", "title": "The nyx API",
            "description": "Read this first. Without arguments: how to work with nyx. `topic` is a guide (nyx-instructions, geometry, optical-modes, troubleshooting) or `commands` (every command with its title). `command` is one command with its Python signature, its documentation and its arguments.",
            "inputSchema": { "type": "object", "properties": {
                "topic": { "type": "string", "enum": ["nyx-instructions", "geometry", "optical-modes", "troubleshooting", "commands"], "description": "which guide; omit for nyx-instructions" },
                "command": { "type": "string", "description": "a command name, e.g. add_rect" }
            } },
            "annotations": { "title": "The nyx API", "readOnlyHint": true, "destructiveHint": false, "idempotentHint": true, "openWorldHint": false }
        }
    ])
}

/// A short id for this process: the session a result belongs to.
fn session_id() -> String {
    let t = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    format!("s_{:08x}", (t ^ ((std::process::id() as u128) << 20)) as u32)
}

/// A saved script's file stem: letters, digits, `-` and `_`, with `.nys` allowed and dropped; nothing that names a folder.
fn script_stem(name: &str) -> Result<String, String> {
    let stem = name.strip_suffix(".nys").unwrap_or(name);
    if stem.is_empty() || stem.len() > 80 || !stem.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_') {
        return Err(format!("`{name}` is not a script name: letters, digits, - and _ (the file is scripts/<name>.nys of the workspace)"));
    }
    Ok(stem.to_string())
}

/// A value too long for a model's context is replaced by its start and its size.
fn cap_value(v: Value) -> Value {
    let s = v.to_string();
    if s.len() <= MAX_VALUE_BYTES {
        return v;
    }
    let mut end = 2000;
    while !s.is_char_boundary(end) {
        end -= 1;
    }
    json!({ "truncated": true, "bytes": s.len(), "head": &s[..end], "hint": "print the numbers you need from the script instead of returning a table" })
}

struct Server {
    engine: Engine,
    kernel: Option<Kernel>,
    dir: PathBuf,
    root: PathBuf,
    proto: &'static str,
    session: String,
    /// counts the scripts that ran a command, so a host can tell that something moved
    revision: u64,
}

impl Server {
    fn new(dir: PathBuf, root: PathBuf) -> Server {
        let mut engine = Engine::new();
        engine.set_workspace(dir.clone(), true);
        engine.set_owner("ai"); // its requests queue behind a person's at the same priority (ADR 87)
        engine.set_self_as_worker(true); // solves run in a child process: running out of memory fails the job, not the server (ADR 86)
        engine.set_history_file(&dir.join("solve_history.jsonl")); // what solves cost, for memory estimates
        Server { engine, kernel: None, dir, root, proto: PROTOCOLS[0], session: session_id(), revision: 0 }
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

    fn kernel(&mut self) -> Result<&mut Kernel, String> {
        if self.kernel.is_none() {
            // the path guard wants its folder to exist; the workspace is made here the first time a script needs it
            for d in [&self.dir, &self.root] {
                std::fs::create_dir_all(d).map_err(|e| format!("cannot create {}: {e}", d.display()))?;
            }
            // what a script prints and returns is what costs a model's context: a smaller cap than a person's `nyx run`
            let limits = Limits { output_bytes: 64 << 10, ..Limits::default() };
            self.kernel = Some(Kernel::new(&self.engine, limits, Options::new(self.root.clone()))?);
        }
        Ok(self.kernel.as_mut().expect("the kernel was just made"))
    }

    fn scripts_dir(&self) -> PathBuf {
        self.dir.join("scripts")
    }

    /// The script a call names: its text and the name its errors and its run records carry.
    fn script_of(&self, args: &Value) -> Result<(String, String), String> {
        match (args.get("source").and_then(Value::as_str), args.get("path").and_then(Value::as_str)) {
            (Some(_), Some(_)) => Err("give `source` or `path`, not both".to_string()),
            (Some(src), None) => {
                let name = match args.get("save_as").and_then(Value::as_str) {
                    Some(n) => format!("{}.nys", script_stem(n)?),
                    None => "script.nys".to_string(),
                };
                Ok((src.to_string(), name))
            }
            (None, Some(p)) => {
                let stem = script_stem(p)?;
                let file = self.scripts_dir().join(format!("{stem}.nys"));
                let text = std::fs::read_to_string(&file).map_err(|_| format!("no saved script `{stem}` (they are the files of {})", self.scripts_dir().display()))?;
                Ok((text, format!("{stem}.nys")))
            }
            (None, None) => Err("give `source` (the script text) or `path` (a saved script)".to_string()),
        }
    }

    fn wait_of(args: &Value) -> Result<Duration, String> {
        match args.get("wait_s") {
            None | Some(Value::Null) => Ok(Duration::from_secs_f64(DEFAULT_WAIT_S)),
            Some(v) => match v.as_f64() {
                Some(w) if (0.0..=MAX_WAIT_S).contains(&w) => Ok(Duration::from_secs_f64(w)),
                _ => Err(format!("`wait_s` is a number of seconds from 0 to {MAX_WAIT_S}")),
            },
        }
    }

    /// The JSON the host gets for one answer of the kernel.
    fn answer(&mut self, r: Result<Step, String>) -> (Value, bool) {
        let (mut v, is_error) = match r {
            Ok(Step::Done(o)) => {
                if o.calls > 0 {
                    self.revision += 1;
                }
                (json!({ "status": "done", "value": cap_value(o.value), "printed": o.printed, "calls": o.calls }), false)
            }
            Ok(Step::Running(p)) => {
                if p.calls > 0 {
                    self.revision += 1;
                }
                (
                    json!({ "status": "running", "waiting_for": p.waiting, "printed": p.printed, "calls": p.calls, "next": "call script_wait to continue, or script_wait with cancel true to abandon it" }),
                    false,
                )
            }
            Err(m) => (json!({ "status": "error", "message": m }), true),
        };
        v["nyx"] = self.header();
        (v, is_error)
    }

    fn script_run(&mut self, args: &Value) -> (Value, bool) {
        let prepared = (|| {
            let (src, name) = self.script_of(args)?;
            let wait = Self::wait_of(args)?;
            if let (Some(_), Some(save)) = (args.get("source"), args.get("save_as").and_then(Value::as_str)) {
                let stem = script_stem(save)?;
                std::fs::create_dir_all(self.scripts_dir()).map_err(|e| format!("cannot create the scripts folder: {e}"))?;
                std::fs::write(self.scripts_dir().join(format!("{stem}.nys")), &src).map_err(|e| format!("cannot save the script: {e}"))?;
            }
            Ok::<_, String>((src, name, wait))
        })();
        let (src, name, wait) = match prepared {
            Ok(x) => x,
            Err(m) => return self.answer(Err(m)),
        };
        let params = args.get("params").cloned().filter(|p| !p.is_null());
        if params.as_ref().is_some_and(|p| !p.is_object()) {
            return self.answer(Err("`params` is an object (a dict)".to_string()));
        }
        let allow = args.get("allow_destructive").and_then(Value::as_bool).unwrap_or(false);
        let r = match self.kernel() {
            // the kernel and the engine are two fields of this struct
            Ok(_) => {
                let (engine, kernel) = (&mut self.engine, self.kernel.as_mut().expect("made above"));
                kernel.run(engine, &name, &src, params.as_ref(), Some(wait), allow)
            }
            Err(m) => Err(m),
        };
        self.answer(r)
    }

    fn script_wait(&mut self, args: &Value) -> (Value, bool) {
        let wait = match Self::wait_of(args) {
            Ok(w) => w,
            Err(m) => return self.answer(Err(m)),
        };
        let cancel = args.get("cancel").and_then(Value::as_bool).unwrap_or(false);
        let Some(kernel) = self.kernel.as_mut() else { return self.answer(Err("no script is running".to_string())) };
        if cancel {
            return match kernel.cancel(&mut self.engine) {
                Ok(()) => {
                    let mut v = json!({ "status": "cancelled" });
                    v["nyx"] = self.header();
                    (v, false)
                }
                Err(m) => self.answer(Err(m)),
            };
        }
        let r = kernel.wait(&mut self.engine, Some(wait));
        self.answer(r)
    }

    fn script_check(&mut self, args: &Value) -> (Value, bool) {
        let (src, name) = match self.script_of(args) {
            Ok(x) => x,
            Err(m) => return self.answer(Err(m)),
        };
        let r = match self.kernel() {
            Ok(_) => {
                let (engine, kernel) = (&mut self.engine, self.kernel.as_mut().expect("made above"));
                kernel.check(engine, &name, &src)
            }
            Err(m) => Err(m),
        };
        match r {
            Ok(()) => {
                let mut v = json!({ "status": "ok" });
                v["nyx"] = self.header();
                (v, false)
            }
            Err(m) => self.answer(Err(m)),
        }
    }

    /// Plain text: a guide card, the list of commands, or one command.
    fn api(&self, args: &Value) -> Result<String, String> {
        if let Some(c) = args.get("command").and_then(Value::as_str) {
            return stubs::command_stub(&self.engine, c).ok_or_else(|| format!("no command `{c}` (api topic \"commands\" lists them)"));
        }
        match args.get("topic").and_then(Value::as_str).unwrap_or("nyx-instructions") {
            "commands" => Ok(stubs::command_list(&self.engine)),
            t => guide::card(t).map(str::to_string).ok_or_else(|| format!("no topic `{t}`: nyx-instructions, geometry, optical-modes, troubleshooting, commands")),
        }
    }

    fn tools_call(&mut self, params: &Value) -> Result<Value, (i64, String)> {
        let name = params.get("name").and_then(Value::as_str).ok_or((-32602, "tools/call needs a string 'name'".to_string()))?.to_string();
        let args = match params.get("arguments") {
            None | Some(Value::Null) => json!({}),
            Some(a) => a.clone(),
        };
        if !args.is_object() {
            return Err((-32602, "arguments is an object".to_string()));
        }
        let out = catch_unwind(AssertUnwindSafe(|| match name.as_str() {
            "script_run" => Some(self.script_run(&args)),
            "script_wait" => Some(self.script_wait(&args)),
            "script_check" => Some(self.script_check(&args)),
            "api" => Some(match self.api(&args) {
                Ok(text) => (json!({ "text": text }), false),
                Err(m) => (json!({ "status": "error", "message": m }), true),
            }),
            _ => None,
        }));
        let Some((value, is_error)) = (match out {
            Ok(x) => x,
            Err(_) => Some((json!({ "status": "error", "message": "nyx hit an internal fault on this call; the session is still usable", "nyx": self.header() }), true)),
        }) else {
            return Err((-32602, format!("unknown tool '{name}'")));
        };
        // one copy of the answer, as text: a guide is plain text, everything else compact JSON
        let text = match (name.as_str(), value.get("text").and_then(Value::as_str)) {
            ("api", Some(t)) if !is_error => t.to_string(),
            _ => value.to_string(),
        };
        Ok(json!({ "content": [{ "type": "text", "text": text }], "isError": is_error }))
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
            "tools/list" => Ok(json!({ "tools": tools() })),
            "tools/call" => self.tools_call(&params),
            m => Err((-32601, format!("method not found: {m}"))),
        };
        Some(match r {
            Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }),
            Err((code, message)) => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } }),
        })
    }
}

fn serve(dir: PathBuf, root: PathBuf) {
    let mut s = Server::new(dir, root);
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
    println!(
        "nyx-mcp {VERSION}\n\nUSAGE\n  nyx-mcp [--dir PATH] [--root PATH]   serve MCP on stdio: the tools script_run, script_wait, script_check and api\n\nWorkflows and saved scripts go to --dir (default {}). The files a script's commands may name must lie in --root (default: --dir when given, else your Documents folder).\nThe program nyx-script-worker must sit next to this one (or be named in NYX_SCRIPT_WORKER).",
        default_dir().display()
    );
}

/// Takes `--flag VALUE` out of the arguments.
fn take_opt(args: &mut Vec<String>, flag: &str) -> Option<PathBuf> {
    let i = args.iter().position(|a| a == flag)?;
    if i + 1 >= args.len() {
        eprintln!("{flag} needs a path");
        std::process::exit(2);
    }
    let v = PathBuf::from(args.remove(i + 1));
    args.remove(i);
    Some(v)
}

fn main() {
    nyx_engine::process::opt_out_of_power_throttling();
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    // started by an engine as its worker (`--job FOLDER`, see `nyx_engine::worker`): solve and answer on stdout, never serve MCP
    if args.first().map(String::as_str) == Some("--job") {
        std::process::exit(nyx_engine::worker::serve(&args));
    }
    let dir_arg = take_opt(&mut args, "--dir");
    let root_arg = take_opt(&mut args, "--root");
    let dir_given = dir_arg.is_some();
    let dir = dir_arg.unwrap_or_else(default_dir);
    let root = root_arg.unwrap_or_else(|| default_root(dir_given, &dir));
    match args.first().map(String::as_str) {
        None => serve(dir, root),
        Some("--version") => println!("{NAME} {VERSION}"),
        Some("--help") | Some("-h") => usage(),
        Some(other) => {
            eprintln!("unknown argument '{other}'");
            usage();
            std::process::exit(2);
        }
    }
}
