// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-mcp: a Model Context Protocol server over stdio for nyx. Newline-delimited JSON-RPC 2.0 on stdin/stdout, logs on stderr.
//! No async runtime and no SDK: the protocol surface nyx needs is five methods. Since S7 (ADR 94) the server offers four tools, not one per engine
//! command: `script_run`, `script_wait`, `script_check` and `api`. The server holds one `nyx_engine::Engine` and one script kernel
//! (`nyx_script::kernel`, created at the first script, so `initialize` and `tools/list` start no process): the agent writes Python that calls the
//! engine commands and only what the script prints or returns comes back. What a host does is stored as a workflow of the workspace folder
//! (origin `ai`), so the person sees it in the window's dashboard, with the script that made each run.
//!
//! Reading and waiting (L1, ADR 97): a reader thread takes the lines of stdin, so the server sees a host's newer request or a cancellation while
//! a script is waiting for a solve. The script itself runs in its own worker process and the engine and the kernel stay on the main thread:
//! a wait is made of slices of 200 ms (`Kernel::wait` with a short wait), and between two slices the server looks at what arrived. A `script_run`
//! or `script_wait` that arrives answers the waiting call `running` at once (the host that sent it gave up on it) and is then served; a
//! `notifications/cancelled` for the waiting call stops the wait and sends no answer; `ping` and `tools/list` are answered during a wait. A wait is
//! never longer than `--max-wait` (default 90 s) whatever the call asked: a host's tool timeout is often 180 s or less.
use std::collections::{HashSet, VecDeque};
use std::io::{BufRead, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::{Duration, Instant};

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
/// The longest wait unless `--max-wait` says otherwise. Hosts end a tool call after 60 to 180 s (Antigravity: 180 s); an answer `running` is cheap,
/// a call the host gave up on is not (L1).
const MAX_WAIT_S: f64 = 90.0;
/// The most `--max-wait` may name.
const LIMIT_WAIT_S: f64 = 900.0;
/// A wait is made of slices of this length; between two the server looks at what the host sent.
const SLICE: Duration = Duration::from_millis(200);
/// A call that carries a progress token gets a `notifications/progress` this often while it waits (some hosts restart their timeout on it).
const PROGRESS_EVERY: Duration = Duration::from_secs(5);
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

fn tools(max_wait: f64) -> Value {
    let default = DEFAULT_WAIT_S.min(max_wait);
    let wait = json!({ "type": "number", "maximum": max_wait, "description": format!("seconds to let the script run before answering `running` (default {default}, at most {max_wait}: a longer value is cut to it; script_wait continues, so keep it under your host's tool timeout)") });
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

/// What the reader thread took from stdin.
enum Incoming {
    Msg(Value),
    /// a line that is not JSON: answered with a parse error
    Bad(String),
}

/// Takes the lines of stdin on a thread of its own, so the main loop can look at them while a script waits.
fn spawn_reader() -> mpsc::Receiver<Incoming> {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        for line in std::io::stdin().lock().lines() {
            let Ok(line) = line else { break };
            if line.trim().is_empty() {
                continue;
            }
            let m = match serde_json::from_str::<Value>(&line) {
                Ok(v) => Incoming::Msg(v),
                Err(e) => Incoming::Bad(e.to_string()),
            };
            if tx.send(m).is_err() {
                break;
            }
        }
    });
    rx
}

/// One line on stdout. The reader of the lines and the main loop both answer, so a line is written under the lock of stdout. `false`: the host went away.
fn send(v: &Value) -> bool {
    let out = std::io::stdout();
    let mut o = out.lock();
    writeln!(o, "{v}").and_then(|_| o.flush()).is_ok()
}

/// The text a request id (or the `requestId` of a cancellation) is matched by.
fn key_of(v: &Value) -> Option<String> {
    (!v.is_null()).then(|| v.to_string())
}

/// What `drive` starts with.
enum Go {
    Run { name: String, src: String, params: Option<Value>, allow: bool },
    Wait,
}

/// What the host sent while a call waited.
#[derive(Default)]
struct Pumped {
    /// a `script_run` or `script_wait` is next in line: the host has given up on the call that waits
    newer: bool,
    /// the host cancelled the call that waits
    cancelled: bool,
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
    /// the longest wait of a call (`--max-wait`)
    max_wait: f64,
    /// what the reader thread took, the messages that wait their turn, the requests the host cancelled before their turn, and whether stdin ended
    rx: mpsc::Receiver<Incoming>,
    backlog: VecDeque<Incoming>,
    cancelled: HashSet<String>,
    eof: bool,
    /// the `tools/call` being served: its id, its progress token, and the notes its answer carries
    current: Option<String>,
    token: Option<Value>,
    note: Option<String>,
    /// the host cancelled the call being served: it gets no answer
    suppress: bool,
}

impl Server {
    fn new(dir: PathBuf, root: PathBuf, queue: PathBuf, max_wait: f64, rx: mpsc::Receiver<Incoming>) -> Server {
        let mut engine = Engine::new();
        engine.set_workspace(dir.clone(), true);
        engine.set_owner("ai"); // its requests queue behind a person's at the same priority (ADR 87)
        engine.set_self_as_worker(true); // solves run in a child process: running out of memory fails the job, not the server (ADR 86)
        engine.set_history_file(&dir.join("solve_history.jsonl")); // what solves cost, for memory estimates
        // one queue with the window: solves are run by the queue worker of that folder, in order, whoever asked (Q3, ADR 95)
        engine.set_shared_queue(&queue, true);
        Server {
            engine,
            kernel: None,
            dir,
            root,
            proto: PROTOCOLS[0],
            session: session_id(),
            revision: 0,
            max_wait,
            rx,
            backlog: VecDeque::new(),
            cancelled: HashSet::new(),
            eof: false,
            current: None,
            token: None,
            note: None,
            suppress: false,
        }
    }

    /// Takes in what the reader thread has for the server and answers what needs no kernel at once (`ping`, `tools/list`); the rest waits its turn
    /// in `backlog`. `waiting` is the id of the call that waits, if one does.
    fn pump(&mut self, waiting: Option<&str>) -> Pumped {
        let mut out = Pumped::default();
        loop {
            let m = match self.rx.try_recv() {
                Ok(m) => m,
                Err(mpsc::TryRecvError::Empty) => return out,
                Err(mpsc::TryRecvError::Disconnected) => {
                    self.eof = true;
                    return out;
                }
            };
            let Incoming::Msg(v) = &m else {
                self.backlog.push_back(m);
                continue;
            };
            let method = v.get("method").and_then(Value::as_str).unwrap_or("");
            let has_id = v.get("id").is_some_and(|i| !i.is_null());
            match method {
                "ping" | "tools/list" if has_id => {
                    if let Some(r) = self.handle(v) {
                        let _ = send(&r);
                    }
                }
                "notifications/cancelled" => {
                    if let Some(k) = v.get("params").and_then(|p| p.get("requestId")).and_then(key_of) {
                        if waiting == Some(k.as_str()) {
                            out.cancelled = true;
                        } else {
                            self.cancelled.insert(k);
                        }
                    }
                }
                "tools/call" => {
                    let tool = v.get("params").and_then(|p| p.get("name")).and_then(Value::as_str).unwrap_or("");
                    out.newer |= matches!(tool, "script_run" | "script_wait");
                    self.backlog.push_back(m);
                }
                _ => self.backlog.push_back(m),
            }
        }
    }

    /// The next message to serve, waiting for one; `None` when stdin has ended and nothing is left. A request the host cancelled while it waited
    /// its turn is dropped (the protocol sends no answer for it).
    fn next_message(&mut self) -> Option<Incoming> {
        loop {
            self.pump(None);
            let m = match self.backlog.pop_front() {
                Some(m) => m,
                None if self.eof => return None,
                None => match self.rx.recv() {
                    Ok(m) => {
                        self.backlog.push_back(m);
                        continue;
                    }
                    Err(_) => return None,
                },
            };
            if let Incoming::Msg(v) = &m {
                if v.get("id").and_then(key_of).is_some_and(|k| self.cancelled.remove(&k)) {
                    continue;
                }
            }
            return Some(m);
        }
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

    /// The wait a call asked for, cut to the longest wait of the server; the second is the value that was asked when it was cut.
    fn wait_of(&self, args: &Value) -> Result<(Duration, Option<f64>), String> {
        match args.get("wait_s") {
            None | Some(Value::Null) => Ok((Duration::from_secs_f64(DEFAULT_WAIT_S.min(self.max_wait)), None)),
            Some(v) => match v.as_f64() {
                Some(w) if w >= 0.0 && w.is_finite() => {
                    let used = w.min(self.max_wait);
                    Ok((Duration::from_secs_f64(used), (w > self.max_wait).then_some(w)))
                }
                _ => Err(format!("`wait_s` is a number of seconds from 0 to {} (a longer wait is cut to it)", self.max_wait)),
            },
        }
    }

    /// The answer of a call that may wait: the note of a wait that was cut or ended early goes with it.
    fn answer_wait(&mut self, r: Result<Step, String>, cut: Option<f64>) -> (Value, bool) {
        let (mut v, is_error) = self.answer(r);
        if let Some(asked) = cut {
            v["wait_s"] = json!({ "asked": asked, "used": self.max_wait, "note": "a wait is cut to the server's longest; script_wait continues" });
        }
        (v, is_error)
    }

    /// Runs a script or continues the one that waits, for at most `wait`, in slices of `SLICE`. Between two slices the host's messages are looked at
    /// (`pump`): a newer `script_run` or `script_wait` or a cancellation ends the wait with the script still running (the job goes on in the
    /// engine's queue, the next `script_wait` picks it up). What was printed in every slice is kept, and a command is counted once.
    fn drive(&mut self, go: Go, wait: Duration) -> Result<Step, String> {
        let started = Instant::now();
        let deadline = started + wait;
        let mut last_note = started;
        let mut printed: Vec<String> = vec![];
        let slice = wait.min(SLICE);
        let mut r = {
            self.kernel()?;
            let (engine, kernel) = (&mut self.engine, self.kernel.as_mut().expect("made above"));
            match go {
                Go::Run { name, src, params, allow } => kernel.run(engine, &name, &src, params.as_ref(), Some(slice), allow),
                Go::Wait => kernel.wait(engine, Some(slice)),
            }
        };
        loop {
            match r {
                Ok(Step::Running(mut p)) => {
                    printed.append(&mut p.printed);
                    let now = Instant::now();
                    let current = self.current.clone();
                    let got = self.pump(current.as_deref());
                    if got.cancelled {
                        self.suppress = true;
                    }
                    if now >= deadline || got.newer || got.cancelled || self.eof {
                        if got.newer {
                            self.note = Some("a newer request arrived: this wait ended; the script is still running".to_string());
                        }
                        p.printed = printed;
                        return Ok(Step::Running(p));
                    }
                    if now.duration_since(last_note) >= PROGRESS_EVERY {
                        last_note = now;
                        if let Some(t) = &self.token {
                            let n = json!({ "jsonrpc": "2.0", "method": "notifications/progress", "params": { "progressToken": t, "progress": now.duration_since(started).as_secs_f64(), "message": format!("nyx waits for {}", p.waiting) } });
                            let _ = send(&n);
                        }
                    }
                    let slice = SLICE.min(deadline.saturating_duration_since(now));
                    let (engine, kernel) = (&mut self.engine, self.kernel.as_mut().expect("the kernel exists while a script waits"));
                    r = kernel.wait(engine, Some(slice));
                }
                Ok(Step::Done(mut o)) => {
                    printed.append(&mut o.printed);
                    o.printed = printed;
                    return Ok(Step::Done(o));
                }
                Err(m) => return Err(m),
            }
        }
    }

    /// The JSON the host gets for one answer of the kernel.
    fn answer(&mut self, r: Result<Step, String>) -> (Value, bool) {
        let note = self.note.take();
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
                let mut v = json!({ "status": "running", "waiting_for": p.waiting, "printed": p.printed, "calls": p.calls, "next": "call script_wait to continue, or script_wait with cancel true to abandon it" });
                if let Some(n) = note {
                    v["note"] = json!(n);
                }
                // what the engine says of the job: its stage, how long it has run, how much of it is done and what such a solve took on this
                // machine before (no number when nothing was learned): enough to decide between waiting and cancelling
                let entry = self.kernel.as_ref().and_then(Kernel::waiting_job).and_then(|id| self.engine.queue_entry(id));
                if let Some(e) = entry {
                    let mut job = serde_json::Map::new();
                    for k in ["id", "state", "stage", "stage_seconds", "seconds", "fraction", "position", "queued_seconds", "estimate"] {
                        if let Some(x) = e.get(k).filter(|x| !x.is_null()) {
                            job.insert(k.to_string(), x.clone());
                        }
                    }
                    v["job"] = Value::Object(job);
                }
                (v, false)
            }
            Err(m) => (json!({ "status": "error", "message": m }), true),
        };
        v["nyx"] = self.header();
        (v, is_error)
    }

    fn script_run(&mut self, args: &Value) -> (Value, bool) {
        let prepared = (|| {
            let (src, name) = self.script_of(args)?;
            let wait = self.wait_of(args)?;
            if let (Some(_), Some(save)) = (args.get("source"), args.get("save_as").and_then(Value::as_str)) {
                let stem = script_stem(save)?;
                std::fs::create_dir_all(self.scripts_dir()).map_err(|e| format!("cannot create the scripts folder: {e}"))?;
                std::fs::write(self.scripts_dir().join(format!("{stem}.nys")), &src).map_err(|e| format!("cannot save the script: {e}"))?;
            }
            Ok::<_, String>((src, name, wait))
        })();
        let (src, name, (wait, cut)) = match prepared {
            Ok(x) => x,
            Err(m) => return self.answer(Err(m)),
        };
        let params = args.get("params").cloned().filter(|p| !p.is_null());
        if params.as_ref().is_some_and(|p| !p.is_object()) {
            return self.answer(Err("`params` is an object (a dict)".to_string()));
        }
        let allow = args.get("allow_destructive").and_then(Value::as_bool).unwrap_or(false);
        let r = self.drive(Go::Run { name, src, params, allow }, wait);
        self.answer_wait(r, cut)
    }

    fn script_wait(&mut self, args: &Value) -> (Value, bool) {
        let (wait, cut) = match self.wait_of(args) {
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
        if kernel.waiting().is_none() {
            return self.answer(Err("no script is running".to_string()));
        }
        let r = self.drive(Go::Wait, wait);
        self.answer_wait(r, cut)
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
            "tools/list" => Ok(json!({ "tools": tools(self.max_wait) })),
            "tools/call" => {
                // the call being served: what a wait looks for in the host's messages, and what it reports progress to
                self.current = key_of(&id);
                self.token = params.get("_meta").and_then(|m| m.get("progressToken")).cloned();
                self.suppress = false;
                let r = self.tools_call(&params);
                self.current = None;
                self.token = None;
                self.note = None;
                if std::mem::take(&mut self.suppress) {
                    return None; // the host cancelled it: the protocol sends no answer
                }
                r
            }
            m => Err((-32601, format!("method not found: {m}"))),
        };
        Some(match r {
            Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }),
            Err((code, message)) => json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } }),
        })
    }
}

fn serve(dir: PathBuf, root: PathBuf, queue: PathBuf, max_wait: f64) {
    let mut s = Server::new(dir, root, queue, max_wait, spawn_reader());
    eprintln!("{NAME} {VERSION} ready");
    while let Some(incoming) = s.next_message() {
        let reply = match incoming {
            Incoming::Msg(Value::Array(batch)) => {
                // 2025-06-18 dropped batching; answering one is harmless
                let out: Vec<Value> = batch.iter().filter_map(|m| s.handle(m)).collect();
                if out.is_empty() { None } else { Some(Value::Array(out)) }
            }
            Incoming::Msg(m) => s.handle(&m),
            Incoming::Bad(e) => Some(json!({ "jsonrpc": "2.0", "id": null, "error": { "code": -32700, "message": format!("parse error: {e}") } })),
        };
        if let Some(r) = reply {
            if !send(&r) {
                break; // the host went away
            }
        }
    }
    eprintln!("{NAME} stdin closed, exiting");
}

fn usage() {
    println!(
        "nyx-mcp {VERSION}\n\nUSAGE\n  nyx-mcp [--dir PATH] [--root PATH] [--max-wait SECONDS]   serve MCP on stdio: the tools script_run, script_wait, script_check and api\n\nWorkflows and saved scripts go to --dir (default {}). The files a script's commands may name must lie in --root (default: --dir when given, else your Documents folder).\nA call never waits longer than --max-wait seconds (default {MAX_WAIT_S}, at most {LIMIT_WAIT_S}) whatever its wait_s says: keep it under your host's tool timeout; script_wait continues.\nThe program nyx-script-worker must sit next to this one (or be named in NYX_SCRIPT_WORKER).",
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
    if matches!(args.first().map(String::as_str), Some("--job" | "--worker")) {
        std::process::exit(nyx_engine::worker::serve(&args));
    }
    let dir_arg = take_opt(&mut args, "--dir");
    let root_arg = take_opt(&mut args, "--root");
    let max_wait = match take_opt(&mut args, "--max-wait") {
        None => MAX_WAIT_S,
        Some(v) => match v.to_str().and_then(|s| s.parse::<f64>().ok()).filter(|w| (1.0..=LIMIT_WAIT_S).contains(w)) {
            Some(w) => w,
            None => {
                eprintln!("--max-wait is a number of seconds from 1 to {LIMIT_WAIT_S}");
                std::process::exit(2);
            }
        },
    };
    let dir_given = dir_arg.is_some();
    let dir = dir_arg.unwrap_or_else(default_dir);
    let root = root_arg.unwrap_or_else(|| default_root(dir_given, &dir));
    match args.first().map(String::as_str) {
        None => {
            // the queue is the window's (`NYX_WORKSPACE`, else Documents/nyx) unless `--dir` named a folder: then the folder has its own
            let queue = if dir_given { dir.clone() } else { nyx_engine::shared_queue::default_workspace() };
            serve(dir, root, queue, max_wait)
        }
        Some("--version") => println!("{NAME} {VERSION}"),
        Some("--help") | Some("-h") => usage(),
        Some(other) => {
            eprintln!("unknown argument '{other}'");
            usage();
            std::process::exit(2);
        }
    }
}
