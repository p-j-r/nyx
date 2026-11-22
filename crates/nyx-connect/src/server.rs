// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `nyx-mcp` program: where it is, what to paste for a host that is not in the table, and a real MCP handshake against it.
use serde_json::{json, Value};
use std::ffi::OsString;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crate::env::Env;
use crate::hosts::hosts;
use crate::ENTRY;

/// The file name of the server program on this operating system.
fn exe_name() -> &'static str {
    if cfg!(windows) {
        "nyx-mcp.exe"
    } else {
        "nyx-mcp"
    }
}

/// Drops the `\\?\` prefix `canonicalize` puts on Windows paths (the hosts start the program by this path).
pub fn plain_path(p: PathBuf) -> PathBuf {
    match p.to_string_lossy().strip_prefix(r"\\?\") {
        Some(rest) => PathBuf::from(rest),
        None => p,
    }
}

/// Looks for `nyx-mcp` in the path `NYX_MCP` and next to the running program; in a cargo build tree also around it (see [`find_server_with`]).
pub fn find_server() -> Option<PathBuf> {
    find_server_with(std::env::var_os("NYX_MCP"), std::env::current_exe().ok().as_deref())
}

/// True when `exe` was built by cargo in place: its folder holds cargo's `deps` and `.fingerprint` folders. A shipped program is never in one.
pub fn is_cargo_build_tree(exe: &Path) -> bool {
    exe.parent().is_some_and(|d| d.join("deps").is_dir() && d.join(".fingerprint").is_dir())
}

/// [`find_server`] with its two inputs given: the value of `NYX_MCP` and the path of the running program. The first existing file wins:
/// 1. the override (if it is set and wrong, nothing is found; there is no fallback);
/// 2. the folder of `running`: this is the whole lookup of a shipped program (the caller then offers "Locate nyx-mcp...");
/// 3. only when `running` is in a cargo build tree ([`is_cargo_build_tree`], so developer builds of any profile, `--release` included):
///    each of the five folders made of the program's own folder and its four parents, and for each of them its `target/release`.
pub fn find_server_with(override_path: Option<OsString>, running: Option<&Path>) -> Option<PathBuf> {
    if let Some(p) = override_path.filter(|v| !v.is_empty()).map(PathBuf::from) {
        return p.is_file().then(|| p.canonicalize().map(plain_path).unwrap_or(p));
    }
    let mut dirs: Vec<PathBuf> = vec![];
    let mut push_dir = |d: PathBuf| {
        if !dirs.contains(&d) {
            dirs.push(d);
        }
    };
    if let Some(me) = running {
        let resolved = me.canonicalize().ok().map(plain_path);
        for cand in [Some(me), resolved.as_deref()] {
            let Some(p) = cand else { continue };
            if let Some(own) = p.parent() {
                push_dir(own.to_path_buf());
            }
            if is_cargo_build_tree(p) {
                for a in p.ancestors().skip(1).take(5) {
                    push_dir(a.to_path_buf());
                    push_dir(a.join("target").join("release"));
                }
            }
        }
    }
    dirs.into_iter().map(|d| d.join(exe_name())).find(|p| p.is_file()).and_then(|p| p.canonicalize().ok()).map(plain_path)
}

/// Where `nyx-mcp` should be, found or not: the path in `NYX_MCP` as given, else what [`find_server`] finds, else next to the running
/// program. For a command line that wants to say "not found" about a definite path.
pub fn server_exe() -> PathBuf {
    if let Some(p) = std::env::var_os("NYX_MCP").filter(|v| !v.is_empty()) {
        return PathBuf::from(p);
    }
    find_server().unwrap_or_else(|| std::env::current_exe().ok().and_then(|e| e.parent().map(|d| d.join(exe_name()))).unwrap_or_else(|| PathBuf::from(exe_name())))
}

/// The JSON a user can paste into any host that is not listed.
pub fn snippet(exe: &Path) -> String {
    serde_json::to_string_pretty(&json!({ "mcpServers": { ENTRY: { "command": exe.to_string_lossy(), "args": [] } } })).unwrap_or_default()
}

/// Text to paste into an AI app that can edit files (Antigravity, Cursor, Claude Code ...) so that it connects nyx itself: which file, which
/// key, which entry, and the rules `connect` follows (a backup first, valid JSON, only the `nyx` entry touched). The paths are this machine's.
pub fn agent_instructions(env: &Env, exe: &Path) -> String {
    let plain = serde_json::to_string(&json!({ "command": exe.to_string_lossy(), "args": [] })).unwrap_or_default();
    let mut s = format!(
        "Connect the nyx MCP server to the app you are running in. Edit that app's MCP configuration file and change nothing else.\n\n\
         The server program is: {exe}\n\
         The entry to add under the servers key, named \"{ENTRY}\": {plain}\n\n\
         Which file (use the one of the app you run in; create it if it is missing):\n",
        exe = exe.display()
    );
    for h in hosts(env) {
        let seen = if h.config.is_file() || h.marker.is_dir() { "" } else { " (not found on this computer)" };
        let entry = h.entry(exe);
        let own = if entry == serde_json::from_str::<Value>(&plain).unwrap_or(Value::Null) { String::new() } else { format!(", with this entry: {entry}") };
        s.push_str(&format!("- {}: {}, under the key \"{}\"{own}{seen}\n", h.name, h.config.display(), h.servers_key()));
    }
    s.push_str(&format!(
        "\nRules:\n\
         1. Keep every other server and setting in the file, in their order.\n\
         2. If there is no copy of the file named like it with \".nyx-backup\" added, make one before you change it.\n\
         3. The file must stay plain JSON. If it has comments or does not parse, stop and tell me; do not rewrite it.\n\
         4. If a \"{ENTRY}\" entry exists, change only its \"command\" to the program above and keep its other keys.\n\
         5. Then tell me to restart the app (or refresh its MCP servers). When it is back, nyx offers four tools: script_run, script_wait, script_check and api.\n"
    ));
    s
}

/// The `claude mcp add` line for Claude Code, for people who prefer the command.
pub fn claude_code_command(exe: &Path) -> String {
    format!("claude mcp add {ENTRY} -- \"{}\"", exe.display())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ServerInfo {
    pub name: String,
    pub version: String,
    pub protocol: String,
    pub tools: Vec<String>,
}

/// Starts the server like a host would, asks `initialize` and `tools/list`, and stops it. Each answer may take 10 s.
pub fn test_server(exe: &Path) -> Result<ServerInfo, String> {
    test_server_within(exe, Duration::from_secs(10))
}

/// [`test_server`] with the time one answer may take.
pub fn test_server_within(exe: &Path, per_answer: Duration) -> Result<ServerInfo, String> {
    test_server_with(Command::new(exe), per_answer)
}

/// [`test_server_within`] for a command the caller prepared (arguments, environment). Its standard streams are replaced.
pub fn test_server_with(mut command: Command, per_answer: Duration) -> Result<ServerInfo, String> {
    let what = command.get_program().to_string_lossy().into_owned();
    let mut child = command.stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::null()).spawn().map_err(|e| format!("cannot start {what}: {e}"))?;
    let mut tx = child.stdin.take().ok_or("no stdin")?;
    let rx = child.stdout.take().ok_or("no stdout")?;
    // read on a thread so a hung server cannot hang the caller
    let (send, recv) = std::sync::mpsc::channel::<String>();
    std::thread::spawn(move || {
        for line in BufReader::new(rx).lines().map_while(Result::ok) {
            if send.send(line).is_err() {
                break;
            }
        }
    });
    let mut ask = |id: u64, method: &str, params: Value| -> Result<Value, String> {
        let req = json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params });
        writeln!(tx, "{req}").and_then(|_| tx.flush()).map_err(|e| format!("the server closed its input: {e}"))?;
        let deadline = Instant::now() + per_answer;
        loop {
            let left = deadline.saturating_duration_since(Instant::now());
            let line = recv.recv_timeout(left).map_err(|_| format!("no answer to {method} within {} s", per_answer.as_secs_f64()))?;
            let v: Value = serde_json::from_str(&line).map_err(|e| format!("not JSON from the server: {e}"))?;
            if v["id"] == json!(id) {
                return v.get("result").cloned().ok_or_else(|| format!("{method} failed: {}", v["error"]));
            }
        }
    };
    let out = (|| {
        let init = ask(1, "initialize", json!({ "protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": { "name": "nyx-connect", "version": env!("CARGO_PKG_VERSION") } }))?;
        let list = ask(2, "tools/list", json!({}))?;
        Ok(ServerInfo {
            name: init["serverInfo"]["name"].as_str().unwrap_or("").to_string(),
            version: init["serverInfo"]["version"].as_str().unwrap_or("").to_string(),
            protocol: init["protocolVersion"].as_str().unwrap_or("").to_string(),
            tools: list["tools"].as_array().into_iter().flatten().filter_map(|t| t["name"].as_str().map(String::from)).collect(),
        })
    })();
    let _ = child.kill();
    let _ = child.wait();
    out
}
