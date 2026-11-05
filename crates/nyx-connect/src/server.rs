//! The `nyx-mcp` program: where it is, what to paste for a host that is not in the table, and a real MCP handshake against it.
use serde_json::{json, Value};
use std::ffi::OsString;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

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
    if let Some(me) = running {
        if let Some(own) = me.parent() {
            dirs.push(own.to_path_buf());
        }
        if is_cargo_build_tree(me) {
            for a in me.ancestors().skip(1).take(5) {
                dirs.push(a.to_path_buf());
                dirs.push(a.join("target").join("release"));
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
