//! Connects the `nyx-mcp` server to the AI hosts installed on this computer, the way a "Connect" button in a design tool does:
//! find each host's MCP config file, add or remove one `nyx` entry, keep a backup, and never overwrite a file it cannot parse.
//!
//! No GUI and no heavy dependencies, so the nyx app can call it from a settings panel and the `nyx-connect` command can call it from a shell.
//! Everything that depends on the machine goes through [`Env`], so the logic is tested against temporary folders.
use serde_json::{json, Map, Value};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// The name of the entry written into every host's config.
pub const ENTRY: &str = "nyx";

/// Where things live on this machine. [`Env::detect`] reads the process environment; tests build one by hand.
#[derive(Clone, Debug)]
pub struct Env {
    pub home: PathBuf,
    /// `%APPDATA%` on Windows.
    pub appdata: Option<PathBuf>,
    /// `%LOCALAPPDATA%` on Windows.
    pub local_appdata: Option<PathBuf>,
    pub os: Os,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Os {
    Windows,
    Mac,
    Linux,
}

impl Env {
    pub fn detect() -> Result<Env, String> {
        let var = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let os = if cfg!(windows) { Os::Windows } else if cfg!(target_os = "macos") { Os::Mac } else { Os::Linux };
        let home = var("USERPROFILE").or_else(|| var("HOME")).ok_or("cannot find the home folder (USERPROFILE or HOME)")?;
        Ok(Env { home, appdata: var("APPDATA"), local_appdata: var("LOCALAPPDATA"), os })
    }
}

/// How a host stores MCP servers in its JSON file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shape {
    /// `{"mcpServers": {"nyx": {"command": ..., "args": []}}}`
    McpServers,
    /// Claude Code additionally records `"type": "stdio"`.
    McpServersTyped,
    /// VS Code: `{"servers": {"nyx": {"type": "stdio", "command": ...}}}`.
    VsCode,
}

impl Shape {
    fn key(self) -> &'static str {
        match self {
            Shape::VsCode => "servers",
            _ => "mcpServers",
        }
    }
    fn entry(self, exe: &Path) -> Value {
        let cmd = exe.to_string_lossy();
        match self {
            Shape::McpServers => json!({ "command": cmd, "args": [] }),
            Shape::McpServersTyped | Shape::VsCode => json!({ "type": "stdio", "command": cmd, "args": [] }),
        }
    }
}

/// A host nyx knows how to connect to.
#[derive(Clone, Debug)]
pub struct Host {
    pub id: &'static str,
    pub name: &'static str,
    /// The config file nyx writes. It may not exist yet.
    pub config: PathBuf,
    /// A folder whose presence says the host is installed (its config folder).
    marker: PathBuf,
    shape: Shape,
}

/// All hosts, with their paths for this machine.
pub fn hosts(env: &Env) -> Vec<Host> {
    let h = &env.home;
    let mut v = vec![];
    v.push(Host { id: "claude-code", name: "Claude Code", config: h.join(".claude.json"), marker: h.join(".claude"), shape: Shape::McpServersTyped });

    // Claude Desktop: the installer and the Microsoft Store package keep their data in different places.
    let desktop_dir = match env.os {
        Os::Windows => store_package_dir(env).or_else(|| env.appdata.clone().map(|a| a.join("Claude"))).unwrap_or_else(|| h.join("AppData").join("Roaming").join("Claude")),
        Os::Mac => h.join("Library").join("Application Support").join("Claude"),
        Os::Linux => h.join(".config").join("Claude"),
    };
    v.push(Host { id: "claude-desktop", name: "Claude Desktop", config: desktop_dir.join("claude_desktop_config.json"), marker: desktop_dir, shape: Shape::McpServers });

    v.push(Host { id: "cursor", name: "Cursor", config: h.join(".cursor").join("mcp.json"), marker: h.join(".cursor"), shape: Shape::McpServers });

    let code_user = match env.os {
        Os::Windows => env.appdata.clone().unwrap_or_else(|| h.join("AppData").join("Roaming")).join("Code").join("User"),
        Os::Mac => h.join("Library").join("Application Support").join("Code").join("User"),
        Os::Linux => h.join(".config").join("Code").join("User"),
    };
    v.push(Host { id: "vscode", name: "VS Code", config: code_user.join("mcp.json"), marker: code_user, shape: Shape::VsCode });

    let ws = h.join(".codeium").join("windsurf");
    v.push(Host { id: "windsurf", name: "Windsurf", config: ws.join("mcp_config.json"), marker: ws, shape: Shape::McpServers });
    v
}

/// The Store build of Claude Desktop keeps its config under `%LOCALAPPDATA%\Packages\Claude_*\LocalCache\Roaming\Claude`.
fn store_package_dir(env: &Env) -> Option<PathBuf> {
    let packages = env.local_appdata.as_ref()?.join("Packages");
    let mut found: Vec<PathBuf> = fs::read_dir(packages)
        .ok()?
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with("Claude_"))
        .map(|e| e.path().join("LocalCache").join("Roaming").join("Claude"))
        .filter(|p| p.is_dir())
        .collect();
    found.sort();
    found.into_iter().next()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Status {
    /// The host's config folder does not exist.
    NotInstalled,
    /// The host is there, nyx is not in its config.
    NotConnected,
    /// nyx is in the config and points at this exe.
    Connected,
    /// nyx is in the config but points somewhere else (an older build, another folder).
    Stale { found: String },
    /// The config exists but is not plain JSON (comments, a syntax error), so nyx will not touch it.
    Unreadable { why: String },
}

#[derive(Clone, Debug)]
pub struct HostInfo {
    pub id: &'static str,
    pub name: &'static str,
    pub config: PathBuf,
    pub status: Status,
}

pub fn status(host: &Host, exe: &Path) -> Status {
    let cfg = match read_config(&host.config) {
        Ok(c) => c,
        Err(why) => return Status::Unreadable { why },
    };
    let entry = cfg.as_ref().and_then(|c| c.get(host.shape.key())).and_then(|s| s.get(ENTRY));
    match entry {
        Some(e) => {
            let found = e.get("command").and_then(Value::as_str).unwrap_or("").to_string();
            if same_path(&found, exe) { Status::Connected } else { Status::Stale { found } }
        }
        None if cfg.is_some() || host.marker.is_dir() => Status::NotConnected,
        None => Status::NotInstalled,
    }
}

/// Every host with its status, for a settings panel.
pub fn list(env: &Env, exe: &Path) -> Vec<HostInfo> {
    hosts(env).into_iter().map(|h| HostInfo { id: h.id, name: h.name, status: status(&h, exe), config: h.config.clone() }).collect()
}

fn same_path(a: &str, b: &Path) -> bool {
    let b = b.to_string_lossy();
    if cfg!(windows) { a.eq_ignore_ascii_case(&b) } else { a == b }
}

/// `Ok(None)`: no file or an empty file. `Err`: present but not plain JSON, or not an object.
fn read_config(path: &Path) -> Result<Option<Value>, String> {
    let text = match fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(format!("cannot read {}: {e}", path.display())),
    };
    let text = text.trim_start_matches('\u{feff}');
    if text.trim().is_empty() {
        return Ok(None);
    }
    match serde_json::from_str::<Value>(text) {
        Ok(v) if v.is_object() => Ok(Some(v)),
        Ok(_) => Err("the file is JSON but not an object".to_string()),
        Err(e) => Err(format!("not plain JSON ({e}); it may contain comments. Add the entry by hand (see the snippet)")),
    }
}

#[derive(Clone, Debug)]
pub struct Done {
    pub host: &'static str,
    pub config: PathBuf,
    /// The copy of the previous file, when there was one.
    pub backup: Option<PathBuf>,
    pub changed: bool,
}

fn find<'a>(hs: &'a [Host], id: &str) -> Result<&'a Host, String> {
    hs.iter().find(|h| h.id == id).ok_or_else(|| format!("no host '{id}'; valid: {}", hs.iter().map(|h| h.id).collect::<Vec<_>>().join(", ")))
}

/// Adds or updates the `nyx` entry. Other servers and settings in the file are kept, in their order.
pub fn connect(env: &Env, id: &str, exe: &Path) -> Result<Done, String> {
    let hs = hosts(env);
    let host = find(&hs, id)?;
    if !exe.is_file() {
        return Err(format!("the nyx-mcp program is not at {}; build it first (cd ai && cargo build --release -p nyx-mcp)", exe.display()));
    }
    let mut cfg = read_config(&host.config)?.unwrap_or_else(|| Value::Object(Map::new()));
    let key = host.shape.key();
    let servers = cfg.as_object_mut().expect("checked to be an object").entry(key).or_insert_with(|| Value::Object(Map::new()));
    let servers = servers.as_object_mut().ok_or_else(|| format!("'{key}' in {} is not an object; not touching it", host.config.display()))?;
    let want = host.shape.entry(exe);
    if servers.get(ENTRY) == Some(&want) {
        return Ok(Done { host: host.id, config: host.config.clone(), backup: None, changed: false });
    }
    servers.insert(ENTRY.to_string(), want);
    let backup = write_with_backup(&host.config, &cfg)?;
    Ok(Done { host: host.id, config: host.config.clone(), backup, changed: true })
}

/// Removes the `nyx` entry and leaves everything else.
pub fn disconnect(env: &Env, id: &str) -> Result<Done, String> {
    let hs = hosts(env);
    let host = find(&hs, id)?;
    let Some(mut cfg) = read_config(&host.config)? else {
        return Ok(Done { host: host.id, config: host.config.clone(), backup: None, changed: false });
    };
    let removed = cfg.get_mut(host.shape.key()).and_then(Value::as_object_mut).map(|s| s.shift_remove(ENTRY).is_some()).unwrap_or(false);
    if !removed {
        return Ok(Done { host: host.id, config: host.config.clone(), backup: None, changed: false });
    }
    let backup = write_with_backup(&host.config, &cfg)?;
    Ok(Done { host: host.id, config: host.config.clone(), backup, changed: true })
}

/// Copies the old file to `<name>.nyx-backup-<unix seconds>`, then replaces it by writing a temp file next to it and renaming.
fn write_with_backup(path: &Path, cfg: &Value) -> Result<Option<PathBuf>, String> {
    let dir = path.parent().ok_or("config path has no folder")?;
    fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let backup = if path.is_file() {
        let b = dir.join(format!("{name}.nyx-backup-{secs}"));
        fs::copy(path, &b).map_err(|e| format!("cannot back up {}: {e}", path.display()))?;
        Some(b)
    } else {
        None
    };
    let tmp = dir.join(format!("{name}.nyx-tmp"));
    let text = serde_json::to_string_pretty(cfg).map_err(|e| e.to_string())? + "\n";
    fs::write(&tmp, text).map_err(|e| format!("cannot write {}: {e}", tmp.display()))?;
    fs::rename(&tmp, path).map_err(|e| {
        let _ = fs::remove_file(&tmp);
        format!("cannot replace {}: {e}", path.display())
    })?;
    Ok(backup)
}

/// The JSON a user can paste into any host that is not listed.
pub fn snippet(exe: &Path) -> String {
    serde_json::to_string_pretty(&json!({ "mcpServers": { ENTRY: Shape::McpServers.entry(exe) } })).unwrap_or_default()
}

/// The `claude mcp add` line for Claude Code, for people who prefer the command.
pub fn claude_code_command(exe: &Path) -> String {
    format!("claude mcp add {ENTRY} -- \"{}\"", exe.display())
}

/// Where `nyx-mcp` should be: next to the running program, or the path in `NYX_MCP`.
pub fn server_exe() -> PathBuf {
    if let Some(p) = std::env::var_os("NYX_MCP").filter(|v| !v.is_empty()) {
        return PathBuf::from(p);
    }
    let name = if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" };
    std::env::current_exe().ok().and_then(|e| e.parent().map(|d| d.join(name))).unwrap_or_else(|| PathBuf::from(name))
}

#[derive(Clone, Debug)]
pub struct ServerInfo {
    pub name: String,
    pub version: String,
    pub protocol: String,
    pub tools: Vec<String>,
}

/// Starts the server like a host would, asks `initialize` and `tools/list`, and stops it.
pub fn test_server(exe: &Path) -> Result<ServerInfo, String> {
    let mut child = Command::new(exe)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("cannot start {}: {e}", exe.display()))?;
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
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            let left = deadline.saturating_duration_since(Instant::now());
            let line = recv.recv_timeout(left).map_err(|_| format!("no answer to {method} within 10 s"))?;
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

#[cfg(test)]
mod tests {
    use super::*;

    struct Sandbox(PathBuf);
    impl Sandbox {
        fn new(tag: &str) -> Sandbox {
            let p = std::env::temp_dir().join(format!("nyx_connect_{tag}_{}", std::process::id()));
            let _ = fs::remove_dir_all(&p);
            fs::create_dir_all(&p).unwrap();
            Sandbox(p)
        }
        fn env(&self, os: Os) -> Env {
            Env { home: self.0.join("home"), appdata: Some(self.0.join("home/AppData/Roaming")), local_appdata: Some(self.0.join("home/AppData/Local")), os }
        }
        fn exe(&self) -> PathBuf {
            let p = self.0.join(if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" });
            fs::write(&p, "x").unwrap();
            p
        }
    }
    impl Drop for Sandbox {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn host<'a>(hs: &'a [Host], id: &str) -> &'a Host {
        hs.iter().find(|h| h.id == id).unwrap()
    }

    #[test]
    fn a_host_that_is_not_there_is_not_installed_and_a_present_one_is_not_connected() {
        let sb = Sandbox::new("detect");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        assert!(list(&env, &exe).iter().all(|h| h.status == Status::NotInstalled));
        fs::create_dir_all(env.home.join(".cursor")).unwrap();
        let cursor = list(&env, &exe).into_iter().find(|h| h.id == "cursor").unwrap();
        assert_eq!(cursor.status, Status::NotConnected);
    }

    #[test]
    fn connect_keeps_other_servers_and_settings_in_order_and_makes_a_backup() {
        let sb = Sandbox::new("merge");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        let cfg = env.home.join(".claude.json");
        fs::create_dir_all(&env.home).unwrap();
        fs::write(&cfg, r#"{"numStartups": 7, "mcpServers": {"paper": {"type": "http", "url": "http://x"}}, "theme": "dark"}"#).unwrap();
        let done = connect(&env, "claude-code", &exe).unwrap();
        assert!(done.changed);
        let backup = done.backup.expect("a backup of the old file");
        assert!(fs::read_to_string(backup).unwrap().contains("numStartups"));
        let v: Value = serde_json::from_str(&fs::read_to_string(&cfg).unwrap()).unwrap();
        assert_eq!(v["numStartups"], 7);
        assert_eq!(v["theme"], "dark");
        assert_eq!(v["mcpServers"]["paper"]["url"], "http://x");
        assert_eq!(v["mcpServers"]["nyx"]["command"], json!(exe.to_string_lossy()));
        assert_eq!(v["mcpServers"]["nyx"]["type"], "stdio");
        let keys: Vec<&String> = v.as_object().unwrap().keys().collect();
        assert_eq!(keys, ["numStartups", "mcpServers", "theme"]);
        let order: Vec<&String> = v["mcpServers"].as_object().unwrap().keys().collect();
        assert_eq!(order, ["paper", "nyx"]);
        assert_eq!(list(&env, &exe).into_iter().find(|h| h.id == "claude-code").unwrap().status, Status::Connected);
    }

    #[test]
    fn connecting_twice_changes_nothing_the_second_time() {
        let sb = Sandbox::new("twice");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        assert!(connect(&env, "cursor", &exe).unwrap().changed);
        let again = connect(&env, "cursor", &exe).unwrap();
        assert!(!again.changed && again.backup.is_none());
    }

    #[test]
    fn an_entry_for_another_exe_is_stale_and_connect_repairs_it() {
        let sb = Sandbox::new("stale");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        let cfg = env.home.join(".cursor").join("mcp.json");
        fs::create_dir_all(cfg.parent().unwrap()).unwrap();
        fs::write(&cfg, r#"{"mcpServers": {"nyx": {"command": "/old/nyx-mcp", "args": []}}}"#).unwrap();
        let hs = hosts(&env);
        assert!(matches!(status(host(&hs, "cursor"), &exe), Status::Stale { .. }));
        assert!(connect(&env, "cursor", &exe).unwrap().changed);
        assert_eq!(status(host(&hs, "cursor"), &exe), Status::Connected);
    }

    #[test]
    fn disconnect_removes_only_nyx() {
        let sb = Sandbox::new("remove");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        let cfg = env.home.join(".cursor").join("mcp.json");
        fs::create_dir_all(cfg.parent().unwrap()).unwrap();
        fs::write(&cfg, r#"{"mcpServers": {"other": {"command": "o"}}}"#).unwrap();
        connect(&env, "cursor", &exe).unwrap();
        assert!(disconnect(&env, "cursor").unwrap().changed);
        let v: Value = serde_json::from_str(&fs::read_to_string(&cfg).unwrap()).unwrap();
        assert!(v["mcpServers"].get("nyx").is_none());
        assert_eq!(v["mcpServers"]["other"]["command"], "o");
        assert!(!disconnect(&env, "cursor").unwrap().changed);
    }

    #[test]
    fn a_file_with_comments_is_never_overwritten() {
        let sb = Sandbox::new("jsonc");
        let env = sb.env(Os::Linux);
        let exe = sb.exe();
        let cfg = env.home.join(".config/Code/User/mcp.json");
        fs::create_dir_all(cfg.parent().unwrap()).unwrap();
        let original = "{\n  // my servers\n  \"servers\": {}\n}\n";
        fs::write(&cfg, original).unwrap();
        let hs = hosts(&env);
        assert!(matches!(status(host(&hs, "vscode"), &exe), Status::Unreadable { .. }));
        assert!(connect(&env, "vscode", &exe).is_err());
        assert_eq!(fs::read_to_string(&cfg).unwrap(), original);
    }

    #[test]
    fn vscode_uses_the_servers_key() {
        let sb = Sandbox::new("vscode");
        let env = sb.env(Os::Mac);
        let exe = sb.exe();
        connect(&env, "vscode", &exe).unwrap();
        let cfg = env.home.join("Library/Application Support/Code/User/mcp.json");
        let v: Value = serde_json::from_str(&fs::read_to_string(cfg).unwrap()).unwrap();
        assert_eq!(v["servers"]["nyx"]["type"], "stdio");
        assert!(v.get("mcpServers").is_none());
    }

    #[test]
    fn the_store_build_of_claude_desktop_is_found_on_windows() {
        let sb = Sandbox::new("store");
        let env = sb.env(Os::Windows);
        let pkg = env.local_appdata.clone().unwrap().join("Packages/Claude_pzs8sxrjxfjjc/LocalCache/Roaming/Claude");
        fs::create_dir_all(&pkg).unwrap();
        let hs = hosts(&env);
        assert_eq!(host(&hs, "claude-desktop").config, pkg.join("claude_desktop_config.json"));
    }

    #[test]
    fn connect_refuses_when_the_server_program_is_missing() {
        let sb = Sandbox::new("noexe");
        let env = sb.env(Os::Linux);
        let e = connect(&env, "cursor", &sb.0.join("missing")).unwrap_err();
        assert!(e.contains("build it first"), "{e}");
        assert!(connect(&env, "nope", &sb.0).unwrap_err().contains("valid:"));
    }
}
