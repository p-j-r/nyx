//! The Connect menu: registers the `nyx-mcp` server (the AI interface in `ai/`) with Claude Code and Claude Desktop in one click.
//!
//! Both clients keep their MCP servers in a JSON file under `mcpServers`; connecting adds a `nyx` entry that runs
//! `nyx-mcp`, disconnecting removes it, and every other key of the file is left as it was. The first time a file is
//! changed a copy is kept next to it as `<name>.nyx-backup`. Claude Code reads `~/.claude.json` (what
//! `claude mcp add --scope user` writes) and starts new servers with its next session; Claude Desktop reads
//! `claude_desktop_config.json` and has to be restarted.

use serde_json::{json, Map, Value};
use std::path::{Path, PathBuf};

/// The name the server is registered under in both clients.
pub const SERVER: &str = "nyx";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Client {
    ClaudeCode,
    ClaudeDesktop,
}

impl Client {
    pub const ALL: [Client; 2] = [Client::ClaudeCode, Client::ClaudeDesktop];

    pub fn name(self) -> &'static str {
        match self {
            Client::ClaudeCode => "Claude Code",
            Client::ClaudeDesktop => "Claude Desktop",
        }
    }

    pub fn after_connect(self) -> &'static str {
        match self {
            Client::ClaudeCode => "Start a new Claude Code session to use it.",
            Client::ClaudeDesktop => "Restart Claude Desktop to use it.",
        }
    }

    /// The config file of the client on this computer (it may not exist yet).
    pub fn config_path(self) -> Option<PathBuf> {
        let home = std::env::var_os("USERPROFILE").or_else(|| std::env::var_os("HOME")).map(PathBuf::from);
        match self {
            Client::ClaudeCode => home.map(|h| h.join(".claude.json")),
            Client::ClaudeDesktop => desktop_config_dir(home).map(|d| d.join("claude_desktop_config.json")),
        }
    }
}

/// The folder of `claude_desktop_config.json`. The Store build of Claude Desktop sees `%APPDATA%` through its package
/// (`%LOCALAPPDATA%\Packages\Claude_*\LocalCache\Roaming`), so that folder is used when the package is installed.
fn desktop_config_dir(home: Option<PathBuf>) -> Option<PathBuf> {
    if let Some(local) = std::env::var_os("LOCALAPPDATA") {
        if let Ok(rd) = std::fs::read_dir(Path::new(&local).join("Packages")) {
            for e in rd.flatten() {
                if e.file_name().to_string_lossy().starts_with("Claude_") {
                    return Some(e.path().join("LocalCache").join("Roaming").join("Claude"));
                }
            }
        }
    }
    if let Some(a) = std::env::var_os("APPDATA") {
        return Some(PathBuf::from(a).join("Claude"));
    }
    home.map(|h| if cfg!(target_os = "macos") { h.join("Library").join("Application Support").join("Claude") } else { h.join(".config").join("Claude") })
}

/// Looks for `nyx-mcp` next to the running program and in the `ai/target` folder of the source tree the program was built from.
pub fn find_server() -> Option<PathBuf> {
    let exe = if cfg!(windows) { "nyx-mcp.exe" } else { "nyx-mcp" };
    let mut dirs: Vec<PathBuf> = vec![];
    if let Ok(me) = std::env::current_exe() {
        for a in me.ancestors().skip(1).take(5) {
            dirs.push(a.to_path_buf());
            dirs.push(a.join("ai").join("target").join("release"));
        }
    }
    dirs.push(Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..").join("ai").join("target").join("release"));
    dirs.into_iter().map(|d| d.join(exe)).find(|p| p.is_file()).and_then(|p| p.canonicalize().ok()).map(plain_path)
}

/// Drops the `\\?\` prefix `canonicalize` puts on Windows paths (the clients start the program by this path).
pub fn plain_path(p: PathBuf) -> PathBuf {
    match p.to_string_lossy().strip_prefix(r"\\?\") {
        Some(rest) => PathBuf::from(rest),
        None => p,
    }
}

/// The `mcpServers` entry that runs `server`.
pub fn entry(server: &Path) -> Value {
    json!({ "type": "stdio", "command": server.display().to_string(), "args": [], "env": {} })
}

/// The JSON a user can paste into a client's config by hand.
pub fn snippet(server: &Path) -> String {
    let mut servers = Map::new();
    servers.insert(SERVER.into(), json!({ "command": server.display().to_string() }));
    serde_json::to_string_pretty(&json!({ "mcpServers": servers })).unwrap_or_default()
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Status {
    /// No config file, or no `nyx` entry in it.
    NotConnected,
    /// The entry runs this program.
    Connected,
    /// There is a `nyx` entry that runs something else (a moved or older build).
    OtherPath(String),
    /// The file could not be read as JSON; it is left alone.
    Unreadable(String),
}

fn read_config(path: &Path) -> Result<Value, String> {
    match std::fs::read_to_string(path) {
        Ok(t) if t.trim().is_empty() => Ok(json!({})),
        Ok(t) => {
            let v: Value = serde_json::from_str(t.trim_start_matches('\u{feff}')).map_err(|e| format!("{} is not valid JSON ({e}); not changed", path.display()))?;
            if v.is_object() { Ok(v) } else { Err(format!("{} does not hold a JSON object; not changed", path.display())) }
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(json!({})),
        Err(e) => Err(format!("cannot read {}: {e}", path.display())),
    }
}

pub fn status(path: &Path, server: &Path) -> Status {
    let cfg = match read_config(path) {
        Ok(c) => c,
        Err(e) => return Status::Unreadable(e),
    };
    match cfg.get("mcpServers").and_then(|m| m.get(SERVER)) {
        None => Status::NotConnected,
        Some(e) => match e.get("command").and_then(|c| c.as_str()) {
            Some(c) if Path::new(c) == server => Status::Connected,
            Some(c) => Status::OtherPath(c.to_string()),
            None => Status::OtherPath(String::new()),
        },
    }
}

/// Writes the file through a temporary file in the same folder, so a failure cannot leave half a config.
fn write_config(path: &Path, cfg: &Value) -> Result<(), String> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    if path.is_file() {
        let mut name = path.file_name().unwrap_or_default().to_os_string();
        name.push(".nyx-backup");
        let backup = path.with_file_name(name);
        if !backup.exists() {
            std::fs::copy(path, &backup).map_err(|e| format!("cannot back up {}: {e}", path.display()))?;
        }
    }
    let mut tmp_name = path.file_name().unwrap_or_default().to_os_string();
    tmp_name.push(".nyx-tmp");
    let tmp = path.with_file_name(tmp_name);
    let text = serde_json::to_string_pretty(cfg).map_err(|e| e.to_string())?;
    std::fs::write(&tmp, text).map_err(|e| format!("cannot write {}: {e}", tmp.display()))?;
    std::fs::rename(&tmp, path).map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        format!("cannot replace {}: {e}", path.display())
    })
}

/// Adds (or points at `server`) the `nyx` entry, keeping every other key of the file.
pub fn connect(path: &Path, server: &Path) -> Result<(), String> {
    let mut cfg = read_config(path)?;
    let root = cfg.as_object_mut().ok_or("config is not a JSON object")?;
    let servers = root.entry("mcpServers").or_insert_with(|| json!({}));
    let servers = servers.as_object_mut().ok_or_else(|| format!("\"mcpServers\" in {} is not an object; not changed", path.display()))?;
    servers.insert(SERVER.into(), entry(server));
    write_config(path, &cfg)
}

/// Removes the `nyx` entry; nothing is written if there is none.
pub fn disconnect(path: &Path) -> Result<(), String> {
    let mut cfg = read_config(path)?;
    let removed = cfg.get_mut("mcpServers").and_then(|m| m.as_object_mut()).and_then(|m| m.remove(SERVER)).is_some();
    if removed { write_config(path, &cfg) } else { Ok(()) }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nyx_connect_{}_{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d.join("config.json")
    }

    #[test]
    fn connecting_keeps_the_other_keys_and_servers_and_backs_up_once() {
        let p = temp("keep");
        std::fs::write(&p, r#"{"zeta":1,"mcpServers":{"paper":{"command":"p"}},"alpha":[1,2]}"#).unwrap();
        let exe = PathBuf::from("C:\\x\\nyx-mcp.exe");
        assert_eq!(status(&p, &exe), Status::NotConnected);
        connect(&p, &exe).unwrap();
        assert_eq!(status(&p, &exe), Status::Connected);
        let v: Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        assert_eq!(v["zeta"], 1);
        assert_eq!(v["alpha"], json!([1, 2]));
        assert_eq!(v["mcpServers"]["paper"]["command"], "p");
        assert_eq!(v["mcpServers"]["nyx"]["command"], "C:\\x\\nyx-mcp.exe");
        // the keys keep their order
        let keys: Vec<_> = v.as_object().unwrap().keys().cloned().collect();
        assert_eq!(keys, ["zeta", "mcpServers", "alpha"]);
        let backup = p.with_file_name("config.json.nyx-backup");
        assert!(std::fs::read_to_string(&backup).unwrap().contains("\"zeta\":1"));
        connect(&p, &exe).unwrap();
        assert!(std::fs::read_to_string(&backup).unwrap().contains("\"zeta\":1"), "the backup is the original, not a later state");
        disconnect(&p).unwrap();
        assert_eq!(status(&p, &exe), Status::NotConnected);
        let v: Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        assert_eq!(v["mcpServers"]["paper"]["command"], "p");
    }

    #[test]
    fn a_missing_file_is_created_and_a_broken_one_is_left_alone() {
        let p = temp("new");
        connect(&p, Path::new("/opt/nyx-mcp")).unwrap();
        assert_eq!(status(&p, Path::new("/opt/nyx-mcp")), Status::Connected);
        assert!(matches!(status(&p, Path::new("/elsewhere")), Status::OtherPath(_)));

        let q = temp("bad");
        std::fs::write(&q, "{ not json").unwrap();
        assert!(connect(&q, Path::new("/opt/nyx-mcp")).is_err());
        assert_eq!(std::fs::read_to_string(&q).unwrap(), "{ not json");
        assert!(matches!(status(&q, Path::new("/opt/nyx-mcp")), Status::Unreadable(_)));
    }
}
