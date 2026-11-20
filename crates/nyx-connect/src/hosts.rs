//! The table of AI hosts nyx can connect to: where each keeps its MCP config file and what an entry in it looks like.
//! Adding a host is one row here and one row in the contract exemplar `connect_hosts.v1.json`.
use crate::env::{Env, Os};
use serde_json::{json, Map, Value};
use std::fs;
use std::path::{Path, PathBuf};

/// How a host stores MCP servers in its JSON file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Shape {
    /// `{"mcpServers": {"nyx": {"command": ..., "args": []}}}`
    McpServers,
    /// Claude Code, as `claude mcp add` writes it: also `"type": "stdio"` and `"env": {}`.
    McpServersTyped,
    /// VS Code: `{"servers": {"nyx": {"type": "stdio", "command": ...}}}`.
    VsCode,
}

impl Shape {
    /// The key of the file's object that holds the servers.
    pub(crate) fn key(self) -> &'static str {
        match self {
            Shape::VsCode => "servers",
            _ => "mcpServers",
        }
    }

    /// The entry for `exe` as a fresh file would get it.
    pub(crate) fn entry(self, exe: &Path) -> Value {
        self.entry_over(None, exe)
    }

    /// The entry for `exe`, built on the one already in the file: the keys nyx owns (`type`, `command`, `args`) are set, every other key the
    /// user added (an `env` with their variables, say) is kept, and `env` is only added to Claude Code's entry when it is missing.
    pub(crate) fn entry_over(self, existing: Option<&Value>, exe: &Path) -> Value {
        let mut m: Map<String, Value> = existing.and_then(Value::as_object).cloned().unwrap_or_default();
        if self != Shape::McpServers {
            m.insert("type".into(), json!("stdio"));
        }
        m.insert("command".into(), json!(exe.to_string_lossy()));
        m.insert("args".into(), json!([]));
        if self == Shape::McpServersTyped {
            m.entry("env").or_insert_with(|| json!({}));
        }
        Value::Object(m)
    }
}

/// A host nyx knows how to connect to.
#[derive(Clone, Debug)]
pub struct Host {
    pub id: &'static str,
    pub name: &'static str,
    /// The config file nyx writes. It may not exist yet.
    pub config: PathBuf,
    /// What to tell the user after connecting.
    pub after_connect: &'static str,
    /// A folder whose presence says the host is installed (its config folder).
    pub(crate) marker: PathBuf,
    pub(crate) shape: Shape,
}

impl Host {
    /// The key of the config file's top-level object under which the servers are kept (`mcpServers`, or `servers` for VS Code).
    pub fn servers_key(&self) -> &'static str {
        self.shape.key()
    }

    /// The entry written for `exe` into an empty config.
    pub fn entry(&self, exe: &Path) -> Value {
        self.shape.entry(exe)
    }
}

/// All hosts, with their paths for this machine, in a fixed order.
pub fn hosts(env: &Env) -> Vec<Host> {
    let h = &env.home;
    let mut v = vec![];
    v.push(Host {
        id: "claude-code",
        name: "Claude Code",
        config: h.join(".claude.json"),
        after_connect: "Start a new Claude Code session to use it.",
        marker: h.join(".claude"),
        shape: Shape::McpServersTyped,
    });

    // Claude Desktop: the installer and the Microsoft Store package keep their data in different places.
    let desktop_dir = match env.os {
        Os::Windows => store_package_dir(env).or_else(|| env.appdata.clone().map(|a| a.join("Claude"))).unwrap_or_else(|| h.join("AppData").join("Roaming").join("Claude")),
        Os::Mac => h.join("Library").join("Application Support").join("Claude"),
        Os::Linux => h.join(".config").join("Claude"),
    };
    v.push(Host {
        id: "claude-desktop",
        name: "Claude Desktop",
        config: desktop_dir.join("claude_desktop_config.json"),
        after_connect: "Restart Claude Desktop to use it.",
        marker: desktop_dir,
        shape: Shape::McpServers,
    });

    v.push(Host {
        id: "cursor",
        name: "Cursor",
        config: h.join(".cursor").join("mcp.json"),
        after_connect: "Restart Cursor to use it.",
        marker: h.join(".cursor"),
        shape: Shape::McpServers,
    });

    let code_user = match env.os {
        Os::Windows => env.appdata.clone().unwrap_or_else(|| h.join("AppData").join("Roaming")).join("Code").join("User"),
        Os::Mac => h.join("Library").join("Application Support").join("Code").join("User"),
        Os::Linux => h.join(".config").join("Code").join("User"),
    };
    v.push(Host {
        id: "vscode",
        name: "VS Code",
        config: code_user.join("mcp.json"),
        after_connect: "Reload the VS Code window to use it.",
        marker: code_user,
        shape: Shape::VsCode,
    });

    let ws = h.join(".codeium").join("windsurf");
    v.push(Host {
        id: "windsurf",
        name: "Windsurf",
        config: ws.join("mcp_config.json"),
        after_connect: "Restart Windsurf to use it.",
        marker: ws,
        shape: Shape::McpServers,
    });

    // Antigravity: one file for the IDE and its command line, `~/.gemini/config/mcp_config.json` (seen on a machine that runs it). The marker is
    // the `config` folder, not `~/.gemini`, which Gemini CLI users have too. A `.agents/mcp_config.json` in a workspace is documented but has
    // been reported as ignored, so it is not used.
    let gemini_config = h.join(".gemini").join("config");
    v.push(Host {
        id: "antigravity",
        name: "Antigravity",
        config: gemini_config.join("mcp_config.json"),
        after_connect: "Restart Antigravity, or refresh its MCP servers, to use it.",
        marker: gemini_config,
        shape: Shape::McpServers,
    });
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
