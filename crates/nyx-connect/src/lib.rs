// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Connects the `nyx-mcp` server to the AI hosts installed on this computer, the way a "Connect" button in a design tool does:
//! find each host's MCP config file, add or remove one `nyx` entry, keep a backup, and never overwrite a file it cannot parse.
//!
//! This is the only copy of that logic (ADR 28). Three callers use it: the app's Connect menu, the `nyx-connect` command, and the
//! MCP server's own tests. No GUI and no heavy dependencies (`serde_json` only).
//!
//! Modules:
//! * [`env`]: where things live on this machine ([`Env`]); tests build one by hand, so nothing here reads the real home folder unless asked.
//! * [`hosts`]: the table of known hosts and the shape of the entry each one wants.
//! * [`config`]: [`status`], [`connect`], [`disconnect`], [`list`]: reading and rewriting a host's JSON file.
//! * [`server`]: finding the `nyx-mcp` program and probing it with a real MCP handshake ([`test_server`]).
//! * [`report`]: the machine-readable form (`nyx-connect ... --json`) and its JSON Schema, versioned (`connect_report.v1`).
//!
//! The contract with the outside world is in `docs/architecture/contracts/connect_*.v1.*`, pinned by `tests/contracts.rs`.
pub mod config;
pub mod env;
pub mod hosts;
pub mod report;
pub mod server;

pub use config::{connect, disconnect, list, status, Done, HostInfo, Status};
pub use env::{Env, Os};
pub use hosts::{hosts, Host};
pub use server::{agent_instructions, claude_code_command, find_server, is_cargo_build_tree, find_server_with, plain_path, server_exe, snippet, test_server, test_server_with, test_server_within, ServerInfo};

/// The name of the entry written into every host's config.
pub const ENTRY: &str = "nyx";
