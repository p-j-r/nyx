# 15. The AI interface is a separate workspace with a hand-written stdio MCP server

Date: 2026-11-01

## Status

Accepted (`docs/MCP_ASSISTANT_PLAN.md` sections 18 to 20; code brought under version control in `36f6b3d`; Connect menu `839efa1`).

Superseded in part by ADR 26: `ai/` is no longer a separate workspace; its crates are in `crates/` and the root workspace.

Superseded in part by ADR 28: the Connect logic exists once, in `nyx-connect`; the app has no copy (R07 resolved).

## Context

An AI host should be able to build a geometry, mesh it, define a problem, solve and pick a mode, with nyx deciding *how* and the AI deciding *what*. The plan was written under the constraint that `app` and `femcore` stay untouched.

## Decision

* The code lives in `ai/`, **its own Cargo workspace** with its own lockfile and target: `nyx-session` (typed tools, knows nothing of MCP), `nyx-mcp` (JSON-RPC 2.0 over stdio, newline-delimited, logs on stderr), `nyx-connect` (detect hosts, write or remove the `nyx` entry with a backup, test the server).
* The server is **hand-written**: no async runtime, no SDK. The subset needed (initialize, tools, resources, progress, cancel) is small, and an SDK would add a dependency tree to a repository whose `geom` is serde-only.
* Each tool has a **parameter table** (`params.rs`) from which the JSON Schema is produced, arguments validated, defaults filled and every assumed value recorded in a **ledger**. Parameters carry a `Matters` rating: high ones are asked for when missing (`status: needs_input`), the others are assumed and listed. `use_defaults: true` accepts the stated fallback.
* **Content-hashed handles**: `g_`, `m_`, `p_`, `r_` plus the FNV-1a digest of the inputs, so repeating a request returns the same handle without recomputing. A journal (`journal.jsonl`), the geom JSON of every geometry (`g_*.geom.json`) and a digest of every result are written to the session folder; meshes and results stay in memory.
* Units are SI everywhere. The reference check is oxide-clad 500 x 300 nm on 2 um BOX: TE0 2.6718, TM0 2.3338.
* "One owner per session": a nyx window and `nyx-mcp` must not both own a session. The agent and the window run in one light process; every solve can run in a detached worker so a host exit does not kill a large solve (settled in plan section 20). Local models: detect and guide only.
* The Connect menu in `app` registers `nyx-mcp` with Claude Code (`~/.claude.json`) and Claude Desktop (`claude_desktop_config.json`), keeps a `.nyx-backup` copy and leaves other keys alone.

## Consequences

* `Project` is still defined in `app`, so the session cannot read or write a project file (R16).
* The Connect logic exists twice, in `app/src/connect.rs` and `ai/nyx-connect` (R07); `nyx-connect` also handles Cursor, VS Code and Windsurf.
* `CLAUDE.md` and `ai/README.md` say the workspace depends on `femcore` and `geom`; `nyx-session` depends on `femcore` only (R17).
* The tool list and schemas are an external contract with every host that has registered `nyx` (contract `mcp_tools`).

## Code

`ai/nyx-session/src/{tools,params,workspace,guide}.rs`, `ai/nyx-mcp/src/main.rs`, `ai/nyx-connect/src/*`, `crates/app/src/connect.rs`.

## Contract

`contracts/mcp_tools` (`ai/nyx-mcp/tests/contracts.rs`): the names, access annotations and input-schema property names of the tools, snapshotted.
