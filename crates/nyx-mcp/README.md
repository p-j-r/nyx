# nyx-mcp, nyx-session, nyx-connect: the MCP server for nyx

Three crates of the root workspace (a separate `ai/` workspace until the restructure). Phase 7 of `docs/architecture/RESTRUCTURE_PLAN.md` turns the tools into engine commands and dissolves `nyx-session`.
Design: `docs/MCP_ASSISTANT_PLAN.md`.

| crate | what |
|---|---|
| `nyx-session` | the tools an AI can call (geometry, mesh, problem, solve, mode pick, checks, one-call recipe), handles, assumptions ledger, structured errors, guide cards. Knows nothing about MCP. |
| `nyx-mcp` | the MCP server: JSON-RPC over stdio, no async runtime, no SDK. |
| `nyx-connect` | one-click connection to AI apps (Claude Code, Claude Desktop, Cursor, VS Code, Windsurf): detects each, writes the `nyx` entry with a backup, removes it, tests the server. A library for the nyx settings panel and a command. |

## Build

```
cargo build --release -p nyx-mcp        # target\release\nyx-mcp.exe
cargo test --release -p nyx-session -p nyx-mcp -p nyx-connect   # unit tests + a stdio test with a real solve
```

## Try it without a host

```
nyx-mcp --list
nyx-mcp --call recipe_run "{\"recipe\":\"wg_modes_2d\",\"thickness\":3e-7,\"width\":5e-7,\"cladding\":\"oxide\",\"box\":2e-6}"
```

Units are SI (metres). A missing `width` returns `status: needs_input` with a question; pass `use_defaults: true` to accept the stated fallback.
Reference check: oxide-clad 500 x 300 nm on 2 um BOX gives TE0 2.6718, TM0 2.3338 (the repository documents 2.6718 / 2.3337).

## Connect a host

Easiest: `nyx-connect list`, then `nyx-connect connect claude-code` (or `all`). It keeps a `.nyx-backup-*` copy of any file it changes and refuses files with comments. `nyx-connect test` starts the server and lists its tools. By hand:

Claude Code: `claude mcp add nyx -- C:\path\to\nyx2\target\release\nyx-mcp.exe`

Claude Desktop, Cursor, and others (JSON config):

```json
{ "mcpServers": { "nyx": { "command": "C:\\path\\to\\nyx\\ai\\target\\release\\nyx-mcp.exe", "args": [] } } }
```

Then ask: "what is the TE0 effective index of a 300 nm silicon-on-insulator waveguide?"

## Files it writes

`--dir PATH` (default `Documents\nyx_runs\mcp`), one folder per server run: `journal.jsonl` (every call and outcome), `g_*.geom.json` (the geom `FeatureList`, loadable in nyx), `r_*.result.json`.
Handles are content hashes, so the same inputs give the same handle in any session.

## Not yet

Detached solve worker for big models, 3D, the typed project-file mirror, recipes beyond the SOI strip, the assistant window. See the plan, sections 22.5 and following.
