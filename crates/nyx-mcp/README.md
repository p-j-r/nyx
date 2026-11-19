# nyx-mcp, nyx-connect: the MCP server for nyx

`nyx-mcp` is a stdio MCP server (JSON-RPC 2.0, no async runtime, no SDK) over one `nyx_engine::Engine`. Since S7 of the productivity plan (ADR 94) it has **four tools**,
not one per engine command: the agent writes a short Python script that calls the engine commands, and only what the script prints or returns comes back.
`nyx-connect` connects the server to the AI apps installed on this computer (Claude Code, Claude Desktop, Cursor, VS Code, Windsurf).

| tool | what |
|---|---|
| `api` | read first. No argument: how to work with nyx. `topic`: a guide (`geometry`, `optical-modes`, `troubleshooting`) or `commands` (every command with its title). `command`: one command with its Python signature, documentation and arguments |
| `script_run` | runs a script (`source`, or `path` of a script saved in the workspace; `params` a dict the script reads as `params`; `save_as`; `wait_s`; `allow_destructive`). Answers `{status, value, printed, calls, nyx}`; `status` is `done`, `running` or `error` |
| `script_wait` | continues a script that answered `running` (a solve is still going), or with `cancel` abandons it |
| `script_check` | parses and type-checks a script against the commands and the variables of the session, runs nothing |

A script is Python syntax (`*.nys`) run by Monty in a worker process (`nyx-script-worker`, which must sit next to `nyx-mcp`; `NYX_SCRIPT_WORKER` names it):

```python
nyx.project_open(path="examples/15_soi_wire.nyx.json")
for w in params.get("widths", [4e-7, 5e-7, 6e-7]):
    nyx.set_param(name="w", value=w)
    nyx.solve_run()
    print(w, nyx.result_info(physics="ems")[0]["solutions"][0]["eigenvalue"][0])
```

* Variables and functions stay from one `script_run` to the next. A script saved with `save_as` is run again with `path` and other `params`: nothing is resent.
* No tool blocks for ever: a solve that outlasts `wait_s` leaves the script paused and `script_run` answers `running`; `script_wait` continues it.
* `project_new`, `project_load`, `project_open`, `store_open` and `remove_*` need `allow_destructive: true`.
* The files a script's commands may name must lie in `--root` (default: `--dir` when given, else your Documents folder). Saved scripts are `<dir>/scripts/<name>.nys`.
* The same API is `nyx api stubs` (the type stub an editor reads) and `nyx run FILE.nys` runs a script from the command line.

## Build and test

```
cargo build --release -p nyx-mcp -p nyx-script     # target\release\nyx-mcp.exe and nyx-script-worker.exe
cargo test --release -p nyx-mcp                    # the tests that run a script need the worker built first
```

## Try it without a host

```
nyx-mcp --help
nyx run docs/scripting-samples/w1_soi_width.nys --workspace .
```

## Connect a host

`nyx-connect list`, then `nyx-connect connect claude-code` (or `all`). It keeps a `.nyx-backup-*` copy of any file it changes and refuses files with comments. `nyx-connect test` starts the server and does a real handshake.

Claude Code: `claude mcp add nyx -- C:\path\to\nyx2\target\release\nyx-mcp.exe`

Other hosts (JSON config):

```json
{ "mcpServers": { "nyx": { "command": "C:\path\to\nyx2\target\release\nyx-mcp.exe", "args": [] } } }
```

## Files it writes

`--dir PATH` (default `Documents\nyx_runs\mcp`): the workspace the work is stored in (workflows and their runs, each run with the script that made it), `scripts/`, and the solve history.
