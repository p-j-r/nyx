+++
title = "Install & Host Setup"
description = "Pure-Rust nyx-setup installer, Cargo release builds on Rust 1.99.1, headless batch CLI, and MCP host configuration via nyx-connect."
template = "install.html"
+++

A complete `nyx` installation consists of six side-by-side executables built from the workspace: `nyx` (desktop application and headless CLI), `nyx-worker` (out-of-process solver worker), `nyx-script-worker` (sandboxed Monty `.nys` script runner), `nyx-mcp` (MCP server), `nyx-connect` (AI host configurator), and `nyx-uninstall`.

## Pure-Rust Installer (`nyx-setup`) & Release Bundles

Cross-platform release bundles (`.zip`, self-extracting `nyx-setup`, `.nyxpack`, and macOS `.dmg`) are built by `nyx-pack` (`crates/nyx-install`). `nyx-setup` supports both interactive installation and non-interactive terminal automation:

```sh
nyx-setup install --user --yes                 # current-user install, PATH & file associations
nyx-setup install --prefix /opt/nyx --cli -y   # custom prefix without GUI wizard
nyx-setup status --json                        # machine-readable install_report.v1 status
nyx-setup doctor                               # verify binary layout, PATH, and MCP links
nyx-setup uninstall --purge-cache              # removes binaries & temp cache; preserves Documents/nyx
```

## Building from Source (Rust `1.99.1`)

The toolchain is pinned to Rust `1.99.1` in `rust-toolchain.toml`. HDF5 results files (`results.h5`, feature `hdf5` enabled by default) are read and written via pure-Rust `hdf5-pure` and `zlib-rs` (ADR 107), so **no external `HDF5_DIR`, C HDF5 library, or `hdf5.dll` is needed**. Always build in `--release` mode:

```sh
cargo build --release -p nyx -p nyx-script -p nyx-mcp -p nyx-connect
cargo run --release -p nyx
```

To run the test suite and local VVUQ verification gate:

```sh
cargo test --release --workspace
cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
bash tools/vvuq_ci.sh --tier 1          # Windows PowerShell: .\tools\vvuq_ci.ps1 -Tier 1
```

## Headless Batch & Script CLI

`nyx` exposes deterministic headless workflows using the exact same engine as the GUI:

```sh
nyx --batch model.nyx.json [--physics ems,sms]
nyx check script.nys [--json]
nyx run script.nys [--workspace DIR] [--seconds 120] [--memory 4096] [--calls 500]
nyx api list
nyx api show figure
nyx queue submit WORKSPACE PROJECT.nyx.json ems --priority 1
```

When running scripts (`nyx run FILE.nys`), all file I/O is confined strictly to `--workspace` (or the script's parent directory if omitted). On Windows, `nyx` opts out of OS power-throttling efficiency-core demotion (`nyx_engine::process::opt_out_of_power_throttling`) so background solves retain full performance-core clock speeds.

## Connecting AI Hosts (`nyx-connect` & `nyx-mcp`)

`nyx-mcp` exposes four tools over stdio JSON-RPC (`script_run`, `script_wait`, `script_check`, `api`). Wire it into Antigravity, Claude Code, Claude Desktop, Cursor, VS Code, or Windsurf using either the **Connect** menu in the GUI or `nyx-connect`:

```sh
nyx-connect list
nyx-connect connect all          # or: antigravity, claude-code, claude-desktop, cursor, vscode, windsurf
nyx-connect test
nyx-connect snippet              # print raw mcpServers JSON
nyx-connect instructions         # print setup prompt for an AI assistant
```

### Latency & Timeout Guidance (`docs/HOST_SETUP.md`)

Experience from instrumented MCP sessions shows that **over 90% of wall time in unoptimised sessions goes to model turns and tool timeouts rather than solving**. `nyx-mcp` therefore enforces four protections:

1. **Bounded Blocking (`--max-wait`):** Calls never block longer than `90 s` by default (`--max-wait SECONDS` in `1..900` inside the host's `args` array). If a host imposes a shorter timeout, set `"args": ["--max-wait", "60"]`.
2. **Pre-emptible Waits & Progress Heartbeats:** Sending `notifications/cancelled`, `ping`, or a follow-up `script_wait` interrupts a blocked wait cleanly without killing the worker job, and `_meta.progressToken` requests receive `notifications/progress` every `5 s`.
3. **One `script_run` per Task:** Combine geometry construction, `mesh_build`, `solve_run`, reductions (`result_info`, `analyse`, `confinement`, `probe`), and `print` into a single script rather than issuing separate round-trips per command.
4. **Inspect `cost` Before Refining:** `mesh_build` reports predicted time and memory from local solve history, and `convergence_check` supports `dry_run = true`, `max_cells`, and `max_seconds` to reject unaffordable meshes before solving.
