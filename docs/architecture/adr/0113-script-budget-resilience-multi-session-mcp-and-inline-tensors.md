# 113. Script per-turn suspension resilience, UTF-8 BOM tolerance, inline tensor properties, and multi-session MCP

Date: 2026-11-21

## Status

Accepted. Amends ADR 25, 89, 93, 94, and 107.

## Context

Long multi-turn `.nys` parameter sweeps in `nyx-mcp` could exhaust `monty-pool`'s cumulative session suspension counter even when each individual turn stayed well under `Limits::calls`, or poison the worker after a single turn hit a call/memory ceiling. In addition, Windows editors frequently prepend a UTF-8 BOM (`\u{feff}`) to `.nys`, `.nyx.json`, `.jsonl`, and `.toml` files; specifying diagonal or symmetric permittivity/conductivity tensors required authoring a full TOML material; and concurrent AI agents sharing one `nyx-mcp` server collided on a single `Engine` and `Kernel`.

## Decision

* **Per-turn Monty suspension budget, proactive recycle, and `reset_session` (`nyx-script::kernel`):**
  * Separates per-turn `calls` (reset on every `start()`) from cumulative `session_suspensions`, setting Monty's internal ceiling to `max(16 * calls, 200_000)` and proactively recycling the checkout via `restore(dump)` before rollover so REPL globals survive arbitrarily many turns.
  * Classifies per-turn call overruns, `MemoryError`, suspension/timeout limits, and `PoolError::Finished` as `Fail::Fatal`, discarding the poisoned checkout and restoring the last clean post-turn `dump` on the next turn; adds `Kernel::reset()` and `script_run(reset_session=true)` to clear REPL globals without touching `Engine`.
  * Extends `_nyx_conv` (`prelude.nys`) and `trace::roots` (`trace.rs`) to recurse into Python dicts so batch named-map `expr` lambdas trace and statically guard cleanly.
* **UTF-8 BOM (`\u{feff}`) stripping across all text entry points:**
  * Strips leading `\u{feff}` in `Kernel::{run,check}`, `host::{run,check}`, `nyx-mcp::script_of`, `app::runfile::read`, `commandlist::parse`, `Project::from_json`, and `MaterialSpec::from_toml_named` / `Library::{load_dir,read_table}`.
* **Inline tensor shorthand (`PropValue`, `MaterialAssign::PropsTensor` in `nyx-physics`):**
  * Adds `PropValue::{Scalar(f64), Diag([f64; 3]), Sym6([f64; 6]), Matrix([[f64; 3]; 3])}` and `MaterialAssign::PropsTensor` (keeping `MaterialAssign::Props` intact), validating `PropKind::Tensor`, finiteness, and $3\times 3$ symmetry ($10^{-12}$ rel) and decoding into `Materials::aniso`.
* **Multi-session `nyx-mcp`, `"cqed"` guide card, and portable binary check:**
  * Replaces the singleton `Engine` + `Kernel` in `nyx-mcp` with `BTreeMap<String, McpSession>` keyed by optional `session` (default `"default"`, rejecting empty or path-separator/traversal names, isolating non-default sessions under `dir/sessions/<name>`).
  * Adds embedded guide card `crates/nyx-engine/cards/cqed.md` (`topic: "cqed"` in `get_guide` and `nyx-mcp`'s `api` tool) and a non-fatal dynamic-symbol check in `tools/package_dist.sh` confirming zero external `hdf5.dll` / `libhdf5` linkage (ADR 107).

## Consequences

* Multi-agent workflows can run isolated concurrent sessions over a single `nyx-mcp` process and execute long multi-turn parameter sweeps without session-budget rollovers or BOM parse failures.

## Code

`crates/nyx-script/src/{kernel,host,prelude.nys,trace}.rs`, `crates/nyx-script/tests/resilience.rs`, `crates/nyx-mcp/src/main.rs`, `crates/nyx-mcp/tests/{reset_session,multi_session}.rs`, `crates/nyx-physics/src/{spec,validate,decode,schema,lib}.rs`, `crates/nyx-physics/tests/inline_tensor.rs`, `crates/femcore/src/material/spec.rs`, `crates/nyx-engine/cards/cqed.md`, `crates/nyx-engine/src/{guide,project,commandlist,store_link,commands}.rs`, `crates/app/src/{runfile,app,tabs}.rs`, `tools/package_dist.sh`.
