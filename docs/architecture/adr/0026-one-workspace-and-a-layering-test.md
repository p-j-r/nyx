# 26. One Cargo workspace for everything, and a test that holds the crate layering

Date: 2026-11-04

## Status

Accepted (foundation of the restructure, following the plan's "`ai/` joins the root workspace"). Supersedes ADR 15 in part (the separate `ai/` workspace).

## Context

`ai/` was a separate Cargo workspace "so that nothing in `../crates` is touched" (ADR 15). That gave it its own lockfile and target directory, two `nyx-connect` copies (R07), and a build graph nobody could check as a whole (R18). Separately, `cargo test --workspace` showed that `serde_json`'s `preserve_order` feature, switched on by the app and `nyx-connect`, changes the key order of JSON objects in the test binaries: seventeen golden tests and the MCP tool-list contract compared text or key lists and failed for that reason alone.

## Decision

* `nyx-session`, `nyx-mcp`, `nyx-connect` live in `crates/` and are members of the root workspace; `ai/` is gone. One lockfile, one `cargo test --workspace`.
* `.cargo/config.toml` sets `target-dir = "../nyx2-target"`, a sibling folder, so every jj workspace of the repository (ADR 27) shares the compiled registry dependencies. Cargo serialises builds in one directory; the cost is waiting for the lock, the gain is disk (a target is 3 GB and the machine had 17 GB free).
* `crates/nyx-engine/tests/layering.rs` holds the dependency table of the target structure (plan section 2) and a short `TEMPORARY` list of edges that exist today with the phase that removes each. A new dependency between workspace crates fails the test until the table is changed on purpose; a `TEMPORARY` row whose edge is gone fails too. Nothing but the engine may depend on a physics crate; `femcore` depends on no nyx crate.
* **Tests never compare JSON text or key order.** Goldens compare sorted key lists; contract exemplars are written and compared with `canonical_json` (sorted keys at every depth) or as `Value`s.
* New crates (`nyx-physics`, the three physics, `nyx-analysis`, `nyx-engine`) are formatted with `rustfmt` (the `crates/geom/rustfmt.toml` settings) and clippy-clean with `-D warnings`; the older crates are not yet (R23, decided in phase 9).

## Consequences

* The MCP build is part of `cargo build --workspace`; `ai/target` and `ai/Cargo.lock` are gone.
* The layering table is the single place that shows what the structure is allowed to be; a phase that removes an edge deletes a row.
* A parallel working session that adds a crate adds a row to the table (a one-line conflict at worst).

## Code

`Cargo.toml`, `.cargo/config.toml`, `crates/nyx-engine/tests/layering.rs`, `crates/nyx-session`, `crates/nyx-mcp`, `crates/nyx-connect`, `crates/femcore/tests/golden_frames.rs` (sorted keys), `crates/nyx-mcp/tests/contracts.rs` (values).

## Contract

None of its own; the layering test is the guard.
