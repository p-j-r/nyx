# 65. Closing the restructure: what the final structure replaces, and the rules left

Date: 2026-11-13

## Status

Accepted (phase 9). Brings the older records in line with the code after phases 1 to 8. Supersedes ADR 4, ADR 22 and ADR 23, and supersedes ADR 15 and ADR 14 in part. It decides R23 (formatting, lints) and R14 (stale design notes).

## Context

The decision records of 2026-10-28 to 2026-11-03 describe the structure *before* the restructure. Several are now wrong in what they tell a reader to do, even though they were right when written.

## Decision

**Superseded or amended.**
* **ADR 4 (SimplexMesh is the one mesh the kernels read): superseded** by ADR 5 and ADR 32 and by the physics interface of ADR 25. `Physics::solve` takes the exchange mesh `femcore::mesh::Mesh`; a kernel numbers its own unknowns with `nyx_physics::discretise` (the old `SimplexMesh` numbering lives on as that discretisation and as the view the window and the results file draw on).
* **ADR 5 (native mixed-cell mesh)** is no longer "not implemented on the solver side": the mesh is native from the mesher to the physics interface (ADR 32) and the structured layers are native (ADR 52).
* **ADR 8** was superseded by ADR 33 (phase 3). **ADR 14 (workflow store)** stays in force for the folder layout and the archives; the spec files are `experiment.nyx.json` and `specs/<hash>.nyx.json` (ADR 63), the store takes spec text and a card function (ADR 39) and knows nothing of a project (ADR 59). **ADR 15** (AI interface in its own workspace) is superseded in part by ADR 26 (one workspace), ADR 28 (one `nyx-connect`) and ADR 57 (the MCP tools are the engine commands; `nyx-session` is dissolved); the MCP server itself and the host connection stand.
* **ADR 22 (use cases and chips live in femcore): superseded** by ADR 24 (chips and Brillouin removed) and ADR 60 (examples are project files; the use cases are gone).
* **ADR 23 (the app is a thin shell over one project): superseded** by ADR 57: the project and every action live in `nyx-engine`; the window and the MCP server are front-ends of one `Engine`.

**Formatting and lints (R23).**
* `cargo clippy --workspace --all-targets` is clean. Where a lint is allowed on purpose it is allowed on the item or the module with a reason (`type_complexity` for local tuples of per-cell arrays, `large_enum_variant` for `Stiffness`, which is made once per material).
* `rustfmt` is applied to the crates that already follow `geom`'s settings (`geom`, `nyx-physics`, `nyx-dfs`, `nyx-sms`, `nyx-ems`, `nyx-analysis`, `nyx-engine`) and to nothing else. `femcore`, `nyx-store`, `nyx-mcp`, `nyx-connect` and `app` stay as written: formatting them would rewrite most of their lines (and the line endings of the older files, which are CRLF) for no change in behaviour, and would bury the history that `jj` and `git blame` are used for. New files in those crates follow the style around them.
* `rust-version` in the manifest is 1.99, the toolchain the code is built and tested with (`rust-toolchain.toml`); the old "1.84" was never tested.

**Records.** The design notes of the time before the restructure (`GEOMETRY_ARCHITECTURE.md`, `GEOM_HANDOFF.md`, `SWEEP_CACHING_PROPOSAL.md`, `WORKFLOWS_DESIGN.md`, `MCP_ASSISTANT_PLAN.md`, `KERNEL_PORT_STATUS.md` and others) are history and are not rewritten (R14). The place to find what is true now is `docs/architecture/` (ADRs, `workspace.dsl`, contracts) and `CLAUDE.md`.

**Parked, as paul asked ("we will see at the very end"):** run-time plug-ins, two-way coupling, Brillouin as an optional physics crate (tag `sbs-last`), and the licence check of the BSD-2-Clause ports against nyx's own licence before anything is published.

## Consequences

* A reader who follows an ADR marked superseded is sent to the one that replaced it.
* The guards that keep the final structure are tests: `tests/layering.rs` (crate layering, and geom's dependencies: R18), `femcore/tests/expr_grammar.rs` (R06), the contract suites (ADR 19) and `check_code_to_adr.sh`.

## Code

`docs/architecture/` (this directory, `workspace.dsl`, `diagrams.md`, `CODE_TO_ADR.md`), `CLAUDE.md`, `Cargo.toml` (`rust-version`).
