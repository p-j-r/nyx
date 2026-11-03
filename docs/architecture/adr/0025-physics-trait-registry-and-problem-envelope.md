# 25. A physics is a crate behind the `Physics` trait; the problem is an envelope; the registry is built from cargo features

Date: 2026-11-04

## Status

Accepted by paul on 2026-11-03 as the target structure (restructure plan, section 2); crate names are gate G3 (default: as written). Implemented as an interface in the foundation phase; the kernels move in phase 4. Supersedes in part ADR 4 (`SimplexMesh` is still the solver contract until phase 2, but physics take it through the alias `Mesh`) and ADR 8 (the ported solver library is wired in phase 3).

## Context

Adding a physics meant touching the shared `Physics` enum, the shared `ProblemSpec` (every physics-specific option lived in it: `out_of_plane`, `thermal`, `plane`, `brillouin`, ...), `solve::dispatch`, the window tabs, the MCP tools, the mesher's sizing and the store's concept card (ROUGH_SPOTS R10, R11, R12, R15). paul's goals 1, 2, 7 and 10: a new physics must not touch the rest of the codebase; window and MCP have the same capabilities; a project loads one to three physics; physics are compiled in.

## Decision

* `nyx-physics` holds the interface and nothing else of any physics: the `Physics` trait, `Ctx` (progress, cancel, imported fields, a `LabelResolver`), the `ProblemSpec` envelope, the descriptors, `SolveResult`, `Error` (with `Cancelled`), the `Registry`, structural validation and generated JSON Schemas.
* A physics is a crate: `nyx-dfs`, `nyx-sms`, `nyx-ems` (three-letter codes, ADR 12). It owns the typed form of its study, condition and option parameters, its kernels (assemble, solve, fields stay public, ADR 9) and its descriptors.
* **The envelope** is `{ physics, order, study: { kind, params }, materials: { region -> material }, bcs: [{ on, kind, params }], params }`. `params` objects are `serde_json::Value` at the interface and are validated against the descriptors (`validate_spec`): unknown keys, wrong types and out-of-range values are errors with plain messages. There is no `format` key until phase 8.
* **Names, not numbers.** Regions and boundaries are `geom` names and selectors (ADR 2); the engine supplies a `LabelResolver` in `Ctx` that turns them into mesh labels (numbered from 1 after phase 2). A material is a name or, until phase 6, inline scalar properties by descriptor name (`MaterialAssign`).
* **Descriptors** (`studies`, `bc_kinds`, `properties`, `params`, `fields`, `exports`, `imports`, `sizing`) are plain serde data; the window builds forms from them, the MCP server builds its schemas from them, the mesher reads the sizing hint. Deviation from the plan's trait: a `params()` method describes the physics-level options of the envelope (plane assumption, thickness, inactive regions); the plan listed `params` in the envelope but no descriptor for it.
* **Couplings** are one-way field links: a physics `exports` fields under a `role` and `imports` a role; `Registry::solve_order` orders a chain and refuses a cycle. Two-way coupling is later.
* **Only plain data crosses the trait** (flat arrays for fields and meshes), so a separate-process plug-in remains possible.
* **The registry is built by the engine from cargo features** (`dfs`, `sms`, `ems`, default all). Nothing outside the physics crates names a physics (checked by the layering test, ADR 26).
* `Mesh` is `femcore::mesh::Mesh`, an alias of `SimplexMesh` today. Phase 2 changes what it is; code outside femcore names it only through `nyx_physics::Mesh`.
* Until phase 4, the physics crates answer `solve` with `Error::NotYetMoved`. Their descriptors are a first version read from the current kernels and `ProblemSpec`, to be refined by the phase 4 owner under the version rules of ADR 19.

## Consequences

* A new physics is one crate plus one line in `nyx-engine`. The window and MCP server list it without code.
* `Ctx` still holds femcore's `Progress`, which kernels reach through a thread-local (R19); phase 4 makes it explicit.
* The sizing hint is `fn sizing(&self, spec)`; to compute wavelengths it needs the materials, which are names until phase 6 resolves them. The phase 4 or 6 owner may add a materials argument; that is a version bump of `physics_api`.
* `SolveResult` is serialisable (R27 resolved at the interface); it does not carry an analysis (ADR 21 is superseded in part by phase 5).
* Descriptor content for the three physics is a first version and not authoritative until phase 4.

## Code

`crates/nyx-physics/**`, `crates/nyx-dfs/**`, `crates/nyx-sms/**`, `crates/nyx-ems/**`, `crates/nyx-engine/src/lib.rs`, `crates/femcore/src/mesh.rs` (the `Mesh` alias).

## Contract

`problem_spec.v2.json`, `solve_result.v2.json`, `physics_descriptors.v1.json`, `physics_dfs.v1.json`, `physics_sms.v1.json`, `physics_ems.v1.json`; the trait surface is pinned by `crates/nyx-physics/tests/contracts.rs` (`the_physics_trait_is_object_safe_and_has_these_methods`).
