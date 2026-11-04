# 37. The engine core is built before the kernels, with two seams

Date: 2026-11-04

## Status

Accepted (restructure, bookmark `restructure-4e-engine-core`). Starts phase 7 early, for the part that does not depend on phases 2 to 4. Supersedes ADR 23 for the project and the jobs (the app still has its own until it is rewired).

## Context

The plan puts `nyx-engine` (project, commands, job runner, store wiring) in phase 7, after the kernels move (phase 4). Checking the plan and the code showed that almost all of it does not touch the kernels: a project is a `geom` feature list, mesh settings and one `ProblemSpec` envelope per physics (all defined by the foundation); a solve is `Physics::solve` behind the registry; the jobs need only `femcore::progress`; persistence is `nyx-store`. Only two calls touch work that is in flight: the **mesher** (`femcore::meshing`, `mesher::extrude_y`: phase 2 replaces the mesh types and moves those modules) and the **solve** itself (`NotYetMoved` until phase 4).

## Decision

* `nyx-engine` holds `Project` (title, feature list, `MeshSettings`, one `ProblemSpec` per physics, the active physics), `Engine` (the project plus the mesh, the results, the jobs and the store binding), `commands` (every action by name, JSON in and out, with a JSON Schema per command and `Engine::catalog()`), `jobs` (a thread per job with `Progress`, cancel, `info`, `try_take`, `wait`; a panic becomes a failed job), `solve` (physics in coupling order, imports filled from the exports of the physics before) and `store_link` (persistence, card, metrics).
* **Seam 1, meshing.** `Mesher` is a trait; `FemcoreMesher` is the only code of the engine that names femcore's mesher. A `MeshedModel` is a mesh plus the names (regions, boundaries) and the `LabelResolver` the physics use; names and labels come from `geom`'s `LabelMap` (domains numbered from 1), not from femcore. Phase 2 rewrites `FemcoreMesher::mesh` and nothing else. `physics_sizing` answers `Unsupported` until the neutral `Sizing` exists.
* **Seam 2, solving** is `Physics::solve`. With the real registry the answer is `NotYetMoved`, which the engine reports (and records as a failed run); the engine tests use stand-in physics registered in a test registry, so the whole path (mesh, resolve, couple, cancel, record) is tested before the kernels move.
* Mesh and results are caches keyed by content (hash of geometry and mesh settings; hash of the project). An edit never invalidates anything by hand; a result whose stamp is not the project's is reported `stale`.
* The project JSON layout is provisional (gate G6: derived in phase 8). It has no `format` key; `to_json` writes sorted keys so the text and hash are the same in every build.
* The window and the MCP server are not rewired here (that is phase 7 proper): `app` and `nyx-session` are untouched and still build.

## Consequences

* Phase 2 and 3 sessions do not conflict with this work: the files they own are not touched; the one coupling is `FemcoreMesher`.
* Phase 4 fills `solve` without touching the engine; phase 7 proper rewires `app` and `nyx-mcp` to `Engine::call` and dissolves `nyx-session`.
* A problem that names a region or selector the model lacks fails with the list of valid names.

## Code

`crates/nyx-engine/src/{project,engine,commands,jobs,meshing,solve,store_link,geometry}.rs`, `crates/nyx-engine/tests/{engine,contracts}.rs`.

## Contract

`engine_commands.v1.json` (the catalog), `engine_project.v1.json` (an example project; provisional, version 0.1 until phase 8).
