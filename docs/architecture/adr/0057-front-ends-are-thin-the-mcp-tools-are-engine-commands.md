# 57. The window and the MCP server are thin front-ends: every tool is an engine command

Date: 2026-11-05

## Status

Accepted. Implements restructure phase 7. Supersedes the parts of ADR 15 (AI interface: separate workspace, stdio MCP) that describe the `nyx-session` tools, handles and recipes, and the parts of ADR 23 (the app is a shell over one `Project`) that describe `app/state.rs` and `worker.rs`. ADR 41's legacy bridge is deleted.

## Context

After phases 4 to 6 the physics lived in crates behind the `Physics` trait and the engine core (ADR 37, 38) existed with 64 commands, but the two front-ends still used the pre-restructure stack: the window held its own `Project` of the old `femcore::problem::ProblemSpec` and called `nyx_engine::legacy` for solves, and `nyx-mcp` served nine hand-written tools of `nyx-session` (a silicon strip on oxide: `geom_template`, `mesh_make`, `problem_define`, `solve_run`, `result_modes`, `check_run`, `recipe_run`) over `femcore::usecases`. A new physics, a new geometry operation or a new analysis needed edits in three places, and the AI could do less than the window.

## Decision

* **One engine, two front-ends.** `nyx-engine::Engine` holds the project, the mesh, the results and the jobs. The window holds an `Engine` and edits it with commands (`set_problem`, `set_geometry`, `set_mesh`, `add_condition`, `set_material`, `solve_run` through `start_solve`, ...); `nyx-mcp` holds an `Engine` and nothing else. The tool list of the MCP server is `Engine::catalog()` and a tool call is `Engine::call`; the geometry tools are the `add_*` commands generated from the descriptors of the geom operations (ADR 38), all 23, none picked by hand (gate G5).
* **Deleted:** `nyx-session` (its handles, assumptions ledger and recipes), `nyx_engine::legacy`, `femcore::problem` (the old `ProblemSpec`, `SolveResult`, `Study`, `BcKind`), `femcore::usecases` (the SOI wire and the Bragg grating as built-in devices), `femcore::xmlio`, `femcore::h5io`, the XML project format in the window, the built-in material presets and the 18 built-in examples as code. What the person used to get from a use case is a project file (ADR 60) or a sequence of commands (the guide card `nyx-instructions`).
* **MCP behaviour.** `get_guide` is a command (cards in `nyx-engine/cards`). Every result carries `nyx.{session, revision, project_hash}` (the revision counts calls that changed the project). An error is `{error: <code>, message}` with `isError`. The server stores what a host does: the first solve binds the project to a new workflow of the server's workspace folder (origin `ai`), so every solve is a row of the run table and the work shows in the dashboard (`Engine::set_workspace(dir, auto_store)`). The tool list is a new contract version (`mcp_tools.v2.json`; the nine old tools are gone: a major change for hosts, announced in `VERSIONS.md`).
* **The window draws forms from descriptors** (`app/src/forms.rs`): the Study tab, the conditions of the Physics tab and the options of a physics are made from the `ParamDesc`s of the physics, so a new physics needs no change in the window (ADR 25). The Materials tab lists the regions of the model and the files of the library (ADR 58).
* **Layering.** The app depends on `nyx-engine` and `nyx-connect` only; what it draws with and edits (`ViewMesh`, `Display`, the geometry kernel, the store, the physics result types) comes through `nyx_engine::facade`, so the dependency table has no other edge. `nyx-mcp` depends on `nyx-engine` only. The `TEMPORARY` list of the layering test is empty.
* **Not in this step:** `sapphire()` and `silicon()` in `femcore::material` (still used by tests of other crates), the final project layout and its `format` key (phase 8), HDF5 export (an optional engine export, not ported), renaming `experiment.nyx.xml` (phase 8).

## Consequences

* A new geometry operation in `geom` appears in the engine, the MCP server and (through the geometry editor's own forms) the window with the descriptor; a new physics appears in the Physics, Study and Materials tabs and as `add_physics`/`describe_physics` for a host.
* The AI interface has the capabilities of the window (geometry, materials, conditions, mesh, solve, analysis, store). It no longer has a one-call recipe for a silicon strip; the guide shows the sequence. Mesh-convergence checks are done by repeating `solve_run` with a finer `set_mesh`.
* Tests that used the use cases are small generic fixtures (`femcore::testkit::soi`, the example files).

## Code

`crates/nyx-mcp/src/main.rs`, `crates/nyx-engine/src/{commands,guide,analyse,pml}.rs`, `crates/nyx-engine/cards/*.md`, `crates/app/src/{app,state,tabs,forms,examples,shell}.rs`, `crates/nyx-engine/tests/layering.rs`.

## Contract

`mcp_tools.v2` (the tool list), `engine_commands.v1.json` (1.1: the new commands).
