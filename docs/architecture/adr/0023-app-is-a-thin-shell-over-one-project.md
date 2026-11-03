# 23. The app is a shell over one saved Project: jobs on worker threads, the same code headless

Date: 2026-10-28

## Status

Accepted by practice (`6de7458`, `2b5b66a`, `37e5871`). The plan to move `Project` out of the app (MCP plan, phase 0) was not carried out (R16).

## Context

A researcher edits a model, meshes, solves, looks at results and iterates; the same model must also run on a workstation or cluster without a window, and be driven by tests and screenshot scripts.

## Decision

* `state::Project` is what is saved: `GeometrySpec` (a geom `FeatureList`, global element size `h`, optional extrusion layers, `adapt_scale`, `physics_sizing`), one `ProblemSpec` per physics, the active physics and the colour map. It is written as XML through `xmlio` or JSON. `ResultView` wraps a `SolveResult` with plot state.
* `App` holds the project, the mesh in use (`Doc`), the results and the viewport. Anything slow (building a mesh, opening a file, solving) is a numbered **job on its own thread**; only one runs at a time, and a message for a job that is no longer current cannot end the busy state. A solve has a `Progress` handle: the toolbar shows its stage and a Cancel button (Esc too), and the stage log is part of the result.
* A solve writes the results file requested in the study and an **autosave** of the run (project XML, stage log, optionally HDF5) to the runs folder as `temp_untitled_YYYYMMDD_HHMMSS`.
* **Headless**: `nyx --batch PROJECT.nyx.xml [--physics ...] [--output serial|parallel FILE.h5]` runs the same code as the window; in a parallel (MPI) run every process solves the whole problem and writes its own block of rows. `nyx --script "cmd; cmd"` drives the window for tests and screenshots (`example N`, `compute`, `shot`, ...). `App::with_context` builds the state without a window for tests.
* The viewport splits drawing by what changes (data, rotation, camera) and caches each part (`scene.rs`), on the shared thread pool.
* The workflow dashboard and workspace live beside this in `shell`, `wfview`, `sheets`, backed by `nyx-store` (ADR 14); the Connect menu registers the MCP server (ADR 15).

## Consequences

* The session layer (`ai/`) cannot load or save a `Project` because the type is in a binary crate (R16).
* The number of the project file's `version` (1) is not read anywhere (R20).
* `app.rs` (about 1.1k lines) and `tabs.rs` hold most behaviour; the workflow design anticipated a refactor of about 2,800 lines (`WORKFLOWS_DESIGN.md` risk 7).
* One job at a time means a long Brillouin solve blocks a mesh rebuild; detached workers are planned for AI sessions (ADR 15).

## Code

`crates/app/src/{main,app,state,worker,batch,script,examples,tabs,geomui,viewport,scene,colormap,theme,kit,dashboard,shell,wfview,sheets,connect}.rs`.

## Contract

`contracts/project` is not written: the exemplar for the saved project would have to live in the app crate. The `ProblemSpec` part is `contracts/problem_spec`.
