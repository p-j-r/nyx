# 92. Feedback while editing is a value, never a job

Date: 2026-11-18

## Status

Accepted as step P3 of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3). **Number 92: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.**

## Context

The plan applies the paper's point (a response inside a second keeps the train of thought) **per interaction**, and says that nothing ever starts a mesh-and-solve chain on its own. P3 had to decide, interaction by interaction, from numbers measured on the merged build (UX stack and P6 included; `docs/PRODUCTIVITY_STATUS.md`, section "P3 done"), which ones need work and which are left alone. The measurements, on this machine (i5-13500), with a frame log that now says where each frame's time went (`NYX_FRAME_LOG` with `spans`):

* **Viewing** (rotate, switch solution or field, 3D cavity of 9,480 cells, 56,226 unknowns, a second solve running): frame p95 18.4 ms, p99 19.3 ms, max 25.8 ms. The data pass of 362 field and solution switches took 2.5 ms at most. The gate (33 ms) holds and nothing is precomputed.
* **Geometry preview, input check, mesh, on every example** (a fresh engine per repetition): geometry 0.05 to 0.51 ms, check 0.01 to 0.09 ms, first mesh build 0.5 to 12.1 ms (2D up to 8,012 cells; 3D swept up to 6,936 cells). Inside the 100 ms of the plan everywhere.
* **Run table at 500 runs:** UI 4.4 ms p50, 6.0 ms p95, frame p95 18 ms. Inside the budget, left alone.
* **Polarization** costs 1.4 to 5 microseconds per cell on the UI thread (9 ms at 9,480 cells in 3D, 15 ms at 10,800 cells of the SOI wire, 41 to 48 ms at the 8,012 cells of the silicon-on-sapphire example): over a frame at moderate sizes.
* **The frame in which a solve's result is taken** costs 53 to 63 ms (once per solve): the engine applying the result 9 to 19 ms, the polarization 9 ms, recording the run (files, digest, manifest) 9 to 18 ms, refreshing the dashboard entries 4 to 6 ms; the frame after it lays out the Results tab for the first time, 29 ms (13.8 ms of it the solver log, 20.7 ms the tab).

## Decision

* **Feedback is a value the front-end asks for, never a job.** Three engine functions answer a question about an edit that is not applied yet, or about a result, without changing the engine, its mesh, its results or its queue (`nyx_engine::preview`):
  * `check_edits(geometry, mesh)`: `Project::check` on a copy with the edits in place (the feature list, the element size, every problem against its physics). Microseconds; runs on every change of the Geometry form.
  * `size_preview(geometry, mesh)`: a `SizeTask`, `Send`, that meshes a copy of the edits (the mesher and sizing of the project, no stage cache) and looks up the memory and time estimate of ADR 86 for that many cells. The caller runs it on a thread of its own. Offered only where the mesh is cheap: 2D, and 3D swept from a 2D section; any other 3D mesh is refused with the reason ("can take seconds: it is built by Build mesh"). The answer says the cells, the dimension, the seconds the trial took and the estimate with its basis, or no estimate when fewer than two runs of different sizes exist (nothing is invented).
  * `polarization_task(physics)`: the analysis of ADR 56 as a function that owns its inputs, with the result it was made for, so a front-end can compute it beside the frame and tell whether the answer still belongs to what is on show.
* **The window** (`feedback.rs`) shows the problem of the edits as a glyph, the word "Problem" and the engine's words, and, if there is none, "About N cells, 2D (a trial mesh of these settings, not built)" with the cost line, under the mesh controls of the Geometry form. The trial mesh starts when an edit has rested for 200 ms and is dropped if the edit moved on. **A trial mesh that takes 100 ms or more turns the automatic preview off** (the rule of plan question 2), says so, and waits for the person ("Estimate size"). The polarization of a new result is computed on a thread and lands in the result view that holds that very result.
* **The mesh preview is numbers, not a drawing.** The plan allowed a small 2D mesh preview if P0 showed the mesher cheap; it is cheap (at most 12 ms), and what the person needs while editing the element size is the count and the cost. Drawing a trial mesh in the viewport needs a viewport content that is not the engine's model, which the UX stack owns; it is not built, and a later decision can add it on top of `size_preview` (the task already meshes the edits).
* **No mesh, solve, run or job without a request** is a test, not a promise: `editing_every_form_starts_nothing` edits every buffered form of the window (element size, sizing from the physics, adaptive, absorbing layers, sweep of the mesh, a refused size, and the Sweep form of the Build page), draws frames until every check and trial mesh has run (and asserts the last one came back), and asserts that no job was submitted, no run recorded, the engine's mesh and project are the same objects, the result on show is still there and nothing is busy. Compute, a queue submission, Run sweep and an MCP command remain the only things that start work.
* **Left alone, with the numbers above:** viewing, the geometry preview of the feature list, the run table, solving (a request), convergence checks, exports and sweeps (background requests). Not changed in this step although they are over a frame once per solve: recording the run, applying the result in the engine and the first layout of the Results tab. Each is named, with its cost, in the status file as follow-up work.
* **Diagnostics:** `NYX_FRAME_LOG` rows gain a `spans` column (`spans.rs`: named timings of `poll`, `take_solved`, `record_ok`, the Results tab and its parts). Nothing is read or stored when the variable is unset. The example `p3measure` times the check, the cold mesh and the polarization on every example project.

## Consequences

* An edit that cannot be built is named at once, not at Compute; the size of an edit is on screen before a mesh is built; the polarization no longer lengthens the frame in which a result arrives (by 9 ms at 9,480 cells, more on finer 2D meshes).
* A trial mesh is a second mesh of the same edits when Build mesh is pressed. That costs milliseconds for the models it is offered for and nothing for the others; the stage cache of the engine is deliberately not used for it, so a trial mesh can never become the mesh a solve reads.
* The estimate is only as good as the runs this machine has done; the line says "no cost estimate yet" until two sizes exist.
* An analysis computed on a thread can arrive after the person has moved on: it is applied only to the result view that holds the very result it was computed from (`Arc::ptr_eq`), otherwise dropped.

## Code

`crates/nyx-engine/src/preview.rs`, `analyse.rs` (`analysis_inputs` is crate-visible), `engine.rs` (`store_results`), `estimate.rs` (the history writer), `crates/nyx-store/src/flow.rs` (one manifest write per run), `crates/app/src/{feedback,spans,app,tabs,script,shell,wfview}.rs`, `crates/nyx-engine/examples/p3measure.rs`, `tools/ux-shots.sh` (two more shots). Tests: `crates/nyx-engine/tests/preview.rs` and the window tests `editing_every_form_starts_nothing`, `the_size_of_an_edit_arrives_from_a_trial_mesh` and `the_polarization_arrives_after_the_result_and_only_for_the_result_it_was_made_from`.

## Contract

None changed: no command was added and no answer changed, so `engine_commands`, `mcp_tools` and the exemplars keep their versions.

## Addendum: the frame that takes a result (same step, second commit)

The first commit left one stutter per solve: the frame that takes a result was 53 to 63 ms and the next one, which draws the Results tab for the first time, 37 to 45 ms. A frame log with spans attributed them, and each part was dealt with by what it is:

* **The engine's apply step copied every result** (`Arc::new(r.clone())`, 6 ms for 56,226 unknowns). Results are now moved into the store after the run is recorded (`Engine::store_results`); what is stored, and under which key, is unchanged. The failed-chain path stores what was solved as before.
* **The history appended to its file on the calling thread** (10 ms for an open and a write on this machine). `History::add` updates the rows in memory at once, which is what estimates read, and a thread of the history writes the line; lines are written in order, `flush` waits for them, and dropping the history flushes. The format of the file is unchanged.
* **The run record is 9 to 17 ms of file operations** (freeze the spec, read the run table twice, append, write the digest and the manifest). The store now writes the manifest once instead of twice per run, and the window writes the record **on the frame after the result is taken**, and the dashboard refresh (4 to 6 ms) on the one after that (`shell::Deferred`, one piece per `poll`, never in the frame that took a result). A solve's run context is taken with its result, so a second solve cannot be recorded under the first one's. Scripts and the test helper wait for the queue to be empty before going on.
* **The solver log starts closed.** Drawing its eleven short monospace lines for the first time cost 13 to 14 ms (one label for the whole log, and laying the glyphs out beforehand, both measured the same or worse and were dropped); the header still says how many lines it holds. This is a visible default changed on purpose, and `default_open` is the one place to change it back.

Measured again with the same 3D case, every frame around a result is now at most 31 ms (the worst: the first Results draw 30.6 ms, the run record 28.8 ms), against 53 to 63 and 37 to 45 ms.
