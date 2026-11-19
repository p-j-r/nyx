# 84. Results are stale per physics: the key is the problem key

Date: 2026-11-16

## Status

Accepted by paul (productivity plan question 1, 2026-11-16). **Number 84: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.** Refines ADR 64 (cache keys); implements step P1 of `docs/NYX2_PRODUCTIVITY_PLAN.md`.

## Context

A result was stamped with a hash of the whole project, so any edit made every result stale, and the window dropped the flag that the engine already returned: a plot of a solve from before an edit looked current. Two solves of one physics could finish out of order (an AI host can start two) and the last to finish became the stored result whatever it was made from.

## Decision

* **The key of a result is the problem key of its physics** (`Engine::problem_key`): a hash of that physics' problem, the materials the project carries, the mesh key (geometry, mesh settings, sizing) and the problem keys of the physics it takes fields from (its imports, found through the `exports` and `imports` of the `Physics` trait, as `solve_order` does). The active physics, the title and everything about how a result is looked at (colormap, camera, field, range) are not in it. A result is *stale* when its key is not the key its physics has now. Staleness follows the content: undoing an edit makes the result current again.
* **It is a cache key in the sense of ADR 64**: computed in the process, never stored in a file, never sent. The replies of the commands carry `stale` and `age_s`, not the key.
* **Late jobs.** A finished solve replaces the stored result of a physics only when its key is the one the physics has now, or when the stored result is stale and the job started later (jobs are numbered in the order they start). Otherwise it is recorded in the run table but not stored; the reply says `current` and `stored` per result. The same rule keeps an old mesh from displacing the mesh of the project. A run is recorded with the mesh of its own job.
* **The window** shows a stale result as stale: a mark on the plot, a band with "Solve again" in the Results tab and panel, "(stale)" in the list of results. A result is stale for the window when the engine says so, when a newer result of the physics replaced it, or when the Geometry or Mesh tab holds edits that were not applied to the engine yet.

## Known limits

* **Not an identity ("same key, same numbers").** `run_chain` hands a physics the fields of the physics that were *in the same job*, so `solve_run` of `sms` alone and of `dfs, sms` have the same `sms` key and different numbers. For staleness this is acceptable: both are results of the current problem. A later use of the key to deduplicate work (the queue, plan step Q1) must add what was chained.
* **The whole `materials` map of the project is hashed**, not only the entries a problem names, so editing a project material that a physics does not use stales it. This errs on the safe side; the library files are loaded once per session and are not in the key.
* **A mesh sized from the physics** (`physics_sizing`) puts the wavelengths of every physics in the mesh key, so a wavelength edit then stales every result on that mesh. That is correct (the mesh is different); with `physics_sizing` off it stales only its own physics.
* The key is computed on each call (no cache): about 0.2 ms for the examples and 5 ms at a project of 1,000 features (measured, `docs/PRODUCTIVITY_STATUS.md`). Caching waits for a model that needs it.

## Consequences

* An edit of one physics no longer marks the others stale; an edit of the geometry still marks all.
* An MCP host sees `stale` on every answer that quotes a result; it cannot quote a stale number without the JSON saying so.

## Code

`crates/nyx-engine/src/engine.rs` (`problem_key`, `problem_keys`, `is_stale`, `apply`), `commands.rs`, `analyse.rs`, `results_file.rs`, `crates/app/src/{app,tabs}.rs`; tests `crates/nyx-engine/tests/stale.rs` and the window tests `a_result_on_show_says_when_the_model_changed`, `a_result_the_engine_did_not_keep_is_not_put_on_show`.

## Contract

`engine_commands` 1.9: descriptions of `result_info`, `result_field`, `analyse`, `list_analyses`; replies gain `stale`, `age_s`, `current`, `stored`.
