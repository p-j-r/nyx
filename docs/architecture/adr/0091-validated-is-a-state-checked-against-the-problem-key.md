# 91. Validated is a state: a result is checked only by a convergence check on its own problem key

Date: 2026-11-18

## Status

Accepted as step P6 of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3). **The number 91 is provisional**: it is the next free one on the productivity line (84 to 87 and 89 and 90 are its own, 88 is the UX decision), and is renumbered when the lines merge.

## Context

`convergence_check` (ADR 61) has been an engine command, an MCP tool and a Study-tab button for a long time, but its answer lived and died with the call. A result that had been checked and one that had not looked the same in `result_info`, in the window and in the run table, so a number could be quoted, or a run kept, without anything saying whether anyone had looked at the mesh. The plan asks for validation to be a state a reader can see, and adds a third state for the case the plan found in its own measurements: in 3D the check solves on a mesh finer by a factor, unknowns grow as `factor^-3` and the memory of a direct factorisation grows faster, so the check may be the one thing the machine cannot do, and a 3D result must not look merely "unchecked".

## Decision

* **Three states, one reason each** (`nyx_engine::validated`). `exploratory`: no convergence check with a stated tolerance has passed on this problem (this is also the state of a result whose last check failed, with the numbers of that check as the reason). `checked`: a check with a tolerance passed. `check_not_affordable`: not checked, and the check is expected not to fit this machine. Every answer that carries `stale` also carries `validated: {state, reason, check?}`: `project_info`, `result_info`, `result_field`, `analyse` and the reductions of ADR 90 (a reduction over several physics reports the weakest state, with one row per physics).
* **A verdict belongs to a problem, not to a result.** The engine keeps the checks it has seen by the problem key of ADR 84 (`Engine::problem_key`, taken when the check starts). A result is checked when a passed check exists for its own key. So a result solved again from the same problem, a file opened, or a project edited back to an earlier state is checked again, and an edited project is not: **a stale result is never checked**, whatever was checked before. The last check of a problem decides (a later failing check turns it back to `exploratory`).
* **A tolerance is what makes a verdict.** A check without `tolerance` records nothing and the state does not move. With one, `passed` means every compared quantity was within it, the rule `convergence_check` always had.
* **What was compared is said.** The verdict names its observable (every eigenvalue, value and numeric scalar of the result, the number of them) and how solutions were matched: **by index**, the order the solver returns them. Matching by overlap or tracking (plan R11) is not built; a 3D check over many near-degenerate modes can therefore pass or fail on a reordering, and the verdict says "matched by index" so nobody reads more into it.
* **Not affordable is an estimate, with its basis.** With no check, the engine takes the cells of the result, the cells of a mesh finer by the planned factor (0.6: `cells / 0.6^dim`) and asks the memory history of ADR 86 for that size (`History::estimate`); when the estimate with its margin exceeds the fraction of this machine's memory that admission refuses (`REFUSE_ABOVE`, 90 %), the state is `check_not_affordable` and the reason names the estimate, the runs it was fitted to and whether it is extrapolated. With no estimate nothing is claimed (the state is `exploratory`: the engine does not know that the check is out of reach). A check that was made always beats the estimate.
* **The run record keeps the verdict.** `nyx_store::runs::RunRecord` gains the optional `convergence: {factor, tolerance, observable, passed, worst_relative}`, written as a patch line (the same mechanism as star, tags and note: nothing is rewritten, the last patch wins), so a run table and a later session know the run was checked. The field is optional and written only when there is a verdict: every earlier `runs.jsonl` and the `store_run.v1` exemplar are unchanged and no format key moves (ADR 63 concerns project, feature list and problem-spec files, not run lines). The engine writes the patch into the run that solved the problem when it records runs (`Stored.run`); the window, which writes its own run records, writes it into the run on show (`App::note_check_in_run`).
* **The window** shows the state as a glyph and a word, never by colour alone (`kit::State`: `Checked` a ring with a check mark, `Exploratory` a ring half filled, `NotAffordable` a ring with a slash; ADR 88): a chip on the plot, and a line with the reason and a "Check convergence" button in the Results tab and under the convergence form of the Study tab. Nothing is drawn for a stale result: the stale band says it.

## Consequences

* An agent or a person can see whether a number rests on a check, and a failed check is as visible as a passed one. A script can refuse to quote an `exploratory` number.
* The check is still two solves, run through the queue like any heavy job (ADR 87); this decision adds no work to a solve.
* What the check compares decides what "checked" is worth: the 04 thermal bar's scalars are its boundary temperatures, so it passes any tolerance. `observable` is in the verdict for that reason; choosing the quantity a check must be about (a mode's effective index, a stress) is a later argument of `convergence_check`, not built here.
* A physics that takes fields from another is checked the way `convergence_check` solves it, on its own; its verdict is stored under its own problem key (which includes its imports' keys), so a change in the physics it imports makes it stale as before.

## Not in this decision

* Cheaper evidence counted as partial support: the Maxwell residual estimator (ADR 53), a two-mesh ladder, a Richardson estimate. None of them leaves a stored number on a result today; when one does, it belongs in `validated.check` as a second kind of evidence, not in a fourth state.
* Tracking modes between the two meshes (plan R11), a convergence check of a reduction (ADR 90) over refinement, a column for the verdict in the runs table, and a filter term for it.

## Code

`crates/nyx-engine/src/validated.rs`, `convergence.rs` (the answer names the worst quantity, the count and the matching), `engine.rs` (`checks`, `Stored.run`, `Pending::Convergence { physics, key }`), `crates/nyx-store/src/runs.rs` (`Convergence`, the patch), `crates/app/src/{activity,kit,shell,app}.rs`. Tests: `crates/nyx-engine/tests/validated.rs`, the store test of the patch, the window test `the_window_shows_the_validated_state_and_a_check_goes_into_the_run`.

## Contract

`engine_commands` 1.13 and `mcp_tools` 2.9 (provisional, additive: `validated` in the answers above, `recorded`, `run`, `worst_quantity`, `quantities` and `matching` in the `convergence_check` answer); `store_run` stays 1.0 (an optional key that is absent unless a check was made).
