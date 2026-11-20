# Latency work (L1 to L3): status

Handoff for the work that came out of the benchmark session of 2026-11-19 (two silicon-on-sapphire waveguides driven from Antigravity; the analysis is `LATENCY_AND_WORKFLOW_ANALYSIS.md` in the benchmark workspace). Decision record: ADR 97. For people and agents using a host: `docs/HOST_SETUP.md`.

Workspace `nyx2-latency` (`C:\Users\robin-adm\Documents\opomech\nyx2-latency`), based on `main` = `productivity-merge-main` 219fe4e1. Three bookmarks, stacked:

| Bookmark | What |
|---|---|
| `latency-1-mcp-reader` | L1. `nyx-mcp` reads stdin on a thread; a wait is cut to `--max-wait` (default 90 s, 1 to 900) and ends for a newer `script_run` or `script_wait` or `notifications/cancelled`; `ping` and `tools/list` are answered during a wait; progress notifications; the job (stage, seconds, estimate) in a `running` answer. `crates/nyx-mcp`, `Kernel::waiting_job`. |
| `latency-2-estimates-noise` | L2, L3. Noise rows in `convergence_check`; `dry_run`, `max_cells`, `max_seconds`; `cost` in the answers of `mesh_build` and `sweep_run`; starred materials in the wave sizing of the mesher; the guide cards; `docs/HOST_SETUP.md`; ADR 97; `engine_commands` 1.22. |
| `latency-3-connect-antigravity` | L4. `nyx-connect`: host `antigravity`, `agent_instructions` (`nyx-connect instructions`, button "Copy instructions for an AI"), the double click, the hint in the Connect menu; ADR 98; `connect_hosts` and `connect_report` 1.1. |

## What was found, in the code

* The wait limit was 900 s and the server handled one line of stdin at a time: confirmed (`crates/nyx-mcp/src/main.rs` before this work).
* `convergence_check` divided by the larger of the two values with a floor of 1e-12 on the value itself: an imaginary part of 1e-10 is a 100 % difference. Reproduced with sapphire at `euler_deg [0, 57.6, 0]` (the SOS wire): `neff 2.685547 + 5.5e-11 i`, `Q 1 = 2.4e10`, and on the finer mesh Q 1.3e10. Test `a_lossless_guide_on_a_rotated_crystal_converges_although_its_imaginary_part_is_rounding`.
* **The starred material was a real bug on `main`:** `set_material region="BOX*"` with `physics_sizing` true made `mesh_build` fail ("the sizing names region 'BOX*', which the model does not have"); the solver expanded the star, the mesher's wave sizing did not. `Sizing::expand_wildcards`.
* `confinement` with `weight "eps_r"` on a tensor crystal is not silently wrong: it is refused with the names of the entries (`eps_r_xx` ...). The card gives the contraction; a test runs it.
* Not checked: the claim that a loop variable reused over different literal lists fails the type checker (`ty`). The card says to use a new name; nothing was reproduced.

## Decisions taken without asking (change them if you disagree)

* The noise floor is `im < 1e-9 * re` on both meshes (a Q above 5e8). A cavity with Q 1e6 is never noise. The advisor's alternative (noise only when the problem has no loss mechanism) needs each physics to say so; not done.
* `wait_s` above the limit is cut, not refused, and the answer says so.
* Only `script_run` and `script_wait` end a wait; `api` and `script_check` wait their turn.
* The history of solves does not know the complex solver from the real one, so an estimate across both is wide; the answer carries the number of runs and `extrapolated`.

## Not done

* Antigravity in Connect: done as L4 (bookmark `latency-3-connect-antigravity`, ADR 98): a host row, `nyx-connect instructions` and the button "Copy instructions for an AI", a double click on `nyx-connect.exe` shows the commands. Not tried on a running Antigravity (whether it reloads the file without a restart).
* The A-benchmarks, `tools/ux-shots.sh`, `tools/baseline.sh` and the `-- --ignored` tests were not run on this build.
* The C4 pages in `docs/architecture/c4` are rebuilt; the copy in `Documents/opomech/nyx-c4` is not.
* No run with a real host: the stdio tests drive the binary the way a host does (a newer request, a cancellation, a ping, a progress token), Antigravity itself was not tried.

## Checks (2026-11-20, Windows)

* `cargo test --release --workspace --no-fail-fast`: 1232 pass, 0 fail, 12 ignored (main was 1220 / 0 / 12). The last two edits (the sweep cost test, a clippy fix) came after that run; `latency`, `sweep`, `api_docs` and `contracts` were run again and pass.
* `cargo clippy --release --workspace --all-targets --no-deps -- -D warnings` clean; `cargo fmt -p nyx-engine --check` clean; `check_code_to_adr.sh` 0 without a row; C4 rebuilt.
* New tests: `crates/nyx-mcp/tests/stdio.rs` (the cut wait, a newer request, a cancellation, a ping and a progress token while a solve runs), `crates/nyx-engine/tests/latency.rs` (starred material, dry run and limits, noise rows, estimates), unit tests in `convergence.rs`, `preview.rs` and `femcore::meshing::sizing`.
