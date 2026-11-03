# Restructure: status

Handoff file between working sessions (see `CLAUDE.md`). The plan is `docs/architecture/RESTRUCTURE_PLAN.md`; read it first. Update this file before finishing every session.

Last updated: 2026-11-03, by the session that wrote the plan. **No phase has started.**

## Phases

| Phase | Name | State | Bookmark | Notes |
|---|---|---|---|---|
| 0 | Baselines | not started | | tags `pre-restructure`, `sbs-last`; record timings and test counts below |
| 1 | Delete sbs, chips | not started | | gate G1 first |
| 2 | femcore: one mesh model, all meshing | not started | | |
| 3 | femcore: one linear algebra | not started | | gate G2 |
| 4 | nyx-physics and physics crates | not started | | gate G3 |
| 5 | nyx-analysis | not started | | |
| 6 | Materials | not started | | gate G4 |
| 7 | nyx-engine, MCP parity, thin app | not started | | gate G5 |
| 8 | File format | not started | | gate G6 |
| 9 | Close-out | not started | | |

## Baselines (fill in during phase 0)

* Test counts (`cargo test --release --workspace --no-fail-fast`): _pending_. Expected on main a047506: 632 pass, 10 ignored, 7 fail (golden_frames).
* Benchmark timings (`bench`, `profile3d`, `mem_profile`, `cavity_bench`, `solve_scaling`): _pending_.
* Contract suites: pass (`cargo test -p femcore -p geom -p nyx-store --test contracts`, `cd ai && cargo test -p nyx-mcp --test contracts`).

## Decisions taken (paul, 2026-11-03)

sbs removed; a few physics per project; separate crates per physics; no built-in devices; femcore = whole FEM core; analysis crate; TOML materials with Voigt tensors and dispersion models; MCP through the store with the same capabilities as the window; XML retired with **no migration**; compile-time physics.

## Open questions for paul

See section 8 (gates) and section 7 (parked recommendations) of the plan. Current open items: the numerical tolerance / reference platform (R22), keeping `sms_kz` and the axisymmetric variants, Structurizr validation of `workspace.dsl`.

## Log

* 2026-11-03: architecture model, 23 ADRs, 28 rough spots, contract tests (bookmark `arch-model-adrs`), restructure plan written. Nothing in `crates/` or `ai/` source changed.
