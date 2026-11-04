# nyx (repository `nyx2`)

Desktop finite-element application for photonics and phononics, in pure Rust. No Python is needed to build or run anything. This repository is the restructure of `nyx`: read `docs/architecture/RESTRUCTURE_PLAN.md` and `docs/RESTRUCTURE_STATUS.md` **before changing crate layout, mesh types, solvers or the problem types**. `README.md` is the feature tour (parts of it describe the pre-restructure state) and `docs/HOW_IT_WORKS.md` explains what the solver does and what each log line means.

## Layout

| Path | What |
|---|---|
| `crates/geom` | geometry kernel (parametric wafer-frame models, selectors, golden files) |
| `crates/femcore` | the FEM core: meshes, meshing, sparse and eigen solvers, the old physics kernels (`dfs`, `sms*`, `ems*`), analysis, use cases. Shrinks as the plan proceeds |
| `crates/nyx-physics` | the `Physics` trait, the problem envelope (`ProblemSpec`), descriptors, `SolveResult`, `Registry`, validation, JSON Schemas (ADR 25) |
| `crates/nyx-dfs`, `nyx-sms`, `nyx-ems` | one physics each. **Skeletons**: descriptors and contract tests now, `solve` answers `NotYetMoved` until phase 4 moves the kernels in |
| `crates/nyx-analysis` | analyses of results (skeleton; polarization moves in at phase 5) |
| `crates/nyx-store` | workspace, workflow, experiment and run storage |
| `crates/nyx-engine` | the registry of compiled-in physics (cargo features `dfs sms ems`) and the layering test; `Project` and commands arrive in phase 7 |
| `crates/app` | the egui application (package and binary `nyx`) |
| `crates/nyx-session`, `nyx-mcp`, `nyx-connect` | the AI interface (stdio MCP server, one-click host connection); formerly the separate `ai/` workspace |
| `docs/` | design notes and per-project status files |
| `docs/architecture/` | C4 model (`workspace.dsl`), ADRs (`adr/`, index `ADR_INDEX.md`), `ROUGH_SPOTS.md`, `CODE_TO_ADR.md`, interface `contracts/`, the restructure plan. Read `README.md` there first |
| `tools/` | `baseline.sh` (tests and benchmark timings recorded per phase) |
| `runs/`, `dist/` | local outputs, git-ignored; never commit results |

Dependency rules are in `crates/nyx-engine/tests/layering.rs` (a table plus the list of temporary edges and the phase that removes each). Nothing outside the physics crates names a physics.

## Build and test

```
cargo run --release -p nyx                                  # the application
cargo test --release --workspace                            # everything; a few minutes
cargo test --release --workspace -- --ignored               # slow tests (Bragg ~1.5 min, solver library)
cargo clippy --release -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine --all-targets --no-deps -- -D warnings   # the new crates are clean
cargo fmt -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine   # the new crates are rustfmt-formatted (settings as geom)
cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering   # interface contracts (ADR 19)
bash docs/architecture/check_code_to_adr.sh                 # every file has an ADR row
bash tools/baseline.sh <dir>                                # test counts and the five benchmark timings (BENCH_ONLY=1 for timings only)
```

Always use `--release`: the solvers are unusably slow unoptimised. The toolchain is pinned in `rust-toolchain.toml`. All jj workspaces of this repository share the build directory `../nyx2-target` (`.cargo/config.toml`); cargo serialises builds in it, so the first build in a new workspace waits for the lock but does not recompile egui, wgpu or faer. If a long build or test is needed, run it in the background.

## Conventions

* Units are SI (metres) everywhere in the library and the MCP tools.
* **Wafer frame**: x lateral, z normal to the wafer, waves propagate along **+y** (`exp(-j beta y)`). `femcore::section` is the one place that knows axis numbers (`plane_tensor`, `local_tensor`, `LOCAL_AXES`, `wafer_point`); do not scatter axis permutations elsewhere. Solid kernels use the Voigt order (x, z, y).
* Physics carry three-letter codes chosen by Paul: `dfs` (diffusion), `sms` (solid mechanics), `ems` (EM waves). Brillouin (`sbs`) is removed (ADR 24). Name new physics the same way; a physics is a crate behind the `Physics` trait (ADR 25).
* Golden files (`crates/femcore/tests/golden_frames/*.json`, `crates/geom/golden/`) are the numerical oracle. `NYX_WRITE_GOLDEN=1` rewrites them; do it only for a kernel known to be right, and say so in the commit. Never regenerate a golden to hide a difference (ADR 10).
* **Tests never compare JSON text or key order** (ADR 26): `serde_json` keeps insertion order when any crate of the build enables `preserve_order` and sorts otherwise, so `--workspace` and `-p` builds differ. Compare `Value`s or `nyx_physics::canonical_json`.
* Contract exemplars (`docs/architecture/contracts/`) change only with a line in `VERSIONS.md` (ADR 19). `NYX_WRITE_CONTRACT=1` rewrites them.
* Every source file of `femcore` and of the new crates starts with a short comment saying what it does; keep that when adding files.
* `*.py` is git-ignored on purpose: nothing Python ships. Any GDS layout is built with gdsfactory, never raw gdstk.
* Determinism: the same input gives the same numbers for any thread count. `NYX_THREADS=N` sets the rayon pool, `NYX_PROFILE=1` echoes the solver log.

## Known state of the checks (foundation, 2026-11-04, Windows, Rust 1.99.1)

* `cargo test --release --workspace --no-fail-fast`: 764 pass, 9 ignored, 2 fail: the 2 tests in `femcore --test golden_frames`: `golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`. The differences are small (about 1e-6 relative in fields, up to 2e-5 in one effective index), identical with `NYX_THREADS=1`, and present at the commit that wrote the golden (35cceba), so they come from the machine the goldens were frozen on, not from a later change (R22, gate G2). Do not regenerate them to hide it.
* `cargo clippy --workspace -D warnings` is **not** clean with 1.99.1 for the old crates (8 errors in `femcore`: `ems_pml.rs`, `sms.rs`, `sms_kz.rs`); the new crates are clean. `cargo fmt --check` fails widely for the old crates (R23, decided in phase 9).
* The numbers and timings of the pre-restructure state and of each phase are in `docs/RESTRUCTURE_STATUS.md`.

## Version control

* `nyx2` is a jj repository (colocated git), cloned from `nyx` with its full history, **with no remote: never push**. Do not touch the old `nyx` repository. `main` is the starting point; never commit to it.
* One jj workspace per working session (`jj workspace add ../nyx2-<name> -r <base>`), one bookmark per step, named `restructure-<n>-<name>` (ADR 27). Commit messages say what changed and why in one line.
* `docs/RESTRUCTURE_STATUS.md` is the handoff between sessions: read it before starting, update your own section before finishing. Where the plan conflicts with what you find in the code, stop, write the conflict into the status file and ask paul.
* Do not leave work only in an uncommitted working tree or a side folder; commit it to your bookmark.
