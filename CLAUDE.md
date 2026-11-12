# nyx (repository `nyx2`)

Desktop finite-element application for photonics and phononics, in pure Rust. No Python is needed to build or run anything. This repository is the restructure of `nyx`: read `docs/architecture/RESTRUCTURE_PLAN.md` and `docs/RESTRUCTURE_STATUS.md` **before changing crate layout, mesh types, solvers or the problem types**. `README.md` is the feature tour (parts of it describe the pre-restructure state) and `docs/HOW_IT_WORKS.md` explains what the solver does and what each log line means.

## Layout

| Path | What |
|---|---|
| `crates/geom` | geometry kernel (parametric wafer-frame models, selectors, golden files) |
| `crates/femcore` | the FEM core: meshes (`mesh`, with `Mesh`, the exchange mesh the physics take), meshing, adaptation, the linear algebra (`linalg`), the finite-element basics (`fem`, `edge`, `pml`, `periodic`, `props` for materials per domain), the material specification and library (`material`, TOML files in `materials/`), `props` (materials per domain), and fixtures for tests (`testkit`, `testkit::soi`). The physics kernels, the old window problem types and the use cases are gone from it |
| `crates/nyx-physics` | the `Physics` trait, the problem envelope (`ProblemSpec`), descriptors, `SolveResult`, `Registry`, validation, JSON Schemas (ADR 25) |
| `crates/nyx-dfs`, `nyx-sms`, `nyx-ems` | one physics each since phase 4 (ADR 41 to 47): its typed `Problem` and `decode` of the envelope, the kernels in three stages, `run`, `Physics::solve`, descriptors, its goldens and tests. |
| `crates/nyx-analysis` | analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a `SolveResult` and its mesh (phase 5, ADR 56). A mode study states its frame in the scalar `field frame` (`nyx_physics::FieldFrame`) and carries the complex field in its cell fields; results carry no analysis |
| `crates/nyx-store` | workspace, workflow, experiment and run storage (JSON specs; no dependency on another crate of the workspace, ADR 59) |
| `crates/nyx-engine` | `Engine`: the project, commands (every action by name, with JSON Schemas; the MCP tools are exactly the catalog, the window uses the same commands), the job runner, the geometry operation descriptors, the materials of a session (library plus the project's own, ADR 58), analyses, the guide cards, store wiring, and the registry of compiled-in physics (cargo features `dfs sms ems`). `facade` re-exports what a front-end draws with. Examples (`examples/`) are the benchmark programs `bench`, `profile3d`, `mem_profile`, `cavity_bench`, `stage_log`, `eig_beam` |
| `crates/app` | the egui application (package and binary `nyx`): a front-end of `nyx-engine`; the Study, Physics and Materials tabs are drawn from the physics descriptors (`forms.rs`) |
| `crates/nyx-mcp`, `nyx-connect` | the AI interface (stdio MCP server whose tools are the engine commands, ADR 57; one-click host connection) |
| `examples/`, `materials/` | example projects (`*.nyx.json`, written by `crates/nyx-engine/tests/make_examples.rs`) and the material library (`*.toml`), both found next to the program or in a parent folder |
| `docs/` | design notes and per-project status files |
| `docs/architecture/` | C4 model (`workspace.dsl`), ADRs (`adr/`, index `ADR_INDEX.md`), `ROUGH_SPOTS.md`, `CODE_TO_ADR.md`, interface `contracts/`, the restructure plan. Read `README.md` there first |
| `tools/` | `baseline.sh` (tests and benchmark timings recorded per phase) |
| `runs/`, `dist/` | local outputs, git-ignored; never commit results |

Dependency rules are in `crates/nyx-engine/tests/layering.rs` (a table; the list of temporary edges is empty since phase 7). Nothing outside the physics crates names a physics. The window and the MCP server depend on `nyx-engine` only (`nyx-connect` for the Connect menu).

## Build and test

```
cargo run --release -p nyx                                  # the application
cargo test --release --workspace                            # everything; a few minutes
cargo test --release --workspace -- --ignored               # slow tests (solver library); `-p nyx-engine --test make_examples -- --ignored` rewrites examples/*.nyx.json
cargo clippy --release -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine --all-targets --no-deps -- -D warnings   # the new crates are clean
cargo fmt -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine   # the new crates are rustfmt-formatted (settings as geom)
cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering   # interface contracts (ADR 19)
bash docs/architecture/check_code_to_adr.sh                 # every file has an ADR row
HDF5_DIR=<prefix of the HDF5 library> cargo test --release -p nyx-engine --features hdf5 --test h5   # results files (needs the library, its bin folder on PATH)
bash tools/baseline.sh <dir>                                # test counts and the five benchmark timings (BENCH_ONLY=1 for timings only)
```

Always use `--release`: the solvers are unusably slow unoptimised. The toolchain is pinned in `rust-toolchain.toml`. Each jj workspace should build in its own target directory (`CARGO_TARGET_DIR`), several sessions share a machine. If a long build or test is needed, run it in the background.

## Conventions

* Units are SI (metres) everywhere in the library and the MCP tools.
* **Wafer frame**: x lateral, z normal to the wafer, waves propagate along **+y** (`exp(-j beta y)`). `femcore::section` is the one place that knows axis numbers (`plane_tensor`, `local_tensor`, `LOCAL_AXES`, `wafer_point`); do not scatter axis permutations elsewhere. Solid kernels use the Voigt order (x, z, y).
* Physics carry three-letter codes chosen by Paul: `dfs` (diffusion), `sms` (solid mechanics), `ems` (EM waves). Brillouin (`sbs`) is removed (ADR 24). Name new physics the same way; a physics is a crate behind the `Physics` trait (ADR 25).
* Golden files (`crates/nyx-{dfs,sms,ems}/tests/golden/*.json`, `crates/geom/golden/`) are the numerical oracle; the harness is `nyx_physics::testkit`. `NYX_WRITE_GOLDEN=1` rewrites them; do it only for a kernel known to be right, and say so in the commit. Never regenerate a golden to hide a difference (ADR 10).
* **Tests never compare JSON text or key order** (ADR 26): `serde_json` keeps insertion order when any crate of the build enables `preserve_order` and sorts otherwise, so `--workspace` and `-p` builds differ. Compare `Value`s or `nyx_physics::canonical_json`.
* Contract exemplars (`docs/architecture/contracts/`) change only with a line in `VERSIONS.md` (ADR 19). `NYX_WRITE_CONTRACT=1` rewrites them.
* Every source file of `femcore` and of the new crates starts with a short comment saying what it does; keep that when adding files.
* `*.py` is git-ignored on purpose: nothing Python ships. Any GDS layout is built with gdsfactory, never raw gdstk.
* Determinism: the same input gives the same numbers for any thread count. `NYX_THREADS=N` sets the rayon pool, `NYX_PROFILE=1` echoes the solver log.
* **File format** (ADR 63): `Project`, `FeatureList` and `ProblemSpec` carry `format` (value 1); files without it or of another format are refused, nothing is converted. Project and spec files are JSON (`*.nyx.json`, `experiment.nyx.json`); results and meshes are HDF5 (feature `hdf5`), never JSON.

## Known state of the checks (phase 7, 2026-11-05, Windows, Rust 1.99.1)

* `cargo test --release --workspace --no-fail-fast`: 894 pass, 10 ignored, 0 fail (phase 8, `restructure-8-format`; 889 / 10 / 0 at phase 7, `restructure-7-engine-frontends`; phase 5 closed at 815 / 11 / 0, phase 4 at 790 / 9 / 0). The two SOI-wire goldens that failed before phase 3 pass at the G2 tolerances (5e-5 relative on n_eff, 1e-3 on field maxima; `nyx-ems/tests/golden`); a geom area property test can be flaky. Do not regenerate a golden to hide a difference.
* `cargo clippy -D warnings` is clean for `nyx-physics`, `nyx-dfs`, `nyx-sms`, `nyx-ems`, `nyx-analysis` and `nyx-engine`; the old `femcore` code has a few warnings left. `cargo fmt --check` fails for the old code (R23, decided in phase 9).
* The numbers and timings of the pre-restructure state and of each phase are in `docs/RESTRUCTURE_STATUS.md`.

## Version control

* `nyx2` is a jj repository (colocated git), cloned from `nyx` with its full history, **with no remote: never push**. Do not touch the old `nyx` repository. `main` is the starting point; never commit to it.
* One jj workspace per working session (`jj workspace add ../nyx2-<name> -r <base>`), one bookmark per step, named `restructure-<n>-<name>` (ADR 27). Commit messages say what changed and why in one line.
* `docs/RESTRUCTURE_STATUS.md` is the handoff between sessions: read it before starting, update your own section before finishing. Where the plan conflicts with what you find in the code, stop, write the conflict into the status file and ask paul.
* Do not leave work only in an uncommitted working tree or a side folder; commit it to your bookmark.
