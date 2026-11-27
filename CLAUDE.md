# nyx

Desktop finite-element application for photonics and phononics, in pure Rust. No Python is needed to build or run anything (tests and CI use Python 3.12 with `h5py` and `numpy` as a second HDF5 reader, and `tools/vvuq_telemetry.py` is a CI script). Read `docs/architecture/README.md`, the ADRs (`docs/architecture/adr/`) and `crates/nyx-engine/tests/layering.rs` **before changing crate layout, mesh types, solvers or the problem types**. The restructure plan and the session handoff logs of the earlier work are kept in `docs/history/` for reference only; they are not current, and an ADR wins over them. This file is the one source of truth for agents (`AGENTS.md` points here). `README.md` is the overview of what nyx is, how to build and drive it, and what it solves, and `docs/HOW_IT_WORKS.md` explains what the solver does and what each log line means.

## Layout

| Path | What |
|---|---|
| `crates/geom` | geometry kernel (parametric wafer-frame models, selectors, golden files) |
| `crates/femcore` | the FEM core: meshes (`mesh`, with `Mesh`, the exchange mesh the physics take), meshing, adaptation, the linear algebra (`linalg`), the finite-element basics (`fem`, `edge`, `pml`, `periodic`, `props` for materials per domain), the material specification and library (`material`, TOML files in `materials/`), `props` (materials per domain), and fixtures for tests (`testkit`, `testkit::soi`). The physics kernels, the old window problem types and the use cases are gone from it |
| `crates/nyx-physics` | the `Physics` trait, the problem envelope (`ProblemSpec`), descriptors, `SolveResult`, `Registry`, validation, JSON Schemas (ADR 25) |
| `crates/nyx-dfs`, `nyx-sms`, `nyx-ems` | one physics each since phase 4 (ADR 41 to 47): its typed `Problem` and `decode` of the envelope, the kernels in three stages, `run`, `Physics::solve`, descriptors, its goldens and tests. `nyx-sms` has four studies: `static`, `eigenfrequency`, `axisymmetric_eigenfrequency` (bodies of revolution on the (r, z) half plane, axis conditions, cylindrical absorbing layers; ADR 68, 69) and `modes` (2.5D elastic waveguide modes `exp(-j q y)` on the 2D (x, z) cross-section; ADR 101); every cross-section mode carries numeric frequency parts and energy shares per component and region (ADR 80, 101). |
| `crates/nyx-analysis` | analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a `SolveResult` and its mesh (phase 5, ADR 56). A mode study states its frame in the scalar `field frame` (`nyx_physics::FieldFrame`) and carries the complex field in its cell fields; results carry no analysis; it also draws the modes of a result as a figure, SVG or PNG (`figure` command, ADR 99) |
| `crates/nyx-store` | workspace, workflow, experiment and run storage (JSON specs; no dependency on another crate of the workspace, ADR 59). Not to be confused with `nyx-storage` |
| `crates/nyx-storage` | the pure-Rust HDF5 results file engine (`hdf5-pure` + `zlib-rs`, packed layout; ADR 107) |
| `crates/nyx-install` | installer and packager: `nyx-setup` (GUI installer), `nyx-uninstall`, `nyx-pack` (`.nyxpack`, zip, dmg; ADR 122) |
| `crates/nyx-engine` | `Engine`: the project, commands (every action by name, with JSON Schemas; the window and the scripts use the same commands; the MCP server has four script tools, ADR 94), the job runner, the geometry operation descriptors, the materials of a session (library plus the project's own, ADR 58), analyses, the guide cards, store wiring, and the registry of compiled-in physics (cargo features `dfs sms ems`). `facade` re-exports what a front-end draws with. Examples (`examples/`) are benchmark and measurement programs (`bench`, `profile3d`, `mem_profile`, `cavity_bench`, `stage_log`, `eig_beam`, `isolate_cost`, `p3measure`, `peakrun`, `latency`, `workflow`, `squadds_bench`) |
| `crates/app` | the egui application (package and binary `nyx`; on Windows a GUI-subsystem program that opens no console, its command-line modes join the terminal that started it, ADR 123; the workflow window is one bar and five tabs with the viewport always on, the script a dock, the physics one dropdown, ADR 124): a front-end of `nyx-engine`; the Study, Physics and Materials tabs are drawn from the physics descriptors (`forms.rs`) |
| `crates/nyx-script` | scripts, a third front-end of `nyx-engine` (ADR 93): Python syntax (`*.nys`, not `.py`, which is git-ignored) run by Monty in a **worker process** (`nyx-script-worker`, built from this package, must sit next to `nyx`; `NYX_SCRIPT_WORKER` names it); its `kernel` is one session that outlives a script (ADR 94). `nyx check FILE.nys`, `nyx run FILE.nys`, `nyx api stubs`. The engine's path guard (`nyx-engine::path_guard`) keeps a script's file arguments inside one folder. `monty-*` crates are pinned with `=` and `Cargo.lock` pins `salsa` to 0.28.2 |
| `crates/nyx-mcp`, `nyx-connect` | the AI interface (stdio MCP server with four script tools, `script_run` `script_wait` `script_check` `api`, over one script kernel of `nyx-script`: no tool per command, ADR 94, amends ADR 57; it needs `nyx-script-worker` next to it, so `cargo build --release -p nyx-script` before `cargo test -p nyx-mcp`; one-click host connection) |
| `examples/`, `materials/` | example projects (`*.nyx.json`, written by `crates/nyx-engine/tests/make_examples.rs`) and the material library (`*.toml`), both found next to the program or in a parent folder |
| `docs/` | `HOW_IT_WORKS.md`, `MATERIALS.md`, `HOST_SETUP.md`, the generated command reference `api/`, `architecture/`, `vvuq/`, `scripting-samples/`, reports; `history/` holds earlier plans, proposals and status logs (read-only, not current) |
| `docs/architecture/` | C4 model (`workspace.dsl`, matches the crate layering test) and the C4 set in `c4/` (the spec, model fragments and tools are tracked; the component and code pages are generated and git-ignored: `bash docs/architecture/c4/tools/build.sh`, needs Perl, Java and the Structurizr CLI), the ADRs (`adr/`, 119 records up to ADR 133, index `ADR_INDEX.md`; some numbers are reserved), `CODE_TO_ADR.md` (every tracked file has a row), interface `contracts/`. Read `README.md` there first |
| `tools/` | `baseline.sh` (tests and benchmark timings), `vvuq_ci.sh` / `vvuq_ci.ps1` (the CI tiers; the workflow runs them), `vvuq_telemetry.py` (CI telemetry, the one tracked Python file), `package_dist.sh`, `bench_mesh3d.sh`, `ux-shots.sh`, `workflows/` |
| `.github/workflows/` | `vvuq-ci.yml` (Tier 1 on every PR, more on `main`, tags and nightly) and `release.yml` (binaries for Windows, Linux and macOS on a `v*` tag) |
| `runs/`, `dist/`, `release-out/` | local outputs, git-ignored; never commit results |

Dependency rules are in `crates/nyx-engine/tests/layering.rs` (a table; the list of temporary edges is empty since phase 7). Nothing outside the physics crates names a physics. The window depends on `nyx-engine`, `nyx-script` and `nyx-connect` (the Connect menu); the MCP server on `nyx-engine` and `nyx-script`.

## Build and test

```
cargo run --release -p nyx                                  # the application
cargo test --release --workspace                            # everything; a few minutes
cargo test --release --workspace -- --ignored               # slow tests (solver library); `-p nyx-engine --test make_examples -- --ignored` rewrites examples/*.nyx.json
cargo clippy --release --workspace --all-targets --no-deps -- -D warnings   # the whole workspace is clean (ADR 65)
cargo fmt -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine   # only these crates are rustfmt-formatted (settings as geom); the others are left as written (ADR 65)
cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering   # interface contracts (ADR 19)
bash docs/architecture/check_code_to_adr.sh                 # every file has an ADR row
cargo test --release -p nyx-engine --test h5                 # results files (pure Rust, no HDF5 library needed, ADR 107; NYX_REQUIRE_H5PY=1 makes the h5py cross-check mandatory: pip install h5py numpy)
bash tools/baseline.sh <dir>                                # test counts and the five benchmark timings (BENCH_ONLY=1 for timings only)
```

Always use `--release`: the solvers are unusably slow unoptimised. The toolchain is pinned in `rust-toolchain.toml`. The build goes to `target/` inside the checkout; several checkouts on one machine can share one directory with `CARGO_TARGET_DIR`. If a long build or test is needed, run it in the background.

## Conventions

* Units are SI (metres) everywhere in the library and the MCP tools.
* **Wafer frame**: x lateral, z normal to the wafer, waves propagate along **+y** (`exp(-j beta y)`). `femcore::section` is the one place that knows axis numbers (`plane_tensor`, `local_tensor`, `LOCAL_AXES`, `wafer_point`); do not scatter axis permutations elsewhere. Solid kernels use the Voigt order (x, z, y).
* Physics carry three-letter codes chosen by Paul: `dfs` (diffusion), `sms` (solid mechanics), `ems` (EM waves). Brillouin (`sbs`) is removed (ADR 24). Name new physics the same way; a physics is a crate behind the `Physics` trait (ADR 25).
* Golden files (`crates/nyx-{dfs,sms,ems}/tests/golden/*.json`, `crates/geom/golden/`) are the numerical oracle; the harness is `nyx_physics::testkit`. `NYX_WRITE_GOLDEN=1` rewrites them; do it only for a kernel known to be right, and say so in the commit. Never regenerate a golden to hide a difference (ADR 10).
* **Tests never compare JSON text or key order** (ADR 26): `serde_json` keeps insertion order when any crate of the build enables `preserve_order` and sorts otherwise, so `--workspace` and `-p` builds differ. Compare `Value`s or `nyx_physics::canonical_json`.
* Contract exemplars (`docs/architecture/contracts/`) change only with a line in `VERSIONS.md` (ADR 19). `NYX_WRITE_CONTRACT=1` rewrites them.
* New Rust files start with the licence header (`// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial`, `// Copyright (C) 2026 p-j-r`; `geom` is BSD-2-Clause, see `NOTICE`) and a short `//!` comment saying what the file does. Some older files in `geom` and `femcore` have no licence header yet; add it when you touch them.
* `*.py` is git-ignored on purpose: the product ships no Python (scripts are `*.nys`). The one tracked exception is `tools/vvuq_telemetry.py` (CI, `!` rule in `.gitignore`). Any GDS layout is built with gdsfactory, never raw gdstk.
* Determinism: the same input gives the same numbers for any thread count. `NYX_THREADS=N` sets the rayon pool, `NYX_PROFILE=1` echoes the solver log.
* **Curved geometry** (ADR 82, 83): the mesh setting `curved` (`set_mesh`, default off) gives the 3D cells on a curved surface quadratic geometry (`CellMesh::edge_nodes`, `femcore::fem::CurvedTet`); only the `ems` `cavity` study uses it, the other physics see straight cells. With it off every number is the one of the straight mesh. `add_size` on a face and `max_dev` act in the 3D lattice mesher.
* **File format** (ADR 63): `Project`, `FeatureList` and `ProblemSpec` carry `format` (value 1); files without it or of another format are refused, nothing is converted. Project and spec files are JSON (`*.nyx.json`, `experiment.nyx.json`); results and meshes are HDF5 (feature `hdf5`), never JSON.

## State of the checks

CI is the record: `.github/workflows/vvuq-ci.yml` (Tier 1 on every pull request: `cargo fmt --check` for the seven formatted crates, clippy with warnings denied, the contract, layering, API-docs and h5 tests, the ADR check; Tier 2 and 3 on `main`, tags and nightly). Test counts are not tracked in this file. The toolchain is Rust 1.99.1. The earlier numbers, and the timings of each restructure phase, are in `docs/history/RESTRUCTURE_STATUS.md`. Two SOI-wire goldens pass at the G2 tolerances (5e-5 relative on n_eff, 1e-3 on field maxima; `nyx-ems/tests/golden`); a geom area property test can be flaky. Do not regenerate a golden to hide a difference.

## Version control

* The repository is `p-j-r/nyx` on GitHub; `main` is protected by CI, not by convention. Work on a branch (`<topic>`), open a pull request against `main`, and merge it when CI is green. Do not push to `main` directly.
* Commit messages say what changed and why in one line, then details if needed. Do not leave work only in an uncommitted working tree.
* Before a pull request: `cargo fmt --check` for the formatted crates, clippy, the `contracts` and `layering` tests, and `bash docs/architecture/check_code_to_adr.sh` (a new or moved file needs a row in `docs/architecture/CODE_TO_ADR.md`). Decisions that change the picture get an ADR (see `docs/architecture/README.md`).
* Where the code conflicts with a doc, the code and the ADRs are right: fix the doc. Where an instruction here conflicts with what you find, ask paul.
