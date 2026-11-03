# nyx

Desktop finite-element application for photonics, phononics and stimulated Brillouin scattering (SBS), in pure Rust. No Python is needed to build or run anything. Read `README.md` for the feature tour and `docs/HOW_IT_WORKS.md` for what the solver does and what each log line means.

## Layout

| Path | What |
|---|---|
| `crates/femcore` | solver library: meshes, sparse and eigen solvers, physics modules, use cases |
| `crates/geom` | geometry kernel (parametric wafer-frame models, selectors, golden files) |
| `crates/nyx-store` | workspace, workflow, experiment and run storage |
| `crates/app` | the egui application (binary `nyx`) |
| `ai/` | a **separate** Cargo workspace: `nyx-session`, `nyx-mcp` (stdio MCP server), `nyx-connect`. Depends on `femcore` and `geom` by path. See `ai/README.md` |
| `docs/` | design notes and per-project status files (`KERNEL_PORT_STATUS.md`, `FRAME_UNIFICATION_STATUS.md`, ...) |
| `docs/architecture/` | C4 model (`workspace.dsl`), ADRs (`adr/`, index `ADR_INDEX.md`), `ROUGH_SPOTS.md`, `CODE_TO_ADR.md`, interface `contracts/`. Read `README.md` there first |
| `runs/`, `dist/` | local outputs, git-ignored; never commit results |

## Build and test

```
cargo run --release -p nyx                                  # the application
cargo test --release --workspace                            # fast suite, a few minutes
cargo test --release --workspace -- --ignored               # slow tests (Brillouin 30-60 s each, Bragg ~1.5 min)
cargo clippy --release --workspace --all-targets -- -D warnings
cd ai && cargo test --release                               # MCP server, own lockfile and target
cargo test --release -p femcore -p geom -p nyx-store --test contracts   # interface contracts (ADR 19); do not need the golden files
bash docs/architecture/check_code_to_adr.sh                 # every file has an ADR row
```

Always use `--release`: the solvers are unusably slow unoptimised. The toolchain is pinned in `rust-toolchain.toml`. The first build compiles egui, wgpu and faer and takes several minutes. If a long build or test is needed, run it in the background.

## Conventions

* Units are SI (metres) everywhere in the library and the MCP tools.
* **Wafer frame**: x lateral, z normal to the wafer, waves propagate along **+y** (`exp(-j beta y)`). `femcore::section` is the one place that knows axis numbers (`plane_tensor`, `local_tensor`, `LOCAL_AXES`, `wafer_point`); do not scatter axis permutations elsewhere. Solid kernels use the Voigt order (x, z, y).
* Physics modules carry three-letter codes chosen by Paul: `dfs` (diffusion), `sms` (solid mechanics), `ems` (EM waves), `sbs` (Brillouin). Name new modules the same way.
* Golden files (`crates/femcore/tests/golden_frames/*.json`, `crates/geom/golden/`) are the numerical oracle. `NYX_WRITE_GOLDEN=1` rewrites them; do it only for a kernel known to be right, and say so in the commit.
* Every `femcore` source file starts with a short comment saying what it does; keep that when adding files.
* `*.py` is git-ignored on purpose: nothing Python ships. Any GDS layout is built with gdsfactory, never raw gdstk.
* Determinism: the same input gives the same numbers. `NYX_THREADS=N` sets the rayon pool, `NYX_PROFILE=1` echoes the solver log.

## Known state of the checks (2026-11-03, Windows, Rust 1.99.1)

* `cargo test --release --workspace --no-fail-fast`: 632 pass, 10 ignored, **7 fail**, all in `femcore --test golden_frames` (`golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`, `golden_brillouin_wire_in_a_cubic_crystal`, `golden_brillouin_wire_in_a_tilted_crystal`, `golden_brillouin_forward_scattering_in_a_tilted_crystal`, `golden_ring_brillouin_cubic_crystal`, `golden_ring_brillouin_tilted_crystal`). The differences are small (about 1e-6 relative in fields, up to 2e-5 in one effective index), identical with `NYX_THREADS=1`, and present at the commit that wrote the golden (35cceba), so they come from the machine the goldens were frozen on, not from a later change. Do not regenerate the goldens to hide them; decide the reference platform or the tolerance first.
* `cargo clippy ... -D warnings` is **not** clean with 1.99.1: 8 errors in `femcore` (`ems_pml.rs`, `sms.rs`, `sms_kz.rs`: needless borrows, a redundant closure, an `unwrap` after `is_none`).
* `cargo fmt --check` fails widely (the code was never rustfmt-formatted; only `crates/geom` has a `rustfmt.toml`).
* `ai/` tests pass.

## Version control

* One branch per task, named for the task; never commit to `main` directly. Commit messages say what changed and why in one line.
* `docs/*_STATUS.md` files are the handoff between sessions: read the relevant one before starting, update it before finishing.
* A restructure is planned but not started: read `docs/architecture/RESTRUCTURE_PLAN.md` and `docs/RESTRUCTURE_STATUS.md` before changing crate layout, mesh types, solvers or the problem types.
* Do not leave work only in an uncommitted working tree or a side folder; commit it to a branch.
