# AGENTS.md

Agent-neutral guide for coding agents working in the `nyx` repository (`nyx2` / `nyx3`). For human-facing overview, build background, and physics summaries see [`README.md`](README.md) and [`docs/HOW_IT_WORKS.md`](docs/HOW_IT_WORKS.md); read [`docs/architecture/RESTRUCTURE_PLAN.md`](docs/architecture/RESTRUCTURE_PLAN.md) and [`docs/RESTRUCTURE_STATUS.md`](docs/RESTRUCTURE_STATUS.md) before changing crate layout, mesh types, solvers, or problem types.

## Project overview

`nyx` is a pure-Rust desktop finite-element application and simulation engine for photonics, phononics, electrostatics/heat diffusion, and electromagnetic cavities/waveguides. No Python is required to build or run `nyx`. One `nyx-engine` crate drives three front-ends identically through JSON-schema-validated commands: the `egui` desktop app (`nyx`), the Monty-backed scripting runtime (`*.nys` via `nyx-script-worker`), and the four-tool stdio MCP server (`nyx-mcp`: `script_run`, `script_wait`, `script_check`, `api`).

## Workspace layout and layering

| Path | Responsibility |
|---|---|
| `crates/geom` | Parametric wafer-frame geometry kernel, selectors, arrangement, goldens (`serde` + optional `rayon` only) |
| `crates/femcore` | FEM core: `Mesh` / `SimplexMesh`, 2D/3D meshing & adaptation, sparse/eigen/AMG linear algebra (`linalg::driver`), finite elements, materials (`material`) & TOML reader |
| `crates/nyx-physics` | `Physics` trait, `ProblemSpec` envelope, descriptors, `SolveResult`, `Registry`, validation, JSON Schemas (ADR 25) |
| `crates/nyx-dfs`, `nyx-sms`, `nyx-ems` | Independent physics crates (`dfs` diffusion/electrostatics, `sms` solid mechanics, `ems` EM waves), each implementing 3-stage kernels and descriptors |
| `crates/nyx-analysis` | Post-solve analyses over `SolveResult` + `Mesh` (Stokes/Poincare polarization, SVG/PNG mode figures) |
| `crates/nyx-store` | Workspace, workflow, experiment, and run storage (JSON manifests; no workspace dependencies) |
| `crates/nyx-engine` | `Engine`, command catalog & schemas, geometry descriptors (`geomdesc`), materials, analyses, job runner/queue, sweeps, HDF5 I/O, guide cards |
| `crates/nyx-script` | `*.nys` script execution in isolated worker (`nyx-script-worker`) over `Engine::call` |
| `crates/nyx-mcp`, `nyx-connect` | Four-tool stdio MCP server (`nyx-mcp`) and AI host configuration (`nyx-connect`) |
| `crates/app` | `egui` desktop shell (`nyx` binary) over `nyx-engine` |
| `examples/`, `materials/`, `docs/`, `tools/` | Project examples (`*.nyx.json`), TOML material library (`format = 1`), ADRs/contracts/C4 docs, `baseline.sh` |

Dependency layering is enforced by `crates/nyx-engine/tests/layering.rs`:
- `geom` depends only on `serde` and optional `rayon`.
- `femcore` depends on no `nyx-*` crate (only `geom` in the workspace).
- Physics crates (`nyx-dfs`, `nyx-sms`, `nyx-ems`) and `nyx-analysis` depend only on `nyx-physics` and `femcore` — never on each other.
- Only `nyx-engine` may depend on the physics crates and `nyx-analysis`.

## Dev environment, build, and test commands

The Rust toolchain is pinned in [`rust-toolchain.toml`](rust-toolchain.toml) (Rust `1.99.1` with `clippy` and `rustfmt`). **Always build and test with `--release`**: solvers are unusably slow unoptimized. Set `CARGO_TARGET_DIR` (or keep `.cargo/config.toml`) so concurrent `jj` workspaces share external registry builds cleanly without polluting workspaces.

```powershell
# Run desktop app or batch/script commands
cargo run --release -p nyx
cargo test --release --workspace                    # full suite
cargo test --release --workspace -- --ignored       # include slow solver-library tests

# nyx-mcp tests require nyx-script-worker next to the test binaries:
cargo build --release -p nyx-script
cargo test --release -p nyx-mcp

# HDF5 results tests (hdf5 feature is on by default; requires HDF5_DIR and bin on PATH)
cargo test --release -p nyx-engine --test h5

# Linting and selective formatting (ADR 65: format ONLY these 7 crates)
cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
cargo fmt -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine

# Architecture, layering, and contract checks
cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering
bash docs/architecture/check_code_to_adr.sh
bash tools/baseline.sh <dir>
```

## Coding, numerical, and contract conventions

- **Units and wafer frame:** All library and engine interfaces use SI units (metres, Hz, F, H, S/m). In the wafer frame, `x` is lateral, `z` is the wafer normal, and `+y` is the propagation direction (`exp(-j beta y)`). Keep axis mappings centralized in `femcore::section`; solid kernels use Voigt order `(x, z, y)`.
- **Three-stage kernels & curved geometry:** Physics kernels split cleanly into assembly, solve, and field extraction (`ADR 9`). Curved edge nodes (`CellMesh::edge_nodes`, `ADR 82`) apply to `ems` `cavity` when `curved` is enabled; straight meshes must preserve exact baseline behavior.
- **File headers & licensing:** New Rust files in `femcore` and `nyx-*` crates must start with the dual-license SPDX header (`// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial`, `// Copyright (C) 2026 p-j-r`) and a `//!` module summary comment.
- **Determinism:** Identical inputs must produce identical outputs regardless of thread count (`NYX_THREADS=N`; fixed-chunk parallel reductions).
- **Testing & JSON rules:**
  - Never compare raw JSON strings or key ordering (`ADR 26`); compare parsed `serde_json::Value` or `nyx_physics::canonical_json`.
  - Golden files (`crates/nyx-{dfs,sms,ems}/tests/golden/`, `crates/geom/golden/`) are the numerical oracle. Never regenerate a golden file (`NYX_WRITE_GOLDEN=1`) to hide a regression.
  - Additive or breaking interface/schema/command changes require updating contract exemplars (`NYX_WRITE_CONTRACT=1`) and bumping [`docs/architecture/contracts/VERSIONS.md`](docs/architecture/contracts/VERSIONS.md) in the same commit (`ADR 19`).
  - Never commit generated `.py` files (`*.py` is git-ignored; ship `*.nys` scripts). Any external chip/photonic GDS layout helpers use `gdsfactory` / `klayout.db`, never raw `gdstk`/`gdspy`.
- **Architecture tracking:** Every tracked file must match a row in [`docs/architecture/CODE_TO_ADR.md`](docs/architecture/CODE_TO_ADR.md) (verified by `docs/architecture/check_code_to_adr.sh`). New architectural decisions receive sequential ADRs in `docs/architecture/adr/` (note: `70..79` and `81` are reserved) and entries in `docs/architecture/ADR_INDEX.md` and `docs/architecture/CODE_TO_ADR.md`.

## Version control (`jj` multi-agent workflow)

- This repository uses colocated Jujutsu (`jj`) with **no remote: never push**.
- Never commit directly to `main` or modify another session's working copy (`default@`, `nyx3@`, etc.). Wait for explicit user acceptance before advancing `main`.
- Every coding agent, subagent, or swarm worker modifying code or running mutating builds must operate in its own dedicated `jj` workspace (e.g., under `.jj-workspaces/<name>` with `.jj-workspaces/` excluded locally in `.git/info/exclude`), describe its commits cleanly (`jj describe -m "..."`), and leave clean verified states for integration.
