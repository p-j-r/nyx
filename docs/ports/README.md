# Ported solver and mesher libraries

`femcore::solvers` and `femcore::meshers` come from the `inspired-rust-port` repository (two staging crates, `femsolve` and
`femmesh`, written by separate working sessions without access to nyx). They were moved into `femcore` as modules; this folder keeps
their documentation, provenance and licence texts.

| Folder | Module | Start here |
|---|---|---|
| [`solvers/`](solvers/) | `femcore::solvers` | `STATUS.md` (what was done and measured), `EXPLAINER.md`, `NYX_INTEGRATION.md` (the author's plan for swapping the solvers into the physics), `PROVENANCE.md` |
| [`meshers/`](meshers/) | `femcore::meshers` | `HANDOFF.md`, `STATUS.md`, `PROVENANCE.md` |

[`NOTICES.md`](NOTICES.md) is the repository-wide origin and licence summary. **Both ports are BSD-2-Clause**; the licence files and the
PETSc/SLEPc notices are next to the provenance files. Whether that is compatible with nyx's own licence has not been checked by
anyone; do that before publishing.

## What changed in the move

* Paths: `crate::x` became `crate::solvers::x` / `crate::meshers::x`; the tests became `tests/solvers/*` (one test binary) and
  `tests/meshers.rs`; the studies and the Matrix Market command-line tools (`fs_eps`, `fs_gen`, `fs_ksp`, `fs_contour`, `fs_snes`) are
  examples (`mesh_*` for the mesher studies). The PETSc/SLEPc comparison harness (C and Python drivers, benchmark results) was left
  out: nyx ships no Python. It stays in the original repository.
* Threads: the `parallel` feature is gone. The kernels always use rayon and run on nyx's pool (`femcore::par`, `NYX_THREADS`); with
  `NYX_THREADS=1` they run serially. Results are bit-identical for every thread count (test `solvers::determinism`).
* Edition: the ports are edition 2024, `femcore` is 2021, so three `if let ... && ...` chains were rewritten and two
  `is_multiple_of` calls (newer than the workspace's rust-version) replaced. Nothing else in the algorithms was touched.
* Meshes: the mesher keeps its own dimension-generic mesh types (they were stand-ins for the nyx types). `meshers::interop` converts to and
  from `crate::mesh::SimplexMesh` / `crate::cellmesh::CellMesh` and turns the PML regions of `build_pml` into `problem::Pml` domains.
* Not moved: `rsolvers` (the earlier real-only Krylov/Krylov-Schur port, superseded by `femsolve`).
* Feature `portable-kernels` of `femsolve` is now a `femcore` feature of the same name.

## What the merge did not do

The new libraries are **not wired into the physics yet**: `emw`, `solid`, `brillouin`, ... still call the old `eigen`, `lobpcg`, `amg` and
`linalg` modules. `solvers/NYX_INTEGRATION.md` lists the order in which to swap them and the checks to run (the author wrote it without
seeing nyx, so treat its nyx statements as unverified). The mesher's adaptation loop is 2D and serial and has not been run against
the nyx solvers.
