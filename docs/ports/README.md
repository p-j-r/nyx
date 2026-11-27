# Ported solver and mesher libraries

The solver and mesher code of `femcore` (today `femcore::linalg`, `femcore::meshing` and `femcore::adapt`) comes from two staging crates, `femsolve` and
`femmesh`, written by separate working sessions without access to nyx. This folder keeps the provenance and licence texts; their staging notes (status, handoff, integration plans, measurements) are archived in [`../history/ports/`](../history/ports) and describe modules under their old names (`femcore::solvers`, `femcore::meshers`).

| Folder | Module | Start here |
|---|---|---|
| [`solvers/`](solvers/) | `femcore::linalg` | `PROVENANCE.md`, `THIRD_PARTY_NOTICES.md`, `LICENSE` |
| [`meshers/`](meshers/) | `femcore::meshing`, `femcore::adapt` | `PROVENANCE.md`, `REFERENCES.md`, `THIRD_PARTY_NOTICES.md`, `LICENSE` |

[`NOTICES.md`](NOTICES.md) is the repository-wide origin and licence summary. **Both ports are BSD-2-Clause**; the licence files and the
PETSc/SLEPc notices are next to the provenance files. `README.md` and `NOTICE` state the carve-out (these parts keep BSD-2-Clause and are not relicensed); whether that is compatible with nyx's own dual licence has not been reviewed by a lawyer, so do that before a commercial release.

## What changed in the move (as of the merge; module names then)

*Historical: the sections below describe the merge when it happened, with the module names of that time.*

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
