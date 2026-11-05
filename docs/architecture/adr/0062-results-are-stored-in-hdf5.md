# 62. Results and meshes are stored in HDF5; JSON is for small things

Date: 2026-11-05

## Status

Accepted on paul's instruction (HDF5 is central to how nyx stores the data of a solve). Restores the results file that phase 7 had dropped with the old result type (ADR 59 said it was a follow-up) and reverses the JSON fallback in the run record.

## Context

A solve produces fields with hundreds of thousands of values per solution and, for mode studies, complex fields per cell. As JSON text these are large, slow to write and parse, and lossy unless every digit is kept. The old `femcore::h5io` wrote them to HDF5; it was bound to the deleted `femcore::problem` types.

## Decision

* `nyx-engine::h5io` (cargo feature `hdf5`, off by default because it needs the HDF5 library at build time; set `HDF5_DIR`) writes **the project, the mesh and the latest result of every solved physics** (all solutions, all fields with their components and units, the per-cell fields, scalars and the stage log) into one file, chunked and deflate-compressed, and reads it back **bit for bit** (test `tests/h5.rs`). Layout version 3 is documented in the module and listed in `VERSIONS.md` (`results_file`). The parallel (MPI-IO) writer is not ported.
* Commands `export_results {path}` and `results_file_info {path}`; a build without the feature answers them with "unsupported", so a host can tell.
* **Every recorded run keeps `results.h5`** in its run folder when the feature is on (engine runs, MCP runs and the window's runs); no JSON copy of the fields is written in that case. Without the feature the engine falls back to `result_<physics>.json`, which is larger; builds that care about size or speed use the feature.
* Measured on two examples, the HDF5 file including the mesh and the project is 2.5 to 2.9 times smaller than the JSON of the result alone (SOI wire 487 kB against 1243 kB; thermal bar 413 kB against 1167 kB).
* **XDMF** (`nyx-engine::xdmf`, recovered from the old writer, XML written with `quick-xml`): next to every `.h5` that `export_results` writes and in every recorded run (`results.xdmf`), a description that names the nodes, cells and fields under `/viz/<physics>`; a scalar field is a hard link to the stored dataset, a vector field has a three-component copy. Each physics is a collection of grids, one per solution. A test resolves every reference against the file.
* Project and spec files, run records and the workflow manifests stay small JSON text. Whether the project file should be XML after all is open with paul (it was JSON by plan constraint 11).

## Consequences

* The window's "Open results file" is not back yet; `results_file_info` and `read_results` exist for it.
* Making `hdf5` a default feature would make the library a build requirement for every developer; that is a decision for paul.

## Code

`crates/nyx-engine/src/{h5io,xdmf,results_file,store_link}.rs`, `crates/nyx-engine/tests/h5.rs`, `crates/app/src/shell.rs`.
