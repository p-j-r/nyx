# 16. `.fmsh`: a lossless, self-checking binary format for meshes and results

Date: 2026-10-28

## Status

Accepted by the mesher port (`docs/ports/meshers/STATUS.md` "Decisions I made"). The format has no consumer in the app yet; HDF5/XDMF (`h5io`) remains the app's results format.

## Context

The original mesher stack used HDF5 plus a small XML descriptor. HDF5 needs the C library or an unproven pure-Rust writer; VTK and text writers produce large files; the requirement was no loss of information.

## Decision

`.fmsh` (version 2) is one self-describing, CRC-checked file holding a mixed-cell mesh (cell types, labels, boundary facets) and any number of result steps with named point or cell fields. Floats round-trip bit for bit. Each array is stored with the smallest of byte planes, XOR-with-previous or dictionary coding, then an adaptive range coder (order 0, order 1 or raw), implemented in `io/entropy.rs` with no dependency. Sizes are 12 to 30 times smaller than plain binary on structured and PML meshes, 2.1 times on an adapted mesh, and within 7 percent of `xz -9e`. A lossy quantised option was deliberately not added.

## Consequences

* Not wired to the MCP plan's "mesh not stored in v1; rebuilt deterministically" (a row of the session table in `MCP_ASSISTANT_PLAN.md`) or to `nyx-store`'s `Meshes` share level, which could use it.
* The connectivity of adapted meshes (5.5 bytes per triangle) could drop about 4 times with a topological coder; not done.
* Three results/mesh encodings exist: `.fmsh`, HDF5/XDMF (`h5io`, attr `version = 2`) and the project XML/JSON; `.fmsh` (version byte 2) and the HDF5 file (attribute `version = 2`) record a version, the project encodings do not.

## Code

`crates/femcore/src/meshers/io/{fmsh,entropy}.rs`, `h5io.rs`.

## Contract

`contracts/mesh` includes an `.fmsh` round-trip and a pinned version byte (`mesh_fmsh_round_trip_and_version`).
