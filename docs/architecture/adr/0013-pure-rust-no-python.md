# 13. Pure Rust; nothing Python ships; HDF5 is optional; the toolchain is pinned

Date: 2026-10-28

## Status

Accepted (`README.md`, `CLAUDE.md`, a047506).

## Context

nyx is a desktop application for researchers on Windows. Python and C dependencies make it hard to build, ship and reproduce. Two parts of the ported work were developed against Python or C tools (the SLEPc comparison harness, HDF5).

## Decision

* The application, library and MCP server are pure Rust. `*.py` is git-ignored on purpose. Any GDS layout is built with gdsfactory, never raw gdstk (Paul's rule; this is a layout-tool rule, not a build dependency).
* HDF5/XDMF results are behind the `hdf5-metno` feature of `femcore` (`hdf5` in `app`; `phdf5` adds MPI-IO). They need the HDF5 library; the default build does not.
* `faer` supplies dense/sparse linear algebra, `spade` the Delaunay triangulation, `quick-xml` the project XML, `serde`/`serde_json` the data model, `rayon` the pool, `egui`/`eframe` the window. Projects save as JSON or XML (`.nyx.json`, `.nyx.xml`) from the same serde value (`xmlio`).
* `rust-toolchain.toml` pins the compiler; `Cargo.toml` says `rust-version = 1.84`. Always build with `--release`: the solvers are unusably slow unoptimised.

## Consequences

* The state of the checks depends on the toolchain: clippy `-D warnings` is not clean with 1.99.1 and `cargo fmt --check` fails widely (R23).
* XML is a derived encoding of a serde value, not a schema of its own: field names of the Rust types are the file format, and nothing versions them (R20).
* The ports' Python harnesses live outside this repository.

## Code

`Cargo.toml`, `rust-toolchain.toml`, `crates/femcore/Cargo.toml`, `crates/femcore/src/xmlio.rs`, `h5io.rs`, `crates/app/src/batch.rs`.

## Contract

`contracts/problem_spec` exercises the JSON and XML encodings of the `ProblemSpec` part of a project; a contract for the whole `Project` is not written (it lives in the `app` crate, R16).
