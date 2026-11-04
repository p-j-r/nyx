# 54. A mesh is asked for and reported in data

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). New contracts `mesh_request`, `mesh_report` and their schemas (ADR 19). Nothing in the engine or the MCP crates is rewired here; this is what they can call.

## Context

Adaptive meshing, wavelength sizing and structured absorbing layers were Rust structs (`MeshPlan`, `Adaptive`, `PhysicsSizing`) that no command could name, and the engine answers "unsupported" for physics sizing. The engine's commands and the MCP tools are generated from JSON Schemas (ADR 37, 38); the mesh needed a request and a report of the same kind, for 2D and 3D.

## Decision

`femcore::meshing::request`:

* **`MeshRequest`** (serde): `h`, optional `max_dev`, `mesher` (`{"kind": "simple"}` or `{"kind": "adaptive", scale, grade, max_cells}`), `structured_pml`, and optional `waves` (`sizing` = the neutral `Sizing` of a physics, `points_per_wavelength`, `pml_relax`, `eigen_headroom`). `plan()` gives the `MeshPlan`.
* **`mesh_request(list, &request)`** meshes either dimension (`mesh_for_sizing` when there are waves, the plain entry point otherwise) and returns the `MeshedGeom` and a **`MeshReport`**: dimension, vertices, cells, facets, per domain the number of cells and the wave size it asked for, a summary of the adaptation (iterations, converged, capped, worst quality, cells at the start and at the end), and the layers per side of structured layers.
* **JSON Schemas** (draft 2020-12) `request_schema()` and `report_schema()`, checked in as `mesh_request.schema.v1.json` and `mesh_report.schema.v1.json`; exemplars `mesh_request.v1.json` (compared as a value) and `mesh_report.v1.json` (compared by structure: the counts belong to the mesher, not to the contract). `additionalProperties` is false: a typo in a request is an error.

## Consequences

The engine can add a `mesh` argument that is this schema and report this struct without knowing a mesher; the 3D features are reachable through it. The window keeps its own `MeshPlan` use.

## Code

`crates/femcore/src/meshing/request.rs`.

## Contract

`mesh_request_and_report_match_their_exemplars_and_schemas` (`tests/contracts.rs`; a small validator of the schema subset used is part of the test file, no schema crate is needed).
