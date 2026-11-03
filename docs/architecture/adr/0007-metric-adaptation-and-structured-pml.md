# 7. Metric-based mesh adaptation and structured PML layers, ported as `femcore::meshers`

Date: 2026-11-02

## Status

Accepted. The port (`7a9a85e`) was written by a separate session without access to nyx (see `docs/ports/meshers/HANDOFF.md`); 2D integration `3fbb24d`, 3D `44db967` to `fdb9f55`.

## Context

Hand-tuned meshes limited the accuracy of the chips and ring studies. Metric-based adaptation (Riemannian metrics from recovered Hessians or from physics) with local operators (split, collapse, flip, smooth) gives meshes sized to the field. Absorbing layers were unstructured triangles of uneven quality.

## Decision

* The mesher port lives in `femcore::meshers` as its own dimension-generic stack (`D = 2` or `3`): metrics (`adapt/metric.rs`), geometry oracles (`geom`, `geom3`), working meshes (`topo`, `topo3`), operators (`ops`, `ops3`), drivers (`drive`, `drive3`), PML generator (`adapt/pml.rs`), quality statistics, fixtures and a small independent Helmholtz validation harness (`validate/`).
* Constrained features (boundaries, material interfaces, corners, frozen faces) are classified from the labels and preserved exactly; periodic pairs stay congruent through frozen faces.
* Structured PML layers are generated around a box-shaped physical mesh: quads (2D) and prisms/hexes (3D), conforming, with `meshers::interop::pml_domains` producing the `Pml` description per layer label, and materials copied from the interface (`pml3`, `meshing_pml`).
* The **solve-estimate-adapt loop** (`estimate`, `estimate3`) recovers Hessians from a solution and feeds `adapt_to_metrics`; at order 1 it does not beat uniform meshes (Maxwell residual estimator needed), at order 2 it does.
* Everything is deterministic and bit-identical across thread counts (ADR 11).

## Consequences

* Adaptation is 2D and serial in the original port; 3D adds edge removal, flat-boundary flips and surface/crease collapses, still serial, with linear-time closest-point queries (`MISSING_FEATURES.md` section 3).
* PML reflection has been measured only for the 2D scalar problem; for the vector solvers only mesh convergence of damped modes and lossless equivalences are tested.
* Attempts at 30 split passes in a row degrade slivers; flips are interleaved. These lessons are in `ports/meshers/HANDOFF.md` section 9 and apply to any future threaded version.
* Provenance and licence (BSD-2-Clause) are recorded in `docs/ports/meshers/PROVENANCE.md` and `THIRD_PARTY_NOTICES.md`; whether BSD-2 is compatible with nyx's own licence has not been checked.

## Code

`crates/femcore/src/meshers/**`, `meshing.rs`, `meshing_pml.rs`, `pml3.rs`, `estimate.rs`, `estimate3.rs`, `refine3d.rs`. History: `7a9a85e`, `3fbb24d`, `154b64e`, `b9a135e`, `44db967`, `c49c4aa`, `281b2f8`, `79936d6`, `62d6cdc`, `ba625ec`, `fdb9f55`.

## Contract

`contracts/mesh`. The mesher's own acceptance tests are `tests/meshers.rs`, `adapt3d.rs`, `meshing.rs`, `physics_meshing.rs`.
