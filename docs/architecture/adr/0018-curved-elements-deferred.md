# 18. Curved elements are deferred; the mesher projects onto exact geometry but stores no curved nodes

Date: 2026-11-02

## Status

Accepted as a constraint (paul, 2026-11-02, quoted in `docs/ports/meshers/FEATURE_DEPENDENCIES.md`). Open: node-storage convention, Jacobian checks and maximum order on the solver side.

## Context

Cells are straight-sided (constant gradients), so any curved boundary carries an O(h^2) geometric error whatever the element order; for P2 and Nedelec-2 this is the accuracy ceiling for discs, rings and spheres (`GEOMETRY_ARCHITECTURE.md` problem 6).

## Decision

* Curved nodes are **not stored yet**. Until the solvers define a convention, the mesher (a) places boundary and interface vertices on the exact geom edge or face (`GeomOracle2` in 2D, `geom3` in 3D) and (b) will do Jacobian validity checks and projection of mid-edge nodes onto the exact geometry when asked.
* The mesher assumes **no maximum element order** and assumes both Nedelec families.
* The curved-element refactor of the solvers belongs to the integration session, benchmarked against closed forms (`GEOM_HANDOFF.md` decision 6); geom's job is to make it possible.
* `matfree.rs` in the ported solver library (sum-factorised high-order operators for curved hexahedra) needs nyx's element geometry and is the intended consumer.

## Consequences

* `MISSING_FEATURES.md` section 1: no stored high-order nodes, no Jacobian checks for curved hex/prism/tet, no projection of mid-edge nodes. The 3D oracle can project but nothing places nodes.
* Smooth analytic surfaces (sphere, SDF projectors) are unit-tested but no end-to-end test adapts a curved body and checks convergence of a solver quantity.
* A format that stores curved nodes will change `SimplexMesh`/`CellMesh` and `.fmsh`; contracts `mesh` and `kernel_stages` are written so such a change fails loudly and is versioned (ADR 19).

## Code

`crates/geom/src/query.rs`, `curve.rs`; `crates/femcore/src/meshing.rs::GeomOracle2`, `meshers/adapt/geom.rs`, `geom3.rs`; `solvers/matfree.rs`.

## Contract

`contracts/mesh` states that meshes carry vertices only (test `mesh_has_no_curved_node_storage`).
