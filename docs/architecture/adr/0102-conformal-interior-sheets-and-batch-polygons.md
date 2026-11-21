# 102. Conformal interior 2D sheets (`Sheet`) and batch polygons (`Polygons`)

Date: 2026-11-21

## Status

Accepted (superconducting qubit milestone M1). Extends the geometry feature list of ADR 2, the 3D mesher of ADR 40 / ADR 55, and the boundary integral reduction of ADR 90.

## Context

Planar superconducting quantum circuits combine $100\ \text{nm}$–$10\ \mu\text{m}$ coplanar gaps, transmon pads, Josephson junction straps, and flux-trap holes at the $z = z_0$ metal interface inside a millimetre-scale 3D chip box. Modelling $35\text{–}200\ \text{nm}$ thick metal as a 3D solid forces extreme aspect ratios or millions of tetrahedra, and an un-snapped background lattice cannot place conforming triangles across zero-thickness planar polygons. Meanwhile, importing multi-island layouts (`gdsfactory` / KLayout polygons) per individual polygon creates excessive feature-list boilerplate and risks opening sub-tolerance cracks at shared vertices during simplification.

## Decision

* **Batch polygon feature (`geom::model::Op::Polygons` / engine `add_polygons`).** Accepts `polygons` (`Rings(2)` outer loops), optional `holes`, `simplify_tol` (default `1e-9` m), and `min_hole_area` (default `0.0` m²). `geom::poly_batch` cleans consecutive duplicates, filters undersized holes, and applies Ramer–Douglas–Peucker simplification while pinning vertices within `pin_tol` of any sibling ring vertex so abutting pads and junction straps never split apart.
* **Interior 2D conductive sheet (`geom::model::Op::Sheet` / engine `add_sheet`).** Embeds a 2D planar feature or boolean tree (`Dims::Three`) at plane `z` (`Expr`) inside a 3D layer stack (`geom::stack`). Adjacent 3D layers sharing `z_0` stay electrically continuous across uncovered gaps, while facets covered by the sheet receive the sheet boundary label (`"<sheet>"`) and its 1D boundary loop receives `"<sheet>:edge"` for localized edge refinement (`size`) and 3D TLS edge-tube regularization.
* **Conformal 3D sheet-stack mesher (`femcore::meshing::sheet3d`).** When a 3D model contains any `Op::Sheet`, `femcore::meshing::mesh_model` dispatches to `sheet3d`:
  1. Projects all 3D layers and 2D sheets onto a single conformal 2D planar straight-line graph (PSLG) and triangulates via `spade` CDT + Ruppert refinement honoring face, region, and `"<sheet>:edge"` size rules with zero T-junctions.
  2. Extrudes graded horizontal z-slabs matching local element sizes at every sheet plane `z_s` and growing smoothly (`1.35x`) toward bulk layers.
  3. Splits each triangular prism into three tetrahedra along deterministic lowest-global-vertex-id (Dompierre) diagonals so shared triangular faces conform bit-identically.
  4. Emits each interior boundary facet (`"<sheet>"` or material interface `"rb|ra"`) exactly once into `SimplexMesh::facets`.
* **Interior-surface `integrate` reduction.** `nyx-engine::reduce` allows `integrate(boundary=..., side=...)` on interior sheet/interface boundaries so surface charge, energy, and TLS integrands can be queried directly from the engine.

## Consequences

* 3D planar superconducting chips with sub-micron junctions and micrometre CPW gaps conform cleanly with thousands—rather than millions—of tetrahedra.
* `geom_operations.v1.json` and `engine_commands.v1.json` gain `Polygons` / `add_polygons` and `Sheet` / `add_sheet` (additive).

## Code

`crates/geom/src/{poly_batch,model,complex,planar,scene,select,stack,topo3,print,lib}.rs`, `crates/femcore/src/meshing/{sheet3d,mod,request}.rs`, `crates/femcore/src/mesh/{cells,mod}.rs`, `crates/nyx-physics/src/desc.rs`, `crates/nyx-engine/src/{geomdesc,geometry,reduce}.rs`, `crates/app/src/geomui.rs`; tests `crates/geom/tests/stack3d.rs`, `crates/femcore/tests/sheet_mesh3d.rs`.
