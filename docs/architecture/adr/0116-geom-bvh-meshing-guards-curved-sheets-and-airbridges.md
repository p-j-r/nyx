# 116. Geometry BVH acceleration, polymorphic selectors, stepped airbridges, sheet3d Z-slab/max_cells guards, Ruppert diagnostics, and curved sheet projection

Date: 2026-11-22

## Status

Accepted. Amends ADR 2, 3, 6, 52, 82, 102, and 114.

## Context

Large superconducting chips (6,000+ meander segments, multi-layer GDSII imports with transformed/patterned sheets, 1000:1 feature-size ratios, and suspended airbridges under structured PML) exposed quadratic geometry and meshing bottlenecks, silent Ruppert Steiner-cap termination, unguarded `sheet3d` tetrahedral blowup when fine lateral sizing coupled directly to Z-slab count, dropped interior sheets under `structured_pml`, and polyhedral error on circular planar pads.

## Decision

* **`geom` spatial acceleration, selectors, DAG pins, and airbridges (`crates/geom`):**
  * Added flat 2D AABB BVH (`crates/geom/src/bvh.rs`: `AabbBvh2`, `ComplexBvh`) rebuilt in curve `(u, v)` coordinates so `Complex::locate`, `point_in_loop_idx`, arrangement side/hole checks, and `stack` lateral matching run in $O(\log N)$ (`< 350 ms` on 6,000-segment meanders, bit-deterministic across `NYX_THREADS`), and added spatial hashing in `poly_batch::simplify_rings_with_pins`.
  * Made `Selector::InBox` polymorphic across `Regions`, `Faces`, `Edges`, and `Vertices` (including 3D `edges3`) with concrete-first evaluation in `Inter`/`Diff`/`Union`/`Not`.
  * Added recursive DAG pin propagation (`Builder::collect_pins`) through `Transform`, `Pattern`, boolean ops, and `Extrude` so transformed/patterned abutting sheets share pinned conformal vertices.
  * Added `Model::stepped_airbridge` / `SteppedAirbridge` lowering airbridges into prismatic `Op::Box` landing pads (`RegionKind::Hole`) + elevated `Op::Sheet` decks (`needs_planar == false`), while making `planar::build_scene` return an explicit error instead of silently dropping `Op::Sheet` features.
* **`femcore::meshing` guards, octave fill grid, PML preservation, and curved sheet projection (`crates/femcore/src/meshing`):**
  * Added optional `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, `max_steiner_multiplier`, and `strict_mesh` to `GeomMeshOpts`, `MeshRequest` (and JSON Schema), `MeshSettings`, `set_mesh`, `scriptout::mesh_args`, and `convergence::refined` (`dz_min` scales linearly when `scale_sizes` is true).
  * `sheet3d` floors vertical slab thickness at `dz_min`, grows slabs by `z_growth_ratio` (default `1.35`), enforces `max_z_slabs` while strictly preserving all stack layer interfaces, and aborts before allocation whenever `projected_tets > max_cells`.
  * Replaced single-level fill lattices and linear facet scans with a hierarchical `OctaveFillGrid` and `SegBvh2`/`Bvh2` (`< 250 ms` at 1000:1 size ratios).
  * Captured `spade` Ruppert `RefinementResult` and sliver audits into `MeshedGeom::warnings` / `MeshReport::warnings`, honoring `max_steiner_multiplier` and failing cleanly when `strict_mesh` is true.
  * Preserved interior `Op::Sheet` facets (`la == lb` and `la != lb`) through `structured_pml_3d` via `interface_facets_with_preserved` in both `z_coarsen` modes.
  * Projected horizontal curved `Circle`/`Ellipse` rim edges of `sheet3d` meshes onto exact mid-edge nodes (`curve_sheet_mesh_edges` / `apply_edge_nodes_validated`) when `curved == true`, keeping `curved == false` bit-identical.

## Consequences

* Large meander/chip layouts arrange and mesh in sub-second time with deterministic thread scaling, hard pre-allocation memory/cell bounds, explicit Ruppert quality diagnostics, PML-safe airbridges, and curved 3D sheet pads.

## Code

`crates/geom/src/{bvh,lib,complex,arrange,stack,poly_batch,select,scene,model,planar}.rs`, `crates/femcore/src/meshing/{geomdrive,request,sheet3d,mesher,adapt3d,pml_geom,pml_geom3,curve3,geomdrive3,mod}.rs`, `crates/nyx-engine/src/{project,meshing,commands,scriptout,convergence}.rs`, `crates/app/src/app.rs`, `crates/nyx-engine/tests/squadds_fr_e2e.rs`, `docs/architecture/contracts/mesh_request.schema.v1.json`.
