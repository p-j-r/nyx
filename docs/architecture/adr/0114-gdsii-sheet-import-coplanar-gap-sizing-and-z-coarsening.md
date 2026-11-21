# 114. Pure-Rust binary GDSII `import_gds`, `sheet3d` PSLG curved-hole guard, coplanar `gap_cells`, and pure-tet `z_coarsen`

Date: 2026-11-21

## Status

Accepted. Amends ADR 3, 54, 89, and 102.

## Context

Superconducting quantum chip layouts are authored as hierarchical binary GDSII streams (`.gds`) with narrow ($2\text{–}20\,\mu\text{m}$) coplanar etched gaps inside centimetre-scale 3D vacuum/substrate boxes. Importing `.gds` without Python dependencies, preventing Ruppert PSLG encroachment slivers across curved concentric pocket notches, resolving coplanar gaps automatically without per-island manual `add_size` boilerplate, and avoiding uniform vertical prism column counts far from $z = z_{\text{sheet}}$ are essential for fast, reliable 3D cQED meshing.

## Decision

* **Zero-dependency binary GDSII parser & `Engine::import_gds` (`geom::gds`, `nyx-engine::gds`):**
  * Implements a pure-Rust big-endian GDSII stream reader and lossless 8-byte IBM hexadecimal float codec (`decode_ibm_f64` / `encode_ibm_f64`) inside `geom` (preserving `geom`'s `serde`-only dependency rule), flattening `SREF`/`AREF` (`STRANS` X-reflection, magnification, rotation, `COLROW`, cycle guard), `BOUNDARY`, `BOX`, and `PATH` (`PATHTYPE` 0, 1, 2, 4) into SI metres.
  * Decomposes keyhole cutlines, boolean-unions each requested `(layer, datatype)` via `geom::arrange` into CCW exterior hulls and CW interior holes, splits or merges islands deterministically, simplifies rings with cross-layer contact-vertex pinning (`simplify_tol`), renames islands from `TEXT` labels on `pin_layer`, and adds `Polygons` + `Sheet` features at height `z` behind `path_guard`.
* **PSLG curved-hole encroachment fix & symmetric coplanar `gap_cells` (`femcore::meshing::sheet3d`):**
  * Fixes cyclic-trapezoid co-circularity slivers on curved pocket holes (`repro_concentric_mesh_fail_P2_fast`) via AABB-pruned symmetric proximity clamping on facing sheet-gap rim edges, deterministic `fine_edges` ordering, and a collar sizing guard.
  * Adds optional `gap_cells: Option<f64>` (`MeshSettings`, `MeshRequest`, `set_mesh`), running a segment-to-segment distance pass across facing sheet-to-gap rims and clamping `ov_edge_h[e] <= d_gap / gap_cells` symmetrically with smooth `0.35` spatial grading.
* **Conforming pure-tetrahedral vertical coarsening (`z_coarsen: bool`):**
  * Keeps exact uncoarsened Dompierre prisms within $|z - z_{\text{sheet}}| \le 4 h_{\text{sheet}}$ ($\ge 2$ slabs per side) and progressively collapses interior 2D edges on outer $z$-slabs via maximal independent-set half-edge collapses (`coarsen_2d_once`: 2D Link Condition, triangle quality $\ge 0.12$, frozen lateral/cap boundaries, positive unit-slab Jacobian check) bridged by 100% conforming crack-free 3-tet/2-tet pure-tetrahedral transition slabs ($\sim 3\times$ 3D cell and P2 DOF reduction with $< 0.2\%$ capacitance shift).

## Consequences

* Foundry `.gds` transmon and resonator layouts import directly into 3D `nyx` projects in pure Rust and mesh into watertight conformal tetrahedral grids with automatic gap resolution and $\sim 3\times$ fewer far-field unknowns.

## Code

`crates/geom/src/{gds,lib}.rs`, `crates/geom/tests/gds.rs`, `crates/femcore/src/meshing/{sheet3d,geomdrive,request}.rs`, `crates/femcore/tests/{sheet3d_qcap,sheet_mesh3d,contracts}.rs`, `crates/nyx-engine/src/{gds,project,meshing,scriptout,commands,lib}.rs`, `crates/nyx-engine/tests/{import_gds,path_guard,cqed_qcap_e2e}.rs`, `docs/architecture/contracts/mesh_request.schema.v1.json`.
