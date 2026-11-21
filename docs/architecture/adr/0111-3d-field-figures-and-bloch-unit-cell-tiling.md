# 111. Un-averaged P2 3D field figures (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`) and zero-copy Bloch tiling

Date: 2026-11-21

## Status

Accepted. Amends ADR 99.

## Context

`figure` previously drew only 2D sections (`gallery` and `components`, cutting a 3D mesh at one `y` plane). Inspecting 3D photonic gratings, cavities, and superconducting chips requires arbitrary planar slices, box-clipped cutaways, translucent 3D isosurfaces, mesh inspection views, and multi-period Bloch wave visualization—while preserving sharp normal electric-field jumps ($\varepsilon_1 E_{n,1} = \varepsilon_2 E_{n,2}$) across high-contrast dielectric interfaces without vertex-averaging smear.

## Decision

* **Un-averaged 8-subtet P2 evaluation (`nyx-analysis::view3d`):**
  * Each quadratic (10-node) tetrahedron is subdivided in nyx's exact barycentric node order (`v0..v3, e01..e23`) into 8 positively oriented linear sub-tetrahedra (`P2_SUBTETS`, each $+1/8$ parent volume) evaluated directly on per-cell `CellFields` without cross-cell vertex welding, preserving sharp material-interface field jumps pixel-for-pixel.
  * Arbitrary scalar formulas (`expr`, e.g. `"re(Ex)"`, `"abs2(E)"`, `"im(Hz)"`) are parsed once via `femcore::formula` and evaluated over complex vector/scalar cell fields; signed expressions default to symmetric zero-centered `Coolwarm` colormaps.
* **3D figure kinds and Z-buffered rasterization (`nyx-analysis::{figure,view3d,raster}`):**
  * Adds `Kind::{Slice3d, Cutaway3d, Isosurface3d, Mesh3d}` plus `plane`, `clip`, `isosurfaces`, `show_edges` (`none`, `regions`, `mesh`), and orthographic `camera` (`azimuth_deg`, `elevation_deg`).
  * `raster::raster_field_3d` renders opaque triangles with a per-pixel `f32` Z-buffer, depth-biased wireframe edge quads, and back-to-front depth-sorted translucent triangles with Porter–Duff `over` compositing across L2-cached horizontal strips (`femcore::par::map_chunks`), producing bit-identical PNG output regardless of `NYX_THREADS`.
* **Zero-copy Bloch unit-cell tiling (`repeat_periodic`):**
  * Replicates geometry along `m * periodic_shift` for $m \in \{0, \dots, \text{repeat\_periodic}-1\}$ and multiplies complex cell fields by $\exp(-i m\,\text{periodic\_phase})$ *before* evaluating `expr`, with `Engine::figure` auto-resolving `periodic_shift` and `periodic_phase` from `periodic_pairs()`, result scalars (`phase N (rad)` / `k_y re N (rad/m)`), or periodic boundary conditions.

## Consequences

* Scripts, the desktop app, and MCP sessions can render publication-grade 3D slices, cutaways, isosurfaces, and multi-period propagating Bloch waves directly to `.png` or `.svg` in tens of milliseconds without external post-processors.

## Code

`crates/nyx-analysis/src/{view3d,figure,raster,lib}.rs`, `crates/nyx-analysis/tests/figure3d.rs`, `crates/nyx-engine/src/{figure,commands}.rs`, `crates/nyx-engine/tests/{figure3d,periodic_photonics_e2e}.rs`.
