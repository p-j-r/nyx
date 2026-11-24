# 125. Finite-die perimeter exclusion, per-terminal `rank1` lumped ports, and deterministic `ems` assembly/solvers

Date: 2026-11-24

## Status

Accepted. Amends ADR 11, 33, 103, 104, 105, 106, 117, and 118.

## Context

Post-SQuADDS benchmark testing exposed three numerical and topological defects on realistic planar cQED chips (`docs/SQUADDS_BUG_FEATURE_REQUESTS.md`, Bugs 1, 2, and 5):
1. **Finite-die boundary segment exclusion (`Bug 1` / `FR-25`):** When conductive sheets (`ground_plane`, CPW launchpads) extend flush to the lateral `xmin`/`xmax`/`ymin`/`ymax` bounding box of a finite die, `detect_sheet_perimeter_3d` classified those outer-die segments as singular interior metal–dielectric knife-edges. For `far_field` / `open` exterior BCs, `E` is nonzero on the outer die boundary, blowing up `I_edge` by $10^2\times$ and breaking `metal_thickness_m` capacitance correction (`non-SPD corrected Maxwell matrix`).
2. **Co-planar ground-split `rank1` ports (`Bug 2` / `FR-26`):** `rank1` `lumped_jj` and `lumped_port` sheets built $\mathbf{g}_{\text{port}}$ from a single shortest edge path between `T_-` and `T_+`. When a sheet bridges a split ground plane (two disconnected PEC islands on `T_+` at the same potential, e.g. CPW left/right ground straps around a resonator or SQUID), a single path shorted only one side and left an unphysical slot-line circulation loop on the other.
3. **Thread-count drift in 3D `ems` (`Bug 5` / `FR-29`):** While `dfs` and `sms` achieved `0.00e0` thread drift across `NYX_THREADS = 1..8`, `ems` `cavity` and `driven` drifted by $\sim 5\times 10^{-11}$ relative frequency because parallel matrix assembly partitioned cells by `threads * 4` (`div_ceil`) and `femcore::linalg` supernodal `LuColamd` / `Lblt` / dense `faer` factorizations switched parallel elimination trees with thread count.

## Decision

* **Outer-bbox perimeter segment exclusion (`FR-25`, `nyx-analysis::{tls,capacitance}`):**
  * `detect_sheet_perimeter_3d` identifies mesh edges lying on the 6 outer bounding-box planes (`ext_edges` via `is_on_outer_bbox`) and excludes `ext_edges` when extracting the 1D singular perimeter `SegmentBvh3` of zero-thickness sheets. True interior etched edges (pocket gaps, CPW slots, interdigital fingers) remain intact while truncated outer die borders no longer masquerade as singular knife-edges.
  * `capacitance::extract_terminal_basis_edge_intensities` adds an explicit physical clamping guard whenever a cylindrical `R_0` tube cannot fit cleanly inside the mesh, preventing non-SPD over-subtraction in `apply_3d_metal_thickness_correction`.
* **Per-terminal multi-component `rank1` port basis (`FR-26`, `nyx-ems::{ibc,ports}`):**
  * `SheetPortGeometry` (`ibc.rs`) and `DrivenPortGeom` (`ports.rs`) run a global multi-source Dijkstra distance field from all `T_-` vertices across the port sheet to all `T_+` vertices, group reached `T_+` vertices by PEC-island connected component, extract one minimum-hop shortest path per distinct `T_+` component, and average the unit-voltage path vectors ($\mathbf{g}_{\text{port}} = \frac{1}{K}\sum_{c=1}^K \mathbf{g}^{(c)}$). Single-island `T_+` geometries ($K = 1$) produce bit-identical results; split co-planar grounds ($K \ge 2$) enforce equal voltage across every ground bank.
* **Thread-count-invariant `ems` assembly and `faer` / supernodal kernels (`FR-29`, `nyx-ems::{cavity,bloch,modes,modes_complex,axisym}`, `femcore::linalg::{faer_glue,lblt,sparse_lu}`):**
  * All `ems` 2D/3D assemblers partition elements into fixed `FIXED_CELL_CHUNK = 2048` ranges independent of `femcore::par::threads()` and merge partial CSR matrices through a deterministic pairwise binary reduction tree (`reduce_tree_ccsr`).
  * `lblt.rs` and `sparse_lu.rs` execute supernodal `Par::Seq` factorizations with deterministic sibling postorder merges, and `faer_glue.rs` runs dense factorizations/EVDs sequentially up to `SMALL_DENSE_SEQ_N = 1024`, eliminating `NYX_THREADS = 1` vs `8` eigenvalue and S-parameter drift (`< 10^{-14}`).

## Consequences

* `metal_thickness_m` 3D capacitance corrections preserve Loewner SPD positivity on finite dies with flush ground planes and `open`/`far_field` BCs, `rank1` and `distributed_sheet` lumped ports agree on split-ground CPW/SQUID layouts, and `ems` `cavity`/`driven` runs are thread-count invariant within `< 10^{-13}` across `NYX_THREADS = 1..8`.

## Code

`crates/nyx-analysis/src/{tls,capacitance}.rs`, `crates/nyx-analysis/tests/squadds_smoke.rs`, `crates/nyx-ems/src/{ibc,ports,cavity,bloch,modes,modes_complex,axisym}.rs`, `crates/femcore/src/linalg/{faer_glue,lblt,sparse_lu}.rs`, `crates/nyx-engine/tests/squadds_engine_smoke.rs`.
