# 126. `z_coarsen` cell-budget enforcement, medial-axis gap seeding, and Stieltjes–Padé capacitance extrapolation

Date: 2026-11-24

## Status

Accepted. Amends ADR 6, 102, 103, 112, 114, 116, 120, and 121.

## Context

Three meshing and electrostatic convergence bottlenecks appeared on dense superconducting multi-element chips (`docs/SQUADDS_BUG_FEATURE_REQUESTS.md`, Bugs 3–4 and Limitation L2):
1. **Unchecked `z_coarsen` budget overflow (`Bug 3` / `FR-27`):** `sheet3d` checked `max_cells` only before building conforming $z$-coarsening transition slabs; when initial fine vertical layers near $z = 0$ already approached or exceeded `max_cells`, `mesh_sheet_stack_3d` returned an oversized mesh instead of failing fast or coarsening earlier.
2. **Single-element cross-gap bridging (`Bug 4` / `FR-28`):** Even with `gap_cells >= 2.0`, `gap_sizing_seeds_2d` clamped vertex target sizes $h(\mathbf{v})$ only *on* conductor boundaries; Ruppert circumradius refinement did not split triangles spanning a narrow dielectric gap $g \ll h_{\text{global}}$ with two vertices on one conductor and one on the opposite conductor, creating direct Dirichlet-to-Dirichlet cross-gap edges that inflated mutual capacitance by $20\text{–}40\%$.
3. **Log-singular 3D knife-edge extrapolation (`Limitation L2` / `FR-31`):** Zero-thickness 3D planar conductors exhibit inverse-logarithmic / fractional-power pre-asymptotic errors $C(h) = C_\infty + a_1 h^{p_1} + a_2 h^{p_2}$ ($p_1 \approx 0.5\text{–}0.85$), causing fixed-order polynomial Richardson extrapolation (`p = 1` or `2`) on a 3-pass AMR sequence to overshoot or stagnate at $1\text{–}2\%$.

## Decision

* **Exact `max_cells` enforcement & adaptive `z_coarsen` stride (`FR-27`, `femcore::meshing::sheet3d`):**
  * When `z_coarsen: true` and `max_cells = Some(limit)` are both active, `extrude_layer_z_coarsened` advances transition coarsening immediately after the first adjacent interface slab (`n_fine_keep = 1`) and adapts coarsening stride when the baseline column count threatens `limit`.
  * A post-assembly guard verifies `total_tets <= limit` across all `sheet3d` paths (uniform and `z_coarsen`), returning a descriptive error naming both the actual cell count and `max_cells` whenever even the minimum 1-slab-per-layer topological floor cannot fit inside `limit`.
* **Medial-axis gap centerline Steiner seeding & cross-gap edge splitting (`FR-28`, `femcore::meshing::{sheet3d,mesher}`):**
  * When `gap_cells = Some(n_g)` with $n_g \ge 2.0$, `gap_sizing_seeds_2d` emits interior medial-axis Steiner seed points (`Complex::seed_points` in `Model2`) at midpoint $\mathbf{m} = \frac{1}{2}(\mathbf{p}_a + \mathbf{p}_b)$ with target size $h_{\text{gap}} = d / n_g$ along detected close-approach conductor pairs ($d < n_g \cdot h_{\text{global}}$).
  * During Ruppert refinement (`mesher.rs`), any interior dielectric triangle edge directly bridging two distinct conductor sheets whose length exceeds $1.5\,h_{\text{gap}}$ is split at its midpoint before termination, guaranteeing $\ge 2$ elements (and $\ge 1$ free interior dielectric vertex) across resolved coplanar gaps.
* **3-level Aitken / Stieltjes–Padé SPD capacitance extrapolation (`FR-31`, `nyx-analysis::capacitance`):**
  * When an electrostatic AMR run records $\ge 3$ passes (`amr.history.len() >= 3`) with `0` curved elements, `capacitance_matrix` fits the effective asymptotic convergence exponent $p_{\text{eff}} \in [0.35, 2.5]$ from the last three `(n_dofs, C)` levels via $h_{\text{eff}} = N_{\text{dof}}^{-1/d}$ and applies a monotone-guarded SPD Loewner Aitken $\Delta^2$ / $[1/1]$ Stieltjes–Padé extrapolation (`extrapolated_maxwell_matrix_F` and `effective_convergence_order`), falling back gracefully when non-monotone or outside Prager–Synge Loewner bounds.

## Consequences

* `set_mesh(max_cells=...)` is a strict upper bound under both uniform extrusion and `z_coarsen: true`, narrow CPW and transmon pocket gaps never suffer single-edge Dirichlet shorting when `gap_cells >= 2`, and 3-pass straight-element 3D AMR capacitance matrices extrapolate within $< 0.5\%$ of ultra-fine reference values without requiring curved elements.

## Code

`crates/femcore/src/meshing/{sheet3d,mesher}.rs`, `crates/nyx-analysis/src/capacitance.rs`, `crates/nyx-engine/tests/squadds_engine_smoke.rs`.
