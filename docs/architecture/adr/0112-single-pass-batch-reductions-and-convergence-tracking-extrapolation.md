# 112. Single-pass bit-exact batch reductions, mode-signature convergence tracking, and SPD Richardson extrapolation

Date: 2026-11-21

## Status

Accepted. Amends ADR 61, 80, 90, 91, and 97.

## Context

Computing energy confinement or overlap-like integrals across 12 modes, 4 regions, and multiple expressions previously required dozens of separate `integrate` / `confinement` round-trips, each re-decoding materials and re-walking the mesh. Meanwhile, `convergence_check` only refined global `h` (`0 < factor < 1`), left local `add_size` rules untouched, could lock onto cladding/PML distractor modes when tracking by frequency alone, and did not perform two-grid Richardson extrapolation on slowly converging $O(h)$ 3D sheet capacitances.

## Decision

* **Single-pass bit-exact batch `integrate` and `confinement` (`nyx-engine::reduce`):**
  * `integrate` accepts `solutions` (`"all"` or `[usize]`), `by_region: true`, and batch `expr` lists or named maps (`BatchExpr::{Single, List, Map}`, distinguished cleanly from single AST node objects via `is_ast_object`); `confinement` accepts `solutions`, `by_region: true`, and `regions: [str]`.
  * Shares non-`Field` leaves (`Coord`, `Material`, `Region`, `Normal`) across solutions in one `Prep`, compiles expressions into a common-subexpression-deduplicated `FlatDag`, and accumulates all `(solution, region, expr)` totals in a single `BLOCK = 2048` pass with exact `to_bits()` floating-point parity against standalone calls across all `NYX_THREADS`.
* **Bidirectional `factor`, proportional `scale_sizes`, and mode-signature `track` (`nyx-engine::convergence`):**
  * `convergence_check` accepts any finite `factor > 0` with `factor != 1` (`< 1` refines, `> 1` coarsens, consistently ordering `(coarse, fine)`).
  * `scale_sizes` (default `true`) scales non-suppressed `Op::Size` (`h`, `max_deviation` by `* factor`; `min_layers` by `/ factor`) and `MeshSettings::gap_cells` (`/ factor`) on the scratch copy without mutating the active `Project`.
  * `TrackFilter` extends `track` with `f` / `n_b` / `n_eff` targets plus `region` + `min_share`, `component` + `min_component`, and `max_layer_fraction`, filtering candidate modes on *both* coarse and fine meshes before nearest-target or mutual-nearest pairing.
* **Log-space SPD Maxwell capacitance Richardson extrapolation (`extrapolate`, `apply_extrapolated`):**
  * Supports `"richardson_log_spd"` (extrapolating off-diagonal mutual couplings $|C_{ij}|$ and diagonal self-excesses $C_{ii} - \sum_{j\ne i}|C_{ij}|$ in log-space with weights $w_1 = r^p/(r^p - 1), w_0 = w_1 - 1$, default $p = 1$ on 3D/sheet capacitance meshes, guaranteeing symmetric positive-definite Maxwell sign structure even under coarse overshoot) and `"richardson_linear"`.
  * When `apply_extrapolated: true`, `Engine::apply` promotes the extrapolated `CapacitanceMatrixData` (and `C_*` / `E_C_*` scalars) onto the stored latest `SolveResult` so downstream `analyse("capacitance_matrix", ...)` consumes the extrapolated matrix directly.

## Consequences

* Multi-mode, multi-region reductions run in a single cache-friendly pass (< 10 ms for 12 modes $\times$ 4 regions $\times$ 4 expressions) with zero rounding drift from standalone calls.
* Mesh convergence checks refine local waveguide/pad features automatically, reject cladding distractors, and promote $O(h^2)$-accurate Richardson-extrapolated capacitance matrices into the cQED analysis pipeline.

## Code

`crates/nyx-engine/src/{reduce,convergence,validated,commands,engine,commandlist}.rs`, `crates/nyx-engine/tests/{reduce,convergence_sizes,convergence_extrapolate,periodic_photonics_e2e,cqed_qcap_e2e}.rs`.
