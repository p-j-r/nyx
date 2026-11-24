# 127. Perturbative Mattis–Bardeen IBC, hybrid quasi-electrostatic modal TLS, and memory-scalable `ems::driven` backends

Date: 2026-11-24

## Status

Accepted. Amends ADR 26, 33, 82, 86, 90, 104, 105, 106, 117, 118, and 119.

## Context

Full-chip 3D superconducting cQED workflows hit three remaining physics and memory-scalability limits (`docs/SQUADDS_BUG_FEATURE_REQUESTS.md`, Limitations L1, L3, and L4):
1. **Ill-conditioned direct thin-film IBC on 3D `cavity` (`Limitation L1` / `FR-30`):** Adding $-i\omega\mu_0 Y_s(\omega_0)$ (`|Y_s| \sim 10^3\text{–}10^5\,\mathrm{S}`) directly into the 3D Nédélec curl–curl operator breaks real symmetry, forces complex non-Hermitian Krylov Arnoldi, and creates $10^7\text{–}10^9$ tangential penalty contrast that stalls shift-invert convergence on large multi-sheet chips.
2. **Coarse-mesh 3D `cavity` TLS & `quarter_point` Jacobian singularity (`Limitation L3` / `FR-32`):** Resolving $r^{-1/2}$ edge field singularities across $2\,\mathrm{nm}$ native oxides requires local sub-micron / `quarter_point` $P_2$ refinement that is fast in scalar `dfs::electrostatic` ($1\,\text{DOF/node}$) but prohibitive in full-chip vector `ems::cavity`, while `CurvedTet::jacobian` had $\det J = 0$ right on Barsoum $\xi_{\text{QP}} = 0.25$ rim vertices.
3. **High-DOF `ems::driven` COLAMD LU memory ceiling & coarse calibration keys (`Limitation L4` / `FR-33`):** `ems::driven` previously factored complex symmetric Maxwell shifts $(K - k_0^2 M + i k_0 D) \mathbf{x} = \mathbf{b}$ with dense-supernode unsymmetric `LuColamd` and keyed engine memory admission solely by `"ems"`, causing legacy COLAMD runs to trigger false `REFUSE_ABOVE` blocks even when switching to low-memory symmetric, BLR, or iterative solvers.

## Decision

* **Perturbative Slater / Wheeler Mattis–Bardeen IBC on `ems::cavity` (`FR-30`, `nyx-ems::{problem,ibc,cavity,lib}`):**
  * `surface_impedance` gains `ibc_mode: "direct" | "perturbative" | "auto"` (default `"direct"` for exact backward compatibility).
  * Under `"perturbative"` (and `"auto"` on `cavity` when $k_0 |L_{k,\square}| / Z_0 < 0.05$), `ems::cavity` treats superconducting sheets as lossless PEC during the real symmetric shift-invert eigensolve, then evaluates the surface magnetic tangential integral $\int_{\Gamma_{\text{sc}}} |\mathbf{n} \times \mathbf{H}_0|^2\,dA$ to apply Wheeler's kinetic-inductance downward frequency shift $\Delta\omega_{\text{re}}$ and Mattis–Bardeen surface-resistance quality factor $Q_{\text{BCS}}$, updating `u_kin_ibc_j` and `CircuitRom` consistently.
* **Curved `quarter_point` cavity guard & hybrid quasi-electrostatic modal TLS (`FR-32`, `nyx-ems::cavity`, `nyx-analysis::tls`, `nyx-engine::{analyse,commands}`):**
  * `cavity::safe_curved_jacobian` pulls degenerate Barsoum $\xi_{\text{QP}} = 0.25$ vertex evaluations slightly into the cell interior (`EPS = 1e-3`) and upgrades curved quadrature to `quad_rule_collapsed(8)` whenever `mesh.has_qp()` is true.
  * `TlsParticipationAnalysis` emits `"edge_singularity_underresolved"` and an actionable `"warning"` when run directly on coarse straight-$P_1$ 3D sheet meshes (`h_edge_mean > 5e-6` without `quarter_point`).
  * `nyx_analysis::tls_participation_quasi_electrostatic` superposes fine-mesh electrostatic terminal basis fields $\mathbf{E}_m(\mathbf{x}) = \sum_i V_{m,i}\,\mathbf{E}^{(i)}(\mathbf{x})$ weighted by modal port voltages to extract converged `p_ma`, `p_ms`, `p_sa`, `q_tls`, and `T1_TLS_s`. `nyx-engine::analyse` resolves each physics from its own `Stored.model` (`ADR 90`) and exposes `electrostatic_physics`, `port_terminals`, and `modal_voltages` / `modal_freqs_hz` while preserving `ADR 26` crate layering.
* **Memory-scalable `ems::driven` solver backends & `(solver, order)` memory calibration (`FR-33`, `femcore::linalg::{lblt,blr,ams,maxwell}`, `nyx-ems::driven`, `nyx-engine::estimate`):**
  * `femcore::linalg::lblt` adds nested-dissection ordering (`LbltNd` when $n \ge 5{,}000$) alongside `LbltAmd`; `ShiftedBlrLu::replace_matrices` reuses `BlrSymbolic` elimination trees across frequency shifts (`symbolic_count == 1`); and `AmsFgmresSolver` supports complex shifted curl–curl updates (`update_matrix` + `solve_many`) for both $P_1$ and $P_2$ Nédélec spaces with floating interior PEC island DSU potentials.
  * `ems::driven` exposes `solver: "auto" | "direct" | "lblt" | "blr" | "ams_fgmres"` and `blr_tol` (default `1e-6`), logging `FactorBackend` (`LbltNd`, `LbltAmd`, `LuColamd`, `BlrNd`, `AmsFgmres`) in `SolveResult::log`.
  * `nyx-engine::estimate::chain_key_for_project` keys `ems:driven` calibration history by `"ems:driven:{solver}:p{order}"` while keeping `"ems"` unchanged for `cavity` / `modes` / `bloch`.

## Consequences

* 3D superconducting cavities with Mattis–Bardeen `surface_impedance` solve at real-PEC speed via `ibc_mode: "perturbative"`, full-chip $P_1$ `CircuitRom` modes pair cleanly with local $P_2$ `quarter_point` electrostatic meshes for converged `T1_TLS_s`, and large 3D `ems::driven` sweeps run within bounded RAM via `lblt`, `blr`, or `ams_fgmres` without false legacy COLAMD admission refusals.

## Code

`crates/femcore/src/linalg/{lblt,blr,ams,maxwell,driver}.rs`, `crates/nyx-ems/src/{problem,ibc,cavity,driven,lib}.rs`, `crates/nyx-ems/tests/driven_prom.rs`, `crates/nyx-analysis/src/{tls,lib}.rs`, `crates/nyx-engine/src/{estimate,analyse,commands,engine,sweep,worker,preview,convergence,validated}.rs`, `crates/nyx-engine/tests/fr25_33_e2e.rs`.
