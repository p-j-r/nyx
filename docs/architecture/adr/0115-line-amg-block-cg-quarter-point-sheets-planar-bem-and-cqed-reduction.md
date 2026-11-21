# 115. Vertical-line P1-iso-P2 Block-CG AMG, quarter-point sheet edges, layered `planar_bem`, and cQED Schur/qubit reduction

Date: 2026-11-21

## Status

Accepted. Amends ADR 8, 11, 82, and 103.

## Context

Extracting multi-terminal Maxwell capacitance matrices to $< 1\%$ on 3D thin-sheet superconducting circuits faces three numerical bottlenecks: (1) high vertical aspect ratios ($100:1$ in thin vacuum/oxide layers) stall point-smoothed AMG while sequential per-terminal PCG traverses the matrix $K$ times; (2) the $r^{-1/2}$ surface-charge singularity at open zero-thickness conductor rims limits standard P2 FEM convergence to $O(h)$; and (3) users need both a sub-second meshless/panel 2D layered-dielectric cross-check and direct analytical reduction of the $N\times N$ Maxwell matrix into grounded/floating transmon $E_C$, $f_{01}$, $\beta$, $g$, and $\chi$.

## Decision

* **Vertical-line P1-iso-P2 AMG smoother & bit-deterministic multi-RHS Block-CG (`femcore::{fem,linalg}`):**
  * `detect_vertical_columns` groups P1/P2 nodes sharing quantized $(x, y)$ coordinates into sorted-$z$ columns; when columns cover $> 50\%$ of vertices, `Amg` (`amg_line_smooth`, default `true`, `AMG_THRESHOLD_3D = 25_000`) runs a pre-factored $LDL^T$ Thomas tridiagonal vertical-line smoother alongside Chebyshev point smoothing.
  * `block_pcg` (`AmgSystem::solve_block_with`) solves all $K$ terminal RHS columns simultaneously in one `[n x K]` block traversal per V-cycle using fixed-chunk (`BLOCK_ROW_CHUNK = 512`) deterministic reductions, paired with thread-independent `assemble_parallel` cell partitioning and `Threads::par_for_factor`, achieving `0.0` bit drift across `NYX_THREADS`.
* **Barsoum quarter-point P2 rim enrichment (`singular_edges = "quarter_point"` in `nyx-dfs::electro`):**
  * Identifies open perimeter vertices of interior conductor sheets on 3D P2 meshes and shifts mid-edge nodes on edges with one rim endpoint (`is_rim[a] ^ is_rim[b]`) to $0.75\,\mathbf{p}_{\text{rim}} + 0.25\,\mathbf{p}_{\text{other}}$, validating positive `CurvedTet` Jacobians across all collapsed quadrature points to capture the $r^{1/2}$ potential / $r^{-1/2}$ charge edge singularity.
* **2D layered-dielectric boundary-element solver (`formulation = "planar_bem"` in `nyx-dfs::planar_bem`):**
  * Evaluates the exact static image-series Green's function $G(\rho; t_{\text{sub}}, \varepsilon_{\text{top}}, \varepsilon_{\text{sub}}, \varepsilon_{\text{bot}})$ (including geometric tail summation and $k_{23} = -1$ grounded backplane), assembles the elastance matrix $P_{mn}$ with Wilton analytical polar integration on self/near triangle panels and 4-point Gauss quadrature on far panels, solves $C = K^T P^{-1} K$ via `dense_eigh`, and supports $1\to 4$ panel subdivision with log-space SPD Richardson extrapolation (`richardson = true`).
* **cQED Schur elimination and lumped qubit/resonator reduction (`nyx-analysis::capacitance`):**
  * `capacitance_matrix_with_params` (`Analysis::run_with_params`) preserves the exact 4-field legacy JSON when `params` is empty/null, and adds `ground` elimination (`"drop"` or `"floating_schur"` $C_{ss} - C_{sg} C_{gg}^{-1} C_{gs}$ with null-mode detection), 1-pad `Xmon` and 2-pad differential `TransmonPocket` reduction ($C_\Sigma$, $E_C/h$, $C_{\text{eff}}$, $\beta$, $E_J/E_C$, $f_{01}$, anharmonicity, $n_{\text{zpf}}$), and quarter-wave/half-wave/lumped resonator $V_{\text{rms}}$, $g/2\pi$, and dispersive shift $\chi/2\pi$.

## Consequences

* 3D multi-terminal superconducting capacitance solves achieve machine-precision agreement between Block-CG line-AMG and direct Cholesky with exact thread-count bit identity, $4\text{–}5\times$ lower single-grid rim error with quarter-point P2 elements, $< 1\%$ cross-verification against `planar_bem`, and one-call extraction of circuit Hamiltonian parameters.

## Code

`crates/femcore/src/fem.rs`, `crates/femcore/src/linalg/amg_elastic/mod.rs`, `crates/femcore/src/linalg/direct.rs`, `crates/nyx-dfs/src/{planar_bem,electro,problem,lib}.rs`, `crates/nyx-dfs/tests/qcap_solvers.rs`, `crates/nyx-analysis/src/capacitance.rs`, `crates/nyx-analysis/tests/cqed_capacitance_reduction.rs`, `crates/nyx-engine/src/commands.rs`, `crates/nyx-engine/tests/cqed_qcap_e2e.rs`, `docs/architecture/contracts/physics_dfs.v1.json`.
