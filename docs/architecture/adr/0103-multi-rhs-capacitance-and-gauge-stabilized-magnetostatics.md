# 103. Multi-RHS capacitance matrices and gauge-stabilized magnetostatics

Date: 2026-11-21

## Status

Accepted (superconducting qubit milestone M2). Extends the diffusion (`dfs`) physics of ADR 41 and the post-solve analyses of ADR 56.

## Context

Lumped-element quantization (`scqubits`, LOM, PyEPR pre-sizing) requires two static electromagnetic circuit matrices: (1) the $N \times N$ Maxwell capacitance matrix $\mathbf{C}_M$, mutual capacitance matrix $\mathbf{C}_{\text{mut}}$, and charging-energy matrix $\mathbf{E}_C / h = \frac{e^2}{2h}\mathbf{C}_M^{-1}$ across floating islands and couplers, and (2) multi-terminal geometric + superconducting kinetic inductance matrices $\mathbf{L}$ across 2D cross-sections (per-unit-length $L'_{ij}$), planar thin-film structures parameterized by Pearl length $\Lambda = 2\lambda_L^2/d$, and 3D Nedelec geometries. Solving terminals sequentially repeats expensive sparse matrix factorizations, and 3D $H(\text{curl})$ magnetostatics suffers from the zero-eigenvalue gradient nullspace unless gauge-regularized.

## Decision

* **Single-factorization multi-terminal electrostatics (`crates/nyx-dfs/src/electro.rs`).** Adds `BcKind::Terminal(name)` to `dfs` `electrostatic` (2D and 3D, P1/P2 including curved P2 3D tetrahedra). For $N$ named terminals, `solve_electrostatics` assembles the free-free stiffness system once, factorizes it **once** via a single `DirectSolver` LU/Cholesky factorization (or iterative CG reuse when requested), computes $N$ unit-potential basis solutions $\mathbf{V}^{(i)}$ (`Terminal i` $= 1\ \text{V}$, other terminals and `Ground` $= 0\ \text{V}$), evaluates the quadratic energy inner product $C_{ij} = \varepsilon_0 (\mathbf{V}^{(i)})^T \mathbf{K}_\varepsilon \mathbf{V}^{(j)}$, symmetrizes and extracts $\mathbf{C}_M$, $\mathbf{C}_{\text{mut}}$, and $\mathbf{E}_C/h$ into `CapacitanceMatrixData`, and stores one named `Solution` per terminal (`solution = 0..N-1`) plus typed payload in `SolveResult::scalars`.
* **Gauge-stabilized London magnetostatics (`crates/nyx-dfs/src/magneto.rs`).** Adds `Study::Magnetostatic` to `dfs` with `formulation` selecting:
  - `cross_section` (2D): out-of-plane vector potential $A_z \in H^1(\Omega_{2\text{D}})$ with $(\mu_0 \mu_r)^{-1}\nabla A_z \cdot \nabla v$ plus London kinetic sheet admittance $L_{k,\square}^{-1} A_z v$, yielding per-unit-length $L'_{ij}$ (H/m);
  - `planar_stream` (2D): Brandt–Clem scalar stream function $g \in H_0^1(\Gamma_{\text{SC}})$ ($\mathbf{K}_s = \nabla\times(g\hat{\mathbf{z}})$) combining kinetic $\int L_{k,\square} \nabla g \cdot \nabla \psi\, dA$ with the exact analytical polar self-term Biot–Savart geometric kernel, yielding $L_{ij}$ (H);
  - `mixed_3d` (3D): coupled lowest-order Nedelec $\mathbf{A} \in H(\text{curl}, \Omega)$ and 2D surface potential $\phi_s \in H^1(\Gamma_{\text{SC}})$ on single-facet deduplicated `london_sheet` boundaries with gradient gauge stabilization $\epsilon_{\text{reg}} = 10^{-6}/L_{\text{domain}}^2$, pinning one reference vertex per connected sheet island and supporting `current_terminal` and `flux_loop` cohomology excitations.
* **Post-solve analyses (`crates/nyx-analysis/src/capacitance.rs`).** Implements `CapacitanceMatrixAnalysis` (`id = "capacitance_matrix"`) and `InductanceMatrixAnalysis` (`id = "inductance_matrix"`) registered in `nyx-engine` to expose both typed JSON reports through `nyx.analyse(...)`.

## Consequences

* A 3D coaxial capacitor matches analytical $2\pi\varepsilon_0\varepsilon_r L / \ln(b/a)$ to `< 0.003%` using one factorization, and 3D planar sheet transmon+coupler capacitance matrices converge to `< 1%`.
* `physics_dfs.v1.json` and `physics_dfs_result.v1.json` gain `terminal`, `magnetostatic`, `london_sheet`, `current_terminal`, and `flux_loop` descriptors and result shapes (additive).

## Code

`crates/nyx-physics/src/cqed.rs` (`CapacitanceMatrixData`, `InductanceMatrixData`), `crates/nyx-dfs/src/{electro,magneto,lib,problem,heat}.rs`, `crates/nyx-analysis/src/capacitance.rs`, `crates/nyx-engine/src/analyse.rs`; test `crates/nyx-dfs/tests/capacitance_qubit.rs`.
