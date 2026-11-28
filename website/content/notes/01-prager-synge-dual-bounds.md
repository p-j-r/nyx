+++
title = "Prager–Synge Dual Energy Bounds & Equilibrated H(div) Flux Recovery"
date = 2026-11-22
description = "How nyx turns primal electrostatic solves into rigorous two-sided Maxwell capacitance intervals [C_low, C_up] via vertex-patch Raviart–Thomas equilibration, and where Barsoum singular cells and Steiner locking set the boundaries."

[extra]
note_number = "01"
crates = "femcore::rt · nyx-dfs::electro"
reading_time = "9 min"
+++

A primal finite-element solve of Laplace's equation $\nabla \cdot (\varepsilon \nabla V_h) = 0$ is variationally one-sided. By Thomson's (Dirichlet) minimum-energy principle, any conforming potential $V_h \in H^1(\Omega)$ that satisfies the prescribed terminal voltages overestimates the true electrostatic energy, so every diagonal capacitance entry $C^{\mathrm{up}}_{ii}$ computed from $\mathbf{D}_h = -\varepsilon \nabla V_h$ sits strictly above the exact value $C^{\mathrm{exact}}_{ii}$. Without an independent lower bound, a single solve cannot certify how tight that overestimate is.{% sidenote(id="sn-thomson", label="1") %}On a fixed mesh without a lower bound, two engineers can run the same problem at P1 and P2 and observe a shift without knowing whether either mesh has converged within 1% or 10% of the continuum limit.{% end %}

In `crates/femcore/src/rt.rs` (`equilibrate_electro_flux_bounds`), `nyx` pairs the primal Lagrange solve ($p \in \{1, 2\}$) with a local $H(\mathrm{div}; \Omega)$ Raviart–Thomas flux reconstruction $\tilde{\mathbf{D}}^{(j)} \in \mathrm{RT}_{p-1}$ (`RT_0` for P1; `RT_1` for P2) to produce certified Loewner-ordered capacitance brackets $C^{\mathrm{low}} \preceq C^{\mathrm{exact}} \preceq C^{\mathrm{up}}$ on every run.

## Local Vertex Half-Star Equilibration

The complementary flux $\tilde{\mathbf{D}}^{(j)}$ for terminal $j$ must be divergence-free in the dielectric bulk ($\nabla \cdot \tilde{\mathbf{D}}^{(j)} = 0$ in $\Omega$), match normal flux continuity across dielectric interfaces, and respect the exterior `FarField` Robin condition $\tilde{\mathbf{D}}^{(j)} \cdot \mathbf{n} = \alpha V_h^{(j)}$ on $\Gamma_\infty$. Rather than solving a costly global mixed problem, `rt.rs` solves small local patch problems around each mesh vertex $v$ using the Braess–Schöberl / Ern–Vohralik equilibration construction.

Crucially, zero-thickness conductive sheets (such as CPW signal traces and ground planes) introduce tangential jumps in $V_h$ and normal jumps in $\mathbf{D}$. When a vertex lies on an interior conductive sheet, `rt.rs` splits its vertex star into disconnected *half-stars* (`n_sheet_splits` in `EquilibratedFluxResult`), solving each side independently so that normal flux jumps across metal sheets are preserved exactly without polluting the global $L^2(\Omega)$ divergence norm `max_div_l2`.

{% callout(title="Raviart–Thomas Local Degrees of Freedom in femcore::rt", kind="note") %}
`rt_nloc(dim, order)` allocates `3` face DOFs (`RT_0`) or `8` DOFs (`6` face + `2` interior cell bubbles, `RT_1`) on 2D triangles, and `4` face DOFs (`RT_0`) or `15` DOFs (`12` face + `3` cell bubbles, `RT_1`) on 3D tetrahedra, mapped from reference elements via the contravariant Piola transform $\frac{1}{|\det J|} J \hat{\phi}$.
{% end %}

## Complementary Energy Matrices & Entrywise Error Bounds

Given the primal potential vectors $\{V_h^{(j)}\}_{j=1}^M$ and equilibrated $H(\mathrm{div})$ fluxes $\{\tilde{\mathbf{D}}^{(j)}\}_{j=1}^M$, `EquilibratedFluxResult` forms two $M \times M$ symmetric matrices:

$$
W^{\mathrm{flux}}_{ij} = \int_\Omega \tilde{\mathbf{D}}^{(i)} \cdot \varepsilon^{-1} \tilde{\mathbf{D}}^{(j)}\, d\Omega + \int_{\Gamma_\infty} \alpha^{-1} \left(\tilde{\mathbf{D}}^{(i)} \cdot \mathbf{n}\right)\left(\tilde{\mathbf{D}}^{(j)} \cdot \mathbf{n}\right) dS
$$

$$
Q_{ij} = -\int_\Omega \nabla V_h^{(i)} \cdot \tilde{\mathbf{D}}^{(j)}\, d\Omega + \int_{\Gamma_\infty} V_h^{(i)} \left(\tilde{\mathbf{D}}^{(j)} \cdot \mathbf{n}\right) dS
$$

From the Cauchy–Schwarz inequality on the mixed pairing $Q$ and complementary flux energy $W^{\mathrm{flux}}$, the exact Maxwell capacitance matrix $C^{\mathrm{exact}}$ is bracketed in the Loewner positive-semidefinite order:

$$
C^{\mathrm{low}} = Q \left(W^{\mathrm{flux}}\right)^+ Q^T \preceq C^{\mathrm{exact}} \preceq C^{\mathrm{up}}
$$

Taking the midpoint $C^{\mathrm{mid}} = \frac{1}{2}\left(C^{\mathrm{up}} + C^{\mathrm{low}}\right)$ provides a certified entrywise bound on every self and mutual capacitance:

$$
\left| C^{\mathrm{mid}}_{ij} - C^{\mathrm{exact}}_{ij} \right| \le \delta C_{ij} = \frac{1}{2}\sqrt{\left(C^{\mathrm{up}}_{ii} - C^{\mathrm{low}}_{ii}\right)\left(C^{\mathrm{up}}_{jj} - C^{\mathrm{low}}_{jj}\right)}
$$

Simultaneously, the cell-local Prager–Synge energy gap

$$
\eta_{T,\mathrm{PS},j}^2 = \left\| \varepsilon^{-1/2}\left(\tilde{\mathbf{D}}^{(j)} - \mathbf{D}_h^{(j)}\right) \right\|_{L^2(T)}^2
$$

feeds goal-oriented `flux_jump_goal` adaptive mesh refinement (`amr_ps`), driving refinement strictly into the edge singularities and narrow gaps that dominate the energy gap.{% sidenote(id="sn-epistemic", label="2") %}In ASME V&V 20 validation (`benchmarks/validation/README.md`), the half-width $\delta_{\mathrm{ep}} = \frac{1}{2}(S_{\mathrm{up}} - S_{\mathrm{low}})$ is kept strictly epistemic and never folded into Gaussian variance $u_{\mathrm{aleatoric}}^2$.{% end %}

## Why `quarter_point` Singular Cells Are Excluded

In planar superconducting circuits, the electric field diverges as $r^{-1/2}$ at the edge of a zero-thickness metal trace. `nyx-dfs` offers two distinct strategies for capturing this singularity:

1. **Goal-oriented `flux_bounds` + `amr_ps` (`AF-1 + AF-4`):** Straight-sided P2 cells adapted around the trace perimeter, paired with `equilibrate_electro_flux_bounds`.
2. **Barsoum `quarter_point` edge grading (`AF-2 + AF-3`):** Quadratic mid-edge nodes on triangles/tetrahedra touching metal rims are shifted to the quarter-point position $\xi_{\mathrm{QP}} \in [0.25, 0.50]$ so the isoparametric Jacobian produces an exact $r^{-1/2}$ gradient singularity on a coarse mesh.

Because $\det J \to 0$ at the collapsed corner of a Barsoum quarter-point cell, `eval_rt_piola` (which divides by $|\det J|$) is singular there. `crates/nyx-engine/src/solve.rs` (`L145–157`) therefore enforces a strict compile-and-runtime separation: `flux_bounds` and `prager_synge` AMR are explicitly forbidden whenever `quarter_point` grading is active, while `qp_graded` meshes are validated via interior 8-point collapsed quadrature (`scaled_jacobian_min_interior_q8`).

## Benchmark Evidence vs. Ansys Q3D (`S1.1–S1.4`)

Across the `418`-run SQuADDS cQED benchmark (`docs/SQUADDS_BENCHMARK_REPORT.md`), both strategies were evaluated on unmodified `HEAD` (`9480093a`) against Ansys Q3D:

- **Certified Bracketing (`amr_ps`):** Across `12` stratified designs (`S1.1–S1.4`), `flux_bounds="midpoint"` strictly bracketed **`68/69` (`98.6%`)** of all Q3D capacitance matrix entries (`100%` across `S1.2 NCap`, `S1.3 CapNInterdigitalTee`, and `S1.4 GenIDC`).
- **Quarter-Point Accuracy (`qp_graded`):** Across all `60/60` `S1.1–S1.3` designs, `dfs::electro` (`order = 2`, `AF-2` + `AF-3` spherical harmonic `far_field(l=2)`) achieved **`1.11%` unweighted (`1.04%` Horvitz–Thompson `[0.78%, 1.35%]`) self-capacitance median error**, `MAE = 1.369 fF`, `-0.06%` signed bias, exact $p_{\mathrm{obs}} = 2.00$ Richardson convergence on `S1.2`, and **bitwise-exact (`0.000e0 < 1e-13`) multi-threaded AMG determinism** (`3.08x–3.67x` speedup at `NYX_THREADS=8`).

{% solver_log(cmd="cargo test --release -p nyx-engine --test validation_cqed_cross_code") %}
coax_cylinder_electrostatic (a=10 um, b=30 um, L=100 um, eps_r=11.7):
  exact C_coax = 5.924750e-14 F
  Prager-Synge [C_low, C_up] brackets exact value with deterministic delta_ep
  dfs::electro NYX_THREADS=1 vs NYX_THREADS=8 max rel drift: 0.000e+00 (bitwise exact)
{% end %}

## Honest Disclosures: Where Baseline `HEAD` Breaks

Four documented defects in `docs/SQUADDS_BENCHMARK_REPORT.md` define the current operating envelope for volumetric electrostatics:

1. **Narrow-Gap PSLG Steiner Locking (Core Bug 4):** On `S1.4` (`GeneralizedCapNInterdigital`) and `Ref.0` (`QCap` GDS) geometries featuring dense `200–300 nm` fillet chords across `2–3 um` gaps inside a large ground pocket, Ruppert 2D refinement exhausts `max_steiner_multiplier = 20` near the short chords, leaving single-prism cross-gap triangles (`365.63%` volumetric error). Switching to the 3-layer Green-function solver **`planar_bem` (`QFR-4`)** eliminates volumetric locking completely and solves `S1.4` **`76x` faster (`1.03 s`) at `4.04%` self-capacitance median error**.
2. **Finite-Die `FR-11` Thickness Guard (Core Bug 1):** When `sub_box` sits inside a larger vacuum box, `crates/nyx-analysis/src/tls.rs` (`L789–797`) mistakes the four outer vertical `Si|Vac` die sidewalls (`|n_z| < 0.99`) for etched trenches and aborts `metal_thickness_m` edge corrections.
3. **Multi-Layer `z_coarsen` Collapse (Core Bug 3):** Enabling `z_coarsen: true` in `sheet3d.rs` collapses exterior cap boundary IDs on finite-die stacks; keep `z_coarsen: false` on multi-layer runs.
4. **`Mode G2` Dual-Box `1/R` Extrapolation Non-M-Matrix:** Independent meshes at outer box radii $R_1$ and $R_2$ caused Richardson $C_\infty = (R_2 C_2 - R_1 C_1)/(R_2 - R_1)$ to produce a small positive off-diagonal entry (`max_offdiag_pos_ratio = 3.75e-3`) on `1/12` ablation runs, whereas single-solve `AF-3` `far_field(l=2)` had `0/249` M-matrix violations.
