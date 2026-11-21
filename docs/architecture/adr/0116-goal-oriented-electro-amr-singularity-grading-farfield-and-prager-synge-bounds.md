# 116. Goal-oriented electrostatic AMR, quarter-point & radical corner singularity grading, asymptotic `far_field` BCs, and Prager–Synge $H(\mathrm{div})$ capacitance bounds

Date: 2026-11-22

## Status

Accepted. Amends ADR 4, 7, 51, 54, 61, 82, 102, 103, 114, and 115.

## Context

Achieving sub-$0.1\%$ capacitance accuracy on 3D superconducting qubit chips without manual `add_size` tuning faces four coupled error sources: (1) uniform or hand-tuned meshes under-resolve coplanar gaps and perimeter charge singularities while wasting tetrahedra in the far field; (2) standard Barsoum quarter-point edges along 1D rims must coexist with sharp $90^\circ / 270^\circ$ 2D sheet corners ($\alpha < 1/2$) without producing `0/0` Jacobian `NaN`s in field recovery; (3) finite simulation boxes truncated with PEC or homogeneous Neumann walls bias self-capacitance by $\mathcal{O}(L/R_{\text{box}})$ unless open-space decay is modeled directly; and (4) standard conforming Lagrange FEM gives only a one-sided upper bound $C^{\text{up}} \ge C^{\text{exact}}$ with no rigorous error certificate.

## Decision

* **Multi-pass goal-oriented electrostatic AMR (`AF-1`, `femcore::adapt::electro`, `femcore::meshing::sheet3d_amr`, `nyx-engine::amr`):**
  * `electro_term_indicators` and `electro_goal_indicators` evaluate cellwise volume, dielectric flux-jump, and Robin/Neumann residuals (including `CurvedTet` cells) in deterministic sorted-face order (`NYX_THREADS` bit-identical) and combine multi-terminal indicators with Cauchy–Schwarz weights $W_{ij} = 1/\max(|C_{ij}|, 10^{-3}\sqrt{C_{ii}C_{jj}})$.
  * `nyx-engine::amr` drives `adapt_passes > 1` (`adapt_rel_tol`, `adapt_min_converged_passes`, `adapt_growth`, `adapt_mode`: `auto | sheet_2d_z | full_3d`, `adapt_estimator`: `flux_jump_goal | prager_synge | hessian`), defaulting on 3D sheet stacks to column-projected `OverlayOracle2` adaptation (`electro_overlay_metrics_2d` with pre-gradation isotropic corner balls) + `extrude_overlay_2d_with_plan`, recording `AdaptPassRecord` in `CapacitanceMatrixData::adapt_history`, persisting the final adapted mesh in HDF5 (`open_results`), and forcing `adapt_passes = 1` inside `convergence_check`.
* **First-class quarter-point edges & Morrison–Lewis radical corner grading (`AF-2`, `femcore::meshing::singularity`, `CellMesh::qp_edges`):**
  * `set_mesh` exposes `singular_edges` (`"none" | "quarter_point"`), `quarter_point_xi` ($\in [0.25, 0.50]$), `corner_grading`, `adapt_grade`, and `adapt_max_cells`. `CellMesh::qp_edges` (and optional HDF5 `/mesh/qp_edge_pairs`) preserves QP edge tags across `MeshedModel` round-trips so `nyx-dfs` deduplicates enrichment.
  * `corner_grading = true` grades overlay edge sizes and coupled vertical `dz` with Morrison–Lewis corner exponents ($\alpha_{90^\circ} \approx 0.2966$, $\alpha_{270^\circ} \approx 0.8147$); `cell_gradient` evaluates singular vertex gradients safely away from collapsed determinants so `V` and `|E|` stay finite everywhere.
* **Asymptotic open/infinite `far_field` boundary condition & `add_chip_stack` (`AF-3`, `femcore::farfield`, `nyx-dfs`, `nyx-engine::gds`):**
  * `BcKind::FarField { order, decay_power, center }` (`"far_field"` on `dfs` `electrostatic`) assembles symmetric positive-definite 2D logarithmic Robin (`order = 1`), 3D radial Robin (`order = 1`), and 3D Bayliss–Turkel symmetric Laplace–Beltrami (`order = 2`, falling back to `order = 1` when `flux_bounds` is active) matrices around area-weighted conductor centroids, enabling $N = 1$ isolated-conductor self-capacitance without an artificial ground wall. `add_chip_stack` and `import_gds(lateral_pad=...)` construct padded 3D `(VacBot, Si, VacTop)` chip boxes directly from layout bounds.
* **Prager–Synge $H(\mathrm{div})$ equilibrated flux & certified Loewner bounds (`AF-4`, `femcore::rt`, `nyx-dfs::electro`, `nyx-analysis::capacitance`):**
  * `equilibrate_electro_flux_bounds` constructs $\mathrm{RT}_{p-1}$ ($\mathrm{RT}_0$ for $p=1$, $\mathrm{RT}_1$ for $p=2$) fluxes via Braess–Schöberl / Ern–Vohralík vertex-patch equilibration split into half-stars across interior conductive sheets, computing $C^{\text{low}} = Q (W^{\text{flux}})^+ Q^T \preceq C^{\text{exact}} \preceq C^{\text{up}}$, midpoint $C^{\text{mid}} = \tfrac{1}{2}(C^{\text{low}} + C^{\text{up}})$, and entrywise error bounds $\delta C_{ij}$.
  * `flux_bounds = "report" | "midpoint"` populates `maxwell_{lower,upper,mid}_matrix_f` and `maxwell_error_bound_f` (promoting $C^{\text{mid}}$ across `mutual_matrix_f`, `ec_matrix_hz`, `Solution::value`, and scalars when `"midpoint"`), propagates rigorous `[C_Sigma_low_F, C_Sigma_up_F]` and `[E_C_low_Hz, E_C_up_Hz]` brackets through `capacitance_matrix`, and rejects incompatible `quarter_point` or non-zero `rho_v` / `surface_charge` inputs with `Error::Invalid`.

## Consequences

* Coarse un-tuned 3D chip meshes automatically refine around coplanar gaps and sharp corners, open-boundary self-capacitances converge in compact $2\text{–}3\,R$ boxes, and every multi-terminal solve can emit certified two-sided capacitance and qubit $E_C$ brackets with $2.5\text{–}5\times$ midpoint error reduction.

## Code

`crates/femcore/src/{farfield,rt,fem}.rs`, `crates/femcore/src/mesh/{cells,exchange,mod}.rs`, `crates/femcore/src/adapt/electro.rs`, `crates/femcore/src/meshing/{singularity,sheet3d_amr,sheet3d,request,mod}.rs`, `crates/femcore/tests/{electro_amr,singularity_qp,rt_equilibration}.rs`, `crates/nyx-physics/src/cqed.rs`, `crates/nyx-dfs/src/{electro,problem,lib}.rs`, `crates/nyx-dfs/tests/farfield_electro.rs`, `crates/nyx-analysis/src/capacitance.rs`, `crates/nyx-analysis/tests/capacitance_bounds.rs`, `crates/nyx-engine/src/{amr,project,meshing,commands,gds,solve,engine,results_file,h5io,convergence,scriptout}.rs`, `crates/nyx-engine/tests/accuracy_features_e2e.rs`, `crates/nyx-script/tests/accuracy_features_script.rs`.
