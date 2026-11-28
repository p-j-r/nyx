# 139. Nonlinear thermal radiation, steady Joule conduction, optical heat exports, and Anderson coupled loops

Date: 2026-11-27

## Status

Accepted (FR-38). Extends ADR 9, ADR 19, ADR 25, ADR 41, ADR 45, ADR 61, ADR 106, ADR 109, and ADR 133.

## Context

Self-consistent electrothermal and photothermal device modeling requires four capabilities working together across `femcore`, `nyx-dfs`, `nyx-ems`, `nyx-physics`, and `nyx-engine` without violating strict crate layering (`nyx-dfs` and `nyx-ems` never depend on each other):

1. **Nonlinear radiative heat transfer and temperature-dependent thermal properties:** High-temperature and vacuum-packaged photonic, MEMS, and RF components radiate heat across boundaries according to the Stefan–Boltzmann $T^4$ law, $-\mathbf{n}\cdot\mathbf{q} = \varepsilon F \sigma_{\mathrm{SB}} (T^4 - T_{\mathrm{amb}}^4)$, while bulk thermal conductivity $k(T)$ and heat capacity $C_p(T)$ vary strongly over cryogenic-to-elevated temperature spans.
2. **Steady DC ohmic conduction and Joule heating:** Resistive heaters, micro-hotplates, and interconnects require solving $-\nabla\cdot(\sigma(T)\nabla V) = 0$ in 2D and 3D with temperature-dependent electrical conductivity $\sigma(T)$ and conservative volumetric Joule dissipation $Q_{\mathrm{joule}} = \nabla V \cdot \sigma(T)\nabla V$.
3. **Electromagnetic volumetric and surface heat generation:** Driven frequency-domain (`ems::driven`), cavity (`ems::cavity`), and complex waveguide (`ems::modes`) solves with dielectric/ohmic loss or resistive `surface_impedance` sheets dissipate power that must map conservatively onto nodal thermal load vectors while excluding artificial absorbing-layer (`pml`) damping and responding to spatially varying temperature fields $T(\mathbf{x})$.
4. **Two-way self-consistent multiphysics coupling:** Once `dfs` exports `temperature` and imports `heat_source` while `ems` exports `heat_source` and imports `temperature`, `Registry::solve_order` rightly detects a bidirectional dependency cycle. Accidental cycles must remain rejected by default while projects that explicitly opt into self-consistent coupling converge rapidly via Anderson acceleration.

## Decision

* **Exact-tangent Newton thermal solver (`nyx-dfs::heat`, `femcore::material`):**
  * Add `BcKind::Radiation { emissivity, t_ambient, view_factor }` (`"radiation"`) on `heat_stationary` and `heat_transient` with Stefan–Boltzmann constant $\sigma_{\mathrm{SB}} = 5.670\,374\,419 \times 10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}$.
  * Add analytical derivative evaluation `TempCurve::eval_deriv` and `TempTensor::eval_matrix_deriv_global` in `femcore::material`, plus optional `Electrical::conductivity_t` (`[electrical] conductivity_T` / curve tables) audited by `nyx_physics::decode::audit_active_materials` and `nyx_engine::material_warnings` (`T_range` / `strict_material_range` aware per ADR 133).
  * When `Radiation` or temperature-dependent $k(T)$ / $C_p(T)$ is active, assemble a fixed union CSR sparsity template once (`symbolic_analyses == 1`) and solve `heat_stationary` and each implicit BDF1/BDF2 step of `heat_transient` via `femcore::linalg::snes::newton` with exact boundary tangent $+4\varepsilon F \sigma_{\mathrm{SB}} T^3$ and volume $k'(T)$ / $C_p'(T)$ tangents, positive-definite linearized warm start, and positive-Kelvin line-search safeguarding. Linear problems keep the existing single-factorization SPD fast path bit-for-bit.
* **Steady conduction and conservative `heat_source` (`nyx-dfs::conduction`):**
  * Add `Study::Conduction` (`"conduction"`) in `nyx-dfs` with `voltage` (and `potential`), `current_density`, and `ground` boundary conditions in 2D and 3D on simplex and native mixed meshes, evaluating $\sigma(T(\mathbf{x}))$ at quadrature points from imported `"temperature"` (`dfs.T`) and exporting `V`, `E`, `J`, scalar `"P_joule (W)"`, and conservative mass-projected `"Q_joule"` (`"W/m^3"`, role `"heat_source"`).
  * Add `FieldLink::new("dfs.Q_joule", "heat_source")` to `Dfs::exports()` and `FieldLink::new("dfs.heat_source", "heat_source")` + `FieldLink::new("dfs.temperature", "temperature")` to `Dfs::imports()`.
* **Optical/RF heat export and spatially varying thermo-optics (`nyx-ems::heat_source`):**
  * Export conservative nodal `"Q_opt"` (`"W/m^3"`, `FieldLink::new("ems.Q_opt", "heat_source")`) and `"P_loss (W)"` on lossy `ems::driven` sweeps (and normalized `"Q_opt"` on lossy `ems::cavity` and `ems::modes` solves), combining physical non-PML dielectric dissipation $Q_{\mathrm{vol}} = \frac{1}{2}\omega\varepsilon_0 \operatorname{Im}(-\mathbf{E}^H \tilde{\boldsymbol{\varepsilon}}_r \mathbf{E})$ with boundary `surface_impedance` sheet loss mapped onto adjacent boundary vertices via lumped $P_1$ mass projection and exact linear $P_2$ prolongation.
  * Evaluate library `[optical]` (`dn_dT`) and `[dielectric]` (`eps_T`, `loss_tangent`) permittivities cell-by-cell against spatially varying imported `"temperature"` vectors across `nyx-ems`.
* **Opt-in Anderson coupled loops (`nyx-physics::registry`, `nyx-engine::{project,solve,commands}`):**
  * Keep `Registry::solve_order` rejecting unconfigured dependency cycles by default.
  * Add `CoupledLoopConfig { max_iters, rtol, anderson_m, auxiliary, source_scale }` (`#[serde(default, skip_serializing_if = "Option::is_none")]` on `Project`, command `set_coupled_loop`) and `nyx_engine::solve::solve_coupled` driven by `femcore::linalg::snes::anderson` (default history depth `m = 5`) enforcing relative $L_\infty$ convergence $\|T^{k+1} - T^k\|_\infty / \max(1, \|T^k\|_\infty) \le \mathrm{rtol}$ for both cross-physics photothermal (`ems <-> dfs`) and intra-physics electrothermal (`dfs::conduction <-> dfs::heat_stationary`) loops.

## Consequences

* Crate layering remains strictly acyclic: `nyx-dfs` and `nyx-ems` exchange `"temperature"` and `"heat_source"` solely through `nyx-physics` role links orchestrated by `nyx-engine::solve`.
* Existing linear thermal solves, lossless electromagnetic solves, and acyclic one-way chains retain identical behavior and golden outputs.

## Code

`crates/femcore/src/material/{state_t,spec}.rs`, `crates/nyx-physics/src/{decode,registry}.rs`, `crates/nyx-dfs/src/{conduction,heat,problem,lib}.rs`, `crates/nyx-ems/src/{heat_source,problem,mixed_ops,driven,cavity,modes,modes_complex,lib}.rs`, `crates/nyx-engine/src/{project,solve,commands,engine,material_warnings,scriptout,lib}.rs`, `crates/nyx-dfs/tests/fr38_thermal.rs`, `crates/nyx-engine/tests/{fr38_coupled_thermal,vvuq_phase4_materials}.rs`.
