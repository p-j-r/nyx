# 144. Curved and mixed-cell Local Time-Stepping Discontinuous Galerkin time-domain EM solver (`Study::Dgtd`)

Date: 2026-11-27

## Status

Accepted (FR-42 Track A). Extends ADR 5, ADR 8, ADR 9, ADR 11, ADR 19, ADR 25, ADR 45, ADR 82, ADR 106, ADR 128, ADR 131, ADR 136, ADR 137, and ADR 138.

## Context

Broadband transient electromagnetic propagation, ultra-wideband S-parameter extraction, pulse scattering on curved and multi-scale sub-wavelength geometries, and time-domain ringdown analysis complement `nyx-ems`'s frequency-domain (`modes`, `cavity`, `bloch`, `driven`, `eme`) solvers:

1. **Unstructured curved and mixed-cell geometries without staircasing:** Sub-wavelength curved dielectric/metallic interfaces and conformal layered meshes (`Tri`, `Quad`, `Tet`, `CurvedTet`, `Prism`, `Hex` via `MixedDiscretisation`, ADR 82, ADR 128) require a discontinuous finite-element time-domain formulation that avoids Cartesian staircasing while keeping element mass matrices strictly block-local so explicit time stepping performs zero global linear solves.
2. **Energy conservation and penetrable zero-thickness 2D material sheets in time domain:** On general curved and mixed meshes, discrete volume differentiation and surface integration must satisfy a Summation-By-Parts (SBP) skew-symmetric identity so that central numerical fluxes conserve discrete electromagnetic energy to machine precision on closed PEC/PMC cavities, while Rankine–Hugoniot upwind fluxes natively incorporate zero-thickness 2D conductive sheet boundary/interface conditions (`sheet_2d`, ADR 137) via $\hat{\mathbf{n}} \times (\mathbf{H}^+ - \mathbf{H}^-) = \boldsymbol{\sigma}_{2\text{D}} \mathbf{E}_t^*$ without auxiliary sheet degrees of freedom or CFL degradation.
3. **Multiscale CFL bottleneck (`h_max / h_min > 50`):** In nanogap resonators, slot waveguides, and sharp metallic corners, the smallest element diameter $h_{\min}$ can be 50–100$\times$ smaller than the bulk element size $h_{\max}$. A uniform global explicit time step $\Delta t = \min_K \Delta t_K$ forces bulk coarse elements to step at the nanogap CFL limit, wasting an order of magnitude of compute work.

## Decision

* **Three-stage high-order curved/mixed-cell DGTD kernel (`nyx-ems::dgtd`, `Study::Dgtd` / `"dgtd"`):**
  * **`assemble_dgtd` (Stage 1):** Supports 2D (`Full`, `Tm`, `Te`) and 3D full-vector Maxwell systems at polynomial basis orders $p \in \{1, 2, 3\}$ across `Tri`, `Quad`, `Tet`, `CurvedTet`, `Prism`, and `Hex` cells. Uses Gauss–Lobatto–Legendre collocation (`femcore::linalg::matfree::Gll`) on `Quad` and `Hex` cells for exact diagonal element mass inversion $M_K^{-1}$, and high-order symmetric quadrature with local $N_p \times N_p$ dense factorization on `Tri`, `Tet`, `CurvedTet`, and `Prism` cells (zero global linear solves). Forms strong-form skew-symmetric SBP differentiation matrices $D^k_K = M_K^{-1}\hat{K}^k_K$ where $\hat{K}^k_K = \frac{1}{2}(K^{\text{vol},k}_K - (K^{\text{vol},k}_K)^T) + \frac{1}{2} F^k_K$, guaranteeing $\hat{K}^k_K + (\hat{K}^k_K)^T = F^k_K$ on straight and curved cells alike.
  * **Riemann fluxes, `sheet_2d` interface jump solver, ADE dispersion, `[nonlinear]` constitutive updates, and FEM-seeded modal ports:** Implements selectable Rankine–Hugoniot upwind ($\alpha = 1$) and energy-conserving central ($\alpha = 0$) numerical fluxes across interior and boundary faces (`Pec`, `Pmc`, `SilverMuller`, `Port`). Integrates penetrable zero-thickness `Sheet2dBc` (`FR-35`) interfaces directly into the Riemann solver via $\big((Y^- + Y^+)\mathbf{I}_2 + \boldsymbol{\sigma}_{2\text{D}}\big)\mathbf{E}_t^* = Y^-\mathbf{E}_t^- + Y^+\mathbf{E}_t^+ + \alpha\,\hat{\mathbf{n}}\times(\mathbf{H}^+ - \mathbf{H}^-)$. Supports multi-pole Drude, Lorentz, and gyrotropic cyclotron-precession (`FR-34`) Auxiliary Differential Equations (ADEs), shared format-2 `[nonlinear]` (`ADR 143`) $\chi^{(2)}$/QPM poling and $\chi^{(3)}$ Kerr + TPA/FCA (`drude_carrier_step`) constitutive updates, uniaxial PML domains, Total-Field/Scattered-Field (TFSF) plane-wave surface injection, point dipoles, and 2D/3D waveguide modal ports solved via `nyx_ems::modes`.
  * **`solve_dgtd` (Stage 2) and power-of-two Local Time-Stepping (`Lts`):** Partitions elements into $L$ power-of-two CFL tiers ($\Delta t_\ell = 2^\ell \Delta t_0$, $\ell \in 0..L$) with 1-ring 2:1 tier grading ($|\ell(K^-) - \ell(K^+)| \le 1$) and advances tiers via a recursive 4th-order multi-rate Runge–Kutta LTS hierarchy (`DgtdTimeStepper::Lts`, plus uniform 5-stage 4th-order Carpenter–Kennedy `Lserk4` and classical 4-stage `Rk4` baselines). Element and face evaluations execute in fixed chunks via `femcore::par`, guaranteeing bitwise-identical fields and monitors across `NYX_THREADS=1` vs `8` (ADR 11).
  * **`dgtd_fields` (Stage 3):** Accumulates on-the-fly running DFTs for `FieldMonitor`, surface Poynting `FluxMonitor`, and waveport modal $(V_p(\omega), I_p(\omega)) \to (a_p(\omega), b_p(\omega))$ power waves to emit broadband `DrivenPortData` (`S`, `Z`, `Y` matrices) alongside `DgtdData` (`DGTD_DATA_KEY`, containing LTS tier histograms, work/runtime speedup metrics, energy history $W(t)$, and `FieldTimeMonitor` probe waveforms).

## Consequences

* `"ems"` remains the single 3-letter electromagnetic physics crate (`crates/nyx-engine/tests/layering.rs`, ADR 12, ADR 25), and `Study::Dgtd` shares the existing `waveport`, `lumped_port`, `sheet_2d`, `pml`, `pec`, `pmc`, and scattered-field excitation descriptors.
* `physics_ems`, `physics_ems_result`, and `kernel_stages` advance additively to include `"dgtd"` and `nyx_ems::dgtd::{assemble_dgtd, solve_dgtd, dgtd_fields, run_dgtd}`.

## Code

`crates/nyx-ems/src/{dgtd,problem,sensitivity,lib}.rs`, `crates/nyx-ems/tests/{fr42_dgtd,contracts,finite_guard}.rs`, `crates/nyx-engine/tests/phase3_integration.rs`.
