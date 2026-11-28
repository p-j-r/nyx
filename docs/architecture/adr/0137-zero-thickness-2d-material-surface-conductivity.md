# 137. Zero-thickness 2D materials and tensor surface conductivity sheets (`FR-35`)

Date: 2026-11-27

## Status

Accepted. Extends ADR 29 & ADR 30 & ADR 108 & ADR 109 (material TOML `format = 2` specification and schema), ADR 45 & ADR 46 (`nyx-ems` 2D/3D electromagnetic kernels), ADR 102 (`Op::Sheet` conformal 2D surface meshes), ADR 104 & ADR 125 (3D surface impedance & sheet current extraction), ADR 106 & ADR 129 (`ImageProm` multi-port frequency-domain sweeps), ADR 128 (`MixedDiscretisation`), and ADR 131–133 (VVUQ `Source` provenance, `T_range`, and determinism).

## Context

Atomically thin two-dimensional materials (monolayer and few-layer graphene, transition-metal dichalcogenides such as $\text{MoS}_2$ and $\text{WS}_2$, and anisotropic van der Waals sheets such as black phosphorus) have physical thicknesses $d \approx 0.34\text{--}0.7\,\text{nm}$, roughly $10^3\text{--}10^4\times$ thinner than infrared or terahertz wavelengths and dielectric waveguide cross-sections. Volumetric 3D meshing of sub-nanometer layers causes severe aspect-ratio locking and unnecessary degree-of-freedom explosion, whereas representing the 2D layer as a zero-thickness interior or exterior surface $\Gamma_{2\text{D}}$ carrying a $2\times 2$ tangential surface conductivity tensor $\boldsymbol{\sigma}_{2\text{D}}(\omega, T, \mu_c, B_0)$ (in Siemens, $\text{S}$) introduces zero auxiliary unknowns and requires no volumetric sub-nanometer cells.

Unlike metallic or superconducting surface-impedance boundaries (`surface_impedance`), a 2D material sheet is penetrable ($\mathbf{n}\times(\mathbf{E}^+ - \mathbf{E}^-) = \mathbf{0}$, so tangential electric-field Nédélec degrees of freedom must remain free rather than clamped to PEC) and induces a tangential magnetic-field jump $\mathbf{n}\times(\mathbf{H}^+ - \mathbf{H}^-) = \mathbf{K}_{2\text{D}} = \boldsymbol{\sigma}_{2\text{D}}(\mathbf{x})\,\mathbf{E}_t$. Furthermore:

1. Graphene requires the analytical intraband + interband Kubo conductivity formula across arbitrary chemical potential $\mu_c$, relaxation time $\tau$, temperature $T$, Fermi velocity $v_F$, and optional out-of-plane magnetostatic bias $B_0$ (producing gyrotropic Hall off-diagonal terms $\sigma_{xy} = -\sigma_{yx}$).
2. Anisotropic 2D crystals require principal-axis Drude/Lorentz oscillator or PCHIP-interpolated tabulated complex sheet conductivities with in-plane rotation $\theta_{2\text{D}}$.
3. Electrostatically gated 2D waveguide modulators require multiphysics coupling from an electrostatic potential field (`dfs.V` $\to$ `ems.potential`) into spatially varying carrier density $n_s(\mathbf{r}) = C_{\text{ox}}(V(\mathbf{r}) - V_{\text{CNP}})/e$ and chemical potential $\mu_c(\mathbf{r}) = \hbar v_F \sqrt{\pi |n_s(\mathbf{r})|}$ at each surface quadrature point.

## Decision

* **2D material constitutive models and TOML `format = 2` `[sheet_2d]` section (`femcore::material::sheet2d`):**
  * Added `Sheet2d`, `Sheet2dModel`, and `Sheet2dPole` in `crates/femcore/src/material/sheet2d.rs`, wired into `MaterialSpec::sheet_2d`, `material_schema()` (`material_schema.v2.json`), and unit conversion (`Dim::SheetConductivity` `"S"`/`"mS"`/`"uS"`, `Dim::Time`, `Dim::Velocity`, `Dim::MagneticFluxDensity`, `Dim::CapacitancePerArea`).
  * Supports three 2D constitutive models under the $\exp(+i\omega t)$ sign convention ($\sigma = \sigma' - i\sigma''$, $\sigma' \ge 0$ passive loss, $\sigma'' > 0$ inductive):
    1. `kubo_graphene`: exact thermal intraband factor $\tilde{\mu}_c(T) = 2 k_B T \ln(2\cosh(\mu_c / (2 k_B T)))$, classical cyclotron frequency $\omega_c = e B_0 v_F^2 / \max(|\mu_c|, \tilde{\mu}_c)$, closed-form $T \to 0$ logarithmic interband conductivity, and nonsingular subtracted 32-point Gauss–Legendre quadrature over $[0, \xi_{\max}]$ with exact analytic asymptotic tail for $T > 0$, scaled by `n_layers`.
    2. `anisotropic`: principal-axis background sheet conductivities $(\sigma_{uu,0}, \sigma_{vv,0})$ plus sums of Drude/Lorentz 2D poles (`poles_u`, `poles_v`), rotated in-plane by $\mathbf{R}(\theta_{2\text{D}})\,\boldsymbol{\sigma}_{\text{prin}}\,\mathbf{R}(\theta_{2\text{D}})^T$.
    3. `tabulated`: monotone PCHIP-interpolated complex principal conductivities $(\sigma_{uu}(\omega), \sigma_{vv}(\omega))$ over ascending angular frequency points.
  * Added electrostatic gate-tuning helpers `gate_oxide_capacitance(eps_ox, d_ox)` and `chemical_potential_from_gate_potential(v_gate, v_cnp, c_ox, v_f)`, plus closed-form verification oracles `sheet2d_fresnel_normal_tra`, `tm_plasmon_wavenumber_nonretarded`, and `tm_plasmon_wavenumber_retarded_symmetric`.
  * Enforces VVUQ provenance (`Source::is_unverified()` emits `unverified_material_provenance` warning) and temperature validity range (`T_range` clamping emits `material_temperature_clamped` warning).
  * Added curated `materials/graphene.toml` (`format = 2`) with literature provenance citations (Falkovsky & Varlamov 2007, Gusynin et al. 2007, Hanson 2008, Liu et al. 2011) and `T_range = [1.0, 1000.0]`.
* **Penetrable variational surface jump integral in `nyx-ems` (`ibc`, `mixed_ops`, `cavity`, `driven`, `bloch`, `modes`, `modes_complex`, `axisym`):**
  * Registered boundary condition `sheet_2d` (`Sheet2dBc`) across `modes`, `axisymmetric_modes`, `cavity`, `bloch`, and `driven` studies, and added optional field link `FieldLink::new("ems.potential", "potential")` on `Ems::imports()`.
  * On 3D `SimplexMesh` and `MixedDiscretisation` (`cavity`, `driven`, `bloch`), `assemble_sheet2d_3d` and `assemble_sheet2d_3d_mixed` deduplicate interior and exterior sheet facets (`Op::Sheet`), keep all Nédélec tangential DOFs free (excluding exterior `sheet_2d` facets from default exterior PEC walls and never collapsing sheet vertices into conductor clusters), interpolate `imported_potential` at surface quadrature points, and assemble the penetrable surface bilinear form:
    $$a_{2\text{D}}(\mathbf{w}, \mathbf{E}) = +i\omega\mu_0 \int_{\Gamma_{2\text{D}}} \mathbf{w}_t \cdot \big(\boldsymbol{\sigma}_{2\text{D}}(\mathbf{x})\,\mathbf{E}_t\big)\,\mathrm{d}S$$
    into $\mathbf{K}_{\text{IBC}} \mathrel{+}= -\omega_{\text{ref}}\mu_0 \operatorname{Im}(\boldsymbol{\sigma}_{2\text{D}})$ and $\mathbf{C}_{\text{IBC}} \mathrel{+}= \mu_0 c_0 \operatorname{Re}(\boldsymbol{\sigma}_{2\text{D}})$, while recording `SheetCellKind::Sheet2D` for Stage-3 surface current $\mathbf{K} = \boldsymbol{\sigma}_{2\text{D}}\mathbf{E}_t$ extraction.
  * When $B_0 \neq 0$ ($\sigma_{uv} \neq \sigma_{vu}$), `DrivenAssembly` sets `has_non_symmetric_sheets = true` so `PromOperator::is_symmetric()` returns `false`, enabling full non-symmetric Galerkin projection in `ImageProm` alongside general sparse LU factorization.
  * On 2D waveguide cross-sections (`modes` / `modes_complex`) and axisymmetric resonators (`axisym`), `add_sheet2d_2d_waveguide[_mixed]` and `add_sheet2d_2d_axisym[_mixed]` assemble transverse $E_\tau$ (`a_coo`), longitudinal $E_y = j\beta\phi$ / azimuthal $E_\theta = -i\Psi/r$ (`b_coo` / `a_coo`), and off-diagonal $(\tau, y)$ coupling (`a1_coo`) on both `SimplexMesh` and `MixedDiscretisation`.

## Consequences

* Graphene, TMD, and anisotropic 2D sheets can be placed on any interior `Op::Sheet` interface or boundary facet in 2D waveguide, 2D axisymmetric, and 3D cavity/Bloch/driven simulations without volumetric sub-nanometer meshing or auxiliary DOFs.
* Electrostatic gate potential fields solved in `nyx-dfs` couple directly into spatially varying Kubo graphene chemical potential $\mu_c(\mathbf{r})$ in `nyx-ems`.

## Code

`crates/femcore/src/material/{sheet2d,mod,spec,schema,units,reader}.rs`, `materials/graphene.toml`, `crates/nyx-ems/src/{lib,problem,ibc,mixed_ops,cavity,driven,bloch,modes,modes_complex,axisym}.rs`, `crates/nyx-ems/tests/sheet2d_graphene.rs`.

## Contract

`material_schema` 2.0 (`material_schema.v2.json`, additive `[sheet_2d]`), `physics_ems` 1.9 (`physics_ems.v1.json`, additive `sheet_2d` BC and `potential` import).
