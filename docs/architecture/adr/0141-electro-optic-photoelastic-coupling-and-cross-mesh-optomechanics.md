# 141. Domain-restricted piezoelectric submesh, cross-mesh electro-optic and photoelastic perturbation, and cross-mesh optomechanics

Date: 2026-11-27

## Status

Accepted (FR-39). Extends ADR 9, ADR 19, ADR 25, ADR 42, ADR 45, ADR 56, ADR 109, ADR 129, and ADR 135.

## Context

Coupled piezoelectric, electro-optic (Pockels), photoelastic (elasto-optic), and cavity/waveguide optomechanical devices span disparate physical domains and spatial scales:

1. **Piezoelectric resonators and actuators with air/vacuum gaps:** Monolithic piezoelectric `(u, phi)` models frequently include surrounding air, vacuum, or dielectric cladding domains where electrostatic potential $\phi$ must be solved for fringe fields and capacitive loading, but mechanical displacement $\mathbf{u}$ has zero physical stiffness. Allocating mechanical displacement unknowns $\mathbf{u}$ on air/vacuum nodes wastes degrees of freedom and requires artificial soft-spring stabilization.
2. **Electro-optic and photoelastic index modulation across independent meshes:** Electro-optic modulators (such as X-cut thin-film $\text{LiNbO}_3$, point group `3m`) and strain-tuned photonic waveguides/cavities perturb the optical relative impermeability tensor $\Delta(\boldsymbol{\varepsilon}_r^{-1})_{IJ} = \sum_k r^S_{Ik} E_{\text{DC},k} + \sum_J p_{IJ} S_J$ from an electrostatic field $\mathbf{E}_{\text{DC}}$ (`nyx-dfs` or `nyx-sms`) and/or mechanical strain tensor $\mathbf{S}$ (`nyx-sms`), often solved on different meshes tailored to each physics (`ADR 135`).
3. **Cross-mesh cavity optomechanics and Brillouin scattering without physics-crate coupling:** Computing vacuum optomechanical coupling rates $g_0 = g_{\text{PE}} + g_{\text{MB}}$ (photoelastic bulk integral $g_{\text{PE}}$ plus Johnson moving-boundary radiation-pressure surface integral $g_{\text{MB}}$) and stimulated Brillouin scattering (SBS) gain $G_B$ requires combining an optical mode (`nyx-ems`, meshed over core + air cladding + optical PML) with a mechanical mode (`nyx-sms`, meshed on the solid core + elastic PML only) inside `nyx-analysis`, which by crate layering (`ADR 26`) must never depend on `nyx-ems` or `nyx-sms`.

## Decision

* **Domain-restricted mechanical submesh in monolithic piezoelectric FEM (`nyx-sms`):**
  * Add optional `u_domains: Vec<i32>` (`"u_domains"`, region names where mechanical displacement $\mathbf{u}$ is active) to `nyx_sms::Problem` and `Sms::params()`. When `u_domains` is omitted, any non-conductor domain with dielectric permittivity (`eps_r` or anisotropic permittivity) and zero elastic stiffness (`E` / `stiffness` absent) automatically sets `mech_active = false, elec_active = true` in `piezo_mat::resolve_domain_mats_with_u_domains`.
  * In `piezo_asm::PiezoConstraints`, build the active mechanical node map `u_nodes` / `node_to_u` over cells with `mech_active == true` and allocate mechanical DOFs strictly on `0..n_mech` (`n_mech = dim * n_u_nodes`) while electrostatic potential $\phi$ spans `n_mech..(n_mech + n_nodes)`. All-solid models (`n_u_nodes == n_nodes`) use the identity numbering `0..n_nodes`, preserving exact 0-ULP parity.
  * Export 6-component wafer-frame solver-Voigt strain `"S"` (`"m/m"`, components `["xx", "yy", "zz", "xy", "yz", "xz"]`, role `"strain"`) across `nyx-sms` static, eigenfrequency, and piezoelectric studies, and export `"E"` (`"V/m"`, role `"electric_field"`) from piezoelectric studies.
* **Cross-mesh electro-optic and photoelastic permittivity perturbation (`femcore`, `nyx-physics`, `nyx-ems`):**
  * Add `ElectroOptic::to_rank3`, `ElectroOptic::rotate_voigt`, `electro_optic_global`, `photoelastic_global`, and `delta_eps_eo_pe` in `femcore::material::coupled`, evaluating the crystal-rotated relative impermeability perturbation in global wafer axes:
    $$\Delta(\varepsilon_r^{-1})_I = \sum_{k=1}^3 r^S_{Ik} E_{\text{DC},k} + \sum_{J=1}^6 p_{IJ} S^{\text{tb}}_J, \qquad \Delta \boldsymbol{\varepsilon}_r = -\boldsymbol{\varepsilon}_r \, \Delta(\boldsymbol{\varepsilon}_r^{-1}) \, \boldsymbol{\varepsilon}_r.$$
  * Preserve `photoelastic` and `electro_optic` tensors on `Anisotropy` during `nyx_physics::decode::materials_by_label` (including optically isotropic cubic crystals such as Silicon and inline `r11..r51`, `p11..p44` properties).
  * Extend `Ems::imports()` with `"electric_field"` (`"ems.electric_field"`) and `"strain"` (`"ems.strain"`). `nyx_ems::Problem::apply_transferred_inputs` samples cross-mesh `"electric_field"` and `"strain"` at cellcentroids strictly filtered by canonical `domain_id` (`femcore::transfer::TransferIndex`), and `Problem::cell_complex_permittivity` (plus `modes_complex`, `axisym`, and `mixed_ops`) adds $\Delta \boldsymbol{\varepsilon}_r$ cell-by-cell.
* **Pure `nyx-analysis` cross-mesh `"optomechanics"` analysis (`nyx-analysis::optomechanics`, `nyx-engine::analyse`):**
  * Add `OptomechanicsAnalysis` (`id: "optomechanics"`) and `optomechanics_with_meshes` in `nyx-analysis` operating purely on optical `(Mesh, SolveResult)` and mechanical `(Mesh, SolveResult)` pairs via `femcore::transfer::{TransferIndex, BoundaryIndex}` with zero dependency on `nyx-ems` or `nyx-sms`.
  * Normalize the selected mechanical eigenmode to unit peak displacement $\max \|\tilde{\mathbf{u}}(\mathbf{x})\|_2 = 1$, compute modal effective mass $m_{\text{eff}} = \int_{\Omega_{\text{mech}}} \rho \|\tilde{\mathbf{u}}\|^2 \,\mathrm{d}V$ and zero-point fluctuation amplitude $x_{\text{zpf}} = \sqrt{\hbar / (2 m_{\text{eff}} \Omega_m)}$.
  * Integrate the cross-mesh photoelastic frequency pull $G_{\text{PE}}$ over optical core cells by locating quadrature points strictly inside matching canonical solid `domain_id`s on the mechanical mesh, and integrate the cross-mesh moving-boundary frequency pull $G_{\text{MB}}$ over dielectric interface facets using Johnson's perturbation theorem:
    $$G_{\text{MB}} = -\frac{\omega_{\text{opt}}}{2 U_{\text{opt}}} \oint_{\partial \Omega_{\text{core}}} (\tilde{\mathbf{u}}\cdot\hat{\mathbf{n}}) \left[ \varepsilon_0 \mathbf{E}_\parallel^* \cdot (\boldsymbol{\varepsilon}_{r,1} - \boldsymbol{\varepsilon}_{r,2}) \mathbf{E}_\parallel - (\varepsilon_1^{-1} - \varepsilon_2^{-1}) |D_\perp|^2 \right] \mathrm{d}A,$$
    reporting $G_{\text{OM}} = G_{\text{PE}} + G_{\text{MB}}$, single-photon coupling rates $g_{\text{PE}}, g_{\text{MB}}, g_0 = g_{\text{PE}} + g_{\text{MB}}$ (in $\text{rad/s}$ and $\text{Hz}$), and waveguide stimulated Brillouin scattering gain $G_B$ ($\text{W}^{-1}\text{m}^{-1}$).

## Consequences

* Piezoelectric models with dielectric air/vacuum domains allocate zero mechanical DOFs outside the active solid submesh while retaining exact 0-ULP parity on all-solid models.
* Electro-optic (`dfs -> ems`, `sms -> ems`), photoelastic (`sms -> ems`), and optomechanical (`ems + sms -> nyx-analysis`) workflows run seamlessly on shared or independent per-physics meshes without breaking crate layering.

## Code

`crates/femcore/src/material/{coupled,spec,mod}.rs`, `crates/nyx-physics/src/{validate,decode,testkit}.rs`, `crates/nyx-sms/src/{problem,piezo_mat,piezo_asm,piezo,piezo_driven,solid,lib}.rs`, `crates/nyx-ems/src/{problem,modes_complex,axisym,mixed_ops,lib}.rs`, `crates/nyx-analysis/src/{optomechanics,lib}.rs`, `crates/nyx-engine/src/analyse.rs`, `crates/nyx-sms/tests/fr39_piezo_submesh.rs`, `crates/nyx-ems/tests/fr39_eo_pe.rs`, `crates/nyx-analysis/tests/fr39_optomechanics.rs`, `crates/nyx-engine/tests/fr39_cross_mesh_eoom.rs`.
