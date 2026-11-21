# 104. BCS superconducting surface impedance, thin-film IBC, and rank-1 lumped Josephson ports

Date: 2026-11-21

## Status

Accepted (superconducting qubit milestone M3 & M4.1). Extends the TOML material library of ADR 29 / ADR 30 and the `ems` kernels of ADR 45 / ADR 67.

## Context

In disordered high-kinetic-inductance films ($35\text{–}100\ \text{nm}$ Nb, NbN, TiN, Ta, granular Al), magnetic energy stored inside the film exceeds vacuum field energy and down-shifts resonator frequencies by tens to hundreds of megahertz with nonlinear temperature dependence up to $T_c$. Treating sheets as PEC misses kinetic inductance, while volumetric meshing of a $35\text{ nm}$ film in a $10\text{ mm}$ die is intractable. Furthermore, simulating linearized Josephson junctions inside 3D Nedelec cavities via distributed edges without a rank-1 circuit port DOF causes galvanic floating-island gradient nullspace pathologies and prevents fast passive reduced-order re-projection (`CircuitRom`).

## Decision

* **`[superconducting]` TOML schema & Zimmermann BCS evaluator (`crates/femcore/src/material/bcs.rs`, `materials/*.toml`).** Extends material format 1 with `[superconducting]` (`tc`, optional `gap_0`, `lambda_l`, `xi_0`, `mean_free_path`, `sigma_n`), adds six shipped superconducting TOMLs (`aluminium`, `niobium`, `tantalum`, `nbn`, `tin`, `granular_aluminium`), implements arbitrary-purity Zimmermann BCS quadrature ($\sigma(\omega, T) = \sigma_1 - i\sigma_2$) plus `kerr` (Park et al. 2026 thin-film IBC including London depth & normal skin correction) and `two_sided` film models, and exposes the read-only engine command `bcs_surface_impedance`.
* **`pec` & `surface_impedance` on single-facet sheets (`crates/nyx-ems/src/ibc.rs`).** Adds `BcKind::Pec` and `BcKind::SurfaceImpedance` (supporting either explicit `L_k_sq`/`R_s` or a named superconducting material + `thickness`, `temperature`, `f_ref`, `film_model`) across 3D `cavity` / `driven` and 2D `modes` / `modes_complex` / `axisym`. On interior sheets (`"<sheet>"`), `ibc.rs` deduplicates triangles shared by two cells and maps series surface impedance $Z_s = R_s + i\omega L_{k,\square}$ to positive parallel admittance $Y_{\text{eff}} = L_{\text{eff}}^{-1}/(i\omega) + G_{\text{eff}}$ integrated cleanly on interior Nedelec traces.
* **Rank-1 `lumped_jj` & `CircuitRom` export (`crates/nyx-ems/src/{cavity,rom,ibc}.rs`).** Adds `BcKind::LumpedJj { name, l_j, e_j_hz, c_j, e_l_hz, phi_ext, n_g, direction, formulation }` (`formulation = "rank1"` default, or `"distributed_sheet"`):
  - `distributed_sheet` contributes $(\ell/w) C_J \mathbf{B}_J$ and $L_J^{-1}$ directly in $H(\text{curl})$.
  - `rank1` (default) projects a uniform sheet current vector $\mathbf{w}_j$ onto each junction sheet (`1` electrical degree of freedom per junction), includes its lumped $C_J$ and $L_J$, augments `g_red` nullspace elimination safely for floating superconducting islands (so floating pads separated by a gap remain clean electrostatic/circuit modes), extracts physical magnetic field $\mathbf{H}$ and sheet surface current $\mathbf{K}_s$ normalized precisely to total electric energy $U_{\text{elec}} = 1\ \text{J}$, computes exact 3D island-tied $H^1$ port capacitance `C_port`, and builds a passive Craig–Bampton [`nyx_physics::CircuitRom`] containing static + dynamic modes, $(\mathbf{K}_r, \mathbf{C}_r, \mathbf{M}_r, \mathbf{B}_r, \mathbf{C}_{\text{port}})$, and per-mode volume magnetic vs. IBC kinetic vs. JJ inductive energies (`TlsMeta` metadata also emitted for `tls_participation`).

## Consequences

* 3D $\lambda/2$ CPW Nb resonator shifts match the $-154\text{ MHz}$ Park et al. (2026) Fig. 1e curve within `< 0.6 MHz` across $T/T_c \in [0.01, 0.85]$.
* Every 3D `ems` `cavity` solve preserves modal energy equipartition to machine precision (`equipartition_error < 1e-12`, typically `1e-13–1e-14`) even with floating islands and kinetic IBC sheets.

## Code

`crates/femcore/src/material/{bcs,spec,reader,schema,units,mod}.rs`, `materials/{aluminium,niobium,tantalum,nbn,tin,granular_aluminium}.toml`, `crates/nyx-physics/src/{cqed,validate,decode}.rs`, `crates/nyx-ems/src/{ibc,rom,cavity,modes,modes_complex,axisym,lib,problem}.rs`, `crates/nyx-engine/src/materials.rs`; tests `crates/femcore/tests/{material_bcs,material_contract,material_library,materials}.rs`, `crates/nyx-ems/tests/kicq_cpw_bcs.rs`.
