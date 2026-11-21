# 109. Material `format = 2`: multi-physics constitutive models, Mattis–Bardeen BCS surface impedance, and coupled-field tensors

Date: 2026-11-21

## Status

Accepted. Amends ADR 29, 30, 34, 35, and 101.

## Context

`format = 1` (ADR 29, 30) represented single-temperature elastic stiffness (`[elastic]`), second-rank thermal and electrical tensors (`[thermal]`, `[electrical]`), and optical dispersion (`[optical]`). The initial multi-source library and importers (ADR 104 and 108) stayed on `format = 1` by storing microwave loss tangents, superconducting parameters ($T_c$, $\Delta_0$, $\lambda_L$, $\xi_0$, $L_k$), interface participation thicknesses, piezoelectric coupling ($e_{iJ}, d_{iJ}$), Brillouin photoelasticity ($p_{IJ}$), Pockels electro-optic coefficients ($r_{IJ}$), and Akhiezer acoustic viscosity ($\eta_{IJ}$) as prose inside `notes`, and by flattening `refractiveindex.info` Formula 4 and orthogonal `k_data` grids during import.

That workaround meant the solvers could not evaluate temperature-dependent cryostat sweeps ($10\text{ mK} \le T \le 400\text{ K}$), complex microwave dielectric loss tangents, complex Mattis–Bardeen superconducting surface impedance $Z_s(\omega, T, d) = R_s + i\omega L_{s,\square}$, or bulk Akhiezer viscous damping directly from a material file.

## Decision

* **`format = 2` canonical schema (`femcore::material`) with `format = 1` read-only shim:**
  * `MaterialSpec` is bumped to `format = 2` (`docs/architecture/contracts/material_schema.v2.json`, `material_kitchen_sink.v2.toml`/`.json`, `material_kitchen_sink_models.v2.toml`/`.json`) while keeping `material_schema.v1.json` and `.v1.json` fixtures immutable.
  * `reader::read` accepts `format = 1` (producing a default-populated `format = 2` `MaterialSpec`) and `format = 2` (emitting canonical `format = 2` on `write`), rejecting all other version integers with `Reason::UnsupportedFormat`.
  * Top-level metadata adds `kind` (`dielectric`, `metal`, `superconductor`, `semiconductor`, `piezoelectric`, `fluid`), crystallographic `point_group` (32 Hermann–Mauguin classes), `crystal_setting` (IEEE Std 176-1987 axis convention), `reference_temperature` (default $293.15\text{ K}$), and `residual_stress` (signed thin-film stress in Pa).
  * `extends` inheritance merges every `format = 2` section (`[dielectric]`, `[superconducting]`, `[interfaces]`, `[piezoelectric]`, `[photoelastic]`, `[electro_optic]`, `[viscoelastic]`) cleanly at the table level.
* **Continuous temperature curves (`state_t.rs`) and lossless optical dispersion (`optics.rs`):**
  * `[thermal]` (`conductivity_T`, `specific_heat_T`, `expansion_T`), `[elastic]` (`c_T`), and `[dielectric]` (`loss_tangent_T`) support `tabulated_T` (PCHIP or linear), `log10_poly` (NIST cryogenic $\log_{10} y = \sum a_k (\log_{10} T)^k$), and `einstein_shift`.
  * `[optical]` adds native `rii_formula_4`, `rii_formula_5`, and `rii_formula_6`, decoupling of `n` and `k_data` wavelength grids, and cryogenic `dn_dT_model = { model = "einstein", ... }`.
* **Microwave `[dielectric]`, `[interfaces]`, and Mattis–Bardeen `[superconducting]` (`bcs.rs`):**
  * `[dielectric]` captures microwave/RF `permittivity`, `loss_tangent`, `loss_tangent_T`, and optional TLS parameters (`tls_tan_delta_0`, `tls_E_c`, `tls_T_0`, `tls_beta`), returning $\tilde{\varepsilon}_r(T) = \varepsilon_r (1 - i\tan\delta(T))$.
  * `[interfaces.{ma,ms,sa}]` stores thin-layer thickness $t$ and loss tangent $\tan\delta$ for metal–air, metal–substrate, and substrate–air participation analyses.
  * `[superconducting]` stores $T_c$, optional $\Delta_0$ (default $\frac{1}{2}\text{bcs\_gap\_ratio}\cdot k_B T_c$), `bcs_gap_ratio` (default $3.528$), $\lambda_{L,0}$, $\xi_0$, $\rho_n$ (or $\sigma_n$), mean free path $\ell$, optional film thickness $d$, `regime` (`local`, `dirty`, `pippard`), and `subgap_loss_floor` (`tan_delta_0` / Dynes $\Gamma$ so $R_s(T\to 0) > 0$ down to $10\text{ mK}$).
  * `bcs::bcs_surface_impedance(sc, f_hz, t_k, d_override)` evaluates the Mattis–Bardeen complex conductivity $\sigma(\omega, T) = \sigma_1 - i\sigma_2$ via Gauss–Legendre quadrature over the thermal distribution and singularity-subtracted sub-gap integral, returning $Z_s = R_s + i\omega L_{s,\square}$ with thin-film $\coth(d/\lambda_{\text{eff}})$ correction.
* **Coupled-field and viscoelastic tensors (`coupled.rs`):**
  * `[piezoelectric]` ($3\times 6$ $e_{iJ}$ in $\text{C/m}^2$ or $d_{iJ}$ in $\text{C/N}$ with Voigt-weighted conversion $e = d\,C^E$ and thermodynamically linked $C^D, C^E, \varepsilon^S, \varepsilon^T$), `[photoelastic]` ($6\times 6$ $p_{IJ}$), `[electro_optic]` ($6\times 3$ $r_{IJ}$ in $\text{pm/V}$), and `[viscoelastic]` ($6\times 6$ $\eta_{IJ}$ in $\text{Pa}\cdot\text{s}$, validated positive semi-definite).
* **Bundled library (`materials/*.toml`, 40 files) and importers (`import.rs`):**
  * Upgraded all 38 existing `materials/*.toml` files to `format = 2` and added `nbn_sc.toml` (NbN thin film, $T_c = 16.0\text{ K}$) and `beta_tin_sc.toml` ($\beta$-Sn, $T_c = 3.722\text{ K}$); kept `materials/examples/silicon_measured.toml` at `format = 1` to test legacy loading continuously.
  * Updated `import.rs` to emit `format = 2` (`rii_formula_4`, separate `k_data`, MP `[piezoelectric]` + `[dielectric]`, NumBAT `[photoelastic]` + `[viscoelastic]`, and AFLOW/MP `point_group`).
* **Solver, engine, MCP, CLI, and GUI integration:**
  * `femcore::props` adds `TensorKind::PermittivityIm` and `Materials::cell_viscosity`; `nyx-physics::decode` evaluates T-dependent `[thermal]`, complex optical/microwave permittivity, and `[viscoelastic]` viscosity at study `(frequency, temperature)`.
  * `nyx-ems` ports `surface_impedance` (`ibc.rs`) onto 2D waveguide (`modes`, `modes_complex`), axisymmetric (`axisym`), and 3D (`cavity`) Maxwell solvers alongside complex dielectric permittivity $\tilde{\varepsilon}_r$.
  * `nyx-sms` routes eigenfrequency solves with `[viscoelastic]` damping (`f_search > 0`) through the quadratic complex eigensolver (`lrb::quadratic_modes`) with $\mathbf{C}_{\text{visc}} = \int \mathbf{B}^T \boldsymbol{\eta} \mathbf{B}\,d\Omega$.
  * `nyx-engine` adds `evaluate_material`, `bcs_surface_impedance`, `search_remote_materials`, and `fetch_material.api_key` (passing `MP_API_KEY` via `curl --config -` stdin without leaking to `argv` or project files); `nyx material` adds `eval`, `bcs`, `search`, and `--api-key`; the GUI Materials tab runs remote search/fetch asynchronously over `mpsc` and provides a live $(T, f, d)$ evaluator.

## Consequences

* Microwave resonator $Q_i$, superconducting kinetic-inductance frequency shifts, and Akhiezer mechanical $Q_m$ are computed directly from the material TOML library across temperature and frequency sweeps.
* All `format = 1` files and existing `.v1.json` contracts continue to load and pass unmodified.

## Code

`crates/femcore/src/material/{mod,spec,schema,units,reader,optics,import,bcs,coupled,state_t}.rs`, `crates/femcore/src/{props,optical}.rs`, `crates/femcore/tests/{material_bcs,material_contract,material_extends,material_import,material_library,materials}.rs`, `crates/nyx-physics/src/{physics,decode}.rs`, `crates/nyx-ems/src/{ibc,problem,lib,modes,modes_complex,axisym,cavity}.rs`, `crates/nyx-sms/src/{lrb,solid}.rs`, `crates/nyx-engine/src/{materials,material_fetch,commands}.rs`, `crates/app/src/{matcmd,state,tabs}.rs`, `materials/*.toml`, `docs/MATERIALS.md`, `docs/architecture/contracts/*v2*`.
