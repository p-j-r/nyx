# 118. cQED open Robin BC, cavity lumped ports, driven PROM/waveports, unified `NetworkAnalysis`, DC-SQUID sweeps, 3D metal-thickness capacitance, and TOML TLS

Date: 2026-11-22

## Status

Accepted. Amends ADR 103, 104, 105, 106, 109, and 115.

## Context

Full SQuADDS superconducting chip characterization requires: (1) accurate far-field electrostatics without artificial ground walls and finite-thickness correction on zero-thickness 3D sheets; (2) fast, accurate multi-port driven PROM sweeps with arbitrary frequency lists, dielectric loss tangents, CPW even-mode lumped ports, contact validation, and modal waveport targeting; (3) passive resistive `lumped_port` coupling and lossless kinetic-inductance IBC inside 3D cavities; (4) unified QEP / Schur-BBQ / Taubin $S_{21}$ circle-fit / Touchstone `.sNp` network extraction alongside cavity $f_{\text{cav}}$, $g$, $\kappa$, Purcell $T_1$, multi-port `sweep_port`, and asymmetric DC-SQUID sweeps; and (5) BVH-accelerated TLS participation auto-loaded from material TOML interfaces with continuous substrate trench depth.

## Decision

* **Electrostatics (`nyx-dfs`):**
  * Added 3D spherical $1/r$ asymptotic Robin boundary condition `kind: "open"` (optional `origin`, default mesh bbox centroid) allowing isolated conductors to solve without ground, plus pre-symmetrization Frobenius asymmetry scalar `maxwell_raw_asymmetry_rel`.
  * Populated `TlsMeta::interfaces` automatically from material TOML `[interfaces.{ma,ms,sa}]` (`nyx_physics::decode::tls_interfaces_from_spec`, including `materials/silicon_4k.toml`).
* **Cavity & driven EM (`nyx-ems`):**
  * Driven study: added `max_hdm_solves`, `strict_prom_tol`, `seed_freqs_hz`, ascending `freqs_hz`, unit-$L_2$ Hermite snapshot normalization, Block-SVD anchor compression, fixed-union CSR pattern with persistent `DirectSolver::update_matrix` (`symbolic_count == 1`), deterministic parallel reduced sweep, complex permittivity via `eps_r_im`/`loss_tangent`, `waveport` `n_eff_guess`/`z0_ohm`, `lumped_port` `mode: "single" | "cpw_even"`, and pre-solve 3D port contact/PEC-short validation.
  * Cavity study: enabled `lumped_port` with optional `complex_damping` (real eigensolve when `false`, `ruiz_unscale_and_rayleigh_damped` when `true`), appended resistive ports into `CircuitRom` (projected strictly from non-PML physical matrices) with `rom_complex_eigenfrequencies`, added `lossless: true` and material-TOML `thickness`/`film_model` fallbacks on `surface_impedance`, added `squid_asymmetry` on `lumped_jj`, fixed complex equipartition scaling `Re(k_0^2)`, and excluded `jj_labels` from `TlsMeta::sheet_labels`.
* **Post-solve cQED analyses (`nyx-analysis` & `nyx-engine::analyse`):**
  * Registered `NetworkAnalysis` (`id: "network"`, `crates/nyx-analysis/src/network.rs`) operating on `CircuitRom` or `DrivenPortData` for $2r\times 2r$ companion QEP, ROM `Z/Y/S` synthesis, dual-Woodbury Schur BBQ (`bbq_modes`), Taubin algebraic circle + arctan phase fit (`s21_fit`), and Touchstone v1.1 (`.sNp`).
  * Extended `HamiltonianAnalysis` with linear Craig–Bampton `f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, `t1_purcell_s` (with Loewner PSD damping deduplication), passive `lumped_port` filtering, `sweep_port` (index or name), and native DC-SQUID `squid_effective_ej`/`squid_asymmetry`/`squid_flux_sweep`/`ej_hz_sweep`.
  * Added 3D matched-asymptotic `I_edge` finite-metal-thickness correction (`metal_thickness_m`, `edge_cutoff_m`) to `CapacitanceMatrixAnalysis` (`< 0.3%` vs 2D finite-thickness CPW, rejecting `planar_bem`), and upgraded `TlsParticipationAnalysis` with flat `SegmentBvh3`, `TlsMeta::interfaces` TOML defaults, and continuous `trench_depth_m`.

## Consequences

* Every cQED simulation stage from open-box finite-thickness capacitance and lossless BCS cavity modes to Purcell $T_1$, asymmetric SQUID flux tuning, BBQ, Taubin $S_{21}$ fitting, and Touchstone export is accessible through uniform `Engine::call` commands with verified contract descriptors.

## Code

`materials/silicon_4k.toml`, `crates/nyx-physics/src/{cqed,decode,lib}.rs`, `crates/nyx-dfs/src/{lib,problem,electro,planar_bem}.rs`, `crates/nyx-ems/src/{lib,problem,driven,ports,ibc,cavity,rom}.rs`, `crates/nyx-analysis/src/{lib,network,epr,hamiltonian,capacitance,tls}.rs`, `crates/nyx-analysis/tests/wave1_cqed_analyses.rs`, `crates/nyx-engine/src/{analyse,commands}.rs`, `crates/nyx-engine/cards/cqed.md`, `crates/nyx-engine/tests/squadds_fr_e2e.rs`, `docs/architecture/contracts/{physics_dfs,physics_ems,engine_commands}.v1.json`.
