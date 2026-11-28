# `analyse`

**Analyse a result** · access `read`: reads the project or the results

Runs an analysis on the latest result of a physics: polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, uq_cqed, far_field, beam_overlap, or optomechanics. For `capacitance_matrix`, optional `params` perform ground projection or floating Schur elimination (`ground`, `eliminate: "drop" | "floating_schur"`), grounded/differential lumped qubit reduction (`qubits: [{ name, pads, E_J_Hz? }]`), resonator ZPF / dispersive coupling extraction (`resonators: [{ name, terminal, f_r_Hz, Z0_ohm?, resonator_type?, C_r_F? }]`), and matched-asymptotic finite-metal-thickness correction (`metal_thickness_m`, `edge_cutoff_m`). For `hamiltonian`, outputs include `f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, and `t1_purcell_s`, and `params` accept `sweep_port` (0-based index or port name), `squid_asymmetry`, `squid_flux_sweep`, and `ej_hz_sweep`. For `tls_participation`, interface thicknesses/permittivities/loss tangents auto-load from the material TOML `[interfaces.{ma,ms,sa}]`, `params` accept `trench_depth_m` (substrate over-etch depth in metres), and passing `electrostatic_physics` (e.g. `"dfs"`, with optional `port_terminals`, `modal_voltages`, `modal_freqs_hz`) evaluates hybrid quasi-electrostatic modal TLS participation and `T1_TLS_s` on the fine electrostatic mesh. For `network`, extracts damped companion-QEP modes, ROM-synthesized Z/Y/S matrices, Schur-complement Black-Box Quantization (`bbq_modes`), Taubin algebraic circle-fit resonator parameters (`s21_fit`), and Touchstone v1.1 (`.sNp`) text from either `CircuitRom` or `DrivenPortData`. For `piezo_resonator`, extracts series/parallel resonance frequencies (`f_s_Hz`, `f_p_Hz`), IEEE coupling factors (`k_t2`, `k_eff2`), unwrapped phase-slope quality factors (`Q_s`, `Q_p`), and bounded 6-parameter mBVD equivalent circuit (`C_0`, `R_0`, `R_s`, `C_m`, `L_m`, `R_m`) from `DrivenPortData`. For `uq_cqed`, propagates fabrication and constitutive 1-sigma `uncertainties` (plus `seed`, `n_samples`, `electrostatic_physics` companion bounds) through zero-solve EPR/Hamiltonian/Capacitance/TLS post-processing while keeping ASME V&V 20 epistemic brackets strictly separated from aleatoric variances. For `optomechanics`, evaluates photoelastic (`g_pe_hz`), moving-boundary (`g_mb_hz`), total vacuum optomechanical coupling (`g_0_hz`, `m_eff_kg`, `x_zpf_m`), and stimulated Brillouin scattering gain (`sbs_gain_w_inv_m_inv`) across shared or independent optical (`optical_physics`) and mechanical (`mechanical_physics`) meshes (ADR 141). The reply carries `stale`, `age_s`, and `validated` of the result it was run on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `analysis` | text | yes |  | analysis id (list_analyses) |
| `params` | object | no |  | optional analysis parameters (capacitance_matrix: ground, eliminate, qubits, resonators, metal_thickness_m, edge_cutoff_m; hamiltonian: basis, fock_trunc, sweep_port, flux_sweep, squid_asymmetry, squid_flux_sweep, ej_hz_sweep; tls_participation: trench_depth_m, t_ma/t_ms/t_sa, electrostatic_physics, port_terminals, modal_voltages, modal_freqs_hz; network: freqs_hz, terminations_ohm, active_ports, s21_ports, resonator_branch, touchstone_format; piezo_resonator: port, terminations, f_window_hz, z_ref_ohm, fix_r_s_ohm; uq_cqed: seed, n_samples, basis, uncertainties, electrostatic_physics, qubits, resonators, metal_thickness_m; optomechanics: optical_physics, mechanical_physics, optical_mode, mechanical_mode, core_domains, interface_boundaries, acoustic_q; far_field/beam_overlap/epr options) |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"analyse","args":{"analysis":"<analysis>","physics":"dfs"}}
```
