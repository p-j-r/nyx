# `analyse`

**Analyse a result** · access `read`: reads the project or the results

Runs an analysis on the latest result of a physics: polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, or network. For `capacitance_matrix`, optional `params` perform ground projection or floating Schur elimination (`ground`, `eliminate: "drop" | "floating_schur"`), grounded/differential lumped qubit reduction (`qubits: [{ name, pads, E_J_Hz? }]`), resonator ZPF / dispersive coupling extraction (`resonators: [{ name, terminal, f_r_Hz, Z0_ohm?, resonator_type?, C_r_F? }]`), and matched-asymptotic finite-metal-thickness correction (`metal_thickness_m`, `edge_cutoff_m`). For `hamiltonian`, outputs include `f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, and `t1_purcell_s`, and `params` accept `sweep_port` (0-based index or port name), `squid_asymmetry`, `squid_flux_sweep`, and `ej_hz_sweep`. For `tls_participation`, interface thicknesses/permittivities/loss tangents auto-load from the material TOML `[interfaces.{ma,ms,sa}]` and `params` accept `trench_depth_m` (substrate over-etch depth in metres). For `network`, extracts damped companion-QEP modes, ROM-synthesized Z/Y/S matrices, Schur-complement Black-Box Quantization (`bbq_modes`), Taubin algebraic circle-fit resonator parameters (`s21_fit`), and Touchstone v1.1 (`.sNp`) text from either `CircuitRom` or `DrivenPortData`. The reply carries `stale`, `age_s`, and `validated` of the result it was run on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `analysis` | text | yes |  | analysis id (list_analyses) |
| `params` | object | no |  | optional analysis parameters (capacitance_matrix: ground, eliminate, qubits, resonators, metal_thickness_m, edge_cutoff_m; hamiltonian: basis, fock_trunc, sweep_port, flux_sweep, squid_asymmetry, squid_flux_sweep, ej_hz_sweep; tls_participation: trench_depth_m, t_ma/t_ms/t_sa, etc.; network: freqs_hz, terminations_ohm, active_ports, s21_ports, resonator_branch, touchstone_format; epr options) |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"analyse","args":{"analysis":"<analysis>","physics":"dfs"}}
```
