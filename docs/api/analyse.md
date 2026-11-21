# `analyse`

**Analyse a result** · access `read`: reads the project or the results

Runs an analysis on the latest result of a physics: polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, or tls_participation. For `capacitance_matrix`, optional `params` perform ground projection or floating Schur elimination (`ground`, `eliminate: "drop" | "floating_schur"`), grounded/differential lumped qubit reduction (`qubits: [{ name, pads, E_J_Hz? }]`), and resonator ZPF / dispersive coupling extraction (`resonators: [{ name, terminal, f_r_Hz, Z0_ohm?, resonator_type?, C_r_F? }]`). The reply carries `stale`, `age_s`, and `validated` of the result it was run on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `analysis` | text | yes |  | analysis id (list_analyses) |
| `params` | object | no |  | optional analysis parameters (for capacitance_matrix: ground, eliminate, qubits, resonators; also epr, hamiltonian, tls_participation options) |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"analyse","args":{"analysis":"<analysis>","physics":"dfs"}}
```
