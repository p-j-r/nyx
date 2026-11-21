# `bcs_surface_impedance`

**Evaluate BCS superconducting complex conductivity and sheet impedance** · access `read`: reads the project or the results

Evaluates Zimmermann arbitrary-purity BCS complex conductivity sigma = sigma_1 - j sigma_2 and thin-film sheet surface impedance Z_s = R_s + j omega L_k_sq (with Dynes pair-breaking Gamma, residual R_res and non-equilibrium x_qp_floor) for a library/project material or explicit superconducting parameters.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `dynes_gamma` | number | no |  | Dynes pair-breaking broadening Gamma in J |
| `ell` | number | no |  | electron mean free length ell in m when material is omitted |
| `film_model` | text: `kerr`, `two_sided`, `thin_limit` | no |  | thin-film sheet impedance model: kerr, two_sided or thin_limit |
| `frequency` | number | yes |  | microwave frequency in Hz |
| `gap_0` | number | no |  | zero-temperature BCS gap Delta(0) in J when material is omitted |
| `lambda_L` | number | no |  | clean-limit London penetration depth in m when material is omitted |
| `material` | text | no |  | optional material name with a [superconducting] section |
| `r_res` | number | no |  | residual sheet resistance floor in Ohm/sq |
| `sigma_n` | number | no |  | normal-state conductivity in S/m when material is omitted |
| `tc` | number | no |  | critical temperature Tc in K when material is omitted |
| `temperature` | number | yes |  | operating temperature in K |
| `thickness` | number | no |  | film thickness in m (default: material nominal_thickness or 100e-9) |
| `x_qp_floor` | number | no |  | non-equilibrium normalized quasiparticle density floor |
| `xi_0` | number | no |  | BCS/Pippard coherence length in m when material is omitted |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"bcs_surface_impedance","args":{"frequency":0.0,"temperature":0.0}}
```
