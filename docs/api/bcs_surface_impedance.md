# `bcs_surface_impedance`

**BCS thin-film surface impedance** · access `read`: reads the project or the results

Evaluates Zimmermann arbitrary-purity BCS complex conductivity and thin-film surface impedance (Z_s = R_s + i omega L_k_sq) for a superconducting material or parameter set.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `dynes_gamma` | number | no |  | Dynes pair-breaking broadening Gamma in J |
| `film_model` | text: `two_sided`, `kerr`, `thin_limit` | no | `kerr` | thin-film boundary convention: two_sided (interior sheet jump), kerr (one-sided backed wall) or thin_limit (Pearl d << lambda_eff) |
| `freq` | number | yes |  | frequency in Hz |
| `gap_0` | number | no |  | zero-temperature energy gap Delta(0) in J |
| `gap_0_eV` | number | no |  | zero-temperature energy gap Delta(0) in eV |
| `lambda_L` | number | no |  | clean zero-temperature London penetration depth in metres |
| `material` | text | no |  | name of a library or project material with a [superconducting] block |
| `mean_free_path` | number | no |  | electron mean free path in metres |
| `r_res` | number | no |  | residual sheet resistance floor in Ohm/sq |
| `sigma_n` | number | no |  | normal-state conductivity in S/m |
| `tc` | number | no |  | critical temperature T_c in K (when material is omitted or overridden) |
| `temperature` | number | no |  | temperature in K, default 0.01 |
| `thickness` | number | yes |  | film thickness in metres |
| `x_qp_floor` | number | no |  | non-equilibrium normalized quasiparticle density floor |
| `xi_0` | number | no |  | BCS coherence length in metres |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"bcs_surface_impedance","args":{"freq":1.0,"thickness":1e-12}}
```
