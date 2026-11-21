# `evaluate_material`

**Evaluate a material at temperature, frequency and wavelength** · access `read`: reads the project or the results

Evaluates a material's temperature-, frequency-, wavelength- and electric-field-dependent properties: thermal curves at T, optical (n, k) and complex permittivity at (wavelength, T), microwave complex permittivity and TLS loss tangent at (frequency, T, e_field), piezoelectric Legendre duals (eps_clamped, eps_free, C^E, C^D), photoelastic, electro-optic, viscoelastic, and BCS superconducting surface impedance.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `e_field` | number | no |  | electric-field amplitude in V/m for saturable TLS dielectric loss |
| `euler_deg` | list of 3 (number) | no |  | Bunge Z-X-Z crystal orientation angles in degrees [phi1, Phi, phi2] |
| `film_model` | text: `kerr`, `two_sided`, `thin_limit` | no |  | thin-film boundary formula for BCS sheet impedance, default kerr |
| `frequency` | number | no |  | evaluation frequency in Hz for microwave dielectric and BCS superconducting models |
| `name` | text | yes |  | material name in the library or active project |
| `temperature` | number | no |  | evaluation temperature in K (default: material reference_temperature) |
| `thickness` | number | no |  | film thickness in m for BCS superconducting sheet impedance |
| `wavelength` | number | no |  | evaluation free-space wavelength in m for optical dispersion models |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"evaluate_material","args":{"name":"<name>"}}
```
