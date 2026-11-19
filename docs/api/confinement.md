# `confinement`

**Fraction of a field in a region** · access `read`: reads the project or the results

The integral of sum_k |c_k|^2 over a region divided by the same over the whole mesh, for the components `fields` (e.g. [Ex, Ey, Ez]). Without `weight` it has no material weight; give `weight` a real formula such as "eps_r" for the energy fraction (the integrand is weight * sum_k |c_k|^2, in the region and over the whole mesh).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `fields` | text or number or object or list of any | yes |  | a field, a formula, an expression, or a list of components |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | yes |  | the region the fraction is of |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | default 0 |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |
| `weight` | text | no |  | a real formula that multiplies the integrand at each point, e.g. eps_r; default none |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"confinement","args":{"fields":"<fields>","region":"<region>"}}
```
