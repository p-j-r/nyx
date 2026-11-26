# `extremum`

**Largest or smallest value of an expression** · access `read`: reads the project or the results

The largest (`op: max`, default) or smallest (`min`) real value of an expression in a region. By default reports the extremum at the nodes while checking P2 candidate cells for larger interior peaks (`interior_extremum`); pass `interior: true` to return the P2 sub-sampled/Newton-refined interior extremum directly as `value`. A complex expression has no order: wrap it in `abs`, `abs2`, `re` or `im`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `interior` | true or false | no |  | for P2 fields, return the sub-sampled and Newton-refined interior extremum directly in `value` instead of the nodal value; default false still checks peak candidates and surfaces `interior_extremum` when higher |
| `op` | text: `max`, `min` | no |  | default max |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `solution` | integer | no |  | default 0 |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"extremum","args":{"expr":"<expr>"}}
```
