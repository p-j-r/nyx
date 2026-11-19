# `probe`

**Value of an expression at points** · access `read`: reads the project or the results

The value of an expression at points, interpolated in the cell that holds each point. Points are in mesh axes: [x, z] on a 2D mesh, [x, y, z] on a 3D one. A point outside the mesh answers `inside: false` and a null value.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `points` | list of list of number | yes |  | a list of points, each a list of coordinates |
| `solution` | integer | no |  | default 0 |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"probe","args":{"expr":"<expr>","points":[]}}
```
