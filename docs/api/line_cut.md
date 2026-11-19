# `line_cut`

**Values of an expression along a line** · access `read`: reads the project or the results

The value of an expression at `n` equally spaced points from `from` to `to` (inclusive), with the distance `s` along the line. Points outside the mesh answer `inside: false`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `from` | list of number | yes |  | start point: [x, z] on a 2D mesh, [x, y, z] on a 3D one |
| `n` | integer | no |  | number of points, default 101 |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `solution` | integer | no |  | default 0 |
| `to` | list of number | yes |  | end point |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"line_cut","args":{"expr":"<expr>","from":[],"to":[]}}
```
