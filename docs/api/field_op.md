# `field_op`

**Arithmetic on handles** · access `read`: reads the project or the results

Combines handles: `add`, `mul` (two or more operands), `sub`, `div`, `neg`, `abs`, `abs2` (|z|^2 with no square root), `re`, `im`, `conj`, `sqrt` and `pow` (with `n`, a whole number 1 to 16). An operand is a handle, a field name of `physics`/`solution`, a number, or another operation `{op, args}`. The operation is checked against the results and returned as one handle.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `args` | list of any | yes |  | the operands |
| `n` | integer | no |  | the power, for pow |
| `op` | text | yes |  | add, sub, mul, div, neg, abs, abs2, re, im, conj, sqrt or pow |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `solution` | integer | no |  | default 0 |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"field_op","args":{"args":[],"op":"<op>"}}
```
