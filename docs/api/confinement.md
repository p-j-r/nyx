# `confinement`

**Fraction of a field in a region** · access `read`: reads the project or the results

The integral of sum_k |c_k|^2 over a region divided by the same over the whole mesh, for the components `fields` (e.g. [Ex, Ey, Ez]). It has no material weight (no permittivity): for the energy or power fraction use an expression that has it.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `fields` | text or number or object or list of any | yes |  | a field, an expression, or a list of components |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | yes |  | the region the fraction is of |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"confinement","args":{"fields":"<fields>","region":"<region>"}}
```
