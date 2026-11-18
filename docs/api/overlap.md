# `overlap`

**Overlap of two fields** · access `read`: reads the project or the results

The overlap integral of two fields (or two lists of components, e.g. [Ex, Ey, Ez]) over a region: the complex integral of sum_k conj(a_k) b_k, the two norms (integrals of sum_k |c_k|^2) and the coefficient |overlap|^2 / (norm_a norm_b), a number from 0 to 1.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `a` | text or number or object or list of any | yes |  | a field, a formula, an expression, or a list of components |
| `b` | text or number or object or list of any | yes |  | a field, a formula, an expression, or a list of components |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"overlap","args":{"a":"<a>","b":"<b>"}}
```
