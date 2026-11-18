# `extremum`

**Largest or smallest value of an expression** · access `read`: reads the project or the results

The largest (`op: max`, default) or smallest (`min`) real value of an expression at the nodes of the cells of a region, with the cell and the point. It is the maximum over the nodes, not over the interpolant between them. A complex expression has no order: wrap it in `abs`, `abs2`, `re` or `im`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `op` | text: `max`, `min` | no |  | default max |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `solution` | integer | no |  | default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"extremum","args":{"expr":"<expr>"}}
```
