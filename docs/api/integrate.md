# `integrate`

**Integrate an expression over a region** · access `read`: reads the project or the results

The integral of an expression of fields over a region (a region name, comma list or `name*`; default the whole mesh), with the degree-4 rule of the finite elements (exact for the product of two P2 fields). A real expression answers a number, a complex one [re, im]. `measure` is the area or volume of the region. `revolution: true` weights by 2 pi r for a body of revolution (a 2D r-z mesh). Answers `stale` and `age_s` of what it read.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a handle, a field name, a number or {op, args} |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"integrate","args":{"expr":"<expr>"}}
```
