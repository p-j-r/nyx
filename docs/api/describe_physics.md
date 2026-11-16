# `describe_physics`

**Describe a physics** · access `read`: reads the project or the results

Studies, conditions, properties, fields and couplings of one physics, and the JSON Schema of its problem.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"describe_physics","args":{"physics":"dfs"}}
```
