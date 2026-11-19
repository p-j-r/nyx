# `result_field`

**Read a field** · access `read`: reads the project or the results

The values of one field of one solution (sampled evenly when there are more than `max_values`), with `stale` and `age_s` of the result it was read from.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `field` | text | yes |  | field name |
| `max_values` | integer | no |  | default 2000 |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `solution` | integer | no |  | solution number from 0, default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"result_field","args":{"field":"<field>","physics":"dfs"}}
```
