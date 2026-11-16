# `remove_condition`

**Remove a condition** · access `edit`: changes the project

By its position (from 0) in the physics' list, as project_info shows.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `index` | integer | yes |  | position from 0 |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"remove_condition","args":{"index":0,"physics":"dfs"}}
```
