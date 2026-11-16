# `remove_physics`

**Remove a physics from the project** · access `edit`: changes the project

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"remove_physics","args":{"physics":"dfs"}}
```
