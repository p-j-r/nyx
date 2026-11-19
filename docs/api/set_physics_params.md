# `set_physics_params`

**Set the options of a physics** · access `edit`: changes the project

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `params` | object | yes |  | physics options |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_physics_params","args":{"params":{},"physics":"dfs"}}
```
