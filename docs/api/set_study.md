# `set_study`

**Set the study of a physics** · access `edit`: changes the project

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `kind` | text | yes |  | study kind |
| `order` | integer | no |  | element order |
| `params` | object | no |  | study parameters |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_study","args":{"kind":"<kind>","physics":"dfs"}}
```
