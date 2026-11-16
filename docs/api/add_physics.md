# `add_physics`

**Add a physics to the project** · access `edit`: changes the project

Adds a problem for a physics (a project loads one to three). The study defaults to the first the physics lists.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `order` | integer | no |  | element order, default 1 |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `study` | text | no |  | study kind (see describe_physics) |
| `study_params` | object | no |  | study parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_physics","args":{"physics":"dfs"}}
```
