# `add_condition`

**Add a boundary or domain condition** · access `edit`: changes the project

`on` is a selector over the geometry (`outer & plane(x = 0)`), a boundary name or a region name.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `kind` | text | yes |  | condition kind (see describe_physics) |
| `on` | text | yes |  | selector, boundary or region |
| `params` | object | no |  | condition parameters |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_condition","args":{"kind":"<kind>","on":"<on>","physics":"dfs"}}
```
