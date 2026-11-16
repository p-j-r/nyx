# `set_material`

**Assign a material to a region** · access `edit`: changes the project

A material by name (library or project: list_materials), a named material placed at an angle (`{"material": "sapphire", "euler_deg": [0, 57.6, 0]}`, Bunge angles in degrees), or scalar properties by descriptor name (`{"E": 1.7e11, "nu": 0.28}`).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `material` | text or object | yes |  | material name or properties |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `region` | text | yes |  | region name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_material","args":{"material":"<material>","physics":"dfs","region":"<region>"}}
```
