# `describe_material`

**Describe a material** · access `read`: reads the project or the results

What a material gives (elastic class, optical model, thermal, electrical) and its TOML text.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | material name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"describe_material","args":{"name":"<name>"}}
```
