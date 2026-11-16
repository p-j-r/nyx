# `add_material`

**Add a material to the project** · access `edit`: changes the project

Adds (or replaces) a material the project carries, as TOML text (units are written in the file; `extends = "silicon"` starts from a library material). Read at once, so a mistake is reported here.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | material name |
| `toml` | text | yes |  | the material file text |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_material","args":{"name":"<name>","toml":"<toml>"}}
```
