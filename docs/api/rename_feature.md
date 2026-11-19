# `rename_feature`

**Rename a feature** · access `edit`: changes the project

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | feature name |
| `new_name` | text | yes |  | new name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"rename_feature","args":{"name":"<name>","new_name":"<new_name>"}}
```
