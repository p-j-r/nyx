# `remove_param`

**Remove a model parameter** · access `edit`: changes the project

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | parameter name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"remove_param","args":{"name":"<name>"}}
```
