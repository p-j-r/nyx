# `update_feature`

**Change a feature** · access `edit`: changes the project

Changes some fields of a feature; give only the fields to change, in the form the add command takes.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `fields` | object | yes |  | fields to change |
| `name` | text | yes |  | feature name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"update_feature","args":{"fields":{},"name":"<name>"}}
```
