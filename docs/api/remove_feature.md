# `remove_feature`

**Remove a feature** · access `edit`: changes the project

Removes a feature that no later feature uses.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | yes |  | feature name |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"remove_feature","args":{"name":"<name>"}}
```
