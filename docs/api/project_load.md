# `project_load`

**Load a project from JSON** · access `edit`: changes the project

Replaces the project by JSON text, after checking it.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `json` | text | yes |  | project JSON text |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"project_load","args":{"json":"<json>"}}
```
