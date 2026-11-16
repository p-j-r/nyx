# `project_open`

**Open a project file** · access `edit`: changes the project

Reads a project JSON file.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `path` | text | yes |  | file path |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"project_open","args":{"path":"<path>"}}
```
