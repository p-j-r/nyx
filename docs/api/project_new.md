# `project_new`

**New project** · access `edit`: changes the project

Starts an empty project (dropping the current one, its mesh and results).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `dim` | integer | no |  | 2 (x-z cross-section, the default) or 3 |
| `h` | number | no |  | global element size in metres |
| `title` | text | no |  | project title |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"project_new","args":{}}
```
