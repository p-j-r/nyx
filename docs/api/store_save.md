# `store_save`

**Save into a workspace** · access `expensive`: starts a mesh build or a solve, or writes files

Binds the project to a workflow folder (created in `workspace` the first time) or stores a new version of it. Later solves are recorded in the run table.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `workspace` | text | no |  | workspace folder, needed the first time |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"store_save","args":{}}
```
