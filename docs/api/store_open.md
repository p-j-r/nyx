# `store_open`

**Open a workflow** · access `edit`: changes the project

Opens a workflow folder and makes its experiment the project.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `experiment` | text | no |  | experiment id, default the first |
| `folder` | text | yes |  | workflow folder |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"store_open","args":{"folder":"<folder>"}}
```
