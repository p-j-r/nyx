# `queue_set_priority`

**Change the priority of a queued job** · access `expensive`: starts a mesh build or a solve, or writes files

1 runs first, 3 last.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |
| `priority` | integer | yes |  | 1, 2 or 3 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_set_priority","args":{"id":0,"priority":0}}
```
