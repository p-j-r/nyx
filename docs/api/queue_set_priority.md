# `queue_set_priority`

**Change the priority of a queued job** · access `expensive`: starts a mesh build or a solve, or writes files

1 runs first, 3 last. Give `id` or `shared_id`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | no |  | job id |
| `priority` | integer | yes |  | 1, 2 or 3 |
| `shared_id` | text | no |  | an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_set_priority","args":{"priority":0}}
```
