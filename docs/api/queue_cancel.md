# `queue_cancel`

**Cancel a queued or running job** · access `expensive`: starts a mesh build or a solve, or writes files

A queued job never runs. A running one stops at its next check (a factorisation in progress finishes first; a solve in a worker process is stopped at once); `cancel_requested` says so until it stops. Give `id` or `shared_id`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | no |  | job id |
| `shared_id` | text | no |  | an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_cancel","args":{}}
```
