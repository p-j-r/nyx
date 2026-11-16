# `queue_cancel`

**Cancel a queued or running job** · access `expensive`: starts a mesh build or a solve, or writes files

A queued job never runs. A running one stops at its next check (a factorisation in progress finishes first); `cancel_requested` says so until it stops.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_cancel","args":{"id":0}}
```
