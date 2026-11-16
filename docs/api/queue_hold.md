# `queue_hold`

**Hold a queued job** · access `expensive`: starts a mesh build or a solve, or writes files

Keeps a queued job from starting until `queue_release`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_hold","args":{"id":0}}
```
