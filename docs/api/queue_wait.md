# `queue_wait`

**Wait for a job, at most a while** · access `expensive`: starts a mesh build or a solve, or writes files

Waits for the job up to `timeout_s` seconds (required: a tool call must not block for ever) and applies its outcome. At the timeout it returns the job's entry with `waiting: true`, state and stage, and the job goes on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |
| `timeout_s` | number | yes |  | longest wait in seconds |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_wait","args":{"id":0,"timeout_s":0.0}}
```
