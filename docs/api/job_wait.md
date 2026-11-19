# `job_wait`

**Wait for a job** · access `expensive`: starts a mesh build or a solve, or writes files

Blocks until the job has finished and applies its outcome. With `timeout_s` it gives up after that many seconds and returns the job's entry with `waiting: true` (the job goes on; wait again).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |
| `timeout_s` | number | no |  | longest wait in seconds; omit to wait until the job has finished |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"job_wait","args":{"id":0}}
```
