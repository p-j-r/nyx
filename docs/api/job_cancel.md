# `job_cancel`

**Cancel a job** · access `expensive`: starts a mesh build or a solve, or writes files

The solver notices at its next check; a factorisation in progress finishes first.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"job_cancel","args":{"id":0}}
```
