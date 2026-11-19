# `job_status`

**Status of a job** · access `read`: reads the project or the results

State, current stage, fraction done and elapsed time.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"job_status","args":{"id":0}}
```
