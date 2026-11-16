# `sweep_cancel`

**Cancel the points of a sweep that have not finished** · access `expensive`: starts a mesh build or a solve, or writes files

Queued points never run; a running one stops at its next check.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `sweep` | integer | yes |  | the sweep's number |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"sweep_cancel","args":{"sweep":0}}
```
