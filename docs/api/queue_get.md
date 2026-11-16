# `queue_get`

**One job of the queue** · access `read`: reads the project or the results

The entry of `queue_list` for one job and the last lines of its stage log.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | yes |  | job id |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_get","args":{"id":0}}
```
