# `queue_get`

**One job of the queue** · access `read`: reads the project or the results

The entry of `queue_list` for one job and the last lines of its stage log. Give `id` (a job of this session) or `shared_id` (an item of the workspace queue).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `id` | integer | no |  | job id |
| `shared_id` | text | no |  | an item of the workspace queue (from queue_list), instead of `id`: the queue is shared, so it can be an item another program wrote |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_get","args":{}}
```
