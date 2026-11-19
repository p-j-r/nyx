# `queue_list`

**List the queue** · access `read`: reads the project or the results

Every job of this session and every item that another program wrote into the workspace queue (`foreign`, with a `shared_id`), oldest first, with its state (`queued`, `running`, `done`, `cancelled`, `failed`), lane, owner, priority, position in the queue, the time it waited and ran, its stage, what it was asked (`kind`, `targets`), `project_changed` and `applied`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `owner` | text | no |  | only this owner: gui, ai or cli |
| `state` | text | no |  | only this state |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_list","args":{}}
```
