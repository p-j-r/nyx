# `sweep_table`

**The table of a sweep** · access `read`: reads the project or the results

Applies the points that have finished and returns one row per point: its label and edits, its state (`queued`, `running`, `done`, `cancelled`, `failed`), the numbers of each physics solved (eigenvalues, values, numeric scalars), the unknowns and cells, the seconds and the run it was recorded as. With `wait_s` it waits up to that many seconds for the points to finish (`complete` says whether they all have). Without `sweep` it lists the sweeps of this session.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `sweep` | integer | no |  | the sweep's number |
| `wait_s` | number | no |  | longest wait in seconds, default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"sweep_table","args":{}}
```
