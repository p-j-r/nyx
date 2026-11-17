# `result_info`

**Describe the results** · access `read`: reads the project or the results

For each physics solved: the study, size, time, solutions with their values, fields (node-averaged, real) and `cell_fields` (the per-cell fields the reductions read: name, unit, complex or not) and scalars, `stale` (whether what this result was made from has changed since: its own problem, the materials, the mesh or a physics it takes fields from) and `age_s` (seconds since it was stored).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"result_info","args":{}}
```
