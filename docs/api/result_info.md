# `result_info`

**Describe the results** · access `read`: reads the project or the results

For each physics solved: the study, size, time, solutions with their values, fields (node-averaged, real) and `cell_fields` (the per-cell fields the reductions read: name, unit, complex or not), scalars, unconditional solver `diagnostics` and numerical `error_budget`, `stale` (whether what this result was made from has changed since: its own problem, the materials, the mesh or a physics it takes fields from) and `age_s` (seconds since it was stored), and `validated`: `{state, tier, reason, check?}` where `state` is `exploratory` (no mesh-convergence check with a tolerance has passed on this problem), `checked` (one did or flux bounds certified it) or `check_not_affordable` (not checked, and a check is expected not to fit this machine). A stale result is never checked; every answer that carries `stale` carries `validated`.

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
