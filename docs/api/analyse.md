# `analyse`

**Analyse a result** · access `read`: reads the project or the results

Runs an analysis on the latest result of a physics: the polarization of the modes gives each mode its class (quasi-TE, quasi-TM, hybrid), Stokes parameters, degree of polarization and the degenerate groups. The reply carries `stale` and `age_s` of the result it was run on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `analysis` | text | yes |  | analysis id (list_analyses) |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"analyse","args":{"analysis":"<analysis>","physics":"dfs"}}
```
