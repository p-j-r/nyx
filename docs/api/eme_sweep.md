# `eme_sweep`

**Run an analytical EME sweep** · access `read`: reads the project or the results

Evaluates an analytical post-solve Eigenmode Expansion sweep (`length`, `periodicity`, or `modes`) on the cached `EmeData` of a solved `eme` study without re-solving 2D cross-section eigenmodes. `length` sweeps the propagation length of segment `segment`, `periodicity` sweeps `PeriodicSubgrid` `subgrid` repetitions in O(log N) Redheffer star products, and `modes` sweeps modal truncation counts across cached interface overlaps. The reply carries `stale`, `age_s`, and `validated` of the result it was run on.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `lengths` | list of number | no |  | segment lengths in metres (>= 0) for mode = 'length' |
| `mode` | text: `length`, `periodicity`, `modes` | yes |  | post-solve EME sweep kind: 'length', 'periodicity', or 'modes' |
| `mode_counts` | list of integer | no |  | modal truncation counts (>= 1) per slice for mode = 'modes' |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `reps` | list of integer | no |  | unit-cell repetition counts (>= 1) for mode = 'periodicity' |
| `segment` | integer | no |  | 0-based EME segment index for mode = 'length', default 0 |
| `subgrid` | integer | no |  | 0-based PeriodicSubgrid index for mode = 'periodicity', default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"eme_sweep","args":{"mode":"length"}}
```
