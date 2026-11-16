# `sweep_run`

**Run a sweep: one recorded run per point** · access `expensive`: starts a mesh build or a solve, or writes files

Queues one solve per point and returns the sweep's number at once. A point is a list of edits (`{ "cmd": "set_param", "args": { "name": "w", "value": 4e-7 } }`: any command that changes the project except those that replace it or start work), applied to a copy of the project: what varies, a geometry parameter, the wavelength in the study, a material, a mesh setting, is the case's to decide. Give `points` (each `{ label, edits }`) or `series` (`{ cmd, args, field, values }`: the command repeated with `field` set to each value). `physics` is the code or list solved at every point (default: what a solve would solve). The points are queued at priority 3 (sweeps and batches) so a person's request and a single request of an agent go before them, and run on the same mesh and solve caches as any request. Every point's edits are applied and checked first: a bad point fails the request with its number and nothing runs. The project, the mesh and the results of the engine are not touched; each point is recorded in the run table when the project is bound to a workflow. Read the table with `sweep_table`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `hold` | true or false | no |  | queue the points held; release each with queue_release |
| `isolate` | true or false | no |  | solve every point in a worker process, as for solve_run |
| `name` | text | no |  | a name for the table |
| `physics` | text or list of text | no |  | code or list of codes to solve at every point |
| `points` | list of object | no |  | the points, each a list of edits |
| `points[].edits` | list of object | yes |  | the edits that make this point, applied in order to a copy of the project |
| `points[].edits[].args` | object | no |  | its arguments |
| `points[].edits[].cmd` | text | yes |  | a command that changes the project |
| `points[].label` | text | no |  | a name for the point in the table; default `point N` |
| `priority` | integer | no |  | 1, 2 or 3; default 3 |
| `series` | object | no |  | one command repeated with `field` set to each value |
| `series.args` | object | no |  | its arguments that stay the same |
| `series.cmd` | text | yes |  | the edit command repeated at every point |
| `series.field` | text | yes |  | the argument that takes each value |
| `series.values` | list of any | yes |  | the values, one point each |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"sweep_run","args":{}}
```
