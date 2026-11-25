# `sweep_run`

**Run a sweep: one recorded run per point** · access `expensive`: starts a mesh build or a solve, or writes files

Queues one solve per point and returns the sweep's number at once. A point is a list of edits (`{ "cmd": "set_param", "args": { "name": "w", "value": 4e-7 } }`: any command that changes the project except those that replace it or start work), applied to a copy of the project: what varies, a geometry parameter, the wavelength in the study, a material, a mesh setting, is the case's to decide. Give `points` (each `{ label, edits }`) or `series` (`{ cmd, args, field, values }`: the command repeated with `field` set to each value). `physics` is the code or list solved at every point (default: what a solve would solve). The points are queued at priority 3 (sweeps and batches) so a person's request and a single request of an agent go before them, and run on the same mesh and solve caches as any request. Every point's edits are applied and checked first: a bad point fails the request with its number and nothing runs. The project, the mesh and the results of the engine are not touched; each point is recorded in the run table when the project is bound to a workflow. Read the table with `sweep_table`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `hold` | true or false | no |  | queue the points held; release each with queue_release |
| `isolate` | true or false | no |  | solve every point in a worker process, as for solve_run |
| `memory_wait_timeout_s` | number | no |  | maximum seconds (>= 0) to wait per sweep point when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default) |
| `name` | text | no |  | a name for the table |
| `physics` | text or list of text | no |  | code or list of codes to solve at every point |
| `points` | list of object | no |  | the points, each a list of edits |
| `points[].edits` | list of object | yes |  | the edits that make this point, applied in order to a copy of the project |
| `points[].edits[].args` | object | no |  | its arguments |
| `points[].edits[].cmd` | text | yes |  | a command that changes the project |
| `points[].label` | text | no |  | a name for the point in the table; default `point N` |
| `post_analyses` | list of object | no |  | read-only post-solve analyses or field reductions run on full in-memory fields before fields are dropped (either { "analysis": "...", "physics"?: "...", "params"?: {...} } shorthand or { "cmd": "analyse" \| "integrate" \| "confinement" \| "extremum" \| "probe" \| "line_cut" \| "overlap", "args": {...} }). each item: one read-only post-solve analysis or field-reduction request |
| `post_analyses[].analysis` | text | no |  | shorthand analysis id to run via `analyse` (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator) |
| `post_analyses[].args` | object | no |  | arguments object passed to `cmd` |
| `post_analyses[].cmd` | text | no |  | explicit read-only post-solve command name (analyse, integrate, confinement, extremum, probe, line_cut, overlap, result_info, result_field, list_analyses, field, field_op) |
| `post_analyses[].label` | text | no |  | optional label attached to this post-analysis entry in the reply |
| `post_analyses[].params` | object | no |  | optional parameters object for shorthand `analysis` |
| `post_analyses[].physics` | text | no |  | three-letter physics code for shorthand `analysis`; defaults to the solved physics that supports the analysis |
| `priority` | integer | no |  | 1, 2 or 3; default 3 |
| `reuse` | text: `none`, `warm` | no |  | `warm`: each point starts its eigensolves from the modes of the previous point instead of a random block (the first point only records them), which shortens the iterations of a sweep that keeps its mesh (a material, a study value, a boundary condition, a frequency with the mesh not sized by the physics; a point whose mesh or number of unknowns differs starts cold). The numbers agree to the solver's tolerance but not bit for bit with a cold solve, so a warm point is marked `exploratory` in its row and run (`reuse`) until it is solved cold. Default `none`: every point is a cold solve |
| `save_fields` | true or false | no |  | when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results at every sweep point after running any post_analyses; default true |
| `series` | object | no |  | one command repeated with `field` set to each value |
| `series.args` | object | no |  | its arguments that stay the same |
| `series.cmd` | text | yes |  | the edit command repeated at every point |
| `series.field` | text | yes |  | the argument that takes each value |
| `series.values` | list of any | yes |  | the values, one point each |
| `shared` | true or false | no |  | hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"sweep_run","args":{}}
```
