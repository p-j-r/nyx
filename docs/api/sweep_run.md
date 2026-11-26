# `sweep_run`

**Run a sweep: one recorded run per point** · access `expensive`: starts a mesh build or a solve, or writes files

Queues one solve per explicit point (`points` or `series`) or streams a deterministic LHS / Sobol QMC ensemble (`sampling`, checking geometry and `SimplexMesh::check()` per sample and recording invalid samples as `skipped_invalid_geometry` without aborting the ensemble). An explicit point is a list of edits (`{ "cmd": "set_param", "args": { "name": "w", "value": 4e-7 } }`) applied to a copy of the project on demand without upfront `Vec<Engine>` cloning. `physics` is the code or list solved at every point (default: what a solve would solve). Points run at priority 3 on the engine's mesh/solve stage caches; optionally pass `checkpoint_path` (`.h5`) to stream points and skipped metadata incrementally to HDF5 `SweepCheckpoint`. Read the table with `sweep_table`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `checkpoint_lzf` | true or false | no |  | compress checkpoint columns with Shuffle+LZF; default false |
| `checkpoint_path` | text | no |  | optional .h5 file path for incremental HDF5 SweepCheckpoint streaming of solved points and skipped invalid geometries |
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
| `post_analyses[].analysis` | text | no |  | shorthand analysis id to run via `analyse` (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, uq_cqed) |
| `post_analyses[].args` | object | no |  | arguments object passed to `cmd` |
| `post_analyses[].cmd` | text | no |  | explicit read-only post-solve command name (analyse, integrate, confinement, extremum, probe, line_cut, overlap, result_info, result_field, list_analyses, field, field_op) |
| `post_analyses[].label` | text | no |  | optional label attached to this post-analysis entry in the reply |
| `post_analyses[].params` | object | no |  | optional parameters object for shorthand `analysis` |
| `post_analyses[].physics` | text | no |  | three-letter physics code for shorthand `analysis`; defaults to the solved physics that supports the analysis |
| `priority` | integer | no |  | 1, 2 or 3; default 3 |
| `reuse` | text: `none`, `warm` | no |  | `warm`: each point starts its eigensolves from the modes of the previous point instead of a random block (the first point only records them), which shortens the iterations of a sweep that keeps its mesh (a material, a study value, a boundary condition, a frequency with the mesh not sized by the physics; a point whose mesh or number of unknowns differs starts cold). The numbers agree to the solver's tolerance but not bit for bit with a cold solve, so a warm point is marked `exploratory` in its row and run (`reuse`) until it is solved cold. Default `none`: every point is a cold solve |
| `sampling` | object | no |  | deterministic LHS / Sobol QMC streaming design-of-experiments (alternative to `points` or `series`) |
| `sampling.method` | text: `lhs`, `sobol` | no |  | QMC generator (`lhs` SplitMix64 Latin Hypercube Sampling, default; or low-discrepancy `sobol`) |
| `sampling.n_samples` | integer | yes |  | number of QMC samples to stream (1..=4096) |
| `sampling.offset` | integer | no |  | starting index in Sobol sequence (default 1 to avoid the origin corner) |
| `sampling.parameters` | list of object | yes |  | sampled parameter dimensions. each item: one parametric margin (uniform `[min, max]` or normal `{ mean, sigma \| rel_tol }`) |
| `sampling.parameters[].args` | object | no |  | base arguments kept constant (e.g. `{ "name": "w" }`) |
| `sampling.parameters[].bounds` | list of 2 (number) | no |  | alternative `[min, max]` uniform range |
| `sampling.parameters[].cmd` | text | no |  | edit command (default `set_param`) |
| `sampling.parameters[].distribution` | text: `uniform`, `normal` | no |  | marginal distribution (default `uniform`, or `normal` when mean/sigma given) |
| `sampling.parameters[].field` | text | no |  | argument key receiving each sampled value (default `value`) |
| `sampling.parameters[].max` | number | no |  | upper bound (uniform range max, or optional normal clamp max) |
| `sampling.parameters[].mean` | number | no |  | nominal mean for `normal` marginal |
| `sampling.parameters[].min` | number | no |  | lower bound (uniform range min, or optional normal clamp min) |
| `sampling.parameters[].name` | text | no |  | column label in checkpoint and table |
| `sampling.parameters[].rel_tol` | number | no |  | relative 1-sigma tolerance `\|mean\| * rel_tol` for `normal` marginal |
| `sampling.parameters[].sigma` | number | no |  | 1-sigma standard deviation for `normal` marginal |
| `sampling.scramble` | true or false | no |  | apply deterministic SplitMix64 XOR digital scramble to Sobol coordinates; default false |
| `sampling.seed` | integer | no |  | deterministic 64-bit SplitMix64 seed |
| `sampling.skip_invalid_geometry` | true or false | no |  | record geometry, meshing, or SimplexMesh::check failures as `skipped_invalid_geometry` without aborting the ensemble; default true |
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
