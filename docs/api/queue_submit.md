# `queue_submit`

**Submit work to the queue** · access `expensive`: starts a mesh build or a solve, or writes files

Queues a `solve` (default; `physics` is a code or a list), a `mesh` build or a mesh `convergence` check (`physics`, `factor`, `grids`, `tolerance`) and returns its entry at once. The heavy lane runs one job at a time in order of priority (1 a person's request, 2 an agent's, the default, 3 sweeps and batches; a waiting job rises one level after a minute), then in turns between owners, then in order of submission. The job carries a frozen copy of the project: later edits do not change it (`project_changed` says they happened). `supersede` replaces a queued request of the same kind and targets; `hold` keeps it from starting until `queue_release`. A job that comes after another of the same project reuses the mesh the first built.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `apply_extrapolated` | true or false | no |  | convergence: promote extrapolated capacitance matrix onto stored result |
| `coarsen` | true or false | no |  | convergence: invert factor < 1 into a coarsening ladder relative to present mesh |
| `detached` | true or false | no |  | solve only, with `shared`: write the item into the queue and return at once, with nothing of this session waiting for it; the answer is the item (`shared_id`, and `results`, the file that `open_results` reads once its state is `done`). Closing the program leaves it running |
| `extrapolate` | text: `richardson_log_spd`, `richardson_linear` | no |  | convergence: optional Richardson extrapolation mode |
| `factor` | number | no |  | convergence: comparison/present element size, default 0.6 |
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `grids` | integer: `2`, `3` | no |  | convergence: number of mesh levels (2 default, or 3 for Roache GCI) |
| `hold` | true or false | no |  | start held |
| `isolate` | true or false | no |  | solve in a worker process (only for kind solve): running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do |
| `kind` | text: `solve`, `mesh`, `convergence` | no |  | default solve |
| `match` | text: `index`, `nearest` | no |  | convergence: solution pairing strategy (index or nearest) |
| `memory_wait_timeout_s` | number | no |  | solve only: maximum seconds (>= 0) to wait when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default) |
| `physics` | text or list of text | no |  | physics code(s); default the active one |
| `post_analyses` | list of object | no |  | read-only post-solve analyses or field reductions run on full in-memory fields before fields are dropped (either { "analysis": "...", "physics"?: "...", "params"?: {...} } shorthand or { "cmd": "analyse" \| "integrate" \| "confinement" \| "extremum" \| "probe" \| "line_cut" \| "overlap", "args": {...} }). each item: one read-only post-solve analysis or field-reduction request |
| `post_analyses[].analysis` | text | no |  | shorthand analysis id to run via `analyse` (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator, uq_cqed) |
| `post_analyses[].args` | object | no |  | arguments object passed to `cmd` |
| `post_analyses[].cmd` | text | no |  | explicit read-only post-solve command name (analyse, integrate, confinement, extremum, probe, line_cut, overlap, result_info, result_field, list_analyses, field, field_op) |
| `post_analyses[].label` | text | no |  | optional label attached to this post-analysis entry in the reply |
| `post_analyses[].params` | object | no |  | optional parameters object for shorthand `analysis` |
| `post_analyses[].physics` | text | no |  | three-letter physics code for shorthand `analysis`; defaults to the solved physics that supports the analysis |
| `priority` | integer | no |  | 1, 2 or 3; lower runs first |
| `quantities` | list of text | no |  | convergence: quantity filter list. each item: quantity filter |
| `retain_amr_data` | true or false | no |  | solve only: retain transient cellwise ElectroAmrData indicators on the stored result instead of stripping them after solve; default false |
| `save_fields` | true or false | no |  | solve only: when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results after running any post_analyses; default true |
| `shared` | true or false | no |  | hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process |
| `strict_material_range` | true or false | no |  | solve only: refuse out-of-range material temperature evaluations with an error instead of clamping and emitting a warning (default false; overrides project strict_material_range) |
| `supersede` | true or false | no |  | replace a queued request of the same kind and targets |
| `tolerance` | number | no |  | convergence: largest accepted relative difference |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_submit","args":{}}
```
