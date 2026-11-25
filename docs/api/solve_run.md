# `solve_run`

**Solve** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the active physics (or the named ones, in coupling order), meshing first when needed. With wait false it starts a job. When the project is bound to a workflow the run is recorded.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `isolate` | true or false | no |  | run the solve in a worker process: running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do |
| `memory_wait_timeout_s` | number | no |  | maximum seconds (>= 0) to wait when free RAM is below the solve estimate (overrides NYX_MEMORY_WAIT_TIMEOUT_S and engine default) |
| `physics` | text or list of text | no |  | physics code(s); default the active one |
| `post_analyses` | list of object | no |  | read-only post-solve analyses or field reductions run on full in-memory fields before fields are dropped (either { "analysis": "...", "physics"?: "...", "params"?: {...} } shorthand or { "cmd": "analyse" \| "integrate" \| "confinement" \| "extremum" \| "probe" \| "line_cut" \| "overlap", "args": {...} }). each item: one read-only post-solve analysis or field-reduction request |
| `post_analyses[].analysis` | text | no |  | shorthand analysis id to run via `analyse` (polarization, capacitance_matrix, inductance_matrix, epr, hamiltonian, tls_participation, network, piezo_resonator) |
| `post_analyses[].args` | object | no |  | arguments object passed to `cmd` |
| `post_analyses[].cmd` | text | no |  | explicit read-only post-solve command name (analyse, integrate, confinement, extremum, probe, line_cut, overlap, result_info, result_field, list_analyses, field, field_op) |
| `post_analyses[].label` | text | no |  | optional label attached to this post-analysis entry in the reply |
| `post_analyses[].params` | object | no |  | optional parameters object for shorthand `analysis` |
| `post_analyses[].physics` | text | no |  | three-letter physics code for shorthand `analysis`; defaults to the solved physics that supports the analysis |
| `save_fields` | true or false | no |  | when false, omit mesh and nodal/cell field arrays from worker HDF5 I/O and stored results after running any post_analyses, keeping scalars and eigenvalues; default true |
| `shared` | true or false | no |  | hand the solve to the queue of the workspace, which every program that uses the workspace shares (the window, the MCP server, `nyx queue`): a queue worker process runs the items one at a time in a worker, and keeps running when the program that asked is closed. Default: what the program was set up to do; implies a worker process |
| `wait` | true or false | no |  | wait for the result, default true |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"solve_run","args":{}}
```
