# `queue_submit`

**Submit work to the queue** · access `expensive`: starts a mesh build or a solve, or writes files

Queues a `solve` (default; `physics` is a code or a list), a `mesh` build or a mesh `convergence` check (`physics`, `factor`, `tolerance`) and returns its entry at once. The heavy lane runs one job at a time in order of priority (1 a person's request, 2 an agent's, the default, 3 sweeps and batches; a waiting job rises one level after a minute), then in turns between owners, then in order of submission. The job carries a frozen copy of the project: later edits do not change it (`project_changed` says they happened). `supersede` replaces a queued request of the same kind and targets; `hold` keeps it from starting until `queue_release`. A job that comes after another of the same project reuses the mesh the first built.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `factor` | number | no |  | convergence: finer/present element size, default 0.6 |
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `hold` | true or false | no |  | start held |
| `isolate` | true or false | no |  | solve in a worker process (only for kind solve): running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do |
| `kind` | text: `solve`, `mesh`, `convergence` | no |  | default solve |
| `physics` | text or list of text | no |  | physics code(s); default the active one |
| `priority` | integer | no |  | 1, 2 or 3; lower runs first |
| `supersede` | true or false | no |  | replace a queued request of the same kind and targets |
| `tolerance` | number | no |  | convergence: largest accepted relative difference |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"queue_submit","args":{}}
```
