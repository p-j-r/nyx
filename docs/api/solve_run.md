# `solve_run`

**Solve** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the active physics (or the named ones, in coupling order), meshing first when needed. With wait false it starts a job. When the project is bound to a workflow the run is recorded.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `force` | true or false | no |  | start even when the memory estimate says this machine cannot hold the solve (an estimate exists once earlier solves of this kind ran in workers); without it such a solve is refused, and one that fits the machine but not the memory free now waits |
| `isolate` | true or false | no |  | run the solve in a worker process: running out of memory or crashing fails this job and not the program, and a cancel stops it at once; the cost is a results round trip. Default: what the program was set up to do |
| `physics` | text or list of text | no |  | physics code(s); default the active one |
| `wait` | true or false | no |  | wait for the result, default true |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"solve_run","args":{}}
```
