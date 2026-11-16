# `solve_run`

**Solve** · access `expensive`: starts a mesh build or a solve, or writes files

Solves the active physics (or the named ones, in coupling order), meshing first when needed. With wait false it starts a job. When the project is bound to a workflow the run is recorded.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text or list of text | no |  | physics code(s); default the active one |
| `wait` | true or false | no |  | wait for the result, default true |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"solve_run","args":{}}
```
