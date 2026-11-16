# `open_results`

**Open an HDF5 results file** · access `edit`: changes the project

Makes the project of a results file the project, builds its mesh again (meshing is deterministic) and checks it is the stored mesh, then takes the stored results as the latest results, so every analysis, plot and export works on them as after a solve. Replaces the present project. Needs a build with HDF5 support. With `wait` false it starts a job (reading and meshing happen on it; the project is replaced when job_wait or poll_jobs applies it).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `path` | text | yes |  | file path, .h5 |
| `wait` | true or false | no |  | wait for it, default true |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"open_results","args":{"path":"<path>"}}
```
