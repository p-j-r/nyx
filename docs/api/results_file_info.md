# `results_file_info`

**Describe an HDF5 results file** · access `read`: reads the project or the results

The project title, the mesh size and for each physics the study, size, time, number of solutions, scalars, solver diagnostics and numerical error budget of a results file.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `path` | text | yes |  | file path |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"results_file_info","args":{"path":"<path>"}}
```
