# `export_results`

**Write the results to an HDF5 file** · access `expensive`: starts a mesh build or a solve, or writes files

Writes the project, the mesh and the latest result of every solved physics (all fields, compressed, exact) into an HDF5 file and an XDMF index of the stored datasets next to it, for nyx and its tools only (not for ParaView or VisIt). Needs a build with HDF5 support; `parallel` also needs `phdf5`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `parallel` | true or false | no |  | write with the MPI-IO driver; every process of an MPI run calls this with the same arguments and writes its own rows (needs a build with the `phdf5` feature) |
| `path` | text | yes |  | file path, .h5 |
| `wait` | true or false | no |  | wait for the file, default true; false starts a job that writes what the engine holds now, under its name only when complete |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"export_results","args":{"path":"<path>"}}
```
