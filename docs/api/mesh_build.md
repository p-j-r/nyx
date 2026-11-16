# `mesh_build`

**Mesh the geometry** · access `expensive`: starts a mesh build or a solve, or writes files

Meshes the geometry with the current settings. With wait false it starts a job and returns its id.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `wait` | true or false | no |  | wait for the mesh, default true |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"mesh_build","args":{}}
```
