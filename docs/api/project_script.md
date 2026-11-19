# `project_script`

**Project as a script** · access `read`: reads the project or the results

The project as a script (Python, `.nys`) that rebuilds it: project_new, set_param (each parameter reads `params`), one add_* per feature, set_mesh, add_material, set_problem, set_active. It does not solve. It is replayed on a scratch engine and offered only when it gives the same project; `form` says how it is written: `commands`, `geometry_json` (the geometry as one set_geometry, for ids and pins that must be kept) or `project_json` (one project_load of the whole project: exact, not readable). Starts with project_new, so over MCP it needs allow_destructive.

## Arguments

None.

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"project_script","args":{}}
```
