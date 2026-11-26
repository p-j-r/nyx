# `vvuq_export`

**Export canonical VVUQ provenance and error-budget report** · access `expensive`: starts a mesh build or a solve, or writes files

Writes an unsigned, canonical JSON VVUQ provenance and error-budget report (`vvuq_report.v1`) linking BLAKE3 digests of canonical project JSON, solver mesh, active/resolved materials, and linked/saved `results.h5` with build metadata, `SimplexMesh::check()` summary, per-physics `SolverDiagnostics`, `NumericalErrorBudget`, `Validated` state/tier, and active material warnings. Pass `export_h5: true` to also write a companion `results.h5` file beside the report when solved and none is linked yet.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `export_h5` | true or false | no |  | also export a companion results.h5 file beside the report if solved and no results.h5 is linked yet, default false |
| `path` | text | yes |  | output JSON file path for the vvuq_report.v1 report |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"vvuq_export","args":{"path":"<path>"}}
```
