# `add_pml_conditions`

**Add the absorbing-layer conditions** · access `edit`: changes the project

Adds a `pml` condition for every absorbing piece the geometry made with add_pml (replacing earlier ones on those pieces), so the physics absorbs where the geometry says. `strength` is the imaginary stretch at the outer edge (default 3).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `strength` | number | no |  | imaginary stretch at the outer edge, default 3 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_pml_conditions","args":{"physics":"dfs"}}
```
