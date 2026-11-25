# `add_pml_conditions`

**Add the absorbing-layer conditions** · access `edit`: changes the project

Adds a `pml` (or `infinite_element` for static studies) condition for every absorbing piece the geometry made with `add_pml` or `add_absorbing_layer` (replacing earlier ones on those pieces), so the physics absorbs where the geometry says. When `strength <= 0` or `auto_scale: true`, parameters are jointly auto-calibrated across materials sharing each PML slab/shell from the study's wavelengths.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `auto_scale` | true or false | no |  | jointly auto-calibrate strength, profile, and kappa_max/curvature from study wavelengths across materials sharing each PML slab/shell |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `profile` | text or null: `auto`, `kappa_polynomial`, `rational`, `exponential`, `null` | no |  | stretching profile ("auto" selects kappa_polynomial for eigenproblems/waveguides and rational for direct driven/static problems) |
| `strength` | number | no |  | imaginary stretch at the outer edge (default 3; <= 0 auto-calibrates from study wavelengths) |
| `target_reflection` | number | no |  | target normal-incidence amplitude reflection R_0 when auto_scale is enabled, default 1e-6 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_pml_conditions","args":{"physics":"dfs"}}
```
