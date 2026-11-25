# `pml_diagnostics`

**Absorbing-layer mesh and stability diagnostics** · access `read`: reads the project or the results

Returns absorbing-layer diagnostics across all PML domains: uncut `element_counts`, `min_pml_scaled_jacobian`, `pml_warnings`, normalized penetration depth `s_dist_bounds` `[0.0, 1.0]`, and pre-solve Christoffel backward-wave check (`has_backward_wave`, `min_projection`, `recommended_mpml_ratio`).

## Arguments

None.

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"pml_diagnostics","args":{}}
```
