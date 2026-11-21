# `confinement`

**Fraction of a field in a region** · access `read`: reads the project or the results

The integral of sum_k |c_k|^2 over one or more regions divided by the same over the whole mesh, for the components `fields` (e.g. [Ex, Ey, Ez]). Without `weight` it has no material weight; give `weight` a real formula such as "eps_r" for the energy fraction (the integrand is weight * sum_k |c_k|^2, in each region and over the whole mesh). Pass `solutions` (`"all"` or a list of 0-based indices), `by_region: true`, or `regions` (`[str]`) to accumulate all requested solutions and regions in a single bit-exact quadrature pass.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `by_region` | true or false | no |  | compute confinement fractions across every region of the mesh (or within `region`) in a single pass, default false |
| `expr` | text or number or object or list of any | no |  | optional alias for `fields`: a field, formula, expression, or list of components |
| `fields` | text or number or object or list of any | yes |  | a field, a formula, an expression, or a list of components |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text or list of text | no |  | the region (or list of regions) the fraction is of; omit when `regions` or `by_region: true` is given |
| `regions` | list of text | no |  | batch list of region names to compute confinement fractions over in a single pass |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | 0-based solution index, default 0 (mutually exclusive with solutions) |
| `solutions` | text or list of integer | no |  | "all" or a list of 0-based solution indices to evaluate in a single pass |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |
| `weight` | text | no |  | a real formula that multiplies the integrand at each point, e.g. eps_r; default none |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"confinement","args":{"fields":"<fields>"}}
```
