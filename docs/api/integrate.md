# `integrate`

**Integrate an expression over a region** · access `read`: reads the project or the results

The integral of an expression of fields over a region (a region name, comma list or `name*`; default the whole mesh), with the degree-4 rule of the finite elements (exact for the product of two P2 fields). A real expression answers a number, a complex one [re, im]. `measure` is the area or volume of the region. `revolution: true` weights by 2 pi r for a body of revolution (a 2D r-z mesh). Answers `stale` and `age_s` of what it read. FORMULAS (every `expr` of the reductions takes one): a text such as `abs(Ex)^2 * eps_r` or `where(x > 0, T, 0)`, parsed once and evaluated at each point in Rust. Numbers (`2j` is imaginary), `pi`, `e`, `+ - * / ^`, comparisons, `and or not`; functions `abs abs2 re im conj sqrt exp sin cos arg log log10 tan sinh cosh tanh atan step min max atan2` and `where(cond, a, b)` (only the branch that applies is evaluated, so `where(x != 0, 1/x, 0)` is safe). Names: the coordinates `x y z`; the fields of the result (`Ex`, `T`; in backticks when the name is not an identifier, `` `|E|` ``; `` `ems:Ex` `` or `` `ems:1:Ex` `` in backticks for another physics or solution); material values by their descriptor name (`eps_r`, `rho`; `eps` means `eps_r`; `material("rho")` to be explicit); `region("Si")` is 1 inside the region and 0 outside. A name like `E` stands for the vector (Ex, Ey, Ez) when the solution has those components: `dot(a, b)`, `inner(a, b)` (conjugates a), `cross`, `norm(v)`, `norm2(v)`, `E.x`. Orders, comparisons, `log` and `tan` read real values: wrap a complex one in `abs`, `re` or `im`. A material value is the one the solver saw in each cell. A tensor property (`eps_r`, `k`, `alpha`) has its entries as `eps_r_xx`, `eps_r_yy`, `eps_r_zz`, `eps_r_xy`, `eps_r_xz`, `eps_r_yz` in wafer axes with the crystal orientation applied (an isotropic material is the scalar times the identity); the plain `eps_r` is refused on a mesh where some material is a tensor. BOUNDARY INTEGRALS: with `boundary` (an outer boundary such as "window.top", or an interface "oxide|silicon" between two regions) the sum is over the facets of the boundary instead of the cells: a line integral on a 2D mesh, a surface integral on a 3D one, with `revolution` weighting by 2 pi r. The formula may use `nx`, `ny`, `nz`, the unit normal in wafer axes, so the flux of a vector F is `dot(F, [nx, ny, nz])`. At an interface the fields and materials differ on the two sides, so give `side`, the region to read them in; the normal points out of it. On an outer boundary no `side` is needed. `wavelength` (m) reads the material values at that wavelength, not the study's; the fields stay those of the solve.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `boundary` | text | no |  | integrate over the facets of this boundary (an outer boundary name, or an interface "a\|b" between two regions) instead of over cells: a line integral on a 2D mesh, a surface integral on 3D. The formula may use nx, ny, nz, the unit normal. Not with region |
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `side` | text | no |  | with boundary: the region an interface is read in (the fields and materials differ across it); the normal points out of it. Not needed on an outer boundary |
| `solution` | integer | no |  | default 0 |
| `wavelength` | number | no |  | metres: where the material values of a formula (eps_r of a dispersive model) are read; default the wavelength of the study. The fields stay those of the solve |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"integrate","args":{"expr":"<expr>"}}
```
