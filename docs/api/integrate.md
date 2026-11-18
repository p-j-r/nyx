# `integrate`

**Integrate an expression over a region** · access `read`: reads the project or the results

The integral of an expression of fields over a region (a region name, comma list or `name*`; default the whole mesh), with the degree-4 rule of the finite elements (exact for the product of two P2 fields). A real expression answers a number, a complex one [re, im]. `measure` is the area or volume of the region. `revolution: true` weights by 2 pi r for a body of revolution (a 2D r-z mesh). Answers `stale` and `age_s` of what it read. FORMULAS (every `expr` of the reductions takes one): a text such as `abs(Ex)^2 * eps_r` or `where(x > 0, T, 0)`, parsed once and evaluated at each point in Rust. Numbers (`2j` is imaginary), `pi`, `e`, `+ - * / ^`, comparisons, `and or not`; functions `abs abs2 re im conj sqrt exp sin cos arg log log10 tan sinh cosh tanh atan step min max atan2` and `where(cond, a, b)` (only the branch that applies is evaluated, so `where(x != 0, 1/x, 0)` is safe). Names: the coordinates `x y z`; the fields of the result (`Ex`, `T`; in backticks when the name is not an identifier, `` `|E|` ``; `` `ems:Ex` `` or `` `ems:1:Ex` `` in backticks for another physics or solution); material values by their descriptor name (`eps_r`, `rho`; `eps` means `eps_r`; `material("rho")` to be explicit); `region("Si")` is 1 inside the region and 0 outside. A name like `E` stands for the vector (Ex, Ey, Ez) when the solution has those components: `dot(a, b)`, `inner(a, b)` (conjugates a), `cross`, `norm(v)`, `norm2(v)`, `E.x`. Orders, comparisons, `log` and `tan` read real values: wrap a complex one in `abs`, `re` or `im`. A material value is a tensor-free scalar of the problem's materials as the solver saw them; an anisotropic `eps_r` is refused.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `expr` | text or number or object or list of any | yes |  | a formula text such as "abs(Ex)^2 * eps_r" (see integrate), a handle, a field name, a number or {op, args} |
| `physics` | text: `dfs`, `ems`, `sms` | no |  | the three-letter code of the physics |
| `region` | text | no |  | region name(s); default all |
| `revolution` | true or false | no |  | weight by 2 pi r, default false |
| `solution` | integer | no |  | default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"integrate","args":{"expr":"<expr>"}}
```
