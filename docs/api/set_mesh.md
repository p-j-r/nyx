# `set_mesh`

**Set the mesh settings** · access `edit`: changes the project

Give only what changes. `adapt_scale`, `adapt_grade`, `adapt_max_cells`, `gap_cells`, `singular_edges`, `quarter_point_xi` or `sweep` set to null return to their defaults. `sweep` extrudes a 2D cross-section or a 3D geometry model along y through layered intervals with optional `fill` or per-layer `regions` remap.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `adapt_grade` | number or null | no |  | adaptive / AMR metric gradation slope (> 0); null uses 0.4 |
| `adapt_max_cells` | integer or null | no |  | adaptive / AMR maximum cell count (>= 1); null uses 2000000 |
| `adapt_scale` | number or null | no |  | adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher |
| `corner_grading` | true or false | no |  | 3D sheet stack: radical Morrison-Lewis corner grading toward sharp 2D sheet corners, default false |
| `curved` | true or false | no |  | 3D: curved (quadratic) isoparametric geometry for cells on curved surfaces and fitted sheet arc rims (ems, dfs, and sheet3d), default false |
| `gap_cells` | number or null | no |  | 3D sheet stack: target minimum number of 2D triangles across narrow coplanar gaps between conductor sheets (positive number, e.g. 2.5); null disables automatic gap clamping |
| `h` | number | no |  | global element size in metres |
| `max_dev` | number or null | no |  | largest chord deviation from a curved edge or face, metres; null: h / 50 for the adaptive mesher and no limit for the simple 3D mesher |
| `physics_sizing` | true or false | no |  | size the elements from the waves of the physics (points per wavelength) |
| `points_per_wavelength` | number or null | no |  | points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2 |
| `quarter_point_xi` | number or null | no |  | 3D sheet stack: Barsoum quarter-point mid-edge parameter xi_QP in [0.25, 0.50]; null uses 0.25 |
| `singular_edges` | text or null: `none`, `quarter_point`, `null` | no |  | 3D sheet stack: edge singularity enrichment on zero-thickness conductor rims ("none" or "quarter_point"); null uses study defaults |
| `structured_pml` | true or false | no |  | rebuild the absorbing layers as structured layers |
| `sweep` | list of object or null | no |  | layers to sweep a 2D cross-section or 3D geometry model along y |
| `sweep[].domain_offset` | integer | no |  | added to the domain numbers of the layer's cells |
| `sweep[].fill` | text | no |  | region name: every cell of the layer belongs to it |
| `sweep[].n` | integer | yes |  | cells across the thickness |
| `sweep[].regions` | object | no |  | optional per-layer region remap { source_region: target_region } in this layer |
| `sweep[].thickness` | number | yes |  | thickness of the layer along y, metres |
| `z_coarsen` | true or false | no |  | 3D sheet stack: progressively coarsen the 2D sheet triangulation away from the conductor plane along z via conforming pure-tetrahedral transition slabs, default false |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_mesh","args":{}}
```
