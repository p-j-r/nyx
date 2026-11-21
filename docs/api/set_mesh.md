# `set_mesh`

**Set the mesh settings** · access `edit`: changes the project

Give only what changes. `adapt_scale`, `gap_cells` or `sweep` set to null return to the simple mesher / explicit-only sheet sizing / an unswept mesh. `sweep` extrudes a 2D cross-section or a 3D geometry model along y through layered intervals with optional `fill` or per-layer `regions` remap.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `adapt_scale` | number or null | no |  | adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher |
| `curved` | true or false | no |  | 3D: curved (quadratic) geometry for the cells on a curved surface; used by the ems cavity and bloch studies, default false |
| `gap_cells` | number or null | no |  | 3D sheet stack: target minimum number of 2D triangles across narrow coplanar gaps between conductor sheets (positive number, e.g. 2.5); null disables automatic gap clamping |
| `h` | number | no |  | global element size in metres |
| `max_dev` | number or null | no |  | largest chord deviation from a curved edge or face, metres; null: h / 50 for the adaptive mesher and no limit for the simple 3D mesher |
| `physics_sizing` | true or false | no |  | size the elements from the waves of the physics (points per wavelength) |
| `points_per_wavelength` | number or null | no |  | points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2 |
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
