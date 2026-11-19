# `set_mesh`

**Set the mesh settings** · access `edit`: changes the project

Give only what changes. `adapt_scale` or `sweep` set to null return to the simple mesher / a 2D mesh.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `adapt_scale` | number or null | no |  | adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher |
| `curved` | true or false | no |  | 3D: curved (quadratic) geometry for the cells on a curved surface; used by the ems cavity study, default false |
| `h` | number | no |  | global element size in metres |
| `max_dev` | number or null | no |  | largest chord deviation from a curved edge or face, metres; null: h / 50 for the adaptive mesher and no limit for the simple 3D mesher |
| `physics_sizing` | true or false | no |  | size the elements from the waves of the physics (points per wavelength) |
| `points_per_wavelength` | number or null | no |  | points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2 |
| `structured_pml` | true or false | no |  | rebuild the absorbing layers as structured layers |
| `sweep` | list of object or null | no |  | layers to sweep the 2D mesh along y |
| `sweep[].domain_offset` | integer | no |  | added to the domain numbers of the layer's cells |
| `sweep[].fill` | text | no |  | region name: every cell of the layer belongs to it |
| `sweep[].n` | integer | yes |  | cells across the thickness |
| `sweep[].thickness` | number | yes |  | thickness of the layer along y, metres |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_mesh","args":{}}
```
