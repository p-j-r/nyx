# `set_mesh`

**Set the mesh settings** · access `edit`: changes the project

Give only what changes. `adapt_scale` or `sweep` set to null return to the simple mesher / a 2D mesh.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `adapt_scale` | number or null | no |  | adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher |
| `h` | number | no |  | global element size in metres |
| `max_dev` | number or null | no |  | largest chord deviation from a curved edge, metres |
| `physics_sizing` | true or false | no |  | size the elements from the waves of the physics (points per wavelength) |
| `points_per_wavelength` | number or null | no |  | points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2 |
| `structured_pml` | true or false | no |  | rebuild the absorbing layers as structured layers |
| `sweep` | list of object or null | no |  | layers to sweep the 2D mesh along y |
| `sweep[].domain_offset` | integer | no |  | added to the domain numbers of the layer's cells |
| `sweep[].n` | integer | yes |  | cells across the thickness |
| `sweep[].thickness` | number | yes |  | thickness of the layer along y, metres |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_mesh","args":{}}
```
