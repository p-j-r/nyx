# `add_polygons`

**Polygon batch** · access `edit`: changes the project

A batch of outer polygon rings and optional hole rings in 2D or 3D, with Ramer-Douglas-Peucker closed-ring simplification (shared corners pinned) and minimum hole area filtering. [shape; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `holes` | list of list of list of 2 (number or text) | no | `[]` | hole polygon rings to subtract |
| `min_hole_area` | number or text | no | `0.0` | drop hole rings with unsigned shoelace area below this threshold in m^2. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `polygons` | list of list of list of 2 (number or text) | yes |  | outer polygon rings, each a list of 2D corners in metres |
| `simplify_tol` | number or text | no | `5e-8` | Ramer-Douglas-Peucker simplification tolerance in metres. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_polygons","args":{"polygons":[[[0.0,0.0],[0.0,0.0],[0.0,0.0]]]}}
```
