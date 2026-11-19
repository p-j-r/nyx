# `add_revolve`

**Revolve** · access `edit`: changes the project

A full turn of a Rect or Polygon profile about the axis through `center` (3D). The profile's (x, y) are the radius (>= 0) and the height. [build; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text: `X`, `Y`, `Z` | no | `Z` | the axis of revolution |
| `center` | list of 3 (number or text) | yes |  | a point on the axis. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `profile` | text or integer | yes |  | the 2D profile feature. the name of an earlier feature (or its id) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_revolve","args":{"center":[0.0,0.0,0.0],"profile":"<profile>"}}
```
