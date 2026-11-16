# `add_box`

**Box** · access `edit`: changes the project

An axis-aligned box (3D). Face roles: xmin xmax ymin ymax zmin zmax. [shape; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `max` | list of 3 (number or text) | yes |  | upper corner (x, y, z). each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `min` | list of 3 (number or text) | yes |  | lower corner (x, y, z). each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_box","args":{"max":[0.0,0.0,0.0],"min":[0.0,0.0,0.0]}}
```
