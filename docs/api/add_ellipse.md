# `add_ellipse`

**Ellipse** · access `edit`: changes the project

An ellipse in the x-z cross-section, rotated by an angle. [shape; dimensions [2]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `angle_deg` | number or text | no | `0.0` | rotation in degrees. a number in metres (or the unit of the field) or an expression over the model parameters |
| `center` | list of 2 (number or text) | yes |  | centre (x, z). each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `radii` | list of 2 (number or text) | yes |  | semi-axes along the rotated x and z. each item: a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_ellipse","args":{"center":[0.0,0.0],"radii":[0.0,0.0]}}
```
