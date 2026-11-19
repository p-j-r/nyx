# `add_cone`

**Cone frustum** · access `edit`: changes the project

A cone frustum (3D) along the axis from `base`; either radius may be 0. Faces side, bottom, top. [shape; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text: `X`, `Y`, `Z` | no | `Z` | the axis of the cone |
| `base` | list of 3 (number or text) | yes |  | centre of the bottom face. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `h` | number or text | yes |  | height. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `r0` | number or text | yes |  | radius at the bottom. a number in metres (or the unit of the field) or an expression over the model parameters |
| `r1` | number or text | yes |  | radius at the top. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_cone","args":{"base":[0.0,0.0,0.0],"h":0.0,"r0":0.0,"r1":0.0}}
```
