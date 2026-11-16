# `add_cylinder`

**Cylinder** · access `edit`: changes the project

A cylinder (3D) from `base` along +z. Face roles: side bottom top. [shape; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text: `X`, `Y`, `Z` | no | `Z` | this build: z only |
| `base` | list of 3 (number or text) | yes |  | centre of the bottom face. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `h` | number or text | yes |  | height. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `r` | number or text | yes |  | radius. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_cylinder","args":{"base":[0.0,0.0,0.0],"h":0.0,"r":0.0}}
```
