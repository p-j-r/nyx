# `add_torus`

**Torus** · access `edit`: changes the project

A torus with axis z (3D). One face, `surface`. [shape; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `center` | list of 3 (number or text) | yes |  | centre. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `major` | number or text | yes |  | major radius. a number in metres (or the unit of the field) or an expression over the model parameters |
| `minor` | number or text | yes |  | minor radius. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_torus","args":{"center":[0.0,0.0,0.0],"major":0.0,"minor":0.0}}
```
