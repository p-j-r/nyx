# `add_circle`

**Circle** · access `edit`: changes the project

A disc in the x-z cross-section. [shape; dimensions [2]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `center` | list of 2 (number or text) | yes |  | centre (x, z). each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `r` | number or text | yes |  | radius. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_circle","args":{"center":[0.0,0.0],"r":0.0}}
```
