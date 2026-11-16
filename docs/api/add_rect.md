# `add_rect`

**Rectangle** · access `edit`: changes the project

A rectangle in the x-z cross-section, given by its two corners. [shape; dimensions [2]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `max` | list of 2 (number or text) | yes |  | upper-right corner (x, z) in metres. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `min` | list of 2 (number or text) | yes |  | lower-left corner (x, z) in metres. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_rect","args":{"max":[0.0,0.0],"min":[0.0,0.0]}}
```
