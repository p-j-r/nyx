# `add_polygon`

**Polygon** · access `edit`: changes the project

A closed polygon in the x-z cross-section; counter-clockwise or clockwise. Edge i goes from point i to point i+1. [shape; dimensions [2]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `names` | list of text | no | `[]` | role names of the edges (default e0, e1, ...), used in boundary names |
| `pts` | list of list of 2 (number or text) | yes |  | the corners (x, z) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_polygon","args":{"pts":[[0.0,0.0],[0.0,0.0],[0.0,0.0]]}}
```
