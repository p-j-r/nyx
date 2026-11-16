# `add_size`

**Mesh size** · access `edit`: changes the project

Sets the element size near the selected boundaries or regions. [naming; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `h` | number or text | no | `null` | element size. a number in metres (or the unit of the field) or an expression over the model parameters |
| `max_deviation` | number or text | no | `null` | largest distance of a chord from its curve. a number in metres (or the unit of the field) or an expression over the model parameters |
| `min_layers` | number or text | no | `null` | fewest element layers across a thin region. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `sel` | text | yes |  | where the size applies. selector text, e.g. `outer & plane(x = 0)` |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_size","args":{"sel":"<sel>"}}
```
