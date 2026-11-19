# `add_cut`

**Cut by a plane** · access `edit`: changes the project

Keeps one side of the plane normal . p = offset. The cut face is named <feature>.cut and tagged `symmetry`. [boundary; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `keep` | text: `Above`, `Below` | yes |  | which side to keep |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `normal` | list of 3 (number or text) | yes |  | plane normal. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `offset` | number or text | yes |  | plane offset along the normal. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_cut","args":{"keep":"Above","normal":[0.0,0.0,0.0],"offset":0.0}}
```
