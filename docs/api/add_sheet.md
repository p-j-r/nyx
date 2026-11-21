# `add_sheet`

**2D sheet or conformal slab in 3D** · access `edit`: changes the project

Embeds a 2D profile or boolean feature at plane z (3D) as a zero-thickness conductive sheet (default) or a conformal finite-thickness conductor slab (`thickness > 0`) with optional substrate over-etch trench (`etch_depth > 0`). [build; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `etch_depth` | number or text | no | `0.0` | substrate over-etch trench depth below z in metres (0 for unetched interface). a number in metres (or the unit of the field) or an expression over the model parameters |
| `input` | text or integer | yes |  | the 2D profile or boolean feature. the name of an earlier feature (or its id) |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `thickness` | number or text | no | `0.0` | conformal metal film thickness above z in metres (0 for zero-thickness sheet). a number in metres (or the unit of the field) or an expression over the model parameters |
| `z` | number or text | no | `0.0` | z coordinate of the sheet plane in metres. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_sheet","args":{"input":"<input>"}}
```
