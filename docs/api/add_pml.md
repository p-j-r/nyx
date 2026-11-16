# `add_pml`

**Absorbing layers** · access `edit`: changes the project

Absorbing layers of `thickness` around a window on some sides; splits every matching region into the physical part and one piece per side combination. [boundary; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `grading` | number or text | yes |  | order of the polynomial profile. a number in metres (or the unit of the field) or an expression over the model parameters |
| `layers` | number or text | yes |  | number of element layers across it. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `regions` | list of text | no | `[]` | regions to split (empty: all) |
| `sides` | list of text | no | `[]` | sides that get a layer (empty: all) |
| `thickness` | number or text | yes |  | layer thickness. a number in metres (or the unit of the field) or an expression over the model parameters |
| `window` | text or integer | yes |  | the feature that bounds the physical window. the name of an earlier feature (or its id) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_pml","args":{"grading":0.0,"layers":0.0,"thickness":0.0,"window":"<window>"}}
```
