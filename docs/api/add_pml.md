# `add_pml`

**Absorbing layers** · access `edit`: changes the project

Absorbing layers around a window (`mode = "Inward"`) or extruded outward from the physical exterior (`mode = "Outward"`, `window` omitted); splits or extrudes matching regions into `<region>.pml_<piece>` sub-domains across Cartesian, Cylindrical, or Spherical topologies. [boundary; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text: `X`, `Y`, `Z` | no | `null` | optional axis (X, Y, or Z) for Cylindrical shells |
| `center` | list of 3 (number or text) | no | `null` | optional centre [x, y, z] for Cylindrical or Spherical shells. each item: a number in metres (or the unit of the field) or an expression over the model parameters |
| `coord` | text: `Auto`, `Cartesian`, `Cylindrical`, `Spherical`, `auto`, `cartesian`, `cylindrical`, `spherical` | no | `Auto` | coordinate topology: Auto, Cartesian, Cylindrical, or Spherical |
| `grading` | number or text | no | `2.0` | order of the polynomial profile. a number in metres (or the unit of the field) or an expression over the model parameters |
| `layer_grading` | number or text | no | `null` | geometric mesh layer grading ratio (>= 1; 1 is uniform). a number in metres (or the unit of the field) or an expression over the model parameters |
| `layers` | number or text | no | `auto` | number of element layers across it (or "auto" / -1). a number in metres (or the unit of the field) or an expression over the model parameters |
| `mode` | text: `Inward`, `Outward`, `inward`, `outward` | no | `Inward` | growth mode: Inward (around window) or Outward (exterior boundary extrusion) |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `regions` | list of text | no | `[]` | regions to split or extrude (empty: all) |
| `sides` | list of text | no | `[]` | sides that get a layer (empty: all) |
| `thickness` | number or text | no | `auto` | layer thickness in metres (or "auto" / -1). a number in metres (or the unit of the field) or an expression over the model parameters |
| `window` | text or integer | no | `null` | the feature that bounds the physical window (required when mode is Inward; omitted for Outward). the name of an earlier feature (or its id) |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_pml","args":{}}
```
