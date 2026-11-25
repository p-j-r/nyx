# `add_absorbing_layer`

**Add an outward absorbing layer** · access `edit`: changes the project

Adds an outward zero-padding absorbing layer (`Op::Pml` with `mode = "Outward"`, `window = null`, `coord = "Auto"`, `layer_grading = 2.0`, `thickness = "auto"`, `layers = "auto"`) around the exterior of the physical domain across Cartesian, Cylindrical, or Spherical topologies.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text or null: `X`, `Y`, `Z`, `null` | no |  | optional cylinder axis (X, Y, or Z) |
| `center` | list of number or text or null | no |  | optional centre [x, y, z] (or [x, z] in 2D) for Cylindrical or Spherical shells |
| `coord` | text: `Auto`, `Cartesian`, `Cylindrical`, `Spherical`, `auto`, `cartesian`, `cylindrical`, `spherical` | no |  | coordinate topology, default "Auto" |
| `grading` | number or text | no |  | polynomial stretch profile order m, default 2.0 |
| `layer_grading` | number or text or null | no |  | geometric mesh layer grading ratio (>= 1; default 2.0) |
| `layers` | number or text | no |  | number of structured element layers across the shell, or "auto" (default "auto") |
| `mode` | text: `Inward`, `Outward`, `inward`, `outward` | no |  | growth mode, default "Outward" |
| `name` | text | no |  | feature name; omit for pml1, pml2, ... |
| `regions` | list of text | no |  | regions to extrude (empty: all) |
| `sides` | list of text | no |  | active sides (empty: all) |
| `thickness` | number or text | no |  | layer thickness in metres, or "auto" (default "auto") |
| `window` | text or integer or null | no |  | optional window feature when mode is Inward |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_absorbing_layer","args":{}}
```
