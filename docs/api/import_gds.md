# `import_gds`

**Import GDSII conductive sheets** · access `edit`: changes the project

Reads a binary GDSII (.gds) stream in a 3D project, flattens SREF/AREF hierarchy and PATH/BOUNDARY/BOX geometry into SI metres, unions each requested [layer, datatype] pair, optionally splits disconnected conductive islands and renames islands from TEXT labels on `pin_layer`, and adds `Polygons` + `Sheet` features at height `z`. Returns `{ sheets, layers, bbox, pins }`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `cell` | text | no |  | top structure cell name to flatten; omit to auto-select the unreferenced root cell |
| `layer_bbox` | true or false | no |  | when true (default), return bbox over imported layers only; when false, return full top-cell bbox |
| `layers` | object | yes |  | map from base sheet name to [layer, datatype] GDSII integer pair, e.g. {"gnd": [1, 0], "sig": [1, 10]} |
| `path` | text | yes |  | file path to the binary GDSII (.gds) stream |
| `pin_layer` | list of 2 (integer) | no |  | [layer, texttype] GDSII pair whose TEXT labels rename matching conductive islands |
| `simplify_tol` | number | no |  | collinear/chord vertex simplification tolerance in metres (shared cross-layer contact vertices stay pinned), default 0 |
| `split_islands` | true or false | no |  | split disconnected polygons on each layer into separate sheets (<name>_0, <name>_1, ... or matched TEXT pin names), default true |
| `z` | number | no |  | wafer-normal z height of the imported conductive sheets in metres, default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"import_gds","args":{"layers":{},"path":"<path>"}}
```
