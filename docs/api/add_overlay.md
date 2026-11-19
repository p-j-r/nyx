# `add_overlay`

**Regions (overlay)** · access `edit`: changes the project

The root of the model: names the regions, in painter order (a later region wins where they overlap). A list has one overlay. Pins fix the integer label of selected boundaries. [naming; dimensions [2, 3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `naming` | text: `Feature`, `Region`, `Numbered` | no | `Feature` | how boundary faces against the void are named |
| `pins` | list of object | no | `[]` | label pins: { sel (selector text), id } |
| `pins[].id` | integer | yes |  |  |
| `pins[].sel` | text | yes |  |  |
| `regions` | list of object | yes |  | regions in painter order: { name, input (feature name), kind (Material \| Hole \| Pml \| Void), pinned_id } |
| `regions[].input` | text or integer | yes |  |  |
| `regions[].kind` | text: `Material`, `Hole`, `Pml`, `Void` | no | `Material` |  |
| `regions[].name` | text | yes |  |  |
| `regions[].pinned_id` | integer or null | no |  |  |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_overlay","args":{"regions":[]}}
```
