# `add_extrude`

**Extrude** · access `edit`: changes the project

Sweeps a 2D profile feature along an axis (default y, the propagation direction) into a 3D solid. Roles side:<profile role> and the caps ymin, ymax. [build; dimensions [3]]

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `axis` | text: `X`, `Y`, `Z` | no | `Y` | the axis of the sweep |
| `from` | number or text | yes |  | start coordinate along the axis. a number in metres (or the unit of the field) or an expression over the model parameters |
| `name` | text | no |  | a name for the feature; default the kind and a number (`rect1`) |
| `profile` | text or integer | yes |  | the 2D profile feature. the name of an earlier feature (or its id) |
| `taper_deg` | number or text | no | `null` | wall angle in degrees (positive narrows towards `to`); makes the solid polyhedral. a number in metres (or the unit of the field) or an expression over the model parameters |
| `to` | number or text | yes |  | end coordinate along the axis. a number in metres (or the unit of the field) or an expression over the model parameters |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_extrude","args":{"from":0.0,"profile":"<profile>","to":0.0}}
```
