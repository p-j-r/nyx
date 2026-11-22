# `add_chip_stack`

**Add a 3D chip dielectric stack** · access `edit`: changes the project

Constructs a 3D chip dielectric stack (`VacBot`, `Si`, `VacTop` boxes and `Overlay`) over a 2D `[xmin, ymin, xmax, ymax]` bounding box (explicitly given or inferred from existing 2D/sheet features in the project), expanded laterally by `lateral_pad >= 0`.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `bbox` | list of 4 (number) | no |  | optional unpadded 2D bounding box [xmin, ymin, xmax, ymax] in metres; omit to infer from existing 2D/sheet features |
| `h_vac_bot` | number | no |  | bottom vacuum thickness below the substrate in metres (0 omits bottom vacuum; aliases: t_vac_bot, vac_bot_thickness), default h_vac_top |
| `h_vac_top` | number | no |  | top vacuum thickness above z in metres (aliases: t_vac_top, vac_top_thickness), default 3 * t_sub |
| `lateral_pad` | number | no |  | non-negative lateral buffer margin in metres added on all 4 horizontal sides, default 0 |
| `overlay` | true or false | no |  | create or update the project's Overlay feature with the stack regions, default true |
| `substrate_box` | text | no |  | feature name of the substrate Box, default "si_sub" |
| `substrate_region` | text | no |  | region name of the substrate in Overlay, default "Si" |
| `substrate_thickness` | number | no |  | alias for t_sub: substrate thickness in metres |
| `t_sub` | number | no |  | substrate thickness in metres (alias: substrate_thickness), default 350e-6 |
| `t_vac_bot` | number | no |  | alias for h_vac_bot: bottom vacuum thickness in metres |
| `t_vac_top` | number | no |  | alias for h_vac_top: top vacuum thickness in metres |
| `vac_bot_box` | text | no |  | feature name of the bottom vacuum Box, default "vac_bot" |
| `vac_bot_region` | text | no |  | region name of the bottom vacuum in Overlay, default "VacBot" |
| `vac_bot_thickness` | number | no |  | alias for h_vac_bot: bottom vacuum thickness in metres |
| `vac_top_box` | text | no |  | feature name of the top vacuum Box, default "vac_top" |
| `vac_top_region` | text | no |  | region name of the top vacuum in Overlay, default "VacTop" |
| `vac_top_thickness` | number | no |  | alias for h_vac_top: top vacuum thickness in metres |
| `z` | number | no |  | wafer-normal z coordinate of the top substrate interface in metres, default 0 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"add_chip_stack","args":{}}
```
