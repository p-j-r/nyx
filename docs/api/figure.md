# `figure`

**Draw the modes of a result** · access `read`: reads the project or the results

Writes a figure of the modes of the latest result of a physics as SVG or PNG (the extension of `path` decides). Each mode is put in a fixed phase first (its in-plane motion largest in Re u, sign fixed, peak 1) and drawn on the section of the mesh (the whole mesh of a 2D result, the first y-plane of a 3D one; a displacement deforms it by Re u), with the outline of the regions and arrows of the in-plane motion. `gallery` shows |u| for each mode; `components` a row per mode with |u|, the in-plane components and the propagation component in quadrature. The reply says, per mode, the angle it was turned by and `quadrature`: the share of the propagation component's energy still in Re after the phase step (near 0 when the mode is in the expected phase). A vector field is found by the names of the cell fields (ux uy uz, Ex Ey Ez).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `arrows` | true or false | no |  | arrows of the in-plane motion, default true |
| `at` | number | no |  | the y of the plane a 3D mesh is cut on, in metres; default its first plane |
| `cmap` | text | no |  | colour map of \|u\|: inferno (default), viridis, cividis, coolwarm |
| `columns` | integer | no |  | panels per row of a gallery, default 4 |
| `core` | list of integer | no |  | domain numbers that carry the mode (the guide, not the substrate): the phase and the scale are taken over them; default all |
| `family` | text | no |  | the vector field to draw, by its prefix (u, E); default the first |
| `kind` | text | no |  | gallery (default) or components |
| `max_modes` | integer | no |  | how many modes when `modes` is not given, default 12 |
| `modes` | list of integer | no |  | the modes to draw, numbered from 1; default the first max_modes |
| `panel_width` | number | no |  | width of a panel in drawing units, default 250 |
| `path` | text | yes |  | file to write, .svg or .png |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `scale` | number | no |  | pixels per drawing unit of a PNG, default 2 |
| `warp` | number | no |  | largest motion of the deformed shape as a fraction of the window width; default 0.05 for a displacement (vector field `u`) and 0 for any other field; 0 draws the undeformed section |
| `window` | list of number | no |  | [x0, x1, z0, z1] in metres; default the extent of the modes |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"figure","args":{"path":"<path>","physics":"dfs"}}
```
