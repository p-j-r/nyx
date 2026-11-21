# `figure`

**Draw the modes or 3D fields of a result** · access `read`: reads the project or the results

Writes a 2D or 3D figure of the latest result of a physics as SVG or PNG (the extension of `path` decides). 2D kinds (`gallery`, `components`) align each mode's phase and draw on the mesh section with region outlines and in-plane arrows. 3D kinds (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`) evaluate `expr` (e.g. `re(Ex)`, `abs2(E)`) directly on un-averaged 8-subtet P2 `CellFields` so normal-field discontinuities across material interfaces stay sharp, render with Z-buffered/translucent orthographic projection (`camera`, `show_edges`), and support zero-copy Bloch unit-cell tiling (`repeat_periodic` with `exp(-i m periodic_phase)` applied before `expr`).

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `arrows` | true or false | no |  | arrows of the in-plane motion, default true |
| `at` | number | no |  | the y of the plane a 3D mesh is cut on (or slice3d axis offset), in metres; default its first plane or midpoint |
| `camera` | object or list of 2 (number) | no |  | 3D orthographic camera angles in degrees: { azimuth_deg, elevation_deg } or [azimuth_deg, elevation_deg] |
| `clip` | object or list of 6 (number or null) | no |  | 3D cutaway clip bounds in metres: { xmin, xmax, ymin, ymax, zmin, zmax } (or x/y/z [min, max] pairs) or [xmin, xmax, ymin, ymax, zmin, zmax] |
| `cmap` | text: `inferno`, `viridis`, `cividis`, `coolwarm` | no |  | colour map: inferno (default for unsigned magnitudes), viridis, cividis, coolwarm (default symmetric map for signed 3D formulas) |
| `columns` | integer | no |  | panels per row of a gallery, default 4 |
| `core` | list of integer | no |  | domain numbers that carry the mode (the guide, not the substrate): the phase and the scale are taken over them; default all |
| `expr` | text | no |  | 3D scalar formula string evaluated per-cell without vertex averaging, e.g. "re(Ex)", "abs2(E)", "im(Hz)"; default \|u\| or \|E\| |
| `family` | text | no |  | the vector field to draw, by its prefix (u, E); default the first |
| `isosurfaces` | list of number | no |  | threshold levels for isosurface3d (default [0.25, 0.6] of peak, or [-0.6, -0.25, 0.25, 0.6] on signed formulas) |
| `kind` | text: `gallery`, `components`, `slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d` | no |  | figure layout: gallery (default), components, slice3d (3D planar slice), cutaway3d (3D clipped cutaway), isosurface3d (3D translucent marching-tets isosurfaces) or mesh3d (3D region/element surface view) |
| `max_modes` | integer | no |  | how many modes when `modes` is not given, default 12 |
| `modes` | list of integer | no |  | the modes to draw, numbered from 1; default the first max_modes |
| `panel_width` | number | no |  | width of a panel in drawing units, default 250 |
| `path` | text | yes |  | file to write, .svg or .png |
| `periodic_phase` | number | no |  | explicit Bloch phase advance per period in radians for repeat_periodic; auto-detected from result scalars or periodic BC when omitted |
| `periodic_shift` | list of 3 (number) | no |  | explicit [dx, dy, dz] lattice translation in metres for repeat_periodic; auto-detected from periodic faces or y-extent when omitted |
| `physics` | text: `dfs`, `ems`, `sms` | yes |  | the three-letter code of the physics |
| `plane` | text or object | no |  | slice plane for slice3d: "x"/"y"/"z", { "axis": "x"\|"y"\|"z", "at": metres }, or { "origin": [x, y, z], "normal": [nx, ny, nz] } |
| `repeat_periodic` | integer | no |  | number of unit-cell copies to tile along periodic_shift with Bloch phase exp(-i m periodic_phase) applied before expr, default 1 |
| `scale` | number | no |  | pixels per drawing unit of a PNG, default 2 |
| `show_edges` | text or true or false | no |  | 3D wireframe overlay: "regions" (default for field views), "mesh" (default for mesh3d, or true) or "none" (or false) |
| `warp` | number | no |  | largest motion of the deformed shape as a fraction of the window width; default 0.05 for a displacement (vector field `u`) and 0 for any other field; 0 draws the undeformed section |
| `window` | list of number | no |  | [x0, x1, z0, z1] in metres; default the extent of the modes |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"figure","args":{"path":"<path>","physics":"dfs"}}
```
