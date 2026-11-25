# `set_mesh`

**Set the mesh settings** · access `edit`: changes the project

Give only what changes. `adapt_scale`, `adapt_grade`, `adapt_max_cells`, `gap_cells`, `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, `max_steiner_multiplier`, `singular_edges`, `quarter_point_xi`, `pml_mode`, `structured_pml`, or `sweep` set to null return to their defaults. `sweep` extrudes a 2D cross-section or a 3D geometry model along y through layered intervals with optional `fill` or per-layer `regions` remap.

## Arguments

| Argument | Type | Required | Default | What it is |
|---|---|---|---|---|
| `adapt_grade` | number or null | no |  | adaptive / AMR metric gradation slope (> 0); null uses 0.4 |
| `adapt_max_cells` | integer or null | no |  | adaptive / AMR maximum cell count (>= 1); null uses 2000000 |
| `adapt_scale` | number or null | no |  | adaptive meshing: the factor on the sizes the model asks for (smaller is finer); null is the simple mesher |
| `corner_grading` | true or false | no |  | 3D sheet stack: radical Morrison-Lewis corner grading toward sharp 2D sheet corners, default false |
| `curved` | true or false | no |  | 3D: curved (quadratic) isoparametric geometry for cells on curved surfaces and fitted sheet arc rims (ems, dfs, and sheet3d), default false |
| `dz_min` | number or null | no |  | 3D sheet stack: floor on initial vertical z-slab thickness near sheet planes in metres; null disables the floor |
| `gap_cells` | number or null | no |  | 3D sheet stack: target minimum number of 2D triangles across narrow coplanar gaps between conductor sheets (positive number, e.g. 2.5); null disables automatic gap clamping |
| `h` | number | no |  | global element size in metres |
| `max_cells` | integer or null | no |  | maximum number of mesh cells allowed before aborting meshing (overrides the default 2,000,000 cap on adaptive meshing and enforces pre-allocation guard on 3D sheet stacks); null restores default |
| `max_dev` | number or null | no |  | largest chord deviation from a curved edge or face, metres; null: h / 50 for the adaptive mesher and no limit for the simple 3D mesher |
| `max_steiner_multiplier` | integer or null | no |  | 2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial PSLG vertices; null uses default max(20 * n_init, 100000) |
| `max_z_slabs` | integer or null | no |  | 3D sheet stack: maximum total number of horizontal z-slabs across active stack layers (must be >= active stack layers); null disables the cap |
| `physics_sizing` | true or false | no |  | size the elements from the waves of the physics (points per wavelength) |
| `pml_mode` | text or null: `auto`, `structured`, `inline`, `null` | no |  | how absorbing layers are meshed: "auto" (default, structured native Quad/Prism/Hex layers when PML regions exist), "structured", or "inline" (legacy unstructured simplices) |
| `points_per_wavelength` | number or null | no |  | points per shortest wavelength when physics_sizing is on; null: 12 for order 1, 6 for order 2 |
| `quarter_point_xi` | number or null | no |  | 3D sheet stack: Barsoum quarter-point mid-edge parameter xi_QP in [0.25, 0.50]; null uses 0.25 |
| `singular_edges` | text or null: `none`, `quarter_point`, `null` | no |  | 3D sheet stack: edge singularity enrichment on zero-thickness conductor rims ("none" or "quarter_point"); null uses study defaults |
| `strict_mesh` | true or false | no |  | fail meshing if 2D Ruppert refinement terminates early before completion or drops unhealed collinear sliver triangles, default false |
| `structured_pml` | true or false or null | no |  | legacy boolean alias for pml_mode (true -> "structured", false -> "inline") |
| `sweep` | list of object or null | no |  | layers to sweep a 2D cross-section or 3D geometry model along y |
| `sweep[].domain_offset` | integer | no |  | added to the domain numbers of the layer's cells |
| `sweep[].fill` | text | no |  | region name: every cell of the layer belongs to it |
| `sweep[].n` | integer | yes |  | cells across the thickness |
| `sweep[].regions` | object | no |  | optional per-layer region remap { source_region: target_region } in this layer |
| `sweep[].thickness` | number | yes |  | thickness of the layer along y, metres |
| `z_coarsen` | true or false | no |  | 3D sheet stack: progressively coarsen the 2D sheet triangulation away from the conductor plane along z via conforming pure-tetrahedral transition slabs, default false |
| `z_growth_ratio` | number or null | no |  | 3D sheet stack: geometric growth ratio (> 1.0) for vertical z-slab thickness away from sheet planes; null uses default 1.35 |

Lengths are in metres and every unit is SI unless the argument says otherwise.

## Example

A line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):

```json
{"cmd":"set_mesh","args":{}}
```
