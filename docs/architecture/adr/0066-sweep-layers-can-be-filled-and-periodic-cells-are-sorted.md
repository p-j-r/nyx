# 66. A swept section can be a periodic cell: sorted vertices and `fill` layers

Date: 2026-11-14

## Status

Accepted (benchmark A3, bookmark `restructure-10-a3-snowflake`). Extends the sweep of ADR 54 and the engine mesher of ADR 57; uses the Floquet pairs of ADR 25.

## Context

A photonic or phononic crystal slab is a periodic cell patterned in the plane of the wafer, with the slab thickness and, for light, a layer of cladding above and below it. In nyx the only way to a 3D cell with Floquet walls that works is a 2D section swept along the surface normal, because the 3D lattice mesher fits curved and slanted walls only to order `alpha h` and does not make opposite faces congruent (`no vertex on the destination boundary matches ...`).

Two things were missing for a slab with a pattern.

1. `extrude_y` cuts each prism into tetrahedra along the vertex numbers. Opposite walls of a periodic section then are cut differently unless the vertices are numbered by their coordinates; `geomdrive::sort_nodes` does that, and the engine mesher did not call it for a sweep.
2. A sweep layer takes the section's domain numbers (plus an offset). A cladding layer above a slab with holes would be partly Si and partly air, so the cladding could not be a layer of air.

## Decision

* `SweepLayer` (project) and `mesher::Layer` (femcore) gain `fill`, a region name (`Option<String>` in the project, the domain number in femcore): every cell of the layer belongs to that region, whatever the region of the section; the domain offset is then not used. An unknown name is refused with the list of regions.
* The engine mesher numbers the section's vertices by coordinates (`sort_nodes`) before it sweeps, so the walls of a periodic section pair. Meshes of swept models that were not periodic change only in their vertex order.
* `set_mesh` takes `fill` in a `sweep` layer (`engine_commands` 1.6, additive).

## Not done

* A primitive hexagonal (rhombic) cell does not work: the 2D mesher leaves degenerate (sub-nanometre) boundary facets on slanted edges, so the two slanted walls cannot be paired. The benchmark uses the rectangular two-site supercell instead (the bands are folded once). Fixing the 2D mesher for slanted periodic walls is recorded in `ROUGH_SPOTS.md`.
* The 3D meshers (lattice, Rivara) still do not pair periodic faces.

## Consequences

* A slab-with-holes cell with cladding is a section, a sweep of `[slab, air, air]` layers (or `[air, slab, air]`), and a pair of Floquet conditions on the section walls; the same project gives the photonic (`ems` cavity) and the phononic (`sms`, air inactive) band structure.
* The vertex renumbering is a pure permutation: results do not change beyond rounding.

## Code

`crates/femcore/src/meshing/mesher.rs` (`Layer::fill`), `crates/nyx-engine/src/meshing.rs`, `crates/nyx-engine/src/project.rs` (`SweepLayer`), tests `a_fill_layer_gives_all_its_cells_one_domain` (`crates/femcore/tests/meshing.rs`) and `crates/nyx-engine/tests/swept_cell.rs` (a vacuum cell with a closed form, and the refusal of a fill that names no region).
