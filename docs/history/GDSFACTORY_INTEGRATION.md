# gdsfactory integration for nyx: design note (2026-11-20, rewritten plugin-first)

Question: how should nyx take layouts from gdsfactory? This is a design proposal, nothing is built. Confidence labels as in `MARKET_NOTE.md`: **sourced** (a page or file says it, read on 2026-11-20), **inferred**, **unchecked**. Items marked **checked in source** were read in this repository today.

## 1. Recommendation in one paragraph

Write a **gplugins plugin for nyx**, not an importer inside nyx. gdsfactory's `gplugins` has one thin plugin per solver (femwell, Tidy3D, Meep, Palace, Elmer, COMSOL, Ansys, Lumerical and others). Each takes a `Component`, a `LayerStack` and a cut line, builds the solver's input, runs it and returns the results as Python objects. A nyx plugin does the same: it slices the layout with gdsfactory's own code, writes a list of **existing nyx engine commands**, runs `nyx run --json`, and reads the replies back. nyx needs no new code for a first version: no new command, no new file format, no Python inside nyx. The earlier idea (a layout-slice JSON and an `import_layout` command) is dropped, because the command list already does that job.

## 2. How the existing plugins connect to a solver (sourced)

The pattern, read from `gplugins/femwell/mode_solver.py` and `gplugins/meshwell/get_meshwell_cross_section.py` (main branch, Nov 2026):

| Step | What femwell's plugin does | Note for nyx |
|---|---|---|
| Entry points | `compute_cross_section_modes(cross_section, layer_stack, wavelength, num_modes, ...)` builds a straight component, then calls `compute_component_slice_modes(component, xsection_bounds, layer_stack, ...)` | Copy both signatures; add `solver="nyx"`-style keywords later |
| Cut line | `xsection_bounds` = two (x, y) points in micrometres, made into a shapely `LineString` | Same |
| Which levels | Loops over `layer_stack.layers`; `layer_level.layer.get_shapes(component)` returns the shapes of each level | **gdsfactory evaluates derived layers (such as core = WG minus etches) itself.** Nothing is redone in nyx or in the plugin |
| Slice | Converts to shapely, intersects each polygon with the line, turns each segment into a rectangle (distance along the line, zmin) to (distance, zmin + thickness) | The reusable part. **It ignores `sidewall_angle`** (rectangles only), so nyx can do better with trapezoids |
| Padding | `add_padding_container` on a wafer layer, `wafer_padding` 2 um | Same idea, or a nyx window |
| Priority | Each surface carries `mesh_order` (from `layer_level.info["mesh_order"]` if present, else `layer_level.mesh_order`) and a name by layer or by material | nyx regions have painter order (section 4) |
| Materials | A hard-coded table (si 3.48, sio2 1.44, sin 2.0) the caller can replace with `material_name_to_index` | nyx has a real library and per-material TOML |
| Solve | In process: meshwell and gmsh make the mesh, femwell solves | nyx is called as a program |
| Result | A femwell `Modes` object | A small `Modes`-like Python object |

Other plugins drive external applications by a client library: the COMSOL plugin through MPh and the Ansys plugin through PyAEDT (both merged 12 Oct 2026, **sourced** from the commit messages), shipped as optional extras such as `gplugins[comsol]`. A plugin that launches a separate program is therefore the accepted pattern; nyx needs no special arrangement.

## 3. How a plugin drives nyx (checked in source)

- **`nyx run FILE --json`** runs a JSON-lines list of engine commands (`{"cmd": ..., "args": {...}}`), checked first, one reply line per command (`crates/app/src/runfile.rs`). `nyx check FILE --json` dry-runs the list and reports the plan without meshing or solving. Exit codes: 0 done, 1 a command failed, 2 the list is wrong.
- **`solve_run` waits for its job** inside a command list (`commandlist.rs`), so no polling is needed.
- **`result_info`** returns, per physics, the study, dofs and a list of solutions with `label`, `value`, `eigenvalue` and field names, plus `scalars`. An `ems` mode study labels its modes "mode N: neff = ..." and also writes scalars "neff N" (as text) and a Q figure (`modes_complex.rs`).
- **`result_field`** returns a named field of one solution. **`export_results`** writes the results file.
- **Geometry commands that already exist:** `Polygon` (corners in the x-z cross-section, with named edges), `Extrude` (with `taper_deg`, positive narrows towards the end), `Union`, `Difference`, `Intersect`, and `Overlay`, which names regions "in painter order (a later region wins where they overlap)" with kinds Material, Hole, Pml and Void (`geomdesc.rs`).
- **Materials:** `add_material` takes TOML (`extends = "silicon"` starts from a library material); `set_material` takes a name, or a name with `euler_deg` (Bunge angles), or scalar properties. **Crystal orientation is therefore already a nyx feature**, which a LayerStack cannot express.
- **Other routes, not needed first:** the MCP server (four script tools over a script kernel), and `.nys` scripts (Python syntax run by Monty).

## 4. Mapping, gplugins to nyx

| gdsfactory / gplugins | nyx command |
|---|---|
| Slice rectangle or trapezoid of one level, in (distance along the line, z) | `Polygon` with `pts` in (x, z). A sidewall angle is a trapezoid, so nyx can honour it where the meshwell path cannot |
| A level made of several segments or with holes | One `Polygon` per piece, `Union` / `Difference` where needed |
| `mesh_order`: which level wins where levels overlap | `Overlay` regions in painter order. **The direction of meshwell's `mesh_order` (lower wins, I believe) is unchecked**; reverse the sort if so |
| Background levels (substrate, cladding) | A Material region first in the overlay, so later regions win over it |
| Material name | `set_material` through a `materials_map`; an unmapped name is an error, never a guess |
| Wavelength, number of modes | the `ems` mode-study parameters (`describe_physics ems`) |
| Wafer padding | the window size; absorbing layers by `Pml` for leaky modes |
| Mesh resolution specs | `add_size` |
| Propagation direction | normal to the cut line in the layout plane; nyx's propagation is +y. The plugin must state which side is +y |

Units: gdsfactory works in micrometres, nyx in metres; the plugin converts once, at the edge.

## 5. Shape of the package

- **Name and place:** `gplugins-nyx` (or `gplugins/nyx/` with a `gplugins[nyx]` extra), a separate repository first. Not in this repo, because `*.py` is git-ignored here on purpose and nothing Python ships (ADR 13). Upstream to gplugins when stable, which also puts nyx beside femwell in the gdsfactory community. The gplugins licence is **unchecked**.
- **API:** the two femwell-shaped functions of section 2, returning a small object with `neff`, mode fields and, where nyx gives them, polarization and Q.
- **Finding nyx:** an environment variable (such as `NYX_BIN`) or the `PATH`. nyx needs its script worker next to it only for `.nys` scripts, not for command lists.
- **Dependencies:** gdsfactory, shapely, the standard library. No nyx crate, no solver library.
- **Inspectable:** the plugin keeps the command list it wrote, so a run can be rerun or opened without Python, and the list is the provenance record (it can also be named in the results file, see `RESULTS_FORMAT_STUDY.md`).
- **Contract:** the file the plugin writes is nyx's own engine command list, versioned by `engine_commands.v1.json`. No new contract line is needed for a first version.

## 6. Two tiers, because sweeps matter

- **Tier 1, polygons from a slice** (all of section 4): works for any layout and any cut. A sweep over a layout parameter reruns the plugin and writes a new list per point.
- **Tier 2, parametric strip, rib and slot:** the plugin writes nyx parameters (width, etch depth, slab) instead of frozen polygons, so nyx's sweeps, warm starts and the shared queue work. This is the agent-loop case; it comes after tier 1 works.
- **3D:** a window of a layout becomes `Extrude` of the polygons (with `taper_deg`) plus the overlay. gplugins has a separate `get_meshwell_3D` for this; reading it is a later step. Whether nyx's 3D mesher copes with arbitrary PIC layouts is **unchecked**.

## 7. Staging

1. **Spike, by hand, no package.** Slice the gdsfactory generic strip waveguide with shapely (or reuse `get_meshwell_cross_section` as is and convert its `PolySurface` polygons), write the command list, run `nyx check` and `nyx run --json`. Oracle: the existing SOI-wire n_eff golden, within its documented tolerance. Do not regenerate a golden (ADR 10). The spike answers the open items in section 8.
2. **The package**, tier 1, tested on the gdsfactory generic PDK stack, including a rib (two levels) and a case with `mesh_order` conflicts.
3. **Tier 2 parametric mapping and the way back:** write an optimised nyx polygon into a gdsfactory `Component`, so a design found in nyx re-enters the PDK flow.
4. **Upstream to gplugins.**
5. **Only if the spike shows a gap**, a nyx change: a numeric reply for n_eff, a quieter `--json` mode, or a command to reload a list. A Rust GDS reader (earlier option B) stays an option for GDS files with no gdsfactory behind them.

## 8. What was unchecked, and where it stands

| Item | Status |
|---|---|
| nyx can be run without a window from a command list | **checked in source**: `nyx run` and `nyx check` |
| `solve_run` waits inside a list | **checked in source** |
| Mode numbers reach the caller | **checked in source** as labels and a text scalar; **unchecked** whether `value` / `eigenvalue` hold n_eff as a number (the spike reads one reply) |
| Polygons with holes and layer priority | **checked in source** (`Difference`, `Overlay` painter order); **unchecked** that a rib or a layout with a hole meshes as intended |
| Crystal orientation | **checked in source**: `set_material` with `euler_deg`; the plugin must supply the angles (the wafer cut against the layout axes) |
| Derived-layer expressions and hierarchy | **resolved by reuse**: `layer_level.layer.get_shapes(component)` in the meshwell slicer, read today; nothing is re-evaluated |
| gdsfactory LayerStack serialisation | **no longer relevant**: the plugin reads the in-memory `LayerStack`, as femwell's does |
| Sidewall angle | **sourced**: the meshwell cross-section path ignores it (rectangles); nyx can use trapezoids |
| Direction of `mesh_order` | **unchecked**; test with a core that overlaps its cladding |
| Finding the `nyx` binary and the HDF5 DLL | **unchecked**; command lists and `result_info` replies need no results file, so the DLL is not on this path unless `export_results` is called |
| Mesh size control from gdsfactory data | **unchecked**; `add_size` exists, the plugin needs a rule for it |
| gdsfactory API drift (notebook APIs differ by release) | **unchecked**; pin a version range and run a nightly test, as other plugins do |
| gplugins licence and contribution rules | **unchecked** |
| The agent needs Python | still true; an agent without Python keeps the MCP route and option D (polygons passed inline) |

## 9. Risks

- Cut-line mapping: which side of the line is +y, and the angle between the cut and the silicon crystal axes. A wrong default would silently give wrong mechanical modes; require the plugin to take the crystal orientation explicitly.
- A GDS grid snaps coordinates to the database unit; a 1 nm grid is coarse next to a few-nanometre sidewall. Record the unit in the command list.
- Rings, couplers and mechanical resonators need a 3D window, not a slice; the 2D path covers cross-section mode problems first.
- Version drift of gdsfactory and gplugins (a changelog line and several recent renames show it moves).
- A plugin that shells out is harder to run on a hosted agent sandbox than one in process; nyx as a pure-Rust program with a static build (see `RESULTS_FORMAT_STUDY.md`) keeps that one binary.

## 10. Sources

- gplugins source tree (list of plugins, last commit 6 Nov 2026): https://github.com/gdsfactory/gplugins/tree/main/gplugins
- femwell plugin, `mode_solver.py`: https://github.com/gdsfactory/gplugins/blob/main/gplugins/femwell/mode_solver.py
- meshwell cross-section slicer: https://github.com/gdsfactory/gplugins/blob/main/gplugins/meshwell/get_meshwell_cross_section.py
- COMSOL and Ansys plugin commits (12 Oct 2026): linked from the tree listing above
- gdsfactory layers and layer stacks notebook: https://gdsfactory.github.io/gdsfactory/notebooks/03_layer_stack/
- gplugins overview (DeepWiki, indexed 2026-05-05): https://deepwiki.com/gdsfactory/gplugins
- nyx, checked in source: `crates/app/src/runfile.rs`, `crates/nyx-engine/src/{commands,commandlist,geomdesc}.rs`, `crates/nyx-ems/src/modes_complex.rs`; rules in `CLAUDE.md`, ADR 10, ADR 13, `docs/GEOMETRY_ARCHITECTURE.md`
