# 99. Figures of modes, and one node order in the results file

Date: 2026-11-20

## Status

Proposed (provisional number: 70 to 79 and 81 are reserved by other plans and the next free number on `main` must be checked at the merge). Built on bookmark `plotting-build`. Follows `docs/PLOTTING_PLAN.md` and `docs/RESULTS_FORMAT_STUDY.md`; amends ADR 62 (layout version 4).

## Context

A script written outside nyx (`plot_sos_mech_modes.py`) drew the modes of a result: `|u|` on the deformed section, arrows of the in-plane motion, the components of a mode, with the phase of each mode fixed by hand. Nothing in nyx could do this without a window, so an assistant or a script could not ask for a picture, and the file layout hid a trap that the script had to find by experiment: the cell fields are stored in the node order of nyx and `/viz/<physics>/cells` in the VTK order, with nothing in the file saying so.

## Decision

* **One node order in the file.** Layout version 4 writes, for each physics, `points` (the coordinates of the nodes) and `cells` (the node of every cell node) in the cell node order of nyx, and says it in the attribute `cell_order` (`v0 v1 v2 e01 e02 e12` for a P2 triangle, `v0 v1 v2 v3 e01 e02 e03 e12 e13 e23` for a P2 tetrahedron, `eij` the midpoint of the edge between vertices i and j). `fields/<name>` is indexed by node, `cell_fields/<name>/{re,im}[c * nloc + a]` is the value at node `cells[c, a]`. A reader gets the nodal complex field without any permutation. Version 3 files are still read. **The `/viz` group is gone** (Paul): it held only a copy of the nodes (identical to `points`) and the cells in another order, 1 to 8 % of a file on seven real files, and no vector copies. **The XDMF file stays** as the light index of the results as stored (Paul: it is for nyx to use the data efficiently, not for ParaView): its Geometry is `/results/<physics>/points`, its Topology `/results/<physics>/cells` with the `Information` element `CellNodeOrder` naming the order, its attributes the stored `fields/<name>` datasets (a vector with its `Components`, shape `[n, m]`). It says in a comment and in an `Information` element `Note` that a viewer such as ParaView must not open it, since it would read the mid-edge nodes in the VTK order and draw them wrongly.
* **Phase, sign and scale of a mode** are an analysis (`nyx_analysis::modeview`, ADR 56): the in-plane components (all but the propagation axis, y in the wafer frame) are rotated by `theta = -1/2 arg(sum u_c^2)` (plain square), then the largest in-plane component at the node of largest in-plane amplitude is made positive, and the peak of `|u|` over the core (the domains asked for, else all) is scaled to 1. The solver's own convention (largest component real, unit modulus) is not changed: it is in the goldens and in the nodal fields. The figure puts every mode in the phase above and reports `quadrature`, the share of the propagation component's energy that is still in `Re` after the phase step. For a guided mode, `d/dy` brings a factor `-iq`, so the propagation component is in quadrature with the in-plane ones when the section has a mirror symmetry along the propagation axis; the number says whether that holds for this mode, it is not assumed. For a real field (a real solver keeps one real amplitude per component) there is no phase and `quadrature` is null.
* **A section** (`nyx_analysis::section`) is the whole mesh of a 2D result, or the face of the cells of a 3D mesh that lies on a plane y = const (the smallest y unless `at` is given). The face of a tetrahedron is a P2 triangle in the native order; P2 triangles are drawn as four linear ones. No interpolation: a plane that no cell face lies on is refused.
* **A figure** is a list of drawing items (`nyx_analysis::canvas`), written as **SVG or PNG**. A field is not a vector drawing: it is rasterised into an image item that the SVG embeds as PNG, while the outline, arrows and text stay vector in the SVG. The PNG is drawn by a small software rasteriser (twice the size, averaged down; text with `ab_glyph` in the font the window carries), with the crates already in the lock file (`png`, `ab_glyph`, `epaint_default_fonts`, `base64`): no system library, no font to install.
* **The command `figure`** (engine, so the window, scripts and the MCP server share it): `physics`, `path` (`.svg` or `.png`), `kind` (`gallery`, `components`), `modes`, `columns`, `cmap`, `warp`, `arrows`, `window`, `at`, `core`, `family`, `panel_width`, `scale`. The reply carries, per mode, the angle the mode was turned by, whether it was flipped, its peak and `quadrature`.
* **No physics is named.** A vector field is found by the names of its cell fields (`<prefix>x`, `<prefix>y`, `<prefix>z`); the title of a mode comes from the label, the frequency scalars and the shares of energy the study wrote (ADR 80); the outline is the border of the section and the edges between cell labels; a displacement (prefix `u`) deforms the section, any other field does not.
* **Where it lives.** In `nyx-analysis` for now (Paul: "analysis is good for now"): no new crate, no new row of the layering table. If the drawing grows (sweep curves, the window drawing the same items), it moves to a crate of its own.

## Measured, and what was found wrong

* On three files of the benches (a 3D Bloch study of a strip, a q = 0 study with absorbing layers, an optical mode study) the figures come out as the script's. The propagation component of the 11 guided modes of the strip is in quadrature to `Re` share 4e-6 to 3e-3; the optical modes to 1e-6. The phase step turns the strip's modes by 0 to 16 degrees (solver convention against the figure's); one mode is turned by 90 degrees and flipped.
* The earlier note that `cell_fields` are "about five times duplicated" was wrong for these files: their meshes are one cell thick, a P2 field has 34980 values per component against 21861 nodes, a factor of 1.6, and the compressed `cell_fields` are 74 % of the file because they carry both parts of three components, not because of duplication. A node-indexed store of continuous cell fields would save at most that factor on these files and is not done.

## Not in this decision

* The window does not draw the figure yet (arrows, the phase of a mode on the screen, a phase animation): the same items can be drawn with egui later.
* Sweep curves (dispersion, mode tracking), the optical next to the mechanical in one figure, scalar results (temperature) and the axisymmetric half plane are not drawn; a result without a vector field is refused with the names it has.
* Nothing of nyx reads the XDMF file yet; it is written, and a test resolves every reference against the .h5. The export reply and the run records list it as before.
* A section inside a 3D mesh, a mesh not periodic along y, a quadrature check for modes with loss or without the mirror symmetry (it reports, it does not correct).

## Consequences

* `engine_commands` gains `figure` (additive, provisional 1.23). The results file is version 4.
* `nyx-analysis` depends on `png`, `ab_glyph`, `epaint_default_fonts` and `base64`.

## Code

`crates/nyx-analysis/src/{section,modeview,canvas,raster,figure}.rs`, `crates/nyx-engine/src/{figure,h5io}.rs`; tests `crates/nyx-engine/tests/{figure,h5}.rs` and the unit tests of the analysis crate.
