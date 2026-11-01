# Geometry in the nyx GUI: inventory and plan

Status (2026-11-01): approved. P0 is committed on branch `geometry-gui`; the legacy kinds are removed (see `crates/app/src/state.rs`); P1 is in `crates/app/src/geomui.rs` (with the live preview); boundary conditions take a selector rule (`Bc.on`) resolved against the built mesh, and the Materials tab shows region names. Open: click-to-rule picking in the viewport, name-keyed materials, 3D preview and the 3D driver (P4). Written against the `geom` crate in the wafer update (about 10.4k lines of Rust, identical to the uncommitted `crates/geom` in the main checkout) and `main` at `54b8829`.

## 1. Where things stand

* `geom` is **not committed anywhere**. It exists as uncommitted files in the main checkout (and as `nyx_wafer_update.tgz`), together with `femcore/src/geomdrive.rs` and edits to `app/{state,tabs,examples}.rs`, `problem.rs`, `mesher.rs`. Every branch and worktree is still at the old geometry. So the first step of any GUI work is to commit that update on a branch.
* What the wafer update already puts in the GUI: a "Model (geom)" button that creates a one-rectangle plate, an element size `h`, the sweep-along-y layers, a read-only table of region and boundary names with their integers, and a **JSON text box** holding the whole feature list. That is the only GUI route to the kernel today, and it re-serialises the list on every frame.
* Meshing a geom model works for 2D (x-z cross-section) and for that section swept along y. A 3D geom model is refused ("cannot be meshed yet"). The layer-stack driver is on the integration side, not in the GUI.

## 2. Capabilities in code that the GUI cannot reach

| Capability (in `geom` / `geomdrive`) | In the GUI today |
|---|---|
| 25 feature kinds: Rect, Circle, Ellipse, Polygon, RoundedRect, Box, Cylinder, Cone, Sphere, Torus, Revolve, Extrude (axis, taper), Transform, Pattern, Union, Difference, Intersect, Overlay, Pml, Cut, Periodic, Size, Tag | JSON only |
| Named parameters and expressions (`w/2+1e-6`), `with()` overrides | JSON only; no table, no sliders |
| Feature list editing: add, remove, reorder, suppress, rename, tags, `origin: Code` grouping | none |
| Overlay: painter order, region kind (Material / Hole / Pml / Void), pinned ints, naming scheme (Feature / Region / Numbered) | JSON only |
| Selectors: the grammar (name, tag, region, interface, outer, plane, normal, extreme, in_box, kind, radius, area, adjacent, boolean algebra), `resolve`, `suggest` | none; boundary conditions still take integers |
| Name-keyed use: `MeshedGeom::domains("sapphire*")`, `boundary("rib\|air")`, `select(text)`, `set_material`, `pml_materials` | code only; `ProblemSpec` has no name-keyed fields, Materials and Physics tabs list integers |
| Preview of the exact geometry (`tessellate2d`, `tessellate`) with face ids | none; the viewport shows only a mesh |
| Pml, Size, Periodic, Cut as data (thickness, layers, grading, h, max_dev, congruence check, symmetry) | JSON only |
| Warnings: selector that resolves to nothing, grazing contact, refusals with feature names | shown only as one red line in "Names and numbers" |
| `to_rust()` (feature list to builder code), `topology_signature`, `diff` | none |
| Legacy lift (`Shapes`, `Solids`, `Rectangle`, `Block` to a feature list) | not written yet |

## 3. Design principles

1. **Code and GUI are two complete routes to one `FeatureList`.** The GUI edits exactly that value and nothing else. Everything the builder can make the GUI can make, and "Show as code" prints it back (`to_rust`). A test builds a model both ways and compares the lists.
2. **No numbers in the interface.** Regions and boundaries are addressed by name or by a selector rule. The integers stay visible, read-only, in a Names panel.
3. **Lean.** One generic form renderer driven by a descriptor table (a few lines per feature kind), not 25 bespoke editors. One selector field widget used everywhere. Estimated 1.5 to 2k lines of app code in one new module; the two legacy editors are left alone and retired later.
4. **Fast.** Evaluation is measured at about 1 ms for a 50-feature model and 6 us for a selector, so the model is re-evaluated on every edit (debounced one frame) on the UI thread, with the result cached by content hash. Meshing and solving stay behind the Build button on the worker thread, as now. The preview never meshes.

## 4. The interface

One Geometry tab, built from four small pieces. It becomes the Concept section if the workflows branch lands.

**Parameters.** A table of `name = expression` with the evaluated value beside it, a drag-to-scrub handle on numeric ones, `+` to add. Features refer to them by name. This is also where a sweep will later attach.

**Features.** An ordered list, one row per feature: icon, name, one-line summary (`Rect 500 nm x 300 nm`), eye (suppress), warning mark if it failed. Click opens its form inline. Forms are generated from the descriptor table: expression fields, a feature picker for inputs, a selector field, choices for enums. `Add` opens one menu grouped as Profiles, Solids, Booleans, Transform and pattern, Simulation (Pml, Cut, Periodic, Size, Tag). Drag to reorder (references only point earlier, so an invalid drop is refused with a message). Features produced by code are grouped and editable through their parameters only, so the next code run does not overwrite hand edits.

**Regions.** The Overlay as an ordered list (later wins): name, source feature, kind, optional material. Naming scheme is one dropdown. Materials are assigned here by name, replacing the per-integer Materials list for geom models.

**Names.** A read-only panel: each region and boundary with its integer, element count after meshing, and the rule that produced it. Zero-hit rules are red.

**Viewport.** A new content mode, "Geometry", beside "Mesh": the exact regions coloured by region, edges drawn, names on hover. Click a face or edge and a small popover lists the simplest selectors that resolve to exactly that set (`suggest`: auto-name, tag, name plus predicate, geometric predicate), each with "resolves to N". Choosing one fills the selector field that asked for it (a boundary condition, a Size rule, a Tag, a Periodic side). Shift-click grows the pick; one-click growths are "same normal", "same interface", "all faces of this feature". Rules are stored, never ids.

**Selector field.** One widget: a text box with live parse, live count, red on error or zero hits, and a Pick button. It is used by Physics (boundary conditions), Size, Tag, Periodic and Pml.

**Physics and Materials for geom models.** Boundary conditions gain `on: <selector text>` next to the integer `labels`, materials gain region names; one `resolve(&mut ProblemSpec, &MeshedGeom)` before `solve` fills the integer fields the solvers already read, so no solver changes. Legacy projects keep their integers.

**Code route.** A "Show as code" button (copy, no editing in place) and "Load model file" for a feature list written by a program (JSON or `.nyx.xml`), with reload on file change. The same file is what the MCP tools and the future workflow store read.

## 5. Phases

| Phase | Content | Done when |
|---|---|---|
| P0 | Commit the wafer update (`geom`, `geomdrive`, app and femcore edits, goldens) on a new branch off `main`; build, test, clippy. No UI change. | `cargo test` and clippy green; JSON box still works |
| P1 | Parameters, Features list with descriptor forms, Regions, live evaluation with inline errors, 2D geometry preview, Show as code. Replace the JSON box. | The 300 nm SOI strip is built entirely in the GUI and equals the code-built `FeatureList`; evaluation stays under a frame |
| P2 | Selector field, pick-to-rule in the viewport, Names panel, warnings. | Click a rib wall, get `name(rib.left)` style rules; zero-hit rule turns red |
| P3 | Name-keyed Materials and boundary conditions, `resolve` before solve, examples SOI, ring chips and Bragg cross-section re-expressed as geom models with no typed integers. New projects default to a geom model; the legacy editors appear only for old projects. | Examples give the same n_eff within the documented mesh band (SOI TE 2.67179) |
| P4 | 3D preview (`tessellate`), Pml / Periodic / Cut drawn in the viewport, 3D meshing once the integration side delivers the layer-stack driver. | Bragg cell shown in 3D with its periodic faces highlighted |

Out of scope here: parameter sweeps and mesh morphing (the workflows F11 item), curved elements, STEP import, retiring the legacy editors.

## 6. Decisions I need

1. **P0 ownership.** The geom update is uncommitted in your main checkout. I would copy it from the tarball (identical, checked by diff) onto a new branch `geometry-gui` and commit it there, leaving your checkout untouched. Is that right, or is another session about to commit it?
2. **Default for new projects.** Recommend: geom model becomes the default in P3; Shapes, Solids, Rectangle and Block stay available but secondary. Fine?
3. **Relation to the workflows branch.** Recommend: all new UI in one module (`geomui.rs`) with a single `show(ui, &mut Project)` entry, so the workflows shell can mount it as Concept without edits to it.

## 7. Risks

* Reorder and suppress can break references or the topology; the evaluator reports it by feature name, and the list shows it on the row instead of blocking edits.
* Editing a feature list is only as fast as `eval` on large models (hundreds of holes): measured 1 ms for 200 holes, but the preview tessellation must be cached per evaluation and rebuilt only on change.
* 3D models can be edited and previewed before they can be meshed. The Build button says so in words.
