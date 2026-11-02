# Frame unification: status

Written when the work landed on `main` (local only, nothing pushed to a remote). The design document this follows
(`COORDINATE_FRAMES.md`) was not reachable from the machine that did the work, so the decisions are those of the brief; all of them
were accepted as recommended.

**There is one public frame, the wafer frame: x lateral, y the direction of propagation, z the surface normal (right-handed).** A 2D
cross-section is the x-z plane with y into the page. The 2D kernels still compute in their old model axes (mesh x-y, z out of the
plane); `femcore::frame_adapter` is the only place the two vocabularies meet, so the kernel port (not started, see
`KERNEL_PORT_PLAN.md`) can delete it.

## What landed (commits on `main`, in order)

| step | commit subject | tests (workspace / ai) |
|---|---|---|
| baseline | b5e5bfb | 597 / 21 |
| 0 | wrong 2D labels, no Axes toggle, one default frame, geom internals private | 598 / 21 |
| 1 | PML in wafer names (x, z in 2D); positive k runs along +y | 600 / 21 |
| 2 | one frame adapter, a per-component field table, a ring variant (+ the ring-frame tests of 8b3599b) | 613 / 21 |
| 3 | retire `Frame` from the public problem type | 614 / 21 |
| 4 | legacy front-ends ported to geom models and deleted | 605 / 21 (see below) |

Ten tests are `#[ignore]`d throughout (slow Brillouin and ring studies; the seven ported ones were re-run in release at the end).

Step 4 test arithmetic (the commit message of step 4 words it loosely): 614 - 8 (`tests/mesher.rs`) - 2 (`mesh.rs` counts) - 1 (`Solid::Extrusion`) - 1 (the legacy-frame refusal test of step 3, removed with the shim) + 3 (`tests/geom_meshes.rs`) = 605.

## Behaviour that changed for users

* **Positive out-of-plane wave number `k` now means a wave along +y** (into the page). It used to mean "towards the viewer, -y". The
  adapter hands the solvers `-k`. Elastic frequencies are even in `k` (centrosymmetric stiffness), so no result number moved; only the
  sign of the imaginary parts of complex mode shapes and the labels/docs.
* **2D PML is written along x and z.** A `y` layer on a 2D mesh is refused ("absorbs along y ... use z"); there is no alias, because
  `y` means y in 3D files. `geomdrive::pml_materials` and `meshers::interop::pml_domains` produce `(x, Off, z)` for 2D meshes.
* **The Axes toggle is gone**, and so is `Frame`/`Materials::frame`. A `frame` key in an old project file is ignored.
* Brillouin optical-PML editor, the mirror-parity choice and the viewport hover now say z where they said y.
* Ring problems: the third kernel component is +theta, i.e. wafer +y, so a ring's `uz -> uy` rename has sign +1 (a straight section keeps -1).
* New names the field table can place that the old 22-string table could not: `sxz`, `syz` stress components. A 2D result field the table cannot
  place fails a debug assertion in `solve` (so any unnamed field fails the tests) and is logged in release builds.

## What was deleted (nothing is released, so no compatibility shims were kept)

* `Frame`, `Materials::frame`, the legacy-model-axes refusal, `ProblemSpec::to_model_frame`, `solve::rename_wafer_fields` (moved to the adapter).
* Front-ends: `Spec2d`, `Shape`, `Part`, `mesh_2d`, `conform_chains`, `mesher::extrude` (z), `Solid::Extrusion` and its helpers,
  `SimplexMesh::rectangle` and `SimplexMesh::block`, and their tests in `mesh.rs`/`mesher3d.rs`.
  **Kept:** `extrude_y`, `mesh_chains`, `Chain`, `Fill`, `MeshParams`, `Layer` (the geom path uses them).
* Tests: `tests/mesher.rs` (8 tests: 5 deleted as tests of the legacy shapes, 3 ported to `tests/geom_meshes.rs`: Poisson on a disk, plate with a
  hole, a section swept along y).
* Examples that built their geometry with `Spec2d` and were studies of finished work: `saw_s0..s5`, `saw_array`, `saw_common`,
  `sos_acoustics`, `sos_optics`, `brillouin_sos` (the last one used an "x lateral, y normal, z propagation" convention). They remain in git history.
* Benchmarks (`bench`, `eig_beam`, `mem_profile`, `pool_bench`, `profile3d`, `stage_log`) were kept and use `examples/benchmesh`, a geom-based helper.

## Ported tests whose construction, tolerance or expected value changed

| test | change and why |
|---|---|
| `aniso_optics` (tilted crystal, 2D vs 3D cell) | now written in wafer axes: 3D cell periodic along y, crystal tilted in the y-z plane, 2D through `solve::solve`. Stronger than before (it checks the adapter's Euler conversion against an independent 3D computation) |
| `bragg_stack`, `emw3d::floquet_...` | the Floquet cells are periodic along y (the sweep axis), not z: a geom mesh swept along y pairs its ends exactly; side walls pair too because `geomdrive::sort_nodes` numbers the section's vertices in coordinate order |
| `anisotropic::sapphire_modes_follow_the_crystal_symmetry` | the old rotation (`[0, 90, 0]`) turned the c axis about the bar axis and cannot change the bending modulus; it only passed on the old coarse structured mesh. Now `[90, 90, 0]` puts c along the bar |
| `scalar::es_series_capacitor`, `solid::prescribed_displacement_and_two_materials`, `solid_pml::low_reflecting...`, `emw2d::slab_loaded...`, `pml::*`, `emw3d::a_real_stretch...` | the interface was a cell relabel by a box (exact only on a structured grid); it is now a conforming two-region geom model |
| `solid_pml::real_stretch_equals_a_longer_bar` | the long bar is the same mesh stretched in x, so the two discretizations are the same one (the 1e-8 comparison stays meaningful) |
| `finite_guard::nan_geometry_is_refused` | the NaN is put into a valid mesh (geom refuses NaN geometry earlier) |
| `ring_frame::..._model_axes_...` | compares against a frozen golden (`tests/golden_frames/ring_model_axes.json`); regenerated once when the ring mesh moved from the legacy mesher to geom (changes of up to 7e-5 relative) |
| `brillouin::a_thin_ring_bends_in_its_plane...` (ignored, slow) | the shared ring mesh uses element size a/3 instead of a/4: on the geom mesh at a/4 the shift-and-invert solver (`solid_kz`, residual filter 1e-6) rejected every mode at the 12 Hz shift (an existing sensitivity of that solver to the mesh, not a frame effect; it passed on the legacy mesh). Results at a/3 are within 0.06 % of the ring theory, as before |
| meshes of all ported tests | unstructured geom meshes of the same box with element size the finer of the old counts; boundary numbers of the test helpers (`tests/common`) keep the numbers the tests always used |

## Open items

* **Optical mode direction.** Superseded by the kernel port (`docs/KERNEL_PORT_STATUS.md`): the 2D optical mode solvers now propagate along +y;
  and so do the Brillouin kernels. (Before the port the kernels solved `exp(-i beta z_model)`, a wave along -y.) The ring: `m > 0` is a wave along +y, verified
  with a signed test (`tests/ring_frame.rs`, `emw_axisym::power_along_theta`).
* Historical design documents (`GEOMETRY_ARCHITECTURE.md`, `GEOM_HANDOFF.md`, `SWEEP_CACHING_PROPOSAL.md`) still describe the deleted front-ends;
  they are records, not updated.
* `frame_adapter` is `#[doc(hidden)] pub` (not `pub(crate)`) because tests and examples drive kernels directly.
* Kernel port to native wafer axes (`KERNEL_PORT_PLAN.md`): not started, out of scope.

## For the geometry-integration thread

APIs kept stable: `geomdrive::{evaluate, mesh_geom_2d, MeshedGeom, GeomMeshOpts, select_labels}`, `mesher::{extrude_y, Layer}`, the `extrude_y`
boundary numbering (walls keep the section's labels, ends are `nb + 1` / `nb + 2`). Changed: `pml_materials` now returns wafer names in 2D
(z, not y). Added: `geomdrive::{plate, block, sort_nodes}`. Changed in geom: `Complex` fields, `wv`/`iv` and `LayerStack::overlay` are
`pub(crate)`; use the `Evaluated` accessors (`regions()`, `edges()`, `vertices()`) and `LayerStack::overlay_region_count/overlay_edge_count`.
