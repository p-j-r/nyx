# Mesher: what is still missing

Companion to `NYX_INTEGRATION_STATUS.md` (what exists). Written at the end of the mesher session on `mesher-work` / `mesher-work-next`.
Items are ordered by how much they matter to the solvers; "depends on" names what must exist first.

## 1. Depends on the new solvers (interface not visible from this repository)

* **Native mixed-cell solving.** `pml3::PmlBox::cells` and `meshing_pml::structured_pml_native` hold the uncut quads / prisms / hexes with
  nyx's tensor-product node order (quad 0,1,3,2 around the cell). Nothing consumes them yet; the solvers run on the simplex split.
  Needed from the solver side: the cell types and node orderings it wants, and whether boundary quads and triangles are both expected.
* **Curved and high-order elements.** No stored high-order nodes, no Jacobian validity checks for curved hex / prism / tet, no projection of
  mid-edge nodes onto the exact geometry. The 3D oracle (`geom3`) can project, but nothing places the nodes. Needs the node-storage
  convention and the maximum order.
* **Nédélec orders above 2 and the second family.** The mesher assumes nothing; sizing uses 12 / 6 points per wavelength for order 1 / 2 only.
  Points per wavelength for higher orders are unmeasured.

## 2. Estimators and loops

* **Lowest-order edge elements.** The recovered-Hessian loop does not beat uniform meshes at order 1 (1.2e-2 on 1490 tets against 5.3e-3 on
  696 for a cavity mode); it does at order 2 (2.6e-6 on 494 tets against 1.8e-4 on 696). Order 1 needs a different estimator: the
  Maxwell residual estimator (Beck, Hiptmair, Hoppe, Wohlmuth) or first-derivative metrics from the curl. Until then use physics sizing.
* **Goal-oriented (dual-weighted residual) adaptivity** with an eigenvalue quantity of interest; adjoint solves. Not started. The brief's
  claims about the adjoint of complex-symmetric PML matrices (forward solve with conjugated data; left eigenvector equals the right one)
  are still unverified numerically.
* **Field transfer** between meshes for time-domain use and warm starts of eigensolves. Not started.
* **Calibration of points per wavelength** beyond smooth modes (sharp corners, thin layers, high-contrast interfaces, Brillouin
  elastic modes). The defaults were measured on a PEC guide, a cavity and a high-index core only.
* **Brillouin and ring (axisymmetric) studies** are covered by the sizing rules on paper only; no test runs them with a sized mesh.

## 3. 3D adaptation limits

* **Surface flips and edge removal.** Only flat boundary patches with two tets on the edge flip; curved surfaces and material interfaces
  keep their triangulation, so surface quality there is whatever the start mesh and the splits give. Edge removal stops at 7 tets around
  an edge; boundary edge removal (more than two tets on a boundary edge) is missing.
* **Curved-surface collapse** refuses to turn a face by more than about 25 degrees. A finer-than-needed curved surface coarsens only where
  it is nearly flat at the scale of the edge.
* **Serial.** No independent-set threading, no deterministic parallel scheme. Closest-point queries are linear in the patch size (no
  bounding-volume tree); a patch with millions of triangles will be slow.
* **Mirror-symmetric models.** Frozen faces keep a periodic pair congruent (tested for a Floquet cell). A symmetry plane is the same
  mechanism but is untested; the result is only symmetric if the start mesh is.
* **Smooth analytic surfaces in the loop.** Sphere and signed-distance projectors exist and are unit-tested, but no end-to-end test adapts a
  curved body (a sphere or a ring) and checks convergence of a solver quantity.
* **No 3D geometry pipeline** in `geom`: 3D adaptation starts from a tet mesh (block, extrusion, import).

## 4. PML

* **Reflection level.** Only mesh convergence of damped modes (ignored study) and lossless equivalences (real stretch, unit stretch) are
  checked. The eigen solvers give no driven response, so the reflection of a lossy layer versus layers, strength and order has not
  been measured. The 2D scalar study in `results/m2_pml_reflection.md` is the only reflection data.
* **PML around non-box shapes** (round windows are refused; non-box 3D bodies are not supported).
* **Adapting after the layers exist.** The layers are rebuilt from the physical mesh each time; the flow adapt, then sweep, is the only one.
* **Native-versus-split agreement** for the vector solver cannot be checked until a native solver exists.

## 5. App

* The `physics_sizing` option is a checkbox on the model editor and is not shown anywhere else (no report of the sizes it used, no
  points-per-wavelength control, no per-domain override).
* `mesh_for_problem` is 2D; the 3D route (`physics_size::adapt_tets`, `pml3::add_pml_3d`, `estimate3`) has no app entry because there is
  no 3D geom model to mesh. The extruded-model path produces tets that could be adapted; that hook is not written.
* No progress or cancel hook inside the 3D adaptation passes (`progress::check()` / `stage()` are not called).

## 6. Process / hygiene

* The older brief's M0 items are still open: the four chip benchmarks with baseline numbers for the hand-tuned meshes, and the offline
  Mmg / Omega_h comparison (question 12 of the brief was never answered). **Dropped (paul, 2026-11-04, decision D1): the Mmg/Omega_h comparison is not going to be done.**
* Provenance and licence notes (`PROVENANCE.md`, `THIRD_PARTY_NOTICES.md`) cover the ported crate; the new files were written from the
  published formulas and the standard algorithms (Ericson for the triangle closest point; Dobrzynski-Frey style operators), with no
  third-party code read. No separate provenance entry was added for them.
* `FEMMESH_CHECK=1` validates the 3D working mesh after every pass, but the validation is not run in the default test run.
* The shared `nyx-jj` working copy was switched by another thread mid-session; further mesher work should use its own jj workspace.

## Update 2026-11-05 (bookmark `restructure-4m-mesh3d-features`, ADR 50 to 55)

* Section 2, **lowest-order edge elements**: done. The Maxwell residual estimator (`adapt::maxwell`) with the equidistributing metric reaches 5.2e-4 on 333 tetrahedra where uniform blocks give 5.3e-3 on 696 (order 1) and 8.0e-6 on 1402 where they give 6.7e-5 on 1470 (order 2). Goal-oriented adaptivity, field transfer and the calibration items remain.
* Section 3, **3D adaptation**: a 3D geom model is adapted with the exact faces as the oracle (ADR 51). Smooth analytic surfaces in the loop: tested end to end (sphere in a box, PEC sphere cavity). Still open: surface flips on curved patches, curved-surface collapse, boundary edge removal, threading, the bounding-volume tree; `adapt3` repeats a full incidence rebuild per sweep, which makes 17 000 strongly graded tetrahedra take about 30 s.
* Section 4, **PML**: structured layers around a box-shaped 3D window of a geom model (ADR 52), with the physical numbering kept; the adapt-then-sweep flow (P3) works. Round windows, non-box bodies and reflection studies remain.
* Section 5, **app**: `mesh_for_problem` is 2D by design; `mesh_for_sizing` and `request::mesh_request` (ADR 50, 54) are the 2D and 3D door for the engine. Not wired.
* The lattice mesher no longer makes zero-volume cells where a model plane lies on a lattice plane, and does not depend on the thread count (ADR 55).
