# femmesh status

Lives in `meshers/femmesh` on `main`. Development history: branch `claude/quirky-euler-20hyv1` (commits `c9486e2` to `9c01565`). Takeover notes: `../../HANDOFF.md`; provenance: `../../PROVENANCE.md`. Build and test: `cargo test --release` (86 tests, ~15 s).
Machine for all numbers below: 4-core Intel Xeon @ 2.8 GHz VM, rustc 1.97.0, release profile
(`opt-level 3`, fat LTO), no `target-cpu=native`. VM timing noise is about ±15 %, so differences under
10 % are ties.

## What was and was not available

The nyx repository, the plan (`nyx_meshing_adaptivity_plan.md`) and the nyx source files named in the
handoff were **not reachable** in this session. paul's instruction was to build the port as a
standalone crate (`femmesh`) that another session moves into `femcore`. The mesh types in
`mesh.rs` / `cellmesh.rs` are therefore **stand-ins** for `SimplexMesh` / `CellMesh`, written from the
handoff's description, and the label scheme for PML regions is configurable instead of copied from
`chips.rs`. Nothing here was checked against the real nyx APIs.

## Milestones

| # | Milestone | State | Notes |
|---|---|---|---|
| M0 | fixtures, statistics, VTK, timing/determinism patterns | **partly done** | Done: box meshes (2D/3D) with boundary labels, analytic fields with exact Hessians, quality and edge-length statistics, `.fmsh` compact lossless format with results (no VTK), per-phase timers, deterministic chunked parallel reductions. **Blocked on nyx:** the four chip benchmarks and the baseline numbers of the hand-tuned chip meshes. **Not done (waiting for answer 12):** Mmg/Omega_h offline export. |
| M1 | `adapt/metric.rs` | **done** | Tier-A algebra tests pass; fused derivative-to-metric kernel; benchmark below. |
| M2 | `adapt/pml.rs` in 2D | **done except the chip comparison** | Generator, labels, split fallback, reflection-vs-layers study done on a stand-in Helmholtz solver. **Stop and report** (this document). The chip eigenvalue comparison needs nyx. |
| M7 | `adapt/pml.rs` in 3D | **done for what can be checked without nyx** | Prisms on faces, hexahedra on edges and corners, conforming, with the simplex split. The scalar harness now runs P1 tets, trilinear hexes and P1xlinear prisms: native and split eigenvalues of the layered mesh converge at second order and agree to 1.2 % at the finest level; 3D PML works for both (`results/m7_pml_3d.md`). **Not done:** the vector `emw3d` cavity check (needs nyx), and the HDF5/XDMF results format of the original wording (see decisions). |
| M3 | `adapt/topo.rs` + `adapt/geom.rs` (2D) | **done** | Classification from cell and facet labels (interior / chain / corner, chains, frozen chains), validity checks, round trip to `SimplexMesh`, polyline geometry from the mesh itself plus analytic projectors (circle, signed distance) for curved boundaries. |
| M4 | 2D operators and driver | **done, stop and report** | split, collapse, flip, smooth with quality gates; Hessian-to-metric-to-mesh end to end; results below. Not compared with hand-tuned nyx meshes (not available). |
| M5, M6, M8–M16 | | not started | |

## M1 result: fused derivative-to-metric kernel

`results/m1_metric_bench.md`, 200 000 vertices, minimum of 7 interleaved runs, 4 threads for the
parallel column:

| | unfused | fused 1 thread | fused 4 threads | complexity miss fused / unfused |
|---|---|---|---|---|
| 2D | 25.8 ms | 26.5 ms | 13.7 ms | 3e-15 / 1e-8 |
| 3D | 196 ms | 134 ms | 40.5 ms | 2e-16 / 1e-5 |

* 3D: 1.46x faster than the multi-pass version and reaches the target complexity exactly. 2D: a tie.
* "Unfused" is the straightforward multi-pass version (separate loops, `|H|`, determinants, scaling,
  complexity, clamp). It scales once and clamps once, so it **misses the target complexity whenever
  clamping is active** (1e-8 in 2D and 1e-5 in 3D on this data). The fused kernel solves for the
  global scale on the *clamped* complexity, so it hits the target to roundoff. The two are therefore not
  equal work; the speed comparison understates the accuracy gain.
* Diagnosis of the first version: it bisected the scale (about 35 evaluations over all vertices) and was
  0.64x in 2D. Replacing the bisection with a safeguarded regula falsi on the monotone log-log function
  brought the scale search to a few evaluations.
* 3D time is 66 % eigen-decomposition (Jacobi) — hypothesis: a closed-form 3x3 solver would help, at some
  cost in eigenvector accuracy. Not tried. 2D is dominated by `powf` and the scale search, not by the eigen-solve (15 %).
* Results are bit-identical for 1, 2, 3 and 8 threads (fixed 512-element chunks, sums in chunk order).
  Assembly of the generated code was not inspected.

## M2 result: structured quad PML vs the same layers split into triangles

`results/m2_pml_reflection.md`, produced by `examples/pml_reflection.rs`. Independent scalar Helmholtz
solver in `validate/` (P1 triangles, Q1 quads, complex coordinate stretching sampled at the quadrature points,
Dirichlet outer boundary, banded LU). Physical domain 3x3 wavelengths, 12 points per wavelength,
Gaussian source, layer thickness 1 wavelength. Error = relative L2 difference over the physical vertices from a
reference on the same interface mesh with a 2-wavelength layer and 48 layers. Quad layers are compared with a quad reference
and split layers with a split reference (the two references differ by 2.6e-3 on the physical vertices — comparing both with a single
reference favoured whichever kind the reference was).

Best error over the design reflection R0 in 1e-1..1e-8 and grading in {1, 4}:

| layers | dof | quad | triangles (split) |
|---|---|---|---|
| 3 | 1849 | 1.0e-1 | 9.8e-2 |
| 4 | 2025 | 4.3e-2 | 4.3e-2 |
| 6 | 2401 | 4.6e-3 | 4.5e-3 |
| 8 | 2809 | 1.2e-3 | 1.8e-3 |
| 12 | 3721 | 1.2e-3 | 1.2e-3 |

* **For scalar P1/Q1 the structured quad layers are not more accurate than the same layers split into
  triangles.** The split has the *same vertices*, hence the same dof count; differences are within the noise of
  the two references. About 6 layers give ~0.5 %, about 8 give ~0.1–0.2 %, and the error floors near 1e-3.
* So the case for quads in 2D rests on what this test cannot show: half the cell count, exact control of layer
  count and grading, no sliver triangles from an unstructured layer mesher, and (hypothesis) Nédélec on quads. I could not compare
  against nyx's current triangular PML mesh because I do not have it.
* The handoff's Tier-A line "native mixed PML mesh and its simplex split give the same eigenvalues" cannot hold
  exactly: Q1/Nédélec on quads and P1/Nédélec on the split triangles are different finite-element spaces. I expect agreement to
  discretisation error only (the table above is that agreement for the scalar problem) and propose the Tier-A test be
  stated that way.

PML generator properties tested (2D and 3D): conforming (every face shared by at most two cells,
the only boundary faces are the outer surface), positive scaled Jacobians, exact total area/volume, deterministic bits,
region boxes with the correct thickness and the stretched directions encoded in the label, materials
crossing the interface continue into PML copies, partial layers (sides with zero thickness), graded interface meshes,
rejection of non-box domains, and the simplex split being conforming and volume-exact (hex into 6 tets, prism into 3,
cone from the smallest vertex id, quad diagonal through the smallest id — this is what makes neighbouring cells agree).

## M3/M4 result: 2D metric-based adaptation

Operators (`adapt/ops.rs`) and driver (`adapt/drive.rs`), written from the standard published rules
(length thresholds 1/sqrt 2 and sqrt 2, link condition for collapses, quality-gated flips, length-equalising
smoothing); no third-party code was read. Everything is serial and deterministic (bit-identical reruns tested).
Constrained features come from the labels: boundaries, material interfaces (cells of different label) and corners are
preserved exactly (tested: interface covered edge for edge, side lengths per label unchanged to 1e-12, corners kept,
L-shaped domain keeps its re-entrant corner, a frozen side keeps its vertices bit for bit).

Independent checks (`tests/adapt2d.rs`, 11 tests): uniform, coarsening, graded isotropic, anisotropic boundary layer
(element count matches the complexity integral to within 40 %), material interface, curved boundary, frozen boundary,
L-shape, determinism, Hessian-driven adaptation.

**Error vs elements** (`results/m4_adapt_convergence.md`, `examples/adapt_convergence.rs`; P1 interpolation error recomputed by
degree-5 quadrature on the output mesh; the comparison mesh is a *uniform structured* mesh with the same number of triangles):

| function | elements | adapted L2 | uniform L2 | gain | q_min | edges in band |
|---|---|---|---|---|---|---|
| anisotropic Gaussian | 8 k / 32 k / 128 k | 3.1e-5 / 8.0e-6 / 2.0e-6 | 2.6e-3 / 7.3e-4 / 1.9e-4 | 85x / 91x / 95x | 0.58 / 0.41 / 0.47 | 99.5 % |
| oblique layer (width 0.01) | 33 k / 95 k / 231 k | 1.5e-6 / 1.0e-6 / 9.0e-7 | 9.4e-3 / 3.4e-3 / 1.4e-3 | 6400x / 3300x / 1600x | 0.46 / 0.42 / 0.49 | 98-99.6 % |
| thin ridge (width 0.005) | 18 k / 53 k / 192 k | 2.0e-5 / 5.3e-6 / 2.1e-6 | 2.9e-2 / 8.8e-3 / 3.1e-3 | 1500x / 1700x / 1400x | 0.49 / 0.47 / 0.60 | 99.3-99.5 % |

(the rows are the graded p=2 variant; the full table has p=inf and ungraded variants too). The adapted error scales as 1/N
(ratio 3.9 and 4.0 per 4x elements), which is the optimal P1 rate; the large gains for the layers are mostly the weakness of the uniform
baseline, which cannot resolve a layer of width 0.01 at these sizes. These are interpolation errors, not solver errors.

**Speed** (serial, this VM): about 15 k triangles/s end to end on the Gaussian (144 k triangles in 9.4 s; phases split 3.7 s, flip 1.5 s,
collapse 1.3 s, smooth 2.4 s for the 128 k-target graded run, 3 iterations). I did not compare with Mmg or Omega_h (question 12 is unanswered); no claim about being
competitive. The first iteration overshoots (splits far more than the final mesh keeps, then collapses), and the passes recompute candidate
lists from scratch, so there is clear room to improve before the threaded milestone (M13).

**Robustness: what I found and fixed, and what remains**

* Fixed: (1) 30 split passes back to back degraded slivers (quality halves per pass), flips between split passes cure it (q_min 0.000 -> 0.41);
  (2) the Hessian metric produced NaN for exactly flat regions (determinant underflow, per-vertex floor): now a global floor and a log-space
  determinant; (3) interior vertices sitting almost on a boundary edge (caps, q 5e-8): smoothing now also tries the ring centroid; (4) edges on a
  frozen boundary that are too long for the metric cannot be split, and their neighbours were split forever: such neighbours are now left alone and
  the count is reported (`unresolved_frozen`); (5) the same flip candidates were re-tested every pass: only edges next to changed triangles now
  (identical output, 2.5x faster flips); the snapshot is linear time.
* **Remaining limitation, important:** adaptation to an *ungraded* metric with sharp features is unreliable. 6 of the 30 ungraded runs in the study do
  not satisfy the stopping rule and one (thin ridge, p=2, 8000) runs away and is stopped by `max_cells`. Every one of the 15 graded runs converged with
  q_min >= 0.28. Use `BackgroundMetric2::new_graded` / `metric::grade` (slope 0.4 used in the study); the cost is more elements than the target
  complexity (1.1x to 3x here).
* **Remaining limitation:** a feature narrower than the starting-mesh spacing can be missed (the ungraded ridge collapsed to 8 triangles from a 7x7 start). Edges are
  also measured through their midpoint in the first two iterations, which helps but is not a guarantee; gradation widens the fine region and
  removes the problem in the study. In the solve-estimate-adapt loop the start mesh is the previous adapted mesh, so this matters mostly for the first cycle.
* The `converged` flag in `AdaptReport` says whether the stopping rule (splits plus collapses under 0.2 % of the triangles) was met; a run that is not
  converged is still a valid mesh (checked on every return), just not fully conforming.

## Decisions I made that you may want to change

* Crate name `femmesh`, no dependencies, BSD-2-Clause (the licence of the other branch; **nyx's licence is unknown to me**).
  It builds in its own workspace (`femmesh/Cargo.toml` has an empty `[workspace]`) so it does not touch PR #1's workspace.
* Adaptation parameters (all in `AdaptOptions`): split above metric length sqrt 2, collapse below 1/sqrt 2, corner angle 40 degrees, smoothing relaxation 0.7, collapse quality floor 0.2.
* PML label = `label_base + region_index * label_stride + material`, region index `sum (code_a + 1) 3^a`; defaults 1000 / 1000.
  nyx's real labels (`chips.rs`) will need an adapter.
* Edge/corner region material rule: material of the lowest-index boundary facet on the first stretched side touching it.
* Gradation uses a linear size-growth law with slope `g` (my parameterisation; not verified equal to Alauzet's).
* Storage is `.fmsh` only. VTK and any text writer were removed on request (the original stack used HDF5 plus a small XML
  descriptor; the point was to avoid large files). `.fmsh` v2 is one self-describing, CRC-checked file holding the mesh (mixed cell
  types, labels, boundary facets) and any number of result steps with named point or cell fields. Everything is lossless
  (floats round-trip bit for bit). Each array is stored with the smallest of: byte planes, XOR-with-previous, or dictionary
  coding, then an adaptive range coder (order-0 or order-1, or raw if neither helps); implemented in `io/entropy.rs`, no dependency.
  Sizes (`results/m0_fmsh_size.md`): 12-30x smaller than plain binary on structured and PML meshes, 2.1x on an adapted
  unstructured mesh, 1.7x with three full-entropy f64 fields; xz -9e is within 7 % either way (smaller by 7 % on the adapted mesh, larger elsewhere).
  Real HDF5 would need the C library or an unproven pure-Rust writer, so I did not use it. If the nyx side must read HDF5/XDMF, a thin
  exporter belongs there (`fmsh::decode_dataset` gives everything it needs). A lossy quantised option was not added, per the no-loss requirement.
  Known gap: connectivity of adapted meshes (5.5 B per triangle) could drop about 4x with a topological coder; not done.

## Still open from the handoff (questions 3, 5, 7–13) and new ones

* Max Nédélec order and families (3), what users draw for 3D geometry (5), thin structured films outside the PML (7),
  target sizes and time per cycle (8), chip export (9), QoIs (10), Floquet/mirror importance (11), offline Mmg/Omega_h
  comparison (12), `brillouin.rs` (13).
* Where is the nyx repo, and is `p-j-r/brillouin` related? (I was told not to worry about nyx; this only matters for the adapter.)
* The plan file and the nyx files were never visible to me; the summary copy "to the project files" could not be made
  (`/mnt/project-files` does not exist here). This document is the summary.
