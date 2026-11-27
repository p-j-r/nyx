> Verbatim copy of the brief (revision 3) that started this work, as it was given to the implementing session
> (it was attached to the first message; the file itself was never on disk, so this copy is recovered from the
> session transcript). The plan it refers to, `nyx_meshing_adaptivity_plan.md`, and the nyx source files were **not**
> available to the implementing session. Status against this brief: [`../HANDOFF.md`](../HANDOFF.md).

# Handoff brief: nyx metric-based mesh adaptation, structured PML meshes and goal-oriented adaptivity (revision 3)

Paste this whole file into the implementing session. It is self-contained. The reasoning, licence table and sources are in `nyx_meshing_adaptivity_plan.md` (same folder, `/mnt/project-files/rust_solvers/`); read it before starting. The nyx source files paul shared are in `/mnt/project-files/` (`mesh.rs`, `mesher.rs`, `mesher3d.rs`, `refine3d.rs`, `problem.rs`, `periodic.rs`, `par.rs`, `solid_pml.rs`, `solve.rs`, and others). Read them before designing anything.

Written 2026-10-29 for paul, revised after paul sent the nyx files. I have not seen the nyx repository itself, only those files, and I read `mesher3d.rs` only in part.

## 0. Decisions already made by paul (do not re-ask)

- nyx is a **2D and 3D** solver with proper **PML in both**; **Nédélec orders higher and lower**; **all domains** (frequency, time, eigen; optical, elastic, EM).
- **Nothing is linked or called that is not ported to Rust.** No Gmsh, Mmg, TetGen or C/C++ libraries in nyx, no subprocess calls. **Geometry is created inside nyx.**
- **Reuse existing Rust** (crates and nyx's own code).
- "Best in class" = **mesh quality, speed and low solver error**, together.
- 3D geometry is **smooth and arbitrary**. Paul has no preferred method (the plan recommends one, section 4).
- **All element types** (triangle, quad, tet, hex, prism) will be in the mesh. The **current PML meshing is wrong**: it must be **quads in 2D, prisms and hexahedra in 3D**.
- Time stepping: **implicit and explicit**, your choice (the plan recommends implicit as the default for adapted meshes, explicit with a metric floor as an option).
- **Curved elements are needed.** The **solvers are being updated** and will grow (orders above 2, curved Jacobians and Nédélec on prism/hex are their work, not yours).
- The plan was approved.

## 1. Task

Add to `femcore` (as a folder of modules `adapt/`, documented in the style of the existing modules; paul to confirm the location) a pure-Rust stack, generic over dimension 2 and 3:
1. a **structured PML generator** (quads, prisms, hexahedra) and mixed-mesh output;
2. a **metric-driven remesher** for triangles and tets with a **geometry oracle** (signed-distance based) for smooth and arbitrary surfaces;
3. **metric sources**: physics-based a priori (from `ProblemSpec`/`Materials`), derivative-based, Maxwell residual, goal-oriented adjoint;
4. an **adaptation loop**, field transfer for time-domain use, and a validation harness.

One milestone per pull request, in the order of section 5. Prove each against analytic problems and nyx's own benchmarks, and report honestly where the Rust version is worse. Three streams can run in parallel (A: PML and mixed meshes; B: adaptation; C: estimators and goal-oriented).

**Do not start coding before you have asked paul the open questions in section 7** (at least 1, 2 and 6). Do M0 while you wait.

nyx repository: **`<owner>/<name>` to be filled in by paul.** Verify every path and API before relying on it.

## 2. Licence and dependency rules (hard)

1. **You may read and port code from:** Omega_h (FreeBSD), PRAgMaTIc (BSD 2-Clause), NASA refine (Apache-2.0 as its repo page shows), deal.II (use its Apache-2.0 with LLVM exception option), MFEM (BSD-3), PETSc (BSD-2). Open each project's `LICENSE` file and confirm it **before** reading code; some files carry other licences (for example under deal.II `bundled/`, `examples/`, `tests/`). Keep notices in `THIRD_PARTY_NOTICES.md`, add a provenance table (file, derived from, upstream file, licence), and put an SPDX header on derived files.
2. **Do not read, copy or translate source from:** Mmg, ParMmg (LGPL-3.0), TetGen (AGPLv3 for 1.5.0 to 1.6.0), Gmsh (GPL), CGAL (GPL/commercial), libMesh (LGPL-2.1), Netgen (LGPL-2.1), rapidmesh (AGPL-3.0). Implement from **published papers** only.
3. **Nothing external is linked or called by nyx or by its build and test.** The only open point (question 12) is an **offline** comparison: export meshes and metrics to files, run Mmg and Omega_h by hand outside the workspace, compare statistics in a separate step. Until paul says yes, build only the Rust side of Tier B.
4. **Rust dependencies:** only MIT, Apache-2.0 or BSD crates, each checked. Already in nyx: `spade`, `faer`, `rayon`, `serde`. No new geometry crates are planned. Do not use `tritet`, `rapidmesh` or any C-backed crate.
5. Unchecked licences: BAMG, EPIC, Goalie, pyadjoint/dolfin-adjoint. Check before reading anything from them.
6. If you are unsure whether something is allowed, stop and ask paul.

## 3. Reference work to read first

1. `/mnt/project-files/rust_solvers/nyx_meshing_adaptivity_plan.md`, especially section 4 (what nyx already has).
2. The nyx files above. Pay attention to: `mesh.rs` (`SimplexMesh`, labels, `from_cells_split`), `mesher.rs` (2D, `spade`), `mesher3d.rs` (stuffing, signed distances, `newton`, crease snapping; most functions are private and must be exposed), `refine3d.rs` (the baseline), `emw_pml.rs` (`cell_stretches`, `stretch`) and `problem.rs` (`Pml`, `Materials`, `ProblemSpec`), `chips.rs` (the PML region labels and the benchmark), `periodic.rs`, `par.rs`.
3. `/mnt/project-files/rust_solvers/nyx_solver_port_review_and_handoff.md` for testing style (tiered tolerances, determinism, fused passes, timing noise) and `BENCHMARKS.md`.
4. Papers and docs (find and read; I did not read all in full): Parallel metric-based mesh adaptation in PETSc using ParMmg (arXiv 2201.02806); PRAgMaTIc thread-parallel paper (arXiv 1308.2480); CDT3D (doi 10.1007/s00366-024-01994-0); Omega_h SCOREC report 2016-25; Goalie docs, "Goal-oriented mesh adaptation" (mesh-adaptation.github.io/goalie). From memory and unchecked: Loseille and Alauzet continuous mesh framework parts I and II; "Very High Order Anisotropic Metric-Based Mesh Adaptation in 3D" (title seen, unread); Dobrzynski and Frey 2008; Becker and Rannacher on DWR; Beck, Hiptmair, Hoppe and Wohlmuth on residual estimators for Maxwell; Demkowicz on projection-based interpolation; Labelle and Shewchuk 2007 (isosurface stuffing, already implemented in nyx). Confirm each citation.
5. The earlier solver port is in `p-j-r/inspired-rust-port`, PR #1 branch `claude/project-thread-x1wdg8` (head `5ad411c`, not merged; `main` has only the initial commit). Reuse its patterns: fixed-chunk reductions (`src/par.rs`), run-time ISA dispatch with bit-identity tests (`src/kernels.rs`), per-phase timers, the `bench/run.py` report style. Do not push to that branch.

## 4. Architecture in one page

```
adapt/metric.rs       SPD algebra 2x2/3x3: length, quality, interpolation, intersection, gradation,
                      derivative-to-metric, Lp complexity normalisation, element-count targeting,
                      isotropic metric from the existing size closures
adapt/geom.rs         geometry oracle: sdf, classify (surface/crease/corner), project; from Spec2d/Spec3d
adapt/topo.rs         working mesh, classification from cell and facet labels, rebuild via SimplexMesh::new
adapt/ops.rs          split, collapse, 2-3/3-2 flips, edge removal, smoothing; quality-gated, new-before-old
adapt/drive.rs        adapt(mesh, metric, geometry, options) -> SimplexMesh; serial first, then independent sets
adapt/pml.rs          structured PML around a box-shaped physical mesh -> CellMesh (quad / prism / hex)
adapt/src_phys.rs     a priori metrics from ProblemSpec/Materials (optical, elastic, PML, thin layers)
adapt/src_deriv.rs    derivative-based metrics;  adapt/est_maxwell.rs  residual estimator
adapt/goal.rs         DWR: adjoint, enrichment, eigenvalue QoI
adapt/transfer.rs     field transfer (time domain, warm starts)
adapt/loop.rs         solve-estimate-adapt loop, stopping rules, history log
```
Design constraints:
- **Match the surrounding code**: `Result<_, String>`, `//!` module docs, `progress::check()` and `progress::stage()` in every long loop, `NYX_*` environment overrides, inline tests, `par::threads()` pool (do not start a second pool).
- Generic over D = 2, 3. Deterministic for a given thread count; thread-count independent for combinatorial operations (independent sets), fixed order for floating-point sums.
- Fusion hooks and per-phase timers from the first commit.
- Metric sources combine by **intersection**, then normalise and grade. Frequency sweeps, modes, time windows and the optical-plus-elastic Brillouin case use the same mechanism.
- **Frequency and eigen problems are solved again on the new mesh; no field transfer is needed.** Only the metric is interpolated. Transfer is for time domain and warm starts.
- **PML is regenerated after adaptation** by sweeping the box-shaped interface. Adapt the physical domain only, with the box faces, edges and corners classified so that vertices stay in their plane, line or corner.
- **PML details (plan 5.4):** each region is its own axis-aligned box with its own domain label and `Pml` directions (faces: one direction, edges: two, corners: three), because `cell_stretches` takes the bounding box of a domain's cells; a structure crossing the layer carries its material into a PML-labelled copy of that material (see `chips.rs`); 2D uses quads, 3D uses prisms on faces and hexahedra on edges and corners, so triangular faces meet the physical tets and quad faces meet each other (no pyramids; verify on the corners); output a `CellMesh`; with today's solvers it goes through `SimplexMesh::from_cells` (the same mesh solved natively and split must agree).
- **Geometry oracle (plan 5.3):** keep the signed-distance representation of `mesher3d`; expose its distance, projection and crease code; surface vertices slide by projection. Extending the primitives and an SDF from triangulated surfaces are later items, with a spike first.
- **Interpolation-error order:** the metric uses the (p+1)-th derivative for a full degree-p space; for lowest-order first-family Nédélec, E converges O(h), so use the first derivative (hypothesis; verify).
- **Adjoint (plan 5.7):** nyx's PML matrices are complex symmetric, so `A^H y = c` is the forward solve with conjugated data, and the left eigenvector of `A x = lambda B x` is `x` under the unconjugated product. Verify both on small cases before relying on them.
- **Constraints:** Floquet-periodic pairs must stay identical up to the translation (`periodic.rs` errors otherwise); mirror symmetry (`mirror_axis`) needs a symmetric mesh. First version: freeze those boundaries or planes and mirror the result.
- **Curved elements (plan 5.9):** needs stored high-order nodes in `SimplexMesh` (today midpoints are computed) and non-affine Jacobians in the solvers. Do the mesh side (projection plus Jacobian checks) only after agreeing the interface with the solver team.

## 5. Milestones (one PR each)

| # | Stream | Milestone | Done when |
|---|---|---|---|
| M0 | all | Analytic test problems (2D, 3D), the four chips from `chips.rs` as benchmarks, mesh statistics and VTK output, timing and determinism patterns | Fixtures and report format exist; the baseline numbers of the hand-tuned chip meshes at several `scale` values are recorded; if paul agrees, an export format for the offline comparison exists |
| M1 | B | `adapt/metric.rs` | Tier-A algebra tests; fused derivative-to-metric kernel; benchmark vs a naive version |
| M2 | A | `adapt/pml.rs` in 2D | Structured quad PML around a box-shaped 2D mesh, labels per region, split fallback; chip eigenvalues compared with the current triangular PML; layers needed for a given reflection level reported; **stop and report** |
| M3 | B | `adapt/topo.rs` + `adapt/geom.rs` (2D) | Classification from labels; validity checks; deterministic |
| M4 | B | 2D adaptation operators and driver | Adapts to analytic metrics; quality statistics; beats the hand-tuned meshes on time to target error or the report says it does not; **stop and report** |
| M5 | C | A priori physics metrics (2D) | Meshes for the four chips from `ProblemSpec` alone; error vs degrees of freedom |
| M6 | B+C | Derivative-based metrics, normalisation, gradation, loop | End to end on the chips, orders 1 and 2, frequency domain |
| M7 | A | `adapt/pml.rs` in 3D | Prisms and hexahedra; conforming; mixed-cell results format (HDF5 `cells`, XDMF) extended; `emw3d` cavities with PML agree between native and split paths |
| M8 | B | 3D oracle from `Spec3d` (expose the private functions) and 3D operators | Tier A and quality statistics in 3D |
| M9 | B+C | 3D end to end; graded metric for `mesher3d`'s lattice | `emw3d` examples; **stop and report** |
| M10 | C | Maxwell residual estimator | Effectivity reported for orders 1 and 2 |
| M11 | C | DWR with eigenvalue QoI | Effectivity index; error per degree of freedom vs derivative-only adaptation |
| M12 | B | Periodic and mirror-constrained adaptation | Floquet and mirror examples solve |
| M13 | B | Threaded adaptation by independent sets | Same bits per thread count; scaling reported; **stop and report** |
| M14 | D | Curved geometry (needs solver support) | Convergence on a sphere or ring at order 2 |
| M15 | D | Field transfer and time-domain mesh sequences | Pulse-through-layers test once a time-domain solver exists |
| M16 | B | More primitives, SDF from triangulated surfaces | Spike first |

Effort is guessed in the plan (about 35 to 50 sessions); re-estimate after M4 and M8.

## 6. How to test

Acceptance tiers:
- **A. Exact properties:** positive volumes (areas in 2D) after every operation; every interior face (edge in 2D) shared by exactly two cells; classification preserved; total volume preserved to roundoff; conforming meshes including prism/hex/tet faces at the PML interface; SPD metric algebra; **native mixed PML mesh and its simplex split give the same eigenvalues**; the same bits for the same thread count; gradient fields stay curl-free through transfer.
- **B. Against reference codes (offline, only if paul agrees):** same input mesh and metric exported and run through Mmg and Omega_h by hand; compare element count, % edges with metric length in [1/√2, √2], mean and minimum mean-ratio quality, histogram; thresholds set after measuring. If declined: analytic quality targets and published statistics.
- **C. End-to-end accuracy:** error vs degrees of freedom and vs wall-clock against analytic or very-fine references, versus uniform refinement, `refine3d` and the hand-tuned nyx meshes.

Benchmarks first: the **four chips** (`neff`, elastic eigenfrequencies, Brillouin gain) with hand-tuned, adapted-triangular-PML and adapted-structured-PML meshes; **3D `emw3d` cavities** with and without PML; **PML reflection level vs layers** for structured and triangulated layers. Analytic problems in 2D and 3D: boundary and interior layers, anisotropic Gaussian, PEC cavity eigenvalues, scattering with PML (Mie series; cylinder series in 2D), thin slab or waveguide, L-shaped or Fichera corner, order-by-order convergence, P and S wave resolution, goal-oriented effectivity index near 1.

Rules (from the solver port):
1. Record machine, compiler flags and (if used) reference-code build configuration with every table.
2. Timing: minimum of repeated runs after warm-up, interleaved, in one session. The VM varies about ±15%, so under 10% is a tie. Never compare a container number to paul's desktop.
3. Report **time to a target error**, not only remesh time; measure the split first (the complex LU of the PML systems will dominate).
4. Check generated assembly for every hot kernel; no `target-cpu=native`. Measure thread wake-up cost on paul's Windows machine before choosing a parallel scheme; nyx already has a shared pool and a smaller factorisation pool, so never add a third.
5. If Rust is worse, diagnose before reporting (assembly, `callgrind`, per-phase timers) and report the cause.
6. Recompute quality and error independently of the code under test.

## 7. Ask paul before coding

Answered: section 0. Still open (do not guess):
1. What will the updated solvers consume (mixed-cell mesh type, element bases on prism/hex, order cap, timing)? **(blocking for M2 and M7)**
2. Who adds non-affine Jacobians and stored high-order nodes; may `SimplexMesh` change?
3. Maximum Nédélec order and families needed.
4. PML parameters (thickness, layers, `strength`, `order`, grading); must the physical domain always be a box, or are non-box PML shapes needed? **(blocking for M2)**
5. What do users draw for "weird" 3D geometry (CSG of the four primitives, imported surfaces, parametric curves)? Is an STL-style import wanted?
6. Modules inside `femcore` or a separate crate? nyx repo owner/name and licence. **(blocking)**
7. Do users ask for structured hex/prism regions outside the PML (thin films)?
8. Target sizes and wall-clock per adaptation cycle; thin-layer aspect ratios.
9. May real chip geometries be exported for tests (the four in `chips.rs` are usable now)?
10. Quantities of interest for goal-oriented work.
11. How important is adapting Floquet-periodic and mirror-symmetric models?
12. Is an offline comparison against Mmg and Omega_h on exported files acceptable (never linked or called by nyx)?
13. `brillouin.rs` did not reach me; ask for it if meshing touches it.

## 8. Commit and PR conventions (as used in PR #1)

- One logical milestone per branch and per PR. Commit message: imperative summary line, a short paragraph on what changed and why, then the attribution trailer lines **your own session instructions** specify. Do not copy mine.
- PR body: the two-line attribution block your session context requires first, copied exactly, then **Summary**, **Result** (statistics table with the tier outcome), **Testing** (exact commands, ISA tiers and thread counts), **Licence** (provenance of anything derived), and the footer your instructions require. Mirror a PR template if the repo has one.
- After creating the PR: assign `p-j-r`, request their review, and subscribe to PR activity so you fix CI and review comments yourself.
- Never force-push over someone else's work, never skip or quarantine a test to get green, never merge.
- Put results under a results folder in the repo and copy a human-readable summary to the project files so paul can read it.

## 9. Working agreement

- Stop and report to paul after M2, M4, M9 and M13, and whenever a licence question appears.
- Say what you verified and what you did not. Tag unmeasured claims as hypotheses.
- If a milestone's result is negative (for example derivative-based metrics do not help for edge elements, or structured PML is not cheaper), report it rather than adding code to make the number look better.
- For any new crate or library: check its licence first, record it, and confirm it is pure Rust.

## 10. Known gaps in what I could tell you

- I have not seen the nyx repository, `Cargo.toml` or `solver_port_list.md`. I read the source files paul attached, not all of them (the plan lists which), and `mesher3d.rs` only in part.
- `brillouin.rs`, `material.rs`, `progress.rs`, `scalar.rs`, `solid*.rs` (except `solid_pml.rs`), `sparse.rs`, `validate.rs`, `view.rs`, `xmlio.rs`, `lobpcg.rs` and `ring_brillouin.rs` were not read.
- Nothing about the updated solvers is known to me beyond paul's message.
- Algorithm details of Omega_h's split/collapse/swap, Mmg's internals, the Nédélec parts of MFEM and deal.II, and the high-order metric and projection-based interpolation theory were not read.
- The claims in section 4 that the adjoint is a conjugated forward solve and that the left eigenvector equals the right one for complex-symmetric pencils are derived by me and unverified numerically.
- Unchecked licences: BAMG, EPIC, Goalie, pyadjoint. Unchecked features: libMesh adjoint estimators, MFEM estimators and Nédélec on non-conforming meshes, deal.II step-14.
- Effort numbers are guesses.
