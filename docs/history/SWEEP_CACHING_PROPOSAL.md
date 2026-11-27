# Parameter sweeps: what nyx and femcore should change to use geometry-side caching

Proposal, 2026-10-31. No code written, nothing measured. Companion to `GEOMETRY_ARCHITECTURE.md` and `GEOM_HANDOFF.md`; `geom` stays a separate crate that
femcore calls (`app -> femcore -> geom`).

## 1. Where the time actually goes (from the repo's own logs and docs)

Reuse only pays where the cost is. From `HOW_IT_WORKS.md` section 11 and the README:

| Case | Reported cost | Dominant part |
|---|---|---|
| SOI wire, 2 modes (2D) | ~0.1 s | everything is small; a 50-point sweep is seconds already |
| Solid cube, 6 modes, LOBPCG | 7.4 s total | assembly 0.19 s + 0.03 s; multigrid setup 0.27 s; **eigensolver 6.4 s** |
| Brillouin chip (example 15) | ~43 s | **complex factorisations 4-10 s each** (one per optical solve, one for the elastic solve); elastic eigenmodes 24.9 s of one pump/Stokes step |
| Bragg cell, coarse | ~40 s each | 3D complex eigenproblem |

So assembly is not the bottleneck. The expensive things are **factorisations, eigensolves and (for adaptive/remeshed runs) meshing**. The proposal targets
those: symbolic reuse, warm starts, preconditioner reuse, avoided remeshing. Sweeps of 2D cross-sections gain little; 3D and Brillouin sweeps gain most.
Step 0 below measures before building.

## 2. What already exists in femcore (and is under-used)

- `solvers::direct::ShiftedFactorization`: symbolic analysis once, numeric factorisation per shift (`replace_matrices`). Used by the frequency sweeps.
- `solvers::sweep::rom_sweep`: reduced-order frequency sweeps with a posteriori residual certificate. Its own doc says it has **no multi-parameter models**.
- `solvers::nep` split form `T(w) = sum_j f_j(w) A_j`: the affine structure that makes parameter changes cheap.
- `track_modes` (Brillouin): follows modes between two `q` by overlap of shapes.
- Shared cancellable progress (`Progress`), a worker thread in the app, `meshers::quality` statistics, edge/dof numbering cached inside `SimplexMesh`.
- Missing: `SolveState` carries only `heat_t`; `solve::solve(mesh, spec, state)` rebuilds everything per call; nothing links geometry changes to the solver.

## 3. Three kinds of sweep parameter, three reuse levels

| Class | Examples | What changes | Reuse |
|---|---|---|---|
| **A: no mesh change** | frequency/wavelength, material constants, BC values, loads, kz, Floquet phase | matrix entries only (pattern identical) | everything geometric; symbolic factorisation; preconditioner structure; warm starts; **parametric ROM** (affine in the parameter) |
| **B: same topology, vertices move** | width, thickness, radius, sidewall angle, period (z only in a layer stack) | node coordinates; connectivity, labels, dof numbering, sparse pattern unchanged | all of A, plus **mesh morphing** instead of remeshing |
| **C: topology or size changes** | a feature appears/disappears, large moves, new holes | everything | only geom-side caches (edge discretisations, 2D overlay mesh when only z changes), fresh symbolic analysis |

Class B is the new opportunity, and it also has an accuracy side effect to test (hypothesis): remeshing from scratch at every point adds mesh noise to a
curve like neff(w), while a morphed mesh varies smoothly with the parameter. That matters for derivatives and optimisation.

## 4. Proposed femcore changes

### 4.1 Solve plan and cache (replaces the hidden rebuild-everything)

```rust
SolveCache {                                  // generalises SolveState (which only holds heat_t today)
    layout_key:  u64,                         // hash of (mesh topology signature, element order, dof layout, constrained-dof set)
    symbolic:    Option<SymbolicFactorization>,   // from ShiftedFactorization, valid while layout_key matches
    amg:         Option<AmgHierarchy>,        // reused as preconditioner; rebuilt when PCG iterations exceed a threshold
    warm:        WarmStart { vectors: Vec<Vec<C64|f64>>, shifts: Vec<f64> /*beta, omega, n_target*/, shapes: ... },
}
solve::solve_cached(mesh, spec, &mut SolveCache, progress)   // same results as solve_with within solver tolerance
```

Rules: the cache key is structural (topology + layout), never geometric, so class A and B hit it; any miss rebuilds that item only and logs it
("symbolic analysis reused", "AMG rebuilt: 31 > 25 iterations"). Warm starts feed LOBPCG its initial block, shift-and-invert Krylov its starting
block and shift (previous `beta`/frequency), and PCG its initial guess; `track_modes` aligns mode order across points so sweeps follow a branch instead of re-sorting.
Fallbacks already in femcore (multigrid -> factorisation, LOBPCG stall -> factorisation) stay and are logged per point.

### 4.2 Affine term assembly for class A

Assemble per-domain (and per-boundary-condition) term matrices once, `A(p) = sum_j f_j(p) A_j`, with all terms on the union pattern (the same trick as
`ShiftedFactorization`). A material or frequency change is then a sparse axpy. Because assembly is cheap next to factorisation (section 1), this is **second
priority**; its main value is feeding the parametric ROM.

### 4.3 Multi-parameter ROM (extends `rom_sweep`)

Generalise the reduced basis from one frequency to several affine parameters (frequency, a few material constants, kz): greedy selection of expansion
points over a parameter grid with the existing a posteriori residual check. Applies to class A only (class B changes the matrices non-affinely). Needs a
review of `docs/ports/solvers/STATUS.md` first, since the doc explicitly lists parametric models as not done.

### 4.4 Mesh morphing for class B (`meshers::morph`)

Inputs: old mesh, the node association (which exact face/edge each boundary node lies on, **with its parametric position**), the old and new geom evaluations.
Steps: (1) new boundary node positions from the new exact entities at the same parametric positions; (2) extend the displacement to interior nodes (harmonic/
linear-elasticity extension solved with our own solvers, or inverse-distance weights as a cheaper first version); (3) check for inverted cells and quality
loss with `meshers::quality`; (4) return the same `SimplexMesh` (same connectivity, labels, edges, dof numbering) with new coordinates, or `NeedsRemesh`.
A periodic anchor (full remesh every N points or when quality drops below a threshold) keeps drift bounded. Morphing only runs when geom reports an unchanged
topology signature. For P2, edge midpoints are recomputed (on exact curves where the association says so).

### 4.5 Staged mesh cache for class C and for layer stacks

Mesh drivers cache by stage with content hashes of that stage's inputs: (1) per-edge discretisation (geom already caches), (2) the 2D overlay mesh,
(3) the extrusion. A sweep that only changes layer thickness or a Bragg period changes stage 3 only: the 2D mesh is reused and the extrusion is rerun,
which is nearly free. A sweep over a width changes stage 2 (class B morph or class C remesh).

### 4.6 `SweepSession` (orchestrator)

```rust
SweepSession::new(model: geom::FeatureList, problem: ProblemSpec) 
session.run_point(&params) -> SweepPoint { scalars, optional fields, report: {reused: [...], remeshed: bool, iterations, times} }
session.run(plan: SweepPlan, progress, cancel) -> SweepResult
```

Per point it asks geom for the evaluation, compares topology signatures with the previous point, picks class A/B/C, runs the cheapest valid path, and records
what it reused. `SweepPlan` orders points for continuity (sorted / nearest neighbour), so warm starts are always from a close point. **Parallelism**: points
in different chains can run concurrently, but femcore has one shared thread pool, a separate pool of 8 for factorisations, and 3D memory is large (about 13 GB for
3M unknowns), so the plan caps concurrency by an estimated per-point memory and defaults to sequential-with-parallel-solver. Measure before enabling.

### 4.7 Correctness rules

- Cached and cold results must agree within the solver tolerance: a test per physics runs a short sweep both ways (`NYX_NO_CACHE=1` forces cold).
- Results are not bitwise independent of sweep order (warm starts change iterations), only within tolerance. Re-running a single point may differ at the
  solver-tolerance level from its value inside a sweep. State this in the docs; the existing "same input gives same numbers" holds for the same sweep plan.
- Morphing never silently degrades: quality thresholds and the periodic anchor are logged, and a morphed result at a checkpoint is compared with a fresh
  remesh in tests (difference must be within the mesh-convergence band already documented per use case).

## 5. Proposed nyx (app) changes

- A named **parameter table** in the project (shared by geometry expressions and by material, BC and study values), saved in `.nyx.xml`.
- A **Sweep** section in the Study tab: pick parameters, ranges (linear, log, list), outputs to record (neff, f, Q, gain, any scalar the study returns), and which
  points keep full fields (memory policy). Runs in the existing worker thread with progress and Cancel.
- Results: a table and `egui_plot` curves per scalar, parameter columns saved with the run (autosave folder; HDF5 builds also store them); partial sweeps
  are resumable from the saved points.
- Per-point log line showing what was reused or rebuilt ("symbolic reused, warm start, morphed"), so speed-ups are visible rather than assumed.
- Batch: `nyx --batch project.nyx.xml --sweep w=400e-9:600e-9:21` prints the table.
- Legacy projects: class A sweeps work immediately; geometry sweeps need the new `geom` path.

## 6. What `geom` needs to add (small; none of it is in the handoff yet)

1. `topology_signature() -> u64`: hash of the combinatorial structure (regions, faces, edges, loops, incidence, surface kinds), independent of parameter values.
2. **Stable ids within a signature**: same signature implies identical id sets and ordering (the handoff asks for determinism and stability under *unrelated*
   edits; this extends it to *parameter* changes).
3. **Parametric positions in the association**: boundary nodes carry `(entity id, parametric fraction)` (one float on an edge, two on a face), and geom offers
   `position_at(entity, fractions)` on the new evaluation. This is the input to morphing (4.4), and it is a small change to the association contract in
   section 10.2 of the architecture doc.
4. `diff(old_eval, new_eval) -> ChangeReport { changed_features, signature_equal, kind: None | ZOnly | Morph | Remesh }` so the session and the staged mesh cache can
   choose a path without recomputing hashes.
5. Edge parameterisation continuous in the parameters (no jumps of seam or start point when a radius changes), documented.

## 7. Order of work

| Step | What | Depends on geom? | Done when |
|---|---|---|---|
| 0 | Measure: run `NYX_PROFILE=1` on a SOI-width sweep, a Bragg cell and a chip step; run `examples/symbolic_share` (how much of a factorisation is symbolic); write a `sweep_bench` baseline (cold, every point) | No | Table of where time goes per physics and a baseline to beat |
| 1 | `SolveCache`: symbolic reuse, warm starts, AMG reuse; class A sweeps; cold-vs-cached tests | No | Class A sweep faster by a measured factor, results within tolerance |
| 2 | `SweepSession`, parameter table, Sweep UI and batch flag, result table/plots | Only for geometry parameters (legacy `Spec2d`/`Spec3d` parameters can be driven by hand meanwhile) | A frequency/material sweep runs end to end in the app |
| 3 | Staged mesh cache (layer stacks: z-only changes reuse the 2D mesh) | Yes (stage hashes, `diff`) | Bragg period sweep skips 2D meshing |
| 4 | `meshers::morph` for class B | Yes (items 1-5 above) | A width sweep stays on one topology, quality log clean, checkpoint remesh agrees |
| 5 | Multi-parameter ROM | No (class A) | Residual-certified multi-parameter sweep matches full solves |
| 6 | Cross-point parallelism with memory cap | No | Throughput gain measured on a 2D sweep without exceeding memory |

Steps 0-2 do not need `geom` at all and give the class A benefit right away; steps 3-4 are where geom's caching and stable ids pay off.

## 8. Risks and open questions

- **Morph quality**: large parameter ranges will need the anchor remesh often; the benefit is mostly for small to moderate steps. Measure the fraction of
  points that morph successfully on the real use cases before investing in the interior-extension solver.
- **Preconditioner reuse** may need a rebuild frequently for eigenproblems whose matrix changes a lot; the rule is "rebuild when iterations grow", not a fixed schedule.
- **Memory** for retaining warm-start vectors and per-point fields in 3D sweeps.
- I have not read `solvers/sens.rs` or `solvers/opt.rs`; if they cover parameter sensitivities, smooth morphed meshes would let them be used across geometry
  parameters. Worth a look before step 4.
- Decision for you: is the first priority 3D/Brillouin sweeps (where the time is) or the app-level sweep workflow (Study tab, plots), which is useful even
  for fast 2D cases? I would start with step 0 and step 1 either way.
