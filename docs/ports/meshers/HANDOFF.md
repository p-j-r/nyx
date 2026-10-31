# Handoff: `femmesh` (meshers/)

For the session that takes this over with proper access to the nyx repository. Written at the end of the implementing session.
Read in this order: this file, [`PROVENANCE.md`](PROVENANCE.md), [`femmesh/docs/STATUS.md`](femmesh/docs/STATUS.md), the original brief
[`docs/nyx_meshing_adaptivity_handoff.md`](docs/nyx_meshing_adaptivity_handoff.md) (sections 2, 5, 6, 7, 9 matter most).

## 1. What this is and how it was built

A standalone, dependency-free Rust crate, `meshers/femmesh`, that implements part of the brief's meshing stack. The author was told "you do
not need to worry about nyx, just build it as `femmesh`, another session moves it inside `femcore`". So:

* **The nyx repository, `nyx_meshing_adaptivity_plan.md` and the nyx source files (`mesh.rs`, `mesher.rs`, `chips.rs`, ...) were never visible.**
  `src/mesh.rs` (`SimplexMesh<D>`) and `src/cellmesh.rs` (`CellMesh<D>`) are **stand-ins written from the brief's description**. The adaptation code
  only relies on the accessors in those two files, so porting means swapping them for the real types.
* The PML label scheme is my own (`label_base + region*stride + material`, defaults 1000/1000) and is **not** the one in `chips.rs`; an adapter is needed.
* The original user answers: placement "inside femcore eventually, another session handles that"; mixed meshes as `CellMesh` plus a split to simplices;
  PML boxes only; no VTK (the old stack used HDF5 plus a small XML descriptor and the user does not want big files).
* Nothing has been run against the real solvers, the four chips, or hand-tuned meshes. Every number is from a stand-in scalar Helmholtz solver or from interpolation-error tests.

Machine for all numbers: 4-core Xeon 2.8 GHz VM, rustc 1.97, release profile (fat LTO), no `target-cpu=native`; timing noise about +-15 %.

## 2. Build, test, run

```bash
cd meshers/femmesh
cargo test --release                       # 86 tests, about 15 s (the 3D solver tests are ~11 s of that)
cargo clippy --release --all-targets       # clean
cargo fmt --check                          # rustfmt.toml: max_width 120
cargo run --release --example pml_reflection   # M2 study (2D)
cargo run --release --example pml_3d           # M7 study
cargo run --release --example adapt_convergence -- 32000   # M4 study; argument caps the largest target (full run is minutes)
cargo run --release --example bench_metric     # M1 kernel benchmark
cargo run --release --example fmsh_size        # storage study (uses gzip/xz/zstd if installed, only for comparison)
```

Gotchas: the 3D banded solver is O(n * bw^2) in time and O(n * bw) in memory; keep 3D harness problems under about 12 000 unknowns. Long runs
(`adapt_convergence` without an argument, `pml_3d`) should be started detached and polled; foreground tool calls have a time limit. `FEMMESH_THREADS`
sets the thread count of the parallel reductions; `FEMMESH_VERBOSE=1` prints adaptation iterations.

## 3. Status against the brief's milestones

| # | State | Where | What is still missing |
|---|---|---|---|
| M0 | partly | `fixtures.rs`, `quality.rs`, `par.rs`, `timer.rs`, `io/fmsh.rs` | The four chip benchmarks and baseline numbers of the hand-tuned chip meshes (need `chips.rs`); offline Mmg/Omega_h export (question 12 unanswered). |
| M1 | done | `adapt/metric.rs` | 3x3 eigen-solve is 66 % of the 3D kernel time (Jacobi); a closed form is untried. |
| M2 | done except chip comparison | `adapt/pml.rs`, `validate/` | Chip eigenvalues vs the current triangular PML (need nyx). Reported at the M2 stop; no reply was received. |
| M3 | done (2D) | `adapt/topo.rs`, `adapt/geom.rs` | |
| M4 | done (2D) | `adapt/ops.rs`, `drive.rs`, `field.rs` | Comparison with hand-tuned nyx meshes and with Mmg/Omega_h. Reported at the M4 stop. |
| M5, M6 | **not started** | | A priori metrics need `ProblemSpec`; the loop needs a solver and an estimator, i.e. nyx. |
| M7 | done as far as possible without nyx | `adapt/pml.rs` (3D), `validate/helmholtz.rs` | Vector `emw3d` cavity agreement between native and split; the HDF5/XDMF results format is replaced by `.fmsh` (decision below). |
| M8-M16 | not started | | |

The brief says: stop and report after M2, M4, M9 and M13. M2 and M4 were reported. The next stop is M9.

## 4. Module map (`femmesh/src`)

| File | Role |
|---|---|
| `mesh.rs` | Stand-in `SimplexMesh<D>`: flat cells, labels, boundary facets; `check()`, edges, facet incidence, orientation. |
| `cellmesh.rs` | Stand-in `CellMesh<D>`: Tri, Quad, Tet, Prism, Hex; `split_to_simplices` (conforming), `check()` (including inverted-cell detection by scaled Jacobian). |
| `adapt/metric.rs` | SPD metric algebra for D = 2, 3 (`Sym<D>`, intersection, log-Euclidean mean, gradation, edge length), Lp Hessian-to-metric kernel with complexity normalisation. **Dimension-generic; reusable for M8 as is.** |
| `adapt/field.rs` | `MetricField<D>` trait; `UniformSize`; `BackgroundMetric2` (metric stored on a background mesh, bucket grid, interpolated). 2D only. |
| `adapt/topo.rs` | 2D working mesh `Work2`, vertex classes (interior / chain / corner), constrained edges, chains, frozen chains, CSR snapshot. 2D only. |
| `adapt/geom.rs` | 2D geometry oracle: polylines from the mesh, analytic projectors (circle, signed distance). |
| `adapt/ops.rs`, `adapt/drive.rs` | 2D split / collapse / flip / smooth and the adaptation loop (`adapt2`, `adapt2_with`, `AdaptOptions`, `AdaptReport`). |
| `adapt/pml.rs` | `PmlSpec<D>`, `build_pml`: structured PML around a box-shaped simplex mesh (quads in 2D; prisms on faces, hexes on edges and corners in 3D). |
| `validate/` | Independent scalar Helmholtz FE solver with PML (2D and 3D, mixed elements), banded LU, RCM, inverse iteration. A test harness, not a solver for production. |
| `io/fmsh.rs`, `io/entropy.rs` | `.fmsh` v2 file: mesh plus result steps, lossless, CRC-32, entropy coded. |
| `fixtures.rs`, `quality.rs`, `par.rs`, `timer.rs` | Test meshes and analytic fields with exact Hessians; quality and edge statistics; deterministic chunked reductions; phase timers. |
| `tests/adapt2d.rs` | 11 integration tests of the 2D adaptation with independently recomputed figures. |

## 5. Results worth knowing (details in `femmesh/results/` and `femmesh/docs/STATUS.md`)

* **M1:** fused Hessian-to-metric kernel is 1.46x faster in 3D (tie in 2D) and hits the target complexity to roundoff; bit-identical for any thread count.
* **M2:** for scalar P1/Q1, structured quad PML layers are **not more accurate** than the same layers split into triangles (same vertices, same dof);
  about 6 layers give 0.5 % and 8 layers 0.1-0.2 %. The case for quads rests on half the cell count and control of layers, not accuracy. The brief's Tier-A
  line "native mixed PML mesh and its split give the same eigenvalues" cannot hold exactly (different spaces); it should be stated as agreement to discretisation error.
* **M7:** 3D scalar harness; native hex/prism/tet and split-to-tet eigenvalues of the layered mesh both converge at second order and agree to 1.2 % at the finest level;
  native is 4.3-4.5x more accurate on the same vertices (cause not isolated). 3D PML works for both. 1-2 layers are unresolved (16-38 % error). The two references differ by 2.4 %.
* **M4:** adapted P1 interpolation error falls as 1/N; 85-95x lower than a uniform mesh with the same triangle count for an anisotropic Gaussian, 1400-6400x for thin
  layers/ridges (mostly the weakness of the uniform baseline). q_min 0.28-0.60, about 99 % of edges in band. About 15 k triangles/s serial. **No comparison with Mmg or Omega_h was made.**
* **`.fmsh`:** 12-30x smaller than plain binary on structured and PML meshes, 2.1x on an adapted unstructured mesh, 1.7x with full-entropy float fields. xz -9e is within 7 %.

## 6. Known limitations (do not rediscover these)

1. **Adapting to an ungraded metric with sharp features is unreliable** (6 of 30 ungraded study runs did not meet the stopping rule; one ran away until `max_cells`). All 15 graded runs converged.
   Use `BackgroundMetric2::new_graded` or `metric::grade` (slope 0.4 in the study). Cost: 1.1-3x more elements than the target.
2. A feature thinner than the starting-mesh spacing can be missed. Midpoint sampling in the first two iterations helps; gradation hides it. Matters most for the first cycle of a solve-estimate-adapt loop.
3. Frozen boundaries must already be resolved by the input mesh (an over-long frozen edge cannot be split; neighbours are left alone and `unresolved_frozen` is reported).
4. Everything in `topo`, `geom`, `ops`, `drive`, `field` is **2D and serial**. M8 needs 3D equivalents; only `metric.rs`, `quality.rs`, `mesh.rs`, `cellmesh.rs` generalise.
5. Adapted-mesh connectivity compresses poorly (5.5 B per triangle) because vertex ids have no spatial locality. A topological coder would cut it about 4x; not done.
6. The first adaptation iteration overshoots (many splits, then collapses); candidate lists are rebuilt each pass. Room for speed before M13.
7. `.fmsh` is new (v2); there are no other files in the wild, but there is no reader for any other format. No HDF5/XDMF export exists.
8. The banded LU does not pivot. It is adequate for the PML systems here; check the returned residual.

## 7. Decisions you may want to change

* Crate name `femmesh`, BSD-2-Clause, own workspace root (`[workspace]` in its `Cargo.toml`), no dependencies. Keep it dependency-free unless a pure-Rust MIT/Apache/BSD crate is justified (brief, section 2).
* Storage: `.fmsh` only, no VTK/text writer, no HDF5 (needs the C library or an unproven pure-Rust writer). If nyx needs HDF5/XDMF, write a thin exporter on the nyx side from `fmsh::decode_dataset`.
  A lossy option was deliberately not added ("no loss of information"). If the user later accepts quantised *fields*, that is where the large saving is.
* Adaptation defaults (`AdaptOptions`): split above metric length sqrt 2, collapse below 1/sqrt 2, corner angle 40 degrees, smoothing relaxation 0.7, collapse quality floor 0.2, stop when splits plus collapses are under 0.2 % of the triangles.
* Edge/corner PML region material rule: material of the lowest-index boundary facet on the first stretched side touching it.
* Gradation: linear size-growth law with slope `g` (my parameterisation).

## 8. What to do next (suggested order)

1. **Get access, then answer the blocking questions** (brief section 7): 1 (what the solvers consume), 4 (PML parameters, box only was assumed), 6 (femcore module vs crate, nyx licence).
   Unanswered and still open: 2, 3, 5, 7-13.
2. **Replace the stand-ins.** Swap `mesh.rs` and `cellmesh.rs` for nyx's types, write the PML label adapter against `chips.rs`, and re-run `cargo test` to find every place the stand-in API differed.
   Check the node orderings in the `cellmesh.rs` header against nyx's before trusting the prism/hex Jacobians.
3. **Finish what was blocked on nyx:** M0 chip benchmarks and hand-tuned baselines; M2 chip-eigenvalue comparison against the current triangular PML; M7 `emw3d` native-vs-split cavity check
   (state the tier-A line as "agree to discretisation error"); M4 comparison with the hand-tuned meshes (time to target error).
4. **M5/M6** with a real `ProblemSpec` and solver: a priori metrics, derivative-based metrics, the solve-estimate-adapt loop. Use graded metrics by default (limitation 1).
5. **M8/M9:** 3D oracle, 3D topology/operators (start from the 2D structure; `metric.rs` is already dimension-generic), then end to end on `emw3d`. **Stop and report at M9** (working agreement).
6. Quality-of-life that would pay off early: closed-form 3x3 symmetric eigen-solve (M1 hypothesis), incremental candidate lists in the 2D passes, a topological connectivity coder for `.fmsh` if file size matters for adapted meshes.
7. Conventions from the brief (section 8): one milestone per branch and PR, PR assigned to `p-j-r`, review requested, never force-push over others, never skip a test, never merge. **No PR was opened by the previous session.**

## 9. Process lessons that cost time

* Compare each discretisation with a reference of its own kind; a single reference favours whichever kind it was built with (M2 and M7).
* Always recompute errors and quality independently of the code under test (the tests do).
* Do not `pkill -f <name>` from the shell that is running the job (killed the session shell once); run long jobs with `setsid nohup ... &` and poll a marker file.
* Interleave flip passes between split passes; 30 split passes in a row gave q_min 0.000.
* A guard that skips splits near the minimum length (to stop split/collapse ping-pong) made things worse; the `converge_frac` stopping rule worked.
* Compression: the first `.fmsh` coordinate coding (XOR) did not help decimal coordinates; per-axis dictionary and byte planes did. Measure on adapted meshes, not only structured ones.

## 10. What was not verified

* Anything against nyx (APIs, labels, chips, `emw3d`, Nédélec elements). Anything about vector Maxwell on hex/prism/tet interfaces.
* The recalled literature citations (see `femmesh/docs/REFERENCES.md`), including whether `grade` equals Alauzet's gradation.
* Assembly of the hot kernels was not inspected; no `callgrind`. Speed numbers are wall-clock on a shared VM.
* Adaptation quality in 3D (nothing exists), threaded adaptation (M13, nothing exists), goal-oriented adaptivity (M10-M11, nothing exists).
