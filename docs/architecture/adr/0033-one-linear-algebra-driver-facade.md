# 33. One linear algebra: the ported library behind one driver

Date: 2026-11-04

## Status

Accepted; implemented in phase 3 (bookmark `restructure-3-linalg`). Supersedes ADR 8 ("ported solver library staged, not wired"). The module rename `solvers` to `linalg` (plan section 2) is left to the last commit of the phase, see the status file.

## Context

femcore carried two linear-algebra stacks (R02, R03): the original one (`sparse`, `csparse`, `linalg`, `amg`, `eigen`, `lobpcg`: real and complex copies of the same code, `usize` indices, public matrix fields) and the ported library `solvers` (one implementation for `f64` and `Complex<f64>`, symbolic-factorisation reuse, block Krylov-Schur, LOBPCG, AMG, Krylov methods, determinism by fixed chunks, ADR 11) that nothing used. ADR 8 staged it; the plan (phase 3) wires it.

## Decision

* **One matrix type**: `solvers::csr::Csr<T>` for `f64` and `Complex<f64>`, with the assembly builder `Triplets<T>` (aliases `Coo`, `CCoo`, `RCsr`, `CCsr`) and the operations the kernels use (`submatrix`, `block_diag`, `sum_all`, `add_scaled`, `scaled`, `get`, `matvec`, ...). Fields are private; kernels use `indptr()`, `indices()`, `values()`, `values_mut()` and `Csr::from_raw`.
* **One door for the kernels**: `solvers::driver`. A kernel's solve stage calls only
  * `DirectSolver<T>`: Cholesky with LU fallback (AMD ordering, the lower triangle built without sorting), factor once, solve many; runs on the factorisation pool;
  * `eigs_near` and its kernel-shaped wrappers `eigs_sym`, `eigs_sym_indefinite`, `eigs_complex(_projected)`: block Krylov-Schur on the shift-and-invert operator, Hermitian or general, optional `B`-orthogonal null-space projector (discrete gradients of the curl-curl operator);
  * `eigs_general`: the real non-symmetric pencil of the 2D waveguide modes, in **real** block Arnoldi (`arnoldi_real`, ported unchanged from the old `eigen.rs`): the library's Krylov-Schur solves a real non-symmetric pencil only in complex arithmetic, where the real eigenvalues of this non-normal pencil acquire imaginary parts of 1e-8 to 1e-5 relative and large meshes do not reach the tolerance. A real Schur form with 2x2 blocks in the library would remove the need (not done);
    **Complex eigenvalues are never dropped** (some modes of a real pencil are leaky or have a complex propagation constant). The Ritz values come from a real matrix: a real Ritz value has an imaginary part of exactly zero, a complex one comes in a conjugate pair (the 2x2 blocks of the real Schur form). A value is taken as complex when `|Im theta| > 10 tol |theta|`, i.e. when its imaginary part exceeds the accuracy `tol` asks of the value itself; below that it is indistinguishable from a real value, and the cutoff only guards two nearly equal real values that the projection merged into a pair. A complex pair counts as a wanted mode only if it has converged (its residual `|OP x - theta x| / |theta x|` is below `tol`) and ranks among the nearest `nev + 2` eigenvalues, otherwise it is a Ritz artefact and is skipped as before. When one is found, `EigPairs::complex_found` is set and the solver stops; `ems::solve_emw` then repeats the study in the complex solver of `ems_pml` (complex vectors and n_eff). Contract test: `a_converged_complex_pair_of_a_real_pencil_is_reported_and_never_dropped`;
  * `lobpcg_guarded`: the library LOBPCG with a closure preconditioner and the stall guard (`LobpcgOptions::stall_gain`);
  * `dense_eigh`, `dense_gen_eigh`;
  * `EigenOptions { tol, block, ncv }` (serialisable; contract `linalg_options`).
* **Progress and cancel stay out of the library** and live in the driver: a stage line per factorisation and eigensolve, the cancel flag checked at every operator application. The library's operator trait cannot return an error, so a cancel leaves the Krylov loop as a panic payload (`CancelUnwind`) that the driver catches and turns into `Err("cancelled")`. No other panic is caught.
* **Elasticity multigrid** (rigid-body near-null-space, node blocks, Chebyshev smoother, P2 to P1 first coarse level) is `solvers::amg_elastic`, moved from the old `amg.rs` and `fem::AmgSystem` without a change of behaviour: the 3D eigenfrequency benchmark (33.7 s) depends on it. The hierarchy is built on the library `Csr<f64>`; only the single-precision copy used by the cycles (`Csr32`, 4-byte indices, `f32` values) is a module-private type.
* **Tolerance (gate G2, answered by paul: option B, loose tolerances everywhere, no reference platform)**: tolerances are options, not constants. `EigenOptions::tol` is what a kernel asks for (1e-9 to 1e-11 on the transformed residual, as before); `NYX_EIG_TOL` overrides it for experiments. Tests compare at tolerances justified by the physics they check (n_eff about 5e-5 relative and fields about 1e-3 of the maximum for the two `ems` silicon-wire goldens), and tighten only where the code needs it. No strict reference platform.
* **Migration** was kernel by kernel behind wrappers with the old signatures, with a switch that ran the old code (`NYX_OLD_LINALG`), and the goldens and closed-form checks as oracle; the switch, the bridge that converted matrix types and the old modules are deleted.

## Consequences

* The determinism rule (ADR 11) covers every kernel's solve stage; the old `par` helpers remain only for assembly.
* A kernel that needs a factorisation no longer sees faer.
* The cancel mechanism by unwinding is a deliberate exception to "no panics"; it is local to `driver.rs`.
* The two PML goldens `emw_pml_x` and `solid_kz_pml_z` freeze a mode the old solver dropped unconverged; the library converges it (3 against 4 and 4 against 5 solutions). Left as they are until paul decides.
* Retired: `sparse`, `csparse`, `linalg`, `amg`, `eigen`, `lobpcg`, and the examples `lu_scaling`, `pool_bench` (studies of the old stack).

## Code

`crates/femcore/src/solvers/driver.rs`, `solvers/csr.rs` and `csr/assembly.rs`, `solvers/arnoldi_real.rs`, `solvers/amg_elastic/`, `solvers/direct.rs`, `solvers/lobpcg.rs`; the kernels `dfs.rs`, `ems*.rs`, `sms*.rs`, `fem.rs`, `periodic.rs`.

## Contract

`contracts/linalg_options.v1.json` and `crates/femcore/tests/linalg_contract.rs`: options shape and defaults, driver surface by name, closed-form solves (Poisson 1D, Laplacian eigenvalues, a mildly non-normal real pencil with known real eigenvalues, dense helpers), cancel returning `"cancelled"`, repeat-run identical bits. `tests/axisym_closed_form.rs` checks the two axisymmetric kernels against closed forms.

## Solver path per stage, old stack against the new driver

Read from the code (old: `restructure-foundation`, new: this bookmark). "Timing" is `bench`/examples, old and new at the same machine load unless marked; **no number here is the required best-of-3 on an idle machine**.

| Stage (file) | Old path | New path | Specialisation lost or changed | Timing, old to new | Regression risk |
|---|---|---|---|---|---|
| dfs static, medium (`fem::solve_dirichlet`, `dfs.rs`) | `linalg::SparseSolver`: faer high-level supernodal Cholesky (lower triangle, AMD), LU if not positive definite, factorisation pool | `driver::DirectSolver`: faer low-level Cholesky through `ShiftedFactorization` (AMD, lower triangle built without sorting), same LU fallback, same pool | none in the method; the wrapper does more work per call (pattern assembly, a copy of the inputs) | electrostatics 2D 0.08 to 0.11 s (before the pattern fix; fix not re-measured) | **generic replacement, small solves slower** |
| dfs static, large (AMG-CG, `fem::AmgSystem`) | old `amg.rs` (smoothed aggregation, rigid-body/constant near-null-space, Chebyshev smoother, P2 to P1 coarse level) + `pcg` | same code in `solvers::amg_elastic`, same options | none; hierarchy built on a private copy of the old matrix type, `pcg` multiplies with the library SpMV (fixed-chunk sums) | solid static 3D 3.1 to 2.8 s; AMG-CG 40x8x8 3.7 to 3.9 s; `mem_profile 12 12 12 p2 static auto` baseline 1.28 s, best of 3 now 0.80 s (load unknown) | none seen |
| dfs transient (`dfs.rs`) | dense generalised eigendecomposition up to 2500 free dofs; above, BDF2 with two Cholesky factorisations or AMG-CG | `dense_gen_eigh` (same code), then `DirectSolver` twice or AMG-CG | none | not timed | as dfs static |
| sms static (`sms::solve_static` through `fem::solve_spd`) | as dfs static, medium or AMG-CG by `spec.solver` and size | as dfs static | none | solid static 2D 0.26 to 0.39 s; 3D direct 12.3 to 13.9 s (before the pattern fix) | **small and 3D direct solves slower** |
| sms eigenfrequencies, shift-invert (`sms::solve_real_modes`, `eigs_sym`) | `eigen::eigs_sym`: block Arnoldi/Lanczos with a growing basis, `SparseSolver` (Cholesky when sigma <= 0), block `min(nev, 8)` | `driver::eigs_sym`: block Krylov-Schur with restarts, same block, Cholesky when sigma <= 0 | algorithm changed (restarted instead of growing basis); convergence test is the same transformed residual | solid 6 eigenfrequencies 2D: 1.73 to 1.06 s (basis 66 to 21 vectors) | none seen |
| sms eigenfrequencies, LOBPCG + multigrid (`sms::solve_real_modes`, large 3D, sigma <= 0) | `lobpcg::lobpcg_guarded` + `AmgSystem::precondition_free`; Auto falls back to the factorisation when multigrid is slow or stalls | library `lobpcg` (soft locking, B-orthonormal basis) + the same `AmgSystem`; same thresholds and fallback | LOBPCG implementation changed; stall guard approximated (window 10, 10 % gain, was best-of-window comparison); residual definition differs slightly | solid 6 eigenfrequencies 3D: 58 to 41 s under load; **old idle baseline 33.7 s, not matched against an idle new run** | **possible regression until measured idle; stall behaviour on slender beams untested** |
| sms_pml, sms_floquet, sms_lrb (quadratic, 2n pencil), axisym (`sms_axisym`, `ems_axisym`) | `csparse::eigs_complex`: complex block Arnoldi with a basis up to `max_dim` (400 to 800), one `CSolver` (faer high-level sparse LU, factorisation pool); multi-RHS solves spread over threads when n > 5000 | `driver::eigs_complex`: block Krylov-Schur (restarts, library default basis), `ShiftedFactorization` LU, solves sequential | **`max_dim` is ignored** (library default basis); parallel multi-RHS solves lost (the library measured faer's parallel multi-RHS sparse solve as pathologically slow, sequential is its default); algorithm changed | not timed (no benchmark covers these); goldens, closed forms and determinism tests pass | **generic replacement of a tuned solver; large complex problems (n > 5000, many solves) may be slower** |
| ems real 2D (`ems::solve_mode_pairs`) | `eigen::eigs_general`: real block Arnoldi, `SparseSolver` LU | `arnoldi_real::eigs_general`: same code on the new matrix type, `DirectSolver` LU | none; complex pairs now reported and re-solved in the complex path (new) | EM modes 2D 2.18 to 2.02 s | none seen |
| ems_pml, tilted crystal (`ems_pml::solve_mode_pairs_complex`, including the fixed-point loop with 1-mode solves) | `eigs_complex` as above, a fresh `CSolver` per refinement step | `driver::eigs_complex`, a fresh factorisation per step | as the complex row | not timed | as the complex row |
| ems3d real (`ems3d::solve_emw_modes`) | `eigs_sym_projected` with the solves of a `SparseSolver` LU of `K - sigma M` and a Cholesky projector of `G^T M G` | `driver::eigs_sym_indefinite` (Krylov-Schur, LU) with the same projector (`DirectSolver`, Cholesky) | algorithm changed as in the shift-invert row | `cavity_bench 0.012 2`: baseline 4.34 s (4.43 at phase 1), best of 3 now 3.32 s (load unknown) | none seen |
| ems3d complex (PML, Floquet) | `csparse::eigs_complex_projected`, `CSolver` projector | `driver::eigs_complex_projected`, `DirectSolver` LU projector | as the complex row | not timed | as the complex row |
| 3D direct (static solid, `fem::solve_dirichlet`) | faer default ordering (AMD) | AMD (the library's `Auto`, which also tries nested dissection, was measured by its authors to cut 3D fill by 30 % but cost a second symbolic analysis that doubled small solves, so it is not used) | **nested dissection for large 3D is available and not used** | 12.3 to 13.9 s | **a size-dependent ordering could win here** |

Not adopted: ILU/ICC, Schwarz, field-split, AMS/ADS, Jacobi-Davidson, contour, mixed precision, BLR/HODLR factorisations of the library; nothing selects them per problem type yet. Open: best-of-3 idle measurement of every row marked "not timed" or "under load"; a benchmark that covers PML/axisym complex problems of realistic size (none exists today); the ordering rule for large 3D direct solves.

## Parallel solves of a block (complex stages), measured 2026-11-04

Why solves were sequential: the library's default `ThreadPolicy::solve` is `Sequential` because faer 0.24.4's *own* parallel triangular solve with more than one right-hand side was measured by the library authors at about 1.4 s per call whatever the block width (64 000-unknown 3D problem), against 12 ms per vector sequentially, and a single-vector solve gains nothing from threads (`direct.rs`, `LEARNINGS.md` 1, 3, 24; measured on a 4-vCPU Linux VM). The old complex eigensolver did something else: it ran the *sequential* solve on groups of columns in separate threads (`CSolver::solve_many`, for n > 5000). That is independent solves against one shared factorisation, not faer's parallel solve, so the library's measurement does not apply to it.

Switch: `EigenOptions::par_solve_min_n` (contract `linalg_options` 1.1; 0 = one thread). The complex wrappers (`eigs_complex`, `eigs_complex_projected`) default to 5000, the old rule; `NYX_PAR_SOLVE_MIN_N=<n>` overrides it (0 turns it off). Real shift-invert stays sequential, as the old `eigs_sym` was.

Measurement (`examples/complex_bench.rs`: `ems_pml` rectangular guide with an absorbing layer, order 2, 4 modes, best of 3 runs inside the process, machine otherwise quiet but not guaranteed idle), seconds:

| unknowns | old stack (foundation) | new, sequential solves | new, column-parallel solves |
|---|---|---|---|
| 6 779 | 0.175 | 0.178 | 0.175 |
| 19 837 | 0.539 | 0.588 | 0.515 |
| 35 465 | 0.963 | 1.033 | 0.914 |

The effective indices agree to the printed 6 digits in all three. Sequential is 6 to 7 % slower than the old stack from about 20 000 unknowns; column-parallel solves are 5 % faster than old and 12 % faster than sequential, and equal below 7 000 unknowns where no thread is used. Decision: keep column-parallel (n >= 5000) as the default for the complex stages. Limits of the evidence: one problem family (`ems_pml`), sizes up to 35 000 unknowns and about one second per solve; no `sms_pml`, axisym or `ems3d` complex case was timed, and nothing larger than 35 000 unknowns.

## Fill-reducing ordering: AMD everywhere (measured 2026-11-04)

The library's `FillOrdering::Auto` (a second symbolic analysis with nested dissection, keep the ordering with fewer factor non-zeros) was tried as a size rule (nested dissection for large 3D, AMD otherwise) and lost. `bench`, `NYX_ORDERING` forcing the ordering for every solve, two runs each, loaded machine (other sessions were running, so absolute times vary by 2x between runs of the same setting):

| case | AMD run 1 / run 2 | Auto run 1 / run 2 |
|---|---|---|
| solid static 3D direct, 205 479 unknowns | 4.68 / 10.92 s | 13.46 / 17.78 s |
| electrostatics 2D, 38 697 | 0.09 / 0.07 s | 0.14 / 0.13 s |
| solid static 2D, 77 394 | 0.37 / 0.23 s | 0.41 / 0.41 s |

AMD is faster in every pair. The library's 30 % less fill and 2.4x faster factorisation were measured on a scalar 3D Laplacian with a sequential factorisation; they do not carry over to P2 elasticity with the parallel factorisation (the library's own notes record a parallel Cholesky with nested dissection slower than sequential). Rule: AMD for every size (as the old stack); `NYX_ORDERING=amd|auto|nd` overrides it for experiments. With AMD the 3D direct solve (4.7 s in the best run) is at the phase 0 idle baseline (4.84 s).
