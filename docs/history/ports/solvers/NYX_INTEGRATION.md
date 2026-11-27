# Moving `femsolve` into nyx `femcore`

**How sure I am.** Everything about nyx here comes from `solver_port_list.md`; the nyx repository
was not available to the session that wrote this (the brief left its name blank). The
`femsolve` side was built and tested in this repository. Treat every nyx assumption below as
unverified until it has been checked against the code.

## What is where

| `femsolve` file | Role | nyx counterpart (per the port list) | Action in nyx |
|---|---|---|---|
| `scalar.rs`, `kernels.rs`, `vecops.rs`, `basis.rs` | one implementation for `f64` and `Complex<f64>`; fixed-chunk deterministic kernels; block Krylov basis | parts of `eigen.rs`, `csparse.rs`, `lobpcg.rs` (real and complex copies) | new code; can replace the duplicated real/complex orthogonalisation loops |
| `par.rs` | deterministic chunk parallelism on rayon | `par.rs` (one shared rayon pool) | **replace the two functions `chunk_map` / `for_chunks_mut` by calls into nyx's pool**; nothing else in the crate talks to rayon |
| `csr.rs` | interchange matrix type, generators | nyx assembly types (`CCoo`, CSC) | write an adapter (`Csr::from_triplets`) or implement `LinearOperator` for nyx's own matrix type; the solvers only need `LinearOperator` |
| `linop.rs` | operator/preconditioner traits with fusion hooks | none | adopt as the solver interface |
| `direct.rs` | shifted sparse factorisation with symbolic reuse and a thread policy | `linalg.rs` (faer Cholesky), `csparse.rs` (faer LU) | use as the `ShiftedSolver` backend, or implement `ShiftedSolver` on top of nyx's existing wrappers |
| `eps/` | block Krylov–Schur | `eigen.rs` (block shift-and-invert Arnoldi) | replace after the validation below |
| `problem.rs` | `ShiftInvertEigenSolver` trait: "eigenpairs near a shift" | none (models call `eigen.rs` directly) | the interface models should call |
| `mm.rs`, `bundle.rs`, `bin/fs_*.rs` | Matrix Market export, problem bundles, CLIs | none | use to export matrices (item V1) |
| `ksp/`, `pc.rs`, `ilu.rs`, `schwarz.rs`, `fieldsplit.rs`, `amg.rs`, `ams.rs`, `mixed.rs`, `bsr.rs` | Krylov solvers (CG, COCG, GMRES, FGMRES, MINRES, BiCGStab, COCR, IDR(s)), preconditioners (Jacobi, ILU(k)/ILUT/IC, additive/restricted Schwarz, field-split/Schur, smoothed-aggregation AMG, AMS for H(curl)), mixed-precision LU, block/symmetric CSR | `pcg.rs`/`amg.rs`/`lobpcg.rs`, FEM drivers | new code behind `LinearOperator`/`Preconditioner`; replace nyx's own AMG-PCG only after the comparison in `STATUS.md` item 2; AMS needs the discrete gradient and edge vectors from `emw3d.rs` |
| `jd.rs`, `lobpcg.rs`, `chebfsi.rs`, `contour.rs`, `nep.rs`, `toar.rs`, `sens.rs` | interior (JD), lowest (LOBPCG, Chebyshev filtering) and region (contour, Beyn) eigensolvers, nonlinear and quadratic eigenproblems, eigenvalue derivatives and mode tracking | `eigen.rs`, `lobpcg.rs`, quadratic linearisations | alternatives selectable per model once validated on the exported matrices |
| `ts/`, `sweep.rs`, `snes.rs`, `opt.rs`, `expm.rs` | time steppers (BDF, Newmark/gen-α, θ, Rosenbrock), reduced-order frequency sweeps, nonlinear solvers (Newton, JFNK, Anderson), optimisers (L-BFGS, CG, MMA), matrix-exponential action | transient and sweep code in `femcore` | use when transient, dispersive or optimisation work starts |
| `nleigs.rs`, `ads.rs`, `matfree.rs`, `lowrank.rs`, `sparse_lu.rs`, `blr.rs`, `ksp/qmr.rs`, `ksp/gcrodr.rs` | NLEIGS for frequency-dependent materials with branch cuts or poles; ADS for H(div); sum-factorised operators for curved high-order hexahedra; HODLR/ACA; LU with a chosen ordering; QMR; recycling GMRES | none | `nleigs`: when `Σ f_j(ω) A_j` has `√` or pole terms and many modes are wanted (needs the singularity locations). `ads`: only for face-element formulations (needs the discrete curl and face-flux interpolation from the mesh). `matfree`: needs nyx's element geometry (replace `mapped`'s coordinate input). `sparse_lu`: when complex-symmetric PML matrices should use the nested-dissection ordering. `blr`: **try first as a drop-in replacement for the faer LU path on structurally symmetric matrices** (`ShiftedBlrLu` with `tol = 0` was 2–16× faster to factor on 2-D/3-D Helmholtz-type matrices); the compression (`tol > 0`) only pays for much larger 3-D problems. `gcrodr`: frequency sweeps and Newton loops with slowly varying matrices |
| `harness/` | SLEPc comparison (C driver, Python runner) | none | runs on Linux/WSL next to the exported bundles |

## faer version

`femsolve` uses faer **0.24.4**. All `faer` calls are in two files:

* `faer_glue.rs`: dense Hermitian eigendecomposition (`Mat::self_adjoint_eigen`) and the
  `FaerScalar` conversions.
* `direct.rs`: the low-level sparse API (`factorize_symbolic_cholesky`, `factorize_symbolic_lu`,
  `factorize_numeric_llt`, `factorize_numeric_lu`, `LltRef::new`, `LuRef::new_unchecked`,
  `solve_in_place_with_conj`, `solve_in_place_scratch`). The low-level API is used because it
  takes the parallelism and the scratch memory per call; the high-level `Lu`/`Llt` types read a
  *global* parallelism setting, which cannot be set per phase (see `LEARNINGS.md`, item 1).

If nyx pins a different faer, adapt those two files only. If it cannot move to 0.24, the
high-level `SymbolicLu::try_new` + `Lu::try_new_with_symbolic` pair (faer 0.19 and later) gives
the same symbolic reuse; the thread-policy part then has to be done by toggling
`faer::set_global_parallelism` around solves (not thread-safe if several threads factorise
concurrently; acceptable for one solve at a time).

## Suggested integration order

1. **Adapters.** Convert one nyx matrix pair to `Csr` and write it with `Bundle::write`
   (`A.mtx`, `B.mtx`, `v0.mtx`, `meta.txt`). Do this for one optical and one elastic problem
   from a Brillouin run; this is also what port-list item 10 asks permission for.
2. **Reproduce.** `fs-eps <bundle> --nev 24 --sigma-re <shift> ...` and compare against the
   existing result: the port list asks for the 24-mode, 102k-unknown elastic result to 1e-9 with
   a basis of about 60 (`--ncv 60`). Note the residual definition: `fs-eps` prints the
   independently recomputed `‖Ax−λBx‖/(|λ|‖x‖)`.
3. **Compare with SLEPc** on the same bundle: `python3 femsolve/harness/run.py` (see its header
   for the one-time setup). Add the bundle as a case in `SUITES`.
4. **Plug in the thread policy.** Check what nyx's block solves do: if they go through faer's
   *global* parallel solve with more than one right-hand side, they may be hitting the
   multi-RHS cliff described in `LEARNINGS.md` item 1 (a fixed cost of seconds per call on a
   64k-unknown 3-D problem). This is cheap to check with `examples/solve_scaling.rs`
   adapted to nyx's own matrices.
5. **Swap the eigensolver** behind `ShiftInvertEigenSolver` and keep the old one selectable until
   the nyx analytic tests (closed forms, rod velocity, bar bending, Q = 5 complex bar, no-PML
   limits) pass with both.
6. **Replace `par.rs` glue** (two functions) so the kernels run on nyx's pool, then run
   `tests/determinism.rs` with `--features parallel`.

## Things that will differ from the existing nyx code

* Restarts, not growth: memory is `n·(ncv + block)` scalars instead of up to 3·n·600 doubles.
* Hermitian generalised problems return **B-orthonormal** eigenvectors (mass-normalised modes);
  general (non-Hermitian) problems return 2-norm unit vectors.
* The default convergence test of `KrylovSchurShiftInvert` is the relative residual of the
  *original* problem, `‖Ax−λBx‖/|λ| < tol`. SLEPc's own default (`‖r‖/|θ|` in the transformed space) is
  available via `theta_convergence` and is what the parity runs use.
* Real non-symmetric operators are not supported by the general kind (Schur form with 2×2 blocks
  is not implemented); promote them to complex.
