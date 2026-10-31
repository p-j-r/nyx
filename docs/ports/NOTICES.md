# Notices

Summary for the whole repository. Details live next to the code; this file only says what applies where.

## Licence

* `solvers/` is released under **BSD-2-Clause** (`solvers/LICENSE`).
* `meshers/` (`femmesh`) is released under **BSD-2-Clause** (`meshers/femmesh/LICENSE`) with its own statements in
  [`meshers/PROVENANCE.md`](meshers/PROVENANCE.md) and [`meshers/femmesh/THIRD_PARTY_NOTICES.md`](meshers/femmesh/THIRD_PARTY_NOTICES.md).
  It was written by a different session; this file summarises only what `solvers/` contains. The meshing session reports no
  code taken from other projects and no dependencies.

## `solvers/`: how the code relates to other projects

The complete, module-by-module statement is [`solvers/PROVENANCE.md`](solvers/PROVENANCE.md); licence texts and references are in
[`solvers/THIRD_PARTY_NOTICES.md`](solvers/THIRD_PARTY_NOTICES.md). In short:

1. **Code taken from other projects: none is recorded.** No source file carries a third-party copyright header, and no
   session that wrote the folder recorded copying source text. Published numerical constants (for example the Padé
   coefficients of the matrix exponential, Higham 2005) are taken from their papers.
2. **Inspired by, after reading (PETSc and SLEPc 3.19, both BSD-2-Clause).** For the earlier real-only crate `rsolvers` and for the
   Krylov–Schur eigensolver, CG/COCG and the convergence tests of `femsolve`, specific PETSc/SLEPc source files were *read*
   (listed by name in `solvers/THIRD_PARTY_NOTICES.md`) and their observable behaviour — defaults, tolerances, stopping,
   restart, locking and selection rules — was reproduced in **new Rust code** so that iteration counts can be compared one to one
   with the C libraries. The data structures, kernels and APIs are new. The PETSc and SLEPc BSD-2-Clause notices are
   reproduced in `solvers/THIRD_PARTY_NOTICES.md` as a matter of course. The earlier sessions' claims could not be re-audited
   line by line by the session that wrote this notice; see `solvers/PROVENANCE.md` for exactly what was and was not checked.
3. **Written from the literature (no third-party source read).** Everything else: the linear solvers (MINRES, BiCGStab, COCR,
   IDR(s), QMR, GCRO-DR), the preconditioners (AMG, AMS, ADS, ILU/IC, Schwarz, field-split), the eigensolvers (Jacobi–Davidson,
   LOBPCG, Chebyshev filtering, contour/Beyn, TOAR, NLEIGS), time steppers, nonlinear solvers, optimisers, matrix functions,
   the low-rank tools (ACA, HODLR), the block low-rank multifrontal LU, the matrix-free spectral-element operators and the
   sparse LU with a user ordering. The papers are cited in each module's header comment and collected in `PROVENANCE.md`.
   hypre, PRIMME, ChASE, BLOPEX, MUMPS and STRUMPACK sources were **not** read.
4. **PETSc and SLEPc as black boxes.** The comparison harness (`solvers/bench`, `solvers/femsolve/harness`) *runs* the installed
   libraries through their public API to produce reference numbers. None of their code is part of the Rust crates.
5. **Libraries called (not included):** `faer` (MIT; it bundles notices for Eigen, LAPACK, AMD and COLAMD in its own files),
   `num-complex`, `rayon` and `proptest` (MIT / Apache-2.0). `faer` supplies the sparse Cholesky/LU, AMD ordering and the dense
   LU, triangular-solve and matrix-product kernels used by the solvers.
6. **Nothing from nyx** was available: no nyx source, data or licence was read or used. Whether nyx's own licence is compatible with
   the BSD-2-Clause code here has **not** been checked.

If a statement above turns out to be wrong for a file, treat the file as covered by the licence of its true source and fix this
notice together with `solvers/PROVENANCE.md`.
