# Provenance of the code in `solvers/`: what was taken, what was only inspired, and how

This file answers, module by module, one question: **where did this code come from?** It complements
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) (licence texts and the full reference list). Read the vocabulary first; the
words are used strictly.

## Vocabulary

| Term | Meaning here |
|---|---|
| **Taken** | Source text (code or comments) copied or mechanically translated from another project. |
| **Read, then re-implemented** | A third-party source file was *read* (only permissively licensed ones: PETSc and SLEPc, BSD-2-Clause), and the behaviour — defaults, tolerances, stopping and restart rules, selection order — was reproduced in new Rust code, so iteration counts can be compared one to one. No text was copied; data structures, kernels and APIs are new. The BSD-2-Clause notice is retained anyway (`THIRD_PARTY_NOTICES.md`). |
| **From the literature** | Written from a published paper or book only; no third-party source code was read. The references are in the module's header comment and in the table below. |
| **Library call** | The code calls a dependency; it does not contain it. |
| **Black box** | PETSc/SLEPc are *run* through their public API by the comparison harness (`bench/`, `femsolve/harness/`) to produce reference numbers; none of their code is part of the Rust crates. |
| **Original** | Design and code are new (kernels, data layouts, test generators, studies). |

## What can and cannot be attested

* **Taken: nothing is recorded.** None of the three sessions that wrote this folder recorded copying third-party source text, and
  no file carries a third-party copyright header. Published numerical constants (for example the degree-13 Padé coefficients and
  `θ₁₃` of the matrix exponential, Higham 2005, in `femsolve/src/expm.rs`; Gauss–Lobatto and quadrature points computed in
  `matfree.rs`) are taken from the papers, as they are in every implementation.
* **Phase 3 (the most recent session: QMR, GCRO-DR, ADS, `matfree`, `lowrank`, `nleigs`, `sparse_lu`, `blr`, and their tests and
  studies)** read **no third-party source code at all**; everything was derived from the papers listed below. The only
  external code *called* is `faer` (MIT): its AMD ordering, dense LU, triangular solves and matrix product.
* **Earlier phases (`rsolvers`, and `femsolve` up to the nonlinear/optimisation items)**: the reading of PETSc 3.19 and
  SLEPc 3.19 files is recorded in `THIRD_PARTY_NOTICES.md` (exact file names). I could not re-audit those sessions line by line;
  the statements "no code was copied" and "not a mechanical translation" are theirs. What I could check in this session: no source
  file, C driver or script in `solvers/` carries a copyright header (`grep -i copyright`), and the structure differs from the
  upstream C (different data layout, fixed-chunk deterministic kernels, generic real/complex scalars). I did **not** diff the
  code against the upstream files. Anyone who needs certainty should compare the Krylov–Schur, CG and convergence-test code with
  `krylovschur.c`, `cg.c` and `iterativ.c`.
* **Nothing from nyx** was available to any session: no nyx source or data was read, copied or used. All test problems are
  synthetic (generated here, or with closed-form answers).
* **Licence rule** for any future reading: third-party sources are read only after confirming MIT / Apache-2.0 / BSD from the
  repository's own licence file; LGPL and CeCILL code (BLOPEX, MUMPS) is never read. Not yet checked: nyx's own licence.

## `rsolvers` (the earlier real-only crate, `solvers/src`, `solvers/tests`)

| Component | Category | Sources |
|---|---|---|
| `ksp::Cg`, `ksp::Gmres` | Read, then re-implemented | PETSc 3.19.6 `src/ksp/ksp/impls/{cg,gmres}`, `src/ksp/ksp/interface/iterativ.c` |
| `pc::Jacobi`, `pc::Sor` | Read, then re-implemented | PETSc `src/ksp/pc/impls/{jacobi,sor}` |
| `pc::Ilu0`, `pc::Ic0` | From the literature | Saad, *Iterative Methods for Sparse Linear Systems*, ch. 10; Meijerink & van der Vorst 1977 |
| `eps::krylov_schur` (Lanczos, thick restart, locking) | Read, then re-implemented | SLEPc 3.19.2 `src/eps/impls/krylov/{epskrylov.c,krylovschur/krylovschur.c}`, `src/eps/interface/epssetup.c` |
| `basis`, `kernels`, `vecops`, `par`, `sparse`, `rng`, tests, `bench/` | Original | fixed-chunk deterministic reductions, fused SpMV+dot, cache-blocked basis operations; `bench/` drives PETSc/SLEPc as a **black box** |

## `femsolve` (`solvers/femsolve/src`)

| Module | Category | Sources |
|---|---|---|
| `eps/` (block Krylov–Schur, complex Schur reordering), `problem.rs` | Read, then re-implemented | SLEPc 3.19.2 `krylovschur.c`, `epskrylov.c` (`EPSKrylovConvergence`), `epssolve.c` (`EPSGetStartVector`), `epsdefault.c` (`EPSConvergedRelative`, purification), `epssetup.c` (`EPSSetUpSort_Default`), `bvkrylov.c`, `dsnhep.c`; block variant and Schur code original (written because faer exposes no Schur reordering) |
| `ksp/cg.rs` (CG, COCG), convergence test in `ksp/mod.rs` | Read, then re-implemented | PETSc 3.19 `cg.c`, `iterativ.c` |
| `ksp/gmres.rs` (GMRES, FGMRES) | From the literature + PETSc's *documented* defaults | textbook Arnoldi/Givens GMRES; PETSc manual for restart 30, CGS |
| `ksp/{minres,bicgstab,cocr,idrs,qmr,gcrodr}.rs` | From the literature | Paige & Saunders 1975, Choi & Saunders; van der Vorst 1992; Sogabe & Zhang 2007; Sonneveld & van Gijzen 2008 and van Gijzen & Sonneveld 2011; Freund & Nachtigal 1991, Freund 1992, Barrett et al. 1994 (*Templates*); Parks, de Sturler, Mackey, Johnson & Maiti 2006 |
| `amg.rs`, `ams.rs`, `ads.rs` | From the literature | Vaněk, Mandel & Brezina 1996 (smoothed aggregation); Hiptmair & Xu 2007, Kolev & Vassilevski 2009 and 2012 (auxiliary-space Maxwell and H(div)). hypre was **not** read (`ams.rs`, `ads.rs` say so) |
| `ilu.rs`, `schwarz.rs`, `fieldsplit.rs`, `mixed.rs` | From the literature | Saad ch. 10, Meijerink–van der Vorst, Manteuffel 1980; Smith–Bjørstad–Gropp, Dolean–Jolivet–Nataf; Murphy–Golub–Wathen 2000; mixed-precision refinement is standard (Higham). PETSc `PCILU/PCICC/PCASM` are **black boxes** in the harness only |
| `jd.rs`, `lobpcg.rs`, `chebfsi.rs`, `contour.rs`, `nep.rs`, `toar.rs`, `nleigs.rs`, `sens.rs` | From the literature | Sleijpen & van der Vorst 1996, Paige–Parlett–van der Vorst 1995, Fokkema–Sleijpen–van der Vorst 1998; Knyazev 2001, Duersch et al. 2018; Zhou & Saad 2007, Winkelmann et al. 2019 (ChASE: described in the paper only); Sakurai & Sugano 2003, Ikegami–Sakurai–Nagashima 2010; Neumaier 1985, Beyn 2012, Asakura et al. 2009; Lu–Su–Bai 2016; Güttel, Van Beeumen, Meerbergen & Michiels 2014; Lancaster 1964, Kuhn–Munkres. SLEPc CISS is a **black box** for comparison |
| `ts/`, `snes.rs`, `opt.rs`, `sweep.rs`, `expm.rs` | From the literature | Shampine & Reichelt 1997, Chung & Hulbert 1993, Verwer et al. 1999; Eisenstat & Walker 1996, Walker & Ni 2011, Brown & Saad 1990 (PETSc `SNES` is a **black box**); Nocedal & Wright, Svanberg 1987, Polak–Ribière; reduced-basis / rational-Krylov idea for sweeps; Higham 2005 (Padé constants taken from the paper) and Sidje 1998 (EXPOKIT, described in the paper only) |
| `sparse_lu.rs` | From the literature + **library call** | Gilbert & Peierls 1988; calls `faer`'s AMD ordering (MIT) |
| `blr.rs`, `lowrank.rs` | From the literature + **library call** | Amestoy, Ashcraft, Boiteau, Buttari, L'Excellent & Weisbecker 2015 (BLR); Duff & Reid 1983, Liu 1992 (multifrontal); Bebendorf 2000 (ACA); Ambikasaran & Darve 2013 (HODLR solver). MUMPS and STRUMPACK were **not** read. Dense LU, triangular solves and matrix products are `faer` calls (MIT) |
| `matfree.rs` | From the literature | Orszag 1980 (sum factorisation); Deville–Fischer–Mund 2002 (Gauss–Lobatto spectral elements) |
| `maxwell.rs` (Nédélec and Raviart–Thomas cube generators) | From the literature + original | lowest-order edge/face elements are textbook (Nédélec 1980, Raviart–Thomas 1977); the mesh, numbering and checks are original |
| `direct.rs`, `faer_glue.rs`, `ordering.rs` | **Library call** + original | `faer` sparse Cholesky/LU (MIT); the nested-dissection ordering (BFS level sets, no refinement) and its tree are original |
| `scalar.rs`, `kernels.rs`, `vecops.rs`, `basis.rs`, `par.rs`, `csr.rs`, `bsr.rs`, `linop.rs`, `rng.rs`, `dense/`, `mm.rs`, `bundle.rs`, `bin/`, `examples/`, `tests/` | Original | the real/complex `Scalar` layer, deterministic kernels, formats, Matrix Market I/O, generators, studies |
| `harness/` (C drivers, Python runner) | Original; **black box** use of PETSc/SLEPc | links the system PETSc/SLEPc (BSD-2-Clause) through public API calls; not part of the library |

## Dependencies (called, not included)

`faer` 0.24.4 (MIT; ships notices for Eigen (MPL-2.0), LAPACK (BSD), AMD and COLAMD (BSD) in its own `COPYING.*`),
`num-complex` (MIT / Apache-2.0), `rayon` (MIT / Apache-2.0, optional `parallel` feature), `proptest` (MIT / Apache-2.0, tests only).

## How to check a claim

Every module's header comment names its sources (`grep -rn "Written from\|from the paper" solvers/femsolve/src`); the reading of
PETSc/SLEPc files is listed by file name in `THIRD_PARTY_NOTICES.md`; `Cargo.lock` lists every dependency.
