# Third-party notices (`solvers/`)

This file holds the licence texts and the reference list. **Which module took code, read code, or was only inspired by a paper is
spelled out per module in [`PROVENANCE.md`](PROVENANCE.md)**; the repository-level summary is in the root `NOTICES.md`.

`rsolvers` is a Rust port of ideas, defaults and control flow from two C
libraries. Their sources were consulted while writing this code (PETSc 3.19.6
`src/ksp/ksp/impls/{cg,gmres}`, `src/ksp/ksp/interface/iterativ.c`,
`src/ksp/pc/impls/{jacobi,sor}`; SLEPc 3.19.2
`src/eps/impls/krylov/{epskrylov.c,krylovschur/krylovschur.c}`,
`src/eps/interface/epssetup.c`), and the algorithms' observable behaviour
(default tolerances, norm and convergence tests, restart rules, selection
order) is reproduced deliberately so results can be compared one-to-one.

Both libraries are released under permissive BSD-2-Clause licenses, which allow
derivative works provided that the copyright notice, the list of conditions and
the disclaimer are retained. They are reproduced verbatim below, and the
`solvers/` folder is released under BSD-2-Clause (see `LICENSE` in this folder).

The code in this repository is not a mechanical translation: data structures,
kernels (fixed-chunk deterministic reductions, fused SpMV+dot, cache-blocked
basis operations) and the Rust API are new.

---

## PETSc

Source: <https://gitlab.com/petsc/petsc/-/blob/main/LICENSE> (BSD-2-Clause)

    Copyright (c) 1991-2021, UChicago Argonne, LLC and the PETSc Development Team
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without modification,
    are permitted provided that the following conditions are met:
    
    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this
      list of conditions and the following disclaimer in the documentation and/or
      other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    This license DOES NOT apply to any software that may be obtained via the
    --download-package option of the PETSc configuration. Each of those packages is
    covered by its own license(s).


---

## SLEPc

Source: <https://gitlab.com/slepc/slepc/-/blob/main/LICENSE.md> (BSD-2-Clause)

    Copyright (c) 2002-2023, Universitat Politecnica de Valencia, Spain
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    

## `femsolve`

`femsolve/` follows the same practice: ideas, defaults and observable behaviour were taken from
SLEPc 3.19.2 after reading `src/eps/impls/krylov/krylovschur/krylovschur.c`,
`src/eps/impls/krylov/epskrylov.c` (`EPSKrylovConvergence`), `src/eps/interface/epssolve.c`
(`EPSGetStartVector`), `src/eps/interface/epsdefault.c` (`EPSConvergedRelative`, purification),
`src/eps/interface/epssetup.c` (`EPSSetUpSort_Default`), `src/sys/classes/bv/interface/bvkrylov.c` and
`src/sys/classes/ds/impls/nhep/dsnhep.c`, to reproduce its restart, locking, start-vector and
convergence rules so iteration counts can be compared one to one. No code was copied; the
BSD-2-Clause notice for SLEPc and PETSc above applies to that derived behaviour.

Dependencies added by `femsolve`: `faer` 0.24.4 (MIT; it bundles notices for Eigen (MPL-2.0),
LAPACK (BSD), AMD and COLAMD (BSD) in its own `COPYING.*` files), `num-complex` (MIT/Apache-2.0),
optional `rayon` (MIT/Apache-2.0). The harness (`femsolve/harness`) links against the system PETSc/SLEPc
(BSD-2-Clause) and is not part of the library.

Phase 2 (K1, P2, E3, E2): for the Krylov solvers the PETSc 3.19 files `src/ksp/ksp/impls/cg/cg.c`
and `src/ksp/ksp/interface/iterativ.c` (default convergence test) were read (BSD-2-Clause, the notice above applies) to
reproduce PETSc's stopping rules and CG variants; GMRES/FGMRES follow the textbook algorithm with PETSc's documented
defaults. SLEPc's CISS, PETSc's GAMG and PETSc's GMRES/FGMRES are used only as black boxes through the public API in
`femsolve/harness` for comparison. The AMG, contour and Jacobi–Davidson code were written from the literature
(smoothed aggregation: Vaněk, Mandel, Brezina 1996; MIS-2 aggregation; block Sakurai–Sugiura/CIRR: Sakurai & Sugano 2003,
Ikegami–Sakurai–Nagashima 2010; Jacobi–Davidson: Sleijpen & Van der Vorst 1996, harmonic extraction: Paige–Parlett–Van der Vorst 1995,
Fokkema–Sleijpen–Van der Vorst 1998; AMS: Hiptmair & Xu 2007, Kolev & Vassilevski 2009; time stepping: Shampine & Reichelt 1997, Chung & Hulbert 1993, Verwer et al. 1999; nonlinear eigenproblems: Neumaier 1985, Beyn 2012, Asakura et al. 2009; TOAR: Lu, Su & Bai 2016; nonlinear solvers: Eisenstat & Walker 1996, Walker & Ni 2011, Brown & Saad 1990; incomplete factorisations: Saad, Iterative Methods for Sparse Linear Systems, ch. 10; Meijerink & van der Vorst 1977; Manteuffel 1980). PETSc `PCILU`/`PCICC`/`PCASM` are used as black boxes in the harness only. Optimisation: Nocedal & Wright (line searches, L-BFGS), Svanberg 1987 (MMA), Polak–Ribière; sensitivities: Lancaster 1964, Kuhn–Munkres (Hungarian algorithm); MINRES: Paige & Saunders 1975, Choi & Saunders; BiCGStab: van der Vorst 1992; LOBPCG: Knyazev 2001, Duersch et al. 2018; Chebyshev filtering: Zhou & Saad 2007, Winkelmann et al. 2019 (ChASE, described in the paper only); matrix exponential: Higham 2005, Sidje 1998 (EXPOKIT, described in the paper only); Schur complement preconditioning: Murphy, Golub & Wathen 2000. additive Schwarz: Smith, Bjørstad & Gropp; Dolean, Jolivet & Nataf. No hypre, PRIMME, BLOPEX, MUMPS or ChASE source has been read.

Licence rule agreed with Paul for the remaining items: read third-party sources only after confirming from the repository
itself that the licence is MIT, Apache-2.0 or BSD (hypre: MIT/Apache-2.0 dual; PRIMME: BSD-3-Clause), and record the check
here; never read LGPL or CeCILL code (BLOPEX, MUMPS). **Not yet checked:** nyx's own licence, and each repository above
(do the check from the repo's `LICENSE` file first).

Phase 3 (QMR, GCRO-DR, ADS, matrix-free operators, HODLR/ACA, NLEIGS, `sparse_lu`): **no third-party source code was read.** Everything
was written from the literature: QMR: Freund & Nachtigal 1991, Freund 1992 (complex symmetric), the Templates book (Barrett et al. 1994);
GCRO-DR: Parks, de Sturler, Mackey, Johnson & Maiti 2006; ADS: Hiptmair & Xu 2007, Kolev & Vassilevski 2012 (the `S, C, Π, C, Sᵀ` cycle
is the standard auxiliary-space one); sum factorisation: Orszag 1980, Deville–Fischer–Mund 2002 (GLL spectral elements); ACA: Bebendorf 2000;
HODLR solver: Ambikasaran & Darve 2013; LU: Gilbert & Peierls 1988; NLEIGS: Güttel, Van Beeumen, Meerbergen & Michiels 2014; block low-rank multifrontal LU: Amestoy, Ashcraft, Boiteau, Buttari, L'Excellent & Weisbecker 2015, with the multifrontal method of Duff & Reid 1983 and Liu 1992 (no MUMPS/STRUMPACK code was read; both are not permissively licensed or were not checked). The only external
code called is faer's AMD ordering (`faer::sparse::linalg::amd`) and its dense LU / triangular-solve / matmul kernels (MIT, already a dependency).
