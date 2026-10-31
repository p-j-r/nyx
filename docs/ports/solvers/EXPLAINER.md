# What `femsolve` does, explained for an undergraduate

This note explains, without assuming more than first-year linear algebra and some programming, what was built in the
`femsolve` crate, why each piece exists, and when you would reach for it. The technical status (what was measured, what is
missing) is in `STATUS.md`; the lessons are in `LEARNINGS.md`. File names below are in `femsolve/src/`.

---

## 1. The big picture: from a physical problem to a matrix

Engineers simulate light in a chip, vibrations of a bridge, or the electric field in a cavity with the **finite element method
(FEM)**: cut the object into small cells, describe the unknown field by a few numbers per cell, and the physics turns into
**big sparse matrices** (most entries are zero because a cell only talks to its neighbours). Everything we need afterwards is
one of three jobs:

| Job | Written as | Question it answers |
|---|---|---|
| **Linear system** | `A x = b` | what is the field when a source `b` pushes on the structure? |
| **Eigenproblem** | `A x = λ B x` | at which frequencies (`λ`) does the structure like to vibrate/resonate, and with which shape `x`? |
| **Nonlinear problem** | `F(x) = 0` or `T(λ) x = 0` | the same, but the material depends on the frequency or the load |

The matrices have millions of rows, so we never invert them. We use **direct** methods (factor the matrix once, like Gaussian
elimination written for sparse matrices) or **iterative** methods (improve a guess step by step, using only matrix–vector
products). Most of the crate is the collection of such methods, written **once** for real *and* complex numbers (`scalar.rs`),
so the same code handles lossy materials, where numbers are complex.

Three recurring ideas to remember:

* **Krylov subspace.** Start with a vector `v` and repeatedly multiply by `A`: `v, Av, A²v, …`. The span of these is a cheap
  subspace that already contains good approximations of the answer. Almost every iterative method below builds one.
* **Preconditioner.** An easy-to-invert approximation `M ≈ A`. Solving `M⁻¹A x = M⁻¹b` instead makes the Krylov method converge
  in far fewer steps. Choosing `M` is where most of the craft is.
* **Shift-and-invert.** To find eigenvalues near a number `σ`, apply `(A − σB)⁻¹` instead of `A`: eigenvalues near `σ` become
  the largest ones, which Krylov methods find first.

---

## 2. Eigenvalue solvers: "find the resonances"

| Method (file) | In one sentence | Use it when |
|---|---|---|
| **Krylov–Schur** (`eps/`) | Build a Krylov subspace of `(A−σB)⁻¹B`, extract the best eigen-estimates, keep the good ones, restart. | The default: a handful to a few dozen modes near a chosen frequency. Matches SLEPc's behaviour step for step. |
| **Jacobi–Davidson** (`jd.rs`) | Instead of the exact inverse, solve a *correction equation* roughly at each step. | You cannot afford a factorisation (matrix-free) but want interior eigenvalues. |
| **Contour integral / Beyn / CIRR** (`contour.rs`, `nep.rs`) | Complex analysis: the integral of the resolvent `(zB−A)⁻¹` around a loop picks out exactly the eigenvalues inside the loop. | "Give me *all* modes in this region of the complex plane", easily run in parallel. |
| **LOBPCG** (`lobpcg.rs`) | Minimise the Rayleigh quotient with a three-term block recurrence. | Many *lowest* modes of a symmetric problem, with a preconditioner. |
| **Chebyshev filtering** (`chebfsi.rs`) | A polynomial that is tiny on the unwanted part of the spectrum and large on the wanted part, applied repeatedly. | Dozens of lowest modes of a symmetric problem, using only matrix products. |
| **TOAR** (`toar.rs`) | A quadratic eigenproblem `(λ²M + λC + K)x = 0` (damping!) solved without doubling the matrix size in memory. | Damped vibrations, absorbing boundaries. |
| **Residual inverse iteration** (`nep.rs`) | Newton-like polishing of one eigenpair of `T(λ)x = 0`. | You have a good guess and want 12 digits. |
| **NLEIGS** (`nleigs.rs`) | Replace `T(λ)` by a *rational interpolant*, which is an ordinary linear eigenproblem in disguise, and run a rational Krylov method. | Materials with square-root or pole terms in frequency (radiation conditions, Drude metals) and many modes wanted. |
| **Derivatives and mode tracking** (`sens.rs`) | `dλ/dp` from one extra solve; match modes between designs by their overlap. | Optimisation, sensitivity studies, following a mode while a parameter changes. |

**NLEIGS in a bit more detail** (new in this phase). Suppose the material makes the matrix depend on frequency in a messy way,
`T(λ) = Σ f_j(λ) A_j` with `f_j` like `√λ` or `1/(λ+iγ)`. We approximate every `f_j` on the region of interest by a rational
function through cleverly chosen points (the "Leja–Bagby" rule puts nodes where the approximation is hardest and poles near the
trouble spots). The approximated problem is equivalent to a bigger *linear* eigenproblem of size `d·n`. The trick that makes it
practical: we never build that big matrix; applying its shift-and-invert operator needs just **one** sparse factorisation of
the original size. We tested it on problems whose eigenvalues are known in closed form (33 of 33 found, errors ~1e-13).

---

## 3. Linear solvers: "find the field"

| Method (file) | Idea | Use it when |
|---|---|---|
| **CG / COCG** (`ksp/cg.rs`) | Conjugate gradients: optimal steps along "conjugate" directions. COCG is the variant for complex *symmetric* (not Hermitian) matrices. | Symmetric positive definite (CG), or lossy/PML Helmholtz-type (COCG). |
| **MINRES** | Like CG but safe for symmetric *indefinite* matrices. | Saddle-point and shifted problems. |
| **GMRES / FGMRES** (`ksp/gmres.rs`) | Pick the best correction in the whole Krylov space (minimal residual); needs memory for every basis vector, hence "restart every `m` steps". FGMRES allows a preconditioner that changes each step. | General nonsymmetric systems; the robust workhorse. |
| **BiCGStab, IDR(s), COCR** | Short-recurrence cousins: little memory, but erratic residuals. | Big nonsymmetric systems where GMRES's memory is too much. |
| **QMR** (`ksp/qmr.rs`, new) | Two-sided Lanczos with a *quasi-minimal* residual: short recurrences like BiCGStab but a smooth residual curve. Our complex-symmetric variant needs one matrix product per step. | Same niche as BiCGStab, when its jumpy residual is a problem; complex symmetric systems. |
| **GCRO-DR** (`ksp/gcrodr.rs`, new) | GMRES that **remembers** the slowest-converging directions (approximate eigenvectors of small eigenvalues) between restarts and between solves. | A *sequence* of similar systems: many right-hand sides, frequency sweeps, Newton steps. Measured: 632 vs 2 095 iterations. |

**Preconditioners** (the "approximate inverses"):

| Preconditioner (file) | Idea | Use it for |
|---|---|---|
| **Jacobi** (`pc.rs`) | Divide by the diagonal. | Baseline, and smoothing. |
| **ILU / incomplete Cholesky** (`ilu.rs`) | A sparse LU that drops small fill-in. | General sparse systems. |
| **AMG** (`amg.rs`) | Algebraic multigrid: solve on coarser and coarser "aggregated" problems and combine. | Elliptic/Helmholtz-type problems; grid-independent iteration counts. |
| **AMS** (`ams.rs`) | Multigrid for the *curl–curl* operator of Maxwell's equations, which has a huge null space (all gradients) that ordinary multigrid cannot see. Fix: also solve on the *nodes* via the discrete gradient. | Electromagnetic cavity/edge-element (Nédélec) problems. |
| **ADS** (`ads.rs`, new) | The same philosophy one dimension up, for the **div–div** operator on face elements (Raviart–Thomas): the null space is all curls, so we add a correction through the discrete curl (using AMS inside) and one through nodal vector fields. | Mixed/H(div) formulations (flux unknowns); iteration counts 15–21 instead of hundreds. |
| **Schwarz** (`schwarz.rs`) | Cut the domain into overlapping pieces, solve each exactly, add up. | Parallel machines; domain decomposition. |
| **Field-split / Schur** (`fieldsplit.rs`) | Treat coupled unknowns (velocity/pressure, E/H) block by block. | Coupled multiphysics. |
| **HODLR** (`lowrank.rs`, new) | Hierarchical low-rank matrices (see §5). | Dense-ish operators; loose-accuracy preconditioning. |
| **Multifrontal LU, optionally low-rank** (`blr.rs`, new) | Direct solver along a nested-dissection tree with dense fronts. | A fast exact solver for structurally symmetric matrices; the compressed version as a preconditioner for very large 3-D problems. |

**Direct solvers.** `direct.rs` wraps the `faer` sparse LU/Cholesky and remembers the *structure* of the factorisation when only
the numbers change (the shift `σ` in shift-and-invert). `mixed.rs` factors in single precision (twice as fast) and repairs the
accuracy by iterative refinement. `sparse_lu.rs` (new) is a small LU of our own that **accepts a user-chosen ordering**.

### Why the ordering matters (and what we found)

Eliminating unknowns in a bad order creates *fill-in* (new non-zeros), and a factorisation can become 5–10× more expensive.
Nested dissection cuts the mesh into halves recursively and numbers the cut last. On a 3-D grid it cut the fill 4.5× compared
with the natural order; on a 2-D grid the classical minimum-degree order is better, so the code picks per matrix. The `faer`
LU does not let us choose its ordering; our `sparse_lu.rs` does, at the price of being 1–2× slower in the factorisation itself
(it lacks the dense "supernode" kernels) but faster per solve.

---

## 4. Nonlinear problems, time, and design

| Piece (file) | Idea | Use it for |
|---|---|---|
| **Newton, JFNK, Picard, Anderson** (`snes.rs`) | Newton: linearise and solve repeatedly (with line search and Jacobian reuse). JFNK: replace the Jacobian by finite differences of the residual. Picard/Anderson: fixed-point iterations. | Nonlinear materials, nonlinear boundary conditions. Matches PETSc's SNES iteration counts. |
| **Time steppers** (`ts/`) | BDF 1–5 (variable step and order), Newmark/generalised-α, θ-method, Rosenbrock. | Transient heat, wave and structural dynamics; checked against known solutions and convergence orders. |
| **Matrix exponential** (`expm.rs`) | Compute `exp(tA)v` with a Krylov method, never forming `exp(tA)`. | Exact time stepping of linear systems. |
| **Frequency sweeps with reduced models** (`sweep.rs`) | Solve at a few well-chosen frequencies, build a small model from those solutions, evaluate it at every frequency. | Hundreds of frequencies for the price of a handful of solves. |
| **Optimisers + adjoints** (`opt.rs`) | L-BFGS, nonlinear CG, MMA; the *adjoint* trick gets the gradient with respect to thousands of parameters for the price of one extra solve. | Topology/shape optimisation. |

---

## 5. The new pieces of this phase, in more detail

### Low-rank structure and the fast direct solver (`lowrank.rs`)
Take the block of a matrix that couples two well-separated parts of a smooth problem: it is a big table of numbers, but its
**numerical rank is tiny**, so it can be stored as `U Vᵀ` with two thin matrices. *Adaptive cross approximation (ACA)* finds
such factors by looking at only a few rows and columns. A **HODLR** matrix applies this recursively: split the matrix in 2×2
blocks, compress the two off-diagonal blocks, recurse into the diagonal blocks. Matrix–vector products and even a **direct
solve** (by the Woodbury identity, applied level by level) then cost about `n log² n` instead of `n²`/`n³`. We measured a
16 384-row kernel matrix stored in 2 % of the dense memory and solved in six milliseconds.
Why it matters for sparse FEM: the dense blocks that appear on the *separators* of a nested-dissection factorisation have exactly
this structure (we measured ranks of only 8–10 for separators of 64–256 unknowns). Building a complete low-rank *sparse* solver
(BLR/HSS multifrontal) is the logical sequel and is listed as future work.

### Low-rank multifrontal sparse solver (`blr.rs`)
A sparse direct solver can be organised as a *tree*: cut the mesh in two with a separator (a line or a surface), solve the two
halves independently, and let the small dense problem on the separator tie them together; repeat recursively. Each tree node
works on a dense "front" (a small dense matrix), which dense linear-algebra libraries handle very fast. That is the
**multifrontal method**, and just using it (with a good nested-dissection ordering) made our solver 2–16× faster than the
`faer` sparse LU on the Helmholtz test matrices. The **low-rank** part adds the idea of the previous section: the big blocks inside
the fronts that couple far-apart pieces of a surface are stored as `U Vᵀ`. At the sizes we could run, only 3–6 % of the memory
was saved, because the blocks are not yet large compared with their rank; the literature reports larger savings for
problems with millions of unknowns. The honest summary is in `LEARNINGS.md` §32.

### Matrix-free high-order operators (`matfree.rs`)
High-order elements (degree `p` polynomials per cell) are accurate, but their element matrices are dense when the cells are
curved. **Sum factorisation** applies such an element operator one coordinate direction at a time — cost `~p⁴` instead of `~p⁶`
in 3-D — without ever storing the matrix. We measured it honestly: on a perfect brick mesh the assembled matrix stays so sparse
that ordinary storage wins; on a curved mesh with `p ≥ 5` matrix-free is faster *and* needs no matrix memory.

### H(div) and the exact sequence (`maxwell.rs`, `ads.rs`)
Gradient, curl and divergence form a chain: `div(curl u) = 0` and `curl(grad φ) = 0`. The discrete versions (matrices `G`, `C`,
`D`) must satisfy `CG = 0` and `DC = 0` exactly, and our generator checks this (plus the strong identity that the Nédélec
curl–curl matrix equals `Cᵀ M C`). Preconditioners for curl and div operators work precisely because they respect this chain.

---

## 6. How we know it works (verification in one page)

* **Same algorithm → same numbers (tier A).** For Krylov–Schur, GMRES, FGMRES, CG, ILU, ICC, Schwarz, BiCGStab and Newton we ran
  PETSc/SLEPc on the same matrices and compared *iteration counts and residual histories*: they agree digit for digit.
* **Different path to the same answer (tier B).** For the methods without a PETSc twin (Jacobi–Davidson, AMS, ADS, NLEIGS, HODLR,
  GCRO-DR…) we compare with closed-form answers (e.g. the exact eigenvalues of a 1-D Laplacian or a scalar-reducible nonlinear
  problem), dense reference solutions, or finite differences, and we recompute residuals independently of the solver.
* **Time-dependent or adaptive methods (tier C).** We check convergence *order*: halving the step must reduce the error by the
  predicted factor (4 for a second-order method).
* **Determinism.** Parallel reductions use fixed chunks and a fixed summation order, so results are bit-identical for 1 or 4 threads.
* **Honest negative results are recorded.** Example: the first matrix-free benchmark was 15–100× *slower* than assembled storage;
  the explanation (GLL box matrices are sparse) is in `LEARNINGS.md` §29 together with the case where it wins.

---

## 7. Cheat sheet: which tool for which problem?

| You have… | Reach for… |
|---|---|
| a few modes near a target frequency, matrix factorisable | `eps` (Krylov–Schur) with shift-and-invert |
| all modes inside a region | `contour` (or `nleigs` if frequency dependence is nonlinear) |
| many lowest modes of a symmetric problem | `lobpcg`, `chebfsi` |
| damped modes | `toar` |
| material with `√ω` or poles in `ω` | `nleigs` (+ `nep::residual_inverse_iteration` to polish) |
| a big SPD system | `Cg` + `Amg` |
| complex symmetric Helmholtz/PML | `Cg` (COCG) or `Qmr::solve_symmetric` + `Amg`/`Ilu` |
| nonsymmetric system | `Gmres` + `Ilu`; `Bicgstab`/`Idrs`/`Qmr` when memory is tight |
| many similar systems in a row | `GcroDr` |
| Maxwell edge elements | `Cg` + `Ams` |
| flux (face) elements | `Cg` + `Ads` |
| curved high-order elements | `matfree::SpectralBox::mapped` |
| a dense or kernel matrix with smooth structure | `lowrank::Hodlr` |
| LU where the ordering matters | `sparse_lu::SparseLu` with `LuOrdering::NestedDissection`, or (faster) `blr::BlrLu` with `tol = 0` |
| a nonlinear system | `snes` |
| time stepping | `ts` |

Not done (and why) is listed in `STATUS.md` → "Next steps": the flop-saving variants of the low-rank multifrontal solver, PETSc `TAO`/`TS`
comparison drivers, GPU and MPI (design notes only), and everything that needs the real nyx matrices, which were never available
to this work.
