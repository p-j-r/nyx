# 117. Indefinite supernodal `LBL^H`/`Sym LBL^T` direct solver and $P_2$ Nédélec Hiptmair–Xu AMS

Date: 2026-11-22

## Status

Accepted. Amends ADR 8, 11, and 33.

## Context

Shifted 3D $P_2$ Nédélec curl-curl eigenproblems (`K - sigma M`) and driven frequency solves (`K + (-omega^2 + i omega gamma) M`) are symmetric indefinite (real-symmetric, complex-Hermitian, or non-conjugate complex-symmetric under PML/lossy dielectrics/IBC). Factorizing them with unsymmetric `faer` COLAMD LU stores both $L$ and $U$ with off-diagonal fill pivoting, consuming $\sim 4\times$ more factor nonzeros and RAM than symmetric factorization. Additionally, scaling beyond direct-solver RAM budgets requires a mesh-independent $P_2$ Nédélec auxiliary-space preconditioner (`AmsP2`).

## Decision

* **Ordering-aware supernodal `LBL^H` / `Sym LBL^T` direct solver (`femcore::linalg::{direct,driver}`):**
  * Added `FactorBackend` (`CholeskyAmd`, `CholeskyNd`, `LbltAmd`, `LbltNd`, `LuColamd`, `LuAmd`, `LuNd`, `BlrNd`) and `factor_nnz` on `FactorStats`.
  * Added $O(\text{nnz})$ transpose-map structural/Hermitian/complex-symmetric detection inside `ShiftedFactorization`.
  * Wired `faer` supernodal intranodal Bunch–Kaufman $LBL^H$ (`SymbolicCholesky::factorize_numeric_intranode_lblt`, lower triangle only, honoring `FillOrdering::{Amd, NestedDissection, Auto}`) for real-symmetric and complex-Hermitian indefinite systems, and implemented non-conjugate supernodal + simplicial Bunch–Kaufman $LBL^T$ over `faer`'s `SymbolicCholesky` and BLAS-3 dense kernels for complex-symmetric indefinite systems.
  * Guarded both factorizations with finiteness and a deterministic 1-vector backward residual check (`<= 1e-8`), falling back transparently to COLAMD LU on singular or inter-supernode pivot loss while preserving `symbolic_count == 1` across shifts.
  * Enabled `FillOrdering::Auto` (honoring `NYX_ORDERING=amd|nd|auto`) in `driver::eigs_near`, automatically benefiting all `ems::cavity` and `ems::driven` solves with $-76\%$ to $-77\%$ factor nonzeros/RAM and $4\text{–}8\times$ faster factorization on 3D $P_2$ Nédélec systems.
* **$P_2$ Nédélec Hiptmair–Xu AMS (`femcore::linalg::{ams,driver,maxwell}`):**
  * Delivered `AmsP2`, `AmsP2Topology`, `AmsP2Options` (symmetric multiplicative `0 - G2 - P1(AMS) - G2 - 0` cycle with $2\times 2$ block Gauss–Seidel smoothing, full $P_2$ scalar gradient AMG, and coarse $P_1$ `Ams`), `solve_curl_curl_p2_iterative`, `eigs_curl_curl_p2_jd`, `IterativeShiftInvertP2`, and `nedelec_p2_cube` cleanly inside `femcore::linalg`.
  * Kept `nyx-ems/src/cavity.rs` on `ShiftedFactorization` (now backed by `LBL^H`/`Sym LBL^T`) rather than wiring `AmsP2` prolongation into `cavity.rs` yet.

## Consequences

* All existing `ems::cavity` and `ems::driven` shift-invert solves automatically cut factorization memory by $\sim 76\%$ and run $4\text{–}8\times$ faster without changing any kernel API, while `femcore::linalg` exposes a verified, mesh-independent $P_2$ Nédélec AMS preconditioner and Jacobi–Davidson / FGMRES drivers for future large-scale iterative EM integration.

## Code

`crates/femcore/src/linalg/{direct,driver,ams,maxwell,sweep,blr,mixed,sparse_lu}.rs`, `crates/femcore/tests/linalg/fr16.rs`.
