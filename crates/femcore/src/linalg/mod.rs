//! # solvers
//!
//! Solver building blocks for finite-element eigenproblems and linear/nonlinear systems, written once for real and
//! complex scalars. Ported from the staging crate `femsolve` of the inspired-rust-port repository (BSD-2-Clause; where
//! each module comes from is in `docs/ports/solvers/PROVENANCE.md`, the notices in `docs/ports/NOTICES.md`).
//!
//! * eigen: thick-restart block Krylov-Schur (`eps`), Jacobi-Davidson (`jd`), LOBPCG (`lobpcg`), Chebyshev filtering
//!   (`chebfsi`), contour integrals (`contour`), nonlinear and quadratic problems (`nep`, `toar`, `nleigs`), eigenvalue
//!   derivatives and mode tracking (`sens`); shift-and-invert with symbolic factorisation reuse (`direct`, `problem`)
//! * linear: Krylov methods (`ksp`: CG, COCG, GMRES, MINRES, BiCGStab, COCR, IDR(s), QMR, GCRO-DR), preconditioners
//!   (`pc`, `ilu`, `schwarz`, `fieldsplit`, `amg`, `ams`, `ads`), mixed precision (`mixed`), block low-rank and sparse LU
//!   (`blr`, `sparse_lu`, `ordering`), low-rank tools (`lowrank`)
//! * time, nonlinear, optimisation: `ts`, `snes`, `opt`, `expm`, reduced-order frequency sweeps (`sweep`)
//! * plumbing: `csr` (interchange matrix), `linop` (operator and preconditioner traits), `scalar`, `kernels`, `vecops`,
//!   `basis`, `par` (deterministic chunking on the shared nyx pool), `mm` and `bundle` (Matrix Market export of problems)
//!
//! These modules are independent of the older modules of `femcore` (`eigen`, `lobpcg`, `amg`, ...), which still run the
//! physics; the new solvers are reached through `femcore::linalg::...`. `docs/ports/solvers/NYX_INTEGRATION.md` is the
//! author's plan for swapping them in one by one; `docs/ports/solvers/STATUS.md` says what was measured. All parallel
//! kernels run on the pool of `crate::par` (`NYX_THREADS`) and give the same bits for any thread count.

pub mod ads;
pub mod amg;
pub mod amg_elastic;
pub mod ams;
pub mod basis;
pub mod arnoldi_real;
pub mod blr;
pub mod bsr;
pub mod bundle;
pub mod chebfsi;
pub mod contour;
pub mod csr;
pub mod dense;
pub mod direct;
pub mod driver;
pub mod eps;
pub mod expm;
pub mod faer_glue;
pub mod fieldsplit;
pub mod ilu;
pub mod jd;
pub mod kernels;
pub mod ksp;
pub mod linop;
pub mod lobpcg;
pub mod lowrank;
pub mod matfree;
pub mod maxwell;
pub mod mixed;
pub mod mm;
pub mod nep;
pub mod nleigs;
pub mod opt;
pub mod ordering;
pub mod par;
pub mod pc;
pub mod problem;
pub mod rng;
pub mod scalar;
pub mod schwarz;
pub mod sens;
pub mod snes;
pub mod sparse_lu;
pub mod sweep;
pub mod toar;
pub mod ts;
pub mod vecops;

pub use ams::{
    Ams, AmsFgmresOptions, AmsFgmresSolver, AmsFgmresStats, AmsOptions, AmsP2, AmsP2Options, AmsP2Topology,
    assemble_shifted_curl_curl, solve_ams_fgmres_many, solve_ams_p2_fgmres_many,
};
pub use csr::Csr;
pub use linop::{FnOperator, Identity, LinearOperator, Preconditioner};
pub use scalar::{C64, Scalar};
