//! Ready-made photonic use cases built on the mode solvers and the polarization analysis, in the way `chips` builds the
//! Brillouin ring chips: a parameter struct, the mesh, the problem and a `solve` (with the solver passed in) that returns the usual `SolveResult`
//! (the Stokes parameters and the Poincare sphere of every mode are `nyx_engine::legacy::analyse` of it).
//!
//! * `soi`: the optical modes of a wire on 300 nm silicon-on-insulator (2D cross-section), with the analytic slab limit as a reference
//! * `bragg`: symmetric and asymmetric Bragg gratings in that wire (3D Floquet cell), stop band and coupling strength
pub mod bragg;
pub mod soi;

/// What the use cases solve with: a function from a mesh and a problem to the old-shape result. The kernels live in `nyx-ems` (ADR 45), which
/// `femcore` cannot name; the window and the tests pass `nyx_engine::legacy::solve`.
pub type Solver<'a> = &'a dyn Fn(&crate::mesh::SimplexMesh, &crate::problem::ProblemSpec) -> Result<crate::problem::SolveResult, String>;
