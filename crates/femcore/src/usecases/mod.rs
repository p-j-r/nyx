//! Ready-made photonic use cases built on the mode solvers and the polarization analysis, in the way `chips` builds the
//! Brillouin ring chips: a parameter struct, the mesh, the problem and a `solve` that returns the usual `SolveResult`
//! (with `analysis` set: Stokes parameters and Poincare sphere of every mode).
//!
//! * `soi`: the optical modes of a wire on 300 nm silicon-on-insulator (2D cross-section), with the analytic slab limit as a reference
//! * `bragg`: symmetric and asymmetric Bragg gratings in that wire (3D Floquet cell), stop band and coupling strength
pub mod bragg;
pub mod soi;
