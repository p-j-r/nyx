//! femcore: finite elements for the nyx application, in pure Rust. docs/HOW_IT_WORKS.md explains the whole chain in plain
//! language; this is the map of the modules.
//!
//! * entry and plumbing: `solve` (one entry point for every physics), `problem` (what is solved: the problem and result
//!   types), `validate` (refuses NaN/infinity), `progress` (stage log, cancel), `par` (the shared thread pools)
//! * meshes: `mesh` (simplex meshes, numbering of unknowns), `geomdrive` (geom models to meshes and named problems), `mesher` (2D shapes to triangles, extrusion), `mesher3d` and
//!   `refine3d` (3D solids to tetrahedra), `view` (surface for drawing), `cellmesh` (mixed-cell meshes)
//! * finite-element basics: `fem` (quadrature, shape functions, scalar matrices, constraints), `fields` (gradients, node
//!   averaging, plotting surface), `material` (anisotropic data, rotations), `expr` (user expressions)
//! * linear algebra: `sparse`, `csparse` (sparse matrices, real and complex), `linalg` (sparse factorisations, dense
//!   eigenproblems), `amg` (multigrid and CG), `eigen` (shift-and-invert Krylov), `lobpcg` (multigrid eigensolver)
//! * physics: `dfs` (electrostatics, heat), `sms` and its variants `sms_kz` (out-of-plane wave number), `sms_pml`,
//!   `sms_lrb` (absorbing boundaries), `sms_floquet` and `periodic` (periodic cells), `sms_axisym` (bodies of
//!   revolution); `ems` (2D waveguide modes), `ems_pml` (complex modes, tilted crystals), `ems_axisym`, `ems3d` (3D cavities)
//! * Brillouin scattering: `sbs` (overlap integrals and gain), `sbs_ring` (rings), `chips` (the four ring chips)
//! * ported libraries: `solvers` (Krylov, eigen, preconditioner, time-stepping and optimisation library for real and complex
//!   systems) and `meshers` (metric mesh adaptation, structured PML meshes, `.fmsh` files), see docs/ports
//! * files: `xmlio` (project XML), `h5io` (HDF5/XDMF results, feature `hdf5`)
//!
//! Lint note: `clippy::needless_range_loop` is allowed crate-wide. The numerical kernels (element matrices, orthogonalisation,
//! multigrid smoothers) index several arrays with one loop variable (`for i in 0..n { y[i] += a[i] * x[idx[i]] }`); the
//! indexed form is the one that matches the formulas, and rewriting 49 such loops as iterator chains would not change
//! the machine code. `too_many_arguments` is allowed on the few assembly/solver entry points that take a full problem
//! description (each marked where it applies).
#![allow(clippy::needless_range_loop)]
pub mod amg;
pub mod analysis;
pub mod sbs;
pub mod chips;
pub mod csparse;
pub mod eigen;
pub mod ems;
pub mod ems3d;
pub mod ems_axisym;
pub mod ems_pml;
pub mod estimate;
pub mod estimate3;
pub mod expr;
pub mod fem;
pub mod fields;
pub mod geomdrive;
#[cfg(feature = "hdf5-metno")]
pub mod h5io;
pub mod linalg;
pub mod material;
pub mod lobpcg;
pub mod mesh;
pub mod meshers;
pub mod mesher;
pub mod mesher3d;
pub mod meshing;
pub mod meshing_pml;
pub mod physics_size;
pub mod pml3;
pub mod cellmesh;
pub mod par;
pub mod periodic;
pub mod polarization;
pub mod problem;
pub mod progress;
pub mod refine3d;
pub mod sbs_ring;
pub mod dfs;
pub mod section;
pub mod sms;
pub mod sms_axisym;
pub mod sms_floquet;
pub mod sms_kz;
pub mod sms_lrb;
pub mod sms_pml;
pub mod solve;
pub mod solvers;
pub mod validate;
pub mod sparse;
pub mod usecases;
pub mod view;
pub mod xmlio;
