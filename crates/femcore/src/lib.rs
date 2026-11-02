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
//! * physics: `scalar` (electrostatics, heat), `solid` and its variants `solid_kz` (out-of-plane wave number), `solid_pml`,
//!   `solid_lrb` (absorbing boundaries), `solid_floquet` and `periodic` (periodic cells), `solid_axisym` (bodies of
//!   revolution); `emw` (2D waveguide modes), `emw_pml` (complex modes, tilted crystals), `emw_axisym`, `emw3d` (3D cavities)
//! * Brillouin scattering: `brillouin` (overlap integrals and gain), `ring_brillouin` (rings), `chips` (the four ring chips)
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
pub mod brillouin;
pub mod chips;
pub mod csparse;
pub mod eigen;
pub mod emw;
pub mod emw3d;
pub mod emw_axisym;
pub mod emw_pml;
pub mod expr;
pub mod fem;
pub mod fields;
#[doc(hidden)]
pub mod frame_adapter;
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
pub mod cellmesh;
pub mod par;
pub mod periodic;
pub mod polarization;
pub mod problem;
pub mod progress;
pub mod refine3d;
pub mod ring_brillouin;
pub mod scalar;
pub mod section;
pub mod solid;
pub mod solid_axisym;
pub mod solid_floquet;
pub mod solid_kz;
pub mod solid_lrb;
pub mod solid_pml;
pub mod solve;
pub mod solvers;
pub mod validate;
pub mod sparse;
pub mod usecases;
pub mod view;
pub mod xmlio;
