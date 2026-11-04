//! femcore: finite elements for the nyx application, in pure Rust. docs/HOW_IT_WORKS.md explains the whole chain in plain
//! language; this is the map of the modules.
//!
//! * entry and plumbing: `solve` (one entry point for every physics), `problem` (what is solved: the problem and result
//!   types), `validate` (refuses NaN/infinity), `progress` (stage log, cancel), `par` (the shared thread pools)
//! * meshes: `mesh` (the meshes: `SimplexMesh` and its numbering of unknowns, the native per-cell-type mesh `native` and `cells`, `view`, `io` for `.fmsh`, `quality`, `fixtures`)
//! * meshing: `meshing` (geom models to meshes: `MeshPlan`, `geomdrive`, `mesher` 2D shapes to triangles, `mesher3d` and `refine3d` 3D solids to tetrahedra, structured PML `pml`, `pml3`, `pml_box`, `physics_size`, `validate`)
//! * adaptation: `adapt` (metric adaptation 2D and 3D, error estimators `estimate`, `estimate3`)
//! * finite-element basics: `fem` (quadrature, shape functions, scalar matrices, constraints), `fields` (gradients, node
//!   averaging, plotting surface), `material` (anisotropic data, rotations), `expr` (user expressions)
//! * linear algebra: `linalg` (the one library: `Csr<T>` matrices, direct and Krylov solvers, eigensolvers, multigrid;
//!   `linalg::driver` is the door the kernels use)
//! * physics: `dfs` (electrostatics, heat), `sms` and its variants `sms_kz` (the shared cross-section solid kernel, today for rings), `sms_pml`,
//!   `sms_lrb` (absorbing boundaries), `sms_floquet` and `periodic` (periodic cells), `sms_axisym` (bodies of
//!   revolution); the electromagnetic kernels are in `nyx-ems` (ADR 45); `edge` (the 2D edge-element basis) and `pml` (coordinate
//!   stretching) are what they share with `optical` and the solid kernels
//! * optical mode sets: `optical` (modes with the data to evaluate their fields; used by `analysis` and `polarization`)
//! * ported libraries: `linalg` (see above; Krylov, eigen, preconditioner, time-stepping and optimisation library for
//!   real and complex systems) and the mesher port in `mesh`, `meshing` and `adapt`, see docs/ports
//! * files: `xmlio` (project XML), `h5io` (HDF5/XDMF results, feature `hdf5`)
//!
//! Lint note: `clippy::needless_range_loop` is allowed crate-wide. The numerical kernels (element matrices, orthogonalisation,
//! multigrid smoothers) index several arrays with one loop variable (`for i in 0..n { y[i] += a[i] * x[idx[i]] }`); the
//! indexed form is the one that matches the formulas, and rewriting 49 such loops as iterator chains would not change
//! the machine code. `too_many_arguments` is allowed on the few assembly/solver entry points that take a full problem
//! description (each marked where it applies).
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]
pub mod analysis;
pub mod optical;
pub mod edge;
pub mod expr;
pub mod fem;
pub mod fields;
#[cfg(feature = "hdf5-metno")]
pub mod h5io;
pub mod material;
pub mod adapt;
pub mod mesh;
pub mod meshing;
pub mod par;
pub mod periodic;
pub mod pml;
pub mod polarization;
pub mod problem;
pub mod props;
pub mod progress;
pub mod section;
pub mod sms;
pub mod sms_axisym;
pub mod sms_floquet;
pub mod sms_kz;
pub mod sms_lrb;
pub mod sms_pml;
pub mod solve;
pub mod testkit;
pub mod linalg;
pub mod validate;
pub mod usecases;
pub mod xmlio;

