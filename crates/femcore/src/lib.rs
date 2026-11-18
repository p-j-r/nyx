//! femcore: finite elements for the nyx application, in pure Rust. docs/HOW_IT_WORKS.md explains the whole chain in plain
//! language; this is the map of the modules.
//!
//! * entry and plumbing: `props` (the data of a problem that is not specific to one physics: materials per domain, Floquet pairs, linear solver),
//!   `validate` (refuses NaN/infinity), `progress` (stage log, cancel), `par` (the shared thread pools)
//! * meshes: `mesh` (the meshes: `SimplexMesh` and its numbering of unknowns, the native per-cell-type mesh `native` and `cells`, `view`, `io` for `.fmsh`, `quality`, `fixtures`)
//! * meshing: `meshing` (geom models to meshes: `MeshPlan`, `geomdrive`, `mesher` 2D shapes to triangles, `mesher3d` and `refine3d` 3D solids to tetrahedra, structured PML `pml`, `pml3`, `pml_box`, `physics_size`, `validate`)
//! * adaptation: `adapt` (metric adaptation 2D and 3D, error estimators `estimate`, `estimate3`)
//! * finite-element basics: `fem` (quadrature, shape functions, scalar matrices, constraints), `fields` (gradients, node
//!   averaging, plotting surface), `material` (anisotropic data, rotations), `expr` (user expressions)
//! * linear algebra: `linalg` (the one library: `Csr<T>` matrices, direct and Krylov solvers, eigensolvers, multigrid;
//!   `linalg::driver` is the door the kernels use)
//! * physics: `dfs`, `sms` and `ems` are crates since phase 4 (ADR 41, 42 to 47); what they share stays here: `periodic` (Floquet pairs, used by the solid and the cavity
//!   kernels), `edge` (the 2D edge-element basis), `pml` (coordinate stretching), `props` (materials per domain), `section` (the wafer axes)
//! * optical mode sets: `optical` (modes with the data to evaluate their fields; the polarization analysis moved to `nyx-analysis` in phase 5)
//! * ported libraries: `linalg` (see above; Krylov, eigen, preconditioner, time-stepping and optimisation library for
//!   real and complex systems) and the mesher port in `mesh`, `meshing` and `adapt`, see docs/ports
//! * files: none: projects are JSON (engine), meshes are `.fmsh` (`mesh::io`); HDF5 export is a later engine export
//!
//! Lint note: `clippy::needless_range_loop` is allowed crate-wide. The numerical kernels (element matrices, orthogonalisation,
//! multigrid smoothers) index several arrays with one loop variable (`for i in 0..n { y[i] += a[i] * x[idx[i]] }`); the
//! indexed form is the one that matches the formulas, and rewriting 49 such loops as iterator chains would not change
//! the machine code. `too_many_arguments` is allowed on the few assembly/solver entry points that take a full problem
//! description (each marked where it applies).
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]
pub mod optical;
pub mod edge;
pub mod expr;
pub mod fieldexpr;
pub mod fem;
pub mod fields;
pub mod material;
pub mod adapt;
pub mod mesh;
pub mod meshing;
pub mod par;
pub mod periodic;
pub mod pml;
pub mod props;
pub mod progress;
pub mod section;
pub mod testkit;
pub mod linalg;
pub mod validate;

