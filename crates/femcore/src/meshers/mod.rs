//! # meshers
//!
//! Metric-based mesh adaptation, structured PML meshes, mixed-cell meshes and a compact lossless mesh/results format,
//! ported from the staging crate `femmesh` of the inspired-rust-port repository (BSD-2-Clause; provenance in
//! `docs/ports/meshers/PROVENANCE.md`, takeover notes in `docs/ports/meshers/HANDOFF.md`).
//!
//! * `adapt`: `metric` (SPD metric algebra, Hessian to metric), `field` (metric fields), `topo`, `geom`, `ops`, `drive`
//!   (2D split/collapse/flip/smooth adaptation loop), `pml` (structured PML meshes around a box, 2D quadrilaterals and 3D
//!   prisms/hexahedra)
//! * `io`: the `.fmsh` mesh and result file (CRC-32, entropy coded, lossless)
//! * `validate`: an independent scalar Helmholtz solver with PML used to test the meshes (not a production solver)
//! * `mesh`, `cellmesh`: the mesher's own dimension-generic mesh types; `interop` translates them from and to
//!   `crate::mesh::SimplexMesh` and `crate::cellmesh::CellMesh` and turns the PML regions into `problem::Pml` domains
//! * `fixtures`, `quality`, `par`, `timer`: test meshes, quality statistics, chunked reductions, phase timers
//!
//! Limits (from the author's notes): the adaptation loop is 2D and serial; the PML mesh is for boxes; nothing was run against
//! the nyx physics solvers before this merge (the new tests in `interop` check areas, volumes and labels only).

// Index loops over the spatial dimension are the clearest form for small fixed-size matrix code, and
// `!(x > 0.0)` is deliberate: it also rejects NaN.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

pub mod adapt;
pub mod cellmesh;
pub mod fixtures;
pub mod interop;
pub mod io;
pub mod mesh;
pub mod par;
pub mod quality;
pub mod timer;
pub mod validate;
