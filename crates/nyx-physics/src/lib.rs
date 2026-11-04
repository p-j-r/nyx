//! nyx-physics: the interface between nyx and a physics (ADR 25). A physics crate (`nyx-dfs`, `nyx-sms`, `nyx-ems`, ...) implements
//! `Physics`; the engine finds it in a `Registry`; the window, the MCP server and the mesher read only its descriptors.
//!
//! * `physics`: the `Physics` trait, `Ctx` (progress, cancel, imported fields) and `LabelResolver` (names to mesh labels)
//! * `spec`: the problem envelope (`ProblemSpec`, `StudySpec`, `BcSpec`, `MaterialAssign`) and `canonical_json` for exemplar files
//! * `desc`: descriptors of studies, conditions, properties, fields, couplings and the sizing hint
//! * `result`: `SolveResult`, `Solution`, `Field`, plain serialisable data
//! * `validate`, `schema`: structural validation against the descriptors and JSON Schemas generated from them
//! * `decode`: typed parameters with the descriptor defaults, names and selectors resolved to mesh numbers (what every physics does first)
//! * `registry`: the physics a build contains; solve order of a one-way chain
//! * `error`: the error type, with `Cancelled`
//! * `contract`: helpers for the contract tests of every physics crate
//!
//! It depends on `femcore` for the mesh type and the progress handle, and on nothing else of nyx.

pub mod contract;
pub mod decode;
pub mod desc;
pub mod error;
pub mod physics;
pub mod registry;
pub mod result;
pub mod schema;
pub mod spec;
pub mod testkit;
pub mod validate;

pub use desc::*;
pub use error::Error;
/// The mesh type physics take. Today `femcore::mesh::SimplexMesh`; restructure phase 2 makes it the one native mesh. Name it through this path.
pub use femcore::mesh::Mesh;
/// The mesh with the numbering of unknowns (edges, P2 dofs) a kernel assembles on. Today the same type as `Mesh`; the removal of
/// `SimplexMesh` (ADR 32, ADR 41) makes `Mesh` the native cell mesh and gives each solve its own `Discretisation` built by `discretise`.
/// Kernels name only these two paths, never `femcore::mesh::SimplexMesh`.
pub type Discretisation = femcore::mesh::SimplexMesh;

/// The `Discretisation` a physics assembles on, made from the `Mesh` it was given.
pub fn discretise(mesh: &Mesh) -> Result<std::borrow::Cow<'_, Discretisation>, Error> {
    Ok(std::borrow::Cow::Borrowed(mesh))
}
pub use decode::{materials_by_label, resolve_conditions, Params, ResolvedBc};
pub use physics::{Ctx, LabelResolver, Physics};
pub use registry::Registry;
pub use result::{cell_fields_back, cell_fields_of, CellField, Field, Scalar, Solution, SolveResult};
pub use spec::{canonical_json, BcSpec, MaterialAssign, ProblemSpec, StudySpec};
pub use validate::{check_params, check_result_finite, validate_spec};
