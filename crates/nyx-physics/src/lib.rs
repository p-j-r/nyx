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
/// The mesh a physics is given: the one native cell mesh, 2D or 3D (restructure phase 4, ADR 32 and 41; it was `SimplexMesh`).
pub use femcore::mesh::Mesh;
/// The numbering of unknowns (edges, P2 dofs) a kernel assembles on. Each solve builds its own from the `Mesh` with `discretise`. The type
/// is still the solver mesh of femcore, which the window uses until phase 7. Kernels name only `Mesh` and `Discretisation`.
pub type Discretisation = femcore::mesh::SimplexMesh;

/// The `Discretisation` a physics assembles on, made from the `Mesh` it was given (non-simplex cells are cut into simplices).
pub fn discretise(mesh: &Mesh) -> Result<Discretisation, Error> {
    mesh.discretise().map_err(Error::invalid)
}

/// Refuses a mesh with a NaN or infinite coordinate, naming the vertex and the axis.
pub fn check_mesh(mesh: &Mesh) -> Result<(), Error> {
    match mesh.first_non_finite() {
        Some((i, k)) => Err(Error::NonFinite(format!(
            "the mesh has a non-finite coordinate (point {i}, axis {k}): NaN or infinity in the geometry"
        ))),
        None => Ok(()),
    }
}
pub use physics::{Ctx, LabelResolver, Physics};
pub use registry::Registry;
pub use result::{cell_fields_back, cell_fields_of, CellField, Field, Scalar, Solution, SolveResult};
pub use spec::{canonical_json, BcSpec, MaterialAssign, ProblemSpec, StudySpec};
pub use validate::{check_params, check_result_finite, validate_spec};
