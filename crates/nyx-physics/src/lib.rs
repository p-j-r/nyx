//! nyx-physics: the interface between nyx and a physics (ADR 25). A physics crate (`nyx-dfs`, `nyx-sms`, `nyx-ems`, ...) implements
//! `Physics`; the engine finds it in a `Registry`; the window, the MCP server and the mesher read only its descriptors.
//!
//! * `physics`: the `Physics` trait, `Ctx` (progress, cancel, imported fields) and `LabelResolver` (names to mesh labels)
//! * `spec`: the problem envelope (`ProblemSpec`, `StudySpec`, `BcSpec`, `MaterialAssign`) and `canonical_json` for exemplar files
//! * `desc`: descriptors of studies, conditions, properties, fields, couplings and the sizing hint
//! * `result`: `SolveResult`, `Solution`, `Field`, plain serialisable data
//! * `validate`, `schema`: structural validation against the descriptors and JSON Schemas generated from them
//! * `registry`: the physics a build contains; solve order of a one-way chain
//! * `error`: the error type, with `Cancelled`
//! * `contract`: helpers for the contract tests of every physics crate
//!
//! It depends on `femcore` for the mesh type and the progress handle, and on nothing else of nyx.

pub mod contract;
pub mod desc;
pub mod error;
pub mod physics;
pub mod registry;
pub mod result;
pub mod schema;
pub mod spec;
pub mod validate;

pub use desc::*;
pub use error::Error;
/// The mesh type physics take. Today `femcore::mesh::SimplexMesh`; restructure phase 2 makes it the one native mesh. Name it through this path.
pub use femcore::mesh::Mesh;
pub use physics::{Ctx, LabelResolver, Physics};
pub use registry::Registry;
pub use result::{CellField, Field, Scalar, Solution, SolveResult};
pub use spec::{canonical_json, BcSpec, MaterialAssign, ProblemSpec, StudySpec};
pub use validate::{check_params, validate_spec};
