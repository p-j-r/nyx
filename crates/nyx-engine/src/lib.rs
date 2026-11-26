// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-engine: the one place that knows how a nyx session works. The window (`app`) and the MCP server (`nyx-mcp`) are thin front-ends
//! over it: every action of either is a command of the engine (restructure phase 7, ADR 23 and 25).
//!
//! * `project`: `Project` (feature list, mesh settings, one problem per physics), saved as JSON
//! * `geometry`: editing a feature list; the operations come from `geomdesc`
//! * `meshing`: `MeshedModel` (mesh and names, the label resolver) and the `Mesher` seam that restructure phase 2 rewrites
//! * `jobs`: the job runner (progress, cancel); `solve`: the coupled chain of physics; `store_link`: persistence through `nyx-store`
//! * `engine`: `Engine`, the state around one project; `commands`: every action by name, JSON in and out, with schemas (`catalog`)
//!
//! Built before phase 4 (the kernels): `solve` goes through `Physics::solve`, which answers `NotYetMoved` for the skeleton physics until
//! phase 4 fills them, and meshing goes through one seam. The rule it enforces: **nothing outside the physics crates names a physics**;
//! code that needs one asks `registry().get(id)`.

/// The analyses of results (restructure phase 5); the window draws the polarization from here until it calls engine commands.
pub use nyx_analysis as analysis;
pub mod amr;
pub mod analyse;
pub mod commands;
pub mod convergence;
pub mod engine;
pub mod estimate;
mod figure;
pub mod gds;
pub mod geomdesc;
pub mod geometry;
pub mod guide;
#[cfg(feature = "hdf5")]
pub mod h5io;
pub mod jobs;
pub mod material_fetch;
pub mod material_warnings;
pub mod materials;
#[cfg(feature = "hdf5")]
pub mod xdmf;

pub mod apiref;
pub mod commandlist;
pub mod meshing;
pub mod path_guard;
pub mod pml;
pub mod preview;
pub mod process;
pub mod project;
pub mod qmc;
pub mod queue;
pub mod reduce;
pub mod results_file;
pub mod reuse;
pub mod scriptout;
pub mod shared_queue;
pub mod solve;
pub mod store_link;
pub mod sweep;
pub mod validated;
pub mod worker;

pub use commands::{Access, CommandDesc};
pub use engine::Engine;
pub use estimate::{Estimate, History};
pub use jobs::{JobId, JobInfo, JobOutput, JobRunner, JobState, Lane, Submit, WorkerRun};
pub use meshing::{FemcoreMesher, MeshedModel, Mesher};
pub use project::{MeshSettings, PmlModeRequest, Project, SweepLayer};
pub use queue::QueueOpts;

use nyx_physics::{Error, Registry};
#[allow(unused_imports)]
use std::sync::Arc;

/// The physics compiled into this build: `dfs`, `sms`, `ems`, each behind the cargo feature of the same name.
pub fn registry() -> Registry {
    #[allow(unused_mut)]
    let mut r = Registry::new();
    #[cfg(feature = "dfs")]
    r.register(Arc::new(nyx_dfs::Dfs)).expect("dfs registers once");
    #[cfg(feature = "sms")]
    r.register(Arc::new(nyx_sms::Sms)).expect("sms registers once");
    #[cfg(feature = "ems")]
    r.register(Arc::new(nyx_ems::Ems)).expect("ems registers once");
    r
}

/// Looks up the physics a problem names; the error text lists what this build has.
pub fn physics_for<'a>(registry: &'a Registry, id: &str) -> Result<&'a dyn nyx_physics::Physics, Error> {
    registry.get(id).ok_or_else(|| {
        Error::Unsupported(format!("this build has no physics '{id}' (it has: {})", registry.ids().join(", ")))
    })
}

/// What a front-end needs from the layers below without depending on them: it draws on these types and edits the geometry with them, so
/// the window and the MCP server name only `nyx-engine` (restructure phase 7; the layering table has no other edge for them).
pub mod facade {
    /// The mesh as the plots read it (vertices, edges, P2 dofs), the surface that is drawn and the nodal fields drawn on it.
    pub use femcore::fields::Display;
    pub use femcore::mesh::view::ViewMesh;
    pub use femcore::mesh::SimplexMesh;
    /// Evaluation of a feature list (what the geometry editor previews).
    pub use femcore::meshing::geomdrive;
    pub use femcore::meshing::request::PmlModeRequest;
    pub use femcore::progress::Progress;
    /// The geometry kernel: the feature list the geometry editor changes, evaluation, selectors.
    pub use geom;
    /// The problem envelope, descriptors and results every physics speaks.
    pub use nyx_physics as physics;
    /// Workspaces, workflows, experiments, runs, versions, search and `.nyx` files.
    pub use nyx_store as store;
}
