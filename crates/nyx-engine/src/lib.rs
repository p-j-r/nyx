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
pub mod commands;
pub mod engine;
pub mod geomdesc;
pub mod geometry;
pub mod jobs;
pub mod legacy;
pub mod meshing;
pub mod project;
pub mod solve;
pub mod store_link;

pub use commands::{Access, CommandDesc};
pub use engine::Engine;
pub use jobs::{JobId, JobInfo, JobOutput, JobRunner, JobState};
pub use meshing::{FemcoreMesher, MeshedModel, Mesher};
pub use project::{MeshSettings, Project, SweepLayer};

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
