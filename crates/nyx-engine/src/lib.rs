//! nyx-engine: the one place that knows how a nyx session works. The window (`app`) and the MCP server (`nyx-mcp`) are thin front-ends
//! over it: every action of either is a command of the engine (restructure phase 7, ADR 23 and 25).
//!
//! SKELETON (restructure foundation). Today it holds only the registry of compiled-in physics, built from cargo features. Phase 7
//! adds `Project` (moved from `app/state.rs`), the command set, the job runner (progress, cancel), persistence through `nyx-store`, the
//! batch runner and the exports. The rule it enforces from now on: **nothing outside the physics crates names a physics**; code that
//! needs one asks `registry().get(id)`.

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
