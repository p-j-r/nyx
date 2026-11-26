// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-store: where nyx keeps its work.
//!
//! A *workspace* folder holds *workflows*; a workflow holds *experiments* (simulation definitions, forks and links of each other);
//! an experiment has an append-only table of *runs*. This crate owns all of it: the manifest and run schemas (`model`, `runs`), the
//! folder (`flow`), the version history (`versions`), the auto-generated concept card (`card`), the dashboard rules and settings
//! (`workspace`), search (`search`) and `.nyx` files (`archive`). It has no window code, so the app, `nyx --batch`, scripts and a
//! future assistant all call the same functions.

pub mod archive;
pub mod card;
pub mod flow;
pub mod model;
pub mod runs;
pub mod search;
pub mod util;
pub mod versions;
pub mod workspace;

pub use flow::WorkflowDir;
pub use model::*;
pub use runs::{BuildProvenance, RunDiagnosticsSummary, RuntimeProvenance};
pub use workspace::{Entry, Sections, Settings, Workspace};
