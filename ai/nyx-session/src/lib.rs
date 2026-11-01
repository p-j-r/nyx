//! nyx session layer: the typed tools an AI host can call to build a geometry, mesh it, define and solve an optical mode problem and read
//! the result, with handles, an assumptions ledger and structured errors. It calls `femcore` and `geom` as they are and knows nothing about MCP.
pub mod guide;
pub mod params;
pub mod tools;
pub mod workspace;

pub use params::{Stop, ToolDef, ToolError};
pub use tools::{call, render, Outcome, TOOLS};
pub use workspace::Workspace;
