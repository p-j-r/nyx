//! Scripts as a front-end of the nyx engine (productivity plan S4, ADR [89]).
//!
//! A script is a text file in an embedded, pure-Rust language. It calls engine commands (`Engine::call`, the same JSON as the window and the MCP
//! server), and what it prints or returns is all that goes back to the caller. The host gives it no file, network, clock or random access and
//! stops it when it runs too long or allocates too much.
//!
//! This is the S4 spike: both candidates (Starlark and Rhai) are bound to the same command catalog so the same sample runs in each.

use serde_json::Value;

pub mod rhai_host;
#[cfg(feature = "starlark")]
pub mod star;

/// What the host allows a script to use. The defaults are strict; a run may raise them (S4 decision).
#[derive(Clone, Debug)]
pub struct Limits {
    /// interpreter steps (Starlark ticks, Rhai operations)
    pub steps: u64,
    /// bytes of script-owned memory (Starlark heap; Rhai bounds strings, arrays and maps instead)
    pub heap_bytes: usize,
    /// nested calls
    pub depth: usize,
}

impl Default for Limits {
    fn default() -> Limits {
        Limits { steps: 50_000_000, heap_bytes: 256 << 20, depth: 50 }
    }
}

/// What a finished script gave back.
#[derive(Debug)]
pub struct Outcome {
    /// the script's value (`result` for Starlark, the last expression for Rhai), as JSON
    pub value: Value,
    /// every line the script printed
    pub printed: Vec<String>,
    /// engine commands the script sent
    pub calls: u64,
}

/// Fixed-point text of a number, for the tables a script prints. Neither language has a precision in its own formatting that suits it
/// (Starlark's `%` has none, Rhai has no printf), so the host gives both the same `fmt(x, digits)`.
pub fn fmt_fixed(x: f64, digits: i64) -> String {
    format!("{:.*}", digits.clamp(0, 17) as usize, x)
}
