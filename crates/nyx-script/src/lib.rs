//! Scripts as a front-end of the nyx engine (productivity plan S4, ADR [89]).
//!
//! A script is a text file in Python syntax (`*.nys`). Monty, a Python interpreter written in Rust, runs it in a *worker process*; each
//! `nyx.<command>(name=value, ...)` the script makes pauses it and is answered here by `Engine::call`, the same JSON the window and the MCP server
//! use. What a script prints or returns is all that goes back to the caller. The host gives the script no file, network, clock or random access
//! outside the engine commands, stops it when it computes too long, allocates too much or sends too many commands, and an engine path guard keeps
//! the file arguments of those commands inside one folder.
//!
//! * `kernel`: one Monty session that outlives a script (ADR 94): `run` a script on it, `wait` for one that is waiting for a job, `check` against
//!   its variables; a command that waits for a job is polled by the host so no call blocks for ever; `nyx.aio` futures; a dump replaces a dead worker
//! * `host`: `run` (execute one script against an engine) and `check` (parse and type-check it against the catalog without running anything),
//!   each on a kernel of its own
//! * `stubs`: the `nyx` module's type stub generated from `Engine::catalog()`: the static check and the file an editor reads
//! * `convert`: JSON to and from script values
//! * `trace`: a script function given as a formula is called once with a symbolic point and becomes formula text (N3, ADR 90); the prelude the
//!   worker is fed, its stub, and the static guard against a branch on a value of the point
//!
//! The worker program is `nyx-script-worker` (built from this package); it has to be found next to the running program, or named by the
//! environment variable `NYX_SCRIPT_WORKER`.

use std::path::PathBuf;

pub mod convert;
pub mod host;
pub mod kernel;
pub mod stubs;
pub mod trace;

/// What the host allows a script to use. The defaults are strict; a run may raise them (S4 decision).
#[derive(Clone, Debug)]
pub struct Limits {
    /// seconds the script itself may compute in all, not counting the time the engine takes to answer its commands (a 40 s solve is not script time)
    pub seconds: f64,
    /// bytes of memory the worker process may use; past it the worker is stopped and the engine goes on
    pub memory_bytes: usize,
    /// nested function calls
    pub depth: usize,
    /// engine commands the script may send
    pub calls: usize,
    /// bytes a script may print
    pub output_bytes: usize,
}

impl Default for Limits {
    fn default() -> Limits {
        Limits { seconds: 30.0, memory_bytes: 256 << 20, depth: 100, calls: 20_000, output_bytes: 1 << 20 }
    }
}

/// How a script is run.
#[derive(Clone, Debug)]
pub struct Options {
    /// the folder the file and folder arguments of the script's commands must lie in (the engine's path guard); relative paths are taken from it
    pub root: PathBuf,
    /// the worker program; `None` looks for `nyx-script-worker` next to the running program
    pub worker: Option<PathBuf>,
}

impl Options {
    pub fn new(root: impl Into<PathBuf>) -> Options {
        Options { root: root.into(), worker: None }
    }
}

/// What a finished script gave back.
#[derive(Debug)]
pub struct Outcome {
    /// the script's value (its last expression), as JSON
    pub value: serde_json::Value,
    /// every line the script printed
    pub printed: Vec<String>,
    /// engine commands the script sent
    pub calls: u64,
}

/// The script worker program: named by `NYX_SCRIPT_WORKER`, else `nyx-script-worker` in the folder of the running program or the one above it
/// (the test programs of cargo sit in `deps/`).
pub fn find_worker() -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("NYX_SCRIPT_WORKER") {
        return Some(PathBuf::from(p));
    }
    let name = format!("nyx-script-worker{}", std::env::consts::EXE_SUFFIX);
    let exe = std::env::current_exe().ok()?;
    let dir = exe.parent()?;
    [dir.to_path_buf(), dir.parent()?.to_path_buf()].into_iter().map(|d| d.join(&name)).find(|p| p.is_file())
}
