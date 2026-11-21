// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Running and checking one script (S4): `nyx run FILE.nys` and `nyx check FILE.nys`. Each is a script on a kernel of its own (`kernel::Kernel`,
//! S7), which starts one worker process, feeds it the script and answers each pause until the script completes or a limit ends it. The engine
//! never leaves this process.
//!
//! What a script can do, and what it gets when it tries more:
//!
//! * engine commands: `nyx.<command>(...)`, answered by `Engine::call`; an engine error is raised in the script as `RuntimeError`, which it can catch
//! * files: none. `open`, `pathlib` and the `os` calls raise `PermissionError`; the commands that read or write files (`project_save`, ...) work
//!   only inside `Options::root`, enforced by the engine's path guard
//! * clock and random numbers: the clock is refused (`RuntimeError`), `random` starts from the fixed seed 0, `sleep` returns at once, so the same
//!   script and inputs give the same numbers (ADR 11)
//! * time, memory, depth, commands, output: `Limits`; the worker process is stopped past them and the engine goes on

use nyx_engine::Engine;

use crate::kernel::{Kernel, Step};
use crate::{Limits, Options, Outcome};

/// Runs a script. Its value is its last expression. The engine's path guard is on for the run and restored after it. A person runs their own
/// file, so the commands that discard the project are allowed.
pub fn run(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<Outcome, String> {
    let src = src.strip_prefix('\u{feff}').unwrap_or(src);
    let mut kernel = Kernel::new(engine, limits.clone(), opts.clone())?;
    kernel.set_script_name(name);
    match kernel.run(engine, name, src, None, None, true)? {
        Step::Done(o) => Ok(o),
        Step::Running(_) => Err("internal error: a script with no time limit is still waiting".to_string()),
    }
}

/// Parses and type-checks a script against the commands of `engine` and runs nothing: a misspelled command or argument, a missing argument, a
/// syntax error. `Ok` when the script is fine.
pub fn check(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<(), String> {
    let src = src.strip_prefix('\u{feff}').unwrap_or(src);
    let mut kernel = Kernel::new(engine, limits.clone(), opts.clone())?;
    kernel.set_script_name(name);
    kernel.check(engine, name, src)
}
