//! What a running solve is doing, and a way to stop it.
//!
//! A `Progress` is shared between the thread that solves and whoever watches it (the window, a script, a
//! test). It holds three things:
//!
//! * a **cancel flag**: the watcher calls `cancel()`; long loops in the solvers call `check()` once per
//!   iteration (assembly stages, each Krylov/LOBPCG/CG iteration, each fixed-point step, each Brillouin
//!   pump/Stokes pair) and return the error `"cancelled"` when the flag is set. A sparse factorisation
//!   cannot be interrupted half way, so a cancel during one takes effect when it finishes;
//! * the **current stage** (`"factorising the complex matrix: 161,625 unknowns"`) with its elapsed time and an
//!   optional fraction done, for the progress text in the window;
//! * the **stage log**: one plain-language line per finished stage with its duration, plus notes such as
//!   `"iteration 30: residual 3.1e-6"`. `solve::solve` appends it to `SolveResult::log`, so the window shows
//!   it and the autosaved run keeps it. With the environment variable `NYX_PROFILE` set every line is also
//!   echoed to stderr as it is written.
//!
//! The solvers do not take a `Progress` argument. `solve::solve_with` installs one for the current thread
//! (`install`) and the free functions below (`stage`, `note`, `check`, `fraction`) use it; without an installed
//! one they do nothing, so library code and tests that call a solver directly are unaffected. Worker threads
//! spawned inside a solver do not see it; only the solving thread logs and checks the flag.

use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// The error text of a cancelled solve.
pub const CANCELLED: &str = "cancelled";

/// Longest stage log kept (further lines are dropped after one "log truncated" line).
const MAX_LOG_LINES: usize = 600;

struct Inner {
    /// stack of the stages currently running (outermost first): text and start time
    stack: Vec<(String, Instant)>,
    fraction: Option<f64>,
    log: Vec<String>,
    truncated: bool,
}

/// Shared state of one job. Create with `Progress::new()`.
pub struct Progress {
    cancel: AtomicBool,
    start: Instant,
    echo: bool,
    inner: Mutex<Inner>,
}

/// A snapshot for display.
#[derive(Clone, Debug, Default)]
pub struct Status {
    /// innermost running stage, empty before the first one
    pub stage: String,
    /// seconds in that stage
    pub stage_seconds: f64,
    /// seconds since the job started
    pub total_seconds: f64,
    pub fraction: Option<f64>,
    pub cancelled: bool,
}

impl Progress {
    pub fn new() -> Arc<Progress> {
        Arc::new(Progress { cancel: AtomicBool::new(false), start: Instant::now(), echo: std::env::var_os("NYX_PROFILE").is_some(), inner: Mutex::new(Inner { stack: vec![], fraction: None, log: vec![], truncated: false }) })
    }

    /// Asks the solver to stop at its next check.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::SeqCst);
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancel.load(Ordering::SeqCst)
    }

    pub fn status(&self) -> Status {
        let g = self.inner.lock().unwrap();
        let (stage, stage_seconds) = g.stack.last().map(|(s, t)| (s.clone(), t.elapsed().as_secs_f64())).unwrap_or_default();
        Status { stage, stage_seconds, total_seconds: self.start.elapsed().as_secs_f64(), fraction: g.fraction, cancelled: self.is_cancelled() }
    }

    /// The finished-stage lines and notes so far.
    pub fn log(&self) -> Vec<String> {
        self.inner.lock().unwrap().log.clone()
    }

    /// Adds a line to the log; returns its index (None when the log is full).
    fn push(&self, depth: usize, text: String) -> Option<usize> {
        let mut g = self.inner.lock().unwrap();
        let line = format!("{}{}", "  ".repeat(depth), text);
        if self.echo {
            eprintln!("  [profile] {line}");
        }
        if g.log.len() >= MAX_LOG_LINES {
            if !g.truncated {
                g.truncated = true;
                g.log.push("... (stage log truncated)".into());
            }
            None
        } else {
            g.log.push(line);
            Some(g.log.len() - 1)
        }
    }

    /// Replaces a line written earlier (a stage line gets its duration when the stage ends).
    fn replace(&self, index: usize, depth: usize, text: String) {
        let mut g = self.inner.lock().unwrap();
        let line = format!("{}{}", "  ".repeat(depth), text);
        if self.echo {
            eprintln!("  [profile] {line}");
        }
        if let Some(l) = g.log.get_mut(index) {
            *l = line;
        }
    }
}

thread_local! {
    static CURRENT: RefCell<Option<Arc<Progress>>> = const { RefCell::new(None) };
}

/// Makes `p` the progress of the current thread until the returned guard is dropped.
pub fn install(p: &Arc<Progress>) -> InstallGuard {
    let prev = CURRENT.with(|c| c.replace(Some(p.clone())));
    InstallGuard { prev }
}

pub struct InstallGuard {
    prev: Option<Arc<Progress>>,
}

impl Drop for InstallGuard {
    fn drop(&mut self) {
        let prev = self.prev.take();
        CURRENT.with(|c| *c.borrow_mut() = prev);
    }
}

/// The progress installed on this thread, if any.
pub fn current() -> Option<Arc<Progress>> {
    CURRENT.with(|c| c.borrow().clone())
}

fn profile_env() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("NYX_PROFILE").is_some())
}

/// A running stage: its line stands in the log from the start ("... running") and gets its duration when the guard is
/// dropped, so the log reads in the order things happened with sub-stages indented below their stage.
pub struct Stage {
    p: Option<Arc<Progress>>,
    line: Option<usize>,
    /// without an installed `Progress` the stage is only echoed to stderr when `NYX_PROFILE` is set
    local: Option<(String, Instant)>,
}

/// Starts a stage with a plain-language description, e.g. `stage(format!("assembling the stiffness matrix: {n} unknowns"))`.
/// Keep the guard alive for the duration of the stage.
pub fn stage(text: impl Into<String>) -> Stage {
    let Some(p) = current() else {
        return Stage { p: None, line: None, local: if profile_env() { Some((text.into(), Instant::now())) } else { None } };
    };
    let text: String = text.into();
    let depth = {
        let mut g = p.inner.lock().unwrap();
        let depth = g.stack.len();
        g.stack.push((text.clone(), Instant::now()));
        g.fraction = None;
        depth
    };
    let line = p.push(depth, format!("{text} ..."));
    Stage { p: Some(p), line, local: None }
}

impl Drop for Stage {
    fn drop(&mut self) {
        if let Some((text, t)) = self.local.take() {
            eprintln!("  [profile] {text}: {:.2} s", t.elapsed().as_secs_f64());
        }
        let Some(p) = self.p.take() else { return };
        let (text, secs, depth) = {
            let mut g = p.inner.lock().unwrap();
            let Some((text, t)) = g.stack.pop() else { return };
            g.fraction = None;
            (text, t.elapsed().as_secs_f64(), g.stack.len())
        };
        let suffix = if p.is_cancelled() { " (cancelled)" } else { "" };
        let done = format!("{text}: {secs:.2} s{suffix}");
        match self.line {
            Some(i) => p.replace(i, depth, done),
            None => {
                p.push(depth, done);
            }
        }
    }
}

/// Adds a line to the stage log (and to the display), e.g. `"iteration 30: residual 3.1e-6"`.
pub fn note(text: impl Into<String>) {
    if let Some(p) = current() {
        let depth = p.inner.lock().unwrap().stack.len();
        p.push(depth, text.into());
    } else if profile_env() {
        eprintln!("  [profile] {}", text.into());
    }
}

/// Timing of a finished piece of work: `what` (the label) took `seconds`. Labels that start with four spaces are
/// fine-grained details (single multigrid set-up steps, mesher phases): they are only echoed to stderr under
/// `NYX_PROFILE` and kept out of the log shown to the user.
pub fn lap(what: &str, seconds: f64) {
    if what.starts_with("    ") {
        if profile_env() {
            eprintln!("  [profile] {what}: {seconds:.3} s");
        }
    } else {
        note(format!("{}: {seconds:.3} s", what.trim()));
    }
}

/// Fraction (0..1) of the current stage that is done, for a progress bar.
pub fn fraction(f: f64) {
    if let Some(p) = current() {
        p.inner.lock().unwrap().fraction = Some(f.clamp(0.0, 1.0));
    }
}

/// `Err("cancelled")` when the watcher asked to stop; call it once per iteration of a long loop.
pub fn check() -> Result<(), String> {
    match current() {
        Some(p) if p.is_cancelled() => Err(CANCELLED.to_string()),
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stages_are_logged_with_nesting_and_cancel_is_seen() {
        let p = Progress::new();
        let _g = install(&p);
        {
            let _outer = stage("outer");
            {
                let _inner = stage("inner");
                note("a note");
                assert_eq!(p.status().stage, "inner");
            }
            assert_eq!(p.status().stage, "outer");
        }
        let log = p.log();
        // in the order things happened: the stage, its sub-stage, the note inside it
        assert!(log[0].starts_with("outer: ") && log[0].contains(" s"), "{log:?}");
        assert!(log[1].starts_with("  inner: "), "{log:?}");
        assert!(log[2].starts_with("    a note"), "{log:?}");
        assert!(check().is_ok());
        p.cancel();
        assert_eq!(check().unwrap_err(), CANCELLED);
    }

    #[test]
    fn without_an_installed_progress_everything_is_a_no_op() {
        let _s = stage("nothing");
        note("nothing");
        assert!(check().is_ok());
    }
}
