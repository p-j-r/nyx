//! Per-phase wall-clock timers, so every long operation reports where its time went.

use std::time::Instant;

/// Ordered list of `(phase, seconds)`; repeated phases accumulate.
#[derive(Clone, Debug, Default)]
pub struct Phases {
    pub entries: Vec<(String, f64)>,
}

impl Phases {
    /// Run `f` and record its duration under `name`.
    pub fn time<R>(&mut self, name: &str, f: impl FnOnce() -> R) -> R {
        let t = Instant::now();
        let r = f();
        self.add(name, t.elapsed().as_secs_f64());
        r
    }

    /// Add `secs` to phase `name`.
    pub fn add(&mut self, name: &str, secs: f64) {
        match self.entries.iter_mut().find(|(n, _)| n == name) {
            Some(e) => e.1 += secs,
            None => self.entries.push((name.to_string(), secs)),
        }
    }

    /// Total of all phases.
    pub fn total(&self) -> f64 {
        self.entries.iter().map(|e| e.1).sum()
    }

    /// One-line summary such as `eig 1.2 ms, scale 0.1 ms`.
    pub fn summary(&self) -> String {
        self.entries.iter().map(|(n, s)| format!("{n} {:.3} ms", s * 1e3)).collect::<Vec<_>>().join(", ")
    }
}
