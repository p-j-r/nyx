// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Warm starts between the points of a sweep (productivity plan P5b, ADR 96): the engine's side of `femcore::warm`. A sweep asked with
//! `reuse: "warm"` links its points in a chain: point i starts its eigensolves from the modes point i-1 found, when point i-1 has finished and
//! left them, and from a random block (as every solve does) when it has not. Nothing else ever carries from one solve to another: a solve
//! that is not a point of a warm sweep has no link and takes the path it always took, bit for bit.
//!
//! A link is a **source** (where the previous modes come from), a **cell** (where this point leaves its own) and the point's index. The cell
//! is filled when the point's job ends, whatever way it ends (`CellGuard`), so a point that waits for its predecessor is never left waiting
//! for one that failed, was cancelled or never ran. A point waits only in the `shared` lane (ADR 95), where it holds no core; in the heavy lane
//! the lane runs the points in order, so a predecessor that has not finished is one that cannot (held, cancelled): the point is cold, and the
//! record says why (a wait there would hold the lane for ever).
//!
//! The vectors go to a worker process, and to the queue worker, as bytes in the job's folder (`warm_in.bin`, and back `warm_out.bin` with
//! `warm_calls.json`), so the three ways a solve can run give one answer.

use crate::jobs::{JobId, JobRunner};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// What a finished point leaves for the next: the modes as bytes (none when it failed, or its eigensolves left none) and what its eigensolves did.
#[derive(Clone, Debug, Default)]
pub struct WarmDone {
    pub bytes: Option<Arc<Vec<u8>>>,
    pub calls: Vec<Value>,
}

/// Where a point leaves its modes. Empty until the point's job has ended.
#[derive(Default)]
pub struct WarmCell(Mutex<Option<WarmDone>>);

impl WarmCell {
    pub fn new() -> Arc<WarmCell> {
        Arc::new(WarmCell::default())
    }

    pub fn get(&self) -> Option<WarmDone> {
        self.0.lock().expect("warm cell").clone()
    }

    pub fn set(&self, d: WarmDone) {
        *self.0.lock().expect("warm cell") = Some(d);
    }

    fn is_set(&self) -> bool {
        self.0.lock().expect("warm cell").is_some()
    }
}

/// Where a point's start block comes from.
#[derive(Clone)]
pub enum WarmSource {
    /// the first point of a sweep: nothing to start from, the modes are only recorded
    Cold,
    /// modes handed over directly (a worker process, which is given the previous point's bytes)
    Bytes(Arc<Vec<u8>>),
    /// the cell of point `index`, solved by job `job`
    Cell { cell: Arc<WarmCell>, job: JobId, index: usize },
}

/// A point's place in a warm sweep.
#[derive(Clone)]
pub struct WarmLink {
    pub source: WarmSource,
    pub mine: Arc<WarmCell>,
    pub index: usize,
}

impl std::fmt::Debug for WarmLink {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "WarmLink(point {})", self.index)
    }
}

/// What a point found to start from.
pub struct Resolved {
    pub bytes: Option<Arc<Vec<u8>>>,
    /// the point whose modes they are
    pub from: Option<usize>,
    /// why the point is cold, when it is
    pub why_cold: Option<&'static str>,
}

/// Sets the point's cell empty-handed if its job ends without having filled it.
pub struct CellGuard(Arc<WarmCell>);

impl Drop for CellGuard {
    fn drop(&mut self) {
        if !self.0.is_set() {
            self.0.set(WarmDone::default());
        }
    }
}

impl WarmLink {
    pub fn guard(&self) -> CellGuard {
        CellGuard(self.mine.clone())
    }

    /// A link for a worker process: it is given the previous point's modes (or none) and leaves its own in `mine`.
    pub fn for_worker(input: Option<Arc<Vec<u8>>>, mine: Arc<WarmCell>) -> WarmLink {
        WarmLink { source: input.map_or(WarmSource::Cold, WarmSource::Bytes), mine, index: 0 }
    }

    /// Finds what to start from. With `wait` a point whose predecessor is still running waits for it (and stops at once on its own cancel);
    /// without, or when the predecessor can no longer finish, it is cold.
    pub fn resolve(
        &self,
        jobs: &JobRunner,
        wait: bool,
        progress: &Arc<femcore::progress::Progress>,
    ) -> Result<Resolved, Error> {
        match &self.source {
            WarmSource::Cold => {
                Ok(Resolved { bytes: None, from: None, why_cold: Some("the first point of the sweep") })
            }
            WarmSource::Bytes(b) => Ok(Resolved { bytes: Some(b.clone()), from: None, why_cold: None }),
            WarmSource::Cell { cell, job, index } => {
                let _installed = femcore::progress::install(progress);
                let mut _stage = None;
                loop {
                    if let Some(d) = cell.get() {
                        return Ok(match d.bytes {
                            Some(b) => Resolved { bytes: Some(b), from: Some(*index), why_cold: None },
                            None => Resolved {
                                bytes: None,
                                from: Some(*index),
                                why_cold: Some("the previous point left no modes (it failed or was cancelled)"),
                            },
                        });
                    }
                    let open = jobs.info(*job).is_some_and(|i| i.state.is_open());
                    if !open {
                        // the job ended after setting the cell, or never ran: look once more, then give up
                        return Ok(match cell.get().and_then(|d| d.bytes) {
                            Some(b) => Resolved { bytes: Some(b), from: Some(*index), why_cold: None },
                            None => Resolved {
                                bytes: None,
                                from: Some(*index),
                                why_cold: Some("the previous point did not run"),
                            },
                        });
                    }
                    if !wait {
                        return Ok(Resolved {
                            bytes: None,
                            from: Some(*index),
                            why_cold: Some("the previous point had not finished when this one started"),
                        });
                    }
                    if progress.is_cancelled() {
                        return Err(Error::Cancelled);
                    }
                    if _stage.is_none() {
                        _stage = Some(femcore::progress::stage(
                            "waiting for the previous point of the sweep (its modes are this point's start)",
                        ));
                    }
                    std::thread::sleep(Duration::from_millis(50));
                }
            }
        }
    }
}

/// The record of what a point did with its start, for the table row and the run: whether it started warm and from which point, how many
/// restarts each eigensolve took, and that a warm result is not confirmed until it is solved again cold.
pub fn reuse_json(link: &WarmLink, r: &Resolved, calls: &[Value]) -> Value {
    let applied = calls.iter().any(|c| c["warm"] == json!(true));
    let mut v = json!({
        "mode": "warm", "point": link.index, "warm": applied, "from_point": if applied { json!(r.from) } else { Value::Null },
        "eigensolves": calls,
    });
    if applied {
        v["exploratory"] = json!(true);
        v["note"] = json!("started from the modes of the previous point; the same numbers to the solver's tolerance, not the same bits as a cold solve. Solve this point cold to confirm it");
    } else {
        v["why_cold"] =
            json!(r.why_cold.unwrap_or("the eigensolves found nothing of the same kind and size to start from"));
    }
    v
}

/// One eigensolve's record as JSON.
pub fn call_json(c: &femcore::warm::Call) -> Value {
    json!({ "kind": c.kind, "unknowns": c.n, "modes": c.nev, "warm": c.warm, "iterations": c.iterations })
}

/// Runs `f` (the in-process solve of one point) under a warm context made from `resolved`, and returns what it returned with the modes the
/// solve left and what its eigensolves did.
pub fn with_context<R>(resolved: &Resolved, f: impl FnOnce() -> R) -> (R, WarmDone) {
    let entries = resolved.bytes.as_ref().and_then(|b| femcore::warm::from_bytes(b).ok()).unwrap_or_default();
    let w = femcore::warm::Warm::new(entries);
    let out = {
        let _g = femcore::warm::install(&w);
        f()
    };
    let left = w.output();
    let bytes = (left.iter().any(Option::is_some)).then(|| Arc::new(femcore::warm::to_bytes(&left)));
    (out, WarmDone { bytes, calls: w.calls().iter().map(call_json).collect() })
}
