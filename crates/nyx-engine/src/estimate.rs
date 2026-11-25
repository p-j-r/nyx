// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What a solve will cost, learned from what solves cost (productivity plan Q2, ADR 86). Every solve that runs in a worker process measures its
//! own peak memory exactly (`process::peak_memory`); the engine keeps one row per such run (the physics of the chain, the dimension, the cells,
//! the unknowns, the peak memory, the seconds) and, before a solve starts, fits the rows of the same chain and dimension with a power law
//! (`value = a * cells^p`, least squares on the logarithms). It never invents a number: **with fewer than two runs of different sizes there is
//! no estimate**, and one beyond the sizes seen says it is extrapolated. The memory estimate gates the start (`Engine` waits for free memory, or
//! refuses what the machine cannot hold unless `force`); the time estimate is only shown. P0 measured memory growing as the 1.5th power of the
//! unknowns and time as the 2nd on the 3D cavity, so a constant per unknown would be wrong, and a fit is only as good as its rows: the answer
//! carries the number of runs and the largest error the fit makes on them, and admission uses the estimate plus that error.

use serde::{Deserialize, Serialize};
use serde_json::json;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// One finished worker run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Row {
    /// the chain that was solved: its physics codes sorted and joined by commas (`dfs,sms`)
    pub chain: String,
    pub dim: u8,
    pub cells: u64,
    pub dofs: u64,
    /// peak working set (bytes in RAM) and peak commit of the worker
    pub peak_working_set: u64,
    pub peak_commit: u64,
    pub seconds: f64,
}

/// An estimate with its basis.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Estimate {
    /// the memory the solve is expected to need, bytes of RAM
    pub bytes: u64,
    /// `bytes` plus the largest error the fit made on the runs it was made from (at least 10 %): what admission compares with free memory
    pub bytes_with_margin: u64,
    /// expected seconds (display only), when the runs have them
    pub seconds: Option<f64>,
    /// how many runs the fit used
    pub runs: usize,
    /// the cells asked for are outside the sizes of those runs
    pub extrapolated: bool,
    /// the exponent of the fit of memory against cells
    pub exponent: f64,
}

/// A power law `ln y = a + p ln x` fitted to `(x, y)`; `None` for fewer than two different x.
fn fit(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    let n = points.len() as f64;
    let (sx, sy) = points.iter().fold((0.0, 0.0), |(a, b), (x, y)| (a + x.ln(), b + y.ln()));
    let (mx, my) = (sx / n, sy / n);
    let sxx: f64 = points.iter().map(|(x, _)| (x.ln() - mx).powi(2)).sum();
    if sxx < 1e-6 {
        return None;
    }
    let sxy: f64 = points.iter().map(|(x, y)| (x.ln() - mx) * (y.ln() - my)).sum();
    let p = sxy / sxx;
    Some((my - p * mx, p))
}

/// The estimate for a solve of `chain` in `dim` dimensions on `cells` cells from `rows`, or `None` when the rows do not support one.
pub fn estimate(rows: &[Row], chain: &str, dim: u8, cells: u64) -> Option<Estimate> {
    let same: Vec<&Row> =
        rows.iter().filter(|r| r.chain == chain && r.dim == dim && r.cells > 0 && r.peak_working_set > 0).collect();
    let mem: Vec<(f64, f64)> = same.iter().map(|r| (r.cells as f64, r.peak_working_set as f64)).collect();
    let (_, p) = fit(&mem)?;
    // a fit through noise can have a silly exponent; memory of a direct solver does not fall with size or grow faster than the cube
    let p = p.clamp(0.5, 3.0);
    // with the exponent clamped the intercept is the mean of ln y - p ln x
    let a = mem.iter().map(|(x, y)| y.ln() - p * x.ln()).sum::<f64>() / mem.len() as f64;
    let predict = |x: f64| (a + p * x.ln()).exp();
    let worst = mem.iter().map(|(x, y)| (predict(*x) / y - 1.0).abs()).fold(0.0, f64::max);
    let x = cells as f64;
    let bytes = predict(x);
    let seconds = {
        let t: Vec<(f64, f64)> = same.iter().filter(|r| r.seconds > 0.0).map(|r| (r.cells as f64, r.seconds)).collect();
        fit(&t).map(|(ta, tp)| (ta + tp.clamp(0.5, 4.0) * x.ln()).exp()).filter(|s| s.is_finite())
    };
    let (lo, hi) = same.iter().fold((u64::MAX, 0), |(lo, hi), r| (lo.min(r.cells), hi.max(r.cells)));
    Some(Estimate {
        bytes: bytes as u64,
        bytes_with_margin: (bytes * (1.0 + worst.max(0.10))) as u64,
        seconds,
        runs: same.len(),
        extrapolated: cells < lo || cells > hi,
        exponent: p,
    })
}

/// What the writer thread of a history is asked to do.
enum Write {
    Line(PathBuf, String),
    Flush(std::sync::mpsc::Sender<()>),
}

/// The thread that appends lines for a `History`, one at a time in the order they were added.
struct Writer {
    tx: std::sync::mpsc::Sender<Write>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl Writer {
    fn start() -> Writer {
        let (tx, rx) = std::sync::mpsc::channel::<Write>();
        let handle = std::thread::spawn(move || {
            while let Ok(w) = rx.recv() {
                match w {
                    Write::Line(p, line) => {
                        use std::io::Write as _;
                        if let Some(d) = p.parent() {
                            let _ = std::fs::create_dir_all(d);
                        }
                        if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&p) {
                            let _ = f.write_all(line.as_bytes());
                        }
                    }
                    Write::Flush(ack) => {
                        let _ = ack.send(());
                    }
                }
            }
        });
        Writer { tx, handle: Some(handle) }
    }
}

/// The rows an engine has learned from, in memory and, when a file is set, appended to it as JSON lines (one `write` per line, like the run
/// table, so two programs can append to one file). The append is made by a thread of the history, not by the caller: a file open and write
/// took 10 ms on the thread that draws (productivity plan P3), and the rows in memory, which the estimates read, are updated at once.
/// `flush` waits until every line added so far is in the file; dropping the history flushes too.
#[derive(Default)]
pub struct History {
    rows: Mutex<Vec<Row>>,
    file: Mutex<Option<PathBuf>>,
    writer: Mutex<Option<Writer>>,
}

impl Drop for History {
    fn drop(&mut self) {
        if let Some(mut w) = self.writer.lock().ok().and_then(|mut g| g.take()) {
            drop(w.tx);
            if let Some(h) = w.handle.take() {
                let _ = h.join();
            }
        }
    }
}

impl History {
    /// Keeps the rows in `path` from now on and reads the ones it holds; a line that cannot be read is skipped.
    pub fn set_file(&self, path: &Path) {
        self.flush();
        let mut rows = vec![];
        if let Ok(text) = std::fs::read_to_string(path) {
            rows.extend(text.lines().filter_map(|l| serde_json::from_str::<Row>(l).ok()));
        }
        *self.rows.lock().expect("history") = rows;
        *self.file.lock().expect("history") = Some(path.to_path_buf());
    }

    /// Reads the file again: another program (the queue worker of the workspace) appends rows to it too.
    pub fn reload(&self) {
        let file = self.file.lock().expect("history").clone();
        if let Some(p) = file {
            self.set_file(&p);
        }
    }

    pub fn rows(&self) -> Vec<Row> {
        self.rows.lock().expect("history").clone()
    }

    /// Waits until every line added so far has been written to the file.
    pub fn flush(&self) {
        let tx = self.writer.lock().expect("history").as_ref().map(|w| w.tx.clone());
        if let Some(tx) = tx {
            let (ack, done) = std::sync::mpsc::channel();
            if tx.send(Write::Flush(ack)).is_ok() {
                let _ = done.recv();
            }
        }
    }

    pub fn add(&self, row: Row) {
        let file = self.file.lock().expect("history").clone();
        if let Some(p) = file {
            let line = format!("{}\n", json!(row));
            let mut w = self.writer.lock().expect("history");
            let _ = w.get_or_insert_with(Writer::start).tx.send(Write::Line(p, line));
        }
        self.rows.lock().expect("history").push(row);
    }

    pub fn estimate(&self, chain: &str, dim: u8, cells: u64) -> Option<Estimate> {
        estimate(&self.rows.lock().expect("history"), chain, dim, cells)
    }
}

/// What the machine says about its memory: `(total, available)` in bytes. Replaceable, so a test can play a machine that is full.
pub type MemoryProbe = std::sync::Arc<dyn Fn() -> Option<(u64, u64)> + Send + Sync>;

/// How much of the machine's memory a solve may be expected to take at most: above this fraction of the physical memory a solve is refused,
/// because it could not run even with nothing else on the machine.
pub const REFUSE_ABOVE: f64 = 0.9;

fn gb(b: u64) -> String {
    format!("{:.1} GB", b as f64 / 1e9)
}

/// Validates a memory-wait timeout in seconds (`finite` and `>= 0.0`).
pub fn validate_memory_wait_timeout(timeout_s: f64, source: &str) -> Result<f64, nyx_physics::Error> {
    if !(timeout_s.is_finite() && timeout_s >= 0.0) {
        return Err(nyx_physics::Error::invalid(format!(
            "{source} must be a finite number of seconds >= 0, got {timeout_s}"
        )));
    }
    Ok(timeout_s)
}

/// Reads `NYX_MEMORY_WAIT_TIMEOUT_S` when set to a non-empty value.
pub fn env_memory_wait_timeout_s() -> Result<Option<f64>, nyx_physics::Error> {
    let Ok(raw) = std::env::var("NYX_MEMORY_WAIT_TIMEOUT_S") else {
        return Ok(None);
    };
    let s = raw.trim();
    if s.is_empty() {
        return Ok(None);
    }
    let val: f64 = s.parse().map_err(|_| {
        nyx_physics::Error::invalid(format!(
            "`NYX_MEMORY_WAIT_TIMEOUT_S` must be a finite number of seconds >= 0, got '{s}'"
        ))
    })?;
    validate_memory_wait_timeout(val, "`NYX_MEMORY_WAIT_TIMEOUT_S`").map(Some)
}

/// Resolves the effective memory-wait timeout in priority order:
/// explicit `timeout_s` > `NYX_MEMORY_WAIT_TIMEOUT_S` env var > `caller_default`.
pub fn resolve_memory_wait_timeout(
    explicit: Option<f64>,
    caller_default: Option<f64>,
) -> Result<Option<f64>, nyx_physics::Error> {
    if let Some(t) = explicit {
        return validate_memory_wait_timeout(t, "`memory_wait_timeout_s`").map(Some);
    }
    if let Some(d) = caller_default {
        validate_memory_wait_timeout(d, "default `memory_wait_timeout_s`")?;
    }
    if let Some(env_t) = env_memory_wait_timeout_s()? {
        return Ok(Some(env_t));
    }
    Ok(caller_default)
}

/// Lets a solve start only when its estimate fits. Above `REFUSE_ABOVE` of the physical memory it fails at once, unless `force`: it cannot be
/// waited for. When it fits the machine but not what is free now (another program, another solve), it waits, shows "waiting for memory:
/// needs about 9.0 GB, 6.0 GB free" as its stage and re-asks every quarter second (clamped to any remaining timeout), until it fits, times out,
/// or is cancelled. Priority of timeout: explicit `timeout_s` > `NYX_MEMORY_WAIT_TIMEOUT_S` > infinite wait (`None`). A machine that does not
/// say how much memory it has admits everything.
pub fn admit(
    est: &Estimate,
    force: bool,
    probe: &MemoryProbe,
    progress: &std::sync::Arc<femcore::progress::Progress>,
    timeout_s: Option<f64>,
) -> Result<(), nyx_physics::Error> {
    let effective_timeout = resolve_memory_wait_timeout(timeout_s, None)?;
    if force {
        return Ok(());
    }
    let Some((total, _)) = probe() else { return Ok(()) };
    let need = est.bytes_with_margin;
    if need as f64 > total as f64 * REFUSE_ABOVE {
        return Err(nyx_physics::Error::invalid(format!(
            "this solve needs about {} (from {} earlier run{}{}), more than this machine can hold ({} of memory in all), so it was not started; use a coarser mesh, or `force: true` to try anyway",
            gb(need),
            est.runs,
            if est.runs == 1 { "" } else { "s" },
            if est.extrapolated { ", extrapolated" } else { "" },
            gb(total)
        )));
    }
    let _installed = femcore::progress::install(progress);
    let started = std::time::Instant::now();
    let deadline = effective_timeout.map(std::time::Duration::from_secs_f64);
    let mut shown = String::new();
    let mut stage = None;
    loop {
        let (total_now, free) = match probe() {
            Some((t, f)) => (t, f),
            None => return Ok(()),
        };
        if need <= free {
            return Ok(());
        }
        if progress.is_cancelled() {
            return Err(nyx_physics::Error::Cancelled);
        }
        if let Some(limit) = deadline {
            if started.elapsed() >= limit {
                let timeout_val = effective_timeout.unwrap_or(0.0);
                return Err(nyx_physics::Error::invalid(format!(
                    "timed out after {timeout_val:.3} s (timeout_s = {timeout_val}) waiting for memory: need {} ({need} B), free {} ({free} B), total {} ({total_now} B); free memory, raise `memory_wait_timeout_s`, use a coarser mesh, or pass `force: true`",
                    gb(need),
                    gb(free),
                    gb(total_now),
                )));
            }
        }
        let text = format!("waiting for memory: needs about {}, {} free", gb(need), gb(free));
        if text != shown {
            drop(stage.take());
            stage = Some(femcore::progress::stage(text.clone()));
            shown = text;
        }
        let poll = std::time::Duration::from_millis(250);
        let sleep_dur = match deadline {
            Some(limit) => poll.min(limit.saturating_sub(started.elapsed())),
            None => poll,
        };
        if !sleep_dur.is_zero() {
            std::thread::sleep(sleep_dur);
        }
    }
}

/// The chain's key: codes sorted and joined.
pub fn chain_key(targets: &[String]) -> String {
    let mut t: Vec<&str> = targets.iter().map(String::as_str).collect();
    t.sort_unstable();
    t.dedup();
    t.join(",")
}

/// Calibration chain key aware of `ems:driven`, `sms:driven`, and coupled `sms:piezo` `(solver_backend, order)` profiles (`FR-33`),
/// keeping pure-elastic `sms` and non-driven `ems` targets identical to [`chain_key`] so existing calibration history stays intact.
pub fn chain_key_for_project(project: &crate::project::Project, targets: &[String]) -> String {
    let mut tokens: Vec<String> = targets
        .iter()
        .map(|t| {
            if let Some(spec) = project.spec(t) {
                if spec.physics == "ems" && spec.study.kind == "driven" {
                    let solver = spec
                        .study
                        .params
                        .get("solver")
                        .and_then(|v| v.as_str())
                        .map(str::trim)
                        .filter(|s| !s.is_empty())
                        .unwrap_or("auto");
                    let order = spec.order.clamp(1, 2);
                    return format!("ems:driven:{solver}:p{order}");
                }
                if spec.physics == "sms" {
                    let order = spec.order.clamp(1, 2);
                    if spec.study.kind == "driven" {
                        let solver = spec
                            .study
                            .params
                            .get("solver")
                            .and_then(|v| v.as_str())
                            .map(str::trim)
                            .filter(|s| !s.is_empty())
                            .unwrap_or("auto");
                        return format!("sms:driven:{solver}:p{order}");
                    }
                    let piezo_param = spec.params.get("piezoelectric").and_then(|v| v.as_str()).unwrap_or("auto");
                    let has_elec_bc = spec.bcs.iter().any(|bc| {
                        matches!(
                            bc.kind.as_str(),
                            "ground" | "potential" | "surface_charge" | "floating_potential" | "conductor" | "terminal"
                        )
                    });
                    if piezo_param == "on" || (piezo_param != "off" && has_elec_bc) {
                        return format!("sms:piezo:p{order}");
                    }
                }
            }
            t.clone()
        })
        .collect();
    tokens.sort_unstable();
    tokens.dedup();
    tokens.join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(cells: u64, gb: f64) -> Row {
        Row {
            chain: "ems".into(),
            dim: 3,
            cells,
            dofs: cells * 6,
            peak_working_set: (gb * 1e9) as u64,
            peak_commit: (gb * 1.3e9) as u64,
            seconds: gb * 3.0,
        }
    }

    #[test]
    fn no_estimate_without_two_sizes() {
        assert!(estimate(&[], "ems", 3, 5000).is_none());
        assert!(estimate(&[row(6672, 2.31)], "ems", 3, 5000).is_none());
        assert!(
            estimate(&[row(6672, 2.31), row(6672, 2.4)], "ems", 3, 5000).is_none(),
            "the same size twice is one size"
        );
        // another chain or dimension is another problem
        let rows = [row(6672, 2.31), row(18120, 10.67)];
        assert!(estimate(&rows, "dfs", 3, 9000).is_none());
        assert!(estimate(&rows, "ems", 2, 9000).is_none());
    }

    #[test]
    fn the_fit_goes_through_the_p0_points_and_says_it_extrapolates() {
        // the two measured points of the 3D cavity (6,672 cells 2.31 GB; 18,120 cells 10.67 GB): exponent about 1.8 in cells
        let rows = [row(6672, 2.31), row(18120, 10.67)];
        let at_a = estimate(&rows, "ems", 3, 6672).unwrap();
        assert!((at_a.bytes as f64 / 2.31e9 - 1.0).abs() < 0.02, "{at_a:?}");
        assert!(!at_a.extrapolated);
        let big = estimate(&rows, "ems", 3, 60000).unwrap();
        assert!(big.extrapolated && big.bytes > 10_670_000_000, "{big:?}");
        assert!(big.bytes_with_margin as f64 >= big.bytes as f64 * 1.099);
        assert_eq!(big.runs, 2);
    }

    #[test]
    fn the_history_is_kept_in_a_file_other_programs_can_append_to() {
        let dir = std::env::temp_dir().join(format!("nyx-history-{}", std::process::id()));
        let file = dir.join("h.jsonl");
        let h = History::default();
        h.set_file(&file);
        h.add(row(1000, 0.1));
        h.add(row(4000, 0.5));
        h.flush(); // the file is written by a thread (the frame that takes a result never waits for it)
        let again = History::default();
        again.set_file(&file);
        assert_eq!(again.rows().len(), 2);
        assert!(again.estimate("ems", 3, 2000).is_some());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn driven_ems_keys_calibration_by_solver_backend_and_order() {
        let mut p = crate::project::Project::new("Driven", 3);
        let mut s = nyx_physics::ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({ "f_min": 4e9, "f_max": 8e9, "solver": "blr" });
        p.physics.push(s);
        assert_eq!(chain_key_for_project(&p, &["ems".into()]), "ems:driven:blr:p1");

        p.physics[0].order = 2;
        p.physics[0].study.params = json!({ "f_min": 4e9, "f_max": 8e9, "solver": "ams_fgmres" });
        assert_eq!(chain_key_for_project(&p, &["ems".into()]), "ems:driven:ams_fgmres:p2");

        // Non-driven ems keeps plain "ems" key
        p.physics[0] = nyx_physics::ProblemSpec::new("ems", 1, "cavity");
        assert_eq!(chain_key_for_project(&p, &["ems".into()]), "ems");

        // Legacy "ems" or "ems:driven:direct:p2" history rows never block "blr" / "ams_fgmres" / "lblt"
        let rows = [
            row(100_000, 40.0),
            row(200_000, 95.0),
            Row { chain: "ems:driven:direct:p2".into(), ..row(100_000, 45.0) },
            Row { chain: "ems:driven:direct:p2".into(), ..row(200_000, 110.0) },
        ];
        assert!(estimate(&rows, "ems:driven:direct:p2", 3, 150_000).is_some());
        assert!(estimate(&rows, "ems:driven:lblt:p1", 3, 150_000).is_none());
        assert!(estimate(&rows, "ems:driven:blr:p1", 3, 150_000).is_none());
        assert!(estimate(&rows, "ems:driven:ams_fgmres:p1", 3, 150_000).is_none());
    }
}
