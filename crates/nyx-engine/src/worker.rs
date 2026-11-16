//! Solving in a child process (productivity plan P4, ADR 86). A solve that needs more memory than the machine has, or a factorisation that
//! nobody can interrupt, must not take the window or the MCP server with it. A *worker* is a small program (`nyx-worker`, built from this crate)
//! that solves one frozen project and writes the results file; the engine that asked runs it as a child, shows its stage log as its own, kills
//! it when the job is cancelled (at once, even in the middle of a factorisation) and reads the results back. A crash or an out-of-memory
//! failure then fails one job, not the process the person is working in.
//!
//! **What the worker is given** is everything the in-process job would have used: the frozen project (JSON), the materials of the library the
//! parent sees (TOML, by name) and the chain of physics to solve. It builds a scratch `Engine` from them and runs the same `start_solve` path,
//! so the numbers cannot differ (ADR 11; `tests/worker.rs` compares them to the bit over the examples). **What it sends back** is on stdout,
//! one line per message, each starting with `@nyx ` and holding JSON (anything else on stdout is ignored): `status` (the stage, how long it has
//! been running, the fraction done and the stage log when it changed) every 100 ms, and one `done` (the error if the solve failed, whether
//! results were written, the peak memory the process used). The results are the HDF5 file `results.h5` in the job's folder (ADR 62, 63), the
//! same round trip as opening a results file. The worker exits when its standard input closes, so it never outlives the process that started it.

use crate::engine::Engine;
use crate::jobs::WorkerRun;
use crate::project::Project;
use femcore::progress::Progress;
use nyx_physics::{Error, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Every line of the worker's standard output that is a message starts with this.
pub const PREFIX: &str = "@nyx ";

/// What a worker is asked to solve.
#[derive(Clone)]
pub(crate) struct RemoteSpec {
    pub worker: PathBuf,
    pub project_json: String,
    /// every material of the library the job sees, as TOML by name
    pub materials: BTreeMap<String, String>,
    pub targets: Vec<String>,
}

/// What a worker handed back.
pub(crate) struct RemoteDone {
    /// physics code and result, as written (sorted by code)
    #[cfg_attr(not(feature = "hdf5"), allow(dead_code))]
    pub results: Vec<(String, SolveResult)>,
    #[cfg_attr(not(feature = "hdf5"), allow(dead_code))]
    pub mesh: Option<femcore::mesh::SimplexMesh>,
    /// why the chain stopped, when a physics failed after others were solved or before any was
    pub failed: Option<Error>,
    pub run: WorkerRun,
}

/// The materials of a library as TOML by name, what a worker rebuilds the library from.
pub fn material_texts(lib: &femcore::material::Library) -> BTreeMap<String, String> {
    lib.names().iter().filter_map(|n| lib.get(n)).map(|m| (m.name.clone(), m.to_toml())).collect()
}

fn error_from(code: &str, message: String) -> Error {
    match code {
        "invalid" => Error::Invalid(message),
        "unsupported" => Error::Unsupported(message),
        "non_finite" => Error::NonFinite(message),
        "cancelled" => Error::Cancelled,
        "not_yet_moved" => Error::NotYetMoved(message),
        _ => Error::Solver(message),
    }
}

// ------------------------------------------------------------------------------------------------------------------ the worker

/// The worker program: `nyx-worker --job FOLDER`. Returns the exit code (0 when it ran, whatever the solve's outcome: the outcome is in the
/// `done` message).
pub fn serve(args: &[String]) -> i32 {
    let dir = match args {
        [flag, dir] if flag == "--job" => PathBuf::from(dir),
        _ => {
            eprintln!("usage: nyx-worker --job FOLDER   (a folder written by an engine: project.nyx.json and job.json; this program is started by nyx, not by hand)");
            return 2;
        }
    };
    // the process that started this one is gone when its end of our standard input closes: stop, so no solve runs on for nobody
    std::thread::spawn(|| {
        let mut sink = [0u8; 256];
        let mut stdin = std::io::stdin();
        while let Ok(n) = std::io::Read::read(&mut stdin, &mut sink) {
            if n == 0 {
                std::process::exit(3);
            }
        }
        std::process::exit(3);
    });
    match serve_job(&dir) {
        Ok(()) => 0,
        Err(message) => {
            say(
                &json!({ "done": { "ok": false, "error": { "code": "solver", "message": message }, "has_results": false } }),
            );
            1
        }
    }
}

fn say(v: &Value) {
    use std::io::Write;
    let mut out = std::io::stdout().lock();
    let _ = writeln!(out, "{PREFIX}{v}");
    let _ = out.flush();
}

#[cfg(feature = "hdf5")]
fn serve_job(dir: &Path) -> Result<(), String> {
    let read =
        |name: &str| std::fs::read_to_string(dir.join(name)).map_err(|e| format!("{}: {e}", dir.join(name).display()));
    let job: Value = serde_json::from_str(&read("job.json")?).map_err(|e| format!("job.json: {e}"))?;
    let targets: Vec<String> = job["targets"]
        .as_array()
        .ok_or("job.json has no targets")?
        .iter()
        .filter_map(|t| t.as_str().map(String::from))
        .collect();
    let texts: BTreeMap<String, String> = job["materials"]
        .as_object()
        .ok_or("job.json has no materials")?
        .iter()
        .filter_map(|(k, v)| v.as_str().map(|t| (k.clone(), t.to_string())))
        .collect();
    let project = Project::from_json(&read("project.nyx.json")?).map_err(|e| e.to_string())?;
    let library = femcore::material::Library::default().with_texts(&texts).map_err(|errs| {
        format!("the materials do not read: {}", errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; "))
    })?;
    let mut e = Engine::with(crate::registry(), Arc::new(crate::meshing::FemcoreMesher));
    e.set_base_library(library);
    e.set_project(project);
    if std::env::var("NYX_WORKER_DRILL").as_deref() == Ok("abort") {
        // a crash drill for the tests: the process dies the way one does when the allocator gives up, with no answer
        std::process::abort();
    }
    let id = e.start_solve_with(Some(targets), &crate::queue::QueueOpts::default()).map_err(|e| e.to_string())?;
    let mut sent_log = 0usize;
    let mut sent_last = String::new();
    let outcome = loop {
        match e.wait_timeout(id, std::time::Duration::from_millis(100)) {
            Ok(v) if v.get("waiting") == Some(&json!(true)) => {
                if let Some(p) = e.jobs().progress(id) {
                    let st = p.status();
                    let log = p.log();
                    let changed = log.len() != sent_log || log.last().map(String::as_str) != Some(sent_last.as_str());
                    let mut m =
                        json!({ "stage": st.stage, "stage_seconds": st.stage_seconds, "fraction": st.fraction });
                    if changed {
                        sent_log = log.len();
                        sent_last = log.last().cloned().unwrap_or_default();
                        m["log"] = json!(log);
                    }
                    say(&json!({ "status": m }));
                }
            }
            other => break other,
        }
    };
    let error = outcome.err();
    let has_results = match e.results_snapshot() {
        Ok(s) => {
            let refs: Vec<(String, &SolveResult)> = s.results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();
            crate::h5io::write_results(&dir.join("results.h5"), &s.project, &s.model.solver_mesh(), &refs)?;
            true
        }
        Err(_) => false,
    };
    let (ws, commit) = crate::process::peak_memory().unzip();
    say(&json!({ "done": {
        "ok": error.is_none(),
        "error": error.as_ref().map(|e| json!({ "code": e.code(), "message": e.to_string() })),
        "has_results": has_results,
        "peak_working_set": ws, "peak_commit": commit,
    } }));
    Ok(())
}

#[cfg(not(feature = "hdf5"))]
fn serve_job(_dir: &Path) -> Result<(), String> {
    Err("this build of nyx-worker has no HDF5 support, so it cannot hand results back".into())
}

// ------------------------------------------------------------------------------------------------------------------ the parent

/// A folder for one job's files, removed when it goes out of scope.
struct JobDir(PathBuf);

impl Drop for JobDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn new_job_dir() -> Result<JobDir, Error> {
    static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    let dir = std::env::temp_dir().join(format!("nyx-job-{}-{n}", std::process::id()));
    std::fs::create_dir_all(&dir)
        .map_err(|e| Error::Solver(format!("could not make the job folder {}: {e}", dir.display())))?;
    Ok(JobDir(dir))
}

fn exit_words(status: &std::process::ExitStatus, stderr: &str) -> String {
    let tail: Vec<&str> = stderr.lines().rev().take(3).collect::<Vec<_>>().into_iter().rev().collect();
    let tail = if tail.is_empty() { String::new() } else { format!("; it said: {}", tail.join(" / ")) };
    match status.code() {
        Some(c) => format!(
            "exit code {c}{}{tail}",
            if c < 0 || c as u32 >= 0xC000_0000 {
                format!(" (0x{:08X}: the process crashed or ran out of memory)", c as u32)
            } else {
                String::new()
            }
        ),
        None => format!("it was ended by a signal{tail}"),
    }
}

/// Solves `spec` in a worker and waits for it, mirroring its stage log into `progress`. A cancel on `progress` kills the worker at once.
#[cfg(feature = "hdf5")]
pub(crate) fn run_remote(spec: &RemoteSpec, progress: &Arc<Progress>) -> Result<RemoteDone, Error> {
    use std::io::BufRead;
    use std::process::{Command, Stdio};
    use std::sync::mpsc;
    use std::time::Duration;
    if progress.is_cancelled() {
        return Err(Error::Cancelled);
    }
    let dir = new_job_dir()?;
    let io = |what: &str, e: std::io::Error| Error::Solver(format!("{what}: {e}"));
    std::fs::write(dir.0.join("project.nyx.json"), &spec.project_json)
        .map_err(|e| io("writing the job's project", e))?;
    std::fs::write(dir.0.join("job.json"), json!({ "targets": spec.targets, "materials": spec.materials }).to_string())
        .map_err(|e| io("writing the job's description", e))?;
    let stderr_path = dir.0.join("stderr.txt");
    let stderr = std::fs::File::create(&stderr_path).map_err(|e| io("making the worker's log file", e))?;
    let mut cmd = Command::new(&spec.worker);
    cmd.arg("--job").arg(&dir.0).stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::from(stderr));
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        // no console window flashing up beside the application
        cmd.creation_flags(0x0800_0000);
    }
    let mut child = cmd.spawn().map_err(|e| {
        Error::Solver(format!(
            "could not start the worker program {}: {e} (the solve was not run; run it in this process instead)",
            spec.worker.display()
        ))
    })?;
    // keeping this end open is what tells the worker its parent is alive
    let _alive = child.stdin.take();
    let stdout = child.stdout.take().ok_or_else(|| Error::Solver("the worker has no output".into()))?;
    let (tx, rx) = mpsc::channel::<String>();
    std::thread::spawn(move || {
        for line in std::io::BufReader::new(stdout).lines().map_while(Result::ok) {
            if tx.send(line).is_err() {
                break;
            }
        }
    });
    let mut done: Option<Value> = None;
    // what the job logged before the worker started (the mesh, a wait for memory) stays above the worker's own log
    let base = progress.log();
    loop {
        match rx.recv_timeout(Duration::from_millis(40)) {
            Ok(line) => {
                if let Some(Ok(m)) = line.strip_prefix(PREFIX).map(serde_json::from_str::<Value>) {
                    if let Some(s) = m.get("status") {
                        let log =
                            s.get("log").and_then(|l| serde_json::from_value::<Vec<String>>(l.clone()).ok()).map(|l| {
                                let mut all = base.clone();
                                all.extend(l);
                                all
                            });
                        progress.mirror(
                            s.get("stage").and_then(Value::as_str).unwrap_or(""),
                            s.get("stage_seconds").and_then(Value::as_f64).unwrap_or(0.0),
                            s.get("fraction").and_then(Value::as_f64),
                            log,
                        );
                    } else if let Some(d) = m.get("done") {
                        done = Some(d.clone());
                    }
                }
            }
            Err(mpsc::RecvTimeoutError::Timeout) => {}
            Err(mpsc::RecvTimeoutError::Disconnected) => break,
        }
        if progress.is_cancelled() {
            // the whole point of a worker: no factorisation has to be waited for
            let _ = child.kill();
            let _ = child.wait();
            return Err(Error::Cancelled);
        }
    }
    let status = child.wait().map_err(|e| io("waiting for the worker", e))?;
    let Some(d) = done else {
        let said = std::fs::read_to_string(&stderr_path).unwrap_or_default();
        return Err(Error::Solver(format!(
            "the worker process stopped without an answer ({}): the solve needed more memory than was free, or the process was ended. Nothing else was affected",
            exit_words(&status, &said)
        )));
    };
    let failed = d.get("error").filter(|e| !e.is_null()).map(|e| {
        error_from(
            e.get("code").and_then(Value::as_str).unwrap_or("solver"),
            e.get("message").and_then(Value::as_str).unwrap_or("the worker failed").to_string(),
        )
    });
    let (results, mesh) = if d.get("has_results").and_then(Value::as_bool) == Some(true) {
        let _stage = femcore::progress::install(progress);
        let _s = femcore::progress::stage("reading the results of the worker");
        let stored = crate::h5io::read_results(&dir.0.join("results.h5")).map_err(Error::Solver)?;
        (stored.results, Some(stored.mesh))
    } else {
        (vec![], None)
    };
    let run = WorkerRun {
        peak_working_set: d.get("peak_working_set").and_then(Value::as_u64),
        peak_commit: d.get("peak_commit").and_then(Value::as_u64),
    };
    Ok(RemoteDone { results, mesh, failed, run })
}

#[cfg(not(feature = "hdf5"))]
pub(crate) fn run_remote(_spec: &RemoteSpec, _progress: &Arc<Progress>) -> Result<RemoteDone, Error> {
    Err(Error::Unsupported(
        "this build of nyx has no HDF5 support, which a worker needs to hand results back; solve in this process"
            .into(),
    ))
}

/// True when the mesh a worker stored is the mesh the parent holds (meshing is deterministic, ADR 55): the check that parent and worker solved
/// the same problem.
pub(crate) fn same_mesh(a: &femcore::mesh::SimplexMesh, b: &femcore::mesh::SimplexMesh) -> bool {
    a.dim == b.dim
        && a.points == b.points
        && a.cells == b.cells
        && a.cell_labels == b.cell_labels
        && a.facets == b.facets
        && a.facet_labels == b.facet_labels
}
