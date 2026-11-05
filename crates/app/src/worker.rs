//! Background jobs: building meshes, opening files and solving, so the window stays responsive.
//!
//! Each job runs on its own thread and sends its result back through a channel together with the job number (`Job`); the
//! window picks it up in `App::poll`. A solve calls `nyx_engine::legacy::solve_with` with the job's `Progress`, so it can be watched and
//! cancelled, then writes the results file requested in the study and the autosave of the run (project XML, stage log and, in builds
//! with HDF5, the results) to the runs folder. Files are named `temp_untitled_YYYYMMDD_HHMMSS` with the local time.

use crate::state::GeometrySpec;
use eframe::egui;
use femcore::fields::Display;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::mesh::view::ViewMesh;
use std::path::{Path, PathBuf};
use std::sync::mpsc::Sender;
use femcore::progress::Progress;
use std::sync::Arc;

pub struct Built {
    pub mesh: Arc<SimplexMesh>,
    pub view: Arc<ViewMesh>,
}

pub struct Solved {
    pub result: SolveResult,
    pub display: Display,
    pub state: SolveState,
    /// outcome of writing the results file, when one was requested
    pub output: Option<Result<String, String>>,
    /// the polarization analysis of the modes of an EM-wave study (a function of the result, not part of it)
    pub analysis: Option<nyx_engine::analysis::ModeAnalysis>,
}

pub enum Msg {
    Built(Result<Built, String>),
    Solved { physics: Physics, res: Result<Solved, String> },
    Script(Result<String, String>),
    H5(String, Result<(Arc<SimplexMesh>, Solved), String>),
}

/// A background job: its number (results of earlier jobs are told apart by it), the channel back to the window and
/// the progress/cancel handle shared with the solver.
#[derive(Clone)]
pub struct Job {
    pub id: u64,
    pub tx: Sender<(u64, Msg)>,
    pub ctx: egui::Context,
    pub progress: Arc<Progress>,
}

fn send(job: &Job, m: Msg) {
    let _ = job.tx.send((job.id, m));
    job.ctx.request_repaint();
}

pub fn spawn_build(spec: GeometrySpec, problem: Option<ProblemSpec>, job: Job) {
    std::thread::spawn(move || {
        let r = spec.build_named_for(problem.as_ref()).map(|g| g.mesh).map(|m| {
            let view = Arc::new(ViewMesh::from_simplex(&m));
            Built { mesh: Arc::new(m), view }
        });
        send(&job, Msg::Built(r));
    });
}

/// Every run is saved to a folder: `temp_untitled_<date>_<time>.hdf5` with its XDMF description and the project as
/// `temp_untitled_<date>_<time>.nyx.xml`.
#[derive(Clone, Debug)]
pub struct AutoSave {
    pub dir: PathBuf,
    /// the project as XML (made on the interface thread)
    pub project_xml: Result<String, String>,
}

/// Local date and time as `YYYYMMDD_HHMMSS`: the Windows clock when available (no helper process is started), else UTC.
fn timestamp() -> String {
    #[cfg(windows)]
    {
        #[repr(C)]
        #[derive(Default)]
        struct SystemTime {
            year: u16,
            month: u16,
            day_of_week: u16,
            day: u16,
            hour: u16,
            minute: u16,
            second: u16,
            milliseconds: u16,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn GetLocalTime(out: *mut SystemTime);
        }
        let mut st = SystemTime::default();
        // SAFETY: GetLocalTime only writes the SYSTEMTIME structure it is given.
        unsafe { GetLocalTime(&mut st) };
        if st.year >= 2000 {
            return format!("{:04}{:02}{:02}_{:02}{:02}{:02}", st.year, st.month, st.day, st.hour, st.minute, st.second);
        }
    }
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    stamp_from_unix(secs)
}

/// `YYYYMMDD_HHMMSS` of a Unix time in UTC (civil date from days since 1970-01-01, Howard Hinnant's algorithm).
fn stamp_from_unix(secs: u64) -> String {
    let (days, rem) = (secs / 86_400, secs % 86_400);
    let z = days as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let (y, doy) = (yoe + era * 400, doe - (365 * yoe + yoe / 4 - yoe / 100));
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{y:04}{m:02}{d:02}_{:02}{:02}{:02}", rem / 3600, rem % 3600 / 60, rem % 60)
}
fn autosave_run(auto: &AutoSave, mesh: &SimplexMesh, spec: &ProblemSpec, result: &SolveResult) -> Result<String, String> {
    std::fs::create_dir_all(&auto.dir).map_err(|e| format!("cannot create {}: {e}", auto.dir.display()))?;
    let stem = format!("temp_untitled_{}", timestamp());
    let mut written: Vec<String> = vec![];
    match &auto.project_xml {
        Ok(xml) => {
            let path = auto.dir.join(format!("{stem}.nyx.xml"));
            std::fs::write(&path, xml).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
            written.push(path.display().to_string());
        }
        Err(e) => return Err(format!("project not saved: {e}")),
    }
    // the plain-language log of the run (what was done, how long each stage took)
    let log_path = auto.dir.join(format!("{stem}.log.txt"));
    std::fs::write(&log_path, format!("{}\n", result.log.join("\n"))).map_err(|e| format!("cannot write {}: {e}", log_path.display()))?;
    written.push(log_path.display().to_string());
    #[cfg(feature = "hdf5")]
    {
        let path = auto.dir.join(format!("{stem}.hdf5"));
        femcore::h5io::write_outputs(&path, femcore::h5io::H5Mode::Serial, mesh, Some(spec), result, true).map_err(|e| format!("results not saved: {e}"))?;
        written.push(path.display().to_string());
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (mesh, spec, result);
        written.push("(no HDF5 in this build: only the project was saved)".into());
    }
    Ok(format!("saved {}", written.join(" and ")))
}

pub fn spawn_solve(mesh: Arc<SimplexMesh>, spec: ProblemSpec, mut state: SolveState, auto: Option<AutoSave>, job: Job) {
    std::thread::spawn(move || {
        let physics = spec.physics;
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| nyx_engine::legacy::solve_with_analysis(&mesh, &spec, &mut state, &job.progress)));
        let res = match r {
            Ok(Ok((result, analysis))) => {
                let display = Display::new(&mesh, result.order);
                let mut output = write_output(&mesh, &spec, &result);
                if let Some(a) = &auto {
                    let saved = autosave_run(a, &mesh, &spec, &result);
                    output = match (output, saved) {
                        (None, s) => Some(s),
                        (Some(Ok(o)), Ok(s)) => Some(Ok(format!("{o}; {s}"))),
                        (Some(Err(e)), _) => Some(Err(e)),
                        (Some(Ok(_)), Err(e)) => Some(Err(e)),
                    };
                }
                Ok(Solved { result, display, state, output, analysis })
            }
            Ok(Err(e)) => Err(e),
            Err(_) => Err("the solver stopped unexpectedly (internal error)".into()),
        };
        send(&job, Msg::Solved { physics, res });
    });
}

/// Runs an optional user Python script: `python script.py project.json [result.csv]`. Python is
/// never needed by the application itself.
pub fn spawn_script(script: PathBuf, project_json: String, csv: Option<String>, job: Job) {
    std::thread::spawn(move || {
        let run = || -> Result<String, String> {
            let dir = std::env::temp_dir().join(format!("nyx-script-{}", std::process::id()));
            std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            let pj = dir.join("project.json");
            std::fs::write(&pj, project_json).map_err(|e| e.to_string())?;
            let mut args: Vec<std::ffi::OsString> = vec![script.clone().into(), pj.clone().into()];
            if let Some(c) = csv {
                let cp = dir.join("result.csv");
                std::fs::write(&cp, c).map_err(|e| e.to_string())?;
                args.push(cp.into());
            }
            let mut candidates: Vec<(String, Vec<String>)> = vec![];
            if let Ok(p) = std::env::var("NYX_PYTHON") {
                candidates.push((p, vec![]));
            }
            candidates.push(("python".into(), vec![]));
            candidates.push(("py".into(), vec!["-3".into()]));
            candidates.push(("python3".into(), vec![]));
            for (exe, pre) in candidates {
                let out = std::process::Command::new(&exe).args(&pre).args(&args).current_dir(script.parent().unwrap_or(Path::new("."))).output();
                match out {
                    Ok(o) => {
                        let mut text = String::from_utf8_lossy(&o.stdout).into_owned();
                        let err = String::from_utf8_lossy(&o.stderr);
                        if !err.trim().is_empty() {
                            text += "\n--- stderr ---\n";
                            text += &err;
                        }
                        if !o.status.success() {
                            text += &format!("\n(exit status {})", o.status);
                        }
                        return Ok(text);
                    }
                    Err(_) => continue,
                }
            }
            Err("Python was not found. It is optional: it is only used to run your own scripts (set NYX_PYTHON to the interpreter to use).".into())
        };
        send(&job, Msg::Script(run()));
    });
}

/// Writes the results file requested by the simulation's output settings.
pub fn write_output(mesh: &SimplexMesh, spec: &ProblemSpec, result: &SolveResult) -> Option<Result<String, String>> {
    use femcore::problem::OutputMode;
    if spec.output.mode == OutputMode::None {
        return None;
    }
    if spec.output.path.trim().is_empty() {
        return Some(Err("an output file is selected but no path is given".into()));
    }
    #[cfg(feature = "hdf5")]
    {
        let mode = if spec.output.mode == OutputMode::Hdf5Parallel { femcore::h5io::H5Mode::Parallel } else { femcore::h5io::H5Mode::Serial };
        let r = femcore::h5io::write_outputs(Path::new(spec.output.path.trim()), mode, mesh, Some(spec), result, true);
        return Some(r.map(|files| {
            if files.is_empty() {
                format!("results written by rank 0 ({})", spec.output.mode.name())
            } else {
                format!("wrote {} ({})", files.iter().map(|f| f.display().to_string()).collect::<Vec<_>>().join(" and "), spec.output.mode.name())
            }
        }));
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (mesh, result);
        Some(Err("this build of nyx was made without HDF5 support (rebuild with --features hdf5)".into()))
    }
}

/// Opens an HDF5 results file written by nyx.
pub fn open_h5(path: &Path) -> Result<(SimplexMesh, SolveResult, Option<String>), String> {
    #[cfg(feature = "hdf5")]
    {
        let (m, spec, r) = femcore::h5io::read_results(path)?;
        return Ok((m, r, spec));
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = path;
        Err("this build of nyx was made without HDF5 support (rebuild with --features hdf5)".into())
    }
}

pub fn spawn_open_h5(path: PathBuf, job: Job) {
    std::thread::spawn(move || {
        let r = open_h5(&path).map(|(mesh, result, _)| {
            let display = Display::new(&mesh, result.order);
            (Arc::new(mesh), Solved { result, display, state: SolveState::default(), output: None, analysis: None })
        });
        send(&job, Msg::H5(path.display().to_string(), r));
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unix_times_are_formatted_as_utc_stamps() {
        assert_eq!(stamp_from_unix(0), "19700101_000000");
        assert_eq!(stamp_from_unix(1_700_000_000), "20231114_221320");
        assert_eq!(stamp_from_unix(951_827_696), "20000229_123456");
    }

    #[test]
    fn the_local_timestamp_is_well_formed_and_costs_microseconds() {
        let t0 = std::time::Instant::now();
        let s = timestamp();
        let us = t0.elapsed().as_secs_f64() * 1e6;
        assert_eq!(s.len(), 15, "{s}");
        assert_eq!(s.as_bytes()[8], b'_');
        assert!(s.bytes().enumerate().all(|(i, b)| i == 8 || b.is_ascii_digit()), "{s}");
        eprintln!("timestamp {s} took {us:.1} us");
        // the old version started PowerShell: 170-1600 ms
        assert!(us < 5_000.0, "took {us} us");
    }
}