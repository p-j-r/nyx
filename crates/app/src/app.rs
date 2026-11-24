// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The application: the state of the window around one `nyx_engine::Engine`, and the window layout.
//!
//! `App` holds the engine (the project, the mesh it built, the results, the jobs), the mesh and results on show and the viewport. Every action
//! of the window is an engine command (`Engine::call`) or one of its job methods; the window keeps only what is about drawing. Anything slow
//! (a mesh build, a solve) is an engine job; only one runs at a time, and the window polls the engine each frame. A solve has a progress
//! handle: the toolbar shows its current stage and a Cancel button (Esc too), and the stage log it writes is part of the result.
//! `with_context` builds the state without a window, which is how the tests drive it.

use crate::colormap::Colormap;
use crate::examples;
use crate::state::*;
use crate::viewport::{Content, Viewport};
use eframe::egui::{self, Color32, RichText};
use nyx_connect::{Env, HostInfo, Status};
use nyx_engine::facade::geom::model::FeatureList;
use nyx_engine::facade::physics::SolveResult;
use nyx_engine::facade::Progress;
use nyx_engine::{Engine, MeshSettings, QueueOpts};
use serde_json::{json, Value};
use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::Arc;
use std::time::Instant;

/// Choices of the File > New dialog.
pub struct NewDialog {
    pub physics: String,
    pub three_d: bool,
}

/// Work waiting for the running job to finish (the latest request wins).
pub enum Queued {
    Open(PathBuf),
    Rebuild,
}

/// What a running engine job is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JobKind {
    Mesh,
    Solve,
    /// the two solves of a mesh convergence check
    Convergence,
    /// reading a results file, meshing its project again and checking the mesh (the draft workflow is made when it ends)
    Open,
}

/// A job of this window that waits in the engine's queue behind the running one (ADR 87).
pub struct WinJob {
    pub id: u64,
    pub kind: JobKind,
    pub label: String,
    pub cancellable: bool,
}

/// Result of running a `.nys` script on a background thread.
pub struct NysRunResult {
    pub engine: Box<Engine>,
    pub model: Option<Arc<nyx_engine::MeshedModel>>,
    pub results: Vec<(String, Arc<SolveResult>)>,
    pub printed: Vec<String>,
    pub value: Value,
    pub calls: u64,
    pub seconds: f64,
}

/// A message from a helper thread (the optional Python script) to the window.
pub enum Msg {
    Script(Result<String, String>),
    NysDone(Result<NysRunResult, String>),
    /// the trial mesh of an edit (the number in the tuple is the edit it was made for, `Feedback::generation`)
    Size(Result<nyx_engine::preview::SizePreview, String>),
    /// the polarization of a result, with the result it was computed from
    Polarization { physics: String, result: Arc<SolveResult>, analysis: Result<nyx_engine::analysis::ModeAnalysis, String> },
}

pub struct App {
    pub engine: Engine,
    pub doc: Doc,
    pub tab: Tab,
    pub vp: Viewport,
    pub results: Vec<ResultView>,
    pub current: Option<usize>,
    pub show_results: bool,
    /// the Reduce form of the Results tab (integrals, extremes and fractions of the field on show: the engine commands of ADR 90)
    pub reduce: crate::state::ReduceForm,
    pub tx: Sender<(u64, Msg)>,
    pub rx: Receiver<(u64, Msg)>,
    pub busy: Option<(String, Instant)>,
    /// the running engine job: its number, its progress/cancel handle and whether it can be cancelled
    pub job: Option<(u64, Arc<Progress>, bool)>,
    pub job_kind: Option<JobKind>,
    /// solves and convergence checks asked for while another job runs: they wait in the engine's queue in this order and become the running job
    /// one after the other
    pub waiting: Vec<WinJob>,
    /// the Queue panel is open
    pub show_queue: bool,
    /// scripted runs: cancel the running job at this time (`cancel_after SECONDS`)
    pub cancel_at: Option<Instant>,
    /// a mesh build or file open asked for while another job was running; it starts when that job is done
    pub queued: Option<Queued>,
    pub status: String,
    pub error: Option<String>,
    pub show_edges: bool,
    pub show_numbers: bool,
    /// the condition (index in the active physics) whose boundaries are being picked on the mesh
    pub picking: Option<usize>,
    /// the boundary numbers picked so far
    pub picked: Vec<i32>,
    pub project_path: Option<PathBuf>,
    pub geometry_dirty: bool,
    /// the geometry and mesh settings as being edited; they go into the engine when the mesh is built or a solve starts
    pub geom_edit: FeatureList,
    pub mesh_edit: MeshSettings,
    pub geom_ui: crate::geomui::State,
    /// why a condition (index in the active physics) does not resolve on the mesh
    pub rule_errors: HashMap<usize, String>,
    /// the boundary numbers each condition of the active physics picks on the mesh on show
    pub bc_labels: Vec<Vec<i32>>,
    /// region numbers with their names, for the Materials tab
    pub domain_names: Vec<(i32, String)>,
    /// the conditions or the mesh changed: resolve them again
    pub rules_dirty: bool,
    /// selectors that pick what was just clicked in the viewport
    pub suggestions: Vec<String>,
    pub bottom_open: bool,
    /// the File > New dialog while it is open
    pub new_dialog: Option<NewDialog>,
    pub script_output: Option<String>,
    pub nys: crate::sheets::ScriptSheet,
    pub script: VecDeque<String>,
    pub settle: u32,
    /// a script with `wait off` carries on while a job runs (to shoot a running job and the queue); by default it waits for the job
    pub script_free: bool,
    /// the sweep started from the window, if any (sweepui.rs)
    pub sweep: Option<crate::sweepui::Active>,
    /// colour-vision simulation applied to the screenshots of a script (`cvd protan|deutan|tritan|none`)
    pub shot_cvd: Option<crate::cvd::Kind>,
    pub shot_pending: Arc<AtomicBool>,
    pub shell: crate::shell::Shell,
    pub colormap: Colormap,
    /// the last convergence check: its answer, whether one is running, and the settings of the next
    pub convergence: Option<Value>,
    /// the validated state of a physics as the engine last said it (the plot asks every frame; the state changes only with a solve, a check or an edit)
    pub validated_cache: std::cell::RefCell<Option<(std::time::Instant, String, nyx_engine::validated::Validated)>>,
    /// the results files being written on engine jobs (job id, path): they do not make the window busy, a failure is shown when they end
    pub exports: HashMap<u64, String>,
    /// the title of the results file being opened on a job, `Some` from the request until the draft workflow is made
    pub opening: Option<String>,
    pub conv_factor: f64,
    pub conv_tol: f64,
    /// the `nyx-mcp` program the Connect menu registers with Claude (found at start, or picked by the user)
    pub mcp_server: Option<PathBuf>,
    /// the Connect menu: what each AI app's config says about nyx (installed apps only) and the last thing done
    pub connect_status: Vec<HostInfo>,
    pub connect_note: Option<String>,
    pub ctx: egui::Context,
    /// the problem and the size of the edits in the Geometry tab, and what a trial mesh said (`feedback.rs`)
    pub feedback: crate::feedback::Feedback,
    /// work that follows a solve, one piece per frame (`shell::Deferred`)
    pub deferred: VecDeque<crate::shell::Deferred>,
    /// form state of the Materials tab's on-demand download / import panel
    pub mat_fetch: crate::state::MaterialFetchForm,
}

impl App {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        crate::theme::install_fonts(&cc.egui_ctx);
        let mut app = Self::with_context(cc.egui_ctx.clone());
        let root = std::env::var_os("NYX_WORKSPACE").map(PathBuf::from).unwrap_or_else(nyx_engine::facade::store::workspace::default_root);
        app.attach_workspace(root);
        app
    }

    /// The application state on an egui context; no window is needed (tests use a default context). It opens with the first example.
    pub fn with_context(ctx: egui::Context) -> Self {
        crate::theme::apply(&ctx);
        let (tx, rx) = channel();
        let mut app = App {
            engine: Engine::new(),
            doc: Doc::default(),
            tab: Tab::Geometry,
            vp: Viewport::default(),
            results: vec![],
            current: None,
            show_results: false,
            reduce: Default::default(),
            tx,
            rx,
            busy: None,
            job: None,
            job_kind: None,
            waiting: vec![],
            show_queue: false,
            cancel_at: None,
            queued: None,
            status: "Ready".into(),
            error: None,
            show_edges: true,
            show_numbers: true,
            picking: None,
            picked: vec![],
            project_path: None,
            geometry_dirty: false,
            geom_edit: FeatureList::new(2),
            mesh_edit: MeshSettings::default(),
            feedback: Default::default(),
            deferred: VecDeque::new(),
            geom_ui: Default::default(),
            rule_errors: HashMap::new(),
            bc_labels: vec![],
            domain_names: vec![],
            rules_dirty: false,
            suggestions: vec![],
            bottom_open: true,
            new_dialog: None,
            script_output: None,
            nys: Default::default(),
            script: VecDeque::new(),
            settle: 6,
            script_free: false,
            sweep: None,
            shot_cvd: None,
            shot_pending: Arc::new(AtomicBool::new(false)),
            shell: crate::shell::Shell::detached(),
            colormap: Colormap::Viridis,
            convergence: None,
            validated_cache: std::cell::RefCell::new(None),
            exports: HashMap::new(),
            opening: None,
            conv_factor: 0.6,
            conv_tol: 1e-3,
            mcp_server: nyx_connect::find_server(),
            connect_status: vec![],
            connect_note: None,
            mat_fetch: Default::default(),
            ctx,
        };
        app.engine.set_owner("gui"); // a person's requests run before an agent's (ADR 87)
        app.refresh_connect();
        let first = examples::list().into_iter().next();
        let started = match first {
            Some(e) => app.engine.call("project_open", &json!({ "path": e.path.to_str().unwrap_or_default() })).map(|_| ()).map_err(|e| e.to_string()),
            None => examples::blank(&mut app.engine, "dfs", false),
        };
        if let Err(e) = started {
            app.set_error(e);
        }
        app.adopt();
        app.rebuild_mesh();
        app
    }

    // ------------------------------------------------------------- the engine

    /// Runs an engine command and reports its error in the status bar. `None` after an error.
    pub fn cmd(&mut self, name: &str, args: Value) -> Option<Value> {
        match self.engine.call(name, &args) {
            Ok(v) => Some(v),
            Err(e) => {
                self.set_error(format!("{name}: {e}"));
                None
            }
        }
    }

    /// The code of the active physics, if the project has one.
    pub fn active(&self) -> Option<String> {
        let p = self.engine.project();
        p.active.clone().or_else(|| p.physics.first().map(|s| s.physics.clone()))
    }

    /// The display name of a physics from its descriptor.
    pub fn physics_name(&self, id: &str) -> String {
        self.engine.registry().get(id).map(|p| p.name().to_string()).unwrap_or_else(|| id.to_string())
    }

    /// The window's copies of what is edited (geometry, mesh settings) are taken from the engine's project: after a project was opened or made.
    pub fn adopt(&mut self) {
        self.geom_edit = self.engine.project().geometry.clone();
        self.mesh_edit = self.engine.project().mesh.clone();
        self.feedback.reset();
        self.geometry_dirty = false;
        self.picking = None;
        self.picked.clear();
        self.suggestions.clear();
    }

    /// Hands the geometry and mesh settings being edited to the engine. False (and an error shown) when the engine refuses them.
    pub fn push_edits(&mut self) -> bool {
        // both edits are checked before either is applied, so a mesh edit that is refused cannot leave the geometry edit applied (a half-pushed
        // project); each command checks its own argument again when it runs
        if self.geom_edit != self.engine.project().geometry {
            if let Err(e) = self.geom_edit.validate() {
                self.set_error(format!("geometry: {e}"));
                return false;
            }
        }
        if self.mesh_edit != self.engine.project().mesh && !(self.mesh_edit.h > 0.0 && self.mesh_edit.h.is_finite()) {
            self.set_error("mesh: the element size must be a positive finite number");
            return false;
        }
        if self.geom_edit != self.engine.project().geometry {
            let g = serde_json::to_value(&self.geom_edit).unwrap_or(Value::Null);
            if self.cmd("set_geometry", json!({ "geometry": g })).is_none() {
                return false;
            }
            self.prune_stale_materials();
        }
        if self.mesh_edit != self.engine.project().mesh {
            let m = &self.mesh_edit;
            let args = json!({
                "h": m.h, "adapt_scale": m.adapt_scale, "physics_sizing": m.physics_sizing, "points_per_wavelength": m.points_per_wavelength,
                "structured_pml": m.structured_pml, "max_dev": m.max_dev, "curved": m.curved,
                "gap_cells": m.gap_cells, "z_coarsen": m.z_coarsen, "max_cells": m.max_cells,
                "dz_min": m.dz_min, "max_z_slabs": m.max_z_slabs, "z_growth_ratio": m.z_growth_ratio,
                "max_steiner_multiplier": m.max_steiner_multiplier, "strict_mesh": m.strict_mesh,
                "sweep": m.sweep,
            });
            if self.cmd("set_mesh", args).is_none() {
                return false;
            }
        }
        true
    }

    // ------------------------------------------------------------- jobs

    pub fn set_error(&mut self, e: impl Into<String>) {
        let e = e.into();
        self.status = "Error".into();
        self.error = Some(e);
    }

    /// A job was submitted to the engine. With nothing running it is the running job; otherwise it waits in the queue behind it.
    pub(crate) fn started(&mut self, id: u64, kind: JobKind, label: String, cancellable: bool) {
        let job = WinJob { id, kind, label, cancellable };
        if self.job.is_some() {
            self.status = format!("{} is queued behind the running job ({} waiting)", job.label, self.waiting.len() + 1);
            self.error = None;
            self.waiting.push(job);
        } else {
            self.make_front(job);
        }
    }

    /// `job` is the running job: the busy label, the progress handle and the run it will record.
    fn make_front(&mut self, job: WinJob) {
        let progress = self.engine.jobs().progress(job.id).unwrap_or_else(Progress::new);
        self.busy = Some((job.label, Instant::now()));
        self.job = Some((job.id, progress, job.cancellable));
        self.job_kind = Some(job.kind);
        self.error = None;
    }

    /// The next waiting job becomes the running one (the engine started it by itself when the lane freed).
    fn promote_waiting(&mut self) {
        if self.job.is_none() {
            // a job held in the queue does not start, so it is not the running one until it is released
            let held = |id: u64| self.engine.jobs().info(id).is_some_and(|i| i.held);
            if let Some(i) = self.waiting.iter().position(|w| !held(w.id)) {
                let next = self.waiting.remove(i);
                self.make_front(next);
            }
        }
    }

    /// Cancels a job that waits in the queue (it never runs).
    pub fn cancel_waiting(&mut self, id: u64) {
        if self.engine.jobs().cancel(id) {
            self.waiting.retain(|w| w.id != id);
            self.status = "A queued job was cancelled before it started".into();
        }
    }

    /// Asks the running solve to stop (it does so at its next check; a factorisation in progress has to finish first).
    /// Solves the active physics again on a finer mesh on its own thread and compares (the engine command `convergence_check`); the answer
    /// arrives as a message and is shown in the Study tab.
    pub fn check_convergence(&mut self) {
        let Some(active) = self.active() else { return };
        if self.refuse_unless_queueable() || !self.push_edits() {
            return;
        }
        // an engine job: it shows in the job list, reports its two solves as stages and can be cancelled; behind a running job it waits in the
        // queue and replaces a check that has not started
        let opts = QueueOpts { supersede: true, ..Default::default() };
        match self.engine.start_convergence_with(&active, self.conv_factor, Some(self.conv_tol), &opts) {
            Ok(id) => {
                self.waiting.retain(|w| w.kind != JobKind::Convergence);
                self.convergence = None;
                self.started(id, JobKind::Convergence, "Checking convergence".into(), true);
            }
            Err(e) => self.set_error(format!("convergence: {e}")),
        }
    }

    /// Whether a Compute can be asked for now: not while a results file is opened (it replaces the project), and in a workspace only when nothing
    /// runs (each Compute there is one run record and the page is locked while it runs).
    pub fn can_queue_compute(&self) -> bool {
        self.job_kind != Some(JobKind::Open) && (!self.shell.active() || self.busy.is_none())
    }

    /// True (and a message in the status bar) when a solve or a check cannot be submitted now: while a results file is being opened, which
    /// replaces the project the request would have been made for. Behind any other job the request waits in the queue.
    fn refuse_unless_queueable(&mut self) -> bool {
        if self.job_kind == Some(JobKind::Open) {
            self.status = "Opening a results file: it replaces the project; wait for it to finish or cancel it".into();
            return true;
        }
        false
    }

    /// True while a mesh convergence check runs.
    pub fn convergence_running(&self) -> bool {
        self.job_kind == Some(JobKind::Convergence)
    }

    /// Why the result at `i` is not up to date, in words, or `None` when it is current (productivity plan P1, ADR 84). It is stale when the engine
    /// says what it was solved from has changed (the problem key of its physics, so an edit of another physics does not count), when a newer result
    /// of the physics replaced it, or when the Geometry or Mesh tab holds edits that were not applied to the engine yet.
    pub fn stale_reason(&self, i: usize) -> Option<&'static str> {
        let rv = self.results.get(i)?;
        if let Some(id) = &rv.physics {
            match self.engine.result(id) {
                None => return Some("this result is no longer the engine's"),
                Some((r, stale)) => {
                    if !Arc::ptr_eq(r, &rv.result) {
                        return Some("a newer result of this physics replaced it");
                    }
                    if stale {
                        return Some("the model changed since this solve");
                    }
                }
            }
        }
        if self.geom_edit != self.engine.project().geometry || self.mesh_edit != self.engine.project().mesh {
            return Some("Geometry or Mesh edits are not applied yet");
        }
        None
    }

    /// Diagnostics only (productivity plan P0): when `NYX_FRAME_LOG` names a file ending in `.csv`, appends one row per frame: wall time of the
    /// frame since the previous row, our layout and drawing time, the viewport time, and what the engine is doing (job kind, stage, seconds in
    /// the stage and in the job). Nothing is written, and nothing is read from the engine, when the variable is unset or is not such a path.
    fn frame_csv(&self, ui_ms: f64, viewport_ms: f64) {
        use std::io::Write;
        use std::sync::{Mutex, OnceLock};
        static LOG: OnceLock<Option<Mutex<(std::fs::File, Instant)>>> = OnceLock::new();
        let log = LOG.get_or_init(|| {
            let path = std::env::var("NYX_FRAME_LOG").ok().filter(|p| p.to_lowercase().ends_with(".csv"))?;
            let mut f = std::fs::File::create(path).ok()?;
            let _ = writeln!(f, "t_s,frame_dt_ms,ui_ms,viewport_ms,busy,job_kind,stage,stage_s,job_s,convergence,show_results,spans");
            Some(Mutex::new((f, Instant::now())))
        });
        let Some(log) = log else { return };
        let Ok(mut g) = log.lock() else { return };
        let now = Instant::now();
        let dt = now.duration_since(g.1).as_secs_f64() * 1e3;
        g.1 = now;
        let kind = match self.job_kind {
            Some(JobKind::Mesh) => "mesh",
            Some(JobKind::Solve) => "solve",
            Some(JobKind::Convergence) => "convergence",
            Some(JobKind::Open) => "open",
            None => "",
        };
        let (stage, stage_s) = match &self.job {
            Some((_, p, _)) => {
                let st = p.status();
                (st.stage.replace([',', '\n'], " "), st.stage_seconds)
            }
            None => (String::new(), 0.0),
        };
        let job_s = self.busy.as_ref().map(|(_, t)| t.elapsed().as_secs_f64()).unwrap_or(0.0);
        let t = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs_f64() % 100000.0).unwrap_or(0.0);
        let _ = writeln!(
            g.0,
            "{t:.3},{dt:.2},{ui_ms:.2},{viewport_ms:.2},{},{kind},{stage},{stage_s:.2},{job_s:.2},{},{},{}",
            self.busy.is_some() as u8,
            self.convergence_running() as u8,
            self.show_results as u8,
            crate::spans::take()
        );
        let _ = g.0.flush();
    }

    pub fn cancel_job(&mut self) {
        if let Some((id, _, true)) = &self.job {
            self.engine.jobs().cancel(*id);
            self.status = "Cancel requested: the job stops at its next check (a factorisation in progress has to finish first)".into();
        }
    }

    /// True (and a message in the status bar) when a job is running and a new one cannot start.
    fn refuse_while_busy(&mut self) -> bool {
        if let Some((what, _)) = &self.busy {
            self.status = format!("Busy ({what}): wait for it to finish{}", if self.job.as_ref().is_some_and(|j| j.2) { " or cancel it" } else { "" });
            return true;
        }
        false
    }

    pub fn rebuild_mesh(&mut self) {
        if self.busy.is_some() {
            self.queued = Some(Queued::Rebuild);
            return;
        }
        if !self.push_edits() {
            return;
        }
        match self.engine.start_mesh() {
            Ok(id) => self.started(id, JobKind::Mesh, "Building mesh".into(), true),
            Err(e) => self.set_error(format!("Geometry: {e}")),
        }
    }

    /// Solves run in a child process (this program started with `--job`, see `nyx_engine::worker`), so a solve that runs out of memory or crashes
    /// fails its job and not the window, and a cancel stops it at once (productivity plan P4, ADR 86). Called by `main`, not by `App::new`: a test
    /// binary must not be its own worker.
    pub fn solve_in_child_processes(&mut self) {
        self.engine.set_self_as_worker(true);
        // what solves cost on this machine is learned in the workspace folder, for the memory estimates of later solves (ADR 86)
        let root = std::env::var_os("NYX_WORKSPACE").map(PathBuf::from).unwrap_or_else(nyx_engine::facade::store::workspace::default_root);
        self.engine.set_history_file(&root.join("solve_history.jsonl"));
        // the queue of the workspace is shared with the MCP server and `nyx queue`; closing the window leaves its items running (Q3, ADR 95)
        self.engine.set_shared_queue(&root, true);
    }

    /// Loads a `.nys` script into the built-in script sheet without executing it (safety rule against untrusted files).
    pub fn open_nys_file(&mut self, path: PathBuf) {
        match std::fs::read_to_string(&path) {
            Ok(text) => {
                self.nys.code = text.trim_start_matches('\u{feff}').to_string();
                self.nys.path = Some(path.clone());
                self.nys.check_msg = None;
                self.nys.output = None;
                self.shell.sheet = crate::shell::Sheet::None;
                self.nys.open = true;
                self.status = format!("Loaded script {}", path.display());
            }
            Err(e) => self.set_error(format!("Could not read script {}: {e}", path.display())),
        }
    }

    pub fn open_path(&mut self, path: PathBuf) {
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
        if ext == "json" {
            if self.busy.is_some() {
                self.queued = Some(Queued::Open(path));
                return;
            }
            if self.shell.active() {
                self.open_received(path);
            } else {
                self.load_project(&path);
            }
            return;
        }
        if ext == "h5" {
            if self.busy.is_some() {
                self.queued = Some(Queued::Open(path));
                return;
            }
            let title = path.file_stem().and_then(|s| s.to_str()).unwrap_or("Results").to_string();
            self.new_draft(&title, crate::shell::Start::Results(path));
            return;
        }
        if ext == "nyx" {
            self.open_received(path);
            return;
        }
        if ext == "nys" {
            self.open_nys_file(path);
            return;
        }
        self.set_error(format!("Don't know how to open '.{ext}' files (supported: .nyx.json, .nyx, .nys, .h5)"));
    }

    /// Solves the active physics (with the physics it takes fields from), meshing first when the geometry or mesh settings changed.
    pub fn compute(&mut self) {
        // in a workspace every Compute is one run record and the workflow page is locked while it runs, so the shell keeps one solve at a time;
        // without one a Compute behind a running job waits in the queue
        if self.refuse_unless_queueable() || (self.shell.active() && self.refuse_while_busy()) {
            return;
        }
        self.sync_rules();
        if let Some(((i, e), _)) = self.rule_errors.iter().map(|(k, v)| ((*k, v.clone()), ())).min_by_key(|((k, _), _)| *k) {
            self.set_error(format!("Condition {}: {e}", i + 1));
            return;
        }
        if !self.push_edits() {
            return;
        }
        let Some(active) = self.active() else {
            self.set_error("The project has no physics: add one in the Physics tab.");
            return;
        };
        // behind a running job the solve waits in the engine's queue; a Compute that has not started yet is replaced by this one
        let opts = QueueOpts { supersede: true, ..Default::default() };
        match self.engine.start_solve_with(None, &opts) {
            Ok(id) => {
                self.waiting.retain(|w| w.kind != JobKind::Solve);
                let label = format!("Solving {}", self.physics_name(&active));
                self.started(id, JobKind::Solve, label, true);
            }
            Err(e) => self.set_error(e.to_string()),
        }
    }

    fn poll(&mut self) {
        self.poll_sweep();
        while let Ok((generation, msg)) = self.rx.try_recv() {
            match msg {
                Msg::Script(Ok(t)) => self.script_output = Some(t),
                Msg::Script(Err(e)) => self.set_error(e),
                Msg::NysDone(res) => self.take_nys_done(res),
                Msg::Size(out) => self.take_size(generation, out),
                Msg::Polarization { physics, result, analysis } => self.take_polarization(&physics, &result, analysis),
            }
        }
        let finished = crate::spans::span("engine_poll", || self.engine.poll());
        let n_finished = finished.len();
        for (id, outcome) in finished {
            // a results file written on a job: it never held the window, only a failure is shown
            if let Some(path) = self.exports.remove(&id) {
                if let Err(e) = outcome {
                    self.error = Some(format!("Results file {path}: {e}"));
                }
                continue;
            }
            // a job that waited in the queue and is already over (the engine started it when the one before it ended) is the running one now; a
            // job that is none of ours (replaced by a newer request, or cancelled in the queue) has nothing to show
            if self.job.as_ref().map(|j| j.0) != Some(id) {
                match self.waiting.iter().position(|w| w.id == id) {
                    Some(i) if self.job.is_none() => {
                        let w = self.waiting.remove(i);
                        self.make_front(w);
                    }
                    _ => continue,
                }
            }
            // only the job that is running clears the busy state
            let started = if self.job.as_ref().map(|j| j.0) == Some(id) {
                self.job = None;
                self.busy.take().map(|b| b.1)
            } else {
                None
            };
            self.forget_validated();
            let kind = self.job_kind.take();
            let took = started.map(|s| s.elapsed().as_secs_f64()).unwrap_or(0.0);
            match (kind, outcome) {
                (Some(JobKind::Mesh), Ok(_)) => {
                    self.take_model();
                    let (c, f) = self.doc.mesh.as_ref().map(|m| (m.n_cells(), m.n_facets())).unwrap_or((0, 0));
                    self.status = format!("Mesh built: {c} cells, {f} boundary facets ({took:.2} s)");
                    self.geometry_dirty = false;
                }
                (Some(JobKind::Mesh), Err(nyx_engine::facade::physics::Error::Cancelled)) => {
                    self.error = None;
                    self.status = format!("Mesh build cancelled after {took:.1} s; the previous mesh is unchanged");
                }
                (Some(JobKind::Mesh), Err(e)) => self.set_error(format!("Geometry: {e}")),
                (Some(JobKind::Convergence), Ok(v)) => {
                    self.convergence = Some(v);
                    self.note_check_in_run();
                    self.status = format!("Convergence checked ({took:.1} s)");
                }
                (Some(JobKind::Convergence), Err(nyx_engine::facade::physics::Error::Cancelled)) => {
                    self.error = None;
                    self.status = format!("Convergence check cancelled after {took:.1} s");
                }
                (Some(JobKind::Convergence), Err(e)) => self.set_error(format!("convergence: {e}")),
                (Some(JobKind::Open), Ok(summary)) => self.finish_open(&summary),
                (Some(JobKind::Open), Err(nyx_engine::facade::physics::Error::Cancelled)) => {
                    self.opening = None;
                    self.status = "Opening the results file cancelled; nothing was changed".into();
                }
                (Some(JobKind::Open), Err(e)) => {
                    self.opening = None;
                    self.set_error(e.to_string());
                }
                (Some(JobKind::Solve), Ok(summary)) => self.take_solved(&summary),
                (Some(JobKind::Solve), Err(nyx_engine::facade::physics::Error::Cancelled)) => {
                    self.error = None;
                    self.record_end(nyx_engine::facade::store::runs::RunStatus::Cancelled, None);
                    self.status = format!("Solve cancelled after {took:.1} s; the previous results are unchanged");
                }
                (Some(JobKind::Solve), Err(e)) => {
                    self.record_end(nyx_engine::facade::store::runs::RunStatus::Failed, Some(e.to_string()));
                    self.set_error(e.to_string())
                }
                (None, _) => {}
            }
            self.promote_waiting();
        }
        self.promote_waiting();
        // one piece of the work that follows a solve, never in the frame that took a result (that frame has enough to do)
        if n_finished == 0 {
            self.run_deferred();
        }
    }

    /// Takes the mesh the engine has now (a new mesh clears the results of the old one).
    fn take_model(&mut self) {
        let Some(m) = self.engine.model().cloned() else { return };
        if self.doc.model.as_ref().is_some_and(|d| Arc::ptr_eq(d, &m)) {
            return;
        }
        self.doc.set_model(m, "generated geometry");
        self.after_mesh_change();
    }

    /// A solve came back: the mesh it used, and one result view per physics it solved.
    pub(crate) fn take_solved(&mut self, summary: &Value) {
        crate::spans::span("take_model", || self.take_model());
        let Some(mesh) = self.doc.mesh.clone() else { return };
        // a job that finished after the project changed again, or after a newer job, can come back with results the engine did not keep (`stored`
        // false): they are not shown and not recorded as the run on screen
        let all: Vec<&Value> = summary["results"].as_array().into_iter().flatten().collect();
        let solved: Vec<String> = all.iter().filter(|r| r["stored"] != json!(false)).filter_map(|r| r["physics"].as_str().map(str::to_string)).collect();
        if solved.is_empty() && !all.is_empty() {
            self.shell.run = None;
            self.status = "A solve finished after the model changed again; its result was not kept (solve again)".into();
            return;
        }
        let active = self.active();
        let mut last_secs = 0.0;
        let mut live: Option<Arc<SolveResult>> = None;
        for id in &solved {
            let Some((result, _)) = self.engine.result(id) else { continue };
            let result = result.clone();
            let name = self.physics_name(id);
            last_secs = result.seconds;
            let rv = crate::spans::span("result_view", || ResultView::new(name, Some(id.clone()), result.clone(), &mesh));
            // the polarization is computed beside the frame and arrives through the channel (`feedback.rs`): the cost grows with the cells
            // (8 to 9 ms at 9,500 cells in 3D, 41 ms for the 8,000 cells of the silicon-on-sapphire example) and the frame has 33 ms
            self.start_polarization(id);
            if active.as_deref() == Some(id.as_str()) {
                live = Some(result);
            }
            self.push_result(rv);
        }
        if let Some(r) = live {
            // the run record is written on a later frame (`run_deferred`): file writes, 9 to 17 ms, are not part of showing the result
            self.shell.live_run = None;
            if let Some(rc) = self.shell.run.take() {
                self.defer(crate::shell::Deferred::Record { rc, result: r });
            }
        }
        // the result on show is the active physics' one
        if let Some(i) = active.and_then(|a| self.results.iter().position(|r| r.physics.as_deref() == Some(a.as_str()))) {
            self.current = Some(i);
        }
        let name = self.current.map(|i| self.results[i].name.clone()).unwrap_or_default();
        self.status = format!("{name}: solved in {last_secs:.2} s");
        self.tab = Tab::Results;
    }

    /// Starts the work that was asked for while a job was running, once nothing is running.
    fn run_queued(&mut self) {
        if self.busy.is_some() {
            return;
        }
        match self.queued.take() {
            Some(Queued::Open(path)) => self.open_path(path),
            Some(Queued::Rebuild) => self.rebuild_mesh(),
            None => {}
        }
    }

    fn push_result(&mut self, rv: ResultView) {
        if let Some(i) = self.results.iter().position(|r| r.name == rv.name) {
            self.results[i] = rv;
            self.current = Some(i);
        } else {
            self.results.push(rv);
            self.current = Some(self.results.len() - 1);
        }
        self.show_results = true;
        self.vp.fit_token += 1;
        // in the window of the shell a new result opens its tab beside the field (ADR 124)
        if self.shell.active() {
            if let Some(w) = self.shell.current_mut() {
                w.tab = Tab::Results;
                w.sub = crate::shell::Sub::Main;
            }
            self.tab = Tab::Results;
        }
    }

    fn after_mesh_change(&mut self) {
        self.results.clear();
        self.current = None;
        self.show_results = false;
        self.picking = None;
        self.picked.clear();
        self.vp.fit_token += 1;
        self.sync_rules();
    }

    /// Resolves the selectors of the conditions of the active physics to the boundary numbers of the mesh on show (to draw them and to
    /// report a selector that selects nothing), and reads the names of the regions.
    pub fn sync_rules(&mut self) {
        self.rules_dirty = false;
        self.rule_errors.clear();
        self.bc_labels.clear();
        let Some(model) = self.doc.model.clone() else { return };
        self.domain_names = model.regions();
        let Some(active) = self.active() else { return };
        let Some(spec) = self.engine.project().spec(&active).cloned() else { return };
        use nyx_engine::facade::physics::LabelResolver;
        for (i, bc) in spec.bcs.iter().enumerate() {
            let applies_to_domain = self.engine.registry().get(&active).is_some_and(|p| p.bc_kinds().iter().any(|k| k.kind == bc.kind && k.applies == nyx_engine::facade::physics::Applies::Domain));
            let r = if applies_to_domain { LabelResolver::domains(model.as_ref(), &bc.on) } else { LabelResolver::boundaries(model.as_ref(), &bc.on) };
            match r {
                Ok(l) => self.bc_labels.push(if applies_to_domain { vec![] } else { l }),
                Err(e) => {
                    self.bc_labels.push(vec![]);
                    self.rule_errors.insert(i, e.to_string());
                }
            }
        }
    }

    // ------------------------------------------------------------- files

    pub fn save_project(&mut self, as_new: bool) {
        let path = if as_new || self.project_path.is_none() {
            rfd::FileDialog::new().add_filter("nyx project (JSON)", &["json"]).set_file_name("project.nyx.json").save_file()
        } else {
            self.project_path.clone()
        };
        let Some(path) = path else { return };
        if !self.push_edits() {
            return;
        }
        if self.cmd("project_save", json!({ "path": path.to_str().unwrap_or_default() })).is_some() {
            self.status = format!("Saved {}", path.display());
            self.project_path = Some(path);
        }
    }


    /// File > Export as script: the project as a script that rebuilds it (`project_script`, ADR 94), written where an agent finds it by name, the
    /// MCP server's `scripts` folder. The window's buffered edits are pushed first, as Save does, so the script is what the person sees.
    pub fn export_script(&mut self) {
        if !self.push_edits() {
            return;
        }
        let Some(r) = self.cmd("project_script", json!({})) else { return };
        let (Some(text), form) = (r["script"].as_str(), r["form"].as_str().unwrap_or("commands")) else { return };
        // a name an agent can give to `script_run path=...`: letters, digits, - and _
        let stem: String = self.engine.project().title.trim().chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' { c.to_ascii_lowercase() } else { '_' }).collect();
        let stem = if stem.trim_matches('_').is_empty() { "project".to_string() } else { stem.trim_matches('_').chars().take(60).collect() };
        let dir = nyx_script::mcp_dir().join("scripts");
        let _ = std::fs::create_dir_all(&dir);
        let Some(path) = rfd::FileDialog::new().add_filter("nyx script", &["nys"]).set_directory(&dir).set_file_name(format!("{stem}.nys")).save_file() else { return };
        match std::fs::write(&path, text) {
            Ok(()) => {
                let how = if form == "commands" { "" } else { " (one block, not command by command)" };
                self.status = format!("Exported {}{how}. An agent runs it with script_run path \"{stem}\" and allow_destructive.", path.display());
            }
            Err(e) => self.set_error(format!("Could not write the file: {e}")),
        }
    }

    pub fn load_project(&mut self, path: &std::path::Path) {
        if self.busy.is_some() {
            self.queued = Some(Queued::Open(path.to_path_buf()));
            return;
        }
        if self.cmd("project_open", json!({ "path": path.to_str().unwrap_or_default() })).is_some() {
            self.adopt();
            self.project_path = Some(path.to_path_buf());
            self.rebuild_mesh();
            self.status = format!("Loaded {}", path.display());
        }
    }

    pub fn load_example(&mut self, i: usize) {
        let Some(e) = examples::list().into_iter().nth(i) else {
            self.set_error("There are no examples: the `examples` folder was not found next to the program (set NYX_EXAMPLES).");
            return;
        };
        if self.cmd("project_open", json!({ "path": e.path.to_str().unwrap_or_default() })).is_some() {
            self.adopt();
            self.project_path = None;
            self.rebuild_mesh();
            self.tab = Tab::Geometry;
            self.status = format!("Loaded example: {}", e.title);
        }
    }

    /// CSV of the current solution: coordinates of every dof and all fields.
    pub fn result_csv(&self) -> Option<String> {
        let (i, mesh) = (self.current?, self.doc.mesh.clone()?);
        let rv = &self.results[i];
        let sol = rv.solution();
        let order = rv.result.order as usize;
        let co = mesh.dof_coordinates(order);
        let d = mesh.dim;
        let mut s = String::new();
        s += &format!("% {} | {} | {}\n", rv.name, sol.label, self.doc.source);
        s += "% x,y";
        if d == 3 {
            s += ",z";
        }
        for f in &sol.fields {
            s += &format!(",{} [{}]", f.name, f.unit);
        }
        s += "\n";
        for n in 0..mesh.n_dofs(order) {
            let mut row: Vec<String> = (0..d).map(|k| format!("{:.12e}", co[n * d + k])).collect();
            for f in &sol.fields {
                row.push(format!("{:.12e}", f.values[n]));
            }
            s += &row.join(",");
            s += "\n";
        }
        Some(s)
    }

    pub fn export_csv(&mut self) {
        let Some(s) = self.result_csv() else {
            self.set_error("Nothing to export: solve something first.");
            return;
        };
        let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("results.csv").save_file() else { return };
        match std::fs::write(&path, s) {
            Ok(()) => self.status = format!("Exported {}", path.display()),
            Err(e) => self.set_error(format!("Could not write the file: {e}")),
        }
    }

    /// Writes the polarization analysis of the result on show: `csv` (one row per mode) or `svg` (the Poincare sphere).
    pub fn export_analysis(&mut self, kind: &str) {
        let Some(a) = self.current.and_then(|i| self.results[i].analysis.as_ref()) else {
            self.set_error("This result has no polarization analysis.");
            return;
        };
        let (text, filter, name) = match kind {
            "svg" => (a.sphere_svg(480.0), "SVG", "poincare_sphere.svg"),
            _ => (a.csv(), "CSV", "polarization.csv"),
        };
        let Some(path) = rfd::FileDialog::new().add_filter(filter, &[kind]).set_file_name(name).save_file() else { return };
        match std::fs::write(&path, text) {
            Ok(()) => self.status = format!("Exported {}", path.display()),
            Err(e) => self.set_error(format!("Could not write the file: {e}")),
        }
    }

    pub fn run_python_script(&mut self) {
        let Some(script) = rfd::FileDialog::new().add_filter("Python script", &["py"]).pick_file() else { return };
        if !self.push_edits() {
            return;
        }
        let json = self.engine.project().to_json();
        if self.refuse_while_busy() {
            return;
        }
        self.status = format!("Running {}", script.file_name().and_then(|n| n.to_str()).unwrap_or("script"));
        crate::pyscript::spawn(script, json, self.result_csv(), self.tx.clone(), self.ctx.clone());
    }

    pub fn fill_nys_from_model(&mut self) {
        let _ = self.push_edits();
        if let Ok(v) = self.engine.call("project_script", &json!({})) {
            if let Some(s) = v["script"].as_str() {
                self.nys.code = s.to_string();
                self.nys.check_msg = None;
            }
        }
    }

    pub fn open_nys_script(&mut self) {
        self.shell.sheet = crate::shell::Sheet::None;
        self.nys.open = true;
        if self.nys.code.trim().is_empty() {
            self.fill_nys_from_model();
        }
    }

    /// The key and the button of the script: opens the dock, or closes it when it is open.
    pub fn toggle_nys_script(&mut self) {
        if self.nys.open {
            self.nys.open = false;
        } else {
            self.open_nys_script();
        }
    }

    fn nys_root(&self) -> PathBuf {
        if let crate::shell::Page::Workflow(id) = &self.shell.page {
            if let Some(d) = self.shell.open.get(id).and_then(|w| w.dir.as_ref()) {
                return d.root.clone();
            }
        }
        if let Some(p) = self.nys.path.as_ref().and_then(|p| p.parent()) {
            if !p.as_os_str().is_empty() {
                return p.to_path_buf();
            }
        }
        if let Some(ws) = &self.shell.ws {
            return ws.root.clone();
        }
        std::env::temp_dir()
    }

    pub fn check_nys_script(&mut self) {
        let _ = self.push_edits();
        let mut eng = Engine::new();
        if let Ok(spec) = self.spec_json() {
            let _ = eng.call("project_load", &json!({ "json": spec }));
        }
        let name = self.nys.path.as_ref().and_then(|p| p.file_name()).and_then(|n| n.to_str()).unwrap_or("script.nys");
        let code = self.nys.code.trim_start_matches('\u{feff}');
        let opts = nyx_script::Options::new(self.nys_root());
        match nyx_script::host::check(&mut eng, name, code, &nyx_script::Limits::default(), &opts) {
            Ok(()) => {
                self.nys.check_msg = Some(Ok("Static check passed: syntax and nyx commands are valid.".into()));
            }
            Err(e) => {
                self.nys.check_msg = Some(Err(e));
            }
        }
    }

    pub fn run_nys_script(&mut self) {
        if self.nys.running || self.refuse_while_busy() {
            return;
        }
        if nyx_script::find_worker().is_none() {
            let msg = "nyx-script-worker was not found next to nyx (build it with `cargo build --release -p nyx-script-worker` or set NYX_SCRIPT_WORKER).".to_string();
            self.nys.output = Some(Err(msg.clone()));
            self.set_error(msg);
            return;
        }
        let _ = self.push_edits();
        let spec = self.spec_json().unwrap_or_else(|_| self.engine.project().to_json());
        let code = self.nys.code.trim_start_matches('\u{feff}').to_string();
        let name = self.nys.path.as_ref().and_then(|p| p.file_name()).and_then(|n| n.to_str()).unwrap_or("script.nys").to_string();
        let root = self.nys_root();
        self.nys.running = true;
        self.nys.check_msg = None;
        self.nys.output = None;
        self.status = format!("Running {name}...");
        let (tx, ctx) = (self.tx.clone(), self.ctx.clone());
        std::thread::spawn(move || {
            let started = Instant::now();
            let mut eng = Engine::new();
            eng.set_self_as_worker(false);
            eng.set_owner("gui");
            let _ = eng.call("project_load", &json!({ "json": spec }));
            let opts = nyx_script::Options::new(root);
            let mut pols: Vec<(String, Arc<SolveResult>, nyx_engine::analysis::ModeAnalysis)> = vec![];
            let res = nyx_script::host::run(&mut eng, &name, &code, &nyx_script::Limits::default(), &opts).map(|out| {
                let model = eng.model().cloned();
                let mut results = vec![];
                for s in &eng.project().physics {
                    if let Some((r, _)) = eng.result(&s.physics) {
                        let r = r.clone();
                        if let Ok(a) = eng.polarization(&s.physics) {
                            pols.push((s.physics.clone(), r.clone(), a));
                        }
                        results.push((s.physics.clone(), r));
                    }
                }
                NysRunResult {
                    engine: Box::new(eng),
                    model,
                    results,
                    printed: out.printed,
                    value: out.value,
                    calls: out.calls,
                    seconds: started.elapsed().as_secs_f64(),
                }
            });
            let _ = tx.send((0, Msg::NysDone(res)));
            for (physics, result, analysis) in pols {
                let _ = tx.send((0, Msg::Polarization { physics, result, analysis: Ok(analysis) }));
            }
            ctx.request_repaint();
        });
    }

    fn take_nys_done(&mut self, res: Result<NysRunResult, String>) {
        self.nys.running = false;
        match res {
            Err(e) => {
                self.nys.output = Some(Err(e.clone()));
                self.set_error(e);
            }
            Ok(out) => {
                let project_json = out.engine.project().to_json();
                self.engine = *out.engine;
                self.adopt();
                self.forget_validated();
                if let Some(m) = out.model {
                    self.doc.set_model(m, "script");
                    self.after_mesh_change();
                } else {
                    self.rebuild_mesh();
                }
                let active = self.active();
                let mut active_res: Option<Arc<SolveResult>> = None;
                if let Some(mesh) = self.doc.mesh.clone() {
                    for (phys, r) in &out.results {
                        let name = self.physics_name(phys);
                        let rv = ResultView::new(name, Some(phys.clone()), r.clone(), &mesh);
                        if active.as_deref() == Some(phys.as_str()) {
                            active_res = Some(r.clone());
                        }
                        self.push_result(rv);
                    }
                }
                if let Some(i) = active.as_ref().and_then(|a| self.results.iter().position(|r| r.physics.as_deref() == Some(a.as_str()))) {
                    self.current = Some(i);
                }
                if let crate::shell::Page::Workflow(wid) = self.shell.page.clone() {
                    if self.ensure_dir(&wid) {
                        self.flush();
                        if let Some(r) = active_res.or_else(|| out.results.first().map(|(_, r)| r.clone())) {
                            let eid = self.shell.open.get(&wid).map(|w| w.exp.clone()).unwrap_or_else(|| "e1".into());
                            let spec = self.spec_json().unwrap_or(project_json);
                            let params = crate::shell::params_of(&spec);
                            let rc = crate::shell::RunCtx {
                                wf_id: wid,
                                exp: eid,
                                spec_json: spec,
                                params,
                                started: Instant::now(),
                            };
                            self.shell.live_run = self.record_ok(rc, &r);
                        }
                    }
                }
                let mut lines = out.printed;
                if !out.value.is_null() {
                    lines.push(serde_json::to_string_pretty(&out.value).unwrap_or_default());
                }
                lines.push(format!("{} engine commands in {:.2} s", out.calls, out.seconds));
                self.nys.output = Some(Ok(lines.join("\n")));
                self.status = format!("Script finished: {} engine commands in {:.2} s", out.calls, out.seconds);
            }
        }
    }

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new()
            .add_filter("nyx files", &["nyx", "nys", "json", "h5"])
            .add_filter("nyx projects (JSON)", &["json"])
            .add_filter("nyx workflow archives", &["nyx"])
            .add_filter("nyx scripts", &["nys"])
            .add_filter("nyx results (HDF5)", &["h5"])
            .pick_file()
        {
            self.open_path(p);
        }
    }

    /// Replaces the project by a new blank one (nothing is kept of the old geometry, conditions or results).
    pub fn new_project(&mut self, physics: &str, three_d: bool, what: &str) {
        match examples::blank(&mut self.engine, physics, three_d) {
            Ok(()) => {
                self.adopt();
                self.project_path = None;
                self.rebuild_mesh();
                self.tab = Tab::Geometry;
                self.status = format!("New simulation: {what}. Set the geometry, materials and conditions, then Compute (F5).");
            }
            Err(e) => self.set_error(e),
        }
    }


    // ------------------------------------------------------------- layout

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("New... (Ctrl+N)").clicked() {
                    ui.close();
                    self.new_dialog = Some(NewDialog { physics: self.active().unwrap_or_else(|| "dfs".into()), three_d: false });
                }
                if ui.button("Open... (Ctrl+O)").clicked() {
                    ui.close();
                    self.open_dialog();
                }
                if ui.button("Save project (Ctrl+S)").clicked() {
                    ui.close();
                    self.save_project(false);
                }
                if ui.button("Save project as...").clicked() {
                    ui.close();
                    self.save_project(true);
                }
                ui.separator();
                if ui.button("Export as script...").clicked() {
                    ui.close();
                    self.export_script();
                }
                if ui.button("Export current result as CSV...").clicked() {
                    ui.close();
                    self.export_csv();
                }
                ui.separator();
                if ui.button("Quit").clicked() {
                    ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });
            ui.menu_button("Examples", |ui| {
                let list = examples::list();
                if list.is_empty() {
                    ui.weak("No examples folder found (set NYX_EXAMPLES).");
                }
                for (i, e) in list.iter().enumerate() {
                    if ui.button(&e.title).clicked() {
                        ui.close();
                        self.load_example(i);
                    }
                }
            });
            let connect_menu = ui.menu_button("Connect", |ui| self.connect_menu(ui));
            if connect_menu.response.clicked() {
                self.refresh_connect();
            }
            ui.menu_button("Python (optional)", |ui| {
                ui.label("Run your own Python script on the current project and result.");
                ui.label("It receives the project JSON and the result CSV as arguments.");
                if ui.button("Run script...").clicked() {
                    ui.close();
                    self.run_python_script();
                }
            });
            ui.menu_button("View", |ui| {
                ui.checkbox(&mut self.show_edges, "Mesh edges");
                ui.checkbox(&mut self.show_numbers, "Boundary numbers (mesh view)");
                ui.checkbox(&mut self.bottom_open, "Results panel");
                if ui.button("Reset view").clicked() {
                    ui.close();
                    self.vp.fit_token += 1;
                }
            });
            ui.menu_button("Help", |ui| {
                ui.label("F5 computes the active physics.");
                ui.label("Pick boundaries on the mesh: Physics tab > a condition > 'Pick on mesh'.");
                ui.label("The numerics are Rust (femcore); no Python is needed.");
            });
        });
    }

    /// Re-reads what Claude Code's and Claude Desktop's config files say about nyx.
    pub(crate) fn refresh_connect(&mut self) {
        self.connect_status = match (Env::detect(), &self.mcp_server) {
            (Ok(env), Some(server)) => nyx_connect::list(&env, server).into_iter().filter(|h| h.status != Status::NotInstalled).collect(),
            _ => vec![],
        };
    }

    /// The Connect menu: one row per Claude app with a Connect/Disconnect button, so the AI interface (`nyx-mcp`)
    /// is registered without editing a config file.
    pub(crate) fn connect_menu(&mut self, ui: &mut egui::Ui) {
        ui.set_min_width(320.0);
        ui.label(RichText::new("Let Claude run nyx").strong());
        ui.weak("Registers the nyx MCP server with the app, so Claude can set up and solve simulations here.");
        ui.separator();
        let Some(server) = self.mcp_server.clone() else {
            ui.colored_label(crate::theme::activity_text(ui.ctx()), "nyx-mcp was not found next to nyx.");
            ui.weak("Build it (cargo build --release -p nyx-mcp) or point to it:");
            if ui.button("Locate nyx-mcp...").clicked() {
                ui.close();
                if let Some(p) = rfd::FileDialog::new().set_title("nyx-mcp").pick_file() {
                    self.mcp_server = Some(nyx_connect::plain_path(p));
                    self.refresh_connect();
                }
            }
            return;
        };
        let Ok(env) = Env::detect() else {
            ui.colored_label(crate::theme::danger(ui.ctx()), "cannot find the home folder");
            return;
        };
        if self.connect_status.is_empty() {
            ui.weak("No supported AI app found (Claude Code, Claude Desktop, Cursor, VS Code, Windsurf, Antigravity).");
            ui.weak("Paste the instructions below into your AI app and it connects nyx itself, or copy the config JSON and paste it by hand.");
        }
        for host in self.connect_status.clone() {
            let status = host.status.clone();
            ui.horizontal(|ui| {
                ui.label(RichText::new(host.name).strong());
                match &status {
                    Status::Connected => ui.colored_label(crate::theme::success(ui.ctx()), "connected"),
                    Status::Stale { .. } => ui.colored_label(crate::theme::activity_text(ui.ctx()), "points to another nyx-mcp"),
                    Status::Unreadable { .. } => ui.colored_label(crate::theme::danger(ui.ctx()), "config unreadable"),
                    Status::NotConnected | Status::NotInstalled => ui.weak("not connected"),
                };
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let unreadable = matches!(status, Status::Unreadable { .. });
                    if status == Status::Connected {
                        if ui.button("Disconnect").clicked() {
                            self.connect_note = Some(match nyx_connect::disconnect(&env, host.id) {
                                Ok(_) => format!("{}: disconnected.", host.name),
                                Err(e) => e,
                            });
                            self.refresh_connect();
                        }
                    } else if ui.add_enabled(!unreadable, egui::Button::new(if matches!(status, Status::Stale { .. }) { "Reconnect" } else { "Connect" })).clicked() {
                        self.connect_note = Some(match nyx_connect::connect(&env, host.id, &server) {
                            Ok(_) => format!("{}: connected. {}", host.name, host.after_connect),
                            Err(e) => e,
                        });
                        self.refresh_connect();
                    }
                });
            });
            ui.weak(host.config.display().to_string());
            if let Status::Unreadable { why } = &status {
                ui.colored_label(crate::theme::danger(ui.ctx()), why);
            }
            ui.add_space(4.0);
        }
        if let Some(n) = &self.connect_note {
            ui.separator();
            ui.label(n);
        }
        ui.separator();
        ui.weak(format!("Server: {}", server.display()));
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy instructions for an AI").on_hover_text("A text for any AI app that can edit files: it tells the app which file to change and what to add.").clicked() {
                ui.ctx().copy_text(nyx_connect::agent_instructions(&env, &server));
                self.connect_note = Some("Instructions copied; paste them into your AI app's chat and let it edit its own MCP configuration.".into());
            }
            if ui.button("Copy config JSON").clicked() {
                ui.ctx().copy_text(nyx_connect::snippet(&server));
                self.connect_note = Some("Config JSON copied; paste it into the mcpServers of any MCP client.".into());
            }
            if ui.button("Locate nyx-mcp...").clicked() {
                if let Some(p) = rfd::FileDialog::new().set_title("nyx-mcp").pick_file() {
                    self.mcp_server = Some(nyx_connect::plain_path(p));
                    self.refresh_connect();
                }
            }
        });
    }

    /// The File > New window: physics, 2D or 3D; Create starts a blank project.
    fn new_dialog_ui(&mut self, ctx: &egui::Context) {
        let Some(mut d) = self.new_dialog.take() else { return };
        let (mut open, mut create, mut cancel) = (true, false, false);
        let ids: Vec<String> = self.engine.registry().ids().iter().map(|s| s.to_string()).collect();
        egui::Window::new("New simulation").open(&mut open).collapsible(false).resizable(false).show(ctx, |ui| {
            egui::Grid::new("new_sim").num_columns(2).spacing([12.0, 8.0]).show(ui, |ui| {
                ui.label("Physics");
                egui::ComboBox::from_id_salt("new_physics").selected_text(self.physics_name(&d.physics)).show_ui(ui, |ui| {
                    for p in &ids {
                        ui.selectable_value(&mut d.physics, p.clone(), self.physics_name(p));
                    }
                });
                ui.end_row();
                ui.label("Space");
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut d.three_d, false, "2D (cross-section)");
                    ui.selectable_value(&mut d.three_d, true, "3D (the section swept along y)");
                });
                ui.end_row();
            });
            ui.add_space(6.0);
            ui.weak("The new project has typical materials and no boundary conditions; everything can be changed in the tabs on the left.");
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                if ui.button(RichText::new("  Create  ").strong()).clicked() {
                    create = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
            });
        });
        if create {
            let what = format!("{}, {}", self.physics_name(&d.physics), if d.three_d { "3D" } else { "2D" });
            self.new_project(&d.physics, d.three_d, &what);
        } else if open && !cancel {
            self.new_dialog = Some(d);
        }
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Physics");
            let mut active = self.active().unwrap_or_default();
            let before = active.clone();
            let ids: Vec<String> = self.engine.project().physics.iter().map(|s| s.physics.clone()).collect();
            egui::ComboBox::from_id_salt("active_physics").selected_text(self.physics_name(&active)).show_ui(ui, |ui| {
                for p in &ids {
                    ui.selectable_value(&mut active, p.clone(), self.physics_name(p));
                }
            });
            if active != before {
                self.cmd("set_active", json!({ "physics": active }));
                self.picking = None;
                self.rules_dirty = true;
            }
            let has_mesh = self.doc.mesh.is_some();
            // behind a running job a Compute waits in the queue (without a workspace; a workspace keeps one run at a time)
            let queue_label = if self.busy.is_some() { format!("  Queue · position {}  ", self.waiting.len() + 1) } else { "  Compute  ".to_string() };
            let b = ui.add_enabled(self.can_queue_compute() && has_mesh, egui::Button::new(RichText::new(queue_label).strong().color(crate::theme::on_accent(ui.ctx()))).fill(crate::theme::accent(ui.ctx())));
            if b.on_hover_text("Solve the active physics on the current mesh (F5)").clicked() {
                self.compute();
            }
            if let Some((what, t)) = &self.busy {
                ui.spinner();
                let st = self.job.as_ref().map(|j| j.1.status()).unwrap_or_default();
                let stage = if st.stage.is_empty() { String::new() } else { format!(": {}", st.stage) };
                let pct = st.fraction.map(|f| format!(" ({:.0} %)", 100.0 * f)).unwrap_or_default();
                ui.label(format!("{what}{stage}{pct}  [{:.1} s in this stage, {:.1} s in all]", st.stage_seconds, t.elapsed().as_secs_f64()));
                if self.job.as_ref().is_some_and(|j| j.2) {
                    let cancelled = st.cancelled;
                    if cancelled {
                        // what the person waits for is the stage, not an estimate: how long it has been running
                        ui.colored_label(crate::theme::activity_text(ui.ctx()), format!("Cancel requested: waiting for the stage to end ({:.1} s so far; it stops at its next check)", st.stage_seconds));
                    }
                    if ui.add_enabled(!cancelled, egui::Button::new("Cancel")).on_hover_text("Stop the job at the next safe point (Esc). A factorisation in progress has to finish first. Earlier results stay as they are.").clicked() {
                        self.cancel_job();
                    }
                }
                if !self.waiting.is_empty() {
                    ui.label(format!("+{} queued", self.waiting.len()));
                }
                if ui.small_button("Queue").on_hover_text("What runs and what waits, with cancel, hold and pause").clicked() {
                    self.show_queue = !self.show_queue;
                }
                ui.ctx().request_repaint();
            }
            ui.separator();
            if !self.results.is_empty() {
                ui.label("Show:");
                if ui.selectable_label(!self.show_results, "mesh").clicked() {
                    self.show_results = false;
                }
                if ui.selectable_label(self.show_results, "result").clicked() {
                    self.show_results = true;
                }
            }
        });
    }

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if let Some(e) = self.error.clone() {
                ui.colored_label(crate::theme::danger(ui.ctx()), format!("⚠ {e}"));
                if ui.small_button("dismiss").clicked() {
                    self.error = None;
                }
            } else {
                ui.label(&self.status);
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let order = self.active().and_then(|a| self.engine.project().spec(&a).map(|s| s.order as usize)).unwrap_or(2);
                if let Some(m) = &self.doc.mesh {
                    ui.weak(format!("{}  |  {}D  |  {} points  {} cells  |  P{} nodes: {}", self.doc.source, m.dim, m.n_points(), m.n_cells(), order, m.n_dofs(order)));
                } else if let Some(v) = &self.doc.view {
                    ui.weak(format!("{}  |  {}D  |  view only, {} cells", self.doc.source, v.dim, v.n_cells));
                }
            });
        });
    }

    /// The colour plan of the model on show (`regions.rs`): the regions of a 2D mesh keyed by material or by their own name (Settings: Colour regions
    /// by), the surfaces of a 3D mesh by boundary number (a surface has no region; the materials belong to the volume).
    ///
    /// The plan reads every physics of the project, in project order, never only the active one (ADR 123): the first physics that gives a region a
    /// material names its colour, and a region that any physics marks as a PML is an absorbing layer, drawn in one grey graded by its strength. So
    /// switching the physics on show never repaints the model.
    pub fn region_paint(&self, view: &nyx_engine::facade::ViewMesh) -> crate::regions::Paint {
        use crate::regions::{assign, assign_with_pml, keys_by_material};
        use nyx_engine::facade::physics::{LabelResolver, MaterialAssign};
        use std::collections::BTreeMap;
        if view.dim == 3 {
            let labels = if view.boundary_labels.is_empty() { &view.tri_label } else { &view.boundary_labels };
            return assign(labels, |l| l.to_string(), false);
        }
        let labels = if view.domain_labels.is_empty() { &view.tri_label } else { &view.domain_labels };
        let name_of = |l: i32| self.domain_names.iter().find(|(i, _)| *i == l).map_or_else(|| l.to_string(), |(_, n)| n.clone());
        let by_material = self.shell.settings().region_colouring != "region";
        let mut materials: BTreeMap<String, Option<String>> = BTreeMap::new();
        let mut layers: BTreeMap<i32, f64> = BTreeMap::new();
        if let Some(model) = self.doc.model.as_ref() {
            for spec in &self.engine.project().physics {
                for (region, m) in &spec.materials {
                    let material = match m {
                        MaterialAssign::Named(n) => Some(n.clone()),
                        MaterialAssign::Oriented { material, .. } => Some(material.clone()),
                        MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_) => None,
                    };
                    let entry = materials.entry(region.clone()).or_insert(None);
                    if entry.is_none() {
                        *entry = material;
                    }
                }
                for bc in spec.bcs.iter().filter(|bc| bc.kind == "pml") {
                    let strength = bc.params.get("strength").and_then(|v| v.as_f64()).unwrap_or(3.0);
                    for l in LabelResolver::domains(model.as_ref(), &bc.on).unwrap_or_default() {
                        let e = layers.entry(l).or_insert(strength);
                        *e = e.max(strength);
                    }
                }
            }
        }
        let key: HashMap<i32, String> = match (by_material, self.doc.model.as_ref()) {
            (true, Some(model)) => keys_by_material(&self.domain_names, &materials, |p| LabelResolver::domains(model.as_ref(), p).unwrap_or_default()),
            _ => HashMap::new(),
        };
        assign_with_pml(labels, |l| key.get(&l).cloned().unwrap_or_else(|| name_of(l)), true, &layers)
    }

    pub fn viewport_ui(&mut self, ui: &mut egui::Ui) {
        // everything on the canvas uses the colours of the canvas: dark by default, also in a light window (Settings: Canvas)
        let ink = crate::theme::ink(ui.ctx(), self.shell.settings().canvas == "theme");
        let mut bc_colors: HashMap<i32, Color32> = HashMap::new();
        for (i, labels) in self.bc_labels.iter().enumerate() {
            for &l in labels {
                bc_colors.insert(l, crate::colormap::palette(i + 3, ink.mode));
            }
        }
        let selected: Vec<i32> = if self.picking.is_some() { self.picked.clone() } else { vec![] };
        let picking = self.picking.is_some();
        let cmap = self.colormap;
        let (show_edges, show_numbers, show_results) = (self.show_edges, self.show_numbers, self.show_results);
        // while the geometry is being edited the exact model is drawn, not the mesh of an earlier version
        let preview = if self.tab == Tab::Geometry && (self.geometry_dirty || self.geom_edit != self.engine.project().geometry) { self.geom_ui.preview.clone() } else { None };
        let show_results = show_results && preview.is_none();
        let view = preview.or_else(|| self.doc.view.clone());

        let stale = if show_results { self.current.and_then(|i| self.stale_reason(i)) } else { None };
        let out = match (self.current, view) {
            (Some(i), Some(view)) if show_results => {
                let rv = &mut self.results[i];
                let values = rv.field_values();
                let range = rv.range();
                let warp = rv.displacement().filter(|_| rv.warp != 0.0).map(|d| (d, rv.warp));
                let sol = rv.sol.min(rv.sols.len() - 1);
                let nf = rv.sols[sol].fields.len();
                let f = &rv.sols[sol].fields[rv.field.min(nf - 1)];
                let (title, unit) = (f.name.clone(), f.unit.clone());
                let rv = &self.results[i];
                self.vp.show(
                    ui,
                    Content::Field {
                        display: &rv.display,
                        view: &view,
                        values: &values,
                        range,
                        cmap,
                        unit: &unit,
                        title: &title,
                        warp: warp.as_ref().map(|(d, s)| (d.as_slice(), *s)),
                        show_edges,
                    },
                    false,
                    &ink,
                )
            }
            (_, Some(view)) => {
                let paint = self.region_paint(&view);
                self.vp.show(ui, Content::Mesh { view: &view, bc_colors: &bc_colors, selected: &selected, show_edges, show_numbers, paint: &paint }, picking, &ink)
            }
            _ => self.vp.show(ui, Content::Empty, false, &ink),
        };
        if let Some(why) = stale {
            // on the plot itself, so a stale field is never read as current; the button to solve again is in the Results tab and its panel.
            // A small chip in the style of the trust chip (an amber ring and the word), not a filled amber band over the data; the reason is its tooltip (ADR 123)
            let ctx = ui.ctx().clone();
            let at = out.rect.left_top() + egui::vec2(10.0, 10.0);
            let galley = ui.painter().layout_no_wrap("Stale".to_string(), egui::FontId::proportional(14.0), crate::theme::text(&ctx));
            let rect = egui::Rect::from_min_size(at, egui::vec2(galley.size().x + 36.0, galley.size().y + 8.0));
            ui.painter().rect(rect, 4.0, crate::theme::card(&ctx), egui::Stroke::new(1.0, crate::theme::edge(&ctx)), egui::StrokeKind::Inside);
            crate::kit::paint_glyph(ui, rect.left_center() + egui::vec2(14.0, 0.0), crate::kit::State::Waiting);
            ui.painter().galley(rect.left_top() + egui::vec2(28.0, 4.0), galley, crate::theme::text(&ctx));
            let _ = ui.interact(rect, egui::Id::unique("nyx-stale-chip"), egui::Sense::hover()).on_hover_text(format!("Stale: {why}. Compute again to refresh this field."));
        } else if show_results {
            // the trust state of the field on show (P6): a chip on the plot, in the colours of the window, so it reads on either canvas and does not
            // depend on the side panel being wide enough to show the same line
            if let Some(physics) = self.current.and_then(|i| self.results.get(i)).and_then(|r| r.physics.clone()) {
                let v = self.validated_of(&physics);
                let state = crate::activity::validation_state(v.state);
                let ctx = ui.ctx().clone();
                let at = out.rect.left_top() + egui::vec2(10.0, 10.0);
                let galley = ui.painter().layout_no_wrap(state.word().to_string(), egui::FontId::proportional(14.0), crate::theme::text(&ctx));
                let rect = egui::Rect::from_min_size(at, egui::vec2(galley.size().x + 36.0, galley.size().y + 8.0));
                ui.painter().rect(rect, 4.0, crate::theme::card(&ctx), egui::Stroke::new(1.0, crate::theme::edge(&ctx)), egui::StrokeKind::Inside);
                crate::kit::paint_glyph(ui, rect.left_center() + egui::vec2(14.0, 0.0), state);
                ui.painter().galley(rect.left_top() + egui::vec2(28.0, 4.0), galley, crate::theme::text(&ctx));
            }
        }
        if let (Some(label), Some(i)) = (out.picked, self.picking) {
            if let Some(p) = self.picked.iter().position(|&l| l == label) {
                self.picked.remove(p);
            } else {
                self.picked.push(label);
                self.picked.sort();
            }
            self.suggestions = self.doc.model.as_ref().map(|m| m.suggest(&self.picked)).unwrap_or_default();
            // the simplest selector that picks exactly this set becomes the condition's; otherwise the person chooses
            if let Some(first) = self.suggestions.first().cloned() {
                self.set_condition_on(i, &first);
            }
        }
    }

    /// Sets the selector of condition `i` of the active physics (an edit of the problem through `set_problem`).
    pub fn set_condition_on(&mut self, i: usize, on: &str) {
        let Some(active) = self.active() else { return };
        let Some(mut spec) = self.engine.project().spec(&active).cloned() else { return };
        if let Some(bc) = spec.bcs.get_mut(i) {
            bc.on = on.to_string();
            self.set_problem(&active, &spec);
        }
    }

    /// Replaces the problem of a physics by `spec` (the Physics and Study tabs edit a copy and hand it back here).
    pub fn set_problem(&mut self, physics: &str, spec: &nyx_engine::facade::physics::ProblemSpec) {
        let v = serde_json::to_value(spec).unwrap_or(Value::Null);
        if self.cmd("set_problem", json!({ "physics": physics, "problem": v })).is_some() {
            self.rules_dirty = true;
        }
    }
}

impl eframe::App for App {
    fn logic(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        crate::spans::span("poll", || self.poll());
        crate::spans::span("run_queued", || self.run_queued());
        if self.rules_dirty {
            crate::spans::span("sync_rules", || self.sync_rules());
        }
        let c = ctx.clone();
        crate::spans::span("script", || self.script_step(&c));
        let dropped: Vec<PathBuf> = ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        if let Some(p) = dropped.into_iter().next() {
            self.open_path(p);
        }
        let (f5, esc, ctrl_shift_r, alt1, alt2, alt3, alt4, alt5, alt6) = ctx.input_mut(|i| (
            i.consume_key(egui::Modifiers::NONE, egui::Key::F5),
            i.key_pressed(egui::Key::Escape),
            i.consume_key(egui::Modifiers::COMMAND | egui::Modifiers::SHIFT, egui::Key::R),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num1),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num2),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num3),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num4),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num5),
            i.consume_key(egui::Modifiers::ALT, egui::Key::Num6),
        ));
        if f5 {
            self.shell_compute();
        }
        if ctrl_shift_r {
            self.toggle_nys_script();
        }
        for (pressed, act) in [
            (alt1, crate::menu::Act::ShowGeometry),
            (alt2, crate::menu::Act::ShowMaterials),
            (alt3, crate::menu::Act::ShowPhysics),
            (alt4, crate::menu::Act::ShowStudy),
            (alt5, crate::menu::Act::ShowSweep),
            (alt6, crate::menu::Act::ShowResults),
        ] {
            if pressed {
                self.do_act(ctx, act);
            }
        }
        if esc && self.busy.is_some() {
            if self.shell.active() {
                // Esc never cancels silently in the shell: it asks first
                self.shell.confirm_cancel = true;
            } else {
                self.cancel_job();
            }
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let t_ui = Instant::now();
        crate::spans::span("queue_panel", || self.queue_panel(ui.ctx()));
        if let Some(text) = self.script_output.clone() {
            let mut open = true;
            egui::Window::new("Python output").open(&mut open).default_size([560.0, 320.0]).show(ui.ctx(), |ui| {
                egui::ScrollArea::both().show(ui, |ui| {
                    ui.add(egui::TextEdit::multiline(&mut text.as_str()).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
                });
            });
            if !open {
                self.script_output = None;
            }
        }
        if self.shell.active() {
            egui::Panel::top("shell_menubar").show(ui, |ui| self.bar_ui(ui));
            crate::spans::span("shell_ui", || self.shell_ui(ui));
            self.frame_csv(t_ui.elapsed().as_secs_f64() * 1e3, f64::NAN);
            return;
        }
        if self.nys.open || self.shell.sheet == crate::shell::Sheet::Open {
            self.sheets_ui(ui.ctx());
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::N)) {
            self.new_dialog = Some(NewDialog { physics: self.active().unwrap_or_else(|| "dfs".into()), three_d: false });
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::O)) {
            self.open_dialog();
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::S)) {
            self.save_project(false);
        }
        egui::Panel::top("menu").show(ui, |ui| self.menu_bar(ui));
        egui::Panel::top("toolbar").show(ui, |ui| {
            ui.add_space(2.0);
            self.toolbar(ui);
            ui.add_space(2.0);
        });
        egui::Panel::bottom("status").show(ui, |ui| self.status_bar(ui));
        egui::Panel::left("inspector").default_size(400.0).size_range(320.0..=700.0).resizable(true).show(ui, |ui| self.inspector(ui));
        if self.bottom_open && self.show_results && self.current.is_some_and(|i| self.results[i].sols.len() > 1) {
            egui::Panel::bottom("results_panel").default_size(190.0).resizable(true).show(ui, |ui| self.results_panel(ui));
        }
        self.new_dialog_ui(ui.ctx());
        let t_vp = Instant::now();
        egui::CentralPanel::default().frame(egui::Frame::NONE).show(ui, |ui| self.viewport_ui(ui));
        if std::env::var_os("NYX_FRAME_LOG").is_some() {
            // diagnostics: time (s, wrapping every 1000) of the frame, our own layout and drawing code, whole viewport, state
            let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs_f64() % 1000.0).unwrap_or(0.0);
            eprintln!("FRAME t={now:.3} ui_ms={:.2} viewport_ms={:.2} busy={} show_results={}", t_ui.elapsed().as_secs_f64() * 1e3, t_vp.elapsed().as_secs_f64() * 1e3, self.busy.is_some(), self.show_results);
            self.frame_csv(t_ui.elapsed().as_secs_f64() * 1e3, t_vp.elapsed().as_secs_f64() * 1e3);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn app() -> App {
        App::with_context(egui::Context::default())
    }

    /// One headless frame. The texture changes of the frame (the font atlas) are cleared on purpose: egui asserts when a frame's output is
    /// dropped with unapplied texture deltas, and no window is there to apply them.
    fn frame(ctx: &egui::Context, run: impl FnMut(&mut egui::Ui)) {
        let mut out = ctx.run_ui(Default::default(), run);
        out.textures_delta.clear();
    }

    /// Processes finished jobs until nothing runs (or panics after `secs`).
    fn wait_idle(a: &mut App, secs: u64) {
        let t0 = Instant::now();
        loop {
            a.poll();
            a.run_queued();
            if a.busy.is_none() && a.queued.is_none() && a.exports.is_empty() && a.deferred.is_empty() {
                return;
            }
            assert!(t0.elapsed() < Duration::from_secs(secs), "still busy after {secs} s: {:?}", a.busy);
            std::thread::sleep(Duration::from_millis(10));
        }
    }

    fn example(name_part: &str) -> usize {
        examples::list().iter().position(|e| e.path.to_string_lossy().contains(name_part)).unwrap_or_else(|| panic!("no example {name_part}"))
    }

    #[test]
    fn a_new_simulation_starts_empty_and_can_be_set_up_and_solved() {
        let mut a = app();
        wait_idle(&mut a, 30);
        for physics in ["dfs", "sms", "ems"] {
            for three_d in [false, true] {
                let mut e = Engine::new();
                examples::blank(&mut e, physics, three_d).unwrap_or_else(|err| panic!("{physics} 3d={three_d}: {err}"));
                let p = e.project();
                assert_eq!(p.active.as_deref(), Some(physics));
                assert!(p.physics.iter().all(|s| s.bcs.is_empty()), "no boundary conditions in a new project");
                e.call("mesh_build", &json!({})).unwrap_or_else(|err| panic!("{physics} 3d={three_d}: {err}"));
                assert_eq!(e.model().unwrap().mesh().dim(), if three_d { 3 } else { 2 });
            }
        }
        a.new_project("dfs", false, "test");
        wait_idle(&mut a, 30);
        assert!(a.results.is_empty() && a.doc.mesh.is_some() && a.project_path.is_none());
        // set up from scratch: plates on the left and right, then compute
        a.cmd("add_condition", json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "potential", "params": { "value": 1.0 } }));
        a.cmd("add_condition", json!({ "physics": "dfs", "on": "outer & plane(x = 1)", "kind": "ground" }));
        a.compute();
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert_eq!(a.results.len(), 1);
    }

    #[test]
    fn an_example_is_built_computed_and_shown() {
        let mut a = app();
        a.load_example(example("01_capacitor"));
        assert!(a.busy.is_some());
        wait_idle(&mut a, 30);
        assert!(a.doc.mesh.is_some(), "{:?}", a.error);
        a.compute();
        assert!(a.busy.is_some(), "compute started: {} / {:?}", a.status, a.error);
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert_eq!(a.results.len(), 1);
        assert!(a.status.contains("solved"), "{}", a.status);
        // the stage log of the run is part of the result
        assert!(a.results[0].result.log.iter().any(|l| l.contains("solving")), "{:?}", a.results[0].result.log);
    }


    #[test]
    fn the_reduce_form_shows_the_numbers_of_the_engine_commands() {
        let mut a = app();
        a.load_example(example("15_soi_wire"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert_eq!(a.results.len(), 1);
        let physics = a.results[0].physics.clone().unwrap();
        let sol = a.results[0].result.solutions[0].clone();
        let names: Vec<String> = sol.cell_fields.iter().map(|f| f.name.clone()).collect();
        let ex = names.iter().position(|n| n == "Ex").expect("a field Ex");
        let regions = a.engine.result_regions(&physics);
        let si = regions.iter().position(|r| r == "Si").expect("a region Si") + 1;
        // the integral of |Ex|^2 over the silicon: the line of the form carries the engine's number
        a.reduce.field = ex;
        a.reduce.quantity = 2;
        a.reduce.region = si;
        a.reduce_now("integrate");
        let want = a.engine.call("integrate", &json!({ "expr": { "op": "abs2", "args": ["Ex"] }, "physics": physics, "region": "Si" })).unwrap();
        assert_eq!(a.reduce.lines[0], format!("integral of |Ex|² over Si = {}", crate::tabs::reduce_number(&want["value"])));
        // the share of the whole vector in the silicon is the engine's `confinement` over the components `field` finds (Ex, Ey, Ez: a mode has the propagation component too)
        a.reduce_now("fraction");
        let conf = a.engine.call("confinement", &json!({ "fields": ["Ex", "Ey", "Ez"], "physics": physics, "region": "Si" })).unwrap();
        assert_eq!(a.reduce.lines[0], format!("share of |Ex, Ey, Ez|² in Si = {:.6}", conf["fraction"].as_f64().unwrap()));
        // the largest squared modulus over the whole mesh, then the largest plain value (the mode of this example is real, so it has an order)
        a.reduce.region = 0;
        a.reduce_now("max");
        assert!(a.reduce.lines[0].starts_with("max of |Ex|² in the whole mesh = "), "{:?}", a.reduce.lines);
        a.reduce.quantity = 0;
        a.reduce_now("min");
        assert!(a.reduce.lines[0].starts_with("min of Ex in the whole mesh = "), "{:?}", a.reduce.lines);
        assert!(a.error.is_none(), "{:?}", a.error);
        // a formula in the box replaces the field and quantity: the line carries the number of the engine for the same text
        a.reduce.region = si;
        a.reduce.formula = "  eps_r * norm2(E)  ".into();
        a.reduce_now("integrate");
        let want = a.engine.call("integrate", &json!({ "expr": "eps_r * norm2(E)", "physics": physics, "region": "Si" })).unwrap();
        assert_eq!(a.reduce.lines[0], format!("integral of `eps_r * norm2(E)` over Si = {}", crate::tabs::reduce_number(&want["value"])));
        // a formula that cannot be read is an error that names it, and no line is added
        let before = a.reduce.lines.len();
        a.reduce.formula = "abs2(Ex) +".into();
        a.reduce_now("integrate");
        assert!(a.error.is_some(), "a broken formula is reported");
        assert_eq!(a.reduce.lines.len(), before);
    }
    #[test]
    fn a_compute_behind_a_running_job_waits_its_turn_and_a_newer_one_replaces_it() {
        let mut a = app();
        a.load_example(example("05_block"));
        assert!(a.busy.is_some());
        let first = a.job.as_ref().map(|j| j.0);
        a.engine.jobs().pause(); // the mesh job may be over already: nothing new starts, so the queue can be looked at
        a.compute(); // a job is running: the solve waits in the engine's queue behind it
        assert_eq!(a.waiting.len(), 1, "{} / {:?}", a.status, a.error);
        assert!(a.status.contains("queued"), "{}", a.status);
        let queued = a.waiting[0].id;
        assert_eq!(a.engine.jobs().info(queued).unwrap().position, Some(1));
        assert_eq!(a.job.as_ref().map(|j| j.0), first, "the running job is not replaced");
        a.compute(); // a Compute that has not started is replaced by the newer one
        assert_eq!(a.waiting.len(), 1);
        assert_ne!(a.waiting[0].id, queued);
        a.rebuild_mesh(); // a mesh build still has its one slot
        assert!(a.queued.is_some());
        // the Queue panel lists what runs and what waits, and draws on a headless context
        a.show_queue = true;
        let ctx = a.ctx.clone();
        frame(&ctx, |ui| a.queue_panel(ui.ctx()));
        a.engine.jobs().resume();
        wait_idle(&mut a, 120);
        assert!(a.doc.mesh.is_some() && a.queued.is_none() && a.waiting.is_empty());
        assert!(!a.results.is_empty(), "the queued Compute ran when the first job ended and its result is on show");
    }

    #[test]
    fn a_refused_mesh_edit_cannot_leave_the_geometry_edit_applied() {
        let mut a = app();
        a.load_example(example("15_soi_wire"));
        wait_idle(&mut a, 30);
        let before = a.engine.project().clone();
        nyx_engine::geometry::set_param(&mut a.geom_edit, "w", &json!(6.0e-7)).unwrap();
        a.mesh_edit.h = -1.0;
        assert!(!a.push_edits(), "a negative element size is refused");
        assert_eq!(*a.engine.project(), before, "neither edit was applied");
        // every mesh field the window sends is accepted by the engine
        a.mesh_edit.h = 2.0e-7;
        a.mesh_edit.max_dev = Some(1.0e-9);
        a.mesh_edit.points_per_wavelength = Some(10.0);
        a.mesh_edit.structured_pml = true;
        assert!(a.push_edits(), "{:?}", a.error);
        assert_eq!(a.engine.project().mesh, a.mesh_edit);
        assert_eq!(a.engine.project().geometry, a.geom_edit);
    }

    #[test]
    fn a_compute_that_waits_in_the_queue_can_be_cancelled_and_never_runs() {
        let mut a = app();
        a.load_example(example("05_block"));
        let running = a.job.as_ref().map(|j| j.0).unwrap();
        a.engine.jobs().pause(); // nothing new starts while the queue is looked at
        a.compute();
        let queued = a.waiting[0].id;
        a.cancel_waiting(queued);
        assert!(a.waiting.is_empty());
        assert_eq!(a.engine.jobs().info(queued).unwrap().state, nyx_engine::JobState::Cancelled);
        assert_eq!(a.job.as_ref().map(|j| j.0), Some(running), "the running job goes on");
        a.engine.jobs().resume();
        wait_idle(&mut a, 60);
        assert!(a.results.is_empty(), "the cancelled Compute never ran");
        assert!(a.engine.jobs().log(queued).is_empty());
    }

    #[test]
    fn a_message_of_a_script_does_not_clear_the_busy_state() {
        let mut a = app();
        a.load_example(example("01_capacitor"));
        let running = a.job.as_ref().map(|j| j.0).unwrap();
        a.tx.send((running + 100, Msg::Script(Ok("from a script".into())))).unwrap();
        a.poll();
        assert_eq!(a.script_output.as_deref(), Some("from a script"));
        assert!(a.busy.is_some() || a.job.is_none(), "only the running job's own end clears it");
        wait_idle(&mut a, 30);
        assert!(a.busy.is_none() && a.job.is_none());
    }

    #[test]
    fn a_long_solve_can_be_cancelled_and_earlier_results_stay() {
        let mut a = app();
        a.load_example(example("01_capacitor"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 60);
        assert_eq!(a.results.len(), 1);
        // a big 3D eigenfrequency problem on the same window
        a.load_example(example("05_block"));
        wait_idle(&mut a, 30);
        a.cmd("set_mesh", json!({ "h": 0.1, "sweep": [{ "thickness": 6.0, "n": 60 }] }));
        a.cmd("set_physics_params", json!({ "physics": "sms", "params": { "solver": "iterative" } }));
        a.cmd("set_study", json!({ "physics": "sms", "kind": "eigenfrequency", "params": { "n_modes": 6, "f_search": 0.0 }, "order": 2 }));
        // a boundary load belongs to the static study: the eigenfrequency study has only the clamp
        a.cmd("remove_condition", json!({ "physics": "sms", "index": 1 }));
        a.adopt();
        wait_idle(&mut a, 60);
        let kept = a.results.len();
        a.compute();
        assert!(a.busy.is_some(), "{} / {:?}", a.status, a.error);
        std::thread::sleep(Duration::from_millis(1500));
        let t = Instant::now();
        a.cancel_job();
        wait_idle(&mut a, 90);
        assert!(t.elapsed() < Duration::from_secs(60), "cancelled after {:?}", t.elapsed());
        assert!(a.status.contains("cancelled"), "{} / {:?}", a.status, a.error);
        assert!(a.error.is_none());
        assert_eq!(a.results.len(), kept, "a cancelled solve adds no result");
    }

    #[test]
    fn a_project_saves_and_loads_again() {
        let mut a = app();
        a.load_example(example("06_wr90"));
        wait_idle(&mut a, 30);
        let path = std::env::temp_dir().join(format!("nyx-test-{}.nyx.json", std::process::id()));
        a.project_path = Some(path.clone());
        a.save_project(false);
        assert!(path.exists(), "{}", a.status);
        let mut b = app();
        b.load_project(&path);
        wait_idle(&mut b, 30);
        assert_eq!(a.engine.project().to_json(), b.engine.project().to_json());
        assert!(b.doc.mesh.is_some());
        let _ = std::fs::remove_file(&path);
    }

    /// A results file written by one window opens in another: the project, the mesh and the stored results are on show, with no solve.
    /// A result on show says when it is stale (productivity plan P1): an unapplied edit of the Geometry tab, an edit the engine has taken, a
    /// Study edit; looking at it differently (colormap, camera, field) is not a change of the model.
    #[test]
    fn a_result_on_show_says_when_the_model_changed() {
        let mut a = app();
        a.load_example(example("15_soi_wire"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 120);
        assert!(a.error.is_none(), "{:?}", a.error);
        let i = a.current.expect("a result on show");
        assert_eq!(a.stale_reason(i), None);
        // the way it is drawn is not in the project, so no key moves
        let keys = a.engine.problem_keys();
        a.colormap = Colormap::Gray;
        a.vp.fit_token += 1;
        a.results[i].field = a.results[i].field.saturating_sub(1);
        assert_eq!(a.engine.problem_keys(), keys);
        assert_eq!(a.stale_reason(i), None);
        // an edit of the Geometry tab that was not applied: the engine does not know it yet, the window does
        nyx_engine::geometry::set_param(&mut a.geom_edit, "w", &json!(6.0e-7)).unwrap();
        assert_eq!(a.stale_reason(i), Some("Geometry or Mesh edits are not applied yet"));
        // applied: now the engine says so
        assert!(a.push_edits(), "{:?}", a.error);
        assert_eq!(a.stale_reason(i), Some("the model changed since this solve"));
        // the band is drawn on a headless context without a panic, with the button to solve again
        let ctx = a.ctx.clone();
        frame(&ctx, |ui| {
            a.results_panel(ui);
            a.viewport_ui(ui);
        });
        // a Study edit (a wavelength) is a change of the problem of its physics
        let mut b = app();
        b.load_example(example("15_soi_wire"));
        wait_idle(&mut b, 30);
        b.compute();
        wait_idle(&mut b, 120);
        let j = b.current.expect("a result on show");
        assert_eq!(b.stale_reason(j), None);
        b.cmd("set_study", json!({ "physics": "ems", "kind": "modes", "params": { "freq": 190.0e12, "n_modes": 2, "n_target": 2.75 } }));
        assert_eq!(b.stale_reason(j), Some("the model changed since this solve"));
    }

    /// A sweep started from a workflow (UX plan U5): the points are jobs of the engine's queue that the chip and the popover show, each finished
    /// point is recorded as a run of the experiment with the metrics of a Compute, and the window's own project is not touched.
    #[test]
    fn a_sweep_from_a_workflow_queues_its_points_and_records_each_as_a_run() {
        let mut a = app();
        let root = std::env::temp_dir().join(format!("nyx-sweep-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        a.attach_workspace(root.clone());
        assert!(a.shell.active(), "a workspace is attached");
        let path = examples::list().into_iter().find(|e| e.path.to_string_lossy().contains("15_soi_wire")).unwrap().path;
        a.new_draft("Sweep", crate::shell::Start::Example(path));
        wait_idle(&mut a, 60);
        let wid = a.shell.current_mut().expect("the workflow is open").id.clone();
        let before = a.engine.project().clone().to_json();
        let mut form = crate::sweepui::Form::default();
        (form.param, form.from, form.to, form.n) = ("w".into(), "4.5e-7".into(), "5.5e-7".into(), 3);
        a.start_sweep(&wid, &form);
        assert!(a.error.is_none(), "{:?}", a.error);
        // the engine's queue holds the points: the activity the chip draws is not idle, and says what runs
        let act = a.activity();
        assert!(!act.is_idle() && act.rows.len() <= 3, "{act:?}");
        assert!(act.chip_text().contains("point"), "{}", act.chip_text());
        let t0 = Instant::now();
        while a.sweep.as_ref().is_some_and(|s| !(s.complete() && s.points.iter().all(|p| p.recorded))) {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(120), "the sweep did not finish: {:?}", a.sweep.as_ref().map(|s| s.finished()));
            std::thread::sleep(Duration::from_millis(10));
        }
        let s = a.sweep.as_ref().expect("the sweep");
        assert!(s.points.iter().all(|p| p.state == "done"), "{:?}", s.points);
        let runs = &a.shell.open[&wid].runs;
        assert_eq!(runs.len(), 3, "one run per point");
        let first = runs[0].metrics["eigenvalue_re"];
        let last = runs[2].metrics["eigenvalue_re"];
        assert!(runs.iter().all(|r| r.status == nyx_engine::facade::store::runs::RunStatus::Ok && r.metrics.contains_key("dofs")));
        assert!((first - last).abs() > 1e-6, "the widths differ, so the numbers do: {first} {last}");
        assert!(a.activity().is_idle());
        assert_eq!(a.engine.project().clone().to_json(), before, "a sweep never changes the window's project");
        let _ = std::fs::remove_dir_all(&root);
    }

    /// ADR 124: in a workflow a new result opens its tab and shows the field.
    #[test]
    fn a_new_result_opens_the_results_tab_with_the_field() {
        let mut a = app();
        let root = std::env::temp_dir().join(format!("nyx-results-tab-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        a.attach_workspace(root.clone());
        let path = examples::list().into_iter().find(|e| e.path.to_string_lossy().contains("02_cantilever")).unwrap().path;
        a.new_draft("Fresh", crate::shell::Start::Example(path));
        wait_idle(&mut a, 30);
        assert!(a.shell.active());
        assert_eq!(a.shell.current_mut().unwrap().tab, Tab::Geometry, "a workflow opens on its geometry");
        a.shell_compute();
        wait_idle(&mut a, 120);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert!(!a.results.is_empty());
        assert!(a.show_results, "the field is on show");
        let w = a.shell.current_mut().unwrap();
        assert_eq!(w.tab, Tab::Results, "the tab of the result is open beside the field");
        assert_eq!(w.sub, crate::shell::Sub::Main);
        let _ = std::fs::remove_dir_all(&root);
    }

    /// Review I-9: a result that was just computed is never stale, whichever way it was asked for (the plain window or a workflow, which
    /// records the run and saves the model around the solve), for the physics of the examples that solve fast.
    #[test]
    fn a_result_just_computed_is_not_stale() {
        for name in ["03_heated_plate", "04_thermal_bar", "02_cantilever", "01_capacitor", "05_block_cantilever"] {
            // the plain window
            let mut a = app();
            a.load_example(example(name));
            wait_idle(&mut a, 30);
            a.compute();
            wait_idle(&mut a, 120);
            assert!(a.error.is_none(), "{name}: {:?}", a.error);
            for i in 0..a.results.len() {
                assert_eq!(a.stale_reason(i), None, "{name}: result {i} of the plain window");
            }
            // a workflow, as the scripted shots of the review ran it (`wf example N; compute; view result`)
            let mut b = app();
            let path = examples::list().into_iter().find(|e| e.path.to_string_lossy().contains(name)).unwrap().path;
            b.new_draft("Fresh", crate::shell::Start::Example(path));
            wait_idle(&mut b, 30);
            b.shell_compute();
            wait_idle(&mut b, 120);
            assert!(b.error.is_none(), "{name}: {:?}", b.error);
            assert!(!b.results.is_empty(), "{name}: the workflow solved");
            b.show_results = true;
            let keys = b.engine.problem_keys();
            // whatever the frames do after the solve (saving the workflow, the run record, the page itself) must not move a key
            let ctx = b.ctx.clone();
            for _ in 0..12 {
                frame(&ctx, |ui| b.shell_ui(ui));
                b.poll();
                std::thread::sleep(Duration::from_millis(5));
            }
            assert_eq!(b.engine.problem_keys(), keys, "{name}: a frame after the solve moved a problem key");
            for i in 0..b.results.len() {
                assert_eq!(b.stale_reason(i), None, "{name}: result {i} of a workflow, after twelve frames");
            }
        }
    }

    /// P6: a result starts as exploratory, a passed convergence check makes it checked, the window draws the state, and the verdict is
    /// written into the run record that solved it (the cantilever: its scalars move with the mesh, so a tight tolerance fails and a loose one
    /// passes).
    #[test]
    fn the_window_shows_the_validated_state_and_a_check_goes_into_the_run() {
        use nyx_engine::validated::Validation;
        let mut b = app();
        let root = std::env::temp_dir().join(format!("nyx-validated-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        b.attach_workspace(root.clone());
        let path =examples::list().into_iter().find(|e| e.path.to_string_lossy().contains("02_cantilever")).unwrap().path;
        b.new_draft("Checked", crate::shell::Start::Example(path));
        wait_idle(&mut b, 30);
        b.shell_compute();
        wait_idle(&mut b, 120);
        assert!(b.error.is_none(), "{:?}", b.error);
        assert_eq!(b.engine.validated("sms").state, Validation::Exploratory);
        let ctx = b.ctx.clone();
        frame(&ctx, |ui| b.validated_line(ui));

        // a tolerance the cantilever misses: still exploratory, and the reason says the check did not pass
        b.conv_tol = 1e-9;
        b.check_convergence();
        wait_idle(&mut b, 120);
        assert!(b.error.is_none(), "{:?}", b.error);
        let v = b.engine.validated("sms");
        assert_eq!(v.state, Validation::Exploratory, "{v:?}");
        assert!(v.reason.contains("did not pass"), "{v:?}");

        b.conv_tol = 0.5;
        b.check_convergence();
        wait_idle(&mut b, 120);
        assert_eq!(b.engine.validated("sms").state, Validation::Checked);
        frame(&ctx, |ui| b.validated_line(ui));

        let live = b.shell.live_run.clone().expect("the run whose result is on show");
        let runs = b.shell.open[&live.wf_id].dir.as_ref().expect("a workflow folder").runs(&live.exp).0;
        let run = runs.iter().find(|r| r.id == live.run).expect("the run");
        let c = run.convergence.as_ref().expect("the verdict is in the run record");
        assert!(c.passed && c.tolerance == 0.5 && c.factor == 0.6, "{c:?}");

        // a mesh edit makes the result stale: nothing is claimed for it
        b.cmd("set_mesh", json!({ "h": 0.015 }));
        assert_eq!(b.engine.validated("sms").state, Validation::Exploratory);
        let _ = std::fs::remove_dir_all(&root);
    }

    /// A solve that comes back after the project changed again and was not kept by the engine is not shown as the run on screen.
    #[test]
    fn a_result_the_engine_did_not_keep_is_not_put_on_show() {
        let mut a = app();
        a.load_example(example("15_soi_wire"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 120);
        let before = a.results.len();
        a.take_solved(&json!({ "results": [{ "physics": "ems", "current": false, "stored": false }] }));
        assert_eq!(a.results.len(), before);
        assert!(a.status.contains("not kept"), "{}", a.status);
    }

    #[test]
    fn a_results_file_opens_with_its_results_on_show() {
        let mut a = app();
        a.load_example(example("04_thermal_bar"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 120);
        assert!(a.error.is_none(), "{:?}", a.error);
        let path = std::env::temp_dir().join(format!("nyx-test-{}.h5", std::process::id()));
        let v = a.cmd("export_results", json!({ "path": path.to_str().unwrap() }));
        assert!(v.is_some(), "{:?}", a.error);
        let mut b = app();
        b.open_path(path.clone());
        // opening returns at once: reading, meshing again and checking run on a job, nothing is on show yet
        assert!(b.busy.is_some(), "{} / {:?}", b.status, b.error);
        assert!(b.doc.mesh.is_none() && b.results.is_empty(), "the file is not opened yet when open_path returns");
        wait_idle(&mut b, 30);
        assert!(b.error.is_none(), "{:?}", b.error);
        assert!(b.doc.mesh.is_some() && !b.results.is_empty(), "the results of the file are on show: mesh {} results {} status {}", b.doc.mesh.is_some(), b.results.len(), b.status);
        assert_eq!(a.engine.project().to_json(), b.engine.project().to_json());
        assert_eq!(b.tab, Tab::Results);
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(path.with_extension("xdmf"));
    }

    /// Draws every tab of the inspector, the toolbar and the viewport for several examples on a headless egui context: a panic in a form
    /// (a descriptor the window cannot draw, a project it cannot show) fails here.
    #[test]
    fn every_tab_draws_for_the_examples() {
        let mut a = app();
        for name in ["01_capacitor", "04_thermal_bar", "09_pml_slab", "15_soi_wire", "13_silicon_on_sapphire"] {
            a.load_example(example(name));
            wait_idle(&mut a, 60);
            if !name.starts_with("13") {
                a.compute();
                wait_idle(&mut a, 120);
                assert!(a.error.is_none(), "{name}: {:?}", a.error);
            }
            for t in Tab::ALL {
                a.tab = t;
                if t == Tab::Results && a.results.is_empty() || t == Tab::Analysis && !a.has_analysis() {
                    continue;
                }
                let ctx = a.ctx.clone();
                frame(&ctx, |ui| {
                    a.inspector(ui);
                    a.toolbar(ui);
                    a.status_bar(ui);
                    a.viewport_ui(ui);
                    if a.show_results {
                        a.results_panel(ui);
                    }
                });
            }
            // the mesh view with a condition being picked
            a.tab = Tab::Physics;
            a.show_results = false;
            a.picking = if a.bc_labels.is_empty() { None } else { Some(0) };
            let ctx = a.ctx.clone();
            frame(&ctx, |ui| {
                a.inspector(ui);
                a.viewport_ui(ui);
            });
        }
    }

    #[test]
    fn the_window_checks_convergence_on_its_own_thread() {
        let mut a = app();
        a.load_example(example("06_wr90"));
        wait_idle(&mut a, 30);
        a.conv_tol = 5e-3;
        a.check_convergence();
        assert!(a.convergence_running());
        assert!(a.busy.is_some(), "the check is the running job");
        let t0 = Instant::now();
        while a.convergence_running() {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(60), "no answer");
            std::thread::sleep(Duration::from_millis(10));
        }
        assert_eq!(a.convergence.as_ref().map(|c| c["converged"].clone()), Some(json!(true)), "{:?} {:?}", a.convergence, a.error);
        a.tab = Tab::Study;
        let ctx = a.ctx.clone();
        frame(&ctx, |ui| a.inspector(ui));
    }

    /// A convergence check is a job of the window: it is the busy job, it can be cancelled and it leaves the results as they were.
    #[test]
    fn the_window_cancels_a_convergence_check() {
        let mut a = app();
        a.load_example(example("14_box_cavity"));
        wait_idle(&mut a, 30);
        a.cmd("set_mesh", json!({ "h": 0.003, "sweep": [{ "domain_offset": 0, "n": 8, "thickness": 0.03 }] }));
        a.mesh_edit = a.engine.project().mesh.clone();
        a.check_convergence();
        assert!(a.convergence_running());
        let t0 = Instant::now();
        while a.job.as_ref().is_some_and(|j| j.1.status().stage.is_empty()) {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(60), "the check never reported a stage");
            std::thread::sleep(Duration::from_millis(5));
        }
        a.cancel_job();
        assert!(a.job.as_ref().is_some_and(|j| j.1.status().cancelled), "the cancel is requested while the stage runs");
        wait_idle(&mut a, 60);
        assert!(a.convergence.is_none(), "a cancelled check gives no answer");
        assert!(a.status.contains("cancelled"), "{}", a.status);
        assert!(a.error.is_none(), "{:?}", a.error);
    }

    /// Draws `n` frames of the shell (the workflow page on show), polling between them, for as long as `secs` at least.
    fn frames(a: &mut App, n: usize, secs: f64) {
        let ctx = a.ctx.clone();
        let t0 = Instant::now();
        for _ in 0..n {
            frame(&ctx, |ui| a.shell_ui(ui));
            a.poll();
            std::thread::sleep(Duration::from_millis(5));
        }
        while t0.elapsed().as_secs_f64() < secs {
            frame(&ctx, |ui| a.shell_ui(ui));
            a.poll();
            std::thread::sleep(Duration::from_millis(10));
        }
    }

    /// Productivity plan P3: nothing starts a mesh, a solve or a run on its own. A workflow with a result is edited in every form the window has
    /// that buffers an edit (the Geometry tab's element size, sizing, adaptive, absorbing layers and sweep of the mesh, and the sweep form of the
    /// Build page) and frames are drawn long enough for every trial mesh and check to run: the engine holds the same mesh, no job was submitted,
    /// no run was recorded, the result on show is still there, and nothing is busy. Run sweep and Compute are the only things that submit.
    #[test]
    fn editing_every_form_starts_nothing() {
        let mut a = app();
        let path = examples::list().into_iter().find(|e| e.path.to_string_lossy().contains("02_cantilever")).unwrap().path;
        a.new_draft("Fresh", crate::shell::Start::Example(path));
        wait_idle(&mut a, 30);
        a.shell_compute();
        wait_idle(&mut a, 120);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert!(!a.results.is_empty());
        {
            let w = a.shell.current_mut().expect("a workflow is open");
            w.sub = crate::shell::Sub::Main;
            // the Geometry tab has the form of the mesh (element size, sizing, adaptive, absorbing layers, sweep of the mesh); the Sweep form of the
            // Study tab comes last
            w.sweep.shown = false;
            w.tab = Tab::Geometry;
        }
        frames(&mut a, 3, 0.0);

        let jobs = a.engine.jobs().list().len();
        let mesh_key = a.engine.mesh_key();
        let model = a.engine.model().cloned().expect("the solve left a mesh");
        let project = a.engine.project().clone().to_json();
        let runs = a.shell.current_mut().unwrap().runs.len();
        let results = a.results.len();

        // every buffered edit of the Geometry tab, one after the other, each followed by frames that let the trial mesh and the check run
        type Edit = Box<dyn Fn(&mut App)>;
        let edits: Vec<Edit> = vec![
            Box::new(|a| a.mesh_edit.h *= 0.5),
            Box::new(|a| a.mesh_edit.physics_sizing = !a.mesh_edit.physics_sizing),
            Box::new(|a| a.mesh_edit.adapt_scale = Some(0.8)),
            Box::new(|a| a.mesh_edit.structured_pml = !a.mesh_edit.structured_pml),
            Box::new(|a| a.mesh_edit.sweep = Some(vec![nyx_engine::SweepLayer { thickness: 1.0, n: 4, domain_offset: 0, fill: None, regions: None }])),
            Box::new(|a| a.mesh_edit.h = -1.0),
            Box::new(|a| a.mesh_edit.h = 0.02),
        ];
        for edit in &edits {
            edit(&mut a);
            a.geometry_dirty = true;
            let ctx = a.ctx.clone();
            a.edits_changed(&ctx);
            frames(&mut a, 3, crate::feedback::REST.as_secs_f64() + 0.15);
            assert_eq!(a.engine.jobs().list().len(), jobs, "an edit submitted a job");
            assert_eq!(a.engine.mesh_key(), mesh_key, "an edit that was not applied moved the mesh key");
            assert!(a.engine.model().is_some_and(|m| Arc::ptr_eq(m, &model)), "an edit replaced the mesh");
            assert_eq!(a.engine.project().clone().to_json(), project, "an edit reached the project without being applied");
            assert!(a.busy.is_none() && a.job.is_none(), "an edit made the window busy");
            assert_eq!(a.results.len(), results, "an edit cleared the result on show");
        }

        // the frames really drew the Geometry tab and the trial mesh really ran (otherwise "nothing started" would prove nothing)
        for _ in 0..40 {
            if matches!(a.feedback.size, crate::feedback::SizeView::Ready(_)) {
                break;
            }
            frames(&mut a, 1, 0.05);
        }
        assert!(matches!(a.feedback.size, crate::feedback::SizeView::Ready(_)), "the last edit's trial mesh did not come back: {:?}", a.feedback.size);

        // the sweep form of the Build page: its values and its count changed, the plan redrawn; only Run sweep submits
        {
            let w = a.shell.current_mut().unwrap();
            w.tab = Tab::Study;
            w.sweep.shown = true;
            w.sweep.n = 9;
            w.sweep.from = "1e-3".into();
            w.sweep.to = "5e-3".into();
        }
        frames(&mut a, 5, 0.1);
        assert!(a.sweep.is_none(), "the sweep form started a sweep");
        assert_eq!(a.engine.jobs().list().len(), jobs);
        assert_eq!(a.shell.current_mut().unwrap().runs.len(), runs, "no run was recorded");
    }

    /// The trial mesh of an edit comes back to the window, says how many cells the edit would mesh to, and is dropped when the edit moved on.
    #[test]
    fn the_size_of_an_edit_arrives_from_a_trial_mesh() {
        let mut a = app();
        a.load_example(example("02_cantilever"));
        wait_idle(&mut a, 30);
        a.mesh_edit.h *= 0.6;
        a.geometry_dirty = true;
        let ctx = a.ctx.clone();
        a.edits_changed(&ctx);
        assert!(a.feedback.problem.is_none());
        assert!(matches!(a.feedback.size, crate::feedback::SizeView::Pending));
        a.start_size_preview();
        let t0 = Instant::now();
        while matches!(a.feedback.size, crate::feedback::SizeView::Pending) {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(10), "no answer from the trial mesh");
            std::thread::sleep(Duration::from_millis(2));
        }
        let crate::feedback::SizeView::Ready(p) = a.feedback.size.clone() else { panic!("{:?}", a.feedback.size) };
        assert!(p.cells > a.doc.mesh.as_ref().unwrap().n_cells() as u64, "a finer h means more cells: {p:?}");
        assert!(crate::feedback::size_line(&p).contains("cells"));
        assert!(crate::feedback::cost_line(&p).contains("No cost estimate"), "no solve history in a test: {}", crate::feedback::cost_line(&p));

        // an answer for an edit that has been replaced is dropped
        let old = a.feedback.generation;
        a.mesh_edit.h *= 0.9;
        let ctx = a.ctx.clone();
        a.edits_changed(&ctx);
        a.take_size(old, Ok(p));
        assert!(matches!(a.feedback.size, crate::feedback::SizeView::Pending), "{:?}", a.feedback.size);

        // a problem is shown at once, and no trial mesh is made of it
        a.mesh_edit.h = -1.0;
        a.edits_changed(&ctx);
        assert!(a.feedback.problem.as_deref().is_some_and(|p| p.contains("element size")), "{:?}", a.feedback.problem);
    }

    /// The polarization of a new result is computed beside the frame and lands in the result view that holds that very result.
    #[test]
    fn the_polarization_arrives_after_the_result_and_only_for_the_result_it_was_made_from() {
        let mut a = app();
        a.load_example(example("06_wr90"));
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 120);
        let t0 = Instant::now();
        while !a.has_analysis() {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(30), "the polarization never arrived: {:?}", a.status);
            std::thread::sleep(Duration::from_millis(5));
        }
        let want = a.engine.polarization("ems").unwrap().csv();
        let got = a.results[a.current.unwrap()].analysis.as_ref().unwrap().csv();
        assert_eq!(got, want, "the window shows exactly what the engine command computes");

        // an analysis that was made for another result is not put on this one
        let other = a.results[a.current.unwrap()].result.clone();
        a.results[a.current.unwrap()].analysis = None;
        let mut stale = (*other).clone();
        stale.seconds += 1.0;
        a.take_polarization("ems", &Arc::new(stale), a.engine.polarization("ems").map_err(|e| e.to_string()));
        assert!(!a.has_analysis(), "an analysis of another result must not be taken");
    }

    /// Plan P3: the frame that takes a solve's result shows it and does nothing else that waits on a disk. The run record is written on a later
    /// frame, one piece per frame (the record, then the dashboard refresh), the result on show is the one that is recorded, and a second solve
    /// taken before the first record is written cannot be recorded under the first one's context.
    #[test]
    fn the_run_record_is_written_on_a_frame_after_the_result() {
        let mut a = app();
        let root = std::env::temp_dir().join(format!("nyx-deferred-record-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        a.attach_workspace(root.clone());
        let path = examples::list().into_iter().find(|e| e.path.to_string_lossy().contains("02_cantilever")).unwrap().path;
        a.new_draft("Fresh", crate::shell::Start::Example(path));
        wait_idle(&mut a, 30);
        a.shell_compute();
        let t0 = Instant::now();
        while a.results.is_empty() {
            a.poll();
            assert!(t0.elapsed() < Duration::from_secs(60), "no result");
            std::thread::sleep(Duration::from_millis(5));
        }
        // the frame that took the result queued the record and did not write it
        assert!(matches!(a.deferred.front(), Some(crate::shell::Deferred::Record { .. })), "the record waits for a later frame (queued: {}, run ctx kept: {}, error {:?}, status {})", a.deferred.len(), a.shell.run.is_some(), a.error, a.status);
        assert!(a.shell.live_run.is_none() && a.shell.run.is_none(), "its context was taken with the result");
        let runs = a.shell.current_mut().unwrap().runs.len();
        assert_eq!(runs, 0, "nothing written yet");
        // the next frames: the record, then the refresh of the dashboard (queued by the record), one piece each
        a.poll();
        assert!(a.shell.live_run.is_some(), "the record was written by the next poll");
        assert!(matches!(a.deferred.front(), Some(crate::shell::Deferred::RefreshEntries)), "the refresh is the next piece");
        a.poll();
        assert!(a.deferred.is_empty());
        wait_idle(&mut a, 30);
        assert_eq!(a.shell.current_mut().unwrap().runs.len(), 1, "exactly one run");
        assert!(a.error.is_none(), "{:?}", a.error);
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn picking_a_boundary_on_the_mesh_gives_a_selector() {
        let mut a = app();
        a.load_example(example("02_cantilever"));
        wait_idle(&mut a, 30);
        a.sync_rules();
        assert_eq!(a.bc_labels.len(), 2, "two conditions resolve on the mesh");
        assert!(a.rule_errors.is_empty(), "{:?}", a.rule_errors);
        let left = a.bc_labels[0].clone();
        let sugg = a.doc.model.as_ref().unwrap().suggest(&left);
        assert!(!sugg.is_empty(), "no selector for the left edge");
        for s in &sugg {
            assert_eq!(nyx_engine::facade::physics::LabelResolver::boundaries(a.doc.model.as_ref().unwrap().as_ref(), s).unwrap(), left, "{s}");
        }
    }

    // ------------------------------------------------------------- simulated-input E2E GUI tests

    const E2E_SCREEN: egui::Vec2 = egui::vec2(1456.0, 819.0);

    fn collect_shape_texts(clip: egui::Rect, shape: &egui::epaint::Shape, out: &mut Vec<(String, egui::Rect)>) {
        match shape {
            egui::epaint::Shape::Text(t) => {
                let r = t.visual_bounding_rect();
                if clip.intersects(r) && !t.galley.job.text.is_empty() {
                    out.push((t.galley.job.text.clone(), r));
                }
            }
            egui::epaint::Shape::Vec(v) => {
                for s in v {
                    collect_shape_texts(clip, s, out);
                }
            }
            _ => {}
        }
    }

    fn gui_frame_sized_scroll(a: &mut App, mut raw: egui::RawInput, screen: egui::Vec2, scroll_dy: f32) -> Vec<(String, egui::Rect)> {
        raw.screen_rect = Some(egui::Rect::from_min_size(egui::Pos2::ZERO, screen));
        a.poll();
        a.run_queued();
        if a.rules_dirty {
            a.sync_rules();
        }
        let ctx = a.ctx.clone();
        let mut out = ctx.run_ui(raw, |ui| {
            if scroll_dy != 0.0 {
                ui.ctx().input_mut(|i| {
                    i.smooth_scroll_delta.y += scroll_dy;
                });
            }
            let (f5, ctrl_shift_r, alt1, alt2, alt3, alt4, alt5, alt6) = ui.ctx().input_mut(|i| (
                i.consume_key(egui::Modifiers::NONE, egui::Key::F5),
                i.consume_key(egui::Modifiers::COMMAND | egui::Modifiers::SHIFT, egui::Key::R),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num1),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num2),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num3),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num4),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num5),
                i.consume_key(egui::Modifiers::ALT, egui::Key::Num6),
            ));
            if f5 {
                a.shell_compute();
            }
            if ctrl_shift_r {
                a.open_nys_script();
            }
            for (pressed, act) in [
                (alt1, crate::menu::Act::ShowGeometry),
                (alt2, crate::menu::Act::ShowMaterials),
                (alt3, crate::menu::Act::ShowPhysics),
                (alt4, crate::menu::Act::ShowStudy),
                (alt5, crate::menu::Act::ShowSweep),
                (alt6, crate::menu::Act::ShowResults),
            ] {
                if pressed {
                    let c = ui.ctx().clone();
                    a.do_act(&c, act);
                }
            }
            if a.shell.active() {
                egui::Panel::top("shell_menubar").show(ui, |ui| a.bar_ui(ui));
                a.shell_ui(ui);
            }
        });
        out.textures_delta.clear();
        let mut texts = Vec::new();
        for cs in &out.shapes {
            collect_shape_texts(cs.clip_rect, &cs.shape, &mut texts);
        }
        texts
    }

    fn gui_frame_sized(a: &mut App, raw: egui::RawInput, screen: egui::Vec2) -> Vec<(String, egui::Rect)> {
        gui_frame_sized_scroll(a, raw, screen, 0.0)
    }

    fn gui_frame(a: &mut App, raw: egui::RawInput) -> Vec<(String, egui::Rect)> {
        gui_frame_sized(a, raw, E2E_SCREEN)
    }

    fn gui_focusable_rects(a: &App) -> Vec<(egui::Id, egui::LayerId, egui::Rect)> {
        let candidates: Vec<(egui::Id, egui::LayerId, egui::Rect)> = a.ctx.viewport(|vp| {
            let mut v = Vec::new();
            for layer in vp.prev_pass.widgets.layer_ids() {
                for w in vp.prev_pass.widgets.get_layer(layer) {
                    if w.enabled && w.sense.is_focusable() {
                        v.push((w.id, layer, w.rect));
                    }
                }
            }
            v
        });
        candidates
            .into_iter()
            .filter(|(id, _, _)| egui::text_edit::TextEditState::load(&a.ctx, *id).is_some())
            .collect()
    }

    fn gui_find_text(texts: &[(String, egui::Rect)], needle: &str) -> Option<egui::Rect> {
        texts
            .iter()
            .rev()
            .find(|(t, _)| t == needle)
            .or_else(|| texts.iter().rev().find(|(t, _)| t.contains(needle)))
            .map(|(_, r)| *r)
    }

    fn gui_click_at(a: &mut App, pos: egui::Pos2) -> Vec<(String, egui::Rect)> {
        let mut down = egui::RawInput::default();
        down.events.push(egui::Event::PointerMoved(pos));
        down.events.push(egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(a, down);

        let mut up = egui::RawInput::default();
        up.events.push(egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(a, up);

        gui_frame(a, egui::RawInput::default());
        gui_frame(a, egui::RawInput::default())
    }

    fn gui_scroll_inspector(a: &mut App, dy: f32) -> Vec<(String, egui::Rect)> {
        let pos = egui::pos2(200.0, 500.0);
        let mut raw = egui::RawInput::default();
        raw.events.push(egui::Event::PointerMoved(pos));
        gui_frame_sized_scroll(a, raw, E2E_SCREEN, dy);
        gui_frame(a, egui::RawInput::default())
    }

    fn gui_click_text(a: &mut App, needle: &str) -> egui::Rect {
        let mut texts = gui_frame(a, egui::RawInput::default());
        if gui_find_text(&texts, needle).is_none() {
            texts = gui_scroll_inspector(a, -420.0);
        }
        if gui_find_text(&texts, needle).is_none() {
            texts = gui_scroll_inspector(a, 840.0);
        }
        let rect = gui_find_text(&texts, needle).unwrap_or_else(|| {
            let sample: Vec<&str> = texts.iter().map(|(t, _)| t.as_str()).collect();
            panic!("GUI text '{needle}' not found on screen; visible texts: {sample:?}");
        });
        gui_click_at(a, rect.center());
        rect
    }

    fn gui_click_combo_item(a: &mut App, item: &str) {
        let mut texts = gui_frame(a, egui::RawInput::default());
        if let Some(rect) = texts.iter().rev().find(|(t, _)| t == item).map(|(_, r)| *r) {
            gui_click_at(a, rect.center());
            return;
        }
        let popup_pos = gui_find_text(&texts, "air")
            .or_else(|| gui_find_text(&texts, "numbers"))
            .or_else(|| gui_find_text(&texts, "aluminium"))
            .unwrap_or_else(|| panic!("open material ComboBox popup not found while looking for '{item}'; texts={texts:?}"))
            .center();
        for _ in 0..8 {
            let mut raw = egui::RawInput::default();
            raw.events.push(egui::Event::PointerMoved(popup_pos));
            gui_frame_sized_scroll(a, raw, E2E_SCREEN, -200.0);
            texts = gui_frame(a, egui::RawInput::default());
            if let Some(rect) = texts.iter().rev().find(|(t, _)| t == item).map(|(_, r)| *r) {
                gui_click_at(a, rect.center());
                return;
            }
        }
        panic!("ComboBox item '{item}' not found after scrolling popup; texts={texts:?}");
    }

    fn gui_type_in_rect(a: &mut App, rect: egui::Rect, text: &str) {
        gui_click_at(a, rect.center());
        let mut raw = egui::RawInput::default();
        raw.events.push(egui::Event::Text(text.to_string()));
        gui_frame(a, raw);

        let mut enter = egui::RawInput::default();
        enter.events.push(egui::Event::Key {
            key: egui::Key::Enter,
            physical_key: Some(egui::Key::Enter),
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        enter.events.push(egui::Event::Key {
            key: egui::Key::Enter,
            physical_key: Some(egui::Key::Enter),
            pressed: false,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(a, enter);
    }

    fn gui_send_key(a: &mut App, key: egui::Key, modifiers: egui::Modifiers) -> Vec<(String, egui::Rect)> {
        let mut raw = egui::RawInput::default();
        raw.events.push(egui::Event::Key {
            key,
            physical_key: Some(key),
            pressed: true,
            repeat: false,
            modifiers,
        });
        gui_frame(a, raw);
        let mut up = egui::RawInput::default();
        up.events.push(egui::Event::Key {
            key,
            physical_key: Some(key),
            pressed: false,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(a, up)
    }

    fn populate_seven_saved_workflows(root: &std::path::Path) -> (Vec<String>, String) {
        use nyx_engine::facade::store::{Workspace, WorkflowDir};
        std::fs::create_dir_all(root).unwrap();
        let mut ids = Vec::new();
        let mut sms_3d_id = String::new();
        let specs: [(&str, &str, bool); 7] = [
            ("wf1_cap_electrostatics", "dfs", false),
            ("wf2_cavity_resonator", "ems", true),
            ("wf3_waveguide_ems", "ems", false),
            ("wf4_thermal_plate", "dfs", false),
            ("wf5_beam_vibration", "sms", false),
            ("wf6_rib_coupler", "ems", false),
            ("wf7_cantilever_solid_3d", "sms", true),
        ];
        for (i, (title, phys, three_d)) in specs.iter().enumerate() {
            let wid = format!("wf0{}", i + 1);
            let spec = crate::examples::blank_json(phys, *three_d).unwrap();
            WorkflowDir::create_with(root, Some(&wid), title, "gui", &spec, nyx_engine::store_link::card_of).unwrap();
            if *phys == "sms" && *three_d {
                sms_3d_id = wid.clone();
            }
            ids.push(wid);
        }
        let mut ws = Workspace::open(root).unwrap();
        ws.settings.open_tabs = ids.clone();
        ws.save_settings().unwrap();
        (ids, sms_3d_id)
    }

    #[test]
    fn e2e_gui_clicks_blank_to_soi_300nm() {
        let root = std::env::temp_dir().join(format!("nyx-e2e-soi300-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let (initial_tabs, _) = populate_seven_saved_workflows(&root);
        assert_eq!(initial_tabs.len(), 7);

        let mut a = app();
        wait_idle(&mut a, 30);
        a.attach_workspace(root.clone());
        assert!(a.shell.active());
        assert_eq!(a.shell.tabs.len(), 7, "all 7 pre-populated tabs opened");

        // Render at 1456 x 819 and verify Settings and Search buttons in top_bar stay within [0.0, 1456.0]
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let settings_rect = gui_find_text(&texts, "Settings").expect("Settings button in top_bar");
        let search_rect = gui_find_text(&texts, "Search or run anything").expect("Search button in top_bar");
        assert!(
            settings_rect.min.x >= 0.0 && settings_rect.max.x <= 1456.0,
            "Settings rect {settings_rect:?} pushed outside [0, 1456] with 7 open tabs"
        );
        assert!(
            search_rect.min.x >= 0.0 && search_rect.max.x <= 1456.0,
            "Search rect {search_rect:?} pushed outside [0, 1456] with 7 open tabs"
        );

        // Also check at 900 x 600 that Settings and Search remain inside [0.0, 900.0]
        let texts_900 = gui_frame_sized(&mut a, egui::RawInput::default(), egui::vec2(900.0, 600.0));
        let settings_900 = gui_find_text(&texts_900, "Settings").expect("Settings at 900x600");
        let search_900 = gui_find_text(&texts_900, "Search or run anything").expect("Search at 900x600");
        assert!(settings_900.min.x >= 0.0 && settings_900.max.x <= 900.0, "Settings at 900x600: {settings_900:?}");
        assert!(search_900.min.x >= 0.0 && search_900.max.x <= 900.0, "Search at 900x600: {search_900:?}");
        gui_frame(&mut a, egui::RawInput::default());

        // Click '+' in top_bar to open the New workflow sheet
        gui_click_text(&mut a, "+");
        assert_eq!(a.shell.sheet, crate::shell::Sheet::New);

        // Click "Blank experiment", select ems and 2D, type "soi300_blank_gui" into the title field, click "Create workflow"
        gui_click_text(&mut a, "Blank experiment");
        gui_click_text(&mut a, "Electromagnetic waves");
        gui_click_text(&mut a, "Electromagnetic waves");
        gui_click_text(&mut a, "2D section");
        let modal_inputs: Vec<egui::Rect> = gui_focusable_rects(&a)
            .into_iter()
            .filter(|(_, layer, _)| layer.order == egui::Order::Foreground)
            .map(|(_, _, r)| r)
            .collect();
        assert!(!modal_inputs.is_empty(), "expected title TextEdit in New workflow modal");
        gui_type_in_rect(&mut a, modal_inputs[0], "soi300_blank_gui");
        assert_eq!(a.shell.new.title, "soi300_blank_gui");

        gui_click_text(&mut a, "Create workflow");
        wait_idle(&mut a, 30);
        let wid = match &a.shell.page {
            crate::shell::Page::Workflow(id) => id.clone(),
            other => panic!("expected Workflow page, got {other:?}"),
        };
        assert_eq!(a.shell.open[&wid].title(), "soi300_blank_gui");
        assert_eq!(a.active().as_deref(), Some("ems"));

        // Click "+ parameter" and inspect frame-1 TextEdit rects & feedback layout stability
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let add_param_rect = gui_find_text(&texts, "+ parameter").expect("+ parameter button");
        // Click "+ parameter" (down + up -> frame 1 of the new parameter row!)
        let mut down = egui::RawInput::default();
        down.events.push(egui::Event::PointerMoved(add_param_rect.center()));
        down.events.push(egui::Event::PointerButton {
            pos: add_param_rect.center(),
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(&mut a, down);
        let mut up = egui::RawInput::default();
        up.events.push(egui::Event::PointerButton {
            pos: add_param_rect.center(),
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        });
        gui_frame(&mut a, up);

        // Frame 1 with the newly added parameter row rendered (after egui Grid's initial sizing pass):
        gui_frame(&mut a, egui::RawInput::default());
        let texts_f1 = gui_frame(&mut a, egui::RawInput::default());
        let params_hdr_f1 = gui_find_text(&texts_f1, "Parameters").expect("Parameters heading");
        let add_param_f1 = gui_find_text(&texts_f1, "+ parameter").expect("+ parameter button");
        let features_hdr_f1 = gui_find_text(&texts_f1, "Features").expect("Features heading");
        let mut param_boxes_f1: Vec<egui::Rect> = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .filter(|r| r.center().y > params_hdr_f1.max.y && r.center().y < add_param_f1.min.y && r.max.x < 700.0)
            .collect();
        param_boxes_f1.sort_by(|a, b| a.min.x.partial_cmp(&b.min.x).unwrap());
        assert_eq!(param_boxes_f1.len(), 2, "expected param_name and param_val boxes on frame 1: {param_boxes_f1:?}");
        assert!(param_boxes_f1[0].width() >= 72.0, "frame-1 param_name width={}", param_boxes_f1[0].width());
        assert!(param_boxes_f1[1].width() >= 104.0, "frame-1 param_val width={}", param_boxes_f1[1].width());

        // Verify Parameters/Features y-positions shift by 0.0 px while feedback_tick transitions from Pending to Ready
        assert!(matches!(a.feedback.size, crate::feedback::SizeView::Pending), "expected Pending feedback: {:?}", a.feedback.size);
        let y_params_pending = params_hdr_f1.min.y;
        let y_features_pending = features_hdr_f1.min.y;
        std::thread::sleep(crate::feedback::REST + Duration::from_millis(30));
        let t0 = Instant::now();
        while !matches!(a.feedback.size, crate::feedback::SizeView::Ready(_)) {
            gui_frame(&mut a, egui::RawInput::default());
            assert!(t0.elapsed() < Duration::from_secs(10), "feedback never became Ready: {:?}", a.feedback.size);
            std::thread::sleep(Duration::from_millis(10));
        }
        let texts_ready = gui_frame(&mut a, egui::RawInput::default());
        let y_params_ready = gui_find_text(&texts_ready, "Parameters").unwrap().min.y;
        let y_features_ready = gui_find_text(&texts_ready, "Features").unwrap().min.y;
        assert_eq!(y_params_ready - y_params_pending, 0.0, "Parameters shifted vertically during feedback transition");
        assert_eq!(y_features_ready - y_features_pending, 0.0, "Features shifted vertically during feedback transition");

        // Set element size h (m) to "220 nm"
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let h_label = gui_find_text(&texts, "element size h (m)").expect("element size h (m) label");
        let h_box = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .find(|r| (r.center().y - h_label.center().y).abs() < 12.0 && r.min.x >= h_label.max.x - 4.0)
            .expect("mesh_h TextEdit");
        gui_type_in_rect(&mut a, h_box, "220 nm");
        assert!((a.mesh_edit.h - 2.2e-7).abs() < 1e-18, "mesh_edit.h = {}", a.mesh_edit.h);

        // Enable physics_sizing + adaptive mesh (0.85)
        gui_click_text(&mut a, "size the mesh from the physics");
        assert!(a.mesh_edit.physics_sizing);
        gui_click_text(&mut a, "adaptive mesh");
        assert_eq!(a.mesh_edit.adapt_scale, Some(1.0));
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let sf_label = gui_find_text(&texts, "size factor").expect("size factor label");
        let sf_box = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .find(|r| (r.center().y - sf_label.center().y).abs() < 12.0 && r.min.x >= sf_label.max.x - 4.0)
            .expect("adapt_s TextEdit");
        gui_type_in_rect(&mut a, sf_box, "0.85");
        assert_eq!(a.mesh_edit.adapt_scale, Some(0.85));

        // Add 2nd parameter and set w = 500 nm (5e-7) and h_si = 300 nm (3e-7)
        gui_click_text(&mut a, "+ parameter");
        assert_eq!(a.geom_edit.params.len(), 2);
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let params_hdr = gui_find_text(&texts, "Parameters").unwrap();
        let add_param = gui_find_text(&texts, "+ parameter").unwrap();
        let mut pboxes: Vec<egui::Rect> = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .filter(|r| r.center().y > params_hdr.max.y && r.center().y < add_param.min.y && r.max.x < 700.0)
            .collect();
        pboxes.sort_by_key(|a| (a.center().y as i32 / 14, a.min.x as i32));
        assert_eq!(pboxes.len(), 4, "expected 4 parameter boxes (2 rows x 2 cols): {pboxes:?}");
        gui_type_in_rect(&mut a, pboxes[0], "w");
        gui_type_in_rect(&mut a, pboxes[1], "500 nm");
        gui_type_in_rect(&mut a, pboxes[2], "h_si");
        gui_type_in_rect(&mut a, pboxes[3], "300 nm");
        assert_eq!(a.geom_edit.params[0].name, "w");
        assert_eq!(a.geom_edit.params[0].value, nyx_engine::facade::geom::expr::Expr::Num(5e-7));
        assert_eq!(a.geom_edit.params[1].name, "h_si");
        assert_eq!(a.geom_edit.params[1].value, nyx_engine::facade::geom::expr::Expr::Num(3e-7));

        // Scroll the Geometry inspector down so the Features section is in the centre of the viewport
        gui_scroll_inspector(&mut a, -420.0);

        // Remove the default 'regions' overlay and 'plate_box' rectangle via GUI clicks
        // (clicking the bottom-most "remove" button first removes 'regions', then 'plate_box')
        for _ in 0..2 {
            let texts = gui_frame(&mut a, egui::RawInput::default());
            let rm_rect = texts
                .iter()
                .filter(|(t, _)| t == "remove")
                .map(|(_, r)| *r)
                .max_by(|a, b| a.center().y.partial_cmp(&b.center().y).unwrap())
                .expect("remove feature button");
            gui_click_at(&mut a, rm_rect.center());
        }
        assert!(a.geom_edit.features.is_empty(), "default features removed");

        // Build the 3 rectangles: r_box, r_clad, r_si via GUI clicks/edits
        let rect_specs = [
            ("r_box", ["-2e-6", "-2e-6"], ["2e-6", "0"]),
            ("r_clad", ["-2e-6", "0"], ["2e-6", "2e-6"]),
            ("r_si", ["-w/2", "0"], ["w/2", "h_si"]),
        ];
        for (fname, min_xz, max_xz) in rect_specs {
            gui_click_text(&mut a, "+ add feature");
            gui_click_text(&mut a, "Rect");

            let texts = gui_frame(&mut a, egui::RawInput::default());
            let active_chk = texts
                .iter()
                .filter(|(t, _)| t == "active")
                .map(|(_, r)| *r)
                .max_by(|a, b| a.center().y.partial_cmp(&b.center().y).unwrap())
                .expect("active checkbox in Rect feature");
            let min_lbl = texts
                .iter()
                .filter(|(t, _)| t == "min (x, z)")
                .map(|(_, r)| *r)
                .max_by(|a, b| a.center().y.partial_cmp(&b.center().y).unwrap())
                .expect("min (x, z) label");
            let max_lbl = texts
                .iter()
                .filter(|(t, _)| t == "max (x, z)")
                .map(|(_, r)| *r)
                .max_by(|a, b| a.center().y.partial_cmp(&b.center().y).unwrap())
                .expect("max (x, z) label");

            let fboxes: Vec<egui::Rect> = gui_focusable_rects(&a).into_iter().map(|(_, _, r)| r).filter(|r| r.max.x < 700.0).collect();
            let name_box = fboxes
                .iter()
                .copied()
                .find(|r| (r.center().y - active_chk.center().y).abs() < 12.0 && r.max.x <= active_chk.min.x + 8.0)
                .expect("feat_name box");
            let mut min_boxes: Vec<egui::Rect> = fboxes
                .iter()
                .copied()
                .filter(|r| (r.center().y - min_lbl.center().y).abs() < 12.0 && r.min.x >= min_lbl.max.x - 4.0)
                .collect();
            min_boxes.sort_by(|a, b| a.min.x.partial_cmp(&b.min.x).unwrap());
            let mut max_boxes: Vec<egui::Rect> = fboxes
                .iter()
                .copied()
                .filter(|r| (r.center().y - max_lbl.center().y).abs() < 12.0 && r.min.x >= max_lbl.max.x - 4.0)
                .collect();
            max_boxes.sort_by(|a, b| a.min.x.partial_cmp(&b.min.x).unwrap());
            assert_eq!(min_boxes.len(), 2, "expected 2 min (x, z) boxes for {fname}");
            assert_eq!(max_boxes.len(), 2, "expected 2 max (x, z) boxes for {fname}");

            gui_type_in_rect(&mut a, name_box, fname);
            gui_type_in_rect(&mut a, min_boxes[0], min_xz[0]);
            gui_type_in_rect(&mut a, min_boxes[1], min_xz[1]);
            gui_type_in_rect(&mut a, max_boxes[0], max_xz[0]);
            gui_type_in_rect(&mut a, max_boxes[1], max_xz[1]);

            // Collapse this Rect feature header so the inspector stays compact for the next feature
            let hdr_needle = format!("{fname}  ·  Rect");
            gui_click_text(&mut a, &hdr_needle);
        }

        // Add Overlay feature and name the 3 regions BOX, Clad, Si
        gui_scroll_inspector(&mut a, -500.0);
        gui_click_text(&mut a, "+ add feature");
        let mut texts = gui_frame(&mut a, egui::RawInput::default());
        if !texts.iter().any(|(t, _)| t == "Overlay") {
            let popup_pos = gui_find_text(&texts, "Rect")
                .or_else(|| gui_find_text(&texts, "Profiles"))
                .expect("ComboBox popup open with Profiles/Rect")
                .center();
            let mut raw = egui::RawInput::default();
            raw.events.push(egui::Event::PointerMoved(popup_pos));
            gui_frame_sized_scroll(&mut a, raw, E2E_SCREEN, -260.0);
            texts = gui_frame(&mut a, egui::RawInput::default());
        }
        let ov_rect = texts
            .iter()
            .rev()
            .find(|(t, _)| t == "Overlay")
            .map(|(_, r)| *r)
            .unwrap_or_else(|| panic!("exact 'Overlay' item not found; texts={texts:?}"));
        gui_click_at(&mut a, ov_rect.center());
        assert_eq!(a.geom_edit.features.len(), 4, "expected 3 Rects + 1 Overlay");
        gui_scroll_inspector(&mut a, -400.0);
        let mut all_boxes: Vec<egui::Rect> = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .filter(|r| r.max.x < 700.0)
            .collect();
        all_boxes.sort_by(|a, b| a.center().y.partial_cmp(&b.center().y).unwrap());
        assert!(all_boxes.len() >= 3, "expected at least 3 region name boxes in Overlay: {all_boxes:?}");
        let rboxes = &all_boxes[all_boxes.len() - 3..];
        gui_type_in_rect(&mut a, rboxes[0], "BOX");
        gui_type_in_rect(&mut a, rboxes[1], "Clad");
        gui_type_in_rect(&mut a, rboxes[2], "Si");

        // Click "Build mesh"
        gui_click_text(&mut a, "Build mesh");
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "mesh build error: {:?}", a.error);
        assert!(a.doc.mesh.is_some());

        // Click "Materials" pill and assign silica / silica / silicon
        gui_click_text(&mut a, "Materials");
        assert_eq!(a.shell.open[&wid].tab, Tab::Materials);
        // Assign "silica" to all regions via "apply material...", then override "Si" to "silicon"
        gui_click_text(&mut a, "apply material...");
        gui_click_combo_item(&mut a, "silica");
        let mut texts = gui_frame(&mut a, egui::RawInput::default());
        if !texts.iter().any(|(t, r)| t == "Si" && r.max.x < 700.0) {
            texts = gui_scroll_inspector(&mut a, -350.0);
        }
        let si_label = texts
            .iter()
            .filter(|(t, r)| t == "Si" && r.max.x < 700.0)
            .map(|(_, r)| *r)
            .next()
            .expect("Si region label in Materials tab");
        let si_combo = texts
            .iter()
            .filter(|(t, r)| t == "silica" && (r.center().y - si_label.center().y).abs() < 12.0 && r.min.x > si_label.max.x)
            .map(|(_, r)| *r)
            .next()
            .expect("Si material combo");
        gui_click_at(&mut a, si_combo.center());
        gui_click_combo_item(&mut a, "silicon");

        // Click the "Study" tab of the left panel
        gui_click_text(&mut a, "Study");
        assert_eq!(a.shell.open[&wid].tab, Tab::Study);

        // Set wavelength 1550 nm (freq ≈ 193.414e12), n_target = 2.85, n_modes = 4, P2
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let lam_lbl = gui_find_text(&texts, "wavelength λ₀ (vacuum) [m]").expect("wavelength helper label");
        let lam_box = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .find(|r| (r.center().y - lam_lbl.center().y).abs() < 12.0 && r.min.x >= lam_lbl.max.x - 4.0)
            .expect("wavelength lambda0 TextEdit");
        gui_type_in_rect(&mut a, lam_box, "1550 nm");
        let freq = a.engine.project().spec("ems").unwrap().study.params["freq"].as_f64().unwrap();
        assert!((freq - 299_792_458.0 / 1.55e-6).abs() < 1e6, "freq={freq}");

        let texts = gui_frame(&mut a, egui::RawInput::default());
        let nt_lbl = gui_find_text(&texts, "Effective index near").expect("Effective index near label");
        let nt_box = gui_focusable_rects(&a)
            .into_iter()
            .map(|(_, _, r)| r)
            .find(|r| (r.center().y - nt_lbl.center().y).abs() < 12.0 && r.min.x >= nt_lbl.max.x - 4.0)
            .expect("n_target TextEdit");
        gui_type_in_rect(&mut a, nt_box, "2.85");
        assert_eq!(a.engine.project().spec("ems").unwrap().study.params["n_target"].as_f64(), Some(2.85));

        // Number of modes starts at default 4 in ems modes study; also click P2
        gui_click_text(&mut a, "P2");
        assert_eq!(a.engine.project().spec("ems").unwrap().order, 2);
        let n_modes = a.engine.project().spec("ems").unwrap().study.params.get("n_modes").and_then(Value::as_i64).unwrap_or(4);
        assert_eq!(n_modes, 4);

        // Press F5 and wait for solve
        gui_send_key(&mut a, egui::Key::F5, egui::Modifiers::NONE);
        assert!(a.busy.is_some(), "F5 started solve: {} / {:?}", a.status, a.error);
        wait_idle(&mut a, 120);
        frames(&mut a, 4, 0.05);
        assert!(a.error.is_none(), "solve error: {:?}", a.error);
        assert_eq!(a.results.len(), 1);
        let sols = &a.results[0].sols;
        assert_eq!(sols.len(), 4, "expected 4 modes");
        let te0 = sols[0].neff.expect("TE0 neff").0;
        let tm0 = sols[1].neff.expect("TM0 neff").0;
        assert!((te0 - 2.67175).abs() < 1e-4, "TE0 n_eff = {te0}, expected 2.67175 ± 1e-4");
        assert!((tm0 - 2.33382).abs() < 1e-4, "TM0 n_eff = {tm0}, expected 2.33382 ± 1e-4");

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn e2e_gui_draft_tab_switch_and_ctrl_w_isolation() {
        let root = std::env::temp_dir().join(format!("nyx-e2e-draft-iso-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let (_, sms_id) = populate_seven_saved_workflows(&root);
        {
            let mut ws = nyx_engine::facade::store::Workspace::open(&root).unwrap();
            ws.settings.open_tabs = vec![sms_id.clone()];
            ws.save_settings().unwrap();
        }

        let mut a = app();
        wait_idle(&mut a, 30);
        a.attach_workspace(root.clone());

        // Open the saved 3D sms workflow tab and record its spec on disk
        a.show_workflow(&sms_id);
        wait_idle(&mut a, 30);
        assert_eq!(a.active().as_deref(), Some("sms"));
        assert_eq!(a.doc.mesh.as_ref().map(|m| m.dim), Some(3));
        let sms_spec_path = a.shell.open[&sms_id].dir.as_ref().unwrap().spec_path("e1");
        let sms_spec_before = std::fs::read_to_string(&sms_spec_path).unwrap();

        // Create an unsaved 2D ems draft via the GUI ('+' in top_bar -> Blank experiment -> Create workflow)
        let texts = gui_frame(&mut a, egui::RawInput::default());
        let top_plus = texts
            .iter()
            .filter(|(t, r)| t == "+" && r.center().y < 85.0)
            .map(|(_, r)| *r)
            .next()
            .expect("'+' button in top_bar");
        gui_click_at(&mut a, top_plus.center());
        assert_eq!(a.shell.sheet, crate::shell::Sheet::New);
        gui_click_text(&mut a, "Blank experiment");
        gui_click_text(&mut a, "2D section");
        let modal_inputs: Vec<egui::Rect> = gui_focusable_rects(&a)
            .into_iter()
            .filter(|(_, layer, _)| layer.order == egui::Order::Foreground)
            .map(|(_, _, r)| r)
            .collect();
        gui_type_in_rect(&mut a, modal_inputs[0], "Draft2DEms");
        gui_click_text(&mut a, "Create workflow");
        wait_idle(&mut a, 30);

        let draft_id = match &a.shell.page {
            crate::shell::Page::Workflow(id) => id.clone(),
            other => panic!("expected draft workflow page, got {other:?}"),
        };
        assert!(a.shell.open[&draft_id].dir.is_none(), "draft has no folder on disk yet");
        assert_eq!(a.active().as_deref(), Some("ems"));
        assert_eq!(a.doc.mesh.as_ref().map(|m| m.dim), Some(2));

        // Click the saved 3D sms tab in top_bar
        gui_click_text(&mut a, "wf7_cantilever_solid_…");
        wait_idle(&mut a, 30);
        assert_eq!(a.shell.page, crate::shell::Page::Workflow(sms_id.clone()));
        assert_eq!(a.shell.loaded, Some((sms_id.clone(), "e1".into())));
        assert_eq!(a.active().as_deref(), Some("sms"));
        assert_eq!(a.doc.mesh.as_ref().map(|m| m.dim), Some(3));

        // Click back to the unsaved 2D ems draft tab in top_bar
        gui_click_text(&mut a, "Draft2DEms");
        wait_idle(&mut a, 30);
        assert_eq!(a.shell.page, crate::shell::Page::Workflow(draft_id.clone()));
        assert_eq!(a.shell.loaded, Some((draft_id.clone(), "e1".into())));
        assert_eq!(a.active().as_deref(), Some("ems"), "2D ems draft must not be overwritten by 3D sms model");
        assert_eq!(a.doc.mesh.as_ref().map(|m| m.dim), Some(2), "doc must not hold the 3D mesh after switching back to 2D draft");
        assert!(a.shell.open[&draft_id].dir.is_none(), "switching back without edits keeps draft unsaved");

        // Modify the draft via GUI click ("+ parameter") and flush: creates its own folder while leaving 3D sms untouched
        gui_click_text(&mut a, "+ parameter");
        a.flush();
        assert!(a.shell.open[&draft_id].dir.is_some(), "modifying the draft creates its own folder");
        let sms_spec_after = std::fs::read_to_string(&sms_spec_path).unwrap();
        assert_eq!(sms_spec_after, sms_spec_before, "saved 3D sms workflow on disk must remain completely untouched");

        // Send Ctrl+W and assert the newly active tab immediately reloads its own project and mesh
        gui_send_key(&mut a, egui::Key::W, egui::Modifiers::COMMAND);
        wait_idle(&mut a, 30);
        assert!(!a.shell.tabs.contains(&draft_id));
        let active_after_close = match &a.shell.page {
            crate::shell::Page::Workflow(id) => id.clone(),
            other => panic!("expected remaining workflow tab active after Ctrl+W, got {other:?}"),
        };
        assert_eq!(a.shell.loaded, Some((active_after_close.clone(), "e1".into())));
        assert!(a.doc.mesh.is_some(), "remaining tab reloaded its mesh after Ctrl+W");

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn e2e_gui_nys_script_sheet_check_run_and_param_units() {
        let root = std::env::temp_dir().join(format!("nyx-e2e-nys-sheet-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);

        let mut a = app();
        wait_idle(&mut a, 30);
        a.attach_workspace(root.clone());
        a.new_draft("soi300_script_wf", crate::shell::Start::Blank { physics: "ems".into(), three_d: false });
        wait_idle(&mut a, 30);
        let wid = match &a.shell.page {
            crate::shell::Page::Workflow(id) => id.clone(),
            other => panic!("expected workflow page, got {other:?}"),
        };

        // Test Ctrl+Shift+R opens the ScriptSheet modal, then close it
        gui_send_key(&mut a, egui::Key::R, egui::Modifiers::COMMAND | egui::Modifiers::SHIFT);
        assert!(a.nys.open, "Ctrl+Shift+R should open ScriptSheet");
        a.nys.open = false;

        // Click "Script" in the bar: the dock opens
        gui_click_text(&mut a, "Script");
        assert!(a.nys.open, "Clicking 'Script' in the bar opens the script dock");

        // Put the exact 15-command soi300_by_script.nys script into the script editor
        let script = "\
c0 = 299792458.0
f0 = c0 / 1.55e-6
nyx.project_new(title=\"soi300_by_script\", dim=2, h=2.2e-7)
nyx.set_param(name=\"w\", value=5.0e-7)
nyx.set_param(name=\"h_si\", value=3.0e-7)
nyx.add_rect(name=\"r_box\", min=[-2.0e-6, -2.0e-6], max=[2.0e-6, 0.0])
nyx.add_rect(name=\"r_clad\", min=[-2.0e-6, 0.0], max=[2.0e-6, 2.0e-6])
nyx.add_rect(name=\"r_si\", min=[\"-w/2\", \"0\"], max=[\"w/2\", \"h_si\"])
nyx.add_overlay(name=\"regions\", regions=[
    {\"name\": \"BOX\", \"input\": \"r_box\"},
    {\"name\": \"Clad\", \"input\": \"r_clad\"},
    {\"name\": \"Si\", \"input\": \"r_si\"},
])
nyx.add_physics(physics=\"ems\", order=2, study=\"modes\", study_params={\"freq\": f0, \"n_target\": 2.85, \"n_modes\": 4})
nyx.set_material(physics=\"ems\", region=\"BOX\", material=\"silica\")
nyx.set_material(physics=\"ems\", region=\"Clad\", material=\"silica\")
nyx.set_material(physics=\"ems\", region=\"Si\", material=\"silicon\")
nyx.set_mesh(h=2.2e-7, physics_sizing=True, adapt_scale=0.85)
nyx.mesh_build()
nyx.solve_run()
nyx.result_info(physics=\"ems\")
";
        a.nys.code = script.to_string();
        gui_frame(&mut a, egui::RawInput::default());

        // Click "Check" and assert check_msg is Ok
        gui_click_text(&mut a, "Check");
        assert!(
            matches!(&a.nys.check_msg, Some(Ok(msg)) if msg.contains("Static check passed")),
            "expected static check Ok, got {:?}",
            a.nys.check_msg
        );

        // Click "Run script" and poll until done
        gui_click_text(&mut a, "Run script");
        let t0 = Instant::now();
        while a.nys.running {
            gui_frame(&mut a, egui::RawInput::default());
            assert!(t0.elapsed() < Duration::from_secs(120), "script run timed out");
            std::thread::sleep(Duration::from_millis(10));
        }
        wait_idle(&mut a, 60);
        gui_frame(&mut a, egui::RawInput::default());
        assert!(a.error.is_none(), "unexpected error after script run: {:?}", a.error);
        let out_text = a.nys.output.as_ref().expect("script output present").as_ref().expect("script run succeeded");
        assert!(out_text.contains("15 engine commands"), "expected 15 engine commands in output: {out_text}");

        // Assert all 4 modes (2.671750, 2.333823, 1.631230, 1.611151) and result freshness
        assert_eq!(a.results.len(), 1);
        assert_eq!(a.stale_reason(0), None, "script result must be current in self.engine");
        let sols = &a.results[0].sols;
        assert_eq!(sols.len(), 4);
        let expected_modes = [2.671750, 2.333823, 1.631230, 1.611151];
        for (i, want) in expected_modes.iter().enumerate() {
            let got = sols[i].neff.expect("mode neff").0;
            assert!((got - want).abs() < 1e-5, "mode {} neff: got {got:.6}, expected {want:.6}", i + 1);
        }

        // Assert run is recorded in w.runs with params["w"] == 5e-7 and params["h_si"] == 3e-7
        let w = &a.shell.open[&wid];
        assert_eq!(w.runs.len(), 1, "expected 1 recorded run in w.runs");
        assert_eq!(w.runs[0].params.get("w").copied(), Some(5e-7));
        assert_eq!(w.runs[0].params.get("h_si").copied(), Some(3e-7));
        assert!(!w.runs[0].params.contains_key("parameter_w"));

        // Assert fmt::key_label and fmt::value never render w as Watts ("500 nW" / "Parameter [W]")
        assert_eq!(crate::fmt::key_label("w"), "W");
        assert_eq!(crate::fmt::value("w", 5e-7), "500 nm");
        assert_ne!(crate::fmt::key_label("parameter_w"), "Parameter [W]");
        assert_eq!(crate::fmt::value("parameter_w", 5e-7), "500 nm");

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn open_path_opens_nys_and_nyx_files_without_running_scripts() {
        let dir = std::env::temp_dir().join(format!("nyx-open-ext-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let nys_path = dir.join("demo.nys");
        std::fs::write(&nys_path, "\u{feff}nyx.project_info()\n").unwrap();
        let nyx_path = dir.join("demo.nyx");
        std::fs::write(&nyx_path, b"not-a-zip").unwrap();

        let mut a = app();
        wait_idle(&mut a, 30);

        a.open_path(nys_path.clone());
        assert!(a.nys.open, "opening .nys should open ScriptSheet");
        assert!(!a.nys.running, "opening .nys must NEVER auto-execute the script");
        assert_eq!(a.nys.path.as_deref(), Some(nys_path.as_path()));
        assert_eq!(a.nys.code, "nyx.project_info()\n");
        assert!(a.error.is_none(), "{:?}", a.error);

        a.open_path(nyx_path.clone());
        assert!(!a.nys.open, "opening .nyx should close ScriptSheet / show Open sheet");
        assert_eq!(a.shell.sheet, crate::shell::Sheet::Open);
        assert!(a.shell.open_file.as_ref().is_some_and(|(p, _)| p == &nyx_path));

        let _ = std::fs::remove_dir_all(&dir);
    }
}
