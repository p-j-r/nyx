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

/// A message from a helper thread (the optional Python script) to the window.
pub enum Msg {
    Script(Result<String, String>),
}

pub struct App {
    pub engine: Engine,
    pub doc: Doc,
    pub tab: Tab,
    pub vp: Viewport,
    pub results: Vec<ResultView>,
    pub current: Option<usize>,
    pub show_results: bool,
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
    pub script: VecDeque<String>,
    pub settle: u32,
    pub shot_pending: Arc<AtomicBool>,
    pub shell: crate::shell::Shell,
    pub colormap: Colormap,
    /// the last convergence check: its answer, whether one is running, and the settings of the next
    pub convergence: Option<Value>,
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
            geom_ui: Default::default(),
            rule_errors: HashMap::new(),
            bc_labels: vec![],
            domain_names: vec![],
            rules_dirty: false,
            suggestions: vec![],
            bottom_open: true,
            new_dialog: None,
            script_output: None,
            script: VecDeque::new(),
            settle: 6,
            shot_pending: Arc::new(AtomicBool::new(false)),
            shell: crate::shell::Shell::detached(),
            colormap: Colormap::Viridis,
            convergence: None,
            exports: HashMap::new(),
            opening: None,
            conv_factor: 0.6,
            conv_tol: 1e-3,
            mcp_server: nyx_connect::find_server(),
            connect_status: vec![],
            connect_note: None,
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
        self.geometry_dirty = false;
        self.picking = None;
        self.picked.clear();
        self.suggestions.clear();
    }

    /// Hands the geometry and mesh settings being edited to the engine. False (and an error shown) when the engine refuses them.
    pub fn push_edits(&mut self) -> bool {
        if self.geom_edit != self.engine.project().geometry {
            let g = serde_json::to_value(&self.geom_edit).unwrap_or(Value::Null);
            if self.cmd("set_geometry", json!({ "geometry": g })).is_none() {
                return false;
            }
        }
        if self.mesh_edit != self.engine.project().mesh {
            let m = &self.mesh_edit;
            let args = json!({
                "h": m.h, "adapt_scale": m.adapt_scale, "physics_sizing": m.physics_sizing, "points_per_wavelength": m.points_per_wavelength,
                "structured_pml": m.structured_pml, "max_dev": m.max_dev, "sweep": m.sweep,
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

    /// The Queue panel: what the engine runs and what waits (the same entries `queue_list` gives an agent), with cancel for each, hold and release for a
    /// waiting job and pause for the queue (ADR 87). It calls the runner's methods that the `queue_*` commands call.
    pub fn queue_panel(&mut self, ctx: &egui::Context) {
        if !self.show_queue {
            return;
        }
        let entries: Vec<Value> = self
            .engine
            .queue_list(None, None)
            .into_iter()
            .filter(|e| matches!(e["state"].as_str(), Some("queued" | "running")))
            .collect();
        let paused = self.engine.jobs().paused();
        let mut open = true;
        let mut act: Vec<(&str, u64)> = vec![];
        egui::Window::new("Queue").open(&mut open).default_width(460.0).show(ctx, |ui| {
            if entries.is_empty() {
                ui.label("Nothing is running or waiting.");
            }
            egui::Grid::new("queue_grid").striped(true).num_columns(4).show(ui, |ui| {
                for e in &entries {
                    let id = e["id"].as_u64().unwrap_or(0);
                    let queued = e["state"] == "queued";
                    ui.label(e["label"].as_str().unwrap_or(""));
                    if queued {
                        let held = if e["held"] == true { ", held" } else { "" };
                        ui.label(format!("waiting, number {}{held}", e["position"].as_u64().unwrap_or(0)));
                    } else {
                        let stage = e["stage"].as_str().unwrap_or("");
                        let stage = if stage.is_empty() { "starting" } else { stage };
                        let cancelled = if e["cancel_requested"] == true { " (cancel requested)" } else { "" };
                        ui.label(format!("{stage}, {:.1} s{cancelled}", e["seconds"].as_f64().unwrap_or(0.0)));
                    }
                    ui.label(format!("{} · priority {}", e["owner"].as_str().unwrap_or(""), e["priority"]));
                    ui.horizontal(|ui| {
                        if ui.small_button("Cancel").clicked() {
                            act.push(("cancel", id));
                        }
                        if queued {
                            if e["held"] == true {
                                if ui.small_button("Release").clicked() {
                                    act.push(("release", id));
                                }
                            } else if ui.small_button("Hold").clicked() {
                                act.push(("hold", id));
                            }
                        }
                    });
                    ui.end_row();
                }
            });
            ui.horizontal(|ui| {
                if ui.button(if paused { "Resume the queue" } else { "Pause the queue" }).on_hover_text("A paused queue starts nothing new; what runs finishes").clicked() {
                    act.push(("pause", 0));
                }
            });
        });
        for (what, id) in act {
            match what {
                "cancel" => {
                    if self.job.as_ref().is_some_and(|j| j.0 == id) {
                        self.cancel_job();
                    } else {
                        self.cancel_waiting(id);
                    }
                }
                "hold" => {
                    self.engine.jobs().hold(id);
                }
                "release" => {
                    self.engine.jobs().release(id);
                }
                _ => {
                    if paused {
                        self.engine.jobs().resume();
                    } else {
                        self.engine.jobs().pause();
                    }
                }
            }
        }
        self.show_queue = open;
        if !entries.is_empty() {
            ctx.request_repaint_after(std::time::Duration::from_millis(250));
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

    /// One line above the result on show when it is stale, with a button to solve again.
    pub fn stale_band(&mut self, ui: &mut egui::Ui) {
        let Some(why) = self.current.and_then(|i| self.stale_reason(i)) else { return };
        ui.horizontal(|ui| {
            ui.colored_label(Color32::from_rgb(235, 170, 60), format!("⚠ Stale: {why}."));
            if self.busy.is_none() && ui.small_button("Solve again").clicked() {
                self.compute();
            }
        });
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

    pub fn open_path(&mut self, path: PathBuf) {
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
        if ext == "json" {
            self.load_project(&path);
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
        self.set_error(format!("Don't know how to open '.{ext}' files (nyx projects are .nyx.json)"));
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
        while let Ok((_, msg)) = self.rx.try_recv() {
            match msg {
                Msg::Script(Ok(t)) => self.script_output = Some(t),
                Msg::Script(Err(e)) => self.set_error(e),
            }
        }
        let finished = self.engine.poll();
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
        self.take_model();
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
            let mut rv = ResultView::new(name, Some(id.clone()), result.clone(), &mesh);
            if self.engine.polarization_applies(id) {
                rv.analysis = self.engine.polarization(id).ok();
            }
            if active.as_deref() == Some(id.as_str()) {
                live = Some(result);
            }
            self.push_result(rv);
        }
        if let Some(r) = live {
            self.shell.live_run = self.record_ok(&r);
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

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new().add_filter("nyx projects", &["json"]).pick_file() {
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
    fn refresh_connect(&mut self) {
        self.connect_status = match (Env::detect(), &self.mcp_server) {
            (Ok(env), Some(server)) => nyx_connect::list(&env, server).into_iter().filter(|h| h.status != Status::NotInstalled).collect(),
            _ => vec![],
        };
    }

    /// The Connect menu: one row per Claude app with a Connect/Disconnect button, so the AI interface (`nyx-mcp`)
    /// is registered without editing a config file.
    fn connect_menu(&mut self, ui: &mut egui::Ui) {
        ui.set_min_width(320.0);
        ui.label(RichText::new("Let Claude run nyx").strong());
        ui.weak("Registers the nyx MCP server with the app, so Claude can set up and solve simulations here.");
        ui.separator();
        let Some(server) = self.mcp_server.clone() else {
            ui.colored_label(Color32::from_rgb(200, 120, 20), "nyx-mcp was not found next to nyx.");
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
            ui.colored_label(Color32::from_rgb(200, 60, 60), "cannot find the home folder");
            return;
        };
        if self.connect_status.is_empty() {
            ui.weak("No supported AI app found (Claude Code, Claude Desktop, Cursor, VS Code, Windsurf). Use 'Copy config JSON' below.");
        }
        for host in self.connect_status.clone() {
            let status = host.status.clone();
            ui.horizontal(|ui| {
                ui.label(RichText::new(host.name).strong());
                match &status {
                    Status::Connected => ui.colored_label(Color32::from_rgb(40, 150, 70), "connected"),
                    Status::Stale { .. } => ui.colored_label(Color32::from_rgb(200, 120, 20), "points to another nyx-mcp"),
                    Status::Unreadable { .. } => ui.colored_label(Color32::from_rgb(200, 60, 60), "config unreadable"),
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
                ui.colored_label(Color32::from_rgb(200, 60, 60), why);
            }
            ui.add_space(4.0);
        }
        if let Some(n) = &self.connect_note {
            ui.separator();
            ui.label(n);
        }
        ui.separator();
        ui.weak(format!("Server: {}", server.display()));
        ui.horizontal(|ui| {
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
            let b = ui.add_enabled(self.can_queue_compute() && has_mesh, egui::Button::new(RichText::new(if self.busy.is_some() { "  Queue  " } else { "  Compute  " }).strong()).fill(Color32::from_rgb(38, 110, 70)));
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
                        ui.colored_label(Color32::from_rgb(235, 170, 60), format!("Cancel requested: waiting for the stage to end ({:.1} s so far; it stops at its next check)", st.stage_seconds));
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
                ui.colored_label(Color32::from_rgb(235, 110, 100), format!("⚠ {e}"));
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

    pub fn viewport_ui(&mut self, ui: &mut egui::Ui) {
        let bg = Color32::from_rgb(24, 26, 30);
        let mut bc_colors: HashMap<i32, Color32> = HashMap::new();
        for (i, labels) in self.bc_labels.iter().enumerate() {
            for &l in labels {
                bc_colors.insert(l, crate::colormap::palette(i + 3));
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
                    bg,
                )
            }
            (_, Some(view)) => self.vp.show(ui, Content::Mesh { view: &view, bc_colors: &bc_colors, selected: &selected, show_edges, show_numbers }, picking, bg),
            _ => self.vp.show(ui, Content::Empty, false, bg),
        };
        if let Some(why) = stale {
            // on the plot itself, so a stale field is never read as current; the button to solve again is in the Results tab and its panel
            let at = ui.max_rect().left_top() + egui::vec2(10.0, 10.0);
            let text = format!("⚠ STALE: {why}");
            let galley = ui.painter().layout_no_wrap(text, egui::FontId::proportional(14.0), Color32::from_rgb(40, 30, 10));
            let rect = egui::Rect::from_min_size(at, galley.size() + egui::vec2(12.0, 6.0));
            ui.painter().rect_filled(rect, 4.0, Color32::from_rgb(235, 170, 60));
            ui.painter().galley(rect.min + egui::vec2(6.0, 3.0), galley, Color32::BLACK);
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
        self.poll();
        self.run_queued();
        if self.rules_dirty {
            self.sync_rules();
        }
        let c = ctx.clone();
        self.script_step(&c);
        let dropped: Vec<PathBuf> = ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        if let Some(p) = dropped.into_iter().next() {
            self.open_path(p);
        }
        let (f5, esc) = ctx.input_mut(|i| (i.consume_key(egui::Modifiers::NONE, egui::Key::F5), i.key_pressed(egui::Key::Escape)));
        if f5 {
            self.shell_compute();
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
        self.queue_panel(ui.ctx());
        if self.shell.active() {
            self.shell_ui(ui);
            return;
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
        let t_vp = Instant::now();
        egui::CentralPanel::default().frame(egui::Frame::NONE).show(ui, |ui| self.viewport_ui(ui));
        if std::env::var_os("NYX_FRAME_LOG").is_some() {
            // diagnostics: time (s, wrapping every 1000) of the frame, our own layout and drawing code, whole viewport, state
            let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs_f64() % 1000.0).unwrap_or(0.0);
            eprintln!("FRAME t={now:.3} ui_ms={:.2} viewport_ms={:.2} busy={} show_results={}", t_ui.elapsed().as_secs_f64() * 1e3, t_vp.elapsed().as_secs_f64() * 1e3, self.busy.is_some(), self.show_results);
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

    /// Processes finished jobs until nothing runs (or panics after `secs`).
    fn wait_idle(a: &mut App, secs: u64) {
        let t0 = Instant::now();
        loop {
            a.poll();
            a.run_queued();
            if a.busy.is_none() && a.queued.is_none() && a.exports.is_empty() {
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
        let _ = ctx.run_ui(Default::default(), |ui| a.queue_panel(ui.ctx()));
        a.engine.jobs().resume();
        wait_idle(&mut a, 120);
        assert!(a.doc.mesh.is_some() && a.queued.is_none() && a.waiting.is_empty());
        assert!(!a.results.is_empty(), "the queued Compute ran when the first job ended and its result is on show");
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
        let _ = ctx.run_ui(Default::default(), |ui| {
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
                let _ = ctx.run_ui(Default::default(), |ui| {
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
            let _ = ctx.run_ui(Default::default(), |ui| {
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
        let _ = ctx.run_ui(Default::default(), |ui| a.inspector(ui));
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
}
