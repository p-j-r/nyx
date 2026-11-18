//! The shell around the simulation window: the tab bar (Dashboard, one tab per open workflow, +), the status chip, the search
//! palette entry, and the glue between the window's single live project and the workspace on disk (`nyx-store`).
//!
//! One project lives in memory at a time (`App::proj`, the mesh and the results). Each open workflow tab remembers which of its experiments
//! is loaded; switching tabs or experiments first stores the current spec (`autosave_tick`, `flush`) and then loads the other one from its
//! `experiment.nyx.spec`. Nothing is saved until something changes: a new workflow is a *draft* in memory and gets its folder on the first
//! edit, run, title or note. Every Compute becomes a run record (`record_ok`, `record_end`). Without a workspace (tests, `--batch`) the shell is
//! switched off and the window behaves as before.

use crate::app::App;
use crate::theme;
use eframe::egui::{self, RichText, Vec2};
use nyx_engine::facade::store::runs::{RunPatch, RunRecord, RunStatus};
use nyx_engine::facade::store::util::{hash_bytes, now_secs};
use nyx_engine::facade::store::{Entry, Settings, Workspace, WorkflowDir};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;
use std::time::{Duration, Instant};

#[derive(Clone, Debug, PartialEq)]
pub enum Page {
    Dashboard,
    Workflow(String),
    Settings,
}

/// Where a new workflow starts from: an example (or any project file), or a blank simulation.
pub enum Start {
    Example(PathBuf),
    /// an HDF5 results file: the project it holds with its mesh and results on show (ADR 62, 63)
    Results(PathBuf),
    Blank { physics: String, three_d: bool },
}

/// The three sections of an experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Part {
    Concept,
    Build,
    Analysis,
}

/// What fills the workflow view's centre besides the section itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sub {
    Main,
    Compare,
    Field,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sheet {
    None,
    New,
    Palette,
    Versions,
    Share,
    Open,
    /// the version, the licences of the fonts and the tagline
    About,
}

pub struct OpenWf {
    pub id: String,
    pub dir: Option<WorkflowDir>,
    /// draft state (no folder yet)
    pub draft_title: String,
    pub spec0: String,
    pub tags: Vec<String>,
    pub tag_input: String,
    pub exp: String,
    pub part: Part,
    pub sub: Sub,
    pub concept_tab: crate::state::Tab,
    pub build_tab: crate::state::Tab,
    pub runs: Vec<RunRecord>,
    pub runs_bad: usize,
    pub filter_text: String,
    pub selected: BTreeSet<String>,
    pub sort: String,
    pub desc: bool,
    pub hidden_cols: BTreeSet<String>,
    pub show_failed: bool,
    pub by_version: bool,
    pub plot_x: String,
    pub plot_y: String,
    pub last_hash: String,
    pub saved_at: Option<Instant>,
    pub notes_open: bool,
    pub notes: String,
    pub function: String,
    pub compare_metric: String,
    pub draft_buf: String,
    pub focus: Option<String>,
    pub card_for: String,
    pub card_summary: String,
    pub sweep: crate::sweepui::Form,
}

impl OpenWf {
    pub fn new(id: String) -> OpenWf {
        OpenWf {
            id,
            dir: None,
            draft_title: String::new(),
            spec0: String::new(),
            tags: vec![],
            tag_input: String::new(),
            exp: "e1".into(),
            part: Part::Concept,
            sub: Sub::Main,
            concept_tab: crate::state::Tab::Geometry,
            build_tab: crate::state::Tab::Physics,
            runs: vec![],
            runs_bad: 0,
            filter_text: String::new(),
            selected: BTreeSet::new(),
            sort: "run".into(),
            desc: true,
            hidden_cols: BTreeSet::new(),
            show_failed: false,
            by_version: false,
            plot_x: String::new(),
            plot_y: String::new(),
            last_hash: String::new(),
            saved_at: None,
            notes_open: false,
            notes: String::new(),
            function: "Metric vs parameter".into(),
            compare_metric: String::new(),
            draft_buf: String::new(),
            focus: None,
            card_for: String::new(),
            card_summary: String::new(),
            sweep: Default::default(),
        }
    }

    pub fn title(&self) -> String {
        match &self.dir {
            Some(d) => d.wf.title.clone(),
            None => if self.draft_title.trim().is_empty() { "New workflow".into() } else { self.draft_title.clone() },
        }
    }

    pub fn reload_runs(&mut self) {
        if let Some(d) = &self.dir {
            let (r, bad) = d.runs(&self.exp);
            self.runs = r;
            self.runs_bad = bad;
            self.selected.retain(|s| self.runs.iter().any(|r| &r.id == s));
        } else {
            self.runs.clear();
        }
    }
}

/// The run that is being solved, so its record can be written whatever the window shows when it ends.
/// The run record whose result is on show: where it lives, and the physics it holds the result of (a convergence verdict is written into it, P6).
#[derive(Clone, Debug)]
pub struct LiveRun {
    pub wf_id: String,
    pub exp: String,
    pub run: String,
    pub physics: String,
}

pub struct RunCtx {
    pub wf_id: String,
    pub exp: String,
    pub spec_json: String,
    pub params: BTreeMap<String, f64>,
    pub started: Instant,
}

pub struct NewState {
    pub title: String,
    pub template: usize,
    pub physics: String,
    pub three_d: bool,
    pub all_examples: bool,
}

impl Default for NewState {
    fn default() -> Self {
        NewState { title: String::new(), template: 0, physics: "ems".into(), three_d: false, all_examples: false }
    }
}

pub struct Shell {
    pub ws: Option<Workspace>,
    pub entries: Vec<Entry>,
    pub spark: HashMap<String, Vec<f64>>,
    pub page: Page,
    pub tabs: Vec<String>,
    pub open: HashMap<String, OpenWf>,
    /// (workflow id, experiment id) whose spec is in `App::proj`
    pub loaded: Option<(String, String)>,
    pub run: Option<RunCtx>,
    /// the id of the run whose fields are still in `App::results`
    pub live_run: Option<LiveRun>,
    pub sheet: Sheet,
    pub new: NewState,
    pub palette_q: String,
    pub dash_q: String,
    pub dash_physics: String,
    pub dash_status: String,
    pub dash_material: String,
    pub earlier_open: bool,
    pub version_sel: Option<u32>,
    pub version_name: String,
    pub share_level: nyx_engine::facade::store::ShareLevel,
    pub share_only_current: bool,
    pub share_msg: String,
    pub open_file: Option<(PathBuf, nyx_engine::facade::store::archive::Inspection)>,
    pub classic: bool,
    pub dark: bool,
    /// the zoom factor last seen: a change is Ctrl +/- (or the View menu) and is saved
    pub last_zoom: f32,
    /// the native menu bar (Windows)
    #[cfg(windows)]
    pub native_menu: Option<crate::menu::Native>,
    pub logo: Option<egui::TextureHandle>,
    pub last_tick: Instant,
    pub last_refresh: Instant,
    pub closed: Vec<String>,
    pub notice: Option<String>,
    /// the "cancel this run?" question is on screen
    pub confirm_cancel: bool,
    /// scripts left out of the file being shared
    pub share_exclude: Vec<String>,
}

impl Shell {
    /// No workspace: the shell is off (tests and `--batch`).
    pub fn detached() -> Shell {
        Shell {
            ws: None,
            entries: vec![],
            spark: HashMap::new(),
            page: Page::Dashboard,
            tabs: vec![],
            open: HashMap::new(),
            loaded: None,
            run: None,
            live_run: None,
            sheet: Sheet::None,
            new: NewState::default(),
            palette_q: String::new(),
            dash_q: String::new(),
            dash_physics: String::new(),
            dash_status: String::new(),
            dash_material: String::new(),
            earlier_open: true,
            version_sel: None,
            version_name: String::new(),
            share_level: nyx_engine::facade::store::ShareLevel::Recipe,
            share_only_current: false,
            share_msg: String::new(),
            open_file: None,
            classic: true,
            dark: false,
            last_zoom: 1.0,
            #[cfg(windows)]
            native_menu: None,
            logo: None,
            last_tick: Instant::now(),
            last_refresh: Instant::now() - Duration::from_secs(60),
            closed: vec![],
            notice: None,
            confirm_cancel: false,
            share_exclude: vec![],
        }
    }

    pub fn active(&self) -> bool {
        self.ws.is_some() && !self.classic
    }

    pub fn settings(&self) -> Settings {
        self.ws.as_ref().map(|w| w.settings.clone()).unwrap_or_default()
    }

    pub fn current_mut(&mut self) -> Option<&mut OpenWf> {
        match &self.page {
            Page::Workflow(id) => {
                let id = id.clone();
                self.open.get_mut(&id)
            }
            _ => None,
        }
    }
}

impl App {
    /// Switches the shell on with the workspace at `root` (the window passes the default one; `NYX_WORKSPACE` overrides it).
    pub fn attach_workspace(&mut self, root: PathBuf) {
        match Workspace::open(&root) {
            Ok(ws) => {
                let zoom = ws.settings.zoom.clamp(0.5, 3.0);
                self.shell.ws = Some(ws);
                self.shell.classic = false;
                self.ctx.set_zoom_factor(zoom);
                self.shell.last_zoom = zoom;
                self.sync_look(&self.ctx.clone());
                self.refresh_entries();
                let tabs = self.shell.settings().open_tabs;
                for id in tabs {
                    if self.shell.ws.as_ref().and_then(|w| w.find(&id)).is_some() {
                        self.open_tab(&id, false);
                    }
                }
                self.shell.page = Page::Dashboard;
            }
            Err(e) => self.set_error(format!("Workspace {}: {e}; using the classic layout", root.display())),
        }
    }

    pub fn refresh_entries(&mut self) {
        let Some(ws) = &self.shell.ws else { return };
        self.shell.entries = ws.list();
        let mut spark = HashMap::new();
        for e in &self.shell.entries {
            if let Some(x) = e.wf.experiments.last() {
                if x.runs.count > 1 {
                    if let Ok(d) = WorkflowDir::open_with(&e.dir, nyx_engine::store_link::card_of) {
                        spark.insert(e.wf.id.clone(), d.sparkline(&x.id));
                    }
                }
            }
        }
        self.shell.spark = spark;
        self.shell.last_refresh = Instant::now();
    }

    fn save_tabs(&mut self) {
        if let Some(ws) = &mut self.shell.ws {
            ws.settings.open_tabs = self.shell.tabs.iter().filter(|t| self.shell.open.get(*t).is_some_and(|o| o.dir.is_some())).cloned().collect();
            let _ = ws.save_settings();
        }
    }

    // ------------------------------------------------------------------------------------------------------ the live project

    /// The project as it would be saved (the geometry being edited included).
    pub fn spec_json(&self) -> Result<String, String> {
        let mut p = self.engine.project().clone();
        p.geometry = self.geom_edit.clone();
        p.mesh = self.mesh_edit.clone();
        Ok(p.to_json())
    }

    /// The engine has a new project (opened, or made): the window takes over its geometry and settings and meshes it.
    fn adopt_project(&mut self) {
        self.adopt();
        self.project_path = None;
        self.rebuild_mesh();
        self.tab = crate::state::Tab::Geometry;
    }

    /// Stores the live project in its workflow when it changed (called about once a second, and before anything that replaces it).
    pub fn flush(&mut self) {
        let Some((wid, eid)) = self.shell.loaded.clone() else { return };
        let Ok(spec) = self.spec_json() else { return };
        let h = hash_bytes(spec.as_bytes());
        let Some(w) = self.shell.open.get_mut(&wid) else { return };
        if h == w.last_hash {
            return;
        }
        if w.dir.is_none() {
            if spec == w.spec0 {
                w.last_hash = h;
                return;
            }
            if self.persist_draft(&wid, &spec).is_err() {
                return;
            }
        }
        let Some(w) = self.shell.open.get_mut(&wid) else { return };
        if let Some(d) = &mut w.dir {
            match d.commit_spec(&eid, &spec) {
                Ok(_) => {
                    w.last_hash = h;
                    w.saved_at = Some(Instant::now());
                }
                Err(e) => {
                    let e = format!("Could not save the workflow: {e}");
                    self.set_error(e);
                }
            }
        }
        self.refresh_entries();
    }

    pub fn persist_draft(&mut self, wid: &str, spec: &str) -> Result<(), String> {
        let Some(root) = self.shell.ws.as_ref().map(|w| w.root.clone()) else { return Err("no workspace".into()) };
        let Some(w) = self.shell.open.get_mut(wid) else { return Err("no such workflow".into()) };
        let title = if w.draft_title.trim().is_empty() { "Untitled workflow".to_string() } else { w.draft_title.trim().to_string() };
        // the very first version is the spec the draft started from; this change is then committed on top of it
        let first = if w.spec0.is_empty() { spec.to_string() } else { w.spec0.clone() };
        match WorkflowDir::create_with(&root, Some(wid), &title, "gui", &first, nyx_engine::store_link::card_of) {
            Ok(mut d) => {
                if !w.tags.is_empty() {
                    let _ = d.set_tags(w.tags.clone());
                }
                if !w.notes.is_empty() {
                    let _ = d.set_notes(&w.notes);
                }
                w.last_hash = hash_bytes(first.as_bytes());
                w.dir = Some(d);
                self.save_tabs();
                Ok(())
            }
            Err(e) => {
                self.set_error(format!("Could not create the workflow folder: {e}"));
                Err(e)
            }
        }
    }

    /// Loads experiment `eid` of workflow `wid` into the window (the current one is stored first).
    pub fn load_experiment(&mut self, wid: &str, eid: &str) {
        self.flush();
        let Some(w) = self.shell.open.get_mut(wid) else { return };
        w.exp = eid.to_string();
        w.sub = Sub::Main;
        w.reload_runs();
        w.selected.clear();
        let Some(d) = &w.dir else { return };
        let path = d.spec_path(eid);
        match self.engine.call("project_open", &serde_json::json!({ "path": path.to_str().unwrap_or_default() })) {
            Ok(_) => {
                self.adopt_project();
                let h = self.spec_json().map(|x| hash_bytes(x.as_bytes())).unwrap_or_default();
                if let Some(w) = self.shell.open.get_mut(wid) {
                    w.last_hash = h;
                    w.notes = w.dir.as_ref().map(|d| d.notes()).unwrap_or_default();
                }
                self.shell.loaded = Some((wid.to_string(), eid.to_string()));
                self.status = format!("Opened {}", path.display());
            }
            Err(e) => self.set_error(format!("{}: {e}", path.display())),
        }
    }

    /// Opens a tab for the workflow `id` of the workspace and (when `show`) switches to it.
    pub fn open_tab(&mut self, id: &str, show: bool) {
        if !self.shell.open.contains_key(id) {
            let Some(entry) = self.shell.ws.as_ref().and_then(|w| w.find(id)) else { return };
            match WorkflowDir::open_with(&entry.dir, nyx_engine::store_link::card_of) {
                Ok(mut d) => {
                    if d.digest_is_stale() {
                        for eid in d.wf.experiments.iter().map(|e| e.id.clone()).collect::<Vec<_>>() {
                            let _ = d.rebuild_digest(&eid);
                        }
                    }
                    let mut w = OpenWf::new(id.to_string());
                    w.exp = d.wf.experiments.last().map(|e| e.id.clone()).unwrap_or_else(|| "e1".into());
                    w.tags = d.wf.tags.clone();
                    w.dir = Some(d);
                    self.shell.open.insert(id.to_string(), w);
                }
                Err(e) => {
                    self.set_error(e);
                    return;
                }
            }
        }
        if !self.shell.tabs.iter().any(|t| t == id) {
            self.shell.tabs.push(id.to_string());
            self.save_tabs();
        }
        if show {
            self.show_workflow(id);
        }
    }

    pub fn show_workflow(&mut self, id: &str) {
        self.flush();
        self.shell.page = Page::Workflow(id.to_string());
        let eid = self.shell.open.get(id).map(|w| w.exp.clone()).unwrap_or_default();
        if self.shell.loaded.as_ref() != Some(&(id.to_string(), eid.clone())) {
            self.load_experiment(id, &eid);
        }
    }

    /// A new draft workflow from a project (nothing is written until the first change).
    pub fn new_draft(&mut self, title: &str, start: Start) {
        self.flush();
        if let Start::Results(path) = &start {
            // reading the file, meshing its project again and checking the mesh run on a job; the draft workflow is made when it ends (`finish_open`)
            match self.engine.start_open_results(path.to_str().unwrap_or_default()) {
                Ok(job) => {
                    self.opening = Some(title.to_string());
                    self.started(job, crate::app::JobKind::Open, format!("Opening {title}"), true);
                    self.status = format!("Opening {title}: reading the file, meshing its project and checking the mesh");
                }
                Err(e) => self.set_error(e.to_string()),
            }
            return;
        }
        let id = nyx_engine::facade::store::util::short_id();
        let mut w = OpenWf::new(id.clone());
        w.draft_title = title.trim().to_string();
        self.shell.open.insert(id.clone(), w);
        self.shell.tabs.push(id.clone());
        let shown = false;
        match start {
            Start::Example(path) => {
                if self.cmd("project_open", serde_json::json!({ "path": path.to_str().unwrap_or_default() })).is_none() {
                    return;
                }
            }
            Start::Results(_) => return,
            Start::Blank { physics, three_d } => {
                if let Err(e) = crate::examples::blank(&mut self.engine, &physics, three_d) {
                    self.set_error(e);
                    return;
                }
            }
        }
        if !shown {
            self.adopt_project();
        }
        self.finish_draft(&id, title, shown);
    }

    /// The results file the window was opening on a job is in: its project is the engine's, so the draft workflow is made around it and the
    /// stored results are put on show.
    pub(crate) fn finish_open(&mut self, summary: &serde_json::Value) {
        let title = self.opening.take().unwrap_or_else(|| "Results".into());
        let id = nyx_engine::facade::store::util::short_id();
        let mut w = OpenWf::new(id.clone());
        w.draft_title = title.trim().to_string();
        self.shell.open.insert(id.clone(), w);
        self.shell.tabs.push(id.clone());
        self.adopt();
        self.project_path = None;
        self.take_solved(summary);
        self.finish_draft(&id, &title, true);
    }

    /// The end of making a draft workflow: remembers its first spec, makes it the loaded one and shows it.
    fn finish_draft(&mut self, id: &str, title: &str, shown: bool) {
        let id = id.to_string();
        let spec = self.spec_json().unwrap_or_default();
        if let Some(w) = self.shell.open.get_mut(&id) {
            w.spec0 = spec.clone();
            w.last_hash = hash_bytes(spec.as_bytes());
        }
        self.shell.loaded = Some((id.clone(), "e1".into()));
        self.shell.page = Page::Workflow(id);
        self.shell.sheet = Sheet::None;
        self.status = if shown { format!("Opened {title}: its project, mesh and results. Nothing is saved until you change something.") } else { "New workflow. Nothing is saved until you change something.".into() };
    }

    pub fn close_tab(&mut self, id: &str) {
        if self.shell.run.as_ref().is_some_and(|r| r.wf_id == id) {
            self.shell.notice = Some("This workflow is still computing. Closing its tab does not stop the run; it keeps going and the dashboard card follows it.".into());
        }
        if self.shell.loaded.as_ref().is_some_and(|l| l.0 == id) {
            self.flush();
        }
        self.shell.tabs.retain(|t| t != id);
        if self.shell.run.as_ref().is_none_or(|r| r.wf_id != id) {
            self.shell.open.remove(id);
        }
        if self.shell.page == Page::Workflow(id.to_string()) {
            self.shell.page = self.shell.tabs.last().map(|t| Page::Workflow(t.clone())).unwrap_or(Page::Dashboard);
        }
        self.shell.closed.push(id.to_string());
        self.save_tabs();
        self.refresh_entries();
    }

    // ------------------------------------------------------------------------------------------------------ runs

    /// Compute with a run record: stores the spec, remembers what is being solved and starts the solve.
    pub fn shell_compute(&mut self) {
        if !self.shell.active() {
            self.compute();
            return;
        }
        let Page::Workflow(wid) = self.shell.page.clone() else {
            self.status = "Open a workflow to compute.".into();
            return;
        };
        if self.busy.is_some() {
            self.compute();
            return;
        }
        if self.shell.loaded.as_ref().map(|l| &l.0) != Some(&wid) {
            self.show_workflow(&wid);
        }
        let Some((_, eid)) = self.shell.loaded.clone() else { return };
        let Ok(spec) = self.spec_json() else { return };
        // a run needs a folder: a draft gets it now (the run is its first change)
        if self.shell.open.get(&wid).is_some_and(|w| w.dir.is_none()) && self.persist_draft(&wid, &spec).is_err() {
            return;
        }
        self.flush();
        let params = params_of(&spec);
        self.shell.run = Some(RunCtx { wf_id: wid.clone(), exp: eid, spec_json: spec, params, started: Instant::now() });
        self.compute();
        if self.busy.is_none() {
            // refused (nothing to solve on, wrong study): there is no job, so no record
            self.shell.run = None;
        }
    }

    /// A solve finished well: writes its run record and keeps the log with it. Returns the run id.
    pub fn record_ok(&mut self, res: &nyx_engine::facade::physics::SolveResult) -> Option<LiveRun> {
        let rc = self.shell.run.take()?;
        let (wf_id, exp) = (rc.wf_id.clone(), rc.exp.clone());
        let rec = RunRecord { sec: res.seconds, status: RunStatus::Ok, params: rc.params.clone(), metrics: nyx_engine::store_link::metrics_of(res), ..Default::default() };
        self.write_run(rc, rec, Some(res.log.join("\n"))).map(|run| LiveRun { wf_id, exp, run, physics: res.physics.clone() })
    }

    /// A convergence check of the result on show has finished: when the engine says the result now rests on a check, the verdict goes into the
    /// run record that solved it (an appended patch line, like a star), so the runs table and a later session know the run was checked.
    pub fn note_check_in_run(&mut self) {
        let Some(live) = self.shell.live_run.clone() else { return };
        let Some(check) = self.engine.validated(&live.physics).check else { return };
        let patch = RunPatch { convergence: Some(check.to_record()), ..Default::default() };
        let Some(w) = self.shell.open.get_mut(&live.wf_id) else { return };
        let Some(d) = w.dir.as_mut() else { return };
        match d.patch_run(&live.exp, &live.run, patch) {
            Ok(()) => {
                if w.exp == live.exp {
                    w.reload_runs();
                }
            }
            Err(e) => self.set_error(format!("Could not record the check in the run: {e}")),
        }
    }

    /// A solve failed or was cancelled.
    pub fn record_end(&mut self, status: RunStatus, error: Option<String>) {
        let Some(rc) = self.shell.run.take() else { return };
        let rec = RunRecord { sec: rc.started.elapsed().as_secs_f64(), status, params: rc.params.clone(), error, ..Default::default() };
        self.write_run(rc, rec, None);
    }

    fn write_run(&mut self, rc: RunCtx, mut rec: RunRecord, log: Option<String>) -> Option<String> {
        let w = self.shell.open.get_mut(&rc.wf_id)?;
        if log.is_some() {
            rec.artifacts.push(nyx_engine::facade::store::runs::Artifact { kind: "log".into(), path: "log.txt".into(), bytes: 0 });
            #[cfg(feature = "hdf5")]
            rec.artifacts.push(nyx_engine::facade::store::runs::Artifact { kind: "result".into(), path: "results.h5".into(), bytes: 0 });
        }
        let res = w.dir.as_mut()?.record_run(&rc.exp, rec, &rc.spec_json);
        match res {
            Ok(r) => {
                let d = w.dir.as_mut()?;
                if let Some(text) = log {
                    let dir = d.run_dir(&rc.exp, &r.id);
                    let _ = std::fs::create_dir_all(&dir);
                    let _ = std::fs::write(dir.join("log.txt"), text);
                    // written on an engine job, from a snapshot of what the engine holds now: the window does not wait for the file
                    #[cfg(feature = "hdf5")]
                    {
                        let target = dir.join("results.h5").to_string_lossy().into_owned();
                        match self.engine.start_export_results(&target, false) {
                            Ok(job) => {
                                self.exports.insert(job, target);
                            }
                            Err(e) => self.error = Some(format!("Results file: {e}")),
                        }
                    }
                }
                if let Some(ws) = &self.shell.ws {
                    // housekeeping after each run, by the person's own settings (the defaults delete nothing)
                    let _ = ws.housekeeping(d, now_secs());
                }
                if w.exp == rc.exp {
                    w.reload_runs();
                }
                w.saved_at = Some(Instant::now());
                self.refresh_entries();
                Some(r.id)
            }
            Err(e) => {
                self.set_error(format!("Could not record the run: {e}"));
                None
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------ frame

    /// The whole window in the shell layout.
    pub fn shell_ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        self.sync_look(&ctx);
        #[cfg(windows)]
        {
            let clicked = self.shell.native_menu.as_ref().map(|m| m.take()).unwrap_or_default();
            for act in clicked {
                self.do_act(&ctx, act);
            }
        }
        self.shell_shortcuts(&ctx);
        if self.shell.logo.is_none() {
            self.shell.logo = load_logo(&ctx);
        }
        if self.shell.last_refresh.elapsed() > Duration::from_secs(5) && self.shell.page == Page::Dashboard {
            self.refresh_entries();
        }
        egui::Panel::top("tabbar").frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::symmetric(14, 8)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| self.top_bar(ui));
        match self.shell.page.clone() {
            Page::Dashboard => self.dashboard_ui(ui),
            Page::Workflow(id) => self.workflow_ui(ui, &id),
            Page::Settings => self.settings_ui(ui),
        }
        self.sheets_ui(&ctx);
        self.cancel_dialog(&ctx);
        if let Some(n) = self.shell.notice.clone() {
            egui::Area::new(egui::Id::unique("notice")).order(egui::Order::Tooltip).anchor(egui::Align2::CENTER_BOTTOM, [0.0, -24.0]).show(&ctx, |ui| {
                crate::kit::card(ui, theme::text(&ctx), theme::text(&ctx), 12, |ui| {
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(n).color(theme::page(&ctx)));
                        if ui.button("OK").clicked() {
                            self.shell.notice = None;
                        }
                    });
                });
            });
        }
        if self.shell.last_tick.elapsed() > Duration::from_millis(1200) {
            self.shell.last_tick = Instant::now();
            if !ctx.input(|i| i.pointer.any_down()) && self.shell.run.is_none() {
                self.flush();
            }
        }
        ctx.request_repaint_after(Duration::from_millis(1200));
    }

    fn shell_shortcuts(&mut self, ctx: &egui::Context) {
        let (k, n, w, o, tab, shift, t) = ctx.input(|i| (i.modifiers.command && i.key_pressed(egui::Key::K), i.modifiers.command && i.key_pressed(egui::Key::N), i.modifiers.command && i.key_pressed(egui::Key::W), i.modifiers.command && i.key_pressed(egui::Key::O), i.modifiers.command && i.key_pressed(egui::Key::Tab), i.modifiers.shift, i.modifiers.command && i.modifiers.shift && i.key_pressed(egui::Key::T)));
        if k {
            self.shell.sheet = if self.shell.sheet == Sheet::Palette { Sheet::None } else { Sheet::Palette };
            self.shell.palette_q.clear();
        }
        if n {
            self.shell.sheet = Sheet::New;
        }
        if o {
            self.pick_nyx_file();
        }
        if t {
            if let Some(id) = self.shell.closed.pop() {
                self.open_tab(&id, true);
            }
        } else if w {
            if let Page::Workflow(id) = self.shell.page.clone() {
                self.close_tab(&id);
            }
        }
        if tab && !self.shell.tabs.is_empty() {
            let pages: Vec<Page> = std::iter::once(Page::Dashboard).chain(self.shell.tabs.iter().map(|t| Page::Workflow(t.clone()))).collect();
            let i = pages.iter().position(|p| *p == self.shell.page).unwrap_or(0);
            let next = if shift { (i + pages.len() - 1) % pages.len() } else { (i + 1) % pages.len() };
            match pages[next].clone() {
                Page::Workflow(id) => self.show_workflow(&id),
                p => {
                    self.flush();
                    self.shell.page = p;
                }
            }
        }
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let running_wf = self.shell.run.as_ref().map(|r| r.wf_id.clone());
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 10.0;
            if let Some(tex) = &self.shell.logo {
                // the whole mark, not clipped: it is a square with its own offset shadow, and a rounded clip cut its corners
                ui.add(egui::Image::new(tex).fit_to_exact_size(Vec2::splat(38.0)));
            } else {
                ui.label(RichText::new("Nyx").font(theme::title(24.0)).color(theme::accent(&ctx)));
            }
            let on_dash = self.shell.page == Page::Dashboard;
            let dash = crate::kit::tab(ui, "dashboard", "Dashboard", on_dash, None, false).select;
            if dash.clicked() {
                self.flush();
                self.shell.page = Page::Dashboard;
                self.refresh_entries();
            }
            let mut close = None;
            let mut show = None;
            for id in self.shell.tabs.clone() {
                let Some(w) = self.shell.open.get(&id) else { continue };
                let on = self.shell.page == Page::Workflow(id.clone());
                let title = crate::fmt::ellipsize(&w.title(), 22);
                let running = running_wf.as_deref() == Some(&id);
                let state = if running { crate::kit::State::Running } else if w.dir.is_none() { crate::kit::State::Idle } else { crate::kit::State::Done };
                let hit = crate::kit::tab(ui, ("tab", &id), &title, on, Some(state), true);
                hit.select.clone().on_hover_text(state.word());
                if hit.select.clicked() {
                    show = Some(id.clone());
                }
                if hit.select.middle_clicked() {
                    close = Some(id.clone());
                }
                if let Some(c) = hit.close {
                    if c.on_hover_text("Close the tab (Ctrl+W). The workflow is already saved.").clicked() {
                        close = Some(id.clone());
                    }
                }
            }
            if ui.add(egui::Button::new(RichText::new("+").size(18.0)).min_size(Vec2::new(34.0, 34.0))).on_hover_text("New workflow (Ctrl+N)").clicked() {
                self.shell.sheet = Sheet::New;
            }
            if let Some(id) = show {
                self.show_workflow(&id);
            }
            if let Some(id) = close {
                self.close_tab(&id);
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.spacing_mut().item_spacing.x = 10.0;
                if ui.add(egui::Button::new("Settings").min_size(Vec2::new(0.0, 32.0))).clicked() {
                    self.flush();
                    self.shell.page = Page::Settings;
                }
                let search = ui.add(egui::Button::new(RichText::new("Search or run anything        Ctrl K").color(theme::muted(&ctx))).min_size(Vec2::new(250.0, 32.0)));
                if search.clicked() {
                    self.shell.sheet = Sheet::Palette;
                    self.shell.palette_q.clear();
                }
                // the busy chip and its queue popover (activity.rs): only while something runs or waits (review I-12: an always-on "Idle" is noise)
                let activity = self.activity();
                self.busy_chip(ui, &activity);
            });
        });
    }

    /// File dialog for a `.nyx` file; shows what is inside before anything is saved.
    pub fn pick_nyx_file(&mut self) {
        let Some(p) = rfd::FileDialog::new().add_filter("nyx workflow", &["nyx"]).add_filter("nyx project (JSON)", &["json"]).add_filter("nyx results (HDF5)", &["h5"]).pick_file() else { return };
        self.open_received(p);
    }

    pub fn open_received(&mut self, p: PathBuf) {
        let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
        if ext == "h5" {
            let title = p.file_stem().and_then(|s| s.to_str()).unwrap_or("Results").to_string();
            self.new_draft(&title, Start::Results(p));
        } else if ext == "nyx" {
            let ins = nyx_engine::facade::store::archive::inspect(&p, &nyx_engine::facade::store::archive::Limits::default());
            self.shell.open_file = Some((p, ins));
            self.shell.sheet = Sheet::Open;
        } else {
            let title = p.file_stem().and_then(|s| s.to_str()).unwrap_or("Imported project").trim_end_matches(".nyx").to_string();
            self.new_draft(&title, Start::Example(p));
        }
    }
}

/// Parameters of a run: the dimensions of the geometry and the wavelength, in SI units, from the concept card of the spec.
pub fn params_of(json: &str) -> BTreeMap<String, f64> {
    let mut m = BTreeMap::new();
    if let Ok(v) = serde_json::from_str::<serde_json::Value>(json) {
        let c = nyx_engine::store_link::card_of(&v);
        for d in &c.geometry.dims {
            let key: String = format!("{}_{}", d.of, d.name).chars().map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' }).collect();
            m.entry(key).or_insert(d.m);
        }
        if let Some(w) = c.wavelength_m {
            m.insert("wavelength".into(), w);
        }
    }
    m
}

fn load_logo(ctx: &egui::Context) -> Option<egui::TextureHandle> {
    let bytes = include_bytes!("../../../docs/brand/nyx-mark.png");
    // the mark is 286 px and drawn at about 38: scale it once with a good filter, or the small drawing shimmers
    let img = image::load_from_memory(bytes).ok()?.to_rgba8();
    let img = image::imageops::resize(&img, 96, 96, image::imageops::FilterType::Lanczos3);
    let size = [img.width() as usize, img.height() as usize];
    Some(ctx.load_texture("nyx-mark", egui::ColorImage::from_rgba_unmultiplied(size, img.as_raw()), egui::TextureOptions::LINEAR))
}

impl App {
    /// "Stop this run?": every way of cancelling a solve (Esc, the Cancel buttons) asks first.
    fn cancel_dialog(&mut self, ctx: &egui::Context) {
        if !self.shell.confirm_cancel {
            return;
        }
        if self.busy.is_none() {
            self.shell.confirm_cancel = false;
            return;
        }
        let (mut stop, mut keep) = (false, false);
        let esc_closed = crate::kit::modal(ctx, "cancel", 460.0, |ui| {
            ui.label(RichText::new("Stop this run?").strong().size(18.0));
            ui.label(crate::kit::muted(ui, "The solver stops at its next safe point; a factorisation in progress has to finish first. Earlier results stay as they are, and the run is recorded as cancelled."));
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if crate::kit::danger(ui, "Stop the run").clicked() {
                    stop = true;
                }
                if crate::kit::primary(ui, "Keep running").clicked() {
                    keep = true;
                }
            });
        });
        if stop {
            self.cancel_job();
        }
        if stop || keep || esc_closed {
            self.shell.confirm_cancel = false;
        }
    }
}

impl App {
    /// Makes the window look the way Settings says: light, dark or the system's (followed live, also when the system changes while nyx runs),
    /// High contrast, reduced motion, and the zoom (a change made with Ctrl +/- or the View menu is saved).
    pub fn sync_look(&mut self, ctx: &egui::Context) {
        let Some(ws) = &self.shell.ws else { return };
        let s = &ws.settings;
        let dark = match s.theme.as_str() {
            "dark" => true,
            "light" => false,
            _ => ctx.system_theme() == Some(egui::Theme::Dark),
        };
        let want = theme::Look { dark, high_contrast: theme::resolve_contrast(&s.contrast), reduce_motion: s.reduce_motion };
        if theme::look(ctx) != want {
            theme::apply_look(ctx, want);
        }
        self.shell.dark = dark;
        let z = ctx.zoom_factor();
        if (z - self.shell.last_zoom).abs() > 1e-4 {
            self.shell.last_zoom = z;
            if let Some(ws) = &mut self.shell.ws {
                ws.settings.zoom = z;
                let _ = ws.save_settings();
            }
        }
    }
}
