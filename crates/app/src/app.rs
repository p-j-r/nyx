//! The application: the state of the window, the background jobs that run in it, and the window layout.
//!
//! `App` holds the project (`state::Project`), the mesh in use (`Doc`), the results shown and the viewport. Anything slow
//! (building a mesh, opening a file, solving) runs as a numbered job on its own thread (`begin_job`, `worker`); only
//! one job runs at a time, and a message coming back for a job that is no longer the running one cannot end the busy state.
//! A solve has a `Progress` handle: the toolbar shows its current stage and a Cancel button (Esc too), and the stage log
//! it writes is part of the result. `with_context` builds the state without a window, which is how the tests drive it.

use crate::colormap::Colormap;
use crate::examples;
use crate::state::*;
use crate::viewport::{Content, Viewport};
use crate::worker::{self, Msg};
use eframe::egui::{self, Color32, RichText};
use femcore::problem::*;
use femcore::progress::Progress;
use std::collections::{HashMap, VecDeque};
use std::sync::atomic::AtomicBool;
use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::Arc;
use std::time::Instant;

/// Work waiting for the running job to finish (the latest request wins).
/// Choices of the File > New dialog.
pub struct NewDialog {
    pub physics: Physics,
    pub three_d: bool,
}

pub enum Queued {
    Open(PathBuf),
    Rebuild,
}

pub struct App {
    pub proj: Project,
    pub doc: Doc,
    pub tab: Tab,
    pub vp: Viewport,
    pub results: Vec<ResultView>,
    pub current: Option<usize>,
    pub show_results: bool,
    pub solve_state: SolveState,
    pub tx: Sender<(u64, Msg)>,
    pub rx: Receiver<(u64, Msg)>,
    pub busy: Option<(String, Instant)>,
    /// number of the last job started
    pub next_job: u64,
    /// the running job: its number, its progress/cancel handle and whether it can be cancelled
    pub job: Option<(u64, Arc<Progress>, bool)>,
    /// scripted runs: cancel the running job at this time (`cancel_after SECONDS`)
    pub cancel_at: Option<Instant>,
    /// a mesh build or file open asked for while another job was running; it starts when that job is done
    pub queued: Option<Queued>,
    pub status: String,
    pub error: Option<String>,
    pub show_edges: bool,
    pub show_numbers: bool,
    pub picking: Option<usize>,
    pub project_path: Option<PathBuf>,
    pub label_buffers: HashMap<(Physics, usize), String>,
    pub geometry_dirty: bool,
    pub geom_edit: GeometrySpec,
    pub geom_ui: crate::geomui::State,
    /// why a boundary rule (physics, index) does not resolve
    pub rule_errors: HashMap<(Physics, usize), String>,
    /// region numbers with their names, for the Materials tab
    pub domain_names: Vec<(i32, String)>,
    /// a rule was edited: resolve them again
    pub rules_dirty: bool,
    /// rules that pick what was just clicked in the viewport
    pub suggestions: Vec<String>,
    pub bottom_open: bool,
    /// the File > New dialog while it is open
    pub new_dialog: Option<NewDialog>,
    /// save every run (results and project) to `autosave_dir`
    pub autosave: bool,
    pub autosave_dir: String,
    pub script_output: Option<String>,
    pub script: VecDeque<String>,
    pub settle: u32,
    pub shot_pending: Arc<AtomicBool>,
    pub shell: crate::shell::Shell,
    pub ctx: egui::Context,
}

/// `Documents\nyx_runs` of the user (the current folder if there is none).
pub fn default_runs_dir() -> String {
    std::env::var("USERPROFILE").or_else(|_| std::env::var("HOME")).map(|h| PathBuf::from(h).join("Documents").join("nyx_runs").display().to_string()).unwrap_or_else(|_| "nyx_runs".into())
}

impl App {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        crate::theme::install_fonts(&cc.egui_ctx);
        let mut app = Self::with_context(cc.egui_ctx.clone());
        let root = std::env::var_os("NYX_WORKSPACE").map(PathBuf::from).unwrap_or_else(nyx_store::workspace::default_root);
        app.attach_workspace(root);
        app
    }

    /// The application state on an egui context; no window is needed (tests use a default context).
    pub fn with_context(ctx: egui::Context) -> Self {
        crate::theme::apply(&ctx);
        let (tx, rx) = channel();
        let proj = Project::default_project();
        let mut app = App {
            geom_edit: proj.geometry.clone(),
            geom_ui: Default::default(),
            rule_errors: HashMap::new(),
            domain_names: vec![],
            rules_dirty: false,
            suggestions: vec![],
            proj,
            doc: Doc::default(),
            tab: Tab::Geometry,
            vp: Viewport::default(),
            results: vec![],
            current: None,
            show_results: false,
            solve_state: SolveState::default(),
            tx,
            rx,
            busy: None,
            next_job: 0,
            job: None,
            cancel_at: None,
            queued: None,
            status: "Ready".into(),
            error: None,
            show_edges: true,
            show_numbers: true,
            picking: None,
            project_path: None,
            label_buffers: HashMap::new(),
            geometry_dirty: false,
            bottom_open: true,
            new_dialog: None,
            autosave: true,
            autosave_dir: default_runs_dir(),
            script_output: None,
            script: VecDeque::new(),
            settle: 6,
            shot_pending: Arc::new(AtomicBool::new(false)),
            shell: crate::shell::Shell::detached(),
            ctx,
        };
        app.rebuild_mesh();
        app
    }

    // ------------------------------------------------------------- jobs

    pub fn set_error(&mut self, e: impl Into<String>) {
        let e = e.into();
        self.status = "Error".into();
        self.error = Some(e);
    }

    /// Marks the window busy and returns the handle for a new background job. Only one job runs at a time; the
    /// caller must have checked `self.busy` first.
    fn begin_job(&mut self, label: String, cancellable: bool) -> worker::Job {
        self.next_job += 1;
        let progress = Progress::new();
        self.busy = Some((label, Instant::now()));
        self.job = Some((self.next_job, progress.clone(), cancellable));
        self.error = None;
        worker::Job { id: self.next_job, tx: self.tx.clone(), ctx: self.ctx.clone(), progress }
    }

    /// Asks the running solve to stop (it does so at its next check; a factorisation in progress has to finish first).
    pub fn cancel_job(&mut self) {
        if let Some((_, p, true)) = &self.job {
            p.cancel();
            self.status = "Cancelling: the solver stops at its next check (a factorisation in progress has to finish first)".into();
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
        self.proj.geometry = self.geom_edit.clone();
        self.proj.source = Source::Generated;
        let job = self.begin_job("Building mesh".into(), false);
        worker::spawn_build(self.proj.geometry.clone(), job);
    }

    pub fn open_path(&mut self, path: PathBuf) {
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
        if ext == "json" || ext == "xml" {
            self.load_project(&path);
            return;
        }
        if self.busy.is_some() {
            self.queued = Some(Queued::Open(path));
            return;
        }
        if ext == "h5" || ext == "hdf5" {
            let job = self.begin_job("Reading HDF5 results".into(), false);
            worker::spawn_open_h5(path, job);
            return;
        }
        self.set_error(format!("Don't know how to open '.{ext}' files (nyx projects .xml and .json, HDF5 results .h5)"));
    }

    pub fn compute(&mut self) {
        if self.refuse_while_busy() {
            return;
        }
        self.sync_rules();
        if let Some(((_, i), e)) = self.rule_errors.iter().find(|((p, _), _)| *p == self.proj.active) {
            let msg = format!("Boundary condition {}: {e}", i + 1);
            self.set_error(msg);
            return;
        }
        let Some(mesh) = self.doc.mesh.clone() else {
            let msg = self.doc.mesh_error.clone().unwrap_or_else(|| "There is no mesh to solve on: build a geometry or open a project first.".into());
            self.set_error(msg);
            return;
        };
        let mut spec = self.proj.active_spec().clone();
        if spec.physics == Physics::Emw && mesh.dim == 2 && !matches!(spec.study, Study::Modes { .. }) {
            self.set_error("On a 2D mesh EM waves runs a mode analysis: open the Study tab to set it up.");
            return;
        }
        if spec.physics == Physics::Emw && mesh.dim == 3 && !matches!(spec.study, Study::Eigenfrequency { .. }) {
            spec.study = Study::Eigenfrequency { n_modes: 6, f_search: 10e9 };
        }
        spec.order = spec.order.clamp(1, 2);
        let job = self.begin_job(format!("Solving {}", spec.physics.name()), true);
        let auto = if self.autosave && !self.autosave_dir.trim().is_empty() {
            self.proj.geometry = self.geom_edit.clone();
            Some(worker::AutoSave { dir: PathBuf::from(self.autosave_dir.trim()), project_xml: self.proj.to_xml() })
        } else {
            None
        };
        worker::spawn_solve(mesh, spec, self.solve_state.clone(), auto, job);
    }

    fn poll(&mut self) {
        while let Ok((job_id, msg)) = self.rx.try_recv() {
            // only the job that is running clears the busy state (a message of an older job must not)
            let started = if self.job.as_ref().map(|j| j.0) == Some(job_id) {
                self.job = None;
                self.busy.take().map(|b| b.1)
            } else {
                None
            };
            let took = started.map(|s| s.elapsed().as_secs_f64()).unwrap_or(0.0);
            match msg {
                Msg::Built(Ok(b)) => {
                    self.doc.mesh = Some(b.mesh.clone());
                    self.doc.view = Some(b.view.clone());
                    self.doc.source = "generated geometry".into();
                    self.doc.mesh_error = None;
                    self.after_mesh_change();
                    self.status = format!("Mesh built: {} cells, {} boundary facets ({took:.2} s)", b.mesh.n_cells(), b.mesh.n_facets());
                    self.geometry_dirty = false;
                }
                Msg::Built(Err(e)) => self.set_error(format!("Geometry: {e}")),
                Msg::Solved { physics, res: Ok(s) } => {
                    if physics == Physics::Heat {
                        self.solve_state = s.state.clone();
                    }
                    let dim = self.doc.mesh.as_ref().map(|m| m.dim).unwrap_or(2);
                    let name = physics.name().to_string();
                    let secs = s.result.seconds;
                    self.shell.live_run = self.record_ok(&s.result);
                    let s_output = s.output;
                    let rv = ResultView::new(name.clone(), Some(physics), s.result, s.display, dim);
                    self.push_result(rv);
                    self.status = format!("{name}: solved in {secs:.2} s");
                    match s_output {
                        Some(Ok(m)) => self.status += &format!("; {m}"),
                        Some(Err(e)) => self.error = Some(format!("Results file: {e}")),
                        None => {}
                    }
                    self.tab = Tab::Results;
                }
                Msg::Solved { res: Err(e), .. } if e == femcore::progress::CANCELLED => {
                    self.error = None;
                    self.record_end(nyx_store::runs::RunStatus::Cancelled, None);
                    self.status = format!("Solve cancelled after {took:.1} s; the previous results are unchanged");
                }
                Msg::Solved { res: Err(e), .. } => {
                    self.record_end(nyx_store::runs::RunStatus::Failed, Some(e.clone()));
                    self.set_error(e)
                }
                Msg::H5(name, Ok((mesh, s))) => {
                    self.doc.view = Some(Arc::new(femcore::view::ViewMesh::from_simplex(&mesh)));
                    self.doc.mesh = Some(mesh.clone());
                    self.doc.source = name.clone();
                    self.after_mesh_change();
                    let rv = ResultView::new(format!("HDF5 {}", std::path::Path::new(&name).file_name().and_then(|n| n.to_str()).unwrap_or("file")), None, s.result, s.display, mesh.dim);
                    self.push_result(rv);
                    self.status = format!("Loaded results from {name}");
                    self.tab = Tab::Results;
                }
                Msg::H5(_, Err(e)) => self.set_error(e),
                Msg::Script(Ok(t)) => self.script_output = Some(t),
                Msg::Script(Err(e)) => self.set_error(e),
            }
        }
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
        self.solve_state = SolveState::default();
        self.picking = None;
        self.vp.fit_token += 1;
        self.sync_rules();
    }

    /// Turns the selector rules of the boundary conditions into the integers of the built mesh, and reads the names of the regions.
    pub fn sync_rules(&mut self) {
        self.rules_dirty = false;
        self.rule_errors.clear();
        let Some(mesh) = self.doc.mesh.clone() else { return };
        let Ok(names) = self.proj.geometry.names(&mesh) else { return };
        self.domain_names = names.domains();
        for spec in &mut self.proj.specs {
            for (i, bc) in spec.bcs.iter_mut().enumerate() {
                if let Some(rule) = bc.on.clone() {
                    match names.select(&rule) {
                        Ok(l) => bc.labels = l,
                        Err(e) => {
                            self.rule_errors.insert((spec.physics, i), e);
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------- files

    pub fn save_project(&mut self, as_new: bool) {
        let path = if as_new || self.project_path.is_none() {
            rfd::FileDialog::new().add_filter("nyx project (XML)", &["xml"]).add_filter("nyx project (JSON)", &["json"]).set_file_name("project.nyx.xml").save_file()
        } else {
            self.project_path.clone()
        };
        let Some(path) = path else { return };
        self.proj.geometry = self.geom_edit.clone();
        let is_xml = path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("xml"));
        let text = if is_xml { self.proj.to_xml() } else { serde_json::to_string_pretty(&self.proj).map_err(|e| e.to_string()) };
        match text.and_then(|s| std::fs::write(&path, s).map_err(|e| e.to_string())) {
            Ok(()) => {
                self.status = format!("Saved {}", path.display());
                self.project_path = Some(path);
            }
            Err(e) => self.set_error(format!("Could not save the project: {e}")),
        }
    }

    pub fn load_project(&mut self, path: &std::path::Path) {
        let is_xml = path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("xml"));
        let parse = |s: String| if is_xml { Project::from_xml(&s) } else { serde_json::from_str::<Project>(&s).map_err(|e| format!("not a nyx project: {e}")) };
        match std::fs::read_to_string(path).map_err(|e| e.to_string()).and_then(parse) {
            Ok(p) => {
                let source = p.source.clone();
                self.geom_edit = p.geometry.clone();
                self.proj = p;
                let frame = self.proj.specs.first().map(|s| s.materials.frame).unwrap_or_default();
                for ph in Physics::ALL {
                    if !self.proj.specs.iter().any(|s| s.physics == ph) {
                        let mut ns = ProblemSpec::new(ph);
                        ns.materials.frame = frame;
                        self.proj.specs.push(ns);
                    }
                }
                self.project_path = Some(path.to_path_buf());
                self.label_buffers.clear();
                match source {
                    Source::Generated => self.rebuild_mesh(),
                }
                self.status = format!("Loaded {}", path.display());
            }
            Err(e) => self.set_error(e),
        }
    }

    pub fn load_example(&mut self, i: usize) {
        let p = examples::by_index(i);
        self.geom_edit = p.geometry.clone();
        self.proj = p;
        self.label_buffers.clear();
        self.project_path = None;
        self.rebuild_mesh();
        self.tab = Tab::Geometry;
        self.status = format!("Loaded example: {}", examples::NAMES[i]);
    }

    /// CSV of the current solution: coordinates of every dof and all fields.
    pub fn result_csv(&self) -> Option<String> {
        let (i, mesh) = (self.current?, self.doc.mesh.clone()?);
        let rv = &self.results[i];
        let sol = rv.solution();
        let order = rv.result.order;
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
        let Some(a) = self.current.and_then(|i| self.results[i].result.analysis.as_ref()) else {
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
        self.proj.geometry = self.geom_edit.clone();
        let json = serde_json::to_string_pretty(&self.proj).unwrap_or_default();
        if self.refuse_while_busy() {
            return;
        }
        let job = self.begin_job(format!("Running {}", script.file_name().and_then(|n| n.to_str()).unwrap_or("script")), false);
        worker::spawn_script(script, json, self.result_csv(), job);
    }

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new().add_filter("nyx projects and HDF5 results", &["json", "xml", "h5", "hdf5"]).pick_file() {
            self.open_path(p);
        }
    }

    // ------------------------------------------------------------- layout

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("New... (Ctrl+N)").clicked() {
                    ui.close();
                    self.new_dialog = Some(NewDialog { physics: self.proj.active, three_d: false });
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
                for (i, n) in examples::NAMES.iter().enumerate() {
                    if ui.button(*n).clicked() {
                        ui.close();
                        self.load_example(i);
                    }
                }
            });
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

    /// Replaces the project by a new one (nothing is kept of the old geometry, conditions or results).
    pub fn new_project(&mut self, p: Project, what: &str) {
        self.geom_edit = p.geometry.clone();
        self.proj = p;
        self.label_buffers.clear();
        self.project_path = None;
        self.rebuild_mesh();
        self.tab = Tab::Geometry;
        self.status = format!("New simulation: {what}. Set the geometry, materials and conditions, then Compute (F5).");
    }

    /// The File > New window: physics, 2D or 3D and how the geometry is made; Create starts an empty project.
    fn new_dialog_ui(&mut self, ctx: &egui::Context) {
        let Some(mut d) = self.new_dialog.take() else { return };
        let (mut open, mut create, mut cancel) = (true, false, false);
        egui::Window::new("New simulation").open(&mut open).collapsible(false).resizable(false).show(ctx, |ui| {
            egui::Grid::new("new_sim").num_columns(2).spacing([12.0, 8.0]).show(ui, |ui| {
                ui.label("Physics");
                egui::ComboBox::from_id_salt("new_physics").selected_text(d.physics.name()).show_ui(ui, |ui| {
                    for p in Physics::ALL {
                        ui.selectable_value(&mut d.physics, p, p.name());
                    }
                });
                ui.end_row();
                if d.physics == Physics::Brillouin {
                    d.three_d = false;
                }
                ui.label("Space");
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut d.three_d, false, "2D (cross-section)");
                    ui.add_enabled_ui(d.physics != Physics::Brillouin, |ui| ui.selectable_value(&mut d.three_d, true, "3D (the section swept along y)"));
                });
                ui.end_row();
            });
            ui.add_space(6.0);
            ui.weak("The new project has default materials and no boundary conditions; everything can be changed in the tabs on the left.");
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
            let what = format!("{}, {}", d.physics.name(), if d.three_d { "3D" } else { "2D" });
            self.new_project(examples::blank(d.physics, d.three_d), &what);
        } else if open && !cancel {
            self.new_dialog = Some(d);
        }
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Physics");
            let mut active = self.proj.active;
            egui::ComboBox::from_id_salt("active_physics").selected_text(active.name()).show_ui(ui, |ui| {
                for p in Physics::ALL {
                    ui.selectable_value(&mut active, p, p.name());
                }
            });
            if active != self.proj.active {
                self.proj.active = active;
                self.picking = None;
            }
            let busy = self.busy.is_some();
            let has_mesh = self.doc.mesh.is_some();
            let b = ui.add_enabled(!busy && has_mesh, egui::Button::new(RichText::new("  Compute  ").strong()).fill(Color32::from_rgb(38, 110, 70)));
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
                    if ui.add_enabled(!cancelled, egui::Button::new("Cancel")).on_hover_text("Stop the solve at the next safe point (Esc). A factorisation in progress has to finish first. Earlier results stay as they are.").clicked() {
                        self.cancel_job();
                    }
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
                if let Some(m) = &self.doc.mesh {
                    ui.weak(format!("{}  |  {}D  |  {} points  {} cells  |  P{} dofs: {}", self.doc.source, m.dim, m.n_points(), m.n_cells(), self.proj.active_spec().order, m.n_dofs(self.proj.active_spec().order) * if self.proj.active == Physics::Solid { m.dim } else { 1 }));
                } else if let Some(v) = &self.doc.view {
                    ui.weak(format!("{}  |  {}D  |  view only, {} cells", self.doc.source, v.dim, v.n_cells));
                }
            });
        });
    }

    pub fn viewport_ui(&mut self, ui: &mut egui::Ui) {
        let bg = Color32::from_rgb(24, 26, 30);
        let mut bc_colors: HashMap<i32, Color32> = HashMap::new();
        let phys = self.proj.active;
        for (i, bc) in self.proj.spec(phys).bcs.iter().enumerate() {
            for &l in &bc.labels {
                bc_colors.insert(l, crate::colormap::palette(i + 3));
            }
        }
        let selected: Vec<i32> = match self.picking {
            Some(i) => self.proj.spec(phys).bcs.get(i).map(|b| b.labels.clone()).unwrap_or_default(),
            None => vec![],
        };
        let picking = self.picking.is_some();
        let cmap = self.proj.colormap;
        let (show_edges, show_numbers, show_results) = (self.show_edges, self.show_numbers, self.show_results);
        // while the geometry is being edited the exact model is drawn, not the mesh of an earlier version
        let preview = if self.tab == Tab::Geometry && (self.geometry_dirty || self.geom_edit != self.proj.geometry) { self.geom_ui.preview.clone() } else { None };
        let show_results = show_results && preview.is_none();
        let view = preview.or_else(|| self.doc.view.clone());
        
        let out = match (self.current, view) {
            (Some(i), Some(view)) if show_results => {
                let rv = &mut self.results[i];
                let values = rv.field_values();
                let range = rv.range();
                let warp = rv.displacement().filter(|_| rv.warp != 0.0).map(|d| (d, rv.warp));
                let sol = rv.sol.min(rv.result.solutions.len() - 1);
                let nf = rv.result.solutions[sol].fields.len();
                let f = &rv.result.solutions[sol].fields[rv.field.min(nf - 1)];
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
            (_, Some(view)) => {
                self.vp.show(ui, Content::Mesh { view: &view, bc_colors: &bc_colors, selected: &selected, show_edges, show_numbers }, picking, bg)
            }
            _ => {
                self.vp.show(ui, Content::Empty, false, bg)
            }
        };
        if let (Some(label), Some(i)) = (out.picked, self.picking) {
            let bcs = &mut self.proj.spec_mut(phys).bcs;
            let mut now = vec![];
            if let Some(bc) = bcs.get_mut(i) {
                if let Some(p) = bc.labels.iter().position(|&l| l == label) {
                    bc.labels.remove(p);
                } else {
                    bc.labels.push(label);
                    bc.labels.sort();
                }
                // a hand-picked set is not a rule any more; the suggestions below write it as one
                bc.on = None;
                now = bc.labels.clone();
            }
            self.suggestions = match self.doc.mesh.clone() {
                Some(mesh) => self.proj.geometry.names(&mesh).map(|n| n.suggest(&now)).unwrap_or_default(),
                None => vec![],
            };
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
        if self.shell.active() {
            self.shell_ui(ui);
            return;
        }
        if ui.input(|i| i.modifiers.command && i.key_pressed(egui::Key::N)) {
            self.new_dialog = Some(NewDialog { physics: self.proj.active, three_d: false });
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
        if self.bottom_open && self.show_results && self.current.is_some_and(|i| self.results[i].result.solutions.len() > 1) {
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

#[allow(dead_code)]
fn _unused(_: Arc<()>, _: Colormap) {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn app() -> App {
        let mut a = App::with_context(egui::Context::default());
        a.autosave = false; // tests must not write into the user's Documents
        a
    }

    /// Processes finished jobs until nothing runs (or panics after `secs`).
    fn wait_idle(a: &mut App, secs: u64) {
        let t0 = Instant::now();
        loop {
            a.poll();
            a.run_queued();
            if a.busy.is_none() && a.queued.is_none() {
                return;
            }
            assert!(t0.elapsed() < Duration::from_secs(secs), "still busy after {secs} s: {:?}", a.busy);
            std::thread::sleep(Duration::from_millis(10));
        }
    }

    #[test]
    fn a_new_simulation_starts_empty_and_can_be_set_up_and_solved() {
        let mut a = app();
        a.load_example(0);
        wait_idle(&mut a, 30);
        for physics in Physics::ALL {
            for three_d in [false, true] {
                {
                    let p = crate::examples::blank(physics, three_d);
                    assert_eq!(p.active, physics);
                    assert!(p.specs.iter().all(|s| s.bcs.is_empty()), "no boundary conditions in a new project");
                    let mesh = p.geometry.build().unwrap_or_else(|e| panic!("{physics:?} 3d={three_d}: {e}"));
                    assert_eq!(mesh.dim, if three_d && physics != Physics::Brillouin { 3 } else { 2 });
                }
            }
        }
        a.new_project(crate::examples::blank(Physics::Electrostatics, false), "test");
        wait_idle(&mut a, 30);
        assert!(a.results.is_empty() && a.doc.mesh.is_some() && a.project_path.is_none());
        // set up from scratch: plates on the left and right, then compute
        let names = a.proj.geometry.build_named().unwrap();
        let (left, right) = (names.select("outer & plane(x = 0)").unwrap(), names.select("outer & plane(x = 1)").unwrap());
        let es = a.proj.spec_mut(Physics::Electrostatics);
        es.bcs = vec![Bc { on: None, labels: left, kind: BcKind::Potential(1.0) }, Bc { on: None, labels: right, kind: BcKind::Ground }];
        a.compute();
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert_eq!(a.results.len(), 1);
    }

    #[test]
    fn an_example_is_built_computed_and_shown() {
        let mut a = app();
        a.load_example(0);
        assert!(a.busy.is_some());
        wait_idle(&mut a, 30);
        assert!(a.doc.mesh.is_some(), "{:?}", a.error);
        a.compute();
        assert!(a.busy.is_some(), "compute started: {}", a.status);
        wait_idle(&mut a, 60);
        assert!(a.error.is_none(), "{:?}", a.error);
        assert_eq!(a.results.len(), 1);
        assert!(a.status.contains("solved"), "{}", a.status);
        // the stage log of the run is part of the result
        assert!(a.results[0].result.log.iter().any(|l| l.contains("solving")), "{:?}", a.results[0].result.log);
    }

    #[test]
    fn a_second_job_is_refused_while_one_runs_and_a_mesh_build_waits_its_turn() {
        let mut a = app();
        a.load_example(4);
        assert!(a.busy.is_some());
        a.compute(); // nothing to solve on yet and a job is running: refused with a message, not started
        assert!(a.status.starts_with("Busy"), "{}", a.status);
        let first = a.job.as_ref().map(|j| j.0);
        a.rebuild_mesh(); // queued
        assert!(a.queued.is_some());
        assert_eq!(a.job.as_ref().map(|j| j.0), first, "the running job is not replaced");
        wait_idle(&mut a, 30);
        assert!(a.doc.mesh.is_some() && a.queued.is_none());
    }

    #[test]
    fn a_message_of_an_older_job_does_not_clear_the_busy_state() {
        let mut a = app();
        a.load_example(0);
        let running = a.job.as_ref().map(|j| j.0).unwrap();
        a.tx.send((running + 100, Msg::Script(Ok("from an old job".into())))).unwrap();
        a.poll();
        assert_eq!(a.script_output.as_deref(), Some("from an old job"));
        assert!(a.busy.is_some() || a.job.is_none(), "only the running job's own message ends it");
        wait_idle(&mut a, 30);
        assert!(a.busy.is_none() && a.job.is_none());
    }

    #[test]
    fn a_long_solve_can_be_cancelled_and_earlier_results_stay() {
        let mut a = app();
        a.load_example(0);
        wait_idle(&mut a, 30);
        a.compute();
        wait_idle(&mut a, 60);
        assert_eq!(a.results.len(), 1);
        // a big 3D eigenfrequency problem on the same window
        a.load_example(4);
        wait_idle(&mut a, 30);
        a.geom_edit = GeometrySpec::plate(1.2, 1.2, 0.1, Some(6.0));
        a.rebuild_mesh();
        wait_idle(&mut a, 60);
        let kept = a.results.len();
        let phys = a.proj.active;
        let spec = a.proj.spec_mut(phys);
        spec.study = Study::Eigenfrequency { n_modes: 6, f_search: 0.0 };
        spec.solver = LinearSolver::Iterative;
        spec.order = 2;
        a.compute();
        assert!(a.busy.is_some(), "{} / {:?}", a.status, a.error);
        std::thread::sleep(Duration::from_millis(1500));
        let t = Instant::now();
        a.cancel_job();
        wait_idle(&mut a, 60);
        assert!(t.elapsed() < Duration::from_secs(40), "cancelled after {:?}", t.elapsed());
        assert!(a.status.contains("cancelled"), "{}", a.status);
        assert!(a.error.is_none());
        assert_eq!(a.results.len(), kept, "a cancelled solve adds no result");
    }

    #[test]
    fn a_project_saves_and_loads_again() {
        let mut a = app();
        a.load_example(5);
        wait_idle(&mut a, 30);
        let path = std::env::temp_dir().join(format!("nyx-test-{}.nyx.xml", std::process::id()));
        a.project_path = Some(path.clone());
        a.save_project(false);
        assert!(path.exists(), "{}", a.status);
        let mut b = app();
        b.load_project(&path);
        wait_idle(&mut b, 30);
        assert_eq!(a.proj.to_xml().unwrap(), b.proj.to_xml().unwrap());
        assert!(b.doc.mesh.is_some());
        let _ = std::fs::remove_file(&path);
    }
}