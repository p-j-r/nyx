//! The workflow view: header, experiments list, the three sections (1 Concept, 2 Build, 3 Analysis) and what hangs off them (the live run,
//! the runs table and plot, compare, the field viewer and the model screen).
//!
//! The editors are the ones of `tabs.rs`; this file only decides where they sit. While a solve of this workflow runs, editing is locked.
//! Anything the app cannot do yet (sweeps, surrogate models, stored fields) is shown as such, never with made-up numbers.

use crate::app::App;
use crate::dashboard::physics_label;
use crate::kit;
use crate::shell::{OpenWf, Page, Part, Sheet, Sub};
use crate::state::Tab;
use crate::theme;
use eframe::egui::{self, Color32, RichText, Sense, Vec2};
use nyx_engine::facade::store::runs::{Filter, RunPatch, RunRecord, RunStatus};
use nyx_engine::facade::store::util::{ago, now_secs, parse_iso};
use std::time::Instant;

pub enum Act {
    Part(Part),
    Sub(Sub),
    LoadExp(String),
    NewExp,
    Fork,
    Link(String, String, String),
    Compute,
    Cancel,
    Star(String, bool),
    Title(String),
    Tags(Vec<String>),
    ExportCsv,
    Notes(String),
    Pin,
    Done,
    OpenFolder,
    Close,
    Script,
    Field(String),
}

/// A length shown in nm, other numbers with four significant digits.
pub fn fmt_value(key: &str, v: f64) -> String {
    let lengthy = v != 0.0 && v.abs() < 1e-3 && (key.contains("width") || key.contains("height") || key.contains("thick") || key.contains("diam") || key.contains("wavelength") || key.contains("layer") || key.contains("box") || key.starts_with("l") || key.contains("gap"));
    if lengthy {
        format!("{} nm", trim(v * 1e9))
    } else if v != 0.0 && (v.abs() < 1e-3 || v.abs() >= 1e6) {
        format!("{v:.3e}")
    } else {
        format!("{v:.4}")
    }
}

fn trim(x: f64) -> String {
    let s = format!("{x:.3}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

impl App {
    fn locked(&self, id: &str) -> bool {
        self.shell.run.as_ref().is_some_and(|r| r.wf_id == id) && self.busy.is_some()
    }

    pub fn workflow_ui(&mut self, ui: &mut egui::Ui, id: &str) {
        let Some(mut w) = self.shell.open.remove(id) else {
            self.shell.page = Page::Dashboard;
            return;
        };
        let mut acts: Vec<Act> = vec![];
        let locked = self.locked(id);
        let ctx = ui.ctx().clone();
        self.wf_header(ui, &mut w, &mut acts);
        self.wf_sidebar(ui, &mut w, &mut acts, locked);
        self.wf_sectionbar(ui, &mut w, &mut acts, locked);
        match (w.sub, w.part) {
            (Sub::Compare, _) => self.wf_compare(ui, &mut w, &mut acts),
            (Sub::Model, _) => self.wf_model(ui, &mut w, &mut acts),
            (Sub::Field, _) => self.wf_field(ui, &mut w, &mut acts),
            (Sub::Main, Part::Concept) => self.wf_concept(ui, &mut w, locked),
            (Sub::Main, Part::Build) => self.wf_build(ui, &mut w, &mut acts, locked),
            (Sub::Main, Part::Analysis) => self.wf_analysis(ui, &mut w, &mut acts),
        }
        if w.notes_open {
            let mut open = true;
            egui::Window::new("Notes").open(&mut open).default_size([460.0, 320.0]).show(&ctx, |ui| {
                let r = ui.add(egui::TextEdit::multiline(&mut w.notes).desired_width(f32::INFINITY).desired_rows(12).hint_text("Free text. Saved with the workflow."));
                if r.lost_focus() {
                    acts.push(Act::Notes(w.notes.clone()));
                }
            });
            if !open {
                w.notes_open = false;
                acts.push(Act::Notes(w.notes.clone()));
            }
        }
        self.shell.open.insert(id.to_string(), w);
        for a in acts {
            self.apply_act(id, a);
        }
    }

    /// Makes sure the workflow of the tab has a folder (a draft gets one on its first change). `false` when that failed.
    fn ensure_dir(&mut self, id: &str) -> bool {
        if self.shell.open.get(id).is_some_and(|w| w.dir.is_some()) {
            return true;
        }
        let Ok(spec) = self.spec_json() else { return false };
        self.persist_draft(id, &spec).is_ok()
    }

    fn apply_act(&mut self, id: &str, a: Act) {
        match a {
            Act::Part(p) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.part = p;
                    w.sub = Sub::Main;
                }
                self.picking = None;
            }
            Act::Sub(s) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.sub = s;
                }
            }
            Act::LoadExp(e) => self.load_experiment(id, &e),
            Act::NewExp => {
                if !self.ensure_dir(id) {
                    return;
                }
                self.flush();
                let spec = crate::examples::blank_json(&self.active().unwrap_or_else(|| "ems".into()), false);
                let n = self.shell.open.get(id).and_then(|w| w.dir.as_ref()).map(|d| d.wf.experiments.len() + 1).unwrap_or(2);
                if let Some(Ok(spec)) = Some(spec) {
                    let res = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| d.add_experiment(&format!("Experiment {n}"), &spec, None));
                    match res {
                        Some(Ok(e)) => {
                            self.refresh_entries();
                            self.load_experiment(id, &e);
                        }
                        Some(Err(e)) => self.set_error(e),
                        None => {}
                    }
                }
            }
            Act::Fork => {
                if !self.ensure_dir(id) {
                    return;
                }
                self.flush();
                let Some(w) = self.shell.open.get_mut(id) else { return };
                let (eid, title) = (w.exp.clone(), w.dir.as_ref().and_then(|d| d.wf.experiment(&w.exp)).map(|e| e.title.clone()).unwrap_or_default());
                let res = w.dir.as_mut().map(|d| d.fork(&eid, &format!("{title} (fork)")));
                match res {
                    Some(Ok(e)) => {
                        self.refresh_entries();
                        self.load_experiment(id, &e);
                    }
                    Some(Err(e)) => self.set_error(e),
                    None => {}
                }
            }
            Act::Link(from, to, label) => {
                let r = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| d.link(&from, &to, &label));
                if let Some(Err(e)) = r {
                    self.set_error(e);
                }
                self.refresh_entries();
            }
            Act::Compute => self.shell_compute(),
            Act::Cancel => self.shell.confirm_cancel = true,
            Act::Star(rid, star) => {
                let r = self.shell.open.get_mut(id).and_then(|w| {
                    let eid = w.exp.clone();
                    w.dir.as_mut().map(|d| d.patch_run(&eid, &rid, RunPatch { star: Some(star), ..Default::default() }))
                });
                if let Some(Err(e)) = r {
                    self.set_error(e);
                }
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.reload_runs();
                }
            }
            Act::Title(t) => {
                if t.trim().is_empty() {
                    return;
                }
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.draft_title = t.trim().to_string();
                }
                if self.ensure_dir(id) {
                    let r = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| d.set_title(&t));
                    if let Some(Err(e)) = r {
                        self.set_error(e);
                    }
                    self.refresh_entries();
                }
            }
            Act::Tags(t) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.tags = t.clone();
                }
                if self.ensure_dir(id) {
                    let r = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| d.set_tags(t));
                    if let Some(Err(e)) = r {
                        self.set_error(e);
                    }
                    self.refresh_entries();
                }
            }
            Act::Notes(text) => {
                if text.trim().is_empty() && self.shell.open.get(id).is_some_and(|w| w.dir.is_none()) {
                    return;
                }
                if self.ensure_dir(id) {
                    let r = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| d.set_notes(&text));
                    if let Some(Err(e)) = r {
                        self.set_error(e);
                    }
                }
            }
            Act::Pin | Act::Done => {
                if !self.ensure_dir(id) {
                    return;
                }
                let done = matches!(a, Act::Done);
                let r = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()).map(|d| if done { d.set_status(if d.wf.status == nyx_engine::facade::store::Status::Done { nyx_engine::facade::store::Status::Active } else { nyx_engine::facade::store::Status::Done }) } else { d.set_pinned(!d.wf.pinned) });
                if let Some(Err(e)) = r {
                    self.set_error(e);
                }
                self.refresh_entries();
            }
            Act::OpenFolder => {
                if let Some(d) = self.shell.open.get(id).and_then(|w| w.dir.as_ref()).map(|d| d.root.clone()) {
                    self.reveal(&d);
                } else {
                    self.status = "Nothing is saved yet: the folder appears with the first change.".into();
                }
            }
            Act::Close => self.close_tab(id),
            Act::ExportCsv => self.export_runs_csv(id),
            Act::Script => self.run_python_script(),
            Act::Field(rid) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.focus = Some(rid);
                    w.sub = Sub::Field;
                }
            }
        }
    }

    fn export_runs_csv(&mut self, id: &str) {
        let Some(w) = self.shell.open.get(id) else { return };
        let f = Filter::parse(&w.filter_text).unwrap_or_default();
        let rows: Vec<&RunRecord> = w.runs.iter().filter(|r| f.matches(r) && (w.show_failed || r.status == RunStatus::Ok)).collect();
        let csv = nyx_engine::facade::store::runs::to_csv(&rows);
        let name = format!("{}-runs.csv", nyx_engine::facade::store::util::slug(&w.title()));
        let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name(name).save_file() else { return };
        match std::fs::write(&path, csv) {
            Ok(()) => self.status = format!("Exported {} runs to {}", rows.len(), path.display()),
            Err(e) => self.set_error(format!("Could not write the file: {e}")),
        }
    }

    // ------------------------------------------------------------------------------------------------------ header

    fn wf_header(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        egui::Panel::top("wf_header").frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::symmetric(18, 10)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.horizontal(|ui| {
                let eid = egui::Id::unique(("title_edit", &w.id));
                let focused = ctx.memory(|m| m.has_focus(eid));
                let mut buf = if focused { w.draft_buf.clone() } else { w.title() };
                let width = buf.chars().count().max(10) as f32 * 14.0;
                let r = ui.add(egui::TextEdit::singleline(&mut buf).id(eid).font(theme::title(30.0)).frame(egui::Frame::NONE).desired_width(width.min(460.0)).text_color(theme::accent(&ctx)));
                if focused || r.has_focus() {
                    w.draft_buf = buf.clone();
                }
                if r.lost_focus() && buf.trim() != w.title() {
                    acts.push(Act::Title(buf.clone()));
                }
                // tags
                let mut tags = w.tags.clone();
                let mut changed = false;
                tags.retain(|t| {
                    let mut keep = true;
                    ui.spacing_mut().item_spacing.x = 2.0;
                    kit::chip(ui, t.as_str()).on_hover_text("Click to remove").clicked().then(|| keep = false);
                    changed |= !keep;
                    keep
                });
                ui.spacing_mut().item_spacing.x = 8.0;
                let t = ui.add(egui::TextEdit::singleline(&mut w.tag_input).hint_text("+ tag").desired_width(70.0).frame(egui::Frame::NONE));
                if t.lost_focus() && ctx.input(|i| i.key_pressed(egui::Key::Enter)) && !w.tag_input.trim().is_empty() {
                    tags.push(w.tag_input.trim().to_lowercase());
                    w.tag_input.clear();
                    changed = true;
                }
                if changed {
                    acts.push(Act::Tags(tags));
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if kit::ghost(ui, "Share").clicked() {
                        self.shell.sheet = Sheet::Share;
                    }
                    let v = w.dir.as_ref().map(|d| d.wf.version).unwrap_or(0);
                    if ui.add_enabled(w.dir.is_some(), egui::Button::new(format!("Versions  v{v}")).min_size(Vec2::new(0.0, 30.0))).on_disabled_hover_text("Versions start with the first change").clicked() {
                        self.shell.sheet = Sheet::Versions;
                        self.shell.version_sel = None;
                    }
                    ui.menu_button("More", |ui| {
                        let (pinned, done) = w.dir.as_ref().map(|d| (d.wf.pinned, d.wf.status == nyx_engine::facade::store::Status::Done)).unwrap_or((false, false));
                        if ui.button(if pinned { "Unpin from the dashboard" } else { "Pin to the dashboard" }).clicked() {
                            acts.push(Act::Pin);
                            ui.close();
                        }
                        if ui.button(if done { "Mark active again" } else { "Mark done (moves to Earlier)" }).clicked() {
                            acts.push(Act::Done);
                            ui.close();
                        }
                        if ui.button("Notes").clicked() {
                            w.notes_open = true;
                            ui.close();
                        }
                        if ui.button("Open folder").clicked() {
                            acts.push(Act::OpenFolder);
                            ui.close();
                        }
                        if ui.button("Close tab (Ctrl+W)").clicked() {
                            acts.push(Act::Close);
                            ui.close();
                        }
                    });
                    let text = match (&w.dir, w.saved_at) {
                        (None, _) => "Draft: nothing is saved until you change something".to_string(),
                        (Some(_), Some(t)) => format!("Saved automatically · {}", ago(now_secs() - t.elapsed().as_secs() as i64, now_secs())),
                        (Some(d), None) => format!("Saved · {}", ago(d.wf.touched_secs(), now_secs())),
                    };
                    ui.label(kit::muted(ui, text));
                    kit::dot(ui, if w.dir.is_some() { theme::GREEN } else { theme::line(&ctx) });
                });
            });
        });
    }

    // ------------------------------------------------------------------------------------------------------ experiments

    fn wf_sidebar(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        let ctx = ui.ctx().clone();
        egui::Panel::left("wf_experiments").exact_size(250.0).resizable(false).frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.horizontal(|ui| {
                kit::label_caps(ui, "Experiments", "");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.menu_button("+", |ui| {
                        if ui.add_enabled(!locked, egui::Button::new("New blank experiment")).clicked() {
                            acts.push(Act::NewExp);
                            ui.close();
                        }
                        if ui.add_enabled(!locked, egui::Button::new("Fork this experiment")).on_hover_text("Copies the set-up and remembers its parent").clicked() {
                            acts.push(Act::Fork);
                            ui.close();
                        }
                        if let Some(d) = &w.dir {
                            let others: Vec<(String, String)> = d.wf.experiments.iter().filter(|e| e.id != w.exp).map(|e| (e.id.clone(), e.title.clone())).collect();
                            if !others.is_empty() {
                                ui.menu_button("Link this to...", |ui| {
                                    for (eid, title) in others {
                                        for label in ["same geometry", "validates", "compares"] {
                                            if ui.button(format!("{label}: {title}")).clicked() {
                                                acts.push(Act::Link(w.exp.clone(), eid.clone(), label.into()));
                                                ui.close();
                                            }
                                        }
                                    }
                                });
                            }
                        }
                    });
                });
            });
            ui.add_space(6.0);
            egui::ScrollArea::vertical().auto_shrink([false, false]).max_height(ui.available_height() - 150.0).show(ui, |ui| {
                #[allow(clippy::type_complexity)]
                let exps: Vec<(String, String, Option<String>, usize, String, Vec<(String, String)>)> = match &w.dir {
                    Some(d) => d.wf.experiments.iter().map(|e| (e.id.clone(), e.title.clone(), e.parent.clone(), e.runs.count, format!("{} · {}D", physics_label(&e.card.physics), e.card.dim), e.links.iter().map(|l| (l.to.clone(), l.label.clone())).collect())).collect(),
                    None => vec![("e1".into(), "Experiment 1".into(), None, 0, format!("{} · {}D", physics_label(&self.active().unwrap_or_default()), self.doc.mesh.as_ref().map(|m| m.dim).unwrap_or(2)), vec![])],
                };
                for (eid, title, parent, runs, what, links) in exps {
                    let on = eid == w.exp;
                    let indent = if parent.is_some() { 14.0 } else { 0.0 };
                    ui.horizontal(|ui| {
                        ui.add_space(indent);
                        let r = egui::Frame::new().fill(if on { theme::card(&ctx) } else { Color32::TRANSPARENT }).stroke(egui::Stroke::new(if on { 1.5 } else { 0.0 }, theme::accent(&ctx))).corner_radius(10).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| ui.vertical(|ui| {
                            ui.set_width(ui.available_width());
                            ui.label(RichText::new(&title).strong());
                            let sub = match &parent {
                                Some(p) => format!("fork of {p} · {runs} runs"),
                                None => format!("{runs} runs · {what}"),
                            };
                            ui.label(kit::muted(ui, sub));
                            for (to, label) in &links {
                                ui.label(kit::muted(ui, format!("{label} {to}")).size(12.0));
                            }
                        }));
                        if ui.interact(r.response.rect, ui.scope_id().with(("exp", &eid)), Sense::click()).clicked() && !on && !locked {
                            acts.push(Act::LoadExp(eid.clone()));
                        }
                    });
                }
                ui.add_space(8.0);
                ui.separator();
                let n_runs: usize = w.dir.as_ref().map(|d| d.wf.experiments.iter().map(|e| e.runs.count).sum()).unwrap_or(0);
                let m = ui.selectable_label(w.sub == Sub::Model, RichText::new(format!("Model m1  (not yet functional)\n{n_runs} runs recorded")).size(13.0));
                if m.clicked() {
                    acts.push(Act::Sub(if w.sub == Sub::Model { Sub::Main } else { Sub::Model }));
                }
                if ui.selectable_label(w.notes_open, "Notes").clicked() {
                    w.notes_open = !w.notes_open;
                }
            });
            ui.with_layout(egui::Layout::bottom_up(egui::Align::Min), |ui| {
                ui.label(kit::muted(ui, "Fork copies an experiment and remembers its parent. Link relates any two.").size(12.5));
            });
        });
    }

    // ------------------------------------------------------------------------------------------------------ section bar

    fn wf_sectionbar(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        let ctx = ui.ctx().clone();
        egui::Panel::top("wf_sections").frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(16, 10)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.horizontal(|ui| {
                for (p, label) in [(Part::Concept, "1 · Concept"), (Part::Build, "2 · Build"), (Part::Analysis, "3 · Analysis")] {
                    let on = w.part == p && w.sub == Sub::Main || (w.sub != Sub::Main && p == Part::Analysis && w.part == p);
                    ui.spacing_mut().item_spacing.x = 2.0;
                    if kit::pill(ui, on, label).clicked() && !locked {
                        acts.push(Act::Part(p));
                    }
                }
                ui.spacing_mut().item_spacing.x = 10.0;
                ui.add_space(8.0);
                if locked {
                    kit::dot(ui, theme::AMBER);
                    ui.label(RichText::new("Running").strong());
                    ui.label(kit::muted(ui, "editing is locked until it finishes or you cancel"));
                } else {
                    ui.label(kit::muted(
                        ui,
                        match w.part {
                            Part::Concept => "What are we simulating? Geometry, materials and the essential physics.",
                            Part::Build => "How is it computed? Mesh, boundaries, study and solver.",
                            Part::Analysis => "What came out? Runs, plots, comparisons and fields.",
                        },
                    ));
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if locked {
                        if kit::danger(ui, "Cancel   Esc").clicked() {
                            acts.push(Act::Cancel);
                        }
                    } else {
                        let can = self.busy.is_none() && self.doc.mesh.is_some();
                        if kit::primary_enabled(ui, can, "Compute   F5").on_disabled_hover_text(if self.busy.is_some() { "Something else is running" } else { "There is no mesh to solve on yet" }).clicked() {
                            acts.push(Act::Compute);
                        }
                        if w.part == Part::Analysis && w.selected.len() >= 2 && kit::ghost(ui, format!("Compare {} selected", w.selected.len())).clicked() {
                            acts.push(Act::Sub(Sub::Compare));
                        }
                        if w.part == Part::Build {
                            let ok: Vec<f64> = w.runs.iter().filter(|r| r.status == RunStatus::Ok).map(|r| r.sec).collect();
                            if !ok.is_empty() {
                                ui.label(kit::muted(ui, format!("{} runs · about {:.1} s each · no mesh stored", ok.len(), ok.iter().sum::<f64>() / ok.len() as f64)));
                            }
                        }
                    }
                });
            });
        });
    }

    // ------------------------------------------------------------------------------------------------------ concept and build

    /// The viewport with its two view buttons (geometry, mesh preview).
    fn wf_viewport(&mut self, ui: &mut egui::Ui, strip: &str) {
        let ctx = ui.ctx().clone();
        egui::Panel::bottom("wf_strip").frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(16, 10))).show(ui, |ui| {
            kit::plain_card(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new(strip).color(theme::text(&ctx)));
            });
        });
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(12))).show(ui, |ui| {
            ui.horizontal(|ui| {
                if kit::pill(ui, !self.show_edges, "Geometry").clicked() {
                    self.show_edges = false;
                    self.show_results = false;
                }
                if kit::pill(ui, self.show_edges, "Mesh preview").clicked() {
                    self.show_edges = true;
                    self.show_results = false;
                }
                if let Some(m) = &self.doc.mesh {
                    ui.label(kit::muted(ui, format!("{}D · {} cells · {} points", m.dim, m.n_cells(), m.n_points())));
                } else if let Some(e) = &self.doc.mesh_error {
                    ui.colored_label(theme::RED, e);
                } else {
                    ui.label(kit::muted(ui, "meshing..."));
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.checkbox(&mut self.show_numbers, "Boundary numbers");
                    if ui.small_button("Fit").clicked() {
                        self.vp.fit_token += 1;
                    }
                });
            });
            egui::Frame::new().corner_radius(12).stroke(egui::Stroke::new(1.0, theme::line(&ctx))).show(ui, |ui| {
                self.viewport_ui(ui);
            });
        });
    }

    fn card_line(&mut self, w: &mut OpenWf) -> String {
        let hash = nyx_engine::facade::store::util::hash_bytes(self.spec_json().unwrap_or_default().as_bytes());
        if w.card_for != hash {
            w.card_summary = self.spec_json().ok().and_then(|x| serde_json::from_str::<serde_json::Value>(&x).ok().map(|v| nyx_engine::store_link::card_of(&v))).map(|c| format!("{}  ·  {}", c.summary, c.materials.join("  ·  "))).unwrap_or_default();
            w.card_for = hash;
        }
        w.card_summary.clone()
    }

    fn wf_concept(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, locked: bool) {
        let ctx = ui.ctx().clone();
        let line = self.card_line(w);
        egui::Panel::right("wf_concept_inspector").default_size(440.0).size_range(380.0..=700.0).resizable(true).frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.add_enabled_ui(!locked, |ui| {
                ui.horizontal_wrapped(|ui| {
                    let active = self.active();
                    let codes: Vec<String> = self.engine.project().physics.iter().map(|s| s.physics.clone()).collect();
                    for p in codes {
                        if kit::pill(ui, active.as_deref() == Some(p.as_str()), short_physics(&p)).clicked() {
                            self.cmd("set_active", serde_json::json!({ "physics": p }));
                            self.rules_dirty = true;
                            self.picking = None;
                        }
                    }
                });
                ui.horizontal(|ui| {
                    for (t, l) in [(Tab::Geometry, "Geometry"), (Tab::Materials, "Materials")] {
                        if kit::pill(ui, w.concept_tab == t, l).clicked() {
                            w.concept_tab = t;
                        }
                    }
                });
                ui.separator();
                egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                    self.tab = w.concept_tab;
                    self.inspector_tab(ui, w.concept_tab);
                    ui.add_space(10.0);
                    ui.label(kit::muted(ui, "Meshing, boundaries and solver options live in Build.").size(12.5));
                });
            });
        });
        self.wf_viewport(ui, &format!("WHAT SEARCH WILL SEE    {line}"));
    }

    fn wf_build(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        let ctx = ui.ctx().clone();
        if locked {
            self.wf_live(ui, w, acts);
            return;
        }
        egui::Panel::right("wf_build_inspector").default_size(440.0).size_range(380.0..=700.0).resizable(true).frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.horizontal(|ui| {
                for (t, l) in [(Tab::Physics, "Boundaries and conditions"), (Tab::Study, "Study and solver")] {
                    if kit::pill(ui, w.build_tab == t, l).clicked() {
                        w.build_tab = t;
                    }
                }
            });
            ui.separator();
            egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                self.tab = w.build_tab;
                self.inspector_tab(ui, w.build_tab);
                ui.add_space(12.0);
                kit::card(ui, theme::tint(&ctx), theme::accent(&ctx), 12, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Sweep").strong().size(16.0));
                        ui.label(kit::muted(ui, "coming soon"));
                    });
                    ui.label(kit::muted(ui, "One run per parameter point needs the shared parameter table (see docs/SWEEP_CACHING_PROPOSAL.md). Today one Compute is one run, recorded in Analysis.").size(12.5));
                });
            });
        });
        let strip = match &self.doc.mesh {
            Some(m) => format!("Mesh    {}D · {} cells · {} boundary facets        Compute solves on exactly this mesh; it is rebuilt from the recipe, never stored.", m.dim, m.n_cells(), m.n_facets()),
            None => "Mesh    not built yet".into(),
        };
        self.wf_viewport(ui, &strip);
    }

    /// The live run screen: where the solver is, the runs so far and what it did.
    fn wf_live(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let _ = acts;
        let (status, log) = match &self.job {
            Some(j) => (j.1.status(), j.1.log()),
            None => (Default::default(), vec![]),
        };
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(16))).show(ui, |ui| {
            kit::card(ui, theme::warm(&ctx), theme::AMBER, 16, |ui| {
                ui.set_width(ui.available_width());
                ui.horizontal(|ui| {
                    if self.doc.mesh.is_some() {
                        ui.label(RichText::new("Mesh ready").color(theme::GREEN).strong());
                    }
                    ui.label(kit::muted(ui, "›"));
                    ui.label(RichText::new(if status.stage.is_empty() { "Starting".to_string() } else { status.stage.clone() }).strong().size(16.0));
                });
                ui.add_space(6.0);
                kit::bar(ui, status.fraction.map(|f| f as f32), ui.available_width().min(760.0));
                ui.add_space(4.0);
                let t = status.total_seconds;
                ui.label(kit::mono(format!("{}   elapsed {:02.0}:{:02.0}   {:.1} s in this stage", status.fraction.map(|f| format!("{:.0} %", f * 100.0)).unwrap_or_else(|| "no percentage reported yet".into()), (t / 60.0).floor(), t % 60.0, status.stage_seconds)));
            });
            ui.add_space(12.0);
            ui.columns(2, |cols| {
                kit::plain_card(&mut cols[0], |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("Runs so far").strong());
                    self.runs_plot(ui, w, 220.0);
                    ui.label(kit::muted(ui, "A point is added when a run finishes.").size(12.5));
                });
                kit::plain_card(&mut cols[1], |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("What the solver did").strong().size(16.0));
                    egui::ScrollArea::vertical().max_height(260.0).stick_to_bottom(true).show(ui, |ui| {
                        for l in log.iter().rev().take(60).rev() {
                            ui.label(kit::mono(l.as_str()).size(12.5));
                        }
                        if log.is_empty() {
                            ui.label(kit::muted(ui, "The first line appears when a stage ends."));
                        }
                    });
                    ui.add_space(6.0);
                    ui.label(kit::muted(ui, "Closing this tab does not stop the run. It keeps going and the dashboard card follows it.").size(12.5));
                });
            });
        });
        ctx.request_repaint_after(std::time::Duration::from_millis(250));
    }

    // ------------------------------------------------------------------------------------------------------ analysis

    fn axes(w: &OpenWf) -> (String, String) {
        let ok: Vec<&RunRecord> = w.runs.iter().filter(|r| r.status == RunStatus::Ok).collect();
        let varies = |k: &str| {
            let v: Vec<f64> = ok.iter().filter_map(|r| r.params.get(k).copied()).collect();
            v.windows(2).any(|p| (p[0] - p[1]).abs() > 1e-15 * p[0].abs().max(1e-300))
        };
        let mut pk: Vec<&String> = ok.iter().flat_map(|r| r.params.keys()).collect();
        pk.sort();
        pk.dedup();
        let x = if !w.plot_x.is_empty() { w.plot_x.clone() } else { pk.iter().find(|k| varies(k)).map(|k| k.to_string()).unwrap_or_else(|| "run".into()) };
        let mut mk: Vec<&String> = ok.iter().flat_map(|r| r.metrics.keys()).filter(|k| !matches!(k.as_str(), "dofs" | "n_solutions")).collect();
        mk.sort();
        mk.dedup();
        let y = if !w.plot_y.is_empty() { w.plot_y.clone() } else if mk.iter().any(|k| *k == "neff_re") { "neff_re".into() } else { mk.first().map(|k| k.to_string()).unwrap_or_default() };
        (x, y)
    }

    fn series(w: &OpenWf, only_selected: bool, filter: &Filter) -> (String, String, Vec<(String, f64, f64)>) {
        let (xk, yk) = Self::axes(w);
        let mut pts: Vec<(String, f64, f64)> = vec![];
        for (i, r) in w.runs.iter().enumerate() {
            if r.status != RunStatus::Ok || !filter.matches(r) || (only_selected && !w.selected.contains(&r.id)) {
                continue;
            }
            let x = if xk == "run" { Some(i as f64 + 1.0) } else { r.params.get(&xk).copied() };
            if let (Some(x), Some(y)) = (x, r.metrics.get(&yk).copied()) {
                pts.push((r.id.clone(), x, y));
            }
        }
        pts.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        (xk, yk, pts)
    }

    fn runs_plot(&mut self, ui: &mut egui::Ui, w: &OpenWf, h: f32) {
        let (xk, yk, pts) = Self::series(w, false, &Filter::default());
        let scale = if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { 1e9 } else { 1.0 };
        let xs = pts.iter().map(|p| [p.1 * scale, p.2]).collect::<Vec<_>>();
        let c = theme::accent(ui.ctx());
        egui_plot::Plot::new(("runs_plot", &w.id)).height(h).x_axis_label(if scale == 1e9 { format!("{xk} (nm)") } else { xk }).y_axis_label(yk).allow_scroll(false).show(ui, |p| {
            if xs.len() > 1 {
                p.line(egui_plot::Line::new("runs", egui_plot::PlotPoints::from(xs.clone())).color(c));
            }
            p.points(egui_plot::Points::new("runs", egui_plot::PlotPoints::from(xs)).radius(4.0).color(c));
        });
    }

    fn wf_analysis(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let filter = Filter::parse(&w.filter_text);
        let active_filter = filter.clone().unwrap_or_default();
        egui::Panel::right("wf_analysis_side").exact_size(400.0).resizable(false).frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    let (xk, yk, pts) = Self::series(w, false, &active_filter);
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(format!("{} against {}", if yk.is_empty() { "metric" } else { &yk }, xk)).strong().size(16.0));
                    });
                    let mut keys_x: Vec<String> = vec!["run".into()];
                    let mut keys_y: Vec<String> = vec![];
                    for r in &w.runs {
                        keys_x.extend(r.params.keys().cloned());
                        keys_y.extend(r.metrics.keys().cloned());
                    }
                    keys_x.sort();
                    keys_x.dedup();
                    keys_y.sort();
                    keys_y.dedup();
                    ui.horizontal(|ui| {
                        combo(ui, "x", &mut w.plot_x, &xk, &keys_x);
                        combo(ui, "y", &mut w.plot_y, &yk, &keys_y);
                    });
                    if pts.is_empty() {
                        ui.add_space(20.0);
                        ui.label(kit::muted(ui, "Nothing to plot yet. Each Compute adds one point."));
                        ui.add_space(20.0);
                    } else {
                        self.runs_plot(ui, w, 210.0);
                    }
                });
                ui.add_space(14.0);
                kit::label_caps(ui, "Run a function on these runs", "");
                ui.add_space(4.0);
                let has_pol = self.has_analysis();
                ui.horizontal_wrapped(|ui| {
                    for f in ["Metric vs parameter", "Best run", "Difference of two runs", "Sensitivity"] {
                        if kit::pill(ui, w.function == f, f).clicked() {
                            w.function = f.into();
                        }
                    }
                    if ui.add_enabled(has_pol, egui::Button::new("Polarization on the Poincaré sphere")).on_disabled_hover_text("Needs a result with a polarization analysis in this session: compute a mode study first.").clicked() {
                        acts.push(Act::Field(String::new()));
                    }
                    kit::soon(ui, "Group velocity and dispersion", "needs a wavelength sweep");
                    if ui.add(egui::Button::new("Your script...")).on_hover_text("Runs a Python script on the project and the result on show").clicked() {
                        acts.push(Act::Script);
                    }
                });
                ui.add_space(8.0);
                let (xk, yk, pts) = Self::series(w, false, &active_filter);
                let out = match w.function.as_str() {
                    "Best run" => pts.iter().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal)).map(|p| format!("Highest {yk}: {} at {xk} = {} ({:.6})", p.0, fmt_value(&xk, p.1), p.2)).unwrap_or_else(|| "No runs to rank.".into()),
                    "Sensitivity" => match slope(&pts) {
                        Some(s) => format!("d {yk} / d {xk} ≈ {:.4e} per {} (least squares over {} runs)", if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { s * 1e-9 } else { s }, if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { "nm" } else { "unit" }, pts.len()),
                        None => "Needs at least two runs with different values of the parameter.".into(),
                    },
                    "Difference of two runs" => {
                        if w.selected.len() == 2 { "Open Compare selected for the table of differences.".into() } else { "Tick two runs in the table, then Compare.".into() }
                    }
                    _ => format!("{} runs plotted.", pts.len()),
                };
                ui.label(RichText::new(out).size(13.0));
                ui.add_space(6.0);
                ui.label(kit::muted(ui, "Click a row to open its fields. Tick rows to compare them.").size(12.5));
            });
        });
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(14))).show(ui, |ui| {
            ui.horizontal(|ui| {
                let hint = "params.w between 450 nm and 550 nm and metrics.neff_re > 2.4";
                let bad = filter.is_err();
                ui.add_sized([(ui.available_width() - 330.0).max(180.0), 30.0], egui::TextEdit::singleline(&mut w.filter_text).hint_text(hint).font(egui::TextStyle::Monospace).text_color(if bad { theme::RED } else { theme::text(&ctx) }));
                let cols_changed = &mut w.hidden_cols;
                let mut keys: Vec<String> = w.runs.iter().flat_map(|r| r.params.keys().map(|k| format!("params.{k}")).chain(r.metrics.keys().map(|k| format!("metrics.{k}")))).collect();
                keys.sort();
                keys.dedup();
                ui.menu_button("Columns ▾", |ui| {
                    for k in &keys {
                        let mut on = !cols_changed.contains(k);
                        if ui.checkbox(&mut on, k).changed() {
                            if on { cols_changed.remove(k); } else { cols_changed.insert(k.clone()); }
                        }
                    }
                    if keys.is_empty() {
                        ui.label("No columns yet.");
                    }
                });
                ui.toggle_value(&mut w.by_version, "Group by version");
                ui.toggle_value(&mut w.show_failed, "Show failed");
                if kit::ghost(ui, "Export CSV").clicked() {
                    acts.push(Act::ExportCsv);
                }
            });
            if let Err(e) = &filter {
                ui.colored_label(theme::RED, e);
            }
            ui.add_space(6.0);
            self.runs_table(ui, w, &active_filter, acts);
        });
    }

    fn runs_table(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, filter: &Filter, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let mut pk: Vec<String> = w.runs.iter().flat_map(|r| r.params.keys().cloned()).collect();
        pk.sort();
        pk.dedup();
        pk.retain(|k| !w.hidden_cols.contains(&format!("params.{k}")));
        let mut mk: Vec<String> = w.runs.iter().flat_map(|r| r.metrics.keys().cloned()).filter(|k| !matches!(k.as_str(), "n_solutions")).collect();
        mk.sort();
        mk.dedup();
        mk.retain(|k| !w.hidden_cols.contains(&format!("metrics.{k}")));
        let mut rows: Vec<(usize, &RunRecord)> = w.runs.iter().enumerate().filter(|(_, r)| filter.matches(r) && (w.show_failed || r.status == RunStatus::Ok)).collect();
        let hidden = w.runs.len() - rows.len();
        let sort = w.sort.clone();
        rows.sort_by(|(_, a), (_, b)| {
            let key = |r: &RunRecord| -> f64 {
                if let Some(k) = sort.strip_prefix("params.") { r.params.get(k).copied().unwrap_or(f64::NAN) } else if let Some(k) = sort.strip_prefix("metrics.") { r.metrics.get(k).copied().unwrap_or(f64::NAN) } else if sort == "version" { r.version as f64 } else if sort == "time" { r.sec } else { parse_iso(&r.t).unwrap_or(0) as f64 }
            };
            let o = if w.by_version { b.version.cmp(&a.version) } else { std::cmp::Ordering::Equal };
            let o = o.then(key(a).partial_cmp(&key(b)).unwrap_or(std::cmp::Ordering::Equal).then(a.id.cmp(&b.id)));
            if w.desc { o.reverse() } else { o }
        });
        if w.by_version {
            // newest version first, runs inside a version in the chosen order
            rows.sort_by_key(|r| std::cmp::Reverse(r.1.version));
        }
        let order: Vec<usize> = rows.iter().map(|r| r.0).collect();
        let n_rows = order.len();
        if w.runs.is_empty() {
            kit::plain_card(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("No runs yet").strong().size(16.0));
                ui.label(kit::muted(ui, "Press Compute (F5). Every Compute is recorded here with its parameters, results and the version of the set-up it used."));
                if w.runs_bad > 0 {
                    ui.colored_label(theme::RED, format!("{} lines of runs.jsonl could not be read and were skipped.", w.runs_bad));
                }
            });
            return;
        }
        let widths = |_k: &str| 96.0;
        let head = |ui: &mut egui::Ui, w: &mut OpenWf, label: &str, key: &str, width: f32| {
            let on = w.sort == key;
            let t = format!("{}{}", label.to_uppercase().replace('_', " "), if on { if w.desc { " ▾" } else { " ▴" } } else { "" });
            if ui.add_sized([width, 22.0], egui::Label::new(RichText::new(t).size(12.0).strong().color(theme::muted(&ctx))).sense(Sense::click())).clicked() {
                if on { w.desc = !w.desc } else { w.sort = key.into(); w.desc = true }
            }
        };
        egui::ScrollArea::horizontal().show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 4.0;
                ui.add_space(52.0);
                head(ui, w, "run", "run", 120.0);
                head(ui, w, "version", "version", 64.0);
                for k in &mk {
                    head(ui, w, k, &format!("metrics.{k}"), widths(k));
                }
                for k in &pk {
                    head(ui, w, k, &format!("params.{k}"), widths(k));
                }
                head(ui, w, "time", "time", 70.0);
            });
            ui.separator();
            let h = 26.0;
            egui::ScrollArea::vertical().auto_shrink([false, false]).max_height(ui.available_height() - 40.0).show_rows(ui, h, n_rows, |ui, range| {
                for &ri in &order[range] {
                    let r = w.runs[ri].clone();
                    let sel = w.selected.contains(&r.id);
                    let fill = if sel { theme::tint(&ctx) } else { Color32::TRANSPARENT };
                    let rect_resp = egui::Frame::new().fill(fill).corner_radius(6).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.spacing_mut().item_spacing.x = 4.0;
                            let mut s = sel;
                            if ui.checkbox(&mut s, "").changed() {
                                if s { w.selected.insert(r.id.clone()); } else { w.selected.remove(&r.id); }
                            }
                            if ui.add(egui::Label::new(RichText::new(if r.star { "*" } else { "·" }).strong().size(16.0).color(if r.star { theme::AMBER } else { theme::muted(&ctx) })).sense(Sense::click())).on_hover_text("Star: kept forever, never cleaned up").clicked() {
                                acts.push(Act::Star(r.id.clone(), !r.star));
                            }
                            let id_text = if r.status == RunStatus::Ok { kit::mono(r.id.split('_').next().unwrap_or(&r.id)) } else { kit::mono(format!("{} ({:?})", r.id.split('_').next().unwrap_or(&r.id), r.status).to_lowercase()).color(theme::RED) };
                            let open = ui.add_sized([120.0, h - 4.0], egui::Label::new(id_text).sense(Sense::click()));
                            if open.clicked() {
                                acts.push(Act::Field(r.id.clone()));
                            }
                            ui.add_sized([64.0, h - 4.0], egui::Label::new(kit::mono(format!("v{}", r.version))));
                            for k in &mk {
                                ui.add_sized([widths(k), h - 4.0], egui::Label::new(kit::mono(r.metrics.get(k).map(|v| fmt_value(k, *v)).unwrap_or_default())));
                            }
                            for k in &pk {
                                ui.add_sized([widths(k), h - 4.0], egui::Label::new(kit::mono(r.params.get(k).map(|v| fmt_value(k, *v)).unwrap_or_default())));
                            }
                            ui.add_sized([70.0, h - 4.0], egui::Label::new(kit::mono(format!("{:.2} s", r.sec))));
                        });
                    });
                    let _ = rect_resp;
                }
            });
        });
        ui.label(kit::muted(ui, format!("{} of {} runs match · {} hidden by the filter{}", n_rows, w.runs.len(), hidden, if w.show_failed { "" } else { " · failed and cancelled runs are hidden unless you ask for them" })).size(12.5));
    }

    // ------------------------------------------------------------------------------------------------------ compare, fields, model

    fn wf_compare(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let runs: Vec<RunRecord> = w.runs.iter().filter(|r| w.selected.contains(&r.id)).cloned().collect();
        let colors = [theme::PLUM, theme::AMBER, Color32::from_rgb(0x3A, 0x4A, 0x9B), Color32::from_rgb(0x2F, 0x7A, 0x66), Color32::from_rgb(0x9B, 0x3A, 0x3A)];
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(18))).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(kit::title("Compare runs", 30.0).color(theme::accent(&ctx)));
                ui.label(kit::muted(ui, format!("{} › {}", w.title(), w.exp)));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if kit::primary(ui, "Back to runs").clicked() {
                        acts.push(Act::Sub(Sub::Main));
                    }
                });
            });
            ui.add_space(8.0);
            if runs.len() < 2 {
                ui.label(kit::muted(ui, "Tick two or more runs in the table first."));
                return;
            }
            let (_, yk) = Self::axes(w);
            if w.compare_metric.is_empty() {
                w.compare_metric = yk.clone();
            }
            ui.columns(2, |cols| {
                let ui = &mut cols[0];
                ui.horizontal_wrapped(|ui| {
                    for (i, r) in runs.iter().enumerate() {
                        kit::card(ui, theme::card(&ctx), colors[i % colors.len()], 10, |ui| {
                            ui.label(RichText::new(r.id.split('_').next().unwrap_or(&r.id)).color(colors[i % colors.len()]).monospace());
                            ui.label(RichText::new(r.metrics.get(&w.compare_metric).map(|v| format!("{v:.4}")).unwrap_or_default()).monospace().size(20.0).strong());
                        });
                    }
                });
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("WHAT DIFFERS").size(12.5).strong());
                    let mut keys: Vec<(String, String)> = vec![];
                    for r in &runs {
                        keys.extend(r.params.keys().map(|k| ("params".to_string(), k.clone())));
                        keys.extend(r.metrics.keys().map(|k| ("metrics".to_string(), k.clone())));
                    }
                    keys.sort();
                    keys.dedup();
                    egui::Grid::new("diff").num_columns(runs.len() + 2).spacing([14.0, 4.0]).show(ui, |ui| {
                        ui.label("");
                        for r in &runs {
                            ui.label(kit::mono(r.id.split('_').next().unwrap_or(&r.id)));
                        }
                        ui.label("");
                        ui.end_row();
                        let row = |ui: &mut egui::Ui, name: &str, vals: Vec<String>| {
                            let same = vals.windows(2).all(|p| p[0] == p[1]);
                            ui.label(kit::mono(name));
                            for v in &vals {
                                ui.label(if same { kit::mono(v.as_str()).color(theme::muted(&ctx)) } else { kit::mono(v.as_str()).strong() });
                            }
                            ui.label(if same { kit::muted(ui, "same") } else { RichText::new("") });
                            ui.end_row();
                        };
                        row(ui, "version", runs.iter().map(|r| format!("v{}", r.version)).collect());
                        for (kind, k) in keys {
                            let vals = runs.iter().map(|r| if kind == "params" { r.params.get(&k).map(|v| fmt_value(&k, *v)).unwrap_or_default() } else { r.metrics.get(&k).map(|v| fmt_value(&k, *v)).unwrap_or_default() }).collect();
                            row(ui, &k, vals);
                        }
                        row(ui, "time", runs.iter().map(|r| format!("{:.2} s", r.sec)).collect());
                    });
                    ui.label(kit::muted(ui, "Values that are the same in every run are greyed.").size(12.5));
                });
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("REPRODUCE").size(12.5).strong());
                    ui.label(kit::muted(ui, "Re-mesh and re-solve these runs on this machine and check the numbers against the recorded ones."));
                    kit::soon(ui, format!("Reproduce {} runs", runs.len()), "needs a re-run queue");
                });
                let ui = &mut cols[1];
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Same quantity from each run").strong().size(16.0));
                        let mut keys: Vec<String> = runs.iter().flat_map(|r| r.metrics.keys().cloned()).collect();
                        keys.sort();
                        keys.dedup();
                        let cur = w.compare_metric.clone();
                        combo(ui, "metric", &mut w.compare_metric, &cur, &keys);
                    });
                    let m = w.compare_metric.clone();
                    egui_plot::Plot::new("compare_plot").height(320.0).allow_scroll(false).y_axis_label(m.clone()).show(ui, |p| {
                        for (i, r) in runs.iter().enumerate() {
                            if let Some(v) = r.metrics.get(&m) {
                                p.points(egui_plot::Points::new(r.id.split('_').next().unwrap_or(&r.id).to_string(), egui_plot::PlotPoints::from(vec![[i as f64, *v]])).radius(7.0).color(colors[i % colors.len()]));
                            }
                        }
                    });
                    ui.label(kit::muted(ui, "Field profiles of past runs are not stored yet (coming soon with stored results), so only recorded numbers are compared. Export PNG: coming soon.").size(12.5));
                });
            });
        });
    }

    fn wf_field(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let live = self.current.is_some();
        egui::Panel::right("wf_field_side").exact_size(380.0).resizable(false).frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            if kit::primary(ui, "Back to runs").clicked() {
                acts.push(Act::Sub(Sub::Main));
            }
            ui.add_space(6.0);
            if live {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                    self.inspector_tab(ui, Tab::Results);
                    if self.has_analysis() {
                        ui.separator();
                        self.inspector_tab(ui, Tab::Analysis);
                    }
                    ui.add_space(6.0);
                    if kit::ghost(ui, "Export the field as CSV").clicked() {
                        self.export_csv();
                    }
                });
            } else {
                ui.label(RichText::new("Fields of this run are not stored").strong());
                ui.label(kit::muted(ui, "Heavy results are kept for the session only; stored fields need the HDF5 build and a results setting (coming soon). Compute again to see the fields."));
            }
        });
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(12))).show(ui, |ui| {
            let label = w.focus.clone().filter(|f| !f.is_empty()).map(|f| format!("Field viewer · {}", f.split('_').next().unwrap_or(&f))).unwrap_or_else(|| "Field viewer".into());
            ui.label(kit::title(label, 26.0).color(theme::accent(&ctx)));
            if live {
                self.show_results = true;
            }
            self.viewport_ui(ui);
        });
    }

    fn wf_model(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(24))).show(ui, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(kit::title("Model m1", 38.0).color(theme::accent(&ctx)));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if kit::ghost(ui, "Back").clicked() {
                            acts.push(Act::Sub(Sub::Main));
                        }
                        kit::soon(ui, "Retrain", "no surrogate trainer exists yet");
                        kit::soon(ui, "Export model", "no model to export");
                    });
                });
                ui.label(kit::muted(ui, "A fast stand-in for the solver, learned from every run in this workflow."));
                ui.add_space(8.0);
                kit::card(ui, theme::warm(&ctx), theme::AMBER, 12, |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("Not functional yet").strong());
                    ui.label("nyx has no surrogate trainer in this version. This screen shows what your recorded runs would hand to one; nothing below is a prediction.");
                });
                ui.add_space(10.0);
                ui.columns(2, |cols| {
                    let ui = &mut cols[0];
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("INPUTS (parameters seen in the runs)").size(12.5).strong());
                        let mut n_runs = 0;
                        let mut n_exp = 0;
                        if let Some(d) = &w.dir {
                            for e in &d.wf.experiments {
                                if e.runs.count > 0 {
                                    n_exp += 1;
                                    n_runs += e.runs.count;
                                }
                                for (k, r) in &e.runs.params {
                                    ui.label(kit::mono(format!("{k:<26}{} – {}", fmt_value(k, r.min), fmt_value(k, r.max))));
                                }
                            }
                            ui.add_space(4.0);
                            ui.label(RichText::new("OUTPUTS (metrics recorded)").size(12.5).strong());
                            for e in &d.wf.experiments {
                                for (k, r) in &e.runs.metrics {
                                    ui.label(kit::mono(format!("{k:<26}{:.4} – {:.4}", r.min, r.max)));
                                }
                            }
                        }
                        ui.add_space(4.0);
                        ui.label(kit::mono(format!("FROM  {n_exp} experiments, {n_runs} runs")));
                    });
                    ui.horizontal(|ui| {
                        for m in ["Gaussian process", "Radial basis", "Small network"] {
                            kit::soon(ui, m, "no model code exists yet");
                        }
                    });
                    let ui = &mut cols[1];
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("Predicted against solved, held-out runs").strong().size(16.0));
                        ui.add_space(70.0);
                        ui.label(kit::muted(ui, "Appears once a model has been trained."));
                        ui.add_space(70.0);
                    });
                });
                ui.add_space(8.0);
                ui.label(kit::muted(ui, "Models will only answer inside the ranges above. Outside them nyx will say so instead of guessing."));
            });
        });
    }
}

fn short_physics(p: &str) -> String {
    match p {
        "dfs" => "Diffusion".into(),
        "sms" => "Solid".into(),
        "ems" => "EM waves".into(),
        other => other.to_string(),
    }
}

fn combo(ui: &mut egui::Ui, label: &str, value: &mut String, shown: &str, options: &[String]) {
    egui::ComboBox::from_id_salt(("combo", label)).selected_text(format!("{label}: {shown}")).show_ui(ui, |ui| {
        for o in options {
            if ui.selectable_label(o == shown, o).clicked() {
                *value = o.clone();
            }
        }
    });
}

/// Least-squares slope of y against x.
fn slope(pts: &[(String, f64, f64)]) -> Option<f64> {
    let n = pts.len() as f64;
    if pts.len() < 2 {
        return None;
    }
    let (mx, my) = (pts.iter().map(|p| p.1).sum::<f64>() / n, pts.iter().map(|p| p.2).sum::<f64>() / n);
    let sxx: f64 = pts.iter().map(|p| (p.1 - mx).powi(2)).sum();
    if sxx <= 0.0 {
        return None;
    }
    Some(pts.iter().map(|p| (p.1 - mx) * (p.2 - my)).sum::<f64>() / sxx)
}

#[allow(dead_code)]
fn _t(_: Instant) {}
