// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The workflow view (ADR 124): the header, the one bar (experiment, physics, what is shown, Script, Compute), the five tabs of the left panel with the viewport always on the right, and what hangs off them (the live run,
//! the page of the runs with its table and plot, compare).
//!
//! The editors are the ones of `tabs.rs`; this file only decides where they sit. While a solve of this workflow runs, editing is locked.
//! Anything the app cannot do yet (sweeps, surrogate models, stored fields) is shown as such, never with made-up numbers.

use crate::app::App;
use crate::kit;
use crate::shell::{OpenWf, Page, Sheet, Sub};
use crate::state::Tab;
use crate::theme;
use eframe::egui::{self, Color32, RichText, Sense};
use nyx_engine::facade::store::runs::{Filter, RunPatch, RunRecord, RunStatus};
use nyx_engine::facade::store::util::{ago, now_secs, parse_iso};
use std::time::Instant;

pub enum Act {
    Tab(Tab),
    Sub(Sub),
    LoadExp(String),
    NewExp,
    Fork,
    Link(String, String, String),
    Compute,
    Sweep(crate::sweepui::Form),
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
    crate::fmt::value(key, v)
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
        self.wf_bar(ui, &mut w, &mut acts, locked);
        match w.sub {
            Sub::Compare => self.wf_compare(ui, &mut w, &mut acts),
            Sub::Runs if !locked => self.wf_runs(ui, &mut w, &mut acts),
            _ => self.wf_model(ui, &mut w, &mut acts, locked),
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
    pub(crate) fn ensure_dir(&mut self, id: &str) -> bool {
        if self.shell.open.get(id).is_some_and(|w| w.dir.is_some()) {
            return true;
        }
        let Ok(spec) = self.spec_json() else { return false };
        self.persist_draft(id, &spec).is_ok()
    }

    fn apply_act(&mut self, id: &str, a: Act) {
        match a {
            Act::Tab(t) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.tab = t;
                    w.sub = Sub::Main;
                }
                self.tab = t;
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
            Act::Sweep(f) => self.start_sweep(id, &f),
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
            Act::Script => self.toggle_nys_script(),
            // a run was clicked (or the polarization asked for): the field of the result on show, in the viewport, with the Results tab beside it
            Act::Field(rid) => {
                if let Some(w) = self.shell.open.get_mut(id) {
                    w.focus = Some(rid);
                    w.tab = Tab::Results;
                    w.sub = Sub::Main;
                }
                self.tab = Tab::Results;
                self.show_results = self.current.is_some();
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
            Ok(()) => self.status = format!("Exported {} to {}", crate::fmt::plural(rows.len(), "run"), path.display()),
            Err(e) => self.set_error(format!("Could not write the file: {e}")),
        }
    }

    // ------------------------------------------------------------------------------------------------------ header

    /// The title and the save state; the tags are chips, and versions, share and the other actions of the workflow are one menu (ADR 124).
    fn wf_header(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        egui::Panel::top("wf_header").frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::symmetric(18, 8)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.horizontal(|ui| {
                let eid = egui::Id::unique(("title_edit", &w.id));
                let focused = ctx.memory(|m| m.has_focus(eid));
                let mut buf = if focused { w.draft_buf.clone() } else { w.title() };
                let width = (buf.chars().count().max(10) as f32 * 11.0).min((ui.available_width() - 340.0).max(200.0));
                let r = ui.add(egui::TextEdit::singleline(&mut buf).id(eid).font(theme::title(30.0)).frame(egui::Frame::NONE).desired_width(width).text_color(theme::text(&ctx)));
                if focused || r.has_focus() {
                    w.draft_buf = buf.clone();
                }
                if r.lost_focus() && buf.trim() != w.title() {
                    acts.push(Act::Title(buf.clone()));
                }
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
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.menu_button("More ▾", |ui| {
                        let v = w.dir.as_ref().map(|d| d.wf.version).unwrap_or(0);
                        if ui.add_enabled(w.dir.is_some(), egui::Button::new(format!("Versions  v{v}"))).on_disabled_hover_text("Versions start with the first change").clicked() {
                            self.shell.sheet = Sheet::Versions;
                            self.shell.version_sel = None;
                            ui.close();
                        }
                        if ui.button("Share / export").clicked() {
                            self.shell.sheet = Sheet::Share;
                            ui.close();
                        }
                        ui.separator();
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
                        ui.separator();
                        ui.label(kit::muted(ui, "Add a tag").size(12.5));
                        let t = ui.add(egui::TextEdit::singleline(&mut w.tag_input).hint_text("name, then Enter").desired_width(150.0));
                        if t.lost_focus() && ctx.input(|i| i.key_pressed(egui::Key::Enter)) && !w.tag_input.trim().is_empty() {
                            tags.push(w.tag_input.trim().to_lowercase());
                            w.tag_input.clear();
                            changed = true;
                        }
                        ui.separator();
                        if ui.button("Close tab (Ctrl+W)").clicked() {
                            acts.push(Act::Close);
                            ui.close();
                        }
                    });
                    // the words of the save state give way to one word when the title leaves little room; the whole sentence is the tooltip
                    let (text, short) = match (&w.dir, w.saved_at) {
                        (None, _) => ("Draft: nothing is saved until you change something".to_string(), "Draft"),
                        (Some(_), Some(t)) => (format!("Saved automatically · {}", ago(now_secs() - t.elapsed().as_secs() as i64, now_secs())), "Saved"),
                        (Some(d), None) => (format!("Saved · {}", ago(d.wf.touched_secs(), now_secs())), "Saved"),
                    };
                    let wide = ui.painter().layout_no_wrap(text.clone(), egui::TextStyle::Body.resolve(ui.style()), theme::muted(&ctx)).size().x;
                    let shown = if wide + 24.0 <= ui.available_width() { text.clone() } else { short.to_string() };
                    ui.label(kit::muted(ui, shown)).on_hover_text(text);
                    kit::glyph(ui, if w.dir.is_some() { kit::State::Done } else { kit::State::Idle });
                });
                if changed {
                    acts.push(Act::Tags(tags));
                }
            });
        });
    }

    // ------------------------------------------------------------------------------------------------------ the bar

    /// The experiments of the workflow as one menu of the bar (ADR 124): choose one, start a blank one, fork this one, link two. It replaces the rail on the left.
    fn experiment_menu(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        let exps: Vec<(String, String, Option<String>, usize)> = match &w.dir {
            Some(d) => d.wf.experiments.iter().map(|e| (e.id.clone(), e.title.clone(), e.parent.clone(), e.runs.count)).collect(),
            None => vec![("e1".into(), "Experiment 1".into(), None, 0)],
        };
        let current = exps.iter().find(|e| e.0 == w.exp).map(|e| e.1.clone()).unwrap_or_else(|| "Experiment 1".into());
        ui.add_enabled_ui(!locked, |ui| {
            ui.menu_button(format!("{current} ▾"), |ui| {
                for (eid, title, parent, runs) in &exps {
                    let label = format!("{}{}  ·  {}", if parent.is_some() { "↳ " } else { "" }, title, crate::fmt::plural(*runs, "run"));
                    if ui.selectable_label(*eid == w.exp, label).clicked() && *eid != w.exp {
                        acts.push(Act::LoadExp(eid.clone()));
                        ui.close();
                    }
                }
                ui.separator();
                if ui.button("New blank experiment").clicked() {
                    acts.push(Act::NewExp);
                    ui.close();
                }
                if ui.button("Fork this experiment").on_hover_text("Copies the set-up and remembers its parent").clicked() {
                    acts.push(Act::Fork);
                    ui.close();
                }
                let others: Vec<(String, String)> = exps.iter().filter(|e| e.0 != w.exp).map(|e| (e.0.clone(), e.1.clone())).collect();
                if w.dir.is_some() && !others.is_empty() {
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
            });
        });
    }

    /// The physics of the project as one dropdown (ADR 124): choose the one on show, add another, remove the one on show. A switch changes what the window shows, not the
    /// workflow (`shell::spec_key`).
    fn physics_menu(&mut self, ui: &mut egui::Ui, locked: bool) {
        let in_project: Vec<String> = self.engine.project().physics.iter().map(|s| s.physics.clone()).collect();
        let free: Vec<String> = self.engine.registry().ids().iter().map(|s| s.to_string()).filter(|i| !in_project.contains(i)).collect();
        let active = self.active().unwrap_or_default();
        let shown = if active.is_empty() { "none".to_string() } else { self.physics_name(&active) };
        ui.add_enabled_ui(!locked, |ui| {
            egui::ComboBox::from_id_salt("wf_physics").selected_text(format!("Physics: {shown}")).show_ui(ui, |ui| {
                for p in &in_project {
                    if ui.selectable_label(active == *p, self.physics_name(p)).clicked() && active != *p {
                        self.cmd("set_active", serde_json::json!({ "physics": p }));
                        self.picking = None;
                        self.rules_dirty = true;
                    }
                }
                if !free.is_empty() && in_project.len() < 3 {
                    ui.separator();
                    for id in &free {
                        if ui.selectable_label(false, format!("Add {}", self.physics_name(id))).clicked() {
                            self.cmd("add_physics", serde_json::json!({ "physics": id }));
                            self.cmd("set_active", serde_json::json!({ "physics": id }));
                            self.picking = None;
                            self.rules_dirty = true;
                        }
                    }
                }
                if !active.is_empty() {
                    ui.separator();
                    if ui.selectable_label(false, format!("Remove {} from the project", self.physics_name(&active))).clicked() {
                        self.cmd("remove_physics", serde_json::json!({ "physics": active }));
                        self.picking = None;
                        self.rules_dirty = true;
                    }
                }
            });
        });
    }

    /// The one bar of the workflow (ADR 124): the experiment and the physics on the left; on the right what is shown (the geometry, the mesh, the result, the runs),
    /// the script and the one filled button, Compute. While a run of this workflow goes, only what it is doing and Cancel.
    fn wf_bar(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        let ctx = ui.ctx().clone();
        egui::Panel::top("wf_bar").frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(16, 8)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            // in a narrow window (or beside the script dock) what is shown goes to a row of its own: nothing overlaps
            let narrow = ui.available_width() < 960.0;
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 10.0;
                self.experiment_menu(ui, w, acts, locked);
                self.physics_menu(ui, locked);
                if locked {
                    kit::glyph(ui, kit::State::Running);
                    ui.label(RichText::new("Running").strong());
                    if !narrow {
                        ui.label(kit::muted(ui, "editing is locked until it finishes or you cancel"));
                    }
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    self.bar_actions(ui, acts, locked);
                    if !narrow && !locked {
                        self.bar_show(ui, w, acts);
                    }
                });
            });
            if narrow && !locked {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 10.0;
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| self.bar_show(ui, w, acts));
                });
            }
        });
    }

    /// Compute (the one filled button) and Script, or Cancel while this workflow runs. Drawn right to left: the last one added is the leftmost.
    fn bar_actions(&mut self, ui: &mut egui::Ui, acts: &mut Vec<Act>, locked: bool) {
        if locked {
            if kit::danger(ui, "Cancel   Esc").clicked() {
                acts.push(Act::Cancel);
            }
            return;
        }
        let can = self.busy.is_none() && self.doc.mesh.is_some();
        if kit::primary_enabled(ui, can, "Compute   F5").on_disabled_hover_text(if self.busy.is_some() { "Something else is running" } else { "There is no mesh to solve on yet" }).clicked() {
            acts.push(Act::Compute);
        }
        if kit::ghost(ui, "Script").on_hover_text("Write, load, check or run a .nys script in a dock on the right (Ctrl+Shift+R)").clicked() {
            acts.push(Act::Script);
        }
    }

    /// What is shown: the three views of the model, and the runs from the second run on. Drawn right to left (the runs first, `Show` last).
    fn bar_show(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        ui.add_space(6.0);
        let runs = w.runs.iter().filter(|r| r.status == RunStatus::Ok).count();
        if runs >= 2 || matches!(w.sub, Sub::Runs | Sub::Compare) {
            let on = matches!(w.sub, Sub::Runs | Sub::Compare);
            if ui.selectable_label(on, format!("Runs ({runs})")).on_hover_text("The table, the plot and the comparison of the runs of this experiment").clicked() {
                acts.push(Act::Sub(Sub::Runs));
            }
        }
        let in_model = w.sub == Sub::Main;
        let has_result = !self.results.is_empty();
        if ui.add_enabled(has_result, egui::Button::selectable(in_model && self.show_results, "Result")).on_disabled_hover_text("Compute first: there is no result to show").clicked() {
            self.show_results = true;
            w.sub = Sub::Main;
            w.tab = Tab::Results;
            self.tab = Tab::Results;
        }
        if ui.selectable_label(in_model && !self.show_results && self.show_edges, "Mesh").on_hover_text("The mesh the study solves on").clicked() {
            self.show_edges = true;
            self.show_results = false;
            w.sub = Sub::Main;
        }
        if ui.selectable_label(in_model && !self.show_results && !self.show_edges, "Geometry").on_hover_text("The regions of the model").clicked() {
            self.show_edges = false;
            self.show_results = false;
            w.sub = Sub::Main;
        }
        ui.label(kit::muted(ui, "Show"));
    }

    // ------------------------------------------------------------------------------------------------------ the model

    /// The model page (ADR 124): the tabs of what is edited and read on the left (Geometry, Materials, Physics, Study, Results), the viewport always on the right.
    fn wf_model(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>, locked: bool) {
        if locked {
            self.wf_live(ui, w, acts);
            return;
        }
        self.wf_inspector(ui, w, acts);
        let strip = match &self.doc.mesh {
            Some(m) => format!("Mesh · {}D · {} points · {} cells · {} boundary facets", m.dim, m.n_points(), m.n_cells(), m.n_facets()),
            None => "Mesh · not built yet".into(),
        };
        self.wf_viewport(ui, &strip);
    }

    /// The left panel: one row of five tabs, then the editor of the tab. The sweep is a part of Study; the polarization is a part of Results.
    fn wf_inspector(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        egui::Panel::left("wf_inspector").default_size(400.0).size_range(320.0..=640.0).resizable(true).frame(egui::Frame::new().fill(theme::card(&ctx)).inner_margin(egui::Margin::same(14)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            const TABS: [Tab; 5] = [Tab::Geometry, Tab::Materials, Tab::Physics, Tab::Study, Tab::Results];
            let on = TABS.iter().position(|t| *t == w.tab).unwrap_or(0);
            let labels: Vec<&str> = TABS.iter().map(|t| t.name()).collect();
            if let Some(i) = kit::tabs(ui, &labels, on) {
                acts.push(Act::Tab(TABS[i]));
            }
            ui.add_space(6.0);
            ui.separator();
            // both axes: a row of a tab that is wider than the panel scrolls; it does not widen the panel (a vertical-only area takes the width of its content)
            egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                self.tab = w.tab;
                match w.tab {
                    Tab::Study => {
                        self.inspector_tab(ui, Tab::Study);
                        ui.add_space(10.0);
                        ui.separator();
                        let head = format!("{}  Sweep", if w.sweep.shown { "▾" } else { "▸" });
                        if ui.add(egui::Label::new(RichText::new(head).strong().size(16.0)).sense(Sense::click())).on_hover_text("Run one parameter over a range").clicked() {
                            w.sweep.shown = !w.sweep.shown;
                        }
                        if w.sweep.shown {
                            self.sweep_section(ui, w, acts);
                        }
                    }
                    Tab::Results => {
                        crate::spans::span("inspector_tab", || self.inspector_tab(ui, Tab::Results));
                        if self.has_analysis() {
                            ui.separator();
                            self.inspector_tab(ui, Tab::Analysis);
                        }
                        if self.current.is_some() {
                            ui.add_space(6.0);
                            if kit::ghost(ui, "Export the field as CSV").clicked() {
                                self.export_csv();
                            }
                        }
                    }
                    t => self.inspector_tab(ui, t),
                }
            });
        });
    }

    /// The viewport, always there. Under it, when the result has several modes, the strip of the solutions; under that the status line with the facts of the mesh.
    fn wf_viewport(&mut self, ui: &mut egui::Ui, strip: &str) {
        let ctx = ui.ctx().clone();
        // a status bar of one line, not a card (ADR 123): the facts of the model on show, muted; what it means is in the tooltip
        let status = match self.current.filter(|_| self.show_results) {
            Some(i) => format!("{strip}  ·  result: {}", self.results[i].name),
            None => strip.to_string(),
        };
        egui::Panel::bottom("wf_strip").frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(16, 4)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| {
            ui.add(egui::Label::new(kit::muted(ui, status).size(12.5)).truncate())
                .on_hover_text("Compute solves on exactly this mesh; it is rebuilt from the recipe, never stored.");
        });
        if self.bottom_open && self.show_results && self.current.is_some_and(|i| self.results[i].sols.len() > 1) {
            egui::Panel::bottom("wf_modes").default_size(180.0).resizable(true).frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(16, 8)).stroke(egui::Stroke::new(1.0, theme::line(&ctx)))).show(ui, |ui| self.results_panel(ui));
        }
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(12))).show(ui, |ui| {
            ui.horizontal(|ui| {
                if ui.small_button("Fit").clicked() {
                    self.vp.fit_token += 1;
                }
                if !self.show_results {
                    ui.checkbox(&mut self.show_numbers, "Boundary numbers");
                }
            });
            egui::Frame::new().corner_radius(12).stroke(egui::Stroke::new(1.0, theme::line(&ctx))).show(ui, |ui| {
                crate::spans::span("viewport_ui", || self.viewport_ui(ui));
            });
        });
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
            kit::card(ui, theme::warm(&ctx), theme::activity(&ctx), 16, |ui| {
                ui.set_width(ui.available_width());
                ui.horizontal(|ui| {
                    if self.doc.mesh.is_some() {
                        ui.label(RichText::new("Mesh ready").color(theme::success(ui.ctx())).strong());
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
        // the default metric is the first that varies between the runs (never a constant such as an imaginary part that is zero, nor the run time)
        let mvaries = |k: &str| {
            let v: Vec<f64> = ok.iter().filter_map(|r| r.metrics.get(k).copied()).collect();
            v.windows(2).any(|p| (p[0] - p[1]).abs() > 1e-9 * p[0].abs().max(p[1].abs()).max(1e-300))
        };
        let y = if !w.plot_y.is_empty() {
            w.plot_y.clone()
        } else if mk.iter().any(|k| *k == "neff_re") {
            "neff_re".into()
        } else {
            mk.iter().find(|k| k.as_str() != "seconds" && mvaries(k)).or(mk.first()).map(|k| k.to_string()).unwrap_or_default()
        };
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

    /// The runs plot of the Analysis page and of the live run screen (UX plan B5): the line joins the runs in order of the parameter, **a marker
    /// shape tells the versions of the set-up apart** (colour is never the only channel; a legend names them when there are two or more), **the
    /// last run is labelled at the end of the line**, and **the runs ticked in the table carry a ring** in the text colour.
    fn runs_plot(&mut self, ui: &mut egui::Ui, w: &OpenWf, h: f32) {
        use egui_plot::MarkerShape as M;
        const SHAPES: [M; 5] = [M::Circle, M::Square, M::Diamond, M::Up, M::Cross];
        let (xk, yk, pts) = Self::series(w, false, &Filter::default());
        let scale = if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { 1e9 } else { 1.0 };
        let c = theme::accent(ui.ctx());
        let ring = theme::text(ui.ctx());
        let version_of = |id: &str| w.runs.iter().find(|r| r.id == id).map(|r| r.version).unwrap_or(0);
        let mut versions: Vec<u32> = pts.iter().map(|p| version_of(&p.0)).collect();
        versions.sort_unstable();
        versions.dedup();
        let line: Vec<[f64; 2]> = pts.iter().map(|p| [p.1 * scale, p.2]).collect();
        let selected: Vec<[f64; 2]> = pts.iter().filter(|p| w.selected.contains(&p.0)).map(|p| [p.1 * scale, p.2]).collect();
        let end = pts.last().map(|p| (egui_plot::PlotPoint::new(p.1 * scale, p.2), format!("{}  {}", p.0.split('_').next().unwrap_or(&p.0), fmt_value(&yk, p.2))));
        let mut plot = egui_plot::Plot::new(("runs_plot", &w.id)).height(h).x_axis_label(if scale == 1e9 { format!("{xk} (nm)") } else { xk }).y_axis_label(yk).allow_scroll(false);
        // headroom above the highest point, so the end label is not cut by the frame of the plot
        let (lo, hi) = pts.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), p| (a.min(p.2), b.max(p.2)));
        if hi.is_finite() && lo.is_finite() {
            plot = plot.include_y(hi + 0.18 * (hi - lo).max(hi.abs() * 1e-3));
        }
        if versions.len() > 1 {
            plot = plot.legend(egui_plot::Legend::default());
        }
        plot.show(ui, |p| {
            if line.len() > 1 {
                p.line(egui_plot::Line::new("runs", egui_plot::PlotPoints::from(line.clone())).color(c));
            }
            for (k, v) in versions.iter().enumerate() {
                let of_v: Vec<[f64; 2]> = pts.iter().filter(|q| version_of(&q.0) == *v).map(|q| [q.1 * scale, q.2]).collect();
                p.points(egui_plot::Points::new(format!("v{v}"), egui_plot::PlotPoints::from(of_v)).radius(4.5).color(c).shape(SHAPES[k % SHAPES.len()]));
            }
            if !selected.is_empty() {
                p.points(egui_plot::Points::new("ticked", egui_plot::PlotPoints::from(selected)).radius(9.0).color(ring).filled(false).shape(M::Circle));
            }
            if let Some((at, text)) = end {
                p.text(egui_plot::Text::new("last run", at, text).anchor(egui::Align2::RIGHT_BOTTOM).color(ring));
            }
        });
    }

    /// The runs of the experiment (ADR 124): the plot from two runs on, the functions on the runs behind a header, and the table. It is a page of its own beside the model
    /// (`Runs` in the bar); the field of a run is the Results tab and the viewport, not a page here.
    fn wf_runs(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let filter = Filter::parse(&w.filter_text);
        let active_filter = filter.clone().unwrap_or_default();
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(14))).show(ui, |ui| {
            ui.horizontal(|ui| {
                if kit::ghost(ui, "← Model").on_hover_text("Back to the model and its viewport").clicked() {
                    acts.push(Act::Sub(Sub::Main));
                }
                ui.label(RichText::new(format!("Runs of {}", w.exp)).strong().size(15.0));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if w.selected.len() >= 2 && kit::ghost(ui, format!("Compare {} selected", w.selected.len())).clicked() {
                        acts.push(Act::Sub(Sub::Compare));
                    }
                });
            });
            ui.add_space(6.0);
            let (xk, yk, pts) = Self::series(w, false, &active_filter);
            let flat = pts.windows(2).all(|p| (p[0].2 - p[1].2).abs() <= 1e-12 * p[0].2.abs().max(1e-300)) && xk == "run";
            if pts.len() >= 2 && !flat {
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
                    ui.label(RichText::new(format!("{} against {}", if yk.is_empty() { "metric" } else { &yk }, xk)).strong());
                    combo(ui, "x", &mut w.plot_x, &xk, &keys_x);
                    combo(ui, "y", &mut w.plot_y, &yk, &keys_y);
                });
                self.runs_plot(ui, w, 200.0);
            } else if pts.len() >= 2 {
                ui.label(kit::muted(ui, "The runs are the same so far, so there is nothing to plot: change a parameter, or sweep one (Study tab)."));
            } else if !w.runs.is_empty() {
                ui.label(kit::muted(ui, "A plot needs two runs. Each Compute adds one."));
            }
            ui.add_space(4.0);
            egui::CollapsingHeader::new("Functions on these runs").default_open(false).show(ui, |ui| {
                let has_pol = self.has_analysis();
                ui.horizontal_wrapped(|ui| {
                    for f in ["Metric vs parameter", "Best run", "Difference of two runs", "Sensitivity"] {
                        if ui.selectable_label(w.function == f, f).clicked() {
                            w.function = f.into();
                        }
                    }
                    if ui.add_enabled(has_pol, egui::Button::new("Polarization on the Poincaré sphere")).on_disabled_hover_text("Needs a result with a polarization analysis in this session: compute a mode study first.").clicked() {
                        acts.push(Act::Field(String::new()));
                    }
                    if ui.add(egui::Button::new("Your script...")).on_hover_text("Write, check or run a .nys script (or an external Python script) on this project and result").clicked() {
                        acts.push(Act::Script);
                    }
                });
                let out = match w.function.as_str() {
                    "Best run" => pts.iter().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal)).map(|p| format!("Highest {yk}: {} at {xk} = {} ({:.6})", p.0, fmt_value(&xk, p.1), p.2)).unwrap_or_else(|| "No runs to rank.".into()),
                    "Sensitivity" => match slope(&pts) {
                        Some(s) => format!("d {yk} / d {xk} ≈ {:.4e} per {} (least squares over {})", if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { s * 1e-9 } else { s }, if xk != "run" && pts.iter().all(|p| p.1.abs() < 1e-3) { "nm" } else { "unit" }, crate::fmt::plural(pts.len(), "run")),
                        None => "Needs at least two runs with different values of the parameter.".into(),
                    },
                    "Difference of two runs" => {
                        if w.selected.len() == 2 { "Open Compare selected for the table of differences.".into() } else { "Tick two runs in the table, then Compare.".into() }
                    }
                    _ => format!("{} plotted.", crate::fmt::plural(pts.len(), "run")),
                };
                ui.label(RichText::new(out).size(13.0));
            });
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                let hint = "params.w between 450 nm and 550 nm and metrics.neff_re > 2.4";
                let bad = filter.is_err();
                ui.add_sized([(ui.available_width() - 330.0).max(180.0), 30.0], egui::TextEdit::singleline(&mut w.filter_text).hint_text(hint).font(egui::TextStyle::Monospace).text_color(if bad { theme::danger(ui.ctx()) } else { theme::text(&ctx) }));
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
                ui.colored_label(theme::danger(ui.ctx()), e);
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
                    ui.colored_label(theme::danger(ui.ctx()), format!("{} lines of runs.jsonl could not be read and were skipped.", w.runs_bad));
                }
            });
            return;
        }
        let widths = |_k: &str| 124.0;
        let head = |ui: &mut egui::Ui, w: &mut OpenWf, label: &str, key: &str, width: f32| {
            let on = w.sort == key;
            let t = format!("{}{}", crate::fmt::key_label(label), if on { if w.desc { " ▾" } else { " ▴" } } else { "" });
            if ui.add_sized([width, 22.0], egui::Label::new(RichText::new(t).size(12.0).strong().color(theme::muted(&ctx))).truncate().sense(Sense::click())).on_hover_text(label).clicked() {
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
                            if ui.add(egui::Label::new(RichText::new(if r.star { "*" } else { "·" }).strong().size(16.0).color(if r.star { theme::activity_text(ui.ctx()) } else { theme::muted(&ctx) })).sense(Sense::click())).on_hover_text("Star: kept forever, never cleaned up").clicked() {
                                acts.push(Act::Star(r.id.clone(), !r.star));
                            }
                            let id_text = if r.status == RunStatus::Ok { kit::mono(r.id.split('_').next().unwrap_or(&r.id)) } else { kit::mono(format!("{} ({:?})", r.id.split('_').next().unwrap_or(&r.id), r.status).to_lowercase()).color(theme::danger(ui.ctx())) };
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
        ui.label(kit::muted(ui, format!("{} of {} match · {} hidden by the filter{}", n_rows, crate::fmt::plural(w.runs.len(), "run"), hidden, if w.show_failed { "" } else { " · failed and cancelled runs are hidden unless you ask for them" })).size(12.5));
    }

    // ------------------------------------------------------------------------------------------------------ compare, fields, model

    fn wf_compare(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<Act>) {
        let ctx = ui.ctx().clone();
        let runs: Vec<RunRecord> = w.runs.iter().filter(|r| w.selected.contains(&r.id)).cloned().collect();
        let colors = [theme::accent(ui.ctx()), theme::activity(ui.ctx()), Color32::from_rgb(0x3A, 0x4A, 0x9B), Color32::from_rgb(0x2F, 0x7A, 0x66), Color32::from_rgb(0x9B, 0x3A, 0x3A)];
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::same(18))).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(kit::title("Compare runs", 30.0));
                ui.label(kit::muted(ui, format!("{} › {}", w.title(), w.exp)));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if kit::primary(ui, "Back to runs").clicked() {
                        acts.push(Act::Sub(Sub::Runs));
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
                    ui.label(kit::muted(ui, "Only the recorded numbers of each run are compared; field profiles of past runs are not stored.").size(12.5));
                });
            });
        });
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
