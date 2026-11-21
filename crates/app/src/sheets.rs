// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The sheets that open over the shell (new workflow, search palette, versions, share, open a received file) and the settings page.

use crate::app::App;
use crate::examples;
use crate::kit;
use crate::shell::{Page, Sheet, Start};
use crate::theme;
use eframe::egui::{self, RichText};
use nyx_engine::facade::store::archive::{ExportOptions, Limits};
use nyx_engine::facade::store::util::{ago, now_secs};
use nyx_engine::facade::store::versions::{self, VersionKind};
use nyx_engine::facade::store::ShareLevel;

/// State of the in-window `.nys` script editor and runner sheet.
#[derive(Default)]
pub struct ScriptSheet {
    pub open: bool,
    pub code: String,
    pub path: Option<std::path::PathBuf>,
    pub check_msg: Option<Result<String, String>>,
    pub output: Option<Result<String, String>>,
    pub running: bool,
}

/// Start-from cards: (title, text, index in `examples::list()`, or `None` when it starts blank)
const TEMPLATES: [(&str, &str, Option<usize>); 7] = [
    ("Optical waveguide", "Strip or rib cross-section, modes and n_eff", Some(14)),
    ("Strip on sapphire", "Optical and acoustic modes of a silicon strip", Some(12)),
    ("Rectangular cavity", "3D electromagnetic eigenmodes of a box", Some(13)),
    ("Vibrating structure", "Solid mechanics, eigenfrequencies", Some(1)),
    ("Heat transfer", "Stationary or transient temperature", Some(2)),
    ("Electrostatics", "Potentials, fields, capacitance", Some(0)),
    ("Blank experiment", "Choose physics, dimension and geometry yourself", None),
];

impl App {
    pub fn sheets_ui(&mut self, ctx: &egui::Context) {
        if self.nys.open {
            if kit::modal(ctx, "nys_script", 920.0, |ui| self.sheet_script(ui)) {
                self.nys.open = false;
            }
            return;
        }
        let close = match self.shell.sheet {
            Sheet::None => return,
            Sheet::New => kit::modal(ctx, "new", 940.0, |ui| self.sheet_new(ui)),
            Sheet::Palette => kit::modal(ctx, "palette", 860.0, |ui| self.sheet_palette(ui)),
            Sheet::Versions => kit::modal(ctx, "versions", 1000.0, |ui| self.sheet_versions(ui)),
            Sheet::Share => kit::modal(ctx, "share", 820.0, |ui| self.sheet_share(ui)),
            Sheet::Open => kit::modal(ctx, "open", 820.0, |ui| self.sheet_open(ui)),
            Sheet::About => kit::modal(ctx, "about", 760.0, |ui| self.sheet_about(ui)),
        };
        if close {
            self.shell.sheet = Sheet::None;
        }
    }

    fn sheet_script(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.label(kit::title("Script (.nys)", 34.0).color(theme::accent(&ctx)));
        ui.label(kit::muted(ui, "Python-syntax script run by nyx-script-worker on the open setup. Changes to geometry, materials, mesh and solves update this workflow directly."));
        ui.add_space(8.0);
        ui.horizontal_wrapped(|ui| {
            if kit::ghost(ui, "Load .nys...").clicked() {
                if let Some(p) = rfd::FileDialog::new().add_filter("nyx script", &["nys"]).pick_file() {
                    match std::fs::read_to_string(&p) {
                        Ok(text) => {
                            self.nys.code = text.trim_start_matches('\u{feff}').to_string();
                            self.nys.path = Some(p);
                            self.nys.check_msg = None;
                        }
                        Err(e) => self.set_error(format!("Could not read script: {e}")),
                    }
                }
            }
            if kit::ghost(ui, "Save .nys...").clicked() {
                let default_name = self.nys.path.as_ref().and_then(|p| p.file_name()).and_then(|n| n.to_str()).unwrap_or("script.nys").to_string();
                if let Some(p) = rfd::FileDialog::new().add_filter("nyx script", &["nys"]).set_file_name(default_name).save_file() {
                    match std::fs::write(&p, &self.nys.code) {
                        Ok(()) => {
                            self.status = format!("Saved {}", p.display());
                            self.nys.path = Some(p);
                        }
                        Err(e) => self.set_error(format!("Could not write script: {e}")),
                    }
                }
            }
            if kit::ghost(ui, "Fill from model").on_hover_text("Replace the editor contents with the .nys commands that build the current model").clicked() {
                self.fill_nys_from_model();
            }
            if kit::ghost(ui, "Check").on_hover_text("Parse and validate every nyx command and argument without running").clicked() {
                self.check_nys_script();
            }
            let can_run = !self.nys.running && self.busy.is_none();
            if kit::primary_enabled(ui, can_run, if self.nys.running { "Running..." } else { "Run .nys script" }).clicked() {
                self.run_nys_script();
            }
            if kit::ghost(ui, "External .py on CSV...").on_hover_text("Run an external Python script on the exported project JSON and result CSV").clicked() {
                self.run_python_script();
            }
        });
        if let Some(p) = &self.nys.path {
            ui.label(kit::mono(p.display().to_string()).size(12.0));
        }
        ui.add_space(6.0);
        egui::ScrollArea::vertical().id_salt("nys_code_scroll").max_height(280.0).show(ui, |ui| {
            ui.add(
                egui::TextEdit::multiline(&mut self.nys.code)
                    .code_editor()
                    .desired_rows(15)
                    .desired_width(f32::INFINITY)
                    .lock_focus(true)
                    .hint_text("# Write .nys commands here, e.g.\n# mesh_build()\n# solve()"),
            );
        });
        if self.nys.running {
            ui.add_space(6.0);
            kit::status(ui, kit::State::Running, "Running script in nyx-script-worker...");
        }
        if let Some(chk) = &self.nys.check_msg {
            ui.add_space(6.0);
            match chk {
                Ok(msg) => kit::status(ui, kit::State::Done, msg.clone()),
                Err(err) => {
                    kit::card(ui, theme::raised(&ctx), theme::danger(&ctx), 10, |ui| {
                        ui.set_width(ui.available_width());
                        ui.colored_label(theme::danger(&ctx), RichText::new("Static check failed").strong());
                        ui.label(kit::mono(err.as_str()).size(12.0));
                    });
                }
            }
        }
        if let Some(out) = &self.nys.output {
            ui.add_space(6.0);
            match out {
                Ok(text) => {
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("Script output").strong().color(theme::success(&ctx)));
                        egui::ScrollArea::vertical().id_salt("nys_out_scroll").max_height(180.0).show(ui, |ui| {
                            ui.add(egui::TextEdit::multiline(&mut text.as_str()).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
                        });
                    });
                }
                Err(err) => {
                    kit::card(ui, theme::raised(&ctx), theme::danger(&ctx), 10, |ui| {
                        ui.set_width(ui.available_width());
                        ui.colored_label(theme::danger(&ctx), RichText::new("Script run failed").strong());
                        egui::ScrollArea::vertical().id_salt("nys_err_scroll").max_height(180.0).show(ui, |ui| {
                            ui.label(kit::mono(err.as_str()).size(12.0));
                        });
                    });
                }
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------ new workflow

    fn sheet_new(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.label(kit::title("Start a workflow", 36.0).color(theme::accent(&ctx)));
                ui.label(kit::muted(ui, "Pick a starting point. You can add more experiments to it later, and nothing is saved until you change something."));
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                ui.add(egui::TextEdit::singleline(&mut self.shell.new.title).hint_text("Name (optional)").desired_width(260.0));
            });
        });
        ui.add_space(8.0);
        kit::label_caps(ui, "Start from", "");
        ui.add_space(4.0);
        let sel = self.shell.new.template;
        let w = (ui.available_width() - 3.0 * 12.0) / 4.0;
        for row in 0..TEMPLATES.len().div_ceil(4) {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 12.0;
                for col in 0..4 {
                    let i = row * 4 + col;
                    let Some(&(title, text, _)) = TEMPLATES.get(i) else { break };
                    let on = sel == i;
                    let r = egui::Frame::new().fill(if on { theme::tint(&ctx) } else { theme::card(&ctx) }).stroke(egui::Stroke::new(if on { 2.0 } else { 1.0 }, if on { theme::accent(&ctx) } else { theme::line(&ctx) })).corner_radius(12).inner_margin(egui::Margin::same(12)).show(ui, |ui| ui.vertical(|ui| {
                        ui.set_width(w - 26.0);
                        ui.set_min_height(86.0);
                        ui.label(RichText::new(title).strong());
                        ui.label(kit::muted(ui, text).size(12.5));
                    }));
                    let resp = ui.interact(r.response.rect, ui.scope_id().with(("tpl", i)), egui::Sense::click());
                    if resp.clicked() {
                        self.shell.new.template = i;
                    }
                }
            });
            ui.add_space(8.0);
        }
        ui.separator();
        ui.horizontal(|ui| {
            let blank = self.shell.new.template == TEMPLATES.len() - 1;
            ui.add_enabled_ui(blank, |ui| {
                ui.label("Physics");
                let ids: Vec<String> = self.engine.registry().ids().iter().map(|s| s.to_string()).collect();
                egui::ComboBox::from_id_salt("new_phys").selected_text(self.physics_name(&self.shell.new.physics)).show_ui(ui, |ui| {
                    for p in ids {
                        let name = self.physics_name(&p);
                        ui.selectable_value(&mut self.shell.new.physics, p, name);
                    }
                });
                ui.label("Dimension");
                ui.selectable_value(&mut self.shell.new.three_d, false, "2D section");
                ui.selectable_value(&mut self.shell.new.three_d, true, "3D");
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if kit::primary(ui, "Create workflow").clicked() {
                    self.create_from_template();
                }
                if ui.button("Cancel").clicked() {
                    self.shell.sheet = Sheet::None;
                }
            });
        });
        ui.separator();
        ui.horizontal(|ui| {
            ui.label(kit::muted(ui, "Or"));
            if ui.link("Open a .nyx someone shared").clicked() {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
            }
            if ui.link("Import a .nyx.json project").clicked() {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
            }
            if ui.link("Write/run a .nys script...").clicked() {
                let title = if self.shell.new.title.trim().is_empty() { "Scripted workflow".to_string() } else { self.shell.new.title.trim().to_string() };
                self.new_draft(&title, Start::Blank { physics: "ems".into(), three_d: false });
                self.open_nys_script();
            }
            if ui.link(format!("Browse all {} examples", examples::list().len())).clicked() {
                self.shell.new.all_examples = !self.shell.new.all_examples;
            }
        });
        if self.shell.new.all_examples {
            egui::ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                for e in examples::list() {
                    if ui.link(&e.title).clicked() {
                        let title = if self.shell.new.title.trim().is_empty() { e.title.split('(').next().unwrap_or(&e.title).trim().to_string() } else { self.shell.new.title.clone() };
                        self.new_draft(&title, Start::Example(e.path.clone()));
                    }
                }
            });
        }
    }

    fn create_from_template(&mut self) {
        let n = &self.shell.new;
        let (title, start) = match TEMPLATES[n.template] {
            (t, _, Some(i)) => match examples::list().get(i) {
                Some(e) => (t.to_string(), Start::Example(e.path.clone())),
                None => {
                    self.set_error("The example for this card was not found: the `examples` folder is missing (set NYX_EXAMPLES).");
                    return;
                }
            },
            (t, _, None) => (t.to_string(), Start::Blank { physics: n.physics.clone(), three_d: n.three_d }),
        };
        let name = if n.title.trim().is_empty() { title } else { n.title.trim().to_string() };
        self.new_draft(&name, start);
        self.shell.new = Default::default();
    }

    // ------------------------------------------------------------------------------------------------------ palette

    fn sheet_palette(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let r = ui.add(egui::TextEdit::singleline(&mut self.shell.palette_q).font(egui::TextStyle::Heading).hint_text("Search workflows, runs and commands").desired_width(f32::INFINITY).margin(egui::Margin::symmetric(12, 10)));
        r.request_focus();
        ui.add_space(8.0);
        let q = nyx_engine::facade::store::search::Query::parse(&self.shell.palette_q);
        let entries = self.shell.entries.clone();
        let hits = nyx_engine::facade::store::search::search(&q, &entries);
        let mut open: Option<String> = None;
        let mut go: Option<crate::menu::Act> = None;
        if !q.is_empty() {
            kit::label_caps(ui, "Workflows and experiments", "");
            if hits.is_empty() {
                ui.label(kit::muted(ui, "No workflow matches. Hidden and done ones are searched too."));
            }
            for h in hits.iter().take(8) {
                let e = &entries[h.index];
                let rect = ui.horizontal(|ui| {
                    ui.vertical(|ui| {
                        ui.label(RichText::new(&e.wf.title).strong());
                        ui.label(kit::muted(ui, format!("matches {} · {}", h.reasons.join(" · "), crate::fmt::plural(e.wf.experiments.len(), "experiment"))).size(12.5));
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, format!("{} · {}", if e.wf.status == nyx_engine::facade::store::Status::Done { "done" } else { "active" }, ago(e.wf.touched_secs(), now_secs())))));
                });
                let resp = ui.interact(rect.response.rect, ui.scope_id().with(("ph", &e.wf.id)), egui::Sense::click());
                if resp.hovered() {
                    ui.painter().rect_filled(rect.response.rect, 6.0, theme::tint(&ctx).gamma_multiply(0.6));
                }
                if resp.clicked() {
                    open = Some(e.wf.id.clone());
                }
            }
            let runs = nyx_engine::facade::store::search::run_hits(&q, &entries, &hits, 5);
            if !runs.is_empty() {
                ui.add_space(6.0);
                kit::label_caps(ui, "Runs", "");
                for r in runs {
                    ui.horizontal(|ui| {
                        ui.label(kit::mono(r.run.split('_').next().unwrap_or(&r.run)));
                        ui.label(&r.what);
                        ui.label(kit::muted(ui, &entries[r.index].wf.title));
                    });
                }
            }
            ui.add_space(8.0);
            ui.label(kit::muted(ui, format!("Understood as: {}", q.understood().join(" · "))));
        }
        ui.add_space(8.0);
        let q_raw = self.shell.palette_q.clone();
        let matching_actions: Vec<&crate::menu::Item> = crate::menu::ACTIONS
            .iter()
            .filter(|it| it.act != crate::menu::Act::Palette && action_matches_query(it, &q_raw))
            .collect();
        if !q_raw.trim().is_empty() && ctx.input(|i| i.key_pressed(egui::Key::Enter)) {
            if let Some(first) = matching_actions.first() {
                go = Some(first.act);
            } else if let Some(h) = hits.first() {
                open = Some(entries[h.index].wf.id.clone());
            }
        }
        ui.columns(2, |cols| {
            kit::label_caps(&mut cols[0], "Do something", "");
            if matching_actions.is_empty() {
                cols[0].label(kit::muted(&cols[0], "No matching command."));
            }
            for it in &matching_actions {
                let label = if it.keys.is_empty() { it.label.to_string() } else { format!("{}  ({})", it.label, it.keys) };
                if cols[0].link(label).clicked() {
                    go = Some(it.act);
                }
            }
            kit::label_caps(&mut cols[1], "Browse what nyx can do", "");
            let ui = &mut cols[1];
            let browse: [(&str, &[(&str, crate::menu::Act)]); 3] = [
                (
                    "CONCEPT",
                    &[
                        ("Use-case templates", crate::menu::Act::New),
                        ("Geometry & parameters", crate::menu::Act::ConceptGeometry),
                        ("Materials & regions", crate::menu::Act::ConceptMaterials),
                    ],
                ),
                (
                    "BUILD",
                    &[
                        ("Mesh preview & rebuild", crate::menu::Act::BuildMesh),
                        ("Boundaries & conditions", crate::menu::Act::BuildConditions),
                        ("Study, wavelength & solver", crate::menu::Act::BuildStudy),
                        ("Parameter sweep", crate::menu::Act::BuildSweep),
                        ("Mesh convergence check", crate::menu::Act::CheckConvergence),
                    ],
                ),
                (
                    "ANALYSIS & SCRIPTS",
                    &[
                        ("Runs table, plot & compare", crate::menu::Act::AnalysisRuns),
                        ("Script (.nys) editor & runner", crate::menu::Act::ScriptEditor),
                        ("Export as .nys script", crate::menu::Act::ExportScript),
                    ],
                ),
            ];
            for (area, items) in browse {
                ui.add_space(4.0);
                ui.label(RichText::new(area).monospace().strong());
                ui.horizontal_wrapped(|ui| {
                    for (idx, (label, act)) in items.iter().enumerate() {
                        if idx > 0 {
                            ui.label(kit::muted(ui, "·"));
                        }
                        if ui.link(*label).clicked() {
                            go = Some(*act);
                        }
                    }
                });
            }
        });
        if let Some(id) = open {
            self.shell.sheet = Sheet::None;
            self.open_tab(&id, true);
        }
        if let Some(act) = go {
            // an entry that opens another sheet replaces the palette; the others close it first
            if !matches!(act, crate::menu::Act::New | crate::menu::Act::About | crate::menu::Act::Versions | crate::menu::Act::Share) {
                self.shell.sheet = Sheet::None;
            }
            self.do_act(&ctx, act);
        }
    }

    // ------------------------------------------------------------------------------------------------------ versions

    fn sheet_versions(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let Page::Workflow(id) = self.shell.page.clone() else {
            self.shell.sheet = Sheet::None;
            return;
        };
        let Some(w) = self.shell.open.get(&id) else { return };
        let Some(d) = &w.dir else { return };
        let hist = d.history();
        let title = d.wf.title.clone();
        let runs_all: Vec<nyx_engine::facade::store::runs::RunRecord> = d.wf.experiments.iter().flat_map(|e| d.runs(&e.id).0).collect();
        let root = d.root.clone();
        let sel = self.shell.version_sel.unwrap_or_else(|| hist.last().map(|v| v.n).unwrap_or(0));
        let mut restore = None;
        let mut fork = None;
        let mut name: Option<(u32, String)> = None;
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.label(kit::title("Versions", 36.0).color(theme::accent(&ctx)));
                ui.label(kit::muted(ui, format!("{title} · every change to the setup is kept. Nothing here is ever lost.")));
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                if kit::primary(ui, "Name this version").clicked() && !self.shell.version_name.trim().is_empty() {
                    name = Some((sel, self.shell.version_name.clone()));
                }
                ui.add(egui::TextEdit::singleline(&mut self.shell.version_name).hint_text("Name for the selected version").desired_width(220.0));
                if kit::ghost(ui, "Open folder").clicked() {
                    self.reveal(&root);
                }
            });
        });
        ui.add_space(8.0);
        ui.columns(2, |cols| {
            let ui = &mut cols[0];
            egui::ScrollArea::vertical().max_height(420.0).show(ui, |ui| {
                for v in hist.iter().rev() {
                    let on = v.n == sel;
                    let nruns = runs_all.iter().filter(|r| r.version == v.n).count();
                    let r = egui::Frame::new().fill(if on { theme::tint(&ctx) } else { theme::card(&ctx) }).stroke(egui::Stroke::new(if on { 1.5 } else { 1.0 }, if on { theme::accent(&ctx) } else { theme::line(&ctx) })).corner_radius(12).inner_margin(egui::Margin::same(12)).show(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.horizontal(|ui| {
                            let head = match (&v.name, v.n == hist.last().map(|x| x.n).unwrap_or(0)) {
                                (Some(n), _) => format!("v{} · {n}", v.n),
                                (None, true) => format!("v{} · now", v.n),
                                _ => format!("v{}", v.n),
                            };
                            ui.label(RichText::new(head).strong().size(16.0));
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, format!("{} · {}", if v.kind == VersionKind::Named { "named" } else { "auto" }, ago(nyx_engine::facade::store::util::parse_iso(&v.t).unwrap_or(0), now_secs())))));
                        });
                        ui.label(v.summary.as_str());
                        ui.label(kit::muted(ui, format!("{} recorded on this version", crate::fmt::plural(nruns, "run"))).size(12.5));
                    });
                    if ui.interact(r.response.rect, ui.scope_id().with(("ver", v.n)), egui::Sense::click()).clicked() {
                        self.shell.version_sel = Some(v.n);
                    }
                    ui.add_space(4.0);
                }
            });
            let ui = &mut cols[1];
            if let Some(v) = hist.iter().find(|v| v.n == sel) {
                let prev = v.n.saturating_sub(1);
                kit::label_caps(ui, &format!("What changed in v{}", v.n), &if prev > 0 { format!("from the version before it in this experiment ({})", v.experiment) } else { String::new() });
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                        if v.changes.is_empty() {
                            ui.label(kit::mono(if v.prev_hash.is_none() { "Created.\nThe first version of this experiment." } else { "No change in the setup." }));
                        } else {
                            ui.label(kit::mono(versions::render(&v.changes)).size(12.5));
                        }
                    });
                    ui.add_space(4.0);
                    ui.label(kit::muted(ui, format!("Runs on this version keep the numbers they were solved with ({}).", crate::fmt::plural(runs_all.iter().filter(|r| r.version == v.n).count(), "run"))).size(12.5));
                });
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    if ui.button(format!("Restore this version (adds v{})", hist.last().map(|x| x.n + 1).unwrap_or(1))).clicked() {
                        restore = Some(v.n);
                    }
                    if ui.button("Fork from here").clicked() {
                        fork = Some(v.clone());
                    }
                });
                ui.add_space(6.0);
                kit::plain_card(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.label(RichText::new("Works with git too").strong());
                    ui.label(kit::muted(ui, "The workflow is a plain folder. A generated .gitignore leaves out meshes and heavy results, so git init in the folder just works. nyx never needs git to run.").size(12.5));
                });
            }
        });
        if let Some((n, text)) = name {
            let r = self.shell.open.get_mut(&id).and_then(|w| w.dir.as_mut()).map(|d| d.name_version(n, &text));
            if let Some(Err(e)) = r {
                self.set_error(e);
            }
            self.shell.version_name.clear();
        }
        if let Some(n) = restore {
            self.flush();
            let r = self.shell.open.get_mut(&id).and_then(|w| w.dir.as_mut()).map(|d| d.restore(n));
            match r {
                Some(Ok((eid, _))) => {
                    self.shell.sheet = Sheet::None;
                    self.shell.loaded = None;
                    self.load_experiment(&id, &eid);
                    self.status = format!("Restored v{n} as a new version");
                }
                Some(Err(e)) => self.set_error(e),
                None => {}
            }
        }
        if let Some(v) = fork {
            let r = self.shell.open.get_mut(&id).and_then(|w| {
                let d = w.dir.as_mut()?;
                let spec = d.read_frozen(&v.spec_hash).ok()?;
                Some(d.add_experiment(&format!("Fork of v{}", v.n), &spec, Some(&v.experiment)))
            });
            match r {
                Some(Ok(eid)) => {
                    self.shell.sheet = Sheet::None;
                    self.refresh_entries();
                    self.load_experiment(&id, &eid);
                }
                Some(Err(e)) => self.set_error(e),
                None => self.set_error("Could not read that version's setup"),
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------ share

    fn sheet_share(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let Page::Workflow(id) = self.shell.page.clone() else {
            self.shell.sheet = Sheet::None;
            return;
        };
        self.flush();
        let Some(w) = self.shell.open.get(&id) else { return };
        let Some(d) = &w.dir else {
            ui.label(kit::title("Share this workflow", 34.0).color(theme::accent(&ctx)));
            ui.label(kit::muted(ui, "Nothing is saved yet. Change something first, then share."));
            return;
        };
        let only = vec![w.exp.clone()];
        let only_opt: Option<&[String]> = if self.shell.share_only_current { Some(&only) } else { None };
        let sizes = nyx_engine::facade::store::archive::sizes(d, only_opt).unwrap_or([(ShareLevel::Recipe, 0), (ShareLevel::Results, 0), (ShareLevel::Meshes, 0)]);
        let n_exp = d.wf.experiments.len();
        let exp_title = d.wf.experiment(&w.exp).map(|e| e.title.clone()).unwrap_or_default();
        let script_names = nyx_engine::facade::store::archive::scripts_of(d);
        let scripts = script_names.iter().filter(|s| !self.shell.share_exclude.contains(s)).count();
        let paths = nyx_engine::facade::store::archive::count_personal_paths(d);
        let file_name = format!("{}.nyx", nyx_engine::facade::store::util::slug(&d.wf.title));
        let strip = self.shell.settings().strip_paths;
        let exclude = self.shell.share_exclude.clone();
        let mut save = None;
        ui.label(kit::title("Share this workflow", 34.0).color(theme::accent(&ctx)));
        ui.label(kit::muted(ui, "One small file. The other person opens it in nyx and it rebuilds the setup on their machine."));
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            if kit::pill(ui, !self.shell.share_only_current, format!("Whole workflow · {n_exp} experiments")).clicked() {
                self.shell.share_only_current = false;
            }
            if kit::pill(ui, self.shell.share_only_current, format!("Only {exp_title}")).clicked() {
                self.shell.share_only_current = true;
            }
        });
        ui.add_space(8.0);
        for (lvl, size) in sizes {
            let on = self.shell.share_level == lvl;
            let r = egui::Frame::new().fill(if on { theme::tint(&ctx) } else { theme::card(&ctx) }).stroke(egui::Stroke::new(if on { 1.5 } else { 1.0 }, if on { theme::accent(&ctx) } else { theme::line(&ctx) })).corner_radius(12).inner_margin(egui::Margin::same(12)).show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.horizontal(|ui| {
                    ui.vertical(|ui| {
                        ui.label(RichText::new(lvl.name()).strong().size(16.0));
                        ui.label(kit::muted(ui, match lvl {
                            ShareLevel::Recipe => "Geometry, materials, mesh and solver settings, run table, versions, notes. Recommended.",
                            ShareLevel::Results => "Adds the saved result files, so they can look without re-solving.",
                            ShareLevel::Meshes => "Adds mesh files. nyx does not store meshes yet, so this equals the level before.",
                        }).size(12.5));
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::mono(human(size))));
                });
            });
            if ui.interact(r.response.rect, ui.scope_id().with(("lvl", lvl.name())), egui::Sense::click()).clicked() {
                self.shell.share_level = lvl;
            }
            ui.add_space(4.0);
        }
        ui.columns(2, |cols| {
            kit::plain_card(&mut cols[0], |ui| {
                ui.set_width(ui.available_width());
                ui.label(kit::mono("INCLUDED").strong());
                for l in ["geometry and materials", "mesh and solver settings", "run table, versions", "notes and tags"] {
                    ui.label(kit::mono(format!("+ {l}")));
                }
                ui.label(kit::mono("- meshes (rebuilt on the other side)"));
                ui.label(kit::mono(if self.shell.share_level == ShareLevel::Recipe { "- heavy result files" } else { "+ heavy result files" }));
            });
            kit::plain_card(&mut cols[1], |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("CHECK BEFORE SENDING").size(12.5).strong());
                if script_names.is_empty() {
                    ui.label("No custom scripts are included.");
                } else {
                    ui.label("Custom scripts. They stay switched off until the recipient approves them; untick any you do not want to send.");
                    for s in &script_names {
                        let mut on = !self.shell.share_exclude.contains(s);
                        if ui.checkbox(&mut on, s.as_str()).changed() {
                            if on { self.shell.share_exclude.retain(|x| x != s); } else { self.shell.share_exclude.push(s.clone()); }
                        }
                    }
                    ui.label(kit::muted(ui, format!("{scripts} of {} will be sent.", script_names.len())).size(12.5));
                }
                ui.label(kit::muted(ui, if strip { format!("Personal paths are reduced to file names ({paths} found).") } else { format!("Personal paths are kept as they are ({paths} found); turn on stripping in Settings.") }));
            });
        });
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.label(kit::mono(&file_name));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if kit::primary(ui, "Save .nyx").clicked() {
                    save = Some(("nyx", self.shell.share_level, strip));
                }
                if kit::ghost(ui, "Zip whole folder").on_hover_text("Everything including result files, as a plain .zip").clicked() {
                    save = Some(("zip", ShareLevel::Meshes, false));
                }
            });
        });
        if !self.shell.share_msg.is_empty() {
            ui.label(RichText::new(&self.shell.share_msg).color(theme::success(ui.ctx())));
        }
        if let Some((ext, level, strip)) = save {
            let Some(path) = rfd::FileDialog::new().add_filter(ext, &[ext]).set_file_name(if ext == "zip" { file_name.replace(".nyx", ".zip") } else { file_name.clone() }).save_file() else { return };
            let d = self.shell.open.get(&id).and_then(|w| w.dir.as_ref());
            let r = d.map(|d| nyx_engine::facade::store::archive::export_file(d, &ExportOptions { level, only: only_opt, strip_paths: strip, exclude: &exclude }, &path));
            match r {
                Some(Ok(n)) => self.shell.share_msg = format!("Saved {} ({})", path.display(), human(n)),
                Some(Err(e)) => self.set_error(e),
                None => {}
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------ open a received file

    fn sheet_open(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let Some((path, ins)) = self.shell.open_file.clone() else {
            self.shell.sheet = Sheet::None;
            return;
        };
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("file.nyx").to_string();
        ui.label(kit::title("Open a shared workflow", 34.0).color(theme::accent(&ctx)));
        ui.label(kit::mono(format!("{name} · {} · {}", human(ins.file_bytes), ins.level.map(|l| l.name()).unwrap_or("unknown level"))));
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            if ins.intact {
                kit::badge(ui, theme::raised(&ctx), theme::success(&ctx), "File is intact");
            } else {
                kit::badge(ui, theme::raised(&ctx), theme::danger(&ctx), "File cannot be opened");
            }
            if !ins.nyx_version.is_empty() {
                kit::badge(ui, theme::raised(&ctx), theme::success(&ctx), format!("Made with nyx {}", ins.nyx_version));
            }
            if !ins.scripts.is_empty() {
                kit::badge(ui, theme::warm(&ctx), theme::activity_text(&ctx), format!("{}, off until you approve", crate::fmt::plural(ins.scripts.len(), "script")));
            }
        });
        ui.add_space(8.0);
        if !ins.problems.is_empty() {
            kit::card(ui, theme::raised(&ctx), theme::danger(&ctx), 12, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("This file was not opened").strong().color(theme::danger(ui.ctx())));
                for p in &ins.problems {
                    ui.label(p);
                }
                ui.label(kit::muted(ui, "Nothing was written to your disk.").size(12.5));
            });
            ui.add_space(8.0);
            if ui.button("Close").clicked() {
                self.shell.sheet = Sheet::None;
            }
            return;
        }
        ui.columns(2, |cols| {
            kit::plain_card(&mut cols[0], |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("WHAT IS INSIDE").size(12.5).strong());
                ui.label(RichText::new(&ins.title).strong().size(17.0));
                ui.label(&ins.summary);
                ui.label(kit::mono(format!("{} · {} recorded", crate::fmt::plural(ins.experiments, "experiment"), crate::fmt::plural(ins.runs, "run"))));
                ui.label(kit::mono(format!("{} versions", ins.versions)));
                ui.label(kit::mono(format!("tags  {}", ins.tags.join("  "))));
            });
            kit::plain_card(&mut cols[1], |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("NOT INCLUDED, REBUILT HERE").size(12.5).strong());
                ui.label("Meshes and fields. nyx meshes and solves them again on this computer from the recipe.");
                ui.label(kit::muted(ui, format!("Unpacked size {}.", human(ins.unpacked_bytes))).size(12.5));
            });
        });
        ui.add_space(8.0);
        if !ins.scripts.is_empty() {
            ui.label(kit::muted(ui, format!("Scripts: {}. They stay off; nothing runs until you approve it in the workflow.", ins.scripts.join(", "))).size(12.5));
        }
        let dest = self.shell.ws.as_ref().map(|w| w.root.display().to_string()).unwrap_or_default();
        ui.horizontal(|ui| {
            ui.label(kit::muted(ui, "Save in"));
            ui.label(kit::mono(dest));
        });
        ui.add_space(8.0);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if kit::primary(ui, "Open as new workflow").clicked() {
                let root = self.shell.ws.as_ref().map(|w| w.root.clone());
                if let Some(root) = root {
                    match nyx_engine::facade::store::archive::import(&path, &root, &Limits::default()) {
                        Ok(d) => {
                            let id = d.wf.id.clone();
                            self.shell.sheet = Sheet::None;
                            self.refresh_entries();
                            self.open_tab(&id, true);
                            self.status = format!("Opened {} as a new workflow", d.wf.title);
                        }
                        Err(e) => self.set_error(e),
                    }
                }
            }
            if ui.button("Just look, do not save").clicked() {
                self.shell.sheet = Sheet::None;
            }
        });
    }

    // ------------------------------------------------------------------------------------------------------ about

    /// The version, what the fonts are and under which licence: the OFL asks that the licence text travels with the font and is shown.
    fn sheet_about(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.label(kit::title("nyx", 42.0).color(theme::accent(&ctx)));
        ui.label(RichText::new("creative chaos.").font(theme::title(20.0)));
        ui.label(kit::muted(ui, format!("Version {}. Finite-element photonics and phononics.", env!("CARGO_PKG_VERSION"))));
        ui.add_space(10.0);
        kit::label_caps(ui, "Licence", "");
        ui.label("nyx is dual licensed: you may use it under the GNU Affero General Public License, version 3, or under a commercial licence, whichever you choose.");
        ui.label("Copyright (C) 2026 p-j-r. nyx comes with ABSOLUTELY NO WARRANTY, to the extent permitted by law. You may convey it under the terms of the AGPL; the full text is the file LICENSE that comes with nyx and https://www.gnu.org/licenses/agpl-3.0.html. NOTICE lists the parts that keep their own licences (the geometry kernel and the ported solvers and meshers: BSD-2-Clause).");
        ui.add_space(10.0);
        kit::label_caps(ui, "Fonts and their licences", "");
        ui.label("Titles: EB Garamond Italic, Copyright 2017 The EB Garamond Project Authors, SIL Open Font License 1.1 (full text below). It is shipped as crates/app/assets/fonts/EBGaramond-Italic[wght].ttf with its licence next to it; when it is not found the window uses Georgia Italic.");
        ui.label("Text and numbers: Inter and JetBrains Mono when installed (not shipped with nyx), otherwise Segoe UI and Cascadia or Consolas.");
        ui.label("Fallback fonts inside egui (Ubuntu-Light, Hack, NotoEmoji, emoji-icon-font) keep their own licences; see the egui project.");
        ui.add_space(10.0);
        kit::label_caps(ui, "SIL Open Font License 1.1 (EB Garamond)", "");
        ui.label(kit::mono(include_str!("../assets/fonts/OFL-EBGaramond.txt")).size(12.0));
    }

    // ------------------------------------------------------------------------------------------------------ settings

    pub fn settings_ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let Some(ws) = &self.shell.ws else { return };
        let mut s = ws.settings.clone();
        let root = ws.root.clone();
        let before = s.clone();
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(40, 24))).show(ui, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.label(kit::title("Settings", 42.0).color(theme::accent(&ctx)));
                ui.label(kit::muted(ui, "Only what you would actually change. Everything else picks a good default."));
                ui.add_space(10.0);
                ui.columns(2, |cols| {
                    let ui = &mut cols[0];
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("WHERE THINGS LIVE").size(12.5).strong());
                        ui.label(kit::mono(root.display().to_string()));
                        ui.label(kit::muted(ui, "One folder per workflow. Move it and nyx follows (set NYX_WORKSPACE). Back it up like any folder.").size(12.5));
                        if ui.link("Open in file manager").clicked() {
                            self.reveal(&root);
                        }
                    });
                    ui.add_space(10.0);
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("COMPUTING").size(12.5).strong());
                        ui.horizontal(|ui| {
                            ui.label("Run at once");
                            ui.label(kit::mono("1 experiment"));
                        });
                        ui.label(kit::muted(ui, "One solve runs at a time: the solvers already use every core.").size(12.5));
                    });
                    ui.add_space(10.0);
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("SHARING").size(12.5).strong());
                        ui.horizontal(|ui| {
                            ui.label("Default level");
                            egui::ComboBox::from_id_salt("share_default").selected_text(s.share_default.name()).show_ui(ui, |ui| {
                                for l in ShareLevel::ALL {
                                    ui.selectable_value(&mut s.share_default, l, l.name());
                                }
                            });
                        });
                        ui.checkbox(&mut s.strip_paths, "Strip personal paths");
                        ui.horizontal(|ui| {
                            ui.label("Run shared scripts");
                            ui.label(kit::mono("ask first"));
                        });
                    });
                    let ui = &mut cols[1];
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("DASHBOARD TIDY-UP").size(12.5).strong());
                        ui.horizontal(|ui| {
                            ui.label("Hide after untouched");
                            ui.add(egui::DragValue::new(&mut s.hide_after_days).range(1..=365).suffix(" days"));
                        });
                        ui.horizontal(|ui| {
                            ui.label("Cards on the dashboard");
                            ui.add(egui::DragValue::new(&mut s.dashboard_cap).range(3..=60));
                        });
                        let mut keep = s.keep_heavy.is_some();
                        ui.horizontal(|ui| {
                            ui.label("Delete old heavy files");
                            if ui.checkbox(&mut keep, "keep the newest").changed() {
                                s.keep_heavy = if keep { Some(5) } else { None };
                            }
                            if let Some(n) = &mut s.keep_heavy {
                                ui.add(egui::DragValue::new(n).range(1..=1000).suffix(" per experiment"));
                            }
                        });
                        ui.label(kit::muted(ui, "Hidden workflows are still found by search. Nothing is deleted without asking; starred runs and every run record are always kept.").size(12.5));
                    });
                    ui.add_space(10.0);
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("UNITS AND LOOK").size(12.5).strong());
                        ui.horizontal(|ui| {
                            ui.label("Theme");
                            for (t, l) in [("system", "Follow system"), ("light", "Light"), ("dark", "Dark")] {
                                if ui.selectable_label(s.theme == t, l).clicked() {
                                    s.theme = t.into();
                                }
                            }
                        });
                        ui.horizontal(|ui| {
                            ui.label("Contrast");
                            for (t, l) in [("auto", "Follow system"), ("standard", "Standard"), ("high", "High")] {
                                if ui.selectable_label(s.contrast == t, l).clicked() {
                                    s.contrast = t.into();
                                }
                            }
                        });
                        ui.label(kit::muted(ui, "High contrast: stronger edges on every control. Follow system reads the Windows High Contrast setting; on other systems it stays standard.").size(12.5));
                        ui.checkbox(&mut s.reduce_motion, "Reduce motion (no sliding progress bar, no easing)");
                        ui.horizontal(|ui| {
                            ui.label("Canvas");
                            for (t, l) in [("dark", "Always dark"), ("theme", "Follow the look")] {
                                if ui.selectable_label(s.canvas == t, l).clicked() {
                                    s.canvas = t.into();
                                }
                            }
                        });
                        ui.label(kit::muted(ui, "The canvas is the background of the model and the plots. Dark is the default also in a light window; the region colours are chosen for the canvas that is drawn.").size(12.5));
                        ui.horizontal(|ui| {
                            ui.label("Colour regions by");
                            for (t, l) in [("material", "Material"), ("region", "Region")] {
                                if ui.selectable_label(s.region_colouring == t, l).clicked() {
                                    s.region_colouring = t.into();
                                }
                            }
                        });
                        ui.label(kit::muted(ui, "Material: regions of one material share a colour. Region: every region has its own. Eight colours are solid, the next eight are hatched, further regions are grey; the number chip and the legend always say which is which.").size(12.5));
                        ui.horizontal(|ui| {
                            ui.label("Length");
                            ui.label(kit::mono("nm"));
                            ui.label("Frequency");
                            ui.label(kit::mono("THz"));
                        });
                        ui.label(kit::muted(ui, "Unit choice is fixed for now; the numbers are always stored in SI.").size(12.5));
                        if ui.checkbox(&mut self.shell.classic, "Use the classic layout (menus and six tabs)").changed() {
                            self.status = "Layout changed.".into();
                        }
                        if ui.link("About nyx and licences").clicked() {
                            self.shell.sheet = Sheet::About;
                        }
                    });
                    ui.add_space(10.0);
                    kit::plain_card(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(RichText::new("ASSISTANT").size(12.5).strong());
                        ui.checkbox(&mut s.assistant, "Let an AI assistant work in nyx");
                        ui.add_enabled_ui(s.assistant, |ui| {
                            egui::ComboBox::from_id_salt("assist_level").selected_text(s.assistant_level.clone()).show_ui(ui, |ui| {
                                for l in ["Read results only", "Read results and queue runs", "Ask before changing a setup"] {
                                    ui.selectable_value(&mut s.assistant_level, l.to_string(), l);
                                }
                            });
                        });
                        ui.label(kit::muted(ui, "Off by default. This records what you allow; no assistant is connected in this version, so nothing can act yet.").size(12.5));
                    });
                });
            });
        });
        if s != before {
            ctx.request_repaint(); // `sync_look` applies the look at the start of the next frame
            if let Some(ws) = &mut self.shell.ws {
                ws.settings = s;
                if let Err(e) = ws.save_settings() {
                    self.set_error(e);
                }
            }
            self.refresh_entries();
        }
    }
}

fn action_matches_query(it: &crate::menu::Item, q_raw: &str) -> bool {
    let q = q_raw.trim().to_lowercase();
    if q.is_empty() {
        return true;
    }
    let hay = format!("{} {} {}", it.menu, it.label, it.keys).to_lowercase();
    let extra = match it.act {
        crate::menu::Act::BuildStudy => "frequency wavelength modes eigenfrequency solver order p1 p2 target neff study",
        crate::menu::Act::BuildConditions => "boundary boundaries condition conditions pml pec port source load constraint physics",
        crate::menu::Act::BuildSweep => "sweep parameter range points batch",
        crate::menu::Act::ConceptGeometry => "geometry dimensions parameters width height thickness rect circle polygon cad step",
        crate::menu::Act::ConceptMaterials => "material materials silicon silica sapphire refractive index epsilon conductivity",
        crate::menu::Act::AnalysisRuns => "analysis runs plot compare table filter export csv field modes polarization poincare",
        crate::menu::Act::ScriptEditor => "script .nys python monty code automate batch editor runfile",
        crate::menu::Act::ExportScript => "export script .nys save code",
        crate::menu::Act::BuildMesh => "mesh remesh elements cells triangulation",
        crate::menu::Act::CheckConvergence => "convergence refine accuracy check error",
        crate::menu::Act::Compute => "compute solve run f5 start simulation",
        _ => "",
    };
    q.split_whitespace().all(|tok| hay.contains(tok) || extra.contains(tok))
}

fn human(b: u64) -> String {
    let b = b as f64;
    if b < 1024.0 {
        format!("{b:.0} B")
    } else if b < 1024.0 * 1024.0 {
        format!("{:.0} KB", b / 1024.0)
    } else if b < 1024.0 * 1024.0 * 1024.0 {
        format!("{:.1} MB", b / 1024.0 / 1024.0)
    } else {
        format!("{:.2} GB", b / 1024.0 / 1024.0 / 1024.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn palette_query_filters_commands_and_synonyms() {
        let find = |q: &str| -> Vec<crate::menu::Act> {
            crate::menu::ACTIONS
                .iter()
                .filter(|it| it.act != crate::menu::Act::Palette && action_matches_query(it, q))
                .map(|it| it.act)
                .collect()
        };
        assert!(find("wavelength").contains(&crate::menu::Act::BuildStudy));
        assert!(find("freq").contains(&crate::menu::Act::BuildStudy));
        assert!(find(".nys").contains(&crate::menu::Act::ScriptEditor));
        assert!(find("material").contains(&crate::menu::Act::ConceptMaterials));
        assert!(find("geom").contains(&crate::menu::Act::ConceptGeometry));
        assert!(find("mesh").contains(&crate::menu::Act::BuildMesh));
    }

    #[test]
    fn nys_script_sheet_prefills_from_project() {
        let ctx = egui::Context::default();
        let mut a = App::with_context(ctx.clone());
        assert!(!a.nys.open);
        a.open_nys_script();
        assert!(a.nys.open);
        assert!(!a.nys.code.trim().is_empty(), "prefilled with project_script");
        let mut out = ctx.run_ui(Default::default(), |ui| a.sheets_ui(ui.ctx()));
        out.textures_delta.clear();
    }
}
