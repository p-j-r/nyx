//! The sheets that open over the shell (new workflow, search palette, versions, share, open a received file) and the settings page.

use crate::app::App;
use crate::examples;
use crate::kit;
use crate::shell::{Page, Sheet};
use crate::theme;
use eframe::egui::{self, Color32, RichText};
use femcore::problem::Physics;
use nyx_store::archive::{ExportOptions, Limits};
use nyx_store::util::{ago, now_secs};
use nyx_store::versions::{self, VersionKind};
use nyx_store::ShareLevel;

/// Start-from cards: (title, text, example index in `examples::NAMES`, or `None` when no template exists yet)
const TEMPLATES: [(&str, &str, Option<usize>); 8] = [
    ("Optical waveguide", "Strip or rib cross-section, modes and n_eff", Some(14)),
    ("Strip on sapphire", "Optical and acoustic modes of a silicon strip", Some(12)),
    ("Bragg grating", "Periodic cell, 3D, band-edge modes", Some(16)),
    ("Ring resonator", "Whispering-gallery modes, axisymmetric", None),
    ("Vibrating structure", "Solid mechanics, eigenfrequencies", Some(1)),
    ("Heat transfer", "Stationary or transient temperature", Some(2)),
    ("Electrostatics", "Potentials, fields, capacitance", Some(0)),
    ("Blank experiment", "Choose physics, dimension and geometry yourself", None),
];

impl App {
    pub fn sheets_ui(&mut self, ctx: &egui::Context) {
        let close = match self.shell.sheet {
            Sheet::None => return,
            Sheet::New => kit::modal(ctx, "new", 940.0, |ui| self.sheet_new(ui)),
            Sheet::Palette => kit::modal(ctx, "palette", 860.0, |ui| self.sheet_palette(ui)),
            Sheet::Versions => kit::modal(ctx, "versions", 1000.0, |ui| self.sheet_versions(ui)),
            Sheet::Share => kit::modal(ctx, "share", 820.0, |ui| self.sheet_share(ui)),
            Sheet::Open => kit::modal(ctx, "open", 820.0, |ui| self.sheet_open(ui)),
        };
        if close {
            self.shell.sheet = Sheet::None;
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
        for row in 0..2 {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 12.0;
                for col in 0..4 {
                    let i = row * 4 + col;
                    let (title, text, ex) = TEMPLATES[i];
                    let on = sel == i;
                    let avail = ex.is_some() || i == 7;
                    let r = egui::Frame::new().fill(if on { theme::tint(&ctx) } else { theme::card(&ctx) }).stroke(egui::Stroke::new(if on { 2.0 } else { 1.0 }, if on { theme::accent(&ctx) } else { theme::line(&ctx) })).corner_radius(12).inner_margin(egui::Margin::same(12)).show(ui, |ui| ui.vertical(|ui| {
                        ui.set_width(w - 26.0);
                        ui.set_min_height(86.0);
                        ui.label(RichText::new(title).strong().color(if avail { theme::text(&ctx) } else { theme::muted(&ctx) }));
                        ui.label(kit::muted(ui, if avail { text.to_string() } else { format!("{text} (template coming soon)") }).size(12.5));
                    }));
                    let resp = ui.interact(r.response.rect, ui.id().with(("tpl", i)), egui::Sense::click());
                    if resp.clicked() && avail {
                        self.shell.new.template = i;
                    }
                }
            });
            ui.add_space(8.0);
        }
        ui.separator();
        ui.horizontal(|ui| {
            let blank = self.shell.new.template == 7;
            ui.add_enabled_ui(blank, |ui| {
                ui.label("Physics");
                egui::ComboBox::from_id_salt("new_phys").selected_text(self.shell.new.physics.name()).show_ui(ui, |ui| {
                    for p in Physics::ALL {
                        ui.selectable_value(&mut self.shell.new.physics, p, p.name());
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
            if ui.link("Import a .nyx.xml project").clicked() {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
            }
            if ui.link(format!("Browse all {} examples", examples::NAMES.len())).clicked() {
                self.shell.new.all_examples = !self.shell.new.all_examples;
            }
        });
        if self.shell.new.all_examples {
            egui::ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                for (i, n) in examples::NAMES.iter().enumerate() {
                    if ui.link(*n).clicked() {
                        let title = if self.shell.new.title.trim().is_empty() { n.split('(').next().unwrap_or(n).trim().to_string() } else { self.shell.new.title.clone() };
                        self.new_draft(&title, examples::by_index(i));
                    }
                }
            });
        }
    }

    fn create_from_template(&mut self) {
        let n = &self.shell.new;
        let (title, project) = match TEMPLATES[n.template] {
            (t, _, Some(i)) => (t.to_string(), examples::by_index(i)),
            (t, _, None) => (t.to_string(), examples::blank(n.physics, n.three_d)),
        };
        let name = if n.title.trim().is_empty() { title } else { n.title.trim().to_string() };
        self.new_draft(&name, project);
        self.shell.new = Default::default();
    }

    // ------------------------------------------------------------------------------------------------------ palette

    fn sheet_palette(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let r = ui.add(egui::TextEdit::singleline(&mut self.shell.palette_q).font(egui::TextStyle::Heading).hint_text("Search workflows, runs and commands").desired_width(f32::INFINITY).margin(egui::Margin::symmetric(12, 10)));
        r.request_focus();
        ui.add_space(8.0);
        let q = nyx_store::search::Query::parse(&self.shell.palette_q);
        let entries = self.shell.entries.clone();
        let hits = nyx_store::search::search(&q, &entries);
        let mut open: Option<String> = None;
        let mut go: Option<&'static str> = None;
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
                        ui.label(kit::muted(ui, format!("matches {} · {} experiments", h.reasons.join(" · "), e.wf.experiments.len())).size(12.5));
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, format!("{} · {}", if e.wf.status == nyx_store::Status::Done { "done" } else { "active" }, ago(e.wf.touched_secs(), now_secs())))));
                });
                let resp = ui.interact(rect.response.rect, ui.id().with(("ph", &e.wf.id)), egui::Sense::click());
                if resp.hovered() {
                    ui.painter().rect_filled(rect.response.rect, 6.0, theme::tint(&ctx).gamma_multiply(0.6));
                }
                if resp.clicked() {
                    open = Some(e.wf.id.clone());
                }
            }
            let runs = nyx_store::search::run_hits(&q, &entries, &hits, 5);
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
        ui.columns(2, |cols| {
            kit::label_caps(&mut cols[0], "Do something", "");
            for (label, act) in [("New workflow", "new"), ("Open a .nyx file", "open"), ("Go to the dashboard", "dash"), ("Compute the open experiment (F5)", "compute"), ("Settings", "settings")] {
                if cols[0].link(label).clicked() {
                    go = Some(act);
                }
            }
            kit::label_caps(&mut cols[1], "Browse what nyx can do", "");
            let ui = &mut cols[1];
            for (area, items) in [("CONCEPT", "Use-case templates · Named regions · Geometry as code (soon)"), ("BUILD", "Mesh preview · Boundary numbers · Solver diagnostics"), ("ANALYSIS", "Runs table and filter · Compare runs · Polarization · Sensitivity · Your script")] {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new(area).monospace().strong());
                    ui.label(kit::muted(ui, items).size(12.5));
                });
            }
            ui.label(kit::muted(ui, "Not yet: fast frequency sweep, optimise a parameter, adaptive meshing, surrogate models.").size(12.5));
        });
        if let Some(id) = open {
            self.shell.sheet = Sheet::None;
            self.open_tab(&id, true);
        }
        match go {
            Some("new") => self.shell.sheet = Sheet::New,
            Some("open") => {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
            }
            Some("dash") => {
                self.shell.sheet = Sheet::None;
                self.flush();
                self.shell.page = Page::Dashboard;
            }
            Some("compute") => {
                self.shell.sheet = Sheet::None;
                self.shell_compute();
            }
            Some("settings") => {
                self.shell.sheet = Sheet::None;
                self.shell.page = Page::Settings;
            }
            _ => {}
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
        let runs_all: Vec<nyx_store::runs::RunRecord> = d.wf.experiments.iter().flat_map(|e| d.runs(&e.id).0).collect();
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
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, format!("{} · {}", if v.kind == VersionKind::Named { "named" } else { "auto" }, ago(nyx_store::util::parse_iso(&v.t).unwrap_or(0), now_secs())))));
                        });
                        ui.label(v.summary.as_str());
                        ui.label(kit::muted(ui, format!("{nruns} run{} recorded on this version", if nruns == 1 { "" } else { "s" })).size(12.5));
                    });
                    if ui.interact(r.response.rect, ui.id().with(("ver", v.n)), egui::Sense::click()).clicked() {
                        self.shell.version_sel = Some(v.n);
                    }
                    ui.add_space(4.0);
                }
            });
            let ui = &mut cols[1];
            if let Some(v) = hist.iter().find(|v| v.n == sel) {
                let prev = v.n.checked_sub(1).unwrap_or(0);
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
                    ui.label(kit::muted(ui, format!("Runs on this version keep the numbers they were solved with ({} run{}).", runs_all.iter().filter(|r| r.version == v.n).count(), if runs_all.iter().filter(|r| r.version == v.n).count() == 1 { "" } else { "s" })).size(12.5));
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
                let xml = d.read_frozen(&v.spec_hash).ok()?;
                Some(d.add_experiment(&format!("Fork of v{}", v.n), &xml, Some(&v.experiment)))
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
        let sizes = nyx_store::archive::sizes(d, only_opt).unwrap_or([(ShareLevel::Recipe, 0), (ShareLevel::Results, 0), (ShareLevel::Meshes, 0)]);
        let n_exp = d.wf.experiments.len();
        let exp_title = d.wf.experiment(&w.exp).map(|e| e.title.clone()).unwrap_or_default();
        let script_names = nyx_store::archive::scripts_of(d);
        let scripts = script_names.iter().filter(|s| !self.shell.share_exclude.contains(s)).count();
        let paths = nyx_store::archive::count_personal_paths(d);
        let file_name = format!("{}.nyx", nyx_store::util::slug(&d.wf.title));
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
            if ui.interact(r.response.rect, ui.id().with(("lvl", lvl.name())), egui::Sense::click()).clicked() {
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
            ui.label(RichText::new(&self.shell.share_msg).color(theme::GREEN));
        }
        if let Some((ext, level, strip)) = save {
            let Some(path) = rfd::FileDialog::new().add_filter(ext, &[ext]).set_file_name(if ext == "zip" { file_name.replace(".nyx", ".zip") } else { file_name.clone() }).save_file() else { return };
            let d = self.shell.open.get(&id).and_then(|w| w.dir.as_ref());
            let r = d.map(|d| nyx_store::archive::export_file(d, &ExportOptions { level, only: only_opt, strip_paths: strip, exclude: &exclude }, &path));
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
                kit::badge(ui, Color32::from_rgb(0xE3, 0xF2, 0xE9), theme::GREEN, "File is intact");
            } else {
                kit::badge(ui, Color32::from_rgb(0xF8, 0xE1, 0xE1), theme::RED, "File cannot be opened");
            }
            if !ins.nyx_version.is_empty() {
                kit::badge(ui, Color32::from_rgb(0xE3, 0xF2, 0xE9), theme::GREEN, format!("Made with nyx {}", ins.nyx_version));
            }
            if !ins.scripts.is_empty() {
                kit::badge(ui, theme::AMBER_TINT, Color32::from_rgb(0x8A, 0x55, 0x00), format!("{} script{}, off until you approve", ins.scripts.len(), if ins.scripts.len() == 1 { "" } else { "s" }));
            }
        });
        ui.add_space(8.0);
        if !ins.problems.is_empty() {
            kit::card(ui, Color32::from_rgb(0xF8, 0xE1, 0xE1), theme::RED, 12, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("This file was not opened").strong().color(theme::RED));
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
                ui.label(kit::mono(format!("{} experiments · {} runs recorded", ins.experiments, ins.runs)));
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
        ui.add_enabled(false, egui::Checkbox::new(&mut false, "Rebuild now and check the answers match the recorded ones (coming soon)"));
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
                    match nyx_store::archive::import(&path, &root, &Limits::default()) {
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
                        ui.label(kit::muted(ui, "One solve runs at a time: the solvers already use every core. A queue, per-run threads and memory limits are coming soon.").size(12.5));
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
                            ui.label("Length");
                            ui.label(kit::mono("nm"));
                            ui.label("Frequency");
                            ui.label(kit::mono("THz"));
                        });
                        ui.label(kit::muted(ui, "Unit choice is fixed for now; the numbers are always stored in SI.").size(12.5));
                        if ui.checkbox(&mut self.shell.classic, "Use the classic layout (menus and six tabs)").changed() {
                            self.status = "Layout changed.".into();
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
            self.shell.dark = s.theme == "dark" || (s.theme == "system" && ctx.system_theme() == Some(egui::Theme::Dark));
            theme::apply_mode(&ctx, self.shell.dark);
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
