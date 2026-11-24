// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The dashboard: running workflows with live status, pinned and recent ones as cards, everything older in a list, and search.
//! The rules (three days, the card cap, "done", nothing deleted) are `nyx_engine::facade::store::Workspace::sections`.

use crate::app::App;
use crate::kit;
use crate::shell::Sheet;
use crate::theme;
use eframe::egui::{self, Color32, RichText, Sense, Vec2};
use nyx_engine::facade::store::util::{ago, now_secs};
use nyx_engine::facade::store::{Entry, Status, WorkflowDir};

pub fn physics_label(p: &str) -> &'static str {
    match p {
        "Electrostatics" => "Electrostatics",
        "Heat" => "Heat transfer",
        "Solid" => "Solid mechanics",
        "Emw" | "ems" => "EM waves",
        "dfs" => "Diffusion",
        "sms" => "Solid mechanics",
        _ => "Other",
    }
}

impl App {
    /// Changes a workflow's pin or status on disk (open tabs are updated too).
    pub(crate) fn edit_workflow(&mut self, id: &str, f: impl FnOnce(&mut WorkflowDir) -> Result<(), String>) {
        let res = if let Some(d) = self.shell.open.get_mut(id).and_then(|w| w.dir.as_mut()) {
            f(d)
        } else if let Some(e) = self.shell.entries.iter().find(|e| e.wf.id == id) {
            WorkflowDir::open_with(&e.dir, nyx_engine::store_link::card_of).and_then(|mut d| f(&mut d))
        } else {
            Ok(())
        };
        if let Err(e) = res {
            self.set_error(e);
        }
        self.refresh_entries();
    }

    pub fn reveal(&mut self, path: &std::path::Path) {
        #[cfg(windows)]
        let r = std::process::Command::new("explorer").arg(path).spawn();
        #[cfg(target_os = "macos")]
        let r = std::process::Command::new("open").arg(path).spawn();
        #[cfg(all(unix, not(target_os = "macos")))]
        let r = std::process::Command::new("xdg-open").arg(path).spawn();
        if let Err(e) = r {
            self.set_error(format!("Could not open {}: {e}", path.display()));
        }
    }

    fn entry_matches_filters(&self, e: &Entry) -> bool {
        let s = &self.shell;
        (s.dash_physics.is_empty() || e.wf.experiments.iter().any(|x| physics_label(&x.card.physics) == s.dash_physics))
            && (s.dash_status.is_empty() || format!("{:?}", e.wf.status).eq_ignore_ascii_case(&s.dash_status))
            && (s.dash_material.is_empty() || e.wf.experiments.iter().any(|x| x.card.keywords.iter().any(|k| k.eq_ignore_ascii_case(&s.dash_material))))
    }

    pub fn dashboard_ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        egui::CentralPanel::default().frame(egui::Frame::new().fill(theme::page(&ctx)).inner_margin(egui::Margin::symmetric(40, 24))).show(ui, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.vertical(|ui| {
                        ui.label(kit::title("Workflows", 46.0));
                        let days = self.shell.settings().hide_after_days;
                        ui.label(kit::muted(ui, format!("Everything you touched in the last {days} days. Older work is one search away.")));
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if kit::primary(ui, "+  New workflow   Ctrl N").clicked() {
                            self.shell.sheet = Sheet::New;
                        }
                        if kit::ghost(ui, "Open .nyx").clicked() {
                            self.pick_nyx_file();
                        }
                    });
                });
                ui.add_space(14.0);
                self.dashboard_search(ui);
                ui.add_space(16.0);
                if !self.shell.dash_q.trim().is_empty() {
                    self.dashboard_results(ui);
                } else {
                    self.dashboard_sections(ui);
                }
            });
        });
    }

    fn dashboard_search(&mut self, ui: &mut egui::Ui) {
        let entries = self.shell.entries.clone();
        let mut physics: Vec<&str> = entries.iter().flat_map(|e| e.wf.experiments.iter().map(|x| physics_label(&x.card.physics))).collect();
        physics.sort();
        physics.dedup();
        let mut materials: Vec<String> = entries.iter().flat_map(|e| e.wf.experiments.iter().flat_map(|x| x.card.keywords.iter().filter(|k| matches!(k.as_str(), "Si" | "SiO2" | "Al2O3" | "Si3N4" | "air")).cloned())).collect();
        materials.sort();
        materials.dedup();
        ui.horizontal(|ui| {
            let w = (ui.available_width() - 330.0).max(200.0);
            ui.add_sized([w, 34.0], egui::TextEdit::singleline(&mut self.shell.dash_q).hint_text("Find a workflow, e.g. \u{201c}plasmon scattering 300 nm sapphire frame\u{201d}").margin(egui::Margin::symmetric(12, 8)));
            filter_menu(ui, "Physics", &mut self.shell.dash_physics, &physics.iter().map(|s| s.to_string()).collect::<Vec<_>>());
            filter_menu(ui, "Material", &mut self.shell.dash_material, &materials);
            filter_menu(ui, "Status", &mut self.shell.dash_status, &["Active".to_string(), "Done".to_string()]);
        });
    }

    fn dashboard_results(&mut self, ui: &mut egui::Ui) {
        let q = nyx_engine::facade::store::search::Query::parse(&self.shell.dash_q);
        let entries = self.shell.entries.clone();
        let hits = nyx_engine::facade::store::search::search(&q, &entries);
        ui.label(kit::muted(ui, format!("Understood as: {}", q.understood().join(" · "))));
        ui.add_space(6.0);
        let mut open = None;
        let mut shown = 0;
        for h in &hits {
            let e = &entries[h.index];
            if !self.entry_matches_filters(e) {
                continue;
            }
            shown += 1;
            let r = kit::plain_card_r(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.horizontal(|ui| {
                    ui.label(RichText::new(&e.wf.title).strong().size(16.0));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, format!("{} · {}", if e.wf.status == Status::Done { "done" } else { "active" }, ago(e.wf.touched_secs(), now_secs())))));
                });
                ui.label(kit::muted(ui, format!("matches {}", h.reasons.join(" · "))));
            });
            if ui.interact(r.rect, ui.scope_id().with(("hit", &e.wf.id)), Sense::click()).clicked() {
                open = Some(e.wf.id.clone());
            }
            ui.add_space(4.0);
        }
        if shown == 0 {
            ui.label(kit::muted(ui, "Nothing matches. Hidden and done workflows are searched too; try fewer words."));
        }
        let runs = nyx_engine::facade::store::search::run_hits(&q, &entries, &hits, 5);
        if !runs.is_empty() {
            ui.add_space(8.0);
            kit::label_caps(ui, "Runs", "");
            for r in runs {
                ui.horizontal(|ui| {
                    ui.label(kit::mono(&r.run));
                    ui.label(r.what);
                    ui.label(kit::muted(ui, &entries[r.index].wf.title));
                });
            }
        }
        if let Some(id) = open {
            self.open_tab(&id, true);
        }
    }

    fn dashboard_sections(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let entries: Vec<Entry> = self.shell.entries.iter().filter(|e| self.entry_matches_filters(e)).cloned().collect();
        let running_ids: Vec<String> = self.shell.run.iter().map(|r| r.wf_id.clone()).collect();
        let Some(ws) = &self.shell.ws else { return };
        let sec = ws.sections(&entries, &running_ids, now_secs());
        let mut open: Option<String> = None;
        let mut menu: Option<(String, &'static str)> = None;

        if entries.is_empty() {
            kit::plain_card(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("No workflows yet").strong().size(17.0));
                ui.label(kit::muted(ui, "Start one from a template, or open a .nyx file somebody shared."));
                ui.add_space(6.0);
                if kit::primary(ui, "New workflow").clicked() {
                    self.shell.sheet = Sheet::New;
                }
            });
            return;
        }

        if !sec.running.is_empty() {
            ui.horizontal(|ui| {
                kit::glyph(ui, kit::State::Running);
                kit::label_caps(ui, "Running", "follow along without opening anything");
            });
            ui.add_space(6.0);
            for &i in &sec.running {
                let e = &entries[i];
                let mut cancel = false;
                kit::card(ui, theme::warm(&ctx), theme::activity(&ctx), 16, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(&e.wf.title).strong().size(18.0));
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if kit::danger(ui, "Cancel").clicked() {
                                cancel = true;
                            }
                        });
                    });
                    ui.label(kit::muted(ui, e.wf.summary()));
                    if let Some(j) = &self.job {
                        let st = j.1.status();
                        ui.add_space(6.0);
                        ui.label(RichText::new(if st.stage.is_empty() { "Starting".to_string() } else { st.stage.clone() }).strong());
                        kit::bar(ui, st.fraction.map(|f| f as f32), 360.0);
                        ui.label(kit::mono(format!("{}  elapsed {:.0}:{:02.0}", st.fraction.map(|f| format!("{:.0} %", 100.0 * f)).unwrap_or_else(|| "no percentage reported yet".into()), (st.total_seconds / 60.0).floor(), st.total_seconds % 60.0)));
                    }
                });
                if cancel {
                    self.shell.confirm_cancel = true;
                }
                ui.add_space(6.0);
            }
            ui.add_space(10.0);
        }

        ui.horizontal(|ui| {
            kit::label_caps(ui, "Pinned and recent", &format!("{} of {} shown", sec.recent.len() + sec.running.len(), self.shell.settings().dashboard_cap.min(entries.len())));
        });
        ui.add_space(6.0);
        if sec.recent.is_empty() {
            ui.label(kit::muted(ui, "Nothing was touched recently. Open one from Earlier below, or start a new workflow."));
        }
        let cols = ((ui.available_width() + 14.0) / 280.0).floor().clamp(1.0, 4.0) as usize;
        let w = (ui.available_width() - 14.0 * (cols as f32 - 1.0)) / cols as f32;
        for row in sec.recent.chunks(cols) {
            ui.horizontal_top(|ui| {
                ui.spacing_mut().item_spacing.x = 14.0;
                for &i in row {
                    let e = &entries[i];
                    let spark = self.shell.spark.get(&e.wf.id).cloned().unwrap_or_default();
                    let resp = workflow_card(ui, e, &spark, w);
                    if resp.clicked() {
                        open = Some(e.wf.id.clone());
                    }
                    resp.context_menu(|ui| {
                        if ui.button(if e.wf.pinned { "Unpin" } else { "Pin to the dashboard" }).clicked() {
                            menu = Some((e.wf.id.clone(), "pin"));
                            ui.close();
                        }
                        if ui.button("Mark done (moves to Earlier)").clicked() {
                            menu = Some((e.wf.id.clone(), "done"));
                            ui.close();
                        }
                        if ui.button("Open folder").clicked() {
                            menu = Some((e.wf.id.clone(), "folder"));
                            ui.close();
                        }
                    });
                }
            });
            ui.add_space(14.0);
        }

        ui.add_space(8.0);
        let head = ui.horizontal(|ui| {
            let arrow = if self.shell.earlier_open { "v" } else { ">" };
            if ui.add(egui::Label::new(RichText::new(arrow).strong()).sense(Sense::click())).clicked() {
                self.shell.earlier_open = !self.shell.earlier_open;
            }
            kit::label_caps(ui, "Earlier", &format!("{} · untouched for {} or marked done · never deleted", crate::fmt::plural(sec.earlier.len(), "workflow"), crate::fmt::plural(self.shell.settings().hide_after_days as usize, "day")));
        });
        if head.response.clicked() {
            self.shell.earlier_open = !self.shell.earlier_open;
        }
        if self.shell.earlier_open {
            ui.add_space(4.0);
            // a table, one row per workflow: name, setup, state (glyph and word), when it was last touched
            let aw = ui.available_width();
            let gap = 18.0;
            let (w_name, w_setup, w_state) = (aw * 0.26, aw * 0.38, 110.0);
            let w_time = (aw - w_name - w_setup - w_state - 3.0 * gap - 4.0).max(60.0);
            // add_sized centres a label in its box; a table wants its text on the left
            let cell = |ui: &mut egui::Ui, w: f32, label: egui::Label| ui.allocate_ui_with_layout(Vec2::new(w, 22.0), egui::Layout::left_to_right(egui::Align::Center), |ui| { ui.set_min_width(w); ui.add(label) }).inner;
            egui::Grid::new("earlier_table").num_columns(4).spacing([gap, 10.0]).show(ui, |ui| {
                for (h, w) in [("Workflow", w_name), ("Setup", w_setup), ("State", w_state), ("Last touched", w_time)] {
                    cell(ui, w, egui::Label::new(RichText::new(h).size(12.0).strong().color(theme::muted(&ctx))).truncate());
                }
                ui.end_row();
                for &i in &sec.earlier {
                    let e = &entries[i];
                    let done = e.wf.status == Status::Done;
                    let mut cells: Vec<egui::Response> = vec![];
                    cells.push(cell(ui, w_name, egui::Label::new(RichText::new(&e.wf.title).strong()).truncate().sense(Sense::click())));
                    cells.push(cell(ui, w_setup, egui::Label::new(kit::muted(ui, e.wf.summary())).truncate().sense(Sense::click())));
                    let state = ui.allocate_ui_with_layout(Vec2::new(w_state, 22.0), egui::Layout::left_to_right(egui::Align::Center), |ui| {
                        ui.set_min_width(w_state);
                        kit::glyph(ui, if done { kit::State::Done } else { kit::State::Idle });
                        ui.label(kit::muted(ui, if done { "Done" } else { "Active" }));
                    });
                    cells.push(ui.interact(state.response.rect, ui.scope_id().with(("earlier_state", &e.wf.id)), Sense::click()));
                    cells.push(cell(ui, w_time, egui::Label::new(kit::muted(ui, ago(e.wf.touched_secs(), now_secs()))).truncate().sense(Sense::click())));
                    ui.end_row();
                    for c in cells {
                        if c.clicked() {
                            open = Some(e.wf.id.clone());
                        }
                        c.context_menu(|ui| {
                            if ui.button("Pin to the dashboard").clicked() {
                                menu = Some((e.wf.id.clone(), "pin"));
                                ui.close();
                            }
                            if ui.button("Mark active again").clicked() {
                                menu = Some((e.wf.id.clone(), "active"));
                                ui.close();
                            }
                            if ui.button("Open folder").clicked() {
                                menu = Some((e.wf.id.clone(), "folder"));
                                ui.close();
                            }
                        });
                    }
                }
            });
        }

        if let Some((id, what)) = menu {
            let pinned = self.shell.entries.iter().find(|e| e.wf.id == id).is_some_and(|e| e.wf.pinned);
            let dir = self.shell.entries.iter().find(|e| e.wf.id == id).map(|e| e.dir.clone());
            match what {
                "pin" => self.edit_workflow(&id, |d| d.set_pinned(!pinned)),
                "done" => self.edit_workflow(&id, |d| d.set_status(Status::Done)),
                "active" => self.edit_workflow(&id, |d| d.set_status(Status::Active)),
                _ => {
                    if let Some(d) = dir {
                        self.reveal(&d);
                    }
                }
            }
        }
        if let Some(id) = open {
            self.open_tab(&id, true);
        }
    }
}

fn filter_menu(ui: &mut egui::Ui, label: &str, value: &mut String, options: &[String]) {
    let text = if value.is_empty() { label.to_string() } else { format!("{label}: {value}") };
    ui.menu_button(text, |ui| {
        if ui.selectable_label(value.is_empty(), "Any").clicked() {
            value.clear();
            ui.close();
        }
        for o in options {
            if ui.selectable_label(value == o, o).clicked() {
                *value = o.clone();
                ui.close();
            }
        }
    });
}

fn workflow_card(ui: &mut egui::Ui, e: &Entry, spark: &[f64], width: f32) -> egui::Response {
    let ctx = ui.ctx().clone();
    let colors = [theme::accent(ui.ctx()), Color32::from_rgb(0x3A, 0x4A, 0x9B), Color32::from_rgb(0xC9, 0x7F, 0x1F), Color32::from_rgb(0x2F, 0x7A, 0x66)];
    let color = colors[e.wf.id.bytes().last().unwrap_or(0) as usize % colors.len()];
    let ok_last = e.wf.experiments.last().is_some_and(|x| x.runs.ok > 0 && x.runs.failed == 0);
    let r = egui::Frame::new().fill(theme::card(&ctx)).stroke(egui::Stroke::new(1.0, theme::line(&ctx))).corner_radius(12).inner_margin(egui::Margin::same(14)).show(ui, |ui| ui.vertical(|ui| {
        ui.set_width(width - 30.0);
        ui.set_min_height(150.0);
        ui.horizontal(|ui| {
            kit::glyph(ui, if ok_last { kit::State::Done } else { kit::State::Idle });
            ui.add(egui::Label::new(RichText::new(&e.wf.title).strong().size(16.0)).truncate());
            if e.wf.pinned {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(RichText::new("pinned").size(12.0).color(theme::accent(&ctx))));
            }
        });
        ui.add(egui::Label::new(kit::muted(ui, e.wf.summary())).truncate());
        ui.add_space(8.0);
        kit::sparkline(ui, spark, color, Vec2::new(width - 30.0, 34.0));
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            let n = e.wf.experiments.len();
            ui.label(kit::muted(ui, format!("{} · {}", crate::fmt::plural(n, "experiment"), crate::fmt::plural(e.wf.run_count(), "run"))));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(kit::muted(ui, ago(e.wf.touched_secs(), now_secs()))));
        });
    }));
    let resp = ui.interact(r.response.rect, ui.scope_id().with(("card", &e.wf.id)), Sense::click());
    if resp.hovered() {
        ui.painter().rect_stroke(r.response.rect, 12.0, egui::Stroke::new(1.5, theme::accent(&ctx)), egui::StrokeKind::Inside);
    }
    resp
}
