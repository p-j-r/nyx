//! What the engine runs and what waits, as the window shows it (UX plan U5, ADR 88): the busy chip of the top bar, the queue popover under it
//! (cancel, priority and pause; no hold, decision D8) and the stale band. The numbers are the engine's own, read from `queue_list` (the entries an
//! assistant gets from the `queue_*` commands), not the window's own list of its jobs: a sweep submits its points straight to the engine and they
//! must show here too. Every state has a glyph and a word; amber means time and attention and nothing else.

use crate::app::App;
use crate::fmt;
use crate::kit::{self, State};
use crate::theme;
use eframe::egui::{self, Align2, Color32, RichText, Stroke, Vec2};
use serde_json::Value;

/// One job of the queue that has not finished.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub id: u64,
    pub label: String,
    pub running: bool,
    /// the stage the solver reports ("" before it reports one)
    pub stage: String,
    pub seconds: f64,
    /// place among the waiting jobs, 1 for the next to start
    pub position: u64,
    pub priority: u8,
    pub owner: String,
    pub cancel_requested: bool,
}

/// The running and waiting jobs, running first, waiting in the order they will start.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Activity {
    pub rows: Vec<Row>,
    pub paused: bool,
}

impl Activity {
    /// From the entries of `queue_list` (finished, failed and cancelled jobs are not activity).
    pub fn of(entries: &[Value], paused: bool) -> Activity {
        let mut rows: Vec<Row> = entries
            .iter()
            .filter_map(|e| {
                let running = match e["state"].as_str()? {
                    "running" => true,
                    "queued" => false,
                    _ => return None,
                };
                Some(Row {
                    id: e["id"].as_u64()?,
                    label: sentence(e["label"].as_str().unwrap_or("")),
                    running,
                    stage: e["stage"].as_str().unwrap_or("").to_string(),
                    seconds: e["seconds"].as_f64().unwrap_or(0.0),
                    position: e["position"].as_u64().unwrap_or(0),
                    priority: e["priority"].as_u64().unwrap_or(2).clamp(1, 3) as u8,
                    owner: e["owner"].as_str().unwrap_or("").to_string(),
                    cancel_requested: e["cancel_requested"] == true,
                })
            })
            .collect();
        rows.sort_by_key(|r| (!r.running, r.position));
        Activity { rows, paused }
    }

    pub fn running(&self) -> Option<&Row> {
        self.rows.iter().find(|r| r.running)
    }

    pub fn waiting(&self) -> usize {
        self.rows.iter().filter(|r| !r.running).count()
    }

    /// Nothing runs, nothing waits.
    pub fn is_idle(&self) -> bool {
        self.rows.is_empty()
    }

    /// The state the chip shows: running if anything runs, else queued.
    pub fn state(&self) -> State {
        if self.running().is_some() {
            State::Running
        } else {
            State::Waiting
        }
    }

    /// The words of the chip, only what the engine reports: what runs and its stage, then how many wait ("Solving · assembling · 3 waiting").
    pub fn chip_text(&self) -> String {
        let mut parts: Vec<String> = vec![];
        match self.running() {
            Some(r) => {
                parts.push(r.label.clone());
                if r.cancel_requested {
                    parts.push("stopping".into());
                } else if !r.stage.is_empty() {
                    parts.push(r.stage.clone());
                }
            }
            None if self.paused => parts.push("Queue paused".into()),
            None => {}
        }
        let w = self.waiting();
        if w > 0 {
            parts.push(format!("{w} waiting"));
        }
        if parts.is_empty() {
            parts.push("Idle".into());
        }
        parts.join(" · ")
    }
}

/// The engine's job label ("solve ems") as a sentence start ("Solve ems").
fn sentence(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(f) => f.to_uppercase().chain(c).collect(),
        None => String::new(),
    }
}

/// Who asked, in words, for a row that is not the person's own (the window's engine only ever holds the person's jobs today; the queue entry
/// carries the owner, so a job of an assistant would say so).
pub fn asked_by(owner: &str) -> Option<&'static str> {
    match owner {
        "ai" => Some("asked by an assistant"),
        "cli" => Some("asked by a script"),
        _ => None,
    }
}

/// What the priority numbers mean to a person.
pub fn priority_word(p: u8) -> &'static str {
    match p {
        1 => "first",
        2 => "normal",
        _ => "later",
    }
}

impl App {
    /// The engine's queue as it is now.
    pub fn activity(&self) -> Activity {
        Activity::of(&self.engine.queue_list(None, None), self.engine.jobs().paused())
    }

    /// The busy chip: glyph, words and the counts, in the activity tint; a click opens the queue popover. Drawn only while something runs or waits.
    pub fn busy_chip(&mut self, ui: &mut egui::Ui, activity: &Activity) {
        if activity.is_idle() && !activity.paused && !self.show_queue {
            return;
        }
        let ctx = ui.ctx().clone();
        let text = activity.chip_text();
        let out = egui::Frame::new()
            .fill(theme::warm(&ctx))
            .stroke(Stroke::new(1.0, theme::activity_text(&ctx)))
            .corner_radius(16)
            .inner_margin(egui::Margin::symmetric(12, 5))
            .show(ui, |ui| {
                // painted by hand: the top bar lays out right to left, which would put the glyph after the words
                let font = egui::TextStyle::Body.resolve(ui.style());
                let galley = ui.painter().layout_no_wrap(text.clone(), font, theme::text(&ctx));
                let (rect, _) = ui.allocate_exact_size(Vec2::new(14.0 + 7.0 + galley.size().x, galley.size().y.max(14.0)), egui::Sense::hover());
                kit::paint_glyph(ui, egui::pos2(rect.left() + 7.0, rect.center().y), if activity.is_idle() { State::Idle } else { activity.state() });
                ui.painter().galley(egui::pos2(rect.left() + 21.0, rect.center().y - galley.size().y / 2.0), galley, theme::text(&ctx));
            });
        let r = ui.interact(out.response.rect, ui.scope_id().with("busy-chip"), egui::Sense::click());
        r.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, format!("Queue: {text}")));
        if r.on_hover_text("What runs and what waits, with cancel, priority and pause").clicked() {
            self.show_queue = !self.show_queue;
        }
    }

    /// The queue popover under the chip: one row per job (state, what it does, its stage or its place, priority, Cancel) and Pause / Resume.
    pub fn queue_panel(&mut self, ctx: &egui::Context) {
        if !self.show_queue {
            return;
        }
        let activity = self.activity();
        let mut act: Vec<(Action, u64)> = vec![];
        let area = egui::Area::new(egui::Id::unique("queue-popover"))
            .order(egui::Order::Foreground)
            .anchor(Align2::RIGHT_TOP, Vec2::new(-12.0, 52.0))
            .show(ctx, |ui| {
                egui::Frame::popup(ui.style()).fill(theme::card(ctx)).stroke(Stroke::new(1.0, theme::edge(ctx))).inner_margin(egui::Margin::same(14)).show(ui, |ui| {
                    ui.set_width(440.0);
                    ui.label(RichText::new("Queue").strong().size(16.0));
                    ui.add_space(4.0);
                    if activity.is_idle() {
                        ui.label(kit::muted(ui, "Nothing is running or waiting."));
                    }
                    for r in &activity.rows {
                        ui.horizontal(|ui| {
                            ui.spacing_mut().item_spacing.x = 8.0;
                            if r.running {
                                kit::status(ui, State::Running, "Running");
                            } else {
                                kit::status(ui, State::Waiting, format!("Waiting, number {}", r.position));
                            }
                            ui.label(RichText::new(&r.label).strong());
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                if ui.small_button("Cancel").clicked() {
                                    act.push((Action::Cancel, r.id));
                                }
                                if !r.running {
                                    let mut p = r.priority;
                                    egui::ComboBox::from_id_salt(("queue-priority", r.id)).selected_text(priority_word(p)).width(70.0).show_ui(ui, |ui| {
                                        for q in 1..=3u8 {
                                            ui.selectable_value(&mut p, q, priority_word(q));
                                        }
                                    });
                                    if p != r.priority {
                                        act.push((Action::Priority(p), r.id));
                                    }
                                }
                            });
                        });
                        let mut detail: Vec<String> = vec![];
                        if r.running {
                            detail.push(if r.stage.is_empty() { "starting".to_string() } else { r.stage.clone() });
                            detail.push(fmt::seconds(r.seconds));
                            if r.cancel_requested {
                                detail.push("cancel requested, it stops at its next check".into());
                            }
                        }
                        if let Some(who) = asked_by(&r.owner) {
                            detail.push(who.into());
                        }
                        if !detail.is_empty() {
                            ui.label(kit::muted(ui, detail.join(" · ")));
                        }
                        ui.add_space(4.0);
                    }
                    ui.separator();
                    ui.horizontal(|ui| {
                        if ui.button(if activity.paused { "Resume the queue" } else { "Pause the queue" }).on_hover_text("A paused queue starts nothing new; what runs finishes").clicked() {
                            act.push((Action::Pause, 0));
                        }
                        if activity.paused {
                            kit::status(ui, State::Waiting, "Paused: nothing new starts");
                        }
                    });
                });
            });
        // a click outside the popover (and not on the chip, which toggles it) or Esc closes it
        let esc = ctx.input(|i| i.key_pressed(egui::Key::Escape));
        let outside = ctx.input(|i| i.pointer.any_pressed()) && !area.response.contains_pointer() && ctx.input(|i| i.pointer.interact_pos()).is_some_and(|p| p.y > 50.0);
        for (what, id) in act {
            match what {
                Action::Cancel => {
                    if self.job.as_ref().is_some_and(|j| j.0 == id) {
                        self.cancel_job();
                    } else {
                        self.cancel_waiting(id);
                    }
                }
                Action::Priority(p) => {
                    self.engine.jobs().set_priority(id, p);
                }
                Action::Pause => {
                    if activity.paused {
                        self.engine.jobs().resume();
                    } else {
                        self.engine.jobs().pause();
                    }
                }
            }
        }
        if esc || outside {
            self.show_queue = false;
        }
        if !activity.is_idle() {
            ctx.request_repaint_after(std::time::Duration::from_millis(250));
        }
    }

    /// The stale band: one line in the activity tint with a glyph, the word and the reason, and a button to solve again when nothing else runs.
    /// Drawn above the result it is about (Results tab) and on the plot (viewport), so a stale field is never read as current.
    pub fn stale_band(&mut self, ui: &mut egui::Ui) {
        let Some(why) = self.current.and_then(|i| self.stale_reason(i)) else { return };
        let ctx = ui.ctx().clone();
        let mut again = false;
        egui::Frame::new().fill(theme::warm(&ctx)).stroke(Stroke::new(1.0, theme::activity_text(&ctx))).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 6)).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 7.0;
                kit::glyph(ui, State::Waiting);
                ui.label(RichText::new("Stale").strong().color(theme::activity_text(&ctx)));
                ui.label(RichText::new(format!("{why}.")).color(theme::text(&ctx)));
                if self.busy.is_none() && ui.small_button("Solve again").clicked() {
                    again = true;
                }
            });
        });
        if again {
            self.compute();
        }
    }
}

enum Action {
    Cancel,
    Priority(u8),
    Pause,
}

/// The colours of the band drawn on the plot itself: amber fill, the text of an amber fill.
pub fn band_colours(ctx: &egui::Context) -> (Color32, Color32) {
    (theme::activity(ctx), theme::on_activity(ctx))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn entry(id: u64, state: &str, label: &str, position: u64, priority: u64) -> Value {
        json!({ "id": id, "state": state, "label": label, "position": position, "priority": priority, "owner": "gui", "stage": "", "seconds": 0.0 })
    }

    #[test]
    fn finished_jobs_are_not_activity_and_running_comes_first() {
        let list = vec![
            entry(1, "done", "Solving", 0, 1),
            entry(2, "queued", "Sweep point 2", 2, 3),
            entry(3, "running", "Solving", 0, 1),
            entry(4, "queued", "Sweep point 1", 1, 3),
            entry(5, "cancelled", "Solving", 0, 1),
        ];
        let a = Activity::of(&list, false);
        assert_eq!(a.rows.iter().map(|r| r.id).collect::<Vec<_>>(), vec![3, 4, 2]);
        assert_eq!(a.waiting(), 2);
        assert_eq!(a.state(), State::Running);
        assert!(!a.is_idle());
        assert!(Activity::of(&list[..1], false).is_idle());
    }

    #[test]
    fn the_chip_says_only_what_the_engine_reports() {
        let mut running = entry(3, "running", "Solving", 0, 1);
        assert_eq!(Activity::of(std::slice::from_ref(&running), false).chip_text(), "Solving");
        running["stage"] = json!("factorising");
        let list = vec![running.clone(), entry(4, "queued", "Sweep point 1", 1, 3), entry(5, "queued", "Sweep point 2", 2, 3)];
        assert_eq!(Activity::of(&list, false).chip_text(), "Solving · factorising · 2 waiting");
        running["cancel_requested"] = json!(true);
        assert_eq!(Activity::of(&[running], false).chip_text(), "Solving · stopping");
        // nothing runs and the queue is paused: the chip says so; with jobs waiting it counts them
        assert_eq!(Activity::of(&[], true).chip_text(), "Queue paused");
        assert_eq!(Activity::of(&[entry(4, "queued", "Solving", 1, 1)], true).chip_text(), "Queue paused · 1 waiting");
        assert_eq!(Activity::of(&[], false).chip_text(), "Idle");
    }

    #[test]
    fn a_job_of_another_owner_says_who_asked() {
        assert_eq!(asked_by("ai"), Some("asked by an assistant"));
        assert_eq!(asked_by("gui"), None);
    }
}
