// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Feedback while editing and while a result arrives (productivity plan P3, ADR 92): the project checked with the edits that are not applied
//! yet, the size and cost of those edits from a trial mesh, and the polarization of a new result computed beside the frame. Everything here
//! answers a person's edit with a value or a line of text: it starts no mesh build, no solve, no job and no run, and it never replaces the
//! mesh or the results the window shows (the test `editing_every_form_starts_nothing` holds that).
//!
//! The trial mesh runs on a thread of its own, after the edit has rested for `REST`, and only for models whose mesh is cheap (the engine
//! refuses the others). When one takes longer than `SLOW` the window stops previewing by itself and says so, until the person asks.

use crate::app::{App, Msg};
use crate::kit::{self, State};
use eframe::egui::{self, RichText, Ui};
use nyx_engine::preview::SizePreview;
use std::time::{Duration, Instant};

/// An edit rests this long before a trial mesh is made of it.
pub const REST: Duration = Duration::from_millis(200);
/// A trial mesh that takes longer than this (seconds) turns the automatic preview off: the rule of plan question 2 (about 100 ms).
pub const SLOW: f64 = 0.1;

/// What the size line shows.
#[derive(Clone, Debug, Default)]
pub enum SizeView {
    /// nothing asked for yet
    #[default]
    None,
    /// a trial mesh is running
    Pending,
    Ready(SizePreview),
    /// the engine does not preview this model (the reason in words)
    Refused(String),
    /// the last trial mesh took `SLOW` or more; previewing waits for the person
    Slow(f64),
}

/// The state of the feedback for the edits in the Geometry tab.
#[derive(Default)]
pub struct Feedback {
    /// the first problem of the project with the edits in place (`Engine::check_edits`), in words
    pub problem: Option<String>,
    pub size: SizeView,
    /// a trial mesh is wanted for the edit as it is now once this time has passed
    due: Option<Instant>,
    /// counts the edits; an answer for an older count is dropped
    pub(crate) generation: u64,
    /// previewing was turned off by a slow mesh; the person turns it on again
    paused: bool,
    /// the edits as they were first seen after a project was adopted have had their first check
    seeded: bool,
}

impl Feedback {
    /// A project was adopted (loaded, new, an example): the edits are the project's own, ask again.
    pub fn reset(&mut self) {
        *self = Feedback { generation: self.generation + 1, ..Default::default() };
    }
}

impl App {
    /// The edits of the Geometry tab changed (or were seen for the first time): check them at once, and ask for a trial mesh once they rest.
    pub fn edits_changed(&mut self, ctx: &egui::Context) {
        let f = &mut self.feedback;
        f.generation += 1;
        f.seeded = true;
        f.problem = self.engine.check_edits(&self.geom_edit, &self.mesh_edit).err().map(|e| e.to_string());
        if f.problem.is_some() {
            f.size = SizeView::None;
            f.due = None;
            return;
        }
        if f.paused {
            return;
        }
        f.size = SizeView::Pending;
        f.due = Some(Instant::now() + REST);
        ctx.request_repaint_after(REST);
    }

    /// Called every frame the Geometry tab is drawn: the first look at a project's edits, and the trial mesh when an edit has rested.
    pub fn feedback_tick(&mut self, ctx: &egui::Context) {
        if !self.feedback.seeded {
            self.edits_changed(ctx);
        }
        let Some(due) = self.feedback.due else { return };
        if Instant::now() < due {
            ctx.request_repaint_after(due - Instant::now());
            return;
        }
        self.feedback.due = None;
        self.start_size_preview();
    }

    /// Starts the trial mesh of the edits on its own thread (also what the "Estimate size" button does after a slow one).
    pub fn start_size_preview(&mut self) {
        let task = match self.engine.size_preview(self.geom_edit.clone(), self.mesh_edit.clone()) {
            Ok(t) => t,
            Err(e) => {
                self.feedback.size = SizeView::Refused(e.to_string());
                return;
            }
        };
        let generation = self.feedback.generation;
        self.feedback.size = SizeView::Pending;
        let (tx, ctx) = (self.tx.clone(), self.ctx.clone());
        std::thread::spawn(move || {
            let out = task.run().map_err(|e| e.to_string());
            let _ = tx.send((generation, Msg::Size(out)));
            ctx.request_repaint();
        });
    }

    /// An answer of a trial mesh came back.
    pub fn take_size(&mut self, generation: u64, out: Result<SizePreview, String>) {
        if generation != self.feedback.generation {
            return; // the edit moved on; its own answer is on the way
        }
        self.feedback.size = match out {
            Ok(p) if p.mesh_seconds >= SLOW => {
                self.feedback.paused = true;
                SizeView::Slow(p.mesh_seconds)
            }
            Ok(p) => SizeView::Ready(p),
            Err(e) => SizeView::Refused(e),
        };
    }

    /// The polarization of the result that was just taken, computed beside the frame and handed back through the window's channel.
    pub fn start_polarization(&mut self, physics: &str) {
        let Some((result, task)) = self.engine.polarization_task(physics) else { return };
        let (tx, ctx, physics) = (self.tx.clone(), self.ctx.clone(), physics.to_string());
        std::thread::spawn(move || {
            let analysis = task().map_err(|e| e.to_string());
            let _ = tx.send((0, Msg::Polarization { physics, result, analysis }));
            ctx.request_repaint();
        });
    }

    /// The analysis arrived: it belongs to the result view that holds the very result it was made from, if that is still on screen.
    pub fn take_polarization(&mut self, physics: &str, result: &std::sync::Arc<nyx_engine::facade::physics::SolveResult>, analysis: Result<nyx_engine::analysis::ModeAnalysis, String>) {
        let Some(rv) = self.results.iter_mut().find(|r| r.physics.as_deref() == Some(physics) && std::sync::Arc::ptr_eq(&r.result, result)) else { return };
        match analysis {
            Ok(a) => rv.analysis = Some(a),
            Err(e) => self.status = format!("Polarization: {e}"),
        }
    }

    /// The problem line and the size line above the Build mesh button.
    pub fn feedback_lines(&mut self, ui: &mut Ui) {
        let mut ask = false;
        let problem = self.feedback.problem.clone();
        let size = self.feedback.size.clone();
        ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), 40.0), egui::Layout::top_down(egui::Align::Min), |ui| {
            ui.set_min_height(40.0);
            if let Some(p) = problem {
                kit::status(ui, State::Failed, format!("Problem: {p}"));
                return;
            }
            let muted = |ui: &Ui, t: String| kit::muted(ui, t).size(12.5);
            match size {
                SizeView::None | SizeView::Pending => {
                    ui.label(muted(ui, "Counting cells…".into()));
                }
                SizeView::Refused(why) => {
                    ui.label(muted(ui, format!("No size preview: {why}")));
                }
                SizeView::Slow(s) => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(muted(ui, format!("The trial mesh took {}: the preview waits for you.", crate::fmt::seconds(s))));
                        ask = ui.small_button("Estimate size").clicked();
                    });
                }
                SizeView::Ready(p) => {
                    ui.label(RichText::new(size_line(&p)).size(12.5));
                    ui.label(muted(ui, cost_line(&p)));
                }
            }
        });
        if ask {
            self.feedback.paused = false;
            self.start_size_preview();
        }
    }
}

/// "About 2,442 cells, 2D".
pub fn size_line(p: &SizePreview) -> String {
    format!("About {} cells, {}D (a trial mesh of these settings, not built)", group(p.cells), p.dim)
}

/// What a solve of that size costs on this machine, from its own history, or why there is no number.
pub fn cost_line(p: &SizePreview) -> String {
    match &p.estimate {
        None => "No cost estimate yet: it needs two earlier solves of different sizes.".into(),
        Some(e) => {
            let mut s = format!("A solve needs about {} of memory", crate::fmt::si(e.bytes as f64, "B"));
            if let Some(t) = e.seconds {
                s += &format!(" and about {}", crate::fmt::seconds(t));
            }
            s += &format!(" (from {}", crate::fmt::plural(e.runs, "earlier solve"));
            if e.extrapolated {
                s += ", beyond the sizes seen";
            }
            s + ")"
        }
    }
}

/// 2442 as "2,442".
fn group(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}
