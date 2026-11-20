// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The Sweep section of a workflow's Build page (UX plan U5, ADR 88): vary one model parameter over a range, see what it costs before it
//! starts, watch the points come in and keep each finished point as a run of the workflow. Everything is the engine's own: the plan is the engine's
//! dry run of the command (`check_list`, the same count `nyx check` prints), the points are `sweep_run` jobs of the queue at priority 3 (so a
//! Compute asked for meanwhile goes before them), the table is `sweep_table`, Cancel is `sweep_cancel`. The window adds only what it must: the
//! engine of the window has no workspace link (the shell writes its own run records), so each finished point is recorded here as a run of the
//! experiment, with the metrics the engine puts on the row (the keys the runs table and its plots use). The sweep varies a *geometry
//! parameter* (a number in the model's parameter list); a study value or a material is a sweep an assistant or a command list can ask for, the
//! window does not offer it yet.

use crate::app::App;
use crate::fmt;
use crate::kit::{self, State};
use crate::shell::{params_of, OpenWf};
use crate::theme;
use eframe::egui::{self, RichText};
use nyx_engine::facade::store::runs::{RunRecord, RunStatus};
use serde_json::{json, Value};
use std::collections::BTreeMap;

/// What the engine's dry run says a sweep will cost.
#[derive(Clone, Debug, PartialEq)]
pub struct Plan {
    pub points: usize,
    pub meshes: usize,
    pub solves: usize,
    /// why it cannot start, in the engine's words
    pub problems: Vec<String>,
}

impl Plan {
    /// "5 points: 5 meshes, 5 solves"
    pub fn words(&self) -> String {
        format!("{}: {}, {}", fmt::plural(self.points, "point"), if self.meshes == 1 { "1 mesh".to_string() } else { format!("{} meshes", self.meshes) }, fmt::plural(self.solves, "solve"))
    }
}

/// The form of the section, one per open workflow.
#[derive(Clone, Debug)]
pub struct Form {
    /// the section is the one on show in the inspector
    pub shown: bool,
    pub param: String,
    pub from: String,
    pub to: String,
    pub n: usize,
    /// the plan for the last values that were planned, with the text it was planned for
    plan: Option<(String, Plan)>,
    /// the metric the chart shows ("" = the one the Analysis page plots)
    pub metric: String,
}

impl Default for Form {
    fn default() -> Self {
        Form { shown: false, param: String::new(), from: String::new(), to: String::new(), n: 5, plan: None, metric: String::new() }
    }
}

/// `n` values from `from` to `to`, evenly spaced, ends included (one value when `n` is 1).
pub fn values(from: f64, to: f64, n: usize) -> Vec<f64> {
    match n {
        0 => vec![],
        1 => vec![from],
        _ => (0..n).map(|i| from + (to - from) * i as f64 / (n - 1) as f64).collect(),
    }
}

/// A number as typed: plain, with an exponent ("450e-9"), or with a length unit ("450 nm", "1.5 um").
pub fn parse_number(s: &str) -> Option<f64> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        return v.is_finite().then_some(v);
    }
    // the exponent is part of the text, so "450 nm" is exactly 450e-9 and not 450 times a rounded 1e-9
    for (unit, exp) in [("nm", -9), ("um", -6), ("µm", -6), ("mm", -3), ("cm", -2), ("m", 0)] {
        if let Some(num) = s.strip_suffix(unit) {
            let num = num.trim();
            if num.parse::<f64>().is_ok() && !num.contains(['e', 'E']) {
                return format!("{num}e{exp}").parse::<f64>().ok().filter(|v| v.is_finite());
            }
        }
    }
    None
}

/// One point of a running or finished sweep.
#[derive(Clone, Debug)]
pub struct Pt {
    pub value: f64,
    /// the engine's state: queued, running, done, failed, cancelled
    pub state: String,
    pub stage: String,
    pub seconds: f64,
    pub metrics: BTreeMap<String, f64>,
    pub error: Option<String>,
    pub recorded: bool,
}

impl Pt {
    fn finished(&self) -> bool {
        matches!(self.state.as_str(), "done" | "failed" | "cancelled")
    }

    pub fn glyph(&self) -> State {
        match self.state.as_str() {
            "done" => State::Done,
            "running" => State::Running,
            "failed" => State::Failed,
            "cancelled" => State::Idle,
            _ => State::Waiting,
        }
    }

    pub fn word(&self) -> &'static str {
        match self.state.as_str() {
            "done" => "Done",
            "running" => "Running",
            "failed" => "Failed",
            "cancelled" => "Cancelled",
            _ => "Waiting",
        }
    }
}

/// A sweep that was started from the window: where its runs go and how its points stand.
pub struct Active {
    pub wf_id: String,
    pub exp: String,
    pub id: u64,
    pub param: String,
    /// the project the sweep started from, as the shell stores a spec (each point's spec is this with the parameter set)
    base: Value,
    pub points: Vec<Pt>,
}

impl Active {
    pub fn finished(&self) -> usize {
        self.points.iter().filter(|p| p.finished()).count()
    }

    pub fn complete(&self) -> bool {
        self.finished() == self.points.len()
    }

    /// The spec a point solved: the sweep's project with the parameter set to the point's value.
    pub fn point_spec(&self, value: f64) -> String {
        let mut v = self.base.clone();
        if let Some(list) = v["geometry"]["params"].as_array_mut() {
            for p in list.iter_mut() {
                if p["name"] == json!(self.param) {
                    p["value"] = json!(value);
                }
            }
        }
        v.to_string()
    }

    /// The run record of a finished point (`None` for a point that was cancelled before it started: it never ran).
    pub fn record_of(&self, i: usize) -> Option<(RunRecord, String)> {
        let p = self.points.get(i)?;
        let status = match p.state.as_str() {
            "done" => RunStatus::Ok,
            "failed" => RunStatus::Failed,
            _ => return None,
        };
        let spec = self.point_spec(p.value);
        let rec = RunRecord { sec: p.seconds, status, params: params_of(&spec), metrics: p.metrics.clone(), error: p.error.clone(), ..Default::default() };
        Some((rec, spec))
    }
}

/// The numeric parameters of the model: name and value (a parameter that is an expression over others is not swept).
pub fn numeric_params(list: &nyx_engine::facade::geom::model::FeatureList) -> Vec<(String, f64)> {
    list.params
        .iter()
        .filter_map(|p| match &p.value {
            nyx_engine::facade::geom::expr::Expr::Num(v) => Some((p.name.clone(), *v)),
            _ => None,
        })
        .collect()
}

impl App {
    fn sweep_args(&self, f: &Form) -> Option<Value> {
        let (from, to) = (parse_number(&f.from)?, parse_number(&f.to)?);
        let vals = values(from, to, f.n);
        if f.param.is_empty() || vals.is_empty() {
            return None;
        }
        let points: Vec<Value> = vals
            .iter()
            .map(|v| json!({ "label": format!("{} = {}", f.param, fmt::si(*v, "m")), "edits": [{ "cmd": "set_param", "args": { "name": f.param, "value": v } }] }))
            .collect();
        Some(json!({ "points": points, "name": format!("{} sweep", f.param) }))
    }

    /// The engine's dry run of the sweep the form describes (planned again only when the form changed).
    fn sweep_plan(&self, f: &mut Form) -> Option<Plan> {
        let args = self.sweep_args(f)?;
        let text = json!({ "cmd": "sweep_run", "args": args }).to_string();
        if let Some((t, p)) = &f.plan {
            if *t == text {
                return Some(p.clone());
            }
        }
        let rep = self.engine.check_list(&text);
        let plan = Plan { points: f.n, meshes: rep.meshes, solves: rep.solves, problems: rep.problems.iter().map(|p| p.message.clone()).collect() };
        f.plan = Some((text, plan.clone()));
        Some(plan)
    }

    /// Starts the sweep of the form: the engine takes the points as jobs of its queue.
    pub fn start_sweep(&mut self, wf_id: &str, f: &Form) {
        // a run needs a folder: a draft gets it now (the sweep is its first change)
        if !self.ensure_dir(wf_id) {
            self.set_error("A sweep records its runs in a workflow folder, and this workflow has none (is there a workspace?).".to_string());
            return;
        }
        self.flush();
        let Some(exp) = self.shell.open.get(wf_id).map(|w| w.exp.clone()) else { return };
        if !self.push_edits() {
            return;
        }
        let Some(args) = self.sweep_args(f) else { return };
        let Some(r) = self.cmd("sweep_run", args.clone()) else { return };
        let Some(id) = r["sweep"].as_u64() else { return };
        let base: Value = serde_json::from_str(&self.spec_json().unwrap_or_default()).unwrap_or(Value::Null);
        let points = args["points"]
            .as_array()
            .map(|ps| ps.iter().map(|p| p["edits"][0]["args"]["value"].as_f64().unwrap_or(f64::NAN)).collect::<Vec<_>>())
            .unwrap_or_default()
            .into_iter()
            .map(|value| Pt { value, state: "queued".into(), stage: String::new(), seconds: 0.0, metrics: BTreeMap::new(), error: None, recorded: false })
            .collect();
        self.sweep = Some(Active { wf_id: wf_id.to_string(), exp, id, param: f.param.clone(), base, points });
        self.status = format!("Sweep of {} started: {} queued", f.param, fmt::plural(f.n, "point"));
    }

    /// Reads the sweep's table from the engine (which applies the points that finished) and records each finished point as a run.
    pub fn poll_sweep(&mut self) {
        let Some(a) = self.sweep.as_ref() else { return };
        if a.complete() && a.points.iter().all(|p| p.recorded || p.state == "cancelled") {
            return;
        }
        let Ok(t) = self.engine.call("sweep_table", &json!({ "sweep": a.id })) else { return };
        let rows = t["points"].as_array().cloned().unwrap_or_default();
        let mut finished_now = vec![];
        if let Some(a) = self.sweep.as_mut() {
            for (i, r) in rows.iter().enumerate() {
                let Some(p) = a.points.get_mut(i) else { continue };
                p.state = r["state"].as_str().unwrap_or("queued").to_string();
                p.stage = r["stage"].as_str().unwrap_or("").to_string();
                p.seconds = r["seconds"].as_f64().unwrap_or(p.seconds);
                p.error = r["error"].as_str().map(String::from);
                if let Some(m) = r["metrics"].as_object() {
                    p.metrics = m.iter().filter_map(|(k, v)| v.as_f64().map(|x| (k.clone(), x))).collect();
                }
                if p.finished() && !p.recorded {
                    finished_now.push(i);
                }
            }
        }
        for i in finished_now {
            let Some(a) = self.sweep.as_ref() else { return };
            let (wf, exp) = (a.wf_id.clone(), a.exp.clone());
            let rec = a.record_of(i);
            if let Some((rec, spec)) = rec {
                self.record_sweep_run(&wf, &exp, rec, &spec);
            }
            if let Some(p) = self.sweep.as_mut().and_then(|a| a.points.get_mut(i)) {
                p.recorded = true;
            }
        }
    }

    fn record_sweep_run(&mut self, wf: &str, exp: &str, rec: RunRecord, spec: &str) {
        let Some(w) = self.shell.open.get_mut(wf) else { return };
        let Some(d) = w.dir.as_mut() else { return };
        match d.record_run(exp, rec, spec) {
            Ok(_) => {
                if w.exp == exp {
                    w.reload_runs();
                }
                self.refresh_entries();
            }
            Err(e) => self.set_error(format!("Could not record a point of the sweep: {e}")),
        }
    }

    /// Cancels the points of the sweep that have not finished (queued ones never run; a running one stops at its next check).
    pub fn cancel_sweep(&mut self) {
        let Some(id) = self.sweep.as_ref().map(|a| a.id) else { return };
        self.cmd("sweep_cancel", json!({ "sweep": id }));
    }

    /// The section in the Build inspector.
    pub fn sweep_section(&mut self, ui: &mut egui::Ui, w: &mut OpenWf, acts: &mut Vec<crate::wfview::Act>) {
        let ctx = ui.ctx().clone();
        let params = numeric_params(&self.geom_edit);
        let f = &mut w.sweep;
        ui.label(RichText::new("Sweep").strong().size(16.0));
        ui.label(kit::muted(ui, "Vary one parameter of the model over a range. Each point is a run of this workflow; the points wait in the queue behind anything you start yourself."));
        ui.add_space(6.0);
        if params.is_empty() {
            ui.label(kit::muted(ui, "This model has no numeric parameter to vary. Add one in the Concept section (a named number the geometry uses)."));
            return;
        }
        if !params.iter().any(|(n, _)| *n == f.param) {
            let (n, v) = &params[0];
            f.param = n.clone();
            f.from = format!("{:.6e}", v * 0.8);
            f.to = format!("{:.6e}", v * 1.2);
        }
        ui.horizontal(|ui| {
            ui.label("Parameter");
            egui::ComboBox::from_id_salt("sweep-param").selected_text(f.param.clone()).show_ui(ui, |ui| {
                for (n, v) in &params {
                    if ui.selectable_label(*n == f.param, format!("{n}  (now {v:e})")).clicked() && *n != f.param {
                        f.param = n.clone();
                        f.from = format!("{:.6e}", v * 0.8);
                        f.to = format!("{:.6e}", v * 1.2);
                    }
                }
            });
        });
        ui.horizontal(|ui| {
            ui.label("From");
            ui.add(egui::TextEdit::singleline(&mut f.from).desired_width(96.0));
            ui.label("to");
            ui.add(egui::TextEdit::singleline(&mut f.to).desired_width(96.0));
            ui.label("in");
            ui.add(egui::DragValue::new(&mut f.n).range(1..=crate::sweepui::MAX_POINTS));
            ui.label("points");
        });
        let mut form = f.clone();
        let plan = self.sweep_plan(&mut form);
        w.sweep = form;
        let f = &w.sweep;
        match (&plan, parse_number(&f.from), parse_number(&f.to)) {
            (_, None, _) | (_, _, None) => {
                ui.label(RichText::new("Type each end as a number (450e-9) or with a unit (450 nm).").color(theme::activity_text(&ctx)));
            }
            (Some(p), _, _) if !p.problems.is_empty() => {
                for m in &p.problems {
                    ui.label(RichText::new(m).color(theme::danger(&ctx)));
                }
            }
            (Some(p), _, _) => {
                ui.label(RichText::new(p.words()).strong());
            }
            _ => {}
        }
        let can = plan.as_ref().is_some_and(|p| p.problems.is_empty()) && self.sweep.as_ref().is_none_or(|a| a.complete());
        let mut start = false;
        ui.horizontal(|ui| {
            if ui.add_enabled(can, egui::Button::new("Run sweep")).on_disabled_hover_text("A sweep is still running, or the plan has a problem").clicked() {
                start = true;
            }
            if self.sweep.as_ref().is_some_and(|a| !a.complete()) && kit::danger(ui, "Cancel the sweep").clicked() {
                self.cancel_sweep();
            }
        });
        if start {
            acts.push(crate::wfview::Act::Sweep(w.sweep.clone()));
        }
        let Some(a) = self.sweep.as_ref().filter(|a| a.wf_id == w.id) else { return };
        ui.add_space(8.0);
        ui.separator();
        ui.label(RichText::new(format!("{} of {} done", a.finished(), a.points.len())).strong());
        egui::Grid::new("sweep-table").striped(true).num_columns(3).show(ui, |ui| {
            for p in &a.points {
                kit::status(ui, p.glyph(), p.word());
                ui.label(kit::mono(fmt::si(p.value, "m")));
                let note = match p.state.as_str() {
                    "running" if !p.stage.is_empty() => p.stage.clone(),
                    "done" => fmt::seconds(p.seconds),
                    "failed" => p.error.clone().map(|e| fmt::ellipsize(&e, 48)).unwrap_or_default(),
                    _ => String::new(),
                };
                ui.label(kit::muted(ui, note));
                ui.end_row();
            }
        });
        if a.finished() > 0 {
            ui.add_space(6.0);
            self.sweep_chart(ui, w);
        }
        if !a.complete() {
            ctx.request_repaint_after(std::time::Duration::from_millis(250));
        }
    }

    /// One chart of a metric against the parameter, a marker for every finished point.
    fn sweep_chart(&self, ui: &mut egui::Ui, w: &mut OpenWf) {
        let Some(a) = self.sweep.as_ref() else { return };
        let done: Vec<&Pt> = a.points.iter().filter(|p| p.state == "done").collect();
        let mut keys: Vec<&String> = done.iter().flat_map(|p| p.metrics.keys()).filter(|k| !matches!(k.as_str(), "dofs" | "n_solutions" | "seconds")).collect();
        keys.sort();
        keys.dedup();
        let varies = |k: &str| {
            let v: Vec<f64> = done.iter().filter_map(|p| p.metrics.get(k).copied()).collect();
            v.windows(2).any(|p| (p[0] - p[1]).abs() > 1e-12 * p[0].abs().max(1e-300))
        };
        if w.sweep.metric.is_empty() || !keys.iter().any(|k| **k == w.sweep.metric) {
            // the metric the Analysis page would plot first (neff, the eigenvalue), else the first that varies
            let preferred = ["neff_re", "eigenvalue_re", "value_1"].iter().find_map(|p| keys.iter().find(|k| k.as_str() == *p));
            w.sweep.metric = preferred.or_else(|| keys.iter().find(|k| varies(k))).or(keys.first()).map(|k| k.to_string()).unwrap_or_default();
        }
        if w.sweep.metric.is_empty() {
            return;
        }
        ui.horizontal(|ui| {
            ui.label("Metric");
            egui::ComboBox::from_id_salt("sweep-metric").selected_text(fmt::key_label(&w.sweep.metric)).show_ui(ui, |ui| {
                for k in &keys {
                    ui.selectable_value(&mut w.sweep.metric, k.to_string(), fmt::key_label(k));
                }
            });
        });
        let pts: Vec<[f64; 2]> = done.iter().filter_map(|p| p.metrics.get(&w.sweep.metric).map(|y| [p.value, *y])).collect();
        let c = theme::accent(ui.ctx());
        egui_plot::Plot::new(("sweep_plot", &w.id)).height(180.0).x_axis_label(a.param.clone()).allow_scroll(false).show(ui, |p| {
            if pts.len() > 1 {
                p.line(egui_plot::Line::new("sweep", egui_plot::PlotPoints::from(pts.clone())).color(c));
            }
            p.points(egui_plot::Points::new("sweep", egui_plot::PlotPoints::from(pts)).radius(4.0).color(c).shape(egui_plot::MarkerShape::Circle));
        });
    }
}

/// The most points one sweep may have (the engine's own limit).
pub const MAX_POINTS: usize = 500;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_are_even_and_include_both_ends() {
        assert_eq!(values(1.0, 3.0, 5), vec![1.0, 1.5, 2.0, 2.5, 3.0]);
        assert_eq!(values(2.0, 9.0, 1), vec![2.0]);
        assert!(values(0.0, 1.0, 0).is_empty());
    }

    #[test]
    fn numbers_are_read_with_exponents_and_units() {
        assert_eq!(parse_number("450e-9"), Some(450e-9));
        assert_eq!(parse_number(" 450 nm "), Some(450e-9));
        assert_eq!(parse_number("1.5um"), Some(1.5e-6));
        assert_eq!(parse_number("2 mm"), Some(2e-3));
        assert_eq!(parse_number("wide"), None);
        assert_eq!(parse_number("inf"), None);
    }

    #[test]
    fn a_plan_says_what_it_costs() {
        let p = Plan { points: 5, meshes: 5, solves: 5, problems: vec![] };
        assert_eq!(p.words(), "5 points: 5 meshes, 5 solves");
        let one = Plan { points: 1, meshes: 1, solves: 1, problems: vec![] };
        assert_eq!(one.words(), "1 point: 1 mesh, 1 solve");
    }

    #[test]
    fn a_point_spec_sets_the_swept_parameter_and_nothing_else() {
        let base = json!({ "geometry": { "params": [{ "name": "w", "value": 5.0e-7 }, { "name": "t", "value": 2.2e-7 }] } });
        let a = Active { wf_id: "w".into(), exp: "e".into(), id: 1, param: "w".into(), base, points: vec![] };
        let spec: Value = serde_json::from_str(&a.point_spec(4.5e-7)).unwrap();
        assert_eq!(spec["geometry"]["params"][0]["value"], json!(4.5e-7));
        assert_eq!(spec["geometry"]["params"][1]["value"], json!(2.2e-7));
    }
}
