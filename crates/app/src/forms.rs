// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Forms drawn from descriptors: a study, a condition or the options of a physics are a list of `ParamDesc`s (name, label, unit, type and
//! range, default), and the form for them is made here, the same for every physics. A new physics needs no change in the window.

use crate::tabs::num;
use eframe::egui::{self, DragValue, Ui};
use nyx_engine::facade::physics::{ParamDesc, ParamKind};
use serde_json::{json, Map, Value};

/// The JSON value a descriptor gives a parameter that was not set.
pub fn default_value(d: &ParamDesc) -> Value {
    match &d.kind {
        ParamKind::Number { default, .. } => json!(default),
        ParamKind::Integer { default, .. } => json!(default),
        ParamKind::Bool { default } => json!(default),
        ParamKind::Choice { default, .. } => json!(default),
        ParamKind::Vector { default, .. } => json!(default),
        ParamKind::Text { default } => json!(default),
        ParamKind::Names { default } => json!(default),
    }
}

/// Makes `value` an object holding every required parameter (with its default): a study or condition that is created must be complete.
pub fn fill_required(descs: &[ParamDesc], value: &mut Value) {
    if !value.is_object() {
        *value = Value::Object(Map::new());
    }
    let m = value.as_object_mut().expect("an object");
    for d in descs {
        if d.required && !m.contains_key(&d.name) {
            m.insert(d.name.clone(), default_value(d));
        }
    }
}

fn label(d: &ParamDesc) -> String {
    if d.unit.is_empty() {
        d.label.clone()
    } else {
        format!("{} [{}]", d.label, d.unit)
    }
}

/// Draws one row per parameter and writes what the person changes into `value` (an object; `Null` becomes one). A parameter that was not
/// touched stays absent, so the physics' default applies. True when something changed.
pub fn param_form(ui: &mut Ui, salt: &str, descs: &[ParamDesc], value: &mut Value) -> bool {
    if descs.is_empty() {
        return false;
    }
    if !value.is_object() {
        *value = Value::Object(Map::new());
    }
    let mut changed = false;
    egui::Grid::new(format!("form_{salt}")).num_columns(2).spacing([10.0, 4.0]).show(ui, |ui| {
        for d in descs {
            let given = value.get(&d.name).cloned();
            let mut set: Option<Value> = None;
            let resp_label = ui.label(label(d));
            if !d.doc.is_empty() {
                resp_label.on_hover_text(&d.doc);
            }
            match &d.kind {
                ParamKind::Number { default, min, max } => {
                    let mut v = given.as_ref().and_then(Value::as_f64).unwrap_or(*default);
                    if num(ui, &mut v).changed() {
                        if let Some(lo) = min {
                            v = v.max(*lo);
                        }
                        if let Some(hi) = max {
                            v = v.min(*hi);
                        }
                        set = Some(json!(v));
                    }
                }
                ParamKind::Integer { default, min, max } => {
                    let mut v = given.as_ref().and_then(Value::as_i64).unwrap_or(*default);
                    let r = ui.add(DragValue::new(&mut v).range(min.unwrap_or(i64::MIN / 2)..=max.unwrap_or(i64::MAX / 2)));
                    if r.changed() {
                        set = Some(json!(v));
                    }
                }
                ParamKind::Bool { default } => {
                    let mut v = given.as_ref().and_then(Value::as_bool).unwrap_or(*default);
                    if ui.checkbox(&mut v, "").changed() {
                        set = Some(json!(v));
                    }
                }
                ParamKind::Choice { options, default } => {
                    let cur = given.as_ref().and_then(Value::as_str).unwrap_or(default).to_string();
                    let mut pick = cur.clone();
                    egui::ComboBox::from_id_salt(format!("choice_{salt}_{}", d.name)).selected_text(&cur).show_ui(ui, |ui| {
                        for o in options {
                            ui.selectable_value(&mut pick, o.clone(), o);
                        }
                    });
                    if pick != cur {
                        set = Some(json!(pick));
                    }
                }
                ParamKind::Vector { len, default } => {
                    let mut v: Vec<f64> = given.as_ref().and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_f64).collect()).filter(|v: &Vec<f64>| v.len() == *len).unwrap_or_else(|| default.clone());
                    let mut edited = false;
                    ui.horizontal(|ui| {
                        for x in v.iter_mut() {
                            edited |= num(ui, x).changed();
                        }
                    });
                    if edited {
                        set = Some(json!(v));
                    }
                }
                ParamKind::Text { default } => {
                    let mut s = given.as_ref().and_then(Value::as_str).unwrap_or(default).to_string();
                    if ui.add(egui::TextEdit::singleline(&mut s).desired_width(180.0)).changed() {
                        set = Some(json!(s));
                    }
                }
                ParamKind::Names { default } => {
                    let cur: Vec<String> = given.as_ref().and_then(Value::as_array).map(|a| a.iter().filter_map(|x| x.as_str().map(str::to_string)).collect()).unwrap_or_else(|| default.clone());
                    let mut s = cur.join(", ");
                    if ui.add(egui::TextEdit::singleline(&mut s).desired_width(180.0).hint_text("region names, comma separated")).changed() {
                        let names: Vec<String> = s.split(',').map(|x| x.trim().to_string()).filter(|x| !x.is_empty()).collect();
                        set = Some(json!(names));
                    }
                }
            }
            if let Some(v) = set {
                value.as_object_mut().expect("an object").insert(d.name.clone(), v);
                changed = true;
            }
            ui.end_row();
        }
    });
    changed
}
