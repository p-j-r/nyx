// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Inspector tabs (geometry, materials, physics, study, results, analysis) and the results panel. The Materials, Physics and Study tabs are
//! drawn from the descriptors of the physics (`forms`), so they are the same for every physics; every change is handed to the engine as a
//! command (`set_problem`, `set_mesh`, `add_condition`, ...).

use crate::app::App;
use crate::colormap::Colormap;
use crate::forms;
use crate::theme;
use crate::state::*;
use eframe::egui::{self, Color32, DragValue, RichText, Ui};
use nyx_engine::facade::geom::model::{FeatureList, Model};
use nyx_engine::facade::physics::{Applies, MaterialAssign, ProblemSpec, PropKind};
use nyx_engine::SweepLayer;
use serde_json::{json, Value};
use std::collections::BTreeMap;

// ------------------------------------------------------------------ helpers

pub fn fmt_num(n: f64) -> String {
    if n == 0.0 {
        "0".into()
    } else if n.abs() >= 1e5 || n.abs() < 1e-3 {
        format!("{n:.5e}")
    } else {
        let s = format!("{n:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        s.to_string()
    }
}

/// Parses a float or scientific-notation number, with optional SI unit suffixes (`nm`, `um`, `µm`, `mm`, `THz`, `GPa`, ...).
pub fn parse_si_f64(s: &str) -> Option<f64> {
    let norm = s.trim().replace(',', ".");
    let s = norm.trim();
    if s.is_empty() {
        return None;
    }
    if let Ok(v) = s.parse::<f64>() {
        return v.is_finite().then_some(v);
    }
    const UNITS: [(&str, i32); 18] = [
        ("THz", 12),
        ("GHz", 9),
        ("MHz", 6),
        ("kHz", 3),
        ("GPa", 9),
        ("MPa", 6),
        ("kPa", 3),
        ("fm", -15),
        ("pm", -12),
        ("nm", -9),
        ("um", -6),
        ("µm", -6),
        ("μm", -6),
        ("mm", -3),
        ("cm", -2),
        ("Hz", 0),
        ("Pa", 0),
        ("m", 0),
    ];
    for (unit, exp) in UNITS {
        if let Some(num) = s.strip_suffix(unit) {
            let num = num.trim();
            if let Ok(base) = num.parse::<f64>() {
                if base.is_finite() && !num.contains(['e', 'E']) {
                    return format!("{num}e{exp}").parse::<f64>().ok().filter(|v| v.is_finite());
                }
            }
        }
    }
    None
}

pub fn select_all_text(ctx: &egui::Context, id: egui::Id, text: &str) {
    let mut state = egui::text_edit::TextEditState::load(ctx, id).unwrap_or_default();
    state.cursor.set_char_range(Some(egui::text::CCursorRange::two(
        egui::text::CCursor::new(0),
        egui::text::CCursor::new(text.chars().count()),
    )));
    state.store(ctx, id);
}

pub fn select_all_on_focus_gain(ctx: &egui::Context, r: &egui::Response, text: &str) {
    let flag = r.id.with("just_focused");
    if r.gained_focus() {
        select_all_text(ctx, r.id, text);
        ctx.data_mut(|d| d.insert_temp(flag, true));
    } else if r.clicked() && ctx.data_mut(|d| d.get_temp::<bool>(flag).unwrap_or(false)) {
        select_all_text(ctx, r.id, text);
        ctx.data_mut(|d| d.remove::<bool>(flag));
    } else if !r.has_focus() {
        ctx.data_mut(|d| d.remove::<bool>(flag));
    }
}

fn num_with_id(ui: &mut Ui, id: egui::Id, v: &mut f64) -> egui::Response {
    let focused = ui.memory(|m| m.has_focus(id));
    let cached = if focused { ui.data_mut(|d| d.get_temp::<String>(id)) } else { None };
    let mut text = cached.clone().unwrap_or_else(|| fmt_num(*v));
    if focused && cached.is_none() {
        select_all_text(ui.ctx(), id, &text);
    }
    let mut r = ui.add(
        egui::TextEdit::singleline(&mut text)
            .id(id)
            .min_size(egui::vec2(104.0, 28.0))
            .desired_width(116.0)
            .clip_text(false)
            .font(egui::TextStyle::Monospace),
    );
    select_all_on_focus_gain(ui.ctx(), &r, &text);
    if r.changed() || r.lost_focus() {
        if let Some(parsed) = parse_si_f64(&text) {
            if parsed != *v {
                *v = parsed;
                r.mark_changed();
            }
        }
    }
    if r.has_focus() {
        ui.data_mut(|d| d.insert_temp(id, text));
    } else if r.lost_focus() {
        ui.data_mut(|d| d.remove::<String>(id));
    }
    r
}

/// A single-click-focusable number field with a persistent ID: accepts scientific notation and SI unit suffixes (`220 nm`, `193.414 THz`).
pub fn num_id(ui: &mut Ui, salt: impl std::hash::Hash + std::fmt::Debug, v: &mut f64) -> egui::Response {
    let id = ui.make_persistent_id(salt);
    num_with_id(ui, id, v)
}

/// A number field: single-click and type (scientific notation and SI suffixes accepted).
#[allow(dead_code)]
pub fn num(ui: &mut Ui, v: &mut f64) -> egui::Response {
    let id = ui.next_auto_id();
    ui.skip_ahead_auto_ids(1);
    num_with_id(ui, id, v)
}

/// A plate `lx` wide (x) and `lz` high (z), one region `plate`: the geometry a person starts from.
fn plate(lx: f64, lz: f64) -> FeatureList {
    let mut m = Model::new_2d();
    let r = m.rect("plate_box", (0.0, 0.0), (lx, lz));
    m.overlay(&[("plate", r)]);
    m.list
}

/// Layers along y (the propagation axis) of a swept cross-section.
fn sweep_editor(ui: &mut Ui, sweep: &mut Option<Vec<SweepLayer>>) {
    let mut on = sweep.is_some();
    if ui.checkbox(&mut on, "Sweep along y into 3D layers").changed() {
        *sweep = if on { Some(vec![SweepLayer { thickness: 0.1, n: 2, domain_offset: 0, fill: None, regions: None }]) } else { None };
    }
    if let Some(layers) = sweep {
        let mut rm = None;
        for (i, l) in layers.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("layer {}", i + 1));
                ui.label("thickness");
                num_id(ui, ("sweep_layer", i), &mut l.thickness);
                ui.label("elements");
                ui.add(DragValue::new(&mut l.n).range(1..=200));
                ui.label("domain +");
                ui.add(DragValue::new(&mut l.domain_offset).range(0..=900));
                if ui.small_button("✕").clicked() {
                    rm = Some(i);
                }
            });
        }
        if let Some(i) = rm {
            if layers.len() > 1 {
                layers.remove(i);
            }
        }
        if ui.small_button("+ layer").clicked() {
            layers.push(SweepLayer { thickness: 0.1, n: 2, domain_offset: 0, fill: None, regions: None });
        }
        ui.weak("The ends of the swept model are the boundaries ymin (at y = 0) and ymax.");
    }
}

/// The key a material is stored under for a region: the name, or `name*` when the geometry cut absorbing pieces `name.pml_*` out of it, so
/// that one assignment covers the region and its layers.
fn material_key(base: &str, all: &[(i32, String)]) -> String {
    if all.iter().any(|(_, n)| n.strip_prefix(base).is_some_and(|r| r.starts_with('.'))) {
        format!("{base}*")
    } else {
        base.to_string()
    }
}

impl App {

    // ------------------------------------------------------------------ inspector

    pub fn inspector(&mut self, ui: &mut Ui) {
        ui.horizontal_wrapped(|ui| {
            for t in Tab::ALL {
                if t == Tab::Results && self.results.is_empty() {
                    continue;
                }
                if t == Tab::Analysis && !self.has_analysis() {
                    continue;
                }
                if ui.selectable_label(self.tab == t, RichText::new(t.name()).strong()).clicked() {
                    self.tab = t;
                    if t != Tab::Physics {
                        self.picking = None;
                    }
                }
            }
        });
        ui.separator();
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| match self.tab {
            Tab::Geometry => self.tab_geometry(ui),
            Tab::Materials => self.tab_materials(ui),
            Tab::Physics => self.tab_physics(ui),
            Tab::Study => self.tab_study(ui),
            Tab::Results => self.tab_results(ui),
            Tab::Analysis => self.tab_analysis(ui),
        });
    }

    /// The heading of a tab: the classic window needs it; in the window of the shell the row of tabs above names the tab already (ADR 124).
    fn tab_heading(&self, ui: &mut Ui, text: &str) {
        if !self.shell.active() {
            ui.heading(text);
        }
    }

    /// One inspector tab on its own: the workflow view lays the six tabs out in its three sections.
    pub fn inspector_tab(&mut self, ui: &mut Ui, t: Tab) {
        match t {
            Tab::Geometry => self.tab_geometry(ui),
            Tab::Materials => self.tab_materials(ui),
            Tab::Physics => self.tab_physics(ui),
            Tab::Study => self.tab_study(ui),
            Tab::Results => self.tab_results(ui),
            Tab::Analysis => self.tab_analysis(ui),
        }
    }

    // ------------------------------------------------------------------ geometry

    fn tab_geometry(&mut self, ui: &mut Ui) {
        self.tab_heading(ui, "Geometry");
        ui.weak("Axes: x lateral, y the direction of propagation, z the surface normal (right-handed). A 2D cross-section is drawn in x-z with y pointing into the page; crystal orientations, constraint and load components and the names of result fields use x, y, z.");
        let mut g = self.geom_edit.clone();
        let mut m = self.mesh_edit.clone();
        ui.horizontal(|ui| {
            ui.label("Start from");
            if ui.button("plate (2D)").clicked() {
                g = plate(1.0, 0.5);
                m.h = 0.05;
                m.sweep = None;
            }
            if ui.button("block (plate swept along y)").clicked() {
                g = plate(1.0, 0.5);
                m.h = 0.1;
                m.sweep = Some(vec![SweepLayer { thickness: 1.0, n: 8, domain_offset: 0, fill: None, regions: None }]);
            }
        });
        ui.horizontal(|ui| {
            ui.label("element size h (m)");
            num_id(ui, "mesh_h", &mut m.h);
        });
        let mesh_hint = if g.dim == 2 && m.sweep.as_ref().is_some_and(|s| !s.is_empty()) {
            "3D swept model (x-z cross-section along y); the size rules of the model refine it"
        } else if g.dim == 2 {
            "2D model, x-z cross-section; the size rules of the model refine it"
        } else {
            "3D model; the size rules of the model refine it"
        };
        ui.add(egui::Label::new(RichText::new(mesh_hint).weak()).wrap());
        ui.checkbox(&mut m.physics_sizing, "size the mesh from the physics")
            .on_hover_text("Elements are made fine enough for the wavelength of the study in each material (points per wavelength by element order). Applies to EM waves and solid-mechanics eigenfrequencies; the sizes of the model still apply.");
        if m.physics_sizing {
            let mut ppw = m.points_per_wavelength.unwrap_or(0.0);
            ui.horizontal(|ui| {
                ui.label("points per wavelength (0 = by element order)");
                if num_id(ui, "ppw", &mut ppw).changed() {
                    m.points_per_wavelength = (ppw > 0.0).then_some(ppw);
                }
            });
        }
        let mut adaptive = m.adapt_scale.is_some();
        ui.horizontal(|ui| {
            if ui.checkbox(&mut adaptive, "adaptive mesh").on_hover_text("Metric adaptation towards the sizes of the model: the boundary and interface vertices lie on the exact curves. The factor scales every size (below one is finer).").changed() {
                m.adapt_scale = adaptive.then_some(1.0);
            }
            if let Some(s) = &mut m.adapt_scale {
                ui.label("size factor");
                num_id(ui, "adapt_s", s);
            }
        });
        let cur_pml_mode = m.resolved_pml_mode();
        ui.horizontal(|ui| {
            ui.label("absorbing layers");
            for (val, label, tip) in [
                (
                    nyx_engine::PmlModeRequest::Auto,
                    "Auto",
                    "Structured native Quad/Prism/Hex layers whenever the geometry has absorbing layers (default).",
                ),
                (
                    nyx_engine::PmlModeRequest::Structured,
                    "Structured",
                    "Always build structured Quad/Prism/Hex absorbing layers.",
                ),
                (
                    nyx_engine::PmlModeRequest::Inline,
                    "Inline",
                    "Legacy unstructured simplex triangulation inside inward CAD-sliced PML domains.",
                ),
            ] {
                if ui.selectable_label(cur_pml_mode == val, label).on_hover_text(tip).clicked() {
                    m.pml_mode = val;
                    m.structured_pml = None;
                }
            }
        });
        sweep_editor(ui, &mut m.sweep);
        ui.add_space(4.0);
        crate::geomui::show(ui, &mut self.geom_ui, &mut g);
        ui.add_space(6.0);
        let changed = g != self.geom_edit || m != self.mesh_edit;
        self.geom_edit = g;
        self.mesh_edit = m;
        if changed {
            self.geometry_dirty = true;
            self.prune_stale_materials();
            let ctx = ui.ctx().clone();
            self.edits_changed(&ctx);
        }
        let ctx = ui.ctx().clone();
        self.feedback_tick(&ctx);
        self.feedback_lines(ui);
        ui.add_space(4.0);
        let dirty = self.geometry_dirty || self.geom_edit != self.engine.project().geometry || self.mesh_edit != self.engine.project().mesh;
        ui.horizontal(|ui| {
            let b = ui.add_enabled(self.busy.is_none(), egui::Button::new(RichText::new("Build mesh").strong().color(if dirty { crate::theme::on_accent(ui.ctx()) } else { crate::theme::text(ui.ctx()) })).fill(if dirty { crate::theme::accent(ui.ctx()) } else { crate::theme::raised(ui.ctx()) }));
            if b.clicked() {
                self.rebuild_mesh();
            }
            if dirty {
                ui.weak("(geometry or mesh settings changed)");
            }
        });

        ui.add_space(8.0);
        ui.separator();
        ui.heading("Mesh");
        if let Some(v) = self.doc.view.clone() {
            ui.label(format!("{}D, {} points, {} cells, {} boundary facets/segments", v.dim, v.points.len(), v.n_cells, if v.dim == 2 { v.segs.len() } else { v.tris.len() }));
            if let Some(model) = self.doc.model.clone() {
                let regions: Vec<String> = model.regions().iter().map(|(i, n)| format!("{n} ({i})")).collect();
                ui.label(format!("regions: {}", regions.join(", ")));
            }
            ui.label(format!("boundaries {:?}", v.boundary_labels));
            ui.label(format!("element types: {}", v.types.join(", ")));
            if let Some(r) = self.doc.model.as_ref().and_then(|m| m.report().cloned()) {
                if let Some(a) = r.adaptation {
                    ui.weak(format!("adaptation: {} iterations, {} to {} cells, worst element quality {:.2}", a.iterations, a.cells_start, a.cells_end, a.q_min));
                }
            }
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.show_edges, "edges");
                ui.checkbox(&mut self.show_numbers, "boundary numbers");
            });
        } else {
            ui.weak("no mesh yet");
        }
        ui.add_space(8.0);
        ui.separator();
        ui.heading("Examples");
        for (i, e) in crate::examples::list().iter().enumerate() {
            if ui.link(&e.title).clicked() {
                self.load_example(i);
            }
        }
    }

    // ------------------------------------------------------------------ materials

    /// The regions of the model by base name (the absorbing pieces `name.pml_*` belong to `name`), with the key a material is stored under.
    pub(crate) fn material_regions(&self) -> Vec<(String, String)> {
        let from_geom: Vec<(i32, String)> = nyx_engine::facade::geomdrive::evaluate(&self.geom_edit)
            .ok()
            .and_then(|ev| ev.label_map().ok())
            .map(|lm| {
                lm.domains
                    .into_iter()
                    .filter(|d| matches!(d.kind, nyx_engine::facade::geom::model::RegionKind::Material | nyx_engine::facade::geom::model::RegionKind::Pml))
                    .map(|d| (d.id, d.name))
                    .collect()
            })
            .unwrap_or_default();
        let all = if !from_geom.is_empty() { &from_geom } else { &self.domain_names };
        let mut seen: Vec<(String, String)> = vec![];
        for (_, n) in all {
            let base = n.split('.').next().unwrap_or(n).to_string();
            if !seen.iter().any(|(b, _)| *b == base) {
                let key = material_key(&base, all);
                seen.push((base, key));
            }
        }
        seen
    }

    /// Removes material assignments for region names that no longer exist in `geom_edit`.
    pub fn prune_stale_materials(&mut self) {
        let regions = self.material_regions();
        if regions.is_empty() {
            return;
        }
        let specs: Vec<(String, ProblemSpec)> = self
            .engine
            .project()
            .physics
            .iter()
            .filter_map(|s| {
                let mut sp = s.clone();
                let before = sp.materials.len();
                sp.materials.retain(|k, _| regions.iter().any(|(base, key)| k == key || k == base));
                (sp.materials.len() != before).then_some((sp.physics.clone(), sp))
            })
            .collect();
        for (phys, sp) in specs {
            self.set_problem(&phys, &sp);
        }
    }

    fn tab_materials(&mut self, ui: &mut Ui) {
        let Some(active) = self.active() else {
            ui.heading("Materials");
            ui.weak("The project has no physics: add one in the Physics tab.");
            return;
        };
        let Some(physics) = self.engine.registry().get(&active) else { return };
        let props = physics.properties();
        let pname = physics.name();
        let Some(spec0) = self.engine.project().spec(&active).cloned() else { return };
        let mut spec = spec0.clone();
        self.tab_heading(ui, "Materials");
        ui.label(format!("What each region of the model is made of, for {pname}. A material is a file of the library (units written in it, optical models evaluated at the wavelength of the study) or numbers typed here."));
        let regions = self.material_regions();
        if regions.is_empty() {
            ui.weak("Build the mesh to list the regions.");
            ui.add_space(6.0);
            self.material_import_panel(ui);
            return;
        }
        let lib_names: Vec<String> = self.engine.library().map(|l| l.names()).unwrap_or_default();
        if lib_names.is_empty() {
            ui.colored_label(crate::theme::activity_text(ui.ctx()), "No material files were found (set NYX_MATERIALS, or keep a `materials` folder next to the program).");
        }
        let mut apply_all: Option<String> = None;
        ui.horizontal(|ui| {
            ui.label("All regions:");
            egui::ComboBox::from_id_salt("mat_all").height(320.0).selected_text("apply material...").show_ui(ui, |ui| {
                for n in &lib_names {
                    if ui.selectable_label(false, n).clicked() {
                        apply_all = Some(n.clone());
                    }
                }
            });
        });
        if let Some(n) = apply_all {
            for (_, key) in &regions {
                spec.materials.insert(key.clone(), MaterialAssign::Named(n.clone()));
            }
        }
        for (base, key) in &regions {
            ui.group(|ui| {
                let cur = spec.materials.get(key).or_else(|| spec.materials.get(base)).cloned();
                ui.horizontal(|ui| {
                    ui.label(RichText::new(if key == base { base.clone() } else { format!("{base}  (and its absorbing layers)") }).strong());
                    let shown = match &cur {
                        Some(MaterialAssign::Named(n)) => n.clone(),
                        Some(MaterialAssign::Oriented { material, .. }) => format!("{material} (rotated)"),
                        Some(MaterialAssign::Props(_)) | Some(MaterialAssign::PropsTensor(_)) => "numbers".into(),
                        None => "none".into(),
                    };
                    egui::ComboBox::from_id_salt(format!("mat_{key}")).height(320.0).selected_text(shown).show_ui(ui, |ui| {
                        if ui.selectable_label(matches!(cur, Some(MaterialAssign::Props(_))), "numbers").clicked() && !matches!(cur, Some(MaterialAssign::Props(_))) {
                            let map: BTreeMap<String, f64> = props.iter().filter(|p| p.kind != PropKind::Matrix3x6).map(|p| (p.name.clone(), p.default.unwrap_or_else(|| crate::examples::typical(&p.name)))).collect();
                            spec.materials.insert(key.clone(), MaterialAssign::Props(map));
                        }
                        ui.separator();
                        for n in &lib_names {
                            if ui.selectable_label(false, n).clicked() {
                                spec.materials.insert(key.clone(), MaterialAssign::Named(n.clone()));
                            }
                        }
                    });
                });
                match spec.materials.get(key).cloned().or_else(|| spec.materials.get(base).cloned()) {
                    Some(MaterialAssign::Named(n)) => {
                        let mut rotated = false;
                        ui.horizontal(|ui| {
                            if ui.checkbox(&mut rotated, "crystal axes rotated in the wafer").changed() && rotated {
                                spec.materials.insert(key.clone(), MaterialAssign::Oriented { material: n.clone(), euler_deg: [0.0; 3] });
                            }
                        });
                        self.material_summary(ui, &n);
                    }
                    Some(MaterialAssign::Oriented { material, mut euler_deg }) => {
                        let mut keep = true;
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut keep, "crystal axes rotated in the wafer");
                        });
                        let mut edited = false;
                        ui.horizontal(|ui| {
                            ui.label("Bunge angles (degrees)");
                            for (i, a) in euler_deg.iter_mut().enumerate() {
                                edited |= num_id(ui, ("bunge_angle", key, i), a).changed();
                            }
                        });
                        if !keep {
                            spec.materials.insert(key.clone(), MaterialAssign::Named(material.clone()));
                        } else if edited {
                            spec.materials.insert(key.clone(), MaterialAssign::Oriented { material: material.clone(), euler_deg });
                        }
                        self.material_summary(ui, &material);
                    }
                    Some(MaterialAssign::Props(mut map)) => {
                        let mut edited = false;
                        egui::Grid::new(format!("props_{key}")).num_columns(3).show(ui, |ui| {
                            for p in &props {
                                if p.kind == PropKind::Matrix3x6 {
                                    continue;
                                }
                                let mut v = map.get(&p.name).copied().unwrap_or_else(|| p.default.unwrap_or_else(|| crate::examples::typical(&p.name)));
                                ui.label(&p.label);
                                if num_id(ui, ("mat_prop", key, &p.name), &mut v).changed() {
                                    map.insert(p.name.clone(), v);
                                    edited = true;
                                }
                                ui.weak(&p.unit);
                                ui.end_row();
                            }
                        });
                        if edited {
                            spec.materials.insert(key.clone(), MaterialAssign::Props(map));
                        }
                    }
                    Some(MaterialAssign::PropsTensor(map)) => {
                        for p in &props {
                            if let Some(v) = map.get(&p.name) {
                                ui.horizontal(|ui| {
                                    ui.label(&p.label);
                                    ui.monospace(format!("{v:?}"));
                                    ui.weak(&p.unit);
                                });
                            }
                        }
                    }
                    None => {
                        ui.weak("No material: choose one (a problem needs a material for every region it solves in).");
                    }
                }
            });
        }
        if spec != spec0 {
            self.set_problem(&active, &spec);
        }
        ui.add_space(6.0);
        ui.weak("Library values are typical room-temperature numbers with their sources written in the files; check them for your application. Anisotropic crystal data replaces the scalar value it covers in that region.");
        ui.add_space(6.0);
        self.material_import_panel(ui);
    }

    /// What a library material gives, in one line, live (T, f, d) evaluation, and its source file text on request.
    fn material_summary(&mut self, ui: &mut Ui, name: &str) {
        match self.engine.describe_material(name) {
            Ok(d) => {
                if let Some(desc) = d["description"].as_str().filter(|s| !s.is_empty()) {
                    ui.weak(desc);
                }
                let mut parts: Vec<String> = vec![];
                if let Some(k) = d["kind"].as_str() {
                    parts.push(k.to_string());
                }
                if let Some(pg) = d["point_group"].as_str() {
                    parts.push(format!("pg: {pg}"));
                }
                if let Some(c) = d["elastic"]["class"].as_str() {
                    parts.push(format!("elastic: {c}"));
                }
                if d["optical"] == true {
                    parts.push("optical".into());
                }
                if d["thermal"] == true {
                    parts.push("thermal".into());
                }
                if d["dielectric"] == true {
                    parts.push("dielectric".into());
                } else if d["electrical"] == true {
                    parts.push("electrical".into());
                }
                if d["superconducting"] == true {
                    parts.push("superconducting (BCS)".into());
                }
                if d["piezoelectric"] == true {
                    parts.push("piezoelectric".into());
                }
                if d["photoelastic"] == true {
                    parts.push("photoelastic".into());
                }
                if d["electro_optic"] == true {
                    parts.push("electro-optic".into());
                }
                if d["viscoelastic"] == true {
                    parts.push("viscoelastic".into());
                }
                if d["interfaces"] == true {
                    parts.push("TLS interfaces".into());
                }
                ui.weak(parts.join("  ·  "));
                egui::CollapsingHeader::new("Live (T, f, d) evaluator")
                    .id_salt(format!("eval_{name}"))
                    .default_open(false)
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label("T (K):");
                            num(ui, &mut self.mat_fetch.eval_temp_k);
                            ui.label("f (GHz):");
                            num(ui, &mut self.mat_fetch.eval_freq_ghz);
                            if d["superconducting"] == true {
                                ui.label("d (nm):");
                                num(ui, &mut self.mat_fetch.eval_thickness_nm);
                            }
                        });
                        let t_k = self.mat_fetch.eval_temp_k.max(0.001);
                        let f_hz = self.mat_fetch.eval_freq_ghz.max(1e-6) * 1.0e9;
                        let d_m = self.mat_fetch.eval_thickness_nm.max(0.1) * 1.0e-9;
                        if let Ok(ev) = self.engine.evaluate_material(name, Some(t_k), Some(f_hz), None, None, Some(d_m), None, None) {
                            if let Some(th) = ev.get("thermal") {
                                if let Some(cp) = th["specific_heat_j_kgk"].as_f64() {
                                    ui.weak(format!("Cp({t_k:.2} K) = {cp:.4e} J/(kg*K)"));
                                }
                            }
                            if let Some(dl) = ev.get("dielectric") {
                                let exx = dl["eps_re"][0][0].as_f64().unwrap_or(1.0);
                                let ezz = dl["eps_re"][2][2].as_f64().unwrap_or(exx);
                                let td = dl["tan_delta_xx"].as_f64().unwrap_or(0.0);
                                ui.weak(format!("eps_r = ({exx:.3}, {ezz:.3}), tan delta = {td:.3e}"));
                            }
                            if let Some(sc) = ev.get("superconducting") {
                                let rs = sc["r_s_ohm_sq"].as_f64().unwrap_or(0.0);
                                let lk = sc["l_k_sq_h_sq"].as_f64().unwrap_or(0.0) * 1.0e12;
                                let gap = sc["delta_t_mev"].as_f64().unwrap_or(0.0);
                                ui.weak(format!("BCS: R_s = {rs:.3e} Ohm/sq, L_k = {lk:.3} pH/sq, Delta(T) = {gap:.3} meV"));
                            }
                        }
                    });
                if let Some(t) = d["toml"].as_str() {
                    egui::CollapsingHeader::new("material file").id_salt(format!("toml_{name}")).default_open(false).show(ui, |ui| {
                        ui.add(egui::TextEdit::multiline(&mut t.to_string().as_str()).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
                    });
                }
            }
            Err(e) => {
                ui.colored_label(crate::theme::danger(ui.ctx()), e.to_string());
            }
        }
    }

    /// Collapsed panel to download curated presets, search open databases asynchronously, or convert pasted database payloads (`refractiveindex.info`, `AFLOW`, `Materials Project`, `NumBAT`).
    fn material_import_panel(&mut self, ui: &mut Ui) {
        use crate::state::MaterialBgMsg;
        if let Some(rx) = &self.mat_fetch.pending {
            ui.ctx().request_repaint();
            if let Ok(msg) = rx.try_recv() {
                self.mat_fetch.pending = None;
                self.mat_fetch.busy_label = None;
                match msg {
                    MaterialBgMsg::Fetched { result } => match result {
                        Ok((nm, toml_text, save_dir)) => match self.engine.add_material(&nm, &toml_text) {
                            Ok(_) => {
                                let mut note = format!("Added '{nm}' to project.");
                                if let Some(dir) = save_dir {
                                    let _ = std::fs::create_dir_all(&dir);
                                    let path = dir.join(format!("{nm}.toml"));
                                    if std::fs::write(&path, &toml_text).is_ok() {
                                        note = format!("Added '{nm}' to project and saved {}.", path.display());
                                    }
                                }
                                self.mat_fetch.is_error = false;
                                self.mat_fetch.message = Some(note);
                            }
                            Err(e) => {
                                self.mat_fetch.is_error = true;
                                self.mat_fetch.message = Some(e.to_string());
                            }
                        },
                        Err(e) => {
                            self.mat_fetch.is_error = true;
                            self.mat_fetch.message = Some(e);
                        }
                    },
                    MaterialBgMsg::Searched { result } => match result {
                        Ok(v) => {
                            let mut rows = v["catalog_matches"].as_array().cloned().unwrap_or_default();
                            if let Some(live) = v["live_matches"].as_array() {
                                rows.extend(live.iter().cloned());
                            }
                            self.mat_fetch.is_error = false;
                            self.mat_fetch.message = Some(format!("Found {} matching record(s).", rows.len()));
                            self.mat_fetch.search_results = rows;
                        }
                        Err(e) => {
                            self.mat_fetch.is_error = true;
                            self.mat_fetch.message = Some(e);
                        }
                    },
                }
            }
        }
        let is_busy = self.mat_fetch.pending.is_some();
        let catalog = nyx_engine::material_fetch::remote_catalog();
        let formats = ["rii", "aflow", "mp", "numbat"];
        let opt_trim = |s: &str| {
            let t = s.trim();
            (!t.is_empty()).then(|| t.to_string())
        };
        egui::CollapsingHeader::new("Download or import material").id_salt("mat_import_panel").default_open(false).show(ui, |ui| {
            ui.weak("Fetch additional entries from open databases on demand (HTTPS, non-blocking background fetch) or paste an offline JSON/YAML/JSONC record.");
            if let Some(busy) = &self.mat_fetch.busy_label {
                ui.weak(format!("Working in background: {busy}..."));
            }
            if !catalog.is_empty() {
                if self.mat_fetch.preset >= catalog.len() {
                    self.mat_fetch.preset = 0;
                }
                let sel = &catalog[self.mat_fetch.preset];
                ui.horizontal(|ui| {
                    ui.label("Preset:");
                    egui::ComboBox::from_id_salt("mat_fetch_preset").selected_text(format!("[{}] {}", sel.domain, sel.id)).show_ui(ui, |ui| {
                        for (idx, entry) in catalog.iter().enumerate() {
                            if ui.selectable_label(self.mat_fetch.preset == idx, format!("[{}] {} ({})", entry.domain, entry.id, entry.format)).clicked() {
                                self.mat_fetch.preset = idx;
                                self.mat_fetch.source = entry.id.to_string();
                                self.mat_fetch.target.clear();
                                self.mat_fetch.name = entry.id.to_string();
                                self.mat_fetch.extends = entry.extends.unwrap_or("").to_string();
                                self.mat_fetch.description = entry.summary.to_string();
                            }
                        }
                    });
                    if ui.add_enabled(!is_busy, egui::Button::new("Download preset")).clicked() {
                        let preset = catalog[self.mat_fetch.preset].clone();
                        let name_ov = opt_trim(&self.mat_fetch.name).unwrap_or_else(|| preset.id.to_string());
                        let ext_ov = opt_trim(&self.mat_fetch.extends).or_else(|| preset.extends.map(str::to_string));
                        let desc_ov = opt_trim(&self.mat_fetch.description).or_else(|| Some(preset.summary.to_string()));
                        let save_dir = self.mat_fetch.save_to_library.then(nyx_engine::material_fetch::default_materials_folder).flatten();
                        let api_key = opt_trim(&self.mat_fetch.mp_api_key);
                        let (tx, rx) = std::sync::mpsc::channel();
                        self.mat_fetch.pending = Some(rx);
                        self.mat_fetch.busy_label = Some(format!("downloading {}", preset.id));
                        std::thread::spawn(move || {
                            let res = (|| -> Result<(String, String, Option<std::path::PathBuf>), String> {
                                let url = nyx_engine::material_fetch::resolve_remote_url(preset.format, preset.target).map_err(|e| e.to_string())?;
                                let raw = nyx_engine::material_fetch::fetch_url_text_with_api_key(&url, api_key.as_deref()).map_err(|e| e.to_string())?;
                                let toml = nyx_engine::material_fetch::convert_to_toml(preset.format, &raw, &name_ov, desc_ov.as_deref(), ext_ov.as_deref(), Some(&url)).map_err(|e| e.to_string())?;
                                Ok((name_ov, toml, save_dir))
                            })();
                            let _ = tx.send(MaterialBgMsg::Fetched { result: res });
                        });
                    }
                });
                ui.weak(sel.summary);
            }

            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.label("Search catalog:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.search_query).hint_text("keyword or formula (e.g. Si, AlN, NbN)").desired_width(150.0));
                ui.checkbox(&mut self.mat_fetch.live_search, "live AFLOW/MP");
                if ui.add_enabled(!is_busy, egui::Button::new("Search")).clicked() {
                    let q = self.mat_fetch.search_query.trim().to_string();
                    let live = self.mat_fetch.live_search;
                    let key = opt_trim(&self.mat_fetch.mp_api_key);
                    let (tx, rx) = std::sync::mpsc::channel();
                    self.mat_fetch.pending = Some(rx);
                    self.mat_fetch.busy_label = Some(format!("searching '{q}'"));
                    std::thread::spawn(move || {
                        let eng = nyx_engine::Engine::new();
                        let res = eng.search_remote_materials(&q, None, None, live, key.as_deref()).map_err(|e| e.to_string());
                        let _ = tx.send(MaterialBgMsg::Searched { result: res });
                    });
                }
            });
            for row in self.mat_fetch.search_results.iter().take(6) {
                let id = row["id"].as_str().unwrap_or("");
                let fmt = row["format"].as_str().unwrap_or("");
                let tgt = row["target"].as_str().unwrap_or("");
                let sum = row["summary"].as_str().unwrap_or("");
                ui.horizontal(|ui| {
                    if ui.small_button("Use").clicked() {
                        self.mat_fetch.source = fmt.to_string();
                        self.mat_fetch.target = tgt.to_string();
                        self.mat_fetch.name = id.to_string();
                        self.mat_fetch.description = sum.to_string();
                    }
                    ui.weak(format!("[{fmt}] {id}: {sum}"));
                });
            }

            ui.add_space(4.0);
            egui::Grid::new("mat_fetch_fields").num_columns(2).show(ui, |ui| {
                ui.label("Source / format:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.source).hint_text("preset id, rii, aflow, mp, or numbat").desired_width(220.0));
                ui.end_row();
                ui.label("Entry / query:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.target).hint_text("e.g. main/Si/nk/Green-2008 or mp-149").desired_width(220.0));
                ui.end_row();
                ui.label("Material name:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.name).hint_text("optional slug").desired_width(220.0));
                ui.end_row();
                ui.label("Extends:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.extends).hint_text("optional base material (e.g. silicon)").desired_width(220.0));
                ui.end_row();
                ui.label("Description:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.description).hint_text("optional one-line comparison").desired_width(220.0));
                ui.end_row();
                ui.label("MP API key:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.mp_api_key).password(true).hint_text("optional Materials Project key (stdin curl)").desired_width(220.0));
                ui.end_row();
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.mat_fetch.save_to_library, "save to materials/ folder");
                if ui.add_enabled(!is_busy, egui::Button::new("Download custom")).clicked() {
                    let src = self.mat_fetch.source.trim().to_string();
                    let tgt = opt_trim(&self.mat_fetch.target);
                    let nm = opt_trim(&self.mat_fetch.name);
                    let ext = opt_trim(&self.mat_fetch.extends);
                    let desc = opt_trim(&self.mat_fetch.description);
                    let key = opt_trim(&self.mat_fetch.mp_api_key);
                    let save_dir = self.mat_fetch.save_to_library.then(nyx_engine::material_fetch::default_materials_folder).flatten();
                    let (tx, rx) = std::sync::mpsc::channel();
                    self.mat_fetch.pending = Some(rx);
                    self.mat_fetch.busy_label = Some(format!("downloading {src}"));
                    std::thread::spawn(move || {
                        let mut eng = nyx_engine::Engine::new();
                        let res = eng
                            .fetch_material_with_key(&src, tgt.as_deref(), nm.as_deref(), desc.as_deref(), ext.as_deref(), false, None, key.as_deref())
                            .map(|v| {
                                let final_nm = v["name"].as_str().unwrap_or("imported").to_string();
                                let toml = v["imported_toml"].as_str().or_else(|| v["toml"].as_str()).unwrap_or("").to_string();
                                (final_nm, toml, save_dir)
                            })
                            .map_err(|e| e.to_string());
                        let _ = tx.send(MaterialBgMsg::Fetched { result: res });
                    });
                }
            });

            ui.add_space(4.0);
            if self.mat_fetch.format >= formats.len() {
                self.mat_fetch.format = 0;
            }
            ui.horizontal(|ui| {
                ui.label("Offline import:");
                egui::ComboBox::from_id_salt("mat_import_fmt").selected_text(formats[self.mat_fetch.format]).show_ui(ui, |ui| {
                    for (idx, f) in formats.iter().enumerate() {
                        ui.selectable_value(&mut self.mat_fetch.format, idx, *f);
                    }
                });
                if ui.button("Import text").clicked() {
                    let nm = opt_trim(&self.mat_fetch.name).unwrap_or_else(|| "imported_material".into());
                    let ext = opt_trim(&self.mat_fetch.extends);
                    let desc = opt_trim(&self.mat_fetch.description);
                    let raw = self.mat_fetch.raw.clone();
                    let save_dir = self.mat_fetch.save_to_library.then(nyx_engine::material_fetch::default_materials_folder).flatten();
                    match self.engine.import_material(&nm, formats[self.mat_fetch.format], &raw, desc.as_deref(), ext.as_deref(), true, save_dir.as_deref()) {
                        Ok(v) => {
                            let added = v["name"].as_str().unwrap_or(&nm);
                            self.mat_fetch.is_error = false;
                            self.mat_fetch.message = Some(format!("Imported '{added}' into project."));
                        }
                        Err(e) => {
                            self.mat_fetch.is_error = true;
                            self.mat_fetch.message = Some(e.to_string());
                        }
                    }
                }
            });
            ui.add(egui::TextEdit::multiline(&mut self.mat_fetch.raw).desired_rows(4).hint_text("Paste YAML (rii), JSON (aflow/mp), or JSONC (numbat)...").font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
            if let Some(msg) = &self.mat_fetch.message {
                if self.mat_fetch.is_error {
                    ui.colored_label(crate::theme::danger(ui.ctx()), msg);
                } else {
                    ui.weak(msg);
                }
            }
        });
    }

    // ------------------------------------------------------------------ physics

    fn tab_physics(&mut self, ui: &mut Ui) {
        let ids: Vec<String> = self.engine.registry().ids().iter().map(|s| s.to_string()).collect();
        let in_project: Vec<String> = self.engine.project().physics.iter().map(|s| s.physics.clone()).collect();
        self.tab_heading(ui, "Physics");
        // in the window of the shell the physics are one dropdown of the toolbar (ADR 124); the classic window keeps the pills here
        let in_shell = self.shell.active();
        ui.horizontal_wrapped(|ui| {
            if in_shell {
                return;
            }
            let active = self.active();
            for id in &in_project {
                if ui.selectable_label(active.as_deref() == Some(id.as_str()), RichText::new(self.physics_name(id)).strong()).clicked() {
                    self.cmd("set_active", json!({ "physics": id }));
                    self.picking = None;
                    self.rules_dirty = true;
                }
            }
            let free: Vec<&String> = ids.iter().filter(|i| !in_project.contains(i)).collect();
            if !free.is_empty() && in_project.len() < 3 {
                egui::ComboBox::from_id_salt("add_physics").selected_text("+ add physics").show_ui(ui, |ui| {
                    for id in free {
                        if ui.selectable_label(false, self.physics_name(id)).clicked() {
                            self.cmd("add_physics", json!({ "physics": id }));
                            self.cmd("set_active", json!({ "physics": id }));
                            self.rules_dirty = true;
                        }
                    }
                });
            }
        });
        let Some(active) = self.active() else {
            ui.weak(if in_shell { "The project has no physics: add one from the Physics menu in the bar above." } else { "The project has no physics: add one above." });
            return;
        };
        let Some(physics) = self.engine.registry().get(&active) else { return };
        let (kinds, opts) = (physics.bc_kinds(), physics.params());
        let (exports, imports) = (physics.exports(), physics.imports());
        let Some(spec0) = self.engine.project().spec(&active).cloned() else { return };
        let mut spec = spec0.clone();
        ui.add_space(4.0);
        if !in_shell {
            ui.horizontal(|ui| {
                if ui.small_button("remove this physics from the project").clicked() {
                    self.cmd("remove_physics", json!({ "physics": active }));
                    self.rules_dirty = true;
                }
            });
        }
        // what physics exchange is said only when there is another physics to exchange with
        if in_project.len() > 1 && !imports.is_empty() {
            let from: Vec<String> = imports.iter().map(|i| format!("{} ({})", i.role, i.field)).collect();
            ui.weak(format!("Takes from other physics of the project: {}.", from.join(", ")));
        }
        if in_project.len() > 1 && !exports.is_empty() {
            ui.weak(format!("Offers to others: {}.", exports.iter().map(|e| e.role.clone()).collect::<Vec<_>>().join(", ")));
        }
        if !opts.is_empty() {
            ui.separator();
            ui.label(RichText::new("Options").strong());
            forms::param_form(ui, &format!("opt_{active}"), &opts, &mut spec.params);
        }
        ui.separator();
        ui.label(RichText::new("Boundary and domain conditions").strong());
        ui.weak("`on` is a selector over the geometry (`outer & plane(x = 0)`), a boundary or region name. 'Pick on mesh' writes the selector for you.");
        let mut remove = None;
        let mut pick: Option<Option<usize>> = None;
        for i in 0..spec.bcs.len() {
            let kind = kinds.iter().find(|k| k.kind == spec.bcs[i].kind).cloned();
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("{}", i + 1)).strong());
                    let shown = kind.as_ref().map(|k| k.name.clone()).unwrap_or_else(|| spec.bcs[i].kind.clone());
                    egui::ComboBox::from_id_salt(format!("bc_kind_{active}_{i}")).selected_text(shown).show_ui(ui, |ui| {
                        for k in &kinds {
                            if ui.selectable_label(spec.bcs[i].kind == k.kind, &k.name).clicked() && spec.bcs[i].kind != k.kind {
                                spec.bcs[i].kind = k.kind.clone();
                                spec.bcs[i].params = Value::Null;
                                forms::fill_required(&k.params, &mut spec.bcs[i].params);
                            }
                        }
                    });
                    if ui.small_button("✕").on_hover_text("remove this condition").clicked() {
                        remove = Some(i);
                    }
                });
                if let Some(k) = &kind {
                    if !k.doc.is_empty() {
                        ui.weak(&k.doc);
                    }
                }
                ui.horizontal(|ui| {
                    ui.label("on");
                    ui.add(egui::TextEdit::singleline(&mut spec.bcs[i].on).desired_width(220.0).hint_text("selector, boundary or region name"));
                    let is_boundary = kind.as_ref().is_none_or(|k| k.applies == Applies::Boundary);
                    if is_boundary && self.doc.mesh.is_some() {
                        let on = self.picking == Some(i);
                        if ui.selectable_label(on, "Pick on mesh").on_hover_text("Click boundaries in the viewport; the simplest selector that picks them is written here.").clicked() {
                            pick = Some(if on { None } else { Some(i) });
                        }
                    }
                });
                match self.rule_errors.get(&i) {
                    Some(e) => {
                        ui.colored_label(crate::theme::danger(ui.ctx()), e);
                    }
                    None => {
                        if let Some(l) = self.bc_labels.get(i).filter(|l| !l.is_empty()) {
                            ui.weak(format!("{} boundary number(s): {:?}", l.len(), l));
                        }
                    }
                }
                if let Some(k) = &kind {
                    forms::param_form(ui, &format!("bc_{active}_{i}"), &k.params, &mut spec.bcs[i].params);
                }
                if self.picking == Some(i) && !self.suggestions.is_empty() {
                    ui.weak("These selectors pick exactly what you clicked; choose one:");
                    for s in self.suggestions.clone() {
                        if ui.small_button(&s).clicked() {
                            spec.bcs[i].on = s;
                        }
                    }
                } else if self.picking == Some(i) && !self.picked.is_empty() {
                    ui.weak("No selector picks exactly this set of boundaries; type one, or click fewer.");
                }
            });
        }
        if let Some(i) = remove {
            spec.bcs.remove(i);
            if self.picking == Some(i) {
                self.picking = None;
            }
        }
        if let Some(p) = pick {
            self.picking = p;
            self.picked.clear();
            self.suggestions.clear();
        }
        let has_geom_pml = self
            .geom_edit
            .features
            .iter()
            .any(|f| !f.suppressed && matches!(f.op, nyx_engine::facade::geom::model::Op::Pml { .. }));
        let supports_pml_bc = kinds
            .iter()
            .any(|k| matches!(k.kind.as_str(), "pml" | "infinite_element"));
        let mut add_from_geom = false;
        if !kinds.is_empty() {
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt(format!("bc_add_{active}")).selected_text("+ add condition").show_ui(ui, |ui| {
                    for k in &kinds {
                        if ui.selectable_label(false, &k.name).clicked() {
                            let mut params = Value::Null;
                            forms::fill_required(&k.params, &mut params);
                            spec.bcs.push(nyx_engine::facade::physics::BcSpec { on: String::new(), kind: k.kind.clone(), params });
                        }
                    }
                });
                if has_geom_pml
                    && supports_pml_bc
                    && ui
                        .button("Add PML Conditions from Geometry")
                        .on_hover_text("Automatically create or update absorbing conditions for every PML domain in the geometry, jointly auto-calibrated from study wavelengths.")
                        .clicked()
                {
                    add_from_geom = true;
                }
            });
        }
        if spec != spec0 {
            self.set_problem(&active, &spec);
        }
        if add_from_geom && self.push_edits() {
            self.cmd(
                "add_pml_conditions",
                json!({ "physics": active, "auto_scale": true }),
            );
            self.rules_dirty = true;
        }
    }

    // ------------------------------------------------------------------ study

    fn tab_study(&mut self, ui: &mut Ui) {
        let Some(active) = self.active() else {
            ui.heading("Study");
            ui.weak("The project has no physics: add one in the Physics tab.");
            return;
        };
        let Some(physics) = self.engine.registry().get(&active) else { return };
        let studies = physics.studies();
        let pname = physics.name();
        let Some(spec0) = self.engine.project().spec(&active).cloned() else { return };
        let mut spec: ProblemSpec = spec0.clone();
        ui.heading(format!("Study: {pname}"));
        let dim = self.doc.dim();
        ui.horizontal_wrapped(|ui| {
            for s in &studies {
                if ui.selectable_label(spec.study.kind == s.kind, &s.name).on_hover_text(&s.doc).clicked() && spec.study.kind != s.kind {
                    spec.study.kind = s.kind.clone();
                    spec.study.params = Value::Null;
                    forms::fill_required(&s.params, &mut spec.study.params);
                }
            }
        });
        if let Some(s) = studies.iter().find(|s| s.kind == spec.study.kind) {
            ui.add(egui::Label::new(RichText::new(&s.doc).weak()).wrap());
            let runs_here = match s.dims {
                nyx_engine::facade::physics::Dims::Two => dim == 2,
                nyx_engine::facade::physics::Dims::Three => dim == 3,
                nyx_engine::facade::physics::Dims::Both => true,
            };
            if !runs_here {
                ui.colored_label(crate::theme::activity_text(ui.ctx()), format!("This study runs on a {}D mesh; the model on show is {}D.", if dim == 2 { 3 } else { 2 }, dim));
            }
            forms::param_form(ui, &format!("study_{active}_{}", s.kind), &s.params, &mut spec.study.params);
            if let Some(fd) = s.params.iter().find(|d| d.name == "freq" && d.unit == "Hz") {
                let freq = spec
                    .study
                    .params
                    .get("freq")
                    .and_then(Value::as_f64)
                    .unwrap_or_else(|| forms::default_value(fd).as_f64().unwrap_or(193.414e12));
                if freq > 0.0 && freq.is_finite() {
                    let lam0 = 299_792_458.0 / freq;
                    let mut lam = lam0;
                    ui.horizontal_wrapped(|ui| {
                        ui.label("wavelength λ₀ (vacuum) [m]")
                            .on_hover_text("Free-space wavelength c₀ / freq (e.g. type 1550 nm); updates Frequency [Hz] above.");
                        if num_id(ui, ("study_lambda0", &active, &s.kind), &mut lam).changed() && (lam - lam0).abs() > 0.0 && lam > 0.0 && lam.is_finite() {
                            if !spec.study.params.is_object() {
                                spec.study.params = Value::Object(serde_json::Map::new());
                            }
                            spec.study.params["freq"] = json!(299_792_458.0 / lam);
                        }
                        ui.weak(crate::fmt::si(lam, "m"));
                    });
                }
            }
        }
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label("element order");
            let mut o = spec.order;
            for v in [1u8, 2] {
                ui.selectable_value(&mut o, v, format!("P{v}"));
            }
            spec.order = o;
        });
        ui.weak("Options of the physics (plane stress or strain, solver, ...) are in the Physics tab.");
        ui.add_space(4.0);
        ui.separator();
        ui.label(RichText::new("Mesh convergence").strong());
        ui.weak("Solves again on a finer mesh and compares every eigenvalue and scalar. Do this before you trust a number.");
        let mut start = false;
        ui.horizontal(|ui| {
            ui.label("finer by");
            num_id(ui, "conv_f", &mut self.conv_factor);
            ui.label("tolerance");
            num_id(ui, "conv_tol", &mut self.conv_tol);
            start = ui.add_enabled(self.job_kind != Some(crate::app::JobKind::Open), egui::Button::new(if self.busy.is_some() { "Queue convergence check" } else { "Check convergence" })).clicked();
        });
        self.conv_factor = self.conv_factor.clamp(0.1, 0.95);
        if start {
            self.check_convergence();
        }
        if self.convergence_running() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("two solves are running (the toolbar shows the stage; Cancel stops it)...");
            });
        }
        if let Some(c) = &self.convergence {
            let ok = c["converged"] == true;
            ui.colored_label(if ok { crate::theme::success(ui.ctx()) } else { crate::theme::activity_text(ui.ctx()) }, format!("{} (worst relative difference {:.2e}; {} to {} cells, {:.1} s)", if ok { "converged" } else { "NOT converged" }, c["worst_relative"].as_f64().unwrap_or(f64::NAN), c["cells"]["coarse"], c["cells"]["refined"], c["seconds"].as_f64().unwrap_or(0.0)));
            egui::Grid::new("conv_rows").num_columns(4).striped(true).show(ui, |ui| {
                for r in c["rows"].as_array().into_iter().flatten().take(40) {
                    ui.label(r["quantity"].as_str().unwrap_or(""));
                    ui.monospace(format!("{:.6e}", r["coarse"].as_f64().unwrap_or(f64::NAN)));
                    ui.monospace(format!("{:.6e}", r["refined"].as_f64().unwrap_or(f64::NAN)));
                    let pass = r["converged"] == true;
                    ui.colored_label(if pass { crate::theme::success(ui.ctx()) } else { crate::theme::activity_text(ui.ctx()) }, format!("{:.1e}", r["relative"].as_f64().unwrap_or(f64::NAN)));
                    ui.end_row();
                }
            });
        }
        // what the check left the result on show as (it is the engine's state, so a check of another problem does not count)
        self.validated_line(ui);
        if spec != spec0 {
            self.set_problem(&active, &spec);
        }
        // in the window of the shell Compute is the one filled button of the bar (ADR 124)
        if self.shell.active() {
            return;
        }
        ui.add_space(8.0);
        let ok = self.can_queue_compute() && self.doc.mesh.is_some();
        if ui.add_enabled(ok, egui::Button::new(RichText::new("  Compute  ").strong().color(crate::theme::on_accent(ui.ctx()))).fill(crate::theme::accent(ui.ctx()))).clicked() {
            self.compute();
        }
    }

    // ------------------------------------------------------------------ results

    fn tab_results(&mut self, ui: &mut Ui) {
        let Some(mut cur) = self.current else {
            ui.weak("No results yet. Press Compute.");
            return;
        };
        self.tab_heading(ui, "Results");
        crate::spans::span("stale_band", || self.stale_band(ui));
        crate::spans::span("validated_line", || self.validated_line(ui));
        if self.results.len() > 1 {
            let stale: Vec<bool> = (0..self.results.len()).map(|i| self.stale_reason(i).is_some()).collect();
            let label = |i: usize, name: &str| if stale[i] { format!("{name} (stale)") } else { name.to_string() };
            egui::ComboBox::from_id_salt("result_sel").selected_text(label(cur, &self.results[cur].name)).show_ui(ui, |ui| {
                for (i, r) in self.results.iter().enumerate() {
                    ui.selectable_value(&mut cur, i, label(i, &r.name));
                }
            });
            if cur != self.current.unwrap() {
                self.current = Some(cur);
                self.show_results = true;
            }
        } else {
            ui.label(self.results[cur].name.clone());
        }
        let mut cmap = self.colormap;
        let mut fit = false;
        {
            let rv = &mut self.results[cur];
            let nsol = rv.sols.len();
            if nsol > 1 {
                ui.horizontal(|ui| {
                    if ui.small_button("<").clicked() && rv.sol > 0 {
                        rv.sol -= 1;
                    }
                    let mut s = rv.sol;
                    ui.spacing_mut().slider_width = (ui.available_width() - 230.0).clamp(60.0, 260.0);
                    // the rail takes the colour of an unpressed widget, which is the panel's own: give it the text colour, faded, so it can be seen
                    let rail = ui.visuals().text_color().gamma_multiply(0.35);
                    ui.visuals_mut().widgets.inactive.bg_fill = rail;
                    ui.add(egui::Slider::new(&mut s, 0..=nsol - 1).show_value(false));
                    rv.sol = s;
                    if ui.small_button(">").clicked() && rv.sol + 1 < nsol {
                        rv.sol += 1;
                    }
                    ui.label(rv.sols[rv.sol].label.clone());
                });
            } else {
                ui.label(rv.sols[0].label.clone());
            }
            let sol = rv.sol.min(nsol - 1);
            let names: Vec<String> = rv.sols[sol].fields.iter().map(|f| format!("{} [{}]", f.name, f.unit)).collect();
            if rv.field >= names.len() {
                rv.field = 0;
            }
            ui.horizontal(|ui| {
                ui.label("Field");
                egui::ComboBox::from_id_salt("field_sel").selected_text(names[rv.field].clone()).show_ui(ui, |ui| {
                    for (i, n) in names.iter().enumerate() {
                        ui.selectable_value(&mut rv.field, i, n);
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("Colours");
                egui::ComboBox::from_id_salt("cmap").selected_text(cmap.name()).show_ui(ui, |ui| {
                    for c in Colormap::ALL {
                        ui.selectable_value(&mut cmap, c, c.name());
                    }
                });
            });
            // range
            let auto = rv.manual_range.is_none();
            let mut auto_sel = auto;
            ui.horizontal(|ui| {
                ui.checkbox(&mut auto_sel, "auto range");
                if !auto_sel {
                    let (lo, hi) = rv.manual_range.unwrap_or_else(|| {
                        rv.manual_range = None;
                        rv.range()
                    });
                    let (mut a, mut b) = (lo as f64, hi as f64);
                    let c1 = num_id(ui, "manual_range_lo", &mut a).changed();
                    ui.label("to");
                    let c2 = num_id(ui, "manual_range_hi", &mut b).changed();
                    if c1 || c2 || rv.manual_range.is_none() {
                        rv.manual_range = Some((a as f32, b as f32));
                    }
                }
            });
            if auto_sel {
                rv.manual_range = None;
            } else if rv.manual_range.is_none() {
                let r = rv.range();
                rv.manual_range = Some(r);
            }
            if rv.sols.iter().any(|s| s.disp.is_some()) {
                ui.label("deformation scale");
                ui.horizontal(|ui| {
                    // the label is above so that the slider has the width of the panel less its box and two buttons; its rail is made visible as for the modes
                    ui.spacing_mut().slider_width = (ui.available_width() - 200.0).clamp(60.0, 200.0);
                    let rail = ui.visuals().text_color().gamma_multiply(0.35);
                    ui.visuals_mut().widgets.inactive.bg_fill = rail;
                    let max = (rv.auto_warp * 10.0).max(1.0);
                    ui.add(egui::Slider::new(&mut rv.warp, 0.0..=max).logarithmic(false).custom_formatter(|v, _| format!("{:.3e}", v)));
                    if ui.small_button("auto").clicked() {
                        rv.warp = rv.auto_warp;
                    }
                    if ui.small_button("0").clicked() {
                        rv.warp = 0.0;
                    }
                });
            }
            ui.checkbox(&mut self.show_edges, "mesh edges");
            if ui.button("Fit view").clicked() {
                fit = true;
            }
        }
        if fit {
            self.vp.fit_token += 1;
        }
        self.colormap = cmap;
        ui.separator();
        ui.label(RichText::new("Summary").strong());
        let rv = &self.results[cur];
        egui::Grid::new("scalars").num_columns(2).striped(true).show(ui, |ui| {
            for (n, v) in &rv.result.scalars {
                ui.label(n);
                match v {
                    nyx_engine::facade::physics::Scalar::Num(x) => ui.monospace(format!("{x:.6e}")),
                    nyx_engine::facade::physics::Scalar::Text(t) => ui.monospace(t.clone()),
                };
                ui.end_row();
            }
        });
        ui.weak(format!("{} unknowns, solved in {:.2} s", rv.result.dofs, rv.result.seconds));
        crate::spans::span("solver_log", || {
            // closed until opened: drawing eleven short monospace lines for the first time cost 14 ms in the frame after a result arrives (P3
            // frame log; neither one label for the whole log nor laying the glyphs out beforehand made it cheaper), and the header says how many
            egui::CollapsingHeader::new(format!("What the solver did ({} lines)", rv.result.log.len())).default_open(false).show(ui, |ui| {
                egui::ScrollArea::vertical().max_height(260.0).id_salt("solver_log").show(ui, |ui| {
                    for l in &rv.result.log {
                        ui.add(egui::Label::new(RichText::new(l).monospace().weak().size(11.0)).wrap());
                    }
                });
            });
        });
        ui.add_space(6.0);
        crate::spans::span("reductions", || self.reductions(ui, cur));
        ui.add_space(6.0);
        if ui.button("Export CSV...").clicked() {
            self.export_csv();
        }
    }

    /// What the Reduce form offers for the result on show: its physics, the solution on show, the fields it can read and the regions of its mesh.
    fn reduce_choices(&self, cur: usize) -> Option<(String, usize, Vec<String>, Vec<String>)> {
        let physics = self.results[cur].physics.clone()?;
        let r = &self.results[cur];
        let sol = r.sol.min(r.result.solutions.len().checked_sub(1)?);
        let s = &r.result.solutions[sol];
        let mut names: Vec<String> = s.cell_fields.iter().map(|f| f.name.clone()).collect();
        for f in &s.fields {
            if f.components.is_empty() && !names.contains(&f.name) {
                names.push(f.name.clone());
            }
        }
        if names.is_empty() {
            return None;
        }
        let regions = self.engine.result_regions(&physics);
        Some((physics, sol, names, regions))
    }

    /// The Reduce form: integrals, extremes and the fraction in a region of the field on show, by the engine commands `integrate`, `extremum` and
    /// `confinement` (ADR 90), so the number is the one a script or the MCP server gets. The answers stay listed, newest first.
    fn reductions(&mut self, ui: &mut Ui, cur: usize) {
        let Some((_, _, names, regions)) = self.reduce_choices(cur) else { return };
        let f = &mut self.reduce;
        f.field = f.field.min(names.len() - 1);
        f.region = f.region.min(regions.len());
        let mut action: Option<&'static str> = None;
        egui::CollapsingHeader::new("Reduce: integral, extremes, fraction").default_open(false).show(ui, |ui| {
            ui.weak("Computed by the engine on the cells the result was solved on, with the rule in ADR 90: the same numbers a script or the MCP server gets.");
            ui.horizontal(|ui| {
                ui.label("Field");
                egui::ComboBox::from_id_salt("reduce_field").selected_text(names[f.field].clone()).show_ui(ui, |ui| {
                    for (i, n) in names.iter().enumerate() {
                        ui.selectable_value(&mut f.field, i, n);
                    }
                });
                egui::ComboBox::from_id_salt("reduce_quantity").selected_text(QUANTITIES[f.quantity.min(2)]).show_ui(ui, |ui| {
                    for (i, n) in QUANTITIES.iter().enumerate() {
                        ui.selectable_value(&mut f.quantity, i, *n);
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("Region");
                egui::ComboBox::from_id_salt("reduce_region").selected_text(if f.region == 0 { "whole mesh".to_string() } else { regions[f.region - 1].clone() }).show_ui(ui, |ui| {
                    ui.selectable_value(&mut f.region, 0, "whole mesh");
                    for (i, r) in regions.iter().enumerate() {
                        ui.selectable_value(&mut f.region, i + 1, r);
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("Formula");
                ui.add(egui::TextEdit::singleline(&mut f.formula).desired_width(ui.available_width().min(320.0)).hint_text("e.g. eps_r * norm2(E), or where(z > 0, abs2(Ex), 0)"))
                    .on_hover_text("When this is not empty, Integral, Largest and Smallest read the formula instead of the field and quantity above. Names: the fields of the result, x y z, material values (eps_r, rho), region(\"Si\"), where(c, a, b), vectors (E: dot, cross, norm, norm2). The same language as the `expr` of a script or the MCP server.");
            });
            ui.horizontal(|ui| {
                if ui.button("Integral").on_hover_text("the integral over the region (a complex field with the value quantity gives re + im i); of the formula when there is one").clicked() {
                    action = Some("integrate");
                }
                if ui.button("Largest").on_hover_text("the largest real value at the nodes of the region").clicked() {
                    action = Some("max");
                }
                if ui.button("Smallest").clicked() {
                    action = Some("min");
                }
                if ui.add_enabled(f.region > 0, egui::Button::new("Fraction in region")).on_hover_text("the integral of the squared moduli of the components of the vector this field belongs to (Ex, Ey, Ez for Ex) over the region, divided by the same over the whole mesh; no material weight").clicked() {
                    action = Some("fraction");
                }
            });
            for l in f.lines.iter().take(8) {
                ui.add(egui::Label::new(RichText::new(l).monospace().size(11.0)).wrap());
            }
        });
        if let Some(a) = action {
            self.reduce_now(a);
        }
    }

    /// Runs one reduction of the Reduce form (`integrate`, `max`, `min` or `fraction`) on the field, quantity and region it holds, through the engine
    /// command of that name, and puts the answer first in its list.
    pub fn reduce_now(&mut self, action: &str) {
        let Some(cur) = self.current else { return };
        let Some((physics, sol, names, regions)) = self.reduce_choices(cur) else { return };
        let name = names[self.reduce.field.min(names.len() - 1)].clone();
        let quantity = self.reduce.quantity.min(2);
        let region = (self.reduce.region > 0 && self.reduce.region <= regions.len()).then(|| regions[self.reduce.region - 1].clone());
        let formula = self.reduce.formula.trim().to_string();
        let expr = if !formula.is_empty() {
            json!(formula)
        } else {
            match quantity {
                1 => json!({ "op": "abs", "args": [name] }),
                2 => json!({ "op": "abs2", "args": [name] }),
                _ => json!(name),
            }
        };
        let label = if !formula.is_empty() {
            format!("`{formula}`")
        } else {
            match quantity {
                1 => format!("|{name}|"),
                2 => format!("|{name}|²"),
                _ => name.clone(),
            }
        };
        let place = region.clone().unwrap_or_else(|| "the whole mesh".into());
        let mut args = json!({ "physics": physics, "solution": sol });
        if let Some(r) = &region {
            args["region"] = json!(r);
        }
        let stale = |v: &Value| if v["stale"] == json!(true) { "  (stale result)" } else { "" };
        let line = match action {
            "integrate" => {
                args["expr"] = expr;
                self.cmd("integrate", args).map(|v| format!("integral of {label} over {place} = {}{}", reduce_number(&v["value"]), stale(&v)))
            }
            "fraction" => {
                if region.is_none() {
                    self.set_error("fraction: choose a region".to_string());
                    return;
                }
                // the components of the vector this field belongs to: the engine knows them (`field` on the name without its axis letter)
                let stem = name.strip_suffix(['x', 'y', 'z']).unwrap_or(&name).to_string();
                let comps: Vec<String> = self
                    .engine
                    .call("field", &json!({ "physics": physics, "solution": sol, "name": stem }))
                    .ok()
                    .and_then(|v| v["components"].as_object().map(|m| m.values().filter_map(|c| c["name"].as_str().map(String::from)).collect()))
                    .unwrap_or_else(|| vec![name.clone()]);
                args["fields"] = json!(comps);
                let shown = comps.join(", ");
                self.cmd("confinement", args).map(|v| format!("share of |{shown}|² in {place} = {:.6}{}", v["fraction"].as_f64().unwrap_or(f64::NAN), stale(&v)))
            }
            op => {
                args["expr"] = expr;
                args["op"] = json!(op);
                self.cmd("extremum", args).map(|v| {
                    let at: Vec<String> = v["point"].as_array().map(|p| p.iter().map(|x| format!("{:.4e}", x.as_f64().unwrap_or(f64::NAN))).collect()).unwrap_or_default();
                    format!("{op} of {label} in {place} = {:.6e} at ({})", v["value"].as_f64().unwrap_or(f64::NAN), at.join(", "))
                })
            }
        };
        if let Some(l) = line {
            self.reduce.lines.insert(0, l);
            self.reduce.lines.truncate(20);
        }
    }

    // ------------------------------------------------------------------ analysis

    /// Whether the result on show carries a polarization analysis.
    pub fn has_analysis(&self) -> bool {
        self.current.is_some_and(|i| self.results[i].analysis.is_some())
    }

    /// The analysis section: Stokes parameters of every optical mode, the Poincare sphere and the table of `nyx-analysis`.
    fn tab_analysis(&mut self, ui: &mut Ui) {
        use nyx_engine::analysis::polarization::{mode_markers, sphere_drawing, SphereState, MODE_COLOURS};
        let Some(cur) = self.current else {
            ui.weak("No results yet. Press Compute.");
            return;
        };
        ui.heading("Analysis");
        let mut export: Option<&'static str> = None;
        {
            let rv = &mut self.results[cur];
            let selectable = rv.physics.as_deref() == Some("ems");
            let Some(a) = rv.analysis.as_ref() else {
                ui.weak("This result has no polarization analysis. Mode studies of the electromagnetic waves physics have one.");
                return;
            };
            ui.label(RichText::new("Polarization of the modes").strong());
            ui.weak("Stokes parameters from the complex field of each mode, summed over the cross-section. s1 > 0 is x (lateral, quasi-TE), s1 < 0 is z (normal, quasi-TM), s2 the diagonals, s3 the sense of rotation (+ counter-clockwise). P is the degree of polarization: below 1 where the polarization changes across the mode.");
            // the picture
            let side = ui.available_width().clamp(220.0, 360.0);
            let (resp, painter) = ui.allocate_painter(egui::vec2(side, side), egui::Sense::hover());
            let dark = ui.visuals().dark_mode;
            let ink = theme::sphere_ink(dark);
            let with_alpha = |c: Color32, o: f64| Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), (o.clamp(0.0, 1.0) * 255.0) as u8);
            let hex = |s: &str| {
                let v = u32::from_str_radix(s.trim_start_matches('#'), 16).unwrap_or(0);
                Color32::from_rgb((v >> 16) as u8, (v >> 8) as u8, v as u8)
            };
            let d = sphere_drawing(side as f64);
            let origin = resp.rect.min;
            let pt = |x: f64, y: f64| origin + egui::vec2(x as f32, y as f32);
            let centre = pt(d.centre.0, d.centre.1);
            let states: Vec<(usize, SphereState)> = a.rows.iter().enumerate().filter_map(|(i, r)| r.report.as_ref().map(|m| (i, SphereState { label: (i + 1).to_string(), s: m.sphere, dashed: r.group.is_some() }))).collect();
            let markers = mode_markers(&states.iter().map(|(_, s)| s.clone()).collect::<Vec<_>>(), side as f64);
            let selected = if selectable { Some(rv.sol) } else { None };
            let draw_markers = |front: bool| {
                for (m, (row, _)) in markers.iter().zip(&states) {
                    if m.front != front {
                        continue;
                    }
                    let col = with_alpha(hex(MODE_COLOURS[m.colour]), if front { 1.0 } else { 0.6 });
                    let (tip, c0) = (pt(m.tip.0, m.tip.1), pt(m.centre.0, m.centre.1));
                    if (tip - c0).length() >= 1.0 {
                        let stroke = egui::Stroke::new(1.8, col);
                        if m.dashed {
                            painter.extend(egui::Shape::dashed_line(&[c0, tip], stroke, 4.0, 3.0));
                        } else {
                            painter.line_segment([c0, tip], stroke);
                        }
                    }
                    let big = selected == Some(*row);
                    painter.circle(tip, side * if big { 0.040 } else { 0.030 }, col, egui::Stroke::new(if big { 2.0 } else { 0.8 }, theme::on_marker()));
                    painter.text(tip, egui::Align2::CENTER_CENTER, &m.label, egui::FontId::proportional(side * 0.036), theme::on_marker());
                }
            };
            let line = |l: &nyx_engine::analysis::polarization::WireLine| -> Vec<egui::Shape> {
                let pts: Vec<egui::Pos2> = l.pts.iter().map(|p| pt(p.0, p.1)).collect();
                if l.front {
                    vec![egui::Shape::line(pts, egui::Stroke::new(l.width as f32, with_alpha(ink, l.opacity)))]
                } else {
                    egui::Shape::dashed_line(&pts, egui::Stroke::new(l.width as f32, with_alpha(ink, l.opacity * 0.6)), 2.0, 3.0)
                }
            };
            // far side, the translucent ball, near side, labels
            for l in d.lines.iter().filter(|l| !l.front) {
                painter.extend(line(l));
            }
            draw_markers(false);
            painter.circle(centre, d.radius as f32, with_alpha(theme::sphere_body(), if dark { 0.10 } else { 0.22 }), egui::Stroke::new(1.0, with_alpha(ink, 0.6)));
            for l in d.lines.iter().filter(|l| l.front) {
                painter.extend(line(l));
            }
            draw_markers(true);
            for lab in &d.labels {
                painter.text(pt(lab.x, lab.y), egui::Align2::CENTER_CENTER, &lab.text, egui::FontId::proportional(side * 0.05), with_alpha(ui.visuals().text_color(), if lab.front { 0.9 } else { 0.5 }));
            }
            // the modes
            ui.add_space(4.0);
            let mut pick = None;
            for (i, r) in a.rows.iter().enumerate() {
                let colour = states.iter().position(|(k, _)| *k == i).map(|p| hex(MODE_COLOURS[p % MODE_COLOURS.len()])).unwrap_or_else(|| ui.visuals().weak_text_color());
                let text = match &r.report {
                    Some(m) => format!("{:>2}  {} {:.5}  {}  P {:.2}  TE {:.2}{}", i + 1, a.value_label, r.neff.0, m.class.label(), m.dop, m.te_fraction, r.group.map_or(String::new(), |g| format!("  G{}", g + 1))),
                    None => format!("{:>2}  {} {:.5}  no intensity", i + 1, a.value_label, r.neff.0),
                };
                ui.horizontal(|ui| {
                    ui.label(RichText::new("\u{25cf}").color(colour));
                    let on = selected == Some(i);
                    let resp = ui.selectable_label(on, RichText::new(text).monospace().size(11.5));
                    if selectable && resp.clicked() {
                        pick = Some(i);
                    }
                });
            }
            if let Some(i) = pick {
                if i < rv.sols.len() {
                    rv.sol = i;
                    rv.values_cache = None;
                }
            }
            for (gi, g) in a.groups.iter().enumerate() {
                let names: Vec<String> = g.members.iter().map(|i| (i + 1).to_string()).collect();
                egui::CollapsingHeader::new(format!("G{}: modes {} are degenerate", gi + 1, names.join(", "))).default_open(true).show(ui, |ui| {
                    ui.weak("The eigensolver returns some mixture of degenerate modes, so their own points on the sphere are not unique. The group state and the resolved pair are.");
                    if let Some(m) = &g.invariant {
                        ui.label(RichText::new(format!("group: {}", m.summary_line())).monospace().size(11.0));
                    }
                    if let Some(r) = &g.resolved {
                        for (k, rep) in r.reports.iter().enumerate() {
                            if let Some(m) = rep {
                                ui.label(RichText::new(format!("resolved {}: {}", k + 1, m.summary_line())).monospace().size(11.0));
                            }
                        }
                    }
                    if let Some(n) = &g.note {
                        ui.weak(n.clone());
                    }
                });
            }
            egui::CollapsingHeader::new("Table").default_open(false).show(ui, |ui| {
                egui::ScrollArea::horizontal().id_salt("pol_table").show(ui, |ui| {
                    ui.add(egui::Label::new(RichText::new(a.table()).monospace().size(10.5)).extend());
                });
            });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("Export CSV...").clicked() {
                    export = Some("csv");
                }
                if ui.button("Export sphere (SVG)...").clicked() {
                    export = Some("svg");
                }
                if ui.button("Copy table").clicked() {
                    ui.ctx().copy_text(a.table());
                }
            });
        }
        if let Some(kind) = export {
            self.export_analysis(kind);
        }
    }

    /// Bottom panel: how the shown field's extremes (or the mode frequencies) vary over the
    /// solutions of the current result.
    pub fn results_panel(&mut self, ui: &mut Ui) {
        let Some(cur) = self.current else { return };
        self.stale_band(ui);
        let rv = &self.results[cur];
        let nsol = rv.sols.len();
        ui.horizontal(|ui| {
            ui.label(RichText::new(rv.name.clone()).strong());
            ui.weak(if nsol > 1 { format!("{nsol} solutions · click the plot to show one") } else { "1 solution".to_string() });
        });
        if nsol < 2 {
            ui.weak("Single solution: use the viewport probe (hover) for values.");
            return;
        }
        let fname = rv.sols[rv.sol.min(nsol - 1)].fields.get(rv.field).map(|f| f.name.clone()).unwrap_or_default();
        let xs: Vec<f64> = rv.sols.iter().enumerate().map(|(i, s)| s.value.unwrap_or(i as f64 + 1.0)).collect();
        let is_modes = rv.sols[0].neff.is_some();
        // eigenfrequencies: solids and the 3D electromagnetic eigenmodes (frequency, no effective index)
        let is_eigen_plot = (rv.physics.as_deref() == Some("sms") || (rv.physics.as_deref() == Some("ems") && rv.sols[0].neff.is_none())) && rv.sols[0].value.is_some();
        let xlabel = if is_eigen_plot || is_modes { "mode number" } else if rv.sols[0].value.is_some() { "time (s)" } else { "solution" };
        let mut max_pts = vec![];
        let mut min_pts = vec![];
        for (i, s) in rv.sols.iter().enumerate() {
            if let Some(f) = s.fields.iter().find(|f| f.name == fname) {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for &v in &f.values {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
                max_pts.push([xs[i], hi]);
                min_pts.push([xs[i], lo]);
            }
        }
        let current = xs[rv.sol.min(nsol - 1)];
        let is_eigen = is_eigen_plot;
        // where each solution sits on the x axis: the mode number for modes, else its own value
        let at: Vec<f64> = if is_eigen || is_modes { (0..nsol).map(|i| i as f64 + 1.0).collect() } else { xs.clone() };
        let chosen = rv.sol.min(nsol - 1);
        let shown_colour = ui.visuals().warn_fg_color;
        // frequencies are drawn in the unit of their size (10 GHz, not 10000000000)
        let top = if is_eigen { xs.iter().fold(0.0f64, |m, x| m.max(x.abs())) } else { 0.0 };
        let (fscale, funit) = if top >= 1e9 { (1e9, "GHz") } else if top >= 1e6 { (1e6, "MHz") } else if top >= 1e3 { (1e3, "kHz") } else { (1.0, "Hz") };
        let ylabel = if is_eigen { format!("frequency ({funit})") } else if is_modes { "|neff|".to_string() } else { fname.clone() };
        let mut plot = egui_plot::Plot::new("summary_plot")
            .set_margin_fraction(egui::vec2(0.08, 0.3))
            .height(ui.available_height().max(90.0))
            .x_axis_label(xlabel)
            .y_axis_label(ylabel)
            .y_axis_formatter(move |m, _| if is_eigen { format!("{}", m.value / fscale) } else { format!("{}", m.value) })
            .legend(egui_plot::Legend::default());
        if is_eigen || is_modes {
            plot = plot
                .x_grid_spacer(integer_mode_marks)
                .x_axis_formatter(format_integer_mode_tick);
        }
        let out = plot.show(ui, |p| {
            if is_eigen {
                p.points(egui_plot::Points::new("f", egui_plot::PlotPoints::from(xs.iter().enumerate().map(|(i, _)| [i as f64 + 1.0, xs[i]]).collect::<Vec<_>>())).radius(4.0));
                p.points(egui_plot::Points::new("", egui_plot::PlotPoints::from(vec![[chosen as f64 + 1.0, xs[chosen]]])).radius(7.0).color(shown_colour));
            } else if is_modes {
                // neff per mode
                let pts: Vec<[f64; 2]> = rv.sols.iter().enumerate().map(|(i, s)| [i as f64 + 1.0, s.neff.map(|n| n.0.hypot(n.1)).unwrap_or(0.0)]).collect();
                let shown = pts[chosen];
                p.points(egui_plot::Points::new("|neff|", egui_plot::PlotPoints::from(pts)).radius(4.0));
                p.points(egui_plot::Points::new("", egui_plot::PlotPoints::from(vec![shown])).radius(7.0).color(shown_colour));
            } else {
                p.line(egui_plot::Line::new(format!("max {fname}"), egui_plot::PlotPoints::from(max_pts.clone())));
                p.line(egui_plot::Line::new(format!("min {fname}"), egui_plot::PlotPoints::from(min_pts.clone())));
                p.vline(egui_plot::VLine::new("now", current).color(shown_colour));
            }
        });
        // a click on the plot shows the solution nearest to it (the same as the slider and the arrows)
        if out.response.clicked() {
            if let Some(pos) = out.response.interact_pointer_pos() {
                let x = out.transform.value_from_position(pos).x;
                if let Some(i) = nearest_solution(&at, x) {
                    let rv = &mut self.results[cur];
                    if i != rv.sol {
                        rv.sol = i;
                        rv.values_cache = None;
                    }
                }
            }
        }
    }
}

/// The index of the solution drawn closest (along x) to `x`; `None` when there is none or `x` is not a number.
pub fn nearest_solution(at: &[f64], x: f64) -> Option<usize> {
    if x.is_nan() {
        return None;
    }
    at.iter().enumerate().filter(|(_, a)| a.is_finite()).min_by(|a, b| (a.1 - x).abs().total_cmp(&(b.1 - x).abs())).map(|(i, _)| i)
}

/// Grid ticks for discrete 1-based mode indices: integer values >= 1 with nice integer steps (no fractional `1.1`, `1.2`, ...).
pub fn integer_mode_marks(input: egui_plot::GridInput) -> Vec<egui_plot::GridMark> {
    let (lo, hi) = input.bounds;
    if !lo.is_finite() || !hi.is_finite() || hi < 1.0 {
        return vec![];
    }
    let raw = input.base_step_size.abs().max(1.0);
    let exp = 10f64.powf(raw.log10().floor());
    let norm = raw / exp;
    let nice = if norm <= 1.0 {
        1.0
    } else if norm <= 2.0 {
        2.0
    } else if norm <= 5.0 {
        5.0
    } else {
        10.0
    };
    let step = (nice * exp).round().max(1.0);
    let first = ((lo / step).ceil() * step).max(step);
    let mut marks = Vec::new();
    let mut x = first;
    for _ in 0..512 {
        if x > hi + 1e-9 {
            break;
        }
        if x >= lo - 1e-9 && x >= 1.0 {
            marks.push(egui_plot::GridMark { value: x, step_size: step });
        }
        x += step;
    }
    marks
}

/// Formats discrete mode-number ticks as integers and hides any non-integer marks.
pub fn format_integer_mode_tick(m: egui_plot::GridMark, _range: &std::ops::RangeInclusive<f64>) -> String {
    let rounded = m.value.round();
    if rounded >= 1.0 && (m.value - rounded).abs() < 1e-6 {
        format!("{rounded:.0}")
    } else {
        String::new()
    }
}

/// What the Reduce form can take of a field.
const QUANTITIES: [&str; 3] = ["value", "modulus", "modulus squared"];

/// A number of a reduction answer as the form shows it: a real one, or `re + im i` for a pair.
pub fn reduce_number(v: &Value) -> String {
    match v {
        Value::Array(a) if a.len() == 2 => format!("{:.6e} + {:.6e} i", a[0].as_f64().unwrap_or(f64::NAN), a[1].as_f64().unwrap_or(f64::NAN)),
        other => format!("{:.6e}", other.as_f64().unwrap_or(f64::NAN)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_click_on_the_mode_plot_picks_the_nearest_solution() {
        let at = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(nearest_solution(&at, 0.2), Some(0));
        assert_eq!(nearest_solution(&at, 2.4), Some(1));
        assert_eq!(nearest_solution(&at, 2.6), Some(2));
        assert_eq!(nearest_solution(&at, 99.0), Some(3));
        // times, not mode numbers
        assert_eq!(nearest_solution(&[0.0, 0.5, 5.0], 1.0), Some(1));
        assert_eq!(nearest_solution(&[], 1.0), None);
        assert_eq!(nearest_solution(&at, f64::NAN), None);
    }

    #[test]
    fn mode_plot_ticks_are_positive_integers_only() {
        let marks = integer_mode_marks(egui_plot::GridInput { bounds: (0.8, 2.2), base_step_size: 0.1 });
        let vals: Vec<f64> = marks.iter().map(|m| m.value).collect();
        assert_eq!(vals, vec![1.0, 2.0]);
        assert_eq!(format_integer_mode_tick(marks[0], &(0.8..=2.2)), "1");
        assert_eq!(format_integer_mode_tick(marks[1], &(0.8..=2.2)), "2");
        assert_eq!(
            format_integer_mode_tick(egui_plot::GridMark { value: 1.5, step_size: 0.5 }, &(0.8..=2.2)),
            ""
        );
        assert_eq!(
            format_integer_mode_tick(egui_plot::GridMark { value: 0.0, step_size: 1.0 }, &(0.0..=3.0)),
            ""
        );

        let wide = integer_mode_marks(egui_plot::GridInput { bounds: (0.5, 25.0), base_step_size: 2.3 });
        assert!(!wide.is_empty());
        assert!(wide.iter().all(|m| m.value >= 1.0 && (m.value - m.value.round()).abs() < 1e-9));
        assert!(wide.iter().all(|m| !format_integer_mode_tick(*m, &(0.5..=25.0)).is_empty()));
    }

    #[test]
    fn parse_si_f64_handles_units_and_rejects_expressions() {
        assert_eq!(parse_si_f64("220 nm"), Some(2.2e-7));
        assert_eq!(parse_si_f64("500nm"), Some(5e-7));
        assert_eq!(parse_si_f64("1.55 µm"), Some(1.55e-6));
        assert_eq!(parse_si_f64("1.55 um"), Some(1.55e-6));
        assert_eq!(parse_si_f64("193.414 THz"), Some(193.414e12));
        assert_eq!(parse_si_f64("2.85"), Some(2.85));
        assert_eq!(parse_si_f64("2,85"), Some(2.85));
        assert_eq!(parse_si_f64("100 GPa"), Some(100e9));
        assert_eq!(parse_si_f64("-w/2"), None);
        assert_eq!(parse_si_f64("NaN"), None);
        assert_eq!(parse_si_f64("inf"), None);
    }

    #[test]
    fn geomui_parameter_boxes_have_readable_widths_and_num_id_accepts_si_text() {
        let ctx = egui::Context::default();
        let mut m = Model::new_2d();
        m.param("w", 5.0e-7);
        m.param("h_si", 2.2e-7);
        let r = m.rect("core", ("-w/2", 0.0), ("w/2", "h_si"));
        m.overlay(&[("Si", r)]);
        let mut list = m.list;
        let mut st = crate::geomui::State::default();

        let mut num_box_id = egui::Id::NULL;
        let mut num_rect = egui::Rect::NOTHING;
        let mut auto_num_rect = egui::Rect::NOTHING;
        let mut v = 0.05_f64;
        let mut v_auto = 1.0_f64;

        for _ in 0..2 {
            let _ = ctx.run_ui(egui::RawInput::default(), |ui| {
                egui::CentralPanel::default().show(ui, |ui| {
                    ui.set_max_width(350.0);
                    crate::geomui::show(ui, &mut st, &mut list);
                    let resp = num_id(ui, "test_si_input", &mut v);
                    num_box_id = resp.id;
                    num_rect = resp.rect;
                    auto_num_rect = num(ui, &mut v_auto).rect;
                });
            });
        }
        assert!(num_rect.width() >= 104.0, "num_id width={}", num_rect.width());
        assert!(auto_num_rect.width() >= 104.0, "num width={}", auto_num_rect.width());
        let focusable: Vec<egui::Rect> = ctx.viewport(|vp| {
            vp.prev_pass
                .widgets
                .get_layer(egui::LayerId::background())
                .filter(|w| w.enabled && w.sense.is_focusable())
                .map(|w| w.rect)
                .collect()
        });
        assert!(focusable.len() >= 2, "expected param_name and param_val focusable widgets");
        let name_rect = focusable[0];
        let val_rect = focusable[1];
        assert!(name_rect.width() >= 72.0, "param_name width={}", name_rect.width());
        assert!(val_rect.width() >= 104.0, "param_val width={}", val_rect.width());

        // Focus num_id and simulate typing "220 nm"
        ctx.memory_mut(|mem| mem.request_focus(num_box_id));
        let mut raw = egui::RawInput::default();
        raw.events.push(egui::Event::Text("220 nm".into()));
        let mut changed = false;
        let _ = ctx.run_ui(raw, |ui| {
            egui::CentralPanel::default().show(ui, |ui| {
                changed = num_id(ui, "test_si_input", &mut v).changed();
            });
        });
        assert!(changed, "typing '220 nm' into focused num_id should mark changed");
        assert!((v - 2.2e-7).abs() < 1e-18, "v={v}");
    }
}
