// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Inspector tabs (geometry, materials, physics, study, results, analysis) and the results panel. The Materials, Physics and Study tabs are
//! drawn from the descriptors of the physics (`forms`), so they are the same for every physics; every change is handed to the engine as a
//! command (`set_problem`, `set_mesh`, `add_condition`, ...).

use crate::app::App;
use crate::colormap::Colormap;
use crate::forms;
use crate::state::*;
use eframe::egui::{self, Color32, DragValue, RichText, Ui};
use nyx_engine::facade::geom::model::{FeatureList, Model};
use nyx_engine::facade::physics::{Applies, MaterialAssign, ProblemSpec};
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

/// A number field: drag to change, or double-click and type (scientific notation accepted).
pub fn num(ui: &mut Ui, v: &mut f64) -> egui::Response {
    let speed = if *v == 0.0 { 0.01 } else { v.abs() * 0.005 };
    ui.add(DragValue::new(v).speed(speed).custom_formatter(|n, _| fmt_num(n)).custom_parser(|s| s.trim().replace(',', ".").parse::<f64>().ok()))
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
        *sweep = if on { Some(vec![SweepLayer { thickness: 0.1, n: 2, domain_offset: 0, fill: None }]) } else { None };
    }
    if let Some(layers) = sweep {
        let mut rm = None;
        for (i, l) in layers.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("layer {}", i + 1));
                ui.label("thickness");
                num(ui, &mut l.thickness);
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
            layers.push(SweepLayer { thickness: 0.1, n: 2, domain_offset: 0, fill: None });
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
        ui.heading("Geometry");
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
                m.sweep = Some(vec![SweepLayer { thickness: 1.0, n: 8, domain_offset: 0, fill: None }]);
            }
        });
        ui.horizontal(|ui| {
            ui.label("element size h (m)");
            num(ui, &mut m.h);
            ui.weak(format!("{}D model, x-z cross-section; the size rules of the model refine it", g.dim));
        });
        ui.checkbox(&mut m.physics_sizing, "size the mesh from the physics")
            .on_hover_text("Elements are made fine enough for the wavelength of the study in each material (points per wavelength by element order). Applies to EM waves and solid-mechanics eigenfrequencies; the sizes of the model still apply.");
        if m.physics_sizing {
            let mut ppw = m.points_per_wavelength.unwrap_or(0.0);
            ui.horizontal(|ui| {
                ui.label("points per wavelength (0 = by element order)");
                if num(ui, &mut ppw).changed() {
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
                num(ui, s);
            }
        });
        ui.checkbox(&mut m.structured_pml, "structured absorbing layers").on_hover_text("Rebuild the absorbing layers of the model as structured layers (quadrilaterals, prisms, hexahedra) cut into simplices.");
        sweep_editor(ui, &mut m.sweep);
        // as you type: the problem of the edits, or their size and cost from a trial mesh (nothing is built: `feedback.rs`); right under the
        // controls it answers, not below the feature list (the lines show the state of the frame before, which is the edit being typed)
        self.feedback_lines(ui);
        ui.add_space(4.0);
        crate::geomui::show(ui, &mut self.geom_ui, &mut g);
        ui.add_space(6.0);
        let changed = g != self.geom_edit || m != self.mesh_edit;
        self.geom_edit = g;
        self.mesh_edit = m;
        if changed {
            self.geometry_dirty = true;
            let ctx = ui.ctx().clone();
            self.edits_changed(&ctx);
        }
        let ctx = ui.ctx().clone();
        self.feedback_tick(&ctx);
        ui.add_space(6.0);
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
    fn material_regions(&self) -> Vec<(String, String)> {
        let all = &self.domain_names;
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
        ui.heading("Materials");
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
            egui::ComboBox::from_id_salt("mat_all").selected_text("apply material...").show_ui(ui, |ui| {
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
                        Some(MaterialAssign::Props(_)) => "numbers".into(),
                        None => "none".into(),
                    };
                    egui::ComboBox::from_id_salt(format!("mat_{key}")).selected_text(shown).show_ui(ui, |ui| {
                        if ui.selectable_label(matches!(cur, Some(MaterialAssign::Props(_))), "numbers").clicked() && !matches!(cur, Some(MaterialAssign::Props(_))) {
                            let map: BTreeMap<String, f64> = props.iter().map(|p| (p.name.clone(), p.default.unwrap_or_else(|| crate::examples::typical(&p.name)))).collect();
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
                            for a in euler_deg.iter_mut() {
                                edited |= num(ui, a).changed();
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
                                let mut v = map.get(&p.name).copied().unwrap_or_else(|| p.default.unwrap_or_else(|| crate::examples::typical(&p.name)));
                                ui.label(&p.label);
                                if num(ui, &mut v).changed() {
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

    /// What a library material gives, in one line, and its source file text on request.
    fn material_summary(&self, ui: &mut Ui, name: &str) {
        match self.engine.describe_material(name) {
            Ok(d) => {
                if let Some(desc) = d["description"].as_str().filter(|s| !s.is_empty()) {
                    ui.weak(desc);
                }
                let mut parts: Vec<String> = vec![];
                if let Some(c) = d["elastic"]["class"].as_str() {
                    parts.push(format!("elastic: {c}"));
                }
                if d["optical"] == true {
                    parts.push("optical model".into());
                }
                if d["thermal"] == true {
                    parts.push("thermal".into());
                }
                if d["electrical"] == true {
                    parts.push("electrical".into());
                }
                ui.weak(parts.join(", "));
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

    /// Collapsed panel to download curated presets or convert pasted database payloads (`refractiveindex.info`, `AFLOW`, `Materials Project`, `NumBAT`).
    fn material_import_panel(&mut self, ui: &mut Ui) {
        let catalog = nyx_engine::material_fetch::remote_catalog();
        let formats = ["rii", "aflow", "mp", "numbat"];
        let opt_trim = |s: &str| {
            let t = s.trim();
            (!t.is_empty()).then(|| t.to_string())
        };
        egui::CollapsingHeader::new("Download or import material").id_salt("mat_import_panel").default_open(false).show(ui, |ui| {
            ui.weak("Fetch additional entries from open databases on demand (HTTPS) or paste an offline JSON/YAML/JSONC record.");
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
                    if ui.button("Download preset").clicked() {
                        let preset = &catalog[self.mat_fetch.preset];
                        let name_ov = opt_trim(&self.mat_fetch.name);
                        let ext_ov = opt_trim(&self.mat_fetch.extends);
                        let desc_ov = opt_trim(&self.mat_fetch.description);
                        let save_dir = self.mat_fetch.save_to_library.then(nyx_engine::material_fetch::default_materials_folder).flatten();
                        match self.engine.fetch_material(preset.id, None, name_ov.as_deref(), desc_ov.as_deref(), ext_ov.as_deref(), true, save_dir.as_deref()) {
                            Ok(v) => {
                                let nm = v["name"].as_str().unwrap_or(preset.id);
                                self.mat_fetch.is_error = false;
                                self.mat_fetch.message = Some(match v["saved_path"].as_str() {
                                    Some(p) => format!("Added '{nm}' to project and saved {p}."),
                                    None => format!("Added '{nm}' to project."),
                                });
                            }
                            Err(e) => {
                                self.mat_fetch.is_error = true;
                                self.mat_fetch.message = Some(e.to_string());
                            }
                        }
                    }
                });
                ui.weak(sel.summary);
            }

            ui.add_space(4.0);
            egui::Grid::new("mat_fetch_fields").num_columns(2).show(ui, |ui| {
                ui.label("Source / format:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.source).hint_text("preset id, rii, aflow, or numbat").desired_width(220.0));
                ui.end_row();
                ui.label("Entry / query:");
                ui.add(egui::TextEdit::singleline(&mut self.mat_fetch.target).hint_text("e.g. main/Si/nk/Green-2008").desired_width(220.0));
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
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.mat_fetch.save_to_library, "save to materials/ folder");
                if ui.button("Download custom").clicked() {
                    let src = self.mat_fetch.source.trim().to_string();
                    let tgt = opt_trim(&self.mat_fetch.target);
                    let nm = opt_trim(&self.mat_fetch.name);
                    let ext = opt_trim(&self.mat_fetch.extends);
                    let desc = opt_trim(&self.mat_fetch.description);
                    let save_dir = self.mat_fetch.save_to_library.then(nyx_engine::material_fetch::default_materials_folder).flatten();
                    match self.engine.fetch_material(&src, tgt.as_deref(), nm.as_deref(), desc.as_deref(), ext.as_deref(), true, save_dir.as_deref()) {
                        Ok(v) => {
                            let added = v["name"].as_str().unwrap_or("imported");
                            self.mat_fetch.is_error = false;
                            self.mat_fetch.message = Some(format!("Added '{added}' to project."));
                        }
                        Err(e) => {
                            self.mat_fetch.is_error = true;
                            self.mat_fetch.message = Some(e.to_string());
                        }
                    }
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
        ui.heading("Physics");
        ui.horizontal_wrapped(|ui| {
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
            ui.weak("The project has no physics: add one above.");
            return;
        };
        let Some(physics) = self.engine.registry().get(&active) else { return };
        let (kinds, opts) = (physics.bc_kinds(), physics.params());
        let (exports, imports) = (physics.exports(), physics.imports());
        let Some(spec0) = self.engine.project().spec(&active).cloned() else { return };
        let mut spec = spec0.clone();
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            if ui.small_button("remove this physics from the project").clicked() {
                self.cmd("remove_physics", json!({ "physics": active }));
                self.rules_dirty = true;
            }
        });
        if !imports.is_empty() {
            let from: Vec<String> = imports.iter().map(|i| format!("{} ({})", i.role, i.field)).collect();
            ui.weak(format!("Takes from other physics of the project: {}.", from.join(", ")));
        }
        if !exports.is_empty() {
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
        if !kinds.is_empty() {
            egui::ComboBox::from_id_salt(format!("bc_add_{active}")).selected_text("+ add condition").show_ui(ui, |ui| {
                for k in &kinds {
                    if ui.selectable_label(false, &k.name).clicked() {
                        let mut params = Value::Null;
                        forms::fill_required(&k.params, &mut params);
                        spec.bcs.push(nyx_engine::facade::physics::BcSpec { on: String::new(), kind: k.kind.clone(), params });
                    }
                }
            });
        }
        if spec != spec0 {
            self.set_problem(&active, &spec);
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
            ui.weak(&s.doc);
            let runs_here = match s.dims {
                nyx_engine::facade::physics::Dims::Two => dim == 2,
                nyx_engine::facade::physics::Dims::Three => dim == 3,
                nyx_engine::facade::physics::Dims::Both => true,
            };
            if !runs_here {
                ui.colored_label(crate::theme::activity_text(ui.ctx()), format!("This study runs on a {}D mesh; the model on show is {}D.", if dim == 2 { 3 } else { 2 }, dim));
            }
            forms::param_form(ui, &format!("study_{active}_{}", s.kind), &s.params, &mut spec.study.params);
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
            num(ui, &mut self.conv_factor);
            ui.label("tolerance");
            num(ui, &mut self.conv_tol);
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
        ui.heading("Results");
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
                    let c1 = num(ui, &mut a).changed();
                    ui.label("to");
                    let c2 = num(ui, &mut b).changed();
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
                ui.horizontal(|ui| {
                    ui.label("deformation scale");
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
            let ink = if dark { Color32::from_rgb(155, 183, 217) } else { Color32::from_rgb(59, 91, 130) };
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
                    painter.circle(tip, side * if big { 0.040 } else { 0.030 }, col, egui::Stroke::new(if big { 2.0 } else { 0.8 }, Color32::WHITE));
                    painter.text(tip, egui::Align2::CENTER_CENTER, &m.label, egui::FontId::proportional(side * 0.036), Color32::WHITE);
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
            painter.circle(centre, d.radius as f32, with_alpha(Color32::from_rgb(141, 179, 221), if dark { 0.10 } else { 0.22 }), egui::Stroke::new(1.0, with_alpha(ink, 0.6)));
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
                let colour = states.iter().position(|(k, _)| *k == i).map(|p| hex(MODE_COLOURS[p % MODE_COLOURS.len()])).unwrap_or(Color32::GRAY);
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
            ui.weak(format!("{} solution(s)", nsol));
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
        egui_plot::Plot::new("summary_plot").height(ui.available_height().max(90.0)).x_axis_label(xlabel).y_axis_label(if is_eigen { "frequency (Hz)" } else if is_modes { "|neff|" } else { fname.as_str() }).legend(egui_plot::Legend::default()).show(ui, |p| {
            if is_eigen {
                p.points(egui_plot::Points::new("f", egui_plot::PlotPoints::from(xs.iter().enumerate().map(|(i, _)| [i as f64 + 1.0, xs[i]]).collect::<Vec<_>>())).radius(4.0));
            } else if is_modes {
                // neff per mode
                let pts: Vec<[f64; 2]> = rv.sols.iter().enumerate().map(|(i, s)| [i as f64 + 1.0, s.neff.map(|n| n.0.hypot(n.1)).unwrap_or(0.0)]).collect();
                p.points(egui_plot::Points::new("|neff|", egui_plot::PlotPoints::from(pts)).radius(4.0));
            } else {
                p.line(egui_plot::Line::new(format!("max {fname}"), egui_plot::PlotPoints::from(max_pts.clone())));
                p.line(egui_plot::Line::new(format!("min {fname}"), egui_plot::PlotPoints::from(min_pts.clone())));
                p.vline(egui_plot::VLine::new("now", current).color(Color32::from_rgb(230, 190, 80)));
            }
        });
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
