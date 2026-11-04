//! The Geometry tab's model editor: named parameters, the ordered feature list with one generated form per feature kind, the names the
//! model gives its regions and boundaries, and the same model as builder code. It edits `geom::model::FeatureList` and nothing else, so
//! what the forms can make the Rust builder can make and the reverse. The list is evaluated when it changes (about a millisecond for a
//! typical model); meshing stays behind the Build button.

use eframe::egui::{self, Color32, ComboBox, RichText, TextEdit, Ui};
use femcore::mesh::view::ViewMesh;
use geom::eval::Evaluated;
use geom::expr::{Env, Expr};
use geom::labels::LabelMap;
use std::sync::Arc;
use geom::model::*;
use geom::select::Selector;

const RED: Color32 = Color32::from_rgb(220, 100, 100);

/// What the editor remembers between frames: the last list it evaluated and the outcome.
#[derive(Default)]
pub struct State {
    checked: Option<FeatureList>,
    error: Option<String>,
    notice: Option<String>,
    warnings: Vec<String>,
    values: Vec<(String, Result<f64, String>)>,
    domains: Vec<(String, i32, RegionKind)>,
    boundaries: Vec<(String, i32)>,
    code: String,
    json: Option<String>,
    ev: Option<Evaluated>,
    /// the exact regions and edges as a display mesh (2D models), for the viewport before any meshing
    pub preview: Option<Arc<ViewMesh>>,
}

impl State {
    /// Evaluates the list if it differs from the one last checked.
    fn refresh(&mut self, list: &FeatureList) {
        if self.checked.as_ref() == Some(list) {
            return;
        }
        self.checked = Some(list.clone());
        self.code = list.to_rust();
        let mut env = Env::default();
        self.values = list
            .params
            .iter()
            .map(|p| {
                let v = env.eval(&p.value, &p.name).map_err(|e| e.to_string());
                if let Ok(x) = v {
                    env.vars.push((p.name.clone(), x));
                }
                (p.name.clone(), v)
            })
            .collect();
        self.domains.clear();
        self.boundaries.clear();
        self.warnings.clear();
        self.error = None;
        if let Err(e) = list.validate() {
            self.error = Some(e.to_string());
            return;
        }
        self.ev = None;
        self.preview = None;
        match femcore::meshing::geomdrive::evaluate(list).and_then(|ev| ev.label_map().map(|lm| (ev, lm)).map_err(|e| e.to_string())) {
            Ok((ev, lm)) => {
                self.warnings = ev.warnings.clone();
                self.preview = preview(&ev, &lm).map(Arc::new);
                self.ev = Some(ev);
                self.domains = lm.domains.iter().map(|d| (d.name.clone(), d.id, d.kind)).collect();
                self.boundaries = lm.boundaries.iter().map(|b| (b.name.clone(), b.id)).collect();
            }
            Err(e) => self.error = Some(e),
        }
    }
}

/// The exact regions and edges of a 2D model as a display mesh, numbered like the compiled labels (region = domain, edge = boundary), so the
/// viewport draws it as it draws a mesh.
fn preview(ev: &Evaluated, lm: &LabelMap) -> Option<ViewMesh> {
    if ev.list.dim != 2 {
        return None;
    }
    let dm = ev.tessellate(ev.bounds().diagonal() / 400.0);
    let points: Vec<[f32; 3]> = dm.points.iter().map(|p| [p[0] as f32, p[2] as f32, 0.0]).collect();
    let domain = |r: u32| lm.domains.get(r as usize).filter(|d| matches!(d.kind, RegionKind::Material | RegionKind::Pml)).map(|d| d.id);
    let (mut tris, mut tri_label) = (vec![], vec![]);
    for (t, &r) in dm.tris.iter().zip(&dm.tri_face) {
        if let Some(d) = domain(r) {
            tris.push(*t);
            tri_label.push(d);
        }
    }
    let (mut segs, mut seg_label) = (vec![], vec![]);
    for (s, &e) in dm.segs.iter().zip(&dm.seg_edge) {
        if let Some(l) = lm.label_of_edge(e) {
            segs.push(*s);
            seg_label.push(l);
        }
    }
    let mut domain_labels = tri_label.clone();
    domain_labels.sort_unstable();
    domain_labels.dedup();
    let mut boundary_labels = seg_label.clone();
    boundary_labels.sort_unstable();
    boundary_labels.dedup();
    let (mut lo, mut hi) = ([f32::INFINITY; 3], [f32::NEG_INFINITY; 3]);
    for p in &points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    if points.is_empty() {
        return None;
    }
    let n_cells = tris.len();
    Some(ViewMesh { dim: 2, points, tris, tri_label, segs, seg_label, domain_labels, boundary_labels, n_cells, bbox: (lo, hi), types: vec!["geometry".into()], view_only: None })
}

/// Edits a number-or-expression field in place. While the field has the keyboard the typed text is kept as it is, so half-typed input like `1e-` survives.
fn expr_edit(ui: &mut Ui, e: &mut Expr) -> bool {
    let id = ui.next_auto_id();
    let focused = ui.memory(|m| m.has_focus(id));
    let mut text = if focused { ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| e.source()) } else { e.source() };
    let r = ui.add(TextEdit::singleline(&mut text).id(id).desired_width(130.0).font(egui::TextStyle::Monospace));
    let mut changed = false;
    if r.changed() {
        *e = match text.trim().parse::<f64>() {
            Ok(v) if v.is_finite() => Expr::Num(v),
            _ => Expr::Text(text.clone()),
        };
        changed = true;
    }
    if r.has_focus() {
        ui.data_mut(|d| d.insert_temp(id, text));
    }
    changed
}

/// A selector as text: edits that parse replace the selector, others are shown in red and leave it as it was.
fn selector_edit(ui: &mut Ui, s: &mut Selector, ev: Option<&Evaluated>) -> bool {
    let id = ui.next_auto_id();
    let focused = ui.memory(|m| m.has_focus(id));
    let mut text = if focused { ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| s.to_string()) } else { s.to_string() };
    let r = ui.add(TextEdit::singleline(&mut text).id(id).desired_width(260.0).font(egui::TextStyle::Monospace));
    let mut changed = false;
    if r.changed() {
        if let Ok(p) = Selector::parse(&text) {
            *s = p;
            changed = true;
        }
    }
    if r.has_focus() {
        if let Err(e) = Selector::parse(&text) {
            ui.colored_label(RED, e.to_string());
        }
        ui.data_mut(|d| d.insert_temp(id, text));
    }
    if let Some(ev) = ev {
        match ev.resolve(s) {
            Ok(r) if r.set.ids.is_empty() => ui.colored_label(RED, "selects nothing"),
            Ok(r) => ui.label(RichText::new(format!("{} {:?}", r.set.ids.len(), r.set.kind).to_lowercase()).weak()),
            Err(e) => ui.colored_label(RED, e.to_string()),
        };
    }
    changed
}

fn row(ui: &mut Ui, label: &str, add: impl FnOnce(&mut Ui) -> bool) -> bool {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).weak());
        add(ui)
    })
    .inner
}

fn exprs<const N: usize>(ui: &mut Ui, label: &str, v: &mut [Expr; N]) -> bool {
    row(ui, label, |ui| {
        let mut c = false;
        for e in v.iter_mut() {
            c |= expr_edit(ui, e);
        }
        c
    })
}

fn one(ui: &mut Ui, label: &str, e: &mut Expr) -> bool {
    row(ui, label, |ui| expr_edit(ui, e))
}

fn feature_combo(ui: &mut Ui, label: &str, id: &mut FeatureId, prior: &[(FeatureId, String)]) -> bool {
    row(ui, label, |ui| {
        let shown = prior.iter().find(|(i, _)| i == id).map(|(_, n)| n.as_str()).unwrap_or("(none)");
        let mut c = false;
        ComboBox::from_id_salt(label).selected_text(shown).show_ui(ui, |ui| {
            for (i, n) in prior {
                c |= ui.selectable_value(id, *i, n).changed();
            }
        });
        c
    })
}

fn axis_combo(ui: &mut Ui, label: &str, a: &mut Axis) -> bool {
    row(ui, label, |ui| {
        let mut c = false;
        for x in [Axis::X, Axis::Y, Axis::Z] {
            c |= ui.selectable_value(a, x, x.name()).changed();
        }
        c
    })
}

fn optional(ui: &mut Ui, label: &str, e: &mut Option<Expr>) -> bool {
    row(ui, label, |ui| {
        let mut on = e.is_some();
        let mut c = ui.checkbox(&mut on, "").changed();
        if c {
            *e = if on { Some(Expr::Num(0.0)) } else { None };
        }
        if let Some(x) = e {
            c |= expr_edit(ui, x);
        }
        c
    })
}

fn id_list(ui: &mut Ui, label: &str, ids: &mut Vec<FeatureId>, prior: &[(FeatureId, String)]) -> bool {
    let mut c = false;
    ui.label(RichText::new(label).weak());
    let mut rm = None;
    for (k, id) in ids.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.push_id(k, |ui| c |= feature_combo(ui, "", id, prior));
            if ui.small_button("✕").clicked() {
                rm = Some(k);
            }
        });
    }
    if let Some(k) = rm {
        ids.remove(k);
        c = true;
    }
    if ui.small_button("+").clicked() {
        if let Some((i, _)) = prior.last() {
            ids.push(*i);
            c = true;
        }
    }
    c
}

fn kind_name(k: RegionKind) -> &'static str {
    match k {
        RegionKind::Material => "material",
        RegionKind::Hole => "hole",
        RegionKind::Pml => "pml",
        RegionKind::Void => "void",
    }
}

fn step_form(ui: &mut Ui, s: &mut TStep) -> bool {
    let mut c = false;
    match s {
        TStep::Translate(v) => c |= exprs(ui, "translate", v),
        TStep::Rotate { axis, deg, about } => {
            c |= axis_combo(ui, "rotate about", axis);
            c |= one(ui, "degrees", deg);
            c |= exprs(ui, "through", about);
        }
        TStep::Mirror { axis, at } => {
            c |= axis_combo(ui, "mirror across", axis);
            c |= one(ui, "at", at);
        }
        TStep::Scale { factor, about } => {
            c |= one(ui, "scale by", factor);
            c |= exprs(ui, "about", about);
        }
    }
    c
}

fn zero3() -> [Expr; 3] {
    [Expr::Num(0.0), Expr::Num(0.0), Expr::Num(0.0)]
}

/// The form of one feature: a field per number, feature reference or selector it has. Returns whether anything changed.
fn op_form(ui: &mut Ui, op: &mut Op, prior: &[(FeatureId, String)], ev: Option<&Evaluated>) -> bool {
    let mut c = false;
    match op {
        Op::Rect { min, max } => {
            c |= exprs(ui, "min (x, z)", min);
            c |= exprs(ui, "max (x, z)", max);
        }
        Op::RoundedRect { min, max, r } => {
            c |= exprs(ui, "min (x, z)", min);
            c |= exprs(ui, "max (x, z)", max);
            c |= one(ui, "corner radius", r);
        }
        Op::Circle { center, r } => {
            c |= exprs(ui, "center (x, z)", center);
            c |= one(ui, "radius", r);
        }
        Op::Ellipse { center, radii, angle_deg } => {
            c |= exprs(ui, "center (x, z)", center);
            c |= exprs(ui, "radii", radii);
            c |= one(ui, "angle (deg)", angle_deg);
        }
        Op::Polygon { pts, names } => {
            let mut rm = None;
            for (k, p) in pts.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("point {k}")).weak());
                    c |= expr_edit(ui, &mut p[0]);
                    c |= expr_edit(ui, &mut p[1]);
                    if ui.small_button("✕").clicked() {
                        rm = Some(k);
                    }
                });
            }
            if let Some(k) = rm {
                if pts.len() > 3 {
                    pts.remove(k);
                    names.clear();
                    c = true;
                }
            }
            if ui.small_button("+ point").clicked() {
                let last = pts.last().cloned().unwrap_or([Expr::Num(0.0), Expr::Num(0.0)]);
                pts.push(last);
                names.clear();
                c = true;
            }
        }
        Op::Box { min, max } => {
            c |= exprs(ui, "min (x, y, z)", min);
            c |= exprs(ui, "max (x, y, z)", max);
        }
        Op::Cylinder { base, axis, r, h } => {
            c |= exprs(ui, "base", base);
            c |= axis_combo(ui, "axis", axis);
            c |= one(ui, "radius", r);
            c |= one(ui, "height", h);
        }
        Op::Sphere { center, r } => {
            c |= exprs(ui, "center", center);
            c |= one(ui, "radius", r);
        }
        Op::Cone { base, axis, r0, r1, h } => {
            c |= exprs(ui, "base", base);
            c |= axis_combo(ui, "axis", axis);
            c |= one(ui, "radius at base", r0);
            c |= one(ui, "radius at top", r1);
            c |= one(ui, "height", h);
        }
        Op::Revolve { profile, axis, center } => {
            c |= feature_combo(ui, "profile", profile, prior);
            c |= axis_combo(ui, "axis", axis);
            c |= exprs(ui, "center", center);
        }
        Op::Torus { center, major, minor } => {
            c |= exprs(ui, "center", center);
            c |= one(ui, "major radius", major);
            c |= one(ui, "minor radius", minor);
        }
        Op::Extrude { profile, axis, from, to, taper_deg } => {
            c |= feature_combo(ui, "profile", profile, prior);
            c |= axis_combo(ui, "along", axis);
            c |= one(ui, "from", from);
            c |= one(ui, "to", to);
            c |= optional(ui, "taper (deg)", taper_deg);
        }
        Op::Transform { input, steps } => {
            c |= feature_combo(ui, "of", input, prior);
            let mut rm = None;
            for (k, s) in steps.iter_mut().enumerate() {
                ui.push_id(k, |ui| {
                    ui.horizontal(|ui| {
                        ui.label(RichText::new(format!("step {}", k + 1)).strong());
                        if ui.small_button("✕").clicked() {
                            rm = Some(k);
                        }
                    });
                    c |= step_form(ui, s);
                });
            }
            if let Some(k) = rm {
                steps.remove(k);
                c = true;
            }
            ui.horizontal(|ui| {
                if ui.small_button("+ translate").clicked() {
                    steps.push(TStep::Translate(zero3()));
                    c = true;
                }
                if ui.small_button("+ rotate").clicked() {
                    steps.push(TStep::Rotate { axis: Axis::Z, deg: Expr::Num(90.0), about: zero3() });
                    c = true;
                }
                if ui.small_button("+ mirror").clicked() {
                    steps.push(TStep::Mirror { axis: Axis::X, at: Expr::Num(0.0) });
                    c = true;
                }
                if ui.small_button("+ scale").clicked() {
                    steps.push(TStep::Scale { factor: Expr::Num(2.0), about: zero3() });
                    c = true;
                }
            });
        }
        Op::Pattern { input, kind, n } => {
            c |= feature_combo(ui, "of", input, prior);
            c |= one(ui, "copies", n);
            let linear = matches!(kind, PatternKind::Linear { .. });
            ui.horizontal(|ui| {
                if ui.selectable_label(linear, "linear").clicked() && !linear {
                    *kind = PatternKind::Linear { step: [Expr::Num(1.0), Expr::Num(0.0), Expr::Num(0.0)] };
                    c = true;
                }
                if ui.selectable_label(!linear, "polar").clicked() && linear {
                    *kind = PatternKind::Polar { axis: Axis::Z, center: zero3(), total_deg: Expr::Num(360.0) };
                    c = true;
                }
            });
            match kind {
                PatternKind::Linear { step } => c |= exprs(ui, "step", step),
                PatternKind::Polar { axis, center, total_deg } => {
                    c |= axis_combo(ui, "about", axis);
                    c |= exprs(ui, "center", center);
                    c |= one(ui, "total angle (deg)", total_deg);
                }
            }
        }
        Op::Union { ids } | Op::Intersect { ids } => c |= id_list(ui, "of", ids, prior),
        Op::Difference { a, b } => {
            c |= feature_combo(ui, "start from", a, prior);
            c |= id_list(ui, "minus", b, prior);
        }
        Op::Overlay { regions, pins, naming } => {
            ui.label(RichText::new("Regions in painting order: where several contain a point, the last one wins.").weak());
            let mut rm = None;
            let mut swap = None;
            let n = regions.len();
            for (k, r) in regions.iter_mut().enumerate() {
                ui.push_id(k, |ui| {
                    ui.horizontal(|ui| {
                        c |= ui.add(TextEdit::singleline(&mut r.name).desired_width(100.0)).changed();
                        let shown = prior.iter().find(|(i, _)| *i == r.input).map(|(_, n)| n.as_str()).unwrap_or("(none)");
                        ComboBox::from_id_salt("input").selected_text(shown).show_ui(ui, |ui| {
                            for (i, nm) in prior {
                                c |= ui.selectable_value(&mut r.input, *i, nm).changed();
                            }
                        });
                        ComboBox::from_id_salt("kind").selected_text(kind_name(r.kind)).show_ui(ui, |ui| {
                            for kd in [RegionKind::Material, RegionKind::Hole, RegionKind::Pml, RegionKind::Void] {
                                c |= ui.selectable_value(&mut r.kind, kd, kind_name(kd)).changed();
                            }
                        });
                        if k > 0 && ui.small_button("↑").clicked() {
                            swap = Some(k);
                        }
                        if k + 1 < n && ui.small_button("↓").clicked() {
                            swap = Some(k + 1);
                        }
                        if ui.small_button("✕").clicked() {
                            rm = Some(k);
                        }
                    });
                });
            }
            if let Some(k) = swap {
                regions.swap(k - 1, k);
                c = true;
            }
            if let Some(k) = rm {
                regions.remove(k);
                c = true;
            }
            if ui.small_button("+ region").clicked() {
                if let Some((i, nm)) = prior.last() {
                    regions.push(RegionDef { name: format!("{nm}_region"), input: *i, kind: RegionKind::Material, pinned_id: None });
                    c = true;
                }
            }
            row(ui, "boundary names after", |ui| {
                let mut ch = false;
                for (v, label) in [(NameScheme::Feature, "the shape"), (NameScheme::Region, "the region"), (NameScheme::Numbered, "region and number")] {
                    ch |= ui.selectable_value(naming, v, label).changed();
                }
                c |= ch;
                ch
            });
            if !pins.is_empty() {
                ui.weak(format!("{} boundary label pin(s) kept from an older project", pins.len()));
            }
        }
        Op::Pml { window, sides, thickness, layers, grading, regions } => {
            c |= feature_combo(ui, "window", window, prior);
            row(ui, "sides (none: all)", |ui| {
                let mut ch = false;
                for s in [Side::XMin, Side::XMax, Side::YMin, Side::YMax, Side::ZMin, Side::ZMax] {
                    let mut on = sides.contains(&s);
                    if ui.checkbox(&mut on, s.name()).changed() {
                        if on {
                            sides.push(s);
                        } else {
                            sides.retain(|x| *x != s);
                        }
                        ch = true;
                    }
                }
                c |= ch;
                ch
            });
            c |= one(ui, "thickness", thickness);
            c |= one(ui, "layers", layers);
            c |= one(ui, "grading order", grading);
            row(ui, "regions (none: all)", |ui| {
                let mut text = regions.join(", ");
                let ch = ui.add(TextEdit::singleline(&mut text).desired_width(200.0)).changed();
                if ch {
                    *regions = text.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
                }
                c |= ch;
                ch
            });
        }
        Op::Cut { normal, offset, keep } => {
            c |= exprs(ui, "plane normal", normal);
            c |= one(ui, "offset", offset);
            row(ui, "keep", |ui| {
                let a = ui.selectable_value(keep, Keep::Above, "side the normal points to").changed();
                let b = ui.selectable_value(keep, Keep::Below, "the other side").changed();
                c |= a | b;
                a | b
            });
        }
        Op::Periodic { sel_a, sel_b, translation } => {
            c |= row(ui, "face A", |ui| selector_edit(ui, sel_a, ev));
            c |= row(ui, "face B", |ui| selector_edit(ui, sel_b, ev));
            c |= exprs(ui, "translation A to B", translation);
        }
        Op::Size { sel, h, max_deviation, min_layers } => {
            c |= row(ui, "applies to", |ui| selector_edit(ui, sel, ev));
            c |= optional(ui, "element size", h);
            c |= optional(ui, "largest chord deviation", max_deviation);
            c |= optional(ui, "layers across thin parts", min_layers);
        }
        Op::Tag { sel, name } => {
            c |= row(ui, "applies to", |ui| selector_edit(ui, sel, ev));
            c |= row(ui, "tag", |ui| ui.add(TextEdit::singleline(name).desired_width(120.0)).changed());
        }
    }
    c
}

/// The kinds the Add menu offers, grouped.
const GROUPS: [(&str, &[&str]); 5] = [
    ("Profiles", &["Rect", "RoundedRect", "Circle", "Ellipse", "Polygon"]),
    ("Solids", &["Box", "Cylinder", "Cone", "Sphere", "Torus", "Extrude", "Revolve"]),
    ("Combine", &["Union", "Difference", "Intersect", "Transform", "Pattern"]),
    ("Regions", &["Overlay"]),
    ("Simulation", &["Pml", "Size", "Tag", "Cut", "Periodic"]),
];

fn e2(a: f64, b: f64) -> [Expr; 2] {
    [Expr::Num(a), Expr::Num(b)]
}

/// A new feature of the given kind with plausible values, referring to the latest earlier features where it needs inputs.
fn default_op(kind: &str, list: &FeatureList) -> Option<Op> {
    let last = list.features.last().map(|f| f.id);
    let before = list.features.iter().rev().nth(1).map(|f| f.id);
    let three = |v: f64| [Expr::Num(v), Expr::Num(v), Expr::Num(v)];
    Some(match kind {
        "Rect" => Op::Rect { min: e2(0.0, 0.0), max: e2(1.0, 1.0) },
        "RoundedRect" => Op::RoundedRect { min: e2(0.0, 0.0), max: e2(1.0, 1.0), r: Expr::Num(0.1) },
        "Circle" => Op::Circle { center: e2(0.5, 0.5), r: Expr::Num(0.25) },
        "Ellipse" => Op::Ellipse { center: e2(0.5, 0.5), radii: e2(0.4, 0.2), angle_deg: Expr::Num(0.0) },
        "Polygon" => Op::Polygon { pts: vec![e2(0.0, 0.0), e2(1.0, 0.0), e2(0.5, 1.0)], names: vec![] },
        "Box" => Op::Box { min: three(0.0), max: three(1.0) },
        "Cylinder" => Op::Cylinder { base: [Expr::Num(0.5), Expr::Num(0.5), Expr::Num(0.0)], axis: Axis::Z, r: Expr::Num(0.25), h: Expr::Num(1.0) },
        "Cone" => Op::Cone { base: [Expr::Num(0.5), Expr::Num(0.5), Expr::Num(0.0)], axis: Axis::Z, r0: Expr::Num(0.3), r1: Expr::Num(0.1), h: Expr::Num(1.0) },
        "Sphere" => Op::Sphere { center: three(0.5), r: Expr::Num(0.25) },
        "Torus" => Op::Torus { center: three(0.5), major: Expr::Num(0.3), minor: Expr::Num(0.1) },
        "Extrude" => Op::Extrude { profile: last?, axis: Axis::Y, from: Expr::Num(0.0), to: Expr::Num(1.0), taper_deg: None },
        "Revolve" => Op::Revolve { profile: last?, axis: Axis::Z, center: zero3() },
        "Union" => Op::Union { ids: vec![before?, last?] },
        "Intersect" => Op::Intersect { ids: vec![before?, last?] },
        "Difference" => Op::Difference { a: before?, b: vec![last?] },
        "Transform" => Op::Transform { input: last?, steps: vec![TStep::Translate(zero3())] },
        "Pattern" => Op::Pattern { input: last?, kind: PatternKind::Linear { step: [Expr::Num(1.0), Expr::Num(0.0), Expr::Num(0.0)] }, n: Expr::Num(3.0) },
        "Overlay" => {
            let regions = list
                .features
                .iter()
                .filter(|f| !matches!(f.op, Op::Overlay { .. } | Op::Pml { .. } | Op::Size { .. } | Op::Tag { .. } | Op::Cut { .. } | Op::Periodic { .. }))
                .map(|f| RegionDef { name: f.name.clone(), input: f.id, kind: RegionKind::Material, pinned_id: None })
                .collect();
            Op::Overlay { regions, pins: vec![], naming: NameScheme::Feature }
        }
        "Pml" => Op::Pml { window: last?, sides: vec![], thickness: Expr::Num(1.0), layers: Expr::Num(4.0), grading: Expr::Num(2.0), regions: vec![] },
        "Size" => Op::Size { sel: Selector::parse("outer").ok()?, h: Some(Expr::Num(0.1)), max_deviation: None, min_layers: None },
        "Tag" => Op::Tag { sel: Selector::parse("outer").ok()?, name: "tag".into() },
        "Cut" => Op::Cut { normal: [Expr::Num(1.0), Expr::Num(0.0), Expr::Num(0.0)], offset: Expr::Num(0.0), keep: Keep::Above },
        "Periodic" => Op::Periodic { sel_a: Selector::parse("plane(x = 0)").ok()?, sel_b: Selector::parse("plane(x = 1)").ok()?, translation: [Expr::Num(1.0), Expr::Num(0.0), Expr::Num(0.0)] },
        _ => return None,
    })
}

fn is_3d_kind(kind: &str) -> bool {
    matches!(kind, "Box" | "Cylinder" | "Cone" | "Sphere" | "Torus" | "Extrude" | "Revolve")
}

fn summary(op: &Op) -> String {
    let v = |e: &Expr| e.source();
    match op {
        Op::Rect { min, max } => format!("({}, {}) to ({}, {})", v(&min[0]), v(&min[1]), v(&max[0]), v(&max[1])),
        Op::Circle { center, r } => format!("r {}  at ({}, {})", v(r), v(&center[0]), v(&center[1])),
        Op::Polygon { pts, .. } => format!("{} points", pts.len()),
        Op::Union { ids } | Op::Intersect { ids } => format!("{} features", ids.len()),
        Op::Overlay { regions, .. } => format!("{} regions", regions.len()),
        Op::Size { h: Some(h), .. } => format!("h {}", v(h)),
        _ => String::new(),
    }
}

/// The model editor. `list` is edited in place; the return value says whether it changed.
pub fn show(ui: &mut Ui, st: &mut State, list: &mut FeatureList) -> bool {
    let before = list.clone();
    st.refresh(list);

    if let Some(e) = &st.error {
        ui.colored_label(RED, e);
    }
    if let Some(n) = &st.notice {
        ui.colored_label(RED, n);
    }
    for w in &st.warnings {
        ui.colored_label(Color32::from_rgb(220, 170, 80), w);
    }

    // ---- parameters
    ui.label(RichText::new("Parameters").strong());
    let mut rm = None;
    egui::Grid::new("geom_params").num_columns(4).show(ui, |ui| {
        for (k, p) in list.params.iter_mut().enumerate() {
            ui.add(TextEdit::singleline(&mut p.name).desired_width(80.0));
            ui.label("=");
            expr_edit(ui, &mut p.value);
            match st.values.get(k).map(|(_, v)| v) {
                Some(Ok(x)) => ui.label(RichText::new(format!("{x:.6e}")).weak().monospace()),
                Some(Err(e)) => ui.colored_label(RED, e),
                None => ui.label(""),
            };
            if ui.small_button("✕").clicked() {
                rm = Some(k);
            }
            ui.end_row();
        }
    });
    if let Some(k) = rm {
        list.params.remove(k);
    }
    if ui.small_button("+ parameter").clicked() {
        let n = list.params.len() + 1;
        list.params.push(Param { name: format!("p{n}"), value: Expr::Num(1.0) });
    }

    // ---- features
    ui.add_space(6.0);
    ui.label(RichText::new("Features").strong());
    let heads: Vec<(FeatureId, String)> = list.features.iter().map(|f| (f.id, f.name.clone())).collect();
    let used: Vec<FeatureId> = list.features.iter().flat_map(|f| f.op.inputs()).collect();
    let mut remove = None;
    let mut move_by: Option<(usize, isize)> = None;
    for (k, f) in list.features.iter_mut().enumerate() {
        let prior = &heads[..k];
        let title = format!("{}  ·  {} {}", f.name, f.op.kind(), summary(&f.op));
        egui::CollapsingHeader::new(RichText::new(title).color(if f.suppressed { Color32::GRAY } else { ui.visuals().text_color() })).id_salt(("feature", f.id)).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.add(TextEdit::singleline(&mut f.name).desired_width(120.0));
                let mut on = !f.suppressed;
                if ui.checkbox(&mut on, "active").changed() {
                    f.suppressed = !on;
                }
                if ui.small_button("↑").clicked() {
                    move_by = Some((k, -1));
                }
                if ui.small_button("↓").clicked() {
                    move_by = Some((k, 1));
                }
                if ui.small_button("remove").clicked() {
                    remove = Some(k);
                }
                if f.origin == Origin::Code {
                    ui.weak("made by code");
                }
            });
            op_form(ui, &mut f.op, prior, st.ev.as_ref());
        });
    }
    st.notice = None;
    if let Some(k) = remove {
        let id = list.features[k].id;
        if used.contains(&id) {
            st.notice = Some(format!("'{}' is used by another feature; remove or change that one first", list.features[k].name));
        } else {
            list.features.remove(k);
        }
    }
    if let Some((k, d)) = move_by {
        let j = k as isize + d;
        if j >= 0 && (j as usize) < list.features.len() {
            list.features.swap(k, j as usize);
            if list.validate().is_err() {
                list.features.swap(k, j as usize);
                st.notice = Some("a feature cannot come before the features it refers to".into());
            }
        }
    }
    let mut add: Option<&str> = None;
    ui.horizontal(|ui| {
        ComboBox::from_id_salt("add_feature").selected_text("+ add feature").show_ui(ui, |ui| {
            for (group, kinds) in GROUPS {
                ui.label(RichText::new(group).weak());
                for kind in kinds.iter().filter(|k| list.dim == 3 || !is_3d_kind(k)) {
                    if *kind == "Overlay" && list.features.iter().any(|f| matches!(f.op, Op::Overlay { .. })) {
                        continue;
                    }
                    if ui.selectable_label(false, *kind).clicked() {
                        add = Some(kind);
                    }
                }
            }
        });
    });
    if let Some(kind) = add {
        match default_op(kind, list) {
            Some(op) => {
                let id = list.features.iter().map(|f| f.id).max().map_or(1, |m| m + 1);
                let mut n = list.features.len() + 1;
                let lower = kind.to_lowercase();
                while list.features.iter().any(|f| f.name == format!("{lower}_{n}")) {
                    n += 1;
                }
                list.features.push(Feature { id, name: format!("{lower}_{n}"), op, suppressed: false, origin: Origin::Gui, tags: vec![] });
            }
            None => st.notice = Some(format!("a {kind} needs an earlier feature to refer to")),
        }
    }

    // ---- names
    ui.add_space(6.0);
    ui.collapsing("Names and numbers", |ui| {
        ui.weak("Problems refer to regions and boundaries by these names or by a rule; the integers are what the solvers store. Read-only.");
        egui::Grid::new("model_domains").num_columns(3).show(ui, |ui| {
            ui.label(RichText::new("region").strong());
            ui.label(RichText::new("kind").strong());
            ui.label(RichText::new("number").strong());
            ui.end_row();
            for (n, id, kd) in &st.domains {
                ui.label(n);
                ui.label(kind_name(*kd));
                ui.label(id.to_string());
                ui.end_row();
            }
        });
        ui.add_space(4.0);
        egui::Grid::new("model_boundaries").num_columns(2).show(ui, |ui| {
            ui.label(RichText::new("boundary").strong());
            ui.label(RichText::new("number").strong());
            ui.end_row();
            for (n, id) in &st.boundaries {
                ui.label(n);
                ui.label(id.to_string());
                ui.end_row();
            }
        });
    });

    // ---- code
    ui.collapsing("Show as code", |ui| {
        ui.weak("The same model as Rust builder code; running it gives this feature list back.");
        if ui.small_button("copy").clicked() {
            ui.ctx().copy_text(st.code.clone());
        }
        let mut code = st.code.clone();
        ui.add(TextEdit::multiline(&mut code).code_editor().desired_rows(8).desired_width(f32::INFINITY).interactive(false));
    });
    ui.collapsing("Feature list (JSON)", |ui| {
        let text = st.json.get_or_insert_with(|| serde_json::to_string_pretty(&*list).unwrap_or_default());
        ui.add(TextEdit::multiline(text).code_editor().desired_rows(10).desired_width(f32::INFINITY));
        let parsed = serde_json::from_str::<FeatureList>(text);
        ui.horizontal(|ui| {
            if ui.add_enabled(parsed.is_ok(), egui::Button::new("Apply")).clicked() {
                if let Ok(l) = &parsed {
                    *list = l.clone();
                }
                st.json = None;
            } else if ui.button("Reset").clicked() {
                st.json = None;
            } else if let Err(e) = &parsed {
                ui.colored_label(RED, e.to_string());
            }
        });
    });
    // keep the JSON text in step with the form once it is closed again
    if *list != before {
        st.json = None;
    }
    *list != before
}


#[cfg(test)]
mod tests {
    use super::*;

    fn draw(list: &mut FeatureList, st: &mut State) {
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            show(ui, st, list);
        });
        out.textures_delta.clear();
    }

    /// Every kind in the Add menu gets a form that draws and a default that the list accepts.
    #[test]
    fn every_kind_can_be_added_and_drawn() {
        for dim in [2u8, 3] {
            let mut list = FeatureList::new(dim);
            let mut st = State::default();
            let mut seed = Model::new_2d();
            let r = seed.rect("base", (0.0, 0.0), (1.0, 1.0));
            let c = seed.circle("hole", (0.5, 0.5), 0.2);
            list.features = seed.list.features;
            let _ = (r, c);
            for (_, kinds) in GROUPS {
                for kind in kinds.iter().filter(|k| dim == 3 || !is_3d_kind(k)) {
                    let op = default_op(kind, &list).unwrap_or_else(|| panic!("no default for {kind}"));
                    let id = list.features.iter().map(|f| f.id).max().unwrap() + 1;
                    list.features.push(Feature { id, name: format!("{}_{id}", kind.to_lowercase()), op, suppressed: false, origin: Origin::Gui, tags: vec![] });
                }
            }
            list.validate().unwrap_or_else(|e| panic!("{dim}D: {e}"));
            draw(&mut list, &mut st);
            assert!(st.code.contains("let mut m"), "{}", st.code);
        }
    }

    /// The model of the SOI wire built with the form defaults and edits equals the one the builder makes.
    #[test]
    fn a_model_made_from_the_defaults_equals_the_builder_model() {
        let mut m = Model::new_2d();
        let a = m.rect("clad", (-1.0, -1.0), (1.0, 1.0));
        let b = m.rect("core", (-0.25, 0.0), (0.25, 0.3));
        m.overlay(&[("clad", a), ("Si", b)]);
        let mut list = FeatureList::new(2);
        for (name, min, max) in [("clad", (-1.0, -1.0), (1.0, 1.0)), ("core", (-0.25, 0.0), (0.25, 0.3))] {
            let mut op = default_op("Rect", &list).unwrap();
            if let Op::Rect { min: lo, max: hi } = &mut op {
                *lo = e2(min.0, min.1);
                *hi = e2(max.0, max.1);
            }
            let id = list.features.len() as u32 + 1;
            list.features.push(Feature { id, name: name.into(), op, suppressed: false, origin: Origin::Gui, tags: vec![] });
        }
        let mut op = default_op("Overlay", &list).unwrap();
        if let Op::Overlay { regions, .. } = &mut op {
            regions[1].name = "Si".into();
        }
        list.features.push(Feature { id: 3, name: "overlay".into(), op, suppressed: false, origin: Origin::Gui, tags: vec![] });
        list.features.iter_mut().for_each(|f| f.origin = Origin::Gui);
        let mut want = m.list.clone();
        want.features.iter_mut().for_each(|f| f.origin = Origin::Gui);
        assert_eq!(list, want);
        let mut st = State::default();
        st.refresh(&list);
        assert!(st.error.is_none(), "{:?}", st.error);
        assert_eq!(st.domains.iter().map(|d| d.0.as_str()).collect::<Vec<_>>(), ["clad", "Si"]);
    }

    /// The preview of a model carries the same numbers the mesh will: one domain per material region, one boundary label per edge.
    #[test]
    fn the_preview_is_numbered_like_the_mesh_and_selectors_report_their_hits() {
        let p = crate::examples::soi_wire();
        let mut st = State::default();
        st.refresh(&p.geometry.list);
        let v = st.preview.clone().expect("a 2D model has a preview");
        let mesh = p.geometry.build().unwrap();
        assert_eq!(v.domain_labels, mesh.domain_labels());
        // the preview also draws the interfaces between regions, which the mesh has no boundary facets on
        assert!(mesh.boundary_labels().iter().all(|l| v.boundary_labels.contains(l)));
        assert!(v.tris.len() > 5 && v.segs.len() > 5);
        let ev = st.ev.as_ref().unwrap();
        assert!(!ev.resolve(&Selector::parse("outer").unwrap()).unwrap().set.ids.is_empty());
        assert!(ev.resolve(&Selector::parse("plane(z = 123)").unwrap()).unwrap().set.ids.is_empty());
    }

    #[test]
    fn expressions_are_numbers_when_they_parse() {
        let mut m = Model::new_2d();
        m.param("w", 0.5);
        let mut list = m.list;
        list.params.push(Param { name: "half".into(), value: Expr::Text("w/2".into()) });
        let mut st = State::default();
        st.refresh(&list);
        assert_eq!(st.values[1].1, Ok(0.25));
        list.params.push(Param { name: "bad".into(), value: Expr::Text("zz+1".into()) });
        st.refresh(&list);
        assert!(st.values[2].1.is_err());
    }
}
