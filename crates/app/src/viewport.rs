// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The 2D/3D viewport: camera, drawing of meshes and result fields, picking.

use crate::colormap::Colormap;
use crate::regions::{self, Paint, Style};
use crate::theme::{with_alpha, Ink};
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Role, Sense, Shape, Stroke, Ui, Vec2, WidgetInfo};
use crate::scene::{Colors, Scene, SceneData};
use eframe::egui::epaint::{Mesh, Vertex, WHITE_UV};
use nyx_engine::facade::Display;
use nyx_engine::facade::ViewMesh;
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Camera {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub pan: Vec2,
    pub center: [f32; 3],
    pub base_scale: f32,
    pub is3d: bool,
}

impl Default for Camera {
    fn default() -> Self {
        Camera { yaw: -0.7, pitch: 0.5, zoom: 1.0, pan: Vec2::ZERO, center: [0.0; 3], base_scale: 100.0, is3d: false }
    }
}

impl Camera {
    pub fn fit(&mut self, bbox: ([f32; 3], [f32; 3]), size: Vec2, is3d: bool, margin_right: f32) {
        let size = Vec2::new((size.x - margin_right).max(50.0), size.y);
        let (lo, hi) = bbox;
        self.is3d = is3d;
        self.center = [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0, (lo[2] + hi[2]) / 2.0];
        let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        let radius = if is3d {
            0.5 * (ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]).sqrt()
        } else {
            0.5 * ext[0].max(ext[1])
        };
        let radius = if radius > 0.0 { radius } else { 1.0 };
        self.base_scale = 0.42 * size.x.min(size.y).max(50.0) / radius;
        if !is3d {
            // fit the rectangle, not the circle
            let sx = 0.86 * size.x / ext[0].max(1e-30);
            let sy = 0.86 * size.y / ext[1].max(1e-30);
            self.base_scale = sx.min(sy);
            self.yaw = 0.0;
            self.pitch = std::f32::consts::FRAC_PI_2;
        } else {
            self.yaw = -0.7;
            self.pitch = 0.5;
        }
        self.zoom = 1.0;
        self.pan = Vec2::new(-margin_right / 2.0, 0.0);
    }

    #[inline]
    pub fn project(&self, p: [f32; 3], c: Pos2) -> (Pos2, f32) {
        let q = [p[0] - self.center[0], p[1] - self.center[1], p[2] - self.center[2]];
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let xr = cy * q[0] - sy * q[1];
        let yr = sy * q[0] + cy * q[1];
        let up = q[2] * cp + yr * sp;
        let depth = yr * cp - q[2] * sp;
        let s = self.base_scale * self.zoom;
        (Pos2::new(c.x + self.pan.x + xr * s, c.y + self.pan.y - up * s), depth)
    }

    /// World (x, y) under a screen position in 2D mode.
    pub fn unproject_2d(&self, s: Pos2, c: Pos2) -> [f32; 2] {
        let k = self.base_scale * self.zoom;
        [self.center[0] + (s.x - c.x - self.pan.x) / k, self.center[1] - (s.y - c.y - self.pan.y) / k]
    }

    pub(crate) fn view_normal_z(&self, a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> f32 {
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-30);
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        // |cos| between the face normal and the view direction: the gradient of
        // depth = yr cp - z sp with yr = sy x + cy y
        let dv = [sy * cp, cy * cp, -sp];
        ((n[0] * dv[0] + n[1] * dv[1] + n[2] * dv[2]) / len).abs()
    }
}

/// What to draw.
pub enum Content<'a> {
    Mesh {
        view: &'a ViewMesh,
        /// colour of boundary numbers that belong to a boundary condition
        bc_colors: &'a HashMap<i32, Color32>,
        selected: &'a [i32],
        show_edges: bool,
        show_numbers: bool,
        /// the colour plan of the regions (2D) or surfaces (3D), see `regions.rs`
        paint: &'a Paint,
    },
    Field {
        display: &'a Display,
        view: &'a ViewMesh,
        values: &'a [f32],
        range: (f32, f32),
        cmap: Colormap,
        unit: &'a str,
        title: &'a str,
        /// displacement of each display point (3 per point) and its scale
        warp: Option<(&'a [f32], f32)>,
        show_edges: bool,
    },
    Empty,
}

pub struct ViewOutput {
    /// boundary number clicked (in picking mode)
    pub picked: Option<i32>,
    pub probe: Option<String>,
    /// where the plot was drawn (overlays that belong to it are placed from here)
    pub rect: eframe::egui::Rect,
}

pub struct Viewport {
    pub cam: Camera,
    /// bumps when the content changes enough to refit the camera
    pub fit_token: u64,
    last_fit: u64,
    /// what was computed for the last frame (see `scene.rs`)
    scene: Scene,
    /// places of the boundary numbers (world position per label), valid for `anchors_key`
    anchors: Vec<(i32, [f32; 3])>,
    anchors_key: u64,
    /// display point of each mesh vertex (result view, for the element edges), valid for `vmap_key`
    vmap: Vec<u32>,
    vmap_key: u64,
    /// the legend row (slot) that is isolated: the other regions are dimmed and this one keeps its outline; cleared when the mesh changes
    pub isolate: Option<usize>,
    isolate_view: u64,
    /// where the number chip of every region sits (2D), valid for `region_anchors_key`
    region_anchors: Vec<(i32, [f32; 3])>,
    region_anchors_key: u64,
    /// the edges of the isolated slot that no other triangle of the slot shares, valid for `outline_key`
    outline: Vec<[u32; 2]>,
    outline_key: u64,
}

impl Default for Viewport {
    fn default() -> Self {
        Viewport {
            cam: Camera::default(),
            fit_token: 1,
            last_fit: 0,
            scene: Scene::default(),
            anchors: vec![],
            anchors_key: 0,
            vmap: vec![],
            vmap_key: 0,
            isolate: None,
            isolate_view: 0,
            region_anchors: vec![],
            region_anchors_key: 0,
            outline: vec![],
            outline_key: 0,
        }
    }
}

/// Shades a colour by the cosine `k` between the face normal and the view direction.
pub(crate) fn shade(c: Color32, k: f32) -> Color32 {
    let f = 0.35 + 0.65 * k;
    Color32::from_rgba_unmultiplied((c.r() as f32 * f) as u8, (c.g() as f32 * f) as u8, (c.b() as f32 * f) as u8, c.a())
}

pub(crate) fn point_in_tri(p: Pos2, a: Pos2, b: Pos2, c: Pos2) -> Option<[f32; 3]> {
    let d = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
    if d.abs() < 1e-12 {
        return None;
    }
    let l0 = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / d;
    let l1 = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / d;
    let l2 = 1.0 - l0 - l1;
    let e = -1e-4;
    if l0 >= e && l1 >= e && l2 >= e {
        Some([l0, l1, l2])
    } else {
        None
    }
}

fn dist_to_segment(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let t = if ab.length_sq() > 0.0 { ((p - a).dot(ab) / ab.length_sq()).clamp(0.0, 1.0) } else { 0.0 };
    (a + ab * t - p).length()
}

fn text_halo(painter: &egui::Painter, pos: Pos2, text: &str, color: Color32, halo: Color32, size: f32) {
    let font = FontId::proportional(size);
    for (dx, dy) in [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
        painter.text(pos + Vec2::new(dx, dy), Align2::CENTER_CENTER, text, font.clone(), halo);
    }
    painter.text(pos, Align2::CENTER_CENTER, text, font, color);
}

fn mix(h: u64, v: u64) -> u64 {
    (h ^ v).wrapping_mul(0x0000_0100_0000_01b3).rotate_left(23)
}

/// A few values spread over a slice, so that a new array at a recycled address does not look like the old one.
fn sample_hash(v: &[f32]) -> u64 {
    let mut h = mix(0xcbf2_9ce4_8422_2325, v.len() as u64);
    let step = (v.len() / 32).max(1);
    for i in (0..v.len()).step_by(step) {
        h = mix(h, v[i].to_bits() as u64);
    }
    if let Some(last) = v.last() {
        h = mix(h, last.to_bits() as u64);
    }
    h
}

/// Edges are drawn only when this many triangles (those on the screen) or fewer are visible: beyond it they are
/// finer than a pixel and egui would spend tens of milliseconds on lines nobody can see.
fn edge_limit(is3d: bool) -> usize {
    // egui turns every line into geometry on every frame (about 100 ns each): 25,000 triangles are 75,000 lines, ~8 ms
    if is3d {
        20_000
    } else {
        25_000
    }
}

impl Viewport {
    pub fn show(&mut self, ui: &mut Ui, content: Content<'_>, picking: bool, ink: &Ink) -> ViewOutput {
        let size = ui.available_size().max(Vec2::splat(50.0));
        let (resp, painter) = ui.allocate_painter(size, Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, ink.bg);
        let mut out = ViewOutput { picked: None, probe: None, rect };
        let (bbox, is3d, has_content) = match &content {
            Content::Mesh { view, .. } | Content::Field { view, .. } => (view.bbox, view.dim == 3, true),
            Content::Empty => (([0.0; 3], [1.0; 3]), false, false),
        };
        if !has_content {
            painter.text(rect.center(), Align2::CENTER_CENTER, "Build a geometry or open a file (File > Open, or drop it here)", FontId::proportional(15.0), ink.muted);
            return out;
        }
        if self.last_fit != self.fit_token || self.cam.is3d != is3d {
            // room on the right for the colour bar (fields) or the legend of the regions (2D meshes)
            let margin = match &content {
                Content::Field { .. } => 110.0,
                Content::Mesh { paint, .. } if !is3d && !paint.rows.is_empty() => 240.0,
                _ => 0.0,
            };
            self.cam.fit(bbox, rect.size(), is3d, margin);
            self.last_fit = self.fit_token;
        }
        let center = rect.center();

        // ---- interaction
        if resp.dragged_by(egui::PointerButton::Primary) && is3d && !picking {
            let d = resp.drag_delta();
            self.cam.yaw -= d.x * 0.01;
            self.cam.pitch = (self.cam.pitch + d.y * 0.01).clamp(-1.55, 1.55);
        } else if resp.dragged_by(egui::PointerButton::Primary) || resp.dragged_by(egui::PointerButton::Middle) || resp.dragged_by(egui::PointerButton::Secondary) {
            self.cam.pan += resp.drag_delta();
        }
        if resp.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                let f = (scroll * 0.0025).exp();
                if let Some(m) = resp.hover_pos() {
                    let rel = m - center;
                    self.cam.pan = rel - (rel - self.cam.pan) * f;
                }
                self.cam.zoom = (self.cam.zoom * f).clamp(0.02, 1e5);
            }
        }
        if resp.double_clicked() {
            self.fit_token += 1;
        }
        let cam = self.cam.clone();
        let hover = resp.hover_pos();
        let log = std::env::var_os("NYX_FRAME_LOG").is_some();
        let t_frame = std::time::Instant::now();

        match content {
            Content::Empty => {}
            Content::Mesh { view, bc_colors, selected, show_edges, show_numbers, paint } => {
                let view_key = mix(mix(0x1, view as *const ViewMesh as usize as u64), (view.tris.len() * 31 + view.points.len()) as u64);
                if self.isolate_view != view_key {
                    // another mesh: the isolated row of the legend belongs to the old one
                    self.isolate = None;
                    self.isolate_view = view_key;
                }
                if self.isolate.is_some_and(|s| s >= paint.rows.len()) {
                    self.isolate = None;
                }
                let iso = self.isolate;
                let mut key = view_key;
                key = mix(key, is3d as u64);
                // the colours of the canvas, the colour plan and the isolated row are part of what the scene colours were made from
                key = mix(mix(mix(key, ink.mode as u64 + 11), paint.hash), iso.map_or(0, |s| s as u64 + 1));
                let mut cs = 0u64;
                for (l, c) in bc_colors {
                    cs = cs.wrapping_add(mix(*l as u64, u32::from_le_bytes(c.to_array()) as u64));
                }
                key = mix(key, cs);
                for l in selected {
                    key = mix(key, *l as u64 + 1000);
                }
                let times = self.scene.update(
                    key,
                    || {
                        let colors: Vec<Color32> = view
                            .tri_label
                            .par_iter()
                            .map(|&lab| {
                                // opaque fills: a blended fill fails the colour-blind and normal-vision floors (plan section 4)
                                if is3d {
                                    match bc_colors.get(&lab) {
                                        Some(c) => *c,
                                        None => regions::fill(paint.style(lab), ink.mode),
                                    }
                                } else {
                                    let c = regions::fill(paint.style(lab), ink.mode);
                                    // a row is isolated: the other regions fade towards the canvas
                                    if iso.is_some() && paint.slot(lab) != iso { dim(c, ink.bg) } else { c }
                                }
                            })
                            .collect();
                        SceneData { pos: view.points.clone(), colors: Colors::PerTri(colors) }
                    },
                    &view.tris,
                    &cam,
                    center,
                );
                painter.add(Shape::mesh(self.scene.mesh()));
                let screen = self.scene.screen();
                // the hatch of slots 9 to 16 and of the absorbing layers (2D): thin diagonal lines over the solid fill, the one cue that does not depend on colour
                let hatched = |lab: i32| matches!(paint.style(lab), Style::Hatched(_) | Style::Pml(_)) && iso.is_none_or(|s| paint.slot(lab) == Some(s));
                let has_hatch = !is3d && view.tri_label.iter().any(|&l| hatched(l));
                if show_edges || has_hatch {
                    let limit = edge_limit(is3d);
                    let vis = self.scene.visible_tris(&view.tris, rect.min, rect.max, limit + 1);
                    if vis.len() <= limit {
                        if has_hatch {
                            let stroke = Stroke::new(1.0, with_alpha(ink.text, 200));
                            for &t in vis.iter().filter(|&&t| hatched(view.tri_label[t])) {
                                let tri = view.tris[t];
                                hatch_triangle(&painter, [screen[tri[0] as usize], screen[tri[1] as usize], screen[tri[2] as usize]], HATCH_SPACING, stroke);
                            }
                        }
                        if show_edges {
                            let stroke = Stroke::new(0.5, with_alpha(ink.text, if is3d { 45 } else { 70 }));
                            for &t in &vis {
                                let tri = view.tris[t];
                                for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                                    painter.line_segment([screen[tri[a] as usize], screen[tri[b] as usize]], stroke);
                                }
                            }
                        }
                    }
                }
                if !is3d {
                    // the isolated row keeps an outline in the text colour, never the accent: dark slot 6 is almost the brand pink (plan section 4)
                    if let Some(s) = iso {
                        let okey = mix(mix(view_key, 0x07), s as u64 + paint.hash);
                        if self.outline_key != okey {
                            self.outline = slot_outline(view, paint, s);
                            self.outline_key = okey;
                        }
                        for e in &self.outline {
                            painter.line_segment([screen[e[0] as usize], screen[e[1] as usize]], Stroke::new(2.0, ink.line));
                        }
                    }
                    for (seg, &lab) in view.segs.iter().zip(&view.seg_label) {
                        let (a, b) = (screen[seg[0] as usize], screen[seg[1] as usize]);
                        let sel = selected.contains(&lab);
                        let width = if sel { 4.0 } else if bc_colors.contains_key(&lab) { 3.0 } else { 1.6 };
                        let col = bc_colors.get(&lab).copied().unwrap_or(ink.muted);
                        // a casing in the canvas colour keeps a condition's line apart from a region fill of the same hue
                        if sel || bc_colors.contains_key(&lab) {
                            painter.line_segment([a, b], Stroke::new(width + 2.0, ink.bg));
                        }
                        painter.line_segment([a, b], Stroke::new(width, col));
                    }
                    // region number chips: with the legend row these are the second cue that is mandatory for the dark set (plan section 3)
                    let rkey = mix(view_key, 0x5e61);
                    if self.region_anchors_key != rkey {
                        self.region_anchors = region_anchors(view);
                        self.region_anchors_key = rkey;
                    }
                    for (lab, p) in &self.region_anchors {
                        let (s, _) = cam.project(*p, center);
                        if rect.expand(8.0).contains(s) {
                            chip(&painter, s, &lab.to_string(), ink);
                        }
                    }
                }
                if show_numbers {
                    let akey = mix(key, 77);
                    if self.anchors_key != akey {
                        self.anchors = label_anchors(view);
                        self.anchors_key = akey;
                    }
                    for (lab, p) in &self.anchors {
                        let (s, _) = cam.project(*p, center);
                        let col = bc_colors.get(lab).copied().unwrap_or(ink.text);
                        text_halo(&painter, s, &lab.to_string(), col, ink.bg, 13.0);
                    }
                }
                if !is3d && !paint.rows.is_empty() {
                    self.legend(ui, &resp, rect, paint, ink);
                }
                if picking && resp.clicked() {
                    if let Some(p) = resp.interact_pointer_pos() {
                        out.picked = pick_boundary(&self.scene, view, p);
                    }
                }
                if let Some(h) = hover {
                    if picking {
                        if let Some(l) = pick_boundary(&self.scene, view, h) {
                            out.probe = Some(format!("boundary {l}"));
                        }
                    } else if !view.tris.is_empty() && view.dim == 2 {
                        let w = cam.unproject_2d(h, center);
                        out.probe = Some(format!("x = {:.4e}, z = {:.4e}", w[0], w[1]));
                    }
                }
                if log {
                    eprintln!("VIEWPORT mesh: {} tris, passes (ms) data {:?} rotation {:?} camera {:?}, total {:.2} ms", view.tris.len(), times.data.map(|x| (x * 10.0).round() / 10.0), times.rotation.map(|x| (x * 10.0).round() / 10.0), times.camera.map(|x| (x * 10.0).round() / 10.0), t_frame.elapsed().as_secs_f64() * 1e3);
                }
            }
            Content::Field { display, view, values, range, cmap, unit, title, warp, show_edges } => {
                let (lo, hi) = range;
                let mut key = mix(mix(0x2, display as *const Display as usize as u64), (display.tris.len() * 31 + display.points.len()) as u64);
                key = mix(key, sample_hash(values));
                key = mix(key, values.as_ptr() as usize as u64);
                key = mix(mix(key, lo.to_bits() as u64), hi.to_bits() as u64);
                key = mix(key, cmap as u64);
                if let Some((d, s)) = warp {
                    key = mix(mix(key, sample_hash(d)), s.to_bits() as u64);
                }
                let times = self.scene.update(
                    key,
                    || {
                        let pos: Vec<[f32; 3]> = match warp {
                            Some((d, s)) => display.points.par_iter().enumerate().map(|(i, p)| [p[0] + s * d[3 * i], p[1] + s * d[3 * i + 1], p[2] + s * d[3 * i + 2]]).collect(),
                            None => display.points.clone(),
                        };
                        let span = if hi > lo { hi - lo } else { 1.0 };
                        let colors: Vec<Color32> = values.par_iter().map(|&v| cmap.sample((v - lo) / span)).collect();
                        SceneData { pos, colors: Colors::PerVertex(colors) }
                    },
                    &display.tris,
                    &cam,
                    center,
                );
                painter.add(Shape::mesh(self.scene.mesh()));
                if show_edges {
                    // coarse element edges (the original mesh), drawn on the deformed positions of its vertices when
                    // they are part of the display (always true for vertex dofs)
                    let vkey = mix(0x3, display as *const Display as usize as u64);
                    if self.vmap_key != vkey || self.vmap.len() != view.points.len() {
                        let mut m = vec![u32::MAX; view.points.len()];
                        for (i, &d) in display.dofs.iter().enumerate() {
                            if d < view.points.len() {
                                m[d] = i as u32;
                            }
                        }
                        self.vmap = m;
                        self.vmap_key = vkey;
                    }
                    let limit = edge_limit(is3d);
                    let screen = self.scene.screen();
                    let stroke = Stroke::new(0.6, Color32::from_black_alpha(110));
                    let mut drawn = 0usize;
                    for tri in &view.tris {
                        let mut s = [Pos2::ZERO; 3];
                        let mut ok = true;
                        for c in 0..3 {
                            match self.vmap[tri[c] as usize] {
                                u32::MAX => {
                                    ok = false;
                                    break;
                                }
                                i => s[c] = screen[i as usize],
                            }
                        }
                        if !ok || s.iter().all(|p| !rect.contains(*p)) {
                            continue;
                        }
                        drawn += 1;
                        if drawn > limit {
                            break;
                        }
                        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                            painter.line_segment([s[a], s[b]], stroke);
                        }
                    }
                }
                draw_colorbar(&painter, rect, range, cmap, unit, title, ink);
                if let Some(h) = hover {
                    // value under the cursor
                    if let Some((t, l, _)) = self.scene.pick(&display.tris, h) {
                        let tri = display.tris[t];
                        let v = l[0] * values[tri[0] as usize] + l[1] * values[tri[1] as usize] + l[2] * values[tri[2] as usize];
                        out.probe = Some(format!("{title} = {v:.5e} {unit}"));
                    }
                }
                if log {
                    eprintln!("VIEWPORT field: {} tris, passes (ms) data {:?} rotation {:?} camera {:?}, total {:.2} ms", display.tris.len(), times.data.map(|x| (x * 10.0).round() / 10.0), times.rotation.map(|x| (x * 10.0).round() / 10.0), times.camera.map(|x| (x * 10.0).round() / 10.0), t_frame.elapsed().as_secs_f64() * 1e3);
                }
            }
        }
        if is3d {
            draw_triad(&painter, rect, &self.cam, ink);
        }
        if let Some(text) = &out.probe {
            if let Some(h) = hover {
                let galley = painter.layout_no_wrap(text.clone(), FontId::monospace(12.0), ink.text);
                let r = Rect::from_min_size(h + Vec2::new(14.0, 14.0), galley.size() + Vec2::new(10.0, 6.0));
                painter.rect_filled(r, 4.0, with_alpha(ink.chip, 235));
                painter.rect_stroke(r, 4.0, Stroke::new(1.0, with_alpha(ink.text, 90)), egui::StrokeKind::Inside);
                painter.galley(r.min + Vec2::new(5.0, 3.0), galley, ink.text);
            }
        }
        painter.text(rect.left_bottom() + Vec2::new(8.0, -6.0), Align2::LEFT_BOTTOM, if is3d { "drag: rotate   right/middle drag: pan   wheel: zoom   double-click: reset" } else { "drag: pan   wheel: zoom   double-click: reset" }, FontId::proportional(11.0), ink.muted);
        out
    }

    /// The legend of the regions at the top right of the canvas (2D): one row per colour slot with the regions that share it. A click on a row
    /// isolates it (the others fade, the row keeps an outline); a click on the same row, or on the empty canvas, brings everything back.
    fn legend(&mut self, ui: &mut Ui, canvas: &egui::Response, rect: Rect, paint: &Paint, ink: &Ink) {
        const ROW: f32 = 22.0;
        const SHOWN: usize = 12;
        let shown = paint.rows.len().min(SHOWN);
        let more = paint.rows.len() - shown;
        let width = 220.0f32.min(rect.width() - 20.0);
        if width < 120.0 {
            return;
        }
        let height = (shown as f32 + if more > 0 { 1.0 } else { 0.0 }) * ROW + 8.0;
        let card = Rect::from_min_size(Pos2::new(rect.right() - width - 10.0, rect.top() + 10.0), Vec2::new(width, height.min(rect.height() - 20.0)));
        let painter = ui.painter_at(rect);
        painter.rect_filled(card, 6.0, with_alpha(ink.chip, 235));
        painter.rect_stroke(card, 6.0, Stroke::new(1.0, with_alpha(ink.text, 90)), egui::StrokeKind::Inside);
        for (i, row) in paint.rows.iter().take(shown).enumerate() {
            let r = Rect::from_min_size(card.min + Vec2::new(4.0, 4.0 + i as f32 * ROW), Vec2::new(width - 8.0, ROW));
            if r.bottom() > card.bottom() {
                break;
            }
            let id = egui::Id::unique("nyx-legend-row").with(row.slot);
            let resp = ui.interact(r, id, Sense::click());
            let on = self.isolate == Some(row.slot);
            let label = if row.labels.len() == 1 { row.name.clone() } else { format!("{} ({} regions)", row.name, row.labels.len()) };
            resp.widget_info(|| WidgetInfo::selected(Role::Button, true, on, format!("Show only {label}")));
            if resp.clicked() {
                self.isolate = if on { None } else { Some(row.slot) };
            }
            if on || resp.hovered() {
                painter.rect_filled(r, 4.0, with_alpha(ink.text, if on { 40 } else { 20 }));
            }
            if on {
                painter.rect_stroke(r, 4.0, Stroke::new(2.0, ink.line), egui::StrokeKind::Inside);
            }
            // the swatch is drawn the way the canvas draws the slot, hatch included
            let sw = Rect::from_min_size(r.min + Vec2::new(6.0, (ROW - 14.0) / 2.0), Vec2::splat(14.0));
            let style = row.style;
            painter.rect_filled(sw, 2.0, regions::fill(style, ink.mode));
            if matches!(style, Style::Hatched(_) | Style::Pml(_)) {
                let p = painter.with_clip_rect(sw);
                hatch_triangle(&p, [sw.left_top(), sw.right_top(), sw.left_bottom()], 4.0, Stroke::new(1.0, ink.text));
                hatch_triangle(&p, [sw.right_top(), sw.right_bottom(), sw.left_bottom()], 4.0, Stroke::new(1.0, ink.text));
            }
            painter.rect_stroke(sw, 2.0, Stroke::new(1.0, with_alpha(ink.text, 120)), egui::StrokeKind::Inside);
            let ids = row.labels.iter().take(4).map(|l| l.to_string()).collect::<Vec<_>>().join(", ") + if row.labels.len() > 4 { "…" } else { "" };
            let name = crate::fmt::ellipsize(&row.name, 16);
            painter.text(Pos2::new(sw.right() + 8.0, r.center().y), Align2::LEFT_CENTER, name, FontId::proportional(12.5), ink.text);
            painter.text(Pos2::new(r.right() - 6.0, r.center().y), Align2::RIGHT_CENTER, ids, FontId::monospace(11.0), ink.muted);
            let _ = resp.on_hover_text(format!("{label}: region {}\nClick to show only this; click again to show all.", row.labels.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(", ")));
        }
        if more > 0 {
            let r = Rect::from_min_size(card.min + Vec2::new(4.0, 4.0 + shown as f32 * ROW), Vec2::new(width - 8.0, ROW));
            if r.bottom() <= card.bottom() {
                let rest: Vec<String> = paint.rows.iter().skip(shown).map(|r| format!("{} ({})", r.name, r.labels.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(", "))).collect();
                let resp = ui.interact(r, egui::Id::unique("nyx-legend-more"), Sense::hover());
                painter.text(Pos2::new(r.left() + 8.0, r.center().y), Align2::LEFT_CENTER, format!("{more} more"), FontId::proportional(12.5), ink.muted);
                let _ = resp.on_hover_text(rest.join("\n"));
            }
        }
        // a click on the empty canvas (not on a row) clears the isolation
        if canvas.clicked() && self.isolate.is_some() && canvas.interact_pointer_pos().is_some_and(|p| !card.contains(p)) {
            self.isolate = None;
        }
    }
}

/// Lines of the hatch are this many pixels apart (plan section 5).
const HATCH_SPACING: f32 = 6.0;

/// A colour faded towards the canvas colour: what a region that is not isolated looks like.
fn dim(c: Color32, bg: Color32) -> Color32 {
    let m = |a: u8, b: u8| (f32::from(a) * 0.22 + f32::from(b) * 0.78).round() as u8;
    Color32::from_rgb(m(c.r(), bg.r()), m(c.g(), bg.g()), m(c.b(), bg.b()))
}

/// A region number on a chip in the colour of the canvas, so that the number reads the same on every region colour.
fn chip(painter: &egui::Painter, at: Pos2, text: &str, ink: &Ink) {
    let galley = painter.layout_no_wrap(text.to_string(), FontId::proportional(12.0), ink.text);
    let r = Rect::from_center_size(at, Vec2::new(galley.size().x.max(10.0) + 10.0, 18.0));
    painter.rect_filled(r, 9.0, ink.chip);
    painter.rect_stroke(r, 9.0, Stroke::new(1.0, with_alpha(ink.text, 140)), egui::StrokeKind::Inside);
    painter.galley(r.center() - galley.size() / 2.0, galley, ink.text);
}

/// 1 px lines at 45 degrees, `spacing` apart on one grid for the whole canvas (so neighbouring triangles join their lines), clipped to the triangle.
pub(crate) fn hatch_triangle(painter: &egui::Painter, tri: [Pos2; 3], spacing: f32, stroke: Stroke) {
    for s in hatch_segments(tri, spacing) {
        painter.line_segment(s, stroke);
    }
}

/// The segments of `hatch_triangle`.
fn hatch_segments(tri: [Pos2; 3], spacing: f32) -> Vec<[Pos2; 2]> {
    let mut out = vec![];
    let f = |p: Pos2| (p.x + p.y) * std::f32::consts::FRAC_1_SQRT_2;
    let fv = [f(tri[0]), f(tri[1]), f(tri[2])];
    let (lo, hi) = (fv[0].min(fv[1]).min(fv[2]), fv[0].max(fv[1]).max(fv[2]));
    if !(lo.is_finite() && hi.is_finite()) {
        return out;
    }
    let (k0, k1) = ((lo / spacing).ceil() as i64, (hi / spacing).floor() as i64);
    if k1 < k0 || k1 - k0 > 4000 {
        return out; // a triangle that spans no line, or one far larger than the screen
    }
    for k in k0..=k1 {
        let u = k as f32 * spacing;
        let mut pts = [Pos2::ZERO; 2];
        let mut n = 0;
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            let (a, b) = (fv[i], fv[j]);
            // half-open on the value, so a line through a vertex is not counted twice
            if ((a <= u && b > u) || (b <= u && a > u)) && n < 2 {
                let t = (u - a) / (b - a);
                pts[n] = tri[i] + (tri[j] - tri[i]) * t;
                n += 1;
            }
        }
        if n == 2 {
            out.push(pts);
        }
    }
    out
}

/// The edges of the triangles of one slot that belong to exactly one of them: the outline of the regions of that slot.
fn slot_outline(view: &ViewMesh, paint: &Paint, slot: usize) -> Vec<[u32; 2]> {
    let mut count: HashMap<(u32, u32), u8> = HashMap::new();
    for (tri, &lab) in view.tris.iter().zip(&view.tri_label) {
        if paint.slot(lab) != Some(slot) {
            continue;
        }
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            let (p, q) = (tri[a].min(tri[b]), tri[a].max(tri[b]));
            let c = count.entry((p, q)).or_insert(0);
            *c = c.saturating_add(1);
        }
    }
    let mut out: Vec<[u32; 2]> = count.into_iter().filter(|(_, c)| *c == 1).map(|((p, q), _)| [p, q]).collect();
    out.sort_unstable();
    out
}

/// Where the number chip of each region (2D) sits: the triangle of the region whose middle is nearest to the mean of all of its middles.
fn region_anchors(view: &ViewMesh) -> Vec<(i32, [f32; 3])> {
    let mut acc: HashMap<i32, Vec<[f32; 3]>> = HashMap::new();
    for (tri, &lab) in view.tris.iter().zip(&view.tri_label) {
        let (a, b, c) = (view.points[tri[0] as usize], view.points[tri[1] as usize], view.points[tri[2] as usize]);
        acc.entry(lab).or_default().push([(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]);
    }
    nearest_to_mean(acc)
}

/// Where the number of every boundary is written: the surface triangle (3D) or segment (2D) whose middle is nearest to
/// the mean of all of that boundary's middles, in world coordinates.
fn label_anchors(view: &ViewMesh) -> Vec<(i32, [f32; 3])> {
    let mut acc: HashMap<i32, Vec<[f32; 3]>> = HashMap::new();
    if view.dim == 3 {
        for (tri, &lab) in view.tris.iter().zip(&view.tri_label) {
            let (a, b, c) = (view.points[tri[0] as usize], view.points[tri[1] as usize], view.points[tri[2] as usize]);
            acc.entry(lab).or_default().push([(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]);
        }
    } else {
        for (seg, &lab) in view.segs.iter().zip(&view.seg_label) {
            let (a, b) = (view.points[seg[0] as usize], view.points[seg[1] as usize]);
            acc.entry(lab).or_default().push([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0]);
        }
    }
    nearest_to_mean(acc)
}

/// For every label the one of its points that is nearest to the mean of all of them, in label order.
fn nearest_to_mean(acc: HashMap<i32, Vec<[f32; 3]>>) -> Vec<(i32, [f32; 3])> {
    let mut out: Vec<(i32, [f32; 3])> = acc
        .into_iter()
        .map(|(lab, pts)| {
            let n = pts.len() as f32;
            let mean = pts.iter().fold([0.0f32; 3], |s, p| [s[0] + p[0] / n, s[1] + p[1] / n, s[2] + p[2] / n]);
            let d2 = |p: &[f32; 3]| (p[0] - mean[0]).powi(2) + (p[1] - mean[1]).powi(2) + (p[2] - mean[2]).powi(2);
            let best = pts.iter().min_by(|a, b| d2(a).total_cmp(&d2(b))).copied().unwrap_or(mean);
            (lab, best)
        })
        .collect();
    out.sort_by_key(|o| o.0);
    out
}

/// Boundary number under the screen position `p`: the nearest segment within 12 pixels (2D) or the front-most surface
/// triangle (3D), using the screen positions of the last drawn frame.
fn pick_boundary(scene: &Scene, view: &ViewMesh, p: Pos2) -> Option<i32> {
    let screen = scene.screen();
    if view.dim == 2 {
        let mut best: Option<(f32, i32)> = None;
        for (seg, &lab) in view.segs.iter().zip(&view.seg_label) {
            let d = dist_to_segment(p, screen[seg[0] as usize], screen[seg[1] as usize]);
            if d < 12.0 && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, lab));
            }
        }
        best.map(|b| b.1)
    } else {
        scene.pick(&view.tris, p).map(|(t, _, _)| view.tri_label[t])
    }
}
fn draw_colorbar(painter: &egui::Painter, rect: Rect, range: (f32, f32), cmap: Colormap, unit: &str, title: &str, ink: &Ink) {
    let w = 16.0;
    let bar = Rect::from_min_max(Pos2::new(rect.right() - 90.0, rect.top() + 40.0), Pos2::new(rect.right() - 90.0 + w, rect.bottom() - 60.0));
    if bar.height() < 40.0 {
        return;
    }
    let n = 64;
    let mut mesh = Mesh::default();
    for i in 0..=n {
        let t = i as f32 / n as f32;
        let y = bar.bottom() - t * bar.height();
        let c = cmap.sample(t);
        mesh.vertices.push(Vertex { pos: Pos2::new(bar.left(), y), uv: WHITE_UV, color: c });
        mesh.vertices.push(Vertex { pos: Pos2::new(bar.right(), y), uv: WHITE_UV, color: c });
        if i > 0 {
            let b = (2 * i) as u32;
            mesh.indices.extend_from_slice(&[b - 2, b - 1, b, b - 1, b + 1, b]);
        }
    }
    painter.add(Shape::mesh(mesh));
    painter.rect_stroke(bar, 0.0, Stroke::new(1.0, ink.muted), egui::StrokeKind::Outside);
    let font = FontId::monospace(11.0);
    let txt = ink.text;
    for (i, t) in [0.0f32, 0.25, 0.5, 0.75, 1.0].iter().enumerate() {
        let v = range.0 + (range.1 - range.0) * t;
        let y = bar.bottom() - t * bar.height();
        painter.line_segment([Pos2::new(bar.right(), y), Pos2::new(bar.right() + 4.0, y)], Stroke::new(1.0, txt));
        let _ = i;
        painter.text(Pos2::new(bar.right() + 6.0, y), Align2::LEFT_CENTER, format_value(v), font.clone(), txt);
    }
    painter.text(Pos2::new(bar.left(), bar.top() - 18.0), Align2::LEFT_BOTTOM, format!("{title} [{unit}]"), FontId::proportional(12.0), ink.text);
}

pub fn format_value(v: f32) -> String {
    if v == 0.0 {
        "0".into()
    } else if v.abs() >= 1e4 || v.abs() < 1e-2 {
        format!("{v:.3e}")
    } else {
        format!("{v:.4}")
    }
}

fn draw_triad(painter: &egui::Painter, rect: Rect, cam: &Camera, ink: &Ink) {
    let o = rect.left_bottom() + Vec2::new(46.0, -46.0);
    let mut c2 = cam.clone();
    c2.center = [0.0; 3];
    c2.pan = Vec2::ZERO;
    c2.zoom = 1.0;
    c2.base_scale = 30.0;
    for (axis, col, name) in [([1.0, 0.0, 0.0], Color32::from_rgb(230, 80, 80), "x"), ([0.0, 1.0, 0.0], Color32::from_rgb(90, 200, 90), "y"), ([0.0, 0.0, 1.0], Color32::from_rgb(90, 140, 240), "z")] {
        // the axis colours are the convention (x red, y green, z blue); on a light canvas they are darkened to stay readable
        let col = if ink.mode == crate::tokens::Mode::Light { shade(col, 0.4) } else { col };
        let (p, _) = c2.project(axis, Pos2::ZERO);
        let tip = o + p.to_vec2();
        painter.line_segment([o, tip], Stroke::new(2.0, col));
        painter.text(tip + (p.to_vec2().normalized() * 8.0), Align2::CENTER_CENTER, name, FontId::proportional(12.0), col);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn picking_finds_boundaries_in_2d_and_3d() {
        let p = femcore::meshing::geomdrive::plate(2.0, 1.0, 0.5).unwrap();
        let view = ViewMesh::from_simplex(&p.mesh);
        let mut cam = Camera::default();
        cam.fit(view.bbox, Vec2::new(800.0, 600.0), false, 0.0);
        let c = Pos2::new(400.0, 300.0);
        let mut sc = Scene::default();
        let data = |v: &ViewMesh| SceneData { pos: v.points.clone(), colors: Colors::PerTri(vec![Color32::GRAY; v.tris.len()]) };
        sc.update(1, || data(&view), &view.tris, &cam, c);
        // a point on the left edge (x = 0, second axis 0.5) is the boundary of x = 0, the bottom edge midpoint the boundary of z = 0
        let left = cam.project([0.0, 0.5, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, left + Vec2::new(3.0, 0.0)), Some(p.xmin[0]));
        let bottom = cam.project([1.0, 0.0, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, bottom + Vec2::new(0.0, -3.0)), Some(p.zmin[0]));
        let far = cam.project([1.0, 0.5, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, far), None);

        let b = femcore::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.5, 2).unwrap();
        let v3 = ViewMesh::from_simplex(&b.mesh);
        let mut cam = Camera::default();
        cam.fit(v3.bbox, Vec2::new(800.0, 600.0), true, 0.0);
        // looking from above (pitch 90 deg) the nearest facet at the centre is the top face (z = lz)
        cam.pitch = std::f32::consts::FRAC_PI_2;
        cam.yaw = 0.0;
        let mut sc = Scene::default();
        sc.update(2, || data(&v3), &v3.tris, &cam, c);
        let top = cam.project([0.5, 0.5, 1.0], c).0;
        assert_eq!(pick_boundary(&sc, &v3, top), Some(b.zmax[0]));
    }

    #[test]
    fn label_anchors_sit_on_the_boundaries() {
        let p = femcore::meshing::geomdrive::plate(2.0, 1.0, 0.5).unwrap();
        let view = ViewMesh::from_simplex(&p.mesh);
        let a = label_anchors(&view);
        let mut labels: Vec<i32> = a.iter().map(|x| x.0).collect();
        labels.sort();
        assert_eq!(labels, p.mesh.boundary_labels());
        // the boundary of x = 0 is anchored at x = 0, the boundary of x = lx at x = lx
        let at = |l: i32| a.iter().find(|x| x.0 == l).unwrap().1;
        assert!(at(p.xmin[0])[0].abs() < 1e-6 && (at(p.xmax[0])[0] - 2.0).abs() < 1e-6);
    }
}

#[cfg(test)]
mod region_tests {
    use super::*;
    use crate::tokens::Mode;

    fn tri() -> [Pos2; 3] {
        [Pos2::new(10.0, 10.0), Pos2::new(70.0, 14.0), Pos2::new(30.0, 60.0)]
    }

    #[test]
    fn hatch_lines_sit_on_one_grid_and_stay_inside_the_triangle() {
        let t = tri();
        let segs = hatch_segments(t, 6.0);
        // across the lines the triangle spans 14.1 to 63.6 px: the grid lines 3 to 10 (6 px apart) cross it
        assert_eq!(segs.len(), 8, "{segs:?}");
        for s in &segs {
            for p in s {
                // every end is on a grid line x + y = k * 6 * sqrt(2)
                let u = (p.x + p.y) * std::f32::consts::FRAC_1_SQRT_2 / 6.0;
                assert!((u - u.round()).abs() < 1e-3, "{p:?} is off the grid ({u})");
                // and on the triangle (barycentric coordinates are not negative)
                assert!(point_in_tri(*p, t[0], t[1], t[2]).is_some(), "{p:?} is outside the triangle");
            }
        }
        // a tiny triangle between two lines has none
        let tiny = [Pos2::new(10.0, 10.0), Pos2::new(11.0, 10.0), Pos2::new(10.0, 11.0)];
        assert!(hatch_segments(tiny, 6.0).len() <= 1);
        // a degenerate or non-finite triangle never panics
        assert!(hatch_segments([Pos2::new(f32::NAN, 0.0); 3], 6.0).is_empty());
    }

    #[test]
    fn two_triangles_that_share_an_edge_join_their_hatch() {
        // a square cut along its diagonal: each line crossing the diagonal ends at the same point in both halves
        let (a, b, c, d) = (Pos2::new(0.0, 0.0), Pos2::new(60.0, 0.0), Pos2::new(60.0, 60.0), Pos2::new(0.0, 60.0));
        let s1 = hatch_segments([a, b, c], 6.0);
        let s2 = hatch_segments([a, c, d], 6.0);
        let near = |p: Pos2, q: Pos2| (p - q).length() < 0.01;
        let on_diag = |p: Pos2| (p.x - p.y).abs() < 0.01;
        for s in &s1 {
            for p in s.iter().filter(|p| on_diag(**p)) {
                assert!(s2.iter().flatten().any(|q| near(*p, *q)), "{p:?} has no partner in the other triangle");
            }
        }
    }

    #[test]
    fn the_outline_of_a_slot_is_the_edges_only_one_of_its_triangles_has() {
        // two triangles of slot 0 forming a square, one triangle of slot 1 beside it
        let view = ViewMesh {
            dim: 2,
            points: vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            tris: vec![[0, 1, 2], [0, 2, 3], [1, 4, 2]],
            tri_label: vec![1, 1, 2],
            segs: vec![],
            seg_label: vec![],
            domain_labels: vec![1, 2],
            boundary_labels: vec![],
            n_cells: 3,
            bbox: ([0.0; 3], [2.0, 1.0, 0.0]),
            types: vec![],
            view_only: None,
        };
        let paint = regions::assign(&view.domain_labels, |l| l.to_string(), true);
        let o = slot_outline(&view, &paint, 0);
        // the square has four outer edges; the diagonal 0-2 is shared, and 1-2 is an outer edge of the square although slot 1 touches it
        assert_eq!(o.len(), 4, "{o:?}");
        assert!(!o.contains(&[0, 2]));
        assert!(o.contains(&[1, 2]));
        let a = region_anchors(&view);
        assert_eq!(a.iter().map(|x| x.0).collect::<Vec<_>>(), vec![1, 2]);
    }

    #[test]
    fn a_dimmed_region_is_closer_to_the_canvas_and_stays_opaque() {
        let bg = Color32::from_rgb(0x18, 0x1A, 0x1E);
        let c = Color32::from_rgb(0xC8, 0x32, 0x20);
        let d = dim(c, bg);
        assert_eq!(d.a(), 255);
        assert!((i32::from(d.r()) - i32::from(bg.r())).abs() < (i32::from(c.r()) - i32::from(bg.r())).abs());
    }

    #[test]
    fn shading_keeps_the_hue_of_a_region() {
        // the 3D faces are shaded by scaling the three channels together: the ratios, and so the hue, stay
        let hue = |c: Color32| {
            let (r, g, b) = (f32::from(c.r()), f32::from(c.g()), f32::from(c.b()));
            let (mx, mn) = (r.max(g).max(b), r.min(g).min(b));
            let d = mx - mn;
            let h = if mx == r { ((g - b) / d).rem_euclid(6.0) } else if mx == g { (b - r) / d + 2.0 } else { (r - g) / d + 4.0 };
            h * 60.0
        };
        for mode in [Mode::Dark, Mode::Light] {
            for i in 0..8 {
                let c = regions::fill(Style::Solid(i), mode);
                for k in [0.2, 0.6, 1.0] {
                    let s = shade(c, k);
                    let dh = (hue(c) - hue(s)).abs();
                    assert!(dh.min(360.0 - dh) < 3.0, "{mode:?} slot {i} at {k}: hue {} -> {}", hue(c), hue(s));
                }
            }
        }
    }

    #[test]
    fn every_region_fill_is_opaque() {
        for mode in [Mode::Dark, Mode::Light] {
            for n in 0..20 {
                assert_eq!(regions::fill(regions::style_of(n, true), mode).a(), 255);
            }
        }
    }
}
