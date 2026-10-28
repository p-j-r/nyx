//! The 2D/3D viewport: camera, drawing of meshes and result fields, picking.

use crate::colormap::{palette, Colormap};
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Ui, Vec2};
use crate::scene::{Colors, Scene, SceneData};
use eframe::egui::epaint::{Mesh, Vertex, WHITE_UV};
use femcore::fields::Display;
use femcore::view::ViewMesh;
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

#[derive(Default)]
pub struct ViewOutput {
    /// boundary number clicked (in picking mode)
    pub picked: Option<i32>,
    pub probe: Option<String>,
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
}

impl Default for Viewport {
    fn default() -> Self {
        Viewport { cam: Camera::default(), fit_token: 1, last_fit: 0, scene: Scene::default(), anchors: vec![], anchors_key: 0, vmap: vec![], vmap_key: 0 }
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

fn text_halo(painter: &egui::Painter, pos: Pos2, text: &str, color: Color32, size: f32) {
    let font = FontId::proportional(size);
    for (dx, dy) in [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
        painter.text(pos + Vec2::new(dx, dy), Align2::CENTER_CENTER, text, font.clone(), Color32::from_black_alpha(220));
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
    pub fn show(&mut self, ui: &mut Ui, content: Content<'_>, picking: bool, bg: Color32) -> ViewOutput {
        let size = ui.available_size().max(Vec2::splat(50.0));
        let (resp, painter) = ui.allocate_painter(size, Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, bg);
        let mut out = ViewOutput::default();
        let (bbox, is3d, has_content) = match &content {
            Content::Mesh { view, .. } | Content::Field { view, .. } => (view.bbox, view.dim == 3, true),
            Content::Empty => (([0.0; 3], [1.0; 3]), false, false),
        };
        if !has_content {
            painter.text(rect.center(), Align2::CENTER_CENTER, "Build a geometry or open a file (File > Open, or drop it here)", FontId::proportional(15.0), Color32::from_gray(140));
            return out;
        }
        if self.last_fit != self.fit_token || self.cam.is3d != is3d {
            self.cam.fit(bbox, rect.size(), is3d, if matches!(content, Content::Field { .. }) { 110.0 } else { 0.0 });
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
            Content::Mesh { view, bc_colors, selected, show_edges, show_numbers } => {
                let mut key = mix(mix(0x1, view as *const ViewMesh as usize as u64), (view.tris.len() * 31 + view.points.len()) as u64);
                key = mix(key, is3d as u64);
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
                                if is3d {
                                    match bc_colors.get(&lab) {
                                        Some(c) => *c,
                                        None => palette((lab.max(1) - 1) as usize).gamma_multiply(0.55),
                                    }
                                } else {
                                    palette((lab.max(1) - 1) as usize).gamma_multiply(0.45)
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
                if show_edges {
                    let limit = edge_limit(is3d);
                    let vis = self.scene.visible_tris(&view.tris, rect.min, rect.max, limit + 1);
                    if vis.len() <= limit {
                        let stroke = Stroke::new(0.5, Color32::from_white_alpha(if is3d { 40 } else { 60 }));
                        for &t in &vis {
                            let tri = view.tris[t];
                            for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                                painter.line_segment([screen[tri[a] as usize], screen[tri[b] as usize]], stroke);
                            }
                        }
                    }
                }
                if !is3d {
                    for (seg, &lab) in view.segs.iter().zip(&view.seg_label) {
                        let (a, b) = (screen[seg[0] as usize], screen[seg[1] as usize]);
                        let sel = selected.contains(&lab);
                        let col = bc_colors.get(&lab).copied().unwrap_or(Color32::from_gray(200));
                        painter.line_segment([a, b], Stroke::new(if sel { 4.0 } else if bc_colors.contains_key(&lab) { 3.0 } else { 1.6 }, col));
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
                        let col = bc_colors.get(lab).copied().unwrap_or(Color32::WHITE);
                        text_halo(&painter, s, &lab.to_string(), col, 13.0);
                    }
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
                        out.probe = Some(format!("x = {:.4e}, y = {:.4e}", w[0], w[1]));
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
                draw_colorbar(&painter, rect, range, cmap, unit, title);
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
            draw_triad(&painter, rect, &self.cam);
        }
        if let Some(text) = &out.probe {
            if let Some(h) = hover {
                let galley = painter.layout_no_wrap(text.clone(), FontId::monospace(12.0), Color32::WHITE);
                let r = Rect::from_min_size(h + Vec2::new(14.0, 14.0), galley.size() + Vec2::new(10.0, 6.0));
                painter.rect_filled(r, 4.0, Color32::from_black_alpha(200));
                painter.galley(r.min + Vec2::new(5.0, 3.0), galley, Color32::WHITE);
            }
        }
        painter.text(rect.left_bottom() + Vec2::new(8.0, -6.0), Align2::LEFT_BOTTOM, if is3d { "drag: rotate   right/middle drag: pan   wheel: zoom   double-click: reset" } else { "drag: pan   wheel: zoom   double-click: reset" }, FontId::proportional(11.0), Color32::from_gray(110));
        out
    }
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
fn draw_colorbar(painter: &egui::Painter, rect: Rect, range: (f32, f32), cmap: Colormap, unit: &str, title: &str) {
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
    painter.rect_stroke(bar, 0.0, Stroke::new(1.0, Color32::from_gray(150)), egui::StrokeKind::Outside);
    let font = FontId::monospace(11.0);
    let txt = Color32::from_gray(210);
    for (i, t) in [0.0f32, 0.25, 0.5, 0.75, 1.0].iter().enumerate() {
        let v = range.0 + (range.1 - range.0) * t;
        let y = bar.bottom() - t * bar.height();
        painter.line_segment([Pos2::new(bar.right(), y), Pos2::new(bar.right() + 4.0, y)], Stroke::new(1.0, txt));
        let _ = i;
        painter.text(Pos2::new(bar.right() + 6.0, y), Align2::LEFT_CENTER, format_value(v), font.clone(), txt);
    }
    painter.text(Pos2::new(bar.left(), bar.top() - 18.0), Align2::LEFT_BOTTOM, format!("{title} [{unit}]"), FontId::proportional(12.0), Color32::from_gray(230));
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

fn draw_triad(painter: &egui::Painter, rect: Rect, cam: &Camera) {
    let o = rect.left_bottom() + Vec2::new(46.0, -46.0);
    let mut c2 = cam.clone();
    c2.center = [0.0; 3];
    c2.pan = Vec2::ZERO;
    c2.zoom = 1.0;
    c2.base_scale = 30.0;
    for (axis, col, name) in [([1.0, 0.0, 0.0], Color32::from_rgb(230, 80, 80), "x"), ([0.0, 1.0, 0.0], Color32::from_rgb(90, 200, 90), "y"), ([0.0, 0.0, 1.0], Color32::from_rgb(90, 140, 240), "z")] {
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
        let m = femcore::mesh::SimplexMesh::rectangle(2.0, 1.0, 4, 2);
        let view = ViewMesh::from_simplex(&m);
        let mut cam = Camera::default();
        cam.fit(view.bbox, Vec2::new(800.0, 600.0), false, 0.0);
        let c = Pos2::new(400.0, 300.0);
        let mut sc = Scene::default();
        let data = |v: &ViewMesh| SceneData { pos: v.points.clone(), colors: Colors::PerTri(vec![Color32::GRAY; v.tris.len()]) };
        sc.update(1, || data(&view), &view.tris, &cam, c);
        // a point on the left edge (x = 0, y = 0.5) is boundary 1, the bottom edge midpoint boundary 2
        let left = cam.project([0.0, 0.5, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, left + Vec2::new(3.0, 0.0)), Some(1));
        let bottom = cam.project([1.0, 0.0, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, bottom + Vec2::new(0.0, -3.0)), Some(2));
        let far = cam.project([1.0, 0.5, 0.0], c).0;
        assert_eq!(pick_boundary(&sc, &view, far), None);

        let b = femcore::mesh::SimplexMesh::block(1.0, 1.0, 1.0, 2, 2, 2);
        let v3 = ViewMesh::from_simplex(&b);
        let mut cam = Camera::default();
        cam.fit(v3.bbox, Vec2::new(800.0, 600.0), true, 0.0);
        // looking from above (pitch 90 deg) the nearest facet at the centre is the top face (4)
        cam.pitch = std::f32::consts::FRAC_PI_2;
        cam.yaw = 0.0;
        let mut sc = Scene::default();
        sc.update(2, || data(&v3), &v3.tris, &cam, c);
        let top = cam.project([0.5, 0.5, 1.0], c).0;
        assert_eq!(pick_boundary(&sc, &v3, top), Some(4));
    }

    #[test]
    fn label_anchors_sit_on_the_boundaries() {
        let m = femcore::mesh::SimplexMesh::rectangle(2.0, 1.0, 4, 2);
        let view = ViewMesh::from_simplex(&m);
        let a = label_anchors(&view);
        assert_eq!(a.iter().map(|x| x.0).collect::<Vec<_>>(), vec![1, 2, 3, 4]);
        // boundary 1 is x = 0, boundary 4 is x = lx
        assert!(a[0].1[0].abs() < 1e-6 && (a[3].1[0] - 2.0).abs() < 1e-6);
    }
}
