// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What the viewport draws, kept between frames.
//!
//! Drawing a mesh or a result field means: colour every vertex, project the vertices to the screen, sort the
//! triangles back to front (3D) and fill an egui `Mesh`. For a million triangles that is far too much to repeat
//! on every frame, so the work is split by what it depends on and each part is redone only when its inputs change:
//!
//! * the **data** (which mesh, which field and colour map and range, which deformation, which boundary colours):
//!   world positions of the vertices and the base colour of every triangle or vertex;
//! * the **rotation** (yaw and pitch, 3D only): depth of every vertex, the back-to-front order of the triangles and
//!   their shading (an orthographic view: panning and zooming do not change the order);
//! * the **camera** (rotation, pan, zoom, window size): the screen position of every vertex and the filled `Mesh`.
//!
//! Every pass runs on the shared thread pool. The finished mesh is handed to egui as an `Arc`, so a frame in which
//! nothing changed copies nothing. Picking and hover use the cached screen positions, in parallel.

use crate::viewport::{shade, Camera};
use eframe::egui::epaint::{Mesh, Vertex, WHITE_UV};
use eframe::egui::{Color32, Pos2};
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;

/// Base colours: one per triangle (mesh view) or one per vertex (result field).
pub enum Colors {
    PerTri(Vec<Color32>),
    PerVertex(Vec<Color32>),
}

/// What the data pass produces.
pub struct SceneData {
    /// world position of every vertex
    pub pos: Vec<[f32; 3]>,
    pub colors: Colors,
}

/// Which passes ran in the last `update`, and how long they took (milliseconds), for the frame log.
#[derive(Clone, Copy, Debug, Default)]
pub struct SceneTimes {
    pub data: Option<f64>,
    pub rotation: Option<f64>,
    pub camera: Option<f64>,
}

pub struct Scene {
    data_key: u64,
    is3d: bool,
    pos: Vec<[f32; 3]>,
    colors: Colors,
    rot_key: [u32; 3],
    depth: Vec<f32>,
    order: Vec<u32>,
    shade: Vec<f32>,
    cam_key: [u32; 10],
    screen: Vec<Pos2>,
    mesh: Arc<Mesh>,
    valid: bool,
}

impl Default for Scene {
    fn default() -> Self {
        Scene { data_key: 0, is3d: false, pos: vec![], colors: Colors::PerTri(vec![]), rot_key: [0; 3], depth: vec![], order: vec![], shade: vec![], cam_key: [0; 10], screen: vec![], mesh: Arc::new(Mesh::default()), valid: false }
    }
}

fn bits(x: f32) -> u32 {
    x.to_bits()
}

impl Scene {
    /// Brings the cache up to date. `make_data` is only called when `data_key` differs from the last call; `cam` and
    /// `center` (the middle of the drawing rectangle) say how the world is seen.
    pub fn update(&mut self, data_key: u64, make_data: impl FnOnce() -> SceneData, tris: &[[u32; 3]], cam: &Camera, center: Pos2) -> SceneTimes {
        let mut times = SceneTimes::default();
        let is3d = cam.is3d;
        let mut rot_dirty = false;
        if !self.valid || self.data_key != data_key || self.is3d != is3d {
            let t = Instant::now();
            let d = make_data();
            self.pos = d.pos;
            self.colors = d.colors;
            self.data_key = data_key;
            self.is3d = is3d;
            self.valid = true;
            self.rot_key = [u32::MAX; 3]; // force the rotation pass
            self.cam_key = [u32::MAX; 10];
            times.data = Some(t.elapsed().as_secs_f64() * 1e3);
            rot_dirty = true;
        }
        let rot_key = [bits(cam.yaw), bits(cam.pitch), is3d as u32];
        if rot_dirty || self.rot_key != rot_key {
            let t = Instant::now();
            self.rotation_pass(tris, cam);
            self.rot_key = rot_key;
            self.cam_key = [u32::MAX; 10];
            times.rotation = Some(t.elapsed().as_secs_f64() * 1e3);
        }
        let cam_key = [bits(cam.yaw), bits(cam.pitch), bits(cam.zoom), bits(cam.pan.x), bits(cam.pan.y), bits(cam.base_scale), bits(cam.center[0] + 3.0 * cam.center[1] + 7.0 * cam.center[2]), bits(center.x), bits(center.y), is3d as u32];
        if self.cam_key != cam_key {
            let t = Instant::now();
            self.camera_pass(tris, cam, center);
            self.cam_key = cam_key;
            times.camera = Some(t.elapsed().as_secs_f64() * 1e3);
        }
        times
    }

    /// Depth of every vertex, triangles sorted back to front, shading of every triangle (3D); identity order (2D).
    fn rotation_pass(&mut self, tris: &[[u32; 3]], cam: &Camera) {
        let n = tris.len();
        if !self.is3d {
            self.depth.clear();
            self.shade.clear();
            if self.order.len() != n || self.order.iter().enumerate().any(|(i, &o)| o as usize != i) {
                self.order = (0..n as u32).collect();
            }
            return;
        }
        let zero = Pos2::ZERO;
        let pos = &self.pos;
        self.depth = pos.par_iter().map(|p| cam.project(*p, zero).1).collect();
        let depth = &self.depth;
        let mut keyed: Vec<(f32, u32)> = tris.par_iter().enumerate().map(|(t, tri)| (depth[tri[0] as usize] + depth[tri[1] as usize] + depth[tri[2] as usize], t as u32)).collect();
        // far triangles first; ties broken by index so the order does not depend on the scheduling
        keyed.par_sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        self.order = keyed.into_iter().map(|k| k.1).collect();
        self.shade = tris.par_iter().map(|tri| cam.view_normal_z(pos[tri[0] as usize], pos[tri[1] as usize], pos[tri[2] as usize])).collect();
    }

    /// Screen positions and the filled mesh.
    fn camera_pass(&mut self, tris: &[[u32; 3]], cam: &Camera, center: Pos2) {
        let pos = &self.pos;
        self.screen = pos.par_iter().map(|p| cam.project(*p, center).0).collect();
        let n = tris.len();
        let (screen, order, shade_k, is3d) = (&self.screen, &self.order, &self.shade, self.is3d);
        let colors = &self.colors;
        // reuse the allocation of the previous mesh when egui no longer holds it
        let unique = Arc::get_mut(&mut self.mesh).is_some();
        let mut local = Mesh::default();
        let mesh: &mut Mesh = if unique { Arc::get_mut(&mut self.mesh).unwrap() } else { &mut local };
        mesh.vertices.resize(3 * n, Vertex { pos: Pos2::ZERO, uv: WHITE_UV, color: Color32::TRANSPARENT });
        if mesh.indices.len() != 3 * n {
            mesh.indices = (0..3 * n as u32).collect();
        }
        mesh.vertices.par_chunks_mut(3 * 2048).enumerate().for_each(|(ci, chunk)| {
            for (j, tv) in chunk.chunks_mut(3).enumerate() {
                let slot = ci * 2048 + j;
                let t = order[slot] as usize;
                let tri = tris[t];
                let k = if is3d { shade_k[t] } else { 1.0 };
                for c in 0..3 {
                    let v = tri[c] as usize;
                    let base = match colors {
                        Colors::PerTri(tc) => tc[t],
                        Colors::PerVertex(vc) => vc[v],
                    };
                    tv[c] = Vertex { pos: screen[v], uv: WHITE_UV, color: if is3d { shade(base, k) } else { base } };
                }
            }
        });
        if !unique {
            // the old mesh is still shared with egui: replace it by the freshly built one
            self.mesh = Arc::new(local);
        }
    }

    /// The mesh for egui (shares the allocation: cloning the `Arc` copies nothing).
    pub fn mesh(&self) -> Arc<Mesh> {
        self.mesh.clone()
    }

    pub fn screen(&self) -> &[Pos2] {
        &self.screen
    }

    #[cfg(test)]
    pub fn depth(&self) -> &[f32] {
        &self.depth
    }

    /// The front-most triangle under `p` (index into `tris`), its barycentric coordinates and its depth; `None` when
    /// `p` is over no triangle. Searches all triangles on the pool.
    pub fn pick(&self, tris: &[[u32; 3]], p: Pos2) -> Option<(usize, [f32; 3], f32)> {
        let screen = &self.screen;
        let depth = &self.depth;
        let is3d = self.is3d;
        tris.par_iter()
            .enumerate()
            .filter_map(|(t, tri)| {
                let (a, b, c) = (screen[tri[0] as usize], screen[tri[1] as usize], screen[tri[2] as usize]);
                if p.x < a.x.min(b.x).min(c.x) || p.x > a.x.max(b.x).max(c.x) || p.y < a.y.min(b.y).min(c.y) || p.y > a.y.max(b.y).max(c.y) {
                    return None;
                }
                let l = crate::viewport::point_in_tri(p, a, b, c)?;
                let d = if is3d { l[0] * depth[tri[0] as usize] + l[1] * depth[tri[1] as usize] + l[2] * depth[tri[2] as usize] } else { 0.0 };
                Some((t, l, d))
            })
            .reduce_with(|x, y| if y.2 < x.2 || (y.2 == x.2 && y.0 < x.0) { y } else { x })
    }

    /// Triangles that touch the rectangle `[lo, hi]` on the screen, at most `limit` of them (for drawing edges).
    pub fn visible_tris(&self, tris: &[[u32; 3]], lo: Pos2, hi: Pos2, limit: usize) -> Vec<usize> {
        let screen = &self.screen;
        let mut v: Vec<usize> = tris
            .par_iter()
            .enumerate()
            .filter_map(|(t, tri)| {
                let (a, b, c) = (screen[tri[0] as usize], screen[tri[1] as usize], screen[tri[2] as usize]);
                let (xmin, xmax) = (a.x.min(b.x).min(c.x), a.x.max(b.x).max(c.x));
                let (ymin, ymax) = (a.y.min(b.y).min(c.y), a.y.max(b.y).max(c.y));
                (xmax >= lo.x && xmin <= hi.x && ymax >= lo.y && ymin <= hi.y).then_some(t)
            })
            .collect();
        v.truncate(limit);
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui::Vec2;

    fn unit_square_tris() -> (Vec<[f32; 3]>, Vec<[u32; 3]>) {
        (vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], vec![[0, 1, 2], [0, 2, 3]])
    }

    #[test]
    fn passes_run_only_when_their_inputs_change() {
        let (pos, tris) = unit_square_tris();
        let mut cam = Camera::default();
        cam.fit(([0.0; 3], [1.0, 1.0, 0.0]), Vec2::new(400.0, 300.0), false, 0.0);
        let mut sc = Scene::default();
        let make = || SceneData { pos: pos.clone(), colors: Colors::PerTri(vec![Color32::RED, Color32::BLUE]) };
        let c = Pos2::new(200.0, 150.0);
        let t = sc.update(1, make, &tris, &cam, c);
        assert!(t.data.is_some() && t.rotation.is_some() && t.camera.is_some());
        assert_eq!(sc.mesh().vertices.len(), 6);
        // nothing changed: nothing runs and make_data is not called
        let t = sc.update(1, || panic!("data must not be rebuilt"), &tris, &cam, c);
        assert!(t.data.is_none() && t.rotation.is_none() && t.camera.is_none());
        // panning only touches the camera pass
        cam.pan += Vec2::new(5.0, 0.0);
        let t = sc.update(1, || panic!("data must not be rebuilt"), &tris, &cam, c);
        assert!(t.data.is_none() && t.rotation.is_none() && t.camera.is_some());
        // new data key rebuilds everything
        let t = sc.update(2, make, &tris, &cam, c);
        assert!(t.data.is_some() && t.rotation.is_some() && t.camera.is_some());
    }

    #[test]
    fn three_d_triangles_come_out_far_first_and_picking_finds_the_front_one() {
        // two parallel squares, the second one 1 unit further from the viewer (larger depth)
        let mut pos = vec![];
        let mut tris = vec![];
        for (k, z) in [0.0f32, 1.0].iter().enumerate() {
            let b = (4 * k) as u32;
            pos.extend([[0.0, 0.0, *z], [1.0, 0.0, *z], [1.0, 1.0, *z], [0.0, 1.0, *z]]);
            tris.extend([[b, b + 1, b + 2], [b, b + 2, b + 3]]);
        }
        let mut cam = Camera::default();
        cam.fit(([0.0; 3], [1.0, 1.0, 1.0]), Vec2::new(400.0, 300.0), true, 0.0);
        // looking along -z from above: pitch 90 degrees, depth grows as z decreases
        cam.yaw = 0.0;
        cam.pitch = std::f32::consts::FRAC_PI_2;
        let mut sc = Scene::default();
        let c = Pos2::new(200.0, 150.0);
        sc.update(7, || SceneData { pos: pos.clone(), colors: Colors::PerTri(vec![Color32::GRAY; 4]) }, &tris, &cam, c);
        let d: Vec<f32> = sc.order.iter().map(|&t| tris[t as usize].iter().map(|&v| sc.depth()[v as usize]).sum()).collect();
        assert!(d.windows(2).all(|w| w[0] >= w[1]), "far first: {d:?}");
        let mid = cam.project([0.7, 0.3, 0.5], c).0;
        let (t, _, _) = sc.pick(&tris, mid).expect("a triangle under the point");
        let nearest_depth = tris[t].iter().map(|&v| sc.depth()[v as usize]).sum::<f32>();
        assert!(d.iter().all(|x| *x >= nearest_depth - 1e-5), "picked the front one");
        assert!(sc.pick(&tris, Pos2::new(5.0, 5.0)).is_none());
    }
}
