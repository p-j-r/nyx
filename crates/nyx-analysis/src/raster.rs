// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The software rasteriser of the figures (ADR 99): Gouraud triangles, flat triangles, segments as thin quads, text from the font of the window,
//! raster images, and Z-buffered 3D triangles/edges with back-to-front translucent shells. A figure is drawn at twice the size and averaged down,
//! which smooths the edges of everything.

use crate::canvas::{Anchor, Canvas, Item, Rgba};
use ab_glyph::{point, Font, FontRef, PxScale, ScaleFont};
use std::sync::OnceLock;

const SS: usize = 2;
/// Height (in output rows for `render`, or pixels for `raster_field_3d`) of one L2-cache-resident horizontal tile strip.
const BAND_ROWS: usize = 64;

fn font() -> Option<&'static FontRef<'static>> {
    static F: OnceLock<Option<FontRef<'static>>> = OnceLock::new();
    F.get_or_init(|| FontRef::try_from_slice(epaint_default_fonts::UBUNTU_LIGHT).ok()).as_ref()
}

/// The width of `text` at font size `size` (canvas units), as the PNG draws it. The SVG uses the same number to place centred text.
pub fn text_width(text: &str, size: f64) -> f64 {
    let Some(f) = font() else { return 0.5 * size * text.chars().count() as f64 };
    let sf = f.as_scaled(PxScale::from(size as f32));
    text.chars().map(|c| f64::from(sf.h_advance(f.glyph_id(c)))).sum()
}

/// A horizontal band `[y0, y0 + h)` of straight-alpha RGBA pixels as floats.
struct Buf {
    w: usize,
    y0: i64,
    h: usize,
    px: Vec<[f32; 4]>,
}

impl Buf {
    fn new_band(w: usize, y0: usize, h: usize, bg: [f32; 4]) -> Buf {
        Buf { w, y0: y0 as i64, h, px: vec![bg; w * h] }
    }

    #[inline]
    fn over(&mut self, x: i64, y: i64, c: [f32; 3], a: f32) {
        let ly = y - self.y0;
        if x < 0 || ly < 0 || x as usize >= self.w || ly as usize >= self.h || a <= 0.0 {
            return;
        }
        let p = &mut self.px[ly as usize * self.w + x as usize];
        let oa = a + p[3] * (1.0 - a);
        if oa <= 0.0 {
            return;
        }
        for k in 0..3 {
            p[k] = (c[k] * a + p[k] * p[3] * (1.0 - a)) / oa;
        }
        p[3] = oa;
    }

    /// A triangle with a colour per corner (the same colour three times for a flat one), blended with `alpha`.
    fn tri(&mut self, p: [[f64; 2]; 3], c: [[f32; 3]; 3], alpha: f32) {
        let area = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
        if area.abs() < 1e-12 {
            return;
        }
        let by0 = self.y0 as f64;
        let by1 = (self.y0 + self.h as i64) as f64;
        let y0_f = p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min).floor().max(by0);
        let y1_f = p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max).ceil().min(by1);
        if y0_f >= y1_f {
            return;
        }
        let x0_f = p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min).floor().max(0.0);
        let x1_f = p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max).ceil().min(self.w as f64);
        if x0_f >= x1_f {
            return;
        }
        let (x0, x1) = (x0_f as i64, x1_f as i64);
        let (y0, y1) = (y0_f as i64, y1_f as i64);
        let inv = 1.0 / area;
        for y in y0..y1 {
            let qy = y as f64 + 0.5;
            for x in x0..x1 {
                let qx = x as f64 + 0.5;
                let w1 = ((qx - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (qy - p[0][1])) * inv;
                let w2 = ((p[1][0] - p[0][0]) * (qy - p[0][1]) - (qx - p[0][0]) * (p[1][1] - p[0][1])) * inv;
                let w0 = 1.0 - w1 - w2;
                // a little outside is inside: neighbouring triangles then leave no seam
                if w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9 {
                    continue;
                }
                let (w0, w1, w2) = (w0 as f32, w1 as f32, w2 as f32);
                let col = [
                    w0 * c[0][0] + w1 * c[1][0] + w2 * c[2][0],
                    w0 * c[0][1] + w1 * c[1][1] + w2 * c[2][1],
                    w0 * c[0][2] + w1 * c[1][2] + w2 * c[2][2],
                ];
                self.over(x, y, col, alpha);
            }
        }
    }

    fn flat(&mut self, p: [[f64; 2]; 3], c: Rgba) {
        let col = rgb(c);
        self.tri(p, [col; 3], f32::from(c[3]) / 255.0);
    }

    fn segment(&mut self, a: [f64; 2], b: [f64; 2], width: f64, c: Rgba) {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = dx.hypot(dy);
        let (nx, ny) = if len > 1e-12 { (-dy / len * 0.5 * width, dx / len * 0.5 * width) } else { (0.0, 0.5 * width) };
        // the ends are squared off a half width further, which makes the corners of a polyline meet
        let (ex, ey) = if len > 1e-12 { (dx / len * 0.5 * width, dy / len * 0.5 * width) } else { (0.0, 0.0) };
        let (a, b) = ([a[0] - ex, a[1] - ey], [b[0] + ex, b[1] + ey]);
        let q = [[a[0] + nx, a[1] + ny], [b[0] + nx, b[1] + ny], [b[0] - nx, b[1] - ny], [a[0] - nx, a[1] - ny]];
        self.flat([q[0], q[1], q[2]], c);
        self.flat([q[0], q[2], q[3]], c);
    }
}

/// Pre-rasterised glyph coverage bitmap so parallel canvas strips never repeat font outlining.
struct GlyphSprite {
    x0: i64,
    y0: i64,
    w: usize,
    h: usize,
    col: [f32; 3],
    a: f32,
    cov: Vec<f32>,
}

fn build_text_sprites(x: f64, y: f64, text: &str, size: f64, c: Rgba, anchor: Anchor, out: &mut Vec<GlyphSprite>) {
    let Some(f) = font() else { return };
    let scale = PxScale::from(size as f32);
    let sf = f.as_scaled(scale);
    let width: f32 = text.chars().map(|ch| sf.h_advance(f.glyph_id(ch))).sum();
    let mut cx = match anchor {
        Anchor::Start => x as f32,
        Anchor::Middle => x as f32 - 0.5 * width,
        Anchor::End => x as f32 - width,
    };
    let col = rgb(c);
    let a = f32::from(c[3]) / 255.0;
    for ch in text.chars() {
        let id = f.glyph_id(ch);
        if let Some(og) = f.outline_glyph(id.with_scale_and_position(scale, point(cx, y as f32))) {
            let b = og.px_bounds();
            let gw = (b.max.x as i64 - b.min.x as i64).max(0) as usize;
            let gh = (b.max.y as i64 - b.min.y as i64).max(0) as usize;
            if gw > 0 && gh > 0 {
                let mut cov = vec![0.0f32; gw * gh];
                og.draw(|gx, gy, v| {
                    let (gx, gy) = (gx as usize, gy as usize);
                    if gx < gw && gy < gh {
                        cov[gy * gw + gx] = v;
                    }
                });
                out.push(GlyphSprite { x0: b.min.x as i64, y0: b.min.y as i64, w: gw, h: gh, col, a, cov });
            }
        }
        cx += sf.h_advance(id);
    }
}

fn rgb(c: Rgba) -> [f32; 3] {
    [f32::from(c[0]) / 255.0, f32::from(c[1]) / 255.0, f32::from(c[2]) / 255.0]
}

/// A triangle in pixel coordinates with an RGB colour per corner.
pub type ColouredTri = ([[f64; 2]; 3], [[f32; 3]; 3]);

/// Rasterises triangles with a colour per corner into `w` by `h` RGBA pixels; pixels no triangle covers stay transparent. Used for the field of
/// a panel, which both the SVG and the PNG carry as an image.
pub fn raster_field(w: u32, h: u32, tris: &[ColouredTri]) -> Vec<u8> {
    let (wu, hu) = (w as usize, h as usize);
    if wu == 0 || hu == 0 {
        return Vec::new();
    }
    let n_bands = hu.div_ceil(BAND_ROWS);
    let mut bins: Vec<Vec<u32>> = vec![Vec::new(); n_bands];
    for (idx, (p, _)) in tris.iter().enumerate() {
        let y0 = p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min).floor().max(0.0);
        let y1 = p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max).ceil().min(hu as f64);
        if y0 >= y1 {
            continue;
        }
        let b0 = (y0 as usize) / BAND_ROWS;
        let b1 = ((y1 as usize).saturating_sub(1) / BAND_ROWS).min(n_bands - 1);
        for b in &mut bins[b0..=b1] {
            b.push(idx as u32);
        }
    }
    let band_ids: Vec<usize> = (0..n_bands).collect();
    let chunks = femcore::par::map_chunks(&band_ids, 1, |bs| {
        let bi = bs[0];
        let y0 = bi * BAND_ROWS;
        let bh = (hu - y0).min(BAND_ROWS);
        let mut buf = Buf::new_band(wu, y0, bh, [0.0, 0.0, 0.0, 0.0]);
        for &ti in &bins[bi] {
            let (p, c) = &tris[ti as usize];
            buf.tri(*p, *c, 1.0);
        }
        let mut out = Vec::with_capacity(wu * bh * 4);
        for p in &buf.px {
            out.extend_from_slice(&[
                (p[0].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[3].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
            ]);
        }
        out
    });
    let mut out = Vec::with_capacity(wu * hu * 4);
    for c in chunks {
        out.extend_from_slice(&c);
    }
    out
}

/// A 3D triangle in pixel coordinates `(x, y)` and depth `z` (smaller `z` is nearer the camera), with an RGB colour per corner and opacity `alpha`.
#[derive(Clone, Copy, Debug)]
pub struct Tri3d {
    pub p: [[f64; 3]; 3],
    pub c: [[f32; 3]; 3],
    pub alpha: f32,
}

#[derive(Clone, Copy)]
struct PreTri3d {
    p: [[f64; 3]; 3],
    c: [[f32; 3]; 3],
    alpha: f32,
    inv_area: f64,
    x0: i64,
    x1: i64,
    y0: i64,
    y1: i64,
}

#[derive(Clone, Copy)]
struct PreEdge3d {
    a: [f64; 3],
    b: [f64; 3],
    dx: f64,
    dy: f64,
    inv_len2: f64,
    q: [[f64; 2]; 4],
    inv_q0: f64,
    inv_q1: f64,
    col: [f32; 3],
    alpha: f32,
    x0: i64,
    x1: i64,
    y0: i64,
    y1: i64,
}

fn prep_tri3d(t: &Tri3d, w: usize, h: usize) -> Option<PreTri3d> {
    if t.alpha <= 0.0 {
        return None;
    }
    let p = t.p;
    let area = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
    if area.abs() < 1e-12 || !area.is_finite() {
        return None;
    }
    let x0_f = p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min).floor().max(0.0);
    let x1_f = p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max).ceil().min(w as f64);
    let y0_f = p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min).floor().max(0.0);
    let y1_f = p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max).ceil().min(h as f64);
    if x0_f >= x1_f || y0_f >= y1_f {
        return None;
    }
    Some(PreTri3d {
        p,
        c: t.c,
        alpha: t.alpha.clamp(0.0, 1.0),
        inv_area: 1.0 / area,
        x0: x0_f as i64,
        x1: x1_f as i64,
        y0: y0_f as i64,
        y1: y1_f as i64,
    })
}

fn prep_edge3d(e: &([f64; 3], [f64; 3], Rgba, f64), w: usize, h: usize) -> Option<PreEdge3d> {
    let (a, b, rgba, width) = *e;
    if rgba[3] == 0 || width <= 0.0 {
        return None;
    }
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len2 = dx * dx + dy * dy;
    let len = len2.sqrt();
    let (nx, ny) = if len > 1e-12 { (-dy / len * 0.5 * width, dx / len * 0.5 * width) } else { (0.0, 0.5 * width) };
    let (ex, ey) = if len > 1e-12 { (dx / len * 0.5 * width, dy / len * 0.5 * width) } else { (0.0, 0.0) };
    let (ae, be) = ([a[0] - ex, a[1] - ey], [b[0] + ex, b[1] + ey]);
    let q = [[ae[0] + nx, ae[1] + ny], [be[0] + nx, be[1] + ny], [be[0] - nx, be[1] - ny], [ae[0] - nx, ae[1] - ny]];
    let area0 = (q[1][0] - q[0][0]) * (q[2][1] - q[0][1]) - (q[2][0] - q[0][0]) * (q[1][1] - q[0][1]);
    let area1 = (q[2][0] - q[0][0]) * (q[3][1] - q[0][1]) - (q[3][0] - q[0][0]) * (q[2][1] - q[0][1]);
    if area0.abs() < 1e-12 && area1.abs() < 1e-12 {
        return None;
    }
    let x0_f = q.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min).floor().max(0.0);
    let x1_f = q.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max).ceil().min(w as f64);
    let y0_f = q.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min).floor().max(0.0);
    let y1_f = q.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max).ceil().min(h as f64);
    if x0_f >= x1_f || y0_f >= y1_f {
        return None;
    }
    Some(PreEdge3d {
        a,
        b,
        dx,
        dy,
        inv_len2: if len2 > 1e-24 { 1.0 / len2 } else { 0.0 },
        q,
        inv_q0: if area0.abs() >= 1e-12 { 1.0 / area0 } else { 0.0 },
        inv_q1: if area1.abs() >= 1e-12 { 1.0 / area1 } else { 0.0 },
        col: rgb(rgba),
        alpha: f32::from(rgba[3]) / 255.0,
        x0: x0_f as i64,
        x1: x1_f as i64,
        y0: y0_f as i64,
        y1: y1_f as i64,
    })
}

/// Z-buffered 3D software rasteriser:
/// 1. `opaque` triangles are drawn with a per-pixel `f32` Z-buffer (smaller `z` is nearer the camera and wins).
/// 2. `edges` are drawn as thin quads with a small depth bias against the Z-buffer so surface edges are crisp while hidden edges are occluded.
/// 3. `translucent` triangles are sorted back-to-front by average `z` and alpha-blended where `z` passes the opaque Z-buffer.
pub fn raster_field_3d(
    w: u32,
    h: u32,
    opaque: &[Tri3d],
    translucent: &[Tri3d],
    edges: &[([f64; 3], [f64; 3], Rgba, f64)],
) -> Vec<u8> {
    let (wu, hu) = (w as usize, h as usize);
    if wu == 0 || hu == 0 {
        return Vec::new();
    }

    // Scene depth span for scale-invariant edge depth bias
    let (mut zmin, mut zmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for t in opaque.iter().chain(translucent.iter()) {
        for q in &t.p {
            if q[2].is_finite() {
                zmin = zmin.min(q[2]);
                zmax = zmax.max(q[2]);
            }
        }
    }
    for (a, b, _, _) in edges {
        for z in [a[2], b[2]] {
            if z.is_finite() {
                zmin = zmin.min(z);
                zmax = zmax.max(z);
            }
        }
    }
    let zspan = if zmax > zmin && (zmax - zmin).is_finite() { zmax - zmin } else { 1.0 };
    let depth_bias = (4e-3 * zspan).max(1e-5 * zmax.abs().max(1.0)) as f32;

    let pre_opaque: Vec<PreTri3d> = opaque.iter().filter_map(|t| prep_tri3d(t, wu, hu)).collect();
    let pre_edges: Vec<PreEdge3d> = edges.iter().filter_map(|e| prep_edge3d(e, wu, hu)).collect();

    // Sort translucent triangles back-to-front (largest average z = farthest first, down to smallest average z = nearest last)
    let mut trans_idx: Vec<(f64, usize, PreTri3d)> = translucent
        .iter()
        .enumerate()
        .filter_map(|(i, t)| {
            let pt = prep_tri3d(t, wu, hu)?;
            let avg_z = (t.p[0][2] + t.p[1][2] + t.p[2][2]) * (1.0 / 3.0);
            Some((avg_z, i, pt))
        })
        .collect();
    trans_idx.sort_by(|a, b| b.0.total_cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

    let n_bands = hu.div_ceil(BAND_ROWS);
    let mut op_bins: Vec<Vec<u32>> = vec![Vec::new(); n_bands];
    for (i, t) in pre_opaque.iter().enumerate() {
        let b0 = (t.y0 as usize) / BAND_ROWS;
        let b1 = ((t.y1 as usize).saturating_sub(1) / BAND_ROWS).min(n_bands - 1);
        for b in &mut op_bins[b0..=b1] {
            b.push(i as u32);
        }
    }
    let mut ed_bins: Vec<Vec<u32>> = vec![Vec::new(); n_bands];
    for (i, e) in pre_edges.iter().enumerate() {
        let b0 = (e.y0 as usize) / BAND_ROWS;
        let b1 = ((e.y1 as usize).saturating_sub(1) / BAND_ROWS).min(n_bands - 1);
        for b in &mut ed_bins[b0..=b1] {
            b.push(i as u32);
        }
    }
    let mut tr_bins: Vec<Vec<u32>> = vec![Vec::new(); n_bands];
    for (i, (_, _, t)) in trans_idx.iter().enumerate() {
        let b0 = (t.y0 as usize) / BAND_ROWS;
        let b1 = ((t.y1 as usize).saturating_sub(1) / BAND_ROWS).min(n_bands - 1);
        for b in &mut tr_bins[b0..=b1] {
            b.push(i as u32);
        }
    }

    let band_ids: Vec<usize> = (0..n_bands).collect();
    let chunks = femcore::par::map_chunks(&band_ids, 1, |bs| {
        let bi = bs[0];
        let by0 = bi * BAND_ROWS;
        let bh = (hu - by0).min(BAND_ROWS);
        let by0_i = by0 as i64;
        let by1_i = (by0 + bh) as i64;
        let mut px = vec![[0.0f32; 4]; wu * bh];
        let mut zbuf = vec![f32::INFINITY; wu * bh];

        let blend = |p: &mut [f32; 4], col: [f32; 3], a: f32| {
            let oa = a + p[3] * (1.0 - a);
            if oa > 0.0 {
                for k in 0..3 {
                    p[k] = (col[k] * a + p[k] * p[3] * (1.0 - a)) / oa;
                }
                p[3] = oa;
            }
        };

        // 1. Opaque pass (nearer depth wins)
        for &ti in &op_bins[bi] {
            let t = &pre_opaque[ti as usize];
            let y0 = t.y0.max(by0_i);
            let y1 = t.y1.min(by1_i);
            let p = &t.p;
            let dx10 = p[1][0] - p[0][0];
            let dy10 = p[1][1] - p[0][1];
            let dx20 = p[2][0] - p[0][0];
            let dy20 = p[2][1] - p[0][1];
            for y in y0..y1 {
                let qy = y as f64 + 0.5;
                let row = (y - by0_i) as usize * wu;
                for x in t.x0..t.x1 {
                    let qx = x as f64 + 0.5;
                    let w1 = ((qx - p[0][0]) * dy20 - dx20 * (qy - p[0][1])) * t.inv_area;
                    let w2 = (dx10 * (qy - p[0][1]) - (qx - p[0][0]) * dy10) * t.inv_area;
                    let w0 = 1.0 - w1 - w2;
                    if w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9 {
                        continue;
                    }
                    let z = (w0 * p[0][2] + w1 * p[1][2] + w2 * p[2][2]) as f32;
                    let idx = row + x as usize;
                    if z < zbuf[idx] {
                        zbuf[idx] = z;
                        let (w0f, w1f, w2f) = (w0 as f32, w1 as f32, w2 as f32);
                        let col = [
                            (w0f * t.c[0][0] + w1f * t.c[1][0] + w2f * t.c[2][0]).clamp(0.0, 1.0),
                            (w0f * t.c[0][1] + w1f * t.c[1][1] + w2f * t.c[2][1]).clamp(0.0, 1.0),
                            (w0f * t.c[0][2] + w1f * t.c[1][2] + w2f * t.c[2][2]).clamp(0.0, 1.0),
                        ];
                        if t.alpha >= 1.0 - 1e-6 {
                            px[idx] = [col[0], col[1], col[2], 1.0];
                        } else {
                            px[idx] = [0.0, 0.0, 0.0, 0.0];
                            blend(&mut px[idx], col, t.alpha);
                        }
                    }
                }
            }
        }

        // 2. Edges pass (depth-biased against Z-buffer)
        for &ei in &ed_bins[bi] {
            let e = &pre_edges[ei as usize];
            let y0 = e.y0.max(by0_i);
            let y1 = e.y1.min(by1_i);
            let q = &e.q;
            for y in y0..y1 {
                let qy = y as f64 + 0.5;
                let row = (y - by0_i) as usize * wu;
                for x in e.x0..e.x1 {
                    let qx = x as f64 + 0.5;
                    let in_t0 = if e.inv_q0 != 0.0 {
                        let w1 =
                            ((qx - q[0][0]) * (q[2][1] - q[0][1]) - (q[2][0] - q[0][0]) * (qy - q[0][1])) * e.inv_q0;
                        let w2 =
                            ((q[1][0] - q[0][0]) * (qy - q[0][1]) - (qx - q[0][0]) * (q[1][1] - q[0][1])) * e.inv_q0;
                        let w0 = 1.0 - w1 - w2;
                        w0 >= -1e-9 && w1 >= -1e-9 && w2 >= -1e-9
                    } else {
                        false
                    };
                    let in_t1 = !in_t0
                        && if e.inv_q1 != 0.0 {
                            let w1 = ((qx - q[0][0]) * (q[3][1] - q[0][1]) - (q[3][0] - q[0][0]) * (qy - q[0][1]))
                                * e.inv_q1;
                            let w2 = ((q[2][0] - q[0][0]) * (qy - q[0][1]) - (qx - q[0][0]) * (q[2][1] - q[0][1]))
                                * e.inv_q1;
                            let w0 = 1.0 - w1 - w2;
                            w0 >= -1e-9 && w1 >= -1e-9 && w2 >= -1e-9
                        } else {
                            false
                        };
                    if !(in_t0 || in_t1) {
                        continue;
                    }
                    let tpar = if e.inv_len2 > 0.0 {
                        (((qx - e.a[0]) * e.dx + (qy - e.a[1]) * e.dy) * e.inv_len2).clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                    let z = (e.a[2] + tpar * (e.b[2] - e.a[2])) as f32;
                    let idx = row + x as usize;
                    if z - depth_bias <= zbuf[idx] {
                        blend(&mut px[idx], e.col, e.alpha);
                    }
                }
            }
        }

        // 3. Translucent pass (sorted back-to-front, tested against opaque Z-buffer)
        for &ti in &tr_bins[bi] {
            let t = &trans_idx[ti as usize].2;
            let y0 = t.y0.max(by0_i);
            let y1 = t.y1.min(by1_i);
            let p = &t.p;
            let dx10 = p[1][0] - p[0][0];
            let dy10 = p[1][1] - p[0][1];
            let dx20 = p[2][0] - p[0][0];
            let dy20 = p[2][1] - p[0][1];
            for y in y0..y1 {
                let qy = y as f64 + 0.5;
                let row = (y - by0_i) as usize * wu;
                for x in t.x0..t.x1 {
                    let qx = x as f64 + 0.5;
                    let w1 = ((qx - p[0][0]) * dy20 - dx20 * (qy - p[0][1])) * t.inv_area;
                    let w2 = (dx10 * (qy - p[0][1]) - (qx - p[0][0]) * dy10) * t.inv_area;
                    let w0 = 1.0 - w1 - w2;
                    if w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9 {
                        continue;
                    }
                    let z = (w0 * p[0][2] + w1 * p[1][2] + w2 * p[2][2]) as f32;
                    let idx = row + x as usize;
                    if z <= zbuf[idx] {
                        let (w0f, w1f, w2f) = (w0 as f32, w1 as f32, w2 as f32);
                        let col = [
                            (w0f * t.c[0][0] + w1f * t.c[1][0] + w2f * t.c[2][0]).clamp(0.0, 1.0),
                            (w0f * t.c[0][1] + w1f * t.c[1][1] + w2f * t.c[2][1]).clamp(0.0, 1.0),
                            (w0f * t.c[0][2] + w1f * t.c[1][2] + w2f * t.c[2][2]).clamp(0.0, 1.0),
                        ];
                        blend(&mut px[idx], col, t.alpha);
                    }
                }
            }
        }

        let mut out = Vec::with_capacity(wu * bh * 4);
        for p in &px {
            out.extend_from_slice(&[
                (p[0].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
                (p[3].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
            ]);
        }
        out
    });

    let mut out = Vec::with_capacity(wu * hu * 4);
    for c in chunks {
        out.extend_from_slice(&c);
    }
    out
}

/// Draws the canvas at `scale` pixels per unit; returns width, height and RGBA bytes.
pub fn render(canvas: &Canvas, scale: f64) -> Result<(u32, u32, Vec<u8>), String> {
    if !(scale > 0.0 && scale <= 8.0) {
        return Err(format!("the scale of a PNG must be between 0 and 8 pixels per unit, not {scale}"));
    }
    let (ow, oh) = ((canvas.width * scale).ceil() as usize, (canvas.height * scale).ceil() as usize);
    if ow == 0 || oh == 0 || ow * oh > 60_000_000 {
        return Err(format!("a figure of {ow} by {oh} pixels is not drawn"));
    }
    let s = scale * SS as f64;
    let sw = ow * SS;
    let sh = oh * SS;
    let bg = [rgb(canvas.background)[0], rgb(canvas.background)[1], rgb(canvas.background)[2], 1.0];

    // Pre-build glyph coverage sprites per Text item so horizontal strips can run in parallel
    let mut text_sprites: Vec<Option<Vec<GlyphSprite>>> = Vec::with_capacity(canvas.items.len());
    for it in &canvas.items {
        if let Item::Text { x, y, text, size, color, anchor, bold } = it {
            let mut sp = Vec::new();
            build_text_sprites(x * s, y * s, text, size * s, *color, *anchor, &mut sp);
            if *bold {
                // the font has one weight: a second strike a fraction of a point to the right
                build_text_sprites(x * s + 0.045 * size * s, y * s, text, size * s, *color, *anchor, &mut sp);
            }
            text_sprites.push(Some(sp));
        } else {
            text_sprites.push(None);
        }
    }

    let out_band = BAND_ROWS / SS;
    let n_bands = oh.div_ceil(out_band);
    let band_ids: Vec<usize> = (0..n_bands).collect();
    let chunks = femcore::par::map_chunks(&band_ids, 1, |bs| {
        let bi = bs[0];
        let oy0 = bi * out_band;
        let obh = (oh - oy0).min(out_band);
        let sy0 = oy0 * SS;
        let sbh = obh * SS;
        let mut b = Buf::new_band(sw, sy0, sbh, bg);
        let m = |p: [f64; 2]| [p[0] * s, p[1] * s];
        let sy0_i = sy0 as i64;
        let sy1_i = (sy0 + sbh) as i64;

        for (idx, it) in canvas.items.iter().enumerate() {
            match it {
                Item::Rect { x, y, w, h, fill, stroke } => {
                    let (a, c) = (m([*x, *y]), m([x + w, y + h]));
                    if let Some(f) = fill {
                        if f[3] == 255
                            && a[0] == a[0].floor()
                            && c[0] == c[0].floor()
                            && a[1] == a[1].floor()
                            && c[1] == c[1].floor()
                        {
                            let rx0 = (a[0] as i64).max(0) as usize;
                            let rx1 = (c[0] as i64).min(sw as i64).max(rx0 as i64) as usize;
                            let ry0 = (a[1] as i64).max(sy0_i);
                            let ry1 = (c[1] as i64).min(sy1_i);
                            if rx0 < rx1 && ry0 < ry1 {
                                let col = rgb(*f);
                                let val = [col[0], col[1], col[2], 1.0];
                                for yy in ry0..ry1 {
                                    let row = (yy - sy0_i) as usize * sw;
                                    b.px[row + rx0..row + rx1].fill(val);
                                }
                            }
                        } else {
                            b.flat([a, [c[0], a[1]], c], *f);
                            b.flat([a, c, [a[0], c[1]]], *f);
                        }
                    }
                    if let Some((col, wd)) = stroke {
                        let q = [a, [c[0], a[1]], c, [a[0], c[1]], a];
                        for k in 0..4 {
                            b.segment(q[k], q[k + 1], wd * s, *col);
                        }
                    }
                }
                Item::Line { pts, color, width } => {
                    let hw = 0.5 * width * s + 1.0;
                    for k in 0..pts.len().saturating_sub(1) {
                        let (pa, pb) = (m(pts[k]), m(pts[k + 1]));
                        if pa[1].max(pb[1]) + hw < sy0 as f64 || pa[1].min(pb[1]) - hw > (sy0 + sbh) as f64 {
                            continue;
                        }
                        b.segment(pa, pb, width * s, *color);
                    }
                }
                Item::Tri { pts, fill } => b.flat([m(pts[0]), m(pts[1]), m(pts[2])], *fill),
                Item::Text { .. } => {
                    if let Some(sprites) = &text_sprites[idx] {
                        for g in sprites {
                            let gy0 = g.y0.max(sy0_i);
                            let gy1 = (g.y0 + g.h as i64).min(sy1_i);
                            if gy0 >= gy1 {
                                continue;
                            }
                            for py in gy0..gy1 {
                                let row = (py - g.y0) as usize * g.w;
                                for gx in 0..g.w {
                                    let cov = g.cov[row + gx];
                                    if cov > 0.0 {
                                        b.over(g.x0 + gx as i64, py, g.col, cov * g.a);
                                    }
                                }
                            }
                        }
                    }
                }
                Item::Image { x, y, w, h, px, py, rgba } => {
                    let (a, c) = (m([*x, *y]), m([x + w, y + h]));
                    if c[1] <= a[1] || c[0] <= a[0] {
                        continue;
                    }
                    let (x0, x1) = (a[0].floor().max(0.0) as i64, (c[0].ceil() as i64).min(sw as i64));
                    let (y0, y1) =
                        (a[1].floor().max(sy0 as f64) as i64, (c[1].ceil() as i64).min(sy1_i).min(sh as i64));
                    if x0 >= x1 || y0 >= y1 {
                        continue;
                    }
                    let inv_w = f64::from(*px) / (c[0] - a[0]);
                    let inv_h = f64::from(*py) / (c[1] - a[1]);
                    for yy in y0..y1 {
                        let v = ((yy as f64 + 0.5 - a[1]) * inv_h).floor() as i64;
                        if v < 0 || v >= i64::from(*py) {
                            continue;
                        }
                        let v_row = v as usize * *px as usize;
                        for xx in x0..x1 {
                            let u = ((xx as f64 + 0.5 - a[0]) * inv_w).floor() as i64;
                            if u < 0 || u >= i64::from(*px) {
                                continue;
                            }
                            let i = (v_row + u as usize) * 4;
                            let q = &rgba[i..i + 4];
                            if q[3] == 0 {
                                continue;
                            }
                            if q[3] == 255 {
                                b.px[(yy - sy0_i) as usize * sw + xx as usize] =
                                    [f32::from(q[0]) / 255.0, f32::from(q[1]) / 255.0, f32::from(q[2]) / 255.0, 1.0];
                            } else {
                                b.over(
                                    xx,
                                    yy,
                                    [f32::from(q[0]) / 255.0, f32::from(q[1]) / 255.0, f32::from(q[2]) / 255.0],
                                    f32::from(q[3]) / 255.0,
                                );
                            }
                        }
                    }
                }
            }
        }

        // average SS x SS pixels for this horizontal band
        let mut out = Vec::with_capacity(ow * obh * 4);
        let n = (SS * SS) as f32;
        for ly in 0..obh {
            for x in 0..ow {
                let mut acc = [0.0f32; 3];
                for dy in 0..SS {
                    let p_row = (ly * SS + dy) * sw + x * SS;
                    for dx in 0..SS {
                        let p = b.px[p_row + dx];
                        acc[0] += p[0];
                        acc[1] += p[1];
                        acc[2] += p[2];
                    }
                }
                out.extend_from_slice(&[
                    (acc[0] / n * 255.0 + 0.5) as u8,
                    (acc[1] / n * 255.0 + 0.5) as u8,
                    (acc[2] / n * 255.0 + 0.5) as u8,
                    255,
                ]);
            }
        }
        out
    });

    let mut out = Vec::with_capacity(ow * oh * 4);
    for c in chunks {
        out.extend_from_slice(&c);
    }
    Ok((ow as u32, oh as u32, out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_triangle_is_filled_with_its_colours_and_a_png_has_the_size_asked_for() {
        let mut c = Canvas::new(40.0, 30.0);
        c.push(Item::Rect {
            x: 5.0,
            y: 5.0,
            w: 30.0,
            h: 20.0,
            fill: Some([255, 0, 0, 255]),
            stroke: Some(([0, 0, 255, 255], 2.0)),
        });
        c.text(20.0, 20.0, "mode", 8.0, [0, 0, 0, 255], Anchor::Middle);
        let (w, h, px) = render(&c, 2.0).unwrap();
        assert_eq!((w, h), (80, 60));
        // the middle of the rectangle is red away from the text, the corner of the canvas is the background
        let at = |x: usize, y: usize| &px[(y * 80 + x) * 4..(y * 80 + x) * 4 + 3];
        assert_eq!(at(20, 20), &[255, 0, 0]);
        assert_eq!(at(1, 1), &[255, 255, 255]);
        // the PNG bytes start with the PNG signature
        assert_eq!(&c.to_png(2.0).unwrap()[..4], &[0x89, b'P', b'N', b'G']);
        assert!(text_width("mode", 8.0) > 8.0);
    }

    #[test]
    fn a_field_raster_interpolates_the_corner_colours_and_leaves_the_rest_transparent() {
        let t = ([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let img = raster_field(10, 10, &[t]);
        let px = |x: usize, y: usize| &img[(y * 10 + x) * 4..(y * 10 + x) * 4 + 4];
        assert_eq!(px(9, 9)[3], 0);
        let near_origin = px(0, 0);
        assert!(near_origin[0] > 200 && near_origin[3] == 255);
    }

    #[test]
    fn z_buffer_keeps_nearer_opaque_occludes_hidden_edges_and_sorts_translucent() {
        let far_blue =
            Tri3d { p: [[0.0, 0.0, 2.0], [16.0, 0.0, 2.0], [0.0, 16.0, 2.0]], c: [[0.0, 0.0, 1.0]; 3], alpha: 1.0 };
        let near_red =
            Tri3d { p: [[0.0, 0.0, 1.0], [16.0, 0.0, 1.0], [0.0, 16.0, 1.0]], c: [[1.0, 0.0, 0.0]; 3], alpha: 1.0 };
        let front_edge = ([1.0, 4.0, 1.0], [7.0, 4.0, 1.0], [0, 255, 0, 255], 2.0);
        let hidden_edge = ([1.0, 2.0, 2.5], [7.0, 2.0, 2.5], [255, 255, 0, 255], 2.0);
        // Translucent green in front of red (z = 0.5) and translucent cyan behind red (z = 1.5, should be occluded)
        let trans_behind =
            Tri3d { p: [[0.0, 0.0, 1.5], [16.0, 0.0, 1.5], [0.0, 16.0, 1.5]], c: [[0.0, 1.0, 1.0]; 3], alpha: 0.5 };
        let trans_front =
            Tri3d { p: [[0.0, 8.0, 0.5], [8.0, 8.0, 0.5], [0.0, 16.0, 0.5]], c: [[0.0, 1.0, 0.0]; 3], alpha: 0.5 };
        let img =
            raster_field_3d(16, 16, &[far_blue, near_red], &[trans_front, trans_behind], &[front_edge, hidden_edge]);
        let px = |x: usize, y: usize| &img[(y * 16 + x) * 4..(y * 16 + x) * 4 + 4];
        // Hidden yellow edge at y=2 is occluded by near_red (z=1.0) -> pure red
        assert_eq!(px(3, 2), &[255, 0, 0, 255]);
        // Coplanar front green edge at y=4 passes depth bias -> pure green
        assert_eq!(px(3, 4), &[0, 255, 0, 255]);
        // Translucent green at (1, 10) blends 50% green over red -> ~(128, 128, 0, 255)
        let blended = px(1, 10);
        assert!(blended[0] > 120 && blended[0] < 135 && blended[1] > 120 && blended[1] < 135 && blended[2] == 0);
    }
}
