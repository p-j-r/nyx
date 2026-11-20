//! The software rasteriser of the figures (ADR 99): Gouraud triangles, flat triangles, segments as thin quads, text from the font of the window,
//! raster images. A figure is drawn at twice the size and averaged down, which smooths the edges of everything.

use crate::canvas::{Anchor, Canvas, Item, Rgba};
use ab_glyph::{point, Font, FontRef, PxScale, ScaleFont};
use std::sync::OnceLock;

const SS: usize = 2;

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

/// A buffer of straight-alpha RGBA pixels as floats.
struct Buf {
    w: usize,
    h: usize,
    px: Vec<[f32; 4]>,
}

impl Buf {
    fn new(w: usize, h: usize, bg: [f32; 4]) -> Buf {
        Buf { w, h, px: vec![bg; w * h] }
    }

    #[inline]
    fn over(&mut self, x: i64, y: i64, c: [f32; 3], a: f32) {
        if x < 0 || y < 0 || x as usize >= self.w || y as usize >= self.h || a <= 0.0 {
            return;
        }
        let p = &mut self.px[y as usize * self.w + x as usize];
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
        let x0 = p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min).floor().max(0.0) as i64;
        let x1 = p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max).ceil().min(self.w as f64) as i64;
        let y0 = p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min).floor().max(0.0) as i64;
        let y1 = p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max).ceil().min(self.h as f64) as i64;
        let inv = 1.0 / area;
        for y in y0..y1 {
            for x in x0..x1 {
                let (qx, qy) = (x as f64 + 0.5, y as f64 + 0.5);
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

    fn text(&mut self, x: f64, y: f64, text: &str, size: f64, c: Rgba, anchor: Anchor) {
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
                og.draw(|gx, gy, cov| {
                    self.over(b.min.x as i64 + i64::from(gx), b.min.y as i64 + i64::from(gy), col, cov * a)
                });
            }
            cx += sf.h_advance(id);
        }
    }
}

fn rgb(c: Rgba) -> [f32; 3] {
    [f32::from(c[0]) / 255.0, f32::from(c[1]) / 255.0, f32::from(c[2]) / 255.0]
}

/// Rasterises triangles with a colour per corner into `w` by `h` RGBA pixels; pixels no triangle covers stay transparent. Used for the field of
/// a panel, which both the SVG and the PNG carry as an image.
/// A triangle in pixel coordinates with an RGB colour per corner.
pub type ColouredTri = ([[f64; 2]; 3], [[f32; 3]; 3]);

pub fn raster_field(w: u32, h: u32, tris: &[ColouredTri]) -> Vec<u8> {
    let mut b = Buf::new(w as usize, h as usize, [0.0, 0.0, 0.0, 0.0]);
    for (p, c) in tris {
        b.tri(*p, *c, 1.0);
    }
    b.px.iter()
        .flat_map(|p| {
            [
                (p[0] * 255.0 + 0.5) as u8,
                (p[1] * 255.0 + 0.5) as u8,
                (p[2] * 255.0 + 0.5) as u8,
                (p[3] * 255.0 + 0.5) as u8,
            ]
        })
        .collect()
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
    let mut b = Buf::new(
        ow * SS,
        oh * SS,
        [rgb(canvas.background)[0], rgb(canvas.background)[1], rgb(canvas.background)[2], 1.0],
    );
    let m = |p: [f64; 2]| [p[0] * s, p[1] * s];
    for it in &canvas.items {
        match it {
            Item::Rect { x, y, w, h, fill, stroke } => {
                let (a, c) = (m([*x, *y]), m([x + w, y + h]));
                if let Some(f) = fill {
                    b.flat([a, [c[0], a[1]], c], *f);
                    b.flat([a, c, [a[0], c[1]]], *f);
                }
                if let Some((col, wd)) = stroke {
                    let q = [a, [c[0], a[1]], c, [a[0], c[1]], a];
                    for k in 0..4 {
                        b.segment(q[k], q[k + 1], wd * s, *col);
                    }
                }
            }
            Item::Line { pts, color, width } => {
                for k in 0..pts.len().saturating_sub(1) {
                    b.segment(m(pts[k]), m(pts[k + 1]), width * s, *color);
                }
            }
            Item::Tri { pts, fill } => b.flat([m(pts[0]), m(pts[1]), m(pts[2])], *fill),
            Item::Text { x, y, text, size, color, anchor, bold } => {
                b.text(x * s, y * s, text, size * s, *color, *anchor);
                if *bold {
                    // the font has one weight: a second strike a fraction of a point to the right
                    b.text(x * s + 0.045 * size * s, y * s, text, size * s, *color, *anchor);
                }
            }
            Item::Image { x, y, w, h, px, py, rgba } => {
                let (a, c) = (m([*x, *y]), m([x + w, y + h]));
                let (x0, x1) = (a[0].floor().max(0.0) as i64, (c[0].ceil() as i64).min(b.w as i64));
                let (y0, y1) = (a[1].floor().max(0.0) as i64, (c[1].ceil() as i64).min(b.h as i64));
                for yy in y0..y1 {
                    let v = ((yy as f64 + 0.5 - a[1]) / (c[1] - a[1]) * f64::from(*py)).floor() as i64;
                    for xx in x0..x1 {
                        let u = ((xx as f64 + 0.5 - a[0]) / (c[0] - a[0]) * f64::from(*px)).floor() as i64;
                        if u < 0 || v < 0 || u >= i64::from(*px) || v >= i64::from(*py) {
                            continue;
                        }
                        let i = (v as usize * *px as usize + u as usize) * 4;
                        let q = &rgba[i..i + 4];
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
    // average SS x SS pixels
    let mut out = Vec::with_capacity(ow * oh * 4);
    for y in 0..oh {
        for x in 0..ow {
            let mut acc = [0.0f32; 3];
            for dy in 0..SS {
                for dx in 0..SS {
                    let p = b.px[(y * SS + dy) * b.w + x * SS + dx];
                    for k in 0..3 {
                        acc[k] += p[k];
                    }
                }
            }
            let n = (SS * SS) as f32;
            out.extend_from_slice(&[
                (acc[0] / n * 255.0 + 0.5) as u8,
                (acc[1] / n * 255.0 + 0.5) as u8,
                (acc[2] / n * 255.0 + 0.5) as u8,
                255,
            ]);
        }
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
}
