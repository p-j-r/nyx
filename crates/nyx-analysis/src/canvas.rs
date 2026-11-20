//! A figure as a list of drawing items (rectangles, polylines, triangles, text, raster images) on a canvas, written as SVG or as PNG (ADR 99).
//!
//! The items are what both outputs can draw. A field on a mesh is not a vector drawing: thousands of coloured triangles make a heavy SVG that
//! renders slowly, so a field is rasterised into an `Image` item (`raster_field`) and the SVG embeds it as a PNG; the axes, outlines, arrows and
//! text stay vector in the SVG. The PNG output draws everything itself (software rasteriser, drawn at twice the size and averaged down, text
//! with `ab_glyph` in the font of the window): no system library, no font to install.

use base64::Engine as _;

pub type Rgba = [u8; 4];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    Start,
    Middle,
    End,
}

#[derive(Clone, Debug)]
pub enum Item {
    /// a rectangle: fill and/or a border of `width`
    Rect { x: f64, y: f64, w: f64, h: f64, fill: Option<Rgba>, stroke: Option<(Rgba, f64)> },
    /// an open polyline or, with two points, a segment
    Line { pts: Vec<[f64; 2]>, color: Rgba, width: f64 },
    /// a filled triangle (an arrow head)
    Tri { pts: [[f64; 2]; 3], fill: Rgba },
    /// text with its baseline at `y`
    Text { x: f64, y: f64, text: String, size: f64, color: Rgba, anchor: Anchor, bold: bool },
    /// a raster image of `px` by `py` pixels, RGBA, drawn in the rectangle
    Image { x: f64, y: f64, w: f64, h: f64, px: u32, py: u32, rgba: Vec<u8> },
}

#[derive(Clone, Debug)]
pub struct Canvas {
    pub width: f64,
    pub height: f64,
    pub background: Rgba,
    pub items: Vec<Item>,
}

pub const WHITE: Rgba = [255, 255, 255, 255];
pub const INK: Rgba = [24, 28, 36, 255];
pub const MUTED: Rgba = [96, 104, 118, 255];

impl Canvas {
    pub fn new(width: f64, height: f64) -> Canvas {
        Canvas { width, height, background: WHITE, items: vec![] }
    }

    pub fn push(&mut self, item: Item) {
        self.items.push(item);
    }

    pub fn text(&mut self, x: f64, y: f64, text: &str, size: f64, color: Rgba, anchor: Anchor) {
        self.items.push(Item::Text { x, y, text: text.to_string(), size, color, anchor, bold: false });
    }

    pub fn bold(&mut self, x: f64, y: f64, text: &str, size: f64, color: Rgba, anchor: Anchor) {
        self.items.push(Item::Text { x, y, text: text.to_string(), size, color, anchor, bold: true });
    }

    /// The figure as SVG text. Fields are embedded as PNG images.
    pub fn to_svg(&self) -> Result<String, String> {
        use std::fmt::Write;
        let mut s = String::new();
        let css = |c: Rgba| format!("rgb({},{},{})", c[0], c[1], c[2]);
        let op =
            |c: Rgba| if c[3] == 255 { String::new() } else { format!(" opacity=\"{:.3}\"", f64::from(c[3]) / 255.0) };
        let w = |r: &mut String, a: std::fmt::Arguments<'_>| {
            let _ = r.write_fmt(a);
        };
        w(&mut s, format_args!("<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"{0}\" height=\"{1}\" viewBox=\"0 0 {0} {1}\">\n", self.width, self.height));
        w(&mut s, format_args!("<rect width=\"100%\" height=\"100%\" fill=\"{}\"/>\n", css(self.background)));
        for it in &self.items {
            match it {
                Item::Rect { x, y, w: rw, h, fill, stroke } => {
                    let f = fill.map_or("none".to_string(), css);
                    let (sc, sw) = stroke.map_or(("none".to_string(), 0.0), |(c, wd)| (css(c), wd));
                    w(&mut s, format_args!("<rect x=\"{x:.2}\" y=\"{y:.2}\" width=\"{rw:.2}\" height=\"{h:.2}\" fill=\"{f}\" stroke=\"{sc}\" stroke-width=\"{sw}\"/>\n"));
                }
                Item::Line { pts, color, width } => {
                    let p: Vec<String> = pts.iter().map(|q| format!("{:.2},{:.2}", q[0], q[1])).collect();
                    w(&mut s, format_args!("<polyline points=\"{}\" fill=\"none\" stroke=\"{}\"{} stroke-width=\"{width}\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n", p.join(" "), css(*color), op(*color)));
                }
                Item::Tri { pts, fill } => {
                    w(
                        &mut s,
                        format_args!(
                            "<polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{}\"{}/>\n",
                            pts[0][0],
                            pts[0][1],
                            pts[1][0],
                            pts[1][1],
                            pts[2][0],
                            pts[2][1],
                            css(*fill),
                            op(*fill)
                        ),
                    );
                }
                Item::Text { x, y, text, size, color, anchor, bold } => {
                    let a = match anchor {
                        Anchor::Start => "start",
                        Anchor::Middle => "middle",
                        Anchor::End => "end",
                    };
                    let esc = text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
                    w(&mut s, format_args!("<text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"Ubuntu, 'Segoe UI', Helvetica, Arial, sans-serif\" font-size=\"{size}\"{} fill=\"{}\" text-anchor=\"{a}\">{esc}</text>\n", if *bold { " font-weight=\"bold\"" } else { "" }, css(*color)));
                }
                Item::Image { x, y, w: iw, h, px, py, rgba } => {
                    let png = encode_png(*px, *py, rgba)?;
                    let b64 = base64::engine::general_purpose::STANDARD.encode(png);
                    w(&mut s, format_args!("<image x=\"{x:.2}\" y=\"{y:.2}\" width=\"{iw:.2}\" height=\"{h:.2}\" preserveAspectRatio=\"none\" xlink:href=\"data:image/png;base64,{b64}\"/>\n"));
                }
            }
        }
        s.push_str("</svg>\n");
        Ok(s)
    }

    /// The figure as PNG bytes, `scale` pixels per unit of the canvas (2 gives a sharp figure on a screen).
    pub fn to_png(&self, scale: f64) -> Result<Vec<u8>, String> {
        let (w, h, rgba) = crate::raster::render(self, scale)?;
        encode_png(w, h, &rgba)
    }
}

/// A PNG of RGBA pixels.
pub fn encode_png(w: u32, h: u32, rgba: &[u8]) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut out, w, h);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        let mut wr = enc.write_header().map_err(|e| format!("PNG: {e}"))?;
        wr.write_image_data(rgba).map_err(|e| format!("PNG: {e}"))?;
    }
    Ok(out)
}
