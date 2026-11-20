// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Figures of the modes of a result (ADR 99): a gallery of modes (`|u|` on the deformed section, the outline of the regions, arrows of the in-plane
//! motion) and the component view of chosen modes (`|u|`, the in-plane components, the propagation component in quadrature). Every mode is put in
//! the phase of `modeview` first. Nothing here names a physics: a vector field is found by the names of its cell fields, the title of a mode
//! comes from its label and from the scalars the study wrote for it, the outline comes from the cell labels of the mesh.

use crate::canvas::{Anchor, Canvas, Item, Rgba, INK, MUTED};
use crate::modeview::{align, AlignOptions, AlignedMode};
use crate::raster::{raster_field, text_width, ColouredTri};
use crate::section::Section;
use nyx_physics::{discretise, Discretisation, Error, Mesh, Scalar, SolveResult};

/// Colour maps of the figures. The tables are those of the window (`crates/app/src/colormap.rs`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cmap {
    Viridis,
    Cividis,
    Inferno,
    Coolwarm,
}

impl Cmap {
    pub fn parse(s: &str) -> Option<Cmap> {
        match s.to_ascii_lowercase().as_str() {
            "viridis" => Some(Cmap::Viridis),
            "cividis" => Some(Cmap::Cividis),
            "inferno" => Some(Cmap::Inferno),
            "coolwarm" | "cool-warm" => Some(Cmap::Coolwarm),
            _ => None,
        }
    }

    pub fn at(self, t: f64) -> [f32; 3] {
        let t = if t.is_finite() { t.clamp(0.0, 1.0) as f32 } else { 0.0 };
        let table: &[[f32; 3]] = match self {
            Cmap::Viridis => &VIRIDIS,
            Cmap::Cividis => &CIVIDIS,
            Cmap::Inferno => &INFERNO,
            Cmap::Coolwarm => &COOLWARM,
        };
        let x = t * (table.len() - 1) as f32;
        let i = (x.floor() as usize).min(table.len() - 2);
        let f = x - i as f32;
        let (a, b) = (table[i], table[i + 1]);
        [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]
    }
}

#[allow(clippy::approx_constant)]
const VIRIDIS: [[f32; 3]; 11] = [
    [0.267, 0.005, 0.329],
    [0.283, 0.141, 0.458],
    [0.254, 0.265, 0.530],
    [0.207, 0.372, 0.553],
    [0.164, 0.471, 0.558],
    [0.128, 0.567, 0.551],
    [0.135, 0.659, 0.518],
    [0.267, 0.749, 0.441],
    [0.478, 0.821, 0.318],
    [0.741, 0.873, 0.150],
    [0.993, 0.906, 0.144],
];
const INFERNO: [[f32; 3]; 11] = [
    [0.001, 0.000, 0.014],
    [0.087, 0.044, 0.224],
    [0.258, 0.039, 0.406],
    [0.416, 0.090, 0.433],
    [0.578, 0.148, 0.404],
    [0.735, 0.216, 0.330],
    [0.865, 0.317, 0.226],
    [0.954, 0.469, 0.103],
    [0.988, 0.652, 0.040],
    [0.964, 0.859, 0.264],
    [0.988, 0.998, 0.645],
];
const COOLWARM: [[f32; 3]; 5] =
    [[0.230, 0.299, 0.754], [0.552, 0.691, 0.996], [0.866, 0.866, 0.866], [0.956, 0.604, 0.486], [0.706, 0.016, 0.150]];
const CIVIDIS: [[f32; 3]; 11] = [
    [0.000000, 0.135112, 0.304751],
    [0.046205, 0.203903, 0.440196],
    [0.208926, 0.272546, 0.424809],
    [0.313287, 0.343120, 0.425120],
    [0.401418, 0.411790, 0.440708],
    [0.488697, 0.485318, 0.471008],
    [0.582087, 0.558670, 0.468118],
    [0.683950, 0.638793, 0.444251],
    [0.785965, 0.720438, 0.399613],
    [0.896818, 0.811030, 0.320982],
    [0.995737, 0.909344, 0.217772],
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// every chosen mode as `|u|` with the outline and arrows
    Gallery,
    /// per chosen mode a row: `|u|`, each in-plane component, the propagation component in quadrature
    Components,
}

#[derive(Clone, Debug)]
pub struct FigureOptions {
    pub kind: Kind,
    /// the modes (0-based solution numbers); the first `max_modes` when `None`
    pub modes: Option<Vec<usize>>,
    pub max_modes: usize,
    pub columns: usize,
    /// width of a panel in canvas units
    pub panel_width: f64,
    /// the colour map of `|u|` (default Inferno); the signed components always use Coolwarm
    pub cmap: Cmap,
    /// the largest motion of the deformed shape as a fraction of the window width (0 draws the undeformed section). Default: 0.05 for a
    /// displacement (a vector field named `u`, as `nyx-sms` names it) and 0 for any other field, which is not a motion of the material
    pub warp: Option<f64>,
    pub arrows: bool,
    /// the window (x0, x1, z0, z1) in metres; the extent of the mode when `None`
    pub window: Option<[f64; 4]>,
    /// the plane a 3D mesh is cut on (y in metres); its first plane when `None`
    pub at: Option<f64>,
    pub core: Option<Vec<i32>>,
    pub family: Option<String>,
}

impl Default for FigureOptions {
    fn default() -> Self {
        FigureOptions {
            kind: Kind::Gallery,
            modes: None,
            max_modes: 12,
            columns: 4,
            panel_width: 250.0,
            cmap: Cmap::Inferno,
            warp: None,
            arrows: true,
            window: None,
            at: None,
            core: None,
            family: None,
        }
    }
}

const PANEL_BG: Rgba = [17, 24, 39, 255];
const OUTLINE: Rgba = [148, 163, 184, 230];
const ARROW: Rgba = [34, 211, 238, 235];

fn bytes(c: [f32; 3]) -> Rgba {
    [(c[0] * 255.0 + 0.5) as u8, (c[1] * 255.0 + 0.5) as u8, (c[2] * 255.0 + 0.5) as u8, 255]
}

/// `17.944 GHz` for 1.7944e10 with the unit Hz.
fn si(v: f64, unit: &str) -> String {
    let a = v.abs();
    let (s, p) = if a == 0.0 {
        (1.0, "")
    } else if a >= 1e12 {
        (1e-12, "T")
    } else if a >= 1e9 {
        (1e-9, "G")
    } else if a >= 1e6 {
        (1e-6, "M")
    } else if a >= 1e3 {
        (1e-3, "k")
    } else if a >= 1.0 {
        (1.0, "")
    } else if a >= 1e-3 {
        (1e3, "m")
    } else if a >= 1e-6 {
        (1e6, "u")
    } else {
        (1e9, "n")
    };
    format!("{:.3} {p}{unit}", v * s).replace(".000 ", " ")
}

/// The two lines under a mode's number: what the study wrote for it (its label, the frequency, Q, the shares of its energy).
fn mode_lines(result: &SolveResult, k: usize) -> (String, String) {
    let sol = &result.solutions[k];
    let scalar = |name: &str| {
        result
            .scalars
            .iter()
            .find(|(n, _)| n == name)
            .and_then(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
    };
    let n = k + 1;
    let mut head = format!("mode {n}");
    let hz = sol.label.contains("Hz");
    if let (Some(re), Some(im)) = (scalar(&format!("f re {n} (Hz)")), scalar(&format!("f im {n} (Hz)"))) {
        head = format!("{head}: {}", si(re, "Hz"));
        if im.abs() > 0.0 {
            head = format!("{head}, Q {:.1}", re / (2.0 * im.abs()));
        }
    } else if let (true, Some(v)) = (hz, sol.value) {
        head = format!("{head}: {}", si(v, "Hz"));
    } else if let Some(e) = sol.eigenvalue {
        let name = if result.study == "modes" { "n_eff = " } else { "" };
        head = format!(
            "{head}: {name}{:.5}{}",
            e[0],
            if e[1].abs() > 1e-9 { format!(" {:+.2e}i", e[1]) } else { String::new() }
        );
    } else {
        head = sol.label.clone();
    }
    // the shares of the energy, as the study wrote them: `<what> share N` and `share in <region> N`
    let suffix = format!(" {n}");
    let mut shares: Vec<(String, f64)> = result
        .scalars
        .iter()
        .filter_map(|(name, v)| {
            let base = name.strip_suffix(&suffix)?;
            let Scalar::Num(x) = v else { return None };
            (base.ends_with("share") || base.starts_with("share in")).then(|| (base.replace(" share", ""), *x))
        })
        .collect();
    shares.sort_by(|a, b| b.1.total_cmp(&a.1));
    let line: Vec<String> =
        shares.iter().take(3).filter(|(_, x)| *x >= 0.005).map(|(b, x)| format!("{b} {:.0}%", 100.0 * x)).collect();
    (head, line.join("   "))
}

/// Everything a panel needs of one mode.
struct Prepared {
    k: usize,
    mode: AlignedMode,
}

fn extent(section: &Section, mode: &AlignedMode, level: f64) -> Option<[f64; 4]> {
    let mag = mode.magnitude();
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    for (i, p) in section.xz.iter().enumerate() {
        if mag[i] >= level {
            b = [b[0].min(p[0]), b[1].max(p[0]), b[2].min(p[1]), b[3].max(p[1])];
        }
    }
    (b[0] <= b[1] && b[2] <= b[3]).then_some(b)
}

struct Frame {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    win: [f64; 4],
}

impl Frame {
    fn map(&self, p: [f64; 2]) -> [f64; 2] {
        [
            self.x + (p[0] - self.win[0]) / (self.win[1] - self.win[0]) * self.w,
            self.y + (self.win[3] - p[1]) / (self.win[3] - self.win[2]) * self.h,
        ]
    }
}

/// What was done to one mode to draw it: the numbers that make the picture reproducible.
#[derive(Clone, Debug)]
pub struct ModeNote {
    /// 1-based, as the study numbers its modes
    pub mode: usize,
    /// the vector field drawn (`u`, `E`)
    pub field: String,
    /// the angle the mode was rotated by, degrees
    pub phase_deg: f64,
    pub flipped: bool,
    /// the peak of `|u|` before scaling
    pub peak: f64,
    /// the share of the propagation component's energy that is in `Re` after the phase step (near 0: in quadrature)
    pub quadrature: Option<f64>,
    /// whether the field has an imaginary part (a real field is drawn as the solver kept it)
    pub complex: bool,
}

#[derive(Clone, Debug)]
pub struct Figure {
    pub canvas: Canvas,
    pub modes: Vec<ModeNote>,
    /// the window drawn (x0, x1, z0, z1), metres
    pub window: [f64; 4],
    /// the plane of the section of a 3D mesh
    pub at: Option<f64>,
}

/// The figure of the modes of `result`, solved on `mesh`.
pub fn mode_figure(mesh: &Mesh, result: &SolveResult, o: &FigureOptions) -> Result<Figure, Error> {
    let disc = discretise(mesh)?;
    mode_figure_on(&disc, result, o)
}

/// The same on the solver mesh.
pub fn mode_figure_on(disc: &Discretisation, result: &SolveResult, o: &FigureOptions) -> Result<Figure, Error> {
    if result.solutions.is_empty() {
        return Err(Error::invalid("the result has no solutions to draw"));
    }
    let section = Section::new(disc, result.order as usize, o.at).map_err(Error::invalid)?;
    let picks: Vec<usize> = match &o.modes {
        Some(m) => m.clone(),
        None => (0..result.solutions.len().min(o.max_modes.max(1))).collect(),
    };
    if let Some(&bad) = picks.iter().find(|&&k| k >= result.solutions.len()) {
        return Err(Error::invalid(format!(
            "mode {} is not in the result, which has {}",
            bad + 1,
            result.solutions.len()
        )));
    }
    let opts = AlignOptions { family: o.family.clone(), core: o.core.clone(), ..Default::default() };
    let mut prep = vec![];
    for &k in &picks {
        prep.push(Prepared { k, mode: align(&section, &result.solutions[k], &opts)? });
    }
    // The window comes from the data: the places where a mode is above 12 % of its peak, with a margin. Modes of similar extent share one window
    // (so that they can be compared); a mode more than twice as wide as the narrowest (a slab mode next to the guided ones) gets its own.
    let exts: Vec<Option<[f64; 4]>> = prep.iter().map(|p| extent(&section, &p.mode, 0.12)).collect();
    let narrowest = exts.iter().flatten().map(|e| e[1] - e[0]).fold(f64::INFINITY, f64::min);
    let shared = |e: &Option<[f64; 4]>| e.is_some_and(|e| e[1] - e[0] <= 2.0 * narrowest);
    let win = match o.window {
        Some(w) => w,
        None => {
            let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
            for e in exts.iter().filter(|e| shared(e)).flatten() {
                b = [b[0].min(e[0]), b[1].max(e[1]), b[2].min(e[2]), b[3].max(e[3])];
            }
            if !(b[0] < b[1] && b[2] < b[3]) {
                let xs = section.xz.iter().map(|p| p[0]);
                let zs = section.xz.iter().map(|p| p[1]);
                b = [
                    xs.clone().fold(f64::INFINITY, f64::min),
                    xs.fold(f64::NEG_INFINITY, f64::max),
                    zs.clone().fold(f64::INFINITY, f64::min),
                    zs.fold(f64::NEG_INFINITY, f64::max),
                ];
            }
            let m = 0.2 * (b[1] - b[0]).max(b[3] - b[2]).max(1e-12);
            [b[0] - m, b[1] + m, b[2] - m, b[3] + m]
        }
    };
    if !(win[1] > win[0] && win[3] > win[2]) {
        return Err(Error::invalid("the window of the figure is empty (x0 < x1 and z0 < z1 are needed)"));
    }
    let aspect = ((win[3] - win[2]) / (win[1] - win[0])).clamp(0.25, 1.6);
    let (pw, ph) = (o.panel_width, o.panel_width * aspect);
    let win = fit_aspect(win, aspect);
    // the window of mode `i`
    let window_of = |i: usize| -> [f64; 4] {
        match (&o.window, &exts[i]) {
            (None, Some(e)) if !shared(&exts[i]) => {
                let m = 0.2 * (e[1] - e[0]).max(e[3] - e[2]).max(1e-12);
                fit_aspect([e[0] - m, e[1] + m, e[2] - m, e[3] + m], aspect)
            }
            _ => win,
        }
    };
    let tris = section.linear_tris();
    let outline = section.outline();
    let columns_per_mode = match o.kind {
        Kind::Gallery => 1,
        Kind::Components => 1 + prep[0].mode.in_plane.len().min(2) + usize::from(prep[0].mode.along.is_some()),
    };
    // the spectrum of the modes goes in the next free panel of a gallery
    let freqs: Vec<Option<f64>> = prep.iter().map(|p| mode_frequency(result, p.k)).collect();
    let spectrum = o.kind == Kind::Gallery && prep.len() >= 3 && freqs.iter().all(Option::is_some);
    let n_panels = prep.len() + usize::from(spectrum);
    let (ncols, nrows) = match o.kind {
        Kind::Gallery => (o.columns.clamp(1, 8).min(n_panels), n_panels.div_ceil(o.columns.clamp(1, 8).min(n_panels))),
        Kind::Components => (columns_per_mode, prep.len()),
    };
    // layout: the tick labels take `lab` on the left of every panel
    let head_h = if o.kind == Kind::Components { 64.0 } else { 46.0 };
    let (left, lab, gap, top_title, foot) = (14.0, 40.0, 14.0, 54.0, 70.0);
    let stride = lab + pw + gap;
    let (note_a, note_b) = (key_text(&prep[0].mode.prefix, "um"), footer_text(&prep[0].mode.prefix));
    let width = (left + ncols as f64 * stride + 6.0)
        .max(left + lab + 180.0 + 12.0 + text_width(&note_a, 11.0) + 10.0)
        .max(left + lab + text_width(&note_b, 11.0) + 10.0);
    let row_h = head_h + ph + 30.0;
    let height = top_title + nrows as f64 * row_h + foot;
    let mut cv = Canvas::new(width, height);
    let place = match section.at {
        Some(y) => format!("section of the 3D mesh at y = {}", si(y, "m")),
        None => "cross-section".to_string(),
    };
    cv.bold(left, 30.0, &format!("{}: {} ({place})", result.physics, result.study), 18.0, INK, Anchor::Start);
    let (unit, uname) = length_unit(win[1] - win[0]);
    let warp_frac = o.warp.unwrap_or(if prep[0].mode.prefix == "u" { 0.07 } else { 0.0 });
    let axis_of = |m: &AlignedMode, a: char| m.axes.iter().position(|&c| c == a);
    for (idx, p) in prep.iter().enumerate() {
        let (r, c0) = match o.kind {
            Kind::Gallery => (idx / ncols, idx % ncols),
            Kind::Components => (idx, 0),
        };
        let m = &p.mode;
        let mw = window_of(idx);
        let warp_len = warp_frac * (mw[1] - mw[0]);
        let (head, sub) = mode_lines(result, p.k);
        let (ix, iz) = (axis_of(m, 'x'), axis_of(m, 'z'));
        let shift = |i: usize| -> [f64; 2] {
            let q = section.xz[i];
            match (ix, iz) {
                (Some(a), Some(b)) if warp_len > 0.0 => [q[0] + warp_len * m.re[a][i], q[1] + warp_len * m.re[b][i]],
                _ => q,
            }
        };
        let pos: Vec<[f64; 2]> = (0..section.n_nodes()).map(shift).collect();
        let mag = m.magnitude();
        struct Col {
            title: String,
            values: Vec<f64>,
            signed: Option<f64>,
            arrows: bool,
        }
        let mut cols: Vec<Col> = vec![Col {
            title: if o.kind == Kind::Components { format!("|{}|", m.prefix) } else { String::new() },
            values: mag.clone(),
            signed: None,
            arrows: o.arrows,
        }];
        if o.kind == Kind::Components {
            for &c in m.in_plane.iter().take(2) {
                let v = m.re[c].clone();
                let mx = v.iter().fold(0.0f64, |a, x| a.max(x.abs())).max(1e-12);
                cols.push(Col {
                    title: format!("Re {}{}", m.prefix, m.axes[c]),
                    values: v,
                    signed: Some(mx),
                    arrows: false,
                });
            }
            if let Some(c) = m.along {
                let (v, part) = m.quadrature_part(c);
                let v = v.to_vec();
                let mx = v.iter().fold(0.0f64, |a, x| a.max(x.abs())).max(1e-12);
                // a real field has no Re or Im: it is the amplitude the solver kept
                let title = if m.complex {
                    format!("{part} {}{}", m.prefix, m.axes[c])
                } else {
                    format!("{}{}", m.prefix, m.axes[c])
                };
                cols.push(Col { title, values: v, signed: Some(mx), arrows: false });
            }
        }
        for (ci, col) in cols.iter().enumerate() {
            let x = left + lab + (c0 + ci) as f64 * stride;
            let y = top_title + r as f64 * row_h + head_h;
            let fr = Frame { x, y, w: pw, h: ph, win: mw };
            cv.push(Item::Rect { x, y, w: pw, h: ph, fill: Some(PANEL_BG), stroke: None });
            // the field, as an image of twice the panel in pixels
            let (iw, ih) = ((pw * 2.0).ceil() as u32, (ph * 2.0).ceil() as u32);
            let to_img = |q: [f64; 2]| -> [f64; 2] {
                let f = fr.map(q);
                [(f[0] - x) * 2.0, (f[1] - y) * 2.0]
            };
            let colour = |v: f64| match col.signed {
                None => o.cmap.at(v),
                Some(mx) => Cmap::Coolwarm.at(0.5 + 0.5 * v / mx),
            };
            let img_tris: Vec<ColouredTri> = tris
                .iter()
                .map(|t| {
                    let n = [t[0] as usize, t[1] as usize, t[2] as usize];
                    (
                        [to_img(pos[n[0]]), to_img(pos[n[1]]), to_img(pos[n[2]])],
                        [colour(col.values[n[0]]), colour(col.values[n[1]]), colour(col.values[n[2]])],
                    )
                })
                .collect();
            cv.push(Item::Image { x, y, w: pw, h: ph, px: iw, py: ih, rgba: raster_field(iw, ih, &img_tris) });
            // the undeformed section: the borders of its regions, from the cell labels of the mesh
            for seg in &outline {
                if let Some(pts) = clip(
                    fr.map(section.xz[seg[0] as usize]),
                    fr.map(section.xz[seg[1] as usize]),
                    [x, y, x + pw, y + ph],
                ) {
                    cv.push(Item::Line { pts: pts.to_vec(), color: OUTLINE, width: 1.0 });
                }
            }
            if col.arrows {
                if let (Some(a), Some(b)) = (ix, iz) {
                    let extent_w = exts[idx].map_or(mw[1] - mw[0], |e| e[1] - e[0]).min(mw[1] - mw[0]);
                    // about nine arrows across the mode, never closer than 1/40 of the window
                    let spacing = (extent_w / 9.0).max((mw[1] - mw[0]) / 40.0);
                    draw_arrows(&mut cv, &fr, &section, &tris, (&m.re[a], &m.re[b]), spacing, warp_len);
                }
            }
            cv.push(Item::Rect { x, y, w: pw, h: ph, fill: None, stroke: Some(([71, 85, 105, 255], 0.8)) });
            draw_ticks(&mut cv, &fr, unit);
            if let Some(mx) = col.signed {
                cv.text(x + pw - 4.0, y + 12.0, &format!("+-{mx:.2}"), 10.5, INK, Anchor::End);
            }
            if o.kind == Kind::Components {
                cv.bold(x + pw / 2.0, y - 8.0, &col.title, 12.5, INK, Anchor::Middle);
            }
        }
        // the heading of the mode (a gallery: over its panel; components: over the row)
        let hx = left + lab + c0 as f64 * stride;
        let hy = top_title + r as f64 * row_h;
        let head_y = if o.kind == Kind::Components { hy + 10.0 } else { hy + 16.0 };
        cv.bold(hx, head_y, &head, 14.0, INK, Anchor::Start);
        let mut note = sub;
        if let Some(q) = m.quadrature {
            if q > 0.05 {
                note = format!(
                    "{note}{}Re of the propagation component {:.0}% (not in quadrature)",
                    if note.is_empty() { "" } else { "   " },
                    100.0 * q
                );
            }
        }
        if !note.is_empty() {
            cv.text(hx, head_y + 15.0, &note, 11.0, MUTED, Anchor::Start);
        }
    }
    if spectrum {
        let (r, c0) = (prep.len() / ncols, prep.len() % ncols);
        let x = left + lab + c0 as f64 * stride;
        let y = top_title + r as f64 * row_h + head_h;
        let numbers: Vec<usize> = prep.iter().map(|p| p.k + 1).collect();
        let f: Vec<f64> = freqs.iter().map(|f| f.unwrap_or(0.0)).collect();
        draw_spectrum(&mut cv, x, y, pw, ph, &numbers, &f);
        cv.bold(x, top_title + r as f64 * row_h + 16.0, "Spectrum of the modes drawn", 14.0, INK, Anchor::Start);
    }
    // the key and what was done to the modes
    let ky = top_title + nrows as f64 * row_h + 14.0;
    let kw = 180.0;
    let steps = 40;
    for i in 0..steps {
        let t = i as f64 / (steps - 1) as f64;
        cv.push(Item::Rect {
            x: left + lab + i as f64 * kw / steps as f64,
            y: ky,
            w: kw / steps as f64 + 0.6,
            h: 10.0,
            fill: Some(bytes(o.cmap.at(t))),
            stroke: None,
        });
    }
    cv.text(left + lab, ky + 24.0, "0", 10.5, MUTED, Anchor::Start);
    cv.text(left + lab + kw, ky + 24.0, "1", 10.5, MUTED, Anchor::End);
    cv.text(left + lab + kw + 12.0, ky + 9.0, &key_text(&prep[0].mode.prefix, uname), 11.0, MUTED, Anchor::Start);
    cv.text(left + lab, ky + 44.0, &footer_text(&prep[0].mode.prefix), 11.0, MUTED, Anchor::Start);
    let modes = prep
        .iter()
        .map(|p| ModeNote {
            mode: p.k + 1,
            field: p.mode.prefix.clone(),
            phase_deg: p.mode.theta.to_degrees(),
            flipped: p.mode.flipped,
            peak: p.mode.peak,
            quadrature: p.mode.quadrature,
            complex: p.mode.complex,
        })
        .collect();
    Ok(Figure { canvas: cv, modes, window: win, at: section.at })
}

fn key_text(prefix: &str, uname: &str) -> String {
    let deform = if prefix == "u" { " The section is deformed by Re u." } else { "" };
    format!("|{prefix}| / its peak. Grey lines: the undeformed regions (axes in {uname}).{deform}")
}

fn footer_text(prefix: &str) -> String {
    format!("Each mode is turned to the phase where its in-plane motion is largest in Re {prefix} (sign: largest component positive) and scaled to a peak of 1.")
}

/// The frequency of mode `k` in Hz when the study states one.
fn mode_frequency(result: &SolveResult, k: usize) -> Option<f64> {
    let n = k + 1;
    let sc = |name: String| {
        result
            .scalars
            .iter()
            .find(|(m, _)| *m == name)
            .and_then(|(_, v)| if let Scalar::Num(x) = v { Some(*x) } else { None })
    };
    sc(format!("f re {n} (Hz)")).or_else(|| sc(format!("f{n} (Hz)"))).or_else(|| {
        let s = &result.solutions[k];
        if s.label.contains("Hz") {
            s.value
        } else {
            None
        }
    })
}

/// A unit of length for a window `w` metres wide: the factor to metres and its name.
fn length_unit(w: f64) -> (f64, &'static str) {
    if w >= 1e-3 {
        (1e-3, "mm")
    } else if w >= 2e-7 {
        (1e-6, "\u{b5}m")
    } else {
        (1e-9, "nm")
    }
}

/// Round numbers (1, 2, 5 times a power of ten) inside `lo..=hi`, about `n` of them.
fn nice_ticks(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let raw = (hi - lo) / n.max(1) as f64;
    if raw <= 0.0 || !raw.is_finite() {
        return vec![];
    }
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0].iter().map(|k| k * mag).find(|s| *s >= raw).unwrap_or(10.0 * mag);
    let first = (lo / step).ceil();
    (0..).map(|i| (first + f64::from(i)) * step).take_while(|v| *v <= hi + 1e-9 * step).collect()
}

fn tick_text(v: f64, step: f64) -> String {
    let s = if step >= 1.0 {
        format!("{v:.0}")
    } else if step >= 0.1 {
        format!("{v:.1}")
    } else {
        format!("{v:.2}")
    };
    match s.strip_prefix('-') {
        Some(rest) if rest.chars().all(|c| c == '0' || c == '.') => rest.to_string(),
        _ => s,
    }
}

/// Tick marks and numbers on the left and bottom edge of a panel, in `unit` metres.
fn draw_ticks(cv: &mut Canvas, fr: &Frame, unit: f64) {
    let (x0, x1, z0, z1) = (fr.win[0] / unit, fr.win[1] / unit, fr.win[2] / unit, fr.win[3] / unit);
    let tick = [71, 85, 105, 255];
    let xs = nice_ticks(x0, x1, 5);
    let step_x = if xs.len() > 1 { xs[1] - xs[0] } else { 1.0 };
    for v in xs {
        let px = fr.x + (v - x0) / (x1 - x0) * fr.w;
        cv.push(Item::Line { pts: vec![[px, fr.y + fr.h], [px, fr.y + fr.h + 3.5]], color: tick, width: 0.8 });
        cv.text(px, fr.y + fr.h + 15.0, &tick_text(v, step_x), 10.0, MUTED, Anchor::Middle);
    }
    let zs = nice_ticks(z0, z1, 4);
    let step_z = if zs.len() > 1 { zs[1] - zs[0] } else { 1.0 };
    for v in zs {
        let py = fr.y + (z1 - v) / (z1 - z0) * fr.h;
        cv.push(Item::Line { pts: vec![[fr.x - 3.5, py], [fr.x, py]], color: tick, width: 0.8 });
        cv.text(fr.x - 6.0, py + 3.5, &tick_text(v, step_z), 10.0, MUTED, Anchor::End);
    }
}

fn dot(cv: &mut Canvas, x: f64, y: f64, r: f64, c: Rgba) {
    let n = 12;
    for i in 0..n {
        let (a, b) = (
            std::f64::consts::TAU * f64::from(i) / f64::from(n),
            std::f64::consts::TAU * f64::from(i + 1) / f64::from(n),
        );
        cv.push(Item::Tri {
            pts: [[x, y], [x + r * a.cos(), y + r * a.sin()], [x + r * b.cos(), y + r * b.sin()]],
            fill: c,
        });
    }
}

/// Frequency against the number of the mode, in the panel at (x, y).
fn draw_spectrum(cv: &mut Canvas, x: f64, y: f64, w: f64, h: f64, numbers: &[usize], f: &[f64]) {
    cv.push(Item::Rect { x, y, w, h, fill: Some([248, 250, 252, 255]), stroke: Some(([71, 85, 105, 255], 0.8)) });
    let (lo, hi) = f.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
    let (scale, name) = if hi >= 1e9 {
        (1e-9, "GHz")
    } else if hi >= 1e6 {
        (1e-6, "MHz")
    } else {
        (1.0, "Hz")
    };
    let pad = 0.12 * (hi - lo).max(hi * 1e-3);
    let (a, b) = ((lo - pad) * scale, (hi + pad) * scale);
    let ticks = nice_ticks(a, b, 4);
    let step = if ticks.len() > 1 { ticks[1] - ticks[0] } else { 1.0 };
    for t in ticks {
        let py = y + (b - t) / (b - a) * h;
        cv.push(Item::Line { pts: vec![[x, py], [x + w, py]], color: [226, 232, 240, 255], width: 0.7 });
        cv.text(x - 6.0, py + 3.5, &tick_text(t, step), 10.0, MUTED, Anchor::End);
    }
    cv.text(x + 6.0, y + 12.0, name, 10.5, MUTED, Anchor::Start);
    let n = numbers.len() as f64;
    for (i, (num, v)) in numbers.iter().zip(f).enumerate() {
        let px = x + (i as f64 + 0.5) / n * w;
        let py = y + (b - v * scale) / (b - a) * h;
        dot(cv, px, py, 4.2, [234, 88, 12, 255]);
        cv.text(px, y + h + 15.0, &num.to_string(), 10.0, MUTED, Anchor::Middle);
    }
}

/// `win` grown around its centre to the shape `aspect` (height over width) of a panel, so that the scale is the same in x and z.
fn fit_aspect(win: [f64; 4], aspect: f64) -> [f64; 4] {
    let (w, h) = (win[1] - win[0], win[3] - win[2]);
    let (cx, cz) = (0.5 * (win[0] + win[1]), 0.5 * (win[2] + win[3]));
    if h / w < aspect {
        let h2 = w * aspect;
        [win[0], win[1], cz - 0.5 * h2, cz + 0.5 * h2]
    } else {
        let w2 = h / aspect;
        [cx - 0.5 * w2, cx + 0.5 * w2, win[2], win[3]]
    }
}

/// The part of the segment `a b` inside `rect` (x0, y0, x1, y1), or `None` when it lies outside (Liang-Barsky).
fn clip(a: [f64; 2], b: [f64; 2], rect: [f64; 4]) -> Option<[[f64; 2]; 2]> {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (p, q) in [(-dx, a[0] - rect[0]), (dx, rect[2] - a[0]), (-dy, a[1] - rect[1]), (dy, rect[3] - a[1])] {
        if p == 0.0 {
            if q < 0.0 {
                return None;
            }
        } else {
            let r = q / p;
            if p < 0.0 {
                t0 = t0.max(r);
            } else {
                t1 = t1.min(r);
            }
            if t0 > t1 {
                return None;
            }
        }
    }
    Some([[a[0] + t0 * dx, a[1] + t0 * dy], [a[0] + t1 * dx, a[1] + t1 * dy]])
}

/// Arrows of the in-plane motion Re (ux, uz) on a regular grid of `spacing` (metres) over the section: the field is read where a grid point
/// lies inside a triangle (linear between the nodes of the drawn triangles), so the arrows are evenly spaced whatever the mesh is. They are
/// centred on the grid point moved by the deformation, as long as the motion is at least 12 % of the largest one of the panel.
fn draw_arrows(
    cv: &mut Canvas,
    fr: &Frame,
    section: &Section,
    lin: &[[u32; 3]],
    u: (&[f64], &[f64]),
    spacing: f64,
    warp_len: f64,
) {
    use std::collections::BTreeMap;
    let (ux, uz) = u;
    let xz = &section.xz;
    let (ox, oz) = (fr.win[0], fr.win[2]);
    let mut samples: BTreeMap<(i64, i64), [f64; 4]> = BTreeMap::new();
    for t in lin {
        let (n0, n1, n2) = (t[0] as usize, t[1] as usize, t[2] as usize);
        let (a, b, c) = (xz[n0], xz[n1], xz[n2]);
        let area = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
        if area.abs() < 1e-40 {
            continue;
        }
        let (x0, x1) = (a[0].min(b[0]).min(c[0]), a[0].max(b[0]).max(c[0]));
        let (z0, z1) = (a[1].min(b[1]).min(c[1]), a[1].max(b[1]).max(c[1]));
        let (i0, i1) = (((x0 - ox) / spacing - 0.5).ceil() as i64, ((x1 - ox) / spacing - 0.5).floor() as i64);
        let (j0, j1) = (((z0 - oz) / spacing - 0.5).ceil() as i64, ((z1 - oz) / spacing - 0.5).floor() as i64);
        for j in j0..=j1 {
            for i in i0..=i1 {
                let p = [ox + (i as f64 + 0.5) * spacing, oz + (j as f64 + 0.5) * spacing];
                let w1 = ((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])) / area;
                let w2 = ((b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])) / area;
                let w0 = 1.0 - w1 - w2;
                if w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9 {
                    continue;
                }
                samples.entry((i, j)).or_insert([
                    p[0],
                    p[1],
                    w0 * ux[n0] + w1 * ux[n1] + w2 * ux[n2],
                    w0 * uz[n0] + w1 * uz[n1] + w2 * uz[n2],
                ]);
            }
        }
    }
    let biggest = samples.values().map(|s| s[2].hypot(s[3])).fold(0.0f64, f64::max);
    if biggest <= 0.0 {
        return;
    }
    let sp_px = spacing / (fr.win[1] - fr.win[0]) * fr.w;
    for s in samples.values() {
        let m = s[2].hypot(s[3]);
        if m < 0.12 * biggest {
            continue;
        }
        let c = fr.map([s[0] + warp_len * s[2], s[1] + warp_len * s[3]]);
        if c[0] < fr.x || c[0] > fr.x + fr.w || c[1] < fr.y || c[1] > fr.y + fr.h {
            continue;
        }
        // the arrow points where Re u points; the screen y runs down
        let len = 0.85 * sp_px * m / biggest;
        let (ex, ey) = (s[2] / m, -s[3] / m);
        let (a, b) = ([c[0] - 0.5 * len * ex, c[1] - 0.5 * len * ey], [c[0] + 0.5 * len * ex, c[1] + 0.5 * len * ey]);
        cv.push(Item::Line { pts: vec![a, b], color: ARROW, width: 1.3 });
        let hl = (0.35 * len).clamp(3.0, 6.0);
        cv.push(Item::Tri {
            pts: [
                [b[0] + ex * hl * 0.7, b[1] + ey * hl * 0.7],
                [b[0] - ey * hl * 0.45, b[1] + ex * hl * 0.45],
                [b[0] + ey * hl * 0.45, b[1] - ex * hl * 0.45],
            ],
            fill: ARROW,
        });
    }
}
