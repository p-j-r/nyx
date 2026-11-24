// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The colour tokens of the window as plain data (ADR 88): every role in a light and a dark value, and the two region palettes. No egui types
//! here: `theme.rs` turns this into `Visuals`, the tests below check the contrast of every text pair and the colour-blind separation of the
//! region palettes. The values equal `docs/ux/tokens.json` (a test compares the two), which is the document Paul accepted.

/// An sRGB colour, 8 bits per channel.
pub type Rgb = [u8; 3];

/// Which set of colours: the window follows the system, the canvas has its own choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Light,
    Dark,
}

/// A role with one value per mode.
#[derive(Clone, Copy, Debug)]
pub struct Pair {
    pub light: Rgb,
    pub dark: Rgb,
}

impl Pair {
    pub const fn get(&self, mode: Mode) -> Rgb {
        match mode {
            Mode::Light => self.light,
            Mode::Dark => self.dark,
        }
    }
}

const fn digit(b: u8) -> u8 {
    match b {
        b'0'..=b'9' => b - b'0',
        b'a'..=b'f' => b - b'a' + 10,
        b'A'..=b'F' => b - b'A' + 10,
        _ => panic!("not a hex digit"),
    }
}

/// `"#RRGGBB"` to bytes, at compile time, so the table below reads like the JSON.
pub const fn hex(s: &str) -> Rgb {
    let b = s.as_bytes();
    assert!(b.len() == 7 && b[0] == b'#');
    [digit(b[1]) * 16 + digit(b[2]), digit(b[3]) * 16 + digit(b[4]), digit(b[5]) * 16 + digit(b[6])]
}

const fn pair(light: &str, dark: &str) -> Pair {
    Pair { light: hex(light), dark: hex(dark) }
}

/// window background
pub const PAGE: Pair = pair("#F8F6F8", "#141218");
/// cards, sheets, fields
pub const CARD: Pair = pair("#FFFFFF", "#1D1A22");
/// hover and pressed neutral fills, popups
pub const RAISED: Pair = pair("#F1EEF2", "#262230");
/// dividers and card outlines only; never a control edge
pub const SEPARATOR: Pair = pair("#E6E4E6", "#2E2A35");
/// edge of buttons, fields, combo boxes, checkboxes (3:1 on page and card)
pub const CONTROL_EDGE: Pair = pair("#8A8091", "#6E6578");
/// the control edge of the High contrast setting (4.5:1 on page, card and raised)
pub const CONTROL_EDGE_HIGH: Pair = pair("#6B6173", "#9A91A5");
/// primary text
pub const TEXT: Pair = pair("#1D1224", "#ECE7EF");
/// secondary text
pub const MUTED: Pair = pair("#5C5062", "#A9A0B1");
/// the brand colour: primary action fill, links, focus ring, selected bar
pub const ACCENT: Pair = pair("#820C4F", "#D96AA4");
/// text on an accent fill (white on the dark pink fails, ink passes)
pub const ON_ACCENT: Pair = pair("#FFFFFF", "#1D1224");
/// selected item background (with a 2 px accent bar)
pub const TINT: Pair = pair("#F1E4EC", "#33202F");
/// running, queued, stale: fill
pub const ACTIVITY: Pair = pair("#FFBD59", "#F2B35E");
/// text on an activity fill
#[allow(dead_code)] // a role of tokens.json (the text on an amber fill); no widget draws an amber fill since the stale band became a chip (ADR 123)
pub const ON_ACTIVITY: Pair = pair("#1D1224", "#1D1224");
/// the activity colour used as text (on `ACTIVITY_TINT`)
pub const ACTIVITY_TEXT: Pair = pair("#6B4300", "#F2B35E");
/// background of the busy chip and the stale band
pub const ACTIVITY_TINT: Pair = pair("#FFF1DC", "#3D2D1A");
/// saved, converged, passed (always with a glyph and a word)
pub const SUCCESS: Pair = pair("#1F6B43", "#5FC08C");
/// failed, destructive (always with a glyph and a word)
pub const DANGER: Pair = pair("#A32020", "#F27A7A");
/// viewport and plot background (default dark in both themes)
pub const CANVAS: Pair = pair("#F1EFF2", "#181A1E");
/// region number chip background: the canvas colour of the same mode, with `TEXT` on it
pub const CHIP: Pair = pair("#F1EFF2", "#181A1E");
/// fill of regions beyond the sixteenth (slot 17 and up): neutral, told apart by the number chip only
pub const REGION_GREY: Pair = pair("#848890", "#B8BCC4");

/// Every role by its name in `tokens.json`.
#[cfg(test)]
pub const ROLES: [(&str, Pair); 20] = [
    ("page", PAGE),
    ("card", CARD),
    ("raised", RAISED),
    ("separator", SEPARATOR),
    ("control_edge", CONTROL_EDGE),
    ("control_edge_high", CONTROL_EDGE_HIGH),
    ("text", TEXT),
    ("muted", MUTED),
    ("accent", ACCENT),
    ("on_accent", ON_ACCENT),
    ("tint", TINT),
    ("activity", ACTIVITY),
    ("on_activity", ON_ACTIVITY),
    ("activity_text", ACTIVITY_TEXT),
    ("activity_tint", ACTIVITY_TINT),
    ("success", SUCCESS),
    ("danger", DANGER),
    ("canvas", CANVAS),
    ("chip", CHIP),
    ("region_grey", REGION_GREY),
];

/// Eight region colours for the dark canvas (the default canvas, also in a light window). Any two can touch: checked with all pairs.
pub const REGIONS_ON_DARK_CANVAS: [Rgb; 8] =
    [hex("#2B60C2"), hex("#109EC6"), hex("#C83220"), hex("#A68729"), hex("#9D71F3"), hex("#C67294"), hex("#2A7D40"), hex("#9D399B")];

/// Eight region colours for the light canvas (only when Settings says the canvas follows a light theme).
pub const REGIONS_ON_LIGHT_CANVAS: [Rgb; 8] =
    [hex("#B69AF8"), hex("#E64F65"), hex("#9B084D"), hex("#2A81C2"), hex("#BD9D0D"), hex("#1F7C13"), hex("#5F36A1"), hex("#46BEBB")];

/// The region palette for the canvas that is drawn.
pub fn regions(canvas: Mode) -> &'static [Rgb; 8] {
    match canvas {
        Mode::Dark => &REGIONS_ON_DARK_CANVAS,
        Mode::Light => &REGIONS_ON_LIGHT_CANVAS,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- colour maths, the same as the dataviz validator (validate_palette.py): Machado, Oliveira and Fernandes (2009) at severity 1.0,
    // OKLab distance times 100, WCAG contrast. Do not swap the simulation model: the thresholds are calibrated to it.

    fn lin(c: Rgb) -> [f64; 3] {
        c.map(|v| {
            let v = f64::from(v) / 255.0;
            if v <= 0.04045 { v / 12.92 } else { ((v + 0.055) / 1.055).powf(2.4) }
        })
    }

    fn luminance(c: Rgb) -> f64 {
        let [r, g, b] = lin(c);
        0.2126 * r + 0.7152 * g + 0.0722 * b
    }

    fn contrast(a: Rgb, b: Rgb) -> f64 {
        let (x, y) = (luminance(a), luminance(b));
        (x.max(y) + 0.05) / (x.min(y) + 0.05)
    }

    fn oklab([r, g, b]: [f64; 3]) -> [f64; 3] {
        let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
        let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
        let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
        [0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s, 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s, 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s]
    }

    const PROTAN: [[f64; 3]; 3] = [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]];
    const DEUTAN: [[f64; 3]; 3] = [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]];

    fn simulate(c: Rgb, m: &[[f64; 3]; 3]) -> [f64; 3] {
        let v = lin(c);
        m.map(|row| (row[0] * v[0] + row[1] * v[1] + row[2] * v[2]).clamp(0.0, 1.0))
    }

    fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
        100.0 * ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    /// OKLab distance times 100 of two colours; `sim` is a colour-vision simulation or `None` for normal vision.
    fn delta_e(a: Rgb, b: Rgb, sim: Option<&[[f64; 3]; 3]>) -> f64 {
        match sim {
            Some(m) => dist(oklab(simulate(a, m)), oklab(simulate(b, m))),
            None => dist(oklab(lin(a)), oklab(lin(b))),
        }
    }

    fn oklch(c: Rgb) -> (f64, f64) {
        let [l, a, b] = oklab(lin(c));
        (l, a.hypot(b))
    }

    /// (worst colour-blind pair under protan and deutan, worst normal-vision pair) over ALL pairs of the palette.
    fn worst(p: &[Rgb]) -> (f64, f64) {
        let (mut cvd, mut normal) = (f64::INFINITY, f64::INFINITY);
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                for m in [&PROTAN, &DEUTAN] {
                    cvd = cvd.min(delta_e(p[i], p[j], Some(m)));
                }
                normal = normal.min(delta_e(p[i], p[j], None));
            }
        }
        (cvd, normal)
    }

    fn pairs_of(mode: Mode) -> Vec<(&'static str, Rgb, &'static str, Rgb, f64)> {
        let g = |p: Pair| p.get(mode);
        // (foreground role, background role, least contrast): text 7, secondary text and statuses 4.5, controls 3 (WCAG 2.2 AA/AAA)
        let table: [(&str, Pair, &str, Pair, f64); 25] = [
            ("text", TEXT, "page", PAGE, 7.0),
            ("text", TEXT, "card", CARD, 7.0),
            ("text", TEXT, "raised", RAISED, 7.0),
            ("text", TEXT, "tint", TINT, 7.0),
            ("muted", MUTED, "page", PAGE, 4.5),
            ("muted", MUTED, "card", CARD, 4.5),
            ("muted", MUTED, "raised", RAISED, 4.5),
            ("accent", ACCENT, "page", PAGE, 4.5),
            ("accent", ACCENT, "card", CARD, 4.5),
            ("accent", ACCENT, "tint", TINT, 4.5),
            ("on_accent", ON_ACCENT, "accent", ACCENT, 4.5),
            ("on_activity", ON_ACTIVITY, "activity", ACTIVITY, 4.5),
            ("activity_text", ACTIVITY_TEXT, "activity_tint", ACTIVITY_TINT, 4.5),
            ("activity_text", ACTIVITY_TEXT, "page", PAGE, 4.5),
            ("success", SUCCESS, "page", PAGE, 4.5),
            ("success", SUCCESS, "card", CARD, 4.5),
            ("danger", DANGER, "page", PAGE, 4.5),
            ("danger", DANGER, "card", CARD, 4.5),
            ("control_edge", CONTROL_EDGE, "page", PAGE, 3.0),
            ("control_edge", CONTROL_EDGE, "card", CARD, 3.0),
            ("control_edge_high", CONTROL_EDGE_HIGH, "page", PAGE, 4.5),
            ("control_edge_high", CONTROL_EDGE_HIGH, "card", CARD, 4.5),
            ("control_edge_high", CONTROL_EDGE_HIGH, "raised", RAISED, 4.5),
            ("text", TEXT, "chip", CHIP, 7.0),
            ("text", TEXT, "canvas", CANVAS, 7.0),
        ];
        table.iter().map(|(fa, a, fb, b, min)| (*fa, g(*a), *fb, g(*b), *min)).collect()
    }

    #[test]
    fn every_text_and_edge_pair_has_its_contrast() {
        let mut bad = vec![];
        for mode in [Mode::Light, Mode::Dark] {
            for (fa, a, fb, b, min) in pairs_of(mode) {
                let c = contrast(a, b);
                if c < min {
                    bad.push(format!("{mode:?}: {fa} on {fb} is {c:.2}:1, needs {min}"));
                }
            }
        }
        assert!(bad.is_empty(), "{}", bad.join("\n"));
    }

    #[test]
    fn region_palettes_reproduce_the_validator_numbers() {
        // recorded with validate_palette.py --pairs all on 2026-11-17: dark CVD 6.7 and normal 15.8, light CVD 8.8 and normal 18.4
        let (cvd, normal) = worst(&REGIONS_ON_DARK_CANVAS);
        assert!((cvd - 6.7).abs() < 0.1 && (normal - 15.8).abs() < 0.1, "dark: {cvd:.2} {normal:.2}");
        let (cvd, normal) = worst(&REGIONS_ON_LIGHT_CANVAS);
        assert!((cvd - 8.8).abs() < 0.1 && (normal - 18.4).abs() < 0.1, "light: {cvd:.2} {normal:.2}");
    }

    #[test]
    fn any_two_regions_can_be_told_apart() {
        // the floors of the validator: colour-blind 6 (the band from 6 to 8 is legal only with a number chip and a legend, which every region has),
        // normal vision 15; all pairs, because any two regions can touch
        for (name, p) in [("dark", &REGIONS_ON_DARK_CANVAS), ("light", &REGIONS_ON_LIGHT_CANVAS)] {
            let (cvd, normal) = worst(p);
            assert!(cvd >= 6.0, "{name}: worst protan or deutan distance {cvd:.2} < 6");
            assert!(normal >= 15.0, "{name}: worst normal-vision distance {normal:.2} < 15");
        }
    }

    #[test]
    fn the_grey_of_the_far_regions_is_apart_from_every_hue_and_from_the_canvas() {
        for (mode, set) in [(Mode::Dark, &REGIONS_ON_DARK_CANVAS), (Mode::Light, &REGIONS_ON_LIGHT_CANVAS)] {
            let grey = REGION_GREY.get(mode);
            assert!(contrast(grey, CANVAS.get(mode)) >= 3.0, "{mode:?}: grey on the canvas is {:.2}:1", contrast(grey, CANVAS.get(mode)));
            for c in set {
                for m in [&PROTAN, &DEUTAN] {
                    assert!(delta_e(grey, *c, Some(m)) >= 6.0, "{mode:?} {c:?}: grey is too near for a colour-blind reader ({:.1})", delta_e(grey, *c, Some(m)));
                }
                assert!(delta_e(grey, *c, None) >= 12.0, "{mode:?} {c:?}: grey too near ({:.1})", delta_e(grey, *c, None));
            }
            assert!(oklch(grey).1 <= 0.03, "grey must stay neutral");
        }
    }

    #[test]
    fn region_colours_stay_in_the_lightness_band_and_keep_their_chroma() {
        for (name, p, lo, hi) in [("dark", &REGIONS_ON_DARK_CANVAS, 0.48, 0.67), ("light", &REGIONS_ON_LIGHT_CANVAS, 0.43, 0.77)] {
            for c in p {
                let (l, ch) = oklch(*c);
                assert!(l >= lo && l <= hi, "{name} {c:?}: L {l:.3} outside {lo}-{hi}");
                assert!(ch >= 0.10, "{name} {c:?}: chroma {ch:.3} < 0.10 reads as grey");
            }
        }
    }

    #[test]
    fn neutrals_keep_their_chroma_low() {
        // the plum cast lives in text, accent and tint; the neutrals stay (almost) neutral (0.026 in the rule; the measured values round to it, so 0.0265)
        for mode in [Mode::Light, Mode::Dark] {
            for (name, p) in [("page", PAGE), ("card", CARD), ("raised", RAISED), ("separator", SEPARATOR), ("canvas", CANVAS)] {
                let (_, ch) = oklch(p.get(mode));
                assert!(ch <= 0.0265, "{mode:?} {name}: chroma {ch:.3} > 0.026");
            }
        }
    }

    #[test]
    fn the_constants_equal_tokens_json() {
        let doc: serde_json::Value = serde_json::from_str(include_str!("../../../docs/ux/tokens.json")).expect("tokens.json is JSON");
        let roles = doc["roles"].as_object().expect("roles");
        assert_eq!(roles.len(), ROLES.len(), "a role was added or removed on one side only");
        let h = |r: Rgb| format!("#{:02X}{:02X}{:02X}", r[0], r[1], r[2]);
        for (name, p) in ROLES {
            let j = &roles[name];
            assert_eq!(j["light"].as_str().map(str::to_uppercase), Some(h(p.light)), "{name} light");
            assert_eq!(j["dark"].as_str().map(str::to_uppercase), Some(h(p.dark)), "{name} dark");
        }
        for (key, p) in [("regions_on_dark_canvas", &REGIONS_ON_DARK_CANVAS), ("regions_on_light_canvas", &REGIONS_ON_LIGHT_CANVAS)] {
            let list: Vec<String> = doc["data"][key].as_array().expect(key).iter().map(|v| v.as_str().unwrap().to_uppercase()).collect();
            let mine: Vec<String> = p.iter().map(|c| h(*c)).collect();
            assert_eq!(list, mine, "{key}");
        }
    }
}
