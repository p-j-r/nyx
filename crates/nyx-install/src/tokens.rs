// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ADR 88 colour tokens for the optional `eframe`/`glow` installer wizard (`feature = "gui"`),
//! kept in exact sync with `docs/ux/tokens.json`.

use eframe::egui::Color32;

pub type Rgb = [u8; 3];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Light,
    Dark,
}

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

    pub fn color32(&self, dark: bool) -> Color32 {
        let [r, g, b] = if dark { self.dark } else { self.light };
        Color32::from_rgb(r, g, b)
    }
}

const fn digit(b: u8) -> u8 {
    match b {
        b'0'..=b'9' => b - b'0',
        b'a'..=b'f' => b - b'a' + 10,
        b'A'..=b'F' => b - b'A' + 10,
        _ => panic!("invalid hex digit"),
    }
}

pub const fn hex(s: &str) -> Rgb {
    let b = s.as_bytes();
    assert!(b.len() == 7 && b[0] == b'#');
    [
        digit(b[1]) * 16 + digit(b[2]),
        digit(b[3]) * 16 + digit(b[4]),
        digit(b[5]) * 16 + digit(b[6]),
    ]
}

const fn pair(light: &str, dark: &str) -> Pair {
    Pair {
        light: hex(light),
        dark: hex(dark),
    }
}

pub const PAGE: Pair = pair("#F8F6F8", "#141218");
pub const CARD: Pair = pair("#FFFFFF", "#1D1A22");
pub const RAISED: Pair = pair("#F1EEF2", "#262230");
pub const SEPARATOR: Pair = pair("#E6E4E6", "#2E2A35");
pub const CONTROL_EDGE: Pair = pair("#8A8091", "#6E6578");
pub const TEXT: Pair = pair("#1D1224", "#ECE7EF");
pub const MUTED: Pair = pair("#5C5062", "#A9A0B1");
pub const ACCENT: Pair = pair("#820C4F", "#D96AA4");
pub const ON_ACCENT: Pair = pair("#FFFFFF", "#1D1224");
pub const TINT: Pair = pair("#F1E4EC", "#33202F");
pub const ACTIVITY_TEXT: Pair = pair("#6B4300", "#F2B35E");
pub const ACTIVITY_TINT: Pair = pair("#FFF1DC", "#3D2D1A");
pub const SUCCESS: Pair = pair("#1F6B43", "#5FC08C");
pub const DANGER: Pair = pair("#A32020", "#F27A7A");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn installer_tokens_match_docs_ux_tokens_json() {
        let doc: serde_json::Value =
            serde_json::from_str(include_str!("../../../docs/ux/tokens.json"))
                .expect("tokens.json");
        let roles = doc["roles"].as_object().expect("roles");
        let h = |r: Rgb| format!("#{:02X}{:02X}{:02X}", r[0], r[1], r[2]);
        for (name, p) in [
            ("page", PAGE),
            ("card", CARD),
            ("raised", RAISED),
            ("separator", SEPARATOR),
            ("control_edge", CONTROL_EDGE),
            ("text", TEXT),
            ("muted", MUTED),
            ("accent", ACCENT),
            ("on_accent", ON_ACCENT),
            ("tint", TINT),
            ("activity_text", ACTIVITY_TEXT),
            ("activity_tint", ACTIVITY_TINT),
            ("success", SUCCESS),
            ("danger", DANGER),
        ] {
            let j = &roles[name];
            assert_eq!(
                j["light"].as_str().map(str::to_uppercase),
                Some(h(p.light)),
                "{name} light"
            );
            assert_eq!(
                j["dark"].as_str().map(str::to_uppercase),
                Some(h(p.dark)),
                "{name} dark"
            );
        }
    }
}
