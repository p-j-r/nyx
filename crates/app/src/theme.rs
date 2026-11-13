//! The nyx look: the brand palette, fonts and the egui visuals. Light by default, dark on request; every colour of the shell comes from here.
//!
//! Brand: amber #FFBD59, plum #820C4F, shadow #531F55, ink #24102A, surface #F7F6F7, line #E6E4E6. Body text is Inter, numbers JetBrains Mono and
//! titles EB Garamond italic when they are installed (next to the executable in `fonts/`, or on the system); otherwise Segoe UI, Cascadia/Consolas and
//! Georgia italic are used, and finally egui's own fonts.

use eframe::egui::{self, Color32, CornerRadius, FontData, FontDefinitions, FontFamily, FontId, TextStyle};
use std::sync::Arc;

pub const AMBER: Color32 = Color32::from_rgb(0xFF, 0xBD, 0x59);
pub const PLUM: Color32 = Color32::from_rgb(0x82, 0x0C, 0x4F);
pub const INK: Color32 = Color32::from_rgb(0x24, 0x10, 0x2A);
pub const SURFACE: Color32 = Color32::from_rgb(0xF7, 0xF6, 0xF7);
pub const LINE: Color32 = Color32::from_rgb(0xE6, 0xE4, 0xE6);
/// derived tints
pub const PLUM_TINT: Color32 = Color32::from_rgb(0xF2, 0xE6, 0xEE);
pub const AMBER_TINT: Color32 = Color32::from_rgb(0xFF, 0xF1, 0xDC);
pub const MUTED: Color32 = Color32::from_rgb(0x6E, 0x62, 0x73);
pub const GREEN: Color32 = Color32::from_rgb(0x2E, 0x7D, 0x52);
pub const RED: Color32 = Color32::from_rgb(0xB4, 0x2B, 0x2B);

/// The serif italic family of page titles.
pub fn title_family() -> FontFamily {
    FontFamily::Name("title".into())
}

pub fn title(size: f32) -> FontId {
    FontId::new(size, title_family())
}

/// True when the colours are the dark ones (`apply_mode`).
pub fn is_dark(ctx: &egui::Context) -> bool {
    ctx.global_style().visuals.dark_mode
}

/// Colours that depend on the mode: (page, card, line, text, muted text).
pub fn page(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0x17, 0x0E, 0x1B) } else { SURFACE }
}
pub fn card(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0x21, 0x15, 0x27) } else { Color32::WHITE }
}
pub fn line(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0x3A, 0x2B, 0x42) } else { LINE }
}
pub fn text(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0xF1, 0xEA, 0xF3) } else { INK }
}
pub fn muted(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0xA8, 0x9B, 0xAD) } else { MUTED }
}
pub fn accent(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0xD9, 0x6A, 0xA4) } else { PLUM }
}
pub fn tint(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0x3A, 0x1B, 0x3B) } else { PLUM_TINT }
}
pub fn warm(ctx: &egui::Context) -> Color32 {
    if is_dark(ctx) { Color32::from_rgb(0x3D, 0x2D, 0x1A) } else { AMBER_TINT }
}

/// Light visuals with the brand palette (what tests and the first frame use).
pub fn apply(ctx: &egui::Context) {
    apply_mode(ctx, false);
}

pub fn apply_mode(ctx: &egui::Context, dark: bool) {
    let mut v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    let (page, card, line, text, accent) = if dark {
        (Color32::from_rgb(0x17, 0x0E, 0x1B), Color32::from_rgb(0x21, 0x15, 0x27), Color32::from_rgb(0x3A, 0x2B, 0x42), Color32::from_rgb(0xF1, 0xEA, 0xF3), Color32::from_rgb(0xD9, 0x6A, 0xA4))
    } else {
        (SURFACE, Color32::WHITE, LINE, INK, PLUM)
    };
    v.panel_fill = page;
    v.window_fill = card;
    v.extreme_bg_color = if dark { Color32::from_rgb(0x10, 0x09, 0x13) } else { Color32::WHITE };
    v.faint_bg_color = if dark { Color32::from_rgb(0x1C, 0x12, 0x21) } else { Color32::from_rgb(0xFB, 0xFA, 0xFB) };
    v.override_text_color = Some(text);
    v.selection.bg_fill = if dark { Color32::from_rgb(0x6B, 0x25, 0x57) } else { Color32::from_rgb(0xE9, 0xCF, 0xDF) };
    v.selection.stroke = egui::Stroke::new(1.0, accent);
    v.hyperlink_color = accent;
    v.window_stroke = egui::Stroke::new(1.0, line);
    v.window_corner_radius = CornerRadius::same(14);
    let r = CornerRadius::same(8);
    for w in [&mut v.widgets.noninteractive, &mut v.widgets.inactive, &mut v.widgets.hovered, &mut v.widgets.active, &mut v.widgets.open] {
        w.corner_radius = r;
    }
    v.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, line);
    v.widgets.inactive.bg_fill = card;
    v.widgets.inactive.weak_bg_fill = card;
    v.widgets.inactive.bg_stroke = egui::Stroke::new(1.0, line);
    v.widgets.hovered.weak_bg_fill = if dark { Color32::from_rgb(0x30, 0x20, 0x38) } else { PLUM_TINT };
    v.widgets.hovered.bg_stroke = egui::Stroke::new(1.0, accent);
    v.widgets.active.weak_bg_fill = if dark { Color32::from_rgb(0x42, 0x28, 0x4A) } else { Color32::from_rgb(0xE9, 0xCF, 0xDF) };
    ctx.set_visuals(v);
    let mut style = (*ctx.global_style()).clone();
    style.spacing.item_spacing = egui::vec2(8.0, 6.0);
    style.spacing.button_padding = egui::vec2(12.0, 5.0);
    style.spacing.interact_size.y = 26.0;
    style.text_styles.insert(TextStyle::Body, FontId::new(14.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Button, FontId::new(14.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Heading, FontId::new(17.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Small, FontId::new(12.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Monospace, FontId::new(13.0, FontFamily::Monospace));
    ctx.set_global_style(style);
}

fn first_font(paths: &[&str]) -> Option<Vec<u8>> {
    let mut dirs: Vec<std::path::PathBuf> = vec![];
    if let Ok(exe) = std::env::current_exe() {
        if let Some(d) = exe.parent() {
            dirs.push(d.join("fonts"));
        }
    }
    dirs.push("crates/app/assets/fonts".into());
    dirs.push("C:/Windows/Fonts".into());
    if let Ok(h) = std::env::var("USERPROFILE") {
        dirs.push(std::path::PathBuf::from(h).join("AppData/Local/Microsoft/Windows/Fonts"));
    }
    dirs.push("/usr/share/fonts/truetype/inter".into());
    dirs.push("/usr/share/fonts/truetype/dejavu".into());
    for p in paths {
        for d in &dirs {
            if let Ok(b) = std::fs::read(d.join(p)) {
                return Some(b);
            }
        }
    }
    None
}

/// Installs the fonts that are available (see the module comment). Anything missing falls back to egui's own fonts.
pub fn install_fonts(ctx: &egui::Context) {
    let mut fonts = FontDefinitions::default();
    let mut add = |name: &str, bytes: Vec<u8>, families: &[FontFamily], first: bool| {
        fonts.font_data.insert(name.to_string(), Arc::new(FontData::from_owned(bytes)));
        for f in families {
            let list = fonts.families.entry(f.clone()).or_default();
            if first {
                list.insert(0, name.to_string());
            } else {
                list.push(name.to_string());
            }
        }
    };
    if let Some(b) = first_font(&["Inter-Regular.ttf", "Inter.ttf", "segoeui.ttf"]) {
        add("ui", b, &[FontFamily::Proportional], true);
    }
    if let Some(b) = first_font(&["JetBrainsMono-Regular.ttf", "JetBrainsMono.ttf", "CascadiaMono.ttf", "consola.ttf"]) {
        add("mono", b, &[FontFamily::Monospace], true);
    }
    match first_font(&["EBGaramond-Italic.ttf", "EBGaramond-Italic[wght].ttf", "georgiai.ttf"]) {
        Some(b) => add("title", b, &[title_family()], true),
        None => {
            fonts.families.insert(title_family(), vec![]);
        }
    }
    // a family must also fall back to the default proportional fonts for symbols the file lacks
    let prop = fonts.families.get(&FontFamily::Proportional).cloned().unwrap_or_default();
    let t = fonts.families.entry(title_family()).or_default();
    for n in prop {
        if !t.contains(&n) {
            t.push(n);
        }
    }
    ctx.set_fonts(fonts);
}
