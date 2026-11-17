//! The nyx look: every colour of the shell comes from `tokens.rs` (ADR 88) through the functions here, in a light and a dark look, both installed
//! once (`apply_mode`) so a switch of the theme only changes which one is active. Body text is Inter, numbers JetBrains Mono and titles EB Garamond
//! italic when they are installed (next to the executable in `fonts/`, or on the system); otherwise Segoe UI, Cascadia/Consolas and Georgia italic
//! are used, and finally egui's own fonts.

use crate::tokens::{self, Mode, Pair, Rgb};
use eframe::egui::{self, Color32, CornerRadius, FontData, FontDefinitions, FontFamily, FontId, Stroke, TextStyle, Theme, Visuals};
use std::sync::Arc;

/// The serif italic family of page titles.
pub fn title_family() -> FontFamily {
    FontFamily::Name("title".into())
}

pub fn title(size: f32) -> FontId {
    FontId::new(size, title_family())
}

/// A token colour as an egui colour.
pub fn rgb(c: Rgb) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

/// True when the colours are the dark ones (`apply_mode`).
pub fn is_dark(ctx: &egui::Context) -> bool {
    ctx.global_style().visuals.dark_mode
}

pub fn mode(ctx: &egui::Context) -> Mode {
    if is_dark(ctx) { Mode::Dark } else { Mode::Light }
}

/// The value of a token in the look that is active.
pub fn role(ctx: &egui::Context, p: Pair) -> Color32 {
    rgb(p.get(mode(ctx)))
}

pub fn page(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::PAGE)
}
pub fn card(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::CARD)
}
/// hover and pressed neutral fills, popups, chips
pub fn raised(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::RAISED)
}
/// Dividers and card outlines. Never the edge of a control: that is `edge`.
pub fn line(ctx: &egui::Context) -> Color32 {
    if look(ctx).high_contrast { edge(ctx) } else { role(ctx, tokens::SEPARATOR) }
}
/// The edge of buttons, fields and checkboxes (3:1 against the page and the card).
pub fn edge(ctx: &egui::Context) -> Color32 {
    role(ctx, if look(ctx).high_contrast { tokens::CONTROL_EDGE_HIGH } else { tokens::CONTROL_EDGE })
}
pub fn text(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::TEXT)
}
pub fn muted(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::MUTED)
}
pub fn accent(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ACCENT)
}
/// Text on an accent fill: white on the plum, ink on the dark pink.
pub fn on_accent(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ON_ACCENT)
}
pub fn tint(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::TINT)
}
/// Amber: time and attention (running, queued, stale) and nothing else.
pub fn activity(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ACTIVITY)
}
#[allow(dead_code)] // text on an amber fill: the busy chip of plan step U5
pub fn on_activity(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ON_ACTIVITY)
}
/// The activity colour as text (on `warm`).
pub fn activity_text(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ACTIVITY_TEXT)
}
/// Background of the busy chip and the stale band.
pub fn warm(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::ACTIVITY_TINT)
}
pub fn success(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::SUCCESS)
}
pub fn danger(ctx: &egui::Context) -> Color32 {
    role(ctx, tokens::DANGER)
}
/// The black scrim behind a sheet and under a shadow (a veil, not a palette colour).
pub fn scrim(alpha: u8) -> Color32 {
    Color32::from_black_alpha(alpha)
}

/// The settings of the look that are not the light or dark colours.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Look {
    pub dark: bool,
    /// control edges at 4.5:1 and 1.5 px, separators as strong as edges
    pub high_contrast: bool,
    /// no easing, no sliding progress block
    pub reduce_motion: bool,
}

fn look_id() -> egui::Id {
    egui::Id::unique("nyx-look")
}

/// The look that was installed last (`apply_look`).
pub fn look(ctx: &egui::Context) -> Look {
    ctx.data(|d| d.get_temp::<Look>(look_id())).unwrap_or_default()
}

pub fn reduce_motion(ctx: &egui::Context) -> bool {
    look(ctx).reduce_motion
}

/// The viewport and plot background: the dark canvas by default (also in a light window), or the canvas of the look when `follow_theme`.
pub fn canvas(ctx: &egui::Context, follow_theme: bool) -> Color32 {
    if follow_theme { role(ctx, tokens::CANVAS) } else { rgb(tokens::CANVAS.dark) }
}

/// True when the system asks for High Contrast: read from Windows (`SPI_GETHIGHCONTRAST`); other platforms have no such flag we read yet.
pub fn system_high_contrast() -> bool {
    platform_high_contrast()
}

#[cfg(not(windows))]
fn platform_high_contrast() -> bool {
    false
}

#[cfg(windows)]
fn platform_high_contrast() -> bool {
    {
        #[repr(C)]
        struct HighContrastW {
            size: u32,
            flags: u32,
            default_scheme: *mut u16,
        }
        #[link(name = "user32")]
        unsafe extern "system" {
            fn SystemParametersInfoW(action: u32, param: u32, data: *mut HighContrastW, win_ini: u32) -> i32;
        }
        const SPI_GETHIGHCONTRAST: u32 = 0x0042;
        const HCF_HIGHCONTRASTON: u32 = 0x0000_0001;
        let mut hc = HighContrastW { size: std::mem::size_of::<HighContrastW>() as u32, flags: 0, default_scheme: std::ptr::null_mut() };
        // SAFETY: `hc` is a valid HIGHCONTRASTW with its size set; the call only writes into it.
        let ok = unsafe { SystemParametersInfoW(SPI_GETHIGHCONTRAST, hc.size, &mut hc, 0) };
        ok != 0 && hc.flags & HCF_HIGHCONTRASTON != 0
    }
}

/// The High contrast setting ("auto", "standard", "high") as a yes or no.
pub fn resolve_contrast(setting: &str) -> bool {
    match setting {
        "high" => true,
        "standard" => false,
        _ => system_high_contrast(),
    }
}

/// Light visuals (what tests and the first frame use).
pub fn apply(ctx: &egui::Context) {
    apply_mode(ctx, false);
}

/// The full `Visuals` of one look. Neutral fills come from the tokens, the focus ring is the 2 px accent stroke of `widgets.active` (egui draws
/// a focused widget with the active visuals), every control has a 3:1 edge and dividers use the separator colour only.
fn build_visuals(m: Mode, high: bool) -> Visuals {
    let g = |p: Pair| rgb(p.get(m));
    let mut v = if m == Mode::Dark { Visuals::dark() } else { Visuals::light() };
    let edge_role = if high { tokens::CONTROL_EDGE_HIGH } else { tokens::CONTROL_EDGE };
    let (page, card, raised, line, edge, text, accent, tint) = (g(tokens::PAGE), g(tokens::CARD), g(tokens::RAISED), if high { g(edge_role) } else { g(tokens::SEPARATOR) }, g(edge_role), g(tokens::TEXT), g(tokens::ACCENT), g(tokens::TINT));
    v.panel_fill = page;
    v.window_fill = card;
    v.extreme_bg_color = card;
    v.faint_bg_color = page.lerp_to_gamma(raised, 0.5);
    v.override_text_color = Some(text);
    v.selection.bg_fill = card.lerp_to_gamma(accent, 0.30);
    v.selection.stroke = Stroke::new(1.0, accent);
    v.hyperlink_color = accent;
    v.warn_fg_color = g(tokens::ACTIVITY_TEXT);
    v.error_fg_color = g(tokens::DANGER);
    v.window_stroke = Stroke::new(1.0, line);
    v.window_corner_radius = CornerRadius::same(14);
    let r = CornerRadius::same(8);
    for w in [&mut v.widgets.noninteractive, &mut v.widgets.inactive, &mut v.widgets.hovered, &mut v.widgets.active, &mut v.widgets.open] {
        w.corner_radius = r;
        w.fg_stroke = Stroke::new(1.0, text);
    }
    v.widgets.noninteractive.bg_fill = page;
    v.widgets.noninteractive.weak_bg_fill = page;
    v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, line);
    v.widgets.inactive.bg_fill = card;
    v.widgets.inactive.weak_bg_fill = card;
    v.widgets.inactive.bg_stroke = Stroke::new(if high { 1.5 } else { 1.0 }, edge);
    v.widgets.hovered.bg_fill = raised;
    v.widgets.hovered.weak_bg_fill = raised;
    v.widgets.hovered.bg_stroke = Stroke::new(1.0, accent);
    v.widgets.active.bg_fill = tint;
    v.widgets.active.weak_bg_fill = tint;
    v.widgets.active.bg_stroke = Stroke::new(2.0, accent);
    v.widgets.open.bg_fill = raised;
    v.widgets.open.weak_bg_fill = raised;
    v.widgets.open.bg_stroke = Stroke::new(1.0, edge);
    v
}

/// Light or dark with the other options of the look as they were last set (`apply_look`).
pub fn apply_mode(ctx: &egui::Context, dark: bool) {
    apply_look(ctx, Look { dark, ..look(ctx) });
}

/// Installs the light and the dark look and makes `look.dark` the active one. Both are installed every time (it is cheap) so nothing depends on
/// which was installed first. The window decides light or dark itself (Settings: follow the system, light, dark) and pins the theme here, so egui
/// never switches on its own.
pub fn apply_look(ctx: &egui::Context, look: Look) {
    ctx.data_mut(|d| d.insert_temp(look_id(), look));
    ctx.set_visuals_of(Theme::Light, build_visuals(Mode::Light, look.high_contrast));
    ctx.set_visuals_of(Theme::Dark, build_visuals(Mode::Dark, look.high_contrast));
    let (dark, motion) = (look.dark, look.reduce_motion);
    // spacing and fonts belong to the style of each theme: set them on both, or a switch of theme would lose them
    ctx.all_styles_mut(|style| {
        style.spacing.item_spacing = egui::vec2(8.0, 6.0);
        style.spacing.button_padding = egui::vec2(12.0, 5.0);
        style.spacing.interact_size.y = 28.0;
        style.animation_time = if motion { 0.0 } else { 1.0 / 12.0 };
        style.text_styles.insert(TextStyle::Body, FontId::new(14.0, FontFamily::Proportional));
        style.text_styles.insert(TextStyle::Button, FontId::new(14.0, FontFamily::Proportional));
        style.text_styles.insert(TextStyle::Heading, FontId::new(18.0, FontFamily::Proportional));
        style.text_styles.insert(TextStyle::Small, FontId::new(12.0, FontFamily::Proportional));
        style.text_styles.insert(TextStyle::Monospace, FontId::new(13.0, FontFamily::Monospace));
    });
    ctx.set_theme(if dark { Theme::Dark } else { Theme::Light });
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
