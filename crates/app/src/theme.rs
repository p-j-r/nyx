//! Dark, compact look.

use eframe::egui::{self, Color32, CornerRadius, FontFamily, FontId, TextStyle};

pub fn apply(ctx: &egui::Context) {
    let mut v = egui::Visuals::dark();
    v.panel_fill = Color32::from_rgb(31, 34, 40);
    v.window_fill = Color32::from_rgb(36, 39, 46);
    v.extreme_bg_color = Color32::from_rgb(22, 24, 28);
    v.faint_bg_color = Color32::from_rgb(38, 41, 48);
    v.selection.bg_fill = Color32::from_rgb(52, 110, 170);
    v.hyperlink_color = Color32::from_rgb(110, 175, 240);
    let r = CornerRadius::same(4);
    v.widgets.noninteractive.corner_radius = r;
    v.widgets.inactive.corner_radius = r;
    v.widgets.hovered.corner_radius = r;
    v.widgets.active.corner_radius = r;
    v.widgets.inactive.bg_fill = Color32::from_rgb(52, 56, 66);
    v.widgets.inactive.weak_bg_fill = Color32::from_rgb(52, 56, 66);
    v.widgets.hovered.weak_bg_fill = Color32::from_rgb(66, 72, 86);
    ctx.set_visuals(v);
    let mut style = (*ctx.global_style()).clone();
    style.spacing.item_spacing = egui::vec2(7.0, 5.0);
    style.spacing.button_padding = egui::vec2(8.0, 3.0);
    style.text_styles.insert(TextStyle::Body, FontId::new(14.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Button, FontId::new(14.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Heading, FontId::new(18.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Monospace, FontId::new(13.0, FontFamily::Monospace));
    ctx.set_global_style(style);
}
