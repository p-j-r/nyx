//! Small widgets of the shell: cards, chips, buttons in the brand style, sparklines and modal sheets.

use crate::theme;
use eframe::egui::{self, Align2, Color32, CornerRadius, FontId, Margin, Response, RichText, Sense, Stroke, Ui, Vec2};

pub fn card<R>(ui: &mut Ui, fill: Color32, stroke: Color32, pad: i8, f: impl FnOnce(&mut Ui) -> R) -> R {
    egui::Frame::new().fill(fill).stroke(Stroke::new(1.0, stroke)).corner_radius(CornerRadius::same(12)).inner_margin(Margin::same(pad)).show(ui, f).inner
}

/// A white card with the line colour.
pub fn plain_card<R>(ui: &mut Ui, f: impl FnOnce(&mut Ui) -> R) -> R {
    let (fill, line) = (theme::card(ui.ctx()), theme::line(ui.ctx()));
    card(ui, fill, line, 14, f)
}

/// Plum filled button with white text.
pub fn primary(ui: &mut Ui, text: impl Into<String>) -> Response {
    let c = theme::accent(ui.ctx());
    ui.add(egui::Button::new(RichText::new(text.into()).color(Color32::WHITE).strong()).fill(c).stroke(Stroke::NONE).min_size(Vec2::new(0.0, 32.0)))
}

pub fn primary_enabled(ui: &mut Ui, enabled: bool, text: impl Into<String>) -> Response {
    ui.add_enabled_ui(enabled, |ui| primary(ui, text)).inner
}

/// A white card; returns the response of the whole card (for clicks).
pub fn plain_card_r(ui: &mut Ui, f: impl FnOnce(&mut Ui)) -> Response {
    let (fill, line) = (theme::card(ui.ctx()), theme::line(ui.ctx()));
    egui::Frame::new().fill(fill).stroke(Stroke::new(1.0, line)).corner_radius(CornerRadius::same(12)).inner_margin(Margin::same(14)).show(ui, f).response
}

/// Outlined button.
pub fn ghost(ui: &mut Ui, text: impl Into<String>) -> Response {
    ui.add(egui::Button::new(text.into()).min_size(Vec2::new(0.0, 30.0)))
}

/// Outlined button in the danger colour (Cancel).
pub fn danger(ui: &mut Ui, text: impl Into<String>) -> Response {
    ui.add(egui::Button::new(RichText::new(text.into()).color(theme::RED)).stroke(Stroke::new(1.0, theme::RED)).min_size(Vec2::new(0.0, 30.0)))
}

/// A selectable pill: plum when on.
pub fn pill(ui: &mut Ui, on: bool, text: impl Into<String>) -> Response {
    let ctx = ui.ctx().clone();
    let (fill, fg, stroke) = if on { (theme::accent(&ctx), Color32::WHITE, theme::accent(&ctx)) } else { (theme::card(&ctx), theme::text(&ctx), theme::line(&ctx)) };
    ui.add(egui::Button::new(RichText::new(text.into()).color(fg)).fill(fill).stroke(Stroke::new(1.0, stroke)).corner_radius(CornerRadius::same(8)))
}

/// A pill that is not available yet: dimmed, with "soon" in the hover text.
pub fn soon(ui: &mut Ui, text: impl Into<String>, why: &str) -> Response {
    ui.add_enabled(false, egui::Button::new(RichText::new(text.into()).color(theme::muted(ui.ctx())))).on_disabled_hover_text(format!("Coming soon: {why}"))
}

/// A small rounded label such as a tag.
pub fn chip(ui: &mut Ui, text: impl Into<String>) -> Response {
    let ctx = ui.ctx().clone();
    let t = RichText::new(text.into()).size(12.5).color(theme::muted(&ctx));
    egui::Frame::new().fill(if theme::is_dark(&ctx) { Color32::from_rgb(0x2A, 0x1C, 0x31) } else { Color32::from_rgb(0xEE, 0xEC, 0xEF) }).corner_radius(CornerRadius::same(10)).inner_margin(Margin::symmetric(9, 3)).show(ui, |ui| ui.label(t)).response
}

/// A chip with a colour (status chip in the top bar, badges on the open-file screen).
pub fn badge(ui: &mut Ui, fill: Color32, fg: Color32, text: impl Into<String>) {
    egui::Frame::new().fill(fill).corner_radius(CornerRadius::same(14)).inner_margin(Margin::symmetric(12, 5)).show(ui, |ui| ui.label(RichText::new(text.into()).color(fg)));
}

/// "RUNNING" style section label with an optional grey hint after it.
pub fn label_caps(ui: &mut Ui, text: &str, hint: &str) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(text.to_uppercase()).size(12.5).strong().extra_letter_spacing(1.2));
        if !hint.is_empty() {
            ui.label(RichText::new(hint).color(theme::muted(ui.ctx())));
        }
    });
}

pub fn dot(ui: &mut Ui, color: Color32) {
    let (r, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
    ui.painter().circle_filled(r.center(), 4.5, color);
}

pub fn mono(text: impl Into<String>) -> RichText {
    RichText::new(text.into()).monospace()
}

pub fn muted(ui: &Ui, text: impl Into<String>) -> RichText {
    RichText::new(text.into()).color(theme::muted(ui.ctx()))
}

pub fn title(text: impl Into<String>, size: f32) -> RichText {
    RichText::new(text.into()).font(theme::title(size))
}

/// A line of values as a polyline in a box (dashboard cards).
pub fn sparkline(ui: &mut Ui, values: &[f64], color: Color32, size: Vec2) {
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    if values.len() < 2 {
        ui.painter().line_segment([rect.left_center(), rect.right_center()], Stroke::new(1.0, theme::line(ui.ctx())));
        return;
    }
    let (lo, hi) = values.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| (l.min(*v), h.max(*v)));
    let span = (hi - lo).max(1e-300);
    let pts: Vec<egui::Pos2> = values
        .iter()
        .enumerate()
        .map(|(i, v)| egui::pos2(rect.left() + rect.width() * i as f32 / (values.len() - 1) as f32, rect.bottom() - 3.0 - (rect.height() - 6.0) * ((v - lo) / span) as f32))
        .collect();
    ui.painter().add(egui::Shape::line(pts, Stroke::new(2.0, color)));
}

/// A determinate progress bar in plum on a pale track; `None` shows a sliding block (no percentage reported).
pub fn bar(ui: &mut Ui, fraction: Option<f32>, width: f32) {
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, 8.0), Sense::hover());
    let track = if theme::is_dark(ui.ctx()) { Color32::from_rgb(0x3A, 0x2B, 0x42) } else { Color32::from_rgb(0xE9, 0xD7, 0xDF) };
    ui.painter().rect_filled(rect, 4.0, track);
    let color = theme::accent(ui.ctx());
    match fraction {
        Some(f) => {
            let mut r = rect;
            r.set_width(rect.width() * f.clamp(0.0, 1.0));
            ui.painter().rect_filled(r, 4.0, color);
        }
        None => {
            let t = ui.input(|i| i.time) as f32;
            let w = rect.width() * 0.25;
            let x = rect.left() + (rect.width() - w) * (0.5 + 0.5 * (t * 2.0).sin());
            ui.painter().rect_filled(egui::Rect::from_min_size(egui::pos2(x, rect.top()), Vec2::new(w, rect.height())), 4.0, color);
            ui.ctx().request_repaint();
        }
    }
}

/// A modal sheet: dims the page, swallows clicks outside and shows `f` in a centred card. Returns `true` when it should close
/// (Esc or a click on the dimmed page).
pub fn modal(ctx: &egui::Context, id: &str, width: f32, f: impl FnOnce(&mut Ui)) -> bool {
    let screen = ctx.content_rect();
    let mut close = false;
    egui::Area::new(egui::Id::new((id, "backdrop"))).order(egui::Order::Foreground).fixed_pos(screen.min).show(ctx, |ui| {
        let r = ui.allocate_response(screen.size(), Sense::click());
        ui.painter().rect_filled(screen, 0.0, Color32::from_black_alpha(110));
        if r.clicked() {
            close = true;
        }
    });
    let frame = egui::Frame::new().fill(theme::card(ctx)).stroke(Stroke::new(1.0, theme::line(ctx))).corner_radius(CornerRadius::same(16)).inner_margin(Margin::same(22)).shadow(egui::Shadow { offset: [0, 8], blur: 40, spread: 0, color: Color32::from_black_alpha(70) });
    egui::Area::new(egui::Id::new((id, "sheet"))).order(egui::Order::Foreground).anchor(Align2::CENTER_CENTER, [0.0, 0.0]).show(ctx, |ui| {
        frame.show(ui, |ui| {
            ui.set_width(width.min(screen.width() - 60.0));
            egui::ScrollArea::vertical().max_height(screen.height() - 120.0).auto_shrink([false, true]).show(ui, f);
        });
    });
    if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
        close = true;
    }
    close
}

pub fn font(size: f32) -> FontId {
    FontId::proportional(size)
}
