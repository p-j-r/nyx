// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Small widgets of the shell: cards, chips, buttons in the brand style, sparklines and modal sheets.

use crate::theme;
use eframe::egui::{self, Align2, Color32, CornerRadius, Margin, Response, RichText, Sense, Stroke, Ui, Vec2, Role, WidgetInfo};

pub fn card<R>(ui: &mut Ui, fill: Color32, stroke: Color32, pad: i8, f: impl FnOnce(&mut Ui) -> R) -> R {
    egui::Frame::new().fill(fill).stroke(Stroke::new(1.0, stroke)).corner_radius(CornerRadius::same(12)).inner_margin(Margin::same(pad)).show(ui, |ui| ui.vertical(f).inner).inner
}

/// A white card with the line colour.
pub fn plain_card<R>(ui: &mut Ui, f: impl FnOnce(&mut Ui) -> R) -> R {
    let (fill, line) = (theme::card(ui.ctx()), theme::line(ui.ctx()));
    card(ui, fill, line, 14, f)
}

/// The one filled accent button of a view, with the text colour that has contrast on it (`on_accent`).
pub fn primary(ui: &mut Ui, text: impl Into<String>) -> Response {
    let (c, fg) = (theme::accent(ui.ctx()), theme::on_accent(ui.ctx()));
    ui.add(egui::Button::new(RichText::new(text.into()).color(fg).strong()).fill(c).stroke(Stroke::NONE).min_size(Vec2::new(0.0, 32.0)))
}

pub fn primary_enabled(ui: &mut Ui, enabled: bool, text: impl Into<String>) -> Response {
    ui.add_enabled_ui(enabled, |ui| primary(ui, text)).inner
}

/// A white card; returns the response of the whole card (for clicks).
pub fn plain_card_r(ui: &mut Ui, f: impl FnOnce(&mut Ui)) -> Response {
    let (fill, line) = (theme::card(ui.ctx()), theme::line(ui.ctx()));
    egui::Frame::new().fill(fill).stroke(Stroke::new(1.0, line)).corner_radius(CornerRadius::same(12)).inner_margin(Margin::same(14)).show(ui, |ui| {
        ui.vertical(f);
    }).response
}

/// Outlined button.
pub fn ghost(ui: &mut Ui, text: impl Into<String>) -> Response {
    ui.add(egui::Button::new(text.into()).min_size(Vec2::new(0.0, 30.0)))
}

/// Outlined button in the danger colour (Cancel).
pub fn danger(ui: &mut Ui, text: impl Into<String>) -> Response {
    let c = theme::danger(ui.ctx());
    ui.add(egui::Button::new(RichText::new(text.into()).color(c)).stroke(Stroke::new(1.0, c)).min_size(Vec2::new(0.0, 30.0)))
}

/// A selectable pill: plum when on.
pub fn pill(ui: &mut Ui, on: bool, text: impl Into<String>) -> Response {
    let ctx = ui.ctx().clone();
    let (fill, fg, stroke) = if on { (theme::accent(&ctx), theme::on_accent(&ctx), theme::accent(&ctx)) } else { (theme::card(&ctx), theme::text(&ctx), theme::edge(&ctx)) };
    ui.add(egui::Button::new(RichText::new(text.into()).color(fg)).fill(fill).stroke(Stroke::new(1.0, stroke)).corner_radius(CornerRadius::same(8)))
}

/// A small rounded label such as a tag.
pub fn chip(ui: &mut Ui, text: impl Into<String>) -> Response {
    let ctx = ui.ctx().clone();
    let t = RichText::new(text.into()).size(12.5).color(theme::muted(&ctx));
    egui::Frame::new().fill(theme::raised(&ctx)).corner_radius(CornerRadius::same(10)).inner_margin(Margin::symmetric(9, 3)).show(ui, |ui| ui.label(t)).response
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

/// What a workflow, a run or the solver is doing. Colour is never the only channel: every state has its own glyph (drawn, so it needs no font)
/// and a word.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum State {
    /// saved, converged, passed: a check mark
    Done,
    /// running: a filled amber disc
    Running,
    /// queued, waiting: an amber ring
    Waiting,
    /// failed: a cross
    #[allow(dead_code)] // drawn by the sweep table
    Failed,
    /// nothing to report, or not saved yet: a hollow grey ring
    Idle,
    /// a result a convergence check passed on: a green ring with a check mark inside (productivity plan P6)
    Checked,
    /// a result nobody has checked: a ring, its left half filled, in the muted text colour
    Exploratory,
    /// a result whose check would not fit this machine: an amber ring with a slash
    NotAffordable,
}

impl State {
    pub fn word(self) -> &'static str {
        match self {
            State::Done => "Done",
            State::Running => "Running",
            State::Waiting => "Queued",
            State::Failed => "Failed",
            State::Idle => "Idle",
            State::Checked => "Checked",
            State::Exploratory => "Exploratory",
            State::NotAffordable => "Check not affordable",
        }
    }
}

/// The glyph of a state, 14 px square.
pub fn glyph(ui: &mut Ui, state: State) -> Response {
    let (rect, r) = ui.allocate_exact_size(Vec2::splat(14.0), Sense::hover());
    paint_glyph(ui, rect.center(), state);
    r.widget_info(|| WidgetInfo::labeled(Role::Label, true, state.word()));
    r
}

pub fn paint_glyph(ui: &Ui, c: egui::Pos2, state: State) {
    let ctx = ui.ctx();
    let p = ui.painter();
    match state {
        State::Done => {
            let col = theme::success(ctx);
            p.add(egui::Shape::line(vec![c + Vec2::new(-4.5, 0.5), c + Vec2::new(-1.5, 3.5), c + Vec2::new(4.5, -3.5)], Stroke::new(2.0, col)));
        }
        State::Running => {
            // an amber fill alone is 1.6:1 on a white card: the outline in the dark amber carries the shape
            p.circle(c, 5.0, theme::activity(ctx), Stroke::new(1.5, theme::activity_text(ctx)));
        }
        State::Waiting => {
            p.circle_stroke(c, 4.5, Stroke::new(2.0, theme::activity_text(ctx)));
        }
        State::Failed => {
            let col = theme::danger(ctx);
            for (a, b) in [((-4.0, -4.0), (4.0, 4.0)), ((-4.0, 4.0), (4.0, -4.0))] {
                p.line_segment([c + Vec2::new(a.0, a.1), c + Vec2::new(b.0, b.1)], Stroke::new(2.0, col));
            }
        }
        State::Idle => {
            p.circle_stroke(c, 4.5, Stroke::new(1.5, theme::edge(ctx)));
        }
        State::Checked => {
            let col = theme::success(ctx);
            p.circle_stroke(c, 6.0, Stroke::new(1.5, col));
            p.add(egui::Shape::line(vec![c + Vec2::new(-3.0, 0.5), c + Vec2::new(-1.0, 2.5), c + Vec2::new(3.0, -2.5)], Stroke::new(1.5, col)));
        }
        State::Exploratory => {
            // the ring whole and the left half of the disc filled: "partly there", told without a colour
            let col = theme::muted(ctx);
            p.circle_stroke(c, 5.5, Stroke::new(1.5, col));
            let half: Vec<egui::Pos2> = (0..=12)
                .map(|i| {
                    let a = std::f32::consts::FRAC_PI_2 + std::f32::consts::PI * i as f32 / 12.0;
                    c + Vec2::new(a.cos(), -a.sin()) * 5.5
                })
                .collect();
            p.add(egui::Shape::convex_polygon(half, col, Stroke::NONE));
        }
        State::NotAffordable => {
            let col = theme::activity_text(ctx);
            p.circle_stroke(c, 6.0, Stroke::new(1.5, col));
            p.line_segment([c + Vec2::new(-4.2, 4.2), c + Vec2::new(4.2, -4.2)], Stroke::new(1.5, col));
        }
    }
}

/// A state with its word (or `label` instead of the word): glyph, then text in the colour that reads on the page.
pub fn status(ui: &mut Ui, state: State, label: impl Into<String>) {
    let ctx = ui.ctx().clone();
    let color = match state {
        State::Done | State::Checked => theme::success(&ctx),
        State::Failed => theme::danger(&ctx),
        State::Running | State::Waiting | State::NotAffordable => theme::activity_text(&ctx),
        State::Idle | State::Exploratory => theme::muted(&ctx),
    };
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 6.0;
        glyph(ui, state);
        ui.label(RichText::new(label.into()).color(color));
    });
}

/// A row of text tabs (ADR 124): the one on show in the text colour and bold, with a 2 px accent line under it; the others muted. Returns the index clicked.
/// Quiet on purpose: the accent is for the selection and the one action, not for every control.
pub fn tabs(ui: &mut Ui, labels: &[&str], on: usize) -> Option<usize> {
    let ctx = ui.ctx().clone();
    let mut hit = None;
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 16.0;
        for (i, l) in labels.iter().enumerate() {
            let sel = i == on;
            let t = RichText::new(*l).size(14.5).color(if sel { theme::text(&ctx) } else { theme::muted(&ctx) });
            let r = ui.add(egui::Label::new(if sel { t.strong() } else { t }).sense(Sense::click()));
            if sel {
                let y = r.rect.bottom() + 4.0;
                ui.painter().line_segment([egui::pos2(r.rect.left(), y), egui::pos2(r.rect.right(), y)], Stroke::new(2.0, theme::accent(&ctx)));
            }
            if r.clicked() {
                hit = Some(i);
            }
            r.on_hover_cursor(egui::CursorIcon::PointingHand);
        }
    });
    hit
}

/// What `tab` reports: a click on the tab and a click on its close target.
pub struct TabHit {
    pub select: Response,
    pub close: Option<Response>,
}

/// One tab of the top bar: 28 px high, a state glyph, the title, and (when `closable`) a 24 x 24 close target. Selected = tint and a 2 px accent
/// bar under it; keyboard focus draws the 2 px accent ring; a screen reader gets the title as a selectable button.
/// The width of a tab of `text_w` points of title (the same sum `tab` allocates).
fn tab_width_of(text_w: f32, has_state: bool, closable: bool) -> f32 {
    let (pad, gap, close_w) = (10.0, 6.0, 24.0);
    pad + if has_state { 14.0 + gap } else { 0.0 } + text_w + if closable { gap + close_w + 2.0 } else { pad }
}

/// How wide `tab` draws this title, so a bar can decide which tabs fit before it draws any.
pub fn tab_width(ui: &Ui, title: &str, has_state: bool, closable: bool) -> f32 {
    let font = egui::TextStyle::Button.resolve(ui.style());
    tab_width_of(ui.painter().layout_no_wrap(title.to_string(), font, Color32::WHITE).size().x, has_state, closable)
}

pub fn tab(ui: &mut Ui, id: impl std::hash::Hash + std::fmt::Debug, title: &str, selected: bool, state: Option<State>, closable: bool) -> TabHit {
    let ctx = ui.ctx().clone();
    let id = ui.make_persistent_id(id);
    let font = egui::TextStyle::Button.resolve(ui.style());
    let col = if selected { theme::accent(&ctx) } else { theme::text(&ctx) };
    let galley = ui.painter().layout_no_wrap(title.to_string(), font, col);
    let (pad, gap, h) = (10.0, 6.0, 28.0);
    let glyph_w = if state.is_some() { 14.0 + gap } else { 0.0 };
    let close_w = if closable { 24.0 } else { 0.0 };
    let w = tab_width_of(galley.size().x, state.is_some(), closable);
    let (rect, _) = ui.allocate_exact_size(Vec2::new(w, h), Sense::hover());
    let close_rect = egui::Rect::from_min_size(egui::pos2(rect.right() - 2.0 - close_w, rect.center().y - 12.0), Vec2::splat(24.0));
    let select_rect = if closable { egui::Rect::from_min_max(rect.min, egui::pos2(close_rect.left(), rect.max.y)) } else { rect };
    let select = ui.interact(select_rect, id.with("select"), Sense::click());
    let close = closable.then(|| ui.interact(close_rect, id.with("close"), Sense::click()));
    select.widget_info(|| WidgetInfo::selected(Role::Button, true, selected, title));
    if let Some(c) = &close {
        c.widget_info(|| WidgetInfo::labeled(Role::Button, true, format!("Close {title}")));
    }
    let hovered = select.hovered() || close.as_ref().is_some_and(|c| c.hovered());
    let fill = if selected { theme::tint(&ctx) } else if hovered { theme::raised(&ctx) } else { Color32::TRANSPARENT };
    let p = ui.painter();
    p.rect_filled(rect, 8.0, fill);
    if selected {
        p.rect_filled(egui::Rect::from_min_max(egui::pos2(rect.left() + 6.0, rect.bottom() - 2.0), egui::pos2(rect.right() - 6.0, rect.bottom())), 1.0, theme::accent(&ctx));
    }
    let mut x = rect.left() + pad;
    if let Some(s) = state {
        paint_glyph(ui, egui::pos2(x + 7.0, rect.center().y), s);
        x += glyph_w;
    }
    ui.painter().galley(egui::pos2(x, rect.center().y - galley.size().y / 2.0), galley, col);
    if let Some(c) = &close {
        let cc = close_rect.center();
        let colr = if c.hovered() { theme::text(&ctx) } else { theme::muted(&ctx) };
        if c.hovered() {
            ui.painter().rect_filled(close_rect.shrink(2.0), 6.0, theme::raised(&ctx));
        }
        for (a, b) in [((-4.0, -4.0), (4.0, 4.0)), ((-4.0, 4.0), (4.0, -4.0))] {
            ui.painter().line_segment([cc + Vec2::new(a.0, a.1), cc + Vec2::new(b.0, b.1)], Stroke::new(1.5, colr));
        }
    }
    // the focus ring: egui draws a focused control with the active visuals; the tab is drawn by hand, so it draws the same ring
    if select.has_focus() || close.as_ref().is_some_and(|c| c.has_focus()) {
        ui.painter().rect_stroke(rect, 8.0, Stroke::new(2.0, theme::accent(&ctx)), egui::StrokeKind::Inside);
    }
    TabHit { select, close }
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
    // an outlined track: the empty part of a bar is part of the information
    ui.painter().rect_filled(rect, 4.0, theme::raised(ui.ctx()));
    ui.painter().rect_stroke(rect, 4.0, Stroke::new(1.0, theme::edge(ui.ctx())), egui::StrokeKind::Inside);
    let color = theme::accent(ui.ctx());
    match fraction {
        Some(f) => {
            let mut r = rect;
            r.set_width(rect.width() * f.clamp(0.0, 1.0));
            ui.painter().rect_filled(r, 4.0, color);
        }
        None => {
            if theme::reduce_motion(ui.ctx()) {
                // reduced motion: a still block in the middle, no animation and no repaint
                let w = rect.width() * 0.25;
                ui.painter().rect_filled(egui::Rect::from_min_size(egui::pos2(rect.center().x - w / 2.0, rect.top()), Vec2::new(w, rect.height())), 4.0, color);
                return;
            }
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
    egui::Area::new(egui::Id::unique((id, "backdrop"))).order(egui::Order::Foreground).fixed_pos(screen.min).show(ctx, |ui| {
        let r = ui.allocate_response(screen.size(), Sense::click());
        ui.painter().rect_filled(screen, 0.0, theme::scrim(110));
        if r.clicked() {
            close = true;
        }
    });
    let frame = egui::Frame::new().fill(theme::card(ctx)).stroke(Stroke::new(1.0, theme::line(ctx))).corner_radius(CornerRadius::same(16)).inner_margin(Margin::same(22)).shadow(egui::Shadow { offset: [0, 8], blur: 40, spread: 0, color: theme::scrim(70) });
    egui::Area::new(egui::Id::unique((id, "sheet"))).order(egui::Order::Foreground).anchor(Align2::CENTER_CENTER, [0.0, 0.0]).show(ctx, |ui| {
        // an anchored area is only offered the size it had last frame, so a sheet could never grow past its first (sizing) frame: offer the height of the window
        ui.set_max_height(screen.height() - 40.0);
        frame.show(ui, |ui| {
            ui.set_width(width.min(screen.width() - 60.0));
            // every sheet can be closed by name, not only by the Esc key or a click outside (a control, so it is also reachable by keyboard)
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                if ui.button("Close  Esc").clicked() {
                    close = true;
                }
            });
            egui::ScrollArea::vertical().max_height(screen.height() - 120.0).auto_shrink([false, true]).show(ui, f);
        });
    });
    if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
        close = true;
    }
    close
}
