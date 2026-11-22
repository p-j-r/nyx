// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Minimalist 640x480 ADR 88 GUI wizard over `eframe 0.37` (`glow` renderer, ~5-6 MB standalone binary).

use crate::bundle::{exe_filename, IntegrityVerifier, PayloadSource};
use crate::cli::{Action, ParsedArgs};
use crate::connect_hook::{self, ConnectPolicy};
use crate::manifest::InstallScope;
use crate::notices::CORRESPONDING_SOURCE_URL;
use crate::platform::{native_backend, HostPaths};
use crate::process::OsProcessInspector;
use crate::report::InstallReport;
use crate::tokens;
use crate::tx::{
    locate_manifest, run_install_or_upgrade, run_repair, run_uninstall, TxOptions, TxProgress,
};
use eframe::egui::{self, CornerRadius, Margin, RichText, Stroke, Vec2, Visuals};
use nyx_connect::{plain_path, Status};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::mpsc::{self, Receiver};
use std::sync::Arc;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Step {
    Welcome,
    Options,
    Progress,
    Done,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum WizardMode {
    InstallOrUpgrade,
    UninstallOrRepair,
}

enum WorkerMsg {
    Progress(TxProgress),
    Log(String),
    Finished(Box<Result<InstallReport, String>>),
}

struct WizardApp {
    mode: WizardMode,
    step: Step,
    is_upgrade: bool,
    existing_version: Option<String>,
    source_version: String,
    scope: InstallScope,
    custom_prefix: String,
    add_path: bool,
    desktop: bool,
    desktop_shortcut: bool,
    file_assoc: bool,
    uninstall_repair: bool,
    purge_cache: bool,
    host_selected: BTreeMap<String, bool>,
    host_rows: Vec<(String, String, Status, PathBuf)>,
    mcp_preview_exe: PathBuf,
    source: Option<Arc<PayloadSource>>,
    progress: f32,
    progress_label: String,
    logs: Vec<String>,
    show_logs: bool,
    rx: Option<Receiver<WorkerMsg>>,
    report: Option<InstallReport>,
    error: Option<String>,
    copied_snippet: Option<String>,
    allow_root_user: bool,
}

fn load_icon() -> Option<egui::IconData> {
    let img = image::load_from_memory(include_bytes!("../../../docs/brand/nyx-mark.png"))
        .ok()?
        .to_rgba8();
    let img = image::imageops::resize(&img, 128, 128, image::imageops::FilterType::Lanczos3);
    Some(egui::IconData {
        rgba: img.into_raw(),
        width: 128,
        height: 128,
    })
}

pub fn run_wizard(
    action: Action,
    parsed: ParsedArgs,
    source: Option<PayloadSource>,
) -> anyhow::Result<()> {
    let backend = native_backend();
    let paths = HostPaths::detect();
    let current_exe = std::env::current_exe().ok();
    let existing = locate_manifest(
        parsed.tx.prefix.as_deref(),
        current_exe.as_deref(),
        backend.as_ref(),
        parsed.tx.effective_scope(),
        &paths,
    );

    let mode = match action {
        Action::Uninstall | Action::Repair => WizardMode::UninstallOrRepair,
        _ => WizardMode::InstallOrUpgrade,
    };

    let source_version = source
        .as_ref()
        .and_then(|s| s.metadata_summary(&IntegrityVerifier).ok())
        .map(|o| o.toc.nyx_version)
        .unwrap_or_else(|| env!("CARGO_PKG_VERSION").to_string());

    let scope = parsed.tx.scope;
    let default_root = parsed
        .tx
        .prefix
        .clone()
        .unwrap_or_else(|| backend.default_install_root(scope, &paths));
    let layout = backend.resolve_layout(&default_root, scope, &paths);
    let mcp_preview_exe = layout.bin_dir.join(exe_filename("nyx-mcp", backend.os()));

    let mut host_selected = BTreeMap::new();
    let mut host_rows = Vec::new();
    if let Some(env) = paths.connect_env(backend.os()) {
        for h in nyx_connect::list(&env, &mcp_preview_exe) {
            let default_on = matches!(h.status, Status::NotConnected | Status::Connected)
                || matches!(&h.status, Status::Stale { found } if !nyx_connect::is_cargo_build_tree(Path::new(found)));
            host_selected.insert(h.id.to_string(), default_on);
            host_rows.push((h.id.to_string(), h.name.to_string(), h.status, h.config));
        }
    }

    let app = WizardApp {
        mode,
        step: Step::Welcome,
        is_upgrade: existing.is_some(),
        existing_version: existing.as_ref().map(|(_, m)| m.nyx_version.clone()),
        source_version,
        scope,
        custom_prefix: default_root.display().to_string(),
        add_path: parsed.tx.add_path,
        desktop: parsed.tx.desktop,
        desktop_shortcut: parsed.tx.desktop_shortcut,
        file_assoc: parsed.tx.file_assoc,
        uninstall_repair: action == Action::Repair,
        purge_cache: parsed.tx.purge_cache,
        host_selected,
        host_rows,
        mcp_preview_exe,
        source: source.map(Arc::new),
        progress: 0.0,
        progress_label: String::new(),
        logs: Vec::new(),
        show_logs: false,
        rx: None,
        report: None,
        error: None,
        copied_snippet: None,
        allow_root_user: parsed.tx.allow_root_user,
    };

    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size([640.0, 480.0])
        .with_min_inner_size([560.0, 420.0])
        .with_title("nyx Setup");
    if let Some(icon) = load_icon() {
        viewport = viewport.with_icon(icon);
    }

    let native_opts = eframe::NativeOptions {
        viewport,
        ..Default::default()
    };
    eframe::run_native("nyx Setup", native_opts, Box::new(|_cc| Ok(Box::new(app))))
        .map_err(|e| anyhow::anyhow!("GUI startup failed: {e}"))
}

impl WizardApp {
    fn sync_default_prefix_for_scope(&mut self) {
        let backend = native_backend();
        let paths = HostPaths::detect();
        self.custom_prefix = backend
            .default_install_root(self.scope, &paths)
            .display()
            .to_string();
    }

    fn start_worker(&mut self, ctx: &egui::Context) {
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        self.step = Step::Progress;
        self.progress = 0.02;
        self.progress_label = "Starting...".into();
        self.logs.clear();
        self.error = None;

        let mode = self.mode;
        let uninstall_repair = self.uninstall_repair;
        let is_upgrade = self.is_upgrade;
        let selected_hosts: Vec<String> = self
            .host_selected
            .iter()
            .filter(|&(_, &on)| on)
            .map(|(k, _)| k.clone())
            .collect();
        let connect_policy = if selected_hosts.is_empty() {
            ConnectPolicy::None
        } else {
            ConnectPolicy::Explicit(selected_hosts)
        };
        let mut opts = TxOptions::default_user();
        opts.scope = self.scope;
        opts.prefix = Some(PathBuf::from(self.custom_prefix.trim()));
        opts.add_path = self.add_path;
        opts.desktop = self.desktop;
        opts.desktop_shortcut = self.desktop_shortcut;
        opts.file_assoc = self.file_assoc;
        opts.purge_cache = self.purge_cache;
        opts.allow_root_user = self.allow_root_user;
        opts.connect_policy = connect_policy;
        let source = self.source.clone();
        let ctx_clone = ctx.clone();

        std::thread::spawn(move || {
            let backend = native_backend();
            let paths = HostPaths::detect();
            let inspector = OsProcessInspector;
            let mut on_prog = |p: TxProgress| {
                let _ = tx.send(WorkerMsg::Log(format!(
                    "[{}/{}] {}",
                    p.phase, p.total_phases, p.message
                )));
                let _ = tx.send(WorkerMsg::Progress(p));
                ctx_clone.request_repaint();
            };

            let res = match mode {
                WizardMode::UninstallOrRepair if !uninstall_repair => {
                    run_uninstall(&opts, backend.as_ref(), &inspector, &paths, &mut on_prog)
                }
                WizardMode::UninstallOrRepair => match source.as_deref() {
                    Some(src) => run_repair(
                        src,
                        &opts,
                        backend.as_ref(),
                        &inspector,
                        &paths,
                        &mut on_prog,
                    ),
                    None => Err(anyhow::anyhow!(
                        "no installer payload available to repair from"
                    )),
                },
                WizardMode::InstallOrUpgrade => match source.as_deref() {
                    Some(src) => {
                        let act = if is_upgrade { "upgrade" } else { "install" };
                        run_install_or_upgrade(
                            act,
                            src,
                            &opts,
                            backend.as_ref(),
                            &inspector,
                            &paths,
                            &mut on_prog,
                        )
                    }
                    None => Err(anyhow::anyhow!("no installer payload found")),
                },
            };

            let _ = tx.send(WorkerMsg::Finished(Box::new(
                res.map_err(|e| e.to_string()),
            )));
            ctx_clone.request_repaint();
        });
    }

    fn apply_theme(ctx: &egui::Context) {
        let dark = ctx.global_style().visuals.dark_mode;
        let mut vis = if dark {
            Visuals::dark()
        } else {
            Visuals::light()
        };
        vis.panel_fill = tokens::PAGE.color32(dark);
        vis.window_fill = tokens::CARD.color32(dark);
        vis.extreme_bg_color = tokens::RAISED.color32(dark);
        vis.widgets.inactive.bg_fill = tokens::CARD.color32(dark);
        vis.widgets.inactive.bg_stroke = Stroke::new(1.0, tokens::CONTROL_EDGE.color32(dark));
        vis.widgets.hovered.bg_fill = tokens::RAISED.color32(dark);
        vis.selection.bg_fill = tokens::TINT.color32(dark);
        vis.selection.stroke = Stroke::new(1.5, tokens::ACCENT.color32(dark));
        ctx.set_visuals(vis);
    }

    fn poll_worker(&mut self) {
        let Some(ref rx) = self.rx else { return };
        let mut done = false;
        while let Ok(msg) = rx.try_recv() {
            match msg {
                WorkerMsg::Progress(p) => {
                    self.progress = p.fraction;
                    self.progress_label = p.message;
                }
                WorkerMsg::Log(l) => self.logs.push(l),
                WorkerMsg::Finished(res) => {
                    match *res {
                        Ok(report) => self.report = Some(report),
                        Err(err) => self.error = Some(err),
                    }
                    self.step = Step::Done;
                    done = true;
                }
            }
        }
        if done {
            self.rx = None;
        }
    }
}

impl eframe::App for WizardApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        Self::apply_theme(&ctx);
        self.poll_worker();
        let dark = ctx.global_style().visuals.dark_mode;

        let card_fill = tokens::CARD.color32(dark);
        let sep = tokens::SEPARATOR.color32(dark);
        let accent = tokens::ACCENT.color32(dark);
        let on_accent = tokens::ON_ACCENT.color32(dark);
        let text_col = tokens::TEXT.color32(dark);
        let muted = tokens::MUTED.color32(dark);
        let success = tokens::SUCCESS.color32(dark);
        let danger = tokens::DANGER.color32(dark);

        egui::Panel::top("setup_header")
            .frame(
                egui::Frame::new()
                    .fill(card_fill)
                    .inner_margin(Margin::symmetric(20, 14))
                    .stroke(Stroke::new(1.0, sep)),
            )
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("nyx").size(26.0).strong().color(accent));
                    ui.add_space(6.0);
                    let subtitle = match self.mode {
                        WizardMode::UninstallOrRepair => format!(
                            "Maintenance · v{}",
                            self.existing_version
                                .as_deref()
                                .unwrap_or(&self.source_version)
                        ),
                        WizardMode::InstallOrUpgrade if self.is_upgrade => format!(
                            "Upgrade {} → {}",
                            self.existing_version.as_deref().unwrap_or("existing"),
                            self.source_version
                        ),
                        WizardMode::InstallOrUpgrade => format!("Setup · v{}", self.source_version),
                    };
                    ui.label(RichText::new(subtitle).size(14.0).color(muted));
                });
            });

        egui::Panel::bottom("setup_footer")
            .frame(
                egui::Frame::new()
                    .fill(card_fill)
                    .inner_margin(Margin::symmetric(20, 12))
                    .stroke(Stroke::new(1.0, sep)),
            )
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.with_layout(
                        egui::Layout::right_to_left(egui::Align::Center),
                        |ui| match self.step {
                            Step::Welcome => {
                                let next_label = if self.mode == WizardMode::UninstallOrRepair {
                                    "Continue"
                                } else if self.is_upgrade {
                                    "Configure upgrade"
                                } else {
                                    "Configure install"
                                };
                                if ui
                                    .add(
                                        egui::Button::new(
                                            RichText::new(next_label).color(on_accent).strong(),
                                        )
                                        .fill(accent)
                                        .min_size(Vec2::new(120.0, 32.0)),
                                    )
                                    .clicked()
                                {
                                    self.step = Step::Options;
                                }
                                if ui
                                    .add(
                                        egui::Button::new("Cancel").min_size(Vec2::new(80.0, 32.0)),
                                    )
                                    .clicked()
                                {
                                    ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                                }
                            }
                            Step::Options => {
                                let act_label = match self.mode {
                                    WizardMode::UninstallOrRepair if self.uninstall_repair => {
                                        "Repair"
                                    }
                                    WizardMode::UninstallOrRepair => "Uninstall",
                                    WizardMode::InstallOrUpgrade if self.is_upgrade => "Upgrade",
                                    WizardMode::InstallOrUpgrade => "Install",
                                };
                                let btn_fill = if self.mode == WizardMode::UninstallOrRepair
                                    && !self.uninstall_repair
                                {
                                    danger
                                } else {
                                    accent
                                };
                                if ui
                                    .add(
                                        egui::Button::new(
                                            RichText::new(act_label).color(on_accent).strong(),
                                        )
                                        .fill(btn_fill)
                                        .min_size(Vec2::new(110.0, 32.0)),
                                    )
                                    .clicked()
                                {
                                    self.start_worker(&ctx);
                                }
                                if ui
                                    .add(egui::Button::new("Back").min_size(Vec2::new(80.0, 32.0)))
                                    .clicked()
                                {
                                    self.step = Step::Welcome;
                                }
                            }
                            Step::Progress => {
                                ui.label(RichText::new("Working... please wait").color(muted));
                            }
                            Step::Done => {
                                if ui
                                    .add(
                                        egui::Button::new(
                                            RichText::new("Close").color(on_accent).strong(),
                                        )
                                        .fill(accent)
                                        .min_size(Vec2::new(90.0, 32.0)),
                                    )
                                    .clicked()
                                {
                                    ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                                }
                                if self.error.is_none() && self.mode == WizardMode::InstallOrUpgrade
                                {
                                    if let Some(ref rep) = self.report {
                                        let nyx_bin = native_backend()
                                            .resolve_layout(
                                                Path::new(&rep.install_root),
                                                rep.scope,
                                                &HostPaths::detect(),
                                            )
                                            .bin_dir
                                            .join(exe_filename("nyx", native_backend().os()));
                                        if nyx_bin.is_file()
                                            && ui
                                                .add(
                                                    egui::Button::new("Launch nyx")
                                                        .min_size(Vec2::new(100.0, 32.0)),
                                                )
                                                .clicked()
                                        {
                                            let _ = Command::new(&nyx_bin).spawn();
                                            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                                        }
                                    }
                                }
                            }
                        },
                    );
                });
            });

        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(tokens::PAGE.color32(dark)).inner_margin(Margin::same(20)))
            .show(ui, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| match self.step {
                    Step::Welcome => {
                        ui.label(
                            RichText::new(match self.mode {
                                WizardMode::UninstallOrRepair => "Uninstall or Repair nyx",
                                WizardMode::InstallOrUpgrade if self.is_upgrade => "Upgrade existing nyx installation",
                                WizardMode::InstallOrUpgrade => "Finite-element photonics, phononics & electromagnetics",
                            })
                            .size(18.0)
                            .strong()
                            .color(text_col),
                        );
                        ui.add_space(10.0);
                        egui::Frame::new()
                            .fill(card_fill)
                            .stroke(Stroke::new(1.0, sep))
                            .corner_radius(CornerRadius::same(10))
                            .inner_margin(Margin::same(14))
                            .show(ui, |ui| {
                                ui.set_width(ui.available_width());
                                ui.label(RichText::new("LICENSING & SOURCE").size(12.0).strong().color(muted));
                                ui.add_space(4.0);
                                ui.label(
                                    "nyx is dual licensed under GNU Affero General Public License v3.0 \
                                     (AGPL-3.0-only) or a Commercial License (LicenseRef-nyx-Commercial).",
                                );
                                ui.label(format!(
                                    "No warranty is provided (AGPL §§15-16). Corresponding Source: {CORRESPONDING_SOURCE_URL}"
                                ));
                                ui.label(
                                    RichText::new(
                                        "Includes bundled BSD-2-Clause geometry/solver modules, SIL OFL 1.1 EB Garamond font, \
                                         and THIRD_PARTY_NOTICES.txt.",
                                    )
                                    .size(12.0)
                                    .color(muted),
                                );
                            });
                        ui.add_space(10.0);
                        ui.label(
                            RichText::new("Your projects and simulation runs remain untouched in Documents/nyx.")
                                .color(success),
                        );
                    }
                    Step::Options => {
                        if self.mode == WizardMode::UninstallOrRepair {
                            egui::Frame::new()
                                .fill(card_fill)
                                .stroke(Stroke::new(1.0, sep))
                                .corner_radius(CornerRadius::same(10))
                                .inner_margin(Margin::same(14))
                                .show(ui, |ui| {
                                    ui.set_width(ui.available_width());
                                    ui.label(RichText::new("ACTION").size(12.0).strong().color(muted));
                                    ui.radio_value(&mut self.uninstall_repair, false, "Uninstall nyx from this machine");
                                    ui.radio_value(&mut self.uninstall_repair, true, "Repair / verify installed files and integrations");
                                    ui.add_space(8.0);
                                    ui.label(
                                        RichText::new("Your simulations and workflows in Documents/nyx will NOT be deleted.")
                                            .color(success)
                                            .strong(),
                                    );
                                    if !self.uninstall_repair {
                                        ui.add_space(6.0);
                                        ui.checkbox(
                                            &mut self.purge_cache,
                                            "Clear temporary cache (%TEMP%/nyx-*) & solve_history.jsonl (--purge-cache)",
                                        );
                                    }
                                });
                            return;
                        }

                        egui::Frame::new()
                            .fill(card_fill)
                            .stroke(Stroke::new(1.0, sep))
                            .corner_radius(CornerRadius::same(10))
                            .inner_margin(Margin::same(14))
                            .show(ui, |ui| {
                                ui.set_width(ui.available_width());
                                ui.label(RichText::new("INSTALL LOCATION").size(12.0).strong().color(muted));
                                ui.horizontal(|ui| {
                                    let prev_scope = self.scope;
                                    ui.radio_value(&mut self.scope, InstallScope::User, "Current user (recommended)");
                                    ui.radio_value(&mut self.scope, InstallScope::System, "All users (requires Admin/root)");
                                    if self.scope != prev_scope {
                                        self.sync_default_prefix_for_scope();
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.add(egui::TextEdit::singleline(&mut self.custom_prefix).desired_width(420.0));
                                    if ui.button("Browse...").clicked() {
                                        if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                                            self.custom_prefix = dir.join("nyx").display().to_string();
                                        }
                                    }
                                });
                            });
                        ui.add_space(10.0);

                        egui::Frame::new()
                            .fill(card_fill)
                            .stroke(Stroke::new(1.0, sep))
                            .corner_radius(CornerRadius::same(10))
                            .inner_margin(Margin::same(14))
                            .show(ui, |ui| {
                                ui.set_width(ui.available_width());
                                ui.label(RichText::new("INTEGRATIONS").size(12.0).strong().color(muted));
                                ui.checkbox(&mut self.add_path, "Add nyx commands to PATH");
                                ui.checkbox(&mut self.desktop, "Create Start Menu / application launcher entry");
                                ui.checkbox(
                                    &mut self.file_assoc,
                                    "Register file associations (.nyx, .nys default; Open With for .nyx.json, .h5)",
                                );
                                #[cfg(windows)]
                                ui.checkbox(&mut self.desktop_shortcut, "Create Desktop shortcut");
                            });
                        ui.add_space(10.0);

                        egui::Frame::new()
                            .fill(card_fill)
                            .stroke(Stroke::new(1.0, sep))
                            .corner_radius(CornerRadius::same(10))
                            .inner_margin(Margin::same(14))
                            .show(ui, |ui| {
                                ui.set_width(ui.available_width());
                                ui.label(RichText::new("AI HOSTS (MCP SERVER)").size(12.0).strong().color(muted));
                                for (id, name, status, _cfg) in &self.host_rows {
                                    match status {
                                        Status::Unreadable { why } => {
                                            ui.horizontal(|ui| {
                                                ui.add_enabled(false, egui::Checkbox::new(&mut false, format!("{name} (custom/commented JSON)")));
                                                if ui.small_button("Copy JSON snippet").on_hover_text(why).clicked() {
                                                    let snip = connect_hook::status_label(status);
                                                    let _ = snip;
                                                    let text = nyx_connect::snippet(&plain_path(self.mcp_preview_exe.clone()));
                                                    ui.ctx().copy_text(text.clone());
                                                    self.copied_snippet = Some(format!("Copied MCP snippet for {name}"));
                                                }
                                            });
                                        }
                                        Status::NotInstalled => {
                                            if let Some(val) = self.host_selected.get_mut(id) {
                                                ui.checkbox(val, format!("{name} (not detected)"));
                                            }
                                        }
                                        _ => {
                                            if let Some(val) = self.host_selected.get_mut(id) {
                                                ui.checkbox(val, format!("{name} ({})", connect_hook::status_label(status)));
                                            }
                                        }
                                    }
                                }
                                if let Some(msg) = &self.copied_snippet {
                                    ui.label(RichText::new(msg).color(success).size(12.0));
                                }
                            });
                    }
                    Step::Progress => {
                        ui.label(RichText::new(&self.progress_label).strong().color(text_col));
                        ui.add_space(8.0);
                        ui.add(egui::ProgressBar::new(self.progress).show_percentage());
                        ui.add_space(8.0);
                        ui.checkbox(&mut self.show_logs, "Show transaction log");
                        if self.show_logs {
                            egui::ScrollArea::vertical().max_height(200.0).show(ui, |ui| {
                                for line in &self.logs {
                                    ui.label(RichText::new(line).monospace().size(12.0).color(muted));
                                }
                            });
                        }
                    }
                    Step::Done => {
                        if let Some(err) = &self.error {
                            ui.label(RichText::new("Operation did not complete").size(18.0).strong().color(danger));
                            ui.add_space(6.0);
                            ui.label(RichText::new(err).monospace().color(danger));
                        } else if let Some(rep) = &self.report {
                            ui.label(
                                RichText::new(format!("nyx {} finished cleanly", rep.action))
                                    .size(18.0)
                                    .strong()
                                    .color(success),
                            );
                            ui.add_space(8.0);
                            ui.label(RichText::new(format!("Location: {}", rep.install_root)).monospace());
                            for h in &rep.ai_hosts {
                                if h.action == "connected" || h.action == "disconnected" {
                                    ui.label(format!("  • {}: {}", h.id, h.action));
                                }
                            }
                            for next in &rep.next_steps {
                                ui.label(RichText::new(next).color(muted));
                            }
                        }
                    }
                });
            });
    }
}
