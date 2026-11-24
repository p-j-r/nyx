// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The actions of the window as one table (ADR 88): the "Do something" list of the command palette and the native menu bar are two views of it, so
//! a command cannot be in one and missing from the other. The native menu (`muda`, Windows for now) shows the shortcuts; the keys themselves are
//! handled by the window (`shell_shortcuts`, and egui's own Ctrl +/- for the zoom), so nothing fires twice.

use crate::app::App;
use crate::shell::{Page, Sheet, Sub};
use crate::state::Tab;
use eframe::egui;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Act {
    New,
    Open,
    ScriptEditor,
    ExportScript,
    CloseTab,
    Quit,
    Palette,
    Dashboard,
    ShowGeometry,
    ShowMaterials,
    ShowPhysics,
    ShowStudy,
    ShowSweep,
    ShowResults,
    ShowRuns,
    Settings,
    ZoomIn,
    ZoomOut,
    ZoomReset,
    Compute,
    BuildMesh,
    CheckConvergence,
    Stop,
    Versions,
    Share,
    Notes,
    OpenFolder,
    About,
}

pub struct Item {
    /// the top-level menu it belongs to
    pub menu: &'static str,
    pub label: &'static str,
    /// shown next to the label ("Ctrl+N"); empty when there is none
    pub keys: &'static str,
    pub act: Act,
}

/// Menus in order; the items of a menu in the order of this table.
pub const ACTIONS: [Item; 28] = [
    Item { menu: "File", label: "New workflow", keys: "Ctrl+N", act: Act::New },
    Item { menu: "File", label: "Open a .nyx file", keys: "Ctrl+O", act: Act::Open },
    Item { menu: "File", label: "Script (.nys)...", keys: "Ctrl+Shift+R", act: Act::ScriptEditor },
    Item { menu: "File", label: "Export as .nys script...", keys: "", act: Act::ExportScript },
    Item { menu: "File", label: "Close tab", keys: "Ctrl+W", act: Act::CloseTab },
    Item { menu: "File", label: "Quit", keys: "", act: Act::Quit },
    Item { menu: "Edit", label: "Search or run anything", keys: "Ctrl+K", act: Act::Palette },
    Item { menu: "View", label: "Dashboard", keys: "", act: Act::Dashboard },
    Item { menu: "View", label: "Geometry", keys: "Alt+1", act: Act::ShowGeometry },
    Item { menu: "View", label: "Materials", keys: "Alt+2", act: Act::ShowMaterials },
    Item { menu: "View", label: "Physics and conditions", keys: "Alt+3", act: Act::ShowPhysics },
    Item { menu: "View", label: "Study and solver", keys: "Alt+4", act: Act::ShowStudy },
    Item { menu: "View", label: "Sweep", keys: "Alt+5", act: Act::ShowSweep },
    Item { menu: "View", label: "Results", keys: "Alt+6", act: Act::ShowResults },
    Item { menu: "View", label: "Runs, plots and compare", keys: "", act: Act::ShowRuns },
    Item { menu: "View", label: "Settings", keys: "", act: Act::Settings },
    Item { menu: "View", label: "Zoom in", keys: "Ctrl+=", act: Act::ZoomIn },
    Item { menu: "View", label: "Zoom out", keys: "Ctrl+-", act: Act::ZoomOut },
    Item { menu: "View", label: "Actual size", keys: "Ctrl+0", act: Act::ZoomReset },
    Item { menu: "Run", label: "Compute the open experiment", keys: "F5", act: Act::Compute },
    Item { menu: "Run", label: "Build mesh", keys: "", act: Act::BuildMesh },
    Item { menu: "Run", label: "Check mesh convergence", keys: "", act: Act::CheckConvergence },
    Item { menu: "Run", label: "Stop the run", keys: "Esc", act: Act::Stop },
    Item { menu: "Workflow", label: "Versions...", keys: "", act: Act::Versions },
    Item { menu: "Workflow", label: "Share / Export .nyx...", keys: "", act: Act::Share },
    Item { menu: "Workflow", label: "Workflow notes", keys: "", act: Act::Notes },
    Item { menu: "Workflow", label: "Open workflow folder", keys: "", act: Act::OpenFolder },
    Item { menu: "Help", label: "About nyx and licences", keys: "", act: Act::About },
];

/// The top-level menus in order.
pub fn menus() -> Vec<&'static str> {
    let mut m: Vec<&'static str> = vec![];
    for i in &ACTIONS {
        if !m.contains(&i.menu) {
            m.push(i.menu);
        }
    }
    m
}

impl App {
    /// Shows a tab of the left panel of the open workflow (the one of the last tab of the window when it was on the Dashboard), with the centre on the model or the runs.
    fn show_tab(&mut self, tab: Tab, sub: Sub, sweep: bool) {
        self.shell.sheet = Sheet::None;
        if let Some(id) = self.ensure_active_workflow() {
            if let Some(w) = self.shell.open.get_mut(&id) {
                w.tab = tab;
                w.sub = sub;
                w.sweep.shown = sweep;
            }
            self.tab = tab;
            self.picking = None;
        }
    }

    fn ensure_active_workflow(&mut self) -> Option<String> {
        if let Page::Workflow(id) = self.shell.page.clone() {
            return Some(id);
        }
        if let Some(id) = self.shell.tabs.last().cloned() {
            self.show_workflow(&id);
            return Some(id);
        }
        self.status = "Open or create a workflow first (Ctrl+N).".into();
        None
    }

    /// In-window themed menu bar generated from `ACTIONS`, plus the Connect menu.
    pub fn bar_ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let mut chosen: Option<Act> = None;
        egui::MenuBar::new().ui(ui, |ui| {
            for name in menus() {
                ui.menu_button(name, |ui| {
                    ui.set_min_width(240.0);
                    for it in ACTIONS.iter().filter(|it| it.menu == name) {
                        let text = if it.keys.is_empty() {
                            it.label.to_string()
                        } else {
                            format!("{:<34}  {}", it.label, it.keys)
                        };
                        if ui.button(text).clicked() {
                            chosen = Some(it.act);
                            ui.close();
                        }
                    }
                });
            }
            let connect_menu = ui.menu_button("Connect", |ui| self.connect_menu(ui));
            if connect_menu.response.clicked() {
                self.refresh_connect();
            }
        });
        if let Some(act) = chosen {
            self.do_act(&ctx, act);
        }
    }

    /// Does what a menu entry or a palette entry says.
    pub fn do_act(&mut self, ctx: &egui::Context, act: Act) {
        match act {
            Act::New => self.shell.sheet = Sheet::New,
            Act::Open => {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
            }
            Act::ScriptEditor => self.open_nys_script(),
            Act::ExportScript => {
                self.shell.sheet = Sheet::None;
                self.export_script();
            }
            Act::CloseTab => {
                if let Page::Workflow(id) = self.shell.page.clone() {
                    self.close_tab(&id);
                }
            }
            Act::Quit => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
            Act::Palette => {
                self.shell.sheet = Sheet::Palette;
                self.shell.palette_q.clear();
            }
            Act::Dashboard => {
                self.shell.sheet = Sheet::None;
                self.flush();
                self.shell.page = Page::Dashboard;
            }
            Act::ShowGeometry => self.show_tab(Tab::Geometry, Sub::Main, false),
            Act::ShowMaterials => self.show_tab(Tab::Materials, Sub::Main, false),
            Act::ShowPhysics => self.show_tab(Tab::Physics, Sub::Main, false),
            Act::ShowStudy => self.show_tab(Tab::Study, Sub::Main, false),
            Act::ShowSweep => self.show_tab(Tab::Study, Sub::Main, true),
            Act::ShowResults => self.show_tab(Tab::Results, Sub::Main, false),
            Act::ShowRuns => self.show_tab(Tab::Results, Sub::Runs, false),
            Act::Settings => {
                self.shell.sheet = Sheet::None;
                self.shell.page = Page::Settings;
            }
            Act::ZoomIn => ctx.set_zoom_factor((ctx.zoom_factor() * 1.1).clamp(0.5, 3.0)),
            Act::ZoomOut => ctx.set_zoom_factor((ctx.zoom_factor() / 1.1).clamp(0.5, 3.0)),
            Act::ZoomReset => ctx.set_zoom_factor(1.0),
            Act::Compute => {
                self.shell.sheet = Sheet::None;
                let _ = self.ensure_active_workflow();
                self.shell_compute();
            }
            Act::BuildMesh => {
                self.shell.sheet = Sheet::None;
                self.rebuild_mesh();
            }
            Act::CheckConvergence => {
                self.shell.sheet = Sheet::None;
                self.check_convergence();
            }
            Act::Stop => {
                if self.busy.is_some() {
                    self.shell.confirm_cancel = true;
                }
            }
            Act::Versions => {
                if let Some(id) = self.ensure_active_workflow() {
                    if self.shell.open.get(&id).is_some_and(|w| w.dir.is_some()) {
                        self.shell.sheet = Sheet::Versions;
                        self.shell.version_sel = None;
                    } else {
                        self.shell.sheet = Sheet::None;
                        self.status = "Versions start with the first saved change.".into();
                    }
                }
            }
            Act::Share => {
                if self.ensure_active_workflow().is_some() {
                    self.shell.sheet = Sheet::Share;
                }
            }
            Act::Notes => {
                self.shell.sheet = Sheet::None;
                if let Some(id) = self.ensure_active_workflow() {
                    if let Some(w) = self.shell.open.get_mut(&id) {
                        w.notes_open = !w.notes_open;
                    }
                }
            }
            Act::OpenFolder => {
                self.shell.sheet = Sheet::None;
                if let Some(id) = self.ensure_active_workflow() {
                    if let Some(d) = self.shell.open.get(&id).and_then(|w| w.dir.as_ref()).map(|d| d.root.clone()) {
                        self.reveal(&d);
                    } else {
                        self.status = "Nothing is saved yet: the folder appears with the first change.".into();
                    }
                }
            }
            Act::About => self.shell.sheet = Sheet::About,
        }
    }
}

/// The native menu bar (Windows). Other platforms use the command palette until their `muda` integration is tried: GTK on Linux and the main
/// thread on macOS are the risks recorded in docs/ux/UX_MERGED_PLAN.md (R2).
#[cfg(windows)]
pub struct Native {
    _menu: muda::Menu,
    /// ids (indices into `ACTIONS`) clicked since the last frame
    clicked: std::sync::Arc<std::sync::Mutex<Vec<usize>>>,
}

#[cfg(windows)]
impl Native {
    /// Builds the menus from `ACTIONS` and attaches them to the window of `cc`. `None` when the window has no Win32 handle or muda refuses.
    pub fn install(cc: &eframe::CreationContext<'_>) -> Option<Native> {
        use muda::{accelerator::Accelerator, Menu, MenuEvent, MenuItem, Submenu};
        use raw_window_handle::{HasWindowHandle, RawWindowHandle};
        let RawWindowHandle::Win32(h) = cc.window_handle().ok()?.as_raw() else { return None };
        let menu = Menu::new();
        for name in menus() {
            let sub = Submenu::new(name, true);
            for (i, it) in ACTIONS.iter().enumerate().filter(|(_, it)| it.menu == name) {
                let accel = if it.keys.is_empty() { None } else { it.keys.parse::<Accelerator>().ok() };
                sub.append(&MenuItem::with_id(i.to_string(), it.label, true, accel)).ok()?;
            }
            menu.append(&sub).ok()?;
        }
        let clicked = std::sync::Arc::new(std::sync::Mutex::new(vec![]));
        let (q, ctx) = (clicked.clone(), cc.egui_ctx.clone());
        MenuEvent::set_event_handler(Some(move |e: MenuEvent| {
            if let Ok(i) = e.id().as_ref().parse::<usize>() {
                if let Ok(mut v) = q.lock() {
                    v.push(i);
                }
                ctx.request_repaint();
            }
        }));
        // SAFETY: the handle comes from the live window of this creation context.
        unsafe { menu.init_for_hwnd_with_theme(h.hwnd.get(), muda::MenuTheme::Auto).ok()? };
        Some(Native { _menu: menu, clicked })
    }

    /// The actions clicked in the menu since the last call.
    pub fn take(&self) -> Vec<Act> {
        let v: Vec<usize> = self.clicked.lock().map(|mut v| std::mem::take(&mut *v)).unwrap_or_default();
        v.into_iter().filter_map(|i| ACTIONS.get(i).map(|a| a.act)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_has_every_menu_and_no_empty_label() {
        assert_eq!(menus(), vec!["File", "Edit", "View", "Run", "Workflow", "Help"]);
        assert!(ACTIONS.iter().all(|a| !a.label.is_empty()));
        // each action appears once: the palette and the menu show one entry per action
        for (i, a) in ACTIONS.iter().enumerate() {
            assert!(ACTIONS.iter().skip(i + 1).all(|b| b.act != a.act), "{:?} twice", a.act);
        }
    }

    #[test]
    fn menu_navigation_actions_switch_sections_and_tabs() {
        let ctx = egui::Context::default();
        let mut a = App::with_context(ctx.clone());
        a.new_draft("NavTest", crate::shell::Start::Blank { physics: "ems".into(), three_d: false });
        let wid = match &a.shell.page {
            Page::Workflow(id) => id.clone(),
            other => panic!("expected workflow page, got {other:?}"),
        };
        a.do_act(&ctx, Act::ShowStudy);
        assert_eq!(a.shell.open[&wid].tab, Tab::Study);
        assert!(!a.shell.open[&wid].sweep.shown);

        a.do_act(&ctx, Act::ShowSweep);
        assert_eq!(a.shell.open[&wid].tab, Tab::Study);
        assert!(a.shell.open[&wid].sweep.shown);

        a.do_act(&ctx, Act::ShowMaterials);
        assert_eq!(a.shell.open[&wid].tab, Tab::Materials);

        a.do_act(&ctx, Act::ShowResults);
        assert_eq!(a.shell.open[&wid].tab, Tab::Results);
        a.do_act(&ctx, Act::ShowRuns);
        assert_eq!(a.shell.open[&wid].sub, Sub::Runs);
        a.do_act(&ctx, Act::ShowPhysics);
        assert_eq!(a.shell.open[&wid].sub, Sub::Main, "a tab of the model leaves the runs");

        a.do_act(&ctx, Act::Dashboard);
        assert_eq!(a.shell.page, Page::Dashboard);
        // Running a workflow navigation action from Dashboard jumps back to the open workflow tab
        a.do_act(&ctx, Act::ShowGeometry);
        assert_eq!(a.shell.page, Page::Workflow(wid.clone()));
        assert_eq!(a.shell.open[&wid].tab, Tab::Geometry);

        a.do_act(&ctx, Act::ScriptEditor);
        assert!(a.nys.open);
        assert!(!a.nys.code.trim().is_empty());
    }
}
