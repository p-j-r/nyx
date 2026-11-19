//! The actions of the window as one table (ADR 88): the "Do something" list of the command palette and the native menu bar are two views of it, so
//! a command cannot be in one and missing from the other. The native menu (`muda`, Windows for now) shows the shortcuts; the keys themselves are
//! handled by the window (`shell_shortcuts`, and egui's own Ctrl +/- for the zoom), so nothing fires twice.

use crate::app::App;
use crate::shell::{Page, Sheet};
use eframe::egui;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Act {
    New,
    Open,
    CloseTab,
    Quit,
    Palette,
    Dashboard,
    Settings,
    ZoomIn,
    ZoomOut,
    ZoomReset,
    Compute,
    Stop,
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
pub const ACTIONS: [Item; 13] = [
    Item { menu: "File", label: "New workflow", keys: "Ctrl+N", act: Act::New },
    Item { menu: "File", label: "Open a .nyx file", keys: "Ctrl+O", act: Act::Open },
    Item { menu: "File", label: "Close tab", keys: "Ctrl+W", act: Act::CloseTab },
    Item { menu: "File", label: "Quit", keys: "", act: Act::Quit },
    Item { menu: "Edit", label: "Search or run anything", keys: "Ctrl+K", act: Act::Palette },
    Item { menu: "View", label: "Dashboard", keys: "", act: Act::Dashboard },
    Item { menu: "View", label: "Settings", keys: "", act: Act::Settings },
    Item { menu: "View", label: "Zoom in", keys: "Ctrl+=", act: Act::ZoomIn },
    Item { menu: "View", label: "Zoom out", keys: "Ctrl+-", act: Act::ZoomOut },
    Item { menu: "View", label: "Actual size", keys: "Ctrl+0", act: Act::ZoomReset },
    Item { menu: "Run", label: "Compute the open experiment", keys: "F5", act: Act::Compute },
    Item { menu: "Run", label: "Stop the run", keys: "Esc", act: Act::Stop },
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
    /// Does what a menu entry or a palette entry says.
    pub fn do_act(&mut self, ctx: &egui::Context, act: Act) {
        match act {
            Act::New => self.shell.sheet = Sheet::New,
            Act::Open => {
                self.shell.sheet = Sheet::None;
                self.pick_nyx_file();
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
            Act::Settings => {
                self.shell.sheet = Sheet::None;
                self.shell.page = Page::Settings;
            }
            Act::ZoomIn => ctx.set_zoom_factor((ctx.zoom_factor() * 1.1).clamp(0.5, 3.0)),
            Act::ZoomOut => ctx.set_zoom_factor((ctx.zoom_factor() / 1.1).clamp(0.5, 3.0)),
            Act::ZoomReset => ctx.set_zoom_factor(1.0),
            Act::Compute => {
                self.shell.sheet = Sheet::None;
                self.shell_compute();
            }
            Act::Stop => {
                if self.busy.is_some() {
                    self.shell.confirm_cancel = true;
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
        assert_eq!(menus(), vec!["File", "Edit", "View", "Run", "Help"]);
        assert!(ACTIONS.iter().all(|a| !a.label.is_empty()));
        // each action appears once: the palette and the menu show one entry per action
        for (i, a) in ACTIONS.iter().enumerate() {
            assert!(ACTIONS.iter().skip(i + 1).all(|b| b.act != a.act), "{:?} twice", a.act);
        }
    }
}
