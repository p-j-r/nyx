// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx: a desktop application for finite-element simulations (window with `eframe`/`egui`, solvers in
//! the engine and the physics crates).
//!
//! Command line: `nyx [FILE]` opens the window (and a project file), `--size WxH` sets the window size,
//! `api list` and `api show COMMAND` print the command reference (see `apicmd`), `check FILE` and `run FILE` check and run a JSON-lines list of engine commands without a window (see `runfile`), `--script "cmd; cmd; ..."` drives the window for tests and screenshots (see `script`), `--batch PROJECT.nyx.json` solves
//! without a window (see `batch`). Modules: `app` (state and layout), the Connect menu (in `app`, backed by the `nyx-connect` crate), `tabs` (the inspector), `viewport` and `scene` (drawing),
//! `pyscript` (the optional Python script), `state` (the mesh and result on show), `examples`, `colormap`, `theme`.
//! On Windows the program is of the GUI subsystem: no console opens with the window; the command-line modes join the terminal that started them (`console`, ADR 123).
#![cfg_attr(windows, windows_subsystem = "windows")]
mod activity;
mod apicmd;
mod app;
mod batch;
mod colormap;
mod console;
mod cvd;
mod dashboard;
mod examples;
mod feedback;
mod fmt;
mod forms;
mod geomui;
mod kit;
mod matcmd;
mod menu;
mod pyscript;
mod regions;
mod runfile;
mod scene;
mod script;
mod sheets;
mod shell;
mod spans;
mod state;
mod sweepui;
mod tabs;
mod theme;
mod tokens;
mod viewport;
mod wfview;

fn main() -> eframe::Result {
    nyx_engine::process::opt_out_of_power_throttling();
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    args.retain(|a| !a.starts_with("-psn_"));
    if console::wants_console(&args) {
        console::attach_parent();
    }
    if matches!(args.first().map(String::as_str), Some("--version" | "-V")) {
        println!("nyx {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if matches!(args.first().map(String::as_str), Some("--help" | "-h")) {
        println!(
            "nyx {} — finite-element photonics, phononics & electromagnetics\n\n\
             Usage:\n\
               nyx [FILE]                          Open the desktop application (optionally loading .nyx.json, .nyx, .nys, .h5)\n\
               nyx --batch FILE [--physics LIST]   Mesh and solve without a window\n\
               nyx check FILE [--json]             Validate a command list or .nys script without running\n\
               nyx run FILE [OPTIONS]              Execute a command list or .nys script headlessly\n\
               nyx api <list|show COMMAND>         Print the engine command catalog reference\n\
               nyx material <SUBCOMMAND>           Inspect, evaluate, or import materials\n\
               nyx queue <SUBCOMMAND>              Manage the shared workspace solver queue\n\
               nyx --version | --help              Print version or this help message",
            env!("CARGO_PKG_VERSION")
        );
        return Ok(());
    }
    // started by an engine as its worker (`--job FOLDER`, see `nyx_engine::worker`): solve, answer on stdout, never open a window
    if matches!(args.first().map(String::as_str), Some("--job" | "--worker")) {
        std::process::exit(nyx_engine::worker::serve(&args));
    }
    if args.first().map(String::as_str) == Some("queue") {
        std::process::exit(nyx_engine::shared_queue::cli(&args[1..]));
    }
    if args.first().map(|a| a.as_str()) == Some("--batch") {
        if let Err(e) = batch::run(&args[1..]) {
            eprintln!("nyx: {e}");
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(String::as_str) == Some("api") {
        if let Err((code, msg)) = apicmd::run(&args[1..]) {
            eprintln!("nyx: {msg}");
            std::process::exit(code);
        }
        return Ok(());
    }
    if args.first().map(String::as_str) == Some("material") {
        if let Err((code, msg)) = matcmd::run(&args[1..]) {
            eprintln!("nyx: {msg}");
            std::process::exit(code);
        }
        return Ok(());
    }
    if let Some(which) = args.first().map(String::as_str).filter(|a| matches!(*a, "check" | "run")) {
        let done = if which == "check" { runfile::check(&args[1..]) } else { runfile::run(&args[1..]) };
        if let Err((code, msg)) = done {
            eprintln!("nyx: {msg}");
            std::process::exit(code);
        }
        return Ok(());
    }
    let mut script = None;
    if let Some(i) = args.iter().position(|a| a == "--script") {
        if i + 1 < args.len() {
            script = Some(args[i + 1].clone());
            args.drain(i..i + 2);
        }
    }
    let mut size = [1400.0, 860.0];
    if let Some(i) = args.iter().position(|a| a == "--size") {
        if i + 1 < args.len() {
            let v: Vec<f32> = args[i + 1].split('x').filter_map(|x| x.parse().ok()).collect();
            if v.len() == 2 {
                size = [v[0], v[1]];
            }
            args.drain(i..i + 2);
        }
    }
    let arg = args.first().cloned();
    let script_has_scale = script
        .as_ref()
        .is_some_and(|s| s.split(';').any(|c| c.trim().starts_with("scale ")));
    let mut viewport = eframe::egui::ViewportBuilder::default()
        .with_inner_size(size)
        .with_min_inner_size([900.0, 600.0])
        .with_title("nyx")
        .with_drag_and_drop(true)
        .with_visible(true)
        .with_active(true)
        .with_transparent(false)
        .with_clamp_size_to_monitor_size(true);
    if let Some(icon) = window_icon() {
        viewport = viewport.with_icon(icon);
    }
    let mut options = eframe::NativeOptions { viewport, ..Default::default() };
    let wslg = shell::is_wsl();
    if wslg {
        // Explicit FIFO presentation with single-frame queue prevents stale swapchain/compositor borders on WSLg resize/maximize.
        options.wgpu_options.surface = eframe::SurfaceConfig {
            present_mode: eframe::wgpu::PresentMode::Fifo,
            desired_maximum_frame_latency: Some(1),
        };
    }
    #[cfg(target_os = "linux")]
    let backend_choice = configure_linux_backend(&mut options, wslg);
    let res = eframe::run_native(
        "nyx",
        options,
        Box::new(move |cc| {
            let mut a = app::App::new(cc);
            if !script_has_scale {
                apply_startup_window_geometry(cc, size, a.shell.last_zoom, script.is_none());
            } else {
                cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::Visible(true));
            }
            if wslg && script.is_none() {
                if let Some(rs) = &cc.wgpu_render_state {
                    let info = rs.adapter.get_info();
                    if info.device_type == eframe::wgpu::DeviceType::Cpu {
                        eprintln!(
                            "nyx: running under WSLg with software rendering ({}); graphics may be slower.",
                            info.name
                        );
                    }
                }
            }
            #[cfg(windows)]
            if script.is_none() && std::env::var_os("NYX_NATIVE_MENU").is_some() {
                a.shell.native_menu = menu::Native::install(cc);
            }
            a.solve_in_child_processes();
            if let Some(p) = &arg {
                a.open_path(std::path::PathBuf::from(p));
            }
            if let Some(s) = &script {
                a.script = s.split(';').map(|c| c.trim().to_string()).filter(|c| !c.is_empty()).collect();
            }
            Ok(Box::new(a))
        }),
    );
    #[cfg(target_os = "linux")]
    if let Err(err) = &res {
        if wslg {
            eprintln!(
                "nyx: GUI window failed on {backend_choice} under WSLg ({err}). Try WINIT_UNIX_BACKEND=x11 if Wayland is unstable."
            );
        }
    }
    res
}

#[cfg(target_os = "linux")]
fn configure_linux_backend(options: &mut eframe::NativeOptions, wslg: bool) -> &'static str {
    use winit::platform::wayland::EventLoopBuilderExtWayland;
    use winit::platform::x11::EventLoopBuilderExtX11;
    let req = std::env::var("WINIT_UNIX_BACKEND").ok().map(|s| s.to_ascii_lowercase());
    let has_x11 = std::env::var_os("DISPLAY").is_some_and(|v| !v.is_empty());
    let has_wayland = std::env::var_os("WAYLAND_DISPLAY").is_some_and(|v| !v.is_empty())
        || std::env::var_os("WAYLAND_SOCKET").is_some_and(|v| !v.is_empty());
    match req.as_deref() {
        Some("x11") => {
            options.event_loop_builder = Some(Box::new(|b| {
                b.with_x11();
            }));
            "x11"
        }
        Some("wayland") => {
            options.event_loop_builder = Some(Box::new(|b| {
                b.with_wayland();
            }));
            "wayland"
        }
        _ if wslg && has_x11 => {
            // WSLg's Wayland socket frequently fails with `Broken pipe (os error 32)`; X11/XWayland is reliable.
            options.event_loop_builder = Some(Box::new(|b| {
                b.with_x11();
            }));
            "x11"
        }
        _ if has_wayland => "wayland",
        _ => "default",
    }
}

fn apply_startup_window_geometry(
    cc: &eframe::CreationContext<'_>,
    req_logical: [f32; 2],
    scale: f32,
    allow_maximize: bool,
) {
    let mon_logical = cc.winit_window().and_then(|w| w.current_monitor()).map(|m| {
        let sf = m.scale_factor().max(0.5) as f32;
        let sz = m.size();
        [sz.width as f32 / sf, sz.height as f32 / sf]
    });
    let (inner, min_inner, exceeds) =
        shell::compute_startup_window_size(req_logical, [900.0, 600.0], scale, mon_logical);
    cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::MinInnerSize(min_inner.into()));
    cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::InnerSize(inner.into()));
    if let Some(w) = cc.winit_window() {
        let _ = w.request_inner_size(winit::dpi::LogicalSize::new(inner[0], inner[1]));
        if let (Some(m), Ok(pos)) = (w.current_monitor(), w.outer_position()) {
            let mpos = m.position();
            let msz = m.size();
            let sf = m.scale_factor().max(0.5) as f32;
            let phys_w = (inner[0] * sf).round() as i32;
            let phys_h = (inner[1] * sf).round() as i32;
            let max_x = mpos.x + (msz.width as i32 - phys_w).max(0);
            let max_y = mpos.y + (msz.height as i32 - phys_h).max(0);
            let nx = pos.x.clamp(mpos.x, max_x);
            let ny = pos.y.clamp(mpos.y, max_y);
            if nx != pos.x || ny != pos.y {
                w.set_outer_position(winit::dpi::PhysicalPosition::new(nx, ny));
            }
        }
        w.set_visible(true);
    }
    cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::Visible(true));
    if allow_maximize {
        cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::Focus);
        if exceeds {
            cc.egui_ctx.send_viewport_cmd(eframe::egui::ViewportCommand::Maximized(true));
        }
    }
}

/// The nyx mark as the icon of the window (title bar, taskbar, Alt+Tab): the same picture as the top bar, scaled once with a good filter.
fn window_icon() -> Option<eframe::egui::IconData> {
    let img = image::load_from_memory(include_bytes!("../../../docs/brand/nyx-mark.png")).ok()?.to_rgba8();
    let img = image::imageops::resize(&img, 256, 256, image::imageops::FilterType::Lanczos3);
    Some(eframe::egui::IconData { rgba: img.into_raw(), width: 256, height: 256 })
}
