//! nyx: a desktop application for finite-element simulations (window with `eframe`/`egui`, solvers in
//! the engine and the physics crates).
//!
//! Command line: `nyx [FILE]` opens the window (and a project file), `--size WxH` sets the window size,
//! `api list` and `api show COMMAND` print the command reference (see `apicmd`), `check FILE` and `run FILE` check and run a JSON-lines list of engine commands without a window (see `runfile`), `--script "cmd; cmd; ..."` drives the window for tests and screenshots (see `script`), `--batch PROJECT.nyx.json` solves
//! without a window (see `batch`). Modules: `app` (state and layout), the Connect menu (in `app`, backed by the `nyx-connect` crate), `tabs` (the inspector), `viewport` and `scene` (drawing),
//! `pyscript` (the optional Python script), `state` (the mesh and result on show), `examples`, `colormap`, `theme`.
mod activity;
mod apicmd;
mod app;
mod batch;
mod colormap;
mod cvd;
mod dashboard;
mod examples;
mod feedback;
mod fmt;
mod forms;
mod geomui;
mod kit;
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
    let mut viewport = eframe::egui::ViewportBuilder::default().with_inner_size(size).with_min_inner_size([900.0, 600.0]).with_title("nyx").with_drag_and_drop(true);
    if let Some(icon) = window_icon() {
        viewport = viewport.with_icon(icon);
    }
    let options = eframe::NativeOptions { viewport, ..Default::default() };
    eframe::run_native(
        "nyx",
        options,
        Box::new(move |cc| {
            let mut a = app::App::new(cc);
            #[cfg(windows)]
            if script.is_none() {
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
    )
}

/// The nyx mark as the icon of the window (title bar, taskbar, Alt+Tab): the same picture as the top bar, scaled once with a good filter.
fn window_icon() -> Option<eframe::egui::IconData> {
    let img = image::load_from_memory(include_bytes!("../../../docs/brand/nyx-mark.png")).ok()?.to_rgba8();
    let img = image::imageops::resize(&img, 256, 256, image::imageops::FilterType::Lanczos3);
    Some(eframe::egui::IconData { rgba: img.into_raw(), width: 256, height: 256 })
}
