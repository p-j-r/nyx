// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Zero-dependency CLI parser, interactive TTY fallback, and subcommand dispatcher shared by `nyx-setup` and `nyx-uninstall`.

use crate::bundle::PayloadSource;
use crate::connect_hook::ConnectPolicy;
use crate::manifest::InstallScope;
#[cfg(feature = "gui")]
use crate::platform::windows::is_inside_dir;
use crate::platform::windows::{evacuate_to_temp_and_relaunch, run_reap_final};
use crate::platform::{native_backend, HostPaths};
use crate::process::OsProcessInspector;
use crate::report::InstallReport;
#[cfg(feature = "gui")]
use crate::tx::locate_manifest;
use crate::tx::{
    exit_code_of, run_doctor, run_install_or_upgrade, run_repair, run_status, run_uninstall,
    ExitError, TxOptions, TxProgress,
};
use anyhow::{bail, Result};
use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Action {
    Install,
    Upgrade,
    Repair,
    Uninstall,
    Status,
    Doctor,
    ReapFinal,
}

#[derive(Clone, Debug)]
pub struct ParsedArgs {
    pub action: Option<Action>,
    pub tx: TxOptions,
    pub bundle: Option<PathBuf>,
    pub from_dir: Option<PathBuf>,
    pub from_app: Option<PathBuf>,
    pub yes: bool,
    pub quiet: bool,
    pub json: bool,
    pub cli_only: bool,
    pub desktop_explicit: bool,
    pub wait_pid: Option<u32>,
    pub reap_root: Option<PathBuf>,
}

pub fn print_help(prog_name: &str) {
    println!(
        "{prog_name} {} — installer and uninstaller for nyx\n\n\
         Usage:\n\
           {prog_name} [install|upgrade|repair|uninstall|status|doctor] [OPTIONS]\n\n\
         Options:\n\
           --user                 Install for current user only (default)\n\
           --system               Install system-wide (requires Admin/root)\n\
           --prefix DIR           Install into custom directory DIR\n\
           --bundle FILE          Install from .nyxpack or .zip bundle instead of embedded payload\n\
           --from-dir DIR         Install from unpacked distribution directory\n\
           --from-app APP         Install macOS app bundle from existing Nyx.app\n\
           --no-path              Do not modify user/system PATH\n\
           --no-desktop           Do not create Start Menu / desktop entry\n\
           --desktop              Force desktop entry creation on headless Linux\n\
           --desktop-shortcut     Also create a Desktop shortcut (Windows)\n\
           --no-file-assoc        Do not register .nyx / .nys / OpenWith file associations\n\
           --file-assoc           Register file associations (default: enabled)\n\
           --connect SPEC         AI hosts: detected (default) | all | none | id1,id2\n\
           --purge-cache          (uninstall) Also clear temp cache & solve_history.jsonl; never touches Documents/nyx\n\
           --dry-run              Preview planned changes without writing anything\n\
           --force                Terminate running nyx processes / overwrite non-empty targets\n\
           --allow-root-user      Allow --user install when running as root on Unix\n\
           --cli                  Force terminal CLI/TTY mode (skip GUI wizard)\n\
           -y, --yes              Non-interactive confirmation\n\
           -q, --quiet            Minimal output\n\
           --json                 Emit machine-readable install_report.v1 JSON on stdout\n\
           -V, --version          Print version and exit\n\
           -h, --help             Print this help message and exit",
        env!("CARGO_PKG_VERSION")
    );
}

pub fn parse_args(args: &[String]) -> Result<ParsedArgs> {
    let mut action: Option<Action> = None;
    let mut tx = TxOptions::default_user();
    let mut bundle = None;
    let mut from_dir = None;
    let mut from_app = None;
    let mut yes = false;
    let mut quiet = false;
    let mut json = false;
    let mut cli_only = false;
    let mut wait_pid = None;
    let mut reap_root = None;

    let need_val = |i: usize, flag: &str| -> Result<String> {
        args.get(i + 1)
            .cloned()
            .ok_or_else(|| ExitError::new(2, format!("{flag} requires a value")).into())
    };

    let mut i = 0usize;
    while i < args.len() {
        let a = args[i].as_str();
        match a {
            "install" if action.is_none() => action = Some(Action::Install),
            "upgrade" if action.is_none() => action = Some(Action::Upgrade),
            "repair" if action.is_none() => action = Some(Action::Repair),
            "uninstall" | "remove" if action.is_none() => action = Some(Action::Uninstall),
            "status" if action.is_none() => action = Some(Action::Status),
            "doctor" if action.is_none() => action = Some(Action::Doctor),
            "--_reap-final" => action = Some(Action::ReapFinal),
            "--wait-pid" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ExitError::new(2, "--wait-pid requires a pid"))?;
                wait_pid = Some(
                    v.parse()
                        .map_err(|_| ExitError::new(2, "--wait-pid must be a u32"))?,
                );
            }
            "--root" => {
                reap_root = Some(PathBuf::from(need_val(i, "--root")?));
                i += 1;
            }
            "--user" => tx.scope = InstallScope::User,
            "--system" => tx.scope = InstallScope::System,
            "--prefix" => {
                tx.prefix = Some(PathBuf::from(need_val(i, "--prefix")?));
                i += 1;
            }
            "--bundle" => {
                bundle = Some(PathBuf::from(need_val(i, "--bundle")?));
                i += 1;
            }
            "--from-dir" => {
                from_dir = Some(PathBuf::from(need_val(i, "--from-dir")?));
                i += 1;
            }
            "--from-app" => {
                from_app = Some(PathBuf::from(need_val(i, "--from-app")?));
                i += 1;
            }
            "--no-path" => tx.add_path = false,
            "--no-desktop" => tx.desktop = false,
            "--desktop" => {
                tx.desktop = true;
                tx.desktop_explicit = true;
            }
            "--desktop-shortcut" => tx.desktop_shortcut = true,
            "--no-file-assoc" => tx.file_assoc = false,
            "--file-assoc" => tx.file_assoc = true,
            "--connect" => {
                tx.connect_policy = ConnectPolicy::parse(&need_val(i, "--connect")?);
                i += 1;
            }
            "--purge-cache" => tx.purge_cache = true,
            "--dry-run" => tx.dry_run = true,
            "--force" => tx.force = true,
            "--allow-root-user" => tx.allow_root_user = true,
            "--cli" => cli_only = true,
            "-y" | "--yes" => yes = true,
            "-q" | "--quiet" => quiet = true,
            "--json" => json = true,
            "--install" if action.is_none() => action = Some(Action::Install),
            "--uninstall" if action.is_none() => action = Some(Action::Uninstall),
            "--repair" if action.is_none() => action = Some(Action::Repair),
            other => {
                bail!(ExitError::new(
                    2,
                    format!("unknown argument '{other}' (run --help)")
                ));
            }
        }
        i += 1;
    }

    let desktop_explicit = tx.desktop_explicit;
    Ok(ParsedArgs {
        action,
        tx,
        bundle,
        from_dir,
        from_app,
        yes,
        quiet,
        json,
        cli_only,
        desktop_explicit,
        wait_pid,
        reap_root,
    })
}

fn invoked_as_uninstaller(argv0: Option<&Path>) -> bool {
    argv0
        .and_then(|p| p.file_stem())
        .and_then(|s| s.to_str())
        .is_some_and(|s| s.to_ascii_lowercase().starts_with("nyx-uninstall"))
}

fn has_gui_display(paths: &HostPaths) -> bool {
    if cfg!(any(windows, target_os = "macos")) {
        true
    } else {
        paths.display.is_some() || paths.wayland_display.is_some()
    }
}

fn print_human_report(r: &InstallReport, quiet: bool) {
    if quiet {
        return;
    }
    println!(
        "nyx {} {}: {} ({})",
        r.action,
        if r.dry_run { "preview" } else { "complete" },
        r.install_root,
        r.nyx_version.as_deref().unwrap_or("unknown")
    );
    if r.files_written > 0 {
        println!("  files installed: {}", r.files_written);
    }
    if r.files_removed > 0 {
        println!("  files removed:   {}", r.files_removed);
    }
    if r.path_updated {
        println!("  PATH updated:    yes");
    }
    if r.desktop_updated {
        println!("  desktop entries: updated");
    }
    if r.file_assoc_updated {
        println!("  file assocs:     updated");
    }
    for h in &r.ai_hosts {
        println!("  AI host {:<15} {} ({})", h.id, h.action, h.status);
    }
    for c in &r.doctor_checks {
        println!("  check {:<20} [{}] {}", c.name, c.status, c.detail);
    }
    for w in &r.warnings {
        eprintln!("warning: {w}");
    }
    for n in &r.next_steps {
        println!("next: {n}");
    }
}

fn confirm_tty(prompt: &str) -> bool {
    if !std::io::stdin().is_terminal() || !std::io::stdout().is_terminal() {
        return false;
    }
    print!("{prompt} [Y/n] ");
    let _ = std::io::stdout().flush();
    let mut line = String::new();
    if std::io::stdin().read_line(&mut line).is_err() {
        return false;
    }
    let ans = line.trim().to_ascii_lowercase();
    ans.is_empty() || ans == "y" || ans == "yes"
}

pub fn entry() -> i32 {
    let raw_args: Vec<String> = std::env::args()
        .skip(1)
        .filter(|a| !a.starts_with("-psn_"))
        .collect();
    let argv0 = std::env::args().next().map(PathBuf::from);
    let prog_name = argv0
        .as_deref()
        .and_then(|p| p.file_name())
        .and_then(|n| n.to_str())
        .unwrap_or("nyx-setup");

    if matches!(
        raw_args.first().map(String::as_str),
        Some("--version" | "-V")
    ) {
        println!("{prog_name} {}", env!("CARGO_PKG_VERSION"));
        return 0;
    }
    if matches!(raw_args.first().map(String::as_str), Some("--help" | "-h")) {
        print_help(prog_name);
        return 0;
    }

    let parsed = match parse_args(&raw_args) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{prog_name}: {e}");
            return exit_code_of(&e);
        }
    };

    if parsed.action == Some(Action::ReapFinal) {
        let Some(pid) = parsed.wait_pid else {
            eprintln!("--_reap-final requires --wait-pid");
            return 2;
        };
        let Some(root) = parsed.reap_root else {
            eprintln!("--_reap-final requires --root");
            return 2;
        };
        if let Err(e) = run_reap_final(pid, &root) {
            eprintln!("{e}");
            return 1;
        }
        return 0;
    }

    let backend = native_backend();
    backend.harden_process();
    let paths = HostPaths::detect();
    let current_exe = std::env::current_exe().unwrap_or_default();

    let source = match PayloadSource::detect(
        &current_exe,
        parsed.bundle.as_deref(),
        parsed.from_dir.as_deref(),
        parsed.from_app.as_deref(),
    ) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{prog_name}: {e}");
            return exit_code_of(&e);
        }
    };

    let uninstaller_mode = invoked_as_uninstaller(argv0.as_deref());
    let resolved_action = match parsed.action {
        Some(a) => a,
        None if uninstaller_mode => Action::Uninstall,
        None if source.is_some() => Action::Install,
        None => {
            print_help(prog_name);
            return 0;
        }
    };

    // Launch GUI wizard when compiled with `--features gui` and invoked interactively with display available.
    #[cfg(feature = "gui")]
    {
        let headless_flags =
            parsed.cli_only || parsed.yes || parsed.quiet || parsed.json || parsed.tx.dry_run;
        let interactive_action = matches!(
            resolved_action,
            Action::Install | Action::Upgrade | Action::Uninstall | Action::Repair
        );
        if !headless_flags && interactive_action && has_gui_display(&paths) {
            // Mode A on Windows: if launched as GUI from inside the installed root, evacuate to %TEMP% first!
            if cfg!(windows)
                && matches!(resolved_action, Action::Uninstall | Action::Repair)
                && locate_manifest(
                    parsed.tx.prefix.as_deref(),
                    Some(&current_exe),
                    backend.as_ref(),
                    parsed.tx.effective_scope(),
                    &paths,
                )
                .is_some_and(|(_, m)| is_inside_dir(&current_exe, Path::new(&m.install_root)))
            {
                match evacuate_to_temp_and_relaunch(&current_exe, &raw_args) {
                    Ok(true) => return 0,
                    Ok(false) => {}
                    Err(e) => {
                        eprintln!("{prog_name}: failed to evacuate uninstaller to TEMP: {e}");
                        return 1;
                    }
                }
            }
            if let Err(e) = crate::gui::run_wizard(resolved_action, parsed, source) {
                eprintln!("{prog_name}: GUI error: {e}");
                return 1;
            }
            return 0;
        }
    }
    let _ = (has_gui_display(&paths), evacuate_to_temp_and_relaunch);

    let inspector = OsProcessInspector;
    let mut progress_printer = |p: TxProgress| {
        if !parsed.quiet && !parsed.json {
            println!("[{}/{}] {}", p.phase, p.total_phases, p.message);
        }
    };

    let result: Result<InstallReport> = match resolved_action {
        Action::Install | Action::Upgrade => {
            let Some(ref src) = source else {
                eprintln!("{prog_name}: no embedded payload or --bundle / --from-dir / --from-app source found");
                return 2;
            };
            if !parsed.yes && !parsed.tx.dry_run {
                let target_preview = parsed.tx.prefix.clone().unwrap_or_else(|| {
                    backend.default_install_root(parsed.tx.effective_scope(), &paths)
                });
                if !confirm_tty(&format!("Install nyx into {}?", target_preview.display())) {
                    eprintln!("{prog_name}: aborted (pass --yes for non-interactive install)");
                    return 2;
                }
            }
            let act_str = if resolved_action == Action::Upgrade {
                "upgrade"
            } else {
                "install"
            };
            run_install_or_upgrade(
                act_str,
                src,
                &parsed.tx,
                backend.as_ref(),
                &inspector,
                &paths,
                &mut progress_printer,
            )
        }
        Action::Repair => {
            let Some(ref src) = source else {
                eprintln!("{prog_name}: repair requires an embedded payload or --bundle / --from-dir / --from-app");
                return 2;
            };
            if !parsed.yes && !parsed.tx.dry_run && !confirm_tty("Repair nyx installation?") {
                eprintln!("{prog_name}: aborted (pass --yes for non-interactive repair)");
                return 2;
            }
            run_repair(
                src,
                &parsed.tx,
                backend.as_ref(),
                &inspector,
                &paths,
                &mut progress_printer,
            )
        }
        Action::Uninstall => {
            if !parsed.yes
                && !parsed.tx.dry_run
                && !confirm_tty(
                    "Uninstall nyx? (Your projects in Documents/nyx will NOT be deleted)",
                )
            {
                eprintln!("{prog_name}: aborted (pass --yes for non-interactive uninstall)");
                return 2;
            }
            run_uninstall(
                &parsed.tx,
                backend.as_ref(),
                &inspector,
                &paths,
                &mut progress_printer,
            )
        }
        Action::Status => run_status(&parsed.tx, backend.as_ref(), &paths),
        Action::Doctor => run_doctor(&parsed.tx, backend.as_ref(), &paths),
        Action::ReapFinal => unreachable!(),
    };

    match result {
        Ok(report) => {
            if parsed.json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&report).unwrap_or_default()
                );
            } else {
                print_human_report(&report, parsed.quiet);
            }
            if report.ok {
                0
            } else {
                1
            }
        }
        Err(e) => {
            let code = exit_code_of(&e);
            if parsed.json {
                let mut err_report = InstallReport::new(
                    match resolved_action {
                        Action::Install => "install",
                        Action::Upgrade => "upgrade",
                        Action::Repair => "repair",
                        Action::Uninstall => "uninstall",
                        Action::Status => "status",
                        Action::Doctor => "doctor",
                        Action::ReapFinal => "uninstall",
                    },
                    parsed.tx.effective_scope(),
                    parsed
                        .tx
                        .prefix
                        .as_ref()
                        .map(|p| p.display().to_string())
                        .unwrap_or_default(),
                    parsed.tx.dry_run,
                );
                err_report.ok = false;
                err_report.errors.push(e.to_string());
                println!(
                    "{}",
                    serde_json::to_string_pretty(&err_report).unwrap_or_default()
                );
            } else {
                eprintln!("{prog_name}: error: {e}");
            }
            code
        }
    }
}
