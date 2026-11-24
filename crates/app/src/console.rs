// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The console of the command-line modes on Windows (ADR 123). `nyx.exe` is a program of the GUI subsystem, so starting it from Explorer, the Start
//! menu or a shortcut opens no console window. The command-line modes (`api`, `check`, `run`, `material`, `queue`, `--batch`, `--version`, `--help`)
//! print to the console of the terminal that started it: `attach_parent` joins that console and points the standard handles at it when the caller
//! did not redirect them. When the caller redirected them (a pipe, a file: the MCP server, a test, a script with `| Out-String`) they are left alone.
//! Other platforms have nothing to do: there the same executable always has a terminal or none, as it is started.

/// Whether the arguments ask for a command-line mode that prints (and so needs the console of the terminal).
pub fn wants_console(args: &[String]) -> bool {
    matches!(args.first().map(String::as_str), Some("api" | "check" | "run" | "material" | "queue" | "--batch" | "--version" | "-V" | "--help" | "-h"))
}

/// Joins the console of the parent process, when it has one, for the command-line modes. Never fails: without a console the output goes nowhere,
/// as it does for any program of the GUI subsystem.
#[cfg(windows)]
pub fn attach_parent() {
    use std::os::windows::io::IntoRawHandle;
    type Handle = *mut core::ffi::c_void;
    const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
    const STD_INPUT_HANDLE: u32 = -10i32 as u32;
    const STD_OUTPUT_HANDLE: u32 = -11i32 as u32;
    const STD_ERROR_HANDLE: u32 = -12i32 as u32;
    #[link(name = "kernel32")]
    extern "system" {
        fn AttachConsole(process_id: u32) -> i32;
        fn GetStdHandle(which: u32) -> Handle;
        fn SetStdHandle(which: u32, handle: Handle) -> i32;
    }
    // SAFETY: plain kernel32 calls with valid constants; the handles opened here are leaked on purpose (they live as long as the process).
    unsafe {
        if AttachConsole(ATTACH_PARENT_PROCESS) == 0 {
            return;
        }
        let missing = |h: Handle| h.is_null() || h as isize == -1;
        if missing(GetStdHandle(STD_OUTPUT_HANDLE)) {
            if let Ok(f) = std::fs::OpenOptions::new().write(true).open("CONOUT$") {
                let h = f.into_raw_handle();
                SetStdHandle(STD_OUTPUT_HANDLE, h as Handle);
                if missing(GetStdHandle(STD_ERROR_HANDLE)) {
                    SetStdHandle(STD_ERROR_HANDLE, h as Handle);
                }
            }
        } else if missing(GetStdHandle(STD_ERROR_HANDLE)) {
            if let Ok(f) = std::fs::OpenOptions::new().write(true).open("CONOUT$") {
                SetStdHandle(STD_ERROR_HANDLE, f.into_raw_handle() as Handle);
            }
        }
        if missing(GetStdHandle(STD_INPUT_HANDLE)) {
            if let Ok(f) = std::fs::OpenOptions::new().read(true).open("CONIN$") {
                SetStdHandle(STD_INPUT_HANDLE, f.into_raw_handle() as Handle);
            }
        }
    }
}

#[cfg(not(windows))]
pub fn attach_parent() {}

#[cfg(test)]
mod tests {
    use super::wants_console;

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_command_line_modes_get_the_console_and_the_window_does_not() {
        for a in [&["api", "list"][..], &["check", "x.nys"], &["run", "x.nys"], &["material", "list"], &["queue", "list"], &["--batch", "p.nyx.json"], &["--version"], &["--help"], &["-V"], &["-h"]] {
            assert!(wants_console(&args(a)), "{a:?}");
        }
        for a in [&[][..], &["model.nyx.json"], &["--size", "1400x860"], &["--script", "quit"], &["--job", "dir"], &["--worker"]] {
            assert!(!wants_console(&args(a)), "{a:?}");
        }
    }
}
