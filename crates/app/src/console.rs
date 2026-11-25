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
    use std::os::windows::io::{IntoRawHandle, RawHandle};
    type Handle = *mut core::ffi::c_void;
    const _: () = assert!(std::mem::size_of::<Handle>() == std::mem::size_of::<RawHandle>());
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

    let missing = |h: Handle| h.is_null() || h as isize == -1;

    // SAFETY: `AttachConsole(ATTACH_PARENT_PROCESS)` uses the Win32 sentinel `0xFFFF_FFFF` with
    // matching `extern "system"` ABI (`BOOL AttachConsole(DWORD)`), takes no pointers, and returns
    // zero if the parent process has no console.
    if unsafe { AttachConsole(ATTACH_PARENT_PROCESS) } == 0 {
        return;
    }

    // SAFETY: `GetStdHandle` is called only with valid Win32 `STD_*_HANDLE` constants. Any
    // returned handle is checked by `missing` (non-null and != `INVALID_HANDLE_VALUE`) before use.
    let out_missing = missing(unsafe { GetStdHandle(STD_OUTPUT_HANDLE) });
    // SAFETY: Same valid Win32 `STD_ERROR_HANDLE` constant as above.
    let err_missing = missing(unsafe { GetStdHandle(STD_ERROR_HANDLE) });

    if out_missing {
        if let Ok(f) = std::fs::OpenOptions::new().write(true).open("CONOUT$") {
            let h = f.into_raw_handle() as Handle;
            if !missing(h) {
                // SAFETY: `h` is a live, non-null, non-`INVALID_HANDLE_VALUE` OS file handle from
                // `OpenOptions::open("CONOUT$")`. Leaking ownership via `into_raw_handle` is
                // intentional so the Win32 standard output slot remains valid for process lifetime.
                let _ = unsafe { SetStdHandle(STD_OUTPUT_HANDLE, h) };
                if err_missing {
                    // SAFETY: Reusing the same live process-lifetime `CONOUT$` handle for stder
                    // when stderr was also unredirected (`STD_ERROR_HANDLE` is valid).
                    let _ = unsafe { SetStdHandle(STD_ERROR_HANDLE, h) };
                }
            }
        }
    } else if err_missing {
        if let Ok(f) = std::fs::OpenOptions::new().write(true).open("CONOUT$") {
            let h = f.into_raw_handle() as Handle;
            if !missing(h) {
                // SAFETY: `h` is a live `CONOUT$` OS handle; leaked intentionally for process lifetime
                // as the Win32 standard error handle.
                let _ = unsafe { SetStdHandle(STD_ERROR_HANDLE, h) };
            }
        }
    }

    // SAFETY: `STD_INPUT_HANDLE` is a valid Win32 `GetStdHandle` selector.
    if missing(unsafe { GetStdHandle(STD_INPUT_HANDLE) }) {
        if let Ok(f) = std::fs::OpenOptions::new().read(true).open("CONIN$") {
            let h = f.into_raw_handle() as Handle;
            if !missing(h) {
                // SAFETY: `h` is a live `CONIN$` OS handle; leaked intentionally for process lifetime
                // as the Win32 standard input handle.
                let _ = unsafe { SetStdHandle(STD_INPUT_HANDLE, h) };
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
