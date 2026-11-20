// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Process-wide settings of the program that hosts the engine. Today one: on Windows, opt the process out of power throttling.
//!
//! Windows 11 puts a process it considers background work (a console program, a window that is hidden or not in front, a program started by
//! another one) into "efficiency mode" and runs its threads on the efficiency cores of a hybrid processor. A solve that starts on the
//! performance cores is then moved to the efficiency cores after its first seconds: on an i5-13500 the shifted factorisation of the 3D box
//! cavity took 2.2 s in the first solve of a process and 6.2 s in every later one (`docs/PRODUCTIVITY_STATUS.md`, finding 1). The call below asks
//! Windows to leave this process alone (`PROCESS_POWER_THROTTLING_EXECUTION_SPEED` switched off). It changes which core runs the work and when,
//! never a number: results do not depend on the core or the thread count (ADR 11). Elsewhere it does nothing.

/// Asks Windows not to throttle this process. True when the request was accepted; false when it was refused (an older Windows) or on another
/// system. Call once at the start of `main` of a program that hosts an `Engine`.
pub fn opt_out_of_power_throttling() -> bool {
    #[cfg(windows)]
    {
        #[repr(C)]
        struct PowerThrottlingState {
            version: u32,
            control_mask: u32,
            state_mask: u32,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn GetCurrentProcess() -> *mut std::ffi::c_void;
            fn SetProcessInformation(
                process: *mut std::ffi::c_void,
                class: i32,
                info: *const PowerThrottlingState,
                size: u32,
            ) -> i32;
        }
        // ProcessPowerThrottling = 4; version 1; control mask 1 = execution speed; state mask 0 = throttling off
        const PROCESS_POWER_THROTTLING: i32 = 4;
        let state = PowerThrottlingState { version: 1, control_mask: 1, state_mask: 0 };
        // SAFETY: the pseudo handle of the current process and a correctly sized, initialised struct that outlives the call.
        unsafe {
            SetProcessInformation(
                GetCurrentProcess(),
                PROCESS_POWER_THROTTLING,
                &state,
                std::mem::size_of::<PowerThrottlingState>() as u32,
            ) != 0
        }
    }
    #[cfg(not(windows))]
    {
        false
    }
}

/// The most memory this process has used so far, in bytes: `(peak working set, peak commit)`. The working set is what was in RAM, the commit what
/// was promised to the process (RAM plus paging file), and the number that decides whether another program fits next to it. Windows keeps both
/// exactly; on Linux the second is the peak of the virtual size. `None` where the system gives neither. A worker process (`worker.rs`) reports it
/// when its solve ends, which is how a run learns its memory: in the process of a window it would be the peak of everything the window did.
pub fn peak_memory() -> Option<(u64, u64)> {
    #[cfg(windows)]
    {
        #[repr(C)]
        #[derive(Default)]
        struct Counters {
            cb: u32,
            page_fault_count: u32,
            peak_working_set: usize,
            working_set: usize,
            quota_peak_paged_pool: usize,
            quota_paged_pool: usize,
            quota_peak_nonpaged_pool: usize,
            quota_nonpaged_pool: usize,
            pagefile_usage: usize,
            peak_pagefile_usage: usize,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn GetCurrentProcess() -> *mut std::ffi::c_void;
            fn K32GetProcessMemoryInfo(process: *mut std::ffi::c_void, counters: *mut Counters, cb: u32) -> i32;
        }
        let mut c = Counters { cb: std::mem::size_of::<Counters>() as u32, ..Counters::default() };
        // SAFETY: the pseudo handle of the current process and a correctly sized, initialised out-parameter that outlives the call.
        let ok = unsafe { K32GetProcessMemoryInfo(GetCurrentProcess(), &mut c, c.cb) };
        (ok != 0).then_some((c.peak_working_set as u64, c.peak_pagefile_usage as u64))
    }
    #[cfg(not(windows))]
    {
        let text = std::fs::read_to_string("/proc/self/status").ok()?;
        let kb = |key: &str| -> Option<u64> {
            text.lines().find(|l| l.starts_with(key))?.split_whitespace().nth(1)?.parse::<u64>().ok().map(|k| k * 1024)
        };
        Some((kb("VmHWM:")?, kb("VmPeak:")?))
    }
}

/// The physical memory of the machine and what is available now without paging, in bytes: `(total, available)`. Admission (`estimate.rs`)
/// compares a solve's estimate with it. `None` where the system does not say.
pub fn system_memory() -> Option<(u64, u64)> {
    #[cfg(windows)]
    {
        #[repr(C)]
        struct Status {
            length: u32,
            memory_load: u32,
            total_phys: u64,
            avail_phys: u64,
            total_page_file: u64,
            avail_page_file: u64,
            total_virtual: u64,
            avail_virtual: u64,
            avail_extended_virtual: u64,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn GlobalMemoryStatusEx(status: *mut Status) -> i32;
        }
        let mut s = Status {
            length: std::mem::size_of::<Status>() as u32,
            memory_load: 0,
            total_phys: 0,
            avail_phys: 0,
            total_page_file: 0,
            avail_page_file: 0,
            total_virtual: 0,
            avail_virtual: 0,
            avail_extended_virtual: 0,
        };
        // SAFETY: a correctly sized struct with its length field set, which outlives the call.
        let ok = unsafe { GlobalMemoryStatusEx(&mut s) };
        (ok != 0).then_some((s.total_phys, s.avail_phys))
    }
    #[cfg(not(windows))]
    {
        let text = std::fs::read_to_string("/proc/meminfo").ok()?;
        let kb = |key: &str| -> Option<u64> {
            text.lines().find(|l| l.starts_with(key))?.split_whitespace().nth(1)?.parse::<u64>().ok().map(|k| k * 1024)
        };
        Some((kb("MemTotal:")?, kb("MemAvailable:")?))
    }
}

/// Starts `program args...` as a process that belongs to nobody: it outlives the program that started it, has no console and no standard
/// streams, and **inherits no handle** of the starter's. The last matters as much as the first: a child that inherited the pipe a host reads the
/// starter's output from would keep that pipe open after the starter ended, and the host would wait for ever for its end. On Windows it is made
/// with `CreateProcessW` (inheritance off, which `std::process::Command` cannot say) outside the job object of the starter where the system lets
/// it (a window or terminal that ends kills the processes of its job); elsewhere in a process group of its own with the standard streams null.
pub fn spawn_detached(program: &std::path::Path, args: &[std::ffi::OsString]) -> std::io::Result<()> {
    #[cfg(windows)]
    {
        use std::ffi::c_void;
        use std::os::windows::ffi::OsStrExt;
        #[repr(C)]
        struct StartupInfoW {
            cb: u32,
            reserved: *mut u16,
            desktop: *mut u16,
            title: *mut u16,
            x: u32,
            y: u32,
            x_size: u32,
            y_size: u32,
            x_chars: u32,
            y_chars: u32,
            fill: u32,
            flags: u32,
            show: u16,
            reserved2: u16,
            reserved2_ptr: *mut u8,
            std_in: *mut c_void,
            std_out: *mut c_void,
            std_err: *mut c_void,
        }
        #[repr(C)]
        struct ProcessInformation {
            process: *mut c_void,
            thread: *mut c_void,
            pid: u32,
            tid: u32,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn CreateProcessW(
                application: *const u16,
                command_line: *mut u16,
                process_attributes: *const c_void,
                thread_attributes: *const c_void,
                inherit_handles: i32,
                flags: u32,
                environment: *const c_void,
                directory: *const u16,
                startup: *const StartupInfoW,
                info: *mut ProcessInformation,
            ) -> i32;
            fn CloseHandle(handle: *mut c_void) -> i32;
        }
        // the command line as the C runtime of the child will split it: every word quoted, quotes escaped, backslashes before a quote doubled
        fn quote(word: &std::ffi::OsStr, out: &mut Vec<u16>) {
            const BACKSLASH: u16 = b'\\' as u16;
            const QUOTE: u16 = b'"' as u16;
            out.push(QUOTE);
            let mut backslashes = 0usize;
            for c in word.encode_wide() {
                if c == BACKSLASH {
                    backslashes += 1;
                    continue;
                }
                if c == QUOTE {
                    out.extend(std::iter::repeat_n(BACKSLASH, backslashes * 2 + 1));
                } else {
                    out.extend(std::iter::repeat_n(BACKSLASH, backslashes));
                }
                backslashes = 0;
                out.push(c);
            }
            out.extend(std::iter::repeat_n(BACKSLASH, backslashes * 2));
            out.push(QUOTE);
        }
        let mut line: Vec<u16> = vec![];
        quote(program.as_os_str(), &mut line);
        for a in args {
            line.push(' ' as u16);
            quote(a, &mut line);
        }
        line.push(0);
        let app: Vec<u16> = program.as_os_str().encode_wide().chain(Some(0)).collect();
        // DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW, and CREATE_BREAKAWAY_FROM_JOB, which a job that does not allow it refuses
        const BASE: u32 = 0x0000_0008 | 0x0000_0200 | 0x0800_0000;
        const BREAKAWAY: u32 = 0x0100_0000;
        for flags in [BASE | BREAKAWAY, BASE] {
            let si = StartupInfoW {
                cb: std::mem::size_of::<StartupInfoW>() as u32,
                reserved: std::ptr::null_mut(),
                desktop: std::ptr::null_mut(),
                title: std::ptr::null_mut(),
                x: 0,
                y: 0,
                x_size: 0,
                y_size: 0,
                x_chars: 0,
                y_chars: 0,
                fill: 0,
                flags: 0,
                show: 0,
                reserved2: 0,
                reserved2_ptr: std::ptr::null_mut(),
                std_in: std::ptr::null_mut(),
                std_out: std::ptr::null_mut(),
                std_err: std::ptr::null_mut(),
            };
            let mut info =
                ProcessInformation { process: std::ptr::null_mut(), thread: std::ptr::null_mut(), pid: 0, tid: 0 };
            let mut buf = line.clone();
            // SAFETY: every pointer is to a live, correctly laid out value or null; the command line is a writable NUL-terminated buffer, as the
            // function requires; the two handles it returns are closed at once.
            let ok = unsafe {
                CreateProcessW(
                    app.as_ptr(),
                    buf.as_mut_ptr(),
                    std::ptr::null(),
                    std::ptr::null(),
                    0,
                    flags,
                    std::ptr::null(),
                    std::ptr::null(),
                    &si,
                    &mut info,
                )
            };
            if ok != 0 {
                // SAFETY: handles returned by the successful call above.
                unsafe {
                    CloseHandle(info.process);
                    CloseHandle(info.thread);
                }
                return Ok(());
            }
        }
        Err(std::io::Error::last_os_error())
    }
    #[cfg(not(windows))]
    {
        use std::os::unix::process::CommandExt;
        std::process::Command::new(program)
            .args(args)
            .stdin(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .process_group(0)
            .spawn()
            .map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn the_machine_says_how_much_memory_it_has_and_how_much_is_free() {
        if cfg!(any(windows, target_os = "linux")) {
            let (total, free) = super::system_memory().expect("the system gives it");
            assert!(total > 1_000_000_000 && free <= total && free > 0, "{total} {free}");
        }
    }

    #[test]
    fn this_process_has_a_peak_memory() {
        if cfg!(any(windows, target_os = "linux")) {
            let (ws, commit) = super::peak_memory().expect("the system gives it");
            assert!(ws > 1_000_000 && commit > 0, "{ws} {commit}");
        }
    }

    #[test]
    fn the_request_is_accepted_on_windows_and_ignored_elsewhere() {
        // a refusal on an old Windows is not a failure of the program; the call must simply not crash
        let accepted = super::opt_out_of_power_throttling();
        if !cfg!(windows) {
            assert!(!accepted);
        }
    }
}
