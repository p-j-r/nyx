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
