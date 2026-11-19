//! Runs a program and reports its wall time and its peak memory (productivity plan P0): the exact peak working set and peak commit that
//! Windows keeps for the child process, read from the child's handle after it has exited, so nothing is sampled and nothing is missed.
//! Run the built program (`target/release/examples/<name>.exe`), not `cargo run`: under cargo the numbers would be cargo's, not the solver's.
//!
//!   cargo run --release -p nyx-engine --example peakrun -- <program> [arguments...]
//!
//! The report goes to stderr as one line (`peakrun: exit 0, 12.34 s wall, peak working set 1.23 GB, peak commit 1.50 GB, cpu 40.1 s user + 0.9 s kernel`), the child's
//! own output is passed through unchanged. Windows only (the memory counters come from `K32GetProcessMemoryInfo`); elsewhere only the time is reported.

#[cfg(windows)]
mod mem {
    use std::ffi::c_void;
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
        fn K32GetProcessMemoryInfo(process: *mut c_void, counters: *mut Counters, cb: u32) -> i32;
    }
    #[repr(C)]
    #[derive(Default, Clone, Copy)]
    struct FileTime {
        low: u32,
        high: u32,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetProcessTimes(
            process: *mut c_void,
            creation: *mut FileTime,
            exit: *mut FileTime,
            kernel: *mut FileTime,
            user: *mut FileTime,
        ) -> i32;
    }
    /// (user, kernel) CPU seconds of the process behind `handle`.
    pub fn cpu(handle: *mut c_void) -> Option<(f64, f64)> {
        let mut t = [FileTime::default(); 4];
        // SAFETY: as in `peak`; the four out-parameters are valid FILETIME-sized structs.
        let ok = unsafe { GetProcessTimes(handle, &mut t[0], &mut t[1], &mut t[2], &mut t[3]) };
        let secs = |f: FileTime| (((f.high as u64) << 32) | f.low as u64) as f64 * 1e-7;
        (ok != 0).then(|| (secs(t[3]), secs(t[2])))
    }
    /// (peak working set, peak commit) in bytes of the process behind `handle`.
    pub fn peak(handle: *mut c_void) -> Option<(usize, usize)> {
        let mut c = Counters { cb: std::mem::size_of::<Counters>() as u32, ..Counters::default() };
        // SAFETY: `handle` is the live handle of a child process we own, `c` is a correctly sized and initialised out-parameter.
        let ok = unsafe { K32GetProcessMemoryInfo(handle, &mut c, c.cb) };
        (ok != 0).then_some((c.peak_working_set, c.peak_pagefile_usage))
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(program) = args.first() else {
        eprintln!("usage: peakrun <program> [arguments...]");
        std::process::exit(2);
    };
    let t = std::time::Instant::now();
    let mut child = match std::process::Command::new(program).args(&args[1..]).spawn() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("peakrun: cannot start {program}: {e}");
            std::process::exit(2);
        }
    };
    let status = child.wait().expect("waiting for the child");
    let wall = t.elapsed().as_secs_f64();
    #[cfg(windows)]
    let memory = {
        use std::os::windows::io::AsRawHandle;
        mem::peak(child.as_raw_handle() as *mut _)
            .map(|(ws, commit)| {
                format!(", peak working set {:.3} GB, peak commit {:.3} GB", ws as f64 / 1e9, commit as f64 / 1e9)
            })
            .unwrap_or_else(|| ", peak memory unavailable".into())
            + &mem::cpu(child.as_raw_handle() as *mut _)
                .map(|(u, k)| format!(", cpu {u:.1} s user + {k:.1} s kernel"))
                .unwrap_or_default()
    };
    #[cfg(not(windows))]
    let memory = String::from(", peak memory only measured on Windows");
    eprintln!("peakrun: exit {}, {wall:.2} s wall{memory}", status.code().unwrap_or(-1));
    std::process::exit(status.code().unwrap_or(1));
}
