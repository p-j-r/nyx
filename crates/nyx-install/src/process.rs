// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Detects running `nyx` binaries inside an install root and cooperates with ADR 95 shared workspace queues.

use anyhow::Result;
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunningProc {
    pub pid: u32,
    pub exe_path: PathBuf,
    pub name: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct QueueStatus {
    pub active_runner: bool,
    pub open_items: usize,
    pub paused: bool,
}

pub trait ProcessInspector {
    fn running_in_root(&self, root: &Path) -> Vec<RunningProc>;
    fn terminate_pid(&self, pid: u32) -> Result<()>;
    fn shared_queue_status(&self, workspace: &Path) -> QueueStatus;
    fn pause_shared_queue(&self, workspace: &Path) -> Result<bool>;
    fn resume_shared_queue(&self, workspace: &Path) -> Result<()>;
}

#[derive(Clone, Copy, Debug, Default)]
pub struct OsProcessInspector;

impl ProcessInspector for OsProcessInspector {
    fn running_in_root(&self, root: &Path) -> Vec<RunningProc> {
        os_running_in_root(root)
    }

    fn terminate_pid(&self, pid: u32) -> Result<()> {
        os_terminate_pid(pid)
    }

    fn shared_queue_status(&self, workspace: &Path) -> QueueStatus {
        inspect_shared_queue(workspace)
    }

    fn pause_shared_queue(&self, workspace: &Path) -> Result<bool> {
        let q = workspace.join("queue");
        if !q.is_dir() {
            return Ok(false);
        }
        let pause = q.join("pause");
        if pause.exists() {
            return Ok(false);
        }
        std::fs::write(&pause, b"")?;
        Ok(true)
    }

    fn resume_shared_queue(&self, workspace: &Path) -> Result<()> {
        let pause = workspace.join("queue").join("pause");
        if pause.is_file() {
            let _ = std::fs::remove_file(pause);
        }
        Ok(())
    }
}

pub fn inspect_shared_queue(workspace: &Path) -> QueueStatus {
    let q = workspace.join("queue");
    if !q.is_dir() {
        return QueueStatus::default();
    }
    let active_runner = std::fs::metadata(q.join("runner.lock"))
        .ok()
        .and_then(|m| m.modified().ok())
        .and_then(|t| SystemTime::now().duration_since(t).ok())
        .is_some_and(|age| age < Duration::from_secs(6));
    let paused = q.join("pause").is_file();
    let mut open_items = 0usize;
    if let Ok(rd) = std::fs::read_dir(q.join("items")) {
        for e in rd.filter_map(Result::ok) {
            let p = e.path();
            let name = e.file_name();
            if name.to_string_lossy().starts_with('.') || !p.is_dir() {
                continue;
            }
            if !p.join("request.json").is_file() {
                continue;
            }
            let finished = std::fs::read_to_string(p.join("state.json"))
                .ok()
                .and_then(|t| serde_json::from_str::<Value>(&t).ok())
                .and_then(|v| v["state"].as_str().map(str::to_string))
                .is_some_and(|s| matches!(s.as_str(), "done" | "failed" | "cancelled"));
            if !finished {
                open_items += 1;
            }
        }
    }
    QueueStatus {
        active_runner,
        open_items,
        paused,
    }
}

/// Removes `%TEMP%/nyx-*` dirs, `<workspace>/solve_history.jsonl`, and finished idle queue items (`--purge-cache`).
/// Never deletes workflows, specs, or `<workspace>/workspace.json`.
pub fn purge_cache_and_history(workspace: Option<&Path>, temp_dir: &Path) -> usize {
    let mut removed = 0usize;
    if let Ok(rd) = std::fs::read_dir(temp_dir) {
        for e in rd.filter_map(Result::ok) {
            let name = e.file_name().to_string_lossy().into_owned();
            if name.starts_with("nyx-")
                && e.path().is_dir()
                && std::fs::remove_dir_all(e.path()).is_ok()
            {
                removed += 1;
            }
        }
    }
    if let Some(ws) = workspace {
        let hist = ws.join("solve_history.jsonl");
        if hist.is_file() && std::fs::remove_file(&hist).is_ok() {
            removed += 1;
        }
        let items_dir = ws.join("queue").join("items");
        if let Ok(rd) = std::fs::read_dir(&items_dir) {
            for e in rd.filter_map(Result::ok) {
                let p = e.path();
                if !p.is_dir() {
                    continue;
                }
                let finished = std::fs::read_to_string(p.join("state.json"))
                    .ok()
                    .and_then(|t| serde_json::from_str::<Value>(&t).ok())
                    .and_then(|v| v["state"].as_str().map(str::to_string))
                    .is_some_and(|s| matches!(s.as_str(), "done" | "failed" | "cancelled"));
                if finished && std::fs::remove_dir_all(&p).is_ok() {
                    removed += 1;
                }
            }
        }
    }
    removed
}

fn path_inside(root: &Path, cand: &Path) -> bool {
    let r = crate::safety::normalize_lexical(root);
    let c = crate::safety::normalize_lexical(cand);
    #[cfg(windows)]
    {
        let rs = r.to_string_lossy().to_ascii_lowercase();
        let cs = c.to_string_lossy().to_ascii_lowercase();
        cs == rs || cs.starts_with(&format!("{rs}\\"))
    }
    #[cfg(not(windows))]
    {
        c == r || c.starts_with(&r)
    }
}

#[cfg(windows)]
fn os_running_in_root(root: &Path) -> Vec<RunningProc> {
    use std::ffi::c_void;
    const TH32CS_SNAPPROCESS: u32 = 0x0000_0002;
    const PROCESS_QUERY_LIMITED_INFORMATION: u32 = 0x1000;
    const INVALID_HANDLE_VALUE: *mut c_void = (-1isize) as *mut c_void;

    #[repr(C)]
    struct ProcessEntry32W {
        dw_size: u32,
        cnt_usage: u32,
        th32_process_id: u32,
        th32_default_heap_id: usize,
        th32_module_id: u32,
        cnt_threads: u32,
        th32_parent_process_id: u32,
        pc_pri_class_base: i32,
        dw_flags: u32,
        sz_exe_file: [u16; 260],
    }

    #[link(name = "kernel32")]
    extern "system" {
        fn CreateToolhelp32Snapshot(dw_flags: u32, th32_process_id: u32) -> *mut c_void;
        fn Process32FirstW(snapshot: *mut c_void, lppe: *mut ProcessEntry32W) -> i32;
        fn Process32NextW(snapshot: *mut c_void, lppe: *mut ProcessEntry32W) -> i32;
        fn OpenProcess(access: u32, inherit: i32, pid: u32) -> *mut c_void;
        fn QueryFullProcessImageNameW(
            proc_h: *mut c_void,
            flags: u32,
            buf: *mut u16,
            size: *mut u32,
        ) -> i32;
        fn CloseHandle(h: *mut c_void) -> i32;
    }

    let me = std::process::id();
    let mut out = Vec::new();
    unsafe {
        let snap = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
        if snap.is_null() || snap == INVALID_HANDLE_VALUE {
            return out;
        }
        let mut pe: ProcessEntry32W = std::mem::zeroed();
        pe.dw_size = std::mem::size_of::<ProcessEntry32W>() as u32;
        let mut ok = Process32FirstW(snap, &mut pe) != 0;
        while ok {
            let pid = pe.th32_process_id;
            if pid > 0 && pid != me {
                let h = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, 0, pid);
                if !h.is_null() {
                    let mut buf = vec![0u16; 1024];
                    let mut len = buf.len() as u32;
                    if QueryFullProcessImageNameW(h, 0, buf.as_mut_ptr(), &mut len) != 0 && len > 0
                    {
                        let exe = PathBuf::from(String::from_utf16_lossy(&buf[..len as usize]));
                        if path_inside(root, &exe) {
                            let name = exe
                                .file_name()
                                .and_then(|n| n.to_str())
                                .unwrap_or("nyx")
                                .to_string();
                            out.push(RunningProc {
                                pid,
                                exe_path: exe,
                                name,
                            });
                        }
                    }
                    CloseHandle(h);
                }
            }
            ok = Process32NextW(snap, &mut pe) != 0;
        }
        CloseHandle(snap);
    }
    out
}

#[cfg(windows)]
fn os_terminate_pid(pid: u32) -> Result<()> {
    use std::ffi::c_void;
    const PROCESS_TERMINATE: u32 = 0x0001;
    const SYNCHRONIZE: u32 = 0x0010_0000;
    #[link(name = "kernel32")]
    extern "system" {
        fn OpenProcess(access: u32, inherit: i32, pid: u32) -> *mut c_void;
        fn TerminateProcess(proc_h: *mut c_void, exit_code: u32) -> i32;
        fn WaitForSingleObject(h: *mut c_void, ms: u32) -> u32;
        fn CloseHandle(h: *mut c_void) -> i32;
    }
    unsafe {
        let h = OpenProcess(PROCESS_TERMINATE | SYNCHRONIZE, 0, pid);
        if h.is_null() {
            return Ok(());
        }
        let ok = TerminateProcess(h, 1);
        if ok != 0 {
            WaitForSingleObject(h, 5000);
        }
        CloseHandle(h);
    }
    Ok(())
}

#[cfg(target_os = "linux")]
fn os_running_in_root(root: &Path) -> Vec<RunningProc> {
    let me = std::process::id();
    let mut out = Vec::new();
    let Ok(rd) = std::fs::read_dir("/proc") else {
        return out;
    };
    for entry in rd.filter_map(Result::ok) {
        let Ok(pid) = entry.file_name().to_string_lossy().parse::<u32>() else {
            continue;
        };
        if pid == me {
            continue;
        }
        if let Ok(exe) = std::fs::read_link(entry.path().join("exe")) {
            let cleaned = PathBuf::from(exe.to_string_lossy().trim_end_matches(" (deleted)"));
            if path_inside(root, &cleaned) {
                let name = cleaned
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("nyx")
                    .to_string();
                out.push(RunningProc {
                    pid,
                    exe_path: cleaned,
                    name,
                });
            }
        }
    }
    out
}

#[cfg(target_os = "macos")]
fn os_running_in_root(root: &Path) -> Vec<RunningProc> {
    let me = std::process::id();
    let mut out = Vec::new();
    let Ok(ps) = std::process::Command::new("ps")
        .args(["-Ao", "pid=,comm="])
        .output()
    else {
        return out;
    };
    let text = String::from_utf8_lossy(&ps.stdout);
    for line in text.lines() {
        let trimmed = line.trim();
        let Some((pid_s, comm)) = trimmed.split_once(char::is_whitespace) else {
            continue;
        };
        let Ok(pid) = pid_s.trim().parse::<u32>() else {
            continue;
        };
        if pid == me {
            continue;
        }
        let exe = PathBuf::from(comm.trim());
        if path_inside(root, &exe) {
            let name = exe
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("nyx")
                .to_string();
            out.push(RunningProc {
                pid,
                exe_path: exe,
                name,
            });
        }
    }
    out
}

#[cfg(not(any(windows, target_os = "linux", target_os = "macos")))]
fn os_running_in_root(_root: &Path) -> Vec<RunningProc> {
    vec![]
}

#[cfg(unix)]
fn os_terminate_pid(pid: u32) -> Result<()> {
    unsafe {
        if libc::kill(pid as libc::pid_t, libc::SIGTERM) != 0 {
            return Ok(());
        }
    }
    for _ in 0..20 {
        std::thread::sleep(Duration::from_millis(50));
        unsafe {
            if libc::kill(pid as libc::pid_t, 0) != 0 {
                return Ok(());
            }
        }
    }
    unsafe {
        let _ = libc::kill(pid as libc::pid_t, libc::SIGKILL);
    }
    Ok(())
}

#[cfg(not(any(windows, unix)))]
fn os_terminate_pid(_pid: u32) -> Result<()> {
    Ok(())
}
