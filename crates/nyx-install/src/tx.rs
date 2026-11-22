// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 7-phase atomic transaction engine with append-only LIFO rollback journal (`install-journal.jsonl`),
//! process/queue quiesce, BLAKE3 staging verification, and safe uninstall.

use crate::bundle::{exe_filename, IntegrityVerifier, PayloadSource};
use crate::connect_hook::{self, ConnectPolicy};
use crate::manifest::{
    discover_manifest_from_root, file_abs_path, new_install_id, utc_now_iso, DesktopEntry,
    FileAssociation, FileBase, FileRole, InstallLayout, InstallManifest, InstallScope,
    ManifestFile, ManifestSymlink, PathEdit, TargetOs, MANIFEST_FORMAT, MANIFEST_SCHEMA,
};
use crate::platform::windows::{is_inside_dir, prepare_self_delete_sync_reaper};
use crate::platform::{HostPaths, LayoutPaths, PlatformBackend, PlatformOpts};
use crate::process::{purge_cache_and_history, ProcessInspector};
use crate::report::{DoctorCheck, HostOutcome, InstallReport};
use crate::safety::{
    file_blake3_hex, safe_remove_manifest_files, validate_install_root, validate_rel_path,
};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Structured error carrying an installer CLI exit code (`0` ok, `1` error, `2` usage, `3` elevation, `4` busy, `5` integrity).
#[derive(Debug)]
pub struct ExitError {
    pub code: i32,
    pub msg: String,
}

impl ExitError {
    pub fn new(code: i32, msg: impl Into<String>) -> Self {
        Self {
            code,
            msg: msg.into(),
        }
    }
}

impl std::fmt::Display for ExitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl std::error::Error for ExitError {}

pub fn exit_code_of(err: &anyhow::Error) -> i32 {
    if let Some(e) = err.downcast_ref::<ExitError>() {
        return e.code;
    }
    let s = err.to_string().to_ascii_lowercase();
    if s.contains("administrator") || s.contains("root privileges") || s.contains("elevation") {
        3
    } else if s.contains("currently running") || s.contains("shared solve queue") {
        4
    } else if s.contains("blake3")
        || s.contains("corrupted")
        || s.contains("tamper")
        || s.contains("trailer")
    {
        5
    } else {
        1
    }
}

#[derive(Clone, Debug)]
pub struct TxOptions {
    pub scope: InstallScope,
    pub prefix: Option<PathBuf>,
    pub add_path: bool,
    pub desktop: bool,
    pub desktop_explicit: bool,
    pub desktop_shortcut: bool,
    pub file_assoc: bool,
    pub connect_policy: ConnectPolicy,
    pub dry_run: bool,
    pub force: bool,
    pub allow_root_user: bool,
    pub purge_cache: bool,
    pub run_mcp_smoke: bool,
    pub uninstaller_stub: Option<PathBuf>,
    pub testing_reg_root: Option<String>,
}

impl TxOptions {
    pub fn default_user() -> Self {
        Self {
            scope: InstallScope::User,
            prefix: None,
            add_path: true,
            desktop: true,
            desktop_explicit: false,
            desktop_shortcut: false,
            file_assoc: true,
            connect_policy: ConnectPolicy::Detected,
            dry_run: false,
            force: false,
            allow_root_user: false,
            purge_cache: false,
            run_mcp_smoke: true,
            uninstaller_stub: None,
            testing_reg_root: None,
        }
    }

    pub fn effective_scope(&self) -> InstallScope {
        if self.prefix.is_some() && self.scope == InstallScope::User {
            InstallScope::Custom
        } else {
            self.scope
        }
    }

    pub fn to_platform_opts(&self, os: TargetOs, layout: InstallLayout) -> PlatformOpts {
        PlatformOpts {
            os,
            scope: self.effective_scope(),
            layout,
            add_path: self.add_path,
            desktop: self.desktop,
            desktop_explicit: self.desktop_explicit,
            desktop_shortcut: self.desktop_shortcut,
            file_assoc: self.file_assoc,
            force: self.force,
            allow_root_user: self.allow_root_user,
            dry_run: self.dry_run,
            testing_reg_root: self.testing_reg_root.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct TxProgress {
    pub phase: u8,
    pub total_phases: u8,
    pub fraction: f32,
    pub message: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum JournalOp {
    TempDir { path: String },
    CreatedFile { path: String },
    ReplacedFile { dest: String, backup: String },
    MovedAsideFile { orig: String, renamed: String },
    CreatedSymlinks { links: Vec<ManifestSymlink> },
    EditedPath { edits: Vec<PathEdit> },
    CreatedDesktop { entries: Vec<DesktopEntry> },
    EditedFileAssoc { assocs: Vec<FileAssociation> },
    ConnectedHosts { mcp_exe: String, ids: Vec<String> },
}

struct Journal {
    path: PathBuf,
    file: File,
    ops: Vec<JournalOp>,
}

impl Journal {
    fn open(path: &Path) -> Result<Self> {
        if let Some(p) = path.parent() {
            std::fs::create_dir_all(p)?;
        }
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            file,
            ops: Vec::new(),
        })
    }

    fn record(&mut self, op: JournalOp) -> Result<()> {
        let mut line = serde_json::to_string(&op)?;
        line.push('\n');
        self.file.write_all(line.as_bytes())?;
        self.file.flush()?;
        let _ = self.file.sync_all();
        self.ops.push(op);
        Ok(())
    }

    fn commit_and_remove(self) {
        let p = self.path.clone();
        drop(self);
        let _ = std::fs::remove_file(p);
    }
}

struct AdvisoryLock {
    path: PathBuf,
}

impl AdvisoryLock {
    fn acquire(lock_dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(lock_dir)?;
        let path = lock_dir.join(".nyx-install.lock");
        if path.is_file() {
            if let Ok(txt) = std::fs::read_to_string(&path) {
                if let Some(pid) = txt
                    .lines()
                    .next()
                    .and_then(|l| l.trim().parse::<u32>().ok())
                {
                    if pid != std::process::id() && pid_is_alive(pid) {
                        bail!(ExitError::new(
                            4,
                            format!(
                                "another nyx installer process (pid {pid}) holds {}",
                                path.display()
                            )
                        ));
                    }
                }
            }
            let _ = std::fs::remove_file(&path);
        }
        let mut f = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)?;
        let _ = writeln!(f, "{}\n{}", std::process::id(), utc_now_iso());
        Ok(Self { path })
    }
}

impl Drop for AdvisoryLock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

fn pid_is_alive(pid: u32) -> bool {
    #[cfg(unix)]
    {
        unsafe { libc::kill(pid as libc::pid_t, 0) == 0 }
    }
    #[cfg(windows)]
    {
        use std::ffi::c_void;
        #[link(name = "kernel32")]
        extern "system" {
            fn OpenProcess(access: u32, inherit: i32, pid: u32) -> *mut c_void;
            fn GetExitCodeProcess(h: *mut c_void, code: *mut u32) -> i32;
            fn CloseHandle(h: *mut c_void) -> i32;
        }
        unsafe {
            let h = OpenProcess(0x1000 /* PROCESS_QUERY_LIMITED_INFORMATION */, 0, pid);
            if h.is_null() {
                return false;
            }
            let mut code: u32 = 0;
            let ok = GetExitCodeProcess(h, &mut code) != 0 && code == 259 /* STILL_ACTIVE */;
            CloseHandle(h);
            ok
        }
    }
    #[cfg(not(any(unix, windows)))]
    {
        let _ = pid;
        false
    }
}

fn replay_rollback_ops(
    ops: &[JournalOp],
    backend: &dyn PlatformBackend,
    plat_opts: &PlatformOpts,
    paths: &HostPaths,
) {
    let connect_env = paths.connect_env(backend.os());
    for op in ops.iter().rev() {
        match op {
            JournalOp::ConnectedHosts { mcp_exe, ids } => {
                if let Some(env) = &connect_env {
                    for id in ids {
                        let _ = nyx_connect::disconnect_matching(env, id, Path::new(mcp_exe));
                    }
                }
            }
            JournalOp::EditedFileAssoc { assocs } => {
                let _ = backend.remove_file_assoc(assocs, plat_opts);
            }
            JournalOp::CreatedDesktop { entries } => {
                let _ = backend.remove_desktop(entries, plat_opts);
            }
            JournalOp::EditedPath { edits } => {
                let _ = backend.remove_path(edits, plat_opts);
            }
            JournalOp::CreatedSymlinks { links } => {
                let _ = backend.remove_symlinks(links, false);
            }
            JournalOp::MovedAsideFile { orig, renamed } => {
                let r = Path::new(renamed);
                let o = Path::new(orig);
                if r.exists() {
                    if o.exists() {
                        let _ = std::fs::remove_file(o);
                    }
                    let _ = std::fs::rename(r, o);
                }
            }
            JournalOp::ReplacedFile { dest, backup } => {
                let d = Path::new(dest);
                let b = Path::new(backup);
                if b.exists() {
                    if let Some(p) = d.parent() {
                        let _ = std::fs::create_dir_all(p);
                    }
                    if d.exists() {
                        let _ = std::fs::remove_file(d);
                    }
                    let _ = std::fs::rename(b, d);
                }
            }
            JournalOp::CreatedFile { path } => {
                let p = Path::new(path);
                let _ = std::fs::remove_file(p);
                if let Some(parent) = p.parent() {
                    let _ = std::fs::remove_dir(parent);
                }
            }
            JournalOp::TempDir { path } => {
                let _ = std::fs::remove_dir_all(path);
            }
        }
    }
}

fn recover_stale_journals(
    journal_dir: &Path,
    backend: &dyn PlatformBackend,
    plat_opts: &PlatformOpts,
    paths: &HostPaths,
) -> Vec<String> {
    let mut notes = Vec::new();
    let Ok(rd) = std::fs::read_dir(journal_dir) else {
        return notes;
    };
    for e in rd.filter_map(Result::ok) {
        let p = e.path();
        let name = e.file_name().to_string_lossy().into_owned();
        if name.starts_with(".nyx-install-") && name.ends_with(".journal.jsonl") && p.is_file() {
            if let Ok(f) = File::open(&p) {
                let mut ops = Vec::new();
                for line in BufReader::new(f).lines().map_while(Result::ok) {
                    if let Ok(op) = serde_json::from_str::<JournalOp>(&line) {
                        ops.push(op);
                    }
                }
                if !ops.is_empty() {
                    replay_rollback_ops(&ops, backend, plat_opts, paths);
                    notes.push(format!(
                        "recovered and rolled back incomplete previous transaction from {}",
                        p.display()
                    ));
                }
            }
            let _ = std::fs::remove_file(&p);
        }
    }
    notes
}

fn check_free_disk_space(target_dir: &Path, needed_bytes: u64) -> Result<()> {
    let mut check_path = target_dir;
    while !check_path.exists() {
        match check_path.parent() {
            Some(p) => check_path = p,
            None => break,
        }
    }
    #[cfg(unix)]
    {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;
        if let Ok(c) = CString::new(check_path.as_os_str().as_bytes()) {
            let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
            if unsafe { libc::statvfs(c.as_ptr(), &mut st) } == 0 {
                let avail = (st.f_bavail as u64).saturating_mul(st.f_frsize as u64);
                if avail < needed_bytes {
                    bail!(
                        "insufficient disk space on {}: {} MiB free, need at least {} MiB",
                        check_path.display(),
                        avail >> 20,
                        needed_bytes >> 20
                    );
                }
            }
        }
    }
    #[cfg(windows)]
    {
        use std::os::windows::ffi::OsStrExt;
        #[link(name = "kernel32")]
        extern "system" {
            fn GetDiskFreeSpaceExW(
                dir: *const u16,
                free_avail: *mut u64,
                total_bytes: *mut u64,
                total_free: *mut u64,
            ) -> i32;
        }
        let w: Vec<u16> = check_path
            .as_os_str()
            .encode_wide()
            .chain(Some(0))
            .collect();
        let mut avail: u64 = 0;
        let mut total: u64 = 0;
        let mut free_total: u64 = 0;
        if unsafe { GetDiskFreeSpaceExW(w.as_ptr(), &mut avail, &mut total, &mut free_total) } != 0
            && avail < needed_bytes
        {
            bail!(
                "insufficient disk space on {}: {} MiB free, need at least {} MiB",
                check_path.display(),
                avail >> 20,
                needed_bytes >> 20
            );
        }
    }
    Ok(())
}

fn check_processes_and_queue(
    roots: &[&Path],
    inspector: &dyn ProcessInspector,
    workspace: Option<&Path>,
    force: bool,
    dry_run: bool,
) -> Result<(bool, Vec<String>)> {
    let mut warnings = Vec::new();
    for root in roots {
        let procs = inspector.running_in_root(root);
        if !procs.is_empty() {
            if !force {
                let list: Vec<String> = procs
                    .iter()
                    .map(|p| format!("{} (pid {})", p.name, p.pid))
                    .collect();
                bail!(ExitError::new(
                    4,
                    format!(
                        "running nyx process(es) in {}: {}. Close them first or re-run with --force.",
                        root.display(),
                        list.join(", ")
                    )
                ));
            }
            if !dry_run {
                for p in procs {
                    warnings.push(format!(
                        "terminated running process {} (pid {})",
                        p.name, p.pid
                    ));
                    inspector.terminate_pid(p.pid)?;
                }
            }
        }
    }

    let mut paused_here = false;
    if let Some(ws) = workspace {
        let qs = inspector.shared_queue_status(ws);
        if qs.active_runner || qs.open_items > 0 {
            if !force {
                bail!(ExitError::new(
                    4,
                    format!(
                        "shared workspace queue at {} has {} open solve item(s){}. Wait for jobs or pass --force.",
                        ws.display(),
                        qs.open_items,
                        if qs.active_runner { " (runner active)" } else { "" }
                    )
                ));
            }
            if !dry_run {
                paused_here = inspector.pause_shared_queue(ws)?;
                if paused_here {
                    warnings.push(format!("paused shared queue in {}", ws.display()));
                }
            }
        }
    }
    Ok((paused_here, warnings))
}

fn journal_parent_dir(layout: &LayoutPaths) -> PathBuf {
    layout
        .manifest_path
        .parent()
        .filter(|_| layout.layout == InstallLayout::MacosApp)
        .map(Path::to_path_buf)
        .unwrap_or_else(|| {
            layout
                .install_root
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_else(|| layout.install_root.clone())
        })
}

fn dest_path_for(layout: &LayoutPaths, base: FileBase, rel_path: &str) -> PathBuf {
    let dir = match base {
        FileBase::Bin => &layout.bin_dir,
        FileBase::Resource => &layout.resource_dir,
        FileBase::Share => &layout.share_dir,
    };
    dir.join(rel_path)
}

fn read_manifest_if_exists(path: &Path) -> Option<InstallManifest> {
    let text = std::fs::read_to_string(path).ok()?;
    serde_json::from_str::<InstallManifest>(&text).ok()
}

pub fn locate_manifest(
    prefix_or_root: Option<&Path>,
    current_exe: Option<&Path>,
    backend: &dyn PlatformBackend,
    scope: InstallScope,
    paths: &HostPaths,
) -> Option<(LayoutPaths, InstallManifest)> {
    let mut candidates = Vec::new();
    if let Some(p) = prefix_or_root {
        candidates.push(p.to_path_buf());
        if backend.os() == TargetOs::Macos && p.extension().and_then(|e| e.to_str()) != Some("app")
        {
            candidates.push(p.join("Nyx.app"));
        }
    }
    if let Some(exe) = current_exe {
        if let Some(bin_dir) = exe.parent() {
            if let Some(root) = bin_dir.parent() {
                if root.file_name().and_then(|n| n.to_str()) == Some("Contents") {
                    if let Some(app) = root.parent() {
                        candidates.push(app.to_path_buf());
                    }
                } else {
                    candidates.push(root.to_path_buf());
                }
            }
        }
    }
    candidates.push(backend.default_install_root(scope, paths));

    for root in candidates {
        let lp = backend.resolve_layout(&root, scope, paths);
        if let Some(m) = read_manifest_if_exists(&lp.manifest_path) {
            return Some((lp, m));
        }
        if let Some(alt) =
            discover_manifest_from_root(backend.os(), scope, &root, paths.home.as_deref())
        {
            if let Some(m) = read_manifest_if_exists(&alt) {
                return Some((lp, m));
            }
        }
    }
    None
}

fn emit_progress(cb: &mut dyn FnMut(TxProgress), phase: u8, fraction: f32, msg: impl Into<String>) {
    cb(TxProgress {
        phase,
        total_phases: 7,
        fraction: fraction.clamp(0.0, 1.0),
        message: msg.into(),
    });
}

fn sweep_old_renamed_files(dir: &Path) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for e in rd.filter_map(Result::ok) {
        let p = e.path();
        let name = e.file_name().to_string_lossy().into_owned();
        if p.is_file() && name.contains(".nyx-old-") {
            let _ = std::fs::remove_file(p);
        } else if p.is_dir() && !name.starts_with('.') {
            sweep_old_renamed_files(&p);
        }
    }
}

/// Executes Phase 1-7 atomic `install`, `upgrade`, or `repair` transaction.
pub fn run_install_or_upgrade(
    action: &str,
    source: &PayloadSource,
    opts: &TxOptions,
    backend: &dyn PlatformBackend,
    inspector: &dyn ProcessInspector,
    paths: &HostPaths,
    on_progress: &mut dyn FnMut(TxProgress),
) -> Result<InstallReport> {
    backend.harden_process();
    let effective_scope = opts.effective_scope();
    backend.check_elevation(effective_scope, opts.allow_root_user)?;

    let install_root = opts
        .prefix
        .clone()
        .map(|p| {
            if backend.os() == TargetOs::Macos
                && p.extension().and_then(|e| e.to_str()) != Some("app")
            {
                p.join("Nyx.app")
            } else {
                p
            }
        })
        .unwrap_or_else(|| backend.default_install_root(effective_scope, paths));

    validate_install_root(
        &install_root,
        backend.os(),
        effective_scope,
        paths.home.as_deref(),
        paths.local_appdata.as_deref(),
        paths.program_files.as_deref(),
        opts.force,
    )?;

    let layout = backend.resolve_layout(&install_root, effective_scope, paths);
    let plat_opts = opts.to_platform_opts(backend.os(), layout.layout);
    let journal_parent = journal_parent_dir(&layout);

    emit_progress(
        on_progress,
        1,
        0.05,
        "Preflight checks and payload verification",
    );
    let verifier = IntegrityVerifier;
    let opened = source.metadata_summary(&verifier)?;
    let old_manifest = read_manifest_if_exists(&layout.manifest_path);

    let effective_action = match (action, old_manifest.as_ref()) {
        ("install", Some(_)) => "upgrade",
        _ => action,
    };

    let mut report = InstallReport::new(
        effective_action,
        effective_scope,
        layout.install_root.display().to_string(),
        opts.dry_run,
    );
    report.nyx_version = Some(opened.toc.nyx_version.clone());
    report.previous_version = old_manifest.as_ref().map(|m| m.nyx_version.clone());
    report.manifest_path = Some(layout.manifest_path.display().to_string());

    let required_bytes = ((opened.total_uncompressed as f64) * 2.2) as u64;
    check_free_disk_space(&install_root, required_bytes)?;

    if opts.dry_run {
        report.files_written = opened.toc.entries.len();
        report.path_updated = opts.add_path;
        report.desktop_updated = opts.desktop;
        report.file_assoc_updated = opts.file_assoc;
        let mcp_exe = layout.bin_dir.join(exe_filename("nyx-mcp", backend.os()));
        let old_mcp = old_manifest
            .as_ref()
            .map(|m| PathBuf::from(&m.bin_dir).join(exe_filename("nyx-mcp", backend.os())));
        let (_, host_outcomes, w, next) = connect_hook::apply_connections(
            paths.connect_env(backend.os()).as_ref(),
            effective_scope,
            &mcp_exe,
            old_mcp.as_deref(),
            old_manifest.as_ref().map_or(&[], |m| &m.ai_hosts),
            &opts.connect_policy,
            true,
        );
        report.ai_hosts = host_outcomes;
        report.warnings.extend(w);
        report.next_steps.extend(next);
        emit_progress(on_progress, 7, 1.0, "Dry run complete");
        return Ok(report);
    }

    let _lock = AdvisoryLock::acquire(&journal_parent)?;
    report.warnings.extend(recover_stale_journals(
        &journal_parent,
        backend,
        &plat_opts,
        paths,
    ));
    if layout.install_root.is_dir() {
        sweep_old_renamed_files(&layout.install_root);
    }

    let install_id = new_install_id();
    let journal_path = journal_parent.join(format!(".nyx-install-{install_id}.journal.jsonl"));
    let mut journal = Journal::open(&journal_path)?;

    let mut check_roots: Vec<&Path> = vec![layout.install_root.as_path()];
    let old_root_buf = old_manifest
        .as_ref()
        .map(|m| PathBuf::from(&m.install_root));
    if let Some(ref or) = old_root_buf {
        if or != &layout.install_root {
            check_roots.push(or.as_path());
        }
    }

    emit_progress(
        on_progress,
        2,
        0.15,
        "Checking active nyx processes and shared queue",
    );
    let (paused_queue, proc_warnings) = check_processes_and_queue(
        &check_roots,
        inspector,
        paths.documents_dir.as_deref(),
        opts.force,
        opts.dry_run,
    )?;
    report.warnings.extend(proc_warnings);

    let run_result = (|| -> Result<InstallManifest> {
        // Phase 3: Stage & Verify BLAKE3
        emit_progress(on_progress, 3, 0.25, "Staging and verifying files");
        let stage_dir = journal_parent.join(format!(".nyx-stage-{install_id}"));
        let backup_dir = journal_parent.join(format!(".nyx-backup-{install_id}"));
        std::fs::create_dir_all(&stage_dir)?;
        journal.record(JournalOp::TempDir {
            path: stage_dir.display().to_string(),
        })?;
        std::fs::create_dir_all(&backup_dir)?;
        journal.record(JournalOp::TempDir {
            path: backup_dir.display().to_string(),
        })?;

        let current_exe = std::env::current_exe().unwrap_or_default();
        let stub = opts
            .uninstaller_stub
            .as_deref()
            .unwrap_or(current_exe.as_path());
        let toc = source.extract_to_stage(
            &stage_dir,
            backend.os(),
            stub,
            &verifier,
            &mut |idx, total, rel| {
                let frac = 0.25 + 0.30 * (idx as f32 / total.max(1) as f32);
                emit_progress(on_progress, 3, frac, format!("Extracting {rel}"));
            },
        )?;

        // Phase 4: Atomic Swap / Promote
        emit_progress(on_progress, 4, 0.60, "Installing verified files");
        let mut manifest_files = Vec::with_capacity(toc.entries.len());
        let mut new_targets = BTreeSet::new();
        let mut bin_names = Vec::new();

        for (idx, entry) in toc.entries.iter().enumerate() {
            validate_rel_path(&entry.rel_path)?;
            let src_sub = match entry.base {
                FileBase::Bin => stage_dir.join("bin"),
                FileBase::Resource => stage_dir.join("resource"),
                FileBase::Share => stage_dir.join("share"),
            };
            let staged_file = src_sub.join(&entry.rel_path);
            let dst = dest_path_for(&layout, entry.base, &entry.rel_path);
            new_targets.insert(dst.clone());
            if matches!(entry.role, FileRole::Binary | FileRole::Uninstaller) {
                bin_names.push(entry.rel_path.clone());
            }
            if let Some(p) = dst.parent() {
                std::fs::create_dir_all(p)?;
            }

            if dst.exists() {
                let bak = backup_dir.join(format!(
                    "{idx:04}_{}",
                    entry.rel_path.replace(['/', '\\'], "_")
                ));
                if std::fs::rename(&dst, &bak).is_ok() {
                    journal.record(JournalOp::ReplacedFile {
                        dest: dst.display().to_string(),
                        backup: bak.display().to_string(),
                    })?;
                } else {
                    // Windows locked file fallback: rename aside before placing new file.
                    let aside = dst.with_extension(format!("nyx-old-{install_id}"));
                    std::fs::rename(&dst, &aside)
                        .with_context(|| format!("replacing locked file {}", dst.display()))?;
                    journal.record(JournalOp::MovedAsideFile {
                        orig: dst.display().to_string(),
                        renamed: aside.display().to_string(),
                    })?;
                }
            } else {
                journal.record(JournalOp::CreatedFile {
                    path: dst.display().to_string(),
                })?;
            }

            if std::fs::rename(&staged_file, &dst).is_err() {
                std::fs::copy(&staged_file, &dst)?;
                let _ = std::fs::remove_file(&staged_file);
            }
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(entry.mode));
            }
            let hash = if entry.blake3.is_empty() {
                file_blake3_hex(&dst)?
            } else {
                entry.blake3.clone()
            };
            manifest_files.push(ManifestFile {
                rel_path: entry.rel_path.clone(),
                base: entry.base,
                role: entry.role,
                size: entry.uncompressed_size,
                blake3: hash,
                mode: entry.mode,
            });
        }

        // Remove obsolete files from previous manifest (verified BLAKE3; backed up for rollback).
        if let Some(ref prev) = old_manifest {
            for old_f in &prev.files {
                let old_abs = file_abs_path(prev, old_f);
                if !new_targets.contains(&old_abs) && old_abs.is_file() {
                    let ok_hash = file_blake3_hex(&old_abs).is_ok_and(|h| h == old_f.blake3);
                    if ok_hash || opts.force {
                        let bak = backup_dir
                            .join(format!("old_{}", old_f.rel_path.replace(['/', '\\'], "_")));
                        if std::fs::rename(&old_abs, &bak).is_ok() {
                            journal.record(JournalOp::ReplacedFile {
                                dest: old_abs.display().to_string(),
                                backup: bak.display().to_string(),
                            })?;
                        }
                    }
                }
            }
        }
        report.files_written = manifest_files.len();

        // Phase 5: OS Integrations & nyx-connect
        emit_progress(
            on_progress,
            5,
            0.78,
            "Configuring PATH, shortcuts, file associations, and AI hosts",
        );
        let symlinks = backend.install_symlinks(
            &layout.bin_dir,
            layout.link_dir.as_deref(),
            &bin_names,
            opts.force,
            false,
        )?;
        if !symlinks.is_empty() {
            journal.record(JournalOp::CreatedSymlinks {
                links: symlinks.clone(),
            })?;
        }

        let old_bin = old_manifest.as_ref().map(|m| PathBuf::from(&m.bin_dir));
        let (path_edits, path_warn) =
            backend.apply_path(&layout.bin_dir, old_bin.as_deref(), &plat_opts, paths)?;
        report.warnings.extend(path_warn);
        if !path_edits.is_empty() {
            report.path_updated = true;
            journal.record(JournalOp::EditedPath {
                edits: path_edits.clone(),
            })?;
        }

        let now = utc_now_iso();
        let mut manifest = InstallManifest {
            format: MANIFEST_FORMAT,
            schema: MANIFEST_SCHEMA.to_string(),
            nyx_version: toc.nyx_version,
            installer_version: env!("CARGO_PKG_VERSION").to_string(),
            install_id: install_id.clone(),
            target_triple: toc.target_triple,
            os: backend.os(),
            scope: effective_scope,
            layout: layout.layout,
            install_root: layout.install_root.display().to_string(),
            bin_dir: layout.bin_dir.display().to_string(),
            resource_dir: layout.resource_dir.display().to_string(),
            manifest_path: layout.manifest_path.display().to_string(),
            installed_at: old_manifest
                .as_ref()
                .map_or_else(|| now.clone(), |m| m.installed_at.clone()),
            updated_at: now,
            previous_version: old_manifest.as_ref().map(|m| m.nyx_version.clone()),
            files: manifest_files,
            symlinks,
            path_edits,
            desktop_entries: vec![],
            file_associations: vec![],
            ai_hosts: vec![],
        };

        let (extra_files, bundle_warn) = backend.finalize_bundle(&manifest, &plat_opts, paths)?;
        for ef in extra_files {
            if !manifest
                .files
                .iter()
                .any(|f| f.base == ef.base && f.rel_path == ef.rel_path)
            {
                let abs = file_abs_path(&manifest, &ef);
                journal.record(JournalOp::CreatedFile {
                    path: abs.display().to_string(),
                })?;
                manifest.files.push(ef);
            }
        }
        report.files_written = manifest.files.len();
        report.warnings.extend(bundle_warn);

        let (desktop_entries, d_warn) = backend.install_desktop(&manifest, &plat_opts, paths)?;
        report.warnings.extend(d_warn);
        if !desktop_entries.is_empty() {
            report.desktop_updated = true;
            journal.record(JournalOp::CreatedDesktop {
                entries: desktop_entries.clone(),
            })?;
        }
        manifest.desktop_entries = desktop_entries;

        let (file_assocs, f_warn) = backend.install_file_assoc(&manifest, &plat_opts, paths)?;
        report.warnings.extend(f_warn);
        if !file_assocs.is_empty() {
            report.file_assoc_updated = true;
            journal.record(JournalOp::EditedFileAssoc {
                assocs: file_assocs.clone(),
            })?;
        }
        manifest.file_associations = file_assocs;

        let mcp_exe = layout.bin_dir.join(exe_filename("nyx-mcp", backend.os()));
        let old_mcp = old_manifest
            .as_ref()
            .map(|m| PathBuf::from(&m.bin_dir).join(exe_filename("nyx-mcp", backend.os())));
        let (ai_records, host_outcomes, c_warn, c_next) = connect_hook::apply_connections(
            paths.connect_env(backend.os()).as_ref(),
            effective_scope,
            &mcp_exe,
            old_mcp.as_deref(),
            old_manifest.as_ref().map_or(&[], |m| &m.ai_hosts),
            &opts.connect_policy,
            false,
        );
        let new_connected_ids: Vec<String> = ai_records
            .iter()
            .filter(|r| r.previous_state == "not_connected")
            .map(|r| r.host_id.clone())
            .collect();
        if !new_connected_ids.is_empty() {
            journal.record(JournalOp::ConnectedHosts {
                mcp_exe: nyx_connect::plain_path(mcp_exe.clone())
                    .display()
                    .to_string(),
                ids: new_connected_ids,
            })?;
        }
        manifest.ai_hosts = ai_records;
        report.ai_hosts = host_outcomes;
        report.warnings.extend(c_warn);
        report.next_steps.extend(c_next);

        // Phase 6: Write Manifest & MCP Smoke
        emit_progress(
            on_progress,
            6,
            0.92,
            "Writing install manifest and verifying installation",
        );
        if let Some(mp) = layout.manifest_path.parent() {
            std::fs::create_dir_all(mp)?;
        }
        let m_tmp = layout
            .manifest_path
            .with_extension(format!("tmp-{install_id}"));
        let m_json = serde_json::to_string_pretty(&manifest)? + "\n";
        std::fs::write(&m_tmp, m_json.as_bytes())?;
        if layout.manifest_path.exists() {
            let _ = std::fs::remove_file(&layout.manifest_path);
        }
        std::fs::rename(&m_tmp, &layout.manifest_path)?;

        if opts.run_mcp_smoke && mcp_exe.is_file() {
            match connect_hook::smoke_test_mcp(&mcp_exe) {
                Ok(info) => report.doctor_checks.push(DoctorCheck {
                    name: "mcp_handshake".into(),
                    tier: "required".into(),
                    status: "ok".into(),
                    detail: format!(
                        "{} {} ({} tools)",
                        info.name,
                        info.version,
                        info.tools.len()
                    ),
                }),
                Err(e) => {
                    report.doctor_checks.push(DoctorCheck {
                        name: "mcp_handshake".into(),
                        tier: "required".into(),
                        status: "warn".into(),
                        detail: e.clone(),
                    });
                    report
                        .warnings
                        .push(format!("MCP verification warning: {e}"));
                }
            }
        }

        let _ = std::fs::remove_dir_all(&stage_dir);
        let _ = std::fs::remove_dir_all(&backup_dir);
        Ok(manifest)
    })();

    if paused_queue {
        if let Some(ws) = paths.documents_dir.as_deref() {
            let _ = inspector.resume_shared_queue(ws);
        }
    }

    match run_result {
        Ok(_manifest) => {
            journal.commit_and_remove();
            if backend.os() != TargetOs::Windows
                && opts.add_path
                && effective_scope == InstallScope::User
            {
                report.next_steps.push(
                    "Restart your terminal or source your shell profile to use `nyx` on PATH."
                        .into(),
                );
            }
            emit_progress(on_progress, 7, 1.0, "Completed");
            Ok(report)
        }
        Err(err) => {
            emit_progress(on_progress, 7, 1.0, "Rolling back failed transaction");
            let ops = journal.ops.clone();
            replay_rollback_ops(&ops, backend, &plat_opts, paths);
            journal.commit_and_remove();
            Err(err)
        }
    }
}

/// Removes an installed `nyx` distribution cleanly without touching `Documents/nyx`.
pub fn run_uninstall(
    opts: &TxOptions,
    backend: &dyn PlatformBackend,
    inspector: &dyn ProcessInspector,
    paths: &HostPaths,
    on_progress: &mut dyn FnMut(TxProgress),
) -> Result<InstallReport> {
    backend.harden_process();
    let effective_scope = opts.effective_scope();
    backend.check_elevation(effective_scope, opts.allow_root_user)?;

    let current_exe = std::env::current_exe().ok();
    let Some((layout, manifest)) = locate_manifest(
        opts.prefix.as_deref(),
        current_exe.as_deref(),
        backend,
        effective_scope,
        paths,
    ) else {
        bail!("no nyx installation manifest found to uninstall (pass --prefix <INSTALL_DIR> if installed in a custom folder)");
    };

    let install_root = PathBuf::from(&manifest.install_root);
    validate_install_root(
        &install_root,
        backend.os(),
        manifest.scope,
        paths.home.as_deref(),
        paths.local_appdata.as_deref(),
        paths.program_files.as_deref(),
        true,
    )?;

    let plat_opts = opts.to_platform_opts(backend.os(), manifest.layout);
    let mut report = InstallReport::new(
        "uninstall",
        manifest.scope,
        manifest.install_root.clone(),
        opts.dry_run,
    );
    report.nyx_version = Some(manifest.nyx_version.clone());
    report.manifest_path = Some(manifest.manifest_path.clone());

    if opts.dry_run {
        report.files_removed = manifest.files.len();
        report.path_updated = !manifest.path_edits.is_empty();
        report.desktop_updated = !manifest.desktop_entries.is_empty();
        report.file_assoc_updated = !manifest.file_associations.is_empty();
        return Ok(report);
    }

    // Disconnect AI hosts FIRST so IDEs do not respawn `nyx-mcp` before we inspect processes!
    emit_progress(on_progress, 1, 0.15, "Disconnecting AI host MCP entries");
    let mcp_exe = PathBuf::from(&manifest.bin_dir).join(exe_filename("nyx-mcp", backend.os()));
    let (host_outcomes, h_warn) = connect_hook::remove_matching_connections(
        paths.connect_env(backend.os()).as_ref(),
        &mcp_exe,
        false,
    );
    report.ai_hosts = host_outcomes;
    report.warnings.extend(h_warn);

    emit_progress(
        on_progress,
        2,
        0.30,
        "Checking running nyx processes and shared queue",
    );
    let (paused_queue, proc_warnings) = check_processes_and_queue(
        &[install_root.as_path()],
        inspector,
        paths.documents_dir.as_deref(),
        opts.force,
        false,
    )?;
    report.warnings.extend(proc_warnings);

    emit_progress(
        on_progress,
        3,
        0.45,
        "Removing OS desktop, PATH, and file associations",
    );
    let _ = backend.prepare_uninstall_hook(&manifest, &plat_opts);
    report.desktop_updated = backend.remove_desktop(&manifest.desktop_entries, &plat_opts)?;
    report.file_assoc_updated =
        backend.remove_file_assoc(&manifest.file_associations, &plat_opts)?;
    report.path_updated = backend.remove_path(&manifest.path_edits, &plat_opts)?;
    let _ = backend.remove_symlinks(&manifest.symlinks, false)?;

    // Change CWD away from install_root BEFORE removing files (essential on Windows & clean Unix rmdir).
    let _ = std::env::set_current_dir(&paths.temp_dir);
    let running_uninstaller_inside = backend.os() == TargetOs::Windows
        && current_exe
            .as_deref()
            .is_some_and(|exe| is_inside_dir(exe, &install_root));

    emit_progress(on_progress, 4, 0.70, "Removing verified installed files");
    let summary = safe_remove_manifest_files(&manifest, opts.force, running_uninstaller_inside)?;
    report.files_removed = summary.removed_files;
    for skipped in summary.skipped_modified {
        report.warnings.push(format!(
            "preserved modified file {} (pass --force to remove)",
            skipped.display()
        ));
    }

    if paused_queue {
        if let Some(ws) = paths.documents_dir.as_deref() {
            let _ = inspector.resume_shared_queue(ws);
        }
    }

    if opts.purge_cache {
        emit_progress(
            on_progress,
            5,
            0.85,
            "Clearing temporary cache and solve history",
        );
        let n = purge_cache_and_history(paths.documents_dir.as_deref(), &paths.temp_dir);
        report
            .warnings
            .push(format!("purged {n} temporary cache/history item(s)"));
    }

    if let Some(ws) = &paths.documents_dir {
        report.next_steps.push(format!(
            "Your workflows and simulation results in {} were preserved.",
            ws.display()
        ));
    }

    if let Some(deferred) = summary.deferred_uninstaller {
        emit_progress(on_progress, 6, 0.95, "Scheduling final uninstaller cleanup");
        prepare_self_delete_sync_reaper(&deferred, &install_root, std::process::id())?;
    }

    let _ = &layout;
    emit_progress(on_progress, 7, 1.0, "Uninstall complete");
    Ok(report)
}

/// Reports current install status, versions, and optional MCP/binary checks without mutating anything.
pub fn run_status(
    opts: &TxOptions,
    backend: &dyn PlatformBackend,
    paths: &HostPaths,
) -> Result<InstallReport> {
    let scope = opts.effective_scope();
    let current_exe = std::env::current_exe().ok();
    let Some((_, manifest)) = locate_manifest(
        opts.prefix.as_deref(),
        current_exe.as_deref(),
        backend,
        scope,
        paths,
    ) else {
        let fallback_root = opts
            .prefix
            .clone()
            .unwrap_or_else(|| backend.default_install_root(scope, paths));
        let mut r = InstallReport::new("status", scope, fallback_root.display().to_string(), false);
        r.ok = false;
        r.errors
            .push("nyx is not installed at this location".into());
        return Ok(r);
    };

    let mut report = InstallReport::new(
        "status",
        manifest.scope,
        manifest.install_root.clone(),
        false,
    );
    report.manifest_path = Some(manifest.manifest_path.clone());
    report.nyx_version = Some(manifest.nyx_version.clone());
    report.previous_version = manifest.previous_version.clone();

    if let Some(env) = paths.connect_env(backend.os()) {
        let mcp_exe = PathBuf::from(&manifest.bin_dir).join(exe_filename("nyx-mcp", backend.os()));
        for h in nyx_connect::list(&env, &mcp_exe) {
            report.ai_hosts.push(HostOutcome {
                id: h.id.to_string(),
                action: "status".to_string(),
                status: connect_hook::status_label(&h.status).to_string(),
                config_path: Some(h.config.display().to_string()),
                backup_path: None,
                detail: None,
            });
        }
    }
    Ok(report)
}

/// Runs verification checks across installed files (BLAKE3), OS integration, MCP handshake, and GUI libs.
pub fn run_doctor(
    opts: &TxOptions,
    backend: &dyn PlatformBackend,
    paths: &HostPaths,
) -> Result<InstallReport> {
    let scope = opts.effective_scope();
    let current_exe = std::env::current_exe().ok();
    let found = locate_manifest(
        opts.prefix.as_deref(),
        current_exe.as_deref(),
        backend,
        scope,
        paths,
    );
    let root_str = found
        .as_ref()
        .map(|(_, m)| m.install_root.clone())
        .unwrap_or_else(|| {
            opts.prefix
                .clone()
                .unwrap_or_else(|| backend.default_install_root(scope, paths))
                .display()
                .to_string()
        });
    let mut report = InstallReport::new("doctor", scope, root_str, false);
    let manifest_ref = found.as_ref().map(|(_, m)| m);

    if let Some(m) = manifest_ref {
        report.manifest_path = Some(m.manifest_path.clone());
        report.nyx_version = Some(m.nyx_version.clone());
        let mut missing = 0usize;
        let mut modified = 0usize;
        for f in &m.files {
            let abs = file_abs_path(m, f);
            if !abs.is_file() {
                missing += 1;
            } else if file_blake3_hex(&abs).is_ok_and(|h| h != f.blake3) {
                modified += 1;
            }
        }
        let ok = missing == 0 && modified == 0;
        if !ok {
            report.ok = false;
        }
        report.doctor_checks.push(DoctorCheck {
            name: "manifest_integrity".into(),
            tier: "required".into(),
            status: if ok { "ok" } else { "fail" }.into(),
            detail: format!(
                "{} files verified ({missing} missing, {modified} modified)",
                m.files.len()
            ),
        });

        let mcp_exe = PathBuf::from(&m.bin_dir).join(exe_filename("nyx-mcp", backend.os()));
        if opts.run_mcp_smoke && mcp_exe.is_file() {
            match connect_hook::smoke_test_mcp(&mcp_exe) {
                Ok(info) => report.doctor_checks.push(DoctorCheck {
                    name: "mcp_handshake".into(),
                    tier: "required".into(),
                    status: "ok".into(),
                    detail: format!(
                        "{} {} ({} tools)",
                        info.name,
                        info.version,
                        info.tools.len()
                    ),
                }),
                Err(e) => {
                    report.ok = false;
                    report.doctor_checks.push(DoctorCheck {
                        name: "mcp_handshake".into(),
                        tier: "required".into(),
                        status: "fail".into(),
                        detail: e,
                    });
                }
            }
        }
    } else {
        report.ok = false;
        report.doctor_checks.push(DoctorCheck {
            name: "manifest_present".into(),
            tier: "required".into(),
            status: "fail".into(),
            detail: "install-manifest.json not found".into(),
        });
    }

    let plat_checks = backend.doctor_checks(manifest_ref, paths);
    for c in &plat_checks {
        if c.tier == "required" && c.status == "fail" {
            report.ok = false;
        }
    }
    report.doctor_checks.extend(plat_checks);

    // Check which hosts currently see nyx.
    if let (Some(env), Some(m)) = (paths.connect_env(backend.os()), manifest_ref) {
        let mcp_exe = PathBuf::from(&m.bin_dir).join(exe_filename("nyx-mcp", backend.os()));
        for h in nyx_connect::list(&env, &mcp_exe) {
            report.ai_hosts.push(HostOutcome {
                id: h.id.to_string(),
                action: "doctor".to_string(),
                status: connect_hook::status_label(&h.status).to_string(),
                config_path: Some(h.config.display().to_string()),
                backup_path: None,
                detail: None,
            });
        }
    }
    Ok(report)
}

/// Repairs an existing installation from `source`, or restores missing/corrupted files cleanly.
pub fn run_repair(
    source: &PayloadSource,
    opts: &TxOptions,
    backend: &dyn PlatformBackend,
    inspector: &dyn ProcessInspector,
    paths: &HostPaths,
    on_progress: &mut dyn FnMut(TxProgress),
) -> Result<InstallReport> {
    let mut repair_opts = opts.clone();
    repair_opts.force = true;
    run_install_or_upgrade(
        "repair",
        source,
        &repair_opts,
        backend,
        inspector,
        paths,
        on_progress,
    )
}

/// Helper used by `tests/hermetic_tx.rs` to verify every field of a live `InstallManifest` against disk.
pub fn collect_installed_file_map(manifest: &InstallManifest) -> BTreeMap<String, PathBuf> {
    manifest
        .files
        .iter()
        .map(|f| (f.rel_path.clone(), file_abs_path(manifest, f)))
        .collect()
}
