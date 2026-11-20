// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Where things live on this machine: the home folder, the Windows application-data folders and the operating system.
//! [`Env::detect`] reads the process environment; tests build an `Env` by hand and never touch the real home folder.
use std::path::PathBuf;

#[derive(Clone, Debug)]
pub struct Env {
    pub home: PathBuf,
    /// `%APPDATA%` on Windows.
    pub appdata: Option<PathBuf>,
    /// `%LOCALAPPDATA%` on Windows.
    pub local_appdata: Option<PathBuf>,
    pub os: Os,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Os {
    Windows,
    Mac,
    Linux,
}

impl Os {
    /// The operating system this program was built for.
    pub fn current() -> Os {
        if cfg!(windows) {
            Os::Windows
        } else if cfg!(target_os = "macos") {
            Os::Mac
        } else {
            Os::Linux
        }
    }

    /// Lower-case name used in the reports and the contract exemplar.
    pub fn name(self) -> &'static str {
        match self {
            Os::Windows => "windows",
            Os::Mac => "mac",
            Os::Linux => "linux",
        }
    }
}

impl Env {
    pub fn detect() -> Result<Env, String> {
        let var = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let home = var("USERPROFILE").or_else(|| var("HOME")).ok_or("cannot find the home folder (USERPROFILE or HOME)")?;
        Ok(Env { home, appdata: var("APPDATA"), local_appdata: var("LOCALAPPDATA"), os: Os::current() })
    }
}
