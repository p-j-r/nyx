// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-user `nyx-setup` / `nyx-uninstall` multi-role entrypoint.
#![cfg_attr(all(windows, feature = "gui"), windows_subsystem = "windows")]

fn main() {
    #[cfg(windows)]
    {
        nyx_install::platform::windows::attach_parent_console_if_needed();
    }
    std::process::exit(nyx_install::cli::entry());
}
