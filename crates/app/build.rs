// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Build script: puts the nyx mark into the Windows executable as its icon (what Explorer and the taskbar show for the file). Other platforms have
//! nothing to do here; the window icon itself is set at run time from `docs/brand/nyx-mark.png` (`main.rs`).

fn main() {
    #[cfg(windows)]
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        // an icon is a nicety: when the resource compiler is missing the build goes on without it
        let _ = embed_resource::compile("assets/nyx.rc", embed_resource::NONE).manifest_optional();
    }
}
