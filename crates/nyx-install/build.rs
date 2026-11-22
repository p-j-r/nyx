// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Build script: embeds the nyx icon and `asInvoker` execution-level manifest into Windows `nyx-setup` / `nyx-pack` executables
//! (matching `crates/app/build.rs` and preventing Windows' filename-based UAC elevation prompt for `setup` binaries).

fn main() {
    #[cfg(windows)]
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        let _ = embed_resource::compile("assets/nyx-setup.rc", embed_resource::NONE)
            .manifest_optional();
    }
}
