// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cross-platform pure-Rust installer, packager, updater, health checker (`doctor`), and uninstaller for `nyx` (ADR 122).
//!
//! Architecture overview:
//! * [`bundle`] — streaming `NYXPK001` archive reader/writer, 96-byte PE-Authenticode-aware `NYXINST1` trailer, `.zip`/dir/`.app` sources, zip-slip/bomb guards
//! * [`manifest`] — versioned `install-manifest.json` (`format = 1`, `install_manifest.v1`) and JSON Schema
//! * [`report`] — machine-readable `--json` report (`install_report.v1`) and JSON Schema
//! * [`safety`] — poison-prefix denylist and BLAKE3-verified non-recursive uninstall cleanup (`remove_dir_all` never called on install roots)
//! * [`process`] — running-binary detection and cooperative ADR 95 shared workspace queue pause/drain
//! * [`connect_hook`] — AI host MCP integration (`nyx-connect` install, upgrade, surgical uninstall, and JSON-RPC probe)
//! * [`platform`] — trait-abstracted Windows, Linux, and macOS desktop/PATH/shortcut/file-association drivers plus in-memory hermetic fake backend
//! * [`tx`] — 7-phase atomic transaction engine with append-only LIFO rollback journal
//! * [`notices`] — offline `cargo metadata` SPDX collector generating `THIRD_PARTY_NOTICES.txt`
//! * [`cli`] — hand-parsed CLI and TTY fallback shared by `nyx-setup` and `nyx-uninstall`
//! * [`tokens`], [`gui`] — optional (`feature = "gui"`) lightweight `eframe`/`glow` wizard following ADR 88 colour tokens

pub mod bundle;
pub mod cli;
pub mod connect_hook;
pub mod manifest;
pub mod notices;
pub mod platform;
pub mod process;
pub mod report;
pub mod safety;
pub mod tx;

#[cfg(feature = "gui")]
pub mod tokens;

#[cfg(feature = "gui")]
pub mod gui;
