# 122. Cross-platform pure-Rust installer, packager, and uninstaller (`nyx-install`)

Date: 2026-11-22

## Status

Accepted.

## Context

`nyx` ships five cooperating native executables (`nyx`, `nyx-mcp`, `nyx-worker`, `nyx-script-worker`, `nyx-connect`) plus read-only runtime asset directories (`materials/`, `examples/`, `fonts/`) and license notices (`LICENSE`, `NOTICE`, `OFL-EBGaramond.txt`, `THIRD_PARTY_NOTICES.txt`). End users previously relied on developer builds or ad-hoc `dist/` folders (`tools/package_dist.sh`), leaving OS `PATH`, desktop launchers, file associations, AI-host MCP setup (`nyx-connect`), upgrades, and uninstallation manual.

Several runtime assumptions in existing crates also assumed a flat developer directory:

* `theme::first_font()` checked only `<exe_dir>/fonts`, missing split `<prefix>/fonts`, macOS `../Resources/fonts`, and Linux `../share/nyx/fonts`.
* `materials::default_dirs()`, `examples::dirs()`, `nyx_script::find_worker()`, and `nyx_connect::find_server_with()` needed `current_exe().canonicalize()` to follow `~/.local/bin` / `/usr/local/bin` symlinks on Unix/macOS and split-layout ancestor directories.
* `materials::load_library()` stopped at the first non-empty directory while `material_fetch::default_materials_folder()` attempted to write fetched `.toml` files back into the bundled `materials/` folder — failing in read-only system locations or breaking macOS bundle signatures.
* `nyx_connect::disconnect()` removed the `"nyx"` entry even when pointing to a different installation.
* `nyx` needed headless `--version` / `--help` flags, macOS `-psn_*` filtering, full file-open handling for `.nyx` archives and `.nys` scripts (opened in `ScriptSheet` without auto-executing), and `set_self_as_worker(false)` in `runfile` / GUI script execution.

## Decision

* **Single self-contained installer crate (`crates/nyx-install`, layering `["nyx-connect"]`):**
  * Exposes a multi-role user binary `nyx-setup` (also bundled inside every install as `nyx-uninstall` with CLI subcommands `install`, `upgrade`, `repair`, `uninstall`, `status`, `doctor`) and release packager `nyx-pack`.
  * Uses only crates already locked in `Cargo.lock` with `default = []` so `--workspace` builds remain fast and free of `wgpu`; enables an optional lightweight 640x480 ADR 88 GUI wizard via `--features gui` (`eframe 0.37` with `glow`, `default_fonts`, `x11`, `wayland`, `rfd`), compiled separately (~5–6 MB).
  * Embeds an `asInvoker` application manifest on Windows (`assets/nyx-setup.manifest`) and hardens DLL lookup via `SetDefaultDllDirectories(LOAD_LIBRARY_SEARCH_SYSTEM32)`.

* **Streaming `NYXPK001` archive & PE Authenticode-aware 96-byte `NYXINST1` trailer (`bundle.rs`):**
  * Per-file streaming `flate2`/`zlib-rs` compression with strict path-traversal / ADS / reserved DOS name guards, per-entry and total decompression caps, full-stream + per-file BLAKE3 verification, and offline SPDX `THIRD_PARTY_NOTICES.txt` generation (`notices.rs`).
  * Self-extracting executables pad `0..7` zero bytes before an 8-byte-aligned 96-byte `NYXINST1` footer so PE `DataDirectory[SECURITY]` certificates sit cleanly after the trailer; runtime reads `cert_offset - 96` when signed or `EOF - 96` otherwise.

* **Canonical OS layouts & non-destructive user-data separation (`platform/`, `safety.rs`, `tx.rs`):**
  * **Windows**: `%LOCALAPPDATA%\Programs\nyx\` (user) / `%ProgramFiles%\nyx\` (`--system`) with `bin\*.exe` on `Path` via raw UTF-16 `RegQueryValueExW`/`RegSetValueExW` (`REG_EXPAND_SZ` preserved, `WM_SETTINGCHANGE` broadcast), unquoted `DisplayIcon` + quoted `UninstallString`/`QuietUninstallString`/`ModifyPath`, COM `.lnk` shortcuts, safe `OpenWithProgids` for `.nyx.json`/`.h5`, and synchronous Mode B `--_reap-final` / GUI Mode A self-deletion.
  * **Linux**: `${XDG_DATA_HOME:-~/.local/share}/nyx/` (user) / `/usr/local/lib/nyx/` (`--system`) with atomic symlinks in `~/.local/bin` or `/usr/local/bin`, marked shell blocks, XDG `.desktop`, `hicolor` PNG icons, `shared-mime-info` XML (`*.nyx.json` `weight="60"`), and two-tier `dlopen` runtime checks.
  * **macOS**: `~/Applications/Nyx.app` / `/Applications/Nyx.app` with generated `.icns` and `Info.plist`, external receipt `install-manifest.json` in `Application Support/nyx/install/` (or `<prefix>/.nyx-install/`) so `Contents/_CodeSignature` is never invalidated, zero Mach-O byte appending, and `lsregister` hooks.
  * **Atomic 7-phase transactions & safe uninstall**: Backed by `install-journal.jsonl` LIFO rollback, cooperative ADR 95 shared-queue pause, `nyx_connect::disconnect_matching` before process termination, and BLAKE3-verified bottom-up `remove_dir` cleanup (`remove_dir_all` banned on install roots; `Documents/nyx` never deleted).

## Consequences

* Fresh installs, upgrades, repairs, AI-host MCP connections, and uninstalls are transactional, scriptable (`install_report.v1` JSON), and verifiable hermetically (`FakePlatformBackend`) as well as end-to-end on clean empty working directories.
* Bundled materials stay cleanly read-only while user custom/fetched materials overlay smoothly from `Documents/nyx/materials` via `Library::with_texts`.

## Code

`crates/nyx-install/*`, `tools/package_dist.sh`, `crates/nyx-engine/src/{materials,material_fetch}.rs`, `crates/nyx-script/src/lib.rs`, `crates/nyx-connect/src/{config,lib,server}.rs`, `crates/app/src/{main,app,examples,theme,shell,sheets,runfile}.rs`, `docs/architecture/contracts/install_{manifest,report}.v1.*`.
