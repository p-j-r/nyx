# 130. WSLg display stability, automatic high-DPI scaling, layout polish, and deterministic shell script driving

Date: 2026-11-25

## Status

Accepted, implemented on `feat/gui-wslg`. Extends ADR 14, ADR 39, ADR 88, ADR 92, and ADR 124.

## Context

Desktop verification on Linux/WSLg high-DPI displays (such as 3840x2160 reporting `native_pixels_per_point = 1.0` over Wayland with software rendering) and interactive/scripted UI passes uncovered seven concrete gaps:

1. **WSLg Wayland startup instability and noisy logs.** `winit 0.30.13` dropped built-in `WINIT_UNIX_BACKEND` inspection and connects to Wayland before X11, frequently hitting `Broken pipe (os error 32)` on WSLg's Wayland socket; first-frame window visibility, focus, and monitor clamping were not enforced either.
2. **Missing high-DPI auto-scaling and broken zoom persistence.** Body text stayed at 1x on 4K WSLg displays without manual `Ctrl+=` presses, `--size WxH` acted in raw OS pixels rather than logical UI points, there was no `NYX_SCALE` env override or Settings "Interface size" preset, and `attach_workspace` invoked `sync_look` before `Context::set_zoom_factor` took effect on the next `begin_pass`, immediately overwriting `Settings.zoom` on disk back to `1.0`.
3. **Stale compositor/swapchain contents on resize/maximize/scale/theme change.** `App` did not override `eframe::App::clear_color` (whose default alpha `180/255` left semi-transparent remnants), did not configure `PresentMode::Fifo` under WSLg, and did not issue multi-frame repaint bursts on viewport transitions.
4. **Top-bar and workflow-bar click reliability.** Tight fixed `right_reserve` / stacking thresholds allowed open tabs or wide left menus to crowd right-aligned buttons (`Settings`, `+`, `+ New workflow`, `Compute F5`, `Script`), while `PointerButton` arrival did not request an explicit repaint for single-frame press/release delivery.
5. **Swept 3D cards.** `nyx_engine::store_link::card_of` ignored non-empty `MeshSettings.sweep` on 2D geometries, displaying `2D` on Dashboard cards for extruded 3D models.
6. **Inspector clipping, low weak-text contrast, card stretching, and fractional mode ticks.** `Visuals::weak_text_color` fell back to translucent `text_color`, the `element size h (m)` muted hint sat on a non-wrapping horizontal row, `summary_plot` emitted fractional ticks (`1.1`, `1.2`, ...) for discrete modes, single Dashboard cards stretched across the full width, and `wf_inspector` stayed fixed at 400/640 max.
7. **Non-deterministic `--script edit` waits and legacy-only script commands.** `script_step` used a fixed 40-frame countdown after `edit` instead of waiting on `Feedback` timers/threads, and `new`/`example`/`tab`/`physics`/`compute` acted on legacy single-project state rather than the ADR 124 tab shell.

## Decision

* **WSLg backend selection, FIFO surface, and first-frame geometry (`app::main`, `shell::is_wsl`).** Detect WSLg (`WSL_DISTRO_NAME`, `WSL_INTEROP`, or `/mnt/wslg` on Linux), honour `WINIT_UNIX_BACKEND` explicitly via `EventLoopBuilderExtX11`/`EventLoopBuilderExtWayland`, default unset WSLg launches with `DISPLAY` to X11/XWayland, configure `PresentMode::Fifo` with single-frame latency (`desired_maximum_frame_latency = Some(1)`) under WSLg, emit one plain-language stderr notice only if WSLg selects a CPU/software adapter or fails backend init, set `ViewportBuilder` visible/active/non-transparent/monitor-clamped, convert logical `--size` to OS dimensions via resolved startup scale, clamp to usable monitor area, and maximize interactively when exceeded.
* **High-DPI auto-scale, `NYX_SCALE`, `scale_mode`, and persistent zoom (`nyx_store::Settings`, `app::{shell,sheets,menu,app}`).** Add backward-compatible `scale_mode: String` (`"auto"`, `"100%"`, `"125%"`, `"150%"`, `"175%"`, `"200%"`, default `"auto"`) to `nyx_store::Settings`. Resolve scale via `NYX_SCALE` (`0.5..=3.0`) > `scale_mode` > display pixel size when `native_ppp <= 1.15` (picking `2.0` on 4K `>= 3200`, `1.5` on `>= 2500`, `1.25` on `>= 2048`), scaled by user `Settings.zoom`. Track `pending_zoom` across passes so `sync_look` never overwrites saved `zoom` before `begin_pass`, intercept `Ctrl+=`/`-`/`0` directly (`zoom_with_keyboard = false`) including French/AZERTY modifiers, add "Interface size" in Settings > Units and look, and document WSL's `/home/<user>/Documents/nyx` root in Settings > Where things live.
* **Opaque clear and surface repaint bursts (`app::{app,theme}`).** Override `eframe::App::clear_color` with opaque `theme::page_for_dark`, paint `theme::page` across `ctx.content_rect()` on `LayerId::background()`, request a repaint immediately on `Event::PointerButton`, and run a 4-frame repaint burst whenever `inner_rect`, `maximized`, `pixels_per_point`, or `dark_mode` changes.
* **Layout overlap fixes and visual polish (`app::{shell,dashboard,wfview,tabs,theme}`).** Dynamically compute top-bar right reserve, wrap Dashboard header controls below 820 pt and clamp workflow cards to `<= 420.0` pt in a responsive grid, base `wf_bar` inline vs stacked rows on remaining width after menus, scale `wf_inspector` default/max width with zoom and available width, set `Visuals::weak_text_color = Some(tokens::MUTED)`, wrap Geometry `element size h (m)` and Study doc/wavelength notes onto clean lines, and attach `integer_mode_marks` + `format_integer_mode_tick` to `summary_plot` for eigenfrequency and waveguide mode results.
* **Swept 3D concept cards (`nyx_engine::store_link::card_of`).** Report `dim = 3`, `"3D (swept)"` in `summary`, and `"3d"`/`"swept"` search keywords when `geometry.dim == 2` and `mesh.sweep` is non-empty.
* **Deterministic script waits and full ADR 124 shell driver (`app::{feedback,script}`).** Add `Feedback::is_busy()`, hold `script_step` while `!script_free && feedback.is_busy()` while ticking `feedback_tick(ctx)`, route `new`, `example`, `open`, `physics`, `tab`, `compute`, `page`, `sheet`, and non-persistent `scale` cleanly through `Shell`, and keep unit tests covering pointer clicks, repaint bursts, scaling/persistence, mode ticks, swept cards, and shell scripting.

## Consequences

* First launch on 4K WSLg displays scales legibly without manual zoom keypresses, `--size` behaves in logical UI points, window resize/maximize leaves no stale borders, right/top toolbar buttons remain clickable across all widths and scales, and `--script` regressions can exercise the full ADR 124 shell deterministically.
* `nyx_store::Settings` gains backward-compatible field `scale_mode` (`"auto"` default); existing `store_workflow.v1` and `store_run.v1` contract exemplars and `docs/ux/tokens.json` are unchanged.

## Code

`crates/app/src/{main,app,shell,dashboard,wfview,sheets,tabs,theme,menu,feedback,script}.rs`, `crates/app/Cargo.toml`, `crates/nyx-store/src/workspace.rs`, `crates/nyx-engine/src/store_link.rs`. Tests: `store_link::tests`, `workspace::tests`, `tabs::tests`, `shell::tests`, `script::tests`, `app::tests`.
