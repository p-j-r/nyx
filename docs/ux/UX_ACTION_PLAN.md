# nyx window: UX action plan (Rust, egui 0.37)

Date: 2026-11-17. Base: `productivity-4-memory` (`fb57860d`). Evidence: `docs/ux/UX_REVIEW.md`. Decision: ADR 90 (provisional number). Tokens: `docs/ux/tokens.json`, checked by `docs/ux/contrast.py`. Mockups: Paper file, pages "nyx UX" (see the end of this file).
Nothing in this plan has been implemented; this branch changes documentation only. The Paper mockups (section 6) are still to be made: Paper Desktop was not running.

## 0. Principles of the plan

* **Smallest change that removes a measured problem**, in the order the skill sets: accessibility, then conventions, then craft, then polish.
* **The engine is the single source of actions.** The window already calls commands by name (`Engine::catalog()`); the palette, the menu, the shortcuts help and the MCP tools should all be generated from it, so no front end can do something another cannot (`CLAUDE.md`, ADR 25, 57).
* **No new dependency unless it earns its place.** Everything in U1 to U4 is plain egui 0.37 and `std`. `egui_extras` is not used (it is not a dependency, and `egui::Grid` is enough for the dashboard table). The one candidate (`muda`, native menu) is a decision for Paul (D1).
* **Each step ends green:** `cargo test --release -p nyx`, `cargo clippy --release -p nyx --all-targets --no-deps -- -D warnings`, `bash docs/architecture/check_code_to_adr.sh`, and, when a source file is added or removed, `bash docs/architecture/c4/tools/build.sh`. `crates/app` is not rustfmt-formatted (ADR 65): do not run `cargo fmt` on it.
* **Baseline first (U0).** The productivity status records two failing `nyx` window tests (`every_tab_draws_for_the_examples`, `the_window_checks_convergence_on_its_own_thread`: "Dropped TexturesDelta with 1 unapplied deltas"). Record the real baseline before changing anything so no step is blamed for them.

## 1. What changes relative to earlier plans

| Earlier document | Said | Now |
|---|---|---|
| `docs/WORKFLOWS_DESIGN.md`, `docs/mockups/*.png` | A light, warm app: amber "Running" cards with a live mini-plot, plum fills, sparklines | The shipping window defaults to dark on a dark Windows and has no Running card; the mockups show light only. The new mockups show both appearances, the quiet accent rule, and the running, queued and stale states that were never drawn. |
| `docs/NYX2_PRODUCTIVITY_PLAN.md` (P1 stale, P2 jobs, Q1 queue, sweeps) | Engine-side time-to-trustworthy-result | The window surfaces are the weak half: the Queue panel and stale band were "not looked at on screen", the window has no sweep UI, and the Sweep card still says "coming soon". U5 draws them. |
| `docs/GEOMETRY_GUI_PLAN.md` | "No numbers in the interface" | Kept. U3 extends it to metric names and formatting. |
| ADR 65 (fmt/clippy) | `app` is not rustfmt-formatted | Kept; new code follows the surrounding style, clippy clean. |

## 2. Phases

Sizes: S = under a day, M = one to three days, L = a week. Files are in `crates/app/src/`.

### U0. Baseline and scaffolding (S)

* Record test counts and the two known failures for `-p nyx`; record a scripted shot set (`tools/ux-shots.sh`, new, see U6) of the current look in both appearances: the "before".
* Add `docs/ux/` files to the status file; branch `ux-0-review` (this commit).
* **Gate:** the baseline is written into `docs/UX_STATUS.md`.

### U1. Accessibility and the token set (M) , the Critical findings C1 to C4

**New `tokens.rs`** (no egui types; pure and testable): `pub struct Tokens { page, card, raised, separator, control_edge, text, muted, accent, on_accent, tint, activity, on_activity, activity_text, activity_tint, success, danger, canvas: [u8; 3] }`, `pub const LIGHT: Tokens`, `pub const DARK: Tokens` with the values of `tokens.json`, `pub fn contrast(a, b) -> f32` (WCAG relative luminance), and the unit test:

```rust
#[test]
fn the_tokens_meet_their_contrast_floors() {
    for (name, t) in [("light", &LIGHT), ("dark", &DARK)] {
        let text_on = |fg: [u8; 3], bg: [u8; 3]| contrast(fg, bg);
        for bg in [t.page, t.card, t.raised] {
            assert!(text_on(t.text, bg) >= 7.0, "{name} text");
            assert!(text_on(t.muted, bg) >= 4.5, "{name} muted");
            assert!(text_on(t.success, bg) >= 4.5 && text_on(t.danger, bg) >= 4.5, "{name} status");
        }
        assert!(contrast(t.on_accent, t.accent) >= 4.5, "{name} accent fill");
        assert!(contrast(t.accent, t.page) >= 4.5, "{name} accent text");
        for bg in [t.page, t.card] { assert!(contrast(t.control_edge, bg) >= 3.0, "{name} edge"); }
        assert!(contrast(t.on_activity, t.activity) >= 4.5 && contrast(t.activity_text, t.activity_tint) >= 4.5);
    }
}
```

**`theme.rs`** keeps its public functions (`page`, `card`, `line`, `text`, `muted`, `accent`, `tint`, `warm`, the `AMBER`... constants) as thin wrappers over `tokens` so the 100-odd call sites do not change in one go; adds `control_edge(ctx)`, `on_accent(ctx)`, `signal(ctx, Signal::{Success, Danger, Activity})`, `canvas(ctx)`. `apply_mode` builds both visuals once and installs them with `ctx.set_visuals_of(Theme::Light, ..)` / `Theme::Dark` (egui 0.37) and sets `widgets.*.bg_stroke` to `control_edge`, `noninteractive.bg_stroke` to `separator`, `selection.stroke` to a 2 px accent (the focus ring).

**`kit.rs`:** `primary` and `pill` use `theme::on_accent`; selected `pill` becomes tint + 2 px leading bar unless it is the one primary; `danger` uses `signal(Danger)`; `dot` is replaced by `status(ui, State, label)` that draws a glyph (filled circle saved, half-filled running, ring draft, triangle stale, cross failed) and the word; new `tab(ui, title, selected, state) -> TabResponse` (>= 28 px high, `ui.interact`, 2 px focus ring on `has_focus()`, `response.widget_info(|| WidgetInfo::labeled(WidgetType::Button, true, title))`, a separate 24 x 24 close target with tooltip).

**`shell.rs`:** the tab strip uses `kit::tab`; `Search or run anything Ctrl K` stays; the status chip keeps its place (U3 trims it).

**Lint test** (`theme.rs` tests): walk `crates/app/src/*.rs`, fail on `Color32::from_rgb`, `from_rgba`, `Color32::WHITE/BLACK/RED/GREEN...` outside `theme.rs`, `tokens.rs`, `colormap.rs` and an explicit allow-list (`viewport.rs` axis colours, until U4). Start the allow-list with today's 56 literals, then shrink it in U3 and U4. This is how "one colour, one meaning" is enforced.

* **Gate:** tokens test, lint test, `kit::tab` reachable with Tab/Shift+Tab/Enter in a window test, scripted shots before/after, clippy.
* **Rollback:** `tokens.rs` can be reverted to the old values with no widget change, because call sites go through wrappers.

### U2. Conventions (M): C-level follow-ups and I-1, I-6, I-10

* **Appearance follows the system (I-1).** `Settings.theme` default `"system"`; startup `ctx.set_theme(ThemePreference::System)`; remove the top-bar Light/Dark `Button` (`shell.rs`); Settings keeps Follow system / Light / Dark. Test: with `ThemePreference::System` and a faked `system_theme`, the visuals switch without restart.
* **Commands in one place (I-6).** New `menu.rs`: `fn commands() -> Vec<Cmd>` built from `Engine::catalog()` plus window-only actions (`New workflow`, `Open .nyx`, `Settings`, `Compute`, `Next/Previous tab`, `Zoom in/out`). The Ctrl+K palette and an in-window menu strip (File, Edit, View, Run, Help) both render from it; every toolbar button has an entry; shortcuts are printed from the same table (`Ctrl+,` Settings, `Ctrl+1/2/3` sections, `Ctrl+Tab` tabs, `F5`, `Ctrl+W`, `Ctrl+K`). The classic layout's menus move onto it, so the menu exists in both layouts. **Decision D1:** in-window strip (no dependency, same on every OS) or native menu bar via `muda` (a Windows menu bar and, later, the macOS one; adds a dependency and an event bridge). Recommendation: in-window strip now, native later behind a feature.
* **Settings** leave the toolbar for the menu (`Ctrl+,`), keep the page.
* **Scale (I-10).** One `type_scale`: Small 13, Body 14, Heading 17, Title 28; `interact_size.y = 28`; remove the 12 / 12.5 / 11 / 11.5 / 10.5 literals; `ctx.set_zoom_factor` bound to `Ctrl +`, `Ctrl -`, `Ctrl 0` and saved in settings (this is the "text scales" answer).
* **Settings additions:** Contrast Standard / High, Reduce motion, Canvas follows theme / always dark.
* **Gate:** palette lists every engine command (test compares to `Engine::catalog()` names); no top-bar theme button (window test); shots.

### U3. Craft and copy (M): I-7, I-8, I-11, I-12, I-13

* **Title** takes the remaining width with an ellipsis and a hover (`wfview.rs`); `+ tag` moves behind it.
* **Copy pass** through one `strings.rs`-style set of helpers: `plural(n, "run")`, `si(value, unit, sig = 4)` ("0.9447", "500 K", "1.55 um"), metric labels from the physics descriptors (`eigenvalue_re` -> "Eigenvalue (real)"), default plot picks the first metric that varies, never a constant; remove "Model m1 (not yet functional)" and the roadmap lines from the palette; the Sweep card becomes real in U5.
* **Dashboard table (I-11):** `egui::Grid` with fixed columns (title, summary, status, when); summaries never say "geometry not recognised" (fall back to the physics and study names).
* **Idle chip and tagline (I-12):** the chip appears only when a job runs or waits; "creative chaos." moves to the About sheet.
* **Sheets (I-13):** a Close control at the top right of every `kit::modal`.
* **Gate:** window tests for the helpers (plural, si), shots, lint allow-list shrinks.

### U4. Data colour and the canvas (M): I-3, I-4

* **`colormap::palette(i)`** becomes the validated eight per appearance (`tokens.json` `data`), `palette_hatch(i)` returns the hatch pattern (diagonal, cross, dots) for regions 9 and up; the viewport draws the region number on every region already, so the light-mode contrast WARN has its relief.
* **Canvas follows the theme:** `viewport.rs` and `app.rs` (`Color32::from_rgb(24, 26, 30)`) take `theme::canvas(ctx)`; the mesh-edge and axis colours are derived from `text` and `muted`; axis x/y/z colours come from the palette, not literals; the "Always dark canvas" setting restores today's look.
* **Colormaps:** add Cividis (table of 11 stops as for Viridis), label Turbo "Rainbow (not recommended)", keep Viridis default.
* **Legend:** a swatch list (colour + hatch + number + name) beside the viewport, for regions and boundaries, so colour is never the only code.
* **Gate:** a unit test that runs the same checks as the dataviz validator on the eight (lightness band, chroma floor, adjacent OKLab dE >= 15 normal-vision, a CVD simulation with the Machado matrices; about 60 lines, the matrices are published constants) so the palette cannot regress; shots in both appearances.

### U5. Productivity surfaces (L): the window half of the productivity plan

Design first (the mockups below), build after Paul agrees the shape.

* **Running state** on the workflow tab (activity glyph + word), in the status chip (`1 running - Solving 14 / 21 points`), and as an amber card on the Dashboard with Cancel (as `docs/mockups/01` intended).
* **Queue panel** (exists in code, never seen): a list with kind, owner, priority, position, Hold / Release / Cancel, Pause queue; a single "Compute" that becomes "Queue" while busy, with the position in the label.
* **Stale** band: one place (top of the viewport), activity colours, a "Solve again" primary only when nothing else is primary. Reproduce I-9 first.
* **Sweep:** a Sweep section in Build that writes `sweep_run` points (series: command, field, values; or a list of edit lines), shows the plan ("6 meshes, 6 solves"; the engine's dry run from `nyx check`), the table of points as runs fill in, and one chart of the metric against the parameter. The plan and the table are the engine's, not new logic.
* **Gate:** window tests for each state with a fake engine job; shots with the job held (use `Engine` hold in the test, since the script cannot shoot while busy).

### U6. Verification harness (S, start in U0, grow every phase)

* `tools/ux-shots.sh`: runs `nyx --script` into a scratch `NYX_WORKSPACE` for 8 views x light/dark x (1440 x 900, 1100 x 700), writes to `target/ux-shots/`; used for the review and as a smoke test (no panic, files non-empty). It never touches the user's workspace (this review's first run did, see `docs/UX_STATUS.md`).
* Review checklist per phase = the skill's quality floor: responsive down to the 900 x 600 minimum, visible keyboard focus on every control, reduced motion and High contrast respected, `Ctrl +` to 200 percent survivable.
* A manual Narrator and High Contrast pass before U2 closes (the review did neither).

## 3. Order, parallelism and workspaces

U0, then U1 alone (everything depends on tokens). After U1, three writers can run at once in their own jj workspaces (`jj workspace add ../nyx2-ux-<name> -r <ux-1 bookmark>`, `CARGO_TARGET_DIR` each), because their files do not overlap:

| Workspace | Phase | Files |
|---|---|---|
| `nyx2-ux-conv` | U2 | `menu.rs` (new), `shell.rs`, `sheets.rs` (Settings), `main.rs`, `app.rs` shortcuts |
| `nyx2-ux-craft` | U3 | `dashboard.rs`, `wfview.rs`, `tabs.rs` strings, `kit.rs` helpers only for `plural`/`si` (new `fmt.rs` instead, to avoid `kit.rs`) |
| `nyx2-ux-data` | U4 | `colormap.rs`, `viewport.rs`, `scene.rs`, `app.rs` canvas lines only |

U5 after U2 to U4 land (it needs the menu, the formats and the legend). Bookmarks: `ux-1-tokens`, `ux-2-conventions`, `ux-3-craft`, `ux-4-data`, `ux-5-surfaces`; one commit line each ("what changed and why"). Never push (no remote); never commit to `main`. After every phase that adds or removes a source file: add its `CODE_TO_ADR.md` row, `bash docs/architecture/check_code_to_adr.sh`, `bash docs/architecture/c4/tools/build.sh`, and update `docs/UX_STATUS.md`. Every new file starts with a one-paragraph comment of what it does.

## 4. Decisions for Paul

* **D1** Menu: in-window strip now (recommended) or native `muda` menu.
* **D2** Canvas default: follows the theme (recommended) or stays dark.
* **D3** Dark default: follow the system (recommended; the app currently starts in whatever was saved last, which on this machine was dark).
* **D4** Accent in dark: `#E58BB8` (recommended, lighter and calmer than `#D96AA4`) or keep brand pink and use ink text only.
* **D5** Contrast "High" as a setting only (no OS detection in egui) or add a Windows registry read for the system High Contrast flag (a `cfg(windows)` block; small, but it reaches into the OS).
* **D6** Is the `STALE` banner on a fresh result a bug (I-9)? Needs a hand repro.

## 5. Risks

* **Token swap changes the look of every screen at once.** Mitigation: wrappers keep the call sites; the before/after shot set is the review; the old values stay in git and in `docs/ux/contrast.py`'s "current" block.
* **egui limits:** no native menu, no OS reduced-motion or increase-contrast query, hit-testing of overlapping responses is by registration order (the tab close target must be allocated after the tab body). Each is handled above or listed as a decision.
* **Colour tests prove contrast, not comfort.** The comfort claim (less visual stress) rests on the measured chroma and area changes in the review (section 4), not on a user study; ask two people who use it daily for one week with the quiet theme, and keep a "Previous colours" switch for that week only (a theme variant, deleted when D4 is settled).
* **Scope creep into the solver docs.** None of this touches physics, goldens or contracts. `engine_commands` and `mcp_tools` contract versions do not change in U1 to U4.

## 6. Mockups (Paper)

Planned, page "nyx UX" in the Paper file: tokens and roles; Dashboard light and dark with a running card; Build with viewport (follows theme) and the legend; Analysis with the quiet selection; states sheet (saved, running, queued, stale, failed with glyph and word); Queue panel; Settings with the new options. Every colour in them is a value of `docs/ux/tokens.json`; egui cannot blur or cast the soft shadows, so none is drawn.
