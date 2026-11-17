# nyx window: design review (apple-design lens)

Date: 2026-11-17. Reviewed: the window of `productivity-4-memory` (`fb57860d`, built binary of 2026-11-17 00:52), Windows 11, Rust 1.99, `eframe`/`egui` 0.37, release build, driven with `nyx --script` into a scratch workspace (`NYX_WORKSPACE`), 1440 x 900.
Method: the `apple-design` skill (Apple's HIG as design principles; nyx is a desktop Rust/egui app, so the platform pages apply as principles, Windows is not macOS). Every citation below was read from the skill's reference pages in this session. Contrast figures are computed from the hex values in `crates/app/src/theme.rs` and `kit.rs` (scripts and numbers in section 7). Palette checks were run with the `dataviz` validator.
The plan that follows from this review is `docs/ux/UX_ACTION_PLAN.md`; the colour decision is ADR 90 (provisional number).

## 0. What was and was not looked at

| Looked at (screenshots in `docs/ux/shots/`) | Not looked at |
|---|---|
| Dashboard, Concept, Build, Analysis, the result viewport (heat, solid), Start-a-workflow sheet, Ctrl+K palette, Share sheet, Settings; light and dark | A job **running** (the script waits while busy, so a shot cannot be taken; the PowerShell route is blocked by the machine's execution policy and was not bypassed). The amber running chip, the Queue panel and the Compute/Queue buttons were therefore read from code and from `docs/mockups/05-live-run.png`, not seen. |
| Code: `theme.rs`, `kit.rs`, `shell.rs` (top bar), `colormap.rs`, `app.rs` and `tabs.rs` colour use | The classic layout (opt-in, `Use the classic layout`). Narrator / screen-reader behaviour, Windows High Contrast, 150 percent display scaling. |

Limits to keep in mind: the sheets (Start a workflow, palette, Share) were captured within a few frames of opening and show egui's fade-in (text of the page behind shows through). That is not a finding. A "STALE" banner was seen on a freshly computed heat result (`light-heat-result.png`); it could be a script artefact and needs a manual repro before it is called a bug (section 3, I-9).

## 1. Summary

**Rating: Critical issues** (by the skill's rule: any accessibility failure is Critical). Underneath, the design is closer to "Good": the structure is calm and specific, and the failures are concentrated in one file (`theme.rs`) plus a handful of widgets in `kit.rs`.

**Thesis.** *A workbench where a simulation is something you can leave and come back to*: Dashboard of workflows, Concept / Build / Analysis as one-line-purpose sections, autosave stated in words, one primary action (Compute, F5) at the trailing end.
**What it will be remembered by:** the italic serif page titles (EB Garamond) over a quiet neutral UI, with plum as the brand. That is a choice rooted in the product (a lab notebook, not a dashboard), not a template.
**What lets it down:** colour is spent in too many places at once (up to four filled magenta objects in one dark view), the same colour means several things, the secondary colours were never made mode-aware, and the viewport is a black rectangle inside a light shell.

## 2. Critical (accessibility)

**C1. White text on the dark accent is 3.21:1.**
* *What:* `kit::primary` and `kit::pill` set `Color32::WHITE` on `theme::accent` = `#D96AA4` in dark mode. 14 pt text (egui sizes are logical points: 1 per pixel at 100 percent scale, as Apple's 1x) needs 4.5:1; it gets 3.21:1. The `.strong()` on the label does not buy the "bold" 3:1 allowance: in egui it changes the colour, not the weight, and only `Inter-Regular` is loaded (`theme::install_fonts`). This is the Compute button, New workflow, the selected Concept / Build / Analysis pill, every selected physics or study pill. In light mode the same pair is 9.93:1, so it only fails in dark, which is the default on a dark-themed Windows.
* *Why:* `accessibility.md › Vision`: contrast "Up to 17 pts, all weights: 4.5:1" (table); "make sure to check the minimum contrast in both light and dark appearances".
* *Fix:* in dark mode the accent *fill* takes ink text, not white: `#E58BB8` fill with `#1D1224` text = 7.49:1 (and the accent as text on the page = 7.72:1). One `kit::on_accent(ctx)` returns the text colour; both widgets call it.

**C2. Control edges are 1.2 to 1.4:1 against their surface.**
* *What:* `theme::LINE` `#E6E4E6` on `#F7F6F7` = 1.17:1, on white 1.26:1; dark `#3A2B42` on page `#170E1B` = 1.44:1, on card 1.34:1. `apply_mode` uses this one colour both as a decorative separator (correct) and as `widgets.inactive.bg_stroke` (every button, field, combo box, checkbox). In the screenshots the Search field, "Columns", the x:/y: combo boxes and the tag chip have no perceivable boundary.
* *Why:* `accessibility.md › Vision` ("enough contrast between foreground text and icons and background colors"); the 3:1 figure for non-text UI is WCAG 1.4.11, not a HIG number (labelled as such).
* *Fix:* two tokens: `separator` (keep) and `control_edge` (3.3:1: light `#8A8091` = 3.50 on page, 3.77 on white; dark `#6E6578` = 3.36 on page, 3.10 on card). Controls use `control_edge`; cards and dividers keep `separator`.

**C3. Status colours are not mode-aware, and status is carried by colour alone.**
* *What:* `theme::RED` `#B42B2B` is 2.98:1 on the dark page (the Cancel button text and stroke, `kit::danger`); `theme::GREEN` `#2E7D52` is 3.74:1 on the dark page (text uses of it fail 4.5). `AMBER` `#FFBD59` is 1.54 to 1.66:1 on light surfaces whenever it is used as a foreground (dot, icon). The tab dot (green = saved, amber = running, grey = draft) and the status chip carry state by hue only.
* *Why:* `color.md › Inclusive color`: "Avoid relying solely on color to differentiate between objects, indicate interactivity, or communicate essential information."
* *Fix:* `signal(ctx, Kind)` returns a mode-specific colour (dark: green `#5FC08C` 8.33:1, red `#F27A7A` 6.94:1, amber `#F2B35E` 10.08:1; light: green `#1F6B43` 6.02:1, red `#A32020` 7.01:1, amber ink `#6B4300` on `#FFF1DC` 7.77:1) and every state has a shape or a word next to it (section 5 of the plan: dot becomes glyph + label).

**C4. Custom tabs and the close button are labels, not controls.**
* *What (code-read, not tested with Narrator):* the workflow tabs and the "x" are `egui::Label` with `Sense::click()` (`shell.rs`, tab strip). They take keyboard focus but draw no focus ring, expose no `WidgetType::Button` role, and the "x" is a 16 px glyph, roughly 9 x 19 px of hit area.
* *Why:* `accessibility.md › Mobility`: macOS default control 28 x 28 pt, minimum 20 x 20 pt; "Let people use the keyboard alone to navigate and interact with your app."
* *Fix:* a `kit::tab()` widget: allocate a >= 28 px high rect, `ui.interact`, paint a 2 px focus ring when `has_focus()`, `response.widget_info(|| WidgetInfo::labeled(WidgetType::Button, true, title))`; the close control is its own 24 x 24 target with a tooltip.

## 3. Improvements

Each is tagged by the skill's severity.

**I-1 (High). App-specific appearance switch.** The top bar has a "Light"/"Dark" button (`shell.rs`, toggles `shell.dark` and writes `settings.theme`), and Settings has Follow system / Light / Dark. *Why:* `dark-mode.md › Best practices`: "Avoid offering an app-specific appearance setting ... they may think your app is broken because it doesn't respond to their systemwide appearance choice." *Fix:* default `theme = "system"` with live following (`ctx.set_theme(egui::ThemePreference::System)`, which exists in egui 0.37, and `ctx.set_visuals_of(Theme, Visuals)` for the two sets); remove the top-bar button; keep the three-way choice in Settings as the escape hatch.

**I-2 (High). Colour is spent on too many things at once, and one colour means several things.** (Section 4 has the full inventory.) In dark mode the Dashboard shows two filled magenta objects (the active tab and New workflow); a workflow view shows four (tab, section pill, sub-pill, Compute). *Why:* `color.md › Best practices`: "Avoid using the same color to mean different things"; `branding.md`: "Apply your app's accent color judiciously ... use it intentionally for primary actions or status indicators."

**I-3 (High). The viewport is always dark.** `app.rs` paints the canvas `#181A1E` in both modes (`Color32::from_rgb(24, 26, 30)`), so a light shell holds a near-black block (L 0.97 beside L 0.22). *Why:* `dark-mode.md`: a permanently dark appearance is for "an app that supports immersive media viewing". A plot canvas is not immersive media. *Trade-off:* a dark canvas makes colormaps pop and is a convention in FEM tools; the fix keeps it as a *choice* (Settings: Canvas follows theme / always dark) and defaults to following the theme.

**I-4 (High). The domain palette fails the validator.** 12 hand-picked colours in `colormap::palette` (`#569AD3`, `#E69138`, ...), drawn at alpha 0.45 to 0.55 over the canvas. Run on the dark canvas (first 8 slots): **FAIL lightness band** (5 of 8 above L 0.67), **FAIL chroma** (`#50BEBE` C 0.10), **FAIL CVD separation: adjacent `#D678AA` vs `#50BEBE`, deutan dE 4.3 (floor 6, target 8)**. Beyond 8 hues the method says never generate a ninth. *Fix:* the validated 8 (dark `#3987E5 #D95926 #199E70 #C98500 #D55181 #0CA30C #9085E9 #E66767`: all checks PASS, worst adjacent CVD 8.4; light set PASS with the contrast WARN, relief = the boundary and domain **numbers** that the viewport already draws) and, past 8, hatch plus number, not a new hue.

**I-5 (Medium). Competing hard-coded colours.** 56 `Color32::from_rgb/rgba/WHITE` outside `theme.rs` (app.rs 13, tabs.rs 13, sheets.rs 5, others): green Compute `(38,110,70)` in the legacy Study tab against plum Compute in the shell; blue Build mesh `(40,90,150)`; four different oranges for warnings `(200,120,20)`, `(235,170,60)`, `(230,190,80)`, `(235,110,100)`; reds `(200,60,60)`. Most are in `tabs.rs`/`app.rs` code that the classic layout uses (opt-in); some tab bodies are reused in the shell, so they can appear in either. *Fix:* no raw colour outside `theme.rs` and `colormap.rs`, enforced by a test (plan U1).

**I-6 (Medium). No menu bar; Settings sits in the toolbar.** *Why:* `the-menu-bar.md › Best practices`: "Provide app-specific menus for custom commands ... Even when commands are available elsewhere in your app, it's important to list them in the menu bar"; `toolbars.md`: "Make every toolbar item available as a command in the menu bar"; `settings.md`: "Include a settings item in the App menu. Avoid adding settings buttons to a window's toolbar". Mitigation that exists: Ctrl+K palette and the classic layout's menus. *Fix:* the engine already has one catalog of commands with schemas; the palette, an in-window menu strip and (optionally) a native Windows menu through `muda` should all be generated from it (plan U2). Decision for Paul: native menu or in-window strip.

**I-7 (Medium). Title is clipped mid-glyph.** `light-build.png`: "WR-90 rectangular waveguide mc" is cut with no ellipsis while the `+ tag` chip sits directly after it. *Fix:* give the title the width left after the right cluster and truncate with an ellipsis, full text on hover.

**I-8 (Medium). Copy that describes the build, not the task.** "Model m1 (not yet functional)"; the Sweep card "coming soon ... see docs/SWEEP_CACHING_PROPOSAL.md" while the engine has `sweep_run` since the productivity branch; "geometry not recognised 2D" in a dashboard card; "29979246 nm" as a wavelength in a card of the user's data; "1 runs plotted."; Analysis columns `EIGENVALUE IM` with 4 decimals and the default plot "eigenvalue_im against run" of a quantity that is zero. *Why:* `writing.md` (names from what people recognise); skill lens 5. *Fix:* labels from the physics descriptors, one `plural()` helper, SI-prefix formatting helper, no roadmap text in the nav (the palette's "Not yet:" list too).

**I-9 (Medium, needs repro). STALE banner on a fresh result.** `light-heat-result.png`: example 03, `compute` then `view result` shows "STALE: the model changed since this solve" with no edit in the script. Productivity status says the stale mark was "not looked at on screen". Reproduce by hand before fixing; if real, it is an `Engine::problem_key` or `push_edits` ordering bug, not a design one.

**I-10 (Low). Control height and secondary type are just under the desktop default.** `interact_size.y = 26` (default 28, minimum 20) and 43 places at 12.5 and 4 at 12 (default 13, minimum 10) in `kit.rs`, `chip`, `label_caps`, tab bodies; `Body` 14 is fine. Nothing is below a minimum; it is under the default, and the tab-strip "x" (C4) is the one real miss. *Why:* `typography.md`: macOS default 13 pt, minimum 10 pt; `accessibility.md`: control 28 x 28 default, 20 x 20 minimum. *Fix:* `Small` 12 -> 13, `interact_size.y` 28, one scale (13 / 14 / 17 / 22 / title 28) instead of the nine sizes in use (10.5, 11, 11.5, 12, 12.5, 13, 14, 16, 17, 18, 20). Text scaling by the user (Ctrl +/-, `ctx.set_zoom_factor`) is the real answer to "text scales with the setting"; it is not wired (plan U2).

**I-11 (Medium). Dashboard "Earlier" rows are mis-aligned.** In the first run (real workspace, not kept) the two rows had centred titles and summaries instead of one left-aligned column; `docs/mockups/01-dashboard.png` shows the intended grid. *Fix:* `egui_extras::TableBuilder` or fixed-width columns in `dashboard.rs`.

**I-12 (Low). Idle chip and tagline.** "Idle" with a grey dot is always on screen; the dashboard ends with the tagline "creative chaos."; the dashboard shows the logo (tab strip) and the large italic title. *Why:* `branding.md`: "Ensure branding always defers to content ... Resist the temptation to display your logo throughout". *Fix:* the chip appears only when something runs or waits; keep the tagline for the About sheet.

**I-13 (Low). A sheet has no visible way out.** The sheets close on Esc and a click on the backdrop but show no close control (`modality.md`: "obvious way out"). Add a "Close" at the sheet's top right.

**I-14 (Low). Continuous animation ignores reduced motion.** `kit::bar(None)` slides a block with `request_repaint()` for as long as it is shown. egui has no OS reduced-motion query; add a Settings switch "Reduce motion" that makes the unknown-length bar a static striped bar and removes the sheet fade (`accessibility.md › Cognitive`: "reduce automatic and repetitive animations").

## 4. The colour question: why it feels stressful, measured

Brand: plum `#820C4F` (OKLCH L 0.40, C 0.158, h 353), amber `#FFBD59` (L 0.84, C 0.138, h 75), shadow `#531F55`, ink `#24102A`. Not the problem in themselves. Stress comes from six things in the way they are applied:

| # | Cause | Measured | Where |
|---|---|---|---|
| 1 | **Large saturated fills in dark mode.** The dark accent `#D96AA4` is L 0.672, C 0.153 on a page of L 0.182: a luminous, high-chroma magenta on near-black violet. Filled pills are 4 per view. | 4 filled accent objects in `dark-build.png`; contrast with white text 3.21:1 | `kit::primary`, `kit::pill`, tab strip |
| 2 | **Tinted neutrals everywhere.** Dark page C 0.029, card C 0.038, line C 0.045 (hue 314): every surface leans violet, so brand colour never rests. | Light neutrals are C 0.002 | `theme::page/card/line` |
| 3 | **One colour, many meanings.** Plum = brand, primary, selected, link, title text, focus. Green = saved, success, legacy Compute. Amber = running, stale, warning. Red = danger, error, boundary. Blue = dirty Build mesh and domain 1. | see I-5 | 56 literals |
| 4 | **Low-contrast edges and high-contrast fills.** Boundaries at 1.3:1 beside fills at 5:1 make the fills look like the only structure; the eye is pulled to colour, not to content. | C2 | `apply_mode` |
| 5 | **Black canvas in a light shell, rainbow in the canvas.** L 0.97 beside L 0.22; 12 categorical hues at half alpha. | I-3, I-4 | `app.rs`, `colormap.rs` |
| 6 | **Colour with no redundancy**, so users try to read it harder (dots, chips). | C3 | tab strip, status chip |

**The way to use the palette** (decision text in ADR 90; tokens in `docs/ux/tokens.json`; checked values in section 7):

1. **Brand rests in four places only:** the logo, the italic titles, the single primary action, and a one-step tint for the selected item. Everything else is neutral.
2. **Neutrals lose their tint:** dark page `#141218` (C 0.012), card `#1D1A22` (C 0.016), raised `#262230`; light page `#F8F6F8`. The violet cast is kept in the *text* (`#ECE7EF`) so the app still reads as plum, not grey.
3. **Selected means tint + bar, not fill.** A selected tab, pill or experiment gets `tint` (`#33202F` dark, `#F1E4EC` light) and a 2 px accent bar or underline. Only Compute is a filled accent. At most **one filled accent object per view**.
4. **Dark accent is lighter and less saturated:** `#E58BB8` (L 0.744, C 0.122), with ink text on fills (7.49:1).
5. **Amber means time and attention** (running, queued, stale) and nothing else; it is a fill with ink text, or, on light, `#6B4300` text on `#FFF1DC`. Plum never means a state; green and red never decorate.
6. **Every state has a shape and a word** (saved, running, draft, stale, failed): glyph + label next to the colour.
7. **Data colour lives in plots and the canvas only,** from the validated sets (section I-4); sequential default stays Viridis; Turbo stays available but is labelled as rainbow; add Cividis for colour-vision deficiency. The canvas follows the theme by default.
8. **Calm mode:** Settings "Contrast: Standard / High". High raises control edges to 4.5:1 and drops the tint fills for outlines. (egui has no OS Increase Contrast query; this is a setting, not detection.)

## 5. Craft notes

* **Point of view (strong):** the italic serif title over a neutral Inter UI with JetBrains Mono for numbers is specific to the subject, and the three-font split follows the skill's rule (system-style type for controls, personality in display text and content). Keep.
* **Signature:** the italic title. Spend boldness there, not on filled pink chips.
* **Remove one accessory:** the "Idle" chip, the tagline on the dashboard, and the second ("Dark") theme button.
* **Structure encodes information:** the numbered "1 · Concept / 2 · Build / 3 · Analysis" is a real sequence (the one place numbering is earned). Keep.
* **Data typography:** mono numbers with four fixed decimals (`0.9447`, `500.0000`) read as raw output; use significant figures with SI prefix and unit.

## 6. What works (keep)

* Concept / Build / Analysis with a one-line purpose each ("How is it computed? Mesh, boundaries, study and solver.") and Compute as the only primary at the trailing end (`toolbars.md › Best practices`: "Only specify one primary action, and put it on the trailing side").
* Autosave stated in words ("Saved automatically · just now"; "Draft: nothing is saved until you change something"); nothing is deleted by the dashboard rules.
* Ctrl+K palette with sections ("Do something", "Browse what nyx can do"); keyboard shortcuts printed on the buttons.
* Stale handling exists at all (banner and "Solve again") and the Viridis default.
* Cards, spacing and the 8/12/16 radii are consistent; the Share sheet explains each level and gives sizes.
* The whole shell is drawn from tokens in one file; that is what makes the fixes cheap.

## 7. Numbers (reproducible)

Contrast (WCAG 2.x relative luminance), current vs proposed. Script: `docs/ux/contrast.py` (pure Python, no dependencies).

| Pair | Now | Proposed |
|---|---|---|
| Dark: white on accent fill | **3.21** | ink on `#E58BB8` 7.49 |
| Dark: accent text on page | 5.87 | 7.72 |
| Dark: muted on card | 6.62 | 6.83 |
| Dark: red on page | **2.98** | 6.94 |
| Dark: green on page | 3.74 | 8.33 |
| Dark: control edge on card | **1.34** | 3.10 |
| Light: muted on page | 5.33 | 7.03 |
| Light: amber on white (as foreground) | **1.66** | amber ink on tint 7.77 |
| Light: control edge on page | **1.17** | 3.50 |
| Light: plum text on page | 9.21 | 9.91 |

Palette validator (dataviz, OKLab dE x 100, CVD = Machado 2009 severity 1.0):

| Set | Surface | Result |
|---|---|---|
| Current domain palette, first 8 | dark canvas `#181A1E` | lightness FAIL, chroma FAIL, CVD FAIL (worst adjacent 4.3, deutan) |
| Proposed, dark steps | `#181A1E` | all PASS (worst adjacent CVD 8.4, normal-vision 19.3) |
| Proposed, light steps | `#F4F2F4` | all PASS, contrast WARN for 4 slots (needs the numbers, which the viewport draws) |

## 8. Platform notes

* nyx is Windows desktop in Rust/egui. Apple's pages are applied as principles; where Windows has its own convention (title bar, Settings in an app menu is a macOS rule, Windows puts it in a menu or gear), the plan says which it follows.
* **Liquid Glass / blur is not recommended** and egui cannot draw it. `liquid-glass.md` limits it to a floating functional layer anyway; nyx has no such layer to restyle. Do not mimic it.
* egui constraints that shape the fixes: no native menu bar (use `muda` or an in-window strip), no OS Increase Contrast / Reduce Motion query (settings instead), AccessKit exists (`eframe` feature `accesskit`) but custom widgets must set their role.
