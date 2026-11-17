# nyx window: merged UX plan

Date: 2026-11-17. Workspace `nyx2-ux-compare`, bookmark `ux-compare` (documentation, one palette-validation script's worth of numbers and one font file; **no Rust source changed**). Decision record: ADR 88 (`docs/architecture/adr/0088-colour-has-roles-and-the-window-follows-the-system.md`, number chosen by Paul). Tokens: `docs/ux/tokens.json` (corrected here). Evidence: `docs/ux/UX_REVIEW.md`.

This file merges two parallel sessions that reviewed the window and drew mockups from the same base (`ux-0-review`):

* **S1** = workspace `nyx2-ux`, bookmarks `ux-0-review`, `ux-0b-mockups`; boards **B1 to B6**, PNGs in `docs/ux/mockups/`.
* **S2** = workspace `nyx2-ux2`, bookmark `ux2-0-mockups`; boards **U-1 to U-6** in the Paper file "nyx - app mockups", page "nyx UX (review mockups)" (the PNGs are not in git).

They are not independent: one review, one plan, one `tokens.json`, one ADR. The comparison is of execution and coverage. **Paul accepted the plan on 2026-11-17**; nothing below has been implemented.

---

## 1. Decisions taken (Paul, 2026-11-17)

| # | Decision | Consequence |
|---|---|---|
| Colour vision | The colours must work as well as possible for people with colour-vision differences | Region palettes replaced and re-validated with `--pairs all`; number chip and legend on every region; states always glyph + word (section 4) |
| D1 | Native menu with `muda` | `menu.rs` builds a native menu from `Engine::catalog()`; the in-window menu strip drawn in U-1, U-2, U-3 and U-6 is **not** built (section 6, U2 and risk R2) |
| D2 | Canvas may default to dark | The viewport and plot canvas is dark in both themes; Settings offers "follows the theme". The dark region palette is therefore the one normally seen, also in a light window |
| D3 | Follow the system theme by default, with easy options | Settings: Theme = System (default) / Light / Dark. The header Light/Dark button is removed |
| D4 | Keep the brand pink | Dark accent stays today's `#D96AA4` (not `#E58BB8`), light accent stays plum `#820C4F`. Only the text on them changes (ink on the dark pink: 5.62:1; white on it fails at 3.21:1) |
| D5 | High contrast: read the Windows setting if easy, and as a Settings toggle everywhere | Toggle on every platform; on Windows, read the system flag at start (section 6, U2). Linux and macOS use the toggle only |
| D7 | Eight solid region colours; beyond eight as in section 5 (confirmed) | Slots 9 to 16 hatched, 17 and up grey with the chip |
| D8 | Queue: cancel, pause, priority | No hold or release in version one |
| Boards | The S2 set (U-1 to U-6) is the spec | With the four things taken from S1 (section 3) and the corrections in section 7 |
| ADR | Number 88 | Done in this workspace; S2's ADR text still says 90 and is renumbered when merged |
| D6 | The stale banner is probably a bug; note it, look later | Section 8: **known problem I-9**, written for the next session |
| Font | Download EB Garamond and display its licence | `crates/app/assets/fonts/` holds the font, `OFL-EBGaramond.txt` and a README; the About/Licences view is plan item U2 |

Nothing is open: Paul confirmed the beyond-eight policy of section 5 on 2026-11-17.

---

## 2. What both sets agree on

* Rating Critical on accessibility (C1 to C4); findings I-1 to I-14 identical.
* 17 colour roles, the five rules (one filled accent per view; amber only for time and attention; glyph and word for every state; no `Color32` literal outside `theme.rs` and `colormap.rs`; neutral chroma at most 0.026).
* System theme, no Light/Dark button; separate `separator` and `control_edge`; a canvas that can follow the theme.
* Phases U0 to U6 and the same order (U1 alone, then U2, U3 and U4 in parallel).
* Verified by re-running the contrast script: every text and edge pair of the old token file passes (the colours that fail are the data colours, section 4). The code citations of the review match the source (`kit.rs:19,47` white on the accent; `shell.rs:755,762` tabs and close as bare labels; `shell.rs:782` the Light/Dark button; `theme.rs:14,91` one colour for separators and control edges; `theme.rs:99` height 26; `app.rs:1240` hard-coded canvas; `viewport.rs:279,282` alpha 0.45 and 0.55 fills).

## 3. What each set contributes

| Topic | S1 (B1 to B6) | S2 (U-1 to U-6) | Merged |
|---|---|---|---|
| Documents (review, plan, ADR text) | Plan wrongly says `selection.stroke` is the focus ring | Correct focus note, "fills are opaque" decided and tabulated | **S2's documents**; S1's PNGs and boards are layered on top |
| Region fills | Plan says "validated 8" but board B3 draws a blended navy against a bright swatch | Opaque fills, numbers in chips, legend that matches (U-3, U-6) | S2 |
| Shell, menu, Dashboard (Running / Pinned / Recent with sparklines, Earlier table) | Dashboard light, Earlier table with a state column (B4) | Dark and light, menu strip, search | S2, with B4's richer Earlier table |
| Queue and sweep | Page with Remove and Cancel (B6) | Popover, stale band, sweep card with the engine's plan (U-3) | S2, **minus hold and release** (D8) |
| Settings | B6 | U-5 | S2 plus the additions in section 7 |
| Build in light | no | U-6 | S2 |
| States sheet: focus ring, selection, disabled | **B2** | no | **S1** |
| Analysis multi-series plot (markers, end labels, ringed selected run) | **B5** | no | **S1** |
| Compute button before and after queueing | **B1** | no | **S1** |
| "Close  Esc" control on sheets | B6 | no | **S1** |

Dropped: S1's `selection.stroke` focus claim; board B3 (superseded by U-3 and U-6); B6's queue page and its "High contrast raises edges to 7:1" (U-5 and the plan say 4.5:1, and 4.5:1 is what `tokens.json` now says).

## 4. What was wrong in both, and the fix

1. **Focus ring.** egui 0.37 draws focus as the `widgets.active` visuals (`egui-0.37.0/src/style.rs:1386`); `selection.stroke` is the text-selection colour. Fix: give `widgets.active` a 2 px accent `bg_stroke`, and paint the ring by hand on `kit::tab` and on clickable cards.
2. **Blended fills fail.** Re-computed with the dataviz validator: dark at alpha 0.45, 0.55 and 0.65 fails the normal-vision floor (11.1 to 14.3, needs 15); light fails at 0.45 and 0.55. **Region fills are opaque** (`fill_opacity: 1.0`), mesh lines on top.
3. **The old data palettes pass only for adjacent pairs.** On a mesh any two regions can touch, and `--pairs all` fails even when opaque: dark `#0CA30C` against `#D95926` under deuteranopia is dE 1.4; `#E66767` against `#D95926` is 7.1 for normal vision; light `#008300` against `#EB6834` is 3.2 (protan). They are replaced:

   | Set | Colours (slot order) | Colour-blind separation, all pairs | Normal vision | Notes |
   |---|---|---|---|---|
   | Dark canvas | `#2B60C2 #109EC6 #C83220 #A68729 #9D71F3 #C67294 #2A7D40 #9D399B` | 6.7 (WARN band) | 15.8 PASS | Two slots 2.9:1 on the canvas |
   | Light canvas | `#B69AF8 #E64F65 #9B084D #2A81C2 #BD9D0D #1F7C13 #5F36A1 #46BEBB` | 8.8 PASS | 18.4 PASS | Three slots under 3:1 on the canvas |

   The dark set is inside the validator's "legal only with a second cue" band, because the dark lightness limits leave too little room for eight colours. The second cue (number chip, legend row, outline for the selected region) is **mandatory**. A longer search or a hand-tuned set may do better. The palettes were found by a random hill-climb and are untested on real colour-blind readers and untested in the app.
4. **Dark slot 6 is almost the brand pink** (dE 4.6, 2.4 for deuteranopes). A selected region is never shown by accent colour alone: a 2 px outline in the text colour plus the highlighted legend row. No single replacement colour kept the other checks passing (searched).
5. **Ctrl +/- and Ctrl 0 already work in egui** (`zoom_with_keyboard` is on by default, `egui-0.37.0/src/memory/mod.rs:336`, `context.rs:2931`). Only persisting `zoom_factor` and a zoom entry in the menu are missing; do not re-bind the keys. `ThemePreference::System` is also egui's default (`memory/theme.rs:75`); `Context::set_visuals_of` exists (`context.rs:2691`). Read in the source, not run.
6. **The workspace table of the plan claims no file overlap; there is some.** The idle chip (review I-12, planned in U3) is in `shell.rs:795`, which U2 owns; U2 and U4 both touch `app.rs`; U1 touches `shell.rs` and `kit.rs`. The idle chip moves into U2.
7. **The mesh-edge colour `from_white_alpha(60)` (`viewport.rs:298`) disappears on a light canvas.** Derive the edge colour from the canvas.
8. Smaller: B5's series lines are 2.46 to 2.8:1 on the light canvas (under 3:1 for graphics); markers plus end labels in the text colour make up for it, as the board says. U-4 names `contrast.py`; the file is `contrast.py.txt` (`*.py` is git-ignored on purpose).

## 5. Colouring when there are more than eight regions or materials (confirmed by Paul, 2026-11-17)

The rule in `tokens.json` is: **never generate a ninth hue**. A ninth hue would not be separable from the other eight for colour-blind readers, and it would break "colour follows the entity". So:

1. **Colour follows the entity, not the rank.** The colour of a region is set by its *material* when it has one (several regions of the same material then share a colour, which is also what a reader wants), otherwise by its region id. Slots are given in order of first appearance in the project, so filtering, hiding or reordering never repaints.
2. **Slots 1 to 8:** the solid colours above.
3. **Slots 9 to 16:** the same eight hues again, drawn with a *hatch* over the solid fill (1 px diagonal lines in the canvas ink colour, 6 px apart). A hatch is a texture of thin lines: it is the one cue that does not depend on colour, so it also separates the second use of a hue from the first for a colour-blind reader. The number chip says which region it is.
4. **Slots 17 and more:** a neutral grey fill with the number chip only. The legend collapses to "N more" and opens the table view (list of regions with their material and colour).
5. A legend row click isolates the region (others dim); the selected region keeps its outline.
6. Settings keeps "Colour regions by: Material / Region" for models where the default is not what the reader wants.

Cost: the viewport draws regions as one coloured triangle mesh, which cannot hatch; slots 9 and up need a second painter pass (a pattern fill per region polygon, or lines clipped to the region). That is the one new drawing feature in this plan (U4). Until it exists, slots 9 and up fall back to point 4 and nothing is wrong, only less colourful.

---

## 6. Ordered Rust action list (egui 0.37, files in `crates/app/src/`)

Gate for every step: `cargo test --release -p nyx`, `cargo clippy --release -p nyx --all-targets --no-deps -- -D warnings`, `bash docs/architecture/check_code_to_adr.sh`; `crates/app` is not rustfmt-formatted (ADR 65). No step needs a golden change.

**U0, baseline (S).** Re-measure the two failing `nyx` window tests (`every_tab_draws_for_the_examples`, `the_window_checks_convergence_on_its_own_thread`: "Dropped TexturesDelta with 1 unapplied deltas") on the current base and record them in `docs/UX_STATUS.md`. Write `tools/ux-shots.sh` (scripted shots in both appearances, scratch `NYX_WORKSPACE`, never the real workspace). To capture a *running* state, hold a job in a window test through the engine's queue instead of waiting on a real solve.

**U1, tokens and accessibility (M; alone first, everything else waits).**
* New `tokens.rs`: pure data from `docs/ux/tokens.json`, no egui types, and the contrast test of the old plan (text >= 7, muted >= 4.5, status >= 4.5, `on_accent` on `accent` >= 4.5, edge >= 3.0). Add the *all-pairs* data test: every pair of region colours separable under protan and deutan at the thresholds the validator uses, and the 15.0 normal-vision floor, so a future palette edit cannot regress.
* `theme.rs`: build both `Visuals` once and install them with `Context::set_visuals_of`; widget `bg_stroke` = `control_edge`, `noninteractive` = `separator`; `widgets.active` gets the 2 px accent stroke (the focus ring). Keep the public functions as thin wrappers so the call sites do not change in one go.
* `kit.rs`: `primary` (line 19) and `pill` (line 47) use `on_accent`; `danger` (line 40) uses `signal(Danger)`; `dot` (line 78) becomes `status(ui, State, label)` that draws a glyph and a word; new `tab()` (28 px high, `WidgetInfo::labeled(WidgetType::Button)`, focus ring, 24 x 24 close target).
* `shell.rs`: the Dashboard button (735) and the tab strip (755, 762) use `kit::tab`.
* A test that fails on any `Color32` literal outside `theme.rs`, `tokens.rs` and `colormap.rs`.

**U2, conventions (M; parallel with U3 and U4 after U1).**
* Remove the Light/Dark button (`shell.rs:782`); `settings.theme` defaults to "system"; Settings keeps System / Light / Dark.
* **Native menu with `muda`** (decision D1): new `menu.rs` built from `Engine::catalog()` (File, Edit, View, Run, Help). The command palette and the menu are two views of the same catalog. Start with a one-day spike on Windows and then on Linux and macOS (risk R2). `muda` is a new dependency: record it in ADR 88 or a new ADR and justify it in the commit.
* Persist `zoom_factor`; add zoom entries to the View menu; do not re-bind the keys.
* The idle chip (`shell.rs:795`) moves here from U3.
* Type scale and `interact_size.y = 28` (replacing 26, `theme.rs:99`).
* Settings additions: **Theme** (System default), **Canvas** (Dark default / Follow theme), **Contrast** (Standard / High; on Windows initialised from the system High Contrast flag read once at start, e.g. `SystemParametersInfoW(SPI_GETHIGHCONTRAST)` behind `cfg(windows)`: add it only if it costs one small function and one dependency already in the tree; otherwise the toggle ships alone), **Reduce motion**, **Colour regions by** (section 5).
* **About / Licences** view that shows the text of `OFL-EBGaramond.txt` and names EB Garamond (`include_str!`), and lists the other font files and their licences. This is what `crates/app/assets/fonts/README.md` promises.

**U3, craft and copy (S to M).** **Logo (Paul, 2026-11-17):** the mark in the top bar looks cropped (`shell.rs:730` draws `docs/brand/nyx-mark.png`, a 286 x 286 square with an offset shadow, at 34 px with `corner_radius(6)`, which clips the corners and the shadow; draw it uncropped, about 36 to 40 px, no clip), and the window and the executable have no icon (`main.rs:77` sets no `with_icon`; add `ViewportBuilder::with_icon` from the same PNG via the `image` crate already in the tree, and a Windows executable icon resource if that costs one build script). Title ellipsis; new `fmt.rs` helpers (`plural`, `si`); descriptor labels in place of raw keys; the Dashboard table as an `egui::Grid` with the state column of B4; a "Close  Esc" control on every sheet; roadmap copy removed from the UI.

**U4, data colour and the canvas (M).**
* `colormap.rs:100` `palette()` returns the validated eight for the canvas actually drawn; fills opaque in `viewport.rs:279,282` (the 3D branch needs the same decision; shading must not change a region's hue); mesh-edge colour derived from the canvas (`viewport.rs:298`).
* `app.rs:1240` takes the canvas from the theme setting (default dark).
* Region number chips (canvas-coloured chip, text colour inside), legend with isolate-on-click, selected outline in the text colour.
* The beyond-eight policy of section 5: material-keyed slots, hatch pass for slots 9 to 16, grey for 17 and more.
* Cividis added; Turbo labelled as rainbow; Viridis stays the default.

**U5, productivity surfaces (L).** Everything the engine already has (verified in the source: `queue_*` commands in `commands.rs:171-175`, owners `gui`/`ai`/`cli` in `jobs.rs:90`, `sweep_run`/`sweep_table`/`sweep_cancel` in `commands.rs:178-185`, the plan of a command list in `commandlist.rs`):
* Busy chip with glyph, word and counts (`shell.rs:795-812`).
* Queue popover through the `queue_*` commands: **Cancel, Pause, priority** (D8); no hold or release. Rows say "asked by an assistant" when the owner is `ai`.
* Running card on the Dashboard; the stale band (after I-9 is reproduced).
* Sweep section where `wfview.rs:585` has a "coming soon" card, driven by `sweep_run`/`sweep_table`/`sweep_cancel` and showing the plan ("21 meshes, 21 solves").
* The Compute button reads "Queue sweep · position n" while busy (B1).
* Gaps, not engine features today: a *time-remaining* estimate ("left about 01:20"; the engine has a memory estimate, `estimate.rs`, not a time one) and the stage wording ("factorising the complex matrix"; check what `job_status` returns). Draw only what the engine reports.

**U6, verification (M).** `tools/ux-shots.sh` as a smoke test; a manual pass with Narrator, Windows High Contrast and 150 and 200 percent scaling (never done in the review); a screenshot of a running job and of the queue panel (never seen on screen); the colour-blind check by viewing the canvas through a protan and a deutan simulation and a person with a colour-vision difference if one can be found.

**Analysis plot (B5).** Marker shapes, end labels and the ringed selected run in the Analysis tab. Window-only; whether the current plot code can take it has not been checked (put it in U4 or U5 once `analysis.rs` has been read).

### Workspaces and file overlap

U1 alone (it touches `theme.rs`, `kit.rs`, `shell.rs`, new `tokens.rs`). Then U2 (`shell.rs`, `menu.rs`, `settings`), U3 (`fmt.rs`, `dashboard.rs`, `sheets.rs`) and U4 (`colormap.rs`, `viewport.rs`, `app.rs`) in three workspaces; **U2 and U4 both edit `app.rs`**: do the `app.rs:1240` canvas line in U2 and let U4 only consume it. U5 after U2 (shell chip, popover) and U4; U6 last. One writer per jj workspace; never run jj in the default workspace.

---

## 7. Boards to correct before they are used as the spec

* U-1, U-2, U-3, U-6: remove the in-window menu strip (the menu is native).
* U-3: the queue popover loses Hold and Release; keep Cancel, Pause and priority.
* All boards: dark accent is `#D96AA4` with ink text on it, not `#E58BB8`.
* U-3 and U-6: swap the region legend and canvas to the new palettes of section 4 and add the number chips to every region.
* U-4: name `contrast.py.txt` and the all-pairs check.
* U-5: add Theme (System), Canvas, Contrast (High, with the Windows note), Reduce motion, Colour regions by, and an About / Licences entry.
* Canvas default is dark in both themes; the light window boards (U-2, U-6) show a dark viewport.

## 8. KNOWN PROBLEM, NOT YET FIXED: stale banner on a fresh result (review I-9)

**Treat this as a bug until proven otherwise. Do not design around it, do not explain it away, and do not ship the busy and stale surfaces (U5) before it is understood.**

* What was seen: example 03 (heat), `compute`, then `view result` showed the banner "STALE: the model changed since this solve" although nothing was edited (`docs/ux/shots/light-heat-result.png` in the review's evidence).
* Why it is suspicious: a result that was just computed must be fresh. If the banner appears here, either `Engine::problem_key` differs between the solve and the check, or `push_edits` runs after the solve and marks it stale. That is an engine-side ordering or identity bug (ADR 64: one persisted hash, the rest are cache keys), not a design one.
* Why it is not yet a confirmed bug: the review's own scripted run may have caused it, and the productivity status says the stale mark was "not looked at on screen". Nobody has reproduced it by hand.
* What to do: reproduce by hand in the running window (compute example 03, open the result, no edit), then in a test: solve a project, then call the same freshness check the window calls, and assert it is fresh. If the test fails, fix the engine, keep the test.
* Paul's instruction (2026-11-17): note it and look at it later. The note is deliberately loud.

## 9. Risks

* R1 **The token swap changes every screen at once.** Do U1 alone, with before and after shots.
* R2 **`muda` on Linux and macOS.** Native menus on Linux go through GTK and need integration with winit/eframe that has not been tried here; on macOS the menu must be created on the main thread. Spike first on all three systems; if Linux fails, ship the egui menu strip there as a fallback. Not verified, from general knowledge of the crate.
* R3 **Hatch** needs a second painter pass and 3D shading under opaque fills needs a look; slots 9 and up stay grey until it exists.
* R4 **Dark region set is at the "needs a second cue" limit** and was found by a simple search; it is not tested on real colour-blind readers. The chips and legend must ship with the colours, not after.
* R5 **Time-remaining estimate and stage wording** are not in the engine; the boards draw them. Show only what the engine reports.
* R6 **Neither session saw a job running on screen or the queue panel.** The running and queued states are drawn from code and earlier mockups.
* R7 **Merge friction between the two branches.** Both edit the same docs. Take S2's documents wholesale (opaque fills, focus note, ADR text, `contrast.py.txt`, renumbered 90 to 88), and layer S1's PNGs, B-boards and `contrast.py` on top. This workspace already carries the renumbering and the corrected `tokens.json`.
* R8 **Comfort claim** (less visual stress) rests on measured chroma and area, not on users.

## 10. What was verified and what was only read

Verified by running: both contrast scripts; the dataviz validator on blended and opaque palettes (adjacent and all pairs); the palette searches; contrast of the brand pink and ink; pixel sampling of the mockup PNGs; Paper node and screenshot reads; the font download (valid TrueType file, licence text read).

Verified by reading: egui 0.37 source (focus, zoom, theme), the app source (`kit.rs`, `shell.rs`, `theme.rs`, `app.rs`, `viewport.rs`, `colormap.rs`, `dashboard.rs`, `wfview.rs`), the engine source (`queue.rs`, `jobs.rs`, `commands.rs`, `commandlist.rs`).

Only read, never run: the application, Ctrl +/- behaviour, any screenshot of the live window of this base, Narrator, High Contrast, scaling, a running-job screenshot, the stale banner, `cargo test` and `clippy`, and `muda` on any platform. The 8-colour palettes are not tested in the app or by colour-blind readers.

## 11. Fonts and licences

EB Garamond Italic is used for the "Nyx" wordmark, the workflow title and the model heading (`theme.rs:152`; `shell.rs:732`, `wfview.rs:272,1042`); it falls back to Georgia Italic. The font file, its licence (SIL Open Font License 1.1, copyright 2017 The EB Garamond Project Authors) and a README with the source, SHA-256 and obligations are in `crates/app/assets/fonts/`. The loader already reads that folder (`theme.rs`, `first_font`). The licence must be shipped next to the font and shown in About / Licences (U2). Inter, JetBrains Mono and the system fonts the loader also tries have their own licences; list them in the same view when they are bundled.
