# nyx window: reset plan (2026-11-24)

Source: the dist build `nyx-main-9480093a` driven with `nyx --script` (shots in `runs/ux/`, git-ignored) plus a read of `crates/app/src`. Verified = seen in a shot or in code. Hypothesis = not yet measured.

The pattern behind most items: earlier passes (`UX_REVIEW.md`, `UX_ACTION_PLAN.md`) added structure and decoration (italic serif titles, a section strip, sub-pills, sheets, banners) and never removed any. This plan mostly removes.

## 1. Console window opens with the exe (verified)
Cause: no `windows_subsystem = "windows"` anywhere in `crates/app` (grep empty), so Windows gives `nyx.exe` a console. Confirmed: launching it creates a `WindowsTerminal` window titled with the exe path. The same exe is also the CLI (`nyx api`, `check`, `run`, `--batch`, `material`, `queue`, `--job`), which is why nobody set the GUI subsystem.
Fix, option A (recommended): a separate GUI-only binary. Keep `nyx.exe` as the console program (CLI, `--batch`, `--job` worker, scripts, MCP all keep exit codes and piped stdout) and add `nyxw.exe`, same code with `windows_subsystem = "windows"`, which is what the shortcut and the installer point at.
Option B: `#![windows_subsystem = "windows"]` on the one exe plus `AttachConsole(ATTACH_PARENT_PROCESS)` in the CLI branches. Cost: PowerShell's `&` returns at once and `$LASTEXITCODE` is lost for `check`, `run`, `api`, `--batch`, and for any test or tool that waits on `nyx.exe` (`tools/`, `crates/app/tests`, shot scripts); the `--job` worker must be re-checked over a piped stdout.
Option C: keep the console subsystem and call `FreeConsole` at start when `GetConsoleProcessList` says the exe is alone in its console (Explorer launch): a brief flash, no other change.
Workers already use `CREATE_NO_WINDOW` (`worker.rs:366`). Test: launch from Explorer, no console; `nyx api list` from PowerShell still prints and sets its exit code.

## 2. Script editor as a side panel, not a modal (verified)
Cause: `sheets.rs:42` opens `kit::modal(ctx, "nys_script", 920.0, ...)`, covering the model; the editor box is capped at 280 px.
Fix: a right-hand dock (or a bottom dock) in the workflow view, toggled by Ctrl+Shift+R, resizable, sharing the area with the inspector (a tab strip "Inspector | Script" on the right side). Output and check messages go under the editor. The modal code in `sheet_script` becomes the dock body.

## 3. Analysis (verified, the worst one)
What the shot shows (`runs/ux/a.png`): a wide run table whose paint overruns the right panel (panel text is cut at the left edge: "IN A FUNCTION...", "etric vs parameter"); a chart of `eigenvalue_im` against run for one run (one dot, meaningless default); four function buttons; a mode list buried at the bottom; the field of the solve is not on the page at all (it is reached only by clicking a row into `Sub::Field`).
The reference is the sibling project `C:\Users\robin-adm\Documents\opomech\simgui` (its README: "in the style of COMSOL"; shots in `simgui/shots/`, e.g. `21_sapphire_mode.png`, `41_pml_mode.png`, `40_pml_materials.png`). Its whole window is one screen:
- menu bar (File, Examples, Python (optional), View, Help); a toolbar with **one Physics dropdown**, a green Compute, and `Show: mesh | result`;
- left: five tabs (Geometry, Materials, Physics, Study, Results); right: the viewport, always visible, with its colour bar;
- Results tab: a **mode slider with < >** and the label of the mode (`mode 1: neff = 1.813075 -2.223e-6i`), `Field` combo, `Colours` combo, auto range, mesh edges, Fit view, then a plain **Summary** list (`neff 1`, `Q 1 = Re/(2|Im|) of neff`, `f1 (Hz)`, ...) and Export CSV;
- under the viewport a strip plot of the solutions against mode number;
- a **status bar** with the mesh facts (`generated geometry | 2D | 729 points 1280 cells | P2 dofs: 2737`) and the last message. Also the answer to item 8.
- PML (`40_pml_materials.png`): drawn as normal regions, with the PML parameters (direction, strength, order) inside the domain's own card in the Materials tab, not as a separate condition with its own colour.
Fix: rebuild Analysis on that. Order:
1. Analysis (and Results) opens on the latest run: the field in the viewport (the same viewport and colour bar Build shows after Compute), the mode slider with < >, Field and Colours combos, the Summary list on the left. The mode strip plot goes under the viewport.
2. "Runs" becomes a secondary view (a tab inside Analysis, opened only when there is more than one run): compact table, columns chosen from the physics (not every metric), clipped to its panel.
3. Plots and "functions on these runs" sit behind the Runs view and need 2 or more runs; default axes come from the sweep parameter, never `eigenvalue_im` against `run`.
4. Nothing is removed from the engine commands: a machine reads the same data; only the default presentation changes.
Hypothesis for the clipping bug (not tested): the central panel is added after the side panel (`wfview.rs:770` then `:861`) and the table does not clip to its rect.
The old simgui is a separate project (no `simgui` in the nyx2 git history), so what its layout did was read from its screenshots, not from code that I ran.

## 4. Italic serif "Workflows" (verified)
Cause: `theme::title(size)` (EB Garamond italic, accent colour) used for 11 headings: the Dashboard, the workflow title bar, every sheet, Compare, field titles.
Fix: drop the italic face and the accent colour from all of them: sans (Inter) semi-bold, text colour. Remove the bundled EB Garamond (and `OFL-EBGaramond.txt`) once nothing uses it. Rename "Workflows" to what it holds ("Projects" or "Recent") and drop the sub-copy ("Everything you touched in the last 3 days...").

## 5. Navigation depth (verified)
Levels now: workflow tabs / Concept-Build-Analysis / Geometry-Mesh preview / Boundaries-Study-Sweep / (Concept only) Geometry-Materials / experiments list. Concept and Build overlap (the Concept inspector already has "Study & solver →"; Build shows BCs beside the mesh). simgui had one row of tabs and one viewport.
Fix: copy simgui's shape, two levels.
- Left panel tabs: Geometry, Materials, Physics, Study, Results (Sweep inside Study). The viewport on the right is always there, with `Show: mesh | result` in the toolbar.
- The physics are one dropdown in the toolbar (not pills in the inspector), so adding a physics does not rearrange the panel.
- Drop the "1 Concept / 2 Build / 3 Analysis" strip and the sub-pills (`Geometry | Mesh preview`).
- Experiments list and Notes move into a collapsible left section or the workflow menu; Versions, Share, More into one menu.
Menu items `Alt+1..6` keep working by mapping to the new tabs.

## 6. Lag when switching physics (not reproduced)
Frame log (`NYX_FRAME_LOG`) while scripting `physics sms; physics ems` on the silicon-on-sapphire example: worst UI frame 11 ms, nothing slow, so the window does not compute on a switch.
Suspect, not confirmed: `Project.active` is serialised into the spec (`project.rs:148`), so a physics switch changes the spec hash, and `flush()` (`shell.rs:383`) would then commit a version to disk and `refresh_entries()` would reread the workspace. In a scratch-workspace run, four switches left the draft unsaved ("Versions v0"), so either `flush` did not run in that path or a draft is exempt; I did not see a version written by a switch. It needs a real pill click on a saved (non-draft) workflow with the frame log on.
Plan: (a) measure that case first; (b) if confirmed, keep "which physics is shown" in the window state, not in the saved project; (c) the simgui single dropdown removes the switching UI cost anyway. If a lag remains it may be machine load: this session had 2 `nyx.exe` and 4 `nyx-mcp.exe` processes running; close stale ones before judging.

## 7. PML colours (verified)
`runs/ux/pml_b.png`: air blue, slab cyan, `pml_left` red, `pml_right` yellow. Cause: `regions::assign` gives each region a slot by material of the **active physics**, else by its own name (`app.rs:region_paint`), and PML regions have no material. So each PML region takes its own hue, and with several physics the key map changes with the active physics, so regions repaint when you switch.
Fix: colour by role, independent of the active physics. Regions that any physics marks as PML (condition kind PML, its `on` selector) share one hue; the fill opacity (or a hatch density) follows the PML `strength`; the legend shows one row "PML" with the strength range. Materials keep the validated 8 slots, keyed by the material of the first physics that assigns one, so a switch never repaints.

## 8. "Mesh 2D · 7188 cells..." strip (verified)
`wf_strip` bottom panel (`wfview.rs:477`) shows mesh counts plus the sentence "Compute solves on exactly this mesh; it is rebuilt from the recipe, never stored."; in Concept it is a different box ("WHAT SEARCH WILL SEE").
Fix: remove the bottom strip. Mesh facts go to a one-line status bar at the bottom, as in simgui (`2D | 729 points 1280 cells | P2 dofs: 2737`, plus the last message such as "Mesh built: 1280 cells, 176 boundary facets (0.02 s)"); a notes section, if you want one, sits in the left panel. The explanatory sentence moves to a tooltip.

## 9. Yellow "Stale: ..." tag (verified)
`app.rs:1497`: an amber filled box painted over the top-left of the plot whenever the result no longer matches the setup; in the shots it also sits on top of the "Exploratory" chip.
Fix: a small amber dot with the word "stale" in the viewport's own toolbar row (next to the field name), with the reason in a tooltip; a stale result is drawn at reduced opacity instead of a banner. No overlay on the data.

## Found while navigating
- Tab strip clips the first tab ("m_modes", "00nm_modes" cut at the left edge) and keeps 5 to 6 tabs open; truncates titles ("300 nm SOI wire with..."). Fix: a scrolling strip with a fixed minimum width per tab, close on middle click, a "Recent" menu instead of a tab per workflow.
- The workflow title is 46 px italic and is cut off ("sloped left wall: i" then "+ tag" overlapping). Fix: one line, ellipsis, rename on click.
- Header carries too much: "Draft: nothing is saved until you change something", "More", "Versions v0", "Share", "Script", "Compute F5", "1 run · about 0.1 s each · no mesh stored". Fix: save state as a small dot in the tab; Versions, Share, More into one menu.
- Everything selected is the same pink (tab, section, sub-pill, physics pill, Compute): five filled pink objects at once (see `UX_REVIEW.md` I-2, never fixed). Fix: one accent, used only for Compute and the current selection in the rail.
- Inspector text noise: "Takes from other physics of the project: temperature (ems.temperature). Offers to others: optical_field." and the paragraph under Boundary conditions. Fix: show only when it is non-empty and actionable.
- "remove this physics from the project" sits beside the add button at the same weight as "add physics": a destructive action next to a common one. Fix: into a menu on the physics pill.
- Sheets fade in over the page with the page showing through (known, from the previous review): not a bug.
- Not looked at: running state, light theme, 150 percent scaling, 3D results. Interactive clicking was not possible: computer-use cannot see nyx (it is not a Start-menu app), so I used `--script` shots.

## Order of work
1. Console (1) and `active` out of the saved spec (6): small, independent, certain.
2. PML and region colouring (7): `regions.rs`, `app.rs:region_paint`, one test per rule.
3. Chrome removal: title face (4), bottom strip (8), stale tag (9), header, tab strip. Mechanical.
4. Script dock (2).
5. Navigation rework (5) and Analysis (3) together: they decide the same layout. Do this after 1 to 4 so the screens being rearranged are already quiet.
Each step: shots before and after with `tools/ux-shots.sh`, `cargo test --release -p nyx` for the app, `cargo clippy` clean, one jj bookmark per step.

Process notes: this plan reverses parts of ADR 88 and `UX_ACTION_PLAN.md` (colour roles, section strip), so a new ADR is needed. The shots cited from `runs/ux/` are git-ignored; copy them to `docs/ux/shots/` when the plan is committed. Nothing in this file is committed yet.

## Decisions I need from Paul
- Console: separate `nyxw.exe` (A, recommended), one exe with `AttachConsole` (B), or `FreeConsole` flash (C)?
- Script dock: right or bottom?
- Navigation and Analysis: rebuild on the simgui layout (one row of tabs, always-on viewport, mode slider and Summary), as written above?
- Keep experiments, versions and the workflow dashboard at all, or cut to one project file at a time as in simgui? (This is the biggest simplification and the biggest product decision.)
- Rename "Workflows": "Projects" or something else?
