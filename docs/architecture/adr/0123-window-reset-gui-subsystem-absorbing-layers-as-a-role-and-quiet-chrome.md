# 123. The window reset: GUI subsystem, absorbing layers as a role, plain headings and quiet bars

Date: 2026-11-24

## Status

Accepted by Paul on 2026-11-24 ("the points look good", reading `docs/ux/UX_RESET_PLAN.md`), implemented in steps on the bookmarks `ux-reset-<n>-<name>`. This record covers what is built (steps 1 to 3). The larger moves of the plan (one row of tabs and an always-on viewport as in the sibling project `simgui`, the Analysis page opened on the field and the mode slider, the script as a dock, the physics as one dropdown) are decisions still to be taken one by one and get their own records. It amends ADR 88 (italic serif titles, the filled amber stale band, "one region, one slot") and ADR 122 (the executable `nyx.exe`).

## Context

Using the window after the UX line (ADR 88, 92, 93) showed a pile of decoration and structure that makes the work harder, not easier (`docs/ux/UX_RESET_PLAN.md`, shots read on 2026-11-24 from the build of `9480093a`):

* Every start of `nyx.exe` from Explorer opened a console window that stayed open: the program had no `windows_subsystem`, because the same executable is also the command line (`api`, `check`, `run`, `material`, `queue`, `--batch`, `--job`).
* An absorbing layer was coloured like a material: with no material of its own each PML region took a slot of its own, so a slab with two layers showed four hues; and the plan was made from the materials of the **active** physics only, so switching the physics on show repainted the model.
* Headings were a 26 to 46 px italic serif (EB Garamond) in the accent colour: the look of generated pages, and a size that cut the workflow title off.
* The bottom of Build was a card with a sentence ("Compute solves on exactly this mesh...") under the viewport, and a filled amber band painted over the data said when a result was stale.

## Decision

* **`nyx.exe` is a program of the GUI subsystem on Windows** (`#![windows_subsystem = "windows"]`): no console window opens with the window. The command-line modes (`api`, `check`, `run`, `material`, `queue`, `--batch`, `--version`, `--help`) call `console::attach_parent()` first: it joins the console of the terminal that started the program and points the standard handles at it **only when the caller did not redirect them**. A caller that pipes or redirects (the MCP server, the tests of `crates/app/tests/cli.rs`, `| Out-String`) gets exactly what it got before. The worker modes (`--job`, `--worker`) and the window attach nothing. Cost, accepted: in an interactive PowerShell or cmd, the prompt returns at once while `nyx check` still prints, and `$LASTEXITCODE` is not set by a GUI-subsystem program; scripts that need the code use `Start-Process -Wait -PassThru`, a pipe, or `cmd /c`. A second executable (`nyxw.exe`) was considered and rejected for now: ADR 122 names `nyx.exe` in the PATH entry, shortcuts, file associations, uninstall, the packager and the process lists (about twenty places), and one executable keeps the five-executable layout of ADR 122.
* **An absorbing layer is a role, not a material.** Every region that any physics of the project marks with a `pml` condition shares one legend row "PML", drawn hatched in the one neutral grey (`REGION_GREY`) graded by the strength of its layer over the strongest layer of the model (55 percent grey over the canvas for the weakest, the full grey for the strongest; opaque, so the ADR 88 floors still hold, and the order of the grades is the order of lightness, which a colour-blind reader keeps). The row is the last, so the materials keep the slots they had without layers.
* **The colour plan reads every physics of the project, never only the active one.** The first physics (in project order) that gives a region a material names its colour; a region is a PML if any physics says so, with the largest strength. Switching the physics on show does not repaint the model.
* **Headings are plain.** `theme::title` is the interface face (Inter) at 62 percent of the size asked for, in the text colour, not the accent: the italic face is no longer loaded or drawn. The accent colour stays for the action, the selection and the focus (ADR 88). The bundled `EBGaramond-Italic[wght].ttf` and its notice (`OFL-EBGaramond.txt`) are now unused; removing them touches the packager and the installer (ADR 122) and is left to a record of its own.
* **Bars, not cards, under the viewport.** The strip under the viewport is a one-line status bar (muted, 12.5 px, a line above it): the facts of the mesh, with what they mean in a tooltip. The stale state of a result on the plot is a small chip (an amber ring and the word "Stale") in the style of the trust chip, with the reason as its tooltip, not a filled band.

## Consequences

* Launching from Explorer, a shortcut or the installer opens no console. `nyx api list` from a terminal still prints. A process that starts `nyx.exe` and waits for it (`Command::output`, `Start-Process -Wait`, Git Bash) still waits and still gets the exit code.
* The legend of a model with layers has one row less than regions, and a layer never shares a hue with a material.
* `crates/app/tests/colour_literals.rs` is unchanged: no colour literal was added.
* No engine command, schema, contract exemplar or file format changed (ADR 19, 63): `Project` and `ProblemSpec` are as before, so `docs/architecture/contracts/` and `VERSIONS.md` have no new line.
* Not decided here, and so not changed: whether choosing a physics writes a version (`Project.active` is in the saved spec; not reproduced as a lag, see the plan, section 6), the tab strip that scrolls its first tab half out of view, the header, and everything about navigation and Analysis.

## Code

`crates/app/src/{console,main,regions,theme,app,viewport,wfview,dashboard,sheets}.rs`, `docs/ux/UX_RESET_PLAN.md`. Tests: `console::tests`, `regions::tests` (layers share a row, the grade follows the strength, the materials keep their slots).
