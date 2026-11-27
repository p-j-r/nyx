# nyx plotting plan (2026-11-20)

Design note, nothing built. It answers: "`plot_sos_mech_modes.py` was written to plot nyx results; it should be extended to all use cases and be part of nyx." Labels: **read** (I read it in the code or the script), **inferred**, **unchecked**.

## 1. Short answer

1. **The script cannot be merged as it is.** `*.py` is git-ignored on purpose and nothing Python ships (CLAUDE.md). What it is good for: a **specification** of what a mode figure must show, and a **visual oracle** to compare the Rust output against. The feature itself is Rust.
2. **nyx already draws part of this** (read): field colours (six maps), the deformed shape with a scale slider, mesh edges, the Poincare sphere as SVG, CSV export, window screenshots in scripted runs (`nyx --script "...; shot out.png"`). It does **not** have: a headless figure an engine command or an AI can ask for, arrows, a gallery of many modes, phase-aligned complex modes, a cross-section of a 3D periodic mesh, spectrum or dispersion plots, the optical-next-to-mechanical layout.
3. **The script found a trap in the results file** (section 3). Every external reader has to rediscover a node permutation. That is a defect of "portable and shareable", fix it first and it is cheap.
4. **Build it as three layers**: analysis functions (phase, normalisation, slice) in `nyx-analysis`; a figure description with an SVG writer in a new crate; an engine command `figure` that the window, scripts and MCP all use. The window draws the same description on screen.
5. **First slice** (section 7): the layout fix, plus one engine command that writes an SVG mode gallery.

## 2. What the script does, and where each part goes in nyx

| Script part | What it does | Where it goes | Hard-coded in the script (must become data) |
|---|---|---|---|
| `load_sms_mode_2d` | reorders P2 cell nodes, cuts the 3D mesh at `y = 0`, splits each 6-node triangle into 4 linear ones, picks the global phase that makes the in-plane motion real, fixes the sign, normalises by the peak inside a "core", takes Im(uy) as the quadrature part | analysis (`nyx-analysis`) and mesh slice (`femcore`) | the `y` cut (`abs(y) < 1e-12`), the core box `abs(x) <= 0.45, z >= -0.05`, the component names `00_ux`, `01_uy`, `02_uz` |
| `load_ems_mode_2d` | reads precomputed `viz/ems/sol_N/_E_`, `Re_Ex` ... | figure layer reading the result, not the file | a different route from sms (sms goes through raw `cell_fields`) |
| `plot_deformed_mode_panel` | colours the field on the deformed mesh (Gouraud), overlays region outline, draws arrows of in-plane motion on a regular grid inside the solid | figure layers: `field`, `outline`, `arrows` | the outline is a hand-drawn rectangle or rib polygon; "in solid" tests are hand-written inequalities; the arrow scale and grid |
| `make_strip_3d_overview` | grid of all bound modes, each titled, plus a summary scatter of frequency against polarisation share | `modes_gallery` figure | frequencies, phase velocities, "99.95% Si" and "53.5% uz" typed in as strings (they are the mode scalars of ADR 80) |
| `make_strip_3d_components` | per mode: `|u|` with arrows, Re ux, Re uz, Im uy | `components` figure | which modes, which columns |
| `make_strip_q0_figure` | q = 0 resonances with a wide window that shows radiation into the substrate, Q from the complex frequency | `field` figure with a window option, Q from the eigenvalue | window boxes, contrast boost x2.5 |
| `make_optical_vs_mechanical_figure` | optical `|E|` and mechanical `|u|` in the same geometry, side by side | `compare` figure over two results | file paths of two result files from two workspaces |
| `export_interactive_widget_data` | samples modes on a 56 x 32 grid, writes JSON for an animation of `Re(u e^(-i w t))` in a browser | a phase slider or play button in the window viewport, no JSON | everything |
| paths (`WS_DIR`, `OPT_DIR`, `ART_DIR`) | a copy into another tool's folder | not needed | all |

## 3. The trap in the results file

**Read** in `h5io.rs` and `xdmf.rs`:

- `results/<p>/solution_k/cell_fields/<name>/{re,im}` are stored per cell node in **nyx's own local order** (`v0 v1 v2 e01 e02 e12` for a P2 triangle, `v0..v3 e01 e02 e03 e12 e13 e23` for a P2 tetrahedron). The layout comment says only "f64 per cell node; attrs: unit, nloc".
- `viz/<p>/cells` is in **VTK/XDMF order** (the permutations `[0,1,2,3,5,4]` and `[0,1,2,3,4,7,5,6,8,9]` of `xdmf.rs::topology`). These are exactly the `PERM_TRI6` and `PERM_TET10` the script needed, and it verified them to 2.2e-16 by hand.
- Nothing in the file says which order `cell_fields` uses. A reader who combines the two (as the script does, to colour a mode with exact complex values on the viz mesh) must find the permutation by experiment.
- The script reads sms from `cell_fields` and ems from the nodal `viz` fields. **Read** (`nyx-sms/src/solid.rs`, `complex_mode_fields`): the nodal sms fields are `|u|` and `Re u{x,y,z}` only, so the imaginary part exists only in `cell_fields`. That is why the script had to go there.
- The phase is **not arbitrary** in the file: the solver rotates each complex mode so its largest component is real and of unit modulus. The script uses a different convention (maximise the real in-plane motion, then fix the sign). Both are conventions; the two figures can differ by a phase.
- Related, from the format study: `cell_fields` look about five times duplicated for sms (unverified).

Consequence for the window: it shows `Re u` in the solver's convention and cannot show `Im u`, so no quadrature view (`Im uy` in the script's figure 2) and no phase animation without `cell_fields`.

Fix, in order of cost:
1. Add an attribute `local_order` on `cell_fields` (a text such as `nyx: v0..vk, e01 e02 e12 ...` or the permutation to VTK order) and say it in the layout comment. Non-breaking; readers that ignore it are unaffected.
2. Or store `cell_fields` in the VTK order (layout 4, a version bump with a line in `VERSIONS.md`).
3. And store the imaginary part of the nodal complex fields under `viz/<p>/sol_k/` (`Im_<name>`; it is missing today, see above), so a reader needs `cell_fields` only for exactness at element level, not for the second half of a complex mode.

## 4. Architecture

**Layers and crates** (rules: `layering.rs` table, "nothing outside the physics crates names a physics"):

| Layer | Content | Where | Cost |
|---|---|---|---|
| Analysis | phase alignment (the phase that minimises the imaginary energy: `theta = -1/2 arg(sum u.u)`), sign convention, normalisation by a region or by the peak, section of a 3D periodic mesh on a plane, P2 to 4-linear-triangles subdivision | `nyx-analysis` (ADR 56: an analysis is a function of a result and its mesh) | no new crate; it already depends on `nyx-physics` and `femcore` |
| Figure | `Figure` = panels, each a stack of layers (`field`, `outline`, `arrows`, `isolines`, `curve`, `scatter`, `text`) with the data inside, plus an SVG writer | new crate, working name `nyx-figure`, depends on `nyx-physics`, `nyx-analysis`, `femcore` | a row in `layering.rs`, an ADR, C4 and `CODE_TO_ADR` rows |
| Command | `figure` (and `figure_kinds`) in `nyx-engine`; returns a path and a JSON summary; MCP and scripts get it for free (ADR 94: four tools over one kernel, no tool per command) | `nyx-engine` | one command, schema, guide card |
| Window | the viewport draws a `Figure` panel; arrows and the phase slider are layers/settings of the viewport; "Export figure..." | `crates/app` | `viewport.rs` and `scene.rs` extended |

Why not only in `crates/app`: the AI route (a script or MCP asks for a mode gallery and gets a picture or a path) has no window. The Poincare sphere SVG (`sphere_svg`) is the precedent for SVG written without a window.

**Rules for the recipes**, so that no physics is named:
- a *vector field with a displacement role* gives a deformed shape and arrows; a *complex field* gives phase alignment and the Re/Im component views; a *periodic 3D mesh* gives a section on the periodic face. The recipe asks the descriptor for the role (components, unit, complex), not the physics code.
- titles come from the mode scalars (ADR 80: frequency parts, energy shares per component and region) and the study's own labels; none are typed.
- outlines come from the region boundaries of the mesh (facets between different cell labels), not drawn by hand.
- the "core" used to normalise is a region selection (a name or a label), default the whole domain.

**Decisions for Paul** (not mine to take):
1. SVG only first, or PNG too? PNG in pure Rust means a rasteriser (`tiny-skia` or `resvg`), a new dependency; **unchecked** which is acceptable under "pure Rust, no C". SVG renders in every browser and ParaView-free tool and diffs as text.
2. A new crate `nyx-figure`, or a module of `nyx-analysis`? A crate keeps drawing apart from analysis; a module saves an ADR row. I would take the crate.
3. Layout fix as attribute only (section 3, option 1) or a version bump.

## 5. Use cases: physics x study

"App" = what the window does today (read). "Headless" = an engine command could produce it (nothing does today).

| Physics, study | What the figure should show | App today | Missing |
|---|---|---|---|
| `dfs` steady / transient | temperature field, isolines, probe curves over time | colour field, time as solution number | isolines, curves, headless |
| `sms` static | deformed shape, von Mises or `|u|`, arrows | colour field, deformation slider | arrows, headless |
| `sms` eigenfrequency (3D periodic, Bloch) | gallery of bound modes at fixed `q`, section on the periodic face, phase-aligned, labelled with f, vp, energy shares; summary of f against polarisation share | none (one mode at a time, solver's phase, real part only) | slice, Im part, gallery, summary |
| `sms` eigenfrequency at `q = 0` with PML | resonance with a wide window that shows radiation into the substrate; Q from the complex frequency | colour field | window, contrast option, Q label |
| `sms` axisymmetric (r, z) | half-plane mode, mirrored to a full body on request, f and shares | colour field | mirror, headless |
| `ems` modes | `|E|`, Re Ex, Re Ez, Re Ey, effective index, loss in dB/cm, polarisation | colour field, Poincare SVG, CSV | component grid, gallery, headless |
| `ems` cavity (3D) | slice of a cavity mode on a plane, Q | colour field in 3D | slice, headless |
| any, sweeps (`sweep_run`, ADR 88/96) | dispersion curves (f or n_eff against a parameter or `q`), mode tracking by overlap | runs table and plot | the curve figure, tracking (analysis) |
| optical next to mechanical in one geometry | same section, `|E|` and `|u|` side by side | none | `compare` over two results of the same project |
| animation | `Re(u e^(-i w t))` over a period | none | phase slider and play in the viewport; frame export as an SVG sequence later |

Out of scope: **overlap and gain between optical and mechanical modes**. `sbs` was removed (ADR 24). The figures only show the two fields side by side.

## 6. Risks

- **Phase alignment is a convention**, and the solver already applies one (largest component real). The figure layer must either reuse it or say which it used. "Maximise the real in-plane motion" is a choice that works for the guided modes of the script; degenerate or near-degenerate modes make any phase arbitrary within the group, as `polarization` already treats it. The figure must say which convention it used.
- **Slicing a 3D mesh** the way the script does (find 6-node faces on `y = 0`) works for the periodic slab with a mesh node on that plane; a general section needs interpolation. Start with the periodic face.
- **A new SVG renderer is code that can be wrong silently.** Mitigation: the script's six figures as a visual oracle (compare by eye, not by pixel), and tests on the numbers inside the `Figure` (the data), not on the SVG text (ADR 26).
- **Colour maps in the figure** must be the validated ones of `tokens.rs`/`colormap.rs` and readable with colour-vision difference (`cvd.rs`); the script's `magma`, `RdBu_r` are not in nyx.
- **Cost of the layering change** is real: ADR, C4 rebuild, `CODE_TO_ADR`, contract check for the command. Nothing is free.

## 7. Staging

| Step | Content | Done when |
|---|---|---|
| P1 | `local_order` attribute and layout comment (and the nodal `Im_` fields if missing); a test that reads `cell_fields` through the stated order and matches `viz` values | an external reader needs no permutation guess |
| P2 | `nyx-analysis`: phase convention (reuse the solver's or the script's, stated), normalisation, 3D-to-2D section, P2 subdivision; tests on a synthetic mode | numbers of the script's `load_sms_mode_2d` reproduced on one of its h5 files (it exists in `nyx2_benches`) |
| P3 | `nyx-figure` and engine command `figure` with kinds `field` and `modes_gallery`, SVG output | the script's figure 1 reproduced by eye from the h5 of the benches |
| P4 | `components`, `compare`, wide window and Q label | figures 2, 3 and 6 |
| P5 | sweep curves and dispersion with mode tracking | a swept `sms` result gives f(q) |
| P6 | window: arrows layer, phase slider and play, "Export figure..." | the animation of the widget |
| P7 | PNG, if decision 1 allows | |

**First slice buildable now: P1 plus P3 with only `field` and a plain gallery**, about the size of the polarisation SVG. P2 comes with it (phase alignment is needed for the gallery to be stable).

## 8. Not checked

- Whether the other physics with complex fields (ems modes) also keep the imaginary part only in `cell_fields`; only sms was read.
- That the script's six figures are correct physics; the note treats them as the intended picture, not as validated results.
- A pure-Rust PNG route and its licence.
- The next free ADR number on current `main` (the last commit message mentions ADR 98; 70 to 79 and 81 are reserved by other plans).
