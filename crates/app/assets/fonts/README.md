# Fonts shipped with the window

The window looks for its fonts (`crates/app/src/theme.rs`, `first_font`) in `fonts/` next to the executable, then in this folder, then in the system folders. A file that is missing falls back to the next candidate.

## EB Garamond Italic (page and workflow titles)

* File: `EBGaramond-Italic[wght].ttf` (variable font; egui draws its default instance, weight 400).
* Used for: the "Nyx" wordmark, the workflow title and the model heading (`theme::title`).
* Designers: Georg Duffner and Octavio Pardo. Copyright 2017 The EB Garamond Project Authors (https://github.com/octaviopardo/EBGaramond12).
* Licence: **SIL Open Font License, Version 1.1**; the full text is `OFL-EBGaramond.txt` in this folder (copied unchanged from the same source).
* Source: https://github.com/google/fonts/tree/main/ofl/ebgaramond, downloaded 2026-11-17. SHA-256 of the font file: `bba2c4499c93c9612b90b9825d32b07da52fce2fe57562a1eb6b833553f93c4e`.

### What the licence asks of us

* The copyright notice and the licence text travel with every copy of the font: **`OFL-EBGaramond.txt` must be shipped next to the font** in `dist/` and in any installer.
* The licence must be **shown to users**: the application's About / Licences view displays the text of `OFL-EBGaramond.txt` and names the font (plan item U2 in `docs/ux/UX_MERGED_PLAN.md`). Until that view exists, this file and the licence file are the notice.
* The font may be bundled and embedded in the program; it may not be sold on its own, and a modified font must not use the name "EB Garamond" if the licence file declares a Reserved Font Name (this file declares none). We do not modify it.

If the font is ever compiled into the binary (`include_bytes!`), the licence text is compiled in with it (`include_str!`) so the notice cannot be separated from it.
