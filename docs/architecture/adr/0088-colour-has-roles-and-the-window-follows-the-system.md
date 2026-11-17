# 88. Colour has roles, and the window follows the system appearance

Date: 2026-11-17

## Status

Proposed (not yet accepted by Paul; no code changed). **The number 88 was chosen by Paul on 2026-11-17**: it is the number the productivity line (84 to 89) left unused. Evidence: `docs/ux/UX_REVIEW.md`; plan: `docs/ux/UX_ACTION_PLAN.md`; tokens: `docs/ux/tokens.json`, checked by `docs/ux/contrast.py`.

## Context

The window draws every colour from `crates/app/src/theme.rs`, which is good, and then breaks its own rule in 56 places (`Color32` literals in `app.rs`, `tabs.rs`, `sheets.rs` and others). A review of `productivity-4-memory` against Apple's HIG (principles, read from the skill's reference pages) and the data-visualisation checks found:

* white on the dark accent is 3.21:1 (needs 4.5), the dark red is 2.98:1, control edges are 1.17 to 1.44:1 against their surface;
* the same colour carries several meanings (plum is brand, primary, selected and link; amber is running, stale and warning; there are two Compute buttons in two greens and plums), and up to four filled accent objects share one dark view, which is what makes the palette tiring;
* the top bar has its own Light/Dark button although Settings already has Follow system (HIG `dark-mode.md`: avoid an app-specific appearance setting);
* the viewport is always near-black, and its 12-colour domain palette fails the validator (adjacent colour-blind separation 4.3, floor 6).

## Decision

* **Colour is a set of roles, not values.** `theme.rs` owns one `Tokens` value per appearance (the table in `docs/ux/tokens.json`): page, card, raised, separator, control edge, text, muted, accent, on-accent, tint, activity (+ text and tint), success, danger, canvas. Widgets ask for a role; nothing else in the app crate writes a colour. A test fails on a `Color32` literal outside `theme.rs` and `colormap.rs`.
* **One colour, one meaning.** Accent = the primary action, links, focus and the selected-item bar. Activity (amber) = time and attention: running, queued, stale. Success and danger = outcomes. Data colours = plots and the canvas only.
* **Quiet by default.** At most one filled accent object per view (Compute). Selection is a tint with a 2 px accent bar. Neutral surfaces carry OKLCH chroma <= 0.026; the plum character lives in the text colour, the accent, the tint, the logo and the italic titles.
* **Text on fills follows the fill.** The dark accent fill is `#E58BB8` with ink text (7.49:1); light keeps white on plum (10.65:1). Amber is a fill with ink text, and as text on light it is `#6B4300`.
* **Nothing by colour alone.** Every state shows a glyph and a word beside its colour.
* **Contrast is tested, not reviewed.** `theme.rs` has a unit test with the table of (foreground, background, minimum) pairs of the tokens: 4.5 for text, 3 for control edges, and fails when a token moves out of range. Control edges and separators are different tokens.
* **The window follows the system appearance.** The default is `ThemePreference::System`, applied live with `Context::set_theme` and `set_visuals_of` for the two sets; the top-bar Light/Dark button is removed; Settings keeps Follow system / Light / Dark.
* **The canvas follows the theme by default** (`canvas` token); Settings may keep it dark. The categorical palette is the validated eight (`tokens.json`, `data`), a ninth region is told apart by hatch and number, never by a new hue. Viridis stays the sequential default.
* **Settings: Contrast Standard / High** (egui has no operating-system query), and **Reduce motion**.

## Not in this decision

* **Type scale, control sizes, menu bar, keyboard map and copy** (plan U2 and U3): they follow the same review but need no ADR of their own unless the menu is made native (a new dependency, `muda`; that is a separate ADR).
* **Blur and translucent materials.** Rejected: egui does not draw them and the HIG limits them to a floating functional layer nyx does not have.
* **A user-editable palette or accent colour.** Not now; the roles make it possible later.
* **Pixel-exact screenshot tests.** The scripted shots are for review and a smoke test, not goldens.
