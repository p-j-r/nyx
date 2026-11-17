# 88. Colour has roles, and the window follows the system appearance

Date: 2026-11-17

## Status

Accepted by Paul on 2026-11-17 (the merged plan, `docs/ux/UX_MERGED_PLAN.md`), implemented step by step: tokens and looks in U1, system theme, settings, menu and About in U2. **The number 88 was chosen by Paul on 2026-11-17**: it is the number the productivity line (84 to 89) left unused. Evidence: `docs/ux/UX_REVIEW.md`; tokens: `docs/ux/tokens.json` (equal to `crates/app/src/tokens.rs`, a test compares them).

## Context

The window draws every colour from `crates/app/src/theme.rs`, which is good, and then breaks its own rule in 56 places (`Color32` literals in `app.rs`, `tabs.rs`, `sheets.rs` and others). A review of `productivity-4-memory` against Apple's HIG (principles, read from the skill's reference pages) and the data-visualisation checks found:

* white on the dark accent is 3.21:1 (needs 4.5), the dark red is 2.98:1, control edges are 1.17 to 1.44:1 against their surface;
* the same colour carries several meanings (plum is brand, primary, selected and link; amber is running, stale and warning; there are two Compute buttons in two greens and plums), and up to four filled accent objects share one dark view, which is what makes the palette tiring;
* the top bar has its own Light/Dark button although Settings already has Follow system (HIG `dark-mode.md`: avoid an app-specific appearance setting);
* the viewport is always near-black, and its 12-colour domain palette fails the validator (adjacent colour-blind separation 4.3, floor 6).

## Decision

* **Colour is a set of roles, not values.** `crates/app/src/tokens.rs` holds every role in a light and a dark value as plain data (the table of `docs/ux/tokens.json`: page, card, raised, separator, control edge, text, muted, accent, on-accent, tint, activity, success, danger, canvas, chip); `theme.rs` turns it into the two `Visuals` and the role functions the window calls.
* **One colour, one meaning.** Accent = the primary action, links, focus and the selected-item bar. Activity (amber) = time and attention: running, queued, stale. Success and danger = outcomes. Data colours = plots and the canvas only.
* **Quiet by default.** At most one filled accent object per view (Compute). Selection is a tint with a 2 px accent bar. Neutral surfaces carry OKLCH chroma <= 0.026; the plum character lives in the text colour, the accent, the tint, the logo and the italic titles.
* **Text on fills follows the fill.** The brand pink is kept (Paul, 2026-11-17): plum `#820C4F` with white text in the light look, `#D96AA4` with ink text in the dark look (white on that pink is 3.21:1 and fails, ink is 5.62:1). Amber is a fill with ink text, and as text on light it is `#6B4300` on its tint.
* **Nothing by colour alone.** Every state shows a glyph and a word beside its colour.
* **Contrast is tested, not reviewed.** `theme.rs` has a unit test with the table of (foreground, background, minimum) pairs of the tokens: 4.5 for text, 3 for control edges, and fails when a token moves out of range. Control edges and separators are different tokens.
* **The window follows the system appearance.** The default is `ThemePreference::System`, applied live with `Context::set_theme` and `set_visuals_of` for the two sets; the top-bar Light/Dark button is removed; Settings keeps Follow system / Light / Dark.
* **The canvas is dark by default, in both looks** (decision D2); Settings can make it follow the look. Regions are drawn opaque with the validated eight colours of `tokens.json` (checked with all pairs, colour-blind and normal vision), every region carries a number chip and a legend row, and beyond eight the same hues are hatched (9 to 16) and then grey with the number (17 and more).
* **Settings: Theme (follow the system by default), Contrast (follow system, standard, high), Reduce motion, Canvas, Colour regions by.** High contrast is read from Windows (`SPI_GETHIGHCONTRAST`, one small FFI call) and is a setting everywhere; zoom set with Ctrl and +/- is kept.
* **One table of actions, two views** (`crates/app/src/menu.rs`): the command palette and the native menu bar list the same actions. The menu is native through `muda` (a new dependency, Windows only for now: it is the platform menu, so it is keyboard-, screen-reader- and High-contrast-correct without drawing it ourselves; the keys stay handled by the window). Linux (GTK) and macOS (main thread) are not tried yet: they keep the palette.

## Not in this decision

* **Type scale, control sizes (28 px), keyboard map and copy** (plan U2 and U3): they follow the same review and need no ADR of their own; the type scale and the 28 px control height are set in `theme.rs`.
* **Blur and translucent materials.** Rejected: egui does not draw them and the HIG limits them to a floating functional layer nyx does not have.
* **A user-editable palette or accent colour.** Not now; the roles make it possible later.
* **Pixel-exact screenshot tests.** The scripted shots are for review and a smoke test, not goldens.
