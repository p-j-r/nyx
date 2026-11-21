# Look and widgets: code

C4 level 4 (code) of the component **Look and widgets** (`app.style`) in the container **app** (crate `nyx`). [Component view of the container](../../3-components/app.md) · [All code pages](../README.md)

> The nyx look (the colour tokens as data, fonts, egui visuals, light and dark, contrast) and the small brand widgets (cards, chips, buttons, glyphs, sparklines, sheets), and how numbers and counts are written.

**Decisions:** [ADR 23](../../../adr/0023-app-is-a-thin-shell-over-one-project.md), [ADR 88](../../../adr/0088-colour-has-roles-and-the-window-follows-the-system.md)  
**Size:** 5 files, 1308 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/app/src/cvd.rs` | 90 |  |
| `crates/app/src/fmt.rs` | 223 |  |
| `crates/app/src/kit.rs` | 334 |  |
| `crates/app/src/theme.rs` | 334 |  |
| `crates/app/src/tokens.rs` | 327 |  |

## What it uses

Nothing outside the component (the standard library and external crates are not listed).

## Used by

| Component | Modules | Items |
|---|---|---|
| [Scripts and batch](./headless.md) | cvd | Kind, simulate_rgba |
| [Inspector tabs](./inspector.md) | fmt, theme | accent, activity_text, danger, on_accent, plural, raised, si, success +1 |
| [Shell and dashboard](./shell.md) | fmt, kit, theme | State, card, danger, ellipsize, key_label, modal, muted, plural +3 |
| [Viewport](./viewport.md) | fmt, theme, tokens | Ink, Mode, ellipsize, rgb, with_alpha |
| [Window](./window.md) | cvd, fmt, kit, theme | Kind, State, accent, activity_text, apply, card, danger, edge +9 |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace cvd_rs {
    class Kind {
        <<enumeration>>
        Protan
        Deutan
        Tritan
        +parse(s) Option~Kind~
    }
    class fns_cvd_rs {
        <<module cvd.rs>>
        +simulate(rgb, kind) slice
        +simulate_rgba(raw, kind)
    }
  }
  namespace fmt_rs {
    class fns_fmt_rs {
        <<module fmt.rs>>
        +plural(n, noun) String
        +ellipsize(s, max) String
        +seconds(s) String
        +si(v, unit) String
        +unit_of(key) Option
        +key_label(key) String
        +value(key, v) String
    }
  }
  namespace kit_rs {
    class State {
        <<enumeration>>
        Done
        Running
        Waiting
        Failed
        Idle
        Checked
        Exploratory
        NotAffordable
        +word() str
    }
    class TabHit {
        <<struct>>
        +select Response
        +close Option~Response~
    }
    class fns_kit_rs {
        <<module kit.rs>>
        +card(ui, fill, stroke, pad, f) R
        +plain_card(ui, f) R
        +primary(ui, text) Response
        +primary_enabled(ui, enabled, text) Response
        +plain_card_r(ui, f) Response
        +ghost(ui, text) Response
        +danger(ui, text) Response
        +pill(ui, on, text) Response
        +chip(ui, text) Response
        +badge(ui, fill, fg, text)
        +more_11_functions()
    }
  }
  namespace theme_rs {
    class Look {
        <<struct>>
        +dark bool
        +high_contrast bool
        +reduce_motion bool
    }
    class Ink {
        <<struct>>
        +mode Mode
        +bg Color32
        +text Color32
        +muted Color32
        +chip Color32
        +line Color32
        +edge Color32
    }
    class fns_theme_rs {
        <<module theme.rs>>
        +title_family() FontFamily
        +title(size) FontId
        +rgb(c) Color32
        +is_dark(ctx) bool
        +mode(ctx) Mode
        +role(ctx, p) Color32
        +page(ctx) Color32
        +card(ctx) Color32
        +raised(ctx) Color32
        +line(ctx) Color32
        +more_25_functions()
    }
  }
  namespace tokens_rs {
    class Mode {
        <<enumeration>>
        Light
        Dark
    }
    class Pair {
        <<struct>>
        +light Rgb
        +dark Rgb
        +get(mode) Rgb
    }
    class fns_tokens_rs {
        <<module tokens.rs>>
        +hex(s) Rgb
        +regions(canvas) slice
        +Rgb type
        +PAGE const
        +CARD const
        +RAISED const
        +more_20_items
    }
  }
  Ink --> Mode
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `Kind` | `cvd.rs:10` | variants: `Protan`, `Deutan`, `Tritan` —  Which cone type is missing. |
| fn | `simulate` | `cvd.rs:51` | `pub fn simulate(rgb: [u8; 3], kind: Kind) -> [u8` —  One colour as the simulated viewer sees it. |
| fn | `simulate_rgba` | `cvd.rs:59` | `pub fn simulate_rgba(raw: &mut [u8], kind: Kind)` —  A whole image, `rgba` bytes in place (alpha is kept). |
| fn | `Kind::parse` | `cvd.rs:17` | `pub fn parse(s: &str) -> Option<Kind>` |
| fn | `plural` | `fmt.rs:8` | `pub fn plural(n: usize, noun: &str) -> String` —  "1 run", "2 runs", "0 runs": the noun gets an "s" unless the count is one. |
| fn | `ellipsize` | `fmt.rs:13` | `pub fn ellipsize(s: &str, max: usize) -> String` —  The text cut to `max` characters with an ellipsis; text that fits is returned as it is. |
| fn | `seconds` | `fmt.rs:21` | `pub fn seconds(s: f64) -> String` —  A time that has passed: "0.4 s", "12 s", "1:05", "1:02:09" (seconds below ten show a decimal; a minute and more show minutes and seconds). |
| fn | `si` | `fmt.rs:38` | `pub fn si(v: f64, unit: &str) -> String` —  A value with an SI prefix and four significant digits: `si(2.017e-4, "m")` is "201.7 µm", `si(6.175e6, "Pa")` is "6.175 MPa". |
| fn | `unit_of` | `fmt.rs:60` | `pub fn unit_of(key: &str) -> Option<&'static str>` —  The unit a metric key ends with ("max_u_m" is in metres), when the last word is a physical unit. |
| fn | `key_label` | `fmt.rs:89` | `pub fn key_label(key: &str) -> String` —  "max_von_mises_pa" as "Max von Mises [Pa]": first word capitalised, a few names spelled right, the unit in brackets. |
| fn | `value` | `fmt.rs:134` | `pub fn value(key: &str, v: f64) -> String` —  A number as the runs table shows it: SI prefixes when the key names a unit, whole numbers without decimals, lengths of a parameter in nm, |
| fn | `card` | `kit.rs:8` | `pub fn card<R>(ui: &mut Ui, fill: Color32, stroke: Color32, pad: i8, f: impl FnOnce(&mut Ui) -> R) -> R` |
| fn | `plain_card` | `kit.rs:13` | `pub fn plain_card<R>(ui: &mut Ui, f: impl FnOnce(&mut Ui) -> R) -> R` —  A white card with the line colour. |
| fn | `primary` | `kit.rs:19` | `pub fn primary(ui: &mut Ui, text: impl Into<String>) -> Response` —  The one filled accent button of a view, with the text colour that has contrast on it (`on_accent`). |
| fn | `primary_enabled` | `kit.rs:24` | `pub fn primary_enabled(ui: &mut Ui, enabled: bool, text: impl Into<String>) -> Response` |
| fn | `plain_card_r` | `kit.rs:29` | `pub fn plain_card_r(ui: &mut Ui, f: impl FnOnce(&mut Ui)) -> Response` —  A white card; returns the response of the whole card (for clicks). |
| fn | `ghost` | `kit.rs:37` | `pub fn ghost(ui: &mut Ui, text: impl Into<String>) -> Response` —  Outlined button. |
| fn | `danger` | `kit.rs:42` | `pub fn danger(ui: &mut Ui, text: impl Into<String>) -> Response` —  Outlined button in the danger colour (Cancel). |
| fn | `pill` | `kit.rs:48` | `pub fn pill(ui: &mut Ui, on: bool, text: impl Into<String>) -> Response` —  A selectable pill: plum when on. |
| fn | `chip` | `kit.rs:55` | `pub fn chip(ui: &mut Ui, text: impl Into<String>) -> Response` —  A small rounded label such as a tag. |
| fn | `badge` | `kit.rs:62` | `pub fn badge(ui: &mut Ui, fill: Color32, fg: Color32, text: impl Into<String>)` —  A chip with a colour (status chip in the top bar, badges on the open-file screen). |
| fn | `label_caps` | `kit.rs:67` | `pub fn label_caps(ui: &mut Ui, text: &str, hint: &str)` —  "RUNNING" style section label with an optional grey hint after it. |
| enum | `State` | `kit.rs:79` | variants: `Done`, `Running`, `Waiting`, `Failed`, `Idle`, `Checked`, `Exploratory`, `NotAffordable` —  What a workflow, a run or the solver is doing. Colour is never the only channel: every state has its own glyph (drawn, so it needs no font) |
| fn | `glyph` | `kit.rs:115` | `pub fn glyph(ui: &mut Ui, state: State) -> Response` —  The glyph of a state, 14 px square. |
| fn | `paint_glyph` | `kit.rs:122` | `pub fn paint_glyph(ui: &Ui, c: egui::Pos2, state: State)` |
| fn | `status` | `kit.rs:172` | `pub fn status(ui: &mut Ui, state: State, label: impl Into<String>)` —  A state with its word (or `label` instead of the word): glyph, then text in the colour that reads on the page. |
| struct | `TabHit` | `kit.rs:188` | fields: `select`, `close` —  What `tab` reports: a click on the tab and a click on its close target. |
| fn | `tab` | `kit.rs:195` | `pub fn tab(ui: &mut Ui, id: impl std::hash::Hash + std::fmt::Debug, title: &str, selected: bool, state: Option<State>, closable: bool) -> TabHit` —  One tab of the top bar: 28 px high, a state glyph, the title, and (when `closable`) a 24 x 24 close target. Selected = tint and a 2 px accent |
| fn | `mono` | `kit.rs:244` | `pub fn mono(text: impl Into<String>) -> RichText` |
| fn | `muted` | `kit.rs:248` | `pub fn muted(ui: &Ui, text: impl Into<String>) -> RichText` |
| fn | `title` | `kit.rs:252` | `pub fn title(text: impl Into<String>, size: f32) -> RichText` |
| fn | `sparkline` | `kit.rs:257` | `pub fn sparkline(ui: &mut Ui, values: &[f64], color: Color32, size: Vec2)` —  A line of values as a polyline in a box (dashboard cards). |
| fn | `bar` | `kit.rs:274` | `pub fn bar(ui: &mut Ui, fraction: Option<f32>, width: f32)` —  A determinate progress bar in plum on a pale track; `None` shows a sliding block (no percentage reported). |
| fn | `modal` | `kit.rs:304` | `pub fn modal(ctx: &egui::Context, id: &str, width: f32, f: impl FnOnce(&mut Ui)) -> bool` —  A modal sheet: dims the page, swallows clicks outside and shows `f` in a centred card. Returns `true` when it should close |
| fn | `State::word` | `kit.rs:100` | `pub fn word(self) -> &'static str` |
| fn | `title_family` | `theme.rs:13` | `pub fn title_family() -> FontFamily` —  The serif italic family of page titles. |
| fn | `title` | `theme.rs:17` | `pub fn title(size: f32) -> FontId` |
| fn | `rgb` | `theme.rs:22` | `pub fn rgb(c: Rgb) -> Color32` —  A token colour as an egui colour. |
| fn | `is_dark` | `theme.rs:27` | `pub fn is_dark(ctx: &egui::Context) -> bool` —  True when the colours are the dark ones (`apply_mode`). |
| fn | `mode` | `theme.rs:31` | `pub fn mode(ctx: &egui::Context) -> Mode` |
| fn | `role` | `theme.rs:36` | `pub fn role(ctx: &egui::Context, p: Pair) -> Color32` —  The value of a token in the look that is active. |
| fn | `page` | `theme.rs:40` | `pub fn page(ctx: &egui::Context) -> Color32` |
| fn | `card` | `theme.rs:43` | `pub fn card(ctx: &egui::Context) -> Color32` |
| fn | `raised` | `theme.rs:47` | `pub fn raised(ctx: &egui::Context) -> Color32` —  hover and pressed neutral fills, popups, chips |
| fn | `line` | `theme.rs:51` | `pub fn line(ctx: &egui::Context) -> Color32` —  Dividers and card outlines. Never the edge of a control: that is `edge`. |
| fn | `edge` | `theme.rs:55` | `pub fn edge(ctx: &egui::Context) -> Color32` —  The edge of buttons, fields and checkboxes (3:1 against the page and the card). |
| fn | `text` | `theme.rs:58` | `pub fn text(ctx: &egui::Context) -> Color32` |
| fn | `muted` | `theme.rs:61` | `pub fn muted(ctx: &egui::Context) -> Color32` |
| fn | `accent` | `theme.rs:64` | `pub fn accent(ctx: &egui::Context) -> Color32` |
| fn | `on_accent` | `theme.rs:68` | `pub fn on_accent(ctx: &egui::Context) -> Color32` —  Text on an accent fill: white on the plum, ink on the dark pink. |
| fn | `tint` | `theme.rs:71` | `pub fn tint(ctx: &egui::Context) -> Color32` |
| fn | `activity` | `theme.rs:75` | `pub fn activity(ctx: &egui::Context) -> Color32` —  Amber: time and attention (running, queued, stale) and nothing else. |
| fn | `on_activity` | `theme.rs:78` | `pub fn on_activity(ctx: &egui::Context) -> Color32` |
| fn | `activity_text` | `theme.rs:82` | `pub fn activity_text(ctx: &egui::Context) -> Color32` —  The activity colour as text (on `warm`). |
| fn | `warm` | `theme.rs:86` | `pub fn warm(ctx: &egui::Context) -> Color32` —  Background of the busy chip and the stale band. |
| fn | `success` | `theme.rs:89` | `pub fn success(ctx: &egui::Context) -> Color32` |
| fn | `danger` | `theme.rs:92` | `pub fn danger(ctx: &egui::Context) -> Color32` |
| fn | `scrim` | `theme.rs:96` | `pub fn scrim(alpha: u8) -> Color32` —  The black scrim behind a sheet and under a shadow (a veil, not a palette colour). |
| struct | `Look` | `theme.rs:102` | fields: `dark`, `high_contrast`, `reduce_motion` —  The settings of the look that are not the light or dark colours. |
| fn | `look` | `theme.rs:115` | `pub fn look(ctx: &egui::Context) -> Look` —  The look that was installed last (`apply_look`). |
| fn | `reduce_motion` | `theme.rs:119` | `pub fn reduce_motion(ctx: &egui::Context) -> bool` |
| fn | `canvas_mode` | `theme.rs:124` | `pub fn canvas_mode(ctx: &egui::Context, follow_theme: bool) -> Mode` —  The mode of the canvas that is drawn: dark by default also in a light window, the mode of the look when `follow_theme`. |
| struct | `Ink` | `theme.rs:131` | fields: `mode`, `bg`, `text`, `muted`, `chip`, `line`, `edge` —  Everything drawn on the canvas, in the colours of the canvas and not of the window (a light window can hold a dark canvas): the viewport and the |
| fn | `ink` | `theme.rs:144` | `pub fn ink(ctx: &egui::Context, follow_theme: bool) -> Ink` |
| fn | `ink_of` | `theme.rs:148` | `pub fn ink_of(m: Mode) -> Ink` |
| fn | `with_alpha` | `theme.rs:154` | `pub fn with_alpha(c: Color32, a: u8) -> Color32` —  A colour with the given opacity (0 to 255), not premultiplied. |
| fn | `system_high_contrast` | `theme.rs:159` | `pub fn system_high_contrast() -> bool` —  True when the system asks for High Contrast: read from Windows (`SPI_GETHIGHCONTRAST`); other platforms have no such flag we read yet. |
| fn | `resolve_contrast` | `theme.rs:191` | `pub fn resolve_contrast(setting: &str) -> bool` —  The High contrast setting ("auto", "standard", "high") as a yes or no. |
| fn | `apply` | `theme.rs:200` | `pub fn apply(ctx: &egui::Context)` —  Light visuals (what tests and the first frame use). |
| fn | `apply_mode` | `theme.rs:247` | `pub fn apply_mode(ctx: &egui::Context, dark: bool)` —  Light or dark with the other options of the look as they were last set (`apply_look`). |
| fn | `apply_look` | `theme.rs:254` | `pub fn apply_look(ctx: &egui::Context, look: Look)` —  Installs the light and the dark look and makes `look.dark` the active one. Both are installed every time (it is cheap) so nothing depends on |
| fn | `install_fonts` | `theme.rs:299` | `pub fn install_fonts(ctx: &egui::Context)` —  Installs the fonts that are available (see the module comment). Anything missing falls back to egui's own fonts. |
| type | `Rgb` | `tokens.rs:8` |  An sRGB colour, 8 bits per channel. |
| enum | `Mode` | `tokens.rs:12` | variants: `Light`, `Dark` —  Which set of colours: the window follows the system, the canvas has its own choice. |
| struct | `Pair` | `tokens.rs:19` | fields: `light`, `dark` —  A role with one value per mode. |
| fn | `hex` | `tokens.rs:43` | `pub const fn hex(s: &str) -> Rgb` —  `"#RRGGBB"` to bytes, at compile time, so the table below reads like the JSON. |
| const | `PAGE` | `tokens.rs:54` | `Pair` —  window background |
| const | `CARD` | `tokens.rs:56` | `Pair` —  cards, sheets, fields |
| const | `RAISED` | `tokens.rs:58` | `Pair` —  hover and pressed neutral fills, popups |
| const | `SEPARATOR` | `tokens.rs:60` | `Pair` —  dividers and card outlines only; never a control edge |
| const | `CONTROL_EDGE` | `tokens.rs:62` | `Pair` —  edge of buttons, fields, combo boxes, checkboxes (3:1 on page and card) |
| const | `CONTROL_EDGE_HIGH` | `tokens.rs:64` | `Pair` —  the control edge of the High contrast setting (4.5:1 on page, card and raised) |
| const | `TEXT` | `tokens.rs:66` | `Pair` —  primary text |
| const | `MUTED` | `tokens.rs:68` | `Pair` —  secondary text |
| const | `ACCENT` | `tokens.rs:70` | `Pair` —  the brand colour: primary action fill, links, focus ring, selected bar |
| const | `ON_ACCENT` | `tokens.rs:72` | `Pair` —  text on an accent fill (white on the dark pink fails, ink passes) |
| const | `TINT` | `tokens.rs:74` | `Pair` —  selected item background (with a 2 px accent bar) |
| const | `ACTIVITY` | `tokens.rs:76` | `Pair` —  running, queued, stale: fill |
| const | `ON_ACTIVITY` | `tokens.rs:78` | `Pair` —  text on an activity fill |
| const | `ACTIVITY_TEXT` | `tokens.rs:80` | `Pair` —  the activity colour used as text (on `ACTIVITY_TINT`) |
| const | `ACTIVITY_TINT` | `tokens.rs:82` | `Pair` —  background of the busy chip and the stale band |
| const | `SUCCESS` | `tokens.rs:84` | `Pair` —  saved, converged, passed (always with a glyph and a word) |
| const | `DANGER` | `tokens.rs:86` | `Pair` —  failed, destructive (always with a glyph and a word) |
| const | `CANVAS` | `tokens.rs:88` | `Pair` —  viewport and plot background (default dark in both themes) |
| const | `CHIP` | `tokens.rs:90` | `Pair` —  region number chip background: the canvas colour of the same mode, with `TEXT` on it |
| const | `REGION_GREY` | `tokens.rs:92` | `Pair` —  fill of regions beyond the sixteenth (slot 17 and up): neutral, told apart by the number chip only |
| const | `ROLES` | `tokens.rs:96` | `[(&str, Pair)` —  Every role by its name in `tokens.json`. |
| const | `REGIONS_ON_DARK_CANVAS` | `tokens.rs:120` | `[Rgb` —  Eight region colours for the dark canvas (the default canvas, also in a light window). Any two can touch: checked with all pairs. |
| const | `REGIONS_ON_LIGHT_CANVAS` | `tokens.rs:124` | `[Rgb` —  Eight region colours for the light canvas (only when Settings says the canvas follows a light theme). |
| fn | `regions` | `tokens.rs:128` | `pub fn regions(canvas: Mode) -> &'static [Rgb` —  The region palette for the canvas that is drawn. |
| fn | `Pair::get` | `tokens.rs:25` | `pub const fn get(&self, mode: Mode) -> Rgb` |
