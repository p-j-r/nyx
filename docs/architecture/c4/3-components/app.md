# Level 3: components of app

Crate `crates/app`. The egui window (binary nyx): geometry editor, mesh view, Study, Physics and Materials tabs drawn from the physics descriptors, results, analysis, workflow dashboard. A front-end of the engine: it holds an Engine and calls commands.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_app.mmd) · [PlantUML](../export/plantuml/structurizr-Components_app.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - app"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 6 ["app"]
        style 6 fill:#ffffff,stroke:#2e6295,color:#2e6295

        10["<div style='font-weight: bold'>Viewport</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The 2D/3D viewport: camera,<br />drawing of meshes and result<br />fields, picking, colour maps,<br />the colour slots of regions<br />(solid, hatched, grey) with<br />their number chips and<br />legend, and what is kept<br />between frames.</div>"]
        style 10 fill:#85bbf0,stroke:#5d82a8,color:#000000
        11["<div style='font-weight: bold'>Look and widgets</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The nyx look (the colour<br />tokens as data, fonts, egui<br />visuals, light and dark,<br />contrast) and the small brand<br />widgets (cards, chips,<br />buttons, glyphs, sparklines,<br />sheets), and how numbers and<br />counts are written.</div>"]
        style 11 fill:#85bbf0,stroke:#5d82a8,color:#000000
        12["<div style='font-weight: bold'>Scripts and batch</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Runs without a person at the<br />window: scripted runs with<br />screenshots (--script), batch<br />runs for clusters, command<br />lists (nyx check, nyx run)<br />and the command reference<br />(nyx api), and the optional<br />user Python script.</div>"]
        style 12 fill:#85bbf0,stroke:#5d82a8,color:#000000
        7["<div style='font-weight: bold'>Window</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The application object around<br />one Engine: window layout,<br />the run-time state (mesh on<br />show, result on show, plot<br />settings) and the frame<br />update loop.</div>"]
        style 7 fill:#85bbf0,stroke:#5d82a8,color:#000000
        8["<div style='font-weight: bold'>Shell and dashboard</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The tab bar, the native menu<br />and the table of actions it<br />shares with the palette, the<br />dashboard of workflows with<br />live status, the busy chip<br />and queue popover read from<br />the engine queue, the<br />workflow view (concept, build<br />with its Sweep section,<br />analysis), the sheets (new,<br />search, versions, share,<br />settings) and the example<br />starting points.</div>"]
        style 8 fill:#85bbf0,stroke:#5d82a8,color:#000000
        9["<div style='font-weight: bold'>Inspector tabs</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The tabs beside the viewport:<br />geometry editor, materials,<br />physics, study, results,<br />analysis. Forms are drawn<br />from the physics descriptors;<br />the geometry editor generates<br />one form per feature kind.</div>"]
        style 9 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      16("<div style='font-weight: bold'>nyx-script</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library and worker binary]</div><div style='font-size: 80%; margin-top:10px'>Scripts, a third front-end of<br />the engine: Python-syntax<br />text (*.nys) run by Monty in<br />a worker process; each<br />nyx.<command>(...) is one<br />Engine::call. Limits, no<br />ambient access, the engine's<br />path guard.</div>")
      style 16 fill:#438dd5,stroke:#2e6295,color:#ffffff
      22("<div style='font-weight: bold'>nyx-connect</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library and binary]</div><div style='font-size: 80%; margin-top:10px'>Finds Claude Code and Claude<br />Desktop installs and<br />registers nyx-mcp with them;<br />the Connect menu and its<br />command-line twin.</div>")
      style 22 fill:#438dd5,stroke:#2e6295,color:#ffffff
      28["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 28 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    16-. "<div>One engine command per script<br />call; the path guard</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    12-. "<div>uses engine: apiref,<br />commandlist +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    12-. "<div>uses script: crate root, host<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->16
    12-. "<div>uses examples, shell</div><div style='font-size: 70%'>[Rust paths]</div>" .->8
    12-. "<div>uses cvd</div><div style='font-size: 70%'>[Rust paths]</div>" .->11
    12-. "<div>uses app, state</div><div style='font-size: 70%'>[Rust paths]</div>" .->7
    9-. "<div>uses engine: crate root,<br />project</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    9-. "<div>uses examples</div><div style='font-size: 70%'>[Rust paths]</div>" .->8
    9-. "<div>uses fmt, theme</div><div style='font-size: 70%'>[Rust paths]</div>" .->11
    9-. "<div>uses colormap</div><div style='font-size: 70%'>[Rust paths]</div>" .->10
    9-. "<div>uses app, spans +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->7
    8-. "<div>uses engine: crate root,<br />engine +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    8-. "<div>uses fmt, kit +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->11
    8-. "<div>uses app, spans +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->7
    10-. "<div>uses engine: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    10-. "<div>uses fmt, theme +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->11
    7-. "<div>uses connect: config, env +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->22
    7-. "<div>uses engine: crate root,<br />engine +8</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    7-. "<div>uses script: crate root, host</div><div style='font-size: 70%'>[Rust paths]</div>" .->16
    7-. "<div>uses pyscript</div><div style='font-size: 70%'>[Rust paths]</div>" .->12
    7-. "<div>uses geomui</div><div style='font-size: 70%'>[Rust paths]</div>" .->9
    7-. "<div>uses activity, examples +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->8
    7-. "<div>uses cvd, fmt +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->11
    7-. "<div>uses colormap, regions +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->10

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Window](../4-code/app/window.md) | 5 | 3851 | 23, 57, 92 | The application object around one Engine: window layout, the run-time state (mesh on show, result on show, plot settings) and the frame update loop. |
| [Shell and dashboard](../4-code/app/shell.md) | 8 | 5092 | 14, 23, 87, 88 | The tab bar, the native menu and the table of actions it shares with the palette, the dashboard of workflows with live status, the busy chip and queue popover read from the engine queue, the workflow view (concept, build with its Sweep section, analysis), the sheets (new, search, versions, share, settings) and the example starting points. |
| [Inspector tabs](../4-code/app/inspector.md) | 3 | 2536 | 23, 57 | The tabs beside the viewport: geometry editor, materials, physics, study, results, analysis. Forms are drawn from the physics descriptors; the geometry editor generates one form per feature kind. |
| [Viewport](../4-code/app/viewport.md) | 4 | 1598 | 23, 88 | The 2D/3D viewport: camera, drawing of meshes and result fields, picking, colour maps, the colour slots of regions (solid, hatched, grey) with their number chips and legend, and what is kept between frames. |
| [Look and widgets](../4-code/app/style.md) | 5 | 1308 | 23, 88 | The nyx look (the colour tokens as data, fonts, egui visuals, light and dark, contrast) and the small brand widgets (cards, chips, buttons, glyphs, sparklines, sheets), and how numbers and counts are written. |
| [Scripts and batch](../4-code/app/headless.md) | 5 | 733 | 23, 89, 93 | Runs without a person at the window: scripted runs with screenshots (--script), batch runs for clusters, command lists (nyx check, nyx run) and the command reference (nyx api), and the optional user Python script. |
