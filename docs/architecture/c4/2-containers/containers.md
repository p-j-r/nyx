# Level 2: containers

One container per crate of the workspace, plus the files they write. The arrows are exactly the layering table of `crates/nyx-engine/tests/layering.rs`. Source: the `Containers` view. Export: [Mermaid](../export/mermaid/structurizr-Containers.mmd) · [PlantUML](../export/plantuml/structurizr-Containers.puml)

```mermaid
graph LR
  linkStyle default fill:#ffffff

  subgraph diagram ["Container View: nyx"]
    style diagram fill:#ffffff,stroke:#ffffff

    1["<div style='font-weight: bold'>Researcher</div><div style='font-size: 70%; margin-top: 0px'>[Person]</div><div style='font-size: 80%; margin-top:10px'>Designs waveguides, rings and<br />phononic cells; runs mode and<br />eigenfrequency studies;<br />compares runs.</div>"]
    style 1 fill:#1d4e89,stroke:#14365f,color:#ffffff
    2["<div style='font-weight: bold'>AI host</div><div style='font-size: 70%; margin-top: 0px'>[Software System]</div><div style='font-size: 80%; margin-top:10px'>Claude Code or Claude Desktop<br />(any MCP host). Drives nyx<br />through the engine commands<br />as MCP tools.</div>"]
    style 2 fill:#8a8a8a,stroke:#606060,color:#ffffff
    3["<div style='font-weight: bold'>Results viewer</div><div style='font-size: 70%; margin-top: 0px'>[Software System]</div><div style='font-size: 80%; margin-top:10px'>ParaView or VisIt. Reads the<br />XDMF description and the HDF5<br />results file.</div>"]
    style 3 fill:#8a8a8a,stroke:#606060,color:#ffffff
    4[("<div style='font-weight: bold'>Workspace folder</div><div style='font-size: 70%; margin-top: 0px'>[Software System]</div><div style='font-size: 80%; margin-top:10px'>Documents/nyx on disk:<br />workflows, runs, .nyx<br />archives, project files.</div>")]
    style 4 fill:#8a8a8a,stroke:#606060,color:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      117["<div style='font-weight: bold'>geom</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Lean CSG-over-B-rep geometry<br />kernel: feature lists, named<br />regions and boundaries,<br />selectors, exact projection<br />oracles. Depends on serde and<br />optional rayon only.</div>"]
      style 117 fill:#438dd5,stroke:#2e6295,color:#ffffff
      125[("<div style='font-weight: bold'>Project file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .nyx.json]</div><div style='font-size: 80%; margin-top:10px'>Everything saved for one<br />model: geometry, materials,<br />boundary conditions, studies.<br />JSON with a format key.</div>")]
      style 125 fill:#438dd5,stroke:#2e6295,color:#ffffff
      126[("<div style='font-weight: bold'>Results file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .h5 and .xdmf]</div><div style='font-size: 80%; margin-top:10px'>Project, mesh and every<br />solution of every solved<br />physics, compressed, exact;<br />with an XDMF description next<br />to it.</div>")]
      style 126 fill:#438dd5,stroke:#2e6295,color:#ffffff
      127[("<div style='font-weight: bold'>Workspace data</div><div style='font-size: 70%; margin-top: 0px'>[Container: Files and folders]</div><div style='font-size: 80%; margin-top:10px'>workspace.json,<br />workflow.json, history.jsonl,<br />runs.jsonl, specs/<br />(experiment.nyx.json,<br />specs/<hash>.nyx.json), run<br />folders with results.h5, .nyx<br />archives.</div>")]
      style 127 fill:#438dd5,stroke:#2e6295,color:#ffffff
      128[("<div style='font-weight: bold'>Material library</div><div style='font-size: 70%; margin-top: 0px'>[Container: Files]</div><div style='font-size: 80%; margin-top:10px'>materials/*.toml next to the<br />program: silicon, silica,<br />sapphire, metals, nitrides,<br />niobate and more, with<br />sources.</div>")]
      style 128 fill:#438dd5,stroke:#2e6295,color:#ffffff
      13("<div style='font-weight: bold'>nyx-mcp</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust binary]</div><div style='font-size: 80%; margin-top:10px'>The MCP server (stdio,<br />JSON-RPC 2.0): four script<br />tools (script_run,<br />script_wait, script_check,<br />api) over one script kernel;<br />the agent writes Python that<br />calls the engine commands and<br />only what it prints comes<br />back.</div>")
      style 13 fill:#438dd5,stroke:#2e6295,color:#ffffff
      16("<div style='font-weight: bold'>nyx-script</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library and worker binary]</div><div style='font-size: 80%; margin-top:10px'>Scripts, a third front-end of<br />the engine: Python-syntax<br />text (*.nys) run by Monty in<br />a worker process; each<br />nyx.<command>(...) is one<br />Engine::call. Limits, no<br />ambient access, the engine's<br />path guard.</div>")
      style 16 fill:#438dd5,stroke:#2e6295,color:#ffffff
      22("<div style='font-weight: bold'>nyx-connect</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library and binary]</div><div style='font-size: 80%; margin-top:10px'>Finds Claude Code and Claude<br />Desktop installs and<br />registers nyx-mcp with them;<br />the Connect menu and its<br />command-line twin.</div>")
      style 22 fill:#438dd5,stroke:#2e6295,color:#ffffff
      28["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 28 fill:#438dd5,stroke:#2e6295,color:#ffffff
      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      48["<div style='font-weight: bold'>nyx-dfs</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Diffusion physics: heat and<br />electrostatics. Typed<br />problem, three-stage kernels,<br />descriptors, goldens.</div>"]
      style 48 fill:#2e7d32,stroke:#205723,color:#ffffff
      55["<div style='font-weight: bold'>nyx-sms</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Solid mechanics: static,<br />eigenfrequency, PML, Floquet,<br />low-reflecting boundaries,<br />axisymmetric eigenfrequency<br />study (bodies of revolution,<br />axis conditions, cylindrical<br />layers).</div>"]
      style 55 fill:#2e7d32,stroke:#205723,color:#ffffff
      6("<div style='font-weight: bold'>app</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust, eframe/egui]</div><div style='font-size: 80%; margin-top:10px'>The egui window (binary nyx):<br />geometry editor, mesh view,<br />Study, Physics and Materials<br />tabs drawn from the physics<br />descriptors, results,<br />analysis, workflow dashboard.<br />A front-end of the engine: it<br />holds an Engine and calls<br />commands.</div>")
      style 6 fill:#438dd5,stroke:#2e6295,color:#ffffff
      65["<div style='font-weight: bold'>nyx-ems</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>EM waves: waveguide modes<br />(real and complex),<br />axisymmetric, cavity in 3D.</div>"]
      style 65 fill:#2e7d32,stroke:#205723,color:#ffffff
      75["<div style='font-weight: bold'>nyx-analysis</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Analyses of results: the<br />polarization (Stokes<br />parameters, Poincare sphere,<br />classes, degenerate groups)<br />as a function of a<br />SolveResult and its mesh.</div>"]
      style 75 fill:#438dd5,stroke:#2e6295,color:#ffffff
      83["<div style='font-weight: bold'>nyx-store</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Workspace, workflow,<br />experiment and run storage:<br />manifests, append-only run<br />log, versions, search, .nyx<br />archives. Knows nothing of a<br />project.</div>"]
      style 83 fill:#438dd5,stroke:#2e6295,color:#ffffff
      92["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 92 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    1-. "<div>Edits models, runs studies,<br />reads results</div><div style='font-size: 70%'></div>" .->6
    1-. "<div>Asks in words</div><div style='font-size: 70%'></div>" .->2
    2-. "<div>Calls tools</div><div style='font-size: 70%'>[MCP over stdio (JSON-RPC 2.0)]</div>" .->13
    3-. "<div>Reads</div><div style='font-size: 70%'>[XDMF / HDF5]</div>" .->126
    6-. "<div>Reads and writes</div><div style='font-size: 70%'></div>" .->4
    6-. "<div>Calls commands, polls jobs,<br />draws results</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    6-. "<div>Connect menu</div><div style='font-size: 70%'>[Rust API]</div>" .->22
    13-. "<div>The one engine the scripts<br />drive; the project header</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    13-. "<div>script_run, script_wait,<br />script_check</div><div style='font-size: 70%'>[Rust API]</div>" .->16
    16-. "<div>One engine command per script<br />call; the path guard</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    6-. "<div>nyx check and nyx run on a<br />script</div><div style='font-size: 70%'>[Rust API]</div>" .->16
    28-. "<div>Registry, ProblemSpec,<br />SolveResult</div><div style='font-size: 70%'>[Rust API]</div>" .->41
    28-. "<div>Compiled-in physics (feature<br />dfs)</div><div style='font-size: 70%'>[Rust API]</div>" .->48
    28-. "<div>Compiled-in physics (feature<br />sms)</div><div style='font-size: 70%'>[Rust API]</div>" .->55
    28-. "<div>Compiled-in physics (feature<br />ems)</div><div style='font-size: 70%'>[Rust API]</div>" .->65
    28-. "<div>Analyses of results</div><div style='font-size: 70%'>[Rust API]</div>" .->75
    28-. "<div>Workflows and runs; passes<br />spec text and a card</div><div style='font-size: 70%'>[Rust API]</div>" .->83
    28-. "<div>Meshing, materials, plots</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    28-. "<div>Feature lists, labels,<br />selectors</div><div style='font-size: 70%'>[Rust API]</div>" .->117
    48-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->41
    55-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->41
    65-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->41
    48-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    55-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    65-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    75-. "<div>Reads SolveResult</div><div style='font-size: 70%'>[Rust API]</div>" .->41
    75-. "<div>Mesh geometry, quadrature</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    92-. "<div>Evaluates models, reads<br />selectors and labels</div><div style='font-size: 70%'>[Rust API]</div>" .->117
    28-. "<div>Saves and loads</div><div style='font-size: 70%'>[serde JSON]</div>" .->125
    28-. "<div>Writes and opens</div><div style='font-size: 70%'>[HDF5 (serial or MPI-IO), XDMF]</div>" .->126
    28-. "<div>Loads</div><div style='font-size: 70%'>[TOML]</div>" .->128
    83-. "<div>Reads and writes</div><div style='font-size: 70%'></div>" .->127

  end```

| Container | Crate | What it is |
|---|---|---|
| [app](../3-components/app.md) | `crates/app` | The egui window (binary nyx): geometry editor, mesh view, Study, Physics and Materials tabs drawn from the physics descriptors, results, analysis, workflow dashboard. A front-end of the engine: it holds an Engine and calls commands. |
| [nyx-mcp](../3-components/mcp.md) | `crates/nyx-mcp` | The MCP server (stdio, JSON-RPC 2.0): four script tools (script_run, script_wait, script_check, api) over one script kernel; the agent writes Python that calls the engine commands and only what it prints comes back. |
| [nyx-script](../3-components/script.md) | `crates/nyx-script` | Scripts, a third front-end of the engine: Python-syntax text (*.nys) run by Monty in a worker process; each nyx.&lt;command&gt;(...) is one Engine::call. Limits, no ambient access, the engine's path guard. |
| [nyx-connect](../3-components/connect.md) | `crates/nyx-connect` | Finds Claude Code and Claude Desktop installs and registers nyx-mcp with them; the Connect menu and its command-line twin. |
| [nyx-engine](../3-components/engine.md) | `crates/nyx-engine` | The Engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them, and the binding to a workflow. Every action is a named command with a JSON Schema. |
| [nyx-physics](../3-components/physics.md) | `crates/nyx-physics` | The Physics trait, the problem envelope (ProblemSpec), descriptors, SolveResult, the registry type, validation, JSON Schemas. |
| [nyx-dfs](../3-components/dfs.md) | `crates/nyx-dfs` | Diffusion physics: heat and electrostatics. Typed problem, three-stage kernels, descriptors, goldens. |
| [nyx-sms](../3-components/sms.md) | `crates/nyx-sms` | Solid mechanics: static, eigenfrequency, PML, Floquet, low-reflecting boundaries, axisymmetric eigenfrequency study (bodies of revolution, axis conditions, cylindrical layers). |
| [nyx-ems](../3-components/ems.md) | `crates/nyx-ems` | EM waves: waveguide modes (real and complex), axisymmetric, cavity in 3D. |
| [nyx-analysis](../3-components/analysis.md) | `crates/nyx-analysis` | Analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a SolveResult and its mesh. |
| [nyx-store](../3-components/store.md) | `crates/nyx-store` | Workspace, workflow, experiment and run storage: manifests, append-only run log, versions, search, .nyx archives. Knows nothing of a project. |
| [geom](../3-components/geom.md) | `crates/geom` | Lean CSG-over-B-rep geometry kernel: feature lists, named regions and boundaries, selectors, exact projection oracles. Depends on serde and optional rayon only. |
| [femcore](../3-components/femcore.md) | `crates/femcore` | The FEM core: meshes, meshing and adaptation, linear algebra behind one driver, finite-element basics, the material specification and library, expressions. |
| Material library | (data) | materials/*.toml next to the program: silicon, silica, sapphire, metals, nitrides, niobate and more, with sources. |
| Project file | (data) | Everything saved for one model: geometry, materials, boundary conditions, studies. JSON with a format key. |
| Results file | (data) | Project, mesh and every solution of every solved physics, compressed, exact; with an XDMF description next to it. |
| Workspace data | (data) | workspace.json, workflow.json, history.jsonl, runs.jsonl, specs/ (experiment.nyx.json, specs/&lt;hash&gt;.nyx.json), run folders with results.h5, .nyx archives. |
