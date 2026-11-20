# Dynamic view: AICommand

An AI host builds a model and solves it through the same commands. Export: [Mermaid](../export/mermaid/structurizr-AICommand.mmd) · [PlantUML](../export/plantuml/structurizr-AICommand.puml)

```mermaid
graph LR
  linkStyle default fill:#ffffff

  subgraph diagram ["Dynamic View: nyx"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      110["<div style='font-weight: bold'>geom</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Lean CSG-over-B-rep geometry<br />kernel: feature lists, named<br />regions and boundaries,<br />selectors, exact projection<br />oracles. Depends on serde and<br />optional rayon only.</div>"]
      style 110 fill:#438dd5,stroke:#2e6295,color:#ffffff
      13("<div style='font-weight: bold'>nyx-mcp</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust binary]</div><div style='font-size: 80%; margin-top:10px'>The MCP server (stdio,<br />JSON-RPC 2.0): four script<br />tools (script_run,<br />script_wait, script_check,<br />api) over one script kernel;<br />the agent writes Python that<br />calls the engine commands and<br />only what it prints comes<br />back.</div>")
      style 13 fill:#438dd5,stroke:#2e6295,color:#ffffff
      28["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 28 fill:#438dd5,stroke:#2e6295,color:#ffffff
      62["<div style='font-weight: bold'>nyx-ems</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>EM waves: waveguide modes<br />(real and complex),<br />axisymmetric, cavity in 3D.</div>"]
      style 62 fill:#2e7d32,stroke:#205723,color:#ffffff
      76["<div style='font-weight: bold'>nyx-store</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Workspace, workflow,<br />experiment and run storage:<br />manifests, append-only run<br />log, versions, search, .nyx<br />archives. Knows nothing of a<br />project.</div>"]
      style 76 fill:#438dd5,stroke:#2e6295,color:#ffffff
      85["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 85 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    2["<div style='font-weight: bold'>AI host</div><div style='font-size: 70%; margin-top: 0px'>[Software System]</div><div style='font-size: 80%; margin-top:10px'>Claude Code or Claude Desktop<br />(any MCP host). Drives nyx<br />through the engine commands<br />as MCP tools.</div>"]
    style 2 fill:#8a8a8a,stroke:#606060,color:#ffffff

    2-. "<div>1. tools/call script_run (a<br />script that builds, meshes<br />and solves)</div><div style='font-size: 70%'>[MCP over stdio (JSON-RPC 2.0)]</div>" .->13
    13-. "<div>2. call(command, arguments)</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    28-. "<div>3. feature list edit and<br />evaluation</div><div style='font-size: 70%'>[Rust API]</div>" .->110
    28-. "<div>4. mesh</div><div style='font-size: 70%'>[Rust API]</div>" .->85
    28-. "<div>5. solve</div><div style='font-size: 70%'>[Rust API]</div>" .->62
    28-. "<div>6. auto_store: the workflow<br />shows in the dashboard</div><div style='font-size: 70%'>[Rust API]</div>" .->76

  end```
