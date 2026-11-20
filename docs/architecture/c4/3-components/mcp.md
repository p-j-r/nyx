# Level 3: components of nyx-mcp

Crate `crates/nyx-mcp`. The MCP server (stdio, JSON-RPC 2.0): four script tools (script_run, script_wait, script_check, api) over one script kernel; the agent writes Python that calls the engine commands and only what it prints comes back.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_mcp.mmd) · [PlantUML](../export/plantuml/structurizr-Components_mcp.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-mcp"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 13 ["nyx-mcp"]
        style 13 fill:#ffffff,stroke:#2e6295,color:#2e6295

        14["<div style='font-weight: bold'>Protocol loop</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Newline-delimited JSON-RPC<br />2.0 over stdio, no async<br />runtime and no SDK;<br />initialize, ping, tools/list<br />and tools/call. `--dir`<br />(workflows and saved<br />scripts), `--root` (the<br />folder a script's file<br />arguments may name); `--job`<br />is the engine's solve worker.</div>"]
        style 14 fill:#85bbf0,stroke:#5d82a8,color:#000000
        15["<div style='font-weight: bold'>Script tools</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Four tools instead of one per<br />command: `script_run` (a<br />script on the kernel;<br />`params`, `path`, `save_as`,<br />`wait_s`,<br />`allow_destructive`),<br />`script_wait` (continue or<br />cancel a script that is<br />waiting for a job),<br />`script_check` (the static<br />check) and `api` (the guide,<br />the command list, one<br />command). One answer, as<br />text, with the project<br />header.</div>"]
        style 15 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      16("<div style='font-weight: bold'>nyx-script</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library and worker binary]</div><div style='font-size: 80%; margin-top:10px'>Scripts, a third front-end of<br />the engine: Python-syntax<br />text (*.nys) run by Monty in<br />a worker process; each<br />nyx.<command>(...) is one<br />Engine::call. Limits, no<br />ambient access, the engine's<br />path guard.</div>")
      style 16 fill:#438dd5,stroke:#2e6295,color:#ffffff
      28["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 28 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    16-. "<div>One engine command per script<br />call; the path guard</div><div style='font-size: 70%'>[Rust API]</div>" .->28
    14-. "<div>Server instructions text</div><div style='font-size: 70%'>[declared in spec]</div>" .->28
    15-. "<div>The one engine the scripts<br />drive, and its project header</div><div style='font-size: 70%'>[declared in spec]</div>" .->28
    15-. "<div>script_run, script_wait and<br />script_check drive the<br />kernel; the kernel is made at<br />the first script</div><div style='font-size: 70%'>[declared in spec]</div>" .->16
    15-. "<div>api: the command list and one<br />command's documentation</div><div style='font-size: 70%'>[declared in spec]</div>" .->16

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Protocol loop](../4-code/mcp/protocol.md) | 1 | 709 | 15, 57, 94 | Newline-delimited JSON-RPC 2.0 over stdio, no async runtime and no SDK; initialize, ping, tools/list and tools/call. `--dir` (workflows and saved scripts), `--root` (the folder a script's file arguments may name); `--job` is the engine's solve worker. |
| [Script tools](../4-code/mcp/tools.md) | 1 | 709 | 94, 57 | Four tools instead of one per command: `script_run` (a script on the kernel; `params`, `path`, `save_as`, `wait_s`, `allow_destructive`), `script_wait` (continue or cancel a script that is waiting for a job), `script_check` (the static check) and `api` (the guide, the command list, one command). One answer, as text, with the project header. |
