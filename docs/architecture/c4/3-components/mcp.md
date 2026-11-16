# Level 3: components of nyx-mcp

Crate `crates/nyx-mcp`. The MCP server (stdio, JSON-RPC 2.0): every tool is one engine command, with the JSON Schema the engine generates.

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

        14["<div style='font-weight: bold'>Protocol loop</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Newline-delimited JSON-RPC<br />2.0 over stdin/stdout:<br />initialize with protocol<br />negotiation, method dispatch,<br />notifications, logs on<br />stderr. No async runtime and<br />no SDK.</div>"]
        style 14 fill:#85bbf0,stroke:#5d82a8,color:#000000
        15["<div style='font-weight: bold'>Tool bridge</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The tool list is the engine's<br />command catalog (with<br />annotations: read-only,<br />destructive); a tool call is<br />Engine::call, its reply<br />rendered as text.</div>"]
        style 15 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      22["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 22 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    14-. "<div>Server instructions text</div><div style='font-size: 70%'>[declared in spec]</div>" .->22
    15-. "<div>One tool per command: catalog<br />and call</div><div style='font-size: 70%'>[declared in spec]</div>" .->22

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Protocol loop](../4-code/mcp/protocol.md) | 1 | 254 | 15, 57 | Newline-delimited JSON-RPC 2.0 over stdin/stdout: initialize with protocol negotiation, method dispatch, notifications, logs on stderr. No async runtime and no SDK. |
| [Tool bridge](../4-code/mcp/tools.md) | 1 | 254 | 15, 57 | The tool list is the engine's command catalog (with annotations: read-only, destructive); a tool call is Engine::call, its reply rendered as text. |
