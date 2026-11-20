# Level 3: components of nyx-connect

Crate `crates/nyx-connect`. Finds Claude Code and Claude Desktop installs and registers nyx-mcp with them; the Connect menu and its command-line twin.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_connect.mmd) · [PlantUML](../export/plantuml/structurizr-Components_connect.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-connect"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 22 ["nyx-connect"]
        style 22 fill:#ffffff,stroke:#2e6295,color:#2e6295

        23["<div style='font-weight: bold'>Host table</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Where each AI host (Claude<br />Code, Claude Desktop, ...)<br />keeps its MCP config and what<br />an entry looks like; the home<br />folder, application-data<br />folders and the operating<br />system.</div>"]
        style 23 fill:#85bbf0,stroke:#5d82a8,color:#000000
        24["<div style='font-weight: bold'>Host config</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Reads and rewrites a host's<br />JSON config: status of the<br />nyx entry, connect,<br />disconnect, list.</div>"]
        style 24 fill:#85bbf0,stroke:#5d82a8,color:#000000
        25["<div style='font-weight: bold'>Server probe</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Finds the nyx-mcp program,<br />builds the snippet for a host<br />that is not in the table, and<br />does a real MCP handshake<br />against it.</div>"]
        style 25 fill:#85bbf0,stroke:#5d82a8,color:#000000
        26["<div style='font-weight: bold'>Report</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The machine-readable form of<br />what nyx-connect finds and<br />does, with its JSON Schema.</div>"]
        style 26 fill:#85bbf0,stroke:#5d82a8,color:#000000
        27["<div style='font-weight: bold'>Command line</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The command-line twin of the<br />Connect menu: list, connect,<br />disconnect, test, optionally<br />as JSON.</div>"]
        style 27 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    27-. "<div>uses config</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    27-. "<div>uses env</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    27-. "<div>uses report</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    27-. "<div>uses server</div><div style='font-size: 70%'>[Rust paths]</div>" .->25
    24-. "<div>uses env</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    24-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    26-. "<div>uses config</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    26-. "<div>uses env</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    26-. "<div>uses server</div><div style='font-size: 70%'>[Rust paths]</div>" .->25
    25-. "<div>uses env</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    25-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->26

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Host table](../4-code/connect/hosts.md) | 2 | 219 | 28 | Where each AI host (Claude Code, Claude Desktop, ...) keeps its MCP config and what an entry looks like; the home folder, application-data folders and the operating system. |
| [Host config](../4-code/connect/config.md) | 1 | 174 | 28 | Reads and rewrites a host's JSON config: status of the nyx entry, connect, disconnect, list. |
| [Server probe](../4-code/connect/server.md) | 1 | 173 | 28 | Finds the nyx-mcp program, builds the snippet for a host that is not in the table, and does a real MCP handshake against it. |
| [Report](../4-code/connect/report.md) | 2 | 189 | 28, 19 | The machine-readable form of what nyx-connect finds and does, with its JSON Schema. |
| [Command line](../4-code/connect/cli.md) | 1 | 172 | 28 | The command-line twin of the Connect menu: list, connect, disconnect, test, optionally as JSON. |
