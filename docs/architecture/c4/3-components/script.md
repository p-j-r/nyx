# Level 3: components of nyx-script

Crate `crates/nyx-script`. Scripts, a third front-end of the engine: Python-syntax text (*.nys) run by Monty in a worker process; each nyx.<command>(...) is one Engine::call. Limits, no ambient access, the engine's path guard.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_script.mmd) · [PlantUML](../export/plantuml/structurizr-Components_script.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-script"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 16 ["nyx-script"]
        style 16 fill:#ffffff,stroke:#2e6295,color:#2e6295

        17["<div style='font-weight: bold'>Script host</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Runs a script (Python syntax,<br />`*.nys`) through a worker<br />process and answers each<br />`nyx.<command>(...)` it makes<br />with Engine::call; the<br />limits, the clock and random<br />policy, the path guard turned<br />on for the run, the script<br />kept with the runs it makes;<br />`check` parses and<br />type-checks without running.<br />JSON to and from script<br />values.</div>"]
        style 17 fill:#85bbf0,stroke:#5d82a8,color:#000000
        18["<div style='font-weight: bold'>Type stub</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The `nyx` object written as<br />Python declarations from<br />Engine::catalog(): the static<br />check, and the file an editor<br />reads.</div>"]
        style 18 fill:#85bbf0,stroke:#5d82a8,color:#000000
        19["<div style='font-weight: bold'>Worker program</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The only place the Monty<br />interpreter runs: Monty's<br />worker mode over stdin and<br />stdout, with the allocator<br />that bounds its memory. A<br />script that exhausts memory<br />or the stack ends this<br />process, not the engine's.</div>"]
        style 19 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      26["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 26 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    17-. "<div>uses engine: engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    17-. "<div>uses stubs</div><div style='font-size: 70%'>[Rust paths]</div>" .->18
    18-. "<div>uses engine: commands, engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    17-. "<div>Frames over a pipe: the<br />script pauses at each command<br />and the host answers</div><div style='font-size: 70%'>[declared in spec]</div>" .->19

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Script host](../4-code/script/host.md) | 3 | 334 | 93, 57 | Runs a script (Python syntax, `*.nys`) through a worker process and answers each `nyx.<command>(...)` it makes with Engine::call; the limits, the clock and random policy, the path guard turned on for the run, the script kept with the runs it makes; `check` parses and type-checks without running. JSON to and from script values. |
| [Type stub](../4-code/script/stubs.md) | 1 | 72 | 93 | The `nyx` object written as Python declarations from Engine::catalog(): the static check, and the file an editor reads. |
| [Worker program](../4-code/script/worker.md) | 1 | 118 | 93, 86 | The only place the Monty interpreter runs: Monty's worker mode over stdin and stdout, with the allocator that bounds its memory. A script that exhausts memory or the stack ends this process, not the engine's. |
