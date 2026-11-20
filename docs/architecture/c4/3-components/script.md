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
        18["<div style='font-weight: bold'>Script kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>One Monty session in a worker<br />process that outlives a<br />script: `run`, `wait`,<br />`cancel` and `check`;<br />variables stay from one<br />script to the next; a command<br />that waits for a job<br />(`solve_run`,<br />`convergence_check`,<br />`job_wait`, `queue_wait`) is<br />started and polled by the<br />host so no call blocks for<br />ever and a script can be left<br />`running`; `nyx.aio` futures;<br />the destructive commands<br />refused unless allowed; a<br />dead worker replaced from a<br />dump.</div>"]
        style 18 fill:#85bbf0,stroke:#5d82a8,color:#000000
        19["<div style='font-weight: bold'>Type stub</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The `nyx` object written as<br />Python declarations from<br />Engine::catalog(): the static<br />check, and the file an editor<br />reads.</div>"]
        style 19 fill:#85bbf0,stroke:#5d82a8,color:#000000
        20["<div style='font-weight: bold'>Tracing of script functions</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>A Python function given as a<br />formula is called once with a<br />symbolic point and becomes<br />formula text (N3): the<br />prelude the worker is fed<br />(the symbolic value, the<br />formula functions, the `nyx`<br />wrapper), the names the stub<br />declares for it, and the<br />static guard that refuses a<br />branch or comparison on a<br />value of the point with its<br />line (Monty does not call<br />`__bool__` or `__lt__`).</div>"]
        style 20 fill:#85bbf0,stroke:#5d82a8,color:#000000
        21["<div style='font-weight: bold'>Worker program</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The only place the Monty<br />interpreter runs: Monty's<br />worker mode over stdin and<br />stdout, with the allocator<br />that bounds its memory. A<br />script that exhausts memory<br />or the stack ends this<br />process, not the engine's.</div>"]
        style 21 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      28["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 28 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    17-. "<div>uses engine: engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    17-. "<div>uses kernel</div><div style='font-size: 70%'>[Rust paths]</div>" .->18
    18-. "<div>uses engine: engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    18-. "<div>uses convert, crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->17
    18-. "<div>uses stubs</div><div style='font-size: 70%'>[Rust paths]</div>" .->19
    18-. "<div>uses trace</div><div style='font-size: 70%'>[Rust paths]</div>" .->20
    19-. "<div>uses engine: commands, engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    19-. "<div>uses kernel</div><div style='font-size: 70%'>[Rust paths]</div>" .->18
    19-. "<div>uses trace</div><div style='font-size: 70%'>[Rust paths]</div>" .->20
    20-. "<div>uses engine: engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    20-. "<div>uses kernel</div><div style='font-size: 70%'>[Rust paths]</div>" .->18
    17-. "<div>Frames over a pipe: the<br />script pauses at each command<br />and the host answers</div><div style='font-size: 70%'>[declared in spec]</div>" .->21

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Script host](../4-code/script/host.md) | 3 | 198 | 93, 57 | Runs a script (Python syntax, `*.nys`) through a worker process and answers each `nyx.<command>(...)` it makes with Engine::call; the limits, the clock and random policy, the path guard turned on for the run, the script kept with the runs it makes; `check` parses and type-checks without running. JSON to and from script values. |
| [Script kernel](../4-code/script/kernel.md) | 1 | 658 | 94, 93 | One Monty session in a worker process that outlives a script: `run`, `wait`, `cancel` and `check`; variables stay from one script to the next; a command that waits for a job (`solve_run`, `convergence_check`, `job_wait`, `queue_wait`) is started and polled by the host so no call blocks for ever and a script can be left `running`; `nyx.aio` futures; the destructive commands refused unless allowed; a dead worker replaced from a dump. |
| [Type stub](../4-code/script/stubs.md) | 1 | 92 | 93 | The `nyx` object written as Python declarations from Engine::catalog(): the static check, and the file an editor reads. |
| [Tracing of script functions](../4-code/script/tracing.md) | 1 | 414 | 90, 93 | A Python function given as a formula is called once with a symbolic point and becomes formula text (N3): the prelude the worker is fed (the symbolic value, the formula functions, the `nyx` wrapper), the names the stub declares for it, and the static guard that refuses a branch or comparison on a value of the point with its line (Monty does not call `__bool__` or `__lt__`). |
| [Worker program](../4-code/script/worker.md) | 1 | 118 | 93, 86 | The only place the Monty interpreter runs: Monty's worker mode over stdin and stdout, with the allocator that bounds its memory. A script that exhausts memory or the stack ends this process, not the engine's. |
