# Level 3: components of nyx-store

Crate `crates/nyx-store`. Workspace, workflow, experiment and run storage: manifests, append-only run log, versions, search, .nyx archives. Knows nothing of a project.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_store.mmd) · [PlantUML](../export/plantuml/structurizr-Components_store.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-store"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 73 ["nyx-store"]
        style 73 fill:#ffffff,stroke:#2e6295,color:#2e6295

        74["<div style='font-weight: bold'>Data model and card</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The workflow manifest,<br />experiments, run digest and<br />the concept card generated<br />from an experiment's spec.</div>"]
        style 74 fill:#85bbf0,stroke:#5d82a8,color:#000000
        75["<div style='font-weight: bold'>Workspace</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The folder behind the<br />dashboard: settings, the list<br />of workflows (manifests only)<br />and the sections of the<br />dashboard.</div>"]
        style 75 fill:#85bbf0,stroke:#5d82a8,color:#000000
        76["<div style='font-weight: bold'>Workflow folder</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>One workflow folder:<br />workflow.json, history.jsonl,<br />notes, specs and one folder<br />per experiment.</div>"]
        style 76 fill:#85bbf0,stroke:#5d82a8,color:#000000
        77["<div style='font-weight: bold'>Versions</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Every change to an<br />experiment's set-up kept in<br />history.jsonl with the hash<br />of its spec.</div>"]
        style 77 fill:#85bbf0,stroke:#5d82a8,color:#000000
        78["<div style='font-weight: bold'>Runs</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Append-only run records in<br />runs.jsonl, the digest<br />derived from them and the<br />filter language of the runs<br />table.</div>"]
        style 78 fill:#85bbf0,stroke:#5d82a8,color:#000000
        79["<div style='font-weight: bold'>Search</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Search over the manifests:<br />plain words with a synonym<br />table and numbers with units<br />matched against card<br />dimensions.</div>"]
        style 79 fill:#85bbf0,stroke:#5d82a8,color:#000000
        80["<div style='font-weight: bold'>Archives</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>.nyx files: a workflow as one<br />deterministic ZIP and the<br />safe way back in.</div>"]
        style 80 fill:#85bbf0,stroke:#5d82a8,color:#000000
        81["<div style='font-weight: bold'>Utilities</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>UTC time without a time<br />crate, time-ordered ids,<br />atomic file writes, hashes<br />and slugs.</div>"]
        style 81 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      117[("<div style='font-weight: bold'>Workspace data</div><div style='font-size: 70%; margin-top: 0px'>[Container: Files and folders]</div><div style='font-size: 80%; margin-top:10px'>workspace.json,<br />workflow.json, history.jsonl,<br />runs.jsonl, specs/<br />(experiment.nyx.json,<br />specs/<hash>.nyx.json), run<br />folders with results.h5, .nyx<br />archives.</div>")]
      style 117 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    80-. "<div>uses flow</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    80-. "<div>uses model</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    80-. "<div>uses util</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    76-. "<div>uses archive</div><div style='font-size: 70%'>[Rust paths]</div>" .->80
    76-. "<div>uses card, model</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    76-. "<div>uses runs</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    76-. "<div>uses util</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    76-. "<div>uses versions</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    74-. "<div>uses util</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    78-. "<div>uses model</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    78-. "<div>uses util</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    79-. "<div>uses model</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    79-. "<div>uses runs</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    79-. "<div>uses workspace</div><div style='font-size: 70%'>[Rust paths]</div>" .->75
    75-. "<div>uses flow</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    75-. "<div>uses model</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    75-. "<div>uses runs</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    75-. "<div>uses util</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    76-. "<div>Reads and writes workflow<br />folders</div><div style='font-size: 70%'>[declared in spec]</div>" .->117
    78-. "<div>Appends runs.jsonl</div><div style='font-size: 70%'>[declared in spec]</div>" .->117

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Data model and card](../4-code/store/model.md) | 2 | 223 | 14, 39 | The workflow manifest, experiments, run digest and the concept card generated from an experiment's spec. |
| [Workspace](../4-code/store/workspace.md) | 2 | 306 | 14 | The folder behind the dashboard: settings, the list of workflows (manifests only) and the sections of the dashboard. |
| [Workflow folder](../4-code/store/flow.md) | 1 | 484 | 14, 59 | One workflow folder: workflow.json, history.jsonl, notes, specs and one folder per experiment. |
| [Versions](../4-code/store/versions.md) | 1 | 200 | 14, 59 | Every change to an experiment's set-up kept in history.jsonl with the hash of its spec. |
| [Runs](../4-code/store/runs.md) | 1 | 448 | 64 | Append-only run records in runs.jsonl, the digest derived from them and the filter language of the runs table. |
| [Search](../4-code/store/search.md) | 1 | 253 | 14 | Search over the manifests: plain words with a synonym table and numbers with units matched against card dimensions. |
| [Archives](../4-code/store/archive.md) | 1 | 703 | 14 | .nyx files: a workflow as one deterministic ZIP and the safe way back in. |
| [Utilities](../4-code/store/util.md) | 1 | 171 | 14 | UTC time without a time crate, time-ordered ids, atomic file writes, hashes and slugs. |
