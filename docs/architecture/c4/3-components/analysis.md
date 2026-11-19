# Level 3: components of nyx-analysis

Crate `crates/nyx-analysis`. Analyses of results: the polarization (Stokes parameters, Poincare sphere, classes, degenerate groups) as a function of a SolveResult and its mesh.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_analysis.mmd) · [PlantUML](../export/plantuml/structurizr-Components_analysis.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-analysis"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 70 ["nyx-analysis"]
        style 70 fill:#ffffff,stroke:#2e6295,color:#2e6295

        71["<div style='font-weight: bold'>Analysis API</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The entry points an engine<br />calls: an analysis is a<br />function of a SolveResult and<br />its mesh.</div>"]
        style 71 fill:#85bbf0,stroke:#5d82a8,color:#000000
        72["<div style='font-weight: bold'>Field samples</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The modes of a SolveResult as<br />weighted samples of the<br />field, which is what the<br />polarization analysis<br />integrates.</div>"]
        style 72 fill:#85bbf0,stroke:#5d82a8,color:#000000
        73["<div style='font-weight: bold'>Polarization</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Stokes parameters, Poincare<br />sphere and the classification<br />of an optical mode.</div>"]
        style 73 fill:#85bbf0,stroke:#5d82a8,color:#000000
        74["<div style='font-weight: bold'>Mode report</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>One report per mode, a table,<br />the Poincare picture and the<br />degenerate groups.</div>"]
        style 74 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      84["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 84 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->84
    71-. "<div>uses physics: crate root,<br />error +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    71-. "<div>uses modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    71-. "<div>uses samples</div><div style='font-size: 70%'>[Rust paths]</div>" .->72
    74-. "<div>uses polarization</div><div style='font-size: 70%'>[Rust paths]</div>" .->73
    72-. "<div>uses femcore: fem, fields +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->84
    72-. "<div>uses physics: crate root,<br />error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    72-. "<div>uses modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->74
    72-. "<div>uses polarization</div><div style='font-size: 70%'>[Rust paths]</div>" .->73

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Analysis API](../4-code/analysis/api.md) | 1 | 91 | 56 | The entry points an engine calls: an analysis is a function of a SolveResult and its mesh. |
| [Field samples](../4-code/analysis/samples.md) | 1 | 151 | 56 | The modes of a SolveResult as weighted samples of the field, which is what the polarization analysis integrates. |
| [Polarization](../4-code/analysis/polarization.md) | 1 | 936 | 56, 21 | Stokes parameters, Poincare sphere and the classification of an optical mode. |
| [Mode report](../4-code/analysis/modes.md) | 1 | 620 | 56 | One report per mode, a table, the Poincare picture and the degenerate groups. |
