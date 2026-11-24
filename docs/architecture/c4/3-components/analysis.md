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

      subgraph 75 ["nyx-analysis"]
        style 75 fill:#ffffff,stroke:#2e6295,color:#2e6295

        76["<div style='font-weight: bold'>Analysis API</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The entry points an engine<br />calls: an analysis is a<br />function of a SolveResult and<br />its mesh.</div>"]
        style 76 fill:#85bbf0,stroke:#5d82a8,color:#000000
        77["<div style='font-weight: bold'>Capacitance and inductance matrices</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Maxwell/mutual capacitance<br />and magnetostatic inductance<br />matrix extraction from<br />multi-RHS dfs results.</div>"]
        style 77 fill:#85bbf0,stroke:#5d82a8,color:#000000
        78["<div style='font-weight: bold'>cQED quantization and TLS participation</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>KICQ energy participation<br />ratios, normal-Fock and<br />hierarchical Craig-Bampton<br />Josephson quantization, and<br />edge-regularized 3D TLS<br />participation ratios.</div>"]
        style 78 fill:#85bbf0,stroke:#5d82a8,color:#000000
        79["<div style='font-weight: bold'>Field samples</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The modes of a SolveResult as<br />weighted samples of the<br />field, which is what the<br />polarization analysis<br />integrates.</div>"]
        style 79 fill:#85bbf0,stroke:#5d82a8,color:#000000
        80["<div style='font-weight: bold'>Polarization</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Stokes parameters, Poincare<br />sphere and the classification<br />of an optical mode.</div>"]
        style 80 fill:#85bbf0,stroke:#5d82a8,color:#000000
        81["<div style='font-weight: bold'>Mode report</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>One report per mode, a table,<br />the Poincare picture and the<br />degenerate groups.</div>"]
        style 81 fill:#85bbf0,stroke:#5d82a8,color:#000000
        82["<div style='font-weight: bold'>Mode figures</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The cross-section of a<br />result, the phase, sign and<br />scale of a mode, and the<br />figure of its modes (items<br />drawn as SVG or by a small<br />software rasteriser as PNG);<br />no physics is named.</div>"]
        style 82 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      98["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 98 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->98
    76-. "<div>uses physics: crate root,<br />error +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    76-. "<div>uses capacitance</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    76-. "<div>uses tls</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    76-. "<div>uses modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    76-. "<div>uses samples</div><div style='font-size: 70%'>[Rust paths]</div>" .->79
    77-. "<div>uses femcore: fem,<br />linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    77-. "<div>uses physics: cqed, crate<br />root +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    77-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    77-. "<div>uses tls</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    78-. "<div>uses femcore: crate root, fem<br />+5</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    78-. "<div>uses physics: cqed, crate<br />root +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    78-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    82-. "<div>uses femcore: fieldexpr,<br />formula +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    82-. "<div>uses physics: crate root,<br />error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    81-. "<div>uses polarization</div><div style='font-size: 70%'>[Rust paths]</div>" .->80
    79-. "<div>uses femcore: fem, fields +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    79-. "<div>uses physics: crate root,<br />error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    79-. "<div>uses modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->81
    79-. "<div>uses polarization</div><div style='font-size: 70%'>[Rust paths]</div>" .->80

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Analysis API](../4-code/analysis/api.md) | 1 | 127 | 56, 103, 105 | The entry points an engine calls: an analysis is a function of a SolveResult and its mesh. |
| [Capacitance and inductance matrices](../4-code/analysis/capacitance.md) | 1 | 1259 | 103 | Maxwell/mutual capacitance and magnetostatic inductance matrix extraction from multi-RHS dfs results. |
| [cQED quantization and TLS participation](../4-code/analysis/cqed.md) | 4 | 6980 | 105 | KICQ energy participation ratios, normal-Fock and hierarchical Craig-Bampton Josephson quantization, and edge-regularized 3D TLS participation ratios. |
| [Field samples](../4-code/analysis/samples.md) | 1 | 161 | 56 | The modes of a SolveResult as weighted samples of the field, which is what the polarization analysis integrates. |
| [Polarization](../4-code/analysis/polarization.md) | 1 | 938 | 56, 21 | Stokes parameters, Poincare sphere and the classification of an optical mode. |
| [Mode report](../4-code/analysis/modes.md) | 1 | 622 | 56 | One report per mode, a table, the Poincare picture and the degenerate groups. |
| [Mode figures](../4-code/analysis/figures.md) | 6 | 5000 | 99, 56 | The cross-section of a result, the phase, sign and scale of a mode, and the figure of its modes (items drawn as SVG or by a small software rasteriser as PNG); no physics is named. |
