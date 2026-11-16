# Level 3: components of nyx-dfs

Crate `crates/nyx-dfs`. Diffusion physics: heat and electrostatics. Typed problem, three-stage kernels, descriptors, goldens.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_dfs.mmd) · [PlantUML](../export/plantuml/structurizr-Components_dfs.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-dfs"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 40 ["nyx-dfs"]
        style 40 fill:#ffffff,stroke:#205723,color:#205723

        41["<div style='font-weight: bold'>Physics adapter</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics implementation of<br />dfs: descriptors, validation,<br />solve entry.</div>"]
        style 41 fill:#85bbf0,stroke:#5d82a8,color:#000000
        42["<div style='font-weight: bold'>Typed problem</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The typed dfs problem the<br />kernels read, decoded from<br />the ProblemSpec envelope.</div>"]
        style 42 fill:#85bbf0,stroke:#5d82a8,color:#000000
        43["<div style='font-weight: bold'>Heat kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Heat transfer on P1/P2<br />simplices, stationary and<br />transient.</div>"]
        style 43 fill:#85bbf0,stroke:#5d82a8,color:#000000
        44["<div style='font-weight: bold'>Electrostatics kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Electrostatics on P1/P2<br />simplices with<br />fixed-potential and<br />surface-charge conditions.</div>"]
        style 44 fill:#85bbf0,stroke:#5d82a8,color:#000000
        45["<div style='font-weight: bold'>Shared helpers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Per-cell gradients and values<br />shared by the heat and<br />electrostatics kernels.</div>"]
        style 45 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      34["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 34 fill:#438dd5,stroke:#2e6295,color:#ffffff
      77["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 77 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    34-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->77
    41-. "<div>uses femcore: material::spec,<br />progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    41-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    41-. "<div>uses electro</div><div style='font-size: 70%'>[Rust paths]</div>" .->44
    41-. "<div>uses heat</div><div style='font-size: 70%'>[Rust paths]</div>" .->43
    41-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->42
    45-. "<div>uses femcore: fields,<br />material +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    45-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    44-. "<div>uses femcore: fem, fields +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    44-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    44-. "<div>uses common</div><div style='font-size: 70%'>[Rust paths]</div>" .->45
    44-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->42
    43-. "<div>uses femcore: expr, fem +6</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    43-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    43-. "<div>uses common</div><div style='font-size: 70%'>[Rust paths]</div>" .->45
    43-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->42
    42-. "<div>uses femcore: props</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    42-. "<div>uses physics: decode, error<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    42-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->41

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics adapter](../4-code/dfs/adapter.md) | 1 | 217 | 41, 25 | The Physics implementation of dfs: descriptors, validation, solve entry. |
| [Typed problem](../4-code/dfs/problem.md) | 1 | 195 | 41 | The typed dfs problem the kernels read, decoded from the ProblemSpec envelope. |
| [Heat kernel](../4-code/dfs/heat.md) | 1 | 273 | 41, 9 | Heat transfer on P1/P2 simplices, stationary and transient. |
| [Electrostatics kernel](../4-code/dfs/electro.md) | 1 | 136 | 41, 9 | Electrostatics on P1/P2 simplices with fixed-potential and surface-charge conditions. |
| [Shared helpers](../4-code/dfs/common.md) | 1 | 56 | 41 | Per-cell gradients and values shared by the heat and electrostatics kernels. |
