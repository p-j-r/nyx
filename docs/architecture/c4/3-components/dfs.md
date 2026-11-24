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

      subgraph 48 ["nyx-dfs"]
        style 48 fill:#ffffff,stroke:#205723,color:#205723

        49["<div style='font-weight: bold'>Physics adapter</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics implementation of<br />dfs: descriptors, validation,<br />solve entry.</div>"]
        style 49 fill:#85bbf0,stroke:#5d82a8,color:#000000
        50["<div style='font-weight: bold'>Typed problem</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The typed dfs problem the<br />kernels read, decoded from<br />the ProblemSpec envelope.</div>"]
        style 50 fill:#85bbf0,stroke:#5d82a8,color:#000000
        51["<div style='font-weight: bold'>Heat kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Heat transfer on P1/P2<br />simplices, stationary and<br />transient.</div>"]
        style 51 fill:#85bbf0,stroke:#5d82a8,color:#000000
        52["<div style='font-weight: bold'>Electrostatics kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Electrostatics on P1/P2<br />simplices with<br />fixed-potential,<br />surface-charge, and multi-RHS<br />terminal conditions.</div>"]
        style 52 fill:#85bbf0,stroke:#5d82a8,color:#000000
        53["<div style='font-weight: bold'>Magnetostatics kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Gauge-stabilized London<br />magnetostatics (2D<br />cross-section, planar stream<br />function, mixed 3D A-phi) and<br />inductance matrix extraction.</div>"]
        style 53 fill:#85bbf0,stroke:#5d82a8,color:#000000
        54["<div style='font-weight: bold'>Shared helpers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Per-cell gradients and values<br />shared by the heat and<br />electrostatics kernels.</div>"]
        style 54 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      98["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 98 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->98
    49-. "<div>uses femcore: material::spec,<br />progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    49-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    49-. "<div>uses electro, planar_bem</div><div style='font-size: 70%'>[Rust paths]</div>" .->52
    49-. "<div>uses heat</div><div style='font-size: 70%'>[Rust paths]</div>" .->51
    49-. "<div>uses magneto</div><div style='font-size: 70%'>[Rust paths]</div>" .->53
    49-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->50
    54-. "<div>uses femcore: fields,<br />material +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    54-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    52-. "<div>uses femcore: adapt::electro,<br />adapt::estimate +14</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    52-. "<div>uses physics: cqed, crate<br />root +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    52-. "<div>uses common</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    52-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->50
    51-. "<div>uses femcore: expr, fem +6</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    51-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    51-. "<div>uses common</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    51-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->50
    53-. "<div>uses femcore: fem, fields +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    53-. "<div>uses physics: cqed, crate<br />root +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    53-. "<div>uses common</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    53-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->50
    50-. "<div>uses femcore: material,<br />material::bcs +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    50-. "<div>uses physics: cqed, decode +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    50-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->49

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics adapter](../4-code/dfs/adapter.md) | 1 | 501 | 41, 25, 103 | The Physics implementation of dfs: descriptors, validation, solve entry. |
| [Typed problem](../4-code/dfs/problem.md) | 1 | 606 | 41, 103 | The typed dfs problem the kernels read, decoded from the ProblemSpec envelope. |
| [Heat kernel](../4-code/dfs/heat.md) | 1 | 277 | 41, 9 | Heat transfer on P1/P2 simplices, stationary and transient. |
| [Electrostatics kernel](../4-code/dfs/electro.md) | 3 | 3451 | 41, 9, 103 | Electrostatics on P1/P2 simplices with fixed-potential, surface-charge, and multi-RHS terminal conditions. |
| [Magnetostatics kernel](../4-code/dfs/magneto.md) | 1 | 1129 | 103 | Gauge-stabilized London magnetostatics (2D cross-section, planar stream function, mixed 3D A-phi) and inductance matrix extraction. |
| [Shared helpers](../4-code/dfs/common.md) | 1 | 58 | 41 | Per-cell gradients and values shared by the heat and electrostatics kernels. |
