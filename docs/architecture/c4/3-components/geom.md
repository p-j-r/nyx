# Level 3: components of geom

Crate `crates/geom`. Lean CSG-over-B-rep geometry kernel: feature lists, named regions and boundaries, selectors, exact projection oracles. Depends on serde and optional rayon only.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_geom.mmd) · [PlantUML](../export/plantuml/structurizr-Components_geom.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - geom"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 123 ["geom"]
        style 123 fill:#ffffff,stroke:#2e6295,color:#2e6295

        124["<div style='font-weight: bold'>Feature list</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The serialisable feature list<br />(single source of truth for<br />code and GUI), its operations<br />and validation, parameter<br />expressions, printing as Rust<br />builder calls, and the error<br />type.</div>"]
        style 124 fill:#85bbf0,stroke:#5d82a8,color:#000000
        125["<div style='font-weight: bold'>Evaluation</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Evaluator turns a feature<br />list into an evaluated model<br />at a chosen tier with a<br />cache; scenes instantiate<br />profiles with transforms and<br />patterns; change detection<br />between evaluations for<br />sweeps.</div>"]
        style 125 fill:#85bbf0,stroke:#5d82a8,color:#000000
        126["<div style='font-weight: bold'>2D kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Primitives, planar curves,<br />the planar arrangement of all<br />boundary pieces, the<br />evaluated region complex, and<br />batch polygon ring cleanup<br />with shared-corner-pinned RDP<br />simplification.</div>"]
        style 126 fill:#85bbf0,stroke:#5d82a8,color:#000000
        127["<div style='font-weight: bold'>3D kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Prismatic models as a layer<br />stack, general planar-faced<br />solids, and the edges and<br />vertices of a 3D model.</div>"]
        style 127 fill:#85bbf0,stroke:#5d82a8,color:#000000
        128["<div style='font-weight: bold'>Queries and selectors</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Exact queries on an evaluated<br />model, serialisable<br />selectors, the LabelMap<br />between names and the integer<br />labels solvers consume, and<br />simulation operations as<br />data.</div>"]
        style 128 fill:#85bbf0,stroke:#5d82a8,color:#000000
        129["<div style='font-weight: bold'>Display tessellation</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Tessellation for the viewport<br />and picking only, never for<br />meshing.</div>"]
        style 129 fill:#85bbf0,stroke:#5d82a8,color:#000000
        130["<div style='font-weight: bold'>Math helpers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Small vector, box and<br />similarity-transform helpers.</div>"]
        style 130 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    125-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    125-. "<div>uses arrange, complex +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->126
    125-. "<div>uses planar, stack</div><div style='font-size: 70%'>[Rust paths]</div>" .->127
    125-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->130
    125-. "<div>uses sim</div><div style='font-size: 70%'>[Rust paths]</div>" .->128
    124-. "<div>uses select</div><div style='font-size: 70%'>[Rust paths]</div>" .->128
    126-. "<div>uses scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->125
    126-. "<div>uses error, model</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    126-. "<div>uses stack</div><div style='font-size: 70%'>[Rust paths]</div>" .->127
    126-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->130
    127-. "<div>uses eval, scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->125
    127-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    127-. "<div>uses arrange, complex +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->126
    127-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->130
    127-. "<div>uses tess</div><div style='font-size: 70%'>[Rust paths]</div>" .->129
    128-. "<div>uses eval, scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->125
    128-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    128-. "<div>uses complex, curve</div><div style='font-size: 70%'>[Rust paths]</div>" .->126
    128-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->130
    129-. "<div>uses eval</div><div style='font-size: 70%'>[Rust paths]</div>" .->125
    129-. "<div>uses complex</div><div style='font-size: 70%'>[Rust paths]</div>" .->126
    129-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->130
    129-. "<div>uses query</div><div style='font-size: 70%'>[Rust paths]</div>" .->128

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Feature list](../4-code/geom/featurelist.md) | 5 | 1468 | 2, 3, 63 | The serialisable feature list (single source of truth for code and GUI), its operations and validation, parameter expressions, printing as Rust builder calls, and the error type. |
| [Evaluation](../4-code/geom/evaluation.md) | 3 | 1456 | 2, 6 | The Evaluator turns a feature list into an evaluated model at a chosen tier with a cache; scenes instantiate profiles with transforms and patterns; change detection between evaluations for sweeps. |
| [2D kernel](../4-code/geom/kernel2d.md) | 5 | 3323 | 3, 18, 102 | Primitives, planar curves, the planar arrangement of all boundary pieces, the evaluated region complex, and batch polygon ring cleanup with shared-corner-pinned RDP simplification. |
| [3D kernel](../4-code/geom/kernel3d.md) | 3 | 2245 | 3, 1 | Prismatic models as a layer stack, general planar-faced solids, and the edges and vertices of a 3D model. |
| [Queries and selectors](../4-code/geom/queries.md) | 4 | 1915 | 2, 3, 7, 18 | Exact queries on an evaluated model, serialisable selectors, the LabelMap between names and the integer labels solvers consume, and simulation operations as data. |
| [Display tessellation](../4-code/geom/tessellation.md) | 1 | 280 | 3 | Tessellation for the viewport and picking only, never for meshing. |
| [Math helpers](../4-code/geom/math.md) | 1 | 340 | 3 | Small vector, box and similarity-transform helpers. |
