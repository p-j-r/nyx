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

      subgraph 117 ["geom"]
        style 117 fill:#ffffff,stroke:#2e6295,color:#2e6295

        118["<div style='font-weight: bold'>Feature list</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The serialisable feature list<br />(single source of truth for<br />code and GUI), its operations<br />and validation, parameter<br />expressions, printing as Rust<br />builder calls, and the error<br />type.</div>"]
        style 118 fill:#85bbf0,stroke:#5d82a8,color:#000000
        119["<div style='font-weight: bold'>Evaluation</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Evaluator turns a feature<br />list into an evaluated model<br />at a chosen tier with a<br />cache; scenes instantiate<br />profiles with transforms and<br />patterns; change detection<br />between evaluations for<br />sweeps.</div>"]
        style 119 fill:#85bbf0,stroke:#5d82a8,color:#000000
        120["<div style='font-weight: bold'>2D kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Primitives, planar curves,<br />the planar arrangement of all<br />boundary pieces and the<br />evaluated region complex.</div>"]
        style 120 fill:#85bbf0,stroke:#5d82a8,color:#000000
        121["<div style='font-weight: bold'>3D kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Prismatic models as a layer<br />stack, general planar-faced<br />solids, and the edges and<br />vertices of a 3D model.</div>"]
        style 121 fill:#85bbf0,stroke:#5d82a8,color:#000000
        122["<div style='font-weight: bold'>Queries and selectors</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Exact queries on an evaluated<br />model, serialisable<br />selectors, the LabelMap<br />between names and the integer<br />labels solvers consume, and<br />simulation operations as<br />data.</div>"]
        style 122 fill:#85bbf0,stroke:#5d82a8,color:#000000
        123["<div style='font-weight: bold'>Display tessellation</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Tessellation for the viewport<br />and picking only, never for<br />meshing.</div>"]
        style 123 fill:#85bbf0,stroke:#5d82a8,color:#000000
        124["<div style='font-weight: bold'>Math helpers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Small vector, box and<br />similarity-transform helpers.</div>"]
        style 124 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    119-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->118
    119-. "<div>uses arrange, complex +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->120
    119-. "<div>uses planar, stack</div><div style='font-size: 70%'>[Rust paths]</div>" .->121
    119-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    119-. "<div>uses sim</div><div style='font-size: 70%'>[Rust paths]</div>" .->122
    118-. "<div>uses select</div><div style='font-size: 70%'>[Rust paths]</div>" .->122
    120-. "<div>uses scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->119
    120-. "<div>uses error, model</div><div style='font-size: 70%'>[Rust paths]</div>" .->118
    120-. "<div>uses stack</div><div style='font-size: 70%'>[Rust paths]</div>" .->121
    120-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    121-. "<div>uses eval, scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->119
    121-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->118
    121-. "<div>uses arrange, complex +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->120
    121-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    121-. "<div>uses tess</div><div style='font-size: 70%'>[Rust paths]</div>" .->123
    122-. "<div>uses eval, scene</div><div style='font-size: 70%'>[Rust paths]</div>" .->119
    122-. "<div>uses error, expr +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->118
    122-. "<div>uses complex, curve</div><div style='font-size: 70%'>[Rust paths]</div>" .->120
    122-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    123-. "<div>uses eval</div><div style='font-size: 70%'>[Rust paths]</div>" .->119
    123-. "<div>uses complex</div><div style='font-size: 70%'>[Rust paths]</div>" .->120
    123-. "<div>uses math</div><div style='font-size: 70%'>[Rust paths]</div>" .->124
    123-. "<div>uses query</div><div style='font-size: 70%'>[Rust paths]</div>" .->122

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Feature list](../4-code/geom/featurelist.md) | 5 | 1410 | 2, 3, 63 | The serialisable feature list (single source of truth for code and GUI), its operations and validation, parameter expressions, printing as Rust builder calls, and the error type. |
| [Evaluation](../4-code/geom/evaluation.md) | 3 | 1306 | 2, 6 | The Evaluator turns a feature list into an evaluated model at a chosen tier with a cache; scenes instantiate profiles with transforms and patterns; change detection between evaluations for sweeps. |
| [2D kernel](../4-code/geom/kernel2d.md) | 4 | 3025 | 3, 18 | Primitives, planar curves, the planar arrangement of all boundary pieces and the evaluated region complex. |
| [3D kernel](../4-code/geom/kernel3d.md) | 3 | 2202 | 3, 1 | Prismatic models as a layer stack, general planar-faced solids, and the edges and vertices of a 3D model. |
| [Queries and selectors](../4-code/geom/queries.md) | 4 | 1905 | 2, 3, 7, 18 | Exact queries on an evaluated model, serialisable selectors, the LabelMap between names and the integer labels solvers consume, and simulation operations as data. |
| [Display tessellation](../4-code/geom/tessellation.md) | 1 | 280 | 3 | Tessellation for the viewport and picking only, never for meshing. |
| [Math helpers](../4-code/geom/math.md) | 1 | 340 | 3 | Small vector, box and similarity-transform helpers. |
