# Level 3: components of nyx-ems

Crate `crates/nyx-ems`. EM waves: waveguide modes (real and complex), axisymmetric, cavity in 3D.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_ems.mmd) · [PlantUML](../export/plantuml/structurizr-Components_ems.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-ems"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 63 ["nyx-ems"]
        style 63 fill:#ffffff,stroke:#205723,color:#205723

        64["<div style='font-weight: bold'>Physics adapter</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics implementation of<br />ems: descriptors, validation,<br />solve entry.</div>"]
        style 64 fill:#85bbf0,stroke:#5d82a8,color:#000000
        65["<div style='font-weight: bold'>Typed problem</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The typed ems problem the<br />kernels read, decoded from<br />the envelope.</div>"]
        style 65 fill:#85bbf0,stroke:#5d82a8,color:#000000
        66["<div style='font-weight: bold'>Mode kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Full-vector mode analysis of<br />a 2D waveguide cross-section<br />(inhomogeneous media<br />allowed).</div>"]
        style 66 fill:#85bbf0,stroke:#5d82a8,color:#000000
        67["<div style='font-weight: bold'>Complex-mode kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Complex-valued mode analysis:<br />lossy media and perfectly<br />matched layers.</div>"]
        style 67 fill:#85bbf0,stroke:#5d82a8,color:#000000
        68["<div style='font-weight: bold'>Axisymmetric kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Electromagnetic eigenmodes of<br />a body of revolution (rings,<br />disks, coaxial cavities).</div>"]
        style 68 fill:#85bbf0,stroke:#5d82a8,color:#000000
        69["<div style='font-weight: bold'>3D cavity kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Electromagnetic eigenmodes of<br />3D cavities and resonators.</div>"]
        style 69 fill:#85bbf0,stroke:#5d82a8,color:#000000
        70["<div style='font-weight: bold'>Optical mode sets</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The mode solvers' results as<br />an OpticalSet, with the<br />material data needed to<br />evaluate their fields.</div>"]
        style 70 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      92["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 92 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->92
    64-. "<div>uses femcore: material,<br />material::spec +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    64-. "<div>uses physics: crate root,<br />decode</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    64-. "<div>uses axisym</div><div style='font-size: 70%'>[Rust paths]</div>" .->68
    64-. "<div>uses cavity</div><div style='font-size: 70%'>[Rust paths]</div>" .->69
    64-. "<div>uses crate root, modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->66
    64-. "<div>uses modes_complex</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    64-. "<div>uses optical</div><div style='font-size: 70%'>[Rust paths]</div>" .->70
    64-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    68-. "<div>uses femcore: edge, fem +7</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    68-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    68-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->66
    68-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    69-. "<div>uses femcore: fem, fields +11</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    69-. "<div>uses physics: crate root,<br />frame +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    69-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->66
    69-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    66-. "<div>uses femcore: edge, fem +9</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    66-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    66-. "<div>uses modes_complex</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    66-. "<div>uses optical</div><div style='font-size: 70%'>[Rust paths]</div>" .->70
    66-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    67-. "<div>uses femcore: edge, fem +7</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    67-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    67-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->66
    67-. "<div>uses optical</div><div style='font-size: 70%'>[Rust paths]</div>" .->70
    67-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    70-. "<div>uses femcore:<br />linalg::csr::assembly,<br />optical +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    70-. "<div>uses physics: crate root,<br />frame</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    70-. "<div>uses modes</div><div style='font-size: 70%'>[Rust paths]</div>" .->66
    70-. "<div>uses modes_complex</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    70-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->65
    65-. "<div>uses femcore: props</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    65-. "<div>uses physics: decode, error<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    65-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->64

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics adapter](../4-code/ems/adapter.md) | 1 | 295 | 45, 25 | The Physics implementation of ems: descriptors, validation, solve entry. |
| [Typed problem](../4-code/ems/problem.md) | 1 | 189 | 45 | The typed ems problem the kernels read, decoded from the envelope. |
| [Mode kernel](../4-code/ems/modes.md) | 1 | 440 | 45, 9 | Full-vector mode analysis of a 2D waveguide cross-section (inhomogeneous media allowed). |
| [Complex-mode kernel](../4-code/ems/modescx.md) | 1 | 466 | 45 | Complex-valued mode analysis: lossy media and perfectly matched layers. |
| [Axisymmetric kernel](../4-code/ems/axisym.md) | 1 | 469 | 46 | Electromagnetic eigenmodes of a body of revolution (rings, disks, coaxial cavities). |
| [3D cavity kernel](../4-code/ems/cavity.md) | 1 | 791 | 45 | Electromagnetic eigenmodes of 3D cavities and resonators. |
| [Optical mode sets](../4-code/ems/optical.md) | 1 | 99 | 45, 21 | The mode solvers' results as an OpticalSet, with the material data needed to evaluate their fields. |
