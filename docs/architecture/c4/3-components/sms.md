# Level 3: components of nyx-sms

Crate `crates/nyx-sms`. Solid mechanics: static, eigenfrequency, PML, Floquet, low-reflecting boundaries, axisymmetric.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_sms.mmd) · [PlantUML](../export/plantuml/structurizr-Components_sms.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-sms"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 46 ["nyx-sms"]
        style 46 fill:#ffffff,stroke:#205723,color:#205723

        47["<div style='font-weight: bold'>Physics adapter</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics implementation of<br />sms: descriptors, validation,<br />solve entry.</div>"]
        style 47 fill:#85bbf0,stroke:#5d82a8,color:#000000
        48["<div style='font-weight: bold'>Typed problem</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The typed sms problem the<br />kernels read, decoded from<br />the envelope.</div>"]
        style 48 fill:#85bbf0,stroke:#5d82a8,color:#000000
        49["<div style='font-weight: bold'>Solid kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Linear elastic solid<br />mechanics on P1/P2 triangles<br />and tetrahedra: static solve<br />and eigenfrequencies.</div>"]
        style 49 fill:#85bbf0,stroke:#5d82a8,color:#000000
        50["<div style='font-weight: bold'>Cross-section kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Stiffness and mass on a 2D<br />mesh with three displacement<br />components (the shared kernel<br />of the elastic cross-section<br />models).</div>"]
        style 50 fill:#85bbf0,stroke:#5d82a8,color:#000000
        51["<div style='font-weight: bold'>Axisymmetric kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Elastic modes of a body of<br />revolution on its half-plane<br />cross-section.</div>"]
        style 51 fill:#85bbf0,stroke:#5d82a8,color:#000000
        52["<div style='font-weight: bold'>Floquet kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Eigenfrequencies of a solid<br />with Floquet-periodic<br />boundaries (phononic<br />crystals).</div>"]
        style 52 fill:#85bbf0,stroke:#5d82a8,color:#000000
        53["<div style='font-weight: bold'>Absorbing-layer kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Eigenfrequencies of elastic<br />bodies with perfectly matched<br />layers.</div>"]
        style 53 fill:#85bbf0,stroke:#5d82a8,color:#000000
        54["<div style='font-weight: bold'>Low-reflecting boundary</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Low-reflecting boundaries for<br />solid mechanics.</div>"]
        style 54 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      34["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 34 fill:#438dd5,stroke:#2e6295,color:#ffffff
      77["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 77 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    34-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->77
    47-. "<div>uses femcore: material,<br />material::spec +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    47-. "<div>uses physics: crate root,<br />decode</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    47-. "<div>uses cross_section</div><div style='font-size: 70%'>[Rust paths]</div>" .->50
    47-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48
    47-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->49
    51-. "<div>uses femcore: fem,<br />linalg::csr +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    51-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    50-. "<div>uses femcore: fem, fields +6</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    50-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    50-. "<div>uses axisym</div><div style='font-size: 70%'>[Rust paths]</div>" .->51
    50-. "<div>uses lrb</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    50-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48
    52-. "<div>uses femcore: linalg::csr,<br />linalg::csr::assembly +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    52-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    52-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48
    52-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->49
    54-. "<div>uses femcore: fem,<br />linalg::csr +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    54-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    54-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48
    54-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->49
    53-. "<div>uses femcore: fem,<br />linalg::csr +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    53-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    53-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48
    53-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->49
    48-. "<div>uses femcore: props, validate</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    48-. "<div>uses physics: decode, error<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    48-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->47
    49-. "<div>uses femcore: fem, fields +9</div><div style='font-size: 70%'>[Rust paths]</div>" .->77
    49-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->34
    49-. "<div>uses floquet</div><div style='font-size: 70%'>[Rust paths]</div>" .->52
    49-. "<div>uses lrb</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    49-. "<div>uses pml</div><div style='font-size: 70%'>[Rust paths]</div>" .->53
    49-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->48

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics adapter](../4-code/sms/adapter.md) | 1 | 264 | 42, 25 | The Physics implementation of sms: descriptors, validation, solve entry. |
| [Typed problem](../4-code/sms/problem.md) | 1 | 314 | 42 | The typed sms problem the kernels read, decoded from the envelope. |
| [Solid kernel](../4-code/sms/solid.md) | 1 | 739 | 42, 9 | Linear elastic solid mechanics on P1/P2 triangles and tetrahedra: static solve and eigenfrequencies. |
| [Cross-section kernel](../4-code/sms/crosssection.md) | 1 | 354 | 43 | Stiffness and mass on a 2D mesh with three displacement components (the shared kernel of the elastic cross-section models). |
| [Axisymmetric kernel](../4-code/sms/axisym.md) | 1 | 232 | 43 | Elastic modes of a body of revolution on its half-plane cross-section. |
| [Floquet kernel](../4-code/sms/floquet.md) | 1 | 128 | 42 | Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals). |
| [Absorbing-layer kernel](../4-code/sms/pml.md) | 1 | 214 | 42 | Eigenfrequencies of elastic bodies with perfectly matched layers. |
| [Low-reflecting boundary](../4-code/sms/lrb.md) | 1 | 296 | 42 | Low-reflecting boundaries for solid mechanics. |
