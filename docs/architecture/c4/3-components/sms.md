# Level 3: components of nyx-sms

Crate `crates/nyx-sms`. Solid mechanics: static, eigenfrequency, PML, Floquet, low-reflecting boundaries, axisymmetric eigenfrequency study (bodies of revolution, axis conditions, cylindrical layers).

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_sms.mmd) · [PlantUML](../export/plantuml/structurizr-Components_sms.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-sms"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 55 ["nyx-sms"]
        style 55 fill:#ffffff,stroke:#205723,color:#205723

        56["<div style='font-weight: bold'>Physics adapter</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics implementation of<br />sms with four studies<br />(static, eigenfrequency,<br />axisymmetric eigenfrequency,<br />waveguide modes):<br />descriptors, validation,<br />solve entry.</div>"]
        style 56 fill:#85bbf0,stroke:#5d82a8,color:#000000
        57["<div style='font-weight: bold'>Typed problem</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The typed sms problem the<br />kernels read, decoded from<br />the envelope.</div>"]
        style 57 fill:#85bbf0,stroke:#5d82a8,color:#000000
        58["<div style='font-weight: bold'>Solid kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Linear elastic solid<br />mechanics on P1/P2 triangles<br />and tetrahedra: static solve<br />and eigenfrequencies.</div>"]
        style 58 fill:#85bbf0,stroke:#5d82a8,color:#000000
        59["<div style='font-weight: bold'>Cross-section study kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Cross-section study kernel<br />(`axisymmetric_eigenfrequency`<br />on the (r, z) half plane and<br />2.5D `modes` on the (x, z)<br />cross-section): constraints,<br />axis and periodic conditions,<br />absorbing layers, complex<br />eigensolve and fields; every<br />mode carries numeric<br />frequency parts, a layer<br />fraction and its energy<br />shares per component and<br />domain.</div>"]
        style 59 fill:#85bbf0,stroke:#5d82a8,color:#000000
        60["<div style='font-weight: bold'>Axisymmetric kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Elastic modes of a body of<br />revolution on its half-plane<br />cross-section, with the axis<br />regularity conditions, the<br />cylindrical PML stretch and<br />the energy weight of a mode<br />per component and class of<br />cells.</div>"]
        style 60 fill:#85bbf0,stroke:#5d82a8,color:#000000
        61["<div style='font-weight: bold'>Waveguide kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Elastic modes of a 2.5D<br />waveguide on its 2D (x, z)<br />cross-section (`exp(-i q<br />y)`), with full 3D<br />anisotropic crystal<br />orientation, planar PML<br />coordinate stretching and the<br />planar kinetic-energy weight<br />of a mode per component and<br />class of cells.</div>"]
        style 61 fill:#85bbf0,stroke:#5d82a8,color:#000000
        62["<div style='font-weight: bold'>Floquet kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Eigenfrequencies of a solid<br />with Floquet-periodic<br />boundaries (phononic<br />crystals).</div>"]
        style 62 fill:#85bbf0,stroke:#5d82a8,color:#000000
        63["<div style='font-weight: bold'>Absorbing-layer kernel</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Eigenfrequencies of elastic<br />bodies with perfectly matched<br />layers and optional<br />Floquet-periodic boundaries.</div>"]
        style 63 fill:#85bbf0,stroke:#5d82a8,color:#000000
        64["<div style='font-weight: bold'>Low-reflecting boundary</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Low-reflecting boundaries for<br />solid mechanics.</div>"]
        style 64 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      41["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 41 fill:#438dd5,stroke:#2e6295,color:#ffffff
      98["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 98 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    41-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->98
    56-. "<div>uses femcore: material,<br />material::spec +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    56-. "<div>uses physics: crate root,<br />decode</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    56-. "<div>uses cross_section</div><div style='font-size: 70%'>[Rust paths]</div>" .->59
    56-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    56-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->58
    60-. "<div>uses femcore: fem,<br />linalg::csr +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    60-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    59-. "<div>uses femcore: fem, fields +8</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    59-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    59-. "<div>uses axisym</div><div style='font-size: 70%'>[Rust paths]</div>" .->60
    59-. "<div>uses lrb</div><div style='font-size: 70%'>[Rust paths]</div>" .->64
    59-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    59-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->58
    59-. "<div>uses waveguide</div><div style='font-size: 70%'>[Rust paths]</div>" .->61
    62-. "<div>uses femcore: linalg::csr,<br />linalg::csr::assembly +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    62-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    62-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    62-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->58
    64-. "<div>uses femcore: fem,<br />linalg::csr +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    64-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    64-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    64-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->58
    63-. "<div>uses femcore: fem,<br />linalg::csr +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    63-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    63-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    63-. "<div>uses solid</div><div style='font-size: 70%'>[Rust paths]</div>" .->58
    57-. "<div>uses femcore: props, validate</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    57-. "<div>uses physics: decode, error<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    57-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->56
    58-. "<div>uses femcore: fem, fields +9</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    58-. "<div>uses physics: crate root,<br />result</div><div style='font-size: 70%'>[Rust paths]</div>" .->41
    58-. "<div>uses cross_section</div><div style='font-size: 70%'>[Rust paths]</div>" .->59
    58-. "<div>uses floquet</div><div style='font-size: 70%'>[Rust paths]</div>" .->62
    58-. "<div>uses lrb</div><div style='font-size: 70%'>[Rust paths]</div>" .->64
    58-. "<div>uses pml</div><div style='font-size: 70%'>[Rust paths]</div>" .->63
    58-. "<div>uses problem</div><div style='font-size: 70%'>[Rust paths]</div>" .->57
    61-. "<div>uses femcore: fem,<br />linalg::csr +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    61-. "<div>uses physics: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->41

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics adapter](../4-code/sms/adapter.md) | 1 | 352 | 42, 25, 68, 101 | The Physics implementation of sms with four studies (static, eigenfrequency, axisymmetric eigenfrequency, waveguide modes): descriptors, validation, solve entry. |
| [Typed problem](../4-code/sms/problem.md) | 1 | 417 | 42, 68, 101 | The typed sms problem the kernels read, decoded from the envelope. |
| [Solid kernel](../4-code/sms/solid.md) | 1 | 768 | 42, 9 | Linear elastic solid mechanics on P1/P2 triangles and tetrahedra: static solve and eigenfrequencies. |
| [Cross-section study kernel](../4-code/sms/crosssection.md) | 1 | 591 | 43, 68, 69, 80, 101 | Cross-section study kernel (`axisymmetric_eigenfrequency` on the (r, z) half plane and 2.5D `modes` on the (x, z) cross-section): constraints, axis and periodic conditions, absorbing layers, complex eigensolve and fields; every mode carries numeric frequency parts, a layer fraction and its energy shares per component and domain. |
| [Axisymmetric kernel](../4-code/sms/axisym.md) | 1 | 390 | 43, 69, 80 | Elastic modes of a body of revolution on its half-plane cross-section, with the axis regularity conditions, the cylindrical PML stretch and the energy weight of a mode per component and class of cells. |
| [Waveguide kernel](../4-code/sms/waveguide.md) | 1 | 158 | 101 | Elastic modes of a 2.5D waveguide on its 2D (x, z) cross-section (`exp(-i q y)`), with full 3D anisotropic crystal orientation, planar PML coordinate stretching and the planar kinetic-energy weight of a mode per component and class of cells. |
| [Floquet kernel](../4-code/sms/floquet.md) | 1 | 130 | 42 | Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals). |
| [Absorbing-layer kernel](../4-code/sms/pml.md) | 1 | 266 | 42, 101 | Eigenfrequencies of elastic bodies with perfectly matched layers and optional Floquet-periodic boundaries. |
| [Low-reflecting boundary](../4-code/sms/lrb.md) | 1 | 298 | 42 | Low-reflecting boundaries for solid mechanics. |
