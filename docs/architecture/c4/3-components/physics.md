# Level 3: components of nyx-physics

Crate `crates/nyx-physics`. The Physics trait, the problem envelope (ProblemSpec), descriptors, SolveResult, the registry type, validation, JSON Schemas.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_physics.mmd) · [PlantUML](../export/plantuml/structurizr-Components_physics.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-physics"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 35 ["nyx-physics"]
        style 35 fill:#ffffff,stroke:#2e6295,color:#2e6295

        36["<div style='font-weight: bold'>Physics trait and registry</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait (everything<br />the rest of nyx may ask of a<br />physics), the Registry of the<br />physics a build contains, and<br />the plain-data error every<br />method returns.</div>"]
        style 36 fill:#85bbf0,stroke:#5d82a8,color:#000000
        37["<div style='font-weight: bold'>Problem envelope</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The one shape every physics<br />is given (ProblemSpec), its<br />decoding into typed<br />parameters with descriptor<br />defaults and name resolution,<br />and its structural<br />validation.</div>"]
        style 37 fill:#85bbf0,stroke:#5d82a8,color:#000000
        38["<div style='font-weight: bold'>Descriptors and schemas</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>What a physics tells nyx<br />about itself (studies,<br />conditions, parameters) and<br />the JSON Schemas generated<br />from it for the window, the<br />MCP server and the<br />documentation.</div>"]
        style 38 fill:#85bbf0,stroke:#5d82a8,color:#000000
        39["<div style='font-weight: bold'>Result and frame</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>SolveResult as plain<br />serialisable data, and the<br />frame a result states for the<br />analysis of its vector<br />fields.</div>"]
        style 39 fill:#85bbf0,stroke:#5d82a8,color:#000000
        40["<div style='font-weight: bold'>Test support</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The golden-file harness and<br />the contract-test helpers<br />shared by every physics<br />crate.</div>"]
        style 40 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      78["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 78 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    38-. "<div>uses femcore: meshing::sizing</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    38-. "<div>uses physics</div><div style='font-size: 70%'>[Rust paths]</div>" .->36
    37-. "<div>uses femcore: material,<br />material::elastic +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    37-. "<div>uses desc</div><div style='font-size: 70%'>[Rust paths]</div>" .->38
    37-. "<div>uses result</div><div style='font-size: 70%'>[Rust paths]</div>" .->39
    37-. "<div>uses error, physics</div><div style='font-size: 70%'>[Rust paths]</div>" .->36
    39-. "<div>uses femcore: adapt, fields</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    40-. "<div>uses femcore: fields</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    40-. "<div>uses schema</div><div style='font-size: 70%'>[Rust paths]</div>" .->38
    40-. "<div>uses spec</div><div style='font-size: 70%'>[Rust paths]</div>" .->37
    40-. "<div>uses result</div><div style='font-size: 70%'>[Rust paths]</div>" .->39
    40-. "<div>uses crate root, physics</div><div style='font-size: 70%'>[Rust paths]</div>" .->36
    36-. "<div>uses femcore: material::spec,<br />mesh +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->78
    36-. "<div>uses crate root, desc</div><div style='font-size: 70%'>[Rust paths]</div>" .->38
    36-. "<div>uses spec, validate</div><div style='font-size: 70%'>[Rust paths]</div>" .->37
    36-. "<div>uses frame, result</div><div style='font-size: 70%'>[Rust paths]</div>" .->39

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Physics trait and registry](../4-code/physics/trait.md) | 4 | 294 | 25 | The Physics trait (everything the rest of nyx may ask of a physics), the Registry of the physics a build contains, and the plain-data error every method returns. |
| [Problem envelope](../4-code/physics/envelope.md) | 3 | 534 | 25, 41, 63 | The one shape every physics is given (ProblemSpec), its decoding into typed parameters with descriptor defaults and name resolution, and its structural validation. |
| [Descriptors and schemas](../4-code/physics/descriptors.md) | 2 | 332 | 25 | What a physics tells nyx about itself (studies, conditions, parameters) and the JSON Schemas generated from it for the window, the MCP server and the documentation. |
| [Result and frame](../4-code/physics/result.md) | 2 | 286 | 25, 56 | SolveResult as plain serialisable data, and the frame a result states for the analysis of its vector fields. |
| [Test support](../4-code/physics/testsupport.md) | 2 | 229 | 10, 19 | The golden-file harness and the contract-test helpers shared by every physics crate. |
