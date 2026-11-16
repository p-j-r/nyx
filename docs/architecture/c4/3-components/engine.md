# Level 3: components of nyx-engine

Crate `crates/nyx-engine`. The Engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them, and the binding to a workflow. Every action is a named command with a JSON Schema.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_engine.mmd) · [PlantUML](../export/plantuml/structurizr-Components_engine.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - nyx-engine"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 22 ["nyx-engine"]
        style 22 fill:#ffffff,stroke:#2e6295,color:#2e6295

        23["<div style='font-weight: bold'>Commands</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The Engine and its catalog:<br />every action by name with<br />JSON arguments, JSON Schema<br />and JSON reply; the<br />dispatcher; the guide cards a<br />host reads on demand. The<br />window and the MCP server use<br />the same entries.</div>"]
        style 23 fill:#85bbf0,stroke:#5d82a8,color:#000000
        24["<div style='font-weight: bold'>Project</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The saved state: feature<br />list, mesh settings, one<br />ProblemSpec per physics,<br />project materials. The format<br />key is checked at the file<br />boundary.</div>"]
        style 24 fill:#85bbf0,stroke:#5d82a8,color:#000000
        25["<div style='font-weight: bold'>Geometry operations</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Every edit of a geom feature<br />list as a function that does<br />it or says why not, and the<br />descriptors of the geom<br />operations (they live here<br />because geom depends on serde<br />only).</div>"]
        style 25 fill:#85bbf0,stroke:#5d82a8,color:#000000
        26["<div style='font-weight: bold'>Job runner</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Mesh and solve jobs on their<br />own threads with progress and<br />cancel, solving physics in<br />the order of the one-way<br />couplings, mesh convergence<br />on a scratch engine; the<br />start-up request that Windows<br />not throttle the process.</div>"]
        style 26 fill:#85bbf0,stroke:#5d82a8,color:#000000
        27["<div style='font-weight: bold'>Meshing seam</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>MeshRequest to MeshedModel<br />(mesh, names, label<br />resolver); the only caller of<br />the femcore mesher; absorbing<br />layers across the two halves<br />of a project.</div>"]
        style 27 fill:#85bbf0,stroke:#5d82a8,color:#000000
        28["<div style='font-weight: bold'>Materials</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The material library plus the<br />materials a project carries,<br />resolved at the study<br />wavelength.</div>"]
        style 28 fill:#85bbf0,stroke:#5d82a8,color:#000000
        29["<div style='font-weight: bold'>Analyses</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Runs the analyses of<br />nyx-analysis on the latest<br />result of a physics.</div>"]
        style 29 fill:#85bbf0,stroke:#5d82a8,color:#000000
        30["<div style='font-weight: bold'>Results files</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>HDF5 results (serial or<br />MPI-IO), the XDMF<br />description, and opening a<br />results file as a project<br />with its results.</div>"]
        style 30 fill:#85bbf0,stroke:#5d82a8,color:#000000
        31["<div style='font-weight: bold'>Store link</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Binds the project to a<br />workflow of the workspace,<br />commits versions, records<br />runs (card and metrics made<br />here).</div>"]
        style 31 fill:#85bbf0,stroke:#5d82a8,color:#000000
        32["<div style='font-weight: bold'>Registry and facade</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The physics compiled in<br />(cargo features dfs, sms,<br />ems) and the facade that<br />re-exports what a front-end<br />draws with.</div>"]
        style 32 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      101["<div style='font-weight: bold'>geom</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Lean CSG-over-B-rep geometry<br />kernel: feature lists, named<br />regions and boundaries,<br />selectors, exact projection<br />oracles. Depends on serde and<br />optional rayon only.</div>"]
      style 101 fill:#438dd5,stroke:#2e6295,color:#ffffff
      109[("<div style='font-weight: bold'>Project file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .nyx.json]</div><div style='font-size: 80%; margin-top:10px'>Everything saved for one<br />model: geometry, materials,<br />boundary conditions, studies.<br />JSON with a format key.</div>")]
      style 109 fill:#438dd5,stroke:#2e6295,color:#ffffff
      110[("<div style='font-weight: bold'>Results file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .h5 and .xdmf]</div><div style='font-size: 80%; margin-top:10px'>Project, mesh and every<br />solution of every solved<br />physics, compressed, exact;<br />with an XDMF description next<br />to it.</div>")]
      style 110 fill:#438dd5,stroke:#2e6295,color:#ffffff
      112[("<div style='font-weight: bold'>Material library</div><div style='font-size: 70%; margin-top: 0px'>[Container: Files]</div><div style='font-size: 80%; margin-top:10px'>materials/*.toml next to the<br />program: silicon, silica,<br />sapphire, metals, nitrides,<br />niobate and more, with<br />sources.</div>")]
      style 112 fill:#438dd5,stroke:#2e6295,color:#ffffff
      33["<div style='font-weight: bold'>nyx-physics</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Physics trait, the<br />problem envelope<br />(ProblemSpec), descriptors,<br />SolveResult, the registry<br />type, validation, JSON<br />Schemas.</div>"]
      style 33 fill:#438dd5,stroke:#2e6295,color:#ffffff
      39["<div style='font-weight: bold'>nyx-dfs</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Diffusion physics: heat and<br />electrostatics. Typed<br />problem, three-stage kernels,<br />descriptors, goldens.</div>"]
      style 39 fill:#2e7d32,stroke:#205723,color:#ffffff
      45["<div style='font-weight: bold'>nyx-sms</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Solid mechanics: static,<br />eigenfrequency, PML, Floquet,<br />low-reflecting boundaries,<br />axisymmetric.</div>"]
      style 45 fill:#2e7d32,stroke:#205723,color:#ffffff
      54["<div style='font-weight: bold'>nyx-ems</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>EM waves: waveguide modes<br />(real and complex),<br />axisymmetric, cavity in 3D.</div>"]
      style 54 fill:#2e7d32,stroke:#205723,color:#ffffff
      62["<div style='font-weight: bold'>nyx-analysis</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Analyses of results: the<br />polarization (Stokes<br />parameters, Poincare sphere,<br />classes, degenerate groups)<br />as a function of a<br />SolveResult and its mesh.</div>"]
      style 62 fill:#438dd5,stroke:#2e6295,color:#ffffff
      67["<div style='font-weight: bold'>nyx-store</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Workspace, workflow,<br />experiment and run storage:<br />manifests, append-only run<br />log, versions, search, .nyx<br />archives. Knows nothing of a<br />project.</div>"]
      style 67 fill:#438dd5,stroke:#2e6295,color:#ffffff
      76["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 76 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    39-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->33
    45-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->33
    54-. "<div>Implements Physics</div><div style='font-size: 70%'>[Rust API]</div>" .->33
    39-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->76
    45-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->76
    54-. "<div>Mesh, assembly, linear<br />algebra</div><div style='font-size: 70%'>[Rust API]</div>" .->76
    62-. "<div>Reads SolveResult</div><div style='font-size: 70%'>[Rust API]</div>" .->33
    62-. "<div>Mesh geometry, quadrature</div><div style='font-size: 70%'>[Rust API]</div>" .->76
    33-. "<div>Mesh, materials</div><div style='font-size: 70%'>[Rust API]</div>" .->76
    76-. "<div>Evaluates models, reads<br />selectors and labels</div><div style='font-size: 70%'>[Rust API]</div>" .->101
    29-. "<div>uses analysis: crate root,<br />modes +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->62
    29-. "<div>uses physics: error, result</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    29-. "<div>uses engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    29-. "<div>uses meshing</div><div style='font-size: 70%'>[Rust paths]</div>" .->27
    23-. "<div>uses femcore: material::spec,<br />progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    23-. "<div>uses geom: model</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    23-. "<div>uses physics: contract, crate<br />root +5</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    23-. "<div>uses store: util</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    23-. "<div>uses geomdesc, geometry</div><div style='font-size: 70%'>[Rust paths]</div>" .->25
    23-. "<div>uses jobs, solve</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    23-. "<div>uses materials</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    23-. "<div>uses meshing</div><div style='font-size: 70%'>[Rust paths]</div>" .->27
    23-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    23-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->32
    23-. "<div>uses store_link</div><div style='font-size: 70%'>[Rust paths]</div>" .->31
    25-. "<div>uses geom: error, eval +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    25-. "<div>uses physics: error</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    26-. "<div>uses femcore: material::spec,<br />progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    26-. "<div>uses physics: error, physics<br />+3</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    26-. "<div>uses engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    26-. "<div>uses materials</div><div style='font-size: 70%'>[Rust paths]</div>" .->28
    26-. "<div>uses meshing</div><div style='font-size: 70%'>[Rust paths]</div>" .->27
    26-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    26-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->32
    28-. "<div>uses femcore: material,<br />material::spec</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    28-. "<div>uses physics: desc, error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    28-. "<div>uses engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    27-. "<div>uses femcore: mesh,<br />meshing::geomdrive +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    27-. "<div>uses geom: eval, labels +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    27-. "<div>uses physics: crate root,<br />error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    27-. "<div>uses engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    27-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    27-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->32
    24-. "<div>uses geom: model</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    24-. "<div>uses physics: error, registry<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    24-. "<div>uses crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->32
    32-. "<div>uses analysis: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->62
    32-. "<div>uses dfs: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->39
    32-. "<div>uses ems: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->54
    32-. "<div>uses femcore: fields, mesh +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    32-. "<div>uses geom: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    32-. "<div>uses physics: crate root,<br />error +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    32-. "<div>uses sms: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->45
    32-. "<div>uses store: crate root</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    32-. "<div>uses commands, engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    32-. "<div>uses jobs</div><div style='font-size: 70%'>[Rust paths]</div>" .->26
    32-. "<div>uses meshing</div><div style='font-size: 70%'>[Rust paths]</div>" .->27
    32-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    30-. "<div>uses femcore: mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    30-. "<div>uses physics: error, result</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    30-. "<div>uses engine</div><div style='font-size: 70%'>[Rust paths]</div>" .->23
    30-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    31-. "<div>uses femcore: mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->76
    31-. "<div>uses geom: expr, model</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    31-. "<div>uses physics: error, result<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->33
    31-. "<div>uses store: flow, model +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->67
    31-. "<div>uses geomdesc</div><div style='font-size: 70%'>[Rust paths]</div>" .->25
    31-. "<div>uses project</div><div style='font-size: 70%'>[Rust paths]</div>" .->24
    31-. "<div>uses h5io, xdmf</div><div style='font-size: 70%'>[Rust paths]</div>" .->30
    28-. "<div>Loads</div><div style='font-size: 70%'>[declared in spec]</div>" .->112
    24-. "<div>Saves and loads</div><div style='font-size: 70%'>[declared in spec]</div>" .->109
    30-. "<div>Writes and opens</div><div style='font-size: 70%'>[declared in spec]</div>" .->110

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Commands](../4-code/engine/commands.md) | 3 | 987 | 57, 38, 15 | The Engine and its catalog: every action by name with JSON arguments, JSON Schema and JSON reply; the dispatcher; the guide cards a host reads on demand. The window and the MCP server use the same entries. |
| [Project](../4-code/engine/project.md) | 1 | 197 | 63 | The saved state: feature list, mesh settings, one ProblemSpec per physics, project materials. The format key is checked at the file boundary. |
| [Geometry operations](../4-code/engine/geometry.md) | 2 | 730 | 38, 3 | Every edit of a geom feature list as a function that does it or says why not, and the descriptors of the geom operations (they live here because geom depends on serde only). |
| [Job runner](../4-code/engine/jobs.md) | 4 | 442 | 61, 11 | Mesh and solve jobs on their own threads with progress and cancel, solving physics in the order of the one-way couplings, mesh convergence on a scratch engine; the start-up request that Windows not throttle the process. |
| [Meshing seam](../4-code/engine/meshing.md) | 2 | 334 | 54, 50 | MeshRequest to MeshedModel (mesh, names, label resolver); the only caller of the femcore mesher; absorbing layers across the two halves of a project. |
| [Materials](../4-code/engine/materials.md) | 1 | 155 | 58 | The material library plus the materials a project carries, resolved at the study wavelength. |
| [Analyses](../4-code/engine/analyse.md) | 1 | 80 | 56 | Runs the analyses of nyx-analysis on the latest result of a physics. |
| [Results files](../4-code/engine/resultsIO.md) | 3 | 940 | 62, 63 | HDF5 results (serial or MPI-IO), the XDMF description, and opening a results file as a project with its results. |
| [Store link](../4-code/engine/storeLink.md) | 1 | 272 | 39, 59 | Binds the project to a workflow of the workspace, commits versions, records runs (card and metrics made here). |
| [Registry and facade](../4-code/engine/registry.md) | 1 | 85 | 25 | The physics compiled in (cargo features dfs, sms, ems) and the facade that re-exports what a front-end draws with. |
