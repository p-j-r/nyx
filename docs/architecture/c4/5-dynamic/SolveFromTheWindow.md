# Dynamic view: SolveFromTheWindow

What happens when the window solves the active physics. Export: [Mermaid](../export/mermaid/structurizr-SolveFromTheWindow.mmd) · [PlantUML](../export/plantuml/structurizr-SolveFromTheWindow.puml)

```mermaid
graph LR
  linkStyle default fill:#ffffff

  subgraph diagram ["Dynamic View: nyx"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      117[("<div style='font-weight: bold'>Results file</div><div style='font-size: 70%; margin-top: 0px'>[Container: File .h5 and .xdmf]</div><div style='font-size: 80%; margin-top:10px'>Project, mesh and every<br />solution of every solved<br />physics, compressed, exact;<br />with an XDMF description next<br />to it.</div>")]
      style 117 fill:#438dd5,stroke:#2e6295,color:#ffffff
      27["<div style='font-weight: bold'>nyx-engine</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The Engine: one project, the<br />mesh built for it, the latest<br />result of each physics, the<br />jobs that produce them, and<br />the binding to a workflow.<br />Every action is a named<br />command with a JSON Schema.</div>"]
      style 27 fill:#438dd5,stroke:#2e6295,color:#ffffff
      6("<div style='font-weight: bold'>app</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust, eframe/egui]</div><div style='font-size: 80%; margin-top:10px'>The egui window (binary nyx):<br />geometry editor, mesh view,<br />Study, Physics and Materials<br />tabs drawn from the physics<br />descriptors, results,<br />analysis, workflow dashboard.<br />A front-end of the engine: it<br />holds an Engine and calls<br />commands.</div>")
      style 6 fill:#438dd5,stroke:#2e6295,color:#ffffff
      61["<div style='font-weight: bold'>nyx-ems</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>EM waves: waveguide modes<br />(real and complex),<br />axisymmetric, cavity in 3D.</div>"]
      style 61 fill:#2e7d32,stroke:#205723,color:#ffffff
      74["<div style='font-weight: bold'>nyx-store</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Workspace, workflow,<br />experiment and run storage:<br />manifests, append-only run<br />log, versions, search, .nyx<br />archives. Knows nothing of a<br />project.</div>"]
      style 74 fill:#438dd5,stroke:#2e6295,color:#ffffff
      83["<div style='font-weight: bold'>femcore</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>The FEM core: meshes, meshing<br />and adaptation, linear<br />algebra behind one driver,<br />finite-element basics, the<br />material specification and<br />library, expressions.</div>"]
      style 83 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    6-. "<div>1. start_solve(targets)</div><div style='font-size: 70%'>[Rust API]</div>" .->27
    27-. "<div>2. mesh the model when none<br />is current</div><div style='font-size: 70%'>[Rust API]</div>" .->83
    27-. "<div>3. Physics::solve(mesh, spec,<br />ctx)</div><div style='font-size: 70%'>[Rust API]</div>" .->61
    61-. "<div>4. assemble, shift-invert<br />eigensolve</div><div style='font-size: 70%'>[Rust API]</div>" .->83
    27-. "<div>5. record the run, freeze the<br />spec</div><div style='font-size: 70%'>[Rust API]</div>" .->74
    27-. "<div>6. write results.h5 and<br />results.xdmf</div><div style='font-size: 70%'>[HDF5 (serial or MPI-IO), XDMF]</div>" .->117
    6-. "<div>7. poll: take the mesh and<br />the results</div><div style='font-size: 70%'>[Rust API]</div>" .->27

  end```
