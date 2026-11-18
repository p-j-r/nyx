# Level 3: components of femcore

Crate `crates/femcore`. The FEM core: meshes, meshing and adaptation, linear algebra behind one driver, finite-element basics, the material specification and library, expressions.

A component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. Export: [Mermaid](../export/mermaid/structurizr-Components_femcore.mmd) · [PlantUML](../export/plantuml/structurizr-Components_femcore.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - femcore"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 82 ["femcore"]
        style 82 fill:#ffffff,stroke:#2e6295,color:#2e6295

        100["<div style='font-weight: bold'>Time, nonlinear, sweeps</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Time integration, Newton and<br />optimisation, matrix<br />exponential action,<br />reduced-order frequency<br />sweeps and matrix-free<br />high-order operators.</div>"]
        style 100 fill:#85bbf0,stroke:#5d82a8,color:#000000
        101["<div style='font-weight: bold'>Finite-element basics</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Lagrange P1/P2 elements and<br />quadrature, the 2D<br />edge-element basis, derived<br />fields, PML coordinate<br />stretching, Floquet pairs and<br />material properties per<br />domain.</div>"]
        style 101 fill:#85bbf0,stroke:#5d82a8,color:#000000
        102["<div style='font-weight: bold'>Material specification</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>TOML material files with<br />explicit units, tensors,<br />dispersion models evaluated<br />at the wavelength, provenance<br />and the library loader.</div>"]
        style 102 fill:#85bbf0,stroke:#5d82a8,color:#000000
        103["<div style='font-weight: bold'>Optical mode sets</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Optical modes of a 2D<br />cross-section as a<br />field-evaluable set with the<br />material data to evaluate<br />them.</div>"]
        style 103 fill:#85bbf0,stroke:#5d82a8,color:#000000
        104["<div style='font-weight: bold'>Frame and sections</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The one place that knows the<br />wafer axes (x lateral, z<br />normal, waves along +y).</div>"]
        style 104 fill:#85bbf0,stroke:#5d82a8,color:#000000
        105["<div style='font-weight: bold'>Expressions</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>A small safe expression<br />evaluator for user input (a<br />copy of geom's kept equal by<br />a corpus test), and the typed<br />tree over complex field<br />values that the reductions of<br />results evaluate, and the<br />text formulas parsed into it<br />(ADR 90).</div>"]
        style 105 fill:#85bbf0,stroke:#5d82a8,color:#000000
        106["<div style='font-weight: bold'>Runtime services</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The shared worker pool, stage<br />log and cancel, and refusal<br />of NaN and infinity.</div>"]
        style 106 fill:#85bbf0,stroke:#5d82a8,color:#000000
        83["<div style='font-weight: bold'>Mesh</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The native mixed-cell mesh,<br />the exchange mesh the physics<br />take, simplex meshes with<br />their numbering of unknowns,<br />draw-only meshes and quality<br />statistics.</div>"]
        style 83 fill:#85bbf0,stroke:#5d82a8,color:#000000
        84["<div style='font-weight: bold'>Mesh files</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The compact lossless .fmsh<br />format for meshes and results<br />and its entropy coder.</div>"]
        style 84 fill:#85bbf0,stroke:#5d82a8,color:#000000
        85["<div style='font-weight: bold'>Fixtures and test kit</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Test meshes and the SOI wire<br />cross-section used by the<br />tests of every physics crate.</div>"]
        style 85 fill:#85bbf0,stroke:#5d82a8,color:#000000
        86["<div style='font-weight: bold'>Meshing facade</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>MeshPlan and mesh_model (the<br />one way to mesh a geom<br />model), the serialisable<br />request and report, sizing<br />from the physics, and the 2D<br />triangulation core.</div>"]
        style 86 fill:#85bbf0,stroke:#5d82a8,color:#000000
        87["<div style='font-weight: bold'>3D meshing</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The lattice tetrahedral<br />mesher on Domain3, Rivara<br />refinement, the adapter from<br />a 3D geom model, and adaptive<br />meshing driven by exact<br />faces.</div>"]
        style 87 fill:#85bbf0,stroke:#5d82a8,color:#000000
        88["<div style='font-weight: bold'>Structured absorbing layers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Structured PML around a<br />box-shaped or geom-model<br />mesh: prisms, hexahedra,<br />quads, labels.</div>"]
        style 88 fill:#85bbf0,stroke:#5d82a8,color:#000000
        89["<div style='font-weight: bold'>Meshing validation harness</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>A small independent scalar<br />Helmholtz solver with PML<br />used to validate meshes.</div>"]
        style 89 fill:#85bbf0,stroke:#5d82a8,color:#000000
        90["<div style='font-weight: bold'>Metric fields</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Riemannian metrics in 2D and<br />3D, metric fields, the<br />simplex working mesh and<br />per-phase timers.</div>"]
        style 90 fill:#85bbf0,stroke:#5d82a8,color:#000000
        91["<div style='font-weight: bold'>Adaptation operators</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Edge split, collapse, flip<br />and smoothing on the working<br />meshes in 2D and 3D, with the<br />vertex classification.</div>"]
        style 91 fill:#85bbf0,stroke:#5d82a8,color:#000000
        92["<div style='font-weight: bold'>Adaptation drivers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Passes of operators until the<br />mesh conforms to the metric,<br />with geometry oracles saying<br />where a boundary vertex may<br />go.</div>"]
        style 92 fill:#85bbf0,stroke:#5d82a8,color:#000000
        93["<div style='font-weight: bold'>Error estimators</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Error estimation and the<br />solve-estimate-adapt loop;<br />the Maxwell residual<br />estimator for edge elements.</div>"]
        style 93 fill:#85bbf0,stroke:#5d82a8,color:#000000
        94["<div style='font-weight: bold'>Matrices and vectors</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The one matrix type Csr over<br />a generic scalar, block<br />storage, multivector kernels,<br />operator interfaces, dense<br />helpers, Matrix Market,<br />problem bundles and<br />test-problem generators.</div>"]
        style 94 fill:#85bbf0,stroke:#5d82a8,color:#000000
        95["<div style='font-weight: bold'>Linear algebra driver</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The door the physics kernels<br />use: a direct solve,<br />eigenpairs near a shift, and<br />the solver interface for<br />models.</div>"]
        style 95 fill:#85bbf0,stroke:#5d82a8,color:#000000
        96["<div style='font-weight: bold'>Direct solvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Shifted sparse direct solves<br />with symbolic reuse, LU with<br />a user ordering, block<br />low-rank factorisation, mixed<br />precision, nested dissection.</div>"]
        style 96 fill:#85bbf0,stroke:#5d82a8,color:#000000
        97["<div style='font-weight: bold'>Krylov solvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>CG, MINRES, GMRES, BiCGStab,<br />COCR, QMR, IDR(s), GCRO-DR<br />for real and complex systems.</div>"]
        style 97 fill:#85bbf0,stroke:#5d82a8,color:#000000
        98["<div style='font-weight: bold'>Eigensolvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Krylov-Schur, LOBPCG,<br />Jacobi-Davidson, quadratic<br />and nonlinear eigenproblems,<br />contour integral, Chebyshev<br />filtering, real Arnoldi for<br />the 2D modes, sensitivities.</div>"]
        style 98 fill:#85bbf0,stroke:#5d82a8,color:#000000
        99["<div style='font-weight: bold'>Preconditioners</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Smoothed-aggregation<br />multigrid (scalar and<br />elastic), auxiliary-space<br />Maxwell and divergence<br />solvers, Schwarz, field-split<br />and incomplete<br />factorisations.</div>"]
        style 99 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

      107["<div style='font-weight: bold'>geom</div><div style='font-size: 70%; margin-top: 0px'>[Container: Rust library]</div><div style='font-size: 80%; margin-top:10px'>Lean CSG-over-B-rep geometry<br />kernel: feature lists, named<br />regions and boundaries,<br />selectors, exact projection<br />oracles. Depends on serde and<br />optional rayon only.</div>"]
      style 107 fill:#438dd5,stroke:#2e6295,color:#ffffff
    end

    92-. "<div>uses adapt::ops, adapt::ops3<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    92-. "<div>uses adapt::field,<br />adapt::metric +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    92-. "<div>uses progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    91-. "<div>uses adapt::drive,<br />adapt::geom +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    91-. "<div>uses mesh::quality</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    91-. "<div>uses adapt::field,<br />adapt::metric +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    93-. "<div>uses geom: model</div><div style='font-size: 70%'>[Rust paths]</div>" .->107
    93-. "<div>uses adapt::drive3</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    93-. "<div>uses adapt::ops3,<br />adapt::topo3</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    93-. "<div>uses fem</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    93-. "<div>uses linalg::csr::assembly</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    93-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    93-. "<div>uses meshing,<br />meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    93-. "<div>uses adapt, adapt::field +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    105-. "<div>uses linalg::csr::assembly</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    105-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    101-. "<div>uses linalg::csr,<br />linalg::csr::assembly</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    101-. "<div>uses linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->95
    101-. "<div>uses linalg::amg_elastic</div><div style='font-size: 70%'>[Rust paths]</div>" .->99
    101-. "<div>uses material</div><div style='font-size: 70%'>[Rust paths]</div>" .->102
    101-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    101-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    101-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    85-. "<div>uses geom: model, select</div><div style='font-size: 70%'>[Rust paths]</div>" .->107
    85-. "<div>uses mesh, mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    85-. "<div>uses meshing,<br />meshing::geomdrive +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    85-. "<div>uses adapt::metric,<br />adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    94-. "<div>uses par</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    96-. "<div>uses linalg::csr,<br />linalg::faer_glue +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    95-. "<div>uses linalg::basis,<br />linalg::csr +5</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    95-. "<div>uses linalg::direct,<br />linalg::ordering</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    95-. "<div>uses linalg::arnoldi_real,<br />linalg::eps +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    95-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    98-. "<div>uses linalg::basis,<br />linalg::csr +7</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    98-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    98-. "<div>uses linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->95
    98-. "<div>uses linalg::ksp,<br />linalg::ksp::gmres</div><div style='font-size: 70%'>[Rust paths]</div>" .->97
    98-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    97-. "<div>uses linalg::basis,<br />linalg::dense +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    100-. "<div>uses linalg::csr,<br />linalg::dense +5</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    100-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    100-. "<div>uses linalg::nep</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    100-. "<div>uses linalg::ksp,<br />linalg::ksp::gmres</div><div style='font-size: 70%'>[Rust paths]</div>" .->97
    99-. "<div>uses linalg::csr,<br />linalg::csr::assembly +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    99-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    99-. "<div>uses linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->95
    99-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    83-. "<div>uses adapt::metric,<br />adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    83-. "<div>uses par</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    89-. "<div>uses mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    89-. "<div>uses adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    89-. "<div>uses meshing::pml_box</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    86-. "<div>uses geom: diff, eval +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->107
    86-. "<div>uses adapt::drive,<br />adapt::drive3 +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    86-. "<div>uses adapt::ops3, adapt::topo<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    86-. "<div>uses props</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    86-. "<div>uses material</div><div style='font-size: 70%'>[Rust paths]</div>" .->102
    86-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    86-. "<div>uses meshing::geomdrive3</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    86-. "<div>uses adapt, adapt::field +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    86-. "<div>uses meshing::pml_geom</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    87-. "<div>uses geom: eval, labels +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->107
    87-. "<div>uses adapt::drive,<br />adapt::drive3 +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    87-. "<div>uses adapt::ops3, adapt::topo<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    87-. "<div>uses mesh, mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    87-. "<div>uses meshing,<br />meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    87-. "<div>uses adapt::field3,<br />adapt::metric +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    87-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    84-. "<div>uses mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    90-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    90-. "<div>uses par</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    103-. "<div>uses edge, fem +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    103-. "<div>uses linalg::csr,<br />linalg::csr::assembly</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    103-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    103-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    104-. "<div>uses props</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    104-. "<div>uses material</div><div style='font-size: 70%'>[Rust paths]</div>" .->102
    104-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    88-. "<div>uses geom: eval, model</div><div style='font-size: 70%'>[Rust paths]</div>" .->107
    88-. "<div>uses props</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    88-. "<div>uses mesh, mesh::cells +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    88-. "<div>uses meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    88-. "<div>uses meshing::adapt3d</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    88-. "<div>uses adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90

  end```

| Component | Files | Lines | Decisions | What it does |
|---|---:|---:|---|---|
| [Mesh](../4-code/femcore/mesh.md) | 5 | 1594 | 4, 5, 32 | The native mixed-cell mesh, the exchange mesh the physics take, simplex meshes with their numbering of unknowns, draw-only meshes and quality statistics. |
| [Mesh files](../4-code/femcore/meshio.md) | 3 | 878 | 16 | The compact lossless .fmsh format for meshes and results and its entropy coder. |
| [Fixtures and test kit](../4-code/femcore/fixtures.md) | 3 | 582 | 41, 10, 60 | Test meshes and the SOI wire cross-section used by the tests of every physics crate. |
| [Meshing facade](../4-code/femcore/meshing.md) | 6 | 1707 | 6, 17, 54 | MeshPlan and mesh_model (the one way to mesh a geom model), the serialisable request and report, sizing from the physics, and the 2D triangulation core. |
| [3D meshing](../4-code/femcore/meshing3d.md) | 4 | 1677 | 40, 50, 51, 55 | The lattice tetrahedral mesher on Domain3, Rivara refinement, the adapter from a 3D geom model, and adaptive meshing driven by exact faces. |
| [Structured absorbing layers](../4-code/femcore/structpml.md) | 4 | 1408 | 7, 52, 5 | Structured PML around a box-shaped or geom-model mesh: prisms, hexahedra, quads, labels. |
| [Meshing validation harness](../4-code/femcore/meshcheck.md) | 4 | 980 | 7 | A small independent scalar Helmholtz solver with PML used to validate meshes. |
| [Metric fields](../4-code/femcore/metric.md) | 6 | 1846 | 7, 5 | Riemannian metrics in 2D and 3D, metric fields, the simplex working mesh and per-phase timers. |
| [Adaptation operators](../4-code/femcore/adaptops.md) | 4 | 2255 | 7, 51 | Edge split, collapse, flip and smoothing on the working meshes in 2D and 3D, with the vertex classification. |
| [Adaptation drivers](../4-code/femcore/adaptdrive.md) | 4 | 720 | 7, 51 | Passes of operators until the mesh conforms to the metric, with geometry oracles saying where a boundary vertex may go. |
| [Error estimators](../4-code/femcore/estimators.md) | 3 | 1125 | 7, 17, 53 | Error estimation and the solve-estimate-adapt loop; the Maxwell residual estimator for edge elements. |
| [Matrices and vectors](../4-code/femcore/lacore.md) | 17 | 4131 | 33, 8, 11 | The one matrix type Csr over a generic scalar, block storage, multivector kernels, operator interfaces, dense helpers, Matrix Market, problem bundles and test-problem generators. |
| [Linear algebra driver](../4-code/femcore/ladriver.md) | 3 | 731 | 33 | The door the physics kernels use: a direct solve, eigenpairs near a shift, and the solver interface for models. |
| [Direct solvers](../4-code/femcore/ladirect.md) | 6 | 3054 | 33, 8 | Shifted sparse direct solves with symbolic reuse, LU with a user ordering, block low-rank factorisation, mixed precision, nested dissection. |
| [Krylov solvers](../4-code/femcore/laksp.md) | 9 | 2131 | 33, 8 | CG, MINRES, GMRES, BiCGStab, COCR, QMR, IDR(s), GCRO-DR for real and complex systems. |
| [Eigensolvers](../4-code/femcore/laeigen.md) | 12 | 4560 | 33, 8, 11 | Krylov-Schur, LOBPCG, Jacobi-Davidson, quadratic and nonlinear eigenproblems, contour integral, Chebyshev filtering, real Arnoldi for the 2D modes, sensitivities. |
| [Preconditioners](../4-code/femcore/laprecond.md) | 8 | 2516 | 33, 8 | Smoothed-aggregation multigrid (scalar and elastic), auxiliary-space Maxwell and divergence solvers, Schwarz, field-split and incomplete factorisations. |
| [Time, nonlinear, sweeps](../4-code/femcore/laother.md) | 10 | 3627 | 33, 8 | Time integration, Newton and optimisation, matrix exponential action, reduced-order frequency sweeps and matrix-free high-order operators. |
| [Finite-element basics](../4-code/femcore/fem.md) | 6 | 1826 | 9, 4, 45 | Lagrange P1/P2 elements and quadrature, the 2D edge-element basis, derived fields, PML coordinate stretching, Floquet pairs and material properties per domain. |
| [Material specification](../4-code/femcore/material.md) | 9 | 3519 | 29, 30, 34, 35 | TOML material files with explicit units, tensors, dispersion models evaluated at the wavelength, provenance and the library loader. |
| [Optical mode sets](../4-code/femcore/optical.md) | 1 | 231 | 21, 9, 45 | Optical modes of a 2D cross-section as a field-evaluable set with the material data to evaluate them. |
| [Frame and sections](../4-code/femcore/section.md) | 1 | 221 | 1 | The one place that knows the wafer axes (x lateral, z normal, waves along +y). |
| [Expressions](../4-code/femcore/expr.md) | 3 | 1562 | 65, 2, 90 | A small safe expression evaluator for user input (a copy of geom's kept equal by a corpus test), and the typed tree over complex field values that the reductions of results evaluate, and the text formulas parsed into it (ADR 90). |
| [Runtime services](../4-code/femcore/runtime.md) | 4 | 538 | 9, 11 | The shared worker pool, stage log and cancel, and refusal of NaN and infinity. |

## Part: Meshing and adaptation

The same components, only those of this part (fewer arrows to read). Export: [Mermaid](../export/mermaid/structurizr-Components_femcore_meshing.mmd) · [PlantUML](../export/plantuml/structurizr-Components_femcore_meshing.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - femcore"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 82 ["femcore"]
        style 82 fill:#ffffff,stroke:#2e6295,color:#2e6295

        83["<div style='font-weight: bold'>Mesh</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The native mixed-cell mesh,<br />the exchange mesh the physics<br />take, simplex meshes with<br />their numbering of unknowns,<br />draw-only meshes and quality<br />statistics.</div>"]
        style 83 fill:#85bbf0,stroke:#5d82a8,color:#000000
        84["<div style='font-weight: bold'>Mesh files</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The compact lossless .fmsh<br />format for meshes and results<br />and its entropy coder.</div>"]
        style 84 fill:#85bbf0,stroke:#5d82a8,color:#000000
        85["<div style='font-weight: bold'>Fixtures and test kit</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Test meshes and the SOI wire<br />cross-section used by the<br />tests of every physics crate.</div>"]
        style 85 fill:#85bbf0,stroke:#5d82a8,color:#000000
        86["<div style='font-weight: bold'>Meshing facade</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>MeshPlan and mesh_model (the<br />one way to mesh a geom<br />model), the serialisable<br />request and report, sizing<br />from the physics, and the 2D<br />triangulation core.</div>"]
        style 86 fill:#85bbf0,stroke:#5d82a8,color:#000000
        87["<div style='font-weight: bold'>3D meshing</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The lattice tetrahedral<br />mesher on Domain3, Rivara<br />refinement, the adapter from<br />a 3D geom model, and adaptive<br />meshing driven by exact<br />faces.</div>"]
        style 87 fill:#85bbf0,stroke:#5d82a8,color:#000000
        88["<div style='font-weight: bold'>Structured absorbing layers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Structured PML around a<br />box-shaped or geom-model<br />mesh: prisms, hexahedra,<br />quads, labels.</div>"]
        style 88 fill:#85bbf0,stroke:#5d82a8,color:#000000
        89["<div style='font-weight: bold'>Meshing validation harness</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>A small independent scalar<br />Helmholtz solver with PML<br />used to validate meshes.</div>"]
        style 89 fill:#85bbf0,stroke:#5d82a8,color:#000000
        90["<div style='font-weight: bold'>Metric fields</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Riemannian metrics in 2D and<br />3D, metric fields, the<br />simplex working mesh and<br />per-phase timers.</div>"]
        style 90 fill:#85bbf0,stroke:#5d82a8,color:#000000
        91["<div style='font-weight: bold'>Adaptation operators</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Edge split, collapse, flip<br />and smoothing on the working<br />meshes in 2D and 3D, with the<br />vertex classification.</div>"]
        style 91 fill:#85bbf0,stroke:#5d82a8,color:#000000
        92["<div style='font-weight: bold'>Adaptation drivers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Passes of operators until the<br />mesh conforms to the metric,<br />with geometry oracles saying<br />where a boundary vertex may<br />go.</div>"]
        style 92 fill:#85bbf0,stroke:#5d82a8,color:#000000
        93["<div style='font-weight: bold'>Error estimators</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Error estimation and the<br />solve-estimate-adapt loop;<br />the Maxwell residual<br />estimator for edge elements.</div>"]
        style 93 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    92-. "<div>uses adapt::ops, adapt::ops3<br />+2</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    92-. "<div>uses adapt::field,<br />adapt::metric +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    91-. "<div>uses adapt::drive,<br />adapt::geom +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    91-. "<div>uses mesh::quality</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    91-. "<div>uses adapt::field,<br />adapt::metric +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    93-. "<div>uses adapt::drive3</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    93-. "<div>uses adapt::ops3,<br />adapt::topo3</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    93-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    93-. "<div>uses meshing,<br />meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    93-. "<div>uses adapt, adapt::field +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    85-. "<div>uses mesh, mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    85-. "<div>uses meshing,<br />meshing::geomdrive +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    85-. "<div>uses adapt::metric,<br />adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    83-. "<div>uses adapt::metric,<br />adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    89-. "<div>uses mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    89-. "<div>uses adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    89-. "<div>uses meshing::pml_box</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    86-. "<div>uses adapt::drive,<br />adapt::drive3 +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    86-. "<div>uses adapt::ops3, adapt::topo<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    86-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    86-. "<div>uses meshing::geomdrive3</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    86-. "<div>uses adapt, adapt::field +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    86-. "<div>uses meshing::pml_geom</div><div style='font-size: 70%'>[Rust paths]</div>" .->88
    87-. "<div>uses adapt::drive,<br />adapt::drive3 +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->92
    87-. "<div>uses adapt::ops3, adapt::topo<br />+1</div><div style='font-size: 70%'>[Rust paths]</div>" .->91
    87-. "<div>uses mesh, mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    87-. "<div>uses meshing,<br />meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    87-. "<div>uses adapt::field3,<br />adapt::metric +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->90
    84-. "<div>uses mesh::cells</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    90-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    88-. "<div>uses mesh, mesh::cells +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    88-. "<div>uses meshing::geomdrive</div><div style='font-size: 70%'>[Rust paths]</div>" .->86
    88-. "<div>uses meshing::adapt3d</div><div style='font-size: 70%'>[Rust paths]</div>" .->87
    88-. "<div>uses adapt::work</div><div style='font-size: 70%'>[Rust paths]</div>" .->90

  end```

## Part: Linear algebra

The same components, only those of this part (fewer arrows to read). Export: [Mermaid](../export/mermaid/structurizr-Components_femcore_linalg.mmd) · [PlantUML](../export/plantuml/structurizr-Components_femcore_linalg.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - femcore"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 82 ["femcore"]
        style 82 fill:#ffffff,stroke:#2e6295,color:#2e6295

        100["<div style='font-weight: bold'>Time, nonlinear, sweeps</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Time integration, Newton and<br />optimisation, matrix<br />exponential action,<br />reduced-order frequency<br />sweeps and matrix-free<br />high-order operators.</div>"]
        style 100 fill:#85bbf0,stroke:#5d82a8,color:#000000
        94["<div style='font-weight: bold'>Matrices and vectors</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The one matrix type Csr over<br />a generic scalar, block<br />storage, multivector kernels,<br />operator interfaces, dense<br />helpers, Matrix Market,<br />problem bundles and<br />test-problem generators.</div>"]
        style 94 fill:#85bbf0,stroke:#5d82a8,color:#000000
        95["<div style='font-weight: bold'>Linear algebra driver</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The door the physics kernels<br />use: a direct solve,<br />eigenpairs near a shift, and<br />the solver interface for<br />models.</div>"]
        style 95 fill:#85bbf0,stroke:#5d82a8,color:#000000
        96["<div style='font-weight: bold'>Direct solvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Shifted sparse direct solves<br />with symbolic reuse, LU with<br />a user ordering, block<br />low-rank factorisation, mixed<br />precision, nested dissection.</div>"]
        style 96 fill:#85bbf0,stroke:#5d82a8,color:#000000
        97["<div style='font-weight: bold'>Krylov solvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>CG, MINRES, GMRES, BiCGStab,<br />COCR, QMR, IDR(s), GCRO-DR<br />for real and complex systems.</div>"]
        style 97 fill:#85bbf0,stroke:#5d82a8,color:#000000
        98["<div style='font-weight: bold'>Eigensolvers</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Krylov-Schur, LOBPCG,<br />Jacobi-Davidson, quadratic<br />and nonlinear eigenproblems,<br />contour integral, Chebyshev<br />filtering, real Arnoldi for<br />the 2D modes, sensitivities.</div>"]
        style 98 fill:#85bbf0,stroke:#5d82a8,color:#000000
        99["<div style='font-weight: bold'>Preconditioners</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Smoothed-aggregation<br />multigrid (scalar and<br />elastic), auxiliary-space<br />Maxwell and divergence<br />solvers, Schwarz, field-split<br />and incomplete<br />factorisations.</div>"]
        style 99 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    96-. "<div>uses linalg::csr,<br />linalg::faer_glue +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    95-. "<div>uses linalg::basis,<br />linalg::csr +5</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    95-. "<div>uses linalg::direct,<br />linalg::ordering</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    95-. "<div>uses linalg::arnoldi_real,<br />linalg::eps +1</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    98-. "<div>uses linalg::basis,<br />linalg::csr +7</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    98-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    98-. "<div>uses linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->95
    98-. "<div>uses linalg::ksp,<br />linalg::ksp::gmres</div><div style='font-size: 70%'>[Rust paths]</div>" .->97
    97-. "<div>uses linalg::basis,<br />linalg::dense +4</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    100-. "<div>uses linalg::csr,<br />linalg::dense +5</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    100-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    100-. "<div>uses linalg::nep</div><div style='font-size: 70%'>[Rust paths]</div>" .->98
    100-. "<div>uses linalg::ksp,<br />linalg::ksp::gmres</div><div style='font-size: 70%'>[Rust paths]</div>" .->97
    99-. "<div>uses linalg::csr,<br />linalg::csr::assembly +3</div><div style='font-size: 70%'>[Rust paths]</div>" .->94
    99-. "<div>uses linalg::direct</div><div style='font-size: 70%'>[Rust paths]</div>" .->96
    99-. "<div>uses linalg::driver</div><div style='font-size: 70%'>[Rust paths]</div>" .->95

  end```

## Part: Finite elements, materials and services

The same components, only those of this part (fewer arrows to read). Export: [Mermaid](../export/mermaid/structurizr-Components_femcore_fem.mmd) · [PlantUML](../export/plantuml/structurizr-Components_femcore_fem.puml)

```mermaid
graph TB
  linkStyle default fill:#ffffff

  subgraph diagram ["Component View: nyx - femcore"]
    style diagram fill:#ffffff,stroke:#ffffff

    subgraph 5 ["nyx"]
      style 5 fill:#ffffff,stroke:#0b4884,color:#0b4884

      subgraph 82 ["femcore"]
        style 82 fill:#ffffff,stroke:#2e6295,color:#2e6295

        101["<div style='font-weight: bold'>Finite-element basics</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Lagrange P1/P2 elements and<br />quadrature, the 2D<br />edge-element basis, derived<br />fields, PML coordinate<br />stretching, Floquet pairs and<br />material properties per<br />domain.</div>"]
        style 101 fill:#85bbf0,stroke:#5d82a8,color:#000000
        102["<div style='font-weight: bold'>Material specification</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>TOML material files with<br />explicit units, tensors,<br />dispersion models evaluated<br />at the wavelength, provenance<br />and the library loader.</div>"]
        style 102 fill:#85bbf0,stroke:#5d82a8,color:#000000
        103["<div style='font-weight: bold'>Optical mode sets</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>Optical modes of a 2D<br />cross-section as a<br />field-evaluable set with the<br />material data to evaluate<br />them.</div>"]
        style 103 fill:#85bbf0,stroke:#5d82a8,color:#000000
        104["<div style='font-weight: bold'>Frame and sections</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The one place that knows the<br />wafer axes (x lateral, z<br />normal, waves along +y).</div>"]
        style 104 fill:#85bbf0,stroke:#5d82a8,color:#000000
        105["<div style='font-weight: bold'>Expressions</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>A small safe expression<br />evaluator for user input (a<br />copy of geom's kept equal by<br />a corpus test), and the typed<br />tree over complex field<br />values that the reductions of<br />results evaluate, and the<br />text formulas parsed into it<br />(ADR 90).</div>"]
        style 105 fill:#85bbf0,stroke:#5d82a8,color:#000000
        106["<div style='font-weight: bold'>Runtime services</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The shared worker pool, stage<br />log and cancel, and refusal<br />of NaN and infinity.</div>"]
        style 106 fill:#85bbf0,stroke:#5d82a8,color:#000000
        83["<div style='font-weight: bold'>Mesh</div><div style='font-size: 70%; margin-top: 0px'>[Component: Rust modules]</div><div style='font-size: 80%; margin-top:10px'>The native mixed-cell mesh,<br />the exchange mesh the physics<br />take, simplex meshes with<br />their numbering of unknowns,<br />draw-only meshes and quality<br />statistics.</div>"]
        style 83 fill:#85bbf0,stroke:#5d82a8,color:#000000
      end

    end

    105-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    101-. "<div>uses material</div><div style='font-size: 70%'>[Rust paths]</div>" .->102
    101-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    101-. "<div>uses par, progress</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    101-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    83-. "<div>uses par</div><div style='font-size: 70%'>[Rust paths]</div>" .->106
    103-. "<div>uses edge, fem +2</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    103-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83
    103-. "<div>uses section</div><div style='font-size: 70%'>[Rust paths]</div>" .->104
    104-. "<div>uses props</div><div style='font-size: 70%'>[Rust paths]</div>" .->101
    104-. "<div>uses material</div><div style='font-size: 70%'>[Rust paths]</div>" .->102
    104-. "<div>uses mesh</div><div style='font-size: 70%'>[Rust paths]</div>" .->83

  end```
