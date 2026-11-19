# nyx architecture diagrams (Mermaid)

Hand-written from `workspace.dsl` on 2026-11-13, at the end of the restructure (phase 9). They are a hand-drawn summary with colour classes; the complete, exported and code-derived set (all four C4 levels, 91 components) is in [`c4/`](c4/README.md). When the DSL changes, change these in the same commit. They render in GitHub, VS Code (Markdown preview with a Mermaid extension) and mermaid.live.

The container diagram is exactly the layering table of `crates/nyx-engine/tests/layering.rs`: a dependency between crates that is not in the table fails the test. Colours: green = a physics (kernel) crate, blue = a front-end or the engine, grey = external, cylinder = files.

## 1. Context

```mermaid
flowchart LR
  classDef ext fill:#8a8a8a,stroke:#555,color:#fff
  classDef sys fill:#1168bd,stroke:#0b4884,color:#fff
  R([Researcher]):::ext
  AI[AI host<br/>Claude Code / Desktop]:::ext
  V[Results viewer<br/>ParaView, VisIt]:::ext
  W[(Workspace folder<br/>Documents/nyx)]:::ext
  N[nyx]:::sys
  R -->|edits models, runs studies| N
  R -->|asks in words| AI
  AI -->|MCP over stdio| N
  V -->|reads XDMF and HDF5| N
  N -->|workflows, runs, projects| W
```

## 2. Containers: the crates and the layering

```mermaid
flowchart TB
  classDef fe fill:#cfe3f7,stroke:#1168bd,color:#000
  classDef eng fill:#a9c9ee,stroke:#1168bd,color:#000
  classDef kern fill:#cfe8cf,stroke:#2e7d32,color:#000
  classDef lib fill:#e8eef5,stroke:#555,color:#000
  classDef file fill:#fff,stroke:#555,color:#000

  APP["app<br/>egui window (binary nyx)"]:::fe
  MCP["nyx-mcp<br/>MCP server, one tool per command"]:::fe
  CON["nyx-connect<br/>registers nyx-mcp with hosts"]:::fe
  ENG["nyx-engine<br/>Engine: project, commands, jobs,<br/>store link, results files"]:::eng
  PHY["nyx-physics<br/>Physics trait, ProblemSpec,<br/>SolveResult, schemas"]:::lib
  DFS["nyx-dfs<br/>heat, electrostatics"]:::kern
  SMS["nyx-sms<br/>solid mechanics"]:::kern
  EMS["nyx-ems<br/>EM waves"]:::kern
  ANA["nyx-analysis<br/>polarization"]:::lib
  STO["nyx-store<br/>workflows, runs, archives"]:::lib
  FEM["femcore<br/>meshes, meshing, adaptation,<br/>linalg, materials"]:::lib
  GEO["geom<br/>feature lists, selectors<br/>(serde + optional rayon)"]:::lib

  PF[("Project file<br/>.nyx.json, format 1")]:::file
  RF[("Results file<br/>.h5 + .xdmf")]:::file
  WS[("Workspace data<br/>workflows, runs, specs")]:::file
  ML[("Material library<br/>materials/*.toml")]:::file

  APP --> ENG
  APP --> CON
  MCP --> ENG
  ENG --> PHY
  ENG --> DFS
  ENG --> SMS
  ENG --> EMS
  ENG --> ANA
  ENG --> STO
  ENG --> FEM
  ENG --> GEO
  DFS --> PHY
  SMS --> PHY
  EMS --> PHY
  DFS --> FEM
  SMS --> FEM
  EMS --> FEM
  ANA --> PHY
  ANA --> FEM
  PHY --> FEM
  FEM --> GEO

  ENG -.-> PF
  ENG -.-> RF
  ENG -.-> ML
  STO -.-> WS
```

Rules the test enforces besides the arrows: `geom` has no workspace crate and no external dependency but `serde` and optional `rayon`; `nyx-store` and `nyx-connect` depend on no workspace crate; no physics depends on another physics, and only the engine names them.

## 3. Inside the engine

```mermaid
flowchart LR
  classDef c fill:#85bbf0,stroke:#1168bd,color:#000
  CMD[commands<br/>catalog + dispatcher]:::c
  PRJ[project<br/>format key checked at the file]:::c
  JOB[jobs<br/>mesh, solve, convergence]:::c
  MSH[meshing seam<br/>MeshRequest to MeshedModel]:::c
  MAT[materials<br/>library + project, at the wavelength]:::c
  REG[registry<br/>dfs, sms, ems compiled in]:::c
  ANL[analyses]:::c
  RIO[results files<br/>h5io, xdmf, open_results]:::c
  STL[store link<br/>versions, run table]:::c

  CMD --> PRJ
  CMD --> JOB
  CMD --> ANL
  CMD --> RIO
  CMD --> STL
  JOB --> MSH
  JOB --> MAT
  JOB --> REG
  JOB --> STL
  STL --> RIO
  MSH -->|MeshRequest| FM[femcore::meshing]
  MAT --> FMat[femcore::material]
  REG --> K[nyx-dfs / sms / ems]
  ANL --> NA[nyx-analysis]
  STL --> S[nyx-store]
```

## 4. femcore components

```mermaid
flowchart TB
  classDef c fill:#85bbf0,stroke:#1168bd,color:#000
  MESH[mesh<br/>native cells, exchange Mesh,<br/>discretisation]:::c
  MG[meshing<br/>lattice mesher, structured PML,<br/>wavelength sizing]:::c
  AD[adapt<br/>metric adaptation, estimators]:::c
  LA[linalg<br/>one driver: direct, Krylov-Schur,<br/>LOBPCG, AMG]:::c
  FE[fem, edge, pml, periodic, props]:::c
  MT[material<br/>TOML, models, library]:::c
  SEC[section<br/>the wafer frame]:::c
  EX[expr]:::c
  G[geom]

  MG --> MESH
  MG --> AD
  MG --> G
  AD --> MESH
  FE --> MESH
  FE --> SEC
  FE --> MT
  FE --> LA
  MT --> EX
```

## 5. A solve from the window (dynamic view)

```mermaid
sequenceDiagram
  participant W as app
  participant E as nyx-engine
  participant F as femcore
  participant K as nyx-ems
  participant S as nyx-store
  participant H as results.h5
  W->>E: start_solve(targets)
  E->>F: mesh the model (when none is current)
  E->>K: Physics::solve(mesh, spec, ctx)
  K->>F: assemble, shift-invert eigensolve
  K-->>E: SolveResult
  E->>S: freeze the spec, add a run to the table
  E->>H: results.h5 and results.xdmf in the run folder
  W->>E: poll: take the mesh and the results
```

## 6. An AI host builds and solves (dynamic view)

```mermaid
sequenceDiagram
  participant A as AI host
  participant M as nyx-mcp
  participant E as nyx-engine
  participant G as geom
  participant K as nyx-ems
  participant S as nyx-store
  A->>M: tools/call add_rect, add_overlay, add_physics
  M->>E: call(command, arguments)
  E->>G: feature list edit and evaluation
  A->>M: tools/call solve_run
  M->>E: call("solve_run")
  E->>K: mesh and solve
  E->>S: auto_store: the workflow shows in the dashboard
  E-->>M: summary (scalars, run id)
  M-->>A: result
```

## 7. Files and the contracts that pin them

```mermaid
flowchart LR
  classDef c fill:#fff,stroke:#555,color:#000
  P["engine_project.v2.json<br/>(project, format 1)"]:::c
  Q["problem_spec.v3.json<br/>(ProblemSpec, format 1)"]:::c
  R["solve_result.v2.json"]:::c
  D["physics_descriptors.v1.json<br/>physics_dfs / sms / ems .v1"]:::c
  K["engine_commands.v1.json<br/>mcp_tools.v3.json"]:::c
  M["mesh_cells, mesh_request,<br/>mesh_report"]:::c
  X["material_schema.v1.json<br/>material_kitchen_sink*"]:::c
  S["store_workflow.v1.json<br/>store_run.v1.json"]:::c
  H["results_file layout 3<br/>(VERSIONS.md, tests/h5.rs)"]:::c
  P --> Q
  Q --> D
  K --> P
  R --> H
  M --> H
```
