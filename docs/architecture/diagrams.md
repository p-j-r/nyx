# nyx architecture diagrams (Mermaid)

Hand-written from `workspace.dsl` on 2026-11-03, not exported by Structurizr (the CLI could not be run, see `README.md`). They render in GitHub, VS Code (Markdown preview with a Mermaid extension) and https://mermaid.live. When `workspace.dsl` changes, update these by hand or replace them with `structurizr-cli export -format mermaid`.

Colours: green = clean boundary or kernel, amber = duplicated or inconsistent (rough spot), red = no production caller, dashed arrow = coupling listed in `ROUGH_SPOTS.md` (the R-number is on the arrow or in the box).

## 1. Context

```mermaid
flowchart LR
  classDef ext fill:#8a8a8a,stroke:#555,color:#fff
  classDef sys fill:#1168bd,stroke:#0b4884,color:#fff
  R([Researcher]):::ext
  AI[AI host<br/>Claude Code / Desktop]:::ext
  V[Results viewer<br/>ParaView, optional HDF5]:::ext
  W[(Workspace folder<br/>Documents/nyx)]:::ext
  N[nyx]:::sys
  R -->|edits models, runs studies| N
  R -->|asks in words| AI
  AI -->|MCP over stdio| N
  N -->|HDF5/XDMF| V
  N -->|workflows, runs, projects| W
```

## 2. Containers

```mermaid
flowchart TB
  classDef ok fill:#cfe8cf,stroke:#2e7d32,color:#000
  classDef rough fill:#fff1cc,stroke:#b8860b,color:#000
  classDef core fill:#d9d2f2,stroke:#5b4b9a,color:#000
  classDef store fill:#e6e6e6,stroke:#777,color:#000
  R([Researcher]) --> APP
  AI([AI host]) --> MCP
  subgraph nyxws["nyx workspace (crates/)"]
    APP["app (egui)<br/>window, --batch, --script<br/>Project type lives here (R16)"]:::rough
    ST["nyx-store<br/>workflows, runs, .nyx archives<br/>card reads removed layout (R25)"]:::rough
    FC["femcore<br/>kernels, meshers, solvers, use cases,<br/>analysis, files — four jobs in one crate (R15)"]:::core
    GM["geom<br/>feature list, selectors<br/>serde only"]:::ok
  end
  subgraph aiws["ai workspace (separate Cargo workspace)"]
    MCP["nyx-mcp<br/>stdio JSON-RPC, 9 tools"]
    SES["nyx-session<br/>tools, hash handles, journal"]
    CON["nyx-connect<br/>copy also in app (R07)"]:::rough
  end
  PF[("Project file<br/>.nyx.json / .nyx.xml")]:::store
  WS[("Workspace data<br/>workflow.json, runs.jsonl")]:::store
  SS[("Session folder<br/>journal.jsonl")]:::store
  APP --> FC
  APP --> GM
  APP --> ST
  APP --- PF
  ST --- WS
  ST -.->|"xmlio + SolveResult only (R15)"| FC
  FC --> GM
  MCP --> SES
  SES --> FC
  SES --- SS
  CON -.->|registers| MCP
  APP -.->|"connect.rs duplicates it (R07)"| CON
```

## 3. femcore components

```mermaid
flowchart TB
  classDef kernel fill:#cfe8cf,stroke:#2e7d32,color:#000
  classDef rough fill:#fff1cc,stroke:#b8860b,color:#000
  classDef dead fill:#f8d7da,stroke:#b02a37,color:#000
  classDef mesh fill:#d9d2f2,stroke:#5b4b9a,color:#000
  classDef base fill:#d6e6fb,stroke:#1f5fa8,color:#000
  classDef misc fill:#e6e6e6,stroke:#777,color:#000

  SOLVE["solve()<br/>one entry point"]:::base
  PROB["problem types<br/>ProblemSpec, SolveResult"]:::base
  SEC["section<br/>wafer frame"]:::base
  subgraph kernels["Kernels: assemble -> solve -> fields (R26: 7 of 10 public)"]
    DFS[dfs]:::kernel
    SMS[sms]:::kernel
    EMS[ems]:::kernel
    SBS[sbs]:::kernel
  end
  OLD["original linear algebra<br/>sparse, amg, eigen, lobpcg (R02)"]:::rough
  PORT["solvers (ported library)<br/>examples only, no caller"]:::rough
  SM["SimplexMesh<br/>what kernels read"]:::rough
  CM["CellMesh (nyx)<br/>named element types"]:::rough
  MM["meshers mesh types<br/>second CellMesh, second SimplexMesh"]:::rough
  GD["geomdrive<br/>geom model to mesh"]:::mesh
  MG["meshing<br/>MeshPlan, sizing, loop"]:::mesh
  MP["meshers (ported)<br/>adaptation, PML, .fmsh"]:::mesh
  M3["mesher3d, refine3d<br/>tests only (R05)"]:::dead
  UC["use cases, chips<br/>app content (R15)"]:::misc
  AN["analysis, polarization<br/>own axis names (R28)"]:::misc
  FI["files: xmlio, h5io"]:::misc

  SOLVE --> PROB
  SOLVE --> DFS & SMS & EMS & SBS
  DFS & SMS & EMS & SBS --> SEC
  SBS --> EMS
  SBS --> SMS
  DFS & SMS & EMS --> OLD
  PORT -.->|intended replacement| OLD
  SOLVE --> SM
  GD --> SM
  MG --> GD
  MG --> MP
  MP <-.->|"interop (R01, R09)"| MM
  SM --- CM
  CM <-.-> MM
  M3 --> SM
  UC --> GD
  UC --> MG
  AN --> PROB
  FI --> PROB
```

## 4. A solve (dynamic view)

```mermaid
sequenceDiagram
  participant App as app worker
  participant S as solve()
  participant V as validate
  participant K as ems kernel
  participant Sec as section
  participant L as original linalg
  App->>S: solve_with(mesh, spec, state, progress)
  S->>V: check_inputs (order 1 or 2, no NaN)
  S->>K: dispatch Emw 2D
  K->>Sec: wafer axes and names
  K->>L: assemble, shift-invert eigensolve
  K-->>S: SolveResult (+ analysis)
  S->>V: check_result (no NaN or inf)
  S-->>App: result + stage log, or Err("cancelled")
```

## 5. An AI recipe (dynamic view)

```mermaid
sequenceDiagram
  participant H as AI host
  participant M as nyx-mcp
  participant S as nyx-session
  participant F as femcore
  participant J as session folder
  H->>M: tools/call recipe_run
  M->>S: dispatch
  S->>F: geom template -> mesh -> ProblemSpec -> solve
  S->>J: journal the call, write g_*.geom.json
  S-->>M: handles g_ m_ p_ r_ + assumption ledger
  M-->>H: result (structuredContent on 2025-06-18)
```

## 6. Packages the contracts pin

```mermaid
flowchart LR
  classDef c fill:#cfe8cf,stroke:#2e7d32,color:#000
  P[problem_spec.v1.json]:::c --- T1[femcore contracts.rs]
  F[frames] --- T1
  Mh[mesh: label bases, quad order, .fmsh v2] --- T1
  KS[kernel_stages: public fns] --- T1
  SR[solve_result structure] --- T1
  G[geom golden lists] --- T2[geom contracts.rs]
  SW[store_workflow.v1.json] --- T3[nyx-store contracts.rs]
  SRn[store_run.v1.json] --- T3
  SC["concept card (ignored, R25)"] -.- T3
  MT[mcp_tools.v1.json] --- T4[nyx-mcp contracts.rs]
```
