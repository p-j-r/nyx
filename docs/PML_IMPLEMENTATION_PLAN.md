# Implementation Plan: Physics-Controlled Structured PML & Open-Boundary Geometry, Native Mixed-Cell Meshing/Solvers, and Curvilinear Coordinate Stretching

- **Companion Specification:** [`docs/PML_FEATURE_REQUEST.md`](PML_FEATURE_REQUEST.md)
- **Target ADR:** ADR 128 (`docs/architecture/adr/0128-structured-curvilinear-pml-and-open-boundaries.md`)
- **Status:** Ready for Phased Implementation

---

## 0. Architectural Mandate & Repository Guardrails

Before editing any file, the implementing model **must** observe all `nyx` workspace and numerical invariants:

1. **Core Architectural Decision — Option B (Native Mixed-Cell `Quad` / `Prism` / `Hex` FEM Kernels, Amending `ADR 4`, Fulfilling `ADR 5` & `ADR 32`):**
   - Keep structured PML/IE cells (`Quad` in 2D; `Prism` on faces/shells and `Hex` on 3D edges, corners, and cylindrical rims—**zero transition pyramids**) **uncut** inside `CellMesh<D>` (`Mesh::D2` / `Mesh::D3`) all the way from `pml_box` / `pml_cyl` / `pml_sphere` through `MeshedModel` into `Physics::solve(&Mesh, ...)` (`crates/nyx-physics/src/physics.rs:101`).
   - Add `MixedDiscretisation` (`crates/femcore/src/mesh/mixed.rs`) and shared mixed-cell element/basis modules (`fem_mixed.rs`, `edge_mixed.rs`) so all physics crates assemble directly on native `Tri`+`Quad` (2D) and `Tet`+`Prism`+`Hex` (3D) meshes:
     - Reduces PML element count by **2×–6×**, slashes PML DOFs by **37%–65%**, and accelerates sparse `LBL^H` / `LU` factorization.
     - Eliminates oblique Nédélec edges and split-diagonal symmetry breaking ($\Delta f / f \sim 10^{-14}$).
     - Eliminates tangential locking under singular poles, unlocking **Rational PMLs** and **Rational Infinite Elements** at both $P_1$ and $P_2$.
   - Existing `Ams` (`crates/femcore/src/linalg/ams.rs:76`) and `AmsP2` (`ams.rs:429`) are already element-shape agnostic (consuming pure algebraic matrices `G`, `p1_to_p2`, and `edge_vectors`) and require **zero preconditioner rewrites**.
   - Keep `SimplexMesh` fast paths intact for pure-simplex meshes and decouple `view3d.rs` / `h5io.rs` via a Stage-3 display/HDF5 field projection helper.
2. **`jj` Workspace Isolation:**
   - Every coding session modifying files or running mutating builds must operate in its own dedicated Jujutsu workspace under `.jj-workspaces/<name>` (with `.jj-workspaces/` excluded locally in `.git/info/exclude`). Never commit directly to `main` or push.
3. **Dependency Layering (`crates/nyx-engine/tests/layering.rs`):**
   - `geom` (`serde` + optional `rayon` only) $\leftarrow$ `femcore` $\leftarrow$ `nyx-physics` $\leftarrow$ `nyx-{dfs, sms, ems, analysis}` (never depend on each other or `geom`) $\leftarrow$ `nyx-engine`.
4. **Release Builds, Licensing, Formatting, Contracts & Goldens:**
   - Always build and test with `--release`.
   - Every new Rust file in `femcore` and `nyx-*` crates **must** start with:
     ```rust
     // SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
     // Copyright (C) 2026 p-j-r
     //! <Module summary>
     ```
   - Selective formatting (`ADR 65`): run `cargo fmt` **only** on:
     `cargo fmt -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine`
     (never run global `cargo fmt` over `femcore` or `app`).
   - Golden tests (`crates/geom/golden/*.json`, `crates/nyx-{dfs,sms,ems}/tests/golden/`): **never** regenerate goldens (`NYX_WRITE_GOLDEN=1`) to mask regressions. New fields on `geom::model::Op` **must** use `#[serde(default, skip_serializing_if = ...)]` so `golden_feature_lists_round_trip_key_for_key` (`crates/geom/tests/contracts.rs:55-67`) passes unchanged.
   - Contract exemplars (`ADR 19`): whenever changing `geom_operations.v1.json`, `mesh_request.v1.json`, `mesh_report.v1.json`, `physics_{ems,sms,dfs}.v1.json`, or `engine_commands.v1.json`, regenerate with `NYX_WRITE_CONTRACT=1` and bump `docs/architecture/contracts/VERSIONS.md` in the same commit.
   - Architecture tracking: add every new file path/glob to `docs/architecture/CODE_TO_ADR.md` and verify `bash docs/architecture/check_code_to_adr.sh`.

---

## Phase 0 — `MixedDiscretisation` Topology, Uncut `CellMesh<D>` Flow & Physics-Controlled Defaults

**Goal:** Establish the canonical mixed-cell solver topology (`MixedDiscretisation`) without splitting structured PML cells, eliminate wasted 2D inline meshing, decouple `layer_grading` from stretch `order`, and make `PmlMode::Auto` the default.

### Files to Create / Modify
- **Create** `crates/femcore/src/mesh/mixed.rs` (SPDX header + `//!` summary)
- Modify `crates/femcore/src/mesh/mod.rs` (`pub mod mixed`, plus fix `from_cells_split_with` at `L520-703` for fallback/display)
- Modify `crates/femcore/src/mesh/cells.rs` (optional `face_nodes` or quadratic face-center helper if needed for $Q_2$ `Quad` / `CurvedHex` face centers)
- Modify `crates/femcore/src/mesh/exchange.rs` (expose native `CellMesh` discretisation `Mesh::discretise_mixed()`)
- Modify `crates/nyx-physics/src/lib.rs` (define `MixedDiscretisation` export + `discretise_mixed(&Mesh)`)
- Modify `crates/femcore/src/meshing/{geomdrive.rs, pml_geom.rs, pml_geom3.rs, request.rs, mod.rs}`
- Modify `crates/nyx-engine/src/{project.rs, meshing.rs, convergence.rs}`

### Exact Design & Algorithms

1. **`MixedDiscretisation` (`crates/femcore/src/mesh/mixed.rs`):**
   - Stores directly from `CellMesh<D>`:
     - `dim: usize`, `points: Vec<f64>` (`D` coords per vertex),
     - Cell CSR: `cell_types: Vec<CellType>`, `cell_offsets: Vec<u32>`, `cell_conn: Vec<u32>`, `cell_labels: Vec<i32>`,
     - Boundary facet CSR: `facet_types: Vec<FacetType>` (`Edge2`, `Tri3`, `Quad4`), `facet_offsets: Vec<u32>`, `facet_conn: Vec<u32>`, `facet_labels: Vec<i32>`, `exterior_facets: Vec<bool>`,
     - Unique sorted topological entities (deterministic lexicographic sort):
       - `edges: Vec<(u32, u32)>` (with `a < b`) + per-cell edge IDs & edge orientation signs (`+1` / `-1`),
       - `tri_faces: Vec<[u32; 3]>` (sorted vertex triple `a < b < c`) + per-cell tri-face IDs & $D_3$ permutation indices (`0..6`),
       - `quad_faces: Vec<[u32; 4]>` (canonical cyclic face starting at minimum vertex ID, oriented by second-minimum neighbor) + per-cell quad-face IDs & $D_4$ dihedral symmetry indices (`0..8`),
     - Curved geometry nodes: `edge_nodes: Vec<(u32, [f64; 3])>` (valid in **both** 2D and 3D!) and optional `quad_face_nodes: Vec<(u32, [f64; 3])>` (for curved 9-node `CurvedQuad`, 18-node `CurvedPrism`, 27-node `CurvedHex` centers), plus `qp_edges: Vec<u32>`.
   - **Zero-pyramid conformity verification helper (`validate_conforming`):** Verify in debug/test builds that every interior facet shared by two cells has identical vertex set (`Tri` $\leftrightarrow$ `Tri` across `Tet`/`Prism`, `Quad` $\leftrightarrow$ `Quad` across `Prism`/`Hex`) and zero mismatched tri-to-quad contacts.
   - **Keep `SimplexMesh::from_cells_split_with` fixed (`mesh/mod.rs:582-702`):** Even though production solvers use `MixedDiscretisation`, fix positive tet determinant orientation (`L641`), topological-parity 2D quad splitting (`L585`), and `edge_nodes` retention (`L695`) inside `from_cells_split_with` for legacy `SimplexMesh` callers, `pml3::add_pml_3d`, and display/export fallbacks.

2. **Keep `CellMesh<D>` Uncut Through `MeshedGeom` & `MeshedModel`:**
   - Update `MeshedGeom` (`crates/femcore/src/meshing/geomdrive.rs`) and `MeshedModel` (`crates/nyx-engine/src/meshing.rs:17-78`) so the primary stored mesh is `Mesh` (`CellMesh<2>` / `CellMesh<3>`) built by `structured_pml_native` / `structured_pml_3d_native`—never prematurely splitting `Quad`/`Prism`/`Hex` into simplices!
   - `MeshedModel::solver` (`OnceLock<Arc<SimplexMesh>>` at `meshing.rs:30`) computes the split `SimplexMesh` lazily *only* when requested by legacy simplex utilities or `h5io`/`view3d`.

3. **Skip 2D Inline PML Cells & Fix `sweep_model_3d`:**
   - In `mesh_geom_2d` (`geomdrive.rs`) and `meshing/mod.rs:116-133`: skip triangulation `Fill` points and Ruppert refinement inside PML/IE regions when `PmlMode::Structured` is active, matching 3D `Sizes3::skip`.
   - In `nyx-engine/src/meshing.rs:403-441`: wire `PmlMode::Structured` into `sweep_model_3d` (skipping PML regions during interior sweep and extruding lateral/endcap `Prism`+`Hex` PML layers around the swept physical box/cylinder).

4. **Physics-Controlled Default `PmlMode::Auto` & `convergence::refined`:**
   - Replace `structured_pml: bool` in `MeshRequest` (`request.rs:74`) and `MeshSettings` (`project.rs:49`) with `pml_mode: PmlModeRequest` (`Auto` [default] | `Structured` | `Inline`), deserializing legacy `"structured_pml": true/false` cleanly for backward compatibility.
   - Remove `pml_relax = 1.5` coarsening (`request.rs:40, 56`, `meshing.rs:359`) inside structured PMLs.
   - Wire `layer_grading` (default `2.0` for new models; `1.0` if omitted on legacy `Op::Pml`) into `PmlSpec::grading`, leaving stretch exponent `order` (`2.0`) strictly for coordinate stretching.
   - Report `element_counts: BTreeMap<String, usize>` (by `Tri/Quad/Tet/Prism/Hex`), `min_pml_scaled_jacobian`, and `pml_warnings` (`< 5` layers warning) in `MeshReport` (`request.rs:169-187`).
   - In `nyx-engine/src/convergence.rs:58-118`: scale absorbing-layer `layers` by `(layers / factor).ceil().max(5.0)` inside `refined()` when `scale_sizes` is enabled.

### Phase 0 Verification
- `cargo test --release -p femcore`
- Unit tests on `MixedDiscretisation`: verify canonical $D_3$/$D_4$ face-orientation consistency across `Tet`–`Prism` and `Prism`–`Hex` interfaces, positive Jacobian on split fallback, and zero mismatched interfaces (`validate_conforming`).

---

## Phase 1 — Native Mixed-Cell Finite Elements (`fem_mixed.rs`, `edge_mixed.rs`), Curvilinear Stretching & Physics Kernels

**Goal:** Implement shared high-efficiency tensor-product Lagrange and Nédélec basis modules in `femcore`, implement quadrature-point curvilinear $\mathbf{J}_s(\mathbf{x}_q)$ with both $\kappa$-Polynomial and Rational profiles, upgrade `nyx-ems` and `nyx-sms` kernels natively to `MixedDiscretisation`, and add Stage-3 `SimplexMesh` display projection.

### Files to Create / Modify
- **Create** `crates/femcore/src/fem_mixed.rs` (SPDX header + `//!` summary)
- **Create** `crates/femcore/src/edge_mixed.rs` (SPDX header + `//!` summary)
- Modify `crates/femcore/src/{lib.rs, props.rs, pml.rs}`
- Modify `crates/femcore/src/meshing/physics_size.rs` (shared Hellmann–Feynman Christoffel helper `backward_wave_mpml_ratio`)
- Modify `crates/nyx-ems/src/{lib.rs, problem.rs, cavity.rs, bloch.rs, driven.rs, modes.rs, modes_complex.rs, axisym.rs}`
- Modify `crates/nyx-sms/src/{lib.rs, problem.rs, solid.rs, pml.rs, axisym.rs, waveguide.rs, lrb.rs}`
- Modify `crates/nyx-physics/src/{desc.rs, result.rs}`

### Exact Element Spaces, Equations & Algorithms

1. **Shared Tensor-Product Lagrange & Curved Elements (`crates/femcore/src/fem_mixed.rs`):**
   - **2D `Quad`** (reference $[-1, 1]^2$, separable tensor Gauss–Legendre rules):
     - $Q_1$ (4 nodes), serendipity $Q_2'$ (8 nodes—default for $P_2$ Lagrange in `nyx-sms` & `nyx-dfs`: 25% fewer DOFs per quad), full biquadratic $Q_2$ (9 nodes—exact de Rham scalar partner for $N_2$ edge elements).
   - **3D `Prism`** (reference $\hat{T} \times [-1, 1]$, separable Dunavant $\times$ Gauss–Legendre rules):
     - $P_1 \otimes Q_1$ (6 nodes), 15-node serendipity (6 vertices + 9 edge midpoints—37.5% fewer DOFs for `nyx-sms` & `nyx-dfs`), 18-node $P_2 \otimes Q_2$ (6 vertices + 9 edge midpoints + 3 quad-face centers—exact de Rham partner for $N_2$).
   - **3D `Hex`** (reference $[-1, 1]^3$, separable 1D Gauss–Legendre tensor rules):
     - $Q_1$ (8 nodes), 20-node serendipity (8 vertices + 12 edge midpoints—50% fewer DOFs for `nyx-sms` & `nyx-dfs`), 27-node $Q_2$ (exact de Rham partner for $N_2$).
   - **Curved Isoparametric Mappings:** `CurvedQuad` (8/9-node), `CurvedPrism` (15/18-node), `CurvedHex` (20/27-node) computing exact $\mathbf{x}(\hat{\boldsymbol{\xi}})$, $\mathbf{J}_K(\hat{\boldsymbol{\xi}})$, $\det\mathbf{J}_K$, and $\nabla\phi_a = \mathbf{J}_K^{-T}\nabla_{\hat{\xi}}\hat{\phi}_a$.
   - **DOF Manager (`MixedLagrangeDofMap`):** Maps vertices (`0..n_v`), unique edges (`n_e`), and quad-face centers (`n_qf`, when tensor $Q_2$ is used) seamlessly across `Tri`+`Quad` and `Tet`+`Prism`+`Hex` interfaces so $C^0$ continuity holds automatically.

2. **Shared Conforming Mixed-Cell Nédélec $H(\mathrm{curl})$ Spaces & Exact Algebraic `AMS` Operators (`crates/femcore/src/edge_mixed.rs`):**
   - **2D `Quad`:**
     - $N_1$ ($Q_{0,1}\times Q_{1,0}$, 4 edge DOFs; tangential trace matches `Tri` $N_1$ edge DOFs identically).
     - $N_2$ ($Q_{1,2}\times Q_{2,1}$, 8 edge + 4 interior DOFs; tangential edge trace matches `Tri` $N_2$ 2-DOF/edge space via canonical edge orientation sign).
   - **3D `Prism`** (Nédélec first-kind on $\hat{T} \times [-1, 1]$):
     - $N_1$ (9 DOFs: 6 triangle-base/top edge DOFs + 3 vertical edge DOFs; 43% fewer DOFs/vertex than 3 split tets).
     - $N_2$ (40 local DOFs: 18 edge [2/edge] + 4 tri-face [2/tri] + 12 quad-face [4/quad] + 6 interior):
       - **Exact trace compatibility:** Base/top tri faces share the exact 2-DOF tangential space $\mathbb{P}_1(\hat{T})\times$ bubble of `Tet` `Ned3` (`cavity.rs:48-215`) via $D_3$ face permutation; side quad faces share the exact 4-DOF $Q_{1,2}\times Q_{2,1}$ tangential space of `Hex` via $D_4$ dihedral sign/axis permutation!
   - **3D `Hex`** (Nédélec first-kind $Q_{k-1,k,k}\times Q_{k,k-1,k}\times Q_{k,k,k-1}$):
     - $N_1$ (12 edge DOFs; 57% fewer DOFs/vertex than 6 split tets).
     - $N_2$ (54 local DOFs: 24 edge + 24 quad-face + 6 interior; 37% fewer DOFs/vertex than split tets).
   - **Covariant Piola Mapping:** $\mathbf{w}_a(\hat{\boldsymbol{\xi}}) = \mathbf{J}_K^{-T}\hat{\mathbf{w}}_a(\hat{\boldsymbol{\xi}})$ and $\nabla\times\mathbf{w}_a(\hat{\boldsymbol{\xi}}) = \frac{1}{\det\mathbf{J}_K}\mathbf{J}_K(\hat{\nabla}\times\hat{\mathbf{w}}_a(\hat{\boldsymbol{\xi}}))$ (handles both straight and curved `CurvedPrism`/`CurvedHex`/`CurvedQuad` identically!).
   - **Exact Discrete Gradient & $P_1 \to P_2$ Transfer (`gradient_p1`, `gradient_p2`, `p1_to_p2`):**
     - Build exact CSR `gradient_p1` (Nédélec $N_1$ edges $\times$ $Q_1/P_1$ vertices), `gradient_p2` (Nédélec $N_2$ DOFs $\times$ full $Q_2$/$P_2$ scalar DOFs—including quad-face center nodes on `Prism`/`Hex`!), and `p1_to_p2` ($N_2 \times N_1$ inclusion) for `MixedDiscretisation`.
     - Pass directly into `Ams::new` and `AmsP2::new` (`crates/femcore/src/linalg/ams.rs:76, 429`): **zero changes to `ams.rs`**, guaranteed exact $\mathbf{K}\mathbf{G}_1 = \mathbf{0}$ and $\mathbf{K}\mathbf{G}_2 = \mathbf{0}$ across `Tet`–`Prism`–`Hex` interfaces!

3. **Coordinate-Free Stretching, Dual Profiles & Auto-Scaling (`crates/femcore/src/{props.rs, pml.rs}`):**
   - Extend `femcore::props::Pml` (all new fields `#[serde(default)]`):
     - `coord: PmlCoord` (`Cartesian`, `Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness }`, `Spherical { center, r_in, thickness }`),
     - `bounds: Option<([f64; 3], [f64; 3])>` (explicit Cartesian domain box),
     - `profile: PmlProfile` (`KappaPolynomial { order: f64, kappa_max: f64 }` vs `Rational { curvature: f64 }`),
     - `auto_scale: bool`, `target_reflection: f64` (`1e-6`),
     - `stabilization: PmlStabilization` (`Off | Auto | Mpml { ratio } | Adiabatic { eta_max }`),
     - `clamp_outer: bool` (default `true` for `nyx-sms`).
   - Implement `pml::auto_scale_params(...)` per [`PML_FEATURE_REQUEST.md` §2 Verdict 4](PML_FEATURE_REQUEST.md) (jointly calibrating multi-material groups to $(\lambda_{\min}, \lambda_{\max})$ for $C^0$ map continuity across interfaces).
   - Implement `pml::pml_jacobian(spec, bounds, x, dim, axisym) -> PmlJacobian { j, j_inv, det, s_theta, s_dist }`:
     - Supports both `KappaPolynomial`: $s(t) = 1 + (\kappa_{\max}-1)t^m - i\psi_{\max} t^m$ and **Rational**:
       $$s(t) = \frac{1}{(1 - c_{\text{rat}} t)^2} - i\,\psi_{\text{rat}}\frac{t^m}{(1 - c_{\text{rat}} t)^2}, \qquad \tilde{u}(u_0, L, t) = u_0 + L\left(\frac{t}{1 - c_{\text{rat}} t} - i\psi_{\text{rat}}\int_0^t \frac{\xi^m}{(1 - c_{\text{rat}}\xi)^2}\,d\xi\right)$$
       (default $c_{\text{rat}} = 0.85\text{--}0.95$ calibrated to target reflection $R_0$ and study type, evaluated analytically for $\tilde{r}(r)$ and $\tilde{R}(R)$).
     - Exact coordinate-free Jacobians for Cartesian, Cylindrical, and Spherical (regular at $r \to 0^+$ in 2D-axi: $s_\theta = \tilde{R}(R)/R = s_\Omega(R)$).
   - Implement `pullback_maxwell(eps: C3, mu_inv: C3, jac: &PmlJacobian) -> (C3, C3)`.

4. **`nyx-ems` & `nyx-sms` Native Mixed-Cell Assembly & Kernel Upgrades:**
   - **All EM kernels (`cavity`, `bloch`, `driven`, `modes`, `modes_complex`, `axisym`):**
     - Dispatch cell loops over `(Tet/Prism/Hex)` in 3D and `(Tri/Quad)` in 2D using `MixedDiscretisation`, `fem_mixed`, and `edge_mixed` (preserving pure-simplex `SimplexMesh` fast loops when no `Quad`/`Prism`/`Hex` cells exist).
     - Widen `mu_inv` from `[C; 3]` to symmetric `[[C; 3]; 3]` (`C3`).
     - Evaluate `pml_jacobian(..., x_q)` and `pullback_maxwell` **at each Gauss quadrature point $\mathbf{x}_q$** inside PML cells.
     - Keep `driven.rs` (`ADR 106`) PROM strictly $\omega$-independent when $f_{\max}/f_{\min} \le 3$ (auto octave sub-banding above $3\times$).
     - In `nyx-ems/src/axisym.rs`: allow meshes touching $r = 0$, separate $r = 0$ symmetry facets from outer PEC walls, apply $m$-dependent axis regularity ($m=0$: $E_\theta=0$, $E_z$ free; $|m|=1$: $E_z=0$, $E_\theta$ free; $|m|\ge 2$: both $0$), and apply full $2\times 2$ $(r, z)$ in-plane coupling for Spherical PMLs.
   - **All SMS kernels (`solid`, `pml`, `axisym`, `waveguide`, `lrb`):**
     - Assemble on `MixedDiscretisation` using serendipity 8-node `Quad` / 15-node `Prism` / 20-node `Hex` at order 2 (and $Q_1$/$P_1\otimes Q_1$/$Q_1$ at order 1), applying full $\mathbf{J}_s(\mathbf{x}_q)^{-T}\nabla\phi_a$ and $\det\mathbf{J}_s(\mathbf{x}_q)$ at quadrature points.
     - Run pre-solve Hellmann–Feynman Christoffel check (`backward_wave_mpml_ratio` in `physics_size.rs`) and wire `Auto` / `Mpml` / `Adiabatic` stabilization.
     - Fix `f.im.abs()` sign masking (`nyx-sms/src/pml.rs:232`, `lrb.rs:292`) to flag `f.im < -1e-9 * f.re` amplifying modes explicitly.
     - Clamp outer PML exterior facets (`u = 0`) by default (`clamp_outer: true`).

5. **Stage-3 Field Projection for `view3d.rs` & `h5io.rs`:**
   - Add a lightweight `MixedToSimplexProjector` in `nyx-physics` / `femcore` that maps nodal fields directly by vertex ID (`0..n_vertices`) and projects per-cell $P_1$/$P_2$ `CellField` values from `Quad`/`Prism`/`Hex` onto the `SplitMode::Minimal` sub-triangles/sub-tets of `MeshedModel::solver` (or compute `layer fraction N` and physical `CircuitRom` integrals directly on native cells in Stage 3 while projecting visualization fields onto `SimplexMesh` for `h5io.rs`, `xdmf.rs`, and `view3d.rs`).

### Phase 1 Verification
- **De Rham identity test:** Assemble $K$ and $G_1, G_2$ on mixed `Tet + Prism + Hex` meshes and assert $\|K G_1\|_\infty < 10^{-13}$ and $\|K G_2\|_\infty < 10^{-13}$.
- **`Ams` & `AmsP2` solver test:** Verify preconditioned COCG / FGMRES convergence on mixed `Tet + Prism + Hex` curl-curl problems matches or beats split-tet baselines.
- **`nyx-ems::axisym` empty sphere ($r=0$) test:** Match analytical $m=0, 1, 2$ spherical-Bessel eigenfrequencies within $< 0.1\%$.
- **`nyx-sms` backward-wave test:** Tilted R-plane sapphire / `LiNbO3` triggers $p_{\min} > 0$, raw PML flags `f.im < 0`, and `mpml`/`adiabatic` eliminates amplifying modes cleanly.

---

## Phase 2 — Outward Boundary-Extruded Absorbing Layers, Legacy `Op::Pml` Adapter & Ray-Invariance Guard (`crates/geom`)

**Goal:** Provide the new zero-boilerplate outward absorbing-layer geometry model across all 7 topologies, keep legacy inward `Op::Pml` 100% backward compatible, and enforce ray-invariance validation.

### Files to Modify
- `crates/geom/src/{model.rs, sim.rs, scene.rs, planar.rs, stack.rs, select.rs}`

### Exact Design & Algorithms

1. **Extend `geom::model::Op::Pml` (`crates/geom/src/model.rs:295-302`) & Optional `Op::AbsorbingLayer` Alias:**
   Add serde-defaulted fields (`#[serde(default, skip_serializing_if = "...")]`) so all existing `crates/geom/golden/*.json` serialize identically byte-for-byte:
   - `mode: PmlGrowMode` (`Inward` [when omitted in legacy JSON] | `Outward` [default for new UI/script helpers])
   - `coord: PmlCoordKind` (`Auto | Cartesian | Cylindrical | Spherical`)
   - `center: Option<[Expr; 3]>`, `axis: Option<Axis>`
   - `layer_grading: Option<Expr>`
   - Support `"auto"` / `None` sentinel on `thickness` and `layers` so `Op::Pml` can run completely physics-driven without hardcoded CAD dimensions, and make `window` optional in `Outward` mode (defaulting to the outer envelope of the `Overlay` when omitted—eliminating the dummy `"window"` rectangle/box entirely!).
   - *(Optional clean ergonomics)* Expose `add_absorbing_layer(...)` in `nyx-engine::geomdesc` / `commands.rs` / `nyx-script` that creates an `Op::Pml { mode: Outward, window: None, ... }`, making the new workflow unmistakable while keeping `add_pml` working identically for legacy scripts.

2. **Enrich `PmlRegion` (`crates/geom/src/sim.rs:10-22`):**
   Add `coord: PmlCoordKind`, `mode: PmlGrowMode`, `center: [f64; 3]`, `axis: [f64; 3]`, `r_in: f64`, `z_range: [f64; 2]`, `layer_grading: f64`, `role: PmlRole` (`Face | Edge | Corner | RadialShell | EndCap | CornerAnnulus | SphericalShell`).

3. **Auto-Detection, All 7 Topologies & Standalone `Sphere` Host (`scene.rs`, `stack.rs`):**
   - `coord == Auto` infers `Cartesian` (`Rect`/`Box`), `Cylindrical` (2D `Circle` off-axis or 3D `Cylinder`), or `Spherical` (2D `Circle` centered at $x=0$ in $x \ge 0$, or 3D `Sphere`).
   - 2D-axisymmetric `Rect` touching $x = 0$ auto-omits `XMin`.
   - 3D `Cylinder` inward/outward supports radial shell + 2 axial endcaps + 2 corner annuli (all prismatic along the stack axis).
   - Lift `stack.rs:261-263` prismatic-host error for standalone 3D `Sphere` outer hosts (supporting interior primitives strictly inside $R_{\text{in}}$) and emit `Surface::Sphere` for `mesher3d`.

4. **Ray-Invariance Validation Guard (`geom` + `femcore::meshing`):**
   - Before meshing/solving, inspect materials and `Op::Sheet` edges touching $\partial\Omega_{\text{phys}}$:
     - Reject `Spherical` if $> 1$ horizontal material layer or an off-center dielectric/metal sheet intersects the outer sphere across varying polar angles $\vartheta$ ($\partial_R\varepsilon \neq 0$), reporting why `Cylindrical` (layered wafer) or `Cartesian` (piercing waveguide) must be used.
     - Reject `Cylindrical` if an extruded waveguide/strip crosses the barrel wall $r = R_{\text{in}}$ obliquely ($\partial_r\varepsilon \neq 0$), recommending `Cartesian`.

### Phase 2 Verification
- `cargo test --release -p geom` (verify 100% golden round-trip preservation + new unit tests for all 7 topologies, outward mode without dummy `window`, and ray-invariance guard diagnostics).

---

## Phase 3 — Native Mixed-Cell Curvilinear Meshers (`pml_cyl.rs`, `pml_sphere.rs`) & Lockstep Curved Shells

**Goal:** Build native `Quad` (2D) and `Prism` + `Hex` (3D) structured sweepers for Cartesian, Cylindrical, and Spherical topologies with zero pyramids, material/`Op::Sheet` continuation, and lockstep curved concentric shells.

### Files to Create / Modify
- **Create** `crates/femcore/src/meshing/pml_cyl.rs` (SPDX header + `//!` summary)
- **Create** `crates/femcore/src/meshing/pml_sphere.rs` (SPDX header + `//!` summary)
- Modify `crates/femcore/src/meshing/{mod.rs, pml_box.rs, pml_geom.rs, pml_geom3.rs, sheet3d.rs, curve3.rs}`

### Exact Algorithms

1. **`crates/femcore/src/meshing/pml_box.rs` (Cartesian Outward & Inward Upgrade):**
   - Keep `build_pml` emitting pure native `Quad` (2D) and `Prism` + `Hex` (3D) layers cleanly with deterministic canonical `(layer_k, sym_key)` vertex ordering.
   - Extend material inheritance and `Op::Sheet` edge extrusion so exterior sheet edges on $\partial\Omega_{\text{phys}}$ become interior conductive `Quad` facets continuing through the lateral `Prism` layers!

2. **`crates/femcore/src/meshing/pml_cyl.rs`:**
   - **`build_pml_annulus_2d`:** Extrude boundary edges on $r = R_{\text{in}}$ radially into $N_r$ rings of `CellType::Quad` (inheriting base facet material `<mat>.pml_radial`), with lockstep `CurvedQuad` mid-edge nodes on every circle $r_k$.
   - **`build_pml_cylinder_3d`:**
     - Barrel triangles ($r = R_{\text{in}}$) $\xrightarrow{\text{sweep along } \hat{\mathbf{r}}}$ $N_r$ layers of `CellType::Prism` (preserving $z$-stratified wafer materials per base triangle!).
     - Bottom/top endcap triangles ($\zeta = z_{\min}, z_{\max}$) $\xrightarrow{\text{sweep along } \mp\hat{\mathbf{a}}}$ $N_z$ layers of `CellType::Prism`.
     - Circular rim edge chains ($r = R_{\text{in}}, \zeta = z_{\min}/\max$) $\xrightarrow{\text{sweep across } (r_k, \zeta_j)}$ $N_r \times N_z$ rings of `CellType::Hex` (corner annuli), inheriting material from the barrel facet touching each rim edge.
     - Extrude interior `Op::Sheet` edges touching $r = R_{\text{in}}$ radially into conductive `Quad` sheet facets between `Prism` layers.
     - **Lockstep cylindrical curving (`edge_nodes` / `quad_face_nodes`):** For every shell $k = 0 \dots N_r$, place circumferential mid-edge nodes on exact circle $r_k$, and place quad-face center nodes on curved barrel quads when `Q_2` `CurvedPrism`/`CurvedHex` is active.

3. **`crates/femcore/src/meshing/pml_sphere.rs`:**
   - **`build_pml_semicircle_2d` (2D Axisymmetric Spherical):** Extrude boundary edges on $R = R_{\text{in}}$ ($r \ge 0$) radially into $N_R$ rings of `CellType::Quad`, clamping $\vartheta = 0, \pi$ endpoints strictly to $r = 0.0$ and tagging extruded axis edges with the `r=0` symmetry facet label (never outer PEC), with lockstep `CurvedQuad` mid-edge projection.
   - **`build_pml_sphere_3d`:** Sweep exterior triangles on $R = R_{\text{in}}$ radially along $\hat{\mathbf{R}}_v$ into $N_R$ concentric shells of `CellType::Prism` (**zero pyramids**), with lockstep `CurvedPrism` mid-edge projection across all shells $R_0 \dots R_{N_R}$.

4. **Unified Dispatch & Curving Isolation (`pml_geom.rs`, `pml_geom3.rs`, `mod.rs`, `curve3.rs`):**
   - Route both `Outward` and legacy `Inward` PMLs to `pml_box` / `pml_cyl` / `pml_sphere` returning uncut `CellMesh<D>`.
   - Ensure `curve3::curve_solver_mesh` and `sheet3d_curve::apply_sheet3d_curved_edges` operate strictly on interior non-PML faces/edges, leaving lockstep PML `edge_nodes` intact (no thin-layer sagitta inversion).

### Phase 3 Verification
- Topology, conformity (`validate_conforming`), and scaled Jacobian checks across all 7 topologies: zero pyramids, zero hanging edges/faces, `min scaled Jacobian > 0.3` at $h_R / R_{\text{in}} = 0.01$ with `curved = true`, and `sDist` bounds $\in [0, 1]$.

---

## Phase 4 — Engine, Script, MCP, GUI & Examples Parity

**Goal:** Expose the zero-boilerplate outward PML workflow, profile/coordinate choices, diagnostics, and one-click physics wiring consistently across all three `nyx` front-ends.

### Files to Modify
- `crates/nyx-engine/src/{pml.rs, meshing.rs, geomdesc.rs, commands.rs, preview.rs, guide.rs}`
- `crates/nyx-engine/cards/geometry.md` (must remain `< 3600` bytes per `guide.rs:31`)
- `crates/nyx-script/src/{prelude.rs, stubs.rs}`
- `crates/app/src/{geomui.rs, tabs.rs}`
- `examples/{09_pml_slab, 11_anchor_loss, 13_silicon_on_sapphire, 17_disk_anchor_loss}.nyx.json` + new canonical examples

### Exact Changes
1. **`nyx-engine::pml::add_pml_conditions`:** Support both legacy inward `Op::Pml` and new outward absorbing layers, auto-calibrating multi-material groups jointly from study `Sizing` (`auto_scale = true`) and setting `profile = "kappa_polynomial"` (eigenfrequency/iterative) or `"rational"` (direct driven/radiation) automatically.
2. **Pre-solve `sDist` & Post-solve $P_{\text{pml}}$ Diagnostics:** Expose `pml_diagnostics` command verifying `sDist ∈ [0, 1]`, `element_counts_by_type`, `min_pml_scaled_jacobian`, and Christoffel $p_{\min}$; report radiated power $P_{\text{pml}}$ in `nyx-ems::driven` and `nyx-sms` scalar outputs alongside `layer fraction N`.
3. **Preview, Cards & Scripts:** Render analytical outward PML shells as translucent wireframe envelopes in preview (`preview.rs`, `ADR 92`); add `add_absorbing_layer` / updated `add_pml` + `set_mesh(pml_mode=...)` to `commands.rs` and `.nys` stubs; update `cards/geometry.md` within its `< 3600` byte budget.
4. **GUI (`crates/app/src/{geomui.rs, tabs.rs}`):**
   - Default new absorbing layers to `mode = Outward`, `coord = Auto`, `layers = 6` (`Rational`) / `8` (`KappaPolynomial`), `layer_grading = 2.0`.
   - Replace `"structured absorbing layers"` checkbox (`tabs.rs:285`) with a 3-way selector (`Auto [default] | Structured | Inline`).
   - Add one-click **"Add PML Conditions from Geometry"** button in the Physics tab (`tabs.rs:918-1004`).
5. **Examples:** Keep all legacy tests working, modernize `09`, `11`, `13`, `17` to outward structured `Quad`/`Prism`/`Hex` setups, and ship canonical 2D Cylindrical, 2D Axisymmetric Spherical, and 3D Cylindrical/Spherical PML examples.

### Phase 4 Verification
- `cargo build --release -p nyx-script && cargo test --release -p nyx-mcp`
- `cargo test --release -p nyx-engine`

---

## Phase 5 — 3-Tier Open-Boundary Completion in `nyx-dfs` (& `nyx-sms::static`)

**Goal:** Deliver complete open-boundary modeling for electrostatics, stationary heat, magnetostatics, and static elasticity using native `Quad`/`Prism`/`Hex` Rational & Exponential Infinite Elements alongside `far_field` Robin BCs.

### Files to Modify
- `crates/femcore/src/pml.rs` (add `ie_rational_jacobian` & `ie_exponential_jacobian`)
- `crates/nyx-dfs/src/{heat.rs, magneto.rs, electro.rs, problem.rs, lib.rs}`
- `crates/nyx-sms/src/{solid.rs, problem.rs, lib.rs}`

### Exact Changes
1. **Tier 2 (`far_field` Robin extensions):**
   - **`nyx-dfs/src/heat.rs`:** Register `"far_field"` (`Applies::Boundary`) on `heat_stationary`, reusing `femcore::farfield::assemble_farfield_robin` with thermal conductivity $\kappa$ around $T_\infty$.
   - **`nyx-dfs/src/magneto.rs`:** Add 1st-order asymptotic curl-Robin BC (`"far_field"`):
     $$\int_\Gamma \frac{1}{\mu_0 r}(\hat{\mathbf{n}}\times\mathbf{A}_a)\cdot(\hat{\mathbf{n}}\times\mathbf{A}_b)\,dS$$
     for `solve_mixed_3d` (native `Tet`+`Prism`+`Hex` boundary faces) and 2D cross-section, eliminating PEC wall flux-compression bias (`magneto.rs:90`).
2. **Tier 3 (Outward Native `Quad`/`Prism`/`Hex` Infinite Elements for Layered 3D Chips & `sms::static`):**
   - Register `"infinite_element"` (`Applies::Domain`) on `nyx-dfs::{electrostatic, magnetostatic, heat_stationary}` and `nyx-sms::static` using native mixed-cell assembly (`fem_mixed` / `edge_mixed`).
   - Support **Rational Infinite Elements** ($s_u(t) = \frac{d_{\text{pole}}(d_{\text{pole}}+d)}{(d_{\text{pole}}+d - t d)^2}$ mapping $[R_{\text{in}}, R_{\text{out}}) \to [R_{\text{in}}, \infty)$ in 5–6 layers—**locking-free** on `Quad`/`Prism`/`Hex` for monopole, dipole, and multipole fields!) and optional **Bounded Exponential + Order-1 `far_field` Robin** ($\chi = e^\beta \in [20, 50]$ for tighter `AMG`/`PCG` conditioning).
   - **Transient Heat Guard:** Explicitly reject `"infinite_element"` on `heat_transient` with an error recommending unstretched ($s_u \equiv 1$) geometrically graded diffusion-length layers $d = \sqrt{4 D_{\max} t_{\text{end}} \ln(1/R_0)}$.

### Phase 5 Verification
- Isolated 3D sphere capacitance (`nyx-dfs::electrostatic`): `< 0.1%` error at $R_{\text{in}} = 2R$.
- Circular loop inductance (`nyx-dfs::magneto`): curl-Robin reduces box-wall bias $> 10\times$ vs PEC walls.
- 5-layer native `Prism`/`Hex` **Rational Infinite Element** monopole, dipole, and layered-chip thick-metal benchmarks: `< 0.1%` error with zero tangential locking at both $P_1$ and $P_2$.
- `heat_transient` guard test confirming `infinite_element` is cleanly rejected.

---

## Phase 6 — End-to-End Verification Matrix, Contracts & ADR 128

### Verification Suite (All Run in `--release`)

1. **Option B Core & Solver Efficiency Benchmarks:**
   - **De Rham Exactness:** $\|K G_1\|_\infty < 10^{-13}$ and $\|K G_2\|_\infty < 10^{-13}$ on mixed `Tet + Prism + Hex` meshes.
   - **`Ams` & `AmsP2` Parity:** Verify `Ams` and `AmsP2` converge on mixed `Tet + Prism + Hex` meshes in $\le$ the iteration count of split-tet meshes.
   - **DOF & Factorization Benchmark:** Measure and log $\ge 35\%$ PML DOF reduction and `LBL^H` factorization speedup on 3D cavity/driven benchmarks relative to split-tet `Option A`.
   - **Exact Symmetry:** Degenerate doublet splitting $|\Delta f / f| < 10^{-12}$ from PML on $C_{4v}$ / circular symmetric cavities/waveguides.
2. **Physics & Topology Benchmarks (All 7 Topologies):**
   - **2D Cartesian & Cylindrical PML:** $R < 10^{-6}$ ($\kappa$-Polynomial) and $R < 10^{-8}$ (Rational, 5–6 layers).
   - **2D Axisymmetric Spherical ($r = 0$ touching) & Cylindrical PML:** Mie / dipole $f$ and $Q$ for $m = 0, 1, 2$ within $< 0.5\%$ of analytical values.
   - **3D Cartesian, Cylindrical (layered wafer $\partial_r\varepsilon = 0$), & Spherical Curved-Shell PML:** `min scaled Jacobian > 0.3` at $h_R/R_{\text{in}} = 0.01$; $P_2$ cavity/driven $R < 10^{-6}\text{--}10^{-8}$; ray-invariance guard rejects Spherical PML on layered wafer; `driven` PROM exact for $f_{\max}/f_{\min} \le 3$.
   - **`nyx-sms` Backward Waves & Clamping:** Tilted `Sapphire` / `Si [110]` and isotropic $S_1$ Lamb membrane show zero amplifying `f.im < 0` modes under `"auto"`/`"mpml"`/`"adiabatic"`, and zero rigid-body modes when outer-clamped.
   - **Determinism & Convergence:** Bit-identical results across `NYX_THREADS=1` vs `4`; `convergence_check` scales PML `layers` and converges monotonically.
3. **Mandatory Repository Gate Commands:**
   ```powershell
   cargo test --release --workspace
   cargo build --release -p nyx-script
   cargo test --release -p nyx-mcp
   cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering
   cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
   cargo fmt -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine -- --check
   bash docs/architecture/check_code_to_adr.sh
   ```
4. **ADR & Contract Finalization:**
   - Write `docs/architecture/adr/0128-structured-curvilinear-pml-and-open-boundaries.md` (explicitly documenting amendment of `ADR 4` via `MixedDiscretisation`, fulfillment of `ADR 5`/`ADR 32`, outward extrusion, and curvilinear PMLs/IEs).
   - Register ADR 128 in `docs/architecture/ADR_INDEX.md` and update `docs/architecture/CODE_TO_ADR.md` rows for `mixed.rs`, `fem_mixed.rs`, `edge_mixed.rs`, `pml_cyl.rs`, `pml_sphere.rs`, `docs/PML_FEATURE_REQUEST.md`, and `docs/PML_IMPLEMENTATION_PLAN.md` to cite ADR `128`.
   - Regenerate modified contract exemplars (`NYX_WRITE_CONTRACT=1`) and bump `docs/architecture/contracts/VERSIONS.md` in the same commit.

---

## Stop & Escalate / Re-Plan Triggers

The implementing model must stop and re-plan (do not improvise broad workarounds) if any of the following occur:

1. **$D_4$ / $D_3$ Trace Orientation Mismatch in Mixed $N_2$ Nédélec (`edge_mixed.rs`):** If $\|K G_2\|_\infty \ge 10^{-12}$ or `AmsP2` stalls across `Tet`–`Prism` or `Prism`–`Hex` interfaces, stop immediately and verify edge/face tangential orientation canonicalization in `mixed.rs` against `Ned3` (`cavity.rs:48-215`) before proceeding to physics sweeps.
2. **Standalone 3D `Sphere` Host Conflicts in `geom::stack`:** If lifting `stack.rs:261-263` for multi-material interior solids inside a spherical shell would require general 3D non-prismatic CSG booleans, stop and restrict 3D Spherical PML to `mode = "outward"` extrusion around a physical `Sphere` (or `Box`/`Cylinder`) without adding 3D CSG booleans to `geom`.
3. **Iterative Preconditioner Conditioning Under Rational Stretching:** If direct solves (`LBL^H`) run cleanly on Rational PMLs/IEs while `AMS`/`AMG` iterative solves stall at $s_{\max} \gg 10^2$, keep `KappaPolynomial` / `BoundedExponential` as the automatic default whenever an iterative solver is chosen and reserve `Rational` strictly for direct factorization paths.
4. **Existing Golden Drift:** If switching `PmlModeRequest::default()` to `Auto` causes any existing numerical golden in `nyx-{ems,sms}` to drift beyond tolerance, preserve the golden test by pinning its historical fixture JSON explicitly to `"pml_mode": "inline"` rather than running `NYX_WRITE_GOLDEN=1`.
