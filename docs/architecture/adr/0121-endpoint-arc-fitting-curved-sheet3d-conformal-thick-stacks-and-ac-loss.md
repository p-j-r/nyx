# 121. 1-DOF endpoint-anchored GDSII arc fitting, curved `sheet3d` projection, conformal thick-metal / trench stacks, and quasi-electrostatic AC loss

Date: 2026-11-22

## Status

Accepted. Amends ADR 2, 3, 82, 102, 103, 105, 109, 114, 115, 116, and 118.

## Context

Two geometric and constitutive approximations remain above the $0.1\%$ capacitance floor on realistic superconducting and RF devices: (1) GDSII layouts approximate circular pads, concentric rings, and rounded corners as $N$-gon polylines whose chordal area deficit $\Delta A/A = 1 - \frac{N}{2\pi}\sin(2\pi/N)$ biases capacitance unless arcs are recovered and mapped onto isoparametric `CurvedTet` edges; and (2) real chips have finite-thickness metal films ($t_m \sim 100\text{–}300\,\mathrm{nm}$), substrate over-etch trenches ($t_e \sim 50\text{–}200\,\mathrm{nm}$), and microwave dielectric loss ($\tan\delta$, $\sigma$) at $f \sim 5\text{–}15\,\mathrm{GHz}$ where evaluating `[optical]` dispersion fails out of range and zero-thickness TLS edge asymptotics would double-count meshed sidewalls.

## Decision

* **1-DOF endpoint-anchored circular arc fitting & curved `sheet3d` projection (`AF-5`, `geom::{poly_batch,gds}`, `femcore::meshing::sheet3d_curve`, `nyx-dfs::electro`):**
  * `Op::Polygons` (`add_polygons`) and `import_gds` accept `arc_tol >= 0` (default `0.0`): maximal runs of $\ge 4$ vertices between pinned/sharp corners are fitted via the exact perpendicular-bisector 1-DOF family $\mathbf{c}(s) = \mathbf{m} + s\,\mathbf{n}$ anchoring exact endpoints $\mathbf{p}_a, \mathbf{p}_b$ (`0 ULP` gap), guarded against self-intersection and gap encroachment, and emitted directly as `Curve::Circle` edges in `planar::from_scene`.
  * When `set_mesh(curved=true)` is active, `sheet3d_curve::project_sheet3d_curved_edges` projects horizontal `EdgeGeom::Horizontal` mid-edge nodes onto exact `Curve::Circle` / `Curve::Ellipse` arcs and relaxes interior incident edges via shared `curve3::relax_edge_nodes`. When combined with `singular_edges = "quarter_point"` (Rule 4), arc projection runs first and Barsoum shifting runs second with $\xi_{\text{QP}}$ clamped upward where needed; `nyx-dfs::electro` integrates curved 3D Laplace tensors with 216-point `quad_rule_collapsed(8)` and validates $\det J \cdot \det J_{\text{straight}} > 0$ at every Gauss point.
* **Conformal finite-thickness metal & substrate over-etch trenches (`AF-6a`, `geom::{model,stack,labels}`, `femcore::meshing::sheet3d`, `nyx-analysis::tls`):**
  * `Op::Sheet` (`add_sheet`) gains `thickness >= 0` (`[z_s, z_s + t_m]` conductor cavity slab, tagged `<name>`, `<name>:top`, `<name>:bottom`, `<name>:side`) and `etch_depth >= 0` (`[z_s - t_e, z_s]` trench slab reassigned to top dielectric), extruded as single-slab vertical intervals (`vec![l.z0, l.z1]`) in `sheet3d` while preserving exact `0 nm` bit-identity when both are zero.
  * Rule 2 automatically skips zero-thickness `quarter_point` shifting on sheets with `thickness > 0`. `LabelMap` / `LabelResolver` / `TlsMeta::sheet_thickness_m` propagate finite sheet thicknesses into `nyx-analysis::tls`, integrating meshed top/bottom/sidewall `MA`/`MS`/`SA` facets directly without zero-thickness knife-edge asymptotic double-counting.
* **Quasi-electrostatic AC dielectric loss & compatibility rules (`AF-6b`, `nyx-dfs::{ac_loss,planar_bem}`, `nyx-engine::materials`):**
  * `dfs` `electrostatic` gains `freq` (Hz, default `0.0`) and tensor properties `tan_delta` and `sigma`; `study_wavelength` returns `None` for `dfs` `electrostatic` while `study_frequency` populates `Ctx::frequency`, evaluating `[dielectric]` / `[electrical]` loss at microwave frequencies without triggering `[optical]` dispersion errors.
  * `nyx-dfs::ac_loss` (and complex-image `planar_bem`) computes `conductance_matrix_s` ($G_{ij}(\omega)$), `q_dielectric_matrix` ($Q_{ij} = \omega |C_{ij}|/|G_{ij}|$), and the strictly positive second-order / complex Schur capacitance shift $\Delta C^{(2)} = +(K_{\text{loss}}\mathbf{v}_i)_f^T K_{ff}^{-1}(K_{\text{loss}}\mathbf{v}_j)_f > 0$ when $\max \tan\delta_{\text{eff}} > 10^{-2}$; `planar_bem` rejects `thickness > 0`, `etch_depth > 0`, `flux_bounds != "none"`, `far_field`, and `adapt_passes > 1` with `Error::Unsupported` (Rule 5).

## Consequences

* Curved GDSII pads and annular rings converge to exact circular geometry without $N$-gon chordal bias, finite-thickness etched superconducting stacks yield physical 3D `MA`/`MS`/`SA` participation and microwave dielectric $Q$ matrices at $13\,\mathrm{GHz}$, and all cross-feature combinations either compose deterministically or fail fast with typed `Invalid` / `Unsupported` diagnostics.

## Code

`crates/geom/src/{model,poly_batch,planar,scene,stack,topo3,labels,gds,print}.rs`, `crates/geom/tests/{arc_fitting,thick_sheet_stack}.rs`, `crates/femcore/src/meshing/{sheet3d_curve,sheet3d,singularity,curve3}.rs`, `crates/femcore/tests/sheet3d_curved_arcs.rs`, `crates/nyx-physics/src/{cqed,physics,decode}.rs`, `crates/nyx-dfs/src/{ac_loss,electro,planar_bem,problem,lib}.rs`, `crates/nyx-dfs/tests/ac_loss_electro.rs`, `crates/nyx-analysis/src/{tls,capacitance}.rs`, `crates/nyx-analysis/tests/tls_thick_sheet.rs`, `crates/nyx-engine/src/{geomdesc,gds,materials,commands,meshing}.rs`, `crates/app/src/geomui.rs`, `crates/nyx-engine/tests/accuracy_features_e2e.rs`, `crates/nyx-script/tests/accuracy_features_script.rs`.
