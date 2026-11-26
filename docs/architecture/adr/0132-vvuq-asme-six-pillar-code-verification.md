# 132. VVUQ Phase 2: ASME 6-Pillar numerical code verification and Rough-Spot closures

Date: 2026-11-26

## Status

Accepted (VVUQ Phase 2 / Layer 2; renumbered from provisional `0126` in `docs/vvuq/VVUQ_VERIFICATION_PLAN_EGM_REFINED.md` because ADRs `0125–0130` were already allocated and Phase 0–1 took `0131`). Extends ADR 10, ADR 19, ADR 41, ADR 42, ADR 45, ADR 82, and ADR 131.

## Context

Layer 2 ("solving the equations right") requires verifying every PDE solver in `nyx-dfs`, `nyx-sms`, and `nyx-ems` against: (a) exact analytical closed forms, (b) true PDE-level Method of Manufactured Solutions (MMS) with non-polynomial solutions, (c) formal observed $h$-, $p$-, and $\Delta t$-convergence orders, (d) conservation laws & operator symmetries, (e) Euclidean/permutation/port-swap/scaling/spectral invariances plus geometry & curved-element correctness, and (f) golden governance (splitting frozen-mesh solver goldens from live-mesher physics checks per R22 and closing R6a/b/d, R29, R31, and R32–R34).

## Decision

* **Zero-overhead MMS & shared rate checker (`femcore::vvuq`):** Provide `SourceField<const D: usize>` and zero-sized `UniformSource` in `femcore::vvuq` so volumetric assemblers (`nyx-dfs::{electro,heat,magneto}`, `nyx-sms::solid`, `nyx-ems::cavity`) share 100% of production quadrature and Piola transforms with zero production overhead, alongside `assert_convergence_rate` for formal OLS and asymptotic $p_{\text{obs}}$ checks.
* **Geometry & meshing correctness closure (`geom`, `femcore`):** Implement analytic quadric flux integrals in `Complex::check3` + 2-manifold shell edge pairing in `topo3`, union-find canonical clustering in `Verts::intern`, exact quadratic curve ray-crossing in `edge_crosses_ray_x`, safe BVH traversal stack spilling, segmented Bernstein positivity for non-QP `CurvedTet` vs interior `q8` validation for quarter-point cells (with periodic twin sync), slanted 2D periodic wall splitting (`R29`), and `ShiftedFactorization` pivot/inertia diagnostics + fallback (`R31`).
* **6-Pillar verification & Rough-Spot closures (`nyx-dfs`, `nyx-sms`, `nyx-ems`):** Cover every study family across Pillars (a)–(f), close `R6a` (`u_r = i m u_theta` at `r = 0`), `R6b` (`r`-weighted `lrb::damping_matrix`), `R6d` (quadrature-point `s_theta`), `R32–R34` curved mesh guards and curved `ems::cavity` post-solve integrals, `R22` canonical `.fmsh` solver-only vs live-mesher split, and deterministic `emw_pml_x` `+4` Krylov / `layer_fraction < 0.25` selection.

## Consequences

* Pure-elastic / non-MMS code paths remain bit-identical via zero-sized `UniformSource`.
* `crates/nyx-ems/tests/curved.rs:251` and `crates/nyx-ems/tests/meshing.rs:80-97` are un-ignored with formal quantitative assertions.

## Code

`crates/femcore/src/vvuq.rs`, `crates/geom/src/{prim,complex,topo3,arrange,bvh}.rs`, `crates/geom/tests/{arrangement,bvh,stack3d}.rs`, `crates/femcore/src/meshing/{curve3,geomdrive,geomdrive3,mesher,request,sheet3d,sheet3d_curve,singularity}.rs`, `crates/femcore/src/linalg/direct.rs`, `crates/femcore/tests/{curved,meshing,sheet3d_curved_arcs,singularity_qp}.rs`, `crates/nyx-dfs/src/{common,electro,heat,lib,magneto}.rs`, `crates/nyx-dfs/tests/{qcap_solvers,heat_bdf2_order,vvuq_closed_forms,vvuq_mms_convergence,vvuq_conservation_invariance,vvuq_golden}.rs`, `crates/nyx-sms/src/{axisym,cross_section,floquet,lib,lrb,piezo_asm,pml,problem,solid,waveguide}.rs`, `crates/nyx-sms/tests/{vvuq_closed_forms,vvuq_mms_convergence,vvuq_conservation_invariance,vvuq_golden}.rs`, `crates/nyx-ems/src/{axisym,cavity,ibc,mixed_ops,modes,modes_complex,optical,ports}.rs`, `crates/nyx-ems/tests/{axisym_pml,closed_forms,conservation_invariance,convergence,curved,golden,meshing,mms_maxwell}.rs`, `crates/nyx-ems/tests/golden/{emw_soi_wire.fmsh,emw_soi_hybrid.fmsh,emw_soi_wire_solver.v1.json,emw_soi_hybrid_solver.v1.json,emw_cavity_3d_box.v1.json,emw_bloch_3d_wire.v1.json,emw_driven_3d_prom.v1.json}`, `crates/nyx-engine/src/xdmf.rs`.
