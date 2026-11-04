# 53. A Maxwell residual estimator for edge elements, and an adaptation loop that takes any estimator

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). Item E1 of `docs/ports/meshers/FEATURE_DEPENDENCIES.md`, section 2 of `MISSING_FEATURES.md` ("Lowest-order edge elements"). Uses ADR 51 for the adaptation of geom models.

## Context

The solve-estimate-adapt loop of `adapt::estimate3` reads plotted nodal fields and recovers a Hessian. For edge (Nedelec) elements it works at order 2 (2.6e-6 on 494 tetrahedra against 1.8e-4 on 696 uniform ones for a cavity mode) and does not at order 1 (1.2e-2 on 1490 against 5.3e-3 on 696): the field is piecewise linear and only tangentially continuous, so its Hessian is noise. The loop also fixed both the estimator (fields by name) and the adaptation (`adapt3` with the facets of the start mesh), so a geom model could not use the exact faces (ADR 51).

## Decision

* **`adapt::maxwell`: the residual estimator of Beck, Hiptmair, Hoppe and Wohlmuth** for `curl (1/mu) curl E - k0^2 eps E = 0`. Per tetrahedron: `h_T^2 ||curl (1/mu) curl E - k0^2 eps E||^2 + h_T^2 ||div (eps E)||^2` plus half the face terms `h_F ||n x [(1/mu) curl E]||^2 + h_F ||[eps E . n]||^2` of its faces (the tangential trace of `E` is continuous). A boundary face is a perfect conductor (nothing) unless its label is declared natural (a magnetic wall or symmetry plane: `h_F ||n x (1/mu) curl E||^2`). `maxwell_indicators` returns `eta_T^2` and the four parts summed.
* **The field comes through one trait**, `CurlField::sample(cell, lam) -> (E, curl E)` at a barycentric point, complex, polynomial of degree at most 2 in the cell (orders 1 and 2) and continued a little outside the cell, because first derivatives are central differences along the barycentric directions (exact for quadratics). The estimator does not know how degrees of freedom are numbered; the physics crates implement the trait for theirs (the test adapter in `nyx-ems` uses `cavity::basis3` and `Ned3`). Nothing in a physics crate is touched by this ADR (`sms` and `ems` are being moved in parallel).
* **`indicator_metrics`: equidistribution.** A cell of size `s` (the edge of the regular tetrahedron of its volume) with indicator `eta` asks for `s (tau / eta)^(1/(order + 3/2))`, kept within a factor 4 of `s` and between `h_min_fraction * h_max` and `h_max`; `tau` is found by bisection so that the predicted cell count `sum (s / s_new)^3` equals the target. A vertex takes the smallest size of its cells; the result is an isotropic metric, intersected with an optional floor (`adapt3d::geometry_floor3` for geom models).
* **`estimate3::adapt_loop3`**: the loop of `solve_adaptively3` with the two choices passed in: what the solution asks of the mesh (`metrics`) and how a mesh is adapted to it (`adapt`). `solve_adaptively3` is now a call of it with the Hessian estimator and the plain `adapt3`; its behaviour is unchanged (cell counts of the existing tests are identical). For a geom model `adapt` is `adapt3d::adapt_to_metrics3`.
* Limits: real isotropic coefficients, no source term, absorbing layers are left out (`skip`: their coefficients are complex tensors), eigenmodes and homogeneous problems. Goal-oriented (dual-weighted) adaptivity is E2 and not done.

## Consequences

Measured on the PEC box cavity 40 x 30 x 20 mm, lowest mode (`crates/nyx-ems/tests/maxwell_estimator.rs`, closed form):

| order | uniform blocks | Hessian loop (earlier) | Maxwell loop |
|---|---|---|---|
| 1 | 5.3e-3 on 696 tets, 6.9e-3 on 1470 | 1.2e-2 on 1490 | **5.2e-4 on 333**, 5.6e-5 on 882 |
| 2 | 1.8e-4 on 696, 6.7e-5 on 1470 | 2.6e-6 on 494 | 1.7e-4 on 345, **8.0e-6 on 1402** |

The total indicator falls by 2.19 per halving of `h` at order 1 and 4.9 at order 2 (rates 1 and 2). On the sphere cavity from geom (loop with `adapt_to_metrics3`, exact sphere as the oracle, curvature floor), the error falls from 1.8e-2 to about 5e-3 in two cycles, where the polyhedral sphere (ADR 18) is the limit for both orders.

## Code

`crates/femcore/src/adapt/{maxwell,estimate3}.rs`.

## Contract

`adapt_maxwell_estimator_contract` (`tests/contracts.rs`: the trait, the four parts adding up, natural boundaries, skipped domains, refusals); the closed-form residual of a field with known curl and the equidistribution in the unit tests of `adapt/maxwell.rs`; `crates/nyx-ems/tests/maxwell_estimator.rs` (needs the cavity kernel: rates, loop against uniform blocks at both orders, the chain on a sphere from geom).
