# Finite-element basics: code

C4 level 4 (code) of the component **Finite-element basics** (`femcore.fem`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> Lagrange P1/P2 elements and quadrature, the 2D edge-element basis, derived fields, PML coordinate stretching, Floquet pairs and material properties per domain.

**Decisions:** [ADR 9](../../../adr/0009-kernels-are-three-stages.md), [ADR 4](../../../adr/0004-simplexmesh-is-the-solver-contract.md), [ADR 45](../../../adr/0045-ems-moves-into-nyx-ems.md), [ADR 69](../../../adr/0069-axis-conditions-and-cylindrical-stretch-of-the-axisymmetric-solid.md)  
**Size:** 6 files, 2025 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/edge.rs` | 81 |  |
| `crates/femcore/src/fem.rs` | 1128 |  |
| `crates/femcore/src/fields.rs` | 220 |  |
| `crates/femcore/src/periodic.rs` | 290 |  |
| `crates/femcore/src/pml.rs` | 98 |  |
| `crates/femcore/src/props.rs` | 208 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Matrices and vectors](./lacore.md) | linalg::csr, linalg::csr::assembly | C, CCoo, CCsr, Coo, RCsr, c |
| [Linear algebra driver](./ladriver.md) | linalg::driver | DirectSolver |
| [Preconditioners](./laprecond.md) | linalg::amg_elastic | Amg, AmgOptions, pcg |
| [Material specification](./material.md) | material | Anisotropy, M3, M6 |
| [Mesh](./mesh.md) | mesh | SimplexMesh |
| [Runtime services](./runtime.md) | par, progress | check, lap_since, stage, threads, threads_up_to |
| [Frame and sections](./section.md) | section | mesh_pml_dirs |

## Used by

| Component | Modules | Items |
|---|---|---|
| [analysis: Field samples](../analysis/samples.md) | fem, fields | geometry, local_bary, shape |
| [dfs: Shared helpers](../dfs/common.md) | fields | gradients_at_nodes |
| [dfs: Electrostatics kernel](../dfs/electro.md) | fem, fields, props | CellFields, TensorKind, node_average |
| [dfs: Heat kernel](../dfs/heat.md) | fem, fields, props | AmgSystem, CellFields, LinearSolver, TensorKind, amg_wanted, node_average |
| [dfs: Typed problem](../dfs/problem.md) | props | LinearSolver, Materials |
| [ems: Physics adapter](../ems/adapter.md) | props | PmlDir |
| [ems: Axisymmetric kernel](../ems/axisym.md) | edge, fem, fields, pml, props | CellData, CellFields, Materials, Pml, PmlDir, cell_data, edge_basis, local_bary +4 |
| [ems: 3D cavity kernel](../ems/cavity.md) | fem, fields, periodic, pml, props | CellFields, Match, Reduction, TensorKind, ccsr_adjoint, ccsr_mul, cell_stretches, congruence +5 |
| [ems: Mode kernel](../ems/modes.md) | edge, fem, fields, props | CellFields, TensorKind, cell_data, edge_basis, local_bary, node_average, tri_g |
| [ems: Complex-mode kernel](../ems/modescx.md) | edge, fem, fields | CellFields, cell_data, edge_basis, local_bary, node_average, tri_g |
| [ems: Optical mode sets](../ems/optical.md) | props | PmlDir, TensorKind |
| [ems: Typed problem](../ems/problem.md) | props | Materials, Periodic, Pml, PmlDir |
| [engine: Commands](../engine/commands.md) | fem, props | RefElement, TensorKind, facet_rule, facets_with, geometry, quad_rule, shape |
| [engine: Registry and facade](../engine/registry.md) | fields | Display |
| [Error estimators](./estimators.md) | fem | quad_rule |
| [Meshing facade](./meshing.md) | props | Pml, PmlDir |
| [3D meshing](./meshing3d.md) | fem | CurvedTet, quad_rule_collapsed |
| [Optical mode sets](./optical.md) | edge, fem, pml, props | CellData, cell_data, edge_basis, stretch, tri_g |
| [Frame and sections](./section.md) | props |  |
| [Structured absorbing layers](./structpml.md) | props | Materials, Pml, PmlDir |
| [physics: Problem envelope](../physics/envelope.md) | props | Materials |
| [physics: Result and frame](../physics/result.md) | fields | CellField, CellFields |
| [physics: Test support](../physics/testsupport.md) | fields | node_average |
| [sms: Axisymmetric kernel](../sms/axisym.md) | fem, pml, props | Pml, PmlDir, cell_stretches, radial_coordinate |
| [sms: Cross-section study kernel](../sms/crosssection.md) | fem, fields, periodic, pml, props | CellFields, PmlDir, Reduction, cell_stretches, check_boundary_labels, congruence, match_vertices, nodal_relations |
| [sms: Floquet kernel](../sms/floquet.md) | periodic | Reduction, congruence, match_vertices, nodal_relations |
| [sms: Low-reflecting boundary](../sms/lrb.md) | fem |  |
| [sms: Absorbing-layer kernel](../sms/pml.md) | fem, periodic, pml | Reduction, cell_stretches, congruence, match_vertices, nodal_relations |
| [sms: Typed problem](../sms/problem.md) | props | LinearSolver, Materials, Periodic, Pml, PmlDir |
| [sms: Solid kernel](../sms/solid.md) | fem, fields, props | AMG_THRESHOLD, AmgSystem, CellFields, FACTORISE_MAX_FREE, LOBPCG_MAX_FACTOR, LOBPCG_MIN_FREE, LinearSolver, Materials +5 |
| [sms: Waveguide kernel](../sms/waveguide.md) | fem |  |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace edge_rs {
    class CellData {
        <<struct>>
        +sign Vec
        +gdofs Vec
        +net usize
    }
    class fns_edge_rs {
        <<module edge.rs>>
        +edge_basis(l, g, sign, order, w, curl)
        +cell_data(mesh, order) CellData
        +tri_g(geo, c) slice
    }
  }
  namespace fem_rs {
    class Quad {
        <<struct>>
        +lam Vec
        +w Vec~f64~
    }
    class RefElement {
        <<struct>>
        +dim usize
        +order usize
        +nloc usize
        +w Vec~f64~
        +n Vec
        +dn Vec
        +new(dim, order, degree) Self
    }
    class Geometry {
        <<struct>>
        +vol Vec~f64~
        +g Vec~f64~
    }
    class CurvedTet {
        <<struct>>
        +x slice
        +new(x) Self
        +straight(v) Self
        +position(lam) slice
        +jacobian_matrix(lam) slice
        +jacobian(lam) tuple
        +min_det_at_nodes() f64
        +volume() f64
    }
    class AmgSystem {
        <<struct>>
        +setup_seconds f64
        +new(mesh, order, ncomp, a, fixed) Result
        +n_levels() usize
        +sizes() Vec~usize~
        +solve_with(matrix, b, fixed_values, tol, max_iter) tuple
        +convergence_factor() f64
        +precondition_free(r_free) Vec~f64~
        +solve_free(matrix, b_free, tol) Result
    }
    class fns_fem_rs {
        <<module fem.rs>>
        +quad_rule(dim, degree) Quad
        +facet_rule(facet_dim) Quad
        +shape(lam, order) tuple
        +geometry(mesh) Geometry
        +quad_rule_collapsed(degree) Quad
        +shape_gradients(re, geo, c, q, out)
        +assemble_parallel(nc, nrows, ncols, f) Csr
        +assemble_laplace(mesh, order, coeff) Csr
        +assemble_laplace_tensor(mesh, order, k) Csr
        +assemble_mass(mesh, order, coeff) Csr
        +more_15_functions()
        +TET_EDGES const
        +AMG_THRESHOLD const
        +AMG_THRESHOLD_3D const
        +LOBPCG_MIN_FREE const
        +more_2_items
    }
  }
  namespace fields_rs {
    class Display {
        <<struct>>
        +dim usize
        +dofs Vec~usize~
        +points Vec
        +tris Vec
        +tri_label Vec~i32~
        +new(mesh, order) Self
        +pick(values) Vec~f32~
    }
    class CellField {
        <<struct>>
        +name String
        +unit String
        +re Vec~f64~
        +im Vec~f64~
    }
    class CellFields {
        <<struct>>
        +order usize
        +nloc usize
        +fields Vec~CellField~
        +new(mesh, order) Self
        +get(name) Option
        +push_real(name, unit, re)
        +push_complex(name, unit, re, im)
        +push_components(prefix, names, unit, g)
        +node_average(mesh, name) Option
    }
    class fns_fields_rs {
        <<module fields.rs>>
        +local_bary(dim, order) Vec
        +gradients_at_nodes(mesh, u, ncomp, order) Vec~f64~
        +node_average(mesh, vals, m, order) Vec~f64~
    }
  }
  namespace periodic_rs {
    class Match {
        <<struct>>
        +map HashMap
        +src_facets Vec~usize~
        +shift Vec~f64~
    }
    class Reduction {
        <<struct>>
        +n usize
        +roots Vec~usize~
        +index_of Vec~usize~
        +ratio Vec~C~
        +build(n, relations, fixed) Reduction
        +n_reduced() usize
        +prolong() CCsr
        +expand(xr) Vec~C~
    }
    class fns_periodic_rs {
        <<module periodic.rs>>
        +phase_factor(phase) C
        +match_vertices(mesh, p) Result
        +nodal_relations(mesh, order, ncomp, m, phase) Result
        +ccsr_mul(a, b) CCsr
        +ccsr_adjoint(a) CCsr
        +congruence(p, a) CCsr
    }
  }
  namespace pml_rs {
    class fns_pml_rs {
        <<module pml.rs>>
        +stretch(p, dirs, centre, lo, hi) tuple
        +radial_coordinate(p, r, lo, hi) C
        +cell_stretches(mesh, pml) VecC
    }
  }
  namespace props_rs {
    class Materials {
        <<struct>>
        +domains BTreeMap
        +default BTreeMap
        +aniso BTreeMap
        +pml BTreeMap
        +has_pml() bool
        +get(domain, name) Option~f64~
        +set(domain, name, v)
        +has_tensor(kind) bool
        +has_stiffness() bool
        +cell_stiffness(mesh) Result
        +cell_stiffness_skipping(mesh, skip) Result
        +cell_tensor(mesh, kind, fallback) Result
        +more_1_methods()
    }
    class PmlDir {
        <<enumeration>>
        Off
        Plus
        Minus
    }
    class Pml {
        <<struct>>
        +x PmlDir
        +y PmlDir
        +z PmlDir
        +strength f64
        +order f64
        +re f64
    }
    class TensorKind {
        <<enumeration>>
        Conductivity
        Permittivity
        Expansion
    }
    class Periodic {
        <<struct>>
        +source Vec~i32~
        +dest Vec~i32~
        +phase f64
    }
    class LinearSolver {
        <<enumeration>>
        Auto
        Direct
        Iterative
        +name() str
    }
  }
  CellFields --> CellField
  Materials --> Pml
  Pml --> PmlDir
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `edge_basis` | `edge.rs:21` | `pub fn edge_basis(l: &[f64], g: &[[f64; 2]; 3], sign: &[f64; 3], order: usize, w: &mut [[f64; 2]], curl: &mut [f64])` —  Nedelec basis at barycentric point `l` of one cell with barycentric gradients `g[i]`: |
| struct | `CellData` | `edge.rs:48` | fields: `sign`, `gdofs`, `net` |
| fn | `cell_data` | `edge.rs:54` | `pub fn cell_data(mesh: &SimplexMesh, order: usize) -> CellData` |
| fn | `tri_g` | `edge.rs:76` | `pub fn tri_g(geo: &Geometry, c: usize) -> [[f64` |
| struct | `Quad` | `fem.rs:26` | fields: `lam`, `w` —  Quadrature on the reference simplex: barycentric points and weights summing to 1 |
| fn | `quad_rule` | `fem.rs:37` | `pub fn quad_rule(dim: usize, degree: usize) -> Quad` |
| fn | `facet_rule` | `fem.rs:82` | `pub fn facet_rule(facet_dim: usize) -> Quad` —  Rule on a boundary facet (Gauss on edges in 2D, degree-4 Dunavant on triangles in 3D). |
| fn | `shape` | `fem.rs:100` | `pub fn shape(lam: &[f64], order: usize) -> (Vec<f64>, Vec<f64>)` —  Shape functions at barycentric point `lam`: values `n[a]` and derivatives with respect to |
| struct | `RefElement` | `fem.rs:129` | fields: `dim`, `order`, `nloc`, `w`, `n`, `dn` —  Quadrature rule with shape functions tabulated at its points. |
| struct | `Geometry` | `fem.rs:157` | fields: `vol`, `g` —  Per-cell volume and gradients of the barycentric coordinates: `g[c*(d+1)*d + b*d + k]`. |
| fn | `geometry` | `fem.rs:162` | `pub fn geometry(mesh: &SimplexMesh) -> Geometry` |
| fn | `quad_rule_collapsed` | `fem.rs:257` | `pub fn quad_rule_collapsed(degree: usize) -> Quad` —  A quadrature rule on the reference tetrahedron that is exact to at least `degree`, built by collapsing the unit cube (Duffy |
| const | `TET_EDGES` | `fem.rs:275` | `[(usize, usize)` —  The local vertex pairs of a tetrahedron in the order of its edge numbering (and of the geometry nodes 4 to 9 of [`CurvedTet`]). |
| struct | `CurvedTet` | `fem.rs:285` | fields: `x` —  A tetrahedron with quadratic (isoparametric) geometry: the map from the reference cell is the quadratic Lagrange interpolation of ten |
| fn | `shape_gradients` | `fem.rs:397` | `pub fn shape_gradients(re: &RefElement, geo: &Geometry, c: usize, q: usize, out: &mut [f64])` —  Physical gradients of the shape functions at quadrature point `q` of cell `c`: |
| fn | `assemble_parallel` | `fem.rs:415` | `pub fn assemble_parallel<F>(nc: usize, nrows: usize, ncols: usize, f: F) -> Csr` —  Runs `f(range, &mut coo)` over cell chunks on all available threads and merges the results. |
| fn | `assemble_laplace` | `fem.rs:461` | `pub fn assemble_laplace(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Csr` —  `K_ij = integral coeff grad(phi_i) . grad(phi_j)`; `coeff` holds one value per cell. |
| fn | `assemble_laplace_tensor` | `fem.rs:495` | `pub fn assemble_laplace_tensor(mesh: &SimplexMesh, order: usize, k: &[[[f64; 3]; 3]]) -> Csr` —  `K_ij = integral grad(phi_i) . K grad(phi_j)` with a symmetric tensor `K` per cell (only |
| fn | `assemble_mass` | `fem.rs:535` | `pub fn assemble_mass(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Csr` —  `M_ij = integral coeff phi_i phi_j`. |
| fn | `assemble_source` | `fem.rs:564` | `pub fn assemble_source(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Vec<f64>` —  Integral of `coeff phi_i` over every cell: the load of a per-cell source density. |
| fn | `facet_measure` | `fem.rs:571` | `pub fn facet_measure(mesh: &SimplexMesh, f: usize) -> f64` |
| fn | `facets_with` | `fem.rs:585` | `pub fn facets_with<'a>(mesh: &'a SimplexMesh, labels: &'a [i32]) -> impl Iterator<Item = usize> + 'a` |
| fn | `assemble_boundary_mass` | `fem.rs:590` | `pub fn assemble_boundary_mass(mesh: &SimplexMesh, labels: &[i32], order: usize, coeff: f64) -> Csr` —  `integral coeff phi_i phi_j` over the facets carrying one of `labels` (convective term). |
| fn | `assemble_boundary_load` | `fem.rs:616` | `pub fn assemble_boundary_load(mesh: &SimplexMesh, labels: &[i32], value: f64, order: usize) -> Vec<f64>` —  Load `integral value phi_i` over the facets carrying one of `labels`. |
| fn | `add_boundary_load` | `fem.rs:624` | `pub fn add_boundary_load(mesh: &SimplexMesh, labels: &[i32], value: &dyn Fn(usize) -> f64, order: usize, out: &mut [f64], ncomp: usize, comp: usize)` —  Adds the boundary load of a per-facet value into component `comp` of a vector of `ncomp` |
| fn | `boundary_dofs` | `fem.rs:647` | `pub fn boundary_dofs(mesh: &SimplexMesh, labels: &[i32], order: usize) -> Vec<usize>` —  Dofs on the boundary facets with any of the given boundary numbers (sorted, unique). |
| fn | `dirichlet_constraints` | `fem.rs:660` | `pub fn dirichlet_constraints(mesh: &SimplexMesh, values: &BTreeMap<i32, f64>, order: usize) -> (Vec<usize>, Vec<f64>)` —  Fixed dofs and values for `{boundary number: value}`. A dof shared by boundaries with |
| fn | `solve_dirichlet` | `fem.rs:679` | `pub fn solve_dirichlet(k: &Csr, f: &[f64], fixed: &[usize], values: &[f64]) -> Result<Vec<f64>, String>` —  Solves `K u = f` with `u[fixed] = values` by elimination. |
| fn | `untouched_dofs` | `fem.rs:706` | `pub fn untouched_dofs(mesh: &SimplexMesh, order: usize, ncomp: usize, inactive: &[i32]) -> Vec<usize>` —  Dofs (of `ncomp` component blocks of `n_dofs(order)`) that no active cell touches: the cells |
| fn | `free_dofs` | `fem.rs:729` | `pub fn free_dofs(n: usize, fixed: &[usize]) -> Vec<usize>` —  Free dof indices (complement of `fixed`). |
| const | `AMG_THRESHOLD` | `fem.rs:743` | `usize` —  Free-dof count above which `LinearSolver::Auto` switches to multigrid in 2D. In 3D multigrid |
| const | `AMG_THRESHOLD_3D` | `fem.rs:744` | `usize` |
| const | `LOBPCG_MIN_FREE` | `fem.rs:748` | `usize` —  Free dofs above which eigenfrequency studies (Auto) try LOBPCG instead of a factorisation. Measured with the shared |
| const | `FACTORISE_MAX_FREE` | `fem.rs:750` | `usize` —  Above this many free dofs a factorisation is not attempted for eigenfrequencies. |
| const | `LOBPCG_MAX_FACTOR` | `fem.rs:752` | `f64` —  Multigrid is used for eigenfrequencies when its convergence factor is at most this. |
| fn | `amg_wanted` | `fem.rs:755` | `pub fn amg_wanted(dim: usize, n_free: usize) -> bool` —  Whether `LinearSolver::Auto` should use multigrid for a system with `n_free` unknowns. |
| struct | `AmgSystem` | `fem.rs:820` | fields: `setup_seconds` —  A symmetric positive definite system over all dofs with some dofs constrained, prepared for |
| fn | `solve_spd` | `fem.rs:1006` | `pub fn solve_spd(mesh: &SimplexMesh, order: usize, ncomp: usize, solver: crate::props::LinearSolver, k: &Csr, f: &[f64], fixed: &[usize], values: &[f64]) -> Result<(Vec<f64>, Option<String>), String>` —  Solves `K u = f` with `u[fixed] = values` choosing between the sparse direct solver and |
| fn | `check_boundary_labels` | `fem.rs:1048` | `pub fn check_boundary_labels(mesh: &SimplexMesh, labels: impl IntoIterator<Item = i32>) -> Result<(), String>` —  Refuses boundary numbers that the mesh does not have, naming the ones it has. Every physics checks its conditions with it. |
| fn | `cell_values_at_nodes` | `fem.rs:1061` | `pub fn cell_values_at_nodes(mesh: &SimplexMesh, u: &[f64], order: usize) -> Vec<f64>` —  The values of a dof field at the local dofs of every cell: `out[c*nloc + a]`. |
| fn | `RefElement::new` | `fem.rs:141` | `pub fn new(dim: usize, order: usize, degree: usize) -> Self` |
| fn | `CurvedTet::new` | `fem.rs:290` | `pub fn new(x: [[f64; 3]; 10]) -> Self` |
| fn | `CurvedTet::straight` | `fem.rs:295` | `pub fn straight(v: [[f64; 3]; 4]) -> Self` —  The straight cell with vertices `v` (edge nodes at the midpoints). |
| fn | `CurvedTet::position` | `fem.rs:328` | `pub fn position(&self, lam: &[f64; 4]) -> [f64` —  Position at the barycentric point `lam`. |
| fn | `CurvedTet::jacobian_matrix` | `fem.rs:346` | `pub fn jacobian_matrix(&self, lam: &[f64; 4]) -> [[f64` —  The Jacobian matrix `J[i][k] = dx_i / dxi_k` with respect to the reference coordinates `(lambda_1, lambda_2, lambda_3)`. |
| fn | `CurvedTet::jacobian` | `fem.rs:359` | `pub fn jacobian(&self, lam: &[f64; 4]) -> ([[f64` —  Physical gradients of the four barycentric functions at `lam`, and the determinant of the Jacobian (negative for an inverted |
| fn | `CurvedTet::min_det_at_nodes` | `fem.rs:376` | `pub fn min_det_at_nodes(&self) -> f64` —  The smallest determinant over the ten nodes of the cell (vertices and edge nodes), a cheap validity probe. |
| fn | `CurvedTet::volume` | `fem.rs:388` | `pub fn volume(&self) -> f64` —  Volume of the curved cell: the integral of `\|det J\|` over the reference tetrahedron, with a rule exact to degree 8. |
| fn | `AmgSystem::new` | `fem.rs:835` | `pub fn new(mesh: &SimplexMesh, order: usize, ncomp: usize, a: &Csr, fixed: &[usize]) -> Result<AmgSystem, String>` —  `a` is the matrix over all `n = ncomp * n_dofs(order)` dofs (entries in constrained rows and |
| fn | `AmgSystem::n_levels` | `fem.rs:912` | `pub fn n_levels(&self) -> usize` |
| fn | `AmgSystem::sizes` | `fem.rs:916` | `pub fn sizes(&self) -> Vec<usize>` |
| fn | `AmgSystem::solve_with` | `fem.rs:923` | `pub fn solve_with(&self, matrix: Option<&Csr>, b: &[f64], fixed_values: &[(usize, f64)], tol: f64, max_iter: usize) -> (Vec<f64>, usize, f64)` —  Solves `A x = b` for the right-hand side over all dofs (entries at constrained dofs are taken |
| fn | `AmgSystem::convergence_factor` | `fem.rs:953` | `pub fn convergence_factor(&self) -> f64` —  Average residual reduction per iteration of multigrid-preconditioned CG on a random |
| fn | `AmgSystem::precondition_free` | `fem.rs:971` | `pub fn precondition_free(&self, r_free: &[f64]) -> Vec<f64>` —  One multigrid cycle applied to a vector over the free dofs: an approximation of the inverse |
| fn | `AmgSystem::solve_free` | `fem.rs:986` | `pub fn solve_free(&self, matrix: Option<&Csr>, b_free: &[f64], tol: f64) -> Result<Vec<f64>, String>` —  Solves for the unconstrained dofs only: `b` and the result are vectors over the free dofs. |
| fn | `local_bary` | `fields.rs:15` | `pub fn local_bary(dim: usize, order: usize) -> Vec<Vec<f64>>` —  Barycentric coordinates of the local dofs (vertices, then edge midpoints). |
| fn | `gradients_at_nodes` | `fields.rs:30` | `pub fn gradients_at_nodes(mesh: &SimplexMesh, u: &[f64], ncomp: usize, order: usize) -> Vec<f64>` —  Gradients of `ncomp` stacked scalar dof fields (each `n_dofs` long) at the local dofs of |
| fn | `node_average` | `fields.rs:74` | `pub fn node_average(mesh: &SimplexMesh, vals: &[f64], m: usize, order: usize) -> Vec<f64>` —  Mean over the cells touching each dof of per-cell nodal values `vals[(c*nloc + a)*m + j]`. |
| struct | `Display` | `fields.rs:104` | fields: `dim`, `dofs`, `points`, `tris`, `tri_label` —  Surface mesh for plotting: the cells in 2D, the boundary facets in 3D, with every P2 |
| struct | `CellField` | `fields.rs:171` | fields: `name`, `unit`, `re`, `im` —  One quantity at the local dofs of every cell, not averaged over the cells that share a node: |
| struct | `CellFields` | `fields.rs:181` | fields: `order`, `nloc`, `fields` —  The un-averaged per-cell fields of one solution. This is what an error estimator needs (the jump of a flux |
| fn | `Display::new` | `fields.rs:115` | `pub fn new(mesh: &SimplexMesh, order: usize) -> Self` |
| fn | `Display::pick` | `fields.rs:161` | `pub fn pick(&self, values: &[f64]) -> Vec<f32>` |
| fn | `CellFields::new` | `fields.rs:189` | `pub fn new(mesh: &SimplexMesh, order: usize) -> Self` |
| fn | `CellFields::get` | `fields.rs:193` | `pub fn get(&self, name: &str) -> Option<&CellField>` |
| fn | `CellFields::push_real` | `fields.rs:197` | `pub fn push_real(&mut self, name: impl Into<String>, unit: impl Into<String>, re: Vec<f64>)` |
| fn | `CellFields::push_complex` | `fields.rs:201` | `pub fn push_complex(&mut self, name: impl Into<String>, unit: impl Into<String>, re: Vec<f64>, im: Vec<f64>)` |
| fn | `CellFields::push_components` | `fields.rs:207` | `pub fn push_components(&mut self, prefix: &str, names: &[&str], unit: &str, g: &[f64])` —  The `ncomp` interleaved per-cell values `g[(c*nloc + a)*ncomp + j]` as one real field per component, named |
| fn | `CellFields::node_average` | `fields.rs:216` | `pub fn node_average(&self, mesh: &SimplexMesh, name: &str) -> Option<Vec<f64>>` —  The node average of a real component (the plot field it must reproduce). |
| fn | `phase_factor` | `periodic.rs:17` | `pub fn phase_factor(phase: f64) -> C` —  The Floquet factor `exp(-i phase)`. |
| struct | `Match` | `periodic.rs:22` | fields: `map`, `src_facets`, `shift` —  Vertex correspondence of a periodic pair. |
| fn | `match_vertices` | `periodic.rs:46` | `pub fn match_vertices(mesh: &SimplexMesh, p: &Periodic) -> Result<Match, String>` —  Matches the vertices of the source boundaries to those of the destination boundaries. |
| fn | `nodal_relations` | `periodic.rs:113` | `pub fn nodal_relations(mesh: &SimplexMesh, order: usize, ncomp: usize, m: &Match, phase: f64) -> Result<Vec<(usize, usize, C)>, String>` —  Relations for nodal fields of `ncomp` component blocks (dof layout: component blocks of |
| struct | `Reduction` | `periodic.rs:184` | fields: `n`, `roots`, `index_of`, `ratio` —  `x_full = P x_red`: the reduced dofs are the roots of the constraint classes that contain no |
| fn | `ccsr_mul` | `periodic.rs:247` | `pub fn ccsr_mul(a: &CCsr, b: &CCsr) -> CCsr` —  Sparse complex product. |
| fn | `ccsr_adjoint` | `periodic.rs:278` | `pub fn ccsr_adjoint(a: &CCsr) -> CCsr` —  Conjugate transpose. |
| fn | `congruence` | `periodic.rs:284` | `pub fn congruence(p: &CCsr, a: &CCsr) -> CCsr` —  `P^H A P`. |
| fn | `Reduction::build` | `periodic.rs:195` | `pub fn build(n: usize, relations: &[(usize, usize, C)], fixed: &[usize]) -> Reduction` |
| fn | `Reduction::n_reduced` | `periodic.rs:224` | `pub fn n_reduced(&self) -> usize` |
| fn | `Reduction::prolong` | `periodic.rs:229` | `pub fn prolong(&self) -> CCsr` —  The prolongation `P` (`n x n_reduced`). |
| fn | `Reduction::expand` | `periodic.rs:241` | `pub fn expand(&self, xr: &[C]) -> Vec<C>` |
| fn | `stretch` | `pml.rs:13` | `pub fn stretch(p: &Pml, dirs: [PmlDir; 2], centre: [f64; 2], lo: [f64; 2], hi: [f64; 2]) -> (C, C)` —  Stretch factors along the two axes of the mesh (`dirs`, the layer directions of each) of a cell centre for a PML specification and the |
| fn | `radial_coordinate` | `pml.rs:31` | `pub fn radial_coordinate(p: &Pml, r: f64, lo: f64, hi: f64) -> C` —  The complex radius `r~ = int s_r dr` (from the interface of the layer) of a point at radius `r` for a layer whose `x` direction is the |
| fn | `cell_stretches` | `pml.rs:49` | `pub fn cell_stretches(mesh: &SimplexMesh, pml: &std::collections::BTreeMap<i32, Pml>) -> Vec<[C` —  Complex stretch factors `(sx, sy, sz)` of every cell of a 2D or 3D mesh for the PML domains of |
| struct | `Materials` | `props.rs:17` | fields: `domains`, `default`, `aniso`, `pml` —  Material properties per domain; `default` applies to domains without their own value. |
| enum | `PmlDir` | `props.rs:31` | variants: `Off`, `Plus`, `Minus` —  Which side of a domain a perfectly matched layer absorbs towards, per axis. |
| struct | `Pml` | `props.rs:43` | fields: `x`, `y`, `z`, `strength`, `order`, `re` —  A perfectly matched layer by coordinate stretching, `s = re - j strength t^order` with `t` |
| enum | `TensorKind` | `props.rs:65` | variants: `Conductivity`, `Permittivity`, `Expansion` |
| struct | `Periodic` | `props.rs:180` | fields: `source`, `dest`, `phase` —  A periodic condition with a Floquet phase between two boundaries of equal meshes: the dofs on the |
| enum | `LinearSolver` | `props.rs:189` | variants: `Auto`, `Direct`, `Iterative` —  How the large symmetric positive definite systems of static problems are solved. |
| fn | `Materials::has_pml` | `props.rs:90` | `pub fn has_pml(&self) -> bool` —  True if any domain has an active perfectly matched layer. |
| fn | `Materials::get` | `props.rs:94` | `pub fn get(&self, domain: i32, name: &str) -> Option<f64>` |
| fn | `Materials::set` | `props.rs:97` | `pub fn set(&mut self, domain: i32, name: &str, v: f64)` |
| fn | `Materials::has_tensor` | `props.rs:101` | `pub fn has_tensor(&self, kind: TensorKind) -> bool` —  True if any domain defines the tensor `kind` anisotropically. |
| fn | `Materials::has_stiffness` | `props.rs:105` | `pub fn has_stiffness(&self) -> bool` |
| fn | `Materials::cell_stiffness` | `props.rs:111` | `pub fn cell_stiffness(&self, mesh: &SimplexMesh) -> Result<Vec<M6>, String>` —  Stiffness matrix (solver Voigt order) of every cell: crystal data where defined, else |
| fn | `Materials::cell_stiffness_skipping` | `props.rs:116` | `pub fn cell_stiffness_skipping(&self, mesh: &SimplexMesh, skip: &[i32]) -> Result<Vec<M6>, String>` —  Stiffness of every cell; domains in `skip` (not part of this physics) get a zero matrix. |
| fn | `Materials::cell_tensor` | `props.rs:146` | `pub fn cell_tensor(&self, mesh: &SimplexMesh, kind: TensorKind, fallback: Option<f64>) -> Result<Vec<M3>, String>` —  Symmetric tensor property of every cell (global frame): tensor data where defined, |
| fn | `Materials::cell_values` | `props.rs:163` | `pub fn cell_values(&self, mesh: &SimplexMesh, name: &str, fallback: Option<f64>) -> Result<Vec<f64>, String>` —  One value of `name` per cell. Without a `fallback`, a domain lacking the property is an error. |
| fn | `LinearSolver::name` | `props.rs:199` | `pub fn name(self) -> &'static str` |

Trait implementations in this component: `Pml: Default`
