# M7: 3D PML, native (hex / prism) against split (tet)

Harness: `validate::helmholtz` generalised to 3D: P1 tetrahedra, trilinear hexahedra, P1 x linear prisms, complex coordinate
stretching `s_a = 1 + i sigma_a / k` with `sigma = sigma_max (d/t)^2`, `R0 = 1e-8`, Dirichlet on the outer box, banded LU.
All meshes come from `build_pml` around a tetrahedral physical box (`box_tet`); the "split" meshes are the same meshes after
`split_to_simplices` (same vertices, hexes and prisms cut into tets).

## 1. Dirichlet eigenvalue of the layered mesh (stretching off)

With the stretching switched off the layers are plain material and the mesh is a box of side `L = 1 + 2t`, `t = 0.5`, so the
lowest eigenvalue is exactly `3 pi^2 / L^2 = 7.402203`. Inverse iteration, test `mixed_pml_mesh_native_and_split_eigenvalues_agree_and_converge`.

| physical cells | layers | native | split | exact | native rel err | split rel err |
|---|---|---|---|---|---|---|
| 4 | 2 | 7.503811 | 7.857685 | 7.402203 | 1.37e-2 | 6.15e-2 |
| 6 | 3 | 7.448654 | 7.604510 | 7.402203 | 6.28e-3 | 2.73e-2 |
| 8 | 4 | 7.428638 | 7.515974 | 7.402203 | 3.57e-3 | 1.54e-2 |

Both converge at second order, and they agree with each other to 1.2 % at the finest level. They are **not** the same
discretisation, so the eigenvalues are not equal: the native elements are 4.3-4.5x more accurate on the same vertices (the likely
cause is the richer trilinear/prism interpolation inside each hex or prism compared with piecewise-linear tets; not isolated further).

## 2. PML truncation error (scattered-field solve, k = 2 pi)

`cargo run --release --example pml_3d` (6^3 physical cells, 343 physical vertices). Error is the relative L2 difference over the physical
vertices against a reference of the same kind (thickness 0.6, 6 layers), so the physical discretisation error cancels.

reference (native): thickness 0.6, 6 layers, 6859 dof, band 955, 8.2s
reference (split): 6859 dof, band 776, 5.2s
native-reference vs split-reference (physical vertices): 2.37e-2

| thickness | layers | dof | native err | native s | split err | split s |
|---|---|---|---|---|---|---|
| 0.3 | 1 | 729 | 1.60e-1 | 0.0 | 2.97e-1 | 0.0 |
| 0.3 | 2 | 1331 | 2.70e-1 | 0.2 | 2.06e-1 | 0.1 |
| 0.3 | 3 | 2197 | 5.27e-2 | 0.3 | 6.73e-2 | 0.2 |
| 0.3 | 5 | 4913 | 1.85e-2 | 2.6 | 2.56e-2 | 2.2 |
| 0.5 | 1 | 729 | 2.72e-1 | 0.0 | 1.68e-1 | 0.0 |
| 0.5 | 2 | 1331 | 3.79e-1 | 0.1 | 2.49e-1 | 0.1 |
| 0.5 | 3 | 2197 | 5.11e-2 | 0.3 | 4.91e-2 | 0.3 |
| 0.5 | 5 | 4913 | 8.13e-3 | 2.7 | 1.06e-2 | 2.0 |

* 1-2 layers are not resolved and give 16-38 % errors, non-monotone in the layer count (a layer cell is then a large part of a wavelength).
* From 3 layers the native and split meshes behave alike (5 % at 3 layers, 0.8-1.9 % at 5 layers); native is slightly better at 5 layers
  (0.81 % vs 1.06 % at t = 0.5) and uses half the cells.
* The two references differ by 2.4 % (physical vertices), the same order as the discretisation error, which is why each kind is
  compared with its own reference.
* Not studied: R0 and grading tuning (done in 2D in M2), non-box PML, vector Maxwell. A scalar Helmholtz harness cannot show
  edge-element behaviour at hex/prism/tet interfaces. That check needs the nyx `emw3d` solver.
