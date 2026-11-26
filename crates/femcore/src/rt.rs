// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! H(div) Raviart-Thomas finite elements (`RT_0` for `order = 1`, `RT_1` for `order = 2`) on 2D triangles and 3D
//! tetrahedra, Braess-Schoberl / Ern-Vohralik vertex-patch flux equilibration across conductive sheets and dielectric
//! interfaces, and Prager-Synge / Thomson complementary capacitance matrix bounds.
//!
//! This module is pure physics-agnostic finite-element mathematics over [`SimplexMesh`], nodal/edge potential vectors,
//! scalar or symmetric-tensor permittivities, and boundary/sheet facet masks (`femcore` depends only on `geom` in the
//! workspace).

use crate::fem::{facet_rule, quad_rule, quad_rule_collapsed, shape, CurvedTet};
use crate::linalg::driver::dense_eigh;
use crate::material::M3;
use crate::mesh::SimplexMesh;
use std::collections::{BTreeMap, VecDeque};

/// Per-cell permittivity `eps` (in F/m, e.g. `eps0 * eps_r` or `eps0 * eps_tensor`).
#[derive(Clone, Copy)]
pub enum CellPermittivity<'a> {
    /// Isotropic scalar permittivity `eps[c] > 0` per cell (`mesh.n_cells()`).
    Scalar(&'a [f64]),
    /// Symmetric positive-definite tensor permittivity `eps[c]` per cell (`mesh.n_cells()`).
    /// On 2D cross-section meshes (`dim == 2`), the leading `2x2` block is used, matching [`crate::fem::assemble_laplace_tensor`].
    Tensor(&'a [M3]),
}

/// Optional exterior Robin (`FarField`) boundary specification `D . n = alpha * V` (`alpha > 0`) on `Gamma_inf`.
///
/// Facets where `alpha <= 0.0` act as homogeneous Neumann boundaries (`D . n = 0`).
#[derive(Clone, Copy, Default)]
pub enum RobinBoundary<'a> {
    #[default]
    None,
    /// Constant `alpha[f] >= 0` per mesh facet (`0..mesh.n_facets()`).
    FacetAlpha(&'a [f64]),
    /// Spatially varying `alpha(facet_idx, x_phys, outward_unit_normal) >= 0` on mesh facets.
    Pointwise(&'a (dyn Fn(usize, [f64; 3], [f64; 3]) -> f64 + Sync)),
}

/// Output of `RT_{p-1}` vertex-patch equilibration and Thomson / Prager-Synge complementary matrix bounds.
#[derive(Clone, Debug)]
pub struct EquilibratedFluxResult {
    /// Primal Lagrange polynomial order (`1` or `2`).
    pub order: usize,
    /// Raviart-Thomas polynomial degree `p - 1` (`0` for `RT_0`, `1` for `RT_1`).
    pub rt_order: usize,
    /// Spatial dimension (`2` or `3`).
    pub dim: usize,
    /// Total number of connected vertex half-star patches solved.
    pub n_patches: usize,
    /// Number of vertices whose star split into `>= 2` disconnected half-stars across interior conductive sheets.
    pub n_sheet_splits: usize,
    /// Maximum global `L^2(Omega)` divergence norm `max_j ||div(D_tilde^(j))||_{L^2}` across all terminals.
    pub max_div_l2: f64,
    /// Cellwise local `RT_{p-1}` basis coefficients for each terminal `j`: length `n_cells * rt_nloc(dim, order)`.
    pub cell_rt_coeffs: Vec<Vec<f64>>,
    /// Dual complementary flux energy matrix `W^flux_{ij} = int D_tilde^(i) . eps^-1 D_tilde^(j) dOmega + int_{Gamma_inf} alpha^-1 (D_tilde^(i).n)(D_tilde^(j).n) dS`.
    pub w_flux_matrix: Vec<Vec<f64>>,
    /// Mixed potential-flux pairing matrix `Q_{ij} = -int grad(V_h^(i)) . D_tilde^(j) dOmega + int_{Gamma_inf} V_h^(i) (D_tilde^(j).n) dS`.
    pub q_pairing_matrix: Vec<Vec<f64>>,
    /// Primal symmetric PSD capacitance upper bound matrix `C^up`.
    pub c_upper_matrix: Vec<Vec<f64>>,
    /// Guaranteed symmetric PSD Loewner lower bound matrix `C^low = Q (W^flux)^+ Q^T <= C^exact <= C^up`.
    pub c_lower_matrix: Vec<Vec<f64>>,
    /// Complementary midpoint estimate `C^mid = 0.5 * (C^up + C^low)`.
    pub c_mid_matrix: Vec<Vec<f64>>,
    /// Rigorous entrywise error bound `delta_C_{ij} = 0.5 * sqrt((C^up_{ii} - C^low_{ii}) * (C^up_{jj} - C^low_{jj}))` on `|C^mid_{ij} - C^exact_{ij}|`.
    pub c_error_bound_matrix: Vec<Vec<f64>>,
    /// Per-terminal, per-cell Prager-Synge error indicator `eta^2_{T,PS,j} = ||eps^{-1/2}(D_tilde^(j) - D_h^(j))||^2_{L^2(T)}`.
    pub terminal_cell_eta2: Vec<Vec<f64>>,
    /// Goal-oriented Cauchy-Schwarz combined per-cell Prager-Synge error indicator `eta^2_{T,PS}` (`len == n_cells`).
    pub cell_eta2: Vec<f64>,
}

/// Number of local degrees of freedom of `RT_{order - 1}` on a `dim`-simplex (`dim` in `{2, 3}`, `order` in `{1, 2}`).
#[inline]
pub fn rt_nloc(dim: usize, order: usize) -> usize {
    match (dim, order) {
        (2, 1) => 3,
        (3, 1) => 4,
        (2, 2) => 8,
        (3, 2) => 15,
        _ => panic!("RT_{order}-1 requires dim in {{2, 3}} and order in {{1, 2}}"),
    }
}

/// Evaluates the reference `RT_{order - 1}(hat_K)` basis vectors and reference divergences at barycentric coordinates `lam` (`len == dim + 1`).
///
/// Basis ordering (with cell vertices sorted in ascending global index order `u_0 < u_1 < ... < u_d`):
/// * `order == 1` (`RT_0`, `d + 1` face functions):
///   `phi_m(hat_x) = (d - 1)! * (hat_x - hat_v_m)` for `m = 0..=d`, with `div_hat(phi_m) = d!`.
/// * `order == 2` (`RT_1`, `d*(d + 1)` face functions followed by `d` interior cell-bubble functions):
///   - For each face `m = 0..=d` and each vertex `a != m` in ascending order (`d` functions per face):
///     `phi_{m,a}(hat_x) = (d - 1)! * lam[a] * (hat_x - hat_v_m)`, with `div_hat(phi_{m,a}) = (d - 1)! * (d + 1) * lam[a]`.
///   - For `m = 1..=d` (`d` interior bubble functions with identically zero normal trace on all `d + 1` faces):
///     `phi^{bub}_m(hat_x) = (d - 1)! * lam[m] * (hat_x - hat_v_m)`, with `div_hat(phi^{bub}_m) = (d - 1)! * ((d + 1) * lam[m] - 1)`.
pub fn eval_rt_ref(
    dim: usize,
    order: usize,
    lam: &[f64],
    out_vec: &mut [[f64; 3]],
    out_div: &mut [f64],
) {
    debug_assert_eq!(lam.len(), dim + 1);
    let fac = if dim == 3 { 2.0 } else { 1.0 }; // (d - 1)!
    let mut x_hat = [0.0; 3];
    x_hat[..dim].copy_from_slice(&lam[1..(dim + 1)]);

    if order == 1 {
        let div_c = fac * dim as f64; // d!
        for m in 0..=dim {
            let mut v = x_hat;
            if m >= 1 {
                v[m - 1] -= 1.0;
            }
            out_vec[m] = [fac * v[0], fac * v[1], fac * v[2]];
            out_div[m] = div_c;
        }
        return;
    }

    let dp1 = (dim + 1) as f64;
    let mut idx = 0;
    for m in 0..=dim {
        let mut diff = x_hat;
        if m >= 1 {
            diff[m - 1] -= 1.0;
        }
        for a in 0..=dim {
            if a == m {
                continue;
            }
            let s = fac * lam[a];
            out_vec[idx] = [s * diff[0], s * diff[1], s * diff[2]];
            out_div[idx] = fac * dp1 * lam[a];
            idx += 1;
        }
    }
    for m in 1..=dim {
        let mut diff = x_hat;
        diff[m - 1] -= 1.0;
        let s = fac * lam[m];
        out_vec[idx] = [s * diff[0], s * diff[1], s * diff[2]];
        out_div[idx] = fac * (dp1 * lam[m] - 1.0);
        idx += 1;
    }
}

/// Evaluates contravariant Piola-mapped physical `RT_{order - 1}(K)` basis vectors `(1 / |det J|) J * phi_hat`
/// and physical divergences `(1 / |det J|) div_hat(phi_hat)` at barycentric point `lam`.
///
/// Here `j_mat[a][b] = dx_a / dxi_b` (`0 <= a, b < dim`) is the physical-by-reference Jacobian matrix and `det_j = det(j_mat)`.
pub fn eval_rt_piola(
    dim: usize,
    order: usize,
    j_mat: &[[f64; 3]; 3],
    det_j: f64,
    lam: &[f64],
    out_vec: &mut [[f64; 3]],
    out_div: &mut [f64],
) {
    let nrt = rt_nloc(dim, order);
    let mut ref_vec = [[0.0; 3]; 15];
    let mut ref_div = [0.0; 15];
    eval_rt_ref(dim, order, lam, &mut ref_vec[..nrt], &mut ref_div[..nrt]);
    let inv_abs_det = 1.0 / det_j.abs();
    for r in 0..nrt {
        let rv = ref_vec[r];
        let mut pv = [0.0; 3];
        for a in 0..dim {
            let mut s = 0.0;
            for b in 0..dim {
                s += j_mat[a][b] * rv[b];
            }
            pv[a] = s * inv_abs_det;
        }
        out_vec[r] = pv;
        out_div[r] = ref_div[r] * inv_abs_det;
    }
}

#[derive(Clone)]
struct SortedCell {
    verts: [u32; 4],
    dofs: [usize; 10],
    curved: Option<CurvedTet>,
    j_aff: [[f64; 3]; 3],
    det_aff: f64,
    grad_lam_aff: [[f64; 3]; 4],
    eps_inv: [[f64; 3]; 3],
    eps_fwd: [[f64; 3]; 3],
}

type Mat3Pair = ([[f64; 3]; 3], [[f64; 3]; 3]);
type AffineJacGrads = ([[f64; 3]; 3], f64, [[f64; 3]; 4]);

fn invert_spd_cell_eps(dim: usize, eps: CellPermittivity<'_>, c: usize) -> Result<Mat3Pair, String> {
    match eps {
        CellPermittivity::Scalar(s) => {
            let e = s[c];
            if !e.is_finite() || e <= 0.0 {
                return Err(format!("cell {c} has non-positive or non-finite permittivity {e}"));
            }
            let inv = 1.0 / e;
            let mut fwd = [[0.0; 3]; 3];
            let mut inv_m = [[0.0; 3]; 3];
            for k in 0..dim {
                fwd[k][k] = e;
                inv_m[k][k] = inv;
            }
            Ok((fwd, inv_m))
        }
        CellPermittivity::Tensor(ts) => {
            let t = ts[c];
            for i in 0..dim {
                for j in 0..dim {
                    if !t[i][j].is_finite() {
                        return Err(format!("cell {c} has non-finite permittivity tensor entry ({i}, {j})"));
                    }
                }
            }
            if dim == 2 {
                let (a00, a01, a11) = (t[0][0], 0.5 * (t[0][1] + t[1][0]), t[1][1]);
                let det = a00 * a11 - a01 * a01;
                if a00 <= 0.0 || det <= 0.0 {
                    return Err(format!("cell {c} 2D permittivity tensor is not positive definite"));
                }
                let fwd = [[a00, a01, 0.0], [a01, a11, 0.0], [0.0, 0.0, 0.0]];
                let inv_m = [[a11 / det, -a01 / det, 0.0], [-a01 / det, a00 / det, 0.0], [0.0, 0.0, 0.0]];
                Ok((fwd, inv_m))
            } else {
                let a00 = t[0][0];
                let a01 = 0.5 * (t[0][1] + t[1][0]);
                let a02 = 0.5 * (t[0][2] + t[2][0]);
                let a11 = t[1][1];
                let a12 = 0.5 * (t[1][2] + t[2][1]);
                let a22 = t[2][2];
                let minor2 = a00 * a11 - a01 * a01;
                let c00 = a11 * a22 - a12 * a12;
                let c01 = a02 * a12 - a01 * a22;
                let c02 = a01 * a12 - a02 * a11;
                let c11 = a00 * a22 - a02 * a02;
                let c12 = a01 * a02 - a00 * a12;
                let c22 = minor2;
                let det = a00 * c00 + a01 * c01 + a02 * c02;
                if a00 <= 0.0 || minor2 <= 0.0 || det <= 0.0 {
                    return Err(format!("cell {c} 3D permittivity tensor is not positive definite"));
                }
                let fwd = [[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]];
                let inv_m = [
                    [c00 / det, c01 / det, c02 / det],
                    [c01 / det, c11 / det, c12 / det],
                    [c02 / det, c12 / det, c22 / det],
                ];
                Ok((fwd, inv_m))
            }
        }
    }
}

fn affine_jacobian_and_grads(dim: usize, pts: &[[f64; 3]; 4]) -> Result<AffineJacGrads, String> {
    let mut j = [[0.0; 3]; 3];
    for b in 0..dim {
        for a in 0..dim {
            j[a][b] = pts[b + 1][a] - pts[0][a];
        }
    }
    let mut g = [[0.0; 3]; 4];
    if dim == 2 {
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() <= 1.0e-300 {
            return Err("degenerate zero-area 2D triangle in mesh".into());
        }
        // grad lam_1 = row 0 of J^-1 = [j11, -j01] / det; grad lam_2 = row 1 of J^-1 = [-j10, j00] / det
        g[1][0] = j[1][1] / det;
        g[1][1] = -j[0][1] / det;
        g[2][0] = -j[1][0] / det;
        g[2][1] = j[0][0] / det;
        g[0][0] = -(g[1][0] + g[2][0]);
        g[0][1] = -(g[1][1] + g[2][1]);
        Ok((j, det, g))
    } else {
        let c00 = j[1][1] * j[2][2] - j[1][2] * j[2][1];
        let c01 = j[1][2] * j[2][0] - j[1][0] * j[2][2];
        let c02 = j[1][0] * j[2][1] - j[1][1] * j[2][0];
        let det = j[0][0] * c00 + j[0][1] * c01 + j[0][2] * c02;
        if det.abs() <= 1.0e-300 {
            return Err("degenerate zero-volume 3D tetrahedron in mesh".into());
        }
        // J^-1 rows: row 0 = (c00, c10, c20)/det, where c_{ab} is cofactor of J[a][b]
        let c10 = j[2][1] * j[0][2] - j[0][1] * j[2][2];
        let c20 = j[0][1] * j[1][2] - j[1][1] * j[0][2];
        let c11 = j[0][0] * j[2][2] - j[0][2] * j[2][0];
        let c21 = j[1][0] * j[0][2] - j[0][0] * j[1][2];
        let c12 = j[2][0] * j[0][1] - j[0][0] * j[2][1];
        let c22 = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        // Wait: (J^-1)_{b, a} = C_{a, b} / det, and grad(lam_{b+1})_a = (J^-1)_{b, a} = C_{a, b} / det!
        // Check: sum_a (dx_a / dxi_k) * grad(lam_{b+1})_a = sum_a J[a][k] * (J^-1)[b][a] = delta_{b, k}.
        g[1] = [c00 / det, c10 / det, c20 / det];
        g[2] = [c01 / det, c11 / det, c21 / det];
        g[3] = [c02 / det, c12 / det, c22 / det];
        for a in 0..3 {
            g[0][a] = -(g[1][a] + g[2][a] + g[3][a]);
        }
        Ok((j, det, g))
    }
}

fn build_sorted_cells(mesh: &SimplexMesh, order: usize, eps: CellPermittivity<'_>) -> Result<Vec<SortedCell>, String> {
    let d = mesh.dim;
    let nc = mesh.n_cells();
    let nv = mesh.n_points();
    let cquad_check = if d == 3 && mesh.has_curved() { Some(quad_rule_collapsed(6)) } else { None };
    let mut out = Vec::with_capacity(nc);

    for c in 0..nc {
        let raw = mesh.cell(c);
        let mut verts = [u32::MAX; 4];
        for (i, &v) in raw.iter().enumerate() {
            verts[i] = v;
        }
        verts[..=d].sort_unstable();

        let mut pts = [[0.0; 3]; 4];
        for i in 0..=d {
            let p = mesh.point(verts[i] as usize);
            pts[i][..d].copy_from_slice(&p[..d]);
        }
        let (j_aff, det_aff, grad_lam_aff) = affine_jacobian_and_grads(d, &pts)?;
        let (eps_fwd, eps_inv) = invert_spd_cell_eps(d, eps, c)?;

        let mut dofs = [0usize; 10];
        for i in 0..=d {
            dofs[i] = verts[i] as usize;
        }
        if order == 2 {
            let mut k = d + 1;
            for i in 0..=d {
                for j in (i + 1)..=d {
                    dofs[k] = nv + mesh.edge_index(verts[i], verts[j]);
                    k += 1;
                }
            }
        }

        let curved = if d == 3 && mesh.cell_is_curved(c) {
            let mut x10 = [[0.0; 3]; 10];
            x10[..4].copy_from_slice(&pts);
            let mut k = 4;
            for i in 0..4 {
                for j in (i + 1)..4 {
                    x10[k] = mesh.edge_node(mesh.edge_index(verts[i], verts[j])).unwrap_or_else(|| {
                        [
                            0.5 * (pts[i][0] + pts[j][0]),
                            0.5 * (pts[i][1] + pts[j][1]),
                            0.5 * (pts[i][2] + pts[j][2]),
                        ]
                    });
                    k += 1;
                }
            }
            let ct = CurvedTet::new(x10);
            // Cross-Feature Compatibility Rule 3: contravariant Piola RT vectors have divergent L^2 energy
            // on singular quarter-point CurvedTets (det J -> 0 at sheet rims).
            let mut node_lams = vec![
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
            for &(i, j) in crate::fem::TET_EDGES.iter() {
                let mut l = [0.0; 4];
                l[i] = 0.5;
                l[j] = 0.5;
                node_lams.push(l);
            }
            for l in &node_lams {
                let (_, det_node) = ct.jacobian(l);
                if det_node * det_aff <= 1.0e-6 * det_aff * det_aff {
                    return Err(
                        "flux_bounds (H(div) Raviart-Thomas equilibration) is mathematically incompatible with singular quarter-point CurvedTets (det(J) -> 0 at sheet rims causes divergent contravariant Piola L2 energy); disable singular_edges=\"quarter_point\" or flux_bounds"
                            .into(),
                    );
                }
            }
            if let Some(cq) = &cquad_check {
                for l in &cq.lam {
                    let (_, det_q) = ct.jacobian(&[l[0], l[1], l[2], l[3]]);
                    if det_q * det_aff <= 0.0 {
                        return Err(format!("curved cell {c} has non-positive Jacobian determinant during H(div) equilibration"));
                    }
                }
            }
            Some(ct)
        } else {
            None
        };

        out.push(SortedCell {
            verts,
            dofs,
            curved,
            j_aff,
            det_aff,
            grad_lam_aff,
            eps_inv,
            eps_fwd,
        });
    }
    Ok(out)
}

#[inline]
fn sorted_facet_key(fv: &[u32]) -> [u32; 3] {
    let mut k = [u32::MAX; 3];
    for (i, &v) in fv.iter().enumerate() {
        k[i] = v;
    }
    k[..fv.len()].sort_unstable();
    k
}

#[inline]
fn cell_face_key(dim: usize, sorted_verts: &[u32; 4], opp_slot: usize) -> [u32; 3] {
    let mut k = [u32::MAX; 3];
    let mut p = 0;
    for a in 0..=dim {
        if a != opp_slot {
            k[p] = sorted_verts[a];
            p += 1;
        }
    }
    k
}

#[derive(Clone, Copy, Debug)]
struct FacetTag {
    is_dirichlet: bool,
    is_sheet: bool,
    robin_facet: Option<usize>,
}

/// Evaluates outward unit normal `n_out`, physical position `x_phys`, and quadrature surface measure `ds_q = |J_F| * w_q`
/// on face `opp_slot` of `sc` at facet barycentric coordinates `mu` (`len == dim`).
fn eval_cell_face_point(
    dim: usize,
    sc: &SortedCell,
    opp_slot: usize,
    mu: &[f64],
    w_q: f64,
    pts_mesh: &[f64],
) -> ([f64; 4], [f64; 3], [f64; 3], f64) {
    let mut lam = [0.0; 4];
    let mut face_slots = [0usize; 3];
    let mut p = 0;
    for a in 0..=dim {
        if a != opp_slot {
            lam[a] = mu[p];
            face_slots[p] = a;
            p += 1;
        }
    }

    if dim == 2 {
        let v0 = sc.verts[face_slots[0]] as usize;
        let v1 = sc.verts[face_slots[1]] as usize;
        let v_opp = sc.verts[opp_slot] as usize;
        let p0 = [pts_mesh[2 * v0], pts_mesh[2 * v0 + 1], 0.0];
        let p1 = [pts_mesh[2 * v1], pts_mesh[2 * v1 + 1], 0.0];
        let popp = [pts_mesh[2 * v_opp], pts_mesh[2 * v_opp + 1], 0.0];
        let x_phys = [
            mu[0] * p0[0] + mu[1] * p1[0],
            mu[0] * p0[1] + mu[1] * p1[1],
            0.0,
        ];
        let tx = p1[0] - p0[0];
        let ty = p1[1] - p0[1];
        let len = (tx * tx + ty * ty).sqrt();
        let mut nx = ty / len;
        let mut ny = -tx / len;
        if nx * (popp[0] - p0[0]) + ny * (popp[1] - p0[1]) > 0.0 {
            nx = -nx;
            ny = -ny;
        }
        (lam, x_phys, [nx, ny, 0.0], len * w_q)
    } else if let Some(ct) = &sc.curved {
        let x_phys = ct.position(&lam);
        // Surface tangents d(position)/d(mu_1) and d(position)/d(mu_2) where mu_0 = 1 - mu_1 - mu_2:
        // d(lam)/d(mu_k) has -1 at face_slots[0] and +1 at face_slots[k] (k = 1, 2).
        // Using j_mat[a][b] = dx_a / d(lam_{b+1}) with lam_0 = 1 - lam_1 - lam_2 - lam_3:
        let j_mat = ct.jacobian_matrix(&lam);
        let col_of_slot = |slot: usize| -> [f64; 3] {
            if slot == 0 {
                [0.0, 0.0, 0.0]
            } else {
                [j_mat[0][slot - 1], j_mat[1][slot - 1], j_mat[2][slot - 1]]
            }
        };
        let c0 = col_of_slot(face_slots[0]);
        let c1 = col_of_slot(face_slots[1]);
        let c2 = col_of_slot(face_slots[2]);
        let t1 = [c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2]];
        let t2 = [c2[0] - c0[0], c2[1] - c0[1], c2[2] - c0[2]];
        let cr = [
            t1[1] * t2[2] - t1[2] * t2[1],
            t1[2] * t2[0] - t1[0] * t2[2],
            t1[0] * t2[1] - t1[1] * t2[0],
        ];
        let norm_cr = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
        // Outward direction matches -grad(lam[opp_slot])!
        let (gq, _) = ct.jacobian(&lam);
        let g_opp = gq[opp_slot];
        let sign = if cr[0] * g_opp[0] + cr[1] * g_opp[1] + cr[2] * g_opp[2] < 0.0 { 1.0 } else { -1.0 };
        let n_out = [sign * cr[0] / norm_cr, sign * cr[1] / norm_cr, sign * cr[2] / norm_cr];
        (lam, x_phys, n_out, 0.5 * norm_cr * w_q)
    } else {
        let v0 = sc.verts[face_slots[0]] as usize;
        let v1 = sc.verts[face_slots[1]] as usize;
        let v2 = sc.verts[face_slots[2]] as usize;
        let v_opp = sc.verts[opp_slot] as usize;
        let pt = |v: usize| [pts_mesh[3 * v], pts_mesh[3 * v + 1], pts_mesh[3 * v + 2]];
        let (p0, p1, p2, popp) = (pt(v0), pt(v1), pt(v2), pt(v_opp));
        let x_phys = [
            mu[0] * p0[0] + mu[1] * p1[0] + mu[2] * p2[0],
            mu[0] * p0[1] + mu[1] * p1[1] + mu[2] * p2[1],
            mu[0] * p0[2] + mu[1] * p1[2] + mu[2] * p2[2],
        ];
        let t1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let t2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let cr = [
            t1[1] * t2[2] - t1[2] * t2[1],
            t1[2] * t2[0] - t1[0] * t2[2],
            t1[0] * t2[1] - t1[1] * t2[0],
        ];
        let norm_cr = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
        let sign = if cr[0] * (popp[0] - p0[0]) + cr[1] * (popp[1] - p0[1]) + cr[2] * (popp[2] - p0[2]) < 0.0 {
            1.0
        } else {
            -1.0
        };
        let n_out = [sign * cr[0] / norm_cr, sign * cr[1] / norm_cr, sign * cr[2] / norm_cr];
        (lam, x_phys, n_out, 0.5 * norm_cr * w_q)
    }
}

fn eval_robin_alpha(
    robin: RobinBoundary<'_>,
    facet_idx: usize,
    x_phys: [f64; 3],
    n_out: [f64; 3],
) -> f64 {
    match robin {
        RobinBoundary::None => 0.0,
        RobinBoundary::FacetAlpha(a) => a.get(facet_idx).copied().unwrap_or(0.0),
        RobinBoundary::Pointwise(f) => f(facet_idx, x_phys, n_out),
    }
}

#[derive(Clone, Copy, Debug)]
struct PatchCellFace {
    /// Active patch DOF start index (`usize::MAX` if `sigma . n = 0` on this face).
    dof_offset: usize,
    /// Orientation sign (`+1.0` or `-1.0`) mapping patch face DOFs to cell outward normal DOFs.
    sign: f64,
    /// Optional exterior Robin facet index when this face has active `alpha > 0` on `Gamma_inf`.
    robin_facet: Option<usize>,
}

#[derive(Clone, Debug)]
struct VertexPatch {
    cells: Vec<usize>,
    /// Slot in `0..=dim` of vertex `v` inside each patch cell's `SortedCell::verts`.
    v_slot_in_cell: Vec<usize>,
    /// Per patch cell (`0..cells.len()`), its `dim + 1` faces (`0..=dim`).
    cell_faces: Vec<[PatchCellFace; 4]>,
    /// Per patch cell, start index in patch vector for `RT_1` interior bubble DOFs (`usize::MAX` for `RT_0`).
    bubble_offset: Vec<usize>,
    n_sigma: usize,
    /// True if at least one face in the patch is an unconstrained Dirichlet/sheet or active Robin face.
    is_open: bool,
}

/// Solves `A X = B` for SPD `n x n` matrix `a` (row-major) and `n x nrhs` right-hand side `b` (row-major `r * nrhs + col`),
/// using unit-diagonal pre-scaling + Cholesky (`L L^T`), falling back to eigenthresholded `dense_eigh` if semi-definite.
fn solve_spd_multi(a: &[f64], b: &[f64], n: usize, nrhs: usize) -> Result<Vec<f64>, String> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut scale = vec![1.0; n];
    for i in 0..n {
        let di = a[i * n + i];
        if di > 1.0e-300 {
            scale[i] = 1.0 / di.sqrt();
        }
    }
    let mut l = vec![0.0; n * n];
    let mut cholesky_ok = true;
    'chol: for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.5 * (a[i * n + j] + a[j * n + i]) * scale[i] * scale[j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 1.0e-12 {
                    cholesky_ok = false;
                    break 'chol;
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }

    if cholesky_ok {
        let mut y = vec![0.0; n * nrhs];
        for i in 0..n {
            let inv_lii = 1.0 / l[i * n + i];
            let si = scale[i];
            for col in 0..nrhs {
                let mut s = b[i * nrhs + col] * si;
                for k in 0..i {
                    s -= l[i * n + k] * y[k * nrhs + col];
                }
                y[i * nrhs + col] = s * inv_lii;
            }
        }
        let mut x = vec![0.0; n * nrhs];
        for i in (0..n).rev() {
            let inv_lii = 1.0 / l[i * n + i];
            for col in 0..nrhs {
                let mut s = y[i * nrhs + col];
                for k in (i + 1)..n {
                    s -= l[k * n + i] * x[k * nrhs + col];
                }
                x[i * nrhs + col] = s * inv_lii;
            }
        }
        for i in 0..n {
            let si = scale[i];
            for col in 0..nrhs {
                x[i * nrhs + col] *= si;
            }
        }
        return Ok(x);
    }

    // Rank-deficient / semi-definite fallback via symmetric eigendecomposition on scaled matrix
    let mut ascaled = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            ascaled[i * n + j] = 0.5 * (a[i * n + j] + a[j * n + i]) * scale[i] * scale[j];
        }
    }
    let (eigs, vecs) = dense_eigh(&ascaled, n)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let mut x = vec![0.0; n * nrhs];
    if lam_max <= 0.0 {
        return Ok(x);
    }
    let thresh = 1.0e-12 * lam_max;
    for (k, &lam) in eigs.iter().enumerate() {
        if lam > thresh {
            let vk = &vecs[k];
            let inv_l = 1.0 / lam;
            for col in 0..nrhs {
                let mut dot = 0.0;
                for i in 0..n {
                    dot += vk[i] * (b[i * nrhs + col] * scale[i]);
                }
                let coeff = dot * inv_l;
                for i in 0..n {
                    x[i * nrhs + col] += coeff * vk[i];
                }
            }
        }
    }
    for i in 0..n {
        let si = scale[i];
        for col in 0..nrhs {
            x[i * nrhs + col] *= si;
        }
    }
    Ok(x)
}

#[inline]
fn matvec_sym_d(dim: usize, m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for i in 0..dim {
        let mut s = 0.0;
        for j in 0..dim {
            s += m[i][j] * v[j];
        }
        out[i] = s;
    }
    out
}

#[inline]
fn dot_d(dim: usize, a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mut s = 0.0;
    for k in 0..dim {
        s += a[k] * b[k];
    }
    s
}

struct PrecomputedRefTab {
    w: Vec<f64>,
    lam: Vec<[f64; 4]>,
    lag_val: Vec<Vec<f64>>,
    lag_dlam: Vec<Vec<f64>>,
    rt_ref_vec: Vec<Vec<[f64; 3]>>,
    rt_ref_div: Vec<Vec<f64>>,
}

impl PrecomputedRefTab {
    #[inline]
    fn piola_at_q(
        &self,
        dim: usize,
        q: usize,
        j_mat: &[[f64; 3]; 3],
        det_j: f64,
        out_vec: &mut [[f64; 3]],
        out_div: &mut [f64],
    ) {
        let inv_abs_det = 1.0 / det_j.abs();
        let rvecs = &self.rt_ref_vec[q];
        let rdivs = &self.rt_ref_div[q];
        for r in 0..rvecs.len() {
            let rv = rvecs[r];
            let mut pv = [0.0; 3];
            for a in 0..dim {
                let mut s = 0.0;
                for b in 0..dim {
                    s += j_mat[a][b] * rv[b];
                }
                pv[a] = s * inv_abs_det;
            }
            out_vec[r] = pv;
            out_div[r] = rdivs[r] * inv_abs_det;
        }
    }
    fn new(dim: usize, order: usize, curved_3d: bool) -> Self {
        let q = if curved_3d {
            quad_rule_collapsed(6)
        } else {
            quad_rule(dim, if order == 1 { 2 } else { 4 })
        };
        let nrt = rt_nloc(dim, order);
        let mut lam = Vec::with_capacity(q.w.len());
        let mut lag_val = Vec::with_capacity(q.w.len());
        let mut lag_dlam = Vec::with_capacity(q.w.len());
        let mut rt_ref_vec = Vec::with_capacity(q.w.len());
        let mut rt_ref_div = Vec::with_capacity(q.w.len());
        for l in &q.lam {
            let mut la = [0.0; 4];
            for (k, &x) in l.iter().enumerate() {
                la[k] = x;
            }
            lam.push(la);
            let (nv, dn) = shape(l, order);
            lag_val.push(nv);
            lag_dlam.push(dn);
            let mut rv = vec![[0.0; 3]; nrt];
            let mut rd = vec![0.0; nrt];
            eval_rt_ref(dim, order, l, &mut rv, &mut rd);
            rt_ref_vec.push(rv);
            rt_ref_div.push(rd);
        }
        Self { w: q.w, lam, lag_val, lag_dlam, rt_ref_vec, rt_ref_div }
    }
}

/// Computes physical `grad(V_h)` and scalar `V_h` at quadrature point `q` of `sc` for potential vector `v_pot`.
#[inline]
fn eval_vh_and_grad(
    dim: usize,
    order: usize,
    sc: &SortedCell,
    grad_lam: &[[f64; 3]; 4],
    tab: &PrecomputedRefTab,
    q: usize,
    v_pot: &[f64],
) -> (f64, [f64; 3]) {
    let nloc = if order == 1 { dim + 1 } else { (dim + 1) + dim * (dim + 1) / 2 };
    let nv = &tab.lag_val[q];
    let dn = &tab.lag_dlam[q];
    let mut val = 0.0;
    let mut grad = [0.0; 3];
    for a in 0..nloc {
        let va = v_pot[sc.dofs[a]];
        val += va * nv[a];
        for b in 0..=dim {
            let dnab = dn[a * (dim + 1) + b];
            if dnab != 0.0 {
                let c = va * dnab;
                for k in 0..dim {
                    grad[k] += c * grad_lam[b][k];
                }
            }
        }
    }
    (val, grad)
}

/// Goal-oriented Cauchy-Schwarz combination of per-terminal Prager-Synge cell indicators `terminal_cell_eta2[j][c]`,
/// matching `AF-1` relative capacitance weights `W_{ij} = 1 / max(|C_{ij}|, eps_floor * sqrt(C_{ii} * C_{jj}))`.
pub fn goal_weighted_cell_eta2(
    terminal_cell_eta2: &[Vec<f64>],
    c_ref: &[Vec<f64>],
    eps_floor: f64,
) -> Vec<f64> {
    let n_term = terminal_cell_eta2.len();
    if n_term == 0 {
        return Vec::new();
    }
    let nc = terminal_cell_eta2[0].len();
    if n_term == 1 {
        return terminal_cell_eta2[0].clone();
    }
    let floor = eps_floor.max(1.0e-12);
    let mut weights = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        let cii = c_ref[i][i].abs();
        for j in i..n_term {
            let cjj = c_ref[j][j].abs();
            let denom = c_ref[i][j].abs().max(floor * (cii * cjj).sqrt()).max(1.0e-30);
            weights[i][j] = 1.0 / denom;
        }
    }
    let mut out = vec![0.0; nc];
    for c in 0..nc {
        let mut s = 0.0;
        for i in 0..n_term {
            let ei = terminal_cell_eta2[i][c].max(0.0).sqrt();
            for j in i..n_term {
                let ej = terminal_cell_eta2[j][c].max(0.0).sqrt();
                s += weights[i][j] * ei * ej;
            }
        }
        out[c] = s;
    }
    out
}

/// High-level electrostatic flux equilibration and complementary capacitance bound computation from boundary label lists.
///
/// Automatically classifies facets carrying `dirichlet_labels` (`Terminal`, `Ground`, `Potential`) into exterior
/// Dirichlet walls vs interior zero-thickness conductive sheets (`!exterior[f]`) so vertex stars split cleanly into
/// upper/lower half-stars across conductive sheets while enforcing normal-flux continuity across dielectric interfaces.
pub fn equilibrate_electro_flux_bounds(
    mesh: &SimplexMesh,
    order: usize,
    potentials: &[Vec<f64>],
    eps: CellPermittivity<'_>,
    dirichlet_labels: &[i32],
    robin: RobinBoundary<'_>,
    c_upper_override: Option<&[Vec<f64>]>,
) -> Result<EquilibratedFluxResult, String> {
    let nf = mesh.n_facets();
    let exterior = mesh.exterior_facets();
    let dir_set: std::collections::BTreeSet<i32> = dirichlet_labels.iter().copied().collect();
    let mut dirichlet_facets = vec![false; nf];
    let mut sheet_facets = vec![false; nf];
    for f in 0..nf {
        if dir_set.contains(&mesh.facet_labels[f]) {
            dirichlet_facets[f] = true;
            if !exterior[f] {
                sheet_facets[f] = true;
            }
        }
    }
    equilibrate_flux_patches(
        mesh,
        order,
        potentials,
        eps,
        &dirichlet_facets,
        &sheet_facets,
        robin,
        c_upper_override,
    )
}

/// Braess-Schoberl / Ern-Vohralik `RT_{order - 1}` vertex-patch flux equilibration and complementary Thomson /
/// Prager-Synge lower, midpoint, and upper capacitance matrix bounds.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_flux_patches(
    mesh: &SimplexMesh,
    order: usize,
    potentials: &[Vec<f64>],
    eps: CellPermittivity<'_>,
    dirichlet_facets: &[bool],
    sheet_facets: &[bool],
    robin: RobinBoundary<'_>,
    c_upper_override: Option<&[Vec<f64>]>,
) -> Result<EquilibratedFluxResult, String> {
    let d = mesh.dim;
    if d != 2 && d != 3 {
        return Err(format!("RT equilibration supports 2D and 3D meshes, got dim = {d}"));
    }
    if order != 1 && order != 2 {
        return Err(format!("RT equilibration supports order = 1 (RT_0) and order = 2 (RT_1), got order = {order}"));
    }
    let nc = mesh.n_cells();
    if nc == 0 {
        return Err("cannot equilibrate flux on an empty mesh".into());
    }
    let nf = mesh.n_facets();
    if dirichlet_facets.len() != nf || sheet_facets.len() != nf {
        return Err(format!(
            "facet mask lengths ({}, {}) do not match mesh.n_facets() = {nf}",
            dirichlet_facets.len(),
            sheet_facets.len()
        ));
    }
    if let RobinBoundary::FacetAlpha(a) = robin {
        if a.len() != nf {
            return Err(format!("robin facet alpha length {} does not match mesh.n_facets() = {nf}", a.len()));
        }
        if a.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err("robin facet alpha values must be finite and >= 0".into());
        }
    }
    let n_term = potentials.len();
    if n_term == 0 {
        return Err("equilibrate_flux_patches requires at least one terminal potential vector".into());
    }
    let ndofs = mesh.n_dofs(order);
    for (j, v) in potentials.iter().enumerate() {
        if v.len() != ndofs {
            return Err(format!("potential vector {j} has length {}, expected {ndofs}", v.len()));
        }
        if v.iter().any(|x| !x.is_finite()) {
            return Err(format!("potential vector {j} contains NaN or infinite entries"));
        }
    }
    let eps_len = match eps {
        CellPermittivity::Scalar(s) => s.len(),
        CellPermittivity::Tensor(t) => t.len(),
    };
    if eps_len != nc {
        return Err(format!("permittivity slice length {eps_len} does not match mesh.n_cells() = {nc}"));
    }
    if let Some(cu) = c_upper_override {
        if cu.len() != n_term || cu.iter().any(|r| r.len() != n_term) {
            return Err("c_upper_override dimensions do not match number of terminals".into());
        }
    }

    // Ensure shared rayon thread pool is initialized according to NYX_THREADS
    let nt = crate::par::threads();

    // 1. Canonical sorted-vertex cell geometry & material inverses (with Rule-3 QP CurvedTet guard)
    let scells = build_sorted_cells(mesh, order, eps)?;

    // 2. Map mesh facets by sorted global vertex key, and build face -> incident cells [(c, opp_slot)]
    let mut facet_tags: BTreeMap<[u32; 3], FacetTag> = BTreeMap::new();
    for f in 0..nf {
        let key = sorted_facet_key(mesh.facet(f));
        let is_sheet = sheet_facets[f];
        let is_dir = dirichlet_facets[f] || is_sheet;
        let rob_f = match robin {
            RobinBoundary::None => None,
            RobinBoundary::FacetAlpha(a) => (a[f] > 0.0 && !is_dir).then_some(f),
            RobinBoundary::Pointwise(_) => (!is_dir).then_some(f),
        };
        let entry = facet_tags.entry(key).or_insert(FacetTag {
            is_dirichlet: false,
            is_sheet: false,
            robin_facet: None,
        });
        entry.is_dirichlet |= is_dir;
        entry.is_sheet |= is_sheet;
        if entry.is_dirichlet {
            entry.robin_facet = None;
        } else if entry.robin_facet.is_none() {
            entry.robin_facet = rob_f;
        }
    }

    let mut face_to_cells: BTreeMap<[u32; 3], Vec<(usize, usize)>> = BTreeMap::new();
    for (c, sc) in scells.iter().enumerate() {
        for m in 0..=d {
            let key = cell_face_key(d, &sc.verts, m);
            face_to_cells.entry(key).or_default().push((c, m));
        }
    }

    // Pre-verify which candidate Robin exterior faces actually have positive alpha > 0 at facet quadrature points
    let frule = facet_rule(d - 1);
    let mut active_robin_faces: BTreeMap<[u32; 3], usize> = BTreeMap::new();
    for (&key, inc) in &face_to_cells {
        if inc.len() == 1 {
            if let Some(tag) = facet_tags.get(&key) {
                if !tag.is_dirichlet {
                    if let Some(f_idx) = tag.robin_facet {
                        let (c0, m0) = inc[0];
                        let mut any_pos = false;
                        for (q, mu) in frule.lam.iter().enumerate() {
                            let (_, xp, n_out, _) = eval_cell_face_point(d, &scells[c0], m0, mu, frule.w[q], &mesh.points);
                            let alpha = eval_robin_alpha(robin, f_idx, xp, n_out);
                            if !alpha.is_finite() || alpha < 0.0 {
                                return Err(format!("non-finite or negative Robin alpha {alpha} on facet {f_idx}"));
                            }
                            if alpha > 0.0 {
                                any_pos = true;
                            }
                        }
                        if any_pos {
                            active_robin_faces.insert(key, f_idx);
                        }
                    }
                }
            }
        }
    }

    // 3. Build vertex stars and split each vertex star into connected half-stars across interior conductive sheets
    let nv = mesh.n_points();
    let mut vert_cells: Vec<Vec<(usize, usize)>> = vec![Vec::new(); nv];
    for (c, sc) in scells.iter().enumerate() {
        for slot in 0..=d {
            vert_cells[sc.verts[slot] as usize].push((c, slot));
        }
    }

    let n_face_dofs = if order == 1 { 1 } else { d };
    let n_bub_dofs = if order == 1 { 0 } else { d };

    let mut patches: Vec<VertexPatch> = Vec::with_capacity(nv);
    // For each cell c and sorted vertex slot a in 0..=d, record (patch_id, local_cell_idx_in_patch)
    let mut cell_slot_to_patch: Vec<[(usize, usize); 4]> = vec![[(usize::MAX, usize::MAX); 4]; nc];
    let mut n_sheet_splits = 0usize;

    for v in 0..nv {
        let star = &vert_cells[v];
        if star.is_empty() {
            continue;
        }
        let n_star = star.len();
        let mut visited = vec![false; n_star];
        let mut v_comps = 0usize;

        // Map global cell index -> index in `star`
        let star_pos_of = |cell_id: usize| -> Option<usize> {
            star.binary_search_by_key(&cell_id, |&(cid, _)| cid).ok()
        };

        for start_i in 0..n_star {
            if visited[start_i] {
                continue;
            }
            v_comps += 1;
            let mut comp_indices = Vec::new();
            let mut queue = VecDeque::new();
            visited[start_i] = true;
            queue.push_back(start_i);

            while let Some(si) = queue.pop_front() {
                comp_indices.push(si);
                let (c_curr, v_slot) = star[si];
                for m in 0..=d {
                    if m == v_slot {
                        // Opposite face does not contain vertex v
                        continue;
                    }
                    let key = cell_face_key(d, &scells[c_curr].verts, m);
                    let is_cut = facet_tags.get(&key).is_some_and(|t| t.is_sheet || t.is_dirichlet);
                    if is_cut {
                        continue;
                    }
                    if let Some(inc) = face_to_cells.get(&key) {
                        if inc.len() == 2 {
                            let other_c = if inc[0].0 == c_curr { inc[1].0 } else { inc[0].0 };
                            if let Some(other_si) = star_pos_of(other_c) {
                                if !visited[other_si] {
                                    visited[other_si] = true;
                                    queue.push_back(other_si);
                                }
                            }
                        }
                    }
                }
            }
            comp_indices.sort_unstable();

            let patch_id = patches.len();
            let n_pcells = comp_indices.len();
            let mut p_cells = Vec::with_capacity(n_pcells);
            let mut p_vslots = Vec::with_capacity(n_pcells);
            for (loc_ci, &si) in comp_indices.iter().enumerate() {
                let (cid, vslot) = star[si];
                p_cells.push(cid);
                p_vslots.push(vslot);
                cell_slot_to_patch[cid][vslot] = (patch_id, loc_ci);
            }

            let init_face = PatchCellFace { dof_offset: usize::MAX, sign: 1.0, robin_facet: None };
            let mut cell_faces = vec![[init_face; 4]; n_pcells];
            let mut next_dof = 0usize;
            let mut is_open = false;

            // Assign shared DOFs on non-cut interior faces containing v, and independent DOFs on Dirichlet/sheet/Robin faces containing v
            let mut assigned_interior: BTreeMap<[u32; 3], usize> = BTreeMap::new();
            for loc_ci in 0..n_pcells {
                let cid = p_cells[loc_ci];
                let vslot = p_vslots[loc_ci];
                for m in 0..=d {
                    if m == vslot {
                        // Outer face of vertex star (opposite v): sigma . n = 0
                        continue;
                    }
                    let key = cell_face_key(d, &scells[cid].verts, m);
                    let inc = &face_to_cells[&key];
                    let tag = facet_tags.get(&key);
                    let is_dir_or_sheet = tag.is_some_and(|t| t.is_dirichlet || t.is_sheet);

                    if inc.len() == 2 && !is_dir_or_sheet {
                        // Interior dielectric face within this connected half-star
                        let (c0, _) = inc[0];
                        let off = *assigned_interior.entry(key).or_insert_with(|| {
                            let o = next_dof;
                            next_dof += n_face_dofs;
                            o
                        });
                        let sign = if cid == c0 { 1.0 } else { -1.0 };
                        cell_faces[loc_ci][m] = PatchCellFace { dof_offset: off, sign, robin_facet: None };
                    } else if is_dir_or_sheet {
                        // Unconstrained Dirichlet wall or one side of an interior conductive sheet
                        let off = next_dof;
                        next_dof += n_face_dofs;
                        is_open = true;
                        cell_faces[loc_ci][m] = PatchCellFace { dof_offset: off, sign: 1.0, robin_facet: None };
                    } else if inc.len() == 1 {
                        if let Some(&rob_f) = active_robin_faces.get(&key) {
                            let off = next_dof;
                            next_dof += n_face_dofs;
                            is_open = true;
                            cell_faces[loc_ci][m] = PatchCellFace { dof_offset: off, sign: 1.0, robin_facet: Some(rob_f) };
                        }
                        // else: homogeneous Neumann wall -> dof_offset stays usize::MAX (sigma . n = 0)
                    }
                }
            }

            let mut bubble_offset = vec![usize::MAX; n_pcells];
            if n_bub_dofs > 0 {
                for bo in &mut bubble_offset {
                    *bo = next_dof;
                    next_dof += n_bub_dofs;
                }
            }

            patches.push(VertexPatch {
                cells: p_cells,
                v_slot_in_cell: p_vslots,
                cell_faces,
                bubble_offset,
                n_sigma: next_dof,
                is_open,
            });
        }

        if v_comps >= 2 {
            n_sheet_splits += 1;
        }
    }

    let n_patches = patches.len();
    let tab_aff = PrecomputedRefTab::new(d, order, false);
    let tab_curv = (d == 3 && mesh.has_curved()).then(|| PrecomputedRefTab::new(3, order, true));

    // 4. Precompute per-cell, per-vertex-slot divergence targets `g0[c][slot]` (and `g_lin[c][slot][m-1]` for order=2)
    //    so any tiny iterative algebraic residual on closed patches can be swept along a deterministic spanning forest
    //    to open patches (`Gamma_D U Gamma_inf`), preserving exact `sum_slot g_slot == 0` in every cell!
    let nq_per_cell = if order == 1 { 1 } else { d + 1 };
    // cell_div_rhs[c][(slot * nq_per_cell + beta) * n_term + j]
    let stride_cell_rhs = (d + 1) * nq_per_cell * n_term;
    let build_cell_rhs = |c: usize| -> Vec<f64> {
        let sc = &scells[c];
        let tab = if sc.curved.is_some() { tab_curv.as_ref().unwrap() } else { &tab_aff };
        let mut out = vec![0.0; stride_cell_rhs];
        let inv_dp1 = 1.0 / (d + 1) as f64;

        for q in 0..tab.w.len() {
            let lam = &tab.lam[q];
            let (grad_lam, w_vol) = if let Some(ct) = &sc.curved {
                let (gq, det_q) = ct.jacobian(lam);
                (gq, tab.w[q] * det_q.abs() / 6.0)
            } else {
                let vol = sc.det_aff.abs() / if d == 2 { 2.0 } else { 6.0 };
                (sc.grad_lam_aff, tab.w[q] * vol)
            };

            for j in 0..n_term {
                let (_, grad_v) = eval_vh_and_grad(d, order, sc, &grad_lam, tab, q, &potentials[j]);
                let neg_gv = [-grad_v[0], -grad_v[1], -grad_v[2]];
                let dh = matvec_sym_d(d, &sc.eps_fwd, &neg_gv);

                for slot in 0..=d {
                    let gpsi_dot_dh = dot_d(d, &grad_lam[slot], &dh);
                    let base = (slot * nq_per_cell) * n_term + j;
                    out[base] += w_vol * gpsi_dot_dh;
                    if order == 2 {
                        for m in 1..=d {
                            out[(slot * nq_per_cell + m) * n_term + j] += w_vol * gpsi_dot_dh * (lam[m] - inv_dp1);
                        }
                    }
                }
            }
        }
        out
    };

    const RT_CELL_CHUNK: usize = 256;
    let cell_rhs_vecs: Vec<Vec<f64>> = crate::par::chunk_map(
        nc.div_ceil(RT_CELL_CHUNK),
        nt > 1 && nc >= 256,
        |cid| {
            let c0 = cid * RT_CELL_CHUNK;
            let c1 = (c0 + RT_CELL_CHUNK).min(nc);
            (c0..c1).map(build_cell_rhs).collect::<Vec<_>>()
        },
    )
    .into_iter()
    .flatten()
    .collect();

    // Assemble patch divergence vectors `patch_g[p]` of shape `(n_pcells * nq_per_cell) x n_term`
    let mut patch_g: Vec<Vec<f64>> = patches
        .iter()
        .map(|p| {
            let mut gp = vec![0.0; p.cells.len() * nq_per_cell * n_term];
            for (loc_ci, (&cid, &vslot)) in p.cells.iter().zip(&p.v_slot_in_cell).enumerate() {
                let src = &cell_rhs_vecs[cid];
                for beta in 0..nq_per_cell {
                    for j in 0..n_term {
                        gp[(loc_ci * nq_per_cell + beta) * n_term + j] =
                            src[(vslot * nq_per_cell + beta) * n_term + j];
                    }
                }
            }
            gp
        })
        .collect();

    // Deterministic BFS spanning forest from open patches (`is_open == true`) to sweep any closed-patch algebraic residual
    {
        // Build patch adjacency `(neighbor_patch, shared_cell_id, my_loc_ci, nbr_loc_ci)` in deterministic order
        let mut patch_adj: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); n_patches];
        for c in 0..nc {
            let slots = &cell_slot_to_patch[c];
            for a in 0..=d {
                let (pa, la) = slots[a];
                for b in (a + 1)..=d {
                    let (pb, lb) = slots[b];
                    if pa != pb {
                        if patch_adj[pa].last().map(|x| x.0) != Some(pb)
                            && !patch_adj[pa].iter().any(|x| x.0 == pb)
                        {
                            patch_adj[pa].push((pb, la, lb));
                        }
                        if patch_adj[pb].last().map(|x| x.0) != Some(pa)
                            && !patch_adj[pb].iter().any(|x| x.0 == pa)
                        {
                            patch_adj[pb].push((pa, lb, la));
                        }
                    }
                }
            }
        }
        for adj in &mut patch_adj {
            adj.sort_unstable_by_key(|x| x.0);
        }

        let mut visited_p = vec![false; n_patches];
        // (child_patch, parent_patch, child_loc_ci, parent_loc_ci)
        let mut tree_edges: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(n_patches);
        let mut bfs = VecDeque::new();
        for (pid, p) in patches.iter().enumerate() {
            if p.is_open {
                visited_p[pid] = true;
                bfs.push_back(pid);
            }
        }
        while let Some(u) = bfs.pop_front() {
            for &(w, lu, lw) in &patch_adj[u] {
                if !visited_p[w] {
                    visited_p[w] = true;
                    tree_edges.push((w, u, lw, lu));
                    bfs.push_back(w);
                }
            }
        }
        // Any remaining completely closed connected component (no Dirichlet and no Robin faces anywhere)
        for root in 0..n_patches {
            if !visited_p[root] {
                visited_p[root] = true;
                let mut comp = vec![root];
                bfs.push_back(root);
                while let Some(u) = bfs.pop_front() {
                    for &(w, lu, lw) in &patch_adj[u] {
                        if !visited_p[w] {
                            visited_p[w] = true;
                            comp.push(w);
                            tree_edges.push((w, u, lw, lu));
                            bfs.push_back(w);
                        }
                    }
                }
                // Zero-sum residual cleanup on root if needed
                let _ = comp;
            }
        }

        for &(child, parent, l_child, l_parent) in tree_edges.iter().rev() {
            let n_child_cells = patches[child].cells.len();
            for j in 0..n_term {
                let mut res = 0.0;
                for ci in 0..n_child_cells {
                    res += patch_g[child][(ci * nq_per_cell) * n_term + j];
                }
                if res != 0.0 {
                    patch_g[child][(l_child * nq_per_cell) * n_term + j] -= res;
                    patch_g[parent][(l_parent * nq_per_cell) * n_term + j] += res;
                }
            }
        }
        // If any isolated closed root remained with roundoff, remove its mean
        for (pid, p) in patches.iter().enumerate() {
            if !p.is_open && !p.cells.is_empty() {
                let npc = p.cells.len();
                for j in 0..n_term {
                    let mut res = 0.0;
                    for ci in 0..npc {
                        res += patch_g[pid][(ci * nq_per_cell) * n_term + j];
                    }
                    if res != 0.0 {
                        let sub = res / npc as f64;
                        for ci in 0..npc {
                            patch_g[pid][(ci * nq_per_cell) * n_term + j] -= sub;
                        }
                    }
                }
            }
        }
    }

    // 5. Solve each vertex half-star patch independently via pre-scaled SPD Schur complement + 1 refinement step
    let nrt = rt_nloc(d, order);
    let f_lag_shapes: Vec<Vec<f64>> = frule.lam.iter().map(|mu| shape(mu, order).0).collect();

    let solve_patch = |pid: usize| -> Result<Vec<f64>, String> {
        let patch = &patches[pid];
        let npc = patch.cells.len();
        let n_sig = patch.n_sigma;
        // Return flat [loc_ci * nrt * n_term + r * n_term + j]
        let mut cell_contrib = vec![0.0; npc * nrt * n_term];
        if n_sig == 0 || npc == 0 {
            return Ok(cell_contrib);
        }

        let n_q_full = npc * nq_per_cell;
        // On closed patches, row 0 (cell 0 constant mode) is redundant because sum_{ci} row_{(ci, 0)} == 0
        // and sum_{ci} g_{(ci, 0)} == 0.
        let drop_row0 = !patch.is_open;
        let n_q = if drop_row0 { n_q_full.saturating_sub(1) } else { n_q_full };
        let map_q_row = |full_r: usize| -> Option<usize> {
            if drop_row0 {
                if full_r == 0 { None } else { Some(full_r - 1) }
            } else {
                Some(full_r)
            }
        };

        let mut m_mat = vec![0.0; n_sig * n_sig];
        let mut b_mat = vec![0.0; n_q * n_sig];
        let mut f_rhs = vec![0.0; n_sig * n_term];
        let mut g_rhs = vec![0.0; n_q * n_term];
        let gp = &patch_g[pid];
        for full_r in 0..n_q_full {
            if let Some(qr) = map_q_row(full_r) {
                for j in 0..n_term {
                    g_rhs[qr * n_term + j] = gp[full_r * n_term + j];
                }
            }
        }

        let inv_dp1 = 1.0 / (d + 1) as f64;
        let mut loc_map = vec![(usize::MAX, 0.0); nrt];
        let mut phi_vec = vec![[0.0; 3]; nrt];
        let mut phi_div = vec![0.0; nrt];

        for loc_ci in 0..npc {
            let cid = patch.cells[loc_ci];
            let vslot = patch.v_slot_in_cell[loc_ci];
            let sc = &scells[cid];
            let tab = if sc.curved.is_some() { tab_curv.as_ref().unwrap() } else { &tab_aff };

            // Build local cell RT basis -> patch active DOF map `(patch_dof, sign)`
            loc_map.fill((usize::MAX, 0.0));
            if order == 1 {
                for m in 0..=d {
                    let pf = patch.cell_faces[loc_ci][m];
                    if pf.dof_offset != usize::MAX {
                        loc_map[m] = (pf.dof_offset, pf.sign);
                    }
                }
            } else {
                let mut idx = 0;
                for m in 0..=d {
                    let pf = patch.cell_faces[loc_ci][m];
                    for _a in 0..d {
                        if pf.dof_offset != usize::MAX {
                            loc_map[idx] = (pf.dof_offset + _a, pf.sign);
                        }
                        idx += 1;
                    }
                }
                let bo = patch.bubble_offset[loc_ci];
                for b in 0..d {
                    loc_map[idx] = (bo + b, 1.0);
                    idx += 1;
                }
            }

            // Volume integrals on cell `cid`
            for q in 0..tab.w.len() {
                let lam = &tab.lam[q];
                let psi_v = lam[vslot];
                let (j_mat, det_q, grad_lam, w_vol) = if let Some(ct) = &sc.curved {
                    let jm = ct.jacobian_matrix(lam);
                    let (gq, dq) = ct.jacobian(lam);
                    (jm, dq, gq, tab.w[q] * dq.abs() / 6.0)
                } else {
                    let vol = sc.det_aff.abs() / if d == 2 { 2.0 } else { 6.0 };
                    (sc.j_aff, sc.det_aff, sc.grad_lam_aff, tab.w[q] * vol)
                };

                tab.piola_at_q(d, q, &j_mat, det_q, &mut phi_vec, &mut phi_div);

                // Assemble M_p and B_p
                for r in 0..nrt {
                    let (gr, sr) = loc_map[r];
                    if gr == usize::MAX {
                        continue;
                    }
                    let pr = [sr * phi_vec[r][0], sr * phi_vec[r][1], sr * phi_vec[r][2]];
                    let einv_pr = matvec_sym_d(d, &sc.eps_inv, &pr);
                    for s in 0..nrt {
                        let (gs, ss) = loc_map[s];
                        if gs == usize::MAX {
                            continue;
                        }
                        let ps = [ss * phi_vec[s][0], ss * phi_vec[s][1], ss * phi_vec[s][2]];
                        m_mat[gr * n_sig + gs] += w_vol * dot_d(d, &einv_pr, &ps);
                    }

                    let div_r = sr * phi_div[r];
                    if let Some(qr0) = map_q_row(loc_ci * nq_per_cell) {
                        b_mat[qr0 * n_sig + gr] += w_vol * div_r;
                    }
                    if order == 2 {
                        for m in 1..=d {
                            if let Some(qrm) = map_q_row(loc_ci * nq_per_cell + m) {
                                b_mat[qrm * n_sig + gr] += w_vol * div_r * (lam[m] - inv_dp1);
                            }
                        }
                    }
                }

                // Assemble RHS F_p = int psi_v * phi_r . eps^-1 D_h^{(j)} = -int psi_v * phi_r . grad(V_h^{(j)})!
                // Wait: notice eps^-1 * D_h^{(j)} = eps^-1 * (-eps * grad V_h^{(j)}) = -grad V_h^{(j)} IDENTICALLY!
                for j in 0..n_term {
                    let (_, grad_v) = eval_vh_and_grad(d, order, sc, &grad_lam, tab, q, &potentials[j]);
                    let target_e = [-psi_v * grad_v[0], -psi_v * grad_v[1], -psi_v * grad_v[2]];
                    for r in 0..nrt {
                        let (gr, sr) = loc_map[r];
                        if gr != usize::MAX {
                            f_rhs[gr * n_term + j] += w_vol * sr * dot_d(d, &phi_vec[r], &target_e);
                        }
                    }
                }
            }

            // Optional Robin boundary integrals on active `Gamma_inf` faces of cell `cid` containing `v`
            for m in 0..=d {
                let pf = patch.cell_faces[loc_ci][m];
                let Some(rob_f) = pf.robin_facet else { continue };
                // Find slot position of `vslot` among the `d` vertices of face `m`
                let mut v_pos_on_face = 0usize;
                let mut p = 0;
                for a in 0..=d {
                    if a != m {
                        if a == vslot {
                            v_pos_on_face = p;
                        }
                        p += 1;
                    }
                }
                // Build face DOF map on `sc` for `V_h` restriction to face `m`
                let mut face_vdofs = [0usize; 6];
                {
                    let mut fv = [0u32; 3];
                    let mut k = 0;
                    for a in 0..=d {
                        if a != m {
                            fv[k] = sc.verts[a];
                            face_vdofs[k] = fv[k] as usize;
                            k += 1;
                        }
                    }
                    if order == 2 {
                        let mut e_k = d;
                        for a in 0..d {
                            for b in (a + 1)..d {
                                face_vdofs[e_k] = nv + mesh.edge_index(fv[a], fv[b]);
                                e_k += 1;
                            }
                        }
                    }
                }
                let n_f_vdofs = if order == 1 { d } else { d + d * (d - 1) / 2 };

                for (fq, mu) in frule.lam.iter().enumerate() {
                    let (lam_f, xp, n_out, ds_q) = eval_cell_face_point(d, sc, m, mu, frule.w[fq], &mesh.points);
                    let alpha = eval_robin_alpha(robin, rob_f, xp, n_out);
                    if alpha <= 0.0 {
                        continue;
                    }
                    let inv_alpha = 1.0 / alpha;
                    let (j_mat, det_q) = if let Some(ct) = &sc.curved {
                        let jm = ct.jacobian_matrix(&lam_f);
                        let (_, dq) = ct.jacobian(&lam_f);
                        (jm, dq)
                    } else {
                        (sc.j_aff, sc.det_aff)
                    };
                    eval_rt_piola(d, order, &j_mat, det_q, &lam_f[..=d], &mut phi_vec, &mut phi_div);

                    let psi_v_f = mu[v_pos_on_face];
                    for r in 0..nrt {
                        let (gr, sr) = loc_map[r];
                        if gr == usize::MAX {
                            continue;
                        }
                        let phin_r = sr * dot_d(d, &phi_vec[r], &n_out);
                        if phin_r == 0.0 {
                            continue;
                        }
                        for s in 0..nrt {
                            let (gs, ss) = loc_map[s];
                            if gs == usize::MAX {
                                continue;
                            }
                            let phin_s = ss * dot_d(d, &phi_vec[s], &n_out);
                            m_mat[gr * n_sig + gs] += ds_q * inv_alpha * phin_r * phin_s;
                        }
                        for j in 0..n_term {
                            let mut vh_f = 0.0;
                            for a in 0..n_f_vdofs {
                                vh_f += potentials[j][face_vdofs[a]] * f_lag_shapes[fq][a];
                            }
                            f_rhs[gr * n_term + j] += ds_q * phin_r * psi_v_f * vh_f;
                        }
                    }
                }
            }
        }

        // Solve KKT system `[M_p, B_p^T; B_p, 0] [sigma; lam] = [F_p; G_p]` via SPD Schur complement
        let sigma_0 = solve_spd_multi(&m_mat, &f_rhs, n_sig, n_term)?;
        let sigma = if n_q == 0 {
            sigma_0
        } else {
            // bt_mat: n_sig x n_q
            let mut bt_mat = vec![0.0; n_sig * n_q];
            for qr in 0..n_q {
                for gr in 0..n_sig {
                    bt_mat[gr * n_q + qr] = b_mat[qr * n_sig + gr];
                }
            }
            // minv_bt = M_p^-1 * B_p^T (n_sig x n_q)
            let minv_bt = solve_spd_multi(&m_mat, &bt_mat, n_sig, n_q)?;
            // s_mat = B_p * minv_bt (n_q x n_q)
            let mut s_mat = vec![0.0; n_q * n_q];
            for r in 0..n_q {
                for c in r..n_q {
                    let mut sum = 0.0;
                    for k in 0..n_sig {
                        sum += b_mat[r * n_sig + k] * minv_bt[k * n_q + c];
                    }
                    s_mat[r * n_q + c] = sum;
                    s_mat[c * n_q + r] = sum;
                }
            }
            // rhs_s = B_p * sigma_0 - G_p
            let mut rhs_s = vec![0.0; n_q * n_term];
            for r in 0..n_q {
                for j in 0..n_term {
                    let mut bs = 0.0;
                    for k in 0..n_sig {
                        bs += b_mat[r * n_sig + k] * sigma_0[k * n_term + j];
                    }
                    rhs_s[r * n_term + j] = bs - g_rhs[r * n_term + j];
                }
            }
            let lam_sol = solve_spd_multi(&s_mat, &rhs_s, n_q, n_term)?;
            let mut sig = sigma_0;
            for k in 0..n_sig {
                for j in 0..n_term {
                    let mut corr = 0.0;
                    for r in 0..n_q {
                        corr += minv_bt[k * n_q + r] * lam_sol[r * n_term + j];
                    }
                    sig[k * n_term + j] -= corr;
                }
            }

            // One iterative refinement pass on the divergence constraint `B_p * sig = G_p`
            let mut rem_g = vec![0.0; n_q * n_term];
            for r in 0..n_q {
                for j in 0..n_term {
                    let mut bs = 0.0;
                    for k in 0..n_sig {
                        bs += b_mat[r * n_sig + k] * sig[k * n_term + j];
                    }
                    rem_g[r * n_term + j] = bs - g_rhs[r * n_term + j];
                }
            }
            let dlam = solve_spd_multi(&s_mat, &rem_g, n_q, n_term)?;
            for k in 0..n_sig {
                for j in 0..n_term {
                    let mut corr = 0.0;
                    for r in 0..n_q {
                        corr += minv_bt[k * n_q + r] * dlam[r * n_term + j];
                    }
                    sig[k * n_term + j] -= corr;
                }
            }
            sig
        };

        // Scatter `sigma` into cell-local RT coefficients for each cell of this patch
        for loc_ci in 0..npc {
            let out_base = loc_ci * nrt * n_term;
            if order == 1 {
                for m in 0..=d {
                    let pf = patch.cell_faces[loc_ci][m];
                    if pf.dof_offset != usize::MAX {
                        for j in 0..n_term {
                            cell_contrib[out_base + m * n_term + j] =
                                pf.sign * sigma[pf.dof_offset * n_term + j];
                        }
                    }
                }
            } else {
                let mut idx = 0;
                for m in 0..=d {
                    let pf = patch.cell_faces[loc_ci][m];
                    for a in 0..d {
                        if pf.dof_offset != usize::MAX {
                            let gd = pf.dof_offset + a;
                            for j in 0..n_term {
                                cell_contrib[out_base + idx * n_term + j] =
                                    pf.sign * sigma[gd * n_term + j];
                            }
                        }
                        idx += 1;
                    }
                }
                let bo = patch.bubble_offset[loc_ci];
                for b in 0..d {
                    let gd = bo + b;
                    for j in 0..n_term {
                        cell_contrib[out_base + idx * n_term + j] = sigma[gd * n_term + j];
                    }
                    idx += 1;
                }
            }
        }
        Ok(cell_contrib)
    };

    const RT_PATCH_CHUNK: usize = 16;
    let patch_chunks: Vec<Result<Vec<Vec<f64>>, String>> = crate::par::chunk_map(
        n_patches.div_ceil(RT_PATCH_CHUNK),
        nt > 1 && n_patches >= 64,
        |cid| {
            let p0 = cid * RT_PATCH_CHUNK;
            let p1 = (p0 + RT_PATCH_CHUNK).min(n_patches);
            (p0..p1).map(solve_patch).collect::<Result<Vec<_>, _>>()
        },
    );
    let mut patch_sols: Vec<Vec<f64>> = Vec::with_capacity(n_patches);
    for chunk in patch_chunks {
        patch_sols.extend(chunk?);
    }

    // 6. Deterministically accumulate global `D_tilde^{(j)} = sum_p sigma_p^{(j)}` in ascending patch order `0..n_patches`
    let mut cell_rt_coeffs = vec![vec![0.0; nc * nrt]; n_term];
    for (pid, patch) in patches.iter().enumerate() {
        let sol = &patch_sols[pid];
        for (loc_ci, &cid) in patch.cells.iter().enumerate() {
            let src_base = loc_ci * nrt * n_term;
            let dst_base = cid * nrt;
            for r in 0..nrt {
                for j in 0..n_term {
                    cell_rt_coeffs[j][dst_base + r] += sol[src_base + r * n_term + j];
                }
            }
        }
    }

    // 7. Integrate `W^flux`, `Q`, `C^up` (volume), `||div D_tilde||^2`, and per-cell Prager-Synge `eta^2_{T,PS,j}`
    struct CellQuadSummary {
        w_flux: Vec<f64>,
        q_pair: Vec<f64>,
        c_up: Vec<f64>,
        eta2: Vec<f64>,
        div2: Vec<f64>,
    }

    let eval_cell_integrals = |c: usize| -> CellQuadSummary {
        let sc = &scells[c];
        let tab = if sc.curved.is_some() { tab_curv.as_ref().unwrap() } else { &tab_aff };
        let mut w_flux = vec![0.0; n_term * n_term];
        let mut q_pair = vec![0.0; n_term * n_term];
        let mut c_up = vec![0.0; n_term * n_term];
        let mut eta2 = vec![0.0; n_term];
        let mut div2 = vec![0.0; n_term];

        let mut phi_vec = vec![[0.0; 3]; nrt];
        let mut phi_div = vec![0.0; nrt];
        let mut d_tilde = vec![[0.0; 3]; n_term];
        let mut grad_v_all = vec![[0.0; 3]; n_term];
        let mut dh_all = vec![[0.0; 3]; n_term];
        let c_base = c * nrt;

        for q in 0..tab.w.len() {
            let lam = &tab.lam[q];
            let (j_mat, det_q, grad_lam, w_vol) = if let Some(ct) = &sc.curved {
                let jm = ct.jacobian_matrix(lam);
                let (gq, dq) = ct.jacobian(lam);
                (jm, dq, gq, tab.w[q] * dq.abs() / 6.0)
            } else {
                let vol = sc.det_aff.abs() / if d == 2 { 2.0 } else { 6.0 };
                (sc.j_aff, sc.det_aff, sc.grad_lam_aff, tab.w[q] * vol)
            };

            eval_rt_piola(d, order, &j_mat, det_q, &lam[..=d], &mut phi_vec, &mut phi_div);

            for j in 0..n_term {
                let coeffs = &cell_rt_coeffs[j][c_base..c_base + nrt];
                let mut dt = [0.0; 3];
                let mut dv = 0.0;
                for r in 0..nrt {
                    let cr = coeffs[r];
                    dt[0] += cr * phi_vec[r][0];
                    dt[1] += cr * phi_vec[r][1];
                    dt[2] += cr * phi_vec[r][2];
                    dv += cr * phi_div[r];
                }
                d_tilde[j] = dt;
                div2[j] += w_vol * dv * dv;

                let (_, gv) = eval_vh_and_grad(d, order, sc, &grad_lam, tab, q, &potentials[j]);
                grad_v_all[j] = gv;
                let neg_gv = [-gv[0], -gv[1], -gv[2]];
                let dh = matvec_sym_d(d, &sc.eps_fwd, &neg_gv);
                dh_all[j] = dh;

                let diff = [dt[0] - dh[0], dt[1] - dh[1], dt[2] - dh[2]];
                let einv_diff = matvec_sym_d(d, &sc.eps_inv, &diff);
                eta2[j] += w_vol * dot_d(d, &diff, &einv_diff).max(0.0);
            }

            for i in 0..n_term {
                let einv_di = matvec_sym_d(d, &sc.eps_inv, &d_tilde[i]);
                for j in i..n_term {
                    let wij = w_vol * dot_d(d, &einv_di, &d_tilde[j]);
                    w_flux[i * n_term + j] += wij;
                    if i != j {
                        w_flux[j * n_term + i] += wij;
                    }
                    // c_up = int grad(V_i) . eps grad(V_j) = int (-grad V_i) . D_h^(j)
                    let cij = -w_vol * dot_d(d, &grad_v_all[i], &dh_all[j]);
                    c_up[i * n_term + j] += cij;
                    if i != j {
                        c_up[j * n_term + i] += cij;
                    }
                }
                for j in 0..n_term {
                    q_pair[i * n_term + j] += -w_vol * dot_d(d, &grad_v_all[i], &d_tilde[j]);
                }
            }
        }
        CellQuadSummary { w_flux, q_pair, c_up, eta2, div2 }
    };

    let cell_summaries: Vec<CellQuadSummary> = crate::par::chunk_map(
        nc.div_ceil(RT_CELL_CHUNK),
        nt > 1 && nc >= 256,
        |cid| {
            let c0 = cid * RT_CELL_CHUNK;
            let c1 = (c0 + RT_CELL_CHUNK).min(nc);
            (c0..c1).map(eval_cell_integrals).collect::<Vec<_>>()
        },
    )
    .into_iter()
    .flatten()
    .collect();

    let mut w_flux_matrix = vec![vec![0.0; n_term]; n_term];
    let mut q_pairing_matrix = vec![vec![0.0; n_term]; n_term];
    let mut c_upper_matrix = vec![vec![0.0; n_term]; n_term];
    let mut terminal_cell_eta2 = vec![vec![0.0; nc]; n_term];
    let mut div_l2_sq = vec![0.0; n_term];

    // Accumulate in exact cell order `0..nc`
    for (c, s) in cell_summaries.iter().enumerate() {
        for i in 0..n_term {
            terminal_cell_eta2[i][c] = s.eta2[i];
            div_l2_sq[i] += s.div2[i];
            for j in 0..n_term {
                w_flux_matrix[i][j] += s.w_flux[i * n_term + j];
                q_pairing_matrix[i][j] += s.q_pair[i * n_term + j];
                c_upper_matrix[i][j] += s.c_up[i * n_term + j];
            }
        }
    }

    // Add Robin `Gamma_inf` boundary terms to `W^flux`, `Q`, and `C^up` (if not overridden) in deterministic sorted face order
    if !active_robin_faces.is_empty() {
        let mut phi_vec = vec![[0.0; 3]; nrt];
        let mut phi_div = vec![0.0; nrt];
        let n_f_vdofs = if order == 1 { d } else { d + d * (d - 1) / 2 };

        for (&key, &rob_f) in &active_robin_faces {
            let (c0, m0) = face_to_cells[&key][0];
            let sc = &scells[c0];
            let c_base = c0 * nrt;

            let mut face_vdofs = [0usize; 6];
            {
                let mut fv = [0u32; 3];
                let mut k = 0;
                for a in 0..=d {
                    if a != m0 {
                        fv[k] = sc.verts[a];
                        face_vdofs[k] = fv[k] as usize;
                        k += 1;
                    }
                }
                if order == 2 {
                    let mut e_k = d;
                    for a in 0..d {
                        for b in (a + 1)..d {
                            face_vdofs[e_k] = nv + mesh.edge_index(fv[a], fv[b]);
                            e_k += 1;
                        }
                    }
                }
            }

            for (fq, mu) in frule.lam.iter().enumerate() {
                let (lam_f, xp, n_out, ds_q) = eval_cell_face_point(d, sc, m0, mu, frule.w[fq], &mesh.points);
                let alpha = eval_robin_alpha(robin, rob_f, xp, n_out);
                if alpha <= 0.0 {
                    continue;
                }
                let inv_alpha = 1.0 / alpha;
                let (j_mat, det_q) = if let Some(ct) = &sc.curved {
                    let jm = ct.jacobian_matrix(&lam_f);
                    let (_, dq) = ct.jacobian(&lam_f);
                    (jm, dq)
                } else {
                    (sc.j_aff, sc.det_aff)
                };
                eval_rt_piola(d, order, &j_mat, det_q, &lam_f[..=d], &mut phi_vec, &mut phi_div);

                let mut dn_vals = vec![0.0; n_term];
                let mut vh_vals = vec![0.0; n_term];
                for j in 0..n_term {
                    let coeffs = &cell_rt_coeffs[j][c_base..c_base + nrt];
                    let mut dt = [0.0; 3];
                    for r in 0..nrt {
                        let cr = coeffs[r];
                        dt[0] += cr * phi_vec[r][0];
                        dt[1] += cr * phi_vec[r][1];
                        dt[2] += cr * phi_vec[r][2];
                    }
                    dn_vals[j] = dot_d(d, &dt, &n_out);

                    let mut vh_f = 0.0;
                    for a in 0..n_f_vdofs {
                        vh_f += potentials[j][face_vdofs[a]] * f_lag_shapes[fq][a];
                    }
                    vh_vals[j] = vh_f;
                }

                for i in 0..n_term {
                    for j in 0..n_term {
                        w_flux_matrix[i][j] += ds_q * inv_alpha * dn_vals[i] * dn_vals[j];
                        q_pairing_matrix[i][j] += ds_q * vh_vals[i] * dn_vals[j];
                        c_upper_matrix[i][j] += ds_q * alpha * vh_vals[i] * vh_vals[j];
                    }
                }
            }
        }
    }

    if let Some(cu) = c_upper_override {
        for i in 0..n_term {
            for j in i..n_term {
                let s = 0.5 * (cu[i][j] + cu[j][i]);
                c_upper_matrix[i][j] = s;
                c_upper_matrix[j][i] = s;
            }
        }
    } else {
        for i in 0..n_term {
            for j in (i + 1)..n_term {
                let s = 0.5 * (c_upper_matrix[i][j] + c_upper_matrix[j][i]);
                c_upper_matrix[i][j] = s;
                c_upper_matrix[j][i] = s;
            }
        }
    }

    for i in 0..n_term {
        for j in (i + 1)..n_term {
            let s = 0.5 * (w_flux_matrix[i][j] + w_flux_matrix[j][i]);
            w_flux_matrix[i][j] = s;
            w_flux_matrix[j][i] = s;
        }
    }

    // 8. Guaranteed symmetric PSD Loewner lower bound `C^low = Q * (W^flux)^+ * Q^T` via eigenthresholded `dense_eigh`
    let mut flat_w = Vec::with_capacity(n_term * n_term);
    for row in &w_flux_matrix {
        flat_w.extend_from_slice(row);
    }
    let (eigs_w, vecs_w) = dense_eigh(&flat_w, n_term)?;
    let lam_max_w = eigs_w.iter().copied().fold(0.0_f64, f64::max);
    let mut c_lower_matrix = vec![vec![0.0; n_term]; n_term];
    if lam_max_w > 0.0 {
        let thresh = 1.0e-13 * lam_max_w;
        for (k, &lam) in eigs_w.iter().enumerate() {
            if lam > thresh {
                let vk = &vecs_w[k];
                // u_i = (Q * vk)_i
                let u: Vec<f64> = (0..n_term)
                    .map(|i| (0..n_term).map(|j| q_pairing_matrix[i][j] * vk[j]).sum())
                    .collect();
                let inv_l = 1.0 / lam;
                for i in 0..n_term {
                    for j in i..n_term {
                        let add = inv_l * u[i] * u[j];
                        c_lower_matrix[i][j] += add;
                        if i != j {
                            c_lower_matrix[j][i] += add;
                        }
                    }
                }
            }
        }
    }

    let mut c_mid_matrix = vec![vec![0.0; n_term]; n_term];
    let mut c_error_bound_matrix = vec![vec![0.0; n_term]; n_term];
    let diag_gaps: Vec<f64> = (0..n_term)
        .map(|i| (c_upper_matrix[i][i] - c_lower_matrix[i][i]).max(0.0))
        .collect();
    for i in 0..n_term {
        for j in i..n_term {
            let mid = 0.5 * (c_upper_matrix[i][j] + c_lower_matrix[i][j]);
            let err = 0.5 * (diag_gaps[i] * diag_gaps[j]).sqrt();
            c_mid_matrix[i][j] = mid;
            c_mid_matrix[j][i] = mid;
            c_error_bound_matrix[i][j] = err;
            c_error_bound_matrix[j][i] = err;
        }
    }

    let max_div_l2 = div_l2_sq.iter().map(|x| x.max(0.0).sqrt()).fold(0.0_f64, f64::max);
    let cell_eta2 = goal_weighted_cell_eta2(&terminal_cell_eta2, &c_mid_matrix, 1.0e-3);

    Ok(EquilibratedFluxResult {
        order,
        rt_order: order - 1,
        dim: d,
        n_patches,
        n_sheet_splits,
        max_div_l2,
        cell_rt_coeffs,
        w_flux_matrix,
        q_pairing_matrix,
        c_upper_matrix,
        c_lower_matrix,
        c_mid_matrix,
        c_error_bound_matrix,
        terminal_cell_eta2,
        cell_eta2,
    })
}

impl EquilibratedFluxResult {
    /// Evaluates physical `(D_tilde^(terminal)(x), div D_tilde^(terminal)(x))` on `cell` at barycentric coordinates
    /// `sorted_lam` (`len == dim + 1`, ordered with respect to `cell`'s ascending global vertex indices).
    pub fn eval_cell_flux(
        &self,
        mesh: &SimplexMesh,
        terminal: usize,
        cell: usize,
        sorted_lam: &[f64],
    ) -> Result<([f64; 3], f64), String> {
        let d = mesh.dim;
        assert_eq!(sorted_lam.len(), d + 1);
        let raw = mesh.cell(cell);
        let mut verts = [u32::MAX; 4];
        for (i, &v) in raw.iter().enumerate() {
            verts[i] = v;
        }
        verts[..=d].sort_unstable();
        let mut pts = [[0.0; 3]; 4];
        for i in 0..=d {
            let p = mesh.point(verts[i] as usize);
            pts[i][..d].copy_from_slice(&p[..d]);
        }
        let (j_mat, det_j) = if d == 3 && mesh.cell_is_curved(cell) {
            let mut x10 = [[0.0; 3]; 10];
            x10[..4].copy_from_slice(&pts);
            let mut k = 4;
            for i in 0..4 {
                for j in (i + 1)..4 {
                    x10[k] = mesh.edge_node(mesh.edge_index(verts[i], verts[j])).unwrap_or_else(|| {
                        [
                            0.5 * (pts[i][0] + pts[j][0]),
                            0.5 * (pts[i][1] + pts[j][1]),
                            0.5 * (pts[i][2] + pts[j][2]),
                        ]
                    });
                    k += 1;
                }
            }
            let ct = CurvedTet::new(x10);
            let lam4 = [sorted_lam[0], sorted_lam[1], sorted_lam[2], sorted_lam[3]];
            (ct.jacobian_matrix(&lam4), ct.jacobian(&lam4).1)
        } else {
            let (ja, da, _) = affine_jacobian_and_grads(d, &pts)?;
            (ja, da)
        };

        let nrt = rt_nloc(d, self.order);
        let mut phi_vec = vec![[0.0; 3]; nrt];
        let mut phi_div = vec![0.0; nrt];
        eval_rt_piola(d, self.order, &j_mat, det_j, sorted_lam, &mut phi_vec, &mut phi_div);
        let coeffs = &self.cell_rt_coeffs[terminal][cell * nrt..(cell + 1) * nrt];
        let mut dt = [0.0; 3];
        let mut dv = 0.0;
        for r in 0..nrt {
            let cr = coeffs[r];
            dt[0] += cr * phi_vec[r][0];
            dt[1] += cr * phi_vec[r][1];
            dt[2] += cr * phi_vec[r][2];
            dv += cr * phi_div[r];
        }
        Ok((dt, dv))
    }

    /// Computes `(sum_{F not in sheet} ||[[D_tilde^(terminal) . n_F]]||^2_{L^2(F)})^{1/2}` over all interior non-sheet faces,
    /// verifying `H(div, Omega \ Gamma_sheet)` normal-flux continuity across dielectric interfaces.
    pub fn interior_non_sheet_normal_jump_l2(
        &self,
        mesh: &SimplexMesh,
        sheet_facets: &[bool],
        terminal: usize,
    ) -> Result<f64, String> {
        let d = mesh.dim;
        let dummy_eps = vec![1.0; mesh.n_cells()];
        let scells = build_sorted_cells(mesh, self.order, CellPermittivity::Scalar(&dummy_eps))?;
        let mut sheet_keys: std::collections::BTreeSet<[u32; 3]> = std::collections::BTreeSet::new();
        for (f, &is_sh) in sheet_facets.iter().enumerate().take(mesh.n_facets()) {
            if is_sh {
                sheet_keys.insert(sorted_facet_key(mesh.facet(f)));
            }
        }
        let mut face_to_cells: BTreeMap<[u32; 3], Vec<(usize, usize)>> = BTreeMap::new();
        for (c, sc) in scells.iter().enumerate() {
            for m in 0..=d {
                face_to_cells.entry(cell_face_key(d, &sc.verts, m)).or_default().push((c, m));
            }
        }

        let frule = facet_rule(d - 1);
        let nrt = rt_nloc(d, self.order);
        let mut phi0 = vec![[0.0; 3]; nrt];
        let mut div0 = vec![0.0; nrt];
        let mut phi1 = vec![[0.0; 3]; nrt];
        let mut div1 = vec![0.0; nrt];
        let mut jump_sq = 0.0;

        for (key, inc) in &face_to_cells {
            if inc.len() != 2 || sheet_keys.contains(key) {
                continue;
            }
            let (c0, m0) = inc[0];
            let (c1, m1) = inc[1];
            let sc0 = &scells[c0];
            let sc1 = &scells[c1];
            for (fq, mu) in frule.lam.iter().enumerate() {
                let (lam0, _, n0, ds_q) = eval_cell_face_point(d, sc0, m0, mu, frule.w[fq], &mesh.points);
                let (lam1, _, _, _) = eval_cell_face_point(d, sc1, m1, mu, frule.w[fq], &mesh.points);
                let (j0, det0) = if let Some(ct) = &sc0.curved {
                    (ct.jacobian_matrix(&lam0), ct.jacobian(&lam0).1)
                } else {
                    (sc0.j_aff, sc0.det_aff)
                };
                let (j1, det1) = if let Some(ct) = &sc1.curved {
                    (ct.jacobian_matrix(&lam1), ct.jacobian(&lam1).1)
                } else {
                    (sc1.j_aff, sc1.det_aff)
                };
                eval_rt_piola(d, self.order, &j0, det0, &lam0[..=d], &mut phi0, &mut div0);
                eval_rt_piola(d, self.order, &j1, det1, &lam1[..=d], &mut phi1, &mut div1);

                let cf0 = &self.cell_rt_coeffs[terminal][c0 * nrt..(c0 + 1) * nrt];
                let cf1 = &self.cell_rt_coeffs[terminal][c1 * nrt..(c1 + 1) * nrt];
                let mut d0 = [0.0; 3];
                let mut d1 = [0.0; 3];
                for r in 0..nrt {
                    for k in 0..d {
                        d0[k] += cf0[r] * phi0[r][k];
                        d1[k] += cf1[r] * phi1[r][k];
                    }
                }
                let jn = dot_d(d, &d0, &n0) - dot_d(d, &d1, &n0);
                jump_sq += ds_q * jn * jn;
            }
        }
        Ok(jump_sq.sqrt())
    }
}
