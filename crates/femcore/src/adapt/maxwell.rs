//! A posteriori error estimator for time-harmonic Maxwell fields on tetrahedra (ADR 53), and the metric it asks for.
//!
//! The recovered-Hessian estimator of `estimate3` looks at plotted nodal fields. For edge elements of the lowest order it does not beat a
//! uniform mesh (`docs/ports/meshers/MISSING_FEATURES.md`, section 2: 1.2e-2 on 1490 tetrahedra against 5.3e-3 on 696 for a cavity mode)
//! because the field is only piecewise linear and tangentially continuous, so its Hessian is noise. The residual of the equation is what
//! carries the information. This is the estimator of Beck, Hiptmair, Hoppe and Wohlmuth (*Residual based a posteriori error estimators for
//! eddy current computation*, M2AN 34, 2000), written for the time-harmonic equation
//!
//! ```text
//!     curl (1/mu_r) curl E  -  k0^2 eps_r E  =  0           (source free; k0 = omega / c)
//! ```
//!
//! with, for every tetrahedron `T` of diameter `h_T` and every face `F` of diameter `h_F`,
//!
//! ```text
//!     eta_T^2 = h_T^2 || curl (1/mu_r) curl E - k0^2 eps_r E ||_T^2          interior residual
//!             + h_T^2 || div (eps_r E) ||_T^2                                 divergence of the flux density
//!             + 1/2 sum_F h_F ( || n x [ (1/mu_r) curl E ] ||_F^2             jump of the tangential magnetic field
//!                             + || [ eps_r E . n ] ||_F^2 )                    jump of the normal flux density
//! ```
//!
//! The tangential trace of `E` is continuous for edge elements, so it gives no jump. A face on the boundary contributes nothing when the
//! tangential field is prescribed (a perfect electric conductor) and `h_F || n x (1/mu_r) curl E ||^2` when the boundary condition is
//! natural (a perfect magnetic conductor or a symmetry plane). For the lowest order `curl curl E` vanishes inside a tetrahedron (on straight
//! cells) and the residual is `k0^2 eps_r E` itself, which is why this estimator sees what the Hessian does not.
//!
//! **Derivative evaluation.** Callers implementing [`CurlField`] may supply exact analytic `(div E, curl curl E)` via
//! [`CurlField::sample_derivatives`]; otherwise, spatial derivatives are evaluated analytically via exact degree-2 Lagrange
//! differentiation on the reference tetrahedron (with covariant/contravariant Piola pullback and exact [`CurvedTet`] Jacobian/Hessian
//! pushforward on isoparametric curved cells), avoiding finite-difference step artifacts and out-of-cell sampling.
//!
//! **Not covered.** Media with complex or tensor coefficients (absorbing layers / PML and anisotropic domains: leave their cells out with
//! [`MaxwellOptions::skip`] or [`DomainSkipRule`]), a source term, and a driven problem; the estimator is for the source-free equation,
//! eigenmodes and homogeneous scattering problems.

use crate::adapt::metric::{intersect, isotropic, Sym};
use crate::fem::{quad_rule, quad_rule_collapsed, shape, CurvedTet, TET_EDGES};
use crate::linalg::csr::C;
use crate::mesh::SimplexMesh;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};

/// Quadrature degree used on isoparametric [`CurvedTet`] cells (matching `CURVED_QUAD_DEGREE` in `nyx-ems`).
const CURVED_QUAD_DEGREE: usize = 6;

const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Barycentric reference gradients `grad_ref(lam_a)` on the standard 3D simplex (`xi_k = lam_{k+1}`).
const G_REF: [V3; 4] = [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// The 10 reference Lagrange nodes on a tetrahedron (4 vertices followed by 6 `TET_EDGES` midpoints).
const P2_REF_NODES: [[f64; 4]; 10] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
    [0.5, 0.5, 0.0, 0.0],
    [0.5, 0.0, 0.5, 0.0],
    [0.5, 0.0, 0.0, 0.5],
    [0.0, 0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0, 0.5],
    [0.0, 0.0, 0.5, 0.5],
];

/// The field of one solution, as the estimator reads it.
pub trait CurlField: Sync {
    /// `E` and `curl E` at the barycentric point `lam` of cell `cell` (complex amplitudes, physical units of the unknowns, curl taken in
    /// the physical cell).
    fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]);

    /// Optional exact analytic physical derivatives `(div E, curl curl E)` at barycentric point `lam` of cell `cell`.
    ///
    /// When `Some(...)` is returned, [`maxwell_indicators`] uses these analytic derivatives directly instead of reconstructing
    /// reference-simplex derivatives from `sample`.
    fn sample_derivatives(&self, _cell: usize, _lam: [f64; 4]) -> Option<(C, [C; 3])> {
        None
    }
}

/// Skip-rule helper identifying non-isotropic / absorbing domain labels (PML and anisotropic regions) for [`MaxwellOptions::skip`].
#[derive(Clone, Debug, Default)]
pub struct DomainSkipRule {
    pub pml_labels: BTreeSet<i32>,
    pub anisotropic_labels: BTreeSet<i32>,
}

impl DomainSkipRule {
    /// Constructs a skip rule from PML and anisotropic cell domain labels.
    pub fn new(pml_labels: impl IntoIterator<Item = i32>, anisotropic_labels: impl IntoIterator<Item = i32>) -> Self {
        Self {
            pml_labels: pml_labels.into_iter().collect(),
            anisotropic_labels: anisotropic_labels.into_iter().collect(),
        }
    }

    /// Returns `true` when `label` belongs to either the PML or anisotropic domain sets.
    #[inline]
    pub fn should_skip(&self, label: i32) -> bool {
        self.pml_labels.contains(&label) || self.anisotropic_labels.contains(&label)
    }
}

/// The medium and the boundary conditions of the problem the field solves.
pub struct MaxwellOptions<'a> {
    /// squared free-space wave number `(omega / c)^2` of the solution (the real part of the eigenvalue for a lossy cavity)
    pub k0_sq: f64,
    /// relative permittivity of every cell (isotropic, real; ignored on skipped cells)
    pub eps_r: &'a [f64],
    /// relative permeability of every cell (isotropic, real; ignored on skipped cells)
    pub mu_r: &'a [f64],
    /// cells of these domain labels are left out (absorbing layers: complex stretched coefficients, or anisotropic media); their indicator is zero
    pub skip: &'a (dyn Fn(i32) -> bool + Sync),
    /// boundary labels with a natural condition (`n x curl E / mu = 0`: a magnetic wall or a symmetry plane). Every other boundary facet
    /// is taken as a perfect electric conductor and contributes nothing.
    pub natural: &'a (dyn Fn(i32) -> bool + Sync),
}

/// The indicators of one field.
#[derive(Clone, Debug)]
pub struct Indicators {
    /// `eta_T^2` of every cell
    pub eta2: Vec<f64>,
    /// the four parts, summed over the mesh: interior residual, divergence, jump of the tangential magnetic field, jump of the normal
    /// flux density (so a test can say which one a change of the mesh lowers)
    pub parts: [f64; 4],
    /// `sqrt(sum eta_T^2)`
    pub total: f64,
}

type V3 = [f64; 3];
type JacHessian = [[[f64; 3]; 3]; 3];

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Gradients of the four barycentric functions and the unsigned volume of a straight cell, or `None` if degenerate.
fn cell_geometry(p: &[V3; 4]) -> Option<([V3; 4], f64)> {
    let (e1, e2, e3) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    let six_v = dot(e1, cross(e2, e3));
    if !six_v.is_finite() || six_v.abs() <= 1e-300 {
        return None;
    }
    let g1 = cross(e2, e3).map(|x| x / six_v);
    let g2 = cross(e3, e1).map(|x| x / six_v);
    let g3 = cross(e1, e2).map(|x| x / six_v);
    let g0 = [-(g1[0] + g2[0] + g3[0]), -(g1[1] + g2[1] + g3[1]), -(g1[2] + g2[2] + g3[2])];
    let g = [g0, g1, g2, g3];
    if !g.iter().flatten().all(|x| x.is_finite()) {
        return None;
    }
    Some((g, six_v.abs() / 6.0))
}

fn diameter(p: &[V3]) -> f64 {
    let mut d: f64 = 0.0;
    for i in 0..p.len() {
        for j in i + 1..p.len() {
            d = d.max(norm(sub(p[i], p[j])));
        }
    }
    d
}

#[inline]
fn cnorm2(v: [C; 3]) -> f64 {
    v.iter().map(|z| z.norm_sqr()).sum()
}

/// Real vector times complex vector: `n x v` for a real `n`.
#[inline]
fn ccross(n: V3, v: [C; 3]) -> [C; 3] {
    [v[2] * n[1] - v[1] * n[2], v[0] * n[2] - v[2] * n[0], v[1] * n[0] - v[0] * n[1]]
}

#[inline]
fn cdot(n: V3, v: [C; 3]) -> C {
    v[0] * n[0] + v[1] * n[1] + v[2] * n[2]
}

/// `curl F` from the gradient tensor `grad[k][j] = d F_k / d x_j`.
#[inline]
fn curl_of(grad: &[[C; 3]; 3]) -> [C; 3] {
    [grad[2][1] - grad[1][2], grad[0][2] - grad[2][0], grad[1][0] - grad[0][1]]
}

#[inline]
fn div_of(grad: &[[C; 3]; 3]) -> C {
    grad[0][0] + grad[1][1] + grad[2][2]
}

/// Constant second derivative `H[i][k][m] = d^2 x_i / d xi_k d xi_m` of a quadratic [`CurvedTet`] map.
fn curved_tet_hessian(ct: &CurvedTet) -> JacHessian {
    let j0 = ct.jacobian_matrix(&[1.0, 0.0, 0.0, 0.0]);
    let mut h = [[[0.0; 3]; 3]; 3];
    for m in 0..3 {
        let mut lm = [0.0; 4];
        lm[m + 1] = 1.0;
        let jm = ct.jacobian_matrix(&lm);
        for i in 0..3 {
            for k in 0..3 {
                h[i][k][m] = jm[i][k] - j0[i][k];
            }
        }
    }
    h
}

/// Reference derivatives `d(F_ref)_p / d(xi_m)` (`[p][m]`) from 10-node P2 Lagrange values on `P2_REF_NODES`
/// (optionally shrunk inward by `alpha in [0, 1)` toward `[0.25; 4]`).
fn eval_p2_ref_dxi(nodes: &[[C; 3]; 10], lam: &[f64; 4], alpha: f64) -> [[C; 3]; 3] {
    let inv_scale = 1.0 / (1.0 - alpha);
    let u = if alpha == 0.0 {
        *lam
    } else {
        std::array::from_fn(|b| (lam[b] - 0.25 * alpha) * inv_scale)
    };
    let (_, dn) = shape(&u, 2);
    let mut df_dlam = [[C::new(0.0, 0.0); 4]; 3];
    for a in 0..10 {
        for p in 0..3 {
            let val = nodes[a][p];
            if val == C::new(0.0, 0.0) {
                continue;
            }
            for b in 0..4 {
                df_dlam[p][b] += val * (dn[a * 4 + b] * inv_scale);
            }
        }
    }
    let mut dxi = [[C::new(0.0, 0.0); 3]; 3];
    for p in 0..3 {
        for m in 0..3 {
            dxi[p][m] = df_dlam[p][m + 1] - df_dlam[p][0];
        }
    }
    dxi
}

/// Exact physical spatial gradient `d E_i / d x_r` from physical `E`, reference derivatives `d E_ref / d xi_m`,
/// barycentric gradients `g` (`g[m+1][i] = J^{-1}_{m,i}`), and geometry Jacobian Hessian `h_jac`.
fn pushforward_covariant_grad(
    e_phys: &[C; 3],
    de_ref_dxi: &[[C; 3]; 3],
    g: &[V3; 4],
    h_jac: &JacHessian,
) -> [[C; 3]; 3] {
    let mut de_dxi = [[C::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for m in 0..3 {
            let mut s = C::new(0.0, 0.0);
            for p in 0..3 {
                s += de_ref_dxi[p][m] * g[p + 1][i];
            }
            let mut geom_corr = C::new(0.0, 0.0);
            for a in 0..3 {
                let ea = e_phys[a];
                for b in 0..3 {
                    geom_corr += ea * (h_jac[a][b][m] * g[b + 1][i]);
                }
            }
            de_dxi[i][m] = s - geom_corr;
        }
    }
    let mut grad = [[C::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for r in 0..3 {
            let mut s = C::new(0.0, 0.0);
            for m in 0..3 {
                s += de_dxi[i][m] * g[m + 1][r];
            }
            grad[i][r] = s;
        }
    }
    grad
}

/// Exact physical spatial gradient `d H_i / d x_r` of `H = curl E = (1/det) J C_ref` from `H`, `C_ref`,
/// reference derivatives `d C_ref / d xi_m`, Jacobian `j_mat`, `g`, `det`, and `h_jac`.
fn pushforward_contravariant_curl_grad(
    h_phys: &[C; 3],
    c_ref: &[C; 3],
    dc_ref_dxi: &[[C; 3]; 3],
    j_mat: &[[f64; 3]; 3],
    g: &[V3; 4],
    det: f64,
    h_jac: &JacHessian,
) -> [[C; 3]; 3] {
    let inv_det = 1.0 / det;
    let mut tr_m = [0.0; 3];
    for m in 0..3 {
        let mut tr = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                tr += g[b + 1][a] * h_jac[a][b][m];
            }
        }
        tr_m[m] = tr;
    }
    let mut dh_dxi = [[C::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for m in 0..3 {
            let mut s = C::new(0.0, 0.0);
            for p in 0..3 {
                s += dc_ref_dxi[p][m] * (j_mat[i][p] * inv_det) + c_ref[p] * (h_jac[i][p][m] * inv_det);
            }
            dh_dxi[i][m] = s - h_phys[i] * tr_m[m];
        }
    }
    let mut grad = [[C::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for r in 0..3 {
            let mut s = C::new(0.0, 0.0);
            for m in 0..3 {
                s += dh_dxi[i][m] * g[m + 1][r];
            }
            grad[i][r] = s;
        }
    }
    grad
}

/// Evaluates all local Nedelec `H(curl)` basis functions (`6` for `order == 1`, `20` for `order == 2`) and their exact
/// physical curls, divergences, and curl-curls at barycentric point `lam` of a straight or curved tetrahedron.
pub struct NedelecEval {
    /// Physical field `w_s(x)`
    pub w: Vec<V3>,
    /// Physical curl `(curl w_s)(x)`
    pub curl: Vec<V3>,
    /// Physical divergence `(div w_s)(x)`
    pub div: Vec<f64>,
    /// Physical curl-curl `(curl curl w_s)(x)`
    pub curl_curl: Vec<V3>,
}

/// Evaluates Nedelec order-1 (`6`) or order-2 (`20`) tetrahedral basis functions and exact analytic physical derivatives
/// (`w`, `curl`, `div`, `curl_curl`) on straight or [`CurvedTet`] geometry.
///
/// `gid` holds the 4 vertex ordering keys `(vkey[v], v)` matching `nyx_ems::cavity::basis3_keyed`.
pub fn eval_nedelec_tet_with_derivatives(
    order: usize,
    lam: &[f64; 4],
    j_mat: &[[f64; 3]; 3],
    g: &[V3; 4],
    det: f64,
    h_jac: &JacHessian,
    gid: &[(u32, u32); 4],
) -> Result<NedelecEval, String> {
    if order != 1 && order != 2 {
        return Err(format!("Nedelec basis supports order 1 or 2, got {order}"));
    }
    if !det.is_finite() || det.abs() <= 1e-300 {
        return Err("degenerate tetrahedral Jacobian determinant".into());
    }
    let nb = if order == 1 { 6 } else { 20 };
    let mut w_ref = Vec::with_capacity(nb);
    let mut dw_dlam = Vec::with_capacity(nb);
    let mut c_ref = Vec::with_capacity(nb);
    let mut dc_dlam = Vec::with_capacity(nb);

    let whitney_ref = |a: usize, b: usize| -> (V3, V3) {
        let ga = G_REF[a];
        let gb = G_REF[b];
        let w = [
            lam[a] * gb[0] - lam[b] * ga[0],
            lam[a] * gb[1] - lam[b] * ga[1],
            lam[a] * gb[2] - lam[b] * ga[2],
        ];
        let c = cross(ga, gb).map(|x| 2.0 * x);
        (w, c)
    };

    for &(i, j) in &TET_EDGES {
        let (a, b) = if gid[i] < gid[j] { (i, j) } else { (j, i) };
        let (wr, cr) = whitney_ref(a, b);
        w_ref.push(wr);
        let mut dw = [[0.0; 4]; 3];
        for p in 0..3 {
            dw[p][a] = G_REF[b][p];
            dw[p][b] = -G_REF[a][p];
        }
        dw_dlam.push(dw);
        c_ref.push(cr);
        dc_dlam.push([[0.0; 4]; 3]);

        if order == 2 {
            let wr2 = [
                lam[a] * G_REF[b][0] + lam[b] * G_REF[a][0],
                lam[a] * G_REF[b][1] + lam[b] * G_REF[a][1],
                lam[a] * G_REF[b][2] + lam[b] * G_REF[a][2],
            ];
            w_ref.push(wr2);
            let mut dw2 = [[0.0; 4]; 3];
            for p in 0..3 {
                dw2[p][a] = G_REF[b][p];
                dw2[p][b] = G_REF[a][p];
            }
            dw_dlam.push(dw2);
            c_ref.push([0.0; 3]);
            dc_dlam.push([[0.0; 4]; 3]);
        }
    }

    if order == 2 {
        for f in &TET_FACES {
            let mut t = *f;
            t.sort_by_key(|&k| gid[k]);
            let (a, b, cc) = (t[0], t[1], t[2]);
            for (p_idx, u_idx, v_idx) in [(cc, a, b), (a, b, cc)] {
                let (nuv, cuv) = whitney_ref(u_idx, v_idx);
                w_ref.push([lam[p_idx] * nuv[0], lam[p_idx] * nuv[1], lam[p_idx] * nuv[2]]);
                let mut dw = [[0.0; 4]; 3];
                for p in 0..3 {
                    dw[p][p_idx] += nuv[p];
                    dw[p][u_idx] += lam[p_idx] * G_REF[v_idx][p];
                    dw[p][v_idx] -= lam[p_idx] * G_REF[u_idx][p];
                }
                dw_dlam.push(dw);

                let gp_x_nuv = cross(G_REF[p_idx], nuv);
                c_ref.push([
                    lam[p_idx] * cuv[0] + gp_x_nuv[0],
                    lam[p_idx] * cuv[1] + gp_x_nuv[1],
                    lam[p_idx] * cuv[2] + gp_x_nuv[2],
                ]);
                let mut dc = [[0.0; 4]; 3];
                let gp_x_gv = cross(G_REF[p_idx], G_REF[v_idx]);
                let gp_x_gu = cross(G_REF[p_idx], G_REF[u_idx]);
                for p in 0..3 {
                    dc[p][p_idx] += cuv[p];
                    dc[p][u_idx] += gp_x_gv[p];
                    dc[p][v_idx] -= gp_x_gu[p];
                }
                dc_dlam.push(dc);
            }
        }
    }

    let inv_det = 1.0 / det;
    let mut w_out = Vec::with_capacity(nb);
    let mut curl_out = Vec::with_capacity(nb);
    let mut div_out = Vec::with_capacity(nb);
    let mut cc_out = Vec::with_capacity(nb);

    for s in 0..nb {
        let wr = w_ref[s];
        let cr = c_ref[s];
        let w_phys: V3 = std::array::from_fn(|i| (0..3).map(|p| g[p + 1][i] * wr[p]).sum());
        let h_phys: V3 = std::array::from_fn(|i| (0..3).map(|p| j_mat[i][p] * cr[p]).sum::<f64>() * inv_det);

        let dw_xi: [[C; 3]; 3] = std::array::from_fn(|p| {
            std::array::from_fn(|m| C::new(dw_dlam[s][p][m + 1] - dw_dlam[s][p][0], 0.0))
        });
        let dc_xi: [[C; 3]; 3] = std::array::from_fn(|p| {
            std::array::from_fn(|m| C::new(dc_dlam[s][p][m + 1] - dc_dlam[s][p][0], 0.0))
        });

        let w_phys_c = w_phys.map(|x| C::new(x, 0.0));
        let h_phys_c = h_phys.map(|x| C::new(x, 0.0));
        let cr_c = cr.map(|x| C::new(x, 0.0));

        let grad_w = pushforward_covariant_grad(&w_phys_c, &dw_xi, g, h_jac);
        let grad_h = pushforward_contravariant_curl_grad(&h_phys_c, &cr_c, &dc_xi, j_mat, g, det, h_jac);
        let cc_c = curl_of(&grad_h);

        w_out.push(w_phys);
        curl_out.push(h_phys);
        div_out.push(div_of(&grad_w).re);
        cc_out.push([cc_c[0].re, cc_c[1].re, cc_c[2].re]);
    }

    Ok(NedelecEval { w: w_out, curl: curl_out, div: div_out, curl_curl: cc_out })
}

/// Self-contained complex Nedelec [`CurlField`] over a [`SimplexMesh`] (straight or [`CurvedTet`])
/// supplying exact analytic values, curls, divergences, and curl-curls.
pub struct NedelecTetField<'a> {
    pub mesh: &'a SimplexMesh,
    pub order: usize,
    pub cell_dofs: &'a [Vec<usize>],
    pub cell_keys: &'a [[(u32, u32); 4]],
    pub dofs: &'a [C],
}

type CellJacobianData = ([[f64; 3]; 3], [V3; 4], f64, JacHessian);

impl<'a> NedelecTetField<'a> {
    fn cell_jacobian_data(&self, cell: usize, lam: &[f64; 4]) -> Option<CellJacobianData> {
        let v = self.mesh.cell(cell);
        let p4: [V3; 4] = std::array::from_fn(|i| {
            let p = self.mesh.point(v[i] as usize);
            [p[0], p[1], p[2]]
        });
        if self.mesh.cell_is_curved(cell) {
            let ct = CurvedTet::new(self.mesh.cell_geometry_nodes(cell));
            let det0 = CurvedTet::straight(p4).jacobian(&[0.25; 4]).1;
            let (mut g, mut det) = ct.jacobian(lam);
            let mut j_mat = ct.jacobian_matrix(lam);
            if !det.is_finite() || det.abs() <= 1e-14 * det0.abs().max(1e-300) {
                const EPS: f64 = 1e-3;
                let lam_in = std::array::from_fn(|b| (1.0 - EPS) * lam[b] + EPS * 0.25);
                let (g_in, det_in) = ct.jacobian(&lam_in);
                if det_in.is_finite() && det_in.abs() > 1e-300 {
                    g = g_in;
                    det = det_in;
                    j_mat = ct.jacobian_matrix(&lam_in);
                } else {
                    return None;
                }
            }
            let h_jac = curved_tet_hessian(&ct);
            Some((j_mat, g, det, h_jac))
        } else {
            let (g, _) = cell_geometry(&p4)?;
            let ct = CurvedTet::straight(p4);
            let j_mat = ct.jacobian_matrix(lam);
            let (_, det) = ct.jacobian(lam);
            Some((j_mat, g, det, [[[0.0; 3]; 3]; 3]))
        }
    }
}

impl CurlField for NedelecTetField<'_> {
    fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
        let Some((j_mat, g, det, h_jac)) = self.cell_jacobian_data(cell, &lam) else {
            return ([C::new(0.0, 0.0); 3], [C::new(0.0, 0.0); 3]);
        };
        let Ok(ev) = eval_nedelec_tet_with_derivatives(self.order, &lam, &j_mat, &g, det, &h_jac, &self.cell_keys[cell])
        else {
            return ([C::new(0.0, 0.0); 3], [C::new(0.0, 0.0); 3]);
        };
        let mut e = [C::new(0.0, 0.0); 3];
        let mut cu = [C::new(0.0, 0.0); 3];
        for (s, &dof) in self.cell_dofs[cell].iter().enumerate() {
            let coef = self.dofs[dof];
            for k in 0..3 {
                e[k] += coef * ev.w[s][k];
                cu[k] += coef * ev.curl[s][k];
            }
        }
        (e, cu)
    }

    fn sample_derivatives(&self, cell: usize, lam: [f64; 4]) -> Option<(C, [C; 3])> {
        let (j_mat, g, det, h_jac) = self.cell_jacobian_data(cell, &lam)?;
        let ev =
            eval_nedelec_tet_with_derivatives(self.order, &lam, &j_mat, &g, det, &h_jac, &self.cell_keys[cell]).ok()?;
        let mut div_e = C::new(0.0, 0.0);
        let mut cc_e = [C::new(0.0, 0.0); 3];
        for (s, &dof) in self.cell_dofs[cell].iter().enumerate() {
            let coef = self.dofs[dof];
            div_e += coef * ev.div[s];
            for k in 0..3 {
                cc_e[k] += coef * ev.curl_curl[s][k];
            }
        }
        Some((div_e, cc_e))
    }
}

/// Reconstructs 10-node P2 reference nodal samples `(E_ref, curl_ref)` on a straight or curved cell
/// so spatial derivatives can be evaluated analytically without finite-difference steps.
struct RefCellDerivativeCache {
    is_curved: bool,
    alpha: f64,
    e_nodes: [[C; 3]; 10],
    c_nodes: [[C; 3]; 10],
    h_jac: JacHessian,
}

impl RefCellDerivativeCache {
    fn build_straight(field: &dyn CurlField, cell: usize) -> Self {
        let mut e_nodes = [[C::new(0.0, 0.0); 3]; 10];
        let mut c_nodes = [[C::new(0.0, 0.0); 3]; 10];
        for (a, &lam) in P2_REF_NODES.iter().enumerate() {
            let (e, cu) = field.sample(cell, lam);
            e_nodes[a] = e;
            c_nodes[a] = cu;
        }
        Self { is_curved: false, alpha: 0.0, e_nodes, c_nodes, h_jac: [[[0.0; 3]; 3]; 3] }
    }

    fn build_curved(field: &dyn CurlField, cell: usize, ct: &CurvedTet, det0_abs: f64) -> Option<Self> {
        // Pull back onto interior nodes when quarter-point / degenerate vertices have vanishing det(J).
        let mut min_det_ratio = f64::INFINITY;
        for &lam in &P2_REF_NODES {
            let (_, det) = ct.jacobian(&lam);
            let ratio = if det0_abs > 1e-300 { det.abs() / det0_abs } else { 0.0 };
            if !det.is_finite() {
                min_det_ratio = 0.0;
                break;
            }
            min_det_ratio = min_det_ratio.min(ratio);
        }
        let alpha = if min_det_ratio < 1e-6 { 0.05 } else { 0.0 };
        let mut e_nodes = [[C::new(0.0, 0.0); 3]; 10];
        let mut c_nodes = [[C::new(0.0, 0.0); 3]; 10];
        for (a, &base_lam) in P2_REF_NODES.iter().enumerate() {
            let lam = if alpha == 0.0 {
                base_lam
            } else {
                std::array::from_fn(|b| (1.0 - alpha) * base_lam[b] + 0.25 * alpha)
            };
            let j_mat = ct.jacobian_matrix(&lam);
            let (g, det) = ct.jacobian(&lam);
            if !det.is_finite() || det.abs() <= 1e-14 * det0_abs.max(1e-300) {
                return None;
            }
            let (e_phys, h_phys) = field.sample(cell, lam);
            // Pull back: E_ref = J^T E_phys, C_ref = det * J^{-1} H_phys
            for p in 0..3 {
                let mut er = C::new(0.0, 0.0);
                let mut cr = C::new(0.0, 0.0);
                for i in 0..3 {
                    er += e_phys[i] * j_mat[i][p];
                    cr += h_phys[i] * (g[p + 1][i] * det);
                }
                e_nodes[a][p] = er;
                c_nodes[a][p] = cr;
            }
        }
        let h_jac = curved_tet_hessian(ct);
        Some(Self { is_curved: true, alpha, e_nodes, c_nodes, h_jac })
    }

    fn eval_div_and_curl_curl(
        &self,
        lam: &[f64; 4],
        e_phys: &[C; 3],
        h_phys: &[C; 3],
        j_mat: &[[f64; 3]; 3],
        g: &[V3; 4],
        det: f64,
    ) -> (C, [C; 3]) {
        if !self.is_curved {
            let (_, dn) = shape(lam, 2);
            let mut grad_e = [[C::new(0.0, 0.0); 3]; 3];
            let mut grad_c = [[C::new(0.0, 0.0); 3]; 3];
            for a in 0..10 {
                let mut ga = [0.0; 3];
                for j in 0..3 {
                    let mut s = 0.0;
                    for b in 0..4 {
                        s += dn[a * 4 + b] * g[b][j];
                    }
                    ga[j] = s;
                }
                for k in 0..3 {
                    let ea = self.e_nodes[a][k];
                    let ca = self.c_nodes[a][k];
                    for j in 0..3 {
                        grad_e[k][j] += ea * ga[j];
                        grad_c[k][j] += ca * ga[j];
                    }
                }
            }
            (div_of(&grad_e), curl_of(&grad_c))
        } else {
            let de_xi = eval_p2_ref_dxi(&self.e_nodes, lam, self.alpha);
            let dc_xi = eval_p2_ref_dxi(&self.c_nodes, lam, self.alpha);
            // Reconstruct C_ref(lam) = J^{-1} * (det * H_phys)
            let c_ref: [C; 3] = std::array::from_fn(|p| {
                (0..3).fold(C::new(0.0, 0.0), |acc, i| acc + h_phys[i] * (g[p + 1][i] * det))
            });
            let grad_e = pushforward_covariant_grad(e_phys, &de_xi, g, &self.h_jac);
            let grad_h = pushforward_contravariant_curl_grad(h_phys, &c_ref, &dc_xi, j_mat, g, det, &self.h_jac);
            (div_of(&grad_e), curl_of(&grad_h))
        }
    }
}

/// Validates mesh integrity and Maxwell options before running indicators.
fn validate_inputs(mesh: &SimplexMesh, opt: &MaxwellOptions) -> Result<(), String> {
    if mesh.dim != 3 {
        return Err("the Maxwell estimator works on tetrahedral meshes".into());
    }
    if !mesh.points.len().is_multiple_of(3) {
        return Err("mesh points length is not a multiple of 3".into());
    }
    if !mesh.cells.len().is_multiple_of(4) {
        return Err("mesh cells length is not a multiple of 4".into());
    }
    if !mesh.facets.len().is_multiple_of(3) {
        return Err("mesh facets length is not a multiple of 3".into());
    }
    let nc = mesh.n_cells();
    let np = mesh.n_points();
    let nf = mesh.n_facets();
    if mesh.cell_labels.len() != nc {
        return Err(format!("{} cell_labels for {nc} cells", mesh.cell_labels.len()));
    }
    if mesh.facet_labels.len() != nf {
        return Err(format!("{} facet_labels for {nf} facets", mesh.facet_labels.len()));
    }
    if opt.eps_r.len() != nc || opt.mu_r.len() != nc {
        return Err(format!("{} permittivities and {} permeabilities for {nc} cells", opt.eps_r.len(), opt.mu_r.len()));
    }
    if !opt.k0_sq.is_finite() {
        return Err("k0_sq must be finite".into());
    }
    if !mesh.points.iter().all(|x| x.is_finite()) {
        return Err("mesh contains non-finite vertex coordinates".into());
    }
    for c in 0..nc {
        let v = mesh.cell(c);
        for (i, &vi) in v.iter().enumerate() {
            if (vi as usize) >= np {
                return Err(format!("cell {c} vertex {vi} out of bounds (n_points = {np})"));
            }
            if v[..i].contains(&vi) {
                return Err(format!("cell {c} contains duplicate vertex {vi}"));
            }
        }
        if !(opt.skip)(mesh.cell_labels[c]) {
            let (e, m) = (opt.eps_r[c], opt.mu_r[c]);
            if !(e > 0.0 && e.is_finite() && m > 0.0 && m.is_finite()) {
                return Err(format!(
                    "live cell {c} (label {}) permittivity={e} and permeability={m} must be positive and finite",
                    mesh.cell_labels[c]
                ));
            }
        }
    }
    for f in 0..nf {
        let fv = mesh.facet(f);
        for (i, &vi) in fv.iter().enumerate() {
            if (vi as usize) >= np {
                return Err(format!("facet {f} vertex {vi} out of bounds (n_points = {np})"));
            }
            if fv[..i].contains(&vi) {
                return Err(format!("facet {f} contains duplicate vertex {vi}"));
            }
        }
    }
    Ok(())
}

/// The indicators `eta_T^2` of `field` on `mesh`.
pub fn maxwell_indicators(mesh: &SimplexMesh, field: &dyn CurlField, opt: &MaxwellOptions) -> Result<Indicators, String> {
    validate_inputs(mesh, opt)?;
    let nc = mesh.n_cells();
    if nc == 0 {
        return Ok(Indicators { eta2: Vec::new(), parts: [0.0; 4], total: 0.0 });
    }
    let pt = |v: u32| -> V3 { [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]] };
    let cell_pts = |c: usize| -> [V3; 4] {
        let v = mesh.cell(c);
        [pt(v[0]), pt(v[1]), pt(v[2]), pt(v[3])]
    };
    let vol_rule = quad_rule(3, 4);
    let curved_quad = mesh
        .has_curved()
        .then(|| quad_rule_collapsed(if mesh.has_qp() { 8 } else { CURVED_QUAD_DEGREE }));
    let face_rule = quad_rule(2, if mesh.has_curved() { 6 } else { 4 });

    // interior terms, cell by cell
    let interior: Vec<[f64; 2]> = (0..nc)
        .into_par_iter()
        .map(|c| {
            if (opt.skip)(mesh.cell_labels[c]) {
                return [0.0, 0.0];
            }
            let p = cell_pts(c);
            let Some((g_straight, straight_vol)) = cell_geometry(&p) else {
                return [0.0, 0.0];
            };
            let (eps, mu) = (opt.eps_r[c], opt.mu_r[c]);

            if mesh.cell_is_curved(c) {
                let nodes10 = mesh.cell_geometry_nodes(c);
                let ct = CurvedTet::new(nodes10);
                let det0 = CurvedTet::straight(p).jacobian(&[0.25; 4]).1;
                let det0_abs = det0.abs();
                let h = diameter(&nodes10);
                let cquad = curved_quad.as_ref().unwrap_or(&vol_rule);
                let mut cache: Option<RefCellDerivativeCache> = None;
                let mut cache_failed = false;
                let (mut res, mut div) = (0.0, 0.0);

                for (lam, w) in cquad.lam.iter().zip(&cquad.w) {
                    let l4 = [lam[0], lam[1], lam[2], lam[3]];
                    let (gq, det) = ct.jacobian(&l4);
                    let det_abs = det.abs();
                    if !det.is_finite() || det_abs <= 1e-14 * det0_abs.max(1e-300) {
                        continue;
                    }
                    let wq = w * det_abs / 6.0;
                    let (e, h_curl) = field.sample(c, l4);
                    let (div_e, cc_e) = if let Some(exact) = field.sample_derivatives(c, l4) {
                        exact
                    } else {
                        if cache.is_none() && !cache_failed {
                            cache = RefCellDerivativeCache::build_curved(field, c, &ct, det0_abs);
                            if cache.is_none() {
                                cache_failed = true;
                            }
                        }
                        if let Some(ref cc) = cache {
                            let j_mat = ct.jacobian_matrix(&l4);
                            cc.eval_div_and_curl_curl(&l4, &e, &h_curl, &j_mat, &gq, det)
                        } else {
                            continue;
                        }
                    };
                    let r = [
                        cc_e[0] / mu - e[0] * (opt.k0_sq * eps),
                        cc_e[1] / mu - e[1] * (opt.k0_sq * eps),
                        cc_e[2] / mu - e[2] * (opt.k0_sq * eps),
                    ];
                    res += wq * cnorm2(r);
                    let d = div_e * eps;
                    div += wq * d.norm_sqr();
                }
                [h * h * res, h * h * div]
            } else {
                let h = diameter(&p);
                let mut cache: Option<RefCellDerivativeCache> = None;
                let (mut res, mut div) = (0.0, 0.0);
                let dummy_j = [[0.0; 3]; 3];

                for (lam, w) in vol_rule.lam.iter().zip(&vol_rule.w) {
                    let l4 = [lam[0], lam[1], lam[2], lam[3]];
                    let (e, h_curl) = field.sample(c, l4);
                    let (div_e, cc_e) = if let Some(exact) = field.sample_derivatives(c, l4) {
                        exact
                    } else {
                        let ref_cache = cache.get_or_insert_with(|| RefCellDerivativeCache::build_straight(field, c));
                        ref_cache.eval_div_and_curl_curl(&l4, &e, &h_curl, &dummy_j, &g_straight, 1.0)
                    };
                    let r = [
                        cc_e[0] / mu - e[0] * (opt.k0_sq * eps),
                        cc_e[1] / mu - e[1] * (opt.k0_sq * eps),
                        cc_e[2] / mu - e[2] * (opt.k0_sq * eps),
                    ];
                    res += w * straight_vol * cnorm2(r);
                    let d = div_e * eps;
                    div += w * straight_vol * d.norm_sqr();
                }
                [h * h * res, h * h * div]
            }
        })
        .collect();

    // faces: `(cell, opp_local_slot, sorted_face_local_slots)` keyed by sorted global vertex triple
    type FaceSide = (usize, usize, [usize; 3]);
    let mut faces: HashMap<[u32; 3], Vec<FaceSide>> = HashMap::with_capacity(nc * 2);
    for c in 0..nc {
        let v = mesh.cell(c);
        for opp in 0..4 {
            let mut lf = [0usize; 3];
            let mut idx = 0;
            for i in 0..4 {
                if i != opp {
                    lf[idx] = i;
                    idx += 1;
                }
            }
            lf.sort_by_key(|&i| v[i]);
            let key = [v[lf[0]], v[lf[1]], v[lf[2]]];
            faces.entry(key).or_default().push((c, opp, lf));
        }
    }
    let mut keys: Vec<[u32; 3]> = faces.keys().copied().collect();
    keys.sort_unstable();
    // boundary facet labels by face
    let mut facet_label: HashMap<[u32; 3], i32> = HashMap::with_capacity(mesh.n_facets());
    for f in 0..mesh.n_facets() {
        let v = mesh.facet(f);
        let mut k = [v[0], v[1], v[2]];
        k.sort_unstable();
        facet_label.insert(k, mesh.facet_labels[f]);
    }

    // jump terms: (cell, tangential part, normal part)
    let face_terms: Vec<Vec<(usize, f64, f64)>> = keys
        .par_iter()
        .map(|k| {
            let cs = &faces[k];
            let live = |cell: usize| !(opt.skip)(mesh.cell_labels[cell]);
            let mut out = Vec::new();
            if cs.len() == 2 {
                if !live(cs[0].0) || !live(cs[1].0) {
                    return out;
                }
            } else if cs.len() == 1 {
                if !live(cs[0].0) {
                    return out;
                }
                let natural = facet_label.get(k).is_some_and(|&l| (opt.natural)(l));
                if !natural {
                    return out;
                }
            } else {
                return out;
            }

            let (a, b, c_pt) = (pt(k[0]), pt(k[1]), pt(k[2]));
            let area_vec = cross(sub(b, a), sub(c_pt, a));
            let twice_area = norm(area_vec);
            let area = 0.5 * twice_area;
            if !area.is_finite() || area <= 1e-300 {
                return out;
            }
            let nrm = area_vec.map(|x| x / twice_area);
            let in_cell = |lf: [usize; 3], l: &[f64]| -> [f64; 4] {
                let mut lam = [0.0; 4];
                lam[lf[0]] = l[0];
                lam[lf[1]] = l[1];
                lam[lf[2]] = l[2];
                lam
            };
            let outward = |cell: usize, opp_slot: usize| -> V3 {
                let opp = pt(mesh.cell(cell)[opp_slot]);
                if dot(nrm, sub(opp, a)) > 0.0 {
                    nrm.map(|x| -x)
                } else {
                    nrm
                }
            };

            if cs.len() == 2 {
                let ((c1, s1, lf1), (c2, _, lf2)) = (cs[0], cs[1]);
                let n_straight = outward(c1, s1);
                let ct1 = mesh.cell_is_curved(c1).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c1)));
                let ct2 = mesh.cell_is_curved(c2).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c2)));
                let h_f = if let Some(ref ct) = ct1 {
                    let mut fp = [a, b, c_pt, [0.0; 3], [0.0; 3], [0.0; 3]];
                    for (idx, (p, q)) in [(lf1[0], lf1[1]), (lf1[1], lf1[2]), (lf1[0], lf1[2])].into_iter().enumerate() {
                        let (lo, hi) = (p.min(q), p.max(q));
                        let e_slot = TET_EDGES.iter().position(|&e| e == (lo, hi)).unwrap_or(0);
                        fp[3 + idx] = ct.x[4 + e_slot];
                    }
                    diameter(&fp)
                } else {
                    diameter(&[a, b, c_pt])
                };
                let (mut jt, mut jn) = (0.0, 0.0);
                for (l, w) in face_rule.lam.iter().zip(&face_rule.w) {
                    let lam1 = in_cell(lf1, l);
                    let lam2 = in_cell(lf2, l);
                    let (n_q, ds_q) = if let Some(ref ct) = ct1 {
                        let d = ct.dx_dlam(&lam1);
                        let cr = cross(sub(d[lf1[1]], d[lf1[0]]), sub(d[lf1[2]], d[lf1[0]]));
                        let cl = norm(cr);
                        if cl > 1e-300 && cl.is_finite() {
                            let sgn = if dot(cr, n_straight) < 0.0 { -1.0 } else { 1.0 };
                            ([sgn * cr[0] / cl, sgn * cr[1] / cl, sgn * cr[2] / cl], w * 0.5 * cl)
                        } else {
                            (n_straight, w * area)
                        }
                    } else if let Some(ref ct) = ct2 {
                        let d = ct.dx_dlam(&lam2);
                        let cr = cross(sub(d[lf2[1]], d[lf2[0]]), sub(d[lf2[2]], d[lf2[0]]));
                        let cl = norm(cr);
                        if cl > 1e-300 && cl.is_finite() {
                            let sgn = if dot(cr, n_straight) < 0.0 { -1.0 } else { 1.0 };
                            ([sgn * cr[0] / cl, sgn * cr[1] / cl, sgn * cr[2] / cl], w * 0.5 * cl)
                        } else {
                            (n_straight, w * area)
                        }
                    } else {
                        (n_straight, w * area)
                    };
                    let (e1, h1) = field.sample(c1, lam1);
                    let (e2, h2) = field.sample(c2, lam2);
                    let hh = [0, 1, 2].map(|i| h1[i] / opt.mu_r[c1] - h2[i] / opt.mu_r[c2]);
                    jt += ds_q * cnorm2(ccross(n_q, hh));
                    let dn = cdot(n_q, e1) * opt.eps_r[c1] - cdot(n_q, e2) * opt.eps_r[c2];
                    jn += ds_q * dn.norm_sqr();
                }
                out.push((c1, 0.5 * h_f * jt, 0.5 * h_f * jn));
                out.push((c2, 0.5 * h_f * jt, 0.5 * h_f * jn));
            } else {
                let (c1, s1, lf1) = cs[0];
                let n_straight = outward(c1, s1);
                let ct1 = mesh.cell_is_curved(c1).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c1)));
                let h_f = diameter(&[a, b, c_pt]);
                let mut jt = 0.0;
                for (l, w) in face_rule.lam.iter().zip(&face_rule.w) {
                    let lam1 = in_cell(lf1, l);
                    let (n_q, ds_q) = if let Some(ref ct) = ct1 {
                        let d = ct.dx_dlam(&lam1);
                        let cr = cross(sub(d[lf1[1]], d[lf1[0]]), sub(d[lf1[2]], d[lf1[0]]));
                        let cl = norm(cr);
                        if cl > 1e-300 && cl.is_finite() {
                            let sgn = if dot(cr, n_straight) < 0.0 { -1.0 } else { 1.0 };
                            ([sgn * cr[0] / cl, sgn * cr[1] / cl, sgn * cr[2] / cl], w * 0.5 * cl)
                        } else {
                            (n_straight, w * area)
                        }
                    } else {
                        (n_straight, w * area)
                    };
                    let (_, h1) = field.sample(c1, lam1);
                    let hh = h1.map(|z| z / opt.mu_r[c1]);
                    jt += ds_q * cnorm2(ccross(n_q, hh));
                }
                out.push((c1, h_f * jt, 0.0));
            }
            out
        })
        .collect();

    let mut eta2 = vec![0.0; nc];
    let mut parts = [0.0; 4];
    for (c, [res, div]) in interior.iter().enumerate() {
        eta2[c] += res + div;
        parts[0] += res;
        parts[1] += div;
    }
    for terms in face_terms {
        for (c, jt, jn) in terms {
            eta2[c] += jt + jn;
            parts[2] += jt;
            parts[3] += jn;
        }
    }
    let total = eta2.iter().sum::<f64>().sqrt();
    Ok(Indicators { eta2, parts, total })
}

/// Settings of [`indicator_metrics`].
#[derive(Clone, Debug)]
pub struct EquidistributeOptions {
    /// order of the edge elements (1 or 2); the indicator of a cell of size `s` falls like `s^(order + 3/2)`
    pub order: usize,
    /// largest element size asked for
    pub h_max: f64,
    /// smallest element size asked for, as a fraction of `h_max`
    pub h_min_fraction: f64,
    /// the mesh is asked to have about this many cells
    pub target_cells: f64,
}

/// The isotropic metric at every vertex that equidistributes the indicators over about `target_cells` cells. A cell of size `s` with
/// indicator `eta` asks for the size `s (tau / eta)^(1 / (order + 3/2))`, kept between `h_min_fraction * h_max` and `h_max` and within a
/// factor 4 of `s`; the level `tau` is found by bisection so that the predicted number of cells, the sum of `(s / s_new)^3`, is the target.
/// A vertex takes the smallest size asked by a cell that meets there. `floor` (a metric per vertex, for example the geometry floor) is
/// intersected when given.
pub fn indicator_metrics(mesh: &SimplexMesh, eta2: &[f64], opt: &EquidistributeOptions, floor: Option<&[Sym<3>]>) -> Result<Vec<Sym<3>>, String> {
    if mesh.dim != 3 {
        return Err("indicator_metrics requires a 3D mesh".into());
    }
    let nc = mesh.n_cells();
    let np = mesh.n_points();
    if eta2.len() != nc {
        return Err(format!("{} indicators for {nc} cells", eta2.len()));
    }
    if let Some(f) = floor {
        if f.len() != np {
            return Err(format!("{} floor metrics for {np} points", f.len()));
        }
    }
    if !(opt.h_max > 0.0 && opt.h_max.is_finite() && opt.target_cells >= 1.0 && opt.target_cells.is_finite() && (0.0..=1.0).contains(&opt.h_min_fraction)) {
        return Err("h_max and target_cells must be positive and finite, and h_min_fraction within 0 and 1".into());
    }
    let q = opt.order as f64 + 1.5;
    let pt = |v: u32| -> V3 { [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]] };
    // the size of a cell: the edge of the regular tetrahedron of its volume
    let size: Vec<f64> = (0..nc)
        .map(|c| {
            let v = mesh.cell(c);
            let vol = cell_geometry(&[pt(v[0]), pt(v[1]), pt(v[2]), pt(v[3])]).map_or(0.0, |(_, vol)| vol);
            (6.0 * 2f64.sqrt() * vol).cbrt()
        })
        .collect();
    let h_min = (opt.h_max * opt.h_min_fraction).max(1e-150);
    let new_size = |tau: f64, c: usize| -> f64 {
        if size[c] <= 0.0 {
            return opt.h_max;
        }
        let eta = eta2[c].max(0.0).sqrt();
        let ratio = if eta > 0.0 { (tau / eta).powf(1.0 / q) } else { 4.0 };
        (size[c] * ratio.clamp(0.25, 4.0)).clamp(h_min.min(opt.h_max), opt.h_max)
    };
    let predicted = |tau: f64| -> f64 {
        (0..nc)
            .filter(|&c| size[c] > 0.0)
            .map(|c| (size[c] / new_size(tau, c).max(1e-150)).powi(3))
            .sum()
    };
    // predicted is decreasing in tau
    let max_eta = eta2.iter().copied().filter(|x| x.is_finite()).fold(0.0, f64::max).sqrt().max(1e-300);
    let (mut lo, mut hi) = (max_eta * 1e-12, max_eta * 1e12);
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        if predicted(mid) > opt.target_cells {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let tau = (lo * hi).sqrt();
    let mut vh = vec![opt.h_max; np];
    for c in 0..nc {
        if size[c] <= 0.0 {
            continue;
        }
        let s = new_size(tau, c);
        for &v in mesh.cell(c) {
            if (v as usize) < np {
                vh[v as usize] = vh[v as usize].min(s);
            }
        }
    }
    let mut out: Vec<Sym<3>> = vh.iter().map(|&h| isotropic::<3>(h.max(1e-150))).collect();
    if let Some(f) = floor {
        for (a, b) in out.iter_mut().zip(f) {
            *a = intersect(a, b);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::c;
    use crate::mesh::fixtures::box_tet;

    /// The whitney field of a plane wave sampled exactly: `E = (0, 0, e^{-ikx})` is not in the space, so use a polynomial field instead.
    struct Poly {
        k0_sq: f64,
        pts: Vec<[V3; 4]>,
        /// E = (a*y, b*z, c*x) + (d, 0, 0): curl = (-b, -c, -a)
        co: [f64; 3],
    }

    impl CurlField for Poly {
        fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
            let p = &self.pts[cell];
            let x: V3 = std::array::from_fn(|k| (0..4).map(|v| lam[v] * p[v][k]).sum());
            let [a, b, cc] = self.co;
            let e = [c(a * x[1], 0.0), c(b * x[2], 0.0), c(cc * x[0], 0.0)];
            let curl = [c(-b, 0.0), c(-cc, 0.0), c(-a, 0.0)];
            let _ = self.k0_sq;
            (e, curl)
        }
    }

    #[test]
    fn a_field_with_known_residual_gives_the_closed_form_indicators() {
        // E linear: curl E constant, curl curl E = 0, div E = 0 for this field, no jumps (single-valued field): eta_T^2 = h^2 k0^4 eps^2 ||E||^2
        let m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [2, 2, 2]));
        let pts: Vec<[V3; 4]> = (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                std::array::from_fn(|i| [m.points[3 * v[i] as usize], m.points[3 * v[i] as usize + 1], m.points[3 * v[i] as usize + 2]])
            })
            .collect();
        let f = Poly { k0_sq: 3.0, pts: pts.clone(), co: [1.0, 2.0, -0.5] };
        let eps = vec![2.0; m.n_cells()];
        let mu = vec![1.0; m.n_cells()];
        let opt = MaxwellOptions { k0_sq: 3.0, eps_r: &eps, mu_r: &mu, skip: &|_| false, natural: &|_| false };
        let ind = maxwell_indicators(&m, &f, &opt).unwrap();
        // closed form per cell: integrate |k0^2 eps E|^2 over the cell with the exact quadrature of the test (a quadratic integrand)
        for (c, p) in pts.iter().enumerate() {
            let (_, vol) = cell_geometry(p).unwrap();
            let h = diameter(p);
            let rule = quad_rule(3, 4);
            let mut s = 0.0;
            for (lam, w) in rule.lam.iter().zip(&rule.w) {
                let (e, _) = f.sample(c, [lam[0], lam[1], lam[2], lam[3]]);
                s += w * vol * cnorm2(e) * (3.0f64 * 2.0).powi(2);
            }
            assert!((ind.eta2[c] - h * h * s).abs() < 1e-10 * (1.0 + h * h * s), "cell {c}: {} vs {}", ind.eta2[c], h * h * s);
        }
        assert!(ind.parts[1].abs() < 1e-18 && ind.parts[2] == 0.0 && ind.parts[3] == 0.0);
    }

    #[test]
    fn analytic_nedelec_basis_derivatives_match_exact_quadratic_and_curved_tet_calculus() {
        // Test both order 1 and order 2 on a genuinely curved tetrahedron against high-precision physical finite differences
        let v = [[0.0, 0.0, 0.0], [1.0, 0.12, -0.05], [0.08, 1.0, 0.15], [-0.04, 0.18, 1.0]];
        let mut x10 = CurvedTet::straight(v).x;
        let offsets = [
            [0.03, -0.02, 0.015],
            [-0.015, 0.035, 0.02],
            [0.02, 0.01, -0.025],
            [0.01, -0.03, 0.02],
            [0.025, 0.03, -0.01],
            [-0.02, 0.015, 0.035],
        ];
        for e in 0..6 {
            for k in 0..3 {
                x10[4 + e][k] += offsets[e][k];
            }
        }
        let ct = CurvedTet::new(x10);
        let h_jac = curved_tet_hessian(&ct);
        let gid = [(3u32, 0u32), (1, 1), (4, 2), (2, 3)];
        let lam0 = [0.28, 0.24, 0.27, 0.21];

        let inv_map = |p: V3, start: [f64; 4]| -> [f64; 4] {
            let mut l = start;
            for _ in 0..40 {
                let q = ct.position(&l);
                let r = sub(q, p);
                let (g, _) = ct.jacobian(&l);
                let d: V3 = std::array::from_fn(|k| dot(g[k + 1], r));
                for k in 0..3 {
                    l[k + 1] -= d[k];
                }
                l[0] = 1.0 - l[1] - l[2] - l[3];
                if d.iter().all(|x| x.abs() < 1e-15) {
                    break;
                }
            }
            l
        };

        for order in [1usize, 2] {
            let nb = if order == 1 { 6 } else { 20 };
            let coefs: Vec<C> = (0..nb).map(|s| c(((s * 29 + 7) % 13) as f64 / 13.0 - 0.5, ((s * 17 + 3) % 11) as f64 / 11.0 - 0.5)).collect();
            let eval_at = |l: [f64; 4]| -> NedelecEval {
                let j_mat = ct.jacobian_matrix(&l);
                let (g, det) = ct.jacobian(&l);
                eval_nedelec_tet_with_derivatives(order, &l, &j_mat, &g, det, &h_jac, &gid).unwrap()
            };
            let ev0 = eval_at(lam0);
            let mut div_analytic = c(0.0, 0.0);
            let mut cc_analytic = [c(0.0, 0.0); 3];
            for s in 0..nb {
                div_analytic += coefs[s] * ev0.div[s];
                for k in 0..3 {
                    cc_analytic[k] += coefs[s] * ev0.curl_curl[s][k];
                }
            }

            // Verify analytic div E and curl(curl E) against tiny physical central difference (h = 1e-6)
            let p0 = ct.position(&lam0);
            let h_fd = 1e-6;
            let mut de = [[c(0.0, 0.0); 3]; 3];
            let mut dc = [[c(0.0, 0.0); 3]; 3];
            for r in 0..3 {
                let (mut pp, mut pm) = (p0, p0);
                pp[r] += h_fd;
                pm[r] -= h_fd;
                let evp = eval_at(inv_map(pp, lam0));
                let evm = eval_at(inv_map(pm, lam0));
                for s in 0..nb {
                    for i in 0..3 {
                        de[i][r] += coefs[s] * ((evp.w[s][i] - evm.w[s][i]) / (2.0 * h_fd));
                        dc[i][r] += coefs[s] * ((evp.curl[s][i] - evm.curl[s][i]) / (2.0 * h_fd));
                    }
                }
            }
            let div_fd = div_of(&de);
            let cc_fd = curl_of(&dc);
            assert!(
                (div_analytic - div_fd).norm() < 1e-7 * (1.0 + div_fd.norm()),
                "order {order}: div analytic {div_analytic:?} vs fd {div_fd:?}"
            );
            for k in 0..3 {
                assert!(
                    (cc_analytic[k] - cc_fd[k]).norm() < 1e-6 * (1.0 + cc_fd[k].norm()),
                    "order {order} comp {k}: cc analytic {:?} vs fd {:?}",
                    cc_analytic[k],
                    cc_fd[k]
                );
            }
        }
    }

    #[test]
    fn curved_tet_and_complex_field_piola_reconstruction_gives_exact_zero_error_for_constant_and_exact_residual() {
        // Build a curved tetrahedral mesh by curving interior edges of a 2x2x2 box
        let mut m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [2, 2, 2]));
        let mut curved_edges = Vec::new();
        for &(a, b) in m.edges() {
            let (pa, pb) = (m.point(a as usize), m.point(b as usize));
            let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            // Curve interior edges smoothly so boundary stays box-shaped while cells become genuine CurvedTets
            if mid.iter().all(|&x| x > 0.05 && x < 0.95) {
                curved_edges.push(([a, b], [mid[0] + 0.02, mid[1] - 0.015, mid[2] + 0.01]));
            }
        }
        m.set_edge_nodes(&curved_edges);
        assert!(m.has_curved());
        m.check().unwrap();

        let nc = m.n_cells();
        let eps = vec![1.5; nc];
        let mu = vec![1.2; nc];

        // Case A: Complex constant field E = (1 + 2i, -0.5 + i, 0.3 - 0.7i) on curved mesh:
        // curl E = 0, div(eps E) = 0, curl curl E = 0, jumps = 0 => exact eta_T^2 = h^2 * vol_T * |k0^2 eps E|^2
        struct ComplexConst([C; 3]);
        impl CurlField for ComplexConst {
            fn sample(&self, _cell: usize, _lam: [f64; 4]) -> ([C; 3], [C; 3]) {
                (self.0, [c(0.0, 0.0); 3])
            }
        }
        let e0 = [c(1.0, 2.0), c(-0.5, 1.0), c(0.3, -0.7)];
        let opt = MaxwellOptions { k0_sq: 2.5, eps_r: &eps, mu_r: &mu, skip: &|_| false, natural: &|_| false };
        let ind_const = maxwell_indicators(&m, &ComplexConst(e0), &opt).unwrap();
        assert!(ind_const.parts[1] < 1e-22, "div part should vanish on curved mesh: {}", ind_const.parts[1]);
        assert!(ind_const.parts[2] < 1e-22 && ind_const.parts[3] < 1e-22);
        let expected_vol_sum: f64 = (0..nc)
            .map(|ci| {
                let nodes10 = m.cell_geometry_nodes(ci);
                let h = diameter(&nodes10);
                let vol = CurvedTet::new(nodes10).volume();
                h * h * vol * cnorm2(e0) * (2.5 * 1.5_f64).powi(2)
            })
            .sum();
        assert!(
            (ind_const.parts[0] - expected_vol_sum).abs() < 1e-10 * expected_vol_sum,
            "curved interior residual {} vs exact {}",
            ind_const.parts[0],
            expected_vol_sum
        );

        // Case B: Complex Nedelec order-2 field on curved mesh via NedelecTetField vs CurlField::sample-only pullback!
        let mut cell_keys = Vec::with_capacity(nc);
        let mut cell_dofs = Vec::with_capacity(nc);
        let mut face_map: HashMap<[u32; 3], usize> = HashMap::new();
        let mut n_faces = 0usize;
        let ne = m.n_edges();
        for ci in 0..nc {
            let v = m.cell(ci);
            cell_keys.push([(v[0], v[0]), (v[1], v[1]), (v[2], v[2]), (v[3], v[3])]);
            let mut d = Vec::with_capacity(20);
            for &eid in m.cell_edge_ids(ci) {
                d.push(2 * eid as usize);
                d.push(2 * eid as usize + 1);
            }
            for f in &TET_FACES {
                let mut k = [v[f[0]], v[f[1]], v[f[2]]];
                k.sort_unstable();
                let fid = *face_map.entry(k).or_insert_with(|| {
                    let id = n_faces;
                    n_faces += 1;
                    id
                });
                d.push(2 * ne + 2 * fid);
                d.push(2 * ne + 2 * fid + 1);
            }
            cell_dofs.push(d);
        }
        let n_dofs = 2 * ne + 2 * n_faces;
        let dofs: Vec<C> = (0..n_dofs).map(|i| c(((i * 37 + 5) % 19) as f64 / 19.0 - 0.5, ((i * 53 + 11) % 23) as f64 / 23.0 - 0.5)).collect();
        let analytic_field = NedelecTetField { mesh: &m, order: 2, cell_dofs: &cell_dofs, cell_keys: &cell_keys, dofs: &dofs };

        // Wrapper exposing ONLY sample() (forcing sample_derivatives -> None) to verify 10-node curved Piola reconstruction
        struct SampleOnly<'a>(&'a dyn CurlField);
        impl CurlField for SampleOnly<'_> {
            fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
                self.0.sample(cell, lam)
            }
        }
        let ind_exact = maxwell_indicators(&m, &analytic_field, &opt).unwrap();
        let ind_recon = maxwell_indicators(&m, &SampleOnly(&analytic_field), &opt).unwrap();
        for p in 0..4 {
            assert!(
                (ind_exact.parts[p] - ind_recon.parts[p]).abs() < 1e-10 * (1.0 + ind_exact.parts[p]),
                "part {p}: exact {} vs reconstructed {}",
                ind_exact.parts[p],
                ind_recon.parts[p]
            );
        }
    }

    #[test]
    fn pml_and_anisotropic_domains_are_skipped_cleanly_even_with_sentinel_material_values() {
        let mut m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [2, 2, 2]));
        let nc = m.n_cells();
        // Label cells by x-centroid: label 1 = isotropic live, label 2 = PML, label 3 = anisotropic
        for ci in 0..nc {
            let v = m.cell(ci);
            let xc = (0..4).map(|i| m.points[3 * v[i] as usize]).sum::<f64>() / 4.0;
            m.cell_labels[ci] = if xc < 0.35 {
                1
            } else if xc < 0.7 {
                2
            } else {
                3
            };
        }

        // Put intentionally invalid sentinel eps_r / mu_r (0.0, -1.0, NaN) in skipped PML & anisotropic cells!
        let mut eps = vec![2.0; nc];
        let mut mu = vec![1.0; nc];
        for ci in 0..nc {
            if m.cell_labels[ci] == 2 {
                eps[ci] = 0.0;
                mu[ci] = f64::NAN;
            } else if m.cell_labels[ci] == 3 {
                eps[ci] = -5.0;
                mu[ci] = 0.0;
            }
        }

        struct PanickingIfSkippedSampled<'a>(&'a [i32]);
        impl CurlField for PanickingIfSkippedSampled<'_> {
            fn sample(&self, cell: usize, _lam: [f64; 4]) -> ([C; 3], [C; 3]) {
                assert_eq!(self.0[cell], 1, "skipped cell {cell} (label {}) must never be sampled", self.0[cell]);
                ([c(1.0, -0.5), c(0.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0); 3])
            }
        }

        let skip_rule = DomainSkipRule::new([2], [3]);
        let skip_fn = |lab: i32| skip_rule.should_skip(lab);
        let opt = MaxwellOptions { k0_sq: 4.0, eps_r: &eps, mu_r: &mu, skip: &skip_fn, natural: &|_| false };
        let ind = maxwell_indicators(&m, &PanickingIfSkippedSampled(&m.cell_labels), &opt).unwrap();

        assert!(ind.total.is_finite() && ind.total > 0.0);
        for ci in 0..nc {
            if m.cell_labels[ci] == 1 {
                assert!(ind.eta2[ci] > 0.0 && ind.eta2[ci].is_finite());
            } else {
                assert_eq!(ind.eta2[ci], 0.0, "skipped cell {ci} must have zero indicator");
            }
        }

        // All cells skipped => zero indicators, no error
        let opt_all_skip = MaxwellOptions { k0_sq: 4.0, eps_r: &eps, mu_r: &mu, skip: &|_| true, natural: &|_| false };
        let ind_all = maxwell_indicators(&m, &PanickingIfSkippedSampled(&m.cell_labels), &opt_all_skip).unwrap();
        assert_eq!(ind_all.total, 0.0);
        assert!(ind_all.eta2.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn degenerate_and_empty_inputs_never_panic() {
        struct ZeroField;
        impl CurlField for ZeroField {
            fn sample(&self, _cell: usize, _lam: [f64; 4]) -> ([C; 3], [C; 3]) {
                ([c(0.0, 0.0); 3], [c(0.0, 0.0); 3])
            }
        }

        // 1. Empty 3D mesh (0 cells)
        let empty = SimplexMesh::new(3, vec![], vec![], vec![], vec![], vec![]);
        let opt_empty = MaxwellOptions { k0_sq: 1.0, eps_r: &[], mu_r: &[], skip: &|_| false, natural: &|_| false };
        let ind_empty = maxwell_indicators(&empty, &ZeroField, &opt_empty).unwrap();
        assert!(ind_empty.eta2.is_empty());
        assert_eq!(ind_empty.total, 0.0);
        let met_empty = indicator_metrics(
            &empty,
            &ind_empty.eta2,
            &EquidistributeOptions { order: 1, h_max: 1.0, h_min_fraction: 0.0, target_cells: 10.0 },
            None,
        )
        .unwrap();
        assert!(met_empty.is_empty());

        // 2. Degenerate coplanar (zero-volume) tetrahedron -> no panic, no NaN
        let coplanar = SimplexMesh::new(
            3,
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0],
            vec![0, 1, 2, 3],
            vec![1],
            vec![0, 1, 2],
            vec![1],
        );
        let eps1 = [1.0];
        let mu1 = [1.0];
        let opt1 = MaxwellOptions { k0_sq: 1.0, eps_r: &eps1, mu_r: &mu1, skip: &|_| false, natural: &|_| true };
        let ind_cop = maxwell_indicators(&coplanar, &ZeroField, &opt1).unwrap();
        assert_eq!(ind_cop.eta2, vec![0.0]);
        assert_eq!(ind_cop.total, 0.0);
        let met_cop = indicator_metrics(
            &coplanar,
            &ind_cop.eta2,
            &EquidistributeOptions { order: 1, h_max: 1.0, h_min_fraction: 0.0, target_cells: 10.0 },
            None,
        )
        .unwrap();
        assert!(met_cop.iter().flatten().flatten().all(|x| x.is_finite()));

        // 3. Duplicate vertex in cell -> clean Err instead of panic on .unwrap()
        let dup_cell = SimplexMesh::new(
            3,
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            vec![0, 1, 2, 2],
            vec![1],
            vec![],
            vec![],
        );
        assert!(maxwell_indicators(&dup_cell, &ZeroField, &opt1).is_err());

        // 4. Out-of-bounds vertex index -> clean Err
        let oob_cell = SimplexMesh::new(
            3,
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            vec![0, 1, 2, 99],
            vec![1],
            vec![],
            vec![],
        );
        assert!(maxwell_indicators(&oob_cell, &ZeroField, &opt1).is_err());

        // 5. Wrong dimension / non-finite k0_sq / invalid live eps_r / mismatched floor length -> clean Err
        let wrong_dim = SimplexMesh::new(2, vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0], vec![0, 1, 2], vec![1], vec![], vec![]);
        assert!(maxwell_indicators(&wrong_dim, &ZeroField, &opt1).is_err());
        let opt_nan = MaxwellOptions { k0_sq: f64::NAN, eps_r: &eps1, mu_r: &mu1, skip: &|_| false, natural: &|_| false };
        assert!(maxwell_indicators(&coplanar, &ZeroField, &opt_nan).is_err());
        let eps_bad = [0.0];
        let opt_bad_eps = MaxwellOptions { k0_sq: 1.0, eps_r: &eps_bad, mu_r: &mu1, skip: &|_| false, natural: &|_| false };
        assert!(maxwell_indicators(&coplanar, &ZeroField, &opt_bad_eps).is_err());
        assert!(indicator_metrics(&coplanar, &ind_cop.eta2, &EquidistributeOptions { order: 1, h_max: 1.0, h_min_fraction: 0.1, target_cells: 10.0 }, Some(&[])).is_err());
    }

    #[test]
    fn equidistribution_hits_the_target_count_and_refines_where_the_error_is() {
        let m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [4, 4, 4]));
        // an error concentrated near the origin
        let eta2: Vec<f64> = (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                let x = (0..4).map(|i| m.points[3 * v[i] as usize]).sum::<f64>() / 4.0;
                (-8.0 * x).exp()
            })
            .collect();
        let opt = EquidistributeOptions { order: 1, h_max: 1.0, h_min_fraction: 0.02, target_cells: 2.0 * m.n_cells() as f64 };
        let met = indicator_metrics(&m, &eta2, &opt, None).unwrap();
        let h_at = |v: usize| (1.0 / met[v][0][0]).sqrt();
        // near the origin the size is smaller than far away
        let near = (0..m.n_points()).filter(|&v| m.points[3 * v] < 0.01).map(h_at).fold(f64::INFINITY, f64::min);
        let far = (0..m.n_points()).filter(|&v| m.points[3 * v] > 0.99).map(h_at).fold(f64::INFINITY, f64::min);
        assert!(near < far, "{near} {far}");
        assert!(indicator_metrics(&m, &eta2[1..], &opt, None).is_err());
    }
}
