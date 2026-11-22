// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Goal-oriented a posteriori flux-jump residual estimator and metric sizing for electrostatics.
//!
//! Pure physics-agnostic finite-element math over [`SimplexMesh`], terminal potential vectors,
//! scalar (`&[f64]`) or symmetric tensor (`&[M3]`) relative permittivity, optional per-facet
//! Robin weights, and Dirichlet/sheet facet masks. Supports 2D triangular and 3D tetrahedral
//! meshes (including isoparametric [`CurvedTet`] cells), deterministic sorted-face parallel
//! reduction (`NYX_THREADS` bit-identical), goal-oriented Cauchy-Schwarz multi-terminal
//! weighting, 2D column-projected metrics for `sheet3d` overlays, and 3D vertex/metric sizing.

use crate::adapt::metric::{
    elements_for_complexity, from_sizes, grade, identity, intersect, isotropic, map_eig, sizes, Sym,
};
use crate::fem::{build_face_to_cells_3d, facet_rule, geometry, quad_rule, quad_rule_collapsed, shape, CurvedTet, TET_EDGES};
use crate::material::M3;
use crate::mesh::SimplexMesh;
use crate::par;
use std::collections::{HashMap, HashSet};

/// Vacuum permittivity `eps_0` in F/m (matching `nyx-dfs` electrostatics).
pub const EPS0: f64 = 8.854_187_812_8e-12;

/// Default relative floor for off-diagonal capacitance weights `W_ij`.
pub const DEFAULT_EPS_FLOOR: f64 = 1e-3;

/// Per-cell relative permittivity specification (dimensionless).
#[derive(Clone, Copy, Debug)]
pub enum Permittivity<'a> {
    /// Isotropic scalar relative permittivity `eps_r[c] > 0` per cell.
    Scalar(&'a [f64]),
    /// Symmetric positive-definite `3x3` relative permittivity tensor per cell (leading `2x2` used in 2D).
    Tensor(&'a [M3]),
}

impl<'a> Permittivity<'a> {
    fn validate(&self, nc: usize, dim: usize) -> Result<(), String> {
        match *self {
            Permittivity::Scalar(s) => {
                if s.len() != nc {
                    return Err(format!("{} permittivities for {nc} cells", s.len()));
                }
                if !s.iter().all(|&x| x > 0.0 && x.is_finite()) {
                    return Err("scalar relative permittivity must be positive and finite".into());
                }
            }
            Permittivity::Tensor(t) => {
                if t.len() != nc {
                    return Err(format!("{} permittivity tensors for {nc} cells", t.len()));
                }
                for (c, m) in t.iter().enumerate() {
                    if !m.iter().flatten().all(|x| x.is_finite()) {
                        return Err(format!("non-finite permittivity tensor on cell {c}"));
                    }
                    for k in 0..dim {
                        if !(m[k][k] > 0.0) {
                            return Err(format!("non-positive diagonal permittivity on cell {c}"));
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[inline]
    fn matrix(&self, c: usize) -> M3 {
        match *self {
            Permittivity::Scalar(s) => {
                let e = s[c];
                [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]]
            }
            Permittivity::Tensor(t) => t[c],
        }
    }

    #[inline]
    fn mean_scalar(&self, c: usize, dim: usize) -> f64 {
        match *self {
            Permittivity::Scalar(s) => s[c],
            Permittivity::Tensor(t) => ((0..dim).map(|k| t[c][k][k]).sum::<f64>() / dim as f64).max(1e-300),
        }
    }

    #[inline]
    fn normal_scalar(&self, c: usize, n: [f64; 3], dim: usize) -> f64 {
        match *self {
            Permittivity::Scalar(s) => s[c],
            Permittivity::Tensor(t) => {
                let m = &t[c];
                let mut s = 0.0;
                for i in 0..dim {
                    for j in 0..dim {
                        s += n[i] * m[i][j] * n[j];
                    }
                }
                s.max(1e-300)
            }
        }
    }
}

/// Options for [`electro_term_indicators`] and [`electro_goal_indicators`].
pub struct ElectroEstimatorOptions<'a> {
    /// Lagrange polynomial order (`1` or `2`).
    pub order: usize,
    /// Vacuum permittivity scale `eps_0` in F/m (default [`EPS0`]).
    pub eps0: f64,
    /// Per-cell scalar or tensor relative permittivity.
    pub permittivity: Permittivity<'a>,
    /// Mask of length `mesh.n_facets()` marking Dirichlet and conductive sheet facets
    /// (`Terminal`, `Ground`, `Potential`), which are excluded from interior flux jumps and
    /// exterior natural/Robin boundary residuals.
    pub dirichlet_facets: &'a [bool],
    /// Optional per-facet Robin weight `alpha_f >= 0` (in `F/m^2`, length `mesh.n_facets()`).
    pub robin_alpha: Option<&'a [f64]>,
    /// Whether exterior non-Dirichlet facets with `alpha_f == 0` contribute homogeneous Neumann
    /// flux residuals `||n . (eps grad V_h)||^2` (default `true`).
    pub include_neumann_exterior: bool,
    /// Optional predicate skipping cells of inactive/absorbing domain labels.
    pub skip_domain: Option<&'a (dyn Fn(i32) -> bool + Sync)>,
}

impl<'a> ElectroEstimatorOptions<'a> {
    /// Constructs default options for `order` (`1` or `2`), `permittivity`, and `dirichlet_facets`.
    pub fn new(order: usize, permittivity: Permittivity<'a>, dirichlet_facets: &'a [bool]) -> Self {
        ElectroEstimatorOptions {
            order,
            eps0: EPS0,
            permittivity,
            dirichlet_facets,
            robin_alpha: None,
            include_neumann_exterior: true,
            skip_domain: None,
        }
    }
}

/// Per-terminal cell error indicators `eta_{T,i}^2` (in Farads for 1 V terminal excitations).
#[derive(Clone, Debug)]
pub struct ElectroTermIndicators {
    /// `eta2_by_term[i][c]` is `eta_{T_c, i}^2` for terminal `i` and cell `c`.
    pub eta2_by_term: Vec<Vec<f64>>,
    /// Sum over all cells of the volume divergence term per terminal.
    pub vol_part: Vec<f64>,
    /// Sum over all cells of the interior dielectric flux-jump term per terminal.
    pub jump_part: Vec<f64>,
    /// Sum over all cells of the exterior Robin/Neumann boundary term per terminal.
    pub robin_part: Vec<f64>,
}

/// Combined goal-oriented cell error indicators `eta_T^2` across all terminal pairs.
#[derive(Clone, Debug)]
pub struct ElectroGoalIndicators {
    /// Combined dimensionless goal-oriented indicator `eta_T^2` per cell.
    pub eta2: Vec<f64>,
    /// Total indicator `sqrt(sum_T eta_T^2)`.
    pub total: f64,
    /// Per-terminal breakdown.
    pub terms: ElectroTermIndicators,
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn pt3(mesh: &SimplexMesh, v: u32) -> [f64; 3] {
    let p = mesh.point(v as usize);
    [p[0], p[1], p[2]]
}

#[inline]
fn pt2(mesh: &SimplexMesh, v: u32) -> [f64; 2] {
    let p = mesh.point(v as usize);
    [p[0], p[1]]
}

/// Kunert minimum altitude `h_{min,T} = 3 * V / max_face_area` of a 3D tetrahedron.
pub fn kunert_h_min_tet(p: &[[f64; 3]; 4], vol: f64) -> f64 {
    let mut max_area = 0.0_f64;
    for [i, j, k] in [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] {
        let a = 0.5 * norm3(cross3(sub3(p[j], p[i]), sub3(p[k], p[i])));
        max_area = max_area.max(a);
    }
    if max_area > 1e-300 { 3.0 * vol.abs() / max_area } else { 0.0 }
}

/// Kunert minimum altitude `h_{min,T} = 2 * A / max_edge_len` of a 2D triangle.
pub fn kunert_h_min_tri(p: &[[f64; 2]; 3], area: f64) -> f64 {
    let mut max_len = 0.0_f64;
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        max_len = max_len.max((p[j][0] - p[i][0]).hypot(p[j][1] - p[i][1]));
    }
    if max_len > 1e-300 { 2.0 * area.abs() / max_len } else { 0.0 }
}

/// Constant reference Hessian `d^2 f / d xi_k d xi_m` (`3x3`, with `xi_k = lambda_{k+1}`) of a P2
/// ten-node function `vals[0..10]` on the reference tetrahedron.
fn p2_ref_hessian_tet(vals: &[f64; 10]) -> [[f64; 3]; 3] {
    let mut b = [[0.0_f64; 4]; 4];
    for a in 0..4 {
        b[a][a] += 4.0 * vals[a];
    }
    for (e, &(i, j)) in TET_EDGES.iter().enumerate() {
        let w = 4.0 * vals[4 + e];
        b[i][j] += w;
        b[j][i] += w;
    }
    let mut h = [[0.0_f64; 3]; 3];
    for k in 0..3 {
        for m in 0..3 {
            h[k][m] = b[k + 1][m + 1] - b[k + 1][0] - b[0][m + 1] + b[0][0];
        }
    }
    h
}

/// Computes per-terminal dual-weighted residual indicators `eta_{T,i}^2` on `mesh`.
///
/// Parallel volume and face loops reduce in deterministic cell and sorted-face order so outputs
/// are bit-identical across all `NYX_THREADS` settings.
pub fn electro_term_indicators(
    mesh: &SimplexMesh,
    potentials: &[Vec<f64>],
    opt: &ElectroEstimatorOptions,
) -> Result<ElectroTermIndicators, String> {
    if mesh.dim != 2 && mesh.dim != 3 {
        return Err(format!("unsupported mesh dimension {}", mesh.dim));
    }
    if opt.order != 1 && opt.order != 2 {
        return Err(format!("electro_term_indicators supports order 1 or 2, got {}", opt.order));
    }
    if !(opt.eps0 > 0.0 && opt.eps0.is_finite()) {
        return Err("eps0 must be positive and finite".into());
    }
    let nc = mesh.n_cells();
    let nf = mesh.n_facets();
    opt.permittivity.validate(nc, mesh.dim)?;
    if opt.dirichlet_facets.len() != nf {
        return Err(format!("dirichlet_facets mask length {} != n_facets {nf}", opt.dirichlet_facets.len()));
    }
    if let Some(ra) = opt.robin_alpha {
        if ra.len() != nf {
            return Err(format!("robin_alpha length {} != n_facets {nf}", ra.len()));
        }
        if !ra.iter().all(|&x| x >= 0.0 && x.is_finite()) {
            return Err("robin_alpha weights must be non-negative and finite".into());
        }
    }
    let ndofs = mesh.n_dofs(opt.order);
    let n_term = potentials.len();
    for (ti, v) in potentials.iter().enumerate() {
        if v.len() != ndofs {
            return Err(format!("potential[{ti}] has {} dofs, expected {ndofs}", v.len()));
        }
    }
    if n_term == 0 || nc == 0 {
        return Ok(ElectroTermIndicators {
            eta2_by_term: vec![vec![0.0; nc]; n_term],
            vol_part: vec![0.0; n_term],
            jump_part: vec![0.0; n_term],
            robin_part: vec![0.0; n_term],
        });
    }

    if mesh.dim == 2 {
        return electro_term_indicators_2d(mesh, potentials, opt);
    }
    electro_term_indicators_3d(mesh, potentials, opt)
}

fn electro_term_indicators_3d(
    mesh: &SimplexMesh,
    potentials: &[Vec<f64>],
    opt: &ElectroEstimatorOptions,
) -> Result<ElectroTermIndicators, String> {
    let nc = mesh.n_cells();
    let n_term = potentials.len();
    let order = opt.order;
    let nloc = mesh.nloc(order);
    let geo = geometry(mesh);
    let cquad = quad_rule_collapsed(6);
    let cshapes_dn: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let live_cell = |c: usize| -> bool { opt.skip_domain.is_none_or(|sk| !sk(mesh.cell_labels[c])) };

    // 1. Precompute Kunert minimum cell heights and cell volumes (curved-aware).
    let cell_idx: Vec<usize> = (0..nc).collect();
    let vol_chunks = par::map_chunks(&cell_idx, 256, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let mut dofs = Vec::with_capacity(nloc);
        for &c in slice {
            let v = mesh.cell(c);
            let p4 = [pt3(mesh, v[0]), pt3(mesh, v[1]), pt3(mesh, v[2]), pt3(mesh, v[3])];
            let straight_vol = geo.vol[c];
            let curved = mesh.cell_is_curved(c);
            let vol = if curved { CurvedTet::new(mesh.cell_geometry_nodes(c)).volume() } else { straight_vol };
            let h_min = kunert_h_min_tet(&p4, vol);
            let mut terms = vec![0.0_f64; n_term];
            if live_cell(c) && vol > 1e-300 && h_min > 0.0 {
                let eps_m = opt.permittivity.matrix(c);
                let eps_bar = opt.permittivity.mean_scalar(c, 3);
                mesh.cell_dofs(c, order, &mut dofs);
                if curved {
                    let x_nodes = mesh.cell_geometry_nodes(c);
                    let ct = CurvedTet::new(x_nodes);
                    let det0 = CurvedTet::straight(p4).jacobian(&[0.25, 0.25, 0.25, 0.25]).1;
                    let mut hx = [[[0.0_f64; 3]; 3]; 3];
                    for coord in 0..3 {
                        let vals: [f64; 10] = std::array::from_fn(|a| x_nodes[a][coord]);
                        hx[coord] = p2_ref_hessian_tet(&vals);
                    }
                    let mut u_ref_hess = Vec::with_capacity(n_term);
                    for pv in potentials {
                        if order == 2 {
                            let u10: [f64; 10] = std::array::from_fn(|a| pv[dofs[a]]);
                            u_ref_hess.push(p2_ref_hessian_tet(&u10));
                        } else {
                            u_ref_hess.push([[0.0; 3]; 3]);
                        }
                    }
                    for q in 0..cquad.w.len() {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        let (gq, det) = ct.jacobian(&lam);
                        let det_abs = det.abs();
                        if !(det_abs > 1e-14 * det0.abs() && det_abs.is_finite()) {
                            continue;
                        }
                        let det_ratio = if det0.abs() > 1e-300 { det_abs / det0.abs() } else { 1.0 };
                        // Guard Barsoum quarter-point cells (`det_ratio < 0.05` near singular rim)
                        // with local Jacobian size scaling while keeping smooth curved cells (`>= 0.10`) exact.
                        let h_loc = if det_ratio < 0.05 { h_min * (det_ratio / 0.05) } else { h_min };
                        let wq = cquad.w[q] * det_abs / 6.0;
                        let dn = &cshapes_dn[q];
                        for ti in 0..n_term {
                            let pv = &potentials[ti];
                            let mut grad_v = [0.0_f64; 3];
                            for a in 0..nloc {
                                let ua = pv[dofs[a]];
                                for k in 0..3 {
                                    let mut s = 0.0;
                                    for b in 0..4 {
                                        s += dn[a * 4 + b] * gq[b][k];
                                    }
                                    grad_v[k] += ua * s;
                                }
                            }
                            let mut m_ref = u_ref_hess[ti];
                            for k in 0..3 {
                                for m in 0..3 {
                                    m_ref[k][m] -= grad_v[0] * hx[0][k][m]
                                        + grad_v[1] * hx[1][k][m]
                                        + grad_v[2] * hx[2][k][m];
                                }
                            }
                            // Hess_x V = R^T m_ref R, where R[k][i] = gq[k + 1][i] = (J^-1)_{k,i}.
                            let mut div_d = 0.0_f64;
                            for i in 0..3 {
                                for j in 0..3 {
                                    let eij = eps_m[i][j];
                                    if eij == 0.0 {
                                        continue;
                                    }
                                    let mut hij = 0.0;
                                    for k in 0..3 {
                                        let rki = gq[k + 1][i];
                                        for m in 0..3 {
                                            hij += rki * m_ref[k][m] * gq[m + 1][j];
                                        }
                                    }
                                    div_d += eij * hij;
                                }
                            }
                            terms[ti] += (h_loc * h_loc * opt.eps0 / eps_bar) * wq * div_d * div_d;
                        }
                    }
                } else if order == 2 {
                    // On straight P2 tetrahedra, Hess_x V is constant across the cell.
                    let base = c * 12;
                    let g = |b: usize, k: usize| geo.g[base + b * 3 + k];
                    for (ti, pv) in potentials.iter().enumerate() {
                        let u10: [f64; 10] = std::array::from_fn(|a| pv[dofs[a]]);
                        let h_ref = p2_ref_hessian_tet(&u10);
                        let mut div_d = 0.0_f64;
                        for i in 0..3 {
                            for j in 0..3 {
                                let eij = eps_m[i][j];
                                if eij == 0.0 {
                                    continue;
                                }
                                let mut hij = 0.0;
                                for k in 0..3 {
                                    let rki = g(k + 1, i);
                                    for m in 0..3 {
                                        hij += rki * h_ref[k][m] * g(m + 1, j);
                                    }
                                }
                                div_d += eij * hij;
                            }
                        }
                        terms[ti] = (h_min * h_min * opt.eps0 / eps_bar) * vol * div_d * div_d;
                    }
                }
            }
            out.push((vol, terms));
        }
        out
    });

    let mut cell_vols = Vec::with_capacity(nc);
    let mut eta2_by_term = vec![vec![0.0_f64; nc]; n_term];
    let mut vol_part = vec![0.0_f64; n_term];
    let mut cid = 0usize;
    for chunk in vol_chunks {
        for (vol, terms) in chunk {
            cell_vols.push(vol);
            for ti in 0..n_term {
                eta2_by_term[ti][cid] += terms[ti];
                vol_part[ti] += terms[ti];
            }
            cid += 1;
        }
    }

    // 2. Build Dirichlet/sheet face set and Robin alpha map keyed by sorted vertex triple.
    let mut dirichlet_keys: HashSet<[u32; 3]> = HashSet::new();
    let mut robin_map: HashMap<[u32; 3], f64> = HashMap::new();
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let mut key = [fv[0], fv[1], fv[2]];
        key.sort_unstable();
        if opt.dirichlet_facets[f] {
            dirichlet_keys.insert(key);
        } else if let Some(ra) = opt.robin_alpha {
            let entry = robin_map.entry(key).or_insert(0.0);
            *entry = entry.max(ra[f]);
        }
    }

    // 3. Sorted unique mesh faces for deterministic parallel jump & Robin integration.
    let face_to_cells = build_face_to_cells_3d(mesh);
    let mut keys: Vec<[u32; 3]> = face_to_cells.keys().copied().collect();
    keys.sort_unstable();

    let face_rule = quad_rule(2, 4);
    let face_shapes: Vec<(Vec<f64>, Vec<f64>)> = face_rule
        .lam
        .iter()
        .map(|l2| {
            let l4 = [l2[0], l2[1], l2[2], 0.0];
            let _ = l4;
            shape(l2, order)
        })
        .collect();
    let _ = face_shapes;

    struct FaceContrib {
        c1: usize,
        c2: Option<usize>,
        vals: Vec<f64>,
    }

    let eval_grad_v = |c: usize,
                       ct: Option<&CurvedTet>,
                       lam4: &[f64; 4],
                       dn4: &[f64],
                       dofs: &[usize],
                       pv: &[f64]|
     -> [f64; 3] {
        let g_cell: [[f64; 3]; 4] = if let Some(ct_ref) = ct {
            ct_ref.jacobian(lam4).0
        } else {
            let base = c * 12;
            std::array::from_fn(|b| [geo.g[base + b * 3], geo.g[base + b * 3 + 1], geo.g[base + b * 3 + 2]])
        };
        let mut gv = [0.0_f64; 3];
        for a in 0..nloc {
            let ua = pv[dofs[a]];
            for k in 0..3 {
                let mut s = 0.0;
                for b in 0..4 {
                    s += dn4[a * 4 + b] * g_cell[b][k];
                }
                gv[k] += ua * s;
            }
        }
        gv
    };

    let face_chunks = par::map_chunks(&keys, 256, |kslice| {
        let mut out: Vec<FaceContrib> = Vec::with_capacity(kslice.len());
        let mut dofs1 = Vec::with_capacity(nloc);
        let mut dofs2 = Vec::with_capacity(nloc);
        for &key in kslice {
            if dirichlet_keys.contains(&key) {
                continue;
            }
            let cs = &face_to_cells[&key];
            let (p0, p1, p2) = (pt3(mesh, key[0]), pt3(mesh, key[1]), pt3(mesh, key[2]));
            let n_straight_raw = cross3(sub3(p1, p0), sub3(p2, p0));
            let straight_len = norm3(n_straight_raw);
            let straight_area = 0.5 * straight_len;
            if straight_area <= 1e-300 {
                continue;
            }
            let (c1, lf1) = cs[0];
            let opp1 = 6 - lf1[0] - lf1[1] - lf1[2];
            let p_opp1 = pt3(mesh, mesh.cell(c1)[opp1]);
            let sign1 = if dot3(n_straight_raw, sub3(p_opp1, p0)) > 0.0 { -1.0 } else { 1.0 };
            let n_straight = [
                sign1 * n_straight_raw[0] / straight_len,
                sign1 * n_straight_raw[1] / straight_len,
                sign1 * n_straight_raw[2] / straight_len,
            ];

            if cs.len() == 2 {
                let (c2, lf2) = cs[1];
                if !live_cell(c1) || !live_cell(c2) {
                    continue;
                }
                let ct1 = mesh.cell_is_curved(c1).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c1)));
                let ct2 = mesh.cell_is_curved(c2).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c2)));
                mesh.cell_dofs(c1, order, &mut dofs1);
                mesh.cell_dofs(c2, order, &mut dofs2);
                let eps1 = opt.permittivity.matrix(c1);
                let eps2 = opt.permittivity.matrix(c2);

                let mut area_acc = 0.0_f64;
                let mut jump_int = vec![0.0_f64; n_term];
                for (q, l2) in face_rule.lam.iter().enumerate() {
                    let mut lam1 = [0.0_f64; 4];
                    let mut lam2 = [0.0_f64; 4];
                    for m in 0..3 {
                        lam1[lf1[m]] = l2[m];
                        lam2[lf2[m]] = l2[m];
                    }
                    let (_, dn1) = shape(&lam1, order);
                    let (_, dn2) = shape(&lam2, order);
                    let (n_q, ds_q) = if let Some(ref ct) = ct1 {
                        let d = ct.dx_dlam(&lam1);
                        let ds = sub3(d[lf1[1]], d[lf1[0]]);
                        let dt = sub3(d[lf1[2]], d[lf1[0]]);
                        let cr = cross3(ds, dt);
                        let cl = norm3(cr);
                        if cl > 1e-300 {
                            let sgn = if dot3(cr, n_straight) < 0.0 { -1.0 } else { 1.0 };
                            ([sgn * cr[0] / cl, sgn * cr[1] / cl, sgn * cr[2] / cl], face_rule.w[q] * 0.5 * cl)
                        } else {
                            (n_straight, face_rule.w[q] * straight_area)
                        }
                    } else {
                        (n_straight, face_rule.w[q] * straight_area)
                    };
                    area_acc += ds_q;
                    for ti in 0..n_term {
                        let gv1 = eval_grad_v(c1, ct1.as_ref(), &lam1, &dn1, &dofs1, &potentials[ti]);
                        let gv2 = eval_grad_v(c2, ct2.as_ref(), &lam2, &dn2, &dofs2, &potentials[ti]);
                        let mut flux1 = 0.0;
                        let mut flux2 = 0.0;
                        for i in 0..3 {
                            let mut d1i = 0.0;
                            let mut d2i = 0.0;
                            for j in 0..3 {
                                d1i += eps1[i][j] * gv1[j];
                                d2i += eps2[i][j] * gv2[j];
                            }
                            flux1 += n_q[i] * d1i;
                            flux2 += n_q[i] * d2i;
                        }
                        let jmp = flux1 - flux2;
                        jump_int[ti] += ds_q * jmp * jmp;
                    }
                }
                let area_f = area_acc.max(1e-300);
                let h_f = 3.0 * (cell_vols[c1] + cell_vols[c2]) / (2.0 * area_f);
                let eps_bar_f = 0.5
                    * (opt.permittivity.normal_scalar(c1, n_straight, 3)
                        + opt.permittivity.normal_scalar(c2, n_straight, 3));
                let scale = 0.5 * h_f * opt.eps0 / eps_bar_f.max(1e-300);
                let vals: Vec<f64> = jump_int.into_iter().map(|v| scale * v).collect();
                out.push(FaceContrib { c1, c2: Some(c2), vals });
            } else if cs.len() == 1 && live_cell(c1) {
                let alpha_f = robin_map.get(&key).copied().unwrap_or(0.0);
                if alpha_f == 0.0 && !opt.include_neumann_exterior {
                    continue;
                }
                let ct1 = mesh.cell_is_curved(c1).then(|| CurvedTet::new(mesh.cell_geometry_nodes(c1)));
                mesh.cell_dofs(c1, order, &mut dofs1);
                let eps1 = opt.permittivity.matrix(c1);
                let alpha_rel = alpha_f / opt.eps0;

                let mut area_acc = 0.0_f64;
                let mut rob_int = vec![0.0_f64; n_term];
                for (q, l2) in face_rule.lam.iter().enumerate() {
                    let mut lam1 = [0.0_f64; 4];
                    for m in 0..3 {
                        lam1[lf1[m]] = l2[m];
                    }
                    let (nv1, dn1) = shape(&lam1, order);
                    let (n_q, ds_q) = if let Some(ref ct) = ct1 {
                        let d = ct.dx_dlam(&lam1);
                        let ds = sub3(d[lf1[1]], d[lf1[0]]);
                        let dt = sub3(d[lf1[2]], d[lf1[0]]);
                        let cr = cross3(ds, dt);
                        let cl = norm3(cr);
                        if cl > 1e-300 {
                            let sgn = if dot3(cr, n_straight) < 0.0 { -1.0 } else { 1.0 };
                            ([sgn * cr[0] / cl, sgn * cr[1] / cl, sgn * cr[2] / cl], face_rule.w[q] * 0.5 * cl)
                        } else {
                            (n_straight, face_rule.w[q] * straight_area)
                        }
                    } else {
                        (n_straight, face_rule.w[q] * straight_area)
                    };
                    area_acc += ds_q;
                    for ti in 0..n_term {
                        let pv = &potentials[ti];
                        let gv1 = eval_grad_v(c1, ct1.as_ref(), &lam1, &dn1, &dofs1, pv);
                        let mut flux1 = 0.0;
                        for i in 0..3 {
                            let mut d1i = 0.0;
                            for j in 0..3 {
                                d1i += eps1[i][j] * gv1[j];
                            }
                            flux1 += n_q[i] * d1i;
                        }
                        let mut u_q = 0.0;
                        if alpha_rel > 0.0 {
                            for a in 0..nloc {
                                u_q += nv1[a] * pv[dofs1[a]];
                            }
                        }
                        let r = flux1 + alpha_rel * u_q;
                        rob_int[ti] += ds_q * r * r;
                    }
                }
                let area_f = area_acc.max(1e-300);
                let h_f = 3.0 * cell_vols[c1] / area_f;
                let eps_bar_t = opt.permittivity.mean_scalar(c1, 3);
                let scale = h_f * opt.eps0 / eps_bar_t.max(1e-300);
                let vals: Vec<f64> = rob_int.into_iter().map(|v| scale * v).collect();
                out.push(FaceContrib { c1, c2: None, vals });
            }
        }
        out
    });

    let mut jump_part = vec![0.0_f64; n_term];
    let mut robin_part = vec![0.0_f64; n_term];
    for chunk in face_chunks {
        for fc in chunk {
            match fc.c2 {
                Some(c2) => {
                    for ti in 0..n_term {
                        let v = fc.vals[ti];
                        eta2_by_term[ti][fc.c1] += v;
                        eta2_by_term[ti][c2] += v;
                        jump_part[ti] += 2.0 * v;
                    }
                }
                None => {
                    for ti in 0..n_term {
                        let v = fc.vals[ti];
                        eta2_by_term[ti][fc.c1] += v;
                        robin_part[ti] += v;
                    }
                }
            }
        }
    }

    Ok(ElectroTermIndicators { eta2_by_term, vol_part, jump_part, robin_part })
}

fn electro_term_indicators_2d(
    mesh: &SimplexMesh,
    potentials: &[Vec<f64>],
    opt: &ElectroEstimatorOptions,
) -> Result<ElectroTermIndicators, String> {
    let nc = mesh.n_cells();
    let n_term = potentials.len();
    let order = opt.order;
    let nloc = mesh.nloc(order);
    let geo = geometry(mesh);
    let live_cell = |c: usize| -> bool { opt.skip_domain.is_none_or(|sk| !sk(mesh.cell_labels[c])) };

    let mut eta2_by_term = vec![vec![0.0_f64; nc]; n_term];
    let mut vol_part = vec![0.0_f64; n_term];
    let mut dofs = Vec::with_capacity(nloc);
    for c in 0..nc {
        if !live_cell(c) || order < 2 {
            continue;
        }
        let v = mesh.cell(c);
        let p3 = [pt2(mesh, v[0]), pt2(mesh, v[1]), pt2(mesh, v[2])];
        let area = geo.vol[c];
        let h_min = kunert_h_min_tri(&p3, area);
        if area <= 1e-300 || h_min <= 0.0 {
            continue;
        }
        let eps_m = opt.permittivity.matrix(c);
        let eps_bar = opt.permittivity.mean_scalar(c, 2);
        mesh.cell_dofs(c, 2, &mut dofs);
        let base = c * 6;
        let g = |b: usize, k: usize| geo.g[base + b * 2 + k];
        for (ti, pv) in potentials.iter().enumerate() {
            // P2 on triangle: vertices 0,1,2, edges (0,1)=3, (0,2)=4, (1,2)=5
            let mut bmat = [[0.0_f64; 3]; 3];
            for a in 0..3 {
                bmat[a][a] += 4.0 * pv[dofs[a]];
            }
            for (e, &(i, j)) in [(0usize, 1usize), (0, 2), (1, 2)].iter().enumerate() {
                let w = 4.0 * pv[dofs[3 + e]];
                bmat[i][j] += w;
                bmat[j][i] += w;
            }
            let mut h_ref = [[0.0_f64; 2]; 2];
            for k in 0..2 {
                for m in 0..2 {
                    h_ref[k][m] = bmat[k + 1][m + 1] - bmat[k + 1][0] - bmat[0][m + 1] + bmat[0][0];
                }
            }
            let mut div_d = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    let mut hij = 0.0;
                    for k in 0..2 {
                        for m in 0..2 {
                            hij += g(k + 1, i) * h_ref[k][m] * g(m + 1, j);
                        }
                    }
                    div_d += eps_m[i][j] * hij;
                }
            }
            let val = (h_min * h_min * opt.eps0 / eps_bar) * area * div_d * div_d;
            eta2_by_term[ti][c] += val;
            vol_part[ti] += val;
        }
    }

    let mut dirichlet_keys: HashSet<[u32; 2]> = HashSet::new();
    let mut robin_map: HashMap<[u32; 2], f64> = HashMap::new();
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        let key = [fv[0].min(fv[1]), fv[0].max(fv[1])];
        if opt.dirichlet_facets[f] {
            dirichlet_keys.insert(key);
        } else if let Some(ra) = opt.robin_alpha {
            let entry = robin_map.entry(key).or_insert(0.0);
            *entry = entry.max(ra[f]);
        }
    }

    let mut edge_to_cells: HashMap<[u32; 2], Vec<(usize, [usize; 2])>> = HashMap::new();
    for ci in 0..nc {
        let v = mesh.cell(ci);
        for mut le in [[1usize, 2], [0, 2], [0, 1]] {
            if v[le[0]] > v[le[1]] {
                le.swap(0, 1);
            }
            let key = [v[le[0]], v[le[1]]];
            edge_to_cells.entry(key).or_default().push((ci, le));
        }
    }
    let mut keys: Vec<[u32; 2]> = edge_to_cells.keys().copied().collect();
    keys.sort_unstable();
    let erule = facet_rule(1);
    let mut dofs1 = Vec::with_capacity(nloc);
    let mut dofs2 = Vec::with_capacity(nloc);
    let mut jump_part = vec![0.0_f64; n_term];
    let mut robin_part = vec![0.0_f64; n_term];

    for key in keys {
        if dirichlet_keys.contains(&key) {
            continue;
        }
        let cs = &edge_to_cells[&key];
        let (p0, p1) = (pt2(mesh, key[0]), pt2(mesh, key[1]));
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let len = dx.hypot(dy);
        if len <= 1e-300 {
            continue;
        }
        let (c1, le1) = cs[0];
        let opp1 = 3 - le1[0] - le1[1];
        let p_opp1 = pt2(mesh, mesh.cell(c1)[opp1]);
        let mut n = [dy / len, -dx / len];
        if n[0] * (p_opp1[0] - p0[0]) + n[1] * (p_opp1[1] - p0[1]) > 0.0 {
            n = [-n[0], -n[1]];
        }
        let n3 = [n[0], n[1], 0.0];

        if cs.len() == 2 {
            let (c2, le2) = cs[1];
            if !live_cell(c1) || !live_cell(c2) {
                continue;
            }
            mesh.cell_dofs(c1, order, &mut dofs1);
            mesh.cell_dofs(c2, order, &mut dofs2);
            let eps1 = opt.permittivity.matrix(c1);
            let eps2 = opt.permittivity.matrix(c2);
            let h_f = (geo.vol[c1] + geo.vol[c2]) / len;
            let eps_bar_f = 0.5
                * (opt.permittivity.normal_scalar(c1, n3, 2) + opt.permittivity.normal_scalar(c2, n3, 2));
            let scale = 0.5 * h_f * opt.eps0 / eps_bar_f.max(1e-300);

            for (q, l1) in erule.lam.iter().enumerate() {
                let mut lam1 = [0.0; 3];
                let mut lam2 = [0.0; 3];
                lam1[le1[0]] = l1[0];
                lam1[le1[1]] = l1[1];
                lam2[le2[0]] = l1[0];
                lam2[le2[1]] = l1[1];
                let (_, dn1) = shape(&lam1, order);
                let (_, dn2) = shape(&lam2, order);
                let ds_q = erule.w[q] * len;
                for ti in 0..n_term {
                    let pv = &potentials[ti];
                    let mut gv1 = [0.0; 2];
                    let mut gv2 = [0.0; 2];
                    for a in 0..nloc {
                        for k in 0..2 {
                            let mut s1 = 0.0;
                            let mut s2 = 0.0;
                            for b in 0..3 {
                                s1 += dn1[a * 3 + b] * geo.g[c1 * 6 + b * 2 + k];
                                s2 += dn2[a * 3 + b] * geo.g[c2 * 6 + b * 2 + k];
                            }
                            gv1[k] += pv[dofs1[a]] * s1;
                            gv2[k] += pv[dofs2[a]] * s2;
                        }
                    }
                    let mut f1 = 0.0;
                    let mut f2 = 0.0;
                    for i in 0..2 {
                        f1 += n[i] * (eps1[i][0] * gv1[0] + eps1[i][1] * gv1[1]);
                        f2 += n[i] * (eps2[i][0] * gv2[0] + eps2[i][1] * gv2[1]);
                    }
                    let contrib = scale * ds_q * (f1 - f2).powi(2);
                    eta2_by_term[ti][c1] += contrib;
                    eta2_by_term[ti][c2] += contrib;
                    jump_part[ti] += 2.0 * contrib;
                }
            }
        } else if cs.len() == 1 && live_cell(c1) {
            let alpha_f = robin_map.get(&key).copied().unwrap_or(0.0);
            if alpha_f == 0.0 && !opt.include_neumann_exterior {
                continue;
            }
            mesh.cell_dofs(c1, order, &mut dofs1);
            let eps1 = opt.permittivity.matrix(c1);
            let alpha_rel = alpha_f / opt.eps0;
            let h_f = 2.0 * geo.vol[c1] / len;
            let eps_bar = opt.permittivity.mean_scalar(c1, 2);
            let scale = h_f * opt.eps0 / eps_bar.max(1e-300);

            for (q, l1) in erule.lam.iter().enumerate() {
                let mut lam1 = [0.0; 3];
                lam1[le1[0]] = l1[0];
                lam1[le1[1]] = l1[1];
                let (nv1, dn1) = shape(&lam1, order);
                let ds_q = erule.w[q] * len;
                for ti in 0..n_term {
                    let pv = &potentials[ti];
                    let mut gv1 = [0.0; 2];
                    let mut u_q = 0.0;
                    for a in 0..nloc {
                        u_q += nv1[a] * pv[dofs1[a]];
                        for k in 0..2 {
                            let mut s1 = 0.0;
                            for b in 0..3 {
                                s1 += dn1[a * 3 + b] * geo.g[c1 * 6 + b * 2 + k];
                            }
                            gv1[k] += pv[dofs1[a]] * s1;
                        }
                    }
                    let mut f1 = 0.0;
                    for i in 0..2 {
                        f1 += n[i] * (eps1[i][0] * gv1[0] + eps1[i][1] * gv1[1]);
                    }
                    let contrib = scale * ds_q * (f1 + alpha_rel * u_q).powi(2);
                    eta2_by_term[ti][c1] += contrib;
                    robin_part[ti] += contrib;
                }
            }
        }
    }

    Ok(ElectroTermIndicators { eta2_by_term, vol_part, jump_part, robin_part })
}

/// Combines per-terminal indicators `eta_{T,i}^2` into goal-oriented relative capacitance error
/// indicators using Cauchy-Schwarz weights:
///
/// `eta_T^2 = sum_{0 <= i <= j < N} W_{ij} eta_{T,i} eta_{T,j}`,
/// `W_{ij} = 1 / max(|C_{ij}|, eps_floor * sqrt(|C_{ii} C_{jj}|))`.
pub fn combine_goal_oriented(
    terms: ElectroTermIndicators,
    cap_matrix: &[Vec<f64>],
    eps_floor: f64,
) -> Result<ElectroGoalIndicators, String> {
    let n_term = terms.eta2_by_term.len();
    if cap_matrix.len() != n_term || cap_matrix.iter().any(|r| r.len() != n_term) {
        return Err(format!("cap_matrix must be {n_term}x{n_term}, got {} rows", cap_matrix.len()));
    }
    if !(eps_floor >= 0.0 && eps_floor.is_finite()) {
        return Err("eps_floor must be non-negative and finite".into());
    }
    let nc = terms.eta2_by_term.first().map_or(0, Vec::len);
    let mut weights: Vec<(usize, usize, f64)> = Vec::with_capacity(n_term * (n_term + 1) / 2);
    for i in 0..n_term {
        let cii = cap_matrix[i][i].abs();
        for j in i..n_term {
            let cjj = cap_matrix[j][j].abs();
            let cij = cap_matrix[i][j].abs();
            let denom = cij.max(eps_floor * (cii * cjj).sqrt()).max(1e-300);
            weights.push((i, j, 1.0 / denom));
        }
    }
    let mut eta2 = vec![0.0_f64; nc];
    for c in 0..nc {
        let mut s = 0.0_f64;
        for &(i, j, wij) in &weights {
            let ei = terms.eta2_by_term[i][c].max(0.0).sqrt();
            let ej = terms.eta2_by_term[j][c].max(0.0).sqrt();
            s += wij * ei * ej;
        }
        eta2[c] = s;
    }
    let total = eta2.iter().sum::<f64>().sqrt();
    Ok(ElectroGoalIndicators { eta2, total, terms })
}

/// Evaluates per-terminal indicators and combines them with `cap_matrix` (`eps_floor` defaults to [`DEFAULT_EPS_FLOOR`]).
pub fn electro_goal_indicators(
    mesh: &SimplexMesh,
    potentials: &[Vec<f64>],
    cap_matrix: &[Vec<f64>],
    eps_floor: f64,
    opt: &ElectroEstimatorOptions,
) -> Result<ElectroGoalIndicators, String> {
    let terms = electro_term_indicators(mesh, potentials, opt)?;
    combine_goal_oriented(terms, cap_matrix, eps_floor)
}

// ------------------------------------------------------------------ 2D column projection & metric sizing

/// Projects 3D cell indicators `eta2_3d` onto the 2D overlay triangles of `m2` along `plane_axes`
/// (`[0, 1]` for standard `(x, y)` `Axis::Z` stacks; `[perm[0], perm[1]]` for general sweep axes).
///
/// Prismatic cells map 1-to-1 to their owning `m2` triangle (`eta_{t,2D}^2 = sum_{T in col(t)} eta_T^2`),
/// while coarsened far-field cells (`z_coarsen`) spanning multiple `m2` triangles distribute their
/// indicator area-proportionally across covered `m2` columns.
pub fn project_indicators_to_overlay_2d(
    mesh3: &SimplexMesh,
    m2: &SimplexMesh,
    eta2_3d: &[f64],
    plane_axes: [usize; 2],
) -> Result<Vec<f64>, String> {
    if mesh3.dim != 3 || m2.dim != 2 {
        return Err("project_indicators_to_overlay_2d requires a 3D mesh and a 2D overlay mesh".into());
    }
    if eta2_3d.len() != mesh3.n_cells() {
        return Err(format!("{} 3D indicators for {} 3D cells", eta2_3d.len(), mesh3.n_cells()));
    }
    let nt2 = m2.n_cells();
    if nt2 == 0 {
        return Ok(Vec::new());
    }
    let [ax, ay] = plane_axes;
    if ax >= 3 || ay >= 3 || ax == ay {
        return Err(format!("invalid plane_axes {plane_axes:?}"));
    }

    let tris2: Vec<[[f64; 2]; 3]> = (0..nt2)
        .map(|t| {
            let v = m2.cell(t);
            [pt2(m2, v[0]), pt2(m2, v[1]), pt2(m2, v[2])]
        })
        .collect();
    let areas2: Vec<f64> = tris2
        .iter()
        .map(|&[a, b, c]| 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs())
        .collect();
    let cens2: Vec<[f64; 2]> = tris2
        .iter()
        .map(|&[a, b, c]| [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0])
        .collect();

    let (lo2, hi2) = m2.bbox();
    let span_x = (hi2[0] - lo2[0]).max(1e-30);
    let span_y = (hi2[1] - lo2[1]).max(1e-30);
    let cell_sz = (2.0 * span_x * span_y / nt2 as f64).sqrt().max(1e-12);
    let nx = ((span_x / cell_sz).ceil() as usize).max(1);
    let ny = ((span_y / cell_sz).ceil() as usize).max(1);
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); nx * ny];
    for (t, &[a, b, c]) in tris2.iter().enumerate() {
        let x0 = a[0].min(b[0]).min(c[0]);
        let x1 = a[0].max(b[0]).max(c[0]);
        let y0 = a[1].min(b[1]).min(c[1]);
        let y1 = a[1].max(b[1]).max(c[1]);
        let ix0 = (((x0 - lo2[0]) / cell_sz).max(0.0) as usize).min(nx - 1);
        let ix1 = (((x1 - lo2[0]) / cell_sz).max(0.0) as usize).min(nx - 1);
        let iy0 = (((y0 - lo2[1]) / cell_sz).max(0.0) as usize).min(ny - 1);
        let iy1 = (((y1 - lo2[1]) / cell_sz).max(0.0) as usize).min(ny - 1);
        for iy in iy0..=iy1 {
            for ix in ix0..=ix1 {
                buckets[iy * nx + ix].push(t);
            }
        }
    }

    let bary2 = |t: usize, p: [f64; 2]| -> [f64; 3] {
        let [a, b, c] = tris2[t];
        let det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if det.abs() <= 1e-300 {
            return [-1.0; 3];
        }
        let l1 = ((p[0] - a[0]) * (c[1] - a[1]) - (p[1] - a[1]) * (c[0] - a[0])) / det;
        let l2 = ((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) / det;
        [1.0 - l1 - l2, l1, l2]
    };

    let locate_tri = |p: [f64; 2]| -> usize {
        let ix = (((p[0] - lo2[0]) / cell_sz).max(0.0) as usize).min(nx - 1);
        let iy = (((p[1] - lo2[1]) / cell_sz).max(0.0) as usize).min(ny - 1);
        let mut best: Option<(f64, usize)> = None;
        for &t in &buckets[iy * nx + ix] {
            let b = bary2(t, p);
            let m = b[0].min(b[1]).min(b[2]);
            if best.is_none_or(|(bm, bt)| m > bm || (m == bm && t < bt)) {
                best = Some((m, t));
            }
        }
        if best.is_none_or(|b| b.0 < -1e-9) {
            for jy in iy.saturating_sub(1)..=(iy + 1).min(ny - 1) {
                for jx in ix.saturating_sub(1)..=(ix + 1).min(nx - 1) {
                    for &t in &buckets[jy * nx + jx] {
                        let b = bary2(t, p);
                        let m = b[0].min(b[1]).min(b[2]);
                        if best.is_none_or(|(bm, bt)| m > bm || (m == bm && t < bt)) {
                            best = Some((m, t));
                        }
                    }
                }
            }
        }
        if let Some((_, t)) = best {
            return t;
        }
        (0..nt2)
            .min_by(|&a, &b| {
                let da = (cens2[a][0] - p[0]).powi(2) + (cens2[a][1] - p[1]).powi(2);
                let db = (cens2[b][0] - p[0]).powi(2) + (cens2[b][1] - p[1]).powi(2);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b))
            })
            .unwrap_or(0)
    };

    let mut eta2_2d = vec![0.0_f64; nt2];
    let mut covered: Vec<usize> = Vec::new();
    for (c, &val) in eta2_3d.iter().enumerate() {
        if val <= 0.0 {
            continue;
        }
        let v = mesh3.cell(c);
        let p: [[f64; 2]; 4] = std::array::from_fn(|i| {
            let q = mesh3.point(v[i] as usize);
            [q[ax], q[ay]]
        });
        let cen = [0.25 * (p[0][0] + p[1][0] + p[2][0] + p[3][0]), 0.25 * (p[0][1] + p[1][1] + p[2][1] + p[3][1])];
        let t_cen = locate_tri(cen);
        let x0 = p[0][0].min(p[1][0]).min(p[2][0]).min(p[3][0]);
        let x1 = p[0][0].max(p[1][0]).max(p[2][0]).max(p[3][0]);
        let y0 = p[0][1].min(p[1][1]).min(p[2][1]).min(p[3][1]);
        let y1 = p[0][1].max(p[1][1]).max(p[2][1]).max(p[3][1]);
        let bbox_area = (x1 - x0) * (y1 - y0);
        // Fast 1-to-1 prismatic column path; only check multi-triangle distribution when a coarsened tet spans > 2.5x area.
        if bbox_area <= 2.5 * areas2[t_cen] {
            eta2_2d[t_cen] += val;
            continue;
        }
        covered.clear();
        let ix0 = (((x0 - lo2[0]) / cell_sz).max(0.0) as usize).min(nx - 1);
        let ix1 = (((x1 - lo2[0]) / cell_sz).max(0.0) as usize).min(nx - 1);
        let iy0 = (((y0 - lo2[1]) / cell_sz).max(0.0) as usize).min(ny - 1);
        let iy1 = (((y1 - lo2[1]) / cell_sz).max(0.0) as usize).min(ny - 1);
        for iy in iy0..=iy1 {
            for ix in ix0..=ix1 {
                for &t in &buckets[iy * nx + ix] {
                    let tc = cens2[t];
                    if tc[0] >= x0 && tc[0] <= x1 && tc[1] >= y0 && tc[1] <= y1 && pt_in_tet_proj2(&p, tc) {
                        covered.push(t);
                    }
                }
            }
        }
        covered.sort_unstable();
        covered.dedup();
        let area_sum: f64 = covered.iter().map(|&t| areas2[t]).sum();
        if covered.len() <= 1 || area_sum <= 1e-300 {
            eta2_2d[t_cen] += val;
        } else {
            for &t in &covered {
                eta2_2d[t] += val * (areas2[t] / area_sum);
            }
        }
    }
    Ok(eta2_2d)
}

fn pt_in_tet_proj2(p: &[[f64; 2]; 4], q: [f64; 2]) -> bool {
    let in_tri = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| -> bool {
        let d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if d.abs() <= 1e-300 {
            return false;
        }
        let l1 = ((q[0] - a[0]) * (c[1] - a[1]) - (q[1] - a[1]) * (c[0] - a[0])) / d;
        let l2 = ((b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])) / d;
        let l0 = 1.0 - l1 - l2;
        l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9
    };
    in_tri(p[0], p[1], p[2]) || in_tri(p[0], p[1], p[3]) || in_tri(p[0], p[2], p[3]) || in_tri(p[1], p[2], p[3])
}

/// Options for 2D column-projected metric generation [`electro_overlay_metrics_2d`].
#[derive(Clone, Debug)]
pub struct ElectroMetric2dOptions<'a> {
    /// Polynomial order `p` (`1` or `2`); smooth convergence exponent is `q = (p + 1) as f64`.
    pub order: usize,
    /// Singular convergence rate near `sheet:edge` rims (default `1.0`).
    pub q_sing: f64,
    /// Largest allowed element size.
    pub h_max: f64,
    /// Smallest allowed element size as a fraction of `h_max` (`h_min = h_max * h_min_fraction`).
    pub h_min_fraction: f64,
    /// Target number of 2D triangles in `m2`.
    pub target_cells: f64,
    /// Maximum allowed directional aspect ratio along straight sheet rims (default `20.0`).
    pub max_aspect: f64,
    /// Requested tangential-to-normal aspect ratio on straight sheet rims away from corners (`1.0..=max_aspect`).
    pub rim_anisotropy: f64,
    /// 2D segments `(p0, p1)` of `sheet:edge` rims in `m2` coordinates.
    pub rim_segments: &'a [([f64; 2], [f64; 2])],
    /// Sharp `90 deg / 270 deg` 2D sheet corner coordinates where isotropic corner-ball intersection is enforced
    /// BEFORE `BackgroundMetric2::new_graded` so orthogonal anisotropic layers never collide.
    pub corners: &'a [[f64; 2]],
    /// Radius multiplier (in units of local triangle size) of the isotropic corner ball (default `2.5`).
    pub corner_ball_factor: f64,
    /// Metric gradation slope `g` accounted for during complexity targeting (default `0.35`).
    pub grade: f64,
}

impl<'a> ElectroMetric2dOptions<'a> {
    pub fn new(order: usize, h_max: f64, h_min_fraction: f64, target_cells: f64) -> Self {
        ElectroMetric2dOptions {
            order,
            q_sing: 1.0,
            h_max,
            h_min_fraction,
            target_cells,
            max_aspect: 20.0,
            rim_anisotropy: 1.0,
            rim_segments: &[],
            corners: &[],
            corner_ball_factor: 2.5,
            grade: 0.35,
        }
    }
}

fn closest_rim_seg(p: [f64; 2], segs: &[([f64; 2], [f64; 2])]) -> Option<(f64, [f64; 2])> {
    let mut best: Option<(f64, [f64; 2])> = None;
    for &(a, b) in segs {
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let l2 = dx * dx + dy * dy;
        if l2 <= 1e-300 {
            continue;
        }
        let t = (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0);
        let q = [a[0] + t * dx, a[1] + t * dy];
        let d = (p[0] - q[0]).hypot(p[1] - q[1]);
        if best.is_none_or(|(bd, _)| d < bd) {
            let len = l2.sqrt();
            // Unit normal [-dy/len, dx/len] perpendicular to rim segment
            best = Some((d, [-dy / len, dx / len]));
        }
    }
    best
}

/// Converts column-projected 2D triangle indicators `eta2_2d` into per-vertex 2D Riemannian metrics
/// `Vec<Sym<2>>` on `m2`:
/// - Uses smooth exponent `q = p + 1` away from singularities and `q_sing` (default `1.0`) near `rim_segments`.
/// - Builds directional metrics via [`from_sizes`] capped at `max_aspect`.
/// - Intersects an isotropic corner-ball metric around every sharp corner in `corners` *before*
///   [`crate::adapt::field::BackgroundMetric2::new_graded`] so orthogonal anisotropic boundary layers
///   never collide at `90 deg / 270 deg` corners.
/// - Clamps all sizes strictly to `[h_max * h_min_fraction, h_max]`.
pub fn electro_overlay_metrics_2d(
    m2: &SimplexMesh,
    eta2_2d: &[f64],
    opt: &ElectroMetric2dOptions,
    floor: Option<&[Sym<2>]>,
) -> Result<Vec<Sym<2>>, String> {
    if m2.dim != 2 {
        return Err("electro_overlay_metrics_2d requires a 2D mesh".into());
    }
    let nt = m2.n_cells();
    let nv = m2.n_points();
    if eta2_2d.len() != nt {
        return Err(format!("{} 2D indicators for {nt} triangles", eta2_2d.len()));
    }
    if !(opt.h_max > 0.0 && opt.h_max.is_finite() && (0.0..=1.0).contains(&opt.h_min_fraction) && opt.target_cells >= 1.0)
    {
        return Err("h_max and target_cells must be positive and h_min_fraction in [0, 1]".into());
    }
    if let Some(fl) = floor {
        if fl.len() != nv {
            return Err(format!("{} floor metrics for {nv} vertices", fl.len()));
        }
    }
    if nv == 0 || nt == 0 {
        return Ok(Vec::new());
    }

    let q_smooth = (opt.order + 1) as f64;
    let q_sing = opt.q_sing.max(0.25);
    let h_min = (opt.h_max * opt.h_min_fraction).max(1e-30 * opt.h_max).min(opt.h_max);
    let max_aspect = opt.max_aspect.max(1.0);
    let rim_aspect_cap = opt.rim_anisotropy.clamp(1.0, max_aspect);

    // Per-triangle size `s_t`, area, local exponent `q_t`, and rim normal `n_t` when inside the rim band.
    struct TriInfo {
        s: f64,
        area: f64,
        q: f64,
        rim_normal: Option<[f64; 2]>,
    }
    let mut tinfo = Vec::with_capacity(nt);
    let mut v_local_s = vec![f64::INFINITY; nv];
    for t in 0..nt {
        let v = m2.cell(t);
        let (a, b, c) = (pt2(m2, v[0]), pt2(m2, v[1]), pt2(m2, v[2]));
        let area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs();
        let s = (4.0 * area / 3.0_f64.sqrt()).sqrt().max(1e-30);
        for &vi in v {
            v_local_s[vi as usize] = v_local_s[vi as usize].min(s);
        }
        let cen = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0];
        let (q_loc, rim_n) = if let Some((d_rim, n_hat)) = closest_rim_seg(cen, opt.rim_segments) {
            if d_rim <= 2.0 * s {
                (q_sing, Some(n_hat))
            } else if d_rim < 4.0 * s {
                let w = (d_rim - 2.0 * s) / (2.0 * s);
                (q_sing * (1.0 - w) + q_smooth * w, None)
            } else {
                (q_smooth, None)
            }
        } else {
            (q_smooth, None)
        };
        tinfo.push(TriInfo { s, area, q: q_loc, rim_normal: rim_n });
    }

    let pts2: Vec<[f64; 2]> = (0..nv as u32).map(|vi| pt2(m2, vi)).collect();
    let edges2: Vec<[u32; 2]> = m2.edges().iter().map(|&(a, b)| [a, b]).collect();
    let eig_lo = 1.0 / (opt.h_max * opt.h_max);
    let eig_hi = 1.0 / (h_min * h_min);
    let id2 = identity::<2>();
    let max_coarsen_ratio = if opt.target_cells >= nt as f64 { 1.25 } else { 4.0 };

    let cell_sizes_at = |tau: f64, t: usize| -> (f64, f64) {
        let ti = &tinfo[t];
        let eta = eta2_2d[t].max(0.0).sqrt();
        let ratio = if eta > 0.0 { (tau / eta).powf(1.0 / ti.q) } else { max_coarsen_ratio };
        let h_n = (ti.s * ratio.clamp(0.2, max_coarsen_ratio)).clamp(h_min, opt.h_max);
        let h_t = if ti.rim_normal.is_some() && rim_aspect_cap > 1.0 {
            (h_n * rim_aspect_cap).min(ti.s.max(h_n)).clamp(h_min, opt.h_max)
        } else {
            h_n
        };
        (h_n, h_t)
    };

    // Build ungraded vertex metric at threshold `tau` (including floor and pre-gradation isotropic corner balls).
    let build_ungraded_at = |tau: f64| -> Vec<Sym<2>> {
        let mut v_met = vec![isotropic::<2>(opt.h_max); nv];
        for t in 0..nt {
            let (hn, ht) = cell_sizes_at(tau, t);
            let mt = match tinfo[t].rim_normal {
                Some(n) => {
                    let dirs = [[n[0], -n[1]], [n[1], n[0]]];
                    from_sizes(&[hn, ht], &dirs)
                }
                None => from_sizes(&[hn, ht], &id2),
            };
            for &vi in m2.cell(t) {
                v_met[vi as usize] = intersect(&v_met[vi as usize], &mt);
            }
        }
        if let Some(fl) = floor {
            for (m, f) in v_met.iter_mut().zip(fl) {
                *m = intersect(m, f);
            }
        }
        // Isotropic corner-ball intersection BEFORE `BackgroundMetric2::new_graded` so orthogonal
        // anisotropic rim layers never collide at 90/270 deg corners.
        if !opt.corners.is_empty() {
            let ball_fac = opt.corner_ball_factor.max(1.0);
            for &cp in opt.corners {
                let mut best_v = (0usize, f64::INFINITY);
                for (vi, &p) in pts2.iter().enumerate() {
                    let d = (p[0] - cp[0]).hypot(p[1] - cp[1]);
                    if d < best_v.1 {
                        best_v = (vi, d);
                    }
                }
                let s_loc = v_local_s[best_v.0];
                if !s_loc.is_finite() {
                    continue;
                }
                let r_ball = ball_fac * s_loc;
                let mut h_corner = sizes(&v_met[best_v.0]).0[0];
                for (vi, &p) in pts2.iter().enumerate() {
                    if (p[0] - cp[0]).hypot(p[1] - cp[1]) <= s_loc {
                        h_corner = h_corner.min(sizes(&v_met[vi]).0[0]);
                    }
                }
                h_corner = h_corner.clamp(h_min, opt.h_max);
                for (vi, &p) in pts2.iter().enumerate() {
                    let d = (p[0] - cp[0]).hypot(p[1] - cp[1]);
                    if d <= r_ball {
                        let h_iso = (h_corner + opt.grade.max(0.1) * d).clamp(h_min, opt.h_max);
                        if d <= 0.5 * r_ball {
                            let cur_min = sizes(&v_met[vi]).0[0].min(h_iso).clamp(h_min, opt.h_max);
                            v_met[vi] = isotropic::<2>(cur_min);
                        } else {
                            v_met[vi] = intersect(&v_met[vi], &isotropic::<2>(h_iso));
                        }
                    }
                }
            }
        }
        for m in &mut v_met {
            let clamped = map_eig(m, |l| l.clamp(eig_lo, eig_hi));
            let (hs, dirs) = sizes(&clamped);
            let h0 = hs[0].clamp(h_min, opt.h_max);
            let h1 = hs[1].clamp(h_min, (h0 * max_aspect).min(opt.h_max));
            *m = from_sizes(&[h0, h1], &dirs);
        }
        v_met
    };

    let qr = quad_rule(2, 4);
    let w_sum: f64 = qr.w.iter().sum::<f64>().max(1e-300);
    let graded_elements_at = |tau: f64| -> f64 {
        let mut m = build_ungraded_at(tau);
        if opt.grade > 0.0 && opt.grade.is_finite() {
            grade(&mut m, &pts2, &edges2, opt.grade, 60);
        }
        // Exact log-Euclidean complexity integral matching `BackgroundMetric2`:
        // det(exp(sum_i lambda_i log M_i))^{1/2} == exp(0.5 * sum_i lambda_i ln det M_i).
        let half_ln_det: Vec<f64> = m
            .iter()
            .map(|a| 0.5 * (a[0][0] * a[1][1] - a[0][1] * a[1][0]).max(1e-300).ln())
            .collect();
        let mut cx = 0.0_f64;
        for (t, ti) in tinfo.iter().enumerate() {
            let v = m2.cell(t);
            let (l0, l1, l2) = (half_ln_det[v[0] as usize], half_ln_det[v[1] as usize], half_ln_det[v[2] as usize]);
            let mut acc = 0.0_f64;
            for (lam, &wq) in qr.lam.iter().zip(&qr.w) {
                acc += wq * (lam[0] * l0 + lam[1] * l1 + lam[2] * l2).exp();
            }
            cx += (ti.area / w_sum) * acc;
        }
        elements_for_complexity(cx, 2)
    };

    let max_eta = eta2_2d.iter().copied().fold(0.0_f64, f64::max).sqrt().max(1e-300);
    let (mut lo, mut hi) = (max_eta * 1e-12, max_eta * 1e12);
    for _ in 0..28 {
        let mid = (lo * hi).sqrt();
        if graded_elements_at(mid) > opt.target_cells {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let tau = (lo * hi).sqrt();
    Ok(build_ungraded_at(tau))
}

// ------------------------------------------------------------------ 3D sizing & metrics for `full_3d`

/// Options for 3D electrostatic equidistribution [`electro_vertex_sizes_3d`] / [`electro_metrics_3d`].
#[derive(Clone, Debug)]
pub struct ElectroMetric3dOptions<'a> {
    /// Lagrange order (`1` or `2`); smooth 3D exponent is `order as f64 + 1.5`.
    pub order: usize,
    /// Singular exponent near 3D sheet rim vertices (default `1.0`).
    pub q_sing: f64,
    /// Largest element size.
    pub h_max: f64,
    /// Smallest element size as a fraction of `h_max`.
    pub h_min_fraction: f64,
    /// Target number of 3D tetrahedra.
    pub target_cells: f64,
    /// Optional mask of length `mesh.n_points()` marking vertices on open `sheet:edge` rims.
    pub rim_vertices: Option<&'a [bool]>,
}

impl<'a> ElectroMetric3dOptions<'a> {
    pub fn new(order: usize, h_max: f64, h_min_fraction: f64, target_cells: f64) -> Self {
        ElectroMetric3dOptions { order, q_sing: 1.0, h_max, h_min_fraction, target_cells, rim_vertices: None }
    }
}

/// Computes target 3D per-vertex element sizes `vh[v]` in `[h_max * h_min_fraction, h_max]` from
/// 3D cell indicators `eta2` for `full_3d` unstructured adaptation.
pub fn electro_vertex_sizes_3d(
    mesh: &SimplexMesh,
    eta2: &[f64],
    opt: &ElectroMetric3dOptions,
) -> Result<Vec<f64>, String> {
    if mesh.dim != 3 {
        return Err("electro_vertex_sizes_3d requires a 3D tetrahedral mesh".into());
    }
    let nc = mesh.n_cells();
    let nv = mesh.n_points();
    if eta2.len() != nc {
        return Err(format!("{} indicators for {nc} cells", eta2.len()));
    }
    if !(opt.h_max > 0.0 && opt.h_max.is_finite() && (0.0..=1.0).contains(&opt.h_min_fraction) && opt.target_cells >= 1.0)
    {
        return Err("h_max and target_cells must be positive and h_min_fraction in [0, 1]".into());
    }
    if let Some(rv) = opt.rim_vertices {
        if rv.len() != nv {
            return Err(format!("rim_vertices length {} != n_points {nv}", rv.len()));
        }
    }
    if nv == 0 || nc == 0 {
        return Ok(vec![opt.h_max; nv]);
    }

    let q_smooth = opt.order as f64 + 1.5;
    let q_sing = opt.q_sing.max(0.25);
    let h_min = (opt.h_max * opt.h_min_fraction).max(1e-30 * opt.h_max).min(opt.h_max);
    let geo = geometry(mesh);
    let size: Vec<f64> = (0..nc).map(|c| (6.0 * 2.0_f64.sqrt() * geo.vol[c]).cbrt().max(1e-30)).collect();
    let cell_q: Vec<f64> = (0..nc)
        .map(|c| {
            let on_rim = opt.rim_vertices.is_some_and(|rv| mesh.cell(c).iter().any(|&v| rv[v as usize]));
            if on_rim { q_sing } else { q_smooth }
        })
        .collect();

    let new_size = |tau: f64, c: usize| -> f64 {
        let eta = eta2[c].max(0.0).sqrt();
        let ratio = if eta > 0.0 { (tau / eta).powf(1.0 / cell_q[c]) } else { 4.0 };
        (size[c] * ratio.clamp(0.25, 4.0)).clamp(h_min, opt.h_max)
    };
    let predicted = |tau: f64| -> f64 { (0..nc).map(|c| (size[c] / new_size(tau, c)).powi(3)).sum() };

    let max_eta = eta2.iter().copied().fold(0.0_f64, f64::max).sqrt().max(1e-300);
    let (mut lo, mut hi) = (max_eta * 1e-14, max_eta * 1e14);
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        if predicted(mid) > opt.target_cells {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let tau = (lo * hi).sqrt();
    let mut vh = vec![opt.h_max; nv];
    for c in 0..nc {
        let s = new_size(tau, c);
        for &v in mesh.cell(c) {
            vh[v as usize] = vh[v as usize].min(s);
        }
    }
    Ok(vh)
}

/// Computes target 3D per-vertex metrics `Vec<Sym<3>>` from 3D cell indicators `eta2`,
/// intersecting with an optional `floor` metric field.
pub fn electro_metrics_3d(
    mesh: &SimplexMesh,
    eta2: &[f64],
    opt: &ElectroMetric3dOptions,
    floor: Option<&[Sym<3>]>,
) -> Result<Vec<Sym<3>>, String> {
    if let Some(fl) = floor {
        if fl.len() != mesh.n_points() {
            return Err(format!("{} floor metrics for {} vertices", fl.len(), mesh.n_points()));
        }
    }
    let vh = electro_vertex_sizes_3d(mesh, eta2, opt)?;
    let mut out: Vec<Sym<3>> = vh.into_iter().map(|h| isotropic::<3>(h.max(1e-300))).collect();
    if let Some(fl) = floor {
        for (a, b) in out.iter_mut().zip(fl) {
            *a = intersect(a, b);
        }
    }
    Ok(out)
}
