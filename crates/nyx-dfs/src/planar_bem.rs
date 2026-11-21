// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 2D layered-dielectric boundary-element method (`formulation = "planar_bem"`) for zero-thickness planar conductor sheets.
//!
//! Extracts planar conductor sheet triangles on `z = z_sheet` (from 3D `sheet3d` interior conductor facets or a 2D mesh),
//! evaluates the layered-dielectric reflection image series Green's function `G(rho)` with exact analytical polar triangle
//! integration on self and near panel pairs, assembles the symmetric positive definite elastance matrix `P`, factors via
//! dense eigendecomposition / Cholesky, computes `C = K^T P^-1 K`, and optionally performs 1->4 panel subdivision with
//! log-space SPD Richardson extrapolation (`richardson = true`).

use crate::common::{components, minmax, norm_of};
use crate::electro::{electrostatics_fields, EPS0};
use crate::problem::{BcKind, Problem};
use femcore::fem::{boundary_dofs, check_boundary_labels};
use femcore::linalg::driver::dense_eigh;
use femcore::material;
use nyx_physics::{cell_fields_of, CapacitanceMatrixData, Discretisation, Scalar, Solution, SolveResult, TlsMeta};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::time::Instant;

const E_CHARGE: f64 = 1.602_176_634e-19;
const PLANCK_H: f64 = 6.626_070_15e-34;
const MAX_IMAGES: usize = 64;
const IMAGE_TOL: f64 = 1e-10;

type QPt = (i64, i64);
type QEdge = (QPt, QPt);
type Matrix = Vec<Vec<f64>>;

/// Three-layer planar dielectric stack around a conductor sheet plane `z = z_sheet`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredStack {
    /// Relative permittivity of the upper half-space `z > z_sheet` (`eps_1`).
    pub eps_top: f64,
    /// Relative permittivity of the substrate layer `z_sheet - t_sub < z < z_sheet` (`eps_2`).
    pub eps_sub: f64,
    /// Relative permittivity of the lower half-space `z < z_sheet - t_sub` (`eps_3`).
    pub eps_bot: f64,
    /// Substrate thickness `t_sub` (m); `f64::INFINITY` for a semi-infinite substrate.
    pub t_sub: f64,
    /// Whether the bottom of the substrate (`z = z_sheet - t_sub`) is a conducting ground plane (`k_23 = -1`).
    pub grounded_backplane: bool,
}

impl LayeredStack {
    /// Effective half-space relative permittivity `eps_eff = (eps_1 + eps_2) / 2`.
    pub fn eps_eff(&self) -> f64 {
        0.5 * (self.eps_top + self.eps_sub)
    }

    /// Reflection coefficient `k_21 = (eps_2 - eps_1) / (eps_2 + eps_1)`.
    pub fn k21(&self) -> f64 {
        (self.eps_sub - self.eps_top) / (self.eps_sub + self.eps_top)
    }

    /// Reflection coefficient `k_23 = -1` for a grounded backplane, or `(eps_2 - eps_3) / (eps_2 + eps_3)`.
    pub fn k23(&self) -> f64 {
        if self.grounded_backplane {
            -1.0
        } else {
            (self.eps_sub - self.eps_bot) / (self.eps_sub + self.eps_bot)
        }
    }

    /// Precomputed image weights and vertical image depths `(w_n, d_n = 2*n*t_sub)` for `n = 1..=N`,
    /// where `G(rho) = 1 / (4*pi*eps0*eps_eff) * [ 1/rho + sum_n w_n / sqrt(rho^2 + d_n^2) ]`.
    pub fn images(&self) -> Vec<(f64, f64)> {
        if !self.t_sub.is_finite() || self.t_sub <= 0.0 {
            return Vec::new();
        }
        let k21 = self.k21();
        let k23 = self.k23();
        if k23.abs() < 1e-15 {
            return Vec::new();
        }
        let scale = 2.0 * self.eps_sub / (self.eps_top + self.eps_sub);
        let mut out = Vec::new();
        let mut coeff = k23;
        let mut last_n = 0usize;
        for n in 1..=MAX_IMAGES {
            if coeff.abs() < IMAGE_TOL {
                break;
            }
            out.push((scale * coeff, 2.0 * (n as f64) * self.t_sub));
            last_n = n;
            coeff *= k21 * k23;
        }
        // If truncated at MAX_IMAGES with a non-negligible geometric tail (e.g. near t_sub -> 0),
        // append the exact geometric series remainder at depth 2*(MAX_IMAGES + 1)*t_sub so the
        // t_sub -> 0 limit converges to machine precision.
        let r = k21 * k23;
        if last_n == MAX_IMAGES && coeff.abs() >= IMAGE_TOL && (1.0 - r).abs() > 1e-12 {
            out.push((scale * coeff / (1.0 - r), 2.0 * ((MAX_IMAGES + 1) as f64) * self.t_sub));
        }
        out
    }
}

/// Evaluates the exact layered-dielectric reflection image Green's function `G(rho)` on `z = z' = z_sheet`
/// for horizontal distance `rho > 0`.
pub fn layered_greens_function(rho: f64, stack: &LayeredStack) -> f64 {
    assert!(rho > 0.0, "rho must be positive");
    let c0 = 1.0 / (4.0 * PI * EPS0 * stack.eps_eff());
    let mut s = 1.0 / rho;
    for (w_n, d_n) in stack.images() {
        s += w_n / (rho * rho + d_n * d_n).sqrt();
    }
    c0 * s
}

/// Exact analytical surface integral `int_T dA' / sqrt(|p - r'|^2 + d^2)` over a counter-clockwise
/// 2D triangle `tri = [v0, v1, v2]` for any observation point `p` (interior, boundary, or exterior)
/// and vertical image offset `d >= 0` (Wilton et al. 1984 polar line-integral formula).
fn triangle_layer_potential_at(p: [f64; 2], tri: [[f64; 2]; 3], d: f64) -> f64 {
    let d_abs = d.abs();
    let d2 = d_abs * d_abs;
    let mut acc = 0.0;
    for i in 0..3 {
        let a = [tri[i][0] - p[0], tri[i][1] - p[1]];
        let b = [tri[(i + 1) % 3][0] - p[0], tri[(i + 1) % 3][1] - p[1]];
        let e = [b[0] - a[0], b[1] - a[1]];
        let len = (e[0] * e[0] + e[1] * e[1]).sqrt();
        if len <= 1e-30 {
            continue;
        }
        let t = [e[0] / len, e[1] / len];
        // Signed perpendicular distance to the CCW edge: positive when `p` is on the interior side
        let h_signed = (a[0] * e[1] - a[1] * e[0]) / len;
        let h_abs = h_signed.abs();
        if h_abs <= 1e-30 {
            continue;
        }
        let s0 = a[0] * t[0] + a[1] * t[1];
        let s1 = b[0] * t[0] + b[1] * t[1];
        let hd = (h_abs * h_abs + d2).sqrt();
        let term_asinh = h_signed * ((s1 / hd).asinh() - (s0 / hd).asinh());
        let term_atan = if d_abs > 1e-30 {
            let r0 = (hd * hd + s0 * s0).sqrt();
            let r1 = (hd * hd + s1 * s1).sqrt();
            let a1 = (d_abs * s1 / (h_abs * r1)).atan() - (s1 / h_abs).atan();
            let a0 = (d_abs * s0 / (h_abs * r0)).atan() - (s0 / h_abs).atan();
            h_signed.signum() * d_abs * (a1 - a0)
        } else {
            0.0
        };
        acc += term_asinh + term_atan;
    }
    acc
}

#[derive(Clone, Debug)]
struct Panel {
    v: [[f64; 2]; 3],
    area: f64,
    center: [f64; 2],
    radius: f64,
    /// `Some(term_idx)` for an excited terminal panel, `None` for a passive grounded sheet panel.
    terminal: Option<usize>,
    /// True if panel has at least one vertex on the outer perimeter rim of its conductor island.
    is_rim: bool,
    /// Effective edge-singularity weight multiplier (> 1 on rim panels to account for `(d)^(-1/2)` charge crowding).
    q_pts_4: [[f64; 2]; 4],
    q_pts_6: [[f64; 2]; 6],
}

// 4-point symmetric degree-3 quadrature on a triangle (barycentric coordinates and weights summing to 1)
const QUAD4_BARY: [[f64; 3]; 4] =
    [[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], [0.6, 0.2, 0.2], [0.2, 0.6, 0.2], [0.2, 0.2, 0.6]];
const QUAD4_W: [f64; 4] = [-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0];

// 6-point Dunavant degree-4 positive-weight quadrature on a triangle
const QUAD6_BARY: [[f64; 3]; 6] = [
    [0.108103018168070, 0.445948490915965, 0.445948490915965],
    [0.445948490915965, 0.108103018168070, 0.445948490915965],
    [0.445948490915965, 0.445948490915965, 0.108103018168070],
    [0.816847572980459, 0.091576213509771, 0.091576213509771],
    [0.091576213509771, 0.816847572980459, 0.091576213509771],
    [0.091576213509771, 0.091576213509771, 0.816847572980459],
];
const QUAD6_W: [f64; 6] =
    [0.223381589678011, 0.223381589678011, 0.223381589678011, 0.109951743655322, 0.109951743655322, 0.109951743655322];

impl Panel {
    fn new(mut v: [[f64; 2]; 3], terminal: Option<usize>, is_rim: bool) -> Option<Self> {
        let det = (v[1][0] - v[0][0]) * (v[2][1] - v[0][1]) - (v[1][1] - v[0][1]) * (v[2][0] - v[0][0]);
        if det.abs() <= 1e-30 {
            return None;
        }
        if det < 0.0 {
            v.swap(1, 2);
        }
        let area = 0.5 * det.abs();
        let center = [(v[0][0] + v[1][0] + v[2][0]) / 3.0, (v[0][1] + v[1][1] + v[2][1]) / 3.0];
        let radius = (0..3).map(|i| (v[i][0] - center[0]).hypot(v[i][1] - center[1])).fold(0.0_f64, f64::max);
        let interp = |l: &[f64; 3]| -> [f64; 2] {
            [l[0] * v[0][0] + l[1] * v[1][0] + l[2] * v[2][0], l[0] * v[0][1] + l[1] * v[1][1] + l[2] * v[2][1]]
        };
        let q_pts_4 = [interp(&QUAD4_BARY[0]), interp(&QUAD4_BARY[1]), interp(&QUAD4_BARY[2]), interp(&QUAD4_BARY[3])];
        let q_pts_6 = [
            interp(&QUAD6_BARY[0]),
            interp(&QUAD6_BARY[1]),
            interp(&QUAD6_BARY[2]),
            interp(&QUAD6_BARY[3]),
            interp(&QUAD6_BARY[4]),
            interp(&QUAD6_BARY[5]),
        ];
        Some(Panel { v, area, center, radius, terminal, is_rim, q_pts_4, q_pts_6 })
    }
}

/// Subdivides each panel 1->4 by connecting its edge midpoints. If a panel edge lies on the open
/// boundary of a nearly circular conductor disk centered at `(cx, cy)` of radius `R`, projects the
/// boundary midpoint onto the circle so curved disk boundaries converge without polygon-chord locking.
fn subdivide_panels_1_to_4(panels: &[Panel]) -> Vec<Panel> {
    // Detect if the outer perimeter vertices of each terminal form a circle: if all perimeter vertices
    // of a terminal have equal distance R from their centroid within 1e-3 * R, record (center, R) for arc midpoint projection.
    let q_scale = panels.iter().map(|p| p.radius).fold(0.0_f64, f64::max).max(1e-30) * 1e-7;
    let qpt = |p: [f64; 2]| -> QPt { ((p[0] / q_scale).round() as i64, (p[1] / q_scale).round() as i64) };
    let mut edge_counts: BTreeMap<QEdge, usize> = BTreeMap::new();
    for p in panels {
        for i in 0..3 {
            let a = qpt(p.v[i]);
            let b = qpt(p.v[(i + 1) % 3]);
            let key = if a <= b { (a, b) } else { (b, a) };
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }
    let mut rim_pts_by_term: BTreeMap<Option<usize>, Vec<[f64; 2]>> = BTreeMap::new();
    for p in panels {
        for i in 0..3 {
            let a = qpt(p.v[i]);
            let b = qpt(p.v[(i + 1) % 3]);
            let key = if a <= b { (a, b) } else { (b, a) };
            if edge_counts.get(&key).copied().unwrap_or(0) == 1 {
                let entry = rim_pts_by_term.entry(p.terminal).or_default();
                entry.push(p.v[i]);
                entry.push(p.v[(i + 1) % 3]);
            }
        }
    }
    let mut circle_by_term: BTreeMap<Option<usize>, ([f64; 2], f64)> = BTreeMap::new();
    for (term, pts) in &rim_pts_by_term {
        if pts.len() >= 16 {
            let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
            let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
            for pt in pts {
                xmin = xmin.min(pt[0]);
                xmax = xmax.max(pt[0]);
                ymin = ymin.min(pt[1]);
                ymax = ymax.max(pt[1]);
            }
            let c = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax)];
            let r_mean = pts.iter().map(|pt| (pt[0] - c[0]).hypot(pt[1] - c[1])).sum::<f64>() / pts.len() as f64;
            if r_mean > 0.0 && pts.iter().all(|pt| ((pt[0] - c[0]).hypot(pt[1] - c[1]) - r_mean).abs() < 2e-3 * r_mean)
            {
                circle_by_term.insert(*term, (c, r_mean));
            }
        }
    }

    let mut out = Vec::with_capacity(panels.len() * 4);
    for p in panels {
        let [v0, v1, v2] = p.v;
        let mid_edge = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] {
            let mut m = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
            let qa = qpt(a);
            let qb = qpt(b);
            let key = if qa <= qb { (qa, qb) } else { (qb, qa) };
            if edge_counts.get(&key).copied().unwrap_or(0) == 1 {
                if let Some(&(c, r)) = circle_by_term.get(&p.terminal) {
                    let dm = (m[0] - c[0]).hypot(m[1] - c[1]);
                    if dm > 1e-30 {
                        m = [c[0] + (m[0] - c[0]) * (r / dm), c[1] + (m[1] - c[1]) * (r / dm)];
                    }
                }
            }
            m
        };
        let m01 = mid_edge(v0, v1);
        let m12 = mid_edge(v1, v2);
        let m20 = mid_edge(v2, v0);
        for sub in [[v0, m01, m20], [m01, v1, m12], [m20, m12, v2], [m01, m12, m20]] {
            if let Some(sp) = Panel::new(sub, p.terminal, p.is_rim) {
                out.push(sp);
            }
        }
    }
    out
}

/// Evaluates `P_{ij} = 1 / (A_i A_j) int_{T_i} int_{T_j} G(|r - r'|) dA' dA` for panel pair `(pi, pj)`.
fn panel_pair_elastance(pi: &Panel, pj: &Panel, c0: f64, images: &[(f64, f64)], same_idx: bool) -> f64 {
    let dist = (pi.center[0] - pj.center[0]).hypot(pi.center[1] - pj.center[1]);
    let h_pair = pi.radius + pj.radius;
    let near = same_idx || dist < 4.0 * h_pair;

    let mut val = 0.0;
    if near {
        // Exact analytical polar integration over inner triangle combined with 6-point Dunavant
        // quadrature over outer triangle, symmetrized between (pi, pj) and (pj, pi).
        let eval_dir = |outer: &Panel, inner: &Panel| -> f64 {
            let mut s0 = 0.0;
            for q in 0..6 {
                let pt = outer.q_pts_6[q];
                s0 += QUAD6_W[q] * triangle_layer_potential_at(pt, inner.v, 0.0);
            }
            s0 / inner.area
        };
        val += if same_idx { eval_dir(pi, pj) } else { 0.5 * (eval_dir(pi, pj) + eval_dir(pj, pi)) };
    } else {
        // 4-point symmetric Gauss quadrature on both distant triangles
        let mut s0 = 0.0;
        for qa in 0..4 {
            let pa = pi.q_pts_4[qa];
            let wa = QUAD4_W[qa];
            for qb in 0..4 {
                let pb = pj.q_pts_4[qb];
                let rho = (pa[0] - pb[0]).hypot(pa[1] - pb[1]).max(1e-30);
                s0 += wa * QUAD4_W[qb] / rho;
            }
        }
        val += s0;
    }

    for &(w_n, d_n) in images {
        let d_eff = (dist * dist + d_n * d_n).sqrt();
        if d_eff < 3.5 * h_pair {
            let eval_img = |outer: &Panel, inner: &Panel| -> f64 {
                let mut sn = 0.0;
                for q in 0..6 {
                    let pt = outer.q_pts_6[q];
                    sn += QUAD6_W[q] * triangle_layer_potential_at(pt, inner.v, d_n);
                }
                sn / inner.area
            };
            let img = if same_idx { eval_img(pi, pj) } else { 0.5 * (eval_img(pi, pj) + eval_img(pj, pi)) };
            val += w_n * img;
        } else {
            let mut sn = 0.0;
            for qa in 0..4 {
                let pa = pi.q_pts_4[qa];
                let wa = QUAD4_W[qa];
                for qb in 0..4 {
                    let pb = pj.q_pts_4[qb];
                    let dx = pa[0] - pb[0];
                    let dy = pa[1] - pb[1];
                    let r = (dx * dx + dy * dy + d_n * d_n).sqrt().max(1e-30);
                    sn += wa * QUAD4_W[qb] / r;
                }
            }
            val += w_n * sn;
        }
    }
    c0 * val
}

/// Assembles the symmetric `n x n` BEM elastance matrix `P`, factors it via `dense_eigh`, and
/// returns `(maxwell_matrix_f, panel_charges_per_terminal)`.
fn assemble_and_solve_bem(panels: &[Panel], n_term: usize, stack: &LayeredStack) -> Result<(Matrix, Matrix), String> {
    let n = panels.len();
    if n == 0 {
        return Err("planar_bem found zero conductor sheet triangles".into());
    }
    let c0 = 1.0 / (4.0 * PI * EPS0 * stack.eps_eff());
    let images = stack.images();

    // Assemble upper triangle of symmetric P deterministically in row chunks on the shared pool
    let par = n >= 48 && femcore::par::threads() > 1;
    let rows: Vec<Vec<f64>> = femcore::linalg::par::chunk_map(n, par, |i| {
        let mut row_upper = Vec::with_capacity(n - i);
        for j in i..n {
            row_upper.push(panel_pair_elastance(&panels[i], &panels[j], c0, &images, i == j));
        }
        row_upper
    });

    let mut p_flat = vec![0.0; n * n];
    for i in 0..n {
        for (off, &v) in rows[i].iter().enumerate() {
            let j = i + off;
            p_flat[i * n + j] = v;
            p_flat[j * n + i] = v;
        }
    }

    let (eigs, vecs) = dense_eigh(&p_flat, n)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    if lam_max <= 0.0 {
        return Err("planar_bem elastance matrix P is not positive definite".into());
    }

    // Compute gamma[r][t] = u_r^T K_{:, t} and C = K^T P^-1 K
    let mut c_mat = vec![vec![0.0; n_term]; n_term];
    let mut q_cols = vec![vec![0.0; n]; n_term];
    for r in 0..n {
        if eigs[r] > 1e-13 * lam_max {
            let inv_l = 1.0 / eigs[r];
            let ur = &vecs[r];
            let mut gamma = vec![0.0; n_term];
            for (i, pan) in panels.iter().enumerate() {
                if let Some(t) = pan.terminal {
                    gamma[t] += ur[i];
                }
            }
            for a in 0..n_term {
                let sc_a = inv_l * gamma[a];
                if sc_a != 0.0 {
                    for i in 0..n {
                        q_cols[a][i] += sc_a * ur[i];
                    }
                }
                for b in a..n_term {
                    let v = sc_a * gamma[b];
                    c_mat[a][b] += v;
                    if a != b {
                        c_mat[b][a] += v;
                    }
                }
            }
        }
    }
    Ok((c_mat, q_cols))
}

/// Applies log-space SPD Richardson extrapolation (`r = 2`, `w_1 = 2.0`) to coarse (`c0`) and
/// 1->4 subdivided fine (`c1`) Maxwell capacitance matrices.
pub fn extrapolate_spd_log_space(c0: &[Vec<f64>], c1: &[Vec<f64>], w1: f64) -> Result<Vec<Vec<f64>>, String> {
    let m = c0.len();
    if m == 0 {
        return Ok(Vec::new());
    }
    if m == 1 {
        let v0 = c0[0][0].max(1e-30);
        let v1 = c1[0][0].max(1e-30);
        let ext = (w1 * v1.ln() + (1.0 - w1) * v0.ln()).exp();
        return Ok(vec![vec![ext]]);
    }

    // Try positive physical mutual/ground component log-space extrapolation first (preserves sign of all off-diagonals
    // and strict diagonal dominance), falling back to Log-Euclidean matrix log/exp via dense_eigh if any ground sum <= 0.
    let mut g0 = vec![0.0; m];
    let mut g1 = vec![0.0; m];
    let mut all_pos = true;
    for i in 0..m {
        g0[i] = c0[i].iter().sum();
        g1[i] = c1[i].iter().sum();
        if g0[i] <= 1e-30 || g1[i] <= 1e-30 {
            all_pos = false;
        }
        for j in (i + 1)..m {
            if -c0[i][j] <= 1e-30 || -c1[i][j] <= 1e-30 {
                all_pos = false;
            }
        }
    }
    if all_pos {
        let mut ext = vec![vec![0.0; m]; m];
        for i in 0..m {
            let gi = (w1 * g1[i].ln() + (1.0 - w1) * g0[i].ln()).exp();
            ext[i][i] += gi;
            for j in (i + 1)..m {
                let m0 = -c0[i][j];
                let m1 = -c1[i][j];
                let mij = (w1 * m1.ln() + (1.0 - w1) * m0.ln()).exp();
                ext[i][j] = -mij;
                ext[j][i] = -mij;
                ext[i][i] += mij;
                ext[j][j] += mij;
            }
        }
        return Ok(ext);
    }

    // Matrix Log-Euclidean SPD extrapolation: exp(w1 * log(C1) + (1 - w1) * log(C0))
    let mat_log = |c: &[Vec<f64>]| -> Result<Vec<f64>, String> {
        let mut flat = Vec::with_capacity(m * m);
        for r in c {
            flat.extend_from_slice(r);
        }
        let (eigs, vecs) = dense_eigh(&flat, m)?;
        let lmax = eigs.iter().copied().fold(0.0_f64, f64::max).max(1e-30);
        let mut out = vec![0.0; m * m];
        for k in 0..m {
            let lk = eigs[k].max(1e-14 * lmax).ln();
            for i in 0..m {
                for j in 0..m {
                    out[i * m + j] += lk * vecs[k][i] * vecs[k][j];
                }
            }
        }
        Ok(out)
    };
    let l0 = mat_log(c0)?;
    let l1 = mat_log(c1)?;
    let l_ext: Vec<f64> = l0.iter().zip(&l1).map(|(&a, &b)| w1 * b + (1.0 - w1) * a).collect();
    let (eigs, vecs) = dense_eigh(&l_ext, m)?;
    let mut out = vec![vec![0.0; m]; m];
    for k in 0..m {
        let ek = eigs[k].exp();
        for i in 0..m {
            for j in 0..m {
                out[i][j] += ek * vecs[k][i] * vecs[k][j];
            }
        }
    }
    for i in 0..m {
        for j in (i + 1)..m {
            let s = 0.5 * (out[i][j] + out[j][i]);
            out[i][j] = s;
            out[j][i] = s;
        }
    }
    Ok(out)
}

fn domain_scalar_eps(mesh: &Discretisation, spec: &Problem) -> BTreeMap<i32, f64> {
    let mut map = BTreeMap::new();
    for d in mesh.domain_labels() {
        let eps = if let Some(t) = spec.materials.aniso.get(&d).and_then(|a| a.permittivity.as_ref()) {
            let m = material::tensor_global(t, spec.materials.aniso[&d].euler_deg);
            // For layered isotropic/uniaxial substrates, the geometric mean sqrt(eps_xy * eps_z) or arithmetic trace
            (m[0][0] + m[1][1] + m[2][2]) / 3.0
        } else {
            spec.materials.get(d, "eps_r").or_else(|| spec.materials.default.get("eps_r").copied()).unwrap_or(1.0)
        };
        map.insert(d, eps);
    }
    map
}

/// Checks whether `mesh` and `spec` represent a planar sheet capacitance problem suitable for
/// automatic selection (`formulation = "auto"`).
pub fn can_use_planar_bem(mesh: &Discretisation, spec: &Problem) -> bool {
    if mesh.dim != 3 {
        return false;
    }
    let has_terminal = spec.bcs.iter().any(|b| matches!(b.kind, BcKind::Terminal(_)));
    if !has_terminal {
        return false;
    }
    extract_3d_panels_and_stack(mesh, spec).is_ok()
}

#[allow(clippy::type_complexity)]
fn extract_3d_panels_and_stack(
    mesh: &Discretisation,
    spec: &Problem,
) -> Result<(Vec<(String, Vec<i32>)>, Vec<Panel>, LayeredStack, f64), String> {
    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();
    let mut ground_labels: BTreeSet<i32> = BTreeSet::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::Terminal(name) => {
                let name = name.trim();
                if name.is_empty() {
                    return Err("terminal condition requires a non-empty 'name'".into());
                }
                if let Some((_, labs)) = terminals.iter_mut().find(|(n, _)| n == name) {
                    for &l in &bc.labels {
                        if !labs.contains(&l) {
                            labs.push(l);
                        }
                    }
                } else {
                    terminals.push((name.to_string(), bc.labels.clone()));
                }
            }
            BcKind::Ground => {
                ground_labels.extend(bc.labels.iter().copied());
            }
            BcKind::Potential(v) if v.abs() < 1e-15 => {
                ground_labels.extend(bc.labels.iter().copied());
            }
            _ => {}
        }
    }
    if terminals.is_empty() {
        return Err("planar_bem requires at least one terminal boundary condition".into());
    }

    let mut label_to_term: BTreeMap<i32, usize> = BTreeMap::new();
    for (idx, (_, labs)) in terminals.iter().enumerate() {
        for &l in labs {
            label_to_term.insert(l, idx);
        }
    }

    let exterior = mesh.exterior_facets();
    let (mut z_min_all, mut z_max_all) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ext_lo, mut ext_hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for v in 0..mesh.n_points() {
        let p = mesh.point(v);
        z_min_all = z_min_all.min(p[2]);
        z_max_all = z_max_all.max(p[2]);
        for k in 0..3 {
            ext_lo[k] = ext_lo[k].min(p[k]);
            ext_hi[k] = ext_hi[k].max(p[k]);
        }
    }
    let extent = (ext_hi[0] - ext_lo[0]).max(ext_hi[1] - ext_lo[1]).max(ext_hi[2] - ext_lo[2]).max(1e-30);
    let z_tol = 1e-5 * extent;

    // Collect interior terminal facets to determine z_sheet
    let mut z_sum = 0.0;
    let mut z_cnt = 0usize;
    let mut term_facets: Vec<(usize, usize)> = Vec::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if let Some(&t_idx) = label_to_term.get(&lab) {
            if !exterior[f] {
                let fv = mesh.facet(f);
                for &v in fv {
                    z_sum += mesh.point(v as usize)[2];
                    z_cnt += 1;
                }
                term_facets.push((f, t_idx));
            }
        }
    }
    if term_facets.is_empty() {
        // Fallback: if all terminal facets lie on a single exterior plane (e.g. a 1-sided half-box mesh), allow them
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if let Some(&t_idx) = label_to_term.get(&lab) {
                let fv = mesh.facet(f);
                for &v in fv {
                    z_sum += mesh.point(v as usize)[2];
                    z_cnt += 1;
                }
                term_facets.push((f, t_idx));
            }
        }
    }
    if term_facets.is_empty() {
        return Err("planar_bem found no facets matching the terminal boundary labels".into());
    }
    let z_sheet = z_sum / (z_cnt as f64);

    // Verify all terminal facets lie on z = z_sheet
    for &(f, _) in &term_facets {
        for &v in mesh.facet(f) {
            let z = mesh.point(v as usize)[2];
            if (z - z_sheet).abs() > z_tol {
                return Err(format!(
                    "planar_bem requires coplanar terminal sheets on a single z plane (found z = {z:.3e} vs z_sheet = {z_sheet:.3e})"
                ));
            }
        }
    }

    // Also include any interior coplanar ground facets lying on z = z_sheet
    let mut coplanar_ground_facets: Vec<usize> = Vec::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if !exterior[f] && ground_labels.contains(&lab) && !label_to_term.contains_key(&lab) {
            let fv = mesh.facet(f);
            if fv.iter().all(|&v| (mesh.point(v as usize)[2] - z_sheet).abs() <= z_tol) {
                coplanar_ground_facets.push(f);
            }
        }
    }

    // Identify perimeter rim vertices of the sheet facets on z = z_sheet
    let mut edge_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for &(f, _) in &term_facets {
        let fv = mesh.facet(f);
        for &(a, b) in &[(fv[0], fv[1]), (fv[1], fv[2]), (fv[2], fv[0])] {
            *edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    for &f in &coplanar_ground_facets {
        let fv = mesh.facet(f);
        for &(a, b) in &[(fv[0], fv[1]), (fv[1], fv[2]), (fv[2], fv[0])] {
            *edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut rim_verts = vec![false; mesh.n_points()];
    for (&(a, b), &cnt) in &edge_count {
        if cnt == 1 {
            rim_verts[a as usize] = true;
            rim_verts[b as usize] = true;
        }
    }

    let mut panels = Vec::with_capacity(term_facets.len() + coplanar_ground_facets.len());
    for &(f, t_idx) in &term_facets {
        let fv = mesh.facet(f);
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        let p2 = mesh.point(fv[2] as usize);
        let is_rim = rim_verts[fv[0] as usize] || rim_verts[fv[1] as usize] || rim_verts[fv[2] as usize];
        if let Some(pan) = Panel::new([[p0[0], p0[1]], [p1[0], p1[1]], [p2[0], p2[1]]], Some(t_idx), is_rim) {
            panels.push(pan);
        }
    }
    for &f in &coplanar_ground_facets {
        let fv = mesh.facet(f);
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        let p2 = mesh.point(fv[2] as usize);
        let is_rim = rim_verts[fv[0] as usize] || rim_verts[fv[1] as usize] || rim_verts[fv[2] as usize];
        if let Some(pan) = Panel::new([[p0[0], p0[1]], [p1[0], p1[1]], [p2[0], p2[1]]], None, is_rim) {
            panels.push(pan);
        }
    }

    // Auto-detect dielectric stack above and below z_sheet from mesh domains
    let dom_eps = domain_scalar_eps(mesh, spec);
    let mut dom_z_range: BTreeMap<i32, (f64, f64, f64, usize)> = BTreeMap::new();
    for c in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[c];
        let cv = mesh.cell(c);
        let (mut z0, mut z1, mut zc) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
        for &v in cv {
            let z = mesh.point(v as usize)[2];
            z0 = z0.min(z);
            z1 = z1.max(z);
            zc += 0.25 * z;
        }
        let e = dom_z_range.entry(dom).or_insert((f64::INFINITY, f64::NEG_INFINITY, 0.0, 0));
        e.0 = e.0.min(z0);
        e.1 = e.1.max(z1);
        e.2 += zc;
        e.3 += 1;
    }

    let mut above_doms: Vec<(i32, f64, f64, f64)> = Vec::new();
    let mut below_doms: Vec<(i32, f64, f64, f64)> = Vec::new();
    for (&dom, &(z0, z1, zsum, cnt)) in &dom_z_range {
        let z_mean = zsum / (cnt.max(1) as f64);
        if z_mean > z_sheet {
            above_doms.push((dom, z0, z1, z_mean));
        } else {
            below_doms.push((dom, z0, z1, z_mean));
        }
    }
    // Sort above domains by increasing z_mean (closest to z_sheet first), below domains by decreasing z_mean
    above_doms.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    below_doms.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));

    let eps_top = above_doms.first().and_then(|(d, _, _, _)| dom_eps.get(d).copied()).unwrap_or(1.0);
    let (sub_dom_opt, z_sub_bot, eps_sub) = match below_doms.first() {
        Some(&(d, z0, _, _)) => (Some(d), z0, dom_eps.get(&d).copied().unwrap_or(1.0)),
        None => (None, z_sheet, eps_top),
    };
    let t_sub = if sub_dom_opt.is_some() { (z_sheet - z_sub_bot).max(0.0) } else { f64::INFINITY };

    // Check if the plane z = z_sub_bot carries Ground facets (and whether top z = z_max_all is also grounded)
    let mut bot_grounded = false;
    let mut top_grounded = false;
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if ground_labels.contains(&lab) {
            let fv = mesh.facet(f);
            if fv.iter().all(|&v| (mesh.point(v as usize)[2] - z_sub_bot).abs() <= z_tol) {
                bot_grounded = true;
            }
            if fv.iter().all(|&v| (mesh.point(v as usize)[2] - z_max_all).abs() <= z_tol) {
                top_grounded = true;
            }
        }
    }

    let (eps_bot, grounded_backplane, t_sub_eff) = if below_doms.len() >= 2 {
        let d_bot = below_doms[1].0;
        let e_bot = dom_eps.get(&d_bot).copied().unwrap_or(1.0);
        (e_bot, bot_grounded, t_sub)
    } else if bot_grounded && !top_grounded {
        // Explicit grounded backplane on bottom of substrate (with open top)
        (1.0, true, t_sub)
    } else if bot_grounded && top_grounded {
        // All exterior walls grounded (3D FEM box enclosure): bottom wall at z_sub_bot is at V = 0
        (1.0, true, t_sub)
    } else {
        // Ungrounded bottom: if domain name mentions halfspace/inf, treat as semi-infinite; else finite slab in eps_top
        let is_halfspace = sub_dom_opt.and_then(|d| spec.domain_names.get(&d)).is_some_and(|name| {
            let low = name.to_ascii_lowercase();
            low.contains("halfspace") || low.contains("inf")
        });
        if is_halfspace {
            (eps_sub, false, f64::INFINITY)
        } else {
            (eps_top, false, t_sub)
        }
    };

    let stack = LayeredStack { eps_top, eps_sub, eps_bot, t_sub: t_sub_eff, grounded_backplane };
    Ok((terminals, panels, stack, z_sheet))
}

#[allow(clippy::type_complexity)]
fn extract_2d_panels_and_stack(
    mesh: &Discretisation,
    spec: &Problem,
) -> Result<(Vec<(String, Vec<i32>)>, Vec<Panel>, LayeredStack), String> {
    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();
    let mut ground_labels: BTreeSet<i32> = BTreeSet::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::Terminal(name) => {
                let name = name.trim();
                if name.is_empty() {
                    return Err("terminal condition requires a non-empty 'name'".into());
                }
                if let Some((_, labs)) = terminals.iter_mut().find(|(n, _)| n == name) {
                    for &l in &bc.labels {
                        if !labs.contains(&l) {
                            labs.push(l);
                        }
                    }
                } else {
                    terminals.push((name.to_string(), bc.labels.clone()));
                }
            }
            BcKind::Ground => {
                ground_labels.extend(bc.labels.iter().copied());
            }
            BcKind::Potential(v) if v.abs() < 1e-15 => {
                ground_labels.extend(bc.labels.iter().copied());
            }
            _ => {}
        }
    }
    if terminals.is_empty() {
        return Err("planar_bem requires at least one terminal boundary condition".into());
    }
    let mut label_to_term: BTreeMap<i32, usize> = BTreeMap::new();
    for (idx, (_, labs)) in terminals.iter().enumerate() {
        for &l in labs {
            label_to_term.insert(l, idx);
        }
    }

    // Map 2D mesh vertices to boundary facet labels, then assign each 2D domain to the terminal or ground
    // touching its boundary facets.
    let mut vert_term: Vec<Option<usize>> = vec![None; mesh.n_points()];
    let mut vert_ground = vec![false; mesh.n_points()];
    let mut rim_verts = vec![false; mesh.n_points()];
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        if let Some(&t_idx) = label_to_term.get(&lab) {
            for &v in fv {
                vert_term[v as usize] = Some(t_idx);
                rim_verts[v as usize] = true;
            }
        } else if ground_labels.contains(&lab) {
            for &v in fv {
                vert_ground[v as usize] = true;
                rim_verts[v as usize] = true;
            }
        }
    }

    let mut dom_term: BTreeMap<i32, Option<usize>> = BTreeMap::new();
    for d in mesh.domain_labels() {
        let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
        let mut g_cnt = 0usize;
        for c in 0..mesh.n_cells() {
            if mesh.cell_labels[c] == d {
                for &v in mesh.cell(c) {
                    if let Some(t) = vert_term[v as usize] {
                        *counts.entry(t).or_insert(0) += 1;
                    }
                    if vert_ground[v as usize] {
                        g_cnt += 1;
                    }
                }
            }
        }
        if let Some((&best_t, _)) = counts.iter().max_by_key(|e| e.1) {
            // Only assign domain `d` to `best_t` if it does not touch multiple distinct terminals
            if counts.len() == 1 {
                dom_term.insert(d, Some(best_t));
            }
        } else if g_cnt > 0 {
            dom_term.insert(d, None);
        }
    }

    let mut panels = Vec::new();
    for c in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[c];
        if let Some(&term_opt) = dom_term.get(&dom) {
            let cv = mesh.cell(c);
            let p0 = mesh.point(cv[0] as usize);
            let p1 = mesh.point(cv[1] as usize);
            let p2 = mesh.point(cv[2] as usize);
            let is_rim = rim_verts[cv[0] as usize] || rim_verts[cv[1] as usize] || rim_verts[cv[2] as usize];
            if let Some(pan) = Panel::new([[p0[0], p0[1]], [p1[0], p1[1]], [p2[0], p2[1]]], term_opt, is_rim) {
                panels.push(pan);
            }
        }
    }

    let dom_eps = domain_scalar_eps(mesh, spec);
    let eps_eff = dom_eps.values().copied().next().unwrap_or(1.0);
    let stack = LayeredStack {
        eps_top: eps_eff,
        eps_sub: eps_eff,
        eps_bot: eps_eff,
        t_sub: f64::INFINITY,
        grounded_backplane: false,
    };
    Ok((terminals, panels, stack))
}

fn invert_spd_ec_matrix(maxwell: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = maxwell.len();
    if n == 0 {
        return Ok(vec![]);
    }
    let mut flat = Vec::with_capacity(n * n);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = dense_eigh(&flat, n)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; n]; n];
    if lam_max <= 0.0 {
        return Ok(ec);
    }
    for k in 0..n {
        if eigs[k] > 1e-13 * lam_max {
            let inv_l = scale / eigs[k];
            for i in 0..n {
                for j in 0..n {
                    ec[i][j] += inv_l * vecs[k][i] * vecs[k][j];
                }
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let s = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = s;
            ec[j][i] = s;
        }
    }
    Ok(ec)
}

fn build_tls_meta(mesh: &Discretisation, spec: &Problem) -> TlsMeta {
    let domain_eps_r = domain_scalar_eps(mesh, spec);
    let mut cond_set = BTreeSet::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Potential(_) | BcKind::Ground | BcKind::Terminal(_)) {
            cond_set.extend(bc.labels.iter().copied());
        }
    }
    let conductor_labels: Vec<i32> = cond_set.iter().copied().collect();
    let exterior = mesh.exterior_facets();
    let mut sheet_set = BTreeSet::new();
    let mut pec_set = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if cond_set.contains(&lab) {
            if exterior[f] {
                pec_set.insert(lab);
            } else {
                sheet_set.insert(lab);
            }
        }
    }
    TlsMeta {
        domain_eps_r,
        domain_names: spec.domain_names.clone(),
        conductor_labels,
        sheet_labels: sheet_set.into_iter().collect(),
        sheet_thickness_m: BTreeMap::new(),
        pec_labels: pec_set.into_iter().collect(),
    }
}

/// Solves a multi-terminal planar sheet capacitance problem with an explicit [`LayeredStack`].
pub fn solve_planar_bem_with_stack(
    mesh: &Discretisation,
    spec: &Problem,
    stack_override: Option<LayeredStack>,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;

    let (terminals, panels, detected_stack) = if mesh.dim == 3 {
        let (t, p, s, _z) = extract_3d_panels_and_stack(mesh, spec)?;
        (t, p, s)
    } else {
        extract_2d_panels_and_stack(mesh, spec)?
    };
    let stack = stack_override.unwrap_or(detected_stack);
    let n_term = terminals.len();
    let n_panels = panels.len();

    let (c_coarse, _q_coarse) = assemble_and_solve_bem(&panels, n_term, &stack)?;
    let mut log = Vec::new();
    let maxwell_matrix_f = if spec.richardson {
        let fine_panels = subdivide_panels_1_to_4(&panels);
        let n_fine = fine_panels.len();
        let (c_fine, _) = assemble_and_solve_bem(&fine_panels, n_term, &stack)?;
        let c_ext = extrapolate_spd_log_space(&c_coarse, &c_fine, 2.0)?;
        log.push(format!(
            "planar_bem (Richardson 1->4, w_1=2.0): {n_panels} -> {n_fine} sheet panels, eps_top={:.3}, eps_sub={:.3}, eps_bot={:.3}, t_sub={:.2e} m, grounded_backplane={}",
            stack.eps_top, stack.eps_sub, stack.eps_bot, stack.t_sub, stack.grounded_backplane
        ));
        c_ext
    } else {
        log.push(format!(
            "planar_bem: {n_panels} sheet panels served {n_term} terminal(s), eps_top={:.3}, eps_sub={:.3}, eps_bot={:.3}, t_sub={:.2e} m, grounded_backplane={}",
            stack.eps_top, stack.eps_sub, stack.eps_bot, stack.t_sub, stack.grounded_backplane
        ));
        c_coarse
    };

    let mut mutual_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        let mut sum_row = 0.0;
        for k in 0..n_term {
            sum_row += maxwell_matrix_f[i][k];
            if i != k {
                mutual_matrix_f[i][k] = -maxwell_matrix_f[i][k];
            }
        }
        mutual_matrix_f[i][i] = sum_row;
    }

    let ec_matrix_hz = invert_spd_ec_matrix(&maxwell_matrix_f)?;
    let tls = build_tls_meta(mesh, spec);

    let n_dofs = mesh.n_dofs(order);
    let term_dofs: Vec<Vec<usize>> = terminals.iter().map(|(_, labs)| boundary_dofs(mesh, labs, order)).collect();

    let mut solutions = Vec::with_capacity(n_term);
    let mut scalars = Vec::new();
    for j in 0..n_term {
        let mut v_vec = vec![0.0; n_dofs];
        for &d in &term_dofs[j] {
            v_vec[d] = 1.0;
        }
        let (fields, cell_fields) = electrostatics_fields(mesh, order, &v_vec);
        if j == 0 {
            let en = norm_of(&components(&fields[1].values, mesh.dim));
            let w = 0.5 * maxwell_matrix_f[0][0];
            let (vmin, vmax) = minmax(&v_vec);
            scalars.push(("stored energy W_e (J; per metre in 2D)".to_string(), Scalar::Num(w)));
            scalars.push(("max |E| (V/m)".into(), Scalar::Num(minmax(&en).1)));
            scalars.push(("V min (V)".into(), Scalar::Num(vmin)));
            scalars.push(("V max (V)".into(), Scalar::Num(vmax)));
            if n_term == 1 {
                scalars.push((
                    "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
                    Scalar::Num(maxwell_matrix_f[0][0]),
                ));
            }
        }
        solutions.push(Solution {
            label: format!("terminal {}", terminals[j].0),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(maxwell_matrix_f[j][j]),
            ..Default::default()
        });
    }

    for i in 0..n_term {
        for j in 0..n_term {
            let name_key = format!("C_{}_{} (F)", terminals[i].0, terminals[j].0);
            let idx_key = format!("C_{i}_{j} (F)");
            scalars.push((name_key.clone(), Scalar::Num(maxwell_matrix_f[i][j])));
            if idx_key != name_key {
                scalars.push((idx_key, Scalar::Num(maxwell_matrix_f[i][j])));
            }
        }
    }
    for i in 0..n_term {
        for j in 0..n_term {
            let name_key = format!("E_C_{}_{} (Hz)", terminals[i].0, terminals[j].0);
            let idx_key = format!("E_C_{i}_{j} (Hz)");
            scalars.push((name_key.clone(), Scalar::Num(ec_matrix_hz[i][j])));
            if idx_key != name_key {
                scalars.push((idx_key, Scalar::Num(ec_matrix_hz[i][j])));
            }
        }
    }

    let cap_data = CapacitanceMatrixData {
        terminals: terminals.into_iter().map(|(n, _)| n).collect(),
        maxwell_matrix_f,
        mutual_matrix_f,
        ec_matrix_hz,
    };
    scalars.push(cap_data.scalar());
    scalars.push(tls.scalar());

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n_panels,
    })
}

/// Solves a multi-terminal planar sheet capacitance problem by auto-detecting the layered stack
/// from the 3D `sheet3d` or 2D mesh.
pub fn solve_planar_bem(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    solve_planar_bem_with_stack(mesh, spec, None)
}
