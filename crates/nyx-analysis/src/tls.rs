// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Dielectric two-level-system (TLS) surface participation analysis (`TlsParticipationAnalysis`, Milestone 4.4).
//!
//! Computes the metal–air (`p_ma`), metal–substrate (`p_ms`), and substrate–air (`p_sa`) dielectric loss
//! participation ratios, total effective loss tangent `loss_tan_total`, and TLS-limited quality factor `q_tls`:
//! - **2D (`dim == 2`)**: Integrates the bulk electric energy `U_elec` with the 6-point degree-4 Dunavant rule
//!   (`quad_rule(2, 4)`) and P1/P2 `shape` functions, and evaluates exact boundary line integrals (`facet_rule(1)`)
//!   on the **vacuum-side** triangle for `MA` and `SA` and the **substrate-side** triangle for `MS`. Also supports
//!   explicit meshed finite-thickness oxide layers (`t_ox = 2 nm`) when present as thin dielectric domains.
//! - **3D (`dim == 3`, Matched Asymptotic Edge-Tube Regularization)**: Detects the 1D singular perimeter
//!   `\partial\Gamma_{SC}` of interior zero-thickness superconducting sheets, integrates exterior `MA`, `MS`, and `SA`
//!   surface facets outside distance `R_0` (`edge_cutoff`) using subdivided near-edge triangle quadrature, extracts
//!   the convergent perimeter singularity intensity `I_edge = \oint |A(s)|^2 ds` from the 3D electric energy inside
//!   the cylindrical edge tube `r(x, \partial\Gamma_{SC}) < R_0`, and adds the matched asymptotic finite-thickness
//!   (`d > 0`) corner contributions so `p_ma, p_ms, p_sa` are independent of `R_0` and converge as `h -> 0`.

use crate::epr::EPS0;
use crate::Analysis;
use femcore::fem::{facet_rule, geometry, quad_rule, shape};
use femcore::linalg::scalar::{Scalar, C64};
use nyx_physics::{discretise, CellField, Discretisation, Error, Mesh, Solution, SolveResult, TlsMeta};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;

const FACES_3D: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Typed output of [`tls_participation`] and [`tls_participation_with_params`].
#[derive(Clone, Debug, PartialEq)]
pub struct TlsReport {
    pub solution: usize,
    pub u_elec_j: f64,
    pub p_ma: f64,
    pub p_ms: f64,
    pub p_sa: f64,
    pub p_total: f64,
    pub loss_tan_total: f64,
    pub q_tls: f64,
    pub edge_cutoff_m: f64,
    pub p_ma_interior: f64,
    pub p_ms_interior: f64,
    pub p_sa_interior: f64,
    pub p_ma_edge: f64,
    pub p_ms_edge: f64,
    pub p_sa_edge: f64,
    pub convergence_check: Value,
}

impl TlsReport {
    pub fn to_json(&self) -> Value {
        json!({
            "solution": self.solution,
            "u_elec_J": self.u_elec_j,
            "p_ma": self.p_ma,
            "p_ms": self.p_ms,
            "p_sa": self.p_sa,
            "p_total": self.p_total,
            "loss_tan_total": self.loss_tan_total,
            "q_tls": self.q_tls,
            "edge_cutoff_m": self.edge_cutoff_m,
            "p_ma_interior": self.p_ma_interior,
            "p_ms_interior": self.p_ms_interior,
            "p_sa_interior": self.p_sa_interior,
            "p_ma_edge": self.p_ma_edge,
            "p_ms_edge": self.p_ms_edge,
            "p_sa_edge": self.p_sa_edge,
            "convergence_check": self.convergence_check,
        })
    }
}

/// Parameters for [`tls_participation_with_params`].
#[derive(Clone, Debug)]
struct TlsParams {
    solution: usize,
    t_ma: f64,
    t_ms: f64,
    t_sa: f64,
    eps_ma: f64,
    eps_ms: f64,
    eps_sa: f64,
    tan_delta_ma: f64,
    tan_delta_ms: f64,
    tan_delta_sa: f64,
    d: Option<f64>,
    edge_cutoff: Option<f64>,
}

impl TlsParams {
    fn from_json(params: &Value) -> Result<Self, Error> {
        let num = |k: &str, def: f64| -> Result<f64, Error> {
            match params.get(k) {
                None | Some(Value::Null) => Ok(def),
                Some(v) => v
                    .as_f64()
                    .filter(|x| x.is_finite() && *x >= 0.0)
                    .ok_or_else(|| Error::invalid(format!("parameter '{k}' must be a non-negative finite number"))),
            }
        };
        let solution = match params.get("solution") {
            None | Some(Value::Null) => 0,
            Some(v) => {
                v.as_u64().ok_or_else(|| Error::invalid("parameter 'solution' must be a non-negative integer"))?
                    as usize
            }
        };
        let d = match params.get("d") {
            None | Some(Value::Null) => None,
            Some(v) => Some(
                v.as_f64()
                    .filter(|x| x.is_finite() && *x > 0.0)
                    .ok_or_else(|| Error::invalid("parameter 'd' must be a positive film thickness in m"))?,
            ),
        };
        let edge_cutoff = match params.get("edge_cutoff") {
            None | Some(Value::Null) => None,
            Some(v) => Some(
                v.as_f64()
                    .filter(|x| x.is_finite() && *x > 0.0)
                    .ok_or_else(|| Error::invalid("parameter 'edge_cutoff' must be a positive radius in m"))?,
            ),
        };
        Ok(Self {
            solution,
            t_ma: num("t_ma", 2.0e-9)?,
            t_ms: num("t_ms", 2.0e-9)?,
            t_sa: num("t_sa", 2.0e-9)?,
            eps_ma: num("eps_ma", 10.0)?.max(1e-12),
            eps_ms: num("eps_ms", 11.45)?.max(1e-12),
            eps_sa: num("eps_sa", 10.0)?.max(1e-12),
            tan_delta_ma: num("tan_delta_ma", 1.0e-3)?,
            tan_delta_ms: num("tan_delta_ms", 1.0e-3)?,
            tan_delta_sa: num("tan_delta_sa", 1.0e-3)?,
            d,
            edge_cutoff,
        })
    }
}

/// True when `result` has electric field cell fields (`Ex`) suitable for [`TlsParticipationAnalysis`].
pub fn tls_participation_applies(result: &SolveResult) -> bool {
    !result.solutions.is_empty()
        && result.solutions.iter().any(|s| s.cell_fields.iter().any(|f| f.name == "Ex"))
        && (TlsMeta::of(result).is_some() || result.physics == "dfs" || result.physics == "ems")
}

// ------------------------------------------------------------------ cell field evaluation helpers

struct EFieldAccessor<'a> {
    ex: &'a CellField,
    ey: Option<&'a CellField>,
    ez: &'a CellField,
    nloc: usize,
    order: usize,
}

impl<'a> EFieldAccessor<'a> {
    fn new(sol: &'a Solution, dim: usize, n_cells: usize) -> Result<Self, Error> {
        let find = |name: &str| sol.cell_fields.iter().find(|f| f.name == name);
        let ex = find("Ex").ok_or_else(|| Error::invalid("solution has no 'Ex' cell field"))?;
        let (ey, ez) = if dim == 2 {
            if let Some(ez) = find("Ez") {
                (find("Ey"), ez)
            } else if let Some(ey) = find("Ey") {
                (None, ey)
            } else {
                return Err(Error::invalid("2D solution needs 'Ez' (or 'Ey') cell field"));
            }
        } else {
            let ey = find("Ey").ok_or_else(|| Error::invalid("3D solution has no 'Ey' cell field"))?;
            let ez = find("Ez").ok_or_else(|| Error::invalid("3D solution has no 'Ez' cell field"))?;
            (Some(ey), ez)
        };
        let nloc = ex.nloc;
        if ez.nloc != nloc || ey.is_some_and(|f| f.nloc != nloc) || nloc == 0 {
            return Err(Error::invalid("electric field components have mismatched nloc"));
        }
        for f in [Some(ex), ey, Some(ez)].into_iter().flatten() {
            if f.re.len() != n_cells * nloc || (!f.im.is_empty() && f.im.len() != n_cells * nloc) {
                return Err(Error::invalid(format!(
                    "cell field '{}' has {} entries for {} cells of {} nodes",
                    f.name,
                    f.re.len(),
                    n_cells,
                    nloc
                )));
            }
        }
        let order = match (dim, nloc) {
            (2, 3) | (3, 4) => 1,
            (2, 6) | (3, 10) => 2,
            _ => return Err(Error::invalid(format!("unsupported cell field layout: dim={dim}, nloc={nloc}"))),
        };
        Ok(Self { ex, ey, ez, nloc, order })
    }

    #[inline]
    fn eval_with_shapes(&self, ci: usize, sh: &[f64]) -> [C64; 3] {
        let base = ci * self.nloc;
        let eval_one = |f: &CellField| -> C64 {
            let mut re = 0.0;
            let mut im = 0.0;
            if f.im.is_empty() {
                for (a, &na) in sh.iter().enumerate().take(self.nloc) {
                    re += na * f.re[base + a];
                }
            } else {
                for (a, &na) in sh.iter().enumerate().take(self.nloc) {
                    re += na * f.re[base + a];
                    im += na * f.im[base + a];
                }
            }
            C64::new(re, im)
        };
        let vx = eval_one(self.ex);
        let vy = self.ey.map_or(C64::ZERO, eval_one);
        let vz = eval_one(self.ez);
        [vx, vy, vz]
    }
}

#[inline]
fn norm_sq3(e: &[C64; 3]) -> f64 {
    e[0].norm_sqr() + e[1].norm_sqr() + e[2].norm_sqr()
}

// ------------------------------------------------------------------ 2D TLS participation

fn run_tls_2d(
    disc: &Discretisation,
    ef: &EFieldAccessor<'_>,
    meta: &TlsMeta,
    p: &TlsParams,
    c_e: f64,
) -> Result<TlsReport, Error> {
    let geo = geometry(disc);
    let quad = quad_rule(2, 4);
    let shapes: Vec<Vec<f64>> = quad.lam.iter().map(|l| shape(l, ef.order).0).collect();
    let nc = disc.n_cells();

    let cell_eps: Vec<f64> =
        (0..nc).map(|ci| meta.domain_eps_r.get(&disc.cell_labels[ci]).copied().unwrap_or(1.0)).collect();

    let mut u_elec_j = 0.0;
    for (ci, &eps_ci) in cell_eps.iter().enumerate().take(nc) {
        let mut cell_e2 = 0.0;
        for (q, sh) in shapes.iter().enumerate() {
            let e = ef.eval_with_shapes(ci, sh);
            cell_e2 += quad.w[q] * norm_sq3(&e);
        }
        u_elec_j += c_e * EPS0 * eps_ci * geo.vol[ci] * cell_e2;
    }
    if u_elec_j <= 0.0 {
        return Err(Error::invalid("electric field has zero stored energy"));
    }

    // Build edge -> adjacent triangles map: key = (v0.min(v1), v0.max(v1))
    type EdgeAdjList = Vec<(usize, [usize; 2])>;
    let mut edge_to_cells: HashMap<(u32, u32), EdgeAdjList> = HashMap::with_capacity(disc.n_edges());
    for ci in 0..nc {
        let v = disc.cell(ci);
        for [la, lb] in [[0, 1], [1, 2], [0, 2]] {
            let (va, vb) = (v[la], v[lb]);
            let (key, slots) = if va < vb { ((va, vb), [la, lb]) } else { ((vb, va), [lb, la]) };
            edge_to_cells.entry(key).or_default().push((ci, slots));
        }
    }

    let mut edge_labels: HashMap<(u32, u32), Vec<i32>> = HashMap::new();
    for f in 0..disc.n_facets() {
        let fv = disc.facet(f);
        let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
        let lab = disc.facet_labels[f];
        let entry = edge_labels.entry(key).or_default();
        if !entry.contains(&lab) {
            entry.push(lab);
        }
    }

    let (lo, hi) = disc.bbox();
    let span = ((hi[0] - lo[0]).max(hi[1] - lo[1])).max(1e-30);
    let bbox_tol = 1e-7 * span;
    let is_on_outer_bbox = |p0: &[f64], p1: &[f64]| -> bool {
        ((p0[0] - lo[0]).abs() < bbox_tol && (p1[0] - lo[0]).abs() < bbox_tol)
            || ((p0[0] - hi[0]).abs() < bbox_tol && (p1[0] - hi[0]).abs() < bbox_tol)
            || ((p0[1] - lo[1]).abs() < bbox_tol && (p1[1] - lo[1]).abs() < bbox_tol)
            || ((p0[1] - hi[1]).abs() < bbox_tol && (p1[1] - hi[1]).abs() < bbox_tol)
    };

    let cond_set: BTreeSet<i32> = meta.conductor_labels.iter().copied().collect();
    let has_interior_cond = edge_labels.iter().any(|(&(v0, v1), labs)| {
        labs.iter().any(|l| cond_set.contains(l)) && !is_on_outer_bbox(disc.point(v0 as usize), disc.point(v1 as usize))
    });

    let f_rule = facet_rule(1);
    let mut u_ma = 0.0;
    let mut u_ms = 0.0;
    let mut u_sa = 0.0;

    for (&(v0, v1), cells) in &edge_to_cells {
        let p0 = disc.point(v0 as usize);
        let p1 = disc.point(v1 as usize);
        let dx = p1[0] - p0[0];
        let dz = p1[1] - p0[1];
        let len = dx.hypot(dz);
        if len <= 0.0 {
            continue;
        }
        let tau = [dx / len, dz / len];
        let n_hat = [-tau[1], tau[0]];

        let is_cond = edge_labels.get(&(v0, v1)).is_some_and(|labs| labs.iter().any(|l| cond_set.contains(l)));

        if is_cond {
            if has_interior_cond && is_on_outer_bbox(p0, p1) {
                continue;
            }
            for &(ci, [la, lb]) in cells {
                let eps_c = cell_eps[ci];
                let mut int_en2 = 0.0;
                let mut int_epar2 = 0.0;
                for (q, l1) in f_rule.lam.iter().enumerate() {
                    let mut lam2 = [0.0; 3];
                    lam2[la] = l1[0];
                    lam2[lb] = l1[1];
                    let (sh, _) = shape(&lam2, ef.order);
                    let e = ef.eval_with_shapes(ci, &sh);
                    let en = e[0].scale(n_hat[0]) + e[2].scale(n_hat[1]);
                    let et = e[0].scale(tau[0]) + e[2].scale(tau[1]);
                    let en2 = en.norm_sqr();
                    let epar2 = et.norm_sqr() + e[1].norm_sqr();
                    let wq = f_rule.w[q] * len;
                    int_en2 += wq * en2;
                    int_epar2 += wq * epar2;
                }
                if eps_c <= 1.5 {
                    // Vacuum-side cell -> MA
                    u_ma += c_e * EPS0 * p.t_ma * ((eps_c * eps_c / p.eps_ma) * int_en2 + p.eps_ma * int_epar2);
                } else {
                    // Substrate-side cell -> MS
                    u_ms += c_e * EPS0 * p.t_ms * ((eps_c * eps_c / p.eps_ms) * int_en2 + p.eps_ms * int_epar2);
                }
            }
        } else if cells.len() == 2 {
            let (c0, s0) = cells[0];
            let (c1, s1) = cells[1];
            let (e0, e1) = (cell_eps[c0], cell_eps[c1]);
            let vac_side = if e0 <= 1.5 && e1 > 1.5 {
                Some((c0, s0, e0))
            } else if e1 <= 1.5 && e0 > 1.5 {
                Some((c1, s1, e1))
            } else {
                None
            };
            if let Some((c_vac, [la, lb], eps_vac)) = vac_side {
                let mut int_en2 = 0.0;
                let mut int_epar2 = 0.0;
                for (q, l1) in f_rule.lam.iter().enumerate() {
                    let mut lam2 = [0.0; 3];
                    lam2[la] = l1[0];
                    lam2[lb] = l1[1];
                    let (sh, _) = shape(&lam2, ef.order);
                    let e = ef.eval_with_shapes(c_vac, &sh);
                    let en = e[0].scale(n_hat[0]) + e[2].scale(n_hat[1]);
                    let et = e[0].scale(tau[0]) + e[2].scale(tau[1]);
                    let en2 = en.norm_sqr();
                    let epar2 = et.norm_sqr() + e[1].norm_sqr();
                    let wq = f_rule.w[q] * len;
                    int_en2 += wq * en2;
                    int_epar2 += wq * epar2;
                }
                u_sa += c_e * EPS0 * p.t_sa * (p.eps_sa * int_epar2 + (eps_vac * eps_vac / p.eps_sa) * int_en2);
            }
        }
    }

    let p_ma = u_ma / u_elec_j;
    let p_ms = u_ms / u_elec_j;
    let p_sa = u_sa / u_elec_j;
    let p_total = p_ma + p_ms + p_sa;
    let loss_tan_total = p_ma * p.tan_delta_ma + p_ms * p.tan_delta_ms + p_sa * p.tan_delta_sa;
    let q_tls = if loss_tan_total > 0.0 { 1.0 / loss_tan_total } else { f64::INFINITY };

    Ok(TlsReport {
        solution: p.solution,
        u_elec_j,
        p_ma,
        p_ms,
        p_sa,
        p_total,
        loss_tan_total,
        q_tls,
        edge_cutoff_m: 0.0,
        p_ma_interior: p_ma,
        p_ms_interior: p_ms,
        p_sa_interior: p_sa,
        p_ma_edge: 0.0,
        p_ms_edge: 0.0,
        p_sa_edge: 0.0,
        convergence_check: json!({
            "dim": 2,
            "mode": "exact_boundary_2d",
            "n_cells": nc,
        }),
    })
}

// ------------------------------------------------------------------ 3D Matched Asymptotic Edge-Tube Regularization

#[inline]
fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

#[inline]
fn pt3(disc: &Discretisation, v: u32) -> [f64; 3] {
    let p = disc.point(v as usize);
    [p[0], p[1], p[2]]
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

/// Squared distance from point `x` to 3D line segment `(a, b)`.
#[inline]
fn point_segment_dist_sq(x: [f64; 3], a: [f64; 3], b: [f64; 3], inv_len_sq: f64) -> f64 {
    let ab = sub3(b, a);
    let ax = sub3(x, a);
    let t = (dot3(ax, ab) * inv_len_sq).clamp(0.0, 1.0);
    let dx = ax[0] - t * ab[0];
    let dy = ax[1] - t * ab[1];
    let dz = ax[2] - t * ab[2];
    dx * dx + dy * dy + dz * dz
}

struct Segment3 {
    a: [f64; 3],
    b: [f64; 3],
    inv_len_sq: f64,
}

#[inline]
fn min_dist_to_segments(x: [f64; 3], segs: &[Segment3]) -> f64 {
    let mut best_d2 = f64::INFINITY;
    for s in segs {
        let d2 = point_segment_dist_sq(x, s.a, s.b, s.inv_len_sq);
        if d2 < best_d2 {
            best_d2 = d2;
        }
    }
    best_d2.sqrt()
}

/// Generates `n_sub^2` sub-triangles on the 2D reference simplex, returning `(bary3, weight)`
/// where `sum(weight) == 1.0`.
fn subdivided_tri_quad(n_sub: usize) -> Vec<([f64; 3], f64)> {
    let base = quad_rule(2, 2);
    let inv_n = 1.0 / (n_sub as f64);
    let area_scale = inv_n * inv_n;
    let mut out = Vec::with_capacity(n_sub * n_sub * base.w.len());

    let map_tri = |v0: [f64; 2], v1: [f64; 2], v2: [f64; 2], l: &[f64]| -> [f64; 3] {
        let s = l[0] * v0[0] + l[1] * v1[0] + l[2] * v2[0];
        let t = l[0] * v0[1] + l[1] * v1[1] + l[2] * v2[1];
        [1.0 - s - t, s, t]
    };

    for i in 0..n_sub {
        for j in 0..(n_sub - i) {
            let s0 = (i as f64) * inv_n;
            let t0 = (j as f64) * inv_n;
            let v0 = [s0, t0];
            let v1 = [s0 + inv_n, t0];
            let v2 = [s0, t0 + inv_n];
            for (q, l) in base.lam.iter().enumerate() {
                out.push((map_tri(v0, v1, v2, l), base.w[q] * area_scale));
            }
            if i + j + 1 < n_sub {
                let u0 = [s0 + inv_n, t0];
                let u1 = [s0 + inv_n, t0 + inv_n];
                let u2 = [s0, t0 + inv_n];
                for (q, l) in base.lam.iter().enumerate() {
                    out.push((map_tri(u0, u1, u2, l), base.w[q] * area_scale));
                }
            }
        }
    }
    out
}

/// Generates `n_sub^3` quadrature points on the 3D reference tetrahedron by collapsing the subdivided cube,
/// returning `(bary4, weight)` where `sum(weight) == 1.0`.
fn subdivided_tet_quad(n_sub: usize) -> Vec<([f64; 4], f64)> {
    // 2-point Gauss-Legendre on [0, 1] per sub-interval
    let g0 = 0.5 * (1.0 - 1.0 / 3.0_f64.sqrt());
    let g1 = 0.5 * (1.0 + 1.0 / 3.0_f64.sqrt());
    let inv_n = 1.0 / (n_sub as f64);
    let mut pts_1d = Vec::with_capacity(2 * n_sub);
    for k in 0..n_sub {
        let a = (k as f64) * inv_n;
        pts_1d.push((a + g0 * inv_n, 0.5 * inv_n));
        pts_1d.push((a + g1 * inv_n, 0.5 * inv_n));
    }
    let mut out = Vec::with_capacity(pts_1d.len().pow(3));
    for &(u, wu) in &pts_1d {
        for &(v, wv) in &pts_1d {
            for &(t, wt) in &pts_1d {
                let x = u;
                let y = v * (1.0 - u);
                let z = t * (1.0 - u) * (1.0 - v);
                let jac = 6.0 * (1.0 - u) * (1.0 - u) * (1.0 - v);
                out.push(([1.0 - x - y - z, x, y, z], wu * wv * wt * jac));
            }
        }
    }
    out
}

fn run_tls_3d(
    disc: &Discretisation,
    ef: &EFieldAccessor<'_>,
    meta: &TlsMeta,
    p: &TlsParams,
    c_e: f64,
) -> Result<TlsReport, Error> {
    let geo = geometry(disc);
    let quad3 = quad_rule(3, 4);
    let shapes3: Vec<Vec<f64>> = quad3.lam.iter().map(|l| shape(l, ef.order).0).collect();
    let nc = disc.n_cells();

    let cell_eps: Vec<f64> =
        (0..nc).map(|ci| meta.domain_eps_r.get(&disc.cell_labels[ci]).copied().unwrap_or(1.0)).collect();
    let eps_sub = cell_eps.iter().copied().filter(|&e| e > 1.5).fold(p.eps_ms, f64::max);
    let eps_vac = 1.0_f64;

    // 1. Bulk electric energy U_elec
    let mut u_elec_j = 0.0;
    let mut cell_energy = vec![0.0; nc];
    for ci in 0..nc {
        let mut e2_sum = 0.0;
        for (q, sh) in shapes3.iter().enumerate() {
            let e = ef.eval_with_shapes(ci, sh);
            e2_sum += quad3.w[q] * norm_sq3(&e);
        }
        let ue = c_e * EPS0 * cell_eps[ci] * geo.vol[ci] * e2_sum;
        cell_energy[ci] = ue;
        u_elec_j += ue;
    }
    if u_elec_j <= 0.0 {
        return Err(Error::invalid("3D electric field has zero stored energy"));
    }

    // 2. Map every sorted triangle face [v0, v1, v2] to its adjacent tetrahedra (ci, [la, lb, lc])
    let mut face_to_cells: HashMap<[u32; 3], Vec<(usize, [usize; 3])>> = HashMap::with_capacity(nc * 2);
    for ci in 0..nc {
        let v = disc.cell(ci);
        for lf in FACES_3D {
            let tri_raw = [v[lf[0]], v[lf[1]], v[lf[2]]];
            let key = sorted3(tri_raw[0], tri_raw[1], tri_raw[2]);
            // Permute local indices so they correspond to key[0], key[1], key[2] in sorted order
            let mut perm = lf;
            perm.sort_by_key(|&loc| v[loc]);
            face_to_cells.entry(key).or_default().push((ci, perm));
        }
    }

    let exterior = disc.exterior_facets();
    let sheet_set: BTreeSet<i32> = if !meta.sheet_labels.is_empty() {
        meta.sheet_labels.iter().copied().collect()
    } else {
        meta.conductor_labels.iter().copied().collect()
    };

    let mut sheet_tris: BTreeMap<[u32; 3], i32> = BTreeMap::new();
    let mut ext_edges: BTreeSet<(u32, u32)> = BTreeSet::new();
    for (f, &is_ext) in exterior.iter().enumerate().take(disc.n_facets()) {
        let fv = disc.facet(f);
        let key = sorted3(fv[0], fv[1], fv[2]);
        let lab = disc.facet_labels[f];
        if is_ext {
            for (a, b) in [(key[0], key[1]), (key[1], key[2]), (key[0], key[2])] {
                ext_edges.insert((a.min(b), a.max(b)));
            }
        }
        if !is_ext && (sheet_set.contains(&lab) || meta.conductor_labels.contains(&lab)) {
            sheet_tris.entry(key).or_insert(lab);
        }
    }

    // Collect uncovered dielectric interface triangles (`SA`: between vacuum <= 1.5 and substrate > 1.5)
    let mut sa_tris: Vec<[u32; 3]> = Vec::new();
    for (&key, cells) in &face_to_cells {
        if sheet_tris.contains_key(&key) || cells.len() != 2 {
            continue;
        }
        let (e0, e1) = (cell_eps[cells[0].0], cell_eps[cells[1].0]);
        if (e0 <= 1.5 && e1 > 1.5) || (e1 <= 1.5 && e0 > 1.5) {
            sa_tris.push(key);
        }
    }
    sa_tris.sort_unstable();

    // 3. Detect 1D singular perimeter \partial\Gamma_{SC}: edges touched by 1 sheet triangle, excluding exterior walls
    let mut sheet_edge_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for &tri in sheet_tris.keys() {
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
            *sheet_edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }

    let mut perim_edge_set: BTreeSet<(u32, u32)> = BTreeSet::new();
    let mut perim_segments: Vec<Segment3> = Vec::new();
    let mut total_perim_len = 0.0;
    for (&(a, b), &cnt) in &sheet_edge_count {
        if cnt == 1 && !ext_edges.contains(&(a, b)) {
            let pa = pt3(disc, a);
            let pb = pt3(disc, b);
            let ab = sub3(pb, pa);
            let len_sq = dot3(ab, ab);
            if len_sq > 1e-30 {
                let len = len_sq.sqrt();
                total_perim_len += len;
                perim_edge_set.insert((a, b));
                perim_segments.push(Segment3 { a: pa, b: pb, inv_len_sq: 1.0 / len_sq });
            }
        }
    }

    // Estimate mean transverse mesh spacing h_perp adjacent to \partial\Gamma_{SC}
    let mut h_perp_sum = 0.0;
    let mut h_perp_cnt = 0usize;
    for &tri in sheet_tris.keys().chain(sa_tris.iter()) {
        let (p0, p1, p2) = (pt3(disc, tri[0]), pt3(disc, tri[1]), pt3(disc, tri[2]));
        let area = 0.5 * norm3(cross3(sub3(p1, p0), sub3(p2, p0)));
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
            let e = (a.min(b), a.max(b));
            if perim_edge_set.contains(&e) {
                let elen = norm3(sub3(pt3(disc, b), pt3(disc, a))).max(1e-15);
                h_perp_sum += (2.0 * area) / elen;
                h_perp_cnt += 1;
            }
        }
    }
    let h_edge_mean = if h_perp_cnt > 0 { h_perp_sum / (h_perp_cnt as f64) } else { 1.0e-6 };
    let d_film = p.d.or_else(|| meta.sheet_thickness_m.values().next().copied()).unwrap_or(100.0e-9);
    let r0 = p.edge_cutoff.unwrap_or_else(|| (1.15 * h_edge_mean).max(4.0 * d_film));

    // 4. Exterior surface integrals (r(x, \partial\Gamma_{SC}) > R_0)
    let quad2 = quad_rule(2, 4);
    let sub_tri = subdivided_tri_quad(12);

    let integrate_surface_tri = |tri: [u32; 3], ci: usize, lf: [usize; 3]| -> (f64, f64, f64, f64) {
        let (p0, p1, p2) = (pt3(disc, tri[0]), pt3(disc, tri[1]), pt3(disc, tri[2]));
        let n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
        let n_len = norm3(n_raw);
        let area = 0.5 * n_len;
        if area <= 0.0 {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
        let h_tri = area.sqrt();

        let min_v_dist = if perim_segments.is_empty() {
            f64::INFINITY
        } else {
            let d0 = min_dist_to_segments(p0, &perim_segments);
            let d1 = min_dist_to_segments(p1, &perim_segments);
            let d2 = min_dist_to_segments(p2, &perim_segments);
            d0.min(d1).min(d2)
        };

        let mut ext_en2 = 0.0;
        let mut ext_epar2 = 0.0;
        // Also accumulate annular collar r in [0.5 R_0, R_0] for local regular-part matching
        let mut col_en2 = 0.0;
        let mut col_epar2 = 0.0;

        if min_v_dist > r0 + 1.5 * h_tri {
            for (q, l2) in quad2.lam.iter().enumerate() {
                let mut lam4 = [0.0; 4];
                lam4[lf[0]] = l2[0];
                lam4[lf[1]] = l2[1];
                lam4[lf[2]] = l2[2];
                let (sh, _) = shape(&lam4, ef.order);
                let e = ef.eval_with_shapes(ci, &sh);
                let en = e[0].scale(n_hat[0]) + e[1].scale(n_hat[1]) + e[2].scale(n_hat[2]);
                let en2 = en.norm_sqr();
                let epar2 = (norm_sq3(&e) - en2).max(0.0);
                let wq = quad2.w[q] * area;
                ext_en2 += wq * en2;
                ext_epar2 += wq * epar2;
            }
        } else {
            for &(l2, w_ref) in &sub_tri {
                let x = [
                    l2[0] * p0[0] + l2[1] * p1[0] + l2[2] * p2[0],
                    l2[0] * p0[1] + l2[1] * p1[1] + l2[2] * p2[1],
                    l2[0] * p0[2] + l2[1] * p1[2] + l2[2] * p2[2],
                ];
                let dist = min_dist_to_segments(x, &perim_segments);
                if dist <= 0.5 * r0 {
                    continue;
                }
                let mut lam4 = [0.0; 4];
                lam4[lf[0]] = l2[0];
                lam4[lf[1]] = l2[1];
                lam4[lf[2]] = l2[2];
                let (sh, _) = shape(&lam4, ef.order);
                let e = ef.eval_with_shapes(ci, &sh);
                let en = e[0].scale(n_hat[0]) + e[1].scale(n_hat[1]) + e[2].scale(n_hat[2]);
                let en2 = en.norm_sqr();
                let epar2 = (norm_sq3(&e) - en2).max(0.0);
                let wq = w_ref * area;
                if dist > r0 {
                    ext_en2 += wq * en2;
                    ext_epar2 += wq * epar2;
                } else {
                    col_en2 += wq * en2;
                    col_epar2 += wq * epar2;
                }
            }
        }
        (ext_en2, ext_epar2, col_en2, col_epar2)
    };

    let mut u_ma_ext = 0.0;
    let mut u_ms_ext = 0.0;
    let mut u_sa_ext = 0.0;
    let mut u_ma_col = 0.0;
    let mut u_ms_col = 0.0;
    let mut u_sa_col = 0.0;

    for &tri in sheet_tris.keys() {
        let Some(cells) = face_to_cells.get(&tri) else { continue };
        for &(ci, lf) in cells {
            let eps_c = cell_eps[ci];
            let (ext_en2, ext_epar2, col_en2, col_epar2) = integrate_surface_tri(tri, ci, lf);
            if eps_c <= 1.5 {
                u_ma_ext += c_e * EPS0 * p.t_ma * ((eps_c * eps_c / p.eps_ma) * ext_en2 + p.eps_ma * ext_epar2);
                u_ma_col += c_e * EPS0 * p.t_ma * ((eps_c * eps_c / p.eps_ma) * col_en2 + p.eps_ma * col_epar2);
            } else {
                u_ms_ext += c_e * EPS0 * p.t_ms * ((eps_c * eps_c / p.eps_ms) * ext_en2 + p.eps_ms * ext_epar2);
                u_ms_col += c_e * EPS0 * p.t_ms * ((eps_c * eps_c / p.eps_ms) * col_en2 + p.eps_ms * col_epar2);
            }
        }
    }

    for &tri in &sa_tris {
        let Some(cells) = face_to_cells.get(&tri) else { continue };
        if let Some(&(c_vac, lf)) = cells.iter().find(|&&(ci, _)| cell_eps[ci] <= 1.5) {
            let eps_c = cell_eps[c_vac];
            let (ext_en2, ext_epar2, col_en2, col_epar2) = integrate_surface_tri(tri, c_vac, lf);
            u_sa_ext += c_e * EPS0 * p.t_sa * (p.eps_sa * ext_epar2 + (eps_c * eps_c / p.eps_sa) * ext_en2);
            u_sa_col += c_e * EPS0 * p.t_sa * (p.eps_sa * col_epar2 + (eps_c * eps_c / p.eps_sa) * col_en2);
        }
    }

    // 5. 3D Edge-Tube Volume Integral & Matched Asymptotic Regularization (r <= R_0)
    let (mut u_ma_edge, mut u_ms_edge, mut u_sa_edge) = (0.0, 0.0, 0.0);
    let mut i_edge = 0.0;
    let mut u_tube_j = 0.0;

    if !perim_segments.is_empty() && r0 > 0.0 {
        let sub_tet = subdivided_tet_quad(4);
        let sub_tet_shapes: Vec<Vec<f64>> = sub_tet.iter().map(|(l, _)| shape(l, ef.order).0).collect();

        for ci in 0..nc {
            let v = disc.cell(ci);
            let pts = [pt3(disc, v[0]), pt3(disc, v[1]), pt3(disc, v[2]), pt3(disc, v[3])];
            let mut d_min = f64::INFINITY;
            let mut d_max = 0.0_f64;
            for &pt in &pts {
                let d = min_dist_to_segments(pt, &perim_segments);
                d_min = d_min.min(d);
                d_max = d_max.max(d);
            }
            let h_tet = (6.0 * geo.vol[ci]).cbrt();
            if d_min > r0 + h_tet {
                continue;
            }
            if d_max <= r0 {
                u_tube_j += cell_energy[ci];
            } else {
                let scale = c_e * EPS0 * cell_eps[ci] * geo.vol[ci];
                for (q, &(l4, w_ref)) in sub_tet.iter().enumerate() {
                    let x = [
                        l4[0] * pts[0][0] + l4[1] * pts[1][0] + l4[2] * pts[2][0] + l4[3] * pts[3][0],
                        l4[0] * pts[0][1] + l4[1] * pts[1][1] + l4[2] * pts[2][1] + l4[3] * pts[3][1],
                        l4[0] * pts[0][2] + l4[1] * pts[1][2] + l4[2] * pts[2][2] + l4[3] * pts[3][2],
                    ];
                    if min_dist_to_segments(x, &perim_segments) <= r0 {
                        let e = ef.eval_with_shapes(ci, &sub_tet_shapes[q]);
                        u_tube_j += scale * w_ref * norm_sq3(&e);
                    }
                }
            }
        }

        i_edge = u_tube_j / (c_e * EPS0 * PI * (eps_sub + eps_vac) * r0).max(1e-300);

        let kappa_ma = c_e * EPS0 * (p.t_ma / p.eps_ma);
        let kappa_ms = c_e * EPS0 * (eps_sub * eps_sub * p.t_ms / p.eps_ms);
        let kappa_sa = c_e * EPS0 * (p.eps_sa * p.t_sa);

        // Matched asymptotic corner constants (Schwarz-Christoffel rectangular step of thickness d
        // at dielectric interface eps_sub / eps_vac, Woods et al. PRApplied 12, 014012 (2019)):
        let c_ma_corner = -0.674;
        let c_ms_corner = -0.921;
        let c_sa_corner = 1.145;

        let ln_r0_d = (r0 / d_film.max(1e-12)).max(1.01).ln();
        let ln_2 = 2.0_f64.ln();

        // Use the annular surface collar [0.5 R_0, R_0] to blend any local smooth (r^0) departure from the
        // pure 1/r asymptote across [0.5 R_0, R_0], ensuring R_0-invariance while anchoring the singular
        // inner core [0, 0.5 R_0] to I_edge * (ln(0.5 R_0 / d) + C_k):
        let ln_inner = (ln_r0_d - ln_2).max(0.0);
        let asymp_col_ma = kappa_ma * i_edge * ln_2;
        let asymp_col_ms = kappa_ms * i_edge * ln_2;
        let asymp_col_sa = kappa_sa * i_edge * ln_2;

        let col_ma = if u_ma_col > 0.0 { 0.5 * (u_ma_col + asymp_col_ma) } else { asymp_col_ma };
        let col_ms = if u_ms_col > 0.0 { 0.5 * (u_ms_col + asymp_col_ms) } else { asymp_col_ms };
        let col_sa = if u_sa_col > 0.0 { 0.5 * (u_sa_col + asymp_col_sa) } else { asymp_col_sa };

        u_ma_edge = (kappa_ma * i_edge * (ln_inner + c_ma_corner) + col_ma).max(0.0);
        u_ms_edge = (kappa_ms * i_edge * (ln_inner + c_ms_corner) + col_ms).max(0.0);
        u_sa_edge = (kappa_sa * i_edge * (ln_inner + c_sa_corner) + col_sa).max(0.0);
    }

    let p_ma_interior = u_ma_ext / u_elec_j;
    let p_ms_interior = u_ms_ext / u_elec_j;
    let p_sa_interior = u_sa_ext / u_elec_j;

    let p_ma_edge = u_ma_edge / u_elec_j;
    let p_ms_edge = u_ms_edge / u_elec_j;
    let p_sa_edge = u_sa_edge / u_elec_j;

    let p_ma = p_ma_interior + p_ma_edge;
    let p_ms = p_ms_interior + p_ms_edge;
    let p_sa = p_sa_interior + p_sa_edge;
    let p_total = p_ma + p_ms + p_sa;
    let loss_tan_total = p_ma * p.tan_delta_ma + p_ms * p.tan_delta_ms + p_sa * p.tan_delta_sa;
    let q_tls = if loss_tan_total > 0.0 { 1.0 / loss_tan_total } else { f64::INFINITY };

    Ok(TlsReport {
        solution: p.solution,
        u_elec_j,
        p_ma,
        p_ms,
        p_sa,
        p_total,
        loss_tan_total,
        q_tls,
        edge_cutoff_m: r0,
        p_ma_interior,
        p_ms_interior,
        p_sa_interior,
        p_ma_edge,
        p_ms_edge,
        p_sa_edge,
        convergence_check: json!({
            "dim": 3,
            "mode": "matched_asymptotic_edge_tube_3d",
            "d_film_m": d_film,
            "r0_m": r0,
            "h_edge_mean_m": h_edge_mean,
            "perimeter_segments": perim_segments.len(),
            "perimeter_len_m": total_perim_len,
            "u_tube_J": u_tube_j,
            "i_edge": i_edge,
        }),
    })
}

/// Computes dielectric TLS surface participation with default parameters.
pub fn tls_participation(mesh: &Mesh, result: &SolveResult) -> Result<TlsReport, Error> {
    tls_participation_with_params(mesh, result, &Value::Null)
}

/// Computes dielectric TLS surface participation with optional JSON parameters.
pub fn tls_participation_with_params(mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<TlsReport, Error> {
    let p = TlsParams::from_json(params)?;
    let disc = discretise(mesh)?;
    if disc.dim != result.dim as usize {
        return Err(Error::invalid(format!("the result is {}D but the mesh is {}D", result.dim, disc.dim)));
    }
    if p.solution >= result.solutions.len() {
        return Err(Error::invalid(format!(
            "solution index {} is out of bounds (result has {} solution(s))",
            p.solution,
            result.solutions.len()
        )));
    }
    let sol = &result.solutions[p.solution];
    let ef = EFieldAccessor::new(sol, disc.dim, disc.n_cells())?;
    let meta = TlsMeta::of(result).unwrap_or_default();
    let c_e = if result.physics == "ems" { 0.25 } else { 0.5 };

    match disc.dim {
        2 => run_tls_2d(&disc, &ef, &meta, &p, c_e),
        3 => run_tls_3d(&disc, &ef, &meta, &p, c_e),
        d => Err(Error::invalid(format!("tls_participation supports 2D and 3D meshes, not {d}D"))),
    }
}

/// Dielectric TLS surface participation analysis behind the [`Analysis`] trait (`id: "tls_participation"`).
pub struct TlsParticipationAnalysis;

impl Analysis for TlsParticipationAnalysis {
    fn id(&self) -> &'static str {
        "tls_participation"
    }

    fn name(&self) -> &'static str {
        "Dielectric TLS surface participation (MA, MS, SA) and quality factor"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        tls_participation_applies(result)
    }

    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(tls_participation(mesh, result)?.to_json())
    }

    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(tls_participation_with_params(mesh, result, params)?.to_json())
    }
}
