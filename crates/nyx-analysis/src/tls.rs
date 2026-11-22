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
//! - **3D (`dim == 3`, Matched Asymptotic Edge-Tube Regularization + Flat AABB BVH)**: Detects the 1D singular
//!   perimeter `\partial\Gamma_{SC}` of interior zero-thickness superconducting sheets, indexes all perimeter
//!   segments in a deterministic flat 3D AABB bounding-volume hierarchy (`SegmentBvh3`) for `O(log N_seg)` distance
//!   queries, integrates exterior `MA`, `MS`, and `SA` surface facets outside distance `R_0` (`edge_cutoff`),
//!   extracts the convergent perimeter singularity intensity `I_edge = \oint |A(s)|^2 ds` from the 3D electric
//!   energy inside `r(x, \partial\Gamma_{SC}) < R_0`, and adds matched asymptotic finite-thickness (`d > 0`) and
//!   optional substrate-trenching (`trench_depth_m >= 0`) corner contributions.

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
    trench_depth_m: Option<f64>,
}

impl TlsParams {
    fn from_json(params: &Value, meta: &TlsMeta) -> Result<Self, Error> {
        let ifaces = meta.interfaces.as_ref();
        let def_t_ma = ifaces.and_then(|i| i.ma.as_ref()).map_or(2.0e-9, |s| s.thickness_m);
        let def_eps_ma = ifaces.and_then(|i| i.ma.as_ref()).map_or(10.0, |s| s.eps_r);
        let def_td_ma = ifaces.and_then(|i| i.ma.as_ref()).map_or(1.0e-3, |s| s.tan_delta);

        let def_t_ms = ifaces.and_then(|i| i.ms.as_ref()).map_or(2.0e-9, |s| s.thickness_m);
        let def_eps_ms = ifaces.and_then(|i| i.ms.as_ref()).map_or(11.45, |s| s.eps_r);
        let def_td_ms = ifaces.and_then(|i| i.ms.as_ref()).map_or(1.0e-3, |s| s.tan_delta);

        let def_t_sa = ifaces.and_then(|i| i.sa.as_ref()).map_or(2.0e-9, |s| s.thickness_m);
        let def_eps_sa = ifaces.and_then(|i| i.sa.as_ref()).map_or(10.0, |s| s.eps_r);
        let def_td_sa = ifaces.and_then(|i| i.sa.as_ref()).map_or(1.0e-3, |s| s.tan_delta);

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
        let trench_depth_m = match params.get("trench_depth_m") {
            None | Some(Value::Null) => None,
            Some(v) => Some(
                v.as_f64()
                    .filter(|x| x.is_finite() && *x >= 0.0)
                    .ok_or_else(|| Error::invalid("parameter 'trench_depth_m' must be a non-negative depth in m"))?,
            ),
        };
        Ok(Self {
            solution,
            t_ma: num("t_ma", def_t_ma)?,
            t_ms: num("t_ms", def_t_ms)?,
            t_sa: num("t_sa", def_t_sa)?,
            eps_ma: num("eps_ma", def_eps_ma)?.max(1e-12),
            eps_ms: num("eps_ms", def_eps_ms)?.max(1e-12),
            eps_sa: num("eps_sa", def_eps_sa)?.max(1e-12),
            tan_delta_ma: num("tan_delta_ma", def_td_ma)?,
            tan_delta_ms: num("tan_delta_ms", def_td_ms)?,
            tan_delta_sa: num("tan_delta_sa", def_td_sa)?,
            d,
            edge_cutoff,
            trench_depth_m,
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

pub(crate) struct EFieldAccessor<'a> {
    ex: &'a CellField,
    ey: Option<&'a CellField>,
    ez: &'a CellField,
    nloc: usize,
    pub(crate) order: usize,
}

impl<'a> EFieldAccessor<'a> {
    pub(crate) fn new(sol: &'a Solution, dim: usize, n_cells: usize) -> Result<Self, Error> {
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
    pub(crate) fn eval_with_shapes(&self, ci: usize, sh: &[f64]) -> [C64; 3] {
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

    #[inline]
    pub(crate) fn eval_real_with_shapes(&self, ci: usize, sh: &[f64]) -> [f64; 3] {
        let base = ci * self.nloc;
        let eval_one = |f: &CellField| -> f64 {
            let mut re = 0.0;
            for (a, &na) in sh.iter().enumerate().take(self.nloc) {
                re += na * f.re[base + a];
            }
            re
        };
        let vx = eval_one(self.ex);
        let vy = self.ey.map_or(0.0, eval_one);
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
                    u_ma += c_e * EPS0 * p.t_ma * ((eps_c * eps_c / p.eps_ma) * int_en2 + p.eps_ma * int_epar2);
                } else {
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

// ------------------------------------------------------------------ 3D Flat AABB BVH & Edge-Tube Helpers

#[inline]
pub(crate) fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

#[inline]
pub(crate) fn pt3(disc: &Discretisation, v: u32) -> [f64; 3] {
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
pub(crate) fn point_segment_dist_sq(x: [f64; 3], a: [f64; 3], b: [f64; 3], inv_len_sq: f64) -> f64 {
    let ab = sub3(b, a);
    let ax = sub3(x, a);
    let t = (dot3(ax, ab) * inv_len_sq).clamp(0.0, 1.0);
    let dx = ax[0] - t * ab[0];
    let dy = ax[1] - t * ab[1];
    let dz = ax[2] - t * ab[2];
    dx * dx + dy * dy + dz * dz
}

/// One 3D singular sheet-perimeter line segment.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Segment3 {
    pub(crate) a: [f64; 3],
    pub(crate) b: [f64; 3],
    pub(crate) inv_len_sq: f64,
}

#[derive(Clone, Copy, Debug)]
struct BvhNode3 {
    lo: [f64; 3],
    hi: [f64; 3],
    /// For internal nodes: left child index; for leaves (`count > 0`): start offset in `ordered_segs`.
    left_or_start: u32,
    /// For internal nodes (`count == 0`): right child index; for leaves: `0` unused.
    right: u32,
    /// `0` for internal nodes, `> 0` segment count for leaf nodes.
    count: u32,
}

#[inline]
fn aabb_dist_sq(lo: [f64; 3], hi: [f64; 3], x: [f64; 3]) -> f64 {
    let dx = if x[0] < lo[0] {
        lo[0] - x[0]
    } else if x[0] > hi[0] {
        x[0] - hi[0]
    } else {
        0.0
    };
    let dy = if x[1] < lo[1] {
        lo[1] - x[1]
    } else if x[1] > hi[1] {
        x[1] - hi[1]
    } else {
        0.0
    };
    let dz = if x[2] < lo[2] {
        lo[2] - x[2]
    } else if x[2] > hi[2] {
        x[2] - hi[2]
    } else {
        0.0
    };
    dx * dx + dy * dy + dz * dz
}

/// Deterministic flat `Vec`-backed 3D AABB bounding-volume hierarchy over 3D line segments (`Segment3`).
pub(crate) struct SegmentBvh3 {
    nodes: Vec<BvhNode3>,
    segs: Vec<Segment3>,
}

impl SegmentBvh3 {
    pub(crate) fn build(segs: Vec<Segment3>) -> Self {
        let n = segs.len();
        if n == 0 {
            return Self { nodes: Vec::new(), segs: Vec::new() };
        }
        struct Item {
            orig_idx: usize,
            lo: [f64; 3],
            hi: [f64; 3],
            cen: [f64; 3],
        }
        let mut items: Vec<Item> = segs
            .iter()
            .enumerate()
            .map(|(idx, s)| Item {
                orig_idx: idx,
                lo: [s.a[0].min(s.b[0]), s.a[1].min(s.b[1]), s.a[2].min(s.b[2])],
                hi: [s.a[0].max(s.b[0]), s.a[1].max(s.b[1]), s.a[2].max(s.b[2])],
                cen: [0.5 * (s.a[0] + s.b[0]), 0.5 * (s.a[1] + s.b[1]), 0.5 * (s.a[2] + s.b[2])],
            })
            .collect();

        let mut nodes = Vec::with_capacity(2 * n);
        let mut ordered_segs = Vec::with_capacity(n);

        fn build_rec(
            items: &mut [Item],
            orig_segs: &[Segment3],
            nodes: &mut Vec<BvhNode3>,
            ordered_segs: &mut Vec<Segment3>,
        ) -> u32 {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            let mut clo = [f64::INFINITY; 3];
            let mut chi = [f64::NEG_INFINITY; 3];
            for it in items.iter() {
                for k in 0..3 {
                    lo[k] = lo[k].min(it.lo[k]);
                    hi[k] = hi[k].max(it.hi[k]);
                    clo[k] = clo[k].min(it.cen[k]);
                    chi[k] = chi[k].max(it.cen[k]);
                }
            }
            let node_idx = nodes.len() as u32;
            if items.len() <= 4 {
                let start = ordered_segs.len() as u32;
                for it in items.iter() {
                    ordered_segs.push(orig_segs[it.orig_idx]);
                }
                nodes.push(BvhNode3 { lo, hi, left_or_start: start, right: 0, count: items.len() as u32 });
                return node_idx;
            }
            let ext = [chi[0] - clo[0], chi[1] - clo[1], chi[2] - clo[2]];
            let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                0
            } else if ext[1] >= ext[2] {
                1
            } else {
                2
            };
            items.sort_by(|a, b| a.cen[axis].total_cmp(&b.cen[axis]).then_with(|| a.orig_idx.cmp(&b.orig_idx)));
            let mid = items.len() / 2;
            nodes.push(BvhNode3 { lo, hi, left_or_start: 0, right: 0, count: 0 });
            let (left_slice, right_slice) = items.split_at_mut(mid);
            let left_idx = build_rec(left_slice, orig_segs, nodes, ordered_segs);
            let right_idx = build_rec(right_slice, orig_segs, nodes, ordered_segs);
            nodes[node_idx as usize].left_or_start = left_idx;
            nodes[node_idx as usize].right = right_idx;
            node_idx
        }

        build_rec(&mut items, &segs, &mut nodes, &mut ordered_segs);
        Self { nodes, segs: ordered_segs }
    }

    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.segs.is_empty()
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.segs.len()
    }

    /// Exact minimum squared Euclidean distance from `x` to all segments in `O(log N_seg)` with zero heap allocations.
    #[inline]
    pub(crate) fn min_dist_sq(&self, x: [f64; 3]) -> f64 {
        if self.nodes.is_empty() {
            return f64::INFINITY;
        }
        let mut best_d2 = f64::INFINITY;
        let mut stack = [0u32; 64];
        let mut sp = 1usize;
        stack[0] = 0;

        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp] as usize];
            if aabb_dist_sq(node.lo, node.hi, x) >= best_d2 {
                continue;
            }
            if node.count > 0 {
                let start = node.left_or_start as usize;
                let end = start + (node.count as usize);
                for s in &self.segs[start..end] {
                    let d2 = point_segment_dist_sq(x, s.a, s.b, s.inv_len_sq);
                    if d2 < best_d2 {
                        best_d2 = d2;
                    }
                }
            } else {
                let l_idx = node.left_or_start;
                let r_idx = node.right;
                let dl = aabb_dist_sq(self.nodes[l_idx as usize].lo, self.nodes[l_idx as usize].hi, x);
                let dr = aabb_dist_sq(self.nodes[r_idx as usize].lo, self.nodes[r_idx as usize].hi, x);
                // Push farther child first so nearer child is popped andTightens best_d2 first
                if dl <= dr {
                    if dr < best_d2 && sp < stack.len() {
                        stack[sp] = r_idx;
                        sp += 1;
                    }
                    if dl < best_d2 && sp < stack.len() {
                        stack[sp] = l_idx;
                        sp += 1;
                    }
                } else {
                    if dl < best_d2 && sp < stack.len() {
                        stack[sp] = l_idx;
                        sp += 1;
                    }
                    if dr < best_d2 && sp < stack.len() {
                        stack[sp] = r_idx;
                        sp += 1;
                    }
                }
            }
        }
        best_d2
    }

    #[inline]
    pub(crate) fn min_dist(&self, x: [f64; 3]) -> f64 {
        self.min_dist_sq(x).sqrt()
    }
}

pub(crate) type FaceToCells3D = HashMap<[u32; 3], Vec<(usize, [usize; 3])>>;

pub(crate) struct SheetPerimeter3D {
    pub(crate) face_to_cells: FaceToCells3D,
    pub(crate) sheet_tris: BTreeMap<[u32; 3], i32>,
    pub(crate) sa_tris: Vec<[u32; 3]>,
    pub(crate) bvh: SegmentBvh3,
    pub(crate) total_perim_len: f64,
    pub(crate) h_edge_mean: f64,
}

/// Detects interior zero-thickness superconducting sheet triangles, uncovered substrate-air triangles,
/// and the 1D singular sheet perimeter `\partial\Gamma_{SC}` indexed in a [`SegmentBvh3`].
pub(crate) fn detect_sheet_perimeter_3d(disc: &Discretisation, meta: &TlsMeta, cell_eps: &[f64]) -> SheetPerimeter3D {
    let nc = disc.n_cells();
    let mut face_to_cells: FaceToCells3D = HashMap::with_capacity(nc * 2);
    for ci in 0..nc {
        let v = disc.cell(ci);
        for lf in FACES_3D {
            let tri_raw = [v[lf[0]], v[lf[1]], v[lf[2]]];
            let key = sorted3(tri_raw[0], tri_raw[1], tri_raw[2]);
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
    let bvh = SegmentBvh3::build(perim_segments);

    SheetPerimeter3D { face_to_cells, sheet_tris, sa_tris, bvh, total_perim_len, h_edge_mean }
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

/// Generates `(2 * n_sub)^3` quadrature points on the 3D reference tetrahedron by collapsing the subdivided cube,
/// returning `(bary4, weight)` where `sum(weight) == 1.0`.
pub(crate) fn subdivided_tet_quad(n_sub: usize) -> Vec<([f64; 4], f64)> {
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

    // 2. Detect sheet triangles, SA triangles, and build flat 3D AABB BVH over singular perimeter
    let perim = detect_sheet_perimeter_3d(disc, meta, &cell_eps);
    let d_film = p.d.or_else(|| meta.sheet_thickness_m.values().next().copied()).unwrap_or(100.0e-9);
    let r0 = p.edge_cutoff.unwrap_or_else(|| (1.15 * perim.h_edge_mean).max(4.0 * d_film));
    let r0_sq = r0 * r0;
    let half_r0_sq = (0.5 * r0) * (0.5 * r0);

    // 3. Exterior surface integrals (r(x, \partial\Gamma_{SC}) > R_0)
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

        let min_v_dist = if perim.bvh.is_empty() {
            f64::INFINITY
        } else {
            let d0 = perim.bvh.min_dist(p0);
            let d1 = perim.bvh.min_dist(p1);
            let d2 = perim.bvh.min_dist(p2);
            d0.min(d1).min(d2)
        };

        let mut ext_en2 = 0.0;
        let mut ext_epar2 = 0.0;
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
                let dist_sq = perim.bvh.min_dist_sq(x);
                if dist_sq <= half_r0_sq {
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
                if dist_sq > r0_sq {
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

    for &tri in perim.sheet_tris.keys() {
        let Some(cells) = perim.face_to_cells.get(&tri) else { continue };
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

    for &tri in &perim.sa_tris {
        let Some(cells) = perim.face_to_cells.get(&tri) else { continue };
        if let Some(&(c_vac, lf)) = cells.iter().find(|&&(ci, _)| cell_eps[ci] <= 1.5) {
            let eps_c = cell_eps[c_vac];
            let (ext_en2, ext_epar2, col_en2, col_epar2) = integrate_surface_tri(tri, c_vac, lf);
            u_sa_ext += c_e * EPS0 * p.t_sa * (p.eps_sa * ext_epar2 + (eps_c * eps_c / p.eps_sa) * ext_en2);
            u_sa_col += c_e * EPS0 * p.t_sa * (p.eps_sa * col_epar2 + (eps_c * eps_c / p.eps_sa) * col_en2);
        }
    }

    // 4. 3D Edge-Tube Volume Integral & Matched Asymptotic Regularization (r <= R_0)
    let (mut u_ma_edge, mut u_ms_edge, mut u_sa_edge) = (0.0, 0.0, 0.0);
    let mut i_edge = 0.0;
    let mut u_tube_j = 0.0;
    let trench_m = p.trench_depth_m.unwrap_or(0.0);

    if !perim.bvh.is_empty() && r0 > 0.0 {
        let sub_tet = subdivided_tet_quad(4);
        let sub_tet_shapes: Vec<Vec<f64>> = sub_tet.iter().map(|(l, _)| shape(l, ef.order).0).collect();

        for ci in 0..nc {
            let v = disc.cell(ci);
            let pts = [pt3(disc, v[0]), pt3(disc, v[1]), pt3(disc, v[2]), pt3(disc, v[3])];
            let mut d_min = f64::INFINITY;
            let mut d_max = 0.0_f64;
            for &pt in &pts {
                let d = perim.bvh.min_dist(pt);
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
                    if perim.bvh.min_dist_sq(x) <= r0_sq {
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
        // with optional substrate gap trenching depth h_tr >= 0, Woods et al. PRApplied 12, 014012 (2019)):
        let (c_ma_corner, c_ms_corner, c_sa_corner) = if trench_m > 0.0 {
            let tau_ln = (1.0 + trench_m / d_film.max(1e-12)).ln();
            (-0.674 + 0.18 * tau_ln, -0.921 - 0.55 * tau_ln, 1.145 - 0.65 * tau_ln)
        } else {
            (-0.674, -0.921, 1.145)
        };

        let ln_r0_d = (r0 / d_film.max(1e-12)).max(1.01).ln();
        let ln_2 = 2.0_f64.ln();

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

    let mut conv = json!({
        "dim": 3,
        "mode": "matched_asymptotic_edge_tube_3d",
        "d_film_m": d_film,
        "r0_m": r0,
        "h_edge_mean_m": perim.h_edge_mean,
        "perimeter_segments": perim.bvh.len(),
        "perimeter_len_m": perim.total_perim_len,
        "u_tube_J": u_tube_j,
        "i_edge": i_edge,
    });
    if let Some(htr) = p.trench_depth_m {
        conv["trench_depth_m"] = json!(htr);
    }

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
        convergence_check: conv,
    })
}

/// Computes dielectric TLS surface participation with default parameters.
pub fn tls_participation(mesh: &Mesh, result: &SolveResult) -> Result<TlsReport, Error> {
    tls_participation_with_params(mesh, result, &Value::Null)
}

/// Computes dielectric TLS surface participation with optional JSON parameters.
pub fn tls_participation_with_params(mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<TlsReport, Error> {
    let meta = TlsMeta::of(result).unwrap_or_default();
    let p = TlsParams::from_json(params, &meta)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use nyx_physics::{TlsInterfaceSpec, TlsInterfaces};
    use std::time::Instant;

    #[test]
    fn segment_bvh3_matches_exact_brute_force_and_queries_4000_segments_fast() {
        // Build a 4,000-segment meandering 3D perimeter
        let n_seg = 4_000usize;
        let mut segs = Vec::with_capacity(n_seg);
        for i in 0..n_seg {
            let t0 = (i as f64) * 0.01;
            let t1 = ((i + 1) as f64) * 0.01;
            let a = [100e-6 * t0.cos() + 2e-6 * (7.0 * t0).sin(), 100e-6 * t0.sin() + 1e-6 * i as f64, 0.0];
            let b = [100e-6 * t1.cos() + 2e-6 * (7.0 * t1).sin(), 100e-6 * t1.sin() + 1e-6 * (i + 1) as f64, 0.0];
            let ab = sub3(b, a);
            let l2 = dot3(ab, ab);
            segs.push(Segment3 { a, b, inv_len_sq: 1.0 / l2 });
        }

        let t0 = Instant::now();
        let bvh = SegmentBvh3::build(segs.clone());
        let mut checksum = 0.0;
        for q in 0..20_000 {
            let fq = q as f64;
            let pt = [95e-6 * (0.003 * fq).cos(), 0.2 * fq * 1e-6, 1.5e-6 * (0.01 * fq).sin()];
            checksum += bvh.min_dist_sq(pt);
        }
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1e3;
        assert!(checksum > 0.0);
        assert!(elapsed_ms < 500.0, "4000-seg BVH build + 20k queries took {elapsed_ms:.2} ms >= 500 ms");

        // Verify bit-for-bit exactness vs brute force on 100 sample points
        for q in 0..100 {
            let fq = (q * 137) as f64;
            let pt = [80e-6 * (0.11 * fq).cos(), 25e-6 * fq * 0.05, 3e-6 * (0.07 * fq).sin()];
            let got = bvh.min_dist_sq(pt);
            let want =
                segs.iter().map(|s| point_segment_dist_sq(pt, s.a, s.b, s.inv_len_sq)).fold(f64::INFINITY, f64::min);
            assert_eq!(got.to_bits(), want.to_bits(), "BVH mismatch at q={q}: got={got:e}, want={want:e}");
        }
    }

    #[test]
    fn tls_params_inherit_meta_interfaces_and_allow_explicit_overrides() {
        let meta = TlsMeta {
            interfaces: Some(TlsInterfaces {
                ma: Some(TlsInterfaceSpec { thickness_m: 3.2e-9, eps_r: 9.5, tan_delta: 2.1e-3 }),
                ms: Some(TlsInterfaceSpec { thickness_m: 1.8e-9, eps_r: 9.8, tan_delta: 1.5e-3 }),
                sa: Some(TlsInterfaceSpec { thickness_m: 2.4e-9, eps_r: 3.9, tan_delta: 1.2e-3 }),
            }),
            ..Default::default()
        };
        let p_def = TlsParams::from_json(&Value::Null, &meta).unwrap();
        assert_eq!(p_def.t_ma, 3.2e-9);
        assert_eq!(p_def.eps_ma, 9.5);
        assert_eq!(p_def.tan_delta_ma, 2.1e-3);
        assert_eq!(p_def.t_ms, 1.8e-9);
        assert_eq!(p_def.eps_ms, 9.8);
        assert_eq!(p_def.tan_delta_ms, 1.5e-3);
        assert_eq!(p_def.t_sa, 2.4e-9);
        assert_eq!(p_def.eps_sa, 3.9);
        assert_eq!(p_def.tan_delta_sa, 1.2e-3);

        // Explicit JSON param overrides TOML meta default while leaving others from meta
        let p_ov = TlsParams::from_json(&json!({ "tan_delta_ms": 5.0e-4, "trench_depth_m": 150.0e-9 }), &meta).unwrap();
        assert_eq!(p_ov.tan_delta_ms, 5.0e-4);
        assert_eq!(p_ov.t_ma, 3.2e-9);
        assert_eq!(p_ov.trench_depth_m, Some(150.0e-9));

        // Negative trench_depth_m is rejected
        assert!(TlsParams::from_json(&json!({ "trench_depth_m": -1.0e-9 }), &meta).is_err());
    }
}
