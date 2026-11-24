// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Barsoum quarter-point (`xi_QP ~ 0.25`) `CurvedTet` sheet-edge singularity enrichment and
//! Morrison-Lewis / Babuska `r^alpha` radical corner grading (`AF-2`).
//!
//! Around any 1D perimeter edge (`"<sheet>:edge"`) of a zero-thickness conductive sheet, the electrostatic
//! potential expands in cylindrical coordinates as `V(r, theta) = V_0 + A_1 r^(1/2) sin(theta/2) + O(r)`,
//! while flat 2D polygonal corners of interior sheet opening angle `omega in (0, 2*pi)` induce 3D vertex
//! singularities `V(rho) ~ rho^alpha` with the Morrison-Lewis corner exponent `alpha(omega)`.
//!
//! This module provides:
//! 1. Piecewise-interpolated Morrison-Lewis corner singularity exponents `alpha(omega, eps_1 / eps_2)` and
//!    helpers measuring `omega` both by summing incident 2D sheet-covered triangle angles on `m2` and
//!    directly from `st.overlay` sector tangents;
//! 2. Babuska's optimal corner grading law `h(d) = max(h_min, h_edge * (d / R_c)^(1 - alpha / (p + 0.5)))`
//!    (`p = 2`, `R_c = 2 * h_edge`, floor `h_min = 0.05 * h_edge`) with coupled vertical `dz` step sizing;
//! 3. Pure `femcore` 3D quarter-point `EdgeNodes` enrichment [`apply_quarter_point_nodes`], shifting only
//!    edges `(a, b)` with `rim_verts[a] ^ rim_verts[b]`, clamping `xi_QP` upward on curved rims with positive
//!    sagitta `delta > 0` (`4 * xi_QP - 1 > 2 * delta / h_y + 0.05`), and validating strict interior
//!    Jacobian positivity (`det(J) * det_straight > 0`) at all 216 Gauss points of `quad_rule_collapsed(8)`.

use crate::mesh::SimplexMesh;
use crate::meshing::curve3::scaled_jacobian_min_interior_q8;
use geom::complex::{Complex, FaceGeom};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{PI, TAU};

/// Default Barsoum quarter-point parameter (`0.25` from the singular rim vertex).
pub const DEFAULT_QUARTER_POINT_XI: f64 = 0.25;

/// Polynomial degree `p = 2` for Babuska corner grading (`p + 0.5 = 2.5`).
pub const CORNER_GRADING_P: f64 = 2.0;

/// Corner grading influence radius `R_c / h_edge = 2.0`.
pub const CORNER_RADIUS_FACTOR: f64 = 2.0;

/// Mandatory corner grading floor `h_min / h_edge = 0.05` (prevents runaway subdivision).
pub const CORNER_H_MIN_FACTOR: f64 = 0.05;

/// Tabulated Morrison-Lewis (1976) 3D Laplace corner singularity exponents `(omega_rad, alpha)` for a flat
/// conducting sector of interior angle `omega in [0, 2*pi]`.
///
/// Key reference points:
/// - `omega = pi / 2` (`90 deg` convex corner): `alpha = 0.2966` (`~0.297`)
/// - `omega = pi` (`180 deg` straight smooth edge): `alpha = 0.5000` (exact `r^(1/2)` half-plane edge)
/// - `omega = 3 * pi / 2` (`270 deg` re-entrant notch): `alpha = 0.8147` (`~0.815`)
/// - `omega -> 2 * pi` (`360 deg` narrow slit): `alpha -> 1.0000`
pub const MORRISON_LEWIS_TABLE: [(f64, f64); 18] = [
    (0.0, 0.1000),
    (PI / 12.0, 0.1602),
    (PI / 6.0, 0.2142),
    (PI / 4.0, 0.2496),
    (PI / 3.0, 0.2703),
    (PI / 2.0, 0.2966),
    (2.0 * PI / 3.0, 0.3554),
    (3.0 * PI / 4.0, 0.3854),
    (5.0 * PI / 6.0, 0.4192),
    (PI, 0.5000),
    (7.0 * PI / 6.0, 0.5905),
    (5.0 * PI / 4.0, 0.6410),
    (4.0 * PI / 3.0, 0.6955),
    (3.0 * PI / 2.0, 0.8147),
    (5.0 * PI / 3.0, 0.8870),
    (7.0 * PI / 4.0, 0.9195),
    (11.0 * PI / 6.0, 0.9500),
    (TAU, 1.0000),
];

/// Piecewise-linearly interpolated Morrison-Lewis 3D corner singularity exponent `alpha(omega, eps_1 / eps_2)`
/// for a flat conducting sheet corner of interior conductor angle `omega_rad in (0, 2*pi)` at a dielectric
/// interface with permittivity ratio `eps_ratio = eps_1 / eps_2`.
///
/// For a zero-thickness conducting sheet lying on a planar interface `z = const` between two half-spaces
/// `eps_1` and `eps_2`, even reflection symmetry in `z` makes `dV/dz = 0` on the dielectric complement
/// outside the sheet, so the Laplace corner eigenvalue is invariant under `eps_1 / eps_2 > 0`.
pub fn morrison_lewis_exponent(omega_rad: f64, eps_ratio: f64) -> f64 {
    let _ = if eps_ratio.is_finite() && eps_ratio > 0.0 { eps_ratio } else { 1.0 };
    if !omega_rad.is_finite() {
        return 0.5;
    }
    let w = omega_rad.clamp(0.0, TAU);
    if w <= MORRISON_LEWIS_TABLE[0].0 {
        return MORRISON_LEWIS_TABLE[0].1;
    }
    for win in MORRISON_LEWIS_TABLE.windows(2) {
        let (w0, a0) = win[0];
        let (w1, a1) = win[1];
        if w <= w1 {
            let t = if (w1 - w0).abs() > 1e-15 { (w - w0) / (w1 - w0) } else { 0.0 };
            return a0 + t * (a1 - a0);
        }
    }
    MORRISON_LEWIS_TABLE.last().unwrap().1
}

/// Alias for [`morrison_lewis_exponent`].
#[inline]
pub fn corner_singularity_exponent(omega_rad: f64, eps_ratio: f64) -> f64 {
    morrison_lewis_exponent(omega_rad, eps_ratio)
}

/// Babuska's optimal radical corner grading law:
/// `h(d) = max(h_min, h_edge * (d / R_c)^(1 - alpha / (p + 0.5)))`
/// with `p = 2`, `R_c = 2 * h_edge`, and mandatory floor `h_min = 0.05 * h_edge`.
/// For `d >= R_c`, returns `h_edge`.
pub fn corner_grading_h(d: f64, h_edge: f64, alpha: f64) -> f64 {
    if !(h_edge > 0.0 && h_edge.is_finite()) {
        return 0.0;
    }
    let h_min = CORNER_H_MIN_FACTOR * h_edge;
    let r_c = CORNER_RADIUS_FACTOR * h_edge;
    if !(d > 0.0 && d.is_finite()) {
        return h_min;
    }
    if d >= r_c {
        return h_edge;
    }
    let alpha_clamped = if alpha.is_finite() { alpha.clamp(0.05, 1.0) } else { 0.5 };
    let exponent = (1.0 - alpha_clamped / (CORNER_GRADING_P + 0.5)).clamp(0.0, 1.0);
    let raw = h_edge * (d / r_c).powf(exponent);
    raw.max(h_min).min(h_edge)
}

/// Alias for [`corner_grading_h`].
#[inline]
pub fn babuska_corner_h(d: f64, h_edge: f64, alpha: f64) -> f64 {
    corner_grading_h(d, h_edge, alpha)
}

/// A detected 2D conductive sheet corner on `st.overlay` or `m2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SheetCorner {
    pub vertex: usize,
    pub pos: [f64; 2],
    /// Interior sheet angle `omega in (0, 2*pi)` (radians).
    pub omega_rad: f64,
    /// Morrison-Lewis singularity exponent `alpha(omega)`.
    pub alpha: f64,
    /// Baseline target edge size `h_edge` at this corner.
    pub h_edge: f64,
}

impl SheetCorner {
    /// True when the interior angle deviates genuinely from a smooth `180 deg` straight or filleted edge
    /// (`|omega - pi| > 0.35 rad`, ~20 deg — preserving 30/45/60/90/135/270 deg physical corners while
    /// ignoring dense polyline chords along smooth circular arcs and corner fillets).
    pub fn is_sharp(&self) -> bool {
        (self.omega_rad - PI).abs() > 0.35 && self.omega_rad > 1e-4 && self.omega_rad < TAU - 1e-4
    }

    /// Graded element size `h(d)` at point `p` induced by this corner.
    pub fn h_at(&self, p: [f64; 2]) -> f64 {
        let d = (p[0] - self.pos[0]).hypot(p[1] - self.pos[1]);
        corner_grading_h(d, self.h_edge, self.alpha)
    }

    /// Minimum floor size `h_min = 0.05 * h_edge` at the corner apex (`d = 0`).
    pub fn h_min(&self) -> f64 {
        CORNER_H_MIN_FACTOR * self.h_edge
    }
}

/// Evaluates the corner-graded 2D element size at point `p` across all sharp corners in `corners`,
/// clamped above by `base_h`.
pub fn corner_graded_h_at(p: [f64; 2], corners: &[SheetCorner], base_h: f64) -> f64 {
    let mut h = base_h;
    for c in corners {
        if !c.is_sharp() {
            continue;
        }
        let r_c = CORNER_RADIUS_FACTOR * c.h_edge;
        let dx = (p[0] - c.pos[0]).abs();
        let dy = (p[1] - c.pos[1]).abs();
        if dx <= r_c && dy <= r_c {
            let d = dx.hypot(dy);
            if d < r_c {
                h = h.min(corner_grading_h(d, c.h_edge, c.alpha));
            }
        }
    }
    h
}

/// Coupled vertical slab step `dz` at a sheet plane when `corner_grading` is active:
/// clamps `base_dz` down to the minimum corner floor `h_min = 0.05 * h_edge` across all active
/// sharp corners on that sheet plane so corner prisms do not become high-aspect-ratio vertical needles.
pub fn coupled_corner_dz(base_dz: f64, corners: &[SheetCorner]) -> f64 {
    let mut dz = base_dz;
    for c in corners {
        if c.is_sharp() && c.h_edge > 0.0 && c.h_edge.is_finite() {
            dz = dz.min(c.h_min());
        }
    }
    dz
}

/// Measures the interior sheet opening angle `omega in (0, 2*pi)` at every vertex of a 2D overlay triangulation `m2`
/// by summing the interior triangle angles of all incident 2D triangles whose 0-based overlay cell index
/// `(m2.cell_labels[ti] - 1) as usize` is marked `true` in `sheet_covered_cells`.
///
/// Returns `Some(omega)` for vertices lying on the 1D sheet boundary (`0 < omega < 2*pi`), and `None` for
/// vertices strictly outside the sheet or strictly in the smooth 2D interior of the sheet (`omega == 2*pi`).
pub fn measure_sheet_corner_angles_m2(m2: &SimplexMesh, sheet_covered_cells: &[bool]) -> Vec<Option<f64>> {
    if m2.dim != 2 {
        return vec![];
    }
    let nv = m2.n_points();
    let mut angle_sum = vec![0.0f64; nv];
    let mut has_sheet = vec![false; nv];
    let mut has_non_sheet = vec![false; nv];
    let mut edge_sheet_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();

    for ti in 0..m2.n_cells() {
        let v = m2.cell(ti);
        let ov_cell = (m2.cell_labels[ti] - 1).max(0) as usize;
        let is_sheet = sheet_covered_cells.get(ov_cell).copied().unwrap_or(false);
        if !is_sheet {
            for &vi in v {
                has_non_sheet[vi as usize] = true;
            }
            continue;
        }
        for &(i, j) in &[(0, 1), (1, 2), (2, 0)] {
            let key = (v[i].min(v[j]), v[i].max(v[j]));
            *edge_sheet_count.entry(key).or_insert(0) += 1;
        }
        let p = [
            [m2.points[2 * v[0] as usize], m2.points[2 * v[0] as usize + 1]],
            [m2.points[2 * v[1] as usize], m2.points[2 * v[1] as usize + 1]],
            [m2.points[2 * v[2] as usize], m2.points[2 * v[2] as usize + 1]],
        ];
        for k in 0..3 {
            let a = p[k];
            let b = p[(k + 1) % 3];
            let c = p[(k + 2) % 3];
            let ab = [b[0] - a[0], b[1] - a[1]];
            let ac = [c[0] - a[0], c[1] - a[1]];
            let lb = ab[0].hypot(ab[1]);
            let lc = ac[0].hypot(ac[1]);
            if lb > 1e-30 && lc > 1e-30 {
                let cos_t = ((ab[0] * ac[0] + ab[1] * ac[1]) / (lb * lc)).clamp(-1.0, 1.0);
                angle_sum[v[k] as usize] += cos_t.acos();
                has_sheet[v[k] as usize] = true;
            }
        }
    }

    let mut on_rim = vec![false; nv];
    for (&(a, b), &cnt) in &edge_sheet_count {
        if cnt == 1 {
            on_rim[a as usize] = true;
            on_rim[b as usize] = true;
        }
    }

    (0..nv)
        .map(|v| {
            if has_sheet[v] && (on_rim[v] || has_non_sheet[v] || angle_sum[v] < TAU - 1e-5) && angle_sum[v] > 1e-6 {
                Some(angle_sum[v].min(TAU))
            } else {
                None
            }
        })
        .collect()
}

/// Measures the interior sheet opening angle `omega in (0, 2*pi)` at every `st.overlay` vertex `0..ov.vertices.len()`
/// on sheet plane `z_s` directly from the 2D planar arrangement `st.overlay` before triangulation.
pub fn measure_overlay_sheet_corners(cx: &Complex, z_s: f64) -> Vec<Option<f64>> {
    let Some(st) = cx.stack.as_ref() else {
        return vec![];
    };
    let ov = &st.overlay;
    let tol = cx.tol;
    let ncell = ov.regions.len();
    let mut cell_sheet = vec![false; ncell];
    for f in &cx.faces {
        if f.is_sheet() {
            if let FaceGeom::Cap { cell, z } = f.geom {
                if (z - z_s).abs() <= tol && (cell as usize) < ncell {
                    cell_sheet[cell as usize] = true;
                }
            }
        }
    }
    let nv = ov.vertices.len();
    type OutgoingRay = (f64, Option<u32>, Option<u32>);
    // For each vertex v, collect outgoing edge rays: (angle_in_-pi_pi, ccw_left_cell, cw_right_cell)
    let mut rays: Vec<Vec<OutgoingRay>> = vec![Vec::new(); nv];
    for ed in &ov.edges {
        if ed.v0 == ed.v1 {
            continue;
        }
        let curve = &ov.curves[ed.curve as usize];
        let t0_tan = curve.tangent(ed.t0);
        let a0 = t0_tan[1].atan2(t0_tan[0]);
        rays[ed.v0 as usize].push((a0, ed.left, ed.right));

        let t1_tan = curve.tangent(ed.t1);
        let out1 = [-t1_tan[0], -t1_tan[1]];
        let a1 = out1[1].atan2(out1[0]);
        rays[ed.v1 as usize].push((a1, ed.right, ed.left));
    }

    (0..nv)
        .map(|v| {
            let rv = &mut rays[v];
            if rv.len() < 2 {
                return None;
            }
            rv.sort_by(|a, b| a.0.total_cmp(&b.0));
            let k = rv.len();
            let mut omega = 0.0f64;
            let mut has_sheet_sector = false;
            let mut has_non_sheet_sector = false;
            for m in 0..k {
                let next = (m + 1) % k;
                let dtheta = if m + 1 < k { rv[next].0 - rv[m].0 } else { rv[0].0 + TAU - rv[m].0 };
                let cell_opt = rv[m].1.or(rv[next].2);
                let is_sheet = cell_opt.is_some_and(|c| cell_sheet.get(c as usize).copied().unwrap_or(false));
                if is_sheet {
                    omega += dtheta;
                    has_sheet_sector = true;
                } else {
                    has_non_sheet_sector = true;
                }
            }
            if has_sheet_sector && has_non_sheet_sector && omega > 1e-6 && omega < TAU - 1e-6 {
                Some(omega)
            } else {
                None
            }
        })
        .collect()
}

/// Collects all active sharp [`SheetCorner`]s on sheet plane `z_s` of `cx`, pairing each corner with its
/// minimum incident overlay edge target size `ov_edge_h[e]` and Morrison-Lewis exponent `alpha(omega, eps_ratio)`.
pub fn collect_overlay_sheet_corners(
    cx: &Complex,
    z_s: f64,
    ov_edge_h: &[f64],
    eps_ratio: f64,
) -> Vec<SheetCorner> {
    let Some(st) = cx.stack.as_ref() else {
        return vec![];
    };
    let ov = &st.overlay;
    let angles = measure_overlay_sheet_corners(cx, z_s);
    let mut corners = Vec::new();
    for (v, omega_opt) in angles.into_iter().enumerate() {
        let Some(omega_rad) = omega_opt else { continue };
        let mut h_edge = f64::INFINITY;
        for (e, ed) in ov.edges.iter().enumerate() {
            if ed.v0 as usize == v || ed.v1 as usize == v {
                if let Some(&eh) = ov_edge_h.get(e) {
                    h_edge = h_edge.min(eh);
                }
            }
        }
        if !h_edge.is_finite() || h_edge <= 0.0 {
            continue;
        }
        let alpha = morrison_lewis_exponent(omega_rad, eps_ratio);
        let p = ov.vertices[v].p;
        let sc = SheetCorner { vertex: v, pos: [p[0], p[1]], omega_rad, alpha, h_edge };
        if sc.is_sharp() {
            corners.push(sc);
        }
    }
    corners
}

/// Exact analytical 3D Jacobian determinant ratio `det(J(lam)) / det(J_straight)` for a tetrahedron whose
/// edges `(v_i, v_j)` with `is_rim[i] ^ is_rim[j]` are shifted to `(1 - xi_qp) * x_rim + xi_qp * x_other`:
///
/// With `alpha = 2 * (1 - 2 * xi_qp)`, `g(lam_out) = (4 * xi_qp - 1) + alpha * lam_out`, and
/// `lam_out = sum_{j notin S_rim} lam_j`:
/// - **Case 1 (`n_rim_in_tet == 2`, edge-adjacent tet)**:
///   `(1 + alpha * lam_out) * g(lam_out) * (g(lam_out) + alpha * lam_out)` (`-> 2 * lam_out^2 * (1 + lam_out)` at `xi_qp = 0.25`)
/// - **Case 2 (`n_rim_in_tet == 1`, vertex-adjacent tet)**:
///   `g(lam_out)^2 * (g(lam_out) + alpha * lam_out)` (`-> 2 * lam_out^3` at `xi_qp = 0.25`)
/// - **Case 3 (`n_rim_in_tet == 3`, corner/gap tet, 3 vertices on rim)**:
///   `(1 + alpha * lam_out)^2 * (g(lam_out) + alpha * lam_out)` (`-> 2 * lam_out * (1 + lam_out)^2` at `xi_qp = 0.25`)
/// - **Unshifted (`n_rim_in_tet == 0` or `4`)**: `1.0`.
pub fn analytical_qp_jacobian_ratio(n_rim_in_tet: usize, lam_out: f64, xi_qp: f64) -> f64 {
    let alpha = 2.0 * (1.0 - 2.0 * xi_qp);
    let g = (4.0 * xi_qp - 1.0) + alpha * lam_out;
    match n_rim_in_tet {
        1 => g * g * (g + alpha * lam_out),
        2 => (1.0 + alpha * lam_out) * g * (g + alpha * lam_out),
        3 => (1.0 + alpha * lam_out).powi(2) * (g + alpha * lam_out),
        _ => 1.0,
    }
}

/// Identifies all 1D open perimeter vertices (`rim_verts[v] == true`) of interior conductive sheet facets
/// carrying any boundary label in `sheet_labels` (edges belonging to exactly 1 interior facet in `sheet_labels`).
pub fn sheet_rim_vertices(mesh: &SimplexMesh, sheet_labels: &[i32]) -> Vec<bool> {
    let nv = mesh.n_points();
    if mesh.dim != 3 || sheet_labels.is_empty() {
        return vec![false; nv];
    }
    let label_set: BTreeSet<i32> = sheet_labels.iter().copied().collect();
    let exterior = mesh.exterior_facets();
    let mut edge_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for f in 0..mesh.n_facets() {
        if !exterior[f] && label_set.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            for &(i, j) in &[(0, 1), (1, 2), (2, 0)] {
                let key = (fv[i].min(fv[j]), fv[i].max(fv[j]));
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
    }
    let mut is_rim = vec![false; nv];
    for (&(a, b), &cnt) in &edge_count {
        if cnt == 1 {
            is_rim[a as usize] = true;
            is_rim[b as usize] = true;
        }
    }
    is_rim
}

/// Summary report returned by [`apply_quarter_point_nodes`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuarterPointReport {
    /// Number of vertices marked `true` in `rim_verts`.
    pub rim_vertices: usize,
    /// Number of mesh edges `(a, b)` with `rim_verts[a] ^ rim_verts[b]` that carry a quarter-point node.
    pub shifted_edges: usize,
    /// Number of tetrahedra touching at least one shifted quarter-point edge.
    pub qp_cells: usize,
    /// Number of quarter-point edges touching a curved rim segment whose `xi_QP` was clamped upward (`4*xi - 1 > 2*delta/h_y + 0.05`).
    pub clamped_curved_edges: usize,
    /// Number of candidate quarter-point edges skipped or unshifted to preserve strict interior Jacobian positivity.
    pub unshifted_edges: usize,
    /// Smallest `det(J) / det(J_straight)` over all 216 Gauss points of `quad_rule_collapsed(8)` across all `qp_cells`.
    pub min_interior_jacobian_ratio: f64,
}

impl Default for QuarterPointReport {
    fn default() -> Self {
        QuarterPointReport {
            rim_vertices: 0,
            shifted_edges: 0,
            qp_cells: 0,
            clamped_curved_edges: 0,
            unshifted_edges: 0,
            min_interior_jacobian_ratio: 1.0,
        }
    }
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Applies 3D Barsoum quarter-point mid-edge geometry nodes (`EdgeNodes`) to `mesh`:
///
/// 1. Shifts **only** edges `(a, b)` with `rim_verts[a] ^ rim_verts[b]` to `(1.0 - xi_eff) * p_rim + xi_eff * p_other`
///    (skipping all edges where both endpoints lie on `rim_verts`, such as rim segments or gap-spanning chords).
/// 2. For edges touching curved rim edges (`curved_rim_edges` or existing curved rim `edge_nodes` with `sagitta > 0`),
///    clamps `xi_eff` upward so `4 * xi_eff - 1 > 2 * delta / h_y + 0.05` (or skips edges where `xi_min >= 0.5`).
/// 3. Validates strict positivity `det(J) * det_straight > 0.0` at all 216 interior Gauss points of
///    `quad_rule_collapsed(8)` ([`scaled_jacobian_min_interior_q8`]), unshifting any edge whose cell fails.
/// 4. Records active quarter-point edges via [`SimplexMesh::set_qp_edges`] so [`crate::meshing::request::MeshReport`]
///    validates quarter-point cells on interior `quad_rule_collapsed(8)` positivity rather than vertex `det J = 0`.
pub fn apply_quarter_point_nodes(
    mesh: &mut SimplexMesh,
    rim_verts: &[bool],
    xi_qp: f64,
    curved_rim_edges: &BTreeSet<[u32; 2]>,
) -> QuarterPointReport {
    let mut report = QuarterPointReport::default();
    if mesh.dim != 3 || rim_verts.len() != mesh.n_points() {
        return report;
    }
    let xi_base = if xi_qp.is_finite() { xi_qp.clamp(0.25, 0.5) } else { DEFAULT_QUARTER_POINT_XI };
    report.rim_vertices = rim_verts.iter().filter(|&&r| r).count();
    if report.rim_vertices == 0 || (xi_base - 0.5).abs() <= 1e-12 {
        mesh.set_qp_edges(&[]);
        return report;
    }

    let nv = mesh.n_points();
    let pt = |v: usize| -> [f64; 3] {
        let p = mesh.point(v);
        [p[0], p[1], p[2]]
    };

    // Existing curved edge nodes (e.g. from AF-5 circular arc projection, which runs before quarter-point shifting)
    let mut edge_map: BTreeMap<[u32; 2], [f64; 3]> = mesh.edge_node_list().into_iter().collect();

    // Detect curved rim segments and their maximum sagitta `delta` at each rim vertex
    let mut touches_curved_rim = vec![false; nv];
    let mut max_rim_sagitta = vec![0.0f64; nv];

    let mut record_curved_rim_edge = |a: u32, b: u32, q_opt: Option<[f64; 3]>| {
        let (ua, ub) = (a as usize, b as usize);
        if ua >= nv || ub >= nv {
            return;
        }
        let (pa, pb) = (pt(ua), pt(ub));
        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
        let delta = q_opt.map_or(0.0, |q| dist3(q, mid));
        touches_curved_rim[ua] = true;
        touches_curved_rim[ub] = true;
        max_rim_sagitta[ua] = max_rim_sagitta[ua].max(delta);
        max_rim_sagitta[ub] = max_rim_sagitta[ub].max(delta);
    };

    for &[a, b] in curved_rim_edges {
        let key = [a.min(b), a.max(b)];
        let q_opt = edge_map.get(&key).copied();
        record_curved_rim_edge(key[0], key[1], q_opt);
    }
    for (&[a, b], &q) in &edge_map {
        let (ua, ub) = (a as usize, b as usize);
        if rim_verts[ua] && rim_verts[ub] {
            let (pa, pb) = (pt(ua), pt(ub));
            let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            if dist3(q, mid) > 1e-12 * dist3(pa, pb) {
                record_curved_rim_edge(a, b, Some(q));
            }
        }
    }

    // For each rim vertex touching a curved rim, also find the shortest transverse edge / altitude `h_y_min`
    // in cells sharing that curved rim so the sagitta guard `4 * xi - 1 > 2 * delta / h_y + 0.05` holds for the cell
    let mut min_transverse_h = vec![f64::INFINITY; nv];
    for &(a, b) in mesh.edges() {
        let (ua, ub) = (a as usize, b as usize);
        if rim_verts[ua] ^ rim_verts[ub] {
            let v_rim = if rim_verts[ua] { ua } else { ub };
            let len = dist3(pt(ua), pt(ub));
            if len > 1e-30 {
                min_transverse_h[v_rim] = min_transverse_h[v_rim].min(len);
            }
        }
    }

    let mut qp_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
    let mut clamped_set: BTreeSet<[u32; 2]> = BTreeSet::new();

    for &(a, b) in mesh.edges() {
        let (ua, ub) = (a as usize, b as usize);
        let key = [a.min(b), a.max(b)];
        if !(rim_verts[ua] ^ rim_verts[ub]) || edge_map.contains_key(&key) {
            continue;
        }
        let (v_rim, p_rim, p_other) = if rim_verts[ua] { (ua, pt(ua), pt(ub)) } else { (ub, pt(ub), pt(ua)) };
        let h_edge = dist3(p_rim, p_other);
        if h_edge <= 1e-30 {
            continue;
        }
        let mut xi_eff = xi_base;
        if touches_curved_rim[v_rim] {
            let delta = max_rim_sagitta[v_rim];
            let h_y = h_edge.min(min_transverse_h[v_rim]).max(1e-30);
            // Enforce strict inequality: 4 * xi_eff - 1 > 2 * delta / h_y + 0.05
            let xi_min = 0.25 * (1.0 + 2.0 * (delta / h_y) + 0.051);
            if xi_min >= 0.499 {
                report.unshifted_edges += 1;
                continue;
            }
            if xi_eff < xi_min {
                xi_eff = xi_min;
                clamped_set.insert(key);
            }
        }
        let qpt = [
            (1.0 - xi_eff) * p_rim[0] + xi_eff * p_other[0],
            (1.0 - xi_eff) * p_rim[1] + xi_eff * p_other[1],
            (1.0 - xi_eff) * p_rim[2] + xi_eff * p_other[2],
        ];
        edge_map.insert(key, qpt);
        qp_edges.insert(key);
    }

    if qp_edges.is_empty() {
        mesh.set_qp_edges(&[]);
        return report;
    }

    let sync_mesh = |mesh: &mut SimplexMesh, edge_map: &BTreeMap<[u32; 2], [f64; 3]>, qp_edges: &BTreeSet<[u32; 2]>| {
        let nodes_vec: Vec<([u32; 2], [f64; 3])> = edge_map.iter().map(|(&k, &v)| (k, v)).collect();
        mesh.set_edge_nodes(&nodes_vec);
        let qp_vec: Vec<[u32; 2]> = qp_edges.iter().copied().collect();
        mesh.set_qp_edges(&qp_vec);
    };
    sync_mesh(mesh, &edge_map, &qp_edges);

    // Validate strict interior positivity `det(J) * det_straight > 0.0` at all 216 Gauss points of `quad_rule_collapsed(8)`
    for _pass in 0..4 {
        let mut bad_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
        for c in 0..mesh.n_cells() {
            if !mesh.cell_is_curved(c) {
                continue;
            }
            let min_q8 = scaled_jacobian_min_interior_q8(mesh, c);
            if !min_q8.is_finite() || min_q8 <= 1e-12 {
                let cv = mesh.cell(c);
                for i in 0..4 {
                    for j in (i + 1)..4 {
                        let key = [cv[i].min(cv[j]), cv[i].max(cv[j])];
                        if qp_edges.contains(&key) {
                            bad_edges.insert(key);
                        }
                    }
                }
            }
        }
        if bad_edges.is_empty() {
            break;
        }
        for k in &bad_edges {
            edge_map.remove(k);
            qp_edges.remove(k);
            clamped_set.remove(k);
            report.unshifted_edges += 1;
        }
        sync_mesh(mesh, &edge_map, &qp_edges);
        if qp_edges.is_empty() {
            break;
        }
    }

    report.shifted_edges = qp_edges.len();
    report.clamped_curved_edges = clamped_set.len();
    let mut min_sj = f64::INFINITY;
    for c in 0..mesh.n_cells() {
        if mesh.cell_is_qp(c) {
            report.qp_cells += 1;
            min_sj = min_sj.min(scaled_jacobian_min_interior_q8(mesh, c));
        }
    }
    report.min_interior_jacobian_ratio = if min_sj.is_finite() { min_sj } else { 1.0 };
    report
}
