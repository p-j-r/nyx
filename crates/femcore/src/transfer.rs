// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Domain-aware cross-mesh field transfer (ADR 135): evaluate a field solved on one mesh at points, cell nodes, quadrature points or
//! boundary facets of another mesh of the same geometry, and project scalars conservatively.
//!
//! Three rules hold everywhere in this module:
//!
//! 1. **Strict domain-label filtering.** A query point carries the canonical domain name (the `geom` region name of `LabelMap`) of the
//!    target cell it belongs to, and it is evaluated only in source cells of that same canonical name. One cell-AABB BVH is built per
//!    canonical domain, so a discontinuous field (a material jump, a stress concentration at an interface) is never smeared across a
//!    material interface, whatever the two meshes look like there.
//! 2. **Direct sampling and conservative projection.** [`TransferIndex::eval_at_points`], [`TransferIndex::sample_cell_nodes`] and
//!    [`TransferIndex::sample_quadrature`] read the source interpolant (`P1` / `P2` nodal, per-cell `CellFields`, per-cell constants or a
//!    native mixed-cell Lagrange vector) directly at the target points. [`project_scalar_conservative`] solves the domain-wise `L2`
//!    projection `M_tgt Q_tgt = int N_tgt Q_src dV` with the right-hand side integrated by the source quadrature (so `sum_i b_i` is the
//!    source integral) and a final constant correction that makes `int Q_tgt dV == int Q_src dV` per domain to rounding.
//! 3. **Boundary-facet transfer.** [`BoundaryIndex`] keeps one facet BVH per canonical boundary name and evaluates a field at the closest
//!    point of the matching facets, in the cell next to the facet (or the cell of a chosen side for an interface).
//!
//! Point location: the candidates of the BVH are tried in ascending cell number and the first cell that contains the point (barycentric
//! or reference coordinates within `inside_tol`) wins, so a point on a face shared by two cells of one domain is read in the lower one.
//! Straight simplices use exact barycentric coordinates; curved (isoparametric) tetrahedra (`CurvedTet`, ADR 82) and the native mixed
//! cells (`Quad`, `Prism`, `Hex`, curved or not) invert the geometric map by Newton iteration. A point outside every cell of its domain
//! (a non-matching curved boundary) is snapped to the nearest cell of the same domain when it lies within `snap_rel * diagonal`, the
//! nearest found by a branch-and-bound BVH search with ties broken by the lower cell number.
//!
//! Determinism: queries are independent; they run in fixed chunks of points on the shared pool and are concatenated in chunk order, and
//! the projection scatters its contributions sequentially in source-cell order, so the numbers do not depend on the thread count.

use crate::fem::{quad_rule, quad_rule_collapsed, shape, CurvedTet};
use crate::fem_mixed::{lagrange_shape, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap};
use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::{Mesh, MixedDiscretisation, SimplexMesh};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap, HashMap};

// ------------------------------------------------------------------------------------------------------------ errors and options

/// Why a transfer could not be made.
#[derive(Clone, Debug, PartialEq)]
pub enum TransferError {
    /// the source and target meshes have different dimensions
    Dimension { source: usize, target: usize },
    /// an array does not fit the mesh it is given with
    Length(String),
    /// a value or coordinate is NaN or infinite
    NonFinite(String),
    /// a combination this module does not handle (a per-cell field on a mixed mesh, an order other than 1 or 2)
    Unsupported(String),
    /// the projection solver did not reach its tolerance
    NotConverged(String),
}

impl std::fmt::Display for TransferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TransferError::Dimension { source, target } => {
                write!(f, "cross-mesh transfer: the source mesh is {source}D and the target mesh {target}D")
            }
            TransferError::Length(m) => write!(f, "cross-mesh transfer: {m}"),
            TransferError::NonFinite(m) => write!(f, "cross-mesh transfer: non-finite value ({m})"),
            TransferError::Unsupported(m) => write!(f, "cross-mesh transfer: {m}"),
            TransferError::NotConverged(m) => write!(f, "cross-mesh transfer: {m}"),
        }
    }
}

impl std::error::Error for TransferError {}

impl From<TransferError> for String {
    fn from(e: TransferError) -> String {
        e.to_string()
    }
}

/// Tolerances of point location and the chunking of the parallel queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransferOptions {
    /// a point is inside a cell when its smallest barycentric (or scaled reference) coordinate is at least `-inside_tol`
    pub inside_tol: f64,
    /// a point outside every cell of its domain is evaluated in the nearest cell of that domain when it lies within `snap_rel` times the
    /// diagonal of the source bounding box; beyond that it is reported missing. `0` disables snapping
    pub snap_rel: f64,
    /// Newton step size (reference units) at which the inversion of a curved or mixed cell has converged
    pub newton_tol: f64,
    /// most Newton steps per cell
    pub newton_max_iter: usize,
    /// points per parallel chunk
    pub chunk: usize,
    /// run the chunks on the shared pool (the numbers are the same either way)
    pub parallel: bool,
}

impl Default for TransferOptions {
    fn default() -> Self {
        TransferOptions { inside_tol: 1e-10, snap_rel: 1e-3, newton_tol: 1e-14, newton_max_iter: 40, chunk: 256, parallel: true }
    }
}

// ------------------------------------------------------------------------------------------------------------ BVH

type Aabb = ([f64; 3], [f64; 3]);

#[derive(Clone, Debug)]
struct Node {
    lo: [f64; 3],
    hi: [f64; 3],
    left: u32,
    right: u32,
    start: u32,
    len: u32,
}

/// A bounding-volume hierarchy over item boxes (median split along the longest axis of the item centres, ties by item number).
#[derive(Clone, Debug, Default)]
struct Bvh {
    nodes: Vec<Node>,
    items: Vec<u32>,
}

const LEAF: usize = 8;

fn centre(b: &Aabb, k: usize) -> f64 {
    0.5 * (b.0[k] + b.1[k])
}

fn box_dist(lo: &[f64; 3], hi: &[f64; 3], p: &[f64; 3]) -> f64 {
    let mut s = 0.0;
    for k in 0..3 {
        let d = (lo[k] - p[k]).max(0.0).max(p[k] - hi[k]);
        s += d * d;
    }
    s.sqrt()
}

fn box_contains(lo: &[f64; 3], hi: &[f64; 3], p: &[f64; 3], pad: f64) -> bool {
    (0..3).all(|k| p[k] >= lo[k] - pad && p[k] <= hi[k] + pad)
}

/// A total order on finite distances for the search heap (smallest first through `Reverse`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Key(f64, u32);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Key {
    fn cmp(&self, o: &Self) -> Ordering {
        // reversed: the max-heap pops the smallest distance, then the lowest node number
        o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
    }
}

impl Bvh {
    fn build(boxes: &[Aabb], mut ids: Vec<u32>) -> Bvh {
        let mut b = Bvh { nodes: Vec::with_capacity(2 * ids.len() / LEAF + 2), items: Vec::new() };
        if !ids.is_empty() {
            b.build_rec(boxes, &mut ids, 0);
        }
        b.items = ids;
        b
    }

    fn build_rec(&mut self, boxes: &[Aabb], ids: &mut [u32], offset: usize) -> u32 {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        let mut clo = [f64::INFINITY; 3];
        let mut chi = [f64::NEG_INFINITY; 3];
        for &i in ids.iter() {
            let b = &boxes[i as usize];
            for k in 0..3 {
                lo[k] = lo[k].min(b.0[k]);
                hi[k] = hi[k].max(b.1[k]);
                let c = centre(b, k);
                clo[k] = clo[k].min(c);
                chi[k] = chi[k].max(c);
            }
        }
        let me = self.nodes.len() as u32;
        self.nodes.push(Node { lo, hi, left: 0, right: 0, start: offset as u32, len: ids.len() as u32 });
        if ids.len() <= LEAF {
            return me;
        }
        let mut axis = 0;
        for k in 1..3 {
            if chi[k] - clo[k] > chi[axis] - clo[axis] {
                axis = k;
            }
        }
        let mid = ids.len() / 2;
        ids.select_nth_unstable_by(mid, |a, b| {
            centre(&boxes[*a as usize], axis).total_cmp(&centre(&boxes[*b as usize], axis)).then(a.cmp(b))
        });
        let (l, r) = ids.split_at_mut(mid);
        let li = self.build_rec(boxes, l, offset);
        let ri = self.build_rec(boxes, r, offset + mid);
        let n = &mut self.nodes[me as usize];
        n.left = li;
        n.right = ri;
        n.len = 0;
        me
    }

    /// The items whose node boxes contain `p` (padded), in ascending item number.
    fn query(&self, p: &[f64; 3], pad: f64, out: &mut Vec<u32>) {
        out.clear();
        if self.nodes.is_empty() {
            return;
        }
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(i) = stack.pop() {
            let n = &self.nodes[i as usize];
            if !box_contains(&n.lo, &n.hi, p, pad) {
                continue;
            }
            if n.len > 0 {
                out.extend_from_slice(&self.items[n.start as usize..(n.start + n.len) as usize]);
            } else {
                stack.push(n.left);
                stack.push(n.right);
            }
        }
        out.sort_unstable();
    }

    /// The item with the smallest `dist` within `max_d` (branch and bound on the node boxes; ties go to the lower item number).
    fn nearest(&self, p: &[f64; 3], max_d: f64, mut dist: impl FnMut(u32) -> Option<f64>) -> Option<(u32, f64)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut heap = BinaryHeap::new();
        heap.push(Key(box_dist(&self.nodes[0].lo, &self.nodes[0].hi, p), 0));
        let mut best: Option<(u32, f64)> = None;
        while let Some(Key(d_box, i)) = heap.pop() {
            if d_box > max_d || best.is_some_and(|(_, bd)| d_box > bd) {
                break;
            }
            let n = &self.nodes[i as usize];
            if n.len > 0 {
                for &it in &self.items[n.start as usize..(n.start + n.len) as usize] {
                    if let Some(d) = dist(it).filter(|d| d.is_finite()) {
                        let better = match best {
                            None => true,
                            Some((bi, bd)) => d < bd || (d == bd && it < bi),
                        };
                        if better {
                            best = Some((it, d));
                        }
                    }
                }
            } else {
                for c in [n.left, n.right] {
                    let m = &self.nodes[c as usize];
                    heap.push(Key(box_dist(&m.lo, &m.hi, p), c));
                }
            }
        }
        best.filter(|(_, d)| *d <= max_d)
    }
}

// ------------------------------------------------------------------------------------------------------------ small geometry

fn sub(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn to3(p: &[f64]) -> [f64; 3] {
    let mut o = [0.0; 3];
    for (k, v) in p.iter().take(3).enumerate() {
        o[k] = *v;
    }
    o
}

fn simplex_point(mesh: &SimplexMesh, v: u32) -> [f64; 3] {
    to3(mesh.point(v as usize))
}

/// Solves `j * d = r` for the leading `dim x dim` block (`j[phys][ref]`); `None` when singular.
fn solve_small(j: &[[f64; 3]; 3], r: &[f64; 3], dim: usize) -> Option<[f64; 3]> {
    if dim == 2 {
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if !(det.is_finite() && det.abs() > 1e-300) {
            return None;
        }
        return Some([(r[0] * j[1][1] - j[0][1] * r[1]) / det, (j[0][0] * r[1] - r[0] * j[1][0]) / det, 0.0]);
    }
    let c0 = [j[0][0], j[1][0], j[2][0]];
    let c1 = [j[0][1], j[1][1], j[2][1]];
    let c2 = [j[0][2], j[1][2], j[2][2]];
    let det = dot(&c0, &cross(&c1, &c2));
    if !(det.is_finite() && det.abs() > 1e-300) {
        return None;
    }
    Some([dot(r, &cross(&c1, &c2)) / det, dot(&c0, &cross(r, &c2)) / det, dot(&c0, &cross(&c1, r)) / det])
}

/// Barycentric coordinates of `p` in the straight simplex `c`.
fn straight_bary(mesh: &SimplexMesh, c: usize, p: &[f64; 3]) -> Option<[f64; 4]> {
    let v = mesh.cell(c);
    let x0 = simplex_point(mesh, v[0]);
    let r = sub(p, &x0);
    if mesh.dim == 2 {
        let e1 = sub(&simplex_point(mesh, v[1]), &x0);
        let e2 = sub(&simplex_point(mesh, v[2]), &x0);
        let det = e1[0] * e2[1] - e1[1] * e2[0];
        if !(det.is_finite() && det.abs() > 1e-300) {
            return None;
        }
        let l1 = (r[0] * e2[1] - r[1] * e2[0]) / det;
        let l2 = (e1[0] * r[1] - e1[1] * r[0]) / det;
        Some([1.0 - l1 - l2, l1, l2, 0.0])
    } else {
        let e1 = sub(&simplex_point(mesh, v[1]), &x0);
        let e2 = sub(&simplex_point(mesh, v[2]), &x0);
        let e3 = sub(&simplex_point(mesh, v[3]), &x0);
        let det = dot(&e1, &cross(&e2, &e3));
        if !(det.is_finite() && det.abs() > 1e-300) {
            return None;
        }
        let l1 = dot(&r, &cross(&e2, &e3)) / det;
        let l2 = dot(&e1, &cross(&r, &e3)) / det;
        let l3 = dot(&e1, &cross(&e2, &r)) / det;
        Some([1.0 - l1 - l2 - l3, l1, l2, l3])
    }
}

/// Barycentric coordinates of `p` in a curved tetrahedron by Newton iteration from `start`.
fn curved_bary(ct: &CurvedTet, start: [f64; 4], p: &[f64; 3], tol: f64, max_iter: usize) -> Option<[f64; 4]> {
    let mut xi = [start[1], start[2], start[3]];
    for _ in 0..max_iter {
        let lam = [1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]];
        let r = sub(&ct.position(&lam), p);
        let d = solve_small(&ct.jacobian_matrix(&lam), &r, 3)?;
        for k in 0..3 {
            xi[k] -= d[k];
        }
        if !xi.iter().all(|x| x.is_finite() && x.abs() < 1e3) {
            return None;
        }
        if d.iter().map(|x| x.abs()).fold(0.0, f64::max) <= tol {
            return Some([1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]]);
        }
    }
    None
}

fn ref_centre(t: CellType) -> [f64; 3] {
    match t {
        CellType::Tri | CellType::Prism => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        CellType::Tet => [0.25, 0.25, 0.25],
        CellType::Quad | CellType::Hex => [0.0, 0.0, 0.0],
    }
}

/// Smallest distance of reference point `xi` to the faces of the reference cell (scaled to unit edge length); positive inside.
fn ref_margin(t: CellType, xi: &[f64; 3]) -> f64 {
    match t {
        CellType::Tri => xi[0].min(xi[1]).min(1.0 - xi[0] - xi[1]),
        CellType::Quad => 0.5 * (1.0 - xi[0].abs()).min(1.0 - xi[1].abs()),
        CellType::Tet => xi[0].min(xi[1]).min(xi[2]).min(1.0 - xi[0] - xi[1] - xi[2]),
        CellType::Prism => xi[0].min(xi[1]).min(1.0 - xi[0] - xi[1]).min(0.5 * (1.0 - xi[2].abs())),
        CellType::Hex => 0.5 * (1.0 - xi[0].abs()).min(1.0 - xi[1].abs()).min(1.0 - xi[2].abs()),
    }
}

fn clamp_tri(x: &mut f64, y: &mut f64) {
    *x = x.max(0.0);
    *y = y.max(0.0);
    let s = *x + *y;
    if s > 1.0 {
        *x /= s;
        *y /= s;
    }
}

fn ref_clamp(t: CellType, xi: &[f64; 3]) -> [f64; 3] {
    let mut o = *xi;
    match t {
        CellType::Tri => {
            let (mut a, mut b) = (o[0], o[1]);
            clamp_tri(&mut a, &mut b);
            o = [a, b, 0.0];
        }
        CellType::Quad => {
            o = [o[0].clamp(-1.0, 1.0), o[1].clamp(-1.0, 1.0), 0.0];
        }
        CellType::Tet => {
            for v in o.iter_mut() {
                *v = v.max(0.0);
            }
            let s = o[0] + o[1] + o[2];
            if s > 1.0 {
                for v in o.iter_mut() {
                    *v /= s;
                }
            }
        }
        CellType::Prism => {
            let (mut a, mut b) = (o[0], o[1]);
            clamp_tri(&mut a, &mut b);
            o = [a, b, o[2].clamp(-1.0, 1.0)];
        }
        CellType::Hex => {
            for v in o.iter_mut() {
                *v = v.clamp(-1.0, 1.0);
            }
        }
    }
    o
}

fn clamp_bary(lam: &[f64; 4], n: usize) -> [f64; 4] {
    let mut o = [0.0; 4];
    let mut s = 0.0;
    for k in 0..n {
        o[k] = lam[k].max(0.0);
        s += o[k];
    }
    if s > 0.0 {
        for v in o.iter_mut().take(n) {
            *v /= s;
        }
    } else {
        for v in o.iter_mut().take(n) {
            *v = 1.0 / n as f64;
        }
    }
    o
}

// ------------------------------------------------------------------------------------------------------------ fields

/// The mesh a [`TransferIndex`] reads from.
#[derive(Clone, Copy, Debug)]
pub enum SourceMesh<'a> {
    /// the simplex discretisation (`Tri`, `Tet`; curved tetrahedra through `edge_nodes`)
    Simplex(&'a SimplexMesh),
    /// the native mixed-cell discretisation (`Tri`, `Quad`, `Tet`, `Prism`, `Hex`, straight or curved)
    Mixed(&'a MixedDiscretisation),
}

/// A field on the source mesh.
#[derive(Clone, Copy, Debug)]
pub enum SourceField<'f> {
    /// `P1` / `P2` values per dof of the simplex source (vertices, then edge midpoints), `ncomp` interleaved per dof
    Nodal { values: &'f [f64], ncomp: usize, order: usize },
    /// un-averaged values at the local nodes of every simplex cell (`CellFields` layout, `nloc(order)` per cell)
    Cell { values: &'f [f64], order: usize },
    /// one value per cell (a material value, a region indicator)
    PerCell { values: &'f [f64] },
    /// a native mixed-cell Lagrange vector with its dof map, `ncomp` interleaved per dof
    MixedNodal { dofs: &'f MixedLagrangeDofMap, space: LagrangeSpace, values: &'f [f64], ncomp: usize },
}

impl SourceField<'_> {
    /// Components per point.
    pub fn ncomp(&self) -> usize {
        match self {
            SourceField::Nodal { ncomp, .. } | SourceField::MixedNodal { ncomp, .. } => (*ncomp).max(1),
            SourceField::Cell { .. } | SourceField::PerCell { .. } => 1,
        }
    }

    fn values(&self) -> &[f64] {
        match self {
            SourceField::Nodal { values, .. }
            | SourceField::Cell { values, .. }
            | SourceField::PerCell { values }
            | SourceField::MixedNodal { values, .. } => values,
        }
    }
}

/// Where a point was found: the cell, its local coordinates (barycentric for a simplex, reference `xi` for a mixed cell) and whether it was
/// snapped from outside the domain (with the distance it was moved).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Located {
    pub cell: usize,
    pub coords: [f64; 4],
    pub snapped: bool,
    pub distance: f64,
}

/// Values sampled at a list of points: `values[i * ncomp + k]`, `found[i]` false where no source cell of the point's domain was found
/// (the values there are 0).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Sampled {
    pub ncomp: usize,
    pub values: Vec<f64>,
    pub found: Vec<bool>,
    /// points evaluated after a snap to the nearest cell of their domain
    pub snapped: usize,
    /// points with no source cell of their domain within the snap distance
    pub missing: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Filter {
    All,
    Domain(usize),
    Absent,
}

/// The canonical domain name of a label: the name given, else `domain <label>` (the fallback the reductions use).
pub fn canonical_name(names: &BTreeMap<i32, String>, label: i32) -> String {
    names.get(&label).cloned().unwrap_or_else(|| format!("domain {label}"))
}

// ------------------------------------------------------------------------------------------------------------ the index

/// A domain-filtered cell locator over a source mesh (one cell-AABB BVH per canonical domain name and one over all cells).
pub struct TransferIndex<'a> {
    src: SourceMesh<'a>,
    dim: usize,
    opts: TransferOptions,
    boxes: Vec<Aabb>,
    cell_domain: Vec<u32>,
    names: Vec<String>,
    by_domain: Vec<Bvh>,
    all: Bvh,
    curved: Vec<Option<Box<CurvedTet>>>,
    iso: Vec<IsoparametricCell>,
    diag: f64,
}

impl<'a> TransferIndex<'a> {
    /// Builds the index; `domain_names` maps the 1-based domain labels of the mesh to their canonical names.
    pub fn new(
        src: SourceMesh<'a>,
        domain_names: &BTreeMap<i32, String>,
        opts: TransferOptions,
    ) -> Result<TransferIndex<'a>, TransferError> {
        let (dim, labels, nc): (usize, &[i32], usize) = match src {
            SourceMesh::Simplex(m) => (m.dim, &m.cell_labels, m.n_cells()),
            SourceMesh::Mixed(m) => (m.dim, &m.cell_labels, m.n_cells()),
        };
        if dim != 2 && dim != 3 {
            return Err(TransferError::Unsupported(format!("a {dim}D mesh")));
        }
        if labels.len() != nc {
            return Err(TransferError::Length(format!("{} cell labels for {nc} cells", labels.len())));
        }
        let mut curved: Vec<Option<Box<CurvedTet>>> = Vec::new();
        let mut iso: Vec<IsoparametricCell> = Vec::new();
        let mut boxes: Vec<Aabb> = Vec::with_capacity(nc);
        let grow = |b: &mut Aabb, p: &[f64; 3]| {
            for k in 0..3 {
                b.0[k] = b.0[k].min(p[k]);
                b.1[k] = b.1[k].max(p[k]);
            }
        };
        let empty: Aabb = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let pad_curved = |b: &mut Aabb| {
            for k in 0..3 {
                let e = 0.1 * (b.1[k] - b.0[k]);
                b.0[k] -= e;
                b.1[k] += e;
            }
        };
        match src {
            SourceMesh::Simplex(m) => {
                if m.has_curved() && dim == 3 {
                    curved = (0..nc).map(|c| m.cell_is_curved(c).then(|| Box::new(CurvedTet::new(m.cell_geometry_nodes(c))))).collect();
                }
                for c in 0..nc {
                    let mut b = empty;
                    match curved.get(c).and_then(|x| x.as_deref()) {
                        Some(ct) => {
                            for p in &ct.x {
                                grow(&mut b, p);
                            }
                            pad_curved(&mut b);
                        }
                        None => {
                            for &v in m.cell(c) {
                                grow(&mut b, &simplex_point(m, v));
                            }
                        }
                    }
                    boxes.push(b);
                }
            }
            SourceMesh::Mixed(m) => {
                iso.reserve(nc);
                for c in 0..nc {
                    let cell = IsoparametricCell::from_mesh(m, c);
                    let mut b = empty;
                    for p in &cell.nodes {
                        grow(&mut b, p);
                    }
                    if cell.curved {
                        pad_curved(&mut b);
                    }
                    boxes.push(b);
                    iso.push(cell);
                }
            }
        }
        if boxes.iter().any(|b| !(b.0.iter().chain(b.1.iter()).all(|x| x.is_finite()))) {
            return Err(TransferError::NonFinite("a source cell has a non-finite coordinate".into()));
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for b in &boxes {
            for k in 0..3 {
                lo[k] = lo[k].min(b.0[k]);
                hi[k] = hi[k].max(b.1[k]);
            }
        }
        let diag = if boxes.is_empty() { 0.0 } else { norm(&sub(&hi, &lo)) };
        // canonical names, sorted; labels without a name keep an identity of their own
        let mut by_label: BTreeMap<i32, String> = BTreeMap::new();
        for &l in labels {
            by_label.entry(l).or_insert_with(|| canonical_name(domain_names, l));
        }
        let mut names: Vec<String> = by_label.values().cloned().collect();
        names.sort();
        names.dedup();
        let index_of: BTreeMap<i32, u32> = by_label
            .iter()
            .filter_map(|(l, n)| names.binary_search(n).ok().map(|i| (*l, i as u32)))
            .collect();
        let cell_domain: Vec<u32> = labels.iter().map(|l| index_of.get(l).copied().unwrap_or(u32::MAX)).collect();
        let mut members: Vec<Vec<u32>> = vec![Vec::new(); names.len()];
        for (c, &d) in cell_domain.iter().enumerate() {
            if let Some(v) = members.get_mut(d as usize) {
                v.push(c as u32);
            }
        }
        let by_domain: Vec<Bvh> = members.into_iter().map(|ids| Bvh::build(&boxes, ids)).collect();
        let all = Bvh::build(&boxes, (0..nc as u32).collect());
        Ok(TransferIndex { src, dim, opts, boxes, cell_domain, names, by_domain, all, curved, iso, diag })
    }

    /// 2 or 3.
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// The canonical domain names of the source, sorted.
    pub fn domain_names(&self) -> &[String] {
        &self.names
    }

    /// The canonical domain name of source cell `c`.
    pub fn cell_domain_name(&self, c: usize) -> Option<&str> {
        self.cell_domain.get(c).and_then(|&d| self.names.get(d as usize)).map(String::as_str)
    }

    /// The diagonal of the bounding box of the source.
    pub fn diagonal(&self) -> f64 {
        self.diag
    }

    fn n_cells(&self) -> usize {
        self.boxes.len()
    }

    fn filter(&self, domain: Option<&str>) -> Filter {
        match domain {
            None => Filter::All,
            Some(n) => match self.names.binary_search_by(|x| x.as_str().cmp(n)) {
                Ok(i) => Filter::Domain(i),
                Err(_) => Filter::Absent,
            },
        }
    }

    /// Local coordinates of `p` in cell `c` (not clamped); `None` when the inversion fails.
    fn coords_of(&self, c: usize, p: &[f64; 3]) -> Option<[f64; 4]> {
        match self.src {
            SourceMesh::Simplex(m) => {
                let lam = straight_bary(m, c, p)?;
                match self.curved.get(c).and_then(|x| x.as_deref()) {
                    Some(ct) => curved_bary(ct, lam, p, self.opts.newton_tol, self.opts.newton_max_iter),
                    None => Some(lam),
                }
            }
            SourceMesh::Mixed(_) => {
                let cell = &self.iso[c];
                let mut xi = ref_centre(cell.cell_type);
                let (mut n, mut dn) = (Vec::new(), Vec::new());
                for _ in 0..self.opts.newton_max_iter {
                    let e = cell.eval_into(&xi[..], &mut n, &mut dn);
                    let r = sub(&e.x, p);
                    let d = solve_small(&e.j, &r, self.dim)?;
                    for k in 0..self.dim {
                        xi[k] -= d[k];
                    }
                    if !xi.iter().all(|x| x.is_finite() && x.abs() < 1e3) {
                        return None;
                    }
                    if d.iter().map(|x| x.abs()).fold(0.0, f64::max) <= self.opts.newton_tol {
                        return Some([xi[0], xi[1], xi[2], 0.0]);
                    }
                }
                None
            }
        }
    }

    fn margin(&self, c: usize, co: &[f64; 4]) -> f64 {
        match self.src {
            SourceMesh::Simplex(_) => co[..=self.dim].iter().copied().fold(f64::INFINITY, f64::min),
            SourceMesh::Mixed(_) => ref_margin(self.iso[c].cell_type, &[co[0], co[1], co[2]]),
        }
    }

    fn clamp(&self, c: usize, co: &[f64; 4]) -> [f64; 4] {
        match self.src {
            SourceMesh::Simplex(_) => clamp_bary(co, self.dim + 1),
            SourceMesh::Mixed(_) => {
                let x = ref_clamp(self.iso[c].cell_type, &[co[0], co[1], co[2]]);
                [x[0], x[1], x[2], 0.0]
            }
        }
    }

    /// The physical position of local coordinates `co` of cell `c`.
    fn map(&self, c: usize, co: &[f64; 4]) -> [f64; 3] {
        match self.src {
            SourceMesh::Simplex(m) => {
                if let Some(ct) = self.curved.get(c).and_then(|x| x.as_deref()) {
                    return ct.position(co);
                }
                let mut x = [0.0; 3];
                for (k, &v) in m.cell(c).iter().enumerate() {
                    let p = simplex_point(m, v);
                    for a in 0..3 {
                        x[a] += co[k] * p[a];
                    }
                }
                x
            }
            SourceMesh::Mixed(_) => self.iso[c].eval(&co[..3]).x,
        }
    }

    fn locate_in(&self, p: &[f64; 3], filter: Filter, snap: f64, cand: &mut Vec<u32>) -> Option<Located> {
        let bvh = match filter {
            Filter::All => &self.all,
            Filter::Domain(i) => &self.by_domain[i],
            Filter::Absent => return None,
        };
        let pad = 1e-12 * self.diag;
        bvh.query(p, pad, cand);
        for &c in cand.iter() {
            let c = c as usize;
            let b = &self.boxes[c];
            if !box_contains(&b.0, &b.1, p, pad) {
                continue;
            }
            if let Some(co) = self.coords_of(c, p) {
                if self.margin(c, &co) >= -self.opts.inside_tol {
                    return Some(Located { cell: c, coords: co, snapped: false, distance: 0.0 });
                }
            }
        }
        if !(snap > 0.0) {
            return None;
        }
        let found = bvh.nearest(p, snap, |c| {
            let c = c as usize;
            let co = self.clamp(c, &self.coords_of(c, p)?);
            Some(norm(&sub(&self.map(c, &co), p)))
        })?;
        let c = found.0 as usize;
        let co = self.clamp(c, &self.coords_of(c, p)?);
        Some(Located { cell: c, coords: co, snapped: true, distance: found.1 })
    }

    /// Finds the cell of `p` (mesh axes, `dim` values) among the source cells of canonical domain `domain` (all cells when `None`).
    pub fn locate(&self, p: &[f64], domain: Option<&str>) -> Option<Located> {
        let mut cand = Vec::new();
        self.locate_in(&to3(p), self.filter(domain), self.opts.snap_rel * self.diag, &mut cand)
    }

    /// Checks that `field` fits the source mesh and holds finite numbers.
    pub fn check_field(&self, field: &SourceField) -> Result<(), TransferError> {
        let nc = self.n_cells();
        let fit = |what: &str, have: usize, want: usize| -> Result<(), TransferError> {
            if have != want {
                return Err(TransferError::Length(format!("{what}: {have} values where the source mesh needs {want}")));
            }
            Ok(())
        };
        match (field, self.src) {
            (SourceField::Nodal { values, ncomp, order }, SourceMesh::Simplex(m)) => {
                if !matches!(order, 1 | 2) {
                    return Err(TransferError::Unsupported(format!("element order {order}")));
                }
                fit("nodal field", values.len(), (*ncomp).max(1) * m.n_dofs(*order))?;
            }
            (SourceField::Cell { values, order }, SourceMesh::Simplex(m)) => {
                if !matches!(order, 1 | 2) {
                    return Err(TransferError::Unsupported(format!("element order {order}")));
                }
                fit("per-cell nodal field", values.len(), nc * m.nloc(*order))?;
            }
            (SourceField::PerCell { values }, _) => fit("per-cell field", values.len(), nc)?,
            (SourceField::MixedNodal { dofs, values, ncomp, .. }, SourceMesh::Mixed(_)) => {
                let n = values.len() / (*ncomp).max(1);
                if !values.len().is_multiple_of((*ncomp).max(1)) {
                    return Err(TransferError::Length("mixed nodal field: length is not a multiple of the components".into()));
                }
                for c in 0..nc {
                    if dofs.cell_dofs(c).iter().any(|&d| d >= n) {
                        return Err(TransferError::Length(format!("mixed nodal field: {n} dofs do not cover cell {c}")));
                    }
                }
            }
            (SourceField::Nodal { .. } | SourceField::Cell { .. }, SourceMesh::Mixed(_)) => {
                return Err(TransferError::Unsupported("a simplex nodal or per-cell nodal field on a mixed-cell mesh".into()))
            }
            (SourceField::MixedNodal { .. }, SourceMesh::Simplex(_)) => {
                return Err(TransferError::Unsupported("a mixed-cell Lagrange field on a simplex mesh".into()))
            }
        }
        if let Some(i) = field.values().iter().position(|v| !v.is_finite()) {
            return Err(TransferError::NonFinite(format!("source value {i}")));
        }
        Ok(())
    }

    /// The field at a located point (`out` has `ncomp` entries). The field must have passed [`TransferIndex::check_field`].
    pub fn eval_located(&self, field: &SourceField, loc: &Located, out: &mut [f64]) {
        let c = loc.cell;
        out.iter_mut().for_each(|o| *o = 0.0);
        match (field, self.src) {
            (SourceField::Nodal { values, ncomp, order }, SourceMesh::Simplex(m)) => {
                let ncomp = (*ncomp).max(1);
                let (n, _) = shape(&loc.coords[..=self.dim], *order);
                let mut dofs = Vec::with_capacity(n.len());
                m.cell_dofs(c, *order, &mut dofs);
                for (a, &d) in dofs.iter().enumerate() {
                    for k in 0..ncomp.min(out.len()) {
                        out[k] += n[a] * values[d * ncomp + k];
                    }
                }
            }
            (SourceField::Cell { values, order }, SourceMesh::Simplex(_)) => {
                let (n, _) = shape(&loc.coords[..=self.dim], *order);
                let base = c * n.len();
                if let Some(o) = out.first_mut() {
                    *o = n.iter().enumerate().map(|(a, w)| w * values[base + a]).sum();
                }
            }
            (SourceField::PerCell { values }, _) => {
                if let Some(o) = out.first_mut() {
                    *o = values[c];
                }
            }
            (SourceField::MixedNodal { dofs, space, values, ncomp }, SourceMesh::Mixed(_)) => {
                let ncomp = (*ncomp).max(1);
                let (n, _) = lagrange_shape(self.iso[c].cell_type, *space, &loc.coords[..3]);
                for (a, &d) in dofs.cell_dofs(c).iter().enumerate() {
                    for k in 0..ncomp.min(out.len()) {
                        out[k] += n[a] * values[d * ncomp + k];
                    }
                }
            }
            _ => {}
        }
    }

    /// Evaluates at `n` points given by `point(i) -> (position, filter)`, in fixed chunks.
    fn eval_many<P>(&self, field: &SourceField, n: usize, point: P) -> Result<Sampled, TransferError>
    where
        P: Fn(usize) -> ([f64; 3], Filter) + Sync + Send,
    {
        self.check_field(field)?;
        let ncomp = field.ncomp();
        let chunk = self.opts.chunk.max(1);
        let nchunks = n.div_ceil(chunk);
        let snap = self.opts.snap_rel * self.diag;
        type Part = (Vec<f64>, Vec<bool>, usize, usize);
        let parts: Vec<Part> = crate::par::chunk_map(nchunks, self.opts.parallel && nchunks > 1, |ci| {
            let (s, e) = (ci * chunk, ((ci + 1) * chunk).min(n));
            let mut vals = vec![0.0; (e - s) * ncomp];
            let mut found = vec![false; e - s];
            let (mut snapped, mut missing) = (0usize, 0usize);
            let mut cand = Vec::new();
            for i in s..e {
                let (p, f) = point(i);
                match self.locate_in(&p, f, snap, &mut cand) {
                    Some(loc) => {
                        self.eval_located(field, &loc, &mut vals[(i - s) * ncomp..(i - s + 1) * ncomp]);
                        found[i - s] = true;
                        snapped += loc.snapped as usize;
                    }
                    None => missing += 1,
                }
            }
            (vals, found, snapped, missing)
        });
        let mut out = Sampled { ncomp, values: Vec::with_capacity(n * ncomp), found: Vec::with_capacity(n), snapped: 0, missing: 0 };
        for (v, f, s, m) in parts {
            out.values.extend(v);
            out.found.extend(f);
            out.snapped += s;
            out.missing += m;
        }
        if let Some(i) = out.values.iter().position(|v| !v.is_finite()) {
            return Err(TransferError::NonFinite(format!("sampled value {i}")));
        }
        Ok(out)
    }

    /// Evaluates `field` at `points` (flat, `dim` per point). `domains` is empty (no filtering) or holds one canonical domain name per point.
    pub fn eval_at_points(
        &self,
        field: &SourceField,
        points: &[f64],
        domains: &[Option<&str>],
    ) -> Result<Sampled, TransferError> {
        let d = self.dim;
        if !points.len().is_multiple_of(d) {
            return Err(TransferError::Length(format!("{} coordinates are not points of dimension {d}", points.len())));
        }
        if let Some(i) = points.iter().position(|v| !v.is_finite()) {
            return Err(TransferError::NonFinite(format!("query coordinate {i}")));
        }
        let n = points.len() / d;
        if !domains.is_empty() && domains.len() != n {
            return Err(TransferError::Length(format!("{} domain names for {n} points", domains.len())));
        }
        let filters: Vec<Filter> =
            if domains.is_empty() { vec![Filter::All; n] } else { domains.iter().map(|x| self.filter(*x)).collect() };
        self.eval_many(field, n, |i| (to3(&points[i * d..(i + 1) * d]), filters[i]))
    }

    /// The scalar `field` at `p` in domain `domain`; `Ok(None)` outside.
    pub fn eval_scalar(&self, field: &SourceField, p: &[f64], domain: Option<&str>) -> Result<Option<f64>, TransferError> {
        Ok(self.eval_vector(field, p, domain)?.and_then(|v| v.first().copied()))
    }

    /// All components of `field` at `p` in domain `domain`; `Ok(None)` outside.
    pub fn eval_vector(&self, field: &SourceField, p: &[f64], domain: Option<&str>) -> Result<Option<Vec<f64>>, TransferError> {
        let s = self.eval_at_points(field, &p[..self.dim.min(p.len())], &[domain])?;
        Ok(s.found.first().copied().unwrap_or(false).then_some(s.values))
    }

    /// Per target label, the filter of its canonical name.
    fn target_filters(&self, target: &SimplexMesh, names: &BTreeMap<i32, String>) -> BTreeMap<i32, Filter> {
        let mut out = BTreeMap::new();
        for &l in &target.cell_labels {
            out.entry(l).or_insert_with(|| {
                let n = canonical_name(names, l);
                self.filter(Some(&n))
            });
        }
        out
    }

    fn check_target(&self, target: &SimplexMesh) -> Result<(), TransferError> {
        if target.dim != self.dim {
            return Err(TransferError::Dimension { source: self.dim, target: target.dim });
        }
        Ok(())
    }

    /// `field` at the local nodes of every target cell (`nloc(order)` per cell, the `CellFields` layout), each node evaluated only in the
    /// source cells of its target cell's canonical domain: a field that jumps across an interface keeps its two values at the interface.
    pub fn sample_cell_nodes(
        &self,
        field: &SourceField,
        target: &SimplexMesh,
        target_names: &BTreeMap<i32, String>,
        order: usize,
    ) -> Result<Sampled, TransferError> {
        self.check_target(target)?;
        if !matches!(order, 1 | 2) {
            return Err(TransferError::Unsupported(format!("element order {order}")));
        }
        let filters = self.target_filters(target, target_names);
        let nloc = target.nloc(order);
        let d = target.dim;
        self.eval_many(field, target.n_cells() * nloc, |i| {
            let (c, a) = (i / nloc, i % nloc);
            let f = filters.get(&target.cell_labels[c]).copied().unwrap_or(Filter::Absent);
            (target_node(target, c, a, d), f)
        })
    }

    /// `field` at the dofs of the target (`P1`/`P2`, `ncomp` interleaved): the mean of the per-cell node samples over the cells around each
    /// dof that found a value. A dof where no cell found one is 0 and `found` false.
    pub fn sample_dofs(
        &self,
        field: &SourceField,
        target: &SimplexMesh,
        target_names: &BTreeMap<i32, String>,
        order: usize,
    ) -> Result<Sampled, TransferError> {
        let s = self.sample_cell_nodes(field, target, target_names, order)?;
        let (ncomp, nloc, n) = (s.ncomp, target.nloc(order), target.n_dofs(order));
        let mut acc = vec![0.0; n * ncomp];
        let mut cnt = vec![0usize; n];
        let mut dofs = Vec::new();
        for c in 0..target.n_cells() {
            target.cell_dofs(c, order, &mut dofs);
            for (a, &dof) in dofs.iter().enumerate() {
                let i = c * nloc + a;
                if !s.found[i] {
                    continue;
                }
                cnt[dof] += 1;
                for k in 0..ncomp {
                    acc[dof * ncomp + k] += s.values[i * ncomp + k];
                }
            }
        }
        let mut missing = 0;
        for dof in 0..n {
            if cnt[dof] > 0 {
                for k in 0..ncomp {
                    acc[dof * ncomp + k] /= cnt[dof] as f64;
                }
            } else {
                missing += 1;
            }
        }
        Ok(Sampled { ncomp, values: acc, found: cnt.iter().map(|&k| k > 0).collect(), snapped: s.snapped, missing })
    }

    /// `field` at the points of the degree-`degree` rule of `fem::quad_rule` in every target cell (`values[(c * nq + q) * ncomp + k]`), the
    /// direct quadrature-point sampling a cross-mesh integral uses. Curved target tetrahedra map their points through the curved geometry.
    pub fn sample_quadrature(
        &self,
        field: &SourceField,
        target: &SimplexMesh,
        target_names: &BTreeMap<i32, String>,
        degree: usize,
    ) -> Result<Sampled, TransferError> {
        self.check_target(target)?;
        let filters = self.target_filters(target, target_names);
        let rule = quad_rule(target.dim, degree);
        let nq = rule.w.len();
        self.eval_many(field, target.n_cells() * nq, |i| {
            let (c, q) = (i / nq, i % nq);
            let f = filters.get(&target.cell_labels[c]).copied().unwrap_or(Filter::Absent);
            (target_point(target, c, &rule.lam[q]), f)
        })
    }
}

/// The position of local node `a` (vertices, then edge nodes in local pair order) of target cell `c`; curved edges give their geometry node.
fn target_node(target: &SimplexMesh, c: usize, a: usize, d: usize) -> [f64; 3] {
    let v = target.cell(c);
    if a <= d {
        return simplex_point(target, v[a]);
    }
    if d == 3 && target.cell_is_curved(c) {
        return target.cell_geometry_nodes(c)[a];
    }
    let mut k = d + 1;
    for i in 0..=d {
        for j in i + 1..=d {
            if k == a {
                let (p, q) = (simplex_point(target, v[i]), simplex_point(target, v[j]));
                return [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])];
            }
            k += 1;
        }
    }
    simplex_point(target, v[0])
}

/// The position of barycentric point `lam` of target cell `c` (through the curved geometry of a curved tetrahedron).
fn target_point(target: &SimplexMesh, c: usize, lam: &[f64]) -> [f64; 3] {
    if target.dim == 3 && target.cell_is_curved(c) {
        return CurvedTet::new(target.cell_geometry_nodes(c)).position(&[lam[0], lam[1], lam[2], lam[3]]);
    }
    let mut x = [0.0; 3];
    for (k, &v) in target.cell(c).iter().enumerate() {
        let p = simplex_point(target, v);
        for a in 0..3 {
            x[a] += lam[k] * p[a];
        }
    }
    x
}

// ------------------------------------------------------------------------------------------------------------ conservative projection

/// The mass matrix of [`project_scalar_conservative`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MassKind {
    /// the consistent mass matrix (solved by Jacobi-preconditioned conjugate gradients)
    Consistent,
    /// the row-scaled diagonal (HRZ) lumped mass matrix, positive for `P2` too
    Lumped,
}

/// The outcome of a conservative projection: per-cell nodal values on the target (`CellFields` layout, continuous inside a domain and
/// discontinuous across domains) and the domain integrals before and after.
#[derive(Clone, Debug, PartialEq)]
pub struct Projection {
    pub order: usize,
    pub nloc: usize,
    pub cell_values: Vec<f64>,
    /// `int Q_src dV` per canonical domain, by the source quadrature
    pub source_integrals: BTreeMap<String, f64>,
    /// `int Q_tgt dV` per canonical domain, by the target quadrature
    pub target_integrals: BTreeMap<String, f64>,
    /// the source integral of domains the target does not have (it cannot be conserved)
    pub unmatched_source_integral: f64,
    /// conjugate-gradient iterations summed over the domains (0 for a lumped mass)
    pub iterations: usize,
}

/// A target cell's rule: barycentric points and absolute weights (curved tetrahedra use the collapsed degree-8 rule with `|det J| / 6`).
fn cell_rule(mesh: &SimplexMesh, c: usize, straight: &crate::fem::Quad, collapsed: &crate::fem::Quad) -> Vec<([f64; 4], f64)> {
    let d = mesh.dim;
    if d == 3 && mesh.cell_is_curved(c) {
        let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
        return collapsed
            .lam
            .iter()
            .zip(&collapsed.w)
            .map(|(l, w)| {
                let lam = [l[0], l[1], l[2], l[3]];
                (lam, w * ct.jacobian(&lam).1.abs() / 6.0)
            })
            .collect();
    }
    let vol = simplex_volume(mesh, c);
    straight
        .lam
        .iter()
        .zip(&straight.w)
        .map(|(l, w)| {
            let mut lam = [0.0; 4];
            lam[..=d].copy_from_slice(&l[..=d]);
            (lam, w * vol)
        })
        .collect()
}

fn simplex_volume(mesh: &SimplexMesh, c: usize) -> f64 {
    let v = mesh.cell(c);
    let x0 = simplex_point(mesh, v[0]);
    let e1 = sub(&simplex_point(mesh, v[1]), &x0);
    let e2 = sub(&simplex_point(mesh, v[2]), &x0);
    if mesh.dim == 2 {
        (e1[0] * e2[1] - e1[1] * e2[0]).abs() / 2.0
    } else {
        let e3 = sub(&simplex_point(mesh, v[3]), &x0);
        dot(&e1, &cross(&e2, &e3)).abs() / 6.0
    }
}

/// A small symmetric sparse matrix in compressed rows.
struct SpMat {
    ptr: Vec<usize>,
    col: Vec<u32>,
    val: Vec<f64>,
}

impl SpMat {
    fn from_triplets(n: usize, mut t: Vec<(u32, u32, f64)>) -> SpMat {
        // stable sort: equal positions keep their insertion order, so the sums are the same every run
        t.sort_by_key(|a| (a.0, a.1));
        let mut ptr = vec![0usize; n + 1];
        let mut col = Vec::with_capacity(t.len());
        let mut val: Vec<f64> = Vec::with_capacity(t.len());
        let mut last: Option<(u32, u32)> = None;
        for (r, c, v) in t {
            if last == Some((r, c)) {
                if let Some(x) = val.last_mut() {
                    *x += v;
                }
            } else {
                col.push(c);
                val.push(v);
                ptr[r as usize + 1] += 1;
                last = Some((r, c));
            }
        }
        for i in 0..n {
            ptr[i + 1] += ptr[i];
        }
        SpMat { ptr, col, val }
    }

    fn mul(&self, x: &[f64], y: &mut [f64]) {
        for (i, yi) in y.iter_mut().enumerate() {
            let mut s = 0.0;
            for k in self.ptr[i]..self.ptr[i + 1] {
                s += self.val[k] * x[self.col[k] as usize];
            }
            *yi = s;
        }
    }
}

/// Jacobi-preconditioned conjugate gradients to `||r|| <= tol ||b||`.
fn pcg(a: &SpMat, diag: &[f64], b: &[f64], tol: f64, max_iter: usize) -> (Vec<f64>, usize, f64) {
    let n = b.len();
    let bn = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut x = vec![0.0; n];
    if bn == 0.0 {
        return (x, 0, 0.0);
    }
    let mut r = b.to_vec();
    let mut z: Vec<f64> = r.iter().zip(diag).map(|(r, d)| r / d).collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    let mut ap = vec![0.0; n];
    let mut rel = 1.0;
    for it in 0..max_iter {
        a.mul(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        if !(pap > 0.0) {
            return (x, it, rel);
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        rel = r.iter().map(|v| v * v).sum::<f64>().sqrt() / bn;
        if rel <= tol {
            return (x, it + 1, rel);
        }
        for i in 0..n {
            z[i] = r[i] / diag[i];
        }
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    (x, max_iter, rel)
}

/// Conservative domain-wise `L2` projection of a scalar source field onto the `P1`/`P2` space of the target, per canonical domain:
/// `M_tgt Q_tgt = int N_tgt Q_src dV`. The right-hand side is integrated on the target cells with the source sampled directly at their
/// quadrature points (each point read only in source cells of its own domain, snapped to the nearest one at a non-matching boundary), so a
/// source that is a polynomial of the target order inside each target cell is reproduced exactly. The source integral of every domain is
/// taken separately by the source quadrature (exact for a `P1`/`P2` source), and after the solve a constant correction makes
/// `int Q_tgt dV` equal to it to rounding. The target space of each domain is continuous inside the domain and independent of the other
/// domains, so a field that jumps at an interface keeps its jump.
pub fn project_scalar_conservative(
    source: &TransferIndex,
    field: &SourceField,
    target: &TransferIndex,
    order: usize,
    mass: MassKind,
) -> Result<Projection, TransferError> {
    let SourceMesh::Simplex(tm) = target.src else {
        return Err(TransferError::Unsupported("a conservative projection onto a mixed-cell target".into()));
    };
    if source.dim != target.dim {
        return Err(TransferError::Dimension { source: source.dim, target: target.dim });
    }
    if !matches!(order, 1 | 2) {
        return Err(TransferError::Unsupported(format!("element order {order}")));
    }
    source.check_field(field)?;
    if field.ncomp() != 1 {
        return Err(TransferError::Unsupported("a conservative projection of a vector field (project each component)".into()));
    }
    let d = source.dim;
    // 1. the source integral of every source domain, by the source quadrature, in chunks of source cells summed in chunk order
    let nsc = source.n_cells();
    let chunk = 256usize;
    let nchunks = nsc.div_ceil(chunk);
    let straight = quad_rule(d, 4);
    let collapsed = quad_rule_collapsed(8);
    let nd = source.names.len();
    let parts: Vec<Vec<f64>> = crate::par::chunk_map(nchunks, source.opts.parallel && nchunks > 1, |ci| {
        let mut integrals = vec![0.0; nd];
        let mut val = [0.0];
        for c in ci * chunk..((ci + 1) * chunk).min(nsc) {
            let dom = source.cell_domain[c] as usize;
            let pts: Vec<([f64; 4], f64)> = match source.src {
                SourceMesh::Simplex(m) => cell_rule(m, c, &straight, &collapsed),
                SourceMesh::Mixed(_) => {
                    let cell = &source.iso[c];
                    let rule = quad_rule_mixed(cell.cell_type, if cell.curved { 6 } else { 4 });
                    rule.pts
                        .iter()
                        .zip(&rule.w)
                        .map(|(xi, w)| ([xi[0], xi[1], xi[2], 0.0], w * cell.eval(&xi[..]).det.abs()))
                        .collect()
                }
            };
            for (co, w) in pts {
                source.eval_located(field, &Located { cell: c, coords: co, snapped: false, distance: 0.0 }, &mut val);
                if let Some(s) = integrals.get_mut(dom) {
                    *s += w * val[0];
                }
            }
        }
        integrals
    });
    let mut source_int = vec![0.0; nd];
    for i in parts {
        for k in 0..nd {
            source_int[k] += i[k];
        }
    }
    let unmatched: f64 = source
        .names
        .iter()
        .zip(&source_int)
        .filter(|(n, _)| target.filter(Some(n)) == Filter::Absent)
        .map(|(_, v)| *v)
        .sum();
    // 2. the source at the quadrature points of every target cell, read in source cells of the target cell's domain
    let ntc = tm.n_cells();
    let sfilter: Vec<Filter> = target.names.iter().map(|n| source.filter(Some(n))).collect();
    let tchunks = ntc.div_ceil(chunk);
    let qparts: Vec<Vec<Vec<f64>>> = crate::par::chunk_map(tchunks, source.opts.parallel && tchunks > 1, |ci| {
        let mut cand = Vec::new();
        let mut val = [0.0];
        (ci * chunk..((ci + 1) * chunk).min(ntc))
            .map(|c| {
                let f = sfilter.get(target.cell_domain[c] as usize).copied().unwrap_or(Filter::Absent);
                cell_rule(tm, c, &straight, &collapsed)
                    .iter()
                    .map(|(lam, _)| {
                        let x = target_point(tm, c, &lam[..=d]);
                        match source.locate_in(&x, f, f64::INFINITY, &mut cand) {
                            Some(loc) => {
                                source.eval_located(field, &loc, &mut val);
                                val[0]
                            }
                            None => 0.0,
                        }
                    })
                    .collect()
            })
            .collect()
    });
    let qvals: Vec<Vec<f64>> = qparts.into_iter().flatten().collect();
    // 3. per target domain: the mass matrix of its own continuous space, the right-hand side, the solve and the correction
    let nloc = tm.nloc(order);
    let mut cell_values = vec![0.0; ntc * nloc];
    let mut target_int: BTreeMap<String, f64> = BTreeMap::new();
    let mut iterations = 0usize;
    let ndof = tm.n_dofs(order);
    let mut local = vec![u32::MAX; ndof];
    let mut dofs = Vec::new();
    for (td, tname) in target.names.iter().enumerate() {
        let cells: Vec<usize> = (0..ntc).filter(|&c| target.cell_domain[c] as usize == td).collect();
        let mut gl: Vec<usize> = Vec::new();
        for &c in &cells {
            tm.cell_dofs(c, order, &mut dofs);
            for &g in &dofs {
                if local[g] == u32::MAX {
                    local[g] = gl.len() as u32;
                    gl.push(g);
                }
            }
        }
        let n = gl.len();
        let mut trip: Vec<(u32, u32, f64)> = Vec::with_capacity(cells.len() * nloc * nloc);
        let mut lumped = vec![0.0; n];
        let mut m1 = vec![0.0; n];
        let mut ldofs = vec![0u32; nloc];
        for &c in &cells {
            tm.cell_dofs(c, order, &mut dofs);
            for (a, &g) in dofs.iter().enumerate() {
                ldofs[a] = local[g];
            }
            let mut me = vec![0.0; nloc * nloc];
            let mut vol = 0.0;
            for (lam, w) in cell_rule(tm, c, &straight, &collapsed) {
                let (nn, _) = shape(&lam[..=d], order);
                vol += w;
                for a in 0..nloc {
                    for b in 0..nloc {
                        me[a * nloc + b] += w * nn[a] * nn[b];
                    }
                }
            }
            let trace: f64 = (0..nloc).map(|a| me[a * nloc + a]).sum();
            for a in 0..nloc {
                let la = ldofs[a] as usize;
                lumped[la] += if trace > 0.0 { me[a * nloc + a] * vol / trace } else { 0.0 };
                for b in 0..nloc {
                    m1[la] += me[a * nloc + b];
                    trip.push((ldofs[a], ldofs[b], me[a * nloc + b]));
                }
            }
        }
        let mut rhs = vec![0.0; n];
        for &c in &cells {
            tm.cell_dofs(c, order, &mut dofs);
            for ((lam, w), qv) in cell_rule(tm, c, &straight, &collapsed).iter().zip(&qvals[c]) {
                let (nn, _) = shape(&lam[..=d], order);
                for (a, &g) in dofs.iter().enumerate() {
                    rhs[local[g] as usize] += w * qv * nn[a];
                }
            }
        }
        // the integral the domain must keep: the source integral of the same canonical domain
        let source_total = match sfilter.get(td) {
            Some(Filter::Domain(si)) => source_int.get(*si).copied(),
            _ => None,
        };
        let mut q = match mass {
            MassKind::Lumped => rhs.iter().zip(&lumped).map(|(b, m)| if *m > 0.0 { b / m } else { 0.0 }).collect::<Vec<f64>>(),
            MassKind::Consistent => {
                let a = SpMat::from_triplets(n, trip);
                let diag: Vec<f64> = (0..n)
                    .map(|i| {
                        (a.ptr[i]..a.ptr[i + 1]).find(|&k| a.col[k] as usize == i).map_or(1.0, |k| a.val[k]).max(1e-300)
                    })
                    .collect();
                let (x, it, rel) = pcg(&a, &diag, &rhs, 1e-14, 20 * n + 200);
                iterations += it;
                if !(rel <= 1e-10) {
                    return Err(TransferError::NotConverged(format!(
                        "the mass solve of domain '{tname}' stopped at relative residual {rel:e} after {it} iterations"
                    )));
                }
                x
            }
        };
        // conservation: int Q dV = m1 . Q must equal the source integral of the domain
        let msum: f64 = m1.iter().sum();
        let integral: f64 = m1.iter().zip(&q).map(|(m, v)| m * v).sum();
        if let (true, Some(total)) = (msum > 0.0, source_total) {
            let shift = (total - integral) / msum;
            for v in q.iter_mut() {
                *v += shift;
            }
        }
        let mut tint = 0.0;
        for &c in &cells {
            tm.cell_dofs(c, order, &mut dofs);
            for (a, &g) in dofs.iter().enumerate() {
                cell_values[c * nloc + a] = q[local[g] as usize];
            }
            for (lam, w) in cell_rule(tm, c, &straight, &collapsed) {
                let (nn, _) = shape(&lam[..=d], order);
                tint += w * nn.iter().zip(&dofs).map(|(nv, &g)| nv * q[local[g] as usize]).sum::<f64>();
            }
        }
        target_int.insert(tname.clone(), tint);
        for &g in &gl {
            local[g] = u32::MAX;
        }
    }
    if let Some(i) = cell_values.iter().position(|v| !v.is_finite()) {
        return Err(TransferError::NonFinite(format!("projected value {i}")));
    }
    let source_integrals: BTreeMap<String, f64> = source.names.iter().cloned().zip(source_int).collect();
    Ok(Projection {
        order,
        nloc,
        cell_values,
        source_integrals,
        target_integrals: target_int,
        unmatched_source_integral: unmatched,
        iterations,
    })
}

// ------------------------------------------------------------------------------------------------------------ boundary facets

fn closest_on_segment(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let l2 = dot(&ab, &ab);
    let t = if l2 > 0.0 { (dot(&sub(p, a), &ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]]
}

/// The closest point of triangle `abc` to `p` (Voronoi regions of the vertices, edges and face).
fn closest_on_triangle(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(&ab, &ap);
    let d2 = dot(&ac, &ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return *a;
    }
    let bp = sub(p, b);
    let d3 = dot(&ab, &bp);
    let d4 = dot(&ac, &bp);
    if d3 >= 0.0 && d4 <= d3 {
        return *b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return [a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]];
    }
    let cp = sub(p, c);
    let d5 = dot(&ab, &cp);
    let d6 = dot(&ac, &cp);
    if d6 >= 0.0 && d5 <= d6 {
        return *c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return [a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]];
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let bc = sub(c, b);
        return [b[0] + w * bc[0], b[1] + w * bc[1], b[2] + w * bc[2]];
    }
    let denom = va + vb + vc;
    if !(denom.abs() > 0.0) {
        return *a;
    }
    let v = vb / denom;
    let w = vc / denom;
    [a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w]
}

/// Where a point was projected onto a boundary: the facet, the closest point, its distance, and the cell it is read in with its
/// barycentric coordinates there.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryHit {
    pub facet: usize,
    pub point: [f64; 3],
    pub distance: f64,
    pub cell: usize,
    pub coords: [f64; 4],
}

/// A closest-point locator over the boundary facets of a simplex mesh, one facet BVH per canonical boundary name (the `LabelMap` boundary
/// names; a name standing for several labels covers all of them).
pub struct BoundaryIndex<'a> {
    mesh: &'a SimplexMesh,
    names: Vec<String>,
    by_name: Vec<Bvh>,
    /// the cells next to each facet, ascending (`u32::MAX` where there is none)
    facet_cells: Vec<[u32; 2]>,
    curved: BTreeMap<usize, CurvedTet>,
    opts: TransferOptions,
}

fn face_key(v: &[u32]) -> [u32; 3] {
    let mut k = [u32::MAX; 3];
    for (i, &x) in v.iter().take(3).enumerate() {
        k[i] = x;
    }
    k[..v.len().min(3)].sort_unstable();
    k
}

impl<'a> BoundaryIndex<'a> {
    /// Builds the index; `boundary_names` maps the 1-based boundary labels to their canonical names (a label without one is `boundary <n>`).
    pub fn new(
        mesh: &'a SimplexMesh,
        boundary_names: &BTreeMap<i32, String>,
        opts: TransferOptions,
    ) -> Result<BoundaryIndex<'a>, TransferError> {
        let d = mesh.dim;
        let nf = mesh.n_facets();
        if mesh.facet_labels.len() != nf {
            return Err(TransferError::Length(format!("{} facet labels for {nf} facets", mesh.facet_labels.len())));
        }
        let mut face_cells: HashMap<[u32; 3], [u32; 2]> = HashMap::with_capacity(nf);
        for f in 0..nf {
            face_cells.insert(face_key(mesh.facet(f)), [u32::MAX; 2]);
        }
        let nloc = d + 1;
        let mut fbuf = [0u32; 3];
        for c in 0..mesh.n_cells() {
            let vs = mesh.cell(c);
            for skip in 0..nloc {
                let mut k = 0;
                for (i, &x) in vs.iter().enumerate() {
                    if i != skip {
                        fbuf[k] = x;
                        k += 1;
                    }
                }
                if let Some(e) = face_cells.get_mut(&face_key(&fbuf[..d])) {
                    if e[0] == u32::MAX {
                        e[0] = c as u32;
                    } else if e[1] == u32::MAX {
                        e[1] = c as u32;
                    }
                }
            }
        }
        let facet_cells: Vec<[u32; 2]> = (0..nf)
            .map(|f| {
                let mut e = face_cells.get(&face_key(mesh.facet(f))).copied().unwrap_or([u32::MAX; 2]);
                e.sort_unstable();
                e
            })
            .collect();
        let name_of = |l: i32| boundary_names.get(&l).cloned().unwrap_or_else(|| format!("boundary {l}"));
        let mut names: Vec<String> = mesh.facet_labels.iter().map(|&l| name_of(l)).collect();
        names.sort();
        names.dedup();
        let mut boxes: Vec<Aabb> = Vec::with_capacity(nf);
        let mut members: Vec<Vec<u32>> = vec![Vec::new(); names.len()];
        for f in 0..nf {
            let mut b: Aabb = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for &v in mesh.facet(f) {
                let p = simplex_point(mesh, v);
                for k in 0..3 {
                    b.0[k] = b.0[k].min(p[k]);
                    b.1[k] = b.1[k].max(p[k]);
                }
            }
            if !b.0.iter().chain(b.1.iter()).all(|x| x.is_finite()) {
                return Err(TransferError::NonFinite(format!("facet {f} has a non-finite coordinate")));
            }
            boxes.push(b);
            if let Ok(i) = names.binary_search(&name_of(mesh.facet_labels[f])) {
                members[i].push(f as u32);
            }
        }
        let by_name = members.into_iter().map(|ids| Bvh::build(&boxes, ids)).collect();
        let mut curved = BTreeMap::new();
        if d == 3 && mesh.has_curved() {
            for cells in &facet_cells {
                for &c in cells.iter().filter(|&&c| c != u32::MAX) {
                    let c = c as usize;
                    if mesh.cell_is_curved(c) {
                        curved.entry(c).or_insert_with(|| CurvedTet::new(mesh.cell_geometry_nodes(c)));
                    }
                }
            }
        }
        Ok(BoundaryIndex { mesh, names, by_name, facet_cells, curved, opts })
    }

    /// The canonical boundary names, sorted.
    pub fn boundary_names(&self) -> &[String] {
        &self.names
    }

    fn closest_on_facet(&self, f: usize, p: &[f64; 3]) -> [f64; 3] {
        let fv = self.mesh.facet(f);
        let a = simplex_point(self.mesh, fv[0]);
        let b = simplex_point(self.mesh, fv[1]);
        if self.mesh.dim == 2 {
            closest_on_segment(p, &a, &b)
        } else {
            closest_on_triangle(p, &a, &b, &simplex_point(self.mesh, fv[2]))
        }
    }

    /// The closest point of the facets named `boundary` to `p`, read in the adjacent cell whose label is in `side_labels` (else the lower
    /// adjacent cell). Ties between facets go to the lower facet number.
    pub fn closest(&self, p: &[f64], boundary: &str, side_labels: &[i32]) -> Option<BoundaryHit> {
        let i = self.names.binary_search_by(|x| x.as_str().cmp(boundary)).ok()?;
        let p3 = to3(p);
        let (f, dist) = self.by_name[i].nearest(&p3, f64::INFINITY, |f| {
            Some(norm(&sub(&self.closest_on_facet(f as usize, &p3), &p3)))
        })?;
        let f = f as usize;
        let point = self.closest_on_facet(f, &p3);
        let cells = self.facet_cells[f];
        let cell = cells
            .iter()
            .copied()
            .filter(|&c| c != u32::MAX)
            .find(|&c| side_labels.contains(&self.mesh.cell_labels[c as usize]))
            .or_else(|| cells.iter().copied().find(|&c| c != u32::MAX))? as usize;
        let lam = straight_bary(self.mesh, cell, &point)?;
        let lam = match self.curved.get(&cell) {
            Some(ct) => curved_bary(ct, lam, &point, self.opts.newton_tol, self.opts.newton_max_iter).unwrap_or(lam),
            None => lam,
        };
        Some(BoundaryHit { facet: f, point, distance: dist, cell, coords: clamp_bary(&lam, self.mesh.dim + 1) })
    }

    /// `field` (a simplex `Nodal`, `Cell` or `PerCell` field of this mesh) at the closest boundary points of `points` (flat, `dim` each) on
    /// the facets named `boundary`.
    pub fn eval_boundary(
        &self,
        field: &SourceField,
        points: &[f64],
        boundary: &str,
        side_labels: &[i32],
    ) -> Result<Sampled, TransferError> {
        let empty = BTreeMap::new();
        let index = TransferIndex::new(SourceMesh::Simplex(self.mesh), &empty, self.opts)?;
        index.check_field(field)?;
        let d = self.mesh.dim;
        if !points.len().is_multiple_of(d) {
            return Err(TransferError::Length(format!("{} coordinates are not points of dimension {d}", points.len())));
        }
        if let Some(i) = points.iter().position(|v| !v.is_finite()) {
            return Err(TransferError::NonFinite(format!("query coordinate {i}")));
        }
        let n = points.len() / d;
        let ncomp = field.ncomp();
        let mut out = Sampled { ncomp, values: vec![0.0; n * ncomp], found: vec![false; n], snapped: 0, missing: 0 };
        for i in 0..n {
            match self.closest(&points[i * d..(i + 1) * d], boundary, side_labels) {
                Some(hit) => {
                    let loc = Located { cell: hit.cell, coords: hit.coords, snapped: hit.distance > 0.0, distance: hit.distance };
                    index.eval_located(field, &loc, &mut out.values[i * ncomp..(i + 1) * ncomp]);
                    out.found[i] = true;
                    out.snapped += loc.snapped as usize;
                }
                None => out.missing += 1,
            }
        }
        Ok(out)
    }
}

// ------------------------------------------------------------------------------------------------------------ meshes

/// Whether two exchange meshes have the same vertices, cells and domain labels (the zero-cost fast path of a transfer).
pub fn same_exchange_mesh(a: &Mesh, b: &Mesh) -> bool {
    fn same<const D: usize>(a: &CellMesh<D>, b: &CellMesh<D>) -> bool {
        a.points == b.points && a.types == b.types && a.offsets == b.offsets && a.conn == b.conn && a.label == b.label
    }
    match (a, b) {
        (Mesh::D2(x), Mesh::D2(y)) => same(x, y),
        (Mesh::D3(x), Mesh::D3(y)) => same(x, y),
        _ => false,
    }
}

/// Whether two simplex discretisations have the same vertices, cells and domain labels.
pub fn same_simplex_mesh(a: &SimplexMesh, b: &SimplexMesh) -> bool {
    a.dim == b.dim && a.points == b.points && a.cells == b.cells && a.cell_labels == b.cell_labels
}

fn restrict_cells<const D: usize>(m: &CellMesh<D>, keep: &[i32]) -> CellMesh<D> {
    let mut used = vec![false; m.n_points()];
    let kept: Vec<usize> = (0..m.n_cells()).filter(|&c| keep.contains(&m.label[c])).collect();
    for &c in &kept {
        for &v in m.cell(c) {
            used[v as usize] = true;
        }
    }
    let mut new_id = vec![u32::MAX; m.n_points()];
    let mut out = CellMesh::<D>::new();
    for (v, &u) in used.iter().enumerate() {
        if u {
            new_id[v] = out.points.len() as u32;
            out.points.push(m.points[v]);
        }
    }
    let key = |vs: &[u32]| -> [u32; 4] {
        let mut k = [u32::MAX; 4];
        for (i, &x) in vs.iter().take(4).enumerate() {
            k[i] = x;
        }
        k[..vs.len().min(4)].sort_unstable();
        k
    };
    let mut faces: std::collections::HashSet<[u32; 4]> = std::collections::HashSet::new();
    let mut vbuf: Vec<u32> = Vec::with_capacity(8);
    for &c in &kept {
        let vs = m.cell(c);
        for f in m.types[c].faces() {
            let fv: Vec<u32> = f.iter().map(|&l| vs[l as usize]).collect();
            faces.insert(key(&fv));
        }
        vbuf.clear();
        vbuf.extend(vs.iter().map(|&v| new_id[v as usize]));
        out.push_cell(m.types[c], &vbuf, m.label[c]);
    }
    for f in 0..m.n_facets() {
        let fv = m.facet(f);
        if fv.iter().all(|&v| used[v as usize]) && faces.contains(&key(fv)) {
            vbuf.clear();
            vbuf.extend(fv.iter().map(|&v| new_id[v as usize]));
            out.push_facet(&vbuf, m.facet_label[f]);
        }
    }
    let remap2 = |p: &[u32; 2]| -> Option<[u32; 2]> {
        let (a, b) = (new_id[p[0] as usize], new_id[p[1] as usize]);
        (a != u32::MAX && b != u32::MAX).then_some([a.min(b), a.max(b)])
    };
    out.set_edge_nodes(m.edge_nodes.iter().filter_map(|(pr, x)| remap2(pr).map(|q| (q, *x))));
    out.set_quad_face_nodes(m.quad_face_nodes.iter().filter_map(|(q, x)| {
        let r = q.map(|v| new_id[v as usize]);
        r.iter().all(|&v| v != u32::MAX).then_some((r, *x))
    }));
    let mut qp: Vec<[u32; 2]> = m.qp_edges.iter().filter_map(remap2).collect();
    qp.sort_unstable();
    qp.dedup();
    out.qp_edges = qp;
    out.retain_existing_edge_nodes();
    out
}

/// The sub-mesh of the cells whose domain label is in `keep` (an active-domain subset): vertices renumbered in their old order, the
/// labelled facets that are faces of a kept cell kept with their labels, curved edges and quad-face nodes carried over. Labels are not
/// renumbered, so the canonical names of the model still apply.
pub fn restrict_to_domains(mesh: &Mesh, keep: &[i32]) -> Result<Mesh, TransferError> {
    let out = match mesh {
        Mesh::D2(m) => Mesh::D2(restrict_cells(m, keep)),
        Mesh::D3(m) => Mesh::D3(restrict_cells(m, keep)),
    };
    if out.n_cells() == 0 {
        return Err(TransferError::Length(format!("no cell carries one of the active domain labels {keep:?}")));
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------------------------ tests

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::{block_split_z, rect_split_x};

    fn names2() -> BTreeMap<i32, String> {
        BTreeMap::from([(1, "Si".to_string()), (2, "SiO2".to_string())])
    }

    fn quad2(p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        1.0 + 2.0 * x - 3.0 * y + x * x + 0.5 * x * y - y * y
    }

    fn quad3(p: &[f64]) -> f64 {
        let (x, y, z) = (p[0], p[1], p[2]);
        0.5 + x - 2.0 * y + 0.3 * z + x * x - 0.7 * y * z + 0.25 * z * z + x * y
    }

    fn nodal(mesh: &SimplexMesh, order: usize, f: impl Fn(&[f64]) -> f64) -> Vec<f64> {
        mesh.dof_coordinates(order).chunks(mesh.dim).map(f).collect()
    }

    #[test]
    fn p2_quadratic_field_transfers_exactly_between_non_matching_2d_meshes() {
        let src = rect_split_x(1.0, 0.5, 0.4, 0.09);
        let tgt = rect_split_x(1.0, 0.5, 0.4, 0.047);
        assert_ne!(src.n_cells(), tgt.n_cells());
        let vals = nodal(&src, 2, quad2);
        let field = SourceField::Nodal { values: &vals, ncomp: 1, order: 2 };
        let idx = TransferIndex::new(SourceMesh::Simplex(&src), &names2(), TransferOptions::default()).unwrap();
        let s = idx.sample_cell_nodes(&field, &tgt, &names2(), 2).unwrap();
        assert_eq!(s.missing, 0);
        let nloc = tgt.nloc(2);
        let (mut all, mut interior) = (0.0f64, 0.0f64);
        for c in 0..tgt.n_cells() {
            for a in 0..nloc {
                let x = target_node(&tgt, c, a, 2);
                let e = (s.values[c * nloc + a] - quad2(&x)).abs();
                all = all.max(e);
                let on_edge = x[0] < 1e-9 || x[0] > 1.0 - 1e-9 || x[1] < 1e-9 || x[1] > 0.5 - 1e-9;
                if !on_edge {
                    interior = interior.max(e);
                }
            }
        }
        assert!(all < 1e-6, "max error {all:e}");
        assert!(interior < 1e-12, "interior error {interior:e}");
        // quadrature points too
        let q = idx.sample_quadrature(&field, &tgt, &names2(), 4).unwrap();
        let rule = quad_rule(2, 4);
        for c in (0..tgt.n_cells()).step_by(7) {
            for (k, l) in rule.lam.iter().enumerate() {
                let x = target_point(&tgt, c, l);
                assert!((q.values[c * rule.w.len() + k] - quad2(&x)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn p2_quadratic_field_transfers_exactly_between_non_matching_3d_meshes() {
        let src = block_split_z(1.0, 0.6, 0.8, 0.4, 0.3, 2);
        let tgt = block_split_z(1.0, 0.6, 0.8, 0.4, 0.21, 3);
        let names = BTreeMap::from([(1, "lower".to_string()), (2, "upper".to_string())]);
        let vals = nodal(&src, 2, quad3);
        let field = SourceField::Nodal { values: &vals, ncomp: 1, order: 2 };
        let idx = TransferIndex::new(SourceMesh::Simplex(&src), &names, TransferOptions::default()).unwrap();
        let s = idx.sample_dofs(&field, &tgt, &names, 2).unwrap();
        assert_eq!(s.missing, 0);
        let want = nodal(&tgt, 2, quad3);
        let err = s.values.iter().zip(&want).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(err < 1e-12, "3D P2 transfer error {err:e}");
    }

    #[test]
    fn a_12_to_1_jump_at_a_si_sio2_interface_is_not_smeared() {
        let src = rect_split_x(1.0, 0.5, 0.4, 0.08);
        let tgt = rect_split_x(1.0, 0.5, 0.4, 0.05);
        // per-cell P1 values: 12 (1 + x) in Si, (1 + x) in SiO2, discontinuous at x = 0.4
        let nloc = src.nloc(1);
        let mut vals = vec![0.0; src.n_cells() * nloc];
        for c in 0..src.n_cells() {
            let k = if src.cell_labels[c] == 1 { 12.0 } else { 1.0 };
            for (a, &v) in src.cell(c).iter().enumerate() {
                vals[c * nloc + a] = k * (1.0 + src.point(v as usize)[0]);
            }
        }
        let field = SourceField::Cell { values: &vals, order: 1 };
        let idx = TransferIndex::new(SourceMesh::Simplex(&src), &names2(), TransferOptions::default()).unwrap();
        let s = idx.sample_cell_nodes(&field, &tgt, &names2(), 2).unwrap();
        assert_eq!(s.missing, 0);
        let tl = tgt.nloc(2);
        let mut on_interface = 0;
        for c in 0..tgt.n_cells() {
            let k = if tgt.cell_labels[c] == 1 { 12.0 } else { 1.0 };
            for a in 0..tl {
                let x = target_node(&tgt, c, a, 2);
                if (x[0] - 0.4).abs() < 1e-12 {
                    on_interface += 1;
                }
                let want = k * (1.0 + x[0]);
                assert!((s.values[c * tl + a] - want).abs() < 1e-12, "cell {c} node {a}: {} vs {want}", s.values[c * tl + a]);
            }
        }
        assert!(on_interface > 0);
        // a projection keeps the jump too
        let tidx = TransferIndex::new(SourceMesh::Simplex(&tgt), &names2(), TransferOptions::default()).unwrap();
        let p = project_scalar_conservative(&idx, &field, &tidx, 1, MassKind::Consistent).unwrap();
        for c in 0..tgt.n_cells() {
            let k = if tgt.cell_labels[c] == 1 { 12.0 } else { 1.0 };
            for (a, &v) in tgt.cell(c).iter().enumerate() {
                let want = k * (1.0 + tgt.point(v as usize)[0]);
                assert!((p.cell_values[c * 3 + a] - want).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn the_conservative_projection_keeps_every_domain_integral() {
        let coarse = rect_split_x(1.0, 0.5, 0.4, 0.1);
        let fine = rect_split_x(1.0, 0.5, 0.4, 0.043);
        let f = |p: &[f64]| (3.0 * p[0]).sin() * (1.0 + p[1] * p[1]) + if p[0] < 0.4 { 11.0 } else { 0.0 };
        for (src, tgt) in [(&coarse, &fine), (&fine, &coarse)] {
            // a per-cell P2 field (discontinuous at the interface)
            let nloc = src.nloc(2);
            let mut vals = vec![0.0; src.n_cells() * nloc];
            for c in 0..src.n_cells() {
                for a in 0..nloc {
                    let mut x = target_node(src, c, a, 2);
                    if src.cell_labels[c] == 2 {
                        x[0] = x[0].max(0.4 + 1e-15);
                    } else {
                        x[0] = x[0].min(0.4 - 1e-15);
                    }
                    vals[c * nloc + a] = f(&x);
                }
            }
            let field = SourceField::Cell { values: &vals, order: 2 };
            let si = TransferIndex::new(SourceMesh::Simplex(src), &names2(), TransferOptions::default()).unwrap();
            let ti = TransferIndex::new(SourceMesh::Simplex(tgt), &names2(), TransferOptions::default()).unwrap();
            for order in [1, 2] {
                for mass in [MassKind::Consistent, MassKind::Lumped] {
                    let p = project_scalar_conservative(&si, &field, &ti, order, mass).unwrap();
                    assert_eq!(p.unmatched_source_integral, 0.0);
                    for (name, s) in &p.source_integrals {
                        let t = p.target_integrals[name];
                        assert!((t - s).abs() <= 1e-12 * s.abs().max(1.0), "{name} order {order} {mass:?}: {s} vs {t}");
                    }
                }
            }
        }
    }

    #[test]
    fn transfers_are_bitwise_identical_serial_and_parallel() {
        let src = rect_split_x(1.0, 0.5, 0.4, 0.07);
        let tgt = rect_split_x(1.0, 0.5, 0.4, 0.031);
        let vals = nodal(&src, 2, |p| (5.0 * p[0]).cos() + p[1]);
        let field = SourceField::Nodal { values: &vals, ncomp: 1, order: 2 };
        let run = |parallel: bool, chunk: usize| {
            let o = TransferOptions { parallel, chunk, ..TransferOptions::default() };
            let si = TransferIndex::new(SourceMesh::Simplex(&src), &names2(), o).unwrap();
            let ti = TransferIndex::new(SourceMesh::Simplex(&tgt), &names2(), o).unwrap();
            let s = si.sample_cell_nodes(&field, &tgt, &names2(), 2).unwrap();
            let p = project_scalar_conservative(&si, &field, &ti, 2, MassKind::Consistent).unwrap();
            (s.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), p.cell_values.iter().map(|v| v.to_bits()).collect::<Vec<_>>())
        };
        let a = run(false, 256);
        assert_eq!(a, run(true, 256));
        assert_eq!(a.0, run(true, 13).0);
    }

    #[test]
    fn boundary_transfer_reads_the_named_facets_at_the_closest_point() {
        let m = rect_split_x(1.0, 0.5, 0.4, 0.06);
        let vals = nodal(&m, 2, quad2);
        let field = SourceField::Nodal { values: &vals, ncomp: 1, order: 2 };
        let bn = BTreeMap::from([(1, "left".to_string()), (2, "bottom".to_string()), (3, "top".to_string()), (4, "right".to_string())]);
        let bi = BoundaryIndex::new(&m, &bn, TransferOptions::default()).unwrap();
        let pts = [1.3, 0.1, 1.05, 0.37, 0.7, 0.2];
        let s = bi.eval_boundary(&field, &pts, "right", &[]).unwrap();
        assert_eq!(s.missing, 0);
        for i in 0..3 {
            let want = quad2(&[1.0, pts[2 * i + 1]]);
            assert!((s.values[i] - want).abs() < 1e-12, "{} vs {want}", s.values[i]);
        }
        assert!(bi.closest(&[0.5, 0.2], "nowhere", &[]).is_none());
    }

    #[test]
    fn a_curved_tetrahedron_is_inverted_by_newton() {
        let pts = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let mut m = SimplexMesh::new(3, pts, vec![0, 1, 2, 3], vec![1], vec![], vec![]);
        m.set_edge_nodes(&[([1, 2], [0.6, 0.6, 0.0]), ([0, 1], [0.5, -0.05, 0.02])]);
        assert!(m.cell_is_curved(0));
        let idx = TransferIndex::new(SourceMesh::Simplex(&m), &BTreeMap::new(), TransferOptions::default()).unwrap();
        let ct = CurvedTet::new(m.cell_geometry_nodes(0));
        for lam in [[0.1, 0.2, 0.3, 0.4], [0.05, 0.6, 0.3, 0.05], [0.25, 0.25, 0.25, 0.25]] {
            let x = ct.position(&lam);
            let loc = idx.locate(&x, None).unwrap();
            assert!(!loc.snapped);
            for k in 0..4 {
                assert!((loc.coords[k] - lam[k]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn mixed_quads_are_read_through_their_bilinear_map() {
        let mut cm = CellMesh::<2>::new();
        let n = 4usize;
        for j in 0..=n {
            for i in 0..=n {
                let (x, y) = (i as f64 / n as f64, j as f64 / n as f64);
                cm.points.push([x + 0.05 * y * x, y]);
            }
        }
        let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
        for j in 0..n {
            for i in 0..n {
                let lab = if i < 2 { 1 } else { 2 };
                cm.push_cell(CellType::Quad, &[id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)], lab);
            }
        }
        let mixed = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        let dm = MixedLagrangeDofMap::new(&mixed, LagrangeSpace::Order1);
        // a field that is linear in x and y is in Q1 on every (bilinear) quad
        let mut vals = vec![0.0; mixed.n_points()];
        for v in 0..mixed.n_points() {
            let p = mixed.point(v);
            vals[v] = 2.0 + p[0] - 3.0 * p[1];
        }
        let field = SourceField::MixedNodal { dofs: &dm, space: LagrangeSpace::Order1, values: &vals, ncomp: 1 };
        let idx = TransferIndex::new(SourceMesh::Mixed(&mixed), &names2(), TransferOptions::default()).unwrap();
        for p in [[0.13, 0.71], [0.62, 0.33], [0.9, 0.95]] {
            let dom = if idx.locate(&p, None).map(|l| mixed.cell_labels[l.cell]) == Some(1) { "Si" } else { "SiO2" };
            let v = idx.eval_scalar(&field, &p, Some(dom)).unwrap().unwrap();
            assert!((v - (2.0 + p[0] - 3.0 * p[1])).abs() < 1e-12);
        }
    }

    #[test]
    fn an_active_domain_subset_keeps_its_cells_and_facets() {
        let m = rect_split_x(1.0, 0.5, 0.4, 0.08);
        let mesh = Mesh::from_solver_mesh(&m).unwrap();
        let sub = restrict_to_domains(&mesh, &[1]).unwrap();
        let s = sub.discretise().unwrap();
        assert_eq!(s.n_cells(), m.cell_labels.iter().filter(|&&l| l == 1).count());
        assert!(s.cell_labels.iter().all(|&l| l == 1));
        s.check().unwrap();
        assert!(s.points.chunks(2).all(|p| p[0] <= 0.4 + 1e-12));
        assert!(restrict_to_domains(&mesh, &[7]).is_err());
        assert!(same_exchange_mesh(&mesh, &mesh.clone()));
        assert!(!same_exchange_mesh(&mesh, &sub));
    }
}
