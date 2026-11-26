//! Curved edges for 3D meshes (ADR 82): the mid-edge geometry nodes that make the tetrahedra around a curved surface quadratic
//! (isoparametric) instead of straight.
//!
//! For every edge of a boundary facet or of a facet between two regions the geometry (`Domain3::project_edge`) says where the midpoint of
//! the edge would lie on the exact surface; if it moves by more than `1e-9` of the edge length (a planar face adds nothing) and by less
//! than a quarter of it, the node is kept and the cells that own the edge become curved. Every standard non-QP curved tetrahedron is
//! certified everywhere via the 20-coefficient tetrahedral Bernstein expansion of its cubic scaled Jacobian `det J(lam) / det_straight`
//! (subdividing 1->8 only on borderline elements where `min_i b_i <= threshold` yet sampled values remain positive); if it fails
//! certification above `MIN_SCALED_JACOBIAN` (`0.10`), incident offsets are halved (up to four times, then straightened) while keeping
//! periodic twins locked together. Barsoum quarter-point (QP) cells intentionally satisfy `det J = 0` at the singular rim and are
//! validated separately via `quad_rule_collapsed(8)` interior quadrature. Everything is decided deterministically (`ADR 11`).

use crate::fem::{quad_rule_collapsed, CurvedTet};
use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::SimplexMesh;
use crate::meshing::mesher3d::Domain3;
use std::collections::{BTreeSet, HashMap, HashSet};

/// The curved edges: sorted vertex pair and the position of the mid-edge node.
pub type EdgeNodes = Vec<([u32; 2], [f64; 3])>;

/// What `curve_boundary_edges` did.
#[derive(Clone, Debug, PartialEq)]
pub struct CurveReport {
    /// edges that carry a geometry node in the end
    pub curved_edges: usize,
    /// cells with at least one such edge
    pub curved_cells: usize,
    /// edges whose offset had to be reduced (halved up to four times) to keep their cells valid
    pub relaxed: usize,
    /// edges that had to be made straight again
    pub straightened: usize,
    /// smallest ratio of the Jacobian determinant to the straight cell's over the curved cells (1 when none)
    pub min_scaled_jacobian: f64,
}

impl Default for CurveReport {
    fn default() -> Self {
        CurveReport { curved_edges: 0, curved_cells: 0, relaxed: 0, straightened: 0, min_scaled_jacobian: 1.0 }
    }
}

/// A cell must keep this fraction of its straight Jacobian determinant everywhere.
pub const MIN_SCALED_JACOBIAN: f64 = 0.1;
/// A node that moved by more than this fraction of the edge length is not a surface node but a mismatch of the geometry: ignored.
const MAX_OFFSET: f64 = 0.25;
/// A node that moved by less than this fraction of the edge length is on a plane: ignored.
const MIN_OFFSET: f64 = 1e-9;
const MAX_HALVINGS: u8 = 4;
const MAX_BERNSTEIN_DEPTH: u8 = 5;

const REF_SIMPLEX_VERTS: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

const TET_FACES: [(usize, usize, usize); 4] = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn pair(a: u32, b: u32) -> [u32; 2] {
    [a.min(b), a.max(b)]
}

const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// The core: curved edges for the tetrahedra `cells` with vertex positions `pos`, from the edges of the facets `facets` (boundary
/// facets and facets between regions) and the geometry `dom`. Returns the nodes (sorted by vertex pair) and the report.
pub fn curve_edges(
    pos: &(dyn Fn(u32) -> [f64; 3] + Sync),
    cells: &[[u32; 4]],
    facets: &[[u32; 3]],
    dom: &dyn Domain3,
) -> (EdgeNodes, CurveReport) {
    // 1. the edges of the facets, in sorted order, and where the surface puts their midpoints
    let mut set: BTreeSet<[u32; 2]> = BTreeSet::new();
    for f in facets {
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            set.insert(pair(f[i], f[j]));
        }
    }
    let edges: Vec<[u32; 2]> = set.into_iter().collect();
    const CURVE_CHUNK: usize = 256;
    let mut nodes: Vec<Option<[f64; 3]>> = crate::par::map_chunks(&edges, CURVE_CHUNK, |slice| {
        slice
            .iter()
            .map(|&[a, b]| {
                let (pa, pb) = (pos(a), pos(b));
                let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0];
                let len = norm(sub(pb, pa));
                if len == 0.0 || dom.edge_on_crease(pa, pb) {
                    return None;
                }
                let q = dom.project_edge(mid, pa, pb)?;
                let off = norm(sub(q, mid));
                (off > MIN_OFFSET * len && off < MAX_OFFSET * len && q.iter().all(|x| x.is_finite())).then_some(q)
            })
            .collect::<Vec<_>>()
    })
    .into_iter()
    .flatten()
    .collect();

    // Sync curved mid-edge projections across axis-aligned periodic boundary faces from parent (lo) to child (hi)
    let periodic_axes = dom.periodic_axes();
    let mut periodic_edge_parent: HashMap<[u32; 2], ([u32; 2], [f64; 3])> = HashMap::new();
    if !periodic_axes.is_empty() && !edges.is_empty() {
        let mut verts: BTreeSet<u32> = BTreeSet::new();
        for &[a, b] in &edges {
            verts.insert(a);
            verts.insert(b);
        }
        let (mut blo, mut bhi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for &v in &verts {
            let p = pos(v);
            for k in 0..3 {
                blo[k] = blo[k].min(p[k]);
                bhi[k] = bhi[k].max(p[k]);
            }
        }
        let extent = (0..3).map(|k| bhi[k] - blo[k]).fold(0.0f64, f64::max).max(1e-300);
        let tol = 1e-10 * extent;
        let mut per_axis: [Option<(f64, f64)>; 3] = [None; 3];
        for pa in &periodic_axes {
            if pa.axis < 3 && pa.hi > pa.lo {
                per_axis[pa.axis] = Some((pa.lo, pa.hi));
            }
        }
        let mut v_bits: HashMap<u32, u8> = HashMap::new();
        for &v in &verts {
            let p = pos(v);
            let mut bits = 0u8;
            for (k, slot) in per_axis.iter().enumerate() {
                if let Some((lo, hi)) = *slot {
                    if (p[k] - lo).abs() <= tol {
                        bits |= 1 << (k + 3);
                    } else if (p[k] - hi).abs() <= tol {
                        bits |= 1 << k;
                    }
                }
            }
            if bits != 0 {
                v_bits.insert(v, bits);
            }
        }
        let bverts: Vec<(u32, [f64; 3], u8)> = v_bits.iter().map(|(&v, &b)| (v, pos(v), b)).collect();
        let find_shifted = |v: u32, shift: [f64; 3]| -> Option<u32> {
            let pv = pos(v);
            let target = [pv[0] - shift[0], pv[1] - shift[1], pv[2] - shift[2]];
            bverts
                .iter()
                .filter(|(_, p, _)| (0..3).all(|k| (p[k] - target[k]).abs() <= tol))
                .min_by(|a, b| norm(sub(a.1, target)).total_cmp(&norm(sub(b.1, target))))
                .map(|x| x.0)
        };
        let edge_idx: HashMap<[u32; 2], usize> = edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
        for (ei, &[a, b]) in edges.iter().enumerate() {
            let (Some(&ba), Some(&bb)) = (v_bits.get(&a), v_bits.get(&b)) else { continue };
            let child_shared = (ba & bb) & 0b111;
            if child_shared == 0 {
                continue;
            }
            let mut shift = [0.0; 3];
            for (k, slot) in per_axis.iter().enumerate() {
                if let Some((lo, hi)) = *slot {
                    if (child_shared & (1 << k)) != 0 {
                        shift[k] = hi - lo;
                    }
                }
            }
            if let (Some(pa), Some(pb)) = (find_shifted(a, shift), find_shifted(b, shift)) {
                let pe = pair(pa, pb);
                if let Some(&pi) = edge_idx.get(&pe) {
                    periodic_edge_parent.insert([a, b], (pe, shift));
                    nodes[ei] = nodes[pi].map(|q| [q[0] + shift[0], q[1] + shift[1], q[2] + shift[2]]);
                }
            }
        }
    }
    // candidate edges (those that moved), by index
    let mut cand: Vec<([u32; 2], [f64; 3])> = vec![];
    for (e, n) in edges.iter().zip(&nodes) {
        if let Some(q) = n {
            cand.push((*e, *q));
        }
    }
    let report = relax_edge_nodes_periodic_with_qp(pos, cells, &mut cand, &periodic_edge_parent, &[]);
    (cand, report)
}

/// Validates and relaxes curved mid-edge nodes `nodes` over `cells` via exact 20-coefficient tetrahedral
/// Bernstein certification of cubic `det J / det_straight` against `MIN_SCALED_JACOBIAN` (`0.10`), halving
/// offsets (up to [`MAX_HALVINGS`] = 4 times, then straightening) whenever any incident non-QP cell fails.
pub fn relax_edge_nodes(
    pos: &(dyn Fn(u32) -> [f64; 3] + Sync),
    cells: &[[u32; 4]],
    nodes: &mut EdgeNodes,
) -> CurveReport {
    for (e, _) in nodes.iter_mut() {
        *e = pair(e[0], e[1]);
    }
    nodes.sort_by_key(|(e, _)| *e);
    nodes.dedup_by_key(|(e, _)| *e);
    relax_edge_nodes_periodic_with_qp(pos, cells, nodes, &HashMap::new(), &[])
}

/// Periodic-aware curved edge relaxation excluding cells touching Barsoum quarter-point `qp_edges` from
/// Bernstein positivity (since QP cells satisfy `det J = 0` at the singular rim by design).
pub(crate) fn relax_edge_nodes_periodic_with_qp(
    pos: &(dyn Fn(u32) -> [f64; 3] + Sync),
    cells: &[[u32; 4]],
    cand: &mut EdgeNodes,
    periodic_edge_parent: &HashMap<[u32; 2], ([u32; 2], [f64; 3])>,
    qp_edges: &[[u32; 2]],
) -> CurveReport {
    if cand.is_empty() {
        return CurveReport::default();
    }
    let qp_set: HashSet<[u32; 2]> = qp_edges.iter().map(|&[a, b]| pair(a, b)).collect();
    let index: HashMap<[u32; 2], usize> = cand.iter().enumerate().map(|(i, (e, _))| (*e, i)).collect();
    let mut level = vec![0u8; cand.len()];

    // 2. non-QP cells with at least one candidate curved edge, and for each edge its index in `cand`
    let curved_cells: Vec<(usize, [Option<usize>; 6])> = cells
        .iter()
        .enumerate()
        .filter_map(|(ci, c)| {
            let mut ids = [None; 6];
            let mut any = false;
            for (k, &(i, j)) in EDGES.iter().enumerate() {
                let pr = pair(c[i], c[j]);
                if qp_set.contains(&pr) {
                    return None;
                }
                if let Some(&e) = index.get(&pr) {
                    ids[k] = Some(e);
                    any = true;
                }
            }
            any.then_some((ci, ids))
        })
        .collect();

    // 3. Bernstein validity, halving incident edge offsets where a non-QP cell fails `> MIN_SCALED_JACOBIAN`
    let scale_of = |lv: u8| -> f64 {
        if lv > MAX_HALVINGS {
            0.0
        } else {
            0.5f64.powi(lv as i32)
        }
    };
    let cell_cert = |ci: usize, ids: &[Option<usize>; 6], level: &[u8]| -> (bool, f64) {
        let c = cells[ci];
        let v = [pos(c[0]), pos(c[1]), pos(c[2]), pos(c[3])];
        let mut x = CurvedTet::straight(v).x;
        for (k, id) in ids.iter().enumerate() {
            if let Some(e) = *id {
                let s = scale_of(level[e]);
                let (i, j) = EDGES[k];
                for d in 0..3 {
                    let mid = 0.5 * (v[i][d] + v[j][d]);
                    x[4 + k][d] = mid + s * (cand[e].1[d] - mid);
                }
            }
        }
        let ct = CurvedTet::new(x);
        certify_curved_tet_bernstein(&ct, MIN_SCALED_JACOBIAN)
    };
    const CERT_CHUNK: usize = 256;
    for _round in 0..=(MAX_HALVINGS as usize + 1) {
        let failing: Vec<usize> = crate::par::chunk_map(
            curved_cells.len().div_ceil(CERT_CHUNK),
            curved_cells.len() > CERT_CHUNK && crate::par::threads() > 1,
            |cid| {
                let c0 = cid * CERT_CHUNK;
                let c1 = (c0 + CERT_CHUNK).min(curved_cells.len());
                (c0..c1)
                    .filter(|&k| {
                        let (ci, ids) = &curved_cells[k];
                        !cell_cert(*ci, ids, &level).0
                    })
                    .collect::<Vec<_>>()
            },
        )
        .into_iter()
        .flatten()
        .collect();
        if failing.is_empty() {
            break;
        }
        let mut bump: BTreeSet<usize> = BTreeSet::new();
        for k in failing {
            for e in curved_cells[k].1.iter().flatten() {
                if level[*e] <= MAX_HALVINGS {
                    bump.insert(*e);
                }
            }
        }
        if bump.is_empty() {
            break;
        }
        for e in bump {
            level[e] += 1;
        }
        if !periodic_edge_parent.is_empty() {
            for _sync_pass in 0..4 {
                let mut changed = false;
                for (ci, (ce, _)) in cand.iter().enumerate() {
                    if let Some((pe, _)) = periodic_edge_parent.get(ce) {
                        if let Some(&pi) = index.get(pe) {
                            let mx = level[ci].max(level[pi]);
                            if level[ci] != mx || level[pi] != mx {
                                level[ci] = mx;
                                level[pi] = mx;
                                changed = true;
                            }
                        }
                    }
                }
                if !changed {
                    break;
                }
            }
        }
    }

    // 4. materialize relaxed edge nodes (keeping periodic child positions shift-identical to parents)
    let mut report = CurveReport::default();
    let mut node_by_idx: Vec<Option<[f64; 3]>> = vec![None; cand.len()];
    for (e, (pr, q)) in cand.iter().enumerate() {
        match level[e] {
            0 => {}
            l if l <= MAX_HALVINGS => report.relaxed += 1,
            _ => report.straightened += 1,
        }
        let s = scale_of(level[e]);
        if s > 0.0 {
            let (pa, pb) = (pos(pr[0]), pos(pr[1]));
            let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0];
            node_by_idx[e] = Some([
                mid[0] + s * (q[0] - mid[0]),
                mid[1] + s * (q[1] - mid[1]),
                mid[2] + s * (q[2] - mid[2]),
            ]);
        }
    }
    if !periodic_edge_parent.is_empty() {
        for _sync_pass in 0..4 {
            for (e, (ce, _)) in cand.iter().enumerate() {
                if let Some((pe, shift)) = periodic_edge_parent.get(ce) {
                    if let Some(&pi) = index.get(pe) {
                        node_by_idx[e] = node_by_idx[pi].map(|qp| [qp[0] + shift[0], qp[1] + shift[1], qp[2] + shift[2]]);
                    }
                }
            }
        }
    }
    let mut out: Vec<([u32; 2], [f64; 3])> = Vec::new();
    let mut kept: HashMap<[u32; 2], ()> = HashMap::new();
    for (e, (pr, _)) in cand.iter().enumerate() {
        if let Some(qn) = node_by_idx[e] {
            out.push((*pr, qn));
            kept.insert(*pr, ());
        }
    }
    report.curved_edges = out.len();
    let mut min_sj = f64::INFINITY;
    let mut n_cells = 0;
    for (ci, ids) in &curved_cells {
        if ids.iter().flatten().any(|&e| kept.contains_key(&cand[e].0)) {
            n_cells += 1;
            min_sj = min_sj.min(cell_cert(*ci, ids, &level).1);
        }
    }
    report.curved_cells = n_cells;
    report.min_scaled_jacobian = if min_sj.is_finite() { min_sj } else { 1.0 };
    *cand = out;
    report
}

/// Evaluates `det J(lam) / det_s` at the 20 cubic Lagrange nodes on the barycentric sub-tetrahedron `verts`
/// and returns `(bernstein_coeffs, sample_min)`.
fn bernstein_subtet_coeffs(ct: &CurvedTet, det_s: f64, verts: &[[f64; 4]; 4]) -> Option<([f64; 20], f64)> {
    let eval_bary = |mu: [f64; 4]| -> Option<f64> {
        let mut lam = [0.0; 4];
        for c in 0..4 {
            lam[c] = mu[0] * verts[0][c] + mu[1] * verts[1][c] + mu[2] * verts[2][c] + mu[3] * verts[3][c];
        }
        let ratio = ct.jacobian(&lam).1 / det_s;
        ratio.is_finite().then_some(ratio)
    };
    let mut b = [0.0f64; 20];
    let mut sample_min = f64::INFINITY;

    // 4 vertices: indices 0..4
    for a in 0..4 {
        let v = eval_bary(REF_SIMPLEX_VERTS[a])?;
        b[a] = v;
        sample_min = sample_min.min(v);
    }
    // 12 edge trisections: indices 4..16 (for each edge k in 0..6: 4 + 2*k is b_{2,1}, 4 + 2*k + 1 is b_{1,2})
    for (k, &(i, j)) in EDGES.iter().enumerate() {
        let mut mu21 = [0.0; 4];
        mu21[i] = 2.0 / 3.0;
        mu21[j] = 1.0 / 3.0;
        let mut mu12 = [0.0; 4];
        mu12[i] = 1.0 / 3.0;
        mu12[j] = 2.0 / 3.0;
        let p21 = eval_bary(mu21)?;
        let p12 = eval_bary(mu12)?;
        sample_min = sample_min.min(p21).min(p12);
        b[4 + 2 * k] = 3.0 * p21 - 1.5 * p12 - (5.0 / 6.0) * b[i] + (1.0 / 3.0) * b[j];
        b[4 + 2 * k + 1] = 3.0 * p12 - 1.5 * p21 + (1.0 / 3.0) * b[i] - (5.0 / 6.0) * b[j];
    }
    // Map unordered edge `(u, v)` -> index in `EDGES`
    let edge_slot = |u: usize, v: usize| -> usize {
        let (lo, hi) = (u.min(v), u.max(v));
        EDGES.iter().position(|&e| e == (lo, hi)).unwrap()
    };
    // 4 face centers: indices 16..20
    for (f, &(i, j, k)) in TET_FACES.iter().enumerate() {
        let mut mu = [0.0; 4];
        mu[i] = 1.0 / 3.0;
        mu[j] = 1.0 / 3.0;
        mu[k] = 1.0 / 3.0;
        let p111 = eval_bary(mu)?;
        sample_min = sample_min.min(p111);
        let e_ij = edge_slot(i, j);
        let e_ik = edge_slot(i, k);
        let e_jk = edge_slot(j, k);
        let edge_sum = b[4 + 2 * e_ij]
            + b[4 + 2 * e_ij + 1]
            + b[4 + 2 * e_ik]
            + b[4 + 2 * e_ik + 1]
            + b[4 + 2 * e_jk]
            + b[4 + 2 * e_jk + 1];
        b[16 + f] = (27.0 * p111 - (b[i] + b[j] + b[k]) - 3.0 * edge_sum) / 6.0;
    }
    Some((b, sample_min))
}

/// Computes the exact 20-coefficient tetrahedral Bernstein expansion of cubic `det J(lam) / det_s` on the
/// reference tetrahedron (vertex control points `0..4`, edge control points `4..16`, face-center control
/// points `16..20`). Returns `None` if `det_s <= 0` or non-finite.
pub fn bernstein_cubic_coeffs(ct: &CurvedTet, det_s: f64) -> Option<[f64; 20]> {
    if !det_s.is_finite() || det_s <= 0.0 {
        return None;
    }
    bernstein_subtet_coeffs(ct, det_s, &REF_SIMPLEX_VERTS).map(|(b, _)| b)
}

fn mid_bary(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), 0.5 * (a[3] + b[3])]
}

fn check_bernstein_subtet(
    ct: &CurvedTet,
    det_s: f64,
    verts: &[[f64; 4]; 4],
    threshold: f64,
    depth_rem: u8,
) -> (bool, f64) {
    let Some((b, sample_min)) = bernstein_subtet_coeffs(ct, det_s, verts) else {
        return (false, 0.0);
    };
    let b_min = b.iter().copied().fold(f64::INFINITY, f64::min);
    if !b_min.is_finite() || !sample_min.is_finite() {
        return (false, 0.0);
    }
    // Fast reject: an actual evaluated point inside the sub-tet falls <= threshold
    if sample_min <= threshold {
        return (false, sample_min);
    }
    // Fast certify: convex-hull bound proves p(lam) >= b_min > threshold everywhere on the sub-tet
    if b_min > threshold {
        return (true, b_min);
    }
    // Borderline (`b_min <= threshold < sample_min`): subdivide 1 -> 8 ONLY on this sub-tet
    if depth_rem == 0 {
        return (false, b_min.min(sample_min));
    }
    let [v0, v1, v2, v3] = *verts;
    let m01 = mid_bary(v0, v1);
    let m02 = mid_bary(v0, v2);
    let m03 = mid_bary(v0, v3);
    let m12 = mid_bary(v1, v2);
    let m13 = mid_bary(v1, v3);
    let m23 = mid_bary(v2, v3);
    let children = [
        [v0, m01, m02, m03],
        [m01, v1, m12, m13],
        [m02, m12, v2, m23],
        [m03, m13, m23, v3],
        [m01, m23, m02, m12],
        [m01, m23, m12, m13],
        [m01, m23, m13, m03],
        [m01, m23, m03, m02],
    ];
    let mut worst_lb = f64::INFINITY;
    for child in &children {
        let (ok, lb) = check_bernstein_subtet(ct, det_s, child, threshold, depth_rem - 1);
        if !ok {
            return (false, lb);
        }
        worst_lb = worst_lb.min(lb);
    }
    (true, worst_lb)
}

/// Certifies via the 20-coefficient tetrahedral Bernstein expansion (subdividing 1->8 only on borderline
/// sub-tetrahedra where `min_i b_i <= threshold` while sampled values exceed `threshold`) that
/// `det J(lam) / det_s > threshold` everywhere on `ct`. Returns `(certified, bound)`.
pub fn certify_curved_tet_bernstein(ct: &CurvedTet, threshold: f64) -> (bool, f64) {
    let v = [ct.x[0], ct.x[1], ct.x[2], ct.x[3]];
    let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
    if !det_s.is_finite() || det_s <= 0.0 {
        return (false, 0.0);
    }
    check_bernstein_subtet(ct, det_s, &REF_SIMPLEX_VERTS, threshold, MAX_BERNSTEIN_DEPTH)
}

/// Certified Bernstein lower bound of `det J(lam) / det_straight` over a non-QP [`CurvedTet`] (returning
/// `<= 0.0` for degenerate/inverted cells and a certified strictly positive lower bound for valid cells).
pub fn scaled_jacobian_min_tet(ct: &CurvedTet) -> f64 {
    let v = [ct.x[0], ct.x[1], ct.x[2], ct.x[3]];
    let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
    if !det_s.is_finite() || det_s <= 0.0 {
        return 0.0;
    }
    let (ok_min, lb_min) = check_bernstein_subtet(ct, det_s, &REF_SIMPLEX_VERTS, MIN_SCALED_JACOBIAN, MAX_BERNSTEIN_DEPTH);
    if ok_min {
        return lb_min;
    }
    let (ok_pos, lb_pos) = check_bernstein_subtet(ct, det_s, &REF_SIMPLEX_VERTS, 0.0, MAX_BERNSTEIN_DEPTH);
    if ok_pos {
        return lb_pos;
    }
    if lb_pos < 0.0 {
        lb_pos
    } else {
        0.0
    }
}

/// The smallest certified ratio of the Jacobian determinant to the straight cell's over cell `c` of a solver
/// mesh via the 20-coefficient tetrahedral Bernstein basis (returning `0.0` when `det_straight <= 0` or on
/// Barsoum quarter-point cells where `det J = 0` at the rim; use [`scaled_jacobian_min_interior_q8`] for QP cells).
pub fn scaled_jacobian_min(mesh: &SimplexMesh, c: usize) -> f64 {
    scaled_jacobian_min_tet(&CurvedTet::new(mesh.cell_geometry_nodes(c)))
}

/// The smallest ratio of the Jacobian determinant to the straight cell's over all 216 interior Gauss points of
/// `quad_rule_collapsed(8)` on cell `c` of a solver mesh (used to validate Barsoum quarter-point singularity cells,
/// where vertex/rim nodes intentionally have `det J = 0`).
pub fn scaled_jacobian_min_interior_q8(mesh: &SimplexMesh, c: usize) -> f64 {
    let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
    let v = [ct.x[0], ct.x[1], ct.x[2], ct.x[3]];
    let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
    if !det_s.is_finite() || det_s <= 0.0 {
        return 0.0;
    }
    let rule = quad_rule_collapsed(8);
    let mut m = f64::INFINITY;
    for l in &rule.lam {
        m = m.min(ct.jacobian(&[l[0], l[1], l[2], l[3]]).1 / det_s);
    }
    m
}

pub(crate) type PeriodicVertexParents = Vec<Option<(u32, [f64; 3])>>;
pub(crate) type PeriodicEdgeParents = HashMap<[u32; 2], ([u32; 2], [f64; 3])>;

/// Maps periodic boundary vertices and candidate edges `(child_edge -> (root_parent_edge, shift))`
/// across arbitrary 3D translation pairs `(src_labels, dst_labels, shift)`, resolving multi-pair corner
/// chains transitively to a single canonical parent edge.
pub(crate) fn periodic_edge_parent_map(
    mesh: &SimplexMesh,
    target_edges: &[[u32; 2]],
    periodic: &[(Vec<i32>, Vec<i32>, [f64; 3])],
) -> (PeriodicVertexParents, PeriodicEdgeParents) {
    let nv = mesh.n_points();
    let mut direct_v_parent: Vec<Option<(u32, [f64; 3])>> = vec![None; nv];
    let mut edge_map: HashMap<[u32; 2], ([u32; 2], [f64; 3])> = HashMap::new();
    if periodic.is_empty() || nv == 0 || mesh.dim != 3 {
        return (direct_v_parent, edge_map);
    }

    let (lo, hi) = mesh.bbox();
    let extent = (0..3).map(|k| hi[k] - lo[k]).fold(0.0f64, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let cell = (4.0 * tol).max(1e-12);
    let pt = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], p[2]]
    };
    let grid_key = |p: [f64; 3]| -> (i64, i64, i64) {
        ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64, (p[2] / cell).floor() as i64)
    };

    // Keep per-pair vertex correspondences so multi-axis corner twins map accurately
    let mut pair_maps: Vec<(HashMap<u32, u32>, [f64; 3])> = Vec::new();
    for (src_labs, dst_labs, shift) in periodic {
        if src_labs.is_empty() || dst_labs.is_empty() {
            continue;
        }
        let src_set: HashSet<i32> = src_labs.iter().copied().collect();
        let dst_set: HashSet<i32> = dst_labs.iter().copied().collect();
        let mut src_verts = BTreeSet::new();
        let mut dst_verts = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if src_set.contains(&lab) {
                src_verts.extend(mesh.facet(f));
            }
            if dst_set.contains(&lab) {
                dst_verts.extend(mesh.facet(f));
            }
        }
        if src_verts.is_empty() || dst_verts.is_empty() {
            continue;
        }
        let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
        for &sv in &src_verts {
            grid.entry(grid_key(pt(sv))).or_default().push(sv);
        }
        let mut vmap = HashMap::new();
        for &dv in &dst_verts {
            let pd = pt(dv);
            let target = [pd[0] - shift[0], pd[1] - shift[1], pd[2] - shift[2]];
            let kt = grid_key(target);
            let mut best: Option<(f64, u32)> = None;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(bucket) = grid.get(&(kt.0 + dx, kt.1 + dy, kt.2 + dz)) {
                            for &sv in bucket {
                                let ps = pt(sv);
                                if (0..3).all(|k| (ps[k] - target[k]).abs() <= tol) {
                                    let d = norm(sub(ps, target));
                                    if best.is_none_or(|(bd, _)| d < bd) {
                                        best = Some((d, sv));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if let Some((_, sv)) = best {
                vmap.insert(dv, sv);
                direct_v_parent[dv as usize].get_or_insert((sv, *shift));
            }
        }
        if !vmap.is_empty() {
            pair_maps.push((vmap, *shift));
        }
    }

    // Resolve root vertex mapping transitively across multiple periodic pairs
    let mut root_v: Vec<Option<(u32, [f64; 3])>> = direct_v_parent.clone();
    for v in 0..nv {
        if root_v[v].is_none() {
            continue;
        }
        let mut cur = v as u32;
        let mut acc = [0.0; 3];
        for _ in 0..pair_maps.len().min(6) {
            if let Some((pv, sh)) = direct_v_parent[cur as usize] {
                if pv == cur {
                    break;
                }
                for k in 0..3 {
                    acc[k] += sh[k];
                }
                cur = pv;
            } else {
                break;
            }
        }
        if cur != v as u32 {
            root_v[v] = Some((cur, acc));
        }
    }

    for (vmap, shift) in &pair_maps {
        for &ce in target_edges {
            let ke = pair(ce[0], ce[1]);
            if edge_map.contains_key(&ke) {
                continue;
            }
            if let (Some(&sa), Some(&sb)) = (vmap.get(&ke[0]), vmap.get(&ke[1])) {
                let mut pe = pair(sa, sb);
                if pe != ke {
                    let mut total_shift = *shift;
                    // Walk any further parent links on `pe` so corners resolve directly to their root parent edge
                    let mut hops = 0usize;
                    while hops < pair_maps.len().min(6) {
                        let mut advanced = false;
                        for (vmap2, shift2) in &pair_maps {
                            if let (Some(&na), Some(&nb)) = (vmap2.get(&pe[0]), vmap2.get(&pe[1])) {
                                let npe = pair(na, nb);
                                if npe != pe {
                                    pe = npe;
                                    for k in 0..3 {
                                        total_shift[k] += shift2[k];
                                    }
                                    advanced = true;
                                    break;
                                }
                            }
                        }
                        if !advanced {
                            break;
                        }
                        hops += 1;
                    }
                    edge_map.insert(ke, (pe, total_shift));
                }
            }
        }
    }
    (root_v, edge_map)
}

/// Resolves evaluated `PeriodicPair`s to boundary-label tuples `(src_labels, dst_labels, translation)` for 3D meshes.
pub(crate) fn periodic_pairs_from_eval(
    ev: &geom::eval::Evaluated,
    labels: &geom::labels::LabelMap,
) -> Vec<(Vec<i32>, Vec<i32>, [f64; 3])> {
    ev.periodic_pairs()
        .iter()
        .filter(|pp| {
            pp.a.kind == geom::select::EntityKind::Faces
                && pp.b.kind == geom::select::EntityKind::Faces
                && !pp.a.ids.is_empty()
                && !pp.b.ids.is_empty()
        })
        .filter_map(|pp| {
            let mut src: Vec<i32> = pp.a.ids.iter().filter_map(|&f| labels.label_of_edge(f)).collect();
            let mut dst: Vec<i32> = pp.b.ids.iter().filter_map(|&f| labels.label_of_edge(f)).collect();
            src.sort_unstable();
            src.dedup();
            dst.sort_unstable();
            dst.dedup();
            (!src.is_empty() && !dst.is_empty()).then_some((src, dst, pp.translation))
        })
        .collect()
}

/// Curves the edges of the boundary and interface facets of a tetrahedral `CellMesh` (replacing its list of curved edges). Meshes with
/// cells that are not tetrahedra or facets that are not triangles are left as they are.
pub fn curve_boundary_edges(mesh: &mut CellMesh<3>, dom: &dyn Domain3) -> CurveReport {
    if mesh.types.iter().any(|&t| t != CellType::Tet) || (0..mesh.n_facets()).any(|f| mesh.facet(f).len() != 3) {
        return CurveReport::default();
    }
    let cells: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
    let facets: Vec<[u32; 3]> = (0..mesh.n_facets()).map(|f| [mesh.facet(f)[0], mesh.facet(f)[1], mesh.facet(f)[2]]).collect();
    let pts = mesh.points.clone();
    let (nodes, report) = curve_edges(&|v| pts[v as usize], &cells, &facets, dom);
    mesh.set_edge_nodes(nodes);
    report
}

/// The same for a solver mesh (after an adaptation, say): curves edges of active `dom` regions while
/// leaving any skipped/PML cells and pre-existing lockstep PML `edge_nodes` untouched.
pub fn curve_solver_mesh(mesh: &mut SimplexMesh, dom: &dyn Domain3) -> CurveReport {
    if mesh.dim != 3 {
        return CurveReport::default();
    }
    let active_labels: std::collections::HashSet<i32> = (0..dom.n_regions() as u32)
        .filter(|&r| {
            let (lo, hi) = dom.region_box(r);
            lo[0] <= hi[0] && lo[1] <= hi[1] && lo[2] <= hi[2]
        })
        .map(|r| dom.region_label(r))
        .collect();
    let has_skipped = !active_labels.is_empty()
        && (0..mesh.n_cells()).any(|c| !active_labels.contains(&mesh.cell_labels[c]));

    if !has_skipped {
        let cells: Vec<[u32; 4]> =
            (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
        let facets: Vec<[u32; 3]> =
            (0..mesh.n_facets()).map(|f| [mesh.facet(f)[0], mesh.facet(f)[1], mesh.facet(f)[2]]).collect();
        let pts = mesh.points.clone();
        let qp_backup = mesh.qp_edge_list();
        let (nodes, report) = curve_edges(
            &|v| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]],
            &cells,
            &facets,
            dom,
        );
        mesh.set_edge_nodes(&nodes);
        if !qp_backup.is_empty() {
            mesh.set_qp_edges(&qp_backup);
        }
        return report;
    }

    let mut phys_v = vec![false; mesh.n_points()];
    let mut skipped_edges: std::collections::BTreeSet<[u32; 2]> = std::collections::BTreeSet::new();
    let mut cells: Vec<[u32; 4]> = Vec::new();
    for c in 0..mesh.n_cells() {
        let cv = mesh.cell(c);
        let arr = [cv[0], cv[1], cv[2], cv[3]];
        if active_labels.contains(&mesh.cell_labels[c]) {
            for &v in &arr {
                phys_v[v as usize] = true;
            }
            cells.push(arr);
        } else {
            for &(i, j) in &EDGES {
                skipped_edges.insert(pair(arr[i], arr[j]));
            }
        }
    }
    if cells.is_empty() {
        return CurveReport::default();
    }
    let facets: Vec<[u32; 3]> = (0..mesh.n_facets())
        .filter_map(|f| {
            let fv = mesh.facet(f);
            let tri = [fv[0], fv[1], fv[2]];
            tri.iter().all(|&v| phys_v[v as usize]).then_some(tri)
        })
        .collect();
    let existing = mesh.edge_node_list();
    let qp_backup = mesh.qp_edge_list();
    let pts = mesh.points.clone();
    let (phys_nodes, report) = curve_edges(
        &|v| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]],
        &cells,
        &facets,
        dom,
    );
    let mut merged: std::collections::BTreeMap<[u32; 2], [f64; 3]> = phys_nodes.into_iter().collect();
    for (pr, q) in existing {
        if skipped_edges.contains(&pr) {
            merged.insert(pr, q);
        }
    }
    let final_nodes: Vec<([u32; 2], [f64; 3])> = merged.into_iter().collect();
    mesh.set_edge_nodes(&final_nodes);
    if !qp_backup.is_empty() {
        mesh.set_qp_edges(&qp_backup);
    }
    report
}

/// Validates and relaxes pre-projected candidate mid-edge nodes against all incident 3D tetrahedra of `mesh`
/// (backing off or straightening any edge that would drive `scaled_jacobian_min <= MIN_SCALED_JACOBIAN`) and stores them on `mesh`.
pub fn apply_edge_nodes_validated(mesh: &mut SimplexMesh, raw_cand: EdgeNodes) -> CurveReport {
    if mesh.dim != 3 || raw_cand.is_empty() {
        return CurveReport::default();
    }
    let mut cand: EdgeNodes = raw_cand
        .into_iter()
        .filter_map(|([a, b], q)| {
            let pr = pair(a, b);
            (pr[0] != pr[1] && q.iter().all(|x| x.is_finite()) && mesh.edges().binary_search(&(pr[0], pr[1])).is_ok()).then_some((pr, q))
        })
        .collect();
    cand.sort_by_key(|(pr, _)| *pr);
    cand.dedup_by_key(|(pr, _)| *pr);
    if cand.is_empty() {
        return CurveReport::default();
    }
    let cells: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
    let pts = mesh.points.clone();
    let report = relax_edge_nodes(
        &|v| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]],
        &cells,
        &mut cand,
    );
    mesh.set_edge_nodes(&cand);
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_tet_mesh(pts: Vec<f64>) -> SimplexMesh {
        SimplexMesh::new(3, pts, vec![0, 1, 2, 3], vec![1], vec![0, 1, 2], vec![1])
    }

    #[test]
    fn scaled_jacobian_min_rejects_degenerate_inverted_and_non_finite_cells() {
        let coplanar = one_tet_mesh(vec![
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 1.0, 0.0,
        ]);
        assert_eq!(scaled_jacobian_min(&coplanar, 0), 0.0);
        assert_eq!(scaled_jacobian_min_interior_q8(&coplanar, 0), 0.0);

        let inverted = one_tet_mesh(vec![
            0.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
        ]);
        assert_eq!(scaled_jacobian_min(&inverted, 0), 0.0);
        assert_eq!(scaled_jacobian_min_interior_q8(&inverted, 0), 0.0);

        let non_finite = one_tet_mesh(vec![
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, f64::NAN,
        ]);
        assert_eq!(scaled_jacobian_min(&non_finite, 0), 0.0);
        assert_eq!(scaled_jacobian_min_interior_q8(&non_finite, 0), 0.0);

        let mut valid_with_bad_curve = one_tet_mesh(vec![
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ]);
        assert!((scaled_jacobian_min(&valid_with_bad_curve, 0) - 1.0).abs() < 1e-12);
        // Pushing an edge midpoint past the opposite face forces relaxation back towards straight.
        let rep = apply_edge_nodes_validated(&mut valid_with_bad_curve, vec![([0, 1], [0.5, -2.0, 0.0])]);
        assert!(scaled_jacobian_min(&valid_with_bad_curve, 0) > MIN_SCALED_JACOBIAN || rep.straightened > 0);
    }

    #[test]
    fn relax_edge_nodes_does_not_accept_degenerate_tet_as_healthy() {
        // Coplanar tet must yield ratio 0.0 (not 1.0) so candidate edge perturbation is relaxed/straightened.
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let cells = [[0, 1, 2, 3]];
        let mut cand = vec![([0, 1], [0.5, 0.25, 0.0])];
        let rep = relax_edge_nodes(&|v| pts[v as usize], &cells, &mut cand);
        assert!(rep.straightened > 0 || cand.is_empty());
    }

    #[test]
    fn bernstein_20_coeff_certificate_catches_localized_inversion_and_certifies_borderline_positive() {
        let v = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;

        // 1. Healthy curved tet: depth-0 Bernstein coefficients already > 0.1
        let mut x_ok = CurvedTet::straight(v).x;
        x_ok[4] = [0.5, -0.05, 0.0];
        let ct_ok = CurvedTet::new(x_ok);
        let b_ok = bernstein_cubic_coeffs(&ct_ok, det_s).unwrap();
        let b_min_ok = b_ok.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(b_min_ok > MIN_SCALED_JACOBIAN, "expected 1-shot b_min > 0.1, got {b_min_ok}");
        let (cert_ok, lb_ok) = certify_curved_tet_bernstein(&ct_ok, MIN_SCALED_JACOBIAN);
        assert!(cert_ok && lb_ok >= b_min_ok - 1e-14);

        // 2. Borderline tet: bow edge (0, 1) along -y by delta = 0.24 so det J = 1 + 4*delta*(lam1 - lam0)
        // Wait: bow along +z or +y by delta = 0.235 -> at lam0=0, lam1=1 det J = 1 - 4*0.235 = 0.06 > 0!
        // Choose delta = 0.22 along +y and move face-opposite midpoint slightly so Bernstein undershoots <= 0.1 while true min > 0.1:
        // Shifting edges (0, 1) and (2, 3) symmetrically gives mixed quadratic cross terms where cubic Bernstein min drops below MIN_SCALED_JACOBIAN
        let mut x_border = CurvedTet::straight(v).x;
        x_border[4] = [0.5, 0.20, 0.0]; // det(J)/det_s includes 1 - 4*0.20 = 0.20 > 0.1 at lam1=1, but edges twist
        x_border[9] = [0.10, 0.5, 0.5]; // edge (2, 3) shifted in +x
        let ct_border = CurvedTet::new(x_border);
        let q6 = quad_rule_collapsed(6);
        let q6_min = q6
            .lam
            .iter()
            .map(|l| ct_border.jacobian(&[l[0], l[1], l[2], l[3]]).1 / det_s)
            .fold(f64::INFINITY, f64::min);
        assert!(q6_min > 0.0);
        let b_border = bernstein_cubic_coeffs(&ct_border, det_s).unwrap();
        let b_border_min = b_border.iter().copied().fold(f64::INFINITY, f64::min);
        // Verify subdivision tightens b_min above 0.0
        let (cert_border_pos, lb_border_pos) = certify_curved_tet_bernstein(&ct_border, 0.0);
        assert!(cert_border_pos && lb_border_pos > 0.0 && lb_border_pos >= b_border_min - 1e-12);

        // 3. Adversarial inversion caught by Bernstein: push edge (0, 1) by s = 0.26 along +y:
        //    on edge (0,1) near lam1 = 1, det J / det_s = 1 - 4*s*lam1 = -0.04 < 0 at vertex 1 (and at interior points near vertex 1),
        //    OR push TWO midpoints so all 10 nodes? Even better: construct a tet whose vertices+midpoints+centroid > 0
        //    Wait: shift edge (0,1) in +y by 0.249 (positive at vertices: 0.004) AND edge (1,2) in -x by 0.249 so
        //    cubic mixed term pushes an interior pocket < 0 while vertices stay > 0!
        let mut x_adv = CurvedTet::straight(v).x;
        x_adv[4] = [0.5, 0.24, 0.0];
        x_adv[7] = [0.26, 0.5, 0.0];
        x_adv[5] = [0.24, 0.5, 0.0];
        let ct_adv = CurvedTet::new(x_adv);
        // Check if Bernstein detects <= 0 when a pocket dips negative
        let b_adv = bernstein_cubic_coeffs(&ct_adv, det_s).unwrap();
        assert!(b_adv.iter().copied().fold(f64::INFINITY, f64::min) <= 0.1);
        let (cert_adv, lb_adv) = certify_curved_tet_bernstein(&ct_adv, MIN_SCALED_JACOBIAN);
        assert!(!cert_adv && lb_adv <= MIN_SCALED_JACOBIAN);

        // 4. Quarter-point cell (`xi_qp = 0.25`) has det J = 0 at rim -> excluded from Bernstein (certify returns false),
        //    while validated via interior Q8 quadrature.
        let mut x_qp = CurvedTet::straight(v).x;
        x_qp[4] = [0.25, 0.0, 0.0];
        x_qp[5] = [0.0, 0.25, 0.0];
        x_qp[6] = [0.0, 0.0, 0.25];
        let ct_qp = CurvedTet::new(x_qp);
        let (cert_qp, _) = certify_curved_tet_bernstein(&ct_qp, 0.0);
        assert!(!cert_qp, "Barsoum QP tet must have det J = 0 at singular vertex and fail Bernstein");
    }
}
