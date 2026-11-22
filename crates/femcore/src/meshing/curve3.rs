//! Curved edges for 3D meshes (ADR 82): the mid-edge geometry nodes that make the tetrahedra around a curved surface quadratic
//! (isoparametric) instead of straight.
//!
//! For every edge of a boundary facet or of a facet between two regions the geometry (`Domain3::project_edge`) says where the midpoint of
//! the edge would lie on the exact surface; if it moves by more than `1e-9` of the edge length (a planar face adds nothing) and by less
//! than a quarter of it, the node is kept and the cells that own the edge become curved. Then every curved cell is checked: the
//! determinant of its Jacobian must stay above a tenth of the straight cell's at the ten nodes and at the points of a degree-6 rule;
//! if it does not, the offsets of that cell's curved edges are halved (every cell around an edge sees the same offset, so the mesh stays
//! conforming), and after four halvings the edges are made straight. Everything is decided per edge and per cell in a fixed order, so the
//! result does not depend on the number of threads (ADR 11).
//!
//! An edge whose two ends lie on two common faces (a crease between two surfaces) is left straight (`Domain3::edge_on_crease`); an edge
//! that is not on a curved surface is never touched. Anything that changes the connectivity afterwards (refinement, adaptation) must run
//! before this, or drop the list.

use crate::fem::{quad_rule_collapsed, CurvedTet};
use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::SimplexMesh;
use crate::meshing::mesher3d::Domain3;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};

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
    let mut nodes: Vec<Option<[f64; 3]>> = edges
        .par_iter()
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
    if cand.is_empty() {
        return (vec![], CurveReport::default());
    }
    validate_and_relax_candidates(pos, cells, cand, &periodic_edge_parent)
}

fn validate_and_relax_candidates(
    pos: &(dyn Fn(u32) -> [f64; 3] + Sync),
    cells: &[[u32; 4]],
    cand: EdgeNodes,
    periodic_edge_parent: &HashMap<[u32; 2], ([u32; 2], [f64; 3])>,
) -> (EdgeNodes, CurveReport) {
    if cand.is_empty() {
        return (vec![], CurveReport::default());
    }
    let index: HashMap<[u32; 2], usize> = cand.iter().enumerate().map(|(i, (e, _))| (*e, i)).collect();
    let mut level = vec![0u8; cand.len()];

    // 2. the cells with a candidate edge, and for each the candidate index of each of its six edges
    let curved_cells: Vec<(usize, [Option<usize>; 6])> = cells
        .iter()
        .enumerate()
        .filter_map(|(ci, c)| {
            let mut ids = [None; 6];
            let mut any = false;
            for (k, &(i, j)) in EDGES.iter().enumerate() {
                if let Some(&e) = index.get(&pair(c[i], c[j])) {
                    ids[k] = Some(e);
                    any = true;
                }
            }
            any.then_some((ci, ids))
        })
        .collect();

    // 3. validity, with the offsets of the edges halved where a cell fails
    let rule = quad_rule_collapsed(6);
    let mut test: Vec<[f64; 4]> = vec![[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]];
    for &(i, j) in EDGES.iter() {
        let mut l = [0.0; 4];
        l[i] = 0.5;
        l[j] = 0.5;
        test.push(l);
    }
    test.extend(rule.lam.iter().map(|l| [l[0], l[1], l[2], l[3]]));
    let scale_of = |lv: u8| -> f64 {
        if lv > MAX_HALVINGS {
            0.0
        } else {
            0.5f64.powi(lv as i32)
        }
    };
    let cell_ratio = |ci: usize, ids: &[Option<usize>; 6], level: &[u8]| -> f64 {
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
        // the straight determinant (constant): the Jacobian of the straight cell
        let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
        if det_s == 0.0 || !det_s.is_finite() {
            return 1.0;
        }
        test.iter().map(|l| ct.jacobian(l).1 / det_s).fold(f64::INFINITY, f64::min)
    };
    for _round in 0..=(MAX_HALVINGS as usize + 1) {
        let failing: Vec<usize> = curved_cells
            .par_iter()
            .enumerate()
            .filter_map(|(k, (ci, ids))| (cell_ratio(*ci, ids, &level) <= MIN_SCALED_JACOBIAN).then_some(k))
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
            for (ci, (ce, _)) in cand.iter().enumerate() {
                if let Some((pe, _)) = periodic_edge_parent.get(ce) {
                    if let Some(&pi) = index.get(pe) {
                        let mx = level[ci].max(level[pi]);
                        level[ci] = mx;
                        level[pi] = mx;
                    }
                }
            }
        }
    }

    // 4. the result
    let mut report = CurveReport::default();
    let mut out: Vec<([u32; 2], [f64; 3])> = Vec::new();
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
            out.push((*pr, [mid[0] + s * (q[0] - mid[0]), mid[1] + s * (q[1] - mid[1]), mid[2] + s * (q[2] - mid[2])]));
        }
    }
    report.curved_edges = out.len();
    let kept: HashMap<[u32; 2], ()> = out.iter().map(|(e, _)| (*e, ())).collect();
    let mut min_sj = f64::INFINITY;
    let mut n_cells = 0;
    for (ci, ids) in &curved_cells {
        if ids.iter().flatten().any(|&e| kept.contains_key(&cand[e].0)) {
            n_cells += 1;
            min_sj = min_sj.min(cell_ratio(*ci, ids, &level));
        }
    }
    report.curved_cells = n_cells;
    report.min_scaled_jacobian = if min_sj.is_finite() { min_sj } else { 1.0 };
    (out, report)
}

/// The smallest ratio of the Jacobian determinant to the straight cell's over cell `c` of a solver mesh (ten nodes and a degree-6 rule).
pub fn scaled_jacobian_min(mesh: &SimplexMesh, c: usize) -> f64 {
    let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
    let v = [ct.x[0], ct.x[1], ct.x[2], ct.x[3]];
    let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
    if det_s == 0.0 || !det_s.is_finite() {
        return 1.0;
    }
    let rule = quad_rule_collapsed(6);
    let mut m = f64::INFINITY;
    for l in rule.lam.iter().map(|l| [l[0], l[1], l[2], l[3]]).chain(node_points()) {
        m = m.min(ct.jacobian(&l).1 / det_s);
    }
    m
}

/// The four vertices and the six edge midpoints in barycentric coordinates.
fn node_points() -> Vec<[f64; 4]> {
    let mut v: Vec<[f64; 4]> = vec![[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]];
    for &(i, j) in EDGES.iter() {
        let mut l = [0.0; 4];
        l[i] = 0.5;
        l[j] = 0.5;
        v.push(l);
    }
    v
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

/// The same for a solver mesh (after an adaptation, say): replaces its curved edges.
pub fn curve_solver_mesh(mesh: &mut SimplexMesh, dom: &dyn Domain3) -> CurveReport {
    if mesh.dim != 3 {
        return CurveReport::default();
    }
    let cells: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
    let facets: Vec<[u32; 3]> = (0..mesh.n_facets()).map(|f| [mesh.facet(f)[0], mesh.facet(f)[1], mesh.facet(f)[2]]).collect();
    let pts = mesh.points.clone();
    let (nodes, report) = curve_edges(&|v| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]], &cells, &facets, dom);
    mesh.set_edge_nodes(&nodes);
    report
}

/// Validates and relaxes pre-projected candidate mid-edge nodes against all incident 3D tetrahedra of `mesh`
/// (backing off or straightening any edge that would drive `scaled_jacobian_min <= MIN_SCALED_JACOBIAN`) and stores them on `mesh`.
pub(crate) fn apply_edge_nodes_validated(mesh: &mut SimplexMesh, raw_cand: EdgeNodes) -> CurveReport {
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
    let empty_per = HashMap::new();
    let (nodes, report) = validate_and_relax_candidates(
        &|v| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]],
        &cells,
        cand,
        &empty_per,
    );
    mesh.set_edge_nodes(&nodes);
    report
}
