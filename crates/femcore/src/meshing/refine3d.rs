//! Local refinement of a tetrahedral mesh by longest-edge (Rivara) bisection.
//!
//! A cell whose longest edge is larger than the wanted element size at its centroid is marked. For a
//! marked cell the algorithm follows the longest edge to a *terminal* edge (one that is the longest
//! edge of every cell around it) and splits all cells around that edge in two. Because every cell
//! that contains the edge is split along it, the mesh stays conforming, and because terminal edges
//! are the longest ones the cells keep reasonable shapes. Boundary facets are split together with the
//! cells, and the new vertex on an edge that lies on the boundary of the solid or on an interface
//! between domains is moved onto the true (curved) surface by the caller's projection.

use crate::mesh::cells::{CellMesh, CellType};
use crate::meshing::mesher3d::PeriodicAxis3;
use std::collections::HashMap;

type P3 = [f64; 3];

/// Where an edge that is being split lies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeKind {
    Interior,
    /// on a boundary facet of the mesh
    Boundary,
    /// inside the mesh, between cells of different domains
    Interface,
}

#[derive(Clone, Debug, Default)]
pub struct RefineReport {
    pub rounds: usize,
    pub bisections: usize,
    pub cells_before: usize,
    pub cells_after: usize,
}

fn key(a: u32, b: u32) -> u64 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    ((lo as u64) << 32) | hi as u64
}

fn dist2(p: P3, q: P3) -> f64 {
    (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)
}

const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

fn shared_periodic_axes(bits_a: u8, bits_b: u8) -> u8 {
    let common = bits_a & bits_b;
    (common | (common >> 3)) & 0b111
}

fn twin_with_child_mask(v: u32, target_child_mask: u8, parent_v: &[u32], plane_bits: &[u8], twins_of: &HashMap<u32, Vec<(u32, u8, P3)>>) -> u32 {
    if (plane_bits[v as usize] & 0b111) == target_child_mask {
        return v;
    }
    let p = parent_v[v as usize];
    if let Some(list) = twins_of.get(&p) {
        for &(tv, mask, _) in list {
            if mask == target_child_mask {
                return tv;
            }
        }
    }
    v
}

/// Canonical key of an edge `(a, b)`: when both endpoints lie on the same periodic plane(s), reduces both endpoints along
/// their shared periodic axes so periodic twin edges share the exact same tie-break key.
fn edge_order_key(a: u32, b: u32, parent_v: &[u32], plane_bits: &[u8], twins_of: &HashMap<u32, Vec<(u32, u8, P3)>>) -> u64 {
    let axes = shared_periodic_axes(plane_bits[a as usize], plane_bits[b as usize]);
    if axes != 0 {
        let pa = twin_with_child_mask(a, (plane_bits[a as usize] & 0b111) & !axes, parent_v, plane_bits, twins_of);
        let pb = twin_with_child_mask(b, (plane_bits[b as usize] & 0b111) & !axes, parent_v, plane_bits, twins_of);
        if pa != pb {
            return key(pa, pb);
        }
    }
    key(a, b)
}

/// The longest edge of a cell under a strict total order (length, then canonical vertex numbers), so two cells
/// that share an edge (or periodic twin edges) always agree on which of two edges is longer.
fn longest(cell: &[u32; 4], pts: &[P3], parent_v: &[u32], plane_bits: &[u8], twins_of: &HashMap<u32, Vec<(u32, u8, P3)>>) -> (u32, u32) {
    let mut best = (0u64, 0u64);
    let mut edge = (cell[0], cell[1]);
    for &(i, j) in &PAIRS {
        let (a, b) = (cell[i], cell[j]);
        let l = (dist2(pts[a as usize], pts[b as usize]).to_bits() >> 8, edge_order_key(a, b, parent_v, plane_bits, twins_of));
        if l > best {
            best = l;
            edge = (a, b);
        }
    }
    edge
}

fn volume6(p: [P3; 4]) -> f64 {
    let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
    let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
    let e3 = [p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2]];
    e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])
}

/// Refines `mesh` until no cell has an edge longer than `ratio` times `size(centroid)`.
///
/// `project(midpoint, a, b, kind)` returns the position of the new vertex on edge `a`-`b`; it is
/// only asked to move the midpoint for edges on the boundary or on an interface.
pub fn refine_to_size(
    mesh: &CellMesh<3>,
    size: &(dyn Fn(P3) -> f64 + Sync),
    project: &(dyn Fn(P3, P3, P3, EdgeKind) -> P3 + Sync),
    ratio: f64,
    max_cells: usize,
) -> Result<(CellMesh<3>, RefineReport), String> {
    refine_to_size_periodic(mesh, size, project, ratio, max_cells, &[])
}

/// Like [`refine_to_size`], keeping opposite periodic faces listed in `periodic_axes` congruent vertex for vertex,
/// edge for edge, and triangle for triangle by bisecting periodic twin boundary edges synchronously from parent to child.
pub fn refine_to_size_periodic(
    mesh: &CellMesh<3>,
    size: &(dyn Fn(P3) -> f64 + Sync),
    project: &(dyn Fn(P3, P3, P3, EdgeKind) -> P3 + Sync),
    ratio: f64,
    max_cells: usize,
    periodic_axes: &[PeriodicAxis3],
) -> Result<(CellMesh<3>, RefineReport), String> {
    if mesh.types.iter().any(|&t| t != CellType::Tet) || (0..mesh.n_facets()).any(|f| mesh.facet(f).len() != 3) {
        return Err("local refinement works on tetrahedral meshes with triangular boundary facets".into());
    }
    let mut pts: Vec<P3> = mesh.points.clone();
    let mut cells: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
    let mut labels = mesh.label.clone();
    let mut facets: Vec<[u32; 3]> = (0..mesh.n_facets()).map(|f| [mesh.facet(f)[0], mesh.facet(f)[1], mesh.facet(f)[2]]).collect();
    let mut flabels = mesh.facet_label.clone();
    let mut report = RefineReport { cells_before: cells.len(), ..Default::default() };

    let mut per_axis: [Option<(f64, f64)>; 3] = [None; 3];
    for pa in periodic_axes {
        if pa.axis < 3 && pa.hi > pa.lo {
            per_axis[pa.axis] = Some((pa.lo, pa.hi));
        }
    }
    let has_periodic = per_axis.iter().any(|x| x.is_some());
    let (mut blo, mut bhi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &pts {
        for k in 0..3 {
            blo[k] = blo[k].min(p[k]);
            bhi[k] = bhi[k].max(p[k]);
        }
    }
    let extent = (0..3).map(|k| bhi[k] - blo[k]).fold(0.0f64, f64::max).max(1e-300);
    let tol = 1e-10 * extent;

    let classify_planes = |p: P3| -> (u8, P3) {
        let mut bits = 0u8;
        let mut shift = [0.0; 3];
        for (k, slot) in per_axis.iter().enumerate() {
            if let Some((lo, hi)) = *slot {
                if (p[k] - lo).abs() <= tol {
                    bits |= 1 << (k + 3);
                } else if (p[k] - hi).abs() <= tol {
                    bits |= 1 << k;
                    shift[k] = hi - lo;
                }
            }
        }
        (bits, shift)
    };

    let mut plane_bits: Vec<u8> = vec![0; pts.len()];
    let mut parent_v: Vec<u32> = (0..pts.len() as u32).collect();
    let mut shift_v: Vec<P3> = vec![[0.0; 3]; pts.len()];
    let mut twins_of: HashMap<u32, Vec<(u32, u8, P3)>> = HashMap::new();

    if has_periodic {
        let cell_sz = 4.0 * tol;
        let gkey = |p: P3| -> [i64; 3] {
            [(p[0] / cell_sz).floor() as i64, (p[1] / cell_sz).floor() as i64, (p[2] / cell_sz).floor() as i64]
        };
        let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        for (v, &p) in pts.iter().enumerate() {
            let (bits, shift) = classify_planes(p);
            plane_bits[v] = bits;
            shift_v[v] = shift;
            if bits != 0 {
                grid.entry(gkey(p)).or_default().push(v as u32);
            }
        }
        for v in 0..pts.len() {
            let child_mask = plane_bits[v] & 0b111;
            if child_mask == 0 {
                continue;
            }
            let shift = shift_v[v];
            let target = [pts[v][0] - shift[0], pts[v][1] - shift[1], pts[v][2] - shift[2]];
            let kt = gkey(target);
            let mut found = None;
            let mut best_d2 = f64::INFINITY;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&[kt[0] + dx, kt[1] + dy, kt[2] + dz]) {
                            for &cand in list {
                                if (plane_bits[cand as usize] & 0b111) != 0 {
                                    continue;
                                }
                                let pc = pts[cand as usize];
                                if (0..3).all(|k| (pc[k] - target[k]).abs() <= tol) {
                                    let d2 = dist2(pc, target);
                                    if d2 < best_d2 {
                                        best_d2 = d2;
                                        found = Some(cand);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if let Some(p) = found {
                parent_v[v] = p;
            }
        }
        for v in 0..pts.len() {
            if plane_bits[v] != 0 {
                let p = parent_v[v];
                twins_of.entry(p).or_default().push((v as u32, plane_bits[v] & 0b111, shift_v[v]));
            }
        }
    }

    let twin_group = |e: (u32, u32),
                      parent_v: &[u32],
                      plane_bits: &[u8],
                      twins_of: &HashMap<u32, Vec<(u32, u8, P3)>>,
                      around: &HashMap<u64, Vec<u32>>|
     -> Vec<(u32, u32, u8, P3)> {
        let axes = shared_periodic_axes(plane_bits[e.0 as usize], plane_bits[e.1 as usize]);
        if !has_periodic || axes == 0 {
            return vec![(e.0, e.1, 0, [0.0; 3])];
        }
        let base_ma = (plane_bits[e.0 as usize] & 0b111) & !axes;
        let base_mb = (plane_bits[e.1 as usize] & 0b111) & !axes;
        let pa0 = twin_with_child_mask(e.0, base_ma, parent_v, plane_bits, twins_of);
        let pb0 = twin_with_child_mask(e.1, base_mb, parent_v, plane_bits, twins_of);
        let (pa0, pb0, base_ma, base_mb) = if pa0 <= pb0 { (pa0, pb0, base_ma, base_mb) } else { (pb0, pa0, base_mb, base_ma) };
        let mut out: Vec<(u32, u32, u8, P3)> = Vec::new();
        for sub in 0..8u8 {
            if (sub & !axes) != 0 {
                continue;
            }
            let va = twin_with_child_mask(pa0, base_ma | sub, parent_v, plane_bits, twins_of);
            let vb = twin_with_child_mask(pb0, base_mb | sub, parent_v, plane_bits, twins_of);
            if va != vb && around.contains_key(&key(va, vb)) {
                let mut shift = [0.0; 3];
                for (k, slot) in per_axis.iter().enumerate() {
                    if let Some((lo, hi)) = *slot {
                        if (sub & (1 << k)) != 0 {
                            shift[k] = hi - lo;
                        }
                    }
                }
                out.push((va, vb, sub, shift));
            }
        }
        if out.is_empty() {
            vec![(e.0, e.1, 0, [0.0; 3])]
        } else {
            out
        }
    };

    for _round in 0..200 {
        // 1. cells that are too large
        let marked: Vec<usize> = (0..cells.len())
            .filter(|&c| {
                let (a, b) = longest(&cells[c], &pts, &parent_v, &plane_bits, &twins_of);
                let cen = {
                    let t = cells[c].map(|v| pts[v as usize]);
                    [(t[0][0] + t[1][0] + t[2][0] + t[3][0]) / 4.0, (t[0][1] + t[1][1] + t[2][1] + t[3][1]) / 4.0, (t[0][2] + t[1][2] + t[2][2] + t[3][2]) / 4.0]
                };
                let s = ratio * size(cen);
                dist2(pts[a as usize], pts[b as usize]) > s * s
            })
            .collect();
        if marked.is_empty() {
            break;
        }
        if cells.len() > max_cells {
            return Err(format!("local refinement would exceed {max_cells} cells: use a larger element size or fewer refined regions"));
        }
        report.rounds += 1;

        // 2. cells around every edge, and boundary facets around every boundary edge
        let mut around: HashMap<u64, Vec<u32>> = HashMap::with_capacity(cells.len() * 2);
        for (c, t) in cells.iter().enumerate() {
            for &(i, j) in &PAIRS {
                around.entry(key(t[i], t[j])).or_default().push(c as u32);
            }
        }
        let mut facet_edges: HashMap<u64, Vec<u32>> = HashMap::with_capacity(facets.len() * 2);
        for (f, t) in facets.iter().enumerate() {
            for (i, j) in [(0, 1), (0, 2), (1, 2)] {
                facet_edges.entry(key(t[i], t[j])).or_default().push(f as u32);
            }
        }

        // 3. terminal edges (or synchronized periodic twin edge groups) whose stars do not overlap
        let mut claimed = vec![false; cells.len()];
        let mut selected: Vec<Vec<(u32, u32, u8, P3)>> = vec![];
        for &c0 in &marked {
            if claimed[c0] {
                continue;
            }
            let mut c = c0;
            let mut steps = 0;
            let group = loop {
                let e = longest(&cells[c], &pts, &parent_v, &plane_bits, &twins_of);
                let grp = twin_group(e, &parent_v, &plane_bits, &twins_of, &around);
                let mut next_cell: Option<u32> = None;
                'search: for &(ga, gb, _, _) in &grp {
                    let gk = key(ga, gb);
                    if let Some(star) = around.get(&gk) {
                        for &s in star {
                            let le = longest(&cells[s as usize], &pts, &parent_v, &plane_bits, &twins_of);
                            if key(le.0, le.1) != gk {
                                next_cell = Some(s);
                                break 'search;
                            }
                        }
                    }
                }
                match next_cell {
                    None => break grp,
                    Some(s) => c = s as usize,
                }
                steps += 1;
                if steps > 100_000 {
                    return Err("local refinement did not find a terminal edge".into());
                }
            };
            let overlaps = group.iter().any(|&(ga, gb, _, _)| {
                around.get(&key(ga, gb)).is_some_and(|star| star.iter().any(|&s| claimed[s as usize]))
            });
            if overlaps {
                continue;
            }
            for &(ga, gb, _, _) in &group {
                if let Some(star) = around.get(&key(ga, gb)) {
                    for &s in star {
                        claimed[s as usize] = true;
                    }
                }
            }
            selected.push(group);
        }

        // 4. split them (synchronously across periodic twins from parent to child)
        for group in selected {
            let (pa, pb, _, _) = group[0];
            let pk = key(pa, pb);
            let pstar = &around[&pk];
            let pids: Vec<u32> = facet_edges.get(&pk).cloned().unwrap_or_default();
            let kind = if !pids.is_empty() {
                EdgeKind::Boundary
            } else if pstar.iter().any(|&s| labels[s as usize] != labels[pstar[0] as usize]) {
                EdgeKind::Interface
            } else {
                EdgeKind::Interior
            };
            let (ppa, ppb) = (pts[pa as usize], pts[pb as usize]);
            let mid_p = [(ppa[0] + ppb[0]) / 2.0, (ppa[1] + ppb[1]) / 2.0, (ppa[2] + ppb[2]) / 2.0];
            let mut mp = if kind == EdgeKind::Interior { mid_p } else { project(mid_p, ppa, ppb, kind) };
            let shared_p_planes = plane_bits[pa as usize] & plane_bits[pb as usize];
            for (k, slot) in per_axis.iter().enumerate() {
                if let Some((lo, hi)) = *slot {
                    if (shared_p_planes & (1 << (k + 3))) != 0 {
                        mp[k] = lo;
                    } else if (shared_p_planes & (1 << k)) != 0 {
                        mp[k] = hi;
                    }
                }
            }
            // the projection must keep every new cell positively oriented and not squash it across all twins
            let ok_all = |mp_cand: P3| -> bool {
                group.iter().all(|&(ga, gb, mask, shift)| {
                    let mut m = [mp_cand[0] + shift[0], mp_cand[1] + shift[1], mp_cand[2] + shift[2]];
                    for (k, slot) in per_axis.iter().enumerate() {
                        if let Some((_, hi)) = *slot {
                            if (mask & (1 << k)) != 0 {
                                m[k] = hi;
                            }
                        }
                    }
                    let star = &around[&key(ga, gb)];
                    star.iter().all(|&s| {
                        let t = cells[s as usize];
                        let vol = volume6(t.map(|v| pts[v as usize]));
                        let child = |drop: u32| -> f64 {
                            let mut q = t.map(|v| pts[v as usize]);
                            for (i, &v) in t.iter().enumerate() {
                                if v == drop {
                                    q[i] = m;
                                }
                            }
                            volume6(q)
                        };
                        let (v1, v2) = (child(gb), child(ga));
                        v1 * vol > 0.05 * vol * vol * 0.5 && v2 * vol > 0.05 * vol * vol * 0.5
                    })
                })
            };
            if mp != mid_p && !ok_all(mp) {
                mp = mid_p;
            }

            let mut parent_mid_id: Option<u32> = None;
            for &(a, b, mask, shift) in &group {
                let k = key(a, b);
                let star = &around[&k];
                let ids: Vec<u32> = facet_edges.get(&k).cloned().unwrap_or_default();
                let mut m = [mp[0] + shift[0], mp[1] + shift[1], mp[2] + shift[2]];
                let shared_bits = plane_bits[a as usize] & plane_bits[b as usize];
                for (axis, slot) in per_axis.iter().enumerate() {
                    if let Some((lo, hi)) = *slot {
                        if (shared_bits & (1 << (axis + 3))) != 0 {
                            m[axis] = lo;
                        } else if (shared_bits & (1 << axis)) != 0 {
                            m[axis] = hi;
                        }
                    }
                }
                let mv = pts.len() as u32;
                pts.push(m);
                let p_id = *parent_mid_id.get_or_insert(mv);
                plane_bits.push(shared_bits);
                parent_v.push(p_id);
                shift_v.push(shift);
                if has_periodic && shared_bits != 0 {
                    twins_of.entry(p_id).or_default().push((mv, mask, shift));
                }

                for &s in star {
                    let t = cells[s as usize];
                    let mut c1 = t;
                    let mut c2 = t;
                    for i in 0..4 {
                        if t[i] == b {
                            c1[i] = mv;
                        }
                        if t[i] == a {
                            c2[i] = mv;
                        }
                    }
                    cells[s as usize] = c1;
                    cells.push(c2);
                    labels.push(labels[s as usize]);
                }
                for f in ids {
                    let t = facets[f as usize];
                    let mut f1 = t;
                    let mut f2 = t;
                    for i in 0..3 {
                        if t[i] == b {
                            f1[i] = mv;
                        }
                        if t[i] == a {
                            f2[i] = mv;
                        }
                    }
                    facets[f as usize] = f1;
                    facets.push(f2);
                    flabels.push(flabels[f as usize]);
                }
                report.bisections += 1;
            }
        }
    }

    report.cells_after = cells.len();
    let mut out = CellMesh::<3>::new();
    out.points = pts;
    for (t, &l) in cells.iter().zip(&labels) {
        out.push_cell(CellType::Tet, t, l);
    }
    for (f, &l) in facets.iter().zip(&flabels) {
        out.push_facet(f, l);
    }
    Ok((out, report))
}
