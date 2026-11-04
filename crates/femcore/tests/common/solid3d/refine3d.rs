//! Local refinement of a tetrahedral mesh by longest-edge (Rivara) bisection.
//!
//! A cell whose longest edge is larger than the wanted element size at its centroid is marked. For a
//! marked cell the algorithm follows the longest edge to a *terminal* edge (one that is the longest
//! edge of every cell around it) and splits all cells around that edge in two. Because every cell
//! that contains the edge is split along it, the mesh stays conforming, and because terminal edges
//! are the longest ones the cells keep reasonable shapes. Boundary facets are split together with the
//! cells, and the new vertex on an edge that lies on the boundary of the solid or on an interface
//! between domains is moved onto the true (curved) surface by the caller's projection.

use femcore::mesh::SimplexMesh;
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

/// The longest edge of a cell under a strict total order (length, then vertex numbers), so two cells
/// that share an edge always agree on which of two edges is longer.
fn longest(cell: &[u32; 4], pts: &[P3]) -> (u32, u32) {
    let mut best = (0u64, 0u64);
    let mut edge = (cell[0], cell[1]);
    for &(i, j) in &PAIRS {
        let (a, b) = (cell[i], cell[j]);
        let l = (dist2(pts[a as usize], pts[b as usize]).to_bits() >> 8, key(a, b));
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
    mesh: &SimplexMesh,
    size: &(dyn Fn(P3) -> f64 + Sync),
    project: &(dyn Fn(P3, P3, P3, EdgeKind) -> P3 + Sync),
    ratio: f64,
    max_cells: usize,
) -> Result<(SimplexMesh, RefineReport), String> {
    if mesh.dim != 3 {
        return Err("local refinement works on tetrahedral meshes".into());
    }
    let mut pts: Vec<P3> = (0..mesh.n_points()).map(|i| [mesh.points[3 * i], mesh.points[3 * i + 1], mesh.points[3 * i + 2]]).collect();
    let mut cells: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cells[4 * c], mesh.cells[4 * c + 1], mesh.cells[4 * c + 2], mesh.cells[4 * c + 3]]).collect();
    let mut labels = mesh.cell_labels.clone();
    let mut facets: Vec<[u32; 3]> = (0..mesh.n_facets()).map(|f| [mesh.facets[3 * f], mesh.facets[3 * f + 1], mesh.facets[3 * f + 2]]).collect();
    let mut flabels = mesh.facet_labels.clone();
    let mut report = RefineReport { cells_before: cells.len(), ..Default::default() };

    for _round in 0..200 {
        // 1. cells that are too large
        let marked: Vec<usize> = (0..cells.len())
            .filter(|&c| {
                let (a, b) = longest(&cells[c], &pts);
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

        // 3. terminal edges whose stars do not overlap
        let mut claimed = vec![false; cells.len()];
        let mut selected: Vec<(u32, u32)> = vec![];
        for &c0 in &marked {
            if claimed[c0] {
                continue;
            }
            let mut c = c0;
            let mut steps = 0;
            let edge = loop {
                let e = longest(&cells[c], &pts);
                let star = &around[&key(e.0, e.1)];
                let next = star.iter().find(|&&s| {
                    let le = longest(&cells[s as usize], &pts);
                    key(le.0, le.1) != key(e.0, e.1)
                });
                match next {
                    None => break e,
                    Some(&s) => c = s as usize,
                }
                steps += 1;
                if steps > 100_000 {
                    return Err("local refinement did not find a terminal edge".into());
                }
            };
            let star = &around[&key(edge.0, edge.1)];
            if star.iter().any(|&s| claimed[s as usize]) {
                continue;
            }
            for &s in star {
                claimed[s as usize] = true;
            }
            selected.push(edge);
        }

        // 4. split them
        for (a, b) in selected {
            let k = key(a, b);
            let star = &around[&k];
            let ids: Vec<u32> = facet_edges.get(&k).cloned().unwrap_or_default();
            let kind = if !ids.is_empty() {
                EdgeKind::Boundary
            } else if star.iter().any(|&s| labels[s as usize] != labels[star[0] as usize]) {
                EdgeKind::Interface
            } else {
                EdgeKind::Interior
            };
            let (pa, pb) = (pts[a as usize], pts[b as usize]);
            let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0];
            let mut m = if kind == EdgeKind::Interior { mid } else { project(mid, pa, pb, kind) };
            // the projection must keep every new cell positively oriented and not squash it
            let ok = |m: P3| -> bool {
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
                    let (v1, v2) = (child(b), child(a));
                    v1 * vol > 0.05 * vol * vol * 0.5 && v2 * vol > 0.05 * vol * vol * 0.5
                })
            };
            if m != mid && !ok(m) {
                m = mid;
            }
            let mv = pts.len() as u32;
            pts.push(m);
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

    report.cells_after = cells.len();
    let points: Vec<f64> = pts.iter().flatten().copied().collect();
    let out = SimplexMesh::new(3, points, cells.iter().flatten().copied().collect(), labels, facets.iter().flatten().copied().collect(), flabels);
    Ok((out, report))
}
