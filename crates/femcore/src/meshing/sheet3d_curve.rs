// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Isoparametric curved mid-edge node projection and Jacobian relaxation for `sheet3d` layer-stack
//! meshes (AF-5, ADR 82).
//!
//! When `plan.opts.curved` is enabled on a 3D layer-stack model containing 2D sheets (`Op::Sheet`)
//! and curved overlay edges (`Curve::Circle` or `Curve::Ellipse`, e.g. from `Op::Circle` or
//! `Op::Polygons` with `arc_tol > 0`), the 2D PSLG triangulation `m2` places vertices on the exact
//! curves while its straight edges are chords. Extruding `m2` into 3D Dompierre prisms creates two
//! families of mesh edges along every vertical curved wall `(a, b) x [z_s, z_{s+1}]`:
//!
//! 1. Horizontal edges `[vid(a, zp), vid(b, zp)]` on `z_planes[zp]`, and
//! 2. Dompierre diagonal edges `[vid(p, f_plane), vid(q, c_plane)]` (which are
//!    `[vid(a.min(b), s), vid(a.max(b), s + 1)]` on uniform-level slabs) across each vertical quad.
//!
//! Projecting the midpoints of both the horizontal edges and the Dompierre diagonal edges onto
//! `ov.project_edge(oe, mid2)` shifts all three non-vertical edges of every vertical prism wall by
//! the same horizontal sagitta vector `\delta_{xy} = q_2 - mid_2`. On all three Dompierre tetrahedra
//! of the prism, the 3D quadratic mapping then factors into a tensor product of the 2D curved
//! triangle mapping and the affine `z` coordinate, eliminating vertical shear on thin near-sheet
//! slabs. Step-3 Jacobian relaxation (`relax_edge_nodes`) then guarantees `min_scaled_jacobian >= 0.10`
//! across all curved tetrahedra in deterministic order (`ADR 11`).

use crate::mesh::SimplexMesh;
use crate::meshing::curve3::{CurveReport, EdgeNodes};
use crate::meshing::geomdrive::MeshedGeom;
use geom::complex::{Complex, FaceGeom};
use geom::curve::Curve;
use geom::math::{dist2, V2};
use std::collections::{BTreeMap, BTreeSet, HashMap};

const MAX_OFFSET: f64 = 0.25;
const MIN_OFFSET: f64 = 1e-9;
const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

fn pair(a: u32, b: u32) -> [u32; 2] {
    [a.min(b), a.max(b)]
}

/// Borrowed transition mapping for one 2D vertical coarsening step `M_m -> M_{m+1}` in `sheet3d`.
#[derive(Clone, Copy, Debug)]
pub struct Sheet3dTransitionRef<'a> {
    pub map_v: &'a [u32],
    pub rank: &'a [(u8, usize, u8)],
}

/// Validates and relaxes candidate mid-edge geometry nodes `cand` (sorted by `[u32; 2]` vertex pair)
/// against `MIN_SCALED_JACOBIAN` (`0.10`) at all 10 Lagrange nodes and the 125 degree-6 collapsed
/// Gauss points of every curved tetrahedron in `cells`, delegating to [`crate::meshing::curve3::relax_edge_nodes`].
pub fn relax_edge_nodes(
    pos: &(dyn Fn(u32) -> [f64; 3] + Sync),
    cells: &[[u32; 4]],
    mut cand: EdgeNodes,
) -> (EdgeNodes, CurveReport) {
    let report = crate::meshing::curve3::relax_edge_nodes(pos, cells, &mut cand);
    (cand, report)
}

/// Projects midpoints of 2D boundary/interface edges `bnd_edges2` lying on `Curve::Circle` or
/// `Curve::Ellipse` onto `st.overlay`, emits 3D `EdgeNodes` on active horizontal `z_planes` and
/// active Dompierre diagonal edges `[vid(a.min(b), s), vid(a.max(b), s+1)]` (guarding
/// `vid != u32::MAX`), relaxes edge offsets to enforce `min_scaled_jacobian >= 0.10`, and applies
/// `mesh.set_edge_nodes(&nodes)`.
#[allow(clippy::too_many_arguments)]
pub fn project_sheet3d_curved_edges(
    cx: &Complex,
    m2: &SimplexMesh,
    bnd_edges2: &[(u32, u32, u32)],
    z_planes: &[f64],
    used_v: &[u32],
    slab_layer: &[Option<usize>],
    level_of_zp: &[usize],
    transitions: &[Sheet3dTransitionRef<'_>],
    mesh: &mut SimplexMesh,
) -> CurveReport {
    project_sheet3d_curved_edges_periodic(
        cx,
        m2,
        bnd_edges2,
        z_planes,
        used_v,
        slab_layer,
        level_of_zp,
        transitions,
        &[],
        mesh,
    )
}

/// Like [`project_sheet3d_curved_edges`], synchronizing curved mid-edge nodes and relaxation levels
/// across all periodic boundary pairs `(src_labels, dst_labels, shift)` (R33).
#[allow(clippy::too_many_arguments)]
pub fn project_sheet3d_curved_edges_periodic(
    cx: &Complex,
    m2: &SimplexMesh,
    bnd_edges2: &[(u32, u32, u32)],
    z_planes: &[f64],
    used_v: &[u32],
    slab_layer: &[Option<usize>],
    level_of_zp: &[usize],
    transitions: &[Sheet3dTransitionRef<'_>],
    periodic: &[(Vec<i32>, Vec<i32>, [f64; 3])],
    mesh: &mut SimplexMesh,
) -> CurveReport {
    if mesh.dim != 3 || z_planes.len() < 2 || bnd_edges2.is_empty() {
        return CurveReport::default();
    }
    let Some(st) = cx.stack.as_ref() else {
        return CurveReport::default();
    };
    let ov = &st.overlay;
    let tol = cx.tol;
    let nz = z_planes.len();
    let nslab = nz - 1;
    let n2 = m2.n_points();
    if used_v.len() < nz * n2 {
        return CurveReport::default();
    }
    let vid = |v2: u32, zp: usize| -> u32 { used_v[zp * n2 + v2 as usize] };
    let pt2 = |v2: u32| -> V2 { [m2.points[2 * v2 as usize], m2.points[2 * v2 as usize + 1]] };

    // Identify 3D lateral faces `(oe, layer_idx)` and active cap/sheet/interface `z_planes`.
    let mut lat_edges: BTreeSet<(u32, usize)> = BTreeSet::new();
    let mut cap_zps: Vec<bool> = vec![false; nz];
    for f in &cx.faces {
        match f.geom {
            FaceGeom::Cap { z, .. } => {
                if let Some(zp) = z_planes.iter().position(|&zp| (zp - z).abs() <= tol) {
                    cap_zps[zp] = true;
                }
            }
            FaceGeom::Lateral { edge, z0, z1 } => {
                for (k, l) in st.layers.iter().enumerate() {
                    if l.z0 >= z0 - tol && l.z1 <= z1 + tol {
                        lat_edges.insert((edge, k));
                    }
                }
            }
            FaceGeom::Quadric(_) | FaceGeom::Planar { .. } => {}
        }
    }
    for l in &st.layers {
        for z_bnd in [l.z0, l.z1] {
            if let Some(zp) = z_planes.iter().position(|&zp| (zp - z_bnd).abs() <= tol) {
                cap_zps[zp] = true;
            }
        }
    }

    // Collect existing tetrahedron edges so we only emit EdgeNodes on edges present in `mesh`.
    let n_cells = mesh.n_cells();
    let cells: Vec<[u32; 4]> =
        (0..n_cells).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
    let mut mesh_edges: HashMap<[u32; 2], ()> = HashMap::with_capacity(n_cells * 2);
    for c in &cells {
        for &(i, j) in EDGES.iter() {
            mesh_edges.insert(pair(c[i], c[j]), ());
        }
    }

    let lvl_at = |zp: usize| -> usize { level_of_zp.get(zp).copied().unwrap_or(0) };

    let mut cand_map: BTreeMap<[u32; 2], [f64; 3]> = BTreeMap::new();
    for &(a, b, oe) in bnd_edges2 {
        let Some(ed) = ov.edges.get(oe as usize) else { continue };
        let Some(curve) = ov.curves.get(ed.curve as usize) else { continue };
        if !matches!(curve, Curve::Circle { .. } | Curve::Ellipse { .. }) {
            continue;
        }
        let (pa, pb) = (pt2(a), pt2(b));
        let mid2 = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
        let len2 = dist2(pa, pb);
        if len2 <= 0.0 {
            continue;
        }
        let q2 = ov.project_edge(oe, mid2);
        let off2 = dist2(q2, mid2);
        if !(off2 > MIN_OFFSET * len2 && off2 < MAX_OFFSET * len2 && q2[0].is_finite() && q2[1].is_finite()) {
            continue;
        }

        // 1. Active slabs `s`: emit horizontal edges at `s` and `s + 1` AND the Dompierre diagonal edge
        // `[vid(a.min(b), s), vid(a.max(b), s + 1)]` (or transition rank-ordered diagonal) so every
        // curved prism wall has a vertical shear-free tensor-product Jacobian.
        for s in 0..nslab {
            let Some(k) = slab_layer.get(s).copied().flatten() else { continue };
            let (ls, lt) = (lvl_at(s), lvl_at(s + 1));
            let is_fine_slab = ls == 0 && lt == 0;
            let is_lat_slab = lat_edges.contains(&(oe, k));
            if !is_fine_slab && !is_lat_slab {
                continue;
            }

            let diag_spec = if ls == lt {
                let (va_s, vb_s) = (vid(a, s), vid(b, s));
                let (va_t, vb_t) = (vid(a, s + 1), vid(b, s + 1));
                if mesh_edges.contains_key(&pair(va_s, vb_t)) {
                    Some((a, b, s, s + 1))
                } else if mesh_edges.contains_key(&pair(vb_s, va_t)) {
                    Some((b, a, s, s + 1))
                } else {
                    Some((a.min(b), a.max(b), s, s + 1))
                }
            } else {
                let (m, fp, cp) = if ls < lt { (ls, s, s + 1) } else { (lt, s + 1, s) };
                transitions.get(m).and_then(|tr| {
                    let ma = tr.map_v.get(a as usize).copied().unwrap_or(a);
                    let mb = tr.map_v.get(b as usize).copied().unwrap_or(b);
                    if ma == a && mb == b {
                        let ra = tr.rank.get(a as usize).copied().unwrap_or((1, a as usize, 0));
                        let rb = tr.rank.get(b as usize).copied().unwrap_or((1, b as usize, 0));
                        let (p, q) = if ra <= rb { (a, b) } else { (b, a) };
                        Some((p, q, fp, cp))
                    } else {
                        None
                    }
                })
            };

            let Some((p, q, f_plane, c_plane)) = diag_spec else { continue };
            let (va_s, vb_s) = (vid(a, s), vid(b, s));
            let (va_t, vb_t) = (vid(a, s + 1), vid(b, s + 1));
            let (vp, vq) = (vid(p, f_plane), vid(q, c_plane));
            if [va_s, vb_s, va_t, vb_t, vp, vq].contains(&u32::MAX) || va_s == vb_s || va_t == vb_t || vp == vq {
                continue;
            }

            let e_bot = pair(va_s, vb_s);
            if mesh_edges.contains_key(&e_bot) {
                cand_map.insert(e_bot, cx.wv([q2[0], q2[1], z_planes[s]]));
            }
            let e_top = pair(va_t, vb_t);
            if mesh_edges.contains_key(&e_top) {
                cand_map.insert(e_top, cx.wv([q2[0], q2[1], z_planes[s + 1]]));
            }
            let e_diag = pair(vp, vq);
            if mesh_edges.contains_key(&e_diag) {
                let z_mid = 0.5 * (z_planes[f_plane] + z_planes[c_plane]);
                cand_map.insert(e_diag, cx.wv([q2[0], q2[1], z_mid]));
            }
        }

        // 2. Active sheet/cap/material-interface `z_planes` (when `level_of_zp[zp] == 0`).
        for (zp, &z) in z_planes.iter().enumerate() {
            if !cap_zps[zp] || lvl_at(zp) != 0 {
                continue;
            }
            let (va, vb) = (vid(a, zp), vid(b, zp));
            if va != u32::MAX && vb != u32::MAX && va != vb {
                let e_horz = pair(va, vb);
                if mesh_edges.contains_key(&e_horz) {
                    cand_map.insert(e_horz, cx.wv([q2[0], q2[1], z]));
                }
            }
        }
    }

    let mut cand: EdgeNodes = cand_map.into_iter().collect();
    let periodic_edge_parent = if !periodic.is_empty() && !cand.is_empty() {
        let mut periodic_target: Vec<[u32; 2]> = cand.iter().map(|&(e, _)| e).collect();
        let all_dst: std::collections::HashSet<i32> =
            periodic.iter().flat_map(|(_, dst, _)| dst.iter().copied()).collect();
        let mut dst_v: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for f in 0..mesh.n_facets() {
            if all_dst.contains(&mesh.facet_labels[f]) {
                dst_v.extend(mesh.facet(f));
            }
        }
        for &e in mesh_edges.keys() {
            if dst_v.contains(&e[0]) && dst_v.contains(&e[1]) {
                periodic_target.push(e);
            }
        }
        periodic_target.sort_unstable();
        periodic_target.dedup();
        let (_, pmap) = crate::meshing::curve3::periodic_edge_parent_map(mesh, &periodic_target, periodic);
        // Ensure child periodic edges always start from exact parent `q + shift` even at tolerance cusps
        let mut cmap: BTreeMap<[u32; 2], [f64; 3]> = cand.into_iter().collect();
        for (&ce, &(pe, shift)) in &pmap {
            if let Some(&pq) = cmap.get(&pe) {
                cmap.insert(ce, [pq[0] + shift[0], pq[1] + shift[1], pq[2] + shift[2]]);
            } else if let Some(&cq) = cmap.get(&ce) {
                cmap.insert(pe, [cq[0] - shift[0], cq[1] - shift[1], cq[2] - shift[2]]);
            }
        }
        cand = cmap.into_iter().collect();
        pmap
    } else {
        HashMap::new()
    };
    if cand.is_empty() {
        return CurveReport::default();
    }

    let pts = mesh.points.clone();
    let pos = |v: u32| -> [f64; 3] {
        let i = 3 * v as usize;
        [pts[i], pts[i + 1], pts[i + 2]]
    };
    let qp_edges = mesh.qp_edge_list();
    let report = crate::meshing::curve3::relax_edge_nodes_periodic_with_qp(
        &pos,
        &cells,
        &mut cand,
        &periodic_edge_parent,
        &qp_edges,
    );
    mesh.set_edge_nodes(&cand);
    if !qp_edges.is_empty() {
        mesh.set_qp_edges(&qp_edges);
    }
    report
}

/// Reconstructs the prismatic 2D column/slab topology of the physical interior of a `MeshedGeom`
/// built by `sheet3d::mesh_sheet_stack_3d` and projects + relaxes curved `EdgeNodes` strictly on
/// physical non-PML cells, preserving any pre-existing lockstep PML `edge_nodes` untouched.
pub fn apply_sheet3d_curved_edges(g: &mut MeshedGeom) -> CurveReport {
    if g.mesh.dim != 3 || g.mesh.n_points() == 0 || g.mesh.n_cells() == 0 {
        return CurveReport::default();
    }
    let cx = g.eval.complex();
    let Some(st) = cx.stack.as_ref() else {
        return CurveReport::default();
    };
    let ov = &st.overlay;
    if !ov.curves.iter().any(|c| matches!(c, Curve::Circle { .. } | Curve::Ellipse { .. })) {
        return CurveReport::default();
    }
    let pml_doms: std::collections::HashSet<i32> = g
        .eval
        .pml_regions()
        .iter()
        .map(|p| g.labels.domains[p.region as usize].id)
        .collect();
    let has_pml_cells = !pml_doms.is_empty()
        && (0..g.mesh.n_cells()).any(|c| pml_doms.contains(&g.mesh.cell_labels[c]));

    let tol = cx.tol.max(1e-12);
    let n3 = g.mesh.n_points();
    let mut phys_vert = vec![!has_pml_cells; n3];
    let mut pml_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
    if has_pml_cells {
        for c in 0..g.mesh.n_cells() {
            let tet = g.mesh.cell(c);
            if pml_doms.contains(&g.mesh.cell_labels[c]) {
                for &(i, j) in &EDGES {
                    pml_edges.insert(pair(tet[i], tet[j]));
                }
            } else {
                for &v in tet {
                    phys_vert[v as usize] = true;
                }
            }
        }
    }
    if !phys_vert.iter().any(|&b| b) {
        return CurveReport::default();
    }

    // 1. Extract unique `z_planes` and 2D `(x, y)` columns in internal stack frame (`cx.iv`) across physical vertices.
    let mut stack_pts: Vec<[f64; 3]> = Vec::with_capacity(n3);
    let mut z_planes: Vec<f64> = Vec::new();
    for v3 in 0..n3 {
        let pw = [g.mesh.points[3 * v3], g.mesh.points[3 * v3 + 1], g.mesh.points[3 * v3 + 2]];
        let ps = cx.iv(pw);
        stack_pts.push(ps);
        if phys_vert[v3] && !z_planes.iter().any(|&z| (z - ps[2]).abs() <= tol) {
            z_planes.push(ps[2]);
        }
    }
    z_planes.sort_by(|a, b| a.total_cmp(b));
    let nz = z_planes.len();
    if nz < 2 {
        return CurveReport::default();
    }

    let grid_h = (16.0 * tol).max(1e-11);
    let mut col_grid: HashMap<(i64, i64), Vec<u32>> = HashMap::new();
    let mut pts2: Vec<V2> = Vec::new();
    let mut v2_of_v3 = vec![u32::MAX; n3];
    let mut zp_of_v3 = vec![0usize; n3];

    for (v3, ps) in stack_pts.iter().enumerate() {
        if !phys_vert[v3] {
            continue;
        }
        let p2 = [ps[0], ps[1]];
        let gx = (p2[0] / grid_h).floor() as i64;
        let gy = (p2[1] / grid_h).floor() as i64;
        let mut found = None;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(bucket) = col_grid.get(&(gx + dx, gy + dy)) {
                    for &c2 in bucket {
                        if dist2(pts2[c2 as usize], p2) <= 4.0 * tol {
                            found = Some(c2);
                            break 'search;
                        }
                    }
                }
            }
        }
        let v2 = match found {
            Some(c2) => c2,
            None => {
                let id = pts2.len() as u32;
                pts2.push(p2);
                col_grid.entry((gx, gy)).or_default().push(id);
                id
            }
        };
        v2_of_v3[v3] = v2;
        let zp = z_planes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (*a - ps[2]).abs().total_cmp(&(*b - ps[2]).abs()))
            .map(|(i, _)| i)
            .unwrap_or(0);
        zp_of_v3[v3] = zp;
    }

    let n2 = pts2.len();
    let mut used_v = vec![u32::MAX; nz * n2];
    for v3 in 0..n3 {
        if phys_vert[v3] && v2_of_v3[v3] != u32::MAX {
            used_v[zp_of_v3[v3] * n2 + v2_of_v3[v3] as usize] = v3 as u32;
        }
    }

    let (zp_fine, max_active_per_plane) = (0..nz)
        .map(|zp| (zp, (0..n2).filter(|&v2| used_v[zp * n2 + v2] != u32::MAX).count()))
        .max_by_key(|&(_, cnt)| cnt)
        .unwrap_or((0, 0));

    // 2. Collect the exact level-0 2D triangulation `m2` from physical tetrahedron faces lying on `zp_fine`.
    let mut tri_set: BTreeSet<[u32; 3]> = BTreeSet::new();
    let mut tris2: Vec<([u32; 3], usize)> = Vec::new();
    for c in 0..g.mesh.n_cells() {
        if has_pml_cells && pml_doms.contains(&g.mesh.cell_labels[c]) {
            continue;
        }
        let tet = g.mesh.cell(c);
        let mut fine_cols = [0u32; 4];
        let mut nf = 0usize;
        for &v3 in tet {
            if zp_of_v3[v3 as usize] == zp_fine {
                fine_cols[nf] = v2_of_v3[v3 as usize];
                nf += 1;
            }
        }
        if nf == 3 {
            let mut key = [fine_cols[0], fine_cols[1], fine_cols[2]];
            key.sort_unstable();
            if key[0] != key[1] && key[1] != key[2] && tri_set.insert(key) {
                let (p0, p1, p2) = (pts2[key[0] as usize], pts2[key[1] as usize], pts2[key[2] as usize]);
                let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0];
                if let Some(ov_c) = ov.locate(cen) {
                    tris2.push((key, ov_c as usize));
                }
            }
        }
    }

    let mut m2_pts = Vec::with_capacity(2 * n2);
    for p in &pts2 {
        m2_pts.push(p[0]);
        m2_pts.push(p[1]);
    }
    let mut m2_cells = Vec::with_capacity(3 * tris2.len());
    let mut m2_labels = Vec::with_capacity(tris2.len());
    for &(v, ov_c) in &tris2 {
        m2_cells.extend_from_slice(&v);
        m2_labels.push(ov_c as i32 + 1);
    }
    let m2 = SimplexMesh::new(2, m2_pts, m2_cells, m2_labels, Vec::new(), Vec::new());

    // 3. Identify 2D boundary/interface edges `bnd_edges2` matching `sheet3d.rs:888-925`.
    let mut edge2_sides: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for &(v, ov_c) in &tris2 {
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            edge2_sides.entry((a.min(b), a.max(b))).or_default().push(ov_c);
        }
    }
    let mut bnd_edges2: Vec<(u32, u32, u32)> = Vec::new();
    let mut edge_list: Vec<((u32, u32), Vec<usize>)> = edge2_sides.into_iter().collect();
    edge_list.sort_by_key(|x| x.0);
    for &((a, b), ref sides) in &edge_list {
        let (c0, c1) = match sides.as_slice() {
            [c0] => (Some(*c0 as u32), None),
            [c0, c1] if c0 != c1 => (Some(*c0 as u32), Some(*c1 as u32)),
            _ => continue,
        };
        let (pa, pb) = (pts2[a as usize], pts2[b as usize]);
        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
        let mut best: Option<(f64, u32)> = None;
        for pass in 0..2 {
            for (ei, ed) in ov.edges.iter().enumerate() {
                if pass == 0 && !((ed.left == c0 && ed.right == c1) || (ed.left == c1 && ed.right == c0)) {
                    continue;
                }
                let d = dist2(ov.project_edge(ei as u32, mid), mid);
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, ei as u32));
                }
            }
            if best.is_some() {
                break;
            }
        }
        if let Some((_, oe)) = best {
            bnd_edges2.push((a, b, oe));
        }
    }

    // 4. Map each slab `[z_planes[s], z_planes[s+1]]` to its 3D stack layer and detect fine level-0 planes.
    let nslab = nz - 1;
    let mut slab_layer: Vec<Option<usize>> = Vec::with_capacity(nslab);
    for s in 0..nslab {
        let zm = 0.5 * (z_planes[s] + z_planes[s + 1]);
        let layer_idx = st.layers.iter().position(|l| zm >= l.z0 - tol && zm <= l.z1 + tol);
        slab_layer.push(layer_idx);
    }

    let level_of_zp: Vec<usize> = (0..nz)
        .map(|zp| {
            let cnt = (0..n2).filter(|&v2| used_v[zp * n2 + v2] != u32::MAX).count();
            if cnt == max_active_per_plane {
                0
            } else {
                1
            }
        })
        .collect();

    let periodic_twins = crate::meshing::curve3::periodic_pairs_from_eval(&g.eval, &g.labels);
    if !has_pml_cells {
        let qp_backup = g.mesh.qp_edge_list();
        let rep = project_sheet3d_curved_edges_periodic(
            cx,
            &m2,
            &bnd_edges2,
            &z_planes,
            &used_v,
            &slab_layer,
            &level_of_zp,
            &[],
            &periodic_twins,
            &mut g.mesh,
        );
        if !qp_backup.is_empty() {
            g.mesh.set_qp_edges(&qp_backup);
        }
        return rep;
    }

    // Build a temporary physical-only SimplexMesh view sharing `g.mesh.points` so relaxation only inspects physical tets!
    let existing_en = g.mesh.edge_node_list();
    let existing_qp = g.mesh.qp_edge_list();
    let mut phys_cells_flat = Vec::new();
    let mut phys_labels = Vec::new();
    for c in 0..g.mesh.n_cells() {
        if !pml_doms.contains(&g.mesh.cell_labels[c]) {
            phys_cells_flat.extend_from_slice(g.mesh.cell(c));
            phys_labels.push(g.mesh.cell_labels[c]);
        }
    }
    let mut phys_facets_flat = Vec::new();
    let mut phys_facet_labels = Vec::new();
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        if fv.iter().all(|&v| phys_vert[v as usize]) {
            phys_facets_flat.extend_from_slice(fv);
            phys_facet_labels.push(g.mesh.facet_labels[f]);
        }
    }
    let mut tmp_phys = SimplexMesh::new(3, g.mesh.points.clone(), phys_cells_flat, phys_labels, phys_facets_flat, phys_facet_labels);
    if !existing_qp.is_empty() {
        let kept_qp: Vec<[u32; 2]> = existing_qp
            .iter()
            .copied()
            .filter(|&[a, b]| !pml_edges.contains(&pair(a, b)))
            .collect();
        tmp_phys.set_qp_edges(&kept_qp);
    }
    let rep = project_sheet3d_curved_edges_periodic(
        cx,
        &m2,
        &bnd_edges2,
        &z_planes,
        &used_v,
        &slab_layer,
        &level_of_zp,
        &[],
        &periodic_twins,
        &mut tmp_phys,
    );
    let mut merged: BTreeMap<[u32; 2], [f64; 3]> = tmp_phys.edge_node_list().into_iter().collect();
    for (pr, q) in existing_en {
        if pml_edges.contains(&pr) {
            merged.insert(pr, q);
        }
    }
    let merged_vec: Vec<([u32; 2], [f64; 3])> = merged.into_iter().collect();
    g.mesh.set_edge_nodes(&merged_vec);
    if !existing_qp.is_empty() {
        g.mesh.set_qp_edges(&existing_qp);
    }
    rep
}
