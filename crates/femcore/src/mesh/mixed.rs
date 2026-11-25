// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Native mixed-cell solver discretisation (`MixedDiscretisation`) built directly from `CellMesh<D>`
//! without splitting `Quad`, `Prism`, or `Hex` cells into simplices.
//!
//! Preserves exact cell and boundary-facet CSR arrays alongside deterministic sorted unique edges,
//! triangular faces (`D_3` permutations), quadrilateral faces (`D_4` dihedral symmetries), curved
//! mid-edge and quad-face center geometry nodes in both 2D and 3D, and conformity verification.

use crate::mesh::cells::{canonical_quad_face, CellMesh, CellType};
use crate::mesh::SimplexMesh;
use std::collections::{BTreeSet, HashMap};

/// Boundary or interior sheet facet shapes supported in mixed meshes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum FacetType {
    /// 2-node line segment (boundary or interior interface in 2D).
    Edge2,
    /// 3-node triangle facet (boundary or interior sheet in 3D).
    Tri3,
    /// 4-node quadrilateral facet in cyclic order (boundary or interior sheet in 3D).
    Quad4,
}

impl FacetType {
    /// Number of corner vertices of the facet.
    pub fn nv(self) -> usize {
        match self {
            FacetType::Edge2 => 2,
            FacetType::Tri3 => 3,
            FacetType::Quad4 => 4,
        }
    }

    /// Ambient mesh dimension in which this facet serves as a codimension-1 facet (`2` or `3`).
    pub fn mesh_dim(self) -> usize {
        match self {
            FacetType::Edge2 => 2,
            FacetType::Tri3 | FacetType::Quad4 => 3,
        }
    }
}

/// Canonical local edges of `CellType::Tri` (matches `SimplexMesh` / `fem.rs` / `edge.rs`).
pub const TRI_EDGES: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

/// Canonical local edges of `CellType::Quad` in counter-clockwise boundary order.
pub const QUAD_EDGES: [(usize, usize); 4] = [(0, 1), (1, 2), (2, 3), (3, 0)];

/// Canonical local edges of `CellType::Tet` (matches `fem::TET_EDGES` and `nyx_ems::cavity`).
pub const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// Canonical local edges of `CellType::Prism`: 3 base-tri edges, 3 top-tri edges, 3 vertical edges.
pub const PRISM_EDGES: [(usize, usize); 9] = [
    (0, 1),
    (0, 2),
    (1, 2),
    (3, 4),
    (3, 5),
    (4, 5),
    (0, 3),
    (1, 4),
    (2, 5),
];

/// Canonical local edges of `CellType::Hex`: 4 bottom-quad edges, 4 top-quad edges, 4 vertical edges.
pub const HEX_EDGES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

/// Canonical local triangular faces of `CellType::Tet` (matches `fem::TET_FACES` and `cavity::FACES`).
pub const TET_TRI_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Canonical local triangular faces of `CellType::Prism` (base `z=-1`, top `z=+1`).
pub const PRISM_TRI_FACES: [[usize; 3]; 2] = [[0, 1, 2], [3, 4, 5]];

/// Canonical local quadrilateral faces of `CellType::Prism` in outward cyclic order (`cells.rs`).
pub const PRISM_QUAD_FACES: [[usize; 4]; 3] = [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];

/// Canonical local quadrilateral faces of `CellType::Hex` in outward cyclic order (`cells.rs`).
pub const HEX_QUAD_FACES: [[usize; 4]; 6] = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
];

/// All 6 permutations of `[0, 1, 2]` (symmetric group `D_3` / `S_3`).
pub const D3_PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// Local edge endpoint table for any [`CellType`].
#[inline]
pub fn cell_local_edges(t: CellType) -> &'static [(usize, usize)] {
    match t {
        CellType::Tri => &TRI_EDGES,
        CellType::Quad => &QUAD_EDGES,
        CellType::Tet => &TET_EDGES,
        CellType::Prism => &PRISM_EDGES,
        CellType::Hex => &HEX_EDGES,
    }
}

/// Local triangular faces for a 3D [`CellType`] (empty in 2D and for `Hex`).
#[inline]
pub fn cell_local_tri_faces(t: CellType) -> &'static [[usize; 3]] {
    match t {
        CellType::Tet => &TET_TRI_FACES,
        CellType::Prism => &PRISM_TRI_FACES,
        _ => &[],
    }
}

/// Local quadrilateral faces for a 3D [`CellType`] (empty in 2D and for `Tet`).
#[inline]
pub fn cell_local_quad_faces(t: CellType) -> &'static [[usize; 4]] {
    match t {
        CellType::Prism => &PRISM_QUAD_FACES,
        CellType::Hex => &HEX_QUAD_FACES,
        _ => &[],
    }
}

/// Sorted vertex triple `[a, b, c]` (`a <= b <= c`).
#[inline]
pub fn sorted_tri_face(mut f: [u32; 3]) -> [u32; 3] {
    if f[0] > f[1] {
        f.swap(0, 1);
    }
    if f[1] > f[2] {
        f.swap(1, 2);
    }
    if f[0] > f[1] {
        f.swap(0, 1);
    }
    f
}

/// Index `p` in `0..6` into [`D3_PERMS`] such that `keys[perm[0]] <= keys[perm[1]] <= keys[perm[2]]`.
#[inline]
pub fn d3_perm_for_keys<K: Ord + Copy>(keys: [K; 3]) -> u8 {
    for (idx, p) in D3_PERMS.iter().enumerate() {
        if keys[p[0]] <= keys[p[1]] && keys[p[1]] <= keys[p[2]] {
            return idx as u8;
        }
    }
    0
}

/// Applies `D_3` permutation `perm` (`0..6`) to a 3-element array.
#[inline]
pub fn d3_apply<T: Copy>(perm: u8, vals: [T; 3]) -> [T; 3] {
    let p = D3_PERMS[perm as usize];
    [vals[p[0]], vals[p[1]], vals[p[2]]]
}

/// Computes the `D_4` dihedral symmetry index `d` in `0..8` (`rot = d & 3`, `flip = (d & 4) != 0`)
/// that maps cyclic quad keys `keys` to canonical cyclic order (starting at minimum key, stepping
/// toward the smaller of its two cyclic neighbors).
#[inline]
pub fn d4_dihedral_for_keys<K: Ord + Copy>(keys: [K; 4]) -> u8 {
    let rot = (0..4).min_by_key(|&i| keys[i]).unwrap();
    let next = keys[(rot + 1) % 4];
    let prev = keys[(rot + 3) % 4];
    let flip = prev < next;
    (rot as u8) | ((flip as u8) << 2)
}

/// Applies `D_4` dihedral index `d` (`0..8`) to cyclic quad values `cyc`.
#[inline]
pub fn d4_apply<T: Copy>(d: u8, cyc: [T; 4]) -> [T; 4] {
    let rot = (d & 3) as usize;
    let flip = (d & 4) != 0;
    if !flip {
        [cyc[rot], cyc[(rot + 1) % 4], cyc[(rot + 2) % 4], cyc[(rot + 3) % 4]]
    } else {
        [cyc[rot], cyc[(rot + 3) % 4], cyc[(rot + 2) % 4], cyc[(rot + 1) % 4]]
    }
}

/// Canonical cyclic quad `[a, b, c, d]` ordered by `keys`.
#[inline]
pub fn canonical_quad_by_keys<K: Ord + Copy, T: Copy>(keys: [K; 4], vals: [T; 4]) -> (u8, [T; 4]) {
    let d = d4_dihedral_for_keys(keys);
    (d, d4_apply(d, vals))
}

/// Native mixed-cell discretisation for 2D (`Tri` + `Quad`) and 3D (`Tet` + `Prism` + `Hex`) FEM assembly.
#[derive(Clone, Debug)]
pub struct MixedDiscretisation {
    /// Spatial dimension (`2` or `3`).
    pub dim: usize,
    /// Flat vertex coordinates (`dim` values per vertex).
    pub points: Vec<f64>,
    /// Cell shape per cell (`len = n_cells`).
    pub cell_types: Vec<CellType>,
    /// CSR offsets into `cell_conn` (`len = n_cells + 1`).
    pub cell_offsets: Vec<u32>,
    /// Flat cell vertex connectivity.
    pub cell_conn: Vec<u32>,
    /// 1-based domain label per cell (`len = n_cells`).
    pub cell_labels: Vec<i32>,
    /// Facet shape per labelled boundary/sheet facet (`len = n_facets`).
    pub facet_types: Vec<FacetType>,
    /// CSR offsets into `facet_conn` (`len = n_facets + 1`).
    pub facet_offsets: Vec<u32>,
    /// Flat facet vertex connectivity.
    pub facet_conn: Vec<u32>,
    /// 1-based boundary label per facet (`len = n_facets`).
    pub facet_labels: Vec<i32>,
    /// Whether facet `f` lies on the exterior boundary of the mesh (belongs to exactly 1 cell).
    pub exterior_facets: Vec<bool>,
    /// Deterministic lexicographically sorted unique edges `(a, b)` with `a < b`.
    pub edges: Vec<(u32, u32)>,
    /// CSR offsets into `cell_edges` and `cell_edge_signs` (`len = n_cells + 1`).
    pub cell_edge_offsets: Vec<u32>,
    /// Global edge ID for each local edge of cell `c` (in [`cell_local_edges`] order).
    pub cell_edges: Vec<u32>,
    /// Default orientation sign (`+1` if `v[la] < v[lb]`, `-1` otherwise) for each local edge of cell `c`.
    pub cell_edge_signs: Vec<i8>,
    /// Deterministic lexicographically sorted unique 3D triangular faces `[a, b, c]` (`a < b < c`).
    pub tri_faces: Vec<[u32; 3]>,
    /// CSR offsets into `cell_tri_faces` and `cell_tri_perms` (`len = n_cells + 1`).
    pub cell_tri_offsets: Vec<u32>,
    /// Global tri-face ID for each local tri face of cell `c` (in [`cell_local_tri_faces`] order).
    pub cell_tri_faces: Vec<u32>,
    /// `D_3` permutation index (`0..6` into [`D3_PERMS`]) sorting each local tri face by vertex ID.
    pub cell_tri_perms: Vec<u8>,
    /// Deterministic lexicographically sorted unique 3D quadrilateral faces in [`canonical_quad_face`] order.
    pub quad_faces: Vec<[u32; 4]>,
    /// CSR offsets into `cell_quad_faces` and `cell_quad_dihedrals` (`len = n_cells + 1`).
    pub cell_quad_offsets: Vec<u32>,
    /// Global quad-face ID for each local quad face of cell `c` (in [`cell_local_quad_faces`] order).
    pub cell_quad_faces: Vec<u32>,
    /// `D_4` dihedral index (`0..8`, see [`d4_apply`]) canonicalising each local quad face by vertex ID.
    pub cell_quad_dihedrals: Vec<u8>,
    /// Curved mid-edge geometry nodes `(edge_id, [x, y, z])` sorted by `edge_id` (in 2D, 3rd entry is `0.0`).
    pub edge_nodes: Vec<(u32, [f64; 3])>,
    /// Curved quad-face center nodes `(id, [x, y, z])` sorted by `id` (`quad_face_id` in 3D, `cell_id` in 2D).
    pub quad_face_nodes: Vec<(u32, [f64; 3])>,
    /// Sorted global edge IDs carrying Barsoum quarter-point singularity shifts.
    pub qp_edges: Vec<u32>,
}

impl MixedDiscretisation {
    /// Builds a [`MixedDiscretisation`] directly from a [`CellMesh<D>`] without splitting any cells.
    pub fn from_cell_mesh<const D: usize>(cm: &CellMesh<D>) -> Result<Self, String> {
        if D != 2 && D != 3 {
            return Err(format!("mesh dimension {D} is not supported"));
        }
        let nc = cm.n_cells();
        let np = cm.n_points();
        let nf = cm.n_facets();
        if nc == 0 {
            return Err("the mesh has no cells".into());
        }

        let mut points = Vec::with_capacity(np * D);
        for p in &cm.points {
            points.extend_from_slice(p);
        }

        let cell_types = cm.types.clone();
        let cell_offsets = if cm.offsets.is_empty() { vec![0] } else { cm.offsets.clone() };
        let cell_conn = cm.conn.clone();
        let cell_labels = cm.label.clone();

        for c in 0..nc {
            let t = cell_types[c];
            if t.dim() != D {
                return Err(format!("cell {c}: shape {t:?} has dimension {}, expected {D}", t.dim()));
            }
            let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
            if vs.len() != t.nv() {
                return Err(format!("cell {c}: expected {} vertices for {t:?}, got {}", t.nv(), vs.len()));
            }
            for &v in vs {
                if v as usize >= np {
                    return Err(format!("cell {c}: vertex {v} out of range (n_points = {np})"));
                }
            }
        }

        let mut facet_types = Vec::with_capacity(nf);
        let facet_offsets = if cm.facet_offsets.is_empty() { vec![0] } else { cm.facet_offsets.clone() };
        let facet_conn = cm.facet_conn.clone();
        let facet_labels = cm.facet_label.clone();

        for f in 0..nf {
            let vs = &facet_conn[facet_offsets[f] as usize..facet_offsets[f + 1] as usize];
            let ft = match (D, vs.len()) {
                (2, 2) => FacetType::Edge2,
                (3, 3) => FacetType::Tri3,
                (3, 4) => FacetType::Quad4,
                _ => {
                    return Err(format!(
                        "facet {f}: unsupported vertex count {} in {D}D mesh",
                        vs.len()
                    ));
                }
            };
            for &v in vs {
                if v as usize >= np {
                    return Err(format!("facet {f}: vertex {v} out of range (n_points = {np})"));
                }
            }
            facet_types.push(ft);
        }

        // 1. Collect & deduplicate edges deterministically.
        let mut raw_edges: Vec<(u32, u32)> = Vec::with_capacity(nc * 6);
        let mut cell_edge_offsets = Vec::with_capacity(nc + 1);
        cell_edge_offsets.push(0);
        for c in 0..nc {
            let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
            let le = cell_local_edges(cell_types[c]);
            for &(i, j) in le {
                let (a, b) = (vs[i], vs[j]);
                if a == b {
                    return Err(format!("cell {c}: degenerate zero-length edge on repeated vertex {a}"));
                }
                raw_edges.push((a.min(b), a.max(b)));
            }
            cell_edge_offsets.push(raw_edges.len() as u32);
        }
        raw_edges.sort_unstable();
        raw_edges.dedup();
        let edges = raw_edges;

        let total_cell_edges = *cell_edge_offsets.last().unwrap() as usize;
        let mut cell_edges = Vec::with_capacity(total_cell_edges);
        let mut cell_edge_signs = Vec::with_capacity(total_cell_edges);
        let mut edge_cell_count = if D == 2 { vec![0u32; edges.len()] } else { Vec::new() };
        for c in 0..nc {
            let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
            for &(i, j) in cell_local_edges(cell_types[c]) {
                let (va, vb) = (vs[i], vs[j]);
                let key = (va.min(vb), va.max(vb));
                let eid = edges.binary_search(&key).expect("edge must exist") as u32;
                cell_edges.push(eid);
                cell_edge_signs.push(if va < vb { 1 } else { -1 });
                if D == 2 {
                    edge_cell_count[eid as usize] += 1;
                }
            }
        }

        // 2. Collect & deduplicate 3D tri_faces and quad_faces deterministically.
        let mut tri_faces: Vec<[u32; 3]> = Vec::new();
        let mut cell_tri_offsets = Vec::with_capacity(nc + 1);
        cell_tri_offsets.push(0);
        let mut cell_tri_faces: Vec<u32> = Vec::new();
        let mut cell_tri_perms: Vec<u8> = Vec::new();

        let mut quad_faces: Vec<[u32; 4]> = Vec::new();
        let mut cell_quad_offsets = Vec::with_capacity(nc + 1);
        cell_quad_offsets.push(0);
        let mut cell_quad_faces: Vec<u32> = Vec::new();
        let mut cell_quad_dihedrals: Vec<u8> = Vec::new();

        let mut tri_cell_count: Vec<u32> = Vec::new();
        let mut quad_cell_count: Vec<u32> = Vec::new();

        if D == 3 {
            let mut raw_tris: Vec<[u32; 3]> = Vec::with_capacity(nc * 2);
            let mut raw_quads: Vec<[u32; 4]> = Vec::with_capacity(nc * 3);
            for c in 0..nc {
                let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
                let tfs = cell_local_tri_faces(cell_types[c]);
                for lf in tfs {
                    raw_tris.push(sorted_tri_face([vs[lf[0]], vs[lf[1]], vs[lf[2]]]));
                }
                cell_tri_offsets.push((cell_tri_offsets[c] as usize + tfs.len()) as u32);

                let qfs = cell_local_quad_faces(cell_types[c]);
                for lq in qfs {
                    raw_quads.push(canonical_quad_face([vs[lq[0]], vs[lq[1]], vs[lq[2]], vs[lq[3]]]));
                }
                cell_quad_offsets.push((cell_quad_offsets[c] as usize + qfs.len()) as u32);
            }
            raw_tris.sort_unstable();
            raw_tris.dedup();
            tri_faces = raw_tris;

            raw_quads.sort_unstable();
            raw_quads.dedup();
            quad_faces = raw_quads;

            tri_cell_count = vec![0u32; tri_faces.len()];
            quad_cell_count = vec![0u32; quad_faces.len()];

            let total_tris = *cell_tri_offsets.last().unwrap() as usize;
            let total_quads = *cell_quad_offsets.last().unwrap() as usize;
            cell_tri_faces.reserve(total_tris);
            cell_tri_perms.reserve(total_tris);
            cell_quad_faces.reserve(total_quads);
            cell_quad_dihedrals.reserve(total_quads);

            for c in 0..nc {
                let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
                for lf in cell_local_tri_faces(cell_types[c]) {
                    let f = [vs[lf[0]], vs[lf[1]], vs[lf[2]]];
                    let perm = d3_perm_for_keys(f);
                    let key = d3_apply(perm, f);
                    let tf_id = tri_faces.binary_search(&key).expect("tri face must exist") as u32;
                    cell_tri_faces.push(tf_id);
                    cell_tri_perms.push(perm);
                    tri_cell_count[tf_id as usize] += 1;
                }
                for lq in cell_local_quad_faces(cell_types[c]) {
                    let q = [vs[lq[0]], vs[lq[1]], vs[lq[2]], vs[lq[3]]];
                    let (dih, key) = canonical_quad_by_keys(q, q);
                    let qf_id = quad_faces.binary_search(&key).expect("quad face must exist") as u32;
                    cell_quad_faces.push(qf_id);
                    cell_quad_dihedrals.push(dih);
                    quad_cell_count[qf_id as usize] += 1;
                }
            }
        } else {
            cell_tri_offsets.resize(nc + 1, 0);
            cell_quad_offsets.resize(nc + 1, 0);
        }

        // 3. Exterior facet flags.
        let mut exterior_facets = Vec::with_capacity(nf);
        for f in 0..nf {
            let vs = &facet_conn[facet_offsets[f] as usize..facet_offsets[f + 1] as usize];
            let is_ext = match facet_types[f] {
                FacetType::Edge2 => {
                    let key = (vs[0].min(vs[1]), vs[0].max(vs[1]));
                    edges
                        .binary_search(&key)
                        .ok()
                        .is_some_and(|eid| edge_cell_count[eid] == 1)
                }
                FacetType::Tri3 => {
                    let key = sorted_tri_face([vs[0], vs[1], vs[2]]);
                    tri_faces
                        .binary_search(&key)
                        .ok()
                        .is_some_and(|tf| tri_cell_count[tf] == 1)
                }
                FacetType::Quad4 => {
                    let key = canonical_quad_face([vs[0], vs[1], vs[2], vs[3]]);
                    quad_faces
                        .binary_search(&key)
                        .ok()
                        .is_some_and(|qf| quad_cell_count[qf] == 1)
                }
            };
            exterior_facets.push(is_ext);
        }

        // 4. Curved mid-edge nodes (both 2D and 3D), quad-face center nodes, and quarter-point edges.
        let mut edge_nodes: Vec<(u32, [f64; 3])> = cm
            .edge_nodes
            .iter()
            .filter_map(|&([a, b], pos)| {
                let key = (a.min(b), a.max(b));
                let eid = edges.binary_search(&key).ok()? as u32;
                let mut p3 = [0.0; 3];
                p3[..D].copy_from_slice(&pos);
                Some((eid, p3))
            })
            .collect();
        edge_nodes.sort_by_key(|e| e.0);
        edge_nodes.dedup_by_key(|e| e.0);

        let mut qp_edges: Vec<u32> = cm
            .qp_edges
            .iter()
            .filter_map(|&[a, b]| {
                let key = (a.min(b), a.max(b));
                let eid = edges.binary_search(&key).ok()? as u32;
                edge_nodes.binary_search_by_key(&eid, |x| x.0).ok().map(|_| eid)
            })
            .collect();
        qp_edges.sort_unstable();
        qp_edges.dedup();

        let mut quad_face_nodes: Vec<(u32, [f64; 3])> = Vec::new();
        if !cm.quad_face_nodes.is_empty() {
            if D == 3 {
                for &(q, pos) in &cm.quad_face_nodes {
                    let key = canonical_quad_face(q);
                    if let Ok(qf_id) = quad_faces.binary_search(&key) {
                        quad_face_nodes.push((qf_id as u32, [pos[0], pos[1], pos[2]]));
                    }
                }
            } else {
                // In 2D, map canonical 4-vertex quad keys to the owning `CellType::Quad` cell index.
                let mut quad_cell_map: HashMap<[u32; 4], u32> = HashMap::new();
                for c in 0..nc {
                    if cell_types[c] == CellType::Quad {
                        let vs = &cell_conn[cell_offsets[c] as usize..cell_offsets[c + 1] as usize];
                        quad_cell_map.insert(canonical_quad_face([vs[0], vs[1], vs[2], vs[3]]), c as u32);
                    }
                }
                for &(q, pos) in &cm.quad_face_nodes {
                    let key = canonical_quad_face(q);
                    if let Some(&cid) = quad_cell_map.get(&key) {
                        quad_face_nodes.push((cid, [pos[0], pos[1], 0.0]));
                    }
                }
            }
            quad_face_nodes.sort_by_key(|x| x.0);
            quad_face_nodes.dedup_by_key(|x| x.0);
        }

        Ok(MixedDiscretisation {
            dim: D,
            points,
            cell_types,
            cell_offsets,
            cell_conn,
            cell_labels,
            facet_types,
            facet_offsets,
            facet_conn,
            facet_labels,
            exterior_facets,
            edges,
            cell_edge_offsets,
            cell_edges,
            cell_edge_signs,
            tri_faces,
            cell_tri_offsets,
            cell_tri_faces,
            cell_tri_perms,
            quad_faces,
            cell_quad_offsets,
            cell_quad_faces,
            cell_quad_dihedrals,
            edge_nodes,
            quad_face_nodes,
            qp_edges,
        })
    }

    /// Builds a [`MixedDiscretisation`] from a [`SimplexMesh`].
    pub fn from_simplex_mesh(sm: &SimplexMesh) -> Result<Self, String> {
        crate::mesh::Mesh::from_solver_mesh(sm)?.discretise_mixed()
    }

    /// Number of vertices.
    #[inline]
    pub fn n_points(&self) -> usize {
        self.points.len() / self.dim
    }

    /// Number of cells.
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.cell_types.len()
    }

    /// Number of labelled facets.
    #[inline]
    pub fn n_facets(&self) -> usize {
        self.facet_types.len()
    }

    /// Number of unique edges.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Number of unique 3D triangular faces (`0` in 2D).
    #[inline]
    pub fn n_tri_faces(&self) -> usize {
        self.tri_faces.len()
    }

    /// Number of unique 3D quadrilateral faces (`0` in 2D).
    #[inline]
    pub fn n_quad_faces(&self) -> usize {
        self.quad_faces.len()
    }

    /// Coordinates of vertex `v` (`len = dim`).
    #[inline]
    pub fn point(&self, v: usize) -> &[f64] {
        &self.points[v * self.dim..(v + 1) * self.dim]
    }

    /// 3D-padded coordinates of vertex `v` (`[x, z, 0.0]` in 2D, `[x, y, z]` in 3D).
    #[inline]
    pub fn point3(&self, v: usize) -> [f64; 3] {
        let p = self.point(v);
        if self.dim == 2 {
            [p[0], p[1], 0.0]
        } else {
            [p[0], p[1], p[2]]
        }
    }

    /// Vertex IDs of cell `c`.
    #[inline]
    pub fn cell(&self, c: usize) -> &[u32] {
        &self.cell_conn[self.cell_offsets[c] as usize..self.cell_offsets[c + 1] as usize]
    }

    /// [`CellType`] of cell `c`.
    #[inline]
    pub fn cell_type(&self, c: usize) -> CellType {
        self.cell_types[c]
    }

    /// Vertex IDs of facet `f`.
    #[inline]
    pub fn facet(&self, f: usize) -> &[u32] {
        &self.facet_conn[self.facet_offsets[f] as usize..self.facet_offsets[f + 1] as usize]
    }

    /// [`FacetType`] of facet `f`.
    #[inline]
    pub fn facet_type(&self, f: usize) -> FacetType {
        self.facet_types[f]
    }

    /// Global edge IDs of cell `c` in [`cell_local_edges`] order.
    #[inline]
    pub fn cell_edge_ids(&self, c: usize) -> &[u32] {
        &self.cell_edges[self.cell_edge_offsets[c] as usize..self.cell_edge_offsets[c + 1] as usize]
    }

    /// Default orientation signs (`+1` / `-1`) of cell `c`'s local edges.
    #[inline]
    pub fn cell_edge_signs(&self, c: usize) -> &[i8] {
        &self.cell_edge_signs[self.cell_edge_offsets[c] as usize..self.cell_edge_offsets[c + 1] as usize]
    }

    /// Global tri-face IDs of cell `c` in [`cell_local_tri_faces`] order.
    #[inline]
    pub fn cell_tri_face_ids(&self, c: usize) -> &[u32] {
        &self.cell_tri_faces[self.cell_tri_offsets[c] as usize..self.cell_tri_offsets[c + 1] as usize]
    }

    /// `D_3` permutation indices (`0..6`) of cell `c`'s local tri faces.
    #[inline]
    pub fn cell_tri_perms(&self, c: usize) -> &[u8] {
        &self.cell_tri_perms[self.cell_tri_offsets[c] as usize..self.cell_tri_offsets[c + 1] as usize]
    }

    /// Global quad-face IDs of cell `c` in [`cell_local_quad_faces`] order.
    #[inline]
    pub fn cell_quad_face_ids(&self, c: usize) -> &[u32] {
        &self.cell_quad_faces[self.cell_quad_offsets[c] as usize..self.cell_quad_offsets[c + 1] as usize]
    }

    /// `D_4` dihedral indices (`0..8`) of cell `c`'s local quad faces.
    #[inline]
    pub fn cell_quad_dihedrals(&self, c: usize) -> &[u8] {
        &self.cell_quad_dihedrals[self.cell_quad_offsets[c] as usize..self.cell_quad_offsets[c + 1] as usize]
    }

    /// Global index of edge `(a, b)`, if present.
    #[inline]
    pub fn edge_index(&self, a: u32, b: u32) -> Option<usize> {
        self.edges.binary_search(&(a.min(b), a.max(b))).ok()
    }

    /// Global index of 3D triangular face `f`, if present.
    #[inline]
    pub fn tri_face_index(&self, f: [u32; 3]) -> Option<usize> {
        self.tri_faces.binary_search(&sorted_tri_face(f)).ok()
    }

    /// Global index of 3D quadrilateral face `cyc` (in cyclic order), if present.
    #[inline]
    pub fn quad_face_index(&self, cyc: [u32; 4]) -> Option<usize> {
        self.quad_faces.binary_search(&canonical_quad_face(cyc)).ok()
    }

    /// Curved mid-edge position of global edge `e`, if curved.
    #[inline]
    pub fn edge_node(&self, e: usize) -> Option<[f64; 3]> {
        self.edge_nodes
            .binary_search_by_key(&(e as u32), |x| x.0)
            .ok()
            .map(|i| self.edge_nodes[i].1)
    }

    /// Curved quad-face center node (by `quad_face_id` in 3D, or `cell_id` in 2D), if curved.
    #[inline]
    pub fn quad_face_node(&self, id: usize) -> Option<[f64; 3]> {
        self.quad_face_nodes
            .binary_search_by_key(&(id as u32), |x| x.0)
            .ok()
            .map(|i| self.quad_face_nodes[i].1)
    }

    /// Whether any edge or quad face in the mesh carries curved geometry nodes.
    #[inline]
    pub fn has_curved(&self) -> bool {
        !self.edge_nodes.is_empty() || !self.quad_face_nodes.is_empty()
    }

    /// Whether cell `c` has at least one curved edge or curved quad face.
    pub fn cell_is_curved(&self, c: usize) -> bool {
        if !self.has_curved() {
            return false;
        }
        if self.cell_edge_ids(c).iter().any(|&e| self.edge_node(e as usize).is_some()) {
            return true;
        }
        if !self.quad_face_nodes.is_empty() {
            if self.dim == 2 {
                if self.cell_types[c] == CellType::Quad && self.quad_face_node(c).is_some() {
                    return true;
                }
            } else if self.cell_quad_face_ids(c).iter().any(|&qf| self.quad_face_node(qf as usize).is_some()) {
                return true;
            }
        }
        false
    }

    /// Whether edge `e` carries a Barsoum quarter-point singularity node.
    #[inline]
    pub fn edge_is_qp(&self, e: usize) -> bool {
        self.qp_edges.binary_search(&(e as u32)).is_ok()
    }

    /// Whether the mesh carries any quarter-point edge.
    #[inline]
    pub fn has_qp(&self) -> bool {
        !self.qp_edges.is_empty()
    }

    /// Whether cell `c` touches at least one quarter-point edge.
    pub fn cell_is_qp(&self, c: usize) -> bool {
        !self.qp_edges.is_empty() && self.cell_edge_ids(c).iter().any(|&e| self.edge_is_qp(e as usize))
    }

    /// Whether every cell is a simplex (`Tri` in 2D, `Tet` in 3D) and every facet is a simplex (`Edge2` / `Tri3`).
    pub fn is_pure_simplex(&self) -> bool {
        let want_cell = if self.dim == 2 { CellType::Tri } else { CellType::Tet };
        let want_facet = if self.dim == 2 { FacetType::Edge2 } else { FacetType::Tri3 };
        self.cell_types.iter().all(|&t| t == want_cell) && self.facet_types.iter().all(|&f| f == want_facet)
    }

    /// Sorted unique domain labels present in `cell_labels`.
    pub fn domain_labels(&self) -> Vec<i32> {
        self.cell_labels.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
    }

    /// Sorted unique boundary labels present in `facet_labels`.
    pub fn boundary_labels(&self) -> Vec<i32> {
        self.facet_labels.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
    }

    /// Axis-aligned bounding box `(lo, hi)` (3rd coordinate `0.0` in 2D).
    pub fn bbox(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in self.points.chunks(self.dim) {
            for k in 0..self.dim {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        for k in self.dim..3 {
            lo[k] = 0.0;
            hi[k] = 0.0;
        }
        (lo, hi)
    }

    /// Verifies that the mixed mesh is 100% conforming with zero transition pyramids and zero
    /// mismatched `Tri`-to-`Quad` or hanging-diagonal contacts:
    /// - every cell has distinct in-range vertices and dimension matching `self.dim`,
    /// - every interior facet (`Edge` in 2D; `Tri` across `Tet`/`Prism` and `Quad` across `Prism`/`Hex` in 3D)
    ///   is shared by at most 2 cells with identical canonical face vertex sets (no figure-8 twisted quads),
    /// - zero `Tri`-to-`Quad` contacts: no 3-vertex subset of any `Quad` face is a `Tri` face, and no
    ///   diagonal of any 2D `Quad` cell or 3D `Quad` face appears in `self.edges`,
    /// - every labelled facet in `facet_conn` matches a valid topological facet of the mesh.
    pub fn validate_conforming(&self) -> Result<(), String> {
        if self.dim != 2 && self.dim != 3 {
            return Err(format!("unsupported dimension {}", self.dim));
        }
        let np = self.n_points();
        let nc = self.n_cells();
        if self.cell_offsets.len() != nc + 1 || self.cell_labels.len() != nc {
            return Err("cell CSR length mismatch".into());
        }

        for c in 0..nc {
            let t = self.cell_types[c];
            if t.dim() != self.dim {
                return Err(format!("cell {c} ({t:?}): dimension mismatch"));
            }
            let vs = self.cell(c);
            if vs.len() != t.nv() {
                return Err(format!("cell {c} ({t:?}): wrong vertex count {}", vs.len()));
            }
            for (i, &v) in vs.iter().enumerate() {
                if v as usize >= np {
                    return Err(format!("cell {c}: vertex {v} out of range"));
                }
                if vs[..i].contains(&v) {
                    return Err(format!("cell {c}: repeated vertex {v}"));
                }
            }
        }

        if self.dim == 2 {
            let mut edge_counts = vec![0u32; self.edges.len()];
            for c in 0..nc {
                for &eid in self.cell_edge_ids(c) {
                    edge_counts[eid as usize] += 1;
                }
                if self.cell_types[c] == CellType::Quad {
                    let vs = self.cell(c);
                    for (a, b) in [(vs[0], vs[2]), (vs[1], vs[3])] {
                        if self.edge_index(a, b).is_some() {
                            return Err(format!(
                                "cell {c} (Quad): diagonal ({a}, {b}) exists as an edge in the mesh (non-conforming Tri-to-Quad split)"
                            ));
                        }
                    }
                }
            }
            for (eid, &cnt) in edge_counts.iter().enumerate() {
                if cnt == 0 || cnt > 2 {
                    return Err(format!(
                        "2D edge {:?} has invalid cell incidence {cnt} (expected 1 or 2)",
                        self.edges[eid]
                    ));
                }
            }
            for f in 0..self.n_facets() {
                if self.facet_types[f] != FacetType::Edge2 {
                    return Err(format!("2D facet {f}: invalid facet type {:?}", self.facet_types[f]));
                }
                let fv = self.facet(f);
                if self.edge_index(fv[0], fv[1]).is_none() {
                    return Err(format!("2D facet {f} ({fv:?}) does not belong to any cell"));
                }
            }
        } else {
            let mut tri_counts = vec![0u32; self.tri_faces.len()];
            let mut quad_counts = vec![0u32; self.quad_faces.len()];
            for c in 0..nc {
                for &tf in self.cell_tri_face_ids(c) {
                    tri_counts[tf as usize] += 1;
                }
                for &qf in self.cell_quad_face_ids(c) {
                    quad_counts[qf as usize] += 1;
                }
            }
            for (tf, &cnt) in tri_counts.iter().enumerate() {
                if cnt == 0 || cnt > 2 {
                    return Err(format!(
                        "3D tri face {:?} has invalid cell incidence {cnt} (expected 1 or 2)",
                        self.tri_faces[tf]
                    ));
                }
            }
            for (qf, &cnt) in quad_counts.iter().enumerate() {
                if cnt == 0 || cnt > 2 {
                    return Err(format!(
                        "3D quad face {:?} has invalid cell incidence {cnt} (expected 1 or 2)",
                        self.quad_faces[qf]
                    ));
                }
            }

            // Detect twisted quad sharing (same 4 vertices in incompatible cyclic order) and Tri-to-Quad contact.
            let mut sorted_quad_sets: HashMap<[u32; 4], [u32; 4]> = HashMap::with_capacity(self.quad_faces.len());
            for &qf in &self.quad_faces {
                let mut set4 = qf;
                set4.sort_unstable();
                if let Some(&prev_cyc) = sorted_quad_sets.get(&set4) {
                    return Err(format!(
                        "twisted quad face sharing: vertex set {set4:?} appears with conflicting cyclic orders {prev_cyc:?} and {qf:?}"
                    ));
                }
                sorted_quad_sets.insert(set4, qf);

                // Check none of the 4 3-vertex subsets of `qf` is a triangular face in `tri_faces`.
                for sub in [
                    [qf[0], qf[1], qf[2]],
                    [qf[0], qf[1], qf[3]],
                    [qf[0], qf[2], qf[3]],
                    [qf[1], qf[2], qf[3]],
                ] {
                    if self.tri_face_index(sub).is_some() {
                        return Err(format!(
                            "non-conforming Tri-to-Quad contact: triangle {:?} lies on quad face {qf:?}",
                            sorted_tri_face(sub)
                        ));
                    }
                }
                // Check neither diagonal `(qf[0], qf[2])` nor `(qf[1], qf[3])` appears as a mesh edge.
                for (a, b) in [(qf[0], qf[2]), (qf[1], qf[3])] {
                    if self.edge_index(a, b).is_some() {
                        return Err(format!(
                            "non-conforming diagonal edge ({a}, {b}) cuts across quad face {qf:?}"
                        ));
                    }
                }
            }

            for f in 0..self.n_facets() {
                let fv = self.facet(f);
                match self.facet_types[f] {
                    FacetType::Edge2 => {
                        return Err(format!("3D facet {f} cannot be Edge2"));
                    }
                    FacetType::Tri3 => {
                        if self.tri_face_index([fv[0], fv[1], fv[2]]).is_none() {
                            return Err(format!("3D Tri3 facet {f} ({fv:?}) does not belong to any cell"));
                        }
                    }
                    FacetType::Quad4 => {
                        if self.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).is_none() {
                            return Err(format!("3D Quad4 facet {f} ({fv:?}) does not belong to any cell"));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mixed_2d_tri_and_quad_conforming_and_curved_edges() {
        // Left cell: Tri (0, 1, 2); Right cell: Quad (1, 3, 4, 2) sharing edge (1, 2).
        let mut cm = CellMesh::<2>::new();
        cm.points = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 0.0], [2.0, 1.0]];
        cm.push_cell(CellType::Tri, &[0, 1, 2], 1);
        cm.push_cell(CellType::Quad, &[1, 3, 4, 2], 2);
        cm.push_facet(&[0, 1], 1);
        cm.push_facet(&[1, 3], 1);
        cm.push_facet(&[3, 4], 2);
        cm.push_facet(&[4, 2], 3);
        cm.push_facet(&[2, 0], 4);
        // Interior interface facet (shared by 2 cells -> not exterior)
        cm.push_facet(&[1, 2], 5);
        cm.set_edge_nodes([([1, 2], [1.05, 0.5])]);
        cm.set_quad_face_nodes([([1, 3, 4, 2], [1.52, 0.5])]);

        assert!(!cm.is_pure_simplex());
        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        md.validate_conforming().unwrap();
        assert!(!md.is_pure_simplex());
        assert_eq!(md.n_cells(), 2);
        assert_eq!(md.n_edges(), 6);
        assert_eq!(md.exterior_facets, vec![true, true, true, true, true, false]);

        let e_shared = md.edge_index(1, 2).unwrap();
        assert_eq!(md.edge_node(e_shared), Some([1.05, 0.5, 0.0]));
        assert_eq!(md.quad_face_node(1), Some([1.52, 0.5, 0.0]));
        assert!(md.cell_is_curved(0));
        assert!(md.cell_is_curved(1));
    }

    #[test]
    fn mixed_3d_tet_prism_hex_d3_and_d4_orientations_and_conformity() {
        // Build a conforming Tet <-> Prism <-> Hex chain:
        // Hex on [0,1]x[0,1]x[0,1] (vertices 0..8),
        // Prism sharing x=1 quad face (1, 2, 6, 5) with triangle base (1, 8, 2) and top (5, 9, 6),
        // Tet sharing top tri face (5, 9, 6) with apex 10 at (1.3, 0.3, 2.0).
        let mut cm = CellMesh::<3>::new();
        cm.points = vec![
            [0.0, 0.0, 0.0], // 0
            [1.0, 0.0, 0.0], // 1
            [1.0, 1.0, 0.0], // 2
            [0.0, 1.0, 0.0], // 3
            [0.0, 0.0, 1.0], // 4
            [1.0, 0.0, 1.0], // 5
            [1.0, 1.0, 1.0], // 6
            [0.0, 1.0, 1.0], // 7
            [2.0, 0.0, 0.0], // 8
            [2.0, 0.0, 1.0], // 9
            [1.3, 0.3, 2.0], // 10
        ];
        cm.push_cell(CellType::Hex, &[0, 1, 2, 3, 4, 5, 6, 7], 1);
        cm.push_cell(CellType::Prism, &[1, 8, 2, 5, 9, 6], 2);
        cm.push_cell(CellType::Tet, &[5, 9, 6, 10], 3);
        cm.push_facet(&[0, 3, 2, 1], 1);
        cm.push_facet(&[5, 9, 10], 2);
        // Shared interior quad face as sheet facet -> exterior_facets == false
        cm.push_facet(&[1, 2, 6, 5], 9);

        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        md.validate_conforming().unwrap();
        assert_eq!(md.n_cells(), 3);
        assert_eq!(md.exterior_facets, vec![true, true, false]);

        // Verify Hex and Prism agree on the shared Quad face ID and canonical cyclic vertices!
        let hex_qf_pos = HEX_QUAD_FACES.iter().position(|f| *f == [1, 2, 6, 5]).unwrap();
        let prism_qf_pos = PRISM_QUAD_FACES.iter().position(|f| *f == [2, 0, 3, 5]).unwrap();
        let qf_hex = md.cell_quad_face_ids(0)[hex_qf_pos];
        let qf_prism = md.cell_quad_face_ids(1)[prism_qf_pos];
        assert_eq!(qf_hex, qf_prism);

        let hex_vs = md.cell(0);
        let prism_vs = md.cell(1);
        let hex_local_q = HEX_QUAD_FACES[hex_qf_pos].map(|i| hex_vs[i]);
        let prism_local_q = PRISM_QUAD_FACES[prism_qf_pos].map(|i| prism_vs[i]);
        let canon_from_hex = d4_apply(md.cell_quad_dihedrals(0)[hex_qf_pos], hex_local_q);
        let canon_from_prism = d4_apply(md.cell_quad_dihedrals(1)[prism_qf_pos], prism_local_q);
        assert_eq!(canon_from_hex, canon_from_prism);
        assert_eq!(canon_from_hex, md.quad_faces[qf_hex as usize]);

        // Verify Prism and Tet agree on the shared Tri face ID and canonical D3 vertices!
        let prism_tf = md.cell_tri_face_ids(1)[1]; // top tri [3, 4, 5] -> [5, 9, 6]
        let tet_tf_pos = TET_TRI_FACES.iter().position(|f| *f == [0, 1, 2]).unwrap();
        let tet_tf = md.cell_tri_face_ids(2)[tet_tf_pos];
        assert_eq!(prism_tf, tet_tf);
        let prism_local_t = PRISM_TRI_FACES[1].map(|i| prism_vs[i]);
        let tet_vs = md.cell(2);
        let tet_local_t = TET_TRI_FACES[tet_tf_pos].map(|i| tet_vs[i]);
        assert_eq!(
            d3_apply(md.cell_tri_perms(1)[1], prism_local_t),
            d3_apply(md.cell_tri_perms(2)[tet_tf_pos], tet_local_t)
        );
    }

    #[test]
    fn validate_conforming_rejects_tri_to_quad_contact() {
        // Non-conforming mesh: a Hex sharing half of its x=1 Quad face (1, 2, 6, 5) with a Tet!
        let mut cm = CellMesh::<3>::new();
        cm.points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [2.0, 0.5, 0.5],
        ];
        cm.push_cell(CellType::Hex, &[0, 1, 2, 3, 4, 5, 6, 7], 1);
        cm.push_cell(CellType::Tet, &[1, 2, 6, 8], 2);
        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        let err = md.validate_conforming().unwrap_err();
        assert!(err.contains("Tri-to-Quad"), "unexpected error: {err}");
    }
}
