//! Structured PML generator around a box-shaped physical mesh.
//!
//! Given a conforming simplex mesh whose boundary lies on the faces of an axis-aligned box, this
//! builds the absorbing layers as *structured* cells that match the physical boundary mesh node
//! for node:
//!
//! * 2D: quadrilaterals. A face region (one stretched direction) extrudes every boundary edge of a
//!   side outward; a corner region (two directions) is a tensor grid of quads.
//! * 3D: prisms on faces (every boundary triangle extruded: triangular faces meet the physical
//!   tetrahedra, quad faces meet each other), hexahedra on box edges (the chain of boundary
//!   vertices along the edge times a 2D layer grid) and hexahedra at corners.
//!
//! No pyramids or other transition elements are needed because the physical side always meets a
//! triangle (or an edge in 2D).
//!
//! Each region is its own axis-aligned box with its own domain label and set of stretched
//! directions (`code`: -1 or +1 on every stretched axis, 0 elsewhere), because a PML solver takes the
//! bounding box of a domain's cells to define the stretching coordinate. A material that crosses
//! the interface continues into a PML copy of itself: the PML label encodes `(region, material)`,
//! see [`PmlSpec::label`]. The material of a face-region cell is that of the physical cell owning
//! the boundary facet it extrudes; for edge and corner regions it is the material of the owner of
//! the lowest-index boundary facet on the first stretched side touching that edge interval or
//! corner.
//!
//! The result is a [`CellMesh`]; [`CellMesh::split_to_simplices`] gives the simplex fallback for
//! solvers that cannot take quads, prisms or hexahedra yet. Physical vertices keep their ids.

use crate::mesh::cells::{CellMesh, CellType};
use crate::adapt::work::{FacetKey, SimplexMesh, bbox_of, facet_key};
use std::collections::HashMap;

/// Layer specification. Axis `a`, side `s` (0 = minimum, 1 = maximum face of the box); a side
/// with zero thickness or zero layers gets no PML.
#[derive(Clone, Debug)]
pub struct PmlSpec<const D: usize> {
    /// Thickness of the layer on each side.
    pub thickness: [[f64; 2]; D],
    /// Number of cell layers across the thickness on each side.
    pub layers: [[usize; 2]; D],
    /// Ratio of the last (outermost) to the first layer width; 1 for uniform layers.
    pub grading: f64,
    /// Label of the PML cell of region index `r` (see [`PmlSpec::region_index`]) and material
    /// `m` is `label_base + r * label_stride + m`. Materials must lie in `0..label_stride`.
    pub label_base: i32,
    pub label_stride: i32,
    /// Label for the boundary facets on the outer surface of the layers (`None`: leave them
    /// unlabelled).
    pub outer_label: Option<i32>,
    /// Relative tolerance (to the box size) for a vertex to count as lying on a box face.
    pub tol: f64,
}

impl<const D: usize> PmlSpec<D> {
    /// The same `thickness` and number of `layers` on every side, uniform spacing.
    pub fn uniform(thickness: f64, layers: usize) -> Self {
        PmlSpec {
            thickness: [[thickness; 2]; D],
            layers: [[layers; 2]; D],
            grading: 1.0,
            label_base: 1000,
            label_stride: 1000,
            outer_label: None,
            tol: 1e-9,
        }
    }

    /// Index of a region code `c` in `0..3^D`: `sum (c_a + 1) 3^a`. The all-zero code (the
    /// physical domain) has index `(3^D - 1) / 2`.
    pub fn region_index(code: &[i8; D]) -> usize {
        let mut r = 0;
        let mut p = 1;
        for a in 0..D {
            r += (code[a] + 1) as usize * p;
            p *= 3;
        }
        r
    }

    /// Inverse of [`PmlSpec::region_index`].
    pub fn region_code(mut index: usize) -> [i8; D] {
        let mut c = [0i8; D];
        for a in 0..D {
            c[a] = (index % 3) as i8 - 1;
            index /= 3;
        }
        c
    }

    /// PML label for region `code` and physical material label `material`.
    pub fn label(&self, code: &[i8; D], material: i32) -> i32 {
        debug_assert!((0..self.label_stride).contains(&material));
        self.label_base + Self::region_index(code) as i32 * self.label_stride + material
    }

    /// Inverse of [`PmlSpec::label`]: `(region code, material)`, or `None` for a physical label.
    pub fn decode(&self, label: i32) -> Option<([i8; D], i32)> {
        let rel = label.checked_sub(self.label_base)?;
        if rel < 0 {
            return None;
        }
        let (r, m) = (rel / self.label_stride, rel % self.label_stride);
        (r < 3i32.pow(D as u32)).then(|| (Self::region_code(r as usize), m))
    }

    /// Positions of the layer interfaces on side `(a, s)` as distances from the box face,
    /// `layers + 1` values from 0 to the thickness. Widths grow geometrically by
    /// `grading^(1 / (layers - 1))` per layer.
    pub fn layer_offsets(&self, a: usize, s: usize) -> Vec<f64> {
        let (n, t) = (self.layers[a][s], self.thickness[a][s]);
        let mut off = vec![0.0; n + 1];
        if n == 0 {
            return off;
        }
        let q = if n > 1 { self.grading.powf(1.0 / (n as f64 - 1.0)) } else { 1.0 };
        let total: f64 = (0..n).map(|k| q.powi(k as i32)).sum();
        let w0 = t / total;
        for k in 0..n {
            off[k + 1] = off[k] + w0 * q.powi(k as i32);
        }
        off[n] = t;
        off
    }

    fn active(&self, a: usize, s: usize) -> bool {
        self.layers[a][s] > 0 && self.thickness[a][s] > 0.0
    }
}

/// Layers needed to resolve `thickness` with cells of normal size about `h` (at least 1).
pub fn layers_for_size(thickness: f64, h: f64) -> usize {
    ((thickness / h).round() as usize).max(1)
}

/// One PML region.
#[derive(Clone, Debug)]
pub struct PmlRegion<const D: usize> {
    /// Stretched directions: `-1` / `+1` for the minimum / maximum side along an axis, 0 for
    /// axes along which the region is not stretched.
    pub code: [i8; D],
    /// Bounding box of the region's cells.
    pub bbox: [[f64; 2]; D],
    /// `(physical material, PML label)` pairs present in this region.
    pub labels: Vec<(i32, i32)>,
}

impl<const D: usize> PmlRegion<D> {
    /// Number of stretched directions (1 on faces, 2 on edges, 3 at corners).
    pub fn n_directions(&self) -> usize {
        self.code.iter().filter(|&&c| c != 0).count()
    }
}

/// Output of [`build_pml`].
#[derive(Clone, Debug)]
pub struct PmlMesh<const D: usize> {
    /// Physical simplices (ids and order unchanged) followed by the PML cells.
    pub mesh: CellMesh<D>,
    pub regions: Vec<PmlRegion<D>>,
    pub n_physical_cells: usize,
    pub n_physical_points: usize,
    /// The box the physical mesh fills.
    pub interface: [[f64; 2]; D],
    pub spec: PmlSpec<D>,
}

fn key_bits<const D: usize>(p: &[f64; D]) -> [u64; D] {
    // +0.0 and -0.0 must hash alike
    std::array::from_fn(|a| (p[a] + 0.0).to_bits())
}

struct Builder<'a, const D: usize> {
    phys: &'a SimplexMesh<D>,
    out: CellMesh<D>,
    lookup: HashMap<[u64; D], u32>,
    /// layer offsets measured outward, indexed `[axis][side]`
    offsets: Vec<[Vec<f64>; 2]>,
    bbox: [[f64; 2]; D],
}

impl<'a, const D: usize> Builder<'a, D> {
    /// Id of the vertex at `p`, created if new.
    fn vertex(&mut self, p: [f64; D]) -> u32 {
        let k = key_bits(&p);
        if let Some(&id) = self.lookup.get(&k) {
            return id;
        }
        let id = self.out.points.len() as u32;
        self.out.points.push(p);
        self.lookup.insert(k, id);
        id
    }

    /// Absolute coordinate of layer interface `k` on side `(a, s)`.
    fn layer_coord(&self, a: usize, s: usize, k: usize) -> f64 {
        let sign = if s == 0 { -1.0 } else { 1.0 };
        self.bbox[a][s] + sign * self.offsets[a][s][k]
    }

    /// Vertex obtained from the physical vertex `base` by moving the stretched axes (`steps` lists
    /// `(axis, side, layer index)`) to their layer coordinates; layer 0 keeps the vertex's own
    /// coordinate so that physical vertices are reused exactly.
    fn shifted(&mut self, base: u32, steps: &[(usize, usize, usize)]) -> u32 {
        let mut p = self.phys.points[base as usize];
        for &(a, s, k) in steps {
            if k > 0 {
                p[a] = self.layer_coord(a, s, k);
            }
        }
        self.vertex(p)
    }

    fn push(&mut self, t: CellType, vs: &[u32], label: i32) {
        let mut v = vs.to_vec();
        orient(t, &mut v, &self.out.points);
        self.out.push_cell(t, &v, label);
    }
}

/// Reorder the nodes of a cell in place so that its scaled Jacobian at node 0 is positive.
fn orient<const D: usize>(t: CellType, v: &mut [u32], pts: &[[f64; D]]) {
    let mut probe = CellMesh::<D>::new();
    probe.points = pts.to_vec();
    probe.push_cell(t, v, 0);
    if crate::mesh::quality::cell_scaled_jacobian_min(&probe, 0) < 0.0 {
        let swaps: &[(usize, usize)] = match t {
            CellType::Tri => &[(1, 2)],
            CellType::Quad => &[(1, 3)],
            CellType::Tet => &[(2, 3)],
            CellType::Prism => &[(1, 2), (4, 5)],
            CellType::Hex => &[(1, 3), (5, 7)],
        };
        for &(i, j) in swaps {
            v.swap(i, j);
        }
    }
}

/// Build the structured PML around `phys` (see the module docs).
///
/// Errors if the physical boundary is not made of facets lying on the faces of its bounding box,
/// or if an edge or corner of the box needed for an active layer is not a vertex chain of the
/// physical mesh.
pub fn build_pml<const D: usize>(phys: &SimplexMesh<D>, spec: &PmlSpec<D>) -> Result<PmlMesh<D>, String> {
    if !(D == 2 || D == 3) {
        return Err("only D = 2 and D = 3 are supported".into());
    }
    if spec.grading <= 0.0 || !spec.grading.is_finite() {
        return Err("grading must be positive".into());
    }
    if spec.label_stride <= 0 {
        return Err("label_stride must be positive".into());
    }
    let bbox = phys.bbox();
    let size = (0..D).map(|a| bbox[a][1] - bbox[a][0]).fold(0.0, f64::max);
    let tol = spec.tol * size;
    let on_plane = |p: &[f64; D], a: usize, s: usize| (p[a] - bbox[a][s]).abs() <= tol;

    // boundary facets of the physical mesh and the plane each one lies on
    let fc = phys.facet_cells();
    let mut bfacets: Vec<(Vec<u32>, usize, usize, i32)> = Vec::new(); // (verts, axis, side, owner label)
    let mut keys: Vec<(&FacetKey, &Vec<u32>)> = fc.iter().filter(|(_, cs)| cs.len() == 1).collect();
    keys.sort_unstable_by_key(|(k, _)| **k);
    for (k, cs) in keys {
        let vs: Vec<u32> = k.iter().take(D).cloned().collect();
        let mut plane = None;
        'outer: for a in 0..D {
            for s in 0..2 {
                if vs.iter().all(|&v| on_plane(&phys.points[v as usize], a, s)) {
                    plane = Some((a, s));
                    break 'outer;
                }
            }
        }
        let (a, s) = plane.ok_or_else(|| {
            format!("physical boundary facet {vs:?} does not lie on a face of the bounding box; the domain must be box-shaped")
        })?;
        bfacets.push((vs, a, s, phys.cell_label[cs[0] as usize]));
    }
    // restore the original facet orientation (the key is sorted) from the owning cell
    for bf in bfacets.iter_mut() {
        let cell = fc[&facet_key(&bf.0)][0] as usize;
        let cf = phys.cell_facets(cell);
        for f in cf.iter().take(D + 1) {
            if facet_key(&f[..D]) == facet_key(&bf.0) {
                bf.0 = f[..D].to_vec();
            }
        }
    }

    let offsets: Vec<[Vec<f64>; 2]> = (0..D).map(|a| [spec.layer_offsets(a, 0), spec.layer_offsets(a, 1)]).collect();
    let mut b = Builder { phys, out: CellMesh::new(), lookup: HashMap::new(), offsets, bbox };
    b.out = CellMesh::from_simplex(phys);
    b.out.facet_offsets = vec![0];
    b.out.facet_conn.clear();
    b.out.facet_label.clear();
    // physical vertices that sit on a box face may be reached by coordinates
    for (i, p) in phys.points.iter().enumerate() {
        if (0..D).any(|a| (0..2).any(|s| on_plane(p, a, s))) {
            b.lookup.insert(key_bits(p), i as u32);
        }
    }

    let mut regions: HashMap<usize, PmlRegion<D>> = HashMap::new();
    let mut note = |code: [i8; D], bbox_pts: &[[f64; D]], material: i32, label: i32| {
        let r = regions.entry(PmlSpec::<D>::region_index(&code)).or_insert_with(|| PmlRegion {
            code,
            bbox: [[f64::INFINITY, f64::NEG_INFINITY]; D],
            labels: Vec::new(),
        });
        let bb = bbox_of(bbox_pts);
        for a in 0..D {
            r.bbox[a][0] = r.bbox[a][0].min(bb[a][0]);
            r.bbox[a][1] = r.bbox[a][1].max(bb[a][1]);
        }
        if !r.labels.contains(&(material, label)) {
            r.labels.push((material, label));
        }
    };

    // enumerate region codes in index order for determinism
    let n_codes = 3usize.pow(D as u32);
    for ri in 0..n_codes {
        let code = PmlSpec::<D>::region_code(ri);
        let axes: Vec<usize> = (0..D).filter(|&a| code[a] != 0).collect();
        if axes.is_empty() {
            continue;
        }
        let side = |a: usize| if code[a] < 0 { 0 } else { 1 };
        if axes.iter().any(|&a| !spec.active(a, side(a))) {
            continue;
        }
        let first = axes[0];

        if axes.len() == 1 {
            // face region: extrude every boundary facet on plane (first, side)
            let (a, s) = (first, side(first));
            let n = spec.layers[a][s];
            for (vs, fa, fs, mat) in bfacets.iter() {
                if *fa != a || *fs != s {
                    continue;
                }
                if !(0..spec.label_stride).contains(mat) {
                    return Err(format!("material label {mat} outside 0..label_stride"));
                }
                let label = spec.label(&code, *mat);
                for k in 0..n {
                    let lo: Vec<u32> = vs.iter().map(|&v| b.shifted(v, &[(a, s, k)])).collect();
                    let hi: Vec<u32> = vs.iter().map(|&v| b.shifted(v, &[(a, s, k + 1)])).collect();
                    let mut nodes = lo.clone();
                    let t = if D == 2 {
                        nodes.extend([hi[1], hi[0]]);
                        CellType::Quad
                    } else {
                        nodes.extend(hi.iter());
                        CellType::Prism
                    };
                    b.push(t, &nodes, label);
                    let pts: Vec<[f64; D]> = nodes.iter().map(|&v| b.out.points[v as usize]).collect();
                    note(code, &pts, *mat, label);
                }
            }
            continue;
        }

        // edge / corner region: structured block over the stretched axes times (in 3D edges) the
        // chain of physical vertices along the free axis
        let free: Vec<usize> = (0..D).filter(|&a| code[a] == 0).collect();
        let mut chain: Vec<u32> = (0..phys.n_points() as u32)
            .filter(|&v| axes.iter().all(|&a| on_plane(&phys.points[v as usize], a, side(a))))
            .collect();
        if free.is_empty() {
            if chain.len() != 1 {
                return Err(format!(
                    "box corner {code:?} is {} vertices of the physical mesh, expected exactly one",
                    chain.len()
                ));
            }
        } else {
            let t = free[0];
            chain.sort_by(|&p, &q| {
                phys.points[p as usize][t].partial_cmp(&phys.points[q as usize][t]).unwrap().then(p.cmp(&q))
            });
            if chain.len() < 2 {
                return Err(format!("box edge {code:?} has fewer than two vertices"));
            }
        }
        // material of the lowest-index boundary facet on the first stretched side touching the
        // given vertices
        let material_of = |vs: &[u32]| -> i32 {
            bfacets
                .iter()
                .find(|(fv, fa, fs, _)| *fa == first && *fs == side(first) && vs.iter().all(|v| fv.contains(v)))
                .map(|x| x.3)
                .unwrap_or(0)
        };
        // layer index ranges over the stretched axes
        let dims: Vec<usize> = axes.iter().map(|&a| spec.layers[a][side(a)]).collect();
        let n_blocks: usize = dims.iter().product();
        let bases: Vec<Vec<u32>> =
            if free.is_empty() { vec![vec![chain[0]]] } else { chain.windows(2).map(|w| vec![w[0], w[1]]).collect() };
        for base in bases {
            let mat = material_of(&base);
            if !(0..spec.label_stride).contains(&mat) {
                return Err(format!("material label {mat} outside 0..label_stride"));
            }
            let label = spec.label(&code, mat);
            for blk in 0..n_blocks {
                // multi-index of the block's lower layer corner
                let mut idx = vec![0usize; axes.len()];
                let mut r = blk;
                for (i, d) in dims.iter().enumerate() {
                    idx[i] = r % d;
                    r /= d;
                }
                // local binary axes: [tangential (if any), stretched axes...]
                let n_local = axes.len() + free.len();
                let mut nodes = Vec::with_capacity(1 << n_local);
                let order: Vec<Vec<usize>> = if n_local == 2 {
                    vec![vec![0, 0], vec![1, 0], vec![1, 1], vec![0, 1]]
                } else {
                    vec![
                        vec![0, 0, 0],
                        vec![1, 0, 0],
                        vec![1, 1, 0],
                        vec![0, 1, 0],
                        vec![0, 0, 1],
                        vec![1, 0, 1],
                        vec![1, 1, 1],
                        vec![0, 1, 1],
                    ]
                };
                for bits in order {
                    let (bv, off) = if free.is_empty() { (base[0], 0) } else { (base[bits[0]], 1) };
                    let steps: Vec<(usize, usize, usize)> =
                        axes.iter().enumerate().map(|(i, &a)| (a, side(a), idx[i] + bits[off + i])).collect();
                    nodes.push(b.shifted(bv, &steps));
                }
                let pts: Vec<[f64; D]> = nodes.iter().map(|&v| b.out.points[v as usize]).collect();
                let t = if D == 2 { CellType::Quad } else { CellType::Hex };
                b.push(t, &nodes, label);
                note(code, &pts, mat, label);
            }
        }
    }

    // boundary facets: surviving physical ones, plus the outer surface of the layers
    let mut out = b.out;
    let n_phys_cells = phys.n_cells();
    let faces = out.boundary_faces_ordered();
    let mut phys_facets: HashMap<FacetKey, i32> = HashMap::new();
    for f in 0..phys.n_facets() {
        phys_facets.insert(facet_key(phys.facet(f)), phys.facet_label[f]);
    }
    for f in faces {
        let key = facet_key(&f);
        match phys_facets.get(&key) {
            Some(&l) => out.push_facet(&f, l),
            None => {
                if let Some(l) = spec.outer_label {
                    out.push_facet(&f, l);
                }
            }
        }
    }

    let mut regs: Vec<PmlRegion<D>> = regions.into_values().collect();
    regs.sort_by_key(|r| PmlSpec::<D>::region_index(&r.code));
    for r in regs.iter_mut() {
        r.labels.sort();
    }
    Ok(PmlMesh {
        mesh: out,
        regions: regs,
        n_physical_cells: n_phys_cells,
        n_physical_points: phys.n_points(),
        interface: bbox,
        spec: spec.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::fixtures::{box_tet, box_tri, rect_tri};
    use crate::mesh::quality::cellmesh_quality;

    fn area_with_pml(w: f64, h: f64, spec: &PmlSpec<2>) -> f64 {
        let ext = |a: usize, len: f64| {
            len + spec.thickness[a][0] * (spec.layers[a][0] > 0) as u8 as f64
                + spec.thickness[a][1] * (spec.layers[a][1] > 0) as u8 as f64
        };
        ext(0, w) * ext(1, h)
    }

    #[test]
    fn labels_roundtrip() {
        let spec = PmlSpec::<3>::uniform(1.0, 3);
        for r in 0..27 {
            let c = PmlSpec::<3>::region_code(r);
            assert_eq!(PmlSpec::<3>::region_index(&c), r);
            for m in [0, 7, 999] {
                if r != 13 {
                    assert_eq!(spec.decode(spec.label(&c, m)), Some((c, m)));
                }
            }
        }
        assert_eq!(spec.decode(5), None);
    }

    #[test]
    fn layer_offsets_uniform_and_graded() {
        let mut spec = PmlSpec::<2>::uniform(2.0, 4);
        let u = spec.layer_offsets(0, 0);
        assert_eq!(u.len(), 5);
        for k in 0..5 {
            assert!((u[k] - 0.5 * k as f64).abs() < 1e-15);
        }
        spec.grading = 4.0;
        let g = spec.layer_offsets(0, 1);
        let w: Vec<f64> = g.windows(2).map(|x| x[1] - x[0]).collect();
        assert!((w[3] / w[0] - 4.0).abs() < 1e-12);
        assert_eq!(g[4], 2.0);
        assert!(w.windows(2).all(|x| x[1] > x[0]));
    }

    #[test]
    fn two_d_full_ring_is_conforming_and_exact() {
        let phys = box_tri([[0.0, 2.0], [0.0, 1.0]], 8, 5);
        let mut spec = PmlSpec::<2>::uniform(0.5, 3);
        spec.outer_label = Some(99);
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        assert_eq!(p.regions.len(), 8);
        let c = p.mesh.counts();
        assert_eq!(c[&CellType::Tri], phys.n_cells());
        // faces: (8 + 8 + 5 + 5) edges x 3 layers, corners 4 x 9
        assert_eq!(c[&CellType::Quad], (8 + 8 + 5 + 5) * 3 + 4 * 9);
        assert!((p.mesh.total_measure() - 3.0 * 2.0).abs() < 1e-12);
        assert_eq!(p.mesh.n_points(), phys.n_points() + count_new_vertices_2d(8, 5, 3));
        // outer boundary: every boundary face of the whole mesh is on the 4 outer sides
        let bk = p.mesh.boundary_keys();
        assert_eq!(bk.len(), 2 * (8 + 6) + 2 * (5 + 6));
        assert_eq!(p.mesh.n_facets(), bk.len());
        assert!(p.mesh.facet_label.iter().all(|&l| l == 99));
        // all PML cells are rectangles: scaled Jacobian 1
        let q = cellmesh_quality(&p.mesh);
        assert!(q.min > 0.4);
        for ci in p.n_physical_cells..p.mesh.n_cells() {
            let v = crate::mesh::quality::cell_scaled_jacobian_min(&p.mesh, ci);
            assert!((v - 1.0).abs() < 1e-12, "cell {ci}: {v}");
        }
        // area by region
        assert!((p.mesh.total_measure() - area_with_pml(2.0, 1.0, &spec)).abs() < 1e-12);
    }

    /// Vertices created by the generator for an `nx x ny` interface boundary with `n` layers
    /// on all four sides (interface vertex count excluded).
    fn count_new_vertices_2d(nx: usize, ny: usize, n: usize) -> usize {
        // each of the four strips adds n layers of its tangential vertices (the corner vertices
        // of the interface are tangential vertices of both adjacent strips), each corner block
        // adds the n x n vertices strictly outside both strips
        2 * n * (nx + 1) + 2 * n * (ny + 1) + 4 * n * n
    }

    #[test]
    fn region_labels_directions_and_boxes() {
        let phys = box_tri([[0.0, 1.0], [0.0, 1.0]], 4, 4);
        let spec = PmlSpec::<2>::uniform(0.25, 2);
        let p = build_pml(&phys, &spec).unwrap();
        for r in &p.regions {
            let n = r.n_directions();
            assert!(n == 1 || n == 2);
            // the region box has the right thickness along each stretched axis and touches the
            // interface box
            for a in 0..2 {
                let w = r.bbox[a][1] - r.bbox[a][0];
                if r.code[a] == 0 {
                    assert!((w - 1.0).abs() < 1e-12, "{:?} {:?}", r.code, r.bbox);
                } else {
                    assert!((w - 0.25).abs() < 1e-12);
                    let edge = if r.code[a] < 0 { r.bbox[a][1] } else { r.bbox[a][0] };
                    assert!((edge - if r.code[a] < 0 { 0.0 } else { 1.0 }).abs() < 1e-12);
                }
            }
            for (m, l) in &r.labels {
                assert_eq!(spec.decode(*l), Some((r.code, *m)));
            }
        }
        // every PML cell carries a label that decodes to a region whose box contains it
        for c in p.n_physical_cells..p.mesh.n_cells() {
            let (code, _) = spec.decode(p.mesh.label[c]).unwrap();
            let r = p.regions.iter().find(|r| r.code == code).unwrap();
            for &v in p.mesh.cell(c) {
                for a in 0..2 {
                    let x = p.mesh.points[v as usize][a];
                    assert!(x >= r.bbox[a][0] - 1e-12 && x <= r.bbox[a][1] + 1e-12);
                }
            }
        }
    }

    #[test]
    fn partial_pml_and_asymmetric_thickness() {
        let phys = box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5);
        let mut spec = PmlSpec::<2>::uniform(0.0, 0);
        spec.thickness[0][1] = 0.4;
        spec.layers[0][1] = 4;
        spec.thickness[1][0] = 0.3;
        spec.layers[1][0] = 2;
        spec.grading = 3.0;
        spec.outer_label = Some(7);
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        // right strip, bottom strip, and the bottom-right corner (both sides active)
        let codes: Vec<[i8; 2]> = p.regions.iter().map(|r| r.code).collect();
        assert_eq!(codes, vec![[0, -1], [1, -1], [1, 0]]);
        assert!((p.mesh.total_measure() - 1.4 * 1.3).abs() < 1e-12);
        // physical sides without PML keep their boundary and labels
        let phys_labels: std::collections::HashSet<i32> = p.mesh.facet_label.iter().cloned().collect();
        assert!(phys_labels.contains(&crate::mesh::fixtures::side_label(0, false)));
        assert!(phys_labels.contains(&crate::mesh::fixtures::side_label(1, true)));
        assert!(phys_labels.contains(&7));
    }

    #[test]
    fn materials_cross_the_interface_into_pml_copies() {
        // two materials: left half label 1, right half label 2; a structure crossing the top PML
        let mut phys = box_tri([[0.0, 2.0], [0.0, 1.0]], 8, 4);
        for c in 0..phys.n_cells() {
            let cx: f64 = phys.cell(c).iter().map(|&v| phys.points[v as usize][0]).sum::<f64>() / 3.0;
            phys.cell_label[c] = if cx < 1.0 { 1 } else { 2 };
        }
        let spec = PmlSpec::<2>::uniform(0.5, 2);
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        let top = p.regions.iter().find(|r| r.code == [0, 1]).unwrap();
        assert_eq!(top.labels.len(), 2);
        assert_eq!(top.labels[0].0, 1);
        assert_eq!(top.labels[1].0, 2);
        // left strip only sees material 1, right strip only 2
        assert_eq!(p.regions.iter().find(|r| r.code == [-1, 0]).unwrap().labels.len(), 1);
        assert_eq!(p.regions.iter().find(|r| r.code == [-1, 0]).unwrap().labels[0].0, 1);
        assert_eq!(p.regions.iter().find(|r| r.code == [1, 0]).unwrap().labels[0].0, 2);
        // each cell above x < 1 carries the PML copy of material 1
        for c in p.n_physical_cells..p.mesh.n_cells() {
            let (code, mat) = spec.decode(p.mesh.label[c]).unwrap();
            if code == [0, 1] {
                let cx: f64 = p.mesh.cell(c).iter().map(|&v| p.mesh.points[v as usize][0]).sum::<f64>() / 4.0;
                assert_eq!(mat, if cx < 1.0 { 1 } else { 2 });
            }
        }
    }

    #[test]
    fn graded_interface_mesh_gives_conforming_graded_pml() {
        let xs: Vec<f64> = (0..=10).map(|i| (i as f64 / 10.0).powi(2)).collect();
        let ys: Vec<f64> = (0..=6).map(|i| i as f64 * 0.2).collect();
        let phys = rect_tri(&xs, &ys);
        let mut spec = PmlSpec::<2>::uniform(0.3, 4);
        spec.grading = 2.0;
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        let split = p.mesh.split_to_simplices().unwrap();
        split.check().unwrap();
        assert!((split.total_measure() - 1.6 * 1.8).abs() < 1e-12);
        assert_eq!(p.mesh.boundary_keys().len(), 2 * (10 + 8) + 2 * (6 + 8));
    }

    #[test]
    fn rejects_non_box_domains() {
        // L-shaped domain: boundary facets on the inner corner are not on the bounding box
        let mut m = box_tri([[0.0, 2.0], [0.0, 2.0]], 4, 4);
        let keep: Vec<usize> = (0..m.n_cells())
            .filter(|&c| {
                let cx: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][0]).sum::<f64>() / 3.0;
                let cy: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][1]).sum::<f64>() / 3.0;
                !(cx > 1.0 && cy > 1.0)
            })
            .collect();
        let mut l = SimplexMesh::<2> { points: m.points.clone(), ..Default::default() };
        for c in keep {
            l.push_cell(m.cell(c), 1);
        }
        m.facets.clear();
        let err = build_pml(&l, &PmlSpec::<2>::uniform(0.5, 2)).unwrap_err();
        assert!(err.contains("box-shaped"), "{err}");
    }

    #[test]
    fn deterministic() {
        let phys = box_tri([[0.0, 1.0], [0.0, 1.0]], 6, 6);
        let spec = PmlSpec::<2>::uniform(0.3, 3);
        let a = build_pml(&phys, &spec).unwrap();
        let b = build_pml(&phys, &spec).unwrap();
        assert_eq!(a.mesh.conn, b.mesh.conn);
        assert_eq!(a.mesh.label, b.mesh.label);
        for (p, q) in a.mesh.points.iter().zip(&b.mesh.points) {
            assert_eq!(p[0].to_bits(), q[0].to_bits());
            assert_eq!(p[1].to_bits(), q[1].to_bits());
        }
    }

    // ---- 3D ----

    #[test]
    fn three_d_full_shell_is_conforming_with_prisms_and_hexes() {
        let phys = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [3, 3, 3]);
        let mut spec = PmlSpec::<3>::uniform(0.4, 2);
        spec.outer_label = Some(50);
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        let c = p.mesh.counts();
        let faces = 6 * (2 * 3 * 3); // triangles on each box side: 2 per quad, 3x3 quads
        assert_eq!(c[&CellType::Prism], faces * 2);
        assert_eq!(c[&CellType::Hex], 12 * 3 * 2 * 2 + 8 * 2 * 2 * 2);
        assert!((p.mesh.total_measure() - 1.8f64.powi(3)).abs() < 1e-11);
        assert_eq!(p.regions.len(), 26);
        for r in &p.regions {
            assert!((1..=3).contains(&r.n_directions()));
        }
        // outer surface: over each box side the prisms end in 2*3*3 triangles; each of the 12 edge
        // regions ends in 2 * (3 intervals * 2 layers) quads and each of the 8 corners in 3 * 2 * 2
        let bk = p.mesh.boundary_keys();
        assert_eq!(bk.len(), 6 * 18 + 12 * 12 + 8 * 12);
        assert_eq!(p.mesh.n_facets(), bk.len());
        for ci in p.n_physical_cells..p.mesh.n_cells() {
            assert!(crate::mesh::quality::cell_scaled_jacobian_min(&p.mesh, ci) > 0.5);
        }
        // simplex fallback is conforming and has the same volume
        let s = p.mesh.split_to_simplices().unwrap();
        s.check().unwrap();
        assert!((s.total_measure() - 1.8f64.powi(3)).abs() < 1e-11);
        assert_eq!(s.n_cells(), phys.n_cells() + 3 * c[&CellType::Prism] + 6 * c[&CellType::Hex]);
    }

    #[test]
    fn three_d_partial_layers() {
        let phys = box_tet([[0.0, 2.0], [0.0, 1.0], [0.0, 1.0]], [4, 2, 2]);
        let mut spec = PmlSpec::<3>::uniform(0.0, 0);
        spec.thickness[0][0] = 0.5;
        spec.layers[0][0] = 2;
        spec.thickness[0][1] = 0.5;
        spec.layers[0][1] = 3;
        spec.thickness[2][1] = 0.25;
        spec.layers[2][1] = 1;
        let p = build_pml(&phys, &spec).unwrap();
        p.mesh.check().unwrap();
        let s = p.mesh.split_to_simplices().unwrap();
        s.check().unwrap();
        assert!((p.mesh.total_measure() - 3.0 * 1.0 * 1.25).abs() < 1e-11);
    }
}

// ---- the nyx Pml description of a PML mesh (was mesh::interop) ----
/// The `Pml` of every domain label of a PML mesh, keyed by label, to merge into `ProblemSpec::optical_pml`.
/// `strength`, `order` and `re` are the stretch parameters of `Pml`.
pub fn pml_domains<const D: usize>(pml: &PmlMesh<D>, strength: f64, order: f64, re: f64) -> std::collections::BTreeMap<i32, crate::problem::Pml> {
    let dir = |c: i8| match c {
        0 => crate::problem::PmlDir::Off,
        c if c > 0 => crate::problem::PmlDir::Plus,
        _ => crate::problem::PmlDir::Minus,
    };
    let mut out = std::collections::BTreeMap::new();
    for r in &pml.regions {
        let axis = |a: usize| if a < D { dir(r.code[a]) } else { crate::problem::PmlDir::Off };
        for &(_, label) in &r.labels {
            // a 2D mesh is the x-z cross-section: its second axis is z
            let (y, z) = if D == 2 { (crate::problem::PmlDir::Off, axis(1)) } else { (axis(1), axis(2)) };
            out.insert(label, crate::problem::Pml { x: axis(0), y, z, strength, order, re });
        }
    }
    out
}

