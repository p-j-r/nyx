//! Working mesh for 3D adaptation: tetrahedra, per-vertex metrics, and the classification of vertices that tells the operators what must
//! be preserved (the 3D counterpart of `topo`).
//!
//! *Surface faces* are the faces of one tetrahedron (the boundary) and the faces between tetrahedra of different labels (material
//! interfaces). Their *signature* (`Sig`) is the label of the boundary facet, or the pair of cell labels. Surface faces of one signature
//! that are connected through edges where the surface does not fold form a *patch*. A surface edge is a *crease* when the surface folds
//! there (the outward normals of its two faces differ by more than `crease_angle_deg`), when the two faces have different signatures, or
//! when the edge does not have exactly two surface faces (three materials meeting). Connected runs of creases of one signature form
//! *curves*. A vertex is
//!
//! * `Interior`: no surface face;
//! * `Surface(p)`: on patch `p` with no crease edge, free to slide over it;
//! * `Crease(c)`: exactly two crease edges of the same signature that are nearly collinear (deflection below `corner_angle_deg`), free to
//!   slide along curve `c`;
//! * `Corner`: anything else; it never moves and is never removed.
//!
//! Everything comes from the labels and the shape of the input mesh, so boundaries, interfaces, creases and corners are preserved
//! without any geometry beyond the mesh itself; `geom3` adds exact projectors for smooth surfaces.

use super::field::MetricField;
use super::metric::Sym;
use super::topo::{ekey, Edge, FxMap, Sig};
use crate::adapt::work::{facet_key, SimplexMesh};
use std::collections::{BTreeSet, HashMap};

/// Sorted vertex triple.
pub type Face = [u32; 3];

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VClass3 {
    Interior,
    Surface(u32),
    Crease(u32),
    Corner,
}

/// A connected piece of surface of one signature.
#[derive(Clone, Copy, Debug)]
pub struct PatchInfo {
    pub sig: Sig,
    pub frozen: bool,
}

/// A connected run of crease edges.
#[derive(Clone, Debug)]
pub struct CurveInfo {
    /// signatures of the surface faces along the curve, sorted
    pub sigs: Vec<Sig>,
    pub frozen: bool,
}

/// Options for building the classification.
#[derive(Clone, Debug)]
pub struct TopoOptions3 {
    /// A surface edge whose two faces turn by more than this many degrees is a crease.
    pub crease_angle_deg: f64,
    /// A vertex whose two crease edges turn by more than this many degrees is a corner.
    pub corner_angle_deg: f64,
    /// Boundary-facet labels whose patches (and the curves around them) are frozen.
    pub frozen_labels: Vec<i32>,
    /// Vertices that are corners whatever their neighbours look like.
    pub fixed_vertices: Vec<u32>,
}

impl Default for TopoOptions3 {
    fn default() -> Self {
        TopoOptions3 { crease_angle_deg: 40.0, corner_angle_deg: 40.0, frozen_labels: Vec::new(), fixed_vertices: Vec::new() }
    }
}

#[derive(Clone, Debug)]
pub struct Work3 {
    pub pts: Vec<[f64; 3]>,
    pub met: Vec<Sym<3>>,
    pub vclass: Vec<VClass3>,
    pub tets: Vec<[u32; 4]>,
    pub tlabel: Vec<i32>,
    pub talive: Vec<bool>,
    pub patches: Vec<PatchInfo>,
    pub curves: Vec<CurveInfo>,
    /// Surface face to its patch.
    pub patch_of: FxMap<Face, u32>,
    /// Crease edge to its curve.
    pub curve_of: FxMap<Edge, u32>,
}

pub fn fkey(a: u32, b: u32, c: u32) -> Face {
    let mut f = [a, b, c];
    f.sort_unstable();
    f
}

pub fn tet_volume(p: &[[f64; 3]], t: &[u32; 4]) -> f64 {
    let (a, b, c, d) = (p[t[0] as usize], p[t[1] as usize], p[t[2] as usize], p[t[3] as usize]);
    let (u, v, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]], [d[0] - a[0], d[1] - a[1], d[2] - a[2]]);
    (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = dot(a, a).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

/// The four faces of a tetrahedron, each with the opposite vertex.
pub fn tet_faces(t: &[u32; 4]) -> [(Face, u32); 4] {
    [(fkey(t[1], t[2], t[3]), t[0]), (fkey(t[0], t[2], t[3]), t[1]), (fkey(t[0], t[1], t[3]), t[2]), (fkey(t[0], t[1], t[2]), t[3])]
}

/// One surface face found while classifying: its vertices, signature and unit normal pointing out of the lower-labelled side (out of
/// the domain on a boundary).
struct SFace {
    key: Face,
    sig: Sig,
    n: [f64; 3],
}

impl Work3 {
    /// Build from a mesh; the vertex metrics are `field` evaluated at the vertices.
    pub fn from_mesh(mesh: &SimplexMesh<3>, field: &dyn MetricField<3>, opt: &TopoOptions3) -> Result<Work3, String> {
        mesh.check().map_err(|e| format!("input mesh invalid: {e}"))?;
        let tets: Vec<[u32; 4]> = (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2], mesh.cell(c)[3]]).collect();
        let mut w = Work3 {
            pts: mesh.points.clone(),
            met: mesh.points.iter().map(|p| field.at(p)).collect(),
            vclass: vec![VClass3::Interior; mesh.n_points()],
            talive: vec![true; tets.len()],
            tlabel: mesh.cell_label.clone(),
            tets,
            patches: vec![],
            curves: vec![],
            patch_of: FxMap::default(),
            curve_of: FxMap::default(),
        };
        // faces with the tetrahedra on them
        let mut on_face: HashMap<Face, Vec<(usize, u32)>> = HashMap::new();
        for (t, tt) in w.tets.iter().enumerate() {
            for (f, opp) in tet_faces(tt) {
                on_face.entry(f).or_default().push((t, opp));
            }
        }
        let mut facet_label: HashMap<Face, i32> = HashMap::new();
        for f in 0..mesh.n_facets() {
            let v = mesh.facet(f);
            facet_label.insert(facet_key(&[v[0], v[1], v[2]]).iter().take(3).copied().collect::<Vec<_>>().try_into().unwrap(), mesh.facet_label[f]);
        }
        // surface faces: signature and outward normal
        let mut sfaces: Vec<SFace> = vec![];
        let mut keys: Vec<&Face> = on_face.keys().collect();
        keys.sort_unstable();
        for k in keys {
            let ts = &on_face[k];
            let (sig, owner_opp) = match ts.len() {
                1 => (Sig::Bnd(*facet_label.get(k).unwrap_or(&i32::MIN)), ts[0].1),
                2 => {
                    let (l0, l1) = (w.tlabel[ts[0].0], w.tlabel[ts[1].0]);
                    if l0 == l1 {
                        continue;
                    }
                    // the normal points out of the lower-labelled side
                    (Sig::Iface(l0.min(l1), l0.max(l1)), if l0 < l1 { ts[0].1 } else { ts[1].1 })
                }
                _ => return Err(format!("face {k:?} shared by {} tetrahedra", ts.len())),
            };
            let (a, b, c) = (w.pts[k[0] as usize], w.pts[k[1] as usize], w.pts[k[2] as usize]);
            let mut n = unit(cross(sub(b, a), sub(c, a)));
            if dot(n, sub(w.pts[owner_opp as usize], a)) > 0.0 {
                n = [-n[0], -n[1], -n[2]];
            }
            sfaces.push(SFace { key: *k, sig, n });
        }
        // surface edges and the faces on them
        let mut on_edge: HashMap<Edge, Vec<usize>> = HashMap::new();
        for (i, f) in sfaces.iter().enumerate() {
            for (a, b) in [(f.key[0], f.key[1]), (f.key[0], f.key[2]), (f.key[1], f.key[2])] {
                on_edge.entry(ekey(a, b)).or_default().push(i);
            }
        }
        let cos_crease = opt.crease_angle_deg.to_radians().cos();
        let mut crease_sigs: HashMap<Edge, Vec<Sig>> = HashMap::new();
        let mut eks: Vec<&Edge> = on_edge.keys().collect();
        eks.sort_unstable();
        for e in eks {
            let fs = &on_edge[e];
            let fold = fs.len() != 2 || sfaces[fs[0]].sig != sfaces[fs[1]].sig || dot(sfaces[fs[0]].n, sfaces[fs[1]].n) < cos_crease;
            if fold {
                let mut sg: Vec<Sig> = fs.iter().map(|&i| sfaces[i].sig).collect();
                sg.sort_by_key(|s| format!("{s:?}"));
                sg.dedup();
                crease_sigs.insert(*e, sg);
            }
        }
        // patches: union of the faces of one signature through non-crease edges
        let mut parent: Vec<usize> = (0..sfaces.len()).collect();
        fn find(p: &mut [usize], mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        for (e, fs) in &on_edge {
            if !crease_sigs.contains_key(e) {
                let (a, b) = (find(&mut parent, fs[0]), find(&mut parent, fs[1]));
                if a != b {
                    parent[a.max(b)] = a.min(b);
                }
            }
        }
        let mut patch_id: HashMap<usize, u32> = HashMap::new();
        for i in 0..sfaces.len() {
            let r = find(&mut parent, i);
            let id = *patch_id.entry(r).or_insert_with(|| {
                let sig = sfaces[r].sig;
                w.patches.push(PatchInfo { sig, frozen: matches!(sig, Sig::Bnd(l) if opt.frozen_labels.contains(&l)) });
                (w.patches.len() - 1) as u32
            });
            w.patch_of.insert(sfaces[i].key, id);
        }
        // vertex classification
        let mut vfaces: Vec<Vec<usize>> = vec![vec![]; w.pts.len()];
        for (i, f) in sfaces.iter().enumerate() {
            for &v in &f.key {
                vfaces[v as usize].push(i);
            }
        }
        let mut vcrease: Vec<Vec<Edge>> = vec![vec![]; w.pts.len()];
        let mut ces: Vec<&Edge> = crease_sigs.keys().collect();
        ces.sort_unstable();
        for e in ces {
            vcrease[e[0] as usize].push(*e);
            vcrease[e[1] as usize].push(*e);
        }
        let cos_corner = opt.corner_angle_deg.to_radians().cos();
        let mut fixed = vec![false; w.pts.len()];
        for &v in &opt.fixed_vertices {
            if let Some(f) = fixed.get_mut(v as usize) {
                *f = true;
            }
        }
        for v in 0..w.pts.len() {
            let es = &vcrease[v];
            w.vclass[v] = if vfaces[v].is_empty() {
                VClass3::Interior
            } else if fixed[v] {
                VClass3::Corner
            } else if es.is_empty() {
                VClass3::Surface(w.patch_of[&sfaces[vfaces[v][0]].key])
            } else if es.len() == 2 && crease_sigs[&es[0]] == crease_sigs[&es[1]] {
                let other = |e: &Edge| if e[0] as usize == v { e[1] } else { e[0] } as usize;
                let (a, b) = (other(&es[0]), other(&es[1]));
                let (d1, d2) = (unit(sub(w.pts[v], w.pts[a])), unit(sub(w.pts[b], w.pts[v])));
                if dot(d1, d2) >= cos_corner { VClass3::Crease(0) } else { VClass3::Corner }
            } else {
                VClass3::Corner
            };
        }
        // curves: connected runs of crease edges through Crease vertices
        let mut cons: Vec<Edge> = crease_sigs.keys().cloned().collect();
        cons.sort_unstable();
        for &e in &cons {
            if w.curve_of.contains_key(&e) {
                continue;
            }
            let id = w.curves.len() as u32;
            let sigs = crease_sigs[&e].clone();
            let frozen = sigs.iter().any(|s| matches!(s, Sig::Bnd(l) if opt.frozen_labels.contains(l)));
            w.curves.push(CurveInfo { sigs, frozen });
            let mut stack = vec![e];
            w.curve_of.insert(e, id);
            while let Some(cur) = stack.pop() {
                for &end in &cur {
                    if matches!(w.vclass[end as usize], VClass3::Crease(_)) {
                        w.vclass[end as usize] = VClass3::Crease(id);
                        for &n in &vcrease[end as usize] {
                            if let std::collections::hash_map::Entry::Vacant(x) = w.curve_of.entry(n) {
                                x.insert(id);
                                stack.push(n);
                            }
                        }
                    }
                }
            }
        }
        Ok(w)
    }

    /// Drops the dead tetrahedra from the arrays. Operators append the tetrahedra they make and only mark the old ones dead, so without
    /// this every pass would walk an array several times the size of the mesh. Only the order of the surviving tetrahedra is kept, which is
    /// all the passes depend on.
    pub fn compact_tets(&mut self) {
        let mut k = 0;
        for t in 0..self.tets.len() {
            if self.talive[t] {
                self.tets[k] = self.tets[t];
                self.tlabel[k] = self.tlabel[t];
                k += 1;
            }
        }
        self.tets.truncate(k);
        self.tlabel.truncate(k);
        self.talive.clear();
        self.talive.resize(k, true);
    }

    pub fn n_alive(&self) -> usize {
        self.talive.iter().filter(|&&a| a).count()
    }

    /// Internal consistency: positive volumes, manifold faces, every surface face registered in a patch and vice versa, vertex classes
    /// agree with their crease edges.
    pub fn check(&self) -> Result<(), String> {
        let mut on_face: HashMap<Face, Vec<usize>> = HashMap::new();
        for (t, tt) in self.tets.iter().enumerate() {
            if !self.talive[t] {
                continue;
            }
            if !(tet_volume(&self.pts, tt) > 0.0) {
                return Err(format!("tetrahedron {t} has non-positive volume"));
            }
            for (f, _) in tet_faces(tt) {
                on_face.entry(f).or_default().push(t);
            }
        }
        for (f, ts) in &on_face {
            if ts.len() > 2 {
                return Err(format!("face {f:?} shared by {} tetrahedra", ts.len()));
            }
            // the two tetrahedra of an interior face lie on opposite sides of it (otherwise they overlap)
            if ts.len() == 2 {
                let opp = |t: usize| *self.tets[t].iter().find(|x| !f.contains(x)).unwrap();
                let side = |t: usize| tet_volume(&self.pts, &[f[0], f[1], f[2], opp(t)]);
                let (s0, s1) = (side(ts[0]), side(ts[1]));
                if s0 * s1 > 0.0 {
                    return Err(format!("tetrahedra {} and {} overlap across face {f:?}", ts[0], ts[1]));
                }
            }
            let surface = ts.len() == 1 || self.tlabel[ts[0]] != self.tlabel[ts[1]];
            if surface != self.patch_of.contains_key(f) {
                return Err(format!("face {f:?}: surface={surface} but patch entry {}", self.patch_of.contains_key(f)));
            }
        }
        for f in self.patch_of.keys() {
            if !on_face.contains_key(f) {
                return Err(format!("patch entry for missing face {f:?}"));
            }
        }
        let mut n_crease = vec![0usize; self.pts.len()];
        for e in self.curve_of.keys() {
            n_crease[e[0] as usize] += 1;
            n_crease[e[1] as usize] += 1;
        }
        let mut on_surface = vec![false; self.pts.len()];
        for f in self.patch_of.keys() {
            for &v in f {
                on_surface[v as usize] = true;
            }
        }
        for (v, c) in self.vclass.iter().enumerate() {
            let ok = match c {
                VClass3::Interior => !on_surface[v] && n_crease[v] == 0,
                VClass3::Surface(_) => on_surface[v] && n_crease[v] == 0,
                VClass3::Crease(_) => n_crease[v] == 2,
                VClass3::Corner => on_surface[v],
            };
            if !ok && on_face.keys().any(|f| f.contains(&(v as u32))) {
                return Err(format!("vertex {v}: class {c:?} but {} crease edges, on surface {}", n_crease[v], on_surface[v]));
            }
        }
        Ok(())
    }

    /// Compact into a simplex mesh (alive tetrahedra, used vertices in index order). Boundary facets carry the patch's facet label;
    /// unlabelled boundary faces produce no facet. Returns the vertex metrics in the new numbering as well.
    pub fn to_mesh(&self) -> (SimplexMesh<3>, Vec<Sym<3>>) {
        let mut used = vec![false; self.pts.len()];
        for (t, tt) in self.tets.iter().enumerate() {
            if self.talive[t] {
                for &v in tt {
                    used[v as usize] = true;
                }
            }
        }
        let mut new_id = vec![u32::MAX; self.pts.len()];
        let mut m = SimplexMesh::<3>::default();
        let mut mets = vec![];
        for v in 0..self.pts.len() {
            if used[v] {
                new_id[v] = m.points.len() as u32;
                m.points.push(self.pts[v]);
                mets.push(self.met[v]);
            }
        }
        for (t, tt) in self.tets.iter().enumerate() {
            if self.talive[t] {
                m.push_cell(&[new_id[tt[0] as usize], new_id[tt[1] as usize], new_id[tt[2] as usize], new_id[tt[3] as usize]], self.tlabel[t]);
            }
        }
        let mut faces: Vec<(&Face, &u32)> = self.patch_of.iter().collect();
        faces.sort_unstable();
        let mut seen = BTreeSet::new();
        for (f, &p) in faces {
            if let Sig::Bnd(l) = self.patches[p as usize].sig {
                if l != i32::MIN && seen.insert(*f) {
                    // orient like a face of its owner so that outward normals are consistent
                    let owner = self.tets.iter().enumerate().find(|(t, tt)| self.talive[*t] && tet_faces(tt).iter().any(|(k, _)| k == f));
                    if let Some((_, tt)) = owner {
                        let (_, opp) = *tet_faces(tt).iter().find(|(k, _)| k == f).unwrap();
                        let mut v = [new_id[f[0] as usize], new_id[f[1] as usize], new_id[f[2] as usize]];
                        let (a, b, c) = (self.pts[f[0] as usize], self.pts[f[1] as usize], self.pts[f[2] as usize]);
                        if dot(cross(sub(b, a), sub(c, a)), sub(self.pts[opp as usize], a)) > 0.0 {
                            v.swap(1, 2);
                        }
                        m.push_facet(&v, l);
                    }
                }
            }
        }
        (m, mets)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapt::field::UniformSize;
    use crate::mesh::fixtures::box_tet;

    fn count(w: &Work3) -> (usize, usize, usize, usize) {
        let mut c = (0, 0, 0, 0);
        for v in &w.vclass {
            match v {
                VClass3::Interior => c.0 += 1,
                VClass3::Surface(_) => c.1 += 1,
                VClass3::Crease(_) => c.2 += 1,
                VClass3::Corner => c.3 += 1,
            }
        }
        c
    }

    #[test]
    fn a_box_has_corners_creases_faces_and_an_interior() {
        let n = 4;
        let m = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [n, n, n]);
        let w = Work3::from_mesh(&m, &UniformSize(0.25), &TopoOptions3::default()).unwrap();
        w.check().unwrap();
        // 8 corners, 12 edges of n - 1 free vertices, 6 faces of (n - 1)^2, (n - 1)^3 inside
        assert_eq!(count(&w), ((n - 1).pow(3), 6 * (n - 1).pow(2), 12 * (n - 1), 8));
        assert_eq!(w.patches.len(), 6);
        assert_eq!(w.curves.len(), 12);
    }

    #[test]
    fn an_interface_inside_the_box_is_a_patch_of_its_own() {
        let mut m = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [4, 4, 4]);
        for c in 0..m.n_cells() {
            let z: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][2]).sum::<f64>() / 4.0;
            if z > 0.5 {
                m.cell_label[c] = 2;
            }
        }
        let w = Work3::from_mesh(&m, &UniformSize(0.25), &TopoOptions3::default()).unwrap();
        w.check().unwrap();
        // the plane z = 0.5 meets the four side faces in creases of a different signature than the box edges: 4 sides split in two
        let iface: Vec<_> = w.patches.iter().filter(|p| matches!(p.sig, Sig::Iface(1, 2))).collect();
        assert_eq!(iface.len(), 1);
        // 5 vertices per side line on the interface plane, minus the 4 corners of the square, are crease vertices
        let creases_on_plane = (0..w.pts.len()).filter(|&v| matches!(w.vclass[v], VClass3::Crease(_)) && (w.pts[v][2] - 0.5).abs() < 1e-12).count();
        assert_eq!(creases_on_plane, (4 * 3));
        let corners_on_plane = (0..w.pts.len()).filter(|&v| matches!(w.vclass[v], VClass3::Corner) && (w.pts[v][2] - 0.5).abs() < 1e-12).count();
        assert_eq!(corners_on_plane, 4);
    }

    #[test]
    fn round_trip_keeps_cells_and_labelled_facets() {
        let m = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [3, 3, 3]);
        let w = Work3::from_mesh(&m, &UniformSize(0.3), &TopoOptions3::default()).unwrap();
        let (back, mets) = w.to_mesh();
        back.check().unwrap();
        assert_eq!(back.n_cells(), m.n_cells());
        assert_eq!(back.n_facets(), m.n_facets());
        assert_eq!(mets.len(), back.n_points());
        assert!((back.total_measure() - 1.0).abs() < 1e-12);
    }
}
