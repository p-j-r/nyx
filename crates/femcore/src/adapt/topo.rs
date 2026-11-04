//! Working mesh for 2D adaptation: triangles, per-vertex metrics, and the classification of
//! vertices and edges that tells the operators what must be preserved.
//!
//! *Constrained edges* are boundary edges and edges between cells of different labels (material
//! interfaces). Connected runs of constrained edges that share a signature (the boundary-facet
//! label, or the pair of cell labels) form *chains*. A vertex is
//!
//! * `Interior`: no constrained edge;
//! * `Chain(c)`: exactly two constrained edges of the same signature that are nearly collinear
//!   (deflection below `corner_angle_deg`): it may slide along chain `c`;
//! * `Corner`: anything else: it never moves and is never removed.
//!
//! The classification and the edge-to-chain map come from the labels of the input mesh, so
//! boundaries, interfaces and corners are preserved without any geometry beyond the mesh itself.

use super::field::MetricField;
use super::metric::Sym;
use crate::mesh::native::SimplexMesh;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Small multiplicative hasher (the well-known Fx scheme) for the integer-keyed maps in the hot
/// loops; `HashMap`'s default SipHash dominated the profile.
#[derive(Default, Clone, Copy)]
pub struct FxHasher(u64);

impl Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(b as u64);
        }
    }
    fn write_u32(&mut self, i: u32) {
        self.write_u64(i as u64);
    }
    fn write_usize(&mut self, i: usize) {
        self.write_u64(i as u64);
    }
    fn write_u64(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
    fn finish(&self) -> u64 {
        self.0
    }
}

pub type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Sorted vertex pair.
pub type Edge = [u32; 2];

pub fn ekey(a: u32, b: u32) -> Edge {
    if a < b { [a, b] } else { [b, a] }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VClass {
    Interior,
    Chain(u32),
    Corner,
}

/// What a chain separates: the label of its boundary facets, or the two cell labels of an
/// interface. Unlabelled boundary edges use `Bnd(i32::MIN)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sig {
    Bnd(i32),
    Iface(i32, i32),
}

#[derive(Clone, Copy, Debug)]
pub struct ChainInfo {
    pub sig: Sig,
    /// Frozen chains keep their vertices and edges exactly (periodic pairs, symmetry planes).
    pub frozen: bool,
}

/// Options for building the classification.
#[derive(Clone, Debug)]
pub struct TopoOptions {
    /// A vertex whose two constrained edges turn by more than this many degrees is a corner.
    pub corner_angle_deg: f64,
    /// Boundary-facet labels whose chains are frozen.
    pub frozen_labels: Vec<i32>,
    /// Vertices that are corners whatever their neighbours look like: the vertices of the exact geometry (a tangent junction of a
    /// line and an arc would otherwise count as a smooth chain vertex and be free to slide off its junction).
    pub fixed_vertices: Vec<u32>,
}

impl Default for TopoOptions {
    fn default() -> Self {
        TopoOptions { corner_angle_deg: 40.0, frozen_labels: Vec::new(), fixed_vertices: Vec::new() }
    }
}

/// Edge-triangle incidence and vertex stars of the current triangles.
pub struct Snap {
    vstart: Vec<u32>,
    vlist: Vec<u32>,
    /// All edges, sorted, with the (up to two) triangles on them (`-1` for none).
    pub edges: Vec<(Edge, [i32; 2])>,
}

impl Snap {
    /// Alive triangles incident to vertex `v`, in increasing order.
    pub fn star(&self, v: u32) -> &[u32] {
        &self.vlist[self.vstart[v as usize] as usize..self.vstart[v as usize + 1] as usize]
    }

    /// Index of an edge in `edges`.
    pub fn index(&self, e: &Edge) -> Option<usize> {
        self.edges.binary_search_by(|x| x.0.cmp(e)).ok()
    }
}

#[derive(Clone, Debug)]
pub struct Work2 {
    pub pts: Vec<[f64; 2]>,
    pub met: Vec<Sym<2>>,
    pub vclass: Vec<VClass>,
    pub tris: Vec<[u32; 3]>,
    pub tlabel: Vec<i32>,
    pub talive: Vec<bool>,
    /// Triangles changed since the last flip pass: only edges next to them can have become flippable.
    pub touched: Vec<bool>,
    pub chains: Vec<ChainInfo>,
    /// Constrained edge to its chain.
    pub chain_of: FxMap<Edge, u32>,
}

pub fn tri_area(p: &[[f64; 2]], t: &[u32; 3]) -> f64 {
    let (a, b, c) = (p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]);
    0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
}

impl Work2 {
    /// Build from a mesh; the vertex metrics are `field` evaluated at the vertices.
    pub fn from_mesh(mesh: &SimplexMesh<2>, field: &dyn MetricField<2>, opt: &TopoOptions) -> Result<Work2, String> {
        mesh.check().map_err(|e| format!("input mesh invalid: {e}"))?;
        let tris: Vec<[u32; 3]> =
            (0..mesh.n_cells()).map(|c| [mesh.cell(c)[0], mesh.cell(c)[1], mesh.cell(c)[2]]).collect();
        let mut w = Work2 {
            pts: mesh.points.clone(),
            met: mesh.points.iter().map(|p| field.at(p)).collect(),
            vclass: vec![VClass::Interior; mesh.n_points()],
            talive: vec![true; tris.len()],
            touched: vec![true; tris.len()],
            tlabel: mesh.cell_label.clone(),
            tris,
            chains: Vec::new(),
            chain_of: FxMap::default(),
        };
        let snap = w.snapshot();
        let mut facet_label: HashMap<Edge, i32> = HashMap::new();
        for f in 0..mesh.n_facets() {
            let v = mesh.facet(f);
            facet_label.insert(ekey(v[0], v[1]), mesh.facet_label[f]);
        }
        // constrained edges and their signatures
        let mut sig_of: HashMap<Edge, Sig> = HashMap::new();
        for (e, ts) in &snap.edges {
            let s = if ts[1] < 0 {
                Some(Sig::Bnd(*facet_label.get(e).unwrap_or(&i32::MIN)))
            } else {
                let (l0, l1) = (w.tlabel[ts[0] as usize], w.tlabel[ts[1] as usize]);
                (l0 != l1).then(|| Sig::Iface(l0.min(l1), l0.max(l1)))
            };
            if let Some(s) = s {
                sig_of.insert(*e, s);
            }
        }
        // constrained edges at each vertex
        let mut vedges: Vec<Vec<Edge>> = vec![Vec::new(); w.pts.len()];
        for (e, _) in &snap.edges {
            if sig_of.contains_key(e) {
                vedges[e[0] as usize].push(*e);
                vedges[e[1] as usize].push(*e);
            }
        }
        let cos_lim = (opt.corner_angle_deg.to_radians()).cos();
        let mut fixed = vec![false; w.pts.len()];
        for &v in &opt.fixed_vertices {
            if let Some(f) = fixed.get_mut(v as usize) {
                *f = true;
            }
        }
        for v in 0..w.pts.len() {
            let es = &vedges[v];
            w.vclass[v] = match es.len() {
                0 => VClass::Interior,
                _ if fixed[v] => VClass::Corner,
                2 if sig_of[&es[0]] == sig_of[&es[1]] => {
                    let other = |e: &Edge| if e[0] as usize == v { e[1] } else { e[0] } as usize;
                    let (a, b) = (other(&es[0]), other(&es[1]));
                    let d1 = [w.pts[v][0] - w.pts[a][0], w.pts[v][1] - w.pts[a][1]];
                    let d2 = [w.pts[b][0] - w.pts[v][0], w.pts[b][1] - w.pts[v][1]];
                    let c = (d1[0] * d2[0] + d1[1] * d2[1]) / (d1[0].hypot(d1[1]) * d2[0].hypot(d2[1]));
                    if c >= cos_lim { VClass::Chain(0) } else { VClass::Corner }
                }
                _ => VClass::Corner,
            };
        }
        // chains: connected runs of constrained edges through Chain vertices, in edge order
        let mut cons: Vec<Edge> = sig_of.keys().cloned().collect();
        cons.sort_unstable();
        for &e in &cons {
            if w.chain_of.contains_key(&e) {
                continue;
            }
            let id = w.chains.len() as u32;
            let sig = sig_of[&e];
            w.chains.push(ChainInfo { sig, frozen: matches!(sig, Sig::Bnd(l) if opt.frozen_labels.contains(&l)) });
            let mut stack = vec![e];
            w.chain_of.insert(e, id);
            while let Some(cur) = stack.pop() {
                for &end in &cur {
                    if matches!(w.vclass[end as usize], VClass::Chain(_)) {
                        w.vclass[end as usize] = VClass::Chain(id);
                        for &n in &vedges[end as usize] {
                            if let std::collections::hash_map::Entry::Vacant(v) = w.chain_of.entry(n) {
                                v.insert(id);
                                stack.push(n);
                            }
                        }
                    }
                }
            }
        }
        Ok(w)
    }

    pub fn n_alive(&self) -> usize {
        self.talive.iter().filter(|&&a| a).count()
    }

    /// Edge and vertex incidence of the alive triangles. Linear time: vertex stars in CSR form, then
    /// for every vertex the edges to its higher-numbered neighbours are grouped with a tiny local
    /// sort, which also yields the edge list already in sorted order.
    pub fn snapshot(&self) -> Snap {
        let nv = self.pts.len();
        let mut vstart = vec![0u32; nv + 1];
        for (t, tr) in self.tris.iter().enumerate() {
            if self.talive[t] {
                for &v in tr {
                    vstart[v as usize + 1] += 1;
                }
            }
        }
        for i in 0..nv {
            vstart[i + 1] += vstart[i];
        }
        let mut fill = vstart.clone();
        let mut vlist = vec![0u32; *vstart.last().unwrap() as usize];
        for (t, tr) in self.tris.iter().enumerate() {
            if self.talive[t] {
                for &v in tr {
                    vlist[fill[v as usize] as usize] = t as u32;
                    fill[v as usize] += 1;
                }
            }
        }
        let mut edges: Vec<(Edge, [i32; 2])> = Vec::with_capacity(self.tris.len() * 3 / 2 + 8);
        let mut local: Vec<(u32, u32)> = Vec::with_capacity(32);
        for a in 0..nv as u32 {
            local.clear();
            for &t in &vlist[vstart[a as usize] as usize..vstart[a as usize + 1] as usize] {
                for &b in &self.tris[t as usize] {
                    if b > a {
                        local.push((b, t));
                    }
                }
            }
            local.sort_unstable();
            for &(b, t) in &local {
                match edges.last_mut() {
                    Some((le, s)) if le[0] == a && le[1] == b => s[1] = t as i32,
                    _ => edges.push(([a, b], [t as i32, -1])),
                }
            }
        }
        Snap { vstart, vlist, edges }
    }

    /// Internal consistency: positive areas, manifold edges, every boundary and interface edge
    /// registered in a chain and vice versa, vertex classes agree with their constrained edges.
    pub fn check(&self) -> Result<(), String> {
        let snap = self.snapshot();
        for (t, tr) in self.tris.iter().enumerate() {
            if self.talive[t] && !(tri_area(&self.pts, tr) > 0.0) {
                return Err(format!("triangle {t} has non-positive area"));
            }
        }
        for (e, ts) in &snap.edges {
            let constrained = ts[1] < 0 || self.tlabel[ts[0] as usize] != self.tlabel[ts[1] as usize];
            if constrained != self.chain_of.contains_key(e) {
                return Err(format!(
                    "edge {e:?}: constrained={constrained} but chain entry {}",
                    self.chain_of.contains_key(e)
                ));
            }
        }
        for e in self.chain_of.keys() {
            if snap.index(e).is_none() {
                return Err(format!("chain entry for missing edge {e:?}"));
            }
        }
        for (v, c) in self.vclass.iter().enumerate() {
            if snap.star(v as u32).is_empty() {
                continue;
            }
            let n = self.chain_of.keys().filter(|e| e[0] as usize == v || e[1] as usize == v).count();
            let ok = match c {
                VClass::Interior => n == 0,
                VClass::Chain(_) => n == 2,
                VClass::Corner => n >= 1,
            };
            if !ok {
                return Err(format!("vertex {v}: class {c:?} but {n} constrained edges"));
            }
        }
        Ok(())
    }

    /// Compact into a simplex mesh (alive triangles, used vertices in index order). Boundary
    /// facets carry the chain's facet label; unlabelled boundary edges produce no facet. Returns
    /// the vertex metrics in the new numbering as well.
    pub fn to_mesh(&self) -> (SimplexMesh<2>, Vec<Sym<2>>) {
        let mut used = vec![false; self.pts.len()];
        for (t, tr) in self.tris.iter().enumerate() {
            if self.talive[t] {
                for &v in tr {
                    used[v as usize] = true;
                }
            }
        }
        let mut newid = vec![u32::MAX; self.pts.len()];
        let mut out = SimplexMesh::<2>::default();
        let mut mets = Vec::new();
        for v in 0..self.pts.len() {
            if used[v] {
                newid[v] = out.points.len() as u32;
                out.points.push(self.pts[v]);
                mets.push(self.met[v]);
            }
        }
        let snap = self.snapshot();
        for (t, tr) in self.tris.iter().enumerate() {
            if self.talive[t] {
                out.push_cell(&[newid[tr[0] as usize], newid[tr[1] as usize], newid[tr[2] as usize]], self.tlabel[t]);
            }
        }
        for (e, ts) in &snap.edges {
            if ts[1] >= 0 {
                continue;
            }
            let label = self.chain_of.get(e).and_then(|&c| match self.chains[c as usize].sig {
                Sig::Bnd(l) if l != i32::MIN => Some(l),
                _ => None,
            });
            if let Some(l) = label {
                // keep the orientation of the owning triangle (counter-clockwise, outward)
                let tr = self.tris[ts[0] as usize];
                let k = (0..3).find(|&k| ekey(tr[(k + 1) % 3], tr[(k + 2) % 3]) == *e).unwrap();
                out.push_facet(&[newid[tr[(k + 1) % 3] as usize], newid[tr[(k + 2) % 3] as usize]], l);
            }
        }
        (out, mets)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapt::field::UniformSize;
    use crate::mesh::fixtures::{box_tri, side_label};

    #[test]
    fn square_classification_and_roundtrip() {
        let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 4, 4);
        let w = Work2::from_mesh(&m, &UniformSize(0.2), &TopoOptions::default()).unwrap();
        w.check().unwrap();
        let corners = w.vclass.iter().filter(|c| **c == VClass::Corner).count();
        let chain = w.vclass.iter().filter(|c| matches!(c, VClass::Chain(_))).count();
        assert_eq!(corners, 4);
        assert_eq!(chain, 4 * 3);
        assert_eq!(w.chains.len(), 4);
        let sigs: Vec<Sig> = w.chains.iter().map(|c| c.sig).collect();
        for a in 0..2 {
            for s in [false, true] {
                assert!(sigs.contains(&Sig::Bnd(side_label(a, s))));
            }
        }
        let (back, mets) = w.to_mesh();
        back.check().unwrap();
        assert_eq!(back.n_cells(), m.n_cells());
        assert_eq!(back.n_facets(), m.n_facets());
        assert_eq!(mets.len(), m.n_points());
        assert_eq!(back.points, m.points);
        assert_eq!(back.cells, m.cells);
        let mut a: Vec<_> = (0..m.n_facets()).map(|f| (crate::mesh::native::facet_key(m.facet(f)), m.facet_label[f])).collect();
        let mut b: Vec<_> =
            (0..back.n_facets()).map(|f| (crate::mesh::native::facet_key(back.facet(f)), back.facet_label[f])).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b);
    }

    #[test]
    fn interface_chain_and_tjunctions() {
        let mut m = box_tri([[0.0, 2.0], [0.0, 1.0]], 8, 4);
        for c in 0..m.n_cells() {
            let cx: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][0]).sum::<f64>() / 3.0;
            m.cell_label[c] = if cx < 1.0 { 1 } else { 2 };
        }
        let w = Work2::from_mesh(&m, &UniformSize(0.2), &TopoOptions::default()).unwrap();
        w.check().unwrap();
        let iface: Vec<usize> = (0..w.chains.len()).filter(|&c| matches!(w.chains[c].sig, Sig::Iface(1, 2))).collect();
        assert_eq!(iface.len(), 1);
        // interface vertices: 3 interior ones slide, the two ends meet the boundary: corners
        let on_if: Vec<VClass> =
            (0..w.pts.len()).filter(|&v| (w.pts[v][0] - 1.0).abs() < 1e-12).map(|v| w.vclass[v]).collect();
        assert_eq!(on_if.len(), 5);
        assert_eq!(on_if.iter().filter(|c| **c == VClass::Corner).count(), 2);
        assert_eq!(on_if.iter().filter(|c| **c == VClass::Chain(iface[0] as u32)).count(), 3);
    }

    #[test]
    fn frozen_labels_mark_chains() {
        let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 3, 3);
        let opt = TopoOptions { frozen_labels: vec![side_label(0, false)], ..Default::default() };
        let w = Work2::from_mesh(&m, &UniformSize(0.3), &opt).unwrap();
        assert_eq!(w.chains.iter().filter(|c| c.frozen).count(), 1);
    }
}
