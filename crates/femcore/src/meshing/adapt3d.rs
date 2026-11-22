//! Adaptive meshing of a 3D `geom` model (ADR 51): the metric adaptation of `adapt::drive3` driven by the exact faces of the model.
//!
//! The 2D adaptive mesher starts from the simple mesh, builds a metric from what the model asks for (region sizes, curvature) and what the
//! physics asks for (wavelengths), and adapts. This is the same in 3D with the tetrahedral lattice mesh as the start mesh and `adapt3` as
//! the engine. Three things are specific to 3D and live here:
//!
//! * **The oracle** (`GeomOracle3`): a vertex that slides over a boundary or interface patch is projected onto the *exact* geom face
//!   (`Complex::facet_project`: plane, sphere, torus, cylinder wall) instead of the triangles of the start mesh, so refining a sphere
//!   converges to the sphere. Creases and corners keep the polylines of the start mesh (`FacetGeometry`): straight creases are exact,
//!   curved ones are as good as the start mesh, like the lattice mesher itself (ADR 40).
//! * **Interface facets.** `Work3::to_mesh` returns the facets of the outer boundary only. The mesher of ADR 40 also writes the faces
//!   between regions as labelled facets, so they are written back here, labelled by the nearest geom face between exactly those two
//!   regions (the rule of `GeomDomain3::boundary_label`).
//! * **Sizes.** `vertex_sizes3` is the 3D counterpart of the 2D `vertex_sizes`: the size of the region of every cell, the sizes at which a
//!   chord of a curved face stays within `max_dev`; `geometry_floor3` is the part the geometry insists on whatever a solution does.

use crate::adapt::drive::AdaptReport;
use crate::adapt::drive3::adapt3_with;
use crate::adapt::field3::BackgroundMetric3;
use crate::adapt::geom3::{FacetGeometry, Geometry3};
use crate::adapt::metric::{isotropic, Sym};
use crate::adapt::ops3::AdaptOptions3;
use crate::adapt::topo::Sig;
use crate::adapt::topo3::{TopoOptions3, Work3};
use crate::adapt::work::{simplex_to_adapt, SimplexMesh as WorkMesh};
use crate::mesh::SimplexMesh;
use crate::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use crate::meshing::geomdrive3::GeomDomain3;
use crate::meshing::Adaptive;
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use std::collections::{BTreeMap, HashMap};

/// The faces of a model as the adaptation needs them: the domain numbers on their two sides, their boundary number and their box.
pub struct FaceTable<'a> {
    ev: &'a Evaluated,
    /// domain numbers of the regions on the two sides of every face (`i32::MIN`: nothing)
    sides: Vec<(i32, i32)>,
    label: Vec<i32>,
    boxes: Vec<geom::math::Aabb>,
}

impl<'a> FaceTable<'a> {
    pub fn new(ev: &'a Evaluated, labels: &LabelMap) -> Result<Self, String> {
        let cx = ev.complex();
        let nf = cx.facet_count();
        let dom = |r: Option<u32>| r.map_or(i32::MIN, |r| labels.domains[r as usize].id);
        let sides = (0..nf as u32).map(|i| { let (a, b) = cx.facet_sides(i); (dom(a), dom(b)) }).collect();
        let boxes = (0..nf as u32).map(|i| cx.facet_bbox(i)).collect();
        let label = (0..nf as u32)
            .map(|i| labels.label_of_edge(i).ok_or_else(|| format!("face '{}' has no boundary label", cx.facet_name(i))))
            .collect::<Result<Vec<i32>, String>>()?;
        Ok(FaceTable { ev, sides, label, boxes })
    }

    /// Distance from `p` to face `i` and the closest point; the box of the face gives a lower bound that skips the projection.
    fn distance(&self, i: usize, p: [f64; 3], cap: f64) -> (f64, Option<[f64; 3]>) {
        let b = &self.boxes[i];
        let d = (0..3).map(|k| (b.lo[k] - p[k]).max(p[k] - b.hi[k]).max(0.0).powi(2)).sum::<f64>().sqrt();
        if d >= cap {
            return (d, None);
        }
        let q = self.ev.complex().facet_project(i as u32, p);
        (((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2)).sqrt(), Some(q))
    }

    fn closest_of(&self, p: [f64; 3], cands: impl Iterator<Item = usize>) -> Option<[f64; 3]> {
        let mut best = (f64::INFINITY, None);
        for i in cands {
            let (d, q) = self.distance(i, p, best.0);
            if let (true, Some(q)) = (d < best.0, q) {
                best = (d, Some(q));
            }
        }
        best.1
    }

    /// The boundary number of the nearest of the faces `cands` to `p`.
    pub fn nearest_of(&self, p: [f64; 3], cands: impl Iterator<Item = usize>) -> Option<i32> {
        let mut best = (f64::INFINITY, None);
        for i in cands {
            let (d, _) = self.distance(i, p, best.0);
            if d < best.0 {
                best = (d, Some(self.label[i]));
            }
        }
        best.1
    }

    /// Faces a boundary patch of label `l` may lie on.
    fn faces_of_label(&self, l: i32) -> Vec<usize> {
        (0..self.label.len()).filter(|&i| self.label[i] == l).collect()
    }

    /// Faces between exactly the domains `a` and `b`.
    fn faces_between(&self, a: i32, b: i32) -> Vec<usize> {
        (0..self.sides.len()).filter(|&i| self.sides[i] == (a, b) || self.sides[i] == (b, a)).collect()
    }

    /// The boundary number of the face a facet with centroid `p` between the domains `a` and `b` lies on: the nearest face between exactly
    /// those two domains, else the nearest of all.
    pub fn interface_label(&self, p: [f64; 3], a: i32, b: i32) -> i32 {
        let mut best = (f64::INFINITY, None::<usize>);
        for pass in 0..2 {
            for i in 0..self.sides.len() {
                let s = self.sides[i];
                if pass == 0 && !(s == (a, b) || s == (b, a)) {
                    continue;
                }
                let (d, _) = self.distance(i, p, best.0);
                if d < best.0 {
                    best = (d, Some(i));
                }
            }
            if best.1.is_some() {
                break;
            }
        }
        best.1.map_or(1, |i| self.label[i])
    }
}

/// The exact geometry as the oracle of the 3D adaptation.
pub struct GeomOracle3<'a> {
    table: &'a FaceTable<'a>,
    /// candidate faces of every patch of the working mesh
    cands: Vec<Vec<usize>>,
    base: FacetGeometry,
}

impl<'a> GeomOracle3<'a> {
    pub fn from_work(table: &'a FaceTable<'a>, w: &Work3) -> Self {
        let cands = w
            .patches
            .iter()
            .map(|p| match p.sig {
                Sig::Bnd(l) => table.faces_of_label(l),
                Sig::Iface(a, b) => table.faces_between(a, b),
            })
            .collect();
        GeomOracle3 { table, cands, base: FacetGeometry::from_work(w) }
    }
}

impl Geometry3 for GeomOracle3<'_> {
    fn project_surface(&self, p: [f64; 3], patch: u32) -> [f64; 3] {
        let c = &self.cands[patch as usize];
        if c.is_empty() || std::env::var_os("NYX_NO_ORACLE").is_some() {
            return self.base.project_surface(p, patch);
        }
        self.table.closest_of(p, c.iter().copied()).unwrap_or(p)
    }
    fn project_crease(&self, p: [f64; 3], curve: u32) -> [f64; 3] {
        self.base.project_crease(p, curve)
    }
}

/// Faces between cells of different domains, labelled, oriented out of the lower-numbered domain.
pub fn interface_facets(point: &dyn Fn(u32) -> [f64; 3], cells: &[u32], cell_label: &[i32], table: &FaceTable) -> (Vec<[u32; 3]>, Vec<i32>) {
    interface_facets_with_preserved(point, cells, cell_label, table, &BTreeMap::new())
}

/// Like [`interface_facets`], while also emitting surviving physical interior facets recorded in `preserved`
/// (such as `Op::Sheet` airbridges where `la == lb`, keeping their exact vertex orientation and `LabelMap` boundary id).
pub fn interface_facets_with_preserved(
    point: &dyn Fn(u32) -> [f64; 3],
    cells: &[u32],
    cell_label: &[i32],
    table: &FaceTable,
    preserved: &BTreeMap<[u32; 3], ([u32; 3], i32)>,
) -> (Vec<[u32; 3]>, Vec<i32>) {
    let mut faces: BTreeMap<[u32; 3], Vec<(usize, u32)>> = BTreeMap::new();
    for c in 0..cell_label.len() {
        let v = &cells[4 * c..4 * c + 4];
        for skip in 0..4 {
            let f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            let mut k = [f[0], f[1], f[2]];
            k.sort_unstable();
            faces.entry(k).or_default().push((c, v[skip]));
        }
    }
    let (mut out, mut labels) = (vec![], vec![]);
    for (k, cs) in faces {
        if cs.len() != 2 {
            continue;
        }
        if let Some(&(tri, lab)) = preserved.get(&k) {
            out.push(tri);
            labels.push(lab);
            continue;
        }
        let (la, lb) = (cell_label[cs[0].0], cell_label[cs[1].0]);
        if la == lb {
            continue;
        }
        let (owner, opp) = if la < lb { cs[0] } else { cs[1] };
        let _ = owner;
        let p = point;
        let (a, b, c) = (p(k[0]), p(k[1]), p(k[2]));
        let centroid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0];
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let o = p(opp);
        let mut tri = k;
        if n[0] * (o[0] - a[0]) + n[1] * (o[1] - a[1]) + n[2] * (o[2] - a[2]) > 0.0 {
            tri.swap(1, 2);
        }
        out.push(tri);
        labels.push(table.interface_label(centroid, la, lb));
    }
    (out, labels)
}

/// Adapts a meshed 3D model to a metric given at every vertex of its mesh: graded, the faces of the model as the oracle, the labels of
/// cells and of boundary and interface facets carried through.
pub fn adapt_to_metrics3(g: MeshedGeom, metrics: &[Sym<3>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    if g.mesh.dim != 3 {
        return Err("3D adaptation works on tetrahedral meshes".into());
    }
    if metrics.len() != g.mesh.n_points() {
        return Err(format!("{} metrics for {} vertices", metrics.len(), g.mesh.n_points()));
    }
    let mut old = simplex_to_adapt::<3>(&g.mesh)?;
    // the adaptation keeps the outer boundary only; the facets between regions are written back from the geometry afterwards
    drop_interior_facets(&mut old);
    let field = BackgroundMetric3::new_graded(&old, metrics, a.grade);
    // Against a curved surface the last splits and collapses never cancel out (the surface cannot conform exactly), so the default of 30
    // iterations to a 0.5 % change rate is spent on churn: stop at 2 %, after at most 12 iterations
    let opt = AdaptOptions3 { max_cells: a.max_cells, max_iter: 12, converge_frac: 0.02, sweep_frac: 0.003, ..AdaptOptions3::default() };
    // The geometry is known exactly, so the dihedral angle of a coarse facetted sphere is no reason to call its edges creases (they would
    // pin the vertices to the chords): faces of different labels still meet in creases whatever the angle, and a vertex of a patch is
    // projected onto the nearest exact face of its label.
    let topo = TopoOptions3 { crease_angle_deg: 120.0, frozen_labels: crate::meshing::periodic_labels(&g), ..TopoOptions3::default() };
    let table = FaceTable::new(&g.eval, &g.labels)?;
    let result = adapt3_with(&old, &field, &topo, &opt, |w| GeomOracle3::from_work(&table, w))?;
    let wm = &result.mesh;
    let (extra, extra_labels) = interface_facets(&|v| wm.points[v as usize], &wm.cells, &wm.cell_label, &table);
    let mut facets = wm.facets.clone();
    let mut facet_labels = wm.facet_label.clone();
    for (t, l) in extra.iter().zip(&extra_labels) {
        facets.extend_from_slice(t);
        facet_labels.push(*l);
    }
    let points = wm.points.iter().flat_map(|p| p.iter().copied()).collect();
    let mesh = SimplexMesh::new(3, points, wm.cells.clone(), wm.cell_label.clone(), facets, facet_labels);
    Ok((MeshedGeom { mesh, labels: g.labels, eval: g.eval, warnings: g.warnings }, result.report))
}

/// The element size every vertex of a meshed 3D model is asked for by the model itself: the size of the regions it belongs to (the
/// global size lowered by `Size` rules on regions and by the layers of an absorbing piece), and on a curved face the size at which a
/// chord stays within `max_dev` of the face. All times `a.scale`. With `geometry_only` a region that asks for no more than the global
/// size contributes the global size, which makes the result a floor for an estimator metric.
pub fn vertex_sizes3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive, geometry_only: bool) -> Result<Vec<f64>, String> {
    if !(a.scale > 0.0 && a.scale.is_finite()) {
        return Err("the adaptation scale must be positive".into());
    }
    let m = &g.mesh;
    let dom = GeomDomain3::new(&g.eval, &g.labels, opts.h)?;
    let region_of: HashMap<i32, usize> = g.labels.domains.iter().map(|d| (d.id, d.region as usize)).collect();
    let mut vh = vec![f64::INFINITY; m.n_points()];
    for c in 0..m.n_cells() {
        let mut h = region_of.get(&m.cell_labels[c]).map_or(opts.h, |&r| dom.region_sizes()[r]);
        if geometry_only && h >= opts.h * 0.999 {
            h = opts.h;
        }
        for &v in m.cell(c) {
            vh[v as usize] = vh[v as usize].min(h);
        }
    }
    // vertices on a boundary or an interface facet: the only ones a curvature limit concerns
    let mut on_facet = vec![false; m.n_points()];
    for f in 0..m.n_facets() {
        for &v in m.facet(f) {
            on_facet[v as usize] = true;
        }
    }
    let (lo, hi) = m.bbox();
    let diag = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let tol = 1e-6 * diag.max(opts.h * 1e-3);
    let max_dev = opts.max_dev.unwrap_or(opts.h / 50.0);
    let cx = g.eval.complex();
    for i in 0..cx.facet_count() as u32 {
        let Some(r) = cx.facet_radius(i) else { continue };
        let h_face = 2.0 * (2.0 * r * max_dev).sqrt();
        let b = cx.facet_bbox(i);
        for v in 0..m.n_points() {
            if !on_facet[v] {
                continue;
            }
            let p = [m.points[3 * v], m.points[3 * v + 1], m.points[3 * v + 2]];
            if (0..3).any(|k| p[k] < b.lo[k] - tol || p[k] > b.hi[k] + tol) {
                continue;
            }
            let q = cx.facet_project(i, p);
            if ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2)).sqrt() <= tol {
                vh[v] = vh[v].min(h_face);
            }
        }
    }
    Ok(vh.iter().map(|&h| h * a.scale).collect())
}

/// The isotropic metric the model itself asks for at every vertex (see [`vertex_sizes3`]).
pub fn a_priori_metrics3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<3>>, String> {
    Ok(vertex_sizes3(g, opts, a, false)?.iter().map(|&h| isotropic::<3>(h.max(1e-300))).collect())
}

/// The sizes the geometry insists on whatever the solution does, as metrics: the floor of an estimator metric.
pub fn geometry_floor3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<3>>, String> {
    Ok(vertex_sizes3(g, opts, a, true)?.iter().map(|&h| isotropic::<3>(h.max(1e-300))).collect())
}

/// Keeps the boundary facets that belong to one cell only (the adaptation takes labelled facets on the boundary, not between cells).
fn drop_interior_facets(m: &mut WorkMesh<3>) {
    let mut count: HashMap<[u32; 3], u32> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            let mut k = [f[0], f[1], f[2]];
            k.sort_unstable();
            *count.entry(k).or_default() += 1;
        }
    }
    let (mut facets, mut labels) = (vec![], vec![]);
    for f in 0..m.facet_label.len() {
        let mut k = [m.facets[3 * f], m.facets[3 * f + 1], m.facets[3 * f + 2]];
        k.sort_unstable();
        if count.get(&k) == Some(&1) {
            facets.extend_from_slice(&m.facets[3 * f..3 * f + 3]);
            labels.push(m.facet_label[f]);
        }
    }
    m.facets = facets;
    m.facet_label = labels;
}
