//! From a 3D `geom` model to a tetrahedral mesh: the adapter between `geom::Evaluated` and the lattice mesher (`mesher3d`), which only knows
//! the trait `Domain3`.
//!
//! * which region a point is in: `Evaluated::classify` (the painter rule over the exact primitives, no tessellation);
//! * where a new vertex on a boundary goes: `Complex::facet_project`, the closest point of the exact face (plane, cylinder wall, sphere,
//!   torus), so refinement lands on the true surface;
//! * the integers: the model's `LabelMap`. Cells carry the domain number of their region (1-based, overlay order), boundary facets (the
//!   outer faces and the faces between regions) the boundary number of the geom face they lie on (1-based, one per face unless pinned).
//!
//! Same axes as the model, nothing converted: x lateral, y the direction of propagation, z the surface normal (ADR 1). Holes and voids give
//! no cells. What geom can evaluate in 3D is a prismatic stack, curved bodies inside a prismatic host, or a polyhedral model; this adapter
//! meshes whatever it evaluates and adds no restriction of its own.

use crate::mesh::cells::CellMesh;
use crate::mesh::SimplexMesh;
use crate::meshing::curve3::curve_boundary_edges;
use crate::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use crate::meshing::mesher3d::{mesh_3d, Domain3, PeriodicAxis3, Report3d, Spec3d};
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use geom::model::RegionKind;
use geom::select::EntityKind;

/// A `geom` model as the geometry of the lattice mesher.
pub struct GeomDomain3<'a> {
    ev: &'a Evaluated,
    labels: &'a LabelMap,
    /// element size per region (the global size, lowered by `Size` rules on regions and by the layers of an absorbing piece)
    region_h: Vec<f64>,
    meshable: Vec<bool>,
    /// per face: the regions on its two sides, its box
    sides: Vec<(Option<u32>, Option<u32>)>,
    boxes: Vec<geom::math::Aabb>,
    face_label: Vec<i32>,
    /// per face: the element size next to it that a `Size` rule or the chord limit asks for (infinite: none), ADR 83
    face_h: Vec<f64>,
    /// axis-aligned periodic face pairs whose planes are aligned and enforced parent->child by the lattice mesher
    periodic_axes: Vec<PeriodicAxis3>,
}

impl<'a> GeomDomain3<'a> {
    pub fn new(ev: &'a Evaluated, labels: &'a LabelMap, h: f64) -> Result<Self, String> {
        Self::new_with(ev, labels, h, &Sizes3::default())
    }

    /// Like `new`, with the extra requests of `sizes`: a finer element size per region (the wavelength of a study, ADR 50) and regions that
    /// are left out of the mesh (the absorbing pieces when the layers are built structured, ADR 52).
    pub fn new_with(ev: &'a Evaluated, labels: &'a LabelMap, h: f64, sizes: &Sizes3) -> Result<Self, String> {
        if ev.list.dim != 3 {
            return Err("the 3D mesher needs a 3D model".into());
        }
        let nr = ev.regions().len();
        let mut region_h = vec![h; nr];
        for p in ev.pml_regions() {
            region_h[p.region as usize] = region_h[p.region as usize].min(p.thickness / p.layers);
        }
        // the sizes of regions lower the size of the region; the sizes of faces (and their chord limit) are kept per face (ADR 83). A rule on
        // the crease edges or vertices of a 3D model has no meaning for the lattice and is reported as a warning by the mesh report
        let mut face_h = vec![f64::INFINITY; ev.complex().facet_count()];
        for rule in ev.size_rules() {
            match rule.entities.kind {
                EntityKind::Regions => {
                    if let Some(rh) = rule.h {
                        for &r in &rule.entities.ids {
                            region_h[r as usize] = region_h[r as usize].min(rh);
                        }
                    }
                }
                EntityKind::Faces => {
                    for &i in &rule.entities.ids {
                        let i = i as usize;
                        if let Some(rh) = rule.h {
                            face_h[i] = face_h[i].min(rh);
                        }
                        if let (Some(md), Some(r)) = (rule.max_dev, ev.complex().facet_radius(i as u32)) {
                            face_h[i] = face_h[i].min(chord_size(r, md));
                        }
                    }
                }
                _ => {}
            }
        }
        // the wave sizes of a study ask for finer elements only; a region that is not meshed needs none
        for (r, &want) in sizes.region_h.iter().enumerate().take(nr) {
            if want.is_finite() && want > 0.0 {
                region_h[r] = region_h[r].min(want);
            }
        }
        let skipped = |r: usize| sizes.skip.get(r).copied().unwrap_or(false);
        let meshable: Vec<bool> = (0..nr).map(|r| !skipped(r) && !matches!(labels.domains[r].kind, RegionKind::Hole | RegionKind::Void)).collect();
        let cx = ev.complex();
        let nf = cx.facet_count();
        // a face between a meshed region and a skipped one is a boundary of the mesh like a face to the void
        let eff = |s: Option<u32>| s.filter(|&r| meshable[r as usize]);
        let sides = (0..nf as u32)
            .map(|i| {
                let (a, b) = cx.facet_sides(i);
                if sizes.skip.iter().any(|&s| s) {
                    (eff(a), eff(b))
                } else {
                    (a, b)
                }
            })
            .collect();
        let boxes: Vec<geom::math::Aabb> = (0..nf as u32).map(|i| cx.facet_bbox(i)).collect();
        let face_label = (0..nf as u32)
            .map(|i| labels.label_of_edge(i).ok_or_else(|| format!("face '{}' has no boundary label", cx.facet_name(i))))
            .collect::<Result<Vec<i32>, String>>()?;
        let mut periodic_axes: Vec<PeriodicAxis3> = Vec::new();
        for pp in ev.periodic_pairs() {
            if let Some(pa) = Self::extract_periodic_axis(pp, &boxes) {
                if !periodic_axes.iter().any(|x| x.axis == pa.axis) {
                    periodic_axes.push(pa);
                }
            }
        }
        periodic_axes.sort_by_key(|p| p.axis);
        Ok(GeomDomain3 { ev, labels, region_h, meshable, sides, boxes, face_label, face_h, periodic_axes })
    }

    fn extract_periodic_axis(pp: &geom::sim::PeriodicPair, boxes: &[geom::math::Aabb]) -> Option<PeriodicAxis3> {
        if pp.a.kind != EntityKind::Faces || pp.b.kind != EntityKind::Faces || pp.a.ids.is_empty() || pp.b.ids.is_empty() {
            return None;
        }
        let t = pp.translation;
        let t_len = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        if !(t_len > 0.0 && t_len.is_finite()) {
            return None;
        }
        let axis = (0..3).find(|&k| (0..3).all(|m| m == k || t[m].abs() <= 1e-9 * t_len))?;
        let union_box = |ids: &[u32]| -> geom::math::Aabb {
            ids.iter().fold(geom::math::Aabb::EMPTY, |acc, &id| acc.union(&boxes[id as usize]))
        };
        let (ba, bb) = (union_box(&pp.a.ids), union_box(&pp.b.ids));
        if (ba.hi[axis] - ba.lo[axis]).abs() > 1e-6 * t_len || (bb.hi[axis] - bb.lo[axis]).abs() > 1e-6 * t_len {
            return None;
        }
        let ca = 0.5 * (ba.lo[axis] + ba.hi[axis]);
        let cb = 0.5 * (bb.lo[axis] + bb.hi[axis]);
        let (lo, hi) = (ca.min(cb), ca.max(cb));
        if ((hi - lo) - t_len).abs() <= 1e-6 * t_len {
            Some(PeriodicAxis3 { axis, lo, hi })
        } else {
            None
        }
    }

    /// Adds the chord limit `max_dev` (metres) of the request to every curved face: the element size at which a chord stays within it
    /// (ADR 83). `None`: no limit, so a mesh without face rules is the mesh it always was.
    pub fn with_max_dev(mut self, max_dev: Option<f64>) -> Self {
        if let Some(md) = max_dev {
            for i in 0..self.face_h.len() {
                if let Some(r) = self.ev.complex().facet_radius(i as u32) {
                    self.face_h[i] = self.face_h[i].min(chord_size(r, md));
                }
            }
        }
        self
    }

    /// A face with a size request that a non-axis-aligned periodic face is within reach of: refused when not handled by axis-aligned
    /// parent->child periodic conformance.
    pub fn check_periodic_faces(&self) -> Result<(), String> {
        if !self.face_h.iter().any(|h| h.is_finite()) {
            return Ok(());
        }
        let h_glob = self.coarsest_meshed_size();
        let gap = |a: &geom::math::Aabb, b: &geom::math::Aabb| (0..3).map(|k| (a.lo[k] - b.hi[k]).max(b.lo[k] - a.hi[k]).max(0.0)).fold(0.0, f64::max);
        for pp in self.ev.periodic_pairs() {
            if Self::extract_periodic_axis(pp, &self.boxes)
                .is_some_and(|pa| self.periodic_axes.iter().any(|x| x.axis == pa.axis && (x.lo - pa.lo).abs() <= 1e-6 * (pa.hi - pa.lo) && (x.hi - pa.hi).abs() <= 1e-6 * (pa.hi - pa.lo)))
            {
                continue;
            }
            for &pf in pp.a.ids.iter().chain(pp.b.ids.iter()) {
                for (i, &hf) in self.face_h.iter().enumerate() {
                    if hf.is_finite() && gap(&self.boxes[pf as usize], &self.boxes[i]) < (h_glob - hf).max(0.0) / 0.3 + 1e-12 {
                        return Err(format!(
                            "the size requested at face '{}' reaches the periodic face '{}': the two faces of a periodic pair would be meshed differently. Keep faces with a size away from periodic faces (a distance of {:.3e} m) or use a uniform size",
                            self.ev.complex().facet_name(i as u32),
                            self.ev.complex().facet_name(pf),
                            (h_glob - hf).max(0.0) / 0.3
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// The element size of each geom region (global size lowered by `Size` rules, absorbing layers and requested sizes).
    pub fn region_sizes(&self) -> &[f64] {
        &self.region_h
    }

    /// Whether a region contributes cells.
    pub fn is_meshed(&self, r: u32) -> bool {
        self.meshable[r as usize]
    }

    /// The largest element size any meshed region asks for.
    pub fn coarsest_meshed_size(&self) -> f64 {
        (0..self.region_h.len()).filter(|&r| self.meshable[r]).map(|r| self.region_h[r]).fold(0.0, f64::max)
    }

    /// Distance from `p` to face `i`, with a cheap lower bound from its box.
    fn face_distance(&self, i: usize, p: [f64; 3], cap: f64) -> f64 {
        let b = &self.boxes[i];
        let d = (0..3).map(|k| (b.lo[k] - p[k]).max(p[k] - b.hi[k]).max(0.0).powi(2)).sum::<f64>().sqrt();
        if d >= cap {
            return d;
        }
        let q = self.ev.complex().facet_project(i as u32, p);
        ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2)).sqrt()
    }
}

impl Domain3 for GeomDomain3<'_> {
    fn region_at(&self, p: [f64; 3]) -> Option<u32> {
        self.ev.classify(p).filter(|&r| self.meshable[r as usize])
    }
    fn n_regions(&self) -> usize {
        self.meshable.len()
    }
    fn region_box(&self, r: u32) -> ([f64; 3], [f64; 3]) {
        if !self.meshable[r as usize] {
            return ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        }
        let b = self.ev.regions()[r as usize].bbox;
        (b.lo, b.hi)
    }
    fn region_size(&self, r: u32) -> Option<f64> {
        Some(self.region_h[r as usize])
    }
    fn region_label(&self, r: u32) -> i32 {
        self.labels.domains[r as usize].id
    }
    fn boundary_label(&self, p: [f64; 3], a: Option<u32>, b: Option<u32>) -> i32 {
        // among the faces between exactly these two regions, the nearest; if the model has none (a facet on the lattice's own cut), the nearest of all
        let mut best = (f64::INFINITY, None::<usize>);
        for pass in 0..2 {
            for i in 0..self.sides.len() {
                let s = self.sides[i];
                if pass == 0 && !(s == (a, b) || s == (b, a)) {
                    continue;
                }
                let d = self.face_distance(i, p, best.0);
                if d < best.0 {
                    best = (d, Some(i));
                }
            }
            if best.1.is_some() {
                break;
            }
        }
        best.1.map_or(1, |i| self.face_label[i])
    }
    fn project_edge(&self, mid: [f64; 3], pa: [f64; 3], pb: [f64; 3]) -> Option<[f64; 3]> {
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2) + (pb[2] - pa[2]).powi(2)).sqrt();
        let tol = 1e-6 * len;
        let cx = self.ev.complex();
        for i in 0..self.sides.len() {
            if self.face_distance(i, pa, tol) <= tol && self.face_distance(i, pb, tol) <= tol {
                return Some(cx.facet_project(i as u32, mid));
            }
        }
        None
    }
    fn asks_size_at(&self) -> bool {
        let h_glob = self.coarsest_meshed_size();
        self.face_h.iter().any(|&h| h.is_finite() && h < 0.999 * h_glob)
    }
    fn size_at(&self, p: [f64; 3]) -> Option<f64> {
        // the smallest size a face asks for, graded away from it with a growth of 0.3 per unit length (a ratio of 1.3 between neighbours)
        let h_glob = self.coarsest_meshed_size();
        let mut best: Option<f64> = None;
        for (i, &hf) in self.face_h.iter().enumerate() {
            if !hf.is_finite() || hf >= h_glob {
                continue;
            }
            let cap = (h_glob - hf) / 0.3;
            let d = self.face_distance(i, p, cap);
            if d < cap {
                let s = hf + 0.3 * d;
                best = Some(best.map_or(s, |b: f64| b.min(s)));
            }
        }
        best
    }
    fn edge_on_crease(&self, pa: [f64; 3], pb: [f64; 3]) -> bool {
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2) + (pb[2] - pa[2]).powi(2)).sqrt();
        let tol = 1e-6 * len;
        let cx = self.ev.complex();
        let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0];
        // the faces that both ends lie on; they are a crease when they put the midpoint in different places
        let mut first: Option<[f64; 3]> = None;
        for i in 0..self.sides.len() {
            if self.face_distance(i, pa, tol) <= tol && self.face_distance(i, pb, tol) <= tol {
                let q = cx.facet_project(i as u32, mid);
                match first {
                    None => first = Some(q),
                    Some(f) => {
                        if ((q[0] - f[0]).powi(2) + (q[1] - f[1]).powi(2) + (q[2] - f[2]).powi(2)).sqrt() > tol {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }
    fn periodic_axes(&self) -> Vec<PeriodicAxis3> {
        self.periodic_axes.clone()
    }
}

/// The element size at which a chord of a surface with radius of curvature `r` stays within `max_dev` of it (the formula of the adaptive
/// mesher, `adapt3d`).
fn chord_size(r: f64, max_dev: f64) -> f64 {
    2.0 * (2.0 * r * max_dev).sqrt()
}

/// What a caller asks of the 3D mesh beyond the model's own sizes.
#[derive(Clone, Debug, Default)]
pub struct Sizes3 {
    /// element size wanted per geom region (index = region); `f64::INFINITY`, a missing entry or zero: no request. A request only
    /// makes the region finer than the global size. When any is given the lattice itself is made at the coarsest size any meshed
    /// region needs (never above the global one) and the finer regions are refined locally, so a study whose waves are short everywhere
    /// is not meshed at a coarse size and refined afterwards.
    pub region_h: Vec<f64>,
    /// regions left out of the mesh (index = region); the lattice is aligned with the box of the regions that are meshed
    pub skip: Vec<bool>,
}

/// Meshes a 3D model with tetrahedra at element size `opts.h`: cells labelled by region (domain numbers), boundary facets by geom face.
pub fn mesh_geom_3d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String> {
    mesh_geom_3d_report(ev, opts).map(|x| x.0)
}

/// Like `mesh_geom_3d`, with the report of the mesher and the native mesh (all tetrahedra) the solver mesh was made from.
pub fn mesh_geom_3d_report(ev: Evaluated, opts: &GeomMeshOpts) -> Result<(MeshedGeom, Report3d, CellMesh<3>), String> {
    mesh_geom_3d_with(ev, opts, &Sizes3::default())
}

/// Like `mesh_geom_3d_report`, with the extra requests of `sizes` (ADR 50, 52).
pub fn mesh_geom_3d_with(ev: Evaluated, opts: &GeomMeshOpts, sizes: &Sizes3) -> Result<(MeshedGeom, Report3d, CellMesh<3>), String> {
    if ev.list.dim != 3 {
        return Err("mesh_geom_3d needs a 3D model".into());
    }
    if !(opts.h > 0.0 && opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let (cells, report) = {
        let dom = GeomDomain3::new_with(&ev, &labels, opts.h, sizes)?.with_max_dev(opts.max_dev);
        dom.check_periodic_faces()?;
        let h_lattice = if sizes.region_h.iter().any(|h| h.is_finite() && *h > 0.0) { dom.coarsest_meshed_size().min(opts.h) } else { opts.h };
        let (mut cells, mut report) = mesh_3d(&dom, &Spec3d::new(h_lattice))?;
        // curved (quadratic) geometry for the cells on a curved surface (ADR 82), after the last change of the connectivity
        if opts.curved {
            report.curve = Some(curve_boundary_edges(&mut cells, &dom));
        }
        (cells, report)
    };
    let mesh = SimplexMesh::from_cells(&cells)?;
    Ok((MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev) }, report, cells))
}

/// A spherical cavity of radius `r` centred at the origin, meshed at element size `h`: one domain (`cavity`) and one boundary, the surface of
/// the sphere (`cavity.surface`). The sphere sits in a host box that is a hole, so nothing but the ball is meshed. The quick problem for the
/// closed-form PEC cavity (lowest mode `k r = 2.7437`), tests and benchmarks.
pub fn sphere_cavity(r: f64, h: f64) -> Result<MeshedGeom, String> {
    sphere_cavity_with(r, &GeomMeshOpts::new(h))
}

/// `sphere_cavity` with the mesh options of the caller (`curved`, `max_dev`).
pub fn sphere_cavity_with(r: f64, opts: &GeomMeshOpts) -> Result<MeshedGeom, String> {
    let mut m = geom::model::Model::new_3d();
    let host = m.box_("host", (-2.0 * r, -2.0 * r, -2.0 * r), (2.0 * r, 2.0 * r, 2.0 * r));
    let ball = m.sphere("cavity", (0.0, 0.0, 0.0), r);
    m.overlay_full(&[("host", host, RegionKind::Hole, None), ("cavity", ball, RegionKind::Material, None)], &[]);
    mesh_geom_3d(crate::meshing::geomdrive::evaluate(&m.list)?, opts)
}
