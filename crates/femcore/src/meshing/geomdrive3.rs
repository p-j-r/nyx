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
use crate::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use crate::meshing::mesher3d::{mesh_3d, Domain3, Report3d, Spec3d};
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
}

impl<'a> GeomDomain3<'a> {
    pub fn new(ev: &'a Evaluated, labels: &'a LabelMap, h: f64) -> Result<Self, String> {
        if ev.list.dim != 3 {
            return Err("the 3D mesher needs a 3D model".into());
        }
        let nr = ev.regions().len();
        let mut region_h = vec![h; nr];
        for p in ev.pml_regions() {
            region_h[p.region as usize] = region_h[p.region as usize].min(p.thickness / p.layers);
        }
        // only the sizes of regions apply to a 3D mesh (a rule on edges has no meaning for the lattice)
        for rule in ev.size_rules() {
            if let (Some(rh), EntityKind::Regions) = (rule.h, rule.entities.kind) {
                for &r in &rule.entities.ids {
                    region_h[r as usize] = region_h[r as usize].min(rh);
                }
            }
        }
        let meshable: Vec<bool> = (0..nr).map(|r| !matches!(labels.domains[r].kind, RegionKind::Hole | RegionKind::Void)).collect();
        let cx = ev.complex();
        let nf = cx.facet_count();
        let sides = (0..nf as u32).map(|i| cx.facet_sides(i)).collect();
        let boxes = (0..nf as u32).map(|i| cx.facet_bbox(i)).collect();
        let face_label = (0..nf as u32)
            .map(|i| labels.label_of_edge(i).ok_or_else(|| format!("face '{}' has no boundary label", cx.facet_name(i))))
            .collect::<Result<Vec<i32>, String>>()?;
        Ok(GeomDomain3 { ev, labels, region_h, meshable, sides, boxes, face_label })
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
}

/// Meshes a 3D model with tetrahedra at element size `opts.h`: cells labelled by region (domain numbers), boundary facets by geom face.
pub fn mesh_geom_3d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String> {
    mesh_geom_3d_report(ev, opts).map(|x| x.0)
}

/// Like `mesh_geom_3d`, with the report of the mesher and the native mesh (all tetrahedra) the solver mesh was made from.
pub fn mesh_geom_3d_report(ev: Evaluated, opts: &GeomMeshOpts) -> Result<(MeshedGeom, Report3d, CellMesh<3>), String> {
    if ev.list.dim != 3 {
        return Err("mesh_geom_3d needs a 3D model".into());
    }
    if !(opts.h > 0.0 && opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let (cells, report) = {
        let dom = GeomDomain3::new(&ev, &labels, opts.h)?;
        mesh_3d(&dom, &Spec3d::new(opts.h))?
    };
    let mesh = SimplexMesh::from_cells(&cells)?;
    Ok((MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev) }, report, cells))
}

/// A spherical cavity of radius `r` centred at the origin, meshed at element size `h`: one domain (`cavity`) and one boundary, the surface of
/// the sphere (`cavity.surface`). The sphere sits in a host box that is a hole, so nothing but the ball is meshed. The quick problem for the
/// closed-form PEC cavity (lowest mode `k r = 2.7437`), tests and benchmarks.
pub fn sphere_cavity(r: f64, h: f64) -> Result<MeshedGeom, String> {
    let mut m = geom::model::Model::new_3d();
    let host = m.box_("host", (-2.0 * r, -2.0 * r, -2.0 * r), (2.0 * r, 2.0 * r, 2.0 * r));
    let ball = m.sphere("cavity", (0.0, 0.0, 0.0), r);
    m.overlay_full(&[("host", host, RegionKind::Hole, None), ("cavity", ball, RegionKind::Material, None)], &[]);
    mesh_geom_3d(crate::meshing::geomdrive::evaluate(&m.list)?, &GeomMeshOpts::new(h))
}
