//! The evaluated region complex of a 2D model: regions (planar faces, possibly with holes), the edges that separate them
//! (each used by exactly two region sides, one of which may be the void), vertices, and the closed boundary loops of every
//! region. Areas, lengths and the self-checks are analytic (Green's theorem on the exact curves).
//!
//! Orientation convention: an edge runs along the canonical direction of its curve; `left` is the region on the left of that
//! direction, `right` the region on the right, and the edge normal points from `left` to `right` (the first listed region to
//! the second). `None` is the void (outside, or a removed hole).

use crate::curve::Curve;
use crate::math::*;
use crate::model::{FeatureId, RegionKind};

/// One origin of an edge: the feature (and pattern instance) whose boundary piece `role` produced it.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct Prov {
    pub feature: FeatureId,
    pub instance: String,
    pub role: String,
}

impl Prov {
    pub fn name(&self) -> String {
        format!("{}.{}", self.instance, self.role)
    }
}

#[derive(Clone, Debug)]
pub struct Region {
    pub name: String,
    pub kind: RegionKind,
    pub feature: FeatureId,
    pub pinned_id: Option<i32>,
    /// index of the region in the scene (before empty generated regions were dropped)
    pub scene_index: usize,
    /// area
    pub measure: f64,
    pub bbox: Aabb,
    /// loops bounding this region (outer boundaries have positive area, holes negative)
    pub loops: Vec<u32>,
    pub parent: Option<u32>,
}

#[derive(Clone, Debug)]
pub struct Edge {
    pub name: String,
    pub curve: u32,
    pub t0: f64,
    pub t1: f64,
    pub v0: u32,
    pub v1: u32,
    pub left: Option<u32>,
    pub right: Option<u32>,
    pub prov: Vec<Prov>,
    /// primitives (scene indices) whose boundary pieces make up this edge
    pub owner_prims: Vec<u32>,
    pub length: f64,
    pub bbox: Aabb,
}

/// A 3D face fragment of a prismatic model.
#[derive(Clone, Debug)]
pub enum FaceGeom {
    /// a 2D overlay edge swept over `[z0, z1]`
    Lateral { edge: u32, z0: f64, z1: f64 },
    /// a cell of the overlay at height `z`
    Cap { cell: u32, z: f64 },
    /// a closed curved surface (sphere, torus) without edges, from a `Body`
    Quadric(crate::prim::Body),
    /// a cell of a plane of a polyhedral model: `plane` indexes `Complex::planar`, `cell` a region of that plane's 2D complex
    Planar { plane: u32, cell: u32 },
}

#[derive(Clone, Debug)]
pub struct Face {
    pub name: String,
    pub geom: FaceGeom,
    /// the normal points from `first` to `second` (lateral: left to right of the swept edge; caps: below to above)
    pub first: Option<u32>,
    pub second: Option<u32>,
    pub prov: Vec<Prov>,
    pub area: f64,
    pub bbox: Aabb,
}

/// A plane of a polyhedral model with the 2D complex of the cells on it (coordinates `(u, v)` from the origin `o`).
#[derive(Clone, Debug)]
pub struct PlaneData {
    pub n: V3,
    pub d: f64,
    pub o: V3,
    pub u: V3,
    pub v: V3,
    pub cx: Complex,
}

/// Where a 3D edge lies: a vertical run above an overlay vertex, or an overlay edge at a level.
#[derive(Clone, Debug)]
pub enum EdgeGeom {
    Vertical {
        vertex: u32,
        z0: f64,
        z1: f64,
    },
    Horizontal {
        edge: u32,
        z: f64,
    },
    /// a full circle in a plane perpendicular to the sweep axis (the rim of a cone), internal frame
    Circle {
        c: V3,
        r: f64,
    },
    /// a straight segment between two points (polyhedral models)
    Segment {
        a: V3,
        b: V3,
    },
}

/// A 3D edge: shared by two or more face fragments (three or more at a junction of regions).
#[derive(Clone, Debug)]
pub struct Edge3 {
    pub name: String,
    pub geom: EdgeGeom,
    pub v0: u32,
    pub v1: u32,
    /// the faces that meet along the edge, ascending
    pub faces: Vec<u32>,
    pub length: f64,
    pub bbox: Aabb,
}

impl Face {
    /// Clamps a fraction to `[0, 1]`.
    pub fn u_clamp(&self, u: f64) -> f64 {
        u.clamp(0.0, 1.0)
    }

    /// True when this face fragment comes from a 2D `Op::Sheet` feature (zero-thickness sheet or thick-sheet cap/sidewall).
    pub fn is_sheet(&self) -> bool {
        self.prov.iter().any(|p| p.role == "sheet")
    }

    /// True when this face fragment bounds the interior cavity of a finite-thickness `Op::Sheet` conductor (`[z_s, z_s + t_m]`).
    pub fn is_conductor_cavity(&self) -> bool {
        self.is_sheet() && (self.first.is_none() != self.second.is_none())
    }
}

#[derive(Clone, Debug)]
pub struct Vertex {
    pub p: V3,
    pub name: String,
    pub edges: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct Loop {
    pub region: u32,
    /// `(edge, forward)`; the region is on the left of the travel direction
    pub edges: Vec<(u32, bool)>,
    /// signed area enclosed (negative for holes)
    pub area: f64,
    /// for a hole: index of the enclosing outer loop
    pub outer: Option<u32>,
}

#[derive(Clone, Debug, Default)]
pub struct Complex {
    pub dim: u8,
    pub tol: f64,
    /// point about which area integrals are taken (centre of the model)
    pub origin: V2,
    pub curves: Vec<Curve>,
    pub regions: Vec<Region>,
    pub edges: Vec<Edge>,
    pub vertices: Vec<Vertex>,
    pub loops: Vec<Loop>,
    /// scene region index -> final region id (`None` for holes and for dropped empty generated regions)
    pub scene_to_region: Vec<Option<u32>>,
    /// 3D only: named face fragments, and the layer stack they were derived from
    pub faces: Vec<Face>,
    /// 3D only: crease edges between face fragments (see `topo3.rs`); `vertices` holds their end points
    pub edges3: Vec<Edge3>,
    /// polyhedral models: the planes that carry faces, each with the 2D complex of its cells
    pub planar: Vec<PlaneData>,
    /// The sweep axis (0, 1, 2) of a 3D stack model, or the out-of-plane axis (y) of a 2D model. Overlay, layers and body data are
    /// stored in the sweep frame of `scene::axis_perm` (profile plane first, sweep axis last); everything the API returns is in world coordinates
    /// (`wv` / `iv` convert). A 2D model's world frame is the wafer cross-section: `(x, 0, z)`, with the curves themselves
    /// parametrised in `(x, z)`.
    pub axis: usize,
    pub stack: Option<Box<crate::stack::LayerStack>>,
}

pub fn curve_kind(c: &Curve) -> &'static str {
    match c {
        Curve::Line { .. } => "line",
        Curve::Circle { .. } => "circle",
        Curve::Ellipse { .. } => "ellipse",
    }
}

impl Complex {
    pub fn curve_of(&self, e: u32) -> (&Curve, f64, f64) {
        let ed = &self.edges[e as usize];
        (&self.curves[ed.curve as usize], ed.t0, ed.t1)
    }

    pub fn edge_point(&self, e: u32, t: f64) -> V3 {
        to3(self.curve_of(e).0.eval(t))
    }

    /// Exact closest point of edge `e` (the arc between its end points, not the whole curve) to `p`, with its parameter
    /// (unwrapped into `[t0, t1]`). Stable at periodic seams: parameters are taken relative to `t0`.
    pub fn edge_closest(&self, e: u32, p: V2) -> (f64, V2) {
        let (c, t0, t1) = self.curve_of(e);
        let (t, q) = c.closest(p);
        if c.in_range(t, t0, t1, 1e-12) {
            let u = c.param_in(t, t0).min(t1);
            return (u, if u == t { q } else { c.eval(u) });
        }
        let (a, b) = (c.eval(t0), c.eval(t1));
        if dist2(p, a) <= dist2(p, b) {
            (t0, a)
        } else {
            (t1, b)
        }
    }

    pub fn project_edge(&self, e: u32, p: V2) -> V2 {
        self.edge_closest(e, p).1
    }

    /// Unit normal at parameter `t` pointing from the left region to the right region.
    pub fn edge_normal(&self, e: u32, t: f64) -> V2 {
        let tg = self.curve_of(e).0.tangent(t);
        [tg[1], -tg[0]]
    }

    /// Outward normal for boundary edges (towards the void), the left-to-right normal for interfaces.
    pub fn outward_normal(&self, e: u32, t: f64) -> V2 {
        let n = self.edge_normal(e, t);
        match self.edges[e as usize].left {
            None => mul2(n, -1.0),
            Some(_) => n,
        }
    }

    /// Internal-frame vector to world coordinates.
    pub fn wv(&self, v: V3) -> V3 {
        let p = crate::scene::axis_perm(self.axis);
        let mut o = [0.0; 3];
        for i in 0..3 {
            o[p[i]] = v[i];
        }
        o
    }
    /// World vector to the internal frame.
    pub fn iv(&self, v: V3) -> V3 {
        let p = crate::scene::axis_perm(self.axis);
        [v[p[0]], v[p[1]], v[p[2]]]
    }
    /// True when the complex contains at least one 2D sheet cap face.
    pub fn has_sheets(&self) -> bool {
        self.faces.iter().any(|f| f.is_sheet())
    }
    /// Sorted, deduplicated names of all 2D sheet features present in the complex.
    pub fn sheet_names(&self) -> Vec<String> {
        let mut v: Vec<String> = self.faces.iter().filter(|f| f.is_sheet()).map(|f| f.name.clone()).collect();
        v.sort();
        v.dedup();
        v
    }
    /// 3D edge indices (`edges3`) forming the exterior 1D perimeter of `<sheet>`'s union of cap cells
    /// at its horizontal cap plane(s) (excluding internal fragment seams shared by two cap cells of the same sheet).
    pub fn sheet_perimeter_edges(&self, sheet: &str) -> Vec<u32> {
        if self.dim != 3 {
            return vec![];
        }
        (0..self.edges3.len() as u32)
            .filter(|&ei| {
                let e = &self.edges3[ei as usize];
                matches!(e.geom, EdgeGeom::Horizontal { .. })
                    && e.faces
                        .iter()
                        .filter(|&&fi| {
                            let f = &self.faces[fi as usize];
                            f.is_sheet() && matches!(f.geom, FaceGeom::Cap { .. }) && f.name == sheet
                        })
                        .count()
                        == 1
            })
            .collect()
    }
    /// Physical metal film thickness `t_m` (m) of a 3D `Op::Sheet` feature by name, if present.
    pub fn sheet_thickness(&self, sheet: &str) -> Option<f64> {
        self.stack.as_ref()?.sheet_thickness_m.get(sheet).copied()
    }
    /// Substrate over-etch trench depth `t_e` (m) of a 3D `Op::Sheet` feature by name, if present.
    pub fn sheet_etch_depth(&self, sheet: &str) -> Option<f64> {
        self.stack.as_ref()?.sheet_etch_depth_m.get(sheet).copied()
    }
    /// Map of sheet name -> metal film thickness `t_m` (m) on a 3D stack complex.
    pub fn sheet_thickness_m(&self) -> std::collections::BTreeMap<String, f64> {
        self.stack.as_ref().map_or_else(Default::default, |s| s.sheet_thickness_m.clone())
    }
    /// Map of sheet name -> substrate over-etch trench depth `t_e` (m) on a 3D stack complex.
    pub fn sheet_etch_depth_m(&self) -> std::collections::BTreeMap<String, f64> {
        self.stack.as_ref().map_or_else(Default::default, |s| s.sheet_etch_depth_m.clone())
    }
    /// True when 3D facet `i` bounds an interior finite-thickness `Op::Sheet` conductor cavity (`[z_s, z_s + t_m]`).
    pub fn facet_is_conductor_cavity(&self, i: u32) -> bool {
        self.dim == 3 && self.faces.get(i as usize).is_some_and(Face::is_conductor_cavity)
    }
    /// Converts the stored boxes and points of a freshly built stack complex to world coordinates and records the axis.
    pub(crate) fn set_world_frame(&mut self, axis: usize) {
        self.axis = axis;
        let conv = |c: &Complex, b: Aabb| {
            if b.is_empty() {
                b
            } else {
                let (lo, hi) = (c.wv(b.lo), c.wv(b.hi));
                Aabb {
                    lo: [lo[0].min(hi[0]), lo[1].min(hi[1]), lo[2].min(hi[2])],
                    hi: [lo[0].max(hi[0]), lo[1].max(hi[1]), lo[2].max(hi[2])],
                }
            }
        };
        let me = self.clone();
        self.regions.iter_mut().for_each(|r| r.bbox = conv(&me, r.bbox));
        self.faces.iter_mut().for_each(|f| f.bbox = conv(&me, f.bbox));
        self.edges3.iter_mut().for_each(|e| e.bbox = conv(&me, e.bbox));
        self.edges.iter_mut().for_each(|e| e.bbox = conv(&me, e.bbox));
        self.vertices.iter_mut().for_each(|v| v.p = me.wv(v.p));
    }

    // ---- uniform "facet" access (the codimension-1 entities: edges of a 2D complex, faces of a 3D one), used by selectors ----

    pub fn facet_count(&self) -> usize {
        if self.dim == 3 {
            self.faces.len()
        } else {
            self.edges.len()
        }
    }
    pub fn facet_name(&self, i: u32) -> &str {
        if self.dim == 3 {
            &self.faces[i as usize].name
        } else {
            &self.edges[i as usize].name
        }
    }
    pub fn facet_prov(&self, i: u32) -> &[Prov] {
        if self.dim == 3 {
            &self.faces[i as usize].prov
        } else {
            &self.edges[i as usize].prov
        }
    }
    /// `(first, second)`: the facet normal points from the first region to the second.
    pub fn facet_sides(&self, i: u32) -> (Option<u32>, Option<u32>) {
        if self.dim == 3 {
            let f = &self.faces[i as usize];
            (f.first, f.second)
        } else {
            let e = &self.edges[i as usize];
            (e.left, e.right)
        }
    }
    pub fn facet_bbox(&self, i: u32) -> Aabb {
        if self.dim == 3 {
            self.faces[i as usize].bbox
        } else {
            self.edges[i as usize].bbox
        }
    }
    /// Length (2D) or area (3D).
    pub fn facet_measure(&self, i: u32) -> f64 {
        if self.dim == 3 {
            self.faces[i as usize].area
        } else {
            self.edges[i as usize].length
        }
    }
    fn overlay(&self) -> &Complex {
        &self.stack.as_ref().expect("a 3D complex has a layer stack").overlay
    }
    /// The curve a facet lies on (2D edge, or the swept overlay edge of a lateral face).
    fn facet_curve(&self, i: u32) -> Option<&Curve> {
        if self.dim == 3 {
            match self.faces[i as usize].geom {
                FaceGeom::Lateral { edge, .. } => {
                    Some(&self.overlay().curves[self.overlay().edges[edge as usize].curve as usize])
                }
                FaceGeom::Cap { .. } | FaceGeom::Quadric(_) | FaceGeom::Planar { .. } => None,
            }
        } else {
            Some(&self.curves[self.edges[i as usize].curve as usize])
        }
    }
    /// `plane`, `cylinder` (3D) or `line`, `circle`, `ellipse` (2D).
    pub fn facet_kind(&self, i: u32) -> &'static str {
        if let (3, FaceGeom::Quadric(b)) =
            (self.dim, &self.faces.get(i as usize).map_or(FaceGeom::Cap { cell: 0, z: 0.0 }, |f| f.geom.clone()))
        {
            return b.kind();
        }
        match (self.dim, self.facet_curve(i)) {
            (3, Some(Curve::Line { .. })) | (3, None) => "plane",
            (3, Some(Curve::Circle { .. })) => "cylinder",
            (3, Some(_)) => "elliptic cylinder",
            (_, Some(c)) => curve_kind(c),
            (_, None) => "plane",
        }
    }
    pub fn facet_radius(&self, i: u32) -> Option<f64> {
        if self.dim == 3 {
            if let FaceGeom::Quadric(crate::prim::Body::Sphere { r, .. }) = self.faces[i as usize].geom {
                return Some(r);
            }
        }
        match self.facet_curve(i) {
            Some(Curve::Circle { r, .. }) => Some(*r),
            _ => None,
        }
    }
    /// Outward normals (towards the void; for interfaces from first to second) sampled at the ends and the middle.
    pub fn facet_normals(&self, i: u32) -> Vec<V3> {
        self.facet_normals_internal(i).into_iter().map(|v| self.wv(v)).collect()
    }
    fn facet_normals_internal(&self, i: u32) -> Vec<V3> {
        if self.dim == 3 {
            let f = &self.faces[i as usize];
            let sign = if f.first.is_none() { -1.0 } else { 1.0 };
            return match f.geom {
                FaceGeom::Lateral { edge, .. } => {
                    let ov = self.overlay();
                    let e = &ov.edges[edge as usize];
                    [e.t0, 0.5 * (e.t0 + e.t1), e.t1]
                        .iter()
                        .map(|&t| to3(mul2(ov.edge_normal(edge, t), sign)))
                        .collect()
                }
                FaceGeom::Cap { .. } => vec![[0.0, 0.0, sign]],
                FaceGeom::Planar { plane, .. } => vec![mul3(self.planar[plane as usize].n, sign)],
                FaceGeom::Quadric(b) => [(0.25, 0.25), (0.5, 0.5), (0.75, 0.75)]
                    .iter()
                    .map(|&(u, v)| mul3(b.normal(b.point_at(u, v)), sign))
                    .collect(),
            };
        }
        let e = &self.edges[i as usize];
        [e.t0, 0.5 * (e.t0 + e.t1), e.t1].iter().map(|&t| to3(self.outward_normal(i, t))).collect()
    }
    pub fn facet_vertices(&self, i: u32) -> Vec<u32> {
        if self.dim == 3 {
            let mut v: Vec<u32> =
                self.edges3.iter().filter(|e| e.faces.contains(&i)).flat_map(|e| [e.v0, e.v1]).collect();
            v.sort_unstable();
            v.dedup();
            return v;
        }
        let e = &self.edges[i as usize];
        if e.v0 == e.v1 {
            vec![e.v0]
        } else {
            vec![e.v0, e.v1]
        }
    }
    /// The point of the facet closest to `p`.
    pub fn facet_project(&self, i: u32, p: V3) -> V3 {
        self.wv(self.facet_project_internal(i, self.iv(p)))
    }
    fn facet_project_internal(&self, i: u32, p: V3) -> V3 {
        if self.dim == 2 {
            return to3(self.project_edge(i, [p[0], p[1]]));
        }
        if let FaceGeom::Planar { plane, cell } = self.faces[i as usize].geom {
            return self.planar_project(plane, cell, p);
        }
        let ov = self.overlay();
        match self.faces[i as usize].geom {
            FaceGeom::Planar { .. } => unreachable!(),
            FaceGeom::Quadric(b) => b.project(p),
            FaceGeom::Lateral { edge, z0, z1 } => {
                let q = ov.project_edge(edge, [p[0], p[1]]);
                [q[0], q[1], p[2].clamp(z0, z1)]
            }
            FaceGeom::Cap { cell, z } => {
                let xy = [p[0], p[1]];
                if ov.locate(xy) == Some(cell) {
                    return [p[0], p[1], z];
                }
                let mut best: Option<(f64, V2)> = None;
                for &l in &ov.regions[cell as usize].loops {
                    for &(e, _) in &ov.loops[l as usize].edges {
                        let q = ov.project_edge(e, xy);
                        let d = dist2(q, xy);
                        if best.is_none_or(|b| d < b.0) {
                            best = Some((d, q));
                        }
                    }
                }
                best.map_or([p[0], p[1], z], |b| [b.1[0], b.1[1], z])
            }
        }
    }
    /// Representative points of the facet (corners, edge middles, an interior point), used for congruence checks.
    pub fn facet_samples(&self, i: u32) -> Vec<V3> {
        self.facet_samples_internal(i).into_iter().map(|p| self.wv(p)).collect()
    }
    fn facet_samples_internal(&self, i: u32) -> Vec<V3> {
        if self.dim == 2 {
            let (c, t0, t1) = self.curve_of(i);
            return (0..=6).map(|k| to3(c.eval(t0 + (t1 - t0) * k as f64 / 6.0))).collect();
        }
        if let FaceGeom::Planar { plane, cell } = self.faces[i as usize].geom {
            let pd = &self.planar[plane as usize];
            let mut v: Vec<V3> = vec![];
            for &l in &pd.cx.regions[cell as usize].loops {
                for &(e, fwd) in &pd.cx.loops[l as usize].edges {
                    let ed = &pd.cx.edges[e as usize];
                    let q = pd.cx.curves[ed.curve as usize].eval(if fwd { ed.t0 } else { ed.t1 });
                    v.push(self.plane_lift(plane, q));
                }
            }
            v.truncate(48);
            if let Some(q) = crate::stack::interior_point(&pd.cx, cell) {
                v.push(self.plane_lift(plane, q));
            }
            return v;
        }
        let ov = self.overlay();
        match self.faces[i as usize].geom {
            FaceGeom::Planar { .. } => unreachable!(),
            FaceGeom::Quadric(b) => (0..4)
                .flat_map(|u| (0..=4).map(move |v| (u, v)))
                .map(|(u, v)| b.point_at(u as f64 / 4.0, v as f64 / 4.0))
                .collect(),
            FaceGeom::Lateral { edge, z0, z1 } => {
                let (c, t0, t1) = ov.curve_of(edge);
                let mut v = vec![];
                for k in 0..=4 {
                    let q = c.eval(t0 + (t1 - t0) * k as f64 / 4.0);
                    for zk in 0..=2 {
                        v.push([q[0], q[1], z0 + (z1 - z0) * zk as f64 / 2.0]);
                    }
                }
                v
            }
            FaceGeom::Cap { cell, z } => {
                let mut v: Vec<V3> = vec![];
                for &l in &ov.regions[cell as usize].loops {
                    for &(e, fwd) in &ov.loops[l as usize].edges {
                        let ed = &ov.edges[e as usize];
                        let c = &ov.curves[ed.curve as usize];
                        for t in [ed.t0, 0.5 * (ed.t0 + ed.t1)] {
                            let q = c.eval(if fwd { t } else { ed.t0 + ed.t1 - t });
                            v.push([q[0], q[1], z]);
                        }
                    }
                }
                v.truncate(48);
                if let Some(q) = crate::stack::interior_point(ov, cell) {
                    v.push([q[0], q[1], z]);
                }
                v
            }
        }
    }

    pub fn plane_lift(&self, plane: u32, q: V2) -> V3 {
        let pd = &self.planar[plane as usize];
        add3(pd.o, add3(mul3(pd.u, q[0]), mul3(pd.v, q[1])))
    }

    /// Closest point of a planar cell to `p`.
    fn planar_project(&self, plane: u32, cell: u32, p: V3) -> V3 {
        let pd = &self.planar[plane as usize];
        let d = sub3(p, pd.o);
        let xy = [dot3(d, pd.u), dot3(d, pd.v)];
        if pd.cx.locate(xy) == Some(cell) {
            return self.plane_lift(plane, xy);
        }
        let mut best: Option<(f64, V2)> = None;
        for &l in &pd.cx.regions[cell as usize].loops {
            for &(e, _) in &pd.cx.loops[l as usize].edges {
                let q = pd.cx.project_edge(e, xy);
                let dd = dist2(q, xy);
                if best.is_none_or(|b| dd < b.0) {
                    best = Some((dd, q));
                }
            }
        }
        best.map_or(p, |b| self.plane_lift(plane, b.1))
    }

    /// 3D self-checks: every face separates two different sides, volumes are positive, and every region is closed
    /// (the sum of its outward area vectors vanishes: the lateral faces and the caps balance in x, y and z).
    fn check3(&self) -> Result<(), String> {
        let ov = self.stack.as_ref().map(|s| &s.overlay);
        let scale = self.regions.iter().map(|r| r.bbox.diagonal()).fold(0.0, f64::max).max(1e-300);
        let mut flux = vec![[0.0f64; 3]; self.regions.len()];
        for f in &self.faces {
            if f.first == f.second && !f.is_sheet() {
                return Err(format!("face '{}' has the same region on both sides", f.name));
            }
            let v = match f.geom {
                FaceGeom::Planar { plane, .. } => mul3(self.planar[plane as usize].n, f.area),
                FaceGeom::Lateral { edge, z0, z1 } => {
                    let ov = ov.expect("a layer stack");
                    let e = &ov.edges[edge as usize];
                    let c = &ov.curves[e.curve as usize];
                    let (a, b) = (c.eval(e.t0), c.eval(e.t1));
                    [(z1 - z0) * (b[1] - a[1]), -(z1 - z0) * (b[0] - a[0]), 0.0]
                }
                FaceGeom::Cap { .. } => [0.0, 0.0, f.area],
                FaceGeom::Quadric(_) => [0.0; 3],
            };
            for (r, s) in [(f.first, 1.0), (f.second, -1.0)] {
                if let Some(r) = r {
                    for k in 0..3 {
                        flux[r as usize][k] += s * v[k];
                    }
                }
            }
        }
        for (r, fl) in self.regions.iter().zip(&flux) {
            if r.measure <= 0.0 {
                continue;
            }
            // vertices are snapped within `tol`, so closure holds to a few tolerances times the model size
            if fl.iter().any(|x| x.abs() > 64.0 * self.tol * scale + 1e-12 * scale * scale) {
                return Err(format!("region '{}' is not closed: its outward area vectors sum to {fl:?}", r.name));
            }
        }
        self.check_topology3()
    }

    /// Region containing `p` by testing the boundary loops (chord-flattened); agrees with `Evaluated::classify` away from the boundary.
    pub fn locate(&self, p: V2) -> Option<u32> {
        (0..self.regions.len() as u32).find(|&r| {
            self.regions[r as usize]
                .loops
                .iter()
                .filter(|&&l| crate::arrange::point_in_loop(&self.loops[l as usize], &self.edges, &self.curves, p))
                .count()
                % 2
                == 1
        })
    }

    pub fn region_area(&self, r: u32) -> f64 {
        self.regions[r as usize].measure
    }

    /// Total length of the interface between two regions.
    pub fn interface_length(&self, a: u32, b: u32) -> f64 {
        (0..self.facet_count() as u32)
            .filter(|&i| {
                let (l, r) = self.facet_sides(i);
                (l == Some(a) && r == Some(b)) || (l == Some(b) && r == Some(a))
            })
            .map(|i| self.facet_measure(i))
            .sum()
    }

    /// Twice-signed-area integral of one edge in its canonical direction.
    pub fn edge_green(&self, e: u32) -> f64 {
        let (c, t0, t1) = self.curve_of(e);
        c.green(t0, t1, self.origin)
    }

    /// Always-on self-checks (cheap): every edge separates two different sides, loops are closed and consistently oriented,
    /// region areas are positive, and the regions tile their union (the sum of region areas equals the area enclosed by the
    /// edges that touch the void).
    pub fn check(&self) -> Result<(), String> {
        if self.dim == 3 {
            return self.check3();
        }
        let scale = self.regions.iter().map(|r| r.bbox.diagonal()).fold(0.0, f64::max).max(1e-300);
        for (i, e) in self.edges.iter().enumerate() {
            if e.left == e.right {
                return Err(format!("edge {i} ('{}') has the same region on both sides", e.name));
            }
            if e.length.is_nan() || e.length <= 0.0 {
                return Err(format!("edge {i} ('{}') has no length", e.name));
            }
        }
        for (i, l) in self.loops.iter().enumerate() {
            for (k, &(e, fwd)) in l.edges.iter().enumerate() {
                let (n_e, n_f) = l.edges[(k + 1) % l.edges.len()];
                let ed = &self.edges[e as usize];
                let end = if fwd { ed.v1 } else { ed.v0 };
                let nd = &self.edges[n_e as usize];
                let start = if n_f { nd.v0 } else { nd.v1 };
                if end != start {
                    return Err(format!(
                        "loop {i} of region '{}' is not closed at edge '{}'",
                        self.regions[l.region as usize].name, ed.name
                    ));
                }
            }
        }
        let mut total = 0.0;
        for r in &self.regions {
            let a: f64 = r.loops.iter().map(|&l| self.loops[l as usize].area).sum();
            if (a - r.measure).abs() > 1e-12 * r.measure.abs().max(1e-300) + 1e-14 * scale * scale {
                return Err(format!("region '{}' area {} disagrees with its loops {a}", r.name, r.measure));
            }
            if r.measure < 0.0 {
                return Err(format!("region '{}' has negative area {}", r.name, r.measure));
            }
            total += r.measure;
        }
        // area of the union from the edges that border the void
        let mut outer = 0.0;
        for (i, e) in self.edges.iter().enumerate() {
            let g = 0.5 * self.edge_green(i as u32);
            if e.left.is_some() && e.right.is_none() {
                outer += g;
            } else if e.left.is_none() && e.right.is_some() {
                outer -= g;
            }
        }
        if (outer - total).abs() > 1e-10 * total.abs().max(1e-300) + 1e-13 * scale * scale {
            return Err(format!("region areas sum to {total} but the outer boundary encloses {outer}"));
        }
        Ok(())
    }
}
