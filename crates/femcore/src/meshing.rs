//! The one way to mesh a `geom` model: a `MeshPlan` says how fine and with which mesher, `mesh_model` returns a `MeshedGeom` with the
//! names, labels and exact geometry of the model, whichever mesher made the cells.
//!
//! * `Mesher::Simple`: chords of every geom edge, constrained Delaunay triangulation, Ruppert refinement, lattice points in finer regions
//!   (`geomdrive::mesh_geom_2d`). Fast, predictable, what the app has always used. Right for most 2D cross-sections.
//! * `Mesher::Adaptive`: the same start, then the metric adaptation of `meshers::adapt` towards a size field built from the model alone
//!   (region sizes, `Size` rules, absorbing layers, curvature of the edges), graded. Vertices that move or are created on a boundary or an
//!   interface land on the *exact* geom edge (`GeomOracle2`), not on the chord polygon, and the vertices of the geometry never move. Labels
//!   of cells and facets are carried through unchanged, so names, selectors and boundary conditions mean the same on either mesh.
//!
//! `adapt_to_metrics` is the hook for an error estimator: give it a metric per vertex of an existing meshed model (from a recovered
//! Hessian, say) and it adapts the same way. Nothing else of the interface changes.

use crate::geomdrive::{evaluate, mesh_geom_2d, size_map, GeomMeshOpts, MeshedGeom};
use crate::mesh::SimplexMesh;
pub use crate::meshing_pml::{structured_pml, PmlMode, PmlReport};
use crate::meshers::adapt::drive::{adapt2_with, AdaptOptions, AdaptReport};
use crate::meshers::adapt::field::{BackgroundMetric2, MetricField};
use crate::meshers::adapt::geom::Geometry2;
use crate::meshers::adapt::metric::{intersect, isotropic, Sym};
use crate::physics_size::{self, PhysicsSizing};
use crate::problem::ProblemSpec;
use crate::meshers::adapt::topo::{Sig, Work2};
use crate::meshers::interop::{simplex_from_adapt, simplex_to_adapt};
use geom::diff::{diff, ChangeKind};
use geom::eval::{Evaluated, Evaluator, Tier};
use geom::labels::LabelMap;
use geom::model::FeatureList;
use geom::select::EntityKind;
use std::collections::HashMap;

/// Settings of the adaptive mesher.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Adaptive {
    /// Factor on every size the model asks for: 1 keeps them, 0.5 makes the mesh about four times as fine. The start mesh is made at the
    /// model's own sizes, so a factor below one is what adaptation adds.
    pub scale: f64,
    /// Largest growth of the element size per unit distance (see `metric::grade`). Adapting to an ungraded field is unreliable; keep it.
    pub grade: f64,
    /// Stop refining beyond this many cells (the mesh stays valid, `AdaptReport::capped` says so).
    pub max_cells: Option<usize>,
}

impl Default for Adaptive {
    fn default() -> Self {
        Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(2_000_000) }
    }
}

/// Which mesher makes the cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mesher {
    Simple,
    Adaptive(Adaptive),
}

/// How a model is meshed: the sizes (as for `mesh_geom_2d`) and the mesher.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshPlan {
    pub opts: GeomMeshOpts,
    pub mesher: Mesher,
    /// how the absorbing layers of the model are meshed
    pub pml: PmlMode,
}

impl MeshPlan {
    /// The simple mesher at element size `h`.
    pub fn simple(h: f64) -> Self {
        MeshPlan { opts: GeomMeshOpts::new(h), mesher: Mesher::Simple, pml: PmlMode::Inline }
    }

    /// The adaptive mesher: start at element size `h`, adapt to `scale * h` (and to the curvature of the edges at that scale).
    pub fn adaptive(h: f64, scale: f64) -> Self {
        MeshPlan { opts: GeomMeshOpts::new(h), mesher: Mesher::Adaptive(Adaptive { scale, ..Adaptive::default() }), pml: PmlMode::Inline }
    }
}

impl From<GeomMeshOpts> for MeshPlan {
    fn from(opts: GeomMeshOpts) -> Self {
        MeshPlan { opts, mesher: Mesher::Simple, pml: PmlMode::Inline }
    }
}

/// Meshes a model.
pub fn mesh_model(list: &FeatureList, plan: &MeshPlan) -> Result<MeshedGeom, String> {
    mesh_evaluated(evaluate(list)?, plan)
}

/// Meshes an evaluated model.
pub fn mesh_evaluated(ev: Evaluated, plan: &MeshPlan) -> Result<MeshedGeom, String> {
    Ok(mesh_evaluated_report(ev, plan)?.0)
}

/// Like `mesh_evaluated`, with the report of the adaptation when there was one.
pub fn mesh_evaluated_report(ev: Evaluated, plan: &MeshPlan) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let start = mesh_geom_2d(ev, &plan.opts)?;
    let (g, report) = match plan.mesher {
        Mesher::Simple => (start, None),
        Mesher::Adaptive(a) => {
            let metrics = a_priori_metrics(&start, &plan.opts, &a)?;
            let (g, r) = adapt_to_metrics(start, &metrics, &a)?;
            (g, Some(r))
        }
    };
    match plan.pml {
        PmlMode::Inline => Ok((g, report)),
        PmlMode::Structured => Ok((structured_pml(g)?.0, report)),
    }
}

/// Meshes a model for a study: the sizes the model asks for, and in addition the sizes the waves of `spec` need
/// (`physics_size`: points per wavelength by element order). The adaptive mesher intersects the two metrics, so each is a floor on
/// the resolution. The simple mesher cannot vary the size by the wavelength, so it takes the smallest physics size as its global size
/// when that is below the plan's.
pub fn mesh_for_problem(list: &FeatureList, plan: &MeshPlan, spec: &ProblemSpec, sz: &PhysicsSizing) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let mut opts = plan.opts;
    let (g, report) = match plan.mesher {
        Mesher::Simple => {
            let probe = mesh_geom_2d(evaluate(list)?, &opts)?;
            let hmin = physics_size::vertex_sizes(spec, &probe.mesh, sz)?.into_iter().fold(f64::INFINITY, f64::min);
            if hmin < opts.h {
                opts.h = hmin;
                (mesh_geom_2d(evaluate(list)?, &opts)?, None)
            } else {
                (probe, None)
            }
        }
        Mesher::Adaptive(a) => {
            let start = mesh_geom_2d(evaluate(list)?, &opts)?;
            let geo = a_priori_metrics(&start, &opts, &a)?;
            let phys = physics_size::metrics(spec, &start.mesh, sz, opts.h)?;
            let metrics: Vec<Sym<2>> = geo
                .iter()
                .zip(&phys)
                .map(|(g, p)| {
                    // the physics size is scaled like the geometry sizes, so `Adaptive::scale` stays one knob for the whole mesh
                    let p = [[p[0][0] / (a.scale * a.scale), p[0][1] / (a.scale * a.scale)], [p[1][0] / (a.scale * a.scale), p[1][1] / (a.scale * a.scale)]];
                    intersect(g, &p)
                })
                .collect();
            let (g, r) = adapt_to_metrics(start, &metrics, &a)?;
            (g, Some(r))
        }
    };
    match plan.pml {
        PmlMode::Inline => Ok((g, report)),
        PmlMode::Structured => Ok((structured_pml(g)?.0, report)),
    }
}

/// The isotropic metric at every vertex of a meshed model that the model itself asks for: the size of the region(s) the vertex belongs
/// to, the `Size` rules of the edges it lies on, and on a curved edge the size at which a chord stays within `max_dev` of the curve. All
/// times `a.scale`.
pub fn a_priori_metrics(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String> {
    Ok(vertex_sizes(g, opts, a, false)?.iter().map(|&h| isotropic::<2>(h.max(1e-300))).collect())
}

/// The sizes the geometry insists on whatever the solution does: `Size` rules on regions and edges, the layers of an absorbing region, and
/// the curvature limit of the edges. Everywhere else the size is the global `opts.h` (a ceiling). Used as the floor of an estimator metric.
pub fn geometry_floor(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String> {
    Ok(vertex_sizes(g, opts, a, true)?.iter().map(|&h| isotropic::<2>(h.max(1e-300))).collect())
}

fn vertex_sizes(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive, geometry_only: bool) -> Result<Vec<f64>, String> {
    if !(a.scale > 0.0 && a.scale.is_finite()) {
        return Err("the adaptation scale must be positive".into());
    }
    let ev = &g.eval;
    let sizes = size_map(ev, opts.h);
    let m = &g.mesh;
    let max_dev = opts.max_dev.unwrap_or(opts.h / 50.0);
    // region index of a domain label
    let region_of: HashMap<i32, usize> = g.labels.domains.iter().map(|d| (d.id, d.region as usize)).collect();
    let mut vh = vec![f64::INFINITY; m.n_points()];
    for c in 0..m.n_cells() {
        let mut h = region_of.get(&m.cell_labels[c]).map_or(opts.h, |&r| sizes.region_h[r]);
        if geometry_only && h >= opts.h * 0.999 {
            h = opts.h;
        }
        for &v in m.cell(c) {
            vh[v as usize] = vh[v as usize].min(h);
        }
    }
    // vertices on a boundary or an interface: the only ones an edge rule or a curvature limit can concern
    let constrained = constrained_vertices(m);
    let (lo, hi) = m.bbox();
    let tol = 1e-6 * (hi[0] - lo[0]).max(hi[1] - lo[1]).max(opts.h * 1e-3);
    for (i, ed) in ev.edges().iter().enumerate() {
        let curved = ev.min_radius_of_curvature(i as u32).is_finite();
        let mut h_edge = sizes.edge_rule[i].unwrap_or(f64::INFINITY);
        if curved {
            let r = ev.min_radius_of_curvature(i as u32);
            h_edge = h_edge.min(2.0 * (2.0 * r * max_dev).sqrt());
        }
        if !h_edge.is_finite() {
            continue;
        }
        let (blo, bhi) = (ed.bbox.lo, ed.bbox.hi);
        for &v in &constrained {
            let p = [m.points[2 * v as usize], m.points[2 * v as usize + 1]];
            if p[0] < blo[0] - tol || p[0] > bhi[0] + tol || p[1] < blo[2] - tol || p[1] > bhi[2] + tol {
                continue;
            }
            let (q, _) = ev.edge_closest(i as u32, [p[0], 0.0, p[1]]);
            if (q[0] - p[0]).hypot(q[2] - p[1]) <= tol {
                vh[v as usize] = vh[v as usize].min(h_edge);
            }
        }
    }
    Ok(vh.iter().map(|&h| h * a.scale).collect())
}

/// Vertices of the mesh that lie on a boundary facet or between cells of different domains.
fn constrained_vertices(m: &SimplexMesh) -> Vec<u32> {
    let mut seen: HashMap<(u32, u32), (u32, i32)> = HashMap::new();
    let mut out = vec![false; m.n_points()];
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for k in 0..3 {
            let (a, b) = (v[k].min(v[(k + 1) % 3]), v[k].max(v[(k + 1) % 3]));
            match seen.get_mut(&(a, b)) {
                None => {
                    seen.insert((a, b), (1, m.cell_labels[c]));
                }
                Some(e) => {
                    if e.1 != m.cell_labels[c] {
                        e.0 = 1;
                        e.1 = i32::MIN;
                    } else {
                        e.0 += 1;
                    }
                }
            }
        }
    }
    for ((a, b), (n, l)) in seen {
        if n == 1 || l == i32::MIN {
            out[a as usize] = true;
            out[b as usize] = true;
        }
    }
    (0..m.n_points() as u32).filter(|&v| out[v as usize]).collect()
}

/// Adapts a meshed model to a metric given at every vertex of its mesh. The metrics are graded and interpolated on the old mesh; the
/// vertices of the geometry stay, boundary and interface vertices are kept on the exact edges, and the labels of cells and facets are
/// carried through.
pub fn adapt_to_metrics(g: MeshedGeom, metrics: &[Sym<2>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    if g.mesh.dim != 2 {
        return Err("adaptation works on 2D meshes".into());
    }
    if metrics.len() != g.mesh.n_points() {
        return Err(format!("{} metrics for {} vertices", metrics.len(), g.mesh.n_points()));
    }
    let old = simplex_to_adapt::<2>(&g.mesh)?;
    let field = BackgroundMetric2::new_graded(&old, metrics, a.grade);
    adapt_to_field(g, &field, a)
}

/// Adapts a meshed model to a metric field defined anywhere in the plane, for example the final metric of an earlier, slightly different
/// model (`BackgroundMetric2` on that model's mesh). Same guarantees as `adapt_to_metrics`.
pub fn adapt_to_field(g: MeshedGeom, field: &dyn MetricField<2>, a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    if g.mesh.dim != 2 {
        return Err("adaptation works on 2D meshes".into());
    }
    let old = simplex_to_adapt::<2>(&g.mesh)?;
    let mut opts = AdaptOptions { max_cells: a.max_cells, ..AdaptOptions::default() };
    opts.topo.fixed_vertices = geometry_vertices(&g);
    opts.topo.frozen_labels = periodic_labels(&g);
    let result = {
        let ev = &g.eval;
        let labels = &g.labels;
        adapt2_with(&old, field, &opts, |w| GeomOracle2::from_work(ev, labels, w))?
    };
    let mesh = simplex_from_adapt(&result.mesh);
    Ok((MeshedGeom { mesh, labels: g.labels, eval: g.eval }, result.report))
}

/// Boundary labels of the edges that take part in a periodic pair of the model. Their chains stay exactly as the start mesh has them, so
/// the two sides of a pair keep matching vertex for vertex (the start mesh cuts congruent edges identically).
fn periodic_labels(g: &MeshedGeom) -> Vec<i32> {
    let mut v: Vec<i32> = vec![];
    for p in g.eval.periodic_pairs() {
        for set in [&p.a, &p.b] {
            if set.kind == EntityKind::Edges {
                v.extend(set.ids.iter().filter_map(|&e| g.labels.label_of_edge(e)));
            }
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

/// Indices of the mesh vertices that sit on a vertex of the geometry.
fn geometry_vertices(g: &MeshedGeom) -> Vec<u32> {
    let m = &g.mesh;
    let (lo, hi) = m.bbox();
    let scale = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    let tol = 1e-7 * scale;
    let mut out = vec![];
    for gv in g.eval.vertices() {
        let q = [gv.p[0], gv.p[2]];
        let best = (0..m.n_points())
            .map(|i| (i, (m.points[2 * i] - q[0]).hypot(m.points[2 * i + 1] - q[1])))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if let Some((i, d)) = best {
            if d <= tol {
                out.push(i as u32);
            }
        }
    }
    out
}

/// The exact geometry as the oracle of the adaptation: where a vertex of a boundary or interface chain may go is the closest point on
/// the geom edge(s) of that chain, so refinement converges to the true circle or ellipse instead of the chord polygon it started from.
pub struct GeomOracle2<'a> {
    ev: &'a Evaluated,
    /// geom edges of every chain of the working mesh
    cands: Vec<Vec<u32>>,
}

impl<'a> GeomOracle2<'a> {
    pub fn from_work(ev: &'a Evaluated, labels: &LabelMap, w: &Work2) -> Self {
        let domain_of: Vec<i32> = {
            let mut v = vec![i32::MIN; ev.regions().len()];
            for d in &labels.domains {
                v[d.region as usize] = d.id;
            }
            v
        };
        let side = |r: Option<u32>| r.map_or(i32::MIN, |r| domain_of[r as usize]);
        let cands = w
            .chains
            .iter()
            .map(|c| match c.sig {
                Sig::Bnd(l) => labels.edges_of_label(l),
                Sig::Iface(a, b) => ev
                    .edges()
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| {
                        let (l, r) = (side(e.left), side(e.right));
                        (l == a && r == b) || (l == b && r == a)
                    })
                    .map(|(i, _)| i as u32)
                    .collect(),
            })
            .collect();
        GeomOracle2 { ev, cands }
    }
}

impl Geometry2 for GeomOracle2<'_> {
    fn project(&self, p: [f64; 2], chain: u32) -> [f64; 2] {
        let mut best = (f64::INFINITY, p);
        for &e in &self.cands[chain as usize] {
            let bb = &self.ev.edges()[e as usize].bbox;
            // distance from the point to the box of the edge is a lower bound of the distance to the edge
            let dx = (bb.lo[0] - p[0]).max(p[0] - bb.hi[0]).max(0.0);
            let dz = (bb.lo[2] - p[1]).max(p[1] - bb.hi[2]).max(0.0);
            if dx.hypot(dz) >= best.0 {
                continue;
            }
            let (q, _) = self.ev.edge_closest(e, [p[0], 0.0, p[1]]);
            let d = (q[0] - p[0]).hypot(q[2] - p[1]);
            if d < best.0 {
                best = (d, [q[0], q[2]]);
            }
        }
        best.1
    }
}

// ------------------------------------------------------------------ caching

/// How a cached request was answered.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Reuse {
    /// the same model and plan as an earlier request: nothing was recomputed
    Hit,
    /// meshed afresh; `change` says how the model differs from the one meshed just before (`None` for the first request), which is what
    /// a mesh morph could use to skip the remeshing
    Meshed { change: Option<ChangeKind> },
}

/// Counters of a `MeshCache`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CacheStats {
    pub requests: usize,
    pub hits: usize,
    pub eval_seconds: f64,
    pub mesh_seconds: f64,
}

/// Meshes a series of related models (a parameter sweep, a refinement study, the same model for several physics) without repeating work.
///
/// * One geom `Evaluator` lives as long as the cache, so the primitives and pair intersections of features that did not change are reused
///   between evaluations.
/// * The 2D meshed model of every recent request is kept under the model and plan that made it; asking again, or extruding the same section
///   through different layers (`mesh_extruded`), costs a copy of the arrays.
///
/// The cache is keyed by content (the serialised feature list and the plan), never by identity, so equal models made separately hit.
pub struct MeshCache {
    evaluator: Evaluator,
    entries: Vec<(String, MeshedGeom)>,
    capacity: usize,
    last_eval: Option<std::sync::Arc<Evaluated>>,
    pub stats: CacheStats,
}

impl Default for MeshCache {
    fn default() -> Self {
        MeshCache::new(8)
    }
}

impl MeshCache {
    /// A cache that keeps the `capacity` most recent distinct requests.
    pub fn new(capacity: usize) -> Self {
        MeshCache { evaluator: Evaluator::new(), entries: vec![], capacity: capacity.max(1), last_eval: None, stats: CacheStats::default() }
    }

    fn key(list: &FeatureList, plan: &MeshPlan) -> Result<String, String> {
        let l = serde_json::to_string(list).map_err(|e| e.to_string())?;
        Ok(format!("{l}\u{1}{plan:?}"))
    }

    /// The meshed model for `list` and `plan`, from the cache when it is there.
    pub fn mesh(&mut self, list: &FeatureList, plan: &MeshPlan) -> Result<(MeshedGeom, Reuse), String> {
        self.stats.requests += 1;
        let key = Self::key(list, plan)?;
        if let Some(i) = self.entries.iter().position(|e| e.0 == key) {
            // most recently used last
            let e = self.entries.remove(i);
            let g = e.1.clone();
            self.entries.push(e);
            self.stats.hits += 1;
            return Ok((g, Reuse::Hit));
        }
        let t = std::time::Instant::now();
        let ev = self.evaluator.eval(list, Tier::Complex).map_err(|e| e.to_string())?;
        self.stats.eval_seconds += t.elapsed().as_secs_f64();
        let change = self.last_eval.as_ref().map(|prev| diff(prev, &ev).kind);
        let t = std::time::Instant::now();
        let g = mesh_evaluated(ev, plan)?;
        self.stats.mesh_seconds += t.elapsed().as_secs_f64();
        self.last_eval = Some(g.eval.clone());
        self.entries.push((key, g.clone()));
        if self.entries.len() > self.capacity {
            self.entries.remove(0);
        }
        Ok((g, Reuse::Meshed { change }))
    }

    /// The section of `list` swept along y through `layers`. The section is the cached part: changing only the layers re-extrudes it.
    pub fn mesh_extruded(&mut self, list: &FeatureList, plan: &MeshPlan, layers: &[crate::mesher::Layer]) -> Result<(MeshedGeom, Reuse), String> {
        let (mut g, reuse) = self.mesh(list, plan)?;
        g.mesh = crate::mesher::extrude_y(&g.mesh, layers, 0.0)?;
        Ok((g, reuse))
    }
}
