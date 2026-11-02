//! From a `geom` model to a mesh and to the integers the solvers take, in the wafer frame throughout: x lateral, y the direction of
//! propagation, z the surface normal. A 2D model is the x-z cross-section and its mesh has the columns `(x, z)`; a 3D model is
//! meshed in the same axes. Nothing is converted on the way.
//!
//! The geometry is exact (lines, circles, ellipses); the mesh is made from it:
//!
//! * every `geom` edge is cut into chords once (by chord deviation and by the local element size) and both neighbouring regions use the
//!   same points, so the triangulation conforms without any repair step;
//! * the element size of a region is the global size, lowered by the `Size` rules of the model and, in an absorbing layer, by its number
//!   of layers; an edge takes the finer of its two neighbours;
//! * cells and boundary facets carry the integers of the model's `LabelMap` (domain = region in overlay order, boundary = one per edge or
//!   the pinned integer), so names, selectors and the solver see the same labels.
//!
//! Problems are written with names: `MeshedGeom::domains("sapphire*")` is the region and its absorbing pieces, `boundary("rib|air")` the
//! labels of an interface, `pml_materials` the stretching of every absorbing piece from the data of the `Pml` feature.

use crate::mesh::SimplexMesh;
use crate::mesher::{mesh_chains, Chain, Fill, MeshParams};
use crate::problem::{Pml, PmlDir, ProblemSpec};
use geom::eval::{Evaluated, Tier};
use geom::labels::LabelMap;
use geom::model::{FeatureList, RegionKind};
use geom::query::DiscretizeOpts;
use geom::select::{EntityKind, Selector};
use std::collections::BTreeMap;

/// Sizes of a mesh made from a model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeomMeshOpts {
    /// global element size
    pub h: f64,
    /// largest distance of a chord from its curve (default `h / 50`)
    pub max_dev: Option<f64>,
    pub min_angle_deg: f64,
}

impl GeomMeshOpts {
    pub fn new(h: f64) -> Self {
        GeomMeshOpts { h, max_dev: None, min_angle_deg: 25.0 }
    }
}

/// A meshed model: the mesh, the compiled names, and the evaluation they came from.
pub struct MeshedGeom {
    pub mesh: SimplexMesh,
    pub labels: LabelMap,
    pub eval: Evaluated,
}

impl MeshedGeom {
    /// The domain label of a region.
    pub fn domain(&self, region: &str) -> Option<i32> {
        self.labels.domain(region)
    }

    /// Domain labels of regions by name. A list is comma separated; `"sapphire*"` is `sapphire` and every region whose name starts with
    /// `sapphire.` (its absorbing pieces). Unknown names give nothing.
    pub fn domains(&self, pattern: &str) -> Vec<i32> {
        let mut v: Vec<i32> = vec![];
        for item in pattern.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            match item.strip_suffix('*') {
                Some(base) => v.extend(
                    self.labels
                        .domains
                        .iter()
                        .filter(|d| d.name == base || d.name.strip_prefix(base).is_some_and(|r| r.starts_with('.')))
                        .map(|d| d.id),
                ),
                None => v.extend(self.labels.domain(item)),
            }
        }
        v.sort_unstable();
        v.dedup();
        v
    }

    /// Boundary labels by name (`"rib.top"`, `"air|rib"`; a face split into pieces has several) or tag/group.
    pub fn boundary(&self, name: &str) -> Vec<i32> {
        self.labels.boundary(name)
    }

    /// Boundary labels of whatever a selector text picks (`"outer & plane(z = -3e-6)"`).
    pub fn select(&self, selector: &str) -> Result<Vec<i32>, String> {
        select_labels(&self.eval, &self.labels, selector)
    }

    /// Sets a material property on a region or family (see `domains`).
    pub fn set_material(&self, spec: &mut ProblemSpec, regions: &str, property: &str, value: f64) {
        for d in self.domains(regions) {
            spec.materials.set(d, property, value);
        }
    }

    /// The absorbing layers of the model as solver stretching, keyed by domain: `strength` is the imaginary stretch at the outer edge, the
    /// grading of the `Pml` feature is the order of the polynomial profile.
    pub fn pml_materials(&self, strength: f64) -> BTreeMap<i32, Pml> {
        let dir = |s: i8| match s {
            0 => PmlDir::Off,
            s if s < 0 => PmlDir::Minus,
            _ => PmlDir::Plus,
        };
        self.eval
            .pml_regions()
            .iter()
            .map(|p| {
                let d = self.labels.domains[p.region as usize].id;
                // a 2D mesh is the x-z plane: the layers are written along x and z (y is out of the plane)
                let (x, y, z) = if self.mesh.dim == 2 {
                    (dir(p.sides[0]), PmlDir::Off, dir(p.sides[2]))
                } else {
                    (dir(p.sides[0]), dir(p.sides[1]), dir(p.sides[2]))
                };
                (d, Pml { x, y, z, strength, order: p.grading, re: 1.0 })
            })
            .collect()
    }
}

/// Evaluates a feature list to the tier the mesher needs.
pub fn evaluate(list: &FeatureList) -> Result<Evaluated, String> {
    geom::eval::Evaluator::new().eval(list, Tier::Complex).map_err(|e| e.to_string())
}

/// Meshes a 2D model (the x-z cross-section) with triangles.
pub fn mesh_geom_2d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String> {
    if ev.list.dim != 2 {
        return Err("mesh_geom_2d needs a 2D model".into());
    }
    if !(opts.h > 0.0 && opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let max_dev = opts.max_dev.unwrap_or(opts.h / 50.0);

    // element size per region: the global size, the Size rules, and thickness / layers in an absorbing piece
    let mut region_h = vec![opts.h; ev.regions().len()];
    for p in ev.pml_regions() {
        region_h[p.region as usize] = region_h[p.region as usize].min(p.thickness / p.layers);
    }
    let mut edge_rule: Vec<Option<f64>> = vec![None; ev.edges().len()];
    for rule in ev.size_rules() {
        let Some(h) = rule.h else { continue };
        match rule.entities.kind {
            EntityKind::Regions => {
                for &r in &rule.entities.ids {
                    region_h[r as usize] = region_h[r as usize].min(h);
                }
            }
            EntityKind::Edges => {
                for &e in &rule.entities.ids {
                    edge_rule[e as usize] = Some(edge_rule[e as usize].map_or(h, |x| x.min(h)));
                }
            }
            _ => {}
        }
    }

    // chains: one per edge, with the end points of the vertices so that neighbours share them exactly
    let w2 = |p: [f64; 3]| [p[0], p[2]];
    let mut chains = Vec::with_capacity(ev.edges().len());
    for (i, ed) in ev.edges().iter().enumerate() {
        let mut h = [ed.left, ed.right].iter().flatten().map(|&r| region_h[r as usize]).fold(opts.h, f64::min);
        if let Some(x) = edge_rule[i] {
            h = h.min(x);
        }
        let mut pts: Vec<[f64; 2]> =
            ev.discretize_edge(i as u32, DiscretizeOpts { max_dev, max_len: h }).into_iter().map(|(_, p)| w2(p)).collect();
        let last = pts.len() - 1;
        pts[0] = w2(ev.vertices()[ed.v0 as usize].p);
        pts[last] = w2(ev.vertices()[ed.v1 as usize].p);
        let label = labels.label_of_edge(i as u32).ok_or_else(|| format!("edge '{}' has no boundary label", ed.name))?;
        chains.push(Chain { pts, label, part: i });
    }

    // domain of a region; holes and the void have none
    let domain_of: Vec<Option<i32>> = labels
        .domains
        .iter()
        .map(|d| if matches!(d.kind, RegionKind::Hole | RegionKind::Void) { None } else { Some(d.id) })
        .collect();
    let b = ev.bounds();
    let (lo, hi) = ([b.lo[0], b.lo[2]], [b.hi[0], b.hi[2]]);

    // finer regions get a lattice of interior points
    let fills: Vec<Fill> = (0..ev.regions().len())
        .filter(|&r| domain_of[r].is_some() && region_h[r] < opts.h * 0.999)
        .map(|r| {
            let bb = ev.regions()[r].bbox;
            let ev = &ev;
            Fill {
                lo: [bb.lo[0], bb.lo[2]],
                hi: [bb.hi[0], bb.hi[2]],
                h: region_h[r],
                owner: Box::new(move |p| ev.classify([p[0], 0.0, p[1]]) == Some(r as u32)),
            }
        })
        .collect();
    let h_fine_min = region_h.iter().copied().fold(opts.h, f64::min);
    let domain_at = |c: [f64; 2]| ev.classify([c[0], 0.0, c[1]]).and_then(|r| domain_of[r as usize]);
    let mesh = mesh_chains(
        &chains,
        &fills,
        &MeshParams { h: opts.h, h_fine_min, h_reach: opts.h, min_angle_deg: opts.min_angle_deg, lo, hi },
        &domain_at,
    )?;
    drop(fills);
    Ok(MeshedGeom { mesh, labels, eval: ev })
}

/// Labels of whatever a selector text picks in an evaluated model: boundary labels for edges, domain labels for regions.
pub fn select_labels(ev: &Evaluated, labels: &LabelMap, selector: &str) -> Result<Vec<i32>, String> {
    let sel = Selector::parse(selector).map_err(|e| e.to_string())?;
    let r = ev.resolve(&sel).map_err(|e| e.to_string())?;
    match r.set.kind {
        EntityKind::Regions => Ok(r.set.ids.iter().map(|&i| labels.domains[i as usize].id).collect()),
        EntityKind::Vertices => Err(format!("'{selector}' selects vertices, which have no label")),
        _ => {
            let mut v: Vec<i32> = r.set.ids.iter().filter_map(|&e| labels.label_of_edge(e)).collect();
            v.sort_unstable();
            v.dedup();
            Ok(v)
        }
    }
}
