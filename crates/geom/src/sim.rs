//! Simulation operations as *data*: absorbing-layer regions, mesh size rules, periodic face pairs (verified congruent), and the
//! projection oracle. geom describes which entities and how; femcore decides what a condition means and applies it.

use crate::error::{GeomError, GeomResult};
use crate::eval::Evaluated;
use crate::math::*;
use crate::model::*;
use crate::select::{EntityKind, EntitySet};

#[derive(Clone, Debug, PartialEq)]
pub struct PmlRegion {
    pub region: u32,
    /// the physical region this layer was cut from
    pub parent: Option<u32>,
    pub window: Aabb,
    /// per axis -1 / 0 / +1: the window side(s) the layer attaches to (a corner piece has two nonzero entries)
    pub sides: [i8; 3],
    pub thickness: f64,
    pub layers: f64,
    pub grading: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SizeRule {
    pub feature: FeatureId,
    pub entities: EntitySet,
    pub h: Option<f64>,
    pub max_dev: Option<f64>,
    pub min_layers: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicPair {
    pub feature: FeatureId,
    pub a: EntitySet,
    pub b: EntitySet,
    pub translation: V3,
    pub checked_congruent: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntityRef {
    Region(u32),
    Face(u32),
    Edge(u32),
    Vertex(u32),
}

#[derive(Default, Clone, Debug)]
pub struct SimData {
    pub pml: Vec<PmlRegion>,
    pub size: Vec<SizeRule>,
    pub periodic: Vec<PeriodicPair>,
}

impl Evaluated {
    pub fn pml_regions(&self) -> &[PmlRegion] {
        &self.sim.pml
    }
    pub fn size_rules(&self) -> &[SizeRule] {
        &self.sim.size
    }
    pub fn periodic_pairs(&self) -> &[PeriodicPair] {
        &self.sim.periodic
    }

    pub(crate) fn compute_sim(&mut self) -> GeomResult<()> {
        let cx = self.complex();
        let scene = self.scene();
        let mut sim = SimData::default();
        let cx_perm = crate::scene::axis_perm(scene.axis);
        for (i, r) in cx.regions.iter().enumerate() {
            if let Some(p) = &scene.regions[r.scene_index].pml {
                sim.pml.push(PmlRegion {
                    region: i as u32,
                    parent: r.parent,
                    window: to_world_box(&cx_perm, p.window),
                    sides: {
                        let mut w = [0i8; 3];
                        for i in 0..3 {
                            w[cx_perm[i]] = p.sides[i];
                        }
                        w
                    },
                    thickness: p.thickness,
                    layers: p.layers,
                    grading: p.grading,
                });
            }
        }
        let env = &scene.env;
        let mut warnings = vec![];
        let ev_opt = |e: &Option<crate::expr::Expr>, what: &str| -> GeomResult<Option<f64>> {
            e.as_ref().map(|e| env.eval(e, what)).transpose()
        };
        for f in self.list.features.iter().filter(|f| !f.suppressed) {
            match &f.op {
                Op::Size { sel, h, max_deviation, min_layers } => {
                    let r = self.resolve(sel)?;
                    warnings.extend(r.warnings.iter().map(|w| format!("feature '{}': {w}", f.name)));
                    let (h, md, ml) = (
                        ev_opt(h, &format!("size of '{}'", f.name))?,
                        ev_opt(max_deviation, &format!("max_deviation of '{}'", f.name))?,
                        ev_opt(min_layers, &format!("min_layers of '{}'", f.name))?,
                    );
                    if h.is_none() && md.is_none() && ml.is_none() {
                        return Err(GeomError(format!(
                            "Size '{}' sets none of h, max_deviation or min_layers",
                            f.name
                        )));
                    }
                    if h.is_some_and(|v| v <= 0.0) || md.is_some_and(|v| v <= 0.0) || ml.is_some_and(|v| v < 1.0) {
                        return Err(GeomError(format!(
                            "Size '{}' needs h > 0, max_deviation > 0 and min_layers >= 1",
                            f.name
                        )));
                    }
                    sim.size.push(SizeRule { feature: f.id, entities: r.set, h, max_dev: md, min_layers: ml });
                }
                Op::Periodic { sel_a, sel_b, translation } => {
                    let (ra, rb) = (self.resolve(sel_a)?, self.resolve(sel_b)?);
                    warnings
                        .extend(ra.warnings.iter().chain(&rb.warnings).map(|w| format!("feature '{}': {w}", f.name)));
                    let t = [
                        env.eval(&translation[0], "translation")?,
                        env.eval(&translation[1], "translation")?,
                        env.eval(&translation[2], "translation")?,
                    ];
                    if ra.set.kind != EntityKind::Edges && ra.set.kind != EntityKind::Faces
                        || ra.set.kind != rb.set.kind
                    {
                        return Err(GeomError(format!(
                            "Periodic '{}' needs two sets of boundary facets of the same kind",
                            f.name
                        )));
                    }
                    if !ra.set.ids.is_empty() && !rb.set.ids.is_empty() {
                        self.check_congruent(&f.name, &ra.set, &rb.set, t)?;
                    }
                    let ok = !ra.set.ids.is_empty() && !rb.set.ids.is_empty();
                    sim.periodic.push(PeriodicPair {
                        feature: f.id,
                        a: ra.set,
                        b: rb.set,
                        translation: t,
                        checked_congruent: ok,
                    });
                }
                Op::Tag { sel, .. } => {
                    let r = self.resolve(sel)?;
                    warnings.extend(r.warnings.iter().map(|w| format!("feature '{}': {w}", f.name)));
                }
                _ => {}
            }
        }
        self.warnings.extend(warnings);
        self.sim = sim;
        Ok(())
    }

    /// The facets of `b` must be the facets of `a` moved by `t`: equal total measure and every sample point of each set lies on the other.
    fn check_congruent(&self, name: &str, a: &EntitySet, b: &EntitySet, t: V3) -> GeomResult<()> {
        let cx = self.complex();
        let tol = 8.0 * cx.tol;
        let len = |s: &EntitySet| s.ids.iter().map(|&e| cx.facet_measure(e)).sum::<f64>();
        let (la, lb) = (len(a), len(b));
        if (la - lb).abs() > 1e-9 * la.max(lb) {
            return Err(GeomError(format!(
                "Periodic '{name}': the two sides are not congruent under the translation ({:.6e}, {:.6e}, {:.6e}): total sizes differ ({la:.9e} versus {lb:.9e})",
                t[0], t[1], t[2]
            )));
        }
        let on = |set: &EntitySet, p: V3| set.ids.iter().any(|&e| norm3(sub3(cx.facet_project(e, p), p)) <= tol);
        for (src, dst, sign) in [(a, b, 1.0), (b, a, -1.0)] {
            for &e in &src.ids {
                for p in cx.facet_samples(e) {
                    if !on(dst, add3(p, mul3(t, sign))) {
                        return Err(GeomError(format!(
                            "Periodic '{name}': '{}' has no counterpart under the translation ({:.6e}, {:.6e}, {:.6e}); the two sides are not congruent",
                            cx.facet_name(e), t[0], t[1], t[2]
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// Projection oracle: closest point on an exact entity (shaped to be wrapped by femcore's `Geometry2::project`).
    pub fn oracle(&self) -> Oracle<'_> {
        Oracle { ev: self }
    }
}

pub struct Oracle<'a> {
    ev: &'a Evaluated,
}

impl Oracle<'_> {
    /// Closest point of the entity to `p`; `None` for an unknown id or an entity kind the model does not have.
    pub fn project(&self, p: V3, e: EntityRef) -> Option<V3> {
        let cx = self.ev.complex();
        // world -> profile frame; the 2D results are mapped back below
        let pi = cx.iv(p);
        let q = [pi[0], pi[1]];
        match e {
            EntityRef::Edge(i) if cx.dim == 2 && (i as usize) < cx.edges.len() => {
                Some(cx.wv(to3(cx.project_edge(i, q))))
            }
            EntityRef::Face(i) if cx.dim == 3 && (i as usize) < cx.faces.len() => Some(cx.facet_project(i, p)),
            EntityRef::Edge(i) if cx.dim == 3 && (i as usize) < cx.edges3.len() => Some(cx.edge3_closest(i, p).0),
            EntityRef::Vertex(i) => cx.vertices.get(i as usize).map(|v| v.p),
            EntityRef::Region(r) if cx.dim == 2 && (r as usize) < cx.regions.len() => {
                // inside: the point itself; outside: the nearest boundary point of the region
                if self.ev.classify(p) == Some(r) {
                    return Some(p);
                }
                let mut best: Option<(f64, V2)> = None;
                for &l in &cx.regions[r as usize].loops {
                    for &(ed, _) in &cx.loops[l as usize].edges {
                        let c = cx.project_edge(ed, q);
                        let d = dist2(c, q);
                        if best.is_none_or(|b| d < b.0) {
                            best = Some((d, c));
                        }
                    }
                }
                best.map(|b| cx.wv(to3(b.1)))
            }
            _ => None,
        }
    }
}

fn to_world_box(perm: &[usize; 3], b: Aabb) -> Aabb {
    if b.is_empty() {
        return b;
    }
    let (mut lo, mut hi) = ([0.0; 3], [0.0; 3]);
    for i in 0..3 {
        lo[perm[i]] = b.lo[i];
        hi[perm[i]] = b.hi[i];
    }
    Aabb { lo, hi }
}
