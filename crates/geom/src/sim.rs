//! Simulation operations as *data*: absorbing-layer regions, mesh size rules, periodic face pairs (verified congruent), and the
//! projection oracle. geom describes which entities and how; femcore decides what a condition means and applies it.

use crate::error::{GeomError, GeomResult};
use crate::eval::Evaluated;
use crate::math::*;
use crate::model::*;
use crate::select::{EntityKind, EntitySet};

pub use crate::model::{PmlCoordKind, PmlGrowMode, PmlRole};

#[derive(Clone, Debug, PartialEq)]
pub struct PmlRegion {
    pub region: u32,
    /// the physical region this layer was cut or extruded from
    pub parent: Option<u32>,
    pub window: Aabb,
    /// per axis -1 / 0 / +1: the window side(s) the layer attaches to (a corner piece has two nonzero entries)
    pub sides: [i8; 3],
    pub thickness: f64,
    pub layers: f64,
    pub grading: f64,
    pub coord: PmlCoordKind,
    pub mode: PmlGrowMode,
    pub center: [f64; 3],
    pub axis: [f64; 3],
    pub r_in: f64,
    pub z_range: [f64; 2],
    pub layer_grading: f64,
    pub role: PmlRole,
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
        self.expand_outward_pmls();
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
                    coord: p.coord,
                    mode: p.mode,
                    center: to_world_vec(&cx_perm, p.center),
                    axis: to_world_vec(&cx_perm, p.axis),
                    r_in: p.r_in,
                    z_range: p.z_range,
                    layer_grading: p.layer_grading,
                    role: p.role,
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
        self.validate_pml_ray_invariance()?;
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

    /// Expands outward boundary-extruded `OutwardPmlSpec` entries against the physical `Complex`,
    /// attaching virtual outward `Region`s to the physical regions that touch each active boundary piece.
    fn expand_outward_pmls(&mut self) {
        let specs = self.scene().outward_pmls.clone();
        if specs.is_empty() {
            return;
        }
        let dim = self.list.dim as usize;
        let scene_axis = self.scene().axis;
        let cx_perm = crate::scene::axis_perm(scene_axis);
        let all = [Side::XMin, Side::XMax, Side::YMin, Side::YMax, Side::ZMin, Side::ZMax];

        for spec in specs {
            let world_win = to_world_box(&cx_perm, spec.window);
            let world_center = to_world_vec(&cx_perm, spec.center);
            let tol = (1e-5 * world_win.diagonal()).max(1e-12);
            let n_existing = self.complex().regions.len();
            let cand_regs: Vec<u32> = (0..n_existing as u32)
                .filter(|&r_id| {
                    let r = &self.complex().regions[r_id as usize];
                    r.parent.is_none()
                        && matches!(r.kind, RegionKind::Material)
                        && (spec.regions.is_empty() || spec.regions.contains(&r.name))
                })
                .collect();

            for r_id in cand_regs {
                type OutwardPiece = ([i8; 3], String, PmlRole, Aabb, [f64; 2]);
                let mut pieces: Vec<OutwardPiece> = vec![];
                if spec.coord == PmlCoordKind::Spherical {
                    if dim == 3 {
                        let touches = self.complex().faces.iter().any(|f| {
                            !f.is_conductor_cavity()
                                && (f.first.is_none() != f.second.is_none())
                                && f.first.or(f.second) == Some(r_id)
                        });
                        if touches {
                            let r_out = spec.r_in + spec.thickness;
                            let bb = Aabb::new(sub3(world_center, [r_out; 3]), add3(world_center, [r_out; 3]));
                            pieces.push((
                                [1, 1, 1],
                                "spherical".to_string(),
                                PmlRole::SphericalShell,
                                bb,
                                spec.z_range,
                            ));
                        }
                    } else {
                        let touches = self.complex().edges.iter().any(|e| {
                            (e.left.is_none() != e.right.is_none())
                                && e.left.or(e.right) == Some(r_id)
                                && e.bbox.hi[0] > tol
                        });
                        if touches {
                            let r_out = spec.r_in + spec.thickness;
                            let bb =
                                Aabb::new([0.0, 0.0, world_center[2] - r_out], [r_out, 0.0, world_center[2] + r_out]);
                            pieces.push((
                                [1, 1, 0],
                                "spherical".to_string(),
                                PmlRole::SphericalShell,
                                bb,
                                spec.z_range,
                            ));
                        }
                    }
                } else if spec.coord == PmlCoordKind::Cylindrical && spec.curved_outer {
                    if dim == 2 {
                        let touches = self
                            .complex()
                            .edges
                            .iter()
                            .any(|e| (e.left.is_none() != e.right.is_none()) && e.left.or(e.right) == Some(r_id));
                        if touches {
                            let r_out = spec.r_in + spec.thickness;
                            let bb = Aabb::new(
                                [world_center[0] - r_out, 0.0, world_center[2] - r_out],
                                [world_center[0] + r_out, 0.0, world_center[2] + r_out],
                            );
                            pieces.push(([1, 0, 0], "radial".to_string(), PmlRole::RadialShell, bb, spec.z_range));
                        }
                    } else {
                        let ax_idx = scene_axis;
                        let (u_ax, v_ax) = (cx_perm[0], cx_perm[1]);
                        let (cap_min_name, cap_max_name) = match ax_idx {
                            0 => ("xmin", "xmax"),
                            1 => ("ymin", "ymax"),
                            _ => ("zmin", "zmax"),
                        };
                        let (min_side, max_side) = (all[2 * ax_idx], all[2 * ax_idx + 1]);
                        let radial_active = spec.sides.iter().any(|&s| s != min_side && s != max_side);
                        let min_cap_active = spec.sides.contains(&min_side);
                        let max_cap_active = spec.sides.contains(&max_side);
                        let (z_min, z_max) = (spec.z_range[0], spec.z_range[1]);

                        let mut touches_barrel = false;
                        let mut r_z0 = f64::INFINITY;
                        let mut r_z1 = f64::NEG_INFINITY;
                        let mut touches_zmin = false;
                        let mut touches_zmax = false;
                        let cx = self.complex();
                        for (fi, f) in cx.faces.iter().enumerate() {
                            if f.is_conductor_cavity()
                                || (f.first.is_none() == f.second.is_none())
                                || f.first.or(f.second) != Some(r_id)
                            {
                                continue;
                            }
                            if (f.bbox.lo[ax_idx] - z_min).abs() <= tol && (f.bbox.hi[ax_idx] - z_min).abs() <= tol {
                                touches_zmin = true;
                            } else if (f.bbox.lo[ax_idx] - z_max).abs() <= tol
                                && (f.bbox.hi[ax_idx] - z_max).abs() <= tol
                            {
                                touches_zmax = true;
                            } else if cx.facet_kind(fi as u32) == "cylinder"
                                && cx.facet_radius(fi as u32).is_some_and(|r| (r - spec.r_in).abs() <= tol)
                            {
                                touches_barrel = true;
                                r_z0 = r_z0.min(f.bbox.lo[ax_idx]);
                                r_z1 = r_z1.max(f.bbox.hi[ax_idx]);
                            }
                        }
                        let cyl_bb = |r_ext: f64, za: f64, zb: f64| -> Aabb {
                            let (mut lo, mut hi) = ([0.0; 3], [0.0; 3]);
                            lo[u_ax] = world_center[u_ax] - r_ext;
                            hi[u_ax] = world_center[u_ax] + r_ext;
                            lo[v_ax] = world_center[v_ax] - r_ext;
                            hi[v_ax] = world_center[v_ax] + r_ext;
                            lo[ax_idx] = za;
                            hi[ax_idx] = zb;
                            Aabb { lo, hi }
                        };
                        if radial_active && touches_barrel {
                            pieces.push((
                                [0, 0, 0],
                                "radial".to_string(),
                                PmlRole::RadialShell,
                                cyl_bb(spec.r_in + spec.thickness, r_z0, r_z1),
                                [r_z0, r_z1],
                            ));
                        }
                        if min_cap_active && touches_zmin {
                            pieces.push((
                                [0, 0, -1],
                                cap_min_name.to_string(),
                                PmlRole::EndCap,
                                cyl_bb(spec.r_in, z_min - spec.thickness, z_min),
                                [z_min, z_max],
                            ));
                        }
                        if max_cap_active && touches_zmax {
                            pieces.push((
                                [0, 0, 1],
                                cap_max_name.to_string(),
                                PmlRole::EndCap,
                                cyl_bb(spec.r_in, z_max, z_max + spec.thickness),
                                [z_min, z_max],
                            ));
                        }
                        if radial_active
                            && min_cap_active
                            && touches_barrel
                            && (touches_zmin || (r_z0 - z_min).abs() <= tol)
                        {
                            pieces.push((
                                [0, 0, -1],
                                format!("radial_{cap_min_name}"),
                                PmlRole::CornerAnnulus,
                                cyl_bb(spec.r_in + spec.thickness, z_min - spec.thickness, z_min),
                                [z_min, z_max],
                            ));
                        }
                        if radial_active
                            && max_cap_active
                            && touches_barrel
                            && (touches_zmax || (r_z1 - z_max).abs() <= tol)
                        {
                            pieces.push((
                                [0, 0, 1],
                                format!("radial_{cap_max_name}"),
                                PmlRole::CornerAnnulus,
                                cyl_bb(spec.r_in + spec.thickness, z_max, z_max + spec.thickness),
                                [z_min, z_max],
                            ));
                        }
                    }
                } else {
                    let mut combos: Vec<[i8; 3]> = vec![];
                    let total = 3usize.pow(dim as u32);
                    for code in 0..total {
                        let mut c = [0i8; 3];
                        let mut k = code;
                        for ci in c.iter_mut().take(dim) {
                            *ci = (k % 3) as i8 - 1;
                            k /= 3;
                        }
                        let ok = c != [0, 0, 0]
                            && (0..dim).all(|a| {
                                c[a] == 0
                                    || spec.sides.contains(&if c[a] < 0 {
                                        all[2 * cx_perm[a]]
                                    } else {
                                        all[2 * cx_perm[a] + 1]
                                    })
                            });
                        if ok {
                            combos.push(c);
                        }
                    }
                    combos.sort_by_key(|c| {
                        let n = c.iter().filter(|&&x| x != 0).count();
                        let single = (0..dim).find(|&a| c[a] != 0).map_or(0, |a| 2 * a + (c[a] > 0) as usize);
                        (n, if n == 1 { single } else { 0 }, (0..dim).rev().map(|a| c[a] + 1).collect::<Vec<_>>())
                    });
                    for c in &combos {
                        let mut w_sides = [0i8; 3];
                        for i in 0..3 {
                            w_sides[cx_perm[i]] = c[i];
                        }
                        let mut toks: Vec<(usize, String)> = (0..dim)
                            .filter(|&a| c[a] != 0)
                            .map(|a| {
                                (
                                    cx_perm[a],
                                    format!("{}{}", ["x", "y", "z"][cx_perm[a]], if c[a] < 0 { "min" } else { "max" }),
                                )
                            })
                            .collect();
                        toks.sort();
                        let piece = toks.into_iter().map(|t| t.1).collect::<Vec<_>>().join("_");
                        let nz = c.iter().filter(|&&x| x != 0).count();
                        let role = if dim == 2 && spec.coord == PmlCoordKind::Cylindrical {
                            match (c[0], c[1]) {
                                (1, 0) => PmlRole::RadialShell,
                                (0, _) => PmlRole::EndCap,
                                _ => PmlRole::CornerAnnulus,
                            }
                        } else if dim == 2 {
                            if nz == 1 {
                                PmlRole::Face
                            } else {
                                PmlRole::Corner
                            }
                        } else {
                            match nz {
                                1 => PmlRole::Face,
                                2 => PmlRole::Edge,
                                _ => PmlRole::Corner,
                            }
                        };

                        if dim == 2 {
                            let ext_edges: Vec<&crate::complex::Edge> = self
                                .complex()
                                .edges
                                .iter()
                                .filter(|e| (e.left.is_none() != e.right.is_none()) && e.left.or(e.right) == Some(r_id))
                                .collect();
                            if nz == 1 {
                                let ax = if w_sides[0] != 0 { 0 } else { 2 };
                                let tan_ax = 2 - ax;
                                let sgn = w_sides[ax];
                                let plane_val = if sgn < 0 { world_win.lo[ax] } else { world_win.hi[ax] };
                                let matching: Vec<&&crate::complex::Edge> = ext_edges
                                    .iter()
                                    .filter(|e| {
                                        (e.bbox.lo[ax] - plane_val).abs() <= tol
                                            && (e.bbox.hi[ax] - plane_val).abs() <= tol
                                    })
                                    .collect();
                                if !matching.is_empty() {
                                    let t_lo = matching.iter().map(|e| e.bbox.lo[tan_ax]).fold(f64::INFINITY, f64::min);
                                    let t_hi =
                                        matching.iter().map(|e| e.bbox.hi[tan_ax]).fold(f64::NEG_INFINITY, f64::max);
                                    let (n_lo, n_hi) = if sgn < 0 {
                                        (plane_val - spec.thickness, plane_val)
                                    } else {
                                        (plane_val, plane_val + spec.thickness)
                                    };
                                    let (mut lo, mut hi) = ([0.0; 3], [0.0; 3]);
                                    lo[ax] = n_lo;
                                    hi[ax] = n_hi;
                                    lo[tan_ax] = t_lo;
                                    hi[tan_ax] = t_hi;
                                    pieces.push((*c, piece, role, Aabb { lo, hi }, spec.z_range));
                                }
                            } else {
                                let cx_val = if w_sides[0] < 0 { world_win.lo[0] } else { world_win.hi[0] };
                                let cz_val = if w_sides[2] < 0 { world_win.lo[2] } else { world_win.hi[2] };
                                let touches = ext_edges.iter().any(|e| {
                                    e.bbox.lo[0] <= cx_val + tol
                                        && e.bbox.hi[0] >= cx_val - tol
                                        && e.bbox.lo[2] <= cz_val + tol
                                        && e.bbox.hi[2] >= cz_val - tol
                                });
                                if touches {
                                    let (x0, x1) = if w_sides[0] < 0 {
                                        (cx_val - spec.thickness, cx_val)
                                    } else {
                                        (cx_val, cx_val + spec.thickness)
                                    };
                                    let (z0, z1) = if w_sides[2] < 0 {
                                        (cz_val - spec.thickness, cz_val)
                                    } else {
                                        (cz_val, cz_val + spec.thickness)
                                    };
                                    pieces.push((
                                        *c,
                                        piece,
                                        role,
                                        Aabb { lo: [x0, 0.0, z0], hi: [x1, 0.0, z1] },
                                        spec.z_range,
                                    ));
                                }
                            }
                        } else {
                            let ext_faces: Vec<&crate::complex::Face> = self
                                .complex()
                                .faces
                                .iter()
                                .filter(|f| {
                                    !f.is_conductor_cavity()
                                        && (f.first.is_none() != f.second.is_none())
                                        && f.first.or(f.second) == Some(r_id)
                                })
                                .collect();
                            let matching: Vec<&&crate::complex::Face> = ext_faces
                                .iter()
                                .filter(|f| {
                                    let touches_active = (0..3).filter(|&a| w_sides[a] != 0).all(|a| {
                                        let val = if w_sides[a] < 0 { world_win.lo[a] } else { world_win.hi[a] };
                                        f.bbox.lo[a] <= val + tol && f.bbox.hi[a] >= val - tol
                                    });
                                    if !touches_active {
                                        return false;
                                    }
                                    if nz == 1 {
                                        let a = (0..3).find(|&a| w_sides[a] != 0).unwrap();
                                        let val = if w_sides[a] < 0 { world_win.lo[a] } else { world_win.hi[a] };
                                        (f.bbox.lo[a] - val).abs() <= tol && (f.bbox.hi[a] - val).abs() <= tol
                                    } else if nz == 2 {
                                        let a_free = (0..3).find(|&a| w_sides[a] == 0).unwrap();
                                        f.bbox.hi[a_free] - f.bbox.lo[a_free] > tol
                                    } else {
                                        true
                                    }
                                })
                                .collect();
                            if !matching.is_empty() {
                                let (mut lo, mut hi) = ([0.0; 3], [0.0; 3]);
                                for a in 0..3 {
                                    if w_sides[a] < 0 {
                                        lo[a] = world_win.lo[a] - spec.thickness;
                                        hi[a] = world_win.lo[a];
                                    } else if w_sides[a] > 0 {
                                        lo[a] = world_win.hi[a];
                                        hi[a] = world_win.hi[a] + spec.thickness;
                                    } else {
                                        lo[a] = matching.iter().map(|f| f.bbox.lo[a]).fold(f64::INFINITY, f64::min);
                                        hi[a] = matching.iter().map(|f| f.bbox.hi[a]).fold(f64::NEG_INFINITY, f64::max);
                                    }
                                }
                                pieces.push((*c, piece, role, Aabb { lo, hi }, spec.z_range));
                            }
                        }
                    }
                }

                let phys_reg = self.complex().regions[r_id as usize].clone();
                for (c, piece, role, world_bbox, piece_z_range) in pieces {
                    let scene_idx = self.scene().regions.len();
                    let reg_name = format!("{}.pml_{piece}", phys_reg.name);
                    self.scene.as_mut().unwrap().regions.push(crate::scene::RegionX {
                        name: reg_name.clone(),
                        kind: RegionKind::Pml,
                        feature: phys_reg.feature,
                        pinned_id: None,
                        expr: crate::scene::Node::Empty,
                        parent: Some(phys_reg.name.clone()),
                        pml: Some(crate::scene::PmlInfo {
                            feature: spec.feature,
                            sides: c,
                            window: spec.window,
                            thickness: spec.thickness,
                            layers: spec.layers,
                            grading: spec.grading,
                            coord: spec.coord,
                            mode: PmlGrowMode::Outward,
                            center: spec.center,
                            axis: spec.axis,
                            r_in: spec.r_in,
                            z_range: piece_z_range,
                            layer_grading: spec.layer_grading,
                            role,
                        }),
                    });
                    let cx = self.complex.as_mut().unwrap();
                    let new_reg_id = cx.regions.len() as u32;
                    cx.regions.push(crate::complex::Region {
                        name: reg_name,
                        kind: RegionKind::Pml,
                        feature: phys_reg.feature,
                        pinned_id: None,
                        scene_index: scene_idx,
                        measure: 0.0,
                        bbox: world_bbox,
                        loops: vec![],
                        parent: Some(r_id),
                    });
                    cx.scene_to_region.push(Some(new_reg_id));
                }
            }
        }
        self.complex.as_mut().unwrap().rebuild_bvh();
    }

    /// Validates Teixeira-Chew / Collino-Monk analytic-continuation ray invariance (`d_u eps == 0`)
    /// across all absorbing-layer boundaries (`Spherical` and `Cylindrical`).
    pub fn validate_pml_ray_invariance(&self) -> GeomResult<()> {
        if self.sim.pml.is_empty() {
            return Ok(());
        }
        let cx = self.complex();
        let scene = self.scene();
        let dim = self.list.dim as usize;
        let tol = (1e-5 * self.bounds().diagonal()).max(1e-12);

        // 1. Validate Spherical PMLs
        let spherical: Vec<&PmlRegion> = self.sim.pml.iter().filter(|p| p.coord == PmlCoordKind::Spherical).collect();
        if !spherical.is_empty() {
            let p0 = spherical[0];
            let center = p0.center;
            let r_in = p0.r_in;
            let mut parents: Vec<u32> = spherical.iter().filter_map(|p| p.parent).collect();
            parents.sort_unstable();
            parents.dedup();

            // Check if any off-center or piercing sheet reaches the spherical boundary
            let sheet_violates = dim == 3
                && scene.sheets.iter().any(|s| {
                    let sb = s.expr.bbox(&scene.prims);
                    if sb.is_empty() {
                        return false;
                    }
                    let dz = s.z - center[2];
                    if dz.abs() > r_in + tol {
                        return false;
                    }
                    let r_slice = (r_in * r_in - dz * dz).max(0.0).sqrt();
                    let d_max = [sb.lo[0], sb.hi[0]]
                        .into_iter()
                        .flat_map(|x| {
                            [sb.lo[1], sb.hi[1]].into_iter().map(move |y| (x - center[0]).hypot(y - center[1]))
                        })
                        .fold(0.0_f64, f64::max);
                    d_max >= r_slice - tol
                });

            // In 2D axisymmetric spherical, check if multiple materials touch the semicircle across non-radial interfaces
            let mut non_radial_2d = false;
            if dim == 2 && parents.len() > 1 {
                for (ei, e) in cx.edges.iter().enumerate() {
                    let (Some(l), Some(r)) = (e.left, e.right) else { continue };
                    if cx.regions[l as usize].kind == RegionKind::Pml || cx.regions[r as usize].kind == RegionKind::Pml
                    {
                        continue;
                    }
                    // Check if this interior material interface touches the semicircle R = r_in
                    let (c_curve, t0, t1) = cx.curve_of(ei as u32);
                    let p_start = cx.wv(to3(c_curve.eval(t0)));
                    let p_end = cx.wv(to3(c_curve.eval(t1)));
                    let d_start = ((p_start[0] - center[0]).hypot(p_start[2] - center[2]) - r_in).abs();
                    let d_end = ((p_end[0] - center[0]).hypot(p_end[2] - center[2]) - r_in).abs();
                    if d_start <= tol || d_end <= tol {
                        let p_hit = if d_start <= tol { p_start } else { p_end };
                        let r_hat = unit3([p_hit[0] - center[0], 0.0, p_hit[2] - center[2]]);
                        for n in cx.facet_normals(ei as u32) {
                            if dot3(n, r_hat).abs() > 1e-3 {
                                non_radial_2d = true;
                            }
                        }
                    }
                }
            }

            // In 3D spherical, > 1 physical material touching R = r_in or > 1 horizontal stack layer cutting the outer sphere
            let multi_mat_3d = dim == 3
                && (parents.len() > 1
                    || cx.faces.iter().any(|f| {
                        let (Some(a), Some(b)) = (f.first, f.second) else { return false };
                        if cx.regions[a as usize].kind == RegionKind::Pml
                            || cx.regions[b as usize].kind == RegionKind::Pml
                        {
                            return false;
                        }
                        // An interior interface whose bbox reaches the outer sphere radius R = r_in violates d_R eps == 0
                        let fb = f.bbox;
                        let d_max = [fb.lo[0], fb.hi[0]]
                            .into_iter()
                            .flat_map(|x| {
                                [fb.lo[1], fb.hi[1]].into_iter().flat_map(move |y| {
                                    [fb.lo[2], fb.hi[2]].into_iter().map(move |z| norm3(sub3([x, y, z], center)))
                                })
                            })
                            .fold(0.0_f64, f64::max);
                        d_max >= r_in - tol && !matches!(f.geom, crate::complex::FaceGeom::Quadric( crate::prim::Body::Sphere { c, .. }) if norm3(sub3(cx.wv(c), center)) <= tol)
                    }));

            if sheet_violates || non_radial_2d || multi_mat_3d {
                return Err(GeomError(format!(
                    "Spherical PML violates analytic-continuation ray invariance (d_R eps != 0): multiple material layers or off-center interfaces/sheets intersect the spherical boundary R = {:.3e} m around center ({:.3e}, {:.3e}, {:.3e}). Use Cylindrical PML (coord = Cylindrical) for horizontally layered wafers/disks, or Cartesian PML (coord = Cartesian) when waveguides/conductors pierce the boundary.",
                    r_in, center[0], center[1], center[2]
                )));
            }
        }

        // 2. Validate Cylindrical PMLs
        let cylindrical: Vec<&PmlRegion> =
            self.sim.pml.iter().filter(|p| p.coord == PmlCoordKind::Cylindrical && p.r_in > 0.0).collect();
        if !cylindrical.is_empty() {
            let p0 = cylindrical[0];
            let center = p0.center;
            let r_in = p0.r_in;

            if dim == 2 {
                // Only check circular 2D planar cylindrical PMLs (where axis is out-of-plane Y: axis[1].abs() > 0.9)
                if p0.axis[1].abs() > 0.9 {
                    let mut non_radial_2d = false;
                    for (ei, e) in cx.edges.iter().enumerate() {
                        let (Some(l), Some(r)) = (e.left, e.right) else { continue };
                        if cx.regions[l as usize].kind == RegionKind::Pml
                            || cx.regions[r as usize].kind == RegionKind::Pml
                        {
                            continue;
                        }
                        let (c_curve, t0, t1) = cx.curve_of(ei as u32);
                        let p_start = cx.wv(to3(c_curve.eval(t0)));
                        let p_end = cx.wv(to3(c_curve.eval(t1)));
                        let d_start = ((p_start[0] - center[0]).hypot(p_start[2] - center[2]) - r_in).abs();
                        let d_end = ((p_end[0] - center[0]).hypot(p_end[2] - center[2]) - r_in).abs();
                        if d_start <= tol || d_end <= tol {
                            let p_hit = if d_start <= tol { p_start } else { p_end };
                            let r_hat = unit3([p_hit[0] - center[0], 0.0, p_hit[2] - center[2]]);
                            for n in cx.facet_normals(ei as u32) {
                                if dot3(n, r_hat).abs() > 1e-3 {
                                    non_radial_2d = true;
                                }
                            }
                        }
                    }
                    if non_radial_2d {
                        return Err(GeomError(format!(
                            "Cylindrical PML violates analytic-continuation ray invariance (d_r eps != 0): a non-radial material interface or layered substrate plane intersects the 2D circular boundary r = {:.3e} m. Use Cartesian PML (coord = Cartesian) for 2D layered-wafer cross-sections or piercing waveguides.",
                            r_in
                        )));
                    }
                }
            } else if let Some(st) = &cx.stack {
                // 3D Cylindrical PML: horizontally stratified z-layers (d_r eps(z) == 0) are valid,
                // but any waveguide, strip, or sheet piercing the barrel wall r = r_in across azimuth angle theta is invalid.
                let c_iv = cx.iv(center);
                let c2 = [c_iv[0], c_iv[1]];
                // Find all overlay edges lying on the barrel circle r = r_in
                let barrel_edges: Vec<(usize, bool)> = st
                    .overlay
                    .edges
                    .iter()
                    .enumerate()
                    .filter_map(|(ei, oe)| {
                        let curve = &st.overlay.curves[oe.curve as usize];
                        let mid = curve.eval(0.5 * (oe.t0 + oe.t1));
                        let p0 = curve.eval(oe.t0);
                        let p1 = curve.eval(oe.t1);
                        if (dist2(mid, c2) - r_in).abs() <= tol
                            && (dist2(p0, c2) - r_in).abs() <= tol
                            && (dist2(p1, c2) - r_in).abs() <= tol
                        {
                            let n2 = st.overlay.edge_normal(ei as u32, 0.5 * (oe.t0 + oe.t1));
                            let r_dir = sub2(mid, c2);
                            // n2 points from left to right; if n2 . r_dir > 0, left is the interior side (< r_in)
                            Some((ei, dot2(n2, r_dir) >= 0.0))
                        } else {
                            None
                        }
                    })
                    .collect();

                if !barrel_edges.is_empty() {
                    let z_min = cylindrical.iter().map(|p| p.z_range[0]).fold(f64::INFINITY, f64::min);
                    let z_max = cylindrical.iter().map(|p| p.z_range[1]).fold(f64::NEG_INFINITY, f64::max);
                    for l in &st.layers {
                        if l.z1 <= z_min + tol || l.z0 >= z_max - tol {
                            continue;
                        }
                        let mut layer_barrel_regs: Vec<Option<u32>> = vec![];
                        for &(ei, left_is_inner) in &barrel_edges {
                            let oe = &st.overlay.edges[ei];
                            let inner_cell = if left_is_inner { oe.left } else { oe.right };
                            let reg = inner_cell.and_then(|c| l.region_of[c as usize]);
                            if !layer_barrel_regs.contains(&reg) {
                                layer_barrel_regs.push(reg);
                            }
                        }
                        if layer_barrel_regs.len() > 1 {
                            return Err(GeomError(format!(
                                "Cylindrical PML violates analytic-continuation ray invariance (d_r eps != 0): a waveguide or non-axisymmetric region pierces the cylindrical barrel wall r = {:.3e} m in layer z = [{:.3e}, {:.3e}] m. Use Cartesian PML (coord = Cartesian) when waveguides or strips pierce the boundary.",
                                r_in, l.z0, l.z1
                            )));
                        }
                    }
                }

                // Check if any Op::Sheet trace/strip pierces the barrel wall r = r_in over a partial azimuth arc
                for s in &scene.sheets {
                    let mut n_in = 0usize;
                    let mut n_out = 0usize;
                    for &(ei, _) in &barrel_edges {
                        let oe = &st.overlay.edges[ei];
                        let curve = &st.overlay.curves[oe.curve as usize];
                        let mid = curve.eval(0.5 * (oe.t0 + oe.t1));
                        let r_dir = sub2(mid, c2);
                        let d = r_dir[0].hypot(r_dir[1]).max(1e-30);
                        let p2 = [mid[0] - 4.0 * tol * r_dir[0] / d, mid[1] - 4.0 * tol * r_dir[1] / d];
                        let hit = s.expr.eval(&|pi| scene.prims[pi].sdf_near(p2, 0.0) < 0.0);
                        if hit {
                            n_in += 1;
                        } else {
                            n_out += 1;
                        }
                    }
                    for k in 0..72 {
                        let th = std::f64::consts::TAU * (k as f64 + 0.5) / 72.0;
                        let p2 = [c2[0] + (r_in - 4.0 * tol) * th.cos(), c2[1] + (r_in - 4.0 * tol) * th.sin()];
                        let hit = s.expr.eval(&|pi| scene.prims[pi].sdf_near(p2, 0.0) < 0.0);
                        if hit {
                            n_in += 1;
                        } else {
                            n_out += 1;
                        }
                    }
                    if n_in > 0 && n_out > 0 {
                        return Err(GeomError(format!(
                            "Cylindrical PML violates analytic-continuation ray invariance (d_r eps != 0): conductive Sheet '{}' at z = {:.3e} m pierces the cylindrical barrel wall r = {:.3e} m over a partial azimuth arc. Use Cartesian PML (coord = Cartesian) when CPW lines or conductive traces pierce the boundary.",
                            s.name, s.z, r_in
                        )));
                    }
                }
            }
        }
        Ok(())
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

fn to_world_vec(perm: &[usize; 3], v: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[perm[i]] = v[i];
    }
    out
}
