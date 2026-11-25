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

use crate::mesh::{Mesh, SimplexMesh};
use crate::meshing::mesher::{mesh_chains, Chain, Fill, MeshParams};
use crate::props::{Pml, PmlCoord, PmlDir};
use geom::eval::{Evaluated, Tier};
use geom::labels::LabelMap;
use geom::model::{FeatureList, PmlCoordKind, PmlGrowMode, PmlRole, RegionKind};
use geom::query::DiscretizeOpts;
use geom::select::{EntityKind, Selector};
use std::collections::BTreeMap;

/// Sizes of a mesh made from a model.
#[derive(Clone, Debug, PartialEq)]
pub struct GeomMeshOpts {
    /// global element size
    pub h: f64,
    /// largest distance of a chord from its curve (default `h / 50`)
    pub max_dev: Option<f64>,
    pub min_angle_deg: f64,
    /// 3D: give the cells on a curved surface quadratic (isoparametric) geometry (ADR 82). Off by default: straight-sided cells (ADR 18).
    pub curved: bool,
    /// 3D sheet mesher: automatic coplanar gap proximity sizing (`d_gap / gap_cells`); off by default (`None`).
    pub gap_cells: Option<f64>,
    /// 3D sheet mesher: conforming pure-tetrahedral vertical coarsening away from sheet planes; off by default (`false`).
    pub z_coarsen: bool,
    /// Maximum number of cells allowed before aborting meshing; off by default (`None`).
    pub max_cells: Option<usize>,
    /// 3D sheet mesher: floor on initial vertical z-slab thickness near sheet planes (m); off by default (`None`).
    pub dz_min: Option<f64>,
    /// 3D sheet mesher: maximum total number of horizontal z-slabs across active stack layers; off by default (`None`).
    pub max_z_slabs: Option<usize>,
    /// 3D sheet mesher: geometric growth ratio for z-slab thickness away from sheet planes (default `1.35` when `None`).
    pub z_growth_ratio: Option<f64>,
    /// 2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial vertices; off by default (`None`).
    pub max_steiner_multiplier: Option<usize>,
    /// Fail meshing if 2D Ruppert refinement terminates early or drops unhealed collinear sliver triangles; off by default (`false`).
    pub strict_mesh: bool,
    /// 3D sheet mesher: edge singularity enrichment (`"none"` or `"quarter_point"`); off by default (`None`).
    pub singular_edges: Option<String>,
    /// 3D sheet mesher: Barsoum quarter-point mid-edge parameter `xi_QP` in `[0.25, 0.50]`; default `0.25` when `None`.
    pub quarter_point_xi: Option<f64>,
    /// 3D sheet mesher: Morrison-Lewis / Babuska radical corner grading at sharp sheet corners; off by default (`false`).
    pub corner_grading: bool,
}

impl GeomMeshOpts {
    pub fn new(h: f64) -> Self {
        GeomMeshOpts {
            h,
            max_dev: None,
            min_angle_deg: 25.0,
            curved: false,
            gap_cells: None,
            z_coarsen: false,
            max_cells: None,
            dz_min: None,
            max_z_slabs: None,
            z_growth_ratio: None,
            max_steiner_multiplier: None,
            strict_mesh: false,
            singular_edges: None,
            quarter_point_xi: None,
            corner_grading: false,
        }
    }
}

/// A meshed model: the mesh, the uncut mixed-cell exchange mesh (when present), the compiled names,
/// and the evaluation they came from.
#[derive(Clone)]
pub struct MeshedGeom {
    pub mesh: SimplexMesh,
    /// Uncut mixed-cell `Mesh` (`CellMesh<2>` or `CellMesh<3>`) preserving `Quad`, `Prism`, and `Hex`
    /// absorbing-layer elements without eager simplex splitting.
    pub cells: Option<Mesh>,
    pub labels: LabelMap,
    /// shared, so that a cached mesh can be handed out many times without copying the evaluation
    pub eval: std::sync::Arc<Evaluated>,
    /// diagnostics and non-fatal warnings emitted during meshing (e.g. Ruppert vertex cap or sliver audit)
    pub warnings: Vec<String>,
}

impl MeshedGeom {
    /// Returns the uncut mixed-cell exchange [`Mesh`] (`CellMesh<2>` or `CellMesh<3>`) when available,
    /// falling back to [`Mesh::from_solver_mesh`] on `self.mesh` otherwise.
    pub fn into_exchange_mesh(&self) -> Result<Mesh, String> {
        match &self.cells {
            Some(m) => Ok(m.clone()),
            None => Mesh::from_solver_mesh(&self.mesh),
        }
    }

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

    /// The absorbing layers of the model as solver stretching, keyed by domain: `strength` is the imaginary stretch at the outer edge, the
    /// `grading` of the `Pml` feature is strictly the order of the polynomial profile (decoupled from geometric `layer_grading`), and
    /// `coord` / `bounds` carry exact `Cartesian`, `Cylindrical`, or `Spherical` geometry in wafer coordinates.
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
                let (coord, bounds) = match p.coord {
                    PmlCoordKind::Spherical => (
                        PmlCoord::Spherical {
                            center: p.center,
                            r_in: p.r_in,
                            thickness: p.thickness,
                        },
                        None,
                    ),
                    PmlCoordKind::Cylindrical if p.r_in > 0.0 => {
                        let ax_len = (p.axis[0] * p.axis[0] + p.axis[1] * p.axis[1] + p.axis[2] * p.axis[2]).sqrt();
                        let a_hat = if ax_len > 1e-30 {
                            [p.axis[0] / ax_len, p.axis[1] / ax_len, p.axis[2] / ax_len]
                        } else if self.mesh.dim == 2 {
                            [0.0, 1.0, 0.0]
                        } else {
                            [0.0, 0.0, 1.0]
                        };
                        let c_dot_a = p.center[0] * a_hat[0] + p.center[1] * a_hat[1] + p.center[2] * a_hat[2];
                        let s_dot = p.sides[0] as f64 * a_hat[0]
                            + p.sides[1] as f64 * a_hat[1]
                            + p.sides[2] as f64 * a_hat[2];
                        let radial_on = matches!(p.role, PmlRole::RadialShell | PmlRole::CornerAnnulus);
                        let axial = if matches!(p.role, PmlRole::EndCap | PmlRole::CornerAnnulus) {
                            if s_dot < -0.5 {
                                PmlDir::Minus
                            } else if s_dot > 0.5 {
                                PmlDir::Plus
                            } else {
                                PmlDir::Off
                            }
                        } else {
                            PmlDir::Off
                        };
                        let z_in = match axial {
                            PmlDir::Plus => p.z_range[1] - c_dot_a,
                            PmlDir::Minus => p.z_range[0] - c_dot_a,
                            PmlDir::Off => 0.0,
                        };
                        let r_thickness = if radial_on { p.thickness } else { 0.0 };
                        let z_thickness = if axial != PmlDir::Off { p.thickness } else { 0.0 };
                        (
                            PmlCoord::Cylindrical {
                                center: p.center,
                                axis: a_hat,
                                r_in: p.r_in,
                                r_thickness,
                                radial_on,
                                axial,
                                z_in,
                                z_thickness,
                            },
                            None,
                        )
                    }
                    _ => {
                        let reg_bb = self.eval.regions()[p.region as usize].bbox;
                        let mut lo = reg_bb.lo;
                        let mut hi = reg_bb.hi;
                        for a in 0..3 {
                            if p.sides[a] < 0
                                && (p.mode == PmlGrowMode::Outward
                                    || ((reg_bb.hi[a] - reg_bb.lo[a]) - p.thickness).abs() <= 1e-5 * p.thickness)
                            {
                                lo[a] = p.window.lo[a] - p.thickness;
                                hi[a] = p.window.lo[a];
                            } else if p.sides[a] > 0
                                && (p.mode == PmlGrowMode::Outward
                                    || ((reg_bb.hi[a] - reg_bb.lo[a]) - p.thickness).abs() <= 1e-5 * p.thickness)
                            {
                                lo[a] = p.window.hi[a];
                                hi[a] = p.window.hi[a] + p.thickness;
                            }
                        }
                        (PmlCoord::Cartesian, Some((lo, hi)))
                    }
                };
                (
                    d,
                    Pml {
                        x,
                        y,
                        z,
                        strength,
                        order: p.grading,
                        re: 1.0,
                        coord,
                        bounds,
                        ..Default::default()
                    },
                )
            })
            .collect()
    }
}

/// Evaluates a feature list to the tier the mesher needs.
pub fn evaluate(list: &FeatureList) -> Result<Evaluated, String> {
    geom::eval::Evaluator::new().eval(list, Tier::Complex).map_err(|e| e.to_string())
}

/// The element sizes a model asks for: per region (the global size, the `Size` rules on regions, and thickness / layers in an absorbing
/// piece) and per edge (a `Size` rule on edges). Shared by the simple mesher and by the adaptive one, which starts from the same sizes.
pub(crate) struct SizeMap {
    pub region_h: Vec<f64>,
    pub edge_rule: Vec<Option<f64>>,
}

pub(crate) fn size_map_with_skip(ev: &Evaluated, h: f64, skip: &[bool]) -> SizeMap {
    let mut region_h = vec![h; ev.regions().len()];
    for p in ev.pml_regions() {
        if !skip.get(p.region as usize).copied().unwrap_or(false) {
            region_h[p.region as usize] = region_h[p.region as usize].min(p.thickness / p.layers);
        }
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
    SizeMap { region_h, edge_rule }
}

/// Meshes a 2D model (the x-z cross-section) with triangles.
pub fn mesh_geom_2d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String> {
    mesh_geom_2d_with_skip(ev, opts, &[])
}

/// Meshes a 2D model (the x-z cross-section) with triangles, leaving out any regions marked in `skip`
/// (e.g. `Inward` PML domains when structured `Quad` layers will be built around the physical window).
pub fn mesh_geom_2d_with_skip(
    ev: Evaluated,
    opts: &GeomMeshOpts,
    skip: &[bool],
) -> Result<MeshedGeom, String> {
    if ev.list.dim != 2 {
        return Err("mesh_geom_2d needs a 2D model".into());
    }
    if !(opts.h > 0.0 && opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    if let Some(mc) = opts.max_cells {
        if mc < 1 {
            return Err("max_cells must be at least 1".into());
        }
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let max_dev = opts.max_dev.unwrap_or(opts.h / 50.0);

    let SizeMap { region_h, edge_rule } = size_map_with_skip(&ev, opts.h, skip);
    let skipped = |r: usize| skip.get(r).copied().unwrap_or(false);
    let meshable: Vec<bool> = labels
        .domains
        .iter()
        .map(|d| !skipped(d.region as usize) && !matches!(d.kind, RegionKind::Hole | RegionKind::Void))
        .collect();
    let has_skip = skip.iter().any(|&s| s);

    // chains: one per active edge, with the end points of the vertices so that neighbours share them exactly
    let w2 = |p: [f64; 3]| [p[0], p[2]];
    let mut chains = Vec::with_capacity(ev.edges().len());
    for (i, ed) in ev.edges().iter().enumerate() {
        if has_skip
            && ![ed.left, ed.right]
                .into_iter()
                .flatten()
                .any(|r| meshable.get(r as usize).copied().unwrap_or(false))
        {
            continue;
        }
        let mut h = [ed.left, ed.right]
            .iter()
            .flatten()
            .filter(|&&r| meshable.get(r as usize).copied().unwrap_or(false))
            .map(|&r| region_h[r as usize])
            .fold(opts.h, f64::min);
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

    // domain of a region; holes, voids, and skipped PML regions have none
    let domain_of: Vec<Option<i32>> = labels
        .domains
        .iter()
        .map(|d| meshable[d.region as usize].then_some(d.id))
        .collect();
    let (lo, hi) = if has_skip {
        let mut bb = geom::math::Aabb::EMPTY;
        for (r, &m) in meshable.iter().enumerate() {
            if m {
                bb = bb.union(&ev.regions()[r].bbox);
            }
        }
        if bb.is_empty() {
            let b = ev.bounds();
            ([b.lo[0], b.lo[2]], [b.hi[0], b.hi[2]])
        } else {
            ([bb.lo[0], bb.lo[2]], [bb.hi[0], bb.hi[2]])
        }
    } else {
        let b = ev.bounds();
        ([b.lo[0], b.lo[2]], [b.hi[0], b.hi[2]])
    };

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
    let h_fine_min = (0..ev.regions().len())
        .filter(|&r| domain_of[r].is_some())
        .map(|r| region_h[r])
        .fold(opts.h, f64::min);
    let domain_at = |c: [f64; 2]| ev.classify([c[0], 0.0, c[1]]).and_then(|r| domain_of[r as usize]);
    let (mesh, mc_rep) = mesh_chains(
        &chains,
        &fills,
        &MeshParams {
            h: opts.h,
            h_fine_min,
            min_angle_deg: opts.min_angle_deg,
            lo,
            hi,
            max_steiner_multiplier: opts.max_steiner_multiplier,
            strict_mesh: opts.strict_mesh,
            max_cells: opts.max_cells,
        },
        &domain_at,
    )?;
    drop(fills);
    if let Some(mc) = opts.max_cells {
        if mesh.n_cells() > mc {
            return Err(format!(
                "2D mesh exceeds max_cells limit: cells={} > max_cells={mc}",
                mesh.n_cells()
            ));
        }
    }
    let cells = Mesh::from_solver_mesh(&mesh)?;
    Ok(MeshedGeom {
        mesh,
        cells: Some(cells),
        labels,
        eval: std::sync::Arc::new(ev),
        warnings: mc_rep.warnings,
    })
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

/// A rectangle `[0, lx] x [0, lz]` of the x-z plane as one domain, meshed with element size `h`, with the boundary labels of its four sides
/// (`xmin`, `xmax`, `zmin` = bottom, `zmax` = top): the quick problem for tests, benchmarks and studies that need no named geometry.
pub struct Plate {
    pub mesh: SimplexMesh,
    pub domain: i32,
    pub xmin: Vec<i32>,
    pub xmax: Vec<i32>,
    pub zmin: Vec<i32>,
    pub zmax: Vec<i32>,
}

pub fn plate(lx: f64, lz: f64, h: f64) -> Result<Plate, String> {
    let mut m = geom::model::Model::new_2d();
    let r = m.rect("plate", (0.0, 0.0), (lx, lz));
    m.overlay(&[("plate", r)]);
    let g = mesh_geom_2d(evaluate(&m.list)?, &GeomMeshOpts::new(h))?;
    let side = |sel: &str| g.select(sel);
    let (xmin, xmax, zmin, zmax) = (side("plane(x = 0)")?, side(&format!("plane(x = {lx:e})"))?, side("plane(z = 0)")?, side(&format!("plane(z = {lz:e})"))?);
    let domain = g.domain("plate").ok_or("the plate has no domain")?;
    Ok(Plate { mesh: g.mesh, domain, xmin, xmax, zmin, zmax })
}

/// A box `[0, lx] x [0, ly] x [0, lz]`: the plate swept along y through `ny` layers (`extrude_y`). The vertices of the section are put in coordinate order first (`sort_nodes`), so that opposite side walls can be paired by a Floquet condition. The side walls keep the labels of the
/// plate's sides; the ends are `ymin` (y = 0) and `ymax`.
pub struct Block {
    pub mesh: SimplexMesh,
    pub domain: i32,
    pub xmin: Vec<i32>,
    pub xmax: Vec<i32>,
    pub zmin: Vec<i32>,
    pub zmax: Vec<i32>,
    pub ymin: Vec<i32>,
    pub ymax: Vec<i32>,
}

pub fn block(lx: f64, ly: f64, lz: f64, h: f64, ny: usize) -> Result<Block, String> {
    let p = plate(lx, lz, h)?;
    let nb = p.mesh.facet_labels.iter().copied().max().unwrap_or(0);
    let mesh = crate::meshing::mesher::extrude_y(&sort_nodes(&p.mesh), &[crate::meshing::mesher::Layer { thickness: ly, n: ny, domain_offset: 0, fill: None }], 0.0)?;
    Ok(Block { mesh, domain: p.domain, xmin: p.xmin, xmax: p.xmax, zmin: p.zmin, zmax: p.zmax, ymin: vec![nb + 1], ymax: vec![nb + 2] })
}

/// The mesh with its vertices renumbered in the order of their coordinates (x first, then the second axis, equal to a tolerance of
/// `1e-9` of the extent). Cells and facets are unchanged. Opposite sides of a section then list their vertices in the same order, so a
/// sweep along y (`extrude_y`, whose tetrahedra follow the vertex numbers) cuts opposite side walls the same way and the two can be paired
/// by a Floquet condition.
pub fn sort_nodes(m: &SimplexMesh) -> SimplexMesh {
    let n = m.n_points();
    let d = m.dim;
    let (mut lo, mut hi) = (vec![f64::MAX; d], vec![f64::MIN; d]);
    for p in m.points.chunks(d) {
        for k in 0..d {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let q = |i: usize, k: usize| ((m.points[i * d + k] - lo[k]) / ((hi[k] - lo[k]).max(1e-300) * 1e-9)).round() as i64;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (0..d).map(|k| q(i, k)).collect::<Vec<_>>());
    let mut new_of = vec![0u32; n];
    for (new, &old) in order.iter().enumerate() {
        new_of[old] = new as u32;
    }
    let points: Vec<f64> = order.iter().flat_map(|&i| m.points[i * d..(i + 1) * d].iter().copied()).collect();
    let map = |v: &[u32]| v.iter().map(|&i| new_of[i as usize]).collect::<Vec<u32>>();
    SimplexMesh::new(d, points, map(&m.cells), m.cell_labels.clone(), map(&m.facets), m.facet_labels.clone())
}
