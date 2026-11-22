//! Structured absorbing layers for a meshed `geom` model (`PmlMode::Structured`).
//!
//! The default way a model gets its absorbing layers is the one geom describes: the `Pml` feature cuts pieces off the regions outside a
//! window, and the mesher fills them with triangles like any other region (`PmlMode::Inline`). This module is the alternative with the
//! rules of the ported `build_pml`: the model is meshed, the triangles of the pieces are dropped, and the layers are rebuilt as
//! quadrilaterals that continue the boundary of the physical window node for node, with exactly the number of layers and the
//! thickness the `Pml` feature asked for, then cut into two triangles each for the solvers.
//!
//! Nothing the solver sees changes: the domain labels of the new cells are the labels geom already gave the pieces (`sapphire.pml.xmax`,
//! the corner pieces too), so `domains("sapphire*")` and `pml_materials` work as before, and the facets on the outer wall of the layers
//! carry the label of the geom edge they lie on. The facets of the physical window that face a layer are gone, because they are interior
//! now (as in the inline mesh, where an interface has no facet).
//!
//! Needs a window that is a box and a physical mesh whose boundary lies on it (`build_pml` says so when it does not).

use crate::meshing::geomdrive::MeshedGeom;
use crate::mesh::cells::CellMesh;
use crate::mesh::SimplexMesh;
use crate::meshing::pml_box::{build_pml, PmlSpec};
use crate::adapt::work::simplex_to_adapt;
use geom::model::RegionKind;
use std::collections::{HashMap, HashSet};

/// How the absorbing layers of a model are meshed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PmlMode {
    /// Triangles inside the model, made by the same mesher as everything else (what the app has always done).
    #[default]
    Inline,
    /// Quadrilateral layers from `build_pml` around the physical window, split into triangles.
    Structured,
}

/// What `structured_pml` built.
#[derive(Clone, Debug, PartialEq)]
pub struct PmlReport {
    /// cells of the physical window
    pub physical_cells: usize,
    /// quadrilaterals of the layers (each is two cells of the final mesh)
    pub layer_quads: usize,
    /// cells the inline layers had before they were replaced
    pub dropped_cells: usize,
    /// layers per side, `[axis][side]` for the mesh axes x and z (0 = no layer)
    pub layers: [[usize; 2]; 2],
    /// layers per side, `[axis][side]` for the axes x, y, z of a 3D model (all zero for a 2D one)
    pub layers3: [[usize; 2]; 3],
}

/// Replaces the absorbing layers of a meshed 2D model by structured ones, cut into triangles for the simplex solvers. A model without a
/// `Pml` comes back unchanged.
pub fn structured_pml(g: MeshedGeom) -> Result<(MeshedGeom, PmlReport), String> {
    if g.mesh.dim == 3 {
        return crate::meshing::pml_geom3::structured_pml_3d(g);
    }
    let ev = g.eval.clone();
    let labels = g.labels.clone();
    let warnings = g.warnings.clone();
    match structured_pml_native(g)? {
        Native::Unchanged(g, report) => Ok((g, report)),
        Native::Cells(cells, report) => {
            let (mesh, _) = SimplexMesh::from_cells_split(&cells)?;
            Ok((MeshedGeom { mesh, labels, eval: ev, warnings }, report))
        }
    }
}

/// What `structured_pml_native` returns.
#[derive(Clone)]
pub enum Native {
    /// the model has no absorbing layer: the meshed model as it was
    Unchanged(MeshedGeom, PmlReport),
    /// the physical triangles and the quadrilateral layers, uncut, with the labels of `structured_pml` (boundary facets are the 1-based
    /// edges of the lower dimension)
    Cells(CellMesh<2>, PmlReport),
}

/// The structured layers without the final cut: a mixed triangle/quadrilateral `CellMesh`, for solvers that take quadrilaterals natively.
/// Same labels, same facets, same requirements as `structured_pml`, whose result is this mesh with every quadrilateral split.
pub fn structured_pml_native(g: MeshedGeom) -> Result<Native, String> {
    let ev = g.eval.clone();
    let pieces = ev.pml_regions();
    if pieces.is_empty() {
        let n = g.mesh.n_cells();
        return Ok(Native::Unchanged(g, PmlReport { physical_cells: n, layer_quads: 0, dropped_cells: 0, layers: [[0; 2]; 2], layers3: [[0; 2]; 3] }));
    }
    if g.mesh.dim != 2 {
        return Err("structured absorbing layers are only built for 2D meshes".into());
    }
    let domain_of_region: Vec<i32> = {
        let mut v = vec![i32::MIN; ev.regions().len()];
        for d in &g.labels.domains {
            v[d.region as usize] = d.id;
        }
        v
    };

    // the layers geom asks for, per side, from the face pieces; mesh axes (x, z) are geom axes (0, 2)
    const GEOM_AXIS: [usize; 2] = [0, 2];
    let mut spec = PmlSpec::<2>::uniform(0.0, 0);
    spec.layers = [[0; 2]; 2];
    spec.thickness = [[0.0; 2]; 2];
    for p in pieces {
        let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
        if nz.len() != 1 {
            continue;
        }
        let Some(a) = GEOM_AXIS.iter().position(|&ga| ga == nz[0]) else {
            return Err("an absorbing layer along y cannot be built in a 2D section".into());
        };
        let s = if p.sides[nz[0]] < 0 { 0 } else { 1 };
        let bb = ev.regions()[p.region as usize].bbox;
        let extent = bb.hi[nz[0]] - bb.lo[nz[0]];
        if (extent - p.thickness).abs() > 1e-6 * p.thickness {
            return Err(format!(
                "an absorbing piece is {extent:.4e} thick where the Pml asks for {:.4e}: the model ends inside the layer, so the structured layers cannot match it",
                p.thickness
            ));
        }
        let n = (p.layers.round() as usize).max(1);
        if spec.layers[a][s] != 0 && (spec.layers[a][s] != n || (spec.thickness[a][s] - p.thickness).abs() > 1e-9 * p.thickness) {
            return Err(format!(
                "the absorbing pieces on side {s} of axis {a} disagree in thickness or number of layers; the structured layers need one"
            ));
        }
        spec.layers[a][s] = n;
        spec.thickness[a][s] = p.thickness;
    }
    let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;
    spec.label_base = stride;
    spec.label_stride = stride;
    spec.outer_label = None;

    // piece of every (stretch direction, material it continues)
    let mut piece_of: HashMap<([i8; 2], i32), i32> = HashMap::new();
    for p in pieces {
        if let Some(parent) = p.parent {
            piece_of.insert(([p.sides[0], p.sides[2]], domain_of_region[parent as usize]), domain_of_region[p.region as usize]);
        }
    }
    let pml_domains: HashSet<i32> = pieces.iter().map(|p| domain_of_region[p.region as usize]).collect();

    // the physical window: the mesh without the cells of the pieces
    let (phys, facet_label) = physical_part(&g.mesh, &pml_domains);
    if phys.n_cells() == 0 {
        return Err("nothing is left of the model outside its absorbing layers".into());
    }
    let dropped_cells = g.mesh.n_cells() - phys.n_cells();
    let work = simplex_to_adapt::<2>(&phys)?;
    let mut built = build_pml(&work, &spec).map_err(|e| format!("structured absorbing layers: {e}"))?;
    let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);

    // domain labels of the layers are geom's
    let mut layer_quads = 0;
    for c in n_phys_cells..built.mesh.n_cells() {
        let (code, material) = built
            .spec
            .decode(built.mesh.label[c])
            .ok_or_else(|| format!("layer cell {c} carries the physical label {}", built.mesh.label[c]))?;
        let id = piece_of.get(&(code, material)).copied().ok_or_else(|| {
            format!(
                "geom has no absorbing piece of {} for the stretch {:?}; the physical cell next to that layer is a different region than the piece was cut from",
                g.labels.name_of_domain(material).unwrap_or("?"),
                code
            )
        })?;
        built.mesh.label[c] = id;
        layer_quads += 1;
    }

    // boundary facets of the whole mesh: the ones the physical mesh already had, and the outer wall of the layers by the geom edge they lie on
    let outer_edges: HashMap<i32, Vec<u32>> = {
        // holes and the void have no cells: an edge next to one is on the outside
        let mut nothing = vec![true; ev.regions().len()];
        for d in &g.labels.domains {
            nothing[d.region as usize] = matches!(d.kind, RegionKind::Hole | RegionKind::Void);
        }
        let mut m: HashMap<i32, Vec<u32>> = HashMap::new();
        for (i, e) in ev.edges().iter().enumerate() {
            let outside = |r: Option<u32>| r.is_none_or(|r| nothing[r as usize]);
            for (r, other) in [(e.left, e.right), (e.right, e.left)] {
                if let Some(r) = r {
                    if outside(other) && domain_of_region[r as usize] != i32::MIN {
                        m.entry(domain_of_region[r as usize]).or_default().push(i as u32);
                    }
                }
            }
        }
        m
    };
    let owners = {
        let fc = built.mesh.face_cells();
        let mut v: HashMap<[u32; 4], u32> = HashMap::new();
        for (k, cs) in fc {
            if cs.len() == 1 {
                v.insert(k, cs[0]);
            }
        }
        v
    };
    let faces = built.mesh.boundary_faces_ordered();
    for f in faces {
        let key = crate::adapt::work::facet_key(&f);
        let owner = owners[&key] as usize;
        let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
        let label = match physical.then(|| facet_label.get(&sorted_pair(f[0], f[1]))).flatten() {
            Some(&l) => l,
            None => {
                let mid = [
                    0.5 * (built.mesh.points[f[0] as usize][0] + built.mesh.points[f[1] as usize][0]),
                    0.5 * (built.mesh.points[f[0] as usize][1] + built.mesh.points[f[1] as usize][1]),
                ];
                nearest_outer_label(&ev, &g, &outer_edges, built.mesh.label[owner], mid)?
            }
        };
        built.mesh.push_facet(&f, label);
    }

    let report = PmlReport { physical_cells: n_phys_cells, layer_quads, dropped_cells, layers: spec.layers, layers3: [[0; 2]; 3] };
    Ok(Native::Cells(built.mesh, report))
}

fn sorted_pair(a: u32, b: u32) -> (u32, u32) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// The label of the geom edge, among the outer edges of the region `domain`, that is closest to `p`.
fn nearest_outer_label(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    outer_edges: &HashMap<i32, Vec<u32>>,
    domain: i32,
    p: [f64; 2],
) -> Result<i32, String> {
    let cands = outer_edges.get(&domain).ok_or_else(|| {
        format!("region '{}' has no edge on the outside of the model, but the mesh has a boundary facet on it", g.labels.name_of_domain(domain).unwrap_or("?"))
    })?;
    let mut best = (f64::INFINITY, None);
    for &e in cands {
        let (q, _) = ev.edge_closest(e, [p[0], 0.0, p[1]]);
        let d = (q[0] - p[0]).hypot(q[2] - p[1]);
        if d < best.0 {
            best = (d, Some(e));
        }
    }
    let e = best.1.ok_or("no candidate edge")?;
    g.labels.label_of_edge(e).ok_or_else(|| format!("edge {e} has no boundary label"))
}

/// The mesh without the cells of the given domains, vertices renumbered in order, and the labels of the facets it keeps by vertex pair.
fn physical_part(m: &SimplexMesh, drop: &HashSet<i32>) -> (SimplexMesh, HashMap<(u32, u32), i32>) {
    let keep: Vec<usize> = (0..m.n_cells()).filter(|&c| !drop.contains(&m.cell_labels[c])).collect();
    let mut used = vec![false; m.n_points()];
    for &c in &keep {
        for &v in m.cell(c) {
            used[v as usize] = true;
        }
    }
    let mut new_id = vec![u32::MAX; m.n_points()];
    let mut n = 0u32;
    let mut points = vec![];
    for v in 0..m.n_points() {
        if used[v] {
            new_id[v] = n;
            n += 1;
            points.extend_from_slice(m.point(v));
        }
    }
    let mut cells = vec![];
    let mut cell_labels = vec![];
    let mut edges: HashSet<(u32, u32)> = HashSet::new();
    for &c in &keep {
        let v: Vec<u32> = m.cell(c).iter().map(|&x| new_id[x as usize]).collect();
        for k in 0..3 {
            edges.insert(sorted_pair(v[k], v[(k + 1) % 3]));
        }
        cells.extend_from_slice(&v);
        cell_labels.push(m.cell_labels[c]);
    }
    let mut facets = vec![];
    let mut facet_labels = vec![];
    let mut by_pair = HashMap::new();
    for f in 0..m.n_facets() {
        let v = m.facet(f);
        let (a, b) = (new_id[v[0] as usize], new_id[v[1] as usize]);
        if a == u32::MAX || b == u32::MAX || !edges.contains(&sorted_pair(a, b)) {
            continue;
        }
        facets.extend_from_slice(&[a, b]);
        facet_labels.push(m.facet_labels[f]);
        by_pair.insert(sorted_pair(a, b), m.facet_labels[f]);
    }
    (SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels), by_pair)
}
