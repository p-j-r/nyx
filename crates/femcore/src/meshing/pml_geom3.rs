//! Structured absorbing layers for a meshed 3D `geom` model (`PmlMode::Structured` in 3D, ADR 52): the 3D counterpart of `pml_geom`.
//!
//! The physical window is meshed alone (the absorbing pieces of the model are left out of the lattice, so the lattice is aligned with
//! the window box and the mesh boundary lies on its faces exactly), optionally adapted, and then `pml_box::build_pml` sweeps prisms on
//! the faces, hexahedra on the edges and on the corners outward, with the number of layers and the thickness the `Pml` feature asked
//! for. The physical vertices and cells keep their numbers (the contract of `build_pml`), so a field computed on the physical mesh can be
//! mapped back by id (feature P3 of the mesher chart).
//!
//! Labels are geom's: the cells of a layer carry the domain number of the absorbing piece they fill (`air.pml.xmax`, edge and corner
//! pieces too); boundary facets are the surviving facets of the physical mesh and the outer wall of the layers, labelled by the geom
//! face they lie on. The native mesh holds boundary facets only; the tetrahedral mesh also has the faces between cells of different domains, with
//! the label of the geom face between those two domains, as in the inline mesh (ADR 40).
//!
//! Needs a window that is a box and a physical mesh whose boundary lies on it (`build_pml` says so when it does not).

use crate::adapt::work::simplex_to_adapt;
use crate::mesh::cells::CellMesh;
use crate::mesh::SimplexMesh;
use crate::meshing::adapt3d::{interface_facets, FaceTable};
use crate::meshing::geomdrive::MeshedGeom;
use crate::meshing::pml_box::{build_pml, PmlSpec};
use crate::meshing::pml_geom::PmlReport;
use geom::model::RegionKind;
use std::collections::{HashMap, HashSet};

/// What `structured_pml_3d_native` returns.
#[derive(Clone)]
pub enum Native3 {
    /// the model has no absorbing layer: the meshed model as it was
    Unchanged(MeshedGeom, PmlReport),
    /// the physical tetrahedra followed by the prisms and hexahedra of the layers, uncut, with the labels described above
    Cells(CellMesh<3>, PmlReport),
}

/// Replaces the absorbing layers of a meshed 3D model by structured ones, cut into tetrahedra for the simplex solvers. A model without a
/// `Pml` comes back unchanged.
pub fn structured_pml_3d(g: MeshedGeom) -> Result<(MeshedGeom, PmlReport), String> {
    let (ev, labels) = (g.eval.clone(), g.labels.clone());
    match structured_pml_3d_native(g)? {
        Native3::Unchanged(g, report) => Ok((g, report)),
        Native3::Cells(cells, report) => {
            let (mesh, _) = SimplexMesh::from_cells_split(&cells)?;
            // the faces between cells of different domains carry the label of the geom face between those two domains, as in the inline
            // mesh (ADR 40); the native mesh holds boundary facets only
            let table = FaceTable::new(&ev, &labels)?;
            let (extra, extra_labels) = interface_facets(&|v| [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]], &mesh.cells, &mesh.cell_labels, &table);
            let (mut facets, mut facet_labels) = (mesh.facets.clone(), mesh.facet_labels.clone());
            for (t, l) in extra.iter().zip(&extra_labels) {
                facets.extend_from_slice(t);
                facet_labels.push(*l);
            }
            // the cut of a prism or hexahedron does not orient its tetrahedra; the lattice mesher does (positive volume), so do the same
            let mut cells = mesh.cells.clone();
            for t in cells.chunks_exact_mut(4) {
                let p = |v: u32| [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]];
                let (a, b, c, d) = (p(t[0]), p(t[1]), p(t[2]), p(t[3]));
                let (u, v, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]], [d[0] - a[0], d[1] - a[1], d[2] - a[2]]);
                if u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]) < 0.0 {
                    t.swap(2, 3);
                }
            }
            let mesh = SimplexMesh::new(3, mesh.points.clone(), cells, mesh.cell_labels.clone(), facets, facet_labels);
            Ok((MeshedGeom { mesh, labels, eval: ev }, report))
        }
    }
}

/// The structured layers without the final cut: a mixed tetrahedron / prism / hexahedron `CellMesh` with triangular and quadrilateral
/// boundary facets. Cells of the absorbing pieces that the mesh still holds are dropped first (the plan leaves them out of the lattice,
/// but any mesh of the model is accepted).
pub fn structured_pml_3d_native(g: MeshedGeom) -> Result<Native3, String> {
    let ev = g.eval.clone();
    let pieces = ev.pml_regions();
    if pieces.is_empty() {
        let n = g.mesh.n_cells();
        return Ok(Native3::Unchanged(g, PmlReport { physical_cells: n, layer_quads: 0, dropped_cells: 0, layers: [[0; 2]; 2], layers3: [[0; 2]; 3] }));
    }
    if g.mesh.dim != 3 {
        return Err("3D structured absorbing layers are built for 3D meshes".into());
    }
    let domain_of_region: Vec<i32> = {
        let mut v = vec![i32::MIN; ev.regions().len()];
        for d in &g.labels.domains {
            v[d.region as usize] = d.id;
        }
        v
    };

    // the layers geom asks for, per side, from the face pieces
    let mut spec = PmlSpec::<3>::uniform(0.0, 0);
    spec.layers = [[0; 2]; 3];
    spec.thickness = [[0.0; 2]; 3];
    for p in pieces {
        let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
        if nz.len() != 1 {
            continue;
        }
        let (a, s) = (nz[0], if p.sides[nz[0]] < 0 { 0 } else { 1 });
        let bb = ev.regions()[p.region as usize].bbox;
        let extent = bb.hi[a] - bb.lo[a];
        if (extent - p.thickness).abs() > 1e-6 * p.thickness {
            return Err(format!(
                "an absorbing piece is {extent:.4e} thick where the Pml asks for {:.4e}: the model ends inside the layer, so the structured layers cannot match it",
                p.thickness
            ));
        }
        let n = (p.layers.round() as usize).max(1);
        if spec.layers[a][s] != 0 && (spec.layers[a][s] != n || (spec.thickness[a][s] - p.thickness).abs() > 1e-9 * p.thickness) {
            return Err(format!("the absorbing pieces on side {s} of axis {a} disagree in thickness or number of layers; the structured layers need one"));
        }
        spec.layers[a][s] = n;
        spec.thickness[a][s] = p.thickness;
    }
    let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;
    spec.label_base = stride;
    spec.label_stride = stride;
    spec.outer_label = None;

    // the piece of every (stretch direction, material it continues)
    let mut piece_of: HashMap<([i8; 3], i32), i32> = HashMap::new();
    for p in pieces {
        if let Some(parent) = p.parent {
            piece_of.insert((p.sides, domain_of_region[parent as usize]), domain_of_region[p.region as usize]);
        }
    }
    let pml_domains: HashSet<i32> = pieces.iter().map(|p| domain_of_region[p.region as usize]).collect();

    // the physical window: the mesh without the cells of the pieces
    let (phys, facet_label) = physical_part(&g.mesh, &pml_domains);
    if phys.n_cells() == 0 {
        return Err("nothing is left of the model outside its absorbing layers".into());
    }
    let dropped_cells = g.mesh.n_cells() - phys.n_cells();
    let work = simplex_to_adapt::<3>(&phys)?;
    let mut built = build_pml(&work, &spec).map_err(|e| format!("structured absorbing layers: {e}"))?;
    let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);

    // domain labels of the layers are geom's
    let mut layer_cells = 0;
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
        layer_cells += 1;
    }

    // boundary facets: the surviving ones of the physical mesh, and the outer wall of the layers by the geom face they lie on
    let table = FaceTable::new(&ev, &g.labels)?;
    let outer_faces: HashMap<i32, Vec<u32>> = {
        // holes and the void have no cells: a face next to one is on the outside
        let mut nothing = vec![true; ev.regions().len()];
        for d in &g.labels.domains {
            nothing[d.region as usize] = matches!(d.kind, RegionKind::Hole | RegionKind::Void);
        }
        let cx = ev.complex();
        let mut m: HashMap<i32, Vec<u32>> = HashMap::new();
        for i in 0..cx.facet_count() as u32 {
            let (l, r) = cx.facet_sides(i);
            let outside = |x: Option<u32>| x.is_none_or(|x| nothing[x as usize]);
            for (a, other) in [(l, r), (r, l)] {
                if let Some(a) = a {
                    if outside(other) && domain_of_region[a as usize] != i32::MIN {
                        m.entry(domain_of_region[a as usize]).or_default().push(i);
                    }
                }
            }
        }
        m
    };
    let face_cells = built.mesh.face_cells();
    let faces = built.mesh.boundary_faces_ordered();
    for f in faces {
        let key = crate::adapt::work::facet_key(&f);
        let owner = face_cells[&key][0] as usize;
        let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
        let label = match physical.then(|| facet_label.get(&key)).flatten() {
            Some(&l) => l,
            None => {
                let pts = f.iter().map(|&v| built.mesh.points[v as usize]).collect::<Vec<_>>();
                let n = pts.len() as f64;
                let mid = [pts.iter().map(|p| p[0]).sum::<f64>() / n, pts.iter().map(|p| p[1]).sum::<f64>() / n, pts.iter().map(|p| p[2]).sum::<f64>() / n];
                nearest_outer_label(&g, &table, &outer_faces, built.mesh.label[owner], mid)?
            }
        };
        built.mesh.push_facet(&f, label);
    }
    let report = PmlReport {
        physical_cells: n_phys_cells,
        layer_quads: layer_cells,
        dropped_cells,
        layers: [spec.layers[0], spec.layers[2]],
        layers3: spec.layers,
    };
    Ok(Native3::Cells(built.mesh, report))
}

/// The label of the geom face, among the outer faces of the region `domain`, that is closest to `p`.
fn nearest_outer_label(g: &MeshedGeom, table: &FaceTable, outer_faces: &HashMap<i32, Vec<u32>>, domain: i32, p: [f64; 3]) -> Result<i32, String> {
    let cands = outer_faces.get(&domain).ok_or_else(|| {
        format!("region '{}' has no face on the outside of the model, but the mesh has a boundary facet on it", g.labels.name_of_domain(domain).unwrap_or("?"))
    })?;
    table.nearest_of(p, cands.iter().map(|&i| i as usize)).ok_or_else(|| "no candidate face".to_string())
}

/// The mesh without the cells of the given domains, vertices renumbered in order (unchanged when nothing is dropped), and the labels of
/// the facets it keeps by face key.
fn physical_part(m: &SimplexMesh, drop: &HashSet<i32>) -> (SimplexMesh, HashMap<crate::adapt::work::FacetKey, i32>) {
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
    let mut faces: HashSet<[u32; 3]> = HashSet::new();
    for &c in &keep {
        let v: Vec<u32> = m.cell(c).iter().map(|&x| new_id[x as usize]).collect();
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            f.sort_unstable();
            faces.insert([f[0], f[1], f[2]]);
        }
        cells.extend_from_slice(&v);
        cell_labels.push(m.cell_labels[c]);
    }
    let mut facets = vec![];
    let mut facet_labels = vec![];
    let mut by_key = HashMap::new();
    for f in 0..m.n_facets() {
        let v = m.facet(f);
        let t = [new_id[v[0] as usize], new_id[v[1] as usize], new_id[v[2] as usize]];
        let mut k = t;
        k.sort_unstable();
        if t.contains(&u32::MAX) || !faces.contains(&k) {
            continue;
        }
        facets.extend_from_slice(&t);
        facet_labels.push(m.facet_labels[f]);
        by_key.insert(crate::adapt::work::facet_key(&t), m.facet_labels[f]);
    }
    (SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels), by_key)
}
