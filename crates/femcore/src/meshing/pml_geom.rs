// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified 2D structured absorbing-layer dispatcher (`PmlMode::{Auto, Structured, Inline}`) and
//! shared mixed-cell PML scaled-Jacobian quality evaluator across Cartesian, Cylindrical, and
//! Axisymmetric Spherical topologies (both `Outward` and legacy `Inward` modes).

use crate::adapt::work::{facet_key, simplex_to_adapt};
use crate::fem_mixed::{quad_rule_mixed, IsoparametricCell, LagrangeSpace};
use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::mixed::MixedDiscretisation;
use crate::mesh::{Mesh, SimplexMesh};
use crate::meshing::geomdrive::MeshedGeom;
use crate::meshing::pml_box::{build_pml, PmlSpec};
use crate::meshing::pml_cyl::{build_pml_annulus_2d, CylPmlSpec2d};
use crate::meshing::pml_sphere::{build_pml_semicircle_2d, SpherePmlSpec2d};
use geom::model::{PmlCoordKind, PmlGrowMode, PmlRole, RegionKind};
use std::collections::{HashMap, HashSet};

/// How the absorbing layers of a model are meshed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PmlMode {
    /// Physics-controlled default: builds structured native `Quad`/`Prism`/`Hex` layers whenever the model
    /// defines absorbing layers (`!ev.pml_regions().is_empty()`), and does nothing otherwise.
    #[default]
    Auto,
    /// Structured mapped/swept absorbing layers (`Quad` in 2D; `Prism` + `Hex` in 3D, zero pyramids).
    Structured,
    /// Legacy unstructured inline simplices inside CAD-sliced absorbing domains.
    Inline,
}

impl PmlMode {
    /// True when structured absorbing-layer generation is active for a model with `has_pml_regions`.
    #[inline]
    pub fn is_structured(self, has_pml_regions: bool) -> bool {
        match self {
            PmlMode::Auto => has_pml_regions,
            PmlMode::Structured => has_pml_regions,
            PmlMode::Inline => false,
        }
    }
}

/// Summary of structured absorbing-layer generation.
#[derive(Clone, Debug, PartialEq)]
pub struct PmlReport {
    /// Number of physical interior cells (`0..physical_cells`).
    pub physical_cells: usize,
    /// Number of structured absorbing-layer cells (`Quad` in 2D; `Prism` + `Hex` in 3D).
    pub layer_quads: usize,
    /// Inline absorbing-layer cells dropped from the input mesh (`0` when skipped upfront or `Outward`).
    pub dropped_cells: usize,
    /// Layers per side `[axis][side]` for 2D axes `(x, z)` (`0` = inactive).
    pub layers: [[usize; 2]; 2],
    /// Layers per side `[axis][side]` for 3D axes `(x, y, z)` (all zero in 2D).
    pub layers3: [[usize; 2]; 3],
    /// Minimum normalized scaled Jacobian across all generated PML cells (curved and straight).
    pub min_pml_scaled_jacobian: Option<f64>,
}

impl PmlReport {
    pub fn empty(physical_cells: usize) -> Self {
        Self {
            physical_cells,
            layer_quads: 0,
            dropped_cells: 0,
            layers: [[0; 2]; 2],
            layers3: [[0; 2]; 3],
            min_pml_scaled_jacobian: None,
        }
    }
}

/// Evaluates the minimum normalized scaled Jacobian in `[-1, 1]` of cell `ci` in `md` (accounting for
/// curved isoparametric `edge_nodes` and `quad_face_nodes` at all Lagrange nodes and degree-4 Gauss points).
pub fn mixed_cell_scaled_jacobian_min<const D: usize>(
    cm: &CellMesh<D>,
    md: &MixedDiscretisation,
    ci: usize,
) -> f64 {
    if !md.cell_is_curved(ci) {
        return crate::mesh::quality::cell_scaled_jacobian_min(cm, ci);
    }
    let iso = IsoparametricCell::from_mesh(md, ci);
    let ct = iso.cell_type;
    let mut pts = crate::fem_mixed::lagrange_ref_nodes(ct, LagrangeSpace::Order2Tensor);
    let qrule = quad_rule_mixed(ct, 4);
    pts.extend_from_slice(&qrule.pts);

    let mut scratch_n = Vec::with_capacity(27);
    let mut scratch_dn = Vec::with_capacity(27);
    let inv_sqrt3 = 1.0 / 3.0_f64.sqrt();
    let mut min_q = f64::INFINITY;

    for xi in &pts {
        let je = iso.eval_into(xi, &mut scratch_n, &mut scratch_dn);
        if !je.det.is_finite() {
            return -1.0;
        }
        let q = match ct {
            CellType::Quad => {
                let n0 = je.j[0][0].hypot(je.j[1][0]);
                let n1 = je.j[0][1].hypot(je.j[1][1]);
                let denom = n0 * n1;
                if denom > 1e-300 {
                    je.det / denom
                } else {
                    0.0
                }
            }
            CellType::Hex => {
                let col_norm = |k: usize| {
                    (je.j[0][k] * je.j[0][k] + je.j[1][k] * je.j[1][k] + je.j[2][k] * je.j[2][k]).sqrt()
                };
                let denom = col_norm(0) * col_norm(1) * col_norm(2);
                if denom > 1e-300 {
                    je.det / denom
                } else {
                    0.0
                }
            }
            CellType::Prism => {
                // Map reference right-triangle columns `(j_0, j_1)` to equilateral orthonormal basis `(t_0, t_1)`:
                // `t_0 = j_0`, `t_1 = (2 * j_1 - j_0) / sqrt(3)`.
                let j0 = [je.j[0][0], je.j[1][0], je.j[2][0]];
                let j1 = [je.j[0][1], je.j[1][1], je.j[2][1]];
                let j2 = [je.j[0][2], je.j[1][2], je.j[2][2]];
                let t1 = [
                    (2.0 * j1[0] - j0[0]) * inv_sqrt3,
                    (2.0 * j1[1] - j0[1]) * inv_sqrt3,
                    (2.0 * j1[2] - j0[2]) * inv_sqrt3,
                ];
                let sum_t2 = (j0[0] * j0[0] + j0[1] * j0[1] + j0[2] * j0[2])
                    + (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
                let n2 = (j2[0] * j2[0] + j2[1] * j2[1] + j2[2] * j2[2]).sqrt();
                let denom = sum_t2 * n2;
                if denom > 1e-300 {
                    2.0 * je.det / denom
                } else {
                    0.0
                }
            }
            CellType::Tri | CellType::Tet => {
                let base_q = crate::mesh::quality::cell_scaled_jacobian_min(cm, ci).max(0.0);
                let straight_meas = cm.measure(ci);
                let ref_meas = match ct {
                    CellType::Tri => 0.5,
                    CellType::Tet => 1.0 / 6.0,
                    _ => 1.0,
                };
                let det_s = straight_meas / ref_meas;
                if det_s.abs() > 1e-300 {
                    (je.det / det_s) * base_q
                } else {
                    0.0
                }
            }
        };
        min_q = min_q.min(q);
    }
    if min_q.is_finite() {
        min_q
    } else {
        0.0
    }
}

/// Computes the minimum normalized scaled Jacobian across all PML cells `pml_start_cell..cm.n_cells()`.
pub fn compute_min_pml_scaled_jacobian<const D: usize>(
    cm: &CellMesh<D>,
    pml_start_cell: usize,
) -> Option<f64> {
    if pml_start_cell >= cm.n_cells() {
        return None;
    }
    let md = MixedDiscretisation::from_cell_mesh(cm).ok()?;
    let min_sj = (pml_start_cell..cm.n_cells())
        .map(|ci| mixed_cell_scaled_jacobian_min(cm, &md, ci))
        .fold(f64::INFINITY, f64::min);
    min_sj.is_finite().then_some(min_sj)
}

/// Replaces or extrudes the absorbing layers of a meshed model (`2D` or `3D`), keeping the uncut
/// `CellMesh` in `MeshedGeom::cells` and splitting `MeshedGeom::mesh` for legacy simplex callers.
pub fn structured_pml(g: MeshedGeom) -> Result<(MeshedGeom, PmlReport), String> {
    structured_pml_with_curved(g, false)
}

/// Like [`structured_pml`], with explicit `curved` control for concentric curvilinear shells.
pub fn structured_pml_with_curved(g: MeshedGeom, curved: bool) -> Result<(MeshedGeom, PmlReport), String> {
    if g.mesh.dim == 3 {
        return crate::meshing::pml_geom3::structured_pml_3d_with_curved(g, curved);
    }
    let ev = g.eval.clone();
    let labels = g.labels.clone();
    let warnings = g.warnings.clone();
    match structured_pml_native_with_curved(g, curved)? {
        Native::Unchanged(g, report) => Ok((g, report)),
        Native::Cells(cells, report) => {
            let (mesh, _) = SimplexMesh::from_cells_split(&cells)?;
            mesh.check()?;
            Ok((
                MeshedGeom {
                    mesh,
                    cells: Some(Mesh::D2(cells)),
                    labels,
                    eval: ev,
                    warnings,
                },
                report,
            ))
        }
    }
}

/// What [`structured_pml_native`] returns.
#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
pub enum Native {
    /// The model has no absorbing layer: the meshed model as it was.
    Unchanged(MeshedGeom, PmlReport),
    /// The physical triangles followed by structured `Quad` absorbing-layer cells, uncut.
    Cells(CellMesh<2>, PmlReport),
}

/// Builds uncut 2D mixed `Tri` + `Quad` structured absorbing layers (`curved = false`).
pub fn structured_pml_native(g: MeshedGeom) -> Result<Native, String> {
    structured_pml_native_with_curved(g, false)
}

/// Builds uncut 2D mixed `Tri` + `Quad` structured absorbing layers across Cartesian, Cylindrical
/// Annulus, and Axisymmetric Spherical Semicircle topologies (both `Outward` and `Inward`).
pub fn structured_pml_native_with_curved(g: MeshedGeom, curved: bool) -> Result<Native, String> {
    let ev = g.eval.clone();
    ev.validate_pml_ray_invariance().map_err(|e| e.to_string())?;
    let pieces = ev.pml_regions();
    if pieces.is_empty() {
        let n = g.mesh.n_cells();
        return Ok(Native::Unchanged(g, PmlReport::empty(n)));
    }
    if g.mesh.dim != 2 {
        return Err("2D structured absorbing layers require a 2D mesh".into());
    }

    let domain_of_region: Vec<i32> = {
        let mut v = vec![i32::MIN; ev.regions().len()];
        for d in &g.labels.domains {
            v[d.region as usize] = d.id;
        }
        v
    };
    let pml_domains: HashSet<i32> = pieces.iter().map(|p| domain_of_region[p.region as usize]).collect();
    let (phys, facet_label, phys_edge_nodes) = physical_part(&g.mesh, &pml_domains);
    if phys.n_cells() == 0 {
        return Err("nothing is left of the model outside its absorbing layers".into());
    }
    let dropped_cells = g.mesh.n_cells() - phys.n_cells();
    let work = simplex_to_adapt::<2>(&phys)?;

    let outer_edges = collect_outer_geom_edges(&ev, &g, &domain_of_region);
    let mut parent_of_pml: HashMap<i32, i32> = HashMap::new();
    for p in pieces {
        if let Some(par) = p.parent {
            parent_of_pml.insert(domain_of_region[p.region as usize], domain_of_region[par as usize]);
        }
    }

    let merge_phys_edges_2d = |cm: &mut CellMesh<2>, phys_nodes: &[([u32; 2], [f64; 2])]| {
        if phys_nodes.is_empty() {
            return;
        }
        let mut map: std::collections::BTreeMap<[u32; 2], [f64; 2]> =
            phys_nodes.iter().map(|&([a, b], p)| ([a.min(b), a.max(b)], p)).collect();
        for &(pr, p) in &cm.edge_nodes {
            map.insert([pr[0].min(pr[1]), pr[0].max(pr[1])], p);
        }
        cm.set_edge_nodes(map);
    };

    let is_spherical_2d = pieces.iter().any(|p| p.coord == PmlCoordKind::Spherical);
    let is_cyl_annulus_2d = !is_spherical_2d
        && pieces
            .iter()
            .any(|p| p.coord == PmlCoordKind::Cylindrical && p.role == PmlRole::RadialShell && p.axis[1].abs() > 0.9 && p.r_in > 0.0);

    // ------------------------------------------------------------------ 1. 2D Axisymmetric Spherical Semicircle
    if is_spherical_2d {
        let p0 = &pieces[0];
        let n_layers = (p0.layers.round() as usize).max(1);
        let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;
        let mut piece_of_mat: HashMap<i32, i32> = HashMap::new();
        for p in pieces {
            if let Some(par) = p.parent {
                piece_of_mat.insert(domain_of_region[par as usize], domain_of_region[p.region as usize]);
            }
        }
        let mut sp = SpherePmlSpec2d::new(p0.center[2], p0.r_in, p0.thickness, n_layers);
        sp.grading = p0.layer_grading.max(1e-6);
        sp.curved = curved;
        sp.label_base = stride;
        sp.label_stride = stride;
        sp.tol = 5e-3;

        let mut built = build_pml_semicircle_2d(&work, &sp).map_err(|e| format!("structured spherical absorbing layers: {e}"))?;
        merge_phys_edges_2d(&mut built.mesh, &phys_edge_nodes);
        let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);
        let mut layer_quads = 0usize;
        for c in n_phys_cells..built.mesh.n_cells() {
            let mat = built
                .spec
                .decode(built.mesh.label[c])
                .ok_or_else(|| format!("layer cell {c} carries physical label {}", built.mesh.label[c]))?;
            let id = piece_of_mat.get(&mat).copied().ok_or_else(|| {
                format!(
                    "geom has no spherical absorbing piece for region '{}'",
                    g.labels.name_of_domain(mat).unwrap_or("?")
                )
            })?;
            built.mesh.label[c] = id;
            layer_quads += 1;
        }

        // In Inward mode, remap outer semicircle (`R = R_out`) facets to the outer geom edge of `pml_dom`,
        // while keeping all `r = 0` axial facets on the symmetry axis label (never outer PEC!).
        if p0.mode == PmlGrowMode::Inward {
            let tol_axis = (1e-6 * p0.r_in).max(1e-12);
            // Determine physical `r = 0` symmetry labels near north/south poles first.
            let mut north_sym: Option<i32> = None;
            let mut south_sym: Option<i32> = None;
            for (&(a, b), &lab) in &facet_label {
                let pa = built.mesh.points[a as usize];
                let pb = built.mesh.points[b as usize];
                if pa[0].abs() <= tol_axis && pb[0].abs() <= tol_axis {
                    let zm = 0.5 * (pa[1] + pb[1]);
                    if zm >= p0.center[2] {
                        north_sym = Some(lab);
                    } else {
                        south_sym = Some(lab);
                    }
                }
            }
            built.mesh.facet_offsets = vec![0];
            built.mesh.facet_conn.clear();
            built.mesh.facet_label.clear();
            let owners = single_cell_owners(&built.mesh);
            for f in built.mesh.boundary_faces_ordered() {
                let key = facet_key(&f);
                let owner = owners[&key] as usize;
                let pa = built.mesh.points[f[0] as usize];
                let pb = built.mesh.points[f[1] as usize];
                let on_axis = pa[0].abs() <= tol_axis && pb[0].abs() <= tol_axis;
                let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
                let lab = if physical {
                    if let Some(&l) = facet_label.get(&sorted_pair(f[0], f[1])) {
                        l
                    } else {
                        nearest_outer_label(&ev, &g, &outer_edges, built.mesh.label[owner], [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])], None)?
                    }
                } else if on_axis {
                    let zm = 0.5 * (pa[1] + pb[1]);
                    if zm >= p0.center[2] {
                        north_sym.or(south_sym).unwrap_or(1)
                    } else {
                        south_sym.or(north_sym).unwrap_or(1)
                    }
                } else {
                    nearest_outer_label(
                        &ev,
                        &g,
                        &outer_edges,
                        built.mesh.label[owner],
                        [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])],
                        Some(outward_edge_normal_2d(&built.mesh, owner, [f[0], f[1]])),
                    )?
                };
                built.mesh.push_facet(&f, lab);
            }
        }

        let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
        let report = PmlReport {
            physical_cells: n_phys_cells,
            layer_quads,
            dropped_cells,
            layers: [[0, n_layers], [n_layers, n_layers]],
            layers3: [[0; 2]; 3],
            min_pml_scaled_jacobian,
        };
        return Ok(Native::Cells(built.mesh, report));
    }

    // ------------------------------------------------------------------ 2. 2D Planar Cylindrical Annulus
    if is_cyl_annulus_2d {
        let p0 = &pieces[0];
        let n_layers = (p0.layers.round() as usize).max(1);
        let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;
        let mut piece_of_mat: HashMap<i32, i32> = HashMap::new();
        for p in pieces {
            if let Some(par) = p.parent {
                piece_of_mat.insert(domain_of_region[par as usize], domain_of_region[p.region as usize]);
            }
        }
        let mut cs = CylPmlSpec2d::new([p0.center[0], p0.center[2]], p0.r_in, p0.thickness, n_layers);
        cs.grading = p0.layer_grading.max(1e-6);
        cs.curved = curved;
        cs.label_base = stride;
        cs.label_stride = stride;
        cs.tol = 5e-3;

        let mut built = build_pml_annulus_2d(&work, &cs).map_err(|e| format!("structured cylindrical absorbing layers: {e}"))?;
        merge_phys_edges_2d(&mut built.mesh, &phys_edge_nodes);
        let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);
        let mut layer_quads = 0usize;
        for c in n_phys_cells..built.mesh.n_cells() {
            let mat = built
                .spec
                .decode(built.mesh.label[c])
                .ok_or_else(|| format!("layer cell {c} carries physical label {}", built.mesh.label[c]))?;
            let id = piece_of_mat.get(&mat).copied().ok_or_else(|| {
                format!(
                    "geom has no cylindrical absorbing piece for region '{}'",
                    g.labels.name_of_domain(mat).unwrap_or("?")
                )
            })?;
            built.mesh.label[c] = id;
            layer_quads += 1;
        }

        if p0.mode == PmlGrowMode::Inward {
            built.mesh.facet_offsets = vec![0];
            built.mesh.facet_conn.clear();
            built.mesh.facet_label.clear();
            let owners = single_cell_owners(&built.mesh);
            for f in built.mesh.boundary_faces_ordered() {
                let key = facet_key(&f);
                let owner = owners[&key] as usize;
                let pa = built.mesh.points[f[0] as usize];
                let pb = built.mesh.points[f[1] as usize];
                let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
                let lab = match physical.then(|| facet_label.get(&sorted_pair(f[0], f[1]))).flatten() {
                    Some(&l) => l,
                    None => nearest_outer_label(
                        &ev,
                        &g,
                        &outer_edges,
                        built.mesh.label[owner],
                        [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])],
                        Some(outward_edge_normal_2d(&built.mesh, owner, [f[0], f[1]])),
                    )?,
                };
                built.mesh.push_facet(&f, lab);
            }
        }

        let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
        let report = PmlReport {
            physical_cells: n_phys_cells,
            layer_quads,
            dropped_cells,
            layers: [[n_layers, n_layers], [n_layers, n_layers]],
            layers3: [[0; 2]; 3],
            min_pml_scaled_jacobian,
        };
        return Ok(Native::Cells(built.mesh, report));
    }

    // ------------------------------------------------------------------ 3. 2D Cartesian & 2D Axisymmetric Cylindrical
    const GEOM_AXIS: [usize; 2] = [0, 2];
    let mut spec = PmlSpec::<2>::uniform(0.0, 0);
    spec.layers = [[0; 2]; 2];
    spec.thickness = [[0.0; 2]; 2];
    spec.grading = pieces[0].layer_grading.max(1e-6);

    for p in pieces {
        let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
        if nz.len() != 1 {
            continue;
        }
        let Some(a) = GEOM_AXIS.iter().position(|&ga| ga == nz[0]) else {
            return Err("an absorbing layer along y cannot be built in a 2D section".into());
        };
        let s = if p.sides[nz[0]] < 0 { 0 } else { 1 };
        if p.mode == PmlGrowMode::Inward {
            let bb = ev.regions()[p.region as usize].bbox;
            let extent = bb.hi[nz[0]] - bb.lo[nz[0]];
            if (extent - p.thickness).abs() > 1e-6 * p.thickness {
                return Err(format!(
                    "an absorbing piece is {extent:.4e} thick where the Pml asks for {:.4e}: the model ends inside the layer, so the structured layers cannot match it",
                    p.thickness
                ));
            }
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

    let phys_bb = work.bbox();
    let win = pieces[0].window;
    spec.interface_box = Some(std::array::from_fn(|a| {
        let ga = GEOM_AXIS[a];
        [
            if spec.layers[a][0] > 0 { win.lo[ga] } else { phys_bb[a][0] },
            if spec.layers[a][1] > 0 { win.hi[ga] } else { phys_bb[a][1] },
        ]
    }));

    let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;
    spec.label_base = stride;
    spec.label_stride = stride;
    spec.outer_label = None;

    let mut piece_of: HashMap<([i8; 2], i32), i32> = HashMap::new();
    for p in pieces {
        if let Some(parent) = p.parent {
            piece_of.insert(
                ([p.sides[0], p.sides[2]], domain_of_region[parent as usize]),
                domain_of_region[p.region as usize],
            );
        }
    }

    let mut built = build_pml(&work, &spec).map_err(|e| format!("structured absorbing layers: {e}"))?;
    merge_phys_edges_2d(&mut built.mesh, &phys_edge_nodes);
    let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);

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

    built.mesh.facet_offsets = vec![0];
    built.mesh.facet_conn.clear();
    built.mesh.facet_label.clear();
    let owners = single_cell_owners(&built.mesh);
    let faces = built.mesh.boundary_faces_ordered();
    for f in faces {
        let key = facet_key(&f);
        let owner = owners[&key] as usize;
        let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
        let label = match physical.then(|| facet_label.get(&sorted_pair(f[0], f[1]))).flatten() {
            Some(&l) => l,
            None => {
                let mid = [
                    0.5 * (built.mesh.points[f[0] as usize][0] + built.mesh.points[f[1] as usize][0]),
                    0.5 * (built.mesh.points[f[0] as usize][1] + built.mesh.points[f[1] as usize][1]),
                ];
                let owner_dom = built.mesh.label[owner];
                let lookup_dom = if outer_edges.contains_key(&owner_dom) {
                    owner_dom
                } else {
                    parent_of_pml.get(&owner_dom).copied().unwrap_or(owner_dom)
                };
                let n_out = outward_edge_normal_2d(&built.mesh, owner, [f[0], f[1]]);
                nearest_outer_label(&ev, &g, &outer_edges, lookup_dom, mid, Some(n_out))?
            }
        };
        built.mesh.push_facet(&f, label);
    }

    let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
    let report = PmlReport {
        physical_cells: n_phys_cells,
        layer_quads,
        dropped_cells,
        layers: spec.layers,
        layers3: [[0; 2]; 3],
        min_pml_scaled_jacobian,
    };
    Ok(Native::Cells(built.mesh, report))
}

fn single_cell_owners(cm: &CellMesh<2>) -> HashMap<[u32; 4], u32> {
    let fc = cm.face_cells();
    let mut v: HashMap<[u32; 4], u32> = HashMap::new();
    for (k, cs) in fc {
        if cs.len() == 1 {
            v.insert(k, cs[0]);
        }
    }
    v
}

fn outward_edge_normal_2d(cm: &CellMesh<2>, cell_idx: usize, edge: [u32; 2]) -> [f64; 2] {
    let pa = cm.points[edge[0] as usize];
    let pb = cm.points[edge[1] as usize];
    let vs = cm.cell(cell_idx);
    let mut cen = [0.0; 2];
    for &v in vs {
        cen[0] += cm.points[v as usize][0] / vs.len() as f64;
        cen[1] += cm.points[v as usize][1] / vs.len() as f64;
    }
    let dx = pb[0] - pa[0];
    let dz = pb[1] - pa[1];
    let len = dx.hypot(dz).max(1e-300);
    let mut n = [dz / len, -dx / len];
    let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
    if n[0] * (mid[0] - cen[0]) + n[1] * (mid[1] - cen[1]) < 0.0 {
        n = [-n[0], -n[1]];
    }
    n
}

fn collect_outer_geom_edges(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    domain_of_region: &[i32],
) -> HashMap<i32, Vec<u32>> {
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
}

fn sorted_pair(a: u32, b: u32) -> (u32, u32) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// The label of the geom edge, among the outer edges of `domain`, closest to `p` (preferring straight edges whose
/// outward normal aligns with `want_normal` when resolving `Outward` corner facets).
fn nearest_outer_label(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    outer_edges: &HashMap<i32, Vec<u32>>,
    domain: i32,
    p: [f64; 2],
    want_normal: Option<[f64; 2]>,
) -> Result<i32, String> {
    let cands = outer_edges.get(&domain).ok_or_else(|| {
        format!(
            "region '{}' has no edge on the outside of the model, but the mesh has a boundary facet on it",
            g.labels.name_of_domain(domain).unwrap_or("?")
        )
    })?;
    let cx = ev.complex();
    let mut best = (f64::INFINITY, None);
    for pass in 0..2 {
        for &e in cands {
            let is_curved = ev.min_radius_of_curvature(e).is_finite();
            if pass == 0 && !is_curved {
                if let Some(wn) = want_normal {
                    let aligned = cx
                        .facet_normals(e)
                        .into_iter()
                        .any(|n| n[0] * wn[0] + n[2] * wn[1] > 0.35);
                    if !aligned {
                        continue;
                    }
                }
            }
            let (q, _) = ev.edge_closest(e, [p[0], 0.0, p[1]]);
            let d = (q[0] - p[0]).hypot(q[2] - p[1]);
            if d < best.0 {
                best = (d, Some(e));
            }
        }
        if best.1.is_some() {
            break;
        }
    }
    let e = best.1.ok_or("no candidate edge")?;
    g.labels
        .label_of_edge(e)
        .ok_or_else(|| format!("edge {e} has no boundary label"))
}

type PhysPart2d = (
    SimplexMesh,
    HashMap<(u32, u32), i32>,
    Vec<([u32; 2], [f64; 2])>,
);

/// The mesh without the cells of `drop`, vertices renumbered in order (`unchanged` when `drop` removes 0 cells),
/// surviving facet labels by sorted vertex pair, and surviving curved 2D `edge_nodes`.
fn physical_part(m: &SimplexMesh, drop: &HashSet<i32>) -> PhysPart2d {
    let keep: Vec<usize> = (0..m.n_cells()).filter(|&c| !drop.contains(&m.cell_labels[c])).collect();
    if keep.len() == m.n_cells() {
        let mut by_pair = HashMap::with_capacity(m.n_facets());
        for f in 0..m.n_facets() {
            let v = m.facet(f);
            by_pair.insert(sorted_pair(v[0], v[1]), m.facet_labels[f]);
        }
        let en: Vec<([u32; 2], [f64; 2])> = m
            .edge_node_list()
            .into_iter()
            .map(|(pr, p)| (pr, [p[0], p[1]]))
            .collect();
        return (m.clone(), by_pair, en);
    }
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
    let mut remapped_en = Vec::new();
    for ([u, v], p) in m.edge_node_list() {
        let (nu, nv) = (new_id[u as usize], new_id[v as usize]);
        if nu != u32::MAX && nv != u32::MAX && edges.contains(&sorted_pair(nu, nv)) {
            remapped_en.push(([nu.min(nv), nu.max(nv)], [p[0], p[1]]));
        }
    }
    let mut phys = SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels);
    if !remapped_en.is_empty() {
        let en3: Vec<([u32; 2], [f64; 3])> = remapped_en
            .iter()
            .map(|&(pr, p)| (pr, [p[0], p[1], 0.0]))
            .collect();
        phys.set_edge_nodes(&en3);
    }
    (phys, by_pair, remapped_en)
}
