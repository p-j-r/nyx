// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified 3D structured absorbing-layer dispatcher (`PmlMode::{Auto, Structured, Inline}`, ADR 52)
//! across Cartesian (`pml_box`), Cylindrical (`pml_cyl`), and Spherical (`pml_sphere`) topologies
//! in both `Outward` and legacy `Inward` modes.
//!
//! Emits uncut `Tet` + `Prism` + `Hex` `CellMesh<3>` meshes (zero pyramids) preserving physical vertex
//! and cell IDs (`0..n_physical`), `z`-stratified wafer material continuation, radial/lateral `Op::Sheet`
//! `Quad` extrusion, interior dielectric interfaces, and lockstep curved concentric shells.

use crate::adapt::work::{facet_key, simplex_to_adapt, FacetKey};
use crate::mesh::cells::{canonical_quad_face, CellMesh};
use crate::mesh::{Mesh, SimplexMesh};
use crate::meshing::adapt3d::FaceTable;
use crate::meshing::geomdrive::MeshedGeom;
use crate::meshing::pml_box::{build_pml_with_sheets, propagate_cartesian_edge_nodes, PmlSpec};
use crate::meshing::pml_cyl::{
    build_pml_cylinder_3d_with_edges_and_sheets, CylPmlSpec3d, CylRegionRole,
};
use crate::meshing::pml_geom::{compute_min_pml_scaled_jacobian, PmlReport};
use crate::meshing::pml_sphere::{build_pml_sphere_3d_with_edges, SpherePmlSpec3d};
use geom::model::{PmlCoordKind, PmlGrowMode, PmlRole, RegionKind};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// What [`structured_pml_3d_native`] returns.
#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
pub enum Native3 {
    /// The model has no absorbing layer: the meshed model as it was.
    Unchanged(MeshedGeom, PmlReport),
    /// The physical tetrahedra followed by the prisms and hexahedra of the layers, uncut.
    Cells(CellMesh<3>, PmlReport),
}

/// Replaces or extrudes the absorbing layers of a meshed 3D model, storing the uncut `CellMesh<3>` in
/// `MeshedGeom::cells` and splitting `MeshedGeom::mesh` into tetrahedra for legacy simplex callers.
pub fn structured_pml_3d(g: MeshedGeom) -> Result<(MeshedGeom, PmlReport), String> {
    let curved = g.mesh.has_curved();
    structured_pml_3d_with_curved(g, curved)
}

/// Like [`structured_pml_3d`], with explicit `curved` control for concentric curvilinear shells.
pub fn structured_pml_3d_with_curved(
    g: MeshedGeom,
    curved: bool,
) -> Result<(MeshedGeom, PmlReport), String> {
    let curved = curved || g.mesh.has_curved();
    let (ev, labels, warnings) = (g.eval.clone(), g.labels.clone(), g.warnings.clone());
    match structured_pml_3d_native_with_curved(g, curved)? {
        Native3::Unchanged(g, report) => Ok((g, report)),
        Native3::Cells(cells, report) => {
            let mesh = split_3d_cells_preserving_curved_quads(&cells)?;
            mesh.check()?;
            Ok((
                MeshedGeom {
                    mesh,
                    cells: Some(Mesh::D3(cells)),
                    labels,
                    eval: ev,
                    warnings,
                },
                report,
            ))
        }
    }
}

/// Builds uncut 3D mixed `Tet` + `Prism` + `Hex` structured absorbing layers, inferring `curved` from `g.mesh`.
pub fn structured_pml_3d_native(g: MeshedGeom) -> Result<Native3, String> {
    let curved = g.mesh.has_curved();
    structured_pml_3d_native_with_curved(g, curved)
}

/// Builds uncut 3D mixed `Tet` + `Prism` + `Hex` structured absorbing layers across Cartesian, Cylindrical,
/// and Spherical topologies (both `Outward` and `Inward`).
pub fn structured_pml_3d_native_with_curved(g: MeshedGeom, curved: bool) -> Result<Native3, String> {
    let curved = curved || g.mesh.has_curved();
    let ev = g.eval.clone();
    ev.validate_pml_ray_invariance().map_err(|e| e.to_string())?;
    let pieces = ev.pml_regions();
    if pieces.is_empty() {
        let n = g.mesh.n_cells();
        return Ok(Native3::Unchanged(g, PmlReport::empty(n)));
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
    let pml_domains: HashSet<i32> = pieces
        .iter()
        .map(|p| domain_of_region[p.region as usize])
        .collect();
    let mut parent_of_pml: HashMap<i32, i32> = HashMap::new();
    for p in pieces {
        if let Some(par) = p.parent {
            parent_of_pml.insert(domain_of_region[p.region as usize], domain_of_region[par as usize]);
        }
    }

    let (phys, facet_label, phys_edge_nodes, phys_qp_edges) = physical_part(&g.mesh, &pml_domains);
    if phys.n_cells() == 0 {
        return Err("nothing is left of the model outside its absorbing layers".into());
    }
    let dropped_cells = g.mesh.n_cells() - phys.n_cells();
    let interior_sheets = collect_interior_sheet_tris(&ev, &g, &phys);
    let work = simplex_to_adapt::<3>(&phys)?;

    let table = FaceTable::new(&ev, &g.labels)?;
    let outer_faces = collect_outer_geom_faces(&ev, &g, &domain_of_region);
    let stride = g.labels.domains.iter().map(|d| d.id).max().unwrap_or(0) + 1;

    let is_spherical = pieces.iter().any(|p| p.coord == PmlCoordKind::Spherical);
    let is_cylindrical = !is_spherical
        && pieces
            .iter()
            .any(|p| p.coord == PmlCoordKind::Cylindrical && p.r_in > 0.0);

    // ------------------------------------------------------------------ 1. 3D Spherical Shells
    if is_spherical {
        let p0 = &pieces[0];
        let n_layers = (p0.layers.round() as usize).max(1);
        let mut piece_of_mat: HashMap<i32, i32> = HashMap::new();
        for p in pieces {
            if let Some(par) = p.parent {
                piece_of_mat.insert(domain_of_region[par as usize], domain_of_region[p.region as usize]);
            }
        }
        let mut sp = SpherePmlSpec3d::new(p0.center, p0.r_in, p0.thickness, n_layers);
        sp.grading = p0.layer_grading.max(1e-6);
        sp.curved = curved;
        sp.label_base = stride;
        sp.label_stride = stride;
        sp.tol = 5e-3;

        let mut built = build_pml_sphere_3d_with_edges(&work, &sp, &phys_edge_nodes)
            .map_err(|e| format!("structured spherical absorbing layers: {e}"))?;
        if !phys_qp_edges.is_empty() {
            let mut qp: Vec<[u32; 2]> =
                phys_qp_edges.into_iter().map(|[a, b]| [a.min(b), a.max(b)]).collect();
            qp.sort_unstable();
            qp.dedup();
            built.mesh.qp_edges = qp;
        }
        let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);
        let mut layer_cells = 0usize;
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
            layer_cells += 1;
        }

        populate_3d_facets(
            &ev,
            &g,
            &table,
            &outer_faces,
            &parent_of_pml,
            p0.mode,
            &mut built.mesh,
            n_phys_points,
            &facet_label,
            &interior_sheets,
            &[],
            p0.mode == PmlGrowMode::Outward,
        )?;

        let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
        let report = PmlReport {
            physical_cells: n_phys_cells,
            layer_quads: layer_cells,
            dropped_cells,
            layers: [[n_layers, n_layers]; 2],
            layers3: [[n_layers, n_layers]; 3],
            min_pml_scaled_jacobian,
        };
        return Ok(Native3::Cells(built.mesh, report));
    }

    // ------------------------------------------------------------------ 2. 3D Cylindrical Shells / Endcaps / Rim Annuli
    if is_cylindrical {
        let p0 = &pieces[0];
        let ax_len = (p0.axis[0] * p0.axis[0] + p0.axis[1] * p0.axis[1] + p0.axis[2] * p0.axis[2])
            .sqrt()
            .max(1e-30);
        let a_hat = [p0.axis[0] / ax_len, p0.axis[1] / ax_len, p0.axis[2] / ax_len];
        let ax_slot = (0..3)
            .max_by(|&a, &b| a_hat[a].abs().total_cmp(&a_hat[b].abs()))
            .unwrap_or(2);

        let mut z_min = pieces.iter().map(|p| p.z_range[0]).fold(f64::INFINITY, f64::min);
        let mut z_max = pieces.iter().map(|p| p.z_range[1]).fold(f64::NEG_INFINITY, f64::max);
        if !(z_max > z_min && z_min.is_finite() && z_max.is_finite()) {
            z_min = work
                .points
                .iter()
                .map(|q| dot3(*q, a_hat))
                .fold(f64::INFINITY, f64::min);
            z_max = work
                .points
                .iter()
                .map(|q| dot3(*q, a_hat))
                .fold(f64::NEG_INFINITY, f64::max);
        }

        let mut r_layers = 0usize;
        let mut r_thickness = 0.0f64;
        let mut z_min_layers = 0usize;
        let mut z_min_thickness = 0.0f64;
        let mut z_max_layers = 0usize;
        let mut z_max_thickness = 0.0f64;
        let mut piece_of_role: HashMap<(CylRegionRole, i32), i32> = HashMap::new();

        for p in pieces {
            let n = (p.layers.round() as usize).max(1);
            let s_dot = p.sides[0] as f64 * a_hat[0]
                + p.sides[1] as f64 * a_hat[1]
                + p.sides[2] as f64 * a_hat[2];
            let role_opt = match p.role {
                PmlRole::RadialShell => {
                    r_layers = n;
                    r_thickness = p.thickness;
                    Some(CylRegionRole::RadialShell)
                }
                PmlRole::EndCap if s_dot < -0.5 => {
                    z_min_layers = n;
                    z_min_thickness = p.thickness;
                    Some(CylRegionRole::EndCapMin)
                }
                PmlRole::EndCap if s_dot > 0.5 => {
                    z_max_layers = n;
                    z_max_thickness = p.thickness;
                    Some(CylRegionRole::EndCapMax)
                }
                PmlRole::CornerAnnulus if s_dot < -0.5 => {
                    if r_layers == 0 {
                        r_layers = n;
                        r_thickness = p.thickness;
                    }
                    if z_min_layers == 0 {
                        z_min_layers = n;
                        z_min_thickness = p.thickness;
                    }
                    Some(CylRegionRole::CornerAnnulusMin)
                }
                PmlRole::CornerAnnulus if s_dot > 0.5 => {
                    if r_layers == 0 {
                        r_layers = n;
                        r_thickness = p.thickness;
                    }
                    if z_max_layers == 0 {
                        z_max_layers = n;
                        z_max_thickness = p.thickness;
                    }
                    Some(CylRegionRole::CornerAnnulusMax)
                }
                _ => None,
            };
            if let (Some(role), Some(par)) = (role_opt, p.parent) {
                piece_of_role.insert((role, domain_of_region[par as usize]), domain_of_region[p.region as usize]);
            }
        }

        let cs = CylPmlSpec3d {
            center: p0.center,
            axis: a_hat,
            r_in: p0.r_in,
            z_range: [z_min, z_max],
            r_thickness,
            r_layers,
            z_min_thickness,
            z_min_layers,
            z_max_thickness,
            z_max_layers,
            grading: p0.layer_grading.max(1e-6),
            curved,
            label_base: stride,
            label_stride: stride,
            outer_label: None,
            tol: 5e-3,
        };

        let mut built = build_pml_cylinder_3d_with_edges_and_sheets(
            &work,
            &cs,
            &phys_edge_nodes,
            &interior_sheets,
        )
        .map_err(|e| format!("structured cylindrical absorbing layers: {e}"))?;
        if !phys_qp_edges.is_empty() {
            let mut qp: Vec<[u32; 2]> =
                phys_qp_edges.into_iter().map(|[a, b]| [a.min(b), a.max(b)]).collect();
            qp.sort_unstable();
            qp.dedup();
            built.mesh.qp_edges = qp;
        }
        let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);
        let mut layer_cells = 0usize;
        for c in n_phys_cells..built.mesh.n_cells() {
            let (role, mat) = built
                .spec
                .decode(built.mesh.label[c])
                .ok_or_else(|| format!("layer cell {c} carries physical label {}", built.mesh.label[c]))?;
            let id = piece_of_role.get(&(role, mat)).copied().ok_or_else(|| {
                format!(
                    "geom has no cylindrical absorbing piece ({role:?}) for region '{}'",
                    g.labels.name_of_domain(mat).unwrap_or("?")
                )
            })?;
            built.mesh.label[c] = id;
            layer_cells += 1;
        }

        populate_3d_facets(
            &ev,
            &g,
            &table,
            &outer_faces,
            &parent_of_pml,
            p0.mode,
            &mut built.mesh,
            n_phys_points,
            &facet_label,
            &interior_sheets,
            &built.extruded_sheets,
            p0.mode == PmlGrowMode::Outward,
        )?;

        let mut layers3 = [[r_layers, r_layers]; 3];
        layers3[ax_slot] = [z_min_layers, z_max_layers];
        let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
        let report = PmlReport {
            physical_cells: n_phys_cells,
            layer_quads: layer_cells,
            dropped_cells,
            layers: [layers3[0], layers3[2]],
            layers3,
            min_pml_scaled_jacobian,
        };
        return Ok(Native3::Cells(built.mesh, report));
    }

    // ------------------------------------------------------------------ 3. 3D Cartesian Box
    let p0 = &pieces[0];
    let mut spec = PmlSpec::<3>::uniform(0.0, 0);
    spec.layers = [[0; 2]; 3];
    spec.thickness = [[0.0; 2]; 3];
    spec.grading = p0.layer_grading.max(1e-6);

    for p in pieces {
        let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
        if nz.len() != 1 {
            continue;
        }
        let (a, s) = (nz[0], if p.sides[nz[0]] < 0 { 0 } else { 1 });
        if p.mode == PmlGrowMode::Inward {
            let bb = ev.regions()[p.region as usize].bbox;
            let extent = bb.hi[a] - bb.lo[a];
            if (extent - p.thickness).abs() > 1e-6 * p.thickness {
                return Err(format!(
                    "an absorbing piece is {extent:.4e} thick where the Pml asks for {:.4e}: the model ends inside the layer, so the structured layers cannot match it",
                    p.thickness
                ));
            }
        }
        let n = (p.layers.round() as usize).max(1);
        if spec.layers[a][s] != 0
            && (spec.layers[a][s] != n || (spec.thickness[a][s] - p.thickness).abs() > 1e-9 * p.thickness)
        {
            return Err(format!(
                "the absorbing pieces on side {s} of axis {a} disagree in thickness or number of layers; the structured layers need one"
            ));
        }
        spec.layers[a][s] = n;
        spec.thickness[a][s] = p.thickness;
    }

    let phys_bb = work.bbox();
    let win = p0.window;
    spec.interface_box = Some(std::array::from_fn(|a| {
        [
            if spec.layers[a][0] > 0 { win.lo[a] } else { phys_bb[a][0] },
            if spec.layers[a][1] > 0 { win.hi[a] } else { phys_bb[a][1] },
        ]
    }));
    spec.label_base = stride;
    spec.label_stride = stride;
    spec.outer_label = None;

    let mut piece_of: HashMap<([i8; 3], i32), i32> = HashMap::new();
    for p in pieces {
        if let Some(parent) = p.parent {
            piece_of.insert(
                (p.sides, domain_of_region[parent as usize]),
                domain_of_region[p.region as usize],
            );
        }
    }

    let (mut built, extruded_sheets) = build_pml_with_sheets(&work, &spec, &interior_sheets)
        .map_err(|e| format!("structured absorbing layers: {e}"))?;
    propagate_cartesian_edge_nodes(&mut built, &work, &spec, &phys_edge_nodes, &phys_qp_edges);
    let (n_phys_cells, n_phys_points) = (built.n_physical_cells, built.n_physical_points);

    let mut layer_cells = 0usize;
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

    populate_3d_facets(
        &ev,
        &g,
        &table,
        &outer_faces,
        &parent_of_pml,
        p0.mode,
        &mut built.mesh,
        n_phys_points,
        &facet_label,
        &interior_sheets,
        &extruded_sheets,
        false,
    )?;

    let min_pml_scaled_jacobian = compute_min_pml_scaled_jacobian(&built.mesh, n_phys_cells);
    let report = PmlReport {
        physical_cells: n_phys_cells,
        layer_quads: layer_cells,
        dropped_cells,
        layers: [spec.layers[0], spec.layers[2]],
        layers3: spec.layers,
        min_pml_scaled_jacobian,
    };
    Ok(Native3::Cells(built.mesh, report))
}

/// Splits an uncut 3D `CellMesh<3>` into a positively-oriented `SimplexMesh`, projecting split diagonals
/// of curved `Prism`/`Hex` quad faces to match their isoparametric face centers so thin curved shells
/// stay non-inverted even on legacy simplex views.
#[allow(clippy::chunks_exact_to_as_chunks)]
fn split_3d_cells_preserving_curved_quads(cm: &CellMesh<3>) -> Result<SimplexMesh, String> {
    let (mut mesh, _) = SimplexMesh::from_cells_split(cm)?;
    for t in mesh.cells.chunks_exact_mut(4) {
        let p = |v: u32| {
            [
                mesh.points[3 * v as usize],
                mesh.points[3 * v as usize + 1],
                mesh.points[3 * v as usize + 2],
            ]
        };
        let (a, b, c, d) = (p(t[0]), p(t[1]), p(t[2]), p(t[3]));
        let (u, v, w) = (sub3(b, a), sub3(c, a), sub3(d, a));
        if dot3(u, cross3(v, w)) < 0.0 {
            t.swap(2, 3);
        }
    }
    let mut mesh = SimplexMesh::new(
        3,
        mesh.points.clone(),
        mesh.cells.clone(),
        mesh.cell_labels.clone(),
        mesh.facets.clone(),
        mesh.facet_labels.clone(),
    );
    if !cm.edge_nodes.is_empty() {
        let mut edge_map: BTreeMap<[u32; 2], [f64; 3]> = cm
            .edge_nodes
            .iter()
            .map(|&([a, b], p)| ([a.min(b), a.max(b)], p))
            .collect();
        let edge_delta = |u: u32, v: u32, emap: &BTreeMap<[u32; 2], [f64; 3]>| -> [f64; 3] {
            let key = [u.min(v), u.max(v)];
            let Some(&q) = emap.get(&key) else { return [0.0; 3] };
            let (pu, pv) = (cm.points[u as usize], cm.points[v as usize]);
            [
                q[0] - 0.5 * (pu[0] + pv[0]),
                q[1] - 0.5 * (pu[1] + pv[1]),
                q[2] - 0.5 * (pu[2] + pv[2]),
            ]
        };
        let mut seen_quads: BTreeSet<[u32; 4]> = BTreeSet::new();
        let mut diag_additions: Vec<([u32; 2], [f64; 3])> = Vec::new();
        for c in 0..cm.n_cells() {
            let vs = cm.cell(c);
            for f in cm.types[c].faces() {
                if f.len() != 4 {
                    continue;
                }
                let q = [
                    vs[f[0] as usize],
                    vs[f[1] as usize],
                    vs[f[2] as usize],
                    vs[f[3] as usize],
                ];
                let cq = canonical_quad_face(q);
                if !seen_quads.insert(cq) {
                    continue;
                }
                let d01 = edge_delta(q[0], q[1], &edge_map);
                let d12 = edge_delta(q[1], q[2], &edge_map);
                let d23 = edge_delta(q[2], q[3], &edge_map);
                let d30 = edge_delta(q[3], q[0], &edge_map);
                let dc = [
                    0.5 * (d01[0] + d12[0] + d23[0] + d30[0]),
                    0.5 * (d01[1] + d12[1] + d23[1] + d30[1]),
                    0.5 * (d01[2] + d12[2] + d23[2] + d30[2]),
                ];
                if dot3(dc, dc) <= 1e-30 {
                    continue;
                }
                for (da, db) in [(q[0], q[2]), (q[1], q[3])] {
                    let dk = (da.min(db), da.max(db));
                    if mesh.edges().binary_search(&dk).is_ok() && !edge_map.contains_key(&[dk.0, dk.1]) {
                        let (pa, pb) = (cm.points[da as usize], cm.points[db as usize]);
                        diag_additions.push((
                            [dk.0, dk.1],
                            [
                                0.5 * (pa[0] + pb[0]) + dc[0],
                                0.5 * (pa[1] + pb[1]) + dc[1],
                                0.5 * (pa[2] + pb[2]) + dc[2],
                            ],
                        ));
                    }
                }
            }
        }
        for (k, p) in diag_additions {
            edge_map.insert(k, p);
        }
        let nodes: Vec<([u32; 2], [f64; 3])> = edge_map.into_iter().collect();
        mesh.set_edge_nodes(&nodes);
        if !cm.qp_edges.is_empty() {
            mesh.set_qp_edges(&cm.qp_edges);
        }
    }
    Ok(mesh)
}

/// Populates exterior boundary facets, interior `Op::Sheet` facets (`Tri3` and `Quad4`), and interior
/// material-interface facets (`Tri3` and `Quad4`, ADR 40) on `cm`.
#[allow(clippy::too_many_arguments)]
fn populate_3d_facets(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    table: &FaceTable,
    outer_faces: &HashMap<i32, Vec<u32>>,
    parent_of_pml: &HashMap<i32, i32>,
    mode: PmlGrowMode,
    cm: &mut CellMesh<3>,
    n_phys_points: usize,
    phys_facet_label: &HashMap<FacetKey, i32>,
    interior_sheets: &[([u32; 3], i32)],
    extruded_sheets: &[([u32; 4], i32)],
    keep_existing_exterior: bool,
) -> Result<(), String> {
    let face_cells = cm.face_cells();
    let mut emitted_keys: HashSet<FacetKey> = HashSet::new();

    if keep_existing_exterior && cm.n_facets() > 0 {
        for f in 0..cm.n_facets() {
            emitted_keys.insert(facet_key(cm.facet(f)));
        }
    } else {
        cm.facet_offsets = vec![0];
        cm.facet_conn.clear();
        cm.facet_label.clear();
        for f in cm.boundary_faces_ordered() {
            let key = facet_key(&f);
            let owner = face_cells[&key][0] as usize;
            let physical = f.iter().all(|&v| (v as usize) < n_phys_points);
            let label = match physical.then(|| phys_facet_label.get(&key)).flatten() {
                Some(&l) => l,
                None => {
                    let mid = facet_centroid(cm, &f);
                    let n_out = outward_face_normal_3d(cm, owner, &f);
                    let owner_dom = cm.label[owner];
                    let lookup_dom = if outer_faces.contains_key(&owner_dom) {
                        owner_dom
                    } else {
                        parent_of_pml.get(&owner_dom).copied().unwrap_or(owner_dom)
                    };
                    nearest_outer_label_3d(ev, g, table, outer_faces, lookup_dom, mid, Some(n_out))?
                }
            };
            cm.push_facet(&f, label);
            emitted_keys.insert(key);
        }
    }

    // 1. Interior conductive `Op::Sheet` facets (`Tri3` in `phys` + `Quad4` extruded into PML).
    for &(tri, lab) in interior_sheets {
        let key = facet_key(&tri);
        if face_cells.get(&key).is_some_and(|cs| cs.len() == 2) && emitted_keys.insert(key) {
            cm.push_facet(&tri, lab);
        }
    }
    for &(quad, lab) in extruded_sheets {
        let key = facet_key(&quad);
        if face_cells.get(&key).is_some_and(|cs| cs.len() == 2) && emitted_keys.insert(key) {
            cm.push_facet(&quad, lab);
        }
    }

    // 2. Interior domain-interface facets (`cs.len() == 2 && la != lb`, ADR 40):
    //    - In `Inward` mode: all `la != lb` interfaces have a `Complex` face (`table.interface_label(mid, la, lb)`).
    //    - In `Outward` mode: virtual PML pieces have no internal `Complex` wall against their own parent
    //      (where `Selector::Outer` belongs strictly to the outer termination wall!), while interfaces between
    //      distinct physical parent materials (`pa != pb`, e.g. `air|glass` and `air.pml_*|glass.pml_*`) DO carry
    //      the physical material-interface label `table.interface_label(mid, pa, pb)`.
    let eff_dom = |d: i32| -> i32 {
        if mode == PmlGrowMode::Outward {
            parent_of_pml.get(&d).copied().unwrap_or(d)
        } else {
            d
        }
    };
    let mut interior_ifaces: Vec<(FacetKey, usize, i32, i32)> = Vec::new();
    for (&key, cs) in &face_cells {
        if cs.len() != 2 || emitted_keys.contains(&key) {
            continue;
        }
        let (c0, c1) = (cs[0] as usize, cs[1] as usize);
        let (l0, l1) = (cm.label[c0], cm.label[c1]);
        if l0 == l1 {
            continue;
        }
        let (e0, e1) = (eff_dom(l0), eff_dom(l1));
        if mode == PmlGrowMode::Outward && e0 == e1 {
            continue;
        }
        let owner = if l0 < l1 { c0 } else { c1 };
        interior_ifaces.push((key, owner, e0.min(e1), e0.max(e1)));
    }
    interior_ifaces.sort_unstable_by_key(|x| x.0);

    for (key, owner, da, db) in interior_ifaces {
        if let Some(ordered) = cell_outward_face(cm, owner, key) {
            let mid = facet_centroid(cm, &ordered);
            let lab = table.interface_label(mid, da, db);
            cm.push_facet(&ordered, lab);
            emitted_keys.insert(key);
        }
    }
    Ok(())
}

fn cell_outward_face(cm: &CellMesh<3>, cell_idx: usize, want_key: FacetKey) -> Option<Vec<u32>> {
    let vs = cm.cell(cell_idx);
    for f in cm.types[cell_idx].faces() {
        let mut ordered: Vec<u32> = f.iter().map(|&k| vs[k as usize]).collect();
        if facet_key(&ordered) == want_key {
            let n_f = if ordered.len() >= 3 {
                let (p0, p1, p2) = (
                    cm.points[ordered[0] as usize],
                    cm.points[ordered[1] as usize],
                    cm.points[ordered[2] as usize],
                );
                cross3(sub3(p1, p0), sub3(p2, p0))
            } else {
                [0.0; 3]
            };
            let cell_cen = facet_centroid(cm, vs);
            let face_cen = facet_centroid(cm, &ordered);
            if dot3(n_f, sub3(face_cen, cell_cen)) < 0.0 {
                if ordered.len() == 3 {
                    ordered.swap(1, 2);
                } else if ordered.len() == 4 {
                    ordered.swap(1, 3);
                }
            }
            return Some(ordered);
        }
    }
    None
}

fn facet_centroid(cm: &CellMesh<3>, vs: &[u32]) -> [f64; 3] {
    let inv = 1.0 / (vs.len().max(1) as f64);
    let mut c = [0.0; 3];
    for &v in vs {
        let p = cm.points[v as usize];
        c[0] += p[0] * inv;
        c[1] += p[1] * inv;
        c[2] += p[2] * inv;
    }
    c
}

fn outward_face_normal_3d(cm: &CellMesh<3>, cell_idx: usize, f: &[u32]) -> [f64; 3] {
    if f.len() < 3 {
        return [0.0, 0.0, 1.0];
    }
    let (p0, p1, p2) = (
        cm.points[f[0] as usize],
        cm.points[f[1] as usize],
        cm.points[f[2] as usize],
    );
    let mut n = cross3(sub3(p1, p0), sub3(p2, p0));
    let len = dot3(n, n).sqrt().max(1e-300);
    n = [n[0] / len, n[1] / len, n[2] / len];
    let cell_cen = facet_centroid(cm, cm.cell(cell_idx));
    let face_cen = facet_centroid(cm, f);
    if dot3(n, sub3(face_cen, cell_cen)) < 0.0 {
        n = [-n[0], -n[1], -n[2]];
    }
    n
}

fn collect_outer_geom_faces(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    domain_of_region: &[i32],
) -> HashMap<i32, Vec<u32>> {
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
}

fn collect_interior_sheet_tris(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    phys: &SimplexMesh,
) -> Vec<([u32; 3], i32)> {
    let cx = ev.complex();
    let sheet_labels: HashSet<i32> = cx
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.is_sheet())
        .filter_map(|(fi, _)| g.labels.label_of_edge(fi as u32))
        .collect();
    let mut phys_face_doms: HashMap<[u32; 3], Vec<i32>> = HashMap::new();
    for c in 0..phys.n_cells() {
        let v = phys.cell(c);
        let dom = phys.cell_labels[c];
        for skip in 0..4 {
            let mut k = [0u32; 3];
            let mut m = 0;
            for (i, &vi) in v.iter().enumerate() {
                if i != skip {
                    k[m] = vi;
                    m += 1;
                }
            }
            k.sort_unstable();
            phys_face_doms.entry(k).or_default().push(dom);
        }
    }
    let mut map: BTreeMap<[u32; 3], ([u32; 3], i32)> = BTreeMap::new();
    for f in 0..phys.n_facets() {
        let fv = phys.facet(f);
        let tri = [fv[0], fv[1], fv[2]];
        let mut k = tri;
        k.sort_unstable();
        let lab = phys.facet_labels[f];
        if let Some(doms) = phys_face_doms.get(&k) {
            if doms.len() == 2 && (sheet_labels.contains(&lab) || doms[0] == doms[1]) {
                map.insert(k, (tri, lab));
            }
        }
    }
    map.into_values().collect()
}

fn nearest_outer_label_3d(
    ev: &geom::eval::Evaluated,
    g: &MeshedGeom,
    table: &FaceTable,
    outer_faces: &HashMap<i32, Vec<u32>>,
    domain: i32,
    p: [f64; 3],
    want_normal: Option<[f64; 3]>,
) -> Result<i32, String> {
    let cands = outer_faces.get(&domain).ok_or_else(|| {
        format!(
            "region '{}' has no face on the outside of the model, but the mesh has a boundary facet on it",
            g.labels.name_of_domain(domain).unwrap_or("?")
        )
    })?;
    let cx = ev.complex();
    if let Some(wn) = want_normal {
        let aligned: Vec<usize> = cands
            .iter()
            .copied()
            .filter(|&fi| {
                if cx.facet_radius(fi).is_some() {
                    true
                } else {
                    cx.facet_normals(fi).into_iter().any(|n| dot3(n, wn) > 0.35)
                }
            })
            .map(|fi| fi as usize)
            .collect();
        if let Some(lab) = table.nearest_of(p, aligned.into_iter()) {
            return Ok(lab);
        }
    }
    table
        .nearest_of(p, cands.iter().map(|&i| i as usize))
        .ok_or_else(|| "no candidate face".to_string())
}

type PhysPart3d = (
    SimplexMesh,
    HashMap<FacetKey, i32>,
    Vec<([u32; 2], [f64; 3])>,
    Vec<[u32; 2]>,
);

/// The mesh without the cells of `drop`, vertices renumbered in order (`unchanged` when `drop` removes 0 cells),
/// surviving facet labels by face key, surviving `edge_nodes`, and surviving `qp_edges`.
fn physical_part(m: &SimplexMesh, drop: &HashSet<i32>) -> PhysPart3d {
    let keep: Vec<usize> = (0..m.n_cells()).filter(|&c| !drop.contains(&m.cell_labels[c])).collect();
    if keep.len() == m.n_cells() {
        let mut by_key = HashMap::with_capacity(m.n_facets());
        for f in 0..m.n_facets() {
            let v = m.facet(f);
            by_key.insert(facet_key(&[v[0], v[1], v[2]]), m.facet_labels[f]);
        }
        return (m.clone(), by_key, m.edge_node_list(), m.qp_edge_list());
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
    let mut faces: HashSet<[u32; 3]> = HashSet::new();
    let mut edges: HashSet<(u32, u32)> = HashSet::new();
    for &c in &keep {
        let v: Vec<u32> = m.cell(c).iter().map(|&x| new_id[x as usize]).collect();
        for i in 0..4 {
            for j in (i + 1)..4 {
                edges.insert((v[i].min(v[j]), v[i].max(v[j])));
            }
        }
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
        by_key.insert(facet_key(&t), m.facet_labels[f]);
    }
    let mut remapped_en = Vec::new();
    for ([u, v], p) in m.edge_node_list() {
        let (nu, nv) = (new_id[u as usize], new_id[v as usize]);
        if nu != u32::MAX && nv != u32::MAX && edges.contains(&(nu.min(nv), nu.max(nv))) {
            remapped_en.push(([nu.min(nv), nu.max(nv)], p));
        }
    }
    let mut remapped_qp = Vec::new();
    for [u, v] in m.qp_edge_list() {
        let (nu, nv) = (new_id[u as usize], new_id[v as usize]);
        if nu != u32::MAX && nv != u32::MAX && edges.contains(&(nu.min(nv), nu.max(nv))) {
            remapped_qp.push([nu.min(nv), nu.max(nv)]);
        }
    }
    let mut phys = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
    if !remapped_en.is_empty() {
        phys.set_edge_nodes(&remapped_en);
        if !remapped_qp.is_empty() {
            phys.set_qp_edges(&remapped_qp);
        }
    }
    (phys, by_key, remapped_en, remapped_qp)
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
