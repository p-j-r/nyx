// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 0B + Phase 3 integration suite: uncut mixed-cell `Quad` / `Prism` / `Hex` flow,
//! `PmlModeRequest::Auto` defaults and legacy backward compatibility, ray-invariance validation,
//! and lockstep curvilinear structured PMLs across all 7 topologies (`2D Cartesian`,
//! `2D-axi Cartesian/Cylindrical`, `2D Cylindrical Annulus`, `2D-axi Spherical Semicircle`,
//! `3D Cartesian` + `Op::Sheet`, `3D Cylindrical` multi-layer wafer + `Op::Sheet`, and
//! `3D Spherical`).

use femcore::fem_mixed::{quad_rule_mixed, IsoparametricCell};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::mixed::{FacetType, MixedDiscretisation};
use femcore::mesh::Mesh;
use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
use femcore::meshing::geomdrive3::{mesh_geom_3d_with, Sizes3};
use femcore::meshing::pml_geom3::{structured_pml_3d_native_with_curved, Native3};
use femcore::meshing::request::{mesh_request, mesh_request_native, MeshRequest, WavesRequest};
use femcore::meshing::{mesh_model, MeshPlan, Sizing};
use femcore::pml::pml_jacobian;
use femcore::props::{Pml, PmlCoord};
use femcore::section::wafer_point;
use geom::expr::Expr;
use geom::model::{Keep, Model, PmlCoordKind, RegionKind, Side};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// Verifies `s_dist` lies strictly inside `[0.0, 1.0]` across all degree-2 quadrature points of every
/// PML cell in `cm`, and spans `(min < 0.25, max > 0.75)` across the active absorbing layer.
fn verify_pml_sdist_bounds<const D: usize>(
    cm: &CellMesh<D>,
    md: &MixedDiscretisation,
    pml_map: &BTreeMap<i32, Pml>,
    axisym: bool,
) {
    let mut s_min = f64::INFINITY;
    let mut s_max = f64::NEG_INFINITY;
    let mut n_checked = 0usize;
    let mut scratch_n = Vec::with_capacity(27);
    let mut scratch_dn = Vec::with_capacity(27);

    for ci in 0..cm.n_cells() {
        let Some(spec) = pml_map.get(&cm.label[ci]) else {
            continue;
        };
        assert!(spec.is_active(), "PML domain {} must be active", cm.label[ci]);
        let iso = IsoparametricCell::from_mesh(md, ci);
        let qrule = quad_rule_mixed(iso.cell_type, 2);
        for xi in &qrule.pts {
            let je = iso.eval_into(xi, &mut scratch_n, &mut scratch_dn);
            assert!(je.det > 0.0, "cell {ci} ({:?}) has non-positive det J = {}", iso.cell_type, je.det);
            let wp = if D == 2 {
                wafer_point(&[je.x[0], je.x[1]], 2)
            } else {
                je.x
            };
            let jac = pml_jacobian(spec, None, wp, D, axisym);
            assert!(
                (0.0..=1.0).contains(&jac.s_dist) && jac.s_dist.is_finite(),
                "cell {ci} s_dist = {} out of [0, 1] at wp = {wp:?}",
                jac.s_dist
            );
            assert!(
                jac.det.re.is_finite() && jac.det.im.is_finite() && jac.det.norm_sqr() > 0.0,
                "cell {ci} singular complex PML det = {:?}",
                jac.det
            );
            s_min = s_min.min(jac.s_dist);
            s_max = s_max.max(jac.s_dist);
            n_checked += 1;
        }
    }
    assert!(n_checked > 0, "expected at least one PML quadrature point");
    assert!(s_min < 0.25, "expected inner quadrature points near s_dist -> 0, got min = {s_min}");
    assert!(s_max > 0.75, "expected outer quadrature points near s_dist -> 1, got max = {s_max}");
}

// ------------------------------------------------------------------ 1. 2D Cartesian (Outward & Inward)

#[test]
fn topology_1_2d_cartesian_outward_and_inward_uncut_quads_and_grading() {
    let mut m = Model::new_2d();
    let air = m.rect("air", (-1.5, 0.0), (1.5, 1.0));
    let sub = m.rect("sub", (-1.5, -0.8), (1.5, 0.0));
    m.overlay(&[("air", air), ("sub", sub)]);
    m.pml_outward("pml", PmlCoordKind::Cartesian, &[], 0.4, 6.0, 2.5);

    let req = MeshRequest {
        h: 0.35,
        ..MeshRequest::default()
    };
    let (ex_mesh, labels, _eval, _, rep, warnings) = mesh_request_native(&m.list, &req).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
    assert!(rep.pml_warnings.is_empty(), "{:?}", rep.pml_warnings);
    assert_eq!(rep.structured_layers, Some(vec![[6, 6], [6, 6]]));
    assert!(rep.min_pml_scaled_jacobian.unwrap() > 0.99, "{:?}", rep.min_pml_scaled_jacobian);
    assert!(rep.element_counts["tri"] > 0 && rep.element_counts["quad"] > 0);
    assert!(!rep.element_counts.contains_key("tet"));

    let Mesh::D2(cm) = &ex_mesh else { panic!("expected 2D CellMesh") };
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();

    // Total area must equal (3.0 + 2*0.4) * (1.8 + 2*0.4) = 3.8 * 2.6 = 9.88
    let total_area: f64 = (0..cm.n_cells()).map(|c| cm.measure(c)).sum();
    assert!((total_area - 3.8 * 2.6).abs() < 1e-11, "area = {total_area}");

    // Verify `layer_grading = 2.5` graded the +x Quad widths while `order = 2.0` remained decoupled!
    let (g, _) = mesh_request(&m.list, &req).unwrap();
    let pml_mats = g.pml_materials(3.5);
    let air_xmax = labels.domain("air.pml_xmax").unwrap();
    assert_eq!(pml_mats[&air_xmax].order, 2.0, "stretch order must stay 2.0 when layer_grading = 2.5");
    assert_eq!(pml_mats[&air_xmax].coord, PmlCoord::Cartesian);

    let mut xs: Vec<f64> = cm.points.iter().map(|p| p[0]).filter(|&x| x >= 1.5 - 1e-12).collect();
    xs.sort_by(|a, b| a.total_cmp(b));
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
    assert_eq!(xs.len(), 7, "expected 6 graded layers (7 x planes): {xs:?}");
    let w_first = xs[1] - xs[0];
    let w_last = xs[6] - xs[5];
    assert!(((w_last / w_first) - 2.5).abs() < 1e-9, "ratio = {}", w_last / w_first);

    verify_pml_sdist_bounds(cm, &md, &pml_mats, false);
}

// ------------------------------------------------------------------ 2. 2D Axisymmetric Cartesian/Cylindrical

#[test]
fn topology_2_2d_axisym_rect_touching_r0_omits_xmin_and_preserves_axis() {
    let mut m = Model::new_2d();
    let air = m.rect("air", (0.0, -1.0), (1.2, 1.0));
    let disk = m.rect("disk", (0.0, -0.2), (0.6, 0.2));
    m.overlay(&[("air", air), ("disk", disk)]);
    m.absorbing_layer("pml", PmlCoordKind::Cylindrical, 0.3, 6.0);

    let req = MeshRequest { h: 0.25, ..MeshRequest::default() };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert_eq!(rep.structured_layers, Some(vec![[0, 6], [6, 6]]), "XMin (r = 0) must stay inactive");
    assert!(g.domain("air.pml_xmin").is_none(), "no PML region may be created on r = 0");
    assert!(g.domain("air.pml_radial").is_some() || g.domain("air.pml_xmax").is_some());

    let Mesh::D2(cm) = g.cells.as_ref().unwrap() else { panic!("2D") };
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();
    assert!(cm.points.iter().all(|p| p[0] >= -1e-15), "all axisymmetric points must keep r >= 0");

    let pml_mats = g.pml_materials(3.0);
    verify_pml_sdist_bounds(cm, &md, &pml_mats, true);
}

// ------------------------------------------------------------------ 3. 2D Planar Cylindrical Annulus

#[test]
fn topology_3_2d_cylindrical_annulus_thin_curved_shells_and_inward_parity() {
    // Ultra-thin shell test: R_in = 1.0, thickness = 0.06, N_r = 6 -> h_R / R_in = 0.01!
    let r_in = 1.0;
    let thickness = 0.06;
    let mut m = Model::new_2d();
    let cladding = m.circle("cladding", (0.0, 0.0), r_in);
    let core = m.circle("core", (0.0, 0.0), 0.35);
    m.overlay(&[("cladding", cladding), ("core", core)]);
    m.absorbing_layer("pml", PmlCoordKind::Auto, thickness, 6.0);

    // Coarse angular chords (`max_dev = 0.04`) so chord sagitta (~0.04) is 4x larger than layer thickness (0.01)!
    let req = MeshRequest {
        h: 0.45,
        max_dev: Some(0.04),
        curved: true,
        ..MeshRequest::default()
    };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert!(rep.pml_warnings.is_empty());
    let min_sj = rep.min_pml_scaled_jacobian.expect("min_pml_scaled_jacobian must be reported");
    assert!(
        min_sj > 0.3,
        "2D curved cylindrical annulus at h_R/R_in = 0.01 must have min_pml_scaled_jacobian > 0.3, got {min_sj}"
    );

    let Mesh::D2(cm) = g.cells.as_ref().unwrap() else { panic!("2D") };
    assert!(!cm.edge_nodes.is_empty(), "curved annulus must carry lockstep 2D edge_nodes");
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();

    let pml_mats = g.pml_materials(3.0);
    let clad_pml = g.domain("cladding.pml_radial").unwrap();
    assert!(matches!(pml_mats[&clad_pml].coord, PmlCoord::Cylindrical { radial_on: true, .. }));
    verify_pml_sdist_bounds(cm, &md, &pml_mats, false);

    // Also verify legacy Inward 2D circular PML works seamlessly!
    let mut m_in = Model::new_2d();
    let outer = m_in.circle("cladding", (0.0, 0.0), 1.3);
    let win = m_in.circle("win", (0.0, 0.0), 1.0);
    m_in.overlay(&[("cladding", outer)]);
    m_in.pml("pml", win, &[], 0.3, 5.0, 2.0);
    let (ex_in, _, _, _, rep_in, _) = mesh_request_native(
        &m_in.list,
        &MeshRequest { h: 0.35, curved: true, ..MeshRequest::default() },
    )
    .unwrap();
    assert!(rep_in.min_pml_scaled_jacobian.unwrap() > 0.3);
    let Mesh::D2(cm_in) = &ex_in else { panic!("2D") };
    MixedDiscretisation::from_cell_mesh(cm_in).unwrap().validate_conforming().unwrap();
}

// ------------------------------------------------------------------ 4. 2D Axisymmetric Spherical Semicircle

#[test]
fn topology_4_2d_axisym_spherical_semicircle_clamps_r0_and_tags_symmetry_axis() {
    // Semicircle in x >= 0 around (0, 0) with R_in = 1.0, thickness = 0.06, N_R = 6 (h_R / R_in = 0.01)
    let r_in = 1.0;
    let mut m = Model::new_2d();
    let disk = m.circle("air", (0.0, 0.0), r_in);
    m.overlay(&[("air", disk)]);
    m.cut("axisym_half", (1.0, 0.0, 0.0), 0.0, Keep::Above);
    m.absorbing_layer("pml", PmlCoordKind::Auto, 0.06, 6.0);

    let req = MeshRequest {
        h: 0.4,
        max_dev: Some(0.035),
        curved: true,
        ..MeshRequest::default()
    };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    let min_sj = rep.min_pml_scaled_jacobian.expect("min_pml_scaled_jacobian");
    assert!(
        min_sj > 0.3,
        "2D-axi spherical semicircle at h_R/R_in = 0.01 must have min_pml_scaled_jacobian > 0.3, got {min_sj}"
    );

    let Mesh::D2(cm) = g.cells.as_ref().unwrap() else { panic!("2D") };
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();

    // Every vertex near r = 0 (including newly extruded north/south pole shells z > 1.0 and z < -1.0)
    // must be clamped to exact IEEE `+0.0`, and every boundary facet on r = 0 must carry the symmetry
    // boundary label (never the outer spherical arc label!).
    let sym_labels: BTreeSet<i32> = g.select("plane(x = 0)").unwrap().into_iter().collect();
    let outer_arc_labels: BTreeSet<i32> = g
        .select("outer")
        .unwrap()
        .into_iter()
        .filter(|l| !sym_labels.contains(l))
        .collect();
    assert!(!sym_labels.is_empty() && !outer_arc_labels.is_empty());

    let mut extruded_axis_facets = 0usize;
    for f in 0..cm.n_facets() {
        let fv = cm.facet(f);
        let pa = cm.points[fv[0] as usize];
        let pb = cm.points[fv[1] as usize];
        if pa[0].abs() < 1e-9 && pb[0].abs() < 1e-9 {
            assert_eq!(pa[0].to_bits(), 0.0_f64.to_bits(), "pole vertex r must be exact IEEE 0.0, got {}", pa[0]);
            assert_eq!(pb[0].to_bits(), 0.0_f64.to_bits(), "pole vertex r must be exact IEEE 0.0, got {}", pb[0]);
            let lab = cm.facet_label[f];
            assert!(
                sym_labels.contains(&lab) && !outer_arc_labels.contains(&lab),
                "axial r=0 facet {f} ({pa:?} -> {pb:?}) got non-symmetry label {lab} (sym={sym_labels:?}, arc={outer_arc_labels:?})"
            );
            if pa[1].abs() > r_in + 1e-6 || pb[1].abs() > r_in + 1e-6 {
                extruded_axis_facets += 1;
            }
        }
    }
    assert_eq!(extruded_axis_facets, 12, "expected 6 north + 6 south extruded r=0 symmetry facets");

    let pml_mats = g.pml_materials(3.0);
    let pml_dom = g.domain("air.pml_spherical").unwrap();
    assert!(matches!(pml_mats[&pml_dom].coord, PmlCoord::Spherical { .. }));
    verify_pml_sdist_bounds(cm, &md, &pml_mats, true);
}

// ------------------------------------------------------------------ 5. 3D Cartesian + Op::Sheet Continuation

#[test]
fn topology_5_3d_cartesian_outward_preserves_ids_and_extrudes_lateral_sheets() {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-1.0, -1.0, -0.5), (1.0, 1.0, 0.0));
    let vac = m.box_("vac", (-1.0, -1.0, 0.0), (1.0, 1.0, 0.6));
    // Microstrip trace spanning y in [-1.0, +1.0] so it touches both YMin and YMax PML walls!
    let trace_2d = m.polygons(
        "trace_poly",
        vec![vec![pt(-0.2, -1.0), pt(0.2, -1.0), pt(0.2, 1.0), pt(-0.2, 1.0)]],
        vec![],
    );
    m.sheet("trace", trace_2d, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);
    m.pml_outward("pml", PmlCoordKind::Cartesian, &[], 0.3, 5.0, 2.0);

    let req = MeshRequest { h: 0.45, ..MeshRequest::default() };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert!(rep.pml_warnings.is_empty());
    assert!(rep.element_counts["tet"] > 0 && rep.element_counts["prism"] > 0 && rep.element_counts["hex"] > 0);
    assert!(rep.min_pml_scaled_jacobian.unwrap() > 0.3);

    let Mesh::D3(cm) = g.cells.as_ref().unwrap() else { panic!("3D") };
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();

    // Verify conductive `trace` sheet has BOTH physical Tri3 facets AND laterally extruded Quad4 facets in +/-y PMLs!
    let trace_labs: BTreeSet<i32> = g.boundary("trace").into_iter().collect();
    assert!(!trace_labs.is_empty());
    let mut n_sheet_tri = 0usize;
    let mut n_sheet_quad = 0usize;
    let mut quad_area = 0.0;
    for f in 0..md.n_facets() {
        if trace_labs.contains(&md.facet_labels[f]) {
            assert!(!md.exterior_facets[f], "sheet facet must be interior (2-cell)");
            match md.facet_types[f] {
                FacetType::Tri3 => n_sheet_tri += 1,
                FacetType::Quad4 => {
                    n_sheet_quad += 1;
                    let fv = md.facet(f);
                    let p0 = cm.points[fv[0] as usize];
                    let p1 = cm.points[fv[1] as usize];
                    let p2 = cm.points[fv[2] as usize];
                    let p3 = cm.points[fv[3] as usize];
                    let c1 = [
                        (p1[1] - p0[1]) * (p2[2] - p0[2]) - (p1[2] - p0[2]) * (p2[1] - p0[1]),
                        (p1[2] - p0[2]) * (p2[0] - p0[0]) - (p1[0] - p0[0]) * (p2[2] - p0[2]),
                        (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]),
                    ];
                    let c2 = [
                        (p2[1] - p0[1]) * (p3[2] - p0[2]) - (p2[2] - p0[2]) * (p3[1] - p0[1]),
                        (p2[2] - p0[2]) * (p3[0] - p0[0]) - (p2[0] - p0[0]) * (p3[2] - p0[2]),
                        (p2[0] - p0[0]) * (p3[1] - p0[1]) - (p2[1] - p0[1]) * (p3[0] - p0[0]),
                    ];
                    quad_area += 0.5 * ((c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]).sqrt()
                        + (c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]).sqrt());
                }
                FacetType::Edge2 => unreachable!(),
            }
        }
    }
    assert!(n_sheet_tri > 0, "expected physical Tri3 sheet facets");
    assert!(n_sheet_quad >= 10, "expected >= 10 extruded Quad4 sheet facets across +/-y PMLs, got {n_sheet_quad}");
    // Expected extruded trace area across -y (0.3) and +y (0.3) of width 0.4 = 2 * 0.4 * 0.3 = 0.24
    assert!((quad_area - 0.24).abs() < 1e-10, "extruded trace Quad4 area = {quad_area} (expected 0.24)");

    // Verify `Selector::Outer` facets are strictly exterior (`md.exterior_facets[f] == true`)!
    let outer_labs: BTreeSet<i32> = g.select("outer").unwrap().into_iter().collect();
    for f in 0..md.n_facets() {
        if outer_labs.contains(&md.facet_labels[f]) {
            assert!(
                md.exterior_facets[f],
                "Selector::Outer label {} appeared on an interior facet!",
                md.facet_labels[f]
            );
        }
    }

    let pml_mats = g.pml_materials(3.0);
    verify_pml_sdist_bounds(cm, &md, &pml_mats, false);
}

// ------------------------------------------------------------------ 6. 3D Cylindrical Wafer + Op::Sheet

#[test]
fn topology_6_3d_cylindrical_multilayer_wafer_thin_curved_shells_and_radial_sheets() {
    // 2-layer cylindrical wafer (r_in = 1.0, sub z in [-0.4, 0], air z in [0, 0.5]) + circular ground sheet at z = 0
    // Thin radial/axial PML: thickness = 0.06, layers = 6 -> h_R / R_in = 0.01!
    let r_in = 1.0;
    let mut m = Model::new_3d();
    let sub = m.cylinder("sub", (0.0, 0.0, -0.4), r_in, 0.4);
    let air = m.cylinder("air", (0.0, 0.0, 0.0), r_in, 0.5);
    let gnd_2d = m.circle("gnd_2d", (0.0, 0.0), r_in);
    m.sheet("gnd", gnd_2d, 0.0);
    m.overlay(&[("sub", sub), ("air", air)]);
    m.absorbing_layer("pml", PmlCoordKind::Auto, 0.06, 6.0);

    let req = MeshRequest {
        h: 0.42,
        max_dev: Some(0.035),
        curved: true,
        ..MeshRequest::default()
    };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert!(rep.pml_warnings.is_empty());
    assert!(rep.element_counts["tet"] > 0 && rep.element_counts["prism"] > 0 && rep.element_counts["hex"] > 0);
    let min_sj = rep.min_pml_scaled_jacobian.expect("min_pml_scaled_jacobian");
    assert!(
        min_sj > 0.3,
        "3D curved cylindrical wafer at h_R/R_in = 0.01 must have min_pml_scaled_jacobian > 0.3, got {min_sj}"
    );
    assert!(rep.curved.as_ref().unwrap().min_scaled_jacobian > 0.1);

    let Mesh::D3(cm) = g.cells.as_ref().unwrap() else { panic!("3D") };
    assert!(!cm.edge_nodes.is_empty() && !cm.quad_face_nodes.is_empty());
    let md = MixedDiscretisation::from_cell_mesh(cm).unwrap();
    md.validate_conforming().unwrap();

    // Verify z-stratified barrel materials AND corner rim Hex materials!
    for name in [
        "sub.pml_radial",
        "air.pml_radial",
        "sub.pml_zmin",
        "air.pml_zmax",
        "sub.pml_radial_zmin",
        "air.pml_radial_zmax",
    ] {
        let dom = g.domain(name).unwrap_or_else(|| panic!("missing domain {name}"));
        assert!(cm.label.contains(&dom), "no cells emitted for {name}");
    }

    // Verify radial extrusion of `gnd` sheet at z = 0 into conductive Quad4 facets inside barrel PML!
    let gnd_labs: BTreeSet<i32> = g.boundary("gnd").into_iter().collect();
    let n_gnd_quads = (0..md.n_facets())
        .filter(|&f| gnd_labs.contains(&md.facet_labels[f]) && md.facet_types[f] == FacetType::Quad4)
        .count();
    assert!(n_gnd_quads >= 36, "expected >= 36 radially extruded gnd Quad4 facets, got {n_gnd_quads}");

    let pml_mats = g.pml_materials(3.0);
    verify_pml_sdist_bounds(cm, &md, &pml_mats, false);
}

// ------------------------------------------------------------------ 7. 3D Spherical (Outward & Inward)

#[test]
fn topology_7_3d_spherical_thin_curved_prism_shells_preserve_ids_and_jacobian() {
    // R_in = 1.0, thickness = 0.06, N_R = 6 -> h_R / R_in = 0.01 with `curved = true`!
    let r_in = 1.0;
    let mut m = Model::new_3d();
    let host = m.box_("host", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let air = m.sphere("air", (0.0, 0.0, 0.0), r_in);
    let core = m.sphere("core", (0.0, 0.0, 0.0), 0.45);
    m.overlay_full(
        &[
            ("host", host, RegionKind::Hole, None),
            ("air", air, RegionKind::Material, None),
            ("core", core, RegionKind::Material, None),
        ],
        &[],
    );
    m.absorbing_layer("pml", PmlCoordKind::Auto, 0.06, 6.0);

    // First verify physical vertex/cell ID preservation via `structured_pml_3d_native_with_curved`:
    let ev = evaluate(&m.list).unwrap();
    let mut skip = vec![false; ev.regions().len()];
    for p in ev.pml_regions() {
        skip[p.region as usize] = true;
    }
    let mut opts = GeomMeshOpts::new(0.45);
    opts.curved = true;
    let (phys, _, _) = mesh_geom_3d_with(ev, &opts, &Sizes3 { region_h: vec![], skip }).unwrap();
    let (np, nc) = (phys.mesh.n_points(), phys.mesh.n_cells());
    let Native3::Cells(cells, pml_rep) = structured_pml_3d_native_with_curved(phys.clone(), true).unwrap() else {
        panic!("expected Native3::Cells");
    };
    assert_eq!((pml_rep.physical_cells, pml_rep.dropped_cells), (nc, 0));
    let min_sj = pml_rep.min_pml_scaled_jacobian.expect("min_pml_scaled_jacobian");
    assert!(
        min_sj > 0.3,
        "3D curved spherical shells at h_R/R_in = 0.01 must have min_pml_scaled_jacobian > 0.3, got {min_sj}"
    );
    // Physical vertices and cells must keep exact IDs 0..np and 0..nc!
    for v in 0..np {
        assert_eq!(cells.points[v], phys.mesh.point(v), "physical vertex {v} shifted");
    }
    for c in 0..nc {
        assert_eq!(cells.cell(c), phys.mesh.cell(c), "physical cell {c} connectivity changed");
        assert_eq!(cells.label[c], phys.mesh.cell_labels[c]);
    }
    // Zero Hex and zero pyramids: pure Tet interior + Prism shells!
    let counts = cells.counts();
    assert_eq!(counts[&CellType::Tet], nc);
    assert_eq!(counts[&CellType::Prism], pml_rep.layer_quads);
    assert_eq!(counts.get(&CellType::Hex).copied().unwrap_or(0), 0);

    let md = MixedDiscretisation::from_cell_mesh(&cells).unwrap();
    md.validate_conforming().unwrap();

    // Also run end-to-end `mesh_request` and `mesh_request_native` and verify `sDist` bounds!
    let req = MeshRequest { h: 0.45, curved: true, ..MeshRequest::default() };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert!(rep.min_pml_scaled_jacobian.unwrap() > 0.3);
    assert!(rep.curved.as_ref().unwrap().min_scaled_jacobian > 0.1);
    let pml_mats = g.pml_materials(3.0);
    verify_pml_sdist_bounds(&cells, &md, &pml_mats, false);

    // Also verify Inward 3D Spherical PML mode!
    let mut m_in = Model::new_3d();
    let host_in = m_in.box_("host", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let outer_ball = m_in.sphere("air", (0.0, 0.0, 0.0), 1.3);
    let win_ball = m_in.sphere("win", (0.0, 0.0, 0.0), 1.0);
    m_in.overlay_full(
        &[
            ("host", host_in, RegionKind::Hole, None),
            ("air", outer_ball, RegionKind::Material, None),
        ],
        &[],
    );
    m_in.pml("pml", win_ball, &[], 0.3, 5.0, 2.0);
    let (ex_in, _, _, _, rep_in, _) = mesh_request_native(&m_in.list, &req).unwrap();
    assert!(rep_in.min_pml_scaled_jacobian.unwrap() > 0.3);
    let Mesh::D3(cm_in) = &ex_in else { panic!("3D") };
    MixedDiscretisation::from_cell_mesh(cm_in).unwrap().validate_conforming().unwrap();
}

// ------------------------------------------------------------------ Ray-invariance, warnings & sizing

#[test]
fn ray_invariance_guards_reject_invalid_curvilinear_topologies() {
    // 1. 3D horizontally stratified wafer with Spherical PML violates d_R eps == 0
    let mut m3 = Model::new_3d();
    let sub = m3.cylinder("sub", (0.0, 0.0, -0.5), 1.0, 0.5);
    let air = m3.cylinder("air", (0.0, 0.0, 0.0), 1.0, 0.5);
    m3.overlay(&[("sub", sub), ("air", air)]);
    m3.absorbing_layer("pml", PmlCoordKind::Spherical, 0.2, 6.0);
    let err3 = mesh_model(&m3.list, &MeshPlan::simple(0.4)).err().expect("must reject spherical PML on layered wafer");
    assert!(err3.contains("Spherical PML") && err3.contains("Cylindrical"), "err3 = {err3}");

    // 2. 2D off-center horizontal substrate plane inside a Circle with Cylindrical PML violates d_r eps == 0
    let mut m2 = Model::new_2d();
    let disk = m2.circle("air", (0.0, 0.0), 1.0);
    let bot_box = m2.rect("bot_box", (-2.0, -2.0), (2.0, -0.3));
    let sub = m2.intersect("sub_chord", &[disk, bot_box]);
    m2.overlay(&[("air", disk), ("sub", sub)]);
    m2.absorbing_layer("pml", PmlCoordKind::Cylindrical, 0.2, 6.0);
    let err2 = mesh_model(&m2.list, &MeshPlan::simple(0.35)).err().expect("must reject 2D circular PML cut by off-center horizontal substrate");
    assert!(err2.contains("Cylindrical PML") && err2.contains("Cartesian"), "err2 = {err2}");
}

#[test]
fn mesh_report_warns_under_5_layers_disables_structured_pml_relax_and_is_deterministic() {
    let mut m = Model::new_2d();
    let air = m.rect("air", (-1.0, -1.0), (1.0, 1.0));
    m.overlay(&[("air", air)]);
    m.pml_outward("pml", PmlCoordKind::Cartesian, &[Side::XMax, Side::ZMax], 0.3, 3.0, 2.0);

    let req = MeshRequest {
        h: 0.4,
        waves: Some(WavesRequest {
            sizing: Sizing {
                order: 2,
                headroom: 1.0,
                wavelengths: [("air*".to_string(), 1.2)].into_iter().collect(),
            },
            points_per_wavelength: Some(6.0), // wave_size = 1.2 / 6 = 0.2 m
            pml_relax: 3.0,                   // must be ignored (forced to 1.0) in structured mode!
            eigen_headroom: 1.0,
        }),
        ..MeshRequest::default()
    };
    let (_, _, _, _, rep, _) = mesh_request_native(&m.list, &req).unwrap();
    assert_eq!(rep.pml_warnings.len(), 1, "expected < 5 layers warning: {:?}", rep.pml_warnings);
    assert!(rep.pml_warnings[0].contains("< 5"));

    // Verify `pml_relax = 3.0` was NOT applied to `air.pml_xmax` (`wave_size` stays `0.2`, not `0.6`)!
    let pml_dom_rep = rep.domains.iter().find(|d| d.name == "air.pml_xmax").unwrap();
    assert!((pml_dom_rep.wave_size.unwrap() - 0.2).abs() < 1e-12, "{pml_dom_rep:?}");

    // Thread determinism (`1` vs `4` threads) on 3D curved cylindrical wafer
    let mut mc = Model::new_3d();
    let c1 = mc.cylinder("air", (0.0, 0.0, -0.4), 1.0, 0.8);
    mc.overlay(&[("air", c1)]);
    mc.absorbing_layer("pml", PmlCoordKind::Auto, 0.15, 5.0);
    let creq = MeshRequest { h: 0.45, curved: true, ..MeshRequest::default() };

    let run = |nt: usize| {
        rayon::ThreadPoolBuilder::new().num_threads(nt).build().unwrap().install(|| {
            let (Mesh::D3(cm), _, _, _, r, _) = mesh_request_native(&mc.list, &creq).unwrap() else {
                panic!("3D");
            };
            (cm.points, cm.conn, cm.edge_nodes, cm.quad_face_nodes, r)
        })
    };
    let a = run(1);
    let b = run(4);
    assert_eq!(a.0, b.0);
    assert_eq!(a.1, b.1);
    assert_eq!(a.2, b.2);
    assert_eq!(a.3, b.3);
    assert_eq!(a.4, b.4);
    let _ = PI;
}
