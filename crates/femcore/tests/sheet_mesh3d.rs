// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for the conformal 3D sheet-stack mesher (`femcore::meshing::sheet3d`):
//! - Touching `pad_1`, `jj_1`, `pad_2` sheets have 0 T-junctions and share exact boundary vertices in the 3D mesh;
//! - Complement interface subtraction (`SA`): `area("substrate|vacuum") + sum(area(sheets)) == die_area` to `< 1e-12`;
//! - Interior sheet triangles appear **exactly once** in `SimplexMesh::facets` / `sheet_facets`, shared by lower and upper tetrahedra;
//! - Sheet mesh refinement respects `"<sheet>:edge"` `add_size` rules and grades z-slabs near sheet planes.

use femcore::meshing::request::{mesh_request, MeshRequest, MesherRequest};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use std::collections::{HashMap, HashSet};

fn tri_area3(pts: &[f64], t: &[u32]) -> f64 {
    let p = |v: u32| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]];
    let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

fn build_transmon_sheets_model(edge_h: Option<f64>) -> Model {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
    let vac = m.box_("vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));

    // pad_1: [-80..-10] in y, [-60..60] in x, with a 10 um x 10 um flux-trap hole Kept
    // and a 1 um x 1 um noise hole Dropped by min_hole_area = 4e-12
    let pad1_outer = vec![pt(-60e-6, -80e-6), pt(60e-6, -80e-6), pt(60e-6, -10e-6), pt(-60e-6, -10e-6)];
    let hole_small = vec![pt(-30e-6, -50e-6), pt(-29e-6, -50e-6), pt(-29e-6, -49e-6), pt(-30e-6, -49e-6)];
    let hole_flux = vec![pt(20e-6, -50e-6), pt(30e-6, -50e-6), pt(30e-6, -40e-6), pt(20e-6, -40e-6)];
    let p1 = m.polygons_with("pad_1_poly", vec![pad1_outer], vec![hole_small, hole_flux], 50e-9, 4e-12);

    // jj_1 touches pad_1 along y = -10e-6 (x in [-4e-6, 4e-6]) and pad_2 along y = 10e-6 (x in [-4e-6, 4e-6])
    let jj_outer = vec![pt(-4e-6, -10e-6), pt(4e-6, -10e-6), pt(4e-6, 10e-6), pt(-4e-6, 10e-6)];
    let jj = m.polygons("jj_1_poly", vec![jj_outer], vec![]);

    // pad_2: [10..80] in y, [-60..60] in x
    let pad2_outer = vec![pt(-60e-6, 10e-6), pt(60e-6, 10e-6), pt(60e-6, 80e-6), pt(-60e-6, 80e-6)];
    let p2 = m.polygons("pad_2_poly", vec![pad2_outer], vec![]);

    m.sheet("pad_1", p1, 0.0);
    m.sheet("jj_1", jj, 0.0);
    m.sheet("pad_2", p2, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    if let Some(h) = edge_h {
        m.size("jj_edge_size", Selector::parse("jj_1:edge").unwrap(), Some(Expr::Num(h)), None, None);
    }
    m
}

#[test]
fn touching_sheets_have_zero_t_junctions_and_exact_complement_subtraction_in_3d_mesh() {
    let m = build_transmon_sheets_model(Some(4e-6));
    let g = mesh_model(&m.list, &MeshPlan::simple(40e-6)).unwrap();
    let mesh = &g.mesh;

    // 1. Complement subtraction on the 3D mesh triangles at z = 0:
    // area("substrate|vacuum") + area("pad_1") + area("jj_1") + area("pad_2") == die_area (400 um x 400 um)
    let die_area = 400e-6 * 400e-6;
    let bnd_area = |name: &str| -> f64 {
        let labs = g.boundary(name);
        assert!(!labs.is_empty(), "missing boundary label for {name}");
        labs.into_iter()
            .flat_map(|l| mesh.sheet_facets(l))
            .map(|tri| tri_area3(&mesh.points, tri))
            .sum()
    };
    let a_iface = bnd_area("substrate|vacuum");
    let a_pad1 = bnd_area("pad_1");
    let a_jj1 = bnd_area("jj_1");
    let a_pad2 = bnd_area("pad_2");

    // Analytic sheet areas:
    // pad_1 = 120x70 um^2 - 10x10 um^2 (flux hole) = 8300 um^2
    // jj_1  = 8x20 um^2 = 160 um^2
    // pad_2 = 120x70 um^2 = 8400 um^2
    assert!((a_pad1 - 8300e-12).abs() < 1e-12 * 8300e-12, "a_pad1 = {a_pad1:e}");
    assert!((a_jj1 - 160e-12).abs() < 1e-12 * 160e-12, "a_jj1 = {a_jj1:e}");
    assert!((a_pad2 - 8400e-12).abs() < 1e-12 * 8400e-12, "a_pad2 = {a_pad2:e}");
    let total_z0 = a_iface + a_pad1 + a_jj1 + a_pad2;
    assert!(
        (total_z0 - die_area).abs() < 1e-12 * die_area,
        "mesh complement subtraction mismatch: {total_z0:e} vs {die_area:e}"
    );

    // 2. Every interior sheet facet appears EXACTLY ONCE in mesh.facets and is shared by
    //    exactly 1 substrate tet (below z=0) and 1 vacuum tet (above z=0).
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();
    let mut tet_faces: HashMap<[u32; 3], Vec<i32>> = HashMap::new();
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let dom = mesh.cell_labels[c];
        for skip in 0..4 {
            let mut f = [0u32; 3];
            let mut k = 0;
            for (i, &vi) in v.iter().enumerate() {
                if i != skip {
                    f[k] = vi;
                    k += 1;
                }
            }
            f.sort_unstable();
            tet_faces.entry(f).or_default().push(dom);
        }
    }

    let mut seen_sheet_tris: HashSet<[u32; 3]> = HashSet::new();
    for sheet_name in ["pad_1", "jj_1", "pad_2"] {
        for lab in g.boundary(sheet_name) {
            let tris = mesh.sheet_facets(lab);
            assert!(!tris.is_empty());
            for tri in tris {
                let mut key = [tri[0], tri[1], tri[2]];
                key.sort_unstable();
                assert!(seen_sheet_tris.insert(key), "sheet triangle {key:?} emitted more than once!");
                let doms = tet_faces.get(&key).expect("sheet triangle must belong to 3D tetrahedra");
                assert_eq!(doms.len(), 2, "interior sheet triangle must be shared by 2 tetrahedra");
                assert!(doms.contains(&sub_dom) && doms.contains(&vac_dom));
            }
        }
    }

    // 3. Verify 0 T-junctions and exact shared boundary vertices on the contact segments
    //    between pad_1 & jj_1 (y = -10e-6, x in [-4e-6, 4e-6]) and jj_1 & pad_2 (y = +10e-6, x in [-4e-6, 4e-6]).
    let edges_of_sheet = |name: &str| -> HashSet<(u32, u32)> {
        let mut set = HashSet::new();
        for lab in g.boundary(name) {
            for tri in mesh.sheet_facets(lab) {
                for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                    let e = (a.min(b), a.max(b));
                    if !set.insert(e) {
                        // Internal edge shared by 2 triangles of the same sheet: remove to leave only the 1D boundary
                        set.remove(&e);
                    }
                }
            }
        }
        set
    };
    let bnd_pad1 = edges_of_sheet("pad_1");
    let bnd_jj1 = edges_of_sheet("jj_1");
    let bnd_pad2 = edges_of_sheet("pad_2");

    let shared_1_jj: Vec<(u32, u32)> = bnd_pad1.intersection(&bnd_jj1).copied().collect();
    let shared_2_jj: Vec<(u32, u32)> = bnd_pad2.intersection(&bnd_jj1).copied().collect();
    assert!(!shared_1_jj.is_empty(), "pad_1 and jj_1 must share exact boundary edges");
    assert!(!shared_2_jj.is_empty(), "pad_2 and jj_1 must share exact boundary edges");

    let edge_len = |(a, b): (u32, u32)| -> f64 {
        let pa = mesh.point(a as usize);
        let pb = mesh.point(b as usize);
        ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2)).sqrt()
    };
    let len_1_jj: f64 = shared_1_jj.iter().copied().map(edge_len).sum();
    let len_2_jj: f64 = shared_2_jj.iter().copied().map(edge_len).sum();
    // Exact contact width is 8 um on both contacts:
    assert!((len_1_jj - 8e-6).abs() < 1e-12, "shared pad_1/jj_1 contact length = {len_1_jj:e} vs 8e-6");
    assert!((len_2_jj - 8e-6).abs() < 1e-12, "shared pad_2/jj_1 contact length = {len_2_jj:e} vs 8e-6");

    // Check global 2D plane conformity at z = 0: every edge of any z=0 triangle (sheets + substrate|vacuum)
    // is shared by exactly 2 triangles at z = 0 (or 1 if on the outer die perimeter x,y = +-200 um).
    let mut z0_edge_count: HashMap<(u32, u32), usize> = HashMap::new();
    for name in ["substrate|vacuum", "pad_1", "jj_1", "pad_2"] {
        for lab in g.boundary(name) {
            for tri in mesh.sheet_facets(lab) {
                for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                    *z0_edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
                }
            }
        }
    }
    for (&(a, b), &cnt) in &z0_edge_count {
        let pa = mesh.point(a as usize);
        let pb = mesh.point(b as usize);
        let on_outer_rim = (pa[0].abs() - 200e-6).abs() < 1e-9 && (pb[0].abs() - 200e-6).abs() < 1e-9
            || (pa[1].abs() - 200e-6).abs() < 1e-9 && (pb[1].abs() - 200e-6).abs() < 1e-9;
        if on_outer_rim {
            assert_eq!(cnt, 1, "outer rim edge {a}-{b} should belong to 1 z=0 triangle");
        } else {
            assert_eq!(cnt, 2, "interior z=0 edge {a}-{b} has count {cnt} (T-junction detected!)");
        }
    }
}

#[test]
fn sheet_edge_size_rule_refines_perimeter_and_grades_z_slabs_without_warnings() {
    let coarse_m = build_transmon_sheets_model(None);
    let fine_m = build_transmon_sheets_model(Some(2e-6));

    let req = MeshRequest {
        h: 40e-6,
        max_dev: None,
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
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
    };
    let (g_coarse, rep_coarse) = mesh_request(&coarse_m.list, &req).unwrap();
    let (g_fine, rep_fine) = mesh_request(&fine_m.list, &req).unwrap();

    assert!(rep_coarse.warnings.is_empty(), "unexpected warnings: {:?}", rep_coarse.warnings);
    assert!(rep_fine.warnings.is_empty(), "edge size rule on 3D sheet must not emit crease warning: {:?}", rep_fine.warnings);

    let jj_tris_coarse: usize = g_coarse.boundary("jj_1").into_iter().map(|l| g_coarse.mesh.sheet_facets(l).len()).sum();
    let jj_tris_fine: usize = g_fine.boundary("jj_1").into_iter().map(|l| g_fine.mesh.sheet_facets(l).len()).sum();
    assert!(
        jj_tris_fine > jj_tris_coarse * 4,
        "jj_1:edge size rule must refine jj_1 facets (coarse = {jj_tris_coarse}, fine = {jj_tris_fine})"
    );
    assert!(g_fine.mesh.n_cells() > g_coarse.mesh.n_cells());
}

#[test]
fn suspended_airbridge_sheet_preserved_under_structured_pml_for_both_z_coarsen_modes() {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    // Full box [-200..200] x [-200..200] x [-100..120] um; PML window [-150..150] x [-150..150] x [-60..80] um (thickness 40..50? uniform 40 um!)
    let sub = m.box_("sub", (-190e-6, -190e-6, -100e-6), (190e-6, 190e-6, 0.0));
    let vac = m.box_("vac", (-190e-6, -190e-6, 0.0), (190e-6, 190e-6, 120e-6));
    let win = m.box_("window", (-150e-6, -150e-6, -60e-6), (150e-6, 150e-6, 80e-6));

    // Coplanar CPW signal trace at z = 0 (between Substrate and Vacuum: la != lb)
    let cpw_sig = m.polygons("cpw_poly", vec![vec![pt(-8e-6, -120e-6), pt(8e-6, -120e-6), pt(8e-6, 120e-6), pt(-8e-6, 120e-6)]], vec![]);
    m.sheet("cpw_sig", cpw_sig, 0.0);

    // Suspended airbridge deck at z = +15 um strictly inside Vacuum (`la == lb`!)
    let bridge_poly = m.polygons("bridge_poly", vec![vec![pt(-40e-6, -12e-6), pt(40e-6, -12e-6), pt(40e-6, 12e-6), pt(-40e-6, 12e-6)]], vec![]);
    m.sheet("airbridge", bridge_poly, 15e-6);

    m.overlay(&[("substrate", sub), ("vacuum", vac)]);
    m.pml("pml", win, &[], 40e-6, 2.0, 2.0);
    m.size("sz_bridge", Selector::parse("name(airbridge)").unwrap(), Some(Expr::Num(10e-6)), None, None);

    let vac_expected_deck_area = 80e-6 * 24e-6; // 1920 um^2
    let cpw_expected_area = 16e-6 * 240e-6; // 3840 um^2
    let total_vol_expected = 380e-6 * 380e-6 * 220e-6;

    for z_coarsen in [false, true] {
        let req = MeshRequest {
            h: 40e-6,
            max_dev: None,
            mesher: MesherRequest::Simple,
            pml_mode: femcore::meshing::request::PmlModeRequest::Structured,
            structured_pml: None,
            waves: None,
            curved: false,
            gap_cells: None,
            z_coarsen,
            max_cells: None,
            dz_min: None,
            max_z_slabs: None,
            z_growth_ratio: None,
            max_steiner_multiplier: None,
            strict_mesh: true,
            singular_edges: None,
            quarter_point_xi: None,
            corner_grading: false,
        };
        let (g, rep) = mesh_request(&m.list, &req).unwrap_or_else(|e| panic!("z_coarsen={z_coarsen}: {e}"));
        assert!(rep.warnings.is_empty(), "z_coarsen={z_coarsen}: {:?}", rep.warnings);
        assert_eq!(rep.structured_layers, Some(vec![[2, 2], [2, 2], [2, 2]]));

        let vac_dom = g.domain("vacuum").unwrap();
        let sub_dom = g.domain("substrate").unwrap();

        // Build tet face incidence map and check strictly positive tet volumes
        let mut total_vol = 0.0;
        let mut tet_faces: HashMap<[u32; 3], Vec<i32>> = HashMap::new();
        for c in 0..g.mesh.n_cells() {
            let v = g.mesh.cell(c);
            let p = |i: u32| [g.mesh.points[3 * i as usize], g.mesh.points[3 * i as usize + 1], g.mesh.points[3 * i as usize + 2]];
            let (p0, p1, p2, p3) = (p(v[0]), p(v[1]), p(v[2]), p(v[3]));
            let (u, w, q) = (
                [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]],
                [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]],
                [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]],
            );
            let det = u[0] * (w[1] * q[2] - w[2] * q[1]) - u[1] * (w[0] * q[2] - w[2] * q[0]) + u[2] * (w[0] * q[1] - w[1] * q[0]);
            assert!(det > 0.0, "z_coarsen={z_coarsen}: non-positive tet {c} det={det:e}");
            total_vol += det / 6.0;
            for skip in 0..4 {
                let mut f = [0u32; 3];
                let mut k = 0;
                for (i, &vi) in v.iter().enumerate() {
                    if i != skip {
                        f[k] = vi;
                        k += 1;
                    }
                }
                f.sort_unstable();
                tet_faces.entry(f).or_default().push(g.mesh.cell_labels[c]);
            }
        }
        assert!(
            (total_vol - total_vol_expected).abs() < 1e-10 * total_vol_expected,
            "z_coarsen={z_coarsen}: vol {total_vol:e} vs {total_vol_expected:e}"
        );

        // Verify suspended `airbridge` (`la == lb == vac_dom`) facets survived `structured_pml = true` intact
        let bridge_labs = g.boundary("airbridge");
        assert!(!bridge_labs.is_empty(), "missing airbridge boundary label");
        let mut bridge_area = 0.0;
        let mut bridge_seen: HashSet<[u32; 3]> = HashSet::new();
        for lab in bridge_labs {
            let tris = g.mesh.sheet_facets(lab);
            assert!(!tris.is_empty(), "z_coarsen={z_coarsen}: airbridge facets lost under structured_pml=true!");
            for tri in tris {
                let mut key = [tri[0], tri[1], tri[2]];
                key.sort_unstable();
                assert!(bridge_seen.insert(key), "z_coarsen={z_coarsen}: duplicate airbridge facet {key:?}");
                let doms = tet_faces.get(&key).expect("airbridge facet must belong to 3D tetrahedra");
                assert_eq!(doms.as_slice(), &[vac_dom, vac_dom], "z_coarsen={z_coarsen}: airbridge must separate vacuum|vacuum (la == lb)");
                // Check +z orientation and exact z = 15 um elevation
                let p = |i: u32| [g.mesh.points[3 * i as usize], g.mesh.points[3 * i as usize + 1], g.mesh.points[3 * i as usize + 2]];
                let (p0, p1, p2) = (p(tri[0]), p(tri[1]), p(tri[2]));
                for pv in [p0, p1, p2] {
                    assert!((pv[2] - 15e-6).abs() < 1e-13, "airbridge vertex z = {:e}", pv[2]);
                }
                let nz = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);
                assert!(nz > 0.0, "airbridge facet must keep +z normal orientation");
                bridge_area += tri_area3(&g.mesh.points, tri);
            }
        }
        assert!(
            (bridge_area - vac_expected_deck_area).abs() < 1e-12 * vac_expected_deck_area,
            "z_coarsen={z_coarsen}: airbridge area {bridge_area:e} vs {vac_expected_deck_area:e}"
        );

        // Also verify coplanar `cpw_sig` (`la != lb`) at z = 0 survived with exact area
        let cpw_area: f64 = g
            .boundary("cpw_sig")
            .into_iter()
            .flat_map(|l| g.mesh.sheet_facets(l))
            .map(|tri| {
                let mut key = [tri[0], tri[1], tri[2]];
                key.sort_unstable();
                let doms = &tet_faces[&key];
                assert_eq!(doms.len(), 2);
                assert!(doms.contains(&sub_dom) && doms.contains(&vac_dom));
                tri_area3(&g.mesh.points, tri)
            })
            .sum();
        assert!(
            (cpw_area - cpw_expected_area).abs() < 1e-12 * cpw_expected_area,
            "z_coarsen={z_coarsen}: cpw_sig area {cpw_area:e} vs {cpw_expected_area:e}"
        );
    }
}

#[test]
fn circular_disk_sheet3d_curved_mid_edge_projection_populates_valid_curved_tets() {
    let r_disk = 80e-6;
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
    let vac = m.box_("vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));
    let disk_2d = m.circle("disk_2d", (0.0, 0.0), r_disk);
    m.sheet("disk", disk_2d, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    // Coarse chord discretization on purpose (`max_dev = 2.0e-6`) so chords sag visibly inside `r_disk`
    let straight_req = MeshRequest {
        h: 45e-6,
        max_dev: Some(2.0e-6),
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: false,
        gap_cells: None,
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: true,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };
    let (g_straight, rep_straight) = mesh_request(&m.list, &straight_req).unwrap();
    assert!(!g_straight.mesh.has_curved());
    assert!(rep_straight.curved.is_none());

    for z_coarsen in [false, true] {
        let curved_req = MeshRequest { curved: true, z_coarsen, ..straight_req.clone() };
        let (g_curved, rep_curved) = mesh_request(&m.list, &curved_req).unwrap();
        assert!(g_curved.mesh.has_curved(), "z_coarsen={z_coarsen}: expected curved edge nodes");
        let cs = rep_curved.curved.expect("MeshReport::curved must be populated");
        assert!(cs.edges > 0 && cs.cells > 0);
        assert!(cs.min_scaled_jacobian > 0.1, "z_coarsen={z_coarsen}: min_scaled_jacobian = {}", cs.min_scaled_jacobian);

        // When z_coarsen == false, straight vertex/cell arrays must remain bit-identical to `curved = false`
        if !z_coarsen {
            assert_eq!(g_curved.mesh.points, g_straight.mesh.points);
            assert_eq!(g_curved.mesh.cells, g_straight.mesh.cells);
            assert_eq!(g_curved.mesh.facets, g_straight.mesh.facets);
        }

        // Verify every 1D perimeter edge of `disk` at z = 0 carries a mid-edge node on the exact circle r = 80 um
        let mut disk_rim: HashSet<(u32, u32)> = HashSet::new();
        for lab in g_curved.boundary("disk") {
            for tri in g_curved.mesh.sheet_facets(lab) {
                for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                    let e = (a.min(b), a.max(b));
                    if !disk_rim.insert(e) {
                        disk_rim.remove(&e);
                    }
                }
            }
        }
        assert!(!disk_rim.is_empty());
        let edge_map: HashMap<[u32; 2], [f64; 3]> = g_curved.mesh.edge_node_list().into_iter().collect();
        for (a, b) in disk_rim {
            let pa = g_curved.mesh.point(a as usize);
            let pb = g_curved.mesh.point(b as usize);
            let chord_mid_r = (0.5 * (pa[0] + pb[0])).hypot(0.5 * (pa[1] + pb[1]));
            assert!(chord_mid_r < r_disk - 1e-7, "expected noticeable chord sagitta, got {chord_mid_r:e}");
            let q = edge_map.get(&[a, b]).unwrap_or_else(|| panic!("missing mid-edge node on disk rim edge ({a}, {b})"));
            let qr = q[0].hypot(q[1]);
            assert!((qr - r_disk).abs() < 1e-12 * r_disk, "mid-edge node radius {qr:e} != {r_disk:e}");
            assert!(q[2].abs() < 1e-15, "mid-edge node z = {:e}", q[2]);
        }

        // Verify all 3D cells have strictly positive scaled Jacobian and thread-count determinism
        for c in 0..g_curved.mesh.n_cells() {
            let sj = femcore::meshing::curve3::scaled_jacobian_min(&g_curved.mesh, c);
            assert!(sj > 0.0, "tet {c} scaled_jacobian_min = {sj}");
        }
        let m1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap().install(|| {
            mesh_request(&m.list, &curved_req).unwrap().0.mesh
        });
        let m8 = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap().install(|| {
            mesh_request(&m.list, &curved_req).unwrap().0.mesh
        });
        assert_eq!(m1.edge_node_list(), m8.edge_node_list());
    }
}
