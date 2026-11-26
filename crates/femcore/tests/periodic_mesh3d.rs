// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Periodic parent->child 3D lattice and Rivara refinement conformance tests (FR-2):
//! verifies vertex-for-vertex (< 1e-13 * extent), edge-for-edge, and triangle-for-triangle
//! congruence on 3D corrugated waveguide gratings with local region sizing across 20 pitches
//! and thread-count determinism (`NYX_THREADS=1` vs `4`).

use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts, MeshedGeom};
use femcore::meshing::geomdrive3::mesh_geom_3d;
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::periodic::match_vertices;
use femcore::props::Periodic;
use geom::model::Model;
use geom::select::Selector;
use std::collections::BTreeSet;

fn pt(m: &SimplexMesh, v: u32) -> [f64; 3] {
    [m.points[3 * v as usize], m.points[3 * v as usize + 1], m.points[3 * v as usize + 2]]
}

fn sorted3(mut t: [u32; 3]) -> [u32; 3] {
    t.sort_unstable();
    t
}

/// Builds a 3D corrugated strip/rib grating unit cell of pitch `lambda`:
/// - `SiO2` cladding box over `[-1.0 um, 1.0 um] x [0, lambda] x [-0.8 um, 0.8 um]`
/// - Continuous `Si_Core` over `[-0.225 um, 0.225 um] x [0, lambda] x [-0.11 um, 0.11 um]`
/// - Asymmetric grating tooth `Si_Teeth` over `[0.225 um, 0.325 um] x [0.18 * lambda, 0.68 * lambda] x [-0.11 um, 0.11 um]`
/// - Local `Size` rule `h_si` on `region(Si_Core) | region(Si_Teeth)` (`Si*`)
/// - Periodic pair on `plane(y = 0)` -> `plane(y = lambda)` with translation `[0, lambda, 0]`.
fn corrugated_grating(lambda: f64, h_si: f64) -> geom::model::FeatureList {
    let mut m = Model::new_3d();
    let clad = m.box_("clad", (-1.0e-6, 0.0, -0.8e-6), (1.0e-6, lambda, 0.8e-6));
    let core = m.box_("core", (-0.225e-6, 0.0, -0.11e-6), (0.225e-6, lambda, 0.11e-6));
    let teeth = m.box_(
        "teeth",
        (0.225e-6, 0.18 * lambda, -0.11e-6),
        (0.325e-6, 0.68 * lambda, 0.11e-6),
    );
    m.overlay(&[("SiO2", clad), ("Si_Core", core), ("Si_Teeth", teeth)]);
    m.size(
        "si_size",
        Selector::parse("region(Si_Core) | region(Si_Teeth)").unwrap(),
        Some(h_si.into()),
        None,
        None,
    );
    let src = Selector::parse("outer & plane(y = 0)").unwrap();
    let dst = Selector::parse(&format!("outer & plane(y = {lambda:e})")).unwrap();
    m.periodic("py", src, dst, (0.0, lambda, 0.0));
    m.list
}

/// Verifies vertex-for-vertex (< 1e-13 * extent), edge-for-edge, and triangle-for-triangle congruence
/// between the `source` and `dest` boundaries of a periodic pair, including the exact `face_id[&sorted3(vt)]`
/// lookup required by `nyx-ems::cavity::ned_relations` (order 2).
fn assert_periodic_congruence(g: &MeshedGeom, src_sel: &str, dst_sel: &str, expected_shift: [f64; 3]) {
    let m = &g.mesh;
    let source = g.select(src_sel).unwrap();
    let dest = g.select(dst_sel).unwrap();
    assert!(!source.is_empty() && !dest.is_empty());

    let (lo, hi) = m.bbox();
    let extent = (0..3).map(|k| hi[k] - lo[k]).fold(0.0f64, f64::max).max(1e-300);

    let per = Periodic { source: source.clone(), dest: dest.clone(), phase: 0.0 };
    let matched = match_vertices(m, &per).expect("match_vertices must succeed on periodic 3D mesh");

    for k in 0..3 {
        assert!(
            (matched.shift[k] - expected_shift[k]).abs() <= 1e-13 * extent,
            "shift mismatch on axis {k}: got {:?}, expected {:?}",
            matched.shift,
            expected_shift
        );
    }

    // 1. Vertex-for-vertex drift < 1e-13 * extent
    let mut max_drift = 0.0f64;
    for (&sv, &dv) in &matched.map {
        let ps = pt(m, sv);
        let pd = pt(m, dv);
        let drift = ((pd[0] - (ps[0] + expected_shift[0])).powi(2)
            + (pd[1] - (ps[1] + expected_shift[1])).powi(2)
            + (pd[2] - (ps[2] + expected_shift[2])).powi(2))
        .sqrt();
        max_drift = max_drift.max(drift);
    }
    assert!(
        max_drift < 1e-13 * extent,
        "periodic vertex drift {max_drift:.3e} exceeds 1e-13 * extent ({:.3e})",
        1e-13 * extent
    );

    // 2. Collect all cell faces (as `Ned3::new` does in `nyx-ems::cavity`) and destination boundary facets
    let mut all_cell_faces: BTreeSet<[u32; 3]> = BTreeSet::new();
    for ci in 0..m.n_cells() {
        let v = m.cell(ci);
        for f in [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] {
            all_cell_faces.insert(sorted3([v[f[0]], v[f[1]], v[f[2]]]));
        }
    }
    let dst_facets: BTreeSet<[u32; 3]> = (0..m.n_facets())
        .filter(|&f| dest.contains(&m.facet_labels[f]))
        .map(|f| {
            let v = m.facet(f);
            sorted3([v[0], v[1], v[2]])
        })
        .collect();
    assert_eq!(
        matched.src_facets.len(),
        dst_facets.len(),
        "source and destination periodic boundaries must have the same number of triangle facets"
    );

    // 3. Edge-for-edge and triangle-for-triangle congruence
    for &fc in &matched.src_facets {
        let v = m.facet(fc);
        let vs = [v[0], v[1], v[2]];
        let vt = [matched.map[&vs[0]], matched.map[&vs[1]], matched.map[&vs[2]]];
        for (i, j) in [(0, 1), (1, 2), (0, 2)] {
            let edge_key = (vt[i].min(vt[j]), vt[i].max(vt[j]));
            assert!(
                m.edges().binary_search(&edge_key).is_ok(),
                "source edge ({}, {}) mapped to ({}, {}) is missing on destination boundary (diagonal mismatch)",
                vs[i],
                vs[j],
                vt[i],
                vt[j]
            );
        }
        let tri_key = sorted3(vt);
        assert!(
            dst_facets.contains(&tri_key),
            "source facet {vs:?} mapped to {vt:?} is missing from destination facets"
        );
        assert!(
            all_cell_faces.contains(&tri_key),
            "source facet {vs:?} mapped to {vt:?} is missing from cell face table (Ned3::face_dofs would panic)"
        );
    }
}

#[test]
fn corrugated_grating_periodic_faces_match_across_20_pitches() {
    let h_glob = 0.30e-6;
    let h_si = 0.09e-6;
    for idx in 0..20 {
        let lambda = 200e-9 + (idx as f64) * (750e-9 / 19.0);
        let list = corrugated_grating(lambda, h_si);
        let g = mesh_model(&list, &MeshPlan::simple(h_glob)).unwrap();
        assert_periodic_congruence(
            &g,
            "outer & plane(y = 0)",
            &format!("outer & plane(y = {lambda:e})"),
            [0.0, lambda, 0.0],
        );
    }
}

#[test]
fn periodic_boundary_rivara_bisection_and_adaptive_3d_keep_congruence() {
    // 1. Fine local sizing (h_si = 45 nm < lattice h_axis = 75 nm) forces Rivara bisections directly on the periodic faces
    let lambda = 370e-9;
    let coarse = mesh_model(&corrugated_grating(lambda, 0.25e-6), &MeshPlan::simple(0.30e-6)).unwrap();
    let fine = mesh_model(&corrugated_grating(lambda, 0.045e-6), &MeshPlan::simple(0.30e-6)).unwrap();
    let si_core = fine.domain("Si_Core").unwrap();
    let count_core = |g: &MeshedGeom| g.mesh.cell_labels.iter().filter(|&&l| l == si_core).count();
    assert!(
        count_core(&fine) > 4 * count_core(&coarse) && fine.mesh.n_cells() > coarse.mesh.n_cells(),
        "core cells {} -> {}, total cells {} -> {}",
        count_core(&coarse),
        count_core(&fine),
        coarse.mesh.n_cells(),
        fine.mesh.n_cells()
    );
    assert_periodic_congruence(
        &fine,
        "outer & plane(y = 0)",
        &format!("outer & plane(y = {lambda:e})"),
        [0.0, lambda, 0.0],
    );

    // 2. Adaptive 3D meshing freezes periodic boundary faces (`periodic_labels` includes `EntityKind::Faces`)
    let adapted = mesh_model(&corrugated_grating(lambda, 0.12e-6), &MeshPlan::adaptive(0.35e-6, 0.7)).unwrap();
    assert_periodic_congruence(
        &adapted,
        "outer & plane(y = 0)",
        &format!("outer & plane(y = {lambda:e})"),
        [0.0, lambda, 0.0],
    );
}

fn push_mesh_bits(buf: &mut Vec<u8>, m: &SimplexMesh) {
    let push_usize = |b: &mut Vec<u8>, v: usize| b.extend_from_slice(&(v as u64).to_le_bytes());
    push_usize(buf, m.points.len());
    for &p in &m.points {
        buf.extend_from_slice(&p.to_bits().to_le_bytes());
    }
    push_usize(buf, m.cells.len());
    for &c in &m.cells {
        buf.extend_from_slice(&c.to_le_bytes());
    }
    push_usize(buf, m.cell_labels.len());
    for &l in &m.cell_labels {
        buf.extend_from_slice(&l.to_le_bytes());
    }
    push_usize(buf, m.facets.len());
    for &f in &m.facets {
        buf.extend_from_slice(&f.to_le_bytes());
    }
    push_usize(buf, m.facet_labels.len());
    for &l in &m.facet_labels {
        buf.extend_from_slice(&l.to_le_bytes());
    }
    let edge_nodes = m.edge_node_list();
    push_usize(buf, edge_nodes.len());
    for (pair, pos) in &edge_nodes {
        buf.extend_from_slice(&pair[0].to_le_bytes());
        buf.extend_from_slice(&pair[1].to_le_bytes());
        for &x in pos {
            buf.extend_from_slice(&x.to_bits().to_le_bytes());
        }
    }
}

fn push_periodic_match_bits(buf: &mut Vec<u8>, g: &MeshedGeom, src_sel: &str, dst_sel: &str) {
    let m = &g.mesh;
    let per = Periodic { source: g.select(src_sel).unwrap(), dest: g.select(dst_sel).unwrap(), phase: 0.0 };
    let matched = match_vertices(m, &per).unwrap();
    for s in matched.shift {
        buf.extend_from_slice(&s.to_bits().to_le_bytes());
    }
    buf.extend_from_slice(&(matched.map.len() as u64).to_le_bytes());
    let mut pairs: Vec<(u32, u32)> = matched.map.iter().map(|(&sv, &dv)| (sv, dv)).collect();
    pairs.sort_unstable();
    for (sv, dv) in pairs {
        buf.extend_from_slice(&sv.to_le_bytes());
        buf.extend_from_slice(&dv.to_le_bytes());
    }
}

fn curved_periodic_cell(lambda: f64) -> geom::model::FeatureList {
    let mut m = Model::new_3d();
    let clad_2d = m.rect("clad_2d", (-0.8e-6, -0.6e-6), (0.8e-6, 0.6e-6));
    let clad = m.extrude("clad", clad_2d, 0.0, lambda);
    let cyl_2d = m.circle("cyl_2d", (0.0, 0.0), 0.22e-6);
    let cyl = m.extrude("cyl", cyl_2d, 0.0, lambda);
    m.overlay(&[("SiO2", clad), ("Si_Cyl", cyl)]);
    let src = Selector::parse("outer & plane(y = 0)").unwrap();
    let dst = Selector::parse(&format!("outer & plane(y = {lambda:e})")).unwrap();
    m.periodic("py", src, dst, (0.0, lambda, 0.0));
    m.list
}

#[test]
fn periodic_3d_meshing_is_deterministic_across_thread_counts() {
    let lambda = 437e-9;
    if let Ok(out_path) = std::env::var("NYX_PERIODIC3D_OUT") {
        let mut buf = Vec::new();

        // 1. Direct `mesh_geom_3d` lattice periodic corrugated grating
        let list_lat = corrugated_grating(lambda, 0.085e-6);
        let ev_lat = evaluate(&list_lat).unwrap();
        let g_lat = mesh_geom_3d(ev_lat, &GeomMeshOpts::new(0.28e-6)).unwrap();
        push_mesh_bits(&mut buf, &g_lat.mesh);
        push_periodic_match_bits(&mut buf, &g_lat, "outer & plane(y = 0)", &format!("outer & plane(y = {lambda:e})"));

        // 2. Rivara bisection on periodic boundaries via `mesh_model(..., Simple)`
        let g_riv = mesh_model(&corrugated_grating(lambda, 0.065e-6), &MeshPlan::simple(0.30e-6)).unwrap();
        push_mesh_bits(&mut buf, &g_riv.mesh);
        push_periodic_match_bits(&mut buf, &g_riv, "outer & plane(y = 0)", &format!("outer & plane(y = {lambda:e})"));

        // 3. Adaptive 3D periodic meshing
        let g_adp = mesh_model(&corrugated_grating(lambda, 0.12e-6), &MeshPlan::adaptive(0.35e-6, 0.7)).unwrap();
        push_mesh_bits(&mut buf, &g_adp.mesh);

        // 4. Curved periodic cylinder cell exercising `curve3` parallel edge projection + periodic twin sync
        let ev_curv = evaluate(&curved_periodic_cell(lambda)).unwrap();
        let mut copts = GeomMeshOpts::new(0.24e-6);
        copts.curved = true;
        let g_curv = mesh_geom_3d(ev_curv, &copts).unwrap();
        assert!(g_curv.mesh.has_curved(), "expected curved edges on periodic cylinder cell");
        push_mesh_bits(&mut buf, &g_curv.mesh);
        push_periodic_match_bits(&mut buf, &g_curv, "outer & plane(y = 0)", &format!("outer & plane(y = {lambda:e})"));

        std::fs::write(out_path, buf).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |nt: &str| -> Vec<u8> {
        let out_file = std::env::temp_dir().join(format!("nyx-periodic3d-det-{}-t{nt}.bin", std::process::id()));
        let status = std::process::Command::new(&exe)
            .args(["--exact", "periodic_3d_meshing_is_deterministic_across_thread_counts", "--nocapture"])
            .env("NYX_THREADS", nt)
            .env("NYX_PERIODIC3D_OUT", &out_file)
            .status()
            .unwrap();
        assert!(status.success(), "child process NYX_THREADS={nt} failed");
        let bytes = std::fs::read(&out_file).unwrap();
        let _ = std::fs::remove_file(&out_file);
        bytes
    };

    let b1 = run_with_threads("1");
    let b4 = run_with_threads("4");
    assert!(!b1.is_empty());
    assert_eq!(b1, b4, "periodic 3D mesh must be 0-ULP bit-identical for NYX_THREADS=1 and 4");
}
