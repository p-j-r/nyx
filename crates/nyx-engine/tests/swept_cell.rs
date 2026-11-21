// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A periodic cell made by sweeping a 2D section: the section is meshed, its vertices are renumbered by coordinates so that opposite walls
//! are cut alike, and layers with a `fill` region take that region everywhere (cladding above and below a patterned slab). Floquet pairs
//! then join the walls. Checked on a vacuum cell with a closed form, and by the refusal of a fill that names no region.

use nyx_engine::Engine;
use serde_json::{json, Value};

const C0: f64 = 299_792_458.0;

/// Section x in [0, lx], z in [0, lz], with a rectangular `Core` region inside a `Host` region (both vacuum), periodic along x and z.
fn project(fill: &str, lx: f64, lz: f64, ly: f64, h: f64) -> Value {
    let feats = json!([
        {"id": 1, "name": "cell", "op": {"Rect": {"min": [0.0, 0.0], "max": [lx, lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 2, "name": "core", "op": {"Rect": {"min": [0.3 * lx, 0.3 * lz], "max": [0.6 * lx, 0.7 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 3, "name": "overlay", "op": {"Overlay": {"naming": "Feature", "pins": [], "regions": [
            {"input": 1, "kind": "Material", "name": "Host", "pinned_id": null},
            {"input": 2, "kind": "Material", "name": "Core", "pinned_id": null}]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 4, "name": "pz", "op": {"Periodic": {"sel_a": {"Name": "cell.bottom"}, "sel_b": {"Name": "cell.top"}, "translation": [0.0, 0.0, lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 5, "name": "px", "op": {"Periodic": {"sel_a": {"Name": "cell.left"}, "sel_b": {"Name": "cell.right"}, "translation": [lx, 0.0, 0.0]}}, "origin": "Code", "suppressed": false, "tags": []}
    ]);
    json!({
        "active": "ems", "format": 1, "title": "swept periodic cell",
        "geometry": {"dim": 2, "features": feats, "format": 1, "params": []},
        "mesh": {"h": h, "sweep": [{"thickness": ly / 4.0, "n": 1, "fill": fill}, {"thickness": ly / 2.0, "n": 2}, {"thickness": ly / 4.0, "n": 1, "fill": fill}]},
        "physics": [{"physics": "ems", "format": 1, "order": 2,
            "bcs": [
                {"kind": "periodic", "on": "name(cell.bottom)", "params": {"dest": "name(cell.top)", "phase": 0}},
                {"kind": "periodic", "on": "name(cell.left)", "params": {"dest": "name(cell.right)", "phase": 0}}],
            "materials": {"Host": {"eps_r": 1.0}, "Core": {"eps_r": 1.0}},
            "study": {"kind": "cavity", "params": {"f_search": 5.0e9, "n_modes": 4}}}]
    })
}

#[test]
fn a_swept_periodic_cell_with_filled_cladding_layers_gives_the_vacuum_modes() {
    let (lx, lz, ly) = (0.04, 0.03, 0.025);
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": project("Host", lx, lz, ly, 0.01).to_string()})).unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    // the Core region exists only in the middle layer
    let model = e.model().unwrap().clone();
    let core = model.regions().iter().find(|r| r.1 == "Core").unwrap().0;
    let m = model.solver_mesh();
    let n_core = m.cell_labels.iter().filter(|&&l| l == core).count();
    assert!(n_core > 0, "the unfilled middle layer keeps the core");
    let (lo, hi) = (ly / 4.0, 0.75 * ly);
    for c in 0..m.n_cells() {
        if m.cell_labels[c] == core {
            let y: Vec<f64> = m.cell(c).iter().map(|&v| m.point(v as usize)[1]).collect();
            assert!(
                y.iter().all(|&y| y >= lo - 1e-12 && y <= hi + 1e-12),
                "a core cell outside the middle layer: {y:?}"
            );
        }
    }
    e.call("solve_run", &json!({"physics": "ems"})).unwrap();
    let (res, _) = e.result("ems").unwrap();
    let f: Vec<f64> = res.solutions.iter().filter_map(|s| s.value).collect();
    // Gamma point of the vacuum cell, PEC ends along y: two polarisations at c / (2 ly)
    let exact = C0 / (2.0 * ly);
    let near = f.iter().filter(|&&v| (v / exact - 1.0).abs() < 5e-3).count();
    assert!(near >= 2, "expected a degenerate pair at {exact:.4e} Hz, got {f:?}");
}

#[test]
fn a_fill_region_that_does_not_exist_is_refused_by_name() {
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": project("Nope", 0.04, 0.03, 0.025, 0.01).to_string()})).unwrap();
    let err = e.call("mesh_build", &json!({})).unwrap_err().to_string();
    assert!(err.contains("Nope") && err.contains("Host"), "{err}");
}

#[test]
fn a_swept_2d_section_remaps_regions_per_layer_and_refuses_unknown_names() {
    let (lx, lz, ly) = (1.0e-6, 0.8e-6, 0.4e-6);
    let feats = json!([
        {"id": 1, "name": "clad", "op": {"Rect": {"min": [0.0, 0.0], "max": [lx, lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 2, "name": "core", "op": {"Rect": {"min": [0.25 * lx, 0.3 * lz], "max": [0.55 * lx, 0.7 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 3, "name": "teeth", "op": {"Rect": {"min": [0.55 * lx, 0.3 * lz], "max": [0.75 * lx, 0.7 * lz]}}, "origin": "Code", "suppressed": false, "tags": []},
        {"id": 4, "name": "overlay", "op": {"Overlay": {"naming": "Feature", "pins": [], "regions": [
            {"input": 1, "kind": "Material", "name": "Host", "pinned_id": null},
            {"input": 2, "kind": "Material", "name": "Core", "pinned_id": null},
            {"input": 3, "kind": "Material", "name": "Teeth", "pinned_id": null}]}}, "origin": "Code", "suppressed": false, "tags": []}
    ]);
    let proj = json!({
        "active": null, "format": 1, "title": "2D per-layer region remap",
        "geometry": {"dim": 2, "features": feats, "format": 1, "params": []},
        "mesh": {
            "h": 0.2e-6,
            "sweep": [
                {"thickness": 0.5 * ly, "n": 2},
                {"thickness": 0.5 * ly, "n": 2, "regions": {"Teeth": "Host"}}
            ]
        },
        "physics": []
    });
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": proj.to_string()})).unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    let model = e.model().unwrap().clone();
    let core_id = model.regions().iter().find(|r| r.1 == "Core").unwrap().0;
    let teeth_id = model.regions().iter().find(|r| r.1 == "Teeth").unwrap().0;
    let m = model.solver_mesh();
    let mut core_first_half = 0;
    let mut core_second_half = 0;
    let mut teeth_first_half = 0;
    let mut teeth_second_half = 0;
    for c in 0..m.n_cells() {
        let ym: f64 = m.cell(c).iter().map(|&v| m.point(v as usize)[1]).sum::<f64>() / 4.0;
        if m.cell_labels[c] == core_id {
            if ym < 0.5 * ly {
                core_first_half += 1;
            } else {
                core_second_half += 1;
            }
        } else if m.cell_labels[c] == teeth_id {
            if ym < 0.5 * ly {
                teeth_first_half += 1;
            } else {
                teeth_second_half += 1;
            }
        }
    }
    assert!(core_first_half > 0 && core_second_half > 0, "Core persists in both halves");
    assert!(teeth_first_half > 0, "Teeth exists in the first half");
    assert_eq!(teeth_second_half, 0, "Teeth is remapped to Host in the second half");

    // Unknown region name in `regions` remap is refused clearly
    let mut bad = proj;
    bad["mesh"]["sweep"][1]["regions"] = json!({"Teeth": "MissingRegion"});
    e.call("project_load", &json!({"json": bad.to_string()})).unwrap();
    let err = e.call("mesh_build", &json!({})).unwrap_err().to_string();
    assert!(err.contains("MissingRegion"), "{err}");
}

#[test]
fn a_3d_corrugated_grating_swept_by_layers_keeps_3d_selectors_and_congruent_periodic_ends() {
    use femcore::periodic::match_vertices;
    use femcore::props::Periodic;
    use geom::model::Model;
    use geom::select::Selector;
    use nyx_physics::LabelResolver;
    use std::collections::BTreeSet;

    let lambda = 0.44e-6;
    let mut gm = Model::new_3d();
    let clad = gm.box_("clad", (-1.0e-6, 0.0, -0.8e-6), (1.0e-6, lambda, 0.8e-6));
    let core = gm.box_("core", (-0.225e-6, 0.0, -0.11e-6), (0.225e-6, lambda, 0.11e-6));
    let teeth = gm.box_("teeth", (0.225e-6, 0.0, -0.11e-6), (0.325e-6, 0.5 * lambda, 0.11e-6));
    gm.overlay(&[("SiO2", clad), ("Si_Core", core), ("Si_Teeth", teeth)]);
    gm.size(
        "si_size",
        Selector::parse("region(Si_Core) | region(Si_Teeth)").unwrap(),
        Some(0.09e-6.into()),
        None,
        None,
    );
    gm.periodic(
        "py",
        Selector::parse("outer & plane(y = 0)").unwrap(),
        Selector::parse(&format!("outer & plane(y = {lambda:e})")).unwrap(),
        (0.0, lambda, 0.0),
    );
    let proj = json!({
        "active": null, "format": 1, "title": "3D swept corrugated grating",
        "geometry": serde_json::to_value(&gm.list).unwrap(),
        "mesh": {
            "h": 0.25e-6,
            "sweep": [
                {"thickness": 0.5 * lambda, "n": 2},
                {"thickness": 0.5 * lambda, "n": 2}
            ]
        },
        "physics": []
    });
    let mut e = Engine::new();
    e.call("project_load", &json!({"json": proj.to_string()})).unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    let model = e.model().unwrap().clone();
    let core_id = model.regions().iter().find(|r| r.1 == "Si_Core").unwrap().0;
    let teeth_id = model.regions().iter().find(|r| r.1 == "Si_Teeth").unwrap().0;
    let m = model.solver_mesh();
    assert_eq!(m.dim, 3);

    let mut teeth_in_first = 0;
    let mut teeth_in_second = 0;
    let mut core_in_first = 0;
    let mut core_in_second = 0;
    for c in 0..m.n_cells() {
        let ym: f64 = m.cell(c).iter().map(|&v| m.point(v as usize)[1]).sum::<f64>() / 4.0;
        if m.cell_labels[c] == teeth_id {
            if ym < 0.5 * lambda {
                teeth_in_first += 1;
            } else {
                teeth_in_second += 1;
            }
        } else if m.cell_labels[c] == core_id {
            if ym < 0.5 * lambda {
                core_in_first += 1;
            } else {
                core_in_second += 1;
            }
        }
    }
    assert!(core_in_first > 0 && core_in_second > 0, "Si_Core exists in both layers");
    assert!(teeth_in_first > 0 && teeth_in_second == 0, "Si_Teeth exists only in first half-period");

    // 3D selectors and ymin/ymax resolve cleanly and match periodic faces
    let src = model.boundaries("outer & plane(y = 0)").unwrap();
    let dst = model.boundaries(&format!("outer & plane(y = {lambda:e})")).unwrap();
    let iface = model.boundaries("Si_Core|SiO2").unwrap();
    assert!(!src.is_empty() && !dst.is_empty() && !iface.is_empty());
    assert_eq!(model.boundaries("ymin").unwrap(), src);
    assert_eq!(model.boundaries("ymax").unwrap(), dst);

    let matched = match_vertices(&m, &Periodic { source: src, dest: dst.clone(), phase: 0.0 })
        .expect("3D swept periodic ends must match vertex for vertex");
    let dst_tris: BTreeSet<[u32; 3]> = (0..m.n_facets())
        .filter(|&f| dst.contains(&m.facet_labels[f]))
        .map(|f| {
            let mut t = [m.facet(f)[0], m.facet(f)[1], m.facet(f)[2]];
            t.sort_unstable();
            t
        })
        .collect();
    for &fc in &matched.src_facets {
        let v = m.facet(fc);
        let mut vt = [matched.map[&v[0]], matched.map[&v[1]], matched.map[&v[2]]];
        vt.sort_unstable();
        assert!(dst_tris.contains(&vt), "triangle {vt:?} missing on destination periodic face");
    }
}
