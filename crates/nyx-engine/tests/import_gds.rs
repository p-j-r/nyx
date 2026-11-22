// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end integration tests for `Engine::import_gds`: synthesizes a pocket-transmon GDSII
//! layout, imports conductive sheets with `TEXT` pins, extrudes `VacBot`/`Si`/`VacTop` boxes,
//! meshes via `sheet3d`, verifies all sheet & `<sheet>:edge` selectors, and solves electrostatics.

use geom::gds::{encode_ibm_f64, STRANS_REFLECT_X};
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

struct GdsBuilder {
    buf: Vec<u8>,
}

impl GdsBuilder {
    fn new(user_unit: f64, dbu_m: f64) -> Self {
        let mut b = Self { buf: Vec::new() };
        b.rec(0x00, 0x02, &600i16.to_be_bytes());
        b.rec(0x01, 0x02, &[0u8; 24]);
        b.rec_ascii(0x02, "transmon_test");
        let mut u = Vec::with_capacity(16);
        u.extend_from_slice(&encode_ibm_f64(user_unit));
        u.extend_from_slice(&encode_ibm_f64(dbu_m));
        b.rec(0x03, 0x05, &u);
        b
    }

    fn rec(&mut self, rec_type: u8, data_type: u8, payload: &[u8]) {
        let len = (4 + payload.len()) as u16;
        self.buf.extend_from_slice(&len.to_be_bytes());
        self.buf.push(rec_type);
        self.buf.push(data_type);
        self.buf.extend_from_slice(payload);
    }

    fn rec_ascii(&mut self, rec_type: u8, s: &str) {
        let mut p = s.as_bytes().to_vec();
        if !p.len().is_multiple_of(2) {
            p.push(0);
        }
        self.rec(rec_type, 0x06, &p);
    }

    fn rec_xy(&mut self, pts: &[[i32; 2]]) {
        let mut p = Vec::with_capacity(pts.len() * 8);
        for &[x, y] in pts {
            p.extend_from_slice(&x.to_be_bytes());
            p.extend_from_slice(&y.to_be_bytes());
        }
        self.rec(0x10, 0x03, &p);
    }

    fn begin_cell(&mut self, name: &str) {
        self.rec(0x05, 0x02, &[0u8; 24]);
        self.rec_ascii(0x06, name);
    }

    fn end_cell(&mut self) {
        self.rec(0x07, 0x00, &[]);
    }

    fn finish(mut self) -> Vec<u8> {
        self.rec(0x04, 0x00, &[]);
        self.buf
    }

    fn boundary(&mut self, layer: i16, datatype: i16, open_pts: &[[i32; 2]]) {
        self.rec(0x08, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x0E, 0x02, &datatype.to_be_bytes());
        let mut closed = open_pts.to_vec();
        if closed.first() != closed.last() {
            closed.push(open_pts[0]);
        }
        self.rec_xy(&closed);
        self.rec(0x11, 0x00, &[]);
    }

    fn path(&mut self, layer: i16, datatype: i16, pathtype: i16, width: i32, pts: &[[i32; 2]]) {
        self.rec(0x09, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x0E, 0x02, &datatype.to_be_bytes());
        self.rec(0x21, 0x02, &pathtype.to_be_bytes());
        self.rec(0x0F, 0x03, &width.to_be_bytes());
        self.rec_xy(pts);
        self.rec(0x11, 0x00, &[]);
    }

    fn sref(&mut self, sname: &str, xy: [i32; 2], reflect_x: bool) {
        self.rec(0x0A, 0x00, &[]);
        self.rec_ascii(0x12, sname);
        if reflect_x {
            self.rec(0x1A, 0x01, &STRANS_REFLECT_X.to_be_bytes());
        }
        self.rec_xy(&[xy]);
        self.rec(0x11, 0x00, &[]);
    }

    fn text(&mut self, layer: i16, texttype: i16, xy: [i32; 2], label: &str) {
        self.rec(0x0C, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x16, 0x02, &texttype.to_be_bytes());
        self.rec_xy(&[xy]);
        self.rec_ascii(0x19, label);
        self.rec(0x11, 0x00, &[]);
    }
}

fn synth_transmon_gds() -> Vec<u8> {
    // 1 nm DBU (1e-9 m)
    let mut g = GdsBuilder::new(1.0e-3, 1.0e-9);

    // Subcell "SHUNT_PAD": [-60..60] um x [10..70] um (area = 120 * 60 = 7200 um^2)
    g.begin_cell("SHUNT_PAD");
    g.boundary(1, 10, &[[-60_000, 10_000], [60_000, 10_000], [60_000, 70_000], [-60_000, 70_000]]);
    g.end_cell();

    // Top cell "TRANSMON_CHIP"
    g.begin_cell("TRANSMON_CHIP");

    // Ground plane on (1, 0): [-200..200] x [-200..200] um with central pocket hole [-100..100] x [-120..120] um
    // constructed via 4 overlapping ground bars:
    g.boundary(1, 0, &[[-200_000, -200_000], [200_000, -200_000], [200_000, -120_000], [-200_000, -120_000]]);
    g.boundary(1, 0, &[[-200_000, 120_000], [200_000, 120_000], [200_000, 200_000], [-200_000, 200_000]]);
    g.boundary(1, 0, &[[-200_000, -200_000], [-100_000, -200_000], [-100_000, 200_000], [-200_000, 200_000]]);
    g.boundary(1, 0, &[[100_000, -200_000], [200_000, -200_000], [200_000, 200_000], [100_000, 200_000]]);

    // Signal pads on (1, 10): top pad via unreflected SREF, bottom pad via X-reflected SREF!
    g.sref("SHUNT_PAD", [0, 0], false); // y in [10..70] um
    g.sref("SHUNT_PAD", [0, 0], true); // y in [-70..-10] um

    // Readout coupler T-claw on (1, 10) inside top of pocket (PATH + BOUNDARY boolean union!):
    // Horizontal bar [-50..50] x [85..95] um + vertical stem PATH (0, 90) -> (0, 110) um of width 10 um
    g.boundary(1, 10, &[[-50_000, 85_000], [50_000, 85_000], [50_000, 95_000], [-50_000, 95_000]]);
    g.path(1, 10, 0, 10_000, &[[0, 90_000], [0, 110_000]]);

    // Josephson junction strap on (20, 0) touching both shunt pads along y = -10 um and y = +10 um:
    g.boundary(20, 0, &[[-4_000, -10_000], [4_000, -10_000], [4_000, 10_000], [-4_000, 10_000]]);

    // TEXT pins on (1, 20):
    g.text(1, 20, [0, -40_000], "Q1_padA");
    g.text(1, 20, [0, 40_000], "Q1_padB");
    g.text(1, 20, [0, 90_000], "coupler");

    // Outer dicing frame marker on unrequested layer (70, 0) at [-300..300] um (tests layer_bbox vs cell_bbox)
    g.boundary(70, 0, &[[-300_000, -300_000], [300_000, -300_000], [300_000, -290_000], [-300_000, -290_000]]);

    g.end_cell();
    g.finish()
}

fn temp_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-import-gds-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn tri_area3(pts: &[f64], t: &[u32]) -> f64 {
    let p = |v: u32| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]];
    let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

#[test]
fn import_gds_builds_transmon_sheets_meshes_via_sheet3d_and_solves_capacitance() {
    let dir = temp_dir("transmon");
    let gds_path = dir.join("transmon.gds");
    std::fs::write(&gds_path, synth_transmon_gds()).unwrap();

    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Value {
        eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"))
    };

    call(&mut e, "project_new", json!({ "title": "GDS Transmon 3D", "dim": 3, "h": 50.0e-6 }));

    let imported = e
        .import_gds(&json!({
            "path": gds_path.to_str().unwrap(),
            "layers": {
                "gnd": [1, 0],
                "jj": [20, 0],
                "sig": [1, 10]
            },
            "z": 0.0,
            "split_islands": true,
            "pin_layer": [1, 20],
            "simplify_tol": 50.0e-9,
            "layer_bbox": true
        }))
        .unwrap();

    assert_eq!(imported["bbox"], json!([-200.0e-6, -200.0e-6, 200.0e-6, 200.0e-6]));
    assert_eq!(imported["layers"]["gnd"], json!(["gnd"]));
    assert_eq!(imported["layers"]["jj"], json!(["jj"]));
    assert_eq!(imported["layers"]["sig"], json!(["Q1_padA", "Q1_padB", "coupler"]));
    assert_eq!(
        imported["pins"],
        json!({
            "Q1_padA": "Q1_padA",
            "Q1_padB": "Q1_padB",
            "coupler": "coupler"
        })
    );

    let bbox = imported["bbox"].as_array().unwrap();
    let (xmin, ymin, xmax, ymax) =
        (bbox[0].as_f64().unwrap(), bbox[1].as_f64().unwrap(), bbox[2].as_f64().unwrap(), bbox[3].as_f64().unwrap());

    // Extrude dielectric stack from GDS bbox: VacBot [-250..-100] um, Si [-100..0] um, VacTop [0..150] um
    call(
        &mut e,
        "add_box",
        json!({ "name": "vac_bot", "min": [xmin, ymin, -250.0e-6], "max": [xmax, ymax, -100.0e-6] }),
    );
    call(&mut e, "add_box", json!({ "name": "si_sub", "min": [xmin, ymin, -100.0e-6], "max": [xmax, ymax, 0.0] }));
    call(&mut e, "add_box", json!({ "name": "vac_top", "min": [xmin, ymin, 0.0], "max": [xmax, ymax, 150.0e-6] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "VacBot", "input": "vac_bot" },
                { "name": "Si", "input": "si_sub" },
                { "name": "VacTop", "input": "vac_top" }
            ]
        }),
    );

    // Verify all sheet names and `<sheet>:edge` selectors appear in `geometry_names` and resolve in `add_size`!
    let names = call(&mut e, "geometry_names", json!({}));
    let bnds: Vec<&str> = names["boundaries"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
    for s in ["gnd", "jj", "Q1_padA", "Q1_padB", "coupler"] {
        assert!(bnds.contains(&s), "missing sheet boundary {s} in {bnds:?}");
        let edge_sel = format!("{s}:edge");
        assert!(bnds.contains(&edge_sel.as_str()), "missing {edge_sel} in {bnds:?}");
    }

    call(&mut e, "add_size", json!({ "sel": "jj:edge", "h": 4.0e-6 }));
    call(&mut e, "add_size", json!({ "sel": "Q1_padA:edge | Q1_padB:edge | coupler:edge", "h": 12.0e-6 }));
    call(&mut e, "mesh_build", json!({ "wait": true }));

    // Verify exact z=0 complement interface area in the 3D sheet3d mesh:
    // die_area = 400 x 400 um^2 = 160_000 um^2
    let model = e.model().expect("built mesh");
    let mesh = model.solver_mesh();
    let bnd_area = |bname: &str| -> f64 {
        let labs: Vec<i32> = model.boundary_names().into_iter().filter(|(_, n)| n == bname).map(|(id, _)| id).collect();
        assert!(!labs.is_empty(), "no mesh facets for {bname}");
        labs.into_iter().flat_map(|l| mesh.sheet_facets(l)).map(|tri| tri_area3(&mesh.points, tri)).sum()
    };
    let a_gnd = bnd_area("gnd"); // 160_000 - 200*240 = 112_000 um^2
    let a_pad_a = bnd_area("Q1_padA"); // 7_200 um^2
    let a_pad_b = bnd_area("Q1_padB"); // 7_200 um^2
    let a_coupler = bnd_area("coupler"); // 100*10 + 10*15 = 1_150 um^2
    let a_jj = bnd_area("jj"); // 8*20 = 160 um^2
    let a_iface = bnd_area("Si|VacTop");

    assert!((a_gnd - 112_000e-12).abs() < 1e-12 * 112_000e-12);
    assert!((a_pad_a - 7_200e-12).abs() < 1e-12 * 7_200e-12);
    assert!((a_pad_b - 7_200e-12).abs() < 1e-12 * 7_200e-12);
    assert!((a_coupler - 1_150e-12).abs() < 1e-12 * 1_150e-12);
    assert!((a_jj - 160e-12).abs() < 1e-12 * 160e-12);
    let total_z0 = a_iface + a_gnd + a_pad_a + a_pad_b + a_coupler + a_jj;
    let die_area = 160_000e-12;
    assert!((total_z0 - die_area).abs() < 1e-12 * die_area, "total_z0 = {total_z0:e} vs {die_area:e}");

    // Run electrostatic capacitance solve on the imported GDS sheets (with `jj` suppressed so pads stay independent):
    call(&mut e, "suppress_feature", json!({ "name": "jj", "suppressed": true }));
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacBot", "material": { "eps_r": 1.0 } }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "gnd", "kind": "ground" }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_top.zmax", "kind": "ground" }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_bot.zmin", "kind": "ground" }));
    for term in ["Q1_padA", "Q1_padB", "coupler"] {
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": term, "kind": "terminal", "params": { "name": term } }),
        );
    }
    call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true }));
    let cap = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    assert_eq!(cap["terminals"], json!(["Q1_padA", "Q1_padB", "coupler"]));
    let c00 = cap["maxwell_matrix_F"][0][0].as_f64().unwrap();
    let c01 = cap["maxwell_matrix_F"][0][1].as_f64().unwrap();
    assert!(c00 > 1.0e-15 && c01 < 0.0, "expected SPD Maxwell matrix: c00={c00:e}, c01={c01:e}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn import_gds_validates_dim3_path_guard_and_cell_bbox_mode() {
    let dir = temp_dir("guard");
    let gds_path = dir.join("chip.gds");
    std::fs::write(&gds_path, synth_transmon_gds()).unwrap();

    // 1. 2D project must be rejected cleanly without mutating features
    let mut e2 = Engine::new();
    let err_2d = e2
        .import_gds(&json!({
            "path": gds_path.to_str().unwrap(),
            "layers": { "gnd": [1, 0] }
        }))
        .unwrap_err()
        .to_string();
    assert!(err_2d.contains("3D") || err_2d.contains("dim: 3"), "unexpected 2D error: {err_2d}");
    assert!(e2.project().geometry.features.is_empty());

    // 2. Path guard confines relative paths inside root and rejects `..` escape
    let mut e3 = Engine::new();
    e3.call("project_new", &json!({ "dim": 3 })).unwrap();
    e3.set_path_root(Some(&dir)).unwrap();
    assert!(e3.import_gds(&json!({ "path": "../outside.gds", "layers": { "gnd": [1, 0] } })).is_err());

    // `layer_bbox: false` returns full cell bbox including dicing marker on (70, 0) at [-300..300] um,
    // and `split_islands: false` keeps all signal islands merged under "sig".
    let res = e3
        .import_gds(&json!({
            "path": "chip.gds",
            "layers": { "sig": [1, 10] },
            "split_islands": false,
            "layer_bbox": false
        }))
        .unwrap();
    assert_eq!(res["sheets"], json!(["sig"]));
    let bb = res["bbox"].as_array().unwrap();
    for (got, want) in bb.iter().zip([-300.0e-6, -300.0e-6, 300.0e-6, 200.0e-6]) {
        assert!((got.as_f64().unwrap() - want).abs() < 1e-18, "bbox mismatch: {bb:?}");
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn import_gds_lateral_pad_and_add_chip_stack_build_padded_collar_and_solve_farfield() {
    let dir = temp_dir("collar");
    let gds_path = dir.join("transmon.gds");
    std::fs::write(&gds_path, synth_transmon_gds()).unwrap();

    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "Padded Chip Stack", "dim": 3, "h": 80.0e-6 })).unwrap();

    // Negative lateral_pad must be rejected cleanly
    assert!(e
        .import_gds(&json!({
            "path": gds_path.to_str().unwrap(),
            "layers": { "sig": [1, 10] },
            "lateral_pad": -10.0e-6
        }))
        .is_err());

    let pad = 150.0e-6;
    let imported = e
        .import_gds(&json!({
            "path": gds_path.to_str().unwrap(),
            "layers": {
                "gnd": [1, 0],
                "sig": [1, 10]
            },
            "pin_layer": [1, 20],
            "lateral_pad": pad
        }))
        .unwrap();

    assert_eq!(imported["unpadded_bbox"], json!([-200.0e-6, -200.0e-6, 200.0e-6, 200.0e-6]));
    assert_eq!(imported["padded_bbox"], json!([-350.0e-6, -350.0e-6, 350.0e-6, 350.0e-6]));
    assert_eq!(imported["bbox"], json!([-350.0e-6, -350.0e-6, 350.0e-6, 350.0e-6]));
    assert_eq!(imported["lateral_pad"], json!(pad));

    // Build padded (VacBot, Si, VacTop) stack via Engine::add_chip_stack (inferring unpadded bbox + lateral_pad)
    assert!(e.add_chip_stack(&json!({ "lateral_pad": -1.0 })).is_err());
    let stack = e
        .add_chip_stack(&json!({
            "t_sub": 100.0e-6,
            "h_vac_top": 150.0e-6,
            "h_vac_bot": 150.0e-6,
            "lateral_pad": pad
        }))
        .unwrap();
    assert_eq!(stack["unpadded_bbox"], json!([-200.0e-6, -200.0e-6, 200.0e-6, 200.0e-6]));
    assert_eq!(stack["padded_bbox"], json!([-350.0e-6, -350.0e-6, 350.0e-6, 350.0e-6]));
    assert_eq!(stack["regions"], json!(["VacBot", "Si", "VacTop"]));

    e.call("mesh_build", &json!({ "wait": true })).unwrap();
    e.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 1 })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "VacBot", "material": { "eps_r": 1.0 } })).unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "gnd", "kind": "ground" })).unwrap();
    e.call(
        "add_condition",
        &json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 1, "decay_power": 1.0 } }),
    )
    .unwrap();
    for term in ["Q1_padA", "Q1_padB"] {
        e.call(
            "add_condition",
            &json!({ "physics": "dfs", "on": term, "kind": "terminal", "params": { "name": term } }),
        )
        .unwrap();
    }
    e.call("solve_run", &json!({ "physics": "dfs", "wait": true })).unwrap();
    let cap = e.call("analyse", &json!({ "analysis": "capacitance_matrix", "physics": "dfs" })).unwrap();
    assert_eq!(cap["terminals"], json!(["Q1_padA", "Q1_padB"]));
    let c00 = cap["maxwell_matrix_F"][0][0].as_f64().unwrap();
    assert!(c00 > 1.0e-15, "c00 = {c00:e}");

    let _ = std::fs::remove_dir_all(&dir);
}
