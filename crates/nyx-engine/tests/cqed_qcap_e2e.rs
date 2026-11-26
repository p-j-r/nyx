// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end integration tests for superconducting qubit capacitance (`QFR-1..QFR-7`) through `Engine::call`:
//! * QFR-1: `set_mesh` `gap_cells` coplanar gap auto-sizing and `z_coarsen` pure-tetrahedral vertical coarsening.
//! * QFR-2: Two-level P1-iso-P2 + Thomas vertical-line AMG (`amg_line_smooth`) and multi-RHS Block-CG (`solver: "iterative"`).
//! * QFR-3: Barsoum `singular_edges: "quarter_point"` P2 enrichment + `convergence_check` `extrapolate: "richardson_log_spd"`
//!   with `apply_extrapolated: true` promoting the extrapolated matrix onto the stored `SolveResult`.
//! * QFR-4: 2D layered-dielectric boundary-element `formulation: "planar_bem"` (`richardson: true`) cross-checked against 3D FEM.
//! * QFR-5: Binary GDSII stream import via `Engine::call("import_gds", ...)` with hierarchy, `PATH`, pocket hole, and `TEXT` pins.
//! * QFR-6: `analyse("capacitance_matrix")` ground/Schur elimination, differential 2-pad `TransmonPocket` reduction, and `g_Hz` / `chi_Hz`.
//! * QFR-7: `get_guide` `"cqed"` topic and `project_script` round-trip of `gap_cells` and `z_coarsen`.

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
        b.rec_ascii(0x02, "cqed_e2e_lib");
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

fn synth_pocket_transmon_gds() -> Vec<u8> {
    let mut g = GdsBuilder::new(1.0e-3, 1.0e-9);

    // Shunt pad subcell: [-65..65] um x [12..68] um (24 um inter-pad gap)
    g.begin_cell("PAD");
    g.boundary(1, 10, &[[-65_000, 12_000], [65_000, 12_000], [65_000, 68_000], [-65_000, 68_000]]);
    g.end_cell();

    g.begin_cell("TOP");
    // Ground plane [-220..220] x [-220..220] um with central pocket hole [-110..110] x [-130..130] um
    g.boundary(1, 0, &[[-220_000, -220_000], [220_000, -220_000], [220_000, -130_000], [-220_000, -130_000]]);
    g.boundary(1, 0, &[[-220_000, 130_000], [220_000, 130_000], [220_000, 220_000], [-220_000, 220_000]]);
    g.boundary(1, 0, &[[-220_000, -220_000], [-110_000, -220_000], [-110_000, 220_000], [-220_000, 220_000]]);
    g.boundary(1, 0, &[[110_000, -220_000], [220_000, -220_000], [220_000, 220_000], [110_000, 220_000]]);

    // Two differential pads via SREF + X-reflection
    g.sref("PAD", [0, 0], false); // Q1_padB in y > 0
    g.sref("PAD", [0, 0], true); // Q1_padA in y < 0

    // Readout T-claw near Q1_padB on (1, 10)
    g.boundary(1, 10, &[[-55_000, 86_000], [55_000, 86_000], [55_000, 98_000], [-55_000, 98_000]]);
    g.path(1, 10, 0, 14_000, &[[0, 92_000], [0, 118_000]]);

    // TEXT pins on (1, 20)
    g.text(1, 20, [0, -40_000], "Q1_padA");
    g.text(1, 20, [0, 40_000], "Q1_padB");
    g.text(1, 20, [0, 92_000], "RO_claw");

    g.end_cell();
    g.finish()
}

fn temp_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-cqed-qcap-e2e-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn mat_f64(v: &Value) -> Vec<Vec<f64>> {
    v.as_array().unwrap().iter().map(|r| r.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()).collect()
}

#[test]
fn cqed_gds_sheet_qcap_solvers_extrapolation_promotion_and_reduction_e2e() {
    let dir = temp_dir("pocket");
    let gds_path = dir.join("pocket_transmon.gds");
    std::fs::write(&gds_path, synth_pocket_transmon_gds()).unwrap();

    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: Value| -> Value {
        eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"))
    };

    // QFR-7: Verify `get_guide` exposes the `"cqed"` card
    let guide = call(&mut e, "get_guide", json!({ "topic": "cqed" }));
    assert_eq!(guide["topic"], json!("cqed"));
    let guide_text = guide["text"].as_str().unwrap();
    assert!(guide_text.len() > 500 && guide_text.len() < 3600);
    assert!(guide_text.contains("import_gds") && guide_text.contains("capacitance_matrix"));

    // QFR-5: Import synthetic binary GDSII via `Engine::call("import_gds", ...)`
    call(&mut e, "project_new", json!({ "title": "cQED Pocket Transmon E2E", "dim": 3, "h": 75.0e-6 }));
    let imported = call(
        &mut e,
        "import_gds",
        json!({
            "path": gds_path.to_str().unwrap(),
            "layers": {
                "gnd": [1, 0],
                "sig": [1, 10]
            },
            "z": 0.0,
            "split_islands": true,
            "pin_layer": [1, 20],
            "simplify_tol": 50.0e-9,
            "layer_bbox": true
        }),
    );
    assert_eq!(imported["layers"]["gnd"], json!(["gnd"]));
    assert_eq!(imported["layers"]["sig"], json!(["Q1_padA", "Q1_padB", "RO_claw"]));
    assert_eq!(imported["pins"]["Q1_padA"], json!("Q1_padA"));
    assert_eq!(imported["pins"]["Q1_padB"], json!("Q1_padB"));
    assert_eq!(imported["pins"]["RO_claw"], json!("RO_claw"));

    let bb = imported["bbox"].as_array().unwrap();
    let (xmin, ymin, xmax, ymax) =
        (bb[0].as_f64().unwrap(), bb[1].as_f64().unwrap(), bb[2].as_f64().unwrap(), bb[3].as_f64().unwrap());

    // Extrude layered Si substrate [-180..0] um with VacBot [-320..-180] um and VacTop [0..220] um
    call(
        &mut e,
        "add_box",
        json!({ "name": "vac_bot", "min": [xmin, ymin, -320.0e-6], "max": [xmax, ymax, -180.0e-6] }),
    );
    call(&mut e, "add_box", json!({ "name": "si_sub", "min": [xmin, ymin, -180.0e-6], "max": [xmax, ymax, 0.0] }));
    call(&mut e, "add_box", json!({ "name": "vac_top", "min": [xmin, ymin, 0.0], "max": [xmax, ymax, 220.0e-6] }));
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

    // QFR-1: Verify `gap_cells` + `z_coarsen` mesh settings reduce cell count while round-tripping in `project_script`
    call(&mut e, "set_mesh", json!({ "h": 65.0e-6, "gap_cells": 2.2, "z_coarsen": false }));
    let mesh_uncoarsened = call(&mut e, "mesh_build", json!({ "wait": true }));
    let cells_uncoarsened = mesh_uncoarsened["mesh"]["cells"].as_u64().unwrap();

    call(&mut e, "set_mesh", json!({ "h": 65.0e-6, "gap_cells": 2.2, "z_coarsen": true }));
    let script_out = call(&mut e, "project_script", json!({}));
    assert_eq!(script_out["form"], json!("commands"));
    let script_text = script_out["script"].as_str().unwrap();
    assert!(script_text.contains("gap_cells") && script_text.contains("z_coarsen"), "{script_text}");

    let mesh_coarsened = call(&mut e, "mesh_build", json!({ "wait": true }));
    let cells_coarsened = mesh_coarsened["mesh"]["cells"].as_u64().unwrap();
    assert!(
        cells_coarsened < cells_uncoarsened,
        "expected z_coarsen=true ({cells_coarsened}) to have fewer 3D cells than z_coarsen=false ({cells_uncoarsened})"
    );

    // Configure 4-terminal electrostatic study (`gnd`, `Q1_padA`, `Q1_padB`, `RO_claw`) + top/bottom box ground
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "dfs",
            "study": "electrostatic",
            "order": 2,
            "study_params": {
                "formulation": "planar_bem",
                "richardson": true
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacBot", "material": { "eps_r": 1.0 } }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_top.zmax", "kind": "ground" }));
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_bot.zmin", "kind": "ground" }));
    for term in ["gnd", "Q1_padA", "Q1_padB", "RO_claw"] {
        call(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": term, "kind": "terminal", "params": { "name": term } }),
        );
    }

    // QFR-4: Solve via `formulation = "planar_bem"` (`richardson = true`) first
    call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true }));
    let cap_bem = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    assert_eq!(cap_bem["terminals"], json!(["gnd", "Q1_padA", "Q1_padB", "RO_claw"]));
    let m_bem = mat_f64(&cap_bem["maxwell_matrix_F"]);
    assert!(m_bem[1][1] > 10.0e-15 && m_bem[1][2] < -1.0e-15, "unexpected planar_bem matrix: {m_bem:?}");

    // QFR-2 & QFR-3: Solve via 3D `volumetric_fem` with `solver: "iterative"` (Block-CG + line-AMG) + `singular_edges: "quarter_point"`
    call(
        &mut e,
        "set_study",
        json!({
            "physics": "dfs",
            "kind": "electrostatic",
            "order": 2,
            "params": {
                "formulation": "volumetric_fem",
                "singular_edges": "quarter_point"
            }
        }),
    );
    call(
        &mut e,
        "set_physics_params",
        json!({
            "physics": "dfs",
            "params": {
                "solver": "iterative",
                "tol": 1.0e-11,
                "amg_line_smooth": true
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "physics": "dfs", "wait": true }));
    let cap_fem_before = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    let m_before = mat_f64(&cap_fem_before["maxwell_matrix_F"]);
    assert!(m_before[1][1] > 10.0e-15 && m_before[1][2] < -1.0e-15);

    // Cross-check inter-pad mutual coupling order-of-magnitude consistency between planar_bem and enclosed-box 3D FEM
    let mut_bem_12 = -m_bem[1][2];
    let mut_fem_12 = -m_before[1][2];
    let rel_12 = (mut_bem_12 - mut_fem_12).abs() / mut_fem_12;
    assert!(
        rel_12 < 0.20,
        "inter-pad coupling C_12 mismatch between planar_bem ({mut_bem_12:e}) and 3D FEM ({mut_fem_12:e}): rel={rel_12:.3}"
    );

    // QFR-3: Run two-grid `convergence_check` with coarsening `factor: 2.1` (achieving effective ratio r_eff >= 1.3 per R11b), `extrapolate: "richardson_log_spd"`, and `apply_extrapolated: true`
    let conv = call(
        &mut e,
        "convergence_check",
        json!({
            "physics": "dfs",
            "factor": 2.1,
            "scale_sizes": true,
            "tolerance": 0.50,
            "extrapolate": "richardson_log_spd",
            "apply_extrapolated": true
        }),
    );
    assert_eq!(conv["extrapolated"]["method"], json!("richardson_log_spd"));
    assert_eq!(conv["extrapolated"]["order"], json!(1.0));
    assert_eq!(conv["extrapolated"]["apply_extrapolated"], json!(true));
    assert_eq!(conv["validated"]["state"], json!("checked"), "convergence verdict: {conv}");
    let m_ext = mat_f64(&conv["extrapolated"]["capacitance_matrix"]["maxwell_matrix_f"]);

    // Verify `apply_extrapolated: true` promoted the Richardson-extrapolated matrix onto the stored `SolveResult`!
    let cap_after = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    assert_eq!(cap_after["stale"], json!(false));
    assert_eq!(cap_after["validated"]["state"], json!("checked"));
    let m_after = mat_f64(&cap_after["maxwell_matrix_F"]);
    for i in 0..4 {
        for j in 0..4 {
            assert!(
                (m_after[i][j] - m_ext[i][j]).abs() <= 1e-27,
                "stored capacitance matrix entry [{i}][{j}] was not promoted: got {:e}, expected {:e}",
                m_after[i][j],
                m_ext[i][j]
            );
        }
    }
    assert!(
        (m_after[1][1] - m_before[1][1]).abs() > 1e-20,
        "extrapolated matrix should refine raw single-grid capacitance"
    );

    // QFR-6: Run cQED ground projection, differential 2-pad TransmonPocket reduction, and quarter-wave resonator coupling
    let cqed = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "capacitance_matrix",
            "physics": "dfs",
            "params": {
                "ground": "gnd",
                "eliminate": "drop",
                "qubits": [
                    { "name": "Q1", "pads": ["Q1_padA", "Q1_padB"], "E_J_Hz": 16.0e9 }
                ],
                "resonators": [
                    { "name": "RO", "terminal": "RO_claw", "f_r_Hz": 6.8e9, "Z0_ohm": 50.0, "resonator_type": "quarter_wave" }
                ]
            }
        }),
    );
    assert_eq!(cqed["reduced_terminals"], json!(["Q1_padA", "Q1_padB", "RO_claw"]));
    assert_eq!(cqed["has_null_mode"], json!(false));
    assert_eq!(cqed["null_modes"], json!(0));

    let q1 = &cqed["qubits"][0];
    let c_sigma = q1["C_Sigma_F"].as_f64().unwrap();
    let ec_hz = q1["E_C_Hz"].as_f64().unwrap();
    let f01_hz = q1["f01_Hz"].as_f64().unwrap();
    let alpha_hz = q1["anharmonicity_Hz"].as_f64().unwrap();
    let beta_ro = q1["beta"]["RO_claw"].as_f64().unwrap();
    assert!(c_sigma > 15.0e-15 && c_sigma < 200.0e-15, "unexpected C_Sigma_F = {c_sigma:e}");
    assert!(ec_hz > 50.0e6 && ec_hz < 1.5e9, "unexpected E_C_Hz = {ec_hz:e}");
    assert!(f01_hz > 2.0e9 && f01_hz < 12.0e9, "unexpected f01_Hz = {f01_hz:e}");
    assert!(alpha_hz < 0.0);
    assert!(beta_ro > 1.0e-3 && beta_ro < 0.5, "unexpected coupling divider beta = {beta_ro}");

    let coup = &cqed["couplings"][0];
    let g_hz = coup["g_Hz"].as_f64().unwrap();
    let chi_hz = coup["chi_Hz"].as_f64().unwrap();
    assert!(g_hz > 1.0e6 && g_hz < 500.0e6, "unexpected vacuum Rabi coupling g_Hz = {g_hz:e}");
    assert!(chi_hz.is_finite() && chi_hz.abs() > 1.0e3, "unexpected dispersive shift chi_Hz = {chi_hz:e}");

    // Also verify floating Schur elimination on the same 4-terminal matrix
    let schur = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "capacitance_matrix",
            "physics": "dfs",
            "params": {
                "ground": "gnd",
                "eliminate": "floating_schur"
            }
        }),
    );
    assert_eq!(schur["eliminate"], json!("floating_schur"));
    assert_eq!(schur["reduced_terminals"], json!(["Q1_padA", "Q1_padB", "RO_claw"]));

    let _ = std::fs::remove_dir_all(&dir);
}
