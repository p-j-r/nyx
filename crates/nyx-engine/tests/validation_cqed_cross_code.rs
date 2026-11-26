// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 cross-code cQED validation suite:
//! - AWS Palace / analytical coaxial capacitor (`[C_low, C_up]` Prager-Synge + GCI)
//! - Binary GDSII `TransmonPocket` cross-code capacitance matrix (3D volumetric FEM with Prager-Synge
//!   flux bounds vs independent 2D layered Green-function `planar_bem`) + SQuADDS Q3D reduction
//! - 3D lossy rectangular cavity (`f_0, Q_0`) + 2-port `CircuitRom` QEP vs `ems::driven` Taubin `S21` fit

mod common;

use common::{
    cmd, inline_target, load_validation_json, u_num_two_grid, verify_asme_vv20, verify_epistemic_vv20, ValidationTarget,
};
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::rt::{equilibrate_electro_flux_bounds, CellPermittivity, RobinBoundary};
use geom::gds::encode_ibm_f64;
use geom::model::{Model, RegionKind};
use nyx_analysis::capacitance::capacitance_matrix_with_params;
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study, EPS0};
use nyx_engine::Engine;
use nyx_physics::{CapacitanceMatrixData, SolveResult};
use serde::Deserialize;
use serde_json::{json, Value};
use std::f64::consts::PI;
use std::path::PathBuf;

#[derive(Deserialize)]
struct CrossCodeDataset {
    coax_cylinder_electrostatic: CoaxSection,
    planar_gds_pocket_cross_code: PlanarPocketSection,
    squadds_q3d_reduction: SquaddsSection,
    cavity_3d_eigen_and_driven: CavitySection,
}

#[derive(Deserialize)]
struct CoaxGeom {
    a_m: f64,
    b_m: f64,
    length_m: f64,
    eps_r: f64,
}

#[derive(Deserialize)]
struct CoaxSection {
    geometry: CoaxGeom,
    target: ValidationTarget,
}

#[derive(Deserialize)]
struct PlanarPocketGeom {
    pad_half_x_m: f64,
    pad_y_lo_m: f64,
    pad_y_hi_m: f64,
    claw_x_m: f64,
    claw_y_lo_m: f64,
    claw_y_hi_m: f64,
    pocket_x_m: f64,
    pocket_y_m: f64,
    gnd_outer_m: f64,
    t_sub_m: f64,
    eps_r_si: f64,
}

#[derive(Deserialize)]
struct PlanarPocketSection {
    geometry: PlanarPocketGeom,
    tol_rel: f64,
    u_input_rel: f64,
    u_d_rel: f64,
}

#[derive(Deserialize)]
struct XmonQ3d {
    terminals: Vec<String>,
    maxwell_matrix_f: Vec<Vec<f64>>,
    e_j_hz: f64,
    f_r_hz: f64,
    z0_ohm: f64,
    target_ec_hz: ValidationTarget,
}

#[derive(Deserialize)]
struct SquaddsSection {
    xmon_q3d_cap_0554: XmonQ3d,
}

#[derive(Deserialize)]
struct CavityGeom {
    a_m: f64,
    b_m: f64,
    c_m: f64,
    sigma_cu_s_per_m: f64,
}

#[derive(Deserialize)]
struct DrivenTol {
    f_r_tol_rel: f64,
    q_l_tol_rel: f64,
}

#[derive(Deserialize)]
struct CavitySection {
    geometry: CavityGeom,
    targets: Vec<ValidationTarget>,
    driven_s21_clt_tolerances: DrivenTol,
}

fn solve_straight_coax_prager_synge(geom: &CoaxGeom, h2d: f64, n_layers: usize) -> (usize, f64, f64, f64) {
    let mut m = Model::new_2d();
    let outer = m.circle("diel", (0.0, 0.0), geom.b_m);
    let inner = m.circle("core", (0.0, 0.0), geom.a_m);
    m.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let mut opts = GeomMeshOpts::new(h2d);
    opts.max_dev = Some(h2d / 40.0);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (geom.a_m + geom.b_m);
    let mut inner_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        if rm < r_mid {
            if !inner_labs.contains(&lab) {
                inner_labs.push(lab);
            }
        } else if !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }

    let mesh3 =
        extrude_y(&g2.mesh, &[Layer { thickness: geom.length_m, n: n_layers, domain_offset: 0, fill: None }], 0.0)
            .unwrap();
    let n_cells = mesh3.n_cells();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(1, "eps_r", geom.eps_r);
    prob.bcs = vec![
        Bc { labels: inner_labs.clone(), kind: BcKind::Potential(1.0) },
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&mesh3, &prob).unwrap();
    let c_primal = res.scalar("capacitance").expect("coax capacitance scalar");

    // Note: straight-sided faceted circle slightly underestimates geometric annulus annulus area/volume -> compute exact
    // Prager-Synge flux bounds on the discrete tetrahedral mesh for rigorous epistemic interval checking
    let v_vals = res.solutions[0].field("V").unwrap().values.clone();
    let eps_cells = vec![geom.eps_r * EPS0; mesh3.n_cells()];
    let mut dir_labs = inner_labs;
    dir_labs.extend(outer_labs);
    let eq = equilibrate_electro_flux_bounds(
        &mesh3,
        2,
        &[v_vals],
        CellPermittivity::Scalar(&eps_cells),
        &dir_labs,
        RobinBoundary::None,
        Some(&[vec![c_primal]]),
    )
    .unwrap();
    let c_low = eq.c_lower_matrix[0][0];
    let c_up = eq.c_upper_matrix[0][0];
    assert!(c_low <= c_primal && c_primal <= c_up);
    (n_cells, c_primal, c_low, c_up)
}

fn solve_rect_annulus_prager_synge(h2d: f64, n_layers: usize, eps_r: f64) -> (f64, f64, f64) {
    let mut m = Model::new_2d();
    let outer = m.rect("diel", (0.0, 0.0), (120.0e-6, 80.0e-6));
    let inner = m.rect("core", (35.0e-6, 25.0e-6), (85.0e-6, 55.0e-6));
    m.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h2d)).unwrap();
    let mut inner_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let (xm, zm) = (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        let on_outer = xm < 1e-8 || xm > 120.0e-6 - 1e-8 || zm < 1e-8 || zm > 80.0e-6 - 1e-8;
        if on_outer {
            if !outer_labs.contains(&lab) {
                outer_labs.push(lab);
            }
        } else if !inner_labs.contains(&lab) {
            inner_labs.push(lab);
        }
    }
    let mesh3 =
        extrude_y(&g2.mesh, &[Layer { thickness: 100.0e-6, n: n_layers, domain_offset: 0, fill: None }], 0.0).unwrap();
    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(1, "eps_r", eps_r);
    prob.bcs = vec![
        Bc { labels: inner_labs.clone(), kind: BcKind::Potential(1.0) },
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&mesh3, &prob).unwrap();
    let c_primal = res.scalar("capacitance").expect("rect coax capacitance scalar");
    let v_vals = res.solutions[0].field("V").unwrap().values.clone();
    let eps_cells = vec![eps_r * EPS0; mesh3.n_cells()];
    let mut dir_labs = inner_labs;
    dir_labs.extend(outer_labs);
    let eq = equilibrate_electro_flux_bounds(
        &mesh3,
        2,
        &[v_vals],
        CellPermittivity::Scalar(&eps_cells),
        &dir_labs,
        RobinBoundary::None,
        Some(&[vec![c_primal]]),
    )
    .unwrap();
    (c_primal, eq.c_lower_matrix[0][0], eq.c_upper_matrix[0][0])
}

#[test]
fn coax_cylinder_3d_prager_synge_epistemic_and_gci_vv20() {
    let ds: CrossCodeDataset = load_validation_json("cqed_palace_q3d_squadds.json");
    let coax = &ds.coax_cylinder_electrostatic;

    let (nc_c, c_c, _lo_c, _up_c) = solve_straight_coax_prager_synge(&coax.geometry, 4.5e-6, 2);
    let (nc_f, c_f, lo_f, up_f) = solve_straight_coax_prager_synge(&coax.geometry, 2.2e-6, 4);
    let r_eff = ((nc_f as f64) / (nc_c as f64)).cbrt().max(1.2);
    let u_num = u_num_two_grid(c_c, c_f, r_eff, 2.0);

    // 1. Coax analytical validation via 2-grid u_num and Prager-Synge midpoint
    verify_asme_vv20(&coax.target, 0.5 * (lo_f + up_f), u_num);

    // 2. Pure polygonal-exact 3D rectangular coax: strict Prager-Synge bracketing [C_low, C_up]
    //    over coarse mesh vs converged fine P2 reference (zero circle-chord geometry discretization error)
    let (_, lo_rect, up_rect) = solve_rect_annulus_prager_synge(10.0e-6, 2, coax.geometry.eps_r);
    let (c_rect_ref, _, _) = solve_rect_annulus_prager_synge(4.5e-6, 4, coax.geometry.eps_r);
    let rect_tgt = inline_target(
        "rect_coax_prager_synge_bracket",
        "Polygonal-exact 3D coaxial capacitor Prager-Synge [C_low, C_up] strictly brackets fine P2 C_ref",
        "C_rect",
        "F",
        c_rect_ref,
        c_rect_ref * 1e-4,
        c_rect_ref * 1e-4,
        0.02,
        0.0,
    );
    verify_epistemic_vv20(&rect_tgt, lo_rect, up_rect, true);
}

struct GdsBuilder {
    buf: Vec<u8>,
}

impl GdsBuilder {
    fn new(user_unit: f64, dbu_m: f64) -> Self {
        let mut b = Self { buf: Vec::new() };
        b.rec(0x00, 0x02, &600i16.to_be_bytes());
        b.rec(0x01, 0x02, &[0u8; 24]);
        b.rec_ascii(0x02, "p4_pocket");
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

    fn finish(mut self) -> Vec<u8> {
        self.rec(0x04, 0x00, &[]);
        self.buf
    }
}

fn temp_dir(tag: &str) -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!(
        "nyx-p4-cross-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn build_pocket_gds(g: &PlanarPocketGeom) -> Vec<u8> {
    let sc = |v: f64| (v * 1.0e9).round() as i32;
    let mut gb = GdsBuilder::new(1.0e-3, 1.0e-9);
    gb.begin_cell("POCKET");
    let (go, px, py) = (sc(g.gnd_outer_m), sc(g.pocket_x_m), sc(g.pocket_y_m));
    // Outer ground frame formed by four boundary rectangles leaving central pocket hole
    gb.boundary(1, 0, &[[-go, -go], [go, -go], [go, -py], [-go, -py]]);
    gb.boundary(1, 0, &[[-go, py], [go, py], [go, go], [-go, go]]);
    gb.boundary(1, 0, &[[-go, -py], [-px, -py], [-px, py], [-go, py]]);
    gb.boundary(1, 0, &[[px, -py], [go, -py], [go, py], [px, py]]);

    let (phx, pyl, pyh) = (sc(g.pad_half_x_m), sc(g.pad_y_lo_m), sc(g.pad_y_hi_m));
    gb.boundary(2, 0, &[[-phx, -pyh], [phx, -pyh], [phx, -pyl], [-phx, -pyl]]);
    gb.boundary(3, 0, &[[-phx, pyl], [phx, pyl], [phx, pyh], [-phx, pyh]]);

    let (cx, cyl, cyh) = (sc(g.claw_x_m), sc(g.claw_y_lo_m), sc(g.claw_y_hi_m));
    gb.boundary(4, 0, &[[-cx, cyl], [cx, cyl], [cx, cyh], [-cx, cyh]]);
    gb.end_cell();
    gb.finish()
}

fn solve_pocket_engine(
    gds_path: &str,
    g: &PlanarPocketGeom,
    formulation: &str,
    h_mesh: f64,
    flux_bounds: Option<&str>,
) -> Value {
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "CrossCode Pocket", "dim": 3, "h": h_mesh }));
    let imp = cmd(
        &mut e,
        "import_gds",
        json!({
            "path": gds_path,
            "layers": {
                "gnd": [1, 0],
                "Q1_padA": [2, 0],
                "Q1_padB": [3, 0],
                "RO_claw": [4, 0]
            },
            "z": 0.0,
            "split_islands": false
        }),
    );
    cmd(
        &mut e,
        "add_chip_stack",
        json!({
            "bbox": imp["unpadded_bbox"],
            "lateral_pad": 320.0e-6,
            "t_sub": g.t_sub_m,
            "h_vac_top": 350.0e-6,
            "h_vac_bot": 350.0e-6
        }),
    );
    if formulation != "planar_bem" {
        for edge in ["Q1_padA:edge", "Q1_padB:edge", "RO_claw:edge"] {
            cmd(&mut e, "add_size", json!({ "sel": edge, "h": h_mesh * 0.50 }));
        }
    }
    cmd(&mut e, "set_mesh", json!({ "h": h_mesh, "gap_cells": 2.0, "max_z_slabs": 7, "z_coarsen": false }));
    let mut study_params = json!({ "formulation": formulation });
    if formulation == "planar_bem" {
        study_params["richardson"] = json!(true);
    }
    if let Some(fb) = flux_bounds {
        study_params["flux_bounds"] = json!(fb);
    }
    cmd(
        &mut e,
        "add_physics",
        json!({ "physics": "dfs", "study": "electrostatic", "order": 2, "study_params": study_params }),
    );
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": g.eps_r_si } }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } }));
    cmd(&mut e, "set_material", json!({ "physics": "dfs", "region": "VacBot", "material": { "eps_r": 1.0 } }));
    cmd(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_top.zmax", "kind": "ground" }));
    cmd(&mut e, "add_condition", json!({ "physics": "dfs", "on": "vac_bot.zmin", "kind": "ground" }));
    for term in ["gnd", "Q1_padA", "Q1_padB", "RO_claw"] {
        cmd(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": term, "kind": "terminal", "params": { "name": term } }),
        );
    }
    cmd(&mut e, "solve_run", json!({}));
    cmd(
        &mut e,
        "analyse",
        json!({
            "analysis": "capacitance_matrix",
            "physics": "dfs",
            "params": {
                "ground": "gnd",
                "qubits": [{ "name": "Q1", "pads": ["Q1_padA", "Q1_padB"], "E_J_Hz": 16.0e9 }]
            }
        }),
    )
}

#[test]
fn planar_gds_pocket_transmon_prager_synge_vs_planar_bem_and_squadds_q3d_vv20() {
    let ds: CrossCodeDataset = load_validation_json("cqed_palace_q3d_squadds.json");
    let pocket = &ds.planar_gds_pocket_cross_code;
    let dir = temp_dir("pocket");
    let gds_path = dir.join("pocket.gds");
    std::fs::write(&gds_path, build_pocket_gds(&pocket.geometry)).unwrap();
    let gds_str = gds_path.to_str().unwrap();

    // 1. Independent layered-Green-function BEM solver reference
    let bem_ana = solve_pocket_engine(gds_str, &pocket.geometry, "planar_bem", 36.0e-6, None);
    let ec_bem = bem_ana["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let c_ab_bem = -bem_ana["maxwell_matrix_F"][1][2].as_f64().unwrap();

    // 2. Coarse & finer 3D volumetric FEM with Prager-Synge flux bounds [low, up] and midpoint promotion
    let fem_c = solve_pocket_engine(gds_str, &pocket.geometry, "volumetric_fem", 50.0e-6, Some("midpoint"));
    let fem_f = solve_pocket_engine(gds_str, &pocket.geometry, "volumetric_fem", 42.0e-6, Some("midpoint"));
    let ec_fem_f = fem_f["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let c_ab_fem_f = -fem_f["maxwell_matrix_F"][1][2].as_f64().unwrap();

    let ec_low = fem_c["qubits"][0]["E_C_low_Hz"].as_f64().unwrap();
    let ec_up = fem_c["qubits"][0]["E_C_up_Hz"].as_f64().unwrap();
    let c_ab_low = (-fem_c["maxwell_upper_matrix_F"][1][2].as_f64().unwrap())
        .min(-fem_c["maxwell_lower_matrix_F"][1][2].as_f64().unwrap());
    let c_ab_up = (-fem_c["maxwell_upper_matrix_F"][1][2].as_f64().unwrap())
        .max(-fem_c["maxwell_lower_matrix_F"][1][2].as_f64().unwrap());

    // Coarse 3D Prager-Synge bounds bracket finer 3D P2 solution for identical BVP
    let ec_ps_target = inline_target(
        "pocket_ec_prager_synge_3d",
        "GDSII TransmonPocket 3D FEM Prager-Synge [E_C_low, E_C_up] vs fine 3D P2",
        "E_C",
        "Hz",
        ec_fem_f,
        ec_fem_f * 1e-3,
        ec_fem_f * 1e-3,
        pocket.tol_rel,
        0.0,
    );
    verify_epistemic_vv20(&ec_ps_target, ec_low, ec_up, true);

    let cab_ps_target = inline_target(
        "pocket_cab_prager_synge_3d",
        "GDSII TransmonPocket 3D FEM Prager-Synge [C_AB_low, C_AB_up] vs fine 3D P2",
        "C_AB",
        "F",
        c_ab_fem_f,
        c_ab_fem_f * 1e-3,
        c_ab_fem_f * 1e-3,
        pocket.tol_rel,
        0.0,
    );
    verify_epistemic_vv20(&cab_ps_target, c_ab_low, c_ab_up, true);

    // Cross-code ASME V&V 20 parity between fine 3D volumetric FEM and independent 2D planar_bem
    let ec_bem_target = inline_target(
        "pocket_ec_bem_vs_fem",
        "Differential 2-pad GDSII TransmonPocket E_C (3D FEM vs planar_bem)",
        "E_C",
        "Hz",
        ec_bem,
        ec_bem * pocket.u_input_rel,
        ec_bem * pocket.u_d_rel,
        0.04,
        0.0,
    );
    let ec_mid_c = 0.5 * (ec_low + ec_up);
    verify_asme_vv20(&ec_bem_target, ec_fem_f, u_num_two_grid(ec_mid_c, ec_fem_f, 50.0 / 42.0, 2.0));

    let cab_bem_target = inline_target(
        "pocket_c_ab_bem_vs_fem",
        "Inter-pad mutual capacitance C_AB (3D FEM vs planar_bem)",
        "C_AB",
        "F",
        c_ab_bem,
        c_ab_bem * pocket.u_input_rel,
        c_ab_bem * pocket.u_d_rel,
        0.05,
        0.0,
    );
    let c_ab_mid_c = 0.5 * (c_ab_low + c_ab_up);
    verify_asme_vv20(&cab_bem_target, c_ab_fem_f, u_num_two_grid(c_ab_mid_c, c_ab_fem_f, 50.0 / 42.0, 2.0));

    // 3. SQuADDS Q3D q3d_cap_0554 Xmon Eq. (17) Maxwell reduction verification
    let q3d = &ds.squadds_q3d_reduction.xmon_q3d_cap_0554;
    let mut mut_m = vec![vec![0.0; 3]; 3];
    for (i, mut_row) in mut_m.iter_mut().enumerate().take(3) {
        let mut s = 0.0;
        for (j, cell) in mut_row.iter_mut().enumerate().take(3) {
            let val = q3d.maxwell_matrix_f[i][j];
            s += val;
            if i != j {
                *cell = -val;
            }
        }
        mut_row[i] = s;
    }
    let res_q3d = SolveResult {
        physics: "dfs".into(),
        study: "electrostatic".into(),
        scalars: vec![CapacitanceMatrixData {
            terminals: q3d.terminals.clone(),
            maxwell_matrix_f: q3d.maxwell_matrix_f.clone(),
            mutual_matrix_f: mut_m,
            ec_matrix_hz: vec![vec![0.0; 3]; 3],
            ..Default::default()
        }
        .scalar()],
        ..Default::default()
    };
    let ana_q3d = capacitance_matrix_with_params(
        &res_q3d,
        &json!({
            "ground": "c0",
            "qubits": [{ "name": "Xmon", "pads": ["c2"], "E_J_Hz": q3d.e_j_hz }],
            "resonators": [{ "name": "RO", "terminal": "c1", "f_r_Hz": q3d.f_r_hz, "Z0_ohm": q3d.z0_ohm, "resonator_type": "quarter_wave" }]
        }),
    )
    .unwrap();
    let ec_xmon = ana_q3d["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    verify_asme_vv20(&q3d.target_ec_hz, ec_xmon, ec_xmon * 1e-4);

    let _ = std::fs::remove_dir_all(&dir);
}

fn solve_copper_cavity(g: &CavityGeom, h_mesh: f64, f_te101: f64) -> (usize, f64, f64) {
    let mu0 = 4.0e-7 * PI;
    let rs_cu = (PI * f_te101 * mu0 / g.sigma_cu_s_per_m).sqrt();

    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "3D Copper Cavity", "dim": 2, "h": h_mesh }));
    cmd(&mut e, "add_rect", json!({ "name": "base", "min": [0.0, 0.0], "max": [g.a_m, g.b_m] }));
    cmd(&mut e, "add_overlay", json!({ "regions": [{ "name": "air", "input": "base" }] }));
    let ny = ((g.c_m / h_mesh).round() as usize).max(4);
    cmd(&mut e, "set_mesh", json!({ "h": h_mesh, "sweep": [{ "thickness": g.c_m, "n": ny }] }));
    cmd(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 2,
            "study_params": { "n_modes": 1, "f_search": 6.0e9 }
        }),
    );
    cmd(
        &mut e,
        "set_material",
        json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0, "mu_r": 1.0 } }),
    );
    for side in ["outer", "ymin", "ymax"] {
        cmd(
            &mut e,
            "add_condition",
            json!({
                "physics": "ems",
                "on": side,
                "kind": "surface_impedance",
                "params": { "R_s": rs_cu, "L_k_sq": 0.0, "ibc_mode": "perturbative" }
            }),
        );
    }
    cmd(&mut e, "solve_run", json!({}));

    let n_cells = e.model().unwrap().mesh().n_cells();
    let res = e.result("ems").unwrap().0;
    let f0 = res.solutions[0].value.unwrap();
    let q0 = res.scalar("Q1").expect("lossy cavity must report perturbative wall Q1");
    (n_cells, f0, q0)
}

#[test]
fn cavity_3d_te101_q0_and_2port_driven_s21_parity_vv20() {
    let ds: CrossCodeDataset = load_validation_json("cqed_palace_q3d_squadds.json");
    let cav = &ds.cavity_3d_eigen_and_driven;
    let f_ref = cav.targets[0].reference_value;

    let (nc_c, f0_c, q0_c) = solve_copper_cavity(&cav.geometry, 0.0085, f_ref);
    let (nc_f, f0_f, q0_f) = solve_copper_cavity(&cav.geometry, 0.0055, f_ref);
    let r_eff = ((nc_f as f64) / (nc_c as f64)).cbrt().max(1.1);

    verify_asme_vv20(&cav.targets[0], f0_f, u_num_two_grid(f0_c, f0_f, r_eff, 2.0));
    verify_asme_vv20(&cav.targets[1], q0_f, u_num_two_grid(q0_c, q0_f, r_eff, 1.5));

    // Cross-code internal parity: CLT quarter-wave hanger QEP (CircuitRom) vs 2-port ems::driven Taubin S21 fit
    let mut e = Engine::new();
    let (wx, wy) = (2.4e-3, 3.8e-3);
    let (h_sub, h_top, h_bot) = (750.0e-6, 890.0e-6, 900.0e-6);
    cmd(&mut e, "project_new", json!({ "title": "CLT Driven CrossCheck", "dim": 3, "h": 420.0e-6 }));
    cmd(
        &mut e,
        "add_box",
        json!({ "name": "vac_bot", "min": [-wx/2.0, -wy/2.0, -h_sub-h_bot], "max": [wx/2.0, wy/2.0, -h_sub] }),
    );
    cmd(
        &mut e,
        "add_box",
        json!({ "name": "si_box", "min": [-wx/2.0, -wy/2.0, -h_sub], "max": [wx/2.0, wy/2.0, 0.0] }),
    );
    cmd(
        &mut e,
        "add_box",
        json!({ "name": "air_box", "min": [-wx/2.0, -wy/2.0, 0.0], "max": [wx/2.0, wy/2.0, h_top] }),
    );
    cmd(
        &mut e,
        "add_polygons",
        json!({
            "name": "gnd_poly",
            "polygons": [[[-wx/2.0, -wy/2.0], [wx/2.0, -wy/2.0], [wx/2.0, wy/2.0], [-wx/2.0, wy/2.0]]],
            "holes": [
                [[-655e-6, -1850e-6], [-525e-6, -1850e-6], [-525e-6, 1850e-6], [-655e-6, 1850e-6]],
                [[-505e-6, -1450e-6], [-385e-6, -1450e-6], [-385e-6, 1630e-6], [-505e-6, 1630e-6]]
            ]
        }),
    );
    cmd(&mut e, "add_sheet", json!({ "name": "ground_plane", "input": "gnd_poly", "z": 0.0 }));
    cmd(&mut e, "add_rect", json!({ "name": "pin_2d", "min": [-620e-6, -1850e-6], "max": [-560e-6, -1750e-6] }));
    cmd(&mut e, "add_sheet", json!({ "name": "port_in", "input": "pin_2d", "z": 0.0 }));
    cmd(&mut e, "add_rect", json!({ "name": "feed_2d", "min": [-620e-6, -1750e-6], "max": [-560e-6, 1750e-6] }));
    cmd(&mut e, "add_sheet", json!({ "name": "feed_trace", "input": "feed_2d", "z": 0.0 }));
    cmd(&mut e, "add_rect", json!({ "name": "pout_2d", "min": [-620e-6, 1750e-6], "max": [-560e-6, 1850e-6] }));
    cmd(&mut e, "add_sheet", json!({ "name": "port_out", "input": "pout_2d", "z": 0.0 }));
    cmd(&mut e, "add_rect", json!({ "name": "res_2d", "min": [-475e-6, -1450e-6], "max": [-415e-6, 1550e-6] }));
    cmd(&mut e, "add_sheet", json!({ "name": "clt_res", "input": "res_2d", "z": 0.0 }));
    cmd(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "vac_bot", "input": "vac_bot" }, { "name": "si", "input": "si_box" }, { "name": "air", "input": "air_box" }] }),
    );
    let excl = "!in_box([-650e-6, -1880e-6, -5e-6], [-530e-6, -1720e-6, 5e-6], touching) & !in_box([-650e-6, 1720e-6, -5e-6], [-530e-6, 1880e-6, 5e-6], touching)";
    cmd(&mut e, "add_size", json!({ "sel": format!("clt_res:edge & {excl}"), "h": 140.0e-6 }));
    cmd(
        &mut e,
        "set_mesh",
        json!({ "h": 450e-6, "dz_min": 140e-6, "max_z_slabs": 5, "z_growth_ratio": 2.2, "z_coarsen": true, "max_cells": 65000, "structured_pml": false }),
    );
    cmd(
        &mut e,
        "add_physics",
        json!({ "physics": "ems", "study": "cavity", "order": 1, "study_params": { "n_modes": 2, "f_search": 8.8e9 } }),
    );
    cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "vac_bot", "material": { "eps_r": 1.0 } }));
    cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
    for sh in ["ground_plane", "feed_trace", "clt_res"] {
        cmd(&mut e, "add_condition", json!({ "physics": "ems", "on": sh, "kind": "pec", "params": {} }));
    }
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "ems", "on": "port_in", "kind": "lumped_port", "params": { "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "complex_damping": false } }),
    );
    cmd(
        &mut e,
        "add_condition",
        json!({ "physics": "ems", "on": "port_out", "kind": "lumped_port", "params": { "name": "port_out", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "complex_damping": false } }),
    );
    cmd(&mut e, "solve_run", json!({}));

    let net_qep = cmd(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": { "s21_ports": ["port_in", "port_out"], "resonator_branch": "hanger" }
        }),
    );
    let res_mode = net_qep["qep_modes"]
        .as_array()
        .expect("qep_modes")
        .iter()
        .filter(|m| {
            let f = m["f_re_Hz"].as_f64().unwrap_or(0.0);
            let k = m["kappa_Hz"].as_f64().unwrap_or(0.0);
            (2.0e9..14.0e9).contains(&f) && k > 5.0e3 && k < 0.2 * f
        })
        .min_by(|a, b| a["f_re_Hz"].as_f64().unwrap().total_cmp(&b["f_re_Hz"].as_f64().unwrap()))
        .unwrap_or_else(|| panic!("no fundamental CLT QEP mode in {net_qep}"));
    let f_qep = res_mode["f_re_Hz"].as_f64().unwrap();
    let kappa_qep = res_mode["kappa_Hz"].as_f64().unwrap();
    let ql_qep = f_qep / kappa_qep.max(1.0);

    // Now run 2-port complex driven sweep around f_qep
    let f_span = (5.0 * kappa_qep).clamp(0.5e6, 150.0e6);
    cmd(
        &mut e,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "driven",
            "order": 1,
            "params": { "f_min": f_qep - f_span, "f_max": f_qep + f_span, "n_points": 51, "prom_tol": 2.0e-4, "max_hdm_solves": 6 }
        }),
    );
    cmd(&mut e, "solve_run", json!({}));
    let net_drv = cmd(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": { "s21_ports": ["port_in", "port_out"], "resonator_branch": "hanger" }
        }),
    );
    let fit = &net_drv["s21_fit"];
    let f_drv = fit["f_r_Hz"].as_f64().unwrap();
    let kappa_drv = fit["kappa_Hz"].as_f64().unwrap();
    let ql_drv = f_drv / kappa_drv.max(1.0);

    let f_tol = inline_target(
        "driven_s21_fr_vs_qep",
        "Driven S21 Taubin fit f_r vs CircuitRom QEP pole f_r",
        "f_r",
        "Hz",
        f_qep,
        f_qep * 1e-4,
        f_qep * 1e-4,
        cav.driven_s21_clt_tolerances.f_r_tol_rel,
        0.0,
    );
    verify_asme_vv20(&f_tol, f_drv, (f_drv - f_qep).abs().max(f_qep * 1e-5));

    let q_tol = inline_target(
        "driven_s21_ql_vs_qep",
        "Driven S21 Taubin fit Q_L vs CircuitRom QEP Q_L",
        "Q_L",
        "1",
        ql_qep,
        ql_qep * 1e-2,
        ql_qep * 1e-2,
        cav.driven_s21_clt_tolerances.q_l_tol_rel,
        0.0,
    );
    verify_asme_vv20(&q_tol, ql_drv, (ql_drv - ql_qep).abs().max(ql_qep * 5e-3));
}
