// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 experimental validation suite for superconducting cQED:
//! Park et al., *npj Quantum Inf.* (2026) 100 nm Nb CPW kinetic-inductance temperature sweep
//! and 35 nm disordered Nb `Q1 + R1` GDSII mask transmon + readout resonator.

mod common;

use common::{
    cmd, inline_target, load_validation_json, scalar_exact, u_num_two_grid, verify_asme_vv20, ValidationTarget,
};
use femcore::material::{bcs_surface_impedance, FilmModel, Library, Superconducting, EV_TO_J};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::gds::encode_ibm_f64;
use geom::model::Model;
use nyx_ems::Ems;
use nyx_engine::Engine;
use nyx_physics::{BcSpec, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec};
use serde::Deserialize;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Deserialize)]
struct ParkDataset {
    cpw_100nm_nb: CpwSection,
    q1_r1_35nm_disordered_nb: Q1R1Section,
}

#[derive(Deserialize)]
struct CpwGeom {
    w_trace_m: f64,
    s_gap_m: f64,
    w_box_m: f64,
    l_res_m: f64,
    h_sub_m: f64,
    h_air_m: f64,
    eps_r_si: f64,
    film_thickness_m: f64,
}

#[derive(Deserialize)]
struct BcsMatSpec {
    name: String,
    tc_k: f64,
    gap_0_j: f64,
    lambda_l_m: f64,
    xi_0_m: f64,
    mean_free_path_m: Option<f64>,
    sigma_n_s_per_m: f64,
    f_ref_hz: f64,
}

#[derive(Deserialize)]
struct CpwSection {
    geometry: CpwGeom,
    bcs_material: BcsMatSpec,
    t_over_tc_sweep: Vec<f64>,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct Q1R1Geom {
    w_box_m: f64,
    l_box_m: f64,
    h_sub_m: f64,
    h_air_m: f64,
    eps_r_sapphire: f64,
    film_thickness_m: f64,
    pad_1_box_m: [f64; 4],
    pad_2_box_m: [f64; 4],
    jj_box_m: [f64; 4],
    res_1_box_m: [f64; 4],
    e_j_hz: f64,
    c_j_f: f64,
    fock_trunc: [usize; 2],
}

#[derive(Deserialize)]
struct Q1R1Section {
    geometry: Q1R1Geom,
    bcs_material: BcsMatSpec,
    targets: Vec<ValidationTarget>,
}

struct GeomResolver<'a>(&'a MeshedGeom);

impl LabelResolver for GeomResolver<'_> {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.domains(region);
        if v.is_empty() {
            Err(Error::invalid(format!("no region '{region}'")))
        } else {
            Ok(v)
        }
    }

    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.boundary(selector);
        if !v.is_empty() {
            return Ok(v);
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

struct GdsWriter {
    buf: Vec<u8>,
}

impl GdsWriter {
    fn new(user_unit: f64, dbu_m: f64) -> Self {
        let mut b = Self { buf: Vec::new() };
        b.rec(0x00, 0x02, &600i16.to_be_bytes());
        b.rec(0x01, 0x02, &[0u8; 24]);
        b.rec_ascii(0x02, "park2026_lib");
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

    fn rect_m(&mut self, layer: i16, box_m: [f64; 4]) {
        let sc = |v: f64| (v * 1.0e9).round() as i32;
        let (x0, y0, x1, y1) = (sc(box_m[0]), sc(box_m[1]), sc(box_m[2]), sc(box_m[3]));
        let pts = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]];
        self.rec(0x08, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x0E, 0x02, &0i16.to_be_bytes());
        self.rec_xy(&pts);
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
        "nyx-p4-park2026-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn target_by_id<'a>(targets: &'a [ValidationTarget], id: &str) -> &'a ValidationTarget {
    targets.iter().find(|t| t.id == id).unwrap_or_else(|| panic!("missing target '{id}'"))
}

fn build_bcs_toml(spec: &BcsMatSpec) -> String {
    let mfp = spec.mean_free_path_m.map(|l| format!("mean_free_path = \"{l:e} m\"\n")).unwrap_or_default();
    format!(
        "format = 2\nname = \"{}\"\ndescription = \"Park 2026 Nb film\"\n\n[superconducting]\ntc = \"{} K\"\ngap_0 = \"{:e} J\"\nlambda_L = \"{:e} m\"\nxi_0 = \"{:e} m\"\n{mfp}sigma_n = \"{:e} S/m\"\n",
        spec.name, spec.tc_k, spec.gap_0_j, spec.lambda_l_m, spec.xi_0_m, spec.sigma_n_s_per_m
    )
}

fn solve_cpw_model(geom: &CpwGeom, bcs: &BcsMatSpec, h_mesh: f64) -> (usize, f64, f64, Vec<(f64, f64, f64)>) {
    let mut m = Model::new_3d();
    let si = m.box_("si_box", (-geom.w_box_m / 2.0, 0.0, -geom.h_sub_m), (geom.w_box_m / 2.0, geom.l_res_m, 0.0));
    let air = m.box_("air_box", (-geom.w_box_m / 2.0, 0.0, 0.0), (geom.w_box_m / 2.0, geom.l_res_m, geom.h_air_m));
    let trace_2d = m.rect("trace_2d", (-geom.w_trace_m / 2.0, 0.0), (geom.w_trace_m / 2.0, geom.l_res_m));
    let gnd_l_2d = m.rect("gnd_l_2d", (-geom.w_box_m / 2.0, 0.0), (-geom.w_trace_m / 2.0 - geom.s_gap_m, geom.l_res_m));
    let gnd_r_2d = m.rect("gnd_r_2d", (geom.w_trace_m / 2.0 + geom.s_gap_m, 0.0), (geom.w_box_m / 2.0, geom.l_res_m));
    m.sheet("cpw_trace", trace_2d, 0.0);
    m.sheet("cpw_gnd_l", gnd_l_2d, 0.0);
    m.sheet("cpw_gnd_r", gnd_r_2d, 0.0);
    m.overlay(&[("si", si), ("air", air)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(h_mesh)).expect("3D CPW mesh");
    let n_cells = g.mesh.n_cells();
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let mut texts = BTreeMap::new();
    texts.insert(bcs.name.clone(), build_bcs_toml(bcs));
    let lib = Arc::new(Library::default().with_texts(&texts).unwrap());

    let mut spec_pec = ProblemSpec::new("ems", 1, "cavity");
    spec_pec.study.params = json!({ "n_modes": 1, "f_search": 7.0e9 });
    spec_pec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".into(), geom.eps_r_si)].into()));
    spec_pec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    for sh in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
        spec_pec.bcs.push(BcSpec { on: sh.into(), kind: "pec".into(), params: json!({}) });
    }
    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(lib.clone());
    let f_pec = Ems.solve(&exchange, &spec_pec, &mut ctx).unwrap().solutions[0].value.unwrap();

    let solve_temp = |temp_k: f64| -> f64 {
        let mut spec = ProblemSpec::new("ems", 1, "cavity");
        spec.study.params = json!({ "n_modes": 1, "f_search": 6.9e9 });
        spec.materials = spec_pec.materials.clone();
        for sh in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
            spec.bcs.push(BcSpec {
                on: sh.into(),
                kind: "surface_impedance".into(),
                params: json!({
                    "material": bcs.name,
                    "thickness": geom.film_thickness_m,
                    "temperature": temp_k,
                    "f_ref": bcs.f_ref_hz,
                    "film_model": "kerr"
                }),
            });
        }
        let mut c = Ctx::new(&resolver);
        c.library = Some(lib.clone());
        Ems.solve(&exchange, &spec, &mut c).unwrap().solutions[0].value.unwrap()
    };

    let f_10mk = solve_temp(0.01);
    let mut sweep = Vec::new();
    let sc = Superconducting {
        source: Vec::new(),
        tc: bcs.tc_k,
        gap_0: Some(bcs.gap_0_j),
        lambda_l: Some(bcs.lambda_l_m),
        xi_0: Some(bcs.xi_0_m),
        mean_free_path: bcs.mean_free_path_m,
        sigma_n: Some(bcs.sigma_n_s_per_m),
        ..Default::default()
    };
    let eval_0 = bcs_surface_impedance(&sc, bcs.f_ref_hz, 0.01, geom.film_thickness_m, FilmModel::Kerr).unwrap();
    let gamma_geom = ((f_pec / f_10mk).powi(2) - 1.0) / eval_0.l_eff;
    for ratio in [0.01, 0.25, 0.50, 0.65, 0.75, 0.85] {
        let tk = ratio * bcs.tc_k;
        let eval_t = bcs_surface_impedance(&sc, bcs.f_ref_hz, tk, geom.film_thickness_m, FilmModel::Kerr).unwrap();
        let f_analytic = f_pec / (1.0 + gamma_geom * eval_t.l_eff).sqrt();
        let f_fem = solve_temp(tk);
        sweep.push((ratio, f_fem, f_analytic));
    }

    (n_cells, f_pec, f_10mk, sweep)
}

#[test]
fn park2026_100nm_nb_cpw_kinetic_inductance_and_bcs_sweep_vv20() {
    let ds: ParkDataset = load_validation_json("cqed_park2026.json");
    let cpw = &ds.cpw_100nm_nb;

    let (nc_c, f_pec_c, f_0_c, _) = solve_cpw_model(&cpw.geometry, &cpw.bcs_material, 1.6e-3);
    let (nc_f, f_pec_f, f_0_f, sweep) = solve_cpw_model(&cpw.geometry, &cpw.bcs_material, 1.1e-3);
    let r_eff = ((nc_f as f64) / (nc_c as f64)).cbrt().max(1.1);

    let u_pec = u_num_two_grid(f_pec_c, f_pec_f, r_eff, 1.0);
    let u_0 = u_num_two_grid(f_0_c, f_0_f, r_eff, 1.0);
    let shift_c = f_0_c - f_pec_c;
    let shift_f = f_0_f - f_pec_f;
    let u_shift = u_num_two_grid(shift_c, shift_f, r_eff, 1.0);

    verify_asme_vv20(target_by_id(&cpw.targets, "park2026_cpw_f_pec"), f_pec_f, u_pec);
    verify_asme_vv20(target_by_id(&cpw.targets, "park2026_cpw_f_10mk"), f_0_f, u_0);
    verify_asme_vv20(target_by_id(&cpw.targets, "park2026_cpw_df_kin"), shift_f, u_shift);

    assert_eq!(sweep.len(), cpw.t_over_tc_sweep.len());
    for (idx, &(ratio, f_fem, f_analytic)) in sweep.iter().enumerate() {
        let u_t = u_0.max(f_fem * 2e-4);
        let t_target = inline_target(
            format!("park2026_bcs_t_{idx}"),
            format!("BCS T/Tc = {ratio:.2}"),
            "f_BCS(T)",
            "Hz",
            f_analytic,
            2.5e6,
            1.0e6,
            0.005,
            0.0,
        );
        verify_asme_vv20(&t_target, f_fem, u_t);
    }
}

fn solve_q1_r1_gds_via_engine(gds_path: &str, sec: &Q1R1Section, h_mesh: f64) -> (usize, f64, f64, f64, f64) {
    let g = &sec.geometry;
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "Park2026 Q1+R1 GDS", "dim": 3, "h": h_mesh }));
    cmd(
        &mut e,
        "import_gds",
        json!({
            "path": gds_path,
            "layers": {
                "pad_1": [1, 0],
                "pad_2": [2, 0],
                "q1_jj": [3, 0],
                "res_1": [4, 0]
            },
            "z": 0.0,
            "split_islands": false,
            "layer_bbox": false
        }),
    );
    cmd(
        &mut e,
        "add_chip_stack",
        json!({
            "bbox": [-g.w_box_m / 2.0, -g.l_box_m / 2.0, g.w_box_m / 2.0, g.l_box_m / 2.0],
            "t_sub": g.h_sub_m,
            "h_vac_top": g.h_air_m,
            "h_vac_bot": 0.0
        }),
    );
    for edge in ["pad_1:edge", "pad_2:edge", "q1_jj:edge"] {
        cmd(&mut e, "add_size", json!({ "sel": edge, "h": h_mesh * 0.25 }));
    }
    cmd(&mut e, "add_size", json!({ "sel": "res_1:edge", "h": h_mesh * 0.65 }));
    cmd(
        &mut e,
        "set_mesh",
        json!({
            "h": h_mesh,
            "max_z_slabs": 5,
            "gap_cells": 2.0,
            "structured_pml": false
        }),
    );
    cmd(&mut e, "add_material", json!({ "name": sec.bcs_material.name, "toml": build_bcs_toml(&sec.bcs_material) }));
    cmd(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 2,
            "study_params": { "n_modes": 2, "f_search": 5.8e9 }
        }),
    );
    cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "Si", "material": { "eps_r": g.eps_r_sapphire } }));
    cmd(&mut e, "set_material", json!({ "physics": "ems", "region": "VacTop", "material": { "eps_r": 1.0 } }));
    for sheet in ["pad_1", "pad_2", "res_1"] {
        cmd(
            &mut e,
            "add_condition",
            json!({
                "physics": "ems",
                "on": sheet,
                "kind": "surface_impedance",
                "params": {
                    "material": sec.bcs_material.name,
                    "thickness": g.film_thickness_m,
                    "temperature": 0.01,
                    "f_ref": sec.bcs_material.f_ref_hz,
                    "film_model": "kerr"
                }
            }),
        );
    }
    cmd(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "q1_jj",
            "kind": "lumped_jj",
            "params": {
                "name": "Q1_JJ",
                "E_j_Hz": g.e_j_hz,
                "C_j": g.c_j_f,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }
        }),
    );

    cmd(&mut e, "solve_run", json!({}));
    let n_cells = e.model().unwrap().mesh().n_cells();
    let res = e.result("ems").unwrap().0;
    let rom = CircuitRom::of(res).expect("Park2026 GDS Q1+R1 cavity must emit CircuitRom");
    let _c_port = rom.c_port_f[0][0];
    let epr = cmd(&mut e, "analyse", json!({ "analysis": "epr", "physics": "ems" }));
    assert!(epr["equipartition_error"].as_f64().unwrap() < 1e-10);
    assert!(
        epr["kinetic_inductance_fraction"].as_array().unwrap().iter().filter_map(Value::as_f64).any(|ak| ak > 1e-5),
        "35 nm disordered Nb sheets must give non-zero kinetic inductance fraction"
    );

    let ham = cmd(
        &mut e,
        "analyse",
        json!({
            "analysis": "hamiltonian",
            "physics": "ems",
            "params": { "basis": "normal_fock", "fock_trunc": g.fock_trunc }
        }),
    );
    let _ = EV_TO_J;
    let _ = scalar_exact;

    let f01 = ham["f01_Hz"].as_array().unwrap();
    let f_q1 = f01[0].as_f64().unwrap();
    let f_r1 = f01[1].as_f64().unwrap();
    let alpha_1 = ham["alpha_Hz"][0].as_f64().unwrap();
    let chi_q1r1 = ham["chi_matrix_Hz"][0][1].as_f64().unwrap();
    (n_cells, f_q1, alpha_1, f_r1, chi_q1r1)
}

#[test]
fn park2026_35nm_disordered_nb_gds_mask_q1_r1_vv20() {
    let ds: ParkDataset = load_validation_json("cqed_park2026.json");
    let sec = &ds.q1_r1_35nm_disordered_nb;
    let g = &sec.geometry;

    let dir = temp_dir("q1r1");
    let gds_path = dir.join("park2026_q1_r1.gds");
    let mut gw = GdsWriter::new(1.0e-3, 1.0e-9);
    gw.begin_cell("PARK_Q1_R1");
    gw.rect_m(1, g.pad_1_box_m);
    gw.rect_m(2, g.pad_2_box_m);
    gw.rect_m(3, g.jj_box_m);
    gw.rect_m(4, g.res_1_box_m);
    gw.end_cell();
    std::fs::write(&gds_path, gw.finish()).unwrap();
    let gds_str = gds_path.to_str().unwrap();

    let (nc_c, fq_c, alpha_c, fr_c, chi_c) = solve_q1_r1_gds_via_engine(gds_str, sec, 780.0e-6);
    let (nc_f, fq_f, alpha_f, fr_f, chi_f) = solve_q1_r1_gds_via_engine(gds_str, sec, 560.0e-6);
    let r_eff = ((nc_f as f64) / (nc_c as f64)).cbrt().max(1.05);

    verify_asme_vv20(target_by_id(&sec.targets, "park2026_q1_f01"), fq_f, u_num_two_grid(fq_c, fq_f, r_eff, 1.5));
    verify_asme_vv20(
        target_by_id(&sec.targets, "park2026_q1_alpha"),
        alpha_f,
        u_num_two_grid(alpha_c, alpha_f, r_eff, 1.5),
    );
    verify_asme_vv20(target_by_id(&sec.targets, "park2026_r1_freq"), fr_f, u_num_two_grid(fr_c, fr_f, r_eff, 1.5));
    verify_asme_vv20(target_by_id(&sec.targets, "park2026_chi_q1r1"), chi_f, u_num_two_grid(chi_c, chi_f, r_eff, 1.5));

    let _ = std::fs::remove_dir_all(&dir);
}
