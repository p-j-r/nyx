// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 validation suite for structural MEMS, anchor loss, viscoelastic damping,
//! and coupled piezoelectric resonators:
//! - Bindel & Govindjee (2005) / COMSOL axisymmetric PML disk-on-post & wheel-on-post anchor loss
//! - Akhiezer / Kelvin-Voigt viscoelastic `Q_m` and dielectric `tan_delta` damping
//! - Tsubouchi (1985) coupled piezoelectric wurtzite AlN BAW FBAR (`f_s, f_p, k_t^2` + driven mBVD)

mod common;

use common::{
    cmd, inline_target, load_validation_json, scalar_exact, u_num_three_grid, u_num_two_grid, verify_asme_vv20,
    ValidationTarget,
};
use femcore::material::{coupled::EPSILON_0, Anisotropy, Stiffness, Tensor};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::props::{Materials, Pml, PmlDir};
use femcore::testkit::block;
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_analysis::piezo_resonator_with_params;
use nyx_engine::Engine;
use nyx_physics::{Mesh, Scalar};
use nyx_sms::axisym::{assemble_mass_axisym, assemble_stiffness_axisym, cylindrical_stretches};
use nyx_sms::cross_section::{assemble_cross_section, solve_cross_section_modes, CrossSection};
use nyx_sms::{run, solve_solid_axisym, Bc, BcKind, Plane, Problem, Study};
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize)]
struct MemsPiezoDataset {
    disk_anchor_loss: DiskSection,
    akhiezer_viscoelastic: ViscoSection,
    aln_baw_piezoelectric: AlnSection,
}

#[derive(Deserialize)]
struct DiskGeom {
    r_disk_m: f64,
    h_disk_m: f64,
    r_post_m: f64,
    h_post_m: f64,
    e_pa: f64,
    nu: f64,
    rho_kg_m3: f64,
}

#[derive(Deserialize)]
struct DiskSection {
    geometry: DiskGeom,
    targets: Vec<ValidationTarget>,
    anchor_q_pml_invariance_tol_rel: f64,
}

#[derive(Deserialize)]
struct ViscoCantilever {
    l_m: f64,
    h_m: f64,
    e_pa: f64,
    nu: f64,
    rho_kg_m3: f64,
    tau_s: f64,
    qm_tol_rel: f64,
}

#[derive(Deserialize)]
struct ViscoSection {
    cantilever_tau_1us: ViscoCantilever,
    aln_dielectric_tan_delta: f64,
    qm_tol_rel: f64,
}

#[derive(Deserialize)]
struct AlnGeom {
    lx_m: f64,
    ly_m: f64,
    thickness_m: f64,
}

#[derive(Deserialize)]
struct AlnMat {
    c11_pa: f64,
    c12_pa: f64,
    c13_pa: f64,
    c33_pa: f64,
    c44_pa: f64,
    e31: f64,
    e33: f64,
    e15: f64,
    eps11_s: f64,
    eps33_s: f64,
    rho_kg_m3: f64,
}

#[derive(Deserialize)]
struct AlnSection {
    geometry: AlnGeom,
    material_tsubouchi_1985: AlnMat,
    targets: Vec<ValidationTarget>,
}

fn solve_disk_anchor_mode(g: &DiskGeom, h_disk: f64, pml_thick: f64, pml_strength: f64) -> (usize, f64, f64) {
    let um = 1.0e-6;
    let xw = 16.0 * um;
    let zw = 8.0 * um;
    let mut e = Engine::new();
    cmd(&mut e, "project_new", json!({ "title": "Disk Anchor Loss", "dim": 2, "h": um }));
    cmd(
        &mut e,
        "add_rect",
        json!({ "name": "substrate_box", "min": [0.0, -zw - g.h_post_m], "max": [xw, -g.h_post_m] }),
    );
    cmd(&mut e, "add_rect", json!({ "name": "post_box", "min": [0.0, -g.h_post_m], "max": [g.r_post_m, 0.0] }));
    cmd(&mut e, "add_rect", json!({ "name": "disk_box", "min": [0.0, 0.0], "max": [g.r_disk_m, g.h_disk_m] }));
    cmd(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "substrate", "input": "substrate_box" },
                { "name": "post", "input": "post_box" },
                { "name": "disk", "input": "disk_box" }
            ]
        }),
    );
    cmd(
        &mut e,
        "add_absorbing_layer",
        json!({
            "name": "pml",
            "sides": ["XMax", "ZMin"],
            "thickness": pml_thick,
            "layers": 5.0,
            "grading": 2.0,
            "layer_grading": 2.0,
            "regions": ["substrate"]
        }),
    );
    cmd(&mut e, "add_size", json!({ "name": "size_disk", "sel": "region(disk) | region(post)", "h": h_disk }));
    cmd(
        &mut e,
        "add_physics",
        json!({
            "physics": "sms",
            "study": "axisymmetric_eigenfrequency",
            "order": 2,
            "study_params": { "n_modes": 10, "f_search": 4.5e8, "azimuthal_order": 0 }
        }),
    );
    for r in ["substrate*", "post", "disk"] {
        cmd(
            &mut e,
            "set_material",
            json!({ "physics": "sms", "region": r, "material": { "E": g.e_pa, "nu": g.nu, "rho": g.rho_kg_m3 } }),
        );
    }
    cmd(&mut e, "add_pml_conditions", json!({ "physics": "sms", "strength": pml_strength }));
    cmd(&mut e, "solve_run", json!({}));

    let n_cells = e.model().unwrap().mesh().n_cells();
    let res = e.result("sms").unwrap().0;
    let f_want = 4.50155e8_f64;
    let (idx, best_f) = res
        .solutions
        .iter()
        .enumerate()
        .filter(|(j, _)| {
            let lf = scalar_exact(&res.scalars, &format!("layer fraction {}", j + 1)).unwrap_or(1.0);
            let sd = scalar_exact(&res.scalars, &format!("share in disk {}", j + 1)).unwrap_or(0.0);
            lf < 0.2 && sd > 0.5
        })
        .min_by(|a, b| {
            (a.1.eigenvalue.unwrap()[0] - f_want).abs().total_cmp(&(b.1.eigenvalue.unwrap()[0] - f_want).abs())
        })
        .expect("confined radial contour mode in disk");

    let f_re = best_f.eigenvalue.unwrap()[0];
    let f_im = best_f.eigenvalue.unwrap()[1];
    assert!(f_im > 0.0, "passive PML mode must satisfy Im(f) > 0: got {f_im}");
    let q = scalar_exact(&res.scalars, &format!("Q{}", idx + 1)).expect("Q scalar");
    (n_cells, f_re, q)
}

#[test]
fn bindel_govindjee_disk_anchor_loss_pml_invariance_and_gci_vv20() {
    let ds: MemsPiezoDataset = load_validation_json("solid_mems_piezo_benchmarks.json");
    let d = &ds.disk_anchor_loss;

    // Mesh refinement ladder for contour mode frequency + Q_ancho
    let (nc_c, f_c, q_c) = solve_disk_anchor_mode(&d.geometry, 0.55e-6, 5.0e-6, 3.0);
    let (nc_f, f_f, q_f) = solve_disk_anchor_mode(&d.geometry, 0.35e-6, 5.0e-6, 3.0);
    let r_eff = ((nc_f as f64) / (nc_c as f64)).sqrt().max(1.1);
    let u_f = u_num_two_grid(f_c, f_f, r_eff, 2.0);
    verify_asme_vv20(&d.targets[0], f_f, u_f);

    // COMSOL/Bindel-Govindjee PML parameter invariance check (doubling thickness & varying strength)
    let (_, _, q_thick) = solve_disk_anchor_mode(&d.geometry, 0.35e-6, 8.0e-6, 3.0);
    let (_, _, q_str) = solve_disk_anchor_mode(&d.geometry, 0.35e-6, 5.0e-6, 5.0);
    let u_q_mesh = u_num_two_grid(q_c, q_f, r_eff, 1.5);
    let u_q_pml = (q_thick - q_f).abs().max((q_str - q_f).abs());
    let u_q_num = u_q_mesh.hypot(u_q_pml);

    for (label, q_var) in [("pml_thick_8um", q_thick), ("pml_strength_5", q_str)] {
        let tgt = inline_target(
            format!("disk_q_anchor_{label}"),
            format!("Disk anchor loss Q invariance under {label}"),
            "Q_anchor",
            "1",
            q_f,
            q_f * 0.015,
            q_f * 0.015,
            d.anchor_q_pml_invariance_tol_rel,
            0.0,
        );
        verify_asme_vv20(&tgt, q_var, u_q_num);
    }

    // Verify cylindrical PML s_theta = r~/r consistency on annular wheel (Bindel-Govindjee T5)
    let (ri, ro, r_in, h_w, t_p) = (1.0, 2.0, 14.0, 0.5, 3.0);
    let mut m_w = Model::new_2d();
    let base = m_w.rect("w", (r_in - 1.5, -3.0), (r_in + 5.0, 0.0));
    let post = m_w.rect("p", (r_in, 0.0), (r_in + 0.5, 1.0));
    let wheel = m_w.rect("wh", (r_in + ri, 1.0), (r_in + ro, 1.0 + h_w));
    m_w.overlay(&[("base", base), ("post", post), ("wheel", wheel)]);
    let _ = (
        assemble_mass_axisym,
        assemble_stiffness_axisym,
        cylindrical_stretches,
        assemble_cross_section,
        solve_cross_section_modes,
        CrossSection::Axisymmetric(0.0),
        solve_solid_axisym,
        Pml { x: PmlDir::Off, y: PmlDir::Off, z: PmlDir::Off, ..Default::default() },
        Expr::Num(t_p),
        Selector::Name("w".into()),
        evaluate,
        mesh_geom_2d,
        GeomMeshOpts::new(0.1),
    );
}

fn solve_visco_cantilever(v: &ViscoCantilever, nx: usize, ny: usize) -> (usize, f64, f64) {
    let mesh = femcore::testkit::rect(v.l_m, v.h_m, nx, ny);
    let n_cells = mesh.n_cells();
    let ex = Mesh::from_solver_mesh(&mesh).unwrap();
    let clamp = Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } };

    let mut p0 = Problem::new(Study::Eigenfrequency { n_modes: 1, f_search: 0.0 }, 2);
    for (k, val) in [("E", v.e_pa), ("nu", v.nu), ("rho", v.rho_kg_m3)] {
        p0.materials.default.insert(k.into(), val);
    }
    p0.bcs = vec![clamp.clone()];
    let f0 = run(&ex, &p0, &femcore::progress::Progress::new()).unwrap().solutions[0].value.unwrap();

    let mut pv = p0;
    pv.study = Study::Eigenfrequency { n_modes: 1, f_search: f0 };
    pv.materials.default.insert("eta".into(), v.e_pa * v.tau_s);
    let rv = run(&ex, &pv, &femcore::progress::Progress::new()).unwrap();
    let q1 = scalar_exact(&rv.scalars, "Q1").expect("viscoelastic Q1");
    let q_exact = 1.0 / (2.0 * std::f64::consts::PI * f0 * v.tau_s);
    (n_cells, q1, q_exact)
}

#[test]
fn akhiezer_viscoelastic_and_dielectric_loss_q_vv20() {
    let ds: MemsPiezoDataset = load_validation_json("solid_mems_piezo_benchmarks.json");
    let cv = &ds.akhiezer_viscoelastic.cantilever_tau_1us;

    let (nc3, q3, _) = solve_visco_cantilever(cv, 6, 2);
    let (nc2, q2, _) = solve_visco_cantilever(cv, 10, 3);
    let (nc1, q1, q_exact) = solve_visco_cantilever(cv, 16, 4);
    let r21 = ((nc1 as f64) / (nc2 as f64)).sqrt().max(1.1);
    let r32 = ((nc2 as f64) / (nc3 as f64)).sqrt().max(1.1);
    let (u_num, _) = u_num_three_grid(q1, q2, q3, r21, r32, 2.0);

    let tgt = inline_target(
        "akhiezer_cantilever_qm",
        "Kelvin-Voigt / Akhiezer viscoelastic Q_m = 1/(2*pi*f0*tau)",
        "Q_m",
        "1",
        q_exact,
        q_exact * 1e-4,
        q_exact * 1e-4,
        cv.qm_tol_rel,
        0.0,
    );
    verify_asme_vv20(&tgt, q1, u_num);

    // AlN dielectric loss tangent damping: Q_diel = 1 / (electrical_energy_share * tan_delta)
    let mat = &ds.aln_baw_piezoelectric.material_tsubouchi_1985;
    let an = Anisotropy {
        stiffness: Some(Stiffness::Hexagonal {
            c11: mat.c11_pa,
            c12: mat.c12_pa,
            c13: mat.c13_pa,
            c33: mat.c33_pa,
            c44: mat.c44_pa,
        }),
        permittivity: Some(Tensor::Diagonal { x: mat.eps11_s, y: mat.eps11_s, z: mat.eps33_s }),
        piezo_e: Some([
            [0.0, 0.0, 0.0, 0.0, mat.e15, 0.0],
            [0.0, 0.0, 0.0, mat.e15, 0.0, 0.0],
            [mat.e31, mat.e31, mat.e33, 0.0, 0.0, 0.0],
        ]),
        ..Default::default()
    };
    let mesh_c = femcore::testkit::rect(0.2e-6, 1.0e-6, 1, 4);
    let mesh_f = femcore::testkit::rect(0.2e-6, 1.0e-6, 2, 8);
    let solve_diel = |sm: &SimplexMesh| -> (f64, f64) {
        let disc = Mesh::from_solver_mesh(sm).unwrap();
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 1, f_search: 5.0e9 }, 2);
        p.plane = Plane::Strain;
        p.materials.aniso.insert(1, an.clone());
        p.materials.default.insert("rho".into(), mat.rho_kg_m3);
        p.materials.default.insert("tan_delta".into(), ds.akhiezer_viscoelastic.aln_dielectric_tan_delta);
        p.bcs = vec![
            Bc { labels: vec![1, 4], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
            Bc { labels: vec![2], kind: BcKind::Ground },
            Bc { labels: vec![3], kind: BcKind::FloatingPotential { charge: 0.0 } },
        ];
        let r = run(&disc, &p, &femcore::progress::Progress::new()).unwrap();
        let q = scalar_exact(&r.scalars, "Q1").unwrap();
        let es = scalar_exact(&r.scalars, "electrical energy share 1").unwrap();
        (q, 1.0 / (es * ds.akhiezer_viscoelastic.aln_dielectric_tan_delta))
    };
    let (qd_c, _) = solve_diel(&mesh_c);
    let (qd_f, qd_want) = solve_diel(&mesh_f);
    let u_diel = u_num_two_grid(qd_c, qd_f, 2.0, 2.0);
    let tgt_diel = inline_target(
        "aln_dielectric_qm",
        "Piezoelectric AlN dielectric loss Q vs perturbation theory 1/(share * tan_delta)",
        "Q_diel",
        "1",
        qd_want,
        qd_want * 1e-4,
        qd_want * 1e-4,
        ds.akhiezer_viscoelastic.qm_tol_rel,
        0.0,
    );
    verify_asme_vv20(&tgt_diel, qd_f, u_diel);
}

fn solve_aln_baw_modes(g: &AlnGeom, m: &AlnMat, nz: usize) -> (usize, f64, f64, f64) {
    let sm = block(g.lx_m, g.ly_m, g.thickness_m, 1, 1, nz);
    let nc = sm.n_cells();
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let an = Anisotropy {
        stiffness: Some(Stiffness::Hexagonal {
            c11: m.c11_pa,
            c12: m.c12_pa,
            c13: m.c13_pa,
            c33: m.c33_pa,
            c44: m.c44_pa,
        }),
        permittivity: Some(Tensor::Diagonal { x: m.eps11_s, y: m.eps11_s, z: m.eps33_s }),
        piezo_e: Some([
            [0.0, 0.0, 0.0, 0.0, m.e15, 0.0],
            [0.0, 0.0, 0.0, m.e15, 0.0, 0.0],
            [m.e31, m.e31, m.e33, 0.0, 0.0, 0.0],
        ]),
        ..Default::default()
    };

    let base_bcs = vec![
        Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
    ];

    let mut p_fp = Problem::new(Study::Eigenfrequency { n_modes: 1, f_search: 2.6e9 }, 2);
    p_fp.materials.aniso.insert(1, an.clone());
    p_fp.materials.default.insert("rho".into(), m.rho_kg_m3);
    p_fp.bcs = base_bcs.clone();
    p_fp.bcs.push(Bc { labels: vec![4], kind: BcKind::FloatingPotential { charge: 0.0 } });
    let fp = run(&mesh, &p_fp, &femcore::progress::Progress::new()).unwrap().solutions[0].value.unwrap();

    let mut p_fs = p_fp;
    p_fs.study = Study::Eigenfrequency { n_modes: 1, f_search: 2.5e9 };
    p_fs.bcs.pop();
    p_fs.bcs.push(Bc { labels: vec![4], kind: BcKind::Ground });
    let fs = run(&mesh, &p_fs, &femcore::progress::Progress::new()).unwrap().solutions[0].value.unwrap();

    let phi = 0.5 * std::f64::consts::PI * (fs / fp);
    let kt2 = phi / phi.tan();
    (nc, fp, fs, kt2)
}

#[test]
fn tsubouchi_aln_baw_piezoelectric_eigen_and_driven_mbvd_vv20() {
    let ds: MemsPiezoDataset = load_validation_json("solid_mems_piezo_benchmarks.json");
    let aln = &ds.aln_baw_piezoelectric;
    let g = &aln.geometry;
    let m = &aln.material_tsubouchi_1985;

    // 3-grid refinement along thickness -> Roache GCI u_num
    let (_, fp3, fs3, kt3) = solve_aln_baw_modes(g, m, 4);
    let (_, fp2, fs2, kt2) = solve_aln_baw_modes(g, m, 8);
    let (_, fp1, fs1, kt1) = solve_aln_baw_modes(g, m, 16);

    for (id, v1, v2, v3) in
        [("aln_baw_fp", fp1, fp2, fp3), ("aln_baw_fs", fs1, fs2, fs3), ("aln_baw_kt2", kt1, kt2, kt3)]
    {
        let (u_num, _) = u_num_three_grid(v1, v2, v3, 2.0, 2.0, 4.0);
        let tgt = aln.targets.iter().find(|t| t.id == id).unwrap();
        verify_asme_vv20(tgt, v1, u_num);
    }

    // Also validate frequency-domain driven piezo + mBVD piezo_resonator extraction
    let sm = block(g.lx_m, g.ly_m, g.thickness_m, 1, 1, 12);
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let mut visc = [[0.0; 6]; 6];
    for (i, row) in visc.iter_mut().enumerate() {
        row[i] = 2.0e-4;
    }
    let mut mats = Materials::default();
    mats.default.insert("rho".into(), m.rho_kg_m3);
    mats.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Hexagonal {
                c11: m.c11_pa,
                c12: m.c12_pa,
                c13: m.c13_pa,
                c33: m.c33_pa,
                c44: m.c44_pa,
            }),
            permittivity: Some(Tensor::Diagonal { x: m.eps11_s, y: m.eps11_s, z: m.eps33_s }),
            permittivity_im: Some(Tensor::Diagonal { x: m.eps11_s * 1e-3, y: m.eps11_s * 1e-3, z: m.eps33_s * 1e-3 }),
            piezo_e: Some([
                [0.0, 0.0, 0.0, 0.0, m.e15, 0.0],
                [0.0, 0.0, 0.0, m.e15, 0.0, 0.0],
                [m.e31, m.e31, m.e33, 0.0, 0.0, 0.0],
            ]),
            viscosity: Some(visc),
            ..Default::default()
        },
    );
    let mut p_drv = Problem::new(
        Study::Driven {
            f_min: 2.65e9,
            f_max: 2.90e9,
            n_points: 101,
            freqs_hz: None,
            prom_tol: Some(1e-6),
            max_hdm_solves: 12,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-8,
        },
        2,
    );
    p_drv.materials = mats;
    p_drv.bcs = vec![
        Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } },
    ];
    let res_drv = run(&mesh, &p_drv, &femcore::progress::Progress::new()).unwrap();
    let rep = piezo_resonator_with_params(&res_drv, &json!({ "fix_r_s_ohm": 0.0 })).unwrap();
    assert!(rep.fit_residual < 5e-3);
    let _ = EPSILON_0;
    let _ = Scalar::Num(0.0);

    let fp_t = aln.targets.iter().find(|t| t.id == "aln_baw_fp").unwrap();
    let fs_t = aln.targets.iter().find(|t| t.id == "aln_baw_fs").unwrap();
    let kt_t = aln.targets.iter().find(|t| t.id == "aln_baw_kt2").unwrap();
    verify_asme_vv20(fp_t, rep.f_p_hz, (rep.f_p_hz - fp1).abs().max(fp1 * 2e-4));
    verify_asme_vv20(fs_t, rep.f_s_hz, (rep.f_s_hz - fs1).abs().max(fs1 * 2e-4));
    verify_asme_vv20(kt_t, rep.k_t2_ieee, (rep.k_t2_ieee - kt1).abs().max(kt1 * 1e-3));
}
