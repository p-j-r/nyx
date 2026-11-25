// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for Milestone M2: coupled piezoelectric `Study::Driven` multi-port
//! admittance (`Y`), impedance (`Z`), scattering (`S`), clamped capacitance (`C_0`),
//! Hermite PROM acceleration, `auto`/`direct`/`lblt`/`blr` solver backends, mixed-cell PML,
//! `conductor` `role = "terminal"`, and negative/conflict validation checks.

use femcore::linalg::csr::{c, C};
use femcore::material::{coupled::EPSILON_0, Anisotropy, Stiffness, Tensor};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::props::{Materials, Pml, PmlDir};
use femcore::testkit::block;
use nyx_physics::{validate_spec, BcSpec, Ctx, DrivenPortData, Mesh, Physics, ProblemSpec};
use nyx_sms::problem::{decode, Bc, BcKind, ConductorRole, Problem, Study};
use nyx_sms::Sms;
use serde_json::json;
use std::collections::BTreeMap;
use std::f64::consts::PI;

#[derive(Clone)]
struct MapResolver {
    domains: BTreeMap<String, Vec<i32>>,
    boundaries: BTreeMap<String, Vec<i32>>,
}

impl nyx_physics::LabelResolver for MapResolver {
    fn domains(&self, region: &str) -> Result<Vec<i32>, nyx_physics::Error> {
        self.domains.get(region).cloned().ok_or_else(|| nyx_physics::Error::invalid(format!("no region '{region}'")))
    }

    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, nyx_physics::Error> {
        self.boundaries
            .get(selector)
            .cloned()
            .ok_or_else(|| nyx_physics::Error::invalid(format!("no boundary '{selector}'")))
    }
}

fn aln_stiffness_and_piezo() -> (Stiffness, [[f64; 6]; 3], Tensor) {
    let stiff = Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: 389e9, c44: 125e9 };
    let e_tb =
        [[0.0, 0.0, 0.0, 0.0, -0.48, 0.0], [0.0, 0.0, 0.0, -0.48, 0.0, 0.0], [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0]];
    let eps = Tensor::Diagonal { x: 9.21, y: 9.21, z: 10.12 };
    (stiff, e_tb, eps)
}

fn build_aln_materials(eta_iso: f64, tan_delta: f64, euler_deg: [f64; 3]) -> Materials {
    let (stiff, e_tb, eps) = aln_stiffness_and_piezo();
    let mut mats = Materials::default();
    mats.default.insert("rho".into(), 3260.0);
    mats.default.insert("eps_r".into(), 10.12);
    if tan_delta > 0.0 {
        mats.default.insert("eps_r_im".into(), 10.12 * tan_delta);
    }
    let eps_im = if tan_delta > 0.0 {
        Some(Tensor::Diagonal { x: 9.21 * tan_delta, y: 9.21 * tan_delta, z: 10.12 * tan_delta })
    } else {
        None
    };
    let visc = if eta_iso > 0.0 {
        let mut v = [[0.0; 6]; 6];
        for (i, row) in v.iter_mut().enumerate() {
            row[i] = eta_iso;
        }
        Some(v)
    } else {
        None
    };
    mats.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(stiff),
            permittivity: Some(eps),
            permittivity_im: eps_im,
            piezo_e: Some(e_tb),
            viscosity: visc,
            euler_deg,
            ..Default::default()
        },
    );
    mats
}

/// Builds a 2D `Tri` + `Quad` mixed-cell rectangle `[0, lx] x [0, dz]` with left `tri_nx` columns
/// triangulated (domain `d_left`) and right `quad_nx` columns as Quads (domain `d_right`).
/// Boundary labels: 1 = left (x=0), 2 = right (x=lx), 3 = bottom (z=0), 4 = top (z=dz).
fn mixed_tri_quad_rect(lx: f64, dz: f64, tri_nx: usize, quad_nx: usize, nz: usize, d_left: i32, d_right: i32) -> Mesh {
    let nx = tri_nx + quad_nx;
    let dx = lx / nx as f64;
    let h_z = dz / nz as f64;
    let mut cm = CellMesh::<2>::new();
    let idx = |ix: usize, iz: usize| -> u32 { (iz * (nx + 1) + ix) as u32 };
    for iz in 0..=nz {
        for ix in 0..=nx {
            cm.points.push([ix as f64 * dx, iz as f64 * h_z]);
        }
    }
    for iz in 0..nz {
        for ix in 0..nx {
            let (p00, p10, p11, p01) = (idx(ix, iz), idx(ix + 1, iz), idx(ix + 1, iz + 1), idx(ix, iz + 1));
            if ix < tri_nx {
                cm.types.push(CellType::Tri);
                cm.conn.extend([p00, p10, p11]);
                cm.offsets.push(cm.conn.len() as u32);
                cm.label.push(d_left);

                cm.types.push(CellType::Tri);
                cm.conn.extend([p00, p11, p01]);
                cm.offsets.push(cm.conn.len() as u32);
                cm.label.push(d_left);
            } else {
                cm.types.push(CellType::Quad);
                cm.conn.extend([p00, p10, p11, p01]);
                cm.offsets.push(cm.conn.len() as u32);
                cm.label.push(d_right);
            }
        }
    }
    let mut push_edge = |a: u32, b: u32, lbl: i32| {
        cm.facet_conn.extend([a, b]);
        cm.facet_offsets.push(cm.facet_conn.len() as u32);
        cm.facet_label.push(lbl);
    };
    for iz in 0..nz {
        push_edge(idx(0, iz), idx(0, iz + 1), 1);
        push_edge(idx(nx, iz), idx(nx, iz + 1), 2);
    }
    for ix in 0..nx {
        push_edge(idx(ix, 0), idx(ix + 1, 0), 3);
        push_edge(idx(ix, nz), idx(ix + 1, nz), 4);
    }
    Mesh::D2(cm)
}

/// Analytical 1D clamped-lateral thickness-extensional Mason admittance `Y(omega)` for AlN along z.
#[allow(clippy::too_many_arguments)]
fn mason_analytical_y(
    area: f64,
    d: f64,
    omega: f64,
    c33_e: f64,
    e33: f64,
    eps33_s: f64,
    rho: f64,
    eta: f64,
    tan_d: f64,
) -> C {
    let ce = c(c33_e, c33_e * omega * eta);
    let eps_c = c(eps33_s, -eps33_s * tan_d);
    let e_c = c(e33, 0.0);
    let c33_d = ce + e_c * e_c / eps_c;
    let va = (c33_d / rho).sqrt();
    let ka = c(omega, 0.0) / va;
    let half_kd = ka * (0.5 * d);
    let tan_over_x = half_kd.tan() / half_kd;
    let kt2 = (e_c * e_c) / (c33_d * eps_c);
    let c0_c = eps_c * (area / d);
    c(0.0, omega) * c0_c / (c(1.0, 0.0) - kt2 * tan_over_x)
}

#[test]
fn mason_1d_baw_analytical_y_and_c0_match_in_3d_and_mixed_prom_vs_hdm() {
    // 1D thickness-extensional AlN plate (lateral rollers ux=0, uy=0; z-free top & bottom)
    // d = 2 um along z -> f_p ~ sqrt(c33_D / rho) / (2d) ~ 2.86 GHz
    let (lx, ly, d) = (0.4e-6, 0.4e-6, 2.0e-6);
    let area = lx * ly;
    let simplex = block(lx, ly, d, 1, 1, 16);
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();
    let mats = build_aln_materials(2.0e-4, 1.0e-3, [0.0, 0.0, 0.0]);

    let freqs_hz: Vec<f64> = (0..17).map(|k| 2.45e9 + (0.60e9) * (k as f64) / 16.0).collect();

    let make_prob = |prom_tol: Option<f64>, solver: &str| {
        let mut prob = Problem::new(
            Study::Driven {
                f_min: freqs_hz[0],
                f_max: *freqs_hz.last().unwrap(),
                n_points: freqs_hz.len(),
                freqs_hz: Some(freqs_hz.clone()),
                prom_tol,
                max_hdm_solves: 12,
                strict_prom_tol: prom_tol.is_some(),
                seed_freqs_hz: vec![],
                solver: solver.into(),
                blr_tol: 1e-8,
            },
            2,
        );
        prob.materials = mats.clone();
        // block() labels: 1=x_lo, 6=x_hi, 2=y_lo, 5=y_hi, 3=z_lo, 4=z_hi
        // Rollers on x/y lateral faces (1,6 and 2,5) to isolate pure 1D longitudinal thickness wave
        prob.bcs.push(Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } });
        prob.bcs.push(Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } });
        prob.bcs.push(Bc { labels: vec![3], kind: BcKind::Ground });
        prob.bcs.push(Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } });
        prob
    };

    let res_prom = nyx_sms::run(&mesh, &make_prob(Some(1e-7), "auto"), &femcore::progress::Progress::new()).unwrap();
    let res_hdm = nyx_sms::run(&mesh, &make_prob(None, "lblt"), &femcore::progress::Progress::new()).unwrap();
    let res_blr = nyx_sms::run(&mesh, &make_prob(Some(1e-7), "blr"), &femcore::progress::Progress::new()).unwrap();

    let data_prom = DrivenPortData::of(&res_prom).expect("DrivenPortData missing");
    let data_hdm = DrivenPortData::of(&res_hdm).expect("DrivenPortData missing");
    let data_blr = DrivenPortData::of(&res_blr).expect("DrivenPortData missing");

    let c0_expected = EPSILON_0 * 10.12 * area / d;
    let c0_got = res_prom.scalar("C_0 (F)").expect("missing C_0 (F)");
    assert!((c0_got - c0_expected).abs() / c0_expected < 1e-6, "C_0 mismatch: {c0_got} vs {c0_expected}");

    // Mason parameters: c33_E = 389e9, e33 = 1.55, eps33 = 10.12*EPS0, rho = 3260, eta = 2e-4 Pa*s, tan_d = 1e-3
    let c33_e = 389e9;
    let c33_eta = 2.0e-4 / c33_e;
    for (k, &f) in freqs_hz.iter().enumerate() {
        let w = 2.0 * PI * f;
        let y_ana = mason_analytical_y(area, d, w, c33_e, 1.55, 10.12 * EPSILON_0, 3260.0, c33_eta, 1.0e-3);
        let y_prom = c(data_prom.y_matrix[k][0][0][0], data_prom.y_matrix[k][0][0][1]);
        let y_hdm = c(data_hdm.y_matrix[k][0][0][0], data_hdm.y_matrix[k][0][0][1]);
        let y_blr = c(data_blr.y_matrix[k][0][0][0], data_blr.y_matrix[k][0][0][1]);

        assert!(y_prom.re >= -1e-12, "passivity violated at {f} Hz: Re Y = {}", y_prom.re);
        let s11 = c(data_prom.s_matrix[k][0][0][0], data_prom.s_matrix[k][0][0][1]);
        assert!(s11.norm() <= 1.0 + 1e-10, "non-passive |S11| = {}", s11.norm());

        let rel_hdm = (y_prom - y_hdm).norm() / y_hdm.norm().max(1e-30);
        assert!(rel_hdm < 1e-6, "PROM vs HDM mismatch at {f} Hz: {rel_hdm:.3e}");
        let rel_blr = (y_blr - y_hdm).norm() / y_hdm.norm().max(1e-30);
        assert!(rel_blr < 1e-5, "BLR vs HDM mismatch at {f} Hz: {rel_blr:.3e}");
        let rel_ana = (y_hdm - y_ana).norm() / y_ana.norm().max(1e-30);
        assert!(
            rel_ana < 5e-3,
            "FEM vs analytical Mason at {f} Hz: y_fem={y_hdm:?}, y_ana={y_ana:?}, err={rel_ana:.3e}"
        );
    }

    assert!(data_prom.hdm_solves < freqs_hz.len(), "PROM should take fewer HDM solves ({})", data_prom.hdm_solves);
    assert!(data_prom.max_residual <= 1e-7);

    // Native mixed Tri+Quad 2D mesh verification on exact same 1D thickness plate
    let mixed_mesh = mixed_tri_quad_rect(lx, d, 2, 2, 16, 1, 1);
    let mut prob_2d = make_prob(Some(1e-7), "direct");
    prob_2d.thickness = ly;
    prob_2d.bcs = vec![
        Bc { labels: vec![1, 2], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } },
    ];
    let res_2d = nyx_sms::run(&mixed_mesh, &prob_2d, &femcore::progress::Progress::new()).unwrap();
    let data_2d = DrivenPortData::of(&res_2d).unwrap();
    for (k, &f) in freqs_hz.iter().enumerate() {
        let y_2d = c(data_2d.y_matrix[k][0][0][0], data_2d.y_matrix[k][0][0][1]);
        let y_3d = c(data_hdm.y_matrix[k][0][0][0], data_hdm.y_matrix[k][0][0][1]);
        let diff = (y_2d - y_3d).norm() / y_3d.norm().max(1e-30);
        assert!(diff < 2e-2, "2D mixed vs 3D Mason mismatch at {f} Hz: diff={diff:.3e}");
    }
}

#[test]
fn two_port_reciprocity_passivity_and_conductor_terminal() {
    // 2-port piezoelectric block (`4 x 2 x 2 um`) with tilted AlN crystal (`euler = [17.0, 31.0, 12.0]`),
    // anisotropic loss + viscosity, boundary `Terminal` ("left_pad" on x <= 1.2 um, z = 2 um) and
    // volumetric `ConductorRole::Terminal` ("right_metal" on cells strictly x >= 2.5 um, z >= 0.8 um).
    let mut simplex = block(4.0e-6, 2.0e-6, 2.0e-6, 8, 2, 4);
    for ci in 0..simplex.n_cells() {
        let cell = simplex.cell(ci);
        if cell.iter().all(|&v| simplex.point(v as usize)[0] >= 2.5e-6 && simplex.point(v as usize)[2] >= 0.8e-6) {
            simplex.cell_labels[ci] = 2;
        }
    }
    for fi in 0..simplex.n_facets() {
        if simplex.facet_labels[fi] == 4 {
            let fv = simplex.facet(fi);
            if fv.iter().all(|&v| simplex.point(v as usize)[0] <= 1.5e-6 + 1e-12) {
                simplex.facet_labels[fi] = 10;
            }
        }
    }
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();
    let mut mats = build_aln_materials(3.0e-4, 2.0e-3, [17.0, 31.0, 12.0]);
    mats.domains.insert(2, BTreeMap::from([("E".into(), 70e9), ("nu".into(), 0.33), ("rho".into(), 2700.0)]));

    let freqs_hz: Vec<f64> = (0..9).map(|k| 1.2e9 + 0.8e9 * (k as f64) / 8.0).collect();
    let mut prob = Problem::new(
        Study::Driven {
            f_min: freqs_hz[0],
            f_max: *freqs_hz.last().unwrap(),
            n_points: freqs_hz.len(),
            freqs_hz: Some(freqs_hz.clone()),
            prom_tol: Some(1e-8),
            max_hdm_solves: 8,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-8,
        },
        1,
    );
    prob.materials = mats;
    prob.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![10], kind: BcKind::Terminal { name: "left_pad".into(), z0_ohm: 50.0 } },
        Bc {
            labels: vec![2],
            kind: BcKind::Conductor {
                role: ConductorRole::Terminal,
                value: 0.0,
                charge: 0.0,
                name: "right_metal".into(),
                z0_ohm: 75.0,
            },
        },
    ];

    let res = nyx_sms::run(&mesh, &prob, &femcore::progress::Progress::new()).unwrap();
    let data = DrivenPortData::of(&res).unwrap();
    assert_eq!(data.ports, vec!["left_pad", "right_metal"]);
    assert_eq!(data.y_matrix.len(), freqs_hz.len());

    for (k, &f) in freqs_hz.iter().enumerate() {
        let y00 = c(data.y_matrix[k][0][0][0], data.y_matrix[k][0][0][1]);
        let y01 = c(data.y_matrix[k][0][1][0], data.y_matrix[k][0][1][1]);
        let y10 = c(data.y_matrix[k][1][0][0], data.y_matrix[k][1][0][1]);
        let y11 = c(data.y_matrix[k][1][1][0], data.y_matrix[k][1][1][1]);

        let recip = (y01 - y10).norm() / y01.norm().max(y00.norm()).max(1e-30);
        assert!(recip < 1e-11, "Y reciprocity broken at {f} Hz: Y01={y01:?} vs Y10={y10:?}");

        // Check Hermitian part of Y >= 0: trace >= 0 and det(Re_H(Y)) >= 0
        let h00 = y00.re;
        let h11 = y11.re;
        let h01 = c(0.5 * (y01.re + y10.conj().re), 0.5 * (y01.im + y10.conj().im));
        let trace = h00 + h11;
        let det = h00 * h11 - h01.norm_sqr();
        let scale = y00.norm().max(y11.norm()).max(1e-30);
        assert!(trace >= -1e-10 * scale, "passivity trace failed at {f} Hz: {trace}");
        assert!(det >= -1e-10 * scale * scale, "passivity det failed at {f} Hz: {det}");

        let z01 = c(data.z_matrix[k][0][1][0], data.z_matrix[k][0][1][1]);
        let z10 = c(data.z_matrix[k][1][0][0], data.z_matrix[k][1][0][1]);
        let z_recip = (z01 - z10).norm() / z01.norm().max(1e-30);
        assert!(z_recip < 1e-11, "Z reciprocity broken at {f} Hz");
    }
}

#[test]
fn mixed_cell_pml_driven_sweep_absorbs_and_matches_hdm() {
    // 2D mixed Tri+Quad slab with AlN piezo domain 1 and right-side PML domain 2 along +x
    let mesh = mixed_tri_quad_rect(4.0e-6, 1.0e-6, 8, 8, 4, 1, 2);
    let mut mats = build_aln_materials(0.0, 5e-4, [0.0, 0.0, 0.0]);
    let aln_an = mats.aniso[&1].clone();
    mats.aniso.insert(2, aln_an);
    mats.pml.insert(
        2,
        Pml {
            x: PmlDir::Plus,
            strength: 2.5,
            order: 2.0,
            re: 1.0,
            auto_scale: true,
            clamp_outer: true,
            ..Default::default()
        },
    );
    let freqs_hz: Vec<f64> = (0..11).map(|k| 1.5e9 + 0.5e9 * (k as f64) / 10.0).collect();

    let build_prob = |prom_tol: Option<f64>| {
        let mut p = Problem::new(
            Study::Driven {
                f_min: freqs_hz[0],
                f_max: *freqs_hz.last().unwrap(),
                n_points: freqs_hz.len(),
                freqs_hz: Some(freqs_hz.clone()),
                prom_tol,
                max_hdm_solves: 10,
                strict_prom_tol: prom_tol.is_some(),
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 1e-8,
            },
            2,
        );
        p.materials = mats.clone();
        p.bcs = vec![
            Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
            Bc { labels: vec![3], kind: BcKind::Ground },
            Bc { labels: vec![4], kind: BcKind::Terminal { name: "gate".into(), z0_ohm: 50.0 } },
        ];
        p
    };

    let res_prom = nyx_sms::run(&mesh, &build_prob(Some(1e-6)), &femcore::progress::Progress::new()).unwrap();
    let res_hdm = nyx_sms::run(&mesh, &build_prob(None), &femcore::progress::Progress::new()).unwrap();
    let d_prom = DrivenPortData::of(&res_prom).unwrap();
    let d_hdm = DrivenPortData::of(&res_hdm).unwrap();
    let lf = res_prom.scalar("layer fraction 1").expect("missing layer fraction 1");
    assert!(lf > 0.0 && lf < 1.0, "unexpected layer fraction {lf}");

    for k in 0..freqs_hz.len() {
        let yp = c(d_prom.y_matrix[k][0][0][0], d_prom.y_matrix[k][0][0][1]);
        let yh = c(d_hdm.y_matrix[k][0][0][0], d_hdm.y_matrix[k][0][0][1]);
        assert!(yh.re >= -1e-12, "PML passive conductance violated: {}", yh.re);
        let rel = (yp - yh).norm() / yh.norm().max(1e-30);
        assert!(rel < 1e-5, "mixed PML PROM vs HDM mismatch: {rel:.3e}");
    }
}

#[test]
fn driven_negative_validation_and_spec_decode() {
    let simplex = block(1.0e-6, 1.0e-6, 1.0e-6, 2, 2, 4);
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();
    let mats = build_aln_materials(1e-4, 5e-4, [0.0, 0.0, 0.0]);

    let resolver = MapResolver {
        domains: BTreeMap::from([("piezo".into(), vec![1])]),
        boundaries: BTreeMap::from([
            ("bot".into(), vec![5]),
            ("top".into(), vec![6]),
            ("left".into(), vec![1]),
            ("right".into(), vec![2]),
        ]),
    };

    let mut spec = ProblemSpec::new("sms", 2, "driven");
    spec.study.params = json!({
        "f_min": 2.4e9,
        "f_max": 2.9e9,
        "n_points": 21,
        "prom_tol": 1e-5,
        "solver": "lblt"
    });
    spec.materials.insert(
        "piezo".into(),
        nyx_physics::spec::MaterialAssign::PropsTensor(BTreeMap::from([
            ("E".into(), nyx_physics::PropValue::Scalar(330e9)),
            ("nu".into(), nyx_physics::PropValue::Scalar(0.24)),
            ("rho".into(), nyx_physics::PropValue::Scalar(3260.0)),
            ("eps_r".into(), nyx_physics::PropValue::Scalar(9.5)),
            ("eta".into(), nyx_physics::PropValue::Scalar(1e-4)),
        ])),
    );
    spec.bcs = vec![
        BcSpec { on: "bot".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }) },
        BcSpec { on: "bot".into(), kind: "ground".into(), params: json!({}) },
        BcSpec { on: "top".into(), kind: "terminal".into(), params: json!({ "name": "port1", "Z_0": 50.0 }) },
    ];
    validate_spec(&Sms, &spec).unwrap();
    let sizing = Sms.sizing(&spec, None);
    assert!(sizing.wavelengths.contains_key("piezo"));
    assert_eq!(sizing.headroom, 1.0);

    // Confirm infinite_element rejected in driven
    let mut bad_ie = spec.clone();
    bad_ie.bcs.push(BcSpec {
        on: "piezo".into(),
        kind: "infinite_element".into(),
        params: json!({ "coord": "cartesian", "x": "plus" }),
    });
    assert!(decode(&bad_ie, &Ctx::new(&resolver)).is_err());

    // Confirm piezoelectric = off rejected in driven
    let mut bad_off = spec.clone();
    bad_off.params = json!({ "piezoelectric": "off" });
    assert!(Sms.solve(&mesh, &bad_off, &mut Ctx::new(&resolver)).is_err());

    // Confirm missing terminal rejected in driven
    let mut bad_no_term = spec.clone();
    bad_no_term.bcs.pop();
    assert!(decode(&bad_no_term, &Ctx::new(&resolver)).is_err());

    // Confirm overlapping terminals rejected
    let mut prob_overlap = Problem::new(
        Study::Driven {
            f_min: 1e9,
            f_max: 2e9,
            n_points: 3,
            freqs_hz: None,
            prom_tol: Some(1e-4),
            max_hdm_solves: 6,
            strict_prom_tol: false,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        1,
    );
    prob_overlap.materials = mats;
    prob_overlap.bcs = vec![
        Bc { labels: vec![5], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Terminal { name: "a".into(), z0_ohm: 50.0 } },
        Bc { labels: vec![6], kind: BcKind::Terminal { name: "b".into(), z0_ohm: 50.0 } },
    ];
    assert!(nyx_sms::run(&mesh, &prob_overlap, &femcore::progress::Progress::new()).is_err());
}
