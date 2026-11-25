// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for FR-1 (`Study::Bloch` fixed-frequency 3D Floquet-Bloch waveguide & grating solver):
//! 1. Uniform 3D silicon wire `bloch` vs 2D `modes` (`< 5e-5` relative error on `n_B` vs `n_eff`, plus `brillouin_zone = 2`).
//! 2. 3D subwavelength grating (SWG) waveguide `bloch` vs 3D `cavity` (`< 1e-4` relative consistency).
//! 3. Subwavelength bound vs radiative leaky grating transition with transverse `z` PML + transverse `x` periodicity.
//! 4. `nyx_analysis::polarization` on 3D `bloch` (unscaled `n_B` + PML `outside` exclusion) and 3D `cavity` (10-node P2 quadrature).

mod common;

use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, sort_nodes, GeomMeshOpts};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir};
use geom::model::Model;
use geom::select::Selector;
use nyx_ems::{assemble_bloch, bloch_fields_cells, solve_bloch_eigen, Ems, Problem, Study, C0};
use nyx_physics::{BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, Scalar, SolveResult};
use serde_json::json;

fn scalar_num(r: &SolveResult, name: &str) -> Option<f64> {
    r.scalars.iter().find_map(|(k, v)| match (k == name, v) {
        (true, Scalar::Num(x)) => Some(*x),
        _ => None,
    })
}

/// Builds a 2D silicon-on-insulator wire cross-section `[-0.85, 0.85] x [-0.65, 0.65] um` with `450 x 220 nm` Si core.
/// Domain 1 = SiO2 cladding, Domain 2 = Si core.
fn soi_wire_2d(h_clad: f64, h_core: f64) -> SimplexMesh {
    let mut m = Model::new_2d();
    let clad = m.rect("clad", (-0.85e-6, -0.65e-6), (0.85e-6, 0.65e-6));
    let core = m.rect("core", (-0.225e-6, -0.11e-6), (0.225e-6, 0.11e-6));
    m.overlay(&[("clad", clad), ("core", core)]);
    m.size("sz_core", Selector::parse("region(core)").unwrap(), Some(h_core.into()), None, None);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h_clad)).unwrap();
    assert_eq!((g.domain("clad"), g.domain("core")), (Some(1), Some(2)));
    sort_nodes(&g.mesh)
}

#[test]
fn uniform_3d_silicon_wire_bloch_matches_2d_modes_and_brillouin_zone_shift() {
    let lambda0 = 1.55e-6;
    let freq = C0 / lambda0;
    let (n_si, n_sio2) = (3.476_f64, 1.444_f64);
    let pitch = 0.08e-6;

    let mesh2 = soi_wire_2d(0.14e-6, 0.036e-6);
    let mut p2 = Problem::new(Study::Modes { freq, n_target: 2.5, evanescent: false, n_modes: 3 }, 2);
    p2.materials.set(1, "eps_r", n_sio2 * n_sio2);
    p2.materials.set(2, "eps_r", n_si * n_si);
    let res2 = common::solve(&mesh2, &p2).unwrap();
    let neff_te0 = res2.solutions[0].eigenvalue.unwrap()[0];
    let neff_tm0 = res2.solutions[1].eigenvalue.unwrap()[0];
    let ng_2d_te0 = scalar_num(&res2, "n_g 1").expect("2D modes n_g 1");

    // Extrude along +y for one uniform period Lambda = 80 nm (n = 1 layer, balanced aspect ratio)
    let nb = mesh2.facet_labels.iter().copied().max().unwrap_or(0);
    let mesh3 = extrude_y(&mesh2, &[Layer { thickness: pitch, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();
    let (y0_label, y1_label) = (nb + 1, nb + 2);

    let mut p3 = Problem::new(Study::Bloch { freq, n_target: 2.5, n_modes: 3, brillouin_zone: 1 }, 2);
    p3.materials.set(1, "eps_r", n_sio2 * n_sio2);
    p3.materials.set(2, "eps_r", n_si * n_si);
    p3.periodic = vec![Periodic { source: vec![y0_label], dest: vec![y1_label], phase: 0.0 }];

    // Exercise public 3-stage API directly
    let asm = assemble_bloch(&mesh3, &p3).unwrap();
    let solved = solve_bloch_eigen(&asm).unwrap();
    let res3 = bloch_fields_cells(&mesh3, &p3, &asm, &solved).unwrap();

    let nb_te0 = res3.solutions[0].eigenvalue.unwrap()[0];
    let nb_tm0 = res3.solutions[1].eigenvalue.unwrap()[0];
    let ng_3d_te0 = scalar_num(&res3, "n_g 1").expect("3D bloch n_g 1");

    let rel_te0 = ((nb_te0 - neff_te0) / neff_te0).abs();
    let rel_tm0 = ((nb_tm0 - neff_tm0) / neff_tm0).abs();
    let rel_ng = ((ng_3d_te0 - ng_2d_te0) / ng_2d_te0).abs();
    assert!(rel_te0 < 5e-5, "TE0 3D Bloch {nb_te0:.8} vs 2D {neff_te0:.8} (rel {rel_te0:.2e})");
    assert!(rel_tm0 < 5e-5, "TM0 3D Bloch {nb_tm0:.8} vs 2D {neff_tm0:.8} (rel {rel_tm0:.2e})");
    assert!(rel_ng < 5e-4, "TE0 n_g 3D Bloch {ng_3d_te0:.8} vs 2D {ng_2d_te0:.8} (rel {rel_ng:.2e})");

    // Check Brillouin zone m = 2 on a 3D periodic plane-wave cell with ny = 24 P2 layers along y
    // (periodic in x and z with phase = 0, resolving exp(+j 2*pi*y/Lambda) in BZ=2 to < 1e-5!)
    let bz_pitch = 0.40e-6;
    let bz_mesh = common::block(0.12e-6, bz_pitch, 0.12e-6, 1, 24, 1);
    let mut bz1_p = Problem::new(Study::Bloch { freq, n_target: n_si + 0.02, n_modes: 2, brillouin_zone: 1 }, 2);
    bz1_p.materials.set(1, "eps_r", n_si * n_si);
    bz1_p.periodic = vec![
        Periodic { source: vec![1], dest: vec![6], phase: 0.0 },
        Periodic { source: vec![3], dest: vec![4], phase: 0.0 },
        Periodic { source: vec![2], dest: vec![5], phase: 0.0 },
    ];
    let r_bz1 = common::solve(&bz_mesh, &bz1_p).unwrap();
    let nb_1 = r_bz1.solutions[0].eigenvalue.unwrap()[0];
    let ng_1 = scalar_num(&r_bz1, "n_g 1").unwrap();
    assert!(((nb_1 - n_si) / n_si).abs() < 1e-6);

    let mut bz2_p = bz1_p.clone();
    bz2_p.study = Study::Bloch { freq, n_target: n_si + 0.02, n_modes: 2, brillouin_zone: 2 };
    let r_bz2 = common::solve(&bz_mesh, &bz2_p).unwrap();
    let nb_2 = r_bz2.solutions[0].eigenvalue.unwrap()[0];
    let want_bz2 = nb_1 + lambda0 / bz_pitch;
    assert!(((nb_2 - want_bz2) / want_bz2).abs() < 5e-5, "BZ=2 n_B {nb_2:.8} vs expected {want_bz2:.8}");
    let ng_2 = scalar_num(&r_bz2, "n_g 1").unwrap();
    assert!(((ng_2 - ng_1) / ng_1).abs() < 1e-4, "n_g is Brillouin-zone invariant: {ng_2} vs {ng_1}");

    // Verify curved midpoint edge nodes reproduce straight Bloch result
    let small3 = common::block(0.6e-6, 0.2e-6, 0.4e-6, 2, 1, 2);
    let mut small_p = Problem::new(Study::Bloch { freq, n_target: 2.8, n_modes: 2, brillouin_zone: 1 }, 2);
    small_p.materials.set(1, "eps_r", n_si * n_si);
    small_p.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase: 0.0 }];
    let r_straight = common::solve(&small3, &small_p).unwrap();
    let mut curved_mesh = SimplexMesh::new(
        3,
        small3.points.clone(),
        small3.cells.clone(),
        small3.cell_labels.clone(),
        small3.facets.clone(),
        small3.facet_labels.clone(),
    );
    let edge_nodes: Vec<([u32; 2], [f64; 3])> = small3
        .edges()
        .iter()
        .map(|&(a, b)| {
            let (pa, pb) = (small3.point(a as usize), small3.point(b as usize));
            ([a, b], [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])])
        })
        .collect();
    curved_mesh.set_edge_nodes(&edge_nodes);
    let r_curved = common::solve(&curved_mesh, &small_p).unwrap();
    assert!(
        ((r_curved.solutions[0].eigenvalue.unwrap()[0] - r_straight.solutions[0].eigenvalue.unwrap()[0])
            / r_straight.solutions[0].eigenvalue.unwrap()[0])
            .abs()
            < 1e-10
    );
}

#[test]
fn swg_3d_waveguide_bloch_agrees_with_3d_cavity_and_envelope_decode() {
    let lambda0 = 1.55e-6;
    let freq = C0 / lambda0;
    let (n_si, n_sio2) = (3.476_f64, 1.444_f64);
    let pitch = 0.24e-6;

    let mesh2 = soi_wire_2d(0.22e-6, 0.07e-6);
    let nb = mesh2.facet_labels.iter().copied().max().unwrap_or(0);
    // 50% duty cycle SWG with n = 2 layers per half-pitch (4 layers along y so cavity P2 resolves exp(-j k_y y) < 3e-5)
    let mesh3 = extrude_y(
        &mesh2,
        &[
            Layer { thickness: 0.5 * pitch, n: 2, domain_offset: 0, fill: None },
            Layer { thickness: 0.5 * pitch, n: 2, domain_offset: 2, fill: None },
        ],
        0.0,
    )
    .unwrap();
    let (y0_label, y1_label) = (nb + 1, nb + 2);
    let ex_mesh = Mesh::from_solver_mesh(&mesh3).unwrap();

    struct SwgTable {
        y0: i32,
        y1: i32,
    }
    impl LabelResolver for SwgTable {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            match r {
                "si_tooth" => Ok(vec![2]),
                "sio2" => Ok(vec![1, 3, 4]),
                _ => Err(Error::invalid(format!("unknown domain {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            match s {
                "ymin" => Ok(vec![self.y0]),
                "ymax" => Ok(vec![self.y1]),
                _ => Err(Error::invalid(format!("unknown boundary {s}"))),
            }
        }
    }
    let table = SwgTable { y0: y0_label, y1: y1_label };

    // 1. Solve via ProblemSpec envelope (`study = "bloch"`)
    let mut spec_bloch = ProblemSpec::new("ems", 2, "bloch");
    spec_bloch.study.params = json!({ "freq": freq, "n_target": 1.95, "n_modes": 3, "brillouin_zone": 1 });
    spec_bloch.materials.insert("si_tooth".into(), MaterialAssign::Props([("eps_r".into(), n_si * n_si)].into()));
    spec_bloch.materials.insert("sio2".into(), MaterialAssign::Props([("eps_r".into(), n_sio2 * n_sio2)].into()));
    spec_bloch.bcs.push(BcSpec {
        on: "ymin".into(),
        kind: "periodic".into(),
        params: json!({ "dest": "ymax", "phase": 0.0 }),
    });
    let res_bloch = Ems.solve(&ex_mesh, &spec_bloch, &mut Ctx::new(&table)).unwrap();
    assert_eq!(res_bloch.study, "bloch");

    let nb_swg = res_bloch.solutions[0].eigenvalue.unwrap()[0];
    let phase_swg = scalar_num(&res_bloch, "phase 1 (rad)").expect("phase 1 (rad)");
    let ng_bloch = scalar_num(&res_bloch, "n_g 1").expect("bloch n_g 1");

    // 2. Solve `cavity` at `phase = phase_swg` and `f_search = freq`
    let mut spec_cav = ProblemSpec::new("ems", 2, "cavity");
    spec_cav.study.params = json!({ "f_search": freq, "n_modes": 3 });
    spec_cav.materials = spec_bloch.materials.clone();
    spec_cav.bcs.push(BcSpec {
        on: "ymin".into(),
        kind: "periodic".into(),
        params: json!({ "dest": "ymax", "phase": phase_swg }),
    });
    let res_cav = Ems.solve(&ex_mesh, &spec_cav, &mut Ctx::new(&table)).unwrap();
    let (idx_cav, f_cav) = res_cav
        .solutions
        .iter()
        .enumerate()
        .map(|(i, s)| (i, s.value.unwrap()))
        .min_by(|a, b| (a.1 - freq).abs().total_cmp(&(b.1 - freq).abs()))
        .unwrap();

    let rel_f = ((f_cav - freq) / freq).abs();
    assert!(rel_f < 1e-4, "3D SWG cavity frequency {f_cav:.6e} vs Bloch {freq:.6e} (rel {rel_f:.2e})");

    let nb_cav = scalar_num(&res_cav, &format!("n_B {}", idx_cav + 1)).expect("cavity n_B");
    let ng_cav = scalar_num(&res_cav, &format!("n_g {}", idx_cav + 1)).expect("cavity n_g");
    assert!(((nb_cav - nb_swg) / nb_swg).abs() < 1e-4, "cavity n_B {nb_cav} vs bloch {nb_swg}");
    assert!(((ng_cav - ng_bloch) / ng_bloch).abs() < 5e-4, "cavity n_g {ng_cav} vs bloch {ng_bloch}");

    // Also verify nyx_analysis::polarization on 3D cavity uses 10-node P2 quadrature and labels f (THz)
    let pol_cav = nyx_analysis::polarization(&ex_mesh, &res_cav).unwrap();
    assert_eq!(pol_cav.value_label, "f (THz)");
    assert!((pol_cav.rows[idx_cav].neff.0 - f_cav / 1e12).abs() < 1e-9);
}

#[test]
fn leaky_vs_lossless_transverse_pml_grating_transition_and_polarization() {
    let lambda0 = 1.55e-6;
    let freq = C0 / lambda0;
    let (n_si, n_sio2) = (3.476_f64, 1.444_f64);

    // 2D (x, z) slab section of width wx = 0.10 um (periodic along x with phase = 0):
    // domains: 1 bottom z-PML, 2 bottom SiO2, 3 Si base (150 nm), 4 Si grating tooth (70 nm), 5 top SiO2, 6 top z-PML
    let wx = 0.10e-6;
    let mut m = Model::new_2d();
    let d1 = m.rect("pml_bot", (0.0, -1.4e-6), (wx, -1.0e-6));
    let d2 = m.rect("clad_bot", (0.0, -1.0e-6), (wx, -0.11e-6));
    let d3 = m.rect("si_base", (0.0, -0.11e-6), (wx, 0.04e-6));
    let d4 = m.rect("tooth", (0.0, 0.04e-6), (wx, 0.11e-6));
    let d5 = m.rect("clad_top", (0.0, 0.11e-6), (wx, 1.0e-6));
    let d6 = m.rect("pml_top", (0.0, 1.0e-6), (wx, 1.4e-6));
    m.overlay(&[("pml_bot", d1), ("clad_bot", d2), ("si_base", d3), ("tooth", d4), ("clad_top", d5), ("pml_top", d6)]);
    m.size("sz_si", Selector::parse("region(si_base)").unwrap(), Some(0.05e-6.into()), None, None);
    m.size("sz_tooth", Selector::parse("region(tooth)").unwrap(), Some(0.05e-6.into()), None, None);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.11e-6)).unwrap();
    let xmin_b = g2.select("plane(x = 0)").unwrap();
    let xmax_b = g2.select(&format!("plane(x = {wx:e})")).unwrap();
    let sorted2 = sort_nodes(&g2.mesh);
    let nb = sorted2.facet_labels.iter().copied().max().unwrap_or(0);

    let solve_grating = |pitch: f64, ny_half: usize| -> (SimplexMesh, SolveResult) {
        let mesh3 = extrude_y(
            &sorted2,
            &[
                Layer { thickness: 0.5 * pitch, n: ny_half, domain_offset: 0, fill: None },
                Layer { thickness: 0.5 * pitch, n: ny_half, domain_offset: 10, fill: None },
            ],
            0.0,
        )
        .unwrap();
        let mut p = Problem::new(Study::Bloch { freq, n_target: 2.65, n_modes: 6, brillouin_zone: 1 }, 2);
        for off in [0, 10] {
            for &dom in &[1 + off, 2 + off, 5 + off, 6 + off] {
                p.materials.set(dom, "eps_r", n_sio2 * n_sio2);
            }
            p.materials.set(3 + off, "eps_r", n_si * n_si);
            let tooth_n = if off == 0 { n_si } else { n_sio2 };
            p.materials.set(4 + off, "eps_r", tooth_n * tooth_n);
            p.materials.pml.insert(
                1 + off,
                Pml {
                    x: PmlDir::Off,
                    y: PmlDir::Off,
                    z: PmlDir::Minus,
                    strength: 3.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
            p.materials.pml.insert(
                6 + off,
                Pml {
                    x: PmlDir::Off,
                    y: PmlDir::Off,
                    z: PmlDir::Plus,
                    strength: 3.0,
                    order: 2.0,
                    re: 1.0,
                    ..Default::default()
                },
            );
        }
        p.periodic = vec![
            Periodic { source: xmin_b.clone(), dest: xmax_b.clone(), phase: 0.0 },
            Periodic { source: vec![nb + 1], dest: vec![nb + 2], phase: 0.0 },
        ];
        let res = common::solve(&mesh3, &p).unwrap();
        (mesh3, res)
    };

    // 1. Subwavelength pitch Lambda = 240 nm (below radiation threshold: all diffraction orders evanescent in SiO2)
    let (_m_sub, r_sub) = solve_grating(0.24e-6, 1);
    let te_sub = (0..r_sub.solutions.len())
        .find(|&i| {
            scalar_num(&r_sub, &format!("Ex share {}", i + 1)).unwrap_or(0.0) > 0.8
                && scalar_num(&r_sub, &format!("layer fraction {}", i + 1)).unwrap_or(1.0) < 0.2
        })
        .expect("guided TE Bloch mode in subwavelength regime");
    let alpha_sub = scalar_num(&r_sub, &format!("alpha {} (dB/cm)", te_sub + 1)).unwrap();
    assert!(alpha_sub.abs() < 0.5, "subwavelength grating should be bound (alpha = {alpha_sub} dB/cm)");

    // 2. Radiative pitch Lambda = 620 nm (m = -1 diffraction order radiates into SiO2 cladding -> absorbed by z-PML)
    let (m_rad, r_rad) = solve_grating(0.62e-6, 2);
    let te_rad = (0..r_rad.solutions.len())
        .find(|&i| {
            scalar_num(&r_rad, &format!("Ex share {}", i + 1)).unwrap_or(0.0) > 0.8
                && scalar_num(&r_rad, &format!("layer fraction {}", i + 1)).unwrap_or(1.0) < 0.4
        })
        .expect("leaky core-guided TE Bloch mode in radiative regime");
    let alpha_rad = scalar_num(&r_rad, &format!("alpha {} (dB/cm)", te_rad + 1)).unwrap();
    let nb_im_rad = scalar_num(&r_rad, &format!("n_B im {}", te_rad + 1)).unwrap();
    assert!(nb_im_rad < -1e-5, "radiative leaky grating must have passive decay Im(n_B) < 0, got {nb_im_rad:e}");
    assert!(
        alpha_rad > 2.0 && alpha_rad > 20.0 * alpha_sub.abs().max(1e-6),
        "radiative alpha {alpha_rad:.3} dB/cm must far exceed subwavelength {alpha_sub:.3e} dB/cm"
    );

    // Check FieldFrame::outside lists all 4 PML domains and nyx_analysis::polarization displays unscaled n_B!
    let frame = nyx_physics::FieldFrame::of(&r_rad).expect("bloch field frame");
    assert_eq!(frame.outside, vec![1, 6, 11, 16]);

    let ex_rad = Mesh::from_solver_mesh(&m_rad).unwrap();
    let pol = nyx_analysis::polarization(&ex_rad, &r_rad).unwrap();
    assert_eq!(pol.value_label, "neff");
    let nb_re_rad = r_rad.solutions[te_rad].eigenvalue.unwrap()[0];
    assert!(
        (pol.rows[te_rad].neff.0 - nb_re_rad).abs() < 1e-12 && nb_re_rad > 2.0,
        "3D bloch polarization must display unscaled n_B ({nb_re_rad}), got {:?}",
        pol.rows[te_rad].neff
    );
    let rep = pol.rows[te_rad].report.as_ref().unwrap();
    assert_eq!(rep.class.label(), "quasi-TE (x)");
    assert!(rep.te_fraction > 0.95, "TE slab grating TE fraction {}", rep.te_fraction);

    // Reject y-absorbing PML on Bloch study
    let mut bad_pml = Problem::new(Study::Bloch { freq, n_target: 2.5, n_modes: 2, brillouin_zone: 1 }, 2);
    bad_pml.periodic = vec![Periodic { source: vec![nb + 1], dest: vec![nb + 2], phase: 0.0 }];
    bad_pml.materials.pml.insert(1, Pml { y: PmlDir::Plus, ..Pml::default() });
    assert!(common::solve(&m_rad, &bad_pml).is_err(), "y-PML must be refused in bloch study");

    // Reject missing y-periodic BC on Bloch study
    let no_per = Problem::new(Study::Bloch { freq, n_target: 2.5, n_modes: 2, brillouin_zone: 1 }, 2);
    assert!(common::solve(&m_rad, &no_per).is_err(), "missing y-periodic pair must be refused in bloch study");
}
