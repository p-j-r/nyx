// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for Milestone M3: 2.5D coupled piezoelectric waveguide modes (`study = "modes"`),
//! covering `(u_x, u_z, u_y, phi)` assembly on both `SimplexMesh` and `MixedDiscretisation`,
//! gauge-pin skipping at `q != 0` vs `q == 0`, explicit `K_{phi u}` with simultaneous `q != 0` and
//! curvilinear PML, Floquet periodic unit cells, dielectric + viscoelastic lossy QEP, AlN Lamb
//! and Y-Z LiNbO3 SAW shorted-vs-open electromechanical coupling `K^2`, and negative validations.

mod common;

use femcore::material::{Anisotropy, Library, Stiffness, Tensor};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::props::{Periodic, Pml, PmlDir};
use femcore::testkit::rect;
use nyx_physics::{validate_spec, BcSpec, Ctx, MaterialAssign, Mesh, ProblemSpec, Scalar, SolveResult};
use nyx_sms::problem::{decode, Bc, BcKind, Plane, Problem, Study};
use nyx_sms::Sms;
use serde_json::json;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;
use std::sync::Arc;

fn num(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(x))) => *x,
        other => panic!("scalar '{name}' is not a number in {:?}: {other:?}", r.scalars),
    }
}

fn shares(r: &SolveResult, j: usize) -> [f64; 3] {
    [
        num(r, &format!("lateral share {j}")),
        num(r, &format!("vertical share {j}")),
        num(r, &format!("longitudinal share {j}")),
    ]
}

fn mat_lib() -> Library {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Library::load_dir(&dir).expect("load materials library")
}

fn aln_aniso(euler_deg: [f64; 3], tan_delta: f64, eta_iso: f64) -> Anisotropy {
    let eps_re = Tensor::Diagonal { x: 9.21, y: 9.21, z: 10.12 };
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
    Anisotropy {
        stiffness: Some(Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: 389e9, c44: 125e9 }),
        permittivity: Some(eps_re),
        permittivity_im: eps_im,
        piezo_e: Some([
            [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
            [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
            [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0],
        ]),
        viscosity: visc,
        euler_deg,
        ..Default::default()
    }
}

fn mixed_tri_quad_rect(lx: f64, lz: f64, tri_nx: usize, quad_nx: usize, nz: usize, d_left: i32, d_right: i32) -> Mesh {
    let nx = tri_nx + quad_nx;
    let dx = lx / nx as f64;
    let dz = lz / nz as f64;
    let mut cm = CellMesh::<2>::new();
    let idx = |ix: usize, iz: usize| -> u32 { (iz * (nx + 1) + ix) as u32 };
    for iz in 0..=nz {
        for ix in 0..=nx {
            cm.points.push([ix as f64 * dx, iz as f64 * dz]);
        }
    }
    for iz in 0..nz {
        for ix in 0..nx {
            let (p00, p10, p11, p01) = (idx(ix, iz), idx(ix + 1, iz), idx(ix + 1, iz + 1), idx(ix, iz + 1));
            if ix < tri_nx {
                cm.push_cell(CellType::Tri, &[p00, p10, p11], d_left);
                cm.push_cell(CellType::Tri, &[p00, p11, p01], d_left);
            } else {
                cm.push_cell(CellType::Quad, &[p00, p10, p11, p01], d_right);
            }
        }
    }
    // Labels match rect(): 1=left (x=0), 2=bottom (z=0), 3=top (z=lz), 4=right (x=lx)
    for iz in 0..nz {
        cm.push_facet(&[idx(0, iz), idx(0, iz + 1)], 1);
        cm.push_facet(&[idx(nx, iz), idx(nx, iz + 1)], 4);
    }
    for ix in 0..nx {
        cm.push_facet(&[idx(ix, 0), idx(ix + 1, 0)], 2);
        cm.push_facet(&[idx(ix, nz), idx(ix + 1, nz)], 3);
    }
    Mesh::D2(cm)
}

#[derive(Clone)]
struct MapResolver {
    doms: HashMap<String, Vec<i32>>,
    bnds: HashMap<String, Vec<i32>>,
}

impl nyx_physics::LabelResolver for MapResolver {
    fn domains(&self, name: &str) -> Result<Vec<i32>, nyx_physics::Error> {
        self.doms.get(name).cloned().ok_or_else(|| nyx_physics::Error::invalid(format!("unknown domain '{name}'")))
    }

    fn boundaries(&self, name: &str) -> Result<Vec<i32>, nyx_physics::Error> {
        self.bnds.get(name).cloned().ok_or_else(|| nyx_physics::Error::invalid(format!("unknown boundary '{name}'")))
    }
}

#[test]
fn q0_aln_plate_decouples_sh_matches_2d_piezo_plane_strain_and_emits_phi_fields() {
    let (w, h) = (80.0e-6, 20.0e-6);
    let mesh = rect(w, h, 8, 4);

    // 1. Grounded bottom electrode (2), clamped side walls (1, 4), q = 0
    let mut p_wg = Problem::new(Study::Modes { n_modes: 6, f_search: 0.0, wavevector: 0.0 }, 2);
    p_wg.materials.default.insert("rho".into(), 3260.0);
    p_wg.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
    p_wg.bcs.push(Bc { labels: vec![1, 4], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    p_wg.bcs.push(Bc { labels: vec![2], kind: BcKind::Ground });

    let r_wg = common::solve(&mesh, &p_wg);
    assert_eq!((r_wg.physics.as_str(), r_wg.study.as_str()), ("sms", "modes"));
    for sol in &r_wg.solutions {
        assert!(sol.fields.iter().any(|f| f.name == "|phi|"), "missing |phi| field");
        assert!(sol.fields.iter().any(|f| f.name == "Re phi"), "missing Re phi field");
        assert!(sol.cell_fields.iter().any(|f| f.name == "phi"), "missing cell field phi");
    }

    let mut p_ps = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    p_ps.plane = Plane::Strain;
    p_ps.materials = p_wg.materials.clone();
    p_ps.bcs = p_wg.bcs.clone();
    let r_ps = common::solve(&mesh, &p_ps);
    let f_ps: Vec<f64> = r_ps.solutions.iter().map(|s| s.value.unwrap()).collect();

    let mut inplane_wg = Vec::new();
    let mut sh_wg = Vec::new();
    for j in 1..=r_wg.solutions.len() {
        let sh = shares(&r_wg, j);
        assert!((sh.iter().sum::<f64>() - 1.0).abs() < 1e-12, "shares sum: {sh:?}");
        let f = num(&r_wg, &format!("f re {j} (Hz)"));
        if sh[2] > 0.999 {
            assert!(sh[0] + sh[1] < 1e-10, "SH mode should decouple at q=0: {sh:?}");
            sh_wg.push(f);
        } else {
            assert!(sh[2] < 1e-10, "in-plane mode should have zero u_y share at q=0: {sh:?}");
            inplane_wg.push(f);
        }
    }
    assert!(!sh_wg.is_empty(), "expected SH mode at q=0");
    assert!(inplane_wg.len() >= 3, "expected at least 3 in-plane modes, got {:?}", inplane_wg);
    for (k, &f_in) in inplane_wg.iter().take(f_ps.len()).enumerate() {
        assert!(
            (f_in - f_ps[k]).abs() < 1e-7 * f_ps[k],
            "mode {k}: waveguide piezo {f_in:.6e} Hz vs 2D Plane::Strain piezo {:.6e} Hz",
            f_ps[k]
        );
    }

    // 2. Completely ungrounded AlN plate at q = 0: gauge BFS pins 1 phi DOF cleanly
    let mut p_ungrounded = Problem::new(Study::Modes { n_modes: 4, f_search: 20e6, wavevector: 0.0 }, 2);
    p_ungrounded.materials = p_wg.materials.clone();
    p_ungrounded
        .bcs
        .push(Bc { labels: vec![1, 4], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    // Trigger piezoelectric activation cleanly via a SurfaceCharge(0.0) Neumann face (no Ground!)
    p_ungrounded.bcs.push(Bc { labels: vec![3], kind: BcKind::SurfaceCharge(0.0) });
    let r_ug = common::solve(&mesh, &p_ungrounded);
    assert_eq!(r_ug.solutions.len(), 4);
}

#[test]
fn aln_lamb_shorted_vs_open_k2_plus_minus_q_and_mixed_tri_quad() {
    // 2 um thick c-axis AlN Lamb plate with periodic lateral boundary x=0 <-> x=lx (infinite slab)
    let (lx, h) = (1.5e-6, 2.0e-6);
    let q = 1.4e6; // q * h = 2.8 rad
    let simplex = rect(lx, h, 3, 8);
    let mixed = mixed_tri_quad_rect(lx, h, 2, 2, 8, 1, 1);

    let make = |wv: f64, shorted: bool| {
        let mut p = Problem::new(Study::Modes { n_modes: 6, f_search: 2.0e9, wavevector: wv }, 2);
        p.materials.default.insert("rho".into(), 3260.0);
        p.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
        p.periodic.push(Periodic { phase: 0.0, source: vec![1], dest: vec![4] });
        if shorted {
            p.bcs.push(Bc { labels: vec![2, 3], kind: BcKind::Ground });
        } else {
            // Open circuit on both faces, NO Ground at all -> tests q != 0 skipping gauge pin
            p.bcs.push(Bc { labels: vec![2, 3], kind: BcKind::SurfaceCharge(0.0) });
        }
        p
    };

    let r_open = common::solve(&simplex, &make(q, false));
    let r_short = common::solve(&simplex, &make(q, true));
    let r_neg = common::solve(&simplex, &make(-q, false));

    let first_lamb = |r: &SolveResult| -> f64 {
        (1..=r.solutions.len())
            .filter(|&j| shares(r, j)[0] < 0.05)
            .map(|j| num(r, &format!("f re {j} (Hz)")))
            .next()
            .expect("expected at least one sagittal Lamb mode")
    };

    let f_o0 = first_lamb(&r_open);
    let f_s0 = first_lamb(&r_short);
    let f_neg0 = first_lamb(&r_neg);
    assert!((f_o0 - f_neg0).abs() < 1e-8 * f_o0, "+q {f_o0} vs -q {f_neg0}");

    let vo = 2.0 * PI * f_o0 / q;
    let vs = 2.0 * PI * f_s0 / q;
    let k2 = 2.0 * (vo - vs) / vo;
    assert!(vo > vs && k2 > 0.005 && k2 < 0.08, "AlN Lamb K^2 out of physical range: vo={vo}, vs={vs}, K^2={k2}");

    // Native MixedDiscretisation (Tri + Quad) agrees closely with Simplex Tri on both open & shorted
    let r_mix_open = nyx_sms::run(&mixed, &make(q, false), &femcore::progress::Progress::new()).unwrap();
    let r_mix_short = nyx_sms::run(&mixed, &make(q, true), &femcore::progress::Progress::new()).unwrap();
    let fm_o0 = first_lamb(&r_mix_open);
    let fm_s0 = first_lamb(&r_mix_short);
    assert!((fm_o0 - f_o0).abs() < 5e-4 * f_o0, "Mixed open {fm_o0} vs Simplex {f_o0}");
    assert!((fm_s0 - f_s0).abs() < 5e-4 * f_s0, "Mixed shorted {fm_s0} vs Simplex {f_s0}");
}

#[test]
fn linbo3_yz_saw_shorted_vs_open_electromechanical_coupling() {
    // Y-Z LiNbO3 SAW: Bunge [0, 90, 0] places crystal Z along wafer +Y (propagation) and crystal Y along wafer -Z (cut normal).
    // Surface-guided Rayleigh SAW at wavelength lambda = 4 um on a 4*lambda = 16 um deep substrate.
    let lib = mat_lib();
    let lno = lib.get("lithium_niobate").expect("materials/lithium_niobate.toml");
    let rho = lno.density.expect("LiNbO3 density");

    let lambda = 4.0e-6;
    let q = 2.0 * PI / lambda;
    let width = 1.0e-6;
    let depth = 16.0e-6;
    let mesh = rect(width, depth, 2, 24);

    let make = |shorted: bool| {
        let mut p = Problem::new(Study::Modes { n_modes: 4, f_search: 3500.0 / lambda, wavevector: q }, 2);
        p.materials.default.insert("rho".into(), rho);
        let mut an = lno.anisotropy([0.0, 90.0, 0.0]);
        an.permittivity_im = None;
        an.viscosity = None;
        p.materials.aniso.insert(1, an);
        p.periodic.push(Periodic { phase: 0.0, source: vec![1], dest: vec![4] });
        // Clamp substrate bottom (label 2, z = 0) and inspect top free surface (label 3, z = depth)
        p.bcs.push(Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
        p.bcs.push(Bc { labels: vec![2], kind: BcKind::Ground });
        if shorted {
            p.bcs.push(Bc { labels: vec![3], kind: BcKind::Ground });
        }
        p
    };

    let r_open = common::solve(&mesh, &make(false));
    let r_short = common::solve(&mesh, &make(true));
    let vo = 2.0 * PI * num(&r_open, "f re 1 (Hz)") / q;
    let vs = 2.0 * PI * num(&r_short, "f re 1 (Hz)") / q;
    let k2_saw = 2.0 * (vo - vs) / vo;
    // Literature Y-Z LiNbO3 Rayleigh SAW: v_open ~ 3488 m/s, K^2 ~ 4.6-4.9%
    assert!(vo > 3350.0 && vo < 3650.0, "Y-Z LiNbO3 SAW v_open expected near 3488 m/s, got {vo:.2} m/s");
    assert!(
        k2_saw > 0.035 && k2_saw < 0.065,
        "Y-Z LiNbO3 SAW K^2 expected ~0.048, got {k2_saw:.5} (vo={vo:.2}, vs={vs:.2})"
    );
}

#[test]
fn explicit_kpu_with_q_nonzero_pml_damping_and_validation_checks() {
    // Strip AlN core (domain 2, x in [-1, 1] um) flanked by AlN Cartesian x-PML wings
    // (domain 1 left [-3, -1] um, domain 3 right [1, 3] um) -> exercises simultaneous q != 0
    // AND complex PML Jacobian where K_{phi u} is neither K_{u phi}^T nor K_{u phi}^H.
    let h = 2.0e-6;
    let q = 1.4e6;
    let f_guess = 2.1e9;
    let parts = [
        common::part(-3.0e-6, 0.0, -1.0e-6, h, 1, Some(0.4e-6)),
        common::part(-1.0e-6, 0.0, 1.0e-6, h, 2, Some(0.3e-6)),
        common::part(1.0e-6, 0.0, 3.0e-6, h, 3, Some(0.4e-6)),
    ];
    let simplex = common::parts_mesh(&parts, 0.4e-6);

    let find_bot = || -> Vec<i32> {
        let mut s = std::collections::BTreeSet::new();
        for f in 0..simplex.n_facets() {
            let fv = simplex.facet(f);
            if fv.iter().all(|&v| simplex.point(v as usize)[1].abs() < 1e-12) {
                s.insert(simplex.facet_labels[f]);
            }
        }
        s.into_iter().collect()
    };
    let bot_labels = find_bot();

    let mut p_pml = Problem::new(Study::Modes { n_modes: 3, f_search: f_guess, wavevector: q }, 2);
    p_pml.materials.default.insert("rho".into(), 3260.0);
    for dom in [1, 2, 3] {
        p_pml.materials.aniso.insert(dom, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
    }
    p_pml.materials.pml.insert(1, Pml { x: PmlDir::Minus, strength: 3.0, auto_scale: true, ..Pml::default() });
    p_pml.materials.pml.insert(3, Pml { x: PmlDir::Plus, strength: 3.0, auto_scale: true, ..Pml::default() });
    p_pml.bcs.push(Bc { labels: bot_labels, kind: BcKind::Ground });

    let r_pml = common::solve(&simplex, &p_pml);
    assert_eq!(r_pml.solutions.len(), 3);
    for j in 1..=3 {
        let lf = num(&r_pml, &format!("layer fraction {j}"));
        let q_val = num(&r_pml, &format!("Q{j}"));
        assert!((0.0..=1.0).contains(&lf), "layer fraction {j} out of range: {lf}");
        assert!(
            q_val.is_finite() && q_val > 0.0,
            "expected passive positive Q for radiating waveguide mode, got {q_val}"
        );
    }

    // Native MixedDiscretisation Tri+Quad core + PML wings also succeeds cleanly
    let mixed = mixed_tri_quad_rect(6.0e-6, h, 3, 6, 4, 1, 2);
    let mut p_mix = Problem::new(Study::Modes { n_modes: 3, f_search: f_guess, wavevector: q }, 2);
    p_mix.materials.default.insert("rho".into(), 3260.0);
    p_mix.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
    p_mix.materials.aniso.insert(2, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
    p_mix.materials.pml.insert(1, Pml { x: PmlDir::Minus, strength: 3.0, auto_scale: true, ..Pml::default() });
    p_mix.bcs.push(Bc { labels: vec![2], kind: BcKind::Ground });
    let r_mix = nyx_sms::run(&mixed, &p_mix, &femcore::progress::Progress::new()).unwrap();
    assert_eq!(r_mix.solutions.len(), 3);
    for j in 1..=3 {
        assert!(num(&r_mix, &format!("Q{j}")) > 0.0);
    }

    // Dielectric loss + viscoelastic damping QEP mode in waveguide piezo
    let mut p_loss = Problem::new(Study::Modes { n_modes: 2, f_search: 2.0e9, wavevector: q }, 2);
    p_loss.materials.default.insert("rho".into(), 3260.0);
    p_loss.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 0.0], 5e-3, 5e-4));
    p_loss.periodic.push(Periodic { phase: 0.0, source: vec![1], dest: vec![4] });
    p_loss.bcs.push(Bc { labels: vec![2], kind: BcKind::Ground });
    let r_loss = common::solve(&rect(1.5e-6, h, 2, 6), &p_loss);
    assert!(num(&r_loss, "Q1") > 10.0 && num(&r_loss, "f im 1 (Hz)") > 0.0);

    // Negative validation: terminal rejected in modes, shorted floating electrode rejected
    let mesh_sm = Mesh::from_solver_mesh(&rect(1.5e-6, h, 2, 4)).unwrap();
    let resolver = MapResolver {
        doms: [("aln".into(), vec![1])].into(),
        bnds: [("bot".into(), vec![2]), ("top".into(), vec![3])].into(),
    };
    let mut ctx = Ctx::new(&resolver);
    ctx.library = Some(Arc::new(mat_lib()));
    let mut spec_term = ProblemSpec::new("sms", 2, "modes");
    spec_term.study.params = json!({ "n_modes": 1, "wavevector": q });
    spec_term.materials.insert("aln".into(), MaterialAssign::Named("aluminum_nitride".into()));
    spec_term.bcs = vec![BcSpec { kind: "terminal".into(), on: "top".into(), params: json!({ "name": "p1" }) }];
    assert!(validate_spec(&Sms, &spec_term).is_err() || decode(&spec_term, &ctx).is_err());

    let mut p_short = Problem::new(Study::Modes { n_modes: 1, f_search: 0.0, wavevector: q }, 2);
    p_short.materials.default.insert("rho".into(), 3260.0);
    p_short.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 0.0], 0.0, 0.0));
    p_short.bcs.push(Bc { labels: vec![2], kind: BcKind::Ground });
    p_short.bcs.push(Bc { labels: vec![2], kind: BcKind::FloatingPotential { charge: 0.0 } });
    let err = nyx_sms::run(&mesh_sm, &p_short, &femcore::progress::Progress::new()).unwrap_err();
    assert!(err.to_string().contains("shorted"), "unexpected error: {err}");
}
