// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for Milestone M4D: coupled piezoelectric axisymmetric modes (`study = "axisymmetric_eigenfrequency"`)
//! across `SimplexMesh` (P1/P2) and `MixedDiscretisation` (`Tri+Quad` & curved quad), verifying:
//! - wurtzite `6mm` AlN disk thickness-extensional `f_p > f_s` stiffening at `m = 0` and floating electrode at `m = 0`,
//! - `r = 0` regularity (`phi = 0` on axis at `|m| >= 1`, `phi` free on axis at `m = 0`),
//! - cylindrical PML open radial disk/annulus passive radiation `Q > 0`,
//! - rotational invariance enforcement (error on trigonal `3m` LiNbO3 and tilted AlN, warning on small tilt).

mod common;

use femcore::material::Library;
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::progress::Progress;
use femcore::props::{Pml, PmlCoord, PmlDir};
use nyx_physics::{validate_spec, BcSpec, Ctx, LabelResolver, MaterialAssign, Mesh, ProblemSpec, Scalar, SolveResult};
use nyx_sms::problem::{decode, Bc, BcKind, ConductorRole, PiezoMode, Problem, Study};
use nyx_sms::Sms;
use serde_json::json;
use std::path::Path;

fn mat_lib() -> Library {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Library::load_dir(&dir).expect("load materials library")
}

fn aln_aniso(euler_deg: [f64; 3]) -> femcore::material::Anisotropy {
    let mut an = mat_lib().get("aluminum_nitride").unwrap().anisotropy(euler_deg);
    an.permittivity_im = None;
    an.viscosity = None;
    an
}

fn linbo3_aniso() -> femcore::material::Anisotropy {
    let mut an = mat_lib().get("lithium_niobate").unwrap().anisotropy([0.0; 3]);
    an.permittivity_im = None;
    an.viscosity = None;
    an
}

fn num_scalar(r: &SolveResult, prefix: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == prefix) {
        Some((_, Scalar::Num(v))) => *v,
        other => panic!("missing numeric scalar '{prefix}', got {other:?} in {:?}", r.scalars),
    }
}

fn field_vals<'a>(r: &'a SolveResult, sol_idx: usize, name: &str) -> &'a [f64] {
    r.solutions[sol_idx]
        .fields
        .iter()
        .find(|f| f.name == name)
        .map(|f| f.values.as_slice())
        .unwrap_or_else(|| panic!("missing field '{name}'"))
}

/// 2D half-plane `[r0, r1] x [0, h]` `Quad` mesh:
/// boundaries 1(`r=r0`), 2(`z=0`), 3(`z=h`), 4(`r=r1`).
fn quad_rz(r0: f64, r1: f64, h: f64, nr: usize, nz: usize) -> Mesh {
    let mut cm = CellMesh::<2>::new();
    let vid = |i: usize, j: usize| (j * (nr + 1) + i) as u32;
    for j in 0..=nz {
        let z = h * j as f64 / nz as f64;
        for i in 0..=nr {
            cm.points.push([r0 + (r1 - r0) * i as f64 / nr as f64, z]);
        }
    }
    for j in 0..nz {
        for i in 0..nr {
            cm.push_cell(CellType::Quad, &[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)], 1);
        }
    }
    for j in 0..nz {
        cm.push_facet(&[vid(0, j), vid(0, j + 1)], 1);
        cm.push_facet(&[vid(nr, j), vid(nr, j + 1)], 4);
    }
    for i in 0..nr {
        cm.push_facet(&[vid(i, 0), vid(i + 1, 0)], 2);
        cm.push_facet(&[vid(i, nz), vid(i + 1, nz)], 3);
    }
    Mesh::D2(cm)
}

/// 2-domain mixed `Tri + Quad` disk `[0, r_inner] (Tri) U [r_inner, r_outer] (Quad, curved/straight) x [0, h]`:
/// boundaries 1(`r=0`), 2(`z=0`), 3(`z=h`), 4(`r=r_outer`).
fn mixed_rz_disk(r_inner: f64, r_outer: f64, h: f64, n1r: usize, n2r: usize, nz: usize, bow_mid: bool) -> Mesh {
    let mut cm = CellMesh::<2>::new();
    let nr = n1r + n2r;
    let vid = |i: usize, j: usize| (j * (nr + 1) + i) as u32;
    for j in 0..=nz {
        let z = h * j as f64 / nz as f64;
        for i in 0..=n1r {
            cm.points.push([r_inner * i as f64 / n1r as f64, z]);
        }
        for i in 1..=n2r {
            let t = i as f64 / n2r as f64;
            let r = r_inner + (r_outer - r_inner) * t;
            let dr = if bow_mid && i < n2r { 0.02 * (r_outer - r_inner) * (t * (1.0 - t)) } else { 0.0 };
            cm.points.push([r + dr, z]);
        }
    }
    for j in 0..nz {
        for i in 0..n1r {
            cm.push_cell(CellType::Tri, &[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)], 1);
            cm.push_cell(CellType::Tri, &[vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)], 1);
        }
        for i in n1r..nr {
            cm.push_cell(CellType::Quad, &[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)], 2);
        }
    }
    for j in 0..nz {
        cm.push_facet(&[vid(0, j), vid(0, j + 1)], 1);
        cm.push_facet(&[vid(nr, j), vid(nr, j + 1)], 4);
    }
    for i in 0..nr {
        cm.push_facet(&[vid(i, 0), vid(i + 1, 0)], 2);
        cm.push_facet(&[vid(i, nz), vid(i + 1, nz)], 3);
    }
    Mesh::D2(cm)
}

#[test]
fn aln_disk_m0_open_vs_shorted_and_mixed_match_analytical_te() {
    // Thin wurtzite c-axis AlN disk (`r = 0..R`, `R >> h`) thickness-extensional frequency:
    // f_p = sqrt(C33^D / rho) / (2h), C33^D = c33 + e33^2 / eps33
    let h = 2.0e-6;
    let r_max = 10.0e-6;
    let an = aln_aniso([0.0, 0.0, 0.0]);
    let rho = 3260.0;
    let eps_0 = nyx_sms::piezo_mat::EPSILON_0;
    let eps33 = eps_0 * 10.7;
    let c33_d = 389e9 + 1.55_f64.powi(2) / eps33;
    let f_p_ana = (c33_d / rho).sqrt() / (2.0 * h);
    let f_s_ana = (389e9 / rho).sqrt() / (2.0 * h);

    let mesh_simplex = common::rect(r_max, h, 24, 6);
    let raw_simplex = Mesh::from_solver_mesh(&mesh_simplex).unwrap();

    // Shorted electrodes (both top & bottom grounded) vs Open/Floating top electrode at m = 0
    let mut p_short =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 8, f_search: 0.92 * f_s_ana, azimuthal_order: 0 }, 2);
    p_short.materials.aniso.insert(1, an.clone());
    p_short.materials.default.insert("rho".into(), rho);
    p_short.bcs = vec![Bc { labels: vec![2], kind: BcKind::Ground }, Bc { labels: vec![3], kind: BcKind::Ground }];

    let mut p_float = p_short.clone();
    p_float.study = Study::AxisymmetricEigenfrequency { n_modes: 8, f_search: 0.92 * f_p_ana, azimuthal_order: 0 };
    p_float.bcs = vec![
        Bc { labels: vec![2], kind: BcKind::Ground },
        Bc { labels: vec![3], kind: BcKind::FloatingPotential { charge: 0.0 } },
    ];

    let r_short = nyx_sms::run(&raw_simplex, &p_short, &Progress::new()).unwrap();
    let r_float = nyx_sms::run(&raw_simplex, &p_float, &Progress::new()).unwrap();

    // Pick thickness-extensional mode (largest axial energy share near analytical TE frequency)
    let pick_te = |res: &SolveResult, target: f64| -> f64 {
        (1..=res.solutions.len())
            .filter(|&j| num_scalar(res, &format!("axial share {j}")) > 0.6)
            .map(|j| num_scalar(res, &format!("f{j} (Hz)")))
            .min_by(|a, b| (a - target).abs().total_cmp(&(b - target).abs()))
            .expect("missing axial TE mode")
    };
    let f_s = pick_te(&r_short, f_s_ana);
    let f_p = pick_te(&r_float, f_p_ana);
    assert!(f_p > f_s, "expected piezoelectric stiffening f_p ({f_p:.4e}) > f_s ({f_s:.4e})");
    let rel_shift = (f_p - f_s) / f_p;
    let rel_shift_ana = (f_p_ana - f_s_ana) / f_p_ana;
    assert!(
        (rel_shift - rel_shift_ana).abs() < 0.015,
        "coupling frequency shift {rel_shift:.4} differs from analytical {rel_shift_ana:.4}"
    );
    assert!(!field_vals(&r_float, 0, "Re phi").is_empty());
    assert!(!field_vals(&r_float, 0, "|phi|").is_empty());

    // Euler twist about c-axis `[0, 0, 35]` keeps m = 0 axisymmetric frequencies invariant
    let mut p_twist = p_short.clone();
    p_twist.materials.aniso.insert(1, aln_aniso([0.0, 0.0, 35.0]));
    let r_twist = nyx_sms::run(&raw_simplex, &p_twist, &Progress::new()).unwrap();
    let f_s_twist = pick_te(&r_twist, f_s_ana);
    assert!((f_s - f_s_twist).abs() / f_s < 1e-8, "twist altered m=0 frequency: {f_s} vs {f_s_twist}");

    // Mixed Tri + Quad disk (order 2 serendipity) matches simplex frequency
    let raw_mixed = mixed_rz_disk(0.5 * r_max, r_max, h, 8, 8, 4, true);
    let mut p_mix = p_short.clone();
    p_mix.materials.aniso.insert(2, an);
    let r_mix = nyx_sms::run(&raw_mixed, &p_mix, &Progress::new()).unwrap();
    let f_s_mix = pick_te(&r_mix, f_s_ana);
    assert!((f_s - f_s_mix).abs() / f_s < 2e-2, "mixed Tri+Quad f_s ({f_s_mix:.4e}) vs simplex ({f_s:.4e})");
}

#[test]
fn axis_regularity_and_azimuthal_order_m_behaviour() {
    let h = 2.0e-6;
    let r_max = 6.0e-6;
    let mesh = common::rect(r_max, h, 10, 4);
    let raw = Mesh::from_solver_mesh(&mesh).unwrap();
    let an = aln_aniso([0.0, 0.0, 0.0]);
    let rho = 3260.0;

    let coords = mesh.dof_coordinates(1);
    let tol = nyx_sms::axisym::axis_tolerance(&mesh);
    let axis_nodes: Vec<usize> = (0..coords.len() / 2).filter(|&i| coords[2 * i] <= tol).collect();
    assert!(!axis_nodes.is_empty());

    // m = 0: phi is FREE on r = 0 (ungrounded on axis -> nonzero |phi| at axis nodes in TE mode)
    let mut p_m0 =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 6, f_search: 2.0e9, azimuthal_order: 0 }, 1);
    p_m0.materials.aniso.insert(1, an.clone());
    p_m0.materials.default.insert("rho".into(), rho);
    p_m0.bcs = vec![Bc { labels: vec![2], kind: BcKind::Ground }];
    let r_m0 = nyx_sms::run(&raw, &p_m0, &Progress::new()).unwrap();
    let max_axis_phi_m0 = (0..r_m0.solutions.len())
        .map(|s| {
            let p = field_vals(&r_m0, s, "|phi|");
            axis_nodes.iter().map(|&idx| p[idx]).fold(0.0_f64, f64::max)
        })
        .fold(0.0_f64, f64::max);
    assert!(max_axis_phi_m0 > 1e-4, "expected free phi at r=0 for m=0, got {max_axis_phi_m0}");

    // |m| >= 1 (m = 1 and m = 2): phi = 0 on r = 0
    for m_ord in [1_i64, 2_i64] {
        let mut p_m =
            Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search: 0.0, azimuthal_order: m_ord }, 1);
        p_m.piezoelectric = PiezoMode::On;
        p_m.materials.aniso.insert(1, an.clone());
        p_m.materials.default.insert("rho".into(), rho);
        let r_m = nyx_sms::run(&raw, &p_m, &Progress::new()).unwrap();
        for s in 0..r_m.solutions.len() {
            let p = field_vals(&r_m, s, "|phi|");
            for &idx in &axis_nodes {
                assert!(p[idx].abs() < 1e-12, "m={m_ord}: expected phi=0 at axis node {idx}, got {}", p[idx]);
            }
        }
    }
}

#[test]
fn cylindrical_pml_piezo_disk_gives_passive_damping() {
    // Solid AlN disk `[0, r_phys] U [r_phys, r_tot]` with cylindrical radial PML `[r_phys, r_tot]`
    let h = 2.0e-6;
    let r_phys = 8.0e-6;
    let r_pml = 14.0e-6;
    let raw_mixed = mixed_rz_disk(r_phys, r_pml, h, 8, 6, 4, false);

    let an = aln_aniso([0.0, 0.0, 0.0]);
    let mut p = Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 6, f_search: 2.5e9, azimuthal_order: 0 }, 1);
    p.materials.aniso.insert(1, an.clone());
    p.materials.aniso.insert(2, an);
    p.materials.default.insert("rho".into(), 3260.0);
    p.materials.pml.insert(
        2,
        Pml {
            coord: PmlCoord::Cylindrical {
                center: [0.0, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                r_in: r_phys,
                r_thickness: r_pml - r_phys,
                radial_on: true,
                axial: PmlDir::Off,
                z_in: 0.0,
                z_thickness: 0.0,
            },
            strength: 2.0,
            order: 2.0,
            clamp_outer: true,
            ..Default::default()
        },
    );
    p.bcs = vec![Bc { labels: vec![2], kind: BcKind::Ground }];

    let res = nyx_sms::run(&raw_mixed, &p, &Progress::new()).unwrap();
    let has_passive = (1..=res.solutions.len()).any(|j| {
        let q = num_scalar(&res, &format!("Q{j}"));
        let lf = num_scalar(&res, &format!("layer fraction {j}"));
        q.is_finite() && q > 0.0 && lf < 0.9
    });
    assert!(has_passive, "expected passive damped interior mode (Q > 0) with cylindrical PML: {:?}", res.scalars);
}

#[test]
fn axisym_rotational_invariance_and_floating_negative_checks() {
    let mesh = common::rect(4.0e-6, 1.0e-6, 6, 3);
    let raw = Mesh::from_solver_mesh(&mesh).unwrap();

    // 1. Trigonal 3m LiNbO3 -> error (> 5% deviation)
    let mut p_ln = Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 2, f_search: 0.0, azimuthal_order: 0 }, 1);
    p_ln.piezoelectric = PiezoMode::On;
    p_ln.materials.aniso.insert(1, linbo3_aniso());
    p_ln.materials.default.insert("rho".into(), 4628.0);
    let err_ln = nyx_sms::run(&raw, &p_ln, &Progress::new()).unwrap_err().to_string();
    assert!(err_ln.contains("5%") || err_ln.contains("rotational invariance"), "got {err_ln}");

    // 2. Tilted AlN by 30 deg -> error (> 5% deviation)
    let mut p_tilt = p_ln.clone();
    p_tilt.materials.aniso.insert(1, aln_aniso([0.0, 30.0, 0.0]));
    let err_tilt = nyx_sms::run(&raw, &p_tilt, &Progress::new()).unwrap_err().to_string();
    assert!(err_tilt.contains("5%") || err_tilt.contains("rotational invariance"), "got {err_tilt}");

    // 3. Small tilt (0.5 deg) -> succeeds with diagnostic warning (`1e-3 < dev <= 5%`)
    let mut p_small =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 2, f_search: 0.0, azimuthal_order: 0 }, 1);
    p_small.piezoelectric = PiezoMode::On;
    p_small.materials.aniso.insert(1, aln_aniso([0.0, 0.5, 0.0]));
    p_small.materials.default.insert("rho".into(), 3260.0);
    let r_small = nyx_sms::run(&raw, &p_small, &Progress::new()).unwrap();
    assert!(
        r_small.log.iter().any(|l| l.contains("azimuthally averaged")),
        "expected small-tilt averaging warning in log: {:?}",
        r_small.log
    );

    // 4. Reject floating_potential and floating conductor at |m| >= 1
    let mut p_m1_flt =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 2, f_search: 0.0, azimuthal_order: 1 }, 1);
    p_m1_flt.piezoelectric = PiezoMode::On;
    p_m1_flt.materials.aniso.insert(1, aln_aniso([0.0; 3]));
    p_m1_flt.materials.default.insert("rho".into(), 3260.0);
    p_m1_flt.bcs = vec![Bc { labels: vec![3], kind: BcKind::FloatingPotential { charge: 0.0 } }];
    let err_flt = nyx_sms::run(&raw, &p_m1_flt, &Progress::new()).unwrap_err().to_string();
    assert!(err_flt.contains("m = 0"), "got {err_flt}");

    let mut p_m1_cond = p_m1_flt.clone();
    p_m1_cond.bcs = vec![Bc {
        labels: vec![1],
        kind: BcKind::Conductor {
            role: ConductorRole::Floating,
            value: 0.0,
            charge: 0.0,
            name: String::new(),
            z0_ohm: 50.0,
        },
    }];
    let err_cond = nyx_sms::run(&raw, &p_m1_cond, &Progress::new()).unwrap_err().to_string();
    assert!(err_cond.contains("m = 0"), "got {err_cond}");

    // 5. Annular Quad mesh (r > 0) with m = 1 solves cleanly without r = 0 nodes
    let raw_ring = quad_rz(2.0e-6, 5.0e-6, 1.0e-6, 4, 2);
    let mut p_ring =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 2, f_search: 0.0, azimuthal_order: 1 }, 1);
    p_ring.piezoelectric = PiezoMode::On;
    p_ring.materials.aniso.insert(1, aln_aniso([0.0; 3]));
    p_ring.materials.default.insert("rho".into(), 3260.0);
    assert!(nyx_sms::run(&raw_ring, &p_ring, &Progress::new()).is_ok());

    // 6. Descriptor smoke validate with floating_potential at m = 0 succeeds and rejects at m = 1 via decode
    struct R;
    impl LabelResolver for R {
        fn domains(&self, _: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            Ok(vec![1])
        }
        fn boundaries(&self, _: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            Ok(vec![3])
        }
    }
    let mut spec = ProblemSpec::new("sms", 2, "axisymmetric_eigenfrequency");
    spec.materials.insert(
        "aln".into(),
        MaterialAssign::Props(
            [("E".into(), 350e9), ("nu".into(), 0.25), ("rho".into(), 3260.0), ("eps_r".into(), 9.5)].into(),
        ),
    );
    spec.bcs.push(BcSpec { on: "top".into(), kind: "floating_potential".into(), params: json!({}) });
    validate_spec(&Sms, &spec).unwrap();
    spec.study.params = json!({ "azimuthal_order": 1 });
    assert!(decode(&spec, &Ctx::new(&R)).is_err());
}
