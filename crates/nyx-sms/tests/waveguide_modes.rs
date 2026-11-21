// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Native 2.5D elastic waveguide modes (`study = "modes"`, displacement `u(x, z) exp(-j q y)` on a 2D cross-section):
//! - At `q = 0`, out-of-plane shear (`SH`, `u_y` only, `longitudinal share = 1`) decouples from in-plane `(u_x, u_z)` motion
//!   (`longitudinal share = 0`), matches the analytical SH cutoff `f = (c_s / 2) sqrt((m / w)^2 + (n / h)^2)` with
//!   `c_s = sqrt(c44 / rho)`, and the in-plane modes match 2D plane-strain `Study::Eigenfrequency`.
//! - `+q` and `-q` are degenerate and all component and domain energy shares sum to 1.
//! - A silicon strip (`800 x 300 nm`) on tilted R-plane sapphire (`euler_deg = [0.0, 57.6, 0.0]`) at `q = 2.3468e7 rad/m`
//!   agrees with a 3D 1-layer swept Floquet solve (`phase = q * L_y`) to `< 5e-4` relative.
//! - Absorbing layers (PML) yield complex eigenfrequencies with `Q = Re f / (2 |Im f|)` and `layer fraction` in `[0, 1]`.

mod common;

use femcore::material::{Anisotropy, Stiffness};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::sort_nodes;
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir};
use nyx_physics::*;
use nyx_sms::{Bc, BcKind, Plane, Problem, Sms, Study};
use serde_json::json;

fn num(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(x))) => *x,
        other => panic!("scalar '{name}' is not a number: {other:?}"),
    }
}

fn shares(r: &SolveResult, j: usize) -> [f64; 3] {
    [
        num(r, &format!("lateral share {j}")),
        num(r, &format!("vertical share {j}")),
        num(r, &format!("longitudinal share {j}")),
    ]
}

#[test]
fn q0_decouples_sh_closed_form_and_inplane_plane_strain_modes() {
    let (w, h) = (1.0e-3, 0.4e-3);
    let mesh = common::rect(w, h, 16, 6);
    let (e, nu, rho) = (169e9_f64, 0.22_f64, 2329.0_f64);
    let mu = e / (2.0 * (1.0 + nu));
    let cs = (mu / rho).sqrt();
    // Hold u_x = u_z = u_y = 0 on left (1) and right (4), leaving bottom (2) and top (3) traction-free:
    // the lowest pure SH mode (u_y = sin(pi x / w), independent of z) has exact frequency f_SH10 = c_s / (2 w).
    let f_sh10 = cs / (2.0 * w);

    let mut wg = Problem::new(Study::Modes { n_modes: 6, f_search: 0.0, wavevector: 0.0 }, 2);
    for (k, v) in [("E", e), ("nu", nu), ("rho", rho)] {
        wg.materials.set(1, k, v);
    }
    wg.bcs.push(Bc { labels: vec![1, 4], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    let r_wg = common::solve(&mesh, &wg);
    assert_eq!((r_wg.physics.as_str(), r_wg.study.as_str()), ("sms", "modes"));

    // Plain 2D plane-strain eigenfrequency on the same mesh and in-plane BCs
    let mut ps = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    ps.plane = Plane::Strain;
    ps.materials = wg.materials.clone();
    ps.bcs = wg.bcs.clone();
    let r_ps = common::solve(&mesh, &ps);
    let f_ps: Vec<f64> = r_ps.solutions.iter().map(|s| s.value.unwrap()).collect();

    let mut sh_freqs = vec![];
    let mut inplane_freqs = vec![];
    for j in 1..=r_wg.solutions.len() {
        let sh = shares(&r_wg, j);
        assert!((sh.iter().sum::<f64>() - 1.0).abs() < 1e-12, "mode {j} shares: {sh:?}");
        let f = num(&r_wg, &format!("f re {j} (Hz)"));
        if sh[2] > 0.999 {
            assert!(sh[0] + sh[1] < 1e-12, "SH mode {j} must be purely longitudinal: {sh:?}");
            sh_freqs.push(f);
        } else {
            assert!(sh[2] < 1e-12, "in-plane mode {j} at q=0 must have zero longitudinal share: {sh:?}");
            inplane_freqs.push(f);
        }
    }
    assert!(!sh_freqs.is_empty(), "expected at least one SH mode among the lowest modes");
    assert!(
        (sh_freqs[0] - f_sh10).abs() < 2e-4 * f_sh10,
        "SH(1,0) frequency {:.6e} Hz vs closed form {:.6e} Hz",
        sh_freqs[0],
        f_sh10
    );
    // Every in-plane mode at q = 0 must match the 2D plane-strain solve to machine precision
    for (k, &f_in) in inplane_freqs.iter().take(f_ps.len()).enumerate() {
        assert!(
            (f_in - f_ps[k]).abs() < 1e-9 * f_ps[k],
            "in-plane mode {k}: waveguide {:.9e} Hz vs plane-strain {:.9e} Hz",
            f_in,
            f_ps[k]
        );
    }
}

#[test]
fn plus_q_and_minus_q_are_degenerate_and_shares_sum_to_one() {
    let parts = [
        common::part(-0.5e-6, -0.4e-6, 0.5e-6, 0.0, 1, Some(0.2e-6)),
        common::part(-0.3e-6, 0.0, 0.3e-6, 0.25e-6, 2, Some(0.15e-6)),
    ];
    let mesh = common::parts_mesh(&parts, 0.2e-6);
    let q = 1.8e7;
    let make = |wv: f64| {
        let mut p = Problem::new(Study::Modes { n_modes: 5, f_search: 5.0e9, wavevector: wv }, 2);
        // Domain 1: tilted trigonal sapphire; Domain 2: cubic silicon
        p.materials.aniso.insert(
            1,
            Anisotropy {
                stiffness: Some(Stiffness::Trigonal {
                    c11: 496.8e9,
                    c12: 163.6e9,
                    c13: 110.9e9,
                    c14: -23.5e9,
                    c33: 498.1e9,
                    c44: 147.4e9,
                }),
                euler_deg: [0.0, 57.6, 0.0],
                ..Default::default()
            },
        );
        p.materials.set(1, "rho", 3980.0);
        p.materials.aniso.insert(
            2,
            Anisotropy {
                stiffness: Some(Stiffness::Cubic { c11: 165.6e9, c12: 63.9e9, c44: 79.5e9 }),
                euler_deg: [0.0, 0.0, 0.0],
                ..Default::default()
            },
        );
        p.materials.set(2, "rho", 2329.0);
        p
    };
    let pos = common::solve(&mesh, &make(q));
    let neg = common::solve(&mesh, &make(-q));
    assert_eq!(pos.solutions.len(), neg.solutions.len());
    for j in 1..=pos.solutions.len() {
        let fp = num(&pos, &format!("f re {j} (Hz)"));
        let fn_ = num(&neg, &format!("f re {j} (Hz)"));
        assert!((fp - fn_).abs() < 1e-9 * fp, "mode {j}: +q {fp} vs -q {fn_}");
        let shp = shares(&pos, j);
        let shn = shares(&neg, j);
        assert!((shp.iter().sum::<f64>() - 1.0).abs() < 1e-12, "component shares sum: {shp:?}");
        for k in 0..3 {
            assert!((shp[k] - shn[k]).abs() < 1e-7, "mode {j} comp {k}: +q {shp:?} vs -q {shn:?}");
        }
        let d1 = num(&pos, &format!("share in domain 1 {j}"));
        let d2 = num(&pos, &format!("share in domain 2 {j}"));
        assert!((d1 + d2 - 1.0).abs() < 1e-12, "domain shares sum: {d1} + {d2}");
        assert!(d1 > 0.0 && d2 > 0.0);
    }
}

#[test]
fn tilted_r_plane_sapphire_si_strip_matches_3d_swept_floquet() {
    // 800 x 300 nm Si strip (domain 2) on a tilted R-plane sapphire block (domain 1, euler_deg = [0, 57.6, 0])
    let parts = [
        common::part(-0.6e-6, -0.35e-6, 0.6e-6, 0.0, 1, Some(0.055e-6)),
        common::part(-0.4e-6, 0.0, 0.4e-6, 0.3e-6, 2, Some(0.038e-6)),
    ];
    let raw_2d = common::parts_mesh(&parts, 0.055e-6);
    let mesh_2d = sort_nodes(&raw_2d);

    // Fix the bottom of the sapphire substrate (z = -0.35 um) so both 2D and 3D models have identical clamping
    let z_bot = -0.35e-6;
    let mut bot_labels = std::collections::BTreeSet::new();
    for f in 0..mesh_2d.n_facets() {
        let fv = mesh_2d.facet(f);
        if fv.iter().all(|&v| (mesh_2d.point(v as usize)[1] - z_bot).abs() < 1e-12) {
            bot_labels.insert(mesh_2d.facet_labels[f]);
        }
    }
    let bot_labels: Vec<i32> = bot_labels.into_iter().collect();
    assert!(!bot_labels.is_empty());

    let q = 2.3468e7;
    let f_search = 7.45e9;
    let n_modes = 4;

    let sapphire = Anisotropy {
        stiffness: Some(Stiffness::Trigonal {
            c11: 496.8e9,
            c12: 163.6e9,
            c13: 110.9e9,
            c14: -23.5e9,
            c33: 498.1e9,
            c44: 147.4e9,
        }),
        euler_deg: [0.0, 57.6, 0.0],
        ..Default::default()
    };
    let silicon = Anisotropy {
        stiffness: Some(Stiffness::Cubic { c11: 165.6e9, c12: 63.9e9, c44: 79.5e9 }),
        euler_deg: [0.0, 0.0, 0.0],
        ..Default::default()
    };

    // 1. Native 2.5D waveguide modes solve on the 2D cross-section
    let mut p2d = Problem::new(Study::Modes { n_modes, f_search, wavevector: q }, 2);
    p2d.materials.aniso.insert(1, sapphire.clone());
    p2d.materials.set(1, "rho", 3980.0);
    p2d.materials.aniso.insert(2, silicon.clone());
    p2d.materials.set(2, "rho", 2329.0);
    p2d.bcs
        .push(Bc { labels: bot_labels.clone(), kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    let r2d = common::solve(&mesh_2d, &p2d);
    let f2d: Vec<f64> = r2d.solutions.iter().map(|s| s.value.unwrap()).collect();
    assert_eq!(f2d.len(), n_modes);

    // Verify crystal orientation is preserved (not rotationally averaged as in axisymmetric):
    // both the assembled stiffness tensor and the eigenfrequencies differ from the unrotated / averaged ones.
    let asm_wg = nyx_sms::cross_section::assemble_cross_section(
        &mesh_2d,
        &p2d,
        nyx_sms::cross_section::CrossSection::Waveguide(q),
    )
    .unwrap();
    let exact_sapphire = femcore::section::plane_stiffness(
        &femcore::material::stiffness_global(sapphire.stiffness.as_ref().unwrap(), sapphire.euler_deg),
        2,
    );
    let sap_cell = (0..mesh_2d.n_cells()).find(|&ci| mesh_2d.cell_labels[ci] == 1).unwrap();
    assert_eq!(asm_wg.c6[sap_cell], exact_sapphire);
    let avg_sapphire = nyx_sms::axisym::average_about_axis(&exact_sapphire, 1, 16);
    assert!(
        (asm_wg.c6[sap_cell][1][4] - avg_sapphire[1][4]).abs() > 1e9,
        "waveguide assembly must not average crystal stiffness about the axis"
    );
    let mut p2d_unrotated = p2d.clone();
    p2d_unrotated.materials.aniso.get_mut(&1).unwrap().euler_deg = [0.0, 0.0, 0.0];
    let r2d_unrot = common::solve(&mesh_2d, &p2d_unrotated);
    let f2d_unrot: Vec<f64> = r2d_unrot.solutions.iter().map(|s| s.value.unwrap()).collect();
    let max_rot_shift = f2d.iter().zip(&f2d_unrot).map(|(a, b)| ((a - b) / a).abs()).fold(0.0, f64::max);
    assert!(max_rot_shift > 1e-5, "crystal orientation must affect waveguide modes, max shift was {max_rot_shift:.3e}");

    // 2. 3D 1-layer swept Floquet solve on the extruded mesh with L_y = 2 nm, phase = q * L_y
    let ly = 2.0e-9;
    let nb = mesh_2d.facet_labels.iter().copied().max().unwrap();
    let mesh_3d = extrude_y(&mesh_2d, &[Layer { thickness: ly, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();
    let (ymin, ymax) = (nb + 1, nb + 2);

    let mut p3d = Problem::new(Study::Eigenfrequency { n_modes, f_search }, 2);
    p3d.materials = p2d.materials.clone();
    p3d.bcs = p2d.bcs.clone();
    p3d.periodic = vec![Periodic { source: vec![ymin], dest: vec![ymax], phase: q * ly }];
    let r3d = common::solve(&mesh_3d, &p3d);
    let f3d: Vec<f64> = r3d.solutions.iter().map(|s| s.value.unwrap()).collect();
    assert_eq!(f3d.len(), n_modes);

    for (j, (&a, &b)) in f2d.iter().zip(&f3d).enumerate() {
        let rel = ((a - b) / b).abs();
        println!("mode {}: 2.5D modes = {a:.6e} Hz, 3D swept Floquet = {b:.6e} Hz (rel diff = {rel:.3e})", j + 1);
        assert!(rel < 5e-4, "mode {}: 2.5D {a} Hz vs 3D Floquet {b} Hz (rel diff {rel:.3e} >= 5e-4)", j + 1);
    }
}

#[test]
fn pml_complex_eigenfrequencies_periodic_lateral_and_inactive_domains() {
    // Core (domain 1) flanked by a lateral PML (domain 2) and an inactive air region (domain 3, no rho)
    let parts = [
        common::part(-0.5e-6, 0.0, 0.5e-6, 0.3e-6, 1, Some(0.15e-6)),
        common::part(0.5e-6, 0.0, 0.9e-6, 0.3e-6, 2, Some(0.15e-6)),
        common::part(-0.5e-6, 0.3e-6, 0.5e-6, 0.5e-6, 3, Some(0.2e-6)),
    ];
    let mesh = common::parts_mesh(&parts, 0.15e-6);
    let mut p = Problem::new(Study::Modes { n_modes: 4, f_search: 4.0e9, wavevector: 1.5e7 }, 2);
    for d in [1, 2] {
        p.materials.set(d, "E", 169e9);
        p.materials.set(d, "nu", 0.22);
        p.materials.set(d, "rho", 2329.0);
    }
    // Domain 3 is inactive and intentionally has NO "rho", "E", or "nu"
    p.inactive = vec![3];
    p.materials
        .pml
        .insert(2, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 });
    let r = common::solve(&mesh, &p);
    let mut seen_lossy = false;
    for (j, sol) in r.solutions.iter().enumerate() {
        let [re, im] = sol.eigenvalue.expect("eigenvalue on every waveguide mode");
        assert_eq!(num(&r, &format!("f re {} (Hz)", j + 1)), re);
        assert_eq!(num(&r, &format!("f im {} (Hz)", j + 1)), im);
        let lf = num(&r, &format!("layer fraction {}", j + 1));
        assert!((0.0..=1.0).contains(&lf), "layer fraction {lf}");
        let sh = shares(&r, j + 1);
        assert!((sh.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // Only active non-PML domain 1 gets a domain share
        assert!((num(&r, &format!("share in domain 1 {}", j + 1)) - 1.0).abs() < 1e-12);
        assert!(!r.scalars.iter().any(|(n, _)| *n == format!("share in domain 2 {}", j + 1)));
        assert!(!r.scalars.iter().any(|(n, _)| *n == format!("share in domain 3 {}", j + 1)));
        if im != 0.0 {
            seen_lossy = true;
            let q_val = num(&r, &format!("Q{}", j + 1));
            assert!((q_val - re / (2.0 * im.abs())).abs() < 1e-9 * q_val);
        }
    }
    assert!(seen_lossy, "expected complex eigenfrequencies with PML");

    // Lateral Floquet-periodic pair on a 2D rectangle in `Study::Modes`
    let rect = common::rect(1.0e-6, 0.3e-6, 8, 3);
    let mut per = Problem::new(Study::Modes { n_modes: 3, f_search: 3.0e9, wavevector: 1.0e7 }, 2);
    for (k, v) in [("E", 169e9), ("nu", 0.22), ("rho", 2329.0)] {
        per.materials.set(1, k, v);
    }
    per.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 0.6 }];
    let r_per = common::solve(&rect, &per);
    assert_eq!(r_per.solutions.len(), 3);

    // Validation refusals through `Sms::solve`
    struct Table;
    impl LabelResolver for Table {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            match r {
                "bar" => Ok(vec![1]),
                _ => Err(Error::invalid(format!("no region {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            match s {
                "left" => Ok(vec![1]),
                "right" => Ok(vec![4]),
                _ => Err(Error::invalid(format!("no selector {s}"))),
            }
        }
    }
    let mut spec = ProblemSpec::new("sms", 2, "modes");
    spec.study.params = json!({ "n_modes": 3, "f_search": 0.0, "wavevector": 1.0e7 });
    spec.materials.insert(
        "bar".into(),
        MaterialAssign::Props([("E".to_string(), 169e9), ("nu".to_string(), 0.22), ("rho".to_string(), 2329.0)].into()),
    );
    // 3D mesh is refused for `modes`
    let block3d = Mesh::from_solver_mesh(&common::block(1.0, 0.5, 0.4, 3, 2, 2)).unwrap();
    assert!(Sms.solve(&block3d, &spec, &mut Ctx::new(&Table)).unwrap_err().to_string().contains("2D"));
    // Non-zero fixed constraint is refused
    let ex_rect = Mesh::from_solver_mesh(&rect).unwrap();
    let mut bad_fixed = spec.clone();
    bad_fixed.bcs.push(BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 1e-6 }) });
    assert!(Sms.solve(&ex_rect, &bad_fixed, &mut Ctx::new(&Table)).unwrap_err().to_string().contains("zero"));
    // Load BC is refused
    let mut bad_load = spec.clone();
    bad_load.bcs.push(BcSpec {
        on: "right".into(),
        kind: "load".into(),
        params: json!({ "traction": [1.0, 0.0, 0.0] }),
    });
    assert!(Sms.solve(&ex_rect, &bad_load, &mut Ctx::new(&Table)).unwrap_err().to_string().contains("static"));
    // PML absorbing along y on 2D cross-section is refused
    let mut bad_pml = spec;
    bad_pml.bcs.push(BcSpec { on: "bar".into(), kind: "pml".into(), params: json!({ "y": "plus" }) });
    assert!(Sms.solve(&ex_rect, &bad_pml, &mut Ctx::new(&Table)).unwrap_err().to_string().contains("out of the plane"));
}

#[allow(dead_code)]
fn _uses_simplex_mesh(_: &SimplexMesh) {}
