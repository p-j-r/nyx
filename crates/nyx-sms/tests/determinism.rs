// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tier A (`<= 60k` DOFs) subprocess `0-ULP` (`to_bits()`) determinism matrix for `nyx-sms`:
//! verifies `NYX_THREADS=1` vs `NYX_THREADS=4` bit-exact equality across multigrid-CG static,
//! thermal-expansion (`nc >= 4000`), LOBPCG eigenfrequency, 3D PML eigenfrequency, and
//! axisymmetric complex shift-invert solves.

mod common;
use femcore::props::*;
use nyx_physics::{Scalar, SolveResult};
use nyx_sms::*;

fn spec(order: usize, study: Study) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.order = order;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) },
    ];
    s.study = study;
    s.solver = LinearSolver::Iterative;
    s
}

fn push_bits(v: &mut Vec<u64>, r: &SolveResult) {
    assert!(r.dofs <= 60_000, "Tier A sms problem must stay <= 60k DOFs, got {}", r.dofs);
    v.push(r.dofs as u64);
    for sol in &r.solutions {
        if let Some(val) = sol.value {
            v.push(val.to_bits());
        }
        if let Some(ev) = sol.eigenvalue {
            v.push(ev[0].to_bits());
            v.push(ev[1].to_bits());
        }
        for f in &sol.fields {
            v.extend(f.values.iter().map(|x| x.to_bits()));
        }
    }
    v.extend(r.scalars.iter().filter_map(|(_, s)| if let Scalar::Num(x) = s { Some(x.to_bits()) } else { None }));
}

fn run_subproc_pair(test_name: &str, mode: &str) {
    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-sms-det-{mode}-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args(["--exact", test_name, "--nocapture"])
            .env("NYX_THREADS", nt)
            .env_remove("NYX_DETERMINISTIC_FACTOR")
            .env("NYX_SMS_DET_MODE", mode)
            .env("NYX_SMS_DET_OUT", &out_path)
            .output()
            .expect("spawn sms determinism subprocess");
        assert!(
            out.status.success(),
            "child {mode} NYX_THREADS={nt} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        bytes
    };
    let b1 = run_child("1");
    let b4 = run_child("4");
    assert!(!b1.is_empty());
    assert_eq!(b1, b4, "nyx-sms {mode} outputs must be 0-ULP bit-identical for NYX_THREADS=1 and 4");
}

fn write_bits(bits: &[u64]) {
    let path = std::env::var("NYX_SMS_DET_OUT").unwrap();
    let bytes: Vec<u8> = bits.iter().flat_map(|w| w.to_le_bytes()).collect();
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn multigrid_cg_gives_identical_results_every_time() {
    if std::env::var("NYX_SMS_DET_MODE").ok().as_deref() == Some("mg_cg") {
        let m = common::block(1.0, 0.5, 0.5, 24, 12, 12);
        let s = spec(1, Study::Static);
        let a = common::solve(&m, &s);
        assert!(a.log.iter().any(|l| l.contains("conjugate gradients") || l.contains("multigrid")), "{:?}", a.log);
        let mut bits = Vec::new();
        push_bits(&mut bits, &a);
        write_bits(&bits);
        return;
    }
    run_subproc_pair("multigrid_cg_gives_identical_results_every_time", "mg_cg");
}

#[test]
fn lobpcg_gives_identical_results_every_time() {
    if std::env::var("NYX_SMS_DET_MODE").ok().as_deref() == Some("lobpcg") {
        let m = common::block(2.0, 0.5, 0.5, 40, 10, 10);
        let s = spec(1, Study::Eigenfrequency { n_modes: 4, f_search: 0.0 });
        let a = common::solve(&m, &s);
        let mut bits = Vec::new();
        push_bits(&mut bits, &a);
        write_bits(&bits);
        return;
    }
    run_subproc_pair("lobpcg_gives_identical_results_every_time", "lobpcg");
}

#[test]
fn the_axisymmetric_study_gives_identical_results_every_time() {
    if std::env::var("NYX_SMS_DET_MODE").ok().as_deref() == Some("axisym") {
        let flat = common::rect(2.0e-3, 0.4e-3, 20, 4);
        let mut pts = flat.points.clone();
        for i in 0..flat.n_points() {
            pts[2 * i] += 5.0e-3;
        }
        let m = femcore::mesh::SimplexMesh::new(
            2,
            pts,
            flat.cells.clone(),
            flat.cell_labels.clone(),
            flat.facets.clone(),
            flat.facet_labels.clone(),
        );
        let mut s =
            Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 4, f_search: 0.0, azimuthal_order: 1 }, 2);
        for (k, v) in [("E", 70e9), ("nu", 0.33), ("rho", 2700.0)] {
            s.materials.default.insert(k.into(), v);
        }
        s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        s.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 3.0, order: 2.0, ..Pml::default() });
        let a = common::solve(&m, &s);
        let mut bits = Vec::new();
        push_bits(&mut bits, &a);
        write_bits(&bits);
        return;
    }
    run_subproc_pair("the_axisymmetric_study_gives_identical_results_every_time", "axisym");
}

#[test]
fn thermal_load_and_3d_pml_are_bit_exact_across_nyx_threads() {
    if std::env::var("NYX_SMS_DET_MODE").ok().as_deref() == Some("thermal_pml") {
        let mut bits = Vec::new();

        // 1. Thermomechanical solve with nc >= 4000 cells (exercises parallel `thermal_load` chunking)
        let m_th = common::block(1.0, 0.5, 0.5, 16, 8, 8);
        assert!(m_th.n_cells() >= 4000, "expected >= 4000 cells for parallel thermal_load");
        let mut p_th = Problem::new(Study::Static, 1);
        for (k, v) in [("E", 169e9), ("nu", 0.28), ("rho", 2329.0), ("alpha", 2.6e-6)] {
            p_th.materials.default.insert(k.into(), v);
        }
        p_th.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(360.0), t_ref: 300.0 });
        p_th.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        p_th.solver = LinearSolver::Iterative;
        let r_th = common::solve(&m_th, &p_th);
        push_bits(&mut bits, &r_th);

        // 2. 3D PML eigenfrequency with nc >= 512 cells (exercises `pml::assemble_solid_pml` chunking)
        let m_pml = common::block_split_z(0.04, 0.03, 0.03, 0.02, 0.005, 4);
        assert!(m_pml.n_cells() >= 512, "expected >= 512 cells for parallel PML assembly");
        let mut p_pml = Problem::new(Study::Eigenfrequency { n_modes: 3, f_search: 60_000.0 }, 1);
        for (k, v) in [("E", 70e9), ("nu", 0.33), ("rho", 2700.0)] {
            p_pml.materials.default.insert(k.into(), v);
        }
        p_pml.bcs = vec![Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
        p_pml.materials.pml.insert(2, Pml { z: PmlDir::Plus, strength: 2.5, order: 2.0, ..Pml::default() });
        let r_pml = common::solve(&m_pml, &p_pml);
        push_bits(&mut bits, &r_pml);

        write_bits(&bits);
        return;
    }
    run_subproc_pair("thermal_load_and_3d_pml_are_bit_exact_across_nyx_threads", "thermal_pml");
}
