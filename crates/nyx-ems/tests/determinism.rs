// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Two-tier subprocess determinism verification for `nyx-ems`:
//! - Tier A (`<= 60k` DOFs): `NYX_THREADS=1` vs `NYX_THREADS=4` bit-exact (`0-ULP` `to_bits()`)
//!   determinism across 2D waveguide modes (real & complex PML), 3D cavity, 3D driven HDM & PROM,
//!   and 3D Bloch solves.
//! - Tier B (`> 60k` DOFs, marked `#[ignore]`): `NYX_FACTOR_THREADS=1` vs `NYX_FACTOR_THREADS=4`
//!   verifying both multi-threaded `faer` runs satisfy algebraic residual `< 1e-10` and agree on
//!   primary physical QoIs within `< 1e-10` relative error.

mod common;

use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir};
use nyx_ems::{Problem, Study, WavePortBc, C0};
use nyx_physics::{DrivenPortData, Scalar, SolveResult};
use serde_json::{json, Value};

fn push_solve_bits(out: &mut Vec<u64>, r: &SolveResult) {
    assert!(r.dofs <= 60_000, "Tier A ems problem must stay <= 60k DOFs, got {}", r.dofs);
    out.push(r.dofs as u64);
    out.push(r.solutions.len() as u64);
    for sol in &r.solutions {
        if let Some(v) = sol.value {
            out.push(v.to_bits());
        }
        if let Some(ev) = sol.eigenvalue {
            out.push(ev[0].to_bits());
            out.push(ev[1].to_bits());
        }
        for f in &sol.fields {
            out.extend(f.values.iter().map(|x| x.to_bits()));
        }
        for cf in &sol.cell_fields {
            out.extend(cf.re.iter().map(|x| x.to_bits()));
            out.extend(cf.im.iter().map(|x| x.to_bits()));
        }
    }
    for (k, s) in &r.scalars {
        if k == "rom_sweep_seconds" || k.ends_with("_seconds") {
            continue;
        }
        if let Scalar::Num(v) = s {
            out.push(v.to_bits());
        }
    }
    if let Some(dp) = DrivenPortData::of(r) {
        out.extend(dp.freqs_hz.iter().map(|f| f.to_bits()));
        for mat in [&dp.s_matrix, &dp.z_matrix, &dp.y_matrix] {
            for f_row in mat {
                for row in f_row {
                    for entry in row {
                        out.push(entry[0].to_bits());
                        out.push(entry[1].to_bits());
                    }
                }
            }
        }
    }
    if let Some(ref d) = r.diagnostics {
        out.push(d.max_algebraic_residual.to_bits());
    }
}

fn tier_a_named_cases() -> Vec<(&'static str, Vec<u64>)> {
    let mut cases: Vec<(&'static str, Vec<u64>)> = Vec::new();

    // 1. Real 2D waveguide modes (WR-90, order 2)
    {
        let mut out = Vec::new();
        let m2 = common::rect(22.86e-3, 10.16e-3, 28, 14);
        let mut p_real = Problem::new(Study::Modes { freq: 12e9, n_target: 1.0, evanescent: false, n_modes: 3 }, 2);
        p_real.materials.default.insert("eps_r".into(), 2.1);
        let r_real = common::solve(&m2, &p_real).unwrap();
        push_solve_bits(&mut out, &r_real);
        cases.push(("2d_modes_real", out));
    }

    // 2. Complex 2D waveguide modes with PML (order 2)
    {
        let mut out = Vec::new();
        let m2_pml = common::rect_split_x(0.03, 0.015, 0.02, 0.0025);
        let mut p_cplx = Problem::new(Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 3 }, 2);
        p_cplx.materials.set(1, "eps_r", 4.0);
        p_cplx.materials.set(2, "eps_r", 1.0);
        p_cplx.materials.pml.insert(2, Pml { x: PmlDir::Plus, strength: 2.5, order: 2.0, ..Pml::default() });
        let r_cplx = common::solve(&m2_pml, &p_cplx).unwrap();
        push_solve_bits(&mut out, &r_cplx);
        cases.push(("2d_modes_complex_pml", out));
    }

    let (a, l, b) = (22.86e-3, 18.0e-3, 10.16e-3);

    // 3. 3D Cavity (`> 2048` cells to exercise `cavity::assemble_complex` multi-chunk tree reduction, `<= 60k` DOFs)
    {
        let mut out = Vec::new();
        let m_cav = common::block(a, l, b, 10, 8, 5);
        assert!(m_cav.n_cells() > 2048, "expected > 2048 cells for parallel cavity assembly");
        let mut p_cav = Problem::new(Study::Cavity { n_modes: 3, f_search: 9.0e9 }, 1);
        p_cav.materials.default.insert("eps_r".into(), 1.0);
        p_cav.materials.default.insert("mu_r".into(), 1.0);
        let r_cav = common::solve(&m_cav, &p_cav).unwrap();
        push_solve_bits(&mut out, &r_cav);
        cases.push(("3d_cavity", out));
    }

    // 4 & 5. 3D Driven HDM + PROM with WavePorts
    let m_drv = common::block(a, 12.0e-3, b, 6, 4, 3);
    for (tag, prom_tol) in [("3d_driven_hdm", None), ("3d_driven_prom", Some(1e-3))] {
        let mut out = Vec::new();
        let mut p_drv = Problem::new(
            Study::Driven {
                f_min: 8.5e9,
                f_max: 10.5e9,
                n_points: 5,
                freqs_hz: None,
                prom_tol,
                max_hdm_solves: 8,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            2,
        );
        p_drv.waveports = vec![
            WavePortBc { labels: vec![2], name: "wp1".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
            WavePortBc { labels: vec![5], name: "wp2".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
        ];
        let r_drv = common::solve(&m_drv, &p_drv).unwrap();
        push_solve_bits(&mut out, &r_drv);
        cases.push((tag, out));
    }

    // 6. 3D Bloch periodic waveguide
    {
        let mut out = Vec::new();
        let m2 = common::rect(a, b, 6, 3);
        let nb = m2.facet_labels.iter().copied().max().unwrap_or(0);
        let m_bl = extrude_y(&m2, &[Layer { thickness: 0.005, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();
        let mut p_bl = Problem::new(Study::Bloch { freq: 10.0e9, n_target: 0.8, n_modes: 2, brillouin_zone: 1 }, 2);
        p_bl.periodic = vec![Periodic { source: vec![nb + 1], dest: vec![nb + 2], phase: 0.0 }];
        let r_bl = common::solve(&m_bl, &p_bl).unwrap();
        push_solve_bits(&mut out, &r_bl);
        cases.push(("3d_bloch", out));
    }

    cases
}

fn encode_named_cases(cases: &[(&'static str, Vec<u64>)]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for (name, words) in cases {
        let nb = name.as_bytes();
        bytes.extend_from_slice(&(nb.len() as u64).to_le_bytes());
        bytes.extend_from_slice(nb);
        bytes.extend_from_slice(&(words.len() as u64).to_le_bytes());
        for w in words {
            bytes.extend_from_slice(&w.to_le_bytes());
        }
    }
    bytes
}

fn decode_named_cases(bytes: &[u8]) -> Vec<(String, Vec<u64>)> {
    let mut off = 0usize;
    let take_u64 = |pos: &mut usize| -> u64 {
        let arr: [u8; 8] = bytes[*pos..*pos + 8].try_into().unwrap();
        *pos += 8;
        u64::from_le_bytes(arr)
    };
    let mut out = Vec::new();
    while off < bytes.len() {
        let nlen = take_u64(&mut off) as usize;
        let name = String::from_utf8(bytes[off..off + nlen].to_vec()).unwrap();
        off += nlen;
        let wlen = take_u64(&mut off) as usize;
        let mut words = Vec::with_capacity(wlen);
        for _ in 0..wlen {
            words.push(take_u64(&mut off));
        }
        out.push((name, words));
    }
    out
}

#[test]
fn tier_a_ems_modes_cavity_and_driven_are_bit_exact_across_nyx_threads() {
    if let Ok(out_path) = std::env::var("NYX_EMS_DET_A_OUT") {
        let bytes = encode_named_cases(&tier_a_named_cases());
        std::fs::write(&out_path, bytes).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<(String, Vec<u64>)> {
        let out_path = std::env::temp_dir().join(format!("nyx-ems-det-a-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args(["--exact", "tier_a_ems_modes_cavity_and_driven_are_bit_exact_across_nyx_threads", "--nocapture"])
            .env("NYX_THREADS", nt)
            .env_remove("NYX_DETERMINISTIC_FACTOR")
            .env("NYX_EMS_DET_A_OUT", &out_path)
            .output()
            .expect("spawn ems Tier A determinism subprocess");
        assert!(
            out.status.success(),
            "child NYX_THREADS={nt} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        decode_named_cases(&bytes)
    };

    let c1 = run_child("1");
    let c4 = run_child("4");
    assert_eq!(c1.len(), c4.len());
    assert!(!c1.is_empty());
    let mut mismatches = Vec::new();
    for ((n1, w1), (n4, w4)) in c1.iter().zip(&c4) {
        assert_eq!(n1, n4);
        if w1 != w4 {
            let first = w1.iter().zip(w4).position(|(a, b)| a != b).unwrap_or(w1.len().min(w4.len()));
            let a = w1.get(first).copied().unwrap_or(0);
            let b = w4.get(first).copied().unwrap_or(0);
            mismatches.push(format!(
                "{n1}: len {} vs {}, first diff word[{first}] 0x{a:016x} ({:?}) vs 0x{b:016x} ({:?})",
                w1.len(),
                w4.len(),
                f64::from_bits(a),
                f64::from_bits(b),
            ));
        }
    }
    assert!(mismatches.is_empty(), "nyx-ems Tier A 0-ULP mismatches:\n{}", mismatches.join("\n"));
}

struct TierBReport {
    cavity_nred: usize,
    cavity_residual: f64,
    cavity_freqs_hz: Vec<f64>,
    driven_dofs: usize,
    driven_residual: f64,
    driven_s21_mag: Vec<f64>,
}

impl TierBReport {
    fn to_json_bytes(&self) -> Vec<u8> {
        json!({
            "cavity_nred": self.cavity_nred,
            "cavity_residual": self.cavity_residual,
            "cavity_freqs_hz": self.cavity_freqs_hz,
            "driven_dofs": self.driven_dofs,
            "driven_residual": self.driven_residual,
            "driven_s21_mag": self.driven_s21_mag,
        })
        .to_string()
        .into_bytes()
    }

    fn from_json_bytes(bytes: &[u8]) -> Self {
        let v: Value = serde_json::from_slice(bytes).unwrap();
        let f64_vec =
            |arr: &Value| -> Vec<f64> { arr.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect() };
        Self {
            cavity_nred: v["cavity_nred"].as_u64().unwrap() as usize,
            cavity_residual: v["cavity_residual"].as_f64().unwrap(),
            cavity_freqs_hz: f64_vec(&v["cavity_freqs_hz"]),
            driven_dofs: v["driven_dofs"].as_u64().unwrap() as usize,
            driven_residual: v["driven_residual"].as_f64().unwrap(),
            driven_s21_mag: f64_vec(&v["driven_s21_mag"]),
        }
    }
}

fn parse_nred_from_log(log: &[String]) -> usize {
    for line in log {
        if let Some(idx) = line.find(" after constraints)") {
            if let Some(open) = line[..idx].rfind('(') {
                if let Ok(v) = line[open + 1..idx].trim().parse::<usize>() {
                    return v;
                }
            }
        }
    }
    0
}

#[test]
#[ignore = "Tier B (>60k DOFs) multi-threaded faer supernodal determinism check; run in Tier 3 with -- --ignored"]
fn tier_b_over_60k_ems_cavity_and_driven_match_residual_and_qoi_across_factor_threads() {
    if let Ok(out_path) = std::env::var("NYX_EMS_DET_B_OUT") {
        // 1. 3D Box cavity (> 60k free Nedelec DOFs so `PAR_FACTOR_MIN_N = 60_000` engages `NYX_FACTOR_THREADS`)
        let (a, l, b): (f64, f64, f64) = (0.04, 0.03, 0.02);
        let f_te101 = 0.5 * C0 * ((1.0 / a).powi(2) + (1.0 / l).powi(2)).sqrt();
        let mesh = common::block(a, l, b, 15, 11, 9);
        let mut p_cav = Problem::new(Study::Cavity { n_modes: 3, f_search: 0.9 * f_te101 }, 2);
        p_cav.materials.default.insert("eps_r".into(), 1.0);
        p_cav.materials.default.insert("mu_r".into(), 1.0);
        let r_cav = common::solve(&mesh, &p_cav).unwrap();
        let cav_nred = parse_nred_from_log(&r_cav.log);
        assert!(cav_nred > 60_000, "Tier B cavity must have > 60k reduced DOFs, got {cav_nred}");
        let cav_res = r_cav.diagnostics.as_ref().unwrap().max_algebraic_residual;
        let cav_freqs: Vec<f64> = r_cav.solutions.iter().map(|s| s.value.unwrap()).collect();

        // 2. 3D Driven waveguide (> 60k DOFs, direct HDM across 2 frequencies to test multi-threaded LU residual + S21 QoI)
        let mut p_drv = Problem::new(
            Study::Driven {
                f_min: 8.5e9,
                f_max: 9.5e9,
                n_points: 2,
                freqs_hz: None,
                prom_tol: None,
                max_hdm_solves: 4,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            2,
        );
        p_drv.waveports = vec![
            WavePortBc { labels: vec![2], name: "wp1".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
            WavePortBc { labels: vec![5], name: "wp2".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
        ];
        let r_drv = common::solve(&mesh, &p_drv).unwrap();
        let drv_res = r_drv.diagnostics.as_ref().unwrap().max_algebraic_residual;
        let dp = DrivenPortData::of(&r_drv).unwrap();
        let s21_mag: Vec<f64> = dp.s_matrix.iter().map(|m| m[1][0][0].hypot(m[1][0][1])).collect();

        let rep = TierBReport {
            cavity_nred: cav_nred,
            cavity_residual: cav_res,
            cavity_freqs_hz: cav_freqs,
            driven_dofs: r_drv.dofs,
            driven_residual: drv_res,
            driven_s21_mag: s21_mag,
        };
        std::fs::write(&out_path, rep.to_json_bytes()).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_child = |ft: &str| -> TierBReport {
        let out_path = std::env::temp_dir().join(format!("nyx-ems-det-b-{}-ft{ft}.json", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args([
                "--exact",
                "tier_b_over_60k_ems_cavity_and_driven_match_residual_and_qoi_across_factor_threads",
                "--ignored",
                "--nocapture",
            ])
            .env("NYX_THREADS", "4")
            .env("NYX_FACTOR_THREADS", ft)
            .env("NYX_EIG_TOL", "1e-12")
            .env_remove("NYX_DETERMINISTIC_FACTOR")
            .env("NYX_EMS_DET_B_OUT", &out_path)
            .output()
            .expect("spawn ems Tier B determinism subprocess");
        assert!(
            out.status.success(),
            "child NYX_FACTOR_THREADS={ft} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        TierBReport::from_json_bytes(&bytes)
    };

    let t0 = std::time::Instant::now();
    let r1 = run_child("1");
    let r4 = run_child("4");
    let elapsed = t0.elapsed();
    eprintln!(
        "Tier B (>60k DOFs, cav_nred={}, drv_dofs={}) finished in {:.2} s: res1=({:.3e}, {:.3e}), res4=({:.3e}, {:.3e})",
        r1.cavity_nred,
        r1.driven_dofs,
        elapsed.as_secs_f64(),
        r1.cavity_residual,
        r1.driven_residual,
        r4.cavity_residual,
        r4.driven_residual
    );

    assert!(r1.cavity_nred > 60_000 && r4.cavity_nred == r1.cavity_nred);
    assert!(r1.cavity_residual < 1e-10, "FT=1 cavity residual {:.3e} >= 1e-10", r1.cavity_residual);
    assert!(r4.cavity_residual < 1e-10, "FT=4 cavity residual {:.3e} >= 1e-10", r4.cavity_residual);
    assert!(r1.driven_residual < 1e-10, "FT=1 driven residual {:.3e} >= 1e-10", r1.driven_residual);
    assert!(r4.driven_residual < 1e-10, "FT=4 driven residual {:.3e} >= 1e-10", r4.driven_residual);

    assert_eq!(r1.cavity_freqs_hz.len(), r4.cavity_freqs_hz.len());
    for (i, (&f1, &f4)) in r1.cavity_freqs_hz.iter().zip(&r4.cavity_freqs_hz).enumerate() {
        let rel = (f1 - f4).abs() / f1.abs().max(1.0);
        assert!(rel < 1e-10, "cavity mode {i} relative drift {rel:.3e} >= 1e-10 ({f1} vs {f4})");
    }

    assert_eq!(r1.driven_s21_mag.len(), r4.driven_s21_mag.len());
    for (k, (&s1, &s4)) in r1.driven_s21_mag.iter().zip(&r4.driven_s21_mag).enumerate() {
        let rel = (s1 - s4).abs() / s1.abs().max(1e-12);
        assert!(rel < 1e-10, "driven |S21|[{k}] relative drift {rel:.3e} >= 1e-10 ({s1} vs {s4})");
    }
}
