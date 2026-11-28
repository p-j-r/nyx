// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-36 verification (ADR 138): scattered-field excitations, Floquet diffraction ports and their invariants.
//!
//! 1. Periodic dielectric slab between two `diffraction_port`s: TE/TM Fresnel reflection and transmission at normal and oblique
//!    incidence against the closed-form Airy formulas; energy conservation and reciprocity of the channel S-matrix.
//! 2. Lamellar dielectric grating with propagating +-1 orders: power conservation over all Rayleigh channels and reciprocity.
//! 3. Scattered-field plane wave on the same periodic slab (ports absorbing): the scattered field above the slab is the reflected
//!    plane wave `|r| E0`, the optical theorem `P_ext - P_abs = P_sca` holds and the PML-free absorbed power vanishes.
//! 4. Mie scattering of dielectric and lossy spheres (curved tetrahedra + spherical PML) against the Lorenz-Mie series.
//! 5. Determinism of the scattered-field solve across thread counts.

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use femcore::linalg::csr::{c, C};
use femcore::mesh::fixtures::brick_tet;
use femcore::props::{Periodic, Pml, PmlCoord};
use nyx_ems::{BeamModel, DiffractionPortBc, Excitation, ExcitationSpec, Problem, Study, C0};
use nyx_physics::{CellField, DrivenPortData, Scalar, SolveResult};
use std::f64::consts::PI;

fn driven_study(freq: f64) -> Study {
    Study::Driven {
        f_min: freq,
        f_max: freq,
        n_points: 1,
        freqs_hz: Some(vec![freq]),
        prom_tol: None,
        max_hdm_solves: 16,
        strict_prom_tol: false,
        seed_freqs_hz: vec![],
        solver: "auto".into(),
        blr_tol: 1e-6,
    }
}

fn scalar(res: &SolveResult, name: &str) -> Option<f64> {
    res.scalars.iter().find_map(|(k, v)| match (k == name, v) {
        (true, Scalar::Num(x)) => Some(*x),
        _ => None,
    })
}

/// Airy reflection / transmission amplitudes of a slab (index `n2`, thickness `d`) in a medium `n1` at incidence angle `theta`.
fn slab_rt(n1: f64, n2: f64, d: f64, k0: f64, theta: f64, te: bool) -> (C, C) {
    let s1 = n1 * theta.sin();
    let cos1 = theta.cos();
    let cos2 = (1.0 - (s1 / n2).powi(2)).sqrt();
    let r12 = if te {
        (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    } else {
        (n2 * cos1 - n1 * cos2) / (n2 * cos1 + n1 * cos2)
    };
    let t12 = if te { 2.0 * n1 * cos1 / (n1 * cos1 + n2 * cos2) } else { 2.0 * n1 * cos1 / (n2 * cos1 + n1 * cos2) };
    let t21 = if te { 2.0 * n2 * cos2 / (n1 * cos1 + n2 * cos2) } else { 2.0 * n2 * cos2 / (n2 * cos1 + n1 * cos2) };
    let delta = k0 * n2 * d * cos2;
    let ph = c(0.0, -2.0 * delta).exp();
    let den = c(1.0, 0.0) - ph * (r12 * r12);
    let r = (c(r12, 0.0) - ph * r12) / den;
    let t = c(0.0, -delta).exp() * (t12 * t21) / den;
    (r, t)
}

/// Unit cell `[0, a] x [0, a] x [0, h]` with the slab `z in [z1, z2]` (label 2), periodic in x and y, ports on z = 0 / z = h.
fn slab_cell(a: f64, zs: &[f64], z1: f64, z2: f64, nxy: usize) -> femcore::mesh::SimplexMesh {
    let xs: Vec<f64> = (0..=nxy).map(|i| a * i as f64 / nxy as f64).collect();
    let mut m = femcore::adapt::work::simplex_from_adapt(&brick_tet(&xs, &xs, zs));
    for ci in 0..m.n_cells() {
        let v = m.cell(ci).to_vec();
        let zc = v.iter().map(|&p| m.point(p as usize)[2]).sum::<f64>() / 4.0;
        if zc > z1 && zc < z2 {
            m.cell_labels[ci] = 2;
        }
    }
    m
}

fn z_grid(breaks: &[f64], h: f64) -> Vec<f64> {
    let mut zs = vec![breaks[0]];
    for w in breaks.windows(2) {
        let n = ((w[1] - w[0]) / h).ceil().max(1.0) as usize;
        for i in 1..=n {
            zs.push(w[0] + (w[1] - w[0]) * i as f64 / n as f64);
        }
    }
    zs
}

fn periodic_problem(freq: f64, eps_slab: f64, phase_x: f64, phase_y: f64) -> Problem {
    let mut p = Problem::new(driven_study(freq), 2);
    p.materials.set(1, "eps_r", 1.0);
    p.materials.set(2, "eps_r", eps_slab);
    p.periodic.push(Periodic { source: vec![1], dest: vec![2], phase: phase_x });
    p.periodic.push(Periodic { source: vec![3], dest: vec![4], phase: phase_y });
    p
}

fn s_entry(d: &DrivenPortData, row: &str, col: &str) -> C {
    let i = d.ports.iter().position(|n| n == row).unwrap_or_else(|| panic!("no channel {row} in {:?}", d.ports));
    let j = d.ports.iter().position(|n| n == col).unwrap_or_else(|| panic!("no channel {col} in {:?}", d.ports));
    let v = d.s_matrix[0][i][j];
    c(v[0], v[1])
}

#[test]
fn periodic_slab_fresnel_reflection_and_transmission_via_diffraction_ports() {
    let lam = 1.0e-6;
    let freq = C0 / lam;
    let k0 = 2.0 * PI / lam;
    let (a, n_slab, d) = (0.4 * lam, 1.5, 0.3 * lam);
    let (z1, z2, h) = (0.5 * lam, 0.5 * lam + d, 1.3 * lam);
    let zs = z_grid(&[0.0, z1, z2, h], lam / 16.0);
    let mesh = slab_cell(a, &zs, z1, z2, 4);
    for &theta_deg in &[0.0_f64, 30.0] {
        let theta = theta_deg.to_radians();
        let kx = k0 * theta.sin();
        let mut p = periodic_problem(freq, n_slab * n_slab, kx * a, 0.0);
        p.diffraction_ports.push(DiffractionPortBc {
            labels: vec![5],
            name: "bot".into(),
            m_max: Some(1),
            n_max: Some(1),
        });
        p.diffraction_ports.push(DiffractionPortBc {
            labels: vec![6],
            name: "top".into(),
            m_max: Some(1),
            n_max: Some(1),
        });
        let res = nyx_ems::solve_driven(&mesh, &p).expect("periodic slab driven solve");
        let data = DrivenPortData::of(&res).expect("DrivenPortData");
        eprintln!("theta {theta_deg} deg channels {:?}", data.ports);
        let mut p_neg = periodic_problem(freq, n_slab * n_slab, -kx * a, 0.0);
        p_neg.diffraction_ports.push(DiffractionPortBc {
            labels: vec![5],
            name: "bot".into(),
            m_max: Some(1),
            n_max: Some(1),
        });
        p_neg.diffraction_ports.push(DiffractionPortBc {
            labels: vec![6],
            name: "top".into(),
            m_max: Some(1),
            n_max: Some(1),
        });
        let data_neg = DrivenPortData::of(&nyx_ems::solve_driven(&mesh, &p_neg).expect("periodic slab -kx")).unwrap();
        let mut worst = 0.0_f64;
        for (pol, te) in [("TE", true), ("TM", false)] {
            let (r, t) = slab_rt(1.0, n_slab, d, k0, theta, te);
            let s11 = s_entry(&data, &format!("top:(0,0){pol}"), &format!("top:(0,0){pol}"));
            let s21 = s_entry(&data, &format!("bot:(0,0){pol}"), &format!("top:(0,0){pol}"));
            let s12_neg = s_entry(&data_neg, &format!("top:(0,0){pol}"), &format!("bot:(0,0){pol}"));
            let (r_fem, t_fem) = (s11.norm_sqr(), s21.norm_sqr());
            eprintln!(
                "theta {theta_deg} {pol}: R fem {r_fem:.7} exact {:.7} | T fem {t_fem:.7} exact {:.7} | R+T {:.9} | |S12(-kx)-S21(+kx)| {:.2e}",
                r.norm_sqr(),
                t.norm_sqr(),
                r_fem + t_fem,
                (s12_neg - s21).norm()
            );
            worst = worst.max((r_fem - r.norm_sqr()).abs()).max((t_fem - t.norm_sqr()).abs());
            let other = if te { "TM" } else { "TE" };
            let xpol_r = s_entry(&data, &format!("top:(0,0){other}"), &format!("top:(0,0){pol}"));
            let xpol_t = s_entry(&data, &format!("bot:(0,0){other}"), &format!("top:(0,0){pol}"));
            let total_col_power = r_fem + t_fem + xpol_r.norm_sqr() + xpol_t.norm_sqr();
            assert!(
                (total_col_power - 1.0).abs() < 1e-12,
                "lossless slab must conserve total channel power to machine precision: sum |S|^2 = {total_col_power}"
            );
            assert!((r_fem + t_fem - 1.0).abs() < 1e-6, "lossless slab co-pol R + T = {}", r_fem + t_fem);
            // Exact Lorentz reciprocity S_{pq}(+k_x) = S_{qp}(-k_x)
            assert!(
                (s12_neg - s21).norm() < 1e-10,
                "Lorentz reciprocity S12(-kx) = S21(+kx) violated by {:.2e}",
                (s12_neg - s21).norm()
            );
            // cross-polarization stays negligible on a Kuhn-split laterally homogeneous slab
            assert!(
                xpol_r.norm() < 1e-3 && xpol_t.norm() < 1e-3,
                "TE/TM coupling {:.2e} on a homogeneous slab",
                xpol_r.norm().max(xpol_t.norm())
            );
        }
        assert!(worst < 1e-4, "theta {theta_deg}: Fresnel efficiencies off by {worst:.3e}");
    }
}

#[test]
fn periodic_slab_diffraction_port_hermite_prom_sweep() {
    let lam0 = 1.0e-6;
    let f0 = C0 / lam0;
    let (a, n_slab, d) = (0.35 * lam0, 1.5, 0.3 * lam0);
    let (z1, z2, h) = (0.4 * lam0, 0.4 * lam0 + d, 1.1 * lam0);
    let zs = z_grid(&[0.0, z1, z2, h], lam0 / 14.0);
    let mesh = slab_cell(a, &zs, z1, z2, 3);
    let mut p = Problem::new(
        Study::Driven {
            f_min: 0.85 * f0,
            f_max: 1.15 * f0,
            n_points: 15,
            freqs_hz: None,
            prom_tol: Some(1e-4),
            max_hdm_solves: 8,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p.materials.set(1, "eps_r", 1.0);
    p.materials.set(2, "eps_r", n_slab * n_slab);
    p.periodic.push(Periodic { source: vec![1], dest: vec![2], phase: 0.0 });
    p.periodic.push(Periodic { source: vec![3], dest: vec![4], phase: 0.0 });
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![5], name: "bot".into(), m_max: Some(1), n_max: Some(1) });
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![6], name: "top".into(), m_max: Some(1), n_max: Some(1) });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("PROM diffraction_port sweep");
    let data = DrivenPortData::of(&res).expect("DrivenPortData");
    assert!(data.hdm_solves <= 4, "expected <= 4 PROM HDM anchors, got {}", data.hdm_solves);
    let i_top = data.ports.iter().position(|n| n == "top:(0,0)TE").unwrap();
    let i_bot = data.ports.iter().position(|n| n == "bot:(0,0)TE").unwrap();
    let mut worst = 0.0_f64;
    for (k, &f) in data.freqs_hz.iter().enumerate() {
        let (r, t) = slab_rt(1.0, n_slab, d, 2.0 * PI * f / C0, 0.0, true);
        let s11 = c(data.s_matrix[k][i_top][i_top][0], data.s_matrix[k][i_top][i_top][1]);
        let s21 = c(data.s_matrix[k][i_bot][i_top][0], data.s_matrix[k][i_bot][i_top][1]);
        worst = worst.max((s11.norm_sqr() - r.norm_sqr()).abs()).max((s21.norm_sqr() - t.norm_sqr()).abs());
    }
    eprintln!("PROM diffraction_port 15-pt sweep: {} HDM solves, max Fresnel error {worst:.3e}", data.hdm_solves);
    assert!(worst < 1e-4, "PROM diffraction_port Fresnel sweep error {worst:.3e}");
}

#[test]
fn lamellar_grating_rayleigh_orders_conserve_power_and_are_reciprocal() {
    let lam = 1.0e-6;
    let freq = C0 / lam;
    let a = 1.6 * lam; // orders m = -1, 0, +1 propagate in air
    let (z1, z2, h) = (0.4 * lam, 0.7 * lam, 1.1 * lam);
    let zs = z_grid(&[0.0, z1, z2, h], lam / 12.0);
    let xs: Vec<f64> = (0..=16).map(|i| a * i as f64 / 16.0).collect();
    let ys: Vec<f64> = (0..=2).map(|i| 0.2 * lam * i as f64 / 2.0).collect();
    let mut mesh = femcore::adapt::work::simplex_from_adapt(&brick_tet(&xs, &ys, &zs));
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci).to_vec();
        let xc = v.iter().map(|&q| mesh.point(q as usize)[0]).sum::<f64>() / 4.0;
        let zc = v.iter().map(|&q| mesh.point(q as usize)[2]).sum::<f64>() / 4.0;
        if zc > z1 && zc < z2 && xc < 0.5 * a {
            mesh.cell_labels[ci] = 2;
        }
    }
    let mut p = periodic_problem(freq, 4.0, 0.0, 0.0);
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![5], name: "bot".into(), m_max: Some(3), n_max: Some(0) });
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![6], name: "top".into(), m_max: Some(3), n_max: Some(0) });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("grating driven solve");
    let data = DrivenPortData::of(&res).expect("DrivenPortData");
    let n = data.ports.len();
    eprintln!("grating channels ({n}): {:?}", data.ports);
    assert!(data.ports.iter().any(|s| s == "top:(1,0)TE") && data.ports.iter().any(|s| s == "top:(-1,0)TM"));
    assert!(!data.ports.iter().any(|s| s.contains("(2,0)")), "order 2 must be evanescent for a = 1.6 lambda");
    let s = |i: usize, j: usize| c(data.s_matrix[0][i][j][0], data.s_matrix[0][i][j][1]);
    let mut worst_cons = 0.0_f64;
    let mut worst_recip = 0.0_f64;
    for q in 0..n {
        let col: f64 = (0..n).map(|pp| s(pp, q).norm_sqr()).sum();
        worst_cons = worst_cons.max((col - 1.0).abs());
        for pp in 0..n {
            // normal incidence: reciprocity pairs channel (m) with channel (-m) of the same port family
            let name_p = &data.ports[pp];
            let name_q = &data.ports[q];
            let mirror = |nm: &str| -> String {
                let (port, rest) = nm.split_once(':').unwrap();
                let inner = &rest[1..rest.find(')').unwrap()];
                let (m, nn) = inner.split_once(',').unwrap();
                let pol = &rest[rest.find(')').unwrap() + 1..];
                format!("{port}:({},{}){pol}", -m.parse::<i32>().unwrap(), -nn.parse::<i32>().unwrap())
            };
            let (ip, iq) = (
                data.ports.iter().position(|x| *x == mirror(name_q)).unwrap(),
                data.ports.iter().position(|x| *x == mirror(name_p)).unwrap(),
            );
            // Floquet reciprocity S_{p<-q}(k_t) = S_{-q<-(-p)}(-k_t); with k_t0 = 0 both sides belong to the same solve (the
            // TE basis of order -m is the negative of that of order m, so the check is on magnitudes)
            worst_recip = worst_recip.max((s(pp, q).norm() - s(ip, iq).norm()).abs());
        }
    }
    let eff: Vec<(String, f64)> = (0..n)
        .map(|pp| {
            (data.ports[pp].clone(), s(pp, data.ports.iter().position(|x| x == "top:(0,0)TE").unwrap()).norm_sqr())
        })
        .collect();
    eprintln!("grating efficiencies for top (0,0)TE incidence: {eff:?}");
    eprintln!("grating power conservation worst |sum |S|^2 - 1| = {worst_cons:.3e}, reciprocity |.| mismatch {worst_recip:.3e}");
    assert!(worst_cons < 1e-8, "Rayleigh channels must conserve power on a lossless grating: {worst_cons:.3e}");
    assert!(worst_recip < 1e-8, "mirrored-order reciprocity violated: {worst_recip:.3e}");
}

fn cell_field<'a>(res: &'a SolveResult, sol: usize, name: &str) -> &'a CellField {
    res.solutions[sol].cell_fields.iter().find(|f| f.name == name).unwrap_or_else(|| panic!("no cell field {name}"))
}

#[test]
fn scattered_field_plane_wave_on_periodic_slab_matches_airy_and_optical_theorem() {
    let lam = 1.0e-6;
    let freq = C0 / lam;
    let k0 = 2.0 * PI / lam;
    let (a, n_slab, d) = (0.4 * lam, 1.5, 0.3 * lam);
    let (z1, z2, h) = (0.5 * lam, 0.5 * lam + d, 1.3 * lam);
    let zs = z_grid(&[0.0, z1, z2, h], lam / 16.0);
    let mesh = slab_cell(a, &zs, z1, z2, 4);
    let mut p = periodic_problem(freq, n_slab * n_slab, 0.0, 0.0);
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![5], name: "bot".into(), m_max: Some(1), n_max: Some(1) });
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![6], name: "top".into(), m_max: Some(1), n_max: Some(1) });
    let e0 = 3.0;
    p.excitations.push(ExcitationSpec {
        name: "pw".into(),
        excitation: Excitation::PlaneWave {
            k_hat: [0.0, 0.0, -1.0],
            pol: [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
            amplitude: e0,
            phase: 0.0,
            origin: [0.0, 0.0, h],
        },
        domains: vec![2],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: true,
    });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("scattered-field slab solve");
    let (r, t) = slab_rt(1.0, n_slab, d, k0, 0.0, true);
    // scattered field above the slab: the reflected wave |r| E0; below: (t - 1) E_inc in magnitude
    let sol = res.solutions.iter().position(|s| s.label.starts_with("excitation (pw)")).expect("excitation solution");
    let (esx, esy) = (cell_field(&res, sol, "Esx"), cell_field(&res, sol, "Esy"));
    let nloc = esx.re.len() / mesh.n_cells();
    let (mut top, mut bot) = (Vec::new(), Vec::new());
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci).to_vec();
        let zc = v.iter().map(|&q| mesh.point(q as usize)[2]).sum::<f64>() / 4.0;
        for k in 0..nloc {
            let i = ci * nloc + k;
            let mag = (esx.re[i].powi(2) + esx.im[i].powi(2) + esy.re[i].powi(2) + esy.im[i].powi(2)).sqrt();
            if zc > z2 + 0.05 * lam {
                top.push(mag);
            } else if zc < z1 - 0.05 * lam {
                bot.push(mag);
            }
        }
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let spread = |v: &[f64]| {
        v.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - v.iter().cloned().fold(f64::INFINITY, f64::min)
    };
    let (rt, rb) = (mean(&top) / e0, mean(&bot) / e0);
    // the slab transmits t relative to an incident wave referenced at the slab's top face; |t e^{-i k0 d} - 1| is the scattered amplitude below
    let t_below = (t * c(0.0, k0 * d).exp() - c(1.0, 0.0)).norm();
    eprintln!(
        "scattered |Es|/E0 above {rt:.6} (exact |r| = {:.6}, spread {:.2e}), below {rb:.6} (exact {t_below:.6}, spread {:.2e})",
        r.norm(),
        spread(&top) / e0,
        spread(&bot) / e0
    );
    assert!((rt - r.norm()).abs() < 2e-3, "reflected scattered amplitude {rt} vs Airy {}", r.norm());
    assert!((rb - t_below).abs() < 2e-3, "forward scattered amplitude {rb} vs Airy {t_below}");
    let p_abs = scalar(&res, "pw: P_abs (W)").expect("P_abs scalar");
    let p_ext = scalar(&res, "pw: P_ext (W)").expect("P_ext scalar");
    let intensity = e0 * e0 / (2.0 * 376.730313668);
    let p_sca_exact = intensity * a * a * (r.norm_sqr() + t_below * t_below);
    eprintln!("P_abs {p_abs:.4e} W, P_ext {p_ext:.6e} W, exact scattered power {p_sca_exact:.6e} W");
    assert!(p_abs.abs() < 1e-9 * p_ext.abs().max(1e-30), "lossless slab absorbs {p_abs:.3e} W");
    assert!(
        ((p_ext - p_sca_exact) / p_sca_exact).abs() < 5e-3,
        "optical theorem: P_ext {p_ext} vs scattered power {p_sca_exact}"
    );
}

/// Lorenz-Mie efficiencies `(Q_ext, Q_sca)` of a sphere of relative index `m` (absorbing: `Im m > 0`) and size parameter `x`.
fn mie_q(m: C, x: f64) -> (f64, f64) {
    let nstop = (x + 4.0 * x.powf(1.0 / 3.0) + 2.0).ceil() as usize;
    let mx = m * x;
    let nmx = (nstop as f64).max(mx.norm()).ceil() as usize + 16;
    let mut dn = vec![c(0.0, 0.0); nmx + 1];
    for n in (1..=nmx).rev() {
        let r = c(n as f64, 0.0) / mx;
        dn[n - 1] = r - c(1.0, 0.0) / (dn[n] + r);
    }
    let (mut psi0, mut psi1) = (x.cos(), x.sin());
    let (mut chi0, mut chi1) = (-x.sin(), x.cos());
    let mut xi1 = c(psi1, -chi1);
    let (mut qext, mut qsca) = (0.0, 0.0);
    for n in 1..=nstop {
        let nf = n as f64;
        let psi = (2.0 * nf - 1.0) * psi1 / x - psi0;
        let chi = (2.0 * nf - 1.0) * chi1 / x - chi0;
        let xi = c(psi, -chi);
        let da = dn[n] / m + nf / x;
        let db = dn[n] * m + nf / x;
        let an = (da * psi - psi1) / (da * xi - xi1);
        let bn = (db * psi - psi1) / (db * xi - xi1);
        qext += (2.0 * nf + 1.0) * (an + bn).re;
        qsca += (2.0 * nf + 1.0) * (an.norm_sqr() + bn.norm_sqr());
        psi0 = psi1;
        psi1 = psi;
        chi0 = chi1;
        chi1 = chi;
        xi1 = c(psi1, -chi1);
    }
    (2.0 * qext / (x * x), 2.0 * qsca / (x * x))
}

#[test]
fn lorenz_mie_reference_series_sanity() {
    // small-sphere Rayleigh limit Q_sca = 8/3 x^4 |(m^2-1)/(m^2+2)|^2
    let (m, x) = (c(1.5, 0.0), 0.05);
    let (qe, qs) = mie_q(m, x);
    let ray = 8.0 / 3.0 * x.powi(4) * ((m * m - 1.0) / (m * m + 2.0)).norm_sqr();
    assert!((qs - ray).abs() / ray < 1e-2 && (qe - qs).abs() < 1e-12);
    // a tabulated value: m = 1.5, x = 1 gives Q_sca = 0.2152 (Bohren-Huffman BHMIE)
    let (_, qs1) = mie_q(c(1.5, 0.0), 1.0);
    assert!((qs1 - 0.2152).abs() < 5e-4, "Q_sca(1.5, 1) = {qs1}");
}

/// Conforming curved-tetrahedral (`CurvedTet`) mesh of a spherical 3-shell domain (`0..r_core` label 1,
/// `r_core..r_air` label 2, `r_air..r_pml` label 3), either restricted to the quarter-space `x >= 0, z >= 0`
/// (`quarter = true`, for `Symmetry::PecXPmcZ`) or covering all 8 octants (`quarter = false`), with geometric
/// radial grading in the exterior air shell so large-wavelength near-fields stay quasi-isotropic.
fn mie_spherical_mesh(
    r_core: f64,
    r_air: f64,
    r_pml: f64,
    n_ang: usize,
    n_core: usize,
    n_air: usize,
    n_pml: usize,
    quarter: bool,
) -> femcore::mesh::SimplexMesh {
    use std::collections::BTreeMap;
    let mut dirs: Vec<[f64; 3]> = Vec::new();
    let mut dir_id: BTreeMap<(i32, i32, i32), usize> = BTreeMap::new();
    let mut get_dir = |ix: i32, iy: i32, iz: i32| -> usize {
        if let Some(&id) = dir_id.get(&(ix, iy, iz)) {
            return id;
        }
        let v = [ix as f64, iy as f64, iz as f64];
        let nrm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let id = dirs.len();
        dirs.push([v[0] / nrm, v[1] / nrm, v[2] / nrm]);
        dir_id.insert((ix, iy, iz), id);
        id
    };
    let x_signs: &[i32] = if quarter { &[1] } else { &[1, -1] };
    let z_signs: &[i32] = if quarter { &[1] } else { &[1, -1] };
    let mut surf_tris: Vec<[usize; 3]> = Vec::new();
    for &sx in x_signs {
        for &sy in &[1_i32, -1_i32] {
            for &sz in z_signs {
                for i in 0..n_ang {
                    for j in 0..(n_ang - i) {
                        let mut p = |ii: usize, jj: usize| -> usize {
                            let kk = n_ang - ii - jj;
                            get_dir(sx * ii as i32, sy * kk as i32, sz * jj as i32)
                        };
                        let mut t = [p(i, j), p(i + 1, j), p(i, j + 1)];
                        t.sort_unstable();
                        surf_tris.push(t);
                        if i + j + 2 <= n_ang {
                            let mut t2 = [p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)];
                            t2.sort_unstable();
                            surf_tris.push(t2);
                        }
                    }
                }
            }
        }
    }
    let mut radii = Vec::with_capacity(n_core + n_air + n_pml);
    let mut layer_dom = Vec::with_capacity(n_core + n_air + n_pml);
    for s in 1..=n_core {
        radii.push(r_core * s as f64 / n_core as f64);
        layer_dom.push(1_i32);
    }
    let ratio_air = r_air / r_core;
    for s in 1..=n_air {
        radii.push(r_core * ratio_air.powf(s as f64 / n_air as f64));
        layer_dom.push(2_i32);
    }
    for s in 1..=n_pml {
        radii.push(r_air + (r_pml - r_air) * s as f64 / n_pml as f64);
        layer_dom.push(3_i32);
    }
    let nd = dirs.len();
    let mut m = femcore::adapt::work::SimplexMesh::<3>::default();
    let mut vr: Vec<(f64, [f64; 3])> = Vec::with_capacity(1 + radii.len() * nd);
    m.points.push([0.0, 0.0, 0.0]);
    vr.push((0.0, [0.0, 0.0, 0.0]));
    for &r in &radii {
        for &u in &dirs {
            m.points.push([r * u[0], r * u[1], r * u[2]]);
            vr.push((r, u));
        }
    }
    let vid = |layer: usize, d: usize| -> u32 { (1 + layer * nd + d) as u32 };
    for &[a, b, cc] in &surf_tris {
        m.push_cell(&[0, vid(0, a), vid(0, b), vid(0, cc)], layer_dom[0]);
    }
    for l in 1..radii.len() {
        let dom = layer_dom[l];
        for &[a, b, cc] in &surf_tris {
            let (a0, b0, c0) = (vid(l - 1, a), vid(l - 1, b), vid(l - 1, cc));
            let (a1, b1, c1) = (vid(l, a), vid(l, b), vid(l, cc));
            m.push_cell(&[a0, b0, c0, c1], dom);
            m.push_cell(&[a0, b0, b1, c1], dom);
            m.push_cell(&[a0, a1, b1, c1], dom);
        }
    }
    m.orient_positive();
    for fk in m.boundary_keys() {
        let vs = [fk[0], fk[1], fk[2]];
        let lbl = if quarter && vs.iter().all(|&v| m.points[v as usize][0].abs() < 1e-12) {
            1
        } else if quarter && vs.iter().all(|&v| m.points[v as usize][2].abs() < 1e-12) {
            2
        } else {
            3
        };
        m.push_facet(&vs, lbl);
    }
    let mut sm = femcore::adapt::work::simplex_from_adapt(&m);
    let mut edge_nodes = Vec::new();
    for &(va, vb) in sm.edges() {
        let (r0, u0) = vr[va as usize];
        let (r1, u1) = vr[vb as usize];
        if r0 > 0.0 && r1 > 0.0 && u0 != u1 {
            let s = [u0[0] + u1[0], u0[1] + u1[1], u0[2] + u1[2]];
            let sn = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
            let rm = 0.5 * (r0 + r1);
            edge_nodes.push(([va, vb], [rm * s[0] / sn, rm * s[1] / sn, rm * s[2] / sn]));
        }
    }
    sm.set_edge_nodes(&edge_nodes);
    assert!(sm.has_curved(), "spherical Mie mesh must contain CurvedTet elements");
    sm
}

fn mie_case(
    ka: f64,
    eps: C,
    n_ang: usize,
    n_core: usize,
    n_air: usize,
    n_pml: usize,
) -> (f64, f64, f64, f64, SolveResult) {
    let a = 1.0;
    let lam = 2.0 * PI * a / ka;
    let r_air = a + (0.55 * lam).max(0.30 * a);
    let r_pml = r_air + (0.65 * lam).max(0.40 * a);
    let mesh = mie_spherical_mesh(a, r_air, r_pml, n_ang, n_core, n_air, n_pml, true);
    let freq = ka * C0 / (2.0 * PI * a);
    let mut p = Problem::new(driven_study(freq), 2);
    p.symmetry = nyx_ems::Symmetry::PecXPmcZ;
    p.materials.set(1, "eps_r", eps.re);
    if eps.im != 0.0 {
        // exp(+i w t): eps = eps_re - i eps_im_loss
        p.materials.set(1, "eps_r_im", -eps.im);
    }
    p.materials.set(2, "eps_r", 1.0);
    p.materials.set(3, "eps_r", 1.0);
    p.materials.pml.insert(
        3,
        Pml {
            coord: PmlCoord::Spherical { center: [0.0; 3], r_in: r_air, thickness: r_pml - r_air },
            strength: 5.0,
            order: 2.0,
            re: 1.0,
            profile: femcore::props::PmlProfile::KappaPolynomial { kappa_max: 2.0 },
            ..Pml::default()
        },
    );
    // x-polarized wave propagating along +y has PEC on x = 0 and PMC on z = 0 (`Symmetry::PecXPmcZ`)
    p.excitations.push(ExcitationSpec {
        name: "mie".into(),
        excitation: Excitation::PlaneWave {
            k_hat: [0.0, 1.0, 0.0],
            pol: [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
            amplitude: 1.0,
            phase: 0.0,
            origin: [0.0; 3],
        },
        domains: vec![1],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: false,
    });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("Mie quarter-space driven solve");
    let geo = PI * a * a;
    let s_ext = scalar(&res, "mie: sigma_ext (m^2)").unwrap() / geo;
    let s_sca = scalar(&res, "mie: sigma_sca (m^2)").unwrap() / geo;
    let s_pml = scalar(&res, "mie: sigma_sca_pml (m^2)").unwrap() / geo;
    let s_abs = scalar(&res, "mie: sigma_abs (m^2)").unwrap() / geo;
    (s_ext, s_sca, s_pml, s_abs, res)
}

#[test]
fn mie_dielectric_and_lossy_spheres_match_lorenz_mie() {
    let cases = [
        (0.5, c(2.25, 0.0), 6, 3, 7, 5),
        (1.0, c(2.25, 0.0), 6, 4, 6, 5),
        (1.0, c(4.0, -1.0), 6, 5, 6, 5),
        (2.0, c(2.25, 0.0), 7, 5, 5, 5),
        (2.0, c(2.5, -0.5), 7, 5, 5, 5),
        (5.0, c(1.69, 0.0), 10, 7, 4, 5),
        (10.0, c(1.44, -0.1), 12, 9, 3, 4),
    ];
    for (ka, eps, n_ang, n_core, n_air, n_pml) in cases {
        let t0 = std::time::Instant::now();
        let m_bh = eps.sqrt().conj();
        let (qe, qs) = mie_q(m_bh, ka);
        let (fe, fs, fp, fa, res) = mie_case(ka, eps, n_ang, n_core, n_air, n_pml);
        let (ee, es) = ((fe - qe) / qe, (fs - qs) / qs);
        eprintln!(
            "Mie ka={ka} eps={eps}: Q_ext fem {fe:.6} exact {qe:.6} ({ee:+.3e}) | Q_sca fem {fs:.6} pml {fp:.6} exact {qs:.6} ({es:+.3e}) | Q_abs {fa:.6} exact {:.6} | {} dofs, {:.2} s",
            qe - qs,
            res.dofs,
            t0.elapsed().as_secs_f64()
        );
        assert!(ee.abs() < 5e-3 && es.abs() < 5e-3, "Mie ka={ka}: Q_ext err {ee:.3e}, Q_sca err {es:.3e} (> 0.5%)");
    }
}

#[test]
fn full_sphere_mie_scattered_field_feeds_stratton_chu_nffft_cross_section() {
    use nyx_physics::Mesh;
    let a = 1.0;
    let ka = 2.0;
    let lam = 2.0 * PI * a / ka;
    let r_air = a + 0.40 * lam;
    let r_pml = r_air + 0.55 * lam;
    let mut mesh = mie_spherical_mesh(a, r_air, r_pml, 5, 4, 4, 4, false);
    // Split air into inner air (label 2, r <= r_mid) and outer air (label 4) so we can also test an interior NFFFT sphere.
    let r_mid = a * (r_air / a).powf(0.5);
    for ci in 0..mesh.n_cells() {
        if mesh.cell_labels[ci] == 2 {
            let v = mesh.cell(ci);
            let rc = v
                .iter()
                .map(|&q| {
                    let pt = mesh.point(q as usize);
                    (pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt()
                })
                .sum::<f64>()
                / 4.0;
            if rc > r_mid {
                mesh.cell_labels[ci] = 4;
            }
        }
    }
    let eps = c(2.25, 0.0);
    let freq = ka * C0 / (2.0 * PI * a);
    let mut p = Problem::new(driven_study(freq), 2);
    p.materials.set(1, "eps_r", eps.re);
    p.materials.set(2, "eps_r", 1.0);
    p.materials.set(4, "eps_r", 1.0);
    p.materials.set(3, "eps_r", 1.0);
    p.materials.pml.insert(
        3,
        Pml {
            coord: PmlCoord::Spherical { center: [0.0; 3], r_in: r_air, thickness: r_pml - r_air },
            strength: 5.0,
            order: 2.0,
            re: 1.0,
            profile: femcore::props::PmlProfile::KappaPolynomial { kappa_max: 2.0 },
            ..Pml::default()
        },
    );
    p.excitations.push(ExcitationSpec {
        name: "mie".into(),
        excitation: Excitation::PlaneWave {
            k_hat: [0.0, 0.0, 1.0],
            pol: [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
            amplitude: 1.0,
            phase: 0.0,
            origin: [0.0; 3],
        },
        domains: vec![1],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: true,
    });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("full-sphere Mie solve");
    let q_vol = scalar(&res, "mie: sigma_sca (m^2)").unwrap() / (PI * a * a);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();
    let ff_auto = nyx_analysis::farfield_with_params(
        &exchange,
        &res,
        &serde_json::json!({
            "solution": 0,
            "field_prefix": "s",
            "incident_amplitude_v_per_m": 1.0,
        }),
    )
    .expect("NFFFT on scattered Mie fields (auto FieldFrame surface)");
    let ff_inner = nyx_analysis::farfield_with_params(
        &exchange,
        &res,
        &serde_json::json!({
            "solution": 0,
            "field_prefix": "s",
            "surface": { "kind": "domains", "inside": [1, 2] },
            "incident_amplitude_v_per_m": 1.0,
        }),
    )
    .expect("NFFFT on scattered Mie fields (inner sphere)");
    let (_, qs_exact) = mie_q(eps.sqrt().conj(), ka);
    let q_auto = ff_auto["scattering_cross_section_m2"].as_f64().unwrap() / (PI * a * a);
    let q_inner = ff_inner["scattering_cross_section_m2"].as_f64().unwrap() / (PI * a * a);
    eprintln!(
        "full-sphere Mie (ka={ka}): Q_vol = {q_vol:.6}, NFFFT inner = {q_inner:.6}, NFFFT outer = {q_auto:.6}, exact = {qs_exact:.6}"
    );
    assert!(((q_vol - qs_exact) / qs_exact).abs() < 5e-3, "volume Q_sca {q_vol} vs exact {qs_exact}");
    // One-sided P1 `H_s = i/(w mu0) curl E_s` trace on an n_ang=5 (18 deg) sphere converges at O(h^2) ~ 1.2%.
    assert!(((q_inner - qs_exact) / qs_exact).abs() < 1.5e-2, "inner NFFFT Q_sca {q_inner} vs exact {qs_exact}");
    assert!(((q_auto - qs_exact) / qs_exact).abs() < 2.5e-2, "outer NFFFT Q_sca {q_auto} vs exact {qs_exact}");
}

#[test]
fn scattered_field_solve_is_deterministic_across_thread_counts() {
    if let Ok(out_path) = std::env::var("NYX_FR36_DET_OUT") {
        let words: Vec<u8> = determinism_case().iter().flat_map(|x| x.to_bits().to_le_bytes()).collect();
        std::fs::write(&out_path, words).unwrap();
        return;
    }
    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-fr36-det-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args(["--exact", "scattered_field_solve_is_deterministic_across_thread_counts", "--nocapture"])
            .env("NYX_THREADS", nt)
            .env_remove("NYX_DETERMINISTIC_FACTOR")
            .env("NYX_FR36_DET_OUT", &out_path)
            .output()
            .expect("spawn FR-36 determinism subprocess");
        assert!(out.status.success(), "child NYX_THREADS={nt} failed: {}", String::from_utf8_lossy(&out.stderr));
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        bytes
    };
    let (a1, a8) = (run_child("1"), run_child("8"));
    assert!(!a1.is_empty());
    assert!(a1 == a8, "scattered-field driven solve differs between NYX_THREADS=1 and 8");
}

fn determinism_case() -> Vec<f64> {
    let lam = 1.0e-6;
    let freq = C0 / lam;
    let (a, d) = (0.4 * lam, 0.3 * lam);
    let (z1, z2, h) = (0.5 * lam, 0.5 * lam + d, 1.3 * lam);
    let zs = z_grid(&[0.0, z1, z2, h], lam / 10.0);
    let mesh = slab_cell(a, &zs, z1, z2, 3);
    let mut p = periodic_problem(freq, 2.25, 0.3, 0.0);
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![5], name: "bot".into(), m_max: Some(1), n_max: Some(1) });
    p.diffraction_ports.push(DiffractionPortBc { labels: vec![6], name: "top".into(), m_max: Some(1), n_max: Some(1) });
    p.excitations.push(ExcitationSpec {
        name: "gb".into(),
        excitation: Excitation::GaussianBeam {
            focus: [0.2 * lam, 0.2 * lam, h],
            w0: 0.8 * lam,
            k_hat: [0.0, 0.0, -1.0],
            pol: [c(0.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)],
            amplitude: 1.0,
            phase: 0.0,
            model: BeamModel::ComplexSource,
        },
        domains: vec![2],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: true,
    });
    let res = nyx_ems::solve_driven(&mesh, &p).expect("determinism solve");
    let mut v: Vec<f64> = res
        .solutions
        .iter()
        .flat_map(|s| s.cell_fields.iter().flat_map(|f| f.re.iter().chain(&f.im).copied()))
        .collect();
    v.push(scalar(&res, "gb: P_ext (W)").unwrap());
    let data = DrivenPortData::of(&res).unwrap();
    v.extend(data.s_matrix[0].iter().flatten().flat_map(|z| [z[0], z[1]]));
    v
}
