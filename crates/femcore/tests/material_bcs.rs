// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of `[superconducting]` TOML parsing, shipped superconducting materials,
//! the self-consistent BCS gap solver, and the Zimmermann arbitrary-purity BCS complex
//! conductivity and thin-film surface impedance evaluator.

use femcore::material::*;
use num_complex::Complex64;
use std::f64::consts::PI;
use std::path::Path;

fn library() -> Library {
    Library::load_dir(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials"))
        .unwrap_or_else(|e| panic!("{}", e.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n")))
}

#[test]
fn superconducting_units_and_validation() {
    let m = MaterialSpec::from_toml(
        "format = 1\nname = \"sc_test\"\n[superconducting]\ntc = \"9200 mK\"\ngap_0 = \"1395 ueV\"\nlambda_L = \"0.0333 um\"\nxi_0 = \"39 nm\"\nmean_free_path = \"20 nm\"\nsigma_n = \"5.5e7 S/m\"\n",
    )
    .unwrap();
    let sc = m.superconducting.as_ref().unwrap();
    assert!((sc.tc - 9.2).abs() < 1e-12);
    assert!((sc.gap_0.unwrap() - 1.395e-3 * EV_TO_J).abs() < 1e-12 * EV_TO_J);
    assert!((sc.lambda_l.unwrap() - 33.3e-9).abs() < 1e-18);
    assert!((sc.xi_0.unwrap() - 39.0e-9).abs() < 1e-18);
    assert!((sc.mean_free_path.unwrap() - 20.0e-9).abs() < 1e-18);
    assert!((sc.sigma_n.unwrap() - 5.5e7).abs() < 1e-5);

    // Omitted gap_0 defaults to 1.7638 * K_B * tc when effective_gap_0() is queried
    let m_default_gap = MaterialSpec::from_toml(
        "format = 1\nname = \"sc_def\"\n[superconducting]\ntc = \"1.2 K\"\nsigma_n = \"1e7 S/m\"\n",
    )
    .unwrap();
    let sc_def = m_default_gap.superconducting.unwrap();
    assert!(sc_def.gap_0.is_none());
    assert!((sc_def.effective_gap_0() - 1.7638 * K_B * 1.2).abs() < 1e-15 * K_B);

    // Missing or non-positive tc is rejected with a clear message
    let e_missing = MaterialSpec::from_toml("format = 1\n[superconducting]\nsigma_n = \"1e7 S/m\"\n").unwrap_err().to_string();
    assert!(e_missing.contains("tc") && e_missing.contains("missing"), "{e_missing}");
    let e_neg = MaterialSpec::from_toml("format = 1\n[superconducting]\ntc = \"-1 K\"\n").unwrap_err().to_string();
    assert!(e_neg.contains("positive"), "{e_neg}");
}

#[test]
fn bcs_gap_accuracy_and_limits() {
    let tc = 9.2;
    let d0 = 1.7638 * K_B * tc;
    assert_eq!(bcs_gap(tc, None, 0.0).unwrap(), d0);
    assert_eq!(bcs_gap(tc, None, 0.01 * tc).unwrap(), d0);
    assert_eq!(bcs_gap(tc, None, tc).unwrap(), 0.0);
    assert_eq!(bcs_gap(tc, None, 1.5 * tc).unwrap(), 0.0);

    // Monotonicity and self-consistency across T / T_c in (0.05, 0.99)
    let mut prev = d0;
    for k in 1..20 {
        let t_ratio = 0.05 * k as f64;
        let temp = t_ratio * tc;
        let d = bcs_gap(tc, None, temp).unwrap();
        assert!(d > 0.0 && d <= prev, "non-monotonic at t={t_ratio}: {d} vs {prev}");
        prev = d;
    }

    // At T / T_c = 0.5, standard Muhlschlegel table value is Delta(0.5 T_c) / Delta(0) = 0.9569
    let d_half = bcs_gap(tc, None, 0.5 * tc).unwrap() / d0;
    assert!((d_half - 0.9569).abs() < 5e-4, "Delta(0.5 T_c)/Delta(0) = {d_half}");

    // Near T_c (t = 0.9999), Ginzburg-Landau asymptotic form Delta(T)/Delta(0) -> 1.7367 * sqrt(1 - T/T_c)
    let t_near = 0.9999;
    let d_near = bcs_gap(tc, None, t_near * tc).unwrap() / d0;
    let gl_want = 1.7367 * (1.0 - t_near).sqrt();
    assert!((d_near - gl_want).abs() < 2e-3 * gl_want, "near Tc: {d_near} vs {gl_want}");
}

#[test]
fn zimmermann_agrees_with_mattis_bardeen_dirty_limit() {
    let tc = 9.2;
    let gap_0 = 1.395e-3 * EV_TO_J;
    let sigma_n = 1.0e6;
    let freq = 1.0e9; // hbar * omega << 2 * Delta_0
    let omega = 2.0 * PI * freq;
    let d = 50e-9;

    // 1. T -> 0 analytic Mattis-Bardeen limit: sigma_2(0) / sigma_n = pi * Delta_0 / (hbar * omega)
    let mb_zero = surface_impedance_from_params(
        tc,
        Some(gap_0),
        None,
        None,
        None,
        Some(sigma_n),
        freq,
        0.01,
        d,
        FilmModel::ThinLimit,
    )
    .unwrap();
    let sigma_2_mb = -mb_zero.sigma.im;
    let want_ratio = PI * gap_0 / (HBAR * omega);
    let got_ratio = sigma_2_mb / sigma_n;
    assert!(
        (got_ratio - want_ratio).abs() < 2e-4 * want_ratio,
        "MB dirty limit sigma_2/sigma_n: got {got_ratio}, want {want_ratio}"
    );
    assert!(mb_zero.sigma.re / sigma_2_mb < 1e-15, "zero sub-gap loss at T -> 0");

    // 2. Zimmermann with ell / xi_0 = 1e-5 must match the direct Mattis-Bardeen dirty limit across T / T_c
    for t_frac in [0.01, 0.3, 0.5, 0.7, 0.85] {
        let temp = t_frac * tc;
        let mb = surface_impedance_from_params(
            tc,
            Some(gap_0),
            None,
            None,
            None,
            Some(sigma_n),
            freq,
            temp,
            d,
            FilmModel::Kerr,
        )
        .unwrap();
        let zim_dirty = surface_impedance_from_params(
            tc,
            Some(gap_0),
            None,
            Some(100e-9),
            Some(1e-12), // ell / xi_0 = 1e-5
            Some(sigma_n),
            freq,
            temp,
            d,
            FilmModel::Kerr,
        )
        .unwrap();
        assert!(
            (mb.sigma - zim_dirty.sigma).norm() <= 1e-4 * mb.sigma.norm(),
            "t={t_frac}: MB {:?} vs Zimmermann dirty {:?}",
            mb.sigma,
            zim_dirty.sigma
        );
        // Check against analytic finite-T low-frequency dirty-limit superfluid density:
        // sigma_2(T) / sigma_n = (pi * Delta(T) / (hbar * omega)) * tanh(Delta(T) / (2 * k_B * T))
        let delta_t = mb.delta_j;
        let want_s2 = sigma_n * (PI * delta_t / (HBAR * omega)) * ((delta_t + 0.5 * HBAR * omega) / (2.0 * K_B * temp)).tanh();
        let got_s2 = -mb.sigma.im;
        assert!(
            (got_s2 - want_s2).abs() <= 1e-3 * want_s2,
            "t={t_frac}: sigma_2(T) got {got_s2}, want {want_s2}"
        );
    }
}

#[test]
fn zimmermann_agrees_with_clean_london_limit() {
    let tc = 9.2;
    let gap_0 = 1.395e-3 * EV_TO_J;
    let lambda_l = 33.3e-9;
    let xi_0 = 39.0e-9;
    let ell_clean = 1e5 * xi_0; // ell >> xi_0
    let freq = 5.0e9;
    let omega = 2.0 * PI * freq;

    let clean = surface_impedance_from_params(
        tc,
        Some(gap_0),
        Some(lambda_l),
        Some(xi_0),
        Some(ell_clean),
        None,
        freq,
        0.01,
        100e-9,
        FilmModel::Kerr,
    )
    .unwrap();

    // In the clean limit at T -> 0 and hbar * omega << 2 * Delta_0, lambda_eff -> lambda_L
    // and sigma_2 -> 1 / (mu_0 * omega * lambda_L^2).
    let want_sigma_2 = 1.0 / (MU0 * omega * lambda_l * lambda_l);
    let got_sigma_2 = -clean.sigma.im;
    assert!(
        (got_sigma_2 - want_sigma_2).abs() < 2e-4 * want_sigma_2,
        "clean limit sigma_2: got {got_sigma_2}, want {want_sigma_2}"
    );
    assert!(
        (clean.lambda_eff - lambda_l).abs() < 1e-4 * lambda_l,
        "clean limit lambda_eff: got {}, want {lambda_l}",
        clean.lambda_eff
    );
}

#[test]
fn park_2026_niobium_film_models_and_temperature_curve() {
    // Park et al. npj Quantum Inf. (2026) Fig. 1d-e Nb parameters:
    // tc = 9.2 K, gap_0 = 1.395 meV, lambda_L = 33.3 nm, xi_0 = 39.0 nm, sigma_n = 5.5e7 S/m
    let lib = library();
    let nb = lib.get("niobium").unwrap().superconducting.as_ref().unwrap();
    let freq = 7.0e9;
    let omega = 2.0 * PI * freq;
    let d = 100e-9;
    let t_base = 10e-3; // 10 mK

    let kerr = bcs_surface_impedance(nb, freq, t_base, d, FilmModel::Kerr).unwrap();
    let two_sided = bcs_surface_impedance(nb, freq, t_base, d, FilmModel::TwoSided).unwrap();
    let thin = bcs_surface_impedance(nb, freq, t_base, d, FilmModel::ThinLimit).unwrap();

    // Intermediate purity (ell ~ 19.9 nm, xi_0 = 39.0 nm) increases lambda_eff above lambda_L = 33.3 nm to ~51.5 nm
    assert!(
        kerr.lambda_eff > 45e-9 && kerr.lambda_eff < 60e-9,
        "Nb effective penetration depth at 10 mK: {} m",
        kerr.lambda_eff
    );

    // At 10 mK, sigma_1 ~ 0, so L_k_sq matches the lossless formulas for each film model:
    let leff = kerr.lambda_eff;
    let want_thin = MU0 * leff * leff / d;
    let want_two = 0.5 * MU0 * leff / (0.5 * d / leff).tanh();
    let want_kerr = MU0 * leff / (d / leff).tanh();
    assert!((thin.l_k_sq - want_thin).abs() < 1e-12 * want_thin);
    assert!((two_sided.l_k_sq - want_two).abs() < 1e-12 * want_two);
    assert!((kerr.l_k_sq - want_kerr).abs() < 1e-12 * want_kerr);
    assert!(
        kerr.l_k_sq > two_sided.l_k_sq && two_sided.l_k_sq > thin.l_k_sq && thin.l_k_sq > 0.0,
        "Kerr ({}) > TwoSided ({}) > ThinLimit ({})",
        kerr.l_k_sq,
        two_sided.l_k_sq,
        thin.l_k_sq
    );

    // In the ultra-thin limit d << lambda_eff, all three models converge
    let d_thin = 0.5e-9;
    let k_ut = bcs_surface_impedance(nb, freq, t_base, d_thin, FilmModel::Kerr).unwrap();
    let ts_ut = bcs_surface_impedance(nb, freq, t_base, d_thin, FilmModel::TwoSided).unwrap();
    let tl_ut = bcs_surface_impedance(nb, freq, t_base, d_thin, FilmModel::ThinLimit).unwrap();
    assert!((k_ut.l_k_sq - tl_ut.l_k_sq).abs() < 1e-4 * tl_ut.l_k_sq);
    assert!((ts_ut.l_k_sq - tl_ut.l_k_sq).abs() < 1e-4 * tl_ut.l_k_sq);

    // Temperature curve over T / T_c in [0.01, 0.85]:
    // kinetic inductance L_k_sq(T) and surface resistance R_s(T) increase monotonically with T,
    // and the series-to-parallel decomposition 1 / Z_s == G_eff + 1 / (i omega L_eff) is exact.
    let mut prev_lk = kerr.l_k_sq;
    let mut prev_rs = kerr.r_s;
    for t_frac in [0.05, 0.2, 0.4, 0.6, 0.75, 0.85] {
        let ev = bcs_surface_impedance(nb, freq, t_frac * nb.tc, d, FilmModel::Kerr).unwrap();
        assert!(ev.l_k_sq >= prev_lk, "L_k_sq must increase with T at t={t_frac}");
        assert!(ev.r_s >= prev_rs, "R_s must increase with T at t={t_frac}");
        prev_lk = ev.l_k_sq;
        prev_rs = ev.r_s;

        // Verify series-to-parallel equivalence: 1 / (R_s + i omega L_k_sq) == G_eff + 1 / (i omega L_eff)
        let y_series = Complex64::new(1.0, 0.0) / ev.z_s;
        let y_parallel = Complex64::new(ev.g_eff, -1.0 / (omega * ev.l_eff));
        assert!(
            (y_series - y_parallel).norm() <= 1e-12 * y_series.norm(),
            "t={t_frac}: series admittance {y_series:?} vs parallel {y_parallel:?}"
        );
    }
    // At T = 0.85 T_c, fractional kinetic inductance shift (L_k(0.85 Tc) - L_k(0.01 Tc)) / L_k(0.01 Tc) > 40%
    let frac_shift = (prev_lk - kerr.l_k_sq) / kerr.l_k_sq;
    assert!(frac_shift > 0.4 && frac_shift < 2.0, "Nb fractional L_k shift at 0.85 T_c: {frac_shift}");
}

#[test]
fn all_shipped_superconducting_materials_evaluate() {
    let lib = library();
    for name in ["aluminium", "niobium", "tantalum", "nbn", "tin", "granular_aluminium"] {
        let m = lib.get(name).unwrap_or_else(|| panic!("missing {name}"));
        let sc = m.superconducting.as_ref().unwrap_or_else(|| panic!("{name} has no [superconducting]"));
        let ev = bcs_surface_impedance(sc, 5.0e9, 20e-3, 50e-9, FilmModel::Kerr)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(ev.l_k_sq > 0.0 && ev.l_k_sq.is_finite(), "{name}: l_k_sq={}", ev.l_k_sq);
        assert!(ev.lambda_eff > 0.0 && ev.lambda_eff.is_finite(), "{name}: lambda_eff={}", ev.lambda_eff);
        assert!(ev.r_s >= 0.0 && ev.r_s.is_finite(), "{name}: r_s={}", ev.r_s);
    }
    // Disordered superinductors (granular_aluminium, nbn) have much higher L_k_sq than elemental superconductors
    let lk_al = bcs_surface_impedance(
        lib.get("aluminium").unwrap().superconducting.as_ref().unwrap(),
        5.0e9,
        20e-3,
        20e-9,
        FilmModel::ThinLimit,
    )
    .unwrap()
    .l_k_sq;
    let lk_gral = bcs_surface_impedance(
        lib.get("granular_aluminium").unwrap().superconducting.as_ref().unwrap(),
        5.0e9,
        20e-3,
        20e-9,
        FilmModel::ThinLimit,
    )
    .unwrap()
    .l_k_sq;
    assert!(lk_gral > 100.0 * lk_al, "grAl L_k_sq ({lk_gral}) >> Al L_k_sq ({lk_al})");
}
