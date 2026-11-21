// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of `[superconducting]` TOML parsing (including mK sub-gap loss extensions),
//! shipped superconducting materials, the self-consistent BCS gap solver, and the Zimmermann
//! arbitrary-purity BCS complex conductivity and thin-film surface impedance evaluator.

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
        "format = 2\nname = \"sc_test\"\nkind = \"superconductor\"\n[superconducting]\ntc = \"9200 mK\"\ngap_0 = \"1395 ueV\"\nlambda_L = \"0.0333 um\"\nxi_0 = \"39 nm\"\nmean_free_path = \"20 nm\"\nsigma_n = \"5.5e7 S/m\"\ndynes_gamma = \"1.2 ueV\"\nr_res = \"12 nOhm/sq\"\nx_qp_floor = 1e-6\nnominal_thickness = \"100 nm\"\nfilm_model = \"kerr\"\n",
    )
    .unwrap();
    let sc = m.superconducting.as_ref().unwrap();
    assert!((sc.tc - 9.2).abs() < 1e-12);
    assert!((sc.gap_0.unwrap() - 1.395e-3 * EV_TO_J).abs() < 1e-12 * EV_TO_J);
    assert!((sc.lambda_l.unwrap() - 33.3e-9).abs() < 1e-18);
    assert!((sc.xi_0.unwrap() - 39.0e-9).abs() < 1e-18);
    assert!((sc.mean_free_path.unwrap() - 20.0e-9).abs() < 1e-18);
    assert!((sc.sigma_n.unwrap() - 5.5e7).abs() < 1e-5);
    assert!((sc.dynes_gamma.unwrap() - 1.2e-6 * EV_TO_J).abs() < 1e-15 * EV_TO_J);
    assert!((sc.r_res.unwrap() - 12e-9).abs() < 1e-22);
    assert!((sc.x_qp_floor.unwrap() - 1e-6).abs() < 1e-18);
    assert!((sc.nominal_thickness.unwrap() - 100e-9).abs() < 1e-18);
    assert_eq!(sc.film_model, Some(FilmModel::Kerr));

    // Omitted gap_0 defaults to 1.7638 * K_B * tc when effective_gap_0() is queried
    let m_default_gap = MaterialSpec::from_toml(
        "format = 1\nname = \"sc_def\"\n[superconducting]\ntc = \"1.2 K\"\nsigma_n = \"1e7 S/m\"\n",
    )
    .unwrap();
    let sc_def = m_default_gap.superconducting.unwrap();
    assert!(sc_def.gap_0.is_none());
    assert!((sc_def.effective_gap_0() - 1.7638 * K_B * 1.2).abs() < 1e-15 * K_B);

    // Missing or non-positive tc is rejected with a clear message
    let e_missing = MaterialSpec::from_toml("format = 2\n[superconducting]\nsigma_n = \"1e7 S/m\"\n").unwrap_err().to_string();
    assert!(e_missing.contains("tc") && e_missing.contains("missing"), "{e_missing}");
    let e_neg = MaterialSpec::from_toml("format = 2\n[superconducting]\ntc = \"-1 K\"\n").unwrap_err().to_string();
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

    let mut prev = d0;
    for k in 1..20 {
        let t_ratio = 0.05 * k as f64;
        let temp = t_ratio * tc;
        let d = bcs_gap(tc, None, temp).unwrap();
        assert!(d > 0.0 && d <= prev, "non-monotonic at t={t_ratio}: {d} vs {prev}");
        prev = d;
    }

    let d_half = bcs_gap(tc, None, 0.5 * tc).unwrap() / d0;
    assert!((d_half - 0.9569).abs() < 5e-4, "Delta(0.5 T_c)/Delta(0) = {d_half}");

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
    let freq = 1.0e9;
    let omega = 2.0 * PI * freq;
    let d = 50e-9;

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
            Some(1e-12),
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
    let ell_clean = 1e5 * xi_0;
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
fn park_2026_niobium_film_models_and_subgap_mk_loss() {
    let nb = Superconducting {
        tc: 9.2,
        gap_0: Some(1.395e-3 * EV_TO_J),
        lambda_l: Some(33.3e-9),
        xi_0: Some(39.0e-9),
        sigma_n: Some(5.5e7),
        dynes_gamma: Some(0.5e-6 * EV_TO_J),
        r_res: Some(8e-9),
        x_qp_floor: Some(1e-7),
        nominal_thickness: Some(100e-9),
        film_model: Some(FilmModel::Kerr),
        ..Default::default()
    };
    let freq = 7.0e9;
    let omega = 2.0 * PI * freq;
    let d = 100e-9;
    let t_base = 10e-3;

    let kerr = bcs_surface_impedance(&nb, freq, t_base, d, FilmModel::Kerr).unwrap();
    let two_sided = bcs_surface_impedance(&nb, freq, t_base, d, FilmModel::TwoSided).unwrap();
    let thin = bcs_surface_impedance(&nb, freq, t_base, d, FilmModel::ThinLimit).unwrap();

    assert!(kerr.lambda_eff > 45e-9 && kerr.lambda_eff < 60e-9);
    assert!(kerr.r_s >= 8e-9, "subgap r_res floor guarantees finite R_s at 10 mK: {}", kerr.r_s);
    assert!(kerr.sigma.re > 0.0, "dynes + x_qp guarantee finite sigma_1 at 10 mK: {}", kerr.sigma.re);
    assert!(kerr.l_k_sq > two_sided.l_k_sq && two_sided.l_k_sq > thin.l_k_sq && thin.l_k_sq > 0.0);

    // Verify series-to-parallel equivalence across T / T_c
    for t_frac in [0.01, 0.2, 0.5, 0.75, 0.85] {
        let ev = bcs_surface_impedance(&nb, freq, t_frac * nb.tc, d, FilmModel::Kerr).unwrap();
        let y_series = Complex64::new(1.0, 0.0) / ev.z_s;
        let y_parallel = Complex64::new(ev.g_eff, -1.0 / (omega * ev.l_eff));
        assert!(
            (y_series - y_parallel).norm() <= 1e-12 * y_series.norm(),
            "t={t_frac}: series admittance {y_series:?} vs parallel {y_parallel:?}"
        );
    }
}

#[test]
fn all_shipped_superconducting_materials_evaluate() {
    let lib = library();
    let mut count = 0usize;
    for name in lib.names() {
        let m = lib.get(&name).unwrap();
        if let Some(sc) = m.superconducting.as_ref() {
            count += 1;
            let ev = bcs_surface_impedance(
                sc,
                5.0e9,
                20e-3,
                sc.nominal_thickness.unwrap_or(50e-9),
                sc.film_model.unwrap_or(FilmModel::Kerr),
            )
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(ev.l_k_sq > 0.0 && ev.l_k_sq.is_finite(), "{name}: l_k_sq={}", ev.l_k_sq);
            assert!(ev.lambda_eff > 0.0 && ev.lambda_eff.is_finite(), "{name}: lambda_eff={}", ev.lambda_eff);
            assert!(ev.r_s >= 0.0 && ev.r_s.is_finite(), "{name}: r_s={}", ev.r_s);
        }
    }
    assert!(count >= 5, "expected at least 5 superconducting materials, got {count}");
}
