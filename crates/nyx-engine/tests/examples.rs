// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The example projects in `examples/` (written by `tests/make_examples.rs`): each file opens, passes the checks, meshes and solves through
//! engine commands, and the ones with a closed form give it. They replace the use-case tests: a device is a project file, not code.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples").join(format!("{name}.nyx.json"))
}

/// Opens an example, checks it, meshes it and solves its active physics (and what it depends on). Returns the engine and the reply.
fn solved(name: &str) -> (Engine, Value) {
    let mut e = Engine::new();
    let path = file(name);
    e.call("project_open", &json!({ "path": path.to_str().unwrap() })).unwrap_or_else(|err| panic!("{name}: {err}"));
    e.call("check_project", &json!({})).unwrap_or_else(|err| panic!("{name}: {err}"));
    let r = e.call("solve_run", &json!({})).unwrap_or_else(|err| panic!("{name}: {err}"));
    (e, r)
}

fn first_eigenvalue(e: &Engine, physics: &str) -> [f64; 2] {
    let r = e.result(physics).expect("a result").0.clone();
    r.solutions[0].eigenvalue.expect("an eigenvalue")
}

fn scalar_num(res: &nyx_physics::SolveResult, name: &str) -> Option<f64> {
    res.scalars.iter().find(|(n, _)| n == name).and_then(|(_, v)| match v {
        nyx_physics::Scalar::Num(x) => Some(*x),
        _ => None,
    })
}

fn probe_re(e: &mut Engine, physics: &str, expr: &str, points: Value) -> Vec<f64> {
    let out = e
        .call("probe", &json!({ "physics": physics, "expr": expr, "points": points }))
        .unwrap_or_else(|err| panic!("probe {physics} {expr}: {err}"));
    out["values"].as_array().unwrap().iter().map(|v| v.as_f64().expect("finite probed scalar")).collect()
}

#[test]
fn every_example_file_opens_and_passes_the_checks() {
    let mut n = 0;
    for entry in std::fs::read_dir(file("x").parent().unwrap()).unwrap() {
        let p = entry.unwrap().path();
        if !p.to_string_lossy().ends_with(".nyx.json") {
            continue;
        }
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": p.to_str().unwrap() }))
            .unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        e.call("check_project", &json!({})).unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        e.call("geometry_names", &json!({})).unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        n += 1;
    }
    assert!(n >= 17, "{n} examples");
}

#[test]
fn capacitor_has_a_potential_between_the_plates() {
    let (e, _) = solved("01_capacitor");
    let r = e.result("dfs").unwrap().0.clone();
    let v = r.solutions[0].field("V").expect("potential");
    let (lo, hi) = v.values.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| (l.min(*x), h.max(*x)));
    assert!((lo - 0.0).abs() < 1e-9 && (hi - 1.0).abs() < 1e-9, "{lo} {hi}");
}

#[test]
fn the_cantilever_bends_down_under_its_load() {
    let (e, _) = solved("02_cantilever");
    let r = e.result("sms").unwrap().0.clone();
    let u = r.solutions[0].field("u").expect("displacement");
    // the tip moves along -z by the Euler-Bernoulli value (F L^3 / 3 E I) to a few percent for this slender beam
    let uz_min = (0..u.values.len() / u.components.len())
        .map(|i| u.values[i * u.components.len() + u.components.iter().position(|c| c == "z").unwrap()])
        .fold(f64::INFINITY, f64::min);
    let (e_mod, l, h, q): (f64, f64, f64, f64) = (200e9, 1.0, 0.1, 1e5);
    let analytic = q * h * l.powi(3) / (3.0 * e_mod * h.powi(3) / 12.0);
    assert!(uz_min < 0.0 && (uz_min.abs() / analytic - 1.0).abs() < 0.03, "tip {uz_min} against {analytic}");
}

#[test]
fn a_heated_bar_expands_with_the_temperature_of_the_heat_physics() {
    let (mut e, r) = solved("04_thermal_bar");
    let physics: Vec<&str> = r["results"].as_array().unwrap().iter().map(|x| x["physics"].as_str().unwrap()).collect();
    assert_eq!(physics, ["dfs", "sms"], "the heat result feeds the thermal expansion: {r}");
    assert!(e.result("sms").is_some());

    // Exact linear steady profile T(x) = 293.15 + 100 * x
    let t_mid = probe_re(&mut e, "dfs", "T", json!([[0.25, 0.05], [0.5, 0.05], [0.75, 0.05]]));
    for (actual, want) in t_mid.iter().zip([318.15, 343.15, 368.15]) {
        assert!((actual - want).abs() < 1e-4, "bar T = {actual} vs exact {want}");
    }
    // Free axial thermal elongation: integral_0^L alpha*(T(x)-T_ref) dx = 0.5 * 23e-6 * 100 * 1.0 = 1.15e-3 m
    let ux_tip = probe_re(&mut e, "sms", "ux", json!([[1.0, 0.05]]))[0];
    let want_ux = 0.5 * 23e-6 * 100.0 * 1.0;
    assert!((ux_tip / want_ux - 1.0).abs() < 0.01, "bar tip ux = {ux_tip} vs {want_ux}");
    // Free thermal expansion carries near-zero axial stress away from clamped root x=0
    let sxx_far = probe_re(&mut e, "sms", "sxx", json!([[0.8, 0.05]]))[0];
    assert!(sxx_far.abs() < 1.0e5, "unconstrained thermal expansion residual sxx = {sxx_far}");
}

#[test]
fn the_wr90_guide_has_the_te10_cutoff() {
    let (e, _) = solved("06_wr90_waveguide");
    let n = first_eigenvalue(&e, "ems")[0];
    // fc = c / 2a = 6.557 GHz; n = sqrt(1 - (fc/f)^2) at 20 GHz
    let fc = 299792458.0 / (2.0 * 22.86e-3);
    let want = (1.0f64 - (fc / 20e9) * (fc / 20e9)).sqrt();
    assert!((n - want).abs() < 2e-3, "{n} against {want}");
}

#[test]
fn the_box_cavity_has_the_te101_frequency() {
    let (e, _) = solved("14_box_cavity");
    let f = first_eigenvalue(&e, "ems")[0];
    // 40 x 30 x 20 mm: f = c/2 sqrt((1/a)^2 + (1/d)^2) with a = 40 mm, d = 30 mm gives 6.2457 GHz
    let want = 299792458.0 / 2.0 * ((1.0 / 0.04f64).powi(2) + (1.0 / 0.03f64).powi(2)).sqrt();
    assert!((f / want - 1.0).abs() < 5e-3, "{f} against {want}");
}

#[test]
fn the_soi_wire_has_the_documented_te0_index_and_is_quasi_te() {
    let (e, _) = solved("15_soi_wire");
    let a = e.analyse("polarization", "ems").unwrap();
    let rows = a["rows"].as_array().unwrap();
    let te0 =
        rows.iter().find(|r| r["report"]["class"].as_str().unwrap_or("").contains("TE")).expect("a quasi-TE mode");
    let neff = te0["value"][0].as_f64().unwrap();
    // 500 x 300 nm wire in oxide on 3 um BOX: the repository documented 2.6718 with fixed indices 3.4757 and 1.444
    assert!((neff - 2.672).abs() < 0.02, "TE0 {neff}");
}

#[test]
fn heated_plate_matches_transient_convective_series_and_cooling_bounds() {
    let (mut e, _) = solved("03_heated_plate");
    let res = e.result("dfs").unwrap().0.clone();
    assert_eq!(res.solutions.len(), 21, "expected t=0 plus n_out=20 steps");
    let mean_t: Vec<f64> = res
        .solutions
        .iter()
        .map(|s| {
            let v = &s.field("T").unwrap().values;
            v.iter().sum::<f64>() / (v.len() as f64)
        })
        .collect();
    assert!(mean_t[0] >= 293.15 && mean_t[0] < 300.0, "initial plate mean T0 (incl. x=0 BC) = {}", mean_t[0]);
    for w in mean_t.windows(2) {
        assert!(w[1] >= w[0] - 1e-6, "transient heating mean T must increase monotonically: {w:?}");
    }
    let dt_first = mean_t[1] - mean_t[0];
    let dt_last = mean_t[20] - mean_t[19];
    assert!(dt_last < 0.25 * dt_first, "plate must approach thermal equilibrium: last={dt_last}, first={dt_first}");

    // Exact 2D separation-of-variables transient series at t = 20000 s (last solution index 20)
    let last_idx = res.solutions.len() - 1;
    let out = e
        .call(
            "probe",
            &json!({ "physics": "dfs", "expr": "T", "points": [[0.5, 0.1], [1.0, 0.1]], "solution": last_idx }),
        )
        .unwrap();
    let t_pts: Vec<f64> = out["values"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    assert!((t_pts[0] / 359.404614 - 1.0).abs() < 0.01, "T(0.5,0.1) = {} vs 359.4046", t_pts[0]);
    assert!((t_pts[1] / 321.798855 - 1.0).abs() < 0.01, "T(1.0,0.1) = {} vs 321.7989", t_pts[1]);
    assert!(500.0 > t_pts[0] && t_pts[0] > t_pts[1] && t_pts[1] > 293.15);
}

#[test]
fn block_cantilever_matches_3d_euler_bernoulli_tip_deflection() {
    let (mut e, _) = solved("05_block_cantilever");
    let ux = probe_re(&mut e, "sms", "ux", json!([[0.05, 1.0, 0.05]]))[0];
    let uz = probe_re(&mut e, "sms", "uz", json!([[0.05, 1.0, 0.05]]))[0];
    // F = 1e5 * 0.1 * 0.1 = 1000 N, I = 0.1^4 / 12, E = 200 GPa, L = 1.0 m -> EB tip uz = -2.0e-4 m
    let want_eb = -(1.0e5_f64 * 0.01) * 1.0_f64.powi(3) / (3.0 * 200.0e9 * (0.1_f64.powi(4) / 12.0));
    assert!(uz < 0.0 && (uz / want_eb - 1.0).abs() < 0.04, "3D tip uz = {uz} vs EB {want_eb}");
    assert!(ux.abs() < 1e-3 * uz.abs(), "lateral drift ux = {ux} should stay << tip deflection {uz}");
}

#[test]
fn loaded_waveguide_matches_exact_lse_root() {
    let (e, _) = solved("07_loaded_waveguide");
    let res = e.result("ems").unwrap().0.clone();
    let neff = res
        .solutions
        .iter()
        .filter_map(|s| s.eigenvalue)
        .filter(|[re, im]| *re > 1.0 && *re < 2.0 && im.abs() < 1e-6)
        .map(|[re, _]| re)
        .fold(f64::NEG_INFINITY, f64::max);
    // Exact half-filled WR-90 LSE_10 root at 8 GHz (a = 22.86 mm, eps_r = 4): 1.57589602
    assert!((neff - 1.57589602).abs() < 2.0e-3, "LSE10 neff = {neff} vs exact 1.57589602");
}

#[test]
fn sapphire_cantilever_matches_anisotropic_crystal_bending_modes() {
    let (e, _) = solved("08_sapphire_cantilever");
    let res = e.result("sms").unwrap().0.clone();
    let mut freqs: Vec<f64> = res
        .solutions
        .iter()
        .filter_map(|s| {
            let re = s.value.or_else(|| s.eigenvalue.map(|e| e[0]))?;
            let im = s.eigenvalue.map_or(0.0, |e| e[1]);
            assert!(im.abs() < 1e-3 * re.abs().max(1.0), "lossless sapphire must have Im(f) ~ 0: {re} + {im}i");
            (re > 1e3).then_some(re)
        })
        .collect();
    freqs.sort_by(f64::total_cmp);
    assert!(freqs.len() >= 2, "expected at least 2 bending modes: {freqs:?}");
    // Euler-deg [0, 90, 0] puts sapphire c-axis along y (E_c = 1/S33 = 460.8535 GPa, f1_EB = 173.828 kHz, f2_EB = 347.656 kHz)
    let f1_eb = 173_828.0_f64;
    let f2_eb = 347_656.0_f64;
    assert!((freqs[0] / f1_eb - 1.0).abs() < 0.06, "z-bending f1 = {} Hz vs EB {f1_eb}", freqs[0]);
    assert!((freqs[1] / f2_eb - 1.0).abs() < 0.12, "x-bending f2 = {} Hz vs EB {f2_eb}", freqs[1]);
    assert!(freqs[1] > 1.7 * freqs[0]);
}

#[test]
fn pml_slab_matches_exact_te0_and_tm0_guided_dispersion() {
    let (e, _) = solved("09_pml_slab");
    let res = e.result("ems").unwrap().0.clone();
    let guided: Vec<f64> = res
        .solutions
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            let [re, im] = s.eigenvalue?;
            let lf = scalar_num(&res, &format!("layer fraction {}", i + 1)).unwrap_or(1.0);
            (re > 1.0 && im.abs() < 1e-4 && lf < 0.05).then_some(re)
        })
        .collect();
    // Exact 10 mm eps_r=4 dielectric slab at 12 GHz with z-PEC height lz=10 mm (k_z = m*pi/lz):
    // TE_0 (m=0): 1.81307421, TE_1 (m=1): sqrt(TE0^2 - (c/(2*f*lz))^2), TM_1 (m=1): sqrt(TM0^2 - (c/(2*f*lz))^2)
    let kz_over_k0 = 299792458.0 / (2.0 * 12.0e9 * 10.0e-3);
    let te0_want = 1.81307421_f64;
    let te1_want = (te0_want * te0_want - kz_over_k0 * kz_over_k0).sqrt();
    let tm0_1d = 1.67671648_f64;
    let tm1_want = (tm0_1d * tm0_1d - kz_over_k0 * kz_over_k0).sqrt();

    for (label, want) in [("TE0", te0_want), ("TE1", te1_want), ("TM1", tm1_want)] {
        let got = guided
            .iter()
            .copied()
            .min_by(|a, b| (*a - want).abs().total_cmp(&(*b - want).abs()))
            .expect("guided slab mode");
        assert!((got - want).abs() < 2.0e-3, "slab {label} = {got} vs exact {want} (all={guided:?})");
    }
}

#[test]
fn plate_with_hole_satisfies_equilibrium_and_howland_stress_concentration() {
    let (mut e, _) = solved("10_plate_with_hole");
    // Divergence theorem: integral sigma_xx dA = L * H * sigma_0 = 0.4 * 0.2 * 1e6 = 80,000 N*m/m
    let int_sxx = e.call("integrate", &json!({ "physics": "sms", "expr": "sxx" })).unwrap()["value"].as_f64().unwrap();
    assert!((int_sxx / 80_000.0 - 1.0).abs() < 0.01, "2D equilibrium int(sxx) = {int_sxx} vs 80000");

    // Howland (1930) finite-width stress concentration (d/H = 0.40): K_{t,gross} * 1 MPa = 3.7261 MPa
    let max_sxx = e.call("extremum", &json!({ "physics": "sms", "expr": "sxx", "kind": "max" })).unwrap()["value"]
        .as_f64()
        .unwrap();
    assert!((max_sxx / 3.7261e6 - 1.0).abs() < 0.06, "2D peak sxx = {max_sxx} vs Howland 3.7261e6");
}

#[test]
fn anchor_loss_cantilever_has_passive_damping_and_confined_beam_energy() {
    let (e, _) = solved("11_anchor_loss");
    let res = e.result("sms").unwrap().0.clone();
    // Planar 2D PML eigenfrequency puts f_re, f_im, Q, layer fraction, P_pml in scalars
    let mut confined = Vec::new();
    for i in 0..res.solutions.len() {
        let f_re = scalar_num(&res, &format!("f re {} (Hz)", i + 1)).unwrap();
        let f_im = scalar_num(&res, &format!("f im {} (Hz)", i + 1)).unwrap();
        let lf = scalar_num(&res, &format!("layer fraction {}", i + 1)).unwrap();
        let q = scalar_num(&res, &format!("Q{}", i + 1)).unwrap();
        let p_pml = scalar_num(&res, &format!("P_pml {} (W)", i + 1)).unwrap();
        if lf < 0.05 && f_re > 50.0 {
            confined.push((i + 1, f_re, f_im, q, lf, p_pml));
        }
    }
    // Rigid-base plane-strain EB gives 189.8 Hz; compliant elastic substrate base softens fundamental bending into [110, 195] Hz
    let (_, f_re, f_im, q, _, p_pml) =
        confined.iter().copied().find(|m| (110.0..=195.0).contains(&m.1)).unwrap_or_else(|| {
            panic!("expected confined fundamental beam bending mode in [110, 195] Hz: {confined:?}")
        });
    assert!(f_im > 0.0 && q.is_finite() && q > 5.0 && p_pml > 0.0, "f={f_re}+{f_im}i, Q={q}, P_pml={p_pml}");
}

#[test]
fn drilled_plate_3d_satisfies_equilibrium_and_stress_concentration() {
    let (mut e, _) = solved("12_drilled_plate");
    // Divergence theorem over 0.4 x 0.2 x 0.04 m block: int sigma_xx dV = 0.4 * 0.2 * 0.04 * 1e6 = 3200 N*m
    let int_sxx = e.call("integrate", &json!({ "physics": "sms", "expr": "sxx" })).unwrap()["value"].as_f64().unwrap();
    assert!((int_sxx / 3_200.0 - 1.0).abs() < 0.02, "3D equilibrium int(sxx) = {int_sxx} vs 3200");

    let max_sxx = e.call("extremum", &json!({ "physics": "sms", "expr": "sxx", "kind": "max" })).unwrap()["value"]
        .as_f64()
        .unwrap();
    assert!(
        (max_sxx / 3.7261e6 - 1.0).abs() < 0.15 && max_sxx > 3.1e6,
        "3D hole peak sxx = {max_sxx} vs Howland 3.7261e6"
    );
}

#[test]
fn silicon_on_sapphire_solves_optical_te0_and_acoustic_eigenmodes() {
    let (mut e, _) = solved("13_silicon_on_sapphire");
    let ems = e.result("ems").unwrap().0.clone();
    let [neff_re, neff_im] = ems.solutions[0].eigenvalue.unwrap();
    // Dispersive silicon (n~3.476 at 1550 nm) + anisotropic r-plane sapphire raises TE0 slightly above fixed-n 2.3175 to ~2.346
    assert!((1.75..3.48).contains(&neff_re) && (neff_re - 2.346).abs() < 0.015, "SOS dispersive TE0 neff = {neff_re}");
    assert!(neff_im.abs() < 1e-4, "SOS TE0 must be lossless: {neff_im}");
    let pol = e.analyse("polarization", "ems").unwrap();
    assert!(pol["rows"][0]["report"]["te_fraction"].as_f64().unwrap() > 0.7);

    // Also solve secondary acoustic eigenfrequency study on r-plane sapphire + Si strip
    e.call("set_active", &json!({ "physics": "sms" })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let sms = e.result("sms").unwrap().0.clone();
    let confined = (0..sms.solutions.len()).any(|i| {
        let f_re = scalar_num(&sms, &format!("f re {} (Hz)", i + 1)).unwrap_or(0.0);
        let f_im = scalar_num(&sms, &format!("f im {} (Hz)", i + 1)).unwrap_or(-1.0);
        let lf = scalar_num(&sms, &format!("layer fraction {}", i + 1)).unwrap_or(1.0);
        (8.0e9..=11.0e9).contains(&f_re) && f_im >= -1e-6 * f_re && lf < 0.5
    });
    assert!(confined, "expected confined acoustic mode in silicon strip near 9.47 GHz");
}

#[test]
fn soi_hybrid_wire_exhibits_tilted_poincare_hybrid_polarization() {
    let (e, _) = solved("16_soi_hybrid_wire");
    let pol = e.analyse("polarization", "ems").unwrap();
    let rows = pol["rows"].as_array().unwrap();
    assert!(rows.len() >= 2);
    let n0 = rows[0]["value"][0].as_f64().unwrap();
    let n1 = rows[1]["value"][0].as_f64().unwrap();
    assert!((1.5..=2.3).contains(&n0) && (1.5..=2.3).contains(&n1), "hybrid bound n_eff = ({n0}, {n1})");
    let r0 = &rows[0]["report"];
    let r1 = &rows[1]["report"];
    let s0 = r0["s"].as_array().unwrap();
    let s1 = r1["s"].as_array().unwrap();
    let (s1_0, s2_0, s3_0) = (s0[0].as_f64().unwrap(), s0[1].as_f64().unwrap(), s0[2].as_f64().unwrap());
    let (s1_1, s2_1, s3_1) = (s1[0].as_f64().unwrap(), s1[1].as_f64().unwrap(), s1[2].as_f64().unwrap());
    assert!(
        s2_0.abs() > 0.30 || s2_1.abs() > 0.30,
        "sloped 75-deg wall must hybridize diagonal Stokes s2: {s2_0}, {s2_1}"
    );
    let p0 = r0["dop"].as_f64().unwrap();
    assert!(p0 > 0.7, "Stokes purity dop = {p0}");
    let dot = s1_0 * s1_1 + s2_0 * s2_1 + s3_0 * s3_1;
    assert!(dot < -0.5, "fundamental hybrid pair must be nearly antipodal on Poincare sphere: dot={dot}");
}

#[test]
fn wgm_cylindrical_pml_matches_exact_step_index_cylinder_modes() {
    let (e, _) = solved("18_wgm_cylindrical_pml");
    let res = e.result("ems").unwrap().0.clone();
    let exact_modes = [2.335860_f64, 2.376073, 2.506467];
    let guided: Vec<f64> = res
        .solutions
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            let [re, im] = s.eigenvalue?;
            let lf = scalar_num(&res, &format!("layer fraction {}", i + 1)).unwrap_or(1.0);
            (re > 1.0 && re < 3.4641 && im.abs() < 1e-4 && lf < 0.05).then_some(re)
        })
        .collect();
    assert!(!guided.is_empty(), "expected confined guided dielectric cylinder modes");
    let best_err =
        guided.iter().flat_map(|&g| exact_modes.iter().map(move |&w| (g - w).abs())).fold(f64::INFINITY, f64::min);
    assert!(
        best_err < 6.0e-3,
        "cylindrical guided modes {guided:?} missed analytical {exact_modes:?} (err={best_err})"
    );
}

#[test]
fn sphere_axisym_pml_matches_exact_mie_whispering_gallery_resonance() {
    // Exact spherical Bessel j_l(x), y_l(x) forward recurrence for l >= 0
    fn sph_bessel_jy(l_max: usize, x: f64) -> (Vec<f64>, Vec<f64>) {
        let mut j = vec![0.0; l_max + 1];
        let mut y = vec![0.0; l_max + 1];
        j[0] = x.sin() / x;
        y[0] = -x.cos() / x;
        if l_max >= 1 {
            j[1] = j[0] / x - x.cos() / x;
            y[1] = y[0] / x - x.sin() / x;
        }
        for l in 1..l_max {
            let c = (2 * l + 1) as f64 / x;
            j[l + 1] = c * j[l] - j[l - 1];
            y[l + 1] = c * y[l] - y[l - 1];
        }
        (j, y)
    }

    // Real-axis Mie characteristic function zeros closely bracket high-Q WGM Re(k_0 * R)
    let mie_zeros = |l: usize, n: f64, tm: bool| -> Vec<f64> {
        let f = |x: f64| -> f64 {
            let (j_in, _) = sph_bessel_jy(l + 1, n * x);
            let (j_out, y_out) = sph_bessel_jy(l + 1, x);
            // Riccati-Bessel derivatives: psi_l'(z) = (l+1)*j_l(z) - z*j_{l+1}(z)
            let psi_in = j_in[l];
            let dpsi_in = (l + 1) as f64 * j_in[l] - (n * x) * j_in[l + 1];
            let dchi_out = (l + 1) as f64 * y_out[l] - x * y_out[l + 1];
            // Wait: exterior outgoing hankel is ~ y_l dominant in WGM resonance (Im -> 0 limit uses Neumann y_l)
            let chi_out = y_out[l];
            let _ = j_out;
            if tm {
                (dpsi_in / (n * psi_in)) - (dchi_out / chi_out)
            } else {
                (n * dpsi_in / psi_in) - (dchi_out / chi_out)
            }
        };
        let mut roots = Vec::new();
        let dx = 0.005;
        let mut x = 1.5;
        let mut prev = f(x);
        while x < 5.0 {
            let next = f(x + dx);
            if prev.is_finite() && next.is_finite() && prev * next < 0.0 && (prev - next).abs() < 10.0 {
                let (mut a, mut b) = (x, x + dx);
                let sa = prev > 0.0;
                for _ in 0..50 {
                    let mid = 0.5 * (a + b);
                    if (f(mid) > 0.0) == sa {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                roots.push(0.5 * (a + b));
            }
            prev = next;
            x += dx;
        }
        roots
    };

    let mut mie_freqs = Vec::new();
    let r_sph = 1.5e-6;
    for l in 4..=7 {
        for tm in [false, true] {
            for xr in mie_zeros(l, 3.0, tm) {
                mie_freqs.push(299792458.0 * xr / (2.0 * std::f64::consts::PI * r_sph));
            }
        }
    }

    let (e, _) = solved("19_sphere_axisym_pml");
    let res = e.result("ems").unwrap().0.clone();
    // Pick the highest-Q confined whispering-gallery mode (minimum radiative loss Im(f) > 0)
    let best = res
        .solutions
        .iter()
        .filter(|s| {
            let [f_re, f_im] = s.eigenvalue.unwrap();
            f_im > 0.0 && f_re / (2.0 * f_im) > 100.0
        })
        .max_by(|a, b| {
            let [ar, ai] = a.eigenvalue.unwrap();
            let [br, bi] = b.eigenvalue.unwrap();
            (ar / (2.0 * ai)).total_cmp(&(br / (2.0 * bi)))
        })
        .expect("high-Q confined sphere WGM");
    let [f_re, f_im] = best.eigenvalue.unwrap();
    let q = f_re / (2.0 * f_im);
    assert!(f_im > 0.0 && q > 500.0, "passive radiative WGM requires high Q: {f_re} + {f_im}i, Q={q}");
    let closest_err = mie_freqs.iter().map(|&fm| (f_re - fm).abs() / fm).fold(f64::INFINITY, f64::min);
    assert!(closest_err < 0.025, "sphere WGM f_re = {f_re:e} missed Mie family (min rel err = {closest_err:e})");
}

#[test]
fn cavity_3d_spherical_pml_exhibits_passive_open_radiation_damping() {
    let (e, _) = solved("20_cavity_3d_spherical_pml");
    let res = e.result("ems").unwrap().0.clone();
    assert!(!res.solutions.is_empty());
    let mut passive_radiative = 0usize;
    for (i, s) in res.solutions.iter().enumerate() {
        let [f_re, f_im] = s.eigenvalue.expect("complex cavity eigenvalue");
        assert!(f_re > 1.0e9 && f_im > 0.0, "open spherical PML mode must be passive (Im f > 0): [{f_re}, {f_im}]");
        let q = f_re / (2.0 * f_im);
        let lf = scalar_num(&res, &format!("layer fraction {}", i + 1)).expect("layer fraction");
        let s_air = scalar_num(&res, &format!("share in air {}", i + 1)).expect("share in air");
        assert!((s_air - 1.0).abs() < 1e-6, "single non-PML region air share must equal 1.0: {s_air}");
        if lf > 0.02 && q < 200.0 {
            passive_radiative += 1;
        }
    }
    assert!(passive_radiative >= 1, "open air sphere in spherical PML must radiate strongly into PML");
}

#[test]
fn a_finer_mesh_confirms_or_refutes_a_number_through_convergence_check() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": file("06_wr90_waveguide").to_str().unwrap() })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let before = e.project().to_json();
    let r = e.call("convergence_check", &json!({ "physics": "ems", "factor": 0.6, "tolerance": 5e-3 })).unwrap();
    assert_eq!(r["converged"], true, "{r}");
    assert!(r["cells"]["refined"].as_u64().unwrap() > r["cells"]["coarse"].as_u64().unwrap(), "{r}");
    assert!(r["rows"].as_array().unwrap().iter().any(|x| x["quantity"] == "eigenvalue[0].re"), "{r}");
    assert_eq!(e.project().to_json(), before, "the project is not changed");
    // a mesh far too coarse is refuted
    e.call("set_mesh", &json!({ "h": 8e-3 })).unwrap();
    let r = e.call("convergence_check", &json!({ "physics": "ems", "factor": 0.5, "tolerance": 1e-6 })).unwrap();
    assert_eq!(r["converged"], false, "{r}");
    assert!(e.call("convergence_check", &json!({ "physics": "ems", "factor": 1.0 })).is_err());
    assert!(e.call("convergence_check", &json!({ "physics": "ems", "factor": 0.0 })).is_err());
    let rc = e.call("convergence_check", &json!({ "physics": "ems", "factor": 1.5 })).unwrap();
    assert!(rc["cells"]["refined"].as_u64().unwrap() > rc["cells"]["coarse"].as_u64().unwrap(), "{rc}");
}

#[test]
fn the_disk_on_a_post_has_a_contour_mode_with_a_quality_factor_from_anchor_loss() {
    let (e, r) = solved("17_disk_anchor_loss");
    assert!(r["results"].as_array().is_some_and(|a| !a.is_empty()), "{r}");
    let res = e.result("sms").expect("a result").0.clone();
    assert_eq!(res.study, "axisymmetric_eigenfrequency");
    let num = |name: &str| {
        res.scalars.iter().find(|(n, _)| n == name).and_then(|(_, v)| match v {
            nyx_physics::Scalar::Num(x) => Some(*x),
            _ => None,
        })
    };
    // every mode carries [Re f, Im f] and a layer fraction; the lossy ones a Q
    for (j, s) in res.solutions.iter().enumerate() {
        let [re, im] = s.eigenvalue.expect("an eigenvalue");
        let lf = num(&format!("layer fraction {}", j + 1)).expect("layer fraction");
        let q = num(&format!("Q{}", j + 1));
        println!("mode {}: f = {re:.5e} {im:+.3e}i Hz, layer fraction {lf:.3}, Q {q:?}", j + 1);
    }
    // the radial contour mode of the free disk (nu = 0.22): root lambda = 1.99765060 of lambda J0 / J1 = 1 - nu
    let c = (160e9_f64 / (2330.0 * (1.0 - 0.22 * 0.22))).sqrt();
    let want = 1.99765060 * c / (2.0 * std::f64::consts::PI * 6e-6);
    let best = res
        .solutions
        .iter()
        .enumerate()
        .filter(|(j, _)| num(&format!("layer fraction {}", j + 1)).unwrap() < 0.2)
        .min_by(|a, b| (a.1.eigenvalue.unwrap()[0] - want).abs().total_cmp(&(b.1.eigenvalue.unwrap()[0] - want).abs()))
        .expect("a mode that sits in the disk");
    let [f, f_im] = best.1.eigenvalue.unwrap();
    let q = num(&format!("Q{}", best.0 + 1)).expect("a quality factor (the anchor loses energy)");
    let radial = num(&format!("radial share {}", best.0 + 1)).expect("radial share");
    let disk = num(&format!("share in disk {}", best.0 + 1)).expect("share in disk");
    println!("contour mode: {f:.5e} {f_im:+.3e}i Hz against the free disk {want:.5e} Hz, Q = {q}");
    assert!((f / want - 1.0).abs() < 0.015, "{f} against {want}");
    assert!(f_im > 0.0 && q.is_finite() && q > 500.0, "Q = {q}, f_im = {f_im}");
    assert!(radial > 0.7 && disk > 0.8, "radial {radial}, disk {disk}");
}

#[test]
fn a_finer_mesh_is_compared_for_the_complex_eigenvalues_of_the_axisymmetric_study() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": file("17_disk_anchor_loss").to_str().unwrap() })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let before = e.project().to_json();
    let r = e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.7, "tolerance": 0.5 })).unwrap();
    println!("{}", serde_json::to_string_pretty(&r).unwrap());
    assert!(r["cells"]["refined"].as_u64().unwrap() > r["cells"]["coarse"].as_u64().unwrap(), "{r}");
    let quantities: Vec<&str> = r["rows"].as_array().unwrap().iter().filter_map(|x| x["quantity"].as_str()).collect();
    assert!(quantities.contains(&"eigenvalue[0].re") && quantities.contains(&"eigenvalue[0].im"), "{quantities:?}");
    assert_eq!(e.project().to_json(), before, "the project is not changed");
}

/// ADR 80: a mode of the axisymmetric study carries numeric frequency parts and energy shares, the shares per region under the region names.
#[test]
fn the_modes_of_the_disk_on_a_post_carry_what_identifies_them() {
    let (e, _) = solved("17_disk_anchor_loss");
    let res = e.result("sms").expect("a result").0.clone();
    let num = |name: &str| match res.scalars.iter().find(|(n, _)| n == name) {
        Some((_, nyx_physics::Scalar::Num(x))) => *x,
        other => panic!("scalar '{name}' is not a number: {other:?}"),
    };
    assert!(!res.scalars.iter().any(|(n, _)| n.contains("domain")), "the regions are named, not numbered");
    for (j, s) in res.solutions.iter().enumerate() {
        let (j1, [re, im]) = (j + 1, s.eigenvalue.unwrap());
        assert_eq!((num(&format!("f re {j1} (Hz)")), num(&format!("f im {j1} (Hz)"))), (re, im));
        let parts = num(&format!("radial share {j1}"))
            + num(&format!("axial share {j1}"))
            + num(&format!("azimuthal share {j1}"));
        let regions = num(&format!("share in disk {j1}"))
            + num(&format!("share in post {j1}"))
            + num(&format!("share in substrate {j1}"));
        assert!((parts - 1.0).abs() < 1e-9 && (regions - 1.0).abs() < 1e-9, "mode {j1}: {parts}, {regions}");
    }
    // the contour mode (a small layer fraction, nearest the free disk's 4.6e8 Hz) is radial and sits in the disk
    let best = (1..=res.solutions.len())
        .filter(|j| num(&format!("layer fraction {j}")) < 0.2 && num(&format!("azimuthal share {j}")) < 0.1)
        .min_by(|a, b| {
            (num(&format!("f re {a} (Hz)")) - 4.6e8).abs().total_cmp(&(num(&format!("f re {b} (Hz)")) - 4.6e8).abs())
        })
        .expect("a mode outside the layers");
    let (radial, disk) = (num(&format!("radial share {best}")), num(&format!("share in disk {best}")));
    println!("contour mode {best}: radial share {radial:.3}, share in disk {disk:.3}");
    assert!(radial > 0.5 && disk > 0.5, "radial {radial}, disk {disk}");
}

/// ADR 80: `convergence_check` pairs the modes of the two meshes by their eigenvalues, or follows one mode; by index nothing changes.
#[test]
fn a_convergence_check_follows_the_modes_when_their_order_changes() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": file("17_disk_anchor_loss").to_str().unwrap() })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let before = e.project().to_json();
    let args = |extra: Value| {
        let mut a = json!({ "physics": "sms", "factor": 0.7, "tolerance": 0.5 });
        a.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
        a
    };
    let by_index = e.call("convergence_check", &args(json!({}))).unwrap();
    assert!(by_index.get("pairs").is_none() && by_index.get("match").is_none(), "{by_index}");

    let near = e.call("convergence_check", &args(json!({ "match": "nearest" }))).unwrap();
    println!("{}", serde_json::to_string_pretty(&near).unwrap());
    assert_eq!(near["match"], "nearest");
    let pairs = near["pairs"].as_array().unwrap();
    assert!(!pairs.is_empty(), "{near}");
    let n_coarse = e.result("sms").unwrap().0.solutions.len();
    assert_eq!(
        pairs.len() + near["unpaired"]["coarse"].as_array().unwrap().len(),
        n_coarse,
        "every mode is paired or listed"
    );
    let worst = near["worst_relative"].as_f64().unwrap();
    assert!(worst.is_finite(), "{near}");
    for p in pairs {
        let i = p["coarse"].as_u64().unwrap();
        let name = format!("eigenvalue[{i}].re");
        let row = near["rows"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["quantity"] == name.as_str())
            .unwrap_or_else(|| panic!("{name}"));
        assert!(
            row["relative"].as_f64().unwrap().is_finite() && p["distance"].as_f64().unwrap().is_finite(),
            "{row} {p}"
        );
    }
    // the rows of a paired mode include its own scalars under its index on the present mesh
    let first = pairs[0]["coarse"].as_u64().unwrap() + 1;
    let quantity = format!("scalar layer fraction {first}");
    assert!(near["rows"].as_array().unwrap().iter().any(|r| r["quantity"] == quantity.as_str()), "{quantity}: {near}");

    let track = e.call("convergence_check", &args(json!({ "track": { "f": 4.5e8 } }))).unwrap();
    assert_eq!(track["match"], "track");
    assert_eq!(track["pairs"].as_array().unwrap().len(), 1, "{track}");
    let i = track["pairs"][0]["coarse"].as_u64().unwrap();
    let rows = track["rows"].as_array().unwrap();
    let re_row = rows.iter().find(|r| r["quantity"].as_str().unwrap().ends_with(".re")).unwrap();
    assert!(re_row["relative"].as_f64().unwrap() < 0.02, "tracked contour f_re drift >= 2%: {re_row}");
    assert!(rows.iter().all(|r| {
        let q = r["quantity"].as_str().unwrap();
        q.contains(&format!("[{i}]")) || q.starts_with("scalar ")
    }));
    assert_eq!(
        rows.iter().filter(|r| r["quantity"].as_str().unwrap().ends_with(".re")).count(),
        1,
        "one mode: {track}"
    );

    assert!(e.call("convergence_check", &args(json!({ "match": "closest" }))).is_err());
    assert!(e.call("convergence_check", &args(json!({ "track": { "g": 1.0 } }))).is_err());
    assert!(e.call("convergence_check", &args(json!({ "match": "index", "track": { "f": 4.6e8 } }))).is_err());
    assert_eq!(e.project().to_json(), before, "the project is not changed");
}

#[test]
fn pml_diagnostics_reports_mixed_elements_jacobian_s_dist_and_christoffel_check() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": file("13_silicon_on_sapphire").to_str().unwrap() })).unwrap();
    let d = e.call("pml_diagnostics", &json!({})).unwrap();
    assert!(d["element_counts"]["quad"].as_u64().unwrap_or(0) > 0, "expected structured Quad PML elements: {d}");
    assert!(d["min_pml_scaled_jacobian"].as_f64().unwrap_or(0.0) > 0.5, "expected positive PML scaled Jacobian: {d}");
    assert_eq!(d["s_dist_bounds"], json!([0.0, 1.0]), "{d}");
    assert!(d["min_projection"].as_f64().is_some(), "{d}");
    assert!(d["recommended_mpml_ratio"].as_f64().is_some(), "{d}");
    assert!(!d["domains"].as_array().unwrap().is_empty(), "{d}");
}

#[test]
fn add_pml_conditions_emits_infinite_element_for_dfs_and_static_sms() {
    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "ie_test", "dim": 2, "h": 0.2 })).unwrap();
    e.call("add_rect", &json!({ "name": "core", "min": [-1.0, -1.0], "max": [1.0, 1.0] })).unwrap();
    e.call("add_overlay", &json!({ "name": "regions", "regions": [{ "name": "diel", "input": "core" }] })).unwrap();
    e.call("add_absorbing_layer", &json!({ "name": "ie", "thickness": 0.5, "layers": 5.0 })).unwrap();
    e.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "diel*", "material": { "eps_r": 4.0 } })).unwrap();
    let r_dfs = e.call("add_pml_conditions", &json!({ "physics": "dfs", "auto_scale": true })).unwrap();
    assert_eq!(r_dfs["kind"], "infinite_element", "{r_dfs}");
    e.call("check_project", &json!({})).unwrap();

    e.call("add_physics", &json!({ "physics": "sms", "study": "static", "order": 2 })).unwrap();
    e.call("set_material", &json!({ "physics": "sms", "region": "diel*", "material": { "E": 100e9, "nu": 0.25 } }))
        .unwrap();
    let r_sms = e.call("add_pml_conditions", &json!({ "physics": "sms", "auto_scale": true })).unwrap();
    assert_eq!(r_sms["kind"], "infinite_element", "{r_sms}");
    e.call("check_project", &json!({})).unwrap();
}
