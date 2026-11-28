// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The shipped material library: provenance of every table, the metals against the published Lorentz-Drude tables, the other published
//! dispersion values, and the Schott, Herzberger and Drude-Lorentz models (derivatives, unit invariance, free-electron limit).
#![allow(clippy::type_complexity)] // tuples of per-cell arrays read better inline than as one-use aliases

use femcore::material::*;
use std::path::Path;

fn library() -> Library {
    Library::load_dir(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials")).unwrap_or_else(|e| panic!("{}", e.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n")))
}

fn parse(body: &str) -> MaterialSpec {
    MaterialSpec::from_toml(&format!("format = 1\n{body}")).unwrap_or_else(|e| panic!("{e}"))
}

// ----- provenance -----

fn check_sources(what: &str, s: &[Source]) {
    assert!(!s.is_empty(), "{what}: no source");
    for x in s {
        assert!(!x.cite.trim().is_empty(), "{what}: empty citation");
        if !x.is_verified() {
            assert!(x.flags_uncertainty(), "{what}: a source without a verification date must say unverified or approximate: {}", x.cite);
        }
    }
}

#[test]
fn every_library_file_and_table_has_a_source() {
    let lib = library();
    assert!(lib.len() >= 10);
    for name in lib.names() {
        let m = lib.get(&name).unwrap();
        check_sources(&format!("{name}: top level"), &m.source);
        if let Some(e) = &m.elastic {
            check_sources(&format!("{name}: [elastic]"), &e.source);
        }
        let t = &m.thermal;
        if t.conductivity.is_some()
            || t.specific_heat.is_some()
            || t.expansion.is_some()
            || t.conductivity_t.is_some()
            || t.specific_heat_t.is_some()
            || t.expansion_t.is_some()
        {
            check_sources(&format!("{name}: [thermal]"), &t.source);
        }
        let e = &m.electrical;
        if e.permittivity.is_some() || e.permeability.is_some() || e.conductivity.is_some() {
            check_sources(&format!("{name}: [electrical]"), &e.source);
        }
        if let Some(dl) = &m.dielectric {
            check_sources(&format!("{name}: [dielectric]"), &dl.source);
        }
        if let Some(sc) = &m.superconducting {
            check_sources(&format!("{name}: [superconducting]"), &sc.source);
        }
        if let Some(pz) = &m.piezoelectric {
            check_sources(&format!("{name}: [piezoelectric]"), &pz.source);
        }
        if let Some(pe) = &m.photoelastic {
            check_sources(&format!("{name}: [photoelastic]"), &pe.source);
        }
        if let Some(eo) = &m.electro_optic {
            check_sources(&format!("{name}: [electro_optic]"), &eo.source);
        }
        if let Some(ve) = &m.viscoelastic {
            check_sources(&format!("{name}: [viscoelastic]"), &ve.source);
        }
        if let Some(s2d) = &m.sheet_2d {
            check_sources(&format!("{name}: [sheet_2d]"), &s2d.source);
        }
        if let Some(nl) = &m.nonlinear {
            if !nl.source.is_empty() {
                check_sources(&format!("{name}: [nonlinear]"), &nl.source);
            }
            if let Some(c2) = &nl.chi2 {
                let srcs = if c2.source.is_empty() { &nl.source } else { &c2.source };
                check_sources(&format!("{name}: [nonlinear.chi2]"), srcs);
            }
            if let Some(c3) = &nl.chi3 {
                let srcs = if c3.source.is_empty() { &nl.source } else { &c3.source };
                check_sources(&format!("{name}: [nonlinear.chi3]"), srcs);
            }
        }
        if let Some(ma) = &m.interfaces.ma {
            check_sources(&format!("{name}: [interfaces.ma]"), &ma.source);
        }
        if let Some(ms) = &m.interfaces.ms {
            check_sources(&format!("{name}: [interfaces.ms]"), &ms.source);
        }
        if let Some(sa) = &m.interfaces.sa {
            check_sources(&format!("{name}: [interfaces.sa]"), &sa.source);
        }
        if let Some(o) = &m.optical {
            // each axis has its own source or shares the one of [optical]
            for (i, ax) in o.axes.models().into_iter().enumerate() {
                if ax.source.is_empty() {
                    check_sources(&format!("{name}: [optical] (axis {i})"), &o.source);
                } else {
                    check_sources(&format!("{name}: [optical] axis {i}"), &ax.source);
                }
            }
        }
        if let Some(sc) = &m.superconducting {
            check_sources(&format!("{name}: [superconducting]"), &sc.source);
        }
    }
}

#[test]
fn sources_survive_writing_and_reading() {
    let lib = library();
    for name in lib.names() {
        let m = lib.get(&name).unwrap();
        let back = MaterialSpec::from_toml(&m.to_toml()).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(back.source, m.source, "{name}");
        assert_eq!(back.elastic.as_ref().map(|e| e.source.clone()), m.elastic.as_ref().map(|e| e.source.clone()), "{name}");
        assert_eq!(back.thermal.source, m.thermal.source, "{name}");
        assert_eq!(back.optical.as_ref().map(|o| o.source.clone()), m.optical.as_ref().map(|o| o.source.clone()), "{name}");
        assert_eq!(back.superconducting.as_ref().map(|s| s.source.clone()), m.superconducting.as_ref().map(|s| s.source.clone()), "{name}");
        assert_eq!(back.nonlinear.as_ref().map(|n| n.source.clone()), m.nonlinear.as_ref().map(|n| n.source.clone()), "{name}");
    }
}

#[test]
fn source_entries_are_validated() {
    let bad = |s: &str| MaterialSpec::from_toml(&format!("format = 1\nsource = {s}\n")).unwrap_err().to_string();
    assert!(bad("{ cite = \"x\", verified = \"04-11-2026\" }").contains("YYYY-MM-DD"));
    assert!(bad("{ cite = \"x\", url = \"ftp://x\" }").contains("http"));
    assert!(bad("{ cite = \"x\", doi = \"https://doi.org/10.1/x\" }").contains("start with 10."));
    assert!(bad("{ url = \"https://x\" }").contains("cite"));
    assert!(bad("{ cite = \"x\", year = 1999 }").contains("unknown entry 'year'"));
    assert!(bad("\"\"").contains("must not be empty"));
    let ok = MaterialSpec::from_toml("format = 1\nsource = [\"a\", { cite = \"b\", doi = \"10.1000/x\", verified = \"2026-11-04\" }]\n").unwrap();
    assert_eq!(ok.source.len(), 2);
    assert!(ok.source[1].is_verified() && !ok.source[0].is_verified());
}

#[test]
fn library_values_that_were_corrected_after_checking() {
    let lib = library();
    let si = lib.get("silicon").unwrap().elastic.as_ref().unwrap().c;
    assert_eq!((si[0][0], si[0][1], si[3][3]), (165.6e9, 63.9e9, 79.5e9), "Hopcroft 2010");
    let sap = lib.get("sapphire").unwrap().elastic.as_ref().unwrap().c;
    assert_eq!((sap[0][0], sap[0][1], sap[0][2], sap[0][3], sap[2][2], sap[3][3]), (496.8e9, 163.6e9, 110.9e9, -23.5e9, 498.1e9, 147.4e9), "Wachtman 1960");
    let silica = lib.get("silica").unwrap().legacy_scalars();
    assert!((silica["E"] - 72.7e9).abs() < 1e-6 * 72.7e9 && silica["nu"] == 0.16 && silica["rho"] == 2201.0 && silica["Cp"] == 770.0 && silica["k"] == 1.38, "Corning 7980: {silica:?}");
}

// ----- published dispersion values -----

/// Rows (wavelength um, n, k) of the Rakic Lorentz-Drude tables of the refractiveindex.info database (CC0), nearest rows to 0.5, 1, 1.55 and 3 um.
const RAKIC: &[(&str, [(f64, f64, f64); 4])] = &[
    ("gold", [(0.49712, 0.81015, 1.8730), (0.99662, 0.25469, 5.9622), (1.5424, 0.57185, 9.6168), (2.9938, 1.7542, 18.725)]),
    ("silver", [(0.50321, 0.13208, 2.7927), (1.0013, 0.21994, 6.4341), (1.5431, 0.39685, 10.126), (3.0107, 1.2442, 19.905)]),
    ("copper", [(0.50054, 0.84204, 2.3418), (1.0075, 0.25557, 6.5372), (1.5520, 0.47077, 10.467), (2.9979, 1.5778, 20.200)]),
    ("aluminium", [(0.49815, 0.75416, 5.8785), (1.0005, 1.4765, 9.2887), (1.5536, 1.5399, 15.368), (2.9935, 4.3326, 29.332)]),
];

#[test]
fn metals_reproduce_the_published_lorentz_drude_tables() {
    let lib = library();
    for (name, rows) in RAKIC {
        let o = lib.get(name).unwrap().optical.as_ref().unwrap();
        for (l, n, k) in rows {
            let (nn, kk, ..) = o.nk(0, l * 1e-6, None).unwrap();
            assert!((nn - n).abs() <= 5e-4 * n.max(0.1), "{name} n at {l} um: {nn} vs {n}");
            assert!((kk - k).abs() <= 5e-4 * k, "{name} k at {l} um: {kk} vs {k}");
        }
    }
}

#[test]
fn metals_refuse_wavelengths_outside_the_fit() {
    let lib = library();
    let au = lib.get("gold").unwrap().optical.as_ref().unwrap();
    assert!(au.nk(0, 0.2e-6, None).is_err() && au.nk(0, 10e-6, None).is_err());
    // the complex index follows n - jk, so the permittivity of a metal has a negative real part and a negative imaginary part
    let e = lib.get("gold").unwrap().optical_permittivity(1.0e-6, None, [0.0; 3]).unwrap();
    assert!(e[0][0].re < 0.0 && e[0][0].im < 0.0, "{:?}", e[0][0]);
}

#[test]
fn silicon_nitride_and_lithium_niobate() {
    let lib = library();
    let sn_spec = lib.get("silicon_nitride").unwrap();
    assert_eq!(sn_spec.density, Some(3100.0));
    let sn_el = sn_spec.elastic.as_ref().unwrap();
    assert_eq!(sn_el.detect_class(), Class::Isotropic);
    let sn_scalars = sn_spec.legacy_scalars();
    assert!((sn_scalars["E"] - 290.0e9).abs() < 1e-6 * 290.0e9 && (sn_scalars["nu"] - 0.27).abs() < 1e-12 && sn_scalars["rho"] == 3100.0);
    let sn = sn_spec.optical.as_ref().unwrap();
    let n = sn.nk(0, 1.55e-6, None).unwrap().0;
    assert!((n - 1.9963).abs() < 5e-4, "{n}");

    let ln_spec = lib.get("lithium_niobate").unwrap();
    assert_eq!(ln_spec.density, Some(4640.0));
    let ln_el = ln_spec.elastic.as_ref().unwrap();
    assert_eq!(ln_el.class, Some(Class::Trigonal));
    assert_eq!(ln_el.detect_class(), Class::Trigonal);
    assert_eq!(
        (ln_el.c[0][0], ln_el.c[0][1], ln_el.c[0][2], ln_el.c[0][3], ln_el.c[2][2], ln_el.c[3][3]),
        (203.0e9, 57.3e9, 75.2e9, 8.5e9, 242.4e9, 59.5e9),
        "Weis and Gaylord 1985"
    );
    let ln = ln_spec.optical.as_ref().unwrap();
    let (no, ne) = (ln.nk(0, 1.55e-6, None).unwrap().0, ln.nk(2, 1.55e-6, None).unwrap().0);
    assert!((no - 2.2111).abs() < 1e-3 && (ne - 2.1379).abs() < 1e-3, "{no} {ne}");
    assert!(no > ne, "lithium niobate is a negative uniaxial crystal");
}

// ----- the new models -----

#[test]
fn drude_matches_the_free_electron_formula() {
    let m = parse(
        "[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 5.0]\nmodel = \"drude\"\nfrequency_unit = \"eV\"\neps_inf = 2.0\nplasma_frequency = 9.0\ndamping = 0.05\n",
    );
    let o = m.optical.unwrap();
    let ev = 1.602_176_634e-19 / 1.054_571_817e-34;
    let (wp, g) = (9.0 * ev, 0.05 * ev);
    for l in [0.6e-6, 1.55e-6, 4.0e-6] {
        let w = 2.0 * std::f64::consts::PI * 299_792_458.0 / l;
        let j = num_complex::Complex64::new(0.0, 1.0);
        let eps = 2.0 - wp * wp / (num_complex::Complex64::new(w * w, 0.0) - j * (g * w));
        let n = eps.sqrt();
        let (nn, kk, ..) = o.nk(0, l, None).unwrap();
        assert!((nn - n.re).abs() < 1e-9 * n.norm() && (kk + n.im).abs() < 1e-9 * n.norm() && kk >= 0.0, "{l}");
        let got = o.permittivity_diag(l, None).unwrap()[0];
        assert!((got - eps).norm() < 1e-9 * eps.norm());
    }
}

#[test]
fn frequency_units_are_equivalent() {
    let ev = |u: &str, scale: f64| {
        parse(&format!(
            "[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 5.0]\nmodel = \"drude_lorentz\"\nfrequency_unit = \"{u}\"\nplasma_frequency = {}\ndamping = {}\noscillators = [[0.3, {}, {}]]\n",
            9.0 * scale,
            0.05 * scale,
            2.0 * scale,
            0.4 * scale
        ))
        .optical
        .unwrap()
        .nk(0, 1.0e-6, None)
        .unwrap()
    };
    let a = ev("eV", 1.0);
    // 1 eV in THz (frequency, not angular) is e / h / 1e12
    let thz = 1.602_176_634e-19 / 6.626_070_15e-34 / 1e12;
    let b = ev("THz", thz);
    let c = ev("rad/s", 1.602_176_634e-19 / 1.054_571_817e-34);
    assert!((a.0 - b.0).abs() < 1e-6 * a.0 && (a.1 - b.1).abs() < 1e-6 * a.1, "{a:?} {b:?}");
    assert!((a.0 - c.0).abs() < 1e-9 * a.0 && (a.1 - c.1).abs() < 1e-9 * a.1, "{a:?} {c:?}");
}

#[test]
fn herzberger_does_not_depend_on_the_wavelength_unit() {
    // coefficients written for um, and the same function rewritten for nm
    let um = parse("[optical]\nwavelength_unit = \"um\"\nrange = [0.4, 1.0]\nmodel = \"herzberger\"\ncoefficients = [1.4, 0.01, 0.0005, -0.0001, 0.00001, 0.0]\n").optical.unwrap();
    let nm = parse("[optical]\nwavelength_unit = \"nm\"\nrange = [400, 1000]\nmodel = \"herzberger\"\ncoefficients = [1.4, 1e4, 5e8, -1e-10, 1e-17, 0.0]\n").optical.unwrap();
    for l in [0.45e-6, 0.7e-6, 0.95e-6] {
        let (a, b) = (um.nk(0, l, None).unwrap(), nm.nk(0, l, None).unwrap());
        assert!((a.0 - b.0).abs() < 1e-12, "{l}: {} vs {}", a.0, b.0);
        assert!((a.2 - b.2).abs() < 1e-6 * a.2.abs(), "derivative per metre: {} vs {}", a.2, b.2);
    }
}

#[test]
fn derivatives_of_the_new_models_match_finite_differences() {
    let models = [
        "[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2.0]\nmodel = \"schott\"\ncoefficients = [2.27, -0.0100, 0.0109, 0.0002, -0.000003, 0.0000002]\n",
        "[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2.0]\nmodel = \"herzberger\"\ncoefficients = [1.45, 0.005, 0.0001, -0.00001, 0.0, 0.0]\n",
        "[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2.0]\nmodel = \"drude_lorentz\"\nfrequency_unit = \"eV\"\neps_inf = 1.5\nplasma_frequency = 9.0\ndamping = 0.05\ndrude_weight = 0.8\noscillators = [[0.1, 3.0, 0.5], [0.4, 6.0, 1.0]]\n",
    ];
    for body in models {
        let o = parse(body).optical.unwrap();
        let (l, h) = (1.3e-6, 1e-11);
        let (_, _, dn, dk) = o.nk(0, l, None).unwrap();
        let (p, m) = (o.nk(0, l + h, None).unwrap(), o.nk(0, l - h, None).unwrap());
        let (fdn, fdk) = ((p.0 - m.0) / (2.0 * h), (p.1 - m.1) / (2.0 * h));
        assert!((dn - fdn).abs() <= 1e-4 * fdn.abs().max(1e3), "{body}: dn {dn} vs {fdn}");
        assert!((dk - fdk).abs() <= 1e-4 * fdk.abs().max(1e3), "{body}: dk {dk} vs {fdk}");
    }
}

#[test]
fn new_models_need_their_range_and_units_and_say_so() {
    let err = |body: &str| MaterialSpec::from_toml(&format!("format = 1\n{body}")).unwrap_err().to_string();
    let e = err("[optical]\nwavelength_unit = \"um\"\nmodel = \"drude\"\nfrequency_unit = \"eV\"\nplasma_frequency = 9\ndamping = 0.05\n");
    assert!(e.contains("needs a validity range"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"drude\"\nplasma_frequency = 9\ndamping = 0.05\n");
    assert!(e.contains("frequency_unit") && e.contains("required"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"drude\"\nfrequency_unit = \"eV\"\nplasma_frequency = 9\n");
    assert!(e.contains("damping") && e.contains("missing"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"drude\"\nfrequency_unit = \"furlong\"\nplasma_frequency = 9\ndamping = 0.1\n");
    assert!(e.contains("not known"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"lorentz\"\nfrequency_unit = \"eV\"\nplasma_frequency = 9\ndamping = 0.1\noscillators = [[1, 2, 0.1]]\n");
    assert!(e.contains("Drude term"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"schott\"\ncoefficients = [1, 2, 3]\n");
    assert!(e.contains("exactly 6 coefficients"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 2]\nmodel = \"lorentz\"\nfrequency_unit = \"eV\"\nplasma_frequency = 9\noscillators = [[1, -2, 0.1]]\n");
    assert!(e.contains("resonance must be positive"), "{e}");
}

#[test]
fn every_library_material_has_a_comparative_description() {
    let lib = library();
    assert_eq!(lib.len(), 46);
    for name in lib.names() {
        let m = lib.get(&name).unwrap();
        assert_eq!(m.format, 2, "{name} should be migrated to format = 2");
        let d = m.description.as_deref().unwrap_or("").trim();
        assert!(d.len() >= 40, "{name} description is missing or too short: {d:?}");
    }
}

#[test]
fn author_process_and_cryogenic_variants_reflect_published_differences() {
    let lib = library();

    // Photonics author variants
    let si_li = lib.get("silicon").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap();
    let si_salz = lib.get("silicon_salzberg").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap();
    let si_palik_vis = lib.get("silicon_palik").unwrap().optical.as_ref().unwrap().nk(0, 0.40e-6, None).unwrap();
    let si_palik_ir = lib.get("silicon_palik").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap();
    assert!((si_li.0 - 3.4757).abs() < 1e-4 && (si_salz.0 - 3.4777).abs() < 5e-4, "{si_li:?} vs {si_salz:?}");
    assert!(si_palik_vis.1 > 0.3 && si_palik_ir.1 == 0.0, "{si_palik_vis:?} {si_palik_ir:?}");

    // Gold Rakić vs Johnson & Christy
    let au_rakic = lib.get("gold").unwrap().optical.as_ref().unwrap().nk(0, 0.80e-6, None).unwrap();
    let au_jc = lib.get("gold_johnson").unwrap().optical.as_ref().unwrap().nk(0, 0.80e-6, None).unwrap();
    assert!(au_jc.0 < au_rakic.0 && au_jc.1 > 4.5 && au_rakic.1 > 4.5, "{au_rakic:?} vs {au_jc:?}");

    // LPCVD vs PECVD SiNx
    let sin_lpcvd = lib.get("silicon_nitride").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0;
    let sin_pecvd = lib.get("silicon_nitride_pecvd").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0;
    let sin_philipp = lib.get("silicon_nitride_philipp").unwrap().optical.as_ref().unwrap().nk(0, 0.633e-6, None).unwrap().0;
    assert!(sin_pecvd > sin_lpcvd && (sin_philipp - 2.011).abs() < 5e-3, "{sin_lpcvd} {sin_pecvd} {sin_philipp}");

    // Congruent LiNbO3 vs MgO:CLN + trigonal mechanics
    let ln = lib.get("lithium_niobate").unwrap();
    let ln_mgo = lib.get("lithium_niobate_mgo").unwrap();
    assert_eq!(ln.elastic.as_ref().unwrap().detect_class(), Class::Trigonal);
    assert!(ln_mgo.optical.as_ref().unwrap().nk(2, 1.55e-6, None).unwrap().0 < ln.optical.as_ref().unwrap().nk(2, 1.55e-6, None).unwrap().0);

    // GaAs Skauli vs Papatryfonos bandgap absorption
    let gaas_ir = lib.get("gallium_arsenide").unwrap().optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0;
    let gaas_vis = lib.get("gallium_arsenide_papatryfonos").unwrap().optical.as_ref().unwrap().nk(0, 0.55e-6, None).unwrap();
    assert!((gaas_ir - 3.3702).abs() < 1e-3 && gaas_vis.1 > 0.25, "{gaas_ir} {gaas_vis:?}");

    // Single-crystal AlN vs sputtered PiezoMUMPs vs Sc0.30Al0.70N
    let aln = lib.get("aluminum_nitride").unwrap();
    let aln_p = lib.get("aluminum_nitride_piezomumps").unwrap();
    let scaln = lib.get("scaln_30").unwrap();
    assert!(aln.elastic.as_ref().unwrap().c[2][2] > scaln.elastic.as_ref().unwrap().c[2][2]);
    assert_eq!(aln_p.elastic.as_ref().unwrap().detect_class(), Class::Isotropic);
    assert!((aln_p.legacy_scalars()["E"] - 330e9).abs() < 1e3);

    // Cryogenic shifts: Si stiffens + lowers index, Sapphire stiffens, Silica anomalously softens, OFHC Cu jumps to >600 W/(m*K)
    let si_4k = lib.get("silicon_4k").unwrap();
    assert!(si_4k.elastic.as_ref().unwrap().c[0][0] > lib.get("silicon").unwrap().elastic.as_ref().unwrap().c[0][0]);
    assert!(si_4k.optical.as_ref().unwrap().nk(0, 1.55e-6, None).unwrap().0 < 3.455);
    assert!(lib.get("silica_4k").unwrap().elastic.as_ref().unwrap().c[0][0] < lib.get("silica").unwrap().elastic.as_ref().unwrap().c[0][0]);
    assert!(lib.get("sapphire_4k").unwrap().elastic.as_ref().unwrap().c[0][0] > lib.get("sapphire").unwrap().elastic.as_ref().unwrap().c[0][0]);
    assert!(lib.get("copper_4k").unwrap().legacy_scalars()["k"] > 600.0);
    assert!(lib.get("niobium_sc").unwrap().notes.as_deref().unwrap_or("").contains("Tc ="));
}

#[test]
fn piezo_inheritance_thermodynamic_consistency_and_piezomumps_d15() {
    let lib = library();

    // 1. lithium_niobate_4k overrides eps_clamped cleanly without inheriting room-temp eps_free
    let lno4k = lib.get("lithium_niobate_4k").unwrap();
    assert!(lno4k.dielectric.as_ref().unwrap().eps_free.is_none());
    let ec = lno4k.dielectric.as_ref().unwrap().eps_clamped.as_ref().unwrap().matrix();
    assert!((ec[0][0] - 44.3).abs() < 1e-12 && (ec[2][2] - 26.5).abs() < 1e-12);

    // 2. PiezoMUMPs d15 has canonical negative sign (-3.6 pC/N)
    let mumps = lib.get("aluminum_nitride_piezomumps").unwrap();
    let pz = mumps.piezoelectric.as_ref().unwrap();
    assert!((pz.matrix[0][4] + 3.6e-12).abs() < 1e-16, "expected d15 = -3.6e-12, got {}", pz.matrix[0][4]);

    // 3. Inconsistent eps_free vs eps_clamped + d e^T / eps_0 (> 5%) is rejected on load
    let bad = r#"
format = 2
name = "bad_piezo"
source = { cite = "synthetic test" }
[elastic]
class = "isotropic"
E = "300 GPa"
nu = 0.25
source = { cite = "synthetic test" }
[dielectric]
eps_clamped = 10.0
eps_free = 25.0
source = { cite = "synthetic test" }
[piezoelectric]
form = "stress_charge"
e33 = "1.5 C/m^2"
source = { cite = "synthetic test" }
"#;
    let err = MaterialSpec::from_toml(bad).unwrap_err();
    assert!(err.message.contains("thermodynamic"), "unexpected error: {err}");
}
