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
        if t.conductivity.is_some() || t.specific_heat.is_some() || t.expansion.is_some() {
            check_sources(&format!("{name}: [thermal]"), &t.source);
        }
        let e = &m.electrical;
        if e.permittivity.is_some() || e.permeability.is_some() || e.conductivity.is_some() {
            check_sources(&format!("{name}: [electrical]"), &e.source);
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
    let sn = lib.get("silicon_nitride").unwrap().optical.as_ref().unwrap();
    let n = sn.nk(0, 1.55e-6, None).unwrap().0;
    assert!((n - 1.9963).abs() < 5e-4, "{n}");
    let ln = lib.get("lithium_niobate").unwrap().optical.as_ref().unwrap();
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
