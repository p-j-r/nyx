//! Behaviour of `femcore::material`: the shipped library against the in-code data it replaces, published dispersion values, the three
//! ways to give a stiffness, round trips from file to matrix and back, and the plain-language error messages.

use femcore::material::*;
use std::path::Path;

fn library() -> Library {
    Library::load_dir(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials")).unwrap_or_else(|e| panic!("{}", e.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n")))
}

fn parse(body: &str) -> Result<MaterialSpec, MaterialError> {
    MaterialSpec::from_toml(&format!("format = 1\nname = \"t\"\n{body}"))
}

fn err(body: &str) -> String {
    match parse(body) {
        Ok(_) => panic!("expected an error for:\n{body}"),
        Err(e) => e.to_string(),
    }
}

fn close6(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6], tol: f64) {
    let s = a.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    for i in 0..6 {
        for j in 0..6 {
            assert!((a[i][j] - b[i][j]).abs() <= tol * s, "entry {i}{j}: {} vs {}", a[i][j], b[i][j]);
        }
    }
}

// ----- the library -----

#[test]
fn the_library_has_the_four_materials() {
    let lib = library();
    assert_eq!(lib.names(), vec!["air", "sapphire", "silica", "silicon"]);
    assert!(lib.get("Silicon").is_some(), "lookup ignores case");
    let e = lib.require("gold").unwrap_err().to_string();
    assert!(e.contains("available: air, sapphire, silica, silicon"), "{e}");
}

#[test]
fn library_stiffness_equals_the_in_code_data_it_replaces() {
    let lib = library();
    let (st, eps, k, alpha, rho) = sapphire();
    let s = lib.get("sapphire").unwrap();
    for euler in [[0.0; 3], [17.0, 33.0, 71.0]] {
        close6(&s.elastic.as_ref().unwrap().global(euler), &stiffness_global(&st, euler), 1e-12);
        close6(&lib.get("silicon").unwrap().elastic.as_ref().unwrap().global(euler), &stiffness_global(&silicon(), euler), 1e-12);
    }
    assert_eq!(s.density, Some(rho));
    let a = s.anisotropy([0.0; 3]);
    for (got, want) in [(a.permittivity.unwrap(), eps), (a.conductivity.unwrap(), k), (a.expansion.unwrap(), alpha)] {
        assert_eq!(got.matrix(), want.matrix());
    }
}

#[test]
fn legacy_scalars_of_silica() {
    let m = library().get("silica").unwrap().legacy_scalars();
    assert!((m["E"] - 73.1e9).abs() < 1e-3 * 73.1e9 * 1e-6 + 1.0);
    assert!((m["nu"] - 0.17).abs() < 1e-12);
    assert_eq!(m["rho"], 2203.0);
    assert_eq!(m["eps_r"], 3.75);
}

// ----- published optical values -----

#[test]
fn silica_matches_malitson() {
    let lib = library();
    let o = lib.get("silica").unwrap().optical.as_ref().unwrap();
    assert!((o.nk(0, 0.5875618e-6, None).unwrap().0 - 1.458464).abs() < 2e-6);
    assert!((o.nk(0, 1.55e-6, None).unwrap().0 - 1.44402).abs() < 1e-4);
}

#[test]
fn sapphire_ordinary_and_extraordinary() {
    let lib = library();
    let o = lib.get("sapphire").unwrap().optical.as_ref().unwrap();
    let (no, ne) = (o.nk(0, 0.5893e-6, None).unwrap().0, o.nk(2, 0.5893e-6, None).unwrap().0);
    assert!((no - 1.7682).abs() < 5e-4 && (ne - 1.7600).abs() < 5e-4, "{no} {ne}");
    assert!(no > ne, "sapphire is a negative uniaxial crystal");
    assert_eq!(o.nk(0, 1.55e-6, None).unwrap().0, o.nk(1, 1.55e-6, None).unwrap().0, "x and y are both ordinary");
    let (no, ne) = (o.nk(0, 1.55e-6, None).unwrap().0, o.nk(2, 1.55e-6, None).unwrap().0);
    assert!((no - 1.7462).abs() < 1e-3 && (ne - 1.7383).abs() < 1e-3, "{no} {ne}");
}

#[test]
fn silicon_table_reproduces_li_and_the_thermo_optic_shift() {
    let lib = library();
    let o = lib.get("silicon").unwrap().optical.as_ref().unwrap();
    assert_eq!(o.nk(0, 1.55e-6, None).unwrap().0, 3.4757);
    let warm = o.nk(0, 1.55e-6, Some(303.15)).unwrap().0;
    assert!((warm - (3.4757 + 1.8e-4 * 10.0)).abs() < 1e-12);
    let ng = o.group_index(0, 1.55e-6, None).unwrap();
    assert!(ng > 3.5 && ng < 3.7, "group index of silicon near 1.55 um is about 3.6, got {ng}");
}

#[test]
fn nothing_extrapolates_silently() {
    let lib = library();
    let si = lib.get("silicon").unwrap().optical.as_ref().unwrap();
    for l in [1.0e-6, 15e-6] {
        let e = si.nk(0, l, None).unwrap_err().to_string();
        assert!(e.contains("outside the range 1.2 to 14 um"), "{e}");
    }
    let silica = lib.get("silica").unwrap().optical.as_ref().unwrap();
    assert!(silica.nk(0, 0.2e-6, None).is_err() && silica.nk(0, 7e-6, None).is_err());
    let e = lib.get("silica").unwrap().optical_permittivity(8e-6, None, [0.0; 3]).unwrap_err().to_string();
    assert!(e.contains("material 'silica'") && e.contains("outside the range"), "{e}");
    assert!(lib.get("silicon").unwrap().optical_permittivity(1.55e-6, None, [0.0; 3]).is_ok());
}

#[test]
fn rotated_sapphire_permittivity_follows_the_optic_axis() {
    let lib = library();
    let s = lib.get("sapphire").unwrap();
    let l = 1.55e-6;
    let o = s.optical.as_ref().unwrap();
    let (no, ne) = (o.nk(0, l, None).unwrap().0, o.nk(2, l, None).unwrap().0);
    // crystal z (optic axis) lies along global -y after a 90 degree rotation about x
    let eps = s.optical_permittivity(l, None, [0.0, 90.0, 0.0]).unwrap();
    assert!((eps[1][1].re - ne * ne).abs() < 1e-12 && (eps[2][2].re - no * no).abs() < 1e-12 && (eps[0][0].re - no * no).abs() < 1e-12);
    assert!(eps[0][1].norm() < 1e-12 && eps[1][2].norm() < 1e-12);
    assert!(eps[0][0].im == 0.0, "lossless material has no imaginary part");
}

#[test]
fn loss_is_n_minus_j_k() {
    let m = parse("[optical]\nmodel = \"constant\"\nn = 2.0\nk = 0.01\n").unwrap();
    let n = m.optical.unwrap().index(0, 1.55e-6, None).unwrap();
    assert_eq!((n.re, n.im), (2.0, -0.01));
}

#[test]
fn cauchy_and_linear_table() {
    let m = parse("[optical]\nwavelength_unit = \"um\"\nrange = [0.4, 1.0]\nmodel = \"cauchy\"\ncoefficients = [1.5, 0.004]\n").unwrap();
    let o = m.optical.unwrap();
    let n = o.nk(0, 0.5e-6, None).unwrap().0;
    assert!((n - (1.5 + 0.004 / 0.25)).abs() < 1e-12);
    let t = parse("[optical]\nwavelength_unit = \"nm\"\nmodel = \"tabulated\"\ninterpolation = \"linear\"\ndata = [[500, 1.5, 0.0], [600, 1.4, 0.2]]\n").unwrap();
    let (n, k, dn, _) = t.optical.unwrap().nk(0, 550e-9, None).unwrap();
    assert!((n - 1.45).abs() < 1e-12 && (k - 0.1).abs() < 1e-12 && (dn + 0.1 / 100e-9).abs() < 1.0);
}

// ----- stiffness in every form -----

const SI_CUBIC: &str = "class = \"cubic\"\nunit = \"GPa\"\nc11 = 165.7\nc12 = 63.9\nc44 = 79.6\n";

fn elastic(body: &str) -> Elastic {
    parse(&format!("[elastic]\n{body}")).unwrap_or_else(|e| panic!("{e}")).elastic.unwrap()
}

#[test]
fn class_matrix_triangle_and_full_agree() {
    let by_class = elastic(SI_CUBIC);
    let tri = elastic(
        "unit = \"GPa\"\nmatrix = [\n [165.7, 63.9, 63.9, 0, 0, 0],\n [165.7, 63.9, 0, 0, 0],\n [165.7, 0, 0, 0],\n [79.6, 0, 0],\n [79.6, 0],\n [79.6],\n]\n",
    );
    let full = elastic(
        "unit = \"GPa\"\nmatrix = [\n [165.7, 63.9, 63.9, 0, 0, 0],\n [63.9, 165.7, 63.9, 0, 0, 0],\n [63.9, 63.9, 165.7, 0, 0, 0],\n [0, 0, 0, 79.6, 0, 0],\n [0, 0, 0, 0, 79.6, 0],\n [0, 0, 0, 0, 0, 79.6],\n]\n",
    );
    close6(&by_class.c, &tri.c, 1e-14);
    close6(&by_class.c, &full.c, 1e-14);
    assert_eq!(tri.detect_class(), Class::Cubic);
    assert_eq!(by_class.c[0][1], 63.9e9);
}

#[test]
fn compliance_matrix_and_constants() {
    let c = elastic(SI_CUBIC);
    let s = c.compliance();
    // cubic compliance constants in 1/TPa
    let t = |x: f64| x * 1e12;
    let by_const = elastic(&format!("class = \"cubic\"\ncompliance_unit = \"1/TPa\"\ns11 = {}\ns12 = {}\ns44 = {}\n", t(s[0][0]), t(s[0][1]), t(s[3][3])));
    close6(&by_const.c, &c.c, 1e-9);
    let rows: Vec<String> = s.iter().map(|r| format!("[{}]", r.iter().map(|x| format!("{:e}", x)).collect::<Vec<_>>().join(", "))).collect();
    let by_matrix = elastic(&format!("compliance_unit = \"1/Pa\"\ncompliance = [{}]\n", rows.join(",\n")));
    close6(&by_matrix.c, &c.c, 1e-9);
}

#[test]
fn trigonal_compliance_constants_carry_the_shear_factors() {
    // s56 = 2 s14 and s66 = 2 (s11 - s12) hold for the engineering compliance; the loader must derive them
    let c = library().get("sapphire").unwrap().elastic.clone().unwrap();
    let s = c.compliance();
    assert!((s[4][5] - 2.0 * s[0][3]).abs() < 1e-9 * s[0][0].abs() && (s[5][5] - 2.0 * (s[0][0] - s[0][1])).abs() < 1e-9 * s[0][0].abs());
    let body = format!(
        "class = \"trigonal\"\ncompliance_unit = \"1/Pa\"\ns11 = {:e}\ns12 = {:e}\ns13 = {:e}\ns14 = {:e}\ns33 = {:e}\ns44 = {:e}\n",
        s[0][0], s[0][1], s[0][2], s[0][3], s[2][2], s[3][3]
    );
    close6(&elastic(&body).c, &c.c, 1e-9);
}

#[test]
fn isotropic_constants_in_three_pairs() {
    let (e, nu) = (73.1e9, 0.17);
    let a = elastic("unit = \"GPa\"\nE = 73.1\nnu = 0.17\n");
    close6(&a.global([0.0; 3]), &isotropic_stiffness(e, nu), 1e-12);
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let b = elastic(&format!("lambda = \"{} Pa\"\nmu = \"{} Pa\"\n", lam, mu));
    close6(&a.c, &b.c, 1e-12);
    let kk = lam + 2.0 * mu / 3.0;
    let c = elastic(&format!("K = \"{kk} Pa\"\nG = \"{mu} Pa\"\n"));
    close6(&a.c, &c.c, 1e-12);
    assert_eq!(a.detect_class(), Class::Isotropic);
}

#[test]
fn monoclinic_and_triclinic_have_13_and_21_constants() {
    let names = Class::Monoclinic(Unique::Y).constants();
    assert_eq!(names.len(), 13);
    assert!(names.contains(&"c15") && !names.contains(&"c16"));
    assert!(Class::Monoclinic(Unique::Z).constants().contains(&"c16"));
    let mut body = String::from("class = \"monoclinic\"\nunit = \"GPa\"\n");
    for (n, v) in names.iter().zip([200.0, 60.0, 55.0, 5.0, 180.0, 50.0, 4.0, 160.0, 3.0, 40.0, 2.0, 45.0, 50.0]) {
        body.push_str(&format!("{n} = {v}\n"));
    }
    let m = elastic(&body);
    assert_eq!(m.c[0][4], m.c[4][0]);
    assert_eq!(m.c[0][4], 5.0e9);
}

// ----- round trips -----

fn approx_json(a: &serde_json::Value, b: &serde_json::Value, path: &str) {
    use serde_json::Value::*;
    match (a, b) {
        (Number(x), Number(y)) => {
            let (x, y) = (x.as_f64().unwrap(), y.as_f64().unwrap());
            assert!((x - y).abs() <= 1e-12 * x.abs().max(y.abs()), "{path}: {x} vs {y}");
        }
        (Array(x), Array(y)) => {
            assert_eq!(x.len(), y.len(), "{path}");
            for (i, (p, q)) in x.iter().zip(y).enumerate() {
                approx_json(p, q, &format!("{path}[{i}]"));
            }
        }
        (Object(x), Object(y)) => {
            assert_eq!(x.keys().collect::<std::collections::BTreeSet<_>>(), y.keys().collect::<std::collections::BTreeSet<_>>(), "{path}");
            for (k, v) in x {
                approx_json(v, &y[k], &format!("{path}.{k}"));
            }
        }
        _ => assert_eq!(a, b, "{path}"),
    }
}

#[test]
fn every_library_material_round_trips_file_to_matrix_and_back() {
    let lib = library();
    for name in lib.names() {
        let m = lib.get(&name).unwrap();
        let text = m.to_toml();
        let back = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
        approx_json(&serde_json::to_value(m).unwrap(), &serde_json::to_value(&back).unwrap(), &name);
        // and the written text reads to the same text again
        assert_eq!(back.to_toml(), MaterialSpec::from_toml(&back.to_toml()).unwrap().to_toml(), "{name}");
    }
}

#[test]
fn a_matrix_only_file_writes_back_as_a_matrix() {
    let m = parse("[elastic]\nunit = \"GPa\"\nmatrix = [[200,60,60,0,0,0],[200,60,0,0,0],[200,0,0,0],[80,0,0],[80,0],[80]]\n").unwrap();
    let text = m.to_toml();
    assert!(text.contains("matrix") && !text.contains("class"), "{text}");
    let back = MaterialSpec::from_toml(&text).unwrap();
    close6(&back.elastic.unwrap().c, &m.elastic.unwrap().c, 1e-14);
}

// ----- error messages -----

#[test]
fn missing_or_wrong_format_is_explained() {
    let e = MaterialSpec::from_toml("name = \"x\"\n").unwrap_err().to_string();
    assert!(e.contains("format = 1"), "{e}");
    let e = MaterialSpec::from_toml("format = 7\n").unwrap_err().to_string();
    assert!(e.contains("format 7") && e.contains("reads format 1"), "{e}");
    let e = MaterialSpec::from_toml("format = 1\nthis is not toml").unwrap_err().to_string();
    assert!(e.contains("not valid TOML"), "{e}");
}

#[test]
fn typos_are_reported_with_a_suggestion() {
    let e = err("densty = \"1 kg/m^3\"\n");
    assert!(e.contains("unknown key 'densty'") && e.contains("did you mean 'density'"), "{e}");
    let e = err("[elastic]\nclass = \"cubic\"\nunit = \"GPa\"\nc11 = 1\nc12 = 0.2\nc44 = 0.3\nc13 = 5\n");
    assert!(e.contains("c13") && e.contains("must equal c12"), "{e}");
}

#[test]
fn units_are_never_assumed() {
    let e = err("[elastic]\nclass = \"cubic\"\nc11 = 165.7\nc12 = 63.9\nc44 = 79.6\n");
    assert!(e.contains("[elastic.c11]") && e.contains("no unit") && e.contains("GPa"), "{e}");
    let e = err("density = \"3 GPa\"\n");
    assert!(e.contains("density") && e.contains("does not fit"), "{e}");
    let e = err("density = 2329\n");
    assert!(e.contains("no unit"), "{e}");
    assert!(parse("density = 2329\ndensity_unit = \"kg/m^3\"\n").is_ok(), "a unit may be declared next to the value");
    let e = err("[optical]\nmodel = \"sellmeier\"\nrange = [1, 2]\nB = [1]\nC = [0.1]\n");
    assert!(e.contains("wavelength_unit"), "{e}");
}

#[test]
fn class_needs_all_its_constants() {
    let e = err("[elastic]\nclass = \"cubic\"\nunit = \"GPa\"\nc11 = 165.7\nc12 = 63.9\n");
    assert!(e.contains("class 'cubic' needs the constants c11, c12, c44") && e.contains("missing: c44"), "{e}");
    let e = err("[elastic]\nunit = \"GPa\"\nc11 = 1\n");
    assert!(e.contains("need a class"), "{e}");
    let e = err("[elastic]\nunit = \"GPa\"\nclass = \"cubic\"\nmatrix = []\nc11 = 1\n");
    assert!(e.contains("one way only"), "{e}");
}

#[test]
fn asymmetric_matrix_is_named_entry_by_entry() {
    let e = err("[elastic]\nunit = \"GPa\"\nmatrix = [[200,60,60,0,0,0],[61,200,60,0,0,0],[60,60,200,0,0,0],[0,0,0,80,0,0],[0,0,0,0,80,0],[0,0,0,0,0,80]]\n");
    assert!(e.contains("not symmetric") && e.contains("c12 = 60 GPa") && e.contains("c21 = 61 GPa"), "{e}");
}

#[test]
fn unstable_stiffness_is_refused() {
    let e = err("[elastic]\nclass = \"cubic\"\nunit = \"GPa\"\nc11 = 10\nc12 = 60\nc44 = 5\n");
    assert!(e.contains("not positive definite") && e.contains("c11 > |c12|"), "{e}");
    let e = err("[elastic]\nunit = \"GPa\"\nE = 70\nnu = 0.6\n");
    assert!(e.contains("Poisson") && e.contains("-1 to 0.5"), "{e}");
    let e = err("[elastic]\nunit = \"GPa\"\nmatrix = [[200,60,60,0,0,0],[200,60,0,0,0],[200,0,0,0],[80,0,0],[80,0],[-80]]\n");
    assert!(e.contains("not positive definite"), "{e}");
}

#[test]
fn declared_class_must_match_the_matrix() {
    let e = err("[elastic]\nclass = \"cubic\"\nunit = \"GPa\"\nmatrix = [[165.7,63.9,63.9,0,0,0],[160,63.9,0,0,0],[165.7,0,0,0],[79.6,0,0],[79.6,0],[79.6]]\n");
    assert!(e.contains("symmetry of the declared class 'cubic'") && e.contains("c22"), "{e}");
    let e = err("[elastic]\nclass = \"hexagonal\"\nunit = \"GPa\"\nc11 = 160\nc12 = 60\nc13 = 50\nc33 = 180\nc44 = 40\nc66 = 70\n");
    assert!(e.contains("c66") && e.contains("(c11 - c12)/2") && e.contains("50 GPa"), "{e}");
}

#[test]
fn tensors_are_checked_in_plain_words() {
    let e = err("[electrical]\npermittivity = [2.0, -1.0, 3.0]\n");
    assert!(e.contains("relative permittivity must be positive in every direction"), "{e}");
    let e = err("[electrical]\npermittivity = [1.0, 2.0]\n");
    assert!(e.contains("2 entries") && e.contains("1 value"), "{e}");
    let e = err("[electrical]\npermittivity = [[1.0, 0.2, 0], [0.3, 1.0, 0], [0, 0, 1]]\n");
    assert!(e.contains("not symmetric"), "{e}");
    let e = err("[thermal]\nconductivity = [1.0, 2.0, 3.0]\n");
    assert!(e.contains("no unit"), "{e}");
    assert!(parse("[thermal]\nexpansion = [-1e-6, 2e-6, 2e-6]\nexpansion_unit = \"1/K\"\n").is_ok(), "negative expansion is physical");
    let e = err("[electrical]\npermeability = -1\n");
    assert!(e.contains("permeability must be positive"), "{e}");
}

#[test]
fn optical_models_are_validated() {
    let e = err("[optical]\nwavelength_unit = \"um\"\nmodel = \"sellmeier\"\nB = [1.0]\nC = [0.01]\n");
    assert!(e.contains("needs a validity range"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.1, 2.0]\nmodel = \"sellmeier\"\nB = [1.0]\nC = [0.25]\n");
    assert!(e.contains("resonance at 0.5"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 1.0]\nmodel = \"sellmeier\"\nB = [1.0, 2.0]\nC = [0.01]\n");
    assert!(e.contains("B has 2 entries but C has 1"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nmodel = \"tabulated\"\ndata = [[1.0, 1.5], [0.9, 1.4]]\n");
    assert!(e.contains("wavelengths must increase"), "{e}");
    let e = err("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 3.0]\nmodel = \"tabulated\"\ndata = [[1.0, 1.5], [2.0, 1.4]]\n");
    assert!(e.contains("wider than the table") && e.contains("extrapolate"), "{e}");
    assert!(parse("[optical]\nwavelength_unit = \"um\"\nrange = [0.5, 3.0]\nmodel = \"tabulated\"\nextrapolate = \"linear\"\ndata = [[1.0, 1.5], [2.0, 1.4]]\n").is_ok());
    let e = err("[optical]\nmodel = \"lorentz\"\n");
    assert!(e.contains("unknown model 'lorentz'"), "{e}");
    let e = err("[optical]\nmodel = \"constant\"\nn = 1.5\n[optical.ordinary]\nmodel = \"constant\"\nn = 1.5\n");
    assert!(e.contains("cannot be combined"), "{e}");
    let e = err("[optical.ordinary]\nmodel = \"constant\"\nn = 1.5\n");
    assert!(e.contains("both [optical.ordinary] and [optical.extraordinary]"), "{e}");
}

#[test]
fn library_reports_every_broken_file() {
    let dir = std::env::temp_dir().join(format!("nyx-materials-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("good.toml"), "format = 1\ndensity = \"1 kg/m^3\"\n").unwrap();
    std::fs::write(dir.join("bad1.toml"), "format = 1\ndensty = 1\n").unwrap();
    std::fs::write(dir.join("bad2.toml"), "name = \"x\"\n").unwrap();
    let errs = Library::load_dir(&dir).unwrap_err();
    std::fs::remove_dir_all(&dir).ok();
    assert_eq!(errs.len(), 2);
    assert!(errs.iter().any(|e| e.material.contains("bad1")) && errs.iter().any(|e| e.material.contains("bad2")));
}
