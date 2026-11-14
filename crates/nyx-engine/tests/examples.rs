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
    let (e, r) = solved("04_thermal_bar");
    let physics: Vec<&str> = r["results"].as_array().unwrap().iter().map(|x| x["physics"].as_str().unwrap()).collect();
    assert_eq!(physics, ["dfs", "sms"], "the heat result feeds the thermal expansion: {r}");
    assert!(e.result("sms").is_some());
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

macro_rules! solves {
    ($($fn_name:ident: $example:literal),* $(,)?) => {$(
        #[test]
        fn $fn_name() {
            let (e, r) = solved($example);
            assert!(r["results"].as_array().is_some_and(|a| !a.is_empty()), "{r}");
            let _ = e;
        }
    )*};
}

solves! {
    heated_plate_solves: "03_heated_plate",
    block_cantilever_solves: "05_block_cantilever",
    loaded_waveguide_solves: "07_loaded_waveguide",
    sapphire_cantilever_solves: "08_sapphire_cantilever",
    pml_slab_solves: "09_pml_slab",
    plate_with_hole_solves: "10_plate_with_hole",
    anchor_loss_solves: "11_anchor_loss",
    drilled_plate_solves: "12_drilled_plate",
    silicon_on_sapphire_solves: "13_silicon_on_sapphire",
    soi_hybrid_wire_solves: "16_soi_hybrid_wire",
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
    assert!(e.call("convergence_check", &json!({ "physics": "ems", "factor": 1.5 })).is_err());
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
    // the radial contour mode of the free disk, lambda = 2.05 of lambda J0 / J1 = 1 - nu: f = lambda c / (2 pi R), c = sqrt(E / (rho (1 - nu^2)))
    let c = (160e9_f64 / (2330.0 * (1.0 - 0.22 * 0.22))).sqrt();
    let want = 2.05 * c / (2.0 * std::f64::consts::PI * 6e-6);
    let best = res
        .solutions
        .iter()
        .enumerate()
        .filter(|(j, _)| num(&format!("layer fraction {}", j + 1)).unwrap() < 0.2)
        .min_by(|a, b| (a.1.eigenvalue.unwrap()[0] - want).abs().total_cmp(&(b.1.eigenvalue.unwrap()[0] - want).abs()))
        .expect("a mode that sits in the disk");
    let f = best.1.eigenvalue.unwrap()[0];
    let q = num(&format!("Q{}", best.0 + 1)).expect("a quality factor (the anchor loses energy)");
    println!("contour mode: {f:.5e} Hz against the free disk {want:.5e} Hz, Q = {q}");
    assert!((f / want - 1.0).abs() < 0.15, "{f} against {want}");
    assert!(q.is_finite() && q > 1.0, "Q = {q}");
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

    let track = e.call("convergence_check", &args(json!({ "track": { "f": 4.6e8 } }))).unwrap();
    assert_eq!(track["match"], "track");
    assert_eq!(track["pairs"].as_array().unwrap().len(), 1, "{track}");
    let i = track["pairs"][0]["coarse"].as_u64().unwrap();
    let rows = track["rows"].as_array().unwrap();
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
