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
        e.call("project_open", &json!({ "path": p.to_str().unwrap() })).unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        e.call("check_project", &json!({})).unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        e.call("geometry_names", &json!({})).unwrap_or_else(|err| panic!("{}: {err}", p.display()));
        n += 1;
    }
    assert!(n >= 16, "{n} examples");
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
    let uz_min = (0..u.values.len() / u.components.len()).map(|i| u.values[i * u.components.len() + u.components.iter().position(|c| c == "z").unwrap()]).fold(f64::INFINITY, f64::min);
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
    let te0 = rows.iter().find(|r| r["report"]["class"].as_str().unwrap_or("").contains("TE")).expect("a quasi-TE mode");
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
