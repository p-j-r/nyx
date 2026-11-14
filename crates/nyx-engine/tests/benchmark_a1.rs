//! Benchmark A1 (Bindel and Govindjee 2005: a polysilicon disk on a post on a substrate, anchor loss of the second radial mode near
//! 715.6 MHz, Q = 6250) as a slow regression test: the model is built by engine commands only, the axisymmetric `sms` study solves it, and the
//! pass criteria are checked. A skeleton until the MCP run of the plan (step 8, `SMS_AXISYM_PLAN.md` section 9) has fixed the material and the
//! numbers: every test here is `#[ignore]`, and `CRITERIA_AGREED` is false, so the criteria are printed and not asserted.
//!
//!   cargo test --release -p nyx-engine --test benchmark_a1 -- --ignored --nocapture
//!
//! Free choices the paper does not give (to be fixed and stated by the A1 report): the substrate half width and depth, the thickness of the
//! layers, the mesh size. Here: about 4 and 4 shear wavelengths at 715 MHz and one wavelength of layer.

use nyx_engine::Engine;
use nyx_physics::Scalar;
use serde_json::{json, Value};

/// Set to true once paul has agreed the criteria below and the material is sourced (plan section 9, decisions 3).
const CRITERIA_AGREED: bool = false;
const F_PAPER: f64 = 715.6e6;
const Q_PAPER: f64 = 6250.0;
/// proposed by the A1 thread, not yet agreed: second radial mode within 2 % of the paper's frequency, Q within 20 %, 5 % under refinement
const TOL_F: f64 = 0.02;
const TOL_Q: f64 = 0.20;
const TOL_CONVERGED: f64 = 0.05;

/// The material is the real risk (plan section 9): the paper gives no E, nu, rho. These are the earlier thread's assumed values, to be replaced
/// by the HiQLab constants with file and line, or reported to paul. Do not tune them to hit 715.6 MHz.
const MATERIAL: (f64, f64, f64) = (160e9, 0.22, 2330.0);

struct Model {
    /// layer thickness, strength of the layers, mesh size of the disk and post, mesh size elsewhere (all SI, metres)
    layer: f64,
    strength: f64,
    h_disk: f64,
    h_bulk: f64,
    /// substrate half width and depth below the post
    half_width: f64,
    depth: f64,
}

impl Model {
    fn reference() -> Model {
        Model { layer: 7.4e-6, strength: 40.0, h_disk: 0.5e-6, h_bulk: 2.0e-6, half_width: 30e-6, depth: 30e-6 }
    }
}

/// Builds the project by commands (the same sequence as `examples/17_disk_anchor_loss.nyx.json`, with the paper's dimensions): disk r 0 to 10 um,
/// 2 um thick; post r 0 to 1 um, 0.5 um tall; a substrate below with layers on `+x` and `-z`; azimuthal order 0.
fn build(e: &mut Engine, m: &Model) {
    let um = 1e-6;
    let (post_h, disk_r, disk_t, post_r) = (0.5 * um, 10.0 * um, 2.0 * um, 1.0 * um);
    let mut c = |cmd: &str, args: Value| {
        e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd} {args}: {err}"));
    };
    c("project_new", json!({ "title": "A1 disk on a post", "dim": 2, "h": m.h_bulk }));
    c(
        "add_rect",
        json!({ "name": "substrate_box", "min": [0.0, -post_h - m.depth - m.layer], "max": [m.half_width + m.layer, -post_h] }),
    );
    c("add_rect", json!({ "name": "post_box", "min": [0.0, -post_h], "max": [post_r, 0.0] }));
    c("add_rect", json!({ "name": "disk_box", "min": [0.0, 0.0], "max": [disk_r, disk_t] }));
    c("add_rect", json!({ "name": "window", "min": [0.0, -post_h - m.depth], "max": [m.half_width, 2.0 * disk_t] }));
    c(
        "add_overlay",
        json!({ "name": "regions", "regions": [
        { "name": "substrate", "input": "substrate_box" }, { "name": "post", "input": "post_box" }, { "name": "disk", "input": "disk_box" } ] }),
    );
    c(
        "add_pml",
        json!({ "name": "pml", "window": "window", "sides": ["XMax", "ZMin"], "thickness": m.layer, "layers": 3.0, "grading": 2.0, "regions": ["substrate"] }),
    );
    c("add_size", json!({ "name": "size_disk", "sel": "region(disk) | region(post)", "h": m.h_disk }));
    c(
        "add_physics",
        json!({ "physics": "sms", "study": "axisymmetric_eigenfrequency", "order": 2,
        "study_params": { "n_modes": 12, "f_search": 715e6, "azimuthal_order": 0 } }),
    );
    let (ym, nu, rho) = MATERIAL;
    for r in ["substrate*", "post", "disk"] {
        c("set_material", json!({ "physics": "sms", "region": r, "material": { "E": ym, "nu": nu, "rho": rho } }));
    }
    c("add_pml_conditions", json!({ "physics": "sms", "strength": m.strength }));
    let fixed = json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 });
    c(
        "add_condition",
        json!({ "physics": "sms", "on": format!("outer & plane(x = {})", m.half_width + m.layer), "kind": "fixed", "params": fixed }),
    );
    c(
        "add_condition",
        json!({ "physics": "sms", "on": format!("outer & plane(z = {})", -post_h - m.depth - m.layer), "kind": "fixed", "params": fixed }),
    );
}

/// `(f, Q, layer fraction)` of the mode nearest `F_PAPER` among those that sit mostly in the resonator (layer fraction below 0.1): the way the
/// A1 report picks the mode.
fn picked_mode(e: &Engine) -> (f64, f64, f64) {
    let res = e.result("sms").expect("a result").0.clone();
    let num = |name: &str| {
        res.scalars.iter().find(|(n, _)| n == name).and_then(|(_, v)| match v {
            Scalar::Num(x) => Some(*x),
            _ => None,
        })
    };
    res.solutions
        .iter()
        .enumerate()
        .filter_map(|(j, s)| {
            let lf = num(&format!("layer fraction {}", j + 1))?;
            let q = num(&format!("Q{}", j + 1))?;
            (lf < 0.1).then(|| (s.eigenvalue.unwrap()[0], q, lf))
        })
        .min_by(|a, b| (a.0 - F_PAPER).abs().total_cmp(&(b.0 - F_PAPER).abs()))
        .expect("a mode in the resonator")
}

fn solve(m: &Model) -> (f64, f64, f64) {
    let mut e = Engine::new();
    build(&mut e, m);
    e.call("solve_run", &json!({})).expect("solve_run");
    picked_mode(&e)
}

fn check(what: &str, ok: bool) {
    println!("{} {what}", if ok { "PASS" } else { "MISS" });
    if CRITERIA_AGREED {
        assert!(ok, "{what}");
    }
}

#[test]
#[ignore = "benchmark A1: slow, and its criteria are not agreed yet"]
fn the_second_radial_mode_matches_the_paper_and_is_converged() {
    let base = Model::reference();
    let (f, q, lf) = solve(&base);
    println!("reference: f = {f:.6e} Hz, Q = {q:.1}, layer fraction {lf:.3} (paper {F_PAPER:.4e} Hz, Q {Q_PAPER})");
    check(&format!("f within {TOL_F} of the paper"), (f / F_PAPER - 1.0).abs() < TOL_F);
    check(&format!("Q within {TOL_Q} of the paper"), (q / Q_PAPER - 1.0).abs() < TOL_Q);
    // nyx's own criteria, which hold whatever the material: finer mesh, thicker layers
    let fine = Model { h_disk: 0.25e-6, h_bulk: 1.0e-6, ..Model::reference() };
    let (_, q_fine, _) = solve(&fine);
    check(
        &format!("Q changes by under {TOL_CONVERGED} between the two meshes"),
        (q_fine / q - 1.0).abs() < TOL_CONVERGED,
    );
    let thick = Model { layer: 2.0 * base.layer, ..Model::reference() };
    let (_, q_thick, _) = solve(&thick);
    check(
        &format!("Q changes by under {TOL_CONVERGED} when the layers are twice as thick"),
        (q_thick / q - 1.0).abs() < TOL_CONVERGED,
    );
}
