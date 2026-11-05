//! Writes the example projects of the program to `examples/*.nyx.json` at the root of the repository, each built through engine commands
//! alone (the same commands a person or a host uses). Run it when an example changes:
//!
//!   cargo test --release -p nyx-engine --test make_examples -- --ignored
//!
//! The examples are data; `tests/examples.rs` loads every file, meshes it and solves the quick ones. The silicon-on-insulator wire takes its
//! geometry from `femcore::usecases::soi`, which phase 7 deletes: the file stays, this generator then builds that geometry from commands.

use nyx_engine::{Engine, Project};
use serde_json::{json, Value};
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples")
}

struct B {
    e: Engine,
}

impl B {
    fn new(title: &str, dim: u8, h: f64) -> B {
        let mut b = B { e: Engine::new() };
        b.c("project_new", json!({ "title": title, "dim": dim, "h": h }));
        b
    }

    fn c(&mut self, cmd: &str, args: Value) -> Value {
        self.e.call(cmd, &args).unwrap_or_else(|err| panic!("{} / {cmd} {args}: {err}", self.e.project().title))
    }

    fn rect(&mut self, name: &str, min: [f64; 2], max: [f64; 2]) {
        self.c("add_rect", json!({ "name": name, "min": min, "max": max }));
    }

    fn overlay(&mut self, regions: &[(&str, &str)]) {
        let r: Vec<Value> = regions.iter().map(|(n, i)| json!({ "name": n, "input": i })).collect();
        self.c("add_overlay", json!({ "name": "regions", "regions": r }));
    }

    /// A plate `lx` wide and `lz` high from the origin, one region `plate`.
    fn plate(&mut self, lx: f64, lz: f64) {
        self.rect("plate_box", [0.0, 0.0], [lx, lz]);
        self.overlay(&[("plate", "plate_box")]);
    }

    fn physics(&mut self, id: &str, study: &str, params: Value, order: u8) {
        self.c("add_physics", json!({ "physics": id, "study": study, "study_params": params, "order": order }));
    }

    fn material(&mut self, physics: &str, region: &str, m: Value) {
        self.c("set_material", json!({ "physics": physics, "region": region, "material": m }));
    }

    fn bc(&mut self, physics: &str, on: &str, kind: &str, params: Value) {
        self.c("add_condition", json!({ "physics": physics, "on": on, "kind": kind, "params": params }));
    }

    fn params(&mut self, physics: &str, params: Value) {
        self.c("set_physics_params", json!({ "physics": physics, "params": params }));
    }

    fn mesh(&mut self, v: Value) {
        self.c("set_mesh", v);
    }

    fn active(&mut self, physics: &str) {
        self.c("set_active", json!({ "physics": physics }));
    }

    fn save(self, file: &str) {
        let p = dir();
        std::fs::create_dir_all(&p).unwrap();
        let mut e = self.e;
        e.call("check_project", &json!({})).unwrap_or_else(|err| panic!("{file}: {err}"));
        std::fs::write(p.join(format!("{file}.nyx.json")), e.project().to_json() + "\n").unwrap();
    }
}

const FIXED: fn() -> Value = || json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 });

fn steel_like(b: &mut B, region: &str) {
    b.material("sms", region, json!({ "E": 200e9, "nu": 0.3, "rho": 7850.0, "alpha": 12e-6 }));
}

fn capacitor() -> B {
    let mut b = B::new("Parallel-plate capacitor with two dielectrics (electrostatics, 2D)", 2, 0.025);
    b.rect("left_box", [0.0, 0.0], [0.5, 0.5]);
    b.rect("right_box", [0.5, 0.0], [1.0, 0.5]);
    b.overlay(&[("left", "left_box"), ("right", "right_box")]);
    b.physics("dfs", "electrostatic", Value::Null, 2);
    b.material("dfs", "left", json!({ "eps_r": 2.0 }));
    b.material("dfs", "right", json!({ "eps_r": 8.0 }));
    b.bc("dfs", "outer & plane(x = 0)", "potential", json!({ "value": 1.0 }));
    b.bc("dfs", "outer & plane(x = 1)", "ground", Value::Null);
    b
}

fn cantilever() -> B {
    let mut b = B::new("Cantilever beam (solid mechanics, static, 2D plane stress)", 2, 0.02);
    b.plate(1.0, 0.1);
    b.physics("sms", "static", Value::Null, 2);
    b.params("sms", json!({ "plane": "stress" }));
    steel_like(&mut b, "plate");
    b.bc("sms", "outer & plane(x = 0)", "fixed", FIXED());
    b.bc("sms", "outer & plane(x = 1)", "load", json!({ "traction": [0.0, 0.0, -1e5] }));
    b
}

fn heated_plate() -> B {
    let mut b = B::new("Heated plate with convective cooling (heat transfer, 2D, transient)", 2, 0.025);
    b.plate(1.0, 0.2);
    b.physics("dfs", "heat_transient", json!({ "t_end": 20000.0, "n_out": 20, "t0": "293.15" }), 2);
    b.material("dfs", "plate", json!({ "k": 50.0, "rho": 7800.0, "Cp": 470.0 }));
    b.bc("dfs", "outer & plane(x = 0)", "temperature", json!({ "value": 500.0 }));
    b.bc("dfs", "outer - plane(x = 0)", "convection", json!({ "h": 25.0, "t_ext": 293.15 }));
    b
}

fn thermal_bar() -> B {
    let mut b = B::new("Thermal expansion of a bar (heat + solid, 2D)", 2, 0.02);
    b.plate(1.0, 0.1);
    b.physics("dfs", "heat_stationary", Value::Null, 2);
    b.material("dfs", "plate", json!({ "k": 160.0 }));
    b.bc("dfs", "outer & plane(x = 0)", "temperature", json!({ "value": 293.15 }));
    b.bc("dfs", "outer & plane(x = 1)", "temperature", json!({ "value": 393.15 }));
    b.physics("sms", "static", Value::Null, 2);
    b.params("sms", json!({ "plane": "stress", "thermal": "imported", "t_ref": 293.15 }));
    b.material("sms", "plate", json!({ "E": 70e9, "nu": 0.33, "alpha": 23e-6, "rho": 2700.0 }));
    b.bc("sms", "outer & plane(x = 0)", "fixed", FIXED());
    b.active("sms");
    b
}

fn block_cantilever() -> B {
    let mut b = B::new("Block cantilever (solid mechanics, 3D, a square section swept along y)", 2, 0.035);
    b.plate(0.1, 0.1);
    b.mesh(json!({ "sweep": [{ "thickness": 1.0, "n": 16 }] }));
    b.physics("sms", "static", Value::Null, 2);
    steel_like(&mut b, "plate");
    b.bc("sms", "ymin", "fixed", FIXED());
    b.bc("sms", "ymax", "load", json!({ "traction": [0.0, 0.0, -1e5] }));
    b
}

fn waveguide() -> B {
    let mut b = B::new("WR-90 rectangular waveguide modes (EM waves, 2D cross-section)", 2, 1.0e-3);
    b.plate(22.86e-3, 10.16e-3);
    b.physics("ems", "modes", json!({ "freq": 20e9, "n_target": 0.9, "n_modes": 4 }), 2);
    b.material("ems", "plate", json!({ "eps_r": 1.0, "mu_r": 1.0 }));
    b
}

fn loaded_waveguide() -> B {
    let (lx, lz) = (22.86e-3, 10.16e-3);
    let mut b = B::new("Dielectric-loaded waveguide (EM waves, 2D cross-section)", 2, 0.7e-3);
    b.rect("air_box", [0.0, 0.0], [lx, lz]);
    b.rect("dielectric_box", [0.0, 0.0], [11.43e-3, lz]);
    b.overlay(&[("air", "air_box"), ("dielectric", "dielectric_box")]);
    b.physics("ems", "modes", json!({ "freq": 8e9, "n_target": 1.7, "n_modes": 3 }), 2);
    b.material("ems", "air", json!({ "eps_r": 1.0, "mu_r": 1.0 }));
    b.material("ems", "dielectric", json!({ "eps_r": 4.0, "mu_r": 1.0 }));
    b
}

fn sapphire_cantilever() -> B {
    let mut b = B::new("Sapphire cantilever, anisotropic crystal (solid mechanics, 3D eigenfrequencies)", 2, 0.05e-3);
    b.plate(0.2e-3, 0.1e-3);
    b.mesh(json!({ "sweep": [{ "thickness": 1.0e-3, "n": 20 }] }));
    b.physics("sms", "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }), 2);
    b.material("sms", "plate", json!({ "material": "sapphire", "euler_deg": [0.0, 90.0, 0.0] }));
    b.bc("sms", "ymin", "fixed", FIXED());
    b
}

fn pml_slab() -> B {
    let (lx, lz) = (40e-3, 10e-3);
    let mut b = B::new("Dielectric slab guide with PML ends (EM waves, complex modes)", 2, 0.8e-3);
    b.rect("air_box", [0.0, 0.0], [lx, lz]);
    b.rect("slab_box", [15e-3, 0.0], [25e-3, lz]);
    b.rect("pml_left_box", [0.0, 0.0], [8e-3, lz]);
    b.rect("pml_right_box", [32e-3, 0.0], [lx, lz]);
    b.overlay(&[("air", "air_box"), ("slab", "slab_box"), ("pml_left", "pml_left_box"), ("pml_right", "pml_right_box")]);
    b.physics("ems", "modes", json!({ "freq": 12e9, "n_target": 1.9, "n_modes": 6 }), 2);
    for (n, eps) in [("air", 1.0), ("slab", 4.0), ("pml_left", 1.0), ("pml_right", 1.0)] {
        b.material("ems", n, json!({ "eps_r": eps, "mu_r": 1.0 }));
    }
    b.bc("ems", "pml_left", "pml", json!({ "x": "minus", "strength": 4.0, "order": 2.0 }));
    b.bc("ems", "pml_right", "pml", json!({ "x": "plus", "strength": 4.0, "order": 2.0 }));
    b
}

fn plate_with_hole() -> B {
    let mut b = B::new("Plate with a hole (solid mechanics, stress concentration)", 2, 0.012);
    b.rect("plate_box", [0.0, 0.0], [0.4, 0.2]);
    b.c("add_circle", json!({ "name": "hole", "center": [0.2, 0.1], "r": 0.04 }));
    b.c("add_difference", json!({ "name": "plate", "a": "plate_box", "b": ["hole"] }));
    b.overlay(&[("plate", "plate")]);
    b.c("add_size", json!({ "name": "size_hole", "sel": "kind(circle)", "h": 0.004 }));
    b.physics("sms", "static", Value::Null, 2);
    b.params("sms", json!({ "plane": "stress" }));
    steel_like(&mut b, "plate");
    b.bc("sms", "outer & plane(x = 0)", "fixed", json!({ "ux": 0.0 }));
    b.bc("sms", "outer & plane(z = 0)", "fixed", json!({ "uz": 0.0 }));
    b.bc("sms", "outer & plane(x = 0.4)", "load", json!({ "traction": [1e6, 0.0, 0.0] }));
    b
}

fn anchor_loss() -> B {
    let mut b = B::new("Beam on a PML substrate: anchor loss (solid eigenfrequencies, complex frequency and Q)", 2, 0.03);
    b.rect("substrate_box", [0.0, 0.0], [1.0, 0.5]);
    b.rect("pml_left_box", [0.0, 0.15], [0.2, 0.5]);
    b.rect("pml_right_box", [0.8, 0.15], [1.0, 0.5]);
    b.rect("pml_bottom_box", [0.2, 0.0], [0.8, 0.15]);
    b.rect("pml_bottom_left_box", [0.0, 0.0], [0.2, 0.15]);
    b.rect("pml_bottom_right_box", [0.8, 0.0], [1.0, 0.15]);
    b.rect("beam_box", [0.49, 0.5], [0.51, 0.8]);
    b.overlay(&[
        ("substrate", "substrate_box"),
        ("pml_left", "pml_left_box"),
        ("pml_right", "pml_right_box"),
        ("pml_bottom", "pml_bottom_box"),
        ("pml_bottom_left", "pml_bottom_left_box"),
        ("pml_bottom_right", "pml_bottom_right_box"),
        ("beam", "beam_box"),
    ]);
    b.physics("sms", "eigenfrequency", json!({ "n_modes": 6, "f_search": 0.0 }), 2);
    b.params("sms", json!({ "plane": "strain" }));
    for r in ["substrate", "pml_left", "pml_right", "pml_bottom", "pml_bottom_left", "pml_bottom_right", "beam"] {
        steel_like(&mut b, r);
    }
    for (name, x, z) in [
        ("pml_left", "minus", "off"),
        ("pml_right", "plus", "off"),
        ("pml_bottom", "off", "minus"),
        ("pml_bottom_left", "minus", "minus"),
        ("pml_bottom_right", "plus", "minus"),
    ] {
        b.bc("sms", name, "pml", json!({ "x": x, "z": z, "strength": 3.0, "order": 2.0 }));
    }
    // clamped on the outer edge of the layers: every boundary below the substrate surface
    b.bc("sms", "(outer & in_box([-1, 0, -1], [2, 0, 0.5000000001], within)) - plane(z = 0.5)", "fixed", FIXED());
    b
}

fn drilled_plate() -> B {
    let mut b = B::new("Plate with a drilled hole in 3D, the section swept along y (solid mechanics, stress concentration)", 2, 0.02);
    b.rect("plate_box", [0.0, 0.0], [0.4, 0.2]);
    b.c("add_circle", json!({ "name": "hole", "center": [0.2, 0.1], "r": 0.04 }));
    b.c("add_difference", json!({ "name": "plate", "a": "plate_box", "b": ["hole"] }));
    b.overlay(&[("plate", "plate")]);
    b.c("add_size", json!({ "name": "size_hole", "sel": "kind(circle)", "h": 0.01 }));
    b.mesh(json!({ "sweep": [{ "thickness": 0.04, "n": 4 }] }));
    b.physics("sms", "static", Value::Null, 2);
    steel_like(&mut b, "plate");
    b.bc("sms", "outer & plane(x = 0)", "fixed", FIXED());
    b.bc("sms", "outer & plane(x = 0.4)", "load", json!({ "traction": [1e6, 0.0, 0.0] }));
    b
}

/// The cross-section of a silicon strip (450 nm x 230 nm) on r-plane sapphire under air at 1550 nm: optical modes and acoustic modes.
fn sos_strip() -> B {
    let um = 1e-6;
    let (xin, hin, ztop, lp) = (6.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let (w_si, h_si, zone) = (450e-9, 230e-9, 1.5 * um);
    let x = xin + lp;
    let mut b = B::new("Silicon strip on sapphire at 1550 nm: optical modes (EM waves) and acoustic modes in plane strain (solid mechanics)", 2, 0.4 * um);
    b.rect("air_box", [-x, -hin - lp], [x, ztop + lp]);
    b.rect("sapphire_box", [-x, -hin - lp], [x, 0.0]);
    b.rect("sapphire_near_box", [-zone, -zone], [zone, 0.0]);
    b.rect("air_near_box", [-zone, 0.0], [zone, zone]);
    b.rect("silicon_box", [-w_si / 2.0, 0.0], [w_si / 2.0, h_si]);
    b.rect("window", [-xin, -hin], [xin, ztop]);
    b.overlay(&[("air", "air_box"), ("sapphire", "sapphire_box"), ("sapphire_near", "sapphire_near_box"), ("air_near", "air_near_box"), ("silicon", "silicon_box")]);
    b.c("add_pml", json!({ "name": "pml", "window": "window", "thickness": lp, "layers": 3.0, "grading": 2.0 }));
    b.c("add_size", json!({ "name": "size_near", "sel": "region(sapphire_near) | region(air_near)", "h": 0.08 * um }));
    b.c("add_size", json!({ "name": "size_si", "sel": "region(silicon)", "h": 40e-9 }));
    // r-plane sapphire: the c axis in the y-z plane, 57.6 degrees from the surface normal
    let t = 57.6f64.to_radians();
    let euler = femcore::material::euler_from_rotation(&[[1.0, 0.0, 0.0], [0.0, t.cos(), -t.sin()], [0.0, t.sin(), t.cos()]]);
    b.physics("ems", "modes", json!({ "freq": 299792458.0 / 1550e-9, "n_target": 2.4, "n_modes": 4 }), 2);
    b.material("ems", "sapphire*", json!({ "material": "sapphire", "euler_deg": euler }));
    b.material("ems", "sapphire_near*", json!({ "material": "sapphire", "euler_deg": euler }));
    b.material("ems", "silicon", json!("silicon"));
    b.material("ems", "air*", json!("air"));
    b.material("ems", "air_near*", json!("air"));
    b.c("add_pml_conditions", json!({ "physics": "ems", "strength": 3.0 }));
    b.physics("sms", "eigenfrequency", json!({ "n_modes": 8, "f_search": 9.47e9 }), 2);
    b.params("sms", json!({ "plane": "strain", "inactive": ["air*", "air_near*"] }));
    b.material("sms", "sapphire*", json!({ "material": "sapphire", "euler_deg": euler }));
    b.material("sms", "sapphire_near*", json!({ "material": "sapphire", "euler_deg": euler }));
    b.material("sms", "silicon", json!("silicon"));
    b.c("add_pml_conditions", json!({ "physics": "sms", "strength": 3.0 }));
    b.bc("sms", "outer & in_box([-1, 0, -1], [1, 0, 1e-12], within)", "fixed", FIXED());
    b.active("ems");
    b
}

fn box_cavity() -> B {
    let mut b = B::new("Rectangular cavity: 3D electromagnetic eigenmodes (40 x 30 x 20 mm, perfect conductor walls)", 2, 0.005);
    b.plate(0.04, 0.02);
    b.mesh(json!({ "sweep": [{ "thickness": 0.03, "n": 6 }] }));
    b.physics("ems", "cavity", json!({ "n_modes": 6, "f_search": 6.0e9 }), 2);
    b.material("ems", "plate", json!({ "eps_r": 1.0, "mu_r": 1.0 }));
    b
}

/// A silicon wire in oxide on a silicon-on-insulator film at 1550 nm; geometry from the old use case, written once as data.
fn soi(title: &str, wire: femcore::testkit::soi::SoiWire) -> B {
    let mut b = B { e: Engine::new() };
    let mut p = Project::new(title, 2);
    p.geometry = wire.model().list;
    p.mesh.h = wire.h_far;
    b.e.set_project(p);
    b.physics("ems", "modes", json!({ "freq": 299792458.0 / wire.wavelength, "n_target": wire.n_target, "n_modes": wire.n_modes }), 2);
    for (regions, mat) in [("BOX", "silica"), ("BOX_near", "silica"), ("Si", "silicon"), ("clad", "silica"), ("clad_near", "silica")] {
        b.material("ems", regions, json!(mat));
    }
    if wire.slab > 0.0 {
        b.material("ems", "Si_slab", json!("silicon"));
    }
    b
}

fn soi_wire() -> B {
    soi("Silicon wire on 300 nm SOI at 1550 nm: TE0 and TM0 modes with Stokes parameters and Poincare sphere (EM waves, 2D cross-section)", femcore::testkit::soi::SoiWire::default())
}

fn soi_hybrid_wire() -> B {
    soi(
        "300 nm SOI wire with a sloped left wall: hybrid modes tilted on the Poincare sphere (EM waves, 2D cross-section)",
        femcore::testkit::soi::SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), ..Default::default() },
    )
}

#[test]
#[ignore = "writes examples/ at the repository root"]
fn write_the_examples() {
    capacitor().save("01_capacitor");
    cantilever().save("02_cantilever");
    heated_plate().save("03_heated_plate");
    thermal_bar().save("04_thermal_bar");
    block_cantilever().save("05_block_cantilever");
    waveguide().save("06_wr90_waveguide");
    loaded_waveguide().save("07_loaded_waveguide");
    sapphire_cantilever().save("08_sapphire_cantilever");
    pml_slab().save("09_pml_slab");
    plate_with_hole().save("10_plate_with_hole");
    anchor_loss().save("11_anchor_loss");
    drilled_plate().save("12_drilled_plate");
    sos_strip().save("13_silicon_on_sapphire");
    box_cavity().save("14_box_cavity");
    soi_wire().save("15_soi_wire");
    soi_hybrid_wire().save("16_soi_hybrid_wire");
}
