//! Frozen reference output of the 2D kernels, taken through `solve::solve` before a kernel is ported to the wafer axes.
//!
//! Each case is a wafer-frame problem. The golden file `tests/golden_frames/<case>.json` holds the names of every field of every solution, the
//! scalars and effective indices, and for every field a few numbers (sum of |v|, max |v|, the values at fixed dofs). A ported kernel must
//! reproduce them to round-off. `NYX_WRITE_GOLDEN=1 cargo test --test golden_frames` writes the files (`NYX_WRITE_GOLDEN=solid` only those named solid*) (do that only on a kernel that is known
//! to be right, i.e. before porting it).
mod common;
use femcore::material::*;
use femcore::problem::*;
use serde_json::{json, Value};
use std::collections::BTreeMap;

const C0: f64 = 299_792_458.0;

fn record(r: &SolveResult) -> Value {
    let mut sols = vec![];
    for s in &r.solutions {
        let mut fields = BTreeMap::new();
        for f in &s.fields {
            let v = &f.values;
            let n = v.len();
            let l1: f64 = v.iter().map(|x| x.abs()).sum();
            let mx = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            fields.insert(f.name.clone(), json!({"n": n, "l1": l1, "max": mx, "probe": [v[0], v[n / 3], v[n / 2], v[n - 1]]}));
        }
        sols.push(json!({"neff": s.neff.map(|(a, b)| json!([a, b])), "value": s.value, "fields": fields}));
    }
    let mut scalars = BTreeMap::new();
    for (k, v) in &r.scalars {
        scalars.insert(k.clone(), match v {
            Scalar::Num(x) => json!(x),
            Scalar::Text(t) => json!(t),
        });
    }
    json!({"solutions": sols, "scalars": scalars})
}

fn close(a: f64, b: f64, scale: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * scale.max(1e-300)
}

fn compare(path: &str, want: &Value, got: &Value, scale: f64, tol: f64, out: &mut Vec<String>) {
    match (want, got) {
        (Value::Number(a), Value::Number(b)) => {
            let (a, b) = (a.as_f64().unwrap(), b.as_f64().unwrap());
            if !close(a, b, scale.max(a.abs()), tol) {
                out.push(format!("{path}: {a:e} against {b:e}"));
            }
        }
        (Value::Array(a), Value::Array(b)) => {
            if a.len() != b.len() {
                out.push(format!("{path}: {} entries against {}", a.len(), b.len()));
            }
            for (i, (x, y)) in a.iter().zip(b).enumerate() {
                compare(&format!("{path}[{i}]"), x, y, scale, tol, out);
            }
        }
        (Value::Object(a), Value::Object(b)) => {
            let (ka, kb): (Vec<_>, Vec<_>) = (a.keys().collect(), b.keys().collect());
            if ka != kb {
                out.push(format!("{path}: names {ka:?} against {kb:?}"));
                return;
            }
            // the scale of a field's probe values is the field's own maximum
            let sc = a.get("max").and_then(|m| m.as_f64()).unwrap_or(scale);
            for k in a.keys() {
                compare(&format!("{path}.{k}"), &a[k], &b[k], sc, tol, out);
            }
        }
        _ => {
            if want != got {
                out.push(format!("{path}: {want} against {got}"));
            }
        }
    }
}

/// Checks (or, with NYX_WRITE_GOLDEN, writes) the record of `r` against the golden file `name`.
fn golden(name: &str, r: &SolveResult, tol: f64) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/golden_frames/{name}.json"));
    let got = record(r);
    // NYX_WRITE_GOLDEN=1 writes every file, NYX_WRITE_GOLDEN=<prefix> only the files whose name starts with the prefix
    if std::env::var("NYX_WRITE_GOLDEN").is_ok_and(|w| w == "1" || name.starts_with(&w)) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, serde_json::to_string_pretty(&got).unwrap() + "\n").unwrap();
        return;
    }
    let want: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {path:?}"))).unwrap();
    // the text of an effective index and Q = Re/(2|Im|) carry the round-off of an imaginary part that is zero: the numbers are in `neff`
    let strip = |v: &Value| -> Value {
        let mut v = v.clone();
        // a field that is zero up to round-off (the longitudinal displacement of a Brillouin mode is a quarter period out of phase, so its
        // real part vanishes) is noise: leave it out
        if name.starts_with("brillouin") {
            if let Some(sols) = v["solutions"].as_array_mut() {
                for sol in sols {
                    if let Some(fs) = sol["fields"].as_object_mut() {
                        let big = fs.values().filter_map(|f| f["max"].as_f64()).fold(0.0, f64::max);
                        fs.retain(|_, f| f["max"].as_f64().is_some_and(|m| m > 1e-9 * big));
                    }
                }
            }
        }
        if let Some(sc) = v["scalars"].as_object_mut() {
            sc.retain(|k, x| !(k.starts_with("Q ") || x.is_string() && k.starts_with("neff")));
            // a Stokes parameter that is zero prints with either sign
            for x in sc.values_mut() {
                if let Some(t) = x.as_str() {
                    let mut t = t.replace("-0.000", "+0.000");
                    // the Q of a lossless ring mode is 1e15 and more: the round-off of an imaginary part that is zero
                    if name.starts_with("brillouin_ring") {
                        if let Some(i) = t.find("| Q ") {
                            let end = t[i + 4..].find(' ').map_or(t.len(), |e| i + 4 + e);
                            t.replace_range(i + 4..end, "-");
                        }
                    }
                    *x = Value::String(t);
                }
            }
        }
        v
    };
    let (want, got) = (strip(&want), strip(&got));
    let mut diffs = vec![];
    compare(name, &want, &got, 1.0, tol, &mut diffs);
    if !diffs.is_empty() {
        // an eigenvector is fixed only up to its sign when two of its entries tie in modulus: accept a solution whose signed fields
        // are all negated (the other fields, scalars and the rest of the record must still agree)
        let mut gauged = got.clone();
        if let (Some(ws), Some(gs)) = (want["solutions"].as_array(), gauged["solutions"].as_array_mut()) {
            for (w, g) in ws.iter().zip(gs.iter_mut()) {
                let mut d = vec![];
                compare("", w, g, 1.0, tol, &mut d);
                let mut flipped = g.clone();
                if !d.is_empty() {
                    if let Some(fs) = flipped["fields"].as_object_mut() {
                        for (fname, f) in fs.iter_mut() {
                            if !fname.starts_with('|') && fname != "von Mises" {
                                if let Some(pr) = f["probe"].as_array_mut() {
                                    pr.iter_mut().for_each(|x| *x = json!(-x.as_f64().unwrap()));
                                }
                            }
                        }
                    }
                    let mut d2 = vec![];
                    compare("", w, &flipped, 1.0, tol, &mut d2);
                    if d2.is_empty() {
                        *g = flipped;
                    }
                }
            }
        }
        diffs.clear();
        compare(name, &want, &gauged, 1.0, tol, &mut diffs);
    }
    assert!(diffs.is_empty(), "{} differences from the golden file {name}:\n{}", diffs.len(), diffs.iter().take(12).cloned().collect::<Vec<_>>().join("\n"));
}

fn solve(mesh: &femcore::mesh::SimplexMesh, spec: &ProblemSpec) -> SolveResult {
    femcore::solve::solve(mesh, spec, &mut SolveState::default()).unwrap()
}

// ---------------------------------------------------------------------- scalar kernels

#[test]
fn golden_electrostatics_anisotropic_plates() {
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    // a generic orientation, so that every entry of the tensor (xz, xy, zy) is non-zero
    spec.materials.aniso.insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }), euler_deg: [30.0, 40.0, 10.0], ..Default::default() });
    spec.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    golden("es_aniso_plates", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_electrostatics_series_with_charge() {
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.materials.set(2, "rho_v", 1e-9);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }, Bc { on: None, labels: vec![3], kind: BcKind::SurfaceCharge(1e-10) }];
    golden("es_series_charge", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_anisotropic_stationary() {
    let mesh = common::rect(1.0, 0.4, 10, 4);
    let mut spec = ProblemSpec::new(Physics::Heat);
    spec.materials.aniso.insert(1, Anisotropy { conductivity: Some(Tensor::Diagonal { x: 30.0, y: 5.0, z: 2.0 }), euler_deg: [10.0, 30.0, 50.0], ..Default::default() });
    spec.materials.default.insert("Q".into(), 1e4);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(300.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Convection { h: 20.0, text: 280.0 } }, Bc { on: None, labels: vec![3], kind: BcKind::Flux(500.0) }];
    golden("heat_aniso_stationary", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_transient_dense() {
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let mut spec = ProblemSpec::new(Physics::Heat);
    spec.materials.default.insert("k".into(), 2.0);
    spec.materials.default.insert("rho".into(), 1.0);
    spec.materials.default.insert("Cp".into(), 1.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    spec.study = Study::Transient { t_end: 0.2, n_out: 4, t0: "sin(pi*x)".into() };
    golden("heat_transient", &solve(&mesh, &spec), 1e-9);
}

// ---------------------------------------------------------------------- optical modes

#[test]
fn golden_soi_wire_modes() {
    let wire = femcore::usecases::soi::SoiWire { n_modes: 2, ..Default::default() };
    golden("emw_soi_wire", &wire.solve().unwrap(), 1e-8);
}

#[test]
fn golden_soi_wire_with_sloped_wall_hybrid_modes() {
    let wire = femcore::usecases::soi::SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), n_modes: 2, ..Default::default() };
    golden("emw_soi_hybrid", &wire.solve().unwrap(), 1e-8);
}

fn crystal(tilt_deg: f64) -> Anisotropy {
    Anisotropy { permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }), euler_deg: [0.0, tilt_deg, 0.0], ..Default::default() }
}

#[test]
fn golden_tilted_crystal_in_a_pec_guide() {
    // the coupled (xz/yz) path of the complex solver, then the uncoupled real path
    for (name, tilt) in [("emw_pec_tilted", 40.0), ("emw_pec_straight", 0.0)] {
        let mesh = common::rect(0.03, 0.015, 12, 6);
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = 2;
        s.materials.aniso.insert(1, crystal(tilt));
        s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 3 };
        golden(name, &solve(&mesh, &s), 1e-7);
    }
}

#[test]
fn golden_pec_guide_with_a_pml() {
    let mesh = common::rect_split_x(0.03, 0.015, 0.02, 0.0025);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.set(1, "eps_r", 4.0);
    s.materials.set(2, "eps_r", 1.0);
    // absorbing layer on the right (domain 2) along x
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Pml::default() });
    s.study = Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 3 };
    let _ = C0;
    golden("emw_pml_x", &solve(&mesh, &s), 1e-7);
}

// ---------------------------------------------------------------------- solid mechanics

fn sapphire_euler() -> [f64; 3] {
    [30.0, 40.0, 10.0]
}

fn sapphire() -> Anisotropy {
    let st = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    Anisotropy { stiffness: Some(st), expansion: Some(Tensor::Diagonal { x: 5.0e-6, y: 5.0e-6, z: 6.0e-6 }), euler_deg: sapphire_euler(), ..Default::default() }
}

fn clamped_left(spec: &mut ProblemSpec) {
    spec.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
}

#[test]
fn golden_solid_static_tilted_sapphire_plane_strain() {
    // a generic crystal orientation turns every entry of the stiffness; loads along x and z, a sliding wall at the bottom
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.plane = Plane::Strain;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    clamped_left(&mut s);
    s.bcs.push(Bc { on: None, labels: vec![4], kind: BcKind::Load([2.0e6, 0.0, -1.5e6]) });
    s.bcs.push(Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0, 0.0, 1.0e-7] } });
    golden("solid_static_sapphire_strain", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_static_plane_stress_with_thermal_expansion() {
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.plane = Plane::Stress;
    s.materials.aniso.insert(1, sapphire());
    s.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
    clamped_left(&mut s);
    golden("solid_static_sapphire_stress_thermal", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_static_isotropic_thermal_from_heat() {
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.order = 2;
    heat.materials.default.insert("k".into(), 5.0);
    heat.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(293.15) }, Bc { on: None, labels: vec![4], kind: BcKind::Temperature(393.15) }, Bc { on: None, labels: vec![2], kind: BcKind::Temperature(293.15) }];
    let mut state = SolveState::default();
    femcore::solve::solve(&mesh, &heat, &mut state).unwrap();
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("alpha", 2.3e-5)] {
        s.materials.default.insert(k.into(), v);
    }
    s.materials.aniso.insert(1, Anisotropy { expansion: Some(Tensor::Diagonal { x: 2.0e-5, y: 3.0e-5, z: 2.6e-5 }), ..Default::default() });
    s.thermal = Some(ThermalExpansion { source: ThermalSource::FromHeat, t_ref: 293.15 });
    clamped_left(&mut s);
    let r = femcore::solve::solve(&mesh, &s, &mut state).unwrap();
    golden("solid_static_thermal_from_heat", &r, 1e-8);
}

#[test]
fn golden_solid_eigenfrequencies_tilted_sapphire() {
    let mesh = common::rect(1.0e-3, 0.2e-3, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    clamped_left(&mut s);
    s.study = Study::Eigenfrequency { n_modes: 5, f_search: 0.0 };
    golden("solid_eigen_sapphire", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_eigenfrequencies_floquet_cell() {
    // a phononic cell: the sides x = 0 and x = lx are periodic with a phase
    let mesh = common::rect(1.0e-3, 0.3e-3, 8, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 1.2 }];
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 1.0e6 };
    golden("solid_eigen_floquet", &solve(&mesh, &s), 1e-7);
}

#[test]
fn golden_solid_eigenfrequencies_with_a_pml() {
    // a clamped isotropic bar whose right end is an absorbing layer along x
    let mesh = common::rect(1.0, 0.1, 30, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.plane = Plane::Strain;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    clamped_left(&mut s);
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    golden("solid_eigen_pml_x", &solve(&mesh, &s), 1e-7);
}

#[test]
fn golden_solid_eigenfrequencies_with_a_pml_along_z() {
    // a plate on a substrate: the layer is the lower part of a tall section and absorbs along the surface normal z
    let mesh = common::rect(0.2, 1.0, 3, 30);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.plane = Plane::Strain;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { on: None, labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.materials.pml.insert(1, Pml { z: PmlDir::Minus, strength: 0.1, order: 0.0, ..Pml::default() });
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    golden("solid_eigen_pml_z", &solve(&mesh, &s), 1e-7);
}

// ---------------------------------------------------------------------- solid mechanics with an out-of-plane wave number

#[test]
fn golden_solid_kz_tilted_sapphire() {
    // the cross-section of a waveguide carries a wave exp(-j k y) along +y: a generic crystal orientation couples it to x and z
    let mesh = common::rect(1.0e-3, 0.2e-3, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    clamped_left(&mut s);
    s.out_of_plane = Some(3.0e3);
    s.study = Study::Eigenfrequency { n_modes: 5, f_search: 0.0 };
    golden("solid_kz_sapphire", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_kz_sapphire_component_constraints() {
    // a wall that holds only x and z (y free) and a wall that holds only y: the components of a constraint are wafer axes
    let mesh = common::rect(1.0e-3, 0.2e-3, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } });
    s.bcs.push(Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } });
    s.out_of_plane = Some(-2.0e3);
    s.study = Study::Eigenfrequency { n_modes: 5, f_search: 0.0 };
    golden("solid_kz_sapphire_constraints", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_kz_with_a_pml_along_x() {
    let mesh = common::rect(1.0, 0.1, 30, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    clamped_left(&mut s);
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    s.out_of_plane = Some(4.0);
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    golden("solid_kz_pml_x", &solve(&mesh, &s), 1e-7);
}

#[test]
fn golden_solid_kz_sapphire_with_a_pml_along_z() {
    let mesh = common::rect(0.2, 1.0, 3, 30);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { on: None, labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.materials.pml.insert(1, Pml { z: PmlDir::Minus, strength: 0.1, order: 0.0, ..Pml::default() });
    s.out_of_plane = Some(3.0);
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    golden("solid_kz_pml_z", &solve(&mesh, &s), 1e-7);
}

// ---------------------------------------------------------------------- Brillouin scattering of a straight waveguide

/// Rotates every crystal of the problem by half a turn about the wafer axis y. Even-rank tensors do not change under an inversion, so this
/// is the mirror image of the crystal across the x-z plane. The golden files of the Brillouin cases were written with the kernels that
/// still propagated along -y: a pump along -y in a crystal is the pump along +y in its mirror image (`NYX_GOLDEN_LEGACY_MIRROR=1`).
fn mirror_crystals_in_y(spec: &mut ProblemSpec) {
    let flip = |e: [f64; 3]| -> [f64; 3] {
        let half_turn_y: M3 = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
        let r = rotation(e);
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    m[i][j] += half_turn_y[i][k] * r[k][j];
                }
            }
        }
        euler_from_rotation(&m)
    };
    for a in spec.materials.aniso.values_mut() {
        a.euler_deg = flip(a.euler_deg);
    }
    if let Some(b) = spec.brillouin.as_mut() {
        for p in b.photoelastic.values_mut() {
            p.euler_deg = flip(p.euler_deg);
        }
        for a in b.viscosity.values_mut() {
            a.euler_deg = flip(a.euler_deg);
        }
    }
}

fn wire_mesh() -> femcore::mesh::SimplexMesh {
    let um = 1e-6;
    let rect = |x0: f64, z0: f64, x1: f64, z1: f64, d: i32, h: Option<f64>| common::part(x0 * um, z0 * um, x1 * um, z1 * um, d, h.map(|v| v * um));
    common::parts_mesh(&[rect(-1.2, -1.0, 1.2, 1.0, 1, None), rect(-0.225, -0.11, 0.225, 0.11, 2, Some(0.04))], 0.3 * um)
}

/// A silicon wire in silica with the elastic, optical and photoelastic data of the crystal in the orientation `euler`.
fn wire_brillouin(mesh: &femcore::mesh::SimplexMesh, euler: [f64; 3], backward: bool, n_optical: usize) -> ProblemSpec {
    let mut sp = ProblemSpec::new(Physics::Brillouin);
    sp.order = 2;
    for d in mesh.domain_labels() {
        if d == 2 {
            sp.materials.aniso.insert(d, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: euler, ..Default::default() });
            sp.materials.set(d, "rho", 2329.0);
        } else {
            sp.materials.set(d, "rho", 2203.0);
        }
        sp.materials.set(d, "E", 73.1e9);
        sp.materials.set(d, "nu", 0.17);
    }
    let mut labels: Vec<i32> = (0..mesh.n_facets()).map(|f| mesh.facet_labels[f]).collect();
    labels.sort();
    labels.dedup();
    sp.bcs = vec![Bc { on: None, labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    sp.study = Study::Eigenfrequency { n_modes: 8, f_search: 16e9 };
    let mut b = BrillouinSpec::default();
    b.wavelength = 1.55e-6;
    b.n_target = 2.35;
    b.n_optical_modes = n_optical;
    b.backward = backward;
    b.pump_modes = vec![1];
    b.stokes_modes = vec![if backward { 1 } else { 2 }];
    b.refractive_index.insert(1, 1.444);
    b.refractive_index.insert(2, 3.476);
    b.photoelastic.insert(2, PhotoSpec { kind: PhotoKind::Cubic { p11: -0.094, p12: 0.017, p44: -0.051 }, euler_deg: [0.0; 3], follow_stiffness: true });
    b.photoelastic.insert(1, PhotoSpec { kind: PhotoKind::Isotropic { p11: 0.121, p12: 0.270 }, euler_deg: [0.0; 3], follow_stiffness: false });
    b.moving = vec![2];
    b.viscosity.insert(2, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 1.0e-3, c12: 0.4e-3, c44: 0.5e-3 }), euler_deg: euler, ..Default::default() });
    sp.brillouin = Some(b);
    sp
}

fn golden_brillouin(name: &str, mut spec: ProblemSpec) {
    if std::env::var("NYX_GOLDEN_LEGACY_MIRROR").is_ok() {
        mirror_crystals_in_y(&mut spec);
    }
    golden(name, &solve(&wire_mesh(), &spec), 1e-6);
}

#[test]
fn golden_brillouin_wire_in_a_cubic_crystal() {
    // silicon along x, y and z: backward scattering of the quasi-TE mode (a crystal this symmetric is its own mirror image)
    golden_brillouin("brillouin_wire_cubic", wire_brillouin(&wire_mesh(), [0.0; 3], true, 1));
}

#[test]
fn golden_brillouin_wire_in_a_tilted_crystal() {
    golden_brillouin("brillouin_wire_tilted", wire_brillouin(&wire_mesh(), [30.0, 40.0, 10.0], true, 1));
}

#[test]
fn golden_brillouin_forward_scattering_in_a_tilted_crystal() {
    golden_brillouin("brillouin_wire_forward", wire_brillouin(&wire_mesh(), [30.0, 40.0, 10.0], false, 2));
}

// ---------------------------------------------------------------------- Brillouin scattering of a ring

/// The cross-section of `wire_mesh` on a ring of radius 20 mm (mesh x = r, mesh z = z).
fn ring_wire_mesh(radius: f64) -> femcore::mesh::SimplexMesh {
    let mesh = wire_mesh();
    let mut pts = mesh.points.clone();
    for i in 0..mesh.n_points() {
        pts[2 * i] += radius;
    }
    femcore::mesh::SimplexMesh::new(2, pts, mesh.cells.clone(), mesh.cell_labels.clone(), mesh.facets.clone(), mesh.facet_labels.clone())
}

fn golden_ring(name: &str, euler: [f64; 3]) {
    let radius = 20e-3;
    let mesh = ring_wire_mesh(radius);
    let mut spec = wire_brillouin(&mesh, euler, true, 1);
    spec.brillouin.as_mut().unwrap().viscosity.clear();
    spec.brillouin.as_mut().unwrap().ring_radius = Some(radius);
    golden(name, &solve(&mesh, &spec), 1e-6);
}

#[test]
fn golden_ring_brillouin_cubic_crystal() {
    golden_ring("brillouin_ring_cubic", [0.0; 3]);
}

#[test]
fn golden_ring_brillouin_tilted_crystal() {
    // a crystal ring is averaged over rotations about the ring axis, so its orientation matters only through the angle to the normal
    golden_ring("brillouin_ring_tilted", [30.0, 40.0, 10.0]);
}
