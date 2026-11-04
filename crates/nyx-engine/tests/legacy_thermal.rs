//! Heat feeding a solid study (thermal expansion from a temperature field), through the window bridge `nyx_engine::legacy` (ADR 41): dfs runs
//! and sms run in their crates. These two tests need both physics, so they live here and not in a physics crate.
//! The coupled chain through the trait (temperature exported by dfs, imported by sms through `Ctx::inputs`) is tested in `nyx-sms/tests/contracts.rs`.
use femcore::material::*;
use femcore::problem::*;
use femcore::testkit::rect;
use nyx_physics as np;
use nyx_physics::testkit::golden;

fn dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// An old-shaped result as a `nyx_physics::SolveResult` (every flat field becomes a scalar field), for the golden harness.
fn to_new(r: &SolveResult) -> np::SolveResult {
    np::SolveResult {
        order: r.order as u8,
        dim: r.dim as u8,
        solutions: r
            .solutions
            .iter()
            .map(|s| np::Solution {
                label: s.label.clone(),
                value: s.value,
                eigenvalue: s.neff.map(|(a, b)| [a, b]),
                fields: s
                    .fields
                    .iter()
                    .map(|f| np::Field::scalar(f.name.clone(), f.unit.clone(), f.values.clone()))
                    .collect(),
                cell_fields: s.cell_fields.as_ref().map(np::cell_fields_of).unwrap_or_default(),
            })
            .collect(),
        scalars: r
            .scalars
            .iter()
            .map(|(n, v)| {
                (
                    n.clone(),
                    match v {
                        Scalar::Num(x) => np::Scalar::Num(*x),
                        Scalar::Text(t) => np::Scalar::Text(t.clone()),
                    },
                )
            })
            .collect(),
        ..Default::default()
    }
}

fn steel() -> Materials {
    let mut m = Materials::default();
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        m.default.insert(k.into(), v);
    }
    m
}

fn cantilever(plane: Plane, study: Study) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.materials = steel();
    s.plane = plane;
    s.study = study;
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s
}

fn clamped_left(spec: &mut ProblemSpec) {
    spec.bcs.push(Bc {
        on: None,
        labels: vec![1],
        kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] },
    });
}

fn maxdiff(a: &[f64], b: &[f64]) -> f64 {
    femcore::testkit::maxdiff(a, b)
}

#[test]
fn heat_result_drives_thermal_expansion() {
    let mesh = rect(1.0, 0.1, 10, 2);
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.materials.default.insert("k".into(), 1.0);
    heat.bcs = vec![Bc { on: None, labels: vec![1, 4], kind: BcKind::Temperature(393.15) }];
    let mut state = SolveState::default();
    nyx_engine::legacy::solve(&mesh, &heat, &mut state).unwrap();
    let mk = |source| {
        let mut spec = cantilever(Plane::Strain, Study::Static);
        for (k, v) in [("E", 70e9), ("nu", 0.33), ("alpha", 23e-6)] {
            spec.materials.default.insert(k.into(), v);
        }
        spec.thermal = Some(ThermalExpansion { source, t_ref: 293.15 });
        spec
    };
    let solve_solid =
        |spec: &ProblemSpec, state: &SolveState| nyx_engine::legacy::solve(&mesh, spec, &mut state.clone());
    let a = solve_solid(&mk(ThermalSource::FromHeat), &state).unwrap();
    let b = solve_solid(&mk(ThermalSource::Uniform(393.15)), &state).unwrap();
    let (ua, ub) = (&a.solutions[0].field("ux").unwrap().values, &b.solutions[0].field("ux").unwrap().values);
    assert!(maxdiff(ua, ub) < 1e-12);
    assert!(solve_solid(&mk(ThermalSource::FromHeat), &SolveState::default()).is_err());
}

#[test]
fn golden_solid_static_isotropic_thermal_from_heat() {
    let mesh = rect(1.0, 0.2, 10, 3);
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.order = 2;
    heat.materials.default.insert("k".into(), 5.0);
    heat.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Temperature(293.15) },
        Bc { on: None, labels: vec![4], kind: BcKind::Temperature(393.15) },
        Bc { on: None, labels: vec![2], kind: BcKind::Temperature(293.15) },
    ];
    let mut state = SolveState::default();
    nyx_engine::legacy::solve(&mesh, &heat, &mut state).unwrap();
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("alpha", 2.3e-5)] {
        s.materials.default.insert(k.into(), v);
    }
    s.materials.aniso.insert(
        1,
        Anisotropy { expansion: Some(Tensor::Diagonal { x: 2.0e-5, y: 3.0e-5, z: 2.6e-5 }), ..Default::default() },
    );
    s.thermal = Some(ThermalExpansion { source: ThermalSource::FromHeat, t_ref: 293.15 });
    clamped_left(&mut s);
    let r = nyx_engine::legacy::solve(&mesh, &s, &mut state).unwrap();
    golden(&dir(), "solid_static_thermal_from_heat", &to_new(&r), 1e-8);
}
