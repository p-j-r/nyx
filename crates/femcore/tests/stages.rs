//! The kernels are three stages: assemble (matrices, constraints), solve (eigenpairs or solution), fields (plot fields and the
//! same quantities per cell). The one-call `solve_*` is their composition; these tests drive the stages on their own.
mod common;
use femcore::material::*;
use femcore::problem::*;

#[test]
fn electrostatics_stages_compose_to_the_one_call_solver() {
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    let whole = femcore::dfs::solve_electrostatics(&mesh, &spec).unwrap();
    let asm = femcore::dfs::assemble_electrostatics(&mesh, &spec).unwrap();
    assert_eq!(asm.k.nrows, mesh.n_dofs(spec.order));
    let (v, _) = femcore::dfs::solve_potential(&mesh, &spec, &asm).unwrap();
    let (fields, cells) = femcore::dfs::electrostatics_fields(&mesh, spec.order, &v);
    let want = &whole.solutions[0];
    assert_eq!(want.field("V").unwrap().values, v);
    assert_eq!(fields.len(), want.fields.len());
    // the flux is discontinuous across the interface between the two dielectrics only in the normal direction: the per-cell
    // field keeps the jump that the node average hides
    let ex = cells.get("Ex").unwrap();
    let avg = cells.node_average(&mesh, "Ex").unwrap();
    assert_eq!(avg, want.field("Ex").unwrap().values);
    assert_eq!(ex.re.len(), mesh.n_cells() * mesh.nloc(spec.order));
    let spread = ex.re.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - ex.re.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(spread > 0.0);
    // no component named after a mesh axis: the cross-section has x and z
    assert!(cells.get("Ey").is_none() && cells.get("Ez").is_some());
}

#[test]
fn the_per_cell_flux_jumps_where_the_node_average_does_not() {
    // two dielectrics side by side across x: D_x is continuous, so E_x jumps by the ratio of the permittivities at the interface
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    let r = femcore::dfs::solve_electrostatics(&mesh, &spec).unwrap();
    let cf = r.solutions[0].cell_fields.as_ref().unwrap();
    let ex = cf.get("Ex").unwrap();
    let nloc = cf.nloc;
    // the largest and the smallest cell mean of |Ex| differ by about the permittivity ratio (4)
    let mean: Vec<f64> = (0..mesh.n_cells()).map(|c| ex.re[c * nloc..(c + 1) * nloc].iter().sum::<f64>() / nloc as f64).collect();
    let hi = mean.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let lo = mean.iter().filter(|x| x.abs() > 1e-9 * hi).fold(f64::INFINITY, |m, x| m.min(x.abs()));
    assert!((hi / lo - 4.0).abs() < 0.2, "ratio {}", hi / lo);
}

#[test]
fn heat_stages_give_the_temperature_and_the_flux_per_cell() {
    let mesh = common::rect(1.0, 0.4, 10, 4);
    let mut spec = ProblemSpec::new(Physics::Heat);
    spec.materials.default.insert("k".into(), 3.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(300.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Temperature(280.0) }];
    let sys = femcore::dfs::assemble_heat(&mesh, &spec).unwrap();
    let st = femcore::dfs::solve_heat_states(&mesh, &spec, &sys).unwrap();
    assert!(st.times.is_empty() && st.temps.len() == 1);
    let (fields, cells) = femcore::dfs::heat_fields(&mesh, &st.temps[0], &sys.ktensor, spec.order);
    assert!(fields.iter().any(|f| f.name == "qx"));
    // a linear temperature drop of 20 K over 1 m in k = 3: q_x = 60 W/m^2 in every cell
    let q = cells.get("qx").unwrap();
    assert!(q.re.iter().all(|v| (v - 60.0).abs() < 1e-6), "{:?}", &q.re[..4]);
}

#[test]
fn optical_mode_stages_compose_to_the_one_call_solver() {
    use femcore::ems::*;
    let mesh = common::rect(1.0e-6, 0.5e-6, 12, 6);
    let nc = mesh.n_cells();
    let eps = vec![EpsCell { xx: 4.0, xy: 0.0, yy: 4.0, zz: 4.0 }; nc];
    let mu = vec![1.0; nc];
    let whole = solve_modes_tensor(&mesh, 2.0e14, &eps, &mu, 1.5, false, 3, 2).unwrap();
    let asm = assemble_modes(&mesh, 2.0e14, &eps, &mu, 2).unwrap();
    let staged = solve_mode_pairs(&asm, 1.5, false, 3).unwrap();
    assert_eq!(whole.neff, staged.neff);
    let (fields, cells) = mode_fields_cells(&mesh, &staged, 0, 2.0e14, 2);
    assert_eq!(fields.len(), 6);
    // the components carry the wafer names: x lateral, z normal, y along the guide
    for name in ["Ex", "Ez", "Ey"] {
        let avg = cells.node_average(&mesh, name).unwrap();
        let plot = fields.iter().find(|f| f.name == name).unwrap();
        let scale = plot.values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(avg.iter().zip(&plot.values).all(|(a, b)| (a - b).abs() < 1e-12 * scale), "{name}");
    }
}

#[test]
fn elastic_cross_section_stages_compose_to_the_one_call_solver() {
    use femcore::sms_kz::*;
    let mesh = common::rect(1.0e-3, 0.2e-3, 8, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    let st = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    s.materials.aniso.insert(1, Anisotropy { stiffness: Some(st), euler_deg: [30.0, 40.0, 10.0], ..Default::default() });
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    let (whole, wm) = solve_cross_section(&mesh, &s, CrossSection::Kz(3.0e3)).unwrap();
    let asm = assemble_cross_section(&mesh, &s, CrossSection::Kz(3.0e3)).unwrap();
    assert_eq!(asm.kf.nrows, asm.free.len());
    let solved = solve_cross_section_modes(&asm).unwrap();
    assert_eq!(solved.modes.freqs, wm.freqs);
    let r = cross_section_fields(&mesh, &asm, &solved).unwrap();
    assert_eq!(r.solutions.len(), whole.solutions.len());
    // the three displacement components come per cell, complex, under wafer names (x, z in the plane; y along the wave)
    let cf = r.solutions[0].cell_fields.as_ref().unwrap();
    let names: Vec<&str> = cf.fields.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["ux", "uz", "uy"]);
    assert!(cf.get("uy").unwrap().im.iter().any(|x| x.abs() > 0.0));
}
