//! The kernels are three stages: assemble (matrices, constraints), solve (eigenpairs or solution), fields (plot fields and the
//! same quantities per cell). The one-call `solve_*` is their composition; these tests drive the stages on their own.
mod common;
use femcore::material::*;
use femcore::problem::*;

#[test]
fn elastic_cross_section_stages_compose_to_the_one_call_solver() {
    use femcore::sms_kz::*;
    // a ring of radius 5 mm with a 1 mm x 0.2 mm section: the mesh is the (r, z) half plane
    let flat = common::rect(1.0e-3, 0.2e-3, 8, 3);
    let mut pts = flat.points.clone();
    for i in 0..flat.n_points() {
        pts[2 * i] += 5.0e-3;
    }
    let mesh = femcore::mesh::SimplexMesh::new(2, pts, flat.cells.clone(), flat.cell_labels.clone(), flat.facets.clone(), flat.facet_labels.clone());
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    let st = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    s.materials.aniso.insert(1, Anisotropy { stiffness: Some(st), euler_deg: [30.0, 40.0, 10.0], ..Default::default() });
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.0 };
    let (whole, wm) = solve_cross_section(&mesh, &s, CrossSection::Axisymmetric(2.0)).unwrap();
    let asm = assemble_cross_section(&mesh, &s, CrossSection::Axisymmetric(2.0)).unwrap();
    assert_eq!(asm.kf.nrows(), asm.free.len());
    let solved = solve_cross_section_modes(&asm).unwrap();
    assert_eq!(solved.modes.freqs, wm.freqs);
    let r = cross_section_fields(&mesh, &asm, &solved).unwrap();
    assert_eq!(r.solutions.len(), whole.solutions.len());
    // the three displacement components come per cell, complex, under wafer names (x, z in the plane; y along the ring)
    let cf = r.solutions[0].cell_fields.as_ref().unwrap();
    let names: Vec<&str> = cf.fields.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["ux", "uz", "uy"]);
    assert!(cf.get("uy").unwrap().im.iter().any(|x| x.abs() > 0.0));
}

fn same_solution_fields(a: &SolveResult, b: &SolveResult) {
    assert_eq!(a.solutions.len(), b.solutions.len());
    for (x, y) in a.solutions.iter().zip(&b.solutions) {
        assert_eq!(x.label, y.label);
        for (f, g) in x.fields.iter().zip(&y.fields) {
            assert_eq!(f.name, g.name);
            assert_eq!(f.values, g.values);
        }
        assert_eq!(x.cell_fields.is_some(), y.cell_fields.is_some());
    }
}

#[test]
fn solid_pml_stages_compose_to_the_one_call_solver() {
    let mesh = common::rect(1.0, 0.1, 20, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    s.plane = Plane::Strain;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    s.study = Study::Eigenfrequency { n_modes: 3, f_search: 0.0 };
    let whole = femcore::sms_pml::solve_solid_pml(&mesh, &s).unwrap();
    let asm = femcore::sms_pml::assemble_solid_pml(&mesh, &s).unwrap();
    assert_eq!(asm.kff.nrows(), asm.free.len());
    let solved = femcore::sms_pml::solve_pml_modes(&asm).unwrap();
    let staged = femcore::sms_pml::pml_fields(&mesh, &asm, &solved).unwrap();
    same_solution_fields(&whole, &staged);
    let cf = staged.solutions[0].cell_fields.as_ref().unwrap();
    assert_eq!(cf.fields.iter().map(|f| f.name.as_str()).collect::<Vec<_>>(), ["ux", "uz"]);
    // the layer makes the frequencies complex
    assert!(solved.freqs.iter().any(|(_, f)| f.im.abs() > 0.0));
}

#[test]
fn solid_static_stages_give_the_displacement_and_the_stress_per_cell() {
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    for (k, v) in [("E", 70e9), ("nu", 0.3)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs.push(Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.bcs.push(Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1.0e6]) });
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let cf = r.solutions[0].cell_fields.as_ref().unwrap();
    for name in ["ux", "uz", "von Mises", "sxx", "szz", "sxz", "syy"] {
        assert!(cf.get(name).is_some(), "{name} in {:?}", cf.fields.iter().map(|f| &f.name).collect::<Vec<_>>());
    }
    assert!(cf.get("uy").is_none());
}

#[test]
fn floquet_stages_compose_to_the_one_call_solver() {
    let mesh = common::rect(1.0e-3, 0.3e-3, 8, 3);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 2;
    for (k, v) in [("E", 70e9), ("nu", 0.3), ("rho", 2700.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 1.2 }];
    s.study = Study::Eigenfrequency { n_modes: 3, f_search: 1.0e6 };
    let whole = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    // the stiffness and mass are the first stage of the plain kernel; drive the Floquet stages with the same matrices
    let st = femcore::sms::setup(&mesh, &s).unwrap();
    let k = femcore::sms::assemble_solid(&mesh, &st.c6, s.order, s.plane, s.thickness);
    let rho = s.materials.cell_values(&mesh, "rho", None).unwrap();
    let rho_t: Vec<f64> = rho.iter().map(|r| r * s.thickness).collect();
    let m1 = femcore::fem::assemble_mass(&mesh, s.order, &rho_t);
    let m = femcore::linalg::csr::RCsr::block_diag(&[&m1, &m1]);
    let asm = femcore::sms_floquet::assemble_floquet(&mesh, &s, &st, &k, &m, 3, 1.0e6).unwrap();
    assert!(asm.nred < asm.n);
    let solved = femcore::sms_floquet::solve_floquet_modes(&asm).unwrap();
    let staged = femcore::sms_floquet::floquet_fields(&mesh, &asm, &solved).unwrap();
    same_solution_fields(&whole, &staged);
}
