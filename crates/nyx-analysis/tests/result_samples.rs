//! The polarization analysis as a function of a `SolveResult`: what it needs from the result (frame, cell fields), what it leaves out (the
//! absorbing domains) and how it refuses a result that does not state a frame. The results here are made by hand on a small mesh, so the
//! numbers are known: the physics crates' own tests (`nyx-ems`: goldens, closed forms) check real modes.
use nyx_analysis::{polarization, polarization_applies, Analysis, PolarizationAnalysis};
use nyx_physics::{CellField, FieldFrame, Mesh, Scalar, Solution, SolveResult};

/// A P1 result of one mode on a 2D mesh: the cell fields `Ex` (lateral) and `Ez` (normal), and `Ey` (propagation) zero, all uniform.
fn result_on(
    mesh: &femcore::mesh::SimplexMesh,
    lat: [f64; 2],
    nor: [f64; 2],
    frame: Option<FieldFrame>,
) -> SolveResult {
    let nc = mesh.n_cells();
    let field = |name: &str, v: [f64; 2]| CellField {
        name: name.into(),
        unit: "a.u.".into(),
        nloc: 3,
        re: vec![v[0]; nc * 3],
        im: vec![v[1]; nc * 3],
    };
    let mut scalars = vec![];
    if let Some(f) = frame {
        scalars.push(f.scalar());
    }
    SolveResult {
        physics: "ems".into(),
        study: "modes".into(),
        order: 1,
        dim: 2,
        solutions: vec![Solution {
            label: "mode 1".into(),
            cell_fields: vec![field("Ex", lat), field("Ez", nor), field("Ey", [0.0, 0.0])],
            eigenvalue: Some([2.5, 0.0]),
            ..Default::default()
        }],
        scalars,
        ..Default::default()
    }
}

fn frame(outside: Vec<i32>) -> FieldFrame {
    FieldFrame { lateral: 0, normal: 2, propagation: 1, left_handed: true, outside }
}

#[test]
fn a_lateral_field_is_quasi_te_and_the_result_must_state_its_frame() {
    let sm = femcore::testkit::rect(1.0, 0.5, 4, 2);
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let r = result_on(&sm, [1.0, 0.0], [0.0, 0.0], Some(frame(vec![])));
    assert!(polarization_applies(&r));
    let a = polarization(&mesh, &r).unwrap();
    let m = a.rows[0].report.as_ref().unwrap();
    assert_eq!(m.class.label(), "quasi-TE (x)");
    assert!((m.te_fraction - 1.0).abs() < 1e-12 && m.long_fraction.abs() < 1e-12 && (m.dop - 1.0).abs() < 1e-12);
    assert_eq!(a.value_label, "neff");
    assert_eq!(a.rows[0].neff, (2.5, 0.0));
    // without a stated frame the analysis does not apply and says why
    let bare = result_on(&sm, [1.0, 0.0], [0.0, 0.0], None);
    assert!(!polarization_applies(&bare));
    let e = polarization(&mesh, &bare).unwrap_err().to_string();
    assert!(e.contains("field frame"), "{e}");
    // the same through the trait: plain JSON, with the table
    let v = PolarizationAnalysis.run(&mesh, &r).unwrap();
    assert_eq!(v["rows"][0]["report"]["class"], "quasi-TE (x)");
    assert!(v["table"].as_str().unwrap().contains("quasi-TE (x)"));
}

#[test]
fn a_90_degree_phase_between_the_components_is_circular_and_the_sense_follows_the_frame() {
    let sm = femcore::testkit::rect(1.0, 0.5, 4, 2);
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let left = polarization(&mesh, &result_on(&sm, [1.0, 0.0], [0.0, 1.0], Some(frame(vec![])))).unwrap();
    let mut right_frame = frame(vec![]);
    right_frame.left_handed = false;
    let right = polarization(&mesh, &result_on(&sm, [1.0, 0.0], [0.0, 1.0], Some(right_frame))).unwrap();
    let (l, r) = (left.rows[0].report.as_ref().unwrap(), right.rows[0].report.as_ref().unwrap());
    assert!(l.sphere[2].abs() > 0.99 && r.sphere[2].abs() > 0.99);
    assert!(
        l.sphere[2] * r.sphere[2] < 0.0,
        "the other handedness reads the sense of rotation the other way: {:?} {:?}",
        l.sphere,
        r.sphere
    );
}

#[test]
fn cells_of_the_absorbing_domains_are_left_out() {
    // two halves: domain 1 carries a lateral field, domain 2 (the absorbing layer) a normal one; leaving domain 2 out gives a pure TE mode
    let sm = femcore::testkit::rect_split_x(1.0, 0.5, 0.5, 0.1);
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let nc = sm.n_cells();
    let mut r = result_on(&sm, [0.0, 0.0], [0.0, 0.0], Some(frame(vec![])));
    for f in r.solutions[0].cell_fields.iter_mut() {
        for c in 0..nc {
            let v = if sm.cell_labels[c] == 2 { (f.name == "Ez") as u8 as f64 } else { (f.name == "Ex") as u8 as f64 };
            f.re[c * 3..c * 3 + 3].fill(v);
            f.im[c * 3..c * 3 + 3].fill(0.0);
        }
    }
    let all = polarization(&mesh, &r).unwrap();
    assert!(all.rows[0].report.as_ref().unwrap().te_fraction < 0.99, "the two halves mix");
    r.scalars.retain(|(n, _)| n != nyx_physics::FRAME_KEY);
    r.scalars.push((nyx_physics::FRAME_KEY.into(), Scalar::Text(frame(vec![2]).to_text())));
    let only = polarization(&mesh, &r).unwrap();
    let m = only.rows[0].report.as_ref().unwrap();
    assert!((m.te_fraction - 1.0).abs() < 1e-12 && m.class.label() == "quasi-TE (x)", "{m:?}");
}

#[test]
fn a_cell_field_of_another_mesh_is_refused_by_name() {
    let sm = femcore::testkit::rect(1.0, 0.5, 4, 2);
    let other = femcore::testkit::rect(1.0, 0.5, 6, 3);
    let r = result_on(&sm, [1.0, 0.0], [0.0, 0.0], Some(frame(vec![])));
    let e = polarization(&Mesh::from_solver_mesh(&other).unwrap(), &r).unwrap_err().to_string();
    assert!(e.contains("cell field 'Ex'") && e.contains("this mesh"), "{e}");
}
