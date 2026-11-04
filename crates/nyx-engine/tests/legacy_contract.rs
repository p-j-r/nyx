//! Contract `solve_result` (ADR 19) of the old result shape the window reads, through the bridge `nyx_engine::legacy` (ADR 41): the structure of a
//! `SolveResult`, the refusal of unsupported orders and non-finite input, and the cancel string. (Moved from the femcore contract test when dfs
//! moved into its crate; the new result shape is pinned by `solve_result.v2.json` in `nyx-physics`.) Deleted with the bridge in phase 7.
use femcore::problem::*;

fn electrostatics_case() -> (femcore::mesh::SimplexMesh, ProblemSpec) {
    let mesh = femcore::testkit::rect(1.0, 1.0, 8, 8);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 4.0);
    spec.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) },
        Bc { on: None, labels: vec![4], kind: BcKind::Ground },
    ];
    (mesh, spec)
}

#[test]
fn solve_result_has_the_documented_structure() {
    // `SolveResult` is not serialisable; h5io, nyx-store::runs::metrics_of and the MCP digests each read it by field, so its structure is the contract
    let (mesh, spec) = electrostatics_case();
    let mut state = SolveState::default();
    let r = nyx_engine::legacy::solve(&mesh, &spec, &mut state).unwrap();
    assert_eq!(r.order, spec.order);
    assert_eq!(r.dim, 2);
    assert_eq!(r.dofs, mesh.n_dofs(spec.order));
    assert!(r.seconds >= 0.0);
    assert!(!r.log.is_empty(), "the stage log is appended to SolveResult::log");
    assert_eq!(r.solutions.len(), 1);
    let s = &r.solutions[0];
    // 2D results are named in wafer axes: a cross-section has x and z, never y (ADR 1)
    let names: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, vec!["V", "|E|", "Ex", "Ez"]);
    for f in &s.fields {
        assert_eq!(f.values.len(), r.dofs, "one value per dof: {}", f.name);
        assert!(f.values.iter().all(|x| x.is_finite()));
    }
    // the per-cell fields: indexed `cell * nloc + a`, node average equals the plot field of the same name (ADR 9)
    let cf = s.cell_fields.as_ref().expect("cell_fields");
    assert_eq!(cf.nloc, mesh.nloc(spec.order));
    assert_eq!(cf.get("Ex").unwrap().re.len(), mesh.n_cells() * cf.nloc);
    assert_eq!(cf.node_average(&mesh, "Ex").unwrap(), s.field("Ex").unwrap().values);
}

#[test]
fn solve_refuses_orders_other_than_one_and_two_and_nan_input() {
    let (mesh, mut spec) = electrostatics_case();
    let mut state = SolveState::default();
    spec.order = 3;
    assert_eq!(nyx_engine::legacy::solve(&mesh, &spec, &mut state).unwrap_err(), "element order must be 1 or 2");
    let (mut mesh, spec) = electrostatics_case();
    mesh.points[0] = f64::NAN;
    assert!(nyx_engine::legacy::solve(&mesh, &spec, &mut state).is_err());
}

#[test]
fn a_cancelled_solve_returns_the_cancelled_string() {
    // callers (the app, the MCP session) compare the error text with `progress::CANCELLED` (ROUGH_SPOTS R11)
    assert_eq!(femcore::progress::CANCELLED, "cancelled");
    let (mesh, spec) = electrostatics_case();
    let p = femcore::progress::Progress::new();
    p.cancel();
    let mut state = SolveState::default();
    let err = nyx_engine::legacy::solve_with(&mesh, &spec, &mut state, &p).unwrap_err();
    assert_eq!(err, femcore::progress::CANCELLED);
}
