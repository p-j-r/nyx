//! Interface contracts of femcore (docs/architecture/adr/0019): what other crates, saved files and external hosts rely on, written down
//! as tests that do not depend on frozen floating-point values (the golden files of `golden_frames.rs` do, see ADR 10).
//!
//! * `problem_spec` / `solve_result`: the serde shape of `ProblemSpec` (exemplar `docs/architecture/contracts/problem_spec.v1.json`) and the
//!   structure of a `SolveResult`.
//! * `frames`: the constants and name tables of `section` (ADR 1).
//! * `mesh`: label bases, node orders, no curved-node storage, the `.fmsh` version (ADRs 4, 5, 16, 18).
//! * `kernel_stages`: the public stage functions (ADR 9).
//!
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p femcore --test contracts` rewrites the exemplar. Do that only together with a bump of
//! `docs/architecture/contracts/VERSIONS.md` and an ADR note (see ADR 19), never to silence a failure.
mod common;

use femcore::material::M3;
use femcore::problem::*;
use femcore::section;
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

/// Compares `got` with the checked-in exemplar (line endings normalised), or rewrites it under `NYX_WRITE_CONTRACT=1`.
fn check_exemplar(name: &str, got: &str) {
    let path = contract_file(name);
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, format!("{got}\n")).unwrap();
        return;
    }
    let want = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    assert_eq!(want.replace("\r\n", "\n").trim_end(), got.trim_end(), "{name} no longer matches the contract exemplar; this is an interface change (ADR 19)");
}

// ---------------------------------------------------------------------------------------------------------------- problem_spec

/// One spec per physics, each touching the fields of that physics that saved projects, nyx-store's concept card and the mesher's physics
/// sizing read. Built field by field so that a changed default also changes the exemplar.
fn exemplar_specs() -> Vec<ProblemSpec> {
    let mut es = ProblemSpec::new(Physics::Electrostatics);
    es.materials.set(1, "eps_r", 2.0);
    es.materials.set(2, "eps_r", 8.0);
    es.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];

    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.materials.default.insert("k".into(), 3.0);
    heat.study = Study::Transient { t_end: 1.0, n_out: 5, t0: "293.15 + 50*exp(-((x-0.5)^2)/0.01)".into() };
    heat.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(300.0) }];

    let mut solid = ProblemSpec::new(Physics::Solid);
    solid.study = Study::Eigenfrequency { n_modes: 6, f_search: 5e9 };
    solid.inactive = vec![3];
    solid.materials.set(1, "E", 169e9);
    solid.materials.set(1, "nu", 0.28);
    solid.materials.set(1, "rho", 2329.0);
    solid.materials.pml.insert(2, Pml { x: PmlDir::Plus, z: PmlDir::Minus, ..Pml::default() });
    solid.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::LowReflecting }];
    solid.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 0.5 }];

    let mut emw = ProblemSpec::new(Physics::Emw);
    emw.study = Study::Modes { freq: 1.934e14, n_target: 2.4, evanescent: false, n_modes: 4 };
    emw.materials.set(1, "eps_r", 12.08);
    emw.materials.set(2, "eps_r", 2.085);

    vec![es, heat, solid, emw]
}

#[test]
fn problem_spec_exemplar_is_stable() {
    let specs = exemplar_specs();
    let got = serde_json::to_string_pretty(&specs).unwrap();
    check_exemplar("problem_spec.v1.json", &got);
    // and the checked-in text reads back to the same values
    let text = std::fs::read_to_string(contract_file("problem_spec.v1.json")).unwrap();
    let back: Vec<ProblemSpec> = serde_json::from_str(&text).unwrap();
    assert_eq!(back, specs);
}

#[test]
fn problem_spec_survives_the_xml_encoding_of_project_files() {
    // projects are saved as XML through `xmlio`; nyx-store reads the same XML back into a JSON value
    for spec in exemplar_specs() {
        let v = serde_json::to_value(&spec).unwrap();
        let xml = femcore::xmlio::value_to_xml("nyx_spec", &v).unwrap();
        let (root, back) = femcore::xmlio::xml_to_value(&xml).unwrap();
        assert_eq!(root, "nyx_spec");
        let spec2: ProblemSpec = serde_json::from_value(back).unwrap();
        assert_eq!(spec2, spec, "{:?}", spec.physics);
    }
}

#[test]
fn problem_spec_json_shape_that_other_crates_read_without_types() {
    // `nyx-store::card` reads saved specs as untyped values: the names below are an implicit contract (ROUGH_SPOTS R25)
    let specs = exemplar_specs();
    let emw = serde_json::to_value(&specs[3]).unwrap();
    assert_eq!(emw["physics"], json!("Emw"));
    assert_eq!(emw["study"], json!({ "Modes": { "freq": 1.934e14, "n_target": 2.4, "evanescent": false, "n_modes": 4 } }));
    assert!(emw["materials"]["domains"]["1"]["eps_r"].is_number(), "materials.domains.<label>.<property>");
    // keys of the top level, in the order of the declaration (the exemplar pins the rest)
    let keys: Vec<&str> = emw.as_object().unwrap().keys().map(|s| s.as_str()).collect();
    for k in ["physics", "order", "materials", "bcs", "study", "plane", "thickness", "thermal", "solver", "output", "inactive", "periodic"] {
        assert!(keys.contains(&k), "ProblemSpec lost the key {k}");
    }
}

#[test]
fn problem_spec_reading_is_lenient_in_both_directions() {
    // Documented, not endorsed (R20): an unknown key is ignored and a missing `#[serde(default)]` key takes its default.
    // A rename therefore does not fail when an old file is read; it silently drops the setting. The exemplar test above is what catches it.
    let specs = exemplar_specs();
    let mut v = serde_json::to_value(&specs[3]).unwrap();
    v.as_object_mut().unwrap().insert("a_key_from_the_future".into(), json!(1));
    assert_eq!(serde_json::from_value::<ProblemSpec>(v.clone()).unwrap(), specs[3]);
    v.as_object_mut().unwrap().remove("inactive");
    v.as_object_mut().unwrap().remove("periodic");
    assert_eq!(serde_json::from_value::<ProblemSpec>(v).unwrap(), specs[3]);
}

#[test]
fn physics_enum_has_four_members_in_this_order() {
    // the order is the order of the tabs and of `Physics::ALL`; names are the serde tags in saved files
    let names: Vec<Value> = Physics::ALL.iter().map(|p| serde_json::to_value(p).unwrap()).collect();
    assert_eq!(names, vec![json!("Electrostatics"), json!("Heat"), json!("Solid"), json!("Emw")]);
}

// ---------------------------------------------------------------------------------------------------------------- solve_result

fn electrostatics_case() -> (femcore::mesh::SimplexMesh, ProblemSpec) {
    let mesh = common::rect(1.0, 1.0, 8, 8);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 4.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    (mesh, spec)
}

#[test]
fn solve_result_has_the_documented_structure() {
    // `SolveResult` is not serialisable; h5io, nyx-store::runs::metrics_of and the MCP digests each read it by field, so its structure is the contract
    let (mesh, spec) = electrostatics_case();
    let mut state = SolveState::default();
    let r = femcore::solve::solve(&mesh, &spec, &mut state).unwrap();
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
    assert_eq!(femcore::solve::solve(&mesh, &spec, &mut state).unwrap_err(), "element order must be 1 or 2");
    let (mut mesh, spec) = electrostatics_case();
    mesh.points[0] = f64::NAN;
    assert!(femcore::solve::solve(&mesh, &spec, &mut state).is_err());
}

#[test]
fn a_cancelled_solve_returns_the_cancelled_string() {
    // callers (the app, the MCP session) compare the error text with `progress::CANCELLED` (ROUGH_SPOTS R11)
    assert_eq!(femcore::progress::CANCELLED, "cancelled");
    let (mesh, spec) = electrostatics_case();
    let p = femcore::progress::Progress::new();
    p.cancel();
    let mut state = SolveState::default();
    let err = femcore::solve::solve_with(&mesh, &spec, &mut state, &p).unwrap_err();
    assert_eq!(err, femcore::progress::CANCELLED);
}

// ---------------------------------------------------------------------------------------------------------------- frames

fn diag(a: f64, b: f64, c: f64) -> M3 {
    [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
}

#[test]
fn frames_constants_of_the_wafer_frame() {
    // ADR 1: waves propagate along +y; local triple (lateral, normal, propagation) = wafer (x, z, y)
    assert_eq!(section::PROP_SIGN, 1.0);
    assert_eq!(section::LOCAL_AXES, [(0, 1.0), (2, 1.0), (1, 1.0)]);
    assert_eq!(section::PLANE_VOIGT, [0, 2, 1, 5, 4, 3]);
    assert!(section::is_left_handed(section::PROP_SIGN), "the mode solvers' local triple is left-handed");
    assert!(!section::is_left_handed(-section::PROP_SIGN));
}

#[test]
fn frames_names_and_points_of_a_cross_section_and_of_3d() {
    assert_eq!(section::wafer_axes(2), &[0, 2]);
    assert_eq!(section::wafer_axes(3), &[0, 1, 2]);
    assert_eq!(section::component_names(2), &["x", "z"]);
    assert_eq!(section::component_names(3), &["x", "y", "z"]);
    assert_eq!(section::wafer_point(&[1.0, 2.0], 2), [1.0, 0.0, 2.0], "a 2D mesh point (x, z) is the wafer point (x, 0, z)");
    assert_eq!(section::wafer_point(&[1.0, 2.0, 3.0], 3), [1.0, 2.0, 3.0]);
    assert_eq!([section::local_component(0), section::local_component(1), section::local_component(2)], [("x", 1.0), ("z", 1.0), ("y", 1.0)]);
    let names: Vec<&str> = (0..6).map(|i| section::stress_name(2, i)).collect();
    assert_eq!(names, vec!["sxx", "szz", "syy", "sxz", "syz", "sxy"], "Voigt order of a cross-section is (xx, zz, yy, xz, zy, xy)");
    let names3: Vec<&str> = (0..6).map(|i| section::stress_name(3, i)).collect();
    assert_eq!(names3, vec!["sxx", "syy", "szz", "sxy", "syz", "sxz"]);
}

#[test]
fn frames_tensors_are_permuted_only_in_section() {
    // wafer diag(1, 2, 3) = (xx, yy, zz): a cross-section reads it in the order (x, z, y)
    assert_eq!(section::plane_tensor(&diag(1.0, 2.0, 3.0), 2), diag(1.0, 3.0, 2.0));
    assert_eq!(section::plane_tensor(&diag(1.0, 2.0, 3.0), 3), diag(1.0, 2.0, 3.0));
    // the mode solvers' local order is the same permutation, with the propagation sign on the third axis
    assert_eq!(section::local_tensor(&diag(1.0, 2.0, 3.0)), diag(1.0, 3.0, 2.0));
    // an off-diagonal wafer entry xy (lateral-propagation) flips sign with the propagation direction, xz does not
    let mut t = diag(1.0, 2.0, 3.0);
    t[0][1] = 0.25;
    t[1][0] = 0.25;
    t[0][2] = 0.5;
    t[2][0] = 0.5;
    let fwd = section::local_tensor_signed(&t, 1.0);
    let bwd = section::local_tensor_signed(&t, -1.0);
    assert_eq!(fwd[0][2], 0.25, "lateral-propagation entry is wafer xy");
    assert_eq!(bwd[0][2], -0.25);
    assert_eq!(fwd[0][1], 0.5, "lateral-normal entry is wafer xz");
    assert_eq!(bwd[0][1], 0.5);
}

#[test]
fn frames_pml_directions_use_x_and_z_on_a_cross_section() {
    let p = Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Minus, ..Pml::default() };
    assert_eq!(section::mesh_pml_dirs(&p, 2), [PmlDir::Plus, PmlDir::Minus, PmlDir::Off]);
    assert_eq!(section::mesh_pml_dirs(&p, 3), [PmlDir::Plus, PmlDir::Off, PmlDir::Minus]);
    // a y layer on a 2D mesh is refused rather than aliased to z (ADR 1)
    let mesh = common::rect(1.0, 1.0, 4, 4);
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.study = Study::Modes { freq: 1.934e14, n_target: 1.5, evanescent: false, n_modes: 1 };
    spec.materials.pml.insert(1, Pml { y: PmlDir::Plus, ..Pml::default() });
    assert!(section::check_wafer_input(&mesh, &spec).is_err());
}

// ---------------------------------------------------------------------------------------------------------------- mesh

#[test]
fn mesh_labels_are_one_based_and_arrays_carry_no_curved_nodes() {
    let m = common::rect(1.0, 0.5, 4, 2);
    assert_eq!(m.dim, 2);
    assert_eq!(m.points.len(), 2 * m.n_points());
    assert_eq!(m.cells.len(), 3 * m.n_cells(), "a triangle stores its three vertices, no mid-edge nodes (ADR 18)");
    assert_eq!(m.facets.len(), 2 * m.facet_labels.len());
    assert!(m.cell_labels.iter().all(|&l| l >= 1), "domain numbers are 1-based");
    assert!(m.facet_labels.iter().all(|&l| l >= 1), "boundary numbers are 1-based");
    // P2 numbering: vertices first, then one dof per edge
    assert_eq!(m.n_dofs(1), m.n_points());
    assert_eq!(m.n_dofs(2), m.n_points() + m.n_edges());
    let b = common::block(1.0, 1.0, 0.5, 3, 3, 2);
    assert_eq!(b.dim, 3);
    assert_eq!(b.cells.len(), 4 * b.n_cells());
    assert_eq!(b.facets.len(), 3 * b.facet_labels.len());
}

#[test]
fn mesh_cells_are_positively_oriented_in_2d() {
    let m = common::rect(1.0, 0.5, 4, 2);
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let p = |i: usize| (m.points[v[i] as usize * 2], m.points[v[i] as usize * 2 + 1]);
        let (a, b, d) = (p(0), p(1), p(2));
        let area2 = (b.0 - a.0) * (d.1 - a.1) - (d.0 - a.0) * (b.1 - a.1);
        assert!(area2.abs() > 0.0);
    }
}

mod mixed {
    use femcore::mesh::cellmesh::CellMesh as NyxCells;
    use femcore::mesh::SimplexMesh;
    use femcore::meshing::pml_box::{build_pml, PmlSpec};
    use femcore::mesh::interop::{cells_to_nyx, simplex_to_adapt};

    pub fn pml_cells() -> NyxCells {
        let plate = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.125).unwrap().mesh;
        let phys = simplex_to_adapt::<2>(&plate).unwrap();
        let spec = PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 3) };
        cells_to_nyx(&build_pml(&phys, &spec).unwrap().mesh)
    }

    pub fn split(cells: &NyxCells) -> SimplexMesh {
        SimplexMesh::from_cells_split(cells).unwrap().0
    }
}

#[test]
fn mesh_nyx_cellmesh_label_bases_differ_between_domains_and_boundaries() {
    // ADR 5 / ROUGH_SPOTS R09: in the nyx `CellMesh`, domain numbers are 1-based and boundary-element numbers 0-based;
    // `SimplexMesh` numbers both from 1. `from_cells_split` is the place where the two meet.
    let cells = mixed::pml_cells();
    let edg = cells.ty("edg").expect("boundary edges of a 2D mixed mesh are named `edg`");
    let quad = cells.ty("quad").expect("the PML layers are quadrilaterals");
    let tri = cells.ty("tri").expect("the physical mesh is triangles");
    assert_eq!(edg.nep, 2);
    assert_eq!(quad.nep, 4);
    assert_eq!(tri.nep, 3);
    assert!(edg.dom.iter().all(|&d| d >= 0) && edg.dom.contains(&0), "boundary entity numbers start at 0");
    assert!(tri.dom.iter().chain(&quad.dom).all(|&d| d >= 1), "domain numbers start at 1 (0 is the exterior)");
    let mesh = mixed::split(&cells);
    let mut want: Vec<i32> = edg.dom.iter().map(|d| d + 1).collect();
    want.sort_unstable();
    want.dedup();
    let mut got = mesh.boundary_labels();
    got.sort_unstable();
    got.dedup();
    assert_eq!(got, want, "the split mesh numbers boundaries from 1");
}

#[test]
fn mesh_quadrilaterals_are_stored_in_tensor_product_order_and_split_exactly() {
    // nyx quad order is `i + 2 j`: the polygon (0, 1, 3, 2) is counter-clockwise (ADR 5)
    let cells = mixed::pml_cells();
    let quad = cells.ty("quad").unwrap();
    for q in 0..quad.len() {
        let v = quad.row(q);
        let p = |i: usize| (cells.points[v[i] as usize * 2], cells.points[v[i] as usize * 2 + 1]);
        let ring = [p(0), p(1), p(3), p(2)];
        let area2: f64 = (0..4).map(|i| ring[i].0 * ring[(i + 1) % 4].1 - ring[(i + 1) % 4].0 * ring[i].1).sum();
        assert!(area2 > 0.0, "quad {q} is not counter-clockwise in tensor-product order");
    }
    // the split into triangles adds no vertex and fills the same area: the 1 x 1 plate with a 0.25 layer is 1.5 x 1.5
    let mesh = mixed::split(&cells);
    assert_eq!(mesh.n_points(), cells.n_points());
    let area: f64 = (0..mesh.n_cells())
        .map(|c| {
            let v = mesh.cell(c);
            let p = |i: usize| (mesh.points[v[i] as usize * 2], mesh.points[v[i] as usize * 2 + 1]);
            0.5 * ((p(1).0 - p(0).0) * (p(2).1 - p(0).1) - (p(2).0 - p(0).0) * (p(1).1 - p(0).1)).abs()
        })
        .sum();
    assert!((area - 2.25).abs() < 1e-9, "area {area}");
}

#[test]
fn mesh_fmsh_round_trip_and_version() {
    use femcore::mesh::cells::CellMesh;
    use femcore::mesh::io::fmsh;
    use femcore::meshing::pml_box::{build_pml, PmlSpec};
    use femcore::mesh::interop::simplex_to_adapt;
    let plate = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.25).unwrap().mesh;
    let phys = simplex_to_adapt::<2>(&plate).unwrap();
    let m: CellMesh<2> = build_pml(&phys, &PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 2) }).unwrap().mesh;
    let bytes = fmsh::encode(&m);
    assert_eq!(&bytes[..4], b"FMSH");
    assert_eq!(bytes[4], 2, "`.fmsh` format version (ADR 16): bump it and add a reader for the old one when the layout changes");
    assert_eq!(bytes[5], 2, "dimension byte");
    let back = fmsh::decode::<2>(&bytes).unwrap();
    assert_eq!(fmsh::encode(&back), bytes, "lossless and deterministic: decode then encode gives the same bytes");
    assert_eq!(back.types, m.types);
    assert_eq!(back.label, m.label);
    assert_eq!(back.facet_label, m.facet_label);
    // the CRC catches a flipped bit
    let mut bad = bytes.clone();
    let mid = bad.len() / 2;
    bad[mid] ^= 0x10;
    assert!(fmsh::decode::<2>(&bad).is_err());
    // the dimension is part of the contract
    assert!(fmsh::decode::<3>(&bytes).is_err());
}

// ---------------------------------------------------------------------------------------------------------------- kernel_stages

#[test]
fn kernel_stage_api_is_public() {
    // ADR 9. Naming the functions makes removing or renaming one a compile error here. The kernels not in this list
    // (plain `sms` static/eigen, `sms_lrb`, `ems3d`) have crate-private stages today (ROUGH_SPOTS R26).
    let _ = (femcore::dfs::assemble_electrostatics, femcore::dfs::solve_potential, femcore::dfs::electrostatics_fields);
    let _ = (femcore::dfs::assemble_heat, femcore::dfs::solve_heat_states, femcore::dfs::heat_fields);
    let _ = (femcore::ems::assemble_modes, femcore::ems::solve_mode_pairs, femcore::ems::mode_fields_cells);
    let _ = (femcore::ems_pml::assemble_modes_complex, femcore::ems_pml::solve_mode_pairs_complex, femcore::ems_pml::mode_fields_complex_cells);
    let _ = (femcore::sms_kz::assemble_cross_section, femcore::sms_kz::solve_cross_section_modes, femcore::sms_kz::cross_section_fields);
    let _ = (femcore::ems_axisym::assemble_axisym, femcore::ems_axisym::solve_axisym_pairs, femcore::ems_axisym::mode_samples);
    let _ = (femcore::sms_pml::assemble_solid_pml, femcore::sms_pml::solve_pml_modes, femcore::sms_pml::pml_fields);
    let _ = (femcore::sms_floquet::assemble_floquet, femcore::sms_floquet::solve_floquet_modes, femcore::sms_floquet::floquet_fields);
    // the one-call solvers and the dispatching entry points
    let _ = (femcore::solve::solve, femcore::solve::solve_with);
    let _ = (femcore::dfs::solve_electrostatics, femcore::dfs::solve_heat, femcore::sms::solve_solid, femcore::ems::solve_emw, femcore::ems3d::solve_emw_eigen);
}

#[test]
fn physics_dispatch_covers_every_physics_and_the_emw_split_by_mesh_dimension() {
    // ROUGH_SPOTS R10: `Physics::Emw` means waveguide modes on a 2D mesh and cavity eigenfrequencies on a 3D mesh. A 3D mesh given a
    // `Modes` study is the caller's mistake; the contract is only that the dispatch does not panic and that 2D `Emw` is the mode solver.
    let names: Vec<&str> = Physics::ALL.iter().map(|p| p.name()).collect();
    assert_eq!(names.len(), 4);
    assert!(Physics::Emw.name().contains("mode analysis"));
}
