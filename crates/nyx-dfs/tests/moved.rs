//! The electrostatics and heat tests that used to sit in the femcore test files (amg_solve, anisotropic, finite_guard, frame, frame_roundtrip,
//! geom_meshes, wafer_names), rewritten against the typed problem. Numbers and tolerances are unchanged.
mod common;
use common::*;
use femcore::material::*;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::props::LinearSolver;
use geom::model::Model;
use nyx_dfs::{Bc, BcKind, Problem, Study};
use std::collections::BTreeSet;

fn set(v: &[&str]) -> BTreeSet<String> {
    v.iter().map(|s| s.to_string()).collect()
}

// ---- from amg_solve.rs
#[test]
fn iterative_poisson_matches() {
    let mesh = block(1.0, 1.0, 1.0, 12, 12, 12);
    let mut s = Problem::new(Study::Electrostatic, 2);
    s.materials.default.insert("eps_r".into(), 4.0);
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    s.solver = LinearSolver::Direct;
    let d = try_solve(&mesh, &s).unwrap();
    s.solver = LinearSolver::Iterative;
    let it = try_solve(&mesh, &s).unwrap();
    assert!(maxdiff(&d.solutions[0].field("V").unwrap().values, &it.solutions[0].field("V").unwrap().values) < 1e-7);
    eprintln!("{}", it.log[0]);
}

#[test]
fn heat_stationary_iterative() {
    let mesh = block(1.0, 0.5, 0.5, 14, 7, 7);
    let mut s = Problem::new(Study::HeatStationary, 2);
    s.materials.default.insert("k".into(), 3.0);
    s.materials.default.insert("Q".into(), 100.0);
    s.bcs = vec![Bc { labels: vec![1, 6], kind: BcKind::Temperature(300.0) }];
    s.solver = LinearSolver::Direct;
    let d = try_solve(&mesh, &s).unwrap();
    s.solver = LinearSolver::Iterative;
    let it = try_solve(&mesh, &s).unwrap();
    assert!(maxdiff(&d.solutions[0].field("T").unwrap().values, &it.solutions[0].field("T").unwrap().values) < 1e-6);
}

#[test]
fn heat_transient_iterative_matches_direct() {
    // large enough to leave the dense eigendecomposition branch and use the BDF2 time stepping
    let mesh = block(1.0, 0.5, 0.5, 12, 6, 6);
    let mut s = Problem::new(Study::HeatStationary, 2);
    for (k, v) in [("k", 3.0), ("rho", 1000.0), ("Cp", 500.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(400.0) },
        Bc { labels: vec![6], kind: BcKind::Temperature(300.0) },
    ];
    s.study = Study::HeatTransient { t_end: 2000.0, n_out: 4, t0: "300".into() };
    s.solver = LinearSolver::Direct;
    let d = try_solve(&mesh, &s).unwrap();
    s.solver = LinearSolver::Iterative;
    let it = try_solve(&mesh, &s).unwrap();
    assert!(it.log.iter().any(|l| l.contains("BDF2")), "{:?}", it.log);
    assert!(it.log.iter().any(|l| l.contains("AMG-CG")), "{:?}", it.log);
    assert_eq!(d.solutions.len(), it.solutions.len());
    for (a, b) in d.solutions.iter().zip(&it.solutions) {
        assert!(maxdiff(&a.field("T").unwrap().values, &b.field("T").unwrap().values) < 1e-5, "{}", a.label);
    }
}

// ---- from anisotropic.rs
#[test]
fn anisotropic_conductivity_is_exact_for_linear_temperature() {
    // any constant tensor: T = T0 + g.x is an exact solution with matching Dirichlet data
    let mesh = rect(1.0, 0.5, 6, 4);
    let k = Tensor::Symmetric { xx: 5.0, yy: 1.0, zz: 1.0, xy: 0.0, yz: 0.0, xz: 1.5 };
    let mut spec = Problem::new(Study::HeatStationary, 2);
    spec.materials
        .aniso
        .insert(1, Anisotropy { conductivity: Some(k.clone()), euler_deg: [0.0; 3], ..Default::default() });
    spec.materials.default.insert("k".into(), 1.0);
    // linear field T = 300 + 20 x + 10 z imposed on all four sides through two values per side
    // (Dirichlet only supports constants, so use a 1D gradient in x and insulated sides: not exact
    // for a tensor with kxz != 0) -> instead check the flux of the 1D x-solution directly
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Temperature(320.0) },
    ];
    let r = nyx_dfs::solve_heat(&mesh, &spec).unwrap();
    let t = &r.solutions[0].field("T").unwrap().values;
    let (tmin, tmax) = (t.iter().cloned().fold(f64::MAX, f64::min), t.iter().cloned().fold(f64::MIN, f64::max));
    assert!((tmin - 300.0).abs() < 1e-9 && (tmax - 320.0).abs() < 1e-9);
    // with kxz != 0 the insulated sides bend the field, so the y flux is not zero in the interior,
    // while for a diagonal tensor the solution is exactly linear and qx = -kxx dT/dx
    let mut diag = spec.clone();
    diag.materials.aniso.insert(
        1,
        Anisotropy {
            conductivity: Some(Tensor::Diagonal { x: 5.0, y: 1.0, z: 1.0 }),
            euler_deg: [0.0; 3],
            ..Default::default()
        },
    );
    let r = nyx_dfs::solve_heat(&mesh, &diag).unwrap();
    let s = &r.solutions[0];
    let x = mesh.dof_coordinates(2);
    let exact_t: Vec<f64> = (0..x.len() / 2).map(|i| 300.0 + 20.0 * x[2 * i]).collect();
    assert!(maxdiff(&s.field("T").unwrap().values, &exact_t) < 1e-10);
    let qx = &s.flat("qx").unwrap();
    assert!(qx.iter().all(|q| (q + 5.0 * 20.0).abs() < 1e-8), "qx should be -100 everywhere");
    assert!(s.flat("qz").unwrap().iter().all(|q| q.abs() < 1e-8));
    // rotating the diagonal tensor by 90 degrees about z swaps kxx and kyy: the x-flux drops to -1*20
    let mut rot = diag.clone();
    rot.materials.aniso.get_mut(&1).unwrap().euler_deg = [90.0, 0.0, 0.0];
    let r = nyx_dfs::solve_heat(&mesh, &rot).unwrap();
    assert!(r.solutions[0].flat("qx").unwrap().iter().all(|q| (q + 20.0).abs() < 1e-8));
}

#[test]
fn anisotropic_permittivity_sets_the_capacitance() {
    // plates at x = 0 and x = L, insulated sides: C = eps0 eps_xx H / L per metre
    let mesh = rect(1.0, 0.5, 8, 4);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.aniso.insert(
        1,
        Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 3.0 }), ..Default::default() },
    );
    spec.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    let r = nyx_dfs::solve_electrostatics(&mesh, &spec).unwrap();
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (nyx_dfs::EPS0 * 7.0 * 0.5) - 1.0).abs() < 1e-9, "{c}");
    // rotated by 90 degrees the roles of x and y swap
    spec.materials.aniso.get_mut(&1).unwrap().euler_deg = [90.0, 0.0, 0.0];
    let r = nyx_dfs::solve_electrostatics(&mesh, &spec).unwrap();
    assert!((r.scalar("capacitance").unwrap() / (nyx_dfs::EPS0 * 3.0 * 0.5) - 1.0).abs() < 1e-9);
}

// ---- from finite_guard.rs
#[test]
fn electrostatics_and_heat_refuse_nan_materials_too() {
    let m = rect(1.0, 1.0, 6, 6);
    let mut es = Problem::new(Study::Electrostatic, 2);
    es.materials.default.insert("eps_r".into(), f64::NAN);
    es.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![3], kind: BcKind::Ground }];
    assert!(try_solve(&m, &es).unwrap_err().contains("'eps_r'"));
    let mut heat = Problem::new(Study::HeatStationary, 2);
    heat.materials.default.insert("k".into(), f64::NAN);
    heat.bcs = vec![Bc { labels: vec![1], kind: BcKind::Temperature(1.0) }];
    assert!(try_solve(&m, &heat).unwrap_err().contains("'k'"));
}

// ---- from frame.rs
#[test]
fn a_permittivity_tensor_in_the_wafer_frame_acts_on_the_vertical_axis_of_the_x_z_plane() {
    // plates on the bottom (z = 0) and top (z = H) of an x-z section of length L: C = eps0 eps_zz L / H per metre
    let mesh = rect(1.0, 0.5, 8, 4);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.aniso.insert(
        1,
        Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }), ..Default::default() },
    );
    spec.bcs = vec![Bc { labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { labels: vec![3], kind: BcKind::Ground }];
    let r = solve(&mesh, &spec.clone());
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (nyx_dfs::EPS0 * 5.0 * 1.0 / 0.5) - 1.0).abs() < 1e-9, "{c}");
    // the field of a uniform vertical field is called Ez
    assert!(r.solutions[0].flat("Ez").is_some() && r.solutions[0].flat("Ey").is_none());
    // plates on the sides see eps_xx
    spec.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    let r = try_solve(&mesh, &spec).unwrap();
    assert!((r.scalar("capacitance").unwrap() / (nyx_dfs::EPS0 * 7.0 * 0.5 / 1.0) - 1.0).abs() < 1e-9);
}

// ---- from frame_roundtrip.rs
#[test]
fn electrostatics_points_from_the_high_plate_along_z() {
    let mesh = rect(1.0, 0.5, 8, 4);
    let mut s = Problem::new(Study::Electrostatic, 2);
    s.materials.default.insert("eps_r".into(), 1.0);
    s.bcs = vec![Bc { labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { labels: vec![3], kind: BcKind::Ground }];
    let r = solve(&mesh, &s);
    assert_eq!(names(&r), set(&["Ex", "Ez", "V", "|E|"]));
    assert!((mean(&r, "Ez") - 2.0).abs() < 1e-6, "E = -grad V = +2 V/m along +z, got {}", mean(&r, "Ez"));
    assert!(mean(&r, "Ex").abs() < 1e-9);
}

#[test]
fn heat_flows_from_the_hot_plate_along_z() {
    let mesh = rect(1.0, 0.5, 8, 4);
    let mut s = Problem::new(Study::HeatStationary, 2);
    s.materials.default.insert("k".into(), 2.0);
    s.bcs = vec![
        Bc { labels: vec![2], kind: BcKind::Temperature(400.0) },
        Bc { labels: vec![3], kind: BcKind::Temperature(300.0) },
    ];
    let r = solve(&mesh, &s);
    assert_eq!(names(&r), set(&["T", "qx", "qz", "|q|"]));
    // q = -k grad T = 2 * 100 / 0.5 along +z
    assert!((mean(&r, "qz") - 400.0).abs() < 1e-6 * 400.0, "{}", mean(&r, "qz"));
}

// ---- from geom_meshes.rs
#[test]
fn poisson_on_a_meshed_disk_matches_the_exact_solution() {
    // -div(eps0 grad V) = eps0 on a disk of radius R, V = 0 on the rim: V = (R^2 - r^2)/4
    let r = 0.5;
    let mut m = Model::new_2d();
    let d = m.circle("disk", (0.0, 0.0), r);
    m.overlay(&[("disk", d)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.05)).unwrap();
    let mut s = Problem::new(Study::Electrostatic, 2);
    s.materials.default.insert("rho_v".into(), nyx_dfs::EPS0);
    s.bcs = vec![Bc { labels: g.select("outer").unwrap(), kind: BcKind::Ground }];
    let res = nyx_dfs::solve_electrostatics(&g.mesh, &s).unwrap();
    let vmax = res.scalar("V max").unwrap();
    assert!((vmax / (r * r / 4.0) - 1.0).abs() < 5e-3, "{vmax} vs {}", r * r / 4.0);
}

#[test]
fn a_section_swept_along_y_solves_exactly() {
    use femcore::meshing::mesher::{extrude_y, Layer};
    let p = femcore::meshing::geomdrive::plate(1.0, 0.4, 0.1).unwrap();
    let layers =
        vec![Layer { thickness: 0.2, n: 2, domain_offset: 0 }, Layer { thickness: 0.1, n: 1, domain_offset: 10 }];
    let m3 = extrude_y(&p.mesh, &layers, 0.0).unwrap();
    // potential 1 on x = 0, ground on x = 1: V = 1 - x exactly
    let mut s = Problem::new(Study::Electrostatic, 2);
    s.bcs = vec![
        Bc { labels: p.xmin.clone(), kind: BcKind::Potential(1.0) },
        Bc { labels: p.xmax.clone(), kind: BcKind::Ground },
    ];
    let r = nyx_dfs::solve_electrostatics(&m3, &s).unwrap();
    let v = &r.solutions[0].field("V").unwrap().values;
    let x = m3.dof_coordinates(2);
    let err = (0..v.len()).map(|i| (v[i] - (1.0 - x[3 * i])).abs()).fold(0.0, f64::max);
    assert!(err < 1e-10, "{err}");
}

// ---- from wafer_names.rs
#[test]
fn an_initial_temperature_on_a_cross_section_is_written_in_x_and_z() {
    // 1 x 0.2 section (x lateral, z vertical), T = 0 on the bottom and the top: T0 = sin(pi z / 0.2) decays as exp(-k (pi / 0.2)^2 t)
    let mesh = rect(1.0, 0.2, 20, 4);
    let mut s = Problem::new(Study::HeatStationary, 2);
    s.materials.default.insert("k".into(), 2.0);
    s.materials.default.insert("rho".into(), 1.0);
    s.materials.default.insert("Cp".into(), 1.0);
    s.bcs = vec![Bc { labels: vec![2, 3], kind: BcKind::Temperature(0.0) }];
    s.study = Study::HeatTransient { t_end: 0.001, n_out: 2, t0: "sin(pi*z/0.2)".into() };
    let r = try_solve(&mesh, &s).unwrap();
    let t = &r.solutions[2].field("T").unwrap().values;
    let xz = mesh.dof_coordinates(2);
    let decay = (-2.0 * (std::f64::consts::PI / 0.2).powi(2) * 0.001).exp();
    let exact: Vec<f64> = (0..t.len()).map(|i| (std::f64::consts::PI * xz[2 * i + 1] / 0.2).sin() * decay).collect();
    assert!(maxdiff(t, &exact) < 5e-3, "{}", maxdiff(t, &exact));
}

// ---- from wafer_names.rs
#[test]
fn scalar_kernels_name_their_components_x_and_z() {
    let mesh = rect(1.0, 0.5, 6, 3);
    let mut es = Problem::new(Study::Electrostatic, 2);
    es.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    assert_eq!(names(&solve(&mesh, &es)), set(&["V", "|E|", "Ex", "Ez"]));
    let mut heat = Problem::new(Study::HeatStationary, 2);
    heat.materials.default.insert("k".into(), 1.0);
    heat.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(1.0) },
        Bc { labels: vec![4], kind: BcKind::Temperature(0.0) },
    ];
    assert_eq!(names(&solve(&mesh, &heat)), set(&["T", "|q|", "qx", "qz"]));
}
