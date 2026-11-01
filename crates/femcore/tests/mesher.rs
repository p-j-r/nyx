//! The 2D mesher and the extrusion to 3D.
use femcore::mesh::SimplexMesh;
use femcore::mesher::*;
use femcore::problem::*;
use std::collections::HashMap;

fn area(m: &SimplexMesh) -> f64 {
    femcore::fem::geometry(m).vol.iter().sum()
}

fn min_angle_deg(m: &SimplexMesh) -> f64 {
    let mut worst: f64 = 180.0;
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let p: Vec<[f64; 2]> = v.iter().map(|&i| [m.point(i as usize)[0], m.point(i as usize)[1]]).collect();
        for k in 0..3 {
            let (a, b, c2) = (p[k], p[(k + 1) % 3], p[(k + 2) % 3]);
            let (u, w) = ([b[0] - a[0], b[1] - a[1]], [c2[0] - a[0], c2[1] - a[1]]);
            let cosv = (u[0] * w[0] + u[1] * w[1]) / ((u[0].hypot(u[1])) * (w[0].hypot(w[1])));
            worst = worst.min(cosv.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    worst
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64, domain: i32) -> Part {
    Part { shape: Shape::Rect { x0, y0, x1, y1 }, domain, h: None }
}

#[test]
fn rectangle_edges_are_numbered_left_bottom_top_right() {
    let m = mesh_2d(&Spec2d { parts: vec![rect(0.0, 0.0, 2.0, 1.0, 1)], h: 0.1, min_angle_deg: 25.0 }).unwrap();
    assert!((area(&m) - 2.0).abs() < 1e-12);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4]);
    assert!(min_angle_deg(&m) > 24.0, "{}", min_angle_deg(&m));
    // facets of boundary 1 (x = 0) have total length 1, boundary 2 (y = 0) length 2
    let len = |lab: i32| -> f64 {
        (0..m.n_facets()).filter(|&f| m.facet_labels[f] == lab).map(|f| femcore::fem::facet_measure(&m, f)).sum()
    };
    assert!((len(1) - 1.0).abs() < 1e-12 && (len(4) - 1.0).abs() < 1e-12);
    assert!((len(2) - 2.0).abs() < 1e-12 && (len(3) - 2.0).abs() < 1e-12);
    for f in 0..m.n_facets() {
        let (p, q) = (m.point(m.facet(f)[0] as usize), m.point(m.facet(f)[1] as usize));
        match m.facet_labels[f] {
            1 => assert!(p[0].abs() < 1e-12 && q[0].abs() < 1e-12),
            4 => assert!((p[0] - 2.0).abs() < 1e-12 && (q[0] - 2.0).abs() < 1e-12),
            _ => {}
        }
    }
}

#[test]
fn hole_and_two_materials() {
    // a plate with a circular hole and a right half of another material
    let n = ((2.0 * std::f64::consts::PI * 0.2 / 0.05f64).ceil() as usize).max(12) as f64;
    let spec = Spec2d {
        parts: vec![rect(0.0, 0.0, 2.0, 1.0, 1), rect(1.0, 0.0, 2.0, 1.0, 2), Part { shape: Shape::Circle { cx: 0.5, cy: 0.5, r: 0.2 }, domain: 0, h: Some(0.05) }],
        h: 0.1,
        min_angle_deg: 25.0,
    };
    let m = mesh_2d(&spec).unwrap();
    let polygon_hole = 0.5 * n * 0.2 * 0.2 * (2.0 * std::f64::consts::PI / n).sin();
    assert!((area(&m) - (2.0 - polygon_hole)).abs() < 1e-9, "{} vs {}", area(&m), 2.0 - polygon_hole);
    assert_eq!(m.domain_labels(), vec![1, 2]);
    let a2: f64 = femcore::fem::geometry(&m).vol.iter().zip(&m.cell_labels).filter(|(_, &l)| l == 2).map(|(v, _)| v).sum();
    assert!((a2 - 1.0).abs() < 1e-9);
    // the hole is boundary 9 (rect 1: labels 1-4, rect 2: 5-8, circle 9)
    assert!(m.boundary_labels().contains(&9), "{:?}", m.boundary_labels());
    let circ: f64 = (0..m.n_facets()).filter(|&f| m.facet_labels[f] == 9).map(|f| femcore::fem::facet_measure(&m, f)).sum();
    assert!((circ - n * 2.0 * 0.2 * (std::f64::consts::PI / n).sin()).abs() < 1e-9, "{circ}");
    assert!(min_angle_deg(&m) > 20.0, "{}", min_angle_deg(&m));
}

#[test]
fn local_size_refines_only_its_shape() {
    let coarse = mesh_2d(&Spec2d { parts: vec![rect(0.0, 0.0, 2.0, 1.0, 1)], h: 0.2, min_angle_deg: 25.0 }).unwrap();
    let mut spec = Spec2d { parts: vec![rect(0.0, 0.0, 2.0, 1.0, 1), Part { shape: Shape::Rect { x0: 0.0, y0: 0.0, x1: 0.5, y1: 1.0 }, domain: 2, h: Some(0.05) }], h: 0.2, min_angle_deg: 25.0 };
    let fine = mesh_2d(&spec).unwrap();
    let cells_in = |m: &SimplexMesh, lab: i32| m.cell_labels.iter().filter(|&&l| l == lab).count();
    assert!(cells_in(&fine, 2) > 8 * (coarse.n_cells() / 8).max(1), "refined region has {} cells", cells_in(&fine, 2));
    assert!(cells_in(&fine, 1) < 2 * coarse.n_cells());
    spec.parts[1].h = None;
    assert!(mesh_2d(&spec).is_ok());
}

#[test]
fn poisson_on_a_meshed_disk_matches_the_exact_solution() {
    // -div(eps0 grad V) = eps0 on a disk of radius R, V = 0 on the rim: V = (R^2 - r^2)/4
    let r = 0.5;
    let m = mesh_2d(&Spec2d { parts: vec![Part { shape: Shape::Circle { cx: 0.0, cy: 0.0, r }, domain: 1, h: None }], h: 0.05, min_angle_deg: 25.0 }).unwrap();
    let mut s = ProblemSpec::new(Physics::Electrostatics);
    s.materials.default.insert("rho_v".into(), femcore::scalar::EPS0);
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Ground }];
    let res = femcore::scalar::solve_electrostatics(&m, &s).unwrap();
    let vmax = res.scalar("V max").unwrap();
    assert!((vmax / (r * r / 4.0) - 1.0).abs() < 5e-3, "{vmax} vs {}", r * r / 4.0);
}

#[test]
fn polygon_with_a_notch() {
    // an L shape: area 3
    let l = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]];
    let m = mesh_2d(&Spec2d { parts: vec![Part { shape: Shape::Polygon(l), domain: 1, h: None }], h: 0.1, min_angle_deg: 25.0 }).unwrap();
    assert!((area(&m) - 3.0).abs() < 1e-9);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4, 5, 6]);
}

fn assert_conforming(m: &SimplexMesh) {
    let mut faces: HashMap<[u32; 3], usize> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            f.sort_unstable();
            *faces.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    assert!(faces.values().all(|&n| n <= 2));
    let mut boundary: Vec<[u32; 3]> = faces.iter().filter(|(_, &n)| n == 1).map(|(f, _)| *f).collect();
    let mut given: Vec<[u32; 3]> = (0..m.n_facets()).map(|f| {
        let mut v = m.facet(f).to_vec();
        v.sort_unstable();
        [v[0], v[1], v[2]]
    }).collect();
    boundary.sort();
    given.sort();
    assert_eq!(boundary, given);
}

#[test]
fn extrusion_is_conforming_and_solves_exactly() {
    let m2 = mesh_2d(&Spec2d { parts: vec![rect(0.0, 0.0, 1.0, 0.4, 1)], h: 0.1, min_angle_deg: 25.0 }).unwrap();
    let layers = vec![Layer { thickness: 0.2, n: 2, domain_offset: 0 }, Layer { thickness: 0.1, n: 1, domain_offset: 10 }];
    let m3 = extrude(&m2, &layers, 0.0).unwrap();
    assert_eq!(m3.n_cells(), m2.n_cells() * 3 * 3);
    assert!((area(&m3) - 1.0 * 0.4 * 0.3).abs() < 1e-12);
    assert_conforming(&m3);
    assert_eq!(m3.domain_labels(), vec![1, 11]);
    assert_eq!(m3.boundary_labels(), vec![1, 2, 3, 4, 5, 6]);
    // potential 1 on x = 0 (boundary 1), ground on x = 1 (boundary 4): V = 1 - x exactly
    let mut s = ProblemSpec::new(Physics::Electrostatics);
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    let r = femcore::scalar::solve_electrostatics(&m3, &s).unwrap();
    let v = &r.solutions[0].field("V").unwrap().values;
    let x = m3.dof_coordinates(2);
    let err = (0..v.len()).map(|i| (v[i] - (1.0 - x[3 * i])).abs()).fold(0.0, f64::max);
    // the shared corner dofs of the two Dirichlet boundaries are not present (they are disjoint), so exact
    assert!(err < 1e-10, "{err}");
}

#[test]
fn spec_round_trips_through_json() {
    let spec = Spec2d { parts: vec![rect(0.0, 0.0, 1.0, 1.0, 1), Part { shape: Shape::Polygon(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), domain: 2, h: Some(0.1) }], h: 0.2, min_angle_deg: 25.0 };
    let s = serde_json::to_string(&spec).unwrap();
    assert_eq!(serde_json::from_str::<Spec2d>(&s).unwrap(), spec);
}

#[test]
fn plate_with_a_hole_has_the_textbook_stress_concentration() {
    // 0.4 x 0.2 plate, hole radius 0.04, 1 MPa gross tension: net-section stress 1e6 * 0.2 / 0.12 = 1.667 MPa and
    // Kt (net, d/w = 0.4) about 2.16 (Peterson), so sigma_max about 3.6 MPa
    let spec = Spec2d {
        parts: vec![rect(0.0, 0.0, 0.4, 0.2, 1), Part { shape: Shape::Circle { cx: 0.2, cy: 0.1, r: 0.04 }, domain: 0, h: Some(0.004) }],
        h: 0.012,
        min_angle_deg: 25.0,
    };
    let m = mesh_2d(&spec).unwrap();
    let mut s = ProblemSpec::new(Physics::Solid);
    s.plane = Plane::Stress;
    s.materials.default.insert("E".into(), 200e9);
    s.materials.default.insert("nu".into(), 0.3);
    s.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![4], kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = femcore::solve::solve(&m, &s, &mut SolveState::default()).unwrap();
    let sxx_max = r.solutions[0].field("sxx").unwrap().values.iter().cloned().fold(0.0, f64::max);
    assert!(sxx_max > 3.4e6 && sxx_max < 3.9e6, "sigma_xx max {sxx_max}");
}
