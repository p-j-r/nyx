//! Physics on meshes made from geom models: what the solvers promise on curved boundaries and holes, and the sweep of a section along y.
mod common;
use femcore::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::mesher::{extrude_y, Layer};
use femcore::problem::*;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use std::collections::HashMap;

#[test]
fn poisson_on_a_meshed_disk_matches_the_exact_solution() {
    // -div(eps0 grad V) = eps0 on a disk of radius R, V = 0 on the rim: V = (R^2 - r^2)/4
    let r = 0.5;
    let mut m = Model::new_2d();
    let d = m.circle("disk", (0.0, 0.0), r);
    m.overlay(&[("disk", d)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.05)).unwrap();
    let mut s = ProblemSpec::new(Physics::Electrostatics);
    s.materials.default.insert("rho_v".into(), femcore::scalar::EPS0);
    s.bcs = vec![Bc { on: None, labels: g.select("outer").unwrap(), kind: BcKind::Ground }];
    let res = femcore::scalar::solve_electrostatics(&g.mesh, &s).unwrap();
    let vmax = res.scalar("V max").unwrap();
    assert!((vmax / (r * r / 4.0) - 1.0).abs() < 5e-3, "{vmax} vs {}", r * r / 4.0);
}

#[test]
fn plate_with_a_hole_has_the_textbook_stress_concentration() {
    // 0.4 x 0.2 plate in the x-z plane, hole radius 0.04, 1 MPa gross tension along x: net-section stress 1e6 * 0.2 / 0.12 = 1.667 MPa and
    // Kt (net, d/w = 0.4) about 2.16 (Peterson), so sigma_max about 3.6 MPa
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (0.4, 0.2));
    let hole = m.circle("hole", (0.2, 0.1), 0.04);
    m.overlay_full(&[("plate", plate, RegionKind::Material, None), ("hole", hole, RegionKind::Hole, None)], &[]);
    m.size("size_hole", Selector::parse("in_box([0.15, 0, 0.05], [0.25, 0, 0.15])").unwrap(), Some(0.004.into()), None, None);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.012)).unwrap();
    let mut s = ProblemSpec::new(Physics::Solid);
    s.plane = Plane::Stress;
    s.materials.default.insert("E".into(), 200e9);
    s.materials.default.insert("nu".into(), 0.3);
    // wafer components (x, y, z): the section is the x-z plane
    s.bcs = vec![
        Bc { on: None, labels: g.select("plane(x = 0)").unwrap(), kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: g.select("plane(z = 0)").unwrap(), kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { on: None, labels: g.select("plane(x = 0.4)").unwrap(), kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = femcore::solve::solve(&g.mesh, &s, &mut SolveState::default()).unwrap();
    let sxx_max = r.solutions[0].field("sxx").unwrap().values.iter().cloned().fold(0.0, f64::max);
    assert!(sxx_max > 3.4e6 && sxx_max < 3.9e6, "sigma_xx max {sxx_max}");
}

fn assert_conforming(m: &femcore::mesh::SimplexMesh) {
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
    let mut given: Vec<[u32; 3]> = (0..m.n_facets())
        .map(|f| {
            let mut v = m.facet(f).to_vec();
            v.sort_unstable();
            [v[0], v[1], v[2]]
        })
        .collect();
    boundary.sort();
    given.sort();
    assert_eq!(boundary, given);
}

#[test]
fn a_section_swept_along_y_is_conforming_and_solves_exactly() {
    let p = femcore::geomdrive::plate(1.0, 0.4, 0.1).unwrap();
    let layers = vec![Layer { thickness: 0.2, n: 2, domain_offset: 0 }, Layer { thickness: 0.1, n: 1, domain_offset: 10 }];
    let m3 = extrude_y(&p.mesh, &layers, 0.0).unwrap();
    assert_eq!(m3.n_cells(), p.mesh.n_cells() * 3 * 3);
    let vol: f64 = femcore::fem::geometry(&m3).vol.iter().sum();
    assert!((vol - 1.0 * 0.4 * 0.3).abs() < 1e-12);
    assert_conforming(&m3);
    assert_eq!(m3.domain_labels(), vec![1, 11]);
    // the walls keep the labels of the section, the ends are one more and two more than the largest
    let nb = p.mesh.facet_labels.iter().copied().max().unwrap();
    let mut want = p.mesh.boundary_labels();
    want.extend([nb + 1, nb + 2]);
    assert_eq!(m3.boundary_labels(), want);
    // potential 1 on x = 0, ground on x = 1: V = 1 - x exactly
    let mut s = ProblemSpec::new(Physics::Electrostatics);
    s.bcs = vec![Bc { on: None, labels: p.xmin.clone(), kind: BcKind::Potential(1.0) }, Bc { on: None, labels: p.xmax.clone(), kind: BcKind::Ground }];
    let r = femcore::scalar::solve_electrostatics(&m3, &s).unwrap();
    let v = &r.solutions[0].field("V").unwrap().values;
    let x = m3.dof_coordinates(2);
    let err = (0..v.len()).map(|i| (v[i] - (1.0 - x[3 * i])).abs()).fold(0.0, f64::max);
    assert!(err < 1e-10, "{err}");
}
