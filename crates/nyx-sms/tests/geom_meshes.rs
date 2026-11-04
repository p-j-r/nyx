//! Physics on meshes made from geom models: what the solid solver promises on curved boundaries and holes. (The sweep of a section along y
//! stays in femcore.)
mod common;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_sms::*;

#[test]
fn plate_with_a_hole_has_the_textbook_stress_concentration() {
    // 0.4 x 0.2 plate in the x-z plane, hole radius 0.04, 1 MPa gross tension along x: net-section stress 1e6 * 0.2 / 0.12 = 1.667 MPa and
    // Kt (net, d/w = 0.4) about 2.16 (Peterson), so sigma_max about 3.6 MPa
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (0.4, 0.2));
    let hole = m.circle("hole", (0.2, 0.1), 0.04);
    m.overlay_full(&[("plate", plate, RegionKind::Material, None), ("hole", hole, RegionKind::Hole, None)], &[]);
    m.size(
        "size_hole",
        Selector::parse("in_box([0.15, 0, 0.05], [0.25, 0, 0.15])").unwrap(),
        Some(0.004.into()),
        None,
        None,
    );
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.012)).unwrap();
    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Stress;
    s.materials.default.insert("E".into(), 200e9);
    s.materials.default.insert("nu".into(), 0.3);
    // wafer components (x, y, z): the section is the x-z plane
    s.bcs = vec![
        Bc {
            labels: g.select("plane(x = 0)").unwrap(),
            kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] },
        },
        Bc {
            labels: g.select("plane(z = 0)").unwrap(),
            kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] },
        },
        Bc { labels: g.select("plane(x = 0.4)").unwrap(), kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = common::solve(&g.mesh, &s);
    let sxx_max = r.solutions[0].flat("sxx").unwrap().iter().cloned().fold(0.0, f64::max);
    assert!(sxx_max > 3.4e6 && sxx_max < 3.9e6, "sigma_xx max {sxx_max}");
}
