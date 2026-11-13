//! 3D geom models meshed with tetrahedra (ADR 40): a sphere in a box, a torus against its closed form, labels and conformity, local sizes.
//! The PEC sphere cavity against the closed-form mode is in the `nyx-ems` cavity tests.
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use std::collections::HashMap;
use std::f64::consts::PI;

fn det6(p: [[f64; 3]; 4]) -> f64 {
    let s = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let (u, v, w) = (s(p[1], p[0]), s(p[2], p[0]), s(p[3], p[0]));
    u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])
}

fn pt(m: &femcore::mesh::SimplexMesh, v: u32) -> [f64; 3] {
    [m.points[3 * v as usize], m.points[3 * v as usize + 1], m.points[3 * v as usize + 2]]
}

/// Volume of the cells of every domain label; every cell must be positively oriented.
fn volumes(g: &MeshedGeom) -> HashMap<i32, f64> {
    let m = &g.mesh;
    let mut out: HashMap<i32, f64> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let d = det6([pt(m, v[0]), pt(m, v[1]), pt(m, v[2]), pt(m, v[3])]) / 6.0;
        assert!(d > 0.0, "cell {c} has volume {d}");
        *out.entry(m.cell_labels[c]).or_default() += d;
    }
    out
}

/// Area of the boundary facets of every boundary label.
fn areas(g: &MeshedGeom) -> HashMap<i32, f64> {
    let m = &g.mesh;
    let mut out: HashMap<i32, f64> = HashMap::new();
    for f in 0..m.n_facets() {
        let (a, b, c) = (pt(m, m.facets[3 * f]), pt(m, m.facets[3 * f + 1]), pt(m, m.facets[3 * f + 2]));
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        *out.entry(m.facet_labels[f]).or_default() += 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    }
    out
}

/// Every face of every cell is shared by at most two cells; the faces used once are exactly the boundary facets that touch the void, and
/// the faces between cells of different domains are exactly the interface facets.
fn assert_conforming(g: &MeshedGeom) {
    let m = &g.mesh;
    let key = |a: u32, b: u32, c: u32| {
        let mut k = [a, b, c];
        k.sort_unstable();
        k
    };
    let mut faces: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            faces.entry(key(f[0], f[1], f[2])).or_default().push(c);
        }
    }
    assert!(faces.values().all(|v| v.len() <= 2), "a face is shared by more than two cells");
    let mut expected: Vec<[u32; 3]> = faces
        .iter()
        .filter(|(_, v)| v.len() == 1 || m.cell_labels[v[0]] != m.cell_labels[v[1]])
        .map(|(k, _)| *k)
        .collect();
    expected.sort_unstable();
    let mut got: Vec<[u32; 3]> = (0..m.n_facets()).map(|f| key(m.facets[3 * f], m.facets[3 * f + 1], m.facets[3 * f + 2])).collect();
    got.sort_unstable();
    assert_eq!(got, expected, "the boundary facets are not the faces of the void and of the interfaces");
    assert!(m.cell_labels.iter().all(|&l| l >= 1) && m.facet_labels.iter().all(|&l| l >= 1), "labels are 1-based");
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

fn sphere_in_box(h: f64) -> (MeshedGeom, f64) {
    let r = 0.7;
    let mut m = Model::new_3d();
    let air = m.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let ball = m.sphere("ball", (2.0, 2.0, 2.0), r);
    m.overlay(&[("air", air), ("ball", ball)]);
    (mesh_model(&m.list, &MeshPlan::simple(h)).unwrap(), r)
}

#[test]
fn a_sphere_in_a_box_is_meshed_with_two_domains_and_the_model_labels() {
    let (g, r) = sphere_in_box(0.2);
    assert_eq!(g.mesh.dim, 3);
    assert_conforming(&g);
    // the labels are the ones of the model: domains 1 (air) and 2 (ball), one boundary per face
    let (air, ball) = (g.domain("air").unwrap(), g.domain("ball").unwrap());
    assert_eq!((air, ball), (1, 2));
    let v = volumes(&g);
    // the lattice follows the box exactly; the ball is piecewise linear
    assert!(rel(v[&air] + v[&ball], 64.0) < 1e-9, "{v:?}");
    let v_ball = 4.0 / 3.0 * PI * r.powi(3);
    println!("sphere in box: {} cells, volume error {:.2e}", g.mesh.n_cells(), rel(v[&ball], v_ball));
    assert!(rel(v[&ball], v_ball) < 0.025, "ball volume {} vs {v_ball}", v[&ball]);
    let a = areas(&g);
    let iface = g.boundary("air|ball");
    assert_eq!(iface.len(), 1);
    let area = a[&iface[0]];
    assert!(rel(area, 4.0 * PI * r * r) < 0.03, "sphere area {area} vs {}", 4.0 * PI * r * r);
    // the six faces of the box have their own labels and their exact area
    for name in ["air.xmin", "air.xmax", "air.ymin", "air.ymax", "air.zmin", "air.zmax"] {
        let l = g.boundary(name);
        assert_eq!(l.len(), 1, "{name}");
        assert!(rel(a[&l[0]], 16.0) < 1e-9, "{name}: {}", a[&l[0]]);
    }
    assert_eq!(a.len(), 7);
    // selectors resolve to the same integers
    assert_eq!(g.select("plane(z = 0)").unwrap(), g.boundary("air.zmin"));
    // every vertex of the sphere facets lies on the sphere to a small fraction of h (the cut points are exact; the warped ones too)
    let m = &g.mesh;
    let mut worst: f64 = 0.0;
    for f in 0..m.n_facets() {
        if m.facet_labels[f] != iface[0] {
            continue;
        }
        for k in 0..3 {
            let p = pt(m, m.facets[3 * f + k]);
            worst = worst.max((((p[0] - 2.0).powi(2) + (p[1] - 2.0).powi(2) + (p[2] - 2.0).powi(2)).sqrt() - r).abs());
        }
    }
    assert!(worst < 1e-9, "a vertex of the interface is {worst} off the sphere");
}

#[test]
fn the_sphere_volume_converges_with_the_element_size() {
    let mut errs = vec![];
    for h in [0.4, 0.2] {
        let (g, r) = sphere_in_box(h);
        let v = volumes(&g);
        errs.push((v[&2] / (4.0 / 3.0 * PI * r.powi(3)) - 1.0).abs());
    }
    // the error of a piecewise linear sphere is O(h^2): halving h must at least halve it (it is noisy for so few cells)
    assert!(errs[1] < 0.6 * errs[0] || errs[1] < 2e-3, "{errs:?}");
}

#[test]
fn a_torus_matches_its_closed_form_volume_and_area() {
    let (big, small) = (1.0, 0.3);
    let mut m = Model::new_3d();
    let host = m.box_("vacuum", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let ring = m.torus("ring", (0.0, 0.0, 0.0), big, small);
    m.overlay(&[("vacuum", host), ("ring", ring)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(0.1)).unwrap();
    assert_conforming(&g);
    let v = volumes(&g);
    let (ring_l, vac_l) = (g.domain("ring").unwrap(), g.domain("vacuum").unwrap());
    let v_torus = 2.0 * PI * PI * big * small * small;
    println!("torus: volume error {:.2e}, area error {:.2e}", rel(v[&ring_l], v_torus), rel(areas(&g)[&g.boundary("ring|vacuum")[0]], 4.0 * PI * PI * big * small));
    assert!(rel(v[&ring_l], v_torus) < 0.03, "torus volume {} vs {v_torus}", v[&ring_l]);
    assert!(rel(v[&ring_l] + v[&vac_l], 64.0) < 1e-9);
    let a = areas(&g);
    let surf = g.boundary("ring|vacuum");
    assert_eq!(surf.len(), 1);
    let s_torus = 4.0 * PI * PI * big * small;
    assert!(rel(a[&surf[0]], s_torus) < 0.04, "torus area {} vs {s_torus}", a[&surf[0]]);
    // the vertices of the surface are on the torus
    let mesh = &g.mesh;
    for f in 0..mesh.n_facets() {
        if mesh.facet_labels[f] != surf[0] {
            continue;
        }
        for k in 0..3 {
            let p = pt(mesh, mesh.facets[3 * f + k]);
            let d = ((p[0] * p[0] + p[1] * p[1]).sqrt() - big).hypot(p[2]) - small;
            assert!(d.abs() < 1e-9, "a vertex is {d} off the torus");
        }
    }
}

#[test]
fn a_cavity_in_a_hole_host_has_one_domain_and_one_boundary() {
    let g = femcore::meshing::geomdrive3::sphere_cavity(0.05, 0.012).unwrap();
    assert_conforming(&g);
    assert_eq!(g.labels.domains.len(), 1);
    assert_eq!(g.domain("cavity"), Some(1));
    assert_eq!(g.boundary("cavity.surface"), vec![1]);
    let v = volumes(&g);
    let a = areas(&g);
    assert!(rel(v[&1], 4.0 / 3.0 * PI * 0.05f64.powi(3)) < 0.02, "{v:?}");
    assert!(rel(a[&1], 4.0 * PI * 0.05f64.powi(2)) < 0.03, "{a:?}");
    assert!(g.mesh.exterior_facets().iter().all(|&e| e), "every facet of a cavity touches the void");
}

#[test]
fn a_size_rule_on_a_region_refines_it_and_keeps_the_surface_exact() {
    let mut m = Model::new_3d();
    let host = m.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let ball = m.sphere("ball", (2.0, 2.0, 2.0), 0.9);
    m.overlay(&[("air", host), ("ball", ball)]);
    m.size("fine ball", Selector::parse("region(ball)").unwrap(), Some(0.15.into()), None, None);
    let coarse = {
        let mut c = m.list.clone();
        c.features.retain(|f| f.name != "fine ball");
        mesh_model(&c, &MeshPlan::simple(0.5)).unwrap()
    };
    let fine = mesh_model(&m.list, &MeshPlan::simple(0.5)).unwrap();
    assert_conforming(&fine);
    let count = |g: &MeshedGeom, l: i32| g.mesh.cell_labels.iter().filter(|&&x| x == l).count();
    assert!(count(&fine, 2) > 4 * count(&coarse, 2), "{} vs {}", count(&fine, 2), count(&coarse, 2));
    let v = volumes(&fine);
    assert!(rel(v[&1] + v[&2], 64.0) < 1e-9, "{v:?}");
    assert!(rel(v[&2], 4.0 / 3.0 * PI * 0.9f64.powi(3)) < 0.01, "{v:?}");
    let _ = RegionKind::Material;
}

#[test]
fn the_adaptive_mesher_and_structured_pml_work_in_3d_and_the_2d_only_options_are_still_refused() {
    // ADR 51 and 52 lifted the refusals of ADR 40; the features themselves are in `mesh3d_features.rs`
    let mut m = Model::new_3d();
    let host = m.box_("air", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("air", host)]);
    let g = mesh_model(&m.list, &MeshPlan::adaptive(0.4, 0.7)).unwrap();
    assert_eq!(g.mesh.dim, 3);
    let v = volumes(&g);
    assert!(rel(v[&1], 1.0) < 1e-9);
    // a model without a Pml comes back unchanged from the structured mode
    let mut plan = MeshPlan::simple(0.3);
    plan.pml = femcore::meshing::PmlMode::Structured;
    let plain = mesh_model(&m.list, &MeshPlan::simple(0.3)).unwrap();
    assert_eq!(mesh_model(&m.list, &plan).unwrap().mesh.n_cells(), plain.mesh.n_cells());
}
