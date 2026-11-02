//! The unified meshing interface: the simple and the adaptive mesher make meshes of the same model with the same names and labels, and
//! the adaptive one puts its boundary and interface vertices on the exact curves of the geometry.
use femcore::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::{mesh_model, Adaptive, MeshPlan, Mesher};
use femcore::meshers::quality::mesh_quality;
use femcore::problem::*;
use geom::model::{Model, RegionKind};
use std::collections::BTreeSet;

fn disk(r: f64) -> Model {
    let mut m = Model::new_2d();
    let d = m.circle("disk", (0.0, 0.0), r);
    m.overlay(&[("disk", d)]);
    m
}

/// A core disk in a square cladding: one curved interface, one straight outer boundary.
fn core_in_clad(r: f64) -> Model {
    let mut m = Model::new_2d();
    let clad = m.rect("clad", (-1.0, -1.0), (1.0, 1.0));
    let core = m.circle("core", (0.0, 0.0), r);
    m.overlay_full(&[("clad", clad, RegionKind::Material, None), ("core", core, RegionKind::Material, None)], &[]);
    m
}

fn area(mesh: &femcore::mesh::SimplexMesh, label: Option<i32>) -> f64 {
    (0..mesh.n_cells())
        .filter(|&c| label.is_none_or(|l| mesh.cell_labels[c] == l))
        .map(|c| {
            let v = mesh.cell(c);
            let p = |i: usize| (mesh.points[2 * v[i] as usize], mesh.points[2 * v[i] as usize + 1]);
            let (a, b, c) = (p(0), p(1), p(2));
            0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs()
        })
        .sum()
}

fn labels(mesh: &femcore::mesh::SimplexMesh) -> (BTreeSet<i32>, BTreeSet<i32>) {
    (mesh.cell_labels.iter().copied().collect(), mesh.facet_labels.iter().copied().collect())
}

#[test]
fn the_simple_mesher_through_the_interface_is_the_old_path() {
    let m = core_in_clad(0.4);
    let a = mesh_model(&m.list, &MeshPlan::simple(0.15)).unwrap();
    let b = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.15)).unwrap();
    assert_eq!(a.mesh.cells, b.mesh.cells);
    assert_eq!(a.mesh.points, b.mesh.points);
    assert_eq!(a.mesh.facet_labels, b.mesh.facet_labels);
}

#[test]
fn adaptation_keeps_names_labels_and_geometry_vertices() {
    let m = core_in_clad(0.4);
    let simple = mesh_model(&m.list, &MeshPlan::simple(0.2)).unwrap();
    let adapted = mesh_model(&m.list, &MeshPlan::adaptive(0.2, 0.5)).unwrap();
    assert_eq!(labels(&simple.mesh), labels(&adapted.mesh), "the same domain and boundary numbers");
    assert_eq!(simple.domain("core"), adapted.domain("core"));
    assert!(adapted.mesh.n_cells() > 2 * simple.mesh.n_cells(), "{} vs {}", adapted.mesh.n_cells(), simple.mesh.n_cells());
    // the four corners of the cladding are still vertices
    for c in [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]] {
        assert!(
            (0..adapted.mesh.n_points()).any(|i| (adapted.mesh.points[2 * i] - c[0]).hypot(adapted.mesh.points[2 * i + 1] - c[1]) < 1e-12),
            "corner {c:?} lost"
        );
    }
    // the outside is still the square
    assert!((area(&adapted.mesh, None) - 4.0).abs() < 1e-12);
}

#[test]
fn adapted_interface_vertices_lie_on_the_circle_and_the_core_area_converges() {
    let r = 0.4;
    let m = core_in_clad(r);
    // a coarse chord polygon on purpose: max_dev h/4
    let mut plan = MeshPlan::adaptive(0.25, 0.4);
    plan.opts.max_dev = Some(0.25 / 4.0);
    let simple_plan = MeshPlan { mesher: Mesher::Simple, ..plan };
    let simple = mesh_model(&m.list, &simple_plan).unwrap();
    let adapted = mesh_model(&m.list, &plan).unwrap();
    let core = adapted.domain("core").unwrap();
    // every vertex between the two domains is on the circle
    let mut on_interface = vec![false; adapted.mesh.n_points()];
    let mut seen = std::collections::HashMap::new();
    for c in 0..adapted.mesh.n_cells() {
        let v = adapted.mesh.cell(c);
        for k in 0..3 {
            let e = (v[k].min(v[(k + 1) % 3]), v[k].max(v[(k + 1) % 3]));
            seen.entry(e).or_insert_with(Vec::new).push(adapted.mesh.cell_labels[c]);
        }
    }
    for ((a, b), ls) in &seen {
        if ls.len() == 2 && ls[0] != ls[1] {
            on_interface[*a as usize] = true;
            on_interface[*b as usize] = true;
        }
    }
    let mut n = 0;
    for (i, on) in on_interface.iter().enumerate() {
        if *on {
            n += 1;
            let rr = adapted.mesh.points[2 * i].hypot(adapted.mesh.points[2 * i + 1]);
            assert!((rr - r).abs() < 1e-9, "interface vertex at radius {rr}");
        }
    }
    assert!(n > 20, "only {n} interface vertices");
    // the core area is closer to the disk than the coarse chord polygon
    let exact = std::f64::consts::PI * r * r;
    let e_simple = (area(&simple.mesh, simple.domain("core")) - exact).abs();
    let e_adapted = (area(&adapted.mesh, Some(core)) - exact).abs();
    assert!(e_adapted < 0.5 * e_simple, "area error {e_adapted} vs {e_simple}");
    let q = mesh_quality_of(&adapted.mesh);
    assert!(q > 0.2, "worst triangle quality {q}");
}

fn mesh_quality_of(mesh: &femcore::mesh::SimplexMesh) -> f64 {
    let a = femcore::meshers::interop::simplex_to_adapt::<2>(mesh).unwrap();
    mesh_quality(&a, None).min
}

#[test]
fn poisson_on_the_adapted_disk_is_at_least_as_accurate_as_on_the_simple_one() {
    // -div(eps0 grad V) = eps0 on a disk of radius R, V = 0 on the rim: V max = R^2/4
    let r = 0.5;
    let m = disk(r);
    let run = |plan: &MeshPlan| {
        let g = mesh_model(&m.list, plan).unwrap();
        let mut s = ProblemSpec::new(Physics::Electrostatics);
        s.materials.default.insert("rho_v".into(), femcore::scalar::EPS0);
        s.bcs = vec![Bc { on: None, labels: g.select("outer").unwrap(), kind: BcKind::Ground }];
        let res = femcore::scalar::solve_electrostatics(&g.mesh, &s).unwrap();
        ((res.scalar("V max").unwrap() / (r * r / 4.0) - 1.0).abs(), g.mesh.n_cells())
    };
    let (e_simple, n_simple) = run(&MeshPlan::simple(0.15));
    let (e_adapt, n_adapt) = run(&MeshPlan::adaptive(0.15, 0.5));
    let (e_fine, n_fine) = run(&MeshPlan::simple(0.15 * 0.5));
    eprintln!("poisson disk: simple h=0.15 {e_simple:.2e} ({n_simple} cells); adaptive 0.5x {e_adapt:.2e} ({n_adapt}); simple h=0.075 {e_fine:.2e} ({n_fine})");
    assert!(n_adapt > n_simple);
    assert!(e_adapt < e_simple, "adapted {e_adapt} ({n_adapt} cells) vs simple {e_simple} ({n_simple} cells)");
}

#[test]
fn the_adaptive_mesher_is_deterministic() {
    let m = core_in_clad(0.4);
    let a = mesh_model(&m.list, &MeshPlan::adaptive(0.2, 0.5)).unwrap();
    let b = mesh_model(&m.list, &MeshPlan::adaptive(0.2, 0.5)).unwrap();
    assert_eq!(a.mesh.cells, b.mesh.cells);
    assert_eq!(a.mesh.points, b.mesh.points);
}

#[test]
fn a_cell_cap_is_reported_and_leaves_a_valid_mesh() {
    let m = disk(0.5);
    let mut plan = MeshPlan::adaptive(0.2, 0.05);
    plan.mesher = Mesher::Adaptive(Adaptive { scale: 0.05, max_cells: Some(3000), ..Adaptive::default() });
    let (g, rep) = femcore::meshing::mesh_evaluated_report(evaluate(&m.list).unwrap(), &plan).unwrap();
    let rep = rep.unwrap();
    assert!(rep.capped);
    assert!(g.mesh.n_cells() > 0);
}

/// The real use case: the quasi-TE index of a 500 x 300 nm SOI wire on the simple mesher, on its a priori graded adaptation, and against a
/// reference on a mesh with every size halved. Prints the table; asserts only what must hold (same labels, close to the reference).
#[test]
fn soi_wire_modes_on_the_simple_and_the_adaptive_mesh() {
    use femcore::usecases::soi::SoiWire;
    let wire = SoiWire::default();
    let te = |plan: &MeshPlan, w: &SoiWire| {
        let g = mesh_model(&w.model().list, plan).unwrap();
        let spec = w.spec(&g);
        let res = femcore::solve::solve(&g.mesh, &spec, &mut SolveState::default()).unwrap();
        let neff = res.solutions.iter().filter_map(|s| s.neff.map(|n| n.0)).fold(0.0, f64::max);
        (neff, g.mesh.n_cells(), labels(&g.mesh))
    };
    let (n_simple, c_simple, l_simple) = te(&MeshPlan::simple(wire.h_far), &wire);
    let (n_ad1, c_ad1, l_ad1) = te(&MeshPlan::adaptive(wire.h_far, 1.0), &wire);
    let (n_ad2, c_ad2, l_ad2) = te(&MeshPlan::adaptive(wire.h_far, 0.6), &wire);
    let fine = SoiWire { h_core: wire.h_core * 0.4, h_near: wire.h_near * 0.4, h_far: wire.h_far * 0.4, ..wire.clone() };
    let (n_ref, c_ref, _) = te(&MeshPlan::simple(fine.h_far), &fine);
    eprintln!("SOI quasi-TE n_eff (reference {n_ref:.6}, {c_ref} cells):");
    eprintln!("  simple            {n_simple:.6}  err {:.1e}  {c_simple} cells", (n_simple - n_ref).abs());
    eprintln!("  adaptive x1.0     {n_ad1:.6}  err {:.1e}  {c_ad1} cells", (n_ad1 - n_ref).abs());
    eprintln!("  adaptive x0.6     {n_ad2:.6}  err {:.1e}  {c_ad2} cells", (n_ad2 - n_ref).abs());
    assert_eq!(l_simple, l_ad1);
    assert_eq!(l_simple, l_ad2);
    assert!((n_ad2 - n_ref).abs() < 2e-3);
}
