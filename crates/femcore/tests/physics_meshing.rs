//! Meshes sized by the physics of the study (`physics_size`, `meshing::mesh_for_problem`): the wavelength per domain sets the element
//! size, and the mode index the solver then finds agrees with the closed form of a PEC guide.
use femcore::ems::{solve_emw, C0};
use femcore::geomdrive::evaluate;
use femcore::meshing::{mesh_for_problem, Adaptive, MeshPlan, Mesher};
use femcore::physics_size::PhysicsSizing;
use femcore::problem::*;
use geom::model::{Model, RegionKind};

const WR90: (f64, f64) = (22.86e-3, 10.16e-3);

fn guide(eps: f64, freq: f64, order: usize, n_target: f64) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = order;
    s.materials.default.insert("eps_r".into(), eps);
    s.study = Study::Modes { freq, n_target, evanescent: false, n_modes: 3 };
    s
}

fn te10_neff(eps: f64, freq: f64) -> f64 {
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    (eps - (std::f64::consts::PI / WR90.0 / k0).powi(2)).sqrt()
}

fn plate_model() -> Model {
    let mut m = Model::new_2d();
    let r = m.rect("plate", (0.0, 0.0), (WR90.0, WR90.1));
    m.overlay(&[("plate", r)]);
    m
}

/// Error of the TE10 index on the mesh the physics asks for, and the number of cells.
fn neff_error(order: usize, ppw: Option<f64>) -> (f64, usize) {
    let (eps, freq) = (2.1, 12e9);
    let exact = te10_neff(eps, freq);
    let spec = guide(eps, freq, order, exact * 1.001);
    // a start size far coarser than the wave needs: the physics must refine it
    let plan = MeshPlan { mesher: Mesher::Adaptive(Adaptive { max_cells: Some(200_000), ..Adaptive::default() }), ..MeshPlan::simple(WR90.1 / 1.5) };
    let sz = PhysicsSizing { points_per_wavelength: ppw, ..PhysicsSizing::default() };
    let (g, _) = mesh_for_problem(&plate_model().list, &plan, &spec, &sz).unwrap();
    let r = solve_emw(&g.mesh, &spec).unwrap();
    let best = r.solutions.iter().map(|s| (s.neff.unwrap().0 - exact).abs()).fold(f64::MAX, f64::min);
    (best, g.mesh.n_cells())
}

#[test]
fn default_points_per_wavelength_resolve_a_guide_mode_at_both_orders() {
    for order in [1, 2] {
        let (err, cells) = neff_error(order, None);
        eprintln!("order {order}: {cells} cells, TE10 neff error {err:.2e}");
        assert!(err < 2e-3, "order {order}: error {err} on {cells} cells");
    }
}

#[test]
fn fewer_points_per_wavelength_give_a_larger_error() {
    let (coarse, nc) = neff_error(1, Some(4.0));
    let (fine, nf) = neff_error(1, Some(16.0));
    eprintln!("order 1: 4 ppw {nc} cells err {coarse:.2e}; 16 ppw {nf} cells err {fine:.2e}");
    assert!(nf > 4 * nc / 2 && fine < coarse / 4.0, "coarse {coarse} ({nc} cells), fine {fine} ({nf} cells)");
}

#[test]
fn a_higher_index_domain_gets_a_finer_mesh() {
    // a core of index 3 in a cladding of index 1: the element size in the core is a third of the cladding's
    let mut m = Model::new_2d();
    let clad = m.rect("clad", (-1e-3, -1e-3), (1e-3, 1e-3));
    let core = m.rect("core", (-0.3e-3, -0.3e-3), (0.3e-3, 0.3e-3));
    m.overlay_full(&[("clad", clad, RegionKind::Material, None), ("core", core, RegionKind::Material, None)], &[]);
    let ev = evaluate(&m.list).unwrap();
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.order = 1;
    spec.study = Study::Modes { freq: 100e9, n_target: 2.0, evanescent: false, n_modes: 2 };
    let g0 = femcore::geomdrive::mesh_geom_2d(ev, &femcore::geomdrive::GeomMeshOpts::new(0.5e-3)).unwrap();
    let (cl, co) = (g0.domain("clad").unwrap(), g0.domain("core").unwrap());
    spec.materials.set(cl, "eps_r", 1.0);
    spec.materials.set(co, "eps_r", 9.0);
    let plan = MeshPlan { mesher: Mesher::Adaptive(Adaptive::default()), ..MeshPlan::simple(0.5e-3) };
    let (g, _) = mesh_for_problem(&m.list, &plan, &spec, &PhysicsSizing::default()).unwrap();
    // the cladding is measured away from the core: the size grows from the core's with the grading slope, so near it the cladding is finer
    let mean_edge = |label: i32| {
        let mesh = &g.mesh;
        let (mut sum, mut n) = (0.0, 0);
        for c in (0..mesh.n_cells()).filter(|&c| mesh.cell_labels[c] == label) {
            let v = mesh.cell(c);
            let cen = v.iter().map(|&i| mesh.point(i as usize)[0].abs().max(mesh.point(i as usize)[1].abs())).sum::<f64>() / 3.0;
            if label == cl && cen < 0.8e-3 {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (mesh.point(v[k] as usize), mesh.point(v[(k + 1) % 3] as usize));
                sum += (a[0] - b[0]).hypot(a[1] - b[1]);
                n += 1;
            }
        }
        sum / n as f64
    };
    let ratio = mean_edge(cl) / mean_edge(co);
    eprintln!("mean edge: clad {:.3e}, core {:.3e}, ratio {ratio:.2}", mean_edge(cl), mean_edge(co));
    assert!((2.4..3.8).contains(&ratio), "edge-length ratio {ratio} (expected about 3, the index ratio)");
}

// ---- vector-field estimator -------------------------------------------------------------------------------------------------------

use femcore::estimate::{estimator_metrics, Estimator};
use femcore::geomdrive::{mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::adapt_to_metrics;

/// A dielectric core (eps 9) in a PEC-walled cladding (eps 1), 1 mm square.
fn loaded_box() -> Model {
    let mut m = Model::new_2d();
    let clad = m.rect("clad", (-0.5e-3, -0.5e-3), (0.5e-3, 0.5e-3));
    let core = m.rect("core", (-0.15e-3, -0.15e-3), (0.15e-3, 0.15e-3));
    m.overlay_full(&[("clad", clad, RegionKind::Material, None), ("core", core, RegionKind::Material, None)], &[]);
    m
}

fn loaded_spec(g: &MeshedGeom) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.set(g.domain("clad").unwrap(), "eps_r", 1.0);
    s.materials.set(g.domain("core").unwrap(), "eps_r", 9.0);
    s.study = Study::Modes { freq: 150e9, n_target: 2.5, evanescent: false, n_modes: 1 };
    s
}

fn first_neff(g: &MeshedGeom) -> (SolveResult, f64) {
    let r = solve_emw(&g.mesh, &loaded_spec(g)).unwrap();
    let n = r.solutions[0].neff.unwrap().0;
    (r, n)
}

#[test]
fn e_field_components_as_one_vector_need_fewer_cells_than_intersected_components() {
    let m = loaded_box();
    let opts = GeomMeshOpts::new(0.12e-3);
    let start = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
    let (res, n0) = first_neff(&start);
    let reference = {
        let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.015e-3)).unwrap();
        first_neff(&g).1
    };
    let target = 2.0 * start.mesh.n_cells() as f64;
    let comps = ["Ex", "Ey", "Ez"].map(String::from).to_vec();
    let ests = [
        ("components, intersected", Estimator { fields: comps.clone(), combine: false, ..Estimator::default() }),
        ("components, one vector", Estimator::e_field()),
        ("|E| (default field)", Estimator::default()),
    ];
    let mut cells = vec![];
    for (name, est) in ests {
        let metrics = estimator_metrics(&start, &res, &est, &opts, target).unwrap();
        let (g, _) = adapt_to_metrics(mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap(), &metrics, &Adaptive::default()).unwrap();
        let n = first_neff(&g).1;
        eprintln!("{name:26} {:6} cells  neff err {:.2e}   (start mesh {} cells err {:.2e})", g.mesh.n_cells(), (n - reference).abs(), start.mesh.n_cells(), (n0 - reference).abs());
        cells.push((g.mesh.n_cells(), (n - reference).abs()));
    }
    assert!(cells[1].0 < cells[0].0, "one vector metric should not ask for more cells than three intersected ones: {cells:?}");
    assert!(cells[1].1 < (n0 - reference).abs() || (n0 - reference).abs() < 1e-6, "adapting to the E components must not worsen the index: {cells:?}");
}

// ---- 3D: a PEC cavity meshed from its physics ----------------------------------------------------------------------------------------

mod common;
use femcore::meshers::adapt::ops3::AdaptOptions3;
use femcore::physics_size::adapt_tets;

fn cavity_modes(a: f64, b: f64, d: f64, fmax: f64) -> Vec<f64> {
    let mut out = vec![];
    for m in 0..5 {
        for n in 0..5 {
            for p in 0..5 {
                let nonzero = (m > 0) as usize + (n > 0) as usize + (p > 0) as usize;
                if nonzero < 2 {
                    continue;
                }
                let f = 0.5 * C0 * ((m as f64 / a).powi(2) + (n as f64 / b).powi(2) + (p as f64 / d).powi(2)).sqrt();
                if f <= fmax {
                    out.push(f);
                    if nonzero == 3 {
                        out.push(f);
                    }
                }
            }
        }
    }
    out.sort_by(|x, y| x.partial_cmp(y).unwrap());
    out
}

#[test]
fn a_cavity_meshed_from_its_physics_reproduces_the_first_modes() {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, d, 20e9);
    for order in [1usize, 2] {
        let start = common::block(a, b, d, 2, 2, 2);
        let mut spec = ProblemSpec::new(Physics::Emw);
        spec.order = order;
        spec.materials.default.insert("eps_r".into(), 1.0);
        spec.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.9 * exact[0] };
        let (mesh, report) = adapt_tets(&start, &spec, &PhysicsSizing::default(), a, 0.4, &AdaptOptions3::default()).unwrap();
        let r = femcore::solve::solve(&mesh, &spec, &mut SolveState::default()).unwrap();
        let got: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
        let worst = got.iter().zip(&exact).map(|(g, e)| (g / e - 1.0).abs()).fold(0.0, f64::max);
        eprintln!("order {order}: {} -> {} tets, worst of 4 modes {worst:.2e}, adaptation iterations {}", start.n_cells(), mesh.n_cells(), report.iters.len());
        assert!(worst < 5e-3, "order {order}: modes {got:?} vs {exact:?}");
        assert!(got[0] > 0.95 * exact[0], "a spurious mode below the first: {got:?}");
    }
}

#[test]
fn a_real_stretch_layer_on_an_adapted_cavity_equals_the_longer_cavity() {
    use femcore::meshers::adapt::pml::PmlSpec;
    use femcore::pml3::{add_pml_3d, apply_materials};
    // PML with the real stretch 2 over a layer of d/2 on the +z end equals a cavity of length d + 2 * d / 2 = 2d. The physical block is
    // first adapted to the physics; its box faces, edges and corners stay exact, so the layers can be swept from it
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, 2.0 * d, 20e9);
    let start = common::block(a, b, d, 2, 2, 2);
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.order = 2;
    spec.materials.default.insert("eps_r".into(), 1.0);
    spec.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.9 * exact[0] };
    let (phys, _) = adapt_tets(&start, &spec, &PhysicsSizing::default(), a, 0.4, &AdaptOptions3::default()).unwrap();
    let mut pspec = PmlSpec::<3>::uniform(0.0, 0);
    pspec.thickness[2][1] = 0.5 * d;
    pspec.layers[2][1] = 3;
    pspec.outer_label = Some(4);
    let pb = add_pml_3d(&phys, &pspec, 0.0, 0.0, 2.0).unwrap();
    assert!(pb.cells.ty("prism").is_some() && pb.split.prisms > 0, "the face layer is prisms");
    apply_materials(&pb, &mut spec.materials);
    // PML cells of the face region carry a z stretch only
    assert!(pb.pml.values().all(|p| p.z != femcore::problem::PmlDir::Off && p.x == femcore::problem::PmlDir::Off));
    let r = femcore::solve::solve(&pb.mesh, &spec, &mut SolveState::default()).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
    let worst = got.iter().zip(&exact).map(|(g, e)| (g / e - 1.0).abs()).fold(0.0, f64::max);
    eprintln!("adapted block {} tets + layers = {} tets ({} prisms cut), worst of 4 modes {worst:.2e}", phys.n_cells(), pb.mesh.n_cells(), pb.split.prisms);
    assert!(worst < 5e-3, "{got:?} vs {exact:?}");
}

#[test]
fn layers_on_all_six_sides_with_a_unit_stretch_equal_the_bigger_cavity() {
    use femcore::meshers::adapt::pml::PmlSpec;
    use femcore::pml3::{add_pml_3d, apply_materials};
    // re = 1 and no loss: the layers are plain vacuum. Faces are prisms, edges and corners hexahedra; all 26 layer domains are in play
    let (a, b, d, t) = (0.04, 0.03, 0.02, 0.01);
    let exact = cavity_modes(a + 2.0 * t, b + 2.0 * t, d + 2.0 * t, 20e9);
    let phys = common::block(a, b, d, 4, 3, 2);
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.order = 2;
    spec.materials.default.insert("eps_r".into(), 1.0);
    spec.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.9 * exact[0] };
    let mut pspec = PmlSpec::<3>::uniform(t, 2);
    pspec.outer_label = Some(1);
    let pb = add_pml_3d(&phys, &pspec, 0.0, 0.0, 1.0).unwrap();
    assert_eq!(pb.pml.len(), 26);
    assert!(pb.split.hexes > 0 && pb.split.prisms > 0);
    apply_materials(&pb, &mut spec.materials);
    let r = femcore::solve::solve(&pb.mesh, &spec, &mut SolveState::default()).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
    let worst = got.iter().zip(&exact).map(|(g, e)| (g / e - 1.0).abs()).fold(0.0, f64::max);
    eprintln!("{} tets ({} hexes, {} prisms cut): worst of 4 modes {worst:.2e}", pb.mesh.n_cells(), pb.split.hexes, pb.split.prisms);
    assert!(worst < 6e-3, "{got:?} vs {exact:?}");
}

/// Lowest three complex frequencies (Hz, as (Re f, Im f)) of the guide a x b closed by a PEC wall at z = 0 and by PML layers beyond
/// z = d, on a block with `n` cells along x and `layers` PML layers (the other directions scale with it).
fn pml_end_modes(n: usize, layers: usize, strength: f64, t: f64) -> Vec<(f64, f64)> {
    use femcore::meshers::adapt::pml::PmlSpec;
    use femcore::pml3::{add_pml_3d, apply_materials};
    let (a, b, d) = (0.04, 0.03, 0.02);
    let phys = common::block(a, b, d, n, n * 3 / 4, n * 3 / 4);
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.order = 2;
    spec.materials.default.insert("eps_r".into(), 1.0);
    let f101 = 0.5 * C0 * ((1.0 / a).powi(2) + (1.0 / (2.0 * d)).powi(2)).sqrt();
    spec.study = Study::Eigenfrequency { n_modes: 3, f_search: f101 };
    let mut pspec = PmlSpec::<3>::uniform(0.0, 0);
    pspec.thickness[2][1] = t;
    pspec.layers[2][1] = layers;
    pspec.outer_label = Some(4);
    let pb = add_pml_3d(&phys, &pspec, strength, 2.0, 1.0).unwrap();
    apply_materials(&pb, &mut spec.materials);
    let r = femcore::solve::solve(&pb.mesh, &spec, &mut SolveState::default()).unwrap();
    (1..=r.solutions.len())
        .map(|j| {
            let f = r.solutions[j - 1].value.unwrap();
            match r.scalar(&format!("Q{j}")) {
                Some(q) => (f, f / (2.0 * q)),
                None => (f, 0.0),
            }
        })
        .collect()
}

#[test]
#[ignore = "study: mesh convergence of the damped modes of a PML-terminated guide; run with --ignored --nocapture"]
fn pml_end_converges_with_the_mesh() {
    // the same continuous problem (PML thickness 2 cm, strength 3) on three meshes: the complex frequencies must converge
    let levels = [(4, 3), (6, 4), (8, 6), (12, 9)];
    let all: Vec<Vec<(f64, f64)>> = levels.iter().map(|&(n, l)| pml_end_modes(n, l, 3.0, 0.02)).collect();
    let fin = all.last().unwrap();
    for (lv, m) in levels.iter().zip(&all) {
        let err: Vec<String> = m.iter().zip(fin).map(|(x, r)| format!("{:.1e}", ((x.0 - r.0).powi(2) + (x.1 - r.1).powi(2)).sqrt() / r.0)).collect();
        eprintln!("n {} layers {}: {:?} error vs finest {err:?}", lv.0, lv.1, m.iter().map(|x| format!("{:.5e}{:+.3e}i", x.0, x.1)).collect::<Vec<_>>());
    }
}
