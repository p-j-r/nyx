//! The unified meshing interface: the simple and the adaptive mesher make meshes of the same model with the same names and labels, and
//! the adaptive one puts its boundary and interface vertices on the exact curves of the geometry.
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::{mesh_model, Adaptive, MeshPlan, Mesher};
use femcore::mesh::quality::mesh_quality;
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
    let a = femcore::adapt::work::simplex_to_adapt::<2>(mesh).unwrap();
    mesh_quality(&a, None).min
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

// ------------------------------------------------------------------ estimator and loop

use femcore::adapt::estimate::recovered_hessians;

#[test]
fn the_recovered_hessian_of_a_quadratic_is_its_constant_hessian() {
    // u = x^2 + 3 z^2 + x z has the Hessian [[2, 1], [1, 6]] (positive definite, so its absolute value is itself)
    let m = Model::new_2d();
    let mut m = m;
    let r = m.rect("plate", (0.0, 0.0), (1.0, 1.0));
    m.overlay(&[("plate", r)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(0.08)).unwrap();
    let u: Vec<f64> = (0..g.mesh.n_points()).map(|i| {
        let (x, z) = (g.mesh.points[2 * i], g.mesh.points[2 * i + 1]);
        x * x + 3.0 * z * z + x * z
    }).collect();
    let h = recovered_hessians(&g.mesh, &u).unwrap();
    let mut worst: f64 = 0.0;
    let mut n = 0;
    for i in 0..g.mesh.n_points() {
        let (x, z) = (g.mesh.points[2 * i], g.mesh.points[2 * i + 1]);
        if x > 0.2 && x < 0.8 && z > 0.2 && z < 0.8 {
            n += 1;
            for (a, b, e) in [(0, 0, 2.0), (0, 1, 1.0), (1, 1, 6.0)] {
                worst = worst.max((h[i][a][b] - e).abs() / e);
            }
        }
    }
    assert!(n > 20);
    assert!(worst < 0.1, "worst relative error of the interior Hessian {worst}");
}

#[test]
fn recovery_does_not_smear_a_kink_over_an_interface() {
    // u = x on the left, 3x on the right of x = 0: continuous, each side linear. One-sided recovery sees no curvature anywhere.
    let mut m = Model::new_2d();
    let l = m.rect("left", (-1.0, 0.0), (0.0, 1.0));
    let r = m.rect("right", (0.0, 0.0), (1.0, 1.0));
    m.overlay(&[("left", l), ("right", r)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(0.1)).unwrap();
    let u: Vec<f64> = (0..g.mesh.n_points()).map(|i| {
        let x = g.mesh.points[2 * i];
        if x > 0.0 { 3.0 * x } else { x }
    }).collect();
    let h = recovered_hessians(&g.mesh, &u).unwrap();
    let on_interface: Vec<usize> = (0..g.mesh.n_points()).filter(|&i| g.mesh.points[2 * i].abs() < 1e-12).collect();
    assert!(on_interface.len() > 5);
    for i in on_interface {
        let n = h[i][0][0].abs() + h[i][1][1].abs() + h[i][0][1].abs();
        assert!(n < 1e-9, "Hessian {:?} at the interface vertex {i}", h[i]);
    }
}

// ------------------------------------------------------------------ caching

use femcore::meshing::{MeshCache, Reuse};

#[test]
fn the_cache_answers_repeats_and_equal_models_without_meshing() {
    let mut cache = MeshCache::new(4);
    let plan = MeshPlan::simple(0.15);
    let (a, r1) = cache.mesh(&core_in_clad(0.4).list, &plan).unwrap();
    assert!(matches!(r1, Reuse::Meshed { change: None }));
    // a separately built, equal model hits
    let (b, r2) = cache.mesh(&core_in_clad(0.4).list, &plan).unwrap();
    assert_eq!(r2, Reuse::Hit);
    assert_eq!(a.mesh.cells, b.mesh.cells);
    assert_eq!(a.mesh.points, b.mesh.points);
    // a changed radius, and a changed plan, mesh again; the change report says what kind of change it was
    let (_, r3) = cache.mesh(&core_in_clad(0.45).list, &plan).unwrap();
    assert!(matches!(r3, Reuse::Meshed { change: Some(_) }));
    let (_, r4) = cache.mesh(&core_in_clad(0.4).list, &MeshPlan::simple(0.2)).unwrap();
    assert!(matches!(r4, Reuse::Meshed { .. }));
    assert_eq!(cache.stats.requests, 4);
    assert_eq!(cache.stats.hits, 1);
}

#[test]
fn extruding_through_other_layers_reuses_the_section() {
    use femcore::meshing::mesher::Layer;
    let mut cache = MeshCache::default();
    let plan = MeshPlan::simple(0.2);
    let list = core_in_clad(0.4).list;
    let layers = |t: f64| vec![Layer { thickness: t, n: 3, domain_offset: 0 }];
    let (a, r1) = cache.mesh_extruded(&list, &plan, &layers(1.0)).unwrap();
    let (b, r2) = cache.mesh_extruded(&list, &plan, &layers(2.0)).unwrap();
    assert!(matches!(r1, Reuse::Meshed { .. }));
    assert_eq!(r2, Reuse::Hit);
    assert_eq!(a.mesh.n_cells(), b.mesh.n_cells());
    let top = |g: &femcore::meshing::geomdrive::MeshedGeom| (0..g.mesh.n_points()).map(|i| g.mesh.points[3 * i + 1]).fold(0.0, f64::max);
    assert!((top(&a) - 1.0).abs() < 1e-12 && (top(&b) - 2.0).abs() < 1e-12);
}

/// Study (run with `--ignored --nocapture`): a 20-point sweep of the SOI wire width, meshed with and without the cache, then the same sweep
/// again as a repeat.
#[test]
#[ignore]
fn cache_sweep_timing() {
    use femcore::usecases::soi::SoiWire;
    let widths: Vec<f64> = (0..20).map(|i| 400e-9 + 10e-9 * i as f64).collect();
    let lists: Vec<_> = widths.iter().map(|&w| SoiWire { width: w, ..SoiWire::default() }.model().list).collect();
    let plan = MeshPlan::simple(0.3e-6);
    let t = std::time::Instant::now();
    for l in &lists {
        mesh_model(l, &plan).unwrap();
    }
    let cold = t.elapsed().as_secs_f64();
    let mut cache = MeshCache::new(32);
    let t = std::time::Instant::now();
    for l in &lists {
        cache.mesh(l, &plan).unwrap();
    }
    let first = t.elapsed().as_secs_f64();
    let t = std::time::Instant::now();
    for l in &lists {
        cache.mesh(l, &plan).unwrap();
    }
    let again = t.elapsed().as_secs_f64();
    eprintln!("20 widths: no cache {cold:.3}s; cache first pass {first:.3}s (eval {:.3}s, mesh {:.3}s); repeat {again:.4}s", cache.stats.eval_seconds, cache.stats.mesh_seconds);
}

// ---------------------------------------------------------------- structured absorbing layers

/// A two-material box with an absorbing frame of thickness 0.5 and four layers, window x in [-1.5, 1.5], z in [-0.5, 0.8].
fn boxed_pml() -> Model {
    let mut m = Model::new_2d();
    let air = m.rect("air", (-2.0, -1.0), (2.0, 1.3));
    let sub = m.rect("sub", (-2.0, -1.0), (2.0, 0.0));
    let win = m.rect("window", (-1.5, -0.5), (1.5, 0.8));
    m.overlay(&[("air", air), ("sub", sub)]);
    m.pml("pml", win, &[], 0.5, 4.0, 2.0);
    m
}

fn label_areas(mesh: &femcore::mesh::SimplexMesh) -> std::collections::BTreeMap<i32, f64> {
    mesh.cell_labels.iter().copied().collect::<BTreeSet<_>>().into_iter().map(|l| (l, area(mesh, Some(l)))).collect()
}

fn structured(plan: MeshPlan) -> MeshPlan {
    MeshPlan { pml: femcore::meshing::PmlMode::Structured, ..plan }
}

#[test]
fn structured_layers_have_the_labels_areas_and_walls_of_the_inline_ones() {
    let m = boxed_pml();
    let inline = mesh_model(&m.list, &MeshPlan::simple(0.25)).unwrap();
    let st = mesh_model(&m.list, &structured(MeshPlan::simple(0.25))).unwrap();
    // the same names mean the same numbers, and each piece fills the same area
    assert_eq!(st.labels, inline.labels);
    let (ai, as_) = (label_areas(&inline.mesh), label_areas(&st.mesh));
    assert_eq!(ai.keys().collect::<Vec<_>>(), as_.keys().collect::<Vec<_>>());
    for (l, a) in &ai {
        assert!((a - as_[l]).abs() < 1e-9 * 16.0, "label {l}: inline {a} structured {}", as_[l]);
    }
    assert!(!st.domains("air*").is_empty() && st.pml_materials(3.0).len() == inline.pml_materials(3.0).len());
    // every outer edge has exactly one facet, interior edges none, and the facets carry the labels of the inline mesh's walls
    let fc = |mesh: &femcore::mesh::SimplexMesh| {
        let mut count = std::collections::HashMap::new();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            for k in 0..3 {
                let (a, b) = (v[k].min(v[(k + 1) % 3]), v[k].max(v[(k + 1) % 3]));
                *count.entry((a, b)).or_insert(0usize) += 1;
            }
        }
        count
    };
    let count = fc(&st.mesh);
    let mut seen = BTreeSet::new();
    for f in 0..st.mesh.n_facets() {
        let v = st.mesh.facet(f);
        let key = (v[0].min(v[1]), v[0].max(v[1]));
        assert_eq!(count[&key], 1, "a facet on an interior edge");
        assert!(seen.insert(key), "a facet twice");
    }
    assert_eq!(seen.len(), count.values().filter(|&&n| n == 1).count(), "an outer edge without a facet");
    assert_eq!(labels(&st.mesh).1, labels(&inline.mesh).1, "the same boundary numbers");
    // four layers across the thickness on the right: the x coordinates beyond the window are 4 distinct values
    let mut xs: Vec<f64> = (0..st.mesh.n_points()).map(|i| st.mesh.points[2 * i]).filter(|&x| x > 1.5 + 1e-9).collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    assert_eq!(xs.len(), 4, "{xs:?}");
    assert!((xs[3] - 2.0).abs() < 1e-12);
    // no cell is inverted
    let q = mesh_quality_of(&st.mesh);
    assert!(q > 0.0, "worst cell quality {q}");
}

#[test]
fn structured_layers_follow_an_adaptive_window() {
    let m = boxed_pml();
    let plan = structured(MeshPlan::adaptive(0.25, 0.5));
    let st = mesh_model(&m.list, &plan).unwrap();
    let inline = mesh_model(&m.list, &MeshPlan::adaptive(0.25, 0.5)).unwrap();
    assert_eq!(labels(&st.mesh).0, labels(&inline.mesh).0);
    let (ai, as_) = (label_areas(&inline.mesh), label_areas(&st.mesh));
    for (l, a) in &ai {
        assert!((a - as_[l]).abs() < 1e-8, "label {l}: {a} vs {}", as_[l]);
    }
    // the physical window kept its refinement; the layers are one quad column per boundary edge, so there are far fewer layer cells
    assert!(st.mesh.n_cells() <= inline.mesh.n_cells(), "{} vs {}", st.mesh.n_cells(), inline.mesh.n_cells());
}

#[test]
fn a_model_without_layers_is_left_alone_and_a_round_window_is_refused() {
    let m = core_in_clad(0.4);
    let plain = mesh_model(&m.list, &MeshPlan::simple(0.2)).unwrap();
    let st = mesh_model(&m.list, &structured(MeshPlan::simple(0.2))).unwrap();
    assert_eq!(plain.mesh.cells, st.mesh.cells);
    // a window with round corners is not a box: the structured builder says so instead of making a bad mesh
    let mut r = Model::new_2d();
    let air = r.rect("air", (-2.0, -2.0), (2.0, 2.0));
    let win = r.circle("window", (0.0, 0.0), 1.2);
    r.overlay(&[("air", air)]);
    r.pml("pml", win, &[], 0.5, 4.0, 2.0);
    let e = mesh_model(&r.list, &structured(MeshPlan::simple(0.25))).err();
    assert!(e.is_some(), "a round window must be refused");
}

/// 3D by sweeping the 2D section: the adapted (and the structured-layer) section extruded along y fills the right volume with positive tets.
#[test]
fn an_adapted_section_extruded_to_3d_has_the_volume_and_the_labels() {
    use femcore::meshing::mesher::Layer;
    let layers = vec![Layer { thickness: 0.6, n: 3, domain_offset: 0 }];
    let vol = |m: &femcore::mesh::SimplexMesh| -> (f64, f64) {
        let (mut total, mut min) = (0.0, f64::INFINITY);
        for c in 0..m.n_cells() {
            let v = m.cell(c);
            let p = |i: usize| -> [f64; 3] { std::array::from_fn(|k| m.points[3 * v[i] as usize + k]) };
            let (a, b, c2, d) = (p(0), p(1), p(2), p(3));
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let w = [c2[0] - a[0], c2[1] - a[1], c2[2] - a[2]];
            let x = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
            let det = u[0] * (w[1] * x[2] - w[2] * x[1]) - u[1] * (w[0] * x[2] - w[2] * x[0]) + u[2] * (w[0] * x[1] - w[1] * x[0]);
            total += det.abs() / 6.0;
            min = f64::min(min, det.abs() / 6.0);
        }
        (total, min)
    };
    for (name, list, plan) in [
        ("adaptive disk", core_in_clad(0.4).list, MeshPlan::adaptive(0.2, 0.5)),
        ("structured layers", boxed_pml().list, structured(MeshPlan::adaptive(0.25, 0.5))),
    ] {
        let mut cache = MeshCache::default();
        let (g2, _) = cache.mesh(&list, &plan).unwrap();
        let (g3, _) = cache.mesh_extruded(&list, &plan, &layers).unwrap();
        let (v, min) = vol(&g3.mesh);
        let expect = area(&g2.mesh, None) * 0.6;
        assert!((v - expect).abs() < 1e-9 * expect.max(1.0), "{name}: volume {v} vs {expect}");
        assert!(min > 0.0, "{name}: smallest tet {min}");
        assert_eq!(g3.mesh.dim, 3);
        assert_eq!(g3.mesh.n_cells(), 3 * 3 * g2.mesh.n_cells(), "{name}: three tets per prism, three layers");
        eprintln!("{name}: {} triangles -> {} tets", g2.mesh.n_cells(), g3.mesh.n_cells());
    }
}

#[test]
fn native_structured_layers_keep_their_quadrilaterals_and_split_to_the_simplex_mesh() {
    use femcore::meshing::pml_geom::{structured_pml, structured_pml_native, Native};
    let m = boxed_pml();
    let start = mesh_model(&m.list, &MeshPlan::simple(0.25)).unwrap();
    let (split, report) = structured_pml(start.clone()).unwrap();
    let Native::Cells(cells, rep2) = structured_pml_native(start).unwrap() else { panic!("the model has a Pml") };
    assert_eq!(report, rep2);
    use femcore::mesh::cells::CellType;
    let quads: Vec<usize> = (0..cells.n_cells()).filter(|&c| cells.types[c] == CellType::Quad).collect();
    let tris: Vec<usize> = (0..cells.n_cells()).filter(|&c| cells.types[c] == CellType::Tri).collect();
    assert_eq!(quads.len(), report.layer_quads, "the layers are quadrilaterals");
    assert_eq!(tris.len(), report.physical_cells, "the window stays triangles");
    // the cut mesh is the native one with each quadrilateral in two
    assert_eq!(split.mesh.n_cells(), tris.len() + 2 * quads.len());
    // the layer quadrilaterals have positive area (native order is counter-clockwise: 0, 1, 2, 3) and the labels geom gave the pieces
    for &q in &quads {
        let v = cells.cell(q);
        let (a, b, c, d) = (cells.points[v[0] as usize], cells.points[v[1] as usize], cells.points[v[2] as usize], cells.points[v[3] as usize]);
        let area = 0.5 * ((a[0] * b[1] - b[0] * a[1]) + (b[0] * c[1] - c[0] * b[1]) + (c[0] * d[1] - d[0] * c[1]) + (d[0] * a[1] - a[0] * d[1]));
        assert!(area > 0.0, "quadrilateral {q} has area {area} (counter-clockwise)");
        assert!(split.mesh.cell_labels.contains(&cells.label[q]));
    }
}
