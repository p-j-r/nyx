// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Wavelength sizing, adaptive meshing and structured absorbing layers for 3D `geom` models (ADR 50, 51, 52): the three features that the
//! 2D mesher has and ADR 40 left out of the lattice mesher.
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, MeshedGeom};
use femcore::meshing::geomdrive3::{mesh_geom_3d_with, Sizes3};
use femcore::meshing::physics_size::PhysicsSizing;
use femcore::meshing::pml_geom3::{structured_pml_3d_native, Native3};
use femcore::meshing::{mesh_for_sizing, mesh_model, MeshPlan, PmlMode, Sizing};
use geom::model::Model;
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;

fn pt(m: &SimplexMesh, v: u32) -> [f64; 3] {
    [m.points[3 * v as usize], m.points[3 * v as usize + 1], m.points[3 * v as usize + 2]]
}

fn det6(p: [[f64; 3]; 4]) -> f64 {
    let s = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let (u, v, w) = (s(p[1], p[0]), s(p[2], p[0]), s(p[3], p[0]));
    u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])
}

/// Volume of the cells of every domain label; every cell must be positively oriented.
fn volumes(m: &SimplexMesh) -> BTreeMap<i32, f64> {
    let mut out: BTreeMap<i32, f64> = BTreeMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let d = det6([pt(m, v[0]), pt(m, v[1]), pt(m, v[2]), pt(m, v[3])]) / 6.0;
        assert!(d > 0.0, "cell {c} has volume {d}");
        *out.entry(m.cell_labels[c]).or_default() += d;
    }
    out
}

/// Longest edge of every cell of domain `label`.
fn longest_edges(m: &SimplexMesh, label: i32) -> Vec<f64> {
    let mut out = vec![];
    for c in 0..m.n_cells() {
        if m.cell_labels[c] != label {
            continue;
        }
        let v = m.cell(c);
        let mut l: f64 = 0.0;
        for i in 0..4 {
            for j in i + 1..4 {
                let (a, b) = (pt(m, v[i]), pt(m, v[j]));
                l = l.max(((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt());
            }
        }
        out.push(l);
    }
    out
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// Every face of every cell is shared by at most two cells; the faces used once are exactly the boundary facets that touch the void, and
/// the faces between cells of different domains are exactly the interface facets (ADR 40).
fn assert_conforming(m: &SimplexMesh) {
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
    let mut expected: Vec<[u32; 3]> = faces.iter().filter(|(_, v)| v.len() == 1 || m.cell_labels[v[0]] != m.cell_labels[v[1]]).map(|(k, _)| *k).collect();
    expected.sort_unstable();
    let mut got: Vec<[u32; 3]> = (0..m.n_facets()).map(|f| key(m.facets[3 * f], m.facets[3 * f + 1], m.facets[3 * f + 2])).collect();
    got.sort_unstable();
    assert_eq!(got, expected, "the boundary facets are not the faces of the void and of the interfaces");
    assert!(m.cell_labels.iter().all(|&l| l >= 1) && m.facet_labels.iter().all(|&l| l >= 1), "labels are 1-based");
}

fn ball_in_box(r: f64) -> Model {
    let mut m = Model::new_3d();
    let air = m.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let ball = m.sphere("ball", (2.0, 2.0, 2.0), r);
    m.overlay(&[("air", air), ("ball", ball)]);
    m
}

fn waves(order: u8, pairs: &[(&str, f64)]) -> Sizing {
    Sizing { order, headroom: 1.0, wavelengths: pairs.iter().map(|(n, l)| (n.to_string(), *l)).collect() }
}

// ---------------------------------------------------------------- wavelength sizing (ADR 50)

#[test]
fn wavelength_sizing_refines_only_the_regions_with_short_waves() {
    let m = ball_in_box(0.9);
    let plan = MeshPlan::simple(0.5);
    let plain = mesh_model(&m.list, &plan).unwrap();
    // order 2 takes 6 points per wavelength: the ball (wavelength 0.6) wants 0.1, the air (2.4) wants 0.4
    let (g, report) = mesh_for_sizing(&m.list, &plan, &waves(2, &[("ball", 0.6), ("air", 2.4)]), &PhysicsSizing::default()).unwrap();
    assert!(report.is_none());
    assert_conforming(&g.mesh);
    let (air, ball) = (g.domain("air").unwrap(), g.domain("ball").unwrap());
    let (le_ball, le_air) = (longest_edges(&g.mesh, ball), longest_edges(&g.mesh, air));
    // Rivara bisection keeps the longest edge within the ratio 1.3 of the wanted size, the stuffed lattice within a factor 2
    assert!(le_ball.iter().cloned().fold(0.0, f64::max) < 0.1 * 1.3 * 1.8, "ball edge {}", le_ball.iter().cloned().fold(0.0, f64::max));
    assert!(mean(&le_ball) < 0.13, "mean ball edge {}", mean(&le_ball));
    assert!(mean(&le_air) > 0.2, "the air stays coarse: {}", mean(&le_air));
    let count = |g: &MeshedGeom, l: i32| g.mesh.cell_labels.iter().filter(|&&x| x == l).count();
    assert!(count(&g, ball) > 20 * count(&plain, ball), "{} vs {}", count(&g, ball), count(&plain, ball));
    // the geometry is as exact as before
    let v = volumes(&g.mesh);
    assert!(rel(v[&air] + v[&ball], 64.0) < 1e-9);
    assert!(rel(v[&ball], 4.0 / 3.0 * PI * 0.9f64.powi(3)) < 0.01, "{v:?}");
}

#[test]
fn wavelength_sizing_is_a_ceiling_never_coarser_than_the_plan_and_names_must_exist() {
    let m = ball_in_box(0.9);
    let plan = MeshPlan::simple(0.5);
    let plain = mesh_model(&m.list, &plan).unwrap();
    // long waves everywhere: the plan's own size stays
    let (g, _) = mesh_for_sizing(&m.list, &plan, &waves(2, &[("ball", 60.0), ("air", 60.0)]), &PhysicsSizing::default()).unwrap();
    assert_eq!(g.mesh.n_cells(), plain.mesh.n_cells());
    // a region the model does not have is an error, not silence
    let e = mesh_for_sizing(&m.list, &plan, &waves(2, &[("core", 1.0)]), &PhysicsSizing::default()).err().unwrap();
    assert!(e.contains("core"), "{e}");
    // short waves everywhere: the lattice itself is made fine (not coarse and refined afterwards), so the cells are uniform in size
    let (g, _) = mesh_for_sizing(&m.list, &plan, &waves(2, &[("ball", 1.2), ("air", 1.2)]), &PhysicsSizing::default()).unwrap();
    let all: Vec<f64> = [1, 2].iter().flat_map(|&l| longest_edges(&g.mesh, l)).collect();
    assert!(all.iter().cloned().fold(0.0, f64::max) < 0.2 * 1.8, "{}", all.iter().cloned().fold(0.0, f64::max));
}

#[test]
fn absorbing_pieces_are_relaxed_by_pml_relax() {
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let win = m.box_("window", (-1.0, -1.0, -1.0), (1.0, 1.0, 1.0));
    m.overlay(&[("air", air)]);
    m.pml("pml", win, &[], 1.0, 2.0, 2.0);
    let plan = MeshPlan { pml: PmlMode::Inline, ..MeshPlan::simple(0.8) };
    let s = waves(1, &[("air", 2.4), ("air.pml_xmax", 2.4)]);
    let strict = PhysicsSizing { pml_relax: 1.0, ..PhysicsSizing::default() };
    let relaxed = PhysicsSizing { pml_relax: 2.0, ..PhysicsSizing::default() };
    let (a, _) = mesh_for_sizing(&m.list, &plan, &s, &strict).unwrap();
    let (b, _) = mesh_for_sizing(&m.list, &plan, &s, &relaxed).unwrap();
    let piece = a.domain("air.pml_xmax").unwrap();
    // the piece is one layer of 0.5 across (thickness / layers): that is its size whatever the waves ask, both ways
    assert!(longest_edges(&a.mesh, piece).iter().cloned().fold(0.0, f64::max) < 0.5 * 1.8);
    assert!(b.mesh.n_cells() <= a.mesh.n_cells());
}

// ---------------------------------------------------------------- adaptive meshing (ADR 51)

#[test]
fn the_adaptive_mesher_refines_a_sphere_and_lands_every_new_vertex_on_it() {
    let r = 0.9;
    let m = ball_in_box(r);
    let simple = mesh_model(&m.list, &MeshPlan::simple(0.6)).unwrap();
    let (g, report) = mesh_for_sizing_free(&m, MeshPlan::adaptive(0.6, 0.5));
    let report = report.expect("an adaptation report");
    assert!(report.iters.last().is_some());
    assert_conforming(&g.mesh);
    assert!(g.mesh.n_cells() > 3 * simple.mesh.n_cells(), "{} vs {}", g.mesh.n_cells(), simple.mesh.n_cells());
    let (air, ball) = (g.domain("air").unwrap(), g.domain("ball").unwrap());
    let iface = g.boundary("air|ball");
    assert_eq!(iface.len(), 1);
    // domains and boundaries have the numbers of the model
    assert_eq!((air, ball), (1, 2));
    let v = volumes(&g.mesh);
    assert!(rel(v[&air] + v[&ball], 64.0) < 1e-9, "{v:?}");
    // the vertices of the sphere are on the sphere, not on the chord polyhedron of the start mesh: this is what the oracle is for
    let mut worst: f64 = 0.0;
    let mut n = 0;
    for f in 0..g.mesh.n_facets() {
        if g.mesh.facet_labels[f] != iface[0] {
            continue;
        }
        n += 1;
        for k in 0..3 {
            let p = pt(&g.mesh, g.mesh.facets[3 * f + k]);
            worst = worst.max((((p[0] - 2.0).powi(2) + (p[1] - 2.0).powi(2) + (p[2] - 2.0).powi(2)).sqrt() - r).abs());
        }
    }
    assert!(n > 0 && worst < 1e-6, "{n} facets, a vertex is {worst} off the sphere");
    // the volume of the ball converges: finer than the start mesh
    let v0 = volumes(&simple.mesh)[&ball];
    let exact = 4.0 / 3.0 * PI * r.powi(3);
    assert!(rel(v[&ball], exact) < rel(v0, exact), "adapted {} start {}", rel(v[&ball], exact), rel(v0, exact));
    // the box faces are exact
    let a = |name: &str| g.boundary(name)[0];
    let area = |l: i32| {
        (0..g.mesh.n_facets())
            .filter(|&f| g.mesh.facet_labels[f] == l)
            .map(|f| {
                let (p, q, s) = (pt(&g.mesh, g.mesh.facets[3 * f]), pt(&g.mesh, g.mesh.facets[3 * f + 1]), pt(&g.mesh, g.mesh.facets[3 * f + 2]));
                let (u, w) = ([q[0] - p[0], q[1] - p[1], q[2] - p[2]], [s[0] - p[0], s[1] - p[1], s[2] - p[2]]);
                let c = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
                0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
            })
            .sum::<f64>()
    };
    for name in ["air.xmin", "air.xmax", "air.zmin", "air.zmax"] {
        assert!(rel(area(a(name)), 16.0) < 1e-9, "{name}: {}", area(a(name)));
    }
}

fn mesh_for_sizing_free(m: &Model, plan: MeshPlan) -> (MeshedGeom, Option<femcore::adapt::drive::AdaptReport>) {
    femcore::meshing::mesh_evaluated_report(evaluate(&m.list).unwrap(), &plan).unwrap()
}

#[test]
fn the_adaptive_mesher_follows_the_waves_of_a_study_in_3d() {
    let m = ball_in_box(0.9);
    let plan = MeshPlan::adaptive(0.8, 1.0);
    let s = waves(2, &[("ball", 0.9), ("air", 5.0)]);
    let t0 = std::time::Instant::now();
    let (g, report) = mesh_for_sizing(&m.list, &plan, &s, &PhysicsSizing::default()).unwrap();
    println!("wave adaptation: {} cells in {:.1} s", g.mesh.n_cells(), t0.elapsed().as_secs_f64());
    assert!(report.is_some());
    assert_conforming(&g.mesh);
    let (air, ball) = (g.domain("air").unwrap(), g.domain("ball").unwrap());
    // 0.9 / 6 = 0.15 inside the ball, graded out into the air
    let (lb, la) = (longest_edges(&g.mesh, ball), longest_edges(&g.mesh, air));
    assert!(mean(&lb) < 0.25, "ball mean edge {}", mean(&lb));
    assert!(mean(&la) > 1.5 * mean(&lb), "air {} ball {}", mean(&la), mean(&lb));
    let v = volumes(&g.mesh);
    println!("total volume error {:e}", rel(v[&air] + v[&ball], 64.0));
    assert!(rel(v[&air] + v[&ball], 64.0) < 1e-9);
}

#[test]
fn adaptation_is_deterministic_for_any_thread_count() {
    let m = ball_in_box(0.9);
    let run = || {
        let (g, _) = mesh_for_sizing_free(&m, MeshPlan::adaptive(0.7, 0.6));
        (g.mesh.points.clone(), g.mesh.cells.clone(), g.mesh.facets.clone(), g.mesh.facet_labels.clone())
    };
    let a = run();
    let b = run();
    assert!(a == b, "two runs of the same adaptation differ");
}

// ---------------------------------------------------------------- structured absorbing layers (ADR 52)

/// Air box [-2, 2]^3 with a window [-1.5, 1.5]^3 and absorbing layers of 0.5 (two layers) on every side; a glass block inside, so the
/// physical window holds two materials.
fn boxed_pml() -> Model {
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let glass = m.box_("glass", (-2.0, -2.0, -2.0), (2.0, 2.0, 0.0));
    let win = m.box_("window", (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5));
    m.overlay(&[("air", air), ("glass", glass)]);
    m.pml("pml", win, &[], 0.5, 2.0, 2.0);
    m
}

fn structured(plan: MeshPlan) -> MeshPlan {
    MeshPlan { pml: PmlMode::Structured, ..plan }
}

fn inline_pml(plan: MeshPlan) -> MeshPlan {
    MeshPlan { pml: PmlMode::Inline, ..plan }
}

#[test]
fn structured_layers_have_the_labels_volumes_and_walls_of_the_inline_ones() {
    let m = boxed_pml();
    let inline = mesh_model(&m.list, &inline_pml(MeshPlan::simple(0.5))).unwrap();
    let st = mesh_model(&m.list, &structured(MeshPlan::simple(0.5))).unwrap();
    // Default `MeshPlan::simple` (`PmlMode::Auto`) selects structured layers automatically
    let auto_st = mesh_model(&m.list, &MeshPlan::simple(0.5)).unwrap();
    assert_eq!(auto_st.mesh.cells, st.mesh.cells);
    assert_eq!(st.labels, inline.labels, "the same names mean the same numbers");
    assert_conforming(&st.mesh);
    let (vi, vs) = (volumes(&inline.mesh), volumes(&st.mesh));
    assert_eq!(vi.keys().collect::<Vec<_>>(), vs.keys().collect::<Vec<_>>());
    for (l, v) in &vi {
        assert!((v - vs[l]).abs() < 1e-9 * 64.0, "{}: inline {v} structured {}", st.labels.name_of_domain(*l).unwrap(), vs[l]);
    }
    // the outer wall is the six faces of the box, each with the label of its geom face
    let outer: Vec<i32> = ["air.xmin", "air.xmax", "air.ymin", "air.ymax", "air.zmin", "air.zmax"].iter().flat_map(|n| st.boundary(n)).collect();
    assert!(outer.len() >= 6);
    for l in &outer {
        assert!((0..st.mesh.n_facets()).any(|f| st.mesh.facet_labels[f] == *l), "no facet carries the boundary {l}");
    }
    assert!(!st.domains("air*").is_empty() && st.pml_materials(3.0).len() == inline.pml_materials(3.0).len());
}

#[test]
fn structured_layers_keep_the_physical_vertex_and_cell_numbers_and_have_the_asked_layers() {
    let m = boxed_pml();
    let ev = evaluate(&m.list).unwrap();
    let opts = femcore::meshing::geomdrive::GeomMeshOpts::new(0.5);
    // the physical window alone
    let mut skip = vec![false; ev.regions().len()];
    for p in ev.pml_regions() {
        skip[p.region as usize] = true;
    }
    let (phys, _, _) = mesh_geom_3d_with(evaluate(&m.list).unwrap(), &opts, &Sizes3 { region_h: vec![], skip }).unwrap();
    let (np, nc) = (phys.mesh.n_points(), phys.mesh.n_cells());
    let Native3::Cells(cells, report) = structured_pml_3d_native(phys.clone()).unwrap() else { panic!("the model has a Pml") };
    assert_eq!(report.layers3, [[2, 2], [2, 2], [2, 2]]);
    assert_eq!((report.physical_cells, report.dropped_cells), (nc, 0));
    // vertices and cells of the physical mesh come first, in the same ids, order, connectivity and label
    for v in 0..np {
        assert_eq!(cells.points[v], pt(&phys.mesh, v as u32), "vertex {v}");
    }
    for c in 0..nc {
        assert_eq!(cells.cell(c), phys.mesh.cell(c), "cell {c}");
        assert_eq!(cells.label[c], phys.mesh.cell_labels[c], "label of cell {c}");
    }
    // the layers: 6 faces + 12 edges + 8 corners of pieces, uncut cells
    use femcore::mesh::cells::CellType;
    let counts = cells.counts();
    assert!(counts[&CellType::Prism] > 0 && counts[&CellType::Hex] > 0);
    // the whole mesh is the box [-2, 2]^3
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &cells.points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    assert!(lo.iter().all(|x| (x + 2.0).abs() < 1e-12) && hi.iter().all(|x| (x - 2.0).abs() < 1e-12), "{lo:?} {hi:?}");
    let total: f64 = (0..cells.n_cells()).map(|c| cells.measure(c)).sum();
    assert!((total - 64.0).abs() < 1e-9, "{total}");
    assert!((0..cells.n_cells()).all(|c| cells.measure(c) > 0.0));
    femcore::mesh::mixed::MixedDiscretisation::from_cell_mesh(&cells)
        .unwrap()
        .validate_conforming()
        .unwrap();
}

#[test]
fn structured_layers_follow_an_adapted_physical_mesh() {
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let ball = m.sphere("ball", (0.0, 0.0, 0.0), 0.6);
    let win = m.box_("window", (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5));
    m.overlay(&[("air", air), ("ball", ball)]);
    m.pml("pml", win, &[], 0.5, 2.0, 2.0);
    let plan = structured(MeshPlan::adaptive(0.6, 0.6));
    let g = mesh_model(&m.list, &plan).unwrap();
    assert_conforming(&g.mesh);
    let v = volumes(&g.mesh);
    assert!(rel(v.values().sum::<f64>(), 64.0) < 1e-9, "{v:?}");
    let ball_l = g.domain("ball").unwrap();
    let err = rel(v[&ball_l], 4.0 / 3.0 * PI * 0.6f64.powi(3));
    println!("ball volume error {err:.3e}, {} cells", g.mesh.n_cells());
    assert!(err < 0.15, "{err}");
    // the vertices of the ball surface are on the sphere (the oracle), whatever the coarse sampling does to the volume
    let iface = g.boundary("air|ball")[0];
    for f in 0..g.mesh.n_facets() {
        if g.mesh.facet_labels[f] != iface {
            continue;
        }
        for k in 0..3 {
            let p = pt(&g.mesh, g.mesh.facets[3 * f + k]);
            let d = ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - 0.6).abs();
            assert!(d < 1e-6, "vertex {p:?} is {d:e} off the sphere, boundary {iface} = {:?}", g.labels.name_of_boundary(iface));
        }
    }
    // the layers got their domain numbers from geom
    assert!(g.domain("air.pml_xmax").is_some() && v.contains_key(&g.domain("air.pml_xmax").unwrap()));
}

#[test]
fn layers_can_replace_the_inline_pieces_of_a_mesh_made_the_ordinary_way() {
    // the same entry point as in 2D: mesh the pieces inline, then replace them. The window planes are surfaces of the model, so the lattice
    // mesher cuts them exactly even where they are off the lattice, and the physical mesh is box-shaped
    let m = boxed_pml();
    let inline = mesh_model(&m.list, &inline_pml(MeshPlan::simple(0.37))).unwrap();
    let (st, report) = femcore::meshing::pml_geom::structured_pml(inline.clone()).unwrap();
    assert!(report.dropped_cells > 0 && report.layer_quads > 0);
    assert_eq!(report.layers3, [[2, 2], [2, 2], [2, 2]]);
    let (vi, vs) = (volumes(&inline.mesh), volumes(&st.mesh));
    for (l, v) in &vi {
        assert!((v - vs[l]).abs() < 1e-9 * 64.0, "label {l}: inline {v} structured {}", vs[l]);
    }
    assert_conforming(&st.mesh);
}

// ---------------------------------------------------------------- the same mesh for any thread count (ADR 11)

fn mesh_hash(plan: MeshPlan) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let ball = m.sphere("ball", (0.0, 0.0, 0.0), 0.7);
    let win = m.box_("window", (-1.5, -1.5, -1.5), (1.5, 1.5, 1.5));
    m.overlay(&[("air", air), ("ball", ball)]);
    m.pml("pml", win, &[], 0.5, 2.0, 2.0);
    let g = mesh_model(&m.list, &plan).unwrap();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    g.mesh.points.iter().for_each(|x| x.to_bits().hash(&mut h));
    g.mesh.cells.hash(&mut h);
    g.mesh.cell_labels.hash(&mut h);
    g.mesh.facets.hash(&mut h);
    g.mesh.facet_labels.hash(&mut h);
    h.finish()
}

/// Run by `meshes_are_identical_for_any_thread_count` in child processes with `NYX_PROBE` set; does nothing otherwise.
#[test]
fn probe_mesh_hashes() {
    if std::env::var_os("NYX_PROBE").is_none() {
        return;
    }
    let plans = [MeshPlan::simple(0.5), MeshPlan::adaptive(0.6, 0.8), structured(MeshPlan::simple(0.5)), structured(MeshPlan::adaptive(0.6, 0.8))];
    let hashes: Vec<String> = plans.iter().map(|p| format!("{:016x}", mesh_hash(p.clone()))).collect();
    println!("HASHES={}", hashes.join(","));
}

#[test]
fn meshes_are_identical_for_any_thread_count() {
    let exe = std::env::current_exe().unwrap();
    let run = |threads: &str| -> String {
        let out = std::process::Command::new(&exe)
            .args(["probe_mesh_hashes", "--exact", "--nocapture", "--test-threads=1"])
            .env("NYX_PROBE", "1")
            .env("NYX_THREADS", threads)
            .output()
            .unwrap();
        let text = String::from_utf8_lossy(&out.stdout).to_string();
        assert!(out.status.success(), "{text}\n{}", String::from_utf8_lossy(&out.stderr));
        text.lines().find_map(|l| l.find("HASHES=").map(|i| l[i + 7..].trim().to_string())).unwrap_or_else(|| panic!("no hashes in {text}"))
    };
    let one = run("1");
    assert_eq!(one, run("3"), "the mesh depends on the number of threads");
    assert_eq!(one, run("8"), "the mesh depends on the number of threads");
}
