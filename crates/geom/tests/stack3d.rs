//! 3D prismatic models: volumes, named faces, closure, layer stacks, periodic pairs, Pml and Cut in 3D, refusals.

use geom::diff::{diff, ChangeKind};
use geom::eval::{Evaluator, Tier};
use geom::model::*;
use geom::select::Selector;
use geom::sim::EntityRef;
use std::f64::consts::PI;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}

#[test]
fn box_with_a_through_cylinder_cavity() {
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
    let c = m.cylinder("cav", (1.0, 1.0, -1.0), 0.3, 3.0);
    let d = m.difference("d", b, &[c]);
    m.overlay(&[("block", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 4.0 - PI * 0.09));
    let mut names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names, ["block.xmax", "block.xmin", "block.ymax", "block.ymin", "block.zmax", "block.zmin", "cav.side"]);
    let area: f64 = ev.faces().iter().map(|f| f.area).sum();
    assert!(close(area, 2.0 * (4.0 - PI * 0.09) + 8.0 + 2.0 * PI * 0.3));
    // selectors work on faces
    let side = ev.resolve(&Selector::parse("name(cav.side)").unwrap()).unwrap();
    assert_eq!(side.set.kind, geom::select::EntityKind::Faces);
    assert_eq!(ev.resolve(&Selector::parse("kind(cylinder)").unwrap()).unwrap().set, side.set);
    assert_eq!(ev.resolve(&Selector::parse("radius(0.3)").unwrap()).unwrap().set, side.set);
    assert_eq!(ev.resolve(&Selector::parse("outer & extreme(+z)").unwrap()).unwrap().set.ids.len(), 1);
    assert_eq!(ev.resolve(&Selector::parse("plane(z = 0)").unwrap()).unwrap().set.ids.len(), 1);
    // classification and the projection oracle
    assert_eq!(ev.classify([0.2, 0.2, 0.5]), Some(0));
    assert_eq!(ev.classify([1.0, 1.0, 0.5]), None);
    let f = side.set.ids[0];
    let q = ev.oracle().project([2.0, 1.0, 0.4], EntityRef::Face(f)).unwrap();
    assert!((q[0] - 1.3).abs() < 1e-12 && (q[1] - 1.0).abs() < 1e-12 && (q[2] - 0.4).abs() < 1e-12);
    let lm = ev.label_map().unwrap();
    assert_eq!(lm.boundary("cav.side").len(), 1);
}

fn bragg(duty: f64, tooth: f64) -> Model {
    let mut m = Model::new_3d();
    m.param("w", 500e-9)
        .param("t", 300e-9)
        .param("period", 300e-9)
        .param("duty", duty)
        .param("tooth", tooth)
        .param("margin", 0.45e-6)
        .param("box_t", 2e-6);
    let e = "(tooth+margin)";
    let clad = m.box_(
        "clad",
        (format!("-w/2-{e}").as_str(), "-box_t", 0.0),
        (format!("w/2+{e}").as_str(), "t+margin", "period"),
    );
    let bx = m.box_("box", (format!("-w/2-{e}").as_str(), "-box_t", 0.0), (format!("w/2+{e}").as_str(), 0.0, "period"));
    let core = m.box_("core", ("-w/2", 0.0, 0.0), ("w/2", "t", "period"));
    let lt = m.box_("left_tooth", ("-w/2-tooth", 0.0, 0.0), ("-w/2", "t", "duty*period"));
    let rt = m.box_("right_tooth", ("w/2", 0.0, 0.0), ("w/2+tooth", "t", "duty*period"));
    m.overlay_defs(
        vec![
            RegionDef { name: "clad".into(), input: clad, kind: RegionKind::Material, pinned_id: Some(3) },
            RegionDef { name: "BOX".into(), input: bx, kind: RegionKind::Material, pinned_id: Some(1) },
            RegionDef { name: "Si".into(), input: core, kind: RegionKind::Material, pinned_id: Some(2) },
            RegionDef { name: "left_tooth".into(), input: lt, kind: RegionKind::Material, pinned_id: Some(4) },
            RegionDef { name: "right_tooth".into(), input: rt, kind: RegionKind::Material, pinned_id: Some(5) },
        ],
        vec![],
    );
    m.periodic(
        "floquet",
        Selector::parse("plane(z = 0)").unwrap(),
        Selector::parse("plane(z = period)").unwrap(),
        (0.0, 0.0, "period"),
    );
    m
}

#[test]
fn bragg_cell_is_a_layer_stack_with_periodic_ends() {
    let m = bragg(0.5, 80e-9);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let (w, t, p, d, tooth, mg, bt) = (500e-9, 300e-9, 300e-9, 0.5, 80e-9, 0.45e-6, 2e-6);
    let wide = w + 2.0 * (tooth + mg);
    let vol = |n: &str| ev.volume(ev.regions().iter().position(|r| r.name == n).unwrap() as u32);
    assert!(close(vol("Si"), w * t * p));
    assert!(close(vol("left_tooth"), tooth * t * d * p) && close(vol("right_tooth"), tooth * t * d * p));
    assert!(close(vol("BOX"), wide * bt * p));
    assert!(close(vol("clad"), wide * (bt + t + mg) * p - wide * bt * p - w * t * p - 2.0 * tooth * t * d * p));
    let st = ev.layer_stack().expect("a stack");
    assert_eq!(st.layers.len(), 2);
    assert!(close(st.layers[0].z1, d * p) && close(st.layers[1].z1, p));
    // the overlay is the common refinement: it has the tooth cells even though the second layer has none
    assert!(st.overlay_region_count() > ev.regions().len() - 2);
    for l in &st.layers {
        assert_eq!(l.region_of.len(), st.overlay_region_count());
    }
    // identical end faces: one Periodic pair, verified congruent
    let pp = &ev.periodic_pairs()[0];
    assert!(pp.checked_congruent && !pp.a.ids.is_empty() && pp.a.ids.len() == pp.b.ids.len());
    // the tooth tops are named after the tooth; the interface between a tooth and the cladding above it is `clad|left_tooth`
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert!(names.contains(&"left_tooth|clad") || names.contains(&"clad|left_tooth"), "{names:?}");
    assert!(names.contains(&"clad.xmin") && names.contains(&"box.zmin") && names.contains(&"clad.zmax"));
    // an unequal tooth makes the same topology: ids stay put while the geometry moves
    let mut e = Evaluator::new();
    let a = e.eval(&m.list, Tier::Complex).unwrap();
    let b = e.eval(&m.list.with(&[("tooth", 90e-9)]), Tier::Complex).unwrap();
    assert_eq!(a.topology_signature(), b.topology_signature());
    assert_eq!(
        a.faces().iter().map(|f| &f.name).collect::<Vec<_>>(),
        b.faces().iter().map(|f| &f.name).collect::<Vec<_>>()
    );
    assert_eq!(diff(&a, &b).kind, ChangeKind::Morph);
    // z-only change: the duty cycle moves the layer boundary, the overlay stays
    let c = e.eval(&m.list.with(&[("duty", 0.4)]), Tier::Complex).unwrap();
    assert_eq!(diff(&a, &c).kind, ChangeKind::ZOnly);
    // removing the teeth changes the topology
    let mut no_teeth = m.list.clone();
    for f in &mut no_teeth.features {
        if f.name.ends_with("_tooth") {
            f.suppressed = true;
        }
    }
    let nt = e.eval(&no_teeth, Tier::Complex).unwrap();
    assert_eq!(diff(&a, &nt).kind, ChangeKind::Remesh);
    assert_eq!(nt.layer_stack().unwrap().layers.len(), 1);
}

#[test]
fn a_wrong_periodic_translation_is_refused() {
    let mut m = bragg(0.5, 80e-9);
    for f in &mut m.list.features {
        if let Op::Periodic { translation, .. } = &mut f.op {
            translation[2] = geom::expr::Expr::Num(250e-9);
        }
    }
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("Periodic 'floquet'"));
}

#[test]
fn extruded_profiles_pml_and_cuts_in_3d() {
    // a rounded plate (an x-z profile) extruded 1 long along y, with PML on the y sides of a window
    let mut m = Model::new_3d();
    let p = m.rounded_rect("outline", (-2.0, -1.0), (2.0, 1.0), 0.5);
    let slab = m.extrude("slab", p, 0.0, 1.0);
    let win = m.box_("window", (-2.0, 0.0, -1.0), (2.0, 1.0, 1.0));
    m.overlay(&[("slab", slab)]);
    m.pml("pml", win, &[Side::YMin, Side::YMax], 0.5, 4.0, 2.0);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), (8.0 - (4.0 - PI) * 0.25) * 1.0));
    // the slab only exists inside the window, so no pad overlaps it: no pml pieces appear
    assert_eq!(ev.regions().len(), 1);
    let mut m = Model::new_3d();
    let a = m.box_("a", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("a", a)]);
    m.cut("sym", (1.0, 0.0, 0.0), 0.5, Keep::Below);
    m.cut("floor", (0.0, 0.0, 1.0), 0.25, Keep::Above);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 0.5 * 1.0 * 0.75));
    assert!(ev.faces().iter().any(|f| f.name == "sym.cut") && ev.faces().iter().any(|f| f.name == "floor.cut"));
    assert_eq!(ev.resolve(&Selector::Tag("symmetry".into())).unwrap().set.ids.len(), 2);
    // a slanted cut of a curved prism: refused with an explanation (a slanted cut of a polyhedron works, see tests/planar.rs)
    let mut m = Model::new_3d();
    let a = m.cylinder("pin", (0.0, 0.0, 0.0), 0.5, 1.0);
    m.overlay(&[("a", a)]);
    m.cut("slant", (1.0, 0.0, 1.0), 0.5, Keep::Below);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("slanted"));
    // PML on all sides of a window around a box: pieces appear for the regions that reach the pad
    let mut m = Model::new_3d();
    let big = m.box_("big", (-1.0, -1.0, -1.0), (2.0, 2.0, 2.0));
    let win = m.box_("window", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("air", big)]);
    m.pml("pml", win, &[], 1.0, 4.0, 2.0);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert_eq!(ev.regions().len(), 27); // the physical core and 26 pad pieces
    assert_eq!(ev.pml_regions().len(), 26);
    assert!(close(ev.regions().iter().map(|r| r.measure).sum::<f64>(), 27.0));
    let corner = ev.pml_regions().iter().find(|p| p.sides == [1, 1, 1]).unwrap();
    assert_eq!(ev.regions()[corner.region as usize].name, "air.pml_xmax_ymax_zmax");
}

#[test]
fn three_d_errors_are_plain() {
    let mut m = Model::new_3d();
    let r = m.rect("flat", (0.0, 0.0), (1.0, 1.0));
    m.overlay(&[("r", r)]);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("must be the input of an Extrude"));
    let mut m = Model::new_2d();
    let b = m.box_("b", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("b", b)]);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("3D solid"));
    let mut m = Model::new_3d();
    let a = m.box_("a", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let b = m.box_("b", (0.0, 0.0, 1.0 + 3e-8), (1.0, 1.0, 2.0)); // z levels 3e-8 apart: undecidable
    m.overlay(&[("a", a), ("b", b)]);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("z levels"));
}

#[test]
fn faces_round_trip_through_fractions() {
    let m = bragg(0.5, 80e-9);
    let ev = m.eval(Tier::Complex).unwrap();
    let l = ev.bounds().diagonal();
    for (i, f) in ev.faces().iter().enumerate() {
        for (u, v) in [(0.0, 0.0), (0.3, 0.7), (1.0, 1.0), (0.5, 0.5)] {
            let p = ev.position_at(EntityRef::Face(i as u32), geom::query::Fractions { u, v }).unwrap();
            let g = ev.fractions_of(EntityRef::Face(i as u32), p).unwrap();
            let q = ev.position_at(EntityRef::Face(i as u32), g).unwrap();
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            assert!(d <= 1e-12 * l, "{} at ({u},{v}): {d}", f.name);
        }
    }
}

#[test]
fn tessellation_and_stuffing_surfaces() {
    let m = bragg(0.5, 80e-9);
    let ev = m.eval(Tier::Display).unwrap();
    let mesh = ev.tessellate(1e-9);
    // triangle area per face matches the analytic face area (all faces here are planar)
    let mut by = vec![0.0; ev.faces().len()];
    for (t, &f) in mesh.tris.iter().zip(&mesh.tri_face) {
        let p = |i: u32| mesh.points[i as usize];
        let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        by[f as usize] += 0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
    }
    for (f, a) in ev.faces().iter().zip(by) {
        assert!((a - f.area).abs() <= 1e-9 * f.area, "{}: {a} vs {}", f.name, f.area);
    }
    // surfaces: every face is on exactly one surface, coincident faces share it, and the extent covers the faces
    let ss = ev.surfaces_for_stuffing();
    let total: usize = ss.iter().map(|s| s.faces.len()).sum();
    assert_eq!(total, ev.faces().len());
    assert!(ss.len() < ev.faces().len());
    for s in &ss {
        for &f in &s.faces {
            let b = ev.faces()[f as usize].bbox;
            assert!((0..3).all(|k| s.face_aabb.lo[k] <= b.lo[k] && s.face_aabb.hi[k] >= b.hi[k]));
        }
    }
    // a cylinder cavity gives a cylinder surface
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
    let c = m.cylinder("cav", (1.0, 1.0, -1.0), 0.3, 3.0);
    let d = m.difference("d", b, &[c]);
    m.overlay(&[("block", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev
        .surfaces_for_stuffing()
        .iter()
        .any(|s| matches!(s.surface, geom::stack::Surface::Cylinder { r, .. } if (r - 0.3).abs() < 1e-15)));
}

#[test]
fn three_d_edges_and_vertices_form_a_closed_topology() {
    // a plain box: 12 edges, 8 vertices, every edge shared by two faces
    let mut m = Model::new_3d();
    let b = m.box_("b", (0.0, 0.0, 0.0), (2.0, 3.0, 4.0));
    m.overlay(&[("b", b)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert_eq!((ev.faces().len(), ev.edges3().len(), ev.vertices().len()), (6, 12, 8));
    assert!(ev.edges3().iter().all(|e| e.faces.len() == 2));
    // Euler: V - E + F = 2 for a closed convex shell
    assert_eq!(ev.vertices().len() as i64 - ev.edges3().len() as i64 + ev.faces().len() as i64, 2);
    assert!(
        ev.edges3().iter().any(|e| e.name == "b.xmin/b.ymin") && ev.edges3().iter().any(|e| e.name == "b.xmax/b.zmax")
    );
    let total: f64 = ev.edges3().iter().map(|e| e.length).sum();
    assert!(close(total, 4.0 * (2.0 + 3.0 + 4.0)));
    // selectors, oracle and fractions on edges and vertices
    let e = ev.resolve(&Selector::parse("name(b.xmin/b.ymin)").unwrap()).unwrap();
    assert_eq!(e.set.kind, geom::select::EntityKind::Edges);
    let id = e.set.ids[0];
    let q = ev.oracle().project([-1.0, -1.0, 1.0], EntityRef::Edge(id)).unwrap();
    assert!((q[0]).abs() < 1e-15 && (q[1]).abs() < 1e-15 && (q[2] - 1.0).abs() < 1e-15);
    let f = ev.fractions_of(EntityRef::Edge(id), [0.0, 0.0, 1.0]).unwrap();
    assert!((f.u - 0.25).abs() < 1e-15);
    let p = ev.position_at(EntityRef::Edge(id), f).unwrap();
    assert!((p[2] - 1.0).abs() < 1e-15);
    // a block with a cylindrical cavity: the circles are edges between the cavity wall and the two caps, the seam is not an edge
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
    let c = m.cylinder("cav", (1.0, 1.0, -1.0), 0.3, 3.0);
    let d = m.difference("d", b, &[c]);
    m.overlay(&[("block", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert_eq!(ev.edges3().len(), 12 + 2);
    assert!(ev.edges3().iter().filter(|e| e.name.contains("cav.side")).count() == 2);
    // a Bragg period: junction edges are shared by three faces (two regions meeting the void, or three regions)
    let ev = bragg(0.5, 80e-9).eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(ev.edges3().iter().any(|e| e.faces.len() >= 3));
    assert!(ev.edges3().iter().all(|e| e.faces.len() >= 2));
    // vertices know their edges and are named after the faces around them
    assert!(ev.vertices().iter().all(|v| !v.edges.is_empty() && !v.name.is_empty()));
}

#[test]
fn curved_bodies_inside_a_host_are_exact_and_conflicts_are_refused() {
    let v_sphere = |r: f64| 4.0 / 3.0 * PI * r.powi(3);
    // sphere in a box: two regions, one extra face, no edges on it
    let mut m = Model::new_3d();
    let air = m.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let ball = m.sphere("ball", (2.0, 2.0, 2.0), 0.7);
    m.overlay(&[("air", air), ("ball", ball)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 64.0 - v_sphere(0.7)) && close(ev.volume(1), v_sphere(0.7)));
    assert_eq!(ev.faces().len(), 7);
    assert_eq!(ev.edges3().len(), 12);
    let f = ev.resolve(&Selector::parse("name(air|ball)").unwrap()).unwrap().set;
    assert_eq!(f.ids.len(), 1);
    assert_eq!(ev.resolve(&Selector::parse("kind(sphere)").unwrap()).unwrap().set, f);
    assert_eq!(ev.resolve(&Selector::parse("radius(0.7)").unwrap()).unwrap().set, f);
    assert_eq!(ev.classify([2.0, 2.0, 2.0]), Some(1));
    assert_eq!(ev.classify([0.5, 0.5, 0.5]), Some(0));
    assert!(close(ev.faces()[f.ids[0] as usize].area, 4.0 * PI * 0.49));
    let q = ev.oracle().project([4.0, 2.0, 2.0], EntityRef::Face(f.ids[0])).unwrap();
    assert!((q[0] - 2.7).abs() < 1e-14 && (q[1] - 2.0).abs() < 1e-14);
    let l = ev.bounds().diagonal();
    for (u, v) in [(0.1, 0.2), (0.9, 0.8), (0.5, 0.5)] {
        let p = ev.position_at(EntityRef::Face(f.ids[0]), geom::query::Fractions { u, v }).unwrap();
        let g = ev.fractions_of(EntityRef::Face(f.ids[0]), p).unwrap();
        let p2 = ev.position_at(EntityRef::Face(f.ids[0]), g).unwrap();
        assert!(norm(p, p2) <= 1e-12 * l);
    }
    let mesh = ev.tessellate(1e-3);
    let area: f64 = mesh
        .tris
        .iter()
        .zip(&mesh.tri_face)
        .filter(|(_, &fc)| fc == f.ids[0])
        .map(|(t, _)| tri_area(&mesh.points, t))
        .sum();
    assert!((area / (4.0 * PI * 0.49) - 1.0).abs() < 1e-2);
    // concentric spheres (nested) and a torus
    let mut m = Model::new_3d();
    let host = m.box_("host", (-3.0, -3.0, -3.0), (3.0, 3.0, 3.0));
    let outer = m.sphere("outer", (0.0, 0.0, 0.0), 2.0);
    let inner = m.sphere("inner", (0.0, 0.0, 0.0), 1.0);
    let tor = m.torus("ring", (0.0, 0.0, 2.6), 0.2, 0.1);
    m.overlay(&[("vacuum", host), ("shell", outer), ("core", inner), ("ring", tor)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let vol = |n: &str| ev.volume(ev.regions().iter().position(|r| r.name == n).unwrap() as u32);
    assert!(close(vol("core"), v_sphere(1.0)) && close(vol("shell"), v_sphere(2.0) - v_sphere(1.0)));
    assert!(close(vol("ring"), 2.0 * PI * PI * 0.2 * 0.01));
    assert!(close(vol("vacuum"), 216.0 - v_sphere(2.0) - 2.0 * PI * PI * 0.2 * 0.01));
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert!(names.contains(&"core|shell") && names.contains(&"shell|vacuum") && names.contains(&"ring|vacuum"));
    // refusals name both features
    let refuse = |build: &dyn Fn(&mut Model) -> Vec<(&'static str, FeatureId)>| {
        let mut m = Model::new_3d();
        let regs = build(&mut m);
        m.overlay(&regs);
        m.eval(Tier::Complex).err().expect("expected a refusal").0
    };
    let e = refuse(&|m| {
        let b = m.box_("box", (0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
        let s = m.sphere("ball", (1.0, 1.0, 1.9), 0.5); // pokes through the top
        vec![("a", b), ("b", s)]
    });
    assert!(e.contains("sphere 'ball'") && e.contains("box.zmax"), "{e}");
    let e = refuse(&|m| {
        let b = m.box_("box", (0.0, 0.0, 0.0), (6.0, 6.0, 6.0));
        let s1 = m.sphere("s1", (2.0, 3.0, 3.0), 1.0);
        let s2 = m.sphere("s2", (3.0, 3.0, 3.0), 1.0); // overlapping spheres
        vec![("a", b), ("b", s1), ("c", s2)]
    });
    assert!(e.contains("'s1'") && e.contains("'s2'") && e.contains("intersect"), "{e}");
    let e = refuse(&|m| {
        let s = m.sphere("alone", (0.0, 0.0, 0.0), 1.0);
        vec![("a", s)]
    });
    assert!(e.contains("prismatic host"), "{e}");
    // printing and the parameter model keep working
    let mut m = Model::new_3d();
    let host = m.box_("host", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    m.param("r", 0.5);
    let s = m.sphere("ball", (2.0, 2.0, 2.0), "r");
    m.overlay(&[("host", host), ("ball", s)]);
    assert_eq!(common::interp::build(&m.list.to_rust()).list, m.list);
    let mut e = Evaluator::new();
    let a = e.eval(&m.list, Tier::Complex).unwrap();
    let b = e.eval(&m.list.with(&[("r", 0.9)]), Tier::Complex).unwrap();
    assert_eq!(a.topology_signature(), b.topology_signature());
    assert_eq!(diff(&a, &b).kind, ChangeKind::Morph);
}

mod common;

fn norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn tri_area(p: &[[f64; 3]], t: &[u32; 3]) -> f64 {
    let (a, b, c) = (p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]);
    let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

#[test]
fn cylinders_and_extrusions_along_x_and_y() {
    // a block with a cavity drilled along each axis: volume, names and projection agree with the z version
    let build = |axis: Axis| {
        let mut m = Model::new_3d();
        let b = m.box_("block", (0.0, 0.0, 0.0), (3.0, 4.0, 5.0));
        let base = match axis {
            Axis::X => (-1.0, 2.0, 2.5),
            Axis::Y => (1.5, -1.0, 2.5),
            Axis::Z => (1.5, 2.0, -1.0),
        };
        let len = match axis {
            Axis::X => 5.0,
            Axis::Y => 6.0,
            Axis::Z => 7.0,
        };
        let c = m.cylinder_along("cav", axis, base, 0.4, len);
        let d = m.difference("d", b, &[c]);
        m.overlay(&[("block", d)]);
        m.eval(Tier::Complex).unwrap()
    };
    let extent = [3.0, 4.0, 5.0];
    for (k, axis) in [Axis::X, Axis::Y, Axis::Z].into_iter().enumerate() {
        let ev = build(axis);
        ev.complex().check().unwrap();
        assert!(close(ev.volume(0), 60.0 - PI * 0.16 * extent[k]), "{axis:?}");
        let mut names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(
            names,
            ["block.xmax", "block.xmin", "block.ymax", "block.ymin", "block.zmax", "block.zmin", "cav.side"],
            "{axis:?}"
        );
        assert_eq!((ev.faces().len(), ev.edges3().len()), (7, 14), "{axis:?}");
        let b = ev.bounds();
        // the cutter is longer than the block along its own axis (5, 6, 7) and the block fixes the other two extents
        let want = [[5.0, 4.0, 5.0], [3.0, 6.0, 5.0], [3.0, 4.0, 7.0]][k];
        assert!(
            close(b.hi[0] - b.lo[0], want[0]) && close(b.hi[1] - b.lo[1], want[1]) && close(b.hi[2] - b.lo[2], want[2]),
            "bounds in world axes {axis:?}: {b:?}"
        );
        // the cavity wall is a cylinder about the right axis, 0.4 from the axis line
        let side = ev.resolve(&Selector::parse("name(cav.side)").unwrap()).unwrap().set.ids[0];
        let probe = [[2.9, 2.0, 2.5], [1.5, 3.9, 2.5], [1.5, 2.0, 4.9]][k];
        let q = ev.oracle().project(probe, EntityRef::Face(side)).unwrap();
        let along = k;
        let off: Vec<usize> = (0..3).filter(|&i| i != along).collect();
        let centre = [[0.0, 2.0, 2.5], [1.5, 0.0, 2.5], [1.5, 2.0, 0.0]][k];
        let r = ((q[off[0]] - centre[off[0]]).powi(2) + (q[off[1]] - centre[off[1]]).powi(2)).sqrt();
        assert!((r - 0.4).abs() < 1e-12 && (q[along] - probe[along]).abs() < 1e-12, "{axis:?}: {q:?}");
        // a point inside the cavity is outside the solid, a point beside it inside
        let mut inside = centre;
        inside[along] = probe[along];
        assert_eq!(ev.classify(inside), None, "{axis:?}");
        let mut beside = inside;
        beside[off[0]] += 0.6;
        assert_eq!(ev.classify(beside), Some(0), "{axis:?}");
        // the stuffing surface is a cylinder with the right axis
        let cyl = ev.surfaces_for_stuffing().into_iter().find_map(|s| match s.surface {
            geom::stack::Surface::Cylinder { axis, r, .. } => Some((axis, r)),
            _ => None,
        });
        let (ax, r) = cyl.unwrap();
        assert!((ax[k].abs() - 1.0).abs() < 1e-12 && (r - 0.4).abs() < 1e-15, "{axis:?}: {ax:?}");
        // fractions and edges in world coordinates
        let l = ev.bounds().diagonal();
        for (u, v) in [(0.2, 0.3), (0.9, 0.6)] {
            let p = ev.position_at(EntityRef::Face(side), geom::query::Fractions { u, v }).unwrap();
            let g = ev.fractions_of(EntityRef::Face(side), p).unwrap();
            let p2 = ev.position_at(EntityRef::Face(side), g).unwrap();
            assert!(norm(p, p2) <= 1e-12 * l, "{axis:?}");
        }
        // the circle edges lie in planes perpendicular to the axis
        for e in ev.edges3().iter().filter(|e| e.name.contains("cav.side")) {
            assert!((e.bbox.hi[k] - e.bbox.lo[k]).abs() < 1e-12, "{axis:?}: edge not flat along its axis");
        }
    }
    // an extrusion along x of a profile given in (y, z); the profile's circle sweeps a cylinder along x
    let mut m = Model::new_3d();
    let prof = m.rect("outline", (0.0, 0.0), (2.0, 3.0));
    let ex = m.extrude_along("bar", Axis::X, prof, 1.0, 5.0);
    m.overlay(&[("bar", ex)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let b = ev.bounds();
    assert!(close(b.lo[0], 1.0) && close(b.hi[0], 5.0) && close(b.hi[1], 2.0) && close(b.hi[2], 3.0));
    assert!(close(ev.volume(0), 24.0));
    assert!(
        ev.faces().iter().any(|f| f.name == "outline.side:left") && ev.faces().iter().any(|f| f.name == "outline.xmax")
    );
    // a Pml around a window in a model swept along y keeps world sides and names
    let mut m = Model::new_3d();
    let big = m.box_("air", (-1.0, -1.0, -1.0), (2.0, 2.0, 2.0));
    let rod = m.cylinder_along("rod", Axis::Y, (0.5, -0.5, 0.5), 0.1, 2.0);
    let win = m.box_("window", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("air", big), ("rod", rod)]);
    m.pml("pml", win, &[], 1.0, 4.0, 2.0);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let corner =
        ev.pml_regions().iter().find(|p| ev.regions()[p.region as usize].name == "air.pml_xmax_ymax_zmax").unwrap();
    assert_eq!(corner.sides, [1, 1, 1]);
    let face = ev.pml_regions().iter().find(|p| ev.regions()[p.region as usize].name == "air.pml_ymin").unwrap();
    assert_eq!(face.sides, [0, -1, 0]);
}

#[test]
fn cone_frustum_bodies() {
    let vol = |r0: f64, r1: f64, h: f64| PI * h * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
    let side_area = |r0: f64, r1: f64, h: f64| PI * (r0 + r1) * h.hypot(r0 - r1);
    // a speaker-like diaphragm: a frustum floating in a box
    let mut m = Model::new_3d();
    let air = m.box_("air", (-3.0, -3.0, -1.0), (3.0, 3.0, 4.0));
    let cone = m.cone("diaphragm", (0.0, 0.0, 0.0), 1.5, 0.5, 2.0);
    m.overlay(&[("air", air), ("membrane", cone)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(1), vol(1.5, 0.5, 2.0)) && close(ev.volume(0), 6.0 * 6.0 * 5.0 - vol(1.5, 0.5, 2.0)));
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names.iter().filter(|n| **n == "air|membrane").count(), 3, "{names:?}");
    assert_eq!((ev.faces().len(), ev.edges3().len()), (6 + 3, 12 + 2));
    let side = ev
        .faces()
        .iter()
        .position(|f| matches!(f.geom, geom::complex::FaceGeom::Quadric(geom::prim::Body::ConeSide { .. })))
        .unwrap();
    assert!(close(ev.faces()[side].area, side_area(1.5, 0.5, 2.0)));
    // exact projection onto the slanted wall: the point (3, 0, 1) projects to the wall at z = 1 where the radius is 1.0
    let q = ev.oracle().project([3.0, 0.0, 1.0], EntityRef::Face(side as u32)).unwrap();
    let want = {
        // foot of the perpendicular on the generator (1.5, 0) - (0.5, 2)
        let (a, d): ([f64; 2], [f64; 2]) = ([1.5, 0.0], [-1.0, 2.0]);
        let t = (((3.0 - a[0]) * d[0] + (1.0 - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1])).clamp(0.0, 1.0);
        [a[0] + t * d[0], a[1] + t * d[1]]
    };
    assert!((q[0] - want[0]).abs() < 1e-14 && q[1].abs() < 1e-14 && (q[2] - want[1]).abs() < 1e-14, "{q:?}");
    assert_eq!(ev.classify([0.0, 0.0, 0.5]), Some(1));
    assert_eq!(ev.classify([1.4, 0.0, 0.1]), Some(1));
    assert_eq!(ev.classify([1.4, 0.0, 1.0]), Some(0));
    // the rim circles are edges named after the two faces they join; their parametrisation round-trips
    let rim = ev
        .edges3()
        .iter()
        .position(|e| matches!(e.geom, geom::complex::EdgeGeom::Circle { .. }) && close(e.length, 2.0 * PI * 1.5))
        .unwrap();
    let f = ev.fractions_of(EntityRef::Edge(rim as u32), [0.0, 1.5, 0.0]).unwrap();
    assert!((f.u - 0.25).abs() < 1e-12);
    let p = ev.position_at(EntityRef::Edge(rim as u32), f).unwrap();
    assert!((p[1] - 1.5).abs() < 1e-12 && p[2].abs() < 1e-12);
    assert!(ev.resolve(&Selector::parse("kind(cone)").unwrap()).unwrap().set.ids == vec![side as u32]);
    // tessellation area of the wall matches the exact area
    let mesh = ev.tessellate(1e-3);
    let area: f64 = mesh
        .tris
        .iter()
        .zip(&mesh.tri_face)
        .filter(|(_, &fc)| fc == side as u32)
        .map(|(t, _)| tri_area(&mesh.points, t))
        .sum();
    assert!((area / side_area(1.5, 0.5, 2.0) - 1.0).abs() < 1e-2);
    // a pointed cone has no top disk and one rim; along x it is the same in world coordinates
    let mut m = Model::new_3d();
    let air = m.box_("air", (-1.0, -3.0, -3.0), (5.0, 3.0, 3.0));
    let c = m.cone_along("horn", Axis::X, (0.0, 0.0, 0.0), 1.0, 0.0, 3.0);
    m.overlay(&[("air", air), ("horn", c)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(1), vol(1.0, 0.0, 3.0)));
    assert_eq!((ev.faces().len(), ev.edges3().len()), (6 + 2, 12 + 1));
    assert_eq!(ev.classify([0.5, 0.0, 0.0]), Some(1));
    assert_eq!(ev.classify([2.8, 0.0, 0.5]), Some(0));
    let disk = ev
        .faces()
        .iter()
        .position(|f| matches!(f.geom, geom::complex::FaceGeom::Quadric(geom::prim::Body::Disk { .. })))
        .unwrap();
    assert!(
        (ev.faces()[disk].bbox.lo[0]).abs() < 1e-15 && (ev.faces()[disk].bbox.hi[0]).abs() < 1e-15,
        "the disk is the plane x = 0"
    );
    // refusals: a cone touching a wall
    let mut m = Model::new_3d();
    let air = m.box_("air", (-3.0, -3.0, 0.0), (3.0, 3.0, 4.0));
    let c = m.cone("base_on_floor", (0.0, 0.0, 0.0), 1.0, 0.5, 1.0);
    m.overlay(&[("air", air), ("c", c)]);
    let e = m.eval(Tier::Complex).err().unwrap().0;
    assert!(e.contains("cone 'base_on_floor'") && e.contains("air.zmin"), "{e}");
    assert_eq!(common::interp::build(&m.list.to_rust()).list, m.list);
}

#[test]
fn revolved_profiles() {
    // a tube (Rect profile off the axis): volume, six faces named after the profile edges, four circle edges
    let mut m = Model::new_3d();
    let air = m.box_("air", (-4.0, -4.0, -2.0), (4.0, 4.0, 6.0));
    let prof = m.rect("section", (1.0, 0.0), (2.0, 3.0));
    let tube = m.revolve("tube", prof, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("steel", tube)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(1), PI * (4.0 - 1.0) * 3.0));
    assert!(close(ev.volume(0), 8.0 * 8.0 * 8.0 - PI * 3.0 * 3.0));
    let kinds: Vec<&str> = ev.faces().iter().filter(|f| f.name == "air|steel").map(|f| f.name.as_str()).collect();
    assert_eq!(kinds.len(), 4); // outer wall, inner wall, bottom ring, top ring
    assert_eq!((ev.faces().len(), ev.edges3().len()), (6 + 4, 12 + 4));
    // the inner wall's normal points towards the axis (out of the metal), the outer wall's away from it
    let inner = ev.faces().iter().position(|f| matches!(f.geom, geom::complex::FaceGeom::Quadric(geom::prim::Body::CylSide { r, .. }) if (r - 1.0).abs() < 1e-12)).unwrap();
    let n_in = ev.oracle().project([0.1, 0.0, 1.5], EntityRef::Face(inner as u32)).unwrap();
    assert!((n_in[0] - 1.0).abs() < 1e-14 && n_in[2] == 1.5);
    assert_eq!(ev.classify([1.5, 0.0, 1.0]), Some(1));
    assert_eq!(ev.classify([0.5, 0.0, 1.0]), Some(0));
    // a bell: stepped profile with a slanted shoulder. Volume by Pappus on the polygon; faces: ring, wall, cone, wall, ring
    let mut m = Model::new_3d();
    let air = m.box_("air", (-5.0, -5.0, -1.0), (5.0, 5.0, 8.0));
    let pts = [(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (1.5, 4.0), (1.5, 6.0), (0.0, 6.0)];
    let bell_p = m.polygon_named("bell_profile", &pts, &["floor", "foot", "shoulder", "neck", "lid", "axis"]);
    let bell = m.revolve("bell", bell_p, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("bronze", bell)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let exact = PI * 9.0 * 2.0 + PI * 2.0 * (9.0 + 4.5 + 2.25) / 3.0 * 1.0 + PI * 2.25 * 2.0;
    // frustum between z = 2 and 4 has radii 3 -> 1.5 over height 2: pi h (r0^2 + r0 r1 + r1^2)/3 with h = 2
    let exact = exact - PI * 2.0 * (9.0 + 4.5 + 2.25) / 3.0 * 1.0 + PI * 2.0 * (9.0 + 4.5 + 2.25) / 3.0;
    assert!(close(ev.volume(1), exact), "{} vs {exact}", ev.volume(1));
    let mut names: Vec<String> = ev.faces().iter().map(|f| f.name.clone()).collect();
    names.sort();
    names.dedup();
    assert!(names.contains(&"air|bronze".to_string()));
    assert_eq!(ev.faces().iter().filter(|f| f.name == "air|bronze").count(), 5);
    // 5 swept faces + the 6 faces of the air box; circle edges at the four corners off the axis plus those of ring edges
    assert_eq!(ev.edges3().iter().filter(|e| matches!(e.geom, geom::complex::EdgeGeom::Circle { .. })).count(), 4);
    // the slanted shoulder: its normal points up and out
    let sh = ev
        .faces()
        .iter()
        .position(|f| matches!(f.geom, geom::complex::FaceGeom::Quadric(geom::prim::Body::ConeSide { .. })))
        .unwrap();
    let nrm = ev.facet_normals(sh as u32);
    assert!(nrm[0][2] > 0.0, "{nrm:?}");
    // a revolve about x of a rect: same volume, world coordinates
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -4.0, -4.0), (8.0, 4.0, 4.0));
    let prof = m.rect("section", (0.0, 1.0), (2.0, 3.0)); // radius 0..2, along the axis 1..3
    let disk = m.revolve_along("disk", Axis::X, prof, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("metal", disk)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(1), PI * 4.0 * 2.0));
    assert_eq!(ev.classify([2.0, 1.0, 0.0]), Some(1));
    assert_eq!(ev.classify([0.5, 1.0, 0.0]), Some(0));
    assert_eq!(ev.classify([2.0, 3.0, 0.0]), Some(0));
    // refusals and round trip
    let mut m = Model::new_3d();
    let air = m.box_("air", (-4.0, -4.0, -2.0), (4.0, 4.0, 6.0));
    let bad = m.rect("crossing", (-1.0, 0.0), (1.0, 2.0));
    let r = m.revolve("spin", bad, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("r", r)]);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("crosses the axis"));
    let mut m = Model::new_3d();
    let air = m.box_("air", (-4.0, -4.0, -2.0), (4.0, 4.0, 6.0));
    let prof = m.rect("touching", (1.0, -2.0), (2.0, 1.0)); // reaches the floor of the air box
    let r = m.revolve("ring", prof, (0.0, 0.0, 0.0));
    m.overlay(&[("air", air), ("r", r)]);
    let e = m.eval(Tier::Complex).err().unwrap().0;
    assert!(e.contains("revolve 'ring'") && e.contains("air.zmin"), "{e}");
    assert_eq!(common::interp::build(&m.list.to_rust()).list, m.list);
}

#[test]
fn polygons_and_sheets_complement_subtraction_and_edge_selectors() {
    use geom::expr::Expr;
    use geom::select::EntityKind;

    // 1. Touching pad_1, jj_1, pad_2 sheets via Op::Polygons + Op::Sheet on substrate|vacuum (z = 0).
    // Include a subtle 20 nm deflection at the pad_1/jj_1 contact corner (0.0, -10e-6 + 20e-9)
    // to verify cross-feature RDP vertex pinning with simplify_tol = 50e-9, plus two holes in pad_1
    // (one 1 um^2 < min_hole_area = 4 um^2 dropped, one 4 um x 4 um = 16 um^2 kept).
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
    let vac = m.box_("vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));
    // A 3D pillar inside vacuum above part of pad_2 (x in [0, 50e-6], y in [20e-6, 60e-6], z in [0, 50e-6])
    // to force an internal overlay seam inside pad_2 and verify that `pad_2:edge` excludes internal seams.
    let pillar = m.box_("pillar", (0.0, 20e-6, 0.0), (50e-6, 60e-6, 50e-6));

    let y_contact1 = -10e-6;
    let y_deflect = y_contact1 + 20e-9; // 20 nm deflection < 50 nm simplify_tol
    let pad1_outer = vec![
        pt(-60e-6, -90e-6),
        pt(0.0, -90e-6 + 10e-9), // unshared 10 nm jitter -> simplified away by RDP
        pt(60e-6, -90e-6),
        pt(60e-6, y_contact1),
        pt(2e-6, y_deflect),  // shared contact corner with jj_1 -> pinned!
        pt(-2e-6, y_deflect), // shared contact corner with jj_1 -> pinned!
        pt(-60e-6, y_contact1),
    ];
    let hole_small = vec![pt(-30e-6, -60e-6), pt(-29e-6, -60e-6), pt(-29e-6, -59e-6), pt(-30e-6, -59e-6)]; // 1e-12 m^2
    let hole_flux = vec![pt(20e-6, -60e-6), pt(24e-6, -60e-6), pt(24e-6, -56e-6), pt(20e-6, -56e-6)]; // 16e-12 m^2
    let p1_poly = m.polygons_with("pad_1_poly", vec![pad1_outer], vec![hole_small, hole_flux], 50e-9, 4e-12);

    let jj_outer = vec![pt(-2e-6, y_deflect), pt(2e-6, y_deflect), pt(2e-6, 10e-6), pt(-2e-6, 10e-6)];
    let jj_poly = m.polygons_with("jj_1_poly", vec![jj_outer], vec![], 50e-9, 0.0);

    let pad2_outer = vec![pt(-60e-6, 10e-6), pt(60e-6, 10e-6), pt(60e-6, 90e-6), pt(-60e-6, 90e-6)];
    let p2_poly = m.polygons("pad_2_poly", vec![pad2_outer], vec![]);

    // Also add an airbridge sheet at z = 25e-6 wholly inside homogeneous vacuum (rb == ra)
    let bridge_rect = m.rect("bridge_rect", (-10e-6, -15e-6), (10e-6, 15e-6));

    m.sheet("pad_1", p1_poly, 0.0);
    m.sheet("jj_1", jj_poly, 0.0);
    m.sheet("pad_2", p2_poly, 0.0);
    m.sheet("bridge", bridge_rect, 25e-6);

    // Note: pillar is painted with the same region name or a dielectric cap; let's use two regions substrate & vacuum,
    // plus an auxiliary cut/feature or third region to split pad_2 into 2 overlay cells.
    let vac_main = m.difference("vac_main", vac, &[pillar]);
    let vac_union = m.union("vac_all", &[vac_main, pillar]);
    m.overlay(&[("substrate", sub), ("vacuum", vac_union)]);
    m.size("jj_edge_sz", Selector::parse("jj_1:edge").unwrap(), Some(Expr::Num(500e-9)), None, None);

    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();

    // Check Complement Interface Subtraction (SA) at z = 0:
    // area("substrate|vacuum") + area("pad_1") + area("jj_1") + area("pad_2") == die_area (400 um x 400 um)
    let die_area = 400e-6 * 400e-6;
    let face_area = |name: &str| -> f64 { ev.faces().iter().filter(|f| f.name == name).map(|f| f.area).sum() };
    let a_iface = face_area("substrate|vacuum");
    let a_pad1 = face_area("pad_1");
    let a_jj1 = face_area("jj_1");
    let a_pad2 = face_area("pad_2");
    let a_bridge = face_area("bridge");

    assert!(a_pad1 > 0.0 && a_jj1 > 0.0 && a_pad2 > 0.0 && a_bridge > 0.0);
    // pad_2 is 120 um x 80 um = 9.6e-9 m^2
    assert!(close(a_pad2, 120e-6 * 80e-6), "a_pad2 = {a_pad2:e}");
    // bridge is 20 um x 30 um = 6e-10 m^2 inside homogeneous vacuum (first == second == Some(1))
    assert!(close(a_bridge, 20e-6 * 30e-6), "a_bridge = {a_bridge:e}");
    for f in ev.faces().iter().filter(|f| f.name == "bridge") {
        assert_eq!(f.first, Some(1));
        assert_eq!(f.second, Some(1));
    }
    // pad_1 outer area has the 16 um^2 flux hole subtracted, but the 1 um^2 small hole filtered out
    // Exact outer area of pad_1: rectangle 120x80 um + triangle/trapezoid deflection on [-60..60] with [-2..2] at +20 nm
    let pad1_outer_exact = 120e-6 * 80e-6 + 0.5 * (120e-6 + 4e-6) * 20e-9;
    let pad1_expected = pad1_outer_exact - 16e-12;
    assert!(close(a_pad1, pad1_expected), "a_pad1 = {a_pad1:e} vs {pad1_expected:e}");
    // jj_1 area: 4 um x (20 um - 20 nm)
    let jj1_expected = 4e-6 * (20e-6 - 20e-9);
    assert!(close(a_jj1, jj1_expected), "a_jj1 = {a_jj1:e} vs {jj1_expected:e}");
    // Complement subtraction invariant to < 1e-12
    let total_z0 = a_iface + a_pad1 + a_jj1 + a_pad2;
    assert!(
        (total_z0 - die_area).abs() < 1e-12 * die_area,
        "complement subtraction mismatch: {total_z0:e} vs {die_area:e} (diff = {:e})",
        (total_z0 - die_area).abs()
    );

    // Check `<sheet>:edge` selectors (both `name(pad_2:edge)` and bare `pad_2:edge`):
    // Even though `pad_2` is split into multiple cap cells by `pillar` and `bridge`,
    // `pad_2:edge` must select ONLY the exterior 1D perimeter of `pad_2` (2*(120 + 80) = 400 um),
    // excluding internal fragment seams!
    let p2_faces_count = ev.faces().iter().filter(|f| f.name == "pad_2").count();
    assert!(p2_faces_count >= 2, "pad_2 should be split into multiple cap fragments by pillar/bridge");
    let p2_edge_sel1 = ev.resolve(&Selector::parse("name(pad_2:edge)").unwrap()).unwrap().set;
    let p2_edge_sel2 = ev.resolve(&Selector::parse("pad_2:edge").unwrap()).unwrap().set;
    assert_eq!(p2_edge_sel1.kind, EntityKind::Edges);
    assert_eq!(p2_edge_sel1, p2_edge_sel2);
    let p2_perim: f64 = p2_edge_sel1.ids.iter().map(|&e| ev.edges3()[e as usize].length).sum();
    assert!(close(p2_perim, 2.0 * (120e-6 + 80e-6)), "pad_2 exterior perimeter = {p2_perim:e}");

    // And `bridge:edge` inside homogeneous vacuum also resolves cleanly to its 100 um perimeter
    let br_edges = ev.resolve(&Selector::parse("bridge:edge").unwrap()).unwrap().set;
    assert_eq!(br_edges.kind, EntityKind::Edges);
    let br_perim: f64 = br_edges.ids.iter().map(|&e| ev.edges3()[e as usize].length).sum();
    assert!(close(br_perim, 2.0 * (20e-6 + 30e-6)), "bridge exterior perimeter = {br_perim:e}");
}

#[test]
fn abutting_sheets_placed_via_transform_and_pattern_pin_shared_vertices_conformally() {
    use geom::expr::Expr;
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];

    // 1. Two abutting sheets where sheet_b is placed via Op::Transform (Translate +50 um in x).
    // Both have a shared contact vertex at world (50e-6, 20e-9) (20 nm deflection < 50 nm simplify_tol)
    // and an unshared 10 nm jitter that should be simplified away.
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
    let vac = m.box_("vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));

    let ring_a = vec![
        pt(0.0, -20e-6),
        pt(25e-6, -20e-6 + 10e-9), // unshared 10 nm jitter -> simplified
        pt(50e-6, -20e-6),
        pt(50e-6, 20e-9), // shared with transformed sheet_b at world (50e-6, 20e-9) -> MUST be pinned!
        pt(50e-6, 20e-6),
        pt(0.0, 20e-6),
    ];
    let poly_a = m.polygons_with("poly_a", vec![ring_a], vec![], 50e-9, 0.0);

    // Local coordinates for sheet_b before +50 um x-translation: contact vertex is at local (0.0, 20e-9)
    let ring_b_local = vec![
        pt(0.0, -20e-6),
        pt(50e-6, -20e-6),
        pt(50e-6, 20e-6),
        pt(0.0, 20e-6),
        pt(0.0, 20e-9), // maps to world (50e-6, 20e-9) after translate -> MUST be pinned!
    ];
    let poly_b_local = m.polygons_with("poly_b_local", vec![ring_b_local], vec![], 50e-9, 0.0);
    let poly_b_xf = m.translate("poly_b_xf", poly_b_local, (50e-6, 0.0, 0.0));

    // 2. Also pattern a unit cell with Op::Pattern (Linear step 50e-6, count = 2) at y in [40e-6, 80e-6],
    // and abut a distinct sheet_c at x = 100e-6 against instance [1]'s right contact vertex (100e-6 + 20e-9, 60e-6)!
    let unit_ring = vec![
        pt(0.0, 40e-6),
        pt(25e-6, 40e-6 + 10e-9), // unshared 10 nm jitter -> simplified
        pt(50e-6, 40e-6),
        pt(50e-6 + 20e-9, 60e-6), // right contact protrusion -> instance [1] places this at (100e-6 + 20e-9, 60e-6)!
        pt(50e-6, 80e-6),
        pt(0.0, 80e-6),
        pt(20e-9, 60e-6), // left contact notch -> instance [1] places this at (50e-6 + 20e-9, 60e-6)
    ];
    let unit_poly = m.polygons_with("unit_poly", vec![unit_ring], vec![], 50e-9, 0.0);
    let pat_poly = m.linear_pattern("pat_poly", unit_poly, (50e-6, 0.0, 0.0), 2);

    let ring_c = vec![
        pt(100e-6, 40e-6),
        pt(150e-6, 40e-6),
        pt(150e-6, 80e-6),
        pt(100e-6, 80e-6),
        pt(100e-6 + 20e-9, 60e-6), // abuts instance [1] of pat_poly at (100e-6 + 20e-9, 60e-6)!
    ];
    let poly_c = m.polygons_with("poly_c", vec![ring_c], vec![], 50e-9, 0.0);

    m.sheet("sheet_a", poly_a, 0.0);
    m.sheet("sheet_b", poly_b_xf, 0.0);
    m.sheet("sheet_pat", pat_poly, 0.0);
    m.sheet("sheet_c", poly_c, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();

    // Verify the shared contact vertex (50e-6, 20e-9, 0.0) and pattern contact vertex (100e-6 + 20e-9, 60e-6, 0.0)
    // both exist in the 3D complex vertices, while the unshared 10 nm jitter (25e-6, -20e-6 + 10e-9, 0.0) was removed!
    let has_v = |wx: f64, wy: f64| -> bool {
        ev.vertices().iter().any(|v| (v.p[0] - wx).abs() < 1e-12 && (v.p[1] - wy).abs() < 1e-12 && v.p[2].abs() < 1e-12)
    };
    assert!(has_v(50e-6, 20e-9), "transformed contact vertex (50 um, 20 nm) must be pinned");
    assert!(has_v(100e-6 + 20e-9, 60e-6), "patterned contact vertex (100.02 um, 60 um) must be pinned");
    assert!(!has_v(25e-6, -20e-6 + 10e-9), "unshared 10 nm jitter must be simplified away");
}

#[test]
fn stepped_airbridge_stays_prismatic_and_evaluates_conformally_in_3d() {
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
    let vac = m.box_("vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));
    let cpw_rect = m.rect("cpw_rect", (-10e-6, -150e-6), (10e-6, 150e-6));
    m.sheet("cpw", cpw_rect, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    // Add stepped airbridge spanning across the CPW from x = -15 um to x = +15 um at z = 3 um
    let ab = m.stepped_airbridge("ab1", (-15e-6, 0.0), (15e-6, 0.0), 10e-6, 8e-6, 0.0, 3e-6);
    assert!(!geom::planar::needs_planar(&m.list), "stepped_airbridge must keep needs_planar == false");
    assert_eq!(ab.hole_region("ab1_pads").kind, RegionKind::Hole);

    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(ev.complex().has_sheets());
    assert_eq!(ev.complex().sheet_names(), vec!["ab1".to_string(), "cpw".to_string()]);

    // Vacuum volume has the two 8 um x 10 um x 3 um landing pad pillars removed
    let vac_vol = 400e-6 * 400e-6 * 100e-6 - 2.0 * (8e-6 * 10e-6 * 3e-6);
    assert!(close(ev.volume(1), vac_vol), "got {:e} vs {:e}", ev.volume(1), vac_vol);

    // Bridge deck sheet area is 30 um x 10 um = 300 um^2, and ab1:edge perimeter is 80 um
    let ab_area: f64 = ev.faces().iter().filter(|f| f.name == "ab1").map(|f| f.area).sum();
    assert!(close(ab_area, 30e-6 * 10e-6));
    let ab_edges = ev.resolve(&Selector::parse("ab1:edge").unwrap()).unwrap().set;
    let ab_perim: f64 = ab_edges.ids.iter().map(|&e| ev.edges3()[e as usize].length).sum();
    assert!(close(ab_perim, 2.0 * (30e-6 + 10e-6)));
}
