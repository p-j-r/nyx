// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Sizes on faces in the 3D lattice mesher (ADR 83): `add_size` on a face and the chord limit `max_dev` refine the mesh next to the face,
//! a mesh without them is unchanged, and the result does not depend on the number of threads.
#![allow(clippy::needless_range_loop)]
use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
use femcore::meshing::geomdrive3::mesh_geom_3d;
use geom::model::{Model, RegionKind};
use geom::select::Selector;

const R: f64 = 1.0;

fn sphere(face_h: Option<f64>) -> geom::model::FeatureList {
    let mut m = Model::new_3d();
    let host = m.box_("host", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let ball = m.sphere("cavity", (0.0, 0.0, 0.0), R);
    m.overlay_full(&[("host", host, RegionKind::Hole, None), ("cavity", ball, RegionKind::Material, None)], &[]);
    if let Some(h) = face_h {
        m.size("surface size", Selector::parse("name(cavity.surface)").unwrap(), Some(h.into()), None, None);
    }
    m.list
}

/// The largest depth of a boundary facet's centroid below the sphere: the chord error of the polyhedron.
fn chord_error(g: &femcore::meshing::geomdrive::MeshedGeom) -> f64 {
    let m = &g.mesh;
    let mut worst = 0.0f64;
    for f in 0..m.n_facets() {
        let mut c = [0.0; 3];
        for &v in m.facet(f) {
            for k in 0..3 {
                c[k] += m.points[3 * v as usize + k] / 3.0;
            }
        }
        worst = worst.max(R - (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt());
    }
    worst
}

#[test]
fn a_size_on_the_sphere_face_refines_it_and_lowers_the_chord_error() {
    let opts = GeomMeshOpts::new(0.5);
    let base = mesh_geom_3d(evaluate(&sphere(None)).unwrap(), &opts).unwrap();
    let fine = mesh_geom_3d(evaluate(&sphere(Some(0.25))).unwrap(), &opts).unwrap();
    let (eb, ef) = (chord_error(&base), chord_error(&fine));
    eprintln!("cells {} -> {}, chord error {eb:.3e} -> {ef:.3e}", base.mesh.n_cells(), fine.mesh.n_cells());
    assert!(fine.mesh.n_cells() > 2 * base.mesh.n_cells());
    // the error of a polyhedron falls as h^2: half the size, about a quarter of the error (the centroid depth of a graded, bisected facet
    // mesh falls somewhat slower than the edge sagitta)
    assert!(ef < eb / 2.0, "{eb} -> {ef}");
}

#[test]
fn max_dev_sets_the_size_of_a_curved_face_and_a_coarser_size_changes_nothing() {
    let mut opts = GeomMeshOpts::new(0.5);
    let base = mesh_geom_3d(evaluate(&sphere(None)).unwrap(), &opts).unwrap();
    opts.max_dev = Some(0.01);
    let g = mesh_geom_3d(evaluate(&sphere(None)).unwrap(), &opts).unwrap();
    // chord limit 0.01 on a radius 1 sphere: h_face = 2 sqrt(2 r 0.01) = 0.283
    assert!(g.mesh.n_cells() > 2 * base.mesh.n_cells());
    // the centroid depth of a facet is a little more than the sagitta of its edges that the limit bounds
    assert!(chord_error(&g) < 0.75 * chord_error(&base), "{} vs {}", chord_error(&g), chord_error(&base));
    // no request, or one that is no finer than the global size: the same mesh
    let same = mesh_geom_3d(evaluate(&sphere(Some(0.8))).unwrap(), &GeomMeshOpts::new(0.5)).unwrap();
    assert_eq!(same.mesh.points, base.mesh.points);
    assert_eq!(same.mesh.cells, base.mesh.cells);
    let loose = GeomMeshOpts { max_dev: Some(1.0), ..GeomMeshOpts::new(0.5) };
    let same = mesh_geom_3d(evaluate(&sphere(None)).unwrap(), &loose).unwrap();
    assert_eq!(same.mesh.cells, base.mesh.cells);
}

#[test]
fn a_face_size_gives_the_same_mesh_for_any_number_of_threads() {
    let make = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        pool.install(|| mesh_geom_3d(evaluate(&sphere(Some(0.25))).unwrap(), &GeomMeshOpts::new(0.5)).unwrap().mesh)
    };
    let (a, b) = (make(1), make(6));
    assert_eq!(a.points, b.points);
    assert_eq!(a.cells, b.cells);
}
