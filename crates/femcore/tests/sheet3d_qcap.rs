// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for `femcore::meshing::sheet3d` cQED enhancements (QFR-1, QFR-5):
//! - Unmodified curved-hole concentric transmon PSLG (`repro_concentric_mesh_fail_P2_fast`) meshes cleanly;
//! - Automatic coplanar gap proximity sizing (`gap_cells`) refines facing sheet rims symmetrically;
//! - Conforming pure-tetrahedral vertical coarsening (`z_coarsen`) produces a watertight, positive-Jacobian,
//!   crack-free mesh with `>= 3x` fewer 3D cells and P2 DOFs while changing P1/P2 capacitance by `< 0.2%`.

use femcore::fem::{assemble_laplace, dirichlet_constraints, solve_spd};
use femcore::mesh::Mesh;
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::request::{mesh_request, MeshRequest, MesherRequest};
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::props::LinearSolver;
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use std::collections::{BTreeMap, HashMap, HashSet};

fn pt(x: f64, y: f64) -> [Expr; 2] {
    [Expr::Num(x), Expr::Num(y)]
}

fn tet_det3(pts: &[f64], v: &[u32]) -> f64 {
    let p = |i: u32| [pts[3 * i as usize], pts[3 * i as usize + 1], pts[3 * i as usize + 2]];
    let (p0, p1, p2, p3) = (p(v[0]), p(v[1]), p(v[2]), p(v[3]));
    let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let c = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn tri_area3(pts: &[f64], t: &[u32]) -> f64 {
    let p = |v: u32| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]];
    let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

/// Verifies that `g.mesh` is a strictly positive-volume, manifold, crack-free pure-tetrahedral mesh
/// whose single-cell faces match the exact outer boundary facets and whose interior facets have 2 adjacent tets.
fn verify_watertight_pure_tet_mesh(g: &MeshedGeom, expected_vol: f64, interior_names: &[&str]) {
    let mesh = &g.mesh;
    assert_eq!(mesh.dim, 3);
    assert!(Mesh::from_solver_mesh(mesh).is_ok(), "SimplexMesh -> Mesh exchange must succeed");

    let mut total_vol = 0.0;
    let mut face_counts: HashMap<[u32; 3], u32> = HashMap::new();
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        assert!(
            v[0] != v[1] && v[0] != v[2] && v[0] != v[3] && v[1] != v[2] && v[1] != v[3] && v[2] != v[3],
            "degenerate repeated vertex in tet {c}: {v:?}"
        );
        let d3 = tet_det3(&mesh.points, v);
        assert!(d3 > 0.0, "tet {c} has non-positive Jacobian det3 = {d3:e}");
        total_vol += d3 / 6.0;
        for skip in 0..4 {
            let mut f = [0u32; 3];
            let mut k = 0;
            for (i, &vi) in v.iter().enumerate() {
                if i != skip {
                    f[k] = vi;
                    k += 1;
                }
            }
            f.sort_unstable();
            *face_counts.entry(f).or_insert(0) += 1;
        }
    }
    assert!(
        (total_vol - expected_vol).abs() < 1e-10 * expected_vol,
        "volume mismatch: got {total_vol:e}, expected {expected_vol:e}"
    );

    for (f, &cnt) in &face_counts {
        assert!(cnt == 1 || cnt == 2, "non-manifold tet face {f:?} shared by {cnt} tets");
    }

    let interior_labels: HashSet<i32> = interior_names.iter().flat_map(|name| g.boundary(name)).collect();

    let ext_flags = mesh.exterior_facets();
    let mut seen_facets: HashSet<[u32; 3]> = HashSet::new();
    let mut outer_facet_keys: HashSet<[u32; 3]> = HashSet::new();
    for (fi, &is_ext) in ext_flags.iter().enumerate().take(mesh.n_facets()) {
        let fv = mesh.facet(fi);
        let mut key = [fv[0], fv[1], fv[2]];
        key.sort_unstable();
        assert!(seen_facets.insert(key), "duplicate boundary/sheet facet {key:?} at index {fi}");
        let cnt = face_counts.get(&key).copied().unwrap_or(0);
        let lab = mesh.facet_labels[fi];
        if interior_labels.contains(&lab) {
            assert_eq!(cnt, 2, "interior sheet/interface facet {fi} (label {lab}) must belong to 2 tets, got {cnt}");
            assert!(!is_ext, "interior facet {fi} must not be marked exterior");
        } else {
            assert_eq!(cnt, 1, "outer boundary facet {fi} (label {lab}) must belong to 1 tet, got {cnt}");
            assert!(is_ext, "outer boundary facet {fi} must be marked exterior");
            outer_facet_keys.insert(key);
        }
    }

    let single_tet_faces: HashSet<[u32; 3]> =
        face_counts.into_iter().filter_map(|(k, cnt)| (cnt == 1).then_some(k)).collect();
    assert_eq!(
        single_tet_faces.len(),
        outer_facet_keys.len(),
        "interior crack detected: {} single-tet faces vs {} outer boundary facets",
        single_tet_faces.len(),
        outer_facet_keys.len()
    );
    assert_eq!(single_tet_faces, outer_facet_keys, "single-tet faces must match outer boundary facets exactly");
}

fn build_repro_concentric_model() -> Model {
    let mut m = Model::new_3d();
    let vac_top = m.box_("vac_top", (-0.00016208, -0.00018958, 0.0), (0.00016208, 0.00018958, 0.001));
    let vac_bot = m.box_("vac_bot", (-0.00016208, -0.00018958, -0.00135), (0.00016208, 0.00018958, -0.00035));
    let si_sub = m.box_("si_sub", (-0.00016208, -0.00018958, -0.00035), (0.00016208, 0.00018958, 0.0));

    let c0_outer = vec![
        pt(-0.00016208, -0.00018958),
        pt(-0.00016208, 0.00018958),
        pt(0.00016208, 0.00018958),
        pt(0.00016208, -0.00018958),
    ];
    let c0_hole: Vec<[Expr; 2]> = [
        [-1.1e-05, -0.00013008], [1.1e-05, -0.00013008], [1.1e-05, -9.7448e-05], [1.2802e-05, -9.7241e-05],
        [1.5977e-05, -9.677e-05], [1.9134e-05, -9.6195e-05], [2.2272e-05, -9.5518e-05], [2.5385e-05, -9.4738e-05],
        [2.8471e-05, -9.3857e-05], [3.1527e-05, -9.2875e-05], [3.4549e-05, -9.1794e-05], [3.7534e-05, -9.0614e-05],
        [4.0478e-05, -8.9338e-05], [4.338e-05, -8.7965e-05], [4.6235e-05, -8.6499e-05], [4.904e-05, -8.494e-05],
        [5.1793e-05, -8.329e-05], [5.449e-05, -8.1551e-05], [5.7129e-05, -7.9724e-05], [5.9707e-05, -7.7812e-05],
        [6.2221e-05, -7.5817e-05], [6.4669e-05, -7.374e-05], [6.7047e-05, -7.1585e-05], [6.9353e-05, -6.9353e-05],
        [7.1585e-05, -6.7047e-05], [7.374e-05, -6.4669e-05], [7.5817e-05, -6.2221e-05], [7.7812e-05, -5.9707e-05],
        [7.9724e-05, -5.7129e-05], [8.1551e-05, -5.449e-05], [8.329e-05, -5.1793e-05], [8.494e-05, -4.904e-05],
        [8.6499e-05, -4.6235e-05], [8.7965e-05, -4.338e-05], [8.9338e-05, -4.0478e-05], [9.0614e-05, -3.7534e-05],
        [9.1794e-05, -3.4549e-05], [9.2875e-05, -3.1527e-05], [9.3857e-05, -2.8471e-05], [9.4738e-05, -2.5385e-05],
        [9.5518e-05, -2.2272e-05], [9.6195e-05, -1.9134e-05], [9.677e-05, -1.5977e-05], [9.7241e-05, -1.2802e-05],
        [9.7608e-05, -9.614e-06], [9.787e-05, -6.415e-06], [9.8027e-05, -3.209e-06], [9.808e-05, 0.0],
        [9.8027e-05, 3.209e-06], [9.787e-05, 6.415e-06], [9.7608e-05, 9.614e-06], [9.7241e-05, 1.2802e-05],
        [9.677e-05, 1.5977e-05], [9.6195e-05, 1.9134e-05], [9.5518e-05, 2.2272e-05], [9.4738e-05, 2.5385e-05],
        [9.3857e-05, 2.8471e-05], [9.2875e-05, 3.1527e-05], [9.1794e-05, 3.4549e-05], [9.0614e-05, 3.7534e-05],
        [8.9338e-05, 4.0478e-05], [8.7965e-05, 4.338e-05], [8.6499e-05, 4.6235e-05], [8.494e-05, 4.904e-05],
        [8.329e-05, 5.1793e-05], [8.1551e-05, 5.449e-05], [7.9724e-05, 5.7129e-05], [7.7812e-05, 5.9707e-05],
        [7.5817e-05, 6.2221e-05], [7.374e-05, 6.4669e-05], [7.1585e-05, 6.7047e-05], [6.9353e-05, 6.9353e-05],
        [6.7047e-05, 7.1585e-05], [6.4669e-05, 7.374e-05], [6.2221e-05, 7.5817e-05], [5.9707e-05, 7.7812e-05],
        [5.7129e-05, 7.9724e-05], [5.449e-05, 8.1551e-05], [5.1793e-05, 8.329e-05], [4.904e-05, 8.494e-05],
        [4.6235e-05, 8.6499e-05], [4.338e-05, 8.7965e-05], [4.0478e-05, 8.9338e-05], [3.7534e-05, 9.0614e-05],
        [3.4549e-05, 9.1794e-05], [3.1527e-05, 9.2875e-05], [2.8471e-05, 9.3857e-05], [2.5385e-05, 9.4738e-05],
        [2.2272e-05, 9.5518e-05], [1.9134e-05, 9.6195e-05], [1.5977e-05, 9.677e-05], [1.2802e-05, 9.7241e-05],
        [1.1e-05, 9.7448e-05], [1.1e-05, 0.00012508], [-1.1e-05, 0.00012508], [-1.1e-05, 9.7448e-05],
        [-1.2802e-05, 9.7241e-05], [-1.5977e-05, 9.677e-05], [-1.9134e-05, 9.6195e-05], [-2.2272e-05, 9.5518e-05],
        [-2.5385e-05, 9.4738e-05], [-2.8471e-05, 9.3857e-05], [-3.1527e-05, 9.2875e-05], [-3.4549e-05, 9.1794e-05],
        [-3.7534e-05, 9.0614e-05], [-4.0478e-05, 8.9338e-05], [-4.338e-05, 8.7965e-05], [-4.6235e-05, 8.6499e-05],
        [-4.904e-05, 8.494e-05], [-5.1793e-05, 8.329e-05], [-5.449e-05, 8.1551e-05], [-5.7129e-05, 7.9724e-05],
        [-5.9707e-05, 7.7812e-05], [-6.2221e-05, 7.5817e-05], [-6.4669e-05, 7.374e-05], [-6.7047e-05, 7.1585e-05],
        [-6.9353e-05, 6.9353e-05], [-7.1585e-05, 6.7047e-05], [-7.374e-05, 6.4669e-05], [-7.5817e-05, 6.2221e-05],
        [-7.7812e-05, 5.9707e-05], [-7.9724e-05, 5.7129e-05], [-8.1551e-05, 5.449e-05], [-8.329e-05, 5.1793e-05],
        [-8.494e-05, 4.904e-05], [-8.6499e-05, 4.6235e-05], [-8.7965e-05, 4.338e-05], [-8.9338e-05, 4.0478e-05],
        [-9.0614e-05, 3.7534e-05], [-9.1794e-05, 3.4549e-05], [-9.2875e-05, 3.1527e-05], [-9.3857e-05, 2.8471e-05],
        [-9.4738e-05, 2.5385e-05], [-9.5518e-05, 2.2272e-05], [-9.6195e-05, 1.9134e-05], [-9.677e-05, 1.5977e-05],
        [-9.7241e-05, 1.2802e-05], [-9.7608e-05, 9.614e-06], [-9.787e-05, 6.415e-06], [-9.8027e-05, 3.209e-06],
        [-9.808e-05, 0.0], [-9.8027e-05, -3.209e-06], [-9.787e-05, -6.415e-06], [-9.7608e-05, -9.614e-06],
        [-9.7241e-05, -1.2802e-05], [-9.677e-05, -1.5977e-05], [-9.6195e-05, -1.9134e-05], [-9.5518e-05, -2.2272e-05],
        [-9.4738e-05, -2.5385e-05], [-9.3857e-05, -2.8471e-05], [-9.2875e-05, -3.1527e-05], [-9.1794e-05, -3.4549e-05],
        [-9.0614e-05, -3.7534e-05], [-8.9338e-05, -4.0478e-05], [-8.7965e-05, -4.338e-05], [-8.6499e-05, -4.6235e-05],
        [-8.494e-05, -4.904e-05], [-8.329e-05, -5.1793e-05], [-8.1551e-05, -5.449e-05], [-7.9724e-05, -5.7129e-05],
        [-7.7812e-05, -5.9707e-05], [-7.5817e-05, -6.2221e-05], [-7.374e-05, -6.4669e-05], [-7.1585e-05, -6.7047e-05],
        [-6.9353e-05, -6.9353e-05], [-6.7047e-05, -7.1585e-05], [-6.4669e-05, -7.374e-05], [-6.2221e-05, -7.5817e-05],
        [-5.9707e-05, -7.7812e-05], [-5.7129e-05, -7.9724e-05], [-5.449e-05, -8.1551e-05], [-5.1793e-05, -8.329e-05],
        [-4.904e-05, -8.494e-05], [-4.6235e-05, -8.6499e-05], [-4.338e-05, -8.7965e-05], [-4.0478e-05, -8.9338e-05],
        [-3.7534e-05, -9.0614e-05], [-3.4549e-05, -9.1794e-05], [-3.1527e-05, -9.2875e-05], [-2.8471e-05, -9.3857e-05],
        [-2.5385e-05, -9.4738e-05], [-2.2272e-05, -9.5518e-05], [-1.9134e-05, -9.6195e-05], [-1.5977e-05, -9.677e-05],
        [-1.2802e-05, -9.7241e-05], [-1.1e-05, -9.7448e-05],
    ].into_iter().map(|[x, y]| pt(x, y)).collect();
    let poly_c0 = m.polygons_with("poly_c0", vec![c0_outer], vec![c0_hole], 7.790615287643187e-07, 0.0);
    m.sheet("c0", poly_c0, 0.0);

    let c1_ring: Vec<[Expr; 2]> = [
        [-1.3e-05, -8.6628e-05], [-1.4271e-05, -8.644e-05], [-1.7092e-05, -8.5927e-05], [-1.9894e-05, -8.5321e-05],
        [-2.2675e-05, -8.4625e-05], [-2.5432e-05, -8.3838e-05], [-2.8161e-05, -8.2961e-05], [-3.0861e-05, -8.1995e-05],
        [-3.3527e-05, -8.0941e-05], [-3.6157e-05, -7.9801e-05], [-3.8749e-05, -7.8575e-05], [-4.1299e-05, -7.7265e-05],
        [-4.3805e-05, -7.5872e-05], [-4.6264e-05, -7.4399e-05], [-4.8674e-05, -7.2845e-05], [-5.1031e-05, -7.1213e-05],
        [-5.3334e-05, -6.9506e-05], [-5.5579e-05, -6.7723e-05], [-5.7765e-05, -6.5869e-05], [-5.989e-05, -6.3943e-05],
        [-6.195e-05, -6.195e-05], [-6.3943e-05, -5.989e-05], [-6.5869e-05, -5.7765e-05], [-6.7723e-05, -5.5579e-05],
        [-6.9506e-05, -5.3334e-05], [-7.1213e-05, -5.1031e-05], [-7.2845e-05, -4.8674e-05], [-7.4399e-05, -4.6264e-05],
        [-7.5872e-05, -4.3805e-05], [-7.7265e-05, -4.1299e-05], [-7.8575e-05, -3.8749e-05], [-7.9801e-05, -3.6157e-05],
        [-8.0941e-05, -3.3527e-05], [-8.1995e-05, -3.0861e-05], [-8.2961e-05, -2.8161e-05], [-8.3838e-05, -2.5432e-05],
        [-8.4625e-05, -2.2675e-05], [-8.5321e-05, -1.9894e-05], [-8.5927e-05, -1.7092e-05], [-8.644e-05, -1.4271e-05],
        [-8.686e-05, -1.1435e-05], [-8.7188e-05, -8.587e-06], [-8.7422e-05, -5.73e-06], [-8.7563e-05, -2.867e-06],
        [-8.761e-05, 0.0], [-8.7563e-05, 2.867e-06], [-8.7422e-05, 5.73e-06], [-8.7188e-05, 8.587e-06],
        [-8.686e-05, 1.1435e-05], [-8.644e-05, 1.4271e-05], [-8.5927e-05, 1.7092e-05], [-8.5321e-05, 1.9894e-05],
        [-8.4625e-05, 2.2675e-05], [-8.3838e-05, 2.5432e-05], [-8.2961e-05, 2.8161e-05], [-8.1995e-05, 3.0861e-05],
        [-8.0941e-05, 3.3527e-05], [-7.9801e-05, 3.6157e-05], [-7.8575e-05, 3.8749e-05], [-7.7265e-05, 4.1299e-05],
        [-7.5872e-05, 4.3805e-05], [-7.4399e-05, 4.6264e-05], [-7.2845e-05, 4.8674e-05], [-7.1213e-05, 5.1031e-05],
        [-6.9506e-05, 5.3334e-05], [-6.7723e-05, 5.5579e-05], [-6.5869e-05, 5.7765e-05], [-6.3943e-05, 5.989e-05],
        [-6.195e-05, 6.195e-05], [-5.989e-05, 6.3943e-05], [-5.7765e-05, 6.5869e-05], [-5.5579e-05, 6.7723e-05],
        [-5.3334e-05, 6.9506e-05], [-5.1031e-05, 7.1213e-05], [-4.8674e-05, 7.2845e-05], [-4.6264e-05, 7.4399e-05],
        [-4.3805e-05, 7.5872e-05], [-4.1299e-05, 7.7265e-05], [-3.8749e-05, 7.8575e-05], [-3.6157e-05, 7.9801e-05],
        [-3.3527e-05, 8.0941e-05], [-3.0861e-05, 8.1995e-05], [-2.8161e-05, 8.2961e-05], [-2.5432e-05, 8.3838e-05],
        [-2.2675e-05, 8.4625e-05], [-1.9894e-05, 8.5321e-05], [-1.7092e-05, 8.5927e-05], [-1.4271e-05, 8.644e-05],
        [-1.1435e-05, 8.686e-05], [-8.587e-06, 8.7188e-05], [-5.73e-06, 8.7422e-05], [-5e-06, 8.7458e-05],
        [-5e-06, 0.00012308], [5e-06, 0.00012308], [5e-06, 8.7458e-05], [5.73e-06, 8.7422e-05],
        [8.587e-06, 8.7188e-05], [1.1435e-05, 8.686e-05], [1.4271e-05, 8.644e-05], [1.7092e-05, 8.5927e-05],
        [1.9894e-05, 8.5321e-05], [2.2675e-05, 8.4625e-05], [2.5432e-05, 8.3838e-05], [2.8161e-05, 8.2961e-05],
        [3.0861e-05, 8.1995e-05], [3.3527e-05, 8.0941e-05], [3.6157e-05, 7.9801e-05], [3.8749e-05, 7.8575e-05],
        [4.1299e-05, 7.7265e-05], [4.3805e-05, 7.5872e-05], [4.6264e-05, 7.4399e-05], [4.8674e-05, 7.2845e-05],
        [5.1031e-05, 7.1213e-05], [5.3334e-05, 6.9506e-05], [5.5579e-05, 6.7723e-05], [5.7765e-05, 6.5869e-05],
        [5.989e-05, 6.3943e-05], [6.195e-05, 6.195e-05], [6.3943e-05, 5.989e-05], [6.5869e-05, 5.7765e-05],
        [6.7723e-05, 5.5579e-05], [6.9506e-05, 5.3334e-05], [7.1213e-05, 5.1031e-05], [7.2845e-05, 4.8674e-05],
        [7.4399e-05, 4.6264e-05], [7.5872e-05, 4.3805e-05], [7.7265e-05, 4.1299e-05], [7.8575e-05, 3.8749e-05],
        [7.9801e-05, 3.6157e-05], [8.0941e-05, 3.3527e-05], [8.1995e-05, 3.0861e-05], [8.2961e-05, 2.8161e-05],
        [8.3838e-05, 2.5432e-05], [8.4625e-05, 2.2675e-05], [8.5321e-05, 1.9894e-05], [8.5927e-05, 1.7092e-05],
        [8.644e-05, 1.4271e-05], [8.686e-05, 1.1435e-05], [8.7188e-05, 8.587e-06], [8.7422e-05, 5.73e-06],
        [8.7563e-05, 2.867e-06], [8.761e-05, 0.0], [8.7563e-05, -2.867e-06], [8.7422e-05, -5.73e-06],
        [8.7188e-05, -8.587e-06], [8.686e-05, -1.1435e-05], [8.644e-05, -1.4271e-05], [8.5927e-05, -1.7092e-05],
        [8.5321e-05, -1.9894e-05], [8.4625e-05, -2.2675e-05], [8.3838e-05, -2.5432e-05], [8.2961e-05, -2.8161e-05],
        [8.1995e-05, -3.0861e-05], [8.0941e-05, -3.3527e-05], [7.9801e-05, -3.6157e-05], [7.8575e-05, -3.8749e-05],
        [7.7265e-05, -4.1299e-05], [7.5872e-05, -4.3805e-05], [7.4399e-05, -4.6264e-05], [7.2845e-05, -4.8674e-05],
        [7.1213e-05, -5.1031e-05], [6.9506e-05, -5.3334e-05], [6.7723e-05, -5.5579e-05], [6.5869e-05, -5.7765e-05],
        [6.3943e-05, -5.989e-05], [6.195e-05, -6.195e-05], [5.989e-05, -6.3943e-05], [5.7765e-05, -6.5869e-05],
        [5.5579e-05, -6.7723e-05], [5.3334e-05, -6.9506e-05], [5.1031e-05, -7.1213e-05], [4.8674e-05, -7.2845e-05],
        [4.6264e-05, -7.4399e-05], [4.3805e-05, -7.5872e-05], [4.1299e-05, -7.7265e-05], [3.8749e-05, -7.8575e-05],
        [3.6157e-05, -7.9801e-05], [3.3527e-05, -8.0941e-05], [3.0861e-05, -8.1995e-05], [2.8161e-05, -8.2961e-05],
        [2.5432e-05, -8.3838e-05], [2.2675e-05, -8.4625e-05], [1.9894e-05, -8.5321e-05], [1.7092e-05, -8.5927e-05],
        [1.4271e-05, -8.644e-05], [1.3e-05, -8.6628e-05], [1.3e-05, -7.4415e-05], [1.4739e-05, -7.4098e-05],
        [1.7156e-05, -7.3576e-05], [1.9554e-05, -7.2976e-05], [2.1931e-05, -7.2297e-05], [2.4285e-05, -7.1541e-05],
        [2.6612e-05, -7.0708e-05], [2.8912e-05, -6.9799e-05], [3.118e-05, -6.8816e-05], [3.3415e-05, -6.7759e-05],
        [3.5614e-05, -6.6629e-05], [3.7775e-05, -6.5428e-05], [3.9896e-05, -6.4157e-05], [4.1973e-05, -6.2818e-05],
        [4.4006e-05, -6.1411e-05], [4.5992e-05, -5.9938e-05], [4.7928e-05, -5.8401e-05], [4.9814e-05, -5.6801e-05],
        [5.1645e-05, -5.5141e-05], [5.3422e-05, -5.3422e-05], [5.5141e-05, -5.1645e-05], [5.6801e-05, -4.9814e-05],
        [5.8401e-05, -4.7928e-05], [5.9938e-05, -4.5992e-05], [6.1411e-05, -4.4006e-05], [6.2818e-05, -4.1973e-05],
        [6.4157e-05, -3.9896e-05], [6.5428e-05, -3.7775e-05], [6.6629e-05, -3.5614e-05], [6.7759e-05, -3.3415e-05],
        [6.8816e-05, -3.118e-05], [6.9799e-05, -2.8912e-05], [7.0708e-05, -2.6612e-05], [7.1541e-05, -2.4285e-05],
        [7.2297e-05, -2.1931e-05], [7.2976e-05, -1.9554e-05], [7.3576e-05, -1.7156e-05], [7.4098e-05, -1.4739e-05],
        [7.4541e-05, -1.2307e-05], [7.4904e-05, -9.861e-06], [7.5186e-05, -7.405e-06], [7.5388e-05, -4.941e-06],
        [7.551e-05, -2.472e-06], [7.555e-05, 0.0], [7.551e-05, 2.472e-06], [7.5388e-05, 4.941e-06],
        [7.5323e-05, 5.73e-06], [7.5186e-05, 7.405e-06], [7.4904e-05, 9.861e-06], [7.4541e-05, 1.2307e-05],
        [7.4098e-05, 1.4739e-05], [7.3576e-05, 1.7156e-05], [7.2976e-05, 1.9554e-05], [7.2297e-05, 2.1931e-05],
        [7.1541e-05, 2.4285e-05], [7.0708e-05, 2.6612e-05], [6.9799e-05, 2.8912e-05], [6.8816e-05, 3.118e-05],
        [6.7759e-05, 3.3415e-05], [6.6629e-05, 3.5614e-05], [6.5428e-05, 3.7775e-05], [6.4157e-05, 3.9896e-05],
        [6.2818e-05, 4.1973e-05], [6.1411e-05, 4.4006e-05], [5.9938e-05, 4.5992e-05], [5.8401e-05, 4.7928e-05],
        [5.6801e-05, 4.9814e-05], [5.5141e-05, 5.1645e-05], [5.3422e-05, 5.3422e-05], [5.1645e-05, 5.5141e-05],
        [4.9814e-05, 5.6801e-05], [4.7928e-05, 5.8401e-05], [4.5992e-05, 5.9938e-05], [4.4006e-05, 6.1411e-05],
        [4.1973e-05, 6.2818e-05], [3.9896e-05, 6.4157e-05], [3.7775e-05, 6.5428e-05], [3.5614e-05, 6.6629e-05],
        [3.3415e-05, 6.7759e-05], [3.118e-05, 6.8816e-05], [2.8912e-05, 6.9799e-05], [2.6612e-05, 7.0708e-05],
        [2.4285e-05, 7.1541e-05], [2.1931e-05, 7.2297e-05], [1.9554e-05, 7.2976e-05], [1.7156e-05, 7.3576e-05],
        [1.4739e-05, 7.4098e-05], [1.2307e-05, 7.4541e-05], [9.861e-06, 7.4904e-05], [7.405e-06, 7.5186e-05],
        [5e-06, 7.5383e-05], [4.941e-06, 7.5388e-05], [2.472e-06, 7.551e-05], [0.0, 7.555e-05],
        [-2.472e-06, 7.551e-05], [-4.941e-06, 7.5388e-05], [-5e-06, 7.5383e-05], [-7.405e-06, 7.5186e-05],
        [-9.861e-06, 7.4904e-05], [-1.2307e-05, 7.4541e-05], [-1.4739e-05, 7.4098e-05], [-1.7156e-05, 7.3576e-05],
        [-1.9554e-05, 7.2976e-05], [-2.1931e-05, 7.2297e-05], [-2.4285e-05, 7.1541e-05], [-2.6612e-05, 7.0708e-05],
        [-2.8912e-05, 6.9799e-05], [-3.118e-05, 6.8816e-05], [-3.3415e-05, 6.7759e-05], [-3.5614e-05, 6.6629e-05],
        [-3.7775e-05, 6.5428e-05], [-3.9896e-05, 6.4157e-05], [-4.1973e-05, 6.2818e-05], [-4.4006e-05, 6.1411e-05],
        [-4.5992e-05, 5.9938e-05], [-4.7928e-05, 5.8401e-05], [-4.9814e-05, 5.6801e-05], [-5.1645e-05, 5.5141e-05],
        [-5.3422e-05, 5.3422e-05], [-5.5141e-05, 5.1645e-05], [-5.6801e-05, 4.9814e-05], [-5.8401e-05, 4.7928e-05],
        [-5.9938e-05, 4.5992e-05], [-6.1411e-05, 4.4006e-05], [-6.2818e-05, 4.1973e-05], [-6.4157e-05, 3.9896e-05],
        [-6.5428e-05, 3.7775e-05], [-6.6629e-05, 3.5614e-05], [-6.7759e-05, 3.3415e-05], [-6.8816e-05, 3.118e-05],
        [-6.9799e-05, 2.8912e-05], [-7.0708e-05, 2.6612e-05], [-7.1541e-05, 2.4285e-05], [-7.2297e-05, 2.1931e-05],
        [-7.2976e-05, 1.9554e-05], [-7.3576e-05, 1.7156e-05], [-7.4098e-05, 1.4739e-05], [-7.4541e-05, 1.2307e-05],
        [-7.4904e-05, 9.861e-06], [-7.5186e-05, 7.405e-06], [-7.5323e-05, 5.73e-06], [-7.5388e-05, 4.941e-06],
        [-7.551e-05, 2.472e-06], [-7.555e-05, 0.0], [-7.551e-05, -2.472e-06], [-7.5388e-05, -4.941e-06],
        [-7.5186e-05, -7.405e-06], [-7.4904e-05, -9.861e-06], [-7.4541e-05, -1.2307e-05], [-7.4098e-05, -1.4739e-05],
        [-7.3576e-05, -1.7156e-05], [-7.2976e-05, -1.9554e-05], [-7.2297e-05, -2.1931e-05], [-7.1541e-05, -2.4285e-05],
        [-7.0708e-05, -2.6612e-05], [-6.9799e-05, -2.8912e-05], [-6.8816e-05, -3.118e-05], [-6.7759e-05, -3.3415e-05],
        [-6.6629e-05, -3.5614e-05], [-6.5428e-05, -3.7775e-05], [-6.4157e-05, -3.9896e-05], [-6.2818e-05, -4.1973e-05],
        [-6.1411e-05, -4.4006e-05], [-5.9938e-05, -4.5992e-05], [-5.8401e-05, -4.7928e-05], [-5.6801e-05, -4.9814e-05],
        [-5.5141e-05, -5.1645e-05], [-5.3422e-05, -5.3422e-05], [-5.1645e-05, -5.5141e-05], [-4.9814e-05, -5.6801e-05],
        [-4.7928e-05, -5.8401e-05], [-4.5992e-05, -5.9938e-05], [-4.4006e-05, -6.1411e-05], [-4.1973e-05, -6.2818e-05],
        [-3.9896e-05, -6.4157e-05], [-3.7775e-05, -6.5428e-05], [-3.5614e-05, -6.6629e-05], [-3.3415e-05, -6.7759e-05],
        [-3.118e-05, -6.8816e-05], [-2.8912e-05, -6.9799e-05], [-2.6612e-05, -7.0708e-05], [-2.4285e-05, -7.1541e-05],
        [-2.1931e-05, -7.2297e-05], [-1.9554e-05, -7.2976e-05], [-1.7156e-05, -7.3576e-05], [-1.4739e-05, -7.4098e-05],
        [-1.3e-05, -7.4415e-05],
    ].into_iter().map(|[x, y]| pt(x, y)).collect();
    let poly_c1 = m.polygons_with("poly_c1", vec![c1_ring], vec![], 7.790615287643187e-07, 0.0);
    m.sheet("c1", poly_c1, 0.0);

    let c2_ring: Vec<[Expr; 2]> = [
        [-5e-06, -0.00012808], [-5e-06, -5.6463e-05], [-5.557e-06, -5.6417e-05], [-7.4e-06, -5.6205e-05],
        [-9.235e-06, -5.5933e-05], [-1.106e-05, -5.5601e-05], [-1.2873e-05, -5.5209e-05], [-1.4672e-05, -5.4758e-05],
        [-1.6456e-05, -5.4249e-05], [-1.8222e-05, -5.3681e-05], [-1.9969e-05, -5.3057e-05], [-2.1694e-05, -5.2375e-05],
        [-2.3396e-05, -5.1637e-05], [-2.5073e-05, -5.0844e-05], [-2.6723e-05, -4.9996e-05], [-2.8345e-05, -4.9095e-05],
        [-2.9936e-05, -4.8141e-05], [-3.1495e-05, -4.7136e-05], [-3.3021e-05, -4.608e-05], [-3.4511e-05, -4.4975e-05],
        [-3.5964e-05, -4.3822e-05], [-3.7378e-05, -4.2622e-05], [-3.8753e-05, -4.1376e-05], [-4.0086e-05, -4.0086e-05],
        [-4.1376e-05, -3.8753e-05], [-4.2622e-05, -3.7378e-05], [-4.3822e-05, -3.5964e-05], [-4.4975e-05, -3.4511e-05],
        [-4.608e-05, -3.3021e-05], [-4.7136e-05, -3.1495e-05], [-4.8141e-05, -2.9936e-05], [-4.9095e-05, -2.8345e-05],
        [-4.9996e-05, -2.6723e-05], [-5.0844e-05, -2.5073e-05], [-5.1637e-05, -2.3396e-05], [-5.2375e-05, -2.1694e-05],
        [-5.3057e-05, -1.9969e-05], [-5.3681e-05, -1.8222e-05], [-5.4249e-05, -1.6456e-05], [-5.4758e-05, -1.4672e-05],
        [-5.5209e-05, -1.2873e-05], [-5.5601e-05, -1.106e-05], [-5.5933e-05, -9.235e-06], [-5.6205e-05, -7.4e-06],
        [-5.6417e-05, -5.557e-06], [-5.6569e-05, -3.708e-06], [-5.666e-05, -1.855e-06], [-5.669e-05, 0.0],
        [-5.666e-05, 1.855e-06], [-5.6569e-05, 3.708e-06], [-5.6417e-05, 5.557e-06], [-5.6205e-05, 7.4e-06],
        [-5.5933e-05, 9.235e-06], [-5.5601e-05, 1.106e-05], [-5.5209e-05, 1.2873e-05], [-5.4758e-05, 1.4672e-05],
        [-5.4249e-05, 1.6456e-05], [-5.3681e-05, 1.8222e-05], [-5.3057e-05, 1.9969e-05], [-5.2375e-05, 2.1694e-05],
        [-5.1637e-05, 2.3396e-05], [-5.0844e-05, 2.5073e-05], [-4.9996e-05, 2.6723e-05], [-4.9095e-05, 2.8345e-05],
        [-4.8141e-05, 2.9936e-05], [-4.7136e-05, 3.1495e-05], [-4.608e-05, 3.3021e-05], [-4.4975e-05, 3.4511e-05],
        [-4.3822e-05, 3.5964e-05], [-4.2622e-05, 3.7378e-05], [-4.1376e-05, 3.8753e-05], [-4.0086e-05, 4.0086e-05],
        [-3.8753e-05, 4.1376e-05], [-3.7378e-05, 4.2622e-05], [-3.5964e-05, 4.3822e-05], [-3.4511e-05, 4.4975e-05],
        [-3.3021e-05, 4.608e-05], [-3.1495e-05, 4.7136e-05], [-2.9936e-05, 4.8141e-05], [-2.8345e-05, 4.9095e-05],
        [-2.6723e-05, 4.9996e-05], [-2.5073e-05, 5.0844e-05], [-2.3396e-05, 5.1637e-05], [-2.1694e-05, 5.2375e-05],
        [-1.9969e-05, 5.3057e-05], [-1.8222e-05, 5.3681e-05], [-1.6456e-05, 5.4249e-05], [-1.4672e-05, 5.4758e-05],
        [-1.2873e-05, 5.5209e-05], [-1.106e-05, 5.5601e-05], [-9.235e-06, 5.5933e-05], [-7.4e-06, 5.6205e-05],
        [-5.557e-06, 5.6417e-05], [-3.708e-06, 5.6569e-05], [-1.855e-06, 5.666e-05], [0.0, 5.669e-05],
        [1.855e-06, 5.666e-05], [3.708e-06, 5.6569e-05], [5.557e-06, 5.6417e-05], [7.4e-06, 5.6205e-05],
        [9.235e-06, 5.5933e-05], [1.106e-05, 5.5601e-05], [1.2873e-05, 5.5209e-05], [1.4672e-05, 5.4758e-05],
        [1.6456e-05, 5.4249e-05], [1.8222e-05, 5.3681e-05], [1.9969e-05, 5.3057e-05], [2.1694e-05, 5.2375e-05],
        [2.3396e-05, 5.1637e-05], [2.5073e-05, 5.0844e-05], [2.6723e-05, 4.9996e-05], [2.8345e-05, 4.9095e-05],
        [2.9936e-05, 4.8141e-05], [3.1495e-05, 4.7136e-05], [3.3021e-05, 4.608e-05], [3.4511e-05, 4.4975e-05],
        [3.5964e-05, 4.3822e-05], [3.7378e-05, 4.2622e-05], [3.8753e-05, 4.1376e-05], [4.0086e-05, 4.0086e-05],
        [4.1376e-05, 3.8753e-05], [4.2622e-05, 3.7378e-05], [4.3822e-05, 3.5964e-05], [4.4975e-05, 3.4511e-05],
        [4.608e-05, 3.3021e-05], [4.7136e-05, 3.1495e-05], [4.8141e-05, 2.9936e-05], [4.9095e-05, 2.8345e-05],
        [4.9996e-05, 2.6723e-05], [5.0844e-05, 2.5073e-05], [5.1637e-05, 2.3396e-05], [5.2375e-05, 2.1694e-05],
        [5.3057e-05, 1.9969e-05], [5.3681e-05, 1.8222e-05], [5.4249e-05, 1.6456e-05], [5.4758e-05, 1.4672e-05],
        [5.5209e-05, 1.2873e-05], [5.5601e-05, 1.106e-05], [5.5933e-05, 9.235e-06], [5.6205e-05, 7.4e-06],
        [5.6417e-05, 5.557e-06], [5.6569e-05, 3.708e-06], [5.666e-05, 1.855e-06], [5.669e-05, 0.0],
        [5.666e-05, -1.855e-06], [5.6569e-05, -3.708e-06], [5.6417e-05, -5.557e-06], [5.6205e-05, -7.4e-06],
        [5.5933e-05, -9.235e-06], [5.5601e-05, -1.106e-05], [5.5209e-05, -1.2873e-05], [5.4758e-05, -1.4672e-05],
        [5.4249e-05, -1.6456e-05], [5.3681e-05, -1.8222e-05], [5.3057e-05, -1.9969e-05], [5.2375e-05, -2.1694e-05],
        [5.1637e-05, -2.3396e-05], [5.0844e-05, -2.5073e-05], [4.9996e-05, -2.6723e-05], [4.9095e-05, -2.8345e-05],
        [4.8141e-05, -2.9936e-05], [4.7136e-05, -3.1495e-05], [4.608e-05, -3.3021e-05], [4.4975e-05, -3.4511e-05],
        [4.3822e-05, -3.5964e-05], [4.2622e-05, -3.7378e-05], [4.1376e-05, -3.8753e-05], [4.0086e-05, -4.0086e-05],
        [3.8753e-05, -4.1376e-05], [3.7378e-05, -4.2622e-05], [3.5964e-05, -4.3822e-05], [3.4511e-05, -4.4975e-05],
        [3.3021e-05, -4.608e-05], [3.1495e-05, -4.7136e-05], [2.9936e-05, -4.8141e-05], [2.8345e-05, -4.9095e-05],
        [2.6723e-05, -4.9996e-05], [2.5073e-05, -5.0844e-05], [2.3396e-05, -5.1637e-05], [2.1694e-05, -5.2375e-05],
        [1.9969e-05, -5.3057e-05], [1.8222e-05, -5.3681e-05], [1.6456e-05, -5.4249e-05], [1.4672e-05, -5.4758e-05],
        [1.2873e-05, -5.5209e-05], [1.106e-05, -5.5601e-05], [9.235e-06, -5.5933e-05], [7.4e-06, -5.6205e-05],
        [5.557e-06, -5.6417e-05], [5e-06, -5.6463e-05], [5e-06, -0.00012808],
    ].into_iter().map(|[x, y]| pt(x, y)).collect();
    let poly_c2 = m.polygons_with("poly_c2", vec![c2_ring], vec![], 7.790615287643187e-07, 0.0);
    m.sheet("c2", poly_c2, 0.0);

    m.overlay(&[("VacTop", vac_top), ("VacBot", vac_bot), ("Si", si_sub)]);
    m.size("sz_c1", Selector::parse("name(c1)").unwrap(), Some(Expr::Num(9.230483429295996e-06)), None, None);
    m.size("sz_c1_e", Selector::parse("c1:edge").unwrap(), Some(Expr::Num(5.429696134879997e-06)), None, None);
    m.size("sz_c2", Selector::parse("name(c2)").unwrap(), Some(Expr::Num(3.431672259133796e-05)), None, None);
    m.size("sz_c2_e", Selector::parse("c2:edge").unwrap(), Some(Expr::Num(2.018630740666939e-05)), None, None);
    m
}

#[test]
fn repro_concentric_mesh_fail_p2_fast_meshes_cleanly() {
    let m = build_repro_concentric_model();
    let g = mesh_model(&m.list, &MeshPlan::simple(0.0001725))
        .expect("unmodified concentric sheet model must mesh cleanly without manual c0:edge workaround");

    let die_area = (2.0 * 0.00016208) * (2.0 * 0.00018958);
    let expected_vol = die_area * (0.001 - (-0.00135));
    verify_watertight_pure_tet_mesh(&g, expected_vol, &["c0", "c1", "c2", "Si|VacTop", "Si|VacBot"]);

    let bnd_area = |name: &str| -> f64 {
        let labs = g.boundary(name);
        assert!(!labs.is_empty(), "missing boundary label for {name}");
        labs.into_iter().flat_map(|l| g.mesh.sheet_facets(l)).map(|t| tri_area3(&g.mesh.points, t)).sum()
    };
    let total_z0 = bnd_area("Si|VacTop") + bnd_area("c0") + bnd_area("c1") + bnd_area("c2");
    assert!(
        (total_z0 - die_area).abs() < 1e-11 * die_area,
        "z=0 complement area mismatch: {total_z0:e} vs {die_area:e}"
    );

    // Also verify `z_coarsen = true` on the same curved-hole concentric model is 100% watertight & reduces DOFs >= 3x
    let req_coarse = MeshRequest {
        h: 0.0001725,
        max_dev: None,
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: false,
        gap_cells: None,
        z_coarsen: true,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: false,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };
    let (g_zc, _) = mesh_request(&m.list, &req_coarse).expect("concentric model with z_coarsen=true must mesh cleanly");
    verify_watertight_pure_tet_mesh(&g_zc, expected_vol, &["c0", "c1", "c2", "Si|VacTop", "Si|VacBot"]);
    eprintln!(
        "[concentric] full cells={}, P2 DOFs={} -> z_coarsen cells={}, P2 DOFs={} ({:.2}x cells, {:.2}x P2 DOFs)",
        g.mesh.n_cells(),
        g.mesh.n_dofs(2),
        g_zc.mesh.n_cells(),
        g_zc.mesh.n_dofs(2),
        g.mesh.n_cells() as f64 / g_zc.mesh.n_cells() as f64,
        g.mesh.n_dofs(2) as f64 / g_zc.mesh.n_dofs(2) as f64,
    );
    assert!(
        g.mesh.n_cells() >= 3 * g_zc.mesh.n_cells(),
        "z_coarsen must reduce concentric cells >= 3x (full={}, coarse={})",
        g.mesh.n_cells(),
        g_zc.mesh.n_cells()
    );
}

#[test]
fn gap_cells_refines_coplanar_gaps_symmetrically_on_both_sides() {
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-250e-6, -250e-6, -150e-6), (250e-6, 250e-6, 0.0));
    let vac = m.box_("vac", (-250e-6, -250e-6, 0.0), (250e-6, 250e-6, 150e-6));

    // Ground plane with [-70..70] x [-70..70] um pocket cutout
    let gnd_outer = vec![pt(-250e-6, -250e-6), pt(250e-6, -250e-6), pt(250e-6, 250e-6), pt(-250e-6, 250e-6)];
    let gnd_pocket = vec![pt(-70e-6, -70e-6), pt(70e-6, -70e-6), pt(70e-6, 70e-6), pt(-70e-6, 70e-6)];
    let poly_gnd = m.polygons("poly_gnd", vec![gnd_outer], vec![gnd_pocket]);

    // Signal island [-60..60] x [-60..60] um inside pocket -> exact d_gap = 10 um on all 4 sides
    let sig_outer = vec![pt(-60e-6, -60e-6), pt(60e-6, -60e-6), pt(60e-6, 60e-6), pt(-60e-6, 60e-6)];
    let poly_sig = m.polygons("poly_sig", vec![sig_outer], vec![]);

    m.sheet("gnd", poly_gnd, 0.0);
    m.sheet("sig", poly_sig, 0.0);
    m.overlay(&[("Si", sub), ("Vac", vac)]);

    let req_base = MeshRequest {
        h: 80e-6,
        max_dev: None,
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: false,
        gap_cells: None,
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: false,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };
    let req_gap = MeshRequest { gap_cells: Some(2.5), ..req_base.clone() };

    let (g_base, _) = mesh_request(&m.list, &req_base).unwrap();
    let (g_gap, _) = mesh_request(&m.list, &req_gap).unwrap();

    // Invalid non-positive gap_cells must return a clean error
    assert!(mesh_request(&m.list, &MeshRequest { gap_cells: Some(0.0), ..req_base.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { gap_cells: Some(-2.0), ..req_base }).is_err());

    // Extract 1D perimeter edges of a sheet inside `|x|, |y| < 100 um` (excluding the outer die rim of `gnd`)
    let inner_rim_stats = |g: &MeshedGeom, sheet_name: &str| -> (usize, f64) {
        let mut set: HashSet<(u32, u32)> = HashSet::new();
        for lab in g.boundary(sheet_name) {
            for tri in g.mesh.sheet_facets(lab) {
                for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                    let e = (a.min(b), a.max(b));
                    if !set.insert(e) {
                        set.remove(&e);
                    }
                }
            }
        }
        let mut n_segs = 0usize;
        let mut max_len = 0.0f64;
        for (a, b) in set {
            let pa = g.mesh.point(a as usize);
            let pb = g.mesh.point(b as usize);
            if pa[0].abs() < 100e-6 && pa[1].abs() < 100e-6 && pb[0].abs() < 100e-6 && pb[1].abs() < 100e-6 {
                let len = ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2)).sqrt();
                n_segs += 1;
                max_len = max_len.max(len);
            }
        }
        (n_segs, max_len)
    };

    let (n_sig_base, max_sig_base) = inner_rim_stats(&g_base, "sig");
    let (n_gnd_base, max_gnd_base) = inner_rim_stats(&g_base, "gnd");
    let (n_sig_gap, max_sig_gap) = inner_rim_stats(&g_gap, "sig");
    let (n_gnd_gap, max_gnd_gap) = inner_rim_stats(&g_gap, "gnd");

    // Target edge length with d_gap = 10 um and gap_cells = 2.5 is 4.0 um (`<= 5.0 um` after 1.25x chord rule)
    assert!(max_sig_base > 20e-6 && max_gnd_base > 20e-6);
    assert!(max_sig_gap <= 5.05e-6, "signal rim max segment {max_sig_gap:e} exceeds d_gap / gap_cells target");
    assert!(max_gnd_gap <= 5.05e-6, "ground pocket rim max segment {max_gnd_gap:e} exceeds d_gap / gap_cells target");
    assert!(n_sig_gap >= 4 * n_sig_base, "signal rim not refined: base={n_sig_base}, gap={n_sig_gap}");
    assert!(n_gnd_gap >= 4 * n_gnd_base, "ground pocket rim not refined: base={n_gnd_base}, gap={n_gnd_gap}");

    // Both facing sides have comparable segment density (~120 um / 4 um vs ~140 um / 4 um per side)
    let avg_sig = (4.0 * 120e-6) / n_sig_gap as f64;
    let avg_gnd = (4.0 * 140e-6) / n_gnd_gap as f64;
    let sym_ratio = avg_sig.max(avg_gnd) / avg_sig.min(avg_gnd);
    assert!(sym_ratio < 1.25, "asymmetric gap refinement: avg_sig={avg_sig:e}, avg_gnd={avg_gnd:e}");
}

fn solve_terminal_capacitance_f(g: &MeshedGeom, order: usize, eps_by_region: &[(&str, f64)], active_sheet: &str, ground_names: &[&str]) -> f64 {
    let eps0 = 8.854_187_812_8e-12;
    let mut dom_eps = HashMap::new();
    for &(rname, eps_r) in eps_by_region {
        let did = g.domain(rname).unwrap_or_else(|| panic!("missing domain {rname}"));
        dom_eps.insert(did, eps_r);
    }
    let coeff: Vec<f64> = g.mesh.cell_labels.iter().map(|d| *dom_eps.get(d).unwrap_or(&1.0)).collect();
    let k_mat = assemble_laplace(&g.mesh, order, &coeff);

    let mut bc_vals: BTreeMap<i32, f64> = BTreeMap::new();
    for &gname in ground_names {
        for lab in g.boundary(gname) {
            bc_vals.insert(lab, 0.0);
        }
    }
    for lab in g.boundary(active_sheet) {
        bc_vals.insert(lab, 1.0);
    }
    let (fixed, values) = dirichlet_constraints(&g.mesh, &bc_vals, order);
    let rhs = vec![0.0; g.mesh.n_dofs(order)];
    let (u, _) = solve_spd(&g.mesh, order, 1, LinearSolver::Iterative, &k_mat, &rhs, &fixed, &values)
        .expect("3D Laplace Dirichlet solve must succeed");
    let ku = k_mat.matvec(&u);
    let energy: f64 = u.iter().zip(&ku).map(|(ui, kui)| ui * kui).sum();
    eps0 * energy
}

#[test]
fn z_coarsen_pure_tet_watertight_3x_dof_reduction_and_sub_0_2_pct_capacitance() {
    let mut m = Model::new_3d();
    let vac_bot = m.box_("vac_bot", (-250e-6, -250e-6, -1500e-6), (250e-6, 250e-6, -350e-6));
    let si_sub = m.box_("si_sub", (-250e-6, -250e-6, -350e-6), (250e-6, 250e-6, 0.0));
    let vac_top = m.box_("vac_top", (-250e-6, -250e-6, 0.0), (250e-6, 250e-6, 1500e-6));

    // Compact pocket transmon on z = 0: ground plane pocket and two pads separated by 2 um gap
    let gnd_outer = vec![pt(-250e-6, -250e-6), pt(250e-6, -250e-6), pt(250e-6, 250e-6), pt(-250e-6, 250e-6)];
    let gnd_hole = vec![pt(-5e-6, -6.5e-6), pt(5e-6, -6.5e-6), pt(5e-6, 6.5e-6), pt(-5e-6, 6.5e-6)];
    let poly_gnd = m.polygons("poly_gnd", vec![gnd_outer], vec![gnd_hole]);

    let pad1_outer = vec![pt(-2.5e-6, -4e-6), pt(2.5e-6, -4e-6), pt(2.5e-6, -1e-6), pt(-2.5e-6, -1e-6)];
    let pad2_outer = vec![pt(-2.5e-6, 1e-6), pt(2.5e-6, 1e-6), pt(2.5e-6, 4e-6), pt(-2.5e-6, 4e-6)];
    let poly_p1 = m.polygons("poly_p1", vec![pad1_outer], vec![]);
    let poly_p2 = m.polygons("poly_p2", vec![pad2_outer], vec![]);

    m.sheet("gnd", poly_gnd, 0.0);
    m.sheet("pad_1", poly_p1, 0.0);
    m.sheet("pad_2", poly_p2, 0.0);
    m.overlay(&[("VacBot", vac_bot), ("Si", si_sub), ("VacTop", vac_top)]);
    m.size("sz_p1", Selector::parse("name(pad_1)").unwrap(), Some(Expr::Num(1.3e-6)), None, None);
    m.size("sz_p2", Selector::parse("name(pad_2)").unwrap(), Some(Expr::Num(1.3e-6)), None, None);

    let req_full = MeshRequest {
        h: 90e-6,
        max_dev: None,
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: false,
        gap_cells: Some(1.5),
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: false,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };
    let req_coarse = MeshRequest { z_coarsen: true, ..req_full.clone() };

    let (g_full, _) = mesh_request(&m.list, &req_full).unwrap();
    let (g_coarse, _) = mesh_request(&m.list, &req_coarse).unwrap();

    let expected_vol = 500e-6 * 500e-6 * (1500e-6 - (-1500e-6));
    let interiors = ["gnd", "pad_1", "pad_2", "Si|VacTop", "Si|VacBot"];
    verify_watertight_pure_tet_mesh(&g_full, expected_vol, &interiors);
    verify_watertight_pure_tet_mesh(&g_coarse, expected_vol, &interiors);

    // Near-sheet plane `z = 0` triangulation must remain 100% uncoarsened (`M_0`)
    for sname in ["gnd", "pad_1", "pad_2", "Si|VacTop"] {
        let n_full: usize = g_full.boundary(sname).into_iter().map(|l| g_full.mesh.sheet_facets(l).len()).sum();
        let n_coarse: usize = g_coarse.boundary(sname).into_iter().map(|l| g_coarse.mesh.sheet_facets(l).len()).sum();
        assert_eq!(n_full, n_coarse, "z=0 sheet/interface {sname} triangle count changed under z_coarsen");
    }

    let cell_ratio = g_full.mesh.n_cells() as f64 / g_coarse.mesh.n_cells() as f64;
    let p2_dof_ratio = g_full.mesh.n_dofs(2) as f64 / g_coarse.mesh.n_dofs(2) as f64;
    assert!(
        cell_ratio >= 3.0,
        "expected >= 3x cell reduction from z_coarsen, got {cell_ratio:.2}x (full={}, coarse={})",
        g_full.mesh.n_cells(),
        g_coarse.mesh.n_cells()
    );
    assert!(
        p2_dof_ratio >= 3.0,
        "expected >= 3x P2 DOF reduction from z_coarsen, got {p2_dof_ratio:.2}x (full={}, coarse={})",
        g_full.mesh.n_dofs(2),
        g_coarse.mesh.n_dofs(2)
    );

    // P1 and P2 Maxwell self-capacitance C_11 must agree to `< 0.2%` (`2e-3` relative)
    let eps_map = [("VacBot", 1.0), ("Si", 11.45), ("VacTop", 1.0)];
    let grounds = ["gnd", "pad_2", "VacTop.zmax", "VacBot.zmin"];

    let c11_p1_full = solve_terminal_capacitance_f(&g_full, 1, &eps_map, "pad_1", &grounds);
    let c11_p1_coarse = solve_terminal_capacitance_f(&g_coarse, 1, &eps_map, "pad_1", &grounds);
    let rel_p1 = (c11_p1_coarse - c11_p1_full).abs() / c11_p1_full;

    let c11_p2_full = solve_terminal_capacitance_f(&g_full, 2, &eps_map, "pad_1", &grounds);
    let c11_p2_coarse = solve_terminal_capacitance_f(&g_coarse, 2, &eps_map, "pad_1", &grounds);
    let rel_p2 = (c11_p2_coarse - c11_p2_full).abs() / c11_p2_full;
    eprintln!(
        "[transmon] cells {} -> {} ({:.2}x), P2 DOFs {} -> {} ({:.2}x), P1 dC/C={:.4}%, P2 dC/C={:.4}%",
        g_full.mesh.n_cells(),
        g_coarse.mesh.n_cells(),
        cell_ratio,
        g_full.mesh.n_dofs(2),
        g_coarse.mesh.n_dofs(2),
        p2_dof_ratio,
        rel_p1 * 100.0,
        rel_p2 * 100.0,
    );
    assert!(
        rel_p1 < 2e-3,
        "P1 capacitance drift {rel_p1:.6e} >= 0.2% (full={c11_p1_full:.6e} F, coarse={c11_p1_coarse:.6e} F)"
    );
    assert!(
        rel_p2 < 2e-3,
        "P2 capacitance drift {rel_p2:.6e} >= 0.2% (full={c11_p2_full:.6e} F, coarse={c11_p2_coarse:.6e} F)"
    );

    // Verify deterministic bit-identical mesh generation across thread pool sizes
    let m1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap().install(|| {
        mesh_request(&m.list, &req_coarse).unwrap().0.mesh
    });
    let m4 = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap().install(|| {
        mesh_request(&m.list, &req_coarse).unwrap().0.mesh
    });
    assert_eq!(m1.points.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), m4.points.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    assert_eq!(m1.cells, m4.cells);
    assert_eq!(m1.cell_labels, m4.cell_labels);
    assert_eq!(m1.facets, m4.facets);
    assert_eq!(m1.facet_labels, m4.facet_labels);
}

#[test]
fn squadds_hfss_v1_z_slab_controls_preserve_interfaces_within_15_slabs_and_abort_max_cells_under_50ms() {
    let mut m = Model::new_3d();
    // SQuADDS HFSS-style 3-layer stack: [-750..-350] um VacBot, [-350..0] um Si, [0..750] um VacTop
    let vac_bot = m.box_("vac_bot", (-300e-6, -300e-6, -750e-6), (300e-6, 300e-6, -350e-6));
    let si_sub = m.box_("si_sub", (-300e-6, -300e-6, -350e-6), (300e-6, 300e-6, 0.0));
    let vac_top = m.box_("vac_top", (-300e-6, -300e-6, 0.0), (300e-6, 300e-6, 750e-6));

    // Ground plane with [-110..110] x [-130..130] um pocket
    let gnd_outer = vec![pt(-300e-6, -300e-6), pt(300e-6, -300e-6), pt(300e-6, 300e-6), pt(-300e-6, 300e-6)];
    let gnd_hole = vec![pt(-110e-6, -130e-6), pt(110e-6, -130e-6), pt(110e-6, 130e-6), pt(-110e-6, 130e-6)];
    let poly_gnd = m.polygons("poly_gnd", vec![gnd_outer], vec![gnd_hole]);

    // Two shunt pads + readout claw + sub-micron 200 nm JJ strap at z = 0
    let pad_a = m.polygons("poly_pa", vec![vec![pt(-65e-6, -68e-6), pt(65e-6, -68e-6), pt(65e-6, -12e-6), pt(-65e-6, -12e-6)]], vec![]);
    let pad_b = m.polygons("poly_pb", vec![vec![pt(-65e-6, 12e-6), pt(65e-6, 12e-6), pt(65e-6, 68e-6), pt(-65e-6, 68e-6)]], vec![]);
    let jj = m.polygons("poly_jj", vec![vec![pt(-0.1e-6, -12e-6), pt(0.1e-6, -12e-6), pt(0.1e-6, 12e-6), pt(-0.1e-6, 12e-6)]], vec![]);
    let claw = m.polygons("poly_claw", vec![vec![pt(-55e-6, 86e-6), pt(55e-6, 86e-6), pt(55e-6, 98e-6), pt(-55e-6, 98e-6)]], vec![]);

    m.sheet("gnd", poly_gnd, 0.0);
    m.sheet("pad_a", pad_a, 0.0);
    m.sheet("pad_b", pad_b, 0.0);
    m.sheet("jj", jj, 0.0);
    m.sheet("claw", claw, 0.0);
    m.overlay(&[("VacBot", vac_bot), ("Si", si_sub), ("VacTop", vac_top)]);
    // Sub-micron 200 nm JJ edge size rule that previously forced dz = 200 nm across the entire die
    m.size("sz_jj", Selector::parse("jj:edge").unwrap(), Some(Expr::Num(0.2e-6)), None, None);

    let count_z_slabs = |g: &MeshedGeom| -> (usize, Vec<f64>) {
        let mut zs: Vec<f64> = g.mesh.points.as_chunks::<3>().0.iter().map(|p| p[2]).collect();
        zs.sort_by(|a, b| a.total_cmp(b));
        zs.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
        let n_slabs = zs.len().saturating_sub(1);
        (n_slabs, zs)
    };

    let base_req = MeshRequest {
        h: 100e-6,
        max_dev: None,
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: false,
        gap_cells: None,
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: true,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };

    // 1. Unconstrained baseline produces > 30 z-slabs due to the 200 nm JJ edge
    let (g_unconstrained, _) = mesh_request(&m.list, &base_req).unwrap();
    let (slabs_unconstrained, _) = count_z_slabs(&g_unconstrained);
    assert!(slabs_unconstrained > 30, "expected > 30 unconstrained z-slabs, got {slabs_unconstrained}");

    // 2. Decoupled Z-slab controls (`dz_min`, `z_growth_ratio`, `max_z_slabs = 15`) cap `N_slabs <= 15`
    //    while strictly preserving all stack interface elevations (-750, -350, 0, +750 um)
    let expected_vol = 600e-6 * 600e-6 * 1500e-6;
    let interiors = ["gnd", "pad_a", "pad_b", "jj", "claw", "Si|VacTop", "Si|VacBot"];
    for z_coarsen in [false, true] {
        let controlled_req = MeshRequest {
            dz_min: Some(6.0e-6),
            z_growth_ratio: Some(1.65),
            max_z_slabs: Some(15),
            z_coarsen,
            ..base_req.clone()
        };
        let (g_ctrl, _) = mesh_request(&m.list, &controlled_req).unwrap();
        let (n_slabs, zs) = count_z_slabs(&g_ctrl);
        assert!(n_slabs <= 15, "z_coarsen={z_coarsen}: expected N_slabs <= 15, got {n_slabs} ({zs:?})");
        for z_iface in [-750e-6, -350e-6, 0.0, 750e-6] {
            assert!(
                zs.iter().any(|&z| (z - z_iface).abs() <= 1e-12),
                "z_coarsen={z_coarsen}: missing stack interface z={z_iface:e} in {zs:?}"
            );
        }
        verify_watertight_pure_tet_mesh(&g_ctrl, expected_vol, &interiors);
    }

    // 3. `max_cells` pre-allocation guard aborts in `< 50 ms` with diagnostic `(n_tri2, n_slabs, projected_tets, max_cells)`
    for z_coarsen in [false, true] {
        let abort_req = MeshRequest {
            max_cells: Some(2_000),
            z_coarsen,
            ..base_req.clone()
        };
        let t0 = std::time::Instant::now();
        let err = mesh_request(&m.list, &abort_req).err().expect("projected_tets > max_cells must abort");
        let abort_ms = t0.elapsed().as_secs_f64() * 1000.0;
        assert!(abort_ms < 50.0, "z_coarsen={z_coarsen}: max_cells abort took {abort_ms:.2} ms (limit 50 ms)");
        assert!(
            err.contains("n_tri2=") && err.contains("n_slabs=") && err.contains("projected_tets=") && err.contains("max_cells=2000"),
            "missing diagnostic fields in error: {err}"
        );
    }

    // Also verify `Adaptive::max_cells` is honoured by the pre-allocation guard
    let adapt_abort = MeshRequest {
        mesher: MesherRequest::Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(1_500) },
        ..base_req.clone()
    };
    let err_adapt = mesh_request(&m.list, &adapt_abort).err().expect("Adaptive::max_cells must abort");
    assert!(err_adapt.contains("max_cells=1500"), "{err_adapt}");

    // 4. Invalid option values are rejected cleanly
    assert!(mesh_request(&m.list, &MeshRequest { dz_min: Some(0.0), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { dz_min: Some(-1e-6), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { z_growth_ratio: Some(1.0), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { z_growth_ratio: Some(0.9), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { max_z_slabs: Some(0), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { max_z_slabs: Some(2), ..base_req.clone() }).is_err());
    assert!(mesh_request(&m.list, &MeshRequest { max_cells: Some(0), ..base_req }).is_err());
}
