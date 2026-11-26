// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME Layer 2 Pillar (a): Exact analytical closed-form verification benchmarks across all `nyx-dfs`
//! formulations (2D/3D electrostatics with curved P2 spheres and `Open`/`FarField` truncation at `R_box = 5a`,
//! magnetostatics `cross_section`/`planar_stream`/`mixed_3d` coax, circular loop + cross-formulation agreement,
//! parallel-strip transmission line, and 1D/3D steady convection & 1D transient eigenmode heat decay).

use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use geom::model::{Model, RegionKind};
use nyx_dfs::{solve_electrostatics, solve_heat, solve_magnetostatics, Bc, BcKind, Problem, Study, EPS0, MU0};
use nyx_physics::InductanceMatrixData;
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn elliptic_k(k: f64) -> f64 {
    let mut a = 1.0_f64;
    let mut b = (1.0 - k * k).sqrt();
    for _ in 0..16 {
        let an = 0.5 * (a + b);
        let bn = (a * b).sqrt();
        a = an;
        b = bn;
        if (a - b).abs() < 1e-15 {
            break;
        }
    }
    0.5 * PI / a
}

/// Builds a pure-`Tet` mesh of the spherical shell `a <= r <= b` with quadratic curved spherical
/// edge midpoints projected onto every concentric radial shell (`r_k`).
fn build_curved_spherical_shell(r_in: f64, r_out: f64, m: usize, n_r: usize) -> SimplexMesh {
    let mut grid_to_surf: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut u_verts: Vec<[f64; 3]> = Vec::new();
    let mut get_surf = |ix: usize, iy: usize, iz: usize| -> usize {
        if let Some(&id) = grid_to_surf.get(&(ix, iy, iz)) {
            return id;
        }
        let id = u_verts.len();
        let q = [
            (0.25 * PI * (2.0 * (ix as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iy as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iz as f64) / (m as f64) - 1.0)).tan(),
        ];
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        u_verts.push([q[0] / len, q[1] / len, q[2] / len]);
        grid_to_surf.insert((ix, iy, iz), id);
        id
    };

    let mut quads: Vec<[usize; 4]> = Vec::with_capacity(6 * m * m);
    for i in 0..m {
        for j in 0..m {
            quads.push([get_surf(m, i, j), get_surf(m, i + 1, j), get_surf(m, i + 1, j + 1), get_surf(m, i, j + 1)]);
            quads.push([get_surf(0, j, i), get_surf(0, j + 1, i), get_surf(0, j + 1, i + 1), get_surf(0, j, i + 1)]);
            quads.push([get_surf(j, m, i), get_surf(j + 1, m, i), get_surf(j + 1, m, i + 1), get_surf(j, m, i + 1)]);
            quads.push([get_surf(i, 0, j), get_surf(i + 1, 0, j), get_surf(i + 1, 0, j + 1), get_surf(i, 0, j + 1)]);
            quads.push([get_surf(i, j, m), get_surf(i + 1, j, m), get_surf(i + 1, j + 1, m), get_surf(i, j + 1, m)]);
            quads.push([get_surf(j, i, 0), get_surf(j + 1, i, 0), get_surf(j + 1, i + 1, 0), get_surf(j, i + 1, 0)]);
        }
    }

    let n_s = u_verts.len();
    let mut cm = CellMesh::<3>::new();
    for k in 0..=n_r {
        let r = r_in + (r_out - r_in) * (k as f64) / (n_r as f64);
        for u in &u_verts {
            cm.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }
    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_r {
        for &[s0, s1, s2, s3] in &quads {
            cm.push_cell(
                CellType::Hex,
                &[
                    vid(k, s0),
                    vid(k, s1),
                    vid(k, s2),
                    vid(k, s3),
                    vid(k + 1, s0),
                    vid(k + 1, s1),
                    vid(k + 1, s2),
                    vid(k + 1, s3),
                ],
                1,
            );
        }
    }
    for &[s0, s1, s2, s3] in &quads {
        cm.push_facet(&[vid(0, s0), vid(0, s3), vid(0, s2), vid(0, s1)], 1);
        cm.push_facet(&[vid(n_r, s0), vid(n_r, s1), vid(n_r, s2), vid(n_r, s3)], 2);
    }

    let mut sm = SimplexMesh::from_cells(&cm).unwrap();
    let mut en = Vec::new();
    for &(va, vb) in sm.edges() {
        let pa = sm.point(va as usize);
        let pb = sm.point(vb as usize);
        let ra = (pa[0] * pa[0] + pa[1] * pa[1] + pa[2] * pa[2]).sqrt();
        let rb = (pb[0] * pb[0] + pb[1] * pb[1] + pb[2] * pb[2]).sqrt();
        if (ra - rb).abs() < 1e-11 {
            let mid: [f64; 3] = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            let rmid = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]).sqrt();
            if rmid > 1e-15 {
                let sc = ra / rmid;
                en.push(([va, vb], [mid[0] * sc, mid[1] * sc, mid[2] * sc]));
            }
        }
    }
    sm.set_edge_nodes(&en);
    sm
}

/// Builds a pure-`Tet` mesh between an inner sphere `r = a` (label 1) and an outer cube `[-R_box, R_box]^3` (label 2).
fn build_curved_sphere_in_box(a: f64, r_box: f64, m: usize, n_r: usize) -> SimplexMesh {
    let mut grid_to_surf: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut box_pts: Vec<[f64; 3]> = Vec::new();
    let mut sph_pts: Vec<[f64; 3]> = Vec::new();
    let mut get_surf = |ix: usize, iy: usize, iz: usize| -> usize {
        if let Some(&id) = grid_to_surf.get(&(ix, iy, iz)) {
            return id;
        }
        let id = box_pts.len();
        let cx = 2.0 * (ix as f64) / (m as f64) - 1.0;
        let cy = 2.0 * (iy as f64) / (m as f64) - 1.0;
        let cz = 2.0 * (iz as f64) / (m as f64) - 1.0;
        box_pts.push([r_box * cx, r_box * cy, r_box * cz]);
        let q = [(0.25 * PI * cx).tan(), (0.25 * PI * cy).tan(), (0.25 * PI * cz).tan()];
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        sph_pts.push([a * q[0] / len, a * q[1] / len, a * q[2] / len]);
        grid_to_surf.insert((ix, iy, iz), id);
        id
    };

    let mut quads: Vec<[usize; 4]> = Vec::with_capacity(6 * m * m);
    for i in 0..m {
        for j in 0..m {
            quads.push([get_surf(m, i, j), get_surf(m, i + 1, j), get_surf(m, i + 1, j + 1), get_surf(m, i, j + 1)]);
            quads.push([get_surf(0, j, i), get_surf(0, j + 1, i), get_surf(0, j + 1, i + 1), get_surf(0, j, i + 1)]);
            quads.push([get_surf(j, m, i), get_surf(j + 1, m, i), get_surf(j + 1, m, i + 1), get_surf(j, m, i + 1)]);
            quads.push([get_surf(i, 0, j), get_surf(i + 1, 0, j), get_surf(i + 1, 0, j + 1), get_surf(i, 0, j + 1)]);
            quads.push([get_surf(i, j, m), get_surf(i + 1, j, m), get_surf(i + 1, j + 1, m), get_surf(i, j + 1, m)]);
            quads.push([get_surf(j, i, 0), get_surf(j + 1, i, 0), get_surf(j + 1, i + 1, 0), get_surf(j, i + 1, 0)]);
        }
    }

    let n_s = sph_pts.len();
    let mut cm = CellMesh::<3>::new();
    for k in 0..=n_r {
        let u = (k as f64) / (n_r as f64);
        let t = u * u;
        for s in 0..n_s {
            cm.points.push([
                (1.0 - t) * sph_pts[s][0] + t * box_pts[s][0],
                (1.0 - t) * sph_pts[s][1] + t * box_pts[s][1],
                (1.0 - t) * sph_pts[s][2] + t * box_pts[s][2],
            ]);
        }
    }
    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_r {
        for &[s0, s1, s2, s3] in &quads {
            cm.push_cell(
                CellType::Hex,
                &[
                    vid(k, s0),
                    vid(k, s1),
                    vid(k, s2),
                    vid(k, s3),
                    vid(k + 1, s0),
                    vid(k + 1, s1),
                    vid(k + 1, s2),
                    vid(k + 1, s3),
                ],
                1,
            );
        }
    }
    for &[s0, s1, s2, s3] in &quads {
        cm.push_facet(&[vid(0, s0), vid(0, s3), vid(0, s2), vid(0, s1)], 1);
        cm.push_facet(&[vid(n_r, s0), vid(n_r, s1), vid(n_r, s2), vid(n_r, s3)], 2);
    }

    let mut sm = SimplexMesh::from_cells(&cm).unwrap();
    let mut en = Vec::new();
    for &(va, vb) in sm.edges() {
        if (va as usize) < n_s && (vb as usize) < n_s {
            let pa = sm.point(va as usize);
            let pb = sm.point(vb as usize);
            let mid: [f64; 3] = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            let rmid = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]).sqrt();
            if rmid > 1e-15 {
                let sc = a / rmid;
                en.push(([va, vb], [mid[0] * sc, mid[1] * sc, mid[2] * sc]));
            }
        }
    }
    sm.set_edge_nodes(&en);
    sm
}

#[test]
fn electro_closed_forms_2d_coax_3d_curved_spheres_and_5a_open_farfield() {
    // 1. 2D coaxial cylinder: C' = 2 pi eps / ln(b/a)
    let (a2, b2, eps2) = (0.25_f64, 1.0_f64, 3.2_f64);
    let mut m2 = Model::new_2d();
    let outer = m2.circle("diel", (0.0, 0.0), b2);
    let inner = m2.circle("core", (0.0, 0.0), a2);
    m2.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let mut opts = GeomMeshOpts::new(0.05);
    opts.max_dev = Some(0.05 / 60.0);
    let g2 = mesh_geom_2d(evaluate(&m2.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (a2 + b2);
    let mut in_labs = Vec::new();
    let mut out_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        if rm < r_mid && !in_labs.contains(&lab) {
            in_labs.push(lab);
        } else if rm >= r_mid && !out_labs.contains(&lab) {
            out_labs.push(lab);
        }
    }
    let mut p2 = Problem::new(Study::Electrostatic, 2);
    p2.materials.default.insert("eps_r".into(), eps2);
    p2.bcs = vec![Bc { labels: in_labs, kind: BcKind::Potential(1.0) }, Bc { labels: out_labs, kind: BcKind::Ground }];
    let res2 = solve_electrostatics(&g2.mesh, &p2).unwrap();
    let c2_num = res2.scalar("capacitance").unwrap();
    let c2_exact = 2.0 * PI * EPS0 * eps2 / (b2 / a2).ln();
    let err2 = (c2_num - c2_exact).abs() / c2_exact;
    assert!(err2 < 2e-3, "2D coax C' err {err2:e} >= 0.2%");

    // 2. 3D concentric spheres on CurvedTet P2 mesh: C = 4 pi eps a b / (b - a)
    let (a3, b3, eps3) = (0.4_f64, 1.2_f64, 2.5_f64);
    let sph_curved = build_curved_spherical_shell(a3, b3, 4, 3);
    assert!(sph_curved.has_curved());
    let mut p3 = Problem::new(Study::Electrostatic, 2);
    p3.materials.default.insert("eps_r".into(), eps3);
    p3.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![2], kind: BcKind::Ground }];
    let c3_exact = 4.0 * PI * EPS0 * eps3 * a3 * b3 / (b3 - a3);
    let c3_curved = solve_electrostatics(&sph_curved, &p3).unwrap().scalar("capacitance").unwrap();
    let err3_curved = (c3_curved - c3_exact).abs() / c3_exact;
    assert!(err3_curved < 2e-3, "3D curved P2 concentric spheres err {err3_curved:e} >= 0.2%");

    let mut sph_straight = sph_curved.clone();
    sph_straight.set_edge_nodes(&[]);
    let c3_straight = solve_electrostatics(&sph_straight, &p3).unwrap().scalar("capacitance").unwrap();
    let err3_straight = (c3_straight - c3_exact).abs() / c3_exact;
    assert!(
        err3_curved < 0.25 * err3_straight,
        "curved P2 ({err3_curved:e}) must sharply beat straight P2 ({err3_straight:e}) on concentric spheres"
    );

    // 3. Isolated sphere capacitance C = 4 pi eps0 a in a R_box = 5a cube with Open and FarField (< 1%)
    let a_iso = 0.2_f64;
    let r_box = 5.0 * a_iso;
    let iso_mesh = build_curved_sphere_in_box(a_iso, r_box, 5, 5);
    let c_iso_exact = 4.0 * PI * EPS0 * a_iso;
    for (name, outer_bc) in [
        ("Open", BcKind::Open { origin: Some([0.0; 3]) }),
        ("FarField", BcKind::FarField { order: 2, decay_power: 1.0, center: Some([0.0; 3]) }),
    ] {
        let mut p_iso = Problem::new(Study::Electrostatic, 2);
        p_iso.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![2], kind: outer_bc }];
        let c_num = solve_electrostatics(&iso_mesh, &p_iso).unwrap().scalar("capacitance").unwrap();
        let rel = (c_num - c_iso_exact).abs() / c_iso_exact;
        eprintln!("Isolated sphere 5a ({name}): rel_err = {:.4}%", rel * 100.0);
        assert!(rel < 0.01, "Isolated sphere 5a ({name}) rel_err {:.3}% exceeds 1%", rel * 100.0);
    }
}

#[test]
fn magneto_closed_forms_coax_loop_parallel_strips_and_cross_formulation_agreement() {
    // 1. 2D coaxial cable L' = (mu0 / 2 pi) ln(b / a)
    let (a, b) = (0.2_f64, 0.8_f64);
    let mut m_coax = Model::new_2d();
    let outer = m_coax.circle("diel", (0.0, 0.0), b);
    let inner = m_coax.circle("core", (0.0, 0.0), a);
    m_coax.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let mut opts = GeomMeshOpts::new(0.035);
    opts.max_dev = Some(0.035 / 80.0);
    let g_coax = mesh_geom_2d(evaluate(&m_coax.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (a + b);
    let (mut in_labs, mut out_labs) = (Vec::new(), Vec::new());
    for f in 0..g_coax.mesh.n_facets() {
        let fv = g_coax.mesh.facet(f);
        let p0 = g_coax.mesh.point(fv[0] as usize);
        let p1 = g_coax.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g_coax.mesh.facet_labels[f];
        if rm < r_mid && !in_labs.contains(&lab) {
            in_labs.push(lab);
        } else if rm >= r_mid && !out_labs.contains(&lab) {
            out_labs.push(lab);
        }
    }
    let mut p_coax = Problem::new(Study::Magnetostatic { formulation: "cross_section".into() }, 2);
    p_coax.bcs = vec![
        Bc {
            labels: in_labs,
            kind: BcKind::CurrentTerminal { name: "center".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc { labels: out_labs, kind: BcKind::Ground },
    ];
    let l_coax = InductanceMatrixData::of(&solve_magnetostatics(&g_coax.mesh, &p_coax).unwrap())
        .unwrap()
        .inductance_matrix_h[0][0];
    let l_coax_exact = (MU0 / (2.0 * PI)) * (b / a).ln();
    assert!((l_coax - l_coax_exact).abs() / l_coax_exact < 3e-3);

    // 2. Thin circular superconducting loop: L ~ mu0 R [ln(8R/a_GMD) - 2], with flat width w = 4 a_GMD
    //    AND direct cross-formulation agreement (< 2%) between `planar_stream` and `mixed_3d`.
    let r_loop = 100.0e-6;
    let w_ring = 20.0e-6;
    let r_in = r_loop - 0.5 * w_ring;
    let r_out = r_loop + 0.5 * w_ring;
    let a_gmd = w_ring / 4.0;
    let l_loop_exact = MU0 * r_loop * ((8.0 * r_loop / a_gmd).ln() - 2.0);

    // 2a. planar_stream annulus
    let mut m_loop2 = Model::new_2d();
    let outer_c = m_loop2.circle("washer", (0.0, 0.0), r_out);
    let inner_c = m_loop2.circle("hole", (0.0, 0.0), r_in);
    m_loop2.overlay_full(
        &[("washer", outer_c, RegionKind::Material, None), ("hole", inner_c, RegionKind::Hole, None)],
        &[],
    );
    let mut opts_l = GeomMeshOpts::new(3.5e-6);
    opts_l.max_dev = Some(3.5e-6 / 60.0);
    let g_loop2 = mesh_geom_2d(evaluate(&m_loop2.list).unwrap(), &opts_l).unwrap();
    let (mut hole_labs, mut outer_labs) = (Vec::new(), Vec::new());
    for f in 0..g_loop2.mesh.n_facets() {
        let fv = g_loop2.mesh.facet(f);
        let p0 = g_loop2.mesh.point(fv[0] as usize);
        let p1 = g_loop2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g_loop2.mesh.facet_labels[f];
        if rm < r_loop && !hole_labs.contains(&lab) {
            hole_labs.push(lab);
        } else if rm >= r_loop && !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }
    let mut p_ps = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_ps.bcs = vec![
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
        Bc { labels: hole_labs.clone(), kind: BcKind::FluxLoop { name: "loop".into(), flux_wb: 1.0 } },
    ];
    let l_ps = InductanceMatrixData::of(&solve_magnetostatics(&g_loop2.mesh, &p_ps).unwrap())
        .unwrap()
        .inductance_matrix_h[0][0];
    let err_ps = (l_ps - l_loop_exact).abs() / l_loop_exact;
    eprintln!("circular loop planar_stream: L={l_ps:.6e}, exact={l_loop_exact:.6e}, rel_err={:.3}%", err_ps * 100.0);
    assert!(err_ps < 0.025, "planar_stream circular loop rel_err {:.3}% too large", err_ps * 100.0);

    // 2b. mixed_3d radial-graded Tet+Prism InfiniteElement spherical shell around z = 0 annulus [r_in, r_out]
    let g_loop3_mesh = {
        use femcore::mesh::cells::{CellMesh, CellType};
        use femcore::props::PmlDir;
        use nyx_dfs::{IeCoordKind, IeProfile, InfiniteElementSpec};
        use nyx_physics::Mesh;
        let n_lat = 18usize; // even so j_eq = 9 sits at exact z = 0, graded toward equator
        let n_lon = 32usize;
        let j_eq = n_lat / 2;
        let n_ring = 5usize;
        let mut radii = vec![0.40 * r_in, 0.78 * r_in, r_in - 2.0e-6];
        let ring_start_k = radii.len(); // index of r_in in `radii`
        for k in 0..=n_ring {
            let u = 0.5 * (1.0 - (PI * k as f64 / n_ring as f64).cos());
            radii.push(r_in + u * (r_out - r_in));
        }
        for &r_extra in &[r_out + 2.0e-6, 1.30 * r_out, 1.85 * r_out, 2.5 * r_out] {
            radii.push(r_extra);
        }
        let n_phys = radii.len() - 1;
        let r_phys = *radii.last().unwrap();
        let r_ie_out = 5.0 * r_out;
        let n_ie = 4usize;
        for k in 1..=n_ie {
            radii.push(r_phys + (k as f64 / n_ie as f64) * (r_ie_out - r_phys));
        }
        let n_tot = n_phys + n_ie;

        let mut u_verts = vec![[0.0, 0.0, -1.0]];
        for j in 1..n_lat {
            let t = -1.0 + 2.0 * (j as f64) / (n_lat as f64);
            let th = 0.5 * PI * t.signum() * t.abs().powf(2.0);
            let (ct, st) = (th.cos(), th.sin());
            for i in 0..n_lon {
                let phi = 2.0 * PI * (i as f64) / (n_lon as f64);
                u_verts.push([ct * phi.cos(), ct * phi.sin(), st]);
            }
        }
        let north_id = u_verts.len();
        u_verts.push([0.0, 0.0, 1.0]);
        let n_s = u_verts.len();
        let ring_v = |j: usize, i: usize| -> usize {
            if j == 0 {
                0
            } else if j == n_lat {
                north_id
            } else {
                1 + (j - 1) * n_lon + (i % n_lon)
            }
        };
        let mut tris = Vec::new();
        for i in 0..n_lon {
            tris.push([0, ring_v(1, i + 1), ring_v(1, i)]);
        }
        for j in 1..(n_lat - 1) {
            for i in 0..n_lon {
                let (a, b, c, d) = (ring_v(j, i), ring_v(j, i + 1), ring_v(j + 1, i + 1), ring_v(j + 1, i));
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
        for i in 0..n_lon {
            tris.push([ring_v(n_lat - 1, i), ring_v(n_lat - 1, i + 1), north_id]);
        }

        let mut cm = CellMesh::<3>::new();
        cm.points.push([0.0, 0.0, 0.0]); // vertex 0 = origin filling core 0..0.35*r_in
        for &r in &radii {
            for u in &u_verts {
                cm.points.push([r * u[0], r * u[1], r * u[2]]);
            }
        }
        let vid = |k: usize, s: usize| -> u32 { (1 + k * n_s + s) as u32 };
        for &[s0, s1, s2] in &tris {
            cm.push_cell(CellType::Tet, &[0, vid(0, s0), vid(0, s1), vid(0, s2)], 1);
        }
        for k in 0..n_tot {
            let dom = if k < n_phys { 1 } else { 2 };
            if k < n_phys {
                for &[s0, s1, s2] in &tris {
                    let mut pairs = [
                        (s0, vid(k, s0), vid(k + 1, s0)),
                        (s1, vid(k, s1), vid(k + 1, s1)),
                        (s2, vid(k, s2), vid(k + 1, s2)),
                    ];
                    pairs.sort_by_key(|p| p.0);
                    let bot = [pairs[0].1, pairs[1].1, pairs[2].1];
                    let top = [pairs[0].2, pairs[1].2, pairs[2].2];
                    cm.push_cell(CellType::Tet, &[bot[0], bot[1], bot[2], top[2]], dom);
                    cm.push_cell(CellType::Tet, &[bot[0], bot[1], top[1], top[2]], dom);
                    cm.push_cell(CellType::Tet, &[bot[0], top[0], top[1], top[2]], dom);
                }
            } else {
                for &[s0, s1, s2] in &tris {
                    cm.push_cell(
                        CellType::Prism,
                        &[vid(k, s0), vid(k, s1), vid(k, s2), vid(k + 1, s0), vid(k + 1, s1), vid(k + 1, s2)],
                        dom,
                    );
                }
            }
        }
        for &[s0, s1, s2] in &tris {
            cm.push_facet(&[vid(n_tot, s0), vid(n_tot, s1), vid(n_tot, s2)], 20);
        }
        for k in ring_start_k..(ring_start_k + n_ring) {
            for i in 0..n_lon {
                let (sa, sb) = (ring_v(j_eq, i), ring_v(j_eq, i + 1));
                let (ua, ub) = if sa < sb { (sa, sb) } else { (sb, sa) };
                let (b0, b1) = (vid(k, ua), vid(k, ub));
                let (t0, t1) = (vid(k + 1, ua), vid(k + 1, ub));
                cm.push_facet(&[b0, b1, t1], 30);
                cm.push_facet(&[b0, t0, t1], 30);
            }
        }
        let mut p_m3 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
        p_m3.bcs = vec![
            Bc { labels: vec![30], kind: BcKind::LondonSheet { l_k_sq: 0.0 } },
            Bc { labels: vec![30], kind: BcKind::FluxLoop { name: "loop".into(), flux_wb: 1.0 } },
            Bc {
                labels: vec![2],
                kind: BcKind::InfiniteElement(InfiniteElementSpec {
                    coord: IeCoordKind::Spherical,
                    profile: IeProfile::Rational,
                    x: PmlDir::Off,
                    y: PmlDir::Off,
                    z: PmlDir::Off,
                    center: Some([0.0; 3]),
                    axis: [0.0, 0.0, 1.0],
                    r_in: Some(r_phys),
                    thickness: Some(r_ie_out - r_phys),
                    pole_dist: Some(r_phys),
                    decay: 30.0,
                    bounds_lo: None,
                    bounds_hi: None,
                }),
            },
        ];
        let ex_mesh = Mesh::D3(cm);
        // Solve the superconducting loop with thin-film kinetic sheet inductance L_k,sq = 60 pH/sq
        // and verify direct cross-formulation agreement (< 2%) between `planar_stream` and `mixed_3d`.
        let l_k_sc = 6.0e-11;
        let mut p_m3_sc = p_m3.clone();
        p_m3_sc.bcs[0] = Bc { labels: vec![30], kind: BcKind::LondonSheet { l_k_sq: l_k_sc } };
        let r_sc = nyx_dfs::magneto::solve_magnetostatics_mixed(&ex_mesh, &p_m3_sc).unwrap();

        let mut p_ps_sc = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
        p_ps_sc.bcs = vec![
            Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
            Bc { labels: hole_labs.clone(), kind: BcKind::FluxLoop { name: "loop".into(), flux_wb: 1.0 } },
            Bc { labels: outer_labs, kind: BcKind::LondonSheet { l_k_sq: l_k_sc } },
        ];
        let l_ps_sc = InductanceMatrixData::of(&solve_magnetostatics(&g_loop2.mesh, &p_ps_sc).unwrap())
            .unwrap()
            .inductance_matrix_h[0][0];
        let l_m3_sc = InductanceMatrixData::of(&r_sc).unwrap().inductance_matrix_h[0][0];
        let sc_cross = (l_m3_sc - l_ps_sc).abs() / l_ps_sc;
        eprintln!(
            "thin SC loop (L_k=60pH/sq): planar_stream={l_ps_sc:.6e}, mixed_3d={l_m3_sc:.6e}, cross={:.3}%",
            sc_cross * 100.0
        );
        sc_cross
    };
    assert!(g_loop3_mesh < 0.02, "planar_stream vs mixed_3d loop disagreement {:.3}% exceeds 2%", g_loop3_mesh * 100.0);

    // 3. Parallel-strip transmission line (`planar_stream` & `mixed_3d` vs L' = mu0 K(k) / K(k'))
    let a_st = 10.0e-6;
    let b_st = 30.0e-6;
    let ell = 500.0e-6;
    let k_mod: f64 = a_st / b_st;
    let kp_mod = (1.0_f64 - k_mod * k_mod).sqrt();
    let l_line_exact = ell * MU0 * elliptic_k(k_mod) / elliptic_k(kp_mod);

    let mut m_st2 = Model::new_2d();
    let st_l = m_st2.rect("s1", (-b_st, -0.5 * ell), (-a_st, 0.5 * ell));
    let st_r = m_st2.rect("s2", (a_st, -0.5 * ell), (b_st, 0.5 * ell));
    m_st2.overlay_full(&[("s1", st_l, RegionKind::Material, None), ("s2", st_r, RegionKind::Material, None)], &[]);
    let g_st2 = mesh_geom_2d(evaluate(&m_st2.list).unwrap(), &GeomMeshOpts::new(6.0e-6)).unwrap();
    let (mut gnd_st2, mut t1_st2, mut t2_st2) = (Vec::new(), Vec::new(), Vec::new());
    for f in 0..g_st2.mesh.n_facets() {
        let fv = g_st2.mesh.facet(f);
        let xm = 0.5 * (g_st2.mesh.point(fv[0] as usize)[0] + g_st2.mesh.point(fv[1] as usize)[0]);
        let lab = g_st2.mesh.facet_labels[f];
        if (xm + b_st).abs() < 1e-8 || (xm - b_st).abs() < 1e-8 {
            if !gnd_st2.contains(&lab) {
                gnd_st2.push(lab);
            }
        } else if (xm + a_st).abs() < 1e-8 {
            if !t1_st2.contains(&lab) {
                t1_st2.push(lab);
            }
        } else if (xm - a_st).abs() < 1e-8 && !t2_st2.contains(&lab) {
            t2_st2.push(lab);
        }
    }
    let mut p_st2 = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_st2.bcs = vec![
        Bc { labels: gnd_st2, kind: BcKind::Ground },
        Bc {
            labels: t1_st2,
            kind: BcKind::CurrentTerminal { name: "s1".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc {
            labels: t2_st2,
            kind: BcKind::CurrentTerminal { name: "s2".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
    ];
    let lmat_st2 =
        InductanceMatrixData::of(&solve_magnetostatics(&g_st2.mesh, &p_st2).unwrap()).unwrap().inductance_matrix_h;
    // Notice opposite x-slopes give +I on s1 and -I on s2 when excited with [+1, +1]
    let l_diff_st2 = lmat_st2[0][0] + lmat_st2[1][1] + 2.0 * lmat_st2[0][1];
    let rel_st2 = (l_diff_st2 - l_line_exact).abs() / l_line_exact;
    eprintln!(
        "parallel strips planar_stream: L={l_diff_st2:.6e}, exact={l_line_exact:.6e}, rel_err={:.2}%",
        rel_st2 * 100.0
    );
    assert!(rel_st2 < 0.06, "planar_stream parallel strips err {:.2}% too large", rel_st2 * 100.0);

    // 3b. mixed_3d extruded 2D isotropic edge-refined cross-section with nested inner refinement box
    let l_diff_st3 = {
        use geom::expr::Expr;
        use geom::select::Selector;
        let mut m_xsec = Model::new_2d();
        let outer = m_xsec.rect("outer", (-500e-6, -400e-6), (500e-6, 400e-6));
        let ib = m_xsec.rect("ib", (-60e-6, -40e-6), (60e-6, 0.0));
        let itl = m_xsec.rect("itl", (-60e-6, 0.0), (-b_st, 40e-6));
        let s1r = m_xsec.rect("s1", (-b_st, 0.0), (-a_st, 40e-6));
        let tm = m_xsec.rect("tm", (-a_st, 0.0), (a_st, 40e-6));
        let s2r = m_xsec.rect("s2", (a_st, 0.0), (b_st, 40e-6));
        let itr = m_xsec.rect("itr", (b_st, 0.0), (60e-6, 40e-6));
        m_xsec.overlay_full(
            &[
                ("outer", outer, RegionKind::Material, Some(10)),
                ("ib", ib, RegionKind::Material, Some(1)),
                ("itl", itl, RegionKind::Material, Some(2)),
                ("s1", s1r, RegionKind::Material, Some(3)),
                ("tm", tm, RegionKind::Material, Some(4)),
                ("s2", s2r, RegionKind::Material, Some(5)),
                ("itr", itr, RegionKind::Material, Some(6)),
            ],
            &[],
        );
        for reg in ["ib", "itl", "s1", "tm", "s2", "itr"] {
            m_xsec.size(&format!("sz_{reg}"), Selector::Region(reg.into()), Some(Expr::Num(2.0e-6)), None, None);
        }
        m_xsec.size(
            "sz_tm_if",
            Selector::Interface("ib".into(), Some("tm".into())),
            Some(Expr::Num(1.2e-6)),
            None,
            None,
        );
        m_xsec.size(
            "sz_if_s1",
            Selector::Interface("ib".into(), Some("s1".into())),
            Some(Expr::Num(0.45e-6)),
            None,
            None,
        );
        m_xsec.size(
            "sz_if_s2",
            Selector::Interface("ib".into(), Some("s2".into())),
            Some(Expr::Num(0.45e-6)),
            None,
            None,
        );
        let g2_tem = mesh_geom_2d(evaluate(&m_xsec.list).unwrap(), &GeomMeshOpts::new(50.0e-6)).unwrap();

        let m2 = &g2_tem.mesh;
        let nv2 = m2.n_points();
        let mut pts3 = Vec::with_capacity(nv2 * 2 * 3);
        for &y_val in &[-0.5 * ell, 0.5 * ell] {
            for v in 0..nv2 {
                let p = m2.point(v);
                pts3.extend_from_slice(&[p[0], y_val, p[1]]);
            }
        }
        let mut cells3 = Vec::with_capacity(m2.n_cells() * 3 * 4);
        let mut clabs3 = Vec::with_capacity(m2.n_cells() * 3);
        let mut facets3 = Vec::new();
        let mut flabs3 = Vec::new();
        let mut seen_strip_edges: std::collections::BTreeSet<(u32, u32)> = std::collections::BTreeSet::new();
        for c in 0..m2.n_cells() {
            let cv = m2.cell(c);
            let mut sv = [cv[0], cv[1], cv[2]];
            sv.sort_unstable();
            let (b0, b1, b2) = (sv[0], sv[1], sv[2]);
            let (t0, t1, t2) = (nv2 as u32 + sv[0], nv2 as u32 + sv[1], nv2 as u32 + sv[2]);
            for tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                cells3.extend_from_slice(&tet);
                clabs3.push(1);
            }
            facets3.extend_from_slice(&[b0, b1, b2, t0, t1, t2]);
            flabs3.extend_from_slice(&[99, 99]);
            for &(eu, ew) in &[(sv[0], sv[1]), (sv[0], sv[2]), (sv[1], sv[2])] {
                let p0 = m2.point(eu as usize);
                let p1 = m2.point(ew as usize);
                let xm = 0.5 * (p0[0] + p1[0]);
                if p0[1].abs() < 1e-10
                    && p1[1].abs() < 1e-10
                    && ((xm > -b_st && xm < -a_st) || (xm > a_st && xm < b_st))
                    && seen_strip_edges.insert((eu, ew))
                {
                    let lab = if xm < 0.0 { 1 } else { 2 };
                    let nu = nv2 as u32 + eu;
                    let nw = nv2 as u32 + ew;
                    facets3.extend_from_slice(&[eu, ew, nw, eu, nu, nw]);
                    flabs3.extend_from_slice(&[lab, lab]);
                }
            }
        }
        let ext2 = m2.exterior_facets();
        for (f, &is_ext) in ext2.iter().enumerate() {
            if is_ext {
                let fv = m2.facet(f);
                let (u, w) = (fv[0].min(fv[1]), fv[0].max(fv[1]));
                let nu = nv2 as u32 + u;
                let nw = nv2 as u32 + w;
                facets3.extend_from_slice(&[u, w, nw, u, nu, nw]);
                flabs3.extend_from_slice(&[99, 99]);
            }
        }
        let mesh_st3 = SimplexMesh::new(3, pts3, cells3, clabs3, facets3, flabs3);
        let mut p_st3 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
        p_st3.bcs = vec![
            Bc { labels: vec![1, 2], kind: BcKind::LondonSheet { l_k_sq: 0.0 } },
            Bc {
                labels: vec![1],
                kind: BcKind::CurrentTerminal { name: "s1".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
            },
            Bc {
                labels: vec![2],
                kind: BcKind::CurrentTerminal { name: "s2".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
            },
            Bc { labels: vec![99], kind: BcKind::Ground },
        ];
        let lmat_st3 =
            InductanceMatrixData::of(&solve_magnetostatics(&mesh_st3, &p_st3).unwrap()).unwrap().inductance_matrix_h;
        lmat_st3[0][0] + lmat_st3[1][1] - 2.0 * lmat_st3[0][1]
    };
    eprintln!(
        "parallel strips mixed_3d: L={l_diff_st3:.6e}, exact={l_line_exact:.6e}, rel_err={:.2}%",
        (l_diff_st3 - l_line_exact).abs() / l_line_exact * 100.0
    );
    assert!((l_diff_st3 - l_line_exact).abs() / l_line_exact < 0.03);
}

#[test]
fn heat_closed_forms_1d_3d_steady_convection_and_1d_transient_eigenmode() {
    // 1. 1D steady conduction along x with uniform source Q, Dirichlet at x=0, Robin convection at x=L:
    //    -k T'' = Q, T(0) = T0, -k T'(L) = h (T(L) - Text)
    let (l_slab, k_cond, q_src, t0, h_conv, t_ext) = (1.5_f64, 3.0_f64, 12.0_f64, 320.0_f64, 4.0_f64, 290.0_f64);
    let c1 =
        (h_conv * (t_ext - t0 + 0.5 * q_src * l_slab * l_slab / k_cond) + q_src * l_slab) / (k_cond + h_conv * l_slab);
    let t_exact_1d = |x: f64| -0.5 * (q_src / k_cond) * x * x + c1 * x + t0;

    let mesh2d = femcore::testkit::rect(l_slab, 0.3, 12, 3);
    let mut p_1d = Problem::new(Study::HeatStationary, 2);
    p_1d.materials.default.insert("k".into(), k_cond);
    p_1d.materials.default.insert("Q".into(), q_src);
    p_1d.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(t0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h: h_conv, text: t_ext } },
    ];
    let r_1d = solve_heat(&mesh2d, &p_1d).unwrap();
    let tv_1d = &r_1d.solutions[0].field("T").unwrap().values;
    let c_1d = mesh2d.dof_coordinates(2);
    let max_err_1d = (0..tv_1d.len()).map(|i| (tv_1d[i] - t_exact_1d(c_1d[2 * i])).abs()).fold(0.0, f64::max);
    assert!(max_err_1d < 1e-11, "1D P2 steady convection max_err {max_err_1d:e} >= 1e-11");

    // 2. 3D steady conduction along z on block [0, 1] x [0, 0.5] x [0, 0.8] with mixed Dirichlet/Flux/Convection
    let height = 0.8_f64;
    let mesh3d = femcore::testkit::block(1.0, 0.5, height, 4, 2, 6);
    let c1_3 =
        (h_conv * (t_ext - t0 + 0.5 * q_src * height * height / k_cond) + q_src * height) / (k_cond + h_conv * height);
    let t_exact_3d = |z: f64| -0.5 * (q_src / k_cond) * z * z + c1_3 * z + t0;
    let mut p_3d = Problem::new(Study::HeatStationary, 2);
    p_3d.materials.default.insert("k".into(), k_cond);
    p_3d.materials.default.insert("Q".into(), q_src);
    p_3d.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Temperature(t0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h: h_conv, text: t_ext } },
        Bc { labels: vec![1, 2, 5, 6], kind: BcKind::Flux(0.0) },
    ];
    let r_3d = solve_heat(&mesh3d, &p_3d).unwrap();
    let tv_3d = &r_3d.solutions[0].field("T").unwrap().values;
    let c_3d = mesh3d.dof_coordinates(2);
    let max_err_3d = (0..tv_3d.len()).map(|i| (tv_3d[i] - t_exact_3d(c_3d[3 * i + 2])).abs()).fold(0.0, f64::max);
    assert!(max_err_3d < 1e-11, "3D P2 steady convection max_err {max_err_3d:e} >= 1e-11");

    // 3. 1D transient eigenmode decay: T(x, t) = sin(pi x / L) exp(-alpha (pi / L)^2 t) (< 1e-3)
    let len = 2.0_f64;
    let (k_tr, rho_tr, cp_tr, t_end) = (1.8_f64, 1.5_f64, 0.8_f64, 0.5_f64);
    let alpha = k_tr / (rho_tr * cp_tr);
    let rod = femcore::testkit::rect(len, 0.2, 20, 2);
    let mut p_tr = Problem::new(Study::HeatTransient { t_end, n_out: 5, t0: "sin(pi * x / 2.0)".into() }, 2);
    p_tr.materials.default.insert("k".into(), k_tr);
    p_tr.materials.default.insert("rho".into(), rho_tr);
    p_tr.materials.default.insert("Cp".into(), cp_tr);
    p_tr.bcs = vec![Bc { labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    let r_tr = solve_heat(&rod, &p_tr).unwrap();
    let crd = rod.dof_coordinates(2);
    for (j, sol) in r_tr.solutions.iter().enumerate() {
        let t_now = t_end * (j as f64) / 5.0;
        let decay = (-alpha * (PI / len).powi(2) * t_now).exp();
        let vals = &sol.field("T").unwrap().values;
        let err = (0..vals.len()).map(|i| (vals[i] - (PI * crd[2 * i] / len).sin() * decay).abs()).fold(0.0, f64::max);
        assert!(err < 1e-3, "1D transient eigenmode step {j} err {err:e} >= 1e-3");
    }
}
