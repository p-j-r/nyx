// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! AF-3 Verification Suite: Open/Infinite Far-Field Boundary Conditions (`BcKind::FarField`).
//!
//! Verifies:
//! 1. Isolated conducting sphere (`R = 50 um`, `C_exact = 4 * pi * eps0 * R ~= 5.56325 fF`) inside a tight
//!    `L = 2R` box (`N = 1` terminal, no ground!): `< 0.25%` error for both Order 1 and Order 2 vs `> 22%`
//!    truncation error without `FarField`.
//! 2. Isolated circular disk (`R = 50 um`, `C_exact = 8 * eps0 * R ~= 3.54168 fF`) inside a tight `L = 2R`
//!    box (`N = 1` terminal, no ground!): `< 0.25%` error for both Order 1 and Order 2 vs `> 22%`
//!    truncation error with lateral Neumann walls.
//! 3. Multi-terminal open-domain capacitance matrix symmetry `|C_ij - C_ji| / C_ii < 1e-14`.
//! 4. Bit-identical results (`0.0` drift) between `NYX_THREADS=1` and `NYX_THREADS=4`.
//! 5. Clean `Error::Unsupported` rejection of `BcKind::FarField` in `planar_bem` and `NonFinite` check
//!    in `Problem::check_finite`.

use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Dfs, Problem, Study, EPS0};
use nyx_physics::{CapacitanceMatrixData, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

const TERM_LABEL: i32 = 1;
const OUTER_LABEL: i32 = 2;
const TOP_BOT_LABEL: i32 = 3;

/// Builds a 3D conforming mesh of the region between an inner conducting sphere of radius `r_sphere`
/// (boundary label `TERM_LABEL = 1`, with isoparametric P2 mid-edge nodes projected onto `r = r_sphere`)
/// and an outer cube `[-l_box, l_box]^3` (`l_box = 2 * r_sphere`, lateral walls label `OUTER_LABEL = 2`,
/// top/bottom `z = +/- l_box` walls label `TOP_BOT_LABEL = 3`).
fn build_sphere_in_box_mesh(r_sphere: f64, l_box: f64, n_side: usize, n_rad: usize) -> SimplexMesh {
    // 1. Enumerate unique surface nodes on the 6 faces of the integer cube [0..=n_side]^3
    let m = n_side;
    let mut grid_to_surf: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut surf_q: Vec<[f64; 3]> = Vec::new();
    let mut get_surf = |ix: usize, iy: usize, iz: usize| -> usize {
        if let Some(&id) = grid_to_surf.get(&(ix, iy, iz)) {
            return id;
        }
        let id = surf_q.len();
        let q = [
            2.0 * (ix as f64) / (m as f64) - 1.0,
            2.0 * (iy as f64) / (m as f64) - 1.0,
            2.0 * (iz as f64) / (m as f64) - 1.0,
        ];
        surf_q.push(q);
        grid_to_surf.insert((ix, iy, iz), id);
        id
    };

    // 6 cube faces with outward counter-clockwise quad orientation [s0, s1, s2, s3]:
    // (face_axis_is_z, [s0, s1, s2, s3])
    let mut quads: Vec<(bool, [usize; 4])> = Vec::with_capacity(6 * m * m);
    for i in 0..m {
        for j in 0..m {
            // +X (ix = m)
            quads.push((
                false,
                [get_surf(m, i, j), get_surf(m, i + 1, j), get_surf(m, i + 1, j + 1), get_surf(m, i, j + 1)],
            ));
            // -X (ix = 0)
            quads.push((
                false,
                [get_surf(0, j, i), get_surf(0, j + 1, i), get_surf(0, j + 1, i + 1), get_surf(0, j, i + 1)],
            ));
            // +Y (iy = m)
            quads.push((
                false,
                [get_surf(j, m, i), get_surf(j + 1, m, i), get_surf(j + 1, m, i + 1), get_surf(j, m, i + 1)],
            ));
            // -Y (iy = 0)
            quads.push((
                false,
                [get_surf(i, 0, j), get_surf(i + 1, 0, j), get_surf(i + 1, 0, j + 1), get_surf(i, 0, j + 1)],
            ));
            // +Z (iz = m)
            quads.push((
                true,
                [get_surf(i, j, m), get_surf(i + 1, j, m), get_surf(i + 1, j + 1, m), get_surf(i, j + 1, m)],
            ));
            // -Z (iz = 0)
            quads.push((
                true,
                [get_surf(j, i, 0), get_surf(j + 1, i, 0), get_surf(j + 1, i + 1, 0), get_surf(j, i + 1, 0)],
            ));
        }
    }

    let n_s = surf_q.len();
    let mut cm = CellMesh::<3>::new();
    cm.points.reserve((n_rad + 1) * n_s);
    for k in 0..=n_rad {
        let t = (k as f64 / n_rad as f64).powf(1.25);
        for q in &surf_q {
            let q_norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
            let d = [q[0] / q_norm, q[1] / q_norm, q[2] / q_norm];
            let p = [
                (1.0 - t) * r_sphere * d[0] + t * l_box * q[0],
                (1.0 - t) * r_sphere * d[1] + t * l_box * q[1],
                (1.0 - t) * r_sphere * d[2] + t * l_box * q[2],
            ];
            cm.points.push(p);
        }
    }

    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_rad {
        for &(_is_z, [s0, s1, s2, s3]) in &quads {
            let hex = [
                vid(k, s0),
                vid(k, s1),
                vid(k, s2),
                vid(k, s3),
                vid(k + 1, s0),
                vid(k + 1, s1),
                vid(k + 1, s2),
                vid(k + 1, s3),
            ];
            cm.push_cell(CellType::Hex, &hex, 1);
        }
    }
    for &(is_z, [s0, s1, s2, s3]) in &quads {
        cm.push_facet(&[vid(0, s0), vid(0, s1), vid(0, s2), vid(0, s3)], TERM_LABEL);
        let out_lab = if is_z { TOP_BOT_LABEL } else { OUTER_LABEL };
        cm.push_facet(&[vid(n_rad, s0), vid(n_rad, s1), vid(n_rad, s2), vid(n_rad, s3)], out_lab);
    }

    let (mut mesh, _) = SimplexMesh::from_cells_split(&cm).unwrap();

    // Project mid-edge geometry nodes of inner-sphere edges (both vertices in 0..n_s) onto r = r_sphere
    let mut edge_nodes = Vec::new();
    for &(u, v) in mesh.edges() {
        if (u as usize) < n_s && (v as usize) < n_s {
            let pu = mesh.point(u as usize);
            let pv = mesh.point(v as usize);
            let mx = 0.5 * (pu[0] + pv[0]);
            let my = 0.5 * (pu[1] + pv[1]);
            let mz = 0.5 * (pu[2] + pv[2]);
            let mr = (mx * mx + my * my + mz * mz).sqrt();
            if mr > 1.0e-18 {
                let sc = r_sphere / mr;
                edge_nodes.push(([u, v], [mx * sc, my * sc, mz * sc]));
            }
        }
    }
    mesh.set_edge_nodes(&edge_nodes);
    mesh
}

/// Builds a 3D mesh of an isolated circular disk of radius `r_disk` at `z = 0` inside a tight
/// box `[-l_box, l_box]^3` (`l_box = 2 * r_disk`) via `geom` + `sheet3d`, curving the circular disk
/// perimeter mid-edge nodes onto `x^2 + y^2 = r_disk^2`.
/// Returns `(mesh, disk_labels, lateral_outer_labels, top_bot_outer_labels)`.
fn build_disk_in_box_mesh(r_disk: f64, l_box: f64) -> (SimplexMesh, Vec<i32>, Vec<i32>, Vec<i32>) {
    let mut m = Model::new_3d();
    let vb = m.box_("vac_bot", (-l_box, -l_box, -l_box), (l_box, l_box, 0.0));
    let vt = m.box_("vac_top", (-l_box, -l_box, 0.0), (l_box, l_box, l_box));
    let dc = m.circle("disk_c", (0.0, 0.0), r_disk);
    m.sheet("disk", dc, 0.0);
    m.overlay(&[("vac_bot", vb), ("vac_top", vt)]);

    let h_global = 28.0e-6;
    let h_rim = 4.5e-6;
    let h_face = 9.0e-6;
    m.size("sz_rim", Selector::parse("disk:edge").unwrap(), Some(Expr::Num(h_rim)), None, None);
    m.size("sz_disk", Selector::Name("disk".into()), Some(Expr::Num(h_face)), None, None);

    let mut plan = MeshPlan::simple(h_global);
    plan.opts.max_dev = Some(r_disk * 2.5e-4);
    let mut g = mesh_model(&m.list, &plan).unwrap();
    let disk_labs = g.boundary("disk");
    let disk_set: BTreeSet<i32> = disk_labs.iter().copied().collect();

    // Project mid-edge nodes of edges whose endpoints lie on the disk perimeter (z == 0, rho == r_disk)
    // onto the exact circle rho = r_disk so there is zero polygon chord area loss.
    let mut rim_v = vec![false; g.mesh.n_points()];
    for (v, is_rim) in rim_v.iter_mut().enumerate() {
        let p = g.mesh.point(v);
        if p[2].abs() < 1.0e-12 && (p[0].hypot(p[1]) - r_disk).abs() < 1.0e-3 * r_disk {
            *is_rim = true;
        }
    }
    let mut edge_nodes = Vec::new();
    for &(u, v) in g.mesh.edges() {
        if rim_v[u as usize] && rim_v[v as usize] {
            let pu = g.mesh.point(u as usize);
            let pv = g.mesh.point(v as usize);
            let mx = 0.5 * (pu[0] + pv[0]);
            let my = 0.5 * (pu[1] + pv[1]);
            let mr = mx.hypot(my);
            if mr > 1.0e-18 {
                edge_nodes.push(([u, v], [mx * (r_disk / mr), my * (r_disk / mr), 0.0]));
            }
        }
    }
    g.mesh.set_edge_nodes(&edge_nodes);

    let ext = g.mesh.exterior_facets();
    let mut lat_labs = Vec::new();
    let mut tb_labs = Vec::new();
    let z_tol = 1.0e-9;
    for (f, &is_ext) in ext.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let lab = g.mesh.facet_labels[f];
        if disk_set.contains(&lab) {
            continue;
        }
        let fv = g.mesh.facet(f);
        let is_tb = fv.iter().all(|&v| (g.mesh.point(v as usize)[2].abs() - l_box).abs() < z_tol);
        if is_tb {
            if !tb_labs.contains(&lab) {
                tb_labs.push(lab);
            }
        } else if !lat_labs.contains(&lab) {
            lat_labs.push(lab);
        }
    }
    (g.mesh, disk_labs, lat_labs, tb_labs)
}

#[test]
fn isolated_sphere_and_disk_in_tight_2r_box_achieve_sub_025_percent_error_and_thread_determinism() {
    let r = 50.0e-6;
    let l_box = 2.0 * r; // tight L = 2R box ([-100 um, 100 um]^3)

    // Subprocess entry for NYX_THREADS=1 vs 4 bit-identity check
    if std::env::var("NYX_AF3_SUBPROC_CHECK").is_ok() {
        let mesh = build_sphere_in_box_mesh(r, l_box, 7, 5);
        let mut prob = Problem::new(Study::Electrostatic, 2);
        prob.materials.default.insert("eps_r".into(), 1.0);
        prob.bcs = vec![
            Bc { labels: vec![TERM_LABEL], kind: BcKind::Terminal("sphere".into()) },
            Bc {
                labels: vec![OUTER_LABEL, TOP_BOT_LABEL],
                kind: BcKind::FarField { order: 2, decay_power: 1.0, center: None },
            },
        ];
        let res = solve_electrostatics(&mesh, &prob).unwrap();
        let cap = CapacitanceMatrixData::of(&res).unwrap();
        let mut bits: Vec<u64> = vec![cap.maxwell_matrix_f[0][0].to_bits()];
        for &v in res.solutions[0].fields[0].values.iter().step_by(11) {
            bits.push(v.to_bits());
        }
        println!("NYX_AF3_BITS:{}", serde_json::to_string(&bits).unwrap());
        return;
    }

    // -------------------------------------------------------------------------
    // 1. Isolated Conducting Sphere (R = 50 um, C_exact = 4 * pi * eps0 * R ~= 5.56325 fF)
    // -------------------------------------------------------------------------
    let c_sphere_exact = 4.0 * PI * EPS0 * r;
    let sphere_mesh = build_sphere_in_box_mesh(r, l_box, 8, 6);

    // Baseline without FarField: lateral walls Neumann (dV/dn = 0), top/bottom grounded at z = +/- 2R
    let mut prob_sphere_neumann = Problem::new(Study::Electrostatic, 2);
    prob_sphere_neumann.materials.default.insert("eps_r".into(), 1.0);
    prob_sphere_neumann.bcs = vec![
        Bc { labels: vec![TERM_LABEL], kind: BcKind::Terminal("sphere".into()) },
        Bc { labels: vec![TOP_BOT_LABEL], kind: BcKind::Ground },
    ];
    let res_sphere_neumann = solve_electrostatics(&sphere_mesh, &prob_sphere_neumann).unwrap();
    let c_sphere_neumann = CapacitanceMatrixData::of(&res_sphere_neumann).unwrap().maxwell_matrix_f[0][0];
    let err_sphere_neumann = (c_sphere_neumann - c_sphere_exact).abs() / c_sphere_exact;

    // Order 1 FarField (N = 1, NO ground!)
    let mut prob_sphere_ff1 = Problem::new(Study::Electrostatic, 2);
    prob_sphere_ff1.materials.default.insert("eps_r".into(), 1.0);
    prob_sphere_ff1.bcs = vec![
        Bc { labels: vec![TERM_LABEL], kind: BcKind::Terminal("sphere".into()) },
        Bc {
            labels: vec![OUTER_LABEL, TOP_BOT_LABEL],
            kind: BcKind::FarField { order: 1, decay_power: 1.0, center: None },
        },
    ];
    let res_sphere_ff1 = solve_electrostatics(&sphere_mesh, &prob_sphere_ff1).unwrap();
    let c_sphere_ff1 = CapacitanceMatrixData::of(&res_sphere_ff1).unwrap().maxwell_matrix_f[0][0];
    let err_sphere_ff1 = (c_sphere_ff1 - c_sphere_exact).abs() / c_sphere_exact;

    // Order 2 FarField (N = 1, NO ground!)
    let mut prob_sphere_ff2 = Problem::new(Study::Electrostatic, 2);
    prob_sphere_ff2.materials.default.insert("eps_r".into(), 1.0);
    prob_sphere_ff2.bcs = vec![
        Bc { labels: vec![TERM_LABEL], kind: BcKind::Terminal("sphere".into()) },
        Bc {
            labels: vec![OUTER_LABEL, TOP_BOT_LABEL],
            kind: BcKind::FarField { order: 2, decay_power: 1.0, center: None },
        },
    ];
    let res_sphere_ff2 = solve_electrostatics(&sphere_mesh, &prob_sphere_ff2).unwrap();
    let c_sphere_ff2 = CapacitanceMatrixData::of(&res_sphere_ff2).unwrap().maxwell_matrix_f[0][0];
    let err_sphere_ff2 = (c_sphere_ff2 - c_sphere_exact).abs() / c_sphere_exact;

    eprintln!(
        "AF-3 Sphere (R=50 um, L=2R): C_exact={:.5} fF | Neumann C={:.5} fF (err {:.2}%) | FF1 C={:.5} fF (err {:.4}%) | FF2 C={:.5} fF (err {:.4}%)",
        c_sphere_exact * 1e15,
        c_sphere_neumann * 1e15,
        err_sphere_neumann * 100.0,
        c_sphere_ff1 * 1e15,
        err_sphere_ff1 * 100.0,
        c_sphere_ff2 * 1e15,
        err_sphere_ff2 * 100.0,
    );
    assert!(
        err_sphere_neumann > 0.22,
        "expected > 22% Neumann truncation error on tight L=2R sphere box, got {:.2}%",
        err_sphere_neumann * 100.0
    );
    assert!(
        err_sphere_ff1 < 0.0025,
        "Order 1 FarField sphere relative error {:.4}% exceeds 0.25% threshold",
        err_sphere_ff1 * 100.0
    );
    assert!(
        err_sphere_ff2 < 0.0025,
        "Order 2 FarField sphere relative error {:.4}% exceeds 0.25% threshold",
        err_sphere_ff2 * 100.0
    );

    // -------------------------------------------------------------------------
    // 2. Isolated Circular Disk (R = 50 um, C_exact = 8 * eps0 * R ~= 3.54168 fF)
    // -------------------------------------------------------------------------
    let c_disk_exact = 8.0 * EPS0 * r;
    let (disk_mesh, disk_labs, lat_labs, tb_labs) = build_disk_in_box_mesh(r, l_box);
    let mut all_outer_labs = lat_labs.clone();
    for &l in &tb_labs {
        if !all_outer_labs.contains(&l) {
            all_outer_labs.push(l);
        }
    }

    // Baseline without FarField: top/bottom Neumann (dV/dn = 0), lateral walls grounded at x,y = +/- 2R
    let mut prob_disk_neumann = Problem::new(Study::Electrostatic, 2).with_singular_edges("quarter_point");
    prob_disk_neumann.materials.default.insert("eps_r".into(), 1.0);
    prob_disk_neumann.bcs = vec![
        Bc { labels: disk_labs.clone(), kind: BcKind::Terminal("disk".into()) },
        Bc { labels: lat_labs.clone(), kind: BcKind::Ground },
    ];
    let res_disk_neumann = solve_electrostatics(&disk_mesh, &prob_disk_neumann).unwrap();
    let c_disk_neumann = CapacitanceMatrixData::of(&res_disk_neumann).unwrap().maxwell_matrix_f[0][0];
    let err_disk_neumann = (c_disk_neumann - c_disk_exact).abs() / c_disk_exact;

    // Order 1 FarField (N = 1, NO ground!)
    let mut prob_disk_ff1 = Problem::new(Study::Electrostatic, 2).with_singular_edges("quarter_point");
    prob_disk_ff1.materials.default.insert("eps_r".into(), 1.0);
    prob_disk_ff1.bcs = vec![
        Bc { labels: disk_labs.clone(), kind: BcKind::Terminal("disk".into()) },
        Bc { labels: all_outer_labs.clone(), kind: BcKind::FarField { order: 1, decay_power: 1.0, center: None } },
    ];
    let res_disk_ff1 = solve_electrostatics(&disk_mesh, &prob_disk_ff1).unwrap();
    let c_disk_ff1 = CapacitanceMatrixData::of(&res_disk_ff1).unwrap().maxwell_matrix_f[0][0];
    let err_disk_ff1 = (c_disk_ff1 - c_disk_exact).abs() / c_disk_exact;

    // Order 2 FarField (N = 1, NO ground!)
    let mut prob_disk_ff2 = Problem::new(Study::Electrostatic, 2).with_singular_edges("quarter_point");
    prob_disk_ff2.materials.default.insert("eps_r".into(), 1.0);
    prob_disk_ff2.bcs = vec![
        Bc { labels: disk_labs.clone(), kind: BcKind::Terminal("disk".into()) },
        Bc { labels: all_outer_labs.clone(), kind: BcKind::FarField { order: 2, decay_power: 1.0, center: None } },
    ];
    let res_disk_ff2 = solve_electrostatics(&disk_mesh, &prob_disk_ff2).unwrap();
    let c_disk_ff2 = CapacitanceMatrixData::of(&res_disk_ff2).unwrap().maxwell_matrix_f[0][0];
    let err_disk_ff2 = (c_disk_ff2 - c_disk_exact).abs() / c_disk_exact;

    eprintln!(
        "AF-3 Disk   (R=50 um, L=2R): C_exact={:.5} fF | Neumann C={:.5} fF (err {:.2}%) | FF1 C={:.5} fF (err {:.4}%) | FF2 C={:.5} fF (err {:.4}%)",
        c_disk_exact * 1e15,
        c_disk_neumann * 1e15,
        err_disk_neumann * 100.0,
        c_disk_ff1 * 1e15,
        err_disk_ff1 * 100.0,
        c_disk_ff2 * 1e15,
        err_disk_ff2 * 100.0,
    );
    assert!(
        err_disk_neumann > 0.22,
        "expected > 22% Neumann truncation error on tight L=2R disk box, got {:.2}%",
        err_disk_neumann * 100.0
    );
    assert!(
        err_disk_ff1 < 0.0025,
        "Order 1 FarField disk relative error {:.4}% exceeds 0.25% threshold",
        err_disk_ff1 * 100.0
    );
    assert!(
        err_disk_ff2 < 0.0025,
        "Order 2 FarField disk relative error {:.4}% exceeds 0.25% threshold",
        err_disk_ff2 * 100.0
    );

    // -------------------------------------------------------------------------
    // 3. Multi-terminal open-domain symmetry |C_ij - C_ji| / C_ii < 1e-14
    // -------------------------------------------------------------------------
    let mut prob_2t = Problem::new(Study::Electrostatic, 2);
    prob_2t.materials.default.insert("eps_r".into(), 1.0);
    prob_2t.bcs = vec![
        Bc { labels: vec![TERM_LABEL], kind: BcKind::Terminal("t1".into()) },
        Bc { labels: vec![TOP_BOT_LABEL], kind: BcKind::Terminal("t2".into()) },
        Bc {
            labels: vec![OUTER_LABEL],
            kind: BcKind::FarField { order: 2, decay_power: 1.0, center: Some([0.0, 0.0, 0.0]) },
        },
    ];
    let res_2t = solve_electrostatics(&sphere_mesh, &prob_2t).unwrap();
    let cap_2t = CapacitanceMatrixData::of(&res_2t).unwrap();
    let sym_err = (cap_2t.maxwell_matrix_f[0][1] - cap_2t.maxwell_matrix_f[1][0]).abs() / cap_2t.maxwell_matrix_f[0][0];
    assert!(sym_err < 1.0e-14, "symmetry error {sym_err:e} exceeds 1e-14");

    // -------------------------------------------------------------------------
    // 4. NYX_THREADS=1 vs NYX_THREADS=4 bit-identity
    // -------------------------------------------------------------------------
    let exe = std::env::current_exe().unwrap();
    let run_subproc = |threads: &str| -> Vec<u64> {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("isolated_sphere_and_disk_in_tight_2r_box_achieve_sub_025_percent_error_and_thread_determinism")
            .arg("--nocapture")
            .env("NYX_THREADS", threads)
            .env("NYX_AF3_SUBPROC_CHECK", "1")
            .output()
            .expect("failed to spawn child process");
        assert!(out.status.success(), "child failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let line = stdout.lines().find_map(|l| l.strip_prefix("NYX_AF3_BITS:")).expect("missing NYX_AF3_BITS");
        serde_json::from_str(line).unwrap()
    };
    let bits_t1 = run_subproc("1");
    let bits_t4 = run_subproc("4");
    assert_eq!(bits_t1, bits_t4, "FarField solution must be bit-identical between NYX_THREADS=1 and NYX_THREADS=4");
}

#[test]
fn farfield_validates_finite_and_planar_bem_rejects_with_unsupported() {
    let mut prob_nan = Problem::new(Study::Electrostatic, 1);
    prob_nan.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("t".into()) },
        Bc {
            labels: vec![2],
            kind: BcKind::FarField { order: 1, decay_power: 1.0, center: Some([0.0, f64::NAN, 0.0]) },
        },
    ];
    assert!(prob_nan.check_finite().is_err(), "check_finite must reject NaN in FarField center");

    // Envelope decode + planar_bem rejection with Error::Unsupported
    struct DummyResolver;
    impl LabelResolver for DummyResolver {
        fn domains(&self, _: &str) -> Result<Vec<i32>, Error> {
            Ok(vec![1])
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            match s {
                "t" => Ok(vec![1]),
                "outer" => Ok(vec![2]),
                _ => Err(Error::invalid("unknown")),
            }
        }
    }
    let mut spec = ProblemSpec::new("dfs", 1, "electrostatic");
    spec.study.params = json!({ "formulation": "planar_bem" });
    spec.materials.insert("vac".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    spec.bcs.push(nyx_physics::BcSpec { on: "t".into(), kind: "terminal".into(), params: json!({ "name": "t" }) });
    spec.bcs.push(nyx_physics::BcSpec {
        on: "outer".into(),
        kind: "far_field".into(),
        params: json!({ "order": 1, "decay_power": 1.0 }),
    });
    let mesh = Mesh::from_solver_mesh(&femcore::testkit::rect(1.0, 1.0, 2, 2)).unwrap();
    let err = Dfs.solve(&mesh, &spec, &mut Ctx::new(&DummyResolver)).unwrap_err();
    assert_eq!(err.code(), "unsupported", "expected Error::Unsupported from planar_bem + far_field, got {err:?}");
}
