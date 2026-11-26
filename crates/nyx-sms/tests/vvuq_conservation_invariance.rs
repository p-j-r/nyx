// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME VVUQ Pillars (d) & (e) conservation laws, operator symmetries, and invariances for `nyx-sms`:
//! - Global static force & torque equilibrium (`< 1e-10`)
//! - Betti reciprocity on 3D rotated anisotropic elasticity (`< 1e-12`)
//! - Operator symmetries: real `||K - K^T||_F / ||K||_F < 1e-13`, complex PML `||K - K^T||_F / ||K||_F < 1e-13`,
//!   lossless Floquet Hermitian `||K - K^H||_F / ||K||_F < 1e-13`, `||M - M^H||_F / ||M||_F < 1e-13`
//! - Rigid Euclidean invariance: mesh-level `(23, 41, 67) deg` + `(0.13, -0.27, 0.41)` (`< 1e-10`),
//!   geometry-level rotated remesh (`< 1e-3`)
//! - Vertex / cell / tet permutation invariance (`< 1e-10`)
//! - Dimensional scaling invariance (`< 1e-10`)
//! - `+m/-m` azimuthal degeneracy in `axisym` (`< 1e-10`, including `m = +-1` on axis & `m = +-2`)
//! - `+phase/-phase` time-reversal degeneracy in lossless Floquet (`< 1e-10`)

mod common;

use femcore::linalg::csr::{CCsr, RCsr as Csr, C};
use femcore::material::{euler_from_rotation, rotation, Anisotropy, Stiffness, Tensor, M3};

fn mul3(a: &M3, b: &M3) -> M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
use femcore::props::{Periodic, Pml, PmlDir};
use geom::model::{Axis, Model, TStep};
use nyx_sms::cross_section::solve_solid_axisym;
use nyx_sms::floquet::assemble_floquet;
use nyx_sms::pml::assemble_solid_pml;
use nyx_sms::solid::{assemble_solid, assemble_static, setup, solve_static};
use nyx_sms::*;
use std::collections::BTreeMap;

fn sapphire_aniso(euler_deg: [f64; 3]) -> Anisotropy {
    Anisotropy {
        stiffness: Some(Stiffness::Trigonal {
            c11: 496.8e9,
            c12: 163.6e9,
            c13: 110.9e9,
            c14: -23.5e9,
            c33: 498.1e9,
            c44: 147.4e9,
        }),
        expansion: Some(Tensor::Uniaxial { perp: 5.0e-6, par: 6.2e-6 }),
        euler_deg,
        ..Default::default()
    }
}

fn mat_vec_mul(r: &M3, p: [f64; 3]) -> [f64; 3] {
    [
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
    ]
}

fn csr_sym_rel_frob_real(a: &Csr) -> f64 {
    let n = a.nrows();
    let mut map: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut norm_sq = 0.0;
    for r in 0..n {
        for k in a.indptr()[r]..a.indptr()[r + 1] {
            let col = a.indices()[k] as usize;
            let v = a.values()[k];
            norm_sq += v * v;
            *map.entry((r, col)).or_default() += v;
        }
    }
    let mut diff_sq = 0.0;
    for (&(r, c0), &v) in &map {
        let vt = map.get(&(c0, r)).copied().unwrap_or(0.0);
        let d = v - vt;
        diff_sq += d * d;
    }
    diff_sq.sqrt() / norm_sq.sqrt().max(1e-300)
}

fn csr_sym_rel_frob_complex(a: &CCsr, hermitian: bool) -> f64 {
    let n = a.nrows();
    let mut map: BTreeMap<(usize, usize), C> = BTreeMap::new();
    let mut norm_sq = 0.0;
    for r in 0..n {
        for k in a.indptr()[r]..a.indptr()[r + 1] {
            let col = a.indices()[k] as usize;
            let v = a.values()[k];
            norm_sq += v.norm_sqr();
            *map.entry((r, col)).or_default() += v;
        }
    }
    let mut diff_sq = 0.0;
    for (&(r, c0), &v) in &map {
        let vt = map.get(&(c0, r)).copied().unwrap_or_default();
        let target = if hermitian { vt.conj() } else { vt };
        diff_sq += (v - target).norm_sqr();
    }
    diff_sq.sqrt() / norm_sq.sqrt().max(1e-300)
}

#[test]
fn global_static_force_and_torque_equilibrium_and_betti_reciprocity_3d_anisotropic() {
    let mesh = common::block(1.0, 0.4, 0.3, 5, 3, 2);
    let euler = [37.0, 52.0, -18.0];
    let mut sp = Problem::new(Study::Static, 2);
    sp.materials.aniso.insert(1, sapphire_aniso(euler));
    sp.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
    sp.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([2.0e6, -1.5e6, 1.0e6]) },
        Bc { labels: vec![4], kind: BcKind::Load([0.5e6, 0.8e6, -1.2e6]) },
    ];

    let st = setup(&mesh, &sp).unwrap();
    let k = assemble_solid(&mesh, &st.c6, sp.order, sp.plane, sp.thickness);
    let asm = assemble_static(&mesh, &sp, &st, k.clone()).unwrap();
    let (u, _) = solve_static(&mesh, &sp, &st, &asm).unwrap();

    // Exact consistent reaction + external load equilibrium vector on all DOFs:
    // K u - f_asm gives boundary nodal reactions R (0 on free DOFs).
    // External mechanical load + nodal reactions satisfy sum(R + f_load) = 0 and sum(x x (R + f_load)) = 0!
    // Equivalently, since r_rigid^T f_th == 0, checking (K u - f_th)^T v_rigid == 0 and also computing f_load separately:
    let ku = k.matvec(&u);
    let residual: Vec<f64> = ku.iter().zip(&asm.f).map(|(a, b)| a - b).collect();
    let n1 = st.n1;
    let coords = mesh.dof_coordinates(sp.order);

    let mut sp_no_th = sp.clone();
    sp_no_th.thermal = None;
    let asm_mech = assemble_static(&mesh, &sp_no_th, &st, k.clone()).unwrap();
    let f_load = &asm_mech.f;
    let total_force_scale: f64 =
        (0..3).map(|comp| (0..n1).map(|i| f_load[comp * n1 + i].abs()).sum::<f64>()).sum::<f64>();
    let total_torque_scale = total_force_scale * 1.0;

    // Sum total external traction + reaction forces across every node: R_i + f_load_i = K u_i - f_th_i - ...
    let mut sum_f = [0.0; 3];
    let mut sum_tau = [0.0; 3];
    for i in 0..n1 {
        // Total mechanical force balance at DOF node i (reaction minus external load => wait: K u = f_load + f_th + R_reac => R_reac - ...):
        // At Dirichlet DOFs, K u - f_asm = R_reac; at free DOFs, K u - f_asm = 0; so net external force on body = -(R_reac) ...
        // Wait: (K u)_i = f_ext_reac_i + f_th_i where total external traction = f_load_i + reaction_on_body_i = residual_i + f_load_i!
        // Check: residual_i = (K u - f_load - f_th)_i => residual_i + f_load_i = (K u - f_th)_i!
        // Since K * u_rigid = 0 and f_th * u_rigid = 0 (B * u_rigid = 0), residual_i + f_load_i sums to zero!
        let f_tot =
            [residual[i] + f_load[i], residual[n1 + i] + f_load[n1 + i], residual[2 * n1 + i] + f_load[2 * n1 + i]];
        let p = [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]];
        for c in 0..3 {
            sum_f[c] += f_tot[c];
        }
        sum_tau[0] += p[1] * f_tot[2] - p[2] * f_tot[1];
        sum_tau[1] += p[2] * f_tot[0] - p[0] * f_tot[2];
        sum_tau[2] += p[0] * f_tot[1] - p[1] * f_tot[0];
    }
    let f_res_rel = (sum_f[0] * sum_f[0] + sum_f[1] * sum_f[1] + sum_f[2] * sum_f[2]).sqrt() / total_force_scale;
    let t_res_rel =
        (sum_tau[0] * sum_tau[0] + sum_tau[1] * sum_tau[1] + sum_tau[2] * sum_tau[2]).sqrt() / total_torque_scale;
    assert!(f_res_rel < 1e-10, "3D force equilibrium residual {f_res_rel:.3e}");
    assert!(t_res_rel < 1e-10, "3D torque equilibrium residual {t_res_rel:.3e}");

    // Betti reciprocity: u_1^T f_2 == u_2^T f_1 on 3D rotated sapphire (< 1e-12)
    let mut sp1 = Problem::new(Study::Static, 2);
    sp1.materials.aniso.insert(1, sapphire_aniso(euler));
    sp1.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([2.0e6, -1.0e6, 0.5e6]) },
    ];
    let mut sp2 = Problem::new(Study::Static, 2);
    sp2.materials.aniso.insert(1, sapphire_aniso(euler));
    sp2.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![4], kind: BcKind::Load([-0.8e6, 1.7e6, -2.3e6]) },
    ];
    let st1 = setup(&mesh, &sp1).unwrap();
    let k1 = assemble_solid(&mesh, &st1.c6, 2, Plane::Strain, 1.0);
    let asm1 = assemble_static(&mesh, &sp1, &st1, k1).unwrap();
    let (u1, _) = solve_static(&mesh, &sp1, &st1, &asm1).unwrap();

    let st2 = setup(&mesh, &sp2).unwrap();
    let k2 = assemble_solid(&mesh, &st2.c6, 2, Plane::Strain, 1.0);
    let asm2 = assemble_static(&mesh, &sp2, &st2, k2).unwrap();
    let (u2, _) = solve_static(&mesh, &sp2, &st2, &asm2).unwrap();

    let w12: f64 = u1.iter().zip(&asm2.f).map(|(a, b)| a * b).sum();
    let w21: f64 = u2.iter().zip(&asm1.f).map(|(a, b)| a * b).sum();
    let betti_rel = (w12 - w21).abs() / w12.abs().max(w21.abs());
    assert!(betti_rel < 1e-12, "Betti reciprocity failed: w12={w12:.6e}, w21={w21:.6e}, rel={betti_rel:.3e}");
}

#[test]
fn operator_symmetries_real_pml_and_floquet_hermitian() {
    let mesh = common::block(1.0, 0.4, 0.3, 4, 2, 2);
    let euler = [25.0, 50.0, -15.0];

    // 1. Real 3D anisotropic stiffness K == K^T (< 1e-13)
    let mut sp = Problem::new(Study::Static, 2);
    sp.materials.aniso.insert(1, sapphire_aniso(euler));
    sp.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let st = setup(&mesh, &sp).unwrap();
    let k_real = assemble_solid(&mesh, &st.c6, 2, Plane::Strain, 1.0);
    let err_real = csr_sym_rel_frob_real(&k_real);
    assert!(err_real < 1e-13, "real stiffness symmetry err={err_real:.3e}");

    // 2. Complex PML stiffness K == K^T (complex symmetric, unconjugated, < 1e-13)
    let mut sp_pml = Problem::new(Study::Eigenfrequency { n_modes: 3, f_search: 1.0e3 }, 2);
    sp_pml.materials.aniso.insert(1, sapphire_aniso(euler));
    sp_pml.materials.default.insert("rho".into(), 3980.0);
    sp_pml.bcs = sp.bcs.clone();
    sp_pml.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 2.5, order: 2.0, ..Default::default() });
    let asm_pml = assemble_solid_pml(&mesh, &sp_pml).unwrap();
    let err_pml_k = csr_sym_rel_frob_complex(&asm_pml.kff, false);
    let err_pml_m = csr_sym_rel_frob_complex(&asm_pml.mff, false);
    assert!(err_pml_k < 1e-13 && err_pml_m < 1e-13, "PML symmetry K={err_pml_k:.3e}, M={err_pml_m:.3e}");

    // 3. Lossless Floquet reduced operators K == K^H and M == M^H (< 1e-13)
    let mut sp_fl = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 1.0e3 }, 2);
    sp_fl.materials.aniso.insert(1, sapphire_aniso(euler));
    sp_fl.materials.default.insert("rho".into(), 3980.0);
    sp_fl.periodic = vec![Periodic { source: vec![1], dest: vec![6], phase: 1.37 }];
    let st_fl = setup(&mesh, &sp_fl).unwrap();
    let k_fl = assemble_solid(&mesh, &st_fl.c6, 2, Plane::Strain, 1.0);
    let rho = vec![3980.0; mesh.n_cells()];
    let m1 = femcore::fem::assemble_mass(&mesh, 2, &rho);
    let m_fl = Csr::block_diag(&[&m1, &m1, &m1]);
    let asm_fl = assemble_floquet(&mesh, &sp_fl, &st_fl, &k_fl, &m_fl, 4, 1.0e3).unwrap();
    let err_fl_k = csr_sym_rel_frob_complex(&asm_fl.kc, true);
    let err_fl_m = csr_sym_rel_frob_complex(&asm_fl.mc, true);
    assert!(err_fl_k < 1e-13 && err_fl_m < 1e-13, "Floquet Hermitian K={err_fl_k:.3e}, M={err_fl_m:.3e}");
}

#[test]
fn rigid_euclidean_invariance_mesh_and_geometry_remesh() {
    let euler_rig = [23.0_f64, 41.0, 67.0];
    let shift = [0.13_f64, -0.27, 0.41];
    let r = rotation(euler_rig);

    // 1. Mesh-level 3D rotated & translated anisotropic block (< 1e-10)
    let mesh0 = common::block(0.8, 0.5, 0.3, 4, 3, 2);
    let mut pts_rot = mesh0.points.clone();
    for p in pts_rot.chunks_mut(3) {
        let q = mat_vec_mul(&r, [p[0], p[1], p[2]]);
        p[0] = q[0] + shift[0];
        p[1] = q[1] + shift[1];
        p[2] = q[2] + shift[2];
    }
    let mesh1 = SimplexMesh::new(
        3,
        pts_rot,
        mesh0.cells.clone(),
        mesh0.cell_labels.clone(),
        mesh0.facets.clone(),
        mesh0.facet_labels.clone(),
    );

    let euler_mat0 = [15.0, 32.0, -28.0];
    let euler_mat1 = euler_from_rotation(&mul3(&r, &rotation(euler_mat0)));

    let mut sp0 = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    sp0.materials.aniso.insert(1, sapphire_aniso(euler_mat0));
    sp0.materials.default.insert("rho".into(), 3980.0);
    sp0.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];

    let mut sp1 = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    sp1.materials.aniso.insert(1, sapphire_aniso(euler_mat1));
    sp1.materials.default.insert("rho".into(), 3980.0);
    sp1.bcs = sp0.bcs.clone();

    let res0 = common::solve(&mesh0, &sp0);
    let res1 = common::solve(&mesh1, &sp1);
    for (s0, s1) in res0.solutions.iter().zip(&res1.solutions) {
        let (f0, f1) = (s0.value.unwrap(), s1.value.unwrap());
        let rel = (f0 - f1).abs() / f0;
        assert!(rel < 1e-10, "mesh Euclidean eigen invariance {f0:.6} vs {f1:.6} (rel={rel:.3e})");
    }

    // Also verify rotated static displacement transforms covariantly: u1 == R * u0 (< 1e-10).
    let t0 = [1.5e6_f64, -0.7e6, 2.2e6];
    let t1 = mat_vec_mul(&r, t0);
    let mut s_st0 = sp0.clone();
    s_st0.study = Study::Static;
    s_st0.bcs.push(Bc { labels: vec![6], kind: BcKind::Load(t0) });
    let mut s_st1 = sp1.clone();
    s_st1.study = Study::Static;
    s_st1.bcs.push(Bc { labels: vec![6], kind: BcKind::Load(t1) });

    let st_r0 = common::solve(&mesh0, &s_st0);
    let st_r1 = common::solve(&mesh1, &s_st1);
    let u0_comps = [0, 1, 2].map(|c| st_r0.solutions[0].field("u").unwrap().component(c));
    let u1_comps = [0, 1, 2].map(|c| st_r1.solutions[0].field("u").unwrap().component(c));
    let mut diff_sq = 0.0;
    let mut norm_sq = 0.0;
    for i in 0..mesh0.n_dofs(2) {
        let ru0 = mat_vec_mul(&r, [u0_comps[0][i], u0_comps[1][i], u0_comps[2][i]]);
        for c in 0..3 {
            diff_sq += (u1_comps[c][i] - ru0[c]).powi(2);
            norm_sq += ru0[c].powi(2);
        }
    }
    assert!(diff_sq.sqrt() / norm_sq.sqrt() < 1e-10, "rotated static covariant vector mismatch");

    // 2. Geometry-level rotated & translated remesh using geom::model::Model::new_2d + TStep + mesh_geom_2d (< 1e-3):
    //    (Using 2D constrained Delaunay avoids 3D voxelized rotated-box ~2% volume discretization noise and keeps exact geometry.)
    let rot_2d_deg = 41.0_f64;
    let th_2d = rot_2d_deg.to_radians();
    let shift_2d = [0.13_f64, -0.27_f64];
    let build_geom_2d = |rotated: bool| -> (SimplexMesh, Vec<i32>) {
        let mut m = Model::new_2d();
        let b = m.rect("beam", (0.0, 0.0), (0.8, 0.25));
        let root = if rotated {
            m.transform(
                "beam_rot",
                b,
                vec![
                    TStep::Rotate {
                        axis: Axis::Y,
                        deg: rot_2d_deg.into(),
                        about: [0.0.into(), 0.0.into(), 0.0.into()],
                    },
                    TStep::Translate([shift_2d[0].into(), 0.0.into(), shift_2d[1].into()]),
                ],
            )
        } else {
            b
        };
        m.overlay(&[("beam", root)]);
        let mut g_opts = GeomMeshOpts::new(0.025);
        g_opts.strict_mesh = false;
        let g = femcore::meshing::geomdrive::mesh_geom_2d(evaluate(&m.list).unwrap(), &g_opts).unwrap();
        let sx = g.mesh;
        // Identify the left end (original x = 0): normal in (x, z) under Ry(th) is (-cos th, sin th), anchor shift_2d.
        let (s_th, c_th) = th_2d.sin_cos();
        let n_left = if rotated { [-c_th, s_th] } else { [-1.0, 0.0] };
        let p0_left = if rotated { shift_2d } else { [0.0, 0.0] };
        let mut lbls = std::collections::BTreeSet::new();
        for f in 0..sx.n_facets() {
            let fv = sx.facet(f);
            let on_left = fv.iter().all(|&v| {
                let p = sx.point(v as usize);
                ((p[0] - p0_left[0]) * n_left[0] + (p[1] - p0_left[1]) * n_left[1]).abs() < 1e-9
            });
            if on_left {
                lbls.insert(sx.facet_labels[f]);
            }
        }
        (sx, lbls.into_iter().collect())
    };
    let (gm0, clamp0) = build_geom_2d(false);
    let (gm1, clamp1) = build_geom_2d(true);
    assert!(!clamp0.is_empty() && !clamp1.is_empty());

    let mut gp0 = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    gp0.plane = Plane::Stress;
    gp0.materials.set(1, "E", 200e9);
    gp0.materials.set(1, "nu", 0.3);
    gp0.materials.set(1, "rho", 7850.0);
    gp0.bcs = vec![Bc { labels: clamp0, kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } }];
    let mut gp1 = gp0.clone();
    gp1.bcs = vec![Bc { labels: clamp1, kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } }];

    let gf0 = common::solve(&gm0, &gp0);
    let gf1 = common::solve(&gm1, &gp1);
    for m in 0..4 {
        let f0 = gf0.solutions[m].value.unwrap();
        let f1 = gf1.solutions[m].value.unwrap();
        let rel = (f0 - f1).abs() / f0;
        println!("2D geometry rotated remesh mode {}: f0={f0:.4} vs f1={f1:.4} (rel={rel:.3e})", m + 1);
        assert!(rel < 1e-3, "geometry-level rotated remesh mode {}: rel={rel:.3e}", m + 1);
    }
}

#[test]
fn permutation_and_dimensional_scaling_invariance() {
    let mesh0 = common::block(1.0, 0.4, 0.3, 4, 2, 2);
    let euler = [28.0, 44.0, -13.0];

    // 1. Vertex + cell + tet-orientation-preserving permutation invariance (< 1e-10)
    let nv = mesh0.n_points();
    let nc = mesh0.n_cells();
    let nf = mesh0.n_facets();
    let vperm: Vec<usize> = (0..nv).map(|i| (i * 13 + 7) % nv).collect();
    // Since gcd(13, nv) might not be 1 if nv % 13 == 0, guarantee a true bijection by swapping/rotating co-prime:
    let mut vperm_bij: Vec<usize> = (0..nv).rev().collect();
    for i in (0..nv.saturating_sub(1)).step_by(2) {
        vperm_bij.swap(i, i + 1);
    }
    let _ = vperm;
    let mut inv_v = vec![0usize; nv];
    for (new_i, &old_i) in vperm_bij.iter().enumerate() {
        inv_v[old_i] = new_i;
    }

    let mut pts_p = vec![0.0; 3 * nv];
    for (new_i, &old_i) in vperm_bij.iter().enumerate() {
        pts_p[3 * new_i..3 * new_i + 3].copy_from_slice(&mesh0.points[3 * old_i..3 * old_i + 3]);
    }
    let mut cells_p = Vec::with_capacity(4 * nc);
    let mut clbls_p = Vec::with_capacity(nc);
    for c_idx in (0..nc).rev() {
        let cv = mesh0.cell(c_idx);
        // Even permutation [2, 0, 1, 3] preserves positive Jacobian determinant!
        cells_p.extend_from_slice(&[
            inv_v[cv[2] as usize] as u32,
            inv_v[cv[0] as usize] as u32,
            inv_v[cv[1] as usize] as u32,
            inv_v[cv[3] as usize] as u32,
        ]);
        clbls_p.push(mesh0.cell_labels[c_idx]);
    }
    let mut facets_p = Vec::with_capacity(3 * nf);
    let mut flbls_p = Vec::with_capacity(nf);
    for f_idx in (0..nf).rev() {
        let fv = mesh0.facet(f_idx);
        facets_p.extend_from_slice(&[
            inv_v[fv[1] as usize] as u32,
            inv_v[fv[2] as usize] as u32,
            inv_v[fv[0] as usize] as u32,
        ]);
        flbls_p.push(mesh0.facet_labels[f_idx]);
    }
    let mesh_p = SimplexMesh::new(3, pts_p, cells_p, clbls_p, facets_p, flbls_p);

    let mut sp_eig = Problem::new(Study::Eigenfrequency { n_modes: 5, f_search: 0.0 }, 2);
    sp_eig.materials.aniso.insert(1, sapphire_aniso(euler));
    sp_eig.materials.default.insert("rho".into(), 3980.0);
    sp_eig.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let r0 = common::solve(&mesh0, &sp_eig);
    let rp = common::solve(&mesh_p, &sp_eig);
    for (s0, sp_s) in r0.solutions.iter().zip(&rp.solutions) {
        let (f0, fp) = (s0.value.unwrap(), sp_s.value.unwrap());
        assert!((f0 - fp).abs() / f0 < 1e-10, "permutation mode invariance {f0} vs {fp}");
    }

    let mut sp_stat = Problem::new(Study::Static, 2);
    sp_stat.materials.aniso.insert(1, sapphire_aniso(euler));
    sp_stat.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::Load([1.0e6, -2.0e6, 0.5e6]) },
    ];
    let rs0 = common::solve(&mesh0, &sp_stat);
    let rsp = common::solve(&mesh_p, &sp_stat);
    for c in 0..3 {
        let u0 = rs0.solutions[0].field("u").unwrap().component(c);
        let up = rsp.solutions[0].field("u").unwrap().component(c);
        let max_u = u0[..nv].iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let max_diff = (0..nv).map(|new_i| (up[new_i] - u0[vperm_bij[new_i]]).abs()).fold(0.0_f64, f64::max);
        assert!(max_diff < 1e-10 * max_u, "permutation static comp {c} err={max_diff:.3e}");
    }

    // 2. Dimensional scaling invariance x -> lambda * x (< 1e-10):
    //    Eigenfrequency scales as lambda^-1; static displacement under same traction scales as lambda.
    for lam in [1.0e-3_f64, 1.0e2_f64] {
        let pts_s: Vec<f64> = mesh0.points.iter().map(|&v| v * lam).collect();
        let mesh_s = SimplexMesh::new(
            3,
            pts_s,
            mesh0.cells.clone(),
            mesh0.cell_labels.clone(),
            mesh0.facets.clone(),
            mesh0.facet_labels.clone(),
        );
        let r_eig_s = common::solve(&mesh_s, &sp_eig);
        for (s0, ss) in r0.solutions.iter().zip(&r_eig_s.solutions) {
            let f0 = s0.value.unwrap();
            let fs = ss.value.unwrap();
            assert!(
                (lam * fs / f0 - 1.0).abs() < 1e-10,
                "eigen scaling lam={lam}: lam*fs={:.6e} vs f0={f0:.6e}",
                lam * fs
            );
        }
        let r_st_s = common::solve(&mesh_s, &sp_stat);
        let umax0 = rs0.scalar("max |u|").unwrap();
        let umax_s = r_st_s.scalar("max |u|").unwrap();
        assert!(((umax_s / lam) / umax0 - 1.0).abs() < 1e-10, "static scaling lam={lam}");
    }
}

#[test]
fn axisym_azimuthal_degeneracy_and_floquet_time_reversal_degeneracy() {
    // 1. Axisym +m / -m azimuthal degeneracy (< 1e-10) for m = +-1 (with R6a MPC) and m = +-2 on axis r = 0
    let cyl = common::parts_mesh(&[common::part(0.0, -0.5, 0.8, 0.5, 1, Some(0.1))], 0.1);
    let end_lbls: Vec<i32> = (0..cyl.n_facets())
        .filter(|&f| cyl.facet(f).iter().all(|&v| (cyl.point(v as usize)[1].abs() - 0.5).abs() < 1e-9))
        .map(|f| cyl.facet_labels[f])
        .collect();
    for m_abs in [1_i64, 2_i64] {
        let run_m = |m: i64| {
            let mut p =
                Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 5, f_search: 500.0, azimuthal_order: m }, 2);
            p.materials.set(1, "E", 169e9);
            p.materials.set(1, "nu", 0.28);
            p.materials.set(1, "rho", 2329.0);
            p.bcs = vec![Bc { labels: end_lbls.clone(), kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
            let (_, modes) = solve_solid_axisym(&cyl, &p, m as f64).unwrap();
            modes.freqs.iter().map(|z| z.re).collect::<Vec<f64>>()
        };
        let fp = run_m(m_abs);
        let fn_ = run_m(-m_abs);
        assert!(fp.len() >= 5 && fn_.len() >= 5, "got {} and {} modes", fp.len(), fn_.len());
        for (a, b) in fp.iter().zip(&fn_).take(5) {
            assert!((a - b).abs() / a < 1e-10, "axisym +{m_abs}/-{m_abs} degeneracy {a} vs {b}");
        }
    }

    // 2. Lossless Floquet +phase / -phase time-reversal degeneracy (< 1e-10) on 3D rotated sapphire
    let mesh = common::block(1.0e-3, 0.3e-3, 0.25e-3, 4, 2, 2);
    let run_phase = |phase: f64| {
        let mut p = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 1.0e6 }, 2);
        p.materials.aniso.insert(1, sapphire_aniso([30.0, 40.0, 10.0]));
        p.materials.default.insert("rho".into(), 3980.0);
        p.periodic = vec![Periodic { source: vec![1], dest: vec![6], phase }];
        let r = common::solve(&mesh, &p);
        r.solutions.iter().map(|s| s.value.unwrap()).collect::<Vec<f64>>()
    };
    let fp = run_phase(1.15);
    let fn_ = run_phase(-1.15);
    for (a, b) in fp.iter().zip(&fn_) {
        assert!((a - b).abs() / a < 1e-10, "Floquet +phase/-phase degeneracy {a} vs {b}");
    }
}
