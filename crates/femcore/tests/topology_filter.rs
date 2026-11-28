// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of `femcore::topology`: unstructured 2D/3D Helmholtz PDE density filter
//! partition-of-unity, volume conservation, exact discrete chain-rule pullback (`<g, F dx> == <F^T g, dx>`),
//! Heaviside + SIMP/RAMP/plasmon-safe optical permittivity derivatives vs central finite differences,
//! harmonic CAD `shape_velocity` and analytic simplex shape derivatives, and multi-thread determinism.

use femcore::fem::{assemble_laplace, assemble_mass};
use femcore::linalg::scalar::C64;
use femcore::meshing::geomdrive::{mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::geomdrive3::mesh_geom_3d;
use femcore::testkit::{block, rect};
use femcore::topology::{
    assemble_laplace_shape_deriv, assemble_mass_shape_deriv, compute_shape_velocity,
    continuation_schedule, displace_simplex_mesh, heaviside_project, optical_eps_interpolate,
    ramp_interpolate, simp_interpolate, HelmholtzFilter, OpticalInterpMode,
};
use geom::eval::Tier;
use geom::model::Model;

#[test]
fn helmholtz_filter_2d_and_3d_partition_of_unity_volume_and_exact_pullback() {
    let mesh2 = rect(1.2, 0.8, 11, 9);
    let mut mesh3 = block(1.0, 0.7, 0.5, 5, 4, 3);
    // Mark an interior box as design domain 2 to test subdomain restriction as well.
    mesh3
        .relabel(&[(vec![0.2, 0.15, 0.1, 0.8, 0.55, 0.4], 2)], &[])
        .unwrap();

    for (label, mesh, doms, lf) in [
        ("2D full", &mesh2, None, 0.15),
        ("3D full", &mesh3, None, 0.18),
        ("3D subdomain", &mesh3, Some(&[2i32][..]), 0.12),
    ] {
        let filter = HelmholtzFilter::new(mesh, lf, doms).unwrap();
        let nc = filter.n_design_cells();
        assert!(nc > 0, "{label}");

        // 1. Exact partition of unity: F(c * 1) == c * 1
        let ones = vec![0.73; nc];
        let filtered_ones = filter.filter_densities(&ones).unwrap();
        for (i, &v) in filtered_ones.iter().enumerate() {
            assert!(
                (v - 0.73).abs() < 1e-12,
                "{label} cell {i}: partition of unity error {:e}",
                (v - 0.73).abs()
            );
        }

        // 2. Pseudo-random oscillatory field: smoothing reduces peak-to-valley variation and conserves volume
        let rho: Vec<f64> = (0..nc)
            .map(|i| 0.5 + 0.45 * ((i as f64 * 1.73205).sin() * (i as f64 * 0.61803).cos()))
            .collect();
        let rho_tilde = filter.filter_densities(&rho).unwrap();
        let vol_before: f64 = rho
            .iter()
            .zip(filter.cell_volumes())
            .map(|(r, v)| r * v)
            .sum();
        let vol_after: f64 = rho_tilde
            .iter()
            .zip(filter.cell_volumes())
            .map(|(r, v)| r * v)
            .sum();
        assert!(
            (vol_before - vol_after).abs() <= 1e-12 * vol_before.abs(),
            "{label}: volume conservation violated: {vol_before} vs {vol_after}"
        );

        let spread_in = rho.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - rho.iter().copied().fold(f64::INFINITY, f64::min);
        let spread_out = rho_tilde.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - rho_tilde.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(spread_out < spread_in, "{label}: Helmholtz filter must smooth spatial oscillations");

        // 3. Exact discrete chain-rule transpose identity: <g, F drho> == <F^T g, drho>
        let g: Vec<f64> = (0..nc).map(|i| ((i as f64 + 0.3) * 2.17).cos()).collect();
        let drho: Vec<f64> = (0..nc).map(|i| ((i as f64 + 1.1) * 1.31).sin()).collect();
        let f_drho = filter.filter_densities(&drho).unwrap();
        let ft_g = filter.pullback_gradient(&g).unwrap();
        let lhs: f64 = g.iter().zip(&f_drho).map(|(a, b)| a * b).sum();
        let rhs: f64 = ft_g.iter().zip(&drho).map(|(a, b)| a * b).sum();
        let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1e-30);
        assert!(rel < 1e-13, "{label}: adjoint transpose mismatch rel = {rel:e}");

        // 4. Volume-weighted self-adjointness: <u, F v>_V == <F u, v>_V
        let f_rho = &rho_tilde;
        let v_lhs: f64 = rho
            .iter()
            .zip(&f_drho)
            .zip(filter.cell_volumes())
            .map(|((u, fv), vol)| vol * u * fv)
            .sum();
        let v_rhs: f64 = f_rho
            .iter()
            .zip(&drho)
            .zip(filter.cell_volumes())
            .map(|((fu, v), vol)| vol * fu * v)
            .sum();
        let v_rel = (v_lhs - v_rhs).abs() / v_lhs.abs().max(v_rhs.abs()).max(1e-30);
        assert!(v_rel < 1e-13, "{label}: V-inner-product self-adjointness rel = {v_rel:e}");

        // 5. End-to-end chain rule with nonlinear objective J(rho) = sum_e V_e * simp(heaviside(F(rho)_e))^2
        let (beta, eta) = (8.0, 0.5);
        let eval_j = |x: &[f64]| -> f64 {
            let rt = filter.filter_densities(x).unwrap();
            let (rb, _) = heaviside_project(&rt, beta, eta);
            rb.iter()
                .zip(filter.cell_volumes())
                .map(|(&r, &vol)| {
                    let (val, _) = simp_interpolate(r, 1.0, 12.0, 3.0);
                    vol * val * val
                })
                .sum()
        };
        let (rb, drb_drt) = heaviside_project(&rho_tilde, beta, eta);
        let dj_drt: Vec<f64> = rb
            .iter()
            .zip(&drb_drt)
            .zip(filter.cell_volumes())
            .map(|((&r, &dr), &vol)| {
                let (val, dval) = simp_interpolate(r, 1.0, 12.0, 3.0);
                2.0 * vol * val * dval * dr
            })
            .collect();
        let grad_adj = filter.pullback_gradient(&dj_drt).unwrap();
        let h = 1e-6;
        let dir_pos: Vec<f64> = (0..nc).map(|i| 0.5 + 0.5 * ((i as f64 + 0.7) * 1.13).cos().abs()).collect();
        let rho_p: Vec<f64> = rho.iter().zip(&dir_pos).map(|(x, d)| x + h * d).collect();
        let rho_m: Vec<f64> = rho.iter().zip(&dir_pos).map(|(x, d)| x - h * d).collect();
        let fd_dir = (eval_j(&rho_p) - eval_j(&rho_m)) / (2.0 * h);
        let adj_dir: f64 = grad_adj.iter().zip(&dir_pos).map(|(g_i, d_i)| g_i * d_i).sum();
        let chain_rel = (fd_dir - adj_dir).abs() / fd_dir.abs().max(1e-30);
        assert!(
            chain_rel < 1e-9,
            "{label}: chain-rule FD vs pullback mismatch rel = {chain_rel:e}"
        );
    }
}

#[test]
fn interpolations_and_plasmon_safe_optical_permittivity_match_finite_differences() {
    let h = 1e-7;
    for rho in [0.08, 0.25, 0.5, 0.77, 0.94] {
        let (_, ds) = simp_interpolate(rho, 1e-3, 200.0, 3.0);
        let fds = (simp_interpolate(rho + h, 1e-3, 200.0, 3.0).0
            - simp_interpolate(rho - h, 1e-3, 200.0, 3.0).0)
            / (2.0 * h);
        assert!((ds - fds).abs() / fds.abs() < 1e-9);

        let (_, dr) = ramp_interpolate(rho, 1.0, 50.0, 4.0);
        let fdr = (ramp_interpolate(rho + h, 1.0, 50.0, 4.0).0
            - ramp_interpolate(rho - h, 1.0, 50.0, 4.0).0)
            / (2.0 * h);
        assert!((dr - fdr).abs() / fdr.abs() < 1e-9);

        // Air-to-metal interpolation: verify zero plasmon crossing (|eps(rho)| bounded well above 0)
        // and exact analytic complex derivative!
        let eps_air = C64::new(1.0, 0.0);
        let eps_metal = C64::new(-12.0, -0.8);
        let mode = OpticalInterpMode::Auto { penalty: 1.0, damping_floor: 0.2 };
        let (eps0, deps) = optical_eps_interpolate(rho, eps_air, eps_metal, mode);
        assert!(
            (eps0.re * eps0.re + eps0.im * eps0.im).sqrt() > 0.5,
            "plasmon-safe interpolation must avoid |eps| -> 0 crossing, got {eps0:?} at rho = {rho}"
        );
        let (eps_p, _) = optical_eps_interpolate(rho + h, eps_air, eps_metal, mode);
        let (eps_m, _) = optical_eps_interpolate(rho - h, eps_air, eps_metal, mode);
        let fd_re = (eps_p.re - eps_m.re) / (2.0 * h);
        let fd_im = (eps_p.im - eps_m.im) / (2.0 * h);
        assert!((deps.re - fd_re).abs() / fd_re.abs().max(1.0) < 1e-9);
        assert!((deps.im - fd_im).abs() / fd_im.abs().max(1.0) < 1e-9);
    }

    let sched = continuation_schedule(1.0, 32.0, 6, 15);
    assert_eq!(sched.len(), 6);
    assert!((sched[0].0 - 1.0).abs() < 1e-12);
    assert!((sched[5].0 - 32.0).abs() < 1e-12);
}

#[test]
fn cad_harmonic_shape_velocity_and_simplex_operator_derivatives_are_exact() {
    // 2D parametric strip waveguide in cladding
    let mut m = Model::new_2d();
    m.param("w", 0.5);
    let clad = m.rect("clad", (-1.5, -1.0), (1.5, 1.0));
    let core = m.rect("core", ("-w/2", -0.2), ("w/2", 0.2));
    m.overlay(&[("clad", clad), ("core", core)]);
    let eval0 = m.eval(Tier::Complex).unwrap();
    let meshed = mesh_geom_2d(eval0, &GeomMeshOpts::new(0.25)).unwrap();
    let mesh = meshed.into_exchange_mesh().unwrap();

    let vel = compute_shape_velocity(&m.list, &meshed.eval, &mesh, "w", None).unwrap();
    let simplex = mesh.discretise().unwrap();
    assert_eq!(vel.vertex_vel.len(), simplex.n_points());

    // Vertices at x = +-0.25, |z| <= 0.2 must have vx = +-0.5, while outer box x = +-1.5 has vx = 0.
    let mut checked_core_wall = 0usize;
    for i in 0..simplex.n_points() {
        let p = simplex.point(i);
        if (p[0].abs() - 0.25).abs() < 1e-9 && p[1].abs() <= 0.2 + 1e-9 {
            let expected_vx = 0.5 * p[0].signum();
            assert!(
                (vel.vertex_vel[i][0] - expected_vx).abs() < 1e-8,
                "core wall vertex {i} at {p:?}: vx = {}, expected {expected_vx}",
                vel.vertex_vel[i][0]
            );
            checked_core_wall += 1;
        }
        if (p[0].abs() - 1.5).abs() < 1e-9 {
            assert!(vel.vertex_vel[i][0].abs() < 1e-12);
        }
    }
    assert!(checked_core_wall >= 4);

    // Analytic `assemble_laplace_shape_deriv` and `assemble_mass_shape_deriv` vs fixed-mesh FD (order 2)
    let coeff: Vec<f64> = simplex
        .cell_labels
        .iter()
        .map(|&l| if l == 2 { 12.0 } else { 2.1 })
        .collect();
    let dk = assemble_laplace_shape_deriv(&simplex, 2, &coeff, &vel).unwrap();
    let dm = assemble_mass_shape_deriv(&simplex, 2, &coeff, &vel).unwrap();

    let h = 1e-6;
    let sp = displace_simplex_mesh(&simplex, h, &vel).unwrap();
    let sm = displace_simplex_mesh(&simplex, -h, &vel).unwrap();
    let kp = assemble_laplace(&sp, 2, &coeff);
    let km = assemble_laplace(&sm, 2, &coeff);
    let mp = assemble_mass(&sp, 2, &coeff);
    let mm = assemble_mass(&sm, 2, &coeff);

    let nd = simplex.n_dofs(2);
    let u: Vec<f64> = (0..nd).map(|i| ((i as f64 + 1.0) * 0.73).sin()).collect();
    let v: Vec<f64> = (0..nd).map(|i| ((i as f64 + 2.0) * 1.19).cos()).collect();
    let quad_form = |mat: &femcore::linalg::csr::RCsr| -> f64 {
        let mu = mat.matvec(&u);
        v.iter().zip(&mu).map(|(vi, mui)| vi * mui).sum()
    };

    let dk_an = quad_form(&dk);
    let dk_fd = (quad_form(&kp) - quad_form(&km)) / (2.0 * h);
    assert!(
        (dk_an - dk_fd).abs() / dk_fd.abs() < 5e-9,
        "2D P2 Laplace shape deriv mismatch: an={dk_an} fd={dk_fd}"
    );

    let dm_an = quad_form(&dm);
    let dm_fd = (quad_form(&mp) - quad_form(&mm)) / (2.0 * h);
    assert!(
        (dm_an - dm_fd).abs() / dm_fd.abs() < 5e-9,
        "2D P2 Mass shape deriv mismatch: an={dm_an} fd={dm_fd}"
    );

    // 3D curved cylinder: radial expansion `d(Vol)/dR == 2 * pi * R * H` (on curved quadratic mesh)
    let mut m3 = Model::new_3d();
    m3.param("R", 0.4);
    let circ = m3.circle("c", (0.0, 0.0), "R");
    let cyl = m3.extrude("cyl", circ, 0.0, 0.6);
    m3.overlay(&[("cav", cyl)]);
    let ev3 = m3.eval(Tier::Complex).unwrap();
    let mut opts3 = GeomMeshOpts::new(0.22);
    opts3.curved = true;
    let meshed3 = mesh_geom_3d(ev3, &opts3).unwrap();
    let ex3 = meshed3.into_exchange_mesh().unwrap();
    let vel3 = compute_shape_velocity(&m3.list, &meshed3.eval, &ex3, "R", None).unwrap();
    assert!(!vel3.edge_mid_vel.is_empty(), "curved 3D cylinder must populate edge_mid_vel");
    let s3 = ex3.discretise().unwrap();
    let ones3 = vec![1.0; s3.n_cells()];
    let dk3 = assemble_laplace_shape_deriv(&s3, 2, &ones3, &vel3).unwrap();
    let dm3 = assemble_mass_shape_deriv(&s3, 1, &ones3, &vel3).unwrap();
    let d_vol_total: f64 = dm3.matvec(&vec![1.0; s3.n_dofs(1)]).iter().sum();

    // Exactness check: 3D analytic operator shape derivatives vs fixed-topology displaced mesh FD
    let s3p = displace_simplex_mesh(&s3, h, &vel3).unwrap();
    let s3m = displace_simplex_mesh(&s3, -h, &vel3).unwrap();
    let vol_fd = (femcore::fem::geometry(&s3p).vol.iter().sum::<f64>()
        - femcore::fem::geometry(&s3m).vol.iter().sum::<f64>())
        / (2.0 * h);
    assert!(
        (d_vol_total - vol_fd).abs() / vol_fd.abs() < 1e-9,
        "3D mass shape derivative must match displaced-mesh FD to < 1e-9: an={d_vol_total} fd={vol_fd}"
    );
    let nd3 = s3.n_dofs(2);
    let u3: Vec<f64> = (0..nd3).map(|i| ((i as f64 + 1.0) * 0.61).sin()).collect();
    let k3p = assemble_laplace(&s3p, 2, &ones3);
    let k3m = assemble_laplace(&s3m, 2, &ones3);
    let q3 = |mat: &femcore::linalg::csr::RCsr| -> f64 {
        u3.iter().zip(mat.matvec(&u3)).map(|(a, b)| a * b).sum()
    };
    let dk3_an = q3(&dk3);
    let dk3_fd = (q3(&k3p) - q3(&k3m)) / (2.0 * h);
    assert!(
        (dk3_an - dk3_fd).abs() / dk3_fd.abs() < 1e-8,
        "3D P2 Laplace shape deriv mismatch: an={dk3_an} fd={dk3_fd}"
    );

    // Physical CAD scaling accuracy on cut-cell 3D cylinder mesh (bevel creases have O(alpha*h^2) cut error)
    let v_simplex: f64 = femcore::fem::geometry(&s3).vol.iter().sum();
    let rel_vol = (d_vol_total - 2.0 * v_simplex / 0.4).abs() / (2.0 * v_simplex / 0.4);
    assert!(
        rel_vol < 1e-3,
        "3D cylinder harmonic velocity must match radial cylinder scaling within cut-cell bevel accuracy (rel={rel_vol:e})"
    );
}

#[test]
fn helmholtz_filter_is_deterministic_across_thread_counts() {
    if std::env::var("NYX_TOPO_SUBPROC").is_ok() {
        let mesh3 = block(1.0, 0.8, 0.6, 6, 5, 4);
        let filter = HelmholtzFilter::new(&mesh3, 0.14, None).unwrap();
        let nc = filter.n_design_cells();
        let rho: Vec<f64> = (0..nc).map(|i| 0.5 + 0.4 * ((i as f64 + 0.5) * 1.37).sin()).collect();
        let rt = filter.filter_densities(&rho).unwrap();
        let pb = filter.pullback_gradient(&rt).unwrap();
        let bits: Vec<String> = pb.iter().take(16).map(|x| format!("{:016x}", x.to_bits())).collect();
        println!("TOPO_BITS={}", bits.join(","));
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_with = |nt: &str| -> String {
        let out = std::process::Command::new(&exe)
            .args(["--exact", "helmholtz_filter_is_deterministic_across_thread_counts", "--nocapture"])
            .env("NYX_TOPO_SUBPROC", "1")
            .env("NYX_THREADS", nt)
            .output()
            .unwrap();
        assert!(out.status.success());
        String::from_utf8(out.stdout)
            .unwrap()
            .lines()
            .find(|l| l.starts_with("TOPO_BITS="))
            .unwrap()
            .to_string()
    };
    let b1 = run_with("1");
    let b8 = run_with("8");
    assert_eq!(b1, b8, "HelmholtzFilter pullback must be bit-identical across NYX_THREADS=1 vs 8");
}
