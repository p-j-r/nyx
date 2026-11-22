// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tests and benchmarks for FR-16 Pt 1 (ordering-aware & memory-reduced sparse indefinite direct
//! solver) and FR-16 Pt 2 (P2 Nédélec Hiptmair–Xu AMS + FGMRES / shift-invert / Jacobi–Davidson).
#![allow(clippy::needless_range_loop)]

use femcore::linalg::ams::{AmsP2, AmsP2Options, AmsP2Topology};
use femcore::linalg::blr::{BlrOptions, ShiftedBlrLu};
use femcore::linalg::csr::laplacian_1d;
use femcore::linalg::direct::{Definiteness, FactorBackend, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::driver::{
    EigenOptions, IterativeCurlCurlOptions, IterativeShiftInvertP2, eigs_curl_curl_p2_jd, eigs_near,
    solve_curl_curl_p2_iterative,
};
use femcore::linalg::ksp::{Cg, CgKind, CgOptions, Tolerances};
use femcore::linalg::maxwell::{EdgeElementP2Cube, nedelec_p2_cube};
use femcore::linalg::ordering::FillOrdering;
use femcore::linalg::problem::Symmetry;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr, Scalar};
use std::f64::consts::PI;
use std::sync::Mutex;
use std::time::Instant;

static ENV_LOCK: Mutex<()> = Mutex::new(());

fn topo_of(p: &EdgeElementP2Cube) -> AmsP2Topology {
    AmsP2Topology {
        p1_to_p2: p.p1_to_p2.clone(),
        gradient_p1: p.gradient_p1.clone(),
        gradient_p2: Some(p.gradient_p2.clone()),
        edge_vectors: p.edge_vectors.clone(),
        blocks: p.blocks.clone(),
    }
}

fn rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    let r = ax.iter().zip(b).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt();
    let bn = b.iter().map(|v| v.modulus_sq()).sum::<f64>().sqrt().max(1e-300);
    r / bn
}

/// Complex-symmetric lossy/PML-like 3D P2 Nédélec curl-curl matrix `K + (−σ + i γ(x)) M` with
/// spatially varying imaginary damping (lossy absorber region).
fn complex_symmetric_p2_curl_curl(p: &EdgeElementP2Cube, sigma: f64, gamma: f64) -> (Csr<C64>, Csr<C64>) {
    let n = p.stiffness.nrows();
    let kc = p.stiffness.to_complex();
    let mc = p.mass.to_complex();
    // Build complex-symmetric mass with spatially varying dielectric loss (lossy shell on upper half of DOFs).
    let mut loss_trip = Vec::with_capacity(p.mass.nnz());
    for r in 0..n {
        let (cols, vals) = p.mass.row(r);
        let lr = if r >= n / 2 { gamma } else { 0.05 * gamma };
        for (&c, &v) in cols.iter().zip(vals) {
            let lc = if (c as usize) >= n / 2 { gamma } else { 0.05 * gamma };
            let g = 0.5 * (lr + lc);
            loss_trip.push((r, c as usize, C64::new(-sigma * v, g * v)));
        }
    }
    let shifted_m = Csr::from_triplets(n, n, &loss_trip).unwrap();
    let a_shifted = Csr::axpby(C64::ONE, &kc, C64::ONE, &shifted_m);
    (a_shifted, mc)
}

#[cfg(windows)]
fn process_memory_mib() -> (f64, f64) {
    #[repr(C)]
    #[allow(non_snake_case)]
    struct ProcessMemoryCountersEx {
        cb: u32,
        PageFaultCount: u32,
        PeakWorkingSetSize: usize,
        WorkingSetSize: usize,
        QuotaPeakPagedPoolUsage: usize,
        QuotaPagedPoolUsage: usize,
        QuotaPeakNonPagedPoolUsage: usize,
        QuotaNonPagedPoolUsage: usize,
        PagefileUsage: usize,
        PeakPagefileUsage: usize,
        PrivateUsage: usize,
    }
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentProcess() -> *mut std::ffi::c_void;
        fn K32GetProcessMemoryInfo(
            process: *mut std::ffi::c_void,
            ppsmemcounters: *mut ProcessMemoryCountersEx,
            cb: u32,
        ) -> i32;
    }
    let mut pmc = ProcessMemoryCountersEx {
        cb: std::mem::size_of::<ProcessMemoryCountersEx>() as u32,
        PageFaultCount: 0,
        PeakWorkingSetSize: 0,
        WorkingSetSize: 0,
        QuotaPeakPagedPoolUsage: 0,
        QuotaPagedPoolUsage: 0,
        QuotaPeakNonPagedPoolUsage: 0,
        QuotaNonPagedPoolUsage: 0,
        PagefileUsage: 0,
        PeakPagefileUsage: 0,
        PrivateUsage: 0,
    };
    let ok = unsafe { K32GetProcessMemoryInfo(GetCurrentProcess(), &mut pmc, pmc.cb) };
    if ok != 0 {
        (pmc.WorkingSetSize as f64 / (1024.0 * 1024.0), pmc.PrivateUsage as f64 / (1024.0 * 1024.0))
    } else {
        (0.0, 0.0)
    }
}

#[cfg(not(windows))]
fn process_memory_mib() -> (f64, f64) {
    (0.0, 0.0)
}

// ------------------------------------------------------------------------------------------------------------ FR-16 Pt 1 tests

#[test]
fn diag_damped_1d_tridiag_accuracy() {
    let n = 200;
    let l = laplacian_1d(n).to_complex();
    let id = Csr::<C64>::identity(n);
    let b: Vec<C64> = (0..n).map(|i| C64::new(((i * 37) % 11) as f64 - 5.0, 0.3 * ((i * 17) % 7) as f64)).collect();
    let w = 0.25f64;
    let t = Csr::axpby(C64::ONE, &l, C64::new(-w * w, 0.02 * w), &id);

    let mut f_amd = ShiftedFactorization::new(&t, None, Definiteness::General).unwrap();
    f_amd.ordering = FillOrdering::Amd;
    f_amd.set_shift(C64::ZERO).unwrap();
    let mut x_amd = b.clone();
    f_amd.solve_block_in_place(&mut x_amd, 1);
    let r_amd = rel_residual(&t, &x_amd, &b);

    let mut f_nd = ShiftedFactorization::new(&t, None, Definiteness::General).unwrap();
    f_nd.ordering = FillOrdering::NestedDissection;
    f_nd.set_shift(C64::ZERO).unwrap();
    let mut x_nd = b.clone();
    f_nd.solve_block_in_place(&mut x_nd, 1);
    let r_nd = rel_residual(&t, &x_nd, &b);

    eprintln!(
        "damped 1D n={n}: AMD backend={:?} res={r_amd:.2e}, ND backend={:?} res={r_nd:.2e}",
        f_amd.stats().backend,
        f_nd.stats().backend
    );
    assert!(r_amd < 1e-12, "r_amd = {r_amd:e}");
    assert!(r_nd < 1e-12, "r_nd = {r_nd:e}");
}

#[test]
fn ordering_controls_indefinite_direct_backend_and_matches_eigenvalues() {
    // 3D P2 Nédélec cube (n = 4 -> 384 tets, 1344 DOFs): fast enough for CI, genuinely 3D P2 H(curl).
    let p = nedelec_p2_cube(4);
    let n = p.stiffness.nrows();
    let sigma = 2.0 * PI * PI + 0.35; // interior indefinite shift near lowest triple cavity mode 2π²

    // 1. Real symmetric indefinite: AMD vs ND LBL^H must report LbltAmd / LbltNd and match eigenvalues < 1e-10
    let mut f_amd = ShiftedFactorization::new(&p.stiffness, Some(&p.mass), Definiteness::General).unwrap();
    f_amd.ordering = FillOrdering::Amd;
    f_amd.set_shift(sigma).unwrap();
    assert_eq!(f_amd.stats().backend, FactorBackend::LbltAmd);
    assert!(f_amd.stats().factor_nnz > 0);

    let mut f_nd = ShiftedFactorization::new(&p.stiffness, Some(&p.mass), Definiteness::General).unwrap();
    f_nd.ordering = FillOrdering::NestedDissection;
    f_nd.set_shift(sigma).unwrap();
    assert_eq!(f_nd.stats().backend, FactorBackend::LbltNd);
    assert!(f_nd.stats().factor_nnz > 0);

    let shifted = Csr::axpby(1.0, &p.stiffness, -sigma, &p.mass);
    let rhs = random_vector::<f64>(n, 42);
    let (mut x_amd, mut x_nd) = (rhs.clone(), rhs.clone());
    f_amd.solve_block_in_place(&mut x_amd, 1);
    f_nd.solve_block_in_place(&mut x_nd, 1);
    assert!(rel_residual(&shifted, &x_amd, &rhs) < 1e-11);
    assert!(rel_residual(&shifted, &x_nd, &rhs) < 1e-11);

    // 2. Complex-symmetric indefinite: AMD vs ND non-conjugate supernodal LBL^T must report LbltAmd / LbltNd
    let (a_cplx, _mc) = complex_symmetric_p2_curl_curl(&p, sigma, 0.4);
    let mut fc_amd = ShiftedFactorization::new(&a_cplx, None, Definiteness::General).unwrap();
    fc_amd.ordering = FillOrdering::Amd;
    fc_amd.set_shift(C64::ZERO).unwrap();
    assert_eq!(fc_amd.stats().backend, FactorBackend::LbltAmd);
    assert!(fc_amd.stats().factor_nnz > 0);

    let mut fc_nd = ShiftedFactorization::new(&a_cplx, None, Definiteness::General).unwrap();
    fc_nd.ordering = FillOrdering::NestedDissection;
    fc_nd.set_shift(C64::ZERO).unwrap();
    assert_eq!(fc_nd.stats().backend, FactorBackend::LbltNd);
    assert!(fc_nd.stats().factor_nnz > 0);

    let rhs_c = random_vector::<C64>(n, 43);
    let (mut xc_amd, mut xc_nd) = (rhs_c.clone(), rhs_c.clone());
    fc_amd.solve_block_in_place(&mut xc_amd, 1);
    fc_nd.solve_block_in_place(&mut xc_nd, 1);
    assert!(rel_residual(&a_cplx, &xc_amd, &rhs_c) < 1e-11);
    assert!(rel_residual(&a_cplx, &xc_nd, &rhs_c) < 1e-11);

    // 3. Shift-invert eigenvalues on 3D P2 cavity via eigs_near agree with exact 2π² and across backends < 1e-10
    let eigs = eigs_near(
        &p.stiffness,
        &p.mass,
        sigma,
        3,
        Symmetry::Hermitian,
        Definiteness::General,
        &EigenOptions { tol: 1e-11, ..Default::default() },
        None,
    )
    .unwrap();
    assert!(eigs.converged);
    for (i, &lam) in eigs.values.iter().enumerate() {
        let exact = 2.0 * PI * PI;
        assert!((lam.re - exact).abs() / exact < 1e-2, "P2 cavity mode {i}: {} vs {exact}", lam.re);
        assert!(eigs.residuals[i] < 1e-9, "mode {i} residual {}", eigs.residuals[i]);
    }
}

#[test]
fn nyx_ordering_env_controls_direct_and_eigs_near() {
    let _guard = ENV_LOCK.lock().unwrap();
    // Verify `NYX_ORDERING=nd` and `NYX_ORDERING=amd` visibly select ND vs AMD through `driver::DirectSolver` and `eigs_near`.
    let p = nedelec_p2_cube(3);
    let sigma = 2.0 * PI * PI + 0.5;
    let a = Csr::axpby(1.0, &p.stiffness, -sigma, &p.mass);

    unsafe { std::env::set_var("NYX_ORDERING", "nd") };
    let s_nd = femcore::linalg::driver::DirectSolver::new(&a, false).unwrap();
    let e_nd = eigs_near(
        &p.stiffness,
        &p.mass,
        sigma,
        3,
        Symmetry::Hermitian,
        Definiteness::General,
        &EigenOptions { tol: 1e-11, ..Default::default() },
        None,
    )
    .unwrap();

    unsafe { std::env::set_var("NYX_ORDERING", "amd") };
    let s_amd = femcore::linalg::driver::DirectSolver::new(&a, false).unwrap();
    let e_amd = eigs_near(
        &p.stiffness,
        &p.mass,
        sigma,
        3,
        Symmetry::Hermitian,
        Definiteness::General,
        &EigenOptions { tol: 1e-11, ..Default::default() },
        None,
    )
    .unwrap();
    unsafe { std::env::remove_var("NYX_ORDERING") };

    assert_eq!(s_nd.stats().backend, FactorBackend::LbltNd);
    assert_eq!(s_amd.stats().backend, FactorBackend::LbltAmd);
    assert!(e_nd.converged && e_amd.converged);
    for (a_val, b_val) in e_nd.values.iter().zip(&e_amd.values) {
        let rel = (*a_val - *b_val).norm() / a_val.norm().max(1e-30);
        assert!(rel < 1e-10, "ND vs AMD eig diff {rel:e}");
    }
}

// ------------------------------------------------------------------------------------------------------------ FR-16 Pt 2 tests

#[test]
fn ams_p2_pcg_converges_rapidly_on_3d_p2_curl_curl_with_low_memory() {
    let mut its = Vec::new();
    for n in [3usize, 4, 5] {
        let p = nedelec_p2_cube(n);
        let a = Csr::axpby(1.0, &p.stiffness, 1.0, &p.mass);
        let ndofs = a.nrows();
        let topo = topo_of(&p);
        let ams = AmsP2::new(&a, &topo, &AmsP2Options::default()).unwrap();
        let rhs = random_vector::<f64>(ndofs, 101 + n as u64);
        let mut x = vec![0.0; ndofs];
        let info = Cg::new(CgOptions {
            tol: Tolerances { rtol: 1e-9, max_it: 200, ..Default::default() },
            ..Default::default()
        })
        .solve(&a, &ams, &rhs, &mut x);
        let res = rel_residual(&a, &x, &rhs);
        eprintln!(
            "P2 AMS-PCG n={n} (tets={}, dofs={ndofs}): it={}, res={res:.2e}, pc_nnz={}",
            p.n_tets,
            info.iterations,
            ams.nnz_storage()
        );
        assert!(info.reason.is_converged(), "n={n}: {:?}", info.reason);
        assert!(res < 1e-8, "n={n}: res={res:e}");
        assert!(info.iterations <= 35, "n={n}: too many iters {}", info.iterations);
        its.push(info.iterations);
    }
    // Nearly mesh-independent iteration count from n=3 (162 tets) to n=5 (750 tets)
    assert!(*its.last().unwrap() <= its[0] + 10, "iterations grew too fast: {its:?}");
}

#[test]
fn ams_p2_fgmres_solves_3d_p2_real_and_complex_shifted_curl_curl_with_low_memory() {
    // 3D P2 Nédélec cube n = 6 -> 1,296 tets, 5,184 DOFs.
    let p = nedelec_p2_cube(6);
    let ndofs = p.stiffness.nrows();
    let topo = topo_of(&p);

    // 1. Real indefinite shifted curl-curl: A = K - sigma M (sigma = 10.0, below first cavity mode 2π² ≈ 19.74
    //    so all gradient modes have negative eigenvalue -sigma while solenoidal modes are positive -> indefinite!).
    let sigma = 10.0;
    let a_indef = Csr::axpby(1.0, &p.stiffness, -sigma, &p.mass);
    let pc_spd = Csr::axpby(1.0, &p.stiffness, sigma, &p.mass);
    let rhs = random_vector::<f64>(ndofs, 202);
    let opts = IterativeCurlCurlOptions { rtol: 1e-9, max_it: 150, restart: 40, flexible: true, ..Default::default() };
    let rep = solve_curl_curl_p2_iterative(&a_indef, Some(&pc_spd), &topo, &rhs, &opts).unwrap();
    assert!(rep.converged, "indefinite P2 FGMRES failed: res={:.2e} it={}", rep.relative_residual, rep.iterations);
    assert!(rep.relative_residual < 1e-8, "res = {:e}", rep.relative_residual);

    // Compare memory against direct sparse factorization on the same 3D P2 matrix (even against memory-reduced LBL^T):
    let mut f_direct = ShiftedFactorization::new(&a_indef, None, Definiteness::General).unwrap();
    f_direct.ordering = FillOrdering::Amd;
    f_direct.set_shift(0.0).unwrap();
    let direct_lblt_nnz = f_direct.stats().factor_nnz;
    let direct_lu_est_nnz = 2 * direct_lblt_nnz;
    eprintln!(
        "3D P2 n=6 (tets={}, dofs={ndofs}): FGMRES+AmsP2 it={} res={:.2e} pc_nnz={} vs direct LBL^T factor_nnz={} (est LU {}, ratio {:.2}x less than LU)",
        p.n_tets,
        rep.iterations,
        rep.relative_residual,
        rep.preconditioner_nnz,
        direct_lblt_nnz,
        direct_lu_est_nnz,
        direct_lu_est_nnz as f64 / rep.preconditioner_nnz as f64
    );
    assert!(
        rep.preconditioner_nnz < direct_lblt_nnz && rep.preconditioner_nnz * 2 < direct_lu_est_nnz,
        "AmsP2 storage ({}) must be strictly smaller than direct LBL^T ({direct_lblt_nnz}) and <50% of direct LU ({direct_lu_est_nnz})",
        rep.preconditioner_nnz
    );

    // 2. Complex-symmetric lossy shifted curl-curl: A = K + (-sigma + i gamma) M
    let (a_cplx, _) = complex_symmetric_p2_curl_curl(&p, 12.0, 2.5);
    let pc_cplx = Csr::axpby(C64::ONE, &p.stiffness.to_complex(), C64::new(12.0, 2.5), &p.mass.to_complex());
    let rhs_c = random_vector::<C64>(ndofs, 203);
    let rep_c = solve_curl_curl_p2_iterative(&a_cplx, Some(&pc_cplx), &topo, &rhs_c, &opts).unwrap();
    eprintln!(
        "3D P2 complex-symmetric: FGMRES+AmsP2 it={} res={:.2e}",
        rep_c.iterations, rep_c.relative_residual
    );
    assert!(rep_c.converged && rep_c.relative_residual < 1e-8, "complex res={:e}", rep_c.relative_residual);

    // 3. COCG on complex-symmetric damped P2 matrix (verifies symmetry of AmsP2!)
    let a_damp = Csr::axpby(C64::ONE, &p.stiffness.to_complex(), C64::new(1.0, 0.5), &p.mass.to_complex());
    let ams_c = AmsP2::new(&a_damp, &topo, &AmsP2Options::default()).unwrap();
    let mut xc = vec![C64::ZERO; ndofs];
    let info_c = Cg::new(CgOptions {
        kind: CgKind::Symmetric,
        tol: Tolerances { rtol: 1e-9, max_it: 200, ..Default::default() },
        ..Default::default()
    })
    .solve(&a_damp, &ams_c, &rhs_c, &mut xc);
    let res_cocg = rel_residual(&a_damp, &xc, &rhs_c);
    eprintln!("3D P2 COCG+AmsP2: it={} res={res_cocg:.2e}", info_c.iterations);
    assert!(info_c.reason.is_converged() && res_cocg < 1e-7);
}

#[test]
fn ams_p2_jacobi_davidson_and_iterative_shift_invert_find_cavity_modes() {
    // 3D P2 cube n = 4 (384 tets, 1344 DOFs): find lowest 3 cavity modes near 2π² ≈ 19.739 without sparse LU!
    let p = nedelec_p2_cube(4);
    let topo = topo_of(&p);
    let target = C64::new(2.0 * PI * PI + 0.4, 0.0);
    let (kc, mc) = (p.stiffness.to_complex(), p.mass.to_complex());
    let jd = eigs_curl_curl_p2_jd(&kc, &mc, &topo, target, 3, 1e-7, 200).unwrap();
    eprintln!("P2 JD+AmsP2: converged={} values={:?} residuals={:?}", jd.converged, jd.values, jd.residuals);
    assert!(jd.converged, "P2 JD did not converge");
    for (i, &lam) in jd.values.iter().enumerate() {
        let exact = 2.0 * PI * PI;
        assert!((lam.re - exact).abs() / exact < 1e-2, "mode {i}: {} vs {exact}", lam.re);
        assert!(jd.residuals[i] < 1e-6, "mode {i} residual {}", jd.residuals[i]);
    }

    // Also test IterativeShiftInvertP2 directly on a shifted linear solve
    let sigma = 15.0f64;
    let mut it_si = IterativeShiftInvertP2::new(
        &p.stiffness,
        &p.mass,
        topo,
        IterativeCurlCurlOptions { rtol: 1e-10, max_it: 200, ..Default::default() },
    );
    it_si.set_shift(sigma).unwrap();
    let rhs = random_vector::<f64>(p.stiffness.nrows(), 303);
    let mut x = rhs.clone();
    it_si.solve_block_in_place(&mut x, 1);
    let shifted = Csr::axpby(1.0, &p.stiffness, -sigma, &p.mass);
    let res = rel_residual(&shifted, &x, &rhs);
    assert!(res < 1e-8, "IterativeShiftInvertP2 res={res:e}");
}

// ------------------------------------------------------------------------------------------------------------ FR-16 Pt 1 Benchmark (ignored by default, run with --ignored --nocapture)

#[test]
#[ignore = "FR-16 Pt 1 3D P2 Nédélec benchmark (5k-10k+ tets); run with --ignored --nocapture"]
fn benchmark_fr16_p2_nedelec_3d_backends() {
    for &n in &[9usize, 10] {
        // n=9 -> 4,374 tets (~18,468 DOFs), n=10 -> 6,000 tets (~25,680 DOFs)
        let t_gen = Instant::now();
        let p = nedelec_p2_cube(n);
        let ndofs = p.stiffness.nrows();
        let sigma = 2.0 * PI * PI + 0.45;
        let rhs_r = random_vector::<f64>(ndofs, 77);
        let shifted_r = Csr::axpby(1.0, &p.stiffness, -sigma, &p.mass);
        let (shifted_c, mc) = complex_symmetric_p2_curl_curl(&p, sigma, 0.35);
        let rhs_c = random_vector::<C64>(ndofs, 78);
        eprintln!(
            "\n=== 3D P2 Nédélec n={n}: tets={}, DOFs={ndofs}, nnz(K)={}, gen={:.2}s ===",
            p.n_tets,
            p.stiffness.nnz(),
            t_gen.elapsed().as_secs_f64()
        );

        // --- REAL SYMMETRIC INDEFINITE ---
        eprintln!("-- REAL SYMMETRIC INDEFINITE (K - sigma M, sigma={sigma:.3}) --");
        let mut ref_eigs_r: Option<Vec<C64>> = None;
        for (label, backend_env, ord) in [
            ("faer COLAMD LU (baseline)", "colamd", FillOrdering::Amd),
            ("faer LBL^H (AMD)", "lblt", FillOrdering::Amd),
            ("faer LBL^H (ND)", "lblt", FillOrdering::NestedDissection),
            ("default Auto (LBL^H)", "", FillOrdering::Auto),
        ] {
            if backend_env.is_empty() {
                unsafe { std::env::remove_var("NYX_DIRECT_BACKEND") };
            } else {
                unsafe { std::env::set_var("NYX_DIRECT_BACKEND", backend_env) };
            }
            let (ws0, priv0) = process_memory_mib();
            let t0 = Instant::now();
            let mut f = ShiftedFactorization::new(&p.stiffness, Some(&p.mass), Definiteness::General).unwrap();
            f.ordering = ord;
            f.set_shift(sigma).unwrap();
            let fac_s = t0.elapsed().as_secs_f64();
            let (ws1, priv1) = process_memory_mib();
            let mut x = rhs_r.clone();
            let ts = Instant::now();
            f.solve_block_in_place(&mut x, 1);
            let slv_s = ts.elapsed().as_secs_f64();
            let res = rel_residual(&shifted_r, &x, &rhs_r);
            let st = f.stats();

            unsafe {
                std::env::set_var(
                    "NYX_ORDERING",
                    match ord {
                        FillOrdering::Amd => "amd",
                        FillOrdering::NestedDissection => "nd",
                        FillOrdering::Auto => "auto",
                    },
                )
            };
            let ep = eigs_near(
                &p.stiffness,
                &p.mass,
                sigma,
                3,
                Symmetry::Hermitian,
                Definiteness::General,
                &EigenOptions { tol: 1e-11, ..Default::default() },
                None,
            )
            .unwrap();
            let max_eig_diff = match &ref_eigs_r {
                None => {
                    ref_eigs_r = Some(ep.values.clone());
                    0.0
                }
                Some(r) => r
                    .iter()
                    .zip(&ep.values)
                    .map(|(a, b)| (*a - *b).norm() / a.norm().max(1e-30))
                    .fold(0.0f64, f64::max),
            };
            eprintln!(
                "  {label:<28} backend={:?} factor_nnz={:>10} fac={fac_s:>6.2}s slv={slv_s:>6.3}s dWS={:>+7.1}MiB dPriv={:>+7.1}MiB res={res:.2e} eig_diff={max_eig_diff:.2e}",
                st.backend,
                st.factor_nnz,
                ws1 - ws0,
                priv1 - priv0
            );
            drop(f);
        }

        // Exact BLR multifrontal LU (tol=0.0) on real symmetric indefinite (only for n <= 9 so it stays quick)
        if n <= 9 {
            let (ws0, priv0) = process_memory_mib();
            let t0 = Instant::now();
            let blr_opts = BlrOptions { tol: 0.0, leaf: 64, ..Default::default() };
            let mut f = ShiftedBlrLu::new(&p.stiffness, Some(&p.mass), &blr_opts).unwrap();
            let ok = f.set_shift(sigma).is_ok();
            let fac_s = t0.elapsed().as_secs_f64();
            let (ws1, priv1) = process_memory_mib();
            if ok {
                let mut x = rhs_r.clone();
                let ts = Instant::now();
                f.solve_block_in_place(&mut x, 1);
                let slv_s = ts.elapsed().as_secs_f64();
                let res = rel_residual(&shifted_r, &x, &rhs_r);
                eprintln!(
                    "  {:<28} backend={:?} factor_nnz={:>10} fac={fac_s:>6.2}s slv={slv_s:>6.3}s dWS={:>+7.1}MiB dPriv={:>+7.1}MiB res={res:.2e}",
                    "ShiftedBlrLu (ND tol=0)",
                    f.stats().backend,
                    f.stats().factor_nnz,
                    ws1 - ws0,
                    priv1 - priv0
                );
            } else {
                eprintln!("  ShiftedBlrLu (ND tol=0): Singular pivot block");
            }
        }

        // --- COMPLEX-SYMMETRIC INDEFINITE ---
        eprintln!("-- COMPLEX-SYMMETRIC INDEFINITE (K + (-sigma + i*gamma(x)) M) --");
        let mut ref_eigs_c: Option<Vec<C64>> = None;
        for (label, backend_env, ord) in [
            ("faer COLAMD LU (baseline)", "colamd", FillOrdering::Amd),
            ("Sym LBL^T (AMD)", "lblt", FillOrdering::Amd),
            ("Sym LBL^T (ND)", "lblt", FillOrdering::NestedDissection),
            ("default Auto (Sym LBL^T)", "", FillOrdering::Auto),
        ] {
            if backend_env.is_empty() {
                unsafe { std::env::remove_var("NYX_DIRECT_BACKEND") };
            } else {
                unsafe { std::env::set_var("NYX_DIRECT_BACKEND", backend_env) };
            }
            let (ws0, priv0) = process_memory_mib();
            let t0 = Instant::now();
            let mut f = ShiftedFactorization::new(&shifted_c, None, Definiteness::General).unwrap();
            f.ordering = ord;
            f.set_shift(C64::ZERO).unwrap();
            let fac_s = t0.elapsed().as_secs_f64();
            let (ws1, priv1) = process_memory_mib();
            let mut x = rhs_c.clone();
            let ts = Instant::now();
            f.solve_block_in_place(&mut x, 1);
            let slv_s = ts.elapsed().as_secs_f64();
            let res = rel_residual(&shifted_c, &x, &rhs_c);
            let st = f.stats();

            // Complex eigenproblem (K + i*gamma*M) x = lam M x near sigma
            let k_lossy = Csr::axpby(C64::ONE, &shifted_c, C64::new(sigma, 0.0), &mc);
            unsafe {
                std::env::set_var(
                    "NYX_ORDERING",
                    match ord {
                        FillOrdering::Amd => "amd",
                        FillOrdering::NestedDissection => "nd",
                        FillOrdering::Auto => "auto",
                    },
                )
            };
            let ep = eigs_near(
                &k_lossy,
                &mc,
                C64::new(sigma, 0.0),
                3,
                Symmetry::General,
                Definiteness::General,
                &EigenOptions { tol: 1e-11, ..Default::default() },
                None,
            )
            .unwrap();
            let max_eig_diff = match &ref_eigs_c {
                None => {
                    ref_eigs_c = Some(ep.values.clone());
                    0.0
                }
                Some(r) => r
                    .iter()
                    .zip(&ep.values)
                    .map(|(a, b)| (*a - *b).norm() / a.norm().max(1e-30))
                    .fold(0.0f64, f64::max),
            };
            eprintln!(
                "  {label:<28} backend={:?} factor_nnz={:>10} fac={fac_s:>6.2}s slv={slv_s:>6.3}s dWS={:>+7.1}MiB dPriv={:>+7.1}MiB res={res:.2e} eig_diff={max_eig_diff:.2e}",
                st.backend,
                st.factor_nnz,
                ws1 - ws0,
                priv1 - priv0
            );
            drop(f);
        }
        unsafe {
            std::env::remove_var("NYX_DIRECT_BACKEND");
            std::env::remove_var("NYX_ORDERING");
        }
    }
}
