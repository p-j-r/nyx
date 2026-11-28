// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-43 Part A verification suite:
//! - `Vector` trait & `CpuVec<T>` 0-ULP bit-for-bit parity with `vecops`, `Cg`, `Gmres`, and `lobpcg`
//! - Atomic-free fixed-order two-level `tree_reduce_sum` (`BLOCK = 256`) and thread determinism
//! - Parallel deterministic Chebyshev `Amg`, mixed-precision `AmgF32` + `AmsF32` driving outer `f64` PCG/FGMRES
//! - `SellCSigma<T>` (`C = 32, sigma = 128`) permutation invariants & 0-ULP SpMV parity
//! - `Csr64<T>` `u64` sparse indexing above `u32::MAX`, `from_triplets64`, checked `Csr <-> Csr64`
//! - Non-finite input guards (`finite_guard`) and `NYX_THREADS=1` vs `8` subprocess determinism
//! - Opt-in `--features gpu` (`wgpu` + WGSL on discrete GPU): `GpuVecF32`, atomic-free tree reductions,
//!   `GpuSellCSigmaF32` SpMV, `GpuAmgF32` Chebyshev V-cycle driving outer `f64` PCG/FGMRES, and bitwise GPU determinism

use femcore::linalg::amg::{Amg, AmgF32, AmgOptions, SmootherKind};
use femcore::linalg::ams::{Ams, AmsF32, AmsOptions};
use femcore::linalg::csr::{Csr, CsrError, laplacian_2d, laplacian_3d};
use femcore::linalg::device::{
    CpuVec, Csr64, DeviceIdentity, HostOp, HostPc, SellCSigma, Vector, cg_vec, fgmres_vec, lobpcg_vec,
    tree_reduce_dot_f32, tree_reduce_dotc, tree_reduce_sum, tree_reduce_sum_f32,
};
use femcore::linalg::gpu::{gpu_available, probe_gpu_device};
use femcore::linalg::ksp::{Cg, CgKind, CgNorm, CgOptions, Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
use femcore::linalg::linop::Identity;
use femcore::linalg::lobpcg::{LobpcgOptions, LobpcgReason, lobpcg};
use femcore::linalg::maxwell::nedelec_cube;
use femcore::linalg::rng::{SplitMix64, random_vector};
use femcore::linalg::scalar::{C64, Scalar};
use femcore::linalg::vecops;

fn assert_f64_bits_eq(a: &[f64], b: &[f64], ctx: &str) {
    assert_eq!(a.len(), b.len(), "{ctx}: length mismatch");
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert_eq!(x.to_bits(), y.to_bits(), "{ctx}: f64 ULP mismatch at index {i}: {x:e} vs {y:e}");
    }
}

fn assert_c64_bits_eq(a: &[C64], b: &[C64], ctx: &str) {
    assert_eq!(a.len(), b.len(), "{ctx}: length mismatch");
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        assert_eq!(x.re.to_bits(), y.re.to_bits(), "{ctx}: re ULP mismatch at {i}: {} vs {}", x.re, y.re);
        assert_eq!(x.im.to_bits(), y.im.to_bits(), "{ctx}: im ULP mismatch at {i}: {} vs {}", x.im, y.im);
    }
}

fn true_rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut r = vec![T::ZERO; b.len()];
    a.spmv(x, &mut r);
    vecops::aypx(-T::ONE, b, &mut r);
    vecops::norm2(&r) / vecops::norm2(b).max(1e-300)
}

fn irregular_spd_csr(n: usize) -> Csr<f64> {
    let mut rng = SplitMix64::new(0x5E11_C516);
    let mut triplets = Vec::new();
    for i in 0..n {
        // Alternate dense hub rows (up to 48 off-diagonals) with sparse rows (2 off-diagonals)
        // inside every 128-row window so sigma-sorting significantly reduces slice padding.
        let deg = if i % 17 == 0 { 24 } else if i % 5 == 0 { 9 } else { 1 };
        let mut diag = 4.0;
        for _ in 0..deg {
            let j = (rng.next_u64() as usize) % n;
            if j != i {
                let w: f64 = -0.05 * (1.5 + 0.5 * rng.next_signed());
                triplets.push((i, j, w));
                triplets.push((j, i, w));
                diag += 2.0 * w.abs();
            }
        }
        triplets.push((i, i, diag + 2.0));
    }
    Csr::from_triplets(n, n, &triplets).expect("valid irregular SPD CSR")
}

#[test]
fn cpu_vec_zero_ulp_parity_with_vecops_krylov_and_lobpcg() {
    let n = 3 * vecops::CHUNK + 137;
    let x_raw = random_vector::<f64>(n, 101);
    let y_raw = random_vector::<f64>(n, 202);
    let xc_raw = random_vector::<C64>(n, 303);
    let yc_raw = random_vector::<C64>(n, 404);

    let xv = CpuVec::new(x_raw.clone());
    let mut yv = CpuVec::new(y_raw.clone());
    assert_eq!(xv.dotc(&yv).to_bits(), vecops::dotc(&x_raw, &y_raw).to_bits());
    assert_eq!(xv.dotu(&yv).to_bits(), vecops::dotu(&x_raw, &y_raw).to_bits());
    assert_eq!(xv.norm2_sq().to_bits(), vecops::norm2_sq(&x_raw).to_bits());
    assert_eq!(xv.norm2().to_bits(), vecops::norm2(&x_raw).to_bits());

    let mut y_ref = y_raw.clone();
    yv.axpy(-0.375, &xv);
    vecops::axpy(-0.375, &x_raw, &mut y_ref);
    assert_f64_bits_eq(&yv, &y_ref, "CpuVec::axpy");

    yv.aypx(1.25, &xv);
    vecops::aypx(1.25, &x_raw, &mut y_ref);
    assert_f64_bits_eq(&yv, &y_ref, "CpuVec::aypx");

    yv.scale(-0.75);
    vecops::scale(-0.75, &mut y_ref);
    assert_f64_bits_eq(&yv, &y_ref, "CpuVec::scale");

    let xcv = CpuVec::new(xc_raw.clone());
    let mut ycv = CpuVec::new(yc_raw.clone());
    let dc_got = xcv.dotc(&ycv);
    let dc_ref = vecops::dotc(&xc_raw, &yc_raw);
    assert_eq!(dc_got.re.to_bits(), dc_ref.re.to_bits());
    assert_eq!(dc_got.im.to_bits(), dc_ref.im.to_bits());
    let mut yc_ref = yc_raw.clone();
    let alpha_c = C64::new(0.4, -0.65);
    ycv.axpy(alpha_c, &xcv);
    vecops::axpy(alpha_c, &xc_raw, &mut yc_ref);
    assert_c64_bits_eq(&ycv, &yc_ref, "CpuVec<C64>::axpy");

    // 1. CG & AMG-PCG 0-ULP parity (real SPD + complex symmetric COCG)
    let a = laplacian_2d(48, 45);
    let m_amg = Amg::new(&a, AmgOptions::default()).expect("Amg setup");
    let b = random_vector::<f64>(a.nrows(), 77);
    let cg_opts = CgOptions {
        tol: Tolerances { rtol: 1e-10, max_it: 200, ..Default::default() },
        record_history: true,
        ..Default::default()
    };
    let mut x_slice = vec![0.0; a.nrows()];
    let info_slice = Cg::new(cg_opts.clone()).solve(&a, &m_amg, &b, &mut x_slice);

    let b_vec = CpuVec::new(b.clone());
    let mut x_vec = CpuVec::zeros(a.nrows());
    let info_vec = cg_vec(&a, &m_amg, &b_vec, &mut x_vec, &cg_opts);
    assert_eq!(info_slice.iterations, info_vec.iterations);
    assert_f64_bits_eq(&info_slice.history, &info_vec.history, "cg_vec history");
    assert_f64_bits_eq(&x_slice, &x_vec, "cg_vec solution");

    // Complex-symmetric COCG parity
    let a_cs = Csr::axpby(
        C64::ONE,
        &laplacian_2d(30, 30).to_complex(),
        C64::new(0.1, 0.35),
        &Csr::identity(900),
    );
    let amg_cs = Amg::new(&a_cs, AmgOptions::default()).expect("complex Amg");
    let b_cs = random_vector::<C64>(900, 88);
    let cocg_opts = CgOptions {
        tol: Tolerances { rtol: 1e-9, max_it: 200, ..Default::default() },
        kind: CgKind::Symmetric,
        record_history: true,
        ..Default::default()
    };
    let mut x_cs_slice = vec![C64::ZERO; 900];
    let info_cs_slice = Cg::new(cocg_opts.clone()).solve(&a_cs, &amg_cs, &b_cs, &mut x_cs_slice);
    let mut x_cs_vec = CpuVec::zeros(900);
    let info_cs_vec = cg_vec(&a_cs, &amg_cs, &CpuVec::new(b_cs), &mut x_cs_vec, &cocg_opts);
    assert_eq!(info_cs_slice.iterations, info_cs_vec.iterations);
    assert_c64_bits_eq(&x_cs_slice, &x_cs_vec, "cocg_vec solution");

    // 2. FGMRES / GMRES 0-ULP parity across orthogonalization strategies
    for orthog in [
        Orthogonalization::ClassicalGramSchmidt,
        Orthogonalization::ClassicalGramSchmidtIfNeeded,
        Orthogonalization::ModifiedGramSchmidt,
    ] {
        let gmres_opts = GmresOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 120, ..Default::default() },
            restart: 15,
            side: PcSide::Right,
            orthog,
            flexible: true,
            record_history: true,
            ..Default::default()
        };
        let mut xg_slice = vec![0.0; a.nrows()];
        let ig_slice = Gmres::new(gmres_opts.clone()).solve(&a, &m_amg, &b, &mut xg_slice);
        let mut xg_vec = CpuVec::zeros(a.nrows());
        let ig_vec = fgmres_vec(&a, &m_amg, &b_vec, &mut xg_vec, &gmres_opts);
        assert_eq!(ig_slice.iterations, ig_vec.iterations);
        assert_f64_bits_eq(&ig_slice.history, &ig_vec.history, "fgmres_vec history");
        assert_f64_bits_eq(&xg_slice, &xg_vec, "fgmres_vec solution");
    }

    // 3. LOBPCG 0-ULP parity (generalised problem with mass matrix, AMG preconditioner, and deflated mode)
    let a_eig = laplacian_2d(14, 14); // 196 unknowns
    let n_eig = a_eig.nrows();
    let b_eig = Csr::axpby(0.5, &Csr::identity(n_eig), 0.05, &a_eig);
    let pc_eig = Amg::new(&a_eig, AmgOptions::default()).expect("2D Amg");
    let mut lopts0 = LobpcgOptions::new(1);
    lopts0.tol = 1e-11;
    let mode0 = lobpcg(&a_eig, Some(&b_eig), Some(&pc_eig), None, &[], &lopts0).expect("mode0");
    assert_eq!(mode0.reason, LobpcgReason::Converged);
    let c0 = mode0.vectors[0].clone();

    let mut lopts = LobpcgOptions::new(3);
    lopts.tol = 1e-9;
    lopts.max_it = 120;

    let res_host = lobpcg(&a_eig, Some(&b_eig), Some(&pc_eig), None, std::slice::from_ref(&c0), &lopts)
        .expect("host lobpcg");
    let proto = CpuVec::zeros(n_eig);
    let cons_vec = [CpuVec::new(c0)];
    let res_dev = lobpcg_vec(
        &HostOp(&a_eig),
        Some(&HostOp(&b_eig)),
        Some(&HostPc(&pc_eig)),
        &proto,
        None,
        &cons_vec,
        &lopts,
    )
    .expect("lobpcg_vec");
    assert_eq!(res_dev.reason, LobpcgReason::Converged);
    assert_eq!(res_host.iterations, res_dev.iterations);
    assert_f64_bits_eq(&res_host.values, &res_dev.values, "lobpcg_vec eigenvalues");
    assert_f64_bits_eq(&res_host.residuals, &res_dev.residuals, "lobpcg_vec residuals");
    for (vh, vd) in res_host.vectors.iter().zip(&res_dev.vectors) {
        assert_f64_bits_eq(vh, vd, "lobpcg_vec eigenvector");
    }

    // Exercise Identity / DeviceIdentity
    let mut xi = CpuVec::zeros(a.nrows());
    let info_id = cg_vec(&a, &DeviceIdentity, &b_vec, &mut xi, &CgOptions::default());
    let mut xi2 = vec![0.0; a.nrows()];
    let info_id2 = Cg::new(CgOptions::default()).solve(&a, &Identity, &b, &mut xi2);
    assert_eq!(info_id.iterations, info_id2.iterations);
}

#[test]
fn tree_reduce_sum_atomic_free_multi_level_and_thread_determinism() {
    // 70_000 > 256 * 256 = 65_536 -> exercises 3 full binary-tree reduction levels!
    for &n in &[0usize, 1, 255, 256, 257, 65_536, 70_123] {
        let xs = random_vector::<f64>(n, 0xABCD + n as u64);
        let ys = random_vector::<f64>(n, 0x1234 + n as u64);
        let got_sum = tree_reduce_sum(&xs);
        let naive_sum: f64 = xs.iter().copied().sum();
        assert!((got_sum - naive_sum).abs() <= 1e-10 * (1.0 + naive_sum.abs()), "n={n}");

        let got_dot = tree_reduce_dotc(&xs, &ys);
        let naive_dot: f64 = xs.iter().zip(&ys).map(|(a, b)| a * b).sum();
        assert!((got_dot - naive_dot).abs() <= 1e-10 * (1.0 + naive_dot.abs()), "n={n}");

        let xs32: Vec<f32> = xs.iter().map(|&v| v as f32).collect();
        let ys32: Vec<f32> = ys.iter().map(|&v| v as f32).collect();
        let got_s32 = tree_reduce_sum_f32(&xs32);
        let got_d32 = tree_reduce_dot_f32(&xs32, &ys32);
        assert!(((got_s32 as f64) - naive_sum).abs() <= 1e-3 * (1.0 + naive_sum.abs()), "f32 sum n={n}");
        assert!(((got_d32 as f64) - naive_dot).abs() <= 1e-3 * (1.0 + naive_dot.abs()), "f32 dot n={n}");
    }

    // 1 vs 8 thread 0-ULP determinism on large vector (> PAR_MIN_ELEMS = 131_072)
    let big_n = 150_019;
    let x_big = random_vector::<C64>(big_n, 42);
    let y_big = random_vector::<C64>(big_n, 43);
    let run = |nt: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(nt).build().expect("pool");
        pool.install(|| (tree_reduce_sum(&x_big), tree_reduce_dotc(&x_big, &y_big)))
    };
    let (s1, d1) = run(1);
    let (s8, d8) = run(8);
    assert_eq!(s1.re.to_bits(), s8.re.to_bits());
    assert_eq!(s1.im.to_bits(), s8.im.to_bits());
    assert_eq!(d1.re.to_bits(), d8.re.to_bits());
    assert_eq!(d1.im.to_bits(), d8.im.to_bits());
}

#[test]
fn chebyshev_amg_and_mixed_precision_amg_f32_ams_f32() {
    assert_eq!(AmgOptions::default().smoother, SmootherKind::GaussSeidel);

    let a = laplacian_3d(20); // 8,000 unknowns
    let b = random_vector::<f64>(a.nrows(), 17);
    let cheb_opts = AmgOptions { smoother: SmootherKind::Chebyshev { degree: 3 }, ..Default::default() };

    // 1. Double-precision Chebyshev AMG is symmetric and converges in O(1) iterations
    let amg_cheb = Amg::new(&a, cheb_opts.clone()).expect("Chebyshev AMG");
    let u = random_vector::<f64>(a.nrows(), 21);
    let v = random_vector::<f64>(a.nrows(), 22);
    let (mut mu, mut mv) = (vec![0.0; a.nrows()], vec![0.0; a.nrows()]);
    femcore::linalg::linop::Preconditioner::apply(&amg_cheb, &u, &mut mu);
    femcore::linalg::linop::Preconditioner::apply(&amg_cheb, &v, &mut mv);
    let uv = vecops::dotc(&u, &mv);
    let vu = vecops::dotc(&v, &mu);
    assert!((uv - vu).abs() < 1e-12 * uv.abs(), "Chebyshev AMG must be symmetric: {uv} vs {vu}");

    let cg_opts = CgOptions {
        tol: Tolerances { rtol: 1e-11, max_it: 150, ..Default::default() },
        norm: CgNorm::Unpreconditioned,
        ..Default::default()
    };
    let mut x_cheb = vec![0.0; a.nrows()];
    let info_cheb = Cg::new(cg_opts.clone()).solve(&a, &amg_cheb, &b, &mut x_cheb);
    assert!(info_cheb.reason.is_converged() && info_cheb.iterations <= 22, "{info_cheb:?}");
    assert!(true_rel_residual(&a, &x_cheb, &b) < 1e-10);

    // 2. Mixed-precision AmgF32 (f32 V-cycle) drives outer f64 PCG and FGMRES to < 1e-10
    let amg_f32 = AmgF32::new_chebyshev(&a, 3).expect("AmgF32 Chebyshev");
    assert_eq!(amg_f32.stats().level_sizes, amg_cheb.stats().level_sizes);

    let mut x_f32_pcg = vec![0.0; a.nrows()];
    let info_f32_pcg = Cg::new(cg_opts.clone()).solve(&a, &amg_f32, &b, &mut x_f32_pcg);
    assert!(info_f32_pcg.reason.is_converged() && info_f32_pcg.iterations <= 24, "{info_f32_pcg:?}");
    let rel_pcg = true_rel_residual(&a, &x_f32_pcg, &b);
    assert!(rel_pcg < 1e-10, "AmgF32 + f64 PCG relative residual {rel_pcg:e}");

    // Verify tiny-norm RHS (1e-28) does not underflow in f32 V-cycle
    let b_tiny: Vec<f64> = b.iter().map(|&v| v * 1e-28).collect();
    let mut x_tiny = vec![0.0; a.nrows()];
    let info_tiny = Cg::new(cg_opts.clone()).solve(&a, &amg_f32, &b_tiny, &mut x_tiny);
    assert!(info_tiny.reason.is_converged());
    assert!(true_rel_residual(&a, &x_tiny, &b_tiny) < 1e-10);

    // Outer f64 FGMRES with AmgF32
    let fg_opts = GmresOptions {
        tol: Tolerances { rtol: 1e-11, max_it: 100, ..Default::default() },
        restart: 25,
        side: PcSide::Right,
        flexible: true,
        orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
        ..Default::default()
    };
    let mut x_fg = CpuVec::zeros(a.nrows());
    let info_fg = fgmres_vec(&a, &amg_f32, &CpuVec::new(b.clone()), &mut x_fg, &fg_opts);
    assert!(info_fg.reason.is_converged());
    assert!(true_rel_residual(&a, &x_fg, &b) < 1e-10);

    // 3. Chebyshev Ams<f64> and mixed-precision AmsF32 on 3D Nedelec H(curl)
    let cube = nedelec_cube(5);
    let a_curl = Csr::axpby(1.0, &cube.stiffness, 1.0, &cube.mass);
    let b_curl = random_vector::<f64>(a_curl.nrows(), 91);
    let ams_cheb_opts = AmsOptions {
        amg: cheb_opts,
        smoother_sweeps: 1,
        smoother: SmootherKind::Chebyshev { degree: 3 },
    };
    let ams_cheb = Ams::new(&a_curl, &cube.gradient, &cube.edge_vectors, &ams_cheb_opts).expect("Chebyshev Ams");
    let ams_f32 = AmsF32::new(&a_curl, &cube.gradient, &cube.edge_vectors, &ams_cheb_opts).expect("AmsF32");

    let mut x_ams1 = vec![0.0; a_curl.nrows()];
    let i_ams1 = Cg::new(cg_opts.clone()).solve(&a_curl, &ams_cheb, &b_curl, &mut x_ams1);
    assert!(i_ams1.reason.is_converged() && true_rel_residual(&a_curl, &x_ams1, &b_curl) < 1e-10);

    let mut x_ams32 = vec![0.0; a_curl.nrows()];
    let i_ams32 = Cg::new(cg_opts).solve(&a_curl, &ams_f32, &b_curl, &mut x_ams32);
    assert!(i_ams32.reason.is_converged() && true_rel_residual(&a_curl, &x_ams32, &b_curl) < 1e-10);
}

#[test]
fn sell_c_sigma_permutation_fill_reduction_and_zero_ulp_spmv() {
    let a = irregular_spd_csr(4_500);
    let sell = SellCSigma::from_csr(&a).expect("SELL-32-128");
    let sell_unsorted = SellCSigma::from_csr_with(&a, 32, 32).expect("SELL-32-32");
    assert_eq!(sell.c(), 32);
    assert_eq!(sell.sigma(), 128);
    assert_eq!(sell.nrows(), a.nrows());
    assert_eq!(sell.nnz(), a.nnz());
    // Sorting over sigma = 128 groups hub rows into dedicated slices, reducing padding vs sigma = 32.
    assert!(
        sell.stored_entries() < sell_unsorted.stored_entries(),
        "sigma=128 stored {} vs sigma=32 stored {}",
        sell.stored_entries(),
        sell_unsorted.stored_entries()
    );
    // Verify bijection and window confinement of `perm` and `inv_perm`.
    for r in 0..a.nrows() {
        let sr = sell.inv_perm()[r] as usize;
        assert_eq!(sell.perm()[sr] as usize, r);
        assert_eq!(sr / 128, r / 128, "permutation escaped sigma=128 window");
    }

    let x = random_vector::<f64>(a.ncols(), 555);
    let (mut y_csr, mut y_sell) = (vec![0.0; a.nrows()], vec![0.0; a.nrows()]);
    a.spmv(&x, &mut y_csr);
    sell.spmv(&x, &mut y_sell);
    assert_f64_bits_eq(&y_csr, &y_sell, "SellCSigma vs Csr SpMV");

    // Complex matrix + parallel threshold check
    let a_big = laplacian_3d(32).to_complex(); // 32,768 rows, > 200k nnz (> PAR_MIN_ELEMS)
    let sell_big = SellCSigma::from_csr(&a_big).expect("complex SELL-32-128");
    let xc = random_vector::<C64>(a_big.ncols(), 777);
    let mut yc_csr = vec![C64::ZERO; a_big.nrows()];
    a_big.spmv(&xc, &mut yc_csr);

    let run = |nt: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(nt).build().expect("pool");
        pool.install(|| {
            let mut y = vec![C64::ZERO; a_big.nrows()];
            sell_big.spmv(&xc, &mut y);
            y
        })
    };
    let y1 = run(1);
    let y8 = run(8);
    assert_c64_bits_eq(&yc_csr, &y1, "SellCSigma<C64> vs Csr<C64>");
    assert_c64_bits_eq(&y1, &y8, "SellCSigma<C64> 1 vs 8 threads");
}

#[test]
fn csr64_u64_indexing_above_u32_max_roundtrip_and_parallel_spmv() {
    let huge_col_0 = (u32::MAX as u64) + 1;
    let huge_col_1 = (u32::MAX as u64) + 4_096;
    let ncols_64 = huge_col_1 + 100;
    let triplets = vec![
        (0u64, huge_col_1, 2.5f64),
        (0u64, 7u64, -1.0f64),
        (0u64, huge_col_1, 0.5f64), // duplicate summed -> 3.0
        (2u64, huge_col_0, 4.25f64),
    ];
    let big = Csr64::from_triplets64(4, ncols_64, &triplets).expect("Csr64 above u32::MAX");
    assert_eq!(big.nrows(), 4);
    assert_eq!(big.ncols(), ncols_64);
    assert_eq!(big.nnz(), 3);
    assert!(big.exceeds_u32_limits());
    assert_eq!(big.row(0), (&[7u64, huge_col_1][..], &[-1.0, 3.0][..]));
    assert_eq!(big.row(1), (&[][..], &[][..]));
    assert_eq!(big.row(2), (&[huge_col_0][..], &[4.25][..]));
    assert_eq!(Csr::try_from(&big), Err(CsrError::DimensionTooLarge));

    // Lossless round-trip and 0-ULP parallel SpMV on 3D Laplacian
    let a32 = laplacian_3d(32); // > 220k nnz
    let a64 = Csr64::from(&a32);
    assert!(!a64.exceeds_u32_limits());
    let back32 = Csr::try_from(&a64).expect("narrowing within u32::MAX");
    assert_eq!(a32, back32);

    let x = random_vector::<f64>(a32.ncols(), 314);
    let mut y32 = vec![0.0; a32.nrows()];
    a32.spmv(&x, &mut y32);

    let run = |nt: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(nt).build().expect("pool");
        pool.install(|| {
            let mut y = vec![0.0; a32.nrows()];
            a64.spmv_checked(&x, &mut y).expect("spmv_checked");
            y
        })
    };
    let y1 = run(1);
    let y8 = run(8);
    assert_f64_bits_eq(&y32, &y1, "Csr64 vs Csr SpMV");
    assert_f64_bits_eq(&y1, &y8, "Csr64 1 vs 8 threads");
}

#[test]
fn finite_guard_rejects_non_finite_matrices_without_panic() {
    let bad_nan = Csr::from_triplets(3, 3, &[(0, 0, 2.0), (1, 1, f64::NAN), (2, 2, 2.0)]).expect("csr");
    let bad_inf = Csr::from_triplets(3, 3, &[(0, 0, 2.0), (1, 1, f64::INFINITY), (2, 2, 2.0)]).expect("csr");
    assert!(Amg::new(&bad_nan, AmgOptions::default()).is_err());
    assert!(Amg::new(&bad_inf, AmgOptions::default()).is_err());
    assert!(AmgF32::new(&bad_nan, AmgOptions::default()).is_err());
    assert!(AmgF32::new(&bad_inf, AmgOptions::default()).is_err());

    let g = Csr::from_triplets(3, 2, &[(0, 0, 1.0), (1, 1, 1.0), (2, 0, -1.0)]).expect("g");
    let ev = [[1.0, 0.0, 0.0]; 3];
    assert!(AmsF32::new(&bad_nan, &g, &ev, &AmsOptions::default()).is_err());
}

/// Subprocess worker entry invoked by `subprocess_nyx_threads_1_vs_8_determinism`.
#[test]
fn subprocess_nyx_threads_1_vs_8_determinism() {
    if std::env::var("NYX_FR43A_SUBPROC").as_deref() == Ok("1") {
        let a = laplacian_3d(26); // 17,576 rows, ~118k nnz; exercise parallel chunks
        let sell = SellCSigma::from_csr(&a).expect("sell");
        let amg = AmgF32::new_chebyshev(&a, 3).expect("amg_f32");
        let b = random_vector::<f64>(a.nrows(), 2026);
        let mut x = CpuVec::zeros(a.nrows());
        let opts = CgOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 80, ..Default::default() },
            ..Default::default()
        };
        let info = cg_vec(&sell, &amg, &CpuVec::new(b), &mut x, &opts);
        let sum_bits = tree_reduce_sum(x.as_slice()).to_bits();
        let nrm_bits = x.norm2_sq().to_bits();
        println!("FR43A_DIGEST={} {:016x} {:016x} {:016x}", info.iterations, sum_bits, nrm_bits, x[123].to_bits());
        return;
    }

    let exe = std::env::current_exe().expect("test binary path");
    let invoke = |threads: &str| -> String {
        let out = std::process::Command::new(&exe)
            .args(["--exact", "subprocess_nyx_threads_1_vs_8_determinism", "--nocapture"])
            .env("NYX_FR43A_SUBPROC", "1")
            .env("NYX_THREADS", threads)
            .output()
            .expect("spawn subprocess");
        assert!(out.status.success(), "subprocess failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        stdout
            .lines()
            .find(|l| l.starts_with("FR43A_DIGEST="))
            .expect("digest line in subprocess stdout")
            .to_string()
    };

    let d1 = invoke("1");
    let d8 = invoke("8");
    assert_eq!(d1, d8, "NYX_THREADS=1 vs 8 digest mismatch");
}

#[cfg(not(feature = "gpu"))]
#[test]
fn default_cpu_build_reports_gpu_disabled() {
    assert!(!gpu_available());
    assert!(probe_gpu_device().is_none());
}

// =======================================================================================
// Real Discrete GPU (`--features gpu`) Compute Kernel & Solver Tests
// =======================================================================================

#[cfg(feature = "gpu")]
mod gpu_tests {
    use super::*;
    use femcore::linalg::gpu::{GpuAmgF32, GpuContext, GpuSellCSigmaF32, GpuSellCSigmaF64, GpuVecF32, GpuVecF64};

    #[test]
    fn gpu_context_selects_discrete_adapter() {
        assert!(gpu_available(), "expected discrete wgpu adapter to be available");
        let info = probe_gpu_device().expect("probe_gpu_device");
        let ctx = GpuContext::new().expect("GpuContext::new");
        println!(
            "Selected GPU adapter: {} (vendor=0x{:04x}, type={}, backend={}, discrete={}, f64={})",
            ctx.info().name,
            ctx.info().vendor,
            ctx.info().device_type,
            ctx.info().backend,
            ctx.info().is_discrete,
            ctx.info().supports_f64
        );
        assert_eq!(ctx.info(), &info);
        assert_eq!(ctx.supports_f64(), info.supports_f64);
        assert!(ctx.info().is_discrete, "expected discrete GPU (e.g. NVIDIA T400), got {:?}", ctx.info());
        // Verify shared device/queue handles are exposed cleanly.
        let _dev = ctx.wgpu_device();
        let _q = ctx.wgpu_queue();
    }

    #[test]
    fn gpu_vec_f32_blas1_and_atomic_free_tree_reduce_parity_and_determinism() {
        let ctx = GpuContext::shared().expect("GpuContext::shared");
        // 70_123 spans 3 levels of the BLOCK = 256 reduction tree (70_123 -> 274 -> 2 -> 1).
        let n = 70_123usize;
        let x64 = random_vector::<f64>(n, 1001);
        let y64 = random_vector::<f64>(n, 1002);
        let x32: Vec<f32> = x64.iter().map(|&v| v as f32).collect();
        let y32: Vec<f32> = y64.iter().map(|&v| v as f32).collect();

        let gx = GpuVecF32::from_f32_slice(&ctx, &x32);
        let mut gy = GpuVecF32::from_f32_slice(&ctx, &y32);

        // Atomic-free tree-reduce dot & norm2 parity vs CPU `tree_reduce_dot_f32`
        let dot_gpu_1 = gx.dot_f32(&gy).expect("dot_f32");
        let dot_gpu_2 = gx.dot_f32(&gy).expect("dot_f32 repeat");
        let dot_cpu = tree_reduce_dot_f32(&x32, &y32);
        assert_eq!(dot_gpu_1.to_bits(), dot_gpu_2.to_bits(), "GPU dot must be 0-ULP deterministic");
        let rel_dot = (dot_gpu_1 - dot_cpu).abs() / dot_cpu.abs().max(1.0);
        assert!(rel_dot < 1e-6, "GPU dot {dot_gpu_1} vs CPU tree_reduce {dot_cpu} (rel {rel_dot:e})");

        let nrm_gpu_1 = gx.norm2_sq_f32().expect("norm2_sq_f32");
        let nrm_gpu_2 = gx.norm2_sq_f32().expect("norm2_sq_f32 repeat");
        let nrm_cpu = tree_reduce_dot_f32(&x32, &x32);
        assert_eq!(nrm_gpu_1.to_bits(), nrm_gpu_2.to_bits(), "GPU norm2_sq must be 0-ULP deterministic");
        let rel_nrm = (nrm_gpu_1 - nrm_cpu).abs() / nrm_cpu.abs().max(1.0);
        assert!(rel_nrm < 1e-6, "GPU norm2_sq {nrm_gpu_1} vs CPU {nrm_cpu} (rel {rel_nrm:e})");

        // Vector trait BLAS-1 ops: axpy, aypx, scale, clone, zeros_like
        let mut y_ref = y32.clone();
        gy.axpy(-0.625, &gx);
        for (yr, &xr) in y_ref.iter_mut().zip(&x32) {
            *yr += -0.625f32 * xr;
        }
        gy.aypx(0.5, &gx);
        for (yr, &xr) in y_ref.iter_mut().zip(&x32) {
            *yr = xr + 0.5f32 * *yr;
        }
        gy.scale_real(1.5);
        for yr in &mut y_ref {
            *yr *= 1.5f32;
        }

        let got = gy.to_host_f32().expect("to_host_f32");
        let max_diff = got.iter().zip(&y_ref).map(|(a, b)| (a - b).abs()).fold(0.0f32, f32::max);
        assert!(max_diff < 1e-6, "GPU BLAS-1 max diff {max_diff:e}");

        let cloned = gy.clone();
        assert_eq!(cloned.to_host_f32().expect("clone readback"), got);
    }

    #[test]
    fn gpu_sell_c_sigma_f32_spmv_accuracy_and_determinism() {
        let ctx = GpuContext::shared().expect("GpuContext::shared");
        for a in [laplacian_3d(24), irregular_spd_csr(6_000)] {
            let n = a.nrows();
            let x = random_vector::<f64>(n, 3030);
            let mut y_ref = vec![0.0; n];
            a.spmv(&x, &mut y_ref);

            let gpu_a = GpuSellCSigmaF32::from_csr(&ctx, &a).expect("GpuSellCSigmaF32");
            let gx = GpuVecF32::from_f64_slice(&ctx, &x);
            let mut gy1 = GpuVecF32::zeros(&ctx, n);
            let mut gy2 = GpuVecF32::zeros(&ctx, n);
            gpu_a.spmv(&gx, &mut gy1).expect("spmv 1");
            gpu_a.spmv(&gx, &mut gy2).expect("spmv 2");

            let h1 = gy1.to_host_f32().expect("readback 1");
            let h2 = gy2.to_host_f32().expect("readback 2");
            for (u, v) in h1.iter().zip(&h2) {
                assert_eq!(u.to_bits(), v.to_bits(), "GPU SpMV must be 0-ULP deterministic");
            }
            let diff_sq: f64 = h1.iter().zip(&y_ref).map(|(&g, &r)| ((g as f64) - r).powi(2)).sum();
            let rel_err = diff_sq.sqrt() / vecops::norm2(&y_ref);
            assert!(rel_err < 1e-5, "GpuSellCSigmaF32 SpMV relative error {rel_err:e}");
        }

        // Finite guard check on GPU upload
        let bad = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, f64::NAN)]).expect("csr");
        assert!(GpuSellCSigmaF32::from_csr(&ctx, &bad).is_err());
    }

    #[test]
    fn gpu_amg_f32_drives_f64_pcg_fgmres_and_device_krylov_deterministically() {
        let ctx = GpuContext::shared().expect("GpuContext::shared");
        let a = laplacian_3d(24); // 13,824 unknowns, multi-level hierarchy
        let n = a.nrows();
        let b = random_vector::<f64>(n, 8080);

        let gpu_amg = GpuAmgF32::new_chebyshev(&ctx, &a, 3).expect("GpuAmgF32");
        let cpu_amg = AmgF32::new_chebyshev(&a, 3).expect("AmgF32");
        assert!(gpu_amg.stats().level_sizes.len() >= 3, "{:?}", gpu_amg.stats());

        let cg_opts = CgOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 100, ..Default::default() },
            record_history: true,
            ..Default::default()
        };

        // 1. Outer f64 PCG preconditioned by GPU Chebyshev V-cycle (two runs for bitwise determinism)
        let mut x_gpu_1 = vec![0.0; n];
        let info_1 = Cg::new(cg_opts.clone()).solve(&a, &gpu_amg, &b, &mut x_gpu_1);
        assert!(info_1.reason.is_converged(), "{info_1:?}");
        let rel_1 = true_rel_residual(&a, &x_gpu_1, &b);
        assert!(rel_1 < 1e-9, "GpuAmgF32 + f64 PCG relative residual {rel_1:e}");

        let mut x_gpu_2 = vec![0.0; n];
        let info_2 = Cg::new(cg_opts.clone()).solve(&a, &gpu_amg, &b, &mut x_gpu_2);
        assert_eq!(info_1.iterations, info_2.iterations);
        assert_f64_bits_eq(&info_1.history, &info_2.history, "GpuAmgF32 PCG history determinism");
        assert_f64_bits_eq(&x_gpu_1, &x_gpu_2, "GpuAmgF32 PCG solution 0-ULP determinism");

        // Compare iteration count with CPU AmgF32 Chebyshev
        let mut x_cpu = vec![0.0; n];
        let info_cpu = Cg::new(cg_opts).solve(&a, &cpu_amg, &b, &mut x_cpu);
        assert!(
            info_1.iterations.abs_diff(info_cpu.iterations) <= 1,
            "GPU AMG iterations {} vs CPU AMG iterations {}",
            info_1.iterations,
            info_cpu.iterations
        );

        // 2. Outer f64 FGMRES via `fgmres_vec` preconditioned by `GpuAmgF32`
        let fg_opts = GmresOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 100, ..Default::default() },
            restart: 25,
            side: PcSide::Right,
            flexible: true,
            orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
            ..Default::default()
        };
        let mut x_fg = CpuVec::zeros(n);
        let info_fg = fgmres_vec(&a, &gpu_amg, &CpuVec::new(b.clone()), &mut x_fg, &fg_opts);
        assert!(info_fg.reason.is_converged());
        assert!(true_rel_residual(&a, &x_fg, &b) < 1e-9);

        // 3. Fully device-resident `cg_vec<GpuVecF32>` and `lobpcg_vec<GpuVecF32>`
        let a_dev = GpuSellCSigmaF32::from_csr(&ctx, &a).expect("a_dev");
        let b_dev = GpuVecF32::from_f64_slice(&ctx, &b);
        let mut x_dev = GpuVecF32::zeros(&ctx, n);
        let dev_cg_opts = CgOptions {
            tol: Tolerances { rtol: 1e-5, max_it: 60, ..Default::default() },
            ..Default::default()
        };
        let dev_info = cg_vec(&a_dev, &gpu_amg, &b_dev, &mut x_dev, &dev_cg_opts);
        assert!(dev_info.reason.is_converged() && dev_info.iterations <= 20, "{dev_info:?}");
        let x_dev_host = x_dev.to_host();
        assert!(true_rel_residual(&a, &x_dev_host, &b) < 5e-5);

        // Small 3D Laplacian eigenproblem with tall vectors kept on the GPU
        let a_eig = laplacian_3d(10); // 1,000 unknowns
        let a_eig_gpu = GpuSellCSigmaF32::from_csr(&ctx, &a_eig).expect("a_eig_gpu");
        let amg_eig_gpu = GpuAmgF32::new_chebyshev(&ctx, &a_eig, 3).expect("amg_eig_gpu");
        let proto_gpu = GpuVecF32::zeros(&ctx, a_eig.nrows());
        let mut lopts = LobpcgOptions::new(2);
        lopts.tol = 2e-4;
        lopts.max_it = 60;
        let eig_gpu = lobpcg_vec::<GpuVecF32, _, GpuSellCSigmaF32, _>(
            &a_eig_gpu,
            None,
            Some(&amg_eig_gpu),
            &proto_gpu,
            None,
            &[],
            &lopts,
        )
        .expect("GPU lobpcg_vec");
        assert_eq!(eig_gpu.reason, LobpcgReason::Converged);
        // Exact lowest eigenvalue of 3D 7-point Laplacian on 10x10x10 interior grid:
        // 3 * (2 - 2 * cos(pi / 11))
        let exact_lam0 = 3.0 * (2.0 - 2.0 * (std::f64::consts::PI / 11.0).cos());
        assert!(
            (eig_gpu.values[0] - exact_lam0).abs() / exact_lam0 < 1e-4,
            "GPU LOBPCG lambda_0 = {} vs exact {}",
            eig_gpu.values[0],
            exact_lam0
        );
    }

    #[test]
    fn gpu_vec_f64_shader_f64_pcg_reaches_1e_10_and_zero_ulp_determinism() {
        let ctx = GpuContext::shared().expect("GpuContext::shared");
        if !ctx.supports_f64() {
            println!("Adapter '{}' lacks SHADER_F64; verifying clear error on GpuVecF64::zeros", ctx.info().name);
            assert!(GpuVecF64::zeros(&ctx, 16).is_err());
            return;
        }

        let a = laplacian_3d(24); // 13,824 unknowns
        let n = a.nrows();
        let b = random_vector::<f64>(n, 9090);
        let a_f64 = GpuSellCSigmaF64::from_csr(&ctx, &a).expect("GpuSellCSigmaF64");
        let gpu_amg = GpuAmgF32::new_chebyshev(&ctx, &a, 3).expect("GpuAmgF32");
        let b_dev = GpuVecF64::from_f64_slice(&ctx, &b).expect("GpuVecF64");

        // Verify dot, norm2_sq, and SpMV accuracy at f64 precision
        let mut y_ref = vec![0.0; n];
        a.spmv(&b, &mut y_ref);
        let mut y_dev = GpuVecF64::zeros(&ctx, n).expect("zeros");
        a_f64.spmv(&b_dev, &mut y_dev).expect("spmv_f64");
        let y_host = y_dev.to_host_f64().expect("readback f64");
        let spmv_rel = y_host.iter().zip(&y_ref).map(|(g, r)| (g - r).powi(2)).sum::<f64>().sqrt() / vecops::norm2(&y_ref);
        assert!(spmv_rel < 1e-14, "GpuSellCSigmaF64 relative error {spmv_rel:e}");

        // Device-resident f64 PCG with fast f32 GPU AMG V-cycle -> < 1e-10 true residual!
        let opts = CgOptions {
            tol: Tolerances { rtol: 1e-11, max_it: 80, ..Default::default() },
            norm: CgNorm::Unpreconditioned,
            record_history: true,
            ..Default::default()
        };
        let mut x1 = GpuVecF64::zeros(&ctx, n).expect("zeros x1");
        let info1 = cg_vec(&a_f64, &gpu_amg, &b_dev, &mut x1, &opts);
        assert!(info1.reason.is_converged() && info1.iterations <= 25, "{info1:?}");
        let x1_host = x1.to_host();
        let rel1 = true_rel_residual(&a, &x1_host, &b);
        assert!(rel1 < 1e-10, "GpuVecF64 + GpuAmgF32 residual {rel1:e}");

        let mut x2 = GpuVecF64::zeros(&ctx, n).expect("zeros x2");
        let info2 = cg_vec(&a_f64, &gpu_amg, &b_dev, &mut x2, &opts);
        assert_eq!(info1.iterations, info2.iterations);
        assert_f64_bits_eq(&info1.history, &info2.history, "GpuVecF64 PCG history determinism");
        assert_f64_bits_eq(&x1_host, &x2.to_host(), "GpuVecF64 PCG solution 0-ULP determinism");
    }
}
