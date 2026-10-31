#![allow(clippy::needless_range_loop)]

use femcore::solvers::ams::{Ams, AmsOptions};
use femcore::solvers::dense::Dense;
use femcore::solvers::ksp::{Cg, CgKind, CgOptions, Tolerances};
use femcore::solvers::linop::{Identity, Preconditioner};
use femcore::solvers::maxwell::{EdgeElementCube, nedelec_cube};
use femcore::solvers::pc::Jacobi;
use femcore::solvers::rng::random_vector;
use femcore::solvers::{C64, Csr};

use std::f64::consts::PI;

/// `K + βM`.
fn shifted(p: &EdgeElementCube, beta: f64) -> Csr<f64> {
    Csr::axpby(1.0, &p.stiffness, beta, &p.mass)
}

fn cg_iterations<A: femcore::solvers::LinearOperator<f64>, M: Preconditioner<f64>>(a: &A, m: &M, n: usize) -> (usize, bool) {
    let b = random_vector::<f64>(n, 11);
    let mut x = vec![0.0; n];
    let info =
        Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, max_it: 2000, ..Default::default() }, ..Default::default() })
            .solve(a, m, &b, &mut x);
    (info.iterations, info.reason.is_converged())
}

#[test]
fn pcg_with_ams_is_nearly_independent_of_the_mesh_while_jacobi_degrades() {
    let (mut ams_its, mut jac_its) = (Vec::new(), Vec::new());
    for n in [4usize, 8, 12] {
        let p = nedelec_cube(n);
        let a = shifted(&p, 1.0); // curl-curl dominated at these mesh sizes (β h² ≤ 0.06)
        let ne = a.nrows();
        let ams = Ams::new(&a, &p.gradient, &p.edge_vectors, &AmsOptions::default()).unwrap();
        let (it, ok) = cg_iterations(&a, &ams, ne);
        assert!(ok, "AMS-PCG did not converge for n = {n}");
        ams_its.push(it);
        let (it, ok) = cg_iterations(&a, &Jacobi::new(&a), ne);
        assert!(ok);
        jac_its.push(it);
    }
    eprintln!("AMS {ams_its:?} Jacobi {jac_its:?}");
    assert!(*ams_its.last().unwrap() <= 30, "{ams_its:?}");
    assert!(
        *ams_its.last().unwrap() as f64 <= 1.6 * ams_its[0] as f64 + 4.0,
        "AMS iterations grow with the mesh: {ams_its:?}"
    );
    assert!(*jac_its.last().unwrap() >= 4 * ams_its.last().unwrap(), "Jacobi {jac_its:?} vs AMS {ams_its:?}");
}

#[test]
fn the_gradient_correction_is_what_removes_the_null_space() {
    // Dropping G from the problem (an identity-like gradient of zero columns is not possible), so
    // compare with a plain scalar-diagonal smoother on a nearly singular system instead:
    // β = 1e-4 makes the gradient kernel almost exact; AMS must still converge in few iterations
    // while Jacobi needs hundreds.
    let p = nedelec_cube(8);
    let a = shifted(&p, 1e-4);
    let ne = a.nrows();
    let ams = Ams::new(&a, &p.gradient, &p.edge_vectors, &AmsOptions::default()).unwrap();
    let (it_ams, ok) = cg_iterations(&a, &ams, ne);
    assert!(ok, "AMS-PCG not converged");
    let (it_jac, _) = cg_iterations(&a, &Jacobi::new(&a), ne);
    let (it_none, _) = cg_iterations(&a, &Identity, ne);
    eprintln!("β = 1e-4: AMS {it_ams}, Jacobi {it_jac}, none {it_none}");
    assert!(it_ams <= 40 && it_jac >= 5 * it_ams && it_none >= 5 * it_ams);
}

#[test]
fn complex_symmetric_edge_matrix_with_cocg() {
    // A = K + (1 + 0.4i) M: complex symmetric, as from a lossy medium.
    let p = nedelec_cube(8);
    let a: Csr<C64> = Csr::axpby(C64::ONE, &p.stiffness.to_complex(), C64::new(1.0, 0.4), &p.mass.to_complex());
    let ne = a.nrows();
    let ams = Ams::new(&a, &p.gradient, &p.edge_vectors, &AmsOptions::default()).unwrap();
    let b = random_vector::<C64>(ne, 3);
    let mut x = vec![C64::ZERO; ne];
    let info = Cg::new(CgOptions {
        kind: CgKind::Symmetric,
        tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() },
        ..Default::default()
    })
    .solve(&a, &ams, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(info.iterations <= 40, "{}", info.iterations);
    // true residual
    let mut ax = vec![C64::ZERO; ne];
    a.spmv(&x, &mut ax);
    let r: f64 = ax.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
    let b0: f64 = b.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
    // the stopping test uses the preconditioned (natural) norm, so allow some slack in the true residual
    assert!(r / b0 < 1e-5, "{}", r / b0);
}

/// Dense reference: all eigenvalues of `K x = λ M x` through the Cholesky factor of `M`
/// (`C = L⁻¹ K L⁻ᵀ`, symmetric).
fn dense_eigenvalues(p: &EdgeElementCube) -> Vec<f64> {
    let n = p.stiffness.nrows();
    let (k, m) = (p.stiffness.to_dense(), p.mass.to_dense());
    let mut l = vec![0.0f64; n * n]; // row-major lower factor
    for i in 0..n {
        for j in 0..=i {
            let s: f64 = (0..j).map(|q| l[i * n + q] * l[j * n + q]).sum();
            l[i * n + j] = if i == j { (m[i * n + i] - s).sqrt() } else { (m[i * n + j] - s) / l[j * n + j] };
        }
    }
    // X = L⁻¹ K (forward substitution on columns), then C = X L⁻ᵀ = (L⁻¹ Xᵀ)ᵀ.
    let solve_lower = |b: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0f64; n * n];
        for c in 0..n {
            for i in 0..n {
                let s: f64 = (0..i).map(|q| l[i * n + q] * out[q * n + c]).sum();
                out[i * n + c] = (b[i * n + c] - s) / l[i * n + i];
            }
        }
        out
    };
    let x = solve_lower(&k);
    let xt: Vec<f64> = (0..n * n).map(|idx| x[(idx % n) * n + idx / n]).collect();
    let c = solve_lower(&xt);
    let dc = Dense::from_fn(n, n, |i, j| 0.5 * (c[i * n + j] + c[j * n + i]));
    let (mut v, _) = femcore::solvers::dense::hermitian_eig(&dc).unwrap();
    v.sort_by(|a, b| a.total_cmp(b));
    v
}

#[test]
fn generator_reproduces_the_unit_cube_cavity_spectrum() {
    // Exact: 2π² (triple), 3π² (double), 5π² (six: (1,2,0)-type), ...
    let mut first_errors = Vec::new();
    for n in [3usize, 5] {
        let p = nedelec_cube(n);
        let ev = dense_eigenvalues(&p);
        let nodes = p.coords.len();
        // the gradient null space: exactly `nodes` zero eigenvalues, then the cavity modes
        assert!(ev[nodes - 1].abs() < 1e-8 * ev[nodes + 2], "n={n}: {} vs {}", ev[nodes - 1], ev[nodes + 2]);
        let modes = &ev[nodes..nodes + 5];
        let exact = [2.0 * PI * PI, 2.0 * PI * PI, 2.0 * PI * PI, 3.0 * PI * PI, 3.0 * PI * PI];
        for (g, w) in modes.iter().zip(exact) {
            assert!((g - w).abs() < 0.12 * w, "n={n}: {g} vs {w}");
        }
        // the Kuhn mesh is not symmetric under axis permutations, so the triple splits a little
        let mean = (modes[0] + modes[1] + modes[2]) / 3.0;
        assert!(modes[2] - modes[0] < 0.15 * mean, "n={n}: triple spread {:?}", &modes[..3]);
        first_errors.push((mean - exact[0]).abs());
    }
    // second-order convergence (n 3 → 5): error ratio about (5/3)² = 2.8
    assert!(first_errors[0] > 1.8 * first_errors[1], "{first_errors:?}");
}

#[test]
fn jacobi_davidson_with_ams_finds_cavity_modes_without_a_factorisation() {
    use femcore::solvers::jd::{JdOptions, jacobi_davidson};
    // Modes near 2π² of the unit cube; the (K − τM) correction equation is preconditioned by
    // AMS built for the SPD operator K + τM. Eigenvalues 0 (gradients) are far from the target.
    let p = nedelec_cube(8);
    let tau = 2.0 * PI * PI + 0.3;
    let (k, m) = (p.stiffness.to_complex(), p.mass.to_complex());
    let spd: Csr<C64> = Csr::axpby(C64::ONE, &k, C64::new(tau, 0.0), &m);
    let ams = Ams::new(&spd, &p.gradient, &p.edge_vectors, &AmsOptions::default()).unwrap();
    let mut opts = JdOptions::new(3, C64::new(tau, 0.0), true);
    opts.max_it = 300;
    // Standard extraction: with the indefinite K − τM, harmonic extraction gets trapped in the
    // 0-eigenspace (gradients) because AMS(K + τM) amplifies exactly those components.
    opts.harmonic = false;
    opts.inner_its = 20;
    let r = jacobi_davidson(&k, Some(&m), &ams, None, &opts).unwrap();
    eprintln!(
        "JD+AMS: converged {} in {} outer / {} inner iterations, values {:?}",
        r.converged, r.iterations, r.inner_iterations, r.values
    );
    assert!(r.converged);
    for v in &r.values {
        assert!(v.im.abs() < 1e-8 && (v.re - 2.0 * PI * PI).abs() < 0.12 * 2.0 * PI * PI, "{v}");
        assert!(v.re > 1.0, "a gradient (zero) mode was returned: {v}");
    }
    assert!(r.residuals.iter().all(|&x| x < 1e-6), "{:?}", r.residuals);
}
