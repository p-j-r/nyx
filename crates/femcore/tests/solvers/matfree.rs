#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::dense::{Dense, hermitian_eig};
use femcore::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::solvers::ksp::{Cg, CgKind, CgOptions, Tolerances};
use femcore::solvers::linop::{Identity, LinearOperator};
use femcore::solvers::matfree::{Boundary, Gll, SpectralBox};
use femcore::solvers::pc::Jacobi;
use femcore::solvers::rng::{SplitMix64, random_vector};
use femcore::solvers::{C64, Csr, Scalar};
use std::f64::consts::PI;

#[test]
fn gll_rule_integrates_and_differentiates_polynomials_exactly() {
    for p in 1..=12 {
        let g = Gll::new(p);
        assert_eq!(g.nodes.len(), p + 1);
        assert!((g.weights.iter().sum::<f64>() - 2.0).abs() < 1e-13);
        for i in 0..=p {
            assert!((g.nodes[i] + g.nodes[p - i]).abs() < 1e-13, "symmetry p={p}");
        }
        // GLL quadrature is exact up to degree 2p−1
        for k in 0..2 * p {
            let q: f64 = g.nodes.iter().zip(&g.weights).map(|(x, w)| w * x.powi(k as i32)).sum();
            let exact = if k % 2 == 0 { 2.0 / (k as f64 + 1.0) } else { 0.0 };
            assert!((q - exact).abs() < 1e-12, "p={p} k={k}: {q} vs {exact}");
        }
        // D differentiates every polynomial of degree ≤ p
        for k in 0..=p {
            for i in 0..=p {
                let d: f64 = (0..=p).map(|j| g.d[i * (p + 1) + j] * g.nodes[j].powi(k as i32)).sum();
                let exact = if k == 0 { 0.0 } else { k as f64 * g.nodes[i].powi(k as i32 - 1) };
                assert!((d - exact).abs() < 1e-9 * (1.0 + exact.abs()), "p={p} k={k} i={i}: {d} vs {exact}");
            }
        }
    }
}

fn random_box<T: Scalar>(p: usize, nel: &[usize], boundary: Boundary, seed: u64, complex: bool) -> SpectralBox<T> {
    let mut rng = SplitMix64::new(seed);
    let widths: Vec<Vec<f64>> = nel.iter().map(|&n| (0..n).map(|_| 0.5 + rng.next_signed().abs()).collect()).collect();
    let total: usize = nel.iter().product();
    let mut coef = |scale: f64| -> Vec<T> {
        (0..total)
            .map(|_| {
                T::from_parts(
                    scale * (1.0 + rng.next_signed().abs()),
                    if complex { 0.3 * rng.next_signed() } else { 0.0 },
                )
            })
            .collect()
    };
    let (kappa, mu) = (coef(1.0), coef(0.7));
    SpectralBox::new(p, &widths, kappa, mu, boundary)
}

fn check_matches_assembled<T: Scalar>(complex: bool) {
    for (dim_nel, ps) in
        [(vec![5usize], vec![1usize, 3, 6]), (vec![3, 2], vec![1, 2, 5]), (vec![2, 2, 3], vec![1, 2, 4])]
    {
        for p in ps {
            for boundary in [Boundary::Dirichlet, Boundary::Natural] {
                let op = random_box::<T>(p, &dim_nel, boundary, 7 + p as u64, complex);
                let a = op.assemble();
                assert_eq!(a.nrows(), op.len());
                let x = random_vector::<T>(op.len(), 11);
                let (mut y1, mut y2) = (vec![T::ZERO; op.len()], vec![T::ZERO; op.len()]);
                op.apply(&x, &mut y1);
                a.spmv(&x, &mut y2);
                let err = y1.iter().zip(&y2).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt() / norm(&y2);
                assert!(err < 1e-13, "dim {} p={p} {boundary:?}: {err}", dim_nel.len());
                // diagonal without assembly
                let d = op.diagonal();
                let da = a.diagonal();
                let derr = d.iter().zip(&da).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt() / norm(&da);
                assert!(derr < 1e-13, "diagonal p={p}: {derr}");
            }
        }
    }
}

#[test]
fn matrix_free_apply_and_diagonal_equal_the_assembled_matrix() {
    check_matches_assembled::<f64>(false);
    check_matches_assembled::<C64>(true);
}

#[test]
fn spectral_accuracy_of_the_laplace_eigenvalues() {
    // −Δ on (0,π)² with Dirichlet conditions: λ = k² + l² = 2, 5, 5, 8, 10, 10, …
    let exact = [2.0, 5.0, 5.0, 8.0, 10.0, 10.0];
    let mut prev = f64::INFINITY;
    for p in [4usize, 6, 8] {
        let h = PI / 2.0;
        let k = SpectralBox::<f64>::uniform(p, &[2, 2], &[h, h], 1.0, 0.0, Boundary::Dirichlet);
        let m = SpectralBox::<f64>::uniform(p, &[2, 2], &[h, h], 0.0, 1.0, Boundary::Dirichlet);
        let n = k.len();
        let interior: Vec<usize> = (0..n).filter(|&g| !k.on_boundary(g)).collect();
        let md = m.diagonal();
        let ni = interior.len();
        // dense  M^{-1/2} K M^{-1/2}  on the interior, built purely from operator applications
        let mut dense = Dense::<f64>::zeros(ni, ni);
        let mut e = vec![0.0; n];
        let mut y = vec![0.0; n];
        for (cj, &gj) in interior.iter().enumerate() {
            e.fill(0.0);
            e[gj] = 1.0;
            k.apply(&e, &mut y);
            for (ci, &gi) in interior.iter().enumerate() {
                dense[(ci, cj)] = y[gi] / (md[gi] * md[gj]).sqrt();
            }
        }
        let (w, _) = hermitian_eig(&dense).unwrap();
        let err = (0..6).map(|i| (w[i] - exact[i]).abs()).fold(0.0, f64::max);
        println!("p = {p}: max error of the 6 lowest eigenvalues {err:.2e}");
        assert!(err < prev, "no improvement with p");
        prev = err;
    }
    assert!(prev < 1e-7, "p=8: {prev}");
}

#[test]
fn poisson_solve_converges_spectrally_with_the_degree() {
    // −Δu = f on (0,π)², u = sin(x) sin(2y);  K u_h = M f at the nodes (GLL collocation)
    let mut errors = Vec::new();
    for p in [2usize, 4, 6, 8] {
        let nel = 3;
        let h = PI / nel as f64;
        let k = SpectralBox::<f64>::uniform(p, &[nel, nel], &[h, h], 1.0, 0.0, Boundary::Dirichlet);
        let m = SpectralBox::<f64>::uniform(p, &[nel, nel], &[h, h], 0.0, 1.0, Boundary::Dirichlet);
        let n = k.len();
        let nx = k.nodes_per_axis()[0];
        let md = m.diagonal();
        let u_exact = |g: usize| (k.coord(0, g % nx)).sin() * (2.0 * k.coord(1, g / nx)).sin();
        let b: Vec<f64> = (0..n).map(|g| if k.on_boundary(g) { 0.0 } else { 5.0 * u_exact(g) * md[g] }).collect();
        let diag = k.diagonal();
        let pre = Jacobi::from_diagonal(diag);
        let mut x = vec![0.0; n];
        let info = Cg::new(CgOptions {
            tol: Tolerances { rtol: 1e-13, max_it: 4000, ..Default::default() },
            ..Default::default()
        })
        .solve(&k, &pre, &b, &mut x);
        assert!(info.reason.is_converged(), "p={p}: {:?}", info.reason);
        let err = (0..n).map(|g| (x[g] - u_exact(g)).abs()).fold(0.0, f64::max);
        println!("p = {p}: nodal max error {err:.2e} ({} CG iterations)", info.iterations);
        errors.push(err);
    }
    assert!(errors.windows(2).all(|w| w[1] < w[0]), "{errors:?}");
    assert!(errors[3] < 1e-8, "{errors:?}");
}

#[test]
fn complex_helmholtz_with_absorption_matches_a_direct_solve_of_the_assembled_matrix() {
    let (p, nel) = (4usize, 4usize);
    let h = 1.0 / nel as f64;
    // −Δ − k² (1 + iα) u:  κ = 1, μ = −k²(1 + iα)
    let mu = C64::new(-20.0, -4.0);
    let op = SpectralBox::<C64>::uniform(p, &[nel, nel], &[h, h], C64::ONE, mu, Boundary::Dirichlet);
    let a = op.assemble();
    let n = op.len();
    let b: Vec<C64> = random_vector::<C64>(n, 3)
        .into_iter()
        .enumerate()
        .map(|(g, v)| if op.on_boundary(g) { C64::ZERO } else { v })
        .collect();
    let mut f = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
    f.set_shift(C64::ZERO).unwrap();
    let mut exact = b.clone();
    f.solve_block_in_place(&mut exact, 1);
    let jac = Jacobi::from_diagonal(op.diagonal());
    let mut x = vec![C64::ZERO; n];
    let info = Cg::new(CgOptions {
        tol: Tolerances { rtol: 1e-11, max_it: 3000, ..Default::default() },
        kind: CgKind::Symmetric,
        ..Default::default()
    })
    .solve(&op, &jac, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    let err = x.iter().zip(&exact).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&exact);
    assert!(err < 1e-8, "{err}");
    let _ = (Identity, Csr::<f64>::identity(1));
}

// ---- mapped (curved) geometry ----------------------------------------------------------------------------

/// Node coordinates of a box mesh warped by a smooth map that keeps the boundary of the box fixed.
fn warped_coords(base: &SpectralBox<f64>, amp: f64, len: f64) -> Vec<Vec<f64>> {
    let nn = base.nodes_per_axis().to_vec();
    let dim = nn.len();
    let n = base.len();
    let mut out = vec![vec![0.0; n]; dim];
    for g in 0..n {
        let mut r = g;
        let mut x = [0.0; 3];
        for d in 0..dim {
            x[d] = base.coord(d, r % nn[d]);
            r /= nn[d];
        }
        for d in 0..dim {
            // x_d + amp·(1 + d/2) Π_e sin(π x_e/len)  (vanishes on every face of the box)
            let bump: f64 = (0..dim).map(|e| (PI * x[e] / len).sin()).product();
            out[d][g] = x[d] + amp * bump * (1.0 + 0.5 * d as f64);
        }
    }
    out
}

#[test]
fn mapped_geometry_reduces_to_the_box_kernel_for_affine_maps() {
    for (nel, p) in [(vec![4usize], 3usize), (vec![3, 3], 4), (vec![2, 2, 2], 3)] {
        let h = 0.7;
        let hs = vec![h; nel.len()];
        let bx = SpectralBox::<f64>::uniform(p, &nel, &hs, 1.3, -2.0, Boundary::Dirichlet);
        let coords: Vec<Vec<f64>> = (0..nel.len())
            .map(|d| {
                let nn = bx.nodes_per_axis().to_vec();
                (0..bx.len()).map(|g| bx.coord(d, (g / nn[..d].iter().product::<usize>()) % nn[d])).collect()
            })
            .collect();
        let total: usize = nel.iter().product();
        let mp = SpectralBox::<f64>::mapped(p, &nel, &coords, vec![1.3; total], vec![-2.0; total], Boundary::Dirichlet);
        let x = random_vector::<f64>(bx.len(), 5);
        let (mut y1, mut y2) = (vec![0.0; bx.len()], vec![0.0; bx.len()]);
        bx.apply(&x, &mut y1);
        mp.apply(&x, &mut y2);
        let err = y1.iter().zip(&y2).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt() / norm(&y1);
        assert!(err < 1e-12, "dim {} p {p}: {err}", nel.len());
    }
}

#[test]
fn mapped_diagonal_matches_the_assembled_matrix_on_a_curved_mesh() {
    for (nel, p) in [(vec![3usize, 3], 3usize), (vec![2, 2, 2], 3)] {
        let hs = vec![1.0; nel.len()];
        let base = SpectralBox::<f64>::uniform(p, &nel, &hs, 1.0, 0.0, Boundary::Natural);
        let coords = warped_coords(&base, 0.12, nel[0] as f64);
        let total: usize = nel.iter().product();
        let op = SpectralBox::<f64>::mapped(p, &nel, &coords, vec![1.0; total], vec![0.5; total], Boundary::Natural);
        let a = op.assemble();
        assert!(a.is_symmetric(1e-12));
        let (d, da) = (op.diagonal(), a.diagonal());
        let err = d.iter().zip(&da).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt() / norm(&da);
        assert!(err < 1e-12, "diagonal {err}");
        // constants are in the kernel of the stiffness part (mu = 0)
        let k = SpectralBox::<f64>::mapped(p, &nel, &coords, vec![1.0; total], vec![0.0; total], Boundary::Natural);
        let ones = vec![1.0; k.len()];
        let mut y = vec![0.0; k.len()];
        k.apply(&ones, &mut y);
        assert!(y.iter().all(|v| v.abs() < 1e-11), "K·1 = {:e}", y.iter().fold(0.0f64, |m, v| m.max(v.abs())));
    }
}

#[test]
fn poisson_on_a_curved_mesh_converges_spectrally() {
    // same exact solution as above, but the interior of the mesh is warped: the geometric tensor G has
    // off-diagonal entries everywhere
    let mut errors = Vec::new();
    for p in [3usize, 5, 7] {
        let nel = 3;
        let base = SpectralBox::<f64>::uniform(p, &[nel, nel], &[PI / nel as f64; 2], 1.0, 0.0, Boundary::Dirichlet);
        let coords = warped_coords(&base, 0.25, PI);
        let n = base.len();
        let ones = vec![1.0; nel * nel];
        let k = SpectralBox::<f64>::mapped(
            p,
            &[nel, nel],
            &coords,
            ones.clone(),
            vec![0.0; nel * nel],
            Boundary::Dirichlet,
        );
        let m = SpectralBox::<f64>::mapped(p, &[nel, nel], &coords, vec![0.0; nel * nel], ones, Boundary::Dirichlet);
        let md = m.diagonal();
        let u_exact = |g: usize| coords[0][g].sin() * (2.0 * coords[1][g]).sin();
        let b: Vec<f64> = (0..n).map(|g| if k.on_boundary(g) { 0.0 } else { 5.0 * u_exact(g) * md[g] }).collect();
        let pre = Jacobi::from_diagonal(k.diagonal());
        let mut x = vec![0.0; n];
        let info = Cg::new(CgOptions {
            tol: Tolerances { rtol: 1e-13, max_it: 5000, ..Default::default() },
            ..Default::default()
        })
        .solve(&k, &pre, &b, &mut x);
        assert!(info.reason.is_converged(), "p={p}");
        let err = (0..n).map(|g| (x[g] - u_exact(g)).abs()).fold(0.0, f64::max);
        println!("curved mesh p = {p}: nodal max error {err:.2e}");
        errors.push(err);
    }
    assert!(errors.windows(2).all(|w| w[1] < w[0]), "{errors:?}");
    assert!(errors[2] < 5e-5, "{errors:?}");
}
