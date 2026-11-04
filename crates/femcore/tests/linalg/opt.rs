//! Optimisers: Rosenbrock and quadratic problems for L-BFGS / nonlinear CG, MMA on a problem with a known KKT
//! solution, adjoint gradients against finite differences, and a 2-D heat-conduction topology optimisation
//! ("volume to point") with MMA and adjoint sensitivities.
#![allow(clippy::needless_range_loop)]

use femcore::linalg::Csr;
use femcore::linalg::csr::laplacian_1d;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::opt::{MinOptions, MmaEval, MmaOptions, StopReason, adjoint_gradient, lbfgs, mma, nonlinear_cg};

fn rosenbrock(x: &[f64], g: &mut [f64]) -> f64 {
    let n = x.len();
    g.fill(0.0);
    let mut f = 0.0;
    for i in 0..n - 1 {
        let (a, b) = (1.0 - x[i], x[i + 1] - x[i] * x[i]);
        f += a * a + 100.0 * b * b;
        g[i] += -2.0 * a - 400.0 * x[i] * b;
        g[i + 1] += 200.0 * b;
    }
    f
}

#[test]
fn lbfgs_and_cg_minimise_rosenbrock() {
    let x0: Vec<f64> = (0..10).map(|i| if i % 2 == 0 { -1.2 } else { 1.0 }).collect();
    let l = lbfgs(&mut rosenbrock, &x0, &MinOptions { gtol: 1e-9, max_it: 500, ..Default::default() });
    assert!(l.reason.is_converged(), "{:?} after {} iterations", l.reason, l.iterations);
    assert!(l.x.iter().all(|v| (v - 1.0).abs() < 1e-6), "{:?}", l.x);
    let c = nonlinear_cg(&mut rosenbrock, &x0, &MinOptions { gtol: 1e-8, max_it: 5000, c2: 0.1, ..Default::default() });
    assert!(c.reason.is_converged(), "{:?} after {} iterations", c.reason, c.iterations);
    assert!(c.x.iter().all(|v| (v - 1.0).abs() < 1e-5), "{:?}", c.x);
    eprintln!(
        "Rosenbrock(10): L-BFGS {} its / {} evals, CG {} its / {} evals",
        l.iterations, l.evaluations, c.iterations, c.evaluations
    );
    assert!(l.iterations < c.iterations);
    assert!(l.history.windows(2).all(|w| w[1] <= w[0] + 1e-12), "line search must give monotone decrease");
}

#[test]
fn lbfgs_solves_an_ill_conditioned_quadratic_in_few_iterations() {
    // f = ½ xᵀ L x − bᵀ x with the 1-D Laplacian (cond ≈ 4n²/π²): n = 100
    let n = 100;
    let l = laplacian_1d(n);
    let b: Vec<f64> = (0..n).map(|i| ((i * 13) % 7) as f64 - 3.0).collect();
    let mut fg = |x: &[f64], g: &mut [f64]| {
        l.spmv(x, g);
        let f = 0.5 * x.iter().zip(g.iter()).map(|(a, c)| a * c).sum::<f64>()
            - x.iter().zip(&b).map(|(a, c)| a * c).sum::<f64>();
        for (gi, bi) in g.iter_mut().zip(&b) {
            *gi -= bi;
        }
        f
    };
    let r = lbfgs(&mut fg, &vec![0.0; n], &MinOptions { gtol: 1e-10, max_it: 1000, memory: 20, ..Default::default() });
    assert!(r.reason.is_converged(), "{:?}", r.reason);
    let mut res = vec![0.0; n];
    l.spmv(&r.x, &mut res);
    assert!(res.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-6));
}

#[test]
fn mma_reproduces_the_kkt_solution_of_a_separable_problem() {
    // min Σ c_j / x_j  s.t.  Σ x_j ≤ V, 0.01 ≤ x ≤ 10  ⇒  x_j = V sqrt(c_j) / Σ sqrt(c_k)
    let c = [1.0, 4.0, 9.0, 16.0, 25.0];
    let v = 10.0;
    let mut ev = |x: &[f64]| MmaEval {
        f0: x.iter().zip(&c).map(|(xi, ci)| ci / xi).sum(),
        g0: x.iter().zip(&c).map(|(xi, ci)| -ci / (xi * xi)).collect(),
        c: vec![x.iter().sum::<f64>() - v],
        gc: vec![vec![1.0; 5]],
    };
    let opts = MmaOptions { max_it: 200, xtol: 1e-9, ..MmaOptions::new(vec![0.01; 5], vec![10.0; 5]) };
    let r = mma(&mut ev, &[2.0; 5], &opts);
    assert!(r.reason.is_converged() || r.iterations > 50, "{:?}", r.reason);
    let s: f64 = c.iter().map(|v| v.sqrt()).sum();
    for (xi, ci) in r.x.iter().zip(&c) {
        let want = v * ci.sqrt() / s;
        assert!((xi - want).abs() < 1e-4 * want.max(1.0), "{xi} vs {want}");
    }
    assert!(r.constraints[0] < 1e-6);
}

#[test]
fn mma_handles_an_initially_infeasible_design_and_two_constraints() {
    // min x0 + x1 + x2  s.t.  1 − x0 x1 ≤ 0 ... use convex form: 1/x0 + 1/x1 ≤ 2 and x2 ≥ 0.5 (i.e. 0.5 − x2 ≤ 0)
    let mut ev = |x: &[f64]| MmaEval {
        f0: x[0] + x[1] + x[2],
        g0: vec![1.0; 3],
        c: vec![1.0 / x[0] + 1.0 / x[1] - 2.0, 0.5 - x[2]],
        gc: vec![vec![-1.0 / (x[0] * x[0]), -1.0 / (x[1] * x[1]), 0.0], vec![0.0, 0.0, -1.0]],
    };
    let opts = MmaOptions { max_it: 300, xtol: 1e-10, ..MmaOptions::new(vec![0.1; 3], vec![5.0; 3]) };
    let r = mma(&mut ev, &[0.3, 0.3, 3.0], &opts); // start violates the first constraint
    assert!((r.x[0] - 1.0).abs() < 1e-3 && (r.x[1] - 1.0).abs() < 1e-3 && (r.x[2] - 0.5).abs() < 1e-3, "{:?}", r.x);
    assert!(r.constraints.iter().all(|&ci| ci < 1e-5), "{:?}", r.constraints);
}

/// 2-D heat conduction on an `nx × ny` grid of bilinear elements with conductivity `ρ_e^p + ε`, heat sink (Dirichlet 0)
/// at the middle of the left edge, uniform heat load.
struct Heat {
    nx: usize,
    ny: usize,
    ke: [[f64; 4]; 4],
    penal: f64,
    eps: f64,
}

impl Heat {
    fn new(nx: usize, ny: usize) -> Self {
        // bilinear element conductivity matrix on a unit square
        let a = [[4.0, -1.0, -2.0, -1.0], [-1.0, 4.0, -1.0, -2.0], [-2.0, -1.0, 4.0, -1.0], [-1.0, -2.0, -1.0, 4.0]];
        let mut ke = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                ke[i][j] = a[i][j] / 6.0;
            }
        }
        Self { nx, ny, ke, penal: 3.0, eps: 1e-3 }
    }
    fn node(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 1) + i
    }
    fn elem_nodes(&self, e: usize) -> [usize; 4] {
        let (i, j) = (e % self.nx, e / self.nx);
        [self.node(i, j), self.node(i + 1, j), self.node(i + 1, j + 1), self.node(i, j + 1)]
    }
    /// Reduced system (sink node removed): returns `K`, the load `f`, and the map from nodes to dofs.
    fn assemble(&self, rho: &[f64]) -> (Csr<f64>, Vec<f64>, Vec<Option<usize>>) {
        let nn = (self.nx + 1) * (self.ny + 1);
        let sink = self.node(0, self.ny / 2);
        let mut map = vec![None; nn];
        let mut k = 0;
        for i in 0..nn {
            if i != sink {
                map[i] = Some(k);
                k += 1;
            }
        }
        let mut trip = Vec::new();
        for e in 0..self.nx * self.ny {
            let c = self.eps + (1.0 - self.eps) * rho[e].powf(self.penal);
            let nodes = self.elem_nodes(e);
            for a in 0..4 {
                for b in 0..4 {
                    if let (Some(ia), Some(ib)) = (map[nodes[a]], map[nodes[b]]) {
                        trip.push((ia, ib, c * self.ke[a][b]));
                    }
                }
            }
        }
        let kmat = Csr::from_triplets(k, k, &trip).unwrap();
        let f: Vec<f64> = vec![1.0 / nn as f64; k]; // total load 1
        (kmat, f, map)
    }
    /// Compliance `fᵀT` and its gradient with respect to the element densities (adjoint = self-adjoint here).
    fn compliance(&self, rho: &[f64]) -> (f64, Vec<f64>) {
        let (kmat, f, map) = self.assemble(rho);
        let ne = self.nx * self.ny;
        let nodes_of = |e: usize| self.elem_nodes(e);
        let dk_u = |e: usize, u: &[f64]| -> Vec<f64> {
            let dc = self.penal * rho[e].powf(self.penal - 1.0) * (1.0 - self.eps);
            let mut out = vec![0.0; u.len()];
            let nodes = nodes_of(e);
            for a in 0..4 {
                for b in 0..4 {
                    if let (Some(ia), Some(ib)) = (map[nodes[a]], map[nodes[b]]) {
                        out[ia] += dc * self.ke[a][b] * u[ib];
                    }
                }
            }
            out
        };
        let (j, g, _) = adjoint_gradient(&kmat, &f, &f, ne, dk_u, |_| None, Definiteness::PositiveDefinite).unwrap();
        (j, g)
    }
}

#[test]
fn adjoint_gradient_matches_finite_differences() {
    let h = Heat::new(8, 6);
    let ne = h.nx * h.ny;
    let rho: Vec<f64> = (0..ne).map(|e| 0.3 + 0.5 * ((e * 7) % 5) as f64 / 5.0).collect();
    let (_, g) = h.compliance(&rho);
    for e in [0usize, 5, 17, 30, 47] {
        let d = 1e-6;
        let (mut rp, mut rm) = (rho.clone(), rho.clone());
        rp[e] += d;
        rm[e] -= d;
        let fd = (h.compliance(&rp).0 - h.compliance(&rm).0) / (2.0 * d);
        assert!(
            (g[e] - fd).abs() < 1e-5 * fd.abs().max(1e-8),
            "element {e}: adjoint {} vs finite difference {fd}",
            g[e]
        );
    }
}

#[test]
fn mma_topology_optimisation_of_a_heat_conduction_design() {
    let h = Heat::new(24, 24);
    let ne = h.nx * h.ny;
    let volfrac = 0.4;
    let mut ev = |x: &[f64]| {
        let (c, g) = h.compliance(x);
        MmaEval {
            f0: c,
            g0: g,
            c: vec![x.iter().sum::<f64>() / ne as f64 - volfrac],
            gc: vec![vec![1.0 / ne as f64; ne]],
        }
    };
    let x0 = vec![volfrac; ne];
    let uniform = h.compliance(&x0).0;
    let opts = MmaOptions { max_it: 120, xtol: 1e-4, ..MmaOptions::new(vec![1e-3; ne], vec![1.0; ne]) };
    let r = mma(&mut ev, &x0, &opts);
    let vol = r.x.iter().sum::<f64>() / ne as f64;
    let grey = r.x.iter().filter(|&&v| v > 0.1 && v < 0.9).count() as f64 / ne as f64;
    eprintln!(
        "topology: compliance {uniform:.4e} → {:.4e} ({:.1}× lower), volume {vol:.4}, grey fraction {grey:.2}, {} iterations",
        r.f0,
        uniform / r.f0,
        r.iterations
    );
    assert!(vol <= volfrac + 1e-4, "volume {vol}");
    assert!(r.f0 < 0.2 * uniform, "compliance should drop by a large factor: {uniform} → {}", r.f0);
    assert!(grey < 0.4, "design should be mostly solid/void (grey fraction {grey})");
    // the history is (nearly) monotone once the constraint is satisfied
    let tail = &r.history[r.history.len() / 2..];
    assert!(tail.windows(2).all(|w| w[1] <= w[0] * 1.05), "{tail:?}");
    let _ = StopReason::MaxIterations;
}

#[test]
fn inverse_design_of_a_potential_to_hit_a_target_eigenvalue_with_eigenvalue_derivatives() {
    // −u'' + V u = λ u on a grid (h = 1): design V (n parameters) so that the lowest eigenvalue equals λ*,
    // with a small smoothness penalty. d λ₁/dV_i = x_i² / ‖x‖² (Hermitian, N3), L-BFGS (O1).
    use femcore::linalg::dense::{Dense, hermitian_eig};
    let n = 40;
    let target = 0.5;
    let reg = 1e-6;
    let l = laplacian_1d(n).to_dense();
    let mut fg = |v: &[f64], g: &mut [f64]| -> f64 {
        let a = Dense::from_fn(n, n, |i, j| l[i * n + j] + if i == j { v[i] } else { 0.0 });
        let (w, q) = hermitian_eig(&a).unwrap();
        let k = (0..n).min_by(|&i, &j| w[i].total_cmp(&w[j])).unwrap();
        let x = q.col(k);
        let xn: f64 = x.iter().map(|t| t * t).sum();
        let mut f = (w[k] - target).powi(2);
        for i in 0..n {
            g[i] = 2.0 * (w[k] - target) * x[i] * x[i] / xn;
        }
        for i in 0..n - 1 {
            let d = v[i + 1] - v[i];
            f += reg * d * d;
            g[i] -= 2.0 * reg * d;
            g[i + 1] += 2.0 * reg * d;
        }
        f
    };
    // gradient check against finite differences at a random point
    let v0: Vec<f64> = (0..n).map(|i| 0.1 + 0.05 * ((i * 11) % 7) as f64).collect();
    let mut g0 = vec![0.0; n];
    fg(&v0, &mut g0);
    let mut scratch = vec![0.0; n];
    for i in [3usize, 17, 31] {
        let (mut vp, mut vm) = (v0.clone(), v0.clone());
        vp[i] += 1e-6;
        vm[i] -= 1e-6;
        let fd = (fg(&vp, &mut scratch) - fg(&vm, &mut scratch)) / 2e-6;
        assert!((g0[i] - fd).abs() < 1e-6 * fd.abs().max(1e-6), "{i}: {} vs {fd}", g0[i]);
    }
    let r = lbfgs(&mut fg, &v0, &MinOptions { gtol: 1e-9, max_it: 200, ..Default::default() });
    // the smoothness penalty keeps the optimum slightly away from λ₁ = λ*
    assert!(
        r.history.last().unwrap() < &(1e-6 * r.history[0]),
        "objective {:?} after {} iterations ({:?})",
        r.history.last(),
        r.iterations,
        r.reason
    );
    let a = Dense::from_fn(n, n, |i, j| l[i * n + j] + if i == j { r.x[i] } else { 0.0 });
    let (w, _) = hermitian_eig(&a).unwrap();
    let lmin = w.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!((lmin - target).abs() < 1e-3, "{lmin}");
}
