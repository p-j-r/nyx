//! **Matrix-free high-order operators** by sum factorisation (port item M2): the spectral-element
//! (Gauss–Lobatto, "SEM") stiffness and mass operators on a structured box mesh, applied without ever
//! assembling a matrix.
//!
//! # Why (and when it does *not* pay)
//! On a **curved or skewed** hexahedral element the geometry couples all reference directions, so the assembled
//! element matrix is dense: `(p+1)^{2d}` entries (`p = 6`, 3-D: 117 649 numbers for 343 unknowns) and the sparse
//! matrix–vector product is bound by memory bandwidth. A tensor-product element can instead be applied in
//! `O(d²(p+1)^{d+1})` operations reading only the vector: **sum factorisation** (Orszag 1980), the technique behind
//! deal.II `MatrixFree`, libCEED, MFEM partial assembly and Nek5000.
//! On an **axis-aligned box with GLL collocation** (the fast path of [`SpectralBox::new`]) this is *not* so: the
//! directions decouple, the element matrices have only about `d(p+1)` non-zeros per row, and a plain CSR matrix is
//! both smaller and faster (measured: `examples/matfree_bench.rs`, `docs/LEARNINGS.md`). The box path exists as the
//! simple, independently checkable reference and for coefficient-heterogeneous box meshes without assembly.
//!
//! # What is here
//! * [`Gll`]: Gauss–Lobatto–Legendre nodes, weights and the collocation derivative matrix.
//! * [`SpectralBox::mapped`]: the same on curved/skewed elements, geometry stored as `d(d+1)/2 + 1` numbers per
//!   quadrature point (memory `7·n` in 3-D, independent of `p`).
//! * [`SpectralBox`]: the operator `A = Σ_e (κ_e K_e + μ_e M_e)` on `nex × ney × nez` boxes of arbitrary
//!   (non-uniform) widths with the same degree `p` in all directions (dimension 1, 2 or 3), per-element
//!   coefficients `κ_e, μ_e` of type `T` (real or complex: Helmholtz with absorption), homogeneous Dirichlet
//!   conditions on the whole boundary or none. `K_e = Σ_d (2/h_d)² detJ · D_dᵀ W D_d` (the 1-D matrices
//!   `DᵀWD` along `d`, the weights `W` along the others) and the *lumped* mass `M_e = detJ · W`, which is
//!   diagonal because quadrature and nodes coincide. It implements [`LinearOperator`].
//! * [`SpectralBox::diagonal`]: the exact diagonal without assembly (Jacobi, Chebyshev smoothing, and the
//!   diagonal weights of the generalised eigenproblem).
//! * [`SpectralBox::assemble`]: an independent route to the same operator as a [`Csr`] through the Kronecker
//!   structure of the element matrices (validation, and the comparison in `examples/matfree_bench.rs`).
//! * Elements are processed in `z`-layers; with the `parallel` feature even layers, then odd layers, run
//!   concurrently on disjoint slices of the output (no atomics, no unsafe code), and the result is bit-identical
//!   for every thread count because each output entry receives its contributions in a fixed order.
//!
//! # Limits
//! Structured (tensor-product) topology only: the mesh is `nex × ney × nez` elements, geometry may be curved
//! ([`SpectralBox::mapped`]) but unstructured hexahedral meshes would need an element-to-node table in the gather/scatter.
//! Collocation (GLL) quadrature with lumped mass; over-integration (Gauss points, exact mass) would add one more
//! 1-D interpolation per axis. The inner kernels are plain loops (about 1 Gflop/s on the test VM): no SIMD
//! specialisation per `p`, which is the main remaining speed-up. No `p`-multigrid or other preconditioner beyond the exact diagonal. Mass lumping is a modelling choice (exact for GLL collocation, `O(h^{p+1})`
//! consistent with the full mass matrix for smooth data).

#![allow(clippy::needless_range_loop)]

use crate::solvers::csr::Csr;
use crate::solvers::linop::LinearOperator;
use crate::solvers::scalar::Scalar;

/// Gauss–Lobatto–Legendre rule of degree `p` (`p + 1` nodes on `[−1, 1]`).
#[derive(Debug, Clone)]
pub struct Gll {
    pub p: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
    /// `d[i·(p+1) + j] = ℓ_j′(x_i)`, row-major.
    pub d: Vec<f64>,
}

fn legendre(p: usize, x: f64) -> (f64, f64) {
    // P_p(x) and P_p'(x) by the three-term recurrence
    if p == 0 {
        return (1.0, 0.0);
    }
    let (mut p0, mut p1) = (1.0, x);
    for k in 2..=p {
        let kf = k as f64;
        let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
        p0 = p1;
        p1 = p2;
    }
    let dp = if (x.abs() - 1.0).abs() < 1e-14 {
        // P'(±1) = ±^{p+1} p(p+1)/2
        let s = if x > 0.0 || p % 2 == 1 { 1.0 } else { -1.0 };
        s * (p * (p + 1)) as f64 / 2.0
    } else {
        p as f64 * (x * p1 - p0) / (x * x - 1.0)
    };
    (p1, dp)
}

impl Gll {
    pub fn new(p: usize) -> Self {
        assert!(p >= 1, "degree must be at least 1");
        let n = p + 1;
        let mut nodes = vec![0.0; n];
        nodes[0] = -1.0;
        nodes[p] = 1.0;
        let pf = (p * (p + 1)) as f64;
        for (j, node) in nodes.iter_mut().enumerate().take(p).skip(1) {
            // Chebyshev–Gauss–Lobatto start, Newton on P_p'
            let mut x = -(std::f64::consts::PI * j as f64 / p as f64).cos();
            for _ in 0..100 {
                let (v, dv) = legendre(p, x);
                let ddv = (2.0 * x * dv - pf * v) / (1.0 - x * x);
                let dx = dv / ddv;
                x -= dx;
                if dx.abs() < 1e-16 {
                    break;
                }
            }
            *node = x;
        }
        nodes.sort_by(f64::total_cmp);
        let pv: Vec<f64> = nodes.iter().map(|&x| legendre(p, x).0).collect();
        let weights: Vec<f64> = pv.iter().map(|v| 2.0 / (pf * v * v)).collect();
        let mut d = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                d[i * n + j] = if i != j {
                    pv[i] / (pv[j] * (nodes[i] - nodes[j]))
                } else if i == 0 {
                    -pf / 4.0
                } else if i == p {
                    pf / 4.0
                } else {
                    0.0
                };
            }
        }
        Self { p, nodes, weights, d }
    }

    pub fn len(&self) -> usize {
        self.p + 1
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// The 1-D stiffness matrix `DᵀWD` on the reference interval (row-major).
    pub fn stiffness_1d(&self) -> Vec<f64> {
        let n = self.len();
        let mut k = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                k[a * n + b] = (0..n).map(|q| self.d[q * n + a] * self.weights[q] * self.d[q * n + b]).sum();
            }
        }
        k
    }
}

/// Homogeneous Dirichlet conditions on the whole boundary, or natural (none).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    Dirichlet,
    Natural,
}

/// The operator `Σ_e (κ_e K_e + μ_e M_e)` on a structured box mesh (see the module documentation).
#[derive(Debug, Clone)]
pub struct SpectralBox<T: Scalar> {
    dim: usize,
    gll: Gll,
    /// Element widths per direction (`widths[d].len()` = number of elements in direction `d`).
    widths: [Vec<f64>; 3],
    nel: [usize; 3],
    /// Global node counts per direction.
    nn: [usize; 3],
    kappa: Vec<T>,
    mu: Vec<T>,
    boundary: Boundary,
    k1: Vec<f64>,
    /// Transpose of the derivative matrix.
    dt: Vec<f64>,
    /// Mapped geometry: per element and local node `d(d+1)/2` geometric factors `G_ab` and the weight `w·|det J|`
    /// (see [`SpectralBox::mapped`]); `None` for axis-aligned boxes.
    geo: Option<Vec<f64>>,
}

impl<T: Scalar> SpectralBox<T> {
    /// `widths[d]` for `d < dim` lists the element widths along axis `d`; `kappa` and `mu` have one entry per
    /// element in x-fastest order.
    pub fn new(p: usize, widths: &[Vec<f64>], kappa: Vec<T>, mu: Vec<T>, boundary: Boundary) -> Self {
        let dim = widths.len();
        assert!((1..=3).contains(&dim), "dimension must be 1, 2 or 3");
        let mut w: [Vec<f64>; 3] = [vec![1.0], vec![1.0], vec![1.0]];
        let mut nel = [1usize; 3];
        let mut nn = [1usize; 3];
        for d in 0..dim {
            assert!(!widths[d].is_empty() && widths[d].iter().all(|&h| h > 0.0));
            w[d] = widths[d].clone();
            nel[d] = widths[d].len();
            nn[d] = nel[d] * p + 1;
        }
        let total: usize = nel.iter().product();
        assert_eq!(kappa.len(), total, "one stiffness coefficient per element");
        assert_eq!(mu.len(), total, "one mass coefficient per element");
        let gll = Gll::new(p);
        let k1 = gll.stiffness_1d();
        let n1 = p + 1;
        let dt: Vec<f64> = (0..n1 * n1).map(|k| gll.d[(k % n1) * n1 + k / n1]).collect();
        Self { dim, gll, widths: w, nel, nn, kappa, mu, boundary, k1, dt, geo: None }
    }

    /// Non-affine (curved, skewed) elements: the same structured topology of `n_el[d]` elements per axis, but the
    /// physical position of every global node is given by `coords[c][g]` (component `c`, node `g`, x-fastest;
    /// `coords.len() == n_el.len()`). The element map is the isoparametric degree-`p` map through the nodes.
    ///
    /// The stiffness becomes `K_e = Σ_{a,b} D_aᵀ (G_ab ⊙ D_b ·)` with the symmetric geometric tensor
    /// `G_ab = w |det J| Σ_c (J⁻¹)_{ac}(J⁻¹)_{bc}` stored per quadrature point: *every* pair of directions is coupled,
    /// so the assembled element matrix is dense (`(p+1)^{2d}` entries) while the sum-factorised kernel still costs
    /// `O(d²(p+1)^{d+1})`. This is the setting in which matrix-free application pays.
    pub fn mapped(
        p: usize,
        n_el: &[usize],
        coords: &[Vec<f64>],
        kappa: Vec<T>,
        mu: Vec<T>,
        boundary: Boundary,
    ) -> Self {
        let dim = n_el.len();
        assert_eq!(coords.len(), dim, "one coordinate array per spatial dimension");
        let widths: Vec<Vec<f64>> = n_el.iter().map(|&n| vec![1.0; n]).collect();
        let mut op = Self::new(p, &widths, kappa, mu, boundary);
        let n1 = p + 1;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let tot = nl[0] * nl[1] * nl[2];
        assert!(coords.iter().all(|c| c.len() == op.len()), "coordinate arrays must have one entry per node");
        let ng = dim * (dim + 1) / 2;
        let mut geo = vec![0.0; op.num_elements() * tot * (ng + 1)];
        let w = &op.gll.weights;
        let mut xl = vec![vec![0.0f64; tot]; dim];
        let mut jac = vec![vec![vec![0.0f64; tot]; dim]; dim]; // jac[c][d][point]
        for ez in 0..op.nel[2] {
            for ey in 0..op.nel[1] {
                for ex in 0..op.nel[0] {
                    let ie = ex + op.nel[0] * (ey + op.nel[1] * ez);
                    for c in 0..nl[2] {
                        for b in 0..nl[1] {
                            for a in 0..nl[0] {
                                let gx = ex * p + a;
                                let gy = if dim > 1 { ey * p + b } else { 0 };
                                let gz = if dim > 2 { ez * p + c } else { 0 };
                                let g = gx + op.nn[0] * (gy + op.nn[1] * gz);
                                let l = a + nl[0] * (b + nl[1] * c);
                                for (k, xk) in xl.iter_mut().enumerate() {
                                    xk[l] = coords[k][g];
                                }
                            }
                        }
                    }
                    for c in 0..dim {
                        for d in 0..dim {
                            SpectralBox::<f64>::axis_apply(&op.gll.d, d, nl, &xl[c], &mut jac[c][d]);
                        }
                    }
                    for l in 0..tot {
                        let (a, b, c) = (l % nl[0], (l / nl[0]) % nl[1], l / (nl[0] * nl[1]));
                        let wp = w[a] * if dim > 1 { w[b] } else { 1.0 } * if dim > 2 { w[c] } else { 1.0 };
                        let mut j = [[0.0f64; 3]; 3];
                        for r in 0..dim {
                            for q in 0..dim {
                                j[r][q] = jac[r][q][l];
                            }
                        }
                        let (det, inv) = invert_small(&j, dim);
                        assert!(det != 0.0 && det.is_finite(), "degenerate element {ie}");
                        let base = (ie * tot + l) * (ng + 1);
                        let wdet = wp * det.abs();
                        let mut idx = 0;
                        for a2 in 0..dim {
                            for b2 in a2..dim {
                                geo[base + idx] = wdet * (0..dim).map(|c2| inv[a2][c2] * inv[b2][c2]).sum::<f64>();
                                idx += 1;
                            }
                        }
                        geo[base + ng] = wdet;
                    }
                }
            }
        }
        op.geo = Some(geo);
        op
    }

    /// Uniform mesh of `n_el[d]` elements of width `h[d]` with constant coefficients.
    pub fn uniform(p: usize, n_el: &[usize], h: &[f64], kappa: T, mu: T, boundary: Boundary) -> Self {
        let widths: Vec<Vec<f64>> = n_el.iter().zip(h).map(|(&n, &w)| vec![w; n]).collect();
        let total: usize = n_el.iter().product();
        Self::new(p, &widths, vec![kappa; total], vec![mu; total], boundary)
    }

    pub fn degree(&self) -> usize {
        self.gll.p
    }
    pub fn spatial_dim(&self) -> usize {
        self.dim
    }
    pub fn num_elements(&self) -> usize {
        self.nel.iter().product()
    }
    /// Global nodes per direction (the dofs are laid out x-fastest).
    pub fn nodes_per_axis(&self) -> &[usize] {
        &self.nn[..self.dim]
    }
    pub fn len(&self) -> usize {
        self.nn.iter().product()
    }
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Physical coordinates of node `(i, j, k)` along each axis: `coord(d, i)`.
    pub fn coord(&self, d: usize, i: usize) -> f64 {
        let p = self.gll.p;
        let (e, a) = if i == self.nn[d] - 1 { (self.nel[d] - 1, p) } else { (i / p, i % p) };
        let x0: f64 = self.widths[d][..e].iter().sum();
        x0 + 0.5 * self.widths[d][e] * (self.gll.nodes[a] + 1.0)
    }

    /// Whether global node `g` lies on the boundary (and so is constrained under [`Boundary::Dirichlet`]).
    pub fn on_boundary(&self, g: usize) -> bool {
        let mut r = g;
        for d in 0..self.dim {
            let i = r % self.nn[d];
            r /= self.nn[d];
            if i == 0 || i == self.nn[d] - 1 {
                return true;
            }
        }
        false
    }

    /// Multiply–add operations of one operator application.
    pub fn flops_per_apply(&self) -> f64 {
        let n1 = (self.gll.p + 1) as f64;
        let per = 2.0 * self.dim as f64 * n1.powi(self.dim as i32 + 1) * 2.0 + 6.0 * n1.powi(self.dim as i32);
        per * self.num_elements() as f64
    }

    /// Applies `κ_e K_e + μ_e M_e` for element `(ex, ey, ez)` to the local vector `u` (x-fastest), *adding* into
    /// `out`. `scratch` must hold `(dim+1) · (p+1)^dim` values.
    fn apply_element(&self, e: [usize; 3], u: &[T], out: &mut [T], scratch: &mut [T]) {
        if self.geo.is_some() {
            return self.apply_element_mapped(e, u, out, scratch);
        }
        let n1 = self.gll.len();
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let tot = nl[0] * nl[1] * nl[2];
        let ie = e[0] + self.nel[0] * (e[1] + self.nel[1] * e[2]);
        let h = [self.widths[0][e[0]], self.widths[1][e[1]], self.widths[2][e[2]]];
        let det: f64 = (0..dim).map(|d| 0.5 * h[d]).product();
        let w = &self.gll.weights;
        let wl = |a: usize, b: usize, c: usize| -> f64 {
            w[a] * if dim > 1 { w[b] } else { 1.0 } * if dim > 2 { w[c] } else { 1.0 }
        };
        let (g, rest) = scratch.split_at_mut(tot);
        let (t, _) = rest.split_at_mut(tot);
        // mass:  out += μ detJ W ⊙ u
        let mu = self.mu[ie];
        for c in 0..nl[2] {
            for b in 0..nl[1] {
                for a in 0..nl[0] {
                    let i = a + nl[0] * (b + nl[1] * c);
                    out[i] += (mu * u[i]).scale(det * wl(a, b, c));
                }
            }
        }
        // stiffness:  out += κ Σ_d (2/h_d)² detJ · D_dᵀ (W ⊙ (D_d u))
        let kappa = self.kappa[ie];
        let dm = &self.gll.d;
        for d in 0..dim {
            let s = (2.0 / h[d]).powi(2) * det;
            Self::axis_apply(dm, d, nl, u, g);
            for c in 0..nl[2] {
                for b in 0..nl[1] {
                    for a in 0..nl[0] {
                        let i = a + nl[0] * (b + nl[1] * c);
                        g[i] = (kappa * g[i]).scale(s * wl(a, b, c));
                    }
                }
            }
            Self::axis_apply(&self.dt, d, nl, g, t);
            for (o, tv) in out.iter_mut().zip(t.iter()) {
                *o += *tv;
            }
        }
    }

    /// Kernel for mapped geometry: `out += μ w|J| u + κ Σ_a D_aᵀ (Σ_b G_ab D_b u)`. `scratch`: `(dim+1)·tot` values.
    fn apply_element_mapped(&self, e: [usize; 3], u: &[T], out: &mut [T], scratch: &mut [T]) {
        let geo = self.geo.as_ref().expect("mapped geometry");
        let n1 = self.gll.len();
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let tot = nl[0] * nl[1] * nl[2];
        let ng = dim * (dim + 1) / 2;
        let ie = e[0] + self.nel[0] * (e[1] + self.nel[1] * e[2]);
        let g = &geo[ie * tot * (ng + 1)..(ie + 1) * tot * (ng + 1)];
        let (mu, kappa) = (self.mu[ie], self.kappa[ie]);
        for i in 0..tot {
            out[i] += (mu * u[i]).scale(g[i * (ng + 1) + ng]);
        }
        let (grads, t) = scratch.split_at_mut(dim * tot);
        for a in 0..dim {
            Self::axis_apply(&self.gll.d, a, nl, u, &mut grads[a * tot..(a + 1) * tot]);
        }
        let mut q = [T::ZERO; 3];
        for i in 0..tot {
            let gi = &g[i * (ng + 1)..i * (ng + 1) + ng];
            let gr = [
                grads[i],
                if dim > 1 { grads[tot + i] } else { T::ZERO },
                if dim > 2 { grads[2 * tot + i] } else { T::ZERO },
            ];
            let mut idx = 0;
            q[..dim].fill(T::ZERO);
            for a in 0..dim {
                for b in a..dim {
                    let gab = gi[idx];
                    idx += 1;
                    q[a] += gr[b].scale(gab);
                    if b != a {
                        q[b] += gr[a].scale(gab);
                    }
                }
            }
            for a in 0..dim {
                grads[a * tot + i] = kappa * q[a];
            }
        }
        for a in 0..dim {
            Self::axis_apply(&self.dt, a, nl, &grads[a * tot..(a + 1) * tot], &mut t[..tot]);
            for (o, tv) in out.iter_mut().zip(t[..tot].iter()) {
                *o += *tv;
            }
        }
    }

    /// `y = (I ⊗ M ⊗ I) x`: the `n×n` row-major matrix `m` (`y_i = Σ_j m[i·n+j] x_j`) acts along axis `d` of the
    /// tensor `x` with extents `nl`. Callers pass `D` or its transpose.
    #[inline]
    fn axis_apply(m: &[f64], d: usize, nl: [usize; 3], x: &[T], y: &mut [T]) {
        let n = nl[d];
        let stride: usize = nl[..d].iter().product();
        let outer: usize = nl[d + 1..].iter().product();
        if stride == 1 {
            // contiguous lines: dot products
            for (xs, ys) in x.chunks_exact(n).zip(y.chunks_exact_mut(n)).take(outer) {
                for (yi, mrow) in ys.iter_mut().zip(m.chunks_exact(n)) {
                    let mut acc = T::ZERO;
                    for (xv, &c) in xs.iter().zip(mrow) {
                        acc += xv.scale(c);
                    }
                    *yi = acc;
                }
            }
            return;
        }
        let block = stride * n;
        for (xb, yb) in x.chunks_exact(block).zip(y.chunks_exact_mut(block)).take(outer) {
            for (yrow, mrow) in yb.chunks_exact_mut(stride).zip(m.chunks_exact(n)) {
                yrow.fill(T::ZERO);
                for (xrow, &c) in xb.chunks_exact(stride).zip(mrow) {
                    for (yv, xv) in yrow.iter_mut().zip(xrow) {
                        *yv += xv.scale(c);
                    }
                }
            }
        }
    }

    /// Adds all elements of `z`-layer `ez` into `ylayer`, the slice of the output covering node planes
    /// `ez·p ..= (ez+1)·p` (the single plane 0 in 1-D/2-D).
    fn apply_layer(&self, ez: usize, x: &[T], ylayer: &mut [T]) {
        let p = self.gll.p;
        let n1 = p + 1;
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let tot = nl[0] * nl[1] * nl[2];
        let plane = self.nn[0] * self.nn[1];
        let z0 = if dim > 2 { ez * p } else { 0 };
        let mut u = vec![T::ZERO; tot];
        let mut o = vec![T::ZERO; tot];
        let mut scratch = vec![T::ZERO; (dim + 1).max(2) * tot];
        for ey in 0..self.nel[1] {
            for ex in 0..self.nel[0] {
                // gather
                for c in 0..nl[2] {
                    for b in 0..nl[1] {
                        let gy = if dim > 1 { ey * p + b } else { 0 };
                        let row = (z0 + c) * plane + gy * self.nn[0] + ex * p;
                        let l = nl[0] * (b + nl[1] * c);
                        u[l..l + nl[0]].copy_from_slice(&x[row..row + nl[0]]);
                    }
                }
                o.fill(T::ZERO);
                self.apply_element([ex, ey, ez], &u, &mut o, &mut scratch);
                // scatter-add into the layer slice (plane index relative to z0)
                for c in 0..nl[2] {
                    for b in 0..nl[1] {
                        let gy = if dim > 1 { ey * p + b } else { 0 };
                        let row = c * plane + gy * self.nn[0] + ex * p;
                        let l = nl[0] * (b + nl[1] * c);
                        for a in 0..nl[0] {
                            ylayer[row + a] += o[l + a];
                        }
                    }
                }
            }
        }
    }

    /// `y = A x` (constrained nodes: identity rows, so the operator stays symmetric).
    pub fn apply_into(&self, x: &[T], y: &mut [T]) {
        let n = self.len();
        assert!(x.len() == n && y.len() == n);
        let xm: Vec<T>;
        let xin: &[T] = if self.boundary == Boundary::Dirichlet {
            xm = x.iter().enumerate().map(|(g, &v)| if self.on_boundary(g) { T::ZERO } else { v }).collect();
            &xm
        } else {
            x
        };
        y.fill(T::ZERO);
        let p = self.gll.p;
        let plane = self.nn[0] * self.nn[1];
        let nlayers = self.nel[2];
        // two passes: even layers, then odd layers; each layer owns a disjoint slice of y within a pass
        for parity in 0..2usize {
            let mut slices: Vec<(usize, &mut [T])> = Vec::new();
            let mut rest: &mut [T] = &mut y[..];
            let mut offset = 0usize; // first plane index of `rest`
            for ez in (parity..nlayers).step_by(2) {
                let (z_lo, z_hi) = if self.dim > 2 { (ez * p, (ez + 1) * p + 1) } else { (0, 1) };
                let skip = (z_lo - offset) * plane;
                let take = (z_hi - z_lo) * plane;
                let (_, tail) = std::mem::take(&mut rest).split_at_mut(skip);
                let (head, tail) = tail.split_at_mut(take);
                slices.push((ez, head));
                rest = tail;
                offset = z_hi;
            }
            crate::solvers::par::for_each_mut(&mut slices, crate::solvers::par::worth_parallel(n), |(ez, ys)| {
                self.apply_layer(*ez, xin, ys)
            });
        }
        if self.boundary == Boundary::Dirichlet {
            for (g, (yv, &xv)) in y.iter_mut().zip(x).enumerate() {
                if self.on_boundary(g) {
                    *yv = xv;
                }
            }
        }
    }

    /// The exact diagonal of the operator (constrained nodes: 1), `O(#elements · (p+1)^dim)`.
    pub fn diagonal(&self) -> Vec<T> {
        let n1 = self.gll.len();
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let p = self.gll.p;
        let w = &self.gll.weights;
        let mut diag = vec![T::ZERO; self.len()];
        for ez in 0..self.nel[2] {
            for ey in 0..self.nel[1] {
                for ex in 0..self.nel[0] {
                    let ie = ex + self.nel[0] * (ey + self.nel[1] * ez);
                    let e = [ex, ey, ez];
                    let h = [self.widths[0][ex], self.widths[1][ey], self.widths[2][ez]];
                    let det: f64 = (0..dim).map(|d| 0.5 * h[d]).product();
                    for c in 0..nl[2] {
                        for b in 0..nl[1] {
                            for a in 0..nl[0] {
                                let loc = [a, b, c];
                                let wprod: f64 = (0..dim).map(|d| w[loc[d]]).product();
                                let mut v = self.mu[ie].scale(det * wprod);
                                if let Some(geo) = &self.geo {
                                    let ng = dim * (dim + 1) / 2;
                                    let tot = nl[0] * nl[1] * nl[2];
                                    let at =
                                        |l: [usize; 3]| (ie * tot + l[0] + nl[0] * (l[1] + nl[1] * l[2])) * (ng + 1);
                                    let dd = &self.gll.d;
                                    let mut k = 0.0;
                                    let mut idx = 0;
                                    for a2 in 0..dim {
                                        for b2 in a2..dim {
                                            if a2 == b2 {
                                                for q in 0..n1 {
                                                    let mut l = loc;
                                                    l[a2] = q;
                                                    let dq = dd[q * n1 + loc[a2]];
                                                    k += geo[at(l) + idx] * dq * dq;
                                                }
                                            } else {
                                                k += 2.0
                                                    * geo[at(loc) + idx]
                                                    * dd[loc[a2] * n1 + loc[a2]]
                                                    * dd[loc[b2] * n1 + loc[b2]];
                                            }
                                            idx += 1;
                                        }
                                    }
                                    v = self.mu[ie].scale(geo[at(loc) + ng]) + self.kappa[ie].scale(k);
                                }
                                for d in 0..dim {
                                    if self.geo.is_some() {
                                        break;
                                    }
                                    let s = (2.0 / h[d]).powi(2) * det;
                                    let wo: f64 = (0..dim).filter(|&q| q != d).map(|q| w[loc[q]]).product();
                                    v += self.kappa[ie].scale(s * wo * self.k1[loc[d] * n1 + loc[d]]);
                                }
                                let gx = e[0] * p + a;
                                let gy = if dim > 1 { e[1] * p + b } else { 0 };
                                let gz = if dim > 2 { e[2] * p + c } else { 0 };
                                diag[gx + self.nn[0] * (gy + self.nn[1] * gz)] += v;
                            }
                        }
                    }
                }
            }
        }
        if self.boundary == Boundary::Dirichlet {
            for (g, d) in diag.iter_mut().enumerate() {
                if self.on_boundary(g) {
                    *d = T::ONE;
                }
            }
        }
        diag
    }

    /// Mapped geometry: the dense element matrices are obtained by applying the element kernel to the unit vectors
    /// (so this path validates the *assembly* and gives the memory/time baseline, not the kernel itself).
    fn assemble_mapped(&self) -> Csr<T> {
        let n1 = self.gll.len();
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let tot = nl[0] * nl[1] * nl[2];
        let p = self.gll.p;
        let mut trip: Vec<(usize, usize, T)> = Vec::new();
        let mut scratch = vec![T::ZERO; (dim + 1).max(2) * tot];
        let (mut u, mut o) = (vec![T::ZERO; tot], vec![T::ZERO; tot]);
        for ez in 0..self.nel[2] {
            for ey in 0..self.nel[1] {
                for ex in 0..self.nel[0] {
                    let e = [ex, ey, ez];
                    let global = |l: usize| -> usize {
                        let (a, b, c) = (l % nl[0], (l / nl[0]) % nl[1], l / (nl[0] * nl[1]));
                        let gx = ex * p + a;
                        let gy = if dim > 1 { ey * p + b } else { 0 };
                        let gz = if dim > 2 { ez * p + c } else { 0 };
                        gx + self.nn[0] * (gy + self.nn[1] * gz)
                    };
                    for lj in 0..tot {
                        u.fill(T::ZERO);
                        u[lj] = T::ONE;
                        o.fill(T::ZERO);
                        self.apply_element_mapped(e, &u, &mut o, &mut scratch);
                        for (li, &v) in o.iter().enumerate() {
                            if v != T::ZERO {
                                trip.push((global(li), global(lj), v));
                            }
                        }
                    }
                }
            }
        }
        let n = self.len();
        if self.boundary == Boundary::Dirichlet {
            trip.retain(|&(i, j, _)| !self.on_boundary(i) && !self.on_boundary(j));
            for g in 0..n {
                if self.on_boundary(g) {
                    trip.push((g, g, T::ONE));
                }
            }
        }
        Csr::from_triplets(n, n, &trip).expect("valid triplets")
    }

    /// Assembles the operator as a CSR matrix from the Kronecker structure of the element matrices
    /// (`κ Σ_d s_d ⊗_e A_e + μ detJ W`), independently of the element kernel.
    pub fn assemble(&self) -> Csr<T> {
        if self.geo.is_some() {
            return self.assemble_mapped();
        }
        let n1 = self.gll.len();
        let dim = self.dim;
        let nl: [usize; 3] = [n1, if dim > 1 { n1 } else { 1 }, if dim > 2 { n1 } else { 1 }];
        let p = self.gll.p;
        let w = &self.gll.weights;
        let loc_count = nl[0] * nl[1] * nl[2];
        let mut trip: Vec<(usize, usize, T)> = Vec::new();
        for ez in 0..self.nel[2] {
            for ey in 0..self.nel[1] {
                for ex in 0..self.nel[0] {
                    let ie = ex + self.nel[0] * (ey + self.nel[1] * ez);
                    let e = [ex, ey, ez];
                    let h = [self.widths[0][ex], self.widths[1][ey], self.widths[2][ez]];
                    let det: f64 = (0..dim).map(|d| 0.5 * h[d]).product();
                    let global = |l: usize| -> usize {
                        let (a, b, c) = (l % nl[0], (l / nl[0]) % nl[1], l / (nl[0] * nl[1]));
                        let gx = e[0] * p + a;
                        let gy = if dim > 1 { e[1] * p + b } else { 0 };
                        let gz = if dim > 2 { e[2] * p + c } else { 0 };
                        gx + self.nn[0] * (gy + self.nn[1] * gz)
                    };
                    for li in 0..loc_count {
                        let ci = [li % nl[0], (li / nl[0]) % nl[1], li / (nl[0] * nl[1])];
                        for lj in 0..loc_count {
                            let cj = [lj % nl[0], (lj / nl[0]) % nl[1], lj / (nl[0] * nl[1])];
                            let mut kval = 0.0;
                            for d in 0..dim {
                                let s = (2.0 / h[d]).powi(2) * det;
                                let mut term = s * self.k1[ci[d] * n1 + cj[d]];
                                for q in (0..dim).filter(|&q| q != d) {
                                    term *= if ci[q] == cj[q] { w[ci[q]] } else { 0.0 };
                                }
                                kval += term;
                            }
                            let mut mval = 0.0;
                            if li == lj {
                                mval = det * (0..dim).map(|d| w[ci[d]]).product::<f64>();
                            }
                            let v = self.kappa[ie].scale(kval) + self.mu[ie].scale(mval);
                            if v != T::ZERO {
                                trip.push((global(li), global(lj), v));
                            }
                        }
                    }
                }
            }
        }
        let n = self.len();
        if self.boundary == Boundary::Dirichlet {
            let keep = |g: usize| !self.on_boundary(g);
            trip.retain(|&(i, j, _)| keep(i) && keep(j));
            for g in 0..n {
                if self.on_boundary(g) {
                    trip.push((g, g, T::ONE));
                }
            }
        }
        Csr::from_triplets(n, n, &trip).expect("valid triplets")
    }
}

/// Determinant and inverse of the leading `dim × dim` block.
fn invert_small(j: &[[f64; 3]; 3], dim: usize) -> (f64, [[f64; 3]; 3]) {
    let mut inv = [[0.0f64; 3]; 3];
    match dim {
        1 => {
            inv[0][0] = 1.0 / j[0][0];
            (j[0][0], inv)
        }
        2 => {
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            inv[0][0] = j[1][1] / det;
            inv[0][1] = -j[0][1] / det;
            inv[1][0] = -j[1][0] / det;
            inv[1][1] = j[0][0] / det;
            (det, inv)
        }
        _ => {
            let c = |a: usize, b: usize, c: usize, d: usize| j[a][c] * j[b][d] - j[a][d] * j[b][c];
            let det = j[0][0] * c(1, 2, 1, 2) - j[0][1] * c(1, 2, 0, 2) + j[0][2] * c(1, 2, 0, 1);
            inv[0][0] = c(1, 2, 1, 2) / det;
            inv[0][1] = -c(0, 2, 1, 2) / det;
            inv[0][2] = c(0, 1, 1, 2) / det;
            inv[1][0] = -c(1, 2, 0, 2) / det;
            inv[1][1] = c(0, 2, 0, 2) / det;
            inv[1][2] = -c(0, 1, 0, 2) / det;
            inv[2][0] = c(1, 2, 0, 1) / det;
            inv[2][1] = -c(0, 2, 0, 1) / det;
            inv[2][2] = c(0, 1, 0, 1) / det;
            (det, inv)
        }
    }
}

impl<T: Scalar> LinearOperator<T> for SpectralBox<T> {
    fn dim(&self) -> usize {
        self.len()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.apply_into(x, y)
    }
}
