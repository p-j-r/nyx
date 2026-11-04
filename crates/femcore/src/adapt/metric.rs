//! Riemannian metrics for mesh adaptation, generic over dimension `D` (2 or 3).
//!
//! A metric is a symmetric positive-definite matrix `M`; the length of a vector `e` is
//! `sqrt(e^T M e)`. An ideal ("unit") mesh has all edges of length 1. Along an eigenvector with
//! eigenvalue `lambda` the requested element size is `h = 1 / sqrt(lambda)`.
//!
//! Contents:
//! * dense symmetric algebra for 2x2 and 3x3 (`sym_eig`, `map_eig`, `det`, `inv`);
//! * lengths: `length`, `edge_length` (geometric size variation along an edge);
//! * `intersect`: simultaneous reduction, the largest ellipsoid inside both unit balls;
//! * `log_euclid_mean`: interpolation;
//! * `grow` / `grade`: anisotropic gradation with a linear size-growth law;
//! * `lp_metrics`: derivative (Hessian) to metric with L^p optimal scaling and complexity
//!   targeting, as a fused kernel, plus `lp_metrics_unfused` as the straightforward multi-pass
//!   reference used for benchmarking and tests.
//!
//! Everything here is written from the published formulae (see `docs/REFERENCES.md`); no third-party
//! code was read or copied.

use crate::par;

/// Symmetric `D x D` matrix (full storage, row-major).
pub type Sym<const D: usize> = [[f64; D]; D];

/// Identity matrix.
pub fn identity<const D: usize>() -> Sym<D> {
    let mut a = [[0.0; D]; D];
    for i in 0..D {
        a[i][i] = 1.0;
    }
    a
}

/// `e^T a e`.
pub fn quad_form<const D: usize>(a: &Sym<D>, e: &[f64; D]) -> f64 {
    let mut s = 0.0;
    for i in 0..D {
        let mut r = 0.0;
        for j in 0..D {
            r += a[i][j] * e[j];
        }
        s += e[i] * r;
    }
    s
}

/// Matrix product `a * b`.
pub fn mul<const D: usize>(a: &Sym<D>, b: &Sym<D>) -> [[f64; D]; D] {
    let mut c = [[0.0; D]; D];
    for i in 0..D {
        for k in 0..D {
            let aik = a[i][k];
            for j in 0..D {
                c[i][j] += aik * b[k][j];
            }
        }
    }
    c
}

/// Transpose.
pub fn transpose<const D: usize>(a: &[[f64; D]; D]) -> [[f64; D]; D] {
    let mut t = [[0.0; D]; D];
    for i in 0..D {
        for j in 0..D {
            t[j][i] = a[i][j];
        }
    }
    t
}

/// Force exact symmetry: `(a + a^T) / 2`.
pub fn symmetrize<const D: usize>(a: &[[f64; D]; D]) -> Sym<D> {
    let mut s = *a;
    for i in 0..D {
        for j in i + 1..D {
            let v = 0.5 * (a[i][j] + a[j][i]);
            s[i][j] = v;
            s[j][i] = v;
        }
    }
    s
}

/// Determinant (closed form for `D = 2, 3`).
pub fn det<const D: usize>(a: &Sym<D>) -> f64 {
    match D {
        2 => a[0][0] * a[1][1] - a[0][1] * a[1][0],
        3 => {
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        }
        _ => panic!("only D = 2 and D = 3 are supported"),
    }
}

/// Inverse by cofactors (`D = 2, 3`).
pub fn inv<const D: usize>(a: &Sym<D>) -> Sym<D> {
    let d = det(a);
    let mut r = [[0.0; D]; D];
    match D {
        2 => {
            r[0][0] = a[1][1] / d;
            r[0][1] = -a[0][1] / d;
            r[1][0] = -a[1][0] / d;
            r[1][1] = a[0][0] / d;
        }
        3 => {
            r[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / d;
            r[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / d;
            r[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / d;
            r[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / d;
            r[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / d;
            r[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / d;
            r[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / d;
            r[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / d;
            r[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / d;
        }
        _ => panic!("only D = 2 and D = 3 are supported"),
    }
    r
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
///
/// Returns `(values, vectors)` with values in descending order and `vectors[i][k]` the `i`-th
/// component of eigenvector `k`. Deterministic: a fixed sweep order and no data-dependent
/// branching beyond the convergence test.
pub fn sym_eig<const D: usize>(a: &Sym<D>) -> ([f64; D], Sym<D>) {
    let mut a = *a;
    let mut v = identity::<D>();
    let scale: f64 = a.iter().flatten().map(|x| x * x).sum::<f64>();
    let tol = 1e-30 * scale;
    for _ in 0..32 {
        let mut off = 0.0;
        for p in 0..D {
            for q in p + 1..D {
                off += a[p][q] * a[p][q];
            }
        }
        if off <= tol {
            break;
        }
        for p in 0..D {
            for q in p + 1..D {
                let apq = a[p][q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..D {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..D {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..D {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut idx: [usize; D] = std::array::from_fn(|i| i);
    idx.sort_by(|&i, &j| a[j][j].partial_cmp(&a[i][i]).unwrap_or(std::cmp::Ordering::Equal).then(i.cmp(&j)));
    let vals = std::array::from_fn(|k| a[idx[k]][idx[k]]);
    let mut vecs = [[0.0; D]; D];
    for k in 0..D {
        for i in 0..D {
            vecs[i][k] = v[i][idx[k]];
        }
    }
    (vals, vecs)
}

/// `V diag(vals) V^T`.
pub fn from_eig<const D: usize>(vals: &[f64; D], vecs: &Sym<D>) -> Sym<D> {
    let mut a = [[0.0; D]; D];
    for k in 0..D {
        for i in 0..D {
            let w = vals[k] * vecs[i][k];
            for j in 0..D {
                a[i][j] += w * vecs[j][k];
            }
        }
    }
    symmetrize(&a)
}

/// Apply a scalar function to the eigenvalues of `a`.
pub fn map_eig<const D: usize>(a: &Sym<D>, f: impl Fn(f64) -> f64) -> Sym<D> {
    let (vals, vecs) = sym_eig(a);
    from_eig(&vals.map(f), &vecs)
}

/// Isotropic metric of size `h` (`M = h^-2 I`).
pub fn isotropic<const D: usize>(h: f64) -> Sym<D> {
    let mut a = [[0.0; D]; D];
    let l = 1.0 / (h * h);
    for i in 0..D {
        a[i][i] = l;
    }
    a
}

/// Metric with the requested sizes `h[k]` along the orthonormal directions `dirs[.][k]`.
pub fn from_sizes<const D: usize>(h: &[f64; D], dirs: &Sym<D>) -> Sym<D> {
    from_eig(&h.map(|x| 1.0 / (x * x)), dirs)
}

/// Requested sizes and directions: `(h, dirs)` with `h[k] = lambda_k^{-1/2}` and sizes ascending.
pub fn sizes<const D: usize>(m: &Sym<D>) -> ([f64; D], Sym<D>) {
    let (vals, vecs) = sym_eig(m);
    (vals.map(|l| 1.0 / l.sqrt()), vecs)
}

/// Length of vector `e` in metric `m`.
pub fn length<const D: usize>(m: &Sym<D>, e: &[f64; D]) -> f64 {
    quad_form(m, e).sqrt()
}

/// Length of the edge `e` between two vertices carrying metrics `ma` and `mb`, assuming the
/// requested size varies geometrically along the edge: with `la`, `lb` the lengths in the two
/// metrics, `l = (la - lb) / ln(la / lb)` (and `la` when they coincide).
pub fn edge_length<const D: usize>(ma: &Sym<D>, mb: &Sym<D>, e: &[f64; D]) -> f64 {
    let la = length(ma, e);
    let lb = length(mb, e);
    let r = la / lb;
    if (r - 1.0).abs() < 1e-6 { 0.5 * (la + lb) } else { (la - lb) / r.ln() }
}

/// Intersection of two metrics by simultaneous reduction: the metric whose unit ball is the
/// largest ellipsoid contained in both unit balls. Also returns whether `m2` was more restrictive
/// than `m1` in some direction by more than `1e-12` relatively (i.e. whether the result differs
/// from `m1`).
pub fn intersect_flag<const D: usize>(m1: &Sym<D>, m2: &Sym<D>) -> (Sym<D>, bool) {
    // R = M1^{1/2}, S = M1^{-1/2}; N = S M2 S = Q L Q^T; result = R Q max(1, L) Q^T R
    let (l1, q1) = sym_eig(m1);
    let r = from_eig(&l1.map(f64::sqrt), &q1);
    let s = from_eig(&l1.map(|x| 1.0 / x.sqrt()), &q1);
    let n = symmetrize(&mul(&mul(&s, m2), &s));
    let (ln, qn) = sym_eig(&n);
    let changed = ln.iter().any(|&x| x > 1.0 + 1e-12);
    if !changed {
        return (*m1, false);
    }
    let mid = from_eig(&ln.map(|x| x.max(1.0)), &qn);
    (symmetrize(&mul(&mul(&r, &mid), &r)), true)
}

/// Intersection of two metrics (see [`intersect_flag`]).
pub fn intersect<const D: usize>(m1: &Sym<D>, m2: &Sym<D>) -> Sym<D> {
    intersect_flag(m1, m2).0
}

/// Matrix logarithm of an SPD matrix.
pub fn log_spd<const D: usize>(m: &Sym<D>) -> Sym<D> {
    map_eig(m, f64::ln)
}

/// Matrix exponential of a symmetric matrix.
pub fn exp_sym<const D: usize>(m: &Sym<D>) -> Sym<D> {
    map_eig(m, f64::exp)
}

/// Log-Euclidean weighted mean `exp(sum w_i log M_i)` (weights should sum to 1).
pub fn log_euclid_mean<const D: usize>(items: &[(f64, Sym<D>)]) -> Sym<D> {
    let mut acc = [[0.0; D]; D];
    for (w, m) in items {
        let l = log_spd(m);
        for i in 0..D {
            for j in 0..D {
                acc[i][j] += w * l[i][j];
            }
        }
    }
    exp_sym(&symmetrize(&acc))
}

/// Metric obtained by letting the requested sizes of `m` grow linearly over a distance `dist`:
/// `h' = h + g * dist` along every eigendirection (eigenvalue `1 / (h + g dist)^2`).
pub fn grow<const D: usize>(m: &Sym<D>, dist: f64, g: f64) -> Sym<D> {
    map_eig(m, |l| {
        let h = 1.0 / l.sqrt() + g * dist;
        1.0 / (h * h)
    })
}

/// Gradation of a vertex-based metric field over the given edges.
///
/// For every edge `(i, j)` the metric at `j` is intersected with the grown metric of `i` and vice
/// versa, so that afterwards the requested size along any direction grows at most by `g` per unit
/// of Euclidean distance (`g` is the size-gradient bound, e.g. `0.3`). Sweeps over the edges in
/// order until no metric changes (at most `max_passes`). The intersection only ever shrinks sizes,
/// so refinement requests are never lost. Returns the number of passes used.
pub fn grade<const D: usize>(
    metrics: &mut [Sym<D>],
    points: &[[f64; D]],
    edges: &[[u32; 2]],
    g: f64,
    max_passes: usize,
) -> usize {
    for pass in 0..max_passes {
        let mut changed = false;
        for e in edges {
            let (i, j) = (e[0] as usize, e[1] as usize);
            let mut d2 = 0.0;
            for a in 0..D {
                let t = points[j][a] - points[i][a];
                d2 += t * t;
            }
            let dist = d2.sqrt();
            let (nj, cj) = intersect_flag(&metrics[j], &grow(&metrics[i], dist, g));
            if cj {
                metrics[j] = nj;
                changed = true;
            }
            let (ni, ci) = intersect_flag(&metrics[i], &grow(&metrics[j], dist, g));
            if ci {
                metrics[i] = ni;
                changed = true;
            }
        }
        if !changed {
            return pass + 1;
        }
    }
    max_passes
}

/// Metric volume of an ideal (unit-edge) simplex: `sqrt(3)/4` in 2D, `sqrt(2)/12` in 3D.
pub fn unit_simplex_measure(d: usize) -> f64 {
    match d {
        2 => 3f64.sqrt() / 4.0,
        3 => 2f64.sqrt() / 12.0,
        _ => panic!("only D = 2 and D = 3 are supported"),
    }
}

/// Complexity (integral of `sqrt(det M)`) that makes a unit mesh of about `n_elements` cells.
pub fn complexity_for_elements(n_elements: f64, d: usize) -> f64 {
    n_elements * unit_simplex_measure(d)
}

/// Approximate number of elements of a unit mesh for a metric field of the given complexity.
pub fn elements_for_complexity(complexity: f64, d: usize) -> f64 {
    complexity / unit_simplex_measure(d)
}

/// Complexity of a vertex-based metric field using lumped vertex volumes.
pub fn complexity<const D: usize>(metrics: &[Sym<D>], vol: &[f64]) -> f64 {
    metrics.iter().zip(vol).map(|(m, v)| v * det(m).sqrt()).sum()
}

/// Options for [`lp_metrics`].
#[derive(Clone, Debug)]
pub struct LpOptions {
    /// Norm exponent `p` of the interpolation-error functional (`f64::INFINITY` for the max norm).
    pub p: f64,
    /// Target complexity, see [`complexity_for_elements`].
    pub target_complexity: f64,
    /// Smallest allowed size.
    pub hmin: f64,
    /// Largest allowed size.
    pub hmax: f64,
    /// Largest allowed ratio between the largest and the smallest size at a vertex
    /// (`f64::INFINITY` for no limit).
    pub max_aspect: f64,
    /// Eigenvalues of `|H|` below `eps_rel` times the largest Hessian norm of the whole field are
    /// raised to that value.
    pub eps_rel: f64,
}

impl Default for LpOptions {
    fn default() -> Self {
        LpOptions {
            p: 2.0,
            target_complexity: 1000.0,
            hmin: 1e-12,
            hmax: 1e12,
            max_aspect: f64::INFINITY,
            eps_rel: 1e-10,
        }
    }
}

/// Result of [`lp_metrics`].
#[derive(Clone, Debug)]
pub struct LpResult<const D: usize> {
    pub metrics: Vec<Sym<D>>,
    /// Global scale applied to the unnormalised field.
    pub scale: f64,
    /// Complexity of the returned field.
    pub complexity: f64,
    /// True when the clamps made the target unreachable and the result sits at a bound.
    pub saturated: bool,
}

/// Eigen-data stored per vertex: unscaled eigenvalues and eigenvectors.
#[derive(Clone, Copy)]
struct Eig<const D: usize> {
    u: [f64; D],
    q: Sym<D>,
}

#[inline]
fn p_exponent(p: f64, d: usize) -> f64 {
    if p.is_finite() { -1.0 / (2.0 * p + d as f64) } else { 0.0 }
}

/// Eigenvalues of the scaled, clamped metric at one vertex.
#[inline]
fn scaled_eigs<const D: usize>(c: f64, u: &[f64; D], lo: f64, hi: f64, max_aspect: f64) -> [f64; D] {
    let mut l: [f64; D] = std::array::from_fn(|k| (c * u[k]).clamp(lo, hi));
    if max_aspect.is_finite() {
        let floor = l.iter().cloned().fold(0.0, f64::max) / (max_aspect * max_aspect);
        for x in l.iter_mut() {
            *x = x.max(floor);
        }
    }
    l
}

/// Complexity of the field at global scale `c`.
fn complexity_at<const D: usize>(eigs: &[Eig<D>], vol: &[f64], c: f64, lo: f64, hi: f64, ar: f64, nt: usize) -> f64 {
    par::sum_range(eigs.len(), par::CHUNK, nt, |i| {
        let l = scaled_eigs(c, &eigs[i].u, lo, hi, ar);
        vol[i] * l.iter().product::<f64>().sqrt()
    })
}

/// Hessian to metric with optimal L^p scaling and complexity targeting.
///
/// For a field `H` of (recovered) Hessians, the metric minimising the L^p linear-interpolation
/// error for a given number of elements is `M = c * det(|H|)^{-1/(2p+D)} |H|`, where `|H|` has the
/// absolute values of the eigenvalues of `H` and `c` is a global constant chosen so that the
/// complexity `sum_i vol_i sqrt(det M_i)` equals the target. Sizes are clamped to
/// `[hmin, hmax]` (and the aspect ratio to `max_aspect`) after scaling; because clamping changes
/// the complexity, `c` is found by bisection on the clamped complexity, which is monotone in `c`.
///
/// Fused kernel: one parallel pass does the eigen-decomposition, the absolute value, the
/// determinant factor and the first complexity partial sum, and keeps the eigen-data; the scale
/// search then runs on the stored eigenvalues only (no further decompositions); a final pass
/// rebuilds the matrices. Results are bit-identical for any thread count.
pub fn lp_metrics<const D: usize>(hess: &[Sym<D>], vol: &[f64], opt: &LpOptions) -> LpResult<D> {
    lp_metrics_with_threads(hess, vol, opt, par::threads())
}

/// [`lp_metrics`] with an explicit thread count.
pub fn lp_metrics_with_threads<const D: usize>(
    hess: &[Sym<D>],
    vol: &[f64],
    opt: &LpOptions,
    nt: usize,
) -> LpResult<D> {
    assert_eq!(hess.len(), vol.len());
    let n = hess.len();
    let expo = p_exponent(opt.p, D);
    let (lo, hi) = (1.0 / (opt.hmax * opt.hmax), 1.0 / (opt.hmin * opt.hmin));
    let ar = opt.max_aspect;

    // Eigenvalues of |H| are raised to a floor that is global, `eps_rel` times the largest Hessian
    // norm: a per-vertex floor would turn an exactly flat region into an arbitrary needle metric.
    let gmax = hess.iter().map(|h| h.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()).fold(0.0f64, f64::max);
    let floor = (opt.eps_rel * gmax).max(1e-300);

    // pass 1: eigen-decompose, |lambda| with floor, determinant factor, unscaled complexity
    let mut eigs = vec![Eig::<D> { u: [0.0; D], q: [[0.0; D]; D] }; n];
    let c1 = par::fold_chunks_mut(&mut eigs, par::CHUNK, nt, |start, chunk| {
        let mut part = 0.0;
        for (k, e) in chunk.iter_mut().enumerate() {
            let i = start + k;
            let (vals, q) = sym_eig(&hess[i]);
            let mut u = vals.map(|v| v.abs().max(floor));
            if expo != 0.0 {
                // det^expo in log space: det of a near-flat Hessian underflows
                let f = (expo * u.iter().map(|x| x.ln()).sum::<f64>()).exp();
                for x in u.iter_mut() {
                    *x *= f;
                }
            }
            *e = Eig { u, q };
            part += vol[i] * u.iter().product::<f64>().sqrt();
        }
        part
    });

    // scale search on the stored eigenvalues
    let target = opt.target_complexity;
    let floor_c: f64 = vol.iter().sum::<f64>() * lo.powf(D as f64 / 2.0);
    let ceil_c: f64 = vol.iter().sum::<f64>() * hi.powf(D as f64 / 2.0);
    let mut saturated = false;
    let c = if n == 0 {
        1.0
    } else if target <= floor_c {
        saturated = true;
        0.0
    } else if target >= ceil_c {
        saturated = true;
        f64::MAX.sqrt()
    } else {
        let c0 = (target / c1).powf(2.0 / D as f64);
        let f = |x: f64| (complexity_at(&eigs, vol, x.exp(), lo, hi, ar, nt) / target).ln();
        // g(x) = ln C(e^x) - ln target is monotone and piecewise linear in x (a power law between
        // clamp events), so a bracketing secant (Illinois) converges in a handful of evaluations.
        let (mut xa, mut xb) = (c0.ln(), c0.ln());
        let (mut fa, mut fb) = (f(xa), f(xb));
        let mut guard = 0;
        while fa > 0.0 && guard < 400 {
            xb = xa;
            fb = fa;
            xa -= 1.4;
            fa = f(xa);
            guard += 1;
        }
        guard = 0;
        while fb < 0.0 && guard < 400 {
            xa = xb;
            fa = fb;
            xb += 1.4;
            fb = f(xb);
            guard += 1;
        }
        let tol = 1e-11;
        let mut x = if fa.abs() < fb.abs() { xa } else { xb };
        let mut side = 0i32;
        for _ in 0..200 {
            if fa.abs() <= tol {
                x = xa;
                break;
            }
            if fb.abs() <= tol {
                x = xb;
                break;
            }
            // regula falsi step, safeguarded by bisection when the bracket is nearly flat
            let mut xm = (xa * fb - xb * fa) / (fb - fa);
            if !(xm > xa && xm < xb) {
                xm = 0.5 * (xa + xb);
            }
            let fm = f(xm);
            x = xm;
            if fm.abs() <= tol {
                break;
            }
            if fm < 0.0 {
                xa = xm;
                fa = fm;
                if side == -1 {
                    fb *= 0.5;
                }
                side = -1;
            } else {
                xb = xm;
                fb = fm;
                if side == 1 {
                    fa *= 0.5;
                }
                side = 1;
            }
            if xb - xa < 1e-15 * (1.0 + xa.abs()) {
                break;
            }
        }
        x.exp()
    };

    // pass 2: rebuild matrices
    let mut metrics = vec![[[0.0; D]; D]; n];
    par::fold_chunks_mut(&mut metrics, par::CHUNK, nt, |start, chunk| {
        for (k, m) in chunk.iter_mut().enumerate() {
            let e = &eigs[start + k];
            *m = from_eig(&scaled_eigs(c, &e.u, lo, hi, ar), &e.q);
        }
        0.0
    });
    let cx = complexity_at(&eigs, vol, c, lo, hi, ar, nt);
    LpResult { metrics, scale: c, complexity: cx, saturated }
}

/// Straightforward multi-pass version of [`lp_metrics`] for benchmarking and tests.
///
/// Each step is its own loop over the vertices with intermediate matrices stored: (1) `|H|`,
/// (2) determinants, (3) scale by the determinant power, (4) complexity, (5) global scaling to the
/// target, (6) clamping by a second eigen-decomposition. The clamp is applied once, so unlike the
/// fused kernel it does not re-normalise afterwards and misses the target when clamping is active.
pub fn lp_metrics_unfused<const D: usize>(hess: &[Sym<D>], vol: &[f64], opt: &LpOptions) -> Vec<Sym<D>> {
    let expo = p_exponent(opt.p, D);
    let (lo, hi) = (1.0 / (opt.hmax * opt.hmax), 1.0 / (opt.hmin * opt.hmin));
    let gmax = hess.iter().map(|h| h.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()).fold(0.0f64, f64::max);
    let floor = (opt.eps_rel * gmax).max(1e-300);
    let abs_h: Vec<Sym<D>> = hess
        .iter()
        .map(|h| {
            let (vals, q) = sym_eig(h);
            from_eig(&vals.map(|v| v.abs().max(floor)), &q)
        })
        .collect();
    let dets: Vec<f64> = abs_h.iter().map(det).collect();
    let scaled: Vec<Sym<D>> = abs_h
        .iter()
        .zip(&dets)
        .map(|(m, d)| {
            let f = if expo != 0.0 { d.ln().mul_add(expo, 0.0).exp() } else { 1.0 };
            let mut r = *m;
            for row in r.iter_mut() {
                for x in row.iter_mut() {
                    *x *= f;
                }
            }
            r
        })
        .collect();
    let c1 = complexity(&scaled, vol);
    let c = (opt.target_complexity / c1).powf(2.0 / D as f64);
    scaled
        .iter()
        .map(|m| {
            let mut r = *m;
            for row in r.iter_mut() {
                for x in row.iter_mut() {
                    *x *= c;
                }
            }
            map_eig(&r, |l| l.clamp(lo, hi))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small deterministic generator for test data.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }
        fn range(&mut self, a: f64, b: f64) -> f64 {
            a + (b - a) * self.next()
        }
        fn sym<const D: usize>(&mut self) -> Sym<D> {
            let mut a = [[0.0; D]; D];
            for i in 0..D {
                for j in i..D {
                    let v = self.range(-3.0, 3.0);
                    a[i][j] = v;
                    a[j][i] = v;
                }
            }
            a
        }
        fn spd<const D: usize>(&mut self) -> Sym<D> {
            let mut q = [[0.0; D]; D];
            let a = self.sym::<D>();
            let (_, v) = sym_eig(&a);
            for i in 0..D {
                for j in 0..D {
                    q[i][j] = v[i][j];
                }
            }
            let vals: [f64; D] = std::array::from_fn(|_| self.range(-3.0, 3.0).exp());
            from_eig(&vals, &q)
        }
        fn vec<const D: usize>(&mut self) -> [f64; D] {
            std::array::from_fn(|_| self.range(-1.0, 1.0))
        }
    }

    fn max_diff<const D: usize>(a: &Sym<D>, b: &Sym<D>) -> f64 {
        let mut m: f64 = 0.0;
        for i in 0..D {
            for j in 0..D {
                m = m.max((a[i][j] - b[i][j]).abs());
            }
        }
        m
    }

    fn eig_roundtrip<const D: usize>() {
        let mut g = Lcg(1);
        for _ in 0..500 {
            let a = g.sym::<D>();
            let (vals, q) = sym_eig(&a);
            assert!(max_diff(&from_eig(&vals, &q), &a) < 1e-12);
            let qtq = mul(&transpose(&q), &q);
            assert!(max_diff(&qtq, &identity::<D>()) < 1e-13);
            assert!(vals.windows(2).all(|w| w[0] >= w[1]));
        }
    }

    #[test]
    fn eigen_decomposition_2d_3d() {
        eig_roundtrip::<2>();
        eig_roundtrip::<3>();
    }

    #[test]
    fn eigen_handles_repeated_and_diagonal() {
        let a = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]];
        let (v, _) = sym_eig(&a);
        assert_eq!(v, [5.0, 2.0, 2.0]);
        let (v, _) = sym_eig(&isotropic::<2>(0.5));
        assert!((v[0] - 4.0).abs() < 1e-15 && (v[1] - 4.0).abs() < 1e-15);
    }

    #[test]
    fn inverse_and_det() {
        let mut g = Lcg(2);
        for _ in 0..100 {
            let a = g.spd::<3>();
            assert!(max_diff(&mul(&a, &inv(&a)), &identity::<3>()) < 1e-11);
            let (v, _) = sym_eig(&a);
            assert!((det(&a) - v.iter().product::<f64>()).abs() < 1e-10 * det(&a));
            let b = g.spd::<2>();
            assert!(max_diff(&mul(&b, &inv(&b)), &identity::<2>()) < 1e-12);
        }
    }

    fn intersection_props<const D: usize>() {
        let mut g = Lcg(3);
        for _ in 0..300 {
            let (a, b) = (g.spd::<D>(), g.spd::<D>());
            let m = intersect(&a, &b);
            // result dominates both (its unit ball is inside both)
            for _ in 0..40 {
                let e = g.vec::<D>();
                let q = quad_form(&m, &e);
                assert!(q >= quad_form(&a, &e) * (1.0 - 1e-10));
                assert!(q >= quad_form(&b, &e) * (1.0 - 1e-10));
            }
            // symmetric, idempotent
            assert!(
                max_diff(&m, &intersect(&b, &a))
                    < 1e-9 * (1.0 + m.iter().flatten().fold(0.0, |s: f64, x| s.max(x.abs())))
            );
            assert!(max_diff(&intersect(&a, &a), &a) < 1e-10);
            // commutes with a congruence: sizes are what matters, not the basis
            // minimality: there is a direction where it equals a or b (a simultaneous eigenvector)
            let (ln, _) = {
                let (l1, q1) = sym_eig(&a);
                let s = from_eig(&l1.map(|x| 1.0 / x.sqrt()), &q1);
                sym_eig(&symmetrize(&mul(&mul(&s, &b), &s)))
            };
            let dm = det(&m);
            let expected = det(&a) * ln.iter().map(|x| x.max(1.0)).product::<f64>();
            assert!((dm - expected).abs() < 1e-8 * expected);
        }
    }

    #[test]
    fn intersection_is_the_simultaneous_reduction() {
        intersection_props::<2>();
        intersection_props::<3>();
    }

    #[test]
    fn isotropic_intersection_takes_smaller_size() {
        let m = intersect(&isotropic::<3>(0.5), &isotropic::<3>(0.2));
        assert!(max_diff(&m, &isotropic::<3>(0.2)) < 1e-12);
    }

    #[test]
    fn anisotropic_intersection_example() {
        // sizes (1, 0.1) along x,y intersected with (0.1, 1): both directions at 0.1
        let a = from_sizes(&[1.0, 0.1], &identity::<2>());
        let b = from_sizes(&[0.1, 1.0], &identity::<2>());
        assert!(max_diff(&intersect(&a, &b), &isotropic::<2>(0.1)) < 1e-9);
    }

    #[test]
    fn edge_length_matches_integral() {
        let mut g = Lcg(4);
        for _ in 0..50 {
            let (a, b) = (g.spd::<2>(), g.spd::<2>());
            let e = [g.range(0.1, 1.0), g.range(-1.0, 1.0)];
            let (la, lb) = (length(&a, &e), length(&b, &e));
            // numeric integral of |e| / h(t) with h(t) = ha^(1-t) hb^t, i.e. sizes interpolated
            // geometrically, which is the model of edge_length
            let n = 20000;
            let mut s = 0.0;
            for k in 0..n {
                let t = (k as f64 + 0.5) / n as f64;
                s += la.powf(1.0 - t) * lb.powf(t);
            }
            s /= n as f64;
            assert!((edge_length(&a, &b, &e) - s).abs() < 1e-8 * s);
        }
        let a = isotropic::<2>(0.3);
        assert!((edge_length(&a, &a, &[0.3, 0.0]) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn log_euclid_properties() {
        let mut g = Lcg(5);
        let a = g.spd::<3>();
        assert!(max_diff(&log_euclid_mean(&[(1.0, a)]), &a) < 1e-9);
        // determinant of the geometric mean is the geometric mean of determinants
        let b = g.spd::<3>();
        let m = log_euclid_mean(&[(0.5, a), (0.5, b)]);
        assert!((det(&m) - (det(&a) * det(&b)).sqrt()).abs() < 1e-9 * det(&m));
        // isotropic: geometric mean of sizes
        let m = log_euclid_mean(&[(0.25, isotropic::<2>(1.0)), (0.75, isotropic::<2>(0.1))]);
        let h = 1.0f64.powf(0.25) * 0.1f64.powf(0.75);
        assert!(max_diff(&m, &isotropic::<2>(h)) < 1e-12);
    }

    #[test]
    fn grow_is_linear_in_size() {
        let m = from_sizes(&[0.1, 0.5], &identity::<2>());
        let gm = grow(&m, 2.0, 0.25);
        let (h, _) = sizes(&gm);
        assert!((h[0] - 0.6).abs() < 1e-12 && (h[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn gradation_bounds_size_growth_on_a_path() {
        // 1D chain along x, tiny size at one end, huge elsewhere
        let n = 60;
        let pts: Vec<[f64; 2]> = (0..n).map(|i| [i as f64 * 0.05, 0.0]).collect();
        let edges: Vec<[u32; 2]> = (0..n - 1).map(|i| [i as u32, i as u32 + 1]).collect();
        let mut ms: Vec<Sym<2>> = vec![isotropic(5.0); n];
        ms[0] = isotropic(0.01);
        let g = 0.3;
        let passes = grade(&mut ms, &pts, &edges, g, 200);
        assert!(passes < 200);
        for e in &edges {
            let (i, j) = (e[0] as usize, e[1] as usize);
            let d = 0.05;
            let (hi, _) = sizes(&ms[i]);
            let (hj, _) = sizes(&ms[j]);
            for k in 0..2 {
                assert!(hj[k] <= hi[k] + g * d + 1e-9);
                assert!(hi[k] <= hj[k] + g * d + 1e-9);
            }
        }
        // refinement is never lost
        assert!((sizes(&ms[0]).0[0] - 0.01).abs() < 1e-9);
    }

    #[test]
    fn anisotropic_gradation_keeps_direction_information() {
        let pts = [[0.0, 0.0], [1.0, 0.0]];
        let mut ms = vec![from_sizes(&[0.01, 1.0], &identity::<2>()), isotropic::<2>(1.0)];
        grade(&mut ms, &pts, &[[0, 1]], 0.5, 50);
        // at x = 1 the x-size can be at most 0.01 + 0.5 = 0.51 but the y-size stays free
        let (h, _) = sizes(&ms[1]);
        assert!((h[0] - 0.51).abs() < 1e-9, "{h:?}"); // smaller size first
        assert!((h[1] - 1.0).abs() < 1e-9);
    }

    fn const_field<const D: usize>(h: &Sym<D>, n: usize, v: f64) -> (Vec<Sym<D>>, Vec<f64>) {
        (vec![*h; n], vec![v; n])
    }

    fn lp_hits_target<const D: usize>(hess: Sym<D>) {
        for p in [1.0, 2.0, f64::INFINITY] {
            let (hs, vol) = const_field(&hess, 100, 0.01);
            let opt = LpOptions { p, target_complexity: 1234.5, ..Default::default() };
            let r = lp_metrics_with_threads(&hs, &vol, &opt, 1);
            assert!(!r.saturated);
            assert!((r.complexity - 1234.5).abs() < 1e-9 * 1234.5, "p={p}: {}", r.complexity);
            assert!((complexity(&r.metrics, &vol) - 1234.5).abs() < 1e-8 * 1234.5);
            // metric proportional to |H|: eigen-ratios equal the Hessian eigenvalue ratios
            let (hv, _) = sym_eig(&hess);
            let (mv, _) = sym_eig(&r.metrics[0]);
            for k in 1..D {
                let a = hv[0].abs() / hv[k].abs();
                let b = mv[0] / mv[k];
                assert!((a - b).abs() < 1e-8 * a.max(1.0));
            }
        }
    }

    #[test]
    fn lp_scaling_hits_target_complexity() {
        lp_hits_target::<2>([[100.0, 0.0], [0.0, 1.0]]);
        lp_hits_target::<3>([[400.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, -1.0]]);
    }

    #[test]
    fn lp_exponent_matches_formula() {
        // 2D, p = 2: M = c * det(|H|)^(-1/6) |H|; two Hessians differing by a uniform factor s
        // must give eigenvalue ratio s^(1 - 2/6) = s^(2/3) before the global normalisation
        let h1: Sym<2> = [[4.0, 0.0], [0.0, 1.0]];
        let h2: Sym<2> = [[4.0 * 8.0, 0.0], [0.0, 8.0]];
        let hs = vec![h1, h2];
        let vol = vec![1.0, 1.0];
        let opt = LpOptions { p: 2.0, target_complexity: 50.0, ..Default::default() };
        let r = lp_metrics_with_threads(&hs, &vol, &opt, 1);
        let ratio = r.metrics[1][0][0] / r.metrics[0][0][0];
        assert!((ratio - 8f64.powf(1.0 - 2.0 / 6.0)).abs() < 1e-10, "{ratio}");
    }

    #[test]
    fn lp_clamps_and_still_hits_target() {
        let mut g = Lcg(6);
        let n = 2000;
        let hs: Vec<Sym<2>> =
            (0..n).map(|_| g.sym::<2>().map(|r| r.map(|x| x * 10f64.powf(g.range(-3.0, 3.0))))).collect();
        let hs: Vec<Sym<2>> = hs.iter().map(symmetrize).collect();
        let vol = vec![1.0 / n as f64; n];
        let opt = LpOptions {
            p: 2.0,
            target_complexity: 400.0,
            hmin: 0.01,
            hmax: 0.5,
            max_aspect: 30.0,
            ..Default::default()
        };
        let r = lp_metrics_with_threads(&hs, &vol, &opt, 1);
        assert!(!r.saturated);
        assert!((r.complexity - 400.0).abs() < 1e-8 * 400.0);
        for m in &r.metrics {
            let (h, _) = sizes(m);
            assert!(h[0] >= 0.01 * (1.0 - 1e-9) && h[1] <= 0.5 * (1.0 + 1e-9));
            assert!(h[1] / h[0] <= 30.0 * (1.0 + 1e-8));
        }
    }

    #[test]
    fn lp_reports_saturation() {
        let hs = vec![[[1.0, 0.0], [0.0, 1.0]]; 10];
        let vol = vec![0.1; 10];
        let opt = LpOptions { target_complexity: 1e9, hmin: 0.1, hmax: 1.0, ..Default::default() };
        let r = lp_metrics_with_threads(&hs, &vol, &opt, 1);
        assert!(r.saturated);
        assert!(r.metrics.iter().all(|m| (sizes(m).0[0] - 0.1).abs() < 1e-9));
    }

    #[test]
    fn fused_equals_unfused_without_clamping_and_is_thread_independent() {
        let mut g = Lcg(7);
        let n = 5000;
        let hs: Vec<Sym<3>> = (0..n).map(|_| g.spd::<3>()).collect();
        let vol: Vec<f64> = (0..n).map(|_| g.range(0.5, 1.5) / n as f64).collect();
        let opt = LpOptions { p: 2.0, target_complexity: 2000.0, ..Default::default() };
        let r1 = lp_metrics_with_threads(&hs, &vol, &opt, 1);
        let un = lp_metrics_unfused(&hs, &vol, &opt);
        for (a, b) in r1.metrics.iter().zip(&un) {
            assert!(max_diff(a, b) < 1e-9 * (1.0 + a[0][0].abs()), "{a:?} {b:?}");
        }
        for t in [2, 3, 8] {
            let rt = lp_metrics_with_threads(&hs, &vol, &opt, t);
            assert_eq!(r1.scale.to_bits(), rt.scale.to_bits());
            for (a, b) in r1.metrics.iter().zip(&rt.metrics) {
                for i in 0..3 {
                    for j in 0..3 {
                        assert_eq!(a[i][j].to_bits(), b[i][j].to_bits());
                    }
                }
            }
        }
    }

    #[test]
    fn unfused_misses_target_when_clamping_is_active_fused_does_not() {
        let mut g = Lcg(8);
        let n = 3000;
        let hs: Vec<Sym<2>> = (0..n)
            .map(|_| {
                let a = g.spd::<2>();
                a.map(|r| r.map(|x| x * 1e3))
            })
            .collect();
        let vol = vec![1.0 / n as f64; n];
        let opt = LpOptions { p: 2.0, target_complexity: 1000.0, hmin: 0.02, hmax: 0.5, ..Default::default() };
        let fused = lp_metrics_with_threads(&hs, &vol, &opt, 1);
        let un = lp_metrics_unfused(&hs, &vol, &opt);
        assert!((fused.complexity - 1000.0).abs() < 1e-8 * 1000.0);
        assert!((complexity(&un, &vol) - 1000.0).abs() > 1e-3 * 1000.0);
    }

    #[test]
    fn exactly_flat_regions_give_finite_metrics() {
        // half of the field has an exactly zero Hessian (as for a saturated tanh layer)
        let mut hs: Vec<Sym<2>> = vec![[[0.0; 2]; 2]; 50];
        hs.extend(vec![[[400.0, 0.0], [0.0, 0.0]]; 50]);
        let vol = vec![0.01; 100];
        for p in [1.0, 2.0, f64::INFINITY] {
            let opt = LpOptions { p, target_complexity: 100.0, hmin: 1e-3, hmax: 1.0, ..Default::default() };
            let r = lp_metrics_with_threads(&hs, &vol, &opt, 1);
            assert!(r.metrics.iter().flatten().flatten().all(|x| x.is_finite()), "p={p}");
            assert!((r.complexity - 100.0).abs() < 1e-8 * 100.0, "p={p}: {}", r.complexity);
            // the curved half is finer than the flat half
            assert!(sizes(&r.metrics[75]).0[0] < sizes(&r.metrics[10]).0[0]);
            let un = lp_metrics_unfused(&hs, &vol, &opt);
            assert!(un.iter().flatten().flatten().all(|x| x.is_finite()), "unfused p={p}");
        }
    }

    #[test]
    fn element_count_conversion() {
        for d in [2, 3] {
            let n = 12345.0;
            assert!((elements_for_complexity(complexity_for_elements(n, d), d) - n).abs() < 1e-9);
        }
    }
}
