//! A posteriori error estimator for time-harmonic Maxwell fields on tetrahedra (ADR 53), and the metric it asks for.
//!
//! The recovered-Hessian estimator of `estimate3` looks at plotted nodal fields. For edge elements of the lowest order it does not beat a
//! uniform mesh (`docs/ports/meshers/MISSING_FEATURES.md`, section 2: 1.2e-2 on 1490 tetrahedra against 5.3e-3 on 696 for a cavity mode)
//! because the field is only piecewise linear and tangentially continuous, so its Hessian is noise. The residual of the equation is what
//! carries the information. This is the estimator of Beck, Hiptmair, Hoppe and Wohlmuth (*Residual based a posteriori error estimators for
//! eddy current computation*, M2AN 34, 2000), written for the time-harmonic equation
//!
//! ```text
//!     curl (1/mu_r) curl E  -  k0^2 eps_r E  =  0           (source free; k0 = omega / c)
//! ```
//!
//! with, for every tetrahedron `T` of diameter `h_T` and every face `F` of diameter `h_F`,
//!
//! ```text
//!     eta_T^2 = h_T^2 || curl (1/mu_r) curl E - k0^2 eps_r E ||_T^2          interior residual
//!             + h_T^2 || div (eps_r E) ||_T^2                                 divergence of the flux density
//!             + 1/2 sum_F h_F ( || n x [ (1/mu_r) curl E ] ||_F^2             jump of the tangential magnetic field
//!                             + || [ eps_r E . n ] ||_F^2 )                    jump of the normal flux density
//! ```
//!
//! The tangential trace of `E` is continuous for edge elements, so it gives no jump. A face on the boundary contributes nothing when the
//! tangential field is prescribed (a perfect electric conductor) and `h_F || n x (1/mu_r) curl E ||^2` when the boundary condition is
//! natural (a perfect magnetic conductor or a symmetry plane). For the lowest order `curl curl E` vanishes inside a tetrahedron and the
//! residual is `k0^2 eps_r E` itself, which is why this estimator sees what the Hessian does not.
//!
//! **What it needs of the field.** Only the trait [`CurlField`]: the field and its curl at a barycentric point of a cell. The physics
//! crates implement it for their degrees of freedom (the edge and face functions of `ems3d`); nothing here knows how they are numbered.
//! The field must be a polynomial of degree at most 2 in the cell (orders 1 and 2), and `sample` must accept barycentric points a little
//! outside the cell (the polynomial continued), because first derivatives are taken by central differences of exact order for quadratics.
//!
//! **Not covered.** Media with complex or tensor coefficients (absorbing layers: leave their cells out with `MaxwellOptions::skip`), a
//! source term, and a driven problem; the estimator is for the source-free equation, eigenmodes and homogeneous scattering problems.

use crate::adapt::metric::{isotropic, intersect, Sym};
use crate::fem::quad_rule;
use crate::linalg::csr::C;
use crate::mesh::SimplexMesh;
use rayon::prelude::*;
use std::collections::HashMap;

/// The field of one solution, as the estimator reads it.
pub trait CurlField: Sync {
    /// `E` and `curl E` at the barycentric point `lam` of cell `cell` (complex amplitudes, physical units of the unknowns, curl taken in
    /// the physical cell).
    fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]);
}

/// The medium and the boundary conditions of the problem the field solves.
pub struct MaxwellOptions<'a> {
    /// squared free-space wave number `(omega / c)^2` of the solution (the real part of the eigenvalue for a lossy cavity)
    pub k0_sq: f64,
    /// relative permittivity of every cell (isotropic, real)
    pub eps_r: &'a [f64],
    /// relative permeability of every cell (isotropic, real)
    pub mu_r: &'a [f64],
    /// cells of these domain labels are left out (absorbing layers: complex stretched coefficients); their indicator is zero
    pub skip: &'a (dyn Fn(i32) -> bool + Sync),
    /// boundary labels with a natural condition (`n x curl E / mu = 0`: a magnetic wall or a symmetry plane). Every other boundary facet
    /// is taken as a perfect electric conductor and contributes nothing.
    pub natural: &'a (dyn Fn(i32) -> bool + Sync),
}

/// The indicators of one field.
#[derive(Clone, Debug)]
pub struct Indicators {
    /// `eta_T^2` of every cell
    pub eta2: Vec<f64>,
    /// the four parts, summed over the mesh: interior residual, divergence, jump of the tangential magnetic field, jump of the normal
    /// flux density (so a test can say which one a change of the mesh lowers)
    pub parts: [f64; 4],
    /// `sqrt(sum eta_T^2)`
    pub total: f64,
}

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Gradients of the four barycentric functions and the volume of a cell.
fn cell_geometry(p: &[V3; 4]) -> ([V3; 4], f64) {
    let (e1, e2, e3) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    let six_v = dot(e1, cross(e2, e3));
    let g1 = cross(e2, e3).map(|x| x / six_v);
    let g2 = cross(e3, e1).map(|x| x / six_v);
    let g3 = cross(e1, e2).map(|x| x / six_v);
    let g0 = [-(g1[0] + g2[0] + g3[0]), -(g1[1] + g2[1] + g3[1]), -(g1[2] + g2[2] + g3[2])];
    ([g0, g1, g2, g3], six_v.abs() / 6.0)
}

fn diameter(p: &[V3]) -> f64 {
    let mut d: f64 = 0.0;
    for i in 0..p.len() {
        for j in i + 1..p.len() {
            d = d.max(norm(sub(p[i], p[j])));
        }
    }
    d
}

fn cnorm2(v: [C; 3]) -> f64 {
    v.iter().map(|z| z.norm_sqr()).sum()
}

/// Real vector times complex vector: `n x v` for a real `n`.
fn ccross(n: V3, v: [C; 3]) -> [C; 3] {
    [v[2] * n[1] - v[1] * n[2], v[0] * n[2] - v[2] * n[0], v[1] * n[0] - v[0] * n[1]]
}

fn cdot(n: V3, v: [C; 3]) -> C {
    v[0] * n[0] + v[1] * n[1] + v[2] * n[2]
}

/// Gradient `d F_k / d x_j` of a field of the cell at `lam`, by central differences along the directions `e_v - (1/4, 1/4, 1/4, 1/4)`:
/// `grad = sum_v (dF/dD_v) grad(lam_v)`. Exact for polynomials of degree 2.
fn gradient(f: &dyn Fn([f64; 4]) -> [C; 3], lam: [f64; 4], g: &[V3; 4]) -> [[C; 3]; 3] {
    const T: f64 = 0.05;
    let mut out = [[C::new(0.0, 0.0); 3]; 3];
    for v in 0..4 {
        let (mut plus, mut minus) = (lam, lam);
        for w in 0..4 {
            let d = if w == v { 0.75 } else { -0.25 };
            plus[w] += T * d;
            minus[w] -= T * d;
        }
        let (fp, fm) = (f(plus), f(minus));
        for k in 0..3 {
            let dk = (fp[k] - fm[k]) / (2.0 * T);
            for j in 0..3 {
                out[k][j] += dk * g[v][j];
            }
        }
    }
    out
}

/// `curl F` from the gradient tensor `grad[k][j] = d F_k / d x_j`.
fn curl_of(grad: &[[C; 3]; 3]) -> [C; 3] {
    [grad[2][1] - grad[1][2], grad[0][2] - grad[2][0], grad[1][0] - grad[0][1]]
}

fn div_of(grad: &[[C; 3]; 3]) -> C {
    grad[0][0] + grad[1][1] + grad[2][2]
}

/// The indicators `eta_T^2` of `field` on `mesh`.
pub fn maxwell_indicators(mesh: &SimplexMesh, field: &dyn CurlField, opt: &MaxwellOptions) -> Result<Indicators, String> {
    if mesh.dim != 3 {
        return Err("the Maxwell estimator works on tetrahedral meshes".into());
    }
    let nc = mesh.n_cells();
    if opt.eps_r.len() != nc || opt.mu_r.len() != nc {
        return Err(format!("{} permittivities and {} permeabilities for {nc} cells", opt.eps_r.len(), opt.mu_r.len()));
    }
    if !(opt.eps_r.iter().chain(opt.mu_r).all(|x| *x > 0.0 && x.is_finite())) {
        return Err("permittivity and permeability must be positive and finite".into());
    }
    let pt = |v: u32| -> V3 { [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]] };
    let cell_pts = |c: usize| -> [V3; 4] {
        let v = mesh.cell(c);
        [pt(v[0]), pt(v[1]), pt(v[2]), pt(v[3])]
    };
    let vol_rule = quad_rule(3, 4);
    let face_rule = quad_rule(2, 4);

    // interior terms, cell by cell
    let interior: Vec<[f64; 2]> = (0..nc)
        .into_par_iter()
        .map(|c| {
            if (opt.skip)(mesh.cell_labels[c]) {
                return [0.0, 0.0];
            }
            let p = cell_pts(c);
            let (g, vol) = cell_geometry(&p);
            let h = diameter(&p);
            let (eps, mu) = (opt.eps_r[c], opt.mu_r[c]);
            // curl curl E is constant in the cell for the orders 1 and 2 (curl E is affine): one gradient of curl E at the centroid
            let centre = [0.25; 4];
            let curl_fn = |lam: [f64; 4]| field.sample(c, lam).1;
            let cc = curl_of(&gradient(&curl_fn, centre, &g)).map(|z| z / mu);
            let (mut res, mut div) = (0.0, 0.0);
            let e_fn = |lam: [f64; 4]| field.sample(c, lam).0;
            for (lam, w) in vol_rule.lam.iter().zip(&vol_rule.w) {
                let l4 = [lam[0], lam[1], lam[2], lam[3]];
                let e = e_fn(l4);
                let r = [cc[0] - e[0] * (opt.k0_sq * eps), cc[1] - e[1] * (opt.k0_sq * eps), cc[2] - e[2] * (opt.k0_sq * eps)];
                res += w * vol * cnorm2(r);
                let d = div_of(&gradient(&e_fn, l4, &g)) * eps;
                div += w * vol * d.norm_sqr();
            }
            [h * h * res, h * h * div]
        })
        .collect();

    // faces: the cells on each side, in sorted-vertex order
    let mut faces: HashMap<[u32; 3], Vec<(usize, usize)>> = HashMap::new();
    for c in 0..nc {
        let v = mesh.cell(c);
        for skip in 0..4 {
            let mut k: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            k.sort_unstable();
            faces.entry([k[0], k[1], k[2]]).or_default().push((c, skip));
        }
    }
    let mut keys: Vec<[u32; 3]> = faces.keys().copied().collect();
    keys.sort_unstable();
    // boundary facet labels by face
    let mut facet_label: HashMap<[u32; 3], i32> = HashMap::new();
    for f in 0..mesh.n_facets() {
        let v = mesh.facet(f);
        let mut k = [v[0], v[1], v[2]];
        k.sort_unstable();
        facet_label.insert(k, mesh.facet_labels[f]);
    }

    // jump terms: (cell, tangential part, normal part)
    let face_terms: Vec<Vec<(usize, f64, f64)>> = keys
        .par_iter()
        .map(|k| {
            let cs = &faces[k];
            let (a, b, c) = (pt(k[0]), pt(k[1]), pt(k[2]));
            let area_vec = cross(sub(b, a), sub(c, a));
            let area = 0.5 * norm(area_vec);
            let h_f = diameter(&[a, b, c]);
            let nrm = area_vec.map(|x| x / (2.0 * area));
            // barycentric point of the face point (l0 l1 l2) in a cell
            let in_cell = |cell: usize, l: &[f64]| -> [f64; 4] {
                let v = mesh.cell(cell);
                let mut lam = [0.0; 4];
                for (q, &fv) in k.iter().enumerate() {
                    let loc = v.iter().position(|&x| x == fv).unwrap();
                    lam[loc] = l[q];
                }
                lam
            };
            // outward normal of the first cell
            let outward = |cell: usize, skip: usize| -> V3 {
                let opp = pt(mesh.cell(cell)[skip]);
                if dot(nrm, sub(opp, a)) > 0.0 {
                    nrm.map(|x| -x)
                } else {
                    nrm
                }
            };
            let live = |cell: usize| !(opt.skip)(mesh.cell_labels[cell]);
            let mut out = vec![];
            if cs.len() == 2 {
                let ((c1, s1), (c2, _)) = (cs[0], cs[1]);
                if !live(c1) && !live(c2) {
                    return out;
                }
                let n = outward(c1, s1);
                let (mut jt, mut jn) = (0.0, 0.0);
                for (l, w) in face_rule.lam.iter().zip(&face_rule.w) {
                    let (e1, h1) = field.sample(c1, in_cell(c1, l));
                    let (e2, h2) = field.sample(c2, in_cell(c2, l));
                    let hh = [0, 1, 2].map(|i| h1[i] / opt.mu_r[c1] - h2[i] / opt.mu_r[c2]);
                    jt += w * area * cnorm2(ccross(n, hh));
                    let dn = cdot(n, e1) * opt.eps_r[c1] - cdot(n, e2) * opt.eps_r[c2];
                    jn += w * area * dn.norm_sqr();
                }
                // a face next to a skipped cell sees only its own side's contribution
                for (cell, other) in [(c1, c2), (c2, c1)] {
                    if live(cell) && live(other) {
                        out.push((cell, 0.5 * h_f * jt, 0.5 * h_f * jn));
                    }
                }
            } else if cs.len() == 1 {
                let (c1, s1) = cs[0];
                let natural = facet_label.get(k).is_some_and(|&l| (opt.natural)(l));
                if natural && live(c1) {
                    let n = outward(c1, s1);
                    let mut jt = 0.0;
                    for (l, w) in face_rule.lam.iter().zip(&face_rule.w) {
                        let (_, h1) = field.sample(c1, in_cell(c1, l));
                        let hh = h1.map(|z| z / opt.mu_r[c1]);
                        jt += w * area * cnorm2(ccross(n, hh));
                    }
                    out.push((c1, h_f * jt, 0.0));
                }
            }
            out
        })
        .collect();

    let mut eta2 = vec![0.0; nc];
    let mut parts = [0.0; 4];
    for (c, [res, div]) in interior.iter().enumerate() {
        eta2[c] += res + div;
        parts[0] += res;
        parts[1] += div;
    }
    for terms in face_terms {
        for (c, jt, jn) in terms {
            eta2[c] += jt + jn;
            parts[2] += jt;
            parts[3] += jn;
        }
    }
    let total = eta2.iter().sum::<f64>().sqrt();
    Ok(Indicators { eta2, parts, total })
}

/// Settings of [`indicator_metrics`].
#[derive(Clone, Debug)]
pub struct EquidistributeOptions {
    /// order of the edge elements (1 or 2); the indicator of a cell of size `s` falls like `s^(order + 3/2)`
    pub order: usize,
    /// largest element size asked for
    pub h_max: f64,
    /// smallest element size asked for, as a fraction of `h_max`
    pub h_min_fraction: f64,
    /// the mesh is asked to have about this many cells
    pub target_cells: f64,
}

/// The isotropic metric at every vertex that equidistributes the indicators over about `target_cells` cells. A cell of size `s` with
/// indicator `eta` asks for the size `s (tau / eta)^(1 / (order + 3/2))`, kept between `h_min_fraction * h_max` and `h_max` and within a
/// factor 4 of `s`; the level `tau` is found by bisection so that the predicted number of cells, the sum of `(s / s_new)^3`, is the target.
/// A vertex takes the smallest size asked by a cell that meets there. `floor` (a metric per vertex, for example the geometry floor) is
/// intersected when given.
pub fn indicator_metrics(mesh: &SimplexMesh, eta2: &[f64], opt: &EquidistributeOptions, floor: Option<&[Sym<3>]>) -> Result<Vec<Sym<3>>, String> {
    let nc = mesh.n_cells();
    if eta2.len() != nc {
        return Err(format!("{} indicators for {nc} cells", eta2.len()));
    }
    if !(opt.h_max > 0.0 && opt.target_cells >= 1.0 && (0.0..=1.0).contains(&opt.h_min_fraction)) {
        return Err("h_max and target_cells must be positive and h_min_fraction within 0 and 1".into());
    }
    let q = opt.order as f64 + 1.5;
    let pt = |v: u32| -> V3 { [mesh.points[3 * v as usize], mesh.points[3 * v as usize + 1], mesh.points[3 * v as usize + 2]] };
    // the size of a cell: the edge of the regular tetrahedron of its volume
    let size: Vec<f64> = (0..nc)
        .map(|c| {
            let v = mesh.cell(c);
            let (_, vol) = cell_geometry(&[pt(v[0]), pt(v[1]), pt(v[2]), pt(v[3])]);
            (6.0 * 2f64.sqrt() * vol).cbrt()
        })
        .collect();
    let h_min = opt.h_max * opt.h_min_fraction;
    let new_size = |tau: f64, c: usize| -> f64 {
        let eta = eta2[c].max(0.0).sqrt();
        let ratio = if eta > 0.0 { (tau / eta).powf(1.0 / q) } else { 4.0 };
        (size[c] * ratio.clamp(0.25, 4.0)).clamp(h_min.min(opt.h_max), opt.h_max)
    };
    let predicted = |tau: f64| -> f64 { (0..nc).map(|c| (size[c] / new_size(tau, c)).powi(3)).sum() };
    // predicted is decreasing in tau
    let max_eta = eta2.iter().cloned().fold(0.0, f64::max).sqrt().max(1e-300);
    let (mut lo, mut hi) = (max_eta * 1e-12, max_eta * 1e12);
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        if predicted(mid) > opt.target_cells {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let tau = (lo * hi).sqrt();
    let mut vh = vec![opt.h_max; mesh.n_points()];
    for c in 0..nc {
        let s = new_size(tau, c);
        for &v in mesh.cell(c) {
            vh[v as usize] = vh[v as usize].min(s);
        }
    }
    let mut out: Vec<Sym<3>> = vh.iter().map(|&h| isotropic::<3>(h.max(1e-300))).collect();
    if let Some(f) = floor {
        for (a, b) in out.iter_mut().zip(f) {
            *a = intersect(a, b);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::c;
    use crate::mesh::fixtures::box_tet;

    /// The whitney field of a plane wave sampled exactly: `E = (0, 0, e^{-ikx})` is not in the space, so use a polynomial field instead.
    struct Poly {
        k0_sq: f64,
        pts: Vec<[V3; 4]>,
        /// E = (a*y, b*z, c*x) + (d, 0, 0): curl = (0 - b, 0 - c, 0 - a)... evaluated from the coordinates
        co: [f64; 3],
    }

    impl CurlField for Poly {
        fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
            let p = &self.pts[cell];
            let x: V3 = std::array::from_fn(|k| (0..4).map(|v| lam[v] * p[v][k]).sum());
            let [a, b, cc] = self.co;
            let e = [c(a * x[1], 0.0), c(b * x[2], 0.0), c(cc * x[0], 0.0)];
            // curl (a y, b z, c x) = (d_y(c x) - d_z(b z), d_z(a y) - d_x(c x), d_x(b z) - d_y(a y)) = (-b, -c, -a)
            let curl = [c(-b, 0.0), c(-cc, 0.0), c(-a, 0.0)];
            let _ = self.k0_sq;
            (e, curl)
        }
    }

    #[test]
    fn a_field_with_known_residual_gives_the_closed_form_indicators() {
        // E linear: curl E constant, curl curl E = 0, div E = 0 for this field, no jumps (single-valued field): eta_T^2 = h^2 k0^4 eps^2 ||E||^2
        let m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [2, 2, 2]));
        let pts: Vec<[V3; 4]> = (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                std::array::from_fn(|i| [m.points[3 * v[i] as usize], m.points[3 * v[i] as usize + 1], m.points[3 * v[i] as usize + 2]])
            })
            .collect();
        let f = Poly { k0_sq: 3.0, pts: pts.clone(), co: [1.0, 2.0, -0.5] };
        let eps = vec![2.0; m.n_cells()];
        let mu = vec![1.0; m.n_cells()];
        let opt = MaxwellOptions { k0_sq: 3.0, eps_r: &eps, mu_r: &mu, skip: &|_| false, natural: &|_| false };
        let ind = maxwell_indicators(&m, &f, &opt).unwrap();
        // closed form per cell: integrate |k0^2 eps E|^2 over the cell with the exact quadrature of the test (a quadratic integrand)
        for (c, p) in pts.iter().enumerate() {
            let (_, vol) = cell_geometry(p);
            let h = diameter(p);
            let rule = quad_rule(3, 4);
            let mut s = 0.0;
            for (lam, w) in rule.lam.iter().zip(&rule.w) {
                let (e, _) = f.sample(c, [lam[0], lam[1], lam[2], lam[3]]);
                s += w * vol * cnorm2(e) * (3.0f64 * 2.0).powi(2);
            }
            assert!((ind.eta2[c] - h * h * s).abs() < 1e-10 * (1.0 + h * h * s), "cell {c}: {} vs {}", ind.eta2[c], h * h * s);
        }
        assert!(ind.parts[1].abs() < 1e-18 && ind.parts[2] == 0.0 && ind.parts[3] == 0.0);
    }

    #[test]
    fn equidistribution_hits_the_target_count_and_refines_where_the_error_is() {
        let m = crate::adapt::work::simplex_from_adapt(&box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [4, 4, 4]));
        // an error concentrated near the origin
        let eta2: Vec<f64> = (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                let x = (0..4).map(|i| m.points[3 * v[i] as usize]).sum::<f64>() / 4.0;
                (-8.0 * x).exp()
            })
            .collect();
        let opt = EquidistributeOptions { order: 1, h_max: 1.0, h_min_fraction: 0.02, target_cells: 2.0 * m.n_cells() as f64 };
        let met = indicator_metrics(&m, &eta2, &opt, None).unwrap();
        let h_at = |v: usize| (1.0 / met[v][0][0]).sqrt();
        // near the origin the size is smaller than far away
        let near = (0..m.n_points()).filter(|&v| m.points[3 * v] < 0.01).map(h_at).fold(f64::INFINITY, f64::min);
        let far = (0..m.n_points()).filter(|&v| m.points[3 * v] > 0.99).map(h_at).fold(f64::INFINITY, f64::min);
        assert!(near < far, "{near} {far}");
        assert!(indicator_metrics(&m, &eta2[1..], &opt, None).is_err());
    }
}
