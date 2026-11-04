//! Scalar Helmholtz equation with PML on mixed meshes: P1 triangles, Q1 quadrilaterals (2D) and
//! P1 tetrahedra, trilinear hexahedra, P1 x linear prisms (3D).
//!
//! Time convention `exp(-i omega t)`, wavenumber `k = omega / c`. In a PML region with stretched
//! axes the coordinate is complexified, `s_a(x) = 1 + i sigma_a(x) / k`, with the polynomial profile
//! `sigma_a = sigma_max (d / t)^m` in the distance `d` from the interface plane and
//! `sigma_max = -(m + 1) ln(R0) / (2 t)`, where `R0` is the design reflection at normal incidence
//! of the continuous layer. With `S = prod_a s_a` the weak form is
//!
//! ```text
//! sum  S/s_a^2 u_a v_a  -  k^2 S u v   =   sum  S f v
//! ```
//!
//! evaluated at the quadrature points of every element (so the profile is sampled, not averaged),
//! with homogeneous Dirichlet conditions on the outer boundary of the mesh.
//!
//! Quadrature: 3-point (edge midpoints) on triangles, 4-point on tetrahedra, 2x2 / 2x2x2 Gauss on
//! quadrilaterals / hexahedra, 3-point x 2-point on prisms. All elements go through one isoparametric
//! routine, so the only per-shape code is the reference shape functions and the rule.

use super::band::{BandLu, BandMatrix, rcm};
use super::cplx::{C64, ONE, ZERO};
use crate::meshing::pml_box::PmlMesh;
use crate::mesh::cells::{CellMesh, CellType};

/// PML profile parameters.
#[derive(Clone, Copy, Debug)]
pub struct Profile {
    /// Polynomial order `m` of the conductivity profile.
    pub order: f64,
    /// Design reflection `R0` of the continuous layer (e.g. `1e-8`).
    pub r0: f64,
    /// Switch the stretching off (the layers then act as plain material): for sanity checks.
    pub enabled: bool,
}

impl Default for Profile {
    fn default() -> Self {
        Profile { order: 2.0, r0: 1e-8, enabled: true }
    }
}

/// Result of a solve.
pub struct Solution {
    pub u: Vec<C64>,
    pub n_dof: usize,
    pub bandwidth: usize,
    pub residual: f64,
}

struct Qp<const D: usize> {
    w: f64,
    x: [f64; D],
    n: [f64; 8],
    dn: [[f64; D]; 8],
}

/// Reference shape functions and gradients at `xi` (unused trailing coordinates are ignored).
fn shape(t: CellType, xi: &[f64; 3]) -> ([f64; 8], [[f64; 3]; 8]) {
    let mut n = [0.0; 8];
    let mut d = [[0.0; 3]; 8];
    const QS: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    match t {
        CellType::Tri => {
            n[..3].copy_from_slice(&[1.0 - xi[0] - xi[1], xi[0], xi[1]]);
            d[..3].copy_from_slice(&[[-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        }
        CellType::Quad => {
            for (i, (sx, sy)) in QS.into_iter().enumerate() {
                n[i] = 0.25 * (1.0 + sx * xi[0]) * (1.0 + sy * xi[1]);
                d[i] = [0.25 * sx * (1.0 + sy * xi[1]), 0.25 * sy * (1.0 + sx * xi[0]), 0.0];
            }
        }
        CellType::Tet => {
            n[..4].copy_from_slice(&[1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]]);
            d[..4].copy_from_slice(&[[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        }
        CellType::Hex => {
            for i in 0..8 {
                let (sx, sy) = QS[i % 4];
                let sz = if i < 4 { -1.0 } else { 1.0 };
                let (fx, fy, fz) = (1.0 + sx * xi[0], 1.0 + sy * xi[1], 1.0 + sz * xi[2]);
                n[i] = 0.125 * fx * fy * fz;
                d[i] = [0.125 * sx * fy * fz, 0.125 * sy * fx * fz, 0.125 * sz * fx * fy];
            }
        }
        CellType::Prism => {
            let l = [1.0 - xi[0] - xi[1], xi[0], xi[1]];
            let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            for i in 0..6 {
                let sz = if i < 3 { -1.0 } else { 1.0 };
                let fz = 0.5 * (1.0 + sz * xi[2]);
                n[i] = l[i % 3] * fz;
                d[i] = [dl[i % 3][0] * fz, dl[i % 3][1] * fz, l[i % 3] * 0.5 * sz];
            }
        }
    }
    (n, d)
}

/// Reference quadrature rule: points and weights.
fn rule(t: CellType) -> Vec<([f64; 3], f64)> {
    let g = 1.0 / 3f64.sqrt();
    match t {
        CellType::Tri => [[0.5, 0.0], [0.5, 0.5], [0.0, 0.5]].iter().map(|p| ([p[0], p[1], 0.0], 1.0 / 6.0)).collect(),
        CellType::Quad => [(-g, -g), (g, -g), (g, g), (-g, g)].iter().map(|p| ([p.0, p.1, 0.0], 1.0)).collect(),
        CellType::Tet => {
            let (a, b) = ((5.0 - 5f64.sqrt()) / 20.0, (5.0 + 3.0 * 5f64.sqrt()) / 20.0);
            [[a, a, a], [b, a, a], [a, b, a], [a, a, b]].iter().map(|p| (*p, 1.0 / 24.0)).collect()
        }
        CellType::Hex => (0..8)
            .map(|i| {
                (
                    [if i & 1 == 0 { -g } else { g }, if i & 2 == 0 { -g } else { g }, if i & 4 == 0 { -g } else { g }],
                    1.0,
                )
            })
            .collect(),
        CellType::Prism => [[0.5, 0.0], [0.5, 0.5], [0.0, 0.5]]
            .iter()
            .flat_map(|p| [-g, g].into_iter().map(move |z| ([p[0], p[1], z], 1.0 / 6.0)))
            .collect(),
    }
}

fn invert3(j: &[[f64; 3]; 3]) -> Option<([[f64; 3]; 3], f64)> {
    let c = [
        [
            j[1][1] * j[2][2] - j[1][2] * j[2][1],
            j[0][2] * j[2][1] - j[0][1] * j[2][2],
            j[0][1] * j[1][2] - j[0][2] * j[1][1],
        ],
        [
            j[1][2] * j[2][0] - j[1][0] * j[2][2],
            j[0][0] * j[2][2] - j[0][2] * j[2][0],
            j[0][2] * j[1][0] - j[0][0] * j[1][2],
        ],
        [
            j[1][0] * j[2][1] - j[1][1] * j[2][0],
            j[0][1] * j[2][0] - j[0][0] * j[2][1],
            j[0][0] * j[1][1] - j[0][1] * j[1][0],
        ],
    ];
    let det = j[0][0] * c[0][0] + j[0][1] * c[1][0] + j[0][2] * c[2][0];
    if !(det > 0.0) {
        return None;
    }
    Some((std::array::from_fn(|r| std::array::from_fn(|q| c[r][q] / det)), det))
}

fn quadrature<const D: usize>(m: &CellMesh<D>, c: usize) -> Result<Vec<Qp<D>>, String> {
    let t = m.types[c];
    if t.dim() != D {
        return Err(format!("cell {c}: {t:?} in a {D}D mesh"));
    }
    let vs = m.cell(c);
    let mut out = Vec::new();
    for (xi, wref) in rule(t) {
        let (n, dref) = shape(t, &xi);
        let mut x = [0.0; D];
        let mut jac = [[0.0; 3]; 3];
        for a in 0..3 {
            jac[a][a] = 1.0;
        }
        for a in 0..D {
            for b in 0..D {
                jac[a][b] = 0.0;
            }
        }
        for (i, &v) in vs.iter().enumerate() {
            let p = m.points[v as usize];
            for a in 0..D {
                x[a] += n[i] * p[a];
                for b in 0..D {
                    jac[a][b] += p[a] * dref[i][b];
                }
            }
        }
        let (inv, det) = invert3(&jac).ok_or_else(|| format!("cell {c} ({t:?}): non-positive Jacobian"))?;
        let mut dn = [[0.0; D]; 8];
        for i in 0..vs.len() {
            for a in 0..D {
                dn[i][a] = (0..D).map(|b| inv[b][a] * dref[i][b]).sum();
            }
        }
        out.push(Qp { w: wref * det, x, n, dn });
    }
    Ok(out)
}

/// Assembled (banded, RCM-ordered) stiffness and mass matrices with Dirichlet rows applied.
pub struct System {
    /// `pos[vertex]` = row of the vertex in the banded matrices.
    pub pos: Vec<usize>,
    pub bandwidth: usize,
    /// Stretched stiffness `sum S/s_a^2 grad N_i . grad N_j`; fixed rows are identity.
    pub stiff: BandMatrix,
    /// Stretched mass `sum S N_i N_j`; fixed rows and columns are zero.
    pub mass: BandMatrix,
    /// Load vector `sum S f N_i`.
    pub rhs: Vec<C64>,
    /// Vertices on the outer boundary of the mesh (homogeneous Dirichlet).
    pub fixed: Vec<bool>,
}

/// Assemble the system for wavenumber `k` (which sets the PML stretching) and source `f`.
pub fn assemble<const D: usize>(
    mesh: &CellMesh<D>,
    pml: Option<(&PmlMesh<D>, Profile)>,
    k: f64,
    f: &dyn Fn(&[f64; D]) -> f64,
) -> Result<System, String> {
    let n = mesh.n_points();
    // stretching factors at a point of a cell with the given label
    let stretch = |label: i32, x: &[f64; D]| -> [C64; D] {
        let mut s = [ONE; D];
        if let Some((pm, prof)) = pml {
            if !prof.enabled {
                return s;
            }
            if let Some((code, _)) = pm.spec.decode(label) {
                for a in 0..D {
                    if code[a] != 0 {
                        let side = if code[a] < 0 { 0 } else { 1 };
                        let t = pm.spec.thickness[a][side];
                        let d = if side == 1 { x[a] - pm.interface[a][1] } else { pm.interface[a][0] - x[a] };
                        let d = d.max(0.0);
                        let sig_max = -(prof.order + 1.0) * prof.r0.ln() / (2.0 * t);
                        let sigma = sig_max * (d / t).powf(prof.order);
                        s[a] = C64::new(1.0, sigma / k);
                    }
                }
            }
        }
        s
    };

    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for c in 0..mesh.n_cells() {
        let vs = mesh.cell(c);
        for &a in vs {
            for &b in vs {
                if a != b {
                    adj[a as usize].push(b);
                }
            }
        }
    }
    for a in adj.iter_mut() {
        a.sort_unstable();
        a.dedup();
    }
    let perm = rcm(&adj);
    let mut pos = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        pos[old as usize] = new;
    }
    let bw = (0..n)
        .flat_map(|u| adj[u].iter().map(move |&v| (u, v as usize)))
        .map(|(u, v)| pos[u].abs_diff(pos[v]))
        .max()
        .unwrap_or(0);
    let mut stiff = BandMatrix::new(n, bw);
    let mut mass = BandMatrix::new(n, bw);
    let mut rhs = vec![ZERO; n];

    let bb = crate::adapt::work::bbox_of(&mesh.points);
    let size = (0..D).map(|a| bb[a][1] - bb[a][0]).fold(0.0, f64::max);
    let on_outer = |p: &[f64; D]| {
        (0..D).any(|ax| (p[ax] - bb[ax][0]).abs() < 1e-9 * size || (p[ax] - bb[ax][1]).abs() < 1e-9 * size)
    };
    let fixed: Vec<bool> = mesh.points.iter().map(on_outer).collect();
    for c in 0..mesh.n_cells() {
        let vs = mesh.cell(c);
        let nn = vs.len();
        let label = mesh.label[c];
        for q in quadrature(mesh, c)? {
            let s = stretch(label, &q.x);
            let big = s.iter().fold(ONE, |p, &z| p * z);
            let wgt = big * q.w;
            let cx: [C64; D] = std::array::from_fn(|a| wgt / (s[a] * s[a]));
            let fv = f(&q.x);
            for i in 0..nn {
                let gi = pos[vs[i] as usize];
                rhs[gi] += wgt * (fv * q.n[i]);
                for j in 0..nn {
                    let gj = pos[vs[j] as usize];
                    let mut kij = ZERO;
                    for a in 0..D {
                        kij += cx[a] * (q.dn[i][a] * q.dn[j][a]);
                    }
                    stiff.add(gi, gj, kij);
                    mass.add(gi, gj, wgt * (q.n[i] * q.n[j]));
                }
            }
        }
    }
    // Dirichlet: identity rows/columns in the stiffness, zero in the mass
    for v in 0..n {
        if fixed[v] {
            let g = pos[v];
            let (lo, hi) = (g.saturating_sub(bw), (g + bw).min(n - 1));
            for j in lo..=hi {
                for (r, c) in [(g, j), (j, g)] {
                    let cur = stiff.get(r, c);
                    stiff.add(r, c, -cur);
                    let cur = mass.get(r, c);
                    mass.add(r, c, -cur);
                }
            }
            stiff.add(g, g, ONE);
            rhs[g] = ZERO;
        }
    }
    Ok(System { pos, bandwidth: bw, stiff, mass, rhs, fixed })
}

/// Solve `-Laplace u - k^2 u = f` with the PML described by `pml` (pass `None` for a plain
/// mesh without stretching). `f` is the real source density.
pub fn solve<const D: usize>(
    mesh: &CellMesh<D>,
    pml: Option<(&PmlMesh<D>, Profile)>,
    k: f64,
    f: &dyn Fn(&[f64; D]) -> f64,
) -> Result<Solution, String> {
    let n = mesh.n_points();
    let sys = assemble(mesh, pml, k, f)?;
    let mut a = sys.stiff;
    a.axpy(C64::real(-k * k), &sys.mass);
    let (x, residual) = a.solve(&sys.rhs)?;
    let u = (0..n).map(|v| x[sys.pos[v]]).collect();
    Ok(Solution { u, n_dof: n, bandwidth: sys.bandwidth, residual })
}

/// Eigenvalue of the (unstretched, real) Dirichlet problem `-Laplace u = lambda u` nearest to
/// `shift`, by inverse iteration with a Rayleigh quotient. Returns `(lambda, iterations)`.
pub fn eigenvalue_near<const D: usize>(
    mesh: &CellMesh<D>,
    shift: f64,
    max_iter: usize,
) -> Result<(f64, usize), String> {
    let sys = assemble(mesh, None, 1.0, &|_: &[f64; D]| 0.0)?;
    let n = mesh.n_points();
    let mut a = BandMatrix::new(n, sys.bandwidth);
    a.axpy(ONE, &sys.stiff);
    a.axpy(C64::real(-shift), &sys.mass);
    let lu: BandLu = a.factor()?;
    // deterministic start vector, zero on the Dirichlet boundary
    let mut x: Vec<C64> = (0..n)
        .map(|v| {
            let row = sys.pos[v];
            let h = (row as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 11;
            C64::real(if sys.fixed[v] { 0.0 } else { (h as f64) / ((1u64 << 53) as f64) - 0.5 })
        })
        .collect();
    let dot = |u: &[C64], v: &[C64]| -> f64 { u.iter().zip(v).map(|(p, q)| p.re * q.re + p.im * q.im).sum() };
    let mut lambda = f64::NAN;
    for it in 1..=max_iter {
        let mx = sys.mass.mul_vec(&x);
        let y = lu.solve(&mx);
        let my = sys.mass.mul_vec(&y);
        let nrm = dot(&y, &my).sqrt();
        x = y.iter().map(|z| *z * (1.0 / nrm)).collect();
        let kx = sys.stiff.mul_vec(&x);
        let mx = sys.mass.mul_vec(&x);
        let new = dot(&x, &kx) / dot(&x, &mx);
        if (new - lambda).abs() < 1e-12 * new.abs() {
            return Ok((new, it));
        }
        lambda = new;
    }
    Ok((lambda, max_iter))
}

/// Relative L2 difference `||u - v|| / ||v||` over the vertex ids `0..n`.
pub fn rel_l2(u: &[C64], v: &[C64], n: usize) -> f64 {
    let num: f64 = (0..n).map(|i| (u[i] - v[i]).norm_sqr()).sum::<f64>().sqrt();
    let den: f64 = (0..n).map(|i| v[i].norm_sqr()).sum::<f64>().sqrt();
    num / den
}

/// Gaussian source of width `w` centred at `x0`.
pub fn gaussian_source<const D: usize>(x0: [f64; D], w: f64) -> impl Fn(&[f64; D]) -> f64 {
    move |x| {
        let r2: f64 = (0..D).map(|a| (x[a] - x0[a]).powi(2)).sum();
        (-r2 / (2.0 * w * w)).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::pml_box::{PmlSpec, build_pml};
    use crate::mesh::fixtures::box_tri;
    use std::f64::consts::PI;

    fn quad_grid(n: usize) -> CellMesh<2> {
        let mut m = CellMesh::<2>::new();
        for j in 0..=n {
            for i in 0..=n {
                m.points.push([i as f64 / n as f64, j as f64 / n as f64]);
            }
        }
        let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
        for j in 0..n {
            for i in 0..n {
                m.push_cell(CellType::Quad, &[id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)], 1);
            }
        }
        m
    }

    fn manufactured_error(mesh: &CellMesh<2>, k: f64) -> f64 {
        let f = move |x: &[f64; 2]| (2.0 * PI * PI - k * k) * (PI * x[0]).sin() * (PI * x[1]).sin();
        let s = solve(mesh, None, k, &f).unwrap();
        assert!(s.residual < 1e-10);
        let mut e: f64 = 0.0;
        for (i, p) in mesh.points.iter().enumerate() {
            let ex = (PI * p[0]).sin() * (PI * p[1]).sin();
            e = e.max((s.u[i] - C64::real(ex)).abs());
        }
        e
    }

    #[test]
    fn p1_and_q1_converge_quadratically() {
        let k = 3.0;
        let p1 = |n| CellMesh::from_simplex(&box_tri([[0.0, 1.0], [0.0, 1.0]], n, n));
        let (a, b) = (manufactured_error(&p1(8), k), manufactured_error(&p1(16), k));
        assert!(a / b > 3.4 && a / b < 4.6, "P1 ratio {}", a / b);
        let (a, b) = (manufactured_error(&quad_grid(8), k), manufactured_error(&quad_grid(16), k));
        assert!(a / b > 3.4 && a / b < 4.6, "Q1 ratio {}", a / b);
    }

    /// Physical-domain error against a thick-PML reference on the same interface mesh.
    fn pml_error(n_layers: usize, thickness: f64, enabled: bool, split: bool) -> f64 {
        let k = 2.0 * PI;
        let phys = box_tri([[0.0, 1.5], [0.0, 1.5]], 15, 15);
        let src = gaussian_source([0.75, 0.75], 0.12);
        let prof = Profile { enabled, ..Default::default() };
        let run = |spec: &PmlSpec<2>, split: bool, prof: Profile| {
            let pm = build_pml(&phys, spec).unwrap();
            let mesh =
                if split { CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap()) } else { pm.mesh.clone() };
            let s = solve(&mesh, Some((&pm, prof)), k, &src).unwrap();
            assert!(s.residual < 1e-8, "residual {}", s.residual);
            s.u
        };
        let reference = run(&PmlSpec::<2>::uniform(1.5, 24), false, Profile::default());
        let u = run(&PmlSpec::<2>::uniform(thickness, n_layers), split, prof);
        rel_l2(&u, &reference, phys.n_points())
    }

    #[test]
    fn pml_absorbs_and_stretching_matters() {
        let with = pml_error(6, 0.75, true, false);
        let without = pml_error(6, 0.75, false, false);
        assert!(with < 0.05, "PML error {with}");
        assert!(without > 5.0 * with, "no-PML error {without} vs {with}");
    }

    #[test]
    fn more_layers_help_and_split_agrees_to_discretisation_error() {
        let coarse = pml_error(2, 0.75, true, false);
        let fine = pml_error(8, 0.75, true, false);
        assert!(fine < coarse, "{fine} vs {coarse}");
        let native = pml_error(6, 0.75, true, false);
        let split = pml_error(6, 0.75, true, true);
        assert!(split < 0.1 && native < 0.1, "{native} {split}");
    }

    // ---- 3D ----

    use crate::mesh::fixtures::box_tet;

    /// Unit-cube hexahedral grid, optionally extruded from triangles into prisms.
    fn hex_grid(n: usize) -> CellMesh<3> {
        let mut m = CellMesh::<3>::new();
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    m.points.push([i as f64 / n as f64, j as f64 / n as f64, k as f64 / n as f64]);
                }
            }
        }
        let id = |i: usize, j: usize, k: usize| (k * (n + 1) * (n + 1) + j * (n + 1) + i) as u32;
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    m.push_cell(
                        CellType::Hex,
                        &[
                            id(i, j, k),
                            id(i + 1, j, k),
                            id(i + 1, j + 1, k),
                            id(i, j + 1, k),
                            id(i, j, k + 1),
                            id(i + 1, j, k + 1),
                            id(i + 1, j + 1, k + 1),
                            id(i, j + 1, k + 1),
                        ],
                        1,
                    );
                }
            }
        }
        m
    }

    fn prism_grid(n: usize) -> CellMesh<3> {
        let tri = box_tri([[0.0, 1.0], [0.0, 1.0]], n, n);
        let np = tri.n_points() as u32;
        let mut m = CellMesh::<3>::new();
        for k in 0..=n {
            for p in &tri.points {
                m.points.push([p[0], p[1], k as f64 / n as f64]);
            }
        }
        for k in 0..n as u32 {
            for c in 0..tri.n_cells() {
                let v = tri.cell(c);
                let (a, b) = (k * np, (k + 1) * np);
                m.push_cell(CellType::Prism, &[a + v[0], a + v[1], a + v[2], b + v[0], b + v[1], b + v[2]], 1);
            }
        }
        m
    }

    fn manufactured_error_3d(mesh: &CellMesh<3>, k: f64) -> f64 {
        mesh.check().unwrap();
        let f = move |x: &[f64; 3]| (3.0 * PI * PI - k * k) * (PI * x[0]).sin() * (PI * x[1]).sin() * (PI * x[2]).sin();
        let s = solve(mesh, None, k, &f).unwrap();
        assert!(s.residual < 1e-9, "residual {}", s.residual);
        let mut e: f64 = 0.0;
        for (i, p) in mesh.points.iter().enumerate() {
            let ex = (PI * p[0]).sin() * (PI * p[1]).sin() * (PI * p[2]).sin();
            e = e.max((s.u[i] - C64::real(ex)).abs());
        }
        e
    }

    #[test]
    fn p1_hex_and_prism_converge_quadratically_in_3d() {
        let k = 2.0;
        let tet = |n| CellMesh::from_simplex(&box_tet([[0.0, 1.0]; 3], [n; 3]));
        for (name, mk) in [
            ("tet", &tet as &dyn Fn(usize) -> CellMesh<3>),
            ("hex", &hex_grid as &dyn Fn(usize) -> CellMesh<3>),
            ("prism", &prism_grid as &dyn Fn(usize) -> CellMesh<3>),
        ] {
            let (a, b) = (manufactured_error_3d(&mk(6), k), manufactured_error_3d(&mk(12), k));
            assert!(a / b > 3.2 && a / b < 5.0, "{name}: ratio {}", a / b);
        }
    }

    /// Dirichlet eigenvalue of the PML-built mesh with the stretching switched off (the layers are
    /// then plain material of a box of side `1 + 2 t`): the lowest eigenvalue is `3 pi^2 / L^2`.
    fn pml_mesh_eigen(n: usize, layers: usize, split: bool) -> (f64, f64) {
        let t = 0.5;
        let phys = box_tet([[0.0, 1.0]; 3], [n; 3]);
        let pm = build_pml(&phys, &PmlSpec::<3>::uniform(t, layers)).unwrap();
        pm.mesh.check().unwrap();
        let mesh = if split { CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap()) } else { pm.mesh.clone() };
        let l = 1.0 + 2.0 * t;
        let exact = 3.0 * PI * PI / (l * l);
        let (lam, _) = eigenvalue_near(&mesh, 0.8 * exact, 30).unwrap();
        (lam, exact)
    }

    #[test]
    fn mixed_pml_mesh_native_and_split_eigenvalues_agree_and_converge() {
        let mut prev: Option<(f64, f64)> = None;
        let mut last = (0.0, 0.0);
        for (n, layers) in [(4usize, 2usize), (6, 3), (8, 4)] {
            let (nat, exact) = pml_mesh_eigen(n, layers, false);
            let (spl, _) = pml_mesh_eigen(n, layers, true);
            let (en, es) = ((nat - exact).abs() / exact, (spl - exact).abs() / exact);
            println!("n={n}: native {nat:.6} split {spl:.6} exact {exact:.6}  rel err {en:.2e} {es:.2e}");
            assert!(en < 0.1 && es < 0.1);
            if let Some((pn, ps)) = prev {
                assert!(en < pn && es < ps, "no convergence");
            }
            prev = Some((en, es));
            last = (nat, spl);
        }
        // the finest native and split eigenvalues agree to the discretisation error of the coarser one
        assert!((last.0 - last.1).abs() / last.0 < 0.02, "{} vs {}", last.0, last.1);
    }

    fn pml_error_3d(n_layers: usize, split: bool, enabled: bool) -> f64 {
        let k = 2.0 * PI;
        let phys = box_tet([[0.0, 1.0]; 3], [6; 3]);
        let src = gaussian_source([0.5, 0.5, 0.5], 0.12);
        let run = |spec: &PmlSpec<3>, split: bool, prof: Profile| {
            let pm = build_pml(&phys, spec).unwrap();
            let mesh =
                if split { CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap()) } else { pm.mesh.clone() };
            let s = solve(&mesh, Some((&pm, prof)), k, &src).unwrap();
            assert!(s.residual < 1e-8, "residual {}", s.residual);
            s.u
        };
        let reference = run(&PmlSpec::<3>::uniform(0.5, 5), split, Profile::default());
        let u = run(&PmlSpec::<3>::uniform(0.5, n_layers), split, Profile { enabled, ..Default::default() });
        rel_l2(&u, &reference, phys.n_points())
    }

    #[test]
    fn pml_works_in_3d_native_and_split() {
        for split in [false, true] {
            let with = pml_error_3d(3, split, true);
            let without = pml_error_3d(3, split, false);
            assert!(with < 0.15, "split={split}: PML error {with}");
            assert!(without > 3.0 * with, "split={split}: no-PML error {without} vs {with}");
        }
    }
}
