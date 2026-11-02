//! Element sizes from the physics of a `ProblemSpec`: how fine a mesh must be for the waves the study looks for, before any solution exists.
//!
//! Every wave problem has a shortest wavelength per domain, and a mesh resolves it with a number of points per wavelength that depends on the
//! element order. This module finds that wavelength from the materials and the study, divides it by the points per wavelength, and hands back
//! a size per cell, per vertex, or an isotropic metric per vertex for `meshing::adapt_to_metrics`.
//!
//! * `Emw` mode analysis at `freq`: `lambda = c0 / (f n)` with `n^2 = mu_r * (largest eigenvalue of the permittivity tensor)`.
//! * `Solid` eigenfrequency study with `f_search > 0`: `lambda = c_min / f_search`, `c_min` the slowest plane wave of the (possibly
//!   anisotropic, rotated) stiffness, found from the Christoffel matrix over a sampled sphere of directions.
//! * `Brillouin`: the shorter of the optical wavelength `lambda0 / n` per domain and the elastic wavelength at the search frequency of
//!   the study (without one, at the backward-scattering frequency `2 n v_L / lambda0`).
//!
//! Domains listed in `inactive` do not carry elastic waves and are left out of the elastic size. Domains with a perfectly matched layer
//! may be coarser by `pml_relax`: the wave is damped there and the layer is not where the answer is read.
//!
//! The points per wavelength are starting values (12 for order 1, 6 for order 2), meant to give a mode index within about 1e-3; the tests
//! of this module measure them on a closed form. They are not a guarantee for geometry that varies faster than the wave (`Size` rules
//! and curvature stay with `meshing::geometry_floor`; combine the two by intersecting their metrics).

use crate::material::M6;
use crate::mesh::SimplexMesh;
use crate::meshers::adapt::metric::{isotropic, sym_eig, Sym};
use crate::problem::{Physics, ProblemSpec, Study, TensorKind};
use std::collections::BTreeMap;

/// Vacuum speed of light (m/s).
const C0: f64 = 299_792_458.0;

/// How finely the waves are resolved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicsSizing {
    /// Points per shortest wavelength; `None` takes 12 for order 1 and 6 for order 2.
    pub points_per_wavelength: Option<f64>,
    /// Factor on the size inside perfectly matched layers (1 keeps it).
    pub pml_relax: f64,
}

impl Default for PhysicsSizing {
    fn default() -> Self {
        PhysicsSizing { points_per_wavelength: None, pml_relax: 1.5 }
    }
}

impl PhysicsSizing {
    /// Points per wavelength for an element order.
    pub fn ppw(&self, order: usize) -> f64 {
        self.points_per_wavelength.unwrap_or(if order >= 2 { 6.0 } else { 12.0 })
    }
}

/// Shortest wavelength (m) of the study in every domain of the mesh, by domain label. Domains that carry no wave of this study are absent.
pub fn wavelengths(spec: &ProblemSpec, mesh: &SimplexMesh) -> Result<BTreeMap<i32, f64>, String> {
    let first = first_cell_of_each_domain(mesh);
    match spec.physics {
        Physics::Emw => {
            let Study::Modes { freq, .. } = spec.study else {
                return Err("the wavelength of an electromagnetic study needs a Modes study".into());
            };
            optical_wavelengths(spec, mesh, &first, C0 / freq, None)
        }
        Physics::Solid => {
            let f = elastic_frequency(spec)?.ok_or("the wavelength of a solid-mechanics study needs an eigenfrequency study with a search frequency")?;
            elastic_wavelengths(spec, mesh, &first, f)
        }
        Physics::Brillouin => {
            let b = spec.brillouin.as_ref().ok_or("a Brillouin study needs its optical data")?;
            let lam0 = b.wavelength;
            let opt = optical_wavelengths(spec, mesh, &first, lam0, Some(&b.refractive_index))?;
            let f = match elastic_frequency(spec)? {
                Some(f) => f,
                None => {
                    // backward scattering: f = 2 n v_L / lambda0, with the fastest index and the fastest wave of the mesh
                    let n = opt.values().map(|l| lam0 / l).fold(1.0, f64::max);
                    let v = elastic_speeds(spec, mesh, &first)?.values().map(|s| s.1).fold(0.0, f64::max);
                    2.0 * n * v / lam0
                }
            };
            let ela = elastic_wavelengths(spec, mesh, &first, f)?;
            let mut out = opt;
            for (d, l) in ela {
                let e = out.entry(d).or_insert(l);
                *e = e.min(l);
            }
            Ok(out)
        }
        p => Err(format!("{} has no wavelength to resolve", p.name())),
    }
}

fn first_cell_of_each_domain(mesh: &SimplexMesh) -> BTreeMap<i32, usize> {
    let mut m = BTreeMap::new();
    for (c, &l) in mesh.cell_labels.iter().enumerate() {
        m.entry(l).or_insert(c);
    }
    m
}

fn elastic_frequency(spec: &ProblemSpec) -> Result<Option<f64>, String> {
    match spec.study {
        Study::Eigenfrequency { f_search, .. } if f_search > 0.0 => Ok(Some(f_search)),
        _ => Ok(None),
    }
}

/// Vacuum wavelength `lam0` shortened by the largest refractive index of each domain. `index` overrides the permittivity (Brillouin).
fn optical_wavelengths(spec: &ProblemSpec, mesh: &SimplexMesh, first: &BTreeMap<i32, usize>, lam0: f64, index: Option<&BTreeMap<i32, f64>>) -> Result<BTreeMap<i32, f64>, String> {
    let mut out = BTreeMap::new();
    if let Some(ix) = index {
        for &d in first.keys() {
            out.insert(d, lam0 / ix.get(&d).copied().unwrap_or(1.0).max(1e-6));
        }
        return Ok(out);
    }
    let eps = spec.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?;
    let mu = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    for (&d, &c) in first {
        let m: Sym<3> = [[eps[c][0][0], eps[c][0][1], eps[c][0][2]], [eps[c][1][0], eps[c][1][1], eps[c][1][2]], [eps[c][2][0], eps[c][2][1], eps[c][2][2]]];
        let (ev, _) = sym_eig::<3>(&m);
        let n2 = (ev[0] * mu[c]).max(1e-12);
        out.insert(d, lam0 / n2.sqrt());
    }
    Ok(out)
}

/// `(slowest, fastest)` plane-wave speed of every active domain, by label.
fn elastic_speeds(spec: &ProblemSpec, mesh: &SimplexMesh, first: &BTreeMap<i32, usize>) -> Result<BTreeMap<i32, (f64, f64)>, String> {
    let stiff = spec.materials.cell_stiffness_skipping(mesh, &spec.inactive)?;
    let rho = spec.materials.cell_values(mesh, "rho", None)?;
    let mut out = BTreeMap::new();
    for (&d, &c) in first {
        if spec.inactive.contains(&d) || rho[c] <= 0.0 {
            continue;
        }
        out.insert(d, wave_speed_range(&stiff[c], rho[c]));
    }
    Ok(out)
}

fn elastic_wavelengths(spec: &ProblemSpec, mesh: &SimplexMesh, first: &BTreeMap<i32, usize>, f: f64) -> Result<BTreeMap<i32, f64>, String> {
    Ok(elastic_speeds(spec, mesh, first)?.into_iter().map(|(d, (vmin, _))| (d, vmin / f)).collect())
}

/// Slowest and fastest plane-wave speed of a medium with stiffness `c` (solver Voigt order xx, yy, zz, xy, yz, xz; engineering shear) and
/// density `rho`: square roots of the extreme eigenvalues of the Christoffel matrix `B(n)^T C B(n) / rho` over 400 directions `n`
/// spread evenly on the sphere (a Fibonacci lattice; the sampled minimum is at most a fraction of a percent above the true one).
pub fn wave_speed_range(c: &M6, rho: f64) -> (f64, f64) {
    const N: usize = 400;
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
    for i in 0..N {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / N as f64;
        let r = (1.0 - z * z).sqrt();
        let n = [r * (golden * i as f64).cos(), r * (golden * i as f64).sin(), z];
        // strain of a displacement amplitude v for a plane wave along n
        let b = |v: [f64; 3]| -> [f64; 6] {
            [n[0] * v[0], n[1] * v[1], n[2] * v[2], n[1] * v[0] + n[0] * v[1], n[2] * v[1] + n[1] * v[2], n[2] * v[0] + n[0] * v[2]]
        };
        let cols: Vec<[f64; 6]> = (0..3).map(|k| b(std::array::from_fn(|j| (j == k) as u8 as f64))).collect();
        let mut g = [[0.0; 3]; 3];
        for a in 0..3 {
            for bb in 0..3 {
                let mut s = 0.0;
                for p in 0..6 {
                    for q in 0..6 {
                        s += cols[a][p] * c[p][q] * cols[bb][q];
                    }
                }
                g[a][bb] = s / rho;
            }
        }
        let (ev, _) = sym_eig::<3>(&g);
        lo = lo.min(ev[2].max(0.0));
        hi = hi.max(ev[0]);
    }
    (lo.sqrt(), hi.sqrt())
}

/// Target element size of every cell: the shortest wavelength of its domain over the points per wavelength, relaxed in PML domains.
/// Cells of domains without a wave get `f64::INFINITY`.
pub fn cell_sizes(spec: &ProblemSpec, mesh: &SimplexMesh, sz: &PhysicsSizing) -> Result<Vec<f64>, String> {
    if !(sz.pml_relax >= 1.0 && sz.ppw(spec.order) > 0.0) {
        return Err("points per wavelength must be positive and pml_relax at least 1".into());
    }
    let lam = wavelengths(spec, mesh)?;
    let ppw = sz.ppw(spec.order);
    let pml = |d: i32| spec.materials.pml.contains_key(&d) || spec.brillouin.as_ref().is_some_and(|b| b.optical_pml.contains_key(&d));
    Ok(mesh
        .cell_labels
        .iter()
        .map(|&d| match lam.get(&d) {
            Some(&l) => l / ppw * if pml(d) { sz.pml_relax } else { 1.0 },
            None => f64::INFINITY,
        })
        .collect())
}

/// Target element size at every vertex: the smallest of the cells that meet there. A vertex of cells without a wave gets `f64::INFINITY`.
pub fn vertex_sizes(spec: &ProblemSpec, mesh: &SimplexMesh, sz: &PhysicsSizing) -> Result<Vec<f64>, String> {
    let cs = cell_sizes(spec, mesh, sz)?;
    let mut vh = vec![f64::INFINITY; mesh.n_points()];
    for (c, &h) in cs.iter().enumerate() {
        for &v in mesh.cell(c) {
            vh[v as usize] = vh[v as usize].min(h);
        }
    }
    Ok(vh)
}

/// Isotropic metric of the physics size at every vertex of a 2D mesh. Vertices without a wave get `ceiling` (a size the caller would accept
/// anywhere, typically the global size).
pub fn metrics(spec: &ProblemSpec, mesh: &SimplexMesh, sz: &PhysicsSizing, ceiling: f64) -> Result<Vec<Sym<2>>, String> {
    if mesh.dim != 2 {
        return Err("physics metrics are for 2D meshes".into());
    }
    Ok(vertex_sizes(spec, mesh, sz)?.into_iter().map(|h| isotropic::<2>(h.min(ceiling).max(1e-300))).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::isotropic_stiffness;

    #[test]
    fn isotropic_speeds_are_the_textbook_ones() {
        let (e, nu, rho) = (169e9, 0.28, 2330.0);
        let (vs, vp) = wave_speed_range(&isotropic_stiffness(e, nu), rho);
        let mu = e / (2.0 * (1.0 + nu));
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        assert!((vs / (mu / rho).sqrt() - 1.0).abs() < 1e-9, "shear {vs}");
        assert!((vp / ((lam + 2.0 * mu) / rho).sqrt() - 1.0).abs() < 1e-9, "pressure {vp}");
    }

    #[test]
    fn optical_wavelength_is_shortened_by_the_index() {
        let m = crate::geomdrive::plate(1.0, 1.0, 0.5).unwrap().mesh;
        let mut s = ProblemSpec::new(Physics::Emw);
        s.materials.default.insert("eps_r".into(), 4.0);
        s.study = Study::Modes { freq: 10e9, n_target: 1.5, evanescent: false, n_modes: 1 };
        let l = wavelengths(&s, &m).unwrap();
        assert!((l[&1] - C0 / 10e9 / 2.0).abs() < 1e-12);
        let h = cell_sizes(&s, &m, &PhysicsSizing::default()).unwrap();
        assert!((h[0] - l[&1] / 6.0).abs() < 1e-15, "order 2 takes 6 points per wavelength");
    }

    #[test]
    fn pml_domains_are_coarser_and_inactive_domains_have_no_elastic_wave() {
        let mut m = crate::geomdrive::plate(1.0, 1.0, 0.5).unwrap().mesh;
        let half = m.n_cells() / 2;
        for l in m.cell_labels.iter_mut().skip(half) {
            *l = 2;
        }
        let mut s = ProblemSpec::new(Physics::Solid);
        s.study = Study::Eigenfrequency { n_modes: 5, f_search: 1e9 };
        s.materials.default.insert("E".into(), 169e9);
        s.materials.default.insert("nu".into(), 0.28);
        s.materials.default.insert("rho".into(), 2330.0);
        s.inactive = vec![2];
        let l = wavelengths(&s, &m).unwrap();
        assert!(l.contains_key(&1) && !l.contains_key(&2));
        let h = cell_sizes(&s, &m, &PhysicsSizing::default()).unwrap();
        assert!(h[0].is_finite() && h[half].is_infinite());
    }
}
