//! Element sizes from the physics: how fine a mesh must be for the waves the study looks for, before any solution exists. The core takes the
//! neutral `meshing::Sizing` a physics reports (`cell_sizes_for`, `vertex_sizes_for`, `metrics_for`, `adapt_tets_for`); the
//! physics crates state their `Sizing` themselves (`Physics::sizing`); `optical_wavelength` and `elastic_wavelength` are the closed forms they use.
//!
//! Every wave problem has a shortest wavelength per domain, and a mesh resolves it with a number of points per wavelength that depends on the
//! element order. This module finds that wavelength from the materials and the study, divides it by the points per wavelength, and hands back
//! a size per cell, per vertex, or an isotropic metric per vertex for `meshing::adapt_to_metrics`.
//!
//! * `Emw` mode analysis at `freq`: `lambda = c0 / (f n)` with `n^2 = mu_r * (largest eigenvalue of the permittivity tensor)`.
//! * `Solid` eigenfrequency study with `f_search > 0`: `lambda = c_min / f_search`, `c_min` the slowest plane wave of the (possibly
//!   anisotropic, rotated) stiffness, found from the Christoffel matrix over a sampled sphere of directions.
//!
//! Domains listed in `inactive` do not carry elastic waves and are left out of the elastic size. Domains with a perfectly matched layer
//! may be coarser by `pml_relax`: the wave is damped there and the layer is not where the answer is read.
//!
//! The points per wavelength are starting values (12 for order 1, 6 for order 2), meant to give a mode index within about 1e-3; the tests
//! of this module measure them on a closed form. They are not a guarantee for geometry that varies faster than the wave (`Size` rules
//! and curvature stay with `meshing::geometry_floor`; combine the two by intersecting their metrics).

use crate::material::M6;
use crate::mesh::SimplexMesh;
use super::sizing::Sizing;
use crate::adapt::metric::{isotropic, sym_eig, Sym};

/// Vacuum speed of light (m/s).
const C0: f64 = 299_792_458.0;

/// How finely the waves are resolved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicsSizing {
    /// Points per shortest wavelength; `None` takes 12 for order 1 and 6 for order 2.
    pub points_per_wavelength: Option<f64>,
    /// Factor on the size inside perfectly matched layers (1 keeps it).
    pub pml_relax: f64,
    /// Eigenfrequency studies look for `n_modes` modes at and above the search frequency; the mesh must resolve the highest of them. The
    /// search frequency is multiplied by this to stand for it.
    pub eigen_headroom: f64,
}

impl Default for PhysicsSizing {
    fn default() -> Self {
        PhysicsSizing { points_per_wavelength: None, pml_relax: 1.5, eigen_headroom: 1.5 }
    }
}

impl PhysicsSizing {
    /// Points per wavelength for an element order.
    pub fn ppw(&self, order: usize) -> f64 {
        self.points_per_wavelength.unwrap_or(if order >= 2 { 6.0 } else { 12.0 })
    }
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

// ---- from a neutral `Sizing` (what a physics reports) to element sizes ----

/// Target element size of every cell: the shortest wavelength of its domain, divided by the headroom and by the points per wavelength, and
/// relaxed in PML domains. `resolve` turns the region names of `sizing` into domain labels (`sizing::label_names` for a mesh with numbers
/// only). Cells of domains without a wave get `f64::INFINITY`.
pub fn cell_sizes_for(
    sizing: &Sizing,
    resolve: &dyn Fn(&str) -> Option<i32>,
    cell_labels: &[i32],
    pml: &dyn Fn(i32) -> bool,
    sz: &PhysicsSizing,
) -> Result<Vec<f64>, String> {
    let order = sizing.order as usize;
    if !(sz.pml_relax >= 1.0 && sz.ppw(order) > 0.0 && sizing.headroom > 0.0) {
        return Err("points per wavelength and headroom must be positive and pml_relax at least 1".into());
    }
    let lam = sizing.resolve(resolve)?;
    let ppw = sz.ppw(order);
    Ok(cell_labels
        .iter()
        .map(|&d| match lam.get(&d) {
            Some(&l) => l / sizing.headroom / ppw * if pml(d) { sz.pml_relax } else { 1.0 },
            None => f64::INFINITY,
        })
        .collect())
}

/// Target element size at every vertex: the smallest of the cells that meet there. A vertex of cells without a wave gets `f64::INFINITY`.
pub fn vertex_sizes_for(
    sizing: &Sizing,
    resolve: &dyn Fn(&str) -> Option<i32>,
    mesh: &SimplexMesh,
    pml: &dyn Fn(i32) -> bool,
    sz: &PhysicsSizing,
) -> Result<Vec<f64>, String> {
    let cs = cell_sizes_for(sizing, resolve, &mesh.cell_labels, pml, sz)?;
    let mut vh = vec![f64::INFINITY; mesh.n_points()];
    for (c, &h) in cs.iter().enumerate() {
        for &v in mesh.cell(c) {
            vh[v as usize] = vh[v as usize].min(h);
        }
    }
    Ok(vh)
}

/// Isotropic metric of the physics size at every vertex of a mesh of any dimension (`D` = 2 or 3). Vertices without a wave get `ceiling` (a
/// size the caller would accept anywhere, typically the global size).
pub fn metrics_for<const D: usize>(
    sizing: &Sizing,
    resolve: &dyn Fn(&str) -> Option<i32>,
    mesh: &SimplexMesh,
    pml: &dyn Fn(i32) -> bool,
    sz: &PhysicsSizing,
    ceiling: f64,
) -> Result<Vec<Sym<D>>, String> {
    if mesh.dim != D {
        return Err(format!("a {}D mesh, but metrics for {D}D were asked", mesh.dim));
    }
    Ok(vertex_sizes_for(sizing, resolve, mesh, pml, sz)?.into_iter().map(|h| isotropic::<D>(h.min(ceiling).max(1e-300))).collect())
}

/// Adapts a tetrahedral mesh to a neutral `Sizing`: the sizes of `vertex_sizes_for` (capped at `ceiling`), graded with slope `grade`,
/// carried to the new vertices by a background field. Domain labels and boundary-facet labels pass through; corners, creases and
/// material interfaces are kept as `adapt::topo3` classifies them. Returns the mesh and the report.
#[allow(clippy::too_many_arguments)]
pub fn adapt_tets_for(
    mesh: &SimplexMesh,
    sizing: &Sizing,
    resolve: &dyn Fn(&str) -> Option<i32>,
    pml: &dyn Fn(i32) -> bool,
    sz: &PhysicsSizing,
    ceiling: f64,
    grade: f64,
    opt: &crate::adapt::ops3::AdaptOptions3,
) -> Result<(SimplexMesh, crate::adapt::drive::AdaptReport), String> {
    use crate::adapt::work::{simplex_from_adapt, simplex_to_adapt};
    use crate::adapt::{drive3::adapt3, field3::BackgroundMetric3, topo3::TopoOptions3};
    let metrics = metrics_for::<3>(sizing, resolve, mesh, pml, sz, ceiling)?;
    let old = simplex_to_adapt::<3>(mesh)?;
    let field = BackgroundMetric3::new_graded(&old, &metrics, grade);
    let out = adapt3(&old, &field, &TopoOptions3::default(), opt)?;
    Ok((simplex_from_adapt(&out.mesh), out.report))
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
        let eps = [[4.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((optical_wavelength(&eps, 1.0, 1.0) - 0.5).abs() < 1e-12, "the largest eigenvalue of the permittivity sets the index");
        assert!((optical_wavelength(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 4.0, speed_of_light() / 10e9) - speed_of_light() / 10e9 / 2.0).abs() < 1e-12);
    }

    #[test]
    fn the_slowest_wave_sets_the_elastic_wavelength() {
        let c = isotropic_stiffness(169e9, 0.28);
        let (vs, _) = wave_speed_range(&c, 2330.0);
        assert!((elastic_wavelength(&c, 2330.0, 1e9) - vs / 1e9).abs() < 1e-18);
    }
}

// ---- wavelengths of a medium, for the physics crates that report their own `Sizing` ----

/// Free-space wavelength `lam0` shortened by the largest refractive index of a medium with the relative permittivity tensor `eps` and the
/// relative permeability `mu`: `n^2 = mu * (largest eigenvalue of eps)`.
pub fn optical_wavelength(eps: &[[f64; 3]; 3], mu: f64, lam0: f64) -> f64 {
    let (ev, _) = sym_eig::<3>(eps);
    lam0 / (ev[0] * mu).max(1e-12).sqrt()
}

/// Vacuum speed of light (m/s).
pub fn speed_of_light() -> f64 {
    C0
}

/// Shortest wavelength of the elastic waves of a medium with stiffness `c` (solver Voigt order) and density `rho` at frequency `f`: the
/// slowest plane wave over the frequency.
pub fn elastic_wavelength(c: &M6, rho: f64, f: f64) -> f64 {
    wave_speed_range(c, rho).0 / f
}
