//! Element sizes from the physics: how fine a mesh must be for the waves the study looks for, before any solution exists. The core takes the
//! neutral `meshing::Sizing` a physics reports (`cell_sizes_for`, `vertex_sizes_for`, `metrics_for`, `adapt_tets_for`); the
//! physics crates state their `Sizing` themselves (`Physics::sizing`); `optical_wavelength` and `elastic_wavelength` are the closed forms they use.
//!
//! Every wave problem has a shortest wavelength per domain, and a mesh resolves it with a number of points per wavelength that depends on the
//! element order. This module finds that wavelength from the materials and the study, divides it by the points per wavelength, and hands back
//! a size per cell, per vertex, or an isotropic metric per vertex for `meshing::adapt_to_metrics`.
//!
//! * `ems` mode analysis at `freq`: `lambda = c0 / (f n)` with `n^2 = mu_r * (largest eigenvalue of the permittivity tensor)`.
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

/// Pre-solve Christoffel stability diagnostic for an anisotropic elastic PML with outward normal `n`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackwardWaveCheck {
    /// Minimum normalized group-velocity projection `(v_g . n_hat) / |v_g|` over outgoing-phase
    /// bulk acoustic branches (`k_hat . n_hat > 0`). Negative whenever a backward wave exists.
    pub min_projection: f64,
    /// True when at least one bulk acoustic branch exhibits antiparallel phase and group normal
    /// components (`(k_hat . n_hat) * (v_g . n_hat) < 0`).
    pub has_backward_wave: bool,
    /// Recommended M-PML auxiliary tangential damping ratio `p = 1.25 * p_min` (`0.0` if stable).
    pub recommended_mpml_ratio: f64,
}

/// Evaluates bulk acoustic group velocities `v_g = d omega / d k` via the Hellmann-Feynman / acoustic
/// Poynting theorem across the 400-point Fibonacci sphere (plus targeted near-grazing rings relative
/// to `outward_normal`) to detect antiparallel backward-wave branches `(k . n_hat)(v_g . n_hat) < 0`
/// and compute the stabilizing M-PML tangential damping ratio (`1.25 * p_min`).
pub fn backward_wave_mpml_ratio(c: &M6, rho: f64, outward_normal: [f64; 3]) -> BackwardWaveCheck {
    let n_norm = (outward_normal[0] * outward_normal[0]
        + outward_normal[1] * outward_normal[1]
        + outward_normal[2] * outward_normal[2])
        .sqrt();
    if !rho.is_finite() || rho <= 0.0 || !n_norm.is_finite() || n_norm <= 1e-30 {
        return BackwardWaveCheck {
            min_projection: 0.0,
            has_backward_wave: false,
            recommended_mpml_ratio: 0.0,
        };
    }
    let n_hat = [outward_normal[0] / n_norm, outward_normal[1] / n_norm, outward_normal[2] / n_norm];

    // Orthonormal tangent basis `(t1, t2)` perpendicular to `n_hat` for high-resolution grazing scan
    let ref_v = if n_hat[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let dot_rn = ref_v[0] * n_hat[0] + ref_v[1] * n_hat[1] + ref_v[2] * n_hat[2];
    let mut t1 = [ref_v[0] - dot_rn * n_hat[0], ref_v[1] - dot_rn * n_hat[1], ref_v[2] - dot_rn * n_hat[2]];
    let t1_len = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    t1 = [t1[0] / t1_len, t1[1] / t1_len, t1[2] / t1_len];
    let t2 = [
        n_hat[1] * t1[2] - n_hat[2] * t1[1],
        n_hat[2] * t1[0] - n_hat[0] * t1[2],
        n_hat[0] * t1[1] - n_hat[1] * t1[0],
    ];

    let mut min_proj = 1.0_f64;
    let mut p_min = 0.0_f64;

    let mut eval_dir = |mut q: [f64; 3]| {
        let mut q_dot_n = q[0] * n_hat[0] + q[1] * n_hat[1] + q[2] * n_hat[2];
        if q_dot_n < 0.0 {
            q = [-q[0], -q[1], -q[2]];
            q_dot_n = -q_dot_n;
        }
        if q_dot_n <= 1e-6 {
            return;
        }
        let b_op = |v: [f64; 3]| -> [f64; 6] {
            [
                q[0] * v[0],
                q[1] * v[1],
                q[2] * v[2],
                q[1] * v[0] + q[0] * v[1],
                q[2] * v[1] + q[1] * v[2],
                q[2] * v[0] + q[0] * v[2],
            ]
        };
        let cols: [[f64; 6]; 3] = [
            b_op([1.0, 0.0, 0.0]),
            b_op([0.0, 1.0, 0.0]),
            b_op([0.0, 0.0, 1.0]),
        ];
        let mut g = [[0.0; 3]; 3];
        for a in 0..3 {
            for bb in 0..3 {
                let mut s = 0.0;
                for p_idx in 0..6 {
                    for q_idx in 0..6 {
                        s += cols[a][p_idx] * c[p_idx][q_idx] * cols[bb][q_idx];
                    }
                }
                g[a][bb] = s / rho;
            }
        }
        let (ev, vecs) = sym_eig::<3>(&g);
        for m in 0..3 {
            if ev[m] <= 1e-30 {
                continue;
            }
            let vp = ev[m].sqrt();
            let p = [vecs[0][m], vecs[1][m], vecs[2][m]];
            let eps = b_op(p);
            let mut sig = [0.0; 6];
            for r_idx in 0..6 {
                let mut acc = 0.0;
                for c_idx in 0..6 {
                    acc += c[r_idx][c_idx] * eps[c_idx];
                }
                sig[r_idx] = acc;
            }
            // Exact Hellmann-Feynman group velocity `v_{g,i} = ((B(e_i) p) . (C B(q) p)) / (rho v_p)`
            // Solver Voigt order: 0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz
            let inv_rvp = 1.0 / (rho * vp);
            let vg = [
                (sig[0] * p[0] + sig[3] * p[1] + sig[5] * p[2]) * inv_rvp,
                (sig[3] * p[0] + sig[1] * p[1] + sig[4] * p[2]) * inv_rvp,
                (sig[5] * p[0] + sig[4] * p[1] + sig[2] * p[2]) * inv_rvp,
            ];
            let vg_norm = (vg[0] * vg[0] + vg[1] * vg[1] + vg[2] * vg[2]).sqrt();
            if vg_norm <= 1e-30 {
                continue;
            }
            let vg_dot_n = vg[0] * n_hat[0] + vg[1] * n_hat[1] + vg[2] * n_hat[2];
            let proj = vg_dot_n / vg_norm;
            min_proj = min_proj.min(proj);

            let kn_vgn = q_dot_n * vg_dot_n;
            if kn_vgn < -1e-9 * vp {
                // Euler homogeneity: `q . v_g == v_p`, so `q_perp . v_{g,perp} == v_p - kn_vgn > v_p > 0`
                let kp_vgp = (vp - kn_vgn).max(1e-30);
                p_min = p_min.max((-kn_vgn) / kp_vgp);
            }
        }
    };

    // 1. Global 400-point Fibonacci sphere (folded into outgoing hemisphere `q . n_hat >= 0`)
    const N_FIB: usize = 400;
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    for i in 0..N_FIB {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / N_FIB as f64;
        let r = (1.0 - z * z).sqrt();
        eval_dir([r * (golden * i as f64).cos(), r * (golden * i as f64).sin(), z]);
    }

    // 2. Near-grazing rings in `(t1, t2, n_hat)` where slowness skewness backward wedges concentrate
    for &sin_el in &[0.015_f64, 0.04, 0.08, 0.14, 0.22, 0.32] {
        let cos_el = (1.0 - sin_el * sin_el).sqrt();
        const N_AZ: usize = 36;
        for j in 0..N_AZ {
            let phi = 2.0 * std::f64::consts::PI * (j as f64) / (N_AZ as f64);
            let (sp, cp) = phi.sin_cos();
            let q = [
                cos_el * (cp * t1[0] + sp * t2[0]) + sin_el * n_hat[0],
                cos_el * (cp * t1[1] + sp * t2[1]) + sin_el * n_hat[1],
                cos_el * (cp * t1[2] + sp * t2[2]) + sin_el * n_hat[2],
            ];
            eval_dir(q);
        }
    }

    let has_backward_wave = p_min > 1e-6 && min_proj < -1e-6;
    let recommended_mpml_ratio = if has_backward_wave {
        (1.25 * p_min).clamp(0.01, 1.0)
    } else {
        0.0
    };
    BackwardWaveCheck {
        min_projection: min_proj,
        has_backward_wave,
        recommended_mpml_ratio,
    }
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

    #[test]
    fn christoffel_backward_wave_check_distinguishes_isotropic_from_skewed_crystals() {
        use crate::material::{sapphire, silicon, stiffness_global, Stiffness};

        // 1. Isotropic media have v_g || k everywhere: never backward-propagating for any normal
        let iso = isotropic_stiffness(72e9, 0.17);
        for n in [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]] {
            let chk = backward_wave_mpml_ratio(&iso, 2200.0, n);
            assert!(!chk.has_backward_wave, "isotropic medium flagged backward wave: {chk:?}");
            assert!(chk.min_projection > 0.0, "isotropic min_projection = {}", chk.min_projection);
            assert_eq!(chk.recommended_mpml_ratio, 0.0);
        }

        // 2. Tilted R-plane sapphire (Al2O3, trigonal c14 = -23.5 GPa, R-cut ~57.6 deg)
        let (saph_stiff, _, _, _, saph_rho) = sapphire();
        let c_saph_rcut = stiffness_global(&saph_stiff, [0.0, 57.6, 0.0]);
        let saph_y = backward_wave_mpml_ratio(&c_saph_rcut, saph_rho, [0.0, 1.0, 0.0]);
        let saph_z = backward_wave_mpml_ratio(&c_saph_rcut, saph_rho, [0.0, 0.0, 1.0]);
        assert!(
            saph_y.has_backward_wave || saph_z.has_backward_wave,
            "tilted R-plane sapphire must exhibit slowness skewness backward waves: y={saph_y:?}, z={saph_z:?}"
        );
        let best_saph = if saph_y.recommended_mpml_ratio >= saph_z.recommended_mpml_ratio { saph_y } else { saph_z };
        assert!(best_saph.min_projection < 0.0);
        assert!(best_saph.recommended_mpml_ratio > 0.01 && best_saph.recommended_mpml_ratio < 0.5);

        // 3. Rotated cubic Silicon (euler = [22.5, 0.0, 0.0], Zener ratio A = 1.564) and rotated LiNbO3
        let c_si_rot = stiffness_global(&silicon(), [22.5, 0.0, 0.0]);
        let si_chk = backward_wave_mpml_ratio(&c_si_rot, 2329.0, [1.0, 0.0, 0.0]);
        assert!(si_chk.has_backward_wave, "rotated Si [22.5 deg] must exhibit qSV backward wave: {si_chk:?}");
        assert!(si_chk.min_projection < 0.0 && si_chk.recommended_mpml_ratio > 0.02);

        let linbo3 = Stiffness::Trigonal {
            c11: 203.0e9,
            c12: 57.3e9,
            c13: 75.2e9,
            c14: 8.5e9,
            c33: 242.4e9,
            c44: 59.5e9,
        };
        let c_ln_ycut = stiffness_global(&linbo3, [0.0, 128.0, 0.0]); // 128-deg Y-X LiNbO3 SAW cut
        let ln_chk = backward_wave_mpml_ratio(&c_ln_ycut, 4640.0, [0.0, 0.0, 1.0]);
        assert!(ln_chk.has_backward_wave, "128Y-X LiNbO3 must exhibit backward wave: {ln_chk:?}");
        assert!(ln_chk.recommended_mpml_ratio > 0.01);
    }
}
