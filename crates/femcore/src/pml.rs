// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Perfectly matched layers (PML) and static Infinite Elements (IE) by coordinate-free complex and
//! real coordinate stretching across Cartesian, Cylindrical, and Spherical topologies in 2D planar,
//! 2D axisymmetric, and 3D.
//!
//! Coordinate conventions:
//! * **Spatial inputs (`x`, `bounds`, `center`, `axis`)** are expressed in the 3D wafer frame
//!   `(x, y, z)` (`femcore::section::wafer_point(pt, dim)` maps a 2D mesh point `[p0, p1]` to
//!   `[p0, 0.0, p1]`).
//! * **Returned Jacobian tensors (`PmlJacobian::j`, `PmlJacobian::j_inv`)** are expressed in the
//!   solver's local axes: `(x, y, z)` in 3D (`dim == 3`), and `section` cross-section axes
//!   `(x/r, z, y/theta)` in 2D (`dim == 2`), matching `section::plane_tensor`, `section::local_tensor`,
//!   `nyx-ems::axisym` `(r, z, theta)`, and 2D mesh shape-function gradients `[d/dx_0, d/dx_1]`.
//! * Under `nyx`'s `exp(+i omega t)` / `exp(-i k u)` convention, outgoing waves are attenuated by
//!   `s(t) = kappa(t) - i * psi(t)` (`kappa >= 1`, `psi >= 0`), and transformation-optics pullbacks
//!   (`pullback_maxwell`) use non-conjugated complex-symmetric transposes (`J_s^{-T}`, `J_s^T`).

use crate::linalg::csr::{c, C};
use crate::mesh::SimplexMesh;
use crate::props::{Pml, PmlCoord, PmlDir, PmlProfile, PmlStabilization};

/// Maximum allowed rational curvature `c_rat < 1` so `1 - c_rat * t` stays strictly positive on `[0, 1]`.
const MAX_RATIONAL_CURVATURE: f64 = 0.999_999;

/// 8-point Gauss-Legendre nodes and weights mapped onto `[0, 1]`.
const GL8_NODES_01: [f64; 8] = [
    0.019_855_071_751_231_884,
    0.101_666_761_293_186_63,
    0.237_233_795_041_835_5,
    0.408_282_678_752_175_1,
    0.591_717_321_247_825,
    0.762_766_204_958_164_5,
    0.898_333_238_706_813_4,
    0.980_144_928_248_768_1,
];

const GL8_WEIGHTS_01: [f64; 8] = [
    0.050_614_268_145_188_13,
    0.111_190_517_226_687_23,
    0.156_853_322_938_943_64,
    0.181_341_891_689_181,
    0.181_341_891_689_181,
    0.156_853_322_938_943_64,
    0.111_190_517_226_687_23,
    0.050_614_268_145_188_13,
];

/// Pointwise coordinate-stretching Jacobian `J_s = d x_tilde / d x` at a quadrature point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PmlJacobian {
    /// Metric Jacobian `J_s` in `(x, y, z)` for `dim == 3`, or `(x/r, z, y/theta)` for `dim == 2`.
    pub j: [[C; 3]; 3],
    /// Exact analytical inverse `J_s^{-1}` in the same basis as `j`.
    pub j_inv: [[C; 3]; 3],
    /// Complex determinant `det(J_s)` (`s_0 * s_1` in 2D planar; `det(J_rz) * s_theta` in 2D axisym; `det(J_3x3)` in 3D).
    pub det: C,
    /// Azimuthal/circumferential stretch (`r_tilde / r` in 2D axisym & Cylindrical; `R_tilde / R` in Spherical; `1` in planar/3D Cartesian).
    pub s_theta: C,
    /// Normalized penetration depth in `[0.0, 1.0]` (`0` at physical interface, `1` at outer PML boundary).
    pub s_dist: f64,
}

impl PmlJacobian {
    /// Unstretched identity Jacobian (`J_s = I`, `det = 1`, `s_theta = 1`, `s_dist = 0`).
    pub fn identity() -> Self {
        let one = c(1.0, 0.0);
        let zero = c(0.0, 0.0);
        let id = [[one, zero, zero], [zero, one, zero], [zero, zero, one]];
        PmlJacobian { j: id, j_inv: id, det: one, s_theta: one, s_dist: 0.0 }
    }
}

impl Default for PmlJacobian {
    fn default() -> Self {
        Self::identity()
    }
}

/// Closed-form dimensionless PML calibration parameters (`PML_FEATURE_REQUEST.md` §2 Verdict 4).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AutoScaleResult {
    /// Recommended physical PML thickness `0.5 * lambda_max`.
    pub d_auto: f64,
    /// Imaginary damping scale `psi_max` (`Pml::strength`).
    pub psi_max: f64,
    /// Real evanescent dilation peak `kappa_max` (`1.0 .. 5.0`, capped at `2.0` for eigenproblems).
    pub kappa_max: f64,
    /// Polynomial profile order `m` (typically `2.0`).
    pub order: f64,
    /// Recommended structured layer count for linear (`P1` / `N1`) elements.
    pub layers_p1: usize,
    /// Recommended structured layer count for quadratic (`P2` / `N2`) elements.
    pub layers_p2: usize,
    /// Geometric layer thickness grading ratio `w_{N_L} / w_1` (`2.0`).
    pub layer_grading: f64,
}

impl AutoScaleResult {
    /// Recommended `curvature` (`c_rat` in `[0.85, 0.95]`) when `PmlProfile::Rational` is used.
    pub fn rational_curvature(&self) -> f64 {
        (1.0 - 1.0 / self.kappa_max.max(6.666_666_666_666_667)).clamp(0.85, 0.95)
    }
}

/// Evaluates `\Phi_m(z) = \int_0^1 s^m / (1 - z s)^2 ds` for `z = c_rat * t` in `[0, 1)` and `m >= 0`.
///
/// Uses an all-positive power series for small `z < 0.05` (eliminating `1 / z^{m+1}` catastrophic
/// cancellation near `c_rat * t -> 0`), exact closed-form antiderivatives for `m = 0, 1, 2, 3`,
/// and a pole-desingularized 8-point Gauss-Legendre rule for general real `m > 0`.
pub fn rational_moment_integral(m: f64, z: f64) -> f64 {
    let m = m.max(0.0);
    let z = z.clamp(0.0, MAX_RATIONAL_CURVATURE);
    if (m - 0.0).abs() <= 1e-14 {
        return 1.0 / (1.0 - z);
    }
    if z < 0.05 {
        let mut sum = 0.0;
        let mut zn = 1.0;
        for n in 0..18 {
            let nf = n as f64;
            sum += ((nf + 1.0) / (m + nf + 1.0)) * zn;
            zn *= z;
        }
        return sum;
    }
    let omz = 1.0 - z;
    if (m - 1.0).abs() <= 1e-14 {
        return (z / omz + omz.ln()) / (z * z);
    }
    if (m - 2.0).abs() <= 1e-14 {
        return (z * (2.0 - z) / omz + 2.0 * omz.ln()) / (z * z * z);
    }
    if (m - 3.0).abs() <= 1e-14 {
        return (z / omz + 3.0 * omz.ln() + 2.0 * z + 0.5 * z * z) / (z * z * z * z);
    }
    // Pole-free Mobius + quadratic substitution: `s = u^2 / (1 - z * (1 - u^2))`, `u` in `[0, 1]`.
    // Geometrically graded panels toward `u = 0` keep the complex branch point `+- i * sqrt((1-z)/z)`
    // outside the Bernstein ellipse `rho >= 7.5` on every panel (`rho^{-16} < 1e-14` up to `z = 0.99`).
    const EDGES: [f64; 6] = [0.0, 0.04, 0.12, 0.28, 0.56, 1.0];
    let mut acc = 0.0;
    for panel in 0..5 {
        let u0 = EDGES[panel];
        let du = EDGES[panel + 1] - u0;
        for q in 0..8 {
            let u = u0 + du * GL8_NODES_01[q];
            let tau = u * u;
            let denom = 1.0 - z * (1.0 - tau);
            acc += du * GL8_WEIGHTS_01[q] * 2.0 * u.powf(2.0 * m + 1.0) / denom.powf(m);
        }
    }
    acc / omz
}

/// Evaluates `(s(t), \int_0^t s(\xi) d\xi, \psi(t))` for normalized PML depth `t` in `[0, 1]`.
#[inline]
fn eval_profile_1d(spec: &Pml, t: f64) -> (C, C, f64) {
    let t = t.clamp(0.0, 1.0);
    let m = spec.order.max(0.0);
    let tm = if spec.order == 0.0 { 1.0 } else { t.powf(m) };
    match spec.profile {
        PmlProfile::KappaPolynomial { kappa_max } => {
            let psi = spec.strength * tm;
            let s_re = if (kappa_max - 1.0).abs() == 0.0 {
                spec.re
            } else {
                1.0 + (kappa_max - 1.0) * tm
            };
            let int_tm = if spec.order == 0.0 { t } else { t.powf(m + 1.0) / (m + 1.0) };
            let int_re = if (kappa_max - 1.0).abs() == 0.0 {
                spec.re * t
            } else {
                t + (kappa_max - 1.0) * int_tm
            };
            let int_im = -spec.strength * int_tm;
            (c(s_re, -psi), c(int_re, int_im), psi)
        }
        PmlProfile::Rational { curvature } => {
            let c_rat = curvature.clamp(0.0, MAX_RATIONAL_CURVATURE);
            let denom = 1.0 - c_rat * t;
            let inv_d2 = 1.0 / (denom * denom);
            let psi = spec.strength * tm * inv_d2;
            let s_re = inv_d2;
            let int_re = t / denom;
            let int_im = if t == 0.0 {
                0.0
            } else {
                let t_pow = if spec.order == 0.0 { t } else { t.powf(m + 1.0) };
                -spec.strength * t_pow * rational_moment_integral(m, c_rat * t)
            };
            (c(s_re, -psi), c(int_re, int_im), psi)
        }
    }
}

/// Adiabatic mechanical stiffness loss factor `1 + i * eta_max * s_dist^2` for `PmlStabilization::Adiabatic`
/// (`1 + 0i` for all other stabilization modes). Multiplying elastic stiffness `C` by `(1 + i * eta)`
/// damps outgoing waves strictly by the sign of group velocity `v_{g,n} > 0` under `exp(+i omega t)`.
pub fn adiabatic_loss_factor(spec: &Pml, s_dist: f64) -> C {
    match spec.stabilization {
        PmlStabilization::Adiabatic { eta_max } if eta_max != 0.0 => {
            let t = s_dist.clamp(0.0, 1.0);
            c(1.0, eta_max * t * t)
        }
        _ => c(1.0, 0.0),
    }
}

/// Normalizes a bounding box into wafer `(x, y, z)` coordinates, automatically converting a raw
/// 2D `SimplexMesh::bbox()` (`lo[2] == +INF, hi[2] == -INF`) into `([lo[0], 0, lo[1]], [hi[0], 0, hi[1]])`.
#[inline]
fn resolve_wafer_bounds(bounds: ([f64; 3], [f64; 3]), dim: usize) -> ([f64; 3], [f64; 3]) {
    let (lo, hi) = bounds;
    if dim == 2 && lo[2] > hi[2] && lo[1] <= hi[1] {
        ([lo[0], 0.0, lo[1]], [hi[0], 0.0, hi[1]])
    } else {
        (lo, hi)
    }
}

/// Permutes a 3x3 complex tensor from wafer `(x, y, z)` into `section` 2D cross-section order `(x, z, y)`.
#[inline]
fn permute_wafer_to_2d(m: [[C; 3]; 3]) -> [[C; 3]; 3] {
    const P: [usize; 3] = [0, 2, 1];
    let mut out = [[c(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = m[P[i]][P[j]];
        }
    }
    out
}

/// Unified coordinate-free curvilinear Jacobian evaluator driven by a 1D ray-stretching law
/// `ray_map(t, length) -> (s_u, delta_u_tilde, psi)`.
#[allow(clippy::too_many_arguments)]
fn evaluate_curvilinear_jacobian(
    coord: PmlCoord,
    cartesian_dirs: [PmlDir; 3],
    bounds_opt: Option<([f64; 3], [f64; 3])>,
    mpml_ratio: f64,
    x: [f64; 3],
    dim: usize,
    axisym: bool,
    ray_map: impl Fn(f64, f64) -> (C, C, f64),
) -> PmlJacobian {
    let one = c(1.0, 0.0);
    let zero = c(0.0, 0.0);

    match coord {
        PmlCoord::Cartesian => {
            let Some(raw_bounds) = bounds_opt else {
                return PmlJacobian::identity();
            };
            let (lo, hi) = resolve_wafer_bounds(raw_bounds, dim);
            // Active wafer axes: in 2D `(x, z)` only (`y` is out-of-plane); in 3D `(x, y, z)`.
            let wafer_dirs = if dim == 2 {
                [cartesian_dirs[0], PmlDir::Off, cartesian_dirs[2]]
            } else {
                cartesian_dirs
            };

            let mut s_prim = [one; 3];
            let mut psi_axis = [0.0_f64; 3];
            let mut t_axis = [0.0_f64; 3];
            let mut delta_x0 = zero;

            for a in 0..3 {
                let len = (hi[a] - lo[a]).max(1e-300);
                let t = match wafer_dirs[a] {
                    PmlDir::Off => continue,
                    PmlDir::Plus => ((x[a] - lo[a]) / len).clamp(0.0, 1.0),
                    PmlDir::Minus => ((hi[a] - x[a]) / len).clamp(0.0, 1.0),
                };
                let (s_val, du, psi) = ray_map(t, len);
                s_prim[a] = s_val;
                psi_axis[a] = psi;
                t_axis[a] = t;
                if a == 0 {
                    delta_x0 = du;
                }
            }

            let s_dist = t_axis[0].max(t_axis[1]).max(t_axis[2]);
            let mut s_diag = s_prim;
            if mpml_ratio != 0.0 {
                if dim == 2 {
                    // Tangential coupling between wafer x (0) and wafer z (2)
                    s_diag[0] += c(0.0, -mpml_ratio * psi_axis[2]);
                    s_diag[2] += c(0.0, -mpml_ratio * psi_axis[0]);
                    if axisym {
                        s_diag[1] += c(0.0, -mpml_ratio * (psi_axis[0] + psi_axis[2]));
                    }
                } else {
                    let psi_sum = psi_axis[0] + psi_axis[1] + psi_axis[2];
                    for a in 0..3 {
                        s_diag[a] += c(0.0, -mpml_ratio * (psi_sum - psi_axis[a]));
                    }
                }
            }

            let mut s_theta = one;
            if dim == 2 && axisym {
                if wafer_dirs[0] != PmlDir::Off {
                    let r = x[0];
                    let len = (hi[0] - lo[0]).max(1e-300);
                    let r_scale = len.max(lo[0].abs()).max(hi[0].abs()).max(1e-300);
                    let s_th_prim = if r.abs() > 1e-15 * r_scale {
                        let r_tilde = match wafer_dirs[0] {
                            PmlDir::Plus => c(lo[0], 0.0) + delta_x0,
                            PmlDir::Minus => c(hi[0], 0.0) - delta_x0,
                            PmlDir::Off => c(r, 0.0),
                        };
                        r_tilde / r
                    } else {
                        s_prim[0]
                    };
                    s_diag[1] = s_th_prim + c(0.0, -mpml_ratio * (psi_axis[0] + psi_axis[2]));
                }
                s_theta = s_diag[1];
            }

            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            for a in 0..3 {
                j_w[a][a] = s_diag[a];
                j_inv_w[a][a] = one / s_diag[a];
            }
            let det = if dim == 2 && !axisym {
                s_diag[0] * s_diag[2]
            } else {
                s_diag[0] * s_diag[1] * s_diag[2]
            };
            let (j, j_inv) = if dim == 2 {
                (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w))
            } else {
                (j_w, j_inv_w)
            };
            PmlJacobian { j, j_inv, det, s_theta, s_dist }
        }

        PmlCoord::Cylindrical {
            center,
            axis,
            r_in,
            r_thickness,
            radial_on,
            axial,
            z_in,
            z_thickness,
        } => {
            let ax_norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            let a_hat = if ax_norm > 1e-30 {
                [axis[0] / ax_norm, axis[1] / ax_norm, axis[2] / ax_norm]
            } else if dim == 2 && !axisym {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            };

            let d_vec = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
            let zeta = d_vec[0] * a_hat[0] + d_vec[1] * a_hat[1] + d_vec[2] * a_hat[2];
            let r_vec = [
                d_vec[0] - zeta * a_hat[0],
                d_vec[1] - zeta * a_hat[1],
                d_vec[2] - zeta * a_hat[2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();

            let (mut s_r, mut s_phi, mut psi_r, mut t_r) = (one, one, 0.0_f64, 0.0_f64);
            if radial_on && r_thickness > 0.0 && r >= r_in {
                t_r = ((r - r_in) / r_thickness).clamp(0.0, 1.0);
                let (sr_val, delta_r, psi_val) = ray_map(t_r, r_thickness);
                s_r = sr_val;
                psi_r = psi_val;
                let r_scale = r_thickness.max(r_in.abs()).max(1e-300);
                s_phi = if r > 1e-15 * r_scale {
                    (c(r_in, 0.0) + delta_r) / r
                } else {
                    sr_val
                };
            }

            let (mut s_zeta, mut psi_z, mut t_z) = (one, 0.0_f64, 0.0_f64);
            if axial != PmlDir::Off && z_thickness > 0.0 {
                let raw_tz = match axial {
                    PmlDir::Plus => (zeta - z_in) / z_thickness,
                    PmlDir::Minus => (z_in - zeta) / z_thickness,
                    PmlDir::Off => 0.0,
                };
                if raw_tz >= 0.0 {
                    t_z = raw_tz.clamp(0.0, 1.0);
                    let (sz_val, _, psi_val) = ray_map(t_z, z_thickness);
                    s_zeta = sz_val;
                    psi_z = psi_val;
                }
            }

            if mpml_ratio != 0.0 {
                if dim == 2 && !axisym {
                    // 2D planar annular ring: axis is out-of-plane `y` (keep `s_zeta = 1`), damp in-plane `phi`.
                    s_phi += c(0.0, -mpml_ratio * psi_r);
                } else {
                    s_r += c(0.0, -mpml_ratio * psi_z);
                    s_zeta += c(0.0, -mpml_ratio * psi_r);
                    s_phi += c(0.0, -mpml_ratio * (psi_r + psi_z));
                }
            }

            let s_dist = t_r.max(t_z);
            let inv_sr = one / s_r;
            let inv_sphi = one / s_phi;
            let inv_szeta = one / s_zeta;

            // Orthogonal rank-1 projectors:
            // P_a = a_hat a_hat^T, P_r = r_hat r_hat^T, P_phi = I - P_a - P_r
            // J = s_phi * (I - P_a) + (s_r - s_phi) * P_r + s_zeta * P_a
            // Note: as r -> 0^+, s_r == s_phi identically, so `(s_r - s_phi) * P_r == 0` with zero division by r = 0!
            let r_hat = if r > 1e-30 {
                [r_vec[0] / r, r_vec[1] / r, r_vec[2] / r]
            } else {
                [0.0, 0.0, 0.0]
            };

            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            let ds_r = s_r - s_phi;
            let dinv_r = inv_sr - inv_sphi;
            for p in 0..3 {
                for q in 0..3 {
                    let id_pq = if p == q { 1.0 } else { 0.0 };
                    let pa = a_hat[p] * a_hat[q];
                    let pr = r_hat[p] * r_hat[q];
                    j_w[p][q] = s_phi * (id_pq - pa) + ds_r * pr + s_zeta * pa;
                    j_inv_w[p][q] = inv_sphi * (id_pq - pa) + dinv_r * pr + inv_szeta * pa;
                }
            }

            let s_theta = if dim == 2 && !axisym { one } else { s_phi };
            let det = if dim == 2 && !axisym {
                // In 2D planar annular ring, out-of-plane `y` (`a_hat = [0, 1, 0]`) has stretch `s_zeta = 1`.
                s_r * s_phi
            } else {
                s_r * s_phi * s_zeta
            };
            let (j, j_inv) = if dim == 2 {
                (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w))
            } else {
                (j_w, j_inv_w)
            };
            PmlJacobian { j, j_inv, det, s_theta, s_dist }
        }

        PmlCoord::Spherical { center, r_in, thickness } => {
            if thickness <= 0.0 {
                return PmlJacobian::identity();
            }
            let r_vec = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
            let big_r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            if big_r < r_in {
                return PmlJacobian::identity();
            }

            let t = ((big_r - r_in) / thickness).clamp(0.0, 1.0);
            let (s_big_r, delta_r, psi_r) = ray_map(t, thickness);
            let r_scale = thickness.max(r_in.abs()).max(1e-300);
            // Exact spherical angular stretch `s_Omega(R) = R_tilde(R) / R`.
            // In 2D axisymmetric `(r >= 0, z)`, `r_tilde / r = (R_tilde * (r / R)) / r = R_tilde / R = s_Omega(R)`,
            // which is completely analytic across the symmetry axis `r -> 0^+` without ever dividing by `r`!
            let mut s_omega = if big_r > 1e-15 * r_scale {
                (c(r_in, 0.0) + delta_r) / big_r
            } else {
                s_big_r
            };
            if mpml_ratio != 0.0 && big_r > r_in {
                s_omega += c(0.0, -mpml_ratio * psi_r);
            }

            let inv_sr = one / s_big_r;
            let inv_som = one / s_omega;
            let r_hat = if big_r > 1e-30 {
                [r_vec[0] / big_r, r_vec[1] / big_r, r_vec[2] / big_r]
            } else {
                [0.0, 0.0, 0.0]
            };

            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            let ds = s_big_r - s_omega;
            let dinv = inv_sr - inv_som;
            for p in 0..3 {
                for q in 0..3 {
                    let id_pq = if p == q { 1.0 } else { 0.0 };
                    let pr = r_hat[p] * r_hat[q];
                    j_w[p][q] = s_omega * id_pq + ds * pr;
                    j_inv_w[p][q] = inv_som * id_pq + dinv * pr;
                }
            }

            if dim == 2 && !axisym {
                // 2D planar radial shell (out-of-plane `y` unstretched):
                j_w[1][1] = one;
                j_inv_w[1][1] = one;
                let (j, j_inv) = (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w));
                return PmlJacobian {
                    j,
                    j_inv,
                    det: s_big_r * s_omega,
                    s_theta: one,
                    s_dist: t,
                };
            }

            let det = s_big_r * s_omega * s_omega;
            let (j, j_inv) = if dim == 2 {
                (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w))
            } else {
                (j_w, j_inv_w)
            };
            PmlJacobian { j, j_inv, det, s_theta: s_omega, s_dist: t }
        }
    }
}

/// Evaluates the pointwise complex coordinate-stretching Jacobian `J_s(x)` of a PML domain.
///
/// * `spec`: domain PML configuration (`Cartesian`, `Cylindrical`, or `Spherical`; `KappaPolynomial` or `Rational`).
/// * `fallback_bounds`: optional `(lo, hi)` domain bounding box for `Cartesian` slabs when `spec.bounds` is `None`.
/// * `x`: evaluation point in 3D wafer coordinates `(x, y, z)` (`femcore::section::wafer_point(pt, dim)`).
/// * `dim`: mesh dimension (`2` or `3`).
/// * `axisym`: `true` for 2D axisymmetric `(r >= 0, z)` models (`x = [r, 0.0, z]`), `false` otherwise.
pub fn pml_jacobian(
    spec: &Pml,
    fallback_bounds: Option<([f64; 3], [f64; 3])>,
    x: [f64; 3],
    dim: usize,
    axisym: bool,
) -> PmlJacobian {
    let mpml_ratio = match spec.stabilization {
        PmlStabilization::Mpml { ratio } => ratio,
        _ => 0.0,
    };
    let bounds = spec.bounds.or(fallback_bounds);
    evaluate_curvilinear_jacobian(
        spec.coord,
        [spec.x, spec.y, spec.z],
        bounds,
        mpml_ratio,
        x,
        dim,
        axisym,
        |t, length| {
            let (s_val, int_s, psi) = eval_profile_1d(spec, t);
            (s_val, int_s * length, psi)
        },
    )
}

/// Transformation-optics pullback of Maxwell permittivity and inverse permeability tensors under
/// complex coordinate stretching `jac`:
///
/// `eps_tilde = det(J_s) * J_s^{-1} * eps * J_s^{-T}`
/// `mu_inv_tilde = (1 / det(J_s)) * J_s^T * mu_inv * J_s`
///
/// Uses non-conjugated transposes so symmetric `eps` and `mu_inv` map to complex-symmetric tensors
/// in 3D `(x, y, z)`, 2D waveguide `(x, z, y)`, and 2D axisymmetric `(r, z, theta)` (where `det`
/// already incorporates `s_theta = r_tilde / r` and `j[2][2] == s_theta`).
pub fn pullback_maxwell(
    eps: [[C; 3]; 3],
    mu_inv: [[C; 3]; 3],
    jac: &PmlJacobian,
) -> ([[C; 3]; 3], [[C; 3]; 3]) {
    let zero = c(0.0, 0.0);
    let inv_det = c(1.0, 0.0) / jac.det;

    // tmp_e = J_inv * eps ; eps_tilde = det * tmp_e * J_inv^T
    let mut tmp_e = [[zero; 3]; 3];
    let mut tmp_m = [[zero; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            let mut acc_e = zero;
            let mut acc_m = zero;
            for l in 0..3 {
                acc_e += jac.j_inv[i][l] * eps[l][k];
                acc_m += jac.j[l][i] * mu_inv[l][k];
            }
            tmp_e[i][k] = acc_e;
            tmp_m[i][k] = acc_m;
        }
    }

    let mut eps_tilde = [[zero; 3]; 3];
    let mut mu_inv_tilde = [[zero; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut acc_e = zero;
            let mut acc_m = zero;
            for k in 0..3 {
                acc_e += tmp_e[i][k] * jac.j_inv[j][k];
                acc_m += tmp_m[i][k] * jac.j[k][j];
            }
            eps_tilde[i][j] = jac.det * acc_e;
            mu_inv_tilde[i][j] = inv_det * acc_m;
        }
    }
    (eps_tilde, mu_inv_tilde)
}

/// Computes closed-form dimensionless PML parameters jointly calibrated across a PML material group
/// (`lambda_min <= lambda_max`) per `docs/PML_FEATURE_REQUEST.md` §2 Verdict 4:
///
/// * `psi_max` and `kappa_max` calibrate to `lambda_max` (fastest wave / lowest refractive index)
///   so continuous reflection satisfies `R_cont <= target_r0` for every mode in the slab/shell.
/// * `layers_p1` and `layers_p2` calibrate to `lambda_min` (slowest wave / highest refractive index)
///   from complex WKB arc-length so the shortest wavelength never aliases.
pub fn auto_scale_params(
    lambda_min: f64,
    lambda_max: f64,
    d_actual_opt: Option<f64>,
    order_m: f64,
    target_r0: f64,
    is_eigen: bool,
    is_rational: bool,
) -> AutoScaleResult {
    let lam_max = lambda_max.max(1e-30);
    let lam_min = lambda_min.clamp(1e-30, lam_max);
    let m = order_m.max(0.0);
    let r0 = target_r0.clamp(1e-15, 0.5);
    let ln_inv_r0 = (1.0 / r0).ln();

    let d_auto = 0.5 * lam_max;
    let d = d_actual_opt.filter(|&v| v > 0.0 && v.is_finite()).unwrap_or(d_auto);
    let d_over_lmax = (d / lam_max).max(1e-12);
    let d_over_lmin = (d / lam_min).max(1e-12);
    let four_pi_dlmax = 4.0 * std::f64::consts::PI * d_over_lmax;

    let bar_kappa = (ln_inv_r0 / four_pi_dlmax).max(1.0);
    let layer_grading = 2.0;

    if is_rational {
        let c_rat = (1.0 - 1.0 / (4.0 * bar_kappa).clamp(6.666_666_666_666_667, 20.0)).clamp(0.85, 0.95);
        let phi_m = rational_moment_integral(m, c_rat).max(1e-12);
        let psi_max = (ln_inv_r0 / (four_pi_dlmax * phi_m)).clamp(0.25, 10.0);
        let kappa_max = 1.0 / (1.0 - c_rat);
        let ratio = (lam_max / lam_min).sqrt();
        let layers_p2 = 6usize.max((6.0 * (d_over_lmin / 0.5).max(1.0) * ratio.min(2.0)).ceil() as usize);
        let layers_p1 = 8usize.max((8.0 * (d_over_lmin / 0.5).max(1.0) * ratio.min(2.0)).ceil() as usize);
        return AutoScaleResult {
            d_auto,
            psi_max,
            kappa_max,
            order: m,
            layers_p1,
            layers_p2,
            layer_grading,
        };
    }

    let psi_max = ((m + 1.0) * ln_inv_r0 / four_pi_dlmax).clamp(1.5, 25.0);
    let kappa_cap = if is_eigen { 2.0 } else { 5.0 };
    let kappa_max = (1.0 + (m + 1.0) * (bar_kappa - 1.0)).clamp(1.0, kappa_cap);

    let wkb_arc = d_over_lmin * (bar_kappa * bar_kappa + (psi_max / (m + 1.0)).powi(2)).sqrt();
    let layers_p2 = 6usize.max((5.0 * wkb_arc).ceil() as usize);
    let layers_p1 = 10usize.max((10.0 * wkb_arc).ceil() as usize);

    AutoScaleResult {
        d_auto,
        psi_max,
        kappa_max,
        order: m,
        layers_p1,
        layers_p2,
        layer_grading,
    }
}

/// 1D rational Infinite Element stretching `(s_u(t), u_tilde(t) - u_0)` for normalized depth `t` in `[0, 1]`,
/// physical layer thickness `d`, and pole distance `d_pole`:
///
/// `u_tilde(t) - u_0 = (d_pole * t * d) / (d_pole + d - t * d)`
/// `s_u(t) = d u_tilde / d u = d_pole * (d_pole + d) / (d_pole + d - t * d)^2`
///
/// Satisfies `(d_pole + u_tilde(t) - u_0) / d_pole == 1 / (1 - t * d / (d_pole + d))` identically.
#[inline]
pub fn ie_rational_stretch_1d(t: f64, d: f64, d_pole: f64) -> (f64, f64) {
    let t = t.clamp(0.0, 1.0);
    let d_pos = d.max(1e-300);
    let dp = d_pole.max(1e-300);
    let denom = (dp + d_pos - t * d_pos).max(1e-300);
    let du = (dp * t * d_pos) / denom;
    let s_u = (dp * (dp + d_pos)) / (denom * denom);
    (s_u, du)
}

/// Pointwise real coordinate-stretching Jacobian for static Rational Infinite Elements (`nyx-dfs`
/// and `nyx-sms::static`) across Cartesian, Cylindrical, and Spherical topologies.
pub fn ie_rational_jacobian(
    spec: &Pml,
    fallback_bounds: Option<([f64; 3], [f64; 3])>,
    d_pole: f64,
    x: [f64; 3],
    dim: usize,
    axisym: bool,
) -> PmlJacobian {
    let bounds = spec.bounds.or(fallback_bounds);
    evaluate_curvilinear_jacobian(
        spec.coord,
        [spec.x, spec.y, spec.z],
        bounds,
        0.0,
        x,
        dim,
        axisym,
        |t, length| {
            let (s_u, du) = ie_rational_stretch_1d(t, length, d_pole);
            (c(s_u, 0.0), c(du, 0.0), 0.0)
        },
    )
}

/// 1D bounded exponential Infinite Element stretching `(s_u(t), u_tilde(t) - u_0)` with stretch
/// ratio `chi = exp(beta) >= 1` (typically `chi` in `[20, 50]` when paired with outer Order-1 Robin):
///
/// `u_tilde(t) - u_0 = d * (exp(beta * t) - 1) / beta`,   `s_u(t) = exp(beta * t)`.
#[inline]
pub fn ie_exponential_stretch_1d(t: f64, d: f64, chi: f64) -> (f64, f64) {
    let t = t.clamp(0.0, 1.0);
    let d_pos = d.max(1e-300);
    let chi = chi.max(1.0);
    let beta = chi.ln();
    if beta < 1e-6 {
        let bt = beta * t;
        let s_u = bt.exp();
        let du = d_pos * t * (1.0 + 0.5 * bt + bt * bt / 6.0 + bt * bt * bt / 24.0);
        (s_u, du)
    } else {
        let e_bt = (beta * t).exp();
        (e_bt, d_pos * (e_bt - 1.0) / beta)
    }
}

/// Pointwise real coordinate-stretching Jacobian for bounded Exponential Infinite Elements
/// (`chi = exp(beta)` typically in `[20, 50]`, gentler conditioning for iterative `AMG`+`PCG`).
pub fn ie_exponential_jacobian(
    spec: &Pml,
    fallback_bounds: Option<([f64; 3], [f64; 3])>,
    chi: f64,
    x: [f64; 3],
    dim: usize,
    axisym: bool,
) -> PmlJacobian {
    let bounds = spec.bounds.or(fallback_bounds);
    evaluate_curvilinear_jacobian(
        spec.coord,
        [spec.x, spec.y, spec.z],
        bounds,
        0.0,
        x,
        dim,
        axisym,
        |t, length| {
            let (s_u, du) = ie_exponential_stretch_1d(t, length, chi);
            (c(s_u, 0.0), c(du, 0.0), 0.0)
        },
    )
}

// ------------------------------------------------------------------ legacy cell-centroid helpers

/// Stretch factors along the two axes of the mesh (`dirs`, the layer directions of each) of a cell centre for a PML specification and the
/// bounding box of its domain. Honors both `KappaPolynomial` (preserving exact legacy `re` behavior) and `Rational` profiles.
pub fn stretch(p: &Pml, dirs: [PmlDir; 2], centre: [f64; 2], lo: [f64; 2], hi: [f64; 2]) -> (C, C) {
    let one = |dir: PmlDir, x: f64, lo: f64, hi: f64| -> C {
        let t = match dir {
            PmlDir::Off => return c(1.0, 0.0),
            PmlDir::Plus => (x - lo) / (hi - lo).max(1e-300),
            PmlDir::Minus => (hi - x) / (hi - lo).max(1e-300),
        };
        eval_profile_1d(p, t).0
    };
    (one(dirs[0], centre[0], lo[0], hi[0]), one(dirs[1], centre[1], lo[1], hi[1]))
}

/// The complex radius `r~ = int s_r dr` (from the interface of the layer) of a point at radius `r` for a layer whose `x` direction is the
/// radial one of an axisymmetric model, `lo` and `hi` being the radial extent of the layer's domain; `r` itself when `p.x` is off.
/// The cylindrical stretch of the axisymmetric kernels is `s_theta = r~ / r` (ADR 46 for `ems`, ADR 69 for `sms`).
pub fn radial_coordinate(p: &Pml, r: f64, lo: f64, hi: f64) -> C {
    let len = (hi - lo).max(1e-300);
    match p.x {
        PmlDir::Off => c(r, 0.0),
        PmlDir::Plus => {
            let t = ((r - lo) / len).clamp(0.0, 1.0);
            if matches!(p.profile, PmlProfile::KappaPolynomial { kappa_max } if (kappa_max - 1.0).abs() == 0.0) {
                c(lo + p.re * (r - lo), -p.strength * len * t.powf(p.order + 1.0) / (p.order + 1.0))
            } else {
                c(lo, 0.0) + eval_profile_1d(p, t).1 * len
            }
        }
        PmlDir::Minus => {
            let t = ((hi - r) / len).clamp(0.0, 1.0);
            if matches!(p.profile, PmlProfile::KappaPolynomial { kappa_max } if (kappa_max - 1.0).abs() == 0.0) {
                c(hi - p.re * (hi - r), p.strength * len * t.powf(p.order + 1.0) / (p.order + 1.0))
            } else {
                c(hi, 0.0) - eval_profile_1d(p, t).1 * len
            }
        }
    }
}

/// Complex stretch factors `(sx, sy, sz)` of every cell of a 2D or 3D mesh for the PML domains of
/// `pml` (cells of other domains get 1). Evaluates at the cell centre against explicit `Pml::bounds`
/// when provided, falling back to the bounding box of the domain's cells.
pub fn cell_stretches(mesh: &SimplexMesh, pml: &std::collections::BTreeMap<i32, Pml>) -> Vec<[C; 3]> {
    cell_stretches_by(mesh, pml, crate::section::mesh_pml_dirs)
}

fn cell_stretches_by(
    mesh: &SimplexMesh,
    pml: &std::collections::BTreeMap<i32, Pml>,
    dirs: impl Fn(&Pml, usize) -> [PmlDir; 3],
) -> Vec<[C; 3]> {
    let d = mesh.dim;
    let mut boxes: std::collections::BTreeMap<i32, ([f64; 3], [f64; 3])> = Default::default();
    for ci in 0..mesh.n_cells() {
        let l = mesh.cell_labels[ci];
        if pml.get(&l).is_some_and(Pml::is_active) {
            let e = boxes.entry(l).or_insert(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
            for &v in mesh.cell(ci) {
                let pt = mesh.point(v as usize);
                for k in 0..d {
                    e.0[k] = e.0[k].min(pt[k]);
                    e.1[k] = e.1[k].max(pt[k]);
                }
            }
        }
    }
    (0..mesh.n_cells())
        .map(|ci| {
            let mut s = [c(1.0, 0.0); 3];
            if let Some(p) = pml.get(&mesh.cell_labels[ci]) {
                let mut cen = [0.0; 3];
                for &v in mesh.cell(ci) {
                    for k in 0..d {
                        cen[k] += mesh.point(v as usize)[k] / (d as f64 + 1.0);
                    }
                }
                match p.coord {
                    PmlCoord::Cartesian => {
                        let mesh_box = p.bounds.map(|(blo, bhi)| {
                            if d == 2 {
                                ([blo[0], blo[2], 0.0], [bhi[0], bhi[2], 0.0])
                            } else {
                                (blo, bhi)
                            }
                        });
                        if let Some((lo, hi)) = mesh_box.or_else(|| boxes.get(&mesh.cell_labels[ci]).copied()) {
                            for (k, dir) in dirs(p, d).into_iter().enumerate().take(d) {
                                s[k] = stretch_one(p, dir, cen[k], lo[k], hi[k]);
                            }
                        }
                    }
                    PmlCoord::Cylindrical { .. } | PmlCoord::Spherical { .. } => {
                        let x_w = crate::section::wafer_point(&cen[..d], d);
                        let jac = pml_jacobian(p, None, x_w, d, false);
                        for k in 0..d {
                            s[k] = jac.j[k][k];
                        }
                    }
                }
            }
            s
        })
        .collect()
}

fn stretch_one(p: &Pml, dir: PmlDir, x: f64, lo: f64, hi: f64) -> C {
    let t = match dir {
        PmlDir::Off => return c(1.0, 0.0),
        PmlDir::Plus => (x - lo) / (hi - lo).max(1e-300),
        PmlDir::Minus => (hi - x) / (hi - lo).max(1e-300),
    };
    eval_profile_1d(p, t).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::rotation;

    const GL16_NODES_01: [f64; 16] = [
        0.005_299_532_504_175_033,
        0.027_712_488_463_383_71,
        0.067_184_398_806_084_12,
        0.122_297_795_822_498_48,
        0.191_061_877_798_678_12,
        0.270_991_611_171_386_3,
        0.359_198_224_610_370_54,
        0.452_493_745_081_181_3,
        0.547_506_254_918_818_7,
        0.640_801_775_389_629_5,
        0.729_008_388_828_613_7,
        0.808_938_122_201_321_9,
        0.877_702_204_177_501_5,
        0.932_815_601_193_915_9,
        0.972_287_511_536_616_3,
        0.994_700_467_495_825,
    ];

    const GL16_WEIGHTS_01: [f64; 16] = [
        0.013_576_229_705_877_047,
        0.031_126_761_969_323_946,
        0.047_579_255_841_246_39,
        0.062_314_485_627_766_94,
        0.074_797_994_408_288_37,
        0.084_578_259_697_501_27,
        0.091_301_707_522_461_79,
        0.094_725_305_227_534_25,
        0.094_725_305_227_534_25,
        0.091_301_707_522_461_79,
        0.084_578_259_697_501_27,
        0.074_797_994_408_288_37,
        0.062_314_485_627_766_94,
        0.047_579_255_841_246_39,
        0.031_126_761_969_323_946,
        0.013_576_229_705_877_047,
    ];

    fn c_abs(z: C) -> f64 {
        z.re.hypot(z.im)
    }

    fn mat_diff_max(a: &[[C; 3]; 3], b: &[[C; 3]; 3]) -> f64 {
        let mut mx = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                mx = mx.max(c_abs(a[i][j] - b[i][j]));
            }
        }
        mx
    }

    fn conj_by_rot(q: &[[f64; 3]; 3], j: &[[C; 3]; 3]) -> [[C; 3]; 3] {
        // Computes Q^T * J * Q
        let mut tmp = [[c(0.0, 0.0); 3]; 3];
        for r in 0..3 {
            for col in 0..3 {
                let mut acc = c(0.0, 0.0);
                for k in 0..3 {
                    acc += j[k][col] * q[k][r];
                }
                tmp[r][col] = acc;
            }
        }
        let mut out = [[c(0.0, 0.0); 3]; 3];
        for r in 0..3 {
            for col in 0..3 {
                let mut acc = c(0.0, 0.0);
                for k in 0..3 {
                    acc += tmp[r][k] * q[k][col];
                }
                out[r][col] = acc;
            }
        }
        out
    }

    fn rot_vec(q: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
        [
            q[0][0] * v[0] + q[0][1] * v[1] + q[0][2] * v[2],
            q[1][0] * v[0] + q[1][1] * v[1] + q[1][2] * v[2],
            q[2][0] * v[0] + q[2][1] * v[1] + q[2][2] * v[2],
        ]
    }

    fn assert_j_times_jinv_is_identity(jac: &PmlJacobian) {
        for i in 0..3 {
            for j in 0..3 {
                let mut acc = c(0.0, 0.0);
                for k in 0..3 {
                    acc += jac.j[i][k] * jac.j_inv[k][j];
                }
                let want = if i == j { c(1.0, 0.0) } else { c(0.0, 0.0) };
                assert!(c_abs(acc - want) < 1e-13, "J * J_inv[{i}][{j}] = {acc:?}");
            }
        }
    }

    #[test]
    fn spherical_and_cylindrical_jacobians_are_rotationally_covariant() {
        let q = rotation([37.0, 53.0, 19.0]);
        let center = [0.25, -0.4, 0.15];
        let dx = [0.6, -0.35, 0.85]; // |dx| = 1.0977 > r_in = 0.8
        let x = [center[0] + dx[0], center[1] + dx[1], center[2] + dx[2]];
        let dx_rot = rot_vec(&q, dx);
        let x_rot = [center[0] + dx_rot[0], center[1] + dx_rot[1], center[2] + dx_rot[2]];

        // 1. Spherical PML with KappaPolynomial and Rational profiles
        for profile in [
            PmlProfile::KappaPolynomial { kappa_max: 3.2 },
            PmlProfile::Rational { curvature: 0.88 },
        ] {
            let spec = Pml {
                coord: PmlCoord::Spherical { center, r_in: 0.8, thickness: 0.5 },
                profile,
                strength: 5.5,
                order: 2.0,
                ..Pml::default()
            };
            let j0 = pml_jacobian(&spec, None, x, 3, false);
            let j_rot = pml_jacobian(&spec, None, x_rot, 3, false);
            assert_j_times_jinv_is_identity(&j0);
            assert_j_times_jinv_is_identity(&j_rot);
            let back = conj_by_rot(&q, &j_rot.j);
            assert!(mat_diff_max(&back, &j0.j) < 1e-13, "spherical Q^T J Q mismatch");
            assert!(c_abs(j_rot.det - j0.det) < 1e-13);
            assert!(c_abs(j_rot.s_theta - j0.s_theta) < 1e-13);
            assert!((j_rot.s_dist - j0.s_dist).abs() < 1e-14);
        }

        // 2. Cylindrical PML (corner annulus: radial + axial active simultaneously) rotated in SO(3)
        let axis0 = [0.0, 0.0, 1.0];
        let axis_rot = rot_vec(&q, axis0);
        let cyl0 = Pml {
            coord: PmlCoord::Cylindrical {
                center,
                axis: axis0,
                r_in: 0.5,
                r_thickness: 0.4,
                radial_on: true,
                axial: PmlDir::Plus,
                z_in: 0.6,
                z_thickness: 0.4,
            },
            profile: PmlProfile::KappaPolynomial { kappa_max: 2.4 },
            strength: 4.2,
            order: 2.0,
            stabilization: PmlStabilization::Mpml { ratio: 0.12 },
            ..Pml::default()
        };
        let cyl_rot = Pml {
            coord: PmlCoord::Cylindrical {
                center,
                axis: axis_rot,
                r_in: 0.5,
                r_thickness: 0.4,
                radial_on: true,
                axial: PmlDir::Plus,
                z_in: 0.6,
                z_thickness: 0.4,
            },
            ..cyl0.clone()
        };
        let jc0 = pml_jacobian(&cyl0, None, x, 3, false);
        let jc_rot = pml_jacobian(&cyl_rot, None, x_rot, 3, false);
        assert_j_times_jinv_is_identity(&jc0);
        let back_c = conj_by_rot(&q, &jc_rot.j);
        assert!(mat_diff_max(&back_c, &jc0.j) < 1e-13, "cylindrical Q^T J Q mismatch");
        assert!(c_abs(jc_rot.det - jc0.det) < 1e-13);

        // 3. Cartesian 90-degree rotation around z swaps x+ <-> y+
        let q90 = rotation([90.0, 0.0, 0.0]);
        let cart_x = Pml {
            x: PmlDir::Plus,
            bounds: Some(([1.0, -1.0, -1.0], [1.5, 1.0, 1.0])),
            profile: PmlProfile::KappaPolynomial { kappa_max: 2.5 },
            strength: 6.0,
            order: 2.0,
            ..Pml::default()
        };
        let cart_y = Pml {
            y: PmlDir::Plus,
            bounds: Some(([-1.0, 1.0, -1.0], [1.0, 1.5, 1.0])),
            profile: PmlProfile::KappaPolynomial { kappa_max: 2.5 },
            strength: 6.0,
            order: 2.0,
            ..Pml::default()
        };
        let pt_x = [1.3, 0.2, -0.4];
        let pt_y = rot_vec(&q90, pt_x);
        let jx = pml_jacobian(&cart_x, None, pt_x, 3, false);
        let jy = pml_jacobian(&cart_y, None, pt_y, 3, false);
        let back_xy = conj_by_rot(&q90, &jy.j);
        assert!(mat_diff_max(&back_xy, &jx.j) < 1e-14, "Cartesian 90-deg rotation mismatch");
    }

    #[test]
    fn axisymmetric_spherical_pml_is_exact_and_regular_on_r_equals_zero_axis() {
        let spec = Pml {
            coord: PmlCoord::Spherical { center: [0.0, 0.0, 0.5], r_in: 1.0, thickness: 0.4 },
            profile: PmlProfile::Rational { curvature: 0.9 },
            strength: 2.5,
            order: 2.0,
            ..Pml::default()
        };
        // Evaluate strictly AT r = 0 on the north pole z = 0.5 + 1.25 = 1.75 (t = 0.625)
        let j_pole = pml_jacobian(&spec, None, [0.0, 0.0, 1.75], 2, true);
        assert!(j_pole.det.re.is_finite() && j_pole.det.im.is_finite());
        assert!(j_pole.s_theta.re.is_finite() && j_pole.s_theta.im.is_finite());
        assert_j_times_jinv_is_identity(&j_pole);

        // At the north pole r = 0, R_hat points along +z (index 1 in 2D `(r, z, theta)`):
        // radial `r` (index 0) and azimuthal `theta` (index 2) are tangential -> both equal `s_Omega(R)`,
        // while axial `z` (index 1) is radial -> equals `s_R(R)`!
        assert!(c_abs(j_pole.j[0][0] - j_pole.s_theta) < 1e-14);
        assert!(c_abs(j_pole.j[2][2] - j_pole.s_theta) < 1e-14);
        assert!(c_abs(j_pole.j[0][1]) < 1e-15 && c_abs(j_pole.j[1][0]) < 1e-15);

        // Verify smooth limit as r -> 0^+ along the sphere R = 1.25
        let r_eps = 1e-8_f64;
        let z_eps = 0.5 + (1.25_f64 * 1.25 - r_eps * r_eps).sqrt();
        let j_near = pml_jacobian(&spec, None, [r_eps, 0.0, z_eps], 2, true);
        assert!(mat_diff_max(&j_near.j, &j_pole.j) < 1e-6);
        assert!(c_abs(j_near.s_theta - j_pole.s_theta) < 1e-14);
    }

    #[test]
    fn rational_antiderivative_matches_16pt_gauss_quadrature_to_machine_precision() {
        for &m in &[0.0, 1.0, 2.0, 3.0, 1.5, 2.5] {
            for &c_rat in &[0.0, 1e-7, 1e-3, 0.25, 0.7, 0.9, 0.95] {
                let spec = Pml {
                    x: PmlDir::Plus,
                    profile: PmlProfile::Rational { curvature: c_rat },
                    strength: 3.7,
                    order: m,
                    ..Pml::default()
                };
                for &t in &[0.0, 0.1, 0.37, 0.73, 1.0] {
                    let (_, int_analytic, _) = eval_profile_1d(&spec, t);
                    // Composite 8-panel x 16-point Gauss-Legendre reference quadrature of `s(xi)` over `[0, t]`
                    // with `xi = t * u^2` (`dxi = 2 * t * u du` on `u in [0, 1]`) so fractional `xi^m` is smooth at `0`.
                    let mut int_gl16 = c(0.0, 0.0);
                    const PANELS: usize = 32;
                    let h = 1.0 / (PANELS as f64);
                    for p in 0..PANELS {
                        let a = (p as f64) * h;
                        for q in 0..16 {
                            let u = a + h * GL16_NODES_01[q];
                            let xi = t * u * u;
                            let (s_xi, _, _) = eval_profile_1d(&spec, xi);
                            int_gl16 += s_xi * (2.0 * t * u * h * GL16_WEIGHTS_01[q]);
                        }
                    }
                    let err = c_abs(int_analytic - int_gl16);
                    let scale = c_abs(int_gl16).max(1.0);
                    assert!(
                        err <= 1e-12 * scale,
                        "m={m}, c_rat={c_rat}, t={t}: analytic={int_analytic:?} vs gl16={int_gl16:?} (err={err:.3e})"
                    );
                }
            }
        }
    }

    #[test]
    fn auto_scale_params_matches_closed_form_calibration() {
        let lam_max = 1.55e-6;
        let lam_min = 1.55e-6 / 3.48; // Si in Air/SiO2/Si group
        let res_open = auto_scale_params(lam_min, lam_max, Some(0.5 * lam_max), 2.0, 1e-6, false, false);
        assert!((res_open.d_auto - 0.5 * lam_max).abs() < 1e-20);
        // psi_max = 3 * ln(1e6) / (2 * pi) = 6.596419553536108
        let want_psi = 3.0 * 1e6_f64.ln() / (2.0 * std::f64::consts::PI);
        let want_bar_kappa = 1e6_f64.ln() / (2.0 * std::f64::consts::PI);
        let want_kappa = 1.0 + 3.0 * (want_bar_kappa - 1.0);
        assert!((res_open.psi_max - want_psi).abs() < 1e-12, "psi_max = {}", res_open.psi_max);
        assert!((res_open.kappa_max - want_kappa).abs() < 1e-12, "kappa_max = {}", res_open.kappa_max);
        assert_eq!(res_open.order, 2.0);
        assert_eq!(res_open.layer_grading, 2.0);
        assert!(res_open.layers_p2 >= 8 && res_open.layers_p1 > res_open.layers_p2);

        // Eigenfrequency mode caps kappa_max <= 2.0; homogeneous lam_min == lam_max gives 8 layers at P2
        let res_eig = auto_scale_params(lam_max, lam_max, None, 2.0, 1e-6, true, false);
        assert_eq!(res_eig.kappa_max, 2.0);
        assert_eq!(res_eig.layers_p2, 8);
        assert_eq!(res_eig.layers_p1, 16);

        // Thick layer d = 1.2 * lam_max -> natural evanescent decay beats R0 so kappa_max == 1.0
        let res_thick = auto_scale_params(lam_max, lam_max, Some(1.2 * lam_max), 2.0, 1e-6, false, false);
        assert_eq!(res_thick.kappa_max, 1.0);
    }

    #[test]
    fn ie_rational_and_exponential_jacobians_match_exact_mappings() {
        let (u0, d, d_pole) = (2.0, 0.5, 2.0);
        let spec = Pml {
            coord: PmlCoord::Spherical { center: [0.0, 0.0, 0.0], r_in: u0, thickness: d },
            ..Pml::default()
        };
        for &t in &[0.0, 0.25, 0.6, 0.9, 1.0] {
            let r = u0 + t * d;
            let (s_u, du) = ie_rational_stretch_1d(t, d, d_pole);
            let mapped_ratio = (d_pole + du) / d_pole;
            let exact_ratio = 1.0 / (1.0 - (t * d) / (d_pole + d));
            assert!((mapped_ratio - exact_ratio).abs() < 1e-15);

            let jac = ie_rational_jacobian(&spec, None, d_pole, [r, 0.0, 0.0], 3, false);
            assert!((jac.j[0][0].re - s_u).abs() < 1e-14 && jac.j[0][0].im == 0.0);
            let r_tilde = u0 + du;
            assert!((jac.s_theta.re - r_tilde / r).abs() < 1e-14);
            assert_j_times_jinv_is_identity(&jac);

            // Exponential IE with chi = 30 in [20, 50]
            let chi = 30.0_f64;
            let jac_exp = ie_exponential_jacobian(&spec, None, chi, [r, 0.0, 0.0], 3, false);
            assert!((jac_exp.j[0][0].re - chi.powf(t)).abs() < 1e-13);
            assert_j_times_jinv_is_identity(&jac_exp);
        }
    }

    #[test]
    fn pullback_maxwell_matches_2d_waveguide_and_axisym_formulas() {
        // 2D waveguide `(x, z, y)` Cartesian stretch
        let p2d = Pml {
            x: PmlDir::Plus,
            z: PmlDir::Minus,
            bounds: Some(([1.0, 0.0, -2.0], [1.5, 0.0, -1.5])),
            strength: 4.0,
            order: 2.0,
            ..Pml::default()
        };
        let jac = pml_jacobian(&p2d, None, [1.2, 0.0, -1.8], 2, false);
        let (sx, sz) = (jac.j[0][0], jac.j[1][1]);
        let eps = [
            [c(2.0, 0.0), c(0.3, 0.0), c(0.5, 0.0)],
            [c(0.3, 0.0), c(3.0, 0.0), c(0.7, 0.0)],
            [c(0.5, 0.0), c(0.7, 0.0), c(4.0, 0.0)],
        ];
        let mu_inv = [
            [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
            [c(0.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)],
            [c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)],
        ];
        let (et, mt) = pullback_maxwell(eps, mu_inv, &jac);
        assert!(c_abs(et[0][0] - (sz / sx) * eps[0][0]) < 1e-14);
        assert!(c_abs(et[1][1] - (sx / sz) * eps[1][1]) < 1e-14);
        assert!(c_abs(et[2][2] - (sx * sz) * eps[2][2]) < 1e-14);
        assert!(c_abs(et[0][1] - eps[0][1]) < 1e-14);
        assert!(c_abs(et[0][2] - sz * eps[0][2]) < 1e-14);
        assert!(c_abs(et[1][2] - sx * eps[1][2]) < 1e-14);
        assert!(c_abs(mt[0][0] - sx / sz) < 1e-14);
        assert!(c_abs(mt[1][1] - sz / sx) < 1e-14);
        assert!(c_abs(mt[2][2] - c(1.0, 0.0) / (sx * sz)) < 1e-14);
    }
}
