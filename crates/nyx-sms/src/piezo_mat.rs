// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Piezoelectric constitutive evaluation, 2D `Plane::Stress` electromechanical static condensation,
//! per-domain material resolution, and block power-of-2 + symmetric Ruiz equilibration.

use crate::problem::{Plane, Problem};
use crate::solid::thermal_strain;
use femcore::adapt::metric::{sym_eig, symmetrize};
use femcore::linalg::csr::{CCsr, Csr, RCsr, C};
use femcore::linalg::equilibrate::{
    ruiz_symmetric_terms, scale_symmetric, scale_vec, unscale_vec, RuizTerm, DEFAULT_RUIZ_PASSES,
};
use femcore::linalg::scalar::Scalar;
use femcore::material::{self, M3, M6};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::mesh::SimplexMesh;
use femcore::props::Materials;
use femcore::section::{plane_piezo_e, plane_stiffness, plane_tensor};
use std::collections::BTreeMap;

/// Vacuum permittivity `epsilon_0` (F/m), shared from `femcore::material::coupled`.
pub const EPSILON_0: f64 = femcore::material::coupled::EPSILON_0;

/// In-plane Voigt indices `[xx, zz, xz]` inside 2D cross-section order `(xx, zz, yy, xz, zy, xy)`.
pub const INPLANE_2D: [usize; 3] = [0, 1, 3];

/// Resolved global/cross-section constitutive properties on a single cell/domain.
#[derive(Clone, Copy, Debug)]
pub struct PiezoCellMat {
    /// Short-circuit elastic stiffness `C^E` (Pa), solver Voigt (3D) or cross-section Voigt (2D).
    pub c6: M6,
    /// Piezoelectric stress-charge coupling `e` (`C/m^2`), solver Voigt (3D) or `plane_piezo_e` (2D).
    pub e: [[f64; 6]; 3],
    /// Real part of clamped permittivity `eps'^S` (F/m), global (3D) or `plane_tensor` (2D).
    pub eps_re: M3,
    /// Imaginary loss part of permittivity `eps''^S` (F/m, `eps = eps_re - i * eps_im`).
    pub eps_im: M3,
    /// Kelvin-Voigt structural viscosity `eta` (Pa*s), solver Voigt (3D) or `plane_stiffness` (2D).
    pub eta: M6,
    /// Mass density `rho` (`kg/m^3`), zero on mechanically inactive domains.
    pub rho: f64,
    /// True if domain participates in mechanical stiffness/mass/coupling (`!inactive`).
    pub mech_active: bool,
    /// True if domain participates in electrical permittivity/coupling (`!conductor`).
    pub elec_active: bool,
}

/// Exact 2D `Plane::Stress` (`T_yy = 0`) electromechanical static condensation of `(C^E, e, eps^S)`.
#[derive(Clone, Copy, Debug)]
pub struct PlaneStressPiezoCond {
    /// Condensed 3x3 row-major in-plane stiffness `C_ps^E` in `(xx, zz, xz)`.
    pub d_uu: [f64; 9],
    /// Condensed 2x3 row-major in-plane coupling `e_ps` (`[x, z] x [xx, zz, xz]`).
    pub e_pu: [f64; 6],
    /// Condensed 2x2 row-major semi-free in-plane real permittivity `eps_ps' = eps' + e_y e_y^T / C22`.
    pub eps_ps_re: [f64; 4],
    /// In-plane 2x2 row-major imaginary loss permittivity `eps_in''`.
    pub eps_ps_im: [f64; 4],
}

/// Performs exact 2D `Plane::Stress` (`T_yy = 0`) electromechanical static condensation.
///
/// Guards against zero out-of-plane stiffness (`c[2][2] <= 1e-30`) so inactive/air cells
/// never evaluate `0.0 / 0.0 = NaN`.
pub fn condense_plane_stress(c: &M6, e: &[[f64; 6]; 3], eps_re: &M3, eps_im: &M3) -> PlaneStressPiezoCond {
    let eps_ps_im = [eps_im[0][0], eps_im[0][1], eps_im[1][0], eps_im[1][1]];
    if c[2][2].abs() <= 1e-30 {
        return PlaneStressPiezoCond {
            d_uu: [0.0; 9],
            e_pu: [0.0; 6],
            eps_ps_re: [eps_re[0][0], eps_re[0][1], eps_re[1][0], eps_re[1][1]],
            eps_ps_im,
        };
    }
    let inv_c22 = 1.0 / c[2][2];
    let mut d_uu = [0.0; 9];
    for (a, &ia) in INPLANE_2D.iter().enumerate() {
        for (b, &ib) in INPLANE_2D.iter().enumerate() {
            d_uu[a * 3 + b] = c[ia][ib] - c[ia][2] * c[2][ib] * inv_c22;
        }
    }
    let mut e_pu = [0.0; 6];
    for p in 0..2 {
        for (b, &ib) in INPLANE_2D.iter().enumerate() {
            e_pu[p * 3 + b] = e[p][ib] - e[p][2] * c[2][ib] * inv_c22;
        }
    }
    let mut eps_ps_re = [0.0; 4];
    for p in 0..2 {
        for q in 0..2 {
            eps_ps_re[p * 2 + q] = eps_re[p][q] + e[p][2] * e[q][2] * inv_c22;
        }
    }
    PlaneStressPiezoCond { d_uu, e_pu, eps_ps_re, eps_ps_im }
}

/// Evaluates in-plane / 3D constitutive blocks `(d_uu, e_pu, eps_re, eps_im)` used for cell assembly.
pub(crate) fn cell_constitutive_blocks(mat: &PiezoCellMat, dim: usize, plane: Plane) -> PlaneStressPiezoCond {
    debug_assert_eq!(dim, 2, "cell_constitutive_blocks 2D helper called with dim={dim}");
    match plane {
        Plane::Stress => condense_plane_stress(&mat.c6, &mat.e, &mat.eps_re, &mat.eps_im),
        Plane::Strain => {
            let mut d_uu = [0.0; 9];
            for (a, &ia) in INPLANE_2D.iter().enumerate() {
                for (b, &ib) in INPLANE_2D.iter().enumerate() {
                    d_uu[a * 3 + b] = mat.c6[ia][ib];
                }
            }
            let mut e_pu = [0.0; 6];
            for p in 0..2 {
                for (b, &ib) in INPLANE_2D.iter().enumerate() {
                    e_pu[p * 3 + b] = mat.e[p][ib];
                }
            }
            PlaneStressPiezoCond {
                d_uu,
                e_pu,
                eps_ps_re: [mat.eps_re[0][0], mat.eps_re[0][1], mat.eps_re[1][0], mat.eps_re[1][1]],
                eps_ps_im: [mat.eps_im[0][0], mat.eps_im[0][1], mat.eps_im[1][0], mat.eps_im[1][1]],
            }
        }
    }
}

/// Pointwise stress (`[xx, yy, zz, xy, yz, xz]` in 3D or `[xx, zz, yy, xz, 0, 0]` in 2D),
/// electric field `E = -grad phi` (`[Ex, Ey, Ez]` or `[Ex, Ez, 0]`), and electric displacement `D`.
#[derive(Clone, Copy, Debug, Default)]
pub struct PiezoPointEval {
    pub stress: [f64; 6],
    /// Elastic strain in 3D wafer solver Voigt order `[xx, yy, zz, xy, yz, xz]` (engineering shear).
    pub strain_wafer: [f64; 6],
    pub e_field: [f64; 3],
    pub d_field: [f64; 3],
    pub elastic_energy_density: f64,
    pub dielectric_energy_density: f64,
}

/// Evaluates pointwise coupled piezoelectric stress, electric field `E`, electric displacement `D`,
/// and stored elastic & dielectric energy densities `(0.5 S_el^T C^E S_el, 0.5 E^T eps^S E)`.
pub fn eval_piezo_point(
    mat: &PiezoCellMat,
    dim: usize,
    plane: Plane,
    grad_u: &[[f64; 3]; 3],
    grad_phi: &[f64; 3],
    thermal: Option<(&M3, f64)>,
) -> PiezoPointEval {
    let e_field = [-grad_phi[0], -grad_phi[1], -grad_phi[2]];
    let eth = thermal.map_or([0.0; 6], |(alpha, dt)| thermal_strain(alpha, dt));

    if !mat.mech_active && !mat.elec_active {
        return PiezoPointEval::default();
    }

    if dim == 3 {
        let s_el = [
            grad_u[0][0] - eth[0],
            grad_u[1][1] - eth[1],
            grad_u[2][2] - eth[2],
            grad_u[0][1] + grad_u[1][0] - eth[3],
            grad_u[1][2] + grad_u[2][1] - eth[4],
            grad_u[0][2] + grad_u[2][0] - eth[5],
        ];
        let mut stress = [0.0; 6];
        let mut u_e = 0.0;
        if mat.mech_active {
            for i in 0..6 {
                let mut c_s = 0.0;
                for j in 0..6 {
                    c_s += mat.c6[i][j] * s_el[j];
                }
                u_e += 0.5 * s_el[i] * c_s;
                let mut e_grad = 0.0;
                if mat.elec_active {
                    for p in 0..3 {
                        e_grad += mat.e[p][i] * grad_phi[p];
                    }
                }
                stress[i] = c_s + e_grad;
            }
        }
        let mut d_field = [0.0; 3];
        let mut u_d = 0.0;
        if mat.elec_active {
            for p in 0..3 {
                let mut es = 0.0;
                if mat.mech_active {
                    for j in 0..6 {
                        es += mat.e[p][j] * (s_el[j] + eth[j]);
                    }
                }
                let mut eps_e = 0.0;
                for q in 0..3 {
                    eps_e += mat.eps_re[p][q] * e_field[q];
                }
                d_field[p] = es + eps_e;
                u_d += 0.5 * e_field[p] * eps_e;
            }
        }
        return PiezoPointEval {
            stress,
            strain_wafer: if mat.mech_active { s_el } else { [0.0; 6] },
            e_field: if mat.elec_active { e_field } else { [0.0; 3] },
            d_field,
            elastic_energy_density: u_e,
            dielectric_energy_density: u_d,
        };
    }

    // 2D cross-section `(x, z)` order: Voigt `[xx=0, zz=1, yy=2, xz=3]`
    let s_xx = grad_u[0][0] - eth[0];
    let s_zz = grad_u[1][1] - eth[1];
    let g_xz = grad_u[0][1] + grad_u[1][0] - eth[3];
    let s_yy_strain = -eth[2];
    let s_yy_eff = match plane {
        Plane::Strain => s_yy_strain,
        Plane::Stress => {
            if !mat.mech_active || mat.c6[2][2].abs() <= 1e-30 {
                0.0
            } else {
                let e_dot_g = if mat.elec_active { mat.e[0][2] * grad_phi[0] + mat.e[1][2] * grad_phi[1] } else { 0.0 };
                -(mat.c6[2][0] * s_xx + mat.c6[2][1] * s_zz + mat.c6[2][3] * g_xz + e_dot_g) / mat.c6[2][2]
            }
        }
    };
    let s4 = [s_xx, s_zz, s_yy_eff, g_xz];
    let s_tot4 = [s_xx + eth[0], s_zz + eth[1], s_yy_eff + eth[2], g_xz + eth[3]];
    let mut stress = [0.0; 6];
    let mut u_e = 0.0;
    if mat.mech_active {
        for &i in &[0usize, 1, 2, 3] {
            let mut c_s = 0.0;
            for (k, &j) in [0usize, 1, 2, 3].iter().enumerate() {
                c_s += mat.c6[i][j] * s4[k];
            }
            let s_i = match i {
                0 => s_xx,
                1 => s_zz,
                2 => s_yy_eff,
                _ => g_xz,
            };
            u_e += 0.5 * s_i * c_s;
            let e_grad = if mat.elec_active { mat.e[0][i] * grad_phi[0] + mat.e[1][i] * grad_phi[1] } else { 0.0 };
            stress[i] = if i == 2 && plane == Plane::Stress { 0.0 } else { c_s + e_grad };
        }
    }
    let mut d_field = [0.0; 3];
    let mut u_d = 0.0;
    if mat.elec_active {
        for p in 0..2 {
            let mut es = 0.0;
            if mat.mech_active {
                for (k, &j) in [0usize, 1, 2, 3].iter().enumerate() {
                    es += mat.e[p][j] * s_tot4[k];
                }
            }
            let eps_e = mat.eps_re[p][0] * e_field[0] + mat.eps_re[p][1] * e_field[1];
            d_field[p] = es + eps_e;
            u_d += 0.5 * e_field[p] * eps_e;
        }
    }
    PiezoPointEval {
        stress,
        strain_wafer: if mat.mech_active { [s_xx, s_yy_eff, s_zz, 0.0, 0.0, g_xz] } else { [0.0; 6] },
        e_field: if mat.elec_active { [e_field[0], e_field[1], 0.0] } else { [0.0; 3] },
        d_field,
        elastic_energy_density: u_e,
        dielectric_energy_density: u_d,
    }
}

fn check_permittivity_spd(dom: i32, eps_re: &M3, eps_im: &M3) -> Result<(), String> {
    let s_re = symmetrize(eps_re);
    let (ev, _) = sym_eig::<3>(&s_re);
    if ev.iter().any(|&v| !v.is_finite() || v <= 0.0) {
        return Err(format!(
            "domain {dom}: clamped permittivity tensor is not positive definite (eigenvalues = {ev:?})"
        ));
    }
    let s_im = symmetrize(eps_im);
    let (ev_im, _) = sym_eig::<3>(&s_im);
    let tol = 1e-12 * ev[0].max(1e-30);
    if ev_im.iter().any(|&v| !v.is_finite() || v < -tol) {
        return Err(format!(
            "domain {dom}: imaginary permittivity loss tensor must be positive semi-definite (eigenvalues = {ev_im:?})"
        ));
    }
    Ok(())
}

/// Resolves per-domain piezoelectric material tables over unique domain labels.
pub fn resolve_domain_mats(
    domain_labels: &[i32],
    materials: &Materials,
    inactive: &[i32],
    conductor_domains: &[i32],
    need_rho: bool,
    dim: usize,
) -> Result<BTreeMap<i32, PiezoCellMat>, String> {
    resolve_domain_mats_with_u_domains(domain_labels, materials, inactive, &[], conductor_domains, need_rho, dim)
}

/// Resolves per-domain piezoelectric material tables over unique domain labels, supporting an
/// optional explicit solid submesh `u_domains` as well as automatic dielectric-only vacuum/air
/// domains (`eps_r` without elastic stiffness `E`/`stiffness`).
pub fn resolve_domain_mats_with_u_domains(
    domain_labels: &[i32],
    materials: &Materials,
    inactive: &[i32],
    u_domains: &[i32],
    conductor_domains: &[i32],
    need_rho: bool,
    dim: usize,
) -> Result<BTreeMap<i32, PiezoCellMat>, String> {
    let any_solid_has_elastic = domain_labels.iter().any(|&dom| {
        !inactive.contains(&dom)
            && (u_domains.is_empty() || u_domains.contains(&dom))
            && (materials.aniso.get(&dom).and_then(|a| a.stiffness.as_ref()).is_some()
                || materials.get(dom, "E").is_some_and(|e| e > 0.0))
    });
    let mut out = BTreeMap::new();
    for &dom in domain_labels {
        let an = materials.aniso.get(&dom);
        let has_elastic =
            an.and_then(|a| a.stiffness.as_ref()).is_some() || materials.get(dom, "E").is_some_and(|e| e > 0.0);
        let in_u_domains = u_domains.is_empty() || u_domains.contains(&dom);
        let mech_active = !inactive.contains(&dom) && in_u_domains && (has_elastic || !any_solid_has_elastic);
        let elec_active = !conductor_domains.contains(&dom);

        let c6 = if !mech_active {
            [[0.0; 6]; 6]
        } else {
            let cw = if let Some(st) = an.and_then(|a| a.stiffness.as_ref()) {
                material::stiffness_global(st, an.unwrap().euler_deg)
            } else {
                let e = materials.get(dom, "E").ok_or_else(|| format!("domain {dom} has no value for 'E'"))?;
                let nu = materials.get(dom, "nu").ok_or_else(|| format!("domain {dom} has no value for 'nu'"))?;
                if nu >= 0.5 || nu <= -1.0 {
                    return Err("Poisson's ratio must be between -1 and 0.5".into());
                }
                if e <= 0.0 {
                    return Err("Young's modulus must be positive".into());
                }
                material::isotropic_stiffness(e, nu)
            };
            if !material::is_positive_definite(&cw) {
                return Err(format!("the stiffness of domain {dom} is not positive definite (check the constants)"));
            }
            plane_stiffness(&cw, dim)
        };

        let eta = if !mech_active {
            [[0.0; 6]; 6]
        } else if let Some(v) = an.and_then(|a| a.viscosity.as_ref()) {
            let ew = material::rotate_m6_global(v, an.unwrap().euler_deg);
            plane_stiffness(&ew, dim)
        } else if let Some(eta_s) = materials.get(dom, "eta") {
            if eta_s < 0.0 {
                return Err(format!("viscosity 'eta' of domain {dom} must be non-negative (found {eta_s})"));
            }
            if eta_s > 0.0 {
                plane_stiffness(&material::isotropic_stiffness(eta_s, 0.0), dim)
            } else {
                [[0.0; 6]; 6]
            }
        } else {
            [[0.0; 6]; 6]
        };

        let rho = if !mech_active || !need_rho {
            0.0
        } else {
            let r = materials.get(dom, "rho").ok_or_else(|| format!("domain {dom} has no value for 'rho'"))?;
            if r < 0.0 {
                return Err(format!("density 'rho' of domain {dom} must be non-negative (found {r})"));
            }
            r
        };

        let e = if mech_active && elec_active {
            if let Some(e_tb) = an.and_then(|a| a.piezo_e.as_ref()) {
                let eg = material::piezo_e_global(e_tb, an.unwrap().euler_deg);
                plane_piezo_e(&eg, dim)
            } else {
                [[0.0; 6]; 3]
            }
        } else {
            [[0.0; 6]; 3]
        };

        let (eps_re, eps_im) = if !elec_active {
            ([[0.0; 3]; 3], [[0.0; 3]; 3])
        } else {
            let re_w = match an.and_then(|a| a.permittivity.as_ref()) {
                Some(t) => material::tensor_global(t, an.unwrap().euler_deg),
                None => {
                    let s = materials.get(dom, "eps_r").unwrap_or(1.0);
                    material::iso_tensor(s)
                }
            };
            let im_w = match an.and_then(|a| a.permittivity_im.as_ref()) {
                Some(t) => material::tensor_global(t, an.unwrap().euler_deg),
                None => {
                    if let Some(s_im) = materials.get(dom, "eps_r_im") {
                        material::iso_tensor(s_im)
                    } else {
                        let td = materials.get(dom, "tan_delta").unwrap_or(0.0);
                        if td > 0.0 {
                            [
                                [re_w[0][0] * td, re_w[0][1] * td, re_w[0][2] * td],
                                [re_w[1][0] * td, re_w[1][1] * td, re_w[1][2] * td],
                                [re_w[2][0] * td, re_w[2][1] * td, re_w[2][2] * td],
                            ]
                        } else {
                            [[0.0; 3]; 3]
                        }
                    }
                }
            };
            check_permittivity_spd(dom, &re_w, &im_w)?;
            let mut er = plane_tensor(&re_w, dim);
            let mut ei = plane_tensor(&im_w, dim);
            for i in 0..3 {
                for j in 0..3 {
                    er[i][j] *= EPSILON_0;
                    ei[i][j] *= EPSILON_0;
                }
            }
            (er, ei)
        };

        out.insert(dom, PiezoCellMat { c6, e, eps_re, eps_im, eta, rho, mech_active, elec_active });
    }
    Ok(out)
}

/// Resolves per-cell [`PiezoCellMat`] entries for a [`SimplexMesh`] and validates domain data.
pub fn resolve_piezo_mats_simplex(
    mesh: &SimplexMesh,
    spec: &Problem,
    need_rho: bool,
) -> Result<(BTreeMap<i32, PiezoCellMat>, Vec<PiezoCellMat>), String> {
    let doms = mesh.domain_labels();
    let conds = spec.conductor_domains();
    let by_dom = resolve_domain_mats_with_u_domains(
        &doms,
        &spec.materials,
        &spec.inactive,
        &spec.u_domains,
        &conds,
        need_rho,
        mesh.dim,
    )?;
    if !by_dom.values().any(|m| m.mech_active) {
        return Err(
            "coupled piezoelectric problem has no mechanically active solid domain (check `inactive` and `u_domains`)"
                .into(),
        );
    }
    let cells = mesh.cell_labels.iter().map(|l| by_dom[l]).collect();
    Ok((by_dom, cells))
}

/// Resolves per-cell [`PiezoCellMat`] entries for a [`MixedDiscretisation`] and validates domain data.
pub fn resolve_piezo_mats_mixed(
    mixed: &MixedDiscretisation,
    spec: &Problem,
    need_rho: bool,
) -> Result<(BTreeMap<i32, PiezoCellMat>, Vec<PiezoCellMat>), String> {
    let mut doms = mixed.cell_labels.clone();
    doms.sort_unstable();
    doms.dedup();
    let conds = spec.conductor_domains();
    let by_dom = resolve_domain_mats_with_u_domains(
        &doms,
        &spec.materials,
        &spec.inactive,
        &spec.u_domains,
        &conds,
        need_rho,
        mixed.dim,
    )?;
    if !by_dom.values().any(|m| m.mech_active) {
        return Err(
            "coupled piezoelectric problem has no mechanically active solid domain (check `inactive` and `u_domains`)"
                .into(),
        );
    }
    let cells = mixed.cell_labels.iter().map(|l| by_dom[l]).collect();
    Ok((by_dom, cells))
}

/// Computes the exact power-of-2 block pre-scaling `d_blk` placing mechanical `[0..n_u)`
/// and electrical `[n_u..n)` blocks on comparable magnitudes before symmetric Ruiz equilibration.
pub fn power_of_two_block_scale<T: Scalar>(a: &Csr<T>, n_u: usize) -> Vec<f64> {
    let n = a.nrows();
    if n == 0 || n_u == 0 || n_u >= n {
        return vec![1.0; n];
    }
    let diag = a.diagonal();
    let mut u_max = 0.0_f64;
    for v in &diag[..n_u] {
        u_max = u_max.max(v.modulus());
    }
    if u_max <= 1e-300 {
        for r in 0..n_u {
            for v in a.row(r).1 {
                u_max = u_max.max(v.modulus());
            }
        }
    }
    let mut p_max = 0.0_f64;
    for v in &diag[n_u..n] {
        p_max = p_max.max(v.modulus());
    }
    if p_max <= 1e-300 {
        for r in n_u..n {
            for v in a.row(r).1 {
                p_max = p_max.max(v.modulus());
            }
        }
    }
    let su = if u_max > 1e-300 { (-(0.5 * u_max.log2()).round()).clamp(-60.0, 60.0).exp2() } else { 1.0 };
    let sp = if p_max > 1e-300 { (-(0.5 * p_max.log2()).round()).clamp(-60.0, 60.0).exp2() } else { 1.0 };
    let mut d = vec![su; n];
    d[n_u..n].fill(sp);
    d
}

/// Applies power-of-2 block scaling followed by 8-pass symmetric Ruiz equilibration to a real system.
pub fn equilibrate_real(a: &RCsr, f: &[f64], n_u: usize) -> (RCsr, Vec<f64>, Vec<f64>) {
    let n = a.nrows();
    let d_blk = power_of_two_block_scale(a, n_u);
    let a_blk = scale_symmetric(a, &d_blk);
    let d_ruiz = ruiz_symmetric_terms(n, &[RuizTerm::Real(&a_blk, 1.0)], DEFAULT_RUIZ_PASSES);
    let d_tot: Vec<f64> = d_blk.iter().zip(&d_ruiz).map(|(&b, &r)| (b * r).clamp(1e-24, 1e24)).collect();
    let a_hat = scale_symmetric(a, &d_tot);
    let f_hat = scale_vec(f, &d_tot);
    (a_hat, f_hat, d_tot)
}

/// Computes power-of-2 block scaling + 8-pass symmetric Ruiz diagonal `d_tot` over weighted complex/real terms.
pub fn compute_piezo_complex_scaling(n: usize, n_u: usize, ref_mat: &CCsr, terms: &[(&CCsr, f64)]) -> Vec<f64> {
    let d_blk = power_of_two_block_scale(ref_mat, n_u);
    let scaled_mats: Vec<(CCsr, f64)> = terms.iter().map(|&(m, w)| (scale_symmetric(m, &d_blk), w)).collect();
    let ruiz_terms: Vec<RuizTerm<'_>> = scaled_mats.iter().map(|(m, w)| RuizTerm::Complex(m, *w)).collect();
    let d_ruiz = ruiz_symmetric_terms(n, &ruiz_terms, DEFAULT_RUIZ_PASSES);
    d_blk.iter().zip(&d_ruiz).map(|(&b, &r)| (b * r).clamp(1e-24, 1e24)).collect()
}

/// Computes power-of-2 block scaling + 8-pass symmetric Ruiz diagonal `d_tot` for real affine terms `(A, w)`.
pub fn compute_piezo_real_scaling(n: usize, n_u: usize, ref_mat: &RCsr, terms: &[(&RCsr, f64)]) -> Vec<f64> {
    let d_blk = power_of_two_block_scale(ref_mat, n_u);
    let scaled_mats: Vec<(RCsr, f64)> = terms.iter().map(|&(m, w)| (scale_symmetric(m, &d_blk), w)).collect();
    let ruiz_terms: Vec<RuizTerm<'_>> = scaled_mats.iter().map(|(m, w)| RuizTerm::Real(m, *w)).collect();
    let d_ruiz = ruiz_symmetric_terms(n, &ruiz_terms, DEFAULT_RUIZ_PASSES);
    d_blk.iter().zip(&d_ruiz).map(|(&b, &r)| (b * r).clamp(1e-24, 1e24)).collect()
}

/// Restores physical solution `x = D x_hat` after equilibrated solve.
#[inline]
pub fn unscale_solution<T: Scalar>(x_hat: &[T], d_tot: &[f64]) -> Vec<T> {
    scale_vec(x_hat, d_tot)
}

/// Converts a reduced solution `x_phys = D x_hat` via symmetric factor `u = d_tot * x_hat` or `f_hat = d_tot * f`.
#[inline]
pub fn scale_rhs<T: Scalar>(f: &[T], d_tot: &[f64]) -> Vec<T> {
    scale_vec(f, d_tot)
}

/// Unscales inverse transform `v / d_tot` (re-export for symmetry checks).
#[inline]
pub fn div_scale<T: Scalar>(v: &[T], d_tot: &[f64]) -> Vec<T> {
    unscale_vec(v, d_tot)
}

/// Suppress unused import warnings if `C` is only referenced under complex tests.
#[allow(dead_code)]
fn _touch_c(_: C) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plane_stress_condensation_matches_zero_out_of_plane_stress_identity() {
        // Wurtzite AlN in cross-section `(xx, zz, yy, xz, zy, xy)`
        let cw = material::stiffness_global(
            &material::Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: 389e9, c44: 125e9 },
            [0.0, 0.0, 0.0],
        );
        let c2 = plane_stiffness(&cw, 2);
        let e_tb =
            [[0.0, 0.0, 0.0, 0.0, -0.48, 0.0], [0.0, 0.0, 0.0, -0.48, 0.0, 0.0], [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0]];
        let e2 = plane_piezo_e(&material::piezo_e_global(&e_tb, [0.0, 0.0, 0.0]), 2);
        let eps_r = material::iso_tensor(9.5);
        let mut eps = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                eps[i][j] = EPSILON_0 * eps_r[i][j];
            }
        }
        let cond = condense_plane_stress(&c2, &e2, &eps, &[[0.0; 3]; 3]);
        // e31_ps = e31 - e32 * c20 / c22 = -0.58 - (-0.58) * 149 / 410
        let expected_e31 = -0.58 - (-0.58) * (149e9 / 410e9);
        assert!((cond.e_pu[3] - expected_e31).abs() < 1e-12);
        // eps33_ps = eps33 + e32^2 / c22
        let expected_eps33 = EPSILON_0 * 9.5 + (-0.58_f64).powi(2) / 410e9;
        assert!((cond.eps_ps_re[3] - expected_eps33).abs() / expected_eps33 < 1e-12);

        // Zero stiffness must not produce NaN
        let zero_c = condense_plane_stress(&[[0.0; 6]; 6], &[[0.0; 6]; 3], &eps, &[[0.0; 3]; 3]);
        assert!(zero_c.d_uu.iter().all(|v| v.is_finite() && *v == 0.0));
        assert!(zero_c.e_pu.iter().all(|v| v.is_finite() && *v == 0.0));
    }
}
