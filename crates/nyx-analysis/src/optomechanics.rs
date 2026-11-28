// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cross-mesh cavity optomechanics and stimulated Brillouin scattering (SBS) coupling analysis (`"optomechanics"`, ADR 141).
//!
//! Operates purely on [`nyx_physics::SolveResult`] and [`nyx_physics::Mesh`] pairs (optical and mechanical)
//! using [`femcore::transfer::TransferIndex`] and [`femcore::transfer::BoundaryIndex`] without depending
//! on `nyx-ems` or `nyx-sms`.
//!
//! Computes:
//! - Modal effective mass $m_{\text{eff}} = \int_{\Omega_{\text{mech}}} \rho(\mathbf{x}) \|\tilde{\mathbf{u}}(\mathbf{x})\|^2 \,\mathrm{d}V$
//!   normalized to peak displacement $\max_{\mathbf{x} \in \Omega_{\text{mech}}} \|\tilde{\mathbf{u}}(\mathbf{x})\|_2 = 1$,
//! - Zero-point fluctuation amplitude $x_{\text{zpf}} = \sqrt{\hbar / (2 m_{\text{eff}} \Omega_m)}$,
//! - Cross-mesh photoelastic (electrostrictive) frequency pull $G_{\text{PE}}$ and single-photon coupling rate $g_{\text{PE}} = G_{\text{PE}} x_{\text{zpf}}$,
//!   integrated over optical dielectric core cells by sampling the mechanical strain tensor $\tilde{\mathbf{S}}(\mathbf{x}_q)$ strictly inside matching canonical `domain_id`s,
//! - Cross-mesh moving-boundary (radiation-pressure) frequency pull $G_{\text{MB}}$ and single-photon coupling rate $g_{\text{MB}} = G_{\text{MB}} x_{\text{zpf}}$,
//!   integrated over optical material interfaces using Johnson's perturbation theorem
//!   $(\varepsilon_1 - \varepsilon_2)|\mathbf{E}_\parallel|^2 - (\varepsilon_1^{-1} - \varepsilon_2^{-1})|D_\perp|^2$,
//! - Total optomechanical coupling $G_{\text{OM}} = G_{\text{PE}} + G_{\text{MB}}$ and $g_0 = g_{\text{PE}} + g_{\text{MB}}$ (in $\text{rad/s}$ and $\text{Hz}$),
//! - Backward/forward stimulated Brillouin scattering (SBS) gain $G_B$ ($\text{W}^{-1}\text{m}^{-1}$) for waveguide modes.

use crate::Analysis;
use femcore::fem::{geometry, quad_rule, shape, RefElement};
use femcore::linalg::csr::{c, C};
use femcore::material::coupled::delta_eps_eo_pe;
use femcore::material::{iso_tensor, M3, M6};
use femcore::optical::{C0, EPS0};
use femcore::transfer::{same_simplex_mesh, BoundaryIndex, SourceField, SourceMesh, TransferIndex, TransferOptions};
use nyx_physics::{discretise, Discretisation, Error, FieldFrame, Mesh, Solution, SolveResult};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// Reduced Planck constant $\hbar$ ($\text{J}\cdot\text{s}$, CODATA 2018 exact).
pub const HBAR: f64 = 1.054_571_817e-34;

/// Per-domain optical, photoelastic, and mass-density properties in the global wafer frame.
#[derive(Clone, Debug, PartialEq)]
pub struct OptomechanicsDomainMat {
    /// Relative permittivity tensor $\boldsymbol{\varepsilon}_r$ in global wafer axes `(x, y, z)`.
    pub eps_r: M3,
    /// Optional photoelastic Pockels tensor $p_{IJ}$ in global wafer axes (**textbook Voigt order** `(xx, yy, zz, yz, xz, xy)`).
    pub photoelastic: Option<M6>,
    /// Mass density $\rho$ ($\text{kg/m}^3$, `0.0` for vacuum/air cladding).
    pub rho: f64,
}

impl Default for OptomechanicsDomainMat {
    fn default() -> Self {
        Self { eps_r: iso_tensor(1.0), photoelastic: None, rho: 0.0 }
    }
}

/// Options and material/domain tables for [`optomechanics_with_meshes`].
#[derive(Clone, Debug)]
pub struct OptomechanicsOptions {
    /// 0-based optical pump mode index in `opt_result.solutions`.
    pub optical_mode_index: usize,
    /// Optional 0-based optical Stokes/scattered mode index (defaults to `optical_mode_index` for self-coupling).
    pub stokes_mode_index: Option<usize>,
    /// 0-based mechanical mode index in `mech_result.solutions`.
    pub mechanical_mode_index: usize,
    /// Out-of-plane integration length $L_y$ ($\text{m}$) for 2D cross-section meshes (ignored in 3D).
    pub length_m: f64,
    /// Optional explicit optical frequency $f_{\text{opt}}$ ($\text{Hz}$); inferred from `opt_result` when `None`.
    pub optical_freq_hz: Option<f64>,
    /// Optional mechanical quality factor $Q_m$ for SBS gain $G_B$; inferred from `mech_result` or defaults to `1000.0` when `None`.
    pub q_mech: Option<f64>,
    /// Optional optical group velocity $v_g$ ($\text{m/s}$) for 2D waveguide SBS gain; inferred from `opt_result` when `None`.
    pub v_g_opt: Option<f64>,
    /// Canonical `domain_id -> region_name` table for the optical mesh.
    pub opt_domain_names: BTreeMap<i32, String>,
    /// Canonical `boundary_id -> boundary_name` table for the optical mesh.
    pub opt_boundary_names: BTreeMap<i32, String>,
    /// Canonical `domain_id -> region_name` table for the mechanical mesh.
    pub mech_domain_names: BTreeMap<i32, String>,
    /// Canonical `boundary_id -> boundary_name` table for the mechanical mesh.
    pub mech_boundary_names: BTreeMap<i32, String>,
    /// Material properties keyed by canonical region name.
    pub domain_materials: BTreeMap<String, OptomechanicsDomainMat>,
    /// Fallback material properties keyed by numeric domain label when region names are omitted.
    pub domain_materials_by_id: BTreeMap<i32, OptomechanicsDomainMat>,
}

impl Default for OptomechanicsOptions {
    fn default() -> Self {
        Self {
            optical_mode_index: 0,
            stokes_mode_index: None,
            mechanical_mode_index: 0,
            length_m: 1.0,
            optical_freq_hz: None,
            q_mech: None,
            v_g_opt: None,
            opt_domain_names: BTreeMap::new(),
            opt_boundary_names: BTreeMap::new(),
            mech_domain_names: BTreeMap::new(),
            mech_boundary_names: BTreeMap::new(),
            domain_materials: BTreeMap::new(),
            domain_materials_by_id: BTreeMap::new(),
        }
    }
}

fn parse_m3_value(v: &Value, name: &str) -> Result<M3, Error> {
    if let Some(s) = v.as_f64() {
        if !s.is_finite() || s <= 0.0 {
            return Err(Error::invalid(format!("{name}: scalar permittivity must be positive and finite")));
        }
        return Ok(iso_tensor(s));
    }
    if let Some(arr) = v.as_array() {
        if arr.len() == 3 && arr.iter().all(Value::is_number) {
            let d0 = arr[0].as_f64().unwrap_or(1.0);
            let d1 = arr[1].as_f64().unwrap_or(1.0);
            let d2 = arr[2].as_f64().unwrap_or(1.0);
            if !d0.is_finite() || !d1.is_finite() || !d2.is_finite() {
                return Err(Error::invalid(format!("{name}: diagonal tensor components must be finite")));
            }
            return Ok([[d0, 0.0, 0.0], [0.0, d1, 0.0], [0.0, 0.0, d2]]);
        }
        if arr.len() == 3 && arr.iter().all(Value::is_array) {
            let mut m = [[0.0; 3]; 3];
            for (r, row) in arr.iter().enumerate() {
                let Some(ra) = row.as_array().filter(|a| a.len() == 3) else {
                    return Err(Error::invalid(format!("{name}: expected 3x3 matrix")));
                };
                for c_idx in 0..3 {
                    let x = ra[c_idx]
                        .as_f64()
                        .filter(|v| v.is_finite())
                        .ok_or_else(|| Error::invalid(format!("{name}[{r}][{c_idx}] must be finite")))?;
                    m[r][c_idx] = x;
                }
            }
            return Ok(m);
        }
    }
    Err(Error::invalid(format!("{name}: expected positive number, [eps_x, eps_y, eps_z], or 3x3 matrix")))
}

fn parse_m6_value(v: &Value, name: &str) -> Result<M6, Error> {
    if let Some(obj) = v.as_object() {
        let p11 = obj.get("p11").or_else(|| obj.get("p_11")).and_then(Value::as_f64).unwrap_or(0.0);
        let p12 = obj.get("p12").or_else(|| obj.get("p_12")).and_then(Value::as_f64).unwrap_or(0.0);
        let p44 = obj.get("p44").or_else(|| obj.get("p_44")).and_then(Value::as_f64).unwrap_or(0.5 * (p11 - p12));
        let p13 = obj.get("p13").or_else(|| obj.get("p_13")).and_then(Value::as_f64).unwrap_or(p12);
        let p31 = obj.get("p31").or_else(|| obj.get("p_31")).and_then(Value::as_f64).unwrap_or(p12);
        let p33 = obj.get("p33").or_else(|| obj.get("p_33")).and_then(Value::as_f64).unwrap_or(p11);
        let mut m = [[0.0; 6]; 6];
        m[0][0] = p11;
        m[1][1] = p11;
        m[2][2] = p33;
        m[0][1] = p12;
        m[1][0] = p12;
        m[0][2] = p13;
        m[1][2] = p13;
        m[2][0] = p31;
        m[2][1] = p31;
        m[3][3] = p44;
        m[4][4] = p44;
        m[5][5] = p44;
        return Ok(m);
    }
    if let Some(arr) = v.as_array().filter(|a| a.len() == 6) {
        let mut m = [[0.0; 6]; 6];
        for (r, row) in arr.iter().enumerate() {
            let Some(ra) = row.as_array().filter(|a| a.len() == 6) else {
                return Err(Error::invalid(format!("{name}: expected 6x6 matrix")));
            };
            for c_idx in 0..6 {
                m[r][c_idx] = ra[c_idx]
                    .as_f64()
                    .filter(|v| v.is_finite())
                    .ok_or_else(|| Error::invalid(format!("{name}[{r}][{c_idx}] must be finite")))?;
            }
        }
        return Ok(m);
    }
    Err(Error::invalid(format!("{name}: expected {{\"p11\": ..., \"p12\": ..., \"p44\": ...}} or 6x6 matrix")))
}

impl OptomechanicsOptions {
    /// Parses [`OptomechanicsOptions`] from a JSON parameter object.
    pub fn from_params(params: &Value) -> Result<Self, Error> {
        let mut opts = Self::default();
        let Some(obj) = params.as_object() else {
            return Ok(opts);
        };
        if let Some(v) = obj.get("optical_mode").or_else(|| obj.get("optical_mode_index")).filter(|v| !v.is_null()) {
            let idx = v.as_u64().ok_or_else(|| Error::invalid("parameter 'optical_mode' must be a positive integer"))?
                as usize;
            opts.optical_mode_index = if obj.contains_key("optical_mode") { idx.saturating_sub(1) } else { idx };
        }
        if let Some(v) = obj.get("stokes_mode").or_else(|| obj.get("stokes_mode_index")).filter(|v| !v.is_null()) {
            let idx =
                v.as_u64().ok_or_else(|| Error::invalid("parameter 'stokes_mode' must be a non-negative integer"))?
                    as usize;
            if obj.contains_key("stokes_mode") {
                if idx > 0 {
                    opts.stokes_mode_index = Some(idx - 1);
                }
            } else {
                opts.stokes_mode_index = Some(idx);
            }
        }
        if let Some(v) =
            obj.get("mechanical_mode").or_else(|| obj.get("mechanical_mode_index")).filter(|v| !v.is_null())
        {
            let idx =
                v.as_u64().ok_or_else(|| Error::invalid("parameter 'mechanical_mode' must be a positive integer"))?
                    as usize;
            opts.mechanical_mode_index = if obj.contains_key("mechanical_mode") { idx.saturating_sub(1) } else { idx };
        }
        if let Some(v) = obj.get("length_m").or_else(|| obj.get("length_y_m")).filter(|v| !v.is_null()) {
            let l = v
                .as_f64()
                .filter(|&x| x > 0.0 && x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'length_m' must be a positive finite length in metres"))?;
            opts.length_m = l;
        }
        if let Some(v) = obj.get("optical_freq_hz").or_else(|| obj.get("freq_hz")).filter(|v| !v.is_null()) {
            let f = v
                .as_f64()
                .filter(|&x| x > 0.0 && x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'optical_freq_hz' must be a positive finite frequency"))?;
            opts.optical_freq_hz = Some(f);
        }
        if let Some(v) = obj.get("q_mech").filter(|v| !v.is_null()) {
            let q = v
                .as_f64()
                .filter(|&x| x >= 0.0 && x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'q_mech' must be a non-negative finite number"))?;
            if q > 0.0 {
                opts.q_mech = Some(q);
            }
        }
        if let Some(v) = obj.get("v_g_opt").filter(|v| !v.is_null()) {
            let vg = v
                .as_f64()
                .filter(|&x| x > 0.0 && x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'v_g_opt' must be a positive finite velocity"))?;
            opts.v_g_opt = Some(vg);
        }
        if let Some(dm) = obj.get("domain_materials").and_then(Value::as_object) {
            for (k, mv) in dm {
                let Some(mobj) = mv.as_object() else {
                    return Err(Error::invalid(format!("domain_materials['{k}'] must be an object")));
                };
                let eps_r = match mobj.get("eps_r") {
                    Some(ev) => parse_m3_value(ev, &format!("domain_materials['{k}'].eps_r"))?,
                    None => iso_tensor(1.0),
                };
                let photoelastic = match mobj.get("photoelastic") {
                    Some(pv) if !pv.is_null() => {
                        Some(parse_m6_value(pv, &format!("domain_materials['{k}'].photoelastic"))?)
                    }
                    _ if mobj.contains_key("p11") || mobj.contains_key("p12") || mobj.contains_key("p44") => {
                        Some(parse_m6_value(mv, &format!("domain_materials['{k}']"))?)
                    }
                    _ => None,
                };
                let rho = mobj.get("rho").and_then(Value::as_f64).unwrap_or(0.0);
                let mat = OptomechanicsDomainMat { eps_r, photoelastic, rho };
                if let Ok(id) = k.parse::<i32>() {
                    opts.domain_materials_by_id.insert(id, mat.clone());
                }
                opts.domain_materials.insert(k.clone(), mat);
            }
        }
        Ok(opts)
    }

    fn resolve_opt_mat(&self, dom_id: i32) -> OptomechanicsDomainMat {
        if let Some(name) = self.opt_domain_names.get(&dom_id) {
            if let Some(m) = self.domain_materials.get(name) {
                return m.clone();
            }
        }
        if let Some(m) = self.domain_materials_by_id.get(&dom_id) {
            return m.clone();
        }
        OptomechanicsDomainMat::default()
    }

    fn resolve_mech_rho(&self, dom_id: i32) -> f64 {
        if let Some(name) = self.mech_domain_names.get(&dom_id) {
            if let Some(m) = self.domain_materials.get(name) {
                if m.rho > 0.0 {
                    return m.rho;
                }
            }
        }
        if let Some(m) = self.domain_materials_by_id.get(&dom_id) {
            if m.rho > 0.0 {
                return m.rho;
            }
        }
        0.0
    }
}

/// Complete output of [`optomechanics_with_meshes`].
#[derive(Clone, Debug, PartialEq)]
pub struct OptomechanicsReport {
    /// 1-based optical mode index.
    pub optical_mode: usize,
    /// 1-based mechanical mode index.
    pub mechanical_mode: usize,
    /// Optical mode frequency $f_{\text{opt}}$ ($\text{Hz}$).
    pub optical_freq_hz: f64,
    /// Mechanical mode frequency $f_m = \Omega_m / (2\pi)$ ($\text{Hz}$).
    pub mechanical_freq_hz: f64,
    /// Modal effective mass $m_{\text{eff}}$ ($\text{kg}$).
    pub m_eff_kg: f64,
    /// Zero-point fluctuation amplitude $x_{\text{zpf}} = \sqrt{\hbar / (2 m_{\text{eff}} \Omega_m)}$ ($\text{m}$).
    pub x_zpf_m: f64,
    /// Normalized optical mode stored electric energy denominator $U_{\text{opt}}$ ($\text{J}$).
    pub u_opt_j: f64,
    /// Photoelastic frequency pull sensitivity $G_{\text{PE}} = \mathrm{d}\omega_{\text{opt,PE}}/\mathrm{d}x$ ($\text{rad}/(\text{s}\cdot\text{m})$).
    pub g_om_pe_rad_s_per_m: f64,
    /// Moving-boundary frequency pull sensitivity $G_{\text{MB}} = \mathrm{d}\omega_{\text{opt,MB}}/\mathrm{d}x$ ($\text{rad}/(\text{s}\cdot\text{m})$).
    pub g_om_mb_rad_s_per_m: f64,
    /// Total optomechanical frequency shift sensitivity $G_{\text{OM}} = G_{\text{PE}} + G_{\text{MB}}$ ($\text{rad}/(\text{s}\cdot\text{m})$).
    pub g_om_rad_s_per_m: f64,
    /// Total optomechanical frequency shift sensitivity $G_{\text{OM}} / (2\pi)$ ($\text{Hz/m}$).
    pub g_om_hz_per_m: f64,
    /// Photoelastic vacuum optomechanical coupling rate $g_{\text{PE}}$ ($\text{rad/s}$).
    pub g_pe_rad_s: f64,
    /// Photoelastic vacuum optomechanical coupling rate $g_{\text{PE}} / (2\pi)$ ($\text{Hz}$).
    pub g_pe_hz: f64,
    /// Moving-boundary vacuum optomechanical coupling rate $g_{\text{MB}}$ ($\text{rad/s}$).
    pub g_mb_rad_s: f64,
    /// Moving-boundary vacuum optomechanical coupling rate $g_{\text{MB}} / (2\pi)$ ($\text{Hz}$).
    pub g_mb_hz: f64,
    /// Total vacuum optomechanical coupling rate $g_0 = g_{\text{PE}} + g_{\text{MB}}$ ($\text{rad/s}$).
    pub g_0_rad_s: f64,
    /// Total vacuum optomechanical coupling rate $g_0 / (2\pi) = (g_{\text{PE}} + g_{\text{MB}}) / (2\pi)$ ($\text{Hz}$).
    pub g_0_hz: f64,
    /// Magnitude $|g_0|$ ($\text{rad/s}$).
    pub g_0_abs_rad_s: f64,
    /// Magnitude $|g_0| / (2\pi)$ ($\text{Hz}$).
    pub g_0_abs_hz: f64,
    /// Mechanical quality factor $Q_m$ used for Brillouin gain.
    pub q_mech: f64,
    /// Stimulated Brillouin scattering gain $G_B$ ($\text{W}^{-1}\text{m}^{-1}$) for waveguide modes.
    pub sbs_gain_w_inv_m_inv: Option<f64>,
    /// Number of optical core cells integrated for $g_{\text{PE}}$.
    pub active_core_cells: usize,
    /// Number of dielectric interface facets integrated for $g_{\text{MB}}$.
    pub interface_facets: usize,
    /// True when the optical and mechanical discretisations are distinct meshes transferred via `TransferIndex` + `BoundaryIndex`.
    pub cross_mesh: bool,
}

impl OptomechanicsReport {
    pub fn to_json(&self) -> Value {
        let mut v = json!({
            "optical_mode": self.optical_mode,
            "mechanical_mode": self.mechanical_mode,
            "optical_freq_hz": self.optical_freq_hz,
            "mechanical_freq_hz": self.mechanical_freq_hz,
            "m_eff_kg": self.m_eff_kg,
            "x_zpf_m": self.x_zpf_m,
            "u_opt_j": self.u_opt_j,
            "g_om_pe_rad_s_per_m": self.g_om_pe_rad_s_per_m,
            "g_om_mb_rad_s_per_m": self.g_om_mb_rad_s_per_m,
            "g_om_rad_s_per_m": self.g_om_rad_s_per_m,
            "g_om_hz_per_m": self.g_om_hz_per_m,
            "g_pe_rad_s": self.g_pe_rad_s,
            "g_pe_hz": self.g_pe_hz,
            "g_mb_rad_s": self.g_mb_rad_s,
            "g_mb_hz": self.g_mb_hz,
            "g_0_rad_s": self.g_0_rad_s,
            "g_0_hz": self.g_0_hz,
            "g_0_abs_rad_s": self.g_0_abs_rad_s,
            "g_0_abs_hz": self.g_0_abs_hz,
            "q_mech": self.q_mech,
            "active_core_cells": self.active_core_cells,
            "interface_facets": self.interface_facets,
            "cross_mesh": self.cross_mesh,
        });
        if let Some(gb) = self.sbs_gain_w_inv_m_inv {
            if let Some(o) = v.as_object_mut() {
                o.insert("sbs_gain_w_inv_m_inv".into(), json!(gb));
            }
        }
        v
    }
}

/// Extracts per-cell complex electric field `[Ex, Ey, Ez]` (`nc * nloc`) in wafer axes from an optical [`Solution`].
fn extract_optical_cell_e(
    disc: &Discretisation,
    order: usize,
    sol: &Solution,
    is_real_waveguide_mode: bool,
) -> Result<Vec<[C; 3]>, Error> {
    let nc = disc.n_cells();
    let nloc = disc.nloc(order);
    let total = nc * nloc;

    let find_cf = |name: &str| -> Option<Vec<C>> {
        let cf = sol.cell_fields.iter().find(|f| f.name == name)?;
        if cf.re.len() != total {
            return None;
        }
        Some((0..total).map(|i| c(cf.re[i], cf.im.get(i).copied().unwrap_or(0.0))).collect())
    };

    let from_nodal = |name: &str, re_name: &str| -> Option<Vec<C>> {
        let vals = sol.flat(name).or_else(|| sol.flat(re_name))?;
        if vals.len() != disc.n_dofs(order) {
            return None;
        }
        let cv = femcore::fem::cell_values_at_nodes(disc, &vals, order);
        Some(cv.into_iter().map(|v| c(v, 0.0)).collect())
    };

    let ex = find_cf("Ex")
        .or_else(|| from_nodal("Ex", "Re Ex"))
        .ok_or_else(|| Error::invalid("optical solution has no 'Ex' field (solve with save_fields: true)"))?;
    let ez = find_cf("Ez")
        .or_else(|| from_nodal("Ez", "Re Ez"))
        .ok_or_else(|| Error::invalid("optical solution has no 'Ez' field (solve with save_fields: true)"))?;
    let ey_raw = find_cf("Ey").or_else(|| from_nodal("Ey", "Re Ey")).unwrap_or_else(|| vec![c(0.0, 0.0); total]);

    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        // In lossless 2D waveguide `modes`, `Ey` is stored as the real amplitude of `j beta phi` (90 deg out of phase with Ex, Ez).
        let ey = if is_real_waveguide_mode && ey_raw[i].im == 0.0 { c(0.0, ey_raw[i].re) } else { ey_raw[i] };
        out.push([ex[i], ey, ez[i]]);
    }
    Ok(out)
}

/// Extracts per-cell complex displacement `[ux, uy, uz]` (`nc * nloc`) in wafer axes from a mechanical [`Solution`].
fn extract_mechanical_cell_u(disc: &Discretisation, order: usize, sol: &Solution) -> Result<Vec<[C; 3]>, Error> {
    let nc = disc.n_cells();
    let nloc = disc.nloc(order);
    let total = nc * nloc;

    let find_cf = |name: &str| -> Option<Vec<C>> {
        let cf = sol.cell_fields.iter().find(|f| f.name == name)?;
        if cf.re.len() != total {
            return None;
        }
        Some((0..total).map(|i| c(cf.re[i], cf.im.get(i).copied().unwrap_or(0.0))).collect())
    };

    let from_nodal = |name: &str, re_name: &str| -> Option<Vec<C>> {
        let vals = sol.flat(name).or_else(|| sol.flat(re_name))?;
        if vals.len() != disc.n_dofs(order) {
            return None;
        }
        let cv = femcore::fem::cell_values_at_nodes(disc, &vals, order);
        Some(cv.into_iter().map(|v| c(v, 0.0)).collect())
    };

    let ux = find_cf("ux")
        .or_else(|| from_nodal("ux", "Re ux"))
        .ok_or_else(|| Error::invalid("mechanical solution has no 'ux' displacement field"))?;
    let uz = find_cf("uz")
        .or_else(|| from_nodal("uz", "Re uz"))
        .ok_or_else(|| Error::invalid("mechanical solution has no 'uz' displacement field"))?;
    let uy = find_cf("uy").or_else(|| from_nodal("uy", "Re uy")).unwrap_or_else(|| vec![c(0.0, 0.0); total]);

    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        out.push([ux[i], uy[i], uz[i]]);
    }
    Ok(out)
}

fn infer_optical_freq_hz(
    opt_result: &SolveResult,
    opt_sol: &Solution,
    opts: &OptomechanicsOptions,
) -> Result<f64, Error> {
    if let Some(f) = opts.optical_freq_hz.filter(|&f| f > 0.0 && f.is_finite()) {
        return Ok(f);
    }
    if matches!(opt_result.study.as_str(), "cavity" | "axisymmetric_modes" | "driven") {
        if let Some(f) = opt_sol.value.filter(|&f| f > 0.0 && f.is_finite()) {
            return Ok(f);
        }
        if let Some([re, _]) = opt_sol.eigenvalue.filter(|e| e[0] > 0.0 && e[0].is_finite()) {
            return Ok(re);
        }
    }
    if let Some(f) = opt_result.scalar("freq").or_else(|| opt_result.scalar("f1")).filter(|&f| f > 0.0 && f.is_finite())
    {
        return Ok(f);
    }
    // Default telecom 1550 nm (193.414489 THz) for 2D waveguide mode results when study frequency is not embedded in SolveResult
    Ok(C0 / 1.55e-6)
}

fn infer_mechanical_freq_hz(mech_result: &SolveResult, mech_sol: &Solution, mode_idx: usize) -> Result<f64, Error> {
    if let Some(f) = mech_sol.value.filter(|&f| f > 0.0 && f.is_finite()) {
        return Ok(f);
    }
    if let Some([re, _]) = mech_sol.eigenvalue.filter(|e| e[0] > 0.0 && e[0].is_finite()) {
        return Ok(re);
    }
    let key_re = format!("f re {}", mode_idx + 1);
    let key_f = format!("f{}", mode_idx + 1);
    if let Some(f) =
        mech_result.scalar(&key_re).or_else(|| mech_result.scalar(&key_f)).filter(|&f| f > 0.0 && f.is_finite())
    {
        return Ok(f);
    }
    Err(Error::invalid(format!(
        "mechanical solution {} has no positive eigenfrequency (expected an eigenfrequency, modes, or driven study)",
        mode_idx + 1
    )))
}

fn infer_q_mech(mech_result: &SolveResult, mode_idx: usize, opts: &OptomechanicsOptions) -> f64 {
    if let Some(q) = opts.q_mech.filter(|&q| q > 0.0 && q.is_finite()) {
        return q;
    }
    let q_key = format!("Q{}", mode_idx + 1);
    if let Some(q) = mech_result.scalar(&q_key).filter(|&q| q > 0.0 && q.is_finite()) {
        return q;
    }
    let f_re = mech_result.scalar(&format!("f re {}", mode_idx + 1)).unwrap_or(0.0);
    let f_im = mech_result.scalar(&format!("f im {}", mode_idx + 1)).unwrap_or(0.0);
    if f_re > 0.0 && f_im.abs() > 1e-30 {
        let q = f_re / (2.0 * f_im.abs());
        if q > 0.0 && q.is_finite() {
            return q;
        }
    }
    1000.0
}

fn infer_v_g(
    opt_result: &SolveResult,
    opt_sol: &Solution,
    mode_idx: usize,
    opts: &OptomechanicsOptions,
) -> Option<f64> {
    if let Some(vg) = opts.v_g_opt.filter(|&v| v > 0.0 && v.is_finite()) {
        return Some(vg);
    }
    let vg_key = format!("v_g {}", mode_idx + 1);
    if let Some(vg) = opt_result.scalar(&vg_key).filter(|&v| v > 0.0 && v.is_finite()) {
        return Some(vg);
    }
    let ng_key = format!("n_g {}", mode_idx + 1);
    if let Some(ng) = opt_result.scalar(&ng_key).filter(|&n| n > 0.0 && n.is_finite()) {
        return Some(C0 / ng);
    }
    if opt_result.study == "modes" {
        if let Some([neff, _]) = opt_sol.eigenvalue.filter(|e| e[0] > 0.0 && e[0].is_finite()) {
            return Some(C0 / neff);
        }
    }
    None
}

/// Computes barycentric coordinates $\boldsymbol{\lambda}$ of physical point `x_q` inside simplex cell `ci`.
fn point_barycentric_in_cell(disc: &Discretisation, geo: &femcore::fem::Geometry, ci: usize, x_q: &[f64]) -> Vec<f64> {
    let d = disc.dim;
    let cv = disc.cell(ci);
    let base = ci * (d + 1) * d;
    let mut lam = vec![0.0; d + 1];
    let p0 = disc.point(cv[0] as usize);
    let mut dx = [0.0; 3];
    for k in 0..d {
        dx[k] = x_q[k] - p0[k];
    }
    let mut sum_other = 0.0;
    for (v_idx, lam_v) in lam.iter_mut().enumerate().take(d + 1).skip(1) {
        let mut lv = 0.0;
        for (k, &dx_k) in dx.iter().enumerate().take(d) {
            lv += geo.g[base + v_idx * d + k] * dx_k;
        }
        *lam_v = lv;
        sum_other += lv;
    }
    lam[0] = 1.0 - sum_other;
    lam
}

/// Evaluates cross-mesh optomechanical coupling rates ($m_{\text{eff}}, x_{\text{zpf}}, G_{\text{OM}}, g_{\text{PE}}, g_{\text{MB}}, g_0, G_B$)
/// between an optical [`SolveResult`] on `opt_mesh` and a mechanical [`SolveResult`] on `mech_mesh`.
pub fn optomechanics_with_meshes(
    opt_mesh: &Mesh,
    opt_result: &SolveResult,
    mech_mesh: &Mesh,
    mech_result: &SolveResult,
    opts: &OptomechanicsOptions,
) -> Result<OptomechanicsReport, Error> {
    let opt_disc = discretise(opt_mesh)?;
    let mech_disc = discretise(mech_mesh)?;
    if opt_disc.dim != mech_disc.dim {
        return Err(Error::invalid(format!(
            "optical mesh dimension ({}D) does not match mechanical mesh dimension ({}D)",
            opt_disc.dim, mech_disc.dim
        )));
    }
    let dim = opt_disc.dim;
    let opt_order = (opt_result.order as usize).clamp(1, 2);
    let mech_order = (mech_result.order as usize).clamp(1, 2);
    let ly = if dim == 2 { opts.length_m } else { 1.0 };

    let opt_sol = opt_result.solutions.get(opts.optical_mode_index).ok_or_else(|| {
        Error::invalid(format!(
            "optical_mode {} is out of bounds (optical result has {} solution(s))",
            opts.optical_mode_index + 1,
            opt_result.solutions.len()
        ))
    })?;
    let stokes_idx = opts.stokes_mode_index.unwrap_or(opts.optical_mode_index);
    let stokes_sol = opt_result.solutions.get(stokes_idx).ok_or_else(|| {
        Error::invalid(format!(
            "stokes_mode {} is out of bounds (optical result has {} solution(s))",
            stokes_idx + 1,
            opt_result.solutions.len()
        ))
    })?;
    let mech_sol = mech_result.solutions.get(opts.mechanical_mode_index).ok_or_else(|| {
        Error::invalid(format!(
            "mechanical_mode {} is out of bounds (mechanical result has {} solution(s))",
            opts.mechanical_mode_index + 1,
            mech_result.solutions.len()
        ))
    })?;

    let f_opt = infer_optical_freq_hz(opt_result, opt_sol, opts)?;
    let omega_opt = 2.0 * PI * f_opt;
    let f_mech = infer_mechanical_freq_hz(mech_result, mech_sol, opts.mechanical_mode_index)?;
    let omega_mech = 2.0 * PI * f_mech;

    let opt_outside: BTreeSet<i32> =
        FieldFrame::of(opt_result).map(|f| f.outside.into_iter().collect()).unwrap_or_default();
    let mech_outside: BTreeSet<i32> =
        FieldFrame::of(mech_result).map(|f| f.outside.into_iter().collect()).unwrap_or_default();

    // 1. Extract and canonicalize mechanical mode displacement u(x) on mech_disc
    let mut u_cell = extract_mechanical_cell_u(&mech_disc, mech_order, mech_sol)?;
    let nc_mech = mech_disc.n_cells();
    let nloc_mech = mech_disc.nloc(mech_order);
    let mech_geo = geometry(&mech_disc);
    let mech_re = RefElement::new(dim, mech_order, 2 * mech_order);

    // Find peak 3D displacement norm u_max outside PML and dominant complex phasor
    let mut u_max_sq = 0.0_f64;
    let mut dom_phasor = c(1.0, 0.0);
    for ci in 0..nc_mech {
        let dom = mech_disc.cell_labels[ci];
        if mech_outside.contains(&dom) || opts.resolve_mech_rho(dom) <= 0.0 {
            continue;
        }
        for a in 0..nloc_mech {
            let u = u_cell[ci * nloc_mech + a];
            let n2 = u[0].norm_sqr() + u[1].norm_sqr() + u[2].norm_sqr();
            if n2 > u_max_sq {
                u_max_sq = n2;
                let best_c = [u[0], u[1], u[2]]
                    .into_iter()
                    .max_by(|x, y| x.norm_sqr().partial_cmp(&y.norm_sqr()).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(c(1.0, 0.0));
                if best_c.norm() > 1e-30 {
                    dom_phasor = best_c / best_c.norm();
                }
            }
        }
    }
    let u_max = u_max_sq.sqrt();
    if !(u_max > 0.0 && u_max.is_finite()) {
        return Err(Error::invalid(
            "mechanical mode has zero displacement on active solid domains (check `rho` in domain_materials)",
        ));
    }

    // Mesh-independent physical sign gauge: compute solid centroid x_cm and radial breathing integral
    // I_div = \int \rho (x - x_cm) . Re(e^{-i theta} u) dV (fallback to weighted moment I_1).
    let mut mass_tot = 0.0_f64;
    let mut x_cm = [0.0_f64; 3];
    for ci in 0..nc_mech {
        let dom = mech_disc.cell_labels[ci];
        if mech_outside.contains(&dom) {
            continue;
        }
        let rho = opts.resolve_mech_rho(dom);
        if rho <= 0.0 {
            continue;
        }
        let cv = mech_disc.cell(ci);
        let mut xc = [0.0; 3];
        for &v in cv {
            let p = mech_disc.point(v as usize);
            for k in 0..dim {
                xc[k] += p[k] / (cv.len() as f64);
            }
        }
        let dm = rho * mech_geo.vol[ci] * ly;
        mass_tot += dm;
        for k in 0..dim {
            x_cm[k] += dm * xc[k];
        }
    }
    if mass_tot > 0.0 {
        for xk in x_cm.iter_mut().take(dim) {
            *xk /= mass_tot;
        }
    }

    let rot0 = dom_phasor.conj() / u_max;
    let mut i_div = 0.0_f64;
    let mut i_mom = 0.0_f64;
    for ci in 0..nc_mech {
        let dom = mech_disc.cell_labels[ci];
        if mech_outside.contains(&dom) {
            continue;
        }
        let rho = opts.resolve_mech_rho(dom);
        if rho <= 0.0 {
            continue;
        }
        let cv = mech_disc.cell(ci);
        let mut xc = [0.0; 3];
        for &v in cv {
            let p = mech_disc.point(v as usize);
            for k in 0..dim {
                xc[k] += p[k] / (cv.len() as f64);
            }
        }
        let xc_w = femcore::section::wafer_point(&xc[..dim], dim);
        let xcm_w = femcore::section::wafer_point(&x_cm[..dim], dim);
        let mut u_avg = [0.0_f64; 3];
        for a in 0..nloc_mech {
            let u = u_cell[ci * nloc_mech + a];
            for k in 0..3 {
                u_avg[k] += (u[k] * rot0).re / (nloc_mech as f64);
            }
        }
        let w = rho * mech_geo.vol[ci] * ly;
        i_div +=
            w * ((xc_w[0] - xcm_w[0]) * u_avg[0] + (xc_w[1] - xcm_w[1]) * u_avg[1] + (xc_w[2] - xcm_w[2]) * u_avg[2]);
        let r_wt = (xc_w[0] - xcm_w[0])
            + std::f64::consts::SQRT_2 * (xc_w[1] - xcm_w[1])
            + 3.0_f64.sqrt() * (xc_w[2] - xcm_w[2]);
        let u_wt = u_avg[0] + std::f64::consts::SQRT_2 * u_avg[1] + 3.0_f64.sqrt() * u_avg[2];
        i_mom += w * r_wt * u_wt;
    }
    let sign_flip = if i_div.abs() > 1e-12 * mass_tot * 1e-6 {
        if i_div < 0.0 {
            -1.0
        } else {
            1.0
        }
    } else if i_mom < 0.0 {
        -1.0
    } else {
        1.0
    };
    let norm_scale = rot0 * sign_flip;
    for u in &mut u_cell {
        u[0] *= norm_scale;
        u[1] *= norm_scale;
        u[2] *= norm_scale;
    }

    // 2. Compute modal effective mass m_eff = \int \rho |\tilde{u}|^2 dV and x_zpf
    let mut m_eff = 0.0_f64;
    for ci in 0..nc_mech {
        let dom = mech_disc.cell_labels[ci];
        if mech_outside.contains(&dom) {
            continue;
        }
        let rho = opts.resolve_mech_rho(dom);
        if rho <= 0.0 {
            continue;
        }
        let vol_scale = rho * mech_geo.vol[ci] * ly;
        for q in 0..mech_re.w.len() {
            let mut u_q = [c(0.0, 0.0); 3];
            for a in 0..nloc_mech {
                let na = mech_re.n[q][a];
                let ua = u_cell[ci * nloc_mech + a];
                u_q[0] += ua[0] * na;
                u_q[1] += ua[1] * na;
                u_q[2] += ua[2] * na;
            }
            let u_norm_sq = u_q[0].norm_sqr() + u_q[1].norm_sqr() + u_q[2].norm_sqr();
            m_eff += vol_scale * mech_re.w[q] * u_norm_sq;
        }
    }
    if !(m_eff > 0.0 && m_eff.is_finite()) {
        return Err(Error::invalid("modal effective mass m_eff is zero or non-finite"));
    }
    let x_zpf = (HBAR / (2.0 * m_eff * omega_mech)).sqrt();

    // 3. Build cell-wise normalized strain \tilde{S} (6-component solver Voigt `[xx, yy, zz, xy, yz, xz]`, real part)
    // and normalized displacement \tilde{u} (`[ux, uy, uz]`, real part) on `mech_disc` for TransferIndex & BoundaryIndex
    let lam_nodes_mech: Vec<Vec<f64>> = match (dim, mech_order) {
        (2, 1) => vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]],
        (2, _) => vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![0.5, 0.5, 0.0],
            vec![0.5, 0.0, 0.5],
            vec![0.0, 0.5, 0.5],
        ],
        (_, 1) => {
            vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0], vec![0.0, 0.0, 1.0, 0.0], vec![0.0, 0.0, 0.0, 1.0]]
        }
        _ => vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
            vec![0.5, 0.5, 0.0, 0.0],
            vec![0.5, 0.0, 0.5, 0.0],
            vec![0.5, 0.0, 0.0, 0.5],
            vec![0.0, 0.5, 0.5, 0.0],
            vec![0.0, 0.5, 0.0, 0.5],
            vec![0.0, 0.0, 0.5, 0.5],
        ],
    };
    let shape_grads_mech: Vec<Vec<f64>> = lam_nodes_mech.iter().map(|l| shape(l, mech_order).1).collect();
    let mut mech_strain_cell_vals = vec![0.0_f64; nc_mech * nloc_mech * 6];
    let mut mech_u_cell_vals = vec![0.0_f64; nc_mech * nloc_mech * 3];

    for ci in 0..nc_mech {
        let base_g = ci * (dim + 1) * dim;
        for a_idx in 0..nloc_mech {
            let u_here = u_cell[ci * nloc_mech + a_idx];
            let u_slot = (ci * nloc_mech + a_idx) * 3;
            mech_u_cell_vals[u_slot] = u_here[0].re;
            mech_u_cell_vals[u_slot + 1] = u_here[1].re;
            mech_u_cell_vals[u_slot + 2] = u_here[2].re;

            let dn = &shape_grads_mech[a_idx];
            let mut g_u = [[0.0_f64; 3]; 3];
            for node_i in 0..nloc_mech {
                let mut gvec = [0.0_f64; 3];
                for (k, gvec_k) in gvec.iter_mut().enumerate().take(dim) {
                    let mut ga = 0.0;
                    for v_idx in 0..=dim {
                        ga += dn[node_i * (dim + 1) + v_idx] * mech_geo.g[base_g + v_idx * dim + k];
                    }
                    *gvec_k = ga;
                }
                let un = u_cell[ci * nloc_mech + node_i];
                if dim == 2 {
                    // In 2D cross-section: mesh axis 0 = wafer x, mesh axis 1 = wafer z
                    g_u[0][0] += un[0].re * gvec[0];
                    g_u[0][2] += un[0].re * gvec[1];
                    g_u[1][0] += un[1].re * gvec[0];
                    g_u[1][2] += un[1].re * gvec[1];
                    g_u[2][0] += un[2].re * gvec[0];
                    g_u[2][2] += un[2].re * gvec[1];
                } else {
                    for comp in 0..3 {
                        for k in 0..3 {
                            g_u[comp][k] += un[comp].re * gvec[k];
                        }
                    }
                }
            }
            let s_wafer =
                [g_u[0][0], g_u[1][1], g_u[2][2], g_u[0][1] + g_u[1][0], g_u[1][2] + g_u[2][1], g_u[0][2] + g_u[2][0]];
            let s_slot = (ci * nloc_mech + a_idx) * 6;
            mech_strain_cell_vals[s_slot..s_slot + 6].copy_from_slice(&s_wafer);
        }
    }

    // Build nodal displacement vector on mech_disc (`n1_mech * 3`, interleaved) for BoundaryIndex
    let n1_mech = mech_disc.n_dofs(mech_order);
    let mut mech_u_nodal_interleaved = vec![0.0_f64; n1_mech * 3];
    let mut mech_u_cnt = vec![0usize; n1_mech];
    let mut dofs_buf = Vec::new();
    for ci in 0..nc_mech {
        let dom = mech_disc.cell_labels[ci];
        if mech_outside.contains(&dom) || opts.resolve_mech_rho(dom) <= 0.0 {
            continue;
        }
        mech_disc.cell_dofs(ci, mech_order, &mut dofs_buf);
        for (a, &gd) in dofs_buf.iter().take(nloc_mech).enumerate() {
            let u = u_cell[ci * nloc_mech + a];
            mech_u_nodal_interleaved[3 * gd] += u[0].re;
            mech_u_nodal_interleaved[3 * gd + 1] += u[1].re;
            mech_u_nodal_interleaved[3 * gd + 2] += u[2].re;
            mech_u_cnt[gd] += 1;
        }
    }
    for gd in 0..n1_mech {
        if mech_u_cnt[gd] > 1 {
            let inv = 1.0 / (mech_u_cnt[gd] as f64);
            mech_u_nodal_interleaved[3 * gd] *= inv;
            mech_u_nodal_interleaved[3 * gd + 1] *= inv;
            mech_u_nodal_interleaved[3 * gd + 2] *= inv;
        }
    }

    let is_cross_mesh = !same_simplex_mesh(&opt_disc, &mech_disc);
    let mech_transfer_idx =
        TransferIndex::new(SourceMesh::Simplex(&mech_disc), &opts.mech_domain_names, TransferOptions::default())
            .map_err(|e| Error::from_kernel(e.to_string()))?;
    let mech_boundary_idx = BoundaryIndex::new(&mech_disc, &opts.mech_boundary_names, TransferOptions::default())
        .map_err(|e| Error::from_kernel(e.to_string()))?;

    let u_nodal_source = SourceField::Nodal { values: &mech_u_nodal_interleaved, ncomp: 3, order: mech_order };

    // 4. Extract optical mode electric fields E_1(x) and E_2(x) on opt_disc and compute U_opt
    let is_real_wg = opt_result.study == "modes"
        && opt_sol.cell_fields.iter().find(|f| f.name == "Ey").is_some_and(|f| f.im.is_empty());
    let e1_cell = extract_optical_cell_e(&opt_disc, opt_order, opt_sol, is_real_wg)?;
    let e2_cell = if stokes_idx == opts.optical_mode_index {
        e1_cell.clone()
    } else {
        extract_optical_cell_e(&opt_disc, opt_order, stokes_sol, is_real_wg)?
    };

    let nc_opt = opt_disc.n_cells();
    let nloc_opt = opt_disc.nloc(opt_order);
    let opt_geo = geometry(&opt_disc);
    let opt_re = RefElement::new(dim, opt_order, 2 * opt_order + 2);
    let qr_opt = quad_rule(dim, 2 * opt_order + 2);

    let mut u_opt_1 = 0.0_f64;
    let mut u_opt_2 = 0.0_f64;
    let mut pe_integral = 0.0_f64;
    let mut active_core_cells = 0usize;

    for ci in 0..nc_opt {
        let dom = opt_disc.cell_labels[ci];
        if opt_outside.contains(&dom) {
            continue;
        }
        let mat = opts.resolve_opt_mat(dom);
        let eps_r = &mat.eps_r;
        let vol_scale = opt_geo.vol[ci] * ly;

        for q in 0..opt_re.w.len() {
            let mut e1_q = [c(0.0, 0.0); 3];
            let mut e2_q = [c(0.0, 0.0); 3];
            for a in 0..nloc_opt {
                let na = opt_re.n[q][a];
                let v1 = e1_cell[ci * nloc_opt + a];
                let v2 = e2_cell[ci * nloc_opt + a];
                for k in 0..3 {
                    e1_q[k] += v1[k] * na;
                    e2_q[k] += v2[k] * na;
                }
            }
            let mut e1_eps_e1 = 0.0_f64;
            let mut e2_eps_e2 = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    e1_eps_e1 += (e1_q[a].conj() * e1_q[b]).re * eps_r[a][b];
                    e2_eps_e2 += (e2_q[a].conj() * e2_q[b]).re * eps_r[a][b];
                }
            }
            u_opt_1 += vol_scale * opt_re.w[q] * EPS0 * e1_eps_e1;
            u_opt_2 += vol_scale * opt_re.w[q] * EPS0 * e2_eps_e2;
        }

        // Cross-mesh photoelastic integral over cells possessing a photoelastic tensor p_IJ
        let Some(ref p_glob) = mat.photoelastic else {
            continue;
        };
        let dom_filter = opts.opt_domain_names.get(&dom).map(String::as_str);
        if let Some(dname) = dom_filter {
            if !opts.mech_domain_names.is_empty() && !opts.mech_domain_names.values().any(|v| v == dname) {
                continue;
            }
        }
        let cv = opt_disc.cell(ci);
        let nq = opt_re.w.len();
        let mut cell_touched = false;
        for q in 0..nq {
            let mut x_q = [0.0_f64; 3];
            let lam = &qr_opt.lam[q];
            for (v_i, &v) in cv.iter().enumerate() {
                let pv = opt_disc.point(v as usize);
                for k in 0..dim {
                    x_q[k] += lam[v_i] * pv[k];
                }
            }
            let Some(loc) = mech_transfer_idx.locate(&x_q[..dim], dom_filter) else {
                continue;
            };
            cell_touched = true;
            let (n_mech, _) = shape(&loc.coords[..=dim], mech_order);
            let mut s_q = [0.0_f64; 6];
            let base = loc.cell * nloc_mech;
            for (a, &na) in n_mech.iter().take(nloc_mech).enumerate() {
                let off = (base + a) * 6;
                for comp in 0..6 {
                    s_q[comp] += na * mech_strain_cell_vals[off + comp];
                }
            }
            let d_eps = delta_eps_eo_pe(eps_r, None, Some(p_glob), None, Some(&s_q));

            let mut e1_q = [c(0.0, 0.0); 3];
            let mut e2_q = [c(0.0, 0.0); 3];
            for a in 0..nloc_opt {
                let na = opt_re.n[q][a];
                let v1 = e1_cell[ci * nloc_opt + a];
                let v2 = e2_cell[ci * nloc_opt + a];
                for k in 0..3 {
                    e1_q[k] += v1[k] * na;
                    e2_q[k] += v2[k] * na;
                }
            }
            let mut e_deps_e = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    e_deps_e += (e1_q[a].conj() * e2_q[b]).re * d_eps[a][b];
                }
            }
            pe_integral += vol_scale * opt_re.w[q] * EPS0 * e_deps_e;
        }
        if cell_touched {
            active_core_cells += 1;
        }
    }

    let u_opt = (u_opt_1 * u_opt_2).sqrt();
    if !(u_opt > 0.0 && u_opt.is_finite()) {
        return Err(Error::invalid("optical mode energy U_opt is zero or non-finite"));
    }

    let g_om_pe = -omega_opt * pe_integral / (2.0 * u_opt);

    // 5. Cross-mesh moving-boundary (radiation-pressure) integral over dielectric interfaces of `opt_disc`
    // Collect interior interface facets `(c_solid, c_clad, facet_verts, opt_bname_if_any)`
    let mut facet_map: BTreeMap<[u32; 3], Vec<usize>> = BTreeMap::new();
    let fkey = |vs: &[u32]| -> [u32; 3] {
        let mut k = [u32::MAX; 3];
        k[..vs.len()].copy_from_slice(vs);
        k[..vs.len()].sort_unstable();
        k
    };
    for ci in 0..nc_opt {
        if opt_outside.contains(&opt_disc.cell_labels[ci]) {
            continue;
        }
        let cv = opt_disc.cell(ci);
        for skip in 0..=dim {
            let fv: Vec<u32> = (0..=dim).filter(|&i| i != skip).map(|i| cv[i]).collect();
            facet_map.entry(fkey(&fv)).or_default().push(ci);
        }
    }
    let mut explicit_facet_labels: BTreeMap<[u32; 3], i32> = BTreeMap::new();
    for fi in 0..opt_disc.n_facets() {
        let lbl = opt_disc.facet_labels[fi];
        if lbl > 0 {
            explicit_facet_labels.insert(fkey(opt_disc.facet(fi)), lbl);
        }
    }

    // 1D (2D mesh) or 2D (3D mesh) facet Gauss quadrature rule (3-point Gauss in 2D, degree-4 triangle rule in 3D)
    let (f_barys, f_weights): (Vec<Vec<f64>>, Vec<f64>) = if dim == 2 {
        let g = (0.6_f64).sqrt();
        let pts = [-g, 0.0, g];
        let ws = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        (pts.iter().map(|&t| vec![0.5 * (1.0 - t), 0.5 * (1.0 + t)]).collect(), ws.iter().map(|&w| 0.5 * w).collect())
    } else {
        let qr_f = quad_rule(2, 4);
        (qr_f.lam, qr_f.w)
    };

    let is_solid_opt_dom = |dom: i32| -> bool {
        if opts.resolve_mech_rho(dom) > 0.0 {
            return true;
        }
        if let Some(name) = opts.opt_domain_names.get(&dom) {
            if opts.domain_materials.get(name).is_some_and(|m| m.rho > 0.0) {
                return true;
            }
            if opts.mech_domain_names.values().any(|v| v == name) {
                return true;
            }
        }
        opts.domain_materials_by_id.get(&dom).is_some_and(|m| m.rho > 0.0)
    };

    let mut mb_integral = 0.0_f64;
    let mut interface_facets = 0usize;

    for (key, cells) in &facet_map {
        if cells.len() != 2 {
            continue;
        }
        let (ca, cb) = (cells[0], cells[1]);
        let (da, db) = (opt_disc.cell_labels[ca], opt_disc.cell_labels[cb]);
        if da == db {
            continue;
        }
        let sa = is_solid_opt_dom(da);
        let sb = is_solid_opt_dom(db);
        if !sa && !sb {
            continue;
        }
        // Order `(c1, d1)` as the solid side and `(c2, d2)` as the exterior/cladding side
        let (c1, d1, c2, d2) = if sa && !sb {
            (ca, da, cb, db)
        } else if sb && !sa {
            (cb, db, ca, da)
        } else if opts.resolve_opt_mat(da).eps_r[0][0] >= opts.resolve_opt_mat(db).eps_r[0][0] {
            (ca, da, cb, db)
        } else {
            (cb, db, ca, da)
        };

        let mat1 = opts.resolve_opt_mat(d1);
        let mat2 = opts.resolve_opt_mat(d2);
        let eps1 = &mat1.eps_r;
        let eps2 = &mat2.eps_r;
        let diff_eps =
            (0..3).flat_map(|i| (0..3).map(move |j| (eps1[i][j] - eps2[i][j]).abs())).fold(0.0_f64, f64::max);
        if diff_eps <= 1e-12 {
            continue;
        }

        let fv: Vec<u32> = key.iter().copied().filter(|&v| v != u32::MAX).collect();
        // Compute facet measure (length in 2D, area in 3D) and outward unit normal n_hat from cell c1 into cell c2 (in wafer axes)
        let (measure, n_wafer) = if dim == 2 {
            let p0 = opt_disc.point(fv[0] as usize);
            let p1 = opt_disc.point(fv[1] as usize);
            let (dx, dz) = (p1[0] - p0[0], p1[1] - p0[1]);
            let len = dx.hypot(dz);
            if len <= 1e-30 {
                continue;
            }
            let mut nx = dz / len;
            let mut nz = -dx / len;
            let cv1 = opt_disc.cell(c1);
            if let Some(&v_opp) = cv1.iter().find(|v| !fv.contains(v)) {
                let po = opt_disc.point(v_opp as usize);
                let mx = 0.5 * (p0[0] + p1[0]);
                let mz = 0.5 * (p0[1] + p1[1]);
                if (mx - po[0]) * nx + (mz - po[1]) * nz < 0.0 {
                    nx = -nx;
                    nz = -nz;
                }
            }
            (len * ly, [nx, 0.0, nz])
        } else {
            let p0 = opt_disc.point(fv[0] as usize);
            let p1 = opt_disc.point(fv[1] as usize);
            let p2 = opt_disc.point(fv[2] as usize);
            let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            let norm = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
            if norm <= 1e-30 {
                continue;
            }
            let area = 0.5 * norm;
            let mut n = [cr[0] / norm, cr[1] / norm, cr[2] / norm];
            let cv1 = opt_disc.cell(c1);
            if let Some(&v_opp) = cv1.iter().find(|v| !fv.contains(v)) {
                let po = opt_disc.point(v_opp as usize);
                let mc = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
                if (mc[0] - po[0]) * n[0] + (mc[1] - po[1]) * n[1] + (mc[2] - po[2]) * n[2] < 0.0 {
                    n = [-n[0], -n[1], -n[2]];
                }
            }
            (area, n)
        };

        let solid_dom_name = opts.opt_domain_names.get(&d1).map(String::as_str);
        let solid_mech_labels: Vec<i32> = match solid_dom_name {
            Some(sname) => {
                opts.mech_domain_names.iter().filter(|(_, n)| n.as_str() == sname).map(|(&id, _)| id).collect()
            }
            None => vec![d1],
        };
        let bname_opt: Option<&str> = explicit_facet_labels
            .get(key)
            .and_then(|lbl| opts.opt_boundary_names.get(lbl))
            .map(String::as_str)
            .filter(|bn| opts.mech_boundary_names.values().any(|mb| mb == *bn))
            .or_else(|| {
                let n1 = opts.opt_domain_names.get(&d1)?;
                let n2 = opts.opt_domain_names.get(&d2)?;
                let pair1 = format!("{n1}|{n2}");
                let pair2 = format!("{n2}|{n1}");
                opts.mech_boundary_names.values().find(|mb| **mb == pair1 || **mb == pair2).map(String::as_str)
            });

        // Normal permittivity projections n . eps_r1 . n and n . eps_r2 . n
        let mut eps1_nn = 0.0_f64;
        let mut eps2_nn = 0.0_f64;
        for a in 0..3 {
            for b in 0..3 {
                eps1_nn += n_wafer[a] * eps1[a][b] * n_wafer[b];
                eps2_nn += n_wafer[a] * eps2[a][b] * n_wafer[b];
            }
        }
        let inv_eps1_minus_inv_eps2 = (1.0 / (EPS0 * eps1_nn.max(1e-30))) - (1.0 / (EPS0 * eps2_nn.max(1e-30)));

        let mut facet_touched = false;
        for (bary_f, &wf) in f_barys.iter().zip(&f_weights) {
            let mut x_q = [0.0_f64; 3];
            for (i, &vid) in fv.iter().enumerate() {
                let pv = opt_disc.point(vid as usize);
                for k in 0..dim {
                    x_q[k] += bary_f[i] * pv[k];
                }
            }

            // Sample mechanical displacement \tilde{u}(x_q) via BoundaryIndex when a named boundary matches,
            // falling back to TransferIndex strictly filtered to the solid domain `solid_dom_name`.
            let u_vec = if let Some(bname) = bname_opt {
                let bev = mech_boundary_idx
                    .eval_boundary(&u_nodal_source, &x_q[..dim], bname, &solid_mech_labels)
                    .map_err(|e| Error::from_kernel(e.to_string()))?;
                if bev.found.first().copied() == Some(true) {
                    Some([bev.values[0], bev.values[1], bev.values[2]])
                } else {
                    None
                }
            } else {
                None
            }
            .or_else(|| {
                let loc = mech_transfer_idx.locate(&x_q[..dim], solid_dom_name)?;
                let (n_mech, _) = shape(&loc.coords[..=dim], mech_order);
                let mut u_q = [0.0_f64; 3];
                let base = loc.cell * nloc_mech;
                for (a, &na) in n_mech.iter().take(nloc_mech).enumerate() {
                    let off = (base + a) * 3;
                    u_q[0] += na * mech_u_cell_vals[off];
                    u_q[1] += na * mech_u_cell_vals[off + 1];
                    u_q[2] += na * mech_u_cell_vals[off + 2];
                }
                Some(u_q)
            });

            let Some(u_q) = u_vec else {
                continue;
            };
            facet_touched = true;
            let u_dot_n = u_q[0] * n_wafer[0] + u_q[1] * n_wafer[1] + u_q[2] * n_wafer[2];

            // Evaluate limiting optical fields E^(1)(x_q) in cell c1 and E^(2)(x_q) in cell c2
            let lam1 = point_barycentric_in_cell(&opt_disc, &opt_geo, c1, &x_q[..dim]);
            let lam2 = point_barycentric_in_cell(&opt_disc, &opt_geo, c2, &x_q[..dim]);
            let (n1_vals, _) = shape(&lam1, opt_order);
            let (n2_vals, _) = shape(&lam2, opt_order);

            let mut e1_side1 = [c(0.0, 0.0); 3];
            let mut e2_side1 = [c(0.0, 0.0); 3];
            let mut e1_side2 = [c(0.0, 0.0); 3];
            let mut e2_side2 = [c(0.0, 0.0); 3];
            for a in 0..nloc_opt {
                let s1 = n1_vals[a];
                let s2 = n2_vals[a];
                let v1_1 = e1_cell[c1 * nloc_opt + a];
                let v2_1 = e2_cell[c1 * nloc_opt + a];
                let v1_2 = e1_cell[c2 * nloc_opt + a];
                let v2_2 = e2_cell[c2 * nloc_opt + a];
                for k in 0..3 {
                    e1_side1[k] += v1_1[k] * s1;
                    e2_side1[k] += v2_1[k] * s1;
                    e1_side2[k] += v1_2[k] * s2;
                    e2_side2[k] += v2_2[k] * s2;
                }
            }

            // Continuous tangential electric field E_par (average of side 1 and side 2 projections)
            let project_par = |e: [C; 3]| -> [C; 3] {
                let en = e[0] * n_wafer[0] + e[1] * n_wafer[1] + e[2] * n_wafer[2];
                [e[0] - en * n_wafer[0], e[1] - en * n_wafer[1], e[2] - en * n_wafer[2]]
            };
            let ep1_a = project_par(e1_side1);
            let ep1_b = project_par(e1_side2);
            let ep2_a = project_par(e2_side1);
            let ep2_b = project_par(e2_side2);
            let e1_par = [(ep1_a[0] + ep1_b[0]) * 0.5, (ep1_a[1] + ep1_b[1]) * 0.5, (ep1_a[2] + ep1_b[2]) * 0.5];
            let e2_par = [(ep2_a[0] + ep2_b[0]) * 0.5, (ep2_a[1] + ep2_b[1]) * 0.5, (ep2_a[2] + ep2_b[2]) * 0.5];

            // Continuous normal displacement D_perp = eps0 * (eps_r E) . n_hat
            // By the two-sided Dyson/Johnson variational derivation, -(eps1^-1 - eps2^-1) (D_perp^(1))^* D_perp^(2)
            // equals (eps1 - eps2) (E_perp^(1))^* E_perp^(2), eliminating O(eps1^2) one-sided FE jump bias.
            let normal_d = |eps: &M3, e: [C; 3]| -> C {
                let mut dn = c(0.0, 0.0);
                for a in 0..3 {
                    for b in 0..3 {
                        dn += e[b] * (EPS0 * eps[a][b] * n_wafer[a]);
                    }
                }
                dn
            };
            let d1_s1 = normal_d(eps1, e1_side1);
            let d1_s2 = normal_d(eps2, e1_side2);
            let d2_s1 = normal_d(eps1, e2_side1);
            let d2_s2 = normal_d(eps2, e2_side2);

            // (eps_1 - eps_2) |E_par|^2 - (eps_1^-1 - eps_2^-1) |D_perp|^2
            let mut delta_eps_epar_sq = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    delta_eps_epar_sq += EPS0 * (eps1[a][b] - eps2[a][b]) * (e1_par[a].conj() * e2_par[b]).re;
                }
            }
            let d_perp_sq = 0.5 * ((d1_s1.conj() * d2_s2).re + (d1_s2.conj() * d2_s1).re);
            let surface_density = delta_eps_epar_sq - inv_eps1_minus_inv_eps2 * d_perp_sq;

            mb_integral += measure * wf * u_dot_n * surface_density;
        }
        if facet_touched {
            interface_facets += 1;
        }
    }

    let g_om_mb = -omega_opt * mb_integral / (2.0 * u_opt);
    let g_om_tot = g_om_pe + g_om_mb;

    let g_pe_rad_s = g_om_pe * x_zpf;
    let g_mb_rad_s = g_om_mb * x_zpf;
    let g_0_rad_s = g_pe_rad_s + g_mb_rad_s;

    let g_pe_hz = g_pe_rad_s / (2.0 * PI);
    let g_mb_hz = g_mb_rad_s / (2.0 * PI);
    let g_0_hz = g_0_rad_s / (2.0 * PI);

    let q_mech = infer_q_mech(mech_result, opts.mechanical_mode_index, opts);
    let vg1 = infer_v_g(opt_result, opt_sol, opts.optical_mode_index, opts);
    let vg2 = infer_v_g(opt_result, stokes_sol, stokes_idx, opts).or(vg1);
    let sbs_gain_w_inv_m_inv = match (vg1, vg2) {
        (Some(v1), Some(v2)) if dim == 2 && v1 > 0.0 && v2 > 0.0 => {
            // G_B = (8 Q_m |G_OM|^2 L_y) / (omega_opt m_eff Omega_m^2 v_{g,1} v_{g,2})
            //     = (16 Q_m |g_0|^2 L_y) / (hbar omega_opt Omega_m v_{g,1} v_{g,2})
            let gb = (16.0 * q_mech * g_0_rad_s * g_0_rad_s * ly) / (HBAR * omega_opt * omega_mech * v1 * v2);
            Some(gb)
        }
        _ => None,
    };

    Ok(OptomechanicsReport {
        optical_mode: opts.optical_mode_index + 1,
        mechanical_mode: opts.mechanical_mode_index + 1,
        optical_freq_hz: f_opt,
        mechanical_freq_hz: f_mech,
        m_eff_kg: m_eff,
        x_zpf_m: x_zpf,
        u_opt_j: u_opt,
        g_om_pe_rad_s_per_m: g_om_pe,
        g_om_mb_rad_s_per_m: g_om_mb,
        g_om_rad_s_per_m: g_om_tot,
        g_om_hz_per_m: g_om_tot / (2.0 * PI),
        g_pe_rad_s,
        g_pe_hz,
        g_mb_rad_s,
        g_mb_hz,
        g_0_rad_s,
        g_0_hz,
        g_0_abs_rad_s: g_0_rad_s.abs(),
        g_0_abs_hz: g_0_hz.abs(),
        q_mech,
        sbs_gain_w_inv_m_inv,
        active_core_cells,
        interface_facets,
        cross_mesh: is_cross_mesh,
    })
}

/// True when [`OptomechanicsAnalysis`] can run on `result`.
pub fn optomechanics_applies(result: &SolveResult) -> bool {
    matches!(result.physics.as_str(), "ems" | "sms") && !result.solutions.is_empty()
}

/// The cross-mesh cavity optomechanics and Brillouin coupling analysis behind [`Analysis`].
pub struct OptomechanicsAnalysis;

impl Analysis for OptomechanicsAnalysis {
    fn id(&self) -> &'static str {
        "optomechanics"
    }

    fn name(&self) -> &'static str {
        "Optomechanical & Brillouin coupling (g_PE, g_MB, g_0, G_B)"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        optomechanics_applies(result)
    }

    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        self.run_with_params(mesh, result, &Value::Null)
    }

    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        let opts = OptomechanicsOptions::from_params(params)?;
        if let Some(mech_res_val) = params.get("mechanical_result").filter(|v| !v.is_null()) {
            let mech_res: SolveResult = serde_json::from_value(mech_res_val.clone())
                .map_err(|e| Error::invalid(format!("invalid 'mechanical_result' in params: {e}")))?;
            return Ok(optomechanics_with_meshes(mesh, result, mesh, &mech_res, &opts)?.to_json());
        }
        if let Some(opt_res_val) = params.get("optical_result").filter(|v| !v.is_null()) {
            let opt_res: SolveResult = serde_json::from_value(opt_res_val.clone())
                .map_err(|e| Error::invalid(format!("invalid 'optical_result' in params: {e}")))?;
            return Ok(optomechanics_with_meshes(mesh, &opt_res, mesh, result, &opts)?.to_json());
        }
        Err(Error::invalid(
            "optomechanics analysis requires both optical ('ems') and mechanical ('sms') results (run via Engine::analyse_with or pass 'mechanical_result' / 'optical_result' in params)",
        ))
    }
}
