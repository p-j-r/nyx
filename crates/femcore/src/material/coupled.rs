// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Coupled multi-physics material tensors and superconducting interface loss models (`format = 2`):
//! - `[dielectric]`: clamped (`eps_clamped`), free (`eps_free`), high-frequency (`eps_inf`),
//!   anisotropic loss tangent (`loss_tangent`) and saturable Two-Level System (`tls_tan_delta_0`, `tls_e_crit`) loss.
//! - `[interfaces.{ma,ms,sa}]`: native surface oxide TLS dielectric loss layers.
//! - `[piezoelectric]`: stress-charge `e_iJ` (`C/m^2`) and strain-charge `d_iJ` (`C/N = m/V`) tensors
//!   plus exact Legendre dual conversions (`e = d C^E`, `C^D = C^E + e^T (eps^S)^-1 e`, `eps^T = eps^S + d e^T`).
//! - `[photoelastic]`: general asymmetric, indefinite Pockels photoelastic `6x6` tensor `p_IJ`.
//! - `[electro_optic]`: clamped (`r^S`) and free (`r^T`) linear electro-optic Pockels `6x3` tensors.
//! - `[viscoelastic]`: symmetric positive-semidefinite phonon viscosity `6x6` tensor `eta_IJ` (`Pa*s`).
//!
//! All tensors stored inside this module use textbook Voigt index order `(xx, yy, zz, yz, xz, xy)`.
//! Conversion helpers [`solver_to_voigt_m6`] and [`voigt_to_solver_m6`] map cleanly to and from
//! `femcore::material`'s solver Voigt order `(xx, yy, zz, xy, yz, xz)`.

use super::elastic::sym_eigenvalues;
use super::reader::{MaterialError, Section, Tree};
use super::source::{read_sources, write_sources, Source};
use super::state_t::{ComponentClamp, PointGroup, TempTensor};
use super::tensor::{Requirement, TensorProp};
use super::units::{self, Dim};
use super::{rotate_tensor, M3, M6, Stiffness};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Vacuum electric permittivity `epsilon_0` in F/m (CODATA 2022).
pub const EPSILON_0: f64 = 8.854_187_818_8e-12;
const KB: f64 = 1.380_649e-23;
const PLANCK_H: f64 = 6.626_070_15e-34;

const VP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];
/// Index in solver Voigt `(xx, yy, zz, xy, yz, xz)` of each textbook Voigt entry `(xx, yy, zz, yz, xz, xy)`.
const TO_SOLVER: [usize; 6] = [0, 1, 2, 4, 5, 3];

/// Converts a `6x6` matrix from solver Voigt `(xx, yy, zz, xy, yz, xz)` to textbook Voigt `(xx, yy, zz, yz, xz, xy)`.
pub fn solver_to_voigt_m6(m: &M6) -> M6 {
    let mut out = [[0.0; 6]; 6];
    for (i, out_row) in out.iter_mut().enumerate() {
        let si = TO_SOLVER[i];
        for (j, slot) in out_row.iter_mut().enumerate() {
            *slot = m[si][TO_SOLVER[j]];
        }
    }
    out
}

/// Converts a `6x6` matrix from textbook Voigt `(xx, yy, zz, yz, xz, xy)` to solver Voigt `(xx, yy, zz, xy, yz, xz)`.
pub fn voigt_to_solver_m6(m: &M6) -> M6 {
    let mut out = [[0.0; 6]; 6];
    for (i, row) in m.iter().enumerate() {
        let si = TO_SOLVER[i];
        for (j, &val) in row.iter().enumerate() {
            out[si][TO_SOLVER[j]] = val;
        }
    }
    out
}

/// Exact `3x3` matrix inverse.
pub fn inv3(a: &M3) -> Result<M3, String> {
    let c00 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    let c01 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    let c02 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    let det = a[0][0] * c00 + a[0][1] * c01 + a[0][2] * c02;
    let scale = a.iter().flatten().fold(0.0_f64, |acc, &x| acc.max(x.abs())).powi(3);
    if !det.is_finite() || det.abs() <= 1e-14 * scale.max(1e-300) {
        return Err("3x3 tensor is singular and cannot be inverted".into());
    }
    let inv = 1.0 / det;
    Ok([
        [c00 * inv, (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv, (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv],
        [c01 * inv, (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv, (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv],
        [c02 * inv, (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv, (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv],
    ])
}

/// Partial-pivoted Gauss-Jordan `6x6` matrix inverse.
pub fn inv6(a: &M6) -> Result<M6, String> {
    let mut m = *a;
    let mut inv = [[0.0; 6]; 6];
    for (i, row) in inv.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let scale = a.iter().flatten().fold(0.0_f64, |acc, &x| acc.max(x.abs())).max(1e-300);
    for col in 0..6 {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for (r, row) in m.iter().enumerate().skip(col + 1) {
            if row[col].abs() > best {
                best = row[col].abs();
                piv = r;
            }
        }
        if !best.is_finite() || best <= 1e-13 * scale {
            return Err("6x6 stiffness matrix is singular and cannot be inverted".into());
        }
        if piv != col {
            m.swap(piv, col);
            inv.swap(piv, col);
        }
        let d = m[col][col];
        for j in 0..6 {
            m[col][j] /= d;
            inv[col][j] /= d;
        }
        for r in 0..6 {
            if r == col {
                continue;
            }
            let f = m[r][col];
            if f != 0.0 {
                for j in 0..6 {
                    m[r][j] -= f * m[col][j];
                    inv[r][j] -= f * inv[col][j];
                }
            }
        }
    }
    Ok(inv)
}

fn symmetrize3(a: &M3) -> M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = 0.5 * (a[i][j] + a[j][i]);
        }
    }
    out
}

fn mul_mat3(a: &M3, b: &M3) -> M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

/// TLS loss tangent factor `tan_delta_tls_0 * tanh(h f / (2 k_B T)) / sqrt(1 + (E / E_c)^2)`.
pub fn tls_saturation_factor(
    tan_delta_tls_0: f64,
    e_crit: Option<f64>,
    freq_hz: Option<f64>,
    temperature_k: Option<f64>,
    e_field_v_per_m: Option<f64>,
) -> f64 {
    if tan_delta_tls_0 <= 0.0 || freq_hz.is_some_and(|f| f <= 0.0) {
        return 0.0;
    }
    let thermal = match (freq_hz, temperature_k.filter(|&t| t > 0.0)) {
        (Some(f), Some(t)) => (PLANCK_H * f / (2.0 * KB * t)).tanh(),
        _ => 1.0,
    };
    let power = match (e_crit.filter(|&ec| ec > 0.0), e_field_v_per_m.filter(|&e| e >= 0.0)) {
        (Some(ec), Some(e)) => 1.0 / (1.0 + (e / ec).powi(2)).sqrt(),
        _ => 1.0,
    };
    tan_delta_tls_0 * thermal * power
}

/// Microwave/RF dielectric permittivity and bulk + TLS dielectric loss (`[dielectric]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Dielectric {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    /// Low-frequency / microwave relative permittivity alias when clamped vs free split is omitted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_r: Option<TensorProp>,
    /// Clamped (constant-strain) static/microwave relative permittivity `eps^S / eps_0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_clamped: Option<TensorProp>,
    /// Free (constant-stress) static/microwave relative permittivity `eps^T / eps_0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_free: Option<TensorProp>,
    /// High-frequency clamped electronic relative permittivity `eps_inf` (above phonons).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_inf: Option<TensorProp>,
    /// Temperature-dependent relative permittivity curve `eps_r(T)`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_t: Option<TempTensor>,
    /// Bulk microwave/RF dielectric loss tangent `tan(delta) >= 0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loss_tangent: Option<TensorProp>,
    /// Optional reference frequency in Hz for `loss_tangent`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loss_ref_freq: Option<f64>,
    /// Zero-field / zero-temperature TLS loss tangent `tan(delta_TLS,0) >= 0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tls_tan_delta_0: Option<f64>,
    /// Critical electric field `E_c` in V/m for TLS saturation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tls_e_crit: Option<f64>,
    /// Optional electric gyrotropy / Voigt axial vector $\mathbf{g}_{\text{re}}$ in the crystal frame.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gyrotropy: Option<[f64; 3]>,
    /// Optional imaginary part $\mathbf{g}_{\text{im}}$ of the electric gyrotropy axial vector.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gyrotropy_im: Option<[f64; 3]>,
    /// Optional relative permeability tensor $\boldsymbol{\mu}_r$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mu_r: Option<TensorProp>,
    /// Optional magnetic gyrotropy / Polder axial vector $\mathbf{g}_{\mu,\text{re}}$ in the crystal frame.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mu_gyrotropy: Option<[f64; 3]>,
    /// Optional imaginary part $\mathbf{g}_{\mu,\text{im}}$ of the magnetic gyrotropy axial vector.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mu_gyrotropy_im: Option<[f64; 3]>,
}

pub const DIELECTRIC_KEYS: &[&str] = &[
    "source",
    "eps_r",
    "eps_clamped",
    "eps_free",
    "eps_inf",
    "eps_T",
    "loss_tangent",
    "loss_ref_freq",
    "loss_ref_freq_unit",
    "tls_tan_delta_0",
    "tls_e_crit",
    "tls_e_crit_unit",
];

const DIELECTRIC_PARSE_KEYS: &[&str] = &[
    "source",
    "eps_r",
    "eps_clamped",
    "eps_free",
    "eps_inf",
    "eps_T",
    "loss_tangent",
    "loss_ref_freq",
    "loss_ref_freq_unit",
    "tls_tan_delta_0",
    "tls_e_crit",
    "tls_e_crit_unit",
    "gyrotropy",
    "gyrotropy_im",
    "mu_r",
    "mu_gyrotropy",
    "mu_gyrotropy_im",
];

impl Dielectric {
    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let eps_r = TensorProp::read(s, "eps_r", Dim::Dimensionless, Requirement::PositiveDefinite, "dielectric.eps_r")?;
        let eps_clamped = TensorProp::read(
            s,
            "eps_clamped",
            Dim::Dimensionless,
            Requirement::PositiveDefinite,
            "dielectric.eps_clamped",
        )?;
        let eps_free = TensorProp::read(
            s,
            "eps_free",
            Dim::Dimensionless,
            Requirement::PositiveDefinite,
            "dielectric.eps_free",
        )?;
        let eps_inf = TensorProp::read(
            s,
            "eps_inf",
            Dim::Dimensionless,
            Requirement::PositiveDefinite,
            "dielectric.eps_inf",
        )?;
        let eps_t = TempTensor::read(s, "eps_T", Dim::Dimensionless, None)?;
        let loss_tangent = TensorProp::read(
            s,
            "loss_tangent",
            Dim::Dimensionless,
            Requirement::NonNegative,
            "dielectric.loss_tangent",
        )?;
        let loss_ref_freq = s.quantity("loss_ref_freq", Dim::FrequencyHz, None)?;
        if let Some(f) = loss_ref_freq {
            if !f.is_finite() || f <= 0.0 {
                return Err(s.err("loss_ref_freq", "loss_ref_freq must be positive"));
            }
        }
        let tls_tan_delta_0 = s.quantity("tls_tan_delta_0", Dim::Dimensionless, None)?;
        if let Some(td) = tls_tan_delta_0 {
            if !td.is_finite() || td < 0.0 {
                return Err(s.err("tls_tan_delta_0", "tls_tan_delta_0 must be non-negative"));
            }
        }
        let tls_e_crit = s.quantity("tls_e_crit", Dim::ElectricField, None)?;
        if let Some(ec) = tls_e_crit {
            if !ec.is_finite() || ec <= 0.0 {
                return Err(s.err("tls_e_crit", "tls_e_crit must be positive"));
            }
        }
        let gyrotropy = super::tensor::read_axial_vec3(s, "gyrotropy")?;
        let gyrotropy_im = super::tensor::read_axial_vec3(s, "gyrotropy_im")?;
        let mu_r = TensorProp::read(s, "mu_r", Dim::Dimensionless, Requirement::PositiveDefinite, "dielectric.mu_r")?;
        let mu_gyrotropy = super::tensor::read_axial_vec3(s, "mu_gyrotropy")?;
        let mu_gyrotropy_im = super::tensor::read_axial_vec3(s, "mu_gyrotropy_im")?;
        s.finish(DIELECTRIC_PARSE_KEYS)?;
        Ok(Self {
            source,
            eps_r,
            eps_clamped,
            eps_free,
            eps_inf,
            eps_t,
            loss_tangent,
            loss_ref_freq,
            tls_tan_delta_0,
            tls_e_crit,
            gyrotropy,
            gyrotropy_im,
            mu_r,
            mu_gyrotropy,
            mu_gyrotropy_im,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(t) = &self.eps_r {
            tbl.insert("eps_r".into(), t.write(1.0));
        }
        if let Some(t) = &self.eps_clamped {
            tbl.insert("eps_clamped".into(), t.write(1.0));
        }
        if let Some(t) = &self.eps_free {
            tbl.insert("eps_free".into(), t.write(1.0));
        }
        if let Some(t) = &self.eps_inf {
            tbl.insert("eps_inf".into(), t.write(1.0));
        }
        if let Some(t) = &self.eps_t {
            tbl.insert("eps_T".into(), t.to_toml(None));
        }
        if let Some(t) = &self.loss_tangent {
            tbl.insert("loss_tangent".into(), t.write(1.0));
        }
        if let Some(f) = self.loss_ref_freq {
            tbl.insert("loss_ref_freq".into(), format!("{f} Hz").into());
        }
        if let Some(td) = self.tls_tan_delta_0 {
            tbl.insert("tls_tan_delta_0".into(), toml::Value::Float(td));
        }
        if let Some(ec) = self.tls_e_crit {
            tbl.insert("tls_e_crit".into(), format!("{ec} V/m").into());
        }
        if let Some(g) = self.gyrotropy {
            tbl.insert("gyrotropy".into(), toml::Value::Array(g.into_iter().map(toml::Value::Float).collect()));
        }
        if let Some(g) = self.gyrotropy_im {
            tbl.insert("gyrotropy_im".into(), toml::Value::Array(g.into_iter().map(toml::Value::Float).collect()));
        }
        if let Some(t) = &self.mu_r {
            tbl.insert("mu_r".into(), t.write(1.0));
        }
        if let Some(g) = self.mu_gyrotropy {
            tbl.insert("mu_gyrotropy".into(), toml::Value::Array(g.into_iter().map(toml::Value::Float).collect()));
        }
        if let Some(g) = self.mu_gyrotropy_im {
            tbl.insert("mu_gyrotropy_im".into(), toml::Value::Array(g.into_iter().map(toml::Value::Float).collect()));
        }
    }

    /// Evaluates `(eps_re, eps_im_loss)` together with any `eps_T` clamping events in the crystal frame.
    pub fn complex_permittivity_checked(
        &self,
        freq_hz: Option<f64>,
        temperature_k: Option<f64>,
        e_field_v_per_m: Option<f64>,
        prefer_free: bool,
    ) -> Option<((M3, M3), Vec<ComponentClamp>)> {
        let mut clamps = Vec::new();
        let base_re = if let (Some(t), Some(curve)) = (temperature_k, &self.eps_t) {
            let (tp, c) = curve.eval_tensor_checked(t);
            clamps = c;
            Some(tp.matrix())
        } else {
            let picked = if prefer_free {
                self.eps_free.as_ref().or(self.eps_clamped.as_ref()).or(self.eps_r.as_ref())
            } else {
                self.eps_clamped.as_ref().or(self.eps_free.as_ref()).or(self.eps_r.as_ref())
            };
            picked.map(TensorProp::matrix)
        }?;

        let mut tan_d = self.loss_tangent.as_ref().map(TensorProp::matrix).unwrap_or([[0.0; 3]; 3]);
        if let Some(tls0) = self.tls_tan_delta_0 {
            let td_tls = tls_saturation_factor(tls0, self.tls_e_crit, freq_hz, temperature_k, e_field_v_per_m);
            for (i, row) in tan_d.iter_mut().enumerate() {
                row[i] += td_tls;
            }
        }
        let loss = symmetrize3(&mul_mat3(&base_re, &tan_d));
        Some(((base_re, loss), clamps))
    }

    /// Evaluates `(eps_re, eps_im_loss)` in the crystal frame where passive loss corresponds to
    /// `eps_tilde = eps_re - i eps_im_loss` (`exp(+i omega t)` convention) and `eps_im_loss >= 0`.
    pub fn complex_permittivity(
        &self,
        freq_hz: Option<f64>,
        temperature_k: Option<f64>,
        e_field_v_per_m: Option<f64>,
        prefer_free: bool,
    ) -> Option<(M3, M3)> {
        self.complex_permittivity_checked(freq_hz, temperature_k, e_field_v_per_m, prefer_free).map(|(v, _)| v)
    }
}

/// Native oxide or interface Two-Level System loss layer (`[interfaces.{ma,ms,sa}]`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InterfaceLayer {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oxide: Option<String>,
    /// Native oxide / interface layer thickness in m.
    pub thickness: f64,
    /// Relative dielectric permittivity of the interface oxide.
    pub eps_r: f64,
    /// Intrinsic TLS loss tangent `tan(delta_TLS,0) >= 0`.
    pub tan_delta_tls_0: f64,
    /// Critical electric field `E_c` in V/m for TLS saturation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub e_crit: Option<f64>,
}

pub const IFACE_KEYS: &[&str] = &[
    "source",
    "oxide",
    "thickness",
    "thickness_unit",
    "eps_r",
    "tan_delta_tls_0",
    "e_crit",
    "e_crit_unit",
];

impl InterfaceLayer {
    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let oxide = s.string("oxide")?;
        let thickness = s
            .quantity("thickness", Dim::Length, None)?
            .ok_or_else(|| s.err("thickness", "interface layer requires thickness, e.g. \"2 nm\""))?;
        if !thickness.is_finite() || thickness <= 0.0 {
            return Err(s.err("thickness", "thickness must be positive"));
        }
        let eps_r = s
            .quantity("eps_r", Dim::Dimensionless, None)?
            .ok_or_else(|| s.err("eps_r", "interface layer requires eps_r > 0"))?;
        if !eps_r.is_finite() || eps_r <= 0.0 {
            return Err(s.err("eps_r", "eps_r must be positive"));
        }
        let tan_delta_tls_0 = s
            .quantity("tan_delta_tls_0", Dim::Dimensionless, None)?
            .ok_or_else(|| s.err("tan_delta_tls_0", "interface layer requires tan_delta_tls_0 >= 0"))?;
        if !tan_delta_tls_0.is_finite() || tan_delta_tls_0 < 0.0 {
            return Err(s.err("tan_delta_tls_0", "tan_delta_tls_0 must be non-negative"));
        }
        let e_crit = s.quantity("e_crit", Dim::ElectricField, None)?;
        if let Some(ec) = e_crit {
            if !ec.is_finite() || ec <= 0.0 {
                return Err(s.err("e_crit", "e_crit must be positive"));
            }
        }
        s.finish(IFACE_KEYS)?;
        Ok(Self {
            source,
            oxide,
            thickness,
            eps_r,
            tan_delta_tls_0,
            e_crit,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(ox) = &self.oxide {
            tbl.insert("oxide".into(), ox.clone().into());
        }
        tbl.insert("thickness".into(), format!("{} m", self.thickness).into());
        tbl.insert("eps_r".into(), toml::Value::Float(self.eps_r));
        tbl.insert("tan_delta_tls_0".into(), toml::Value::Float(self.tan_delta_tls_0));
        if let Some(ec) = self.e_crit {
            tbl.insert("e_crit".into(), format!("{ec} V/m").into());
        }
    }

    pub fn effective_tan_delta(&self, freq_hz: Option<f64>, temperature_k: Option<f64>, e_field_v_per_m: Option<f64>) -> f64 {
        tls_saturation_factor(self.tan_delta_tls_0, self.e_crit, freq_hz, temperature_k, e_field_v_per_m)
    }
}

/// Superconducting / dielectric surface oxide interfaces (`[interfaces]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Interfaces {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ma: Option<InterfaceLayer>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ms: Option<InterfaceLayer>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sa: Option<InterfaceLayer>,
}

impl Interfaces {
    pub fn is_empty(&self) -> bool {
        self.ma.is_none() && self.ms.is_none() && self.sa.is_none()
    }

    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let ma = s.sub("ma")?.map(|mut sub| InterfaceLayer::read(&mut sub)).transpose()?;
        let ms = s.sub("ms")?.map(|mut sub| InterfaceLayer::read(&mut sub)).transpose()?;
        let sa = s.sub("sa")?.map(|mut sub| InterfaceLayer::read(&mut sub)).transpose()?;
        s.finish(&["ma", "ms", "sa"])?;
        Ok(Self { ma, ms, sa })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        for (name, layer) in [("ma", &self.ma), ("ms", &self.ms), ("sa", &self.sa)] {
            if let Some(l) = layer {
                let mut sub = toml::Table::new();
                l.write(&mut sub);
                tbl.insert(name.into(), toml::Value::Table(sub));
            }
        }
    }
}

pub const INTERFACES_KEYS: &[&str] = &["ma", "ms", "sa"];
pub const PIEZOELECTRIC_KEYS: &[&str] = &["source", "form", "unit", "matrix"];
pub const PHOTOELASTIC_KEYS: &[&str] = &["source", "wavelength", "wavelength_unit", "matrix"];
pub const ELECTRO_OPTIC_KEYS: &[&str] = &["source", "wavelength", "wavelength_unit", "unit", "clamped", "free", "matrix", "state"];
pub const VISCOELASTIC_KEYS: &[&str] = &["source", "unit", "matrix", "cryo_qm_floor"];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PiezoForm {
    StressCharge,
    StrainCharge,
}

/// Piezoelectric coupling tensor (`[piezoelectric]`), stored in textbook Voigt `(xx, yy, zz, yz, xz, xy)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Piezoelectric {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    pub form: PiezoForm,
    /// `3x6` coupling matrix in SI (`C/m^2` for `StressCharge`, `C/N = m/V` for `StrainCharge`).
    pub matrix: [[f64; 6]; 3],
}

fn expand_piezo_named(
    form: PiezoForm,
    pg: Option<PointGroup>,
    raw: &[[Option<f64>; 6]; 3],
) -> [[f64; 6]; 3] {
    let mut m = [[0.0; 6]; 3];
    for i in 0..3 {
        for j in 0..6 {
            if let Some(v) = raw[i][j] {
                m[i][j] = v;
            }
        }
    }
    match pg {
        Some(PointGroup::C6v | PointGroup::C4v | PointGroup::C6 | PointGroup::C4) => {
            // Hexagonal/tetragonal 6mm / 4mm: e31 = e32, e15 = e24; optional chiral e14 = -e25
            if raw[2][1].is_none() && raw[2][0].is_some() {
                m[2][1] = m[2][0];
            }
            if raw[1][3].is_none() && raw[0][4].is_some() {
                m[1][3] = m[0][4];
            }
            if matches!(pg, Some(PointGroup::C6 | PointGroup::C4)) && raw[1][4].is_none() && raw[0][3].is_some() {
                m[1][4] = -m[0][3];
            }
        }
        Some(PointGroup::C3v) => {
            // Trigonal 3m (IEEE-1987, mirror perp X): e21 = -e22, e16 = -e22 (or d16 = -2*d22), e24 = e15, e32 = e31
            if raw[1][0].is_none() && raw[1][1].is_some() {
                m[1][0] = -m[1][1];
            }
            if raw[0][5].is_none() && raw[1][1].is_some() {
                let f = if form == PiezoForm::StrainCharge { -2.0 } else { -1.0 };
                m[0][5] = f * m[1][1];
            }
            if raw[1][3].is_none() && raw[0][4].is_some() {
                m[1][3] = m[0][4];
            }
            if raw[2][1].is_none() && raw[2][0].is_some() {
                m[2][1] = m[2][0];
            }
        }
        Some(PointGroup::D3) => {
            // Trigonal 32 (alpha-quartz): d12 = -d11, d26 = -2*d11 (e26 = -e11), d25 = -d14
            if raw[0][1].is_none() && raw[0][0].is_some() {
                m[0][1] = -m[0][0];
            }
            if raw[1][5].is_none() && raw[0][0].is_some() {
                let f = if form == PiezoForm::StrainCharge { -2.0 } else { -1.0 };
                m[1][5] = f * m[0][0];
            }
            if raw[1][4].is_none() && raw[0][3].is_some() {
                m[1][4] = -m[0][3];
            }
        }
        Some(PointGroup::Td | PointGroup::T23) => {
            // Cubic -43m / 23: e14 = e25 = e36
            if let Some(v) = raw[0][3] {
                if raw[1][4].is_none() {
                    m[1][4] = v;
                }
                if raw[2][5].is_none() {
                    m[2][5] = v;
                }
            }
        }
        Some(PointGroup::D2d) if raw[1][4].is_none() && raw[0][3].is_some() => {
            // Tetragonal -42m: e14 = e25, independent e36
            m[1][4] = m[0][3];
        }
        _ => {}
    }
    m
}

impl Piezoelectric {
    /// Expands `[piezoelectric]` textbook Voigt `3x6` matrix to a full rank-3 tensor `t[i][j][k]`.
    pub fn to_rank3(&self) -> [[[f64; 3]; 3]; 3] {
        let mut t = [[[0.0; 3]; 3]; 3];
        for i in 0..3 {
            for (voigt_j, &(a, b)) in VP.iter().enumerate() {
                let scale = if self.form == PiezoForm::StrainCharge && voigt_j >= 3 {
                    0.5
                } else {
                    1.0
                };
                let v = self.matrix[i][voigt_j] * scale;
                t[i][a][b] = v;
                t[i][b][a] = v;
            }
        }
        t
    }

    /// Reconstructs a `3x6` textbook Voigt matrix from a rank-3 tensor `t[i][j][k]` (`t[i][j][k] == t[i][k][j]`).
    pub fn from_rank3(form: PiezoForm, t: &[[[f64; 3]; 3]; 3]) -> [[f64; 6]; 3] {
        let mut m = [[0.0; 6]; 3];
        for i in 0..3 {
            for (voigt_j, &(a, b)) in VP.iter().enumerate() {
                let scale = if form == PiezoForm::StrainCharge && voigt_j >= 3 {
                    2.0
                } else {
                    1.0
                };
                m[i][voigt_j] = t[i][a][b] * scale;
            }
        }
        m
    }

    /// Rotates the `3x6` piezoelectric coupling matrix from crystal frame into global frame via `R`.
    pub fn rotate(&self, r: &M3) -> [[f64; 6]; 3] {
        let t = self.to_rank3();
        let mut out = [[[0.0; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let mut s = 0.0;
                    for a in 0..3 {
                        for b in 0..3 {
                            for c in 0..3 {
                                s += r[i][a] * r[j][b] * r[k][c] * t[a][b][c];
                            }
                        }
                    }
                    out[i][j][k] = s;
                }
            }
        }
        Self::from_rank3(self.form, &out)
    }

    /// Converts to stress-charge `e_iJ` (`C/m^2`) in textbook Voigt, multiplying `d * C^E` if stored as strain-charge.
    pub fn stress_charge(&self, stiffness_voigt: Option<&M6>) -> Result<[[f64; 6]; 3], String> {
        match self.form {
            PiezoForm::StressCharge => Ok(self.matrix),
            PiezoForm::StrainCharge => {
                let ce = stiffness_voigt.ok_or_else(|| {
                    "converting strain-charge [piezoelectric] (d) to stress-charge (e) requires [elastic] stiffness C^E".to_string()
                })?;
                let mut e = [[0.0; 6]; 3];
                for i in 0..3 {
                    for j in 0..6 {
                        let mut sum = 0.0;
                        for k in 0..6 {
                            sum += self.matrix[i][k] * ce[k][j];
                        }
                        e[i][j] = sum;
                    }
                }
                Ok(e)
            }
        }
    }

    /// Converts to strain-charge `d_iJ` (`C/N = m/V`) in textbook Voigt, multiplying `e * (C^E)^-1` if stored as stress-charge.
    pub fn strain_charge(&self, stiffness_voigt: Option<&M6>) -> Result<[[f64; 6]; 3], String> {
        match self.form {
            PiezoForm::StrainCharge => Ok(self.matrix),
            PiezoForm::StressCharge => {
                let ce = stiffness_voigt.ok_or_else(|| {
                    "converting stress-charge [piezoelectric] (e) to strain-charge (d) requires [elastic] stiffness C^E".to_string()
                })?;
                let s_e = inv6(ce)?;
                let mut d = [[0.0; 6]; 3];
                for i in 0..3 {
                    for j in 0..6 {
                        let mut sum = 0.0;
                        for k in 0..6 {
                            sum += self.matrix[i][k] * s_e[k][j];
                        }
                        d[i][j] = sum;
                    }
                }
                Ok(d)
            }
        }
    }

    pub fn read(s: &mut Section, pg: Option<PointGroup>) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let form_str = s.string("form")?;
        let unit_str = s.string("unit")?;

        let form = match form_str.as_deref() {
            Some("stress_charge") | Some("e") => Some(PiezoForm::StressCharge),
            Some("strain_charge") | Some("d") => Some(PiezoForm::StrainCharge),
            Some(other) => {
                return Err(s.err(
                    "form",
                    format!("unknown piezoelectric form '{other}'; expected \"stress_charge\" or \"strain_charge\""),
                ))
            }
            None => {
                // Infer from named keys eij / dij if present
                let has_d = (1..=3).any(|i| (1..=6).any(|j| s.has(&format!("d{i}{j}"))));
                let has_e = (1..=3).any(|i| (1..=6).any(|j| s.has(&format!("e{i}{j}"))));
                if has_e && !has_d {
                    Some(PiezoForm::StressCharge)
                } else if has_d && !has_e {
                    Some(PiezoForm::StrainCharge)
                } else if let Some(u) = unit_str.as_deref() {
                    if units::to_si(1.0, u, Dim::PiezoelectricStress).is_ok() {
                        Some(PiezoForm::StressCharge)
                    } else if units::to_si(1.0, u, Dim::PiezoelectricStrain).is_ok() {
                        Some(PiezoForm::StrainCharge)
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
        }
        .ok_or_else(|| s.err("form", "specify form = \"stress_charge\" (e) or \"strain_charge\" (d)"))?;

        let dim = match form {
            PiezoForm::StressCharge => Dim::PiezoelectricStress,
            PiezoForm::StrainCharge => Dim::PiezoelectricStrain,
        };

        let mut known: Vec<String> = vec!["source".into(), "form".into(), "unit".into(), "matrix".into()];
        let prefix = if form == PiezoForm::StressCharge { 'e' } else { 'd' };
        for i in 1..=3 {
            for j in 1..=6 {
                known.push(format!("{prefix}{i}{j}"));
            }
        }
        let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();

        let matrix = if s.has("matrix") {
            let tree = s
                .quantities("matrix", dim, unit_str.as_deref())?
                .ok_or_else(|| s.err("matrix", "must be a 3x6 array"))?;
            read_3x6_tree(s, "matrix", &tree)?
        } else {
            let mut raw = [[None; 6]; 3];
            let mut any = false;
            for i in 0..3 {
                for j in 0..6 {
                    let key = format!("{prefix}{}{}", i + 1, j + 1);
                    if let Some(v) = s.quantity(&key, dim, unit_str.as_deref())? {
                        raw[i][j] = Some(v);
                        any = true;
                    }
                }
            }
            if !any {
                return Err(s.err("", format!("write matrix = [...] or named constants like {prefix}31, {prefix}33, {prefix}15")));
            }
            expand_piezo_named(form, pg, &raw)
        };

        s.finish(&known_refs)?;
        Ok(Self { source, form, matrix })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        let (form_s, unit_s) = match self.form {
            PiezoForm::StressCharge => ("stress_charge", "C/m^2"),
            PiezoForm::StrainCharge => ("strain_charge", "C/N"),
        };
        tbl.insert("form".into(), form_s.into());
        tbl.insert("unit".into(), unit_s.into());
        let rows = self
            .matrix
            .iter()
            .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
            .collect();
        tbl.insert("matrix".into(), toml::Value::Array(rows));
    }
}

/// Photoelastic Pockels tensor `p_IJ` (`[photoelastic]`) in textbook Voigt `(xx, yy, zz, yz, xz, xy)`.
/// Note: photoelastic matrices can be asymmetric (`p13 != p31`, `p14 != p41`) and indefinite (`p < 0` allowed).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Photoelastic {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    /// Optional measurement wavelength in m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength: Option<f64>,
    /// `6x6` dimensionless Pockels matrix in textbook Voigt order.
    pub matrix: M6,
}

impl Photoelastic {
    /// Expands textbook Voigt `6x6` matrix `p_IJ` directly into full rank-4 tensor `p_{abcd}` (`p_{abcd} = p_{IJ}` without factor 2).
    pub fn to_rank4(&self) -> [[[[f64; 3]; 3]; 3]; 3] {
        let mut t = [[[[0.0; 3]; 3]; 3]; 3];
        for (i, &(a, b)) in VP.iter().enumerate() {
            for (j, &(c, d)) in VP.iter().enumerate() {
                let v = self.matrix[i][j];
                for (p, q) in [(a, b), (b, a)] {
                    for (r, s) in [(c, d), (d, c)] {
                        t[p][q][r][s] = v;
                    }
                }
            }
        }
        t
    }

    /// Rotates photoelastic tensor `p_{IJ}` into the global frame using `R`, returning textbook Voigt `6x6`.
    pub fn rotate_voigt(&self, r: &M3) -> M6 {
        let t = self.to_rank4();
        let mut rot = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut sum = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                for c in 0..3 {
                                    for d in 0..3 {
                                        sum += r[i][a] * r[j][b] * r[k][c] * r[l][d] * t[a][b][c][d];
                                    }
                                }
                            }
                        }
                        rot[i][j][k][l] = sum;
                    }
                }
            }
        }
        let mut out = [[0.0; 6]; 6];
        for (i, &(a, b)) in VP.iter().enumerate() {
            for (j, &(c, d)) in VP.iter().enumerate() {
                out[i][j] = rot[a][b][c][d];
            }
        }
        out
    }

    /// Rotates photoelastic tensor `p_{IJ}` into solver Voigt order `(xx, yy, zz, xy, yz, xz)`.
    pub fn rotate_solver(&self, r: &M3) -> M6 {
        voigt_to_solver_m6(&self.rotate_voigt(r))
    }

    pub fn read(s: &mut Section, pg: Option<PointGroup>) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let wavelength = s.quantity("wavelength", Dim::Length, None)?;
        if let Some(w) = wavelength {
            if !w.is_finite() || w <= 0.0 {
                return Err(s.err("wavelength", "wavelength must be positive"));
            }
        }

        let mut known: Vec<String> = vec!["source".into(), "wavelength".into(), "wavelength_unit".into(), "matrix".into()];
        for i in 1..=6 {
            for j in 1..=6 {
                known.push(format!("p{i}{j}"));
            }
        }
        let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();

        let matrix = if s.has("matrix") {
            let tree = s
                .quantities("matrix", Dim::Dimensionless, None)?
                .ok_or_else(|| s.err("matrix", "must be a 6x6 array"))?;
            read_6x6_tree(s, "matrix", &tree)?
        } else {
            let mut raw = [[None; 6]; 6];
            let mut any = false;
            for i in 0..6 {
                for j in 0..6 {
                    let key = format!("p{}{}", i + 1, j + 1);
                    if let Some(v) = s.quantity(&key, Dim::Dimensionless, None)? {
                        raw[i][j] = Some(v);
                        any = true;
                    }
                }
            }
            if !any {
                return Err(s.err("", "write matrix = [...] or named constants like p11, p12, p44"));
            }
            expand_photoelastic_named(pg, &raw)
        };

        s.finish(&known_refs)?;
        Ok(Self {
            source,
            wavelength,
            matrix,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(w) = self.wavelength {
            tbl.insert("wavelength".into(), format!("{w} m").into());
        }
        let rows = self
            .matrix
            .iter()
            .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
            .collect();
        tbl.insert("matrix".into(), toml::Value::Array(rows));
    }
}

fn expand_photoelastic_named(pg: Option<PointGroup>, raw: &[[Option<f64>; 6]; 6]) -> M6 {
    let mut m = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            if let Some(v) = raw[i][j] {
                m[i][j] = v;
            }
        }
    }
    match pg {
        Some(PointGroup::Isotropic) | None => {
            let p11 = raw[0][0].unwrap_or(0.0);
            let p12 = raw[0][1].unwrap_or(0.0);
            let p44 = raw[3][3].unwrap_or(0.5 * (p11 - p12));
            if raw[0][0].is_some() || raw[0][1].is_some() {
                for i in 0..3 {
                    m[i][i] = raw[i][i].unwrap_or(p11);
                    for j in 0..3 {
                        if i != j {
                            m[i][j] = raw[i][j].unwrap_or(p12);
                        }
                    }
                    m[i + 3][i + 3] = raw[i + 3][i + 3].unwrap_or(p44);
                }
            }
        }
        Some(PointGroup::Oh | PointGroup::Td | PointGroup::O432) => {
            let p11 = raw[0][0].unwrap_or(0.0);
            let p12 = raw[0][1].unwrap_or(0.0);
            let p44 = raw[3][3].unwrap_or(0.0);
            for i in 0..3 {
                m[i][i] = raw[i][i].unwrap_or(p11);
                for j in 0..3 {
                    if i != j {
                        m[i][j] = raw[i][j].unwrap_or(p12);
                    }
                }
                m[i + 3][i + 3] = raw[i + 3][i + 3].unwrap_or(p44);
            }
        }
        Some(
            PointGroup::C6v
            | PointGroup::D6
            | PointGroup::D3h
            | PointGroup::D6h
            | PointGroup::TransverseIsotropic
            | PointGroup::C4v
            | PointGroup::D4
            | PointGroup::D2d
            | PointGroup::D4h,
        ) => {
            if raw[1][1].is_none() && raw[0][0].is_some() {
                m[1][1] = m[0][0];
            }
            if raw[1][0].is_none() && raw[0][1].is_some() {
                m[1][0] = m[0][1];
            }
            if raw[1][2].is_none() && raw[0][2].is_some() {
                m[1][2] = m[0][2];
            }
            if raw[2][1].is_none() && raw[2][0].is_some() {
                m[2][1] = m[2][0];
            }
            if raw[4][4].is_none() && raw[3][3].is_some() {
                m[4][4] = m[3][3];
            }
            if raw[5][5].is_none()
                && matches!(
                    pg,
                    Some(PointGroup::C6v | PointGroup::D6 | PointGroup::D3h | PointGroup::D6h | PointGroup::TransverseIsotropic)
                )
            {
                m[5][5] = 0.5 * (m[0][0] - m[0][1]);
            }
        }
        Some(PointGroup::C3v | PointGroup::D3d | PointGroup::D3) => {
            // Trigonal 3m / -3m / 32 (IEEE-1987, asymmetric): p13 != p31, p14 != p41;
            // p24 = -p14, p65 = p14, p42 = -p41, p56 = p41, p66 = 0.5*(p11 - p12).
            if raw[1][1].is_none() && raw[0][0].is_some() {
                m[1][1] = m[0][0];
            }
            if raw[1][0].is_none() && raw[0][1].is_some() {
                m[1][0] = m[0][1];
            }
            if raw[1][2].is_none() && raw[0][2].is_some() {
                m[1][2] = m[0][2];
            }
            if raw[2][1].is_none() && raw[2][0].is_some() {
                m[2][1] = m[2][0];
            }
            if raw[4][4].is_none() && raw[3][3].is_some() {
                m[4][4] = m[3][3];
            }
            if raw[1][3].is_none() && raw[0][3].is_some() {
                m[1][3] = -m[0][3];
            }
            if raw[5][4].is_none() && raw[0][3].is_some() {
                m[5][4] = m[0][3];
            }
            if raw[3][1].is_none() && raw[3][0].is_some() {
                m[3][1] = -m[3][0];
            }
            if raw[4][5].is_none() && raw[3][0].is_some() {
                m[4][5] = m[3][0];
            }
            if raw[5][5].is_none() {
                m[5][5] = 0.5 * (m[0][0] - m[0][1]);
            }
        }
        _ => {}
    }
    m
}

/// Linear electro-optic (Pockels) tensor `r_Ij` (`[electro_optic]`) in m/V and textbook Voigt order.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ElectroOptic {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength: Option<f64>,
    /// Clamped (constant strain) `6x3` Pockels matrix `r^S_Ij` in m/V.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clamped: Option<[[f64; 3]; 6]>,
    /// Free (constant stress) `6x3` Pockels matrix `r^T_Ij` in m/V.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub free: Option<[[f64; 3]; 6]>,
}

fn expand_eo_named(pg: Option<PointGroup>, raw: &[[Option<f64>; 3]; 6]) -> [[f64; 3]; 6] {
    let mut m = [[0.0; 3]; 6];
    for i in 0..6 {
        for j in 0..3 {
            if let Some(v) = raw[i][j] {
                m[i][j] = v;
            }
        }
    }
    match pg {
        Some(PointGroup::C6v | PointGroup::C4v) => {
            if raw[1][2].is_none() && raw[0][2].is_some() {
                m[1][2] = m[0][2];
            }
            if raw[3][1].is_none() && raw[4][0].is_some() {
                m[3][1] = m[4][0];
            }
        }
        Some(PointGroup::C3v) => {
            // Trigonal 3m: r13 = r23, r22 = -r12 = -r61, r51 = r42
            if raw[1][2].is_none() && raw[0][2].is_some() {
                m[1][2] = m[0][2];
            }
            if raw[0][1].is_none() && raw[1][1].is_some() {
                m[0][1] = -m[1][1];
            }
            if raw[5][0].is_none() && raw[1][1].is_some() {
                m[5][0] = -m[1][1];
            }
            if raw[3][1].is_none() && raw[4][0].is_some() {
                m[3][1] = m[4][0];
            }
        }
        Some(PointGroup::Td | PointGroup::T23) => {
            if let Some(v) = raw[3][0] {
                if raw[4][1].is_none() {
                    m[4][1] = v;
                }
                if raw[5][2].is_none() {
                    m[5][2] = v;
                }
            }
        }
        _ => {}
    }
    m
}

impl ElectroOptic {
    /// Derives free `r^T` from clamped `r^S` (`r^T = r^S + p * d^T`), or returns explicit `free`.
    pub fn resolved_free(&self, pe: Option<&Photoelastic>, d_strain_charge: Option<&[[f64; 6]; 3]>) -> Option<[[f64; 3]; 6]> {
        if let Some(f) = self.free {
            return Some(f);
        }
        let s = self.clamped?;
        let (Some(p), Some(d)) = (pe, d_strain_charge) else {
            return Some(s);
        };
        let mut out = s;
        for i in 0..6 {
            for j in 0..3 {
                let mut piezo_opt = 0.0;
                for k in 0..6 {
                    piezo_opt += p.matrix[i][k] * d[j][k];
                }
                out[i][j] += piezo_opt;
            }
        }
        Some(out)
    }

    /// Derives clamped `r^S` from free `r^T` (`r^S = r^T - p * d^T`), or returns explicit `clamped`.
    pub fn resolved_clamped(&self, pe: Option<&Photoelastic>, d_strain_charge: Option<&[[f64; 6]; 3]>) -> Option<[[f64; 3]; 6]> {
        if let Some(c) = self.clamped {
            return Some(c);
        }
        let f = self.free?;
        let (Some(p), Some(d)) = (pe, d_strain_charge) else {
            return Some(f);
        };
        let mut out = f;
        for i in 0..6 {
            for j in 0..3 {
                let mut piezo_opt = 0.0;
                for k in 0..6 {
                    piezo_opt += p.matrix[i][k] * d[j][k];
                }
                out[i][j] -= piezo_opt;
            }
        }
        Some(out)
    }

    /// Expands a 6x3 textbook-Voigt electro-optic matrix `r_Ik` (`(xx, yy, zz, yz, xz, xy)`) to rank-3 tensor `r_abk`.
    pub fn to_rank3(r_voigt: &[[f64; 3]; 6]) -> [[[f64; 3]; 3]; 3] {
        let mut t = [[[0.0; 3]; 3]; 3];
        for (i, &(a, b)) in VP.iter().enumerate() {
            for k in 0..3 {
                let v = r_voigt[i][k];
                t[a][b][k] = v;
                t[b][a][k] = v;
            }
        }
        t
    }

    /// Rotates a 6x3 textbook-Voigt electro-optic matrix `r_Ik` from the crystal frame to the global frame
    /// via `rot` (`x_global = rot * x_crystal`), returning a 6x3 matrix in textbook Voigt order `(xx, yy, zz, yz, xz, xy)`.
    pub fn rotate_voigt(r_voigt: &[[f64; 3]; 6], rot: &M3) -> [[f64; 3]; 6] {
        let t = Self::to_rank3(r_voigt);
        let mut out = [[0.0; 3]; 6];
        for (ii, &(i, j)) in VP.iter().enumerate() {
            for k in 0..3 {
                let mut sum = 0.0;
                for a in 0..3 {
                    let ria = rot[i][a];
                    if ria == 0.0 {
                        continue;
                    }
                    for b in 0..3 {
                        let rjb = rot[j][b];
                        if rjb == 0.0 {
                            continue;
                        }
                        let ria_rjb = ria * rjb;
                        for c in 0..3 {
                            sum += ria_rjb * rot[k][c] * t[a][b][c];
                        }
                    }
                }
                out[ii][k] = sum;
            }
        }
        out
    }
}

/// Rotates a 6x3 textbook-Voigt electro-optic tensor `r_voigt` by Bunge Z-X-Z `euler_deg` into the global wafer frame.
pub fn electro_optic_global(r_voigt: &[[f64; 3]; 6], euler_deg: [f64; 3]) -> [[f64; 3]; 6] {
    if euler_deg == [0.0; 3] {
        return *r_voigt;
    }
    let rot = super::rotation(euler_deg);
    ElectroOptic::rotate_voigt(r_voigt, &rot)
}

/// Rotates a 6x6 textbook-Voigt photoelastic tensor `p_textbook` by Bunge Z-X-Z `euler_deg` into the global wafer frame.
pub fn photoelastic_global(p_textbook: &M6, euler_deg: [f64; 3]) -> M6 {
    if euler_deg == [0.0; 3] {
        return *p_textbook;
    }
    let rot = super::rotation(euler_deg);
    let pe = Photoelastic { source: Vec::new(), wavelength: None, matrix: *p_textbook };
    pe.rotate_voigt(&rot)
}

/// Evaluates the exact relative permittivity perturbation tensor:
///
/// `\Delta(\varepsilon_r^{-1})_I = \sum_{k=0}^2 r^S_{Ik} E_{\text{DC},k} + \sum_{J=0}^5 p_{IJ} S^{\text{tb}}_J`,
/// `\Delta \boldsymbol{\varepsilon}_r = -\boldsymbol{\varepsilon}_r \, \Delta(\boldsymbol{\varepsilon}_r^{-1}) \, \boldsymbol{\varepsilon}_r`,
///
/// where `r_global_tb` (`6x3`, m/V) and `p_global_tb` (`6x6`, dimensionless) are in the global wafer frame
/// in **textbook Voigt order** `(xx, yy, zz, yz, xz, xy)`, `e_dc_wafer` is `[Ex, Ey, Ez]` (V/m) in wafer axes,
/// and `strain_wafer_solver` is `[Sxx, Syy, Szz, Sxy, Syz, Sxz]` (m/m) in **solver Voigt order**.
pub fn delta_eps_eo_pe(
    eps_r: &M3,
    r_global_tb: Option<&[[f64; 3]; 6]>,
    p_global_tb: Option<&M6>,
    e_dc_wafer: Option<&[f64; 3]>,
    strain_wafer_solver: Option<&[f64; 6]>,
) -> M3 {
    let mut d_eta_v = [0.0; 6];
    let mut active = false;
    if let (Some(r), Some(e)) = (r_global_tb, e_dc_wafer) {
        for i in 0..6 {
            d_eta_v[i] += r[i][0] * e[0] + r[i][1] * e[1] + r[i][2] * e[2];
        }
        active = true;
    }
    if let (Some(p), Some(s_sol)) = (p_global_tb, strain_wafer_solver) {
        let s_tb = [s_sol[0], s_sol[1], s_sol[2], s_sol[4], s_sol[5], s_sol[3]];
        for i in 0..6 {
            let mut acc = 0.0;
            for j in 0..6 {
                acc += p[i][j] * s_tb[j];
            }
            d_eta_v[i] += acc;
        }
        active = true;
    }
    if !active {
        return [[0.0; 3]; 3];
    }
    let d_eta = [
        [d_eta_v[0], d_eta_v[5], d_eta_v[4]],
        [d_eta_v[5], d_eta_v[1], d_eta_v[3]],
        [d_eta_v[4], d_eta_v[3], d_eta_v[2]],
    ];
    let mut tmp = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            tmp[i][j] = eps_r[i][0] * d_eta[0][j] + eps_r[i][1] * d_eta[1][j] + eps_r[i][2] * d_eta[2][j];
        }
    }
    let mut d_eps = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            d_eps[i][j] = -(tmp[i][0] * eps_r[0][j] + tmp[i][1] * eps_r[1][j] + tmp[i][2] * eps_r[2][j]);
        }
    }
    for i in 0..3 {
        for j in (i + 1)..3 {
            let s = 0.5 * (d_eps[i][j] + d_eps[j][i]);
            d_eps[i][j] = s;
            d_eps[j][i] = s;
        }
    }
    d_eps
}

impl ElectroOptic {
    pub fn read(s: &mut Section, pg: Option<PointGroup>) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let wavelength = s.quantity("wavelength", Dim::Length, None)?;
        if let Some(w) = wavelength {
            if !w.is_finite() || w <= 0.0 {
                return Err(s.err("wavelength", "wavelength must be positive"));
            }
        }
        let unit_str = s.string("unit")?;

        let read_sub_6x3 = |sub: &mut Section| -> Result<[[f64; 3]; 6], MaterialError> {
            let sub_unit = sub.string("unit")?.or_else(|| unit_str.clone());
            let mut known: Vec<String> = vec!["unit".into(), "matrix".into()];
            for i in 1..=6 {
                for j in 1..=3 {
                    known.push(format!("r{i}{j}"));
                }
            }
            let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();
            let out = if sub.has("matrix") {
                let tree = sub
                    .quantities("matrix", Dim::ElectroOptic, sub_unit.as_deref())?
                    .ok_or_else(|| sub.err("matrix", "must be a 6x3 array"))?;
                read_6x3_tree(sub, "matrix", &tree)?
            } else {
                let mut raw = [[None; 3]; 6];
                let mut any = false;
                for i in 0..6 {
                    for j in 0..3 {
                        let key = format!("r{}{}", i + 1, j + 1);
                        if let Some(v) = sub.quantity(&key, Dim::ElectroOptic, sub_unit.as_deref())? {
                            raw[i][j] = Some(v);
                            any = true;
                        }
                    }
                }
                if !any {
                    return Err(sub.err("", "write matrix = [...] or named Pockels constants like r33, r13, r51"));
                }
                expand_eo_named(pg, &raw)
            };
            sub.finish(&known_refs)?;
            Ok(out)
        };

        let clamped = s.sub("clamped")?.map(|mut sub| read_sub_6x3(&mut sub)).transpose()?;
        let free = s.sub("free")?.map(|mut sub| read_sub_6x3(&mut sub)).transpose()?;

        let mut known: Vec<String> = vec![
            "source".into(),
            "wavelength".into(),
            "wavelength_unit".into(),
            "unit".into(),
            "clamped".into(),
            "free".into(),
            "matrix".into(),
            "state".into(),
        ];
        for i in 1..=6 {
            for j in 1..=3 {
                known.push(format!("r{i}{j}"));
            }
        }
        let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();

        let (clamped, free) = if clamped.is_some() || free.is_some() {
            let _ = s.string("state")?;
            (clamped, free)
        } else {
            let state = s.string("state")?.unwrap_or_else(|| "clamped".to_string());
            let m = if s.has("matrix") {
                let tree = s
                    .quantities("matrix", Dim::ElectroOptic, unit_str.as_deref())?
                    .ok_or_else(|| s.err("matrix", "must be a 6x3 array"))?;
                read_6x3_tree(s, "matrix", &tree)?
            } else {
                let mut raw = [[None; 3]; 6];
                let mut any = false;
                for i in 0..6 {
                    for j in 0..3 {
                        let key = format!("r{}{}", i + 1, j + 1);
                        if let Some(v) = s.quantity(&key, Dim::ElectroOptic, unit_str.as_deref())? {
                            raw[i][j] = Some(v);
                            any = true;
                        }
                    }
                }
                if !any {
                    return Err(s.err("", "specify [electro_optic.clamped], [electro_optic.free], or r_Ij constants"));
                }
                expand_eo_named(pg, &raw)
            };
            match state.as_str() {
                "clamped" | "S" => (Some(m), None),
                "free" | "T" => (None, Some(m)),
                other => return Err(s.err("state", format!("unknown state '{other}'; expected \"clamped\" or \"free\""))),
            }
        };

        s.finish(&known_refs)?;
        Ok(Self {
            source,
            wavelength,
            clamped,
            free,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(w) = self.wavelength {
            tbl.insert("wavelength".into(), format!("{w} m").into());
        }
        tbl.insert("unit".into(), "m/V".into());
        for (name, mat) in [("clamped", &self.clamped), ("free", &self.free)] {
            if let Some(m) = mat {
                let mut sub = toml::Table::new();
                let rows = m
                    .iter()
                    .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
                    .collect();
                sub.insert("matrix".into(), toml::Value::Array(rows));
                tbl.insert(name.into(), toml::Value::Table(sub));
            }
        }
    }
}

/// Phononic / viscoelastic damping tensor `eta_IJ` (`[viscoelastic]`) in `Pa*s` and textbook Voigt order.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Viscoelastic {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    /// Symmetric positive-semidefinite viscosity matrix in `Pa*s` (textbook Voigt).
    pub matrix: M6,
    /// Optional cryogenic mechanical quality factor floor `Q_m`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cryo_qm_floor: Option<f64>,
}

impl Viscoelastic {
    pub fn to_solver_voigt(&self) -> M6 {
        voigt_to_solver_m6(&self.matrix)
    }

    pub fn read(s: &mut Section, pg: Option<PointGroup>) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let unit_str = s.string("unit")?;
        let cryo_qm_floor = s.quantity("cryo_qm_floor", Dim::Dimensionless, None)?;
        if let Some(q) = cryo_qm_floor {
            if !q.is_finite() || q <= 0.0 {
                return Err(s.err("cryo_qm_floor", "cryo_qm_floor must be positive"));
            }
        }

        let mut known: Vec<String> = vec!["source".into(), "unit".into(), "matrix".into(), "cryo_qm_floor".into()];
        for i in 1..=6 {
            for j in i..=6 {
                known.push(format!("eta{i}{j}"));
            }
        }
        let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();

        let matrix = if s.has("matrix") {
            let tree = s
                .quantities("matrix", Dim::DynamicViscosity, unit_str.as_deref())?
                .ok_or_else(|| s.err("matrix", "must be a symmetric 6x6 array"))?;
            let m = read_6x6_tree(s, "matrix", &tree)?;
            symmetrize_checked_6x6(s, "matrix", &m)?
        } else {
            let mut raw = [[None; 6]; 6];
            let mut any = false;
            for (i, row) in raw.iter_mut().enumerate() {
                for (j, slot) in row.iter_mut().enumerate().skip(i) {
                    let key = format!("eta{}{}", i + 1, j + 1);
                    if let Some(v) = s.quantity(&key, Dim::DynamicViscosity, unit_str.as_deref())? {
                        *slot = Some(v);
                        any = true;
                    }
                }
            }
            if !any {
                return Err(s.err("", "write matrix = [...] or named constants like eta11, eta12, eta44"));
            }
            expand_viscous_named(pg, &raw)
        };

        // Passivity / positive-semidefiniteness check
        let evals = sym_eigenvalues(matrix);
        let scale = matrix.iter().flatten().fold(0.0_f64, |a, &v| a.max(v.abs())).max(f64::MIN_POSITIVE);
        if evals[0] < -1e-9 * scale {
            return Err(s.err(
                "",
                format!(
                    "viscosity tensor must be positive semidefinite (passive), but smallest eigenvalue is {}",
                    evals[0]
                ),
            ));
        }

        s.finish(&known_refs)?;
        Ok(Self {
            source,
            matrix,
            cryo_qm_floor,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        tbl.insert("unit".into(), "Pa*s".into());
        let rows = self
            .matrix
            .iter()
            .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
            .collect();
        tbl.insert("matrix".into(), toml::Value::Array(rows));
        if let Some(q) = self.cryo_qm_floor {
            tbl.insert("cryo_qm_floor".into(), toml::Value::Float(q));
        }
    }
}

fn expand_viscous_named(pg: Option<PointGroup>, raw: &[[Option<f64>; 6]; 6]) -> M6 {
    let mut m = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in i..6 {
            if let Some(v) = raw[i][j] {
                m[i][j] = v;
                m[j][i] = v;
            }
        }
    }
    let sym = matches!(
        pg,
        Some(PointGroup::Isotropic) | Some(PointGroup::Oh | PointGroup::Td | PointGroup::O432) | None
    );
    if sym {
        let e11 = raw[0][0].unwrap_or(0.0);
        let e12 = raw[0][1].unwrap_or(0.0);
        let e44 = raw[3][3].unwrap_or(if matches!(pg, Some(PointGroup::Isotropic) | None) {
            0.5 * (e11 - e12)
        } else {
            0.0
        });
        for i in 0..3 {
            m[i][i] = raw[i][i].unwrap_or(e11);
            for j in (i + 1)..3 {
                let v = raw[i][j].unwrap_or(e12);
                m[i][j] = v;
                m[j][i] = v;
            }
            m[i + 3][i + 3] = raw[i + 3][i + 3].unwrap_or(e44);
        }
    } else if matches!(
        pg,
        Some(PointGroup::C6v | PointGroup::D6 | PointGroup::D6h | PointGroup::TransverseIsotropic | PointGroup::C3v | PointGroup::D3d | PointGroup::D3)
    ) {
        if raw[1][1].is_none() && raw[0][0].is_some() {
            m[1][1] = m[0][0];
        }
        if raw[1][2].is_none() && raw[0][2].is_some() {
            m[1][2] = m[0][2];
            m[2][1] = m[0][2];
        }
        if raw[4][4].is_none() && raw[3][3].is_some() {
            m[4][4] = m[3][3];
        }
        if raw[5][5].is_none() {
            m[5][5] = 0.5 * (m[0][0] - m[0][1]);
        }
    }
    m
}

/// Exact piezoelectric Legendre duals in textbook Voigt order `(xx, yy, zz, yz, xz, xy)`.
/// Given clamped `C^E` (`Pa`), piezoelectric `e` or `d`, and either `eps_clamped` or `eps_free` (relative),
/// returns `(eps_clamped_r, eps_free_r, C^E, C^D)`.
pub fn resolve_piezo_legendre_duals(
    ce_voigt: &M6,
    piezo: &Piezoelectric,
    eps_clamped_r: Option<&M3>,
    eps_free_r: Option<&M3>,
) -> Result<(M3, M3, M6, M6), String> {
    let e = piezo.stress_charge(Some(ce_voigt))?;
    let d = piezo.strain_charge(Some(ce_voigt))?;

    // Piezoelectric dielectric increment d * e^T / eps_0 (dimensionless relative permittivity jump)
    let mut de_t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut sum = 0.0;
            for k in 0..6 {
                sum += d[i][k] * e[j][k];
            }
            de_t[i][j] = sum / EPSILON_0;
        }
    }

    let (ec, ef) = match (eps_clamped_r, eps_free_r) {
        (Some(c), Some(f)) => (*c, *f),
        (Some(c), None) => {
            let mut f = *c;
            for i in 0..3 {
                for j in 0..3 {
                    f[i][j] += de_t[i][j];
                }
            }
            (*c, f)
        }
        (None, Some(f)) => {
            let mut c = *f;
            for i in 0..3 {
                for j in 0..3 {
                    c[i][j] -= de_t[i][j];
                }
            }
            // Sanity check positivity of derived clamped permittivity
            if sym_eigenvalues(c)[0] <= 0.0 {
                return Err("derived clamped permittivity eps^S = eps^T - d e^T is not positive definite; check tensor units".into());
            }
            (c, *f)
        }
        (None, None) => {
            return Err("Legendre dual requires at least eps_clamped or eps_free relative permittivity".into());
        }
    };

    // Open-circuit stiffness C^D = C^E + e^T (eps^S)^{-1} e
    let mut eps_s_si = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            eps_s_si[i][j] = ec[i][j] * EPSILON_0;
        }
    }
    let beta_s = inv3(&eps_s_si)?;
    let mut cd = *ce_voigt;
    for i in 0..6 {
        for j in 0..6 {
            let mut stiff_inc = 0.0;
            for p in 0..3 {
                for q in 0..3 {
                    stiff_inc += e[p][i] * beta_s[p][q] * e[q][j];
                }
            }
            cd[i][j] += stiff_inc;
        }
    }

    Ok((symmetrize3(&ec), symmetrize3(&ef), *ce_voigt, cd))
}

/// Rotates a stiffness `Stiffness` enum into textbook Voigt `(xx, yy, zz, yz, xz, xy)` without shear conversion.
pub fn stiffness_to_textbook_voigt(s: &Stiffness) -> M6 {
    let mut c = [[0.0; 6]; 6];
    let mut set = |i: usize, j: usize, v: f64| {
        c[i][j] = v;
        c[j][i] = v;
    };
    match *s {
        Stiffness::Cubic { c11, c12, c44 } => {
            for i in 0..3 {
                set(i, i, c11);
                set(i + 3, i + 3, c44);
                for j in i + 1..3 {
                    set(i, j, c12);
                }
            }
        }
        Stiffness::Hexagonal { c11, c12, c13, c33, c44 } => {
            set(0, 0, c11);
            set(1, 1, c11);
            set(2, 2, c33);
            set(0, 1, c12);
            set(0, 2, c13);
            set(1, 2, c13);
            set(3, 3, c44);
            set(4, 4, c44);
            set(5, 5, 0.5 * (c11 - c12));
        }
        Stiffness::Trigonal { c11, c12, c13, c14, c33, c44 } => {
            set(0, 0, c11);
            set(1, 1, c11);
            set(2, 2, c33);
            set(0, 1, c12);
            set(0, 2, c13);
            set(1, 2, c13);
            set(0, 3, c14);
            set(1, 3, -c14);
            set(4, 5, c14);
            set(3, 3, c44);
            set(4, 4, c44);
            set(5, 5, 0.5 * (c11 - c12));
        }
        Stiffness::Orthotropic {
            c11,
            c22,
            c33,
            c12,
            c13,
            c23,
            c44,
            c55,
            c66,
        } => {
            set(0, 0, c11);
            set(1, 1, c22);
            set(2, 2, c33);
            set(0, 1, c12);
            set(0, 2, c13);
            set(1, 2, c23);
            set(3, 3, c44);
            set(4, 4, c55);
            set(5, 5, c66);
        }
        Stiffness::General(m) => return m,
    }
    c
}

/// Rotates a `3x3` permittivity/conductivity/loss tensor by `R`.
pub fn rotate_m3(m: &M3, r: &M3) -> M3 {
    rotate_tensor(m, r)
}

/// Computes `R_z(angle_deg)` for crystallographic rotation checks.
pub fn rot_z_deg(angle_deg: f64) -> M3 {
    let (s, c) = (angle_deg * PI / 180.0).sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

fn tree_num(t: &Tree) -> Option<f64> {
    match t {
        Tree::Num(x) => Some(*x),
        _ => None,
    }
}

fn read_3x6_tree(s: &Section, key: &str, tree: &Tree) -> Result<[[f64; 6]; 3], MaterialError> {
    let Tree::List(rows) = tree else {
        return Err(s.err(key, "must be a 3x6 matrix of 3 rows of 6 numbers"));
    };
    if rows.len() != 3 {
        return Err(s.err(key, format!("expected 3 rows, got {}", rows.len())));
    }
    let mut out = [[0.0; 6]; 3];
    for (i, row) in rows.iter().enumerate() {
        let Tree::List(cols) = row else {
            return Err(s.err(key, "each row of the 3x6 matrix must be a list of 6 numbers"));
        };
        let vals: Option<Vec<f64>> = cols.iter().map(tree_num).collect();
        match vals {
            Some(v) if v.len() == 6 => out[i].copy_from_slice(&v),
            _ => return Err(s.err(key, format!("row {} must contain 6 finite numbers", i + 1))),
        }
    }
    Ok(out)
}

fn read_6x3_tree(s: &Section, key: &str, tree: &Tree) -> Result<[[f64; 3]; 6], MaterialError> {
    let Tree::List(rows) = tree else {
        return Err(s.err(key, "must be a 6x3 matrix of 6 rows of 3 numbers"));
    };
    if rows.len() != 6 {
        return Err(s.err(key, format!("expected 6 rows, got {}", rows.len())));
    }
    let mut out = [[0.0; 3]; 6];
    for (i, row) in rows.iter().enumerate() {
        let Tree::List(cols) = row else {
            return Err(s.err(key, "each row of the 6x3 matrix must be a list of 3 numbers"));
        };
        let vals: Option<Vec<f64>> = cols.iter().map(tree_num).collect();
        match vals {
            Some(v) if v.len() == 3 => out[i].copy_from_slice(&v),
            _ => return Err(s.err(key, format!("row {} must contain 3 finite numbers", i + 1))),
        }
    }
    Ok(out)
}

fn read_6x6_tree(s: &Section, key: &str, tree: &Tree) -> Result<M6, MaterialError> {
    let Tree::List(rows) = tree else {
        return Err(s.err(key, "must be a 6x6 matrix of 6 rows of 6 numbers"));
    };
    if rows.len() != 6 {
        return Err(s.err(key, format!("expected 6 rows, got {}", rows.len())));
    }
    let mut out = [[0.0; 6]; 6];
    for (i, row) in rows.iter().enumerate() {
        let Tree::List(cols) = row else {
            return Err(s.err(key, "each row of the 6x6 matrix must be a list of 6 numbers"));
        };
        let vals: Option<Vec<f64>> = cols.iter().map(tree_num).collect();
        match vals {
            Some(v) if v.len() == 6 => out[i].copy_from_slice(&v),
            _ => return Err(s.err(key, format!("row {} must contain 6 finite numbers", i + 1))),
        }
    }
    Ok(out)
}

fn symmetrize_checked_6x6(s: &Section, key: &str, m: &M6) -> Result<M6, MaterialError> {
    let sc = m.iter().flatten().fold(0.0_f64, |a, &v| a.max(v.abs())).max(f64::MIN_POSITIVE);
    let mut out = *m;
    for i in 0..6 {
        for j in (i + 1)..6 {
            if (m[i][j] - m[j][i]).abs() > 1e-9 * sc {
                return Err(s.err(key, format!("matrix entry ({},{}) = {} differs from ({},{}) = {}", i + 1, j + 1, m[i][j], j + 1, i + 1, m[j][i])));
            }
            let avg = 0.5 * (m[i][j] + m[j][i]);
            out[i][j] = avg;
            out[j][i] = avg;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_diff_6x6(a: &M6, b: &M6) -> f64 {
        let mut d = 0.0_f64;
        for i in 0..6 {
            for j in 0..6 {
                d = d.max((a[i][j] - b[i][j]).abs());
            }
        }
        d
    }

    fn max_diff_3x6(a: &[[f64; 6]; 3], b: &[[f64; 6]; 3]) -> f64 {
        let mut d = 0.0_f64;
        for i in 0..3 {
            for j in 0..6 {
                d = d.max((a[i][j] - b[i][j]).abs());
            }
        }
        d
    }

    #[test]
    fn voigt_roundtrip_is_exact_inverse() {
        let mut m = [[0.0; 6]; 6];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, slot) in row.iter_mut().enumerate() {
                *slot = (i * 10 + j + 1) as f64;
            }
        }
        assert_eq!(solver_to_voigt_m6(&voigt_to_solver_m6(&m)), m);
    }

    #[test]
    fn isotropic_photoelastic_is_rotation_invariant() {
        // Silica: p11 = 0.121, p12 = 0.270, p44 = (p11 - p12)/2 = -0.0745
        let toml_text = "p11 = 0.121\np12 = 0.270\n";
        let tbl: toml::Table = toml_text.parse().unwrap();
        let mut s = Section::new("silica", "photoelastic", &tbl);
        let pe = Photoelastic::read(&mut s, Some(PointGroup::Isotropic)).unwrap();
        assert!((pe.matrix[3][3] - (-0.0745)).abs() < 1e-12);
        let r = super::super::rotation([37.0, 52.0, -19.0]);
        let rot = pe.rotate_voigt(&r);
        assert!(max_diff_6x6(&rot, &pe.matrix) < 1e-13, "isotropic photoelastic not invariant");
    }

    #[test]
    fn trigonal_3m_linbo3_asymmetric_photoelastic_and_piezo_symmetry() {
        // LiNbO3 (trigonal 3m): Weis & Gaylord (1985) photoelastic & Smith & Welsh (1971) piezo + C^E
        let pe_toml = "p11 = -0.026\np12 = 0.090\np13 = 0.133\np31 = 0.179\np33 = 0.071\np14 = -0.075\np41 = 0.151\np44 = 0.146\n";
        let tbl: toml::Table = pe_toml.parse().unwrap();
        let mut s = Section::new("linbo3", "photoelastic", &tbl);
        let pe = Photoelastic::read(&mut s, Some(PointGroup::C3v)).unwrap();
        // Verify asymmetry preserved: p13 != p31 and p14 != p41
        assert!((pe.matrix[0][2] - 0.133).abs() < 1e-12);
        assert!((pe.matrix[2][0] - 0.179).abs() < 1e-12);
        assert!((pe.matrix[0][3] - (-0.075)).abs() < 1e-12);
        assert!((pe.matrix[3][0] - 0.151).abs() < 1e-12);

        // Must be strictly invariant under 120 deg C3z rotation, and NOT invariant under 60 deg!
        let rot120 = pe.rotate_voigt(&rot_z_deg(120.0));
        let rot60 = pe.rotate_voigt(&rot_z_deg(60.0));
        assert!(max_diff_6x6(&rot120, &pe.matrix) < 1e-12, "3m photoelastic failed 120 deg C3z invariance");
        assert!(max_diff_6x6(&rot60, &pe.matrix) > 0.05, "3m photoelastic should change p14/p41 under 60 deg");

        // Piezo strain-charge d_iJ with factor-of-2 shear rule d16 = -2 * d22:
        let d_toml = "unit = \"pC/N\"\nd15 = 69.2\nd22 = 20.8\nd31 = -0.85\nd33 = 6.0\n";
        let dtbl: toml::Table = d_toml.parse().unwrap();
        let mut ds = Section::new("linbo3", "piezoelectric", &dtbl);
        let d_pz = Piezoelectric::read(&mut ds, Some(PointGroup::C3v)).unwrap();
        assert!((d_pz.matrix[0][5] - (-41.6e-12)).abs() < 1e-22, "d16 = -2*d22");
        assert!(max_diff_3x6(&d_pz.rotate(&rot_z_deg(120.0)), &d_pz.matrix) < 1e-22);
        assert!(max_diff_3x6(&d_pz.rotate(&rot_z_deg(60.0)), &d_pz.matrix) > 10e-12);

        // LiNbO3 Legendre duals: Weis & Gaylord (1985) C^E [GPa] + clamped eps^S = (44, 44, 29)
        let ce = stiffness_to_textbook_voigt(&Stiffness::Trigonal {
            c11: 203e9,
            c12: 57.3e9,
            c13: 75.2e9,
            c14: 8.5e9,
            c33: 242.4e9,
            c44: 59.5e9,
        });
        let e_stress = d_pz.stress_charge(Some(&ce)).unwrap();
        // Exact e = d C^E values for Weis & Gaylord (1985) LiNbO3 inputs:
        // e22 = d22*(C11-C12) - d15*C14 = 2.44236, e15 = d15*C44 - 2*d22*C14 = 3.7638, e16 = -e22
        assert!((e_stress[0][5] + e_stress[1][1]).abs() < 1e-10, "e16 must equal -e22 for 3m");
        assert!((e_stress[1][1] - 2.44236).abs() < 1e-4, "LiNbO3 e22 = {}", e_stress[1][1]);
        assert!((e_stress[0][4] - 3.7638).abs() < 1e-4, "LiNbO3 e15 = {}", e_stress[0][4]);

        let eps_s = TensorProp::Uniaxial { perp: 44.3, par: 27.9 }.matrix();
        let (ec, ef, _, cd) = resolve_piezo_legendre_duals(&ce, &d_pz, Some(&eps_s), None).unwrap();
        assert!((ec[0][0] - 44.3).abs() < 1e-12);
        // Strong electromechanical coupling in LiNbO3 boosts eps^T_11 (~ 83..86) well above eps^S_11 (44.3)
        assert!(ef[0][0] > 78.0 && ef[0][0] < 90.0, "LiNbO3 eps^T_11 = {}", ef[0][0]);
        // Acoustic stiffening raises C44^D (~ 86..96 GPa) well above C44^E (59.5 GPa)
        assert!(cd[3][3] > 80e9 && cd[3][3] < 100e9, "LiNbO3 C44^D = {} GPa", cd[3][3] * 1e-9);
    }

    #[test]
    fn scaln_6mm_legendre_dual_c33d_boost() {
        // Sc_0.3 Al_0.7 N (wurtzite 6mm): strongly enhanced e33 = 2.45 C/m^2, C^E, clamped eps_r = (14.5, 14.5, 15.8)
        let pz_toml = "form = \"stress_charge\"\nunit = \"C/m^2\"\ne31 = -0.78\ne33 = 2.45\ne15 = -0.35\n";
        let tbl: toml::Table = pz_toml.parse().unwrap();
        let mut s = Section::new("scaln", "piezoelectric", &tbl);
        let pz = Piezoelectric::read(&mut s, Some(PointGroup::C6v)).unwrap();

        let ce = stiffness_to_textbook_voigt(&Stiffness::Hexagonal {
            c11: 318e9,
            c12: 141e9,
            c13: 122e9,
            c33: 262e9,
            c44: 105e9,
        });
        let eps_s = TensorProp::Uniaxial { perp: 14.5, par: 15.8 }.matrix();
        let (_, ef, _, cd) = resolve_piezo_legendre_duals(&ce, &pz, Some(&eps_s), None).unwrap();
        // C33^D = C33^E + e33^2 / (eps0 * eps_s_33) = 262 GPa + 2.45^2 / (8.854e-12 * 15.8) = 262 + 42.9 = 304.9 GPa
        assert!((cd[2][2] - 304.9e9).abs() < 1e9, "Sc0.3Al0.7N C33^D = {} GPa", cd[2][2] * 1e-9);
        assert!(ef[2][2] > 18.0, "Sc0.3Al0.7N eps^T_33 > eps^S_33: {}", ef[2][2]);
    }

    #[test]
    fn tls_dielectric_loss_saturates_with_temperature_and_power() {
        let dl = Dielectric {
            eps_clamped: Some(TensorProp::Scalar(11.45)),
            loss_tangent: Some(TensorProp::Scalar(1e-7)),
            tls_tan_delta_0: Some(3e-6),
            tls_e_crit: Some(10.0),
            ..Default::default()
        };
        // Low power, 5 GHz, 10 mK: h f / (2 k_B T) = 12 >> 1 -> full TLS loss (~3.1e-6 * 11.45)
        let (_, loss_mk_low) = dl.complex_permittivity(Some(5e9), Some(0.01), Some(0.01), false).unwrap();
        // High power E = 1000 V/m >> E_c = 10 V/m: TLS quenched by ~100x, leaving ~1e-7 bulk floor
        let (_, loss_mk_high) = dl.complex_permittivity(Some(5e9), Some(0.01), Some(1000.0), false).unwrap();
        // High temperature 2 K, low power: h f / (2 k_B T) = 0.06 -> tanh ~ 0.06
        let (_, loss_2k_low) = dl.complex_permittivity(Some(5e9), Some(2.0), Some(0.01), false).unwrap();

        assert!(loss_mk_low[0][0] > 10.0 * loss_mk_high[0][0]);
        assert!(loss_mk_low[0][0] > 10.0 * loss_2k_low[0][0]);
        assert!((loss_mk_high[0][0] / 11.45 - 1.3e-7).abs() < 2e-8);
    }

    #[test]
    fn tls_saturation_factor_is_zero_at_dc() {
        assert_eq!(tls_saturation_factor(3e-6, Some(10.0), Some(0.0), Some(0.01), Some(0.0)), 0.0);
        let dl = Dielectric {
            eps_clamped: Some(TensorProp::Scalar(11.45)),
            tls_tan_delta_0: Some(3e-6),
            tls_e_crit: Some(10.0),
            ..Default::default()
        };
        let (_, loss_dc) = dl.complex_permittivity(Some(0.0), Some(0.01), None, false).unwrap();
        assert_eq!(loss_dc, [[0.0; 3]; 3]);
    }
}
