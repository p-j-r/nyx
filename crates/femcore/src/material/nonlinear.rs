// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Second- and third-order optical nonlinearities (`[nonlinear]` TOML v2 section, ADR 143):
//!
//! 1. **Second-order susceptibility (`[nonlinear.chi2]`, `Chi2`)**:
//!    - Contracted $3 \times 6$ $d_{iJ} = \frac{1}{2}\chi^{(2)}_{iJ}$ tensor in $\text{m/V}$ (accepting `"pm/V"`),
//!      stored in textbook Voigt order `(1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy)` where
//!      $$P_i^{(2)} = 2\varepsilon_0 \sum_{j,k=0}^2 d_{ijk} E_j E_k
//!      = 2\varepsilon_0 \big(d_{i1} E_x^2 + d_{i2} E_y^2 + d_{i3} E_z^2
//!      + 2 d_{i4} E_y E_z + 2 d_{i5} E_x E_z + 2 d_{i6} E_x E_y\big).$$
//!    - Crystallographic `PointGroup` symmetry expansion (`3m`, `6mm`, `4mm`, `-43m`, `23`, `-42m`,
//!      `32`, `-6m2`, `mm2`) and Neumann centrosymmetry validation: centrosymmetric groups (plus chiral
//!      cubic `432`) reject non-zero $d_{iJ}$ unless `broken_centrosymmetry = true` is explicitly set
//!      (e.g. DC-field-induced or strain-poled $\chi^{(2)}$).
//!    - Exact rank-3 tensor expansion (`to_rank3`) and Bunge Z-X-Z Euler rotation (`d_global`).
//!    - Optional quasi-phase-matching (`[nonlinear.qpm]`, `QpmSpec`) supporting uniform domain poling
//!      sign (`+1` / `-1`) and spatial periodic poling (`Lambda_qpm`, `axis`, `duty_cycle`, `phase_offset`)
//!      evaluated at Gauss quadrature points via [`QpmSpec::poling_sign_at`].
//!
//! 2. **Third-order Kerr, TPA, and free-carrier nonlinearities (`[nonlinear.chi3]`, `Chi3`)**:
//!    - Isotropic Kerr refractive index `n_2` ($\text{m}^2/\text{W}$) or explicit tensor susceptibility
//!      components `chi3_xxxx`, `chi3_xxyy`, `chi3_xyxy`, `chi3_xyyx` ($\text{m}^2/\text{V}^2$) and
//!      cross-polarization ratio `b_over_a` ($B / A$, default `1.0` for electronic Kleinman symmetry).
//!    - Under `nyx-ems`'s $\tilde{\mathbf{E}}(\mathbf{r},t) = \operatorname{Re}[\mathbf{E}(\mathbf{r}) e^{+i\omega t}]$
//!      peak-phasor convention ($I = \frac{1}{2} c_0 \varepsilon_0 n_0 \|\mathbf{E}\|^2$), the monochromatic
//!      cross-polarization vector Kerr polarization is
//!      $\mathbf{P}^{(3,\text{Kerr})}(\omega) = \varepsilon_0 [ A (\mathbf{E}\cdot\mathbf{E}^*)\mathbf{E} + \frac{1}{2} B (\mathbf{E}\cdot\mathbf{E})\mathbf{E}^* ]$
//!      (or $\frac{3}{4}\varepsilon_0 \sum_{j,k,l} \chi^{(3),\text{glob}}_{ijkl} E_j E_k E_l^*$),
//!      with $A + \frac{1}{2}B = \frac{3}{4}\chi^{(3)}_{xxxx} = c_0 \varepsilon_0 n_0^2 n_2$ (for electronic
//!      Kleinman symmetry $A = B = \frac{2}{3} c_0 \varepsilon_0 n_0^2 n_2$, reproducing the Maker–Terhune
//!      $2/3$ circular-polarization ratio and cross-polarized SPM/XPM birefringence).
//!    - Two-photon absorption `beta_tpa` ($\text{m/W}$), free-carrier absorption cross-section
//!      `sigma_fca` ($\text{m}^2$), free-carrier dispersion `k_fcd` ($\text{m}^3$, $\Delta n_{\text{FCD}} = k_{\text{FCD}} N_c$),
//!      and carrier recombination lifetime `tau_c` ($\text{s}$), yielding steady-state TPA-generated
//!      carrier density $N_c(\mathbf{E}) = \frac{\tau_c \beta_{\text{TPA}}}{2\hbar\omega} I(\mathbf{E})^2$.
//!
//! 3. **Shared Time-Domain Evaluators**:
//!    - [`eval_chi2_polarization_time`], [`eval_chi3_kerr_time`], [`eval_chi3_envelope_time`], and
//!      unconditionally stable exponential-integrator [`drude_carrier_step`] so time-domain Maxwell
//!      solvers (`dgtd`, `fdtd`) share the exact same `[nonlinear]` TOML material definitions.

use super::bcs::HBAR;
use super::coupled::EPSILON_0;
use super::reader::{MaterialError, Section, Tree};
use super::source::{read_sources, write_sources, Source};
use super::state_t::{clamp_temp_range, CheckedTempEval, PointGroup};
use super::units::Dim;
use super::{rotation, M3};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Speed of light in vacuum $c_0$ ($\text{m/s}$).
pub const C0_NL: f64 = 299_792_458.0;

const VP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];

/// Allowed TOML keys in `[nonlinear]`.
pub const NONLINEAR_KEYS: &[&str] = &[
    "source",
    "wavelength",
    "wavelength_unit",
    "T_range",
    "T_unit",
    "broken_centrosymmetry",
    "chi2",
    "chi3",
    "qpm",
];

/// Allowed non-indexed TOML keys in `[nonlinear.chi2]`.
pub const CHI2_KEYS: &[&str] = &[
    "source",
    "wavelength",
    "wavelength_unit",
    "unit",
    "kleinman",
    "broken_centrosymmetry",
    "matrix",
    "qpm",
];

/// Allowed TOML keys in `[nonlinear.qpm]`.
pub const QPM_KEYS: &[&str] = &[
    "sign",
    "Lambda_qpm",
    "Lambda_qpm_unit",
    "period",
    "period_unit",
    "axis",
    "duty_cycle",
    "phase_offset",
    "phase_offset_unit",
];

/// Allowed TOML keys in `[nonlinear.chi3]`.
pub const CHI3_KEYS: &[&str] = &[
    "source",
    "wavelength",
    "wavelength_unit",
    "n_2",
    "n_2_unit",
    "unit",
    "chi3_xxxx",
    "chi3_xxyy",
    "chi3_xyxy",
    "chi3_xyyx",
    "b_over_a",
    "beta_tpa",
    "beta_tpa_unit",
    "sigma_fca",
    "sigma_fca_unit",
    "k_fcd",
    "k_fcd_unit",
    "tau_c",
    "tau_c_unit",
];

fn default_one() -> f64 {
    1.0
}

fn default_half() -> f64 {
    0.5
}

fn default_y_axis() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

fn is_zero_flag(b: &bool) -> bool {
    !*b
}

/// Quasi-phase-matching (QPM) domain poling specification (`[nonlinear.qpm]`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct QpmSpec {
    /// Base domain orientation sign (`+1.0` or `-1.0`).
    #[serde(default = "default_one")]
    pub sign: f64,
    /// Optional spatial poling period $\Lambda_{\text{QPM}}$ in metres (`None` = uniform `sign` across the domain).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub period: Option<f64>,
    /// Unit poling wavevector direction in global wafer coordinates (defaults to `[0, 1, 0]`, propagation `+y`).
    #[serde(default = "default_y_axis")]
    pub axis: [f64; 3],
    /// Positive-domain duty cycle $D \in (0, 1)$ within each period $\Lambda_{\text{QPM}}$ (default `0.5`).
    #[serde(default = "default_half")]
    pub duty_cycle: f64,
    /// Spatial coordinate offset $s_0$ in metres along `axis` where the `+sign` half-period starts.
    #[serde(default)]
    pub phase_offset: f64,
}

impl Default for QpmSpec {
    fn default() -> Self {
        Self {
            sign: 1.0,
            period: None,
            axis: [0.0, 1.0, 0.0],
            duty_cycle: 0.5,
            phase_offset: 0.0,
        }
    }
}

impl QpmSpec {
    /// Evaluates the local ferroelectric domain poling sign (`+1.0` or `-1.0`) at physical position `pos` (m).
    #[inline]
    pub fn poling_sign_at(&self, pos: [f64; 3]) -> f64 {
        let base = if self.sign < 0.0 { -1.0 } else { 1.0 };
        let Some(lam) = self.period.filter(|&p| p > 0.0 && p.is_finite()) else {
            return base;
        };
        let s = self.axis[0] * pos[0] + self.axis[1] * pos[1] + self.axis[2] * pos[2] - self.phase_offset;
        let u = s.rem_euclid(lam);
        if u < self.duty_cycle * lam {
            base
        } else {
            -base
        }
    }

    /// Validates that all fields of `QpmSpec` are finite and physically well-formed.
    pub fn validate(&self) -> Result<(), String> {
        if !self.sign.is_finite() || (self.sign.abs() - 1.0).abs() > 1e-12 {
            return Err(format!("qpm.sign must be +1 or -1 (got {})", self.sign));
        }
        if let Some(p) = self.period {
            if !p.is_finite() || p <= 0.0 {
                return Err(format!("qpm.Lambda_qpm (period) must be positive and finite (got {p})"));
            }
        }
        if self.axis.iter().any(|v| !v.is_finite()) {
            return Err("qpm.axis components must be finite".into());
        }
        let nrm = (self.axis[0] * self.axis[0] + self.axis[1] * self.axis[1] + self.axis[2] * self.axis[2]).sqrt();
        if nrm <= 1e-30 || !nrm.is_finite() {
            return Err("qpm.axis must be a non-zero vector".into());
        }
        if !self.duty_cycle.is_finite() || self.duty_cycle <= 0.0 || self.duty_cycle >= 1.0 {
            return Err(format!("qpm.duty_cycle must lie strictly in (0, 1) (got {})", self.duty_cycle));
        }
        if !self.phase_offset.is_finite() {
            return Err("qpm.phase_offset must be finite".into());
        }
        Ok(())
    }

    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let sign = s.number("sign")?.unwrap_or(1.0);
        if !sign.is_finite() || (sign.abs() - 1.0).abs() > 1e-12 {
            return Err(s.err("sign", "poling sign must be +1 or -1"));
        }
        let sign = if sign < 0.0 { -1.0 } else { 1.0 };

        let lam = s.quantity("Lambda_qpm", Dim::Length, None)?;
        let per = s.quantity("period", Dim::Length, None)?;
        if lam.is_some() && per.is_some() {
            return Err(s.err("period", "specify either Lambda_qpm or period, not both"));
        }
        let period = lam.or(per);
        if let Some(p) = period {
            if !p.is_finite() || p <= 0.0 {
                return Err(s.err("Lambda_qpm", "poling period must be positive"));
            }
        }

        let axis = match s.numbers("axis")? {
            Some(v) => {
                if v.len() != 3 {
                    return Err(s.err("axis", "axis must have 3 components [ax, ay, az]"));
                }
                let nrm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if !nrm.is_finite() || nrm <= 1e-30 {
                    return Err(s.err("axis", "axis must be a non-zero vector"));
                }
                [v[0] / nrm, v[1] / nrm, v[2] / nrm]
            }
            None => [0.0, 1.0, 0.0],
        };

        let duty_cycle = s.number("duty_cycle")?.unwrap_or(0.5);
        if !duty_cycle.is_finite() || duty_cycle <= 0.0 || duty_cycle >= 1.0 {
            return Err(s.err("duty_cycle", "duty_cycle must lie strictly between 0 and 1"));
        }
        let phase_offset = s.quantity("phase_offset", Dim::Length, Some("m"))?.unwrap_or(0.0);
        if !phase_offset.is_finite() {
            return Err(s.err("phase_offset", "phase_offset must be finite"));
        }

        s.finish(QPM_KEYS)?;
        Ok(Self {
            sign,
            period,
            axis,
            duty_cycle,
            phase_offset,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        tbl.insert("sign".into(), toml::Value::Float(self.sign));
        if let Some(p) = self.period {
            tbl.insert("Lambda_qpm".into(), format!("{p} m").into());
        }
        tbl.insert(
            "axis".into(),
            toml::Value::Array(self.axis.iter().copied().map(toml::Value::Float).collect()),
        );
        if (self.duty_cycle - 0.5).abs() > 1e-12 {
            tbl.insert("duty_cycle".into(), toml::Value::Float(self.duty_cycle));
        }
        if self.phase_offset.abs() > 0.0 {
            tbl.insert("phase_offset".into(), format!("{} m", self.phase_offset).into());
        }
    }
}

/// Second-order nonlinear optical tensor `d_iJ = 0.5 * chi^(2)_iJ` (`[nonlinear.chi2]`) in `m/V`
/// and textbook Voigt order `(xx, yy, zz, yz, xz, xy)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Chi2 {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength: Option<f64>,
    #[serde(default, skip_serializing_if = "is_zero_flag")]
    pub kleinman: bool,
    #[serde(default, skip_serializing_if = "is_zero_flag")]
    pub broken_centrosymmetry: bool,
    /// Contracted $3 \times 6$ $d_{iJ}$ tensor in $\text{m/V}$ (crystal frame, textbook Voigt).
    pub d: [[f64; 6]; 3],
}

fn expand_chi2_named(pg: Option<PointGroup>, kleinman: bool, raw: &[[Option<f64>; 6]; 3]) -> [[f64; 6]; 3] {
    let mut m = [[0.0; 6]; 3];
    for i in 0..3 {
        for j in 0..6 {
            if let Some(v) = raw[i][j] {
                m[i][j] = v;
            }
        }
    }
    match pg {
        Some(PointGroup::C3v) => {
            // Trigonal 3m (IEEE-1987, mirror m perp X, e.g. LiNbO3, BBO):
            // d21 = -d22, d16 = -d22 (no factor 2 in contracted optical d_iJ!), d24 = d15, d32 = d31.
            if raw[1][0].is_none() && raw[1][1].is_some() {
                m[1][0] = -m[1][1];
            }
            if raw[0][5].is_none() && raw[1][1].is_some() {
                m[0][5] = -m[1][1];
            }
            if kleinman && raw[0][4].is_none() && raw[2][0].is_some() {
                m[0][4] = m[2][0];
            } else if kleinman && raw[2][0].is_none() && raw[0][4].is_some() {
                m[2][0] = m[0][4];
            }
            if raw[1][3].is_none() && m[0][4] != 0.0 {
                m[1][3] = m[0][4];
            } else if raw[0][4].is_none() && raw[1][3].is_some() {
                m[0][4] = m[1][3];
                if kleinman && raw[2][0].is_none() {
                    m[2][0] = m[0][4];
                }
            }
            if raw[2][1].is_none() && m[2][0] != 0.0 {
                m[2][1] = m[2][0];
            }
        }
        Some(PointGroup::C6v | PointGroup::C4v | PointGroup::C6 | PointGroup::C4) => {
            // Wurtzite / polar hexagonal/tetragonal 6mm, 4mm (AlN, GaN, ZnO, BaTiO3):
            // d32 = d31, d24 = d15; chiral 6 / 4 also allow d25 = -d14.
            if kleinman && raw[0][4].is_none() && raw[2][0].is_some() {
                m[0][4] = m[2][0];
            } else if kleinman && raw[2][0].is_none() && raw[0][4].is_some() {
                m[2][0] = m[0][4];
            }
            if raw[2][1].is_none() && m[2][0] != 0.0 {
                m[2][1] = m[2][0];
            }
            if raw[1][3].is_none() && m[0][4] != 0.0 {
                m[1][3] = m[0][4];
            }
            if matches!(pg, Some(PointGroup::C6 | PointGroup::C4)) && raw[1][4].is_none() && raw[0][3].is_some() {
                m[1][4] = -m[0][3];
            }
        }
        Some(PointGroup::Td | PointGroup::T23) => {
            // Zincblende cubic -43m / 23 (GaAs, InP, GaP, AlGaAs): d14 = d25 = d36
            let v_opt = raw[0][3].or(raw[1][4]).or(raw[2][5]);
            if let Some(v) = v_opt {
                if raw[0][3].is_none() {
                    m[0][3] = v;
                }
                if raw[1][4].is_none() {
                    m[1][4] = v;
                }
                if raw[2][5].is_none() {
                    m[2][5] = v;
                }
            }
        }
        Some(PointGroup::D2d) => {
            // Chalcopyrite / KDP tetragonal -42m: d14 = d25, independent d36 (equal under Kleinman)
            if raw[1][4].is_none() && raw[0][3].is_some() {
                m[1][4] = m[0][3];
            } else if raw[0][3].is_none() && raw[1][4].is_some() {
                m[0][3] = m[1][4];
            }
            if kleinman {
                if raw[2][5].is_none() && m[0][3] != 0.0 {
                    m[2][5] = m[0][3];
                } else if raw[0][3].is_none() && raw[1][4].is_none() && raw[2][5].is_some() {
                    m[0][3] = m[2][5];
                    m[1][4] = m[2][5];
                }
            }
        }
        Some(PointGroup::D3) => {
            // Trigonal 32 (alpha-quartz): d12 = -d11, d26 = -d11, d25 = -d14
            if raw[0][1].is_none() && raw[0][0].is_some() {
                m[0][1] = -m[0][0];
            }
            if raw[1][5].is_none() && raw[0][0].is_some() {
                m[1][5] = -m[0][0];
            }
            if raw[1][4].is_none() && raw[0][3].is_some() {
                m[1][4] = -m[0][3];
            }
        }
        Some(PointGroup::D3h) => {
            // Hexagonal -6m2 (monolayer MoS2/hBN bulk 2H/3R, GaSe): d21 = -d22, d16 = -d22 (mirror perp x)
            // or d12 = -d11, d26 = -d11 (mirror perp y)
            if raw[1][1].is_some() {
                if raw[1][0].is_none() {
                    m[1][0] = -m[1][1];
                }
                if raw[0][5].is_none() {
                    m[0][5] = -m[1][1];
                }
            } else if raw[0][0].is_some() {
                if raw[0][1].is_none() {
                    m[0][1] = -m[0][0];
                }
                if raw[1][5].is_none() {
                    m[1][5] = -m[0][0];
                }
            }
        }
        Some(PointGroup::C2v) if kleinman => {
            // Orthorhombic mm2 (KTP, LBO, KNbO3): independent d15, d24, d31, d32, d33; under Kleinman d15=d31, d24=d32
            if raw[0][4].is_none() && raw[2][0].is_some() {
                m[0][4] = m[2][0];
            }
            if raw[1][3].is_none() && raw[2][1].is_some() {
                m[1][3] = m[2][1];
            }
        }
        _ => {}
    }
    m
}

/// Validates a $3 \times 6$ $d_{iJ}$ tensor against crystallographic `PointGroup` selection rules
/// unless `broken_centrosymmetry` is explicitly enabled.
pub fn validate_chi2_point_group(
    d: &[[f64; 6]; 3],
    pg: Option<PointGroup>,
    broken_centrosymmetry: bool,
) -> Result<(), String> {
    if d.iter().flatten().any(|v| !v.is_finite()) {
        return Err("chi2 tensor contains NaN or infinite entries".into());
    }
    let Some(pg) = pg else {
        return Ok(());
    };
    if broken_centrosymmetry {
        return Ok(());
    }
    let max_abs = d.iter().flatten().fold(0.0_f64, |acc, &x| acc.max(x.abs()));
    if pg.forbids_piezo_and_pockels() {
        if max_abs > 0.0 {
            return Err(format!(
                "point_group \"{}\" is centrosymmetric (or chiral cubic 432) and forbids bulk second-order susceptibility chi2 (d_iJ = 0) by Neumann's principle; set broken_centrosymmetry = true only for field/strain-broken symmetry",
                pg.as_str()
            ));
        }
        return Ok(());
    }
    if max_abs == 0.0 {
        return Ok(());
    }
    let tol = 1e-9 * max_abs;
    let check_zero = |slots: &[(usize, usize)], cls: &str| -> Result<(), String> {
        for &(i, j) in slots {
            if d[i][j].abs() > tol {
                return Err(format!(
                    "point_group \"{cls}\" forbids non-zero d{}{ } (found {:e} m/V)",
                    i + 1,
                    j + 1,
                    d[i][j]
                ));
            }
        }
        Ok(())
    };
    match pg {
        PointGroup::C3v => {
            // Allowed non-zeros in 3m (mirror perp X): d15(0,4), d16(0,5), d21(1,0), d22(1,1), d24(1,3), d31(2,0), d32(2,1), d33(2,2)
            check_zero(
                &[
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 3),
                    (1, 2),
                    (1, 4),
                    (1, 5),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                ],
                pg.as_str(),
            )?;
            if (d[1][0] + d[1][1]).abs() > tol || (d[0][5] + d[1][1]).abs() > tol {
                return Err(format!(
                    "point_group \"3m\" requires d21 = d16 = -d22 (got d21={:e}, d16={:e}, d22={:e})",
                    d[1][0], d[0][5], d[1][1]
                ));
            }
            if (d[1][3] - d[0][4]).abs() > tol {
                return Err(format!(
                    "point_group \"3m\" requires d24 = d15 (got d24={:e}, d15={:e})",
                    d[1][3], d[0][4]
                ));
            }
            if (d[2][1] - d[2][0]).abs() > tol {
                return Err(format!(
                    "point_group \"3m\" requires d32 = d31 (got d32={:e}, d31={:e})",
                    d[2][1], d[2][0]
                ));
            }
        }
        PointGroup::C6v | PointGroup::C4v => {
            check_zero(
                &[
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 3),
                    (0, 5),
                    (1, 0),
                    (1, 1),
                    (1, 2),
                    (1, 4),
                    (1, 5),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                ],
                pg.as_str(),
            )?;
            if (d[1][3] - d[0][4]).abs() > tol || (d[2][1] - d[2][0]).abs() > tol {
                return Err(format!("point_group \"{}\" requires d24 = d15 and d32 = d31", pg.as_str()));
            }
        }
        PointGroup::Td | PointGroup::T23 => {
            check_zero(
                &[
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 4),
                    (0, 5),
                    (1, 0),
                    (1, 1),
                    (1, 2),
                    (1, 3),
                    (1, 5),
                    (2, 0),
                    (2, 1),
                    (2, 2),
                    (2, 3),
                    (2, 4),
                ],
                pg.as_str(),
            )?;
            if (d[0][3] - d[1][4]).abs() > tol || (d[0][3] - d[2][5]).abs() > tol {
                return Err(format!("point_group \"{}\" requires d14 = d25 = d36", pg.as_str()));
            }
        }
        PointGroup::D2d => {
            check_zero(
                &[
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 4),
                    (0, 5),
                    (1, 0),
                    (1, 1),
                    (1, 2),
                    (1, 3),
                    (1, 5),
                    (2, 0),
                    (2, 1),
                    (2, 2),
                    (2, 3),
                    (2, 4),
                ],
                pg.as_str(),
            )?;
            if (d[0][3] - d[1][4]).abs() > tol {
                return Err("point_group \"-42m\" requires d14 = d25".into());
            }
        }
        PointGroup::C2v => {
            check_zero(
                &[
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 3),
                    (0, 5),
                    (1, 0),
                    (1, 1),
                    (1, 2),
                    (1, 4),
                    (1, 5),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                ],
                pg.as_str(),
            )?;
        }
        _ => {}
    }
    Ok(())
}

fn read_3x6_chi2(s: &Section, key: &str, tree: &Tree) -> Result<[[f64; 6]; 3], MaterialError> {
    let Tree::List(rows) = tree else {
        return Err(s.err(key, "must be a 3x6 array (3 rows of 6 entries)"));
    };
    if rows.len() != 3 {
        return Err(s.err(key, format!("expected 3 rows in 3x6 chi2 matrix, found {}", rows.len())));
    }
    let mut out = [[0.0; 6]; 3];
    for (i, r) in rows.iter().enumerate() {
        let Tree::List(cols) = r else {
            return Err(s.err(key, format!("row {} must be an array of 6 numbers", i + 1)));
        };
        if cols.len() != 6 {
            return Err(s.err(key, format!("row {} has {} entries; expected 6", i + 1, cols.len())));
        }
        for (j, c) in cols.iter().enumerate() {
            let Tree::Num(v) = c else {
                return Err(s.err(key, "matrix entries must be scalar quantities"));
            };
            out[i][j] = *v;
        }
    }
    Ok(out)
}

impl Chi2 {
    /// Expands the contracted $3 \times 6$ $d_{iJ}$ matrix (`(xx, yy, zz, yz, xz, xy)`) into a rank-3 tensor
    /// $d_{ijk}$ symmetric in its last two indices ($d_{ijk} = d_{ikj} = d_{iJ}$ directly, with NO $1/2$ shear factor).
    pub fn to_rank3(&self) -> [[[f64; 3]; 3]; 3] {
        d_voigt_to_rank3(&self.d)
    }

    /// Reconstructs a $3 \times 6$ contracted $d_{iJ}$ matrix from a rank-3 tensor $d_{ijk}$ ($d_{ijk} == d_{ikj}$).
    pub fn from_rank3(t: &[[[f64; 3]; 3]; 3]) -> [[f64; 6]; 3] {
        d_rank3_to_voigt(t)
    }

    /// Rotates $d_{iJ}$ from the crystal frame into a global frame using rotation matrix `r` ($x_{\text{glob}} = r x_{\text{cryst}}$).
    pub fn rotate(&self, r: &M3) -> [[f64; 6]; 3] {
        let t = self.to_rank3();
        let mut out = [[[0.0; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let mut sum = 0.0;
                    for a in 0..3 {
                        let ria = r[i][a];
                        if ria == 0.0 {
                            continue;
                        }
                        for b in 0..3 {
                            let rjb = r[j][b];
                            if rjb == 0.0 {
                                continue;
                            }
                            let ria_rjb = ria * rjb;
                            for c in 0..3 {
                                sum += ria_rjb * r[k][c] * t[a][b][c];
                            }
                        }
                    }
                    out[i][j][k] = sum;
                }
            }
        }
        Self::from_rank3(&out)
    }

    /// Evaluates the global-frame contracted $3 \times 6$ $d_{iJ}$ tensor rotated by Bunge Z-X-Z `euler_deg` (degrees).
    pub fn d_global(&self, euler_deg: [f64; 3]) -> [[f64; 6]; 3] {
        if euler_deg == [0.0; 3] {
            return self.d;
        }
        self.rotate(&rotation(euler_deg))
    }

    pub fn read(
        s: &mut Section,
        pg: Option<PointGroup>,
        parent_broken: bool,
    ) -> Result<(Self, Option<QpmSpec>), MaterialError> {
        let source = read_sources(s, "source")?;
        let wavelength = s.quantity("wavelength", Dim::Length, None)?;
        if let Some(w) = wavelength {
            if !w.is_finite() || w <= 0.0 {
                return Err(s.err("wavelength", "wavelength must be positive"));
            }
        }
        let unit_str = s.string("unit")?;
        let kleinman = s.bool("kleinman")?.unwrap_or(false);
        let broken_centrosymmetry = s.bool("broken_centrosymmetry")?.unwrap_or(parent_broken);
        let nested_qpm = s.sub("qpm")?.map(|mut qs| QpmSpec::read(&mut qs)).transpose()?;

        let mut known: Vec<String> = CHI2_KEYS.iter().map(|&k| k.to_string()).collect();
        for i in 1..=3 {
            for j in 1..=6 {
                known.push(format!("d{i}{j}"));
            }
        }
        let known_refs: Vec<&str> = known.iter().map(String::as_str).collect();

        let d = if s.has("matrix") {
            let tree = s
                .quantities("matrix", Dim::ElectroOptic, unit_str.as_deref())?
                .ok_or_else(|| s.err("matrix", "must be a 3x6 array"))?;
            read_3x6_chi2(s, "matrix", &tree)?
        } else {
            let mut raw = [[None; 6]; 3];
            let mut any = false;
            for i in 0..3 {
                for j in 0..6 {
                    let key = format!("d{}{}", i + 1, j + 1);
                    if let Some(v) = s.quantity(&key, Dim::ElectroOptic, unit_str.as_deref())? {
                        raw[i][j] = Some(v);
                        any = true;
                    }
                }
            }
            if !any {
                return Err(s.err("", "write matrix = [...] or named chi^(2) constants like d33, d31, d15, d22, d14"));
            }
            expand_chi2_named(pg, kleinman, &raw)
        };

        validate_chi2_point_group(&d, pg, broken_centrosymmetry).map_err(|msg| s.err("", msg))?;
        s.finish(&known_refs)?;
        Ok((
            Self {
                source,
                wavelength,
                kleinman,
                broken_centrosymmetry,
                d,
            },
            nested_qpm,
        ))
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(w) = self.wavelength {
            tbl.insert("wavelength".into(), format!("{w} m").into());
        }
        if self.kleinman {
            tbl.insert("kleinman".into(), toml::Value::Boolean(true));
        }
        if self.broken_centrosymmetry {
            tbl.insert("broken_centrosymmetry".into(), toml::Value::Boolean(true));
        }
        tbl.insert("unit".into(), "m/V".into());
        let rows = self
            .d
            .iter()
            .map(|r| toml::Value::Array(r.iter().copied().map(toml::Value::Float).collect()))
            .collect();
        tbl.insert("matrix".into(), toml::Value::Array(rows));
    }
}

/// Expands a $3 \times 6$ contracted optical $d_{iJ}$ matrix (`(xx, yy, zz, yz, xz, xy)`) to rank-3 $d_{ijk}$.
#[inline]
pub fn d_voigt_to_rank3(d: &[[f64; 6]; 3]) -> [[[f64; 3]; 3]; 3] {
    let mut t = [[[0.0; 3]; 3]; 3];
    for i in 0..3 {
        for (vj, &(a, b)) in VP.iter().enumerate() {
            let v = d[i][vj];
            t[i][a][b] = v;
            t[i][b][a] = v;
        }
    }
    t
}

/// Contracts a symmetric-in-`(j,k)` rank-3 tensor $d_{ijk}$ to a $3 \times 6$ optical $d_{iJ}$ matrix.
#[inline]
pub fn d_rank3_to_voigt(t: &[[[f64; 3]; 3]; 3]) -> [[f64; 6]; 3] {
    let mut m = [[0.0; 6]; 3];
    for i in 0..3 {
        for (vj, &(a, b)) in VP.iter().enumerate() {
            m[i][vj] = 0.5 * (t[i][a][b] + t[i][b][a]);
        }
    }
    m
}

/// Rotates a crystal-frame $3 \times 6$ contracted $d_{iJ}$ tensor into the global wafer frame by Bunge Z-X-Z `euler_deg`.
pub fn chi2_d_global(d_crystal: &[[f64; 6]; 3], euler_deg: [f64; 3]) -> [[f64; 6]; 3] {
    let c2 = Chi2 {
        source: Vec::new(),
        wavelength: None,
        kleinman: false,
        broken_centrosymmetry: true,
        d: *d_crystal,
    };
    c2.d_global(euler_deg)
}

/// Third-order Kerr, two-photon absorption (TPA), and free-carrier absorption/dispersion (`[nonlinear.chi3]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Chi3 {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength: Option<f64>,
    /// Isotropic nonlinear Kerr index $n_2$ in $\text{m}^2/\text{W}$ ($\Delta n = n_2 I$).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub n_2: Option<f64>,
    /// Diagonal third-order susceptibility $\chi^{(3)}_{xxxx}$ in $\text{m}^2/\text{V}^2$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi3_xxxx: Option<f64>,
    /// Cross-polarized susceptibility $\chi^{(3)}_{xxyy}$ in $\text{m}^2/\text{V}^2$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi3_xxyy: Option<f64>,
    /// Cross-polarized susceptibility $\chi^{(3)}_{xyxy}$ in $\text{m}^2/\text{V}^2$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi3_xyxy: Option<f64>,
    /// Conjugate cross-polarized susceptibility $\chi^{(3)}_{xyyx}$ in $\text{m}^2/\text{V}^2$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi3_xyyx: Option<f64>,
    /// Maker–Terhune cross-polarization ratio $B / A$ when parameterized via `n_2` or `chi3_xxxx` alone (`1.0` = electronic Kleinman).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub b_over_a: Option<f64>,
    /// Two-photon absorption coefficient $\beta_{\text{TPA}} \ge 0$ in $\text{m/W}$ ($\alpha_{\text{TPA}} = \beta_{\text{TPA}} I$).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beta_tpa: Option<f64>,
    /// Free-carrier absorption cross-section $\sigma_{\text{FCA}} \ge 0$ in $\text{m}^2$ ($\alpha_{\text{FCA}} = \sigma_{\text{FCA}} N_c$).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sigma_fca: Option<f64>,
    /// Signed free-carrier dispersion coefficient $k_{\text{FCD}}$ in $\text{m}^3$ ($\Delta n_{\text{FCD}} = k_{\text{FCD}} N_c$).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_fcd: Option<f64>,
    /// Free-carrier recombination lifetime $\tau_c \ge 0$ in $\text{s}$.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tau_c: Option<f64>,
}

impl Chi3 {
    /// Returns `true` when at least one Kerr, TPA, or free-carrier coefficient is non-zero.
    pub fn is_active(&self) -> bool {
        let nz = |v: Option<f64>| v.is_some_and(|x| x != 0.0);
        nz(self.n_2)
            || nz(self.chi3_xxxx)
            || nz(self.chi3_xxyy)
            || nz(self.chi3_xyxy)
            || nz(self.chi3_xyyx)
            || nz(self.beta_tpa)
            || (nz(self.tau_c) && (nz(self.sigma_fca) || nz(self.k_fcd)))
    }

    /// True if cubic (`m-3m` / `-43m`) crystalline Kerr anisotropy $\sigma_{\text{cub}} = \chi_{xxxx} - (\chi_{xxyy} + \chi_{xyxy} + \chi_{xyyx})$ is non-zero.
    pub fn has_cubic_anisotropy(&self, n0: f64) -> bool {
        self.cubic_anisotropy_term(n0).abs() > 0.0
    }

    /// Evaluates Maker–Terhune coefficients `(A, B)` in $\text{m}^2/\text{V}^2$ for
    /// $\mathbf{P}^{(3)}(\omega) = \varepsilon_0 \big(A (\mathbf{E}\cdot\mathbf{E}^*)\mathbf{E} + \frac{1}{2} B (\mathbf{E}\cdot\mathbf{E})\mathbf{E}^*\big)$
    /// where linear polarization sees $A + \frac{1}{2}B = \frac{3}{4}\chi^{(3)}_{xxxx} = c_0 \varepsilon_0 n_0^2 n_2$.
    pub fn isotropic_ab(&self, n0: f64) -> (f64, f64) {
        if self.chi3_xxyy.is_some() || self.chi3_xyxy.is_some() || self.chi3_xyyx.is_some() {
            let xxyy = self.chi3_xxyy.or(self.chi3_xyxy).unwrap_or(0.0);
            let xyxy = self.chi3_xyxy.unwrap_or(xxyy);
            let default_xyyx = if let Some(xxxx) = self.chi3_xxxx {
                xxxx - xxyy - xyxy
            } else {
                xxyy
            };
            let xyyx = self.chi3_xyyx.unwrap_or(default_xyyx);
            let a = 0.75 * (xxyy + xyxy);
            let b = 1.5 * xyyx;
            return (a, b);
        }
        let chi_eff = if let Some(xxxx) = self.chi3_xxxx {
            0.75 * xxxx
        } else if let Some(n2) = self.n_2 {
            C0_NL * EPSILON_0 * n0 * n0 * n2
        } else {
            return (0.0, 0.0);
        };
        let r = self.b_over_a.unwrap_or(1.0);
        let denom = 1.0 + 0.5 * r;
        if denom.abs() > 1e-12 {
            let a = chi_eff / denom;
            (a, r * a)
        } else {
            (chi_eff, 0.0)
        }
    }

    /// Cubic excess susceptibility $\delta\chi^{(3)} = \chi^{(3)}_{xxxx} - (\chi^{(3)}_{xxyy} + \chi^{(3)}_{xyxy} + \chi^{(3)}_{xyyx})$ in $\text{m}^2/\text{V}^2$.
    pub fn cubic_anisotropy_term(&self, _n0: f64) -> f64 {
        let Some(xxxx) = self.chi3_xxxx else {
            return 0.0;
        };
        if self.chi3_xxyy.is_none() && self.chi3_xyxy.is_none() && self.chi3_xyyx.is_none() {
            return 0.0;
        }
        let xxyy = self.chi3_xxyy.or(self.chi3_xyxy).unwrap_or(0.0);
        let xyxy = self.chi3_xyxy.unwrap_or(xxyy);
        let xyyx = self.chi3_xyyx.unwrap_or(xxyy);
        xxxx - (xxyy + xyxy + xyyx)
    }

    /// Builds the full rank-4 susceptibility tensor $\chi^{(3),\text{glob}}_{ijkl}(-\omega; \omega, \omega, -\omega)$
    /// in global wafer coordinates (symmetric in `(j, k)`), rotated from crystal axes by Bunge Z-X-Z `euler_deg`,
    /// such that $P_i^{(3,\text{Kerr})}(\omega) = \frac{3}{4}\varepsilon_0 \sum_{j,k,l=0}^2 \chi^{\text{glob}}_{ijkl} E_j E_k E_l^*$.
    pub fn rank4_global(&self, n0: f64, euler_deg: [f64; 3]) -> [[[[f64; 3]; 3]; 3]; 3] {
        let (a_mt, b_mt) = self.isotropic_ab(n0);
        let cub = self.cubic_anisotropy_term(n0);
        // Since (3/4) * sum_{j,k,l} [c1*(d_ij d_kl + d_ik d_jl) + c2 * d_il d_jk] E_j E_k E_l^*
        // = (3/2)*c1 * (E.E^*) E_i + (3/4)*c2 * (E.E) E_i^*, we have c1 = (2/3)*A and c2 = (2/3)*B.
        let c1 = (2.0 / 3.0) * a_mt;
        let c2 = (2.0 / 3.0) * b_mt;
        let mut cryst = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut v = 0.0;
                        if i == j && k == l {
                            v += c1;
                        }
                        if i == k && j == l {
                            v += c1;
                        }
                        if i == l && j == k {
                            v += c2;
                        }
                        if i == j && j == k && k == l {
                            v += cub;
                        }
                        cryst[i][j][k][l] = v;
                    }
                }
            }
        }
        if cub.abs() <= 1e-36 || euler_deg == [0.0; 3] {
            return cryst;
        }
        let r = rotation(euler_deg);
        let mut glob = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0;
                        for a in 0..3 {
                            let ria = r[i][a];
                            if ria == 0.0 {
                                continue;
                            }
                            for b in 0..3 {
                                let rjb = r[j][b];
                                if rjb == 0.0 {
                                    continue;
                                }
                                let r_ij = ria * rjb;
                                for c in 0..3 {
                                    let rkc = r[k][c];
                                    if rkc == 0.0 {
                                        continue;
                                    }
                                    let r_ijk = r_ij * rkc;
                                    for d in 0..3 {
                                        s += r_ijk * r[l][d] * cryst[a][b][c][d];
                                    }
                                }
                            }
                        }
                        glob[i][j][k][l] = s;
                    }
                }
            }
        }
        glob
    }

    /// Evaluates steady-state TPA-generated free-carrier density
    /// $N_c(\mathbf{E}) = \frac{\tau_c \beta_{\text{TPA}}}{2 \hbar \omega} I(\mathbf{E})^2$ ($\text{m}^{-3}$)
    /// for optical angular frequency `omega` ($\text{rad/s}$) and local refractive index `n0`.
    #[inline]
    pub fn steady_carrier_density(&self, omega: f64, n0: f64, e_peak: &[Complex64; 3]) -> f64 {
        let beta = self.beta_tpa.unwrap_or(0.0);
        let tau = self.tau_c.unwrap_or(0.0);
        if beta <= 0.0 || tau <= 0.0 || omega <= 0.0 {
            return 0.0;
        }
        let e_norm2 = e_peak[0].norm_sqr() + e_peak[1].norm_sqr() + e_peak[2].norm_sqr();
        let intensity = 0.5 * C0_NL * EPSILON_0 * n0 * e_norm2;
        (tau * beta / (2.0 * HBAR * omega)) * intensity * intensity
    }

    pub fn validate(&self) -> Result<(), String> {
        for (name, val) in [
            ("wavelength", self.wavelength),
            ("n_2", self.n_2),
            ("chi3_xxxx", self.chi3_xxxx),
            ("chi3_xxyy", self.chi3_xxyy),
            ("chi3_xyxy", self.chi3_xyxy),
            ("chi3_xyyx", self.chi3_xyyx),
            ("b_over_a", self.b_over_a),
            ("beta_tpa", self.beta_tpa),
            ("sigma_fca", self.sigma_fca),
            ("k_fcd", self.k_fcd),
            ("tau_c", self.tau_c),
        ] {
            if let Some(v) = val {
                if !v.is_finite() {
                    return Err(format!("chi3.{name} must be a finite number"));
                }
            }
        }
        if let Some(w) = self.wavelength {
            if w <= 0.0 {
                return Err("chi3.wavelength must be positive".into());
            }
        }
        if let Some(b) = self.beta_tpa {
            if b < 0.0 {
                return Err(format!("chi3.beta_tpa must be non-negative (passive TPA), got {b}"));
            }
        }
        if let Some(s) = self.sigma_fca {
            if s < 0.0 {
                return Err(format!("chi3.sigma_fca must be non-negative (passive FCA), got {s}"));
            }
        }
        if let Some(t) = self.tau_c {
            if t < 0.0 {
                return Err(format!("chi3.tau_c must be non-negative, got {t}"));
            }
        }
        Ok(())
    }

    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let wavelength = s.quantity("wavelength", Dim::Length, None)?;
        let n_2 = s.quantity("n_2", Dim::NonlinearIndex, None)?;
        let chi_unit = s.string("unit")?;
        let chi3_xxxx = s.quantity("chi3_xxxx", Dim::Chi3Susceptibility, chi_unit.as_deref())?;
        let chi3_xxyy = s.quantity("chi3_xxyy", Dim::Chi3Susceptibility, chi_unit.as_deref())?;
        let chi3_xyxy = s.quantity("chi3_xyxy", Dim::Chi3Susceptibility, chi_unit.as_deref())?;
        let chi3_xyyx = s.quantity("chi3_xyyx", Dim::Chi3Susceptibility, chi_unit.as_deref())?;
        let b_over_a = s.number("b_over_a")?;
        let beta_tpa = s.quantity("beta_tpa", Dim::TwoPhotonAbsorption, None)?;
        let sigma_fca = s.quantity("sigma_fca", Dim::Area, None)?;
        let k_fcd = s.quantity("k_fcd", Dim::Volume, None)?;
        let tau_c = s.quantity("tau_c", Dim::Time, None)?;

        let out = Self {
            source,
            wavelength,
            n_2,
            chi3_xxxx,
            chi3_xxyy,
            chi3_xyxy,
            chi3_xyyx,
            b_over_a,
            beta_tpa,
            sigma_fca,
            k_fcd,
            tau_c,
        };
        if out.n_2.is_none()
            && out.chi3_xxxx.is_none()
            && out.chi3_xxyy.is_none()
            && out.beta_tpa.is_none()
            && out.sigma_fca.is_none()
            && out.k_fcd.is_none()
        {
            return Err(s.err(
                "",
                "specify at least one chi^(3) parameter (n_2, chi3_xxxx, chi3_xxyy, beta_tpa, sigma_fca, or k_fcd)",
            ));
        }
        out.validate().map_err(|msg| s.err("", msg))?;
        s.finish(CHI3_KEYS)?;
        Ok(out)
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(w) = self.wavelength {
            tbl.insert("wavelength".into(), format!("{w} m").into());
        }
        if let Some(n2) = self.n_2 {
            tbl.insert("n_2".into(), format!("{n2} m^2/W").into());
        }
        if self.chi3_xxxx.is_some()
            || self.chi3_xxyy.is_some()
            || self.chi3_xyxy.is_some()
            || self.chi3_xyyx.is_some()
        {
            tbl.insert("unit".into(), "m^2/V^2".into());
        }
        if let Some(v) = self.chi3_xxxx {
            tbl.insert("chi3_xxxx".into(), toml::Value::Float(v));
        }
        if let Some(v) = self.chi3_xxyy {
            tbl.insert("chi3_xxyy".into(), toml::Value::Float(v));
        }
        if let Some(v) = self.chi3_xyxy {
            tbl.insert("chi3_xyxy".into(), toml::Value::Float(v));
        }
        if let Some(v) = self.chi3_xyyx {
            tbl.insert("chi3_xyyx".into(), toml::Value::Float(v));
        }
        if let Some(v) = self.b_over_a {
            tbl.insert("b_over_a".into(), toml::Value::Float(v));
        }
        if let Some(v) = self.beta_tpa {
            tbl.insert("beta_tpa".into(), format!("{v} m/W").into());
        }
        if let Some(v) = self.sigma_fca {
            tbl.insert("sigma_fca".into(), format!("{v} m^2").into());
        }
        if let Some(v) = self.k_fcd {
            tbl.insert("k_fcd".into(), format!("{v} m^3").into());
        }
        if let Some(v) = self.tau_c {
            tbl.insert("tau_c".into(), format!("{v} s").into());
        }
    }
}

/// Complete optical nonlinearity material block (`[nonlinear]`, `format = 2`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Nonlinear {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength: Option<f64>,
    /// Valid temperature interval `[T_min, T_max]` in Kelvin.
    #[serde(rename = "T_range", default, skip_serializing_if = "Option::is_none")]
    pub t_range: Option<[f64; 2]>,
    /// Allows non-zero `chi2` on a nominally centrosymmetric `point_group` when inversion symmetry is broken by DC field, strain, or poling.
    #[serde(default, skip_serializing_if = "is_zero_flag")]
    pub broken_centrosymmetry: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi2: Option<Chi2>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chi3: Option<Chi3>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub qpm: Option<QpmSpec>,
}

impl Nonlinear {
    /// Checks `temperature_k` against `self.t_range` (`[T_min, T_max]`), returning clamping metadata if out of bounds.
    pub fn check_temperature(&self, temperature_k: f64) -> CheckedTempEval<f64> {
        match self.t_range {
            Some([lo, hi]) => {
                let (clamped, info) = clamp_temp_range(temperature_k, lo, hi);
                CheckedTempEval { value: clamped, clamp: info }
            }
            None => CheckedTempEval { value: temperature_k, clamp: None },
        }
    }

    pub fn read(s: &mut Section, pg: Option<PointGroup>) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let wavelength = s.quantity("wavelength", Dim::Length, None)?;
        if let Some(w) = wavelength {
            if !w.is_finite() || w <= 0.0 {
                return Err(s.err("wavelength", "wavelength must be positive"));
            }
        }
        let t_unit = s.string("T_unit")?.unwrap_or_else(|| "K".to_string());
        let t_range = match s.numbers("T_range")? {
            Some(v) => {
                if v.len() != 2 {
                    return Err(s.err("T_range", "T_range must have 2 elements [T_min, T_max]"));
                }
                let t0 = super::units::to_si(v[0], &t_unit, Dim::Temperature).map_err(|m| s.err("T_range", m))?;
                let t1 = super::units::to_si(v[1], &t_unit, Dim::Temperature).map_err(|m| s.err("T_range", m))?;
                if !t0.is_finite() || !t1.is_finite() || t0 < 0.0 || t1 <= t0 {
                    return Err(s.err("T_range", "T_range must satisfy 0 <= T_min < T_max"));
                }
                Some([t0, t1])
            }
            None => None,
        };
        let broken_centrosymmetry = s.bool("broken_centrosymmetry")?.unwrap_or(false);

        let (chi2, nested_qpm) = match s.sub("chi2")? {
            Some(mut sub) => {
                let (c2, q) = Chi2::read(&mut sub, pg, broken_centrosymmetry)?;
                (Some(c2), q)
            }
            None => (None, None),
        };
        let chi3 = s.sub("chi3")?.map(|mut sub| Chi3::read(&mut sub)).transpose()?;
        let top_qpm = s.sub("qpm")?.map(|mut sub| QpmSpec::read(&mut sub)).transpose()?;
        if top_qpm.is_some() && nested_qpm.is_some() {
            return Err(s.err("qpm", "specify [nonlinear.qpm] in one place (not both [nonlinear.qpm] and [nonlinear.chi2.qpm])"));
        }
        let qpm = top_qpm.or(nested_qpm);

        if chi2.is_none() && chi3.is_none() {
            return Err(s.err("", "[nonlinear] must define [nonlinear.chi2], [nonlinear.chi3], or both"));
        }
        s.finish(NONLINEAR_KEYS)?;
        Ok(Self {
            source,
            wavelength,
            t_range,
            broken_centrosymmetry,
            chi2,
            chi3,
            qpm,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        if let Some(w) = self.wavelength {
            tbl.insert("wavelength".into(), format!("{w} m").into());
        }
        if let Some([lo, hi]) = self.t_range {
            tbl.insert(
                "T_range".into(),
                toml::Value::Array(vec![toml::Value::Float(lo), toml::Value::Float(hi)]),
            );
        }
        if self.broken_centrosymmetry {
            tbl.insert("broken_centrosymmetry".into(), toml::Value::Boolean(true));
        }
        if let Some(c2) = &self.chi2 {
            let mut sub = toml::Table::new();
            c2.write(&mut sub);
            tbl.insert("chi2".into(), toml::Value::Table(sub));
        }
        if let Some(c3) = &self.chi3 {
            let mut sub = toml::Table::new();
            c3.write(&mut sub);
            tbl.insert("chi3".into(), toml::Value::Table(sub));
        }
        if let Some(q) = &self.qpm {
            let mut sub = toml::Table::new();
            q.write(&mut sub);
            tbl.insert("qpm".into(), toml::Value::Table(sub));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Shared Frequency-Domain & Instantaneous Time-Domain Nonlinear Polarization Evaluators
// ---------------------------------------------------------------------------------------------

/// Evaluates the instantaneous real time-domain second-order polarization $\mathbf{P}^{(2)}(t)$ ($\text{C/m}^2$)
/// for instantaneous real electric field `e_t = [Ex(t), Ey(t), Ez(t)]` ($\text{V/m}$), global-frame contracted
/// $3 \times 6$ $d_{iJ}$ tensor ($\text{m/V}$, textbook Voigt), and local QPM poling sign `qpm_sign` (`+1` or `-1`).
#[inline]
pub fn eval_chi2_polarization_time(d_global: &[[f64; 6]; 3], qpm_sign: f64, e_t: [f64; 3]) -> [f64; 3] {
    let [ex, ey, ez] = e_t;
    let voigt_e2 = [ex * ex, ey * ey, ez * ez, 2.0 * ey * ez, 2.0 * ex * ez, 2.0 * ex * ey];
    let pref = 2.0 * EPSILON_0 * qpm_sign;
    let mut p = [0.0; 3];
    for i in 0..3 {
        let row = &d_global[i];
        p[i] = pref
            * (row[0] * voigt_e2[0]
                + row[1] * voigt_e2[1]
                + row[2] * voigt_e2[2]
                + row[3] * voigt_e2[3]
                + row[4] * voigt_e2[4]
                + row[5] * voigt_e2[5]);
    }
    p
}

/// Evaluates the instantaneous tangent susceptibility $\frac{1}{\varepsilon_0}\frac{\partial \mathbf{P}^{(2)}(t)}{\partial \mathbf{E}(t)}$
/// ($3 \times 3$ dimensionless matrix) for implicit local constitutive updates in time-domain solvers (`dgtd`, `fdtd`).
#[inline]
pub fn eval_chi2_tangent_time(d_global: &[[f64; 6]; 3], qpm_sign: f64, e_t: [f64; 3]) -> M3 {
    let [ex, ey, ez] = e_t;
    let s = 4.0 * qpm_sign;
    let mut jac = [[0.0; 3]; 3];
    for i in 0..3 {
        let r = &d_global[i];
        jac[i][0] = s * (r[0] * ex + r[4] * ez + r[5] * ey);
        jac[i][1] = s * (r[1] * ey + r[3] * ez + r[5] * ex);
        jac[i][2] = s * (r[2] * ez + r[3] * ey + r[4] * ex);
    }
    jac
}

/// Evaluates the instantaneous real time-domain third-order Kerr polarization $\mathbf{P}^{(3)}(t)$ ($\text{C/m}^2$)
/// for instantaneous real electric field `e_t = [Ex(t), Ey(t), Ez(t)]` ($\text{V/m}$).
///
/// Uses an $O(1)$ isotropic fast path when cubic anisotropy is absent, and exact rotated rank-4 contraction
/// when cubic (`m-3m`) anisotropy is present.
#[inline]
pub fn eval_chi3_kerr_time(chi3: &Chi3, n0: f64, euler_deg: [f64; 3], e_t: [f64; 3]) -> [f64; 3] {
    if !chi3.has_cubic_anisotropy(n0) {
        let (a_mt, b_mt) = chi3.isotropic_ab(n0);
        let chi_xxxx = (4.0 / 3.0) * (a_mt + 0.5 * b_mt);
        let e2 = e_t[0] * e_t[0] + e_t[1] * e_t[1] + e_t[2] * e_t[2];
        let scale = EPSILON_0 * chi_xxxx * e2;
        return [scale * e_t[0], scale * e_t[1], scale * e_t[2]];
    }
    let rank4 = chi3.rank4_global(n0, euler_deg);
    let mut p = [0.0; 3];
    for i in 0..3 {
        let mut s = 0.0;
        for j in 0..3 {
            for k in 0..3 {
                let ejk = e_t[j] * e_t[k];
                for l in 0..3 {
                    s += rank4[i][j][k][l] * ejk * e_t[l];
                }
            }
        }
        p[i] = EPSILON_0 * s;
    }
    p
}

/// Evaluates the complex envelope Kerr polarization $\mathbf{P}^{(3)}(\omega)$ ($\text{C/m}^2$) for complex peak-phasor
/// electric field `e_env = [Ex, Ey, Ez]` ($\text{V/m}$):
/// $$\mathbf{P}^{(3)}(\omega) = \varepsilon_0 \left[ A (\mathbf{E}\cdot\mathbf{E}^*)\mathbf{E} + \frac{1}{2} B (\mathbf{E}\cdot\mathbf{E})\mathbf{E}^* \right]$$
/// (plus rotated cubic `m-3m` anisotropy when present).
#[inline]
pub fn eval_chi3_envelope_time(
    chi3: &Chi3,
    n0: f64,
    euler_deg: [f64; 3],
    e_env: [Complex64; 3],
) -> [Complex64; 3] {
    if !chi3.has_cubic_anisotropy(n0) {
        let (a_mt, b_mt) = chi3.isotropic_ab(n0);
        let e_dot_ec = e_env[0].norm_sqr() + e_env[1].norm_sqr() + e_env[2].norm_sqr();
        let e_dot_e = e_env[0] * e_env[0] + e_env[1] * e_env[1] + e_env[2] * e_env[2];
        let c_a = EPSILON_0 * a_mt * e_dot_ec;
        let c_b = EPSILON_0 * 0.5 * b_mt * e_dot_e;
        return [
            e_env[0] * c_a + c_b * e_env[0].conj(),
            e_env[1] * c_a + c_b * e_env[1].conj(),
            e_env[2] * c_a + c_b * e_env[2].conj(),
        ];
    }
    let rank4 = chi3.rank4_global(n0, euler_deg);
    let pref = 0.75 * EPSILON_0;
    let mut out = [Complex64::new(0.0, 0.0); 3];
    for i in 0..3 {
        let mut acc = Complex64::new(0.0, 0.0);
        for j in 0..3 {
            for k in 0..3 {
                let ejk = e_env[j] * e_env[k];
                for l in 0..3 {
                    acc += ejk * e_env[l].conj() * rank4[i][j][k][l];
                }
            }
        }
        out[i] = acc * pref;
    }
    out
}

/// Unconditionally stable exact exponential-integrator step for the TPA-generated free-carrier density ODE
/// $$\frac{\mathrm{d}N_c}{\mathrm{d}t} = \frac{\beta_{\text{TPA}}\,I(t)^2}{2\hbar\omega_0} - \frac{N_c}{\tau_c},$$
/// advancing `n_c` ($\text{m}^{-3}$) across `dt` ($\text{s}$) at carrier angular frequency `omega_0` ($\text{rad/s}$)
/// and local optical intensity `intensity_w_m2` ($\text{W/m}^2$).
///
/// Converges smoothly to `n_c + G * dt` as `dt / tau_c -> 0` and to the exact steady-state density
/// $N_c = \frac{\tau_c \beta_{\text{TPA}}}{2\hbar\omega_0} I^2$ as `dt / tau_c -> infinity`.
#[inline]
pub fn drude_carrier_step(
    n_c: f64,
    dt: f64,
    omega_0: f64,
    beta_tpa: f64,
    tau_c: f64,
    intensity_w_m2: f64,
) -> f64 {
    let n_prev = n_c.max(0.0);
    if dt <= 0.0 || !dt.is_finite() {
        return n_prev;
    }
    let gen = if beta_tpa > 0.0 && omega_0 > 0.0 && intensity_w_m2 > 0.0 {
        (beta_tpa * intensity_w_m2 * intensity_w_m2) / (2.0 * HBAR * omega_0)
    } else {
        0.0
    };
    if tau_c > 0.0 && tau_c.is_finite() {
        let x = dt / tau_c;
        if x > 700.0 {
            gen * tau_c
        } else {
            let decay = (-x).exp();
            let one_minus_decay = -(-x).exp_m1();
            (n_prev * decay + gen * tau_c * one_minus_decay).max(0.0)
        }
    } else {
        (n_prev + gen * dt).max(0.0)
    }
}
