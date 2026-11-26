// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Crystallographic metadata (`MaterialKind`, `PointGroup`, `CrystalSetting`) and temperature-dependent
//! scalar/uniaxial evaluation (`TempCurve`, `TempTensor`) for format 2 material definitions.

use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
use super::tensor::TensorProp;
use super::units::{self, Dim};
use serde::{Deserialize, Serialize};

/// Macroscopic structural classification of a material.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaterialKind {
    Crystalline,
    Polycrystalline,
    Amorphous,
    Superconductor,
    Semiconductor,
    Dielectric,
    Metal,
    Polymer,
    Fluid,
}

impl MaterialKind {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "crystalline" | "single_crystal" | "single-crystal" | "crystal" => Some(Self::Crystalline),
            "polycrystalline" | "poly" | "ceramic" => Some(Self::Polycrystalline),
            "amorphous" | "glass" => Some(Self::Amorphous),
            "superconductor" | "superconducting" => Some(Self::Superconductor),
            "semiconductor" => Some(Self::Semiconductor),
            "dielectric" | "insulator" => Some(Self::Dielectric),
            "metal" => Some(Self::Metal),
            "polymer" => Some(Self::Polymer),
            "fluid" | "liquid" | "gas" => Some(Self::Fluid),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Crystalline => "crystalline",
            Self::Polycrystalline => "polycrystalline",
            Self::Amorphous => "amorphous",
            Self::Superconductor => "superconductor",
            Self::Semiconductor => "semiconductor",
            Self::Dielectric => "dielectric",
            Self::Metal => "metal",
            Self::Polymer => "polymer",
            Self::Fluid => "fluid",
        }
    }
}

/// Crystal axis frame convention for trigonal/orthorhombic/piezoelectric tensors.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrystalSetting {
    #[default]
    #[serde(rename = "IEEE-1987")]
    Ieee1987,
    #[serde(rename = "IRE-1949")]
    Ire1949,
}

impl CrystalSetting {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_uppercase().replace('_', "-").as_str() {
            "IEEE-1987" | "IEEE1987" => Some(Self::Ieee1987),
            "IRE-1949" | "IRE1949" => Some(Self::Ire1949),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ieee1987 => "IEEE-1987",
            Self::Ire1949 => "IRE-1949",
        }
    }
}

/// All 32 Hermann–Mauguin crystallographic point groups plus isotropic and transverse-isotropic continuous groups.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PointGroup {
    // Triclinic
    C1,
    Ci,
    // Monoclinic
    C2,
    Cs,
    C2h,
    // Orthorhombic
    D2,
    C2v,
    D2h,
    // Tetragonal
    C4,
    S4,
    C4h,
    D4,
    C4v,
    D2d,
    D4h,
    // Trigonal
    C3,
    C3i,
    D3,
    C3v,
    D3d,
    // Hexagonal
    C6,
    C3h,
    C6h,
    D6,
    C6v,
    D3h,
    D6h,
    // Cubic
    T23,
    Th,
    O432,
    Td,
    Oh,
    // Continuous
    Isotropic,
    TransverseIsotropic,
}

impl PointGroup {
    pub const ALL: &'static [PointGroup] = &[
        Self::C1,
        Self::Ci,
        Self::C2,
        Self::Cs,
        Self::C2h,
        Self::D2,
        Self::C2v,
        Self::D2h,
        Self::C4,
        Self::S4,
        Self::C4h,
        Self::D4,
        Self::C4v,
        Self::D2d,
        Self::D4h,
        Self::C3,
        Self::C3i,
        Self::D3,
        Self::C3v,
        Self::D3d,
        Self::C6,
        Self::C3h,
        Self::C6h,
        Self::D6,
        Self::C6v,
        Self::D3h,
        Self::D6h,
        Self::T23,
        Self::Th,
        Self::O432,
        Self::Td,
        Self::Oh,
        Self::Isotropic,
        Self::TransverseIsotropic,
    ];

    pub fn parse(s: &str) -> Option<Self> {
        let t = s.trim();
        let lower = t.to_ascii_lowercase();
        match lower.as_str() {
            "1" => Some(Self::C1),
            "-1" => Some(Self::Ci),
            "2" => Some(Self::C2),
            "m" => Some(Self::Cs),
            "2/m" => Some(Self::C2h),
            "222" => Some(Self::D2),
            "mm2" | "2mm" | "m2m" => Some(Self::C2v),
            "mmm" | "2/m2/m2/m" => Some(Self::D2h),
            "4" => Some(Self::C4),
            "-4" => Some(Self::S4),
            "4/m" => Some(Self::C4h),
            "422" => Some(Self::D4),
            "4mm" => Some(Self::C4v),
            "-42m" | "-4m2" => Some(Self::D2d),
            "4/mmm" | "4/m2/m2/m" => Some(Self::D4h),
            "3" => Some(Self::C3),
            "-3" => Some(Self::C3i),
            "32" => Some(Self::D3),
            "3m" => Some(Self::C3v),
            "-3m" | "-32/m" => Some(Self::D3d),
            "6" => Some(Self::C6),
            "-6" => Some(Self::C3h),
            "6/m" => Some(Self::C6h),
            "622" => Some(Self::D6),
            "6mm" => Some(Self::C6v),
            "-6m2" | "-62m" => Some(Self::D3h),
            "6/mmm" | "6/m2/m2/m" => Some(Self::D6h),
            "23" => Some(Self::T23),
            "m-3" | "m3" | "2/m-3" => Some(Self::Th),
            "432" => Some(Self::O432),
            "-43m" => Some(Self::Td),
            "m-3m" | "m3m" | "4/m-32/m" => Some(Self::Oh),
            "isotropic" | "iso" => Some(Self::Isotropic),
            "transverse_isotropic" | "transverse-isotropic" | "uniaxial" => Some(Self::TransverseIsotropic),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::C1 => "1",
            Self::Ci => "-1",
            Self::C2 => "2",
            Self::Cs => "m",
            Self::C2h => "2/m",
            Self::D2 => "222",
            Self::C2v => "mm2",
            Self::D2h => "mmm",
            Self::C4 => "4",
            Self::S4 => "-4",
            Self::C4h => "4/m",
            Self::D4 => "422",
            Self::C4v => "4mm",
            Self::D2d => "-42m",
            Self::D4h => "4/mmm",
            Self::C3 => "3",
            Self::C3i => "-3",
            Self::D3 => "32",
            Self::C3v => "3m",
            Self::D3d => "-3m",
            Self::C6 => "6",
            Self::C3h => "-6",
            Self::C6h => "6/m",
            Self::D6 => "622",
            Self::C6v => "6mm",
            Self::D3h => "-6m2",
            Self::D6h => "6/mmm",
            Self::T23 => "23",
            Self::Th => "m-3",
            Self::O432 => "432",
            Self::Td => "-43m",
            Self::Oh => "m-3m",
            Self::Isotropic => "isotropic",
            Self::TransverseIsotropic => "transverse_isotropic",
        }
    }

    /// True for the 11 centrosymmetric Laue point groups and continuous isotropic/transverse-isotropic media.
    pub fn is_centrosymmetric(self) -> bool {
        matches!(
            self,
            Self::Ci
                | Self::C2h
                | Self::D2h
                | Self::C4h
                | Self::D4h
                | Self::C3i
                | Self::D3d
                | Self::C6h
                | Self::D6h
                | Self::Th
                | Self::Oh
                | Self::Isotropic
                | Self::TransverseIsotropic
        )
    }

    /// True if Neumann's principle forces all linear piezoelectric and Pockels (rank-3 polar) components to zero.
    /// Note that cubic chiral group `432` lacks an inversion center yet its rank-3 symmetric polar tensor also vanishes identically.
    pub fn forbids_piezo_and_pockels(self) -> bool {
        self.is_centrosymmetric() || matches!(self, Self::O432)
    }
}

impl Serialize for PointGroup {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for PointGroup {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        Self::parse(&s).ok_or_else(|| serde::de::Error::custom(format!("unknown crystallographic point group '{s}'")))
    }
}

fn default_t_ref() -> f64 {
    293.15
}

/// Einstein oscillator thermal shift helper `g(T) = theta_e / (exp(theta_e / T) - 1)`,
/// smoothly saturating to `0` with `dg/dT -> 0` as `T -> 0+` (Third Law of Thermodynamics).
pub fn einstein_thermal_coordinate(temperature_k: f64, theta_e_k: f64) -> f64 {
    if !temperature_k.is_finite() || temperature_k <= 0.0 || !theta_e_k.is_finite() || theta_e_k <= 0.0 {
        return temperature_k.max(0.0);
    }
    let ratio = theta_e_k / temperature_k;
    if ratio > 700.0 {
        return 0.0;
    }
    if ratio < 1e-6 {
        // Taylor expansion: theta / (theta/T + 0.5*(theta/T)^2) = T - 0.5*theta_e + ...
        return temperature_k / (1.0 + 0.5 * ratio + ratio * ratio / 6.0);
    }
    theta_e_k / ratio.exp_m1()
}

/// Allowed TOML keys inside a temperature-curve table (`$defs.temp_curve` in `material_schema`).
pub const TEMP_CURVE_KEYS: &[&str] = &[
    "form",
    "source",
    "T",
    "T_unit",
    "T_range",
    "values",
    "values_unit",
    "coeffs",
    "ref_value",
    "t_ref",
    "theta_e",
    "alpha",
    "unit",
];

/// Metadata recorded when a temperature-dependent constitutive curve or axis model is clamped to its valid domain `[T_min, T_max]`.
#[derive(Clone, Debug, PartialEq)]
pub struct TempClampInfo {
    /// Requested evaluation temperature in Kelvin.
    pub requested_k: f64,
    /// Clamped evaluation temperature actually used in Kelvin.
    pub clamped_k: f64,
    /// Valid temperature interval `[T_min, T_max]` in Kelvin.
    pub bounds_k: [f64; 2],
}

/// Result of evaluating a temperature-dependent quantity together with optional out-of-range clamping metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct CheckedTempEval<T> {
    pub value: T,
    pub clamp: Option<TempClampInfo>,
}

/// Optional tensor-component label (`"perp"` / `"par"` or `None` for scalar) paired with clamping metadata.
pub type ComponentClamp = (Option<&'static str>, TempClampInfo);

/// Optical crystal axis index (`0`, `1`, `2`) paired with thermo-optic clamping metadata.
pub type AxisClamp = (usize, TempClampInfo);

/// Clamps `temperature_k` to `[lo, hi]` with sub-ULP rounding slack on unit conversions, returning clamping metadata when outside bounds.
pub(crate) fn clamp_temp_range(temperature_k: f64, lo: f64, hi: f64) -> (f64, Option<TempClampInfo>) {
    if temperature_k.is_nan() {
        return (temperature_k, None);
    }
    let tol_lo = 1e-12 * lo.abs().max(1.0);
    let tol_hi = 1e-12 * hi.abs().max(1.0);
    let clamped_k = temperature_k.clamp(lo, hi);
    let clamp = if temperature_k < lo - tol_lo || temperature_k > hi + tol_hi {
        Some(TempClampInfo { requested_k: temperature_k, clamped_k, bounds_k: [lo, hi] })
    } else {
        None
    };
    (clamped_k, clamp)
}

/// A temperature-dependent scalar property curve in SI units.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "form", rename_all = "snake_case")]
pub enum TempCurve {
    Constant {
        value: f64,
    },
    #[serde(rename = "tabulated_T")]
    Tabulated {
        #[serde(rename = "T")]
        t: Vec<f64>,
        values: Vec<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        t_ref: Option<f64>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        source: Vec<Source>,
    },
    Log10Poly {
        #[serde(rename = "T_range")]
        t_range: [f64; 2],
        coeffs: Vec<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        t_ref: Option<f64>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        source: Vec<Source>,
    },
    EinsteinShift {
        ref_value: f64,
        #[serde(default = "default_t_ref")]
        t_ref: f64,
        theta_e: f64,
        alpha: f64,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        source: Vec<Source>,
    },
}

impl TempCurve {
    /// Valid temperature interval `[T_min, T_max]` in Kelvin when the curve is calibrated on a finite interval.
    pub fn t_range(&self) -> Option<[f64; 2]> {
        match self {
            Self::Constant { .. } | Self::EinsteinShift { .. } => None,
            Self::Tabulated { t, .. } => (t.len() >= 2).then(|| [t[0], t[t.len() - 1]]),
            Self::Log10Poly { t_range, .. } => {
                let lo = t_range[0].max(1e-6);
                let hi = t_range[1].max(lo);
                Some([lo, hi])
            }
        }
    }

    /// Citations attached directly to this temperature curve.
    pub fn sources(&self) -> &[Source] {
        match self {
            Self::Constant { .. } => &[],
            Self::Tabulated { source, .. }
            | Self::Log10Poly { source, .. }
            | Self::EinsteinShift { source, .. } => source,
        }
    }

    /// Evaluate the scalar quantity (in SI units) at `temperature_k` [K], returning clamping metadata when outside `T_range` / tabulated bounds.
    pub fn eval_checked(&self, temperature_k: f64) -> CheckedTempEval<f64> {
        match self {
            Self::Constant { value } => CheckedTempEval { value: *value, clamp: None },
            Self::Tabulated { t, values, .. } => {
                if t.is_empty() || values.is_empty() {
                    return CheckedTempEval { value: 0.0, clamp: None };
                }
                if temperature_k.is_nan() {
                    return CheckedTempEval { value: f64::NAN, clamp: None };
                }
                let last = t.len().min(values.len()) - 1;
                let (tc, clamp) = clamp_temp_range(temperature_k, t[0], t[last]);
                let value = if tc <= t[0] {
                    values[0]
                } else if tc >= t[last] {
                    values[last]
                } else {
                    let idx = t[..=last].partition_point(|&ti| ti <= tc);
                    let i0 = idx.saturating_sub(1).min(last.saturating_sub(1));
                    let i1 = (i0 + 1).min(last);
                    let dt = t[i1] - t[i0];
                    if dt <= 0.0 {
                        values[i0]
                    } else {
                        let w = (tc - t[i0]) / dt;
                        values[i0] + w * (values[i1] - values[i0])
                    }
                };
                CheckedTempEval { value, clamp }
            }
            Self::Log10Poly { t_range, coeffs, .. } => {
                if coeffs.is_empty() {
                    return CheckedTempEval { value: 0.0, clamp: None };
                }
                if temperature_k.is_nan() {
                    return CheckedTempEval { value: f64::NAN, clamp: None };
                }
                let lo = t_range[0].max(1e-6);
                let hi = t_range[1].max(lo);
                let (tc, clamp) = clamp_temp_range(temperature_k, lo, hi);
                let l = tc.log10();
                let mut acc = 0.0;
                for &c in coeffs.iter().rev() {
                    acc = acc * l + c;
                }
                CheckedTempEval { value: 10.0_f64.powf(acc), clamp }
            }
            Self::EinsteinShift {
                ref_value,
                t_ref,
                theta_e,
                alpha,
                ..
            } => {
                let g_t = einstein_thermal_coordinate(temperature_k, *theta_e);
                let g_ref = einstein_thermal_coordinate(*t_ref, *theta_e);
                CheckedTempEval { value: ref_value + alpha * (g_t - g_ref), clamp: None }
            }
        }
    }

    /// Evaluate the scalar quantity (in SI units) at `temperature_k` [K], refusing out-of-range clamping.
    pub fn eval_strict(&self, temperature_k: f64) -> Result<f64, TempClampInfo> {
        let checked = self.eval_checked(temperature_k);
        match checked.clamp {
            Some(info) => Err(info),
            None => Ok(checked.value),
        }
    }

    /// Evaluate the scalar quantity (in SI units) at `temperature_k` [K].
    pub fn eval(&self, temperature_k: f64) -> f64 {
        self.eval_checked(temperature_k).value
    }

    pub fn at_ref(&self) -> f64 {
        let t_ref = match self {
            Self::Constant { value } => return *value,
            Self::Tabulated { t_ref, .. } | Self::Log10Poly { t_ref, .. } => t_ref.unwrap_or(293.15),
            Self::EinsteinShift { t_ref, .. } => *t_ref,
        };
        self.eval(t_ref)
    }

    pub fn read(s: &mut Section, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<Self>, MaterialError> {
        let Some(v) = s.value(key) else {
            return Ok(None);
        };
        if let toml::Value::Table(tbl) = v {
            if tbl.contains_key("form") || tbl.contains_key("T") || tbl.contains_key("T_range") || tbl.contains_key("theta_e") {
                let path = if s.path.is_empty() { key.to_string() } else { format!("{}.{key}", s.path) };
                let mut sub = Section::new(&s.material, &path, tbl);
                return Self::read_table(&mut sub, dim, default_unit).map(Some);
            }
        }
        let unit = s.string(&format!("{key}_unit"))?.or_else(|| default_unit.map(str::to_string));
        let Some(val) = Self::parse_scalar_val(s, key, v, dim, unit.as_deref())? else {
            return Err(s.err(key, "must be a number, a quantity string, or a temperature-curve table"));
        };
        Ok(Some(Self::Constant { value: val }))
    }

    fn parse_scalar_val(s: &Section, key: &str, v: &toml::Value, dim: Dim, unit: Option<&str>) -> Result<Option<f64>, MaterialError> {
        match v {
            toml::Value::Float(f) => {
                if !f.is_finite() {
                    return Err(s.err(key, "must be finite"));
                }
                if dim == Dim::Dimensionless {
                    Ok(Some(*f))
                } else if let Some(u) = unit {
                    units::to_si(*f, u, dim).map(Some).map_err(|m| s.err(key, m))
                } else {
                    Err(s.err(key, "needs a unit"))
                }
            }
            toml::Value::Integer(i) => Self::parse_scalar_val(s, key, &toml::Value::Float(*i as f64), dim, unit),
            toml::Value::String(text) => {
                let (x, u) = units::split_quantity(text).map_err(|m| s.err(key, m))?;
                let eff_u = if u.is_empty() { unit } else { Some(u.as_str()) };
                Self::parse_scalar_val(s, key, &toml::Value::Float(x), dim, eff_u)
            }
            toml::Value::Table(t) if t.contains_key("value") => {
                let u = t.get("unit").and_then(|u| u.as_str()).or(unit);
                Self::parse_scalar_val(s, key, &t["value"], dim, u)
            }
            _ => Ok(None),
        }
    }

    fn read_table(s: &mut Section, dim: Dim, default_unit: Option<&str>) -> Result<Self, MaterialError> {
        let form = s
            .string("form")?
            .ok_or_else(|| s.err("form", "must specify form = \"tabulated_T\", \"log10_poly\", or \"einstein_shift\""))?;
        let source = read_sources(s, "source")?;
        let unit_str = s.string("unit")?.or_else(|| s.string("values_unit").ok().flatten()).or_else(|| default_unit.map(str::to_string));
        match form.as_str() {
            "tabulated_T" | "tabulated" => {
                let t_unit = s.string("T_unit")?.unwrap_or_else(|| "K".to_string());
                let t_raw = s.numbers("T")?.ok_or_else(|| s.err("T", "tabulated_T requires a temperature array T = [...]"))?;
                let mut t = Vec::with_capacity(t_raw.len());
                for (idx, x) in t_raw.into_iter().enumerate() {
                    let tk = units::to_si(x, &t_unit, Dim::Temperature).map_err(|m| s.err("T", m))?;
                    if idx > 0 && tk <= t[idx - 1] {
                        return Err(s.err("T", "temperatures in T must be strictly increasing"));
                    }
                    t.push(tk);
                }
                if t.len() < 2 {
                    return Err(s.err("T", "tabulated_T requires at least 2 temperature points"));
                }
                let val_unit = s.string("values_unit")?.or(unit_str);
                let v_raw = s.value("values").ok_or_else(|| s.err("values", "tabulated_T requires a values array"))?;
                let values = match v_raw {
                    toml::Value::Array(arr) => {
                        let mut out = Vec::with_capacity(arr.len());
                        for elem in arr {
                            let y = Self::parse_scalar_val(s, "values", elem, dim, val_unit.as_deref())?
                                .ok_or_else(|| s.err("values", "each entry must be a scalar quantity"))?;
                            out.push(y);
                        }
                        out
                    }
                    _ => return Err(s.err("values", "must be an array matching T")),
                };
                if values.len() != t.len() {
                    return Err(s.err("values", format!("length {} must match T length {}", values.len(), t.len())));
                }
                let t_ref = s.quantity("t_ref", Dim::Temperature, Some("K"))?;
                s.finish(&["form", "T", "T_unit", "values", "values_unit", "unit", "t_ref", "source"])?;
                Ok(Self::Tabulated { t, values, t_ref, source })
            }
            "log10_poly" => {
                let t_unit = s.string("T_unit")?.unwrap_or_else(|| "K".to_string());
                let tr = s
                    .numbers("T_range")?
                    .ok_or_else(|| s.err("T_range", "log10_poly requires T_range = [T_min, T_max]"))?;
                if tr.len() != 2 {
                    return Err(s.err("T_range", "T_range must have 2 elements [T_min, T_max]"));
                }
                let t0 = units::to_si(tr[0], &t_unit, Dim::Temperature).map_err(|m| s.err("T_range", m))?;
                let t1 = units::to_si(tr[1], &t_unit, Dim::Temperature).map_err(|m| s.err("T_range", m))?;
                if t0 <= 0.0 || t1 <= t0 {
                    return Err(s.err("T_range", "T_range must satisfy 0 < T_min < T_max"));
                }
                let mut coeffs = s
                    .numbers("coeffs")?
                    .ok_or_else(|| s.err("coeffs", "log10_poly requires coeffs = [c0, c1, ...]"))?;
                if coeffs.is_empty() {
                    return Err(s.err("coeffs", "coeffs must not be empty"));
                }
                if dim != Dim::Dimensionless {
                    let scale = match unit_str.as_deref() {
                        Some(u) => units::to_si(1.0, u, dim).map_err(|m| s.err("unit", m))?,
                        None => 1.0,
                    };
                    if (scale - 1.0).abs() > 1e-12 {
                        coeffs[0] += scale.log10();
                    }
                }
                let t_ref = s.quantity("t_ref", Dim::Temperature, Some("K"))?;
                s.finish(&["form", "T_range", "T_unit", "coeffs", "unit", "values_unit", "t_ref", "source"])?;
                Ok(Self::Log10Poly {
                    t_range: [t0, t1],
                    coeffs,
                    t_ref,
                    source,
                })
            }
            "einstein_shift" => {
                let ref_value = s
                    .quantity("ref_value", dim, unit_str.as_deref())?
                    .ok_or_else(|| s.err("ref_value", "einstein_shift requires ref_value"))?;
                let t_ref = s.quantity("t_ref", Dim::Temperature, Some("K"))?.unwrap_or(293.15);
                let theta_e = s
                    .quantity("theta_e", Dim::Temperature, Some("K"))?
                    .ok_or_else(|| s.err("theta_e", "einstein_shift requires theta_e (Einstein temperature in K)"))?;
                if theta_e <= 0.0 {
                    return Err(s.err("theta_e", "theta_e must be positive"));
                }
                let alpha_raw = s
                    .number("alpha")?
                    .ok_or_else(|| s.err("alpha", "einstein_shift requires high-T slope alpha in SI units per K"))?;
                let scale = if dim == Dim::Dimensionless {
                    1.0
                } else if let Some(u) = unit_str.as_deref() {
                    units::to_si(1.0, u, dim).map_err(|m| s.err("unit", m))?
                } else {
                    1.0
                };
                s.finish(&["form", "ref_value", "unit", "values_unit", "t_ref", "theta_e", "alpha", "source"])?;
                Ok(Self::EinsteinShift {
                    ref_value,
                    t_ref,
                    theta_e,
                    alpha: alpha_raw * scale,
                    source,
                })
            }
            other => Err(s.err("form", format!("unknown temperature-curve form '{other}'; expected tabulated_T, log10_poly or einstein_shift"))),
        }
    }

    pub fn to_toml(&self, unit: Option<&str>) -> toml::Value {
        let mut tbl = toml::Table::new();
        match self {
            Self::Constant { value } => {
                if let Some(u) = unit {
                    return toml::Value::String(format!("{value} {u}"));
                }
                return toml::Value::Float(*value);
            }
            Self::Tabulated { t, values, t_ref, source } => {
                tbl.insert("form".into(), "tabulated_T".into());
                tbl.insert("T".into(), toml::Value::Array(t.iter().copied().map(toml::Value::Float).collect()));
                tbl.insert("values".into(), toml::Value::Array(values.iter().copied().map(toml::Value::Float).collect()));
                if let Some(u) = unit {
                    tbl.insert("unit".into(), u.into());
                }
                if let Some(tr) = t_ref {
                    tbl.insert("t_ref".into(), toml::Value::Float(*tr));
                }
                if let Some(v) = write_sources(source) {
                    tbl.insert("source".into(), v);
                }
            }
            Self::Log10Poly { t_range, coeffs, t_ref, source } => {
                tbl.insert("form".into(), "log10_poly".into());
                tbl.insert(
                    "T_range".into(),
                    toml::Value::Array(vec![toml::Value::Float(t_range[0]), toml::Value::Float(t_range[1])]),
                );
                tbl.insert("coeffs".into(), toml::Value::Array(coeffs.iter().copied().map(toml::Value::Float).collect()));
                if let Some(u) = unit {
                    tbl.insert("unit".into(), u.into());
                }
                if let Some(tr) = t_ref {
                    tbl.insert("t_ref".into(), toml::Value::Float(*tr));
                }
                if let Some(v) = write_sources(source) {
                    tbl.insert("source".into(), v);
                }
            }
            Self::EinsteinShift {
                ref_value,
                t_ref,
                theta_e,
                alpha,
                source,
            } => {
                tbl.insert("form".into(), "einstein_shift".into());
                tbl.insert("ref_value".into(), toml::Value::Float(*ref_value));
                tbl.insert("t_ref".into(), toml::Value::Float(*t_ref));
                tbl.insert("theta_e".into(), toml::Value::Float(*theta_e));
                tbl.insert("alpha".into(), toml::Value::Float(*alpha));
                if let Some(u) = unit {
                    tbl.insert("unit".into(), u.into());
                }
                if let Some(v) = write_sources(source) {
                    tbl.insert("source".into(), v);
                }
            }
        }
        toml::Value::Table(tbl)
    }
}

/// Temperature-dependent scalar or uniaxial (perp/par) property in SI units.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum TempTensor {
    Scalar(TempCurve),
    Uniaxial { perp: TempCurve, par: TempCurve },
}

impl TempTensor {
    /// Returns temperature bounds `[T_min, T_max]` when the underlying scalar curve defines a finite range.
    pub fn t_range(&self) -> Option<[f64; 2]> {
        match self {
            Self::Scalar(c) => c.t_range(),
            Self::Uniaxial { perp, par } => match (perp.t_range(), par.t_range()) {
                (Some([a0, a1]), Some([b0, b1])) => Some([a0.max(b0), a1.min(b1)]),
                (Some(r), None) | (None, Some(r)) => Some(r),
                (None, None) => None,
            },
        }
    }

    /// All source entries attached to the underlying temperature curves.
    pub fn sources(&self) -> Vec<&Source> {
        match self {
            Self::Scalar(c) => c.sources().iter().collect(),
            Self::Uniaxial { perp, par } => perp.sources().iter().chain(par.sources().iter()).collect(),
        }
    }

    /// Evaluates the scalar/uniaxial tensor property at `temperature_k` [K], returning any clamping events with component label.
    pub fn eval_tensor_checked(&self, temperature_k: f64) -> (TensorProp, Vec<ComponentClamp>) {
        match self {
            Self::Scalar(c) => {
                let r = c.eval_checked(temperature_k);
                let clamps = r.clamp.into_iter().map(|info| (None, info)).collect();
                (TensorProp::Scalar(r.value), clamps)
            }
            Self::Uniaxial { perp, par } => {
                let rp = perp.eval_checked(temperature_k);
                let rpar = par.eval_checked(temperature_k);
                let mut clamps = Vec::new();
                if let Some(info) = rp.clamp {
                    clamps.push((Some("perp"), info));
                }
                if let Some(info) = rpar.clamp {
                    clamps.push((Some("par"), info));
                }
                (TensorProp::Uniaxial { perp: rp.value, par: rpar.value }, clamps)
            }
        }
    }

    /// Evaluates the isotropic / orientation-averaged scalar property at `temperature_k` [K], returning any clamping events.
    pub fn eval_scalar_checked(&self, temperature_k: f64) -> (f64, Vec<ComponentClamp>) {
        match self {
            Self::Scalar(c) => {
                let r = c.eval_checked(temperature_k);
                (r.value, r.clamp.into_iter().map(|info| (None, info)).collect())
            }
            Self::Uniaxial { perp, par } => {
                let rp = perp.eval_checked(temperature_k);
                let rpar = par.eval_checked(temperature_k);
                let mut clamps = Vec::new();
                if let Some(info) = rp.clamp {
                    clamps.push((Some("perp"), info));
                }
                if let Some(info) = rpar.clamp {
                    clamps.push((Some("par"), info));
                }
                ((2.0 * rp.value + rpar.value) / 3.0, clamps)
            }
        }
    }

    pub fn eval_tensor(&self, temperature_k: f64) -> TensorProp {
        self.eval_tensor_checked(temperature_k).0
    }

    pub fn eval_scalar(&self, temperature_k: f64) -> f64 {
        self.eval_scalar_checked(temperature_k).0
    }

    pub fn read(s: &mut Section, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<Self>, MaterialError> {
        let Some(v) = s.value(key) else {
            return Ok(None);
        };
        if let toml::Value::Table(tbl) = v {
            if tbl.contains_key("perp") && tbl.contains_key("par") && tbl.len() == 2 {
                let base_path = if s.path.is_empty() { key.to_string() } else { format!("{}.{key}", s.path) };
                let mut sub = Section::new(&s.material, &base_path, tbl);
                let unit = s.string(&format!("{key}_unit"))?.or_else(|| default_unit.map(str::to_string));
                let perp = TempCurve::read(&mut sub, "perp", dim, unit.as_deref())?
                    .ok_or_else(|| sub.err("perp", "uniaxial temperature curve requires 'perp'"))?;
                let par = TempCurve::read(&mut sub, "par", dim, unit.as_deref())?
                    .ok_or_else(|| sub.err("par", "uniaxial temperature curve requires 'par'"))?;
                sub.finish(&["perp", "par"])?;
                return Ok(Some(Self::Uniaxial { perp, par }));
            }
        }
        TempCurve::read(s, key, dim, default_unit).map(|o| o.map(Self::Scalar))
    }

    pub fn to_toml(&self, unit: Option<&str>) -> toml::Value {
        match self {
            Self::Scalar(c) => c.to_toml(unit),
            Self::Uniaxial { perp, par } => {
                let mut tbl = toml::Table::new();
                tbl.insert("perp".into(), perp.to_toml(unit));
                tbl.insert("par".into(), par.to_toml(unit));
                toml::Value::Table(tbl)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn einstein_saturates_flat_at_zero_kelvin_and_hits_ref() {
        let c = TempCurve::EinsteinShift {
            ref_value: 3.4757,
            t_ref: 293.15,
            theta_e: 400.0,
            alpha: 1.8e-4,
            source: Vec::new(),
        };
        assert!((c.eval(293.15) - 3.4757).abs() < 1e-12);
        // Near 0 K, slope dX/dT must approach zero (Third Law)
        let n_1k = c.eval(1.0);
        let n_10k = c.eval(10.0);
        assert!((n_10k - n_1k).abs() < 1e-18, "slope near 0 K should be exponentially flat: {} vs {}", n_10k, n_1k);
        assert!(c.eval(4.0) < c.eval(293.15));
    }

    #[test]
    fn nist_log10_poly_oxygen_free_copper_cp_shape() {
        // Sample 4..300 K NIST Cryogenic DB style log10 polynomial
        let toml_text = r#"
[cp]
form = "log10_poly"
T_range = [4.0, 300.0]
coeffs = [-1.91844, -0.15973, 8.61013, -18.9964, 21.9661, -12.7328, 3.54322, -0.3797]
unit = "J/(kg*K)"
"#;
        let tbl: toml::Table = toml_text.parse().unwrap();
        let mut s = Section::new("test", "", &tbl);
        let curve = TempCurve::read(&mut s, "cp", Dim::SpecificHeat, None).unwrap().unwrap();
        let cp_4 = curve.eval(4.0);
        let cp_300 = curve.eval(300.0);
        assert!(cp_4 > 0.01 && cp_4 < 1.0, "cp(4 K) = {cp_4}");
        assert!(cp_300 > 300.0 && cp_300 < 450.0, "cp(300 K) = {cp_300}");
        // Clamping beyond T_range keeps finite bounded outputs
        assert_eq!(curve.eval(0.5), cp_4);
        assert_eq!(curve.eval(500.0), cp_300);
    }

    #[test]
    fn centrosymmetric_point_groups_flagged_properly() {
        assert!(PointGroup::parse("m-3m").unwrap().forbids_piezo_and_pockels());
        assert!(PointGroup::parse("432").unwrap().forbids_piezo_and_pockels());
        assert!(!PointGroup::parse("432").unwrap().is_centrosymmetric());
        assert!(!PointGroup::parse("3m").unwrap().forbids_piezo_and_pockels());
        assert!(!PointGroup::parse("6mm").unwrap().forbids_piezo_and_pockels());
        assert!(!PointGroup::parse("-43m").unwrap().forbids_piezo_and_pockels());
    }
}
