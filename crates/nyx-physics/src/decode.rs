// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Decoding the problem envelope into what a kernel reads: typed parameters with the descriptor defaults, region names resolved to
//! domain numbers (`materials_by_label`), condition selectors resolved to boundary or domain numbers (`resolve_conditions`). Every physics
//! crate uses these, so that the envelope is read the same way in all of them (ADR 41).

use crate::cqed::{TlsInterfaceSpec, TlsInterfaces};
use crate::desc::{Applies, BcDesc, ParamDesc, ParamKind, PropDesc, PropKind};
use crate::error::Error;
use crate::physics::Ctx;
use crate::result::{DiagnosticSeverity, DiagnosticWarning};
use crate::spec::{BcSpec, MaterialAssign, ProblemSpec, PropValue};
use femcore::material::{MaterialSpec, Source, TempClampInfo};
use femcore::props::Materials;
use serde_json::Value;
use std::collections::BTreeSet;

/// Typed access to the parameters of a study, a condition or the physics: returns the given value or the descriptor default,
/// rejecting unknown names or mismatched JSON types cleanly with [`Error::Invalid`].
pub struct Params<'a> {
    value: &'a Value,
    descs: &'a [ParamDesc],
}

impl<'a> Params<'a> {
    pub fn new(value: &'a Value, descs: &'a [ParamDesc]) -> Params<'a> {
        Params { value, descs }
    }

    fn desc(&self, name: &str) -> Result<&ParamDesc, Error> {
        self.descs.iter().find(|d| d.name == name).ok_or_else(|| Error::invalid(format!("unknown parameter '{name}'")))
    }

    fn given(&self, name: &str) -> Option<&Value> {
        self.value.get(name).filter(|v| !v.is_null())
    }

    /// True if the parameter key was explicitly passed with JSON `null`.
    pub fn is_null(&self, name: &str) -> bool {
        self.value.get(name).is_some_and(Value::is_null)
    }

    /// Optional numeric parameter: returns `None` when explicitly passed as `null`,
    /// otherwise the given value or the descriptor default.
    pub fn optional_number(&self, name: &str) -> Result<Option<f64>, Error> {
        if self.is_null(name) {
            Ok(None)
        } else {
            self.number(name).map(Some)
        }
    }

    /// Whether the parameter was given with a non-null value in the problem (not defaulted).
    pub fn has(&self, name: &str) -> bool {
        self.given(name).is_some()
    }

    pub fn number(&self, name: &str) -> Result<f64, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v.as_f64().ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a number"))),
            (None, ParamKind::Number { default, .. }) => Ok(*default),
            _ => Err(Error::invalid(format!("parameter '{name}' is not a number"))),
        }
    }

    pub fn integer(&self, name: &str) -> Result<i64, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v.as_i64().ok_or_else(|| Error::invalid(format!("parameter '{name}' must be an integer"))),
            (None, ParamKind::Integer { default, .. }) => Ok(*default),
            _ => Err(Error::invalid(format!("parameter '{name}' is not an integer"))),
        }
    }

    pub fn boolean(&self, name: &str) -> Result<bool, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v.as_bool().ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a boolean"))),
            (None, ParamKind::Bool { default }) => Ok(*default),
            _ => Err(Error::invalid(format!("parameter '{name}' is not a boolean"))),
        }
    }

    /// Text or a choice.
    pub fn text(&self, name: &str) -> Result<String, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => {
                v.as_str().map(str::to_string).ok_or_else(|| Error::invalid(format!("parameter '{name}' must be text")))
            }
            (None, ParamKind::Text { default }) | (None, ParamKind::Choice { default, .. }) => Ok(default.clone()),
            _ => Err(Error::invalid(format!("parameter '{name}' is not text"))),
        }
    }

    pub fn vector(&self, name: &str) -> Result<Vec<f64>, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v
                .as_array()
                .ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a vector of numbers")))?
                .iter()
                .map(|item| {
                    item.as_f64().ok_or_else(|| Error::invalid(format!("parameter '{name}' entries must be numbers")))
                })
                .collect(),
            (None, ParamKind::Vector { default, .. }) => Ok(default.clone()),
            _ => Err(Error::invalid(format!("parameter '{name}' is not a vector"))),
        }
    }

    pub fn names(&self, name: &str) -> Result<Vec<String>, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v
                .as_array()
                .ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a list of names")))?
                .iter()
                .map(|item| {
                    item.as_str()
                        .map(str::to_string)
                        .ok_or_else(|| Error::invalid(format!("parameter '{name}' entries must be strings")))
                })
                .collect(),
            (None, ParamKind::Names { default }) => Ok(default.clone()),
            _ => Err(Error::invalid(format!("parameter '{name}' is not a list of names"))),
        }
    }

    pub fn numbers(&self, name: &str) -> Result<Vec<f64>, Error> {
        let desc = self.desc(name)?;
        match (self.given(name), &desc.kind) {
            (Some(v), _) => v
                .as_array()
                .ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a list of numbers")))?
                .iter()
                .map(|item| {
                    item.as_f64().ok_or_else(|| Error::invalid(format!("parameter '{name}' entries must be numbers")))
                })
                .collect(),
            (None, ParamKind::Numbers { default }) | (None, ParamKind::Vector { default, .. }) => Ok(default.clone()),
            _ => Err(Error::invalid(format!("parameter '{name}' is not a list of numbers"))),
        }
    }

    /// Optional number-list parameter: returns `None` when explicitly passed as `null` or when the resolved list is empty.
    pub fn optional_numbers(&self, name: &str) -> Result<Option<Vec<f64>>, Error> {
        if self.is_null(name) {
            Ok(None)
        } else {
            let v = self.numbers(name)?;
            Ok((!v.is_empty()).then_some(v))
        }
    }
}

/// Out-of-range temperature clamping recorded while evaluating an active material table.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialClampEvent {
    pub region: String,
    pub material: String,
    pub section: String,
    pub clamp: TempClampInfo,
}

impl MaterialClampEvent {
    /// Converts this clamping event into a [`DiagnosticWarning`] with code `"temperature_clamped"`.
    pub fn to_warning(&self) -> DiagnosticWarning {
        DiagnosticWarning {
            code: "temperature_clamped".into(),
            severity: DiagnosticSeverity::Warning,
            message: format!(
                "region '{}' (material '{}' [{}]): temperature {:.4} K is outside calibrated range [{:.4}, {:.4}] K; clamped to {:.4} K",
                self.region,
                self.material,
                self.section,
                self.clamp.requested_k,
                self.clamp.bounds_k[0],
                self.clamp.bounds_k[1],
                self.clamp.clamped_k,
            ),
        }
    }

    fn describe_out_of_range(&self) -> String {
        format!(
            "region '{}' (material '{}' [{}]): T = {:.4} K outside [{:.4}, {:.4}] K",
            self.region,
            self.material,
            self.section,
            self.clamp.requested_k,
            self.clamp.bounds_k[0],
            self.clamp.bounds_k[1],
        )
    }

    /// Formats a refusal error for `strict_material_range = true`.
    pub fn to_strict_error(&self) -> Error {
        Error::invalid(format!("strict_material_range: {}", self.describe_out_of_range()))
    }
}

/// Active-study material provenance uncertainty and temperature-clamping findings.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MaterialAudit {
    pub warnings: Vec<DiagnosticWarning>,
    pub clamps: Vec<MaterialClampEvent>,
}

impl MaterialAudit {
    /// Returns all warnings (both provenance uncertainty and `temperature_clamped`).
    pub fn all_warnings(&self) -> Vec<DiagnosticWarning> {
        let mut out = self.warnings.clone();
        out.extend(self.clamps.iter().map(MaterialClampEvent::to_warning));
        out
    }

    /// Returns the first strict refusal error if any active temperature evaluation was clamped.
    pub fn strict_error(&self) -> Option<Error> {
        if self.clamps.is_empty() {
            return None;
        }
        let details: Vec<String> = self.clamps.iter().map(MaterialClampEvent::describe_out_of_range).collect();
        Some(Error::invalid(format!(
            "strict_material_range refused out-of-range temperature evaluation: {}",
            details.join("; ")
        )))
    }
}

/// Resolves `spec.materials` per domain label and returns both [`Materials`] and the active-study [`MaterialAudit`].
pub fn materials_by_label_checked(
    spec: &ProblemSpec,
    ctx: &Ctx,
    props: &[PropDesc],
) -> Result<(Materials, MaterialAudit), Error> {
    let mut m = Materials::default();
    let is_quasi_static_electro =
        (spec.physics == "dfs" && spec.study.kind == "electrostatic") || spec.physics == "sms";
    let has_eps_r = props.iter().any(|d| d.name == "eps_r");
    let find_prop_kind = |name: &str| -> Option<PropKind> {
        if let Some(d) = props.iter().find(|d| d.name == name) {
            Some(if name == "loss_tangent" { PropKind::Tensor } else { d.kind })
        } else if has_eps_r && matches!(name, "loss_tangent" | "tan_delta" | "sigma" | "eps_r_im") {
            Some(PropKind::Tensor)
        } else {
            None
        }
    };
    let freq_hz = effective_freq_hz(spec, ctx, is_quasi_static_electro);
    let temp_k = effective_temp_k(spec, ctx);
    for (region, assign) in &spec.materials {
        let labels = ctx.resolver.domains(region)?;
        match assign {
            MaterialAssign::Props(map) => {
                for (name, v) in map {
                    let Some(prop_kind) = find_prop_kind(name) else {
                        return Err(Error::invalid(format!("region '{region}': unknown property '{name}'")));
                    };
                    if prop_kind == PropKind::Matrix3x6 {
                        return Err(Error::invalid(format!(
                            "region '{region}': property '{name}' must be a 3x6 matrix, not a scalar number"
                        )));
                    }
                    if !v.is_finite() {
                        return Err(Error::invalid(format!(
                            "region '{region}': property '{name}' is not a finite number"
                        )));
                    }
                    if name == "loss_tangent" && *v < 0.0 {
                        return Err(Error::invalid(format!(
                            "region '{region}': property 'loss_tangent' must be non-negative, got {v}"
                        )));
                    }
                    for &l in &labels {
                        m.set(l, name, *v);
                    }
                }
                if let Some(&td) = map.get("loss_tangent") {
                    if td > 0.0 && !map.contains_key("eps_r_im") {
                        let eps_re = map
                            .get("eps_r")
                            .copied()
                            .or_else(|| props.iter().find(|d| d.name == "eps_r").and_then(|d| d.default))
                            .unwrap_or(1.0);
                        let eps_im = eps_re * td;
                        if eps_im > 0.0 {
                            for &l in &labels {
                                m.set(l, "eps_r_im", eps_im);
                            }
                        }
                    }
                }
            }
            MaterialAssign::PropsTensor(map) => {
                let mut inline_td_tensor: Option<femcore::material::Tensor> = None;
                for (name, v) in map {
                    let Some(prop_kind) = find_prop_kind(name) else {
                        return Err(Error::invalid(format!("region '{region}': unknown property '{name}'")));
                    };
                    match v {
                        PropValue::Scalar(s) => {
                            if prop_kind == PropKind::Matrix3x6 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a 3x6 matrix, not a scalar number"
                                )));
                            }
                            if !s.is_finite() {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' is not a finite number"
                                )));
                            }
                            if name == "loss_tangent" && *s < 0.0 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property 'loss_tangent' must be non-negative, got {s}"
                                )));
                            }
                            for &l in &labels {
                                m.set(l, name, *s);
                            }
                        }
                        PropValue::Diag([x, y, z]) => {
                            if prop_kind != PropKind::Tensor {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                                )));
                            }
                            if !x.is_finite() || !y.is_finite() || !z.is_finite() {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite tensor component"
                                )));
                            }
                            if name == "loss_tangent" && (*x < 0.0 || *y < 0.0 || *z < 0.0) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property 'loss_tangent' must be non-negative in every direction"
                                )));
                            }
                            let diag_m3 = [[*x, 0.0, 0.0], [0.0, *y, 0.0], [0.0, 0.0, *z]];
                            if is_isotropic_m3(&diag_m3) {
                                for &l in &labels {
                                    m.set(l, name, *x);
                                }
                            } else if name != "loss_tangent" {
                                let tensor = femcore::material::Tensor::Diagonal { x: *x, y: *y, z: *z };
                                if name == "tan_delta" {
                                    inline_td_tensor = Some(tensor);
                                } else {
                                    set_inline_tensor(&mut m, region, &labels, name, tensor)?;
                                }
                            }
                        }
                        PropValue::Sym6([xx, yy, zz, xy, xz, yz]) => {
                            if prop_kind != PropKind::Tensor {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                                )));
                            }
                            if [*xx, *yy, *zz, *xy, *xz, *yz].iter().any(|v| !v.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite tensor component"
                                )));
                            }
                            let sym_m3 = [[*xx, *xy, *xz], [*xy, *yy, *yz], [*xz, *yz, *zz]];
                            if name == "loss_tangent" {
                                femcore::material::tensor::TensorProp::Symmetric([*xx, *yy, *zz, *yz, *xz, *xy])
                                    .check(femcore::material::tensor::Requirement::NonNegative, "loss_tangent")
                                    .map_err(|msg| Error::invalid(format!("region '{region}': {msg}")))?;
                            }
                            if is_isotropic_m3(&sym_m3) {
                                for &l in &labels {
                                    m.set(l, name, *xx);
                                }
                            } else if name != "loss_tangent" {
                                let tensor = femcore::material::Tensor::Symmetric {
                                    xx: *xx,
                                    yy: *yy,
                                    zz: *zz,
                                    xy: *xy,
                                    yz: *yz,
                                    xz: *xz,
                                };
                                if name == "tan_delta" {
                                    inline_td_tensor = Some(tensor);
                                } else {
                                    set_inline_tensor(&mut m, region, &labels, name, tensor)?;
                                }
                            }
                        }
                        PropValue::Matrix(mat) => {
                            if prop_kind != PropKind::Tensor {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                                )));
                            }
                            if mat.iter().flatten().any(|x| !x.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite tensor component"
                                )));
                            }
                            let max_abs = mat.iter().flatten().fold(0.0_f64, |acc, x| acc.max(x.abs()));
                            let tol = 1e-12 * max_abs.max(1.0);
                            if (mat[0][1] - mat[1][0]).abs() > tol
                                || (mat[0][2] - mat[2][0]).abs() > tol
                                || (mat[1][2] - mat[2][1]).abs() > tol
                            {
                                return Err(Error::invalid(format!(
                                    "region '{region}': tensor property '{name}' must be symmetric"
                                )));
                            }
                            if name == "loss_tangent" {
                                femcore::material::tensor::TensorProp::Symmetric([
                                    mat[0][0],
                                    mat[1][1],
                                    mat[2][2],
                                    0.5 * (mat[1][2] + mat[2][1]),
                                    0.5 * (mat[0][2] + mat[2][0]),
                                    0.5 * (mat[0][1] + mat[1][0]),
                                ])
                                .check(femcore::material::tensor::Requirement::NonNegative, "loss_tangent")
                                .map_err(|msg| Error::invalid(format!("region '{region}': {msg}")))?;
                            }
                            if is_isotropic_m3(mat) {
                                for &l in &labels {
                                    m.set(l, name, mat[0][0]);
                                }
                            } else if name != "loss_tangent" {
                                let tensor = femcore::material::Tensor::Symmetric {
                                    xx: mat[0][0],
                                    yy: mat[1][1],
                                    zz: mat[2][2],
                                    xy: 0.5 * (mat[0][1] + mat[1][0]),
                                    yz: 0.5 * (mat[1][2] + mat[2][1]),
                                    xz: 0.5 * (mat[0][2] + mat[2][0]),
                                };
                                if name == "tan_delta" {
                                    inline_td_tensor = Some(tensor);
                                } else {
                                    set_inline_tensor(&mut m, region, &labels, name, tensor)?;
                                }
                            }
                        }
                        PropValue::Matrix3x6(mat) => {
                            if prop_kind != PropKind::Matrix3x6 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' cannot use a 3x6 matrix"
                                )));
                            }
                            if mat.iter().flatten().any(|x| !x.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite 3x6 matrix component"
                                )));
                            }
                            if name == "piezo_e" {
                                for &l in &labels {
                                    m.aniso.entry(l).or_default().piezo_e = Some(*mat);
                                }
                            } else {
                                return Err(Error::invalid(format!(
                                    "region '{region}': unsupported 3x6 property '{name}'"
                                )));
                            }
                        }
                    }
                }
                if let Some(td_val) = map.get("loss_tangent") {
                    if !map.contains_key("eps_r_im") {
                        let eps_re_m3 = match map.get("eps_r") {
                            Some(ev) => prop_value_to_m3(ev),
                            None => {
                                let d = props.iter().find(|d| d.name == "eps_r").and_then(|d| d.default).unwrap_or(1.0);
                                femcore::material::iso_tensor(d)
                            }
                        };
                        let td_m3 = prop_value_to_m3(td_val);
                        let eps_im_m3 = symmetrize_m3(&mul_m3(&eps_re_m3, &td_m3));
                        let has_im = eps_im_m3.iter().flatten().any(|&v| v.abs() > 0.0);
                        if has_im {
                            if is_isotropic_m3(&eps_im_m3) {
                                for &l in &labels {
                                    m.set(l, "eps_r_im", eps_im_m3[0][0]);
                                    if m.aniso.get(&l).is_some_and(|a| a.permittivity.is_some()) {
                                        m.aniso.entry(l).or_default().permittivity_im =
                                            Some(m3_to_sym_tensor(&eps_im_m3));
                                    }
                                }
                            } else {
                                set_inline_tensor(&mut m, region, &labels, "eps_r_im", m3_to_sym_tensor(&eps_im_m3))?;
                            }
                        }
                    }
                }
                apply_inline_tan_delta_tensors(&mut m, &labels, inline_td_tensor.as_ref());
            }
            MaterialAssign::Named(n) => named(
                &mut m,
                region,
                &labels,
                n,
                [0.0; 3],
                ctx,
                freq_hz,
                temp_k,
                props,
                is_quasi_static_electro,
                &spec.physics,
            )?,
            MaterialAssign::Oriented { material, euler_deg } => named(
                &mut m,
                region,
                &labels,
                material,
                *euler_deg,
                ctx,
                freq_hz,
                temp_k,
                props,
                is_quasi_static_electro,
                &spec.physics,
            )?,
        }
    }
    let temps: Vec<f64> = temp_k.into_iter().collect();
    let audit = audit_active_materials(spec, ctx, props, &temps)?;
    Ok((m, audit))
}

/// The material properties of every region as numbers per domain label: a region name or pattern (`sapphire*`) can match several
/// domains. A property outside `props` is refused for inline properties (`validate_spec` already does that; this keeps a kernel safe when
/// called directly). A named material comes from the library of `ctx` (`MaterialAssign::Named`, or `Oriented` with Bunge angles): the
/// scalar properties the physics lists are taken under their descriptor names, the crystal data (stiffness, tensors) as an `Anisotropy`
/// where the material is anisotropic or oriented, and `eps_r` from the optical model at `ctx.wavelength` where the material has one.
pub fn materials_by_label(spec: &ProblemSpec, ctx: &Ctx, props: &[PropDesc]) -> Result<Materials, Error> {
    materials_by_label_checked(spec, ctx, props).map(|(m, _)| m)
}

fn effective_freq_hz(spec: &ProblemSpec, ctx: &Ctx, is_quasi_static_electro: bool) -> Option<f64> {
    let spec_freq = spec
        .study
        .params
        .get("freq")
        .or_else(|| spec.study.params.get("frequency"))
        .and_then(Value::as_f64)
        .or_else(|| spec.study.params.get("f_search").and_then(Value::as_f64))
        .filter(|&f| f > 0.0 && f.is_finite());
    ctx.frequency
        .or_else(|| {
            (!is_quasi_static_electro)
                .then_some(ctx.wavelength)
                .flatten()
                .filter(|&l| l > 0.0)
                .map(|l| 299_792_458.0 / l)
        })
        .or(spec_freq)
}

fn effective_temp_k(spec: &ProblemSpec, ctx: &Ctx) -> Option<f64> {
    ctx.temperature
        .or_else(|| spec.params.get("temperature").and_then(Value::as_f64).filter(|&t| t > 0.0 && t.is_finite()))
}

fn uses_prop_in_study(spec: &ProblemSpec, props: &[PropDesc], name: &str) -> bool {
    let study = spec.study.kind.as_str();
    if spec.physics == "sms" && matches!(name, "rho" | "eta") {
        return props.iter().any(|d| d.name == name) && study != "static";
    }
    props
        .iter()
        .find(|d| d.name == name)
        .is_some_and(|d| d.needed_by.is_empty() || d.needed_by.iter().any(|s| s == study))
}

#[derive(Default)]
struct AuditCollector {
    seen_warnings: BTreeSet<(String, String, String, &'static str)>,
    seen_clamps: BTreeSet<(String, String, String, u64, u64, u64)>,
    audit: MaterialAudit,
}

impl AuditCollector {
    fn record_sources(&mut self, region: &str, material: &str, section: &str, primary: &[Source], fallback: &[Source]) {
        let sources = if primary.is_empty() { fallback } else { primary };
        for src in sources {
            if !src.flags_uncertainty() {
                continue;
            }
            if src.is_unverified() || !src.is_approximate() {
                self.push_warning(region, material, section, "material_unverified", &src.cite);
            }
            if src.is_approximate() {
                self.push_warning(region, material, section, "material_approximate", &src.cite);
            }
        }
    }

    fn push_warning(&mut self, region: &str, material: &str, section: &str, code: &'static str, cite: &str) {
        let key = (region.to_string(), material.to_string(), section.to_string(), code);
        if self.seen_warnings.insert(key) {
            self.audit.warnings.push(DiagnosticWarning {
                code: code.to_string(),
                severity: DiagnosticSeverity::Warning,
                message: format!("region '{region}' (material '{material}' [{section}]): {cite}"),
            });
        }
    }

    fn record_clamp(&mut self, region: &str, material: &str, section: impl Into<String>, clamp: TempClampInfo) {
        let section = section.into();
        let key = (
            region.to_string(),
            material.to_string(),
            section.clone(),
            clamp.requested_k.to_bits(),
            clamp.bounds_k[0].to_bits(),
            clamp.bounds_k[1].to_bits(),
        );
        if self.seen_clamps.insert(key) {
            self.audit.clamps.push(MaterialClampEvent {
                region: region.to_string(),
                material: material.to_string(),
                section,
                clamp,
            });
        }
    }
}

fn clamp_section(base: &str, axis_or_component: Option<&str>) -> String {
    match axis_or_component {
        Some(tag) => format!("{base} ({tag})"),
        None => base.to_string(),
    }
}

/// Audits active-study material provenance (`Source::flags_uncertainty()`) and out-of-range temperature clamping across `spec.materials` and boundary conditions.
pub fn audit_active_materials(
    spec: &ProblemSpec,
    ctx: &Ctx,
    props: &[PropDesc],
    temps_k: &[f64],
) -> Result<MaterialAudit, Error> {
    let mut collector = AuditCollector::default();
    let mut check_temps: Vec<f64> = Vec::new();
    for &t in effective_temp_k(spec, ctx).iter().chain(temps_k.iter()) {
        if !t.is_nan() && !check_temps.iter().any(|x| x.to_bits() == t.to_bits()) {
            check_temps.push(t);
        }
    }

    let is_quasi_static_electro =
        (spec.physics == "dfs" && spec.study.kind == "electrostatic") || spec.physics == "sms";
    let freq_hz = effective_freq_hz(spec, ctx, is_quasi_static_electro);
    let thermal_load_active = spec.physics != "sms"
        || matches!(spec.params.get("thermal").and_then(Value::as_str), Some("uniform" | "imported"));
    let has_elec_bc = spec.bcs.iter().any(|bc| {
        matches!(
            bc.kind.as_str(),
            "ground" | "potential" | "surface_charge" | "floating_potential" | "terminal" | "conductor"
        )
    });
    let piezo_active = spec.physics != "sms"
        || spec.study.kind == "driven"
        || spec.params.get("piezoelectric").and_then(Value::as_str) == Some("on")
        || (spec.params.get("piezoelectric").and_then(Value::as_str) != Some("off") && has_elec_bc);
    let wants_k = uses_prop_in_study(spec, props, "k");
    let wants_cp = uses_prop_in_study(spec, props, "Cp");
    let wants_rho = uses_prop_in_study(spec, props, "rho");
    let wants_alpha = uses_prop_in_study(spec, props, "alpha") && thermal_load_active;
    let wants_e = uses_prop_in_study(spec, props, "E");
    let wants_eta = uses_prop_in_study(spec, props, "eta");
    let wants_eps_r = uses_prop_in_study(spec, props, "eps_r") && piezo_active;
    let wants_mu_r = uses_prop_in_study(spec, props, "mu_r");
    let wants_sigma = uses_prop_in_study(spec, props, "sigma")
        && (spec.physics == "dfs" || !wants_eps_r)
        && (!wants_eps_r || freq_hz.is_some_and(|f| f > 0.0));
    let tls_active = (spec.physics == "dfs" && spec.study.kind == "electrostatic")
        || (spec.physics == "ems" && matches!(spec.study.kind.as_str(), "cavity" | "driven"));

    for (region, assign) in &spec.materials {
        let (mat_name, euler) = match assign {
            MaterialAssign::Named(n) => (n.as_str(), [0.0; 3]),
            MaterialAssign::Oriented { material, euler_deg } => (material.as_str(), *euler_deg),
            MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_) => continue,
        };
        let Some(lib) = ctx.library.as_ref() else { continue };
        let Ok(mat) = lib.require(mat_name) else { continue };
        let mat: &MaterialSpec = mat;

        if wants_rho {
            collector.record_sources(region, mat_name, "density", &mat.source, &[]);
        }
        if wants_e {
            if let Some(el) = &mat.elastic {
                collector.record_sources(region, mat_name, "elastic", &el.source, &mat.source);
            } else {
                collector.record_sources(region, mat_name, "elastic", &[], &mat.source);
            }
        }
        if wants_eta {
            if let Some(ve) = &mat.viscoelastic {
                collector.record_sources(region, mat_name, "viscoelastic", &ve.source, &mat.source);
            }
        }
        if spec.physics == "sms" && piezo_active && uses_prop_in_study(spec, props, "piezo_e") {
            if let Some(pz) = &mat.piezoelectric {
                collector.record_sources(region, mat_name, "piezoelectric", &pz.source, &mat.source);
            }
        }
        if wants_k {
            collector.record_sources(region, mat_name, "thermal", &mat.thermal.source, &mat.source);
            if let Some(curve) = &mat.thermal.conductivity_t {
                let curve_srcs: Vec<Source> = curve.sources().into_iter().cloned().collect();
                collector.record_sources(region, mat_name, "thermal.conductivity_T", &curve_srcs, &[]);
                for &t in &check_temps {
                    if let Some((_, clamps)) = mat.thermal_conductivity_at_checked(Some(t)) {
                        for (comp, info) in clamps {
                            collector.record_clamp(
                                region,
                                mat_name,
                                clamp_section("thermal.conductivity_T", comp),
                                info,
                            );
                        }
                    }
                }
            }
        }
        if wants_cp {
            collector.record_sources(region, mat_name, "thermal", &mat.thermal.source, &mat.source);
            if let Some(curve) = &mat.thermal.specific_heat_t {
                collector.record_sources(region, mat_name, "thermal.specific_heat_T", curve.sources(), &[]);
                for &t in &check_temps {
                    if let Some((_, Some(info))) = mat.specific_heat_at_checked(Some(t)) {
                        collector.record_clamp(region, mat_name, "thermal.specific_heat_T", info);
                    }
                }
            }
        }
        if wants_alpha {
            collector.record_sources(region, mat_name, "thermal", &mat.thermal.source, &mat.source);
            if let Some(curve) = &mat.thermal.expansion_t {
                let curve_srcs: Vec<Source> = curve.sources().into_iter().cloned().collect();
                collector.record_sources(region, mat_name, "thermal.expansion_T", &curve_srcs, &[]);
                for &t in &check_temps {
                    if let Some((_, clamps)) = mat.thermal_expansion_at_checked(Some(t)) {
                        for (comp, info) in clamps {
                            collector.record_clamp(region, mat_name, clamp_section("thermal.expansion_T", comp), info);
                        }
                    }
                }
            }
        }

        let optical_active = wants_eps_r
            && !is_quasi_static_electro
            && mat.optical.is_some()
            && ctx.wavelength.is_some_and(|l| l <= 3e-4);
        if optical_active {
            if let Some(opt) = &mat.optical {
                let mut opt_sources = opt.source.clone();
                for ax in opt.axes.models() {
                    opt_sources.extend(ax.sources().iter().cloned());
                }
                collector.record_sources(region, mat_name, "optical", &opt_sources, &mat.source);
                if let Some(lam) = ctx.wavelength {
                    for &t in &check_temps {
                        if let Ok((_, clamps)) = mat.optical_permittivity_checked(lam, Some(t), euler) {
                            for (axis, info) in clamps {
                                let label = match &opt.axes {
                                    femcore::material::optics::Axes::Isotropic(_) => {
                                        if axis != 0 {
                                            continue;
                                        }
                                        "optical".to_string()
                                    }
                                    femcore::material::optics::Axes::Uniaxial { .. } => {
                                        if axis == 1 {
                                            continue;
                                        }
                                        if axis == 2 {
                                            "optical.extraordinary".to_string()
                                        } else {
                                            "optical.ordinary".to_string()
                                        }
                                    }
                                    femcore::material::optics::Axes::Biaxial(_) => {
                                        const AX: [&str; 3] = ["optical.x", "optical.y", "optical.z"];
                                        AX[axis.min(2)].to_string()
                                    }
                                };
                                collector.record_clamp(region, mat_name, label, info);
                            }
                        }
                    }
                }
            }
            if mat.electrical.permeability.is_some() {
                collector.record_sources(region, mat_name, "electrical", &mat.electrical.source, &mat.source);
            }
        } else if wants_eps_r {
            if let Some(dl) = &mat.dielectric {
                collector.record_sources(region, mat_name, "dielectric", &dl.source, &mat.source);
                if let Some(curve) = &dl.eps_t {
                    let curve_srcs: Vec<Source> = curve.sources().into_iter().cloned().collect();
                    collector.record_sources(region, mat_name, "dielectric.eps_T", &curve_srcs, &[]);
                    for &t in &check_temps {
                        if let Some((_, clamps)) =
                            mat.complex_microwave_permittivity_checked(freq_hz, Some(t), None, false, euler)
                        {
                            for (comp, info) in clamps {
                                collector.record_clamp(region, mat_name, clamp_section("dielectric.eps_T", comp), info);
                            }
                        }
                    }
                }
            } else if mat.electrical.permittivity.is_some() || mat.electrical.loss_tangent.is_some() {
                collector.record_sources(region, mat_name, "electrical", &mat.electrical.source, &mat.source);
            } else {
                collector.record_sources(region, mat_name, "dielectric", &[], &mat.source);
            }
            if (!is_quasi_static_electro && mat.electrical.permeability.is_some())
                || (wants_sigma && mat.electrical.conductivity.is_some())
            {
                collector.record_sources(region, mat_name, "electrical", &mat.electrical.source, &mat.source);
            }
        } else if (wants_mu_r || wants_sigma)
            && (mat.electrical.permeability.is_some()
                || mat.electrical.conductivity.is_some()
                || mat.electrical.permittivity.is_some())
        {
            collector.record_sources(region, mat_name, "electrical", &mat.electrical.source, &mat.source);
        }

        if tls_active {
            for (iface_name, layer) in [
                ("interfaces.ma", &mat.interfaces.ma),
                ("interfaces.ms", &mat.interfaces.ms),
                ("interfaces.sa", &mat.interfaces.sa),
            ] {
                if let Some(l) = layer {
                    collector.record_sources(region, mat_name, iface_name, &l.source, &mat.source);
                }
            }
        }
    }

    if let Some(lib) = ctx.library.as_ref() {
        for bc in &spec.bcs {
            let Some(mat_name) =
                bc.params.get("material").and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty())
            else {
                continue;
            };
            let Ok(mat) = lib.require(mat_name) else { continue };
            let scope = format!("bc:{}", bc.kind);
            if matches!(bc.kind.as_str(), "london_sheet" | "surface_impedance") {
                if let Some(sc) = &mat.superconducting {
                    collector.record_sources(&scope, mat_name, "superconducting", &sc.source, &mat.source);
                }
            }
            if tls_active {
                for (iface_name, layer) in [
                    ("interfaces.ma", &mat.interfaces.ma),
                    ("interfaces.ms", &mat.interfaces.ms),
                    ("interfaces.sa", &mat.interfaces.sa),
                ] {
                    if let Some(l) = layer {
                        collector.record_sources(&scope, mat_name, iface_name, &l.source, &mat.source);
                    }
                }
            }
        }
    }

    Ok(collector.audit)
}

/// Resolves each material assignment region key in `spec.materials` to its domain labels.
pub fn domain_names_by_label(spec: &ProblemSpec, ctx: &Ctx) -> Result<std::collections::BTreeMap<i32, String>, Error> {
    let mut out = std::collections::BTreeMap::new();
    for region in spec.materials.keys() {
        for l in ctx.resolver.domains(region)? {
            out.entry(l).or_insert_with(|| region.clone());
        }
    }
    Ok(out)
}

fn map_interface_layer(layer: &femcore::material::InterfaceLayer) -> TlsInterfaceSpec {
    TlsInterfaceSpec { thickness_m: layer.thickness, eps_r: layer.eps_r, tan_delta: layer.tan_delta_tls_0 }
}

/// Resolves `[interfaces.{ma,ms,sa}]` TLS surface-loss parameters for a named material from `ctx.library`.
pub fn tls_interface_from_material_name(name: &str, ctx: &Ctx) -> Option<TlsInterfaces> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return None;
    }
    let lib = ctx.library.as_ref()?;
    let mat = lib.get(trimmed)?;
    let out = TlsInterfaces {
        ma: mat.interfaces.ma.as_ref().map(map_interface_layer),
        ms: mat.interfaces.ms.as_ref().map(map_interface_layer),
        sa: mat.interfaces.sa.as_ref().map(map_interface_layer),
    };
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// Resolves `[interfaces.{ma,ms,sa}]` TLS surface-loss parameters from named/oriented domain materials in
/// `spec.materials` (and any BC carrying a `"material"` parameter, such as `surface_impedance` or `london_sheet`)
/// via `ctx.library`.
pub fn tls_interfaces_from_spec(spec: &ProblemSpec, ctx: &Ctx) -> TlsInterfaces {
    let mut out = TlsInterfaces::default();
    let mut merge = |ifaces: TlsInterfaces| {
        if out.ma.is_none() {
            out.ma = ifaces.ma;
        }
        if out.ms.is_none() {
            out.ms = ifaces.ms;
        }
        if out.sa.is_none() {
            out.sa = ifaces.sa;
        }
    };
    for assign in spec.materials.values() {
        let mat_name = match assign {
            MaterialAssign::Named(n) => Some(n.as_str()),
            MaterialAssign::Oriented { material, .. } => Some(material.as_str()),
            MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_) => None,
        };
        if let Some(ifaces) = mat_name.and_then(|n| tls_interface_from_material_name(n, ctx)) {
            merge(ifaces);
        }
    }
    for bc in &spec.bcs {
        if let Some(ifaces) =
            bc.params.get("material").and_then(Value::as_str).and_then(|n| tls_interface_from_material_name(n, ctx))
        {
            merge(ifaces);
        }
    }
    out
}

fn prop_value_to_m3(v: &PropValue) -> femcore::material::M3 {
    match v {
        PropValue::Scalar(s) => femcore::material::iso_tensor(*s),
        PropValue::Diag([x, y, z]) => [[*x, 0.0, 0.0], [0.0, *y, 0.0], [0.0, 0.0, *z]],
        PropValue::Sym6([xx, yy, zz, xy, xz, yz]) => [[*xx, *xy, *xz], [*xy, *yy, *yz], [*xz, *yz, *zz]],
        PropValue::Matrix(mat) => symmetrize_m3(mat),
        PropValue::Matrix3x6(_) => [[0.0; 3]; 3],
    }
}

fn mul_m3(a: &femcore::material::M3, b: &femcore::material::M3) -> femcore::material::M3 {
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

fn symmetrize_m3(a: &femcore::material::M3) -> femcore::material::M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = 0.5 * (a[i][j] + a[j][i]);
        }
    }
    out
}

fn set_inline_tensor(
    m: &mut Materials,
    region: &str,
    labels: &[i32],
    name: &str,
    tensor: femcore::material::Tensor,
) -> Result<(), Error> {
    for &l in labels {
        let an = m.aniso.entry(l).or_default();
        match name {
            "eps_r" => an.permittivity = Some(tensor.clone()),
            "eps_r_im" => an.permittivity_im = Some(tensor.clone()),
            "k" | "sigma" => an.conductivity = Some(tensor.clone()),
            "alpha" => an.expansion = Some(tensor.clone()),
            _ => {
                return Err(Error::invalid(format!(
                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                )));
            }
        }
    }
    Ok(())
}

fn apply_inline_tan_delta_tensors(
    m: &mut Materials,
    labels: &[i32],
    inline_td_tensor: Option<&femcore::material::Tensor>,
) {
    let Some(td_t) = inline_td_tensor else { return };
    let td_m = femcore::material::tensor_global(td_t, [0.0; 3]);
    for &l in labels {
        let eps_m = m
            .aniso
            .get(&l)
            .and_then(|an| an.permittivity.as_ref())
            .map(|t| femcore::material::tensor_global(t, [0.0; 3]))
            .unwrap_or_else(|| {
                let s = m.get(l, "eps_r").unwrap_or(1.0);
                [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]]
            });
        let mut im = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s1 = 0.0;
                let mut s2 = 0.0;
                for k in 0..3 {
                    s1 += eps_m[i][k] * td_m[k][j];
                    s2 += td_m[i][k] * eps_m[k][j];
                }
                im[i][j] = 0.5 * (s1 + s2);
            }
        }
        let an = m.aniso.entry(l).or_default();
        an.permittivity_im = Some(m3_to_sym_tensor(&im));
    }
}

fn m3_to_sym_tensor(m: &femcore::material::M3) -> femcore::material::Tensor {
    femcore::material::Tensor::Symmetric {
        xx: m[0][0],
        yy: m[1][1],
        zz: m[2][2],
        xy: m[0][1],
        yz: m[1][2],
        xz: m[0][2],
    }
}

fn is_isotropic_m3(m: &femcore::material::M3) -> bool {
    let scale = m[0][0].abs().max(m[1][1]).abs().max(m[2][2].abs()).max(1.0);
    (m[0][0] - m[1][1]).abs() + (m[1][1] - m[2][2]).abs() + m[0][1].abs() + m[0][2].abs() + m[1][2].abs()
        < 1e-12 * scale
}

#[allow(clippy::too_many_arguments)]
fn named(
    m: &mut Materials,
    region: &str,
    labels: &[i32],
    name: &str,
    euler: [f64; 3],
    ctx: &Ctx,
    freq_hz: Option<f64>,
    temp_k: Option<f64>,
    props: &[PropDesc],
    is_quasi_static_electro: bool,
    physics: &str,
) -> Result<(), Error> {
    use femcore::material::{Class, Tensor};
    let lib = ctx.library.as_ref().ok_or_else(|| {
        Error::Unsupported(format!(
            "region '{region}': the material '{name}' needs a material library, and this run has none; give the properties as numbers"
        ))
    })?;
    let mat = lib.require(name).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
    let wants = |p: &str| props.iter().any(|d| d.name == p);
    let mut legacy = mat.legacy_scalars();
    if temp_k.is_some() {
        if let Some(cp) = mat.specific_heat_at(temp_k) {
            legacy.insert("Cp".into(), cp);
        }
        if let Some(kt) = mat.thermal_conductivity_at(temp_k) {
            if let Some(s) = kt.as_scalar() {
                legacy.insert("k".into(), s);
            } else {
                legacy.remove("k");
            }
        }
        if let Some(at) = mat.thermal_expansion_at(temp_k) {
            if let Some(s) = at.as_scalar() {
                legacy.insert("alpha".into(), s);
            } else {
                legacy.remove("alpha");
            }
        }
    }
    for (k, v) in &legacy {
        if wants(k)
            || (wants("eps_r")
                && (matches!(k.as_str(), "tan_delta" | "eps_r_im") || (k == "sigma" && physics != "sms")))
        {
            for &l in labels {
                m.set(l, k, *v);
            }
        }
    }
    let mut an = mat.anisotropy(euler);
    if temp_k.is_some() {
        if let Some(kt) = mat.thermal_conductivity_at(temp_k) {
            an.conductivity = Some(m3_to_sym_tensor(&kt.matrix()));
        }
        if let Some(at) = mat.thermal_expansion_at(temp_k) {
            an.expansion = Some(m3_to_sym_tensor(&at.matrix()));
        }
    }
    // what the kernels read as scalars stays scalar: crystal data only where the material has a direction in it or is placed at an angle
    let tilted = euler.iter().any(|a| a.abs() > 1e-12);
    let elastic_aniso = mat.elastic.as_ref().is_some_and(|e| e.detect_class() != Class::Isotropic);
    let keep = |t: &Option<Tensor>, scalar: Option<f64>| -> Option<Tensor> {
        match (t, scalar) {
            (Some(_), Some(_)) if !tilted => None,
            (t, _) => t.clone(),
        }
    };
    an.conductivity = if wants("k") {
        keep(&an.conductivity, legacy.get("k").copied())
    } else if wants("sigma") || (wants("eps_r") && !wants("k") && physics != "sms") {
        let sig_t = mat.electrical.conductivity.as_ref().map(|tp| m3_to_sym_tensor(&tp.matrix()));
        keep(&sig_t, legacy.get("sigma").copied())
    } else {
        None
    };
    an.expansion = if wants("alpha") { keep(&an.expansion, legacy.get("alpha").copied()) } else { None };
    an.permittivity = if wants("eps_r") { keep(&an.permittivity, legacy.get("eps_r").copied()) } else { None };
    if !wants("E") {
        an.stiffness = None;
        an.viscosity = None;
        an.piezo_e = None;
    } else if !elastic_aniso && !tilted {
        an.stiffness = None;
    }
    if !wants("eps_r") {
        an.permittivity = None;
        an.permittivity_im = None;
    } else {
        // Optical dispersion applies when an optical/THz wavelength (<= 300 um, f >= 1 THz) is given and the material has [optical],
        // never for quasi-static dfs electrostatic or sms solid/piezoelectric studies.
        if let (false, Some(_), Some(lambda)) =
            (is_quasi_static_electro, &mat.optical, ctx.wavelength.filter(|&lam| lam <= 3e-4))
        {
            an.permittivity_im = None;
            let eps_cryst = mat
                .optical_permittivity(lambda, temp_k, [0.0; 3])
                .map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            let eps_glob = mat
                .optical_permittivity(lambda, temp_k, euler)
                .map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            let re_cryst = [
                [eps_cryst[0][0].re, eps_cryst[0][1].re, eps_cryst[0][2].re],
                [eps_cryst[1][0].re, eps_cryst[1][1].re, eps_cryst[1][2].re],
                [eps_cryst[2][0].re, eps_cryst[2][1].re, eps_cryst[2][2].re],
            ];
            let re_glob = [
                [eps_glob[0][0].re, eps_glob[0][1].re, eps_glob[0][2].re],
                [eps_glob[1][0].re, eps_glob[1][1].re, eps_glob[1][2].re],
                [eps_glob[2][0].re, eps_glob[2][1].re, eps_glob[2][2].re],
            ];
            if is_isotropic_m3(&re_glob) {
                for &l in labels {
                    m.set(l, "eps_r", re_glob[0][0]);
                }
                an.permittivity = None;
            } else {
                an.permittivity = Some(m3_to_sym_tensor(&re_cryst));
            }
            // In exp(+i omega t), passive optical loss (n - i k)^2 has Im(eps) <= 0, so loss magnitude is -Im(eps).
            let lossy = eps_cryst.iter().flatten().any(|z| z.im.abs() > 1e-12 * (1.0 + z.re.abs()));
            if lossy {
                let im_cryst = [
                    [-eps_cryst[0][0].im, -eps_cryst[0][1].im, -eps_cryst[0][2].im],
                    [-eps_cryst[1][0].im, -eps_cryst[1][1].im, -eps_cryst[1][2].im],
                    [-eps_cryst[2][0].im, -eps_cryst[2][1].im, -eps_cryst[2][2].im],
                ];
                let im_glob = [
                    [-eps_glob[0][0].im, -eps_glob[0][1].im, -eps_glob[0][2].im],
                    [-eps_glob[1][0].im, -eps_glob[1][1].im, -eps_glob[1][2].im],
                    [-eps_glob[2][0].im, -eps_glob[2][1].im, -eps_glob[2][2].im],
                ];
                if is_isotropic_m3(&im_glob) && an.permittivity.is_none() {
                    for &l in labels {
                        m.set(l, "eps_r_im", im_glob[0][0]);
                    }
                    an.permittivity_im = None;
                } else {
                    an.permittivity_im = Some(m3_to_sym_tensor(&im_cryst));
                }
            }
        } else if let Some((re_cryst, im_cryst)) =
            mat.complex_microwave_permittivity(freq_hz, temp_k, None, false, [0.0; 3])
        {
            let re_glob = femcore::material::tensor_global(&m3_to_sym_tensor(&re_cryst), euler);
            if is_isotropic_m3(&re_glob) {
                for &l in labels {
                    m.set(l, "eps_r", re_glob[0][0]);
                }
                an.permittivity = None;
            } else {
                an.permittivity = Some(m3_to_sym_tensor(&re_cryst));
            }
            let has_im = im_cryst.iter().flatten().any(|&v| v.abs() > 1e-15);
            if has_im {
                let im_glob = femcore::material::tensor_global(&m3_to_sym_tensor(&im_cryst), euler);
                if is_isotropic_m3(&im_glob) && an.permittivity.is_none() {
                    for &l in labels {
                        m.set(l, "eps_r_im", im_glob[0][0]);
                        if re_glob[0][0].abs() > 0.0 {
                            m.set(l, "tan_delta", im_glob[0][0] / re_glob[0][0]);
                        }
                    }
                    an.permittivity_im = None;
                } else {
                    an.permittivity_im = Some(m3_to_sym_tensor(&im_cryst));
                }
            } else {
                an.permittivity_im = None;
            }
        }
    }
    if an.stiffness.is_some()
        || an.viscosity.is_some()
        || an.piezo_e.is_some()
        || an.conductivity.is_some()
        || an.permittivity.is_some()
        || an.permittivity_im.is_some()
        || an.expansion.is_some()
    {
        for &l in labels {
            m.aniso.insert(l, an.clone());
        }
    }
    Ok(())
}

/// A condition of the problem with its selector resolved to numbers: boundary numbers when the descriptor says `Boundary`, domain
/// numbers when it says `Domain`.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBc {
    pub kind: String,
    pub labels: Vec<i32>,
    pub applies: Applies,
    pub params: Value,
}

/// Resolves the `on` selector of every condition. A selector that matches nothing is an error, naming the selector.
pub fn resolve_conditions(bcs: &[BcSpec], kinds: &[BcDesc], ctx: &Ctx) -> Result<Vec<ResolvedBc>, Error> {
    let mut out = vec![];
    for bc in bcs {
        let d = kinds
            .iter()
            .find(|k| k.kind == bc.kind)
            .ok_or_else(|| Error::Unsupported(format!("no condition '{}'", bc.kind)))?;
        let labels = match d.applies {
            Applies::Boundary => ctx.resolver.boundaries(&bc.on)?,
            Applies::Domain => ctx.resolver.domains(&bc.on)?,
        };
        if labels.is_empty() {
            return Err(Error::invalid(format!("condition '{}' on '{}' selects nothing", bc.kind, bc.on)));
        }
        out.push(ResolvedBc { kind: bc.kind.clone(), labels, applies: d.applies, params: bc.params.clone() });
    }
    Ok(out)
}

/// Names the regions a problem gives materials to, for the physics' `sizing`: region `i` of the result is label `i + 1` of the returned
/// `Materials`. Needs no mesh. A region key that is a pattern (`sapphire*`) stays one entry; the mesher resolves the pattern.
pub fn materials_by_region(
    spec: &ProblemSpec,
    library: Option<&std::sync::Arc<femcore::material::Library>>,
    wavelength: Option<f64>,
    props: &[PropDesc],
) -> Result<(Vec<String>, Materials), Error> {
    struct Index(Vec<String>);
    impl crate::physics::LabelResolver for Index {
        fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
            self.0.iter().position(|r| r == region).map(|i| vec![i as i32 + 1]).ok_or_else(|| Error::invalid(region))
        }
        fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
            Err(Error::invalid(selector))
        }
    }
    let names: Vec<String> = spec.materials.keys().cloned().collect();
    let resolver = Index(names.clone());
    let mut ctx = Ctx::new(&resolver);
    ctx.library = library.cloned();
    ctx.wavelength = wavelength;
    let m = materials_by_label(spec, &ctx, props)?;
    Ok((names, m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::LabelResolver;
    use crate::schema::params_schema;
    use crate::validate::check_params;
    use serde_json::json;
    use std::path::Path;
    use std::sync::Arc;

    struct DummyResolver;
    impl LabelResolver for DummyResolver {
        fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
            match region {
                "sub" => Ok(vec![1]),
                "vac" => Ok(vec![2]),
                other => Err(Error::invalid(format!("unknown region {other}"))),
            }
        }
        fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
            match selector {
                "metal" => Ok(vec![10]),
                other => Err(Error::invalid(format!("unknown boundary {other}"))),
            }
        }
    }

    #[test]
    fn param_kind_numbers_validates_decodes_and_generates_schema() {
        let descs = vec![ParamDesc::numbers("freqs", "Frequencies", "Hz", "explicit frequency list")];
        let schema = params_schema(&descs);
        assert_eq!(schema["properties"]["freqs"]["type"], json!("array"));
        assert_eq!(schema["properties"]["freqs"]["items"], json!({ "type": "number" }));

        let val_ok = json!({ "freqs": [5.0e9, 6.0e9, 7.5e9] });
        check_params(&val_ok, &descs, "driven").expect("valid numbers array");
        let p = Params::new(&val_ok, &descs);
        assert!(p.has("freqs"));
        assert_eq!(p.numbers("freqs").unwrap(), vec![5.0e9, 6.0e9, 7.5e9]);
        assert_eq!(p.optional_numbers("freqs").unwrap(), Some(vec![5.0e9, 6.0e9, 7.5e9]));

        let val_empty = json!({});
        check_params(&val_empty, &descs, "driven").expect("omitted numbers array uses default");
        let p_empty = Params::new(&val_empty, &descs);
        assert!(!p_empty.has("freqs"));
        assert!(p_empty.numbers("freqs").unwrap().is_empty());
        assert_eq!(p_empty.optional_numbers("freqs").unwrap(), None);

        let val_bad = json!({ "freqs": [5.0e9, "oops"] });
        assert!(check_params(&val_bad, &descs, "driven").is_err());
        assert!(matches!(Params::new(&val_bad, &descs).numbers("freqs"), Err(Error::Invalid(_))));
        assert!(matches!(Params::new(&val_ok, &descs).number("missing"), Err(Error::Invalid(_))));
        assert!(matches!(Params::new(&val_ok, &descs).text("freqs"), Err(Error::Invalid(_))));
        assert!(matches!(
            Params::new(&json!({ "freqs": "not-an-array" }), &descs).numbers("freqs"),
            Err(Error::Invalid(_))
        ));
    }

    #[test]
    fn inline_loss_tangent_populates_complex_permittivity_and_rejects_negative() {
        let props = vec![PropDesc::new("eps_r", "Permittivity", "", PropKind::Tensor, &[], Some(1.0))];
        let mut spec = ProblemSpec::new("ems", 2, "cavity");
        spec.materials.insert(
            "sub".into(),
            MaterialAssign::Props([("eps_r".into(), 11.45), ("loss_tangent".into(), 1.0e-4)].into()),
        );
        let ctx = Ctx::new(&DummyResolver);
        let m = materials_by_label(&spec, &ctx, &props).expect("valid inline loss_tangent");
        assert!(m.has_complex_permittivity());
        let eps_im = m.get(1, "eps_r_im").expect("eps_r_im populated");
        assert!((eps_im - 11.45 * 1.0e-4).abs() < 1e-15);

        // Anisotropic tensor eps_r + tensor loss_tangent
        let mut spec_t = ProblemSpec::new("ems", 2, "cavity");
        spec_t.materials.insert(
            "sub".into(),
            MaterialAssign::PropsTensor(
                [
                    ("eps_r".into(), PropValue::Diag([9.4, 9.4, 11.6])),
                    ("loss_tangent".into(), PropValue::Diag([1e-5, 1e-5, 2e-5])),
                ]
                .into(),
            ),
        );
        let mt = materials_by_label(&spec_t, &ctx, &props).expect("valid tensor loss_tangent");
        assert!(mt.has_complex_permittivity());
        assert!(mt.aniso.get(&1).and_then(|a| a.permittivity_im.as_ref()).is_some());

        // Negative loss_tangent is rejected cleanly
        let mut spec_neg = ProblemSpec::new("ems", 2, "cavity");
        spec_neg.materials.insert(
            "sub".into(),
            MaterialAssign::Props([("eps_r".into(), 11.45), ("loss_tangent".into(), -1.0e-6)].into()),
        );
        assert!(materials_by_label(&spec_neg, &ctx, &props).is_err());
    }

    #[test]
    fn tls_interfaces_resolve_from_material_library() {
        let lib_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
        let lib = Arc::new(femcore::material::Library::load_dir(&lib_dir).expect("load material library"));
        let mut ctx = Ctx::new(&DummyResolver);
        ctx.library = Some(lib);

        let si_ifaces = tls_interface_from_material_name("silicon_4k", &ctx).expect("silicon_4k has interfaces");
        assert!(si_ifaces.sa.is_some(), "silicon_4k [interfaces.sa]");
        assert!(si_ifaces.ms.is_some(), "silicon_4k [interfaces.ms]");
        assert_eq!(si_ifaces.ms.as_ref().unwrap().eps_r, 9.8);

        let al_ifaces =
            tls_interface_from_material_name("aluminium_sc_emft", &ctx).expect("aluminium_sc_emft has interfaces");
        assert!(al_ifaces.ma.is_some(), "aluminium_sc_emft [interfaces.ma]");

        let mut spec = ProblemSpec::new("ems", 2, "cavity");
        spec.materials.insert("sub".into(), MaterialAssign::Named("silicon_4k".into()));
        spec.bcs.push(BcSpec {
            on: "metal".into(),
            kind: "surface_impedance".into(),
            params: json!({ "material": "aluminium_sc_emft" }),
        });
        let combined = tls_interfaces_from_spec(&spec, &ctx);
        assert!(combined.ma.is_some());
        assert!(combined.ms.is_some());
        assert!(combined.sa.is_some());
    }

    #[test]
    fn sms_decode_keeps_piezo_on_isotropic_host_and_never_loads_optical_or_sigma() {
        let lib_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
        let lib = Arc::new(femcore::material::Library::load_dir(&lib_dir).expect("load material library"));
        let mut ctx = Ctx::new(&DummyResolver);
        ctx.library = Some(lib);
        ctx.wavelength = Some(1.55e-6); // Must be ignored for sms!

        let sms_props = vec![
            PropDesc::new("E", "Young's modulus", "Pa", PropKind::Scalar, &[], None),
            PropDesc::new("nu", "Poisson's ratio", "", PropKind::Scalar, &[], None),
            PropDesc::new("rho", "Density", "kg/m^3", PropKind::Scalar, &[], None),
            PropDesc::new("eps_r", "Relative permittivity", "", PropKind::Tensor, &[], Some(1.0)),
        ];
        let mut spec = ProblemSpec::new("sms", 1, "eigenfrequency");
        spec.materials.insert("sub".into(), MaterialAssign::Named("aluminum_nitride_piezomumps".into()));
        spec.materials.insert("vac".into(), MaterialAssign::Named("gold".into()));
        let m = materials_by_label(&spec, &ctx, &sms_props).expect("decode sms materials");
        assert!(m.has_piezoelectric(), "piezo_e must survive even when E/nu and eps_r are isotropic");
        // Optical n^2 at 1.55 um is ~4.07, whereas dielectric eps_clamped is 9.5
        assert!((m.get(1, "eps_r").unwrap() - 9.5).abs() < 1e-12);
        // Electrical conductivity sigma from gold must not leak into sms
        assert!(m.get(2, "sigma").is_none());
        assert!(m.aniso.get(&2).and_then(|a| a.conductivity.as_ref()).is_none());

        // Non-mechanical physics (!wants("E")) strips piezo_e
        let dfs_props = vec![PropDesc::new("eps_r", "Permittivity", "", PropKind::Tensor, &[], Some(1.0))];
        let mut dfs_spec = ProblemSpec::new("dfs", 1, "electrostatic");
        dfs_spec.materials.insert("sub".into(), MaterialAssign::Named("aluminum_nitride_piezomumps".into()));
        let m_dfs = materials_by_label(&dfs_spec, &ctx, &dfs_props).expect("decode dfs materials");
        assert!(!m_dfs.has_piezoelectric());
    }
}
