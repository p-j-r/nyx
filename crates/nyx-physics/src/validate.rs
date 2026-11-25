// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Structural validation of a problem against the descriptors of its physics: the study and condition kinds exist, parameters have the
//! declared type and range, unknown keys are reported, material property names are ones the physics reads (ROUGH_SPOTS R12). What
//! needs the mesh or the physics itself stays in `Physics::validate`.

use crate::desc::{ParamDesc, ParamKind, PropKind};
use crate::error::Error;
use crate::physics::Physics;
use crate::spec::{MaterialAssign, ProblemSpec, PropValue};
use serde_json::Value;

pub fn validate_spec(p: &dyn Physics, spec: &ProblemSpec) -> Result<(), Error> {
    if spec.format != crate::spec::FORMAT {
        return Err(Error::invalid(format!(
            "problem format {} is not supported (this program reads format {})",
            spec.format,
            crate::spec::FORMAT
        )));
    }
    if spec.physics != p.id() {
        return Err(Error::invalid(format!("the problem is for '{}' but was given to '{}'", spec.physics, p.id())));
    }
    // the orders a physics supports are its own business (`Physics::validate`); the envelope only refuses nonsense
    if spec.order == 0 || spec.order > 9 {
        return Err(Error::invalid(format!("element order must be a whole number from 1 to 9, got {}", spec.order)));
    }
    let studies = p.studies();
    let study = studies.iter().find(|s| s.kind == spec.study.kind).ok_or_else(|| {
        Error::Unsupported(format!(
            "{} has no study '{}' (it has: {})",
            p.name(),
            spec.study.kind,
            studies.iter().map(|s| s.kind.as_str()).collect::<Vec<_>>().join(", ")
        ))
    })?;
    check_params(&spec.study.params, &study.params, &format!("study '{}'", study.kind))?;
    check_params(&spec.params, &p.params(), "the physics options")?;
    let kinds = p.bc_kinds();
    for (i, bc) in spec.bcs.iter().enumerate() {
        let d = kinds.iter().find(|k| k.kind == bc.kind).ok_or_else(|| {
            Error::Unsupported(format!(
                "{} has no condition '{}' (it has: {})",
                p.name(),
                bc.kind,
                kinds.iter().map(|k| k.kind.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        if bc.on.trim().is_empty() {
            return Err(Error::invalid(format!(
                "condition {} ('{}') does not say where it applies (`on`)",
                i + 1,
                bc.kind
            )));
        }
        check_params(&bc.params, &d.params, &format!("condition '{}'", bc.kind))?;
    }
    let props = p.properties();
    for (region, m) in &spec.materials {
        match m {
            MaterialAssign::Props(map) => {
                for (name, v) in map {
                    let desc = props.iter().find(|d| &d.name == name);
                    let allowed = desc.is_some() || (name == "loss_tangent" && props.iter().any(|d| d.name == "eps_r"));
                    if !allowed {
                        return Err(Error::invalid(format!(
                            "region '{region}': {} does not read a property '{name}' (it reads: {})",
                            p.name(),
                            props.iter().map(|d| d.name.as_str()).collect::<Vec<_>>().join(", ")
                        )));
                    }
                    if desc.is_some_and(|d| d.kind == PropKind::Matrix3x6) {
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
                }
            }
            MaterialAssign::PropsTensor(map) => {
                for (name, v) in map {
                    let desc_kind = if let Some(desc) = props.iter().find(|d| &d.name == name) {
                        if name == "loss_tangent" {
                            PropKind::Tensor
                        } else {
                            desc.kind
                        }
                    } else if name == "loss_tangent" && props.iter().any(|d| d.name == "eps_r") {
                        PropKind::Tensor
                    } else {
                        return Err(Error::invalid(format!(
                            "region '{region}': {} does not read a property '{name}' (it reads: {})",
                            p.name(),
                            props.iter().map(|d| d.name.as_str()).collect::<Vec<_>>().join(", ")
                        )));
                    };
                    match v {
                        PropValue::Scalar(x) => {
                            if desc_kind == PropKind::Matrix3x6 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a 3x6 matrix, not a scalar number"
                                )));
                            }
                            if !x.is_finite() {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' is not a finite number"
                                )));
                            }
                            if name == "loss_tangent" && *x < 0.0 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property 'loss_tangent' must be non-negative, got {x}"
                                )));
                            }
                        }
                        PropValue::Diag(d) => {
                            if desc_kind != PropKind::Tensor {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                                )));
                            }
                            if d.iter().any(|x| !x.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite tensor component"
                                )));
                            }
                            if name == "loss_tangent" && d.iter().any(|&x| x < 0.0) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property 'loss_tangent' must be non-negative in every direction"
                                )));
                            }
                        }
                        PropValue::Sym6(s6) => {
                            if desc_kind != PropKind::Tensor {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' must be a scalar number, not a tensor"
                                )));
                            }
                            if s6.iter().any(|x| !x.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite tensor component"
                                )));
                            }
                            if name == "loss_tangent" {
                                femcore::material::tensor::TensorProp::Symmetric([
                                    s6[0], s6[1], s6[2], s6[5], s6[4], s6[3],
                                ])
                                .check(femcore::material::tensor::Requirement::NonNegative, "loss_tangent")
                                .map_err(|msg| Error::invalid(format!("region '{region}': {msg}")))?;
                            }
                        }
                        PropValue::Matrix(mat) => {
                            if desc_kind != PropKind::Tensor {
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
                        }
                        PropValue::Matrix3x6(mat) => {
                            if desc_kind != PropKind::Matrix3x6 {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' cannot use a 3x6 matrix"
                                )));
                            }
                            if mat.iter().flatten().any(|x| !x.is_finite()) {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' has a non-finite 3x6 matrix component"
                                )));
                            }
                        }
                    }
                }
            }
            MaterialAssign::Named(_) | MaterialAssign::Oriented { .. } => {}
        }
    }
    Ok(())
}

/// Checks `value` (null or an object) against `descs`.
pub fn check_params(value: &Value, descs: &[ParamDesc], what: &str) -> Result<(), Error> {
    let empty = serde_json::Map::new();
    let obj = match value {
        Value::Null => &empty,
        Value::Object(m) => m,
        _ => return Err(Error::invalid(format!("{what}: parameters must be an object"))),
    };
    for k in obj.keys() {
        if !descs.iter().any(|d| &d.name == k) {
            let known: Vec<&str> = descs.iter().map(|d| d.name.as_str()).collect();
            return Err(Error::invalid(format!(
                "{what}: unknown parameter '{k}' (known: {})",
                if known.is_empty() { "none".into() } else { known.join(", ") }
            )));
        }
    }
    for d in descs {
        match obj.get(&d.name) {
            None | Some(Value::Null) => {
                if d.required {
                    return Err(Error::invalid(format!("{what}: parameter '{}' ({}) is required", d.name, d.label)));
                }
            }
            Some(v) => check_value(v, d).map_err(|m| Error::invalid(format!("{what}: parameter '{}': {m}", d.name)))?,
        }
    }
    Ok(())
}

fn check_value(v: &Value, d: &ParamDesc) -> Result<(), String> {
    match &d.kind {
        ParamKind::Number { min, max, .. } => {
            let x = v.as_f64().filter(|x| x.is_finite()).ok_or("must be a finite number")?;
            if min.is_some_and(|m| x < m) || max.is_some_and(|m| x > m) {
                return Err(format!(
                    "{x} is outside the range {}..{}",
                    min.map_or("-inf".into(), |m| m.to_string()),
                    max.map_or("inf".into(), |m| m.to_string())
                ));
            }
        }
        ParamKind::Integer { min, max, .. } => {
            let x = v.as_i64().ok_or("must be an integer")?;
            if min.is_some_and(|m| x < m) || max.is_some_and(|m| x > m) {
                return Err(format!(
                    "{x} is outside the range {}..{}",
                    min.map_or("-inf".into(), |m| m.to_string()),
                    max.map_or("inf".into(), |m| m.to_string())
                ));
            }
        }
        ParamKind::Bool { .. } => {
            v.as_bool().ok_or("must be true or false")?;
        }
        ParamKind::Choice { options, .. } => {
            let s = v.as_str().ok_or("must be text")?;
            if !options.iter().any(|o| o == s) {
                return Err(format!("'{s}' is not one of {}", options.join(", ")));
            }
        }
        ParamKind::Vector { len, .. } => {
            let a = v.as_array().ok_or("must be a list of numbers")?;
            if a.len() != *len || a.iter().any(|x| !x.as_f64().is_some_and(f64::is_finite)) {
                return Err(format!("must be {len} finite numbers"));
            }
        }
        ParamKind::Text { .. } => {
            v.as_str().ok_or("must be text")?;
        }
        ParamKind::Names { .. } => {
            let a = v.as_array().ok_or("must be a list of region names")?;
            if a.iter().any(|x| !x.is_string()) {
                return Err("must be a list of region names".into());
            }
        }
        ParamKind::Numbers { .. } => {
            let a = v.as_array().ok_or("must be a list of numbers")?;
            if a.iter().any(|x| !x.as_f64().is_some_and(f64::is_finite)) {
                return Err("must be a list of finite numbers".into());
            }
        }
    }
    Ok(())
}

/// Every field of the result must be finite: NaN or infinity does not stop a factorisation or a Krylov iteration, it flows into every
/// field. The message names the solution, the field and how many values are bad. Scalar results (for example Q = infinity of a lossless
/// mode) may legitimately be infinite and are not checked.
pub fn check_result_finite(result: &crate::result::SolveResult) -> Result<(), Error> {
    for (k, sol) in result.solutions.iter().enumerate() {
        for f in &sol.fields {
            let bad = f.values.iter().filter(|v| !v.is_finite()).count();
            if bad > 0 {
                return Err(Error::NonFinite(format!(
                    "the solution is not finite: field '{}' of solution {} ('{}') has {bad} NaN/infinite values out of {}; check materials, boundary conditions and the mesh",
                    f.name,
                    k + 1,
                    sol.label,
                    f.values.len()
                )));
            }
        }
    }
    Ok(())
}
