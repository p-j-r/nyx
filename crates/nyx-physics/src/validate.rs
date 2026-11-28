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
    let has_eps_r = props.iter().any(|d| d.name == "eps_r");
    let is_gyro_key = |name: &str| {
        matches!(
            name,
            "gyrotropy"
                | "gyrotropy_im"
                | "gyro_e"
                | "gyro_e_im"
                | "mu_gyrotropy"
                | "mu_gyrotropy_im"
                | "gyro_m"
                | "gyro_m_im"
        )
    };
    let is_eo_pe_key = |name: &str| {
        p.id() == "ems"
            && matches!(
                name,
                "r11"
                    | "r_11"
                    | "r13"
                    | "r_13"
                    | "r22"
                    | "r_22"
                    | "r33"
                    | "r_33"
                    | "r41"
                    | "r_41"
                    | "r42"
                    | "r_42"
                    | "r51"
                    | "r_51"
                    | "p11"
                    | "p_11"
                    | "p12"
                    | "p_12"
                    | "p13"
                    | "p_13"
                    | "p31"
                    | "p_31"
                    | "p33"
                    | "p_33"
                    | "p14"
                    | "p_14"
                    | "p41"
                    | "p_41"
                    | "p44"
                    | "p_44"
            )
    };
    let is_aux_em_tensor = |name: &str| {
        has_eps_r
            && (matches!(name, "loss_tangent" | "tan_delta" | "sigma" | "eps_r_im" | "mu_r_im")
                || is_gyro_key(name)
                || is_eo_pe_key(name)
                || (p.id() == "ems" && name == "mu_r"))
    };
    for (region, m) in &spec.materials {
        match m {
            MaterialAssign::Props(map) => {
                for (name, v) in map {
                    let desc = props.iter().find(|d| &d.name == name);
                    let allowed = desc.is_some() || is_aux_em_tensor(name);
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
                    let desc_kind = if is_aux_em_tensor(name) {
                        PropKind::Tensor
                    } else if let Some(desc) = props.iter().find(|d| &d.name == name) {
                        desc.kind
                    } else {
                        return Err(Error::invalid(format!(
                            "region '{region}': {} does not read a property '{name}' (it reads: {})",
                            p.name(),
                            props.iter().map(|d| d.name.as_str()).collect::<Vec<_>>().join(", ")
                        )));
                    };
                    let allow_non_symmetric = (p.id() == "ems"
                        && matches!(name.as_str(), "eps_r" | "eps_r_im" | "mu_r" | "mu_r_im"))
                        || is_gyro_key(name);
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
                            if !allow_non_symmetric
                                && ((mat[0][1] - mat[1][0]).abs() > tol
                                    || (mat[0][2] - mat[2][0]).abs() > tol
                                    || (mat[1][2] - mat[2][1]).abs() > tol)
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

/// Every field and scalar solution metric of the result must be finite, with one physical
/// exception: `SolveResult::scalars` allows `+inf` (for example lossless cavity `Q = inf`), while
/// still rejecting `NaN` in scalars and any non-finite value in fields, cell fields, eigenvalues,
/// or solution scalar values.
pub fn check_result_finite(result: &crate::result::SolveResult) -> Result<(), Error> {
    for (k, sol) in result.solutions.iter().enumerate() {
        if let Some(v) = sol.value {
            if !v.is_finite() {
                return Err(Error::NonFinite(format!(
                    "the solution is not finite: value of solution {} ('{}') is {v}; check materials, boundary conditions and the mesh",
                    k + 1,
                    sol.label,
                )));
            }
        }
        if let Some(ev) = sol.eigenvalue {
            if !ev[0].is_finite() || !ev[1].is_finite() {
                return Err(Error::NonFinite(format!(
                    "the solution is not finite: eigenvalue of solution {} ('{}') is {ev:?}; check materials, boundary conditions and the mesh",
                    k + 1,
                    sol.label,
                )));
            }
        }
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
        for cf in &sol.cell_fields {
            let bad_re = cf.re.iter().filter(|v| !v.is_finite()).count();
            let bad_im = cf.im.iter().filter(|v| !v.is_finite()).count();
            if bad_re > 0 || bad_im > 0 {
                return Err(Error::NonFinite(format!(
                    "the solution is not finite: cell field '{}' of solution {} ('{}') has {bad_re}/{} bad real and {bad_im}/{} bad imaginary values; check materials, boundary conditions and the mesh",
                    cf.name,
                    k + 1,
                    sol.label,
                    cf.re.len(),
                    cf.im.len(),
                )));
            }
        }
    }
    for (name, scalar) in &result.scalars {
        if let crate::result::Scalar::Num(v) = scalar {
            if v.is_nan() {
                return Err(Error::NonFinite(format!(
                    "the solution is not finite: scalar '{name}' is NaN; check materials, boundary conditions and the mesh"
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result::{CellField, Field, Scalar, Solution, SolveResult};

    fn base_result() -> SolveResult {
        SolveResult {
            physics: "ems".into(),
            study: "cavity".into(),
            order: 1,
            dim: 3,
            dofs: 1,
            seconds: 0.1,
            solutions: vec![Solution {
                label: "mode 1".into(),
                fields: vec![Field::scalar("V", "V", vec![1.0, 2.0])],
                value: Some(1.5e9),
                eigenvalue: Some([2.4, -0.01]),
                cell_fields: vec![CellField {
                    name: "E".into(),
                    unit: "V/m".into(),
                    nloc: 1,
                    re: vec![0.1, 0.2, 0.3],
                    im: vec![0.0, -0.1, 0.0],
                }],
            }],
            scalars: vec![("f1".into(), Scalar::Num(1.5e9)), ("note".into(), Scalar::Text("ok".into()))],
            log: vec![],
            diagnostics: None,
            error_budget: None,
        }
    }

    #[test]
    fn check_result_finite_allows_lossless_q_infinity() {
        let mut r = base_result();
        r.scalars.push(("Q1".into(), Scalar::Num(f64::INFINITY)));
        check_result_finite(&r).unwrap();
    }

    #[test]
    fn check_result_finite_rejects_every_non_finite_slot() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut r = base_result();
            r.solutions[0].fields[0].values[1] = bad;
            assert!(check_result_finite(&r).is_err());

            let mut r = base_result();
            r.solutions[0].cell_fields[0].re[0] = bad;
            assert!(check_result_finite(&r).is_err());

            let mut r = base_result();
            r.solutions[0].cell_fields[0].im[2] = bad;
            assert!(check_result_finite(&r).is_err());

            let mut r = base_result();
            r.solutions[0].eigenvalue = Some([bad, 0.0]);
            assert!(check_result_finite(&r).is_err());

            let mut r = base_result();
            r.solutions[0].eigenvalue = Some([1.0, bad]);
            assert!(check_result_finite(&r).is_err());

            let mut r = base_result();
            r.solutions[0].value = Some(bad);
            assert!(check_result_finite(&r).is_err());
        }

        let mut r = base_result();
        r.scalars.push(("Q1".into(), Scalar::Num(f64::NAN)));
        assert!(check_result_finite(&r).is_err());
    }
}
