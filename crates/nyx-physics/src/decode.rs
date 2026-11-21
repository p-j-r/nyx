// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Decoding the problem envelope into what a kernel reads: typed parameters with the descriptor defaults, region names resolved to
//! domain numbers (`materials_by_label`), condition selectors resolved to boundary or domain numbers (`resolve_conditions`). Every physics
//! crate uses these, so that the envelope is read the same way in all of them (ADR 41).

use crate::desc::{Applies, BcDesc, ParamDesc, ParamKind, PropDesc, PropKind};
use crate::error::Error;
use crate::physics::Ctx;
use crate::spec::{BcSpec, MaterialAssign, ProblemSpec, PropValue};
use femcore::props::Materials;
use serde_json::Value;

/// Typed access to the parameters of a study, a condition or the physics: a value that was given, else the default of the descriptor.
/// `Params::new` is called after `validate_spec` has checked types and ranges, so a getter on a name the descriptors do not list is a
/// bug of the physics crate and panics with the name.
pub struct Params<'a> {
    value: &'a Value,
    descs: &'a [ParamDesc],
}

impl<'a> Params<'a> {
    pub fn new(value: &'a Value, descs: &'a [ParamDesc]) -> Params<'a> {
        Params { value, descs }
    }

    fn desc(&self, name: &str) -> &ParamDesc {
        self.descs
            .iter()
            .find(|d| d.name == name)
            .unwrap_or_else(|| panic!("parameter '{name}' is not in the descriptors"))
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
    pub fn optional_number(&self, name: &str) -> Option<f64> {
        if self.is_null(name) {
            None
        } else {
            Some(self.number(name))
        }
    }

    /// Whether the parameter was given with a non-null value in the problem (not defaulted).
    pub fn has(&self, name: &str) -> bool {
        self.given(name).is_some()
    }

    pub fn number(&self, name: &str) -> f64 {
        match (self.given(name).and_then(Value::as_f64), &self.desc(name).kind) {
            (Some(x), _) => x,
            (None, ParamKind::Number { default, .. }) => *default,
            _ => panic!("parameter '{name}' is not a number"),
        }
    }

    pub fn integer(&self, name: &str) -> i64 {
        match (self.given(name).and_then(Value::as_i64), &self.desc(name).kind) {
            (Some(x), _) => x,
            (None, ParamKind::Integer { default, .. }) => *default,
            _ => panic!("parameter '{name}' is not an integer"),
        }
    }

    pub fn boolean(&self, name: &str) -> bool {
        match (self.given(name).and_then(Value::as_bool), &self.desc(name).kind) {
            (Some(x), _) => x,
            (None, ParamKind::Bool { default }) => *default,
            _ => panic!("parameter '{name}' is not a boolean"),
        }
    }

    /// Text or a choice.
    pub fn text(&self, name: &str) -> String {
        match (self.given(name).and_then(Value::as_str), &self.desc(name).kind) {
            (Some(x), _) => x.to_string(),
            (None, ParamKind::Text { default }) | (None, ParamKind::Choice { default, .. }) => default.clone(),
            _ => panic!("parameter '{name}' is not text"),
        }
    }

    pub fn vector(&self, name: &str) -> Vec<f64> {
        match (self.given(name).and_then(Value::as_array), &self.desc(name).kind) {
            (Some(a), _) => a.iter().filter_map(Value::as_f64).collect(),
            (None, ParamKind::Vector { default, .. }) => default.clone(),
            _ => panic!("parameter '{name}' is not a vector"),
        }
    }

    pub fn names(&self, name: &str) -> Vec<String> {
        match (self.given(name).and_then(Value::as_array), &self.desc(name).kind) {
            (Some(a), _) => a.iter().filter_map(|v| v.as_str().map(str::to_string)).collect(),
            (None, ParamKind::Names { default }) => default.clone(),
            _ => panic!("parameter '{name}' is not a list of names"),
        }
    }
}

/// The material properties of every region as numbers per domain label: a region name or pattern (`sapphire*`) can match several
/// domains. A property outside `props` is refused for inline properties (`validate_spec` already does that; this keeps a kernel safe when
/// called directly). A named material comes from the library of `ctx` (`MaterialAssign::Named`, or `Oriented` with Bunge angles): the
/// scalar properties the physics lists are taken under their descriptor names, the crystal data (stiffness, tensors) as an `Anisotropy`
/// where the material is anisotropic or oriented, and `eps_r` from the optical model at `ctx.wavelength` where the material has one.
pub fn materials_by_label(spec: &ProblemSpec, ctx: &Ctx, props: &[PropDesc]) -> Result<Materials, Error> {
    let mut m = Materials::default();
    let is_quasi_static_electro = spec.physics == "dfs" && spec.study.kind == "electrostatic";
    let has_eps_r = props.iter().any(|d| d.name == "eps_r");
    let find_prop_kind =
        |name: &str| -> Option<PropKind> {
            props.iter().find(|d| d.name == name).map(|d| d.kind).or_else(|| {
                (has_eps_r && matches!(name, "tan_delta" | "sigma" | "eps_r_im")).then_some(PropKind::Tensor)
            })
        };
    let spec_freq = spec
        .study
        .params
        .get("freq")
        .or_else(|| spec.study.params.get("frequency"))
        .and_then(Value::as_f64)
        .or_else(|| spec.study.params.get("f_search").and_then(Value::as_f64))
        .filter(|&f| f > 0.0 && f.is_finite());
    let freq_hz = ctx
        .frequency
        .or_else(|| {
            (!is_quasi_static_electro)
                .then_some(ctx.wavelength)
                .flatten()
                .filter(|&l| l > 0.0)
                .map(|l| 299_792_458.0 / l)
        })
        .or(spec_freq);
    let temp_k = ctx
        .temperature
        .or_else(|| spec.params.get("temperature").and_then(Value::as_f64).filter(|&t| t > 0.0 && t.is_finite()));
    for (region, assign) in &spec.materials {
        let labels = ctx.resolver.domains(region)?;
        match assign {
            MaterialAssign::Props(map) => {
                for (name, v) in map {
                    if find_prop_kind(name).is_none() {
                        return Err(Error::invalid(format!("region '{region}': unknown property '{name}'")));
                    }
                    for &l in &labels {
                        m.set(l, name, *v);
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
                            if !s.is_finite() {
                                return Err(Error::invalid(format!(
                                    "region '{region}': property '{name}' is not a finite number"
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
                            let diag_m3 = [[*x, 0.0, 0.0], [0.0, *y, 0.0], [0.0, 0.0, *z]];
                            if is_isotropic_m3(&diag_m3) {
                                for &l in &labels {
                                    m.set(l, name, *x);
                                }
                            } else {
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
                            if is_isotropic_m3(&sym_m3) {
                                for &l in &labels {
                                    m.set(l, name, *xx);
                                }
                            } else {
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
                            if is_isotropic_m3(mat) {
                                for &l in &labels {
                                    m.set(l, name, mat[0][0]);
                                }
                            } else {
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
                    }
                }
                apply_inline_tan_delta_tensors(&mut m, &labels, inline_td_tensor.as_ref());
            }
            MaterialAssign::Named(n) => {
                named(&mut m, region, &labels, n, [0.0; 3], ctx, freq_hz, temp_k, props, is_quasi_static_electro)?
            }
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
            )?,
        }
    }
    Ok(m)
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
        if wants(k) || (wants("eps_r") && matches!(k.as_str(), "tan_delta" | "sigma" | "eps_r_im")) {
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
    } else if wants("sigma") || (wants("eps_r") && !wants("k")) {
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
    } else if !elastic_aniso && !tilted {
        an.stiffness = None;
    }
    if !wants("eps_r") {
        an.permittivity = None;
        an.permittivity_im = None;
    } else {
        // Optical dispersion applies when an optical/THz wavelength (<= 300 um, f >= 1 THz) is given and the material has [optical],
        // never for quasi-static dfs electrostatic studies.
        let use_optical =
            !is_quasi_static_electro && mat.optical.is_some() && ctx.wavelength.is_some_and(|lam| lam <= 3e-4);
        if use_optical {
            let lambda = ctx.wavelength.unwrap();
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
