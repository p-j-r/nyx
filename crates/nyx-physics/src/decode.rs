//! Decoding the problem envelope into what a kernel reads: typed parameters with the descriptor defaults, region names resolved to
//! domain numbers (`materials_by_label`), condition selectors resolved to boundary or domain numbers (`resolve_conditions`). Every physics
//! crate uses these, so that the envelope is read the same way in all of them (ADR 41).

use crate::desc::{Applies, BcDesc, ParamDesc, ParamKind, PropDesc};
use crate::error::Error;
use crate::physics::Ctx;
use crate::spec::{BcSpec, MaterialAssign, ProblemSpec};
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
        self.value.get(name)
    }

    /// Whether the parameter was given in the problem (not defaulted).
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

/// The scalar material properties of every region as numbers per domain label: a region name or pattern (`sapphire*`) can match several
/// domains. A property outside `props` is refused (`validate_spec` already does that; this keeps a kernel safe when called directly).
/// A named material is not available before the material files exist (phase 6).
pub fn materials_by_label(spec: &ProblemSpec, ctx: &Ctx, props: &[PropDesc]) -> Result<Materials, Error> {
    let mut m = Materials::default();
    for (region, assign) in &spec.materials {
        let labels = ctx.resolver.domains(region)?;
        match assign {
            MaterialAssign::Props(map) => {
                for (name, v) in map {
                    if !props.iter().any(|d| &d.name == name) {
                        return Err(Error::invalid(format!("region '{region}': unknown property '{name}'")));
                    }
                    for &l in &labels {
                        m.set(l, name, *v);
                    }
                }
            }
            MaterialAssign::Named(n) => {
                return Err(Error::Unsupported(format!(
                    "region '{region}': the material '{n}' cannot be read yet, material files arrive with restructure phase 6; give the properties as numbers"
                )))
            }
        }
    }
    Ok(m)
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
