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

/// The material properties of every region as numbers per domain label: a region name or pattern (`sapphire*`) can match several
/// domains. A property outside `props` is refused for inline properties (`validate_spec` already does that; this keeps a kernel safe when
/// called directly). A named material comes from the library of `ctx` (`MaterialAssign::Named`, or `Oriented` with Bunge angles): the
/// scalar properties the physics lists are taken under their descriptor names, the crystal data (stiffness, tensors) as an `Anisotropy`
/// where the material is anisotropic or oriented, and `eps_r` from the optical model at `ctx.wavelength` where the material has one.
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
            MaterialAssign::Named(n) => named(&mut m, region, &labels, n, [0.0; 3], ctx, props)?,
            MaterialAssign::Oriented { material, euler_deg } => {
                named(&mut m, region, &labels, material, *euler_deg, ctx, props)?
            }
        }
    }
    Ok(m)
}

fn named(
    m: &mut Materials,
    region: &str,
    labels: &[i32],
    name: &str,
    euler: [f64; 3],
    ctx: &Ctx,
    props: &[PropDesc],
) -> Result<(), Error> {
    use femcore::material::{Class, Tensor};
    let lib = ctx.library.as_ref().ok_or_else(|| {
        Error::Unsupported(format!(
            "region '{region}': the material '{name}' needs a material library, and this run has none; give the properties as numbers"
        ))
    })?;
    let mat = lib.require(name).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
    let wants = |p: &str| props.iter().any(|d| d.name == p);
    for (k, v) in mat.legacy_scalars() {
        if wants(&k) {
            for &l in labels {
                m.set(l, &k, v);
            }
        }
    }
    let mut an = mat.anisotropy(euler);
    // what the kernels read as scalars stays scalar: crystal data only where the material has a direction in it or is placed at an angle
    let tilted = euler.iter().any(|a| a.abs() > 1e-12);
    let elastic_aniso = mat.elastic.as_ref().is_some_and(|e| e.detect_class() != Class::Isotropic);
    let keep = |t: &Option<Tensor>, scalar: Option<f64>| -> Option<Tensor> {
        match (t, scalar) {
            (Some(_), Some(_)) if !tilted => None,
            (t, _) => t.clone(),
        }
    };
    let legacy = mat.legacy_scalars();
    an.conductivity = keep(&an.conductivity, legacy.get("k").copied());
    an.expansion = keep(&an.expansion, legacy.get("alpha").copied());
    an.permittivity = keep(&an.permittivity, legacy.get("eps_r").copied());
    if !elastic_aniso && !tilted {
        an.stiffness = None;
    }
    if wants("eps_r") && mat.optical.is_some() {
        if let Some(lambda) = ctx.wavelength {
            let eps = mat.optical_permittivity(lambda, None, euler).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            let lossy = eps.iter().flatten().any(|z| z.im.abs() > 1e-9 * (1.0 + z.re.abs()));
            if lossy {
                return Err(Error::Unsupported(format!(
                    "region '{region}': the material '{name}' is lossy at {:.0} nm (permittivity with an imaginary part); the mode solver takes real permittivities, use an absorbing layer for loss",
                    lambda * 1e9
                )));
            }
            let g = |i: usize, j: usize| eps[i][j].re;
            let iso = (g(0, 0) - g(1, 1)).abs() + (g(1, 1) - g(2, 2)).abs() + g(0, 1).abs() + g(0, 2).abs() + g(1, 2).abs() < 1e-12;
            if iso {
                for &l in labels {
                    m.set(l, "eps_r", g(0, 0));
                }
                an.permittivity = None;
            } else {
                // the optical tensor is already in the wafer frame
                an.permittivity = Some(Tensor::Symmetric { xx: g(0, 0), yy: g(1, 1), zz: g(2, 2), xy: g(0, 1), yz: g(1, 2), xz: g(0, 2) });
                an.euler_deg = [0.0; 3];
            }
        }
    }
    if an.stiffness.is_some() || an.conductivity.is_some() || an.permittivity.is_some() || an.expansion.is_some() {
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
