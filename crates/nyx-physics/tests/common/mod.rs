//! Shared helpers: a mock physics, a table resolver, and the exemplar file check used by the contract tests.
#![allow(dead_code)]

use nyx_physics::*;
use serde_json::Value;
use std::path::PathBuf;

pub struct Mock {
    pub id: &'static str,
    pub exports: Vec<FieldLink>,
    pub imports: Vec<FieldLink>,
}

impl Mock {
    pub fn new(id: &'static str) -> Mock {
        Mock { id, exports: vec![], imports: vec![] }
    }
}

impl Physics for Mock {
    fn id(&self) -> &'static str {
        self.id
    }
    fn name(&self) -> &'static str {
        "Mock physics"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc {
            kind: "modes".into(),
            name: "Modes".into(),
            doc: "a study".into(),
            dims: Dims::Two,
            params: vec![
                ParamDesc::number("freq", "Frequency", "Hz", 1e9, "").required().range(Some(1.0), Some(1e18)),
                ParamDesc::integer("n_modes", "Number of modes", 4, 1, ""),
                ParamDesc::boolean("evanescent", "Evanescent", false, ""),
                ParamDesc::choice("kind", "Kind", &["a", "b"], "a", ""),
                ParamDesc::vector("v", "Vector", "m", 3, ""),
                ParamDesc::text("expr", "Expression", "x", ""),
                ParamDesc::names("skip", "Skipped regions", ""),
            ],
        }]
    }
    fn bc_kinds(&self) -> Vec<BcDesc> {
        vec![BcDesc {
            kind: "fixed".into(),
            name: "Fixed".into(),
            doc: "".into(),
            applies: Applies::Boundary,
            params: vec![ParamDesc::number("ux", "Displacement", "m", 0.0, "")],
        }]
    }
    fn properties(&self) -> Vec<PropDesc> {
        vec![
            PropDesc::new("E", "Young's modulus", "Pa", PropKind::Scalar, &[], None),
            PropDesc::new("eps_r", "Permittivity", "", PropKind::Tensor, &["modes"], Some(1.0)),
        ]
    }
    fn params(&self) -> Vec<ParamDesc> {
        vec![ParamDesc::choice("plane", "Plane", &["strain", "stress"], "strain", "")]
    }
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![FieldDesc::new("E", "V/m", &["x", "y", "z"], "field")]
    }
    fn exports(&self) -> Vec<FieldLink> {
        self.exports.clone()
    }
    fn imports(&self) -> Vec<FieldLink> {
        self.imports.clone()
    }
    fn validate(&self, _spec: &ProblemSpec, _mesh: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, spec: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }
    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        ctx.check()?;
        Ok(SolveResult { physics: self.id.into(), ..Default::default() })
    }
}

pub struct Table;

impl LabelResolver for Table {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "core" => Ok(vec![2]),
            "clad" => Ok(vec![1]),
            other => Err(Error::invalid(format!("no region '{other}'"))),
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        if selector == "outer" {
            Ok(vec![1, 2, 3, 4])
        } else {
            Err(Error::invalid(format!("cannot read the selector '{selector}'")))
        }
    }
}

pub fn contract_file(name: &str) -> PathBuf {
    nyx_physics::contract::contract_dir(env!("CARGO_MANIFEST_DIR")).join(name)
}

pub fn check_exemplar(name: &str, got: &Value) {
    nyx_physics::contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), name, got);
}
