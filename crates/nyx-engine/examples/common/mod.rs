//! Shared by the benchmark programs: a plain generated mesh with numbered boundaries (`femcore/examples/benchmesh`) is driven through the
//! `Physics` trait, with a resolver that maps the boundary names used here to those numbers. The same kernels, on the same meshes, as the
//! window and the MCP server reach through `Engine`; only the geometry and meshing steps are left out so that the times are solver times.
#![allow(dead_code)]

use nyx_engine::facade::physics::{
    BcSpec, Ctx, Error, LabelResolver, MaterialAssign, Mesh, ProblemSpec, SolveResult, StudySpec,
};
use nyx_engine::facade::{Progress, SimplexMesh};
use serde_json::{json, Value};
use std::sync::Arc;

/// Boundary names and the numbers they stand for; every region name is the one domain 1 of the generated mesh.
pub struct Table(pub Vec<(&'static str, Vec<i32>)>);

impl LabelResolver for Table {
    fn domains(&self, _region: &str) -> Result<Vec<i32>, Error> {
        Ok(vec![1])
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        self.0
            .iter()
            .find(|(n, _)| *n == selector)
            .map(|(_, l)| l.clone())
            .ok_or_else(|| Error::invalid(format!("no boundary '{selector}'")))
    }
}

/// `left`, `bottom`, `top`, `right` of `benchmesh::rect`, and `xmin`, `ymin`, `zmin`, `zmax`, `ymax`, `xmax` of `benchmesh::block`.
pub fn table() -> Table {
    Table(vec![
        ("left", vec![1]),
        ("bottom", vec![2]),
        ("top", vec![3]),
        ("right", vec![4]),
        ("xmin", vec![1]),
        ("ymin", vec![2]),
        ("zmin", vec![3]),
        ("zmax", vec![4]),
        ("ymax", vec![5]),
        ("xmax", vec![6]),
    ])
}

/// A problem of physics `id` (`dfs`, `sms`, `ems`) with the study `study` and its parameters.
pub fn spec(id: &str, order: u8, study: &str, params: Value) -> ProblemSpec {
    let mut s = ProblemSpec::new(id, order, study);
    s.study.params = params;
    s
}

pub fn study(s: &mut ProblemSpec, kind: &str, params: Value) {
    s.study = StudySpec { kind: kind.into(), params };
}

/// Properties on the whole mesh.
pub fn props(s: &mut ProblemSpec, p: &[(&str, f64)]) {
    s.materials.insert("all".into(), MaterialAssign::Props(p.iter().map(|(k, v)| (k.to_string(), *v)).collect()));
}

/// Solid steel properties on the whole mesh.
pub fn steel(s: &mut ProblemSpec) {
    props(s, &[("E", 200e9), ("nu", 0.3), ("rho", 7850.0)]);
}

pub fn condition(s: &mut ProblemSpec, on: &str, kind: &str, params: Value) {
    s.bcs.push(BcSpec { on: on.into(), kind: kind.into(), params });
}

pub fn fixed(s: &mut ProblemSpec, on: &str) {
    condition(s, on, "fixed", json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }));
}

pub fn load(s: &mut ProblemSpec, on: &str, t: [f64; 3]) {
    condition(s, on, "load", json!({ "traction": t }));
}

/// The linear solver option of the solid and diffusion physics: `auto`, `direct` or `iterative`.
pub fn solver(s: &mut ProblemSpec, name: &str) {
    s.params = json!({ "solver": name });
}

/// Solves `spec` on `mesh` under `progress` (a fresh one when none is given).
pub fn solve(mesh: &SimplexMesh, spec: &ProblemSpec, progress: Option<Arc<Progress>>) -> Result<SolveResult, String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    let reg = nyx_engine::registry();
    let physics = reg.get(&spec.physics).ok_or_else(|| format!("no physics {}", spec.physics))?;
    let table = table();
    let mut ctx = Ctx::new(&table);
    if let Some(p) = progress {
        ctx.progress = p;
    }
    physics.solve(&exchange, spec, &mut ctx).map_err(|e| e.to_string())
}
