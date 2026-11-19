//! Whole-array operations on results (productivity plan N1, ADR 90): `field`, `field_op`, `integrate`, `extremum`, `probe`, `line_cut`,
//! `overlap` and `confinement` as engine commands. A script, the MCP server and the window's Results tab call these same commands.
//!
//! **Handles are values.** A handle is a small JSON tree, not an id the engine keeps: a leaf `{"field": "Ex", "physics": "ems", "solution": 0}`
//! (what `field` returns), a number, or `{"op": "mul", "args": [a, b]}` with the operations `add sub mul div neg abs abs2 re im conj sqrt pow`
//! (`pow` takes `"n"`, a whole number 1 to 16). A bare string stands for a leaf of the request's `physics` and `solution`. The tree is parsed
//! into `femcore::fieldexpr::Tree` and evaluated there, so a formula of N2 and a handle are one request with one evaluator. Because nothing is
//! kept, a handle cannot go stale in the engine; it is resolved against the stored results each time it is used, and the answer says whether
//! the result it read is `stale`.
//!
//! **What a leaf reads.** The un-averaged per-cell fields of a solution (`cell_fields`: `Ex`, `uz`, ..., complex when the solve is) when one has
//! the name; else a node-averaged field by its flat name (`T`, `V`, `|E|`; real, one value per dof). Both are Lagrange P1 or P2 interpolants
//! on the cells of the mesh the result was solved on (the result keeps that mesh, so a result is read on its own cells whatever the project's
//! mesh is now).
//!
//! **Quadrature (the rule S5 must reproduce bit for bit).** Integrals use the degree-4 symmetric rule of `femcore::fem::quad_rule` (6 points
//! on a triangle, 14 on a tetrahedron), exact for the product of two P2 fields, so `|E|^2` of a P2 field is exact up to rounding; an
//! integrand of higher degree (a product of three fields, `sqrt`, a quotient) carries the quadrature error. At each point every leaf is the
//! interpolant `sum_a N_a(q) v_a` (a in local node order, starting from 0.0) and the tree is evaluated on those values. A cell adds
//! `vol * sum_q (w_q * f_q)`; cells are summed in blocks of `BLOCK` consecutive cell numbers in cell order and the blocks in block order, so the
//! sum has the same bits for any thread count. With `revolution` the integrand at each point is multiplied by `2 pi r` (r is mesh axis 0; the
//! mesh must be 2D); the default is the plane measure, and an answer says which was used.
//!
//! **Maxima** are taken over the values at the nodes of the cells (the vertices, and the edge midpoints of a P2 mesh), not over the interpolant:
//! the true maximum of a P2 field between nodes can be larger. A node shared by cells has one value per cell (the fields are not averaged);
//! all of them count. **Points** (`probe`, `line_cut`) are in mesh axes (a 2D mesh: x then z of the wafer; 3D: x, y, z); a point outside the
//! mesh gives `inside: false` and a null value, not an error; a point on the face between cells is read in the cell with the lowest number.

use crate::engine::Engine;
use femcore::fem::{geometry, quad_rule, shape, RefElement};
use femcore::fieldexpr::{Fun, Op, Tree, MAX_POWER};
use femcore::formula::{self, Named, Names};
use femcore::linalg::csr::C;
use femcore::mesh::SimplexMesh;
use nyx_physics::{Error, LabelResolver, SolveResult};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::sync::Arc;

use crate::meshing::MeshedModel;

type R = Result<Value, Error>;

/// Consecutive cells summed together before the blocks are added in order.
pub const BLOCK: usize = 2048;
/// Barycentric slack for "the point is inside the cell".
const INSIDE: f64 = 1e-10;
const MAX_DEPTH: usize = 64;
const MAX_POINTS: usize = 200_000;

/// What a leaf reads (N2): a field, a coordinate, a material value or the indicator of a region. All of them live on the mesh of one result.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum What {
    Field,
    /// a wafer coordinate, 0 x, 1 y, 2 z; `name` is its letter
    Coord(usize),
    /// the value of a property of the material of each cell, under the descriptor name in `name`, at the wavelength given (m), else the study's
    Material(Option<f64>),
    /// 1 in the cells of the region `name`, else 0
    Region,
}

/// One leaf of an expression: a named field of one solution of one physics (or, with `what`, a coordinate, a material value, a region).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Spec {
    pub physics: String,
    pub solution: usize,
    pub name: String,
    pub what: What,
}

/// Where the arguments say a bare name comes from.
pub(crate) struct Defaults {
    physics: Option<String>,
    solution: usize,
    /// `wavelength` of the request: where the materials are read for a dispersive model (m); None is the wavelength of the study
    wavelength: Option<f64>,
}

enum Kind {
    /// index into `Solution::cell_fields`
    Cell(usize),
    /// the node-averaged values of the field, one per dof
    Dof(Vec<f64>),
    /// a wafer coordinate (0 x, 1 y, 2 z), interpolated from the coordinates of the nodes
    Coord(usize),
    /// one value per cell, the same at every point of the cell (a material value, a region indicator): not interpolated
    Const(Vec<f64>),
}

struct Leaf {
    res: Arc<SolveResult>,
    sol: usize,
    kind: Kind,
    complex: bool,
    unit: String,
}

/// The leaves of one request resolved against the stored results, and the mesh they live on.
pub(crate) struct Prep {
    leaves: Vec<Leaf>,
    model: Arc<MeshedModel>,
    mesh: Arc<SimplexMesh>,
    order: usize,
    nloc: usize,
    physics: Vec<String>,
}

fn invalid<T>(text: impl Into<String>) -> Result<T, Error> {
    Err(Error::invalid(text))
}

fn slot(specs: &mut Vec<Spec>, s: Spec) -> usize {
    match specs.iter().position(|x| *x == s) {
        Some(i) => i,
        None => {
            specs.push(s);
            specs.len() - 1
        }
    }
}

fn same_mesh(a: &SimplexMesh, b: &SimplexMesh) -> bool {
    a.dim == b.dim && a.points == b.points && a.cells == b.cells && a.cell_labels == b.cell_labels
}

fn interp(n: &[f64], loc: &[C]) -> C {
    let (mut re, mut im) = (0.0, 0.0);
    for a in 0..n.len() {
        re += n[a] * loc[a].re;
        im += n[a] * loc[a].im;
    }
    C::new(re, im)
}

/// What the names of a formula mean for this engine: coordinates, the fields of a result, the material values of its problem.
struct Namer<'a> {
    engine: &'a Engine,
    d: &'a Defaults,
    specs: &'a mut Vec<Spec>,
    /// an error that is not about the text (no result yet, no physics given), returned as it is
    missing: Option<Error>,
    /// a name that nothing in the result or the materials has
    unknown: bool,
}

/// `ems:Ex`, `ems:1:Ex` or `Ex`: the physics and solution a name is read from, and the name.
fn split_target(text: &str, d: &Defaults) -> Result<(Option<String>, usize, String), String> {
    let parts: Vec<&str> = text.split(':').collect();
    match parts.as_slice() {
        [n] => Ok((d.physics.clone(), d.solution, n.to_string())),
        [p, n] => Ok((Some(p.to_string()), d.solution, n.to_string())),
        [p, s, n] => {
            let s = s.trim().parse::<usize>().map_err(|_| format!("'{s}' in '{text}' is not a solution number"))?;
            Ok((Some(p.to_string()), s, n.to_string()))
        }
        _ => Err(format!("'{text}' has more than two ':' (it is physics:solution:name)")),
    }
}

/// The material property a formula name means (`eps` and `mu` are `eps_r` and `mu_r`).
fn prop_name(name: &str) -> &str {
    match name {
        "eps" => "eps_r",
        "mu" => "mu_r",
        other => other,
    }
}

impl Namer<'_> {
    fn slot(&mut self, physics: &str, solution: usize, name: &str, what: What) -> Tree {
        Tree::Leaf(slot(self.specs, Spec { physics: physics.to_string(), solution, name: name.to_string(), what }))
    }

    /// The physics and solution of a name, with the result they live in checked.
    fn target(&mut self, text: &str) -> Result<(String, usize, String), String> {
        let (physics, solution, name) = split_target(text, self.d)?;
        let Some(physics) = physics else {
            let e = Error::invalid(format!(
                "the name '{name}' needs a physics: give `physics` (the results are: {})",
                self.engine.results.keys().cloned().collect::<Vec<_>>().join(", ")
            ));
            let m = e.to_string();
            self.missing = Some(e);
            return Err(m);
        };
        let Some(stored) = self.engine.results.get(&physics) else {
            let e = Error::invalid(format!(
                "there is no result for '{physics}' yet; solve first{}",
                if self.engine.results.is_empty() {
                    String::new()
                } else {
                    format!(
                        " (the results are: {})",
                        self.engine.results.keys().cloned().collect::<Vec<_>>().join(", ")
                    )
                }
            ));
            let m = e.to_string();
            self.missing = Some(e);
            return Err(m);
        };
        if stored.result.solutions.get(solution).is_none() {
            let e = Error::invalid(format!(
                "the result of '{physics}' has {} solutions; there is no {solution}",
                stored.result.solutions.len()
            ));
            let m = e.to_string();
            self.missing = Some(e);
            return Err(m);
        }
        Ok((physics, solution, name))
    }

    fn props_of(&self, physics: &str) -> Vec<String> {
        crate::physics_for(&self.engine.registry, physics)
            .map(|p| p.properties().into_iter().map(|d| d.name).collect())
            .unwrap_or_default()
    }

    /// `eps_r_xz` and the like: a component of a tensor property of this physics, as the canonical name (`eps_r_xz`; `eps_zx` is the same).
    fn component(&self, physics: &str, name: &str) -> Option<String> {
        let (base, i, j) = tensor_component(name)?;
        let base = prop_name(base);
        self.props_of(physics).iter().any(|p| p == base).then(|| format!("{base}_{}{}", AXES[i.min(j)], AXES[i.max(j)]))
    }
}

const AXES: [&str; 3] = ["x", "y", "z"];

/// `eps_r_xz`: the property `eps_r` and the wafer axes of the entry, 0 x, 1 y, 2 z (the tensor is symmetric, so the order of the two letters
/// does not matter to the value). Only the tensor properties have components: the caller checks the base against the physics.
fn tensor_component(name: &str) -> Option<(&str, usize, usize)> {
    let (base, ij) = name.rsplit_once('_')?;
    let mut it = ij.chars();
    let (a, b) = (it.next()?, it.next()?);
    if it.next().is_some() {
        return None;
    }
    let axis = |c: char| "xyz".find(c);
    Some((base, axis(a)?, axis(b)?)).filter(|(b, _, _)| matches!(prop_name(b), "eps_r" | "k" | "alpha"))
}

impl Names for Namer<'_> {
    fn name(&mut self, text: &str) -> Result<Option<Named>, String> {
        let (physics, solution, name) = self.target(text)?;
        let sol = &self.engine.results[&physics].result.solutions[solution];
        if let Some(axis) = ["x", "y", "z"].iter().position(|a| *a == name) {
            return Ok(Some(Named::Scalar(self.slot(&physics, solution, &name, What::Coord(axis)))));
        }
        if Engine::field_names(sol).contains(&name) {
            return Ok(Some(Named::Scalar(self.slot(&physics, solution, &name, What::Field))));
        }
        let comps: Vec<(usize, String)> = ["x", "y", "z"]
            .iter()
            .enumerate()
            .map(|(a, c)| (a, format!("{name}{c}")))
            .filter(|(_, n)| sol.cell_fields.iter().any(|f| f.name == *n))
            .collect();
        if comps.len() >= 2 {
            return Ok(Some(Named::Vector(
                comps.into_iter().map(|(a, n)| (a, self.slot(&physics, solution, &n, What::Field))).collect(),
            )));
        }
        let prop = prop_name(&name);
        if self.props_of(&physics).iter().any(|p| p == prop) {
            return Ok(Some(Named::Scalar(self.slot(&physics, solution, prop, What::Material(self.d.wavelength)))));
        }
        if let Some(c) = self.component(&physics, &name) {
            return Ok(Some(Named::Scalar(self.slot(&physics, solution, &c, What::Material(self.d.wavelength)))));
        }
        self.unknown = true;
        Ok(None)
    }

    fn material(&mut self, text: &str) -> Result<Tree, String> {
        let (physics, solution, name) = self.target(text)?;
        let prop = prop_name(&name).to_string();
        let props = self.props_of(&physics);
        if let Some(c) = self.component(&physics, &name) {
            return Ok(self.slot(&physics, solution, &c, What::Material(self.d.wavelength)));
        }
        if !props.contains(&prop) {
            return Err(format!("'{physics}' has no material property '{name}' (it has: {})", props.join(", ")));
        }
        Ok(self.slot(&physics, solution, &prop, What::Material(self.d.wavelength)))
    }

    fn region(&mut self, text: &str) -> Result<Tree, String> {
        let (physics, solution, name) = self.target(text)?;
        Ok(self.slot(&physics, solution, &name, What::Region))
    }

    fn known(&self) -> String {
        let mut v = vec!["x".to_string(), "y".to_string(), "z".to_string()];
        if let Some(sol) = self
            .d
            .physics
            .as_deref()
            .and_then(|p| self.engine.results.get(p))
            .and_then(|s| s.result.solutions.get(self.d.solution))
        {
            v.extend(Engine::field_names(sol));
        }
        if let Some(p) = self.d.physics.as_deref() {
            for n in self.props_of(p) {
                if matches!(n.as_str(), "eps_r" | "k" | "alpha") {
                    v.push(format!("{n} ({n}_xx ... {n}_zz for the tensor)"));
                } else {
                    v.push(n);
                }
            }
        }
        v.join(", ")
    }
}

// ------------------------------------------------------------------------------------------------------------- the engine side

impl Engine {
    fn defaults(&self, args: &Value) -> Result<Defaults, Error> {
        let physics = match args.get("physics") {
            None | Some(Value::Null) => {
                (self.results.len() == 1).then(|| self.results.keys().next().cloned()).flatten()
            }
            Some(Value::String(s)) => Some(s.clone()),
            Some(o) => return invalid(format!("`physics` must be text, not {o}")),
        };
        let solution = match args.get("solution") {
            None | Some(Value::Null) => 0,
            Some(v) => v.as_u64().ok_or_else(|| Error::invalid("`solution` must be a whole number"))? as usize,
        };
        let wavelength = match args.get("wavelength") {
            None | Some(Value::Null) => None,
            Some(v) => match v.as_f64() {
                Some(w) if w > 0.0 && w.is_finite() => Some(w),
                _ => return invalid("`wavelength` must be a positive number of metres"),
            },
        };
        Ok(Defaults { physics, solution, wavelength })
    }

    fn parse_expr(&self, v: &Value, d: &Defaults, specs: &mut Vec<Spec>, depth: usize) -> Result<Tree, Error> {
        if depth > MAX_DEPTH {
            return invalid("the expression is nested too deeply");
        }
        let leaf = |name: &str, physics: Option<&str>, solution: usize, specs: &mut Vec<Spec>| -> Result<Tree, Error> {
            let physics = physics.ok_or_else(|| {
                Error::invalid(format!(
                    "the field '{name}' needs a physics: give `physics` (the results are: {})",
                    self.results.keys().cloned().collect::<Vec<_>>().join(", ")
                ))
            })?;
            Ok(Tree::Leaf(slot(
                specs,
                Spec { physics: physics.to_string(), solution, name: name.to_string(), what: What::Field },
            )))
        };
        match v {
            Value::Number(n) => Ok(Tree::num(n.as_f64().ok_or_else(|| Error::invalid("a number is out of range"))?)),
            Value::String(s) => self.parse_text(s, d, specs),
            Value::Object(m) => {
                if let Some(f) = m.get("field") {
                    let name = f.as_str().ok_or_else(|| Error::invalid("`field` must be the name of a field"))?;
                    let physics = match m.get("physics") {
                        Some(Value::String(p)) => Some(p.as_str()),
                        _ => d.physics.as_deref(),
                    };
                    let solution = match m.get("solution") {
                        Some(s) => {
                            s.as_u64().ok_or_else(|| Error::invalid("`solution` must be a whole number"))? as usize
                        }
                        None => d.solution,
                    };
                    return leaf(name, physics, solution, specs);
                }
                if let Some(n) = m.get("num") {
                    return match n {
                        Value::Number(x) => Ok(Tree::num(x.as_f64().unwrap_or(f64::NAN))),
                        Value::Array(a) if a.len() == 2 => match (a[0].as_f64(), a[1].as_f64()) {
                            (Some(re), Some(im)) => Ok(Tree::Num(C::new(re, im))),
                            _ => invalid("`num` pairs are [re, im] numbers"),
                        },
                        _ => invalid("`num` is a number or a pair [re, im]"),
                    };
                }
                let op = m
                    .get("op")
                    .and_then(Value::as_str)
                    .ok_or_else(|| Error::invalid("an expression object has `field`, `num` or `op`"))?;
                let args: &[Value] = m
                    .get("args")
                    .and_then(Value::as_array)
                    .map(Vec::as_slice)
                    .ok_or_else(|| Error::invalid(format!("the operation '{op}' needs `args`, a list")))?;
                let mut sub = Vec::with_capacity(args.len());
                for a in args {
                    sub.push(self.parse_expr(a, d, specs, depth + 1)?);
                }
                let one = |sub: Vec<Tree>| -> Result<Tree, Error> {
                    if sub.len() != 1 {
                        return invalid(format!("'{op}' takes one argument, not {}", sub.len()));
                    }
                    Ok(sub.into_iter().next().expect("one argument"))
                };
                let two = |sub: Vec<Tree>| -> Result<(Tree, Tree), Error> {
                    if sub.len() != 2 {
                        return invalid(format!("'{op}' takes two arguments, not {}", sub.len()));
                    }
                    let mut it = sub.into_iter();
                    Ok((it.next().expect("two arguments"), it.next().expect("two arguments")))
                };
                let fold = |op_: Op, sub: Vec<Tree>| -> Result<Tree, Error> {
                    if sub.len() < 2 {
                        return invalid(format!("'{op}' takes at least two arguments"));
                    }
                    let mut it = sub.into_iter();
                    let first = it.next().expect("an argument");
                    Ok(it.fold(first, |acc, t| Tree::Bin(op_, Box::new(acc), Box::new(t))))
                };
                match op {
                    "add" => fold(Op::Add, sub),
                    "mul" => fold(Op::Mul, sub),
                    "sub" | "div" => {
                        let (a, b) = two(sub)?;
                        Ok(Tree::Bin(if op == "sub" { Op::Sub } else { Op::Div }, Box::new(a), Box::new(b)))
                    }
                    "neg" => Ok(Tree::Neg(Box::new(one(sub)?))),
                    "abs" | "abs2" | "re" | "im" | "conj" | "sqrt" => {
                        let f = match op {
                            "abs" => Fun::Abs,
                            "abs2" => Fun::Abs2,
                            "re" => Fun::Re,
                            "im" => Fun::Im,
                            "conj" => Fun::Conj,
                            _ => Fun::Sqrt,
                        };
                        Ok(Tree::Call(f, Box::new(one(sub)?)))
                    }
                    "pow" => {
                        let n = m
                            .get("n")
                            .and_then(Value::as_u64)
                            .filter(|n| (1..=MAX_POWER as u64).contains(n))
                            .ok_or_else(|| {
                                Error::invalid(format!("`pow` needs `n`, a whole number from 1 to {MAX_POWER}"))
                            })?;
                        Ok(Tree::Pow(Box::new(one(sub)?), n as u32))
                    }
                    other => invalid(format!(
                        "no operation '{other}' (they are: add, sub, mul, div, neg, abs, abs2, re, im, conj, sqrt, pow)"
                    )),
                }
            }
            other => invalid(format!("an expression is a field name, a number or an object, not {other}")),
        }
    }

    /// A string in an expression: the exact name of a field (`Ex`, `|E|`, `Re Ex`: what a bare string always meant), else a formula
    /// (`femcore::formula`, N2). A string that is not a formula either and reads like a name is reported as a missing field.
    fn parse_text(&self, text: &str, d: &Defaults, specs: &mut Vec<Spec>) -> Result<Tree, Error> {
        let plain = |s: &str| {
            !s.chars().any(|c| "+-*/^()<>=!,[]\"'`.:".contains(c)) && !s.contains(" and ") && !s.contains(" or ")
        };
        if let Some(p) = d.physics.as_deref() {
            let named = self
                .results
                .get(p)
                .and_then(|st| st.result.solutions.get(d.solution))
                .is_some_and(|sol| Self::field_names(sol).iter().any(|n| n == text));
            if named {
                return Ok(Tree::Leaf(slot(
                    specs,
                    Spec { physics: p.to_string(), solution: d.solution, name: text.to_string(), what: What::Field },
                )));
            }
        }
        let mut namer = Namer { engine: self, d, specs, missing: None, unknown: false };
        match formula::parse(text, &mut namer) {
            Ok(t) => Ok(t),
            Err(e) => {
                if let Some(m) = namer.missing.take() {
                    return Err(m);
                }
                if namer.unknown && plain(text.trim()) {
                    // a name that is no field and no formula: say which fields there are, as a missing field always did
                    let physics = d.physics.as_deref().ok_or_else(|| {
                        Error::invalid(format!(
                            "the field '{text}' needs a physics: give `physics` (the results are: {})",
                            self.results.keys().cloned().collect::<Vec<_>>().join(", ")
                        ))
                    })?;
                    namer.specs.push(Spec {
                        physics: physics.to_string(),
                        solution: d.solution,
                        name: text.to_string(),
                        what: What::Field,
                    });
                    return Ok(Tree::Leaf(namer.specs.len() - 1));
                }
                Err(Error::invalid(format!("the formula '{text}' does not read: {e}")))
            }
        }
    }

    /// The value of a material property in each cell of the mesh a result was solved on: the problem's materials as the solver saw them
    /// (library or inline numbers, the optical model at the wavelength of the study), under the descriptor name (`eps_r`, `rho`, ...).
    /// A property that is a tensor in some material of the mesh is refused rather than read as its scalar default.
    fn material_values(
        &self,
        s: &Spec,
        wavelength: Option<f64>,
        model: &Arc<MeshedModel>,
        mesh: &SimplexMesh,
    ) -> Result<(Vec<f64>, String), Error> {
        let p = crate::physics_for(&self.registry, &s.physics)?;
        let spec = self.project.spec(&s.physics).ok_or_else(|| {
            Error::invalid(format!("the project has no problem for '{}', so its materials are not known", s.physics))
        })?;
        let props = p.properties();
        let component = tensor_component(&s.name);
        let base = component.map_or(s.name.as_str(), |(b, _, _)| b);
        let desc = props.iter().find(|d| d.name == base).ok_or_else(|| {
            Error::invalid(format!(
                "'{}' has no material property '{}' (it has: {})",
                s.physics,
                s.name,
                props.iter().map(|d| d.name.clone()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let mut ctx = nyx_physics::Ctx::new(model.as_ref());
        ctx.library = Some(self.library()?);
        ctx.wavelength = wavelength.or_else(|| crate::materials::study_wavelength(p, spec));
        let m = nyx_physics::decode::materials_by_label(spec, &ctx, &props)?;
        let kind = match base {
            "eps_r" => Some(femcore::props::TensorKind::Permittivity),
            "k" => Some(femcore::props::TensorKind::Conductivity),
            "alpha" => Some(femcore::props::TensorKind::Expansion),
            _ => None,
        };
        if let (Some((_, i, j)), Some(kind)) = (component, kind) {
            // the tensor the solver used, in wafer axes (x, y, z) with the crystal orientation applied; an isotropic material is `v I`
            let t = m.cell_tensor(mesh, kind, desc.default).map_err(Error::invalid)?;
            return Ok((t.iter().map(|c| c[i][j]).collect(), desc.unit.clone()));
        }
        let tensor = |a: &femcore::material::Anisotropy| match s.name.as_str() {
            "eps_r" => a.permittivity.is_some(),
            "k" => a.conductivity.is_some(),
            "alpha" => a.expansion.is_some(),
            _ => false,
        };
        if mesh.domain_labels().iter().any(|l| m.aniso.get(l).is_some_and(tensor)) {
            return invalid(format!(
                "'{n}' is a tensor in a material of this mesh (an anisotropic or rotated crystal): read its entries {n}_xx, {n}_yy, {n}_zz, {n}_xy, {n}_xz, {n}_yz (wafer axes)",
                n = s.name
            ));
        }
        let v = m.cell_values(mesh, &s.name, desc.default).map_err(Error::invalid)?;
        Ok((v, desc.unit.clone()))
    }

    /// The names a solution offers: its cell fields and its flat fields.
    fn field_names(sol: &nyx_physics::Solution) -> Vec<String> {
        let mut v: Vec<String> = sol.cell_fields.iter().map(|f| f.name.clone()).collect();
        for (n, _, _) in sol.flat_fields() {
            if !v.contains(&n) {
                v.push(n);
            }
        }
        v
    }

    /// Resolves the leaves against the stored results.
    pub(crate) fn prepare(&self, specs: &[Spec]) -> Result<Prep, Error> {
        let mut base: Option<(Arc<MeshedModel>, Arc<SimplexMesh>, usize)> = None;
        let mut leaves = Vec::with_capacity(specs.len());
        let mut physics: Vec<String> = vec![];
        for s in specs {
            let stored = self.results.get(&s.physics).ok_or_else(|| {
                Error::invalid(format!(
                    "there is no result for '{}' yet; solve first{}",
                    s.physics,
                    if self.results.is_empty() {
                        String::new()
                    } else {
                        format!(" (the results are: {})", self.results.keys().cloned().collect::<Vec<_>>().join(", "))
                    }
                ))
            })?;
            if !physics.contains(&s.physics) {
                physics.push(s.physics.clone());
            }
            let res = &stored.result;
            let sol = res.solutions.get(s.solution).ok_or_else(|| {
                Error::invalid(format!(
                    "the result of '{}' has {} solutions; there is no {}",
                    s.physics,
                    res.solutions.len(),
                    s.solution
                ))
            })?;
            let mesh = stored.model.solver_mesh();
            let order = res.order as usize;
            let nloc = mesh.nloc(order);
            let (kind, complex, unit) = match &s.what {
                What::Coord(a) => (Kind::Coord(*a), false, "m".to_string()),
                What::Region => {
                    let labels = stored.model.domains(&s.name)?;
                    let v: Vec<f64> =
                        mesh.cell_labels.iter().map(|l| if labels.contains(l) { 1.0 } else { 0.0 }).collect();
                    if !v.contains(&1.0) {
                        return invalid(format!("the region '{}' has no cells", s.name));
                    }
                    (Kind::Const(v), false, String::new())
                }
                What::Material(wavelength) => {
                    let (v, unit) = self.material_values(s, *wavelength, &stored.model, &mesh)?;
                    (Kind::Const(v), false, unit)
                }
                What::Field => {
                    if let Some(i) = sol.cell_fields.iter().position(|f| f.name == s.name) {
                        let f = &sol.cell_fields[i];
                        if f.nloc != nloc
                            || f.re.len() != mesh.n_cells() * nloc
                            || (!f.im.is_empty() && f.im.len() != f.re.len())
                        {
                            return invalid(format!(
                                "the field '{}' of '{}' does not fit its mesh ({} values for {} cells of {} nodes)",
                                s.name,
                                s.physics,
                                f.re.len(),
                                mesh.n_cells(),
                                nloc
                            ));
                        }
                        (Kind::Cell(i), !f.im.is_empty(), f.unit.clone())
                    } else if let Some((_, unit, vals)) = sol.flat_fields().into_iter().find(|(n, _, _)| *n == s.name) {
                        if vals.len() != mesh.n_dofs(order) {
                            return invalid(format!(
                                "the field '{}' of '{}' has {} values for {} dofs",
                                s.name,
                                s.physics,
                                vals.len(),
                                mesh.n_dofs(order)
                            ));
                        }
                        (Kind::Dof(vals), false, unit)
                    } else {
                        return invalid(format!(
                            "'{}' has no field '{}' in solution {} (the fields are: {})",
                            s.physics,
                            s.name,
                            s.solution,
                            Self::field_names(sol).join(", ")
                        ));
                    }
                }
            };
            match &base {
                None => base = Some((stored.model.clone(), mesh, order)),
                Some((m, b, o)) => {
                    if *o != order {
                        return invalid(format!(
                            "'{}' is solved at order {order} and an earlier field at order {o}; an expression reads fields of one order",
                            s.physics
                        ));
                    }
                    if !Arc::ptr_eq(m, &stored.model) && !same_mesh(b, &mesh) {
                        return invalid(format!(
                            "'{}' was solved on a different mesh from the other fields of the expression; solve them on one mesh",
                            s.physics
                        ));
                    }
                }
            }
            leaves.push(Leaf { res: res.clone(), sol: s.solution, kind, complex, unit });
        }
        let (model, mesh, order) = base.ok_or_else(|| Error::invalid("the expression reads no field"))?;
        let nloc = mesh.nloc(order);
        Ok(Prep { leaves, model, mesh, order, nloc, physics })
    }

    /// Parses expressions that share one leaf table (each a handle, a name or a number), resolved against the stored results.
    pub(crate) fn expressions(&self, args: &Value, exprs: &[&Value]) -> Result<(Vec<Tree>, Prep), Error> {
        let d = self.defaults(args)?;
        let mut specs = vec![];
        let mut trees = vec![];
        for e in exprs {
            trees.push(self.parse_expr(e, &d, &mut specs, 0)?);
        }
        // an expression of numbers only (`1` for an area) still needs a mesh to integrate over: the one of the default result
        if specs.is_empty() {
            if let Some(p) = &d.physics {
                specs.push(Spec { physics: p.clone(), solution: d.solution, name: "x".into(), what: What::Coord(0) });
            }
        }
        let prep = self.prepare(&specs)?;
        let flags = prep.complex_flags();
        for (t, e) in trees.iter().zip(exprs) {
            t.check(&flags).map_err(|m| match e {
                Value::String(s) => Error::invalid(format!("the formula '{s}': {m}")),
                _ => Error::invalid(m),
            })?;
        }
        Ok((trees, prep))
    }

    fn one_expr(&self, args: &Value, key: &str) -> Result<(Tree, Prep), Error> {
        let e = args
            .get(key)
            .ok_or_else(|| Error::invalid(format!("`{key}` is required: a field name, a number or an expression")))?;
        let (mut t, p) = self.expressions(args, &[e])?;
        Ok((t.remove(0), p))
    }

    /// Adds the staleness of what an answer read.
    fn decorate(&self, v: &mut Value, p: &Prep) {
        self.decorate_physics(v, &p.physics);
    }

    fn decorate_physics(&self, v: &mut Value, physics: &[String]) {
        let stale = physics.iter().any(|k| self.is_stale(k));
        let age = physics
            .iter()
            .filter_map(|k| self.result_age(k))
            .fold(None, |m: Option<f64>, a| Some(m.map_or(a, |m| m.max(a))));
        v["physics"] = json!(physics);
        v["stale"] = json!(stale);
        v["age_s"] = json!(age);
        v["validated"] = self.validated_of(physics);
    }

    /// The names of the regions of the mesh a result was solved on (what the Results tab offers to reduce over).
    pub fn result_regions(&self, physics: &str) -> Vec<String> {
        self.results.get(physics).map(|s| s.model.regions().into_iter().map(|(_, n)| n).collect()).unwrap_or_default()
    }
    // ------------------------------------------------------------------------------------------------------------- commands

    /// `field`: a handle to a field of a solution, with what it holds. A vector name (`E` when `Ex`, `Ez` exist) answers one handle per component.
    pub(crate) fn cmd_field(&self, args: &Value) -> R {
        let d = self.defaults(args)?;
        let name = args
            .get("name")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid("`name` is required: the field, e.g. \"Ex\""))?;
        let physics = d.physics.clone().ok_or_else(|| {
            Error::invalid(format!(
                "give `physics` (the results are: {})",
                self.results.keys().cloned().collect::<Vec<_>>().join(", ")
            ))
        })?;
        let handle = |n: &str| json!({ "field": n, "physics": physics, "solution": d.solution });
        let describe = |n: &str| -> Result<Value, Error> {
            let spec = Spec { physics: physics.clone(), solution: d.solution, name: n.to_string(), what: What::Field };
            let p = self.prepare(std::slice::from_ref(&spec))?;
            let l = &p.leaves[0];
            Ok(json!({
                "handle": handle(n), "name": n, "unit": l.unit, "complex": l.complex,
                "source": match l.kind { Kind::Cell(_) => "cell", Kind::Dof(_) => "dof", _ => "formula" },
                "order": p.order, "cells": p.mesh.n_cells(), "nodes_per_cell": p.nloc,
            }))
        };
        let stored = self
            .results
            .get(&physics)
            .ok_or_else(|| Error::invalid(format!("there is no result for '{physics}' yet; solve first")))?;
        let sol = stored.result.solutions.get(d.solution).ok_or_else(|| {
            Error::invalid(format!(
                "the result has {} solutions; there is no {}",
                stored.result.solutions.len(),
                d.solution
            ))
        })?;
        let names = Self::field_names(sol);
        let mut out = if names.iter().any(|n| n == name) {
            describe(name)?
        } else {
            let comps: Vec<(String, String)> = ["x", "y", "z"]
                .iter()
                .map(|c| (c.to_string(), format!("{name}{c}")))
                .filter(|(_, n)| sol.cell_fields.iter().any(|f| f.name == *n))
                .collect();
            if comps.len() < 2 {
                return invalid(format!(
                    "'{physics}' has no field '{name}' in solution {} (the fields are: {})",
                    d.solution,
                    names.join(", ")
                ));
            }
            let mut m = serde_json::Map::new();
            for (c, n) in &comps {
                m.insert(c.clone(), describe(n)?);
            }
            json!({ "name": name, "components": Value::Object(m) })
        };
        self.decorate_physics(&mut out, std::slice::from_ref(&physics));
        Ok(out)
    }

    /// `field_op`: arithmetic on handles. The operation is checked against the results (every leaf must resolve) and the handle comes back.
    pub(crate) fn cmd_field_op(&self, args: &Value) -> R {
        let op = args.get("op").cloned().ok_or_else(|| Error::invalid("`op` is required"))?;
        let list = args.get("args").cloned().ok_or_else(|| Error::invalid("`args` is required: the operands"))?;
        let mut e = json!({ "op": op, "args": list });
        if let Some(n) = args.get("n") {
            e["n"] = n.clone();
        }
        let (tree, prep) = self.expressions(args, &[&e])?;
        let flags: Vec<bool> = prep.leaves.iter().map(|l| l.complex).collect();
        let mut out = json!({ "handle": e, "complex": !tree[0].is_real_valued(&flags) });
        self.decorate(&mut out, &prep);
        Ok(out)
    }
}

// ------------------------------------------------------------------------------------------------------------- evaluation

struct Region {
    mask: Vec<bool>,
    cells: usize,
}

struct Sum {
    re: f64,
    im: f64,
    measure: f64,
}

impl Prep {
    fn region(&self, region: Option<&str>) -> Result<Region, Error> {
        let n = self.mesh.n_cells();
        let Some(text) = region.filter(|r| !r.trim().is_empty()) else {
            return Ok(Region { mask: vec![true; n], cells: n });
        };
        let labels = self.model.domains(text)?;
        let mask: Vec<bool> = self.mesh.cell_labels.iter().map(|l| labels.contains(l)).collect();
        let cells = mask.iter().filter(|m| **m).count();
        if cells == 0 {
            return invalid(format!("the region '{text}' has no cells"));
        }
        Ok(Region { mask, cells })
    }

    /// The values of every leaf at the nodes of cell `c`: `loc[i * nloc + a]`.
    fn gather(&self, c: usize, dofs: &mut Vec<usize>, loc: &mut [C]) {
        let nloc = self.nloc;
        let need_dofs = self.leaves.iter().any(|l| matches!(l.kind, Kind::Dof(_)));
        if need_dofs {
            self.mesh.cell_dofs(c, self.order, dofs);
        }
        for (i, l) in self.leaves.iter().enumerate() {
            let out = &mut loc[i * nloc..(i + 1) * nloc];
            match &l.kind {
                Kind::Cell(k) => {
                    let f = &l.res.solutions[l.sol].cell_fields[*k];
                    let re = &f.re[c * nloc..(c + 1) * nloc];
                    if f.im.is_empty() {
                        for a in 0..nloc {
                            out[a] = C::new(re[a], 0.0);
                        }
                    } else {
                        let im = &f.im[c * nloc..(c + 1) * nloc];
                        for a in 0..nloc {
                            out[a] = C::new(re[a], im[a]);
                        }
                    }
                }
                Kind::Dof(vals) => {
                    for a in 0..nloc {
                        out[a] = C::new(vals[dofs[a]], 0.0);
                    }
                }
                Kind::Coord(axis) => {
                    for (a, o) in out.iter_mut().enumerate() {
                        *o = C::new(self.node_xyz(c, a)[*axis], 0.0);
                    }
                }
                Kind::Const(vals) => {
                    for o in out.iter_mut() {
                        *o = C::new(vals[c], 0.0);
                    }
                }
            }
        }
    }

    /// The value of every leaf at one point of the cell: the interpolant of its nodal values, or the cell's own value for a constant one.
    fn at_point(&self, tree: &Tree, loc: &[C], n: &[f64], work: &mut Vec<C>) -> C {
        work.clear();
        for (i, l) in self.leaves.iter().enumerate() {
            let nodes = &loc[i * self.nloc..(i + 1) * self.nloc];
            work.push(if matches!(l.kind, Kind::Const(_)) { nodes[0] } else { interp(n, nodes) });
        }
        tree.eval(work)
    }

    /// The wafer coordinates `(x, y, z)` of local node `a` of cell `c` (the corners, then the edge midpoints of a P2 cell).
    fn node_xyz(&self, c: usize, a: usize) -> [f64; 3] {
        let dim = self.mesh.dim;
        let v = self.mesh.cell(c);
        let mut p = [0.0; 3];
        if a <= dim {
            p[..dim].copy_from_slice(self.mesh.point(v[a] as usize));
        } else {
            let mut k = dim + 1;
            'find: for i in 0..=dim {
                for j in i + 1..=dim {
                    if k == a {
                        let (q, r) = (self.mesh.point(v[i] as usize), self.mesh.point(v[j] as usize));
                        for d in 0..dim {
                            p[d] = 0.5 * (q[d] + r[d]);
                        }
                        break 'find;
                    }
                    k += 1;
                }
            }
        }
        femcore::section::wafer_point(&p[..dim], dim)
    }

    fn complex_flags(&self) -> Vec<bool> {
        self.leaves.iter().map(|l| l.complex).collect()
    }

    /// `sum over the region of tree dV` (see the module comment for the rule and the order of the sums).
    fn integrate(&self, tree: &Tree, region: &Region, revolution: bool) -> Result<Sum, Error> {
        let dim = self.mesh.dim;
        if revolution && dim != 2 {
            return invalid("`revolution` applies to a 2D (r, z) mesh");
        }
        femcore::par::threads();
        let el = RefElement::new(dim, self.order, 4);
        let quad = quad_rule(dim, 4);
        let geo = geometry(&self.mesh);
        let nq = el.w.len();
        let nloc = self.nloc;
        let nb = self.mesh.n_cells().div_ceil(BLOCK);
        let block = |b: usize| -> Sum {
            let mut s = Sum { re: 0.0, im: 0.0, measure: 0.0 };
            let mut loc = vec![C::new(0.0, 0.0); self.leaves.len() * nloc];
            let mut dofs = vec![];
            let mut work = Vec::with_capacity(self.leaves.len());
            for c in b * BLOCK..((b + 1) * BLOCK).min(self.mesh.n_cells()) {
                if !region.mask[c] {
                    continue;
                }
                self.gather(c, &mut dofs, &mut loc);
                let (mut cre, mut cim, mut cm) = (0.0, 0.0, 0.0);
                let v = self.mesh.cell(c);
                for q in 0..nq {
                    let mut f = self.at_point(tree, &loc, &el.n[q], &mut work);
                    let mut w = el.w[q];
                    if revolution {
                        let mut r = 0.0;
                        for (k, &vk) in v.iter().enumerate() {
                            r += quad.lam[q][k] * self.mesh.point(vk as usize)[0];
                        }
                        let rho = 2.0 * std::f64::consts::PI * r;
                        f = C::new(f.re * rho, f.im * rho);
                        w *= rho;
                    }
                    cre += el.w[q] * f.re;
                    cim += el.w[q] * f.im;
                    cm += w;
                }
                s.re += geo.vol[c] * cre;
                s.im += geo.vol[c] * cim;
                s.measure += geo.vol[c] * cm;
            }
            s
        };
        let parts: Vec<Sum> = (0..nb).into_par_iter().map(block).collect();
        let mut t = Sum { re: 0.0, im: 0.0, measure: 0.0 };
        for p in parts {
            t.re += p.re;
            t.im += p.im;
            t.measure += p.measure;
        }
        if !(t.re.is_finite() && t.im.is_finite()) {
            return Err(Error::NonFinite(
                "the integral is not finite (a division by zero, or a field with NaN)".into(),
            ));
        }
        Ok(t)
    }

    /// The coordinates of local node `a` of cell `c`.
    fn node_point(&self, c: usize, a: usize) -> Vec<f64> {
        let dim = self.mesh.dim;
        let v = self.mesh.cell(c);
        if a <= dim {
            return self.mesh.point(v[a] as usize).to_vec();
        }
        let mut k = dim + 1;
        for i in 0..=dim {
            for j in i + 1..=dim {
                if k == a {
                    let (p, q) = (self.mesh.point(v[i] as usize), self.mesh.point(v[j] as usize));
                    return (0..dim).map(|d| 0.5 * (p[d] + q[d])).collect();
                }
                k += 1;
            }
        }
        unreachable!("local node {a} of a cell with {} nodes", self.nloc)
    }

    /// The largest or smallest value at the nodes of the cells of the region: (value, cell, local node).
    fn extremum(&self, tree: &Tree, region: &Region, max: bool) -> Result<(f64, usize, usize), Error> {
        let nloc = self.nloc;
        let nb = self.mesh.n_cells().div_ceil(BLOCK);
        let better = |a: f64, b: f64| if max { a > b } else { a < b };
        let block = |b: usize| -> Option<(f64, usize, usize)> {
            let mut best: Option<(f64, usize, usize)> = None;
            let mut loc = vec![C::new(0.0, 0.0); self.leaves.len() * nloc];
            let mut dofs = vec![];
            let mut at = vec![C::new(0.0, 0.0); self.leaves.len()];
            for c in b * BLOCK..((b + 1) * BLOCK).min(self.mesh.n_cells()) {
                if !region.mask[c] {
                    continue;
                }
                self.gather(c, &mut dofs, &mut loc);
                for a in 0..nloc {
                    for (i, x) in at.iter_mut().enumerate() {
                        *x = loc[i * nloc + a];
                    }
                    let v = tree.eval(&at).re;
                    if v.is_nan() {
                        return Some((f64::NAN, c, a));
                    }
                    if best.is_none_or(|(bv, _, _)| better(v, bv)) {
                        best = Some((v, c, a));
                    }
                }
            }
            best
        };
        let parts: Vec<Option<(f64, usize, usize)>> = (0..nb).into_par_iter().map(block).collect();
        let mut best: Option<(f64, usize, usize)> = None;
        for p in parts.into_iter().flatten() {
            if p.0.is_nan() {
                return Err(Error::NonFinite("the expression is NaN at a node".into()));
            }
            if best.is_none_or(|(bv, _, _)| better(p.0, bv)) {
                best = Some(p);
            }
        }
        best.ok_or_else(|| Error::invalid("there is nothing to compare"))
    }
}

/// A uniform grid of the cells' bounding boxes, to find the cell that holds a point.
struct Locator {
    dim: usize,
    lo: [f64; 3],
    h: [f64; 3],
    n: [usize; 3],
    start: Vec<u32>,
    items: Vec<u32>,
}

impl Locator {
    fn new(mesh: &SimplexMesh) -> Locator {
        let dim = mesh.dim;
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in mesh.points.chunks(dim) {
            for k in 0..dim {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let nc = mesh.n_cells();
        let per_axis =
            (((nc as f64) / 2.0).powf(1.0 / dim as f64).ceil() as usize).clamp(1, if dim == 3 { 160 } else { 2048 });
        let mut n = [1usize; 3];
        let mut h = [1.0f64; 3];
        for k in 0..dim {
            n[k] = per_axis;
            let ext = hi[k] - lo[k];
            h[k] = if ext > 0.0 { ext / per_axis as f64 } else { 1.0 };
        }
        let total: usize = n.iter().product();
        let range = |c: usize| -> ([usize; 3], [usize; 3]) {
            let (mut a, mut b) = ([0usize; 3], [0usize; 3]);
            let (mut cl, mut ch) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for &v in mesh.cell(c) {
                let p = mesh.point(v as usize);
                for k in 0..dim {
                    cl[k] = cl[k].min(p[k]);
                    ch[k] = ch[k].max(p[k]);
                }
            }
            for k in 0..dim {
                a[k] = (((cl[k] - lo[k]) / h[k]).floor().max(0.0) as usize).min(n[k] - 1);
                b[k] = (((ch[k] - lo[k]) / h[k]).floor().max(0.0) as usize).min(n[k] - 1);
            }
            (a, b)
        };
        let bins_of = |a: [usize; 3], b: [usize; 3], f: &mut dyn FnMut(usize)| {
            for i in a[0]..=b[0] {
                for j in a[1]..=b[1] {
                    for l in a[2]..=b[2] {
                        f((i * n[1] + j) * n[2] + l);
                    }
                }
            }
        };
        let mut count = vec![0u32; total + 1];
        for c in 0..nc {
            let (a, b) = range(c);
            bins_of(a, b, &mut |bin| count[bin + 1] += 1);
        }
        for i in 0..total {
            count[i + 1] += count[i];
        }
        let mut cursor = count.clone();
        let mut items = vec![0u32; count[total] as usize];
        for c in 0..nc {
            let (a, b) = range(c);
            bins_of(a, b, &mut |bin| {
                items[cursor[bin] as usize] = c as u32;
                cursor[bin] += 1;
            });
        }
        Locator { dim, lo, h, n, start: count, items }
    }

    /// The cells whose boxes touch the bin of `p` (ascending), empty when `p` is outside the grid.
    fn candidates(&self, p: &[f64]) -> &[u32] {
        let mut idx = [0usize; 3];
        for k in 0..self.dim {
            let t = (p[k] - self.lo[k]) / self.h[k];
            let top = self.n[k] as f64;
            if t < -1e-9 || t > top + 1e-9 {
                return &[];
            }
            idx[k] = (t.floor().max(0.0) as usize).min(self.n[k] - 1);
        }
        let bin = (idx[0] * self.n[1] + idx[1]) * self.n[2] + idx[2];
        &self.items[self.start[bin] as usize..self.start[bin + 1] as usize]
    }
}

/// Barycentric coordinates of `p` in cell `c`.
fn barycentric(mesh: &SimplexMesh, c: usize, p: &[f64]) -> Option<[f64; 4]> {
    let dim = mesh.dim;
    let v = mesh.cell(c);
    let p0 = mesh.point(v[0] as usize);
    let mut t = [[0.0f64; 3]; 3];
    for j in 0..dim {
        let pj = mesh.point(v[j + 1] as usize);
        for k in 0..dim {
            t[k][j] = pj[k] - p0[k];
        }
    }
    let mut d = [0.0f64; 3];
    for k in 0..dim {
        d[k] = p[k] - p0[k];
    }
    let l = if dim == 2 {
        let det = t[0][0] * t[1][1] - t[0][1] * t[1][0];
        if det == 0.0 {
            return None;
        }
        [(t[1][1] * d[0] - t[0][1] * d[1]) / det, (t[0][0] * d[1] - t[1][0] * d[0]) / det, 0.0]
    } else {
        let det3 = |m: &[[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let det = det3(&t);
        if det == 0.0 {
            return None;
        }
        let mut out = [0.0; 3];
        for j in 0..3 {
            let mut m = t;
            for k in 0..3 {
                m[k][j] = d[k];
            }
            out[j] = det3(&m) / det;
        }
        out
    };
    let mut lam = [0.0; 4];
    lam[0] = 1.0 - l[..dim].iter().sum::<f64>();
    lam[1..=dim].copy_from_slice(&l[..dim]);
    Some(lam)
}

struct Located {
    cell: usize,
    lam: Vec<f64>,
}

fn locate(mesh: &SimplexMesh, loc: &Locator, p: &[f64]) -> Option<Located> {
    for &c in loc.candidates(p) {
        if let Some(l) = barycentric(mesh, c as usize, p) {
            if l[..=mesh.dim].iter().all(|x| *x >= -INSIDE) {
                return Some(Located { cell: c as usize, lam: l[..=mesh.dim].to_vec() });
            }
        }
    }
    None
}

// ------------------------------------------------------------------------------------------------------------- the commands

fn number_array(v: &Value, key: &str, dim: usize) -> Result<Vec<f64>, Error> {
    let a = v.as_array().ok_or_else(|| Error::invalid(format!("`{key}` is a list of {dim} numbers")))?;
    if a.len() != dim {
        return invalid(format!("`{key}` has {} numbers; this is a {dim}D mesh (x then z; 3D: x, y, z)", a.len()));
    }
    a.iter()
        .map(|x| {
            x.as_f64()
                .filter(|f| f.is_finite())
                .ok_or_else(|| Error::invalid(format!("`{key}` holds a number that is not finite")))
        })
        .collect()
}

fn pair(re: f64, im: f64) -> Value {
    json!([re, im])
}

impl Engine {
    fn revolution(args: &Value) -> bool {
        args.get("revolution").and_then(Value::as_bool).unwrap_or(false)
    }

    fn region_arg(args: &Value) -> Option<&str> {
        args.get("region").and_then(Value::as_str)
    }

    /// A note when a result of a study named `axisymmetric...` is integrated without the revolution weight.
    fn axisym_note(&self, p: &Prep, revolution: bool, out: &mut Value) {
        if !revolution
            && p.physics.iter().any(|k| self.results.get(k).is_some_and(|s| s.result.study.contains("axisym")))
        {
            out["note"] = json!(
                "this is an axisymmetric result: pass `revolution: true` for the 2 pi r weight of a body of revolution"
            );
        }
    }

    pub(crate) fn cmd_integrate(&self, args: &Value) -> R {
        let (tree, p) = self.one_expr(args, "expr")?;
        let region = p.region(Self::region_arg(args))?;
        let rev = Self::revolution(args);
        let s = p.integrate(&tree, &region, rev)?;
        let real = tree.is_real_valued(&p.complex_flags());
        let mut out = json!({
            "value": if real { json!(s.re) } else { pair(s.re, s.im) }, "complex": !real,
            "measure": s.measure, "cells": region.cells, "revolution": rev, "region": Self::region_arg(args),
        });
        self.axisym_note(&p, rev, &mut out);
        self.decorate(&mut out, &p);
        Ok(out)
    }

    pub(crate) fn cmd_extremum(&self, args: &Value) -> R {
        let (tree, p) = self.one_expr(args, "expr")?;
        if !tree.is_real_valued(&p.complex_flags()) {
            return invalid("the expression is complex and has no order: wrap it in abs, abs2, re or im");
        }
        let op = args.get("op").and_then(Value::as_str).unwrap_or("max");
        if op != "max" && op != "min" {
            return invalid("`op` is \"max\" or \"min\"");
        }
        let region = p.region(Self::region_arg(args))?;
        let (v, c, a) = p.extremum(&tree, &region, op == "max")?;
        let mut out = json!({
            "op": op, "value": v, "cell": c, "node": a, "point": p.node_point(c, a), "region": Self::region_arg(args),
            "at": "the nodes of the cells (not between them)",
        });
        self.decorate(&mut out, &p);
        Ok(out)
    }

    fn sample(&self, p: &Prep, tree: &Tree, points: &[Vec<f64>]) -> Value {
        let loc = Locator::new(&p.mesh);
        let nloc = p.nloc;
        let flags = p.complex_flags();
        let real = tree.is_real_valued(&flags);
        let (mut re, mut im, mut inside, mut cells) = (vec![], vec![], vec![], vec![]);
        let mut local = vec![C::new(0.0, 0.0); p.leaves.len() * nloc];
        let mut dofs = vec![];
        let mut work = Vec::with_capacity(p.leaves.len());
        for pt in points {
            match locate(&p.mesh, &loc, pt) {
                Some(l) => {
                    p.gather(l.cell, &mut dofs, &mut local);
                    let (n, _) = shape(&l.lam, p.order);
                    let v = p.at_point(tree, &local, &n, &mut work);
                    re.push(json!(v.re));
                    im.push(json!(v.im));
                    inside.push(true);
                    cells.push(json!(l.cell));
                }
                None => {
                    re.push(Value::Null);
                    im.push(Value::Null);
                    inside.push(false);
                    cells.push(Value::Null);
                }
            }
        }
        let mut out = json!({ "values": re, "inside": inside, "cells": cells, "complex": !real });
        if !real {
            out["im"] = json!(im);
        }
        out
    }

    pub(crate) fn cmd_probe(&self, args: &Value) -> R {
        let (tree, p) = self.one_expr(args, "expr")?;
        let dim = p.mesh.dim;
        let pts = args
            .get("points")
            .and_then(Value::as_array)
            .ok_or_else(|| Error::invalid("`points` is required: a list of points"))?;
        if pts.is_empty() || pts.len() > MAX_POINTS {
            return invalid(format!("`points` holds 1 to {MAX_POINTS} points"));
        }
        let points: Vec<Vec<f64>> = pts.iter().map(|x| number_array(x, "points", dim)).collect::<Result<_, _>>()?;
        let mut out = self.sample(&p, &tree, &points);
        out["points"] = json!(points);
        self.decorate(&mut out, &p);
        Ok(out)
    }

    pub(crate) fn cmd_line_cut(&self, args: &Value) -> R {
        let (tree, p) = self.one_expr(args, "expr")?;
        let dim = p.mesh.dim;
        let from = number_array(args.get("from").unwrap_or(&Value::Null), "from", dim)?;
        let to = number_array(args.get("to").unwrap_or(&Value::Null), "to", dim)?;
        let n = args.get("n").and_then(Value::as_u64).unwrap_or(101) as usize;
        if !(2..=MAX_POINTS).contains(&n) {
            return invalid(format!("`n`, the number of points, is 2 to {MAX_POINTS}"));
        }
        let length = (0..dim).map(|k| (to[k] - from[k]) * (to[k] - from[k])).sum::<f64>().sqrt();
        let mut s = Vec::with_capacity(n);
        let points: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                s.push(length * t);
                (0..dim).map(|k| from[k] + (to[k] - from[k]) * t).collect()
            })
            .collect();
        let mut out = self.sample(&p, &tree, &points);
        out["s"] = json!(s);
        out["from"] = json!(from);
        out["to"] = json!(to);
        self.decorate(&mut out, &p);
        Ok(out)
    }

    /// The expressions of a list argument (`a`, `b`, `fields`): one expression, or a list of components.
    fn components<'a>(args: &'a Value, key: &str) -> Result<Vec<&'a Value>, Error> {
        match args.get(key) {
            None | Some(Value::Null) => invalid(format!("`{key}` is required: a field, or a list of components")),
            Some(Value::Array(a)) if a.is_empty() => invalid(format!("`{key}` is empty")),
            Some(Value::Array(a)) => Ok(a.iter().collect()),
            Some(v) => Ok(vec![v]),
        }
    }

    pub(crate) fn cmd_overlap(&self, args: &Value) -> R {
        let (a, b) = (Self::components(args, "a")?, Self::components(args, "b")?);
        if a.len() != b.len() {
            return invalid(format!("`a` has {} components and `b` {}; they must match", a.len(), b.len()));
        }
        let all: Vec<&Value> = a.iter().chain(b.iter()).copied().collect();
        let (trees, p) = self.expressions(args, &all)?;
        let (ta, tb) = trees.split_at(a.len());
        let sum = |mut it: Vec<Tree>| -> Tree {
            let first = it.remove(0);
            it.into_iter().fold(first, |acc, t| Tree::Bin(Op::Add, Box::new(acc), Box::new(t)))
        };
        let cross = sum(ta
            .iter()
            .zip(tb)
            .map(|(x, y)| Tree::Bin(Op::Mul, Box::new(Tree::Call(Fun::Conj, Box::new(x.clone()))), Box::new(y.clone())))
            .collect());
        let na = sum(ta.iter().map(|x| Tree::Call(Fun::Abs2, Box::new(x.clone()))).collect());
        let nb = sum(tb.iter().map(|x| Tree::Call(Fun::Abs2, Box::new(x.clone()))).collect());
        let region = p.region(Self::region_arg(args))?;
        let rev = Self::revolution(args);
        let x = p.integrate(&cross, &region, rev)?;
        let (xa, xb) = (p.integrate(&na, &region, rev)?, p.integrate(&nb, &region, rev)?);
        let coefficient = (x.re * x.re + x.im * x.im) / (xa.re * xb.re);
        let mut out = json!({
            "overlap": pair(x.re, x.im), "norm_a": xa.re, "norm_b": xb.re, "coefficient": coefficient,
            "cells": region.cells, "revolution": rev, "region": Self::region_arg(args),
            "meaning": "overlap = the integral of sum_k conj(a_k) b_k; norm = the integral of sum_k |c_k|^2; coefficient = |overlap|^2 / (norm_a norm_b), 0 to 1",
        });
        self.axisym_note(&p, rev, &mut out);
        self.decorate(&mut out, &p);
        Ok(out)
    }

    pub(crate) fn cmd_confinement(&self, args: &Value) -> R {
        let key = if args.get("fields").is_some() { "fields" } else { "expr" };
        let mut comps = Self::components(args, key)?;
        // an optional real weight (a formula: `eps_r` for the energy fraction) multiplies the integrand at every point; without one the tree is as it always was
        let weight_arg = args.get("weight").filter(|v| !v.is_null());
        if let Some(w) = weight_arg {
            comps.push(w);
        }
        let (mut trees, p) = self.expressions(args, &comps)?;
        let weight = weight_arg.map(|_| trees.pop().expect("the weight"));
        if weight.as_ref().is_some_and(|w| !w.is_real_valued(&p.complex_flags())) {
            return invalid("`weight` must be real: wrap a complex formula in abs, abs2, re or im");
        }
        let mut it = trees.into_iter().map(|t| Tree::Call(Fun::Abs2, Box::new(t)));
        let first = it.next().expect("a component");
        let mut power = it.fold(first, |acc, t| Tree::Bin(Op::Add, Box::new(acc), Box::new(t)));
        if let Some(w) = weight {
            power = Tree::Bin(Op::Mul, Box::new(w), Box::new(power));
        }
        let rev = Self::revolution(args);
        let all = p.region(None)?;
        let part = p.region(Some(
            Self::region_arg(args)
                .ok_or_else(|| Error::invalid("`region` is required: the region the fraction is of"))?,
        ))?;
        let (inside, total) = (p.integrate(&power, &part, rev)?, p.integrate(&power, &all, rev)?);
        let mut out = json!({
            "fraction": inside.re / total.re, "inside": inside.re, "total": total.re, "cells": part.cells,
            "revolution": rev, "region": Self::region_arg(args),
            "weighted": weight_arg.is_some(),
            "meaning": if weight_arg.is_some() {
                "the integral of weight * sum_k |c_k|^2 over the region divided by the same over the whole mesh"
            } else {
                "the integral of sum_k |c_k|^2 over the region divided by the same over the whole mesh (no material weight)"
            },
        });
        self.axisym_note(&p, rev, &mut out);
        self.decorate(&mut out, &p);
        Ok(out)
    }
}
