// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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
//! **Boundary integrals (N2b).** `integrate` with `boundary` sums over the faces of a named boundary instead of the cells: the facet rule of
//! `fem::facet_rule` (3-point Gauss on an edge, the degree-4 triangle rule on a face) at points mapped into the barycentric coordinates of the
//! cell each face is read in, `measure * sum_q (w_q f_q)` per face, in blocks of `BLOCK` faces in face order. A face next to two cells (an
//! interface) is read in the cell of the region `side`; the normal `nx ny nz` points out of that cell. See ADR 90, N2b.
//!
//! **Maxima** are taken over the values at the nodes of the cells (the vertices, and the edge midpoints of a P2 mesh), not over the interpolant:
//! the true maximum of a P2 field between nodes can be larger. A node shared by cells has one value per cell (the fields are not averaged);
//! all of them count. **Points** (`probe`, `line_cut`) are in mesh axes (a 2D mesh: x then z of the wafer; 3D: x, y, z); a point outside the
//! mesh gives `inside: false` and a null value, not an error; a point on the face between cells is read in the cell with the lowest number.

use crate::engine::Engine;
use femcore::fem::{facet_rule, facets_with, geometry, quad_rule, shape, RefElement};
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
    /// a component of the unit normal of a boundary facet (wafer axis 0 x, 1 y, 2 z); only in an `integrate` over a `boundary`
    Normal(usize),
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
#[derive(Clone)]
pub(crate) struct Defaults {
    physics: Option<String>,
    solution: usize,
    /// `wavelength` of the request: where the materials are read for a dispersive model (m); None is the wavelength of the study
    wavelength: Option<f64>,
    /// the request integrates over a `boundary`, so `nx ny nz` mean something
    boundary: bool,
    /// the request evaluates over a batch of `solutions`, so bare field handles follow `solution`
    batch_solutions: bool,
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
    /// a component of the normal of the facet being integrated (set per facet by the boundary integral, 0 elsewhere)
    Normal(usize),
}

struct Leaf {
    res: Arc<SolveResult>,
    sol: usize,
    kind: Kind,
    complex: bool,
    unit: String,
}

/// How the side of an interior boundary or sheet is chosen in `integrate`.
enum SideSpec {
    /// Domain label(s) of the region on one side of the interface/sheet.
    Domains(Vec<i32>),
    /// Signed wafer axis `(axis: 0=x, 1=y, 2=z, positive: true for +axis side, false for -axis side)`
    /// for a sheet embedded inside a single homogeneous region (`rb == ra`).
    Axis(usize, bool),
}

/// A face of a boundary: its vertices (the first `dim` used) and the cells next to it as (cell, local number of its vertex opposite the face).
struct Face {
    v: [u32; 3],
    cells: Vec<(usize, usize)>,
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

#[inline]
fn cell_vol(mesh: &SimplexMesh, c: usize) -> f64 {
    let v = mesh.cell(c);
    let x0 = mesh.point(v[0] as usize);
    if mesh.dim == 2 {
        let x1 = mesh.point(v[1] as usize);
        let x2 = mesh.point(v[2] as usize);
        let det = (x1[0] - x0[0]) * (x2[1] - x0[1]) - (x1[1] - x0[1]) * (x2[0] - x0[0]);
        det.abs() / 2.0
    } else {
        let x1 = mesh.point(v[1] as usize);
        let x2 = mesh.point(v[2] as usize);
        let x3 = mesh.point(v[3] as usize);
        let m0 = [x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]];
        let m1 = [x2[0] - x0[0], x2[1] - x0[1], x2[2] - x0[2]];
        let m2 = [x3[0] - x0[0], x3[1] - x0[1], x3[2] - x0[2]];
        let c00 = m1[1] * m2[2] - m1[2] * m2[1];
        let c01 = m1[2] * m2[0] - m1[0] * m2[2];
        let c02 = m1[0] * m2[1] - m1[1] * m2[0];
        let det = m0[0] * c00 + m0[1] * c01 + m0[2] * c02;
        det.abs() / 6.0
    }
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
        let sol = if matches!(what, What::Field) { solution } else { 0 };
        Tree::Leaf(slot(self.specs, Spec { physics: physics.to_string(), solution: sol, name: name.to_string(), what }))
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
        if let Some(axis) = ["nx", "ny", "nz"].iter().position(|a| *a == name) {
            if !self.d.boundary {
                return Err(format!(
                    "'{name}' is a component of the normal of a boundary: it exists in `integrate` with `boundary`"
                ));
            }
            return Ok(Some(Named::Scalar(self.slot(&physics, solution, &name, What::Normal(axis)))));
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

/// Whether a JSON object is a single expression AST node (`op`, `field`, `num`, `const`, `coord`, `material`, `region`, `normal`)
/// rather than a batch map of named expressions (`{"energy": "...", "ex2": "..."}`).
fn is_ast_object(m: &serde_json::Map<String, Value>) -> bool {
    if m.is_empty() {
        return false;
    }
    if m.contains_key("op") && m.keys().all(|k| matches!(k.as_str(), "op" | "args" | "n")) {
        return true;
    }
    if m.contains_key("field") && m.keys().all(|k| matches!(k.as_str(), "field" | "physics" | "solution")) {
        return true;
    }
    for key in ["num", "const"] {
        if let Some(v) = m.get(key) {
            if m.len() == 1
                && (v.is_number() || v.as_array().is_some_and(|a| a.len() == 2 && a[0].is_number() && a[1].is_number()))
            {
                return true;
            }
        }
    }
    let plain_name = |v: &Value| {
        v.as_str().is_some_and(|s| {
            !s.trim().is_empty()
                && !s.chars().any(|c| "+-*/^()<>=!,[]\"'`.:".contains(c))
                && !s.contains(" and ")
                && !s.contains(" or ")
        })
    };
    if let Some(v) = m.get("coord") {
        if m.keys().all(|k| matches!(k.as_str(), "coord" | "physics" | "solution"))
            && (v.as_str().is_some_and(|s| matches!(s, "x" | "y" | "z")) || v.as_u64().is_some_and(|a| a < 3))
        {
            return true;
        }
    }
    if let Some(v) = m.get("normal") {
        if m.keys().all(|k| matches!(k.as_str(), "normal" | "physics" | "solution"))
            && (v.as_str().is_some_and(|s| matches!(s, "nx" | "ny" | "nz" | "x" | "y" | "z"))
                || v.as_u64().is_some_and(|a| a < 3))
        {
            return true;
        }
    }
    if let Some(v) = m.get("material") {
        if m.keys().all(|k| matches!(k.as_str(), "material" | "physics" | "solution" | "wavelength")) && plain_name(v) {
            return true;
        }
    }
    if let Some(v) = m.get("region") {
        if m.keys().all(|k| matches!(k.as_str(), "region" | "physics" | "solution")) && plain_name(v) {
            return true;
        }
    }
    false
}

/// The shape of `expr` in `integrate`: a single expression, a batch list, or a batch map of named expressions.
enum BatchExpr<'a> {
    Single(&'a Value),
    List(Vec<&'a Value>),
    Map(Vec<(String, &'a Value)>),
}

impl<'a> BatchExpr<'a> {
    fn from_value(v: &'a Value) -> Result<Self, Error> {
        match v {
            Value::Array(arr) => {
                if arr.is_empty() {
                    return invalid("`expr` list cannot be empty");
                }
                Ok(BatchExpr::List(arr.iter().collect()))
            }
            Value::Object(m) if !is_ast_object(m) => {
                if m.is_empty() {
                    return invalid("`expr` map cannot be empty");
                }
                Ok(BatchExpr::Map(m.iter().map(|(k, val)| (k.clone(), val)).collect()))
            }
            other => Ok(BatchExpr::Single(other)),
        }
    }

    fn is_multi(&self) -> bool {
        !matches!(self, BatchExpr::Single(_))
    }

    fn values(&self) -> Vec<&'a Value> {
        match self {
            BatchExpr::Single(v) => vec![*v],
            BatchExpr::List(vs) => vs.clone(),
            BatchExpr::Map(entries) => entries.iter().map(|(_, v)| *v).collect(),
        }
    }
}

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
        let boundary = args.get("boundary").is_some_and(|v| !v.is_null());
        Ok(Defaults { physics, solution, wavelength, boundary, batch_solutions: false })
    }

    /// Resolves `solutions` (`"all"` or `[0, 1, ...]`) or falls back to `vec![d.solution]`.
    /// Refuses if both `solution` and `solutions` are explicitly passed, or if `solutions` is empty / out-of-range.
    fn parse_solutions(&self, args: &Value, d: &Defaults) -> Result<(Vec<usize>, bool), Error> {
        let Some(sols_val) = args.get("solutions").filter(|v| !v.is_null()) else {
            return Ok((vec![d.solution], false));
        };
        if args.get("solution").is_some_and(|v| !v.is_null()) {
            return invalid("give `solution` or `solutions`, not both");
        }
        let physics = d.physics.as_deref().ok_or_else(|| {
            Error::invalid(format!(
                "give `physics` with `solutions` (the results are: {})",
                self.results.keys().cloned().collect::<Vec<_>>().join(", ")
            ))
        })?;
        let stored = self.results.get(physics).ok_or_else(|| {
            Error::invalid(format!(
                "there is no result for '{physics}' yet; solve first{}",
                if self.results.is_empty() {
                    String::new()
                } else {
                    format!(" (the results are: {})", self.results.keys().cloned().collect::<Vec<_>>().join(", "))
                }
            ))
        })?;
        let n_sol = stored.result.solutions.len();
        match sols_val {
            Value::String(s) if s == "all" => {
                if n_sol == 0 {
                    return invalid(format!("the result of '{physics}' has no solutions"));
                }
                Ok(((0..n_sol).collect(), true))
            }
            Value::Array(arr) => {
                if arr.is_empty() {
                    return invalid("`solutions` cannot be empty");
                }
                let mut out = Vec::with_capacity(arr.len());
                for v in arr {
                    let idx = v.as_u64().ok_or_else(|| {
                        Error::invalid("`solutions` must be \"all\" or a list of whole solution numbers")
                    })? as usize;
                    if idx >= n_sol {
                        return invalid(format!("the result of '{physics}' has {n_sol} solutions; there is no {idx}"));
                    }
                    out.push(idx);
                }
                Ok((out, true))
            }
            _ => invalid("`solutions` must be \"all\" or a list of whole solution numbers"),
        }
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
                    let solution = if d.batch_solutions {
                        d.solution
                    } else {
                        match m.get("solution") {
                            Some(s) => {
                                s.as_u64().ok_or_else(|| Error::invalid("`solution` must be a whole number"))? as usize
                            }
                            None => d.solution,
                        }
                    };
                    return leaf(name, physics, solution, specs);
                }
                if let Some(n) = m.get("num").or_else(|| m.get("const")) {
                    return match n {
                        Value::Number(x) => Ok(Tree::num(x.as_f64().unwrap_or(f64::NAN))),
                        Value::Array(a) if a.len() == 2 => match (a[0].as_f64(), a[1].as_f64()) {
                            (Some(re), Some(im)) => Ok(Tree::Num(C::new(re, im))),
                            _ => invalid("`num` pairs are [re, im] numbers"),
                        },
                        _ => invalid("`num` is a number or a pair [re, im]"),
                    };
                }
                if let Some(c) = m.get("coord") {
                    let axis = match c {
                        Value::String(s) => ["x", "y", "z"]
                            .iter()
                            .position(|a| a == s)
                            .ok_or_else(|| Error::invalid("`coord` must be \"x\", \"y\" or \"z\""))?,
                        Value::Number(n) => n
                            .as_u64()
                            .filter(|&a| a < 3)
                            .ok_or_else(|| Error::invalid("`coord` axis must be 0, 1 or 2"))?
                            as usize,
                        _ => return invalid("`coord` must be \"x\", \"y\" or \"z\""),
                    };
                    let physics = match m.get("physics").and_then(Value::as_str) {
                        Some(p) => p,
                        None => d.physics.as_deref().ok_or_else(|| Error::invalid("give `physics`"))?,
                    };
                    return Ok(Tree::Leaf(slot(
                        specs,
                        Spec {
                            physics: physics.to_string(),
                            solution: 0,
                            name: AXES[axis].to_string(),
                            what: What::Coord(axis),
                        },
                    )));
                }
                if let Some(n) = m.get("normal") {
                    if !d.boundary {
                        return invalid(
                            "`nx`, `ny` and `nz` are the normal of a boundary: they exist in `integrate` with `boundary`",
                        );
                    }
                    let axis = match n {
                        Value::String(s) => ["nx", "ny", "nz"]
                            .iter()
                            .position(|a| a == s)
                            .or_else(|| ["x", "y", "z"].iter().position(|a| a == s))
                            .ok_or_else(|| Error::invalid("`normal` must be \"nx\", \"ny\" or \"nz\""))?,
                        Value::Number(n) => n
                            .as_u64()
                            .filter(|&a| a < 3)
                            .ok_or_else(|| Error::invalid("`normal` axis must be 0, 1 or 2"))?
                            as usize,
                        _ => return invalid("`normal` must be \"nx\", \"ny\" or \"nz\""),
                    };
                    let physics = match m.get("physics").and_then(Value::as_str) {
                        Some(p) => p,
                        None => d.physics.as_deref().ok_or_else(|| Error::invalid("give `physics`"))?,
                    };
                    return Ok(Tree::Leaf(slot(
                        specs,
                        Spec {
                            physics: physics.to_string(),
                            solution: 0,
                            name: format!("n{}", AXES[axis]),
                            what: What::Normal(axis),
                        },
                    )));
                }
                if let Some(mat) = m.get("material").and_then(Value::as_str) {
                    let mut namer = Namer { engine: self, d, specs, missing: None, unknown: false };
                    return namer.material(mat).map_err(|e| namer.missing.take().unwrap_or_else(|| Error::invalid(e)));
                }
                if let Some(reg) = m.get("region").and_then(Value::as_str) {
                    let mut namer = Namer { engine: self, d, specs, missing: None, unknown: false };
                    return namer.region(reg).map_err(|e| namer.missing.take().unwrap_or_else(|| Error::invalid(e)));
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
                What::Normal(a) => (Kind::Normal(*a), false, String::new()),
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
        self.expressions_for_solutions(&d, &[d.solution], false, exprs)
    }

    /// Parses `exprs` across `solutions` (`sol_idx` outer, `expr_idx` inner) against one deduplicated `Prep`.
    pub(crate) fn expressions_for_solutions(
        &self,
        d: &Defaults,
        solutions: &[usize],
        batch_solutions: bool,
        exprs: &[&Value],
    ) -> Result<(Vec<Tree>, Prep), Error> {
        fn can_remap(v: &Value, batch_solutions: bool) -> bool {
            match v {
                Value::Number(_) => true,
                Value::String(s) => !s.contains(':'),
                Value::Object(m) => {
                    if m.contains_key("field") {
                        return batch_solutions || !m.contains_key("solution");
                    }
                    if let Some(args) = m.get("args").and_then(Value::as_array) {
                        return args.iter().all(|a| can_remap(a, batch_solutions));
                    }
                    true
                }
                _ => false,
            }
        }
        fn remap_tree(t: &Tree, sol: usize, specs: &mut Vec<Spec>) -> Tree {
            match t {
                Tree::Num(z) => Tree::Num(*z),
                Tree::Leaf(i) => {
                    let mut s = specs[*i].clone();
                    if matches!(s.what, What::Field) {
                        s.solution = sol;
                    }
                    Tree::Leaf(slot(specs, s))
                }
                Tree::Neg(a) => Tree::Neg(Box::new(remap_tree(a, sol, specs))),
                Tree::Bin(op, a, b) => {
                    Tree::Bin(*op, Box::new(remap_tree(a, sol, specs)), Box::new(remap_tree(b, sol, specs)))
                }
                Tree::Call(f, a) => Tree::Call(*f, Box::new(remap_tree(a, sol, specs))),
                Tree::Call2(f, a, b) => {
                    Tree::Call2(*f, Box::new(remap_tree(a, sol, specs)), Box::new(remap_tree(b, sol, specs)))
                }
                Tree::Pow(a, n) => Tree::Pow(Box::new(remap_tree(a, sol, specs)), *n),
                Tree::Cmp(c, a, b) => {
                    Tree::Cmp(*c, Box::new(remap_tree(a, sol, specs)), Box::new(remap_tree(b, sol, specs)))
                }
                Tree::Where(c, a, b) => Tree::Where(
                    Box::new(remap_tree(c, sol, specs)),
                    Box::new(remap_tree(a, sol, specs)),
                    Box::new(remap_tree(b, sol, specs)),
                ),
            }
        }

        let mut specs = vec![];
        let mut trees = Vec::with_capacity(solutions.len() * exprs.len());
        let remappable = solutions.len() > 1 && exprs.iter().all(|e| can_remap(e, batch_solutions));
        if remappable {
            let sol0 = solutions[0];
            let mut d_sol = d.clone();
            d_sol.solution = sol0;
            d_sol.batch_solutions = batch_solutions;
            for e in exprs {
                trees.push(self.parse_expr(e, &d_sol, &mut specs, 0)?);
            }
            let n_exprs = exprs.len();
            for &sol in &solutions[1..] {
                for e_idx in 0..n_exprs {
                    let remapped = remap_tree(&trees[e_idx], sol, &mut specs);
                    trees.push(remapped);
                }
            }
        } else {
            for &sol in solutions {
                let mut d_sol = d.clone();
                d_sol.solution = sol;
                d_sol.batch_solutions = batch_solutions;
                for e in exprs {
                    trees.push(self.parse_expr(e, &d_sol, &mut specs, 0)?);
                }
            }
        }
        // an expression of numbers only (`1` for an area) still needs a mesh to integrate over: the one of the default result
        if specs.is_empty() {
            if let Some(p) = &d.physics {
                let sol = solutions.first().copied().unwrap_or(d.solution);
                specs.push(Spec { physics: p.clone(), solution: sol, name: "x".into(), what: What::Coord(0) });
            }
        }
        let prep = self.prepare(&specs)?;
        let flags = prep.complex_flags();
        let check_count = if remappable { exprs.len() } else { trees.len() };
        let n_exprs = exprs.len();
        for (idx, t) in trees.iter().take(check_count).enumerate() {
            let e = exprs[idx % n_exprs];
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

/// Per-cell domain index map for `by_region: true` single-pass reductions.
struct DomainMap {
    names: Vec<String>,
    cells: Vec<usize>,
    cell_dom: Vec<usize>,
}

/// Output of a single-pass batch integral over `trees`.
struct BatchSum {
    total: Vec<(f64, f64)>,
    measure: f64,
    region_totals: Vec<Vec<(f64, f64)>>,
    region_measures: Vec<f64>,
}

/// Output of a single-pass batch confinement integral over `power_trees` and `regions`.
struct ConfinementBatch {
    total: Vec<f64>,
    inside: Vec<Vec<f64>>,
}

/// Pre-borrowed leaf view avoiding per-cell `loc` buffer copies for `Cell` and `Const` leaves.
enum FastLeaf<'a> {
    CellReal(&'a [f64]),
    CellComplex(&'a [f64], &'a [f64]),
    Const(&'a [f64]),
    Slow,
}

#[inline(always)]
fn dot_n(n: &[f64], s: &[f64]) -> f64 {
    if n.len() == 6 && s.len() == 6 {
        let mut r = 0.0;
        r += n[0] * s[0];
        r += n[1] * s[1];
        r += n[2] * s[2];
        r += n[3] * s[3];
        r += n[4] * s[4];
        r += n[5] * s[5];
        r
    } else if n.len() == 10 && s.len() == 10 {
        let mut r = 0.0;
        for a in 0..10 {
            r += n[a] * s[a];
        }
        r
    } else {
        let mut r = 0.0;
        for a in 0..n.len() {
            r += n[a] * s[a];
        }
        r
    }
}

#[inline(always)]
fn dot_n_complex(n: &[f64], sr: &[f64], si: &[f64]) -> C {
    if n.len() == 6 && sr.len() == 6 && si.len() == 6 {
        let (mut r, mut im_v) = (0.0, 0.0);
        r += n[0] * sr[0];
        im_v += n[0] * si[0];
        r += n[1] * sr[1];
        im_v += n[1] * si[1];
        r += n[2] * sr[2];
        im_v += n[2] * si[2];
        r += n[3] * sr[3];
        im_v += n[3] * si[3];
        r += n[4] * sr[4];
        im_v += n[4] * si[4];
        r += n[5] * sr[5];
        im_v += n[5] * si[5];
        C::new(r, im_v)
    } else if n.len() == 10 && sr.len() == 10 && si.len() == 10 {
        let (mut r, mut im_v) = (0.0, 0.0);
        for a in 0..10 {
            r += n[a] * sr[a];
            im_v += n[a] * si[a];
        }
        C::new(r, im_v)
    } else {
        let (mut r, mut im_v) = (0.0, 0.0);
        for a in 0..n.len() {
            r += n[a] * sr[a];
            im_v += n[a] * si[a];
        }
        C::new(r, im_v)
    }
}

/// Flattened, common-subexpression-deduplicated instruction step for batch tree evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
enum FlatStep {
    Num(C),
    Add(usize, usize),
    Sub(usize, usize),
    Mul(usize, usize),
    Div(usize, usize),
    Neg(usize),
    Abs(usize),
    Abs2(usize),
    Re(usize),
    Im(usize),
    Conj(usize),
    Pow(usize, u32),
}

struct FlatDag {
    n_leaves: usize,
    steps: Vec<FlatStep>,
    roots: Vec<usize>,
}

impl FlatDag {
    fn compile(trees: &[Tree], n_leaves: usize) -> Option<Self> {
        fn add_node(t: &Tree, n_leaves: usize, steps: &mut Vec<FlatStep>) -> Option<usize> {
            let step = match t {
                Tree::Leaf(i) => return Some(*i),
                Tree::Num(z) => FlatStep::Num(*z),
                Tree::Neg(a) => FlatStep::Neg(add_node(a, n_leaves, steps)?),
                Tree::Bin(Op::Add, a, b) => FlatStep::Add(add_node(a, n_leaves, steps)?, add_node(b, n_leaves, steps)?),
                Tree::Bin(Op::Sub, a, b) => FlatStep::Sub(add_node(a, n_leaves, steps)?, add_node(b, n_leaves, steps)?),
                Tree::Bin(Op::Mul, a, b) => FlatStep::Mul(add_node(a, n_leaves, steps)?, add_node(b, n_leaves, steps)?),
                Tree::Bin(Op::Div, a, b) => FlatStep::Div(add_node(a, n_leaves, steps)?, add_node(b, n_leaves, steps)?),
                Tree::Call(Fun::Abs, a) => FlatStep::Abs(add_node(a, n_leaves, steps)?),
                Tree::Call(Fun::Abs2, a) => FlatStep::Abs2(add_node(a, n_leaves, steps)?),
                Tree::Call(Fun::Re, a) => FlatStep::Re(add_node(a, n_leaves, steps)?),
                Tree::Call(Fun::Im, a) => FlatStep::Im(add_node(a, n_leaves, steps)?),
                Tree::Call(Fun::Conj, a) => FlatStep::Conj(add_node(a, n_leaves, steps)?),
                Tree::Pow(a, n) => FlatStep::Pow(add_node(a, n_leaves, steps)?, *n),
                _ => return None,
            };
            if let Some(pos) = steps.iter().position(|s| *s == step) {
                Some(n_leaves + pos)
            } else {
                steps.push(step);
                Some(n_leaves + steps.len() - 1)
            }
        }

        let mut steps = Vec::new();
        let mut roots = Vec::with_capacity(trees.len());
        for t in trees {
            roots.push(add_node(t, n_leaves, &mut steps)?);
        }
        Some(FlatDag { n_leaves, steps, roots })
    }

    #[inline]
    fn reg_len(&self) -> usize {
        self.n_leaves + self.steps.len()
    }

    #[inline(always)]
    fn eval_all(&self, reg: &mut [C]) {
        for (k, &step) in self.steps.iter().enumerate() {
            let v = match step {
                FlatStep::Num(z) => z,
                FlatStep::Neg(a) => {
                    let z = reg[a];
                    C::new(-z.re, -z.im)
                }
                FlatStep::Add(a, b) => {
                    let (x, y) = (reg[a], reg[b]);
                    C::new(x.re + y.re, x.im + y.im)
                }
                FlatStep::Sub(a, b) => {
                    let (x, y) = (reg[a], reg[b]);
                    C::new(x.re - y.re, x.im - y.im)
                }
                FlatStep::Mul(a, b) => {
                    let (x, y) = (reg[a], reg[b]);
                    C::new(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re)
                }
                FlatStep::Div(a, b) => {
                    let (x, y) = (reg[a], reg[b]);
                    let d = y.re * y.re + y.im * y.im;
                    C::new((x.re * y.re + x.im * y.im) / d, (x.im * y.re - x.re * y.im) / d)
                }
                FlatStep::Abs(a) => {
                    let z = reg[a];
                    C::new((z.re * z.re + z.im * z.im).sqrt(), 0.0)
                }
                FlatStep::Abs2(a) => {
                    let z = reg[a];
                    C::new(z.re * z.re + z.im * z.im, 0.0)
                }
                FlatStep::Re(a) => C::new(reg[a].re, 0.0),
                FlatStep::Im(a) => C::new(reg[a].im, 0.0),
                FlatStep::Conj(a) => {
                    let z = reg[a];
                    C::new(z.re, -z.im)
                }
                FlatStep::Pow(a, n) => {
                    let z = reg[a];
                    let mut r = z;
                    for _ in 1..n {
                        r = C::new(r.re * z.re - r.im * z.im, r.re * z.im + r.im * z.re);
                    }
                    r
                }
            };
            reg[self.n_leaves + k] = v;
        }
    }
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

    /// Maps each cell in `region` to its domain index in `self.model.regions()`, optionally filtered by `region_filter`.
    fn build_domain_map(&self, region: &Region, region_filter: Option<&str>) -> Result<DomainMap, Error> {
        use std::collections::HashMap;
        let allowed = match region_filter.filter(|r| !r.trim().is_empty()) {
            Some(text) => Some(self.model.domains(text)?),
            None => None,
        };
        let mut names: Vec<String> = Vec::new();
        let mut cells: Vec<usize> = Vec::new();
        let mut id_to_dom: HashMap<i32, usize> = HashMap::new();
        for (id, name) in self.model.regions() {
            if allowed.as_ref().is_some_and(|a| !a.contains(&id)) {
                continue;
            }
            let d_idx = match names.iter().position(|n| *n == name) {
                Some(i) => i,
                None => {
                    names.push(name);
                    cells.push(0);
                    names.len() - 1
                }
            };
            id_to_dom.insert(id, d_idx);
        }
        let n_cells = self.mesh.n_cells();
        let mut cell_dom = vec![usize::MAX; n_cells];
        for (c, dom_out) in cell_dom.iter_mut().enumerate().take(n_cells) {
            if !region.mask[c] {
                continue;
            }
            let lbl = self.mesh.cell_labels[c];
            let d_idx = match id_to_dom.get(&lbl) {
                Some(&d) => d,
                None => {
                    let d = names.len();
                    names.push(format!("domain {lbl}"));
                    cells.push(0);
                    id_to_dom.insert(lbl, d);
                    d
                }
            };
            *dom_out = d_idx;
            cells[d_idx] += 1;
        }
        Ok(DomainMap { names, cells, cell_dom })
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
                Kind::Normal(_) => out.fill(C::new(0.0, 0.0)),
            }
        }
    }

    /// The value of every leaf at one point of the cell: the interpolant of its nodal values, or the cell's own value for a constant one.
    fn at_point(&self, tree: &Tree, loc: &[C], n: &[f64], normal: Option<&[f64; 3]>, work: &mut Vec<C>) -> C {
        work.clear();
        for (i, l) in self.leaves.iter().enumerate() {
            let nodes = &loc[i * self.nloc..(i + 1) * self.nloc];
            work.push(match &l.kind {
                Kind::Const(_) => nodes[0],
                Kind::Normal(a) => C::new(normal.map_or(f64::NAN, |v| v[*a]), 0.0),
                _ => interp(n, nodes),
            });
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
        self.reject_normal()?;
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
                    let mut f = self.at_point(tree, &loc, &el.n[q], None, &mut work);
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

    fn fast_leaves(&self) -> (Vec<FastLeaf<'_>>, bool) {
        let mut has_slow = false;
        let v = self
            .leaves
            .iter()
            .map(|l| match &l.kind {
                Kind::Cell(k) => {
                    let f = &l.res.solutions[l.sol].cell_fields[*k];
                    if f.im.is_empty() {
                        FastLeaf::CellReal(&f.re)
                    } else {
                        FastLeaf::CellComplex(&f.re, &f.im)
                    }
                }
                Kind::Const(vals) => FastLeaf::Const(vals),
                _ => {
                    has_slow = true;
                    FastLeaf::Slow
                }
            })
            .collect();
        (v, has_slow)
    }

    #[inline]
    fn set_const_leaves(fast: &[FastLeaf<'_>], c: usize, work: &mut [C]) {
        for (i, fl) in fast.iter().enumerate() {
            if let FastLeaf::Const(vals) = fl {
                work[i] = C::new(vals[c], 0.0);
            }
        }
    }

    #[inline(always)]
    fn interp_fast_leaves(
        &self,
        fast: &[FastLeaf<'_>],
        c: usize,
        n: &[f64],
        loc: &[C],
        normal: Option<&[f64; 3]>,
        work: &mut [C],
    ) {
        let nloc = self.nloc;
        let base = c * nloc;
        for (i, fl) in fast.iter().enumerate() {
            match fl {
                FastLeaf::CellReal(re) => {
                    work[i] = C::new(dot_n(n, &re[base..base + nloc]), 0.0);
                }
                FastLeaf::CellComplex(re, im) => {
                    work[i] = dot_n_complex(n, &re[base..base + nloc], &im[base..base + nloc]);
                }
                FastLeaf::Const(_) => {}
                FastLeaf::Slow => {
                    work[i] = match &self.leaves[i].kind {
                        Kind::Normal(a) => C::new(normal.map_or(f64::NAN, |v| v[*a]), 0.0),
                        _ => interp(n, &loc[i * nloc..(i + 1) * nloc]),
                    };
                }
            }
        }
    }

    /// Single-pass bit-exact batch volume integral over all `trees` and optional `dom_map` (`by_region`).
    /// Preserves the exact `BLOCK = 2048` summation order and `geo.vol[c] * cre` grouping so every result matches standalone
    /// `integrate` bit-for-bit (`to_bits()`) across any thread count.
    fn integrate_batch(
        &self,
        trees: &[Tree],
        region: &Region,
        dom_map: Option<&DomainMap>,
        revolution: bool,
    ) -> Result<BatchSum, Error> {
        self.reject_normal()?;
        let dim = self.mesh.dim;
        if revolution && dim != 2 {
            return invalid("`revolution` applies to a 2D (r, z) mesh");
        }
        let nthreads = femcore::par::threads();
        let el = RefElement::new(dim, self.order, 4);
        let quad = quad_rule(dim, 4);
        let nq = el.w.len();
        let nloc = self.nloc;
        let n_trees = trees.len();
        let n_doms = dom_map.map_or(0, |m| m.names.len());
        let nb = self.mesh.n_cells().div_ceil(BLOCK);
        let (fast_leaves, has_slow) = self.fast_leaves();
        let dag = FlatDag::compile(trees, self.leaves.len());
        let work_len = dag.as_ref().map_or(self.leaves.len(), FlatDag::reg_len);
        let block = |b: usize| -> BatchSum {
            let mut total = vec![(0.0, 0.0); n_trees];
            let mut measure = 0.0;
            let mut region_totals = vec![vec![(0.0, 0.0); n_trees]; n_doms];
            let mut region_measures = vec![0.0; n_doms];
            let mut loc = if has_slow { vec![C::new(0.0, 0.0); self.leaves.len() * nloc] } else { vec![] };
            let mut dofs = vec![];
            let mut work = vec![C::new(0.0, 0.0); work_len];
            let mut cell_acc = vec![(0.0, 0.0); n_trees];
            for c in b * BLOCK..((b + 1) * BLOCK).min(self.mesh.n_cells()) {
                if !region.mask[c] {
                    continue;
                }
                if has_slow {
                    self.gather(c, &mut dofs, &mut loc);
                }
                Self::set_const_leaves(&fast_leaves, c, &mut work);
                cell_acc.fill((0.0, 0.0));
                let mut cm = 0.0;
                let v = self.mesh.cell(c);
                for q in 0..nq {
                    let n = &el.n[q];
                    self.interp_fast_leaves(&fast_leaves, c, n, &loc, None, &mut work);
                    let mut w = el.w[q];
                    let rho = if revolution {
                        let mut r = 0.0;
                        for (k, &vk) in v.iter().enumerate() {
                            r += quad.lam[q][k] * self.mesh.point(vk as usize)[0];
                        }
                        let rho = 2.0 * std::f64::consts::PI * r;
                        w *= rho;
                        rho
                    } else {
                        1.0
                    };
                    let wq = el.w[q];
                    cm += w;
                    if let Some(d) = &dag {
                        d.eval_all(&mut work);
                        if revolution {
                            for (t, &root) in d.roots.iter().enumerate() {
                                let f = work[root];
                                cell_acc[t].0 += wq * (f.re * rho);
                                cell_acc[t].1 += wq * (f.im * rho);
                            }
                        } else {
                            for (t, &root) in d.roots.iter().enumerate() {
                                let f = work[root];
                                cell_acc[t].0 += wq * f.re;
                                cell_acc[t].1 += wq * f.im;
                            }
                        }
                    } else {
                        for (t, tree) in trees.iter().enumerate() {
                            let mut f = tree.eval(&work);
                            if revolution {
                                f = C::new(f.re * rho, f.im * rho);
                            }
                            cell_acc[t].0 += wq * f.re;
                            cell_acc[t].1 += wq * f.im;
                        }
                    }
                }
                let vol = cell_vol(&self.mesh, c);
                let dm = vol * cm;
                measure += dm;
                if let Some(dm_map) = dom_map {
                    let d = dm_map.cell_dom[c];
                    region_measures[d] += dm;
                    let rt = &mut region_totals[d];
                    for t in 0..n_trees {
                        let dv0 = vol * cell_acc[t].0;
                        let dv1 = vol * cell_acc[t].1;
                        total[t].0 += dv0;
                        total[t].1 += dv1;
                        rt[t].0 += dv0;
                        rt[t].1 += dv1;
                    }
                } else {
                    for t in 0..n_trees {
                        total[t].0 += vol * cell_acc[t].0;
                        total[t].1 += vol * cell_acc[t].1;
                    }
                }
            }
            BatchSum { total, measure, region_totals, region_measures }
        };
        let parts: Vec<BatchSum> = if nthreads <= 1 || nb <= 1 {
            (0..nb).map(block).collect()
        } else {
            (0..nb).into_par_iter().map(block).collect()
        };
        let mut out = BatchSum {
            total: vec![(0.0, 0.0); n_trees],
            measure: 0.0,
            region_totals: vec![vec![(0.0, 0.0); n_trees]; n_doms],
            region_measures: vec![0.0; n_doms],
        };
        for p in parts {
            out.measure += p.measure;
            for t in 0..n_trees {
                out.total[t].0 += p.total[t].0;
                out.total[t].1 += p.total[t].1;
            }
            for d in 0..n_doms {
                out.region_measures[d] += p.region_measures[d];
                for t in 0..n_trees {
                    out.region_totals[d][t].0 += p.region_totals[d][t].0;
                    out.region_totals[d][t].1 += p.region_totals[d][t].1;
                }
            }
        }
        if out.total.iter().any(|(re, im)| !(re.is_finite() && im.is_finite()))
            || out.region_totals.iter().flatten().any(|(re, im)| !(re.is_finite() && im.is_finite()))
        {
            return Err(Error::NonFinite(
                "the integral is not finite (a division by zero, or a field with NaN)".into(),
            ));
        }
        Ok(out)
    }

    /// Single-pass bit-exact confinement integral: accumulates `total` over all cells and `inside[r]` over each target region `r`
    /// in ONE `BLOCK = 2048` pass with exact `to_bits()` parity to standalone `integrate`.
    fn confinement_batch(
        &self,
        power_trees: &[Tree],
        regions: &[Region],
        revolution: bool,
    ) -> Result<ConfinementBatch, Error> {
        self.reject_normal()?;
        let dim = self.mesh.dim;
        if revolution && dim != 2 {
            return invalid("`revolution` applies to a 2D (r, z) mesh");
        }
        let nthreads = femcore::par::threads();
        let el = RefElement::new(dim, self.order, 4);
        let quad = quad_rule(dim, 4);
        let nq = el.w.len();
        let nloc = self.nloc;
        let n_sols = power_trees.len();
        let n_regs = regions.len();
        let nb = self.mesh.n_cells().div_ceil(BLOCK);
        let (fast_leaves, has_slow) = self.fast_leaves();
        let dag = FlatDag::compile(power_trees, self.leaves.len());
        let work_len = dag.as_ref().map_or(self.leaves.len(), FlatDag::reg_len);
        let block = |b: usize| -> ConfinementBatch {
            let mut total = vec![0.0; n_sols];
            let mut inside = vec![vec![0.0; n_sols]; n_regs];
            let mut loc = if has_slow { vec![C::new(0.0, 0.0); self.leaves.len() * nloc] } else { vec![] };
            let mut dofs = vec![];
            let mut work = vec![C::new(0.0, 0.0); work_len];
            let mut cell_acc = vec![0.0; n_sols];
            for c in b * BLOCK..((b + 1) * BLOCK).min(self.mesh.n_cells()) {
                if has_slow {
                    self.gather(c, &mut dofs, &mut loc);
                }
                Self::set_const_leaves(&fast_leaves, c, &mut work);
                cell_acc.fill(0.0);
                let v = self.mesh.cell(c);
                for q in 0..nq {
                    let n = &el.n[q];
                    self.interp_fast_leaves(&fast_leaves, c, n, &loc, None, &mut work);
                    let rho = if revolution {
                        let mut r = 0.0;
                        for (k, &vk) in v.iter().enumerate() {
                            r += quad.lam[q][k] * self.mesh.point(vk as usize)[0];
                        }
                        2.0 * std::f64::consts::PI * r
                    } else {
                        1.0
                    };
                    let wq = el.w[q];
                    if let Some(d) = &dag {
                        d.eval_all(&mut work);
                        if revolution {
                            for (s, &root) in d.roots.iter().enumerate() {
                                let f = work[root];
                                cell_acc[s] += wq * (f.re * rho);
                            }
                        } else {
                            for (s, &root) in d.roots.iter().enumerate() {
                                let f = work[root];
                                cell_acc[s] += wq * f.re;
                            }
                        }
                    } else {
                        for (s, tree) in power_trees.iter().enumerate() {
                            let mut f = tree.eval(&work);
                            if revolution {
                                f = C::new(f.re * rho, f.im * rho);
                            }
                            cell_acc[s] += wq * f.re;
                        }
                    }
                }
                let vol = cell_vol(&self.mesh, c);
                for s in 0..n_sols {
                    total[s] += vol * cell_acc[s];
                }
                for (r, reg) in regions.iter().enumerate() {
                    if reg.mask[c] {
                        let br = &mut inside[r];
                        for s in 0..n_sols {
                            br[s] += vol * cell_acc[s];
                        }
                    }
                }
            }
            ConfinementBatch { total, inside }
        };
        let parts: Vec<ConfinementBatch> = if nthreads <= 1 || nb <= 1 {
            (0..nb).map(block).collect()
        } else {
            (0..nb).into_par_iter().map(block).collect()
        };
        let mut out = ConfinementBatch { total: vec![0.0; n_sols], inside: vec![vec![0.0; n_sols]; n_regs] };
        for p in parts {
            for s in 0..n_sols {
                out.total[s] += p.total[s];
            }
            for r in 0..n_regs {
                for s in 0..n_sols {
                    out.inside[r][s] += p.inside[r][s];
                }
            }
        }
        if out.total.iter().any(|x| !x.is_finite()) || out.inside.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite(
                "the integral is not finite (a division by zero, or a field with NaN)".into(),
            ));
        }
        Ok(out)
    }

    /// A leaf that is the normal of a boundary has no meaning in a volume integral, an extremum or a probe.
    fn reject_normal(&self) -> Result<(), Error> {
        if self.leaves.iter().any(|l| matches!(l.kind, Kind::Normal(_))) {
            return invalid(
                "`nx`, `ny` and `nz` are the normal of a boundary: they exist in `integrate` with `boundary`",
            );
        }
        Ok(())
    }
    /// The faces of a boundary, each read in one cell. `labels` are the boundary numbers of the name; the facets of the solver mesh that carry
    /// them are the outer boundaries. An interface between two regions (`"oxide|silicon"`) has no facets in the solver mesh, so when the name
    /// carries none and reads `a|b` with `a` and `b` regions, the faces are those shared by a cell of `a` and a cell of `b`, which is what the
    /// geometry names `a|b`. A face with one cell is read there. A face with two cells is read in the one whose domain is in `side`; with no
    /// `side`, or with both cells in it, that is an error, because the fields and the materials are not continuous across an interface.
    fn boundary_faces(&self, labels: &[i32], side: Option<&SideSpec>, name: &str) -> Result<Vec<Face>, Error> {
        use std::collections::HashMap;
        let dim = self.mesh.dim;
        let key = |v: &[u32]| -> [u32; 3] {
            let mut k = [u32::MAX; 3];
            k[..v.len()].copy_from_slice(v);
            k[..v.len()].sort_unstable();
            k
        };
        let pair: Option<(Vec<i32>, Vec<i32>)> = match name.split_once('|') {
            Some((a, b)) if facets_with(&self.mesh, labels).next().is_none() && !b.contains('|') => {
                match (self.model.domains(a.trim()), self.model.domains(b.trim())) {
                    (Ok(a), Ok(b)) if !a.is_empty() && !b.is_empty() => Some((a, b)),
                    _ => None,
                }
            }
            _ => None,
        };
        // every face of the cells that may belong to the boundary, with the cells next to it, in the order the faces are first met
        let mut faces: Vec<Face> = vec![];
        let mut index: HashMap<[u32; 3], usize> = HashMap::new();
        if pair.is_none() {
            let chosen: Vec<usize> = facets_with(&self.mesh, labels).collect();
            if chosen.is_empty() {
                return invalid(format!(
                    "the boundary '{name}' has no facets in the mesh of this result: an outer boundary, or an interface `a|b` between two regions, can be integrated over"
                ));
            }
            for f in chosen {
                let v = self.mesh.facet(f);
                let mut vs = [0u32; 3];
                vs[..dim].copy_from_slice(v);
                if index.insert(key(v), faces.len()).is_none() {
                    faces.push(Face { v: vs, cells: vec![] });
                }
            }
        }
        let mut face = [0u32; 3];
        for c in 0..self.mesh.n_cells() {
            let v = self.mesh.cell(c);
            for m in 0..=dim {
                let mut n = 0;
                for (a, &va) in v.iter().enumerate() {
                    if a != m {
                        face[n] = va;
                        n += 1;
                    }
                }
                let k = key(&face[..dim]);
                match index.get(&k) {
                    Some(&i) => faces[i].cells.push((c, m)),
                    None if pair.is_some() => {
                        index.insert(k, faces.len());
                        faces.push(Face { v: face, cells: vec![(c, m)] });
                    }
                    None => {}
                }
            }
        }
        if let Some((a, b)) = &pair {
            let label = |c: usize| self.mesh.cell_labels[c];
            faces.retain(|f| match f.cells.as_slice() {
                [(c, _), (d, _)] => {
                    (a.contains(&label(*c)) && b.contains(&label(*d)))
                        || (b.contains(&label(*c)) && a.contains(&label(*d)))
                }
                _ => false,
            });
            if faces.is_empty() {
                return invalid(format!("the regions of '{name}' do not touch in the mesh of this result"));
            }
        }
        let regions = self.model.regions();
        let domain =
            |l: i32| regions.iter().find(|(id, _)| *id == l).map_or_else(|| format!("domain {l}"), |(_, n)| n.clone());
        let mut out = Vec::with_capacity(faces.len());
        for f in faces {
            let near = || -> String {
                f.cells.iter().map(|(c, _)| domain(self.mesh.cell_labels[*c])).collect::<Vec<_>>().join(" and ")
            };
            let picked = match side {
                None if f.cells.len() == 1 => f.cells[0],
                None => {
                    return invalid(format!(
                        "the boundary '{name}' is an interface between {}: the fields and the materials differ across it, so give `side`, the region to read them in (the normal points out of it)",
                        near()
                    ))
                }
                Some(SideSpec::Axis(ax, pos)) => {
                    // A cell on the +ax side of the facet has its outward normal across the facet pointing in -ax (n[ax] < 0).
                    let want_sign = if *pos { -1.0 } else { 1.0 };
                    let pick: Vec<(usize, usize)> = f
                        .cells
                        .iter()
                        .copied()
                        .filter(|&(c, opp)| self.face_geometry(&f.v, c, opp).1[*ax] * want_sign > 1e-6)
                        .collect();
                    match pick.as_slice() {
                        [one, ..] => *one,
                        [] => {
                            return invalid(format!(
                                "the boundary '{name}' has no adjacent cell on the `{}{}` side",
                                if *pos { "+" } else { "-" },
                                ["x", "y", "z"][*ax]
                            ))
                        }
                    }
                }
                Some(SideSpec::Domains(s)) => {
                    let pick: Vec<(usize, usize)> = f
                        .cells
                        .iter()
                        .copied()
                        .filter(|(c, _)| s.contains(&self.mesh.cell_labels[*c]))
                        .collect();
                    match pick.as_slice() {
                        [one] => *one,
                        [] => {
                            return invalid(format!(
                                "the boundary '{name}' has a face next to {} only, not to the `side` region you gave",
                                near()
                            ))
                        }
                        [a, b] if self.mesh.cell_labels[a.0] == self.mesh.cell_labels[b.0] => {
                            // Sheet embedded inside a single homogeneous region (rb == ra): pick the -z cell (outward normal +z)
                            if self.face_geometry(&f.v, a.0, a.1).1[2] >= 0.0 {
                                *a
                            } else {
                                *b
                            }
                        }
                        _ => {
                            return invalid(format!(
                                "the boundary '{name}' has a face with {} both in `side`: name the region on one side only",
                                near()
                            ))
                        }
                    }
                }
            };
            out.push(Face { v: f.v, cells: vec![picked] });
        }
        Ok(out)
    }

    /// The measure of a face (an edge's length in 2D, a triangle's area in 3D) and its unit normal in wafer axes, pointing out of the cell
    /// whose vertex `opposite` (local number) is not on the face.
    fn face_geometry(&self, v: &[u32; 3], cell: usize, opposite: usize) -> (f64, [f64; 3]) {
        let dim = self.mesh.dim;
        let p0 = self.mesh.point(v[0] as usize);
        let n = if dim == 2 {
            let p1 = self.mesh.point(v[1] as usize);
            [p1[1] - p0[1], -(p1[0] - p0[0]), 0.0]
        } else {
            let (p1, p2) = (self.mesh.point(v[1] as usize), self.mesh.point(v[2] as usize));
            let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        };
        let len = n[..dim].iter().map(|x| x * x).sum::<f64>().sqrt();
        let far = self.mesh.point(self.mesh.cell(cell)[opposite] as usize);
        let toward: f64 = (0..dim).map(|k| n[k] * (far[k] - p0[k])).sum();
        let s = if toward > 0.0 { -1.0 } else { 1.0 } / len;
        let m = [n[0] * s, n[1] * s, n[2] * s];
        let measure = if dim == 2 { len } else { 0.5 * len };
        (measure, femcore::section::wafer_point(&m[..dim], dim))
    }

    /// `sum over the faces of tree dS` (module comment, "Boundary integrals"): the facet rule at points mapped into the barycentric
    /// coordinates of the cell each face is read in; blocks of `BLOCK` faces in face order, so the bits do not depend on the threads.
    fn integrate_boundary(&self, tree: &Tree, faces: &[Face], revolution: bool) -> Result<Sum, Error> {
        let dim = self.mesh.dim;
        if revolution && dim != 2 {
            return invalid("`revolution` applies to a 2D (r, z) mesh");
        }
        femcore::par::threads();
        let rule = facet_rule(dim - 1);
        let nloc = self.nloc;
        let nb = faces.len().div_ceil(BLOCK);
        let block = |b: usize| -> Sum {
            let mut s = Sum { re: 0.0, im: 0.0, measure: 0.0 };
            let mut loc = vec![C::new(0.0, 0.0); self.leaves.len() * nloc];
            let mut dofs = vec![];
            let mut work = Vec::with_capacity(self.leaves.len());
            let mut lam = vec![0.0; dim + 1];
            for face in &faces[b * BLOCK..((b + 1) * BLOCK).min(faces.len())] {
                let (c, opposite) = face.cells[0];
                self.gather(c, &mut dofs, &mut loc);
                let (measure, normal) = self.face_geometry(&face.v, c, opposite);
                let fv = &face.v[..dim];
                let v = self.mesh.cell(c);
                let (mut cre, mut cim, mut cm) = (0.0, 0.0, 0.0);
                for (q, lf) in rule.lam.iter().enumerate() {
                    lam.fill(0.0);
                    for (j, &vj) in fv.iter().enumerate() {
                        let at = v.iter().position(|&x| x == vj).expect("a face vertex is a vertex of its cell");
                        lam[at] = lf[j];
                    }
                    let (n, _) = shape(&lam, self.order);
                    let mut g = self.at_point(tree, &loc, &n, Some(&normal), &mut work);
                    let mut w = rule.w[q];
                    if revolution {
                        let r: f64 =
                            fv.iter().enumerate().map(|(j, &vj)| lf[j] * self.mesh.point(vj as usize)[0]).sum();
                        let rho = 2.0 * std::f64::consts::PI * r;
                        g = C::new(g.re * rho, g.im * rho);
                        w *= rho;
                    }
                    cre += rule.w[q] * g.re;
                    cim += rule.w[q] * g.im;
                    cm += w;
                }
                s.re += measure * cre;
                s.im += measure * cim;
                s.measure += measure * cm;
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

    /// Single-pass bit-exact batch boundary integral over all `trees` across `faces` in `BLOCK = 2048` face chunks.
    fn integrate_boundary_batch(&self, trees: &[Tree], faces: &[Face], revolution: bool) -> Result<BatchSum, Error> {
        let dim = self.mesh.dim;
        if revolution && dim != 2 {
            return invalid("`revolution` applies to a 2D (r, z) mesh");
        }
        femcore::par::threads();
        let rule = facet_rule(dim - 1);
        let nloc = self.nloc;
        let n_trees = trees.len();
        let nb = faces.len().div_ceil(BLOCK);
        let block = |b: usize| -> BatchSum {
            let mut total = vec![(0.0, 0.0); n_trees];
            let mut measure = 0.0;
            let mut loc = vec![C::new(0.0, 0.0); self.leaves.len() * nloc];
            let mut dofs = vec![];
            let mut work = vec![C::new(0.0, 0.0); self.leaves.len()];
            let mut lam = vec![0.0; dim + 1];
            let mut face_acc = vec![(0.0, 0.0); n_trees];
            for face in &faces[b * BLOCK..((b + 1) * BLOCK).min(faces.len())] {
                let (c, opposite) = face.cells[0];
                self.gather(c, &mut dofs, &mut loc);
                let (f_measure, normal) = self.face_geometry(&face.v, c, opposite);
                let fv = &face.v[..dim];
                let v = self.mesh.cell(c);
                face_acc.fill((0.0, 0.0));
                let mut cm = 0.0;
                for (q, lf) in rule.lam.iter().enumerate() {
                    lam.fill(0.0);
                    for (j, &vj) in fv.iter().enumerate() {
                        let at = v.iter().position(|&x| x == vj).expect("a face vertex is a vertex of its cell");
                        lam[at] = lf[j];
                    }
                    let (n, _) = shape(&lam, self.order);
                    for (i, l) in self.leaves.iter().enumerate() {
                        let nodes = &loc[i * nloc..(i + 1) * nloc];
                        work[i] = match &l.kind {
                            Kind::Const(_) => nodes[0],
                            Kind::Normal(a) => C::new(normal[*a], 0.0),
                            _ => interp(&n, nodes),
                        };
                    }
                    let mut w = rule.w[q];
                    let rho = if revolution {
                        let r: f64 =
                            fv.iter().enumerate().map(|(j, &vj)| lf[j] * self.mesh.point(vj as usize)[0]).sum();
                        let rho = 2.0 * std::f64::consts::PI * r;
                        w *= rho;
                        rho
                    } else {
                        1.0
                    };
                    let wq = rule.w[q];
                    cm += w;
                    for (t, tree) in trees.iter().enumerate() {
                        let mut g = tree.eval(&work);
                        if revolution {
                            g = C::new(g.re * rho, g.im * rho);
                        }
                        face_acc[t].0 += wq * g.re;
                        face_acc[t].1 += wq * g.im;
                    }
                }
                measure += f_measure * cm;
                for t in 0..n_trees {
                    total[t].0 += f_measure * face_acc[t].0;
                    total[t].1 += f_measure * face_acc[t].1;
                }
            }
            BatchSum { total, measure, region_totals: vec![], region_measures: vec![] }
        };
        let parts: Vec<BatchSum> = (0..nb).into_par_iter().map(block).collect();
        let mut out =
            BatchSum { total: vec![(0.0, 0.0); n_trees], measure: 0.0, region_totals: vec![], region_measures: vec![] };
        for p in parts {
            out.measure += p.measure;
            for t in 0..n_trees {
                out.total[t].0 += p.total[t].0;
                out.total[t].1 += p.total[t].1;
            }
        }
        if out.total.iter().any(|(re, im)| !(re.is_finite() && im.is_finite())) {
            return Err(Error::NonFinite(
                "the integral is not finite (a division by zero, or a field with NaN)".into(),
            ));
        }
        Ok(out)
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
        self.reject_normal()?;
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
        let raw_expr = args
            .get("expr")
            .ok_or_else(|| Error::invalid("`expr` is required: a field name, a number or an expression"))?;
        let batch_expr = BatchExpr::from_value(raw_expr)?;
        let by_region = match args.get("by_region") {
            None | Some(Value::Null) => false,
            Some(Value::Bool(b)) => *b,
            Some(_) => return invalid("`by_region` must be true or false"),
        };
        let d = self.defaults(args)?;
        let (sol_indices, has_solutions_arg) = self.parse_solutions(args, &d)?;
        let is_batch = has_solutions_arg || by_region || batch_expr.is_multi();
        let rev = Self::revolution(args);
        if !is_batch {
            let (tree, p) = self.one_expr(args, "expr")?;
            if let Some(b) = args.get("boundary").filter(|v| !v.is_null()) {
                return self.integrate_over_boundary(args, &tree, &p, b, rev);
            }
            if args.get("side").is_some_and(|v| !v.is_null()) {
                return invalid("`side` goes with `boundary`: the region an interface is read in");
            }
            let region = p.region(Self::region_arg(args))?;
            let s = p.integrate(&tree, &region, rev)?;
            let real = tree.is_real_valued(&p.complex_flags());
            let mut out = json!({
                "value": if real { json!(s.re) } else { pair(s.re, s.im) }, "complex": !real,
                "measure": s.measure, "cells": region.cells, "revolution": rev, "region": Self::region_arg(args),
            });
            self.axisym_note(&p, rev, &mut out);
            self.decorate(&mut out, &p);
            return Ok(out);
        }

        let expr_vals = batch_expr.values();
        let (trees, p) = self.expressions_for_solutions(&d, &sol_indices, has_solutions_arg, &expr_vals)?;
        if let Some(b) = args.get("boundary").filter(|v| !v.is_null()) {
            if by_region {
                return invalid(
                    "give `boundary` or `by_region`, not both: a boundary integral sums over facets, not regions",
                );
            }
            return self.integrate_over_boundary_batch(args, &batch_expr, &sol_indices, &trees, &p, b);
        }
        if args.get("side").is_some_and(|v| !v.is_null()) {
            return invalid("`side` goes with `boundary`: the region an interface is read in");
        }
        let reg_arg = Self::region_arg(args);
        let region = p.region(reg_arg)?;
        let dom_map = if by_region { Some(p.build_domain_map(&region, reg_arg)?) } else { None };
        let batch = p.integrate_batch(&trees, &region, dom_map.as_ref(), rev)?;
        let flags = p.complex_flags();
        let real_flags: Vec<bool> = trees.iter().map(|t| t.is_real_valued(&flags)).collect();
        let n_exprs = expr_vals.len();
        let val_of = |t_idx: usize, (re, im): (f64, f64)| -> Value {
            if real_flags[t_idx] {
                json!(re)
            } else {
                pair(re, im)
            }
        };
        let pack_values = |base: usize, slice: &[(f64, f64)]| -> Value {
            match &batch_expr {
                BatchExpr::Single(_) => val_of(base, slice[0]),
                BatchExpr::List(_) => Value::Array((0..n_exprs).map(|e| val_of(base + e, slice[e])).collect()),
                BatchExpr::Map(entries) => {
                    let mut m = serde_json::Map::with_capacity(n_exprs);
                    for (e, (k, _)) in entries.iter().enumerate() {
                        m.insert(k.clone(), val_of(base + e, slice[e]));
                    }
                    Value::Object(m)
                }
            }
        };

        let mut sol_entries = Vec::with_capacity(sol_indices.len());
        for (s_idx, &sol_id) in sol_indices.iter().enumerate() {
            let base = s_idx * n_exprs;
            let vals = pack_values(base, &batch.total[base..base + n_exprs]);
            let mut sol_obj = serde_json::Map::new();
            sol_obj.insert("solution".into(), json!(sol_id));
            if matches!(batch_expr, BatchExpr::Single(_)) {
                sol_obj.insert("value".into(), vals.clone());
            }
            sol_obj.insert("values".into(), vals);
            if let Some(dm) = &dom_map {
                let mut reg_map = serde_json::Map::with_capacity(dm.names.len());
                for (d_idx, dom_name) in dm.names.iter().enumerate() {
                    let r_vals = pack_values(base, &batch.region_totals[d_idx][base..base + n_exprs]);
                    reg_map.insert(dom_name.clone(), r_vals);
                }
                sol_obj.insert("regions".into(), Value::Object(reg_map));
            }
            sol_obj.insert("measure".into(), json!(batch.measure));
            sol_obj.insert("cells".into(), json!(region.cells));
            sol_entries.push(Value::Object(sol_obj));
        }

        let mut out = json!({
            "batch": true,
            "by_region": by_region,
            "solutions": sol_entries,
            "measure": batch.measure,
            "cells": region.cells,
            "revolution": rev,
            "region": reg_arg,
        });
        if let Some(dm) = &dom_map {
            let mut rm = serde_json::Map::with_capacity(dm.names.len());
            let mut rc = serde_json::Map::with_capacity(dm.names.len());
            for (d_idx, dom_name) in dm.names.iter().enumerate() {
                rm.insert(dom_name.clone(), json!(batch.region_measures[d_idx]));
                rc.insert(dom_name.clone(), json!(dm.cells[d_idx]));
            }
            out["region_measures"] = Value::Object(rm);
            out["region_cells"] = Value::Object(rc);
        }
        if out["solutions"].as_array().is_some_and(|a| a.len() == 1) {
            let first = out["solutions"][0].clone();
            out["values"] = first["values"].clone();
            if matches!(batch_expr, BatchExpr::Single(_)) {
                out["value"] = first["value"].clone();
                out["complex"] = json!(!real_flags[0]);
            }
            if by_region {
                out["regions"] = first["regions"].clone();
            }
        }
        self.axisym_note(&p, rev, &mut out);
        self.decorate(&mut out, &p);
        Ok(out)
    }

    /// `integrate` over the facets of a named boundary (ADR 90, N2b): a line integral on a 2D mesh, a surface integral on a 3D one.
    fn integrate_over_boundary(&self, args: &Value, tree: &Tree, p: &Prep, boundary: &Value, rev: bool) -> R {
        let (name, side_arg, facets) = Self::resolve_boundary_facets(args, p, boundary)?;
        let s = p.integrate_boundary(tree, &facets, rev)?;
        let real = tree.is_real_valued(&p.complex_flags());
        let mut out = json!({
            "value": if real { json!(s.re) } else { pair(s.re, s.im) }, "complex": !real,
            "measure": s.measure, "facets": facets.len(), "revolution": rev, "boundary": name, "side": side_arg,
            "on": if p.mesh.dim == 2 { "a line (the edges of the boundary)" } else { "a surface (the faces of the boundary)" },
            "normal": "nx, ny, nz: the unit normal in wafer axes, out of the cell each facet is read in",
        });
        self.axisym_note(p, rev, &mut out);
        self.decorate(&mut out, p);
        Ok(out)
    }

    fn integrate_over_boundary_batch(
        &self,
        args: &Value,
        batch_expr: &BatchExpr,
        sol_indices: &[usize],
        trees: &[Tree],
        p: &Prep,
        boundary: &Value,
    ) -> R {
        let rev = Self::revolution(args);
        let (name, side_arg, facets) = Self::resolve_boundary_facets(args, p, boundary)?;
        let batch = p.integrate_boundary_batch(trees, &facets, rev)?;
        let flags = p.complex_flags();
        let real_flags: Vec<bool> = trees.iter().map(|t| t.is_real_valued(&flags)).collect();
        let n_exprs = batch_expr.values().len();
        let val_of = |t_idx: usize, (re, im): (f64, f64)| -> Value {
            if real_flags[t_idx] {
                json!(re)
            } else {
                pair(re, im)
            }
        };
        let pack_values = |base: usize, slice: &[(f64, f64)]| -> Value {
            match batch_expr {
                BatchExpr::Single(_) => val_of(base, slice[0]),
                BatchExpr::List(_) => Value::Array((0..n_exprs).map(|e| val_of(base + e, slice[e])).collect()),
                BatchExpr::Map(entries) => {
                    let mut m = serde_json::Map::with_capacity(n_exprs);
                    for (e, (k, _)) in entries.iter().enumerate() {
                        m.insert(k.clone(), val_of(base + e, slice[e]));
                    }
                    Value::Object(m)
                }
            }
        };
        let mut sol_entries = Vec::with_capacity(sol_indices.len());
        for (s_idx, &sol_id) in sol_indices.iter().enumerate() {
            let base = s_idx * n_exprs;
            let vals = pack_values(base, &batch.total[base..base + n_exprs]);
            let mut sol_obj = serde_json::Map::new();
            sol_obj.insert("solution".into(), json!(sol_id));
            if matches!(batch_expr, BatchExpr::Single(_)) {
                sol_obj.insert("value".into(), vals.clone());
            }
            sol_obj.insert("values".into(), vals);
            sol_obj.insert("measure".into(), json!(batch.measure));
            sol_obj.insert("facets".into(), json!(facets.len()));
            sol_entries.push(Value::Object(sol_obj));
        }
        let mut out = json!({
            "batch": true,
            "by_region": false,
            "solutions": sol_entries,
            "measure": batch.measure,
            "facets": facets.len(),
            "revolution": rev,
            "boundary": name,
            "side": side_arg,
            "on": if p.mesh.dim == 2 { "a line (the edges of the boundary)" } else { "a surface (the faces of the boundary)" },
            "normal": "nx, ny, nz: the unit normal in wafer axes, out of the cell each facet is read in",
        });
        if out["solutions"].as_array().is_some_and(|a| a.len() == 1) {
            let first = out["solutions"][0].clone();
            out["values"] = first["values"].clone();
            if matches!(batch_expr, BatchExpr::Single(_)) {
                out["value"] = first["value"].clone();
                out["complex"] = json!(!real_flags[0]);
            }
        }
        self.axisym_note(p, rev, &mut out);
        self.decorate(&mut out, p);
        Ok(out)
    }

    fn resolve_boundary_facets<'a>(
        args: &'a Value,
        p: &Prep,
        boundary: &'a Value,
    ) -> Result<(&'a str, Option<&'a str>, Vec<Face>), Error> {
        let name = boundary.as_str().filter(|s| !s.trim().is_empty()).ok_or_else(|| {
            Error::invalid("`boundary` is the name of a boundary, or a selector that picks boundaries (e.g. \"outer & plane(z = 0)\")")
        })?;
        if Self::region_arg(args).is_some_and(|r| !r.trim().is_empty()) {
            return invalid(
                "give `boundary` or `region`, not both: a boundary integral has no region (an interface takes `side`)",
            );
        }
        let labels = p.model.boundaries(name)?;
        let side_arg = match args.get("side") {
            None | Some(Value::Null) => None,
            Some(Value::String(s)) if !s.trim().is_empty() => Some(s.as_str()),
            Some(_) => {
                return invalid("`side` is the name of a region (or \"+z\" / \"-z\" for a sheet inside one region)")
            }
        };
        let side = match side_arg.map(str::trim) {
            None => None,
            Some("+x") => Some(SideSpec::Axis(0, true)),
            Some("-x") => Some(SideSpec::Axis(0, false)),
            Some("+y") => Some(SideSpec::Axis(1, true)),
            Some("-y") => Some(SideSpec::Axis(1, false)),
            Some("+z") => Some(SideSpec::Axis(2, true)),
            Some("-z") => Some(SideSpec::Axis(2, false)),
            Some(s) => Some(SideSpec::Domains(p.model.domains(s)?)),
        };
        let facets = p.boundary_faces(&labels, side.as_ref(), name)?;
        Ok((name, side_arg, facets))
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
                    let v = p.at_point(tree, &local, &n, None, &mut work);
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
        p.reject_normal()?;
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
        p.reject_normal()?;
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
        let by_region = match args.get("by_region") {
            None | Some(Value::Null) => false,
            Some(Value::Bool(b)) => *b,
            Some(_) => return invalid("`by_region` must be true or false"),
        };
        if args.get("region").is_some_and(|v| !v.is_null()) && args.get("regions").is_some_and(|v| !v.is_null()) {
            return invalid("give `region` or `regions`, not both");
        }
        let explicit_regions: Option<Vec<String>> = match args.get("regions").filter(|v| !v.is_null()) {
            Some(Value::Array(arr)) => {
                if arr.is_empty() {
                    return invalid("`regions` cannot be empty");
                }
                let mut v = Vec::with_capacity(arr.len());
                for item in arr {
                    let s = item
                        .as_str()
                        .filter(|s| !s.trim().is_empty())
                        .ok_or_else(|| Error::invalid("`regions` must be a list of region names"))?;
                    v.push(s.to_string());
                }
                Some(v)
            }
            Some(_) => return invalid("`regions` must be a list of region names"),
            None => match args.get("region") {
                Some(Value::Array(arr)) => {
                    if arr.is_empty() {
                        return invalid("`region` list cannot be empty");
                    }
                    let mut v = Vec::with_capacity(arr.len());
                    for item in arr {
                        let s = item
                            .as_str()
                            .filter(|s| !s.trim().is_empty())
                            .ok_or_else(|| Error::invalid("`region` list must contain region names"))?;
                        v.push(s.to_string());
                    }
                    Some(v)
                }
                _ => None,
            },
        };
        let d = self.defaults(args)?;
        let (sol_indices, has_solutions_arg) = self.parse_solutions(args, &d)?;
        let is_batch = has_solutions_arg || by_region || explicit_regions.is_some();
        let (all_trees, p) = self.expressions_for_solutions(&d, &sol_indices, has_solutions_arg, &comps)?;
        let flags = p.complex_flags();
        let n_comps = comps.len();
        let mut power_trees = Vec::with_capacity(sol_indices.len());
        for s_idx in 0..sol_indices.len() {
            let mut trees: Vec<Tree> = all_trees[s_idx * n_comps..(s_idx + 1) * n_comps].to_vec();
            let weight = weight_arg.map(|_| trees.pop().expect("the weight"));
            if weight.as_ref().is_some_and(|w| !w.is_real_valued(&flags)) {
                return invalid("`weight` must be real: wrap a complex formula in abs, abs2, re or im");
            }
            let mut it = trees.into_iter().map(|t| Tree::Call(Fun::Abs2, Box::new(t)));
            let first = it.next().expect("a component");
            let mut power = it.fold(first, |acc, t| Tree::Bin(Op::Add, Box::new(acc), Box::new(t)));
            if let Some(w) = weight {
                power = Tree::Bin(Op::Mul, Box::new(w), Box::new(power));
            }
            power_trees.push(power);
        }
        let rev = Self::revolution(args);
        let meaning = if weight_arg.is_some() {
            "the integral of weight * sum_k |c_k|^2 over the region divided by the same over the whole mesh"
        } else {
            "the integral of sum_k |c_k|^2 over the region divided by the same over the whole mesh (no material weight)"
        };
        if !is_batch {
            let part = p.region(Some(
                Self::region_arg(args)
                    .ok_or_else(|| Error::invalid("`region` is required: the region the fraction is of"))?,
            ))?;
            let batch = p.confinement_batch(&power_trees, std::slice::from_ref(&part), rev)?;
            let inside = batch.inside[0][0];
            let total = batch.total[0];
            let mut out = json!({
                "fraction": inside / total, "inside": inside, "total": total, "cells": part.cells,
                "revolution": rev, "region": Self::region_arg(args),
                "weighted": weight_arg.is_some(),
                "meaning": meaning,
            });
            self.axisym_note(&p, rev, &mut out);
            self.decorate(&mut out, &p);
            return Ok(out);
        }

        let target_regions: Vec<(String, Region)> = if let Some(names) = &explicit_regions {
            let mut v = Vec::with_capacity(names.len());
            for name in names {
                v.push((name.clone(), p.region(Some(name))?));
            }
            v
        } else if by_region {
            let reg_arg = Self::region_arg(args);
            let filter_reg = p.region(reg_arg)?;
            let dm = p.build_domain_map(&filter_reg, reg_arg)?;
            let n_cells = p.mesh.n_cells();
            dm.names
                .into_iter()
                .enumerate()
                .map(|(d_idx, name)| {
                    let mask: Vec<bool> = (0..n_cells).map(|c| dm.cell_dom[c] == d_idx).collect();
                    let cells = dm.cells[d_idx];
                    (name, Region { mask, cells })
                })
                .collect()
        } else {
            let r_name = Self::region_arg(args)
                .ok_or_else(|| Error::invalid("`region` is required: the region the fraction is of"))?;
            vec![(r_name.to_string(), p.region(Some(r_name))?)]
        };

        let regs_only: Vec<Region> =
            target_regions.iter().map(|(_, r)| Region { mask: r.mask.clone(), cells: r.cells }).collect();
        let batch = p.confinement_batch(&power_trees, &regs_only, rev)?;
        let multi_region = by_region || explicit_regions.is_some();
        let mut sol_entries = Vec::with_capacity(sol_indices.len());
        for (s_idx, &sol_id) in sol_indices.iter().enumerate() {
            let tot = batch.total[s_idx];
            let mut sol_obj = serde_json::Map::new();
            sol_obj.insert("solution".into(), json!(sol_id));
            let mut reg_fracs = serde_json::Map::with_capacity(target_regions.len());
            let mut reg_insides = serde_json::Map::with_capacity(target_regions.len());
            for (r_idx, (r_name, _)) in target_regions.iter().enumerate() {
                let ins = batch.inside[r_idx][s_idx];
                reg_fracs.insert(r_name.clone(), json!(ins / tot));
                reg_insides.insert(r_name.clone(), json!(ins));
            }
            if !multi_region {
                let ins = batch.inside[0][s_idx];
                sol_obj.insert("fraction".into(), json!(ins / tot));
                sol_obj.insert("inside".into(), json!(ins));
                sol_obj.insert("total".into(), json!(tot));
                sol_obj.insert("cells".into(), json!(target_regions[0].1.cells));
                sol_obj.insert("region".into(), json!(target_regions[0].0));
                sol_obj.insert("regions".into(), Value::Object(reg_fracs.clone()));
                sol_obj.insert("fractions".into(), Value::Object(reg_fracs));
            } else {
                sol_obj.insert("regions".into(), Value::Object(reg_fracs.clone()));
                sol_obj.insert("fractions".into(), Value::Object(reg_fracs));
                sol_obj.insert("inside".into(), Value::Object(reg_insides));
                sol_obj.insert("total".into(), json!(tot));
            }
            sol_entries.push(Value::Object(sol_obj));
        }

        let mut region_cells = serde_json::Map::with_capacity(target_regions.len());
        for (r_name, reg) in &target_regions {
            region_cells.insert(r_name.clone(), json!(reg.cells));
        }
        let mut out = json!({
            "batch": true,
            "by_region": multi_region,
            "solutions": sol_entries,
            "region_cells": Value::Object(region_cells),
            "revolution": rev,
            "region": Self::region_arg(args),
            "weighted": weight_arg.is_some(),
            "meaning": meaning,
        });
        if out["solutions"].as_array().is_some_and(|a| a.len() == 1) {
            let first = out["solutions"][0].clone();
            out["regions"] = first["regions"].clone();
            out["fractions"] = first["fractions"].clone();
            out["inside"] = first["inside"].clone();
            out["total"] = first["total"].clone();
            if !multi_region {
                out["fraction"] = first["fraction"].clone();
                out["cells"] = first["cells"].clone();
            }
        }
        self.axisym_note(&p, rev, &mut out);
        self.decorate(&mut out, &p);
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrate_works_on_both_sides_of_interior_sheets_and_geometry_names_exposes_sheet_edges() {
        let mut e = Engine::new();
        let mut call =
            |cmd: &str, args: Value| -> Value { e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}")) };
        call("project_new", json!({ "title": "3D sheet integrate", "dim": 3, "h": 40e-6 }));
        call("add_box", json!({ "name": "sub", "min": [-100e-6, -100e-6, -80e-6], "max": [100e-6, 100e-6, 0.0] }));
        call("add_box", json!({ "name": "vac", "min": [-100e-6, -100e-6, 0.0], "max": [100e-6, 100e-6, 80e-6] }));
        call(
            "add_polygons",
            json!({
                "name": "pad_1_poly",
                "polygons": [[[-40e-6, -40e-6], [40e-6, -40e-6], [40e-6, 40e-6], [-40e-6, 40e-6]]]
            }),
        );
        call("add_sheet", json!({ "name": "pad_1", "input": "pad_1_poly", "z": 0.0 }));
        call("add_rect", json!({ "name": "bridge_rect", "min": [-20e-6, -20e-6], "max": [20e-6, 20e-6] }));
        call("add_sheet", json!({ "name": "bridge", "input": "bridge_rect", "z": 30e-6 }));
        call(
            "add_overlay",
            json!({
                "regions": [
                    { "name": "substrate", "input": "sub" },
                    { "name": "vacuum", "input": "vac" }
                ]
            }),
        );

        let names = call("geometry_names", json!({}));
        let bnds: Vec<&str> = names["boundaries"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
        assert!(bnds.contains(&"pad_1") && bnds.contains(&"pad_1:edge"), "boundaries: {bnds:?}");
        assert!(bnds.contains(&"bridge") && bnds.contains(&"bridge:edge"), "boundaries: {bnds:?}");

        call("add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 1 }));
        call("set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "k": 10.0 } }));
        call("set_material", json!({ "physics": "dfs", "region": "vacuum", "material": { "k": 1.0 } }));
        call(
            "add_condition",
            json!({ "physics": "dfs", "on": "pad_1", "kind": "temperature", "params": { "value": 1.0 } }),
        );
        call(
            "add_condition",
            json!({ "physics": "dfs", "on": "outer", "kind": "temperature", "params": { "value": 0.0 } }),
        );
        call("solve_run", json!({}));

        // 1. Sheet on substrate|vacuum interface (z = 0): side="substrate" vs side="vacuum"
        let pad_area = 80e-6 * 80e-6;
        let r_sub = call("integrate", json!({ "expr": "nz", "boundary": "pad_1", "side": "substrate" }));
        let r_vac = call("integrate", json!({ "expr": "nz", "boundary": "pad_1", "side": "vacuum" }));
        // Out of substrate (-z side) points in +z (+pad_area); out of vacuum (+z side) points in -z (-pad_area)
        assert!((r_sub["measure"].as_f64().unwrap() - pad_area).abs() < 1e-12 * pad_area);
        assert!((r_vac["measure"].as_f64().unwrap() - pad_area).abs() < 1e-12 * pad_area);
        assert!((r_sub["value"].as_f64().unwrap() - pad_area).abs() < 1e-11 * pad_area);
        assert!((r_vac["value"].as_f64().unwrap() + pad_area).abs() < 1e-11 * pad_area);

        // Material property k differs across pad_1 (10 in substrate, 1 in vacuum)
        let k_sub = call("integrate", json!({ "expr": "k", "boundary": "pad_1", "side": "substrate" }));
        let k_vac = call("integrate", json!({ "expr": "k", "boundary": "pad_1", "side": "vacuum" }));
        assert!((k_sub["value"].as_f64().unwrap() - 10.0 * pad_area).abs() < 1e-11 * pad_area);
        assert!((k_vac["value"].as_f64().unwrap() - 1.0 * pad_area).abs() < 1e-11 * pad_area);

        // 2. Sheet embedded inside homogeneous vacuum (z = 30e-6, rb == ra): side="-z", side="+z", and side="vacuum"
        let br_area = 40e-6 * 40e-6;
        let br_mz = call("integrate", json!({ "expr": "nz", "boundary": "bridge", "side": "-z" }));
        let br_pz = call("integrate", json!({ "expr": "nz", "boundary": "bridge", "side": "+z" }));
        let br_vac = call("integrate", json!({ "expr": "nz", "boundary": "bridge", "side": "vacuum" }));
        assert!((br_mz["value"].as_f64().unwrap() - br_area).abs() < 1e-11 * br_area);
        assert!((br_pz["value"].as_f64().unwrap() + br_area).abs() < 1e-11 * br_area);
        assert!((br_vac["value"].as_f64().unwrap() - br_area).abs() < 1e-11 * br_area);
    }
}
