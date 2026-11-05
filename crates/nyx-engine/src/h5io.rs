//! Results in HDF5 files (cargo feature `hdf5`): the project, the mesh and every solution of every solved physics with all its fields and
//! cell fields, written compressed and read back exactly (ADR 62). This is the efficient store of the data a solve produces; the JSON files
//! of a project and a run record stay small.
//!
//! Layout, version 3 (numbers as stored by HDF5, little-endian):
//!
//! ```text
//! /                          attrs: format = "nyx-results", version = 3
//! /project                   the project as JSON text (u8 dataset)
//! /mesh/points               [n_points, dim] f64         /mesh/cells        [n_cells, dim+1] u32
//! /mesh/cell_labels          [n_cells] i32               /mesh/facets       [n_facets, dim] u32
//! /mesh/facet_labels         [n_facets] i32
//! /results/<physics>         attrs: study, order, dim, dofs, seconds, n_solutions
//!   scalars, log             JSON text (u8 datasets)
//!   solution_<k>             attrs: label, value (when present), eig_re, eig_im (when present)
//!     fields/<name>          f64, `n_dofs * components` values interleaved by dof; attrs: unit, components (JSON), name
//!     cell_fields/<name>/re, im   f64 per cell node; attrs on the group: unit, nloc
//! ```
//!
//! Large datasets are chunked and deflate-compressed (level 3); a result written and read back is bit-identical.

use crate::project::Project;
use femcore::mesh::SimplexMesh;
use hdf5_metno as h5;
use h5::types::VarLenUnicode;
use nyx_physics::{CellField, Field, Scalar, Solution, SolveResult};
use serde_json::{json, Value};
use std::path::Path;

pub(crate) fn e<T: std::fmt::Display>(what: &str) -> impl Fn(T) -> String + '_ {
    move |x| format!("HDF5 {what}: {x}")
}

fn text(s: &str) -> Result<VarLenUnicode, String> {
    s.parse::<VarLenUnicode>().map_err(|_| "text could not be stored (it contains a NUL character)".to_string())
}

fn attr<T: h5::H5Type>(loc: &h5::Location, name: &str, v: &T) -> Result<(), String> {
    loc.new_attr::<T>().create(name).map_err(e("attribute"))?.write_scalar(v).map_err(e("attribute"))
}

pub(crate) fn attr_text(loc: &h5::Location, name: &str, s: &str) -> Result<(), String> {
    attr(loc, name, &text(s)?)
}

/// A dataset of shape `[rows, cols]` (`[rows]` when `cols` is 0) holding `data` row by row, chunked and compressed when it is large.
pub(crate) fn put<T: h5::H5Type>(g: &h5::Group, name: &str, data: &[T], cols: usize) -> Result<h5::Dataset, String> {
    let rows = if cols == 0 { data.len() } else { data.len() / cols };
    let big = data.len() >= 4096 && rows > 0;
    let ds = if cols == 0 {
        let b = g.new_dataset::<T>().shape([rows]);
        if big { b.chunk([rows.min(32768)]).deflate(3).create(name) } else { b.create(name) }
    } else {
        let b = g.new_dataset::<T>().shape([rows, cols]);
        if big { b.chunk([rows.min((32768 / cols).max(1)), cols]).deflate(3).create(name) } else { b.create(name) }
    }
    .map_err(e(name))?;
    if !data.is_empty() {
        ds.write_raw(data).map_err(e(name))?;
    }
    Ok(ds)
}

fn put_text(g: &h5::Group, name: &str, s: &str) -> Result<(), String> {
    put(g, name, s.as_bytes(), 0).map(|_| ())
}

/// Name usable as an HDF5 link: `/` is the path separator, everything else is allowed.
pub(crate) fn link_name(name: &str) -> String {
    let n: String = name.chars().map(|c| if c == '/' { '_' } else { c }).collect();
    if n.is_empty() || n == "." { format!("_{n}") } else { n }
}

/// Writes `results` (physics code and result) solved on `mesh` for `project` into a new file at `path`.
pub fn write_results(path: &Path, project: &Project, mesh: &SimplexMesh, results: &[(String, &SolveResult)]) -> Result<(), String> {
    let file = h5::File::create(path).map_err(e("create file"))?;
    attr_text(&file, "format", "nyx-results")?;
    attr(&file, "version", &3i32)?;
    put_text(&file, "project", &project.to_json())?;
    let d = mesh.dim;
    let m = file.create_group("mesh").map_err(e("group"))?;
    put(&m, "points", &mesh.points, d)?;
    put(&m, "cells", &mesh.cells, d + 1)?;
    put(&m, "cell_labels", mesh.cell_labels.as_slice(), 0)?;
    put(&m, "facets", &mesh.facets, d)?;
    put(&m, "facet_labels", mesh.facet_labels.as_slice(), 0)?;
    let root = file.create_group("results").map_err(e("group"))?;
    for (code, r) in results {
        let g = root.create_group(&link_name(code)).map_err(e("group"))?;
        attr_text(&g, "physics", code)?;
        attr_text(&g, "study", &r.study)?;
        attr(&g, "order", &(r.order as i64))?;
        attr(&g, "dim", &(r.dim as i64))?;
        attr(&g, "dofs", &(r.dofs as i64))?;
        attr(&g, "seconds", &r.seconds)?;
        attr(&g, "n_solutions", &(r.solutions.len() as i64))?;
        let scalars: Vec<Value> = r
            .scalars
            .iter()
            .map(|(n, v)| match v {
                Scalar::Num(x) => json!([n, x]),
                Scalar::Text(t) => json!([n, t]),
            })
            .collect();
        put_text(&g, "scalars", &Value::Array(scalars).to_string())?;
        put_text(&g, "log", &serde_json::to_string(&r.log).map_err(|x| x.to_string())?)?;
        for (k, s) in r.solutions.iter().enumerate() {
            let sg = g.create_group(&format!("solution_{k}")).map_err(e("group"))?;
            attr_text(&sg, "label", &s.label)?;
            if let Some(v) = s.value {
                attr(&sg, "value", &v)?;
            }
            if let Some([re, im]) = s.eigenvalue {
                attr(&sg, "eig_re", &re)?;
                attr(&sg, "eig_im", &im)?;
            }
            let fg = sg.create_group("fields").map_err(e("group"))?;
            let order: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
            attr_text(&fg, "order", &serde_json::to_string(&order).map_err(|x| x.to_string())?)?;
            for (i, f) in s.fields.iter().enumerate() {
                let ds = put(&fg, &format!("{i:02}_{}", link_name(&f.name)), &f.values, 0)?;
                attr_text(&ds, "name", &f.name)?;
                attr_text(&ds, "unit", &f.unit)?;
                attr_text(&ds, "components", &serde_json::to_string(&f.components).map_err(|x| x.to_string())?)?;
            }
            if !s.cell_fields.is_empty() {
                let cg = sg.create_group("cell_fields").map_err(e("group"))?;
                for (i, f) in s.cell_fields.iter().enumerate() {
                    let fgp = cg.create_group(&format!("{i:02}_{}", link_name(&f.name))).map_err(e("group"))?;
                    attr_text(&fgp, "name", &f.name)?;
                    attr_text(&fgp, "unit", &f.unit)?;
                    attr(&fgp, "nloc", &(f.nloc as i64))?;
                    put(&fgp, "re", &f.re, 0)?;
                    if !f.im.is_empty() {
                        put(&fgp, "im", &f.im, 0)?;
                    }
                }
            }
        }
        crate::xdmf::write_viz(&file, code, mesh, r)?;
    }
    Ok(())
}

fn read_text(g: &h5::Group, name: &str) -> Result<String, String> {
    let bytes = g.dataset(name).map_err(e(name))?.read_raw::<u8>().map_err(e(name))?;
    String::from_utf8(bytes).map_err(|x| format!("HDF5 {name}: {x}"))
}

fn attr_f64(g: &h5::Location, name: &str) -> Option<f64> {
    g.attr(name).ok().and_then(|a| a.read_scalar::<f64>().ok())
}

fn attr_i64(g: &h5::Location, name: &str) -> usize {
    g.attr(name).ok().and_then(|a| a.read_scalar::<i64>().ok()).unwrap_or(0) as usize
}

fn attr_str(g: &h5::Location, name: &str) -> Option<String> {
    g.attr(name).ok().and_then(|a| a.read_scalar::<VarLenUnicode>().ok()).map(|v| v.as_str().to_string())
}

/// What a results file holds.
pub struct Stored {
    pub mesh: SimplexMesh,
    /// the project as JSON text
    pub project: String,
    /// physics code and result, in the order they were written
    pub results: Vec<(String, SolveResult)>,
}

/// Reads a results file back.
pub fn read_results(path: &Path) -> Result<Stored, String> {
    let file = h5::File::open(path).map_err(e("open file"))?;
    if attr_str(&file, "format").as_deref() != Some("nyx-results") {
        return Err("this HDF5 file was not written by nyx (no format attribute)".into());
    }
    if attr_i64(&file, "version") != 3 {
        return Err("this results file has another version than this nyx reads (3)".into());
    }
    let project = read_text(&file, "project")?;
    let m = file.group("mesh").map_err(e("mesh group"))?;
    let pts = m.dataset("points").map_err(e("points"))?;
    let dim = pts.shape().get(1).copied().ok_or("bad points shape")?;
    let mesh = SimplexMesh::new(
        dim,
        pts.read_raw::<f64>().map_err(e("points"))?,
        m.dataset("cells").map_err(e("cells"))?.read_raw::<u32>().map_err(e("cells"))?,
        m.dataset("cell_labels").map_err(e("cell labels"))?.read_raw::<i32>().map_err(e("cell labels"))?,
        m.dataset("facets").map_err(e("facets"))?.read_raw::<u32>().map_err(e("facets"))?,
        m.dataset("facet_labels").map_err(e("facet labels"))?.read_raw::<i32>().map_err(e("facet labels"))?,
    );
    let root = file.group("results").map_err(e("results group"))?;
    let mut codes = root.member_names().map_err(e("results"))?;
    codes.sort();
    let mut results = vec![];
    for name in codes {
        let g = root.group(&name).map_err(e("physics group"))?;
        let mut r = SolveResult {
            physics: attr_str(&g, "physics").unwrap_or_else(|| name.clone()),
            study: attr_str(&g, "study").unwrap_or_default(),
            order: attr_i64(&g, "order") as u8,
            dim: attr_i64(&g, "dim") as u8,
            dofs: attr_i64(&g, "dofs"),
            seconds: attr_f64(&g, "seconds").unwrap_or(0.0),
            ..Default::default()
        };
        if let Ok(Value::Array(items)) = serde_json::from_str::<Value>(&read_text(&g, "scalars")?) {
            for it in items {
                if let (Some(n), Some(v)) = (it.get(0).and_then(Value::as_str), it.get(1)) {
                    r.scalars.push((n.to_string(), match v.as_f64() {
                        Some(x) => Scalar::Num(x),
                        None => Scalar::Text(v.as_str().unwrap_or("").to_string()),
                    }));
                }
            }
        }
        r.log = serde_json::from_str(&read_text(&g, "log")?).unwrap_or_default();
        for k in 0..attr_i64(&g, "n_solutions") {
            let sg = g.group(&format!("solution_{k}")).map_err(e("solution"))?;
            let mut s = Solution { label: attr_str(&sg, "label").unwrap_or_default(), value: attr_f64(&sg, "value"), ..Default::default() };
            if let (Some(re), Some(im)) = (attr_f64(&sg, "eig_re"), attr_f64(&sg, "eig_im")) {
                s.eigenvalue = Some([re, im]);
            }
            let fg = sg.group("fields").map_err(e("fields"))?;
            let mut names = fg.member_names().map_err(e("fields"))?;
            names.sort();
            for n in names {
                let ds = fg.dataset(&n).map_err(e(&n))?;
                s.fields.push(Field {
                    name: attr_str(&ds, "name").unwrap_or_else(|| n.clone()),
                    unit: attr_str(&ds, "unit").unwrap_or_default(),
                    components: attr_str(&ds, "components").and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default(),
                    values: ds.read_raw::<f64>().map_err(e(&n))?,
                });
            }
            if let Ok(cg) = sg.group("cell_fields") {
                let mut names = cg.member_names().map_err(e("cell fields"))?;
                names.sort();
                for n in names {
                    let fgp = cg.group(&n).map_err(e(&n))?;
                    s.cell_fields.push(CellField {
                        name: attr_str(&fgp, "name").unwrap_or_else(|| n.clone()),
                        unit: attr_str(&fgp, "unit").unwrap_or_default(),
                        nloc: attr_i64(&fgp, "nloc"),
                        re: fgp.dataset("re").map_err(e("re"))?.read_raw::<f64>().map_err(e("re"))?,
                        im: match fgp.dataset("im") {
                            Ok(d) => d.read_raw::<f64>().map_err(e("im"))?,
                            Err(_) => vec![],
                        },
                    });
                }
            }
            r.solutions.push(s);
        }
        results.push((name, r));
    }
    Ok(Stored { mesh, project, results })
}

/// The size of a results file in bytes and how it divides: for the Share screen and the tests that compare with JSON.
pub fn file_bytes(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}
