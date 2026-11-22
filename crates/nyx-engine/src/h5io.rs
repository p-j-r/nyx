// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Results in HDF5 files (cargo feature `hdf5`): the project, the mesh and every solution of every solved physics with all its fields and
//! cell fields, written compressed and read back to 0 ULP (ADR 62, ADR 102) with pure-Rust `hdf5-pure` and `zlib-rs` (no C HDF5 library,
//! `HDF5_DIR`, or runtime DLL). This is the efficient store of the data a solve produces; the JSON files of a project and a run record stay
//! small.
//!
//! Layout, version 4 (numbers as stored by HDF5, little-endian; version 3 files, which lack `points`, `cells` and `cell_order` of a physics,
//! are still read):
//!
//! ```text
//! /                          attrs: format = "nyx-results", version = 4 (i32), writer, processes
//! /project                   the project as JSON text (u8 dataset)
//! /mesh/points               [n_points, dim] f64         /mesh/cells        [n_cells, dim+1] u32
//! /mesh/cell_labels          [n_cells] i32               /mesh/facets       [n_facets, dim] u32
//! /mesh/facet_labels         [n_facets] i32
//! /mesh/edge_node_pairs      [n_curved, 2] u32           /mesh/edge_node_points  [n_curved, 3] f64   (optional: curved edges, ADR 82; absent = straight)
//! /results/<physics>         attrs: physics, study, order, dim, dofs, seconds, n_solutions, cell_order, field_layout
//!   points                   [n_dofs, dim] f64         coordinates of the nodes (the dofs): vertices first, then edge midpoints
//!   cells                    [n_cells, nloc] u32       the node (dof) of every cell node, in the cell node order of `cell_order`
//!   scalars, log             JSON text (u8 datasets)
//!   eigenvalues              [n_solutions, 2] f64      [eig_re, eig_im] (or [value, 0.0]) per solution
//!   cell_fields_packed       [M, C, n_cells, nloc, R] f64   (L2 Shuffle+LZF chunked per mode, or L3 contiguous; R=1 real-only, R=2 [re, im])
//!   cell_fields_meta         JSON text (u8 dataset) of Vec<PackedComponentMeta> (name, unit, nloc, has_im)
//!   solution_<k>             attrs: label, value (when present), eig_re, eig_im (when present)
//!     fields/<name>          f64, `n_dofs * components` values interleaved by dof; attrs: unit, components (JSON), name
//!     cell_fields/<name>/re, im   f64 per cell node (legacy L1 / heterogeneous fallback); attrs on the group: name, unit, nloc
//! ```
//!
//! **One node order.** `fields/<name>` is indexed by node (dof), as `points`. `cell_fields` (`cell_fields_packed[k, c, cell, a, 0..2]` in
//! `L2`/`L3`, or `cell_fields/<name>/re[cell * nloc + a]` in `L1`) holds the value at node `cells[cell, a]` as seen from `cell`, so any
//! reader gets the nodal complex field by `nodal[cells[cell, a]] = cell_field[cell, a]` (for a continuous field every cell that holds the
//! node gives the same value). The cell node order is the one of nyx (`cell_order` says it): a P1 cell lists its vertices; a P2 triangle
//! `v0 v1 v2 e01 e02 e12`; a P2 tetrahedron `v0 v1 v2 v3 e01 e02 e03 e12 e13 e23`, where `eij` is the midpoint of the edge between vertices
//! `i` and `j`. `/mesh/cells` lists the vertices of a cell (the first `dim + 1` entries of its row of `cells`). There is no second copy of
//! the nodes or cells: the XDMF file next to the `.h5` (`xdmf.rs`) is a small XML index that names `points`, `cells` and the stored fields
//! as they are, for nyx and its tools only (not for ParaView or other viewers).

use crate::project::Project;
use femcore::mesh::SimplexMesh;
use hdf5_pure::{
    make_f64_type, make_i32_type, make_u32_type, make_u8_type, AttrValue, Dataset, DatasetBuilder, DeflateLevel, File,
    FileBuilder, Group,
};
use nyx_physics::{CellField, Field, Scalar, Solution, SolveResult};
pub use nyx_storage::{FieldStorageLayout, PackedComponentMeta, ReadBudgetConfig};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::path::Path;

pub(crate) fn e<T: std::fmt::Display>(what: &str) -> impl Fn(T) -> String + '_ {
    move |x| format!("HDF5 {what}: {x}")
}

fn text_attr(s: &str) -> Result<AttrValue, String> {
    if s.contains('\0') {
        return Err("text could not be stored (it contains a NUL character)".into());
    }
    Ok(AttrValue::VarLenString(s.to_string()))
}

fn chunk_dims(rows: usize, cols: usize) -> Vec<u64> {
    match 32768usize.checked_div(cols) {
        None => vec![rows.min(32768) as u64],
        Some(q) => vec![rows.min(q.max(1)) as u64, cols as u64],
    }
}

fn apply_numeric_filter(ds: &mut DatasetBuilder, len: usize, rows: usize, cols: usize, layout: FieldStorageLayout) {
    if rows == 0 {
        return;
    }
    match layout {
        FieldStorageLayout::L1LegacyGzip3 if len >= 4096 => {
            let chunks = chunk_dims(rows, cols);
            ds.with_chunks(&chunks).with_deflate(DeflateLevel::Three);
        }
        FieldStorageLayout::L2ShuffleLzfPacked if len >= 128 => {
            let chunks = chunk_dims(rows, cols);
            ds.with_chunks(&chunks).with_shuffle().with_deflate(DeflateLevel::Three);
        }
        _ => {}
    }
}

fn fill_f64<'a>(
    ds: &'a mut DatasetBuilder,
    data: &[f64],
    cols: usize,
    layout: FieldStorageLayout,
) -> &'a mut DatasetBuilder {
    if data.is_empty() {
        ds.with_dtype(make_f64_type());
        if cols == 0 {
            ds.with_shape(&[0]);
        } else {
            ds.with_shape(&[0, cols as u64]);
        }
        return ds;
    }
    let rows = data.len().checked_div(cols).unwrap_or(data.len());
    ds.with_f64_data(data);
    if cols == 0 {
        ds.with_shape(&[rows as u64]);
    } else {
        ds.with_shape(&[rows as u64, cols as u64]);
    }
    apply_numeric_filter(ds, data.len(), rows, cols, layout);
    ds
}

fn fill_u32<'a>(
    ds: &'a mut DatasetBuilder,
    data: &[u32],
    cols: usize,
    layout: FieldStorageLayout,
) -> &'a mut DatasetBuilder {
    if data.is_empty() {
        ds.with_dtype(make_u32_type());
        if cols == 0 {
            ds.with_shape(&[0]);
        } else {
            ds.with_shape(&[0, cols as u64]);
        }
        return ds;
    }
    let rows = data.len().checked_div(cols).unwrap_or(data.len());
    ds.with_u32_data(data);
    if cols == 0 {
        ds.with_shape(&[rows as u64]);
    } else {
        ds.with_shape(&[rows as u64, cols as u64]);
    }
    apply_numeric_filter(ds, data.len(), rows, cols, layout);
    ds
}

fn fill_i32<'a>(
    ds: &'a mut DatasetBuilder,
    data: &[i32],
    cols: usize,
    layout: FieldStorageLayout,
) -> &'a mut DatasetBuilder {
    if data.is_empty() {
        ds.with_dtype(make_i32_type());
        if cols == 0 {
            ds.with_shape(&[0]);
        } else {
            ds.with_shape(&[0, cols as u64]);
        }
        return ds;
    }
    let rows = data.len().checked_div(cols).unwrap_or(data.len());
    ds.with_i32_data(data);
    if cols == 0 {
        ds.with_shape(&[rows as u64]);
    } else {
        ds.with_shape(&[rows as u64, cols as u64]);
    }
    apply_numeric_filter(ds, data.len(), rows, cols, layout);
    ds
}

fn fill_u8<'a>(ds: &'a mut DatasetBuilder, data: &[u8], layout: FieldStorageLayout) -> &'a mut DatasetBuilder {
    if data.is_empty() {
        return ds.with_dtype(make_u8_type()).with_shape(&[0]);
    }
    ds.with_u8_data(data).with_shape(&[data.len() as u64]);
    if data.len() >= 4096 && layout != FieldStorageLayout::L3ContiguousPacked {
        ds.with_chunks(&[data.len().min(32768) as u64]).with_deflate(DeflateLevel::Three);
    }
    ds
}

/// Name usable as an HDF5 link: `/` is the path separator, everything else is allowed.
pub fn link_name(name: &str) -> String {
    let n: String = name.chars().map(|c| if c == '/' { '_' } else { c }).collect();
    if n.is_empty() || n == "." {
        format!("_{n}")
    } else {
        n
    }
}

/// How a results file is written: with the sequential pure-Rust HDF5 writer, or (stub) parallel MPI-IO.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H5Mode {
    Serial,
    Parallel,
}

/// True if this build can write with the mode: serial always, parallel never (`hdf5-pure` has no MPI-IO driver; `phdf5` is a zero-C stub).
pub fn mode_available(mode: H5Mode) -> bool {
    match mode {
        H5Mode::Serial => true,
        H5Mode::Parallel => false,
    }
}

const NO_PARALLEL: &str = "parallel HDF5 (MPI-IO) is not available in this build: pure-Rust `hdf5-pure` has no MPI-IO driver (`phdf5` is a zero-C stub without MPI); use the serial output";

/// Starts what the mode needs before a solve. Returns `(0, 1)` for `Serial` and an explanatory error for `Parallel`.
pub fn start_runtime(mode: H5Mode) -> Result<(usize, usize), String> {
    match mode {
        H5Mode::Serial => Ok((0, 1)),
        H5Mode::Parallel => Err(NO_PARALLEL.into()),
    }
}

/// Process `rank` of `size` that write one results file together. Each process owns a contiguous block of rows of every dataset; a serial
/// run is the single process 0 of 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Partition {
    pub rank: usize,
    pub size: usize,
}

impl Partition {
    pub const SOLO: Partition = Partition { rank: 0, size: 1 };

    /// The rows `[start, end)` of `n` this process writes; the blocks of all processes tile `0..n` exactly, with sizes differing by at most one.
    pub fn rows(&self, n: usize) -> std::ops::Range<usize> {
        let size = self.size.max(1);
        let (base, rem) = (n / size, n % size);
        let start = self.rank * base + self.rank.min(rem);
        start..start + base + usize::from(self.rank < rem)
    }
}

fn packable_cell_fields(r: &SolveResult, n_cells: usize) -> Option<(usize, Vec<PackedComponentMeta>)> {
    if r.solutions.is_empty() || n_cells == 0 {
        return None;
    }
    let first = &r.solutions[0].cell_fields;
    if first.is_empty() {
        return None;
    }
    let nloc = first[0].nloc;
    if nloc == 0 {
        return None;
    }
    let field_len = n_cells.checked_mul(nloc)?;
    let meta: Vec<PackedComponentMeta> = first
        .iter()
        .map(|f| PackedComponentMeta {
            name: f.name.clone(),
            unit: f.unit.clone(),
            nloc: f.nloc,
            has_im: !f.im.is_empty(),
        })
        .collect();
    for s in &r.solutions {
        if s.cell_fields.len() != meta.len() {
            return None;
        }
        for (f, m) in s.cell_fields.iter().zip(&meta) {
            if f.name != m.name || f.unit != m.unit || f.nloc != nloc || f.re.len() != field_len {
                return None;
            }
            if m.has_im {
                if f.im.len() != field_len {
                    return None;
                }
            } else if !f.im.is_empty() {
                return None;
            }
        }
    }
    Some((nloc, meta))
}

fn build_packed_cell_fields(r: &SolveResult, meta: &[PackedComponentMeta], field_len: usize) -> (usize, Vec<f64>) {
    let any_im = meta.iter().any(|m| m.has_im);
    let ri = if any_im { 2 } else { 1 };
    let mut packed = Vec::with_capacity(r.solutions.len() * meta.len() * field_len * ri);
    if any_im {
        for s in &r.solutions {
            for (f, m) in s.cell_fields.iter().zip(meta) {
                if m.has_im {
                    for (&re, &im) in f.re.iter().zip(&f.im) {
                        packed.push(re);
                        packed.push(im);
                    }
                } else {
                    for &re in &f.re {
                        packed.push(re);
                        packed.push(0.0);
                    }
                }
            }
        }
    } else {
        for s in &r.solutions {
            for f in &s.cell_fields {
                packed.extend_from_slice(&f.re);
            }
        }
    }
    (ri, packed)
}

fn serialize_scalars(scalars: &[(String, Scalar)]) -> String {
    let items: Vec<Value> = scalars
        .iter()
        .map(|(n, v)| match v {
            Scalar::Num(x) => json!([n, x]),
            Scalar::Text(t) => json!([n, t]),
        })
        .collect();
    Value::Array(items).to_string()
}

fn parse_scalars(text: &str) -> Vec<(String, Scalar)> {
    let mut out = Vec::new();
    if let Ok(Value::Array(items)) = serde_json::from_str::<Value>(text) {
        for it in items {
            if let (Some(n), Some(v)) = (it.get(0).and_then(Value::as_str), it.get(1)) {
                out.push((
                    n.to_string(),
                    match v.as_f64() {
                        Some(x) => Scalar::Num(x),
                        None => Scalar::Text(v.as_str().unwrap_or("").to_string()),
                    },
                ));
            }
        }
    }
    out
}

/// Writes `results` (physics code and result) solved on `mesh` for `project` into a new file at `path`, serially using `L2ShuffleLzfPacked`.
pub fn write_results(
    path: &Path,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<(), String> {
    write_results_mode(path, H5Mode::Serial, project, mesh, results).map(|_| ())
}

/// The same in the chosen mode (`Serial` uses `L2ShuffleLzfPacked`; `Parallel` returns an error because `hdf5-pure` has no MPI-IO driver).
/// Returns whether this process is the one that writes the XDMF file (`true` on `Serial`).
pub fn write_results_mode(
    path: &Path,
    mode: H5Mode,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<bool, String> {
    match mode {
        H5Mode::Serial => {
            write_results_with_layout(path, project, mesh, results, FieldStorageLayout::L2ShuffleLzfPacked)?;
            Ok(true)
        }
        H5Mode::Parallel => Err(NO_PARALLEL.into()),
    }
}

/// Writes `results` solved on `mesh` for `project` into `path` using the requested [`FieldStorageLayout`].
pub fn write_results_with_layout(
    path: &Path,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
    layout: FieldStorageLayout,
) -> Result<(), String> {
    write_results_with_options(path, project, mesh, results, layout, true)
}

/// Writes `results` solved on `mesh` for `project` into `path` using `layout`, optionally omitting all mesh and field arrays when
/// `save_fields` is `false` (preserving `project`, `scalars`, `log`, `eigenvalues`, and per-solution metadata with `@fields_saved = 0`).
pub fn write_results_with_options(
    path: &Path,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
    layout: FieldStorageLayout,
    save_fields: bool,
) -> Result<(), String> {
    let mut fb = FileBuilder::new();
    fb.set_attr("writer", text_attr("serial")?);
    fb.set_attr("processes", AttrValue::I64(1));
    fb.set_attr("format", text_attr("nyx-results")?);
    fb.set_attr("version", AttrValue::I32(4));
    if !save_fields {
        fb.set_attr("fields_saved", AttrValue::I64(0));
        fb.set_attr("mesh_dim", AttrValue::I64(mesh.dim as i64));
        fb.set_attr("mesh_points", AttrValue::I64(mesh.n_points() as i64));
        fb.set_attr("mesh_cells", AttrValue::I64(mesh.n_cells() as i64));
    }
    let proj_json = project.to_json();
    fill_u8(fb.create_dataset("project"), proj_json.as_bytes(), layout);

    let d = mesh.dim;
    if save_fields {
        let mut m = fb.create_group("mesh");
        fill_f64(m.create_dataset("points"), &mesh.points, d, layout);
        fill_u32(m.create_dataset("cells"), &mesh.cells, d + 1, layout);
        fill_i32(m.create_dataset("cell_labels"), mesh.cell_labels.as_slice(), 0, layout);
        fill_u32(m.create_dataset("facets"), &mesh.facets, d, layout);
        fill_i32(m.create_dataset("facet_labels"), mesh.facet_labels.as_slice(), 0, layout);
        // curved edges (ADR 82): written only when there are some, so a straight mesh gives the file it always gave
        let curved = mesh.edge_node_list();
        if !curved.is_empty() {
            let pairs: Vec<u32> = curved.iter().flat_map(|e| e.0).collect();
            let pos: Vec<f64> = curved.iter().flat_map(|e| e.1).collect();
            fill_u32(m.create_dataset("edge_node_pairs"), &pairs, 2, layout);
            fill_f64(m.create_dataset("edge_node_points"), &pos, 3, layout);
        }
        let qp = mesh.qp_edge_list();
        if !qp.is_empty() {
            let qp_pairs: Vec<u32> = qp.iter().flat_map(|e| *e).collect();
            fill_u32(m.create_dataset("qp_edge_pairs"), &qp_pairs, 2, layout);
        }
        fb.add_group(m.finish());
    }

    let mut root = fb.create_group("results");
    let n_cells = mesh.n_cells();
    for (code, r) in results {
        let mut g = root.create_group(&link_name(code));
        g.set_attr("physics", text_attr(code)?);
        g.set_attr("study", text_attr(&r.study)?);
        g.set_attr("order", AttrValue::I64(r.order as i64));
        g.set_attr("dim", AttrValue::I64(r.dim as i64));
        g.set_attr("dofs", AttrValue::I64(r.dofs as i64));
        g.set_attr("seconds", AttrValue::F64(r.seconds));
        g.set_attr("n_solutions", AttrValue::I64(r.solutions.len() as i64));
        if !save_fields {
            g.set_attr("fields_saved", AttrValue::I64(0));
        }

        // the nodes of this physics and the node of every cell node, in nyx's own cell node order (the order of every cell field)
        let ord = if r.order == 2 { 2 } else { 1 };
        if save_fields {
            let nloc_mesh = mesh.nloc(ord);
            let mut conn: Vec<u32> = Vec::with_capacity(n_cells * nloc_mesh);
            let mut dofs = Vec::new();
            for c in 0..n_cells {
                mesh.cell_dofs(c, ord, &mut dofs);
                conn.extend(dofs.iter().take(nloc_mesh).map(|&idx| idx as u32));
            }
            let dof_pts = mesh.dof_coordinates(ord);
            fill_f64(g.create_dataset("points"), &dof_pts, mesh.dim, layout);
            fill_u32(g.create_dataset("cells"), &conn, nloc_mesh, layout);
        }
        g.set_attr("cell_order", text_attr(cell_order_text(mesh.dim, ord))?);

        let scalars_text = serialize_scalars(&r.scalars);
        fill_u8(g.create_dataset("scalars"), scalars_text.as_bytes(), layout);
        let log_text = serde_json::to_string(&r.log).map_err(|x| x.to_string())?;
        fill_u8(g.create_dataset("log"), log_text.as_bytes(), layout);

        let m_sol = r.solutions.len();
        if m_sol == 0 {
            g.create_dataset("eigenvalues").with_dtype(make_f64_type()).with_shape(&[0, 2]);
        } else {
            let mut eig_flat = Vec::with_capacity(m_sol * 2);
            for s in &r.solutions {
                let [re, im] = s.eigenvalue.unwrap_or([s.value.unwrap_or(0.0), 0.0]);
                eig_flat.push(re);
                eig_flat.push(im);
            }
            fill_f64(g.create_dataset("eigenvalues"), &eig_flat, 2, layout);
        }

        let packed_info = if save_fields && layout != FieldStorageLayout::L1LegacyGzip3 {
            packable_cell_fields(r, n_cells)
        } else {
            None
        };
        let effective_layout =
            if !save_fields || packed_info.is_some() || r.solutions.iter().all(|s| s.cell_fields.is_empty()) {
                layout
            } else {
                FieldStorageLayout::L1LegacyGzip3
            };
        g.set_attr("field_layout", text_attr(effective_layout.as_str())?);

        if save_fields {
            if let Some((nloc_cf, ref meta)) = packed_info {
                let field_len = n_cells * nloc_cf;
                let (ri, packed) = build_packed_cell_fields(r, meta, field_len);
                let c_comp = meta.len() as u64;
                let ri_u64 = ri as u64;
                let ds_packed = g.create_dataset("cell_fields_packed").with_f64_data(&packed).with_shape(&[
                    m_sol as u64,
                    c_comp,
                    n_cells as u64,
                    nloc_cf as u64,
                    ri_u64,
                ]);
                if layout == FieldStorageLayout::L2ShuffleLzfPacked {
                    ds_packed
                        .with_chunks(&[1, c_comp, n_cells as u64, nloc_cf as u64, ri_u64])
                        .with_shuffle()
                        .with_lzf();
                }
                let meta_bytes = serde_json::to_vec(meta).map_err(|x| x.to_string())?;
                fill_u8(g.create_dataset("cell_fields_meta"), &meta_bytes, layout);
            }
        }

        for (k, s) in r.solutions.iter().enumerate() {
            let mut sg = g.create_group(&format!("solution_{k}"));
            sg.set_attr("label", text_attr(&s.label)?);
            if let Some(v) = s.value {
                sg.set_attr("value", AttrValue::F64(v));
            }
            if let Some([re, im]) = s.eigenvalue {
                sg.set_attr("eig_re", AttrValue::F64(re));
                sg.set_attr("eig_im", AttrValue::F64(im));
            }

            if save_fields {
                let mut fg = sg.create_group("fields");
                let order: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
                fg.set_attr("order", text_attr(&serde_json::to_string(&order).map_err(|x| x.to_string())?)?);
                for (i, f) in s.fields.iter().enumerate() {
                    let comp_json = serde_json::to_string(&f.components).map_err(|x| x.to_string())?;
                    let ds =
                        fill_f64(fg.create_dataset(&format!("{i:02}_{}", link_name(&f.name))), &f.values, 0, layout);
                    ds.set_attr("name", text_attr(&f.name)?)
                        .set_attr("unit", text_attr(&f.unit)?)
                        .set_attr("components", text_attr(&comp_json)?);
                }
                sg.add_group(fg.finish());

                if packed_info.is_none() && !s.cell_fields.is_empty() {
                    let mut cg = sg.create_group("cell_fields");
                    for (i, f) in s.cell_fields.iter().enumerate() {
                        let mut fgp = cg.create_group(&format!("{i:02}_{}", link_name(&f.name)));
                        fgp.set_attr("name", text_attr(&f.name)?);
                        fgp.set_attr("unit", text_attr(&f.unit)?);
                        fgp.set_attr("nloc", AttrValue::I64(f.nloc as i64));
                        fill_f64(fgp.create_dataset("re"), &f.re, 0, FieldStorageLayout::L1LegacyGzip3);
                        if !f.im.is_empty() {
                            fill_f64(fgp.create_dataset("im"), &f.im, 0, FieldStorageLayout::L1LegacyGzip3);
                        }
                        cg.add_group(fgp.finish());
                    }
                    sg.add_group(cg.finish());
                }
            }

            g.add_group(sg.finish());
        }
        root.add_group(g.finish());
    }
    fb.add_group(root.finish());

    fb.write(path).map_err(e("write file"))
}

/// The text of the `cell_order` attribute: the cell node order of the cell fields and of `cells`.
pub fn cell_order_text(dim: usize, order: usize) -> &'static str {
    match (dim, order) {
        (2, 1) => "v0 v1 v2",
        (3, 1) => "v0 v1 v2 v3",
        (2, _) => "v0 v1 v2 e01 e02 e12",
        _ => "v0 v1 v2 v3 e01 e02 e03 e12 e13 e23",
    }
}

fn check_budget(ds: &Dataset, elem_bytes: u64, budget: &ReadBudgetConfig, name: &str) -> Result<Vec<u64>, String> {
    let shape = ds.shape().map_err(e(name))?;
    let elems = shape
        .iter()
        .try_fold(1u64, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| format!("HDF5 {name}: dataset shape overflow ({shape:?})"))?;
    let bytes =
        elems.checked_mul(elem_bytes).ok_or_else(|| format!("HDF5 {name}: dataset byte size overflow ({shape:?})"))?;
    if bytes > budget.max_dataset_bytes {
        return Err(format!("HDF5 {name}: dataset size {bytes} B exceeds read budget {} B", budget.max_dataset_bytes));
    }
    if bytes > budget.max_total_result_bytes {
        return Err(format!(
            "HDF5 {name}: dataset size {bytes} B exceeds read budget {} B",
            budget.max_total_result_bytes
        ));
    }
    Ok(shape)
}

/// Enforces per-mode slice (`shape[1..] <= budget.max_dataset_bytes`) and total dataset
/// (`shape <= budget.max_total_result_bytes`) bounds on 5D `cell_fields_packed` read row-by-row.
fn check_packed_budget(
    ds: &Dataset,
    elem_bytes: u64,
    budget: &ReadBudgetConfig,
    name: &str,
) -> Result<Vec<u64>, String> {
    let shape = ds.shape().map_err(e(name))?;
    if shape.is_empty() {
        return Err(format!("HDF5 {name}: empty dataset shape"));
    }
    let slice_bytes = shape[1..]
        .iter()
        .try_fold(elem_bytes, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| format!("HDF5 {name}: per-mode slice byte size overflow ({shape:?})"))?;
    if slice_bytes > budget.max_dataset_bytes {
        return Err(format!(
            "HDF5 {name}: per-mode slice size {slice_bytes} B exceeds read budget {} B",
            budget.max_dataset_bytes
        ));
    }
    let total_bytes = shape
        .iter()
        .try_fold(elem_bytes, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| format!("HDF5 {name}: dataset byte size overflow ({shape:?})"))?;
    if total_bytes > budget.max_total_result_bytes {
        return Err(format!(
            "HDF5 {name}: total dataset size {total_bytes} B exceeds read budget {} B",
            budget.max_total_result_bytes
        ));
    }
    Ok(shape)
}

fn read_text(g: &Group, name: &str, budget: &ReadBudgetConfig) -> Result<String, String> {
    let ds = g.dataset(name).map_err(e(name))?;
    check_budget(&ds, 1, budget, name)?;
    let bytes = ds.read_u8().map_err(e(name))?;
    String::from_utf8(bytes).map_err(|x| format!("HDF5 {name}: {x}"))
}

fn attr_f64(attrs: &HashMap<String, AttrValue>, name: &str) -> Option<f64> {
    attrs.get(name).and_then(AttrValue::as_f64)
}

fn attr_i64(attrs: &HashMap<String, AttrValue>, name: &str) -> usize {
    attrs.get(name).and_then(|v| v.as_i64().or_else(|| v.as_u64().map(|u| u as i64))).unwrap_or(0).max(0) as usize
}

fn attr_bool(attrs: &HashMap<String, AttrValue>, name: &str) -> Option<bool> {
    let v = attrs.get(name)?;
    if let Some(i) = v.as_i64().or_else(|| v.as_u64().map(|u| u as i64)) {
        return Some(i != 0);
    }
    if let Some(s) = v.as_str() {
        return match s.trim().to_ascii_lowercase().as_str() {
            "false" | "0" | "no" => Some(false),
            "true" | "1" | "yes" => Some(true),
            _ => None,
        };
    }
    None
}

fn attr_str(attrs: &HashMap<String, AttrValue>, name: &str) -> Option<String> {
    let v = attrs.get(name)?;
    if let Some(s) = v.as_str() {
        return Some(s.to_string());
    }
    if let Some(arr) = v.as_strings() {
        if let Some(first) = arr.first() {
            return Some(first.clone());
        }
    }
    None
}

fn project_dim(project_json: &str) -> usize {
    serde_json::from_str::<Value>(project_json)
        .ok()
        .and_then(|v| v.get("geometry")?.get("dim")?.as_u64())
        .filter(|&d| d == 2 || d == 3)
        .unwrap_or(3) as usize
}

fn check_nyx_header(file: &File) -> Result<(), String> {
    let root_attrs = file.root().attrs().map_err(e("root attributes"))?;
    if attr_str(&root_attrs, "format").as_deref() != Some("nyx-results") {
        return Err("this HDF5 file was not written by nyx (no format attribute)".into());
    }
    if !matches!(attr_i64(&root_attrs, "version"), 3 | 4) {
        return Err("this results file has another version than this nyx reads (3 and 4)".into());
    }
    Ok(())
}

/// Summary of one solved physics in a results file, read without loading heavy field or mesh arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct PhysicsSummary {
    pub physics: String,
    pub study: String,
    pub order: u8,
    pub dim: u8,
    pub dofs: usize,
    pub seconds: f64,
    pub solutions: usize,
    pub scalars: Vec<(String, Scalar)>,
}

/// Lightweight metadata summary of a results file (`results_file_info`), read via streaming metadata only.
#[derive(Clone, Debug, PartialEq)]
pub struct ResultsFileInfo {
    pub project: String,
    pub dim: usize,
    pub points: usize,
    pub cells: usize,
    pub results: Vec<(String, PhysicsSummary)>,
}

/// Reads the header summary of a results file via `File::open_streaming` without loading mesh or field arrays.
pub fn read_results_summary(path: &Path) -> Result<ResultsFileInfo, String> {
    let budget = ReadBudgetConfig::from_env();
    let file = File::open_streaming(path).map_err(e("open file"))?;
    check_nyx_header(&file)?;
    let root_grp = file.root();
    let root_attrs = root_grp.attrs().map_err(e("root attributes"))?;
    let fields_saved = attr_bool(&root_attrs, "fields_saved").unwrap_or(true);
    let project = read_text(&root_grp, "project", &budget)?;

    let (dim, points, cells) = if fields_saved {
        let m = file.group("mesh").map_err(e("mesh group"))?;
        let pts_shape = m.dataset("points").map_err(e("points"))?.shape().map_err(e("points"))?;
        let cells_shape = m.dataset("cells").map_err(e("cells"))?.shape().map_err(e("cells"))?;
        let dim = pts_shape.get(1).copied().ok_or("bad points shape")? as usize;
        let points = pts_shape.first().copied().unwrap_or(0) as usize;
        let cells = cells_shape.first().copied().unwrap_or(0) as usize;
        (dim, points, cells)
    } else {
        let d = match attr_i64(&root_attrs, "mesh_dim") {
            2 | 3 => attr_i64(&root_attrs, "mesh_dim"),
            _ => project_dim(&project),
        };
        (d, attr_i64(&root_attrs, "mesh_points"), attr_i64(&root_attrs, "mesh_cells"))
    };

    let res_root = file.group("results").map_err(e("results group"))?;
    let mut codes = res_root.groups().map_err(e("results"))?;
    codes.sort();
    let mut results = Vec::with_capacity(codes.len());
    for name in codes {
        let g = res_root.group(&name).map_err(e("physics group"))?;
        let g_attrs = g.attrs().map_err(e("physics attributes"))?;
        let scalars = parse_scalars(&read_text(&g, "scalars", &budget)?);
        let summary = PhysicsSummary {
            physics: attr_str(&g_attrs, "physics").unwrap_or_else(|| name.clone()),
            study: attr_str(&g_attrs, "study").unwrap_or_default(),
            order: attr_i64(&g_attrs, "order") as u8,
            dim: attr_i64(&g_attrs, "dim") as u8,
            dofs: attr_i64(&g_attrs, "dofs"),
            seconds: attr_f64(&g_attrs, "seconds").unwrap_or(0.0),
            solutions: attr_i64(&g_attrs, "n_solutions"),
            scalars,
        };
        results.push((name, summary));
    }
    Ok(ResultsFileInfo { project, dim, points, cells, results })
}

/// What a results file holds.
pub struct Stored {
    pub mesh: SimplexMesh,
    /// the project as JSON text
    pub project: String,
    /// physics code and result, in the order they were written
    pub results: Vec<(String, SolveResult)>,
    /// false when the file was written with `save_fields == false` (mesh and solution field datasets omitted)
    pub fields_saved: bool,
}

/// Reads a results file back using `File::open_streaming` and [`ReadBudgetConfig::from_env()`].
pub fn read_results(path: &Path) -> Result<Stored, String> {
    read_results_with_budget(path, ReadBudgetConfig::from_env())
}

/// Reads a results file back using `File::open_streaming` and an explicit [`ReadBudgetConfig`].
pub fn read_results_with_budget(path: &Path, budget: ReadBudgetConfig) -> Result<Stored, String> {
    let file = File::open_streaming(path).map_err(e("open file"))?;
    check_nyx_header(&file)?;
    let root_grp = file.root();
    let root_attrs = root_grp.attrs().map_err(e("root attributes"))?;
    let fields_saved = attr_bool(&root_attrs, "fields_saved").unwrap_or(true);
    let project = read_text(&root_grp, "project", &budget)?;

    let mut total_read_bytes = project.len() as u64;
    let mut add_total = |bytes: u64, name: &str| -> Result<(), String> {
        total_read_bytes = total_read_bytes
            .checked_add(bytes)
            .ok_or_else(|| format!("HDF5 {name}: cumulative read byte count overflow"))?;
        if total_read_bytes > budget.max_total_result_bytes {
            return Err(format!(
                "HDF5 {name}: cumulative result size {total_read_bytes} B exceeds read budget {} B",
                budget.max_total_result_bytes
            ));
        }
        Ok(())
    };
    add_total(0, "project")?;

    let mesh = if fields_saved {
        let m = file.group("mesh").map_err(e("mesh group"))?;
        let pts = m.dataset("points").map_err(e("points"))?;
        let pts_shape = check_budget(&pts, 8, &budget, "points")?;
        let dim = pts_shape.get(1).copied().ok_or("bad points shape")? as usize;

        let ds_cells = m.dataset("cells").map_err(e("cells"))?;
        check_budget(&ds_cells, 4, &budget, "cells")?;
        let ds_cl = m.dataset("cell_labels").map_err(e("cell labels"))?;
        check_budget(&ds_cl, 4, &budget, "cell labels")?;
        let ds_facets = m.dataset("facets").map_err(e("facets"))?;
        check_budget(&ds_facets, 4, &budget, "facets")?;
        let ds_fl = m.dataset("facet_labels").map_err(e("facet labels"))?;
        check_budget(&ds_fl, 4, &budget, "facet labels")?;

        let pts_v = pts.read_f64().map_err(e("points"))?;
        add_total((pts_v.len() as u64) * 8, "points")?;
        let cells_v = ds_cells.read_u32().map_err(e("cells"))?;
        add_total((cells_v.len() as u64) * 4, "cells")?;
        let cl_v = ds_cl.read_i32().map_err(e("cell labels"))?;
        add_total((cl_v.len() as u64) * 4, "cell labels")?;
        let facets_v = ds_facets.read_u32().map_err(e("facets"))?;
        add_total((facets_v.len() as u64) * 4, "facets")?;
        let fl_v = ds_fl.read_i32().map_err(e("facet labels"))?;
        add_total((fl_v.len() as u64) * 4, "facet labels")?;

        let mut mesh = SimplexMesh::new(dim, pts_v, cells_v, cl_v, facets_v, fl_v);
        // optional curved edges: a file without them is a straight mesh
        if let (Ok(p), Ok(q)) = (m.dataset("edge_node_pairs"), m.dataset("edge_node_points")) {
            check_budget(&p, 4, &budget, "edge node pairs")?;
            check_budget(&q, 8, &budget, "edge node points")?;
            let pairs = p.read_u32().map_err(e("edge node pairs"))?;
            add_total((pairs.len() as u64) * 4, "edge node pairs")?;
            let pos = q.read_f64().map_err(e("edge node points"))?;
            add_total((pos.len() as u64) * 8, "edge node points")?;
            let nodes: Vec<([u32; 2], [f64; 3])> =
                pairs.as_chunks::<2>().0.iter().copied().zip(pos.as_chunks::<3>().0.iter().copied()).collect();
            mesh.set_edge_nodes(&nodes);
        }
        if let Ok(qp_ds) = m.dataset("qp_edge_pairs") {
            check_budget(&qp_ds, 4, &budget, "qp edge pairs")?;
            let qp_raw = qp_ds.read_u32().map_err(e("qp edge pairs"))?;
            add_total((qp_raw.len() as u64) * 4, "qp edge pairs")?;
            let qp_edges: Vec<[u32; 2]> = qp_raw.as_chunks::<2>().0.to_vec();
            mesh.set_qp_edges(&qp_edges);
        }
        mesh
    } else {
        let d = match attr_i64(&root_attrs, "mesh_dim") {
            2 | 3 => attr_i64(&root_attrs, "mesh_dim"),
            _ => project_dim(&project),
        };
        SimplexMesh::new(d, vec![], vec![], vec![], vec![], vec![])
    };

    let root = file.group("results").map_err(e("results group"))?;
    let mut codes = root.groups().map_err(e("results"))?;
    codes.sort();
    let mut results = vec![];
    for name in codes {
        let g = root.group(&name).map_err(e("physics group"))?;
        let g_attrs = g.attrs().map_err(e("physics attributes"))?;
        let g_fields_saved = attr_bool(&g_attrs, "fields_saved").unwrap_or(fields_saved);
        let n_solutions = attr_i64(&g_attrs, "n_solutions");
        let mut r = SolveResult {
            physics: attr_str(&g_attrs, "physics").unwrap_or_else(|| name.clone()),
            study: attr_str(&g_attrs, "study").unwrap_or_default(),
            order: attr_i64(&g_attrs, "order") as u8,
            dim: attr_i64(&g_attrs, "dim") as u8,
            dofs: attr_i64(&g_attrs, "dofs"),
            seconds: attr_f64(&g_attrs, "seconds").unwrap_or(0.0),
            scalars: parse_scalars(&read_text(&g, "scalars", &budget)?),
            log: serde_json::from_str(&read_text(&g, "log", &budget)?).unwrap_or_default(),
            solutions: Vec::with_capacity(n_solutions),
        };

        let packed_cf = if g_fields_saved {
            if let Ok(ds_packed) = g.dataset("cell_fields_packed") {
                let p_shape = check_packed_budget(&ds_packed, 8, &budget, "cell_fields_packed")?;
                if p_shape.len() != 5 || (p_shape[4] != 1 && p_shape[4] != 2) || p_shape[0] as usize != n_solutions {
                    return Err(format!("HDF5 cell_fields_packed: bad shape {p_shape:?} for {n_solutions} solutions"));
                }
                let packed_total = p_shape
                    .iter()
                    .try_fold(8u64, |acc, &d| acc.checked_mul(d))
                    .ok_or("HDF5 cell_fields_packed: byte size overflow")?;
                add_total(packed_total, "cell_fields_packed")?;
                let ri = p_shape[4] as usize;
                let meta_json = read_text(&g, "cell_fields_meta", &budget)?;
                let meta: Vec<PackedComponentMeta> =
                    serde_json::from_str(&meta_json).map_err(|x| format!("HDF5 cell_fields_meta: {x}"))?;
                if meta.len() != p_shape[1] as usize {
                    return Err(format!(
                        "HDF5 cell_fields_meta: {} components, packed shape has {}",
                        meta.len(),
                        p_shape[1]
                    ));
                }
                if ri == 1 && meta.iter().any(|m| m.has_im) {
                    return Err(
                        "HDF5 cell_fields_packed: shape has R=1 (real-only) but metadata marks has_im=true".into()
                    );
                }
                let field_len = (p_shape[2] as usize)
                    .checked_mul(p_shape[3] as usize)
                    .ok_or("HDF5 cell_fields_packed: field length overflow")?;
                Some((ds_packed, meta, field_len, ri))
            } else {
                None
            }
        } else {
            None
        };

        for k in 0..n_solutions {
            let sg = g.group(&format!("solution_{k}")).map_err(e("solution"))?;
            let sg_attrs = sg.attrs().map_err(e("solution attributes"))?;
            let mut s = Solution {
                label: attr_str(&sg_attrs, "label").unwrap_or_default(),
                value: attr_f64(&sg_attrs, "value"),
                ..Default::default()
            };
            if let (Some(re), Some(im)) = (attr_f64(&sg_attrs, "eig_re"), attr_f64(&sg_attrs, "eig_im")) {
                s.eigenvalue = Some([re, im]);
            }

            if g_fields_saved {
                let fg = sg.group("fields").map_err(e("fields"))?;
                let mut names = fg.datasets().map_err(e("fields"))?;
                names.sort();
                for n in names {
                    let ds = fg.dataset(&n).map_err(e(&n))?;
                    check_budget(&ds, 8, &budget, &n)?;
                    let ds_attrs = ds.attrs().map_err(e(&n))?;
                    let values = ds.read_f64().map_err(e(&n))?;
                    add_total((values.len() as u64) * 8, &n)?;
                    s.fields.push(Field {
                        name: attr_str(&ds_attrs, "name").unwrap_or_else(|| n.clone()),
                        unit: attr_str(&ds_attrs, "unit").unwrap_or_default(),
                        components: attr_str(&ds_attrs, "components")
                            .and_then(|t| serde_json::from_str(&t).ok())
                            .unwrap_or_default(),
                        values,
                    });
                }

                if let Some((ref ds_packed, ref meta, field_len, ri)) = packed_cf {
                    let row = ds_packed.read_f64_rows(k as u64, 1).map_err(e("cell_fields_packed"))?;
                    for (c_idx, cm) in meta.iter().enumerate() {
                        let (re, im) = if ri == 1 {
                            let comp = &row[c_idx * field_len..(c_idx + 1) * field_len];
                            (comp.to_vec(), Vec::new())
                        } else {
                            let comp = &row[c_idx * field_len * 2..(c_idx + 1) * field_len * 2];
                            let mut re = Vec::with_capacity(field_len);
                            let mut im = if cm.has_im { Vec::with_capacity(field_len) } else { Vec::new() };
                            if cm.has_im {
                                for &[re_v, im_v] in comp.as_chunks::<2>().0 {
                                    re.push(re_v);
                                    im.push(im_v);
                                }
                            } else {
                                for &[re_v, _] in comp.as_chunks::<2>().0 {
                                    re.push(re_v);
                                }
                            }
                            (re, im)
                        };
                        s.cell_fields.push(CellField {
                            name: cm.name.clone(),
                            unit: cm.unit.clone(),
                            nloc: cm.nloc,
                            re,
                            im,
                        });
                    }
                } else if let Ok(cg) = sg.group("cell_fields") {
                    let mut names = cg.groups().map_err(e("cell fields"))?;
                    names.sort();
                    for n in names {
                        let fgp = cg.group(&n).map_err(e(&n))?;
                        let fgp_attrs = fgp.attrs().map_err(e(&n))?;
                        let ds_re = fgp.dataset("re").map_err(e("re"))?;
                        check_budget(&ds_re, 8, &budget, "re")?;
                        let re = ds_re.read_f64().map_err(e("re"))?;
                        add_total((re.len() as u64) * 8, "re")?;
                        let im = match fgp.dataset("im") {
                            Ok(ds_im) => {
                                check_budget(&ds_im, 8, &budget, "im")?;
                                let im_v = ds_im.read_f64().map_err(e("im"))?;
                                add_total((im_v.len() as u64) * 8, "im")?;
                                im_v
                            }
                            Err(_) => vec![],
                        };
                        s.cell_fields.push(CellField {
                            name: attr_str(&fgp_attrs, "name").unwrap_or_else(|| n.clone()),
                            unit: attr_str(&fgp_attrs, "unit").unwrap_or_default(),
                            nloc: attr_i64(&fgp_attrs, "nloc"),
                            re,
                            im,
                        });
                    }
                }
            }
            r.solutions.push(s);
        }
        results.push((name, r));
    }
    Ok(Stored { mesh, project, results, fields_saved })
}

/// The size of a results file in bytes and how it divides: for the Share screen and the tests that compare with JSON.
pub fn file_bytes(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::Engine;
    use serde_json::json;

    fn solved(name: &str) -> Engine {
        let mut e = Engine::new();
        let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(format!("{name}.nyx.json"));
        e.call("project_open", &json!({ "path": p.to_str().unwrap() })).unwrap();
        e.call("solve_run", &json!({})).unwrap();
        e
    }

    #[test]
    fn the_blocks_of_the_processes_tile_the_rows() {
        for size in [1usize, 2, 3, 4, 7, 16] {
            for n in [0usize, 1, 5, 16, 17, 1000] {
                let (mut next, mut lo, mut hi) = (0, usize::MAX, 0);
                for rank in 0..size {
                    let r = Partition { rank, size }.rows(n);
                    assert_eq!(r.start, next, "contiguous (size {size}, n {n})");
                    next = r.end;
                    lo = lo.min(r.len());
                    hi = hi.max(r.len());
                }
                assert_eq!(next, n, "all rows covered (size {size}, n {n})");
                assert!(hi - lo <= 1, "balanced (size {size}, n {n})");
            }
        }
    }

    fn write_slice_f64(session: &File, path: &str, data: &[f64], cols: usize, part: Partition) -> Result<(), String> {
        let rows = data.len().checked_div(cols).unwrap_or(data.len());
        let r = part.rows(rows);
        if r.is_empty() {
            return Ok(());
        }
        let stride = cols.max(1);
        let mut ds = session.dataset(path).map_err(e(path))?;
        let shape = ds.shape().map_err(e(path))?;
        let mut buf = ds.read_f64().map_err(e(path))?;
        buf[r.start * stride..r.end * stride].copy_from_slice(&data[r.start * stride..r.end * stride]);
        ds.write_staged(|b| {
            b.with_f64_data(&buf).with_shape(&shape);
        })
        .map_err(e(path))
    }

    fn write_slice_u32(session: &File, path: &str, data: &[u32], cols: usize, part: Partition) -> Result<(), String> {
        let rows = data.len().checked_div(cols).unwrap_or(data.len());
        let r = part.rows(rows);
        if r.is_empty() {
            return Ok(());
        }
        let stride = cols.max(1);
        let mut ds = session.dataset(path).map_err(e(path))?;
        let shape = ds.shape().map_err(e(path))?;
        let mut buf = ds.read_u32().map_err(e(path))?;
        buf[r.start * stride..r.end * stride].copy_from_slice(&data[r.start * stride..r.end * stride]);
        ds.write_staged(|b| {
            b.with_u32_data(&buf).with_shape(&shape);
        })
        .map_err(e(path))
    }

    fn write_slice_i32(session: &File, path: &str, data: &[i32], part: Partition) -> Result<(), String> {
        let r = part.rows(data.len());
        if r.is_empty() {
            return Ok(());
        }
        let mut ds = session.dataset(path).map_err(e(path))?;
        let shape = ds.shape().map_err(e(path))?;
        let mut buf = ds.read_i32().map_err(e(path))?;
        buf[r.clone()].copy_from_slice(&data[r]);
        ds.write_staged(|b| {
            b.with_i32_data(&buf).with_shape(&shape);
        })
        .map_err(e(path))
    }

    fn write_slice_u8(session: &File, path: &str, data: &[u8], part: Partition) -> Result<(), String> {
        let r = part.rows(data.len());
        if r.is_empty() {
            return Ok(());
        }
        let mut ds = session.dataset(path).map_err(e(path))?;
        let shape = ds.shape().map_err(e(path))?;
        let mut buf = ds.read_u8().map_err(e(path))?;
        buf[r.clone()].copy_from_slice(&data[r]);
        ds.write_staged(|b| {
            b.with_u8_data(&buf).with_shape(&shape);
        })
        .map_err(e(path))
    }

    /// Emulates a multi-process row-partitioned write using `hdf5-pure`: rank 0 writes the zero-initialised scaffold of datasets,
    /// and each rank `0..part.size` opens the file via `File::open_rw`, writes only its `Partition::rows` slice, and commits.
    fn write_partitioned_emulated(
        path: &Path,
        part: Partition,
        project: &Project,
        mesh: &SimplexMesh,
        results: &[(String, &SolveResult)],
    ) -> Result<(), String> {
        let layout = FieldStorageLayout::L2ShuffleLzfPacked;
        let proj_json = project.to_json();
        let d = mesh.dim;
        let curved = mesh.edge_node_list();
        let n_cells = mesh.n_cells();

        if part.rank == 0 {
            let mut fb = FileBuilder::new();
            fb.set_attr("writer", text_attr("parallel")?);
            fb.set_attr("processes", AttrValue::I64(part.size as i64));
            fb.set_attr("format", text_attr("nyx-results")?);
            fb.set_attr("version", AttrValue::I32(4));
            fill_u8(fb.create_dataset("project"), &vec![0u8; proj_json.len()], layout);

            let mut m = fb.create_group("mesh");
            fill_f64(m.create_dataset("points"), &vec![0.0; mesh.points.len()], d, layout);
            fill_u32(m.create_dataset("cells"), &vec![0u32; mesh.cells.len()], d + 1, layout);
            fill_i32(m.create_dataset("cell_labels"), &vec![0i32; mesh.cell_labels.len()], 0, layout);
            fill_u32(m.create_dataset("facets"), &vec![0u32; mesh.facets.len()], d, layout);
            fill_i32(m.create_dataset("facet_labels"), &vec![0i32; mesh.facet_labels.len()], 0, layout);
            if !curved.is_empty() {
                fill_u32(m.create_dataset("edge_node_pairs"), &vec![0u32; curved.len() * 2], 2, layout);
                fill_f64(m.create_dataset("edge_node_points"), &vec![0.0; curved.len() * 3], 3, layout);
            }
            fb.add_group(m.finish());

            let mut root = fb.create_group("results");
            for (code, r) in results {
                let mut g = root.create_group(&link_name(code));
                g.set_attr("physics", text_attr(code)?);
                g.set_attr("study", text_attr(&r.study)?);
                g.set_attr("order", AttrValue::I64(r.order as i64));
                g.set_attr("dim", AttrValue::I64(r.dim as i64));
                g.set_attr("dofs", AttrValue::I64(r.dofs as i64));
                g.set_attr("seconds", AttrValue::F64(r.seconds));
                g.set_attr("n_solutions", AttrValue::I64(r.solutions.len() as i64));

                let ord = if r.order == 2 { 2 } else { 1 };
                let nloc_mesh = mesh.nloc(ord);
                let dof_pts = mesh.dof_coordinates(ord);
                fill_f64(g.create_dataset("points"), &vec![0.0; dof_pts.len()], mesh.dim, layout);
                fill_u32(g.create_dataset("cells"), &vec![0u32; n_cells * nloc_mesh], nloc_mesh, layout);
                g.set_attr("cell_order", text_attr(cell_order_text(mesh.dim, ord))?);
                g.set_attr("field_layout", text_attr(layout.as_str())?);

                let scalars_text = serialize_scalars(&r.scalars);
                fill_u8(g.create_dataset("scalars"), &vec![0u8; scalars_text.len()], layout);
                let log_text = serde_json::to_string(&r.log).map_err(|x| x.to_string())?;
                fill_u8(g.create_dataset("log"), &vec![0u8; log_text.len()], layout);
                fill_f64(g.create_dataset("eigenvalues"), &vec![0.0; r.solutions.len() * 2], 2, layout);

                if let Some((nloc_cf, ref meta)) = packable_cell_fields(r, n_cells) {
                    let c_comp = meta.len();
                    let ri = if meta.iter().any(|m| m.has_im) { 2usize } else { 1usize };
                    let zeros = vec![0.0; r.solutions.len() * c_comp * n_cells * nloc_cf * ri];
                    g.create_dataset("cell_fields_packed")
                        .with_f64_data(&zeros)
                        .with_shape(&[
                            r.solutions.len() as u64,
                            c_comp as u64,
                            n_cells as u64,
                            nloc_cf as u64,
                            ri as u64,
                        ])
                        .with_chunks(&[1, c_comp as u64, n_cells as u64, nloc_cf as u64, ri as u64])
                        .with_shuffle()
                        .with_lzf();
                    let meta_bytes = serde_json::to_vec(meta).map_err(|x| x.to_string())?;
                    fill_u8(g.create_dataset("cell_fields_meta"), &vec![0u8; meta_bytes.len()], layout);
                }

                for (k, s) in r.solutions.iter().enumerate() {
                    let mut sg = g.create_group(&format!("solution_{k}"));
                    sg.set_attr("label", text_attr(&s.label)?);
                    if let Some(v) = s.value {
                        sg.set_attr("value", AttrValue::F64(v));
                    }
                    if let Some([re, im]) = s.eigenvalue {
                        sg.set_attr("eig_re", AttrValue::F64(re));
                        sg.set_attr("eig_im", AttrValue::F64(im));
                    }
                    let mut fg = sg.create_group("fields");
                    let order: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
                    fg.set_attr("order", text_attr(&serde_json::to_string(&order).map_err(|x| x.to_string())?)?);
                    for (i, f) in s.fields.iter().enumerate() {
                        let comp_json = serde_json::to_string(&f.components).map_err(|x| x.to_string())?;
                        let ds = fill_f64(
                            fg.create_dataset(&format!("{i:02}_{}", link_name(&f.name))),
                            &vec![0.0; f.values.len()],
                            0,
                            layout,
                        );
                        ds.set_attr("name", text_attr(&f.name)?)
                            .set_attr("unit", text_attr(&f.unit)?)
                            .set_attr("components", text_attr(&comp_json)?);
                    }
                    sg.add_group(fg.finish());
                    g.add_group(sg.finish());
                }
                root.add_group(g.finish());
            }
            fb.add_group(root.finish());
            fb.write(path).map_err(e("write scaffold"))?;
        }

        let session = File::open_rw(path).map_err(e("open_rw"))?;
        write_slice_u8(&session, "project", proj_json.as_bytes(), part)?;
        write_slice_f64(&session, "mesh/points", &mesh.points, d, part)?;
        write_slice_u32(&session, "mesh/cells", &mesh.cells, d + 1, part)?;
        write_slice_i32(&session, "mesh/cell_labels", mesh.cell_labels.as_slice(), part)?;
        write_slice_u32(&session, "mesh/facets", &mesh.facets, d, part)?;
        write_slice_i32(&session, "mesh/facet_labels", mesh.facet_labels.as_slice(), part)?;
        if !curved.is_empty() {
            let pairs: Vec<u32> = curved.iter().flat_map(|e| e.0).collect();
            let pos: Vec<f64> = curved.iter().flat_map(|e| e.1).collect();
            write_slice_u32(&session, "mesh/edge_node_pairs", &pairs, 2, part)?;
            write_slice_f64(&session, "mesh/edge_node_points", &pos, 3, part)?;
        }

        for (code, r) in results {
            let base = format!("results/{}", link_name(code));
            let ord = if r.order == 2 { 2 } else { 1 };
            let nloc_mesh = mesh.nloc(ord);
            let mut conn: Vec<u32> = Vec::with_capacity(n_cells * nloc_mesh);
            let mut dofs = Vec::new();
            for c in 0..n_cells {
                mesh.cell_dofs(c, ord, &mut dofs);
                conn.extend(dofs.iter().take(nloc_mesh).map(|&idx| idx as u32));
            }
            let dof_pts = mesh.dof_coordinates(ord);
            write_slice_f64(&session, &format!("{base}/points"), &dof_pts, mesh.dim, part)?;
            write_slice_u32(&session, &format!("{base}/cells"), &conn, nloc_mesh, part)?;

            let scalars_text = serialize_scalars(&r.scalars);
            write_slice_u8(&session, &format!("{base}/scalars"), scalars_text.as_bytes(), part)?;
            let log_text = serde_json::to_string(&r.log).map_err(|x| x.to_string())?;
            write_slice_u8(&session, &format!("{base}/log"), log_text.as_bytes(), part)?;

            let mut eig_flat = Vec::with_capacity(r.solutions.len() * 2);
            for s in &r.solutions {
                let [re, im] = s.eigenvalue.unwrap_or([s.value.unwrap_or(0.0), 0.0]);
                eig_flat.push(re);
                eig_flat.push(im);
            }
            write_slice_f64(&session, &format!("{base}/eigenvalues"), &eig_flat, 2, part)?;

            if let Some((nloc_cf, ref meta)) = packable_cell_fields(r, n_cells) {
                let field_len = n_cells * nloc_cf;
                let (ri, packed) = build_packed_cell_fields(r, meta, field_len);
                let mode_stride = meta.len() * field_len * ri;
                write_slice_f64(&session, &format!("{base}/cell_fields_packed"), &packed, mode_stride, part)?;
                let meta_bytes = serde_json::to_vec(meta).map_err(|x| x.to_string())?;
                write_slice_u8(&session, &format!("{base}/cell_fields_meta"), &meta_bytes, part)?;
            }

            for (k, s) in r.solutions.iter().enumerate() {
                for (i, f) in s.fields.iter().enumerate() {
                    let fp = format!("{base}/solution_{k}/fields/{i:02}_{}", link_name(&f.name));
                    write_slice_f64(&session, &fp, &f.values, 0, part)?;
                }
            }
        }

        session.commit().map_err(e("commit"))
    }

    /// Three emulated processes each write only their `Partition::rows` slice into one file via `hdf5-pure` staged appends,
    /// producing bit-for-bit what one serial process writes (same mesh, same fields, same cell fields, same project).
    #[test]
    fn three_emulated_processes_write_what_one_serial_process_writes() {
        let e = solved("15_soi_wire");
        let mesh = e.model().unwrap().solver_mesh();
        let results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = e
            .project()
            .physics
            .iter()
            .filter_map(|s| e.result(&s.physics).map(|(r, _)| (s.physics.clone(), r.clone())))
            .collect();
        let refs: Vec<(String, &nyx_physics::SolveResult)> =
            results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();
        let tmp = |n: &str| std::env::temp_dir().join(format!("nyx-part-{n}-{}.h5", std::process::id()));
        let (serial, three) = (tmp("serial"), tmp("three"));
        write_results(&serial, e.project(), &mesh, &refs).unwrap();
        for rank in 0..3 {
            write_partitioned_emulated(&three, Partition { rank, size: 3 }, e.project(), &mesh, &refs).unwrap();
        }
        let (a, b) = (read_results(&serial).unwrap(), read_results(&three).unwrap());
        assert_eq!(a.project, b.project);
        assert_eq!(a.mesh.points, b.mesh.points);
        assert_eq!(a.mesh.cells, b.mesh.cells);
        assert_eq!(a.mesh.facets, b.mesh.facets);
        assert_eq!(a.results.len(), b.results.len());
        for ((ca, ra), (cb, rb)) in a.results.iter().zip(&b.results) {
            assert_eq!(ca, cb);
            assert_eq!(ra, rb, "{ca}: every solution, field and cell field is the same");
        }
        // the nodes and cells of a physics, which the XDMF file names, are partitioned too
        let (fa, fb) = (File::open_streaming(&serial).unwrap(), File::open_streaming(&three).unwrap());
        let code = &results[0].0;
        for leaf in ["points", "cells"] {
            let p = format!("results/{}/{leaf}", link_name(code));
            let (da, db) = (fa.dataset(&p).unwrap(), fb.dataset(&p).unwrap());
            assert_eq!(da.shape().unwrap(), db.shape().unwrap(), "{p}");
            if leaf == "cells" {
                assert_eq!(da.read_u32().unwrap(), db.read_u32().unwrap(), "{p}");
            } else {
                assert_eq!(da.read_f64().unwrap(), db.read_f64().unwrap(), "{p}");
            }
        }
        let fb_attrs = fb.root().attrs().unwrap();
        assert_eq!(attr_str(&fb_attrs, "writer").as_deref(), Some("parallel"));
        assert_eq!(attr_i64(&fb_attrs, "processes"), 3);
        let _ = std::fs::remove_file(&serial);
        let _ = std::fs::remove_file(&three);
    }

    #[test]
    fn parallel_mode_says_what_is_missing_when_the_build_has_no_mpi() {
        assert!(!mode_available(H5Mode::Parallel));
        let e = solved("04_thermal_bar");
        let mesh = e.model().unwrap().solver_mesh();
        let err = write_results_mode(&std::env::temp_dir().join("never.h5"), H5Mode::Parallel, e.project(), &mesh, &[])
            .unwrap_err();
        assert!(err.contains("phdf5") && err.contains("MPI"), "{err}");
        assert!(start_runtime(H5Mode::Parallel).is_err());
        assert_eq!(start_runtime(H5Mode::Serial), Ok((0, 1)));
    }
}
