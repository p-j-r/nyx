// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Results in HDF5 files (cargo feature `hdf5`): the project, the mesh and every solution of every solved physics with all its fields and
//! cell fields, written compressed and read back exactly (ADR 62). This is the efficient store of the data a solve produces; the JSON files
//! of a project and a run record stay small.
//!
//! Layout, version 4 (numbers as stored by HDF5, little-endian; version 3 files, which lack `points`, `cells` and `cell_order` of a physics,
//! are still read):
//!
//! ```text
//! /                          attrs: format = "nyx-results", version = 4
//! /project                   the project as JSON text (u8 dataset)
//! /mesh/points               [n_points, dim] f64         /mesh/cells        [n_cells, dim+1] u32
//! /mesh/cell_labels          [n_cells] i32               /mesh/facets       [n_facets, dim] u32
//! /mesh/facet_labels         [n_facets] i32
//! /mesh/edge_node_pairs      [n_curved, 2] u32           /mesh/edge_node_points  [n_curved, 3] f64   (optional: curved edges, ADR 82; absent = straight)
//! /results/<physics>         attrs: study, order, dim, dofs, seconds, n_solutions, cell_order
//!   points                   [n_dofs, dim] f64         coordinates of the nodes (the dofs): vertices first, then edge midpoints
//!   cells                    [n_cells, nloc] u32       the node (dof) of every cell node, in the cell node order of `cell_order`
//!   scalars, log             JSON text (u8 datasets)
//!   solution_<k>             attrs: label, value (when present), eig_re, eig_im (when present)
//!     fields/<name>          f64, `n_dofs * components` values interleaved by dof; attrs: unit, components (JSON), name
//!     cell_fields/<name>/re, im   f64 per cell node; attrs on the group: unit, nloc
//! ```
//!
//! **One node order.** `fields/<name>` is indexed by node (dof), as `points`. `cell_fields/<name>/re[c * nloc + a]` is the value at the
//! node `cells[c, a]` as seen from cell `c`, so any reader gets the nodal complex field by `nodal[cells[c, a]] = cell_field[c * nloc + a]`
//! (for a continuous field every cell that holds the node gives the same value). The cell node order is the one of nyx (`cell_order`
//! says it): a P1 cell lists its vertices; a P2 triangle `v0 v1 v2 e01 e02 e12`; a P2 tetrahedron `v0 v1 v2 v3 e01 e02 e03 e12 e13 e23`,
//! where `eij` is the midpoint of the edge between vertices `i` and `j`. `/mesh/cells` lists the vertices of a cell (the first `dim + 1`
//! entries of its row of `cells`). There is no second copy of the nodes or cells: the XDMF file next to the `.h5` (`xdmf.rs`) is a small XML
//! index that names `points`, `cells` and the stored fields as they are, for nyx and its tools only (not for ParaView or other viewers).
//!
//! Large datasets are chunked and deflate-compressed (level 3); a result written and read back is bit-identical.

use crate::project::Project;
use femcore::mesh::SimplexMesh;
use h5::types::VarLenUnicode;
use hdf5_metno as h5;
use nyx_physics::{CellField, Field, Scalar, Solution, SolveResult};
use serde_json::{json, Value};
use std::path::Path;

pub(crate) fn e<T: std::fmt::Display>(what: &str) -> impl Fn(T) -> String + '_ {
    move |x| format!("HDF5 {what}: {x}")
}

fn text(s: &str) -> Result<VarLenUnicode, String> {
    s.parse::<VarLenUnicode>().map_err(|_| "text could not be stored (it contains a NUL character)".to_string())
}

/// A dataset of shape `[rows, cols]` (`[rows]` when `cols` is 0) holding `data` row by row, chunked and compressed when it is large.
fn put_serial<T: h5::H5Type>(g: &h5::Group, name: &str, data: &[T], cols: usize) -> Result<h5::Dataset, String> {
    let rows = data.len().checked_div(cols).unwrap_or(data.len());
    let big = data.len() >= 4096 && rows > 0;
    let ds = if cols == 0 {
        let b = g.new_dataset::<T>().shape([rows]);
        if big {
            b.chunk([rows.min(32768)]).deflate(3).create(name)
        } else {
            b.create(name)
        }
    } else {
        let b = g.new_dataset::<T>().shape([rows, cols]);
        if big {
            b.chunk([rows.min((32768 / cols).max(1)), cols]).deflate(3).create(name)
        } else {
            b.create(name)
        }
    }
    .map_err(e(name))?;
    if !data.is_empty() {
        ds.write_raw(data).map_err(e(name))?;
    }
    Ok(ds)
}

/// Name usable as an HDF5 link: `/` is the path separator, everything else is allowed.
pub(crate) fn link_name(name: &str) -> String {
    let n: String = name.chars().map(|c| if c == '/' { '_' } else { c }).collect();
    if n.is_empty() || n == "." {
        format!("_{n}")
    } else {
        n
    }
}

/// How a results file is written: with the ordinary sequential HDF5 driver, or with the MPI-IO driver by several processes at once.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H5Mode {
    Serial,
    Parallel,
}

/// True if this build can write with the mode: serial always, parallel only when it was built with the `phdf5` feature (an MPI-enabled
/// HDF5 library and an MPI SDK).
pub fn mode_available(mode: H5Mode) -> bool {
    match mode {
        H5Mode::Serial => true,
        H5Mode::Parallel => cfg!(feature = "phdf5"),
    }
}

#[cfg(not(feature = "phdf5"))]
const NO_PARALLEL: &str = "parallel HDF5 (MPI-IO) is not available in this build: it was made without the `phdf5` feature (an MPI-enabled HDF5 library and an MPI SDK are needed to build it); use the serial output";

/// Starts what the mode needs before a solve: the MPI runtime for parallel output (a run under `mpiexec` has one process per rank).
/// Returns this process's rank and the number of processes.
pub fn start_runtime(mode: H5Mode) -> Result<(usize, usize), String> {
    match mode {
        H5Mode::Serial => Ok((0, 1)),
        #[cfg(feature = "phdf5")]
        H5Mode::Parallel => {
            let w = mpi_rt::world(true)?.ok_or("MPI is not available")?;
            Ok((w.rank, w.size))
        }
        #[cfg(not(feature = "phdf5"))]
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

/// The MPI process group of a parallel build (raw `mpi-sys` calls: initialise, rank, size, barrier).
#[cfg(feature = "phdf5")]
mod mpi_rt {
    use mpi_sys as ffi;
    use std::os::raw::c_int;

    pub struct World {
        pub comm: ffi::MPI_Comm,
        pub rank: usize,
        pub size: usize,
    }

    /// MPI_COMM_WORLD. With `init` the library is initialised when this process has not done so (a run under `mpiexec` or on a cluster does it
    /// here); without it a process that never started MPI gets `None`.
    pub fn world(init: bool) -> Result<Option<World>, String> {
        unsafe {
            let mut flag: c_int = 0;
            ffi::MPI_Initialized(&mut flag);
            if flag == 0 {
                if !init {
                    return Ok(None);
                }
                if ffi::MPI_Init(std::ptr::null_mut(), std::ptr::null_mut()) != ffi::MPI_SUCCESS as c_int {
                    return Err("MPI_Init failed".into());
                }
            }
            let comm = ffi::RSMPI_COMM_WORLD;
            let (mut r, mut n): (c_int, c_int) = (0, 1);
            ffi::MPI_Comm_rank(comm, &mut r);
            ffi::MPI_Comm_size(comm, &mut n);
            Ok(Some(World { comm, rank: r as usize, size: n as usize }))
        }
    }

    pub fn barrier(w: &World) {
        unsafe {
            ffi::MPI_Barrier(w.comm);
        }
    }
}

/// An open results file and the share of it this process writes. Every process runs the same sequence of group, dataset and attribute
/// creations (HDF5 wants metadata operations to be collective in parallel mode) and writes only its own rows of each dataset. The serial
/// path is the same code with one process (and then chunks and compresses large datasets; the parallel path writes contiguous datasets).
pub(crate) struct Io {
    pub file: h5::File,
    pub part: Partition,
    /// tests only: run the writer again for another rank against a file that exists, opening what is there instead of creating it, to
    /// check that the row blocks tile the data (the MPI-IO driver itself cannot run without an MPI-enabled HDF5 library)
    pub replay: bool,
}

impl Io {
    pub fn group(&self, parent: &h5::Group, name: &str) -> Result<h5::Group, String> {
        if self.replay {
            if let Ok(g) = parent.group(name) {
                return Ok(g);
            }
        }
        parent.create_group(name).map_err(e("group"))
    }

    pub fn attr<T: h5::H5Type>(&self, loc: &h5::Location, name: &str, v: &T) -> Result<(), String> {
        if self.replay && loc.attr(name).is_ok() {
            return Ok(());
        }
        loc.new_attr::<T>().create(name).map_err(e("attribute"))?.write_scalar(v).map_err(e("attribute"))
    }

    pub fn attr_text(&self, loc: &h5::Location, name: &str, s: &str) -> Result<(), String> {
        self.attr(loc, name, &text(s)?)
    }

    /// A dataset of shape `[rows, cols]` (`[rows]` when `cols` is 0) holding `data` row by row.
    pub fn put<T: h5::H5Type>(
        &self,
        g: &h5::Group,
        name: &str,
        data: &[T],
        cols: usize,
    ) -> Result<h5::Dataset, String> {
        if self.part == Partition::SOLO && !self.replay {
            return put_serial(g, name, data, cols);
        }
        let rows = data.len().checked_div(cols).unwrap_or(data.len());
        let existing = if self.replay { g.dataset(name).ok() } else { None };
        let ds = match existing {
            Some(d) => d,
            None if cols == 0 => g.new_dataset::<T>().shape([rows]).create(name).map_err(e(name))?,
            None => g.new_dataset::<T>().shape([rows, cols]).create(name).map_err(e(name))?,
        };
        let r = self.part.rows(rows);
        if r.is_empty() {
            return Ok(ds);
        }
        if cols == 0 {
            let v = ndarray::ArrayView1::from(&data[r.clone()]);
            ds.write_slice(v, ndarray::s![r.start..r.end]).map_err(e(name))?;
        } else {
            let v = ndarray::ArrayView2::from_shape((r.len(), cols), &data[r.start * cols..r.end * cols])
                .map_err(|x| format!("HDF5 {name}: {x}"))?;
            ds.write_slice(v, ndarray::s![r.start..r.end, ..]).map_err(e(name))?;
        }
        Ok(ds)
    }

    /// Text as UTF-8 bytes in a 1-D `u8` dataset (variable-length strings cannot be written in parallel).
    pub fn put_text(&self, g: &h5::Group, name: &str, s: &str) -> Result<(), String> {
        self.put(g, name, s.as_bytes(), 0).map(|_| ())
    }
}

/// Opens the results file for writing in the mode. `Ok(None)`: this process has nothing to write (serial output chosen in a multi-process
/// MPI run: only rank 0 writes).
fn open_io(path: &Path, mode: H5Mode) -> Result<Option<Io>, String> {
    match mode {
        H5Mode::Serial => {
            #[cfg(feature = "phdf5")]
            if let Some(w) = mpi_rt::world(false)? {
                if w.rank != 0 {
                    return Ok(None);
                }
            }
            Ok(Some(Io {
                file: h5::File::create(path).map_err(e("create file"))?,
                part: Partition::SOLO,
                replay: false,
            }))
        }
        #[cfg(feature = "phdf5")]
        H5Mode::Parallel => {
            let w = mpi_rt::world(true)?.ok_or("MPI is not available")?;
            let file = h5::File::with_options()
                .with_fapl(|p| p.mpio(w.comm, None))
                .create(path)
                .map_err(e("create file (MPI-IO)"))?;
            Ok(Some(Io { file, part: Partition { rank: w.rank, size: w.size }, replay: false }))
        }
        #[cfg(not(feature = "phdf5"))]
        H5Mode::Parallel => {
            let _ = path;
            Err(NO_PARALLEL.into())
        }
    }
}

/// Writes `results` (physics code and result) solved on `mesh` for `project` into a new file at `path`, serially.
pub fn write_results(
    path: &Path,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<(), String> {
    write_results_mode(path, H5Mode::Serial, project, mesh, results).map(|_| ())
}

/// The same in the chosen mode. In a parallel run every process calls this with the same arguments; the file is the same as a serial one.
/// Returns whether this process is the one that writes the XDMF file (process 0, after every process is done).
pub fn write_results_mode(
    path: &Path,
    mode: H5Mode,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<bool, String> {
    let Some(io) = open_io(path, mode)? else { return Ok(false) };
    let root = io.part.rank == 0;
    write_into(&io, mode, project, mesh, results)?;
    drop(io);
    #[cfg(feature = "phdf5")]
    if mode == H5Mode::Parallel {
        if let Some(w) = mpi_rt::world(false)? {
            mpi_rt::barrier(&w);
        }
    }
    Ok(root)
}

fn write_into(
    io: &Io,
    mode: H5Mode,
    project: &Project,
    mesh: &SimplexMesh,
    results: &[(String, &SolveResult)],
) -> Result<(), String> {
    let file = &io.file;
    io.attr_text(file, "writer", if mode == H5Mode::Serial { "serial" } else { "parallel" })?;
    io.attr(file, "processes", &(io.part.size as i64))?;
    io.attr_text(file, "format", "nyx-results")?;
    io.attr(file, "version", &4i32)?;
    io.put_text(file, "project", &project.to_json())?;
    let d = mesh.dim;
    let m = io.group(file, "mesh")?;
    io.put(&m, "points", &mesh.points, d)?;
    io.put(&m, "cells", &mesh.cells, d + 1)?;
    io.put(&m, "cell_labels", mesh.cell_labels.as_slice(), 0)?;
    io.put(&m, "facets", &mesh.facets, d)?;
    io.put(&m, "facet_labels", mesh.facet_labels.as_slice(), 0)?;
    // curved edges (ADR 82): written only when there are some, so a straight mesh gives the file it always gave
    let curved = mesh.edge_node_list();
    if !curved.is_empty() {
        let pairs: Vec<u32> = curved.iter().flat_map(|e| e.0).collect();
        let pos: Vec<f64> = curved.iter().flat_map(|e| e.1).collect();
        io.put(&m, "edge_node_pairs", &pairs, 2)?;
        io.put(&m, "edge_node_points", &pos, 3)?;
    }
    let root = io.group(file, "results")?;
    for (code, r) in results {
        let g = io.group(&root, &link_name(code))?;
        io.attr_text(&g, "physics", code)?;
        io.attr_text(&g, "study", &r.study)?;
        io.attr(&g, "order", &(r.order as i64))?;
        io.attr(&g, "dim", &(r.dim as i64))?;
        io.attr(&g, "dofs", &(r.dofs as i64))?;
        io.attr(&g, "seconds", &r.seconds)?;
        io.attr(&g, "n_solutions", &(r.solutions.len() as i64))?;
        // the nodes of this physics and the node of every cell node, in nyx's own cell node order (the order of every cell field)
        {
            let ord = if r.order == 2 { 2 } else { 1 };
            let nloc = mesh.nloc(ord);
            let mut conn: Vec<u32> = Vec::with_capacity(mesh.n_cells() * nloc);
            let mut dofs = Vec::new();
            for c in 0..mesh.n_cells() {
                mesh.cell_dofs(c, ord, &mut dofs);
                conn.extend(dofs.iter().take(nloc).map(|&d| d as u32));
            }
            io.put(&g, "points", &mesh.dof_coordinates(ord), mesh.dim)?;
            io.put(&g, "cells", &conn, nloc)?;
            io.attr_text(&g, "cell_order", cell_order_text(mesh.dim, ord))?;
        }
        let scalars: Vec<Value> = r
            .scalars
            .iter()
            .map(|(n, v)| match v {
                Scalar::Num(x) => json!([n, x]),
                Scalar::Text(t) => json!([n, t]),
            })
            .collect();
        io.put_text(&g, "scalars", &Value::Array(scalars).to_string())?;
        io.put_text(&g, "log", &serde_json::to_string(&r.log).map_err(|x| x.to_string())?)?;
        for (k, s) in r.solutions.iter().enumerate() {
            let sg = io.group(&g, &format!("solution_{k}"))?;
            io.attr_text(&sg, "label", &s.label)?;
            if let Some(v) = s.value {
                io.attr(&sg, "value", &v)?;
            }
            if let Some([re, im]) = s.eigenvalue {
                io.attr(&sg, "eig_re", &re)?;
                io.attr(&sg, "eig_im", &im)?;
            }
            let fg = io.group(&sg, "fields")?;
            let order: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
            io.attr_text(&fg, "order", &serde_json::to_string(&order).map_err(|x| x.to_string())?)?;
            for (i, f) in s.fields.iter().enumerate() {
                let ds = io.put(&fg, &format!("{i:02}_{}", link_name(&f.name)), &f.values, 0)?;
                io.attr_text(&ds, "name", &f.name)?;
                io.attr_text(&ds, "unit", &f.unit)?;
                io.attr_text(&ds, "components", &serde_json::to_string(&f.components).map_err(|x| x.to_string())?)?;
            }
            if !s.cell_fields.is_empty() {
                let cg = io.group(&sg, "cell_fields")?;
                for (i, f) in s.cell_fields.iter().enumerate() {
                    let fgp = io.group(&cg, &format!("{i:02}_{}", link_name(&f.name)))?;
                    io.attr_text(&fgp, "name", &f.name)?;
                    io.attr_text(&fgp, "unit", &f.unit)?;
                    io.attr(&fgp, "nloc", &(f.nloc as i64))?;
                    io.put(&fgp, "re", &f.re, 0)?;
                    if !f.im.is_empty() {
                        io.put(&fgp, "im", &f.im, 0)?;
                    }
                }
            }
        }
    }
    Ok(())
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
    if !matches!(attr_i64(&file, "version"), 3 | 4) {
        return Err("this results file has another version than this nyx reads (3 and 4)".into());
    }
    let project = read_text(&file, "project")?;
    let m = file.group("mesh").map_err(e("mesh group"))?;
    let pts = m.dataset("points").map_err(e("points"))?;
    let dim = pts.shape().get(1).copied().ok_or("bad points shape")?;
    let mut mesh = SimplexMesh::new(
        dim,
        pts.read_raw::<f64>().map_err(e("points"))?,
        m.dataset("cells").map_err(e("cells"))?.read_raw::<u32>().map_err(e("cells"))?,
        m.dataset("cell_labels").map_err(e("cell labels"))?.read_raw::<i32>().map_err(e("cell labels"))?,
        m.dataset("facets").map_err(e("facets"))?.read_raw::<u32>().map_err(e("facets"))?,
        m.dataset("facet_labels").map_err(e("facet labels"))?.read_raw::<i32>().map_err(e("facet labels"))?,
    );
    // optional curved edges: a file without them is a straight mesh
    if let (Ok(p), Ok(q)) = (m.dataset("edge_node_pairs"), m.dataset("edge_node_points")) {
        let pairs = p.read_raw::<u32>().map_err(e("edge node pairs"))?;
        let pos = q.read_raw::<f64>().map_err(e("edge node points"))?;
        let nodes: Vec<([u32; 2], [f64; 3])> =
            pairs.chunks(2).zip(pos.chunks(3)).map(|(a, b)| ([a[0], a[1]], [b[0], b[1], b[2]])).collect();
        mesh.set_edge_nodes(&nodes);
    }
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
                    r.scalars.push((
                        n.to_string(),
                        match v.as_f64() {
                            Some(x) => Scalar::Num(x),
                            None => Scalar::Text(v.as_str().unwrap_or("").to_string()),
                        },
                    ));
                }
            }
        }
        r.log = serde_json::from_str(&read_text(&g, "log")?).unwrap_or_default();
        for k in 0..attr_i64(&g, "n_solutions") {
            let sg = g.group(&format!("solution_{k}")).map_err(e("solution"))?;
            let mut s = Solution {
                label: attr_str(&sg, "label").unwrap_or_default(),
                value: attr_f64(&sg, "value"),
                ..Default::default()
            };
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
                    components: attr_str(&ds, "components")
                        .and_then(|t| serde_json::from_str(&t).ok())
                        .unwrap_or_default(),
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

    /// The MPI-IO driver cannot run without an MPI-enabled HDF5 library, but the logic can: three emulated processes each write only
    /// their own rows into one file, which must read back as what one serial process writes (same mesh, same fields, same project).
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
        {
            let file = h5::File::create(&three).unwrap();
            for rank in 0..3 {
                let io = Io { file: file.clone(), part: Partition { rank, size: 3 }, replay: rank > 0 };
                write_into(&io, H5Mode::Parallel, e.project(), &mesh, &refs).unwrap();
            }
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
        let (fa, fb) = (h5::File::open(&serial).unwrap(), h5::File::open(&three).unwrap());
        let code = &results[0].0;
        for leaf in ["points", "cells"] {
            let p = format!("results/{}/{leaf}", link_name(code));
            let (da, db) = (fa.dataset(&p).unwrap(), fb.dataset(&p).unwrap());
            assert_eq!(da.shape(), db.shape(), "{p}");
            if leaf == "cells" {
                assert_eq!(da.read_raw::<u32>().unwrap(), db.read_raw::<u32>().unwrap(), "{p}");
            } else {
                assert_eq!(da.read_raw::<f64>().unwrap(), db.read_raw::<f64>().unwrap(), "{p}");
            }
        }
        assert_eq!(fb.attr("writer").unwrap().read_scalar::<VarLenUnicode>().unwrap().as_str(), "parallel");
        assert_eq!(fb.attr("processes").unwrap().read_scalar::<i64>().unwrap(), 3);
        let _ = std::fs::remove_file(&serial);
        let _ = std::fs::remove_file(&three);
    }

    #[test]
    fn parallel_mode_says_what_is_missing_when_the_build_has_no_mpi() {
        if !mode_available(H5Mode::Parallel) {
            let e = solved("04_thermal_bar");
            let mesh = e.model().unwrap().solver_mesh();
            let err =
                write_results_mode(&std::env::temp_dir().join("never.h5"), H5Mode::Parallel, e.project(), &mesh, &[])
                    .unwrap_err();
            assert!(err.contains("phdf5") && err.contains("MPI"), "{err}");
            assert!(start_runtime(H5Mode::Parallel).is_err());
        }
        assert_eq!(start_runtime(H5Mode::Serial), Ok((0, 1)));
    }
}
