//! Results in HDF5 files (feature `hdf5`): the mesh, the problem description and every solution
//! with all fields, written and read back.
//!
//! Layout (all numbers little-endian doubles/ints as stored by HDF5):
//!
//! ```text
//! /                       attrs: format = "nyx-results", version = 2, writer = "serial" | "parallel", processes
//! /spec                   JSON of the ProblemSpec as UTF-8 bytes (u8 dataset), optional
//! /mesh/points            [n_points, dim]  f64
//! /mesh/cells             [n_cells, dim+1] u32          /mesh/cell_labels   [n_cells]  i32
//! /mesh/facets            [n_facets, dim]  u32          /mesh/facet_labels  [n_facets] i32
//! /results                attrs: order, dim, dofs, seconds, n_solutions
//! /results/scalars        JSON array of [name, value-or-text] (UTF-8 bytes)
//! /results/log            JSON array of strings (UTF-8 bytes)
//! /results/solution_<k>   attrs: label, value (time / frequency), neff_re, neff_im (when present)
//!   fields/<name>         [n_dofs] f64, attr unit
//!   disp                  [n_dofs, dim] f64 (solid mechanics)
//! ```
//!
//! /viz/...                 the same data in the shape XDMF needs (quadratic cells in VTK node order,
//!                         three-component vectors), referenced by the `.xdmf` file
//!
//! Dual-mode I/O: one writer, two drivers, chosen per simulation at run time.
//! * `H5Mode::Serial` uses the standard sequential HDF5 driver and needs nothing but the HDF5 library.
//! * `H5Mode::Parallel` uses the MPI-IO driver (cargo feature `phdf5`, which needs an MPI-enabled
//!   HDF5 library and an MPI SDK to build). Every process makes the same collective calls and writes
//!   its own block of rows of each dataset (`Partition::rows`), so the file is the same as a serial one.
//!
//! Builds without `phdf5` contain no MPI code at all; asking them for parallel output returns an error.

use crate::mesh::SimplexMesh;
use crate::problem::*;
use hdf5_metno as h5;
use h5::types::VarLenUnicode;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H5Mode {
    Serial,
    Parallel,
}

fn e<T: std::fmt::Display>(what: &str) -> impl Fn(T) -> String + '_ {
    move |x| format!("HDF5 {what}: {x}")
}

fn text(s: &str) -> Result<VarLenUnicode, String> {
    s.parse::<VarLenUnicode>().map_err(|_| "text could not be stored (contains a NUL character)".to_string())
}

/// True if this build can write with the given mode: serial always, parallel only when the
/// program was built with the `phdf5` feature (MPI-enabled HDF5 library and an MPI SDK).
pub fn mode_available(mode: H5Mode) -> bool {
    match mode {
        H5Mode::Serial => true,
        H5Mode::Parallel => cfg!(feature = "phdf5"),
    }
}

/// Starts what the mode needs before the solve: the MPI runtime for parallel output (a run under
/// `mpiexec` has one process per rank). Returns this process's rank and the number of processes.
pub fn start_runtime(mode: H5Mode) -> Result<(usize, usize), String> {
    match mode {
        H5Mode::Serial => Ok((0, 1)),
        H5Mode::Parallel => {
            #[cfg(feature = "phdf5")]
            {
                let w = mpi_rt::world(true)?.ok_or("MPI is not available")?;
                return Ok((w.rank, w.size));
            }
            #[cfg(not(feature = "phdf5"))]
            Err("parallel HDF5 (MPI-IO) is not available in this build: it was made without the phdf5 feature".into())
        }
    }
}

/// Process `rank` of `size` that write one results file together. Each process owns a contiguous
/// block of rows of every dataset; a serial run is the single process 0 of 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Partition {
    pub rank: usize,
    pub size: usize,
}

impl Partition {
    pub const SOLO: Partition = Partition { rank: 0, size: 1 };

    /// The rows `[start, end)` of `n` that this process writes; the blocks of all processes
    /// tile `0..n` exactly, with sizes differing by at most one.
    pub fn rows(&self, n: usize) -> std::ops::Range<usize> {
        let size = self.size.max(1);
        let (base, rem) = (n / size, n % size);
        let start = self.rank * base + self.rank.min(rem);
        start..start + base + usize::from(self.rank < rem)
    }
}

/// MPI process group for the parallel build (raw `mpi-sys` calls: initialise, rank, size, barrier).
#[cfg(feature = "phdf5")]
mod mpi_rt {
    use mpi_sys as ffi;
    use std::os::raw::c_int;

    pub struct World {
        pub comm: ffi::MPI_Comm,
        pub rank: usize,
        pub size: usize,
    }

    /// MPI_COMM_WORLD. With `init` the library is initialised when this process has not done so
    /// (a run under `mpiexec` or on a cluster does it here); without it a process that never
    /// started MPI gets `None`.
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

/// An open results file and the share of it this process writes.
///
/// Every process runs the same sequence of group, dataset and attribute creations (HDF5 requires
/// metadata operations to be collective in parallel mode) and then writes only its own rows of each
/// dataset (`Partition::rows`). The serial path is the same code with one process.
struct Io {
    file: h5::File,
    part: Partition,
    /// tests only: re-run the writer for another rank against a file that already exists by
    /// opening what is there instead of creating it, to check that the row blocks tile the data
    replay: bool,
}

impl Io {
    fn group(&self, parent: &h5::Group, name: &str) -> Result<h5::Group, String> {
        if self.replay {
            if let Ok(g) = parent.group(name) {
                return Ok(g);
            }
        }
        parent.create_group(name).map_err(e("group"))
    }

    fn attr<T: h5::H5Type>(&self, loc: &h5::Location, name: &str, v: &T) -> Result<(), String> {
        if self.replay && loc.attr(name).is_ok() {
            return Ok(());
        }
        loc.new_attr::<T>().create(name).map_err(e("attribute"))?.write_scalar(v).map_err(e("attribute"))
    }

    fn attr_text(&self, loc: &h5::Location, name: &str, s: &str) -> Result<(), String> {
        self.attr(loc, name, &text(s)?)
    }

    fn link(&self, g: &h5::Group, target: &str, name: &str) -> Result<(), String> {
        if self.replay && g.link_exists(name) {
            return Ok(());
        }
        g.link_hard(target, name).map_err(e("link"))
    }

    /// A dataset of shape `[rows, cols]` (`[rows]` when `cols` is 0) holding `data` row by row.
    fn put<T: h5::H5Type>(&self, g: &h5::Group, name: &str, data: &[T], cols: usize) -> Result<h5::Dataset, String> {
        let rows = if cols == 0 { data.len() } else { data.len() / cols };
        let existing = if self.replay { g.dataset(name).ok() } else { None };
        let ds = match existing {
            Some(d) => d,
            None if cols == 0 => g.new_dataset::<T>().shape([rows]).create(name).map_err(e(name))?,
            None => g.new_dataset::<T>().shape([rows, cols]).create(name).map_err(e(name))?,
        };
        if self.part == Partition::SOLO && !self.replay {
            ds.write_raw(data).map_err(e(name))?;
            return Ok(ds);
        }
        let r = self.part.rows(rows);
        if r.is_empty() {
            return Ok(ds);
        }
        if cols == 0 {
            let v = ndarray::ArrayView1::from(&data[r.clone()]);
            ds.write_slice(v, ndarray::s![r.start..r.end]).map_err(e(name))?;
        } else {
            let v = ndarray::ArrayView2::from_shape((r.len(), cols), &data[r.start * cols..r.end * cols]).map_err(|x| format!("HDF5 {name}: {x}"))?;
            ds.write_slice(v, ndarray::s![r.start..r.end, ..]).map_err(e(name))?;
        }
        Ok(ds)
    }

    /// Text as UTF-8 bytes in a 1-D `u8` dataset (variable-length strings cannot be written in parallel).
    fn put_text(&self, g: &h5::Group, name: &str, s: &str) -> Result<(), String> {
        self.put(g, name, s.as_bytes(), 0)?;
        Ok(())
    }
}

/// Opens the results file for writing in the requested mode. `Ok(None)` means this process has
/// nothing to write (serial output chosen in a multi-process MPI run: only rank 0 writes).
fn open_io(path: &Path, mode: H5Mode) -> Result<Option<Io>, String> {
    match mode {
        H5Mode::Serial => {
            #[cfg(feature = "phdf5")]
            if let Some(w) = mpi_rt::world(false)? {
                if w.rank != 0 {
                    return Ok(None);
                }
            }
            Ok(Some(Io { file: h5::File::create(path).map_err(e("create file"))?, part: Partition::SOLO, replay: false }))
        }
        H5Mode::Parallel => {
            #[cfg(feature = "phdf5")]
            {
                let w = mpi_rt::world(true)?.ok_or("MPI is not available")?;
                let file = h5::File::with_options().with_fapl(|p| p.mpio(w.comm, None)).create(path).map_err(e("create file (MPI-IO)"))?;
                return Ok(Some(Io { file, part: Partition { rank: w.rank, size: w.size }, replay: false }));
            }
            #[cfg(not(feature = "phdf5"))]
            {
                let _ = path;
                Err("parallel HDF5 (MPI-IO) is not available in this build: it was made without the phdf5 feature (an MPI-enabled HDF5 library and an MPI SDK are needed to build it). Choose the serial HDF5 output".into())
            }
        }
    }
}

/// Writes the results file in the chosen mode, plus (on the process that writes the XML) the
/// XDMF description next to it when `xdmf` is set. Returns the paths written by this process.
pub fn write_outputs(path: &Path, mode: H5Mode, mesh: &SimplexMesh, spec: Option<&ProblemSpec>, result: &SolveResult, xdmf: bool) -> Result<Vec<std::path::PathBuf>, String> {
    let Some(io) = open_io(path, mode)? else { return Ok(vec![]) };
    let root = io.part.rank == 0;
    write_into(&io, mode, mesh, spec, result)?;
    drop(io);
    #[cfg(feature = "phdf5")]
    if mode == H5Mode::Parallel {
        if let Some(w) = mpi_rt::world(false)? {
            mpi_rt::barrier(&w);
        }
    }
    let mut out = vec![];
    if root {
        out.push(path.to_path_buf());
        if xdmf {
            out.push(write_xdmf(path, mesh, result)?);
        }
    }
    Ok(out)
}

pub fn write_results(path: &Path, mode: H5Mode, mesh: &SimplexMesh, spec: Option<&ProblemSpec>, result: &SolveResult) -> Result<(), String> {
    write_outputs(path, mode, mesh, spec, result, false).map(|_| ())
}

fn write_into(io: &Io, mode: H5Mode, mesh: &SimplexMesh, spec: Option<&ProblemSpec>, result: &SolveResult) -> Result<(), String> {
    let file = &io.file;
    let d = mesh.dim;
    io.attr_text(file, "format", "nyx-results")?;
    io.attr(file, "version", &2i32)?;
    io.attr_text(file, "writer", if mode == H5Mode::Serial { "serial" } else { "parallel" })?;
    io.attr(file, "processes", &(io.part.size as i64))?;
    if let Some(s) = spec {
        io.put_text(file, "spec", &serde_json::to_string(s).map_err(|x| x.to_string())?)?;
    }
    let m = io.group(file, "mesh")?;
    io.put(&m, "points", &mesh.points, d)?;
    io.put(&m, "cells", &mesh.cells, d + 1)?;
    io.put(&m, "cell_labels", mesh.cell_labels.as_slice(), 0)?;
    io.put(&m, "facets", &mesh.facets, d)?;
    io.put(&m, "facet_labels", mesh.facet_labels.as_slice(), 0)?;
    write_viz_mesh(io, mesh, result.order)?;
    let r = io.group(file, "results")?;
    io.attr(&r, "order", &(result.order as i64))?;
    io.attr(&r, "dim", &(result.dim as i64))?;
    io.attr(&r, "dofs", &(result.dofs as i64))?;
    io.attr(&r, "seconds", &result.seconds)?;
    io.attr(&r, "n_solutions", &(result.solutions.len() as i64))?;
    let scalars: Vec<serde_json::Value> = result
        .scalars
        .iter()
        .map(|(n, v)| match v {
            Scalar::Num(x) => serde_json::json!([n, x]),
            Scalar::Text(t) => serde_json::json!([n, t]),
        })
        .collect();
    io.put_text(&r, "scalars", &serde_json::Value::Array(scalars).to_string())?;
    io.put_text(&r, "log", &serde_json::to_string(&result.log).map_err(|x| x.to_string())?)?;
    for (k, s) in result.solutions.iter().enumerate() {
        let g = io.group(&r, &format!("solution_{k}"))?;
        io.attr_text(&g, "label", &s.label)?;
        if let Some(v) = s.value {
            io.attr(&g, "value", &v)?;
        }
        if let Some((re, im)) = s.neff {
            io.attr(&g, "neff_re", &re)?;
            io.attr(&g, "neff_im", &im)?;
        }
        let fg = io.group(&g, "fields")?;
        let order: Vec<&str> = s.fields.iter().map(|f| f.name.as_str()).collect();
        io.attr_text(&fg, "order", &serde_json::to_string(&order).map_err(|x| x.to_string())?)?;
        for f in &s.fields {
            let ds = io.put(&fg, f.name.as_str(), f.values.as_slice(), 0)?;
            io.attr_text(&ds, "unit", &f.unit)?;
        }
        if let Some(disp) = &s.disp {
            io.put(&g, "disp", disp, result.dim)?;
        }
        write_viz_solution(io, k, s, result.dim)?;
    }
    Ok(())
}

/// Name usable in a plain HDF5/XDMF path: letters, digits and underscores.
pub fn sanitize(name: &str) -> String {
    let s: String = name.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
    if s.chars().next().map_or(true, |c| c.is_ascii_digit()) { format!("f_{s}") } else { s }
}

/// XDMF topology type and the local node order of the cells for a mesh of the given order.
/// Quadratic cells use the VTK/XDMF numbering (corners, then edge midpoints), so the nodes are
/// the dofs of the P2 field and the attributes are nodal.
fn viz_topology(dim: usize, order: usize) -> (&'static str, Vec<usize>) {
    match (dim, order) {
        (2, 1) => ("Triangle", vec![0, 1, 2]),
        (3, 1) => ("Tetrahedron", vec![0, 1, 2, 3]),
        // ours: v0 v1 v2 e01 e02 e12  ->  v0 v1 v2 e01 e12 e02
        (2, _) => ("Triangle_6", vec![0, 1, 2, 3, 5, 4]),
        // ours: v0..v3 e01 e02 e03 e12 e13 e23  ->  e01 e12 e02 e03 e13 e23
        (_, _) => ("Tetrahedron_10", vec![0, 1, 2, 3, 4, 7, 5, 6, 8, 9]),
    }
}

fn write_viz_mesh(io: &Io, mesh: &SimplexMesh, order: usize) -> Result<(), String> {
    let order = if order == 2 { 2 } else { 1 };
    let viz = io.group(&io.file, "viz")?;
    let (topo, perm) = viz_topology(mesh.dim, order);
    let mut conn: Vec<u32> = Vec::with_capacity(mesh.n_cells() * perm.len());
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        conn.extend(perm.iter().map(|&i| dofs[i] as u32));
    }
    let points = mesh.dof_coordinates(order);
    io.put(&viz, "points", &points, mesh.dim)?;
    io.put(&viz, "cells", &conn, perm.len())?;
    io.attr_text(&viz, "topology", topo)?;
    Ok(())
}

fn write_viz_solution(io: &Io, k: usize, s: &Solution, dim: usize) -> Result<(), String> {
    let viz = io.file.group("viz").map_err(e("viz group"))?;
    let g = io.group(&viz, &format!("sol_{k}"))?;
    for f in &s.fields {
        io.link(&g, &format!("/results/solution_{k}/fields/{}", f.name), &sanitize(&f.name))?;
    }
    if let Some(d) = &s.disp {
        // vectors in XDMF have three components
        let n = d.len() / dim.max(1);
        let mut v3 = vec![0.0; n * 3];
        for i in 0..n {
            for c in 0..dim {
                v3[i * 3 + c] = d[i * dim + c];
            }
        }
        io.put(&g, "disp3", &v3, 3)?;
    }
    Ok(())
}

/// Writes the XDMF description (XML, `quick-xml`) of the results in `h5_path` next to it
/// (same name, extension `.xdmf`). The heavy data stay in the HDF5 file; ParaView and VisIt open
/// the `.xdmf`. Every solution (time step, mode) is one grid of a temporal collection.
pub fn write_xdmf(h5_path: &Path, mesh: &SimplexMesh, result: &SolveResult) -> Result<std::path::PathBuf, String> {
    use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
    use quick_xml::Writer;
    let order = if result.order == 2 { 2 } else { 1 };
    let (topo, perm) = viz_topology(mesh.dim, order);
    let np = mesh.n_dofs(order);
    let nc = mesh.n_cells();
    let h5name = h5_path.file_name().and_then(|n| n.to_str()).ok_or("bad file name")?.to_string();
    let xe = |x: std::io::Error| x.to_string();
    let mut w = Writer::new_with_indent(Vec::new(), b' ', 2);
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None))).map_err(xe)?;
    w.write_event(Event::DocType(BytesText::from_escaped("Xdmf SYSTEM \"Xdmf.dtd\" []"))).map_err(xe)?;
    let mut x = BytesStart::new("Xdmf");
    x.push_attribute(("Version", "3.0"));
    w.write_event(Event::Start(x)).map_err(xe)?;
    w.write_event(Event::Start(BytesStart::new("Domain"))).map_err(xe)?;
    let item = |w: &mut Writer<Vec<u8>>, dims: &str, number: &str, prec: &str, path: &str| -> Result<(), String> {
        let mut d = BytesStart::new("DataItem");
        d.push_attribute(("Dimensions", dims));
        d.push_attribute(("NumberType", number));
        d.push_attribute(("Precision", prec));
        d.push_attribute(("Format", "HDF"));
        w.write_event(Event::Start(d)).map_err(|x| x.to_string())?;
        w.write_event(Event::Text(BytesText::new(&format!("{h5name}:{path}")))).map_err(|x| x.to_string())?;
        w.write_event(Event::End(BytesEnd::new("DataItem"))).map_err(|x| x.to_string())
    };
    let times_ok = result.solutions.len() > 1 && result.solutions.iter().all(|s| s.value.is_some()) && result.solutions.windows(2).all(|p| p[1].value >= p[0].value) && result.solutions[0].neff.is_none() && result.solutions.iter().any(|s| s.label.starts_with("t ="));
    let collection = result.solutions.len() > 1;
    if collection {
        let mut c = BytesStart::new("Grid");
        c.push_attribute(("Name", "results"));
        c.push_attribute(("GridType", "Collection"));
        c.push_attribute(("CollectionType", "Temporal"));
        w.write_event(Event::Start(c)).map_err(xe)?;
    }
    for (k, s) in result.solutions.iter().enumerate() {
        let mut g = BytesStart::new("Grid");
        g.push_attribute(("Name", format!("solution_{k}").as_str()));
        g.push_attribute(("GridType", "Uniform"));
        w.write_event(Event::Start(g)).map_err(xe)?;
        if collection {
            let mut t = BytesStart::new("Time");
            let tv = if times_ok { s.value.unwrap_or(k as f64) } else { k as f64 };
            t.push_attribute(("Value", format!("{tv:e}").as_str()));
            w.write_event(Event::Empty(t)).map_err(xe)?;
        }
        let mut t = BytesStart::new("Topology");
        t.push_attribute(("TopologyType", topo));
        t.push_attribute(("NumberOfElements", nc.to_string().as_str()));
        w.write_event(Event::Start(t)).map_err(xe)?;
        item(&mut w, &format!("{nc} {}", perm.len()), "UInt", "4", "/viz/cells")?;
        w.write_event(Event::End(BytesEnd::new("Topology"))).map_err(xe)?;
        let mut ge = BytesStart::new("Geometry");
        ge.push_attribute(("GeometryType", if mesh.dim == 2 { "XY" } else { "XYZ" }));
        w.write_event(Event::Start(ge)).map_err(xe)?;
        item(&mut w, &format!("{np} {}", mesh.dim), "Float", "8", "/viz/points")?;
        w.write_event(Event::End(BytesEnd::new("Geometry"))).map_err(xe)?;
        for f in &s.fields {
            let mut a = BytesStart::new("Attribute");
            a.push_attribute(("Name", f.name.as_str()));
            a.push_attribute(("AttributeType", "Scalar"));
            a.push_attribute(("Center", "Node"));
            w.write_event(Event::Start(a)).map_err(xe)?;
            item(&mut w, &np.to_string(), "Float", "8", &format!("/viz/sol_{k}/{}", sanitize(&f.name)))?;
            w.write_event(Event::End(BytesEnd::new("Attribute"))).map_err(xe)?;
        }
        if s.disp.is_some() {
            let mut a = BytesStart::new("Attribute");
            a.push_attribute(("Name", "displacement"));
            a.push_attribute(("AttributeType", "Vector"));
            a.push_attribute(("Center", "Node"));
            w.write_event(Event::Start(a)).map_err(xe)?;
            item(&mut w, &format!("{np} 3"), "Float", "8", &format!("/viz/sol_{k}/disp3"))?;
            w.write_event(Event::End(BytesEnd::new("Attribute"))).map_err(xe)?;
        }
        w.write_event(Event::End(BytesEnd::new("Grid"))).map_err(xe)?;
    }
    if collection {
        w.write_event(Event::End(BytesEnd::new("Grid"))).map_err(xe)?;
    }
    w.write_event(Event::End(BytesEnd::new("Domain"))).map_err(xe)?;
    w.write_event(Event::End(BytesEnd::new("Xdmf"))).map_err(xe)?;
    let xml = String::from_utf8(w.into_inner()).map_err(|x| x.to_string())?;
    let out = h5_path.with_extension("xdmf");
    std::fs::write(&out, xml).map_err(|x| format!("{}: {x}", out.display()))?;
    Ok(out)
}

fn read_text(g: &h5::Group, name: &str) -> Result<String, String> {
    let ds = g.dataset(name).map_err(e(name))?;
    // version 2 files hold the text as UTF-8 bytes, version 1 as one variable-length string
    if let Ok(bytes) = ds.read_raw::<u8>() {
        return String::from_utf8(bytes).map_err(|x| format!("HDF5 {name}: {x}"));
    }
    let v = ds.read_raw::<VarLenUnicode>().map_err(e(name))?;
    Ok(v.first().map(|x| x.as_str().to_string()).unwrap_or_default())
}

fn attr_f64(g: &h5::Group, name: &str) -> Option<f64> {
    g.attr(name).ok().and_then(|a| a.read_scalar::<f64>().ok())
}

/// Reads a results file back: the mesh, the stored problem description (JSON) and the result.
pub fn read_results(path: &Path) -> Result<(SimplexMesh, Option<String>, SolveResult), String> {
    let file = h5::File::open(path).map_err(e("open file"))?;
    let fmt = file.attr("format").ok().and_then(|a| a.read_scalar::<VarLenUnicode>().ok()).map(|v| v.as_str().to_string());
    if fmt.as_deref() != Some("nyx-results") {
        return Err("this HDF5 file was not written by nyx (no format attribute)".into());
    }
    let spec = read_text(&file, "spec").ok();
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
    let r = file.group("results").map_err(e("results group"))?;
    let geti = |n: &str| r.attr(n).ok().and_then(|a| a.read_scalar::<i64>().ok()).unwrap_or(0) as usize;
    let nsol = geti("n_solutions");
    let mut result = SolveResult { order: geti("order"), dim: geti("dim"), dofs: geti("dofs"), seconds: attr_f64(&r, "seconds").unwrap_or(0.0), ..Default::default() };
    if let Ok(t) = read_text(&r, "scalars") {
        if let Ok(serde_json::Value::Array(items)) = serde_json::from_str::<serde_json::Value>(&t) {
            for it in items {
                if let (Some(n), Some(v)) = (it.get(0).and_then(|x| x.as_str()), it.get(1)) {
                    result.scalars.push((n.to_string(), match v.as_f64() {
                        Some(x) => Scalar::Num(x),
                        None => Scalar::Text(v.as_str().unwrap_or("").to_string()),
                    }));
                }
            }
        }
    }
    if let Ok(t) = read_text(&r, "log") {
        result.log = serde_json::from_str(&t).unwrap_or_default();
    }
    for k in 0..nsol {
        let g = r.group(&format!("solution_{k}")).map_err(e("solution"))?;
        let mut s = Solution::default();
        s.label = g.attr("label").ok().and_then(|a| a.read_scalar::<VarLenUnicode>().ok()).map(|v| v.as_str().to_string()).unwrap_or_default();
        s.value = attr_f64(&g, "value");
        if let (Some(re), Some(im)) = (attr_f64(&g, "neff_re"), attr_f64(&g, "neff_im")) {
            s.neff = Some((re, im));
        }
        let fg = g.group("fields").map_err(e("fields"))?;
        let mut names = fg.member_names().map_err(e("fields"))?;
        names.sort();
        if let Some(o) = fg.attr("order").ok().and_then(|a| a.read_scalar::<VarLenUnicode>().ok()).and_then(|v| serde_json::from_str::<Vec<String>>(v.as_str()).ok()) {
            if o.len() == names.len() {
                names = o;
            }
        }
        for name in names {
            let ds = fg.dataset(&name).map_err(e(&name))?;
            let unit = ds.attr("unit").ok().and_then(|a| a.read_scalar::<VarLenUnicode>().ok()).map(|v| v.as_str().to_string()).unwrap_or_default();
            s.fields.push(Field { name, unit, values: ds.read_raw::<f64>().map_err(e("field"))? });
        }
        if let Ok(ds) = g.dataset("disp") {
            s.disp = Some(ds.read_raw::<f64>().map_err(e("displacement"))?);
        }
        result.solutions.push(s);
    }
    Ok((mesh, spec, result))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let mesh = crate::geomdrive::plate(1.0, 0.5, 0.2).unwrap().mesh;
        let mut spec = ProblemSpec::new(Physics::Heat);
        spec.materials.default.insert("k".into(), 2.0);
        spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(300.0) }];
        let r = crate::dfs::solve_heat(&mesh, &spec, &mut SolveState::default()).unwrap();
        let path = std::env::temp_dir().join(format!("nyx-h5-test-{}.h5", std::process::id()));
        write_results(&path, H5Mode::Serial, &mesh, Some(&spec), &r).unwrap();
        let (m2, spec2, r2) = read_results(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!((m2.points.clone(), m2.cells.clone(), m2.cell_labels.clone(), m2.facets.clone(), m2.facet_labels.clone()), (mesh.points.clone(), mesh.cells.clone(), mesh.cell_labels.clone(), mesh.facets.clone(), mesh.facet_labels.clone()));
        assert_eq!(serde_json::from_str::<ProblemSpec>(&spec2.unwrap()).unwrap(), spec);
        assert_eq!(r2.solutions.len(), r.solutions.len());
        let (na, nb): (Vec<&str>, Vec<&str>) = (r.solutions[0].fields.iter().map(|f| f.name.as_str()).collect(), r2.solutions[0].fields.iter().map(|f| f.name.as_str()).collect());
        assert_eq!(na, nb, "field order is kept");
        let ta = &r.solutions[0].field("T").unwrap().values;
        let tb = &r2.solutions[0].field("T").unwrap().values;
        assert_eq!(ta, tb);
        assert_eq!(r.scalars.len(), r2.scalars.len());
        assert_eq!(r2.order, r.order);
    }

    #[test]
    fn xdmf_describes_the_h5_data_consistently() {
        // quadratic triangles and tetrahedra: the midside nodes must sit at the midpoints the VTK/XDMF
        // numbering says, and every dimension the XML states must match the HDF5 datasets
        for (mesh, name) in [(crate::geomdrive::plate(1.0, 0.5, 0.4).unwrap().mesh, "tri"), (crate::geomdrive::block(1.0, 0.5, 0.4, 0.5, 1).unwrap().mesh, "tet")] {
            let mut spec = ProblemSpec::new(if mesh.dim == 2 { Physics::Solid } else { Physics::Heat });
            let r = if mesh.dim == 2 {
                spec.materials.default.insert("E".into(), 1e9);
                spec.materials.default.insert("nu".into(), 0.3);
                spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![4], kind: BcKind::Load([1e5, 0.0, 0.0]) }];
                crate::sms::solve_solid(&mesh, &spec, &SolveState::default()).unwrap()
            } else {
                spec.materials.default.insert("k".into(), 1.0);
                spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(1.0) }, Bc { on: None, labels: vec![6], kind: BcKind::Temperature(0.0) }];
                crate::dfs::solve_heat(&mesh, &spec, &mut SolveState::default()).unwrap()
            };
            let path = std::env::temp_dir().join(format!("nyx-xdmf-{name}-{}.h5", std::process::id()));
            write_results(&path, H5Mode::Serial, &mesh, Some(&spec), &r).unwrap();
            let x = write_xdmf(&path, &mesh, &r).unwrap();
            let xml = std::fs::read_to_string(&x).unwrap();
            let file = h5::File::open(&path).unwrap();
            let pts = file.dataset("viz/points").unwrap().read_raw::<f64>().unwrap();
            let conn = file.dataset("viz/cells").unwrap().read_raw::<u32>().unwrap();
            let d = mesh.dim;
            let per = if d == 2 { 6 } else { 10 };
            assert_eq!(conn.len(), mesh.n_cells() * per);
            assert_eq!(pts.len(), mesh.n_dofs(2) * d);
            let at = |i: u32| &pts[i as usize * d..(i as usize + 1) * d];
            let mid = |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).map(|(p, q)| 0.5 * (p + q)).collect() };
            // (midside node index, corner a, corner b) per the VTK numbering
            let edges: Vec<(usize, usize, usize)> = if d == 2 { vec![(3, 0, 1), (4, 1, 2), (5, 2, 0)] } else { vec![(4, 0, 1), (5, 1, 2), (6, 2, 0), (7, 0, 3), (8, 1, 3), (9, 2, 3)] };
            for cell in conn.chunks(per) {
                for &(m, a, b) in &edges {
                    let want = mid(at(cell[a]), at(cell[b]));
                    let got = at(cell[m]);
                    assert!(want.iter().zip(got).all(|(p, q)| (p - q).abs() < 1e-12), "{name}: midside node {m} is not the midpoint of {a} and {b}");
                }
            }
            assert!(xml.contains(if d == 2 { "Triangle_6" } else { "Tetrahedron_10" }));
            assert!(xml.contains(&format!("NumberOfElements=\"{}\"", mesh.n_cells())));
            assert!(xml.contains(&format!("Dimensions=\"{} {per}\"", mesh.n_cells())));
            assert!(xml.contains(&format!("Dimensions=\"{} {d}\"", mesh.n_dofs(2))));
            // every attribute path exists in the file
            for line in xml.lines().filter(|l| l.contains(".h5:/viz/sol_")) {
                let p = line.trim().trim_start_matches(|c| c != '/').split('<').next().unwrap().to_string();
                assert!(file.dataset(&p).is_ok(), "missing {p}");
            }
            // the XML parses as XML
            let mut rd = quick_xml::Reader::from_str(&xml);
            let mut depth = 0i32;
            loop {
                match rd.read_event().unwrap() {
                    quick_xml::events::Event::Start(_) => depth += 1,
                    quick_xml::events::Event::End(_) => depth -= 1,
                    quick_xml::events::Event::Eof => break,
                    _ => {}
                }
            }
            assert_eq!(depth, 0);
            let _ = std::fs::remove_file(&path);
            let _ = std::fs::remove_file(&x);
        }
    }

    #[test]
    fn row_blocks_tile_every_dataset() {
        for size in 1..=9 {
            for n in [0, 1, 2, 7, 8, 100, 101] {
                let mut next = 0;
                let (mut lo, mut hi) = (usize::MAX, 0);
                for rank in 0..size {
                    let r = Partition { rank, size }.rows(n);
                    assert_eq!(r.start, next, "blocks are contiguous (size {size}, n {n})");
                    next = r.end;
                    lo = lo.min(r.len());
                    hi = hi.max(r.len());
                }
                assert_eq!(next, n, "blocks cover all rows (size {size}, n {n})");
                assert!(hi - lo <= 1, "blocks are balanced (size {size}, n {n})");
            }
        }
    }

    /// The parallel writer cannot run here (no MPI), but its logic can: three emulated processes each
    /// write only their own rows into one file, which must equal what one serial process writes.
    #[test]
    fn partitioned_writes_equal_a_serial_write() {
        let b = crate::geomdrive::block(1.0, 0.5, 0.4, 0.34, 2).unwrap();
        let mesh = b.mesh.clone();
        let mut spec = ProblemSpec::new(Physics::Heat);
        spec.materials.default.insert("k".into(), 1.0);
        spec.bcs = vec![Bc { on: None, labels: b.xmin.clone(), kind: BcKind::Temperature(1.0) }, Bc { on: None, labels: b.xmax.clone(), kind: BcKind::Temperature(0.0) }];
        let r = crate::dfs::solve_heat(&mesh, &spec, &mut SolveState::default()).unwrap();
        let tmp = |n: &str| std::env::temp_dir().join(format!("nyx-part-{n}-{}.h5", std::process::id()));
        let (pa, pb) = (tmp("serial"), tmp("three"));
        write_results(&pa, H5Mode::Serial, &mesh, Some(&spec), &r).unwrap();
        {
            let file = h5::File::create(&pb).unwrap();
            for rank in 0..3 {
                let io = Io { file: file.clone(), part: Partition { rank, size: 3 }, replay: rank > 0 };
                write_into(&io, H5Mode::Parallel, &mesh, Some(&spec), &r).unwrap();
            }
        }
        let (fa, fb) = (h5::File::open(&pa).unwrap(), h5::File::open(&pb).unwrap());
        for path in ["mesh/points", "mesh/cells", "mesh/facets", "viz/points", "viz/cells"] {
            let (a, b) = (fa.dataset(path).unwrap(), fb.dataset(path).unwrap());
            assert_eq!(a.shape(), b.shape(), "{path}");
            if path.ends_with("cells") || path.ends_with("facets") {
                assert_eq!(a.read_raw::<u32>().unwrap(), b.read_raw::<u32>().unwrap(), "{path}");
            } else {
                assert_eq!(a.read_raw::<f64>().unwrap(), b.read_raw::<f64>().unwrap(), "{path}");
            }
        }
        let (ra, rb) = (read_results(&pa).unwrap(), read_results(&pb).unwrap());
        assert_eq!(ra.1, rb.1, "problem description");
        assert_eq!(ra.2.solutions[0].field("T").unwrap().values, rb.2.solutions[0].field("T").unwrap().values);
        assert_eq!(ra.2.scalars.len(), rb.2.scalars.len());
        let _ = std::fs::remove_file(&pa);
        let _ = std::fs::remove_file(&pb);
    }

    #[test]
    fn serial_output_is_readable_from_the_c_library_layout() {
        // text is stored as UTF-8 bytes: the JSON of the spec must parse straight from the dataset
        let mesh = crate::geomdrive::plate(1.0, 1.0, 0.5).unwrap().mesh;
        let mut spec = ProblemSpec::new(Physics::Heat);
        spec.materials.default.insert("k".into(), 1.0);
        let path = std::env::temp_dir().join(format!("nyx-bytes-{}.h5", std::process::id()));
        write_results(&path, H5Mode::Serial, &mesh, Some(&spec), &SolveResult::default()).unwrap();
        let bytes = h5::File::open(&path).unwrap().dataset("spec").unwrap().read_raw::<u8>().unwrap();
        assert_eq!(serde_json::from_slice::<ProblemSpec>(&bytes).unwrap(), spec);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn parallel_mode_reports_what_is_missing() {
        if !mode_available(H5Mode::Parallel) {
            let mesh = crate::geomdrive::plate(1.0, 1.0, 0.5).unwrap().mesh;
            let r = SolveResult::default();
            let err = write_results(&std::env::temp_dir().join("never.h5"), H5Mode::Parallel, &mesh, None, &r).unwrap_err();
            assert!(err.contains("MPI"));
        }
    }
}
