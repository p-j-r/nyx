//! Writing and reading the results file as engine commands. The HDF5 writer lives in `h5io` (cargo feature `hdf5`); a build without the
//! feature answers these commands with a plain refusal, so a host can tell.

use crate::engine::Engine;
use nyx_physics::Error;
#[allow(unused_imports)]
use serde_json::{json, Value};

impl Engine {
    /// Writes the project, the mesh and the latest result of every solved physics into an HDF5 file. With `parallel` the file is written with
    /// the MPI-IO driver by all the processes of an MPI run, each calling this with the same arguments (each writes its own block of rows; the
    /// XDMF file is written by process 0); only a build with the `phdf5` feature can, any other answers with what is missing.
    pub fn export_results(&self, path: &str, parallel: bool) -> Result<Value, Error> {
        #[cfg(feature = "hdf5")]
        {
            let model = self.model().ok_or_else(|| Error::invalid("there is no mesh; solve first"))?;
            let mut results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)> = vec![];
            for s in &self.project.physics {
                if let Some((r, _)) = self.result(&s.physics) {
                    results.push((s.physics.clone(), r.clone()));
                }
            }
            if results.is_empty() {
                return Err(Error::invalid("there are no results yet; solve first"));
            }
            let refs: Vec<(String, &nyx_physics::SolveResult)> = results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();
            let p = std::path::Path::new(path);
            let mesh = model.solver_mesh();
            let mode = if parallel { crate::h5io::H5Mode::Parallel } else { crate::h5io::H5Mode::Serial };
            let root = crate::h5io::write_results_mode(p, mode, &self.project, &mesh, &refs).map_err(Error::invalid)?;
            let physics = results.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>();
            if !root {
                // another process of an MPI run writes the XDMF file and reports
                return Ok(json!({ "path": path, "writer": "parallel", "written_by_this_process": false, "physics": physics }));
            }
            let xdmf = crate::xdmf::write_xdmf(p, &mesh, &refs).map_err(Error::invalid)?;
            Ok(json!({ "path": path, "xdmf": xdmf, "bytes": crate::h5io::file_bytes(p), "physics": physics, "writer": if parallel { "parallel" } else { "serial" } }))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = (path, parallel);
            let _ = path;
            Err(Error::Unsupported("this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature; the HDF5 library must be installed)".into()))
        }
    }

    /// Opens a results file: its project becomes the project, the project is meshed again and the stored results become the latest results
    /// (not stale), so plots, analyses and exports work on them as after a solve. Meshing is deterministic (ADR 55), so the mesh made again
    /// is the stored one; when it is not (a file of another version of nyx) nothing is changed and the message says so.
    pub fn open_results(&mut self, path: &str) -> Result<Value, Error> {
        #[cfg(feature = "hdf5")]
        {
            let stored = crate::h5io::read_results(std::path::Path::new(path)).map_err(Error::invalid)?;
            let project = crate::project::Project::from_json(&stored.project)?;
            project.check(&self.registry)?;
            let before = (self.project.clone(), self.model.clone(), std::mem::take(&mut self.results), self.link.take());
            self.set_project(project);
            let restore = |e: &mut Engine, b: (crate::project::Project, _, _, _)| {
                e.project = b.0;
                e.model = b.1;
                e.results = b.2;
                e.link = b.3;
            };
            if let Err(e) = self.mesh_now() {
                restore(self, before);
                return Err(e);
            }
            let made = self.model().expect("the mesh was just built").solver_mesh();
            let s = &stored.mesh;
            if made.dim != s.dim || made.points != s.points || made.cells != s.cells || made.cell_labels != s.cell_labels || made.facets != s.facets || made.facet_labels != s.facet_labels {
                restore(self, before);
                return Err(Error::invalid(format!(
                    "{path}: the project of this file meshes to something else than the mesh stored in it ({} cells made, {} stored); it was written by another version of nyx, so its results cannot be shown on the geometry",
                    made.n_cells(),
                    s.n_cells()
                )));
            }
            let stamp = self.stamp();
            let mut physics = vec![];
            for (code, r) in stored.results {
                physics.push(json!({
                    "physics": code, "study": r.study, "dofs": r.dofs, "seconds": r.seconds, "solutions": r.solutions.len(), "scalars": r.scalars,
                }));
                self.results.insert(code, crate::engine::Stored { result: std::sync::Arc::new(r), stamp: stamp.clone() });
            }
            Ok(json!({ "title": self.project.title, "path": path, "cells": made.n_cells(), "results": physics }))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = path;
            Err(Error::Unsupported("this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature)".into()))
        }
    }

    /// What a results file holds: the project title, the mesh size and, per physics, the study, the solutions and the scalars.
    pub fn results_file_info(&self, path: &str) -> Result<Value, Error> {
        #[cfg(feature = "hdf5")]
        {
            let s = crate::h5io::read_results(std::path::Path::new(path)).map_err(Error::invalid)?;
            let title = serde_json::from_str::<Value>(&s.project).ok().and_then(|v| v.get("title").cloned()).unwrap_or(Value::Null);
            Ok(json!({
                "title": title, "dim": s.mesh.dim, "cells": s.mesh.n_cells(), "points": s.mesh.n_points(),
                "results": s.results.iter().map(|(c, r)| json!({
                    "physics": c, "study": r.study, "order": r.order, "dofs": r.dofs, "seconds": r.seconds,
                    "solutions": r.solutions.len(), "scalars": r.scalars,
                })).collect::<Vec<_>>(),
            }))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = path;
            Err(Error::Unsupported("this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature)".into()))
        }
    }
}
