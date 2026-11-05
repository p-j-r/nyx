//! Writing and reading the results file as engine commands. The HDF5 writer lives in `h5io` (cargo feature `hdf5`); a build without the
//! feature answers these commands with a plain refusal, so a host can tell.

use crate::engine::Engine;
use nyx_physics::Error;
#[allow(unused_imports)]
use serde_json::{json, Value};

impl Engine {
    /// Writes the project, the mesh and the latest result of every solved physics into an HDF5 file.
    pub fn export_results(&self, path: &str) -> Result<Value, Error> {
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
            crate::h5io::write_results(p, &self.project, &mesh, &refs).map_err(Error::invalid)?;
            let xdmf = crate::xdmf::write_xdmf(p, &mesh, &refs).map_err(Error::invalid)?;
            Ok(json!({ "path": path, "xdmf": xdmf, "bytes": crate::h5io::file_bytes(p), "physics": results.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>() }))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = path;
            Err(Error::Unsupported("this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature; the HDF5 library must be installed)".into()))
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
