// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Writing and reading the results file as engine commands. The HDF5 writer lives in `h5io` (cargo feature `hdf5`); a build without the
//! feature answers these commands with a plain refusal, so a host can tell.

use crate::engine::Engine;
use nyx_physics::Error;
#[allow(unused_imports)]
use serde_json::{json, Value};

/// What a results file is written from: the project, the mesh and the latest result of every physics, as they were when the write was asked for.
/// Owned and `Send`, so the write can run on a job while the engine goes on.
pub struct Snapshot {
    pub project: crate::project::Project,
    pub model: std::sync::Arc<crate::meshing::MeshedModel>,
    pub results: Vec<(String, std::sync::Arc<nyx_physics::SolveResult>)>,
}

/// Writes `snapshot` into the HDF5 file `path` and its XDMF description. With `atomic` the file is written under another name and renamed when it
/// is complete, so an interrupted or failed write never leaves a file that looks finished (not with `parallel`: all processes write one file).
pub fn write_snapshot(snapshot: &Snapshot, path: &str, parallel: bool, atomic: bool) -> Result<Value, Error> {
    #[cfg(feature = "hdf5")]
    {
        let refs: Vec<(String, &nyx_physics::SolveResult)> =
            snapshot.results.iter().map(|(c, r)| (c.clone(), r.as_ref())).collect();
        let p = std::path::Path::new(path);
        let mesh = snapshot.model.solver_mesh();
        let mode = if parallel { crate::h5io::H5Mode::Parallel } else { crate::h5io::H5Mode::Serial };
        let temporary = std::path::PathBuf::from(format!("{path}.partial"));
        let atomic = atomic && !parallel;
        let target = if atomic { temporary.as_path() } else { p };
        let written =
            crate::h5io::write_results_mode(target, mode, &snapshot.project, &mesh, &refs).map_err(Error::invalid);
        let root = match written {
            Ok(r) => r,
            Err(e) => {
                if atomic {
                    let _ = std::fs::remove_file(&temporary);
                }
                return Err(e);
            }
        };
        if atomic {
            std::fs::rename(&temporary, p).map_err(|e| {
                let _ = std::fs::remove_file(&temporary);
                Error::invalid(format!("{path}: {e}"))
            })?;
        }
        let physics = snapshot.results.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>();
        if !root {
            // another process of an MPI run writes the XDMF file and reports
            return Ok(
                json!({ "path": path, "writer": "parallel", "written_by_this_process": false, "physics": physics }),
            );
        }
        let xdmf = crate::xdmf::write_xdmf(p, &mesh, &refs).map_err(Error::invalid)?;
        Ok(
            json!({ "path": path, "xdmf": xdmf, "bytes": crate::h5io::file_bytes(p), "physics": physics, "writer": if parallel { "parallel" } else { "serial" } }),
        )
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (snapshot, path, parallel, atomic);
        Err(Error::Unsupported(
            "this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature)".into(),
        ))
    }
}

/// Reads the results file `path`, meshes its project again on the calling thread (under `progress`) and checks that the mesh is the stored one.
/// Touches no engine: the outcome is applied by `Engine::poll` (`Opened`).
pub(crate) fn read_and_mesh(
    registry: &nyx_physics::Registry,
    mesher: &std::sync::Arc<dyn crate::meshing::Mesher>,
    library: &std::sync::Arc<femcore::material::Library>,
    path: &str,
    progress: &std::sync::Arc<femcore::progress::Progress>,
) -> Result<crate::jobs::Opened, Error> {
    #[cfg(feature = "hdf5")]
    {
        let stored = {
            let _stage = femcore::progress::stage("reading the results file");
            crate::h5io::read_results(std::path::Path::new(path)).map_err(Error::invalid)?
        };
        let project = crate::project::Project::from_json(&stored.project)?;
        project.check(registry)?;
        project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        if progress.is_cancelled() {
            return Err(Error::Cancelled);
        }
        let mut scratch = Engine::with(registry.clone(), mesher.clone());
        scratch.set_base_library_arc(library.clone());
        scratch.set_project(project.clone());
        let sizing = scratch.sizing();
        let model = {
            let _installed = femcore::progress::install(progress);
            let _stage = femcore::progress::stage("meshing the project of the file again");
            std::sync::Arc::new(mesher.mesh(&project.geometry, &project.mesh, sizing.as_ref())?)
        };
        let model = if !stored.fields_saved {
            model
        } else {
            let made = model.solver_mesh();
            let s = &stored.mesh;
            if made.dim != s.dim
                || made.points != s.points
                || made.cells != s.cells
                || made.cell_labels != s.cell_labels
                || made.facets != s.facets
                || made.facet_labels != s.facet_labels
                || made.edge_node_list() != s.edge_node_list()
                || made.qp_edge_list() != s.qp_edge_list()
            {
                if s.n_cells() > 0 && crate::amr::project_has_electro_amr(&project) && made.dim == s.dim {
                    std::sync::Arc::new(
                        crate::meshing::MeshedModel::new(
                            stored.mesh,
                            model.evaluated().clone(),
                            model.labels().clone(),
                            project.mesh.sweep.is_some(),
                        )?
                        .with_amr_adapted(),
                    )
                } else {
                    return Err(Error::invalid(format!(
                        "{path}: the project of this file meshes to something else than the mesh stored in it ({} cells made, {} stored); it was written by another version of nyx, so its results cannot be shown on the geometry",
                        made.n_cells(),
                        s.n_cells()
                    )));
                }
            } else {
                model
            }
        };
        Ok(crate::jobs::Opened { project, model, results: stored.results, path: path.to_string() })
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (registry, mesher, library, path, progress);
        Err(Error::Unsupported(
            "this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature)".into(),
        ))
    }
}

impl Engine {
    /// The project, the mesh and the latest result of every solved physics, to be written to a file (`write_snapshot`).
    pub fn results_snapshot(&self) -> Result<Snapshot, Error> {
        let model = self.model().ok_or_else(|| Error::invalid("there is no mesh; solve first"))?.clone();
        let mut results = vec![];
        for s in &self.project.physics {
            if let Some((r, _)) = self.result(&s.physics) {
                results.push((s.physics.clone(), r.clone()));
            }
        }
        if results.is_empty() {
            return Err(Error::invalid("there are no results yet; solve first"));
        }
        Ok(Snapshot { project: self.project.clone(), model, results })
    }

    /// Writes the project, the mesh and the latest result of every solved physics into an HDF5 file, and waits for it. With `parallel` the file
    /// is written with the MPI-IO driver by all the processes of an MPI run, each calling this with the same arguments (each writes its own block
    /// of rows; the XDMF file is written by process 0); only a build with the `phdf5` feature can, any other answers with what is missing.
    /// A window or a host that must not wait uses `start_export_results`.
    pub fn export_results(&self, path: &str, parallel: bool) -> Result<Value, Error> {
        write_snapshot(&self.results_snapshot()?, path, parallel, false)
    }

    /// Opens a results file: its project becomes the project, the project is meshed again and the stored results become the latest results
    /// (not stale), so plots, analyses and exports work on them as after a solve. Meshing is deterministic (ADR 55), so the mesh made again
    /// is the stored one; when it is not (a file of another version of nyx) nothing is changed and the message says so. Waits for it; a window
    /// uses `start_open_results`.
    pub fn open_results(&mut self, path: &str) -> Result<Value, Error> {
        let id = self.start_open_results(path)?;
        self.wait(id)
    }

    /// What a results file holds: the project title, the mesh size and, per physics, the study, the solutions and the scalars.
    pub fn results_file_info(&self, path: &str) -> Result<Value, Error> {
        #[cfg(feature = "hdf5")]
        {
            let s = crate::h5io::read_results_summary(std::path::Path::new(path)).map_err(Error::invalid)?;
            let title = serde_json::from_str::<Value>(&s.project)
                .ok()
                .and_then(|v| v.get("title").cloned())
                .unwrap_or(Value::Null);
            Ok(json!({
                "title": title, "dim": s.dim, "cells": s.cells, "points": s.points,
                "results": s.results.iter().map(|(c, r)| json!({
                    "physics": c, "study": r.study, "order": r.order, "dofs": r.dofs, "seconds": r.seconds,
                    "solutions": r.solutions, "scalars": r.scalars,
                })).collect::<Vec<_>>(),
            }))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = path;
            Err(Error::Unsupported(
                "this build of nyx was made without HDF5 support (rebuild nyx-engine with the `hdf5` feature)".into(),
            ))
        }
    }
}
