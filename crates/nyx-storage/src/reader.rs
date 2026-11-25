// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Bounded-memory streaming HDF5 reader (`File::open_streaming`), mode slicer, and untrusted fuzz walker.

use crate::{assert_f64_bits_eq, MeshData, ModeHeader, PackedComponentMeta, ReadBudgetConfig};
use anyhow::{anyhow, bail, Result};
use hdf5_pure::{AttrValue, Dataset, File, Group};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// Extracts a UTF-8 string attribute, handling both scalar strings (`as_str()`) and 1-element string arrays (`as_strings()`).
pub fn attr_str(attrs: &HashMap<String, AttrValue>, key: &str) -> Result<String> {
    let v = attrs.get(key).ok_or_else(|| anyhow!("missing string attr '{}'", key))?;
    if let Some(s) = v.as_str() {
        return Ok(s.to_string());
    }
    if let Some(arr) = v.as_strings() {
        if let Some(first) = arr.first() {
            return Ok(first.clone());
        }
    }
    bail!("attr '{}' is not a string: {:?}", key, v)
}

/// Extracts an integer attribute as `i64` (accepting `I8`/`I16`/`I32`/`I64`/`U8`/`U16`/`U32`/`U64`).
pub fn attr_i64(attrs: &HashMap<String, AttrValue>, key: &str) -> Result<i64> {
    let v = attrs.get(key).ok_or_else(|| anyhow!("missing int attr '{}'", key))?;
    v.as_i64()
        .or_else(|| v.as_u64().map(|u| u as i64))
        .ok_or_else(|| anyhow!("attr '{}' is not an integer: {:?}", key, v))
}

/// Extracts a scalar `f64` attribute.
pub fn attr_f64(attrs: &HashMap<String, AttrValue>, key: &str) -> Result<f64> {
    let v = attrs.get(key).ok_or_else(|| anyhow!("missing f64 attr '{}'", key))?;
    v.as_f64().ok_or_else(|| anyhow!("attr '{}' is not an f64: {:?}", key, v))
}

/// Tuple `(name, unit, nloc, re, im)` for one legacy v3 cell-field component.
pub type LegacyCellFieldEntry = (String, String, i64, Vec<f64>, Vec<f64>);

/// Tuple `(name, unit, components_json, data)` for one legacy v3 nodal field.
pub type LegacyNodalFieldEntry = (String, String, String, Vec<f64>);

/// Per-solution data in a legacy `nyx-results` v3 file.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NyxSolutionData {
    pub eig_re: f64,
    pub eig_im: f64,
    pub label: String,
    /// Map of component group name (e.g. `"00_Ex"`) -> `(name, unit, nloc, re, im)`.
    pub cell_fields: BTreeMap<String, LegacyCellFieldEntry>,
    pub fields_order_json: String,
    /// Map of field dataset name (e.g. `"00_|E|"`) -> `(name, unit, components, data)`.
    pub fields: BTreeMap<String, LegacyNodalFieldEntry>,
    /// Map of `/viz/<phys>/sol_<k>/<name>` dataset -> data.
    pub viz_fields: BTreeMap<String, Vec<f64>>,
}

/// Complete representation of a `nyx-results` v3 file for Gate 2 (`EC-01`) & Gate 3 (`EC-02`) 0-ULP verification.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NyxV3Summary {
    pub format: String,
    pub version: i64,
    pub writer: String,
    pub processes: i64,
    pub project_bytes: Vec<u8>,
    pub mesh_cell_labels: Vec<i32>,
    pub mesh_cells_shape: [usize; 2],
    pub mesh_cells: Vec<u32>,
    pub mesh_facet_labels: Vec<i32>,
    pub mesh_facets_shape: [usize; 2],
    pub mesh_facets: Vec<u32>,
    pub mesh_points_shape: [usize; 2],
    pub mesh_points: Vec<f64>,
    pub physics: String,
    pub study: String,
    pub dim: i64,
    pub dofs: i64,
    pub n_solutions: i64,
    pub order: i64,
    pub seconds: f64,
    pub log_bytes: Vec<u8>,
    pub scalars_bytes: Vec<u8>,
    pub viz_topology: String,
    pub viz_cells_shape: [usize; 2],
    pub viz_cells: Vec<u32>,
    pub viz_points_shape: [usize; 2],
    pub viz_points: Vec<f64>,
    pub solutions: Vec<NyxSolutionData>,
}

impl NyxV3Summary {
    /// Verifies exact 0-ULP (`to_bits()`) bitwise equality against another [`NyxV3Summary`].
    pub fn verify_0_ulp(&self, expected: &Self, check_viz_links: bool) -> Result<()> {
        if self.format != expected.format
            || self.version != expected.version
            || self.writer != expected.writer
            || self.processes != expected.processes
        {
            bail!(
                "root attrs mismatch: actual=({},{},{},{}) vs expected=({},{},{},{})",
                self.format,
                self.version,
                self.writer,
                self.processes,
                expected.format,
                expected.version,
                expected.writer,
                expected.processes
            );
        }
        if self.project_bytes != expected.project_bytes {
            bail!("project bytes mismatch");
        }
        if self.mesh_cell_labels != expected.mesh_cell_labels {
            bail!("mesh/cell_labels mismatch");
        }
        if self.mesh_cells_shape != expected.mesh_cells_shape || self.mesh_cells != expected.mesh_cells {
            bail!("mesh/cells mismatch");
        }
        if self.mesh_facet_labels != expected.mesh_facet_labels {
            bail!("mesh/facet_labels mismatch");
        }
        if self.mesh_facets_shape != expected.mesh_facets_shape || self.mesh_facets != expected.mesh_facets {
            bail!("mesh/facets mismatch");
        }
        if self.mesh_points_shape != expected.mesh_points_shape {
            bail!("mesh/points shape mismatch");
        }
        assert_f64_bits_eq("mesh/points", &self.mesh_points, &expected.mesh_points)?;

        if self.physics != expected.physics
            || self.study != expected.study
            || self.dim != expected.dim
            || self.dofs != expected.dofs
            || self.n_solutions != expected.n_solutions
            || self.order != expected.order
            || self.seconds.to_bits() != expected.seconds.to_bits()
        {
            bail!("results/{} attributes mismatch", self.physics);
        }
        if self.log_bytes != expected.log_bytes || self.scalars_bytes != expected.scalars_bytes {
            bail!("results/{} log or scalars bytes mismatch", self.physics);
        }
        if self.viz_topology != expected.viz_topology
            || self.viz_cells_shape != expected.viz_cells_shape
            || self.viz_cells != expected.viz_cells
            || self.viz_points_shape != expected.viz_points_shape
        {
            bail!("viz/{} topology or cells mismatch", self.physics);
        }
        assert_f64_bits_eq("viz/points", &self.viz_points, &expected.viz_points)?;

        if self.solutions.len() != expected.solutions.len() {
            bail!("solutions count mismatch: {} vs {}", self.solutions.len(), expected.solutions.len());
        }
        for (k, (a_sol, e_sol)) in self.solutions.iter().zip(expected.solutions.iter()).enumerate() {
            if a_sol.eig_re.to_bits() != e_sol.eig_re.to_bits()
                || a_sol.eig_im.to_bits() != e_sol.eig_im.to_bits()
                || a_sol.label != e_sol.label
            {
                bail!("solution_{} attrs mismatch", k);
            }
            if a_sol.cell_fields.len() != e_sol.cell_fields.len() {
                bail!("solution_{} cell_fields count mismatch", k);
            }
            for (cf_key, (a_name, a_unit, a_nloc, a_re, a_im)) in &a_sol.cell_fields {
                let (e_name, e_unit, e_nloc, e_re, e_im) =
                    e_sol.cell_fields.get(cf_key).ok_or_else(|| anyhow!("missing cell_field {}", cf_key))?;
                if a_name != e_name || a_unit != e_unit || a_nloc != e_nloc {
                    bail!("solution_{}/cell_fields/{} attrs mismatch", k, cf_key);
                }
                assert_f64_bits_eq(&format!("solution_{}/cell_fields/{}/re", k, cf_key), a_re, e_re)?;
                assert_f64_bits_eq(&format!("solution_{}/cell_fields/{}/im", k, cf_key), a_im, e_im)?;
            }
            if a_sol.fields_order_json != e_sol.fields_order_json || a_sol.fields.len() != e_sol.fields.len() {
                bail!("solution_{}/fields metadata mismatch", k);
            }
            for (f_key, (a_name, a_unit, a_comp, a_data)) in &a_sol.fields {
                let (e_name, e_unit, e_comp, e_data) =
                    e_sol.fields.get(f_key).ok_or_else(|| anyhow!("missing field {}", f_key))?;
                if a_name != e_name || a_unit != e_unit || a_comp != e_comp {
                    bail!("solution_{}/fields/{} attrs mismatch", k, f_key);
                }
                assert_f64_bits_eq(&format!("solution_{}/fields/{}", k, f_key), a_data, e_data)?;
            }
            if check_viz_links {
                if a_sol.viz_fields.len() != e_sol.viz_fields.len() {
                    bail!("viz/sol_{} field count mismatch", k);
                }
                for (v_key, a_data) in &a_sol.viz_fields {
                    let e_data = e_sol.viz_fields.get(v_key).ok_or_else(|| anyhow!("missing viz field {}", v_key))?;
                    assert_f64_bits_eq(&format!("viz/sol_{}/{}", k, v_key), a_data, e_data)?;
                }
            }
        }
        Ok(())
    }
}

/// Bounded-memory streaming HDF5 reader for `nyx-results` files (Layouts L1, L2, and L3).
pub struct NyxReader {
    file: File,
    budget: ReadBudgetConfig,
    file_bytes: u64,
}

impl NyxReader {
    /// Opens an `.h5` file in streaming (`pread`) mode without buffering datasets into RAM.
    pub fn open(path: &Path, budget: ReadBudgetConfig) -> Result<Self> {
        let file_bytes = std::fs::metadata(path)?.len();
        let file = File::open_streaming(path)?;
        Ok(Self { file, budget, file_bytes })
    }

    /// Opens an in-memory `.h5` image with [`ReadBudgetConfig`] protection (used by untrusted fuzzing).
    pub fn from_bytes(bytes: Vec<u8>, budget: ReadBudgetConfig) -> Result<Self> {
        let file_bytes = bytes.len() as u64;
        let file = File::from_bytes(bytes)?;
        Ok(Self { file, budget, file_bytes })
    }

    /// Configured read budget for dataset and group traversal operations.
    pub fn budget(&self) -> ReadBudgetConfig {
        self.budget
    }

    /// Underlying `hdf5-pure` [`File`] handle.
    pub fn file(&self) -> &File {
        &self.file
    }

    fn account_bytes(&self, total_bytes: &mut u64, add_bytes: u64, context: &str) -> Result<()> {
        account_budget_bytes(total_bytes, add_bytes, self.budget.max_total_result_bytes, context)
    }

    /// Enforces [`ReadBudgetConfig::max_dataset_bytes`] using checked multiplication across `ds.shape()`.
    pub fn check_dataset_budget(&self, ds: &Dataset, elem_bytes: u64) -> Result<u64> {
        let shape = ds.shape()?;
        check_shape_budget(&shape, elem_bytes, self.budget.max_dataset_bytes)
    }

    fn check_dataset_budget_tracked(
        &self,
        ds: &Dataset,
        elem_bytes: u64,
        total_bytes: &mut u64,
        context: &str,
    ) -> Result<u64> {
        let elems = self.check_dataset_budget(ds, elem_bytes)?;
        let bytes = elems
            .checked_mul(elem_bytes)
            .ok_or_else(|| anyhow!("dataset byte size overflow at {context}"))?;
        self.account_bytes(total_bytes, bytes, context)?;
        Ok(elems)
    }

    fn validate_n_solutions(&self, raw_n: i64, min_bytes_per_solution: u64, context: &str) -> Result<usize> {
        validate_solution_count(raw_n, min_bytes_per_solution, self.budget.max_solutions, self.file_bytes, context)
    }

    fn capped_solution_capacity(&self, n_solutions: usize) -> usize {
        capped_solution_capacity_for(n_solutions, self.budget.max_solutions, self.file_bytes)
    }

    /// Reads raw dataset bytes only after checking datatype size and dataspace dimensions against [`ReadBudgetConfig`].
    pub fn read_raw_checked(&self, ds: &Dataset) -> Result<Vec<u8>> {
        let elem_bytes = u64::from(ds.datatype()?.type_size()).max(1);
        self.check_dataset_budget(ds, elem_bytes)?;
        Ok(ds.read_raw()?)
    }

    /// Fast O(1) metadata scan of eigenmode headers for `physics` (`"ems"`, `"sms"`, etc.).
    /// Supports both packed `eigenvalues` (`[M, 2]` or streaming `[2 * M]`) and legacy per-solution
    /// attributes (`eig_re` / `eig_im` with mechanical `value` fallback and `eig_im = 0.0`).
    pub fn scan_mode_headers(&self, physics: &str) -> Result<Vec<ModeHeader>> {
        let mut total_bytes = 0u64;
        self.scan_mode_headers_with_budget(physics, &mut total_bytes)
    }

    /// Same as [`Self::scan_mode_headers`] while charging cumulative bytes against `total_bytes`.
    pub fn scan_mode_headers_with_budget(&self, physics: &str, total_bytes: &mut u64) -> Result<Vec<ModeHeader>> {
        let phys_path = format!("results/{physics}");
        let phys_grp = self.file.group(&phys_path)?;
        let phys_attrs = phys_grp.attrs()?;
        let raw_n = attr_i64(&phys_attrs, "n_solutions")?;

        let packed_eig_path = format!("{phys_path}/eigenvalues");
        if let Ok(ds_eig) = self.file.dataset(&packed_eig_path) {
            let n_solutions = self.validate_n_solutions(raw_n, 0, &phys_path)?;
            self.check_dataset_budget_tracked(&ds_eig, 8, total_bytes, &packed_eig_path)?;
            let raw = ds_eig.read_f64()?;
            let need = n_solutions
                .checked_mul(2)
                .ok_or_else(|| anyhow!("n_solutions * 2 overflow at {}", packed_eig_path))?;
            if raw.len() < need {
                bail!(
                    "{} length {} < n_solutions * 2 ({})",
                    packed_eig_path,
                    raw.len(),
                    need
                );
            }
            let labels_opt: Option<Vec<String>> =
                attr_str(&phys_attrs, "mode_labels_json").ok().and_then(|j| serde_json::from_str(&j).ok());
            let mut headers = Vec::with_capacity(self.capped_solution_capacity(n_solutions));
            for k in 0..n_solutions {
                let label = labels_opt
                    .as_ref()
                    .and_then(|v| v.get(k).cloned())
                    .or_else(|| {
                        self.file
                            .group(&format!("{phys_path}/solution_{k}"))
                            .ok()
                            .and_then(|g| g.attrs().ok())
                            .and_then(|a| attr_str(&a, "label").ok())
                    })
                    .unwrap_or_else(|| format!("mode_{k}"));
                let hdr_bytes = (std::mem::size_of::<ModeHeader>() as u64).saturating_add(label.len() as u64);
                self.account_bytes(total_bytes, hdr_bytes, &packed_eig_path)?;
                headers.push(ModeHeader {
                    mode_index: k,
                    eig_re: raw[k * 2],
                    eig_im: raw[k * 2 + 1],
                    label,
                });
            }
            return Ok(headers);
        }

        let n_solutions = self.validate_n_solutions(raw_n, 8, &phys_path)?;
        let mut headers = Vec::with_capacity(self.capped_solution_capacity(n_solutions));
        for k in 0..n_solutions {
            let sol_path = format!("{phys_path}/solution_{k}");
            let sol_grp = self.file.group(&sol_path)?;
            let sattrs = sol_grp.attrs()?;
            let eig_re = attr_f64(&sattrs, "eig_re").or_else(|_| attr_f64(&sattrs, "value"))?;
            let eig_im = attr_f64(&sattrs, "eig_im").unwrap_or(0.0);
            let label = attr_str(&sattrs, "label").unwrap_or_default();
            let hdr_bytes = (std::mem::size_of::<ModeHeader>() as u64).saturating_add(label.len() as u64);
            self.account_bytes(total_bytes, hdr_bytes, &sol_path)?;
            headers.push(ModeHeader { mode_index: k, eig_re, eig_im, label });
        }
        Ok(headers)
    }

    /// Reads `/mesh` topology, connectivity, coordinates, and optional curved-edge geometry (`edge_node_pairs`, `edge_node_points`).
    /// Transparently handles both explicit `/mesh` attributes and legacy `nyx-results` v3/v4 files without `/mesh` attributes.
    pub fn read_mesh(&self) -> Result<MeshData> {
        let mut total_bytes = 0u64;
        self.read_mesh_with_budget(&mut total_bytes)
    }

    fn read_mesh_with_budget(&self, total_bytes: &mut u64) -> Result<MeshData> {
        let mesh_grp = self.file.group("mesh")?;
        let mattrs = mesh_grp.attrs()?;

        let ds_cl = self.file.dataset("mesh/cell_labels")?;
        self.check_dataset_budget_tracked(&ds_cl, 4, total_bytes, "mesh/cell_labels")?;
        let cell_labels = ds_cl.read_i32()?;

        let ds_c = self.file.dataset("mesh/cells")?;
        self.check_dataset_budget_tracked(&ds_c, 4, total_bytes, "mesh/cells")?;
        let cs = ds_c.shape()?;
        let cells_shape = shape_2d(&cs, "mesh/cells")?;
        let cells = ds_c.read_u32()?;

        let ds_fl = self.file.dataset("mesh/facet_labels")?;
        self.check_dataset_budget_tracked(&ds_fl, 4, total_bytes, "mesh/facet_labels")?;
        let facet_labels = ds_fl.read_i32()?;

        let ds_f = self.file.dataset("mesh/facets")?;
        self.check_dataset_budget_tracked(&ds_f, 4, total_bytes, "mesh/facets")?;
        let fs = ds_f.shape()?;
        let facets_shape = shape_2d(&fs, "mesh/facets")?;
        let facets = ds_f.read_u32()?;

        let ds_p = self.file.dataset("mesh/points")?;
        self.check_dataset_budget_tracked(&ds_p, 8, total_bytes, "mesh/points")?;
        let ps = ds_p.shape()?;
        let points_shape = shape_2d(&ps, "mesh/points")?;
        let points = ds_p.read_f64()?;

        let mut edge_node_pairs = Vec::new();
        if let Ok(ds_ep) = self.file.dataset("mesh/edge_node_pairs") {
            self.check_dataset_budget_tracked(&ds_ep, 4, total_bytes, "mesh/edge_node_pairs")?;
            let raw_pairs = ds_ep.read_u32()?;
            let (chunks, rem) = raw_pairs.as_chunks::<2>();
            if !rem.is_empty() {
                bail!("mesh/edge_node_pairs length {} is not a multiple of 2", raw_pairs.len());
            }
            edge_node_pairs.extend_from_slice(chunks);
        }

        let mut edge_node_points = Vec::new();
        if let Ok(ds_ept) = self.file.dataset("mesh/edge_node_points") {
            self.check_dataset_budget_tracked(&ds_ept, 8, total_bytes, "mesh/edge_node_points")?;
            let raw_pts = ds_ept.read_f64()?;
            let (chunks, rem) = raw_pts.as_chunks::<3>();
            if !rem.is_empty() {
                bail!("mesh/edge_node_points length {} is not a multiple of 3", raw_pts.len());
            }
            edge_node_points.extend_from_slice(chunks);
        }

        if edge_node_pairs.len() != edge_node_points.len() {
            bail!(
                "mesh curved edge count mismatch: {} pairs vs {} points",
                edge_node_pairs.len(),
                edge_node_points.len()
            );
        }

        let space_dimension = attr_i64(&mattrs, "space_dimension").unwrap_or(points_shape[1] as i64);
        let element_family = attr_str(&mattrs, "element_family").unwrap_or_else(|_| "Lagrange".to_string());
        let element_degree = attr_i64(&mattrs, "element_degree").unwrap_or_else(|_| {
            if (space_dimension == 2 && cells_shape[1] == 6) || (space_dimension == 3 && cells_shape[1] == 10) {
                2
            } else {
                1
            }
        });
        let cell_topology = attr_str(&mattrs, "cell_topology").unwrap_or_else(|_| {
            match (space_dimension, cells_shape[1], element_degree) {
                (2, 6, _) | (2, _, 2) => "Triangle_6".to_string(),
                (2, 3, _) => "Triangle".to_string(),
                (3, 10, _) | (3, _, 2) => "Tetrahedron_10".to_string(),
                (3, 4, _) => "Tetrahedron".to_string(),
                (d, nv, _) => format!("Simplex_{d}_{nv}"),
            }
        });

        Ok(MeshData {
            cell_topology,
            element_family,
            element_degree,
            space_dimension,
            cell_labels,
            cells_shape,
            cells,
            facet_labels,
            facets_shape,
            facets,
            points_shape,
            points,
            edge_node_pairs,
            edge_node_points,
        })
    }

    /// Reads `/results/<physics>/cell_fields_meta` component descriptors if present.
    pub fn read_packed_components_meta(&self, physics: &str) -> Result<Option<Vec<PackedComponentMeta>>> {
        let meta_path = format!("results/{physics}/cell_fields_meta");
        let Ok(ds) = self.file.dataset(&meta_path) else {
            return Ok(None);
        };
        self.check_dataset_budget(&ds, 1)?;
        let bytes = ds.read_u8()?;
        let parsed: Vec<PackedComponentMeta> = serde_json::from_slice(&bytes)?;
        Ok(Some(parsed))
    }

    /// Returns the logical 5D shape `[M, C, N_c, N_v, 2]` of `/results/<physics>/cell_fields_packed`
    /// (for both whole-file 5D datasets and rank-1 streaming datasets), or `Ok(None)` for legacy L1 v3 files.
    pub fn packed_logical_shape(&self, physics: &str) -> Result<Option<[usize; 5]>> {
        let packed_path = format!("results/{physics}/cell_fields_packed");
        let Ok(ds) = self.file.dataset(&packed_path) else {
            return Ok(None);
        };
        let shape = ds.shape()?;
        if shape.len() == 5 {
            return Ok(Some(shape_5d(&shape)?));
        }
        if shape.len() == 1 {
            let attrs = ds.attrs()?;
            let tail = attrs
                .get("logical_shape_tail")
                .and_then(AttrValue::to_i64s)
                .ok_or_else(|| anyhow!("rank-1 cell_fields_packed missing logical_shape_tail attribute"))?;
            let (partial, needs_nsol) = packed_shape_from_tail(shape[0], &tail)?;
            if needs_nsol {
                let phys_path = format!("results/{physics}");
                let pattrs = self.file.group(&phys_path)?.attrs()?;
                let m = self.validate_n_solutions(attr_i64(&pattrs, "n_solutions")?, 0, &phys_path)?;
                return Ok(Some([m, partial[1], partial[2], partial[3], partial[4]]));
            }
            return Ok(Some(partial));
        }
        bail!("unexpected rank {} for {}", shape.len(), packed_path)
    }

    /// Bounded-memory single-mode slice reader (Gate 4 `EC-09` & `BM-2`).
    /// Dispatches to `ds.read_f64_rows` on `L2`/`L3`/`L1-packed` datasets (decompressing only chunk `mode_idx`),
    /// or reconstructs `[C, N_c, N_v, 2]` (`[re, im]` interleaved) from legacy `/results/<phys>/solution_<k>/cell_fields`.
    pub fn read_mode_cell_fields(&self, physics: &str, mode_idx: usize) -> Result<Vec<f64>> {
        let mut total_bytes = 0u64;
        self.read_mode_cell_fields_with_budget(physics, mode_idx, &mut total_bytes)
    }

    /// Same as [`Self::read_mode_cell_fields`] while charging cumulative bytes against `total_bytes`.
    pub fn read_mode_cell_fields_with_budget(
        &self,
        physics: &str,
        mode_idx: usize,
        total_bytes: &mut u64,
    ) -> Result<Vec<f64>> {
        let packed_path = format!("results/{physics}/cell_fields_packed");
        if let Ok(ds) = self.file.dataset(&packed_path) {
            let shape = ds.shape()?;
            if shape.len() == 5 {
                let (row_elems, row_bytes) =
                    packed_5d_slice_bytes(&shape, mode_idx, self.budget.max_dataset_bytes)?;
                self.account_bytes(total_bytes, row_bytes, &packed_path)?;
                if row_elems == 0 {
                    return Ok(Vec::new());
                }
                return Ok(ds.read_f64_rows(mode_idx as u64, 1)?);
            }
            if shape.len() == 1 {
                let attrs = ds.attrs()?;
                let tail = attrs
                    .get("logical_shape_tail")
                    .and_then(AttrValue::to_i64s)
                    .ok_or_else(|| anyhow!("rank-1 cell_fields_packed missing logical_shape_tail attribute"))?;
                let (mode_elems, mode_bytes, start) =
                    packed_1d_slice_range(shape[0], &tail, mode_idx, self.budget.max_dataset_bytes)?;
                self.account_bytes(total_bytes, mode_bytes, &packed_path)?;
                if mode_elems == 0 {
                    let phys_path = format!("results/{physics}");
                    let pattrs = self.file.group(&phys_path)?.attrs()?;
                    let n_sol = self.validate_n_solutions(attr_i64(&pattrs, "n_solutions")?, 0, &phys_path)?;
                    if mode_idx >= n_sol {
                        bail!("mode_idx {} out of bounds for n_solutions={}", mode_idx, n_sol);
                    }
                    return Ok(Vec::new());
                }
                return Ok(ds.read_f64_rows(start, mode_elems)?);
            }
            bail!("unsupported cell_fields_packed rank {} (shape {:?})", shape.len(), shape);
        }

        // Legacy L1/v3/v4 fallback: read `/results/<phys>/solution_<k>/cell_fields/<comp>/{re,im}`
        let cf_path = format!("results/{physics}/solution_{mode_idx}/cell_fields");
        let cf_grp = self.file.group(&cf_path)?;
        let mut comp_names = cf_grp.groups()?;
        comp_names.sort();
        let mut out = Vec::new();
        for cname in comp_names {
            let re_path = format!("{cf_path}/{cname}/re");
            let ds_re = self.file.dataset(&re_path)?;
            self.check_dataset_budget_tracked(&ds_re, 8, total_bytes, &re_path)?;
            let re = ds_re.read_f64()?;

            let im_path = format!("{cf_path}/{cname}/im");
            let im = if let Ok(ds_im) = self.file.dataset(&im_path) {
                self.check_dataset_budget_tracked(&ds_im, 8, total_bytes, &im_path)?;
                ds_im.read_f64()?
            } else {
                let im_bytes = (re.len() as u64)
                    .checked_mul(8)
                    .ok_or_else(|| anyhow!("{}: im byte size overflow", im_path))?;
                self.account_bytes(total_bytes, im_bytes, &im_path)?;
                vec![0.0; re.len()]
            };
            if re.len() != im.len() {
                bail!("{}/{}: re.len() {} != im.len() {}", cf_path, cname, re.len(), im.len());
            }
            let inter_elems = re
                .len()
                .checked_mul(2)
                .ok_or_else(|| anyhow!("{}/{}: interleaved length overflow", cf_path, cname))?;
            out.try_reserve(inter_elems)
                .map_err(|e| anyhow!("{}/{}: allocation failed for {} elements: {}", cf_path, cname, inter_elems, e))?;
            for (&r, &i) in re.iter().zip(im.iter()) {
                out.push(r);
                out.push(i);
            }
        }
        Ok(out)
    }

    /// Reads the complete `nyx-results` v3 structure for the first alphabetically sorted physics group under `/results`.
    pub fn read_nyx_v3_summary(&self, read_viz_links: bool) -> Result<NyxV3Summary> {
        let res_root = self.file.group("results")?;
        let mut phys_groups = res_root.groups()?;
        phys_groups.sort();
        let physics = phys_groups.first().ok_or_else(|| anyhow!("no physics group under /results"))?.clone();
        self.read_nyx_v3_summary_for_physics(&physics, read_viz_links)
    }

    /// Reads the complete `nyx-results` v3 structure for a specific `physics` group under `/results`.
    pub fn read_nyx_v3_summary_for_physics(&self, physics: &str, read_viz_links: bool) -> Result<NyxV3Summary> {
        let mut total_bytes = 0u64;
        self.read_nyx_v3_summary_for_physics_with_budget(physics, read_viz_links, &mut total_bytes)
    }

    /// Same as [`Self::read_nyx_v3_summary_for_physics`] while charging cumulative bytes against `total_bytes`.
    pub fn read_nyx_v3_summary_for_physics_with_budget(
        &self,
        physics: &str,
        read_viz_links: bool,
        total_bytes: &mut u64,
    ) -> Result<NyxV3Summary> {
        let root_attrs = self.file.root().attrs()?;
        let format = attr_str(&root_attrs, "format")?;
        let version = attr_i64(&root_attrs, "version")?;
        let writer = attr_str(&root_attrs, "writer")?;
        let processes = attr_i64(&root_attrs, "processes")?;

        let ds_proj = self.file.dataset("project")?;
        self.check_dataset_budget_tracked(&ds_proj, 1, total_bytes, "project")?;
        let project_bytes = ds_proj.read_u8()?;

        let ds_cell_labels = self.file.dataset("mesh/cell_labels")?;
        self.check_dataset_budget_tracked(&ds_cell_labels, 4, total_bytes, "mesh/cell_labels")?;
        let mesh_cell_labels = ds_cell_labels.read_i32()?;

        let ds_cells = self.file.dataset("mesh/cells")?;
        self.check_dataset_budget_tracked(&ds_cells, 4, total_bytes, "mesh/cells")?;
        let mesh_cells_shape = shape_2d(&ds_cells.shape()?, "mesh/cells")?;
        let mesh_cells = ds_cells.read_u32()?;

        let ds_facet_labels = self.file.dataset("mesh/facet_labels")?;
        self.check_dataset_budget_tracked(&ds_facet_labels, 4, total_bytes, "mesh/facet_labels")?;
        let mesh_facet_labels = ds_facet_labels.read_i32()?;

        let ds_facets = self.file.dataset("mesh/facets")?;
        self.check_dataset_budget_tracked(&ds_facets, 4, total_bytes, "mesh/facets")?;
        let mesh_facets_shape = shape_2d(&ds_facets.shape()?, "mesh/facets")?;
        let mesh_facets = ds_facets.read_u32()?;

        let ds_points = self.file.dataset("mesh/points")?;
        self.check_dataset_budget_tracked(&ds_points, 8, total_bytes, "mesh/points")?;
        let mesh_points_shape = shape_2d(&ds_points.shape()?, "mesh/points")?;
        let mesh_points = ds_points.read_f64()?;

        let phys_path = format!("results/{physics}");
        let phys_grp = self.file.group(&phys_path)?;
        let phys_attrs = phys_grp.attrs()?;
        let study = attr_str(&phys_attrs, "study")?;
        let dim = attr_i64(&phys_attrs, "dim")?;
        let dofs = attr_i64(&phys_attrs, "dofs")?;
        let n_solutions = attr_i64(&phys_attrs, "n_solutions")?;
        let order = attr_i64(&phys_attrs, "order")?;
        let seconds = attr_f64(&phys_attrs, "seconds")?;

        let ds_log = phys_grp.dataset("log")?;
        self.check_dataset_budget_tracked(&ds_log, 1, total_bytes, &format!("{phys_path}/log"))?;
        let log_bytes = ds_log.read_u8()?;

        let ds_scalars = phys_grp.dataset("scalars")?;
        self.check_dataset_budget_tracked(&ds_scalars, 1, total_bytes, &format!("{phys_path}/scalars"))?;
        let scalars_bytes = ds_scalars.read_u8()?;

        let viz_phys = self.file.group(&format!("viz/{physics}"))?;
        let viz_attrs = viz_phys.attrs()?;
        let viz_topology = attr_str(&viz_attrs, "topology")?;

        let ds_vcells = viz_phys.dataset("cells")?;
        self.check_dataset_budget_tracked(&ds_vcells, 4, total_bytes, &format!("viz/{physics}/cells"))?;
        let viz_cells_shape = shape_2d(&ds_vcells.shape()?, "viz/cells")?;
        let viz_cells = ds_vcells.read_u32()?;

        let ds_vpoints = viz_phys.dataset("points")?;
        self.check_dataset_budget_tracked(&ds_vpoints, 8, total_bytes, &format!("viz/{physics}/points"))?;
        let viz_points_shape = shape_2d(&ds_vpoints.shape()?, "viz/points")?;
        let viz_points = ds_vpoints.read_f64()?;

        let n_sol = self.validate_n_solutions(n_solutions, 8, &phys_path)?;
        let mut solutions = Vec::with_capacity(self.capped_solution_capacity(n_sol));
        for k in 0..n_sol {
            let sol_grp = phys_grp.group(&format!("solution_{k}"))?;
            let sol_attrs = sol_grp.attrs()?;
            let eig_re = attr_f64(&sol_attrs, "eig_re").or_else(|_| attr_f64(&sol_attrs, "value")).unwrap_or(0.0);
            let eig_im = attr_f64(&sol_attrs, "eig_im").unwrap_or(0.0);
            let label = attr_str(&sol_attrs, "label")?;

            let cf_grp = sol_grp.group("cell_fields")?;
            let mut cf_names = cf_grp.groups()?;
            cf_names.sort();
            let mut cell_fields = BTreeMap::new();
            for cf_key in cf_names {
                let comp_grp = cf_grp.group(&cf_key)?;
                let comp_attrs = comp_grp.attrs()?;
                let c_name = attr_str(&comp_attrs, "name")?;
                let c_unit = attr_str(&comp_attrs, "unit")?;
                let c_nloc = attr_i64(&comp_attrs, "nloc")?;

                let re_ctx = format!("{phys_path}/solution_{k}/cell_fields/{cf_key}/re");
                let ds_re = comp_grp.dataset("re")?;
                self.check_dataset_budget_tracked(&ds_re, 8, total_bytes, &re_ctx)?;
                let re = ds_re.read_f64()?;

                let im_ctx = format!("{phys_path}/solution_{k}/cell_fields/{cf_key}/im");
                let im = if let Ok(ds_im) = comp_grp.dataset("im") {
                    self.check_dataset_budget_tracked(&ds_im, 8, total_bytes, &im_ctx)?;
                    ds_im.read_f64()?
                } else {
                    let im_bytes = (re.len() as u64)
                        .checked_mul(8)
                        .ok_or_else(|| anyhow!("{}: im byte size overflow", im_ctx))?;
                    self.account_bytes(total_bytes, im_bytes, &im_ctx)?;
                    vec![0.0; re.len()]
                };
                if re.len() != im.len() {
                    bail!("{phys_path}/solution_{k}/cell_fields/{cf_key}: re.len() {} != im.len() {}", re.len(), im.len());
                }
                cell_fields.insert(cf_key, (c_name, c_unit, c_nloc, re, im));
            }

            let f_grp = sol_grp.group("fields")?;
            let f_attrs = f_grp.attrs()?;
            let fields_order_json = attr_str(&f_attrs, "order")?;
            let mut f_names = f_grp.datasets()?;
            f_names.sort();
            let mut fields = BTreeMap::new();
            for f_key in f_names {
                let ds = f_grp.dataset(&f_key)?;
                self.check_dataset_budget_tracked(
                    &ds,
                    8,
                    total_bytes,
                    &format!("{phys_path}/solution_{k}/fields/{f_key}"),
                )?;
                let ds_attrs = ds.attrs()?;
                let f_name = attr_str(&ds_attrs, "name")?;
                let f_unit = attr_str(&ds_attrs, "unit")?;
                let f_comp = attr_str(&ds_attrs, "components")?;
                let data = ds.read_f64()?;
                fields.insert(f_key, (f_name, f_unit, f_comp, data));
            }

            let mut viz_fields = BTreeMap::new();
            if read_viz_links {
                let vsol_grp = viz_phys.group(&format!("sol_{k}"))?;
                let mut v_names = vsol_grp.datasets()?;
                v_names.sort();
                for v_key in v_names {
                    let v_ds = vsol_grp.dataset(&v_key)?;
                    self.check_dataset_budget_tracked(
                        &v_ds,
                        8,
                        total_bytes,
                        &format!("viz/{physics}/sol_{k}/{v_key}"),
                    )?;
                    let v_data = v_ds.read_f64()?;
                    viz_fields.insert(v_key, v_data);
                }
            }

            solutions.push(NyxSolutionData {
                eig_re,
                eig_im,
                label,
                cell_fields,
                fields_order_json,
                fields,
                viz_fields,
            });
        }

        Ok(NyxV3Summary {
            format,
            version,
            writer,
            processes,
            project_bytes,
            mesh_cell_labels,
            mesh_cells_shape,
            mesh_cells,
            mesh_facet_labels,
            mesh_facets_shape,
            mesh_facets,
            mesh_points_shape,
            mesh_points,
            physics: physics.to_string(),
            study,
            dim,
            dofs,
            n_solutions,
            order,
            seconds,
            log_bytes,
            scalars_bytes,
            viz_topology,
            viz_cells_shape,
            viz_cells,
            viz_points_shape,
            viz_points,
            solutions,
        })
    }
}

fn check_shape_budget(shape: &[u64], elem_bytes: u64, max_dataset_bytes: u64) -> Result<u64> {
    let elems = shape
        .iter()
        .try_fold(1u64, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| anyhow!("dataset shape overflow: {:?}", shape))?;
    let bytes = elems.checked_mul(elem_bytes).ok_or_else(|| anyhow!("dataset byte size overflow"))?;
    if bytes > max_dataset_bytes {
        bail!("dataset size {} B exceeds configured read budget {} B", bytes, max_dataset_bytes);
    }
    Ok(elems)
}

fn account_budget_bytes(
    total_bytes: &mut u64,
    add_bytes: u64,
    max_total_result_bytes: u64,
    context: &str,
) -> Result<()> {
    let next = total_bytes
        .checked_add(add_bytes)
        .ok_or_else(|| anyhow!("cumulative result byte count overflow at {context}"))?;
    if next > max_total_result_bytes {
        bail!(
            "cumulative result size {} B exceeds configured total read budget {} B at {context}",
            next,
            max_total_result_bytes
        );
    }
    *total_bytes = next;
    Ok(())
}

fn validate_solution_count(
    raw_n: i64,
    min_bytes_per_solution: u64,
    max_solutions: usize,
    file_bytes: u64,
    context: &str,
) -> Result<usize> {
    if raw_n < 0 {
        bail!("negative n_solutions {} at {}", raw_n, context);
    }
    let n_solutions = usize::try_from(raw_n)
        .map_err(|_| anyhow!("n_solutions {} overflows usize at {}", raw_n, context))?;
    if n_solutions > max_solutions {
        bail!(
            "n_solutions {} exceeds configured max_solutions budget {} at {}",
            n_solutions,
            max_solutions,
            context
        );
    }
    if min_bytes_per_solution > 0 {
        let min_needed = (n_solutions as u64).saturating_mul(min_bytes_per_solution);
        if min_needed > file_bytes {
            bail!(
                "declared n_solutions {} requires at least {} B, exceeding file size {} B at {}",
                n_solutions,
                min_needed,
                file_bytes,
                context
            );
        }
    }
    Ok(n_solutions)
}

fn capped_solution_capacity_for(n_solutions: usize, max_solutions: usize, file_bytes: u64) -> usize {
    let file_cap = (file_bytes / 16).min(usize::MAX as u64) as usize;
    n_solutions.min(max_solutions).min(4096).min(file_cap)
}

fn shape_2d(shape: &[u64], label: &str) -> Result<[usize; 2]> {
    match shape {
        [r, c] => Ok([
            usize::try_from(*r).map_err(|_| anyhow!("2D shape row {} overflows usize for {}", r, label))?,
            usize::try_from(*c).map_err(|_| anyhow!("2D shape col {} overflows usize for {}", c, label))?,
        ]),
        [0] => Ok([0, 0]),
        _ => bail!("expected 2D shape for {}, got {:?}", label, shape),
    }
}

fn shape_5d(shape: &[u64]) -> Result<[usize; 5]> {
    if shape.len() != 5 {
        bail!("expected 5D shape, got {:?}", shape);
    }
    Ok([
        usize::try_from(shape[0]).map_err(|_| anyhow!("5D dim 0 overflows usize"))?,
        usize::try_from(shape[1]).map_err(|_| anyhow!("5D dim 1 overflows usize"))?,
        usize::try_from(shape[2]).map_err(|_| anyhow!("5D dim 2 overflows usize"))?,
        usize::try_from(shape[3]).map_err(|_| anyhow!("5D dim 3 overflows usize"))?,
        usize::try_from(shape[4]).map_err(|_| anyhow!("5D dim 4 overflows usize"))?,
    ])
}

fn packed_shape_from_tail(flat_len: u64, tail: &[i64]) -> Result<([usize; 5], bool)> {
    if tail.len() != 4 {
        bail!("invalid logical_shape_tail length {}: expected 4", tail.len());
    }
    let c = usize::try_from(tail[0].max(0)).map_err(|_| anyhow!("logical_shape_tail[0] overflows usize"))?;
    let n_c = usize::try_from(tail[1].max(0)).map_err(|_| anyhow!("logical_shape_tail[1] overflows usize"))?;
    let n_v = usize::try_from(tail[2].max(0)).map_err(|_| anyhow!("logical_shape_tail[2] overflows usize"))?;
    let ri = usize::try_from(tail[3].max(0)).map_err(|_| anyhow!("logical_shape_tail[3] overflows usize"))?;
    let mode_elems = c
        .checked_mul(n_c)
        .and_then(|x| x.checked_mul(n_v))
        .and_then(|x| x.checked_mul(ri))
        .ok_or_else(|| anyhow!("logical_shape_tail overflow"))?;
    let flat_usize = usize::try_from(flat_len).map_err(|_| anyhow!("packed flat_len overflows usize"))?;
    if let Some(m) = flat_usize.checked_div(mode_elems) {
        Ok(([m, c, n_c, n_v, ri], false))
    } else {
        Ok(([0, c, n_c, n_v, ri], true))
    }
}

fn packed_5d_slice_bytes(shape: &[u64], mode_idx: usize, max_dataset_bytes: u64) -> Result<(u64, u64)> {
    if shape.len() != 5 {
        bail!("expected 5D shape, got {:?}", shape);
    }
    let mode_u64 = u64::try_from(mode_idx).map_err(|_| anyhow!("mode_idx overflows u64"))?;
    if mode_u64 >= shape[0] {
        bail!("mode_idx {} out of bounds for shape {:?}", mode_idx, shape);
    }
    let row_elems = shape[1..]
        .iter()
        .try_fold(1u64, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| anyhow!("row shape overflow: {:?}", shape))?;
    let row_bytes = row_elems.checked_mul(8).ok_or_else(|| anyhow!("row byte size overflow"))?;
    if row_bytes > max_dataset_bytes {
        bail!("single mode slice ({} B) exceeds read budget {} B", row_bytes, max_dataset_bytes);
    }
    Ok((row_elems, row_bytes))
}

fn packed_1d_slice_range(
    flat_len: u64,
    tail: &[i64],
    mode_idx: usize,
    max_dataset_bytes: u64,
) -> Result<(u64, u64, u64)> {
    if tail.len() != 4 {
        bail!("invalid logical_shape_tail length {}: expected 4", tail.len());
    }
    let mode_elems = tail
        .iter()
        .try_fold(1u64, |acc, &d| acc.checked_mul(d.max(0) as u64))
        .ok_or_else(|| anyhow!("mode_elems overflow"))?;
    let mode_bytes = mode_elems.checked_mul(8).ok_or_else(|| anyhow!("mode_bytes overflow"))?;
    if mode_bytes > max_dataset_bytes {
        bail!("single mode slice ({} B) exceeds read budget {} B", mode_bytes, max_dataset_bytes);
    }
    if mode_elems == 0 {
        return Ok((0, mode_bytes, 0));
    }
    let mode_u64 = u64::try_from(mode_idx).map_err(|_| anyhow!("mode_idx overflows u64"))?;
    let start = mode_u64.checked_mul(mode_elems).ok_or_else(|| anyhow!("mode slice start offset overflow"))?;
    let end = start.checked_add(mode_elems).ok_or_else(|| anyhow!("mode slice end offset overflow"))?;
    if end > flat_len {
        bail!("mode_idx {} (elems {}..{}) out of bounds for len {}", mode_idx, start, end, flat_len);
    }
    Ok((mode_elems, mode_bytes, start))
}

fn walk_group_checked(reader: &NyxReader, grp: &Group, depth: usize, total_bytes: &mut u64) -> Result<usize> {
    if depth > reader.budget.max_group_depth {
        bail!("group recursion depth {} exceeds configured limit {}", depth, reader.budget.max_group_depth);
    }
    let mut local_bytes = 0usize;
    let _ = grp.attrs()?;
    for ds_name in grp.datasets()? {
        let ds = grp.dataset(&ds_name)?;
        let _ = ds.attrs()?;
        let elem_bytes = u64::from(ds.datatype()?.type_size()).max(1);
        reader.check_dataset_budget_tracked(&ds, elem_bytes, total_bytes, &ds_name)?;
        let raw = ds.read_raw()?;
        local_bytes = local_bytes.saturating_add(raw.len());
    }
    for sub_name in grp.groups()? {
        let sub = grp.group(&sub_name)?;
        local_bytes = local_bytes.saturating_add(walk_group_checked(reader, &sub, depth + 1, total_bytes)?);
    }
    Ok(local_bytes)
}

/// Walks an untrusted HDF5 byte buffer enforcing [`ReadBudgetConfig`] (`checked_mul` byte ceiling
/// and recursion depth limit) without panicking or OOM-aborting (`Gate 5 / EC-11`).
pub fn fuzz_walk_untrusted(bytes: &[u8], budget: ReadBudgetConfig) -> Result<usize> {
    let reader = NyxReader::from_bytes(bytes.to_vec(), budget)?;
    let root = reader.file().root();
    let mut total_bytes = 0u64;
    walk_group_checked(&reader, &root, 0, &mut total_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use hdf5_pure::FileBuilder;

    fn build_packed_bytes(n_sol_attr: i64, eig_pairs: usize) -> Result<Vec<u8>> {
        let mut fb = FileBuilder::new();
        let mut res = fb.create_group("results");
        let mut ems = res.create_group("ems");
        ems.set_attr("n_solutions", AttrValue::I64(n_sol_attr));
        let eigs: Vec<f64> = (0..eig_pairs * 2).map(|i| i as f64).collect();
        ems.create_dataset("eigenvalues").with_f64_data(&eigs).with_shape(&[eig_pairs as u64, 2]);
        res.add_group(ems.finish());
        fb.add_group(res.finish());
        Ok(fb.finish()?)
    }

    fn build_legacy_v3_bytes(
        n_sol_attr: i64,
        actual_sols: usize,
        n_comp: usize,
        n_vals: usize,
        with_im: bool,
    ) -> Result<Vec<u8>> {
        let mut fb = FileBuilder::new();
        fb.set_attr("format", AttrValue::VarLenString("nyx-results".into()));
        fb.set_attr("version", AttrValue::I64(3));
        fb.set_attr("writer", AttrValue::VarLenString("nyx".into()));
        fb.set_attr("processes", AttrValue::I64(1));
        fb.create_dataset("project").with_u8_data(b"{}");

        let mut mg = fb.create_group("mesh");
        mg.create_dataset("cell_labels").with_i32_data(&[1]);
        mg.create_dataset("cells").with_u32_data(&[0, 1, 2]).with_shape(&[1, 3]);
        mg.create_dataset("facet_labels").with_i32_data(&[1]);
        mg.create_dataset("facets").with_u32_data(&[0, 1]).with_shape(&[1, 2]);
        mg.create_dataset("points").with_f64_data(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]).with_shape(&[3, 2]);
        fb.add_group(mg.finish());

        let mut res = fb.create_group("results");
        let mut ems = res.create_group("ems");
        ems.set_attr("physics", AttrValue::VarLenString("ems".into()));
        ems.set_attr("study", AttrValue::VarLenString("modes".into()));
        ems.set_attr("dim", AttrValue::I64(2));
        ems.set_attr("dofs", AttrValue::I64(3));
        ems.set_attr("n_solutions", AttrValue::I64(n_sol_attr));
        ems.set_attr("order", AttrValue::I64(1));
        ems.set_attr("seconds", AttrValue::F64(0.1));
        ems.create_dataset("log").with_u8_data(b"ok");
        ems.create_dataset("scalars").with_u8_data(b"{}");

        let mut viz = fb.create_group("viz");
        let mut vems = viz.create_group("ems");
        vems.set_attr("topology", AttrValue::VarLenString("Triangle".into()));
        vems.create_dataset("cells").with_u32_data(&[0, 1, 2]).with_shape(&[1, 3]);
        vems.create_dataset("points").with_f64_data(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]).with_shape(&[3, 2]);

        let vals = vec![1.0f64; n_vals];
        for k in 0..actual_sols {
            let mut sg = ems.create_group(&format!("solution_{k}"));
            sg.set_attr("eig_re", AttrValue::F64(2.5));
            sg.set_attr("eig_im", AttrValue::F64(0.0));
            sg.set_attr("label", AttrValue::VarLenString(format!("mode_{k}")));

            let mut cf = sg.create_group("cell_fields");
            for c in 0..n_comp {
                let mut cg = cf.create_group(&format!("{c:02}_E"));
                cg.set_attr("name", AttrValue::VarLenString(format!("E{c}")));
                cg.set_attr("unit", AttrValue::VarLenString("V/m".into()));
                cg.set_attr("nloc", AttrValue::I64(3));
                cg.create_dataset("re").with_f64_data(&vals);
                if with_im {
                    cg.create_dataset("im").with_f64_data(&vals);
                }
                cf.add_group(cg.finish());
            }
            sg.add_group(cf.finish());

            let mut fg = sg.create_group("fields");
            fg.set_attr("order", AttrValue::VarLenString("[\"E\"]".into()));
            fg.create_dataset("00_E")
                .with_f64_data(&vals)
                .set_attr("name", AttrValue::VarLenString("E".into()))
                .set_attr("unit", AttrValue::VarLenString("V/m".into()))
                .set_attr("components", AttrValue::VarLenString("[\"x\"]".into()));
            sg.add_group(fg.finish());
            ems.add_group(sg.finish());

            let mut vsg = vems.create_group(&format!("sol_{k}"));
            vsg.create_dataset("E").with_f64_data(&vals);
            vems.add_group(vsg.finish());
        }
        res.add_group(ems.finish());
        fb.add_group(res.finish());
        viz.add_group(vems.finish());
        fb.add_group(viz.finish());
        Ok(fb.finish()?)
    }

    #[test]
    fn oversized_n_solutions_fails_cleanly_without_oom() -> Result<()> {
        let bomb_packed = build_packed_bytes(i64::MAX, 1)?;
        let r = NyxReader::from_bytes(bomb_packed, ReadBudgetConfig::default())?;
        let err = r.scan_mode_headers("ems").err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err.contains("max_solutions"), "unexpected packed bomb error: {err}");

        let over_custom = build_packed_bytes(50, 50)?;
        let r_cap = NyxReader::from_bytes(
            over_custom,
            ReadBudgetConfig { max_solutions: 16, ..ReadBudgetConfig::default() },
        )?;
        let err_cap = r_cap.scan_mode_headers("ems").err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err_cap.contains("max_solutions"), "unexpected custom max_solutions error: {err_cap}");

        // Legacy layout with forged n_solutions exceeding file size / max_solutions fails before looping
        let bomb_legacy = build_legacy_v3_bytes(50_000, 1, 1, 3, true)?;
        let r_leg = NyxReader::from_bytes(bomb_legacy, ReadBudgetConfig::default())?;
        let err_leg_hdr = r_leg.scan_mode_headers("ems").err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err_leg_hdr.contains("exceeding file size"), "unexpected legacy header error: {err_leg_hdr}");
        let err_leg_sum = r_leg
            .read_nyx_v3_summary_for_physics("ems", true)
            .err()
            .ok_or_else(|| anyhow!("expected error"))?
            .to_string();
        assert!(err_leg_sum.contains("exceeding file size"), "unexpected legacy summary error: {err_leg_sum}");
        Ok(())
    }

    #[test]
    fn scan_mode_headers_enforces_cumulative_byte_budget() -> Result<()> {
        let packed = build_packed_bytes(10, 10)?;
        // 10 * 2 * 8 = 160 B for eigenvalues alone, plus ModeHeader accounting
        let budget = ReadBudgetConfig {
            max_dataset_bytes: 1024,
            max_total_result_bytes: 120,
            ..ReadBudgetConfig::default()
        };
        let r = NyxReader::from_bytes(packed, budget)?;
        let err = r.scan_mode_headers("ems").err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err.contains("cumulative result size"), "unexpected header budget error: {err}");
        Ok(())
    }

    #[test]
    fn read_mode_cell_fields_and_v3_summary_enforce_cumulative_byte_budget() -> Result<()> {
        // 4 components * (64 re + 64 im) * 8 bytes = 4096 B total; each individual dataset is 512 B.
        let bytes_with_im = build_legacy_v3_bytes(2, 2, 4, 64, true)?;
        let budget = ReadBudgetConfig {
            max_dataset_bytes: 1024,
            max_total_result_bytes: 2048,
            ..ReadBudgetConfig::default()
        };
        let r = NyxReader::from_bytes(bytes_with_im, budget)?;
        let err_cf = r.read_mode_cell_fields("ems", 0).err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err_cf.contains("cumulative result size"), "unexpected cell_fields budget error: {err_cf}");

        let err_sum = r
            .read_nyx_v3_summary_for_physics("ems", true)
            .err()
            .ok_or_else(|| anyhow!("expected error"))?
            .to_string();
        assert!(err_sum.contains("cumulative result size"), "unexpected summary budget error: {err_sum}");

        // Without `im` dataset: synthetic zero-im still counts toward cumulative `re.len() * 2 * 8` budget.
        let bytes_no_im = build_legacy_v3_bytes(1, 1, 4, 64, false)?;
        let r_no_im = NyxReader::from_bytes(bytes_no_im, budget)?;
        let err_no_im = r_no_im.read_mode_cell_fields("ems", 0).err().ok_or_else(|| anyhow!("expected error"))?.to_string();
        assert!(err_no_im.contains("cumulative result size"), "unexpected no-im budget error: {err_no_im}");
        Ok(())
    }
}

#[cfg(kani)]
mod kani_proofs {
    use super::*;

    #[kani::proof]
    #[kani::unwind(16)]
    fn verify_read_budget_no_overflow() {
        let max_ds: u64 = kani::any();
        let max_tot: u64 = kani::any();
        let max_sol: usize = kani::any();
        let cfg = ReadBudgetConfig {
            max_dataset_bytes: max_ds,
            max_total_result_bytes: max_tot,
            max_group_depth: 24,
            max_solutions: max_sol,
        };

        // 1. `check_shape_budget` across arbitrary rank <= 5 and element sizes.
        let dims: [u64; 5] = kani::any();
        let rank: usize = kani::any();
        kani::assume(rank <= 5);
        let elem_bytes: u64 = kani::any();
        if let Ok(elems) = check_shape_budget(&dims[..rank], elem_bytes, cfg.max_dataset_bytes) {
            if let Some(bytes) = elems.checked_mul(elem_bytes) {
                assert!(bytes <= cfg.max_dataset_bytes);
            }
        }

        // 2. Cumulative byte accounting (`account_budget_bytes`).
        let mut total: u64 = kani::any();
        let add: u64 = kani::any();
        if account_budget_bytes(&mut total, add, cfg.max_total_result_bytes, "kani").is_ok() {
            assert!(total <= cfg.max_total_result_bytes);
        }

        // 3. Solution count and capacity capping.
        let raw_n: i64 = kani::any();
        let min_per_sol: u64 = kani::any();
        let file_bytes: u64 = kani::any();
        if let Ok(n_sol) = validate_solution_count(raw_n, min_per_sol, cfg.max_solutions, file_bytes, "kani") {
            assert!(n_sol <= cfg.max_solutions);
            let cap = capped_solution_capacity_for(n_sol, cfg.max_solutions, file_bytes);
            assert!(cap <= cfg.max_solutions && cap <= 4096);
        }

        // 4. `shape_2d` & `packed_logical_shape` tail arithmetic.
        let s2: [u64; 3] = kani::any();
        let s2_len: usize = kani::any();
        kani::assume(s2_len <= 3);
        let _ = shape_2d(&s2[..s2_len], "kani");

        let flat_len: u64 = kani::any();
        let tail: [i64; 4] = kani::any();
        if let Ok(([m, c, n_c, n_v, ri], _needs_nsol)) = packed_shape_from_tail(flat_len, &tail) {
            let per_mode = c.checked_mul(n_c).and_then(|x| x.checked_mul(n_v)).and_then(|x| x.checked_mul(ri));
            assert!(per_mode.is_some());
            if let Some(pm) = per_mode {
                if pm > 0 {
                    assert!(m.checked_mul(pm).is_some());
                }
            }
        }

        // 5. `read_mode_cell_fields` 5D & 1D index arithmetic.
        let mode_idx: usize = kani::any();
        if let Ok((_re, rb)) = packed_5d_slice_bytes(&dims, mode_idx, cfg.max_dataset_bytes) {
            assert!(rb <= cfg.max_dataset_bytes);
            assert!((mode_idx as u64) < dims[0]);
        }
        if let Ok((mode_elems, mode_bytes, start)) =
            packed_1d_slice_range(flat_len, &tail, mode_idx, cfg.max_dataset_bytes)
        {
            assert!(mode_bytes <= cfg.max_dataset_bytes);
            if mode_elems > 0 {
                let end = start.checked_add(mode_elems);
                assert!(end.is_some());
                assert!(end.unwrap() <= flat_len);
            }
        }
    }
}
