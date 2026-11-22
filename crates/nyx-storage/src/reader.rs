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
}

impl NyxReader {
    /// Opens an `.h5` file in streaming (`pread`) mode without buffering datasets into RAM.
    pub fn open(path: &Path, budget: ReadBudgetConfig) -> Result<Self> {
        let file = File::open_streaming(path)?;
        Ok(Self { file, budget })
    }

    /// Opens an in-memory `.h5` image with [`ReadBudgetConfig`] protection (used by untrusted fuzzing).
    pub fn from_bytes(bytes: Vec<u8>, budget: ReadBudgetConfig) -> Result<Self> {
        let file = File::from_bytes(bytes)?;
        Ok(Self { file, budget })
    }

    /// Configured read budget for dataset and group traversal operations.
    pub fn budget(&self) -> ReadBudgetConfig {
        self.budget
    }

    /// Underlying `hdf5-pure` [`File`] handle.
    pub fn file(&self) -> &File {
        &self.file
    }

    /// Enforces [`ReadBudgetConfig::max_dataset_bytes`] using checked multiplication across `ds.shape()`.
    pub fn check_dataset_budget(&self, ds: &Dataset, elem_bytes: u64) -> Result<u64> {
        let shape = ds.shape()?;
        let elems = shape
            .iter()
            .try_fold(1u64, |acc, &d| acc.checked_mul(d))
            .ok_or_else(|| anyhow!("dataset shape overflow: {:?}", shape))?;
        let bytes = elems.checked_mul(elem_bytes).ok_or_else(|| anyhow!("dataset byte size overflow"))?;
        if bytes > self.budget.max_dataset_bytes {
            bail!("dataset size {} B exceeds configured read budget {} B", bytes, self.budget.max_dataset_bytes);
        }
        Ok(elems)
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
        let phys_grp = self.file.group(&format!("results/{physics}"))?;
        let phys_attrs = phys_grp.attrs()?;
        let n_solutions = attr_i64(&phys_attrs, "n_solutions")?.max(0) as usize;

        let packed_eig_path = format!("results/{physics}/eigenvalues");
        if let Ok(ds_eig) = self.file.dataset(&packed_eig_path) {
            self.check_dataset_budget(&ds_eig, 8)?;
            let raw = ds_eig.read_f64()?;
            let labels_opt: Option<Vec<String>> =
                attr_str(&phys_attrs, "mode_labels_json").ok().and_then(|j| serde_json::from_str(&j).ok());
            let mut headers = Vec::with_capacity(n_solutions);
            for k in 0..n_solutions {
                let label = labels_opt
                    .as_ref()
                    .and_then(|v| v.get(k).cloned())
                    .or_else(|| {
                        self.file
                            .group(&format!("results/{physics}/solution_{k}"))
                            .ok()
                            .and_then(|g| g.attrs().ok())
                            .and_then(|a| attr_str(&a, "label").ok())
                    })
                    .unwrap_or_else(|| format!("mode_{k}"));
                headers.push(ModeHeader {
                    mode_index: k,
                    eig_re: raw.get(k * 2).copied().unwrap_or(0.0),
                    eig_im: raw.get(k * 2 + 1).copied().unwrap_or(0.0),
                    label,
                });
            }
            return Ok(headers);
        }

        let mut headers = Vec::with_capacity(n_solutions);
        for k in 0..n_solutions {
            let sol_grp = self.file.group(&format!("results/{physics}/solution_{k}"))?;
            let sattrs = sol_grp.attrs()?;
            let eig_re = attr_f64(&sattrs, "eig_re").or_else(|_| attr_f64(&sattrs, "value"))?;
            let eig_im = attr_f64(&sattrs, "eig_im").unwrap_or(0.0);
            let label = attr_str(&sattrs, "label").unwrap_or_default();
            headers.push(ModeHeader { mode_index: k, eig_re, eig_im, label });
        }
        Ok(headers)
    }

    /// Reads `/mesh` topology, connectivity, coordinates, and optional curved-edge geometry (`edge_node_pairs`, `edge_node_points`).
    /// Transparently handles both explicit `/mesh` attributes and legacy `nyx-results` v3/v4 files without `/mesh` attributes.
    pub fn read_mesh(&self) -> Result<MeshData> {
        let mesh_grp = self.file.group("mesh")?;
        let mattrs = mesh_grp.attrs()?;

        let ds_cl = self.file.dataset("mesh/cell_labels")?;
        self.check_dataset_budget(&ds_cl, 4)?;
        let cell_labels = ds_cl.read_i32()?;

        let ds_c = self.file.dataset("mesh/cells")?;
        self.check_dataset_budget(&ds_c, 4)?;
        let cs = ds_c.shape()?;
        let cells_shape = shape_2d(&cs, "mesh/cells")?;
        let cells = ds_c.read_u32()?;

        let ds_fl = self.file.dataset("mesh/facet_labels")?;
        self.check_dataset_budget(&ds_fl, 4)?;
        let facet_labels = ds_fl.read_i32()?;

        let ds_f = self.file.dataset("mesh/facets")?;
        self.check_dataset_budget(&ds_f, 4)?;
        let fs = ds_f.shape()?;
        let facets_shape = shape_2d(&fs, "mesh/facets")?;
        let facets = ds_f.read_u32()?;

        let ds_p = self.file.dataset("mesh/points")?;
        self.check_dataset_budget(&ds_p, 8)?;
        let ps = ds_p.shape()?;
        let points_shape = shape_2d(&ps, "mesh/points")?;
        let points = ds_p.read_f64()?;

        let mut edge_node_pairs = Vec::new();
        if let Ok(ds_ep) = self.file.dataset("mesh/edge_node_pairs") {
            self.check_dataset_budget(&ds_ep, 4)?;
            let raw_pairs = ds_ep.read_u32()?;
            let (chunks, rem) = raw_pairs.as_chunks::<2>();
            if !rem.is_empty() {
                bail!("mesh/edge_node_pairs length {} is not a multiple of 2", raw_pairs.len());
            }
            edge_node_pairs.extend_from_slice(chunks);
        }

        let mut edge_node_points = Vec::new();
        if let Ok(ds_ept) = self.file.dataset("mesh/edge_node_points") {
            self.check_dataset_budget(&ds_ept, 8)?;
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
            return Ok(Some([
                shape[0] as usize,
                shape[1] as usize,
                shape[2] as usize,
                shape[3] as usize,
                shape[4] as usize,
            ]));
        }
        if shape.len() == 1 {
            let attrs = ds.attrs()?;
            let tail = attrs
                .get("logical_shape_tail")
                .and_then(AttrValue::to_i64s)
                .ok_or_else(|| anyhow!("rank-1 cell_fields_packed missing logical_shape_tail attribute"))?;
            if tail.len() != 4 {
                bail!("invalid logical_shape_tail length {}: expected 4", tail.len());
            }
            let c = tail[0].max(0) as usize;
            let n_c = tail[1].max(0) as usize;
            let n_v = tail[2].max(0) as usize;
            let ri = tail[3].max(0) as usize;
            let mode_elems = c
                .checked_mul(n_c)
                .and_then(|x| x.checked_mul(n_v))
                .and_then(|x| x.checked_mul(ri))
                .ok_or_else(|| anyhow!("logical_shape_tail overflow"))?;
            let m = if let Some(div) = (shape[0] as usize).checked_div(mode_elems) {
                div
            } else {
                let pattrs = self.file.group(&format!("results/{physics}"))?.attrs()?;
                attr_i64(&pattrs, "n_solutions")?.max(0) as usize
            };
            return Ok(Some([m, c, n_c, n_v, ri]));
        }
        bail!("unexpected rank {} for {}", shape.len(), packed_path)
    }

    /// Bounded-memory single-mode slice reader (Gate 4 `EC-09` & `BM-2`).
    /// Dispatches to `ds.read_f64_rows` on `L2`/`L3`/`L1-packed` datasets (decompressing only chunk `mode_idx`),
    /// or reconstructs `[C, N_c, N_v, 2]` (`[re, im]` interleaved) from legacy `/results/<phys>/solution_<k>/cell_fields`.
    pub fn read_mode_cell_fields(&self, physics: &str, mode_idx: usize) -> Result<Vec<f64>> {
        let packed_path = format!("results/{physics}/cell_fields_packed");
        if let Ok(ds) = self.file.dataset(&packed_path) {
            let shape = ds.shape()?;
            if shape.len() == 5 {
                if (mode_idx as u64) >= shape[0] {
                    bail!("mode_idx {} out of bounds for shape {:?}", mode_idx, shape);
                }
                let row_elems = shape[1..]
                    .iter()
                    .try_fold(1u64, |acc, &d| acc.checked_mul(d))
                    .ok_or_else(|| anyhow!("row shape overflow: {:?}", shape))?;
                let row_bytes = row_elems.checked_mul(8).ok_or_else(|| anyhow!("row byte size overflow"))?;
                if row_bytes > self.budget.max_dataset_bytes {
                    bail!(
                        "single mode slice ({} B) exceeds read budget {} B",
                        row_bytes,
                        self.budget.max_dataset_bytes
                    );
                }
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
                let mode_elems = tail
                    .iter()
                    .try_fold(1u64, |acc, &d| acc.checked_mul(d.max(0) as u64))
                    .ok_or_else(|| anyhow!("mode_elems overflow"))?;
                let mode_bytes = mode_elems.checked_mul(8).ok_or_else(|| anyhow!("mode_bytes overflow"))?;
                if mode_bytes > self.budget.max_dataset_bytes {
                    bail!(
                        "single mode slice ({} B) exceeds read budget {} B",
                        mode_bytes,
                        self.budget.max_dataset_bytes
                    );
                }
                if mode_elems == 0 {
                    let pattrs = self.file.group(&format!("results/{physics}"))?.attrs()?;
                    let n_sol = attr_i64(&pattrs, "n_solutions")?.max(0) as usize;
                    if mode_idx >= n_sol {
                        bail!("mode_idx {} out of bounds for n_solutions={}", mode_idx, n_sol);
                    }
                    return Ok(Vec::new());
                }
                let start = (mode_idx as u64)
                    .checked_mul(mode_elems)
                    .ok_or_else(|| anyhow!("mode slice start offset overflow"))?;
                let end = start.checked_add(mode_elems).ok_or_else(|| anyhow!("mode slice end offset overflow"))?;
                if end > shape[0] {
                    bail!("mode_idx {} (elems {}..{}) out of bounds for len {}", mode_idx, start, end, shape[0]);
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
            let ds_re = self.file.dataset(&format!("{cf_path}/{cname}/re"))?;
            self.check_dataset_budget(&ds_re, 8)?;
            let re = ds_re.read_f64()?;

            let im = if let Ok(ds_im) = self.file.dataset(&format!("{cf_path}/{cname}/im")) {
                self.check_dataset_budget(&ds_im, 8)?;
                ds_im.read_f64()?
            } else {
                vec![0.0; re.len()]
            };
            if re.len() != im.len() {
                bail!("{}/{}: re.len() {} != im.len() {}", cf_path, cname, re.len(), im.len());
            }
            out.reserve(re.len() * 2);
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
        let root_attrs = self.file.root().attrs()?;
        let format = attr_str(&root_attrs, "format")?;
        let version = attr_i64(&root_attrs, "version")?;
        let writer = attr_str(&root_attrs, "writer")?;
        let processes = attr_i64(&root_attrs, "processes")?;

        let ds_proj = self.file.dataset("project")?;
        self.check_dataset_budget(&ds_proj, 1)?;
        let project_bytes = ds_proj.read_u8()?;

        let ds_cell_labels = self.file.dataset("mesh/cell_labels")?;
        self.check_dataset_budget(&ds_cell_labels, 4)?;
        let mesh_cell_labels = ds_cell_labels.read_i32()?;

        let ds_cells = self.file.dataset("mesh/cells")?;
        self.check_dataset_budget(&ds_cells, 4)?;
        let mesh_cells_shape = shape_2d(&ds_cells.shape()?, "mesh/cells")?;
        let mesh_cells = ds_cells.read_u32()?;

        let ds_facet_labels = self.file.dataset("mesh/facet_labels")?;
        self.check_dataset_budget(&ds_facet_labels, 4)?;
        let mesh_facet_labels = ds_facet_labels.read_i32()?;

        let ds_facets = self.file.dataset("mesh/facets")?;
        self.check_dataset_budget(&ds_facets, 4)?;
        let mesh_facets_shape = shape_2d(&ds_facets.shape()?, "mesh/facets")?;
        let mesh_facets = ds_facets.read_u32()?;

        let ds_points = self.file.dataset("mesh/points")?;
        self.check_dataset_budget(&ds_points, 8)?;
        let mesh_points_shape = shape_2d(&ds_points.shape()?, "mesh/points")?;
        let mesh_points = ds_points.read_f64()?;

        let phys_grp = self.file.group(&format!("results/{physics}"))?;
        let phys_attrs = phys_grp.attrs()?;
        let study = attr_str(&phys_attrs, "study")?;
        let dim = attr_i64(&phys_attrs, "dim")?;
        let dofs = attr_i64(&phys_attrs, "dofs")?;
        let n_solutions = attr_i64(&phys_attrs, "n_solutions")?;
        let order = attr_i64(&phys_attrs, "order")?;
        let seconds = attr_f64(&phys_attrs, "seconds")?;

        let ds_log = phys_grp.dataset("log")?;
        self.check_dataset_budget(&ds_log, 1)?;
        let log_bytes = ds_log.read_u8()?;

        let ds_scalars = phys_grp.dataset("scalars")?;
        self.check_dataset_budget(&ds_scalars, 1)?;
        let scalars_bytes = ds_scalars.read_u8()?;

        let viz_phys = self.file.group(&format!("viz/{physics}"))?;
        let viz_attrs = viz_phys.attrs()?;
        let viz_topology = attr_str(&viz_attrs, "topology")?;

        let ds_vcells = viz_phys.dataset("cells")?;
        self.check_dataset_budget(&ds_vcells, 4)?;
        let viz_cells_shape = shape_2d(&ds_vcells.shape()?, "viz/cells")?;
        let viz_cells = ds_vcells.read_u32()?;

        let ds_vpoints = viz_phys.dataset("points")?;
        self.check_dataset_budget(&ds_vpoints, 8)?;
        let viz_points_shape = shape_2d(&ds_vpoints.shape()?, "viz/points")?;
        let viz_points = ds_vpoints.read_f64()?;

        let n_sol = n_solutions.max(0) as usize;
        let mut solutions = Vec::with_capacity(n_sol);
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

                let ds_re = comp_grp.dataset("re")?;
                self.check_dataset_budget(&ds_re, 8)?;
                let re = ds_re.read_f64()?;

                let im = if let Ok(ds_im) = comp_grp.dataset("im") {
                    self.check_dataset_budget(&ds_im, 8)?;
                    ds_im.read_f64()?
                } else {
                    vec![0.0; re.len()]
                };
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
                self.check_dataset_budget(&ds, 8)?;
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
                    self.check_dataset_budget(&v_ds, 8)?;
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

fn shape_2d(shape: &[u64], label: &str) -> Result<[usize; 2]> {
    match shape {
        [r, c] => Ok([*r as usize, *c as usize]),
        [0] => Ok([0, 0]),
        _ => bail!("expected 2D shape for {}, got {:?}", label, shape),
    }
}

fn walk_group_checked(reader: &NyxReader, grp: &Group, depth: usize) -> Result<usize> {
    if depth > reader.budget.max_group_depth {
        bail!("group recursion depth {} exceeds configured limit {}", depth, reader.budget.max_group_depth);
    }
    let mut total_bytes = 0usize;
    let _ = grp.attrs()?;
    for ds_name in grp.datasets()? {
        let ds = grp.dataset(&ds_name)?;
        let _ = ds.attrs()?;
        let raw = reader.read_raw_checked(&ds)?;
        total_bytes = total_bytes.saturating_add(raw.len());
        if (total_bytes as u64) > reader.budget.max_total_result_bytes {
            bail!(
                "cumulative result size {} B exceeds configured total read budget {} B",
                total_bytes,
                reader.budget.max_total_result_bytes
            );
        }
    }
    for sub_name in grp.groups()? {
        let sub = grp.group(&sub_name)?;
        total_bytes = total_bytes.saturating_add(walk_group_checked(reader, &sub, depth + 1)?);
        if (total_bytes as u64) > reader.budget.max_total_result_bytes {
            bail!(
                "cumulative result size {} B exceeds configured total read budget {} B",
                total_bytes,
                reader.budget.max_total_result_bytes
            );
        }
    }
    Ok(total_bytes)
}

/// Walks an untrusted HDF5 byte buffer enforcing [`ReadBudgetConfig`] (`checked_mul` byte ceiling
/// and recursion depth limit) without panicking or OOM-aborting (`Gate 5 / EC-11`).
pub fn fuzz_walk_untrusted(bytes: &[u8], budget: ReadBudgetConfig) -> Result<usize> {
    let reader = NyxReader::from_bytes(bytes.to_vec(), budget)?;
    let root = reader.file().root();
    walk_group_checked(&reader, &root, 0)
}
