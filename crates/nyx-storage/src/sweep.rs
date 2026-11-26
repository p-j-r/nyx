// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Incremental in-place parameter-sweep HDF5 checkpointing and atomic Windows NTFS file replacement.

use crate::reader::{attr_i64, attr_str};
use crate::ReadBudgetConfig;
use anyhow::{anyhow, bail, Context, Result};
use hdf5_pure::{make_f64_type, AttrValue, Dataset, File, FileBuilder, MaxExtent};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Computes the sibling `.tmp` staging path in the same directory as `target_path`
/// (e.g. `foo.h5` -> `foo.h5.tmp`, `foo.xdmf` -> `foo.xdmf.tmp`).
pub(crate) fn staging_path_for(target_path: &Path) -> Result<PathBuf> {
    let fname = target_path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| anyhow!("target path {:?} has no valid UTF-8 filename", target_path))?;
    Ok(target_path.with_file_name(format!("{fname}.tmp")))
}

/// Atomically renames `staging_path` onto `target_path` (`MoveFileExW` semantics on Windows NTFS).
pub fn atomic_rename(staging_path: &Path, target_path: &Path) -> Result<()> {
    if let Some(parent) = target_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| format!("failed to create directory {:?}", parent))?;
        }
    }
    std::fs::rename(staging_path, target_path)
        .with_context(|| format!("atomic rename {:?} -> {:?} failed", staging_path, target_path))
}

/// Writes `bytes` to a sibling `.tmp` file, flushes OS buffers, closes the handle,
/// and atomically renames into `target_path` so concurrent readers never see a torn file.
pub fn atomic_write_bytes(target_path: &Path, bytes: &[u8]) -> Result<()> {
    if let Some(parent) = target_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| format!("failed to create directory {:?}", parent))?;
        }
    }
    let tmp_path = staging_path_for(target_path)?;
    {
        let mut f = std::fs::File::create(&tmp_path)
            .with_context(|| format!("failed to create staging file {:?}", tmp_path))?;
        f.write_all(bytes).with_context(|| format!("failed to write staging file {:?}", tmp_path))?;
        f.sync_all().with_context(|| format!("failed to sync staging file {:?}", tmp_path))?;
    }
    atomic_rename(&tmp_path, target_path)
}

/// One parameter-sweep evaluation point (`params` + complex mode eigenvalues).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SweepPoint {
    pub params: Vec<f64>,
    pub eig_re: Vec<f64>,
    pub eig_im: Vec<f64>,
}

/// Metadata for a parametric sweep sample skipped due to invalid geometry, meshing failure, or bad elements.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SkippedSweepPoint {
    pub index: usize,
    pub label: String,
    pub params: Vec<f64>,
    pub reason: String,
}

/// Incremental HDF5 parameter-sweep checkpoint store backed by unlimited rank-1 chunked datasets.
#[derive(Clone, Debug)]
pub struct SweepCheckpoint {
    path: PathBuf,
    param_names: Vec<String>,
    n_modes: usize,
    use_lzf: bool,
}

impl SweepCheckpoint {
    /// Default chunk size (in sweep points) along axis 0 for each 1D column dataset.
    pub const DEFAULT_CHUNK_POINTS: u64 = 16;

    /// Initializes a new `.h5` sweep checkpoint file atomically with unlimited chunked 1D columns:
    /// one dataset per entry in `param_names` plus `eig_re_<m>` and `eig_im_<m>` for `m in 0..n_modes`.
    pub fn create(path: &Path, param_names: &[&str], n_modes: usize) -> Result<Self> {
        Self::create_with_options(path, param_names, n_modes, false)
    }

    /// Initializes a new `.h5` sweep checkpoint file with optional `Shuffle + LZF` chunk compression.
    pub fn create_with_options(path: &Path, param_names: &[&str], n_modes: usize, use_lzf: bool) -> Result<Self> {
        if param_names.is_empty() {
            bail!("SweepCheckpoint requires at least one parameter column name");
        }
        let owned_names: Vec<String> = param_names.iter().map(|&s| s.to_string()).collect();
        let params_json = serde_json::to_string(&owned_names)?;

        let mut fb = FileBuilder::new();
        fb.set_attr("format", AttrValue::VarLenString("nyx-sweep".into()));
        fb.set_attr("version", AttrValue::I64(1));
        fb.set_attr("param_names_json", AttrValue::VarLenString(params_json));
        fb.set_attr("n_modes", AttrValue::I64(n_modes as i64));
        fb.set_attr("n_points", AttrValue::I64(0));
        fb.set_attr("use_lzf", AttrValue::I64(i64::from(use_lzf)));

        let configure_col = |fb: &mut FileBuilder, col_name: &str| {
            let ds = fb
                .create_dataset(col_name)
                .with_dtype(make_f64_type())
                .with_shape(&[0])
                .with_maxshape(&[MaxExtent::Unlimited])
                .with_chunks(&[Self::DEFAULT_CHUNK_POINTS]);
            if use_lzf {
                ds.with_shuffle().with_lzf();
            }
        };

        for pname in &owned_names {
            configure_col(&mut fb, pname);
        }
        for m in 0..n_modes {
            configure_col(&mut fb, &format!("eig_re_{m}"));
            configure_col(&mut fb, &format!("eig_im_{m}"));
        }

        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let tmp_path = staging_path_for(path)?;
        fb.write(&tmp_path)?;
        atomic_rename(&tmp_path, path)?;

        Ok(Self { path: path.to_path_buf(), param_names: owned_names, n_modes, use_lzf })
    }

    /// Opens an existing `.h5` sweep checkpoint file and reads its column schema.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open_streaming(path)?;
        let attrs = file.root().attrs()?;
        let format = attr_str(&attrs, "format")?;
        if format != "nyx-sweep" {
            bail!("unexpected sweep checkpoint format {:?}", format);
        }
        let params_json = attr_str(&attrs, "param_names_json")?;
        let param_names: Vec<String> = serde_json::from_str(&params_json)?;
        let n_modes = attr_i64(&attrs, "n_modes")?.max(0) as usize;
        let use_lzf = attr_i64(&attrs, "use_lzf").unwrap_or(0) != 0;
        Ok(Self { path: path.to_path_buf(), param_names, n_modes, use_lzf })
    }

    /// Path to the backing `.h5` file.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Parameter column names in order.
    pub fn param_names(&self) -> &[String] {
        &self.param_names
    }

    /// Number of tracked eigenmodes per sweep point.
    pub fn n_modes(&self) -> usize {
        self.n_modes
    }

    /// Appends a single solved sweep point in place via `File::open_rw`.
    pub fn append_point(&self, point: &SweepPoint) -> Result<()> {
        self.append_points(std::slice::from_ref(point))
    }

    /// Appends one or more solved sweep points in place via `File::open_rw`.
    /// Uses immediate `Dataset::append` for unfiltered columns when no staged metadata update is needed,
    /// or transactional `append_staged` + `commit()` (updating `@n_points`) for all configurations.
    pub fn append_points(&self, points: &[SweepPoint]) -> Result<()> {
        if points.is_empty() {
            return Ok(());
        }
        for (idx, pt) in points.iter().enumerate() {
            if pt.params.len() != self.param_names.len() {
                bail!("sweep point {} has {} params, expected {}", idx, pt.params.len(), self.param_names.len());
            }
            if pt.eig_re.len() != self.n_modes || pt.eig_im.len() != self.n_modes {
                bail!(
                    "sweep point {} has ({}, {}) eigenvalues, expected {}",
                    idx,
                    pt.eig_re.len(),
                    pt.eig_im.len(),
                    self.n_modes
                );
            }
        }

        let session = File::open_rw(&self.path)?;
        let cur_len = session.dataset(&self.param_names[0])?.shape()?.first().copied().unwrap_or(0);
        let new_len = cur_len.checked_add(points.len() as u64).ok_or_else(|| anyhow!("sweep point count overflow"))?;

        if !self.use_lzf {
            // Immediate crash-atomic Extensible-Array append on each unfiltered column,
            // followed by staged `@n_points` attribute update + commit.
            for (p_idx, pname) in self.param_names.iter().enumerate() {
                let vals: Vec<f64> = points.iter().map(|pt| pt.params[p_idx]).collect();
                let mut ds = session.dataset(pname)?;
                ds.append(&vals)?;
            }
            for m in 0..self.n_modes {
                let re_vals: Vec<f64> = points.iter().map(|pt| pt.eig_re[m]).collect();
                let im_vals: Vec<f64> = points.iter().map(|pt| pt.eig_im[m]).collect();
                session.dataset(&format!("eig_re_{m}"))?.append(&re_vals)?;
                session.dataset(&format!("eig_im_{m}"))?.append(&im_vals)?;
            }
            session.root().set_attr("n_points", AttrValue::I64(new_len as i64))?;
            session.commit()?;
        } else {
            for (p_idx, pname) in self.param_names.iter().enumerate() {
                let vals: Vec<f64> = points.iter().map(|pt| pt.params[p_idx]).collect();
                let mut ds = session.dataset(pname)?;
                ds.append_staged(|b| {
                    b.append_f64(&vals);
                })?;
            }
            for m in 0..self.n_modes {
                let re_vals: Vec<f64> = points.iter().map(|pt| pt.eig_re[m]).collect();
                let im_vals: Vec<f64> = points.iter().map(|pt| pt.eig_im[m]).collect();
                session.dataset(&format!("eig_re_{m}"))?.append_staged(|b| {
                    b.append_f64(&re_vals);
                })?;
                session.dataset(&format!("eig_im_{m}"))?.append_staged(|b| {
                    b.append_f64(&im_vals);
                })?;
            }
            session.root().set_attr("n_points", AttrValue::I64(new_len as i64))?;
            session.commit()?;
        }
        drop(session);
        Ok(())
    }

    /// Reads all stored sweep points back in streaming mode with [`ReadBudgetConfig`] enforcement.
    pub fn read_points(&self) -> Result<Vec<SweepPoint>> {
        self.read_points_with_budget(ReadBudgetConfig::default())
    }

    /// Reads all stored sweep points back under an explicit [`ReadBudgetConfig`].
    pub fn read_points_with_budget(&self, budget: ReadBudgetConfig) -> Result<Vec<SweepPoint>> {
        let file = File::open_streaming(&self.path)?;
        let read_col = |name: &str| -> Result<Vec<f64>> {
            let ds: Dataset = file.dataset(name)?;
            let shape = ds.shape()?;
            let elems = shape
                .iter()
                .try_fold(1u64, |acc, &d| acc.checked_mul(d))
                .ok_or_else(|| anyhow!("shape overflow on sweep column '{}'", name))?;
            let bytes = elems.checked_mul(8).ok_or_else(|| anyhow!("byte size overflow"))?;
            if bytes > budget.max_dataset_bytes {
                bail!("sweep column '{}' ({} B) exceeds read budget {} B", name, bytes, budget.max_dataset_bytes);
            }
            Ok(ds.read_f64()?)
        };

        let mut param_cols = Vec::with_capacity(self.param_names.len());
        for pname in &self.param_names {
            param_cols.push(read_col(pname)?);
        }
        let n_pts = param_cols.first().map_or(0, Vec::len);
        for (i, col) in param_cols.iter().enumerate() {
            if col.len() != n_pts {
                bail!("sweep column '{}' length {} != {}", self.param_names[i], col.len(), n_pts);
            }
        }

        let mut re_cols = Vec::with_capacity(self.n_modes);
        let mut im_cols = Vec::with_capacity(self.n_modes);
        for m in 0..self.n_modes {
            let rc = read_col(&format!("eig_re_{m}"))?;
            let ic = read_col(&format!("eig_im_{m}"))?;
            if rc.len() != n_pts || ic.len() != n_pts {
                bail!("eigenvalue column {} length mismatch ({}, {}) vs {}", m, rc.len(), ic.len(), n_pts);
            }
            re_cols.push(rc);
            im_cols.push(ic);
        }

        let mut out = Vec::with_capacity(n_pts);
        for row in 0..n_pts {
            let params = param_cols.iter().map(|c| c[row]).collect();
            let eig_re = re_cols.iter().map(|c| c[row]).collect();
            let eig_im = im_cols.iter().map(|c| c[row]).collect();
            out.push(SweepPoint { params, eig_re, eig_im });
        }
        Ok(out)
    }

    /// Returns the committed number of solved points (`@n_points`).
    pub fn n_points(&self) -> Result<u64> {
        let file = File::open_streaming(&self.path)?;
        let attrs = file.root().attrs()?;
        Ok(attr_i64(&attrs, "n_points").unwrap_or(0).max(0) as u64)
    }

    /// Returns the committed number of skipped points (`@n_skipped`, or `0` on legacy files).
    pub fn n_skipped(&self) -> Result<u64> {
        let file = File::open_streaming(&self.path)?;
        let attrs = file.root().attrs()?;
        Ok(attr_i64(&attrs, "n_skipped").unwrap_or(0).max(0) as u64)
    }

    /// Appends one skipped sweep point (`skipped_invalid_geometry`) to root metadata cleanly in place.
    pub fn record_skipped(&self, skipped: &SkippedSweepPoint) -> Result<()> {
        self.record_skipped_batch(std::slice::from_ref(skipped))
    }

    /// Appends one or more skipped sweep points atomically via `File::open_rw` + `commit()`.
    pub fn record_skipped_batch(&self, items: &[SkippedSweepPoint]) -> Result<()> {
        if items.is_empty() {
            return Ok(());
        }
        for (i, item) in items.iter().enumerate() {
            if item.params.len() != self.param_names.len() {
                bail!("skipped point {} has {} params, expected {}", i, item.params.len(), self.param_names.len());
            }
        }
        let mut existing = self.read_skipped()?;
        existing.extend(items.iter().cloned());
        let encoded = serde_json::to_string(&existing)?;
        let n_skipped = existing.len() as i64;

        let session = File::open_rw(&self.path)?;
        session.root().set_attr("n_skipped", AttrValue::I64(n_skipped))?;
        session.root().set_attr("skipped_json", AttrValue::VarLenString(encoded))?;
        session.commit()?;
        drop(session);
        Ok(())
    }

    /// Reads all recorded skipped sweep points (empty vector on checkpoints with no skipped entries).
    pub fn read_skipped(&self) -> Result<Vec<SkippedSweepPoint>> {
        let file = File::open_streaming(&self.path)?;
        let attrs = file.root().attrs()?;
        let Ok(raw) = attr_str(&attrs, "skipped_json") else {
            return Ok(Vec::new());
        };
        let list: Vec<SkippedSweepPoint> = serde_json::from_str(&raw)?;
        Ok(list)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sweep_checkpoint_streams_points_and_skipped_invalid_geometry_cleanly() {
        let dir = std::env::temp_dir().join(format!("nyx-storage-sweep-uq-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("qmc_ckpt.h5");

        let ckpt = SweepCheckpoint::create(&path, &["gap_um", "eps_r"], 0).unwrap();
        assert_eq!(ckpt.n_points().unwrap(), 0);
        assert_eq!(ckpt.n_skipped().unwrap(), 0);

        ckpt.append_point(&SweepPoint { params: vec![5.0, 11.45], eig_re: vec![], eig_im: vec![] }).unwrap();
        ckpt.record_skipped(&SkippedSweepPoint {
            index: 1,
            label: "qmc-0001".into(),
            params: vec![0.0, 11.45],
            reason: "gap_um <= 0 invalid geometry".into(),
        })
        .unwrap();
        ckpt.append_point(&SweepPoint { params: vec![7.5, 11.70], eig_re: vec![], eig_im: vec![] }).unwrap();

        let reopened = SweepCheckpoint::open(&path).unwrap();
        assert_eq!(reopened.n_points().unwrap(), 2);
        assert_eq!(reopened.n_skipped().unwrap(), 1);
        let pts = reopened.read_points().unwrap();
        assert_eq!(pts.len(), 2);
        assert_eq!(pts[0].params, vec![5.0, 11.45]);
        assert_eq!(pts[1].params, vec![7.5, 11.70]);
        let skipped = reopened.read_skipped().unwrap();
        assert_eq!(skipped.len(), 1);
        assert_eq!(skipped[0].index, 1);
        assert_eq!(skipped[0].label, "qmc-0001");
        assert!(skipped[0].reason.contains("invalid geometry"));

        let _ = std::fs::remove_dir_all(&dir);
    }
}
