// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-Rust HDF5 storage engine types, read-budget guards, and layout configuration for `nyx`.

pub mod reader;
pub mod sweep;
pub mod writer;
pub mod xdmf;

use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

pub use reader::{fuzz_walk_untrusted, NyxReader, NyxSolutionData, NyxV3Summary};
pub use sweep::{atomic_rename, atomic_write_bytes, SweepCheckpoint, SweepPoint};
pub use writer::{append_solved_mode_in_place, init_nyx_packed_streaming, write_nyx_packed, write_nyx_v3};
pub use xdmf::{generate_direct_xdmf, viz_link_name, write_direct_xdmf};

/// Dataset storage layout for FEM cell and nodal fields.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum FieldStorageLayout {
    /// Legacy `nyx-results` per-mode, per-component datasets compressed with `gzip-3`.
    /// Use only when strict structural compatibility with legacy per-field scripts is required.
    L1LegacyGzip3,
    /// Recommended Default: Packed 5D `[M, C, N_c, N_v, 2]` chunked per mode `[1, C, N_c, N_v, 2]`
    /// with `Shuffle (id=2) + LZF (id=32000)`.
    /// Delivers high compression on smooth P2 modes, fast writes, and native `h5py` compatibility.
    #[default]
    L2ShuffleLzfPacked,
    /// Maximum-speed scratch/hot-loop layout: Packed 5D `[M, C, N_c, N_v, 2]` contiguous uncompressed.
    /// Enables single-positional-`pread` single-mode slicing.
    L3ContiguousPacked,
}

impl FieldStorageLayout {
    /// Stable string identifier stored in HDF5 `/results/<physics>` `@field_layout` attributes.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::L1LegacyGzip3 => "L1LegacyGzip3",
            Self::L2ShuffleLzfPacked => "L2ShuffleLzfPacked",
            Self::L3ContiguousPacked => "L3ContiguousPacked",
        }
    }

    /// Parses a `@field_layout` attribute string back into [`FieldStorageLayout`].
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "L1LegacyGzip3" => Some(Self::L1LegacyGzip3),
            "L2ShuffleLzfPacked" => Some(Self::L2ShuffleLzfPacked),
            "L3ContiguousPacked" => Some(Self::L3ContiguousPacked),
            _ => None,
        }
    }
}

/// Security & memory bounds enforced when opening untrusted or external `.h5` files (Gate EC-11).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadBudgetConfig {
    /// Maximum allowed bytes for a single whole-dataset read call (default: 512 MiB).
    /// Prevents corrupted v1 dataspace dimensions from triggering multi-TB allocations.
    pub max_dataset_bytes: u64,
    /// Maximum cumulative bytes across all datasets when materializing a full result file in memory
    /// (default: 8 GiB, overridable via `NYX_H5_MAX_BYTES`).
    pub max_total_result_bytes: u64,
    /// Maximum group traversal recursion depth (default: 24).
    pub max_group_depth: usize,
}

impl Default for ReadBudgetConfig {
    fn default() -> Self {
        Self {
            max_dataset_bytes: 512 * 1024 * 1024,
            max_total_result_bytes: 8 * 1024 * 1024 * 1024,
            max_group_depth: 24,
        }
    }
}

impl ReadBudgetConfig {
    /// Returns [`ReadBudgetConfig::default()`], overriding `max_total_result_bytes` when the
    /// `NYX_H5_MAX_BYTES` environment variable is set to a positive integer byte count
    /// (and clamping `max_dataset_bytes` to `max_total_result_bytes` when `max_dataset_bytes > max_total_result_bytes`).
    pub fn from_env() -> Self {
        let mut cfg = Self::default();
        if let Some(limit) = std::env::var("NYX_H5_MAX_BYTES")
            .ok()
            .and_then(|s| s.trim().parse::<u64>().ok())
            .filter(|&n| n > 0)
        {
            cfg.max_total_result_bytes = limit;
            if cfg.max_dataset_bytes > limit {
                cfg.max_dataset_bytes = limit;
            }
        }
        cfg
    }
}

/// Complete mesh topology, connectivity, coordinates, and optional curved-edge geometry.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MeshData {
    pub cell_topology: String,
    pub element_family: String,
    pub element_degree: i64,
    pub space_dimension: i64,
    pub cell_labels: Vec<i32>,
    pub cells_shape: [usize; 2],
    pub cells: Vec<u32>,
    pub facet_labels: Vec<i32>,
    pub facets_shape: [usize; 2],
    pub facets: Vec<u32>,
    pub points_shape: [usize; 2],
    pub points: Vec<f64>,
    #[serde(default)]
    pub edge_node_pairs: Vec<[u32; 2]>,
    #[serde(default)]
    pub edge_node_points: Vec<[f64; 3]>,
}

/// Eigenvalue and label summary for a single mode or solution.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModeHeader {
    pub mode_index: usize,
    pub eig_re: f64,
    pub eig_im: f64,
    pub label: String,
}

/// Metadata for a single packed cell-field component stored in `/results/<physics>/cell_fields_meta`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PackedComponentMeta {
    pub name: String,
    pub unit: String,
    pub nloc: usize,
    pub has_im: bool,
}

/// Exact 0-ULP (`to_bits()`) comparison of two `f64` slices.
pub fn assert_f64_bits_eq(context: &str, actual: &[f64], expected: &[f64]) -> Result<()> {
    if actual.len() != expected.len() {
        bail!("{}: length mismatch (actual {}, expected {})", context, actual.len(), expected.len());
    }
    for (i, (&a, &e)) in actual.iter().zip(expected.iter()).enumerate() {
        if a.to_bits() != e.to_bits() {
            bail!(
                "{}[{}]: 0-ULP mismatch: actual={} (0x{:016x}) vs expected={} (0x{:016x})",
                context,
                i,
                a,
                a.to_bits(),
                e,
                e.to_bits()
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_string_round_trip() {
        for l in [
            FieldStorageLayout::L1LegacyGzip3,
            FieldStorageLayout::L2ShuffleLzfPacked,
            FieldStorageLayout::L3ContiguousPacked,
        ] {
            assert_eq!(FieldStorageLayout::parse(l.as_str()), Some(l));
        }
        assert_eq!(FieldStorageLayout::default(), FieldStorageLayout::L2ShuffleLzfPacked);
        assert_eq!(FieldStorageLayout::parse("unknown"), None);
    }

    #[test]
    fn f64_bits_eq_distinguishes_signed_zero_and_nan_payloads() {
        assert!(assert_f64_bits_eq("pos_zero", &[0.0], &[0.0]).is_ok());
        assert!(assert_f64_bits_eq("sign_zero", &[0.0], &[-0.0]).is_err());
        let nan1 = f64::from_bits(0x7ff8_0000_0000_1234);
        let nan2 = f64::from_bits(0x7ff8_0000_0000_5678);
        assert!(assert_f64_bits_eq("same_nan", &[nan1], &[nan1]).is_ok());
        assert!(assert_f64_bits_eq("diff_nan", &[nan1], &[nan2]).is_err());
    }

    #[test]
    fn viz_link_sanitization() {
        assert_eq!(viz_link_name("00_|E|"), "_E_");
        assert_eq!(viz_link_name("03_Re Ex"), "Re_Ex");
        assert_eq!(viz_link_name("plain_name"), "plain_name");
    }

    #[test]
    fn read_budget_config_defaults_and_env() {
        let def = ReadBudgetConfig::default();
        assert_eq!(def.max_dataset_bytes, 512 * 1024 * 1024);
        assert_eq!(def.max_total_result_bytes, 8 * 1024 * 1024 * 1024);
        assert_eq!(def.max_group_depth, 24);
        let env_cfg = ReadBudgetConfig::from_env();
        assert!(env_cfg.max_total_result_bytes > 0);
    }
}
